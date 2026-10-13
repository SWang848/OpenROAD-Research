// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include <chrono>

#include "tdp/Timingdp.h"
#include "infrastructure/Architecture3D.h"
#include "infrastructure/Network3D.h"
#include "infrastructure/Common3D.h"
#include "infrastructure/Net3D.h"
#include "infrastructure/SegmentManager3D.h"
#include "infrastructure/Segment3D.h"
#include "infrastructure/Node3D.h"
#include "eva/eva.h"
#include "eva/TimingEngine3D.h"
#include "hbta/HBTACostModel.h"
#include "placement/placement.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include "grt/GRoute.h"
#include "grt/GlobalRouter.h"
#include "sta/Sta.hh"
#include "sta/Search.hh"
#include "sta/Graph.hh"
#include "sta/Path.hh"
#include "sta/Network.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/TimingArc.hh"
#include "sta/TimingRole.hh"
#include "sta/Delay.hh"
#include "sta/MinMax.hh"
#include "rsz/Resizer.hh"
#include "odb/db.h"
#include "odb/defout.h"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "utl/Logger.h"
#include "lgista/LocalGraphIsta.hh"
#include "PtGraph.hh"

using utl::TDP;

namespace tdp {

// ============================================================================
// Timingdp实现（使用新的转换系统）
// ============================================================================

Timingdp::Timingdp(odb::dbDatabase* db,
                   utl::Logger* logger, 
                   sta::dbSta* sta, 
                   rsz::Resizer* resizer, 
                   grt::GlobalRouter* global_router,
                   dpl::Opendp* opendp) 
    : db_(db), sta_(sta), logger_(logger), 
      rs_(resizer), global_router_(global_router),
      estimate_parasitics_(resizer->getEstimateParasitics()), 
      opendp_(opendp),
      segment_manager_(nullptr) {
    // 初始化3D架构和网络（类似DPL的Opendp）
    timing_engine_ = std::make_unique<TimingEngine3D>(sta_, this);
    arch_3d_ = std::make_unique<Architecture3D>();
    network_3d_ = std::make_unique<Network3D>();
}

void Timingdp::setSta(sta::dbSta* sta) {
    sta_ = sta;
}

Timingdp::~Timingdp() {
    if (segment_manager_) {
        delete segment_manager_;
        segment_manager_ = nullptr;
    }
    destroyDatabase();
}

// ============================================================================
// 工艺组合配置（-tech）
// ============================================================================
Timingdp::TechProfile Timingdp::parseProfile(const std::string& tech) {
    // 命名规则：<upper工艺>_<bottom工艺>_3D；同构时合并为 <工艺>_3D。
    if (tech == "asap7_3D")           return {Tech::ASAP7,     Tech::ASAP7};
    if (tech == "nangate45_3D")       return {Tech::NANGATE45, Tech::NANGATE45};
    if (tech == "asap7_nangate45_3D") return {Tech::ASAP7,     Tech::NANGATE45};
    if (tech == "nangate45_asap7_3D") return {Tech::NANGATE45, Tech::ASAP7};
    // 默认 = 现状（上 asap7 / 下 nangate45）。
    return {Tech::ASAP7, Tech::NANGATE45};
}

const char* Timingdp::techSiteName(Tech t) {
    return (t == Tech::ASAP7) ? "asap7sc7p5t"
                              : "FreePDK45_38x28_10R_NP_162NW_34O";
}

const char* Timingdp::tierSite(const std::string& tier, const std::string& tech) {
    TechProfile p = parseProfile(tech);
    Tech t = (tier == "upper") ? p.upper : p.bottom;
    return techSiteName(t);
}

// ============================================================================
// 数据库管理（类似DPL的Opendp）
// ============================================================================
void Timingdp::detailedPlacement(const int max_displacement_x,
                               const int max_displacement_y,
                               const std::string& report_file_name,
                               double margin_um,
                               const std::string& tech,
                               const std::string& hbt)
{
  // 解析 -tech，确定上/下层各自工艺（供 legalizeTier 选 site）。
  profile_ = parseProfile(tech);

  auto tp0 = std::chrono::steady_clock::now();
  importDb();
  auto tp1 = std::chrono::steady_clock::now();
//   adjustNodesOrient();
//   for (const auto& node : network_3d_->getNode()) {
//     if (node->getType() == Node::CELL && !node->isFixed()) {
//       node->setPlaced(false);
//     }
//   }

    if (max_displacement_x == 0 || max_displacement_y == 0) {
    // defaults
    max_displacement_x_ = 500;
    max_displacement_y_ = 100;
    } else {
    max_displacement_x_ = max_displacement_x;
    max_displacement_y_ = max_displacement_y;
    }

    logger_->info(TDP, 6, "import right");

    std::string hbt_array_json = hbt;
    if (!initializeHBTAllocator(hbt_array_json, block_->getDbUnitsPerMicron())) {
        std::cerr << "WARNING: Failed to initialize HBT allocator, but continuing..." << std::endl;
    } else {
        std::cout << "✓ HBT allocator initialized successfully" << std::endl;

        // 执行HBT分配
        auto allocation_result = allocateHBTsForDesign();
        
        std::cout << "   - Total inter-layer nets: " << allocation_result.total_nets << std::endl;
        std::cout << "   - Successful allocations: " << allocation_result.successful_allocations << std::endl;
        std::cout << "   - Failed allocations: " << allocation_result.failed_allocations << std::endl;
        
        if (allocation_result.total_nets > 0) {
            std::cout << "   - Success rate: " 
                        << (static_cast<double>(allocation_result.successful_allocations) / allocation_result.total_nets * 100.0)
                        << "%" << std::endl;
            std::cout << "   - Total estimated wirelength: " << allocation_result.total_wirelength << " microns" << std::endl;
        }
    }
    rs_->resizeSlackPreamble();

    auto tp2 = std::chrono::steady_clock::now();

    // 先初始化 GRT 的 grid 和 db_net_map_
    int min_layer, max_layer;
    global_router_->getMinMaxLayer(min_layer, max_layer);
    // TDP uses the strict max routing layer from its source implementation.
    // The explicit zero keeps this policy local to TDP; ordinary GRT callers
    // retain TaiWei's default two-layer pin extension.
    global_router_->initFastRoute(min_layer, max_layer, true, 0);
    SteinerTreeGenerator steiner_generator(this);
    steiner_generator.createSteninerRoute();
    executeTimingDriven(true);
    auto stats = timing_engine_->getTimingStats();

    auto tp3 = std::chrono::steady_clock::now();

    Placement3D placetdp(this);
    logger_->info(TDP, 112, "start rundetailedplacement");
    placetdp.rundetailedplacement();

    auto tp4 = std::chrono::steady_clock::now();

    steiner_generator.createSteninerRoute();
    sta_->updateTiming(true);  
    auto stats2 = timing_engine_->getTimingStats();
    logger_->info(TDP, 1007,
        "Early: WNS={:.20f}ns TNS={:.6f}ns -> WNS={:.6f}ns (Δ={:.6f}) TNS={:.6f}ns (Δ={:.6f})",
        stats.early_wns * 1e9, stats.early_tns * 1e9,
        stats2.early_wns * 1e9, (stats2.early_wns - stats.early_wns)/stats.early_wns,
        stats2.early_tns * 1e9 , (stats2.early_tns - stats.early_tns)/stats.early_tns);
    logger_->info(TDP, 1008,
        "Late : WNS={:.6f}ns TNS={:.6f}ns -> WNS={:.6f}ns (Δ={:.6f}) TNS={:.6f}ns (Δ={:.6f})",
        stats.late_wns * 1e9, stats.late_tns * 1e9,
        stats2.late_wns * 1e9, (stats2.late_wns - stats.late_wns)/stats.late_wns,
        stats2.late_tns * 1e9, (stats2.late_tns - stats.late_tns)/stats.late_tns);

    legalizeTier("bottom", margin_um, max_displacement_x_, max_displacement_y_);
    legalizeTier("upper", margin_um, max_displacement_x_, max_displacement_y_);

    auto tp5 = std::chrono::steady_clock::now();

    steiner_generator.createSteninerRoute();
    sta_->updateTiming(true);  
    auto stats3 = timing_engine_->getTimingStats();
    logger_->info(TDP, 1009,
        "Late : WNS={:.20f} TNS={:.20f} -> WNS={:.20f} (Δ={:.6f}) TNS={:.20f} (Δ={:.6f})",
        stats.late_wns, stats.late_tns,
        stats3.late_wns, (stats3.late_wns - stats.late_wns)/stats.late_wns,
        stats3.late_tns, (stats3.late_tns - stats.late_tns)/stats.late_tns);

#define DUR_MS(a,b) std::chrono::duration<double>(b-a).count()
    logger_->info(TDP, 3002,
        "[PHASE] importDb={:.1f}s HBT={:.1f}s Steiner+STA={:.1f}s DLS={:.1f}s postDLS+legalize={:.1f}s",
        DUR_MS(tp0,tp1), DUR_MS(tp1,tp2), DUR_MS(tp2,tp3), DUR_MS(tp3,tp4), DUR_MS(tp4,tp5));
#undef DUR_MS
}

// P2 diagnostic (env TDP_DIAG_P2=<N>, default N=10, default OFF): dump the
// worst-N setup paths arc-by-arc, splitting gate(cell) vs wire(net) delay via
// Edge::isWire(). For each gate arc prints lib cell / arc / input slew / output
// load cap / delay / out slew so it can be cross-checked against the .lib NLDM
// tables; for each wire arc prints net / driver / load cap / delay. Log only,
// does not alter timing or placement.
static void diagWorstSetupPaths(TimingEngine3D* te, utl::Logger* logger)
{
  const char* env = std::getenv("TDP_DIAG_P2");
  if (!env) return;
  int topN = std::atoi(env);
  if (topN <= 0) topN = 10;

  sta::dbSta* sta = te->getSta();
  sta::dbNetwork* nw = sta->getDbNetwork();
  sta::GraphDelayCalc* gdc = sta->graphDelayCalc();
  sta::Scene* scene = te->getDefaultScene();
  const sta::MinMax* mm = sta::MinMax::max();
  if (!sta || !nw || !gdc || !scene) return;

  auto paths = te->getViolatingPaths();  // setup, sorted worst-first
  const int np = std::min<int>(topN, static_cast<int>(paths.size()));
  logger->info(TDP, 2075, "[P2] dumping worst {} setup paths (units: ns, fF)", np);

  auto pin_name = [&](sta::Pin* p) -> std::string {
    return p ? nw->pathName(p) : std::string("?");
  };

  for (int i = 0; i < np; ++i) {
    sta::Path* endp = paths[i];
    if (!endp) continue;

    // Collect the chain endpoint..startpoint, then reverse to print start->end.
    std::vector<sta::Path*> chain;
    for (sta::Path* s = endp; s != nullptr; s = s->prevPath())
      chain.push_back(s);
    std::reverse(chain.begin(), chain.end());
    if (chain.size() < 2) continue;

    logger->info(TDP, 2076,
      "[P2] path#{}  start={}  end={}  arrival={:.6f}ns required={:.6f}ns slack={:.6f}ns",
      i, pin_name(chain.front()->pin(sta)), pin_name(endp->pin(sta)),
      sta::delayAsFloat(endp->arrival()) * 1e9,
      sta::delayAsFloat(endp->required()) * 1e9,
      sta::delayAsFloat(endp->slack(sta)) * 1e9);

    double gate_tot = 0.0, wire_tot = 0.0;
    for (size_t k = 1; k < chain.size(); ++k) {
      sta::Path* prev = chain[k - 1];
      sta::Path* step = chain[k];
      const double d = sta::delayAsFloat(step->arrival())
                     - sta::delayAsFloat(prev->arrival());
      sta::Edge* e = step->prevEdge(sta);
      sta::TimingArc* arc = step->prevArc(sta);
      bool is_wire = e ? e->isWire()
                       : (arc && arc->role() && arc->role()->isWire());
      const sta::RiseFall* rf = step->transition(sta);
      const std::string rf_s = rf ? rf->to_string(true) : std::string("?");
      const double in_slew = sta::delayAsFloat(prev->slew(sta)) * 1e9;   // ns
      const double out_slew = sta::delayAsFloat(step->slew(sta)) * 1e9;  // ns
      sta::Pin* from_pin = prev->pin(sta);
      sta::Pin* to_pin = step->pin(sta);

      if (is_wire) {
        sta::Net* net = to_pin ? nw->net(to_pin) : nullptr;
        const std::string net_name = net ? nw->pathName(net) : std::string("?");
        const float load_cap = from_pin
            ? gdc->loadCap(from_pin, rf, scene, mm) : 0.0f;
        logger->info(TDP, 2077,
          "[P2]   WIRE net={} drvr={} load={} load_cap={:.4f}fF delay={:.6f}ns",
          net_name, pin_name(from_pin), pin_name(to_pin),
          load_cap * 1e15, d * 1e9);
        wire_tot += d;
      } else {
        sta::Instance* inst = to_pin ? nw->instance(to_pin) : nullptr;
        sta::LibertyCell* lc = inst ? nw->libertyCell(inst) : nullptr;
        const std::string cell_name = lc ? lc->name() : std::string("?");
        const std::string from_port =
            (arc && arc->from()) ? arc->from()->name() : std::string("?");
        const std::string to_port =
            (arc && arc->to()) ? arc->to()->name() : std::string("?");
        const float load_cap = to_pin
            ? gdc->loadCap(to_pin, rf, scene, mm) : 0.0f;
        logger->info(TDP, 2078,
          "[P2]   GATE cell={} arc={}->{} rf={} in_slew={:.6f}ns "
          "load_cap={:.4f}fF delay={:.6f}ns out_slew={:.6f}ns",
          cell_name, from_port, to_port, rf_s, in_slew,
          load_cap * 1e15, d * 1e9, out_slew);
        gate_tot += d;
      }
    }
    const double tot = gate_tot + wire_tot;
    logger->info(TDP, 2079,
      "[P2] path#{} totals: gate={:.6f}ns wire={:.6f}ns wire_frac={:.1f}%",
      i, gate_tot * 1e9, wire_tot * 1e9,
      tot != 0.0 ? 100.0 * wire_tot / tot : 0.0);
  }
}

void Timingdp::reportTimingStats(const int max_displacement_x,
                                 const int max_displacement_y,
                                 double margin_um,
                                 const std::string& tech,
                                 const std::string& hbt)
{
  profile_ = parseProfile(tech);
  importDb();

    if (max_displacement_x == 0 || max_displacement_y == 0) {
    // defaults
    max_displacement_x_ = 500;
    max_displacement_y_ = 100;
    } else {
    max_displacement_x_ = max_displacement_x;
    max_displacement_y_ = max_displacement_y;
    }

    logger_->info(TDP, 1012, "import right");

    std::string hbt_array_json = hbt;
    if (!initializeHBTAllocator(hbt_array_json, block_->getDbUnitsPerMicron())) {
        std::cerr << "WARNING: Failed to initialize HBT allocator, but continuing..." << std::endl;
    } else {
        std::cout << "✓ HBT allocator initialized successfully" << std::endl;

        // 执行HBT分配
        auto allocation_result = allocateHBTsForDesign();

        std::cout << "   - Total inter-layer nets: " << allocation_result.total_nets << std::endl;
        std::cout << "   - Successful allocations: " << allocation_result.successful_allocations << std::endl;
        std::cout << "   - Failed allocations: " << allocation_result.failed_allocations << std::endl;

        if (allocation_result.total_nets > 0) {
            std::cout << "   - Success rate: "
                        << (static_cast<double>(allocation_result.successful_allocations) / allocation_result.total_nets * 100.0)
                        << "%" << std::endl;
            std::cout << "   - Total estimated wirelength: " << allocation_result.total_wirelength << " microns" << std::endl;
        }
    }
    rs_->resizeSlackPreamble();

    // 先初始化 GRT 的 grid 和 db_net_map_
    int min_layer, max_layer;
    global_router_->getMinMaxLayer(min_layer, max_layer);
    global_router_->initFastRoute(min_layer, max_layer, true, 0);
    SteinerTreeGenerator steiner_generator(this);
    steiner_generator.createSteninerRoute();
    executeTimingDriven(true);
    auto stats = timing_engine_->getTimingStats();

    logger_->info(TDP, 1010, "Early: WNS={:.6f}ns TNS={:.6f}ns",
                  stats.early_wns * 1e9, stats.early_tns * 1e9);
    logger_->info(TDP, 1011, "Late : WNS={:.6f}ns TNS={:.6f}ns",
                  stats.late_wns * 1e9, stats.late_tns * 1e9);

    // P2 diagnostic: worst-N setup path arc-by-arc gate/wire dump (env TDP_DIAG_P2).
    diagWorstSetupPaths(timing_engine_.get(), logger_);
}

void Timingdp::legalizeTier(const std::string& tier,
                            double margin_um,
                            int max_displacement_x,
                            int max_displacement_y)
{
    // 按本次 -tech 决定的每层工艺选取 core site。
    Tech tier_tech = (tier == "upper") ? profile_.upper : profile_.bottom;
    std::string site_name = techSiteName(tier_tech);
    odb::dbSite* target_site = nullptr;
    for (auto* lib : db_->getLibs()) {
        target_site = lib->findSite(site_name.c_str());
        if (target_site) break;
    }
    if (!target_site) {
        logger_->error(TDP, 700, "Cannot find site '{}'", site_name);
    }

    for (auto* row : block_->getRows()) {
        odb::dbRow::destroy(row);
    }

    int margin_dbu = block_->micronsToDbu(margin_um);
    odb::Rect die = block_->getDieArea();
    odb::Rect core(die.xMin() + margin_dbu,
                   die.yMin() + margin_dbu,
                   die.xMax() - margin_dbu,
                   die.yMax() - margin_dbu);

    sta::dbNetwork* network = sta_->getDbNetwork();
    ifp::InitFloorplan floorplan(block_, logger_, network);
    floorplan.makeRows(core, target_site);

    std::string other_tier = (tier == "bottom") ? "upper" : "bottom";
    std::vector<std::pair<odb::dbMaster*, odb::dbMasterType>> saved_types;
    for (auto* lib : db_->getLibs()) {
        for (auto* master : lib->getMasters()) {
            std::string name = master->getName();
            if (name.find(other_tier) != std::string::npos) {
                saved_types.emplace_back(master, master->getType());
                master->setType(odb::dbMasterType::COVER);
            }
        }
    }

    opendp_->detailedPlacement(0, 0, "", false);

    for (auto& [master, type] : saved_types) {
        master->setType(type);
    }
}

bool Timingdp::destroyDatabase() {
    importClear();  // 清除转换数据（类似DPL的importClear）
    
    // if (db_) {
    //     odb::dbDatabase::destroy(db_);
    //     db_ = nullptr;
    //     chip_ = nullptr;
    //     block_ = nullptr;
    // }
    
    //logger_->info("Destroyed 3D placement database");
    return true;
}

// ============================================================================
// 3D PDK参数提取（仿照DPL从数据库获取）
// ============================================================================

bool Timingdp::extractPDKParametersFromDatabase(const std::string& layer) {
    if (!block_) {
        //logger_->error("No block available for PDK parameter extraction");
        return false;
    }
    
    try {
        logger_->info(TDP, 13, "Extracting PDK parameters from database...");

        if (!layer.empty()) {
            // createnewrow 阶段：此时 OpenDB 仅含该层的 row（newrow 逐层
            // "建一层 site -> 抓一层"）。层归属由调用方显式给定，避免同构
            // （上下层同名 site）下的误判。
            Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(layer);

            // 清空现有配置并仅填充当前层。
            pdk_config_.clear();
            for (odb::dbRow* row : block_->getRows()) {
                if (!row->getSite()) {
                    continue;
                }
                odb::dbSite* site = row->getSite();
                DBU height = site->getHeight();

                bool exists = false;
                for (auto& lc : pdk_config_.layers) {
                    if (lc.height == height) { lc.row_count++; exists = true; break; }
                }
                if (!exists) {
                    PDKLayerConfig layer_config;
                    layer_config.layer_id =
                        (layer_type == Architecture3D::LAYER_BOTTOM) ? 0 : 1;
                    layer_config.layer_name = layer + "_layer";
                    layer_config.height = height;
                    layer_config.site_width = site->getWidth();
                    layer_config.site_name = site->getName();
                    layer_config.tech_lib = site->getLib()->getName();
                    layer_config.row_count = 1;
                    pdk_config_.layers.push_back(layer_config);
                }
            }

            // 间距提取（不清空 layers）。
            extract3DSpacingParametersFromDatabase();

            if (arch_3d_) {
                arch_3d_->setHBTPitch(pdk_config_.hbt_pitch);
                arch_3d_->setInterLayerSpacing(pdk_config_.inter_layer_spacing);
                if (!pdk_config_.layers.empty()) {
                    arch_3d_->setSiteWidth(
                        layer_type, DbuX{pdk_config_.layers.front().site_width});
                    logger_->info(TDP, 8, "Configured layer " + layer +
                        " (height=" + std::to_string(pdk_config_.layers.front().height) +
                        ", width=" + std::to_string(pdk_config_.layers.front().site_width) +
                        ", tech_lib=" + pdk_config_.layers.front().tech_lib + ")");
                }
            }
        } else {
            // importDb 阶段：仅刷新层间/HBT 间距，不重新标层或改 site 宽度，
            // 以免覆盖 createnewrow 阶段已为各层设置好的参数。
            extract3DSpacingParametersFromDatabase();
            if (arch_3d_) {
                arch_3d_->setHBTPitch(pdk_config_.hbt_pitch);
                arch_3d_->setInterLayerSpacing(pdk_config_.inter_layer_spacing);
            }
        }

        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to extract PDK parameters: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 数据访问
// ============================================================================

int64_t Timingdp::getDieWidth() const {
    if (!block_) return 0;
    
    odb::Rect die_area = block_->getDieArea();
    return die_area.xMax() - die_area.xMin();
}

int64_t Timingdp::getDieHeight() const {
    if (!block_) return 0;
    
    odb::Rect die_area = block_->getDieArea();
    return die_area.yMax() - die_area.yMin();
}

std::vector<odb::dbInst*>& Timingdp::getMovableInstances() {
    return movable_insts_;
}

std::vector<odb::dbInst*>& Timingdp::getFixedInstances() {
    return fixed_insts_;
}

std::vector<odb::dbInst*>& Timingdp::getAllInstances() {
    return all_insts_;
}

std::vector<odb::dbNet*>& Timingdp::getNets() {
    return nets_;
}

// ============================================================================
// 3D特定操作
// ============================================================================


bool Timingdp::setInstanceLocation3D(odb::dbInst* inst, int64_t x, int64_t y, int64_t z) {
    if (!inst) return false;
    
    // 设置XY位置（Z位置目前只用于内部表示）
    inst->setLocation(x, y);
    
    // 更新3D网络中的位置
    if (network_3d_) {
        auto* node = network_3d_->getCell(inst);
        if (node) {
            // 确保Z坐标是0或1（bottom层=0, upper层=1）
            DBU clamped_z = (z > 0) ? 1 : 0;
            node->location = Pointdp3D(DbuX{x}, DbuY{y}, DbuZ{clamped_z});
        }
    }
    
    return true;
}


// ============================================================================
// 日志管理
// ============================================================================

utl::Logger* Timingdp::getLogger() const {
    return logger_;
}
// bool Timingdp::addHBTConnection(const std::string& bottom_inst, const std::string& upper_inst) {
//     // TODO: 实现HBT连接添加
//     logger_->warn("HBT connection addition not yet implemented");
//     return true;
// }

// std::vector<std::pair<std::string, std::string>> Timingdp::getHBTConnections() const {
//     // TODO: 实现HBT连接获取
//     return {};
// }

// ============================================================================
// 验证和调试
// ============================================================================

bool Timingdp::check3DPlacementLegality() const {
    if (!arch_3d_ || !network_3d_) {
        return false;
    }
    
    return network_3d_->checkPlacementLegality(arch_3d_.get());
}

/////

bool Timingdp::setupLayerArchitecture(const std::string& layer) {
    if (!arch_3d_ || !block_) {
        return false;
    }
    
    // 提取该层 PDK 参数（site 宽度等）；层归属由传入 layer 决定。
    if (!extractPDKParametersFromDatabase(layer)) {
        // logger_->error("Failed to extract PDK parameters for layer architecture");
        return false;
    }
    
    // 添加当前层（若尚未添加）；site 宽度已在 extractPDK 内设置。
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(layer);
    if (!arch_3d_->hasLayer(layer_type)) {
        arch_3d_->addLayer(layer_type, layer);
    }
    
    // 防御：确保两层都存在（另一层在另一次 createnewrow 调用中配置）。
    if (!arch_3d_->hasLayer(Architecture3D::LAYER_BOTTOM)) {
        arch_3d_->addLayer(Architecture3D::LAYER_BOTTOM, "bottom");
    }
    if (!arch_3d_->hasLayer(Architecture3D::LAYER_UPPER)) {
        arch_3d_->addLayer(Architecture3D::LAYER_UPPER, "upper");
    }
    
    return true;
}

void Timingdp::cacheDesignData() {
    logger_ ->info (TDP, 14, "Importing database to cacheDesignDat ......");
    cacheInstances();
    cacheNets();
}

void Timingdp::cacheInstances() {
    if (!block_) return;
    
    movable_insts_.clear();
    fixed_insts_.clear();
    all_insts_.clear();
    
    for (auto* inst : block_->getInsts()) {
        all_insts_.push_back(inst);
        
        if (inst->isFixed()) {
            fixed_insts_.push_back(inst);
        } else {
            movable_insts_.push_back(inst);
        }
    }
}

void Timingdp::cacheNets() {
    if (!block_) return;
    
    nets_.clear();
    
    for (auto* net : block_->getNets()) {
        nets_.push_back(net);
    }
}

void Timingdp::cleanup() {
    movable_insts_.clear();
    fixed_insts_.clear();
    all_insts_.clear();
    nets_.clear();
}

// ============================================================================
// 3D间距参数提取（仿照DPL从LEF文件获取）
// ============================================================================

void Timingdp::extract3DSpacingParametersFromDatabase() {
    if (!db_) {
        //logger_->warn("Database or chip not available for 3D spacing extraction");
        return;
    }
    
    try {
        //logger_->info("Extracting 3D spacing parameters from database...");
        
        // 重置间距参数
        pdk_config_.hbt_pitch = hbt_pitch_;  // 使用用户设置的HBT间距
        pdk_config_.inter_layer_spacing = inter_layer_spacing_;  // 使用用户设置的层间间距
        
        // 获取所有技术库
        std::vector<odb::dbLib*> libs;
        for (auto lib : db_->getLibs()) {
            libs.push_back(lib);
        }
        
        // 查找hb_layer层（从LEF文件中）- 仅用于信息记录，不覆盖用户设置
        bool found_hb_layer = false;
        for (odb::dbLib* lib : libs) {
            odb::dbTech* tech = lib->getTech();
            if (!tech) continue;
            
            odb::dbTechLayer* hb_layer = tech->findLayer("hb_layer");
            if (hb_layer) {
                DBU lef_hbt_pitch = hb_layer->getSpacing();
                found_hb_layer = true;
                //logger_->info("Found hb_layer in LEF: spacing=" + std::to_string(lef_hbt_pitch) + " DBU");
                
                // 如果用户没有设置HBT间距，使用LEF文件中的值
                if (hbt_pitch_ == 0) {
                    pdk_config_.hbt_pitch = lef_hbt_pitch;
                    hbt_pitch_ = lef_hbt_pitch;  // 更新用户变量
                    //logger_->info("Using HBT pitch from LEF: " + std::to_string(lef_hbt_pitch) + " DBU");
                } else {
                    //logger_->info("Using user-specified HBT pitch: " + std::to_string(hbt_pitch_) + " DBU (LEF value: " + std::to_string(lef_hbt_pitch) + " DBU)");
                }
                break;
            }
        }
        
        // 如果没有找到hb_layer，使用用户设置的值
        if (!found_hb_layer) {
            if (hbt_pitch_ == 0) {
                //logger_->warn("hb_layer not found in LEF files and no user HBT pitch specified. Using default 0.");
            } else {
                //logger_->info("Using user-specified HBT pitch: " + std::to_string(hbt_pitch_) + " DBU");
            }
        }
        
        // logger_->info("3D spacing parameters: HBT pitch=" + std::to_string(pdk_config_.hbt_pitch) +
        //               " DBU, inter-layer spacing=" + std::to_string(pdk_config_.inter_layer_spacing) + " DBU");
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to extract 3D spacing parameters: " + std::string(e.what()));
    }
}

// ============================================================================
// HBT分配接口实现（新增）
// ============================================================================
bool Timingdp::initializeHBTAllocator(const std::string& hbt_array_json, int dbu_per_micron) {
    //logger_->info("Initializing HBT allocator with JSON: " + hbt_array_json);

    // 读取 HBT array，并创建 HBTA allocator。
    hbt_interface_ = std::make_unique<hbta::HBTDatabaseInterface>();
    
    // 这里只初始化接口；真正分配在 allocateHBTsForDesign() 中触发。
    if (!hbt_interface_->initialize(this, hbt_array_json, dbu_per_micron)) {
        //logger_->error("Failed to initialize HBT database interface");
        hbt_interface_.reset();
        return false;
    }
    
    //logger_->info("HBT allocator initialized successfully");
    return true;
}

std::vector<hbta::InterLayerConnection> Timingdp::extractInterLayerConnections() {

    if (!hbt_interface_) {
        return {};
    }
    
    // 将需要 HBT 的跨层 net 转成 HBTA 使用的连接描述。
    auto connections = hbt_interface_->extractInterLayerConnections();
    
    // 打印提取统计，便于检查跨层 net 是否正确识别。
    auto stats = hbt_interface_->getExtractionStats();
    stats.print();
    
    return connections;
}

hbta::HBTAllocator* Timingdp::getHBTAllocator() {
    if (!hbt_interface_) {
        return nullptr;
    }

    // allocator 的所有权仍在 hbt_interface_ 中。
    return &hbt_interface_->getAllocator();
}

const hbta::HBTAllocator* Timingdp::getHBTAllocator() const {
    if (!hbt_interface_) {
        return nullptr;
    }

    // const 版本用于只读访问 allocator 状态。
    return &hbt_interface_->getAllocator();
}

hbta::HBTAllocator::BatchNetAllocationResult Timingdp::allocateHBTsForDesign() {
    hbta::HBTAllocator::BatchNetAllocationResult result;
    
    // 获取HBT分配器
    auto* hbt_allocator = getHBTAllocator();
    if (!hbt_allocator) {
        //logger_->error("HBT allocator not available");
        result.total_nets = 0;
        result.failed_allocations = 0;
        return result;
    }
    
    // 提取跨层 net，并缓存给后续 timing-driven 重匹配复用。
    auto connections = extractInterLayerConnections();
    last_hbt_connections_ = connections;
    
    if (connections.empty()) {
        result.total_nets = 0;
        result.successful_allocations = 0;
        result.failed_allocations = 0;
        return result;
    }

    // 初始分配：HBTA 内部使用四叉树分区 + 局部 Hungarian。
    result = hbt_allocator->allocateHBTsForNets(connections);

    // 同步回 Network3D，后续 routing/STA 才能看到 HBT 连接。
    if (result.successful_allocations > 0) {
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
    }

    return result;
}

bool Timingdp::runTimingDrivenHBTReallocation(
    const std::vector<hbta::InterLayerConnection>& connections) {
    if (!hbt_interface_ || connections.empty()) {
        return false;
    }

    return hbt_interface_->runTimingDrivenReallocation(connections);
}

/////////////////////
///////Segment///////
////////////////////

void Timingdp::initializeSegmentManager() {
    if (arch_3d_) {
        // 动态分配段管理器
        if (segment_manager_) {
            delete segment_manager_;  // 清理旧的（如果存在）
        }
        logger_->info(TDP, 15, "Segment manager initialized");
        segment_manager_ = new SegmentManager3D(arch_3d_.get());
        segment_manager_->initializeSegments();
        //logger_->info("Segment manager initialized");
    } else {
        //logger_->error("Cannot initialize segment manager: Architecture3D not available");
    }
}

// 实现从 .h 文件移出的方法
std::vector<Segment3D*> Timingdp::getCellSegments(const Node3D* cell) const {
    return segment_manager_ ? segment_manager_->getCellSegments(cell) : std::vector<Segment3D*>();
}

bool Timingdp::addCellToSegment(Node3D* cell, int seg_id) {
    return segment_manager_ ? segment_manager_->addCellToSegment(cell, seg_id) : false;
}

bool Timingdp::removeCellFromSegment(const Node3D* cell, int seg_id) {
    return segment_manager_ ? segment_manager_->removeCellFromSegment(cell, seg_id) : false;
}

void Timingdp::initializeSegments() {
    if (segment_manager_) {
        segment_manager_->initializeSegments();
    }
}

void Timingdp::createnewrow(const std::string& layer){
    //设置层架构（仿照DPL的createArchitecture）
    if (!setupLayerArchitecture(layer)) {
        logger_->error(TDP, 591, "Failed to setup layer architecture");
    }
    
    //仿照DPL的dbToOpendp::createArchitecture）
    if (!convertRows(layer)) {
        logger_->error(TDP, 592, "Failed to convert rows");
    }

}

void Timingdp::reportLegalizationStats() const
{
  logger_->report("Placement Analysis");
  logger_->report("---------------------------------");
//   logger_->report("total displacement   {:10.1f} u",
//                   block_->dbuToMicrons(displacement_sum_));
//   logger_->metric("design__instance__displacement__total",
//                   block_->dbuToMicrons(displacement_sum_));
//   logger_->report("average displacement {:10.1f} u",
//                   block_->dbuToMicrons(displacement_avg_));
//   logger_->metric("design__instance__displacement__mean",
//                   block_->dbuToMicrons(displacement_avg_));
//   logger_->report("max displacement     {:10.1f} u",
//                   block_->dbuToMicrons(displacement_max_));
//   logger_->metric("design__instance__displacement__max",
//                   block_->dbuToMicrons(displacement_max_));
//   logger_->report("original HPWL        {:10.1f} u",
//                   block_->dbuToMicrons(hpwl_before_));
//   odb::WireLengthEvaluator eval(block_);
//   const double hpwl_legal = eval.hpwl();
//   logger_->report("legalized HPWL       {:10.1f} u",
//                   block_->dbuToMicrons(hpwl_legal));
//   logger_->metric("route__wirelength__estimated",
//                   block_->dbuToMicrons(hpwl_legal));
//   const int hpwl_delta
//       = (hpwl_before_ == 0.0)
//             ? 0.0
//             : round((hpwl_legal - hpwl_before_) / hpwl_before_ * 100);
//   logger_->report("delta HPWL           {:10} %", hpwl_delta);
  logger_->report("");
}


//////////////
/////test sta
/////////////
bool Timingdp::executeTimingDriven(bool)
{
    sta_->ensureGraph();
    // 更新时序
    sta_->updateTiming(true);
    sta_->findRequireds();

    return true;
}
void Timingdp::executeTimingDrivenI(Node3D* cell, int depth)
{
    if (!cell || !cell->db_inst) return;
    sta::Instance* sta_inst = sta_->getDbNetwork()->dbToSta(cell->db_inst);
    if (!sta_inst) return;

    if (!ista_)
        ista_ = std::make_unique<lgista::LocalGraphIsta>(sta_);

    static bool compare_env = (std::getenv("LGISTA_COMPARE") != nullptr);
    if (compare_env) {
        struct PinSlew { float rise; float fall; };
        std::map<std::string, PinSlew> old_slew;
        sta::DcalcAPIndex ap0 = 0;
        {
            sta_->updateTiming(false);
            sta::InstancePinIterator* pit = sta_->network()->pinIterator(sta_inst);
            while (pit->hasNext()) {
                sta::Pin* pin = pit->next();
                sta::Vertex *v, *bidir;
                sta_->graph()->pinVertices(pin, v, bidir);
                if (v) {
                    std::string name = sta_->network()->name(pin);
                    old_slew[name].rise = sta::delayAsFloat(
                        sta_->graph()->slew(v, sta::RiseFall::rise(), ap0));
                    old_slew[name].fall = sta::delayAsFloat(
                        sta_->graph()->slew(v, sta::RiseFall::fall(), ap0));
                }
            }
            delete pit;
        }

        ista_->analyzeAndCommit(sta_inst, depth);

        double max_slew_delta = 0;
        int pin_count = 0;
        sta::InstancePinIterator* pit2 = sta_->network()->pinIterator(sta_inst);
        while (pit2->hasNext()) {
            sta::Pin* pin = pit2->next();
            sta::Vertex *v, *bidir;
            sta_->graph()->pinVertices(pin, v, bidir);
            if (!v) continue;
            std::string name = sta_->network()->name(pin);
            auto it = old_slew.find(name);
            if (it == old_slew.end()) continue;
            double new_r = sta::delayAsFloat(sta_->graph()->slew(v, sta::RiseFall::rise(), ap0));
            double new_f = sta::delayAsFloat(sta_->graph()->slew(v, sta::RiseFall::fall(), ap0));
            max_slew_delta = std::max(max_slew_delta,
                std::max(std::abs(new_r - (double)it->second.rise),
                         std::abs(new_f - (double)it->second.fall)));
            ++pin_count;
        }
        delete pit2;
        printf("[LGISTA_COMPARE] cell=%s pins=%d max_slew_delta=%.4e\n",
               sta_->network()->name(sta_inst), pin_count, max_slew_delta);
        fflush(stdout);
        return;
    }

    ista_->analyzeAndCommit(sta_inst, depth);
}

lgista::PtGraph* Timingdp::analyzeOnly(Node3D* cell, int depth)
{
    if (!cell || !cell->db_inst) return nullptr;
    sta::Instance* sta_inst = sta_->getDbNetwork()->dbToSta(cell->db_inst);
    if (!sta_inst) return nullptr;

    if (!ista_)
        ista_ = std::make_unique<lgista::LocalGraphIsta>(sta_);

    // Transfer ownership of the local PtGraph to the caller.
    return ista_->analyzeOnly(sta_inst, depth).release();
}

bool Timingdp::isIdealClockNet(odb::dbNet* db_net) const {
    if (!db_net || !sta_) {
        return false;
    }
    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Net* sta_net = db_network->dbToSta(db_net);
    if (!sta_net) {
        return false;
    }
    sta::Network* network = sta_->network();
    sta::PinSet* drivers = network->drivers(sta_net);
    if (drivers && !drivers->empty()) {
        const sta::Pin* drvr_pin = *drivers->begin();
        for (sta::Mode* mode : sta_->modes()) {
            if (!sta_->isIdealClock(drvr_pin, mode)) {
                return false;
            }
        }
    }
    return true;
}

bool Timingdp::isValidSigType(odb::dbNet* db_net) const
{
    if (!db_net)
        return false;
    auto sig = db_net->getSigType();
    return (sig == odb::dbSigType::SIGNAL
            || sig == odb::dbSigType::CLOCK);
}


} // namespace tdp
