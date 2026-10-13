// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "HBTDatabaseInterface.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "HBTACostModel.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "eva/TimingEngine3D.h"
#include "eva/eva.h"
#include "PtGraph.hh"
#include "search/Tag.hh"
#include "search/TagGroup.hh"
#include "sta/Graph.hh"
#include "sta/MinMax.hh"
#include "sta/Network.hh"
#include "sta/Path.hh"
#include "utl/Logger.h"

// 包含必要的3D-TDP头文件以获取完整类型定义
#include "tdp/Timingdp.h"
#include "infrastructure/Coordinates.h"
#include "infrastructure/Network3D.h"
#include "infrastructure/Net3D.h"
#include "infrastructure/Node3D.h"

namespace {

bool ptVertexLocalLateSlack(lgista::PtGraph* pg,
                            lgista::PtVertex& ptv,
                            sta::StaState* sta_state,
                            double& late_slack_out)
{
    sta::Path* paths = ptv.paths();
    if (paths == nullptr) {
        return false;
    }

    sta::TagGroup* tg = pg->tagGroup(ptv);
    if (tg == nullptr) {
        return false;
    }

    const size_t count = tg->pathCount();
    bool have_late = false;
    float worst_late_slack = 0.0f;

    for (size_t i = 0; i < count; ++i) {
        sta::Tag* tag = paths[i].tag(sta_state);
        if (tag == nullptr || tag->isGenClkSrcPath()) {
            continue;
        }
        if (tag->minMax() != sta::MinMax::max()) {
            continue;
        }

        const float slack = sta::delayAsFloat(paths[i].slack(sta_state));
        if (!have_late || slack < worst_late_slack) {
            worst_late_slack = slack;
            have_late = true;
        }
    }

    if (!have_late) {
        return false;
    }

    late_slack_out = worst_late_slack;
    return true;
}

}  // namespace

namespace hbta {

// ============================================================================
// HBTDatabaseInterface 实现
// ============================================================================

HBTDatabaseInterface::HBTDatabaseInterface() {
    // 初始化统计信息
    extraction_stats_ = ExtractionStats();
}

bool HBTDatabaseInterface::initialize(tdp::Timingdp* timingdp,
                                     const std::string& hbt_array_json, 
                                     int dbu_per_micron) {
    if (!timingdp) {
        std::cerr << "ERROR: Timingdp pointer is null" << std::endl;
        return false;
    }
    
    timingdp_ = timingdp;
    hbt_array_json_path_ = hbt_array_json;
    dbu_per_micron_ = dbu_per_micron;
    
    // 初始化HBT分配器
    if (!allocator_.initialize(hbt_array_json, dbu_per_micron)) {
        std::cerr << "ERROR: Failed to initialize HBT allocator" << std::endl;
        return false;
    }
    
    // 重置统计信息
    extraction_stats_ = ExtractionStats();
    cell_layer_cache_.clear();
    
    std::cout << "HBTDatabaseInterface initialized successfully" << std::endl;
    std::cout << "  - HBT array JSON: " << hbt_array_json << std::endl;
    std::cout << "  - DBU per micron: " << dbu_per_micron << std::endl;
    
    return true;
}

// 从 Network3D 提取所有跨层信号网络，并转换成 HBTA 输入。
std::vector<InterLayerConnection> HBTDatabaseInterface::extractInterLayerConnections() {
    std::vector<InterLayerConnection> connections;
    
    if (!timingdp_) {
        std::cerr << "ERROR: Timingdp not initialized" << std::endl;
        return connections;
    }
    
    // 获取3D网络
    auto* network_3d = timingdp_->getNetwork3D();
    if (!network_3d) {
        std::cerr << "ERROR: Network3D not available" << std::endl;
        return connections;
    }
    
    // 重置统计信息
    extraction_stats_ = ExtractionStats();
    
    // 获取数据库中的所有网络
    auto& db_nets = timingdp_->getNets();
    extraction_stats_.total_nets = db_nets.size();
    
    std::cout << "Extracting inter-layer connections from " << db_nets.size() << " nets..." << std::endl;
    
    // 遍历所有数据库网络
    for (auto* db_net : db_nets) {
        if (!db_net) continue;
        
        // 跳过电源网络
        if (db_net->getSigType().isSupply()) {
            continue;
        }
        
        // 跳过时钟网络
        if (timingdp_->isIdealClockNet(db_net)) {
            continue;
        }

        // 获取对应的Net3D对象
        auto* net_3d = network_3d->getNet(db_net);
        if (!net_3d) {
            // 如果没有对应的Net3D，跳过
            continue;
        }
        
        // 检查是否是跨层网络
        if (isInterLayerNetwork(net_3d)) {
            extraction_stats_.inter_layer_nets++;
            
            // 转换为InterLayerConnection
            auto connection = convertNetToConnection(net_3d);
            if (!connection.net_name.empty()) {
                connections.push_back(connection);
            } else {
                extraction_stats_.failed_conversions++;
            }
        }
    }
    
    std::cout << "Found " << connections.size() << " inter-layer connections" << std::endl;
    return connections;
}

// 将单个 Net3D 转成 InterLayerConnection，补齐 driver 和各层 load 信息。
InterLayerConnection HBTDatabaseInterface::convertNetToConnection(tdp::Net3D* net) {
    InterLayerConnection connection;
    
    if (!net || !net->db_net) {
        return connection;
    }
    
    connection.net_name = net->name;
    
    // 通过 STA 获取 driver pin
    auto* sta = timingdp_->getSta();
    auto* db_network = sta->getDbNetwork();
    auto* network = sta->network();
    
    sta::Net* sta_net = db_network->dbToSta(net->db_net);
    sta::PinSet* drivers = network->drivers(sta_net);
    if (!drivers || drivers->empty()) {
        return connection;
    }
    const sta::Pin* drvr_pin = *drivers->begin();
    
    // STA driver pin → db pin
    odb::dbITerm* drvr_iterm = nullptr;
    odb::dbBTerm* drvr_bterm = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    db_network->staToDb(drvr_pin, drvr_iterm, drvr_bterm, moditerm);
    
    // 通过 Network3D 映射 O(1) 查找 driver cell
    auto* network_3d = timingdp_->getNetwork3D();
    tdp::Node3D* driver_cell = nullptr;
    if (drvr_iterm) {
        driver_cell = network_3d->getCell(drvr_iterm->getInst());
    } else if (drvr_bterm) {
        driver_cell = network_3d->getNode(drvr_bterm);
    }
    
    if (!driver_cell) {
        return connection;
    }
    
    // driver 信息（使用真实 layer）
    connection.driver_cell_name = driver_cell->name;
    connection.driver_location = convertToHBTAPoint(driver_cell->location);
    connection.driver_layer = driver_cell->layer;
    
    // 其余 cell 均为 load
    for (auto* cell : net->connected_cells) {
        if (!cell || cell == driver_cell) continue;
        std::string layer = getCellLayer(cell);
        if (layer == "bottom") extraction_stats_.bottom_layer_cells++;
        else if (layer == "upper") extraction_stats_.upper_layer_cells++;
        connection.load_cell_names.push_back(cell->name);
        connection.load_locations.push_back(convertToHBTAPoint(cell->location));
        connection.load_layers.push_back(layer);
    }
    
    // 计算建议的HBT位置
    connection.suggested_hbt_location = calculateSuggestedHBTLocation(net);
    
    return connection;
}

Point HBTDatabaseInterface::calculateSuggestedHBTLocation(tdp::Net3D* net) {
    
    return Point(0.0, 0.0);
}

// 收集重分配前的 baseline slack，供候选筛选和打分使用。
std::unordered_map<std::string, double> HBTDatabaseInterface::collectNetSlacks(
    const std::vector<InterLayerConnection>& connections) const {
    std::unordered_map<std::string, double> slack_by_net;
    if (timingdp_ == nullptr) {
        return slack_by_net;
    }

    auto* timing_engine = timingdp_->gettimingengine();
    if (timing_engine == nullptr) {
        return slack_by_net;
    }

    std::unordered_map<std::string, odb::dbNet*> db_net_by_name;
    auto& nets = timingdp_->getNets();
    db_net_by_name.reserve(nets.size());
    for (auto* db_net : nets) {
        if (db_net != nullptr) {
            db_net_by_name.emplace(db_net->getName(), db_net);
        }
    }

    for (const auto& connection : connections) {
        auto net_it = db_net_by_name.find(connection.net_name);
        if (net_it == db_net_by_name.end()) {
            continue;
        }

        const double slack_value = static_cast<double>(
            timing_engine->getNetSlack(connection.net_name, sta::MinMax::max()));
        if (!std::isfinite(slack_value)) {
            continue;
        }

        slack_by_net[connection.net_name] = slack_value;

        auto* network_3d = timingdp_->getNetwork3D();
        if (network_3d != nullptr) {
            if (auto* net_3d = network_3d->getNet(connection.net_name)) {
                net_3d->updateTimingInfo(slack_value);
            }
        }
    }

    return slack_by_net;
}

std::optional<double> HBTDatabaseInterface::collectNetSlack(
    const InterLayerConnection& connection) const {
    if (timingdp_ == nullptr) {
        return std::nullopt;
    }

    auto* timing_engine = timingdp_->gettimingengine();
    if (timing_engine == nullptr) {
        return std::nullopt;
    }

    const double slack_value = static_cast<double>(
        timing_engine->getNetSlack(connection.net_name, sta::MinMax::max()));

    if (!std::isfinite(slack_value)) {
        return std::nullopt;
    }

    auto* network_3d = timingdp_->getNetwork3D();
    if (network_3d != nullptr) {
        if (auto* net_3d = network_3d->getNet(connection.net_name)) {
            net_3d->updateTimingInfo(slack_value);
        }
    }

    return slack_value;
}

std::optional<double> HBTDatabaseInterface::collectNetSlackFromPtGraph(
    const InterLayerConnection& connection,
    lgista::PtGraph* pt_graph) const
{
    if (timingdp_ == nullptr || pt_graph == nullptr) {
        return std::nullopt;
    }

    auto* sta = timingdp_->getSta();
    auto* db_network = sta ? sta->getDbNetwork() : nullptr;
    auto* network = sta ? sta->network() : nullptr;
    auto* graph = sta ? sta->graph() : nullptr;
    if (sta == nullptr || db_network == nullptr || network == nullptr || graph == nullptr) {
        return std::nullopt;
    }

    sta::Net* net = db_network->findNet(connection.net_name.c_str());
    if (net == nullptr) {
        return std::nullopt;
    }

    std::unique_ptr<sta::NetConnectedPinIterator> pin_iter(
        network->connectedPinIterator(net));
    auto* timing_engine = timingdp_->gettimingengine();

    bool saw_load = false;
    bool have_slack = false;
    double worst_net_slack = 0.0;

    while (pin_iter->hasNext()) {
        const sta::Pin* pin = pin_iter->next();
        if (!network->isLoad(pin)) {
            continue;
        }
        saw_load = true;

        double pin_late_slack = 0.0;
        bool have_pin_slack = false;

        sta::Vertex* load_vertex = graph->pinLoadVertex(pin);
        if (load_vertex != nullptr) {
            lgista::PtVertex* ptv = pt_graph->ptVertex(load_vertex);
            if (ptv != nullptr
                && ptVertexLocalLateSlack(pt_graph, *ptv, sta, pin_late_slack)
                && std::isfinite(pin_late_slack)) {
                have_pin_slack = true;
            }
        }

        // Fallback: sink outside the local PtGraph depth -> reuse frozen
        // global baseline late slack, mirroring LSP3D's PtGraph refresh path.
        if (!have_pin_slack && timing_engine != nullptr) {
            const auto timing_info = timing_engine->getPinTimingInfoDirect(
                const_cast<sta::Pin*>(pin));
            if (timing_info.is_valid && std::isfinite(timing_info.late_slack)) {
                pin_late_slack = timing_info.late_slack;
                have_pin_slack = true;
            }
        }

        if (!have_pin_slack) {
            continue;
        }

        if (!have_slack || pin_late_slack < worst_net_slack) {
            worst_net_slack = pin_late_slack;
            have_slack = true;
        }
    }

    if (!saw_load || !have_slack || !std::isfinite(worst_net_slack)) {
        return std::nullopt;
    }

    return worst_net_slack;
}

bool HBTDatabaseInterface::isInterLayerNetwork(tdp::Net3D* net) const {
    if (!net || net->connected_cells.empty()) {
        return false;
    }
    // 只要出现两种及以上层标签，就视为跨层网络。
    std::set<std::string> layers;
    for (auto* cell : net->connected_cells) {
        if (!cell) continue;
        layers.insert(getCellLayer(cell));
    }
    
    return layers.size() > 1;
}

std::string HBTDatabaseInterface::getCellLayer(tdp::Node3D* cell) const {
    if (!cell) {
        return "unknown";
    }
    
    // 先查缓存，避免在大网表中重复字符串处理。
    auto it = cell_layer_cache_.find(cell->name);
    if (it != cell_layer_cache_.end()) {
        return it->second;
    }
    
    // 直接从 Node3D 读取层标签。
    std::string layer = cell->layer;
    
    // 标准化层名称
    if (layer.empty()) {
        layer = "unknown";
    }
    
    // 更新缓存（使用const_cast，因为缓存应该是mutable的）
    // 在实际实现中，应该将cell_layer_cache_声明为mutable
    cell_layer_cache_[cell->name] = layer;
    
    return layer;
}

Point HBTDatabaseInterface::convertToHBTAPoint(const tdp::Pointdp3D& point) const {
    // 将DBU坐标转换为微米坐标
    // double x_um = static_cast<double>(point.x.v) / dbu_per_micron_;
    // double y_um = static_cast<double>(point.y.v) / dbu_per_micron_;
    // return Point(x_um, y_um);

    // 统一使用 DBU，避免 HBTArray(DBU) 与网络坐标(um) 混用。
    double x_dbu = static_cast<double>(point.x.v);
    double y_dbu = static_cast<double>(point.y.v);
    
    return Point(x_dbu, y_dbu);
}

tdp::Point HBTDatabaseInterface::convertToTDPPoint(const Point& point) const {
    return tdp::Point(tdp::DbuX(convertToTDPDBU(point.x)),
                      tdp::DbuY(convertToTDPDBU(point.y)));
}

int64_t HBTDatabaseInterface::convertToTDPDBU(double value) const {
    return static_cast<int64_t>(std::llround(value));
}

// 按 allocator 的当前分配结果，重建单个 net 在 Network3D 中的 HBT 状态。
void HBTDatabaseInterface::syncHBTAssignmentToNetwork3D(
    const InterLayerConnection& connection) {
    if (timingdp_ == nullptr) {
        return;
    }

    auto* network_3d = timingdp_->getNetwork3D();
    if (network_3d == nullptr) {
        return;
    }

    const tdp::DBU pitch = static_cast<tdp::DBU>(
        convertToTDPDBU(allocator_.getHBTArray().getPitchDBU()));

    auto* net_3d = network_3d->getNet(connection.net_name);
    if (net_3d == nullptr) {
        return;
    }

    // 先清空旧的 HBT 连接，避免多次同步时把历史状态残留在 Net3D 中。
    net_3d->hbt_connections.clear();

    const auto* hbt = allocator_.getHBTForNet(connection.net_name);
    if (hbt == nullptr) {
        return;
    }

    tdp::Net3D::HBTConnection hbt_conn;
    hbt_conn.location = convertToTDPPoint(hbt->location);
    hbt_conn.pitch = pitch;
    hbt_conn.aligned = false;

    // 根据 connection 中记录的层信息，把 driver/load 映射回 Net3D 需要的上下层端点。
    auto set_endpoint = [&](const std::string& cell_name,
                            const std::string& layer) {
        auto* cell = network_3d->getCell(cell_name);
        if (cell == nullptr) {
            return;
        }
        if (layer == "bottom" && hbt_conn.bottom_cell == nullptr) {
            hbt_conn.bottom_cell = cell;
        } else if (layer == "upper" && hbt_conn.upper_cell == nullptr) {
            hbt_conn.upper_cell = cell;
        }
    };

    set_endpoint(connection.driver_cell_name, connection.driver_layer);
    for (size_t i = 0; i < connection.load_cell_names.size(); ++i) {
        const std::string layer =
            i < connection.load_layers.size()
                ? connection.load_layers[i]
                : std::string();
        set_endpoint(connection.load_cell_names[i], layer);
    }

    net_3d->addHBTConnection(hbt_conn);
    net_3d->updateHBTAlignment();
}

// 批量同步所有跨层网的 HBT 分配结果，通常在一次试探、回滚或正式提交后调用。
void HBTDatabaseInterface::syncHBTAssignmentsToNetwork3D(
    const std::vector<InterLayerConnection>& connections) {
    for (const auto& connection : connections) {
        syncHBTAssignmentToNetwork3D(connection);
    }
}

// 通过增量 timing 评估局部 HBT 候选，并决定提交还是回滚本轮重分配。
bool HBTDatabaseInterface::runTimingDrivenReallocation(
    const std::vector<InterLayerConnection>& connections) {
    if (timingdp_ == nullptr || connections.empty()) {
        return false;
    }

    auto* allocator = &allocator_;
    auto* timing_engine = timingdp_->gettimingengine();
    if (timing_engine == nullptr) {
        return false;
    }

    auto* logger = timingdp_->getLogger();
    const auto old_stats = timing_engine->getTimingStats();

    auto slack_by_net = collectNetSlacks(connections);
    if (slack_by_net.empty()) {
        logger->warn(utl::TDP,
                     1070,
                     "Timing-driven HBT reallocation skipped: no STA slack was available for inter-layer nets.");
        return false;
    }

    // 当前策略：把所有跨层网都纳入 timing-driven 流程，但候选 HBT 只在局部窗口内搜索。
    HBTAllocator::TimingDrivenReallocationConfig config;
    config.critical_slack_threshold = 10e-12;
    config.critical_nets_percent = 100.0;
    config.local_bonding_box_radius_in_hbt_pitch = 5.0;

    if (allocator->hasPendingTimingDrivenReallocation()) {
        allocator->rollbackPendingTimingDrivenReallocation();
    }

    const auto problem = allocator->prepareTimingDrivenReallocationProblem(
        connections, slack_by_net, config);

    if (problem.empty()) {
        syncHBTAssignmentsToNetwork3D(connections);
        logger->info(utl::TDP,
                     1071,
                     "Timing-driven HBT reallocation produced no candidate.");
        return false;
    }

    // baseline 记录进入本轮优化前的分配状态，后续每个候选都基于它独立试探。
    const auto baseline_assignment = allocator->captureAssignmentSnapshot();
    const double failed_sta_cost = 1e20;
    std::vector<std::vector<double>> sta_cost_matrix(
        problem.connections.size(),
        std::vector<double>(problem.candidate_hbts.size(), failed_sta_cost));
    std::ofstream rematch_debug_log;
    if (const char* debug_log_path = std::getenv("HBTA_REMATCH_DEBUG_LOG")) {
        rematch_debug_log.open(debug_log_path, std::ios::app);
        rematch_debug_log << std::setprecision(8);
    }

    // 先把 slack 转成 criticality，后面候选打分时直接与 trial slack 结合。
    std::unordered_map<std::string, double> criticality_by_net;
    criticality_by_net.reserve(problem.connections.size());
    for (const auto& connection : problem.connections) {
        double criticality = 0.0;
        auto slack_it = slack_by_net.find(connection.net_name);
        if (slack_it != slack_by_net.end()) {
            criticality = HBTACostModel::criticalityFromSlack(
                slack_it->second,
                problem.worst_slack,
                config.critical_slack_threshold);
        }
        criticality_by_net[connection.net_name] = criticality;
    }

    tdp::SteinerTreeGenerator trial_steiner_generator(timingdp_);


    auto rebuild_connection_route_for_trial =
        [&](const InterLayerConnection& connection,
            tdp::Node3D*& driver_cell) {
            driver_cell = nullptr;
            auto* network_3d = timingdp_->getNetwork3D();
            if (network_3d == nullptr) {
                return false;
            }
            auto* net_3d = network_3d->getNet(connection.net_name);
            if (net_3d == nullptr) {
                return false;
            }
            if (!trial_steiner_generator.updateNetRouteAndParasitics(net_3d)) {
                return false;
            }

            if (!connection.driver_cell_name.empty()) {
                driver_cell = network_3d->getCell(connection.driver_cell_name);
            }
            return true;
        };

    // 某个候选试探失败后，恢复该 net 对应的 baseline 分配
    auto restore_connection_to_baseline =
        [&](const InterLayerConnection& connection) {
            allocator->restoreAssignmentSnapshot(baseline_assignment);
            syncHBTAssignmentToNetwork3D(connection);
            tdp::Node3D* driver_cell = nullptr;
            if (!rebuild_connection_route_for_trial(connection, driver_cell)) {
                return false;
            }
            if (driver_cell == nullptr) {
                return timingdp_->executeTimingDriven(true);
            }
            timingdp_->executeTimingDrivenI(driver_cell, 2);
            return true;
        };

    auto restore_allocator_to_baseline = [&]() {
        allocator->restoreAssignmentSnapshot(baseline_assignment);
    };

    auto restore_full_baseline = [&]() {
        allocator->restoreAssignmentSnapshot(baseline_assignment);
        syncHBTAssignmentsToNetwork3D(connections);
        tdp::SteinerTreeGenerator restore_steiner_generator(timingdp_);
        restore_steiner_generator.createSteninerRoute();
        return timingdp_->executeTimingDriven(true);
    };

    std::vector<size_t> active_row_indices;
    active_row_indices.reserve(problem.connections.size());
    for (size_t net_idx = 0; net_idx < problem.connections.size(); ++net_idx) {
        bool row_has_valid_sta_cost = false;
        bool row_restore_failed = false;
        const auto& candidate_connection = problem.connections[net_idx];
        const std::vector<size_t>* candidate_indices =
            net_idx < problem.candidate_indices_by_connection.size()
                ? &problem.candidate_indices_by_connection[net_idx]
                : nullptr;

        if (candidate_indices == nullptr || candidate_indices->empty()) {
            logger->warn(utl::TDP,
                         1077,
                         "Timing-driven HBT reallocation skipping net {}: no local HBT candidates.",
                         candidate_connection.net_name);
            continue;
        }

        if (rematch_debug_log) {
            rematch_debug_log << "[hbta-rematch-eval] step=row_begin net="
                              << candidate_connection.net_name
                              << " candidate_count=" << candidate_indices->size()
                              << " candidates=";
            bool first_candidate = true;
            for (size_t hbt_idx : *candidate_indices) {
                if (hbt_idx >= problem.candidate_hbts.size()) {
                    continue;
                }
                const auto& candidate_hbt = problem.candidate_hbts[hbt_idx];
                if (!first_candidate) {
                    rematch_debug_log << ';';
                }
                first_candidate = false;
                rematch_debug_log << candidate_hbt.hbt_id << "@("
                                  << candidate_hbt.location.x << ","
                                  << candidate_hbt.location.y << ")";
            }
            rematch_debug_log << '\n';
        }

        std::string current_hbt_id;
        for (const auto& [net_name, hbt_id] : problem.previous_assignments) {
            if (net_name == candidate_connection.net_name) {
                current_hbt_id = hbt_id;
                break;
            }
        }
        bool has_current_cost = false;
        double current_slack = 0.0;
        double current_cost = failed_sta_cost;
        bool has_best_cost = false;
        std::string best_hbt_id;
        double best_slack = 0.0;
        double best_cost = failed_sta_cost;

        for (size_t hbt_idx : *candidate_indices) {
            if (hbt_idx >= problem.candidate_hbts.size()) {
                continue;
            }
            restore_allocator_to_baseline();

            // 把当前 net 强制绑定到候选 HBT，再同步回 Network3D 做一次增量评估。
            const auto& candidate_hbt = problem.candidate_hbts[hbt_idx];
            if (!allocator->forceAssignHBTToNetForEvaluation(
                    candidate_connection, candidate_hbt.hbt_id)) {
                if (!restore_connection_to_baseline(candidate_connection)) {
                    row_restore_failed = true;
                    break;
                }
                continue;
            }

            syncHBTAssignmentToNetwork3D(candidate_connection);
            tdp::Node3D* driver_cell = nullptr;
            if (!rebuild_connection_route_for_trial(candidate_connection, driver_cell)) {
                if (!restore_connection_to_baseline(candidate_connection)) {
                    row_restore_failed = true;
                    break;
                }
                continue;
            }
            if (driver_cell == nullptr) {
                if (!restore_connection_to_baseline(candidate_connection)) {
                    row_restore_failed = true;
                    break;
                }
                continue;
            }

            std::unique_ptr<lgista::PtGraph> pt_graph(
                timingdp_->analyzeOnly(driver_cell, 2));
            if (!pt_graph) {
                if (!restore_connection_to_baseline(candidate_connection)) {
                    row_restore_failed = true;
                    break;
                }
                continue;
            }

            const auto trial_slack =
                collectNetSlackFromPtGraph(candidate_connection, pt_graph.get());
            if (!trial_slack.has_value() || !std::isfinite(*trial_slack)) {
                if (!restore_connection_to_baseline(candidate_connection)) {
                    row_restore_failed = true;
                    break;
                }
                continue;
            }

            auto criticality_it = criticality_by_net.find(candidate_connection.net_name);
            const double criticality = criticality_it == criticality_by_net.end()
                                           ? 0.0
                                           : criticality_it->second;
            sta_cost_matrix[net_idx][hbt_idx] = -criticality * (*trial_slack);
            if (rematch_debug_log) {
                rematch_debug_log << "[hbta-rematch-eval] step=trial_ok net="
                                  << candidate_connection.net_name
                                  << " hbt=" << candidate_hbt.hbt_id
                                  << " hbt_xy=(" << candidate_hbt.location.x
                                  << "," << candidate_hbt.location.y << ")"
                                  << " slack=" << *trial_slack
                                  << " criticality=" << criticality
                                  << " sta_cost=" << sta_cost_matrix[net_idx][hbt_idx]
                                  << '\n';
            }
            if (candidate_hbt.hbt_id == current_hbt_id) {
                has_current_cost = true;
                current_slack = *trial_slack;
                current_cost = sta_cost_matrix[net_idx][hbt_idx];
            }
            if (!has_best_cost || sta_cost_matrix[net_idx][hbt_idx] < best_cost) {
                has_best_cost = true;
                best_hbt_id = candidate_hbt.hbt_id;
                best_slack = *trial_slack;
                best_cost = sta_cost_matrix[net_idx][hbt_idx];
            }
            row_has_valid_sta_cost = true;

            if (!restore_connection_to_baseline(candidate_connection)) {
                row_restore_failed = true;
                break;
            }
        }

        if (rematch_debug_log) {
            rematch_debug_log << "[hbta-rematch-eval] step=row_result net="
                              << candidate_connection.net_name
                              << " has_valid_cost=" << row_has_valid_sta_cost
                              << " current_hbt="
                              << (current_hbt_id.empty() ? "none" : current_hbt_id)
                              << " current_slack="
                              << (has_current_cost ? current_slack : 0.0)
                              << " current_cost="
                              << (has_current_cost ? current_cost : failed_sta_cost)
                              << " best_hbt="
                              << (has_best_cost ? best_hbt_id : "none")
                              << " best_slack="
                              << (has_best_cost ? best_slack : 0.0)
                              << " best_cost="
                              << (has_best_cost ? best_cost : failed_sta_cost)
                              << " improvement="
                              << ((has_current_cost && has_best_cost)
                                      ? current_cost - best_cost
                                      : 0.0)
                              << '\n';
        }

        if (row_restore_failed) {
            if (!restore_full_baseline()) {
                logger->warn(
                    utl::TDP,
                    1084,
                    "Timing-driven HBT reallocation aborted: failed to restore full baseline after trial restore failed on net {}.",
                    candidate_connection.net_name);
                return false;
            }
            logger->warn(utl::TDP,
                         1076,
                         "Timing-driven HBT reallocation skipping net {}: failed to restore baseline after candidate trial.",
                         candidate_connection.net_name);
            continue;
        }

        if (!row_has_valid_sta_cost) {
            if (!restore_full_baseline()) {
                logger->warn(
                    utl::TDP,
                    1083,
                    "Timing-driven HBT reallocation aborted: failed to restore full baseline after skipping net {} with no valid STA candidate cost.",
                    candidate_connection.net_name);
                return false;
            }
            logger->warn(utl::TDP,
                         1075,
                         "Timing-driven HBT reallocation skipping net {}: no valid STA candidate cost.",
                         candidate_connection.net_name);
            continue;
        }

        active_row_indices.push_back(net_idx);
    }

    if (active_row_indices.empty()) {
        if (!restore_full_baseline()) {
            logger->warn(
                utl::TDP,
                1085,
                "Timing-driven HBT reallocation aborted: failed to restore full baseline after all candidate rows were skipped.");
            return false;
        }
        logger->info(utl::TDP,
                     1087,
                     "Timing-driven HBT reallocation produced no valid assignment.");
        return false;
    }

    std::unordered_set<std::string> active_net_names;
    active_net_names.reserve(active_row_indices.size());
    for (size_t net_idx : active_row_indices) {
        active_net_names.insert(problem.connections[net_idx].net_name);
    }

    std::vector<char> keep_candidate(problem.candidate_hbts.size(), false);
    for (size_t hbt_idx = 0; hbt_idx < problem.candidate_hbts.size(); ++hbt_idx) {
        const auto& candidate_hbt = problem.candidate_hbts[hbt_idx];
        const std::string owner_net = allocator->getNetForHBT(candidate_hbt.hbt_id);
        keep_candidate[hbt_idx] = owner_net.empty()
                                  || active_net_names.find(owner_net) != active_net_names.end();
    }

    std::vector<int> remapped_candidate_index(problem.candidate_hbts.size(), -1);
    HBTAllocator::TimingDrivenReallocationProblem active_problem;
    active_problem.worst_slack = problem.worst_slack;
    active_problem.candidate_hbts.reserve(problem.candidate_hbts.size());
    for (size_t hbt_idx = 0; hbt_idx < problem.candidate_hbts.size(); ++hbt_idx) {
        if (!keep_candidate[hbt_idx]) {
            continue;
        }
        remapped_candidate_index[hbt_idx] =
            static_cast<int>(active_problem.candidate_hbts.size());
        active_problem.candidate_hbts.push_back(problem.candidate_hbts[hbt_idx]);
    }

    active_problem.connections.reserve(active_row_indices.size());
    active_problem.candidate_indices_by_connection.reserve(active_row_indices.size());

    std::vector<std::vector<double>> active_sta_cost_matrix;
    active_sta_cost_matrix.reserve(active_row_indices.size());
    std::unordered_set<std::string> retained_active_net_names;
    retained_active_net_names.reserve(active_row_indices.size());
    for (size_t net_idx : active_row_indices) {
        std::vector<size_t> remapped_indices;
        remapped_indices.reserve(problem.candidate_indices_by_connection[net_idx].size());
        std::vector<double> remapped_costs(active_problem.candidate_hbts.size(), failed_sta_cost);
        bool row_has_valid_remapped_sta_cost = false;

        for (size_t original_hbt_idx : problem.candidate_indices_by_connection[net_idx]) {
            if (original_hbt_idx >= remapped_candidate_index.size()) {
                continue;
            }
            const int new_hbt_idx = remapped_candidate_index[original_hbt_idx];
            if (new_hbt_idx < 0) {
                continue;
            }
            remapped_indices.push_back(static_cast<size_t>(new_hbt_idx));
            const double remapped_cost = sta_cost_matrix[net_idx][original_hbt_idx];
            remapped_costs[static_cast<size_t>(new_hbt_idx)] = remapped_cost;
            if (remapped_cost != failed_sta_cost) {
                row_has_valid_remapped_sta_cost = true;
            }
        }

        if (remapped_indices.empty()) {
            logger->warn(utl::TDP,
                         1081,
                         "Timing-driven HBT reallocation skipping net {}: all legal HBT candidates are still owned by skipped nets.",
                         problem.connections[net_idx].net_name);
            continue;
        }
        if (!row_has_valid_remapped_sta_cost) {
            logger->warn(utl::TDP,
                         1082,
                         "Timing-driven HBT reallocation skipping net {}: all remaining HBT candidates have invalid STA cost after filtering skipped-owner HBTs.",
                         problem.connections[net_idx].net_name);
            continue;
        }

        active_problem.connections.push_back(problem.connections[net_idx]);
        active_problem.candidate_indices_by_connection.push_back(std::move(remapped_indices));
        active_sta_cost_matrix.push_back(std::move(remapped_costs));
        retained_active_net_names.insert(problem.connections[net_idx].net_name);
    }

    if (active_problem.connections.empty() || active_problem.candidate_hbts.empty()) {
        if (!restore_full_baseline()) {
            logger->warn(
                utl::TDP,
                1086,
                "Timing-driven HBT reallocation aborted: failed to restore full baseline after filtering active rows and candidate HBTs.");
            return false;
        }
        logger->info(utl::TDP,
                     1080,
                     "Timing-driven HBT reallocation produced no valid assignment.");
        return false;
    }

    active_problem.previous_assignments.reserve(problem.previous_assignments.size());
    for (const auto& [net_name, hbt_id] : problem.previous_assignments) {
        if (retained_active_net_names.find(net_name) != retained_active_net_names.end()) {
            active_problem.previous_assignments.emplace_back(net_name, hbt_id);
        }
    }

    allocator->restoreAssignmentSnapshot(baseline_assignment);
    syncHBTAssignmentsToNetwork3D(connections);

    auto result = allocator->commitTimingDrivenReallocationWithCostMatrix(
        active_problem, active_sta_cost_matrix);

    if (rematch_debug_log) {
        std::unordered_map<std::string, size_t> active_row_by_net;
        active_row_by_net.reserve(active_problem.connections.size());
        for (size_t row_idx = 0; row_idx < active_problem.connections.size(); ++row_idx) {
            active_row_by_net.emplace(active_problem.connections[row_idx].net_name, row_idx);
        }
        std::unordered_map<std::string, size_t> active_hbt_index_by_id;
        active_hbt_index_by_id.reserve(active_problem.candidate_hbts.size());
        for (size_t hbt_idx = 0; hbt_idx < active_problem.candidate_hbts.size(); ++hbt_idx) {
            active_hbt_index_by_id.emplace(active_problem.candidate_hbts[hbt_idx].hbt_id, hbt_idx);
        }
        for (const auto& allocation : result.results) {
            if (!allocation.success || allocation.allocated_hbt == nullptr) {
                continue;
            }
            const auto row_it = active_row_by_net.find(allocation.net_name);
            const auto hbt_it = active_hbt_index_by_id.find(allocation.allocated_hbt->id);
            if (row_it == active_row_by_net.end() || hbt_it == active_hbt_index_by_id.end()) {
                continue;
            }
            const double criticality = criticality_by_net.count(allocation.net_name)
                                           ? criticality_by_net[allocation.net_name]
                                           : 0.0;
            const double final_cost = active_sta_cost_matrix[row_it->second][hbt_it->second];
            const double final_slack = criticality > 0.0 ? -final_cost / criticality : 0.0;
            rematch_debug_log << "[hbta-rematch-eval] step=final_assign net="
                              << allocation.net_name
                              << " hbt=" << allocation.allocated_hbt->id
                              << " hbt_xy=(" << allocation.allocated_hbt->location.x
                              << "," << allocation.allocated_hbt->location.y << ")"
                              << " slack=" << final_slack
                              << " criticality=" << criticality
                              << " sta_cost=" << final_cost
                              << '\n';
        }
    }

    if (!allocator->hasPendingTimingDrivenReallocation()
        || result.total_nets == 0
        || !result.allSuccessful()) {
        syncHBTAssignmentsToNetwork3D(connections);
        tdp::SteinerTreeGenerator restore_steiner_generator(timingdp_);
        restore_steiner_generator.createSteninerRoute();
        timingdp_->executeTimingDriven(true);
        logger->info(utl::TDP,
                     1078,
                     "Timing-driven HBT reallocation produced no valid assignment.");
        return false;
    }

    // 对求解后的完整分配做一次全局 timing
    syncHBTAssignmentsToNetwork3D(connections);
    tdp::SteinerTreeGenerator steiner_generator(timingdp_);
    steiner_generator.createSteninerRoute();
    if (!timingdp_->executeTimingDriven(true)) {
        allocator->rollbackPendingTimingDrivenReallocation();
        syncHBTAssignmentsToNetwork3D(connections);
        steiner_generator.createSteninerRoute();
        timingdp_->executeTimingDriven(true);
        logger->warn(utl::TDP,
                     1072,
                     "Timing-driven HBT reallocation rolled back: tentative timing update failed.");
        return false;
    }

    const auto new_stats = timing_engine->getTimingStats();
    const bool accept_reallocation =
        HBTAllocator::shouldAcceptTimingDrivenReallocation(old_stats.late_wns,
                                                           old_stats.late_tns,
                                                           new_stats.late_wns,
                                                           new_stats.late_tns);

    // 确认这次重分配值得保留，还是应该回滚到 baseline 状态。
    if (!accept_reallocation) {
        allocator->rollbackPendingTimingDrivenReallocation();
        syncHBTAssignmentsToNetwork3D(connections);
        steiner_generator.createSteninerRoute();
        timingdp_->executeTimingDriven(true);
        logger->info(utl::TDP,
                     1073,
                     "Timing-driven HBT reallocation rejected: WNS {:.20f}->{:.20f}, TNS {:.20f}->{:.20f}.",
                     old_stats.late_wns,
                     new_stats.late_wns,
                     old_stats.late_tns,
                     new_stats.late_tns);
        return false;
    }

    allocator->acceptPendingTimingDrivenReallocation();

    logger->info(utl::TDP,
                 1074,
                 "Timing-driven HBT reallocation accepted for {} nets: WNS {:.20f}->{:.20f}, TNS {:.20f}->{:.20f}.",
                 result.successful_allocations,
                 old_stats.late_wns,
                 new_stats.late_wns,
                 old_stats.late_tns,
                 new_stats.late_tns);
    return true;
}

void HBTDatabaseInterface::reset() {
    timingdp_ = nullptr;
    allocator_.resetAllAllocations();
    extraction_stats_ = ExtractionStats();
    cell_layer_cache_.clear();
}

// ============================================================================
// ExtractionStats 实现
// ============================================================================

void HBTDatabaseInterface::ExtractionStats::print() const {
    std::cout << "=== HBT Network Extraction Statistics ===" << std::endl;
    std::cout << "  Total nets analyzed: " << total_nets << std::endl;
    std::cout << "  Inter-layer nets: " << inter_layer_nets << std::endl;
    std::cout << "  Bottom layer cells: " << bottom_layer_cells << std::endl;
    std::cout << "  Upper layer cells: " << upper_layer_cells << std::endl;
    std::cout << "  Failed conversions: " << failed_conversions << std::endl;
    
    if (total_nets > 0) {
        double inter_layer_ratio = static_cast<double>(inter_layer_nets) / total_nets * 100.0;
        std::cout << "  Inter-layer ratio: " << inter_layer_ratio << "%" << std::endl;
    }
}

} // namespace hbta
