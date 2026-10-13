// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
// #include "../../src/infrastructure/SegmentManager3D.h" 
// #include "../../src/infrastructure/Architecture3D.h"
// #include "../../src/infrastructure/Network3D.h"
#include "../../src/infrastructure/Common3D.h"
// HBT分配器头文件
#include "../../src/hbta/HBTDatabaseInterface.h"

#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "db_sta/dbSta.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"
#include "rsz/Resizer.hh"
#include "grt/GRoute.h"
#include "grt/GlobalRouter.h"
#include "est/EstimateParasitics.h"
#include "dpl/Opendp.h"
#include "ifp/InitFloorplan.hh"
// Tcl解释器前向声明
// struct Tcl_Interp;

// Forward declarations
// namespace sta {
//     class dbSta;
//     class LibertyLibrary;
// }

// namespace odb {
//     class dbDatabase;
//     class dbChip;
//     class dbBlock;
//     class dbInst;
//     class dbNet;
//     class dbSite;
//     class dbRow;
// }

// namespace utl {
//     class Logger;
// }

namespace lgista {
class LocalGraphIsta;
class PtGraph;
}

namespace sta {
class ArcDelayCalc;
}

namespace tdp {

    class Architecture3D;      
    class Network3D;           
    class SegmentManager3D;    
    class Node3D;              
    class Segment3D;
    class TimingEngine3D;

    // using odb::Rect;
/**
 * @class Timingdp
 * @brief 3D placement database following DPL design pattern
 * 
 * Similar to DPL's Opendp class, this class manages:
 * 1. 3D architecture information (Architecture3D)
 * 2. 3D network/cell information (Network3D)
 * 3. Database loading and conversion
 * 4. Basic operations for placement algorithms
 * 
 * Algorithm logic should be implemented in separate classes.
 */
class Timingdp
{
public:
    Timingdp(odb::dbDatabase* db, 
             utl::Logger* logger, 
             sta::dbSta* sta, 
             rsz::Resizer* resizer,
             grt::GlobalRouter* global_router,
             dpl::Opendp* opendp);
    ~Timingdp();
    sta::dbSta* getSta() const { return sta_; }
    grt::GlobalRouter* getGlobalRouter() { return global_router_; }
    rsz::Resizer* getResizer() { return rs_; }
    est::EstimateParasitics* getEstimateParasitics() { return estimate_parasitics_; }
    TimingEngine3D* gettimingengine() { return timing_engine_.get(); }
    void setSta(sta::dbSta* sta);
    bool isIdealClockNet(odb::dbNet* db_net) const;
    bool isValidSigType(odb::dbNet* db_net) const;
    // 禁止拷贝（仿照DPL的Opendp）
    Timingdp(const Timingdp&) = delete;
    Timingdp& operator=(const Timingdp&) = delete;

    // ========================================================================
    // 工艺组合配置（由 -tech 决定上/下层各自工艺）
    // ========================================================================
    enum class Tech { ASAP7, NANGATE45 };
    struct TechProfile { Tech upper; Tech bottom; };

    // 解析 -tech 字符串（asap7_3D / nangate45_3D / asap7_nangate45_3D /
    // nangate45_asap7_3D）为 (upper, bottom) 工艺组合。
    static TechProfile parseProfile(const std::string& tech);
    // 单工艺 -> core site 名。
    static const char* techSiteName(Tech t);
    // 供 Tcl/SWIG 调用：给定层（"upper"/"bottom"）与 -tech，返回该层 site 名。
    static const char* tierSite(const std::string& tier, const std::string& tech);

    
    // ========================================================================
    // 数据库管理（类似DPL的Opendp）
    // ========================================================================
    
    bool destroyDatabase();
    
    // 3D PDK参数提取（仿照DPL从数据库获取）。layer 非空（"upper"/"bottom"）时按该层
    // 标定并设置 site 宽度；layer 为空时仅刷新层间/HBT 间距。
    bool extractPDKParametersFromDatabase(const std::string& layer);
    
    // ========================================================================
    // 转换控制（类似DPL的importDb/importClear）
    // ========================================================================
    
    bool importDb();      // 从OpenDB导入并转换到3D内部表示
    void importClear();   // 清除转换数据
    bool designset();
    void createnewrow(const std::string& layer);
    void detailedPlacement(int max_displacement_x,
                        int max_displacement_y,
                        const std::string& report_file_name = std::string(""),
                        double margin_um = 0,
                        const std::string& tech = std::string("asap7_nangate45_3D"),
                        const std::string& hbt = std::string(""));
    void reportTimingStats(int max_displacement_x,
                        int max_displacement_y,
                        double margin_um = 0,
                        const std::string& tech = std::string("asap7_nangate45_3D"),
                        const std::string& hbt = std::string(""));
    void legalizeTier(const std::string& tier,
                      double margin_um,
                      int max_displacement_x,
                      int max_displacement_y);
    void reportLegalizationStats() const;
    
    // ========================================================================
    // 数据访问（类似DPL的getter方法）
    // ========================================================================
    
    // 数据访问
    odb::dbDatabase* getDatabase() const { return db_; }
    odb::dbChip* getChip() const { return chip_; }
    odb::dbBlock* getBlock() const { return block_; }
    
    Architecture3D* getArchitecture3D() { return arch_3d_.get(); }
    const Architecture3D* getArchitecture3D() const { return arch_3d_.get(); }
    
    Network3D* getNetwork3D() { return network_3d_.get(); }
    const Network3D* getNetwork3D() const { return network_3d_.get(); }
    
    // 设计信息
    int64_t getDieWidth() const;
    int64_t getDieHeight() const;

    // 日志管理
    utl::Logger* getLogger() const;
    
    // 实例管理
    std::vector<odb::dbInst*>& getMovableInstances();
    std::vector<odb::dbInst*>& getFixedInstances();
    std::vector<odb::dbInst*>& getAllInstances();
    
    // 网络管理
    std::vector<odb::dbNet*>& getNets();
    
    // 3D特定操作
    bool setInstanceLocation3D(odb::dbInst* inst, int64_t x, int64_t y, int64_t z);
    
    // 层间连接
    // bool addHBTConnection(const std::string& bottom_inst, const std::string& upper_inst);
    // std::vector<std::pair<std::string, std::string>> getHBTConnections() const;
    
    // 验证和调试
    bool check3DPlacementLegality() const;
    void adjustNodesOrient();
    
    // ========================================================================
    // HBT分配接口（新增）
    // ========================================================================
    
    // 初始化HBT分配器
    bool initializeHBTAllocator(const std::string& hbt_array_json, int dbu_per_micron = 1000);
    
    // 提取跨层连接
    std::vector<hbta::InterLayerConnection> extractInterLayerConnections();
    
    // 获取HBT分配器（需要先初始化）
    hbta::HBTAllocator* getHBTAllocator();
    const hbta::HBTAllocator* getHBTAllocator() const;
    
    // 执行HBT分配
    hbta::HBTAllocator::BatchNetAllocationResult allocateHBTsForDesign();
    // 基于 STA late slack 对已分配 HBT 做局部重匹配。
    bool runTimingDrivenHBTReallocation(
        const std::vector<hbta::InterLayerConnection>& connections);
    
    // 添加段管理接口
    void initializeSegmentManager();
    SegmentManager3D* getSegmentManager() { return segment_manager_; }
    const SegmentManager3D* getSegmentManager() const { return segment_manager_; }
    
    // 段管理包装方法
    std::vector<Segment3D*> getCellSegments(const Node3D* cell) const;
    bool addCellToSegment(Node3D* cell, int seg_id);
    bool removeCellFromSegment(const Node3D* cell, int seg_id);
    void initializeSegments();

    odb::dbDatabase* db_ = nullptr;
    odb::dbChip* chip_ = nullptr;
    odb::dbBlock* block_ = nullptr;
    odb::Rect core_;

    int max_displacement_x_ = 0;  // sites
    int max_displacement_y_ = 0;  // sites

    bool executeTimingDriven(bool run_journal_restore);
    void executeTimingDrivenI(Node3D* cell, int depth = 1);
    // Compute-only local STA: build a local timing graph around the cell and
    // propagate, WITHOUT committing back to the global sta::Graph. Returns a
    // newly-allocated PtGraph and TRANSFERS OWNERSHIP to the caller (wrap it in
    // a unique_ptr / delete it). Returns nullptr if the cell has no usable
    // driver. Raw pointer (not unique_ptr) so this header needs only a forward
    // declaration of lgista::PtGraph (keeps the SWIG wrappers happy).
    lgista::PtGraph* analyzeOnly(Node3D* cell, int depth = 2);
private:
    
    // 3D架构和网络（类似DPL的Architecture和Network）
    std::unique_ptr<Architecture3D> arch_3d_;
    std::unique_ptr<Network3D> network_3d_;
    std::unique_ptr<TimingEngine3D> timing_engine_;
    sta::dbSta* sta_ = nullptr; 
    rsz::Resizer* rs_ = nullptr;
    utl::Logger* logger_;
    grt::GlobalRouter* global_router_ = nullptr;
    est::EstimateParasitics* estimate_parasitics_ = nullptr;
    dpl::Opendp* opendp_ = nullptr;
    
    // 缓存数据
    std::vector<odb::dbInst*> movable_insts_;
    std::vector<odb::dbInst*> fixed_insts_;
    std::vector<odb::dbInst*> all_insts_;
    std::vector<odb::dbNet*> nets_;
    
    // 设计参数
    DBU die_width_ = 0;
    DBU die_height_ = 0;
    DBU die_depth_ = 0;
    DBU inter_layer_spacing_ = 0;
    DBU hbt_pitch_ = 0;
    
    // 工艺库信息
    std::string bottom_tech_lib_ = "NangateOpenCellLibrary.macro.mod.bottom.processed";
    std::string upper_tech_lib_ = "asap7sc7p5t_28_R_1x_220121a.upper.processed";
    std::string hbt_tech_lib_;

    // 本次运行的工艺组合（由 detailedPlacement 的 -tech 解析得到）。
    TechProfile profile_ = {Tech::ASAP7, Tech::NANGATE45};
    
    // PDK配置（仿照DPL从数据库动态获取）
    struct LayerInfo {
        DBU height = 0;
        DBU site_width = 0;
        std::string site_name;
        std::string tech_lib;
        std::string layer_type;  // "bottom", "upper", "unknown"
        int row_count = 0;
    };
    
    struct PDKLayerConfig {
        int layer_id = 0;
        std::string layer_name;
        DBU height = 0;
        DBU site_width = 0;
        std::string site_name;
        std::string tech_lib;
        int row_count = 0;
    };
    
    struct PDKConfig {
        std::vector<PDKLayerConfig> layers;
        DBU hbt_pitch = 0;
        DBU inter_layer_spacing = 0;
        
        void clear() {
            layers.clear();
            hbt_pitch = 0;
            inter_layer_spacing = 0;
        }
    };
    
    PDKConfig pdk_config_;
    
    // ========================================================================
    // 私有方法（类似DPL的私有方法）
    // ========================================================================
    
    // 数据库转换方法（类似DPL的dbToOpendp功能）
    bool convertRows(const std::string& layer); // 转换行（按传入层标定）
    bool convertInstances();                  // 转换实例
    bool convertNets();                       // 转换网络
    bool setUpPlacementGroups();              // 设置布局组
    
    // 辅助方法
    bool setupLayerArchitecture(const std::string& layer);
    
    // 数据缓存（类似DPL的缓存）
    void cacheDesignData();
    void cacheInstances();
    void cacheNets();
    
    // 内存管理（类似DPL的importClear）
    void cleanup();
    
    // 状态检查
    bool isDatabaseReady() const;
    bool isConverted() const;
    
    // 转换状态
    bool is_converted_ = false;
    bool is_loaded_ = false;
    
    // ========================================================================
    // 新增的私有辅助方法
    // ========================================================================
    
    // 从LEF文件名提取工艺库名称
    std::string extractLibraryNameFromLEF(const std::string& lef_file);
    
    // 从site信息推断层类型
    std::string determineLayerTypeFromSiteInfo(odb::dbSite* site);
    
    // 提取3D间距参数
    void extract3DSpacingParametersFromDatabase();

    
    // sta::dbSta* sta_ = nullptr;
    // 确定row属于哪个层
    // Architecture3D::LayerType determineLayerForRow(odb::dbRow* row);
    
    // 确定instance属于哪个层
    // Architecture3D::LayerType determineLayerForInstance(odb::dbInst* inst);

    
    // 添加段管理器
    SegmentManager3D* segment_manager_;

    // HBT 分配状态。保持为 Timingdp 成员，避免 initialize/extract/allocate
    // 之间落到不同的静态 interface 实例。
    std::unique_ptr<hbta::HBTDatabaseInterface> hbt_interface_;
    std::vector<hbta::InterLayerConnection> last_hbt_connections_;

    std::unordered_map<std::string, double> collectHBTNetSlacks(
        const std::vector<hbta::InterLayerConnection>& connections);

    // lgista local timing engine
    std::unique_ptr<lgista::LocalGraphIsta> ista_;

    //////
    ///test


};

} // namespace tdp