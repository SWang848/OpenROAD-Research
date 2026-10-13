// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "sta/Sta.hh"
#include "../eva/TimingEngine3D.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sta { class Parasitics; }

namespace tdp {

class SteinerTreeGenerator;

struct LRParams {
    double lambda_wirelength = 1.0;
    double lambda_timing = 1.0;
    double lambda_alignment = 1.0;
    double beta = 0.01;
    int max_iterations = 20;
    double convergence_threshold = 1e-3;
    int    window_y = 4;
    int    window_x = 4;
};
/**
 * @class Placement3D
 * @brief Base class for 3D placement algorithms
 * 
 * Following DPL design pattern, algorithms should be implemented
 * in separate classes that use Timingdp for data access.
 */
class Placement3D
{
public:
    Placement3D(Timingdp* timingdp);
    virtual ~Placement3D();
    
    static Architecture3D::Row3D* findCellRow(Node3D* cell, 
                                            const std::vector<Architecture3D::Row3D*>& rows);
    // 倒序遍历
    void buildTraverseOrder(Timingdp* timingdp);
    void rundetailedplacement();
    void verifyTraverseOrder() const;

    // 刷新单个单元的时序缓存（从 STA 拉取最新值）
    void refreshPinTiming(Node3D* cell);
    void refreshPinTiming2(Node3D* cell);
    void compareLocalVsFullSTA(Node3D* cell);
    bool isSequentialCell(const Node3D* cell) const;
    bool isPrimaryInput(const Node3D* cell) const;
    bool isPrimaryOutput(const Node3D* cell) const;

protected:
    Timingdp* timingdp_ = nullptr;
    sta::dbSta* sta_ = nullptr;
    utl::Logger* logger_ = nullptr;
    Architecture3D* arch_ = nullptr;
    Network3D* network_ = nullptr;
    TimingEngine3D* timing_engine_ = nullptr;
    std::list<Node3D*> traverse_list;
    double clock_period_T_ = 0.0;

    bool isInputPin(const Pin3D* pin);
    bool isOutputPin(const Pin3D* pin);
    bool isPrimaryInput(const Pin3D* pin);
    bool isPrimaryOutput(const Pin3D* pin);
};

/**
 * @class LSP3D
 * @brief Timing-driven 3D placement algorithm using LR method
 * 
 * This algorithm:
 * 1. Considers both upper and lower layers simultaneously
 * 2. Moves cells based on cost sensitivity minimization
 * 3. Uses LR (Lagrangian Relaxation) method
 * 4. Handles inter-layer timing optimization
 */
class LSP3D : public Placement3D
{
public:

    
    LSP3D(Timingdp* timingdp, const LRParams& params = LRParams());
    ~LSP3D() override;
    
    // 主要布局方法
    void place();

    // 窗口寻找相关方法
    std::vector<Pointdp3D> timingDrivenWindowSearch(
        Node3D* cell,
        int num_candidates = 33,
        int window_width = 0,
        int window_height = 0);

    // 行占用空间索引 (row occupancy index): [layer][geometric row idx] -> cells.
    // 用于把 findAvailableSlotsInRow 从 O(全设计cell) 全扫降到 O(行内cell)。
    using RowIndex = std::vector<std::vector<std::vector<Node3D*>>>;
    void buildRowIndexInto(RowIndex& dst) const;
    int  rowIdxForY(int layer, DBU y_bottom) const;   // 几何行号；越界返回 -1
    void updateRowIndexMove(Node3D* cell,
                            const Pointdp3D& old_loc,
                            const Pointdp3D& new_loc);
    RowIndex row_index_;

    // 时序权重初始化
    void initWeight(Timingdp* timingdp, TimingEngine3D* timing_engine);
    void collectW2Stats();
    void collectW3Stats(
        const std::unordered_set<Node3D*>& crit_cells,
        const std::unordered_map<Pin3D*, std::pair<double,double>>& snap);
    void computeWl();
    void propagateWp();
    // void combineWeights();  // commented: now computed in calculateCellTimingCost

    // 成本计算方法
    double calculateCellCost(Node3D* cell, double alpha, double omega,
                             TimingEngine3D* timing_engine = nullptr) const;
    
    double calculateCellTimingCost(Node3D* cell, double alpha, double omega,
                                   TimingEngine3D* timing_engine) const;
    
    // 局部搜索方法
    bool discreteLocalSearch(Node3D* cell, int search_windows, double alpha, double omega,
                             TimingEngine3D* timing_engine = nullptr);

    // No-commit variant of discreteLocalSearch: evaluate every candidate with
    // a compute-only local STA (Timingdp::analyzeOnly) read from a PtGraph,
    // WITHOUT committing to the global graph; commit only the chosen best.
    bool discreteLocalSearchNoCommit(Node3D* cell, int search_windows,
                                     double alpha, double omega,
                                     TimingEngine3D* timing_engine = nullptr);

    // Phase-2 parallel variant: evaluate all candidates concurrently using
    // per-thread {engine, Steiner, parasitic container}; odb is never moved
    // during evaluation (moved-pin override). Commits only the best candidate.
    bool discreteLocalSearchParallel(Node3D* cell, int search_windows,
                                     double alpha, double omega,
                                     TimingEngine3D* timing_engine = nullptr);

    void discreteLocalSearchAll(TimingEngine3D* timing_engine = nullptr);

    // Spike 1: verify concurrent read-only safety of the local-graph STA.
    void runSearchConcurrencySpike();

    // Spike 2: verify per-thread (non-scene) trial parasitics feed the local
    // STA identically to the global scene path.
    void runTrialParasiticSpike();

    // Spike 3: verify the combined inner loop (per-thread Steiner + trial
    // parasitics + local STA) runs concurrently, bit-identically to serial.
    void runConcurrentTrialSpike();

    // Spike 2B (R2): verify the moved-pin override trial geometry matches
    // physically moving the cell in odb.
    void runTrialMovedSpike();

    // Extractor compare: at current geometry, scene extractor
    // (updateSteinerTrees/estimateTdpWireParasitic) vs local extractor
    // (buildLocalNetParasiticInto), both through analyzeOnly.
    void runExtractorCompareSpike();

    // Cost attribution: on the SAME trial PtGraph (moved candidate), decompose
    // the par-vs-no-commit candidate-cost gap into its sub-steps:
    //   driver-frame (③), arc frozen-vs-live (②), baseline source (⑥).
    void runCostAttributionSpike();

    // 构建关键路径优先遍历顺序（top-k 违反约束路径上的 cell 优先）
    std::vector<Node3D*> buildCriticalPathPriorityOrder(
        TimingEngine3D* timing_engine,
        double top_k);

    // 辅助函数
    Pin3D* findDriverPin(Net3D* net, Pin3D* exclude_pin) const;
    std::vector<Pin3D*> getCellInputPins(Node3D* cell) const;
    void flipCell(Node3D* cell) const;
    std::string orientToString(const odb::dbOrientType& orient) const;
    // 权重更新
    void updatePinWeights();
    void balanceCellWeights(Node3D* cell);
    void distributeCellWeights(Node3D* cell);
    // cell Arc delay 计算
    void cellArcDelays(Node3D* cell) const;
    // 辅助函数
    std::vector<Pin3D*> getCellOutputPins(Node3D* cell) const;

    // 刷新受影响的 pin 时序缓存（cell + 其 net 上所有 pin）
    void refreshAffectedPinTiming(Node3D* cell);
    void refreshAffectedPinTiming2(Node3D* cell);

    // 同 refreshAffectedPinTiming，但到达时间从局部 PtGraph 读取（不读全局图）；
    // PtGraph 中缺失的 pin 回退到全局基线值。
    void refreshAffectedPinTimingFromPtGraph(Node3D* cell, lgista::PtGraph* pt_graph);

private:
    LRParams params_;
    int current_iteration_ = 0;
    double current_cost_ = 0.0;
    double previous_cost_ = 0.0;

    // Phase-2 parallel pool + thread-local cost helpers
    void ensureParallelPool(int num_threads);
    void buildLocalArrivalMap(
        Node3D* cell, lgista::PtGraph* pt_graph,
        std::unordered_map<Pin3D*, std::pair<double, double>>& out) const;
    // Per-input-pin internal gate arc delays read from the (moved-geometry)
    // local PtGraph, mirroring cellArcDelays but at the trial position. Keyed by
    // the input Pin3D: {early(min from 0), late(max from 0)}.
    void buildLocalArcMap(
        Node3D* cell, lgista::PtGraph* pt_graph,
        std::unordered_map<Pin3D*, std::pair<double, double>>& out) const;
    double calculateCellTimingCostLocal(
        Node3D* cell, double alpha, double omega,
        const std::unordered_map<Pin3D*, std::pair<double, double>>& arr) const;
    double calculateCellTimingCostLocal(
        Node3D* cell, double alpha, double omega,
        const std::unordered_map<Pin3D*, std::pair<double, double>>& arr,
        const std::unordered_map<Pin3D*, std::pair<double, double>>& arc_map) const;

    std::vector<std::unique_ptr<lgista::LocalGraphIsta>> par_engines_;
    std::vector<std::unique_ptr<SteinerTreeGenerator>>   par_steiners_;
    std::vector<sta::Parasitics*>                        par_local_par_;
    int par_nthreads_ = 0;

    // Timingdp* timingdp_ = nullptr;
    // sta::dbSta* sta_ = nullptr;
    // utl::Logger* logger_;
    
    // 敏感度数据
    struct CellSensitivity {
        double timing_sensitivity = 0.0;
        double wirelength_sensitivity = 0.0;
        double alignment_sensitivity = 0.0;
        double total_sensitivity = 0.0;
    };
    
    std::unordered_map<Node3D*, CellSensitivity> sensitivities_;
    // 辅助数据结构
    struct TimingVector {
        double x;
        double y;
        double magnitude;
        TimingVector() : x(0.0), y(0.0), magnitude(0.0) {}
        TimingVector(double x_val, double y_val, double mag) 
            : x(x_val), y(y_val), magnitude(mag) {}
    };
    
    struct DirectionComponents {
        double left;
        double right;
        double up;
        double down;
        DirectionComponents() : left(0.0), right(0.0), up(0.0), down(0.0) {}
    };
    
    struct SearchWindow {
        DBU x;
        DBU y;
        DBU width;
        DBU height;
        SearchWindow() : x(0), y(0), width(0), height(0) {}
        SearchWindow(DBU x_val, DBU y_val, DBU w, DBU h)
            : x(x_val), y(y_val), width(w), height(h) {}
    };
    
    // 私有辅助方法
    void computeTimingVectors(
        Node3D* cell_j,
        Node3D* connected_cell,
        bool is_fan_in,
        double late_lm,
        double early_lm,
        TimingVector& late_vector,
        TimingVector& early_vector);
    
    Pointdp3D getCellCenter(Node3D* cell) const;
    
    bool computeDirectionComponents(
        Node3D* cell,
        DirectionComponents& late_components,
        DirectionComponents& early_components);
    
    void decomposeToComponents(
        const TimingVector& vector,
        DirectionComponents& components);
    
    SearchWindow determineSearchWindow(
        Node3D* cell,
        DBU window_width,
        DBU window_height);
    
    std::vector<Pointdp3D> generateCandidatePositions(
        const SearchWindow& window,
        int num_candidates,
        DBU step_x = 0,
        DBU step_y = 0);
    
    SearchWindow clampSearchWindow(
        const SearchWindow& window,
        Node3D* cell);
    
    // 候选位置生成相关方法
    std::vector<Pointdp3D> generateCandidatePositions(
        Node3D* cell,
        const SearchWindow& window,
        int num_candidates,
        DBU step_x = 0);

    std::vector<Pointdp3D> generateCandidatePositions2(
        Node3D* cell,
        const SearchWindow& window,
        int num_candidates);
    
    std::vector<std::pair<DBU, DBU>> findAvailableSlotsInRow(
        Architecture3D::Row3D* row,
        Node3D* exclude_cell,
        DBU range_start,
        DBU range_end);
    
    std::vector<std::pair<DBU, DBU>> findAvailableSlotsInRowAligned(
        Architecture3D::Row3D* row,
        Node3D* exclude_cell,
        DBU range_start,
        DBU range_end,
        DBU site_width = 0,
        DBU site_spacing = 0);
    
    DBU alignToSite(DBU coord, DBU site_spacing, bool round_up) const;
};

} // namespace tdp