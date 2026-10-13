#pragma once

#include "HBTTypes.h"

#include <string>
#include <vector>

namespace hbta {

// 规划结果里的单条分配记录
struct HBTQuadTreePlannedAssignment {
    const InterLayerConnection* connection = nullptr;
    const HBTInfo* hbt = nullptr;
    // 物理 routing 估计，单位 DBU
    double cost = 0.0;
    // 匹配优化目标值。
    double matching_cost = 0.0;
    int depth = 0;
};

// 四叉树求解阶段统计
struct HBTQuadTreeAllocationStats {
    // HBT 四叉树空间索引统计
    int quadtree_nodes = 0;
    int quadtree_leaf_nodes = 0;
    int quadtree_max_depth = 0;
    int quadtree_max_node_candidates = 0;
    long long quadtree_total_node_candidates = 0;

    // 求解阶段统计
    int visited_nodes = 0;
    int split_attempts = 0;
    int successful_split_nodes = 0;
    int partial_successful_split_nodes = 0;
    int local_solves = 0;
    int region_fallback_solves = 0;
    int coarse_sampled_solves = 0;
    int coarse_fallback_solves = 0;

    // 局部子问题统计
    int local_solve_samples = 0;
    long long total_coarse_input_candidates = 0;
    long long total_coarse_output_candidates = 0;
    long long total_local_candidate_pool_size = 0;
    long long total_local_connection_count = 0;
    long long total_local_candidate_entries = 0;
    int max_local_matrix_rows = 0;
    int max_local_matrix_cols = 0;
    long long max_local_matrix_entries = 0;
};

// 四叉树规划器返回值
struct HBTQuadTreeAllocationResult {
    bool success = false;
    std::string error_message;
    std::vector<HBTQuadTreePlannedAssignment> planned_assignments;
    HBTQuadTreeAllocationStats stats;
};

enum class HBTMatchingCostMode {
    // 基础匹配代价：routing-aware physical cost。
    Physical
};

// 四叉树分配器配置。
struct HBTQuadTreeAllocatorConfig {
    // 区域宽或高小于该阈值时，不再继续分裂。
    double quad_tree_min_wh_dbu = 50000.0;
    //采样因子
    size_t coarse_sample_target_factor = 6;

    HBTMatchingCostMode matching_cost_mode = HBTMatchingCostMode::Physical;
};

class HBTQuadTreeAllocator {
public:
    // 使用默认配置执行规划
    static HBTQuadTreeAllocationResult planAssignments(
        const std::vector<const InterLayerConnection*>& connections,
        const std::vector<const HBTInfo*>& hbts,
        double pitch_dbu);

    // 使用外部配置执行规划，便于实验批量调参
    static HBTQuadTreeAllocationResult planAssignments(
        const std::vector<const InterLayerConnection*>& connections,
        const std::vector<const HBTInfo*>& hbts,
        double pitch_dbu,
        const HBTQuadTreeAllocatorConfig& config);

    // 使用外部提供的真实 STA cost matrix 执行一对一匹配。
    static HBTQuadTreeAllocationResult planAssignmentsWithCostMatrix(
        const std::vector<const InterLayerConnection*>& connections,
        const std::vector<const HBTInfo*>& hbts,
        const std::vector<std::vector<double>>& cost_matrix);

    // 基于 net-HBT 候选二分图连通分量分解的局部 Hungarian。
    static HBTQuadTreeAllocationResult planAssignmentsWithCostMatrixByComponents(
        const std::vector<const InterLayerConnection*>& connections,
        const std::vector<const HBTInfo*>& hbts,
        const std::vector<std::vector<double>>& cost_matrix,
        const std::vector<std::vector<size_t>>& candidate_indices_by_connection);
};

} // namespace hbta
