#include "HBTQuadTreeAllocator.h"
#include "HBTACostModel.h"
#include "HBTQuadTree.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace hbta {

namespace {

using ReservedHBTSet = std::unordered_set<const HBTInfo*>;
using HBTQuadTreeNode = HBTQuadTree::Node;
constexpr double kInvalidAssignmentCost = 1e20;

double computePhysicalConnectionCost(const InterLayerConnection& connection,
                                     const Point& hbt_location) {
    // 物理线长统计采用 routing-aware MST 估计，用于结果记录与统计输出。
    return HBTACostModel::computeRoutingAwareMSTCost(connection, hbt_location);
}

double computeMatchingConnectionCost(const InterLayerConnection& connection,
                                     const Point& hbt_location) {
    // 初始匹配阶段的边权采用基础物理匹配代价。
    return HBTACostModel::computeMatchingCost(connection, hbt_location);
}

std::vector<std::vector<double>> buildPhysicalCostMatrix(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& candidate_hbts) {
    // 构建基础匹配代价矩阵：矩阵行为 net，矩阵列为候选 HBT。
    std::vector<std::vector<double>> cost_matrix(
        connections.size(), std::vector<double>(candidate_hbts.size(), 0.0));

    for (size_t i = 0; i < connections.size(); ++i) {
        for (size_t j = 0; j < candidate_hbts.size(); ++j) {
            cost_matrix[i][j] =
                computeMatchingConnectionCost(*connections[i], candidate_hbts[j]->location);
        }
    }

    return cost_matrix;
}

std::vector<std::vector<double>> buildCostMatrix(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& candidate_hbts,
    const HBTQuadTreeAllocatorConfig& config) {
    (void)config;
    return buildPhysicalCostMatrix(connections, candidate_hbts);
}

bool solveAssignment(const std::vector<std::vector<double>>& cost,
                     std::vector<int>* assignment,
                     double* min_cost) {
    if (!assignment || !min_cost) {
        return false;
    }

    // n 为 net 数量，m 为候选 HBT 数量；矩形 Hungarian 要求 n <= m。
    const int n = static_cast<int>(cost.size());
    if (n == 0) {
        assignment->clear();
        *min_cost = 0.0;
        return true;
    }

    const int m = static_cast<int>(cost[0].size());
    if (m == 0 || n > m) {
        return false;
    }

    for (int i = 1; i < n; ++i) {
        if (static_cast<int>(cost[i].size()) != m) {
            return false;
        }
    }

    // 在矩形代价矩阵上运行 Hungarian，求满足一对一约束的最小总代价匹配。
    const double inf = std::numeric_limits<double>::max() / 4.0;
    std::vector<double> u(n + 1, 0.0);  // 行势能。
    std::vector<double> v(m + 1, 0.0);  // 列势能。
    // p[j] = i 表示第 j 列当前匹配第 i 行；0 表示该列未匹配。
    std::vector<int> p(m + 1, 0);
    // way 用于记录增广路径，最终回溯并更新匹配关系。
    std::vector<int> way(m + 1, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, inf);
        std::vector<char> used(m + 1, false);

        do {
            used[j0] = true;
            const int i0 = p[j0];
            double delta = inf;
            int j1 = 0;

            for (int j = 1; j <= m; ++j) {
                if (used[j]) {
                    continue;
                }

                const double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }

            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        do {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    // 将列到行的内部匹配结果转换为 assignment[row] = column 的外部表示。
    assignment->assign(n, -1);
    for (int j = 1; j <= m; ++j) {
        if (p[j] != 0) {
            (*assignment)[p[j] - 1] = j - 1;
        }
    }

    *min_cost = -v[0];
    return true;
}

bool hasValidCostMatrixShape(const std::vector<const InterLayerConnection*>& connections,
                             const std::vector<const HBTInfo*>& hbts,
                             const std::vector<std::vector<double>>& cost_matrix,
                             std::string* error_message) {
    if (cost_matrix.size() != connections.size()) {
        if (error_message) {
            *error_message = "STA cost matrix row count does not match connection count";
        }
        return false;
    }

    for (const auto& row : cost_matrix) {
        if (row.size() != hbts.size()) {
            if (error_message) {
                *error_message =
                    "STA cost matrix column count does not match HBT candidate count";
            }
            return false;
        }
    }

    return true;
}

bool isInvalidAssignmentCost(double cost) {
    return !std::isfinite(cost) || cost >= (kInvalidAssignmentCost * 1.0);
}

void recordLocalMatrixStats(size_t rows, size_t cols, HBTQuadTreeAllocationStats* stats) {
    if (stats == nullptr) {
        return;
    }

    const long long row_count = static_cast<long long>(rows);
    const long long col_count = static_cast<long long>(cols);
    const long long entries = row_count * col_count;
    stats->local_solve_samples++;
    stats->total_local_connection_count += row_count;
    stats->total_local_candidate_pool_size += col_count;
    stats->total_local_candidate_entries += entries;
    if (row_count > stats->max_local_matrix_rows) {
        stats->max_local_matrix_rows = static_cast<int>(row_count);
    }
    if (col_count > stats->max_local_matrix_cols) {
        stats->max_local_matrix_cols = static_cast<int>(col_count);
    }
    if (entries > stats->max_local_matrix_entries) {
        stats->max_local_matrix_entries = entries;
    }
}

bool hasCandidateIndex(const std::vector<size_t>& candidate_indices, size_t candidate_index) {
    return std::binary_search(
        candidate_indices.begin(), candidate_indices.end(), candidate_index);
}

bool buildCandidateConnectedComponents(
    size_t net_count,
    size_t hbt_count,
    const std::vector<std::vector<size_t>>& candidate_indices_by_connection,
    std::vector<std::vector<size_t>>* net_components,
    std::vector<std::vector<size_t>>* hbt_components,
    std::string* error_message) {
    if (candidate_indices_by_connection.size() != net_count) {
        if (error_message) {
            *error_message = "Candidate index row count does not match connection count";
        }
        return false;
    }

    std::vector<std::vector<size_t>> net_indices_by_hbt(hbt_count);
    for (size_t net_idx = 0; net_idx < net_count; ++net_idx) {
        if (candidate_indices_by_connection[net_idx].empty()) {
            if (error_message) {
                *error_message = "A connection has no legal HBT candidates";
            }
            return false;
        }

        for (size_t hbt_idx : candidate_indices_by_connection[net_idx]) {
            if (hbt_idx >= hbt_count) {
                if (error_message) {
                    *error_message = "Candidate HBT index is out of range";
                }
                return false;
            }
            net_indices_by_hbt[hbt_idx].push_back(net_idx);
        }
    }

    std::vector<char> visited_nets(net_count, false);
    std::vector<char> visited_hbts(hbt_count, false);
    std::deque<std::pair<bool, size_t>> pending_nodes;

    net_components->clear();
    hbt_components->clear();
    net_components->reserve(net_count);
    hbt_components->reserve(net_count);

    for (size_t start_net_idx = 0; start_net_idx < net_count; ++start_net_idx) {
        if (visited_nets[start_net_idx]) {
            continue;
        }

        visited_nets[start_net_idx] = true;
        pending_nodes.emplace_back(true, start_net_idx);
        net_components->push_back({});
        hbt_components->push_back({});
        auto& component_nets = net_components->back();
        auto& component_hbts = hbt_components->back();

        while (!pending_nodes.empty()) {
            const auto [is_net_node, node_index] = pending_nodes.front();
            pending_nodes.pop_front();

            if (is_net_node) {
                component_nets.push_back(node_index);
                for (size_t hbt_idx : candidate_indices_by_connection[node_index]) {
                    if (visited_hbts[hbt_idx]) {
                        continue;
                    }
                    visited_hbts[hbt_idx] = true;
                    pending_nodes.emplace_back(false, hbt_idx);
                }
                continue;
            }

            component_hbts.push_back(node_index);
            for (size_t net_idx : net_indices_by_hbt[node_index]) {
                if (visited_nets[net_idx]) {
                    continue;
                }
                visited_nets[net_idx] = true;
                pending_nodes.emplace_back(true, net_idx);
            }
        }
    }

    return true;
}

HBTRegion computeConnectionBounds(const InterLayerConnection& connection) {
    // 使用真实 pin 位置构造 net 的空间包围盒
    Point seed = connection.suggested_hbt_location;
    size_t first_load_index = 0;
    if (!connection.driver_cell_name.empty()) {
        seed = connection.driver_location;
    } else if (!connection.load_locations.empty()) {
        seed = connection.load_locations.front();
        first_load_index = 1;
    }

    double min_x = seed.x;
    double max_x = seed.x;
    double min_y = seed.y;
    double max_y = seed.y;

    auto update_bounds = [&](const Point& point) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    };

    if (!connection.driver_cell_name.empty()) {
        update_bounds(connection.driver_location);
    }
    for (size_t i = first_load_index; i < connection.load_locations.size(); ++i) {
        update_bounds(connection.load_locations[i]);
    }

    return HBTRegion(min_x, min_y, max_x, max_y);
}

int findContainingChild(const HBTRegion& connection_bounds,
                        const std::array<HBTRegion, 4>& child_regions) {
    // 只有当 net 包围盒完全落入某个子区域时，才允许该 net 下沉到子节点处理。
    for (int i = 0; i < static_cast<int>(child_regions.size()); ++i) {
        if (child_regions[i].containsRegion(connection_bounds)) {
            return i;
        }
    }
    return -1;
}

std::vector<const HBTInfo*> collectAvailableCandidates(const HBTQuadTreeNode& node,
                                                       const ReservedHBTSet& reserved_hbts) {
    // 收集当前四叉树节点内尚未被全局或局部预留的 HBT。
    std::vector<const HBTInfo*> available_hbts;
    available_hbts.reserve(node.candidate_hbts.size());
    for (const HBTInfo* hbt : node.candidate_hbts) {
        if (reserved_hbts.find(hbt) == reserved_hbts.end()) {
            available_hbts.push_back(hbt);
        }
    }
    return available_hbts;
}

size_t countAvailableCandidates(const HBTQuadTreeNode& node,
                                const ReservedHBTSet& reserved_hbts) {
    // 仅统计当前节点中仍可用于本轮局部匹配的 HBT 数量。
    size_t available_count = 0;
    for (const HBTInfo* hbt : node.candidate_hbts) {
        if (reserved_hbts.find(hbt) == reserved_hbts.end()) {
            ++available_count;
        }
    }
    return available_count;
}

std::vector<const HBTInfo*> buildCoarseSampledCandidateSet(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& region_candidates,
    HBTQuadTreeAllocationStats* stats,
    const HBTQuadTreeAllocatorConfig& config) {
    // 对局部区域内候选 HBT 做规则采样，以控制 Hungarian 代价矩阵规模。
    if (connections.empty() || region_candidates.empty()) {
        return region_candidates;
    }

    // 采样预算约为 net 数量乘以 coarse_sample_target_factor。
    const size_t sample_budget = std::min(
        region_candidates.size(),
        connections.size() * config.coarse_sample_target_factor);
    if (sample_budget >= region_candidates.size()
        || sample_budget < connections.size()
        || sample_budget == 0) {
        return region_candidates;
    }

    // 根据 HBT 阵列行列范围估计采样步长。
    int min_row = std::numeric_limits<int>::max();
    int max_row = std::numeric_limits<int>::min();
    int min_col = std::numeric_limits<int>::max();
    int max_col = std::numeric_limits<int>::min();
    for (const HBTInfo* hbt : region_candidates) {
        if (hbt == nullptr) {
            continue;
        }
        min_row = std::min(min_row, hbt->row);
        max_row = std::max(max_row, hbt->row);
        min_col = std::min(min_col, hbt->col);
        max_col = std::max(max_col, hbt->col);
    }
    if (min_row > max_row || min_col > max_col) {
        return region_candidates;
    }

    const int row_span = std::max(1, max_row - min_row + 1);
    const int col_span = std::max(1, max_col - min_col + 1);
    const double total_grid_slots = static_cast<double>(row_span) * static_cast<double>(col_span);
    int skip_col = std::max(
        static_cast<int>(std::lround(
            total_grid_slots
            / std::max(1.0,
                       static_cast<double>(connections.size())
                           * static_cast<double>(config.coarse_sample_target_factor)))),
        1);
    if (skip_col <= 1) {
        return region_candidates;
    }

    // 按行列坐标进行斜向规则采样，尽量保持候选在空间上的覆盖性。
    auto sample_with_skip = [&](int stride) {
        std::vector<const HBTInfo*> sampled;
        sampled.reserve(std::min(sample_budget, region_candidates.size()));

        for (const HBTInfo* hbt : region_candidates) {
            if (hbt == nullptr) {
                continue;
            }

            const int local_row = hbt->row - min_row;
            const int local_col = hbt->col - min_col;
            const int row_offset = local_row % stride;
            if (((local_col - row_offset) % stride) == 0) {
                sampled.push_back(hbt);
            }
        }

        if (sampled.size() > sample_budget) {
            sampled.resize(sample_budget);
        }
        return sampled;
    };

    std::vector<const HBTInfo*> coarse_candidates = sample_with_skip(skip_col);
    // 若采样后候选数不足以形成一对一匹配，则放宽采样步长。
    if (coarse_candidates.size() < connections.size() && skip_col > 1) {
        skip_col = std::max(skip_col / 2, 1);
        if (skip_col > 1) {
            coarse_candidates = sample_with_skip(skip_col);
        }
    }

    if (coarse_candidates.size() < connections.size()
        || coarse_candidates.size() >= region_candidates.size()) {
        // 采样无效或无法保证可行性时，回退为完整区域候选集合。
        return region_candidates;
    }

    std::sort(coarse_candidates.begin(), coarse_candidates.end(), compareHBTCandidatePointers);
    if (stats) {
        stats->coarse_sampled_solves++;
        stats->total_coarse_input_candidates += static_cast<long long>(region_candidates.size());
        stats->total_coarse_output_candidates += static_cast<long long>(coarse_candidates.size());
    }
    return coarse_candidates;
}

bool solveSubproblemHungarian(const std::vector<const InterLayerConnection*>& connections,
                              const std::vector<const HBTInfo*>& candidate_hbts,
                              const HBTQuadTreeAllocatorConfig& config,
                              std::vector<int>* assignment,
                              double* min_cost,
                              std::vector<double>* assigned_matching_costs) {
    if (!assignment || !min_cost || !assigned_matching_costs) {
        return false;
    }
    if (connections.empty()) {
        assignment->clear();
        assigned_matching_costs->clear();
        *min_cost = 0.0;
        return true;
    }
    if (candidate_hbts.size() < connections.size()) {
        // 候选 HBT 少于 net 数量时，不存在合法的一对一匹配。
        return false;
    }

    // 构建稠密代价矩阵：行是 net，列是候选 HBT。
    const auto cost_matrix = buildCostMatrix(connections, candidate_hbts, config);
    if (!solveAssignment(cost_matrix, assignment, min_cost)) {
        return false;
    }

    // 保存每个 net 实际选中边的匹配代价，供结果记录使用。
    assigned_matching_costs->assign(connections.size(), 0.0);
    for (size_t i = 0; i < assignment->size(); ++i) {
        const int candidate_index = (*assignment)[i];
        if (candidate_index < 0
            || candidate_index >= static_cast<int>(cost_matrix[i].size())) {
            return false;
        }
        (*assigned_matching_costs)[i] = cost_matrix[i][candidate_index];
    }

    return true;
}

bool solveLocalSubproblem(const HBTQuadTreeNode& node,
                          const std::vector<const InterLayerConnection*>& connections,
                          int depth,
                          const HBTQuadTreeAllocatorConfig& config,
                          ReservedHBTSet* reserved_hbts,
                          std::vector<HBTQuadTreePlannedAssignment>* planned_assignments,
                          HBTQuadTreeAllocationStats* stats) {
    if (!reserved_hbts || !planned_assignments) {
        return false;
    }
    if (connections.empty()) {
        return true;
    }

    if (stats) {
        stats->local_solves++;
    }

    // 局部求解仅使用当前四叉树节点内可用候选。
    const auto region_candidates = collectAvailableCandidates(node, *reserved_hbts);
    if (region_candidates.size() < connections.size()) {
        return false;
    }

    auto candidate_pool =
        buildCoarseSampledCandidateSet(connections, region_candidates, stats, config);

    // 优先使用采样候选求解；若采样候选不可行，再回退到完整区域候选。
    std::vector<int> assignment;
    std::vector<double> assigned_matching_costs;
    double min_cost = 0.0;
    bool solved = solveSubproblemHungarian(
        connections, candidate_pool, config, &assignment, &min_cost, &assigned_matching_costs);
    if (!solved) {
        if (candidate_pool.size() != region_candidates.size()) {
            if (stats) {
                stats->coarse_fallback_solves++;
            }
            candidate_pool = region_candidates;
            solved = solveSubproblemHungarian(
                connections, candidate_pool, config, &assignment, &min_cost, &assigned_matching_costs);
        }
    }

    if (!solved) {
        return false;
    }

    // 两阶段提交：先在局部集合中校验所有选择，确认无冲突后再合并到共享结果。
    ReservedHBTSet local_reservations;
    std::vector<HBTQuadTreePlannedAssignment> local_assignments;
    local_assignments.reserve(connections.size());

    for (size_t i = 0; i < connections.size(); ++i) {
        const int candidate_index = assignment[i];
        if (candidate_index < 0 || candidate_index >= static_cast<int>(candidate_pool.size())) {
            return false;
        }

        const HBTInfo* selected_hbt = candidate_pool[candidate_index];
        // 防止已占用 HBT、已预留 HBT 或局部重复选择进入最终规划结果。
        if (selected_hbt == nullptr
            || selected_hbt->used
            || reserved_hbts->find(selected_hbt) != reserved_hbts->end()
            || !local_reservations.insert(selected_hbt).second) {
            return false;
        }

        HBTQuadTreePlannedAssignment local_assignment;
        local_assignment.connection = connections[i];
        local_assignment.hbt = selected_hbt;
        // cost 记录物理线长估计，matching_cost 记录本轮 Hungarian 使用的优化目标值。
        local_assignment.cost =
            computePhysicalConnectionCost(*connections[i], selected_hbt->location);
        local_assignment.matching_cost =
            i < assigned_matching_costs.size()
                ? assigned_matching_costs[i]
                : computeMatchingConnectionCost(*connections[i], selected_hbt->location);
        local_assignment.depth = depth;
        local_assignments.push_back(local_assignment);
    }

    if (stats) {
        // 记录局部矩阵规模，便于评估四叉树分解与候选采样的效果。
        const long long rows = static_cast<long long>(connections.size());
        const long long cols = static_cast<long long>(candidate_pool.size());
        const long long entries = rows * cols;
        stats->local_solve_samples++;
        stats->total_local_candidate_pool_size += cols;
        stats->total_local_connection_count += rows;
        stats->total_local_candidate_entries += entries;
        if (rows > stats->max_local_matrix_rows) {
            stats->max_local_matrix_rows = static_cast<int>(rows);
        }
        if (cols > stats->max_local_matrix_cols) {
            stats->max_local_matrix_cols = static_cast<int>(cols);
        }
        if (entries > stats->max_local_matrix_entries) {
            stats->max_local_matrix_entries = entries;
        }
    }

    // 局部校验通过后，统一提交预留 HBT 和计划分配结果。
    reserved_hbts->insert(local_reservations.begin(), local_reservations.end());
    planned_assignments->insert(
        planned_assignments->end(), local_assignments.begin(), local_assignments.end());
    return true;
}

bool solveNodeRecursively(const HBTQuadTreeNode& node,
                          const std::vector<const InterLayerConnection*>& connections,
                          int depth,
                          const HBTQuadTreeAllocatorConfig& config,
                          ReservedHBTSet* reserved_hbts,
                          std::vector<HBTQuadTreePlannedAssignment>* planned_assignments,
                          HBTQuadTreeAllocationStats* stats) {
    if (!reserved_hbts || !planned_assignments || !stats) {
        return false;
    }

    stats->visited_nodes++;

    if (connections.empty()) {
        return true;
    }

    // 可行性检查：当前节点可用 HBT 数不足时，该区域无法完成一对一匹配。
    if (countAvailableCandidates(node, *reserved_hbts) < connections.size()) {
        return false;
    }

    if (!node.has_children) {
        // 叶子节点不再继续分裂，直接在当前区域内求解局部匹配。
        return solveLocalSubproblem(
            node, connections, depth, config, reserved_hbts, planned_assignments, stats);
    }

    // 非叶子节点先尝试将 net 下沉到子区域；不能安全下沉的 net 保留在父节点处理。
    std::array<std::vector<const InterLayerConnection*>, 4> child_connections;
    std::vector<const InterLayerConnection*> parent_connections;
    std::array<size_t, 4> child_available_counts = {0, 0, 0, 0};

    std::array<HBTRegion, 4> child_regions = {
        node.children[0]->region,
        node.children[1]->region,
        node.children[2]->region,
        node.children[3]->region,
    };

    // 预先统计各子区域可用 HBT 数量，用于判断子问题是否具备基本可行性。
    for (size_t i = 0; i < child_regions.size(); ++i) {
        child_available_counts[i] =
            countAvailableCandidates(*node.children[i], *reserved_hbts);
    }

    for (const InterLayerConnection* connection : connections) {
        // 仅当 net 包围盒被子区域完全包含时才下推，否则留在父节点处理跨区关系。
        const HBTRegion connection_bounds = computeConnectionBounds(*connection);
        const int child_index = findContainingChild(connection_bounds, child_regions);

        if (child_index < 0 || child_available_counts[child_index] == 0) {
            parent_connections.push_back(connection);
            continue;
        }

        child_connections[child_index].push_back(connection);
    }

    for (size_t i = 0; i < child_connections.size(); ++i) {
        // 子区候选不足时，将该子区 net 整体上提到父节点统一求解。
        if (child_connections[i].size() > child_available_counts[i]) {
            parent_connections.insert(
                parent_connections.end(),
                child_connections[i].begin(),
                child_connections[i].end());
            child_connections[i].clear();
        }
    }

    bool has_child_subproblem = false;
    for (const auto& child_group : child_connections) {
        if (!child_group.empty()) {
            has_child_subproblem = true;
            break;
        }
    }

    if (has_child_subproblem) {
        stats->split_attempts++;

        // 子树求解失败时只回滚该子树结果，并将相关 net 上提到父节点。
        ReservedHBTSet working_reserved = *reserved_hbts;
        auto working_assignments = *planned_assignments;
        std::vector<const InterLayerConnection*> promoted_parent_connections = parent_connections;
        int successful_child_solves = 0;
        int failed_child_solves = 0;

        for (size_t i = 0; i < child_connections.size(); ++i) {
            if (child_connections[i].empty()) {
                continue;
            }

            ReservedHBTSet child_reserved_snapshot = working_reserved;
            auto child_assignments_snapshot = working_assignments;

            // 每个子树独立递归求解；失败不会污染其他已成功子树的正式状态。
            if (!solveNodeRecursively(*node.children[i],
                                      child_connections[i],
                                      depth + 1,
                                      config,
                                      &working_reserved,
                                      &working_assignments,
                                      stats)) {
                // 子树失败时恢复快照，并把该子树 net 交给父节点后续统一分配。
                working_reserved = std::move(child_reserved_snapshot);
                working_assignments = std::move(child_assignments_snapshot);
                promoted_parent_connections.insert(promoted_parent_connections.end(),
                                                   child_connections[i].begin(),
                                                   child_connections[i].end());
                failed_child_solves++;
            } else {
                successful_child_solves++;
            }
        }

        // 子区处理后，父节点继续求解跨区 net 与失败上提的 net。
        if (solveLocalSubproblem(node,
                                 promoted_parent_connections,
                                 depth,
                                 config,
                                 &working_reserved,
                                 &working_assignments,
                                 stats)) {
            *reserved_hbts = std::move(working_reserved);
            *planned_assignments = std::move(working_assignments);
            if (successful_child_solves > 0) {
                stats->successful_split_nodes++;
                if (failed_child_solves > 0) {
                    stats->partial_successful_split_nodes++;
                }
            }
            return true;
        }

        // 若“子区结果 + 父区剩余 net”组合不可行，则回退为当前节点整体局部求解。
        stats->region_fallback_solves++;
    }

    // 无有效子问题或区域级 fallback 时，在当前节点完整候选集合上求解。
    return solveLocalSubproblem(
        node, connections, depth, config, reserved_hbts, planned_assignments, stats);
}

void appendUniqueCandidates(const std::vector<const HBTInfo*>& candidates,
                            std::unordered_set<const HBTInfo*>* seen_hbts,
                            std::vector<const HBTInfo*>* reduced_hbts) {
    if (!seen_hbts || !reduced_hbts) {
        return;
    }

    // 对多个路径收集到的 HBT 去重。
    for (const HBTInfo* hbt : candidates) {
        if (hbt != nullptr && seen_hbts->insert(hbt).second) {
            reduced_hbts->push_back(hbt);
        }
    }
}

void appendCandidatesToConnections(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& candidates,
    const std::unordered_map<const InterLayerConnection*, size_t>& connection_to_index,
    std::vector<std::vector<const HBTInfo*>>* candidate_hbts_by_connection) {
    if (candidates.empty() || candidate_hbts_by_connection == nullptr) {
        return;
    }

    for (const InterLayerConnection* connection : connections) {
        if (connection == nullptr) {
            continue;
        }
        const auto index_it = connection_to_index.find(connection);
        if (index_it == connection_to_index.end()) {
            continue;
        }
        auto& per_connection_candidates =
            (*candidate_hbts_by_connection)[index_it->second];
        per_connection_candidates.insert(per_connection_candidates.end(),
                                         candidates.begin(),
                                         candidates.end());
    }
}

bool collectReducedCandidatesRecursively(
    const HBTQuadTreeNode& node,
    const std::vector<const InterLayerConnection*>& connections,
    const HBTQuadTreeAllocatorConfig& config,
    const std::unordered_map<const InterLayerConnection*, size_t>& connection_to_index,
    std::unordered_set<const HBTInfo*>* seen_hbts,
    std::vector<const HBTInfo*>* reduced_hbts,
    std::vector<std::vector<const HBTInfo*>>* candidate_hbts_by_connection,
    HBTQuadTreeAllocationStats* stats) {
    if (!seen_hbts || !reduced_hbts || !stats) {
        return false;
    }

    stats->visited_nodes++;

    if (connections.empty()) {
        return true;
    }
    if (node.candidate_hbts.empty()) {
        return false;
    }

    if (!node.has_children) {
        // 叶子节点只做候选压缩，不做 assignment。
        appendUniqueCandidates(node.candidate_hbts, seen_hbts, reduced_hbts);
        appendCandidatesToConnections(connections,
                                      node.candidate_hbts,
                                      connection_to_index,
                                      candidate_hbts_by_connection);
        return !node.candidate_hbts.empty();
    }

    std::array<std::vector<const InterLayerConnection*>, 4> child_connections;
    std::vector<const InterLayerConnection*> parent_connections;
    std::array<size_t, 4> child_candidate_counts = {0, 0, 0, 0};
    const std::array<HBTRegion, 4> child_regions = {
        node.children[0]->region,
        node.children[1]->region,
        node.children[2]->region,
        node.children[3]->region,
    };

    for (size_t i = 0; i < child_regions.size(); ++i) {
        // reduction 路径按子树原始候选规模判断是否下沉。
        child_candidate_counts[i] = node.children[i]->candidate_hbts.size();
    }

    for (const InterLayerConnection* connection : connections) {
        // 只有完整落入子区的 net 才继续下钻。
        const HBTRegion connection_bounds = computeConnectionBounds(*connection);
        const int child_index = findContainingChild(connection_bounds, child_regions);
        if (child_index < 0 || child_candidate_counts[child_index] == 0) {
            parent_connections.push_back(connection);
            continue;
        }

        child_connections[child_index].push_back(connection);
    }

    for (size_t i = 0; i < child_connections.size(); ++i) {
        // 子区候选不足时，上提回父区处理。
        if (child_connections[i].size() > child_candidate_counts[i]) {
            parent_connections.insert(parent_connections.end(),
                                      child_connections[i].begin(),
                                      child_connections[i].end());
            child_connections[i].clear();
        }
    }

    bool has_child_subproblem = false;
    int successful_child_reductions = 0;
    int failed_child_reductions = 0;
    for (size_t i = 0; i < child_connections.size(); ++i) {
        if (child_connections[i].empty()) {
            continue;
        }

        has_child_subproblem = true;
        // 子树失败时，把这部分 net 上提给父区。
        if (collectReducedCandidatesRecursively(*node.children[i],
                                                child_connections[i],
                                                config,
                                                connection_to_index,
                                                seen_hbts,
                                                reduced_hbts,
                                                candidate_hbts_by_connection,
                                                stats)) {
            successful_child_reductions++;
        } else {
            parent_connections.insert(parent_connections.end(),
                                      child_connections[i].begin(),
                                      child_connections[i].end());
            failed_child_reductions++;
        }
    }

    if (has_child_subproblem) {
        stats->split_attempts++;
        if (successful_child_reductions > 0) {
            stats->successful_split_nodes++;
            if (failed_child_reductions > 0) {
                stats->partial_successful_split_nodes++;
            }
        }
    }

    if (!parent_connections.empty()) {
        // 父区 residual net 直接继承当前区域完整候选，不再做额外采样。
        appendUniqueCandidates(node.candidate_hbts, seen_hbts, reduced_hbts);
        appendCandidatesToConnections(parent_connections,
                                      node.candidate_hbts,
                                      connection_to_index,
                                      candidate_hbts_by_connection);
        return !node.candidate_hbts.empty();
    }

    return true;
}

} // namespace

HBTQuadTreeAllocationResult HBTQuadTreeAllocator::planAssignments(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& hbts,
    double pitch_dbu) {
    // 默认配置入口，供基础批量分配流程调用。
    return planAssignments(connections, hbts, pitch_dbu, HBTQuadTreeAllocatorConfig{});
}

HBTQuadTreeAllocationResult HBTQuadTreeAllocator::planAssignments(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& hbts,
    double pitch_dbu,
    const HBTQuadTreeAllocatorConfig& config) {
    HBTQuadTreeAllocationResult result;

    if (connections.empty()) {
        result.success = true;
        return result;
    }

    if (hbts.empty()) {
        result.error_message = "HBT candidate set is empty";
        return result;
    }

    // 基于候选 HBT 建立四叉树空间索引，后续递归过程只在局部区域内求解匹配。
    HBTQuadTreeConfig tree_config;
    tree_config.min_region_wh_dbu = config.quad_tree_min_wh_dbu;
    HBTQuadTree hbt_tree(hbts, tree_config);
    if (hbt_tree.empty() || hbt_tree.root() == nullptr) {
        result.error_message = "Failed to build HBT quadtree";
        return result;
    }

    // 保存四叉树结构统计信息，便于分析区域划分规模。
    const auto& tree_stats = hbt_tree.stats();
    result.stats.quadtree_nodes = tree_stats.nodes;
    result.stats.quadtree_leaf_nodes = tree_stats.leaf_nodes;
    result.stats.quadtree_max_depth = tree_stats.max_depth;
    result.stats.quadtree_max_node_candidates = tree_stats.max_node_candidates;
    result.stats.quadtree_total_node_candidates = tree_stats.total_node_candidates;

    ReservedHBTSet reserved_hbts;
    // 从根节点开始递归求解；reserved_hbts 用于保证不同局部子问题不会选择同一个 HBT。
    if (!solveNodeRecursively(*hbt_tree.root(),
                              connections,
                              0,
                              config,
                              &reserved_hbts,
                              &result.planned_assignments,
                              &result.stats)) {
        result.error_message = "Quadtree allocation failed";
        return result;
    }

    result.success = true;
    return result;
}

HBTQuadTreeAllocationResult HBTQuadTreeAllocator::planAssignmentsWithCostMatrix(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& hbts,
    const std::vector<std::vector<double>>& cost_matrix) {
    HBTQuadTreeAllocationResult result;

    if (connections.empty()) {
        result.success = true;
        return result;
    }

    if (hbts.empty()) {
        result.error_message = "HBT candidate set is empty";
        return result;
    }

    if (!hasValidCostMatrixShape(connections, hbts, cost_matrix, &result.error_message)) {
        return result;
    }

    std::vector<int> assignment;
    double min_cost = 0.0;
    if (!solveAssignment(cost_matrix, &assignment, &min_cost)) {
        result.error_message = "STA cost matrix assignment failed";
        return result;
    }

    result.stats.local_solves = 1;
    recordLocalMatrixStats(connections.size(), hbts.size(), &result.stats);

    result.planned_assignments.reserve(connections.size());
    for (size_t i = 0; i < connections.size(); ++i) {
        const int candidate_index = assignment[i];
        if (candidate_index < 0 || static_cast<size_t>(candidate_index) >= hbts.size()) {
            result.error_message = "STA cost matrix produced an invalid assignment";
            result.planned_assignments.clear();
            return result;
        }

        HBTQuadTreePlannedAssignment planned_assignment;
        planned_assignment.connection = connections[i];
        planned_assignment.hbt = hbts[candidate_index];
        planned_assignment.cost =
            computePhysicalConnectionCost(*connections[i], hbts[candidate_index]->location);
        planned_assignment.matching_cost = cost_matrix[i][candidate_index];
        planned_assignment.depth = 0;
        result.planned_assignments.push_back(planned_assignment);
    }

    result.success = true;
    return result;
}

HBTQuadTreeAllocationResult HBTQuadTreeAllocator::planAssignmentsWithCostMatrixByComponents(
    const std::vector<const InterLayerConnection*>& connections,
    const std::vector<const HBTInfo*>& hbts,
    const std::vector<std::vector<double>>& cost_matrix,
    const std::vector<std::vector<size_t>>& candidate_indices_by_connection) {
    HBTQuadTreeAllocationResult result;

    if (connections.empty()) {
        result.success = true;
        return result;
    }

    if (hbts.empty()) {
        result.error_message = "HBT candidate set is empty";
        return result;
    }

    if (!hasValidCostMatrixShape(connections, hbts, cost_matrix, &result.error_message)) {
        return result;
    }

    std::vector<std::vector<size_t>> net_components;
    std::vector<std::vector<size_t>> hbt_components;
    if (!buildCandidateConnectedComponents(connections.size(),
                                           hbts.size(),
                                           candidate_indices_by_connection,
                                           &net_components,
                                           &hbt_components,
                                           &result.error_message)) {
        return result;
    }

    result.planned_assignments.resize(connections.size());
    std::vector<char> assigned_nets(connections.size(), false);
    std::vector<char> assigned_hbts(hbts.size(), false);

    for (size_t component_idx = 0; component_idx < net_components.size(); ++component_idx) {
        const auto& component_nets = net_components[component_idx];
        const auto& component_hbts = hbt_components[component_idx];
        if (component_nets.empty()) {
            continue;
        }
        if (component_nets.size() > component_hbts.size()) {
            result.error_message =
                "A candidate component has fewer HBTs than connections";
            result.planned_assignments.clear();
            return result;
        }

        result.stats.local_solves++;
        std::unordered_map<size_t, size_t> local_hbt_index_by_global;
        local_hbt_index_by_global.reserve(component_hbts.size());
        for (size_t local_col = 0; local_col < component_hbts.size(); ++local_col) {
            local_hbt_index_by_global.emplace(component_hbts[local_col], local_col);
        }

        std::vector<std::vector<double>> local_cost_matrix(
            component_nets.size(),
            std::vector<double>(component_hbts.size(), kInvalidAssignmentCost));

        for (size_t local_row = 0; local_row < component_nets.size(); ++local_row) {
            const size_t global_net_idx = component_nets[local_row];
            bool row_has_finite_legal_cost = false;
            for (size_t global_hbt_idx : candidate_indices_by_connection[global_net_idx]) {
                const auto local_col_it =
                    local_hbt_index_by_global.find(global_hbt_idx);
                if (local_col_it == local_hbt_index_by_global.end()) {
                    continue;
                }

                const double candidate_cost = cost_matrix[global_net_idx][global_hbt_idx];
                local_cost_matrix[local_row][local_col_it->second] = candidate_cost;
                if (!isInvalidAssignmentCost(candidate_cost)) {
                    row_has_finite_legal_cost = true;
                }
            }

            if (!row_has_finite_legal_cost) {
                result.error_message =
                    "A connection has no finite STA cost among its legal HBT candidates";
                result.planned_assignments.clear();
                return result;
            }
        }

        std::vector<int> local_assignment;
        double component_min_cost = 0.0;
        if (!solveAssignment(local_cost_matrix, &local_assignment, &component_min_cost)) {
            result.error_message = "STA cost matrix component assignment failed";
            result.planned_assignments.clear();
            return result;
        }

        recordLocalMatrixStats(
            component_nets.size(), component_hbts.size(), &result.stats);

        for (size_t local_row = 0; local_row < component_nets.size(); ++local_row) {
            const int local_col = local_assignment[local_row];
            if (local_col < 0
                || static_cast<size_t>(local_col) >= component_hbts.size()) {
                result.error_message =
                    "STA cost matrix component assignment produced an invalid column";
                result.planned_assignments.clear();
                return result;
            }

            const double selected_cost = local_cost_matrix[local_row][local_col];
            if (isInvalidAssignmentCost(selected_cost)) {
                result.error_message =
                    "STA cost matrix component assignment selected an illegal edge";
                result.planned_assignments.clear();
                return result;
            }

            const size_t global_net_idx = component_nets[local_row];
            const size_t global_hbt_idx = component_hbts[local_col];
            if (assigned_nets[global_net_idx] || assigned_hbts[global_hbt_idx]) {
                result.error_message =
                    "STA cost matrix component assignment is not one-to-one";
                result.planned_assignments.clear();
                return result;
            }
            if (!hasCandidateIndex(candidate_indices_by_connection[global_net_idx],
                                   global_hbt_idx)) {
                result.error_message =
                    "STA cost matrix component assignment violated candidate constraints";
                result.planned_assignments.clear();
                return result;
            }

            assigned_nets[global_net_idx] = true;
            assigned_hbts[global_hbt_idx] = true;

            HBTQuadTreePlannedAssignment planned_assignment;
            planned_assignment.connection = connections[global_net_idx];
            planned_assignment.hbt = hbts[global_hbt_idx];
            planned_assignment.cost = computePhysicalConnectionCost(
                *connections[global_net_idx], hbts[global_hbt_idx]->location);
            planned_assignment.matching_cost = selected_cost;
            planned_assignment.depth = 0;
            result.planned_assignments[global_net_idx] = planned_assignment;
        }
    }

    for (size_t net_idx = 0; net_idx < connections.size(); ++net_idx) {
        if (!assigned_nets[net_idx]
            || result.planned_assignments[net_idx].connection == nullptr
            || result.planned_assignments[net_idx].hbt == nullptr) {
            result.error_message =
                "STA cost matrix component assignment missed a connection";
            result.planned_assignments.clear();
            return result;
        }
    }

    result.success = true;
    return result;
}

} // namespace hbta
