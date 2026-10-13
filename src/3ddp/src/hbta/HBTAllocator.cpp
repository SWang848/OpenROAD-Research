#include "HBTAllocator.h"
#include "HBTACostModel.h"
#include "HBTQuadTree.h"
#include "HBTQuadTreeAllocator.h"
#include "HBTRegion.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace hbta {

namespace {

double computeConnectionCost(const InterLayerConnection& connection, const Point& hbt_location) {
    // 统一采用 routing-aware MST 估计单个 net 经过候选 HBT 后的物理连线代价。
    return HBTACostModel::computeRoutingAwareMSTCost(connection, hbt_location);
}

HBTRegion computePinBoundingBox(const InterLayerConnection& connection) {
    double min_x = connection.suggested_hbt_location.x;
    double max_x = connection.suggested_hbt_location.x;
    double min_y = connection.suggested_hbt_location.y;
    double max_y = connection.suggested_hbt_location.y;

    auto include_point = [&](const Point& point) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    };

    if (!connection.driver_cell_name.empty()) {
        include_point(connection.driver_location);
    }
    for (const auto& load_location : connection.load_locations) {
        include_point(load_location);
    }

    return HBTRegion(min_x, min_y, max_x, max_y);
}

HBTRegion makeLocalBondingBox(const InterLayerConnection& connection,
                              const HBTInfo* current_hbt,
                              double half_size_dbu) {
    const HBTRegion pin_box = computePinBoundingBox(connection);
    const Point center = current_hbt != nullptr
                             ? current_hbt->location
                             : Point(pin_box.centerX(), pin_box.centerY());
    return HBTRegion(center.x - half_size_dbu,
                     center.y - half_size_dbu,
                     center.x + half_size_dbu,
                     center.y + half_size_dbu);
}

} // 匿名命名空间

bool HBTAllocator::initialize(const std::string& hbt_array_json, int dbu_per_micron) {
    dbu_per_micron_ = dbu_per_micron;

    // 加载 HBT 阵列，并建立后续分配所需的物理坐标与可用性状态。
    if (!hbt_array_.loadFromJSON(hbt_array_json, dbu_per_micron)) {
        std::cerr << "Failed to initialize HBTAllocator: cannot load HBT array" << std::endl;
        return false;
    }

    // 初始化后清空历史映射，保证新的 HBT 阵列状态与 net/HBT 映射一致。
    net_to_hbt_.clear();
    hbt_to_net_.clear();
    allocation_history_.clear();
    
    std::cout << "HBTAllocator initialized with " << hbt_array_.getTotalCount() 
              << " HBTs, " << hbt_array_.getAvailableCount() << " available" << std::endl;
    
    return true;
}

HBTAllocator::NetAllocationResult HBTAllocator::allocateHBTForNet(
    const InterLayerConnection& net_connection) {
    // 单网分配转为单元素批处理，保证单网与批量分配使用相同的四叉树 + Hungarian 主流程。
    const std::vector<InterLayerConnection> single_connection_batch = {net_connection};
    auto batch_result = allocateHBTsForNets(single_connection_batch);

    if (!batch_result.results.empty()) {
        return batch_result.results.front();
    }

    NetAllocationResult fallback_result;
    fallback_result.net_name = net_connection.net_name;
    fallback_result.success = false;
    fallback_result.error_message = "Single-net allocation returned no result";
    return fallback_result;
}

HBTAllocator::BatchNetAllocationResult HBTAllocator::allocateHBTsForNets(
    const std::vector<InterLayerConnection>& net_connections) {
    return allocateHBTsForNetsWithConfig(net_connections, HBTQuadTreeAllocatorConfig{});
}

HBTAllocator::BatchNetAllocationResult HBTAllocator::allocateHBTsForNetsWithConfig(
    const std::vector<InterLayerConnection>& net_connections,
    const HBTQuadTreeAllocatorConfig& quad_config) {
    return allocateHBTsForNetsWithConfigAndCandidates(
        net_connections, quad_config, nullptr);
}

HBTAllocator::BatchNetAllocationResult HBTAllocator::allocateHBTsForNetsWithConfigAndCandidates(
    const std::vector<InterLayerConnection>& net_connections,
    const HBTQuadTreeAllocatorConfig& quad_config,
    const std::vector<const HBTInfo*>* candidate_hbt_override) {
    BatchNetAllocationResult batch_result;
    batch_result.total_nets = net_connections.size();

    // 统一写入批处理结果；历史记录只保存真实发生的新分配，避免重复统计复用结果。
    auto append_result = [&](const NetAllocationResult& result, bool record_history = true) {
        batch_result.results.push_back(result);
        if (record_history) {
            allocation_history_.push_back(result);
        }
    };

    // 预检查：连接合法性与批内 net_name 唯一性；失败时整批返回错误，避免部分状态被写入。
    std::unordered_set<std::string> seen_net_names;
    seen_net_names.reserve(net_connections.size());
    for (const auto& conn : net_connections) {
        if (!validateNetConnection(conn)) {
            batch_result.failed_allocations = net_connections.size();
            for (const auto& c : net_connections) {
                NetAllocationResult result;
                result.success = false;
                result.net_name = c.net_name;
                result.error_message = "Invalid network connection";
                append_result(result);
            }
            std::cerr << "ERROR: Invalid network connection in batch" << std::endl;
            return batch_result;
        }

        if (!seen_net_names.insert(conn.net_name).second) {
            batch_result.failed_allocations = net_connections.size();
            const std::string error_message =
                "Duplicate net_name in batch: " + conn.net_name;
            for (const auto& c : net_connections) {
                NetAllocationResult result;
                result.success = false;
                result.net_name = c.net_name;
                result.error_message = error_message;
                append_result(result);
            }
            std::cerr << "ERROR: " << error_message << std::endl;
            return batch_result;
        }
    }

    // 已分配网络直接复用现有 HBT；未分配网络进入四叉树批量匹配。
    std::vector<const InterLayerConnection*> to_assign;
    to_assign.reserve(net_connections.size());
    for (const auto& conn : net_connections) {
        if (isNetAllocated(conn.net_name)) {
            const HBTInfo* existing_hbt = getHBTForNet(conn.net_name);
            NetAllocationResult result;
            result.net_name = conn.net_name;

            if (existing_hbt) {
                result.success = true;
                result.allocated_hbt = existing_hbt;
                result.estimated_wirelength = computeConnectionCost(conn, existing_hbt->location);
                result.driver_cell = conn.driver_cell_name;
                result.load_cells = conn.load_cell_names;
                batch_result.successful_allocations++;
                batch_result.total_wirelength += result.estimated_wirelength;
            } else {
                result.success = false;
                result.error_message = "Net already allocated but HBT not found";
                batch_result.failed_allocations++;
            }

            // 复用分配结果不写入历史记录，避免重复统计。
            append_result(result, false);
        } else {
            to_assign.push_back(&conn);
        }
    }

    // 默认只把当前未被占用的 HBT 作为候选；二次局部重匹配可传入候选集合
    std::vector<const HBTInfo*> available_hbts;
    std::unordered_set<std::string> seen_candidate_hbts;
    if (candidate_hbt_override != nullptr) {
        available_hbts.reserve(candidate_hbt_override->size());
        for (const HBTInfo* hbt : *candidate_hbt_override) {
            if (hbt != nullptr && !hbt->used && seen_candidate_hbts.insert(hbt->id).second) {
                available_hbts.push_back(hbt);
            }
        }
    } else {
        available_hbts.reserve(hbt_array_.getAvailableCount());
        for (const auto& hbt : hbt_array_.getAllHBTs()) {
            if (!hbt.used && seen_candidate_hbts.insert(hbt.id).second) {
                available_hbts.push_back(&hbt);
            }
        }
    }

    std::cout << "Quadtree matching input sizes: nets=" << to_assign.size()
              << ", available_hbts=" << available_hbts.size();
    if (candidate_hbt_override != nullptr) {
        std::cout << " (local override)";
    }
    std::cout << std::endl;

    // Hungarian 要求候选 HBT 数量不少于待分配 net 数量
    if (to_assign.size() > available_hbts.size()) {
        batch_result.failed_allocations += static_cast<int>(to_assign.size());
        for (const auto* connection : to_assign) {
            NetAllocationResult result;
            result.success = false;
            result.net_name = connection->net_name;
            result.error_message = "Insufficient HBTs. Required: "
                                 + std::to_string(to_assign.size())
                                 + ", Available: " + std::to_string(available_hbts.size());
            append_result(result);
        }
        std::cerr << "ERROR: Connection count (" << to_assign.size()
                  << ") exceeds available HBT count (" << available_hbts.size() << ")" << std::endl;
        return batch_result;
    }

    if (to_assign.empty()) {
        std::cout << "Batch allocation completed: " << batch_result.successful_allocations
                  << " successful, " << batch_result.failed_allocations << " failed" << std::endl;
        return batch_result;
    }

    // 四叉树负责缩小局部候选集合，局部子问题再通过 Hungarian 求最小总代价匹配。
    const auto quad_plan = HBTQuadTreeAllocator::planAssignments(
        to_assign, available_hbts, hbt_array_.getPitchDBU(), quad_config);
    if (!quad_plan.success) {
        batch_result.failed_allocations += static_cast<int>(to_assign.size());
        for (const auto* connection : to_assign) {
            NetAllocationResult result;
            result.success = false;
            result.net_name = connection->net_name;
            result.error_message = quad_plan.error_message.empty()
                                 ? "Quadtree allocation failed"
                                 : quad_plan.error_message;
            append_result(result);
        }
        std::cerr << "ERROR: "
                  << (quad_plan.error_message.empty()
                          ? "Quadtree allocation failed"
                          : quad_plan.error_message)
                  << std::endl;
        return batch_result;
    }

    // 将规划结果按 net_name 建索引，便于后续逐个提交到 HBTArray 与双向映射表。
    std::unordered_map<std::string, HBTQuadTreePlannedAssignment> plan_by_net;
    for (const auto& plan : quad_plan.planned_assignments) {
        if (plan.connection != nullptr) {
            plan_by_net[plan.connection->net_name] = plan;
        }
    }

    for (const auto* connection : to_assign) {
        NetAllocationResult result;
        result.net_name = connection->net_name;

        auto plan_it = plan_by_net.find(connection->net_name);
        if (plan_it == plan_by_net.end() || plan_it->second.hbt == nullptr) {
            result.success = false;
            result.error_message = "No planned HBT assignment found";
            batch_result.failed_allocations++;
            append_result(result);
            continue;
        }

        const HBTInfo* selected_hbt = plan_it->second.hbt;
        // 提交阶段再次标记 HBT 占用，防止规划结果与当前全局状态不一致。
        if (!hbt_array_.markHBTUsed(selected_hbt->id, connection->net_name)) {
            result.success = false;
            result.error_message = "Failed to mark HBT " + selected_hbt->id + " as used";
            batch_result.failed_allocations++;
            append_result(result);
            continue;
        }

        // 维护 net->HBT 与 HBT->net 双向映射，供查询、释放和二次重匹配回滚使用。
        net_to_hbt_[connection->net_name] = selected_hbt->id;
        hbt_to_net_[selected_hbt->id] = connection->net_name;

        result.success = true;
        result.allocated_hbt = selected_hbt;
        result.estimated_wirelength = plan_it->second.cost;
        result.driver_cell = connection->driver_cell_name;
        result.load_cells = connection->load_cell_names;

        batch_result.successful_allocations++;
        batch_result.total_wirelength += result.estimated_wirelength;
        append_result(result);
    }

    const double coarse_keep_ratio = quad_plan.stats.total_coarse_input_candidates > 0
                                         ? static_cast<double>(quad_plan.stats.total_coarse_output_candidates)
                                               / quad_plan.stats.total_coarse_input_candidates
                                         : 1.0;
    const double avg_local_candidates = quad_plan.stats.local_solve_samples > 0
                                            ? static_cast<double>(quad_plan.stats.total_local_candidate_pool_size)
                                                  / quad_plan.stats.local_solve_samples
                                            : 0.0;
    const double avg_candidates_per_connection = quad_plan.stats.total_local_connection_count > 0
                                                     ? static_cast<double>(quad_plan.stats.total_local_candidate_entries)
                                                           / quad_plan.stats.total_local_connection_count
                                                     : 0.0;
    std::cout << "Quadtree summary: nodes=" << quad_plan.stats.quadtree_nodes
              << ", local_solves=" << quad_plan.stats.local_solves
              << ", max_matrix=" << quad_plan.stats.max_local_matrix_rows
              << "x" << quad_plan.stats.max_local_matrix_cols
              << ", avg_local_candidates=" << avg_local_candidates
              << ", avg_candidates_per_connection=" << avg_candidates_per_connection
              << ", coarse_keep_ratio=" << coarse_keep_ratio << std::endl;

    const double total_wirelength_um = batch_result.total_wirelength
                                     / static_cast<double>(dbu_per_micron_);
    std::cout << "Batch allocation completed: " << batch_result.successful_allocations
              << " successful, " << batch_result.failed_allocations << " failed"
              << ", total estimated wirelength: " << total_wirelength_um << " um" << std::endl;

    return batch_result;
}

HBTAllocator::TimingDrivenReallocationProblem
HBTAllocator::prepareTimingDrivenReallocationProblem(
    const std::vector<InterLayerConnection>& net_connections,
    const std::unordered_map<std::string, double>& slack_by_net) const {
    return prepareTimingDrivenReallocationProblem(
        net_connections, slack_by_net, TimingDrivenReallocationConfig{});
}

HBTAllocator::TimingDrivenReallocationProblem
HBTAllocator::prepareTimingDrivenReallocationProblem(
    const std::vector<InterLayerConnection>& net_connections,
    const std::unordered_map<std::string, double>& slack_by_net,
    const TimingDrivenReallocationConfig& config) const {
    TimingDrivenReallocationProblem problem;

    // 根据 slack 找 critical nets，收集局部候选 HBT
    struct CriticalConnection {
        const InterLayerConnection* connection = nullptr;
        double slack = 0.0;
        double criticality = 0.0;
    };

    std::vector<CriticalConnection> critical_connections;
    critical_connections.reserve(net_connections.size());

    // 从 STA slack 表中筛选严格小于阈值且连接合法的网络；未出现在 slack 表中的 net 不参与本轮重匹配。
    for (const auto& connection : net_connections) {
        auto slack_it = slack_by_net.find(connection.net_name);
        if (slack_it == slack_by_net.end()) {
            continue;
        }
        if (slack_it->second >= config.critical_slack_threshold) {
            continue;
        }
        if (!validateNetConnection(connection)) {
            continue;
        }
        critical_connections.push_back({&connection, slack_it->second});
    }

    if (critical_connections.empty()) {
        std::cout << "Timing-driven HBT reallocation: no critical nets selected" << std::endl;
        return problem;
    }

    // 使用当前 critical 集合中的 worst slack 计算 criticality
    const auto worst_it = std::min_element(
        critical_connections.begin(),
        critical_connections.end(),
        [](const CriticalConnection& lhs, const CriticalConnection& rhs) {
            return lhs.slack < rhs.slack;
        });
    const double all_worst_slack = worst_it->slack;

    for (auto& critical : critical_connections) {
        critical.criticality = HBTACostModel::criticalityFromSlack(
            critical.slack,
            all_worst_slack,
            config.critical_slack_threshold);
    }

    std::sort(critical_connections.begin(),
              critical_connections.end(),
              [](const CriticalConnection& lhs, const CriticalConnection& rhs) {
                  if (lhs.criticality != rhs.criticality) {
                      return lhs.criticality > rhs.criticality;
                  }
                  if (lhs.slack != rhs.slack) {
                      return lhs.slack < rhs.slack;
                  }
                  return lhs.connection->net_name < rhs.connection->net_name;
              });

    size_t selection_limit = critical_connections.size();
    if (config.critical_nets_percent > 0.0 && config.critical_nets_percent < 100.0) {
        const size_t percent_limit = std::max<size_t>(
            1,
            static_cast<size_t>(std::ceil(net_connections.size()
                                          * config.critical_nets_percent / 100.0)));
        selection_limit = std::min(selection_limit, percent_limit);
    }
    if (config.max_critical_nets > 0) {
        selection_limit = std::min(selection_limit, config.max_critical_nets);
    }
    critical_connections.resize(selection_limit);

    const double worst_slack = critical_connections.front().slack;

    // 把 bonding box 内已占用 HBT 的owner net 同步纳入 selected nets。
    std::unordered_map<std::string, const InterLayerConnection*> connection_by_net;
    connection_by_net.reserve(net_connections.size());
    for (const auto& connection : net_connections) {
        if (validateNetConnection(connection)) {
            connection_by_net.emplace(connection.net_name, &connection);
        }
    }

    // candidate_hbts 是 critical net 局部 bonding box 内可参与本轮匹配的 HBT。
    std::vector<const HBTInfo*> candidate_hbts;
    candidate_hbts.reserve(hbt_array_.getAllHBTs().size());
    std::unordered_set<std::string> candidate_hbt_ids;
    std::vector<const InterLayerConnection*> selected_connections;
    selected_connections.reserve(critical_connections.size());
    std::unordered_set<std::string> selected_net_names;

    auto add_selected_connection = [&](const InterLayerConnection* connection) {
        if (connection == nullptr) {
            return false;
        }
        if (!selected_net_names.insert(connection->net_name).second) {
            return false;
        }
        selected_connections.push_back(connection);
        return true;
    };

    for (const auto& critical : critical_connections) {
        add_selected_connection(critical.connection);
    }

    const double hbt_pitch = std::max(hbt_array_.getPitchDBU(), 1.0);
    const double bonding_box_half_size =
        std::max(config.local_bonding_box_radius_in_hbt_pitch, 0.0) * hbt_pitch;

    for (const auto& critical : critical_connections) {
        const HBTInfo* current_hbt = getHBTForNet(critical.connection->net_name);
        const HBTRegion bonding_box =
            makeLocalBondingBox(*critical.connection, current_hbt, bonding_box_half_size);

        for (const auto& hbt : hbt_array_.getAllHBTs()) {
            if (!bonding_box.containsPoint(hbt.location)) {
                continue;
            }

            if (!hbt.used) {
                if (candidate_hbt_ids.insert(hbt.id).second) {
                    candidate_hbts.push_back(&hbt);
                }
                continue;
            }

            const std::string owner_net = getNetForHBT(hbt.id);
            auto owner_it = connection_by_net.find(owner_net);
            if (owner_it == connection_by_net.end()) {
                std::cerr << "WARNING: occupied HBT " << hbt.id
                          << " is inside a timing-driven bonding box but owner net "
                          << owner_net << " is not available for local rematch"
                          << std::endl;
                continue;
            }

            add_selected_connection(owner_it->second);
            if (candidate_hbt_ids.insert(hbt.id).second) {
                candidate_hbts.push_back(&hbt);
            }
        }
    }

    if (candidate_hbts.empty()) {
        std::cout << "Timing-driven HBT reallocation: no local candidate HBTs found"
                  << std::endl;
        return problem;
    }

    std::vector<const HBTInfo*> reduced_candidate_hbts;
    std::vector<std::vector<const HBTInfo*>> candidate_hbts_by_connection;
    reduced_candidate_hbts.reserve(candidate_hbts.size());
    candidate_hbts_by_connection.resize(selected_connections.size());

    std::unordered_set<std::string> reduced_candidate_hbt_ids;
    reduced_candidate_hbt_ids.reserve(candidate_hbts.size());

    for (size_t connection_idx = 0; connection_idx < selected_connections.size();
         ++connection_idx) {
        const auto* selected_connection = selected_connections[connection_idx];
        if (selected_connection == nullptr) {
            continue;
        }

        const HBTInfo* current_hbt = getHBTForNet(selected_connection->net_name);
        const HBTRegion bonding_box =
            makeLocalBondingBox(*selected_connection, current_hbt, bonding_box_half_size);
        auto& local_candidates = candidate_hbts_by_connection[connection_idx];

        for (const auto& hbt : hbt_array_.getAllHBTs()) {
            if (!bonding_box.containsPoint(hbt.location)) {
                continue;
            }

            bool allow_candidate = !hbt.used;
            if (!allow_candidate) {
                const std::string owner_net = getNetForHBT(hbt.id);
                allow_candidate = !owner_net.empty()
                                  && selected_net_names.find(owner_net) != selected_net_names.end();
            }
            if (!allow_candidate) {
                continue;
            }

            local_candidates.push_back(&hbt);
            if (reduced_candidate_hbt_ids.insert(hbt.id).second) {
                reduced_candidate_hbts.push_back(&hbt);
            }
        }
    }

    // 保留旧 HBT，确保压缩后仍允许 no-op assignment 和完整回滚基线。
    std::unordered_set<std::string> reduced_hbt_ids;
    reduced_hbt_ids.reserve(reduced_candidate_hbts.size());
    for (const auto* hbt : reduced_candidate_hbts) {
        if (hbt != nullptr) {
            reduced_hbt_ids.insert(hbt->id);
        }
    }
    for (size_t connection_idx = 0; connection_idx < selected_connections.size(); ++connection_idx) {
        const auto* selected_connection = selected_connections[connection_idx];
        auto old_hbt_it = net_to_hbt_.find(selected_connection->net_name);
        if (old_hbt_it == net_to_hbt_.end()) {
            continue;
        }
        const HBTInfo* old_hbt = hbt_array_.getHBT(old_hbt_it->second);
        if (old_hbt != nullptr && reduced_hbt_ids.insert(old_hbt->id).second) {
            reduced_candidate_hbts.push_back(old_hbt);
        }
        if (old_hbt != nullptr) {
            candidate_hbts_by_connection[connection_idx].push_back(old_hbt);
        }
    }
    std::sort(reduced_candidate_hbts.begin(),
              reduced_candidate_hbts.end(),
              compareHBTCandidatePointers);

    for (auto& per_connection_candidates : candidate_hbts_by_connection) {
        std::sort(per_connection_candidates.begin(),
                  per_connection_candidates.end(),
                  compareHBTCandidatePointers);
        per_connection_candidates.erase(
            std::unique(per_connection_candidates.begin(), per_connection_candidates.end()),
            per_connection_candidates.end());
    }

    problem.connections.reserve(selected_connections.size());
    problem.candidate_indices_by_connection.reserve(selected_connections.size());
    // previous_assignments 记录本轮被局部释放且已经有旧 HBT 的 net，用于失败回滚。
    problem.previous_assignments.reserve(selected_connections.size());

    // 记录旧 HBT 位置，真实 STA 候选边打分由 Timingdp 层完成。
    for (const auto* selected_connection : selected_connections) {
        InterLayerConnection problem_connection = *selected_connection;

        auto old_hbt_it = net_to_hbt_.find(selected_connection->net_name);
        if (old_hbt_it != net_to_hbt_.end()) {
            problem.previous_assignments.emplace_back(selected_connection->net_name,
                                                      old_hbt_it->second);
            const HBTInfo* old_hbt = hbt_array_.getHBT(old_hbt_it->second);
            if (old_hbt) {
                problem_connection.has_previous_hbt_location = true;
                problem_connection.previous_hbt_location = old_hbt->location;
            }
        }

        problem.connections.push_back(problem_connection);
    }

    problem.candidate_hbts.reserve(reduced_candidate_hbts.size());
    for (const auto* hbt : reduced_candidate_hbts) {
        if (hbt == nullptr) {
            continue;
        }
        problem.candidate_hbts.push_back({hbt->id, hbt->location});
    }

    std::unordered_map<std::string, size_t> candidate_index_by_hbt_id;
    candidate_index_by_hbt_id.reserve(problem.candidate_hbts.size());
    for (size_t hbt_idx = 0; hbt_idx < problem.candidate_hbts.size(); ++hbt_idx) {
        candidate_index_by_hbt_id.emplace(problem.candidate_hbts[hbt_idx].hbt_id, hbt_idx);
    }

    for (size_t connection_idx = 0; connection_idx < selected_connections.size(); ++connection_idx) {
        std::vector<size_t> candidate_indices;
        candidate_indices.reserve(candidate_hbts_by_connection[connection_idx].size());
        for (const auto* hbt : candidate_hbts_by_connection[connection_idx]) {
            if (hbt == nullptr) {
                continue;
            }
            const auto index_it = candidate_index_by_hbt_id.find(hbt->id);
            if (index_it != candidate_index_by_hbt_id.end()) {
                candidate_indices.push_back(index_it->second);
            }
        }
        if (candidate_indices.empty()) {
            std::cout << "Timing-driven HBT reallocation: net "
                      << selected_connections[connection_idx]->net_name
                      << " has no local HBT candidates after filtering"
                      << std::endl;
            return problem;
        } else {
            std::sort(candidate_indices.begin(), candidate_indices.end());
            candidate_indices.erase(
                std::unique(candidate_indices.begin(), candidate_indices.end()),
                candidate_indices.end());
        }
        problem.candidate_indices_by_connection.push_back(std::move(candidate_indices));
    }

    problem.worst_slack = worst_slack;

    std::cout << "Timing-driven HBT reallocation: selected "
              << critical_connections.size() << " critical nets, "
              << problem.connections.size() << " local selected nets"
              << ", local_candidate_hbts=" << candidate_hbts.size()
              << ", reduced_candidate_hbts=" << reduced_candidate_hbts.size()
              << ", selected_worst_slack=" << worst_slack
	              << ", bonding_box_radius_pitch="
	              << config.local_bonding_box_radius_in_hbt_pitch << std::endl;

    return problem;
}

HBTAllocator::AssignmentSnapshot HBTAllocator::captureAssignmentSnapshot() const {
    AssignmentSnapshot snapshot;
    snapshot.hbt_array = hbt_array_;
    snapshot.net_to_hbt = net_to_hbt_;
    snapshot.hbt_to_net = hbt_to_net_;
    snapshot.allocation_history = allocation_history_;
    snapshot.pending_active = pending_timing_reallocation_.active;
    snapshot.pending_history_size_before_reassign =
        pending_timing_reallocation_.history_size_before_reassign;
    snapshot.pending_selected_net_names = pending_timing_reallocation_.selected_net_names;
    snapshot.pending_previous_assignments = pending_timing_reallocation_.previous_assignments;
    return snapshot;
}

void HBTAllocator::restoreAssignmentSnapshot(const AssignmentSnapshot& snapshot) {
    hbt_array_ = snapshot.hbt_array;
    net_to_hbt_ = snapshot.net_to_hbt;
    hbt_to_net_ = snapshot.hbt_to_net;
    allocation_history_ = snapshot.allocation_history;
    pending_timing_reallocation_.active = snapshot.pending_active;
    pending_timing_reallocation_.history_size_before_reassign =
        snapshot.pending_history_size_before_reassign;
    pending_timing_reallocation_.selected_net_names = snapshot.pending_selected_net_names;
    pending_timing_reallocation_.previous_assignments =
        snapshot.pending_previous_assignments;
}

bool HBTAllocator::forceAssignHBTToNetForEvaluation(
    const InterLayerConnection& connection,
    const std::string& hbt_id,
    double estimated_wirelength) {
    const auto owner_it = hbt_to_net_.find(hbt_id);
    const bool candidate_owned_by_other =
        owner_it != hbt_to_net_.end() && owner_it->second != connection.net_name;

    if (isNetAllocated(connection.net_name)) {
        releaseHBTForNet(connection.net_name);
    }

    if (!candidate_owned_by_other) {
        if (!hbt_array_.markHBTUsed(hbt_id, connection.net_name)) {
            return false;
        }
        hbt_to_net_[hbt_id] = connection.net_name;
    }

    net_to_hbt_[connection.net_name] = hbt_id;

    NetAllocationResult result;
    result.success = true;
    result.net_name = connection.net_name;
    result.allocated_hbt = hbt_array_.getHBT(hbt_id);
    result.estimated_wirelength = estimated_wirelength;
    result.driver_cell = connection.driver_cell_name;
    result.load_cells = connection.load_cell_names;
    allocation_history_.push_back(result);

    return true;
}

HBTAllocator::BatchNetAllocationResult
HBTAllocator::commitTimingDrivenReallocationWithCostMatrix(
    const TimingDrivenReallocationProblem& problem,
    const std::vector<std::vector<double>>& cost_matrix) {
    BatchNetAllocationResult batch_result;
    batch_result.total_nets = problem.connections.size();

    if (pending_timing_reallocation_.active) {
        std::cerr << "WARNING: previous timing-driven reallocation is still pending; "
                  << "rolling it back before starting a new one" << std::endl;
        rollbackPendingTimingDrivenReallocation();
    }

    if (problem.empty()) {
        return batch_result;
    }

    PendingTimingReallocation snapshot;
    snapshot.active = true;
    snapshot.history_size_before_reassign = allocation_history_.size();
    snapshot.previous_assignments = problem.previous_assignments;
    snapshot.selected_net_names.reserve(problem.connections.size());
    for (const auto& connection : problem.connections) {
        snapshot.selected_net_names.push_back(connection.net_name);
    }

    // 释放 selected nets： critical nets 以及 critical bonding box 内被占用HBT 的 owner nets。
    for (const auto& [net_name, ignored_hbt_id] : problem.previous_assignments) {
        releaseHBTForNet(net_name);
    }

    std::vector<const InterLayerConnection*> connection_ptrs;
    connection_ptrs.reserve(problem.connections.size());
    for (const auto& connection : problem.connections) {
        connection_ptrs.push_back(&connection);
    }

    std::vector<const HBTInfo*> candidate_hbts;
    candidate_hbts.reserve(problem.candidate_hbts.size());
    for (const auto& candidate : problem.candidate_hbts) {
        const auto* hbt = hbt_array_.getHBT(candidate.hbt_id);
        if (hbt == nullptr) {
            batch_result.failed_allocations = batch_result.total_nets;
            for (const auto& connection : problem.connections) {
                NetAllocationResult result;
                result.net_name = connection.net_name;
                result.error_message = "Candidate HBT disappeared before STA-cost commit";
                batch_result.results.push_back(result);
            }
            restoreTimingDrivenReallocationSnapshot(snapshot);
            return batch_result;
        }
        candidate_hbts.push_back(hbt);
    }

    const auto quad_plan =
        HBTQuadTreeAllocator::planAssignmentsWithCostMatrixByComponents(
            connection_ptrs,
            candidate_hbts,
            cost_matrix,
            problem.candidate_indices_by_connection);
    if (!quad_plan.success
        || quad_plan.planned_assignments.size() != problem.connections.size()) {
        for (const auto& connection : problem.connections) {
            NetAllocationResult result;
            result.net_name = connection.net_name;
            result.error_message = quad_plan.error_message.empty()
                                       ? "STA-cost assignment failed"
                                       : quad_plan.error_message;
            batch_result.results.push_back(result);
        }
        batch_result.failed_allocations = batch_result.total_nets;
        restoreTimingDrivenReallocationSnapshot(snapshot);
        return batch_result;
    }

    std::unordered_map<std::string, size_t> candidate_index_by_hbt_id;
    candidate_index_by_hbt_id.reserve(problem.candidate_hbts.size());
    for (size_t hbt_idx = 0; hbt_idx < problem.candidate_hbts.size(); ++hbt_idx) {
        candidate_index_by_hbt_id.emplace(problem.candidate_hbts[hbt_idx].hbt_id, hbt_idx);
    }

    for (size_t net_idx = 0; net_idx < quad_plan.planned_assignments.size(); ++net_idx) {
        const auto& assignment = quad_plan.planned_assignments[net_idx];
        if (assignment.connection == nullptr || assignment.hbt == nullptr) {
            continue;
        }

        if (net_idx >= problem.candidate_indices_by_connection.size()) {
            batch_result.failed_allocations = batch_result.total_nets;
            for (const auto& connection : problem.connections) {
                NetAllocationResult result;
                result.net_name = connection.net_name;
                result.error_message = "STA-cost assignment lost local candidate constraints";
                batch_result.results.push_back(result);
            }
            restoreTimingDrivenReallocationSnapshot(snapshot);
            return batch_result;
        }

        const auto candidate_it = candidate_index_by_hbt_id.find(assignment.hbt->id);
        if (candidate_it == candidate_index_by_hbt_id.end()
            || !std::binary_search(problem.candidate_indices_by_connection[net_idx].begin(),
                                   problem.candidate_indices_by_connection[net_idx].end(),
                                   candidate_it->second)) {
            batch_result.results.clear();
            batch_result.failed_allocations = batch_result.total_nets;
            for (const auto& connection : problem.connections) {
                NetAllocationResult result;
                result.net_name = connection.net_name;
                result.error_message =
                    "STA-cost assignment violated quadtree-local HBT constraints";
                batch_result.results.push_back(result);
            }
            restoreTimingDrivenReallocationSnapshot(snapshot);
            return batch_result;
        }
    }

    for (const auto& assignment : quad_plan.planned_assignments) {
        NetAllocationResult result;
        if (assignment.connection == nullptr || assignment.hbt == nullptr) {
            result.error_message = "STA-cost assignment produced a null entry";
            batch_result.failed_allocations++;
            batch_result.results.push_back(result);
            continue;
        }

        result.net_name = assignment.connection->net_name;
        result.driver_cell = assignment.connection->driver_cell_name;
        result.load_cells = assignment.connection->load_cell_names;
        result.estimated_wirelength = assignment.cost;

        if (hbt_array_.markHBTUsed(assignment.hbt->id, assignment.connection->net_name)) {
            net_to_hbt_[assignment.connection->net_name] = assignment.hbt->id;
            hbt_to_net_[assignment.hbt->id] = assignment.connection->net_name;
            result.success = true;
            result.allocated_hbt = hbt_array_.getHBT(assignment.hbt->id);
            batch_result.successful_allocations++;
            batch_result.total_wirelength += result.estimated_wirelength;
            allocation_history_.push_back(result);
        } else {
            result.error_message = "Failed to commit STA-cost HBT assignment";
            batch_result.failed_allocations++;
        }

        batch_result.results.push_back(result);
    }

    if (batch_result.allSuccessful()) {
        pending_timing_reallocation_ = std::move(snapshot);
        return batch_result;
    }

    std::cerr << "Timing-driven HBT reallocation failed; restoring previous assignments"
              << std::endl;

    restoreTimingDrivenReallocationSnapshot(snapshot);

    return batch_result;
}

bool HBTAllocator::shouldAcceptTimingDrivenReallocation(
    double old_late_wns,
    double old_late_tns,
    double new_late_wns,
    double new_late_tns) {
    if (new_late_wns > old_late_wns) {
        return true;
    }

    if (new_late_wns == old_late_wns && new_late_tns > old_late_tns) {
        return true;
    }

    return false;
}

bool HBTAllocator::hasPendingTimingDrivenReallocation() const {
    return pending_timing_reallocation_.active;
}

void HBTAllocator::acceptPendingTimingDrivenReallocation() {
    clearPendingTimingDrivenReallocation();
}

bool HBTAllocator::rollbackPendingTimingDrivenReallocation() {
    if (!pending_timing_reallocation_.active) {
        return false;
    }

    const auto snapshot = pending_timing_reallocation_;
    const bool restored = restoreTimingDrivenReallocationSnapshot(snapshot);
    clearPendingTimingDrivenReallocation();
    return restored;
}

bool HBTAllocator::releaseHBT(const std::string& hbt_id) {
    // 根据 HBT->net 映射查找当前占用该 HBT 的网络。
    auto it = hbt_to_net_.find(hbt_id);
    if (it == hbt_to_net_.end()) {
        return false;
    }
    
    std::string net_name = it->second;
    
    // 先更新 HBTArray 的占用状态，再同步清理双向映射。
    if (!hbt_array_.markHBTAvailable(hbt_id)) {
        return false;
    }

    net_to_hbt_.erase(net_name);
    hbt_to_net_.erase(hbt_id);
    
    return true;
}

bool HBTAllocator::releaseHBTForNet(const std::string& net_name) {
    // 根据 net->HBT 映射定位待释放 HBT，再复用 releaseHBT 统一清理状态。
    auto it = net_to_hbt_.find(net_name);
    if (it == net_to_hbt_.end()) {
        return false;
    }
    
    std::string hbt_id = it->second;
    
    return releaseHBT(hbt_id);
}

const HBTInfo* HBTAllocator::getHBTForNet(const std::string& net_name) const {
    auto it = net_to_hbt_.find(net_name);
    if (it == net_to_hbt_.end()) {
        return nullptr;
    }
    
    return hbt_array_.getHBT(it->second);
}

std::string HBTAllocator::getNetForHBT(const std::string& hbt_id) const {
    auto it = hbt_to_net_.find(hbt_id);
    if (it == hbt_to_net_.end()) {
        return "";
    }
    
    return it->second;
}

bool HBTAllocator::isNetAllocated(const std::string& net_name) const {
    return net_to_hbt_.find(net_name) != net_to_hbt_.end();
}

HBTAllocator::Statistics HBTAllocator::getStatistics() const {
    Statistics stats;
    
    stats.total_hbts = hbt_array_.getTotalCount();
    stats.used_hbts = stats.total_hbts - hbt_array_.getAvailableCount();
    stats.available_hbts = hbt_array_.getAvailableCount();
    stats.utilization_rate = stats.total_hbts > 0 ? 
        static_cast<double>(stats.used_hbts) / stats.total_hbts : 0.0;
    stats.allocated_nets = net_to_hbt_.size();
    
    return stats;
}

std::string HBTAllocator::exportAllocationsToJSON() const {
    std::stringstream ss;
    ss << "{\n";
    ss << "  \"allocations\": [\n";
    
    bool first = true;
    for (const auto& [net_name, hbt_id] : net_to_hbt_) {
        if (!first) {
            ss << ",\n";
        }
        first = false;
        
        const HBTInfo* hbt = hbt_array_.getHBT(hbt_id);
        if (hbt) {
            ss << "    {\n";
            ss << "      \"net_name\": \"" << net_name << "\",\n";
            ss << "      \"hbt_id\": \"" << hbt_id << "\",\n";
            ss << "      \"hbt_row\": " << hbt->row << ",\n";
            ss << "      \"hbt_col\": " << hbt->col << ",\n";
            ss << "      \"hbt_x\": " << hbt->location.x << ",\n";
            ss << "      \"hbt_y\": " << hbt->location.y << "\n";
            ss << "    }";
        }
    }
    
    ss << "\n  ],\n";
    
    // 导出当前映射与统计信息，便于外部流程记录最终 HBTA 结果。
    auto stats = getStatistics();
    ss << "  \"statistics\": {\n";
    ss << "    \"total_hbts\": " << stats.total_hbts << ",\n";
    ss << "    \"used_hbts\": " << stats.used_hbts << ",\n";
    ss << "    \"available_hbts\": " << stats.available_hbts << ",\n";
    ss << "    \"utilization_rate\": " << stats.utilization_rate << ",\n";
    ss << "    \"allocated_nets\": " << stats.allocated_nets << "\n";
    ss << "  }\n";
    ss << "}\n";
    
    return ss.str();
}

void HBTAllocator::resetAllAllocations() {
    // 重置 HBT 阵列占用状态。
    hbt_array_.resetAllHBTs();

    // 清空双向映射。
    net_to_hbt_.clear();
    hbt_to_net_.clear();

    // 清空分配记录缓存。
    allocation_history_.clear();
    clearPendingTimingDrivenReallocation();
    
    std::cout << "All allocations reset" << std::endl;
}

bool HBTAllocator::validateNetConnection(const InterLayerConnection& connection) const {
    if (connection.net_name.empty()) {
        return false;
    }

    // 校验负载数组长度一致性，确保索引访问安全。
    if (connection.load_cell_names.size() != connection.load_locations.size()) {
        return false;
    }
    if (!connection.load_layers.empty()
        && connection.load_layers.size() != connection.load_cell_names.size()) {
        return false;
    }
    
    if (!connection.hasInterLayerConnection()) {
        return false;
    }

    // (0, 0) 作为合法坐标参与计算。
    return true;
}

void HBTAllocator::clearNetMapping(const std::string& net_name) {
    // 只清理双向映射，不会把 HBTArray 中的 used 状态改回 available。
    auto it = net_to_hbt_.find(net_name);
    if (it != net_to_hbt_.end()) {
        std::string hbt_id = it->second;
        hbt_to_net_.erase(hbt_id);
        net_to_hbt_.erase(it);
    }
}

void HBTAllocator::clearPendingTimingDrivenReallocation() {
    pending_timing_reallocation_ = PendingTimingReallocation{};
}

bool HBTAllocator::restoreTimingDrivenReallocationSnapshot(
    const PendingTimingReallocation& snapshot) {
    if (allocation_history_.size() > snapshot.history_size_before_reassign) {
        allocation_history_.resize(snapshot.history_size_before_reassign);
    }

    bool restored_all = true;

    // 释放 tentative assignment，包括本轮新加入但原先没有 HBT 的 selected net。
    for (const auto& net_name : snapshot.selected_net_names) {
        releaseHBTForNet(net_name);
    }

    // 恢复本轮释放前已有的 assignment。
    for (const auto& [net_name, hbt_id] : snapshot.previous_assignments) {
        if (hbt_array_.markHBTUsed(hbt_id, net_name)) {
            net_to_hbt_[net_name] = hbt_id;
            hbt_to_net_[hbt_id] = net_name;
        } else {
            restored_all = false;
            std::cerr << "WARNING: failed to restore HBT " << hbt_id
                      << " for net " << net_name << std::endl;
        }
    }

    return restored_all;
}

} // 命名空间 hbta
