#pragma once

#include "HBTTypes.h"
#include "HBTArray.h"
#include <cstddef>
#include <string>
#include <vector>
#include <unordered_map>

namespace hbta {

struct HBTQuadTreeAllocatorConfig;

class HBTAllocator {
public:
    HBTAllocator() = default;
    ~HBTAllocator() = default;
    
    // 初始化
    bool initialize(const std::string& hbt_array_json, 
                    int dbu_per_micron = 1000);

    // 分配HBT给跨层网络
    struct NetAllocationResult {
        bool success = false;
        std::string net_name;
        const HBTInfo* allocated_hbt = nullptr;
        double estimated_wirelength = 0.0;
        std::string error_message;
        
        // 连接信息
        std::string driver_cell;
        std::vector<std::string> load_cells;
    };
    
    NetAllocationResult allocateHBTForNet(const InterLayerConnection& net_connection);
    
    // 批量分配
    struct BatchNetAllocationResult {
        int total_nets = 0;
        int successful_allocations = 0;
        int failed_allocations = 0;
        double total_wirelength = 0.0;
        std::vector<NetAllocationResult> results;
        
        // 检查是否所有网络都分配成功
        bool allSuccessful() const {
            return successful_allocations == total_nets;
        }
        
        // 获取失败的网络列表
        std::vector<std::string> getFailedNetworks() const {
            std::vector<std::string> failed;
            for (const auto& result : results) {
                if (!result.success) {
                    failed.push_back(result.net_name);
                }
            }
            return failed;
        }
    };
    
    BatchNetAllocationResult allocateHBTsForNets(
        const std::vector<InterLayerConnection>& net_connections);

    struct TimingDrivenReallocationConfig {
        // 临界 slack 阈值，满足 slack < threshold 的网络参与重分配。
        double critical_slack_threshold = 0.0;
        double critical_nets_percent = 100.0;
        size_t max_critical_nets = 0;
        double local_bonding_box_radius_in_hbt_pitch = 5.0;
    };

    struct TimingDrivenCandidateHBT {
        std::string hbt_id;
        Point location;
    };

    struct TimingDrivenReallocationProblem {
        std::vector<InterLayerConnection> connections;
        std::vector<TimingDrivenCandidateHBT> candidate_hbts;
        std::vector<std::vector<size_t>> candidate_indices_by_connection;
        std::vector<std::pair<std::string, std::string>> previous_assignments;
        double worst_slack = 0.0;

        bool empty() const {
            return connections.empty() || candidate_hbts.empty();
        }
    };

    struct AssignmentSnapshot {
        HBTArray hbt_array;
        std::unordered_map<std::string, std::string> net_to_hbt;
        std::unordered_map<std::string, std::string> hbt_to_net;
        std::vector<NetAllocationResult> allocation_history;
        bool pending_active = false;
        size_t pending_history_size_before_reassign = 0;
        std::vector<std::string> pending_selected_net_names;
        std::vector<std::pair<std::string, std::string>> pending_previous_assignments;
    };

    TimingDrivenReallocationProblem prepareTimingDrivenReallocationProblem(
        const std::vector<InterLayerConnection>& net_connections,
        const std::unordered_map<std::string, double>& slack_by_net) const;

    TimingDrivenReallocationProblem prepareTimingDrivenReallocationProblem(
        const std::vector<InterLayerConnection>& net_connections,
        const std::unordered_map<std::string, double>& slack_by_net,
        const TimingDrivenReallocationConfig& config) const;

    BatchNetAllocationResult commitTimingDrivenReallocationWithCostMatrix(
        const TimingDrivenReallocationProblem& problem,
        const std::vector<std::vector<double>>& cost_matrix);

    AssignmentSnapshot captureAssignmentSnapshot() const;
    void restoreAssignmentSnapshot(const AssignmentSnapshot& snapshot);
    bool forceAssignHBTToNetForEvaluation(
        const InterLayerConnection& connection,
        const std::string& hbt_id,
        double estimated_wirelength = 0.0);

    static bool shouldAcceptTimingDrivenReallocation(
        double old_late_wns,
        double old_late_tns,
        double new_late_wns,
        double new_late_tns);

    bool hasPendingTimingDrivenReallocation() const;
    void acceptPendingTimingDrivenReallocation();
    bool rollbackPendingTimingDrivenReallocation();
    
    // 释放HBT
    bool releaseHBT(const std::string& hbt_id);
    bool releaseHBTForNet(const std::string& net_name);
    
    // 查询
    const HBTInfo* getHBTForNet(const std::string& net_name) const;
    std::string getNetForHBT(const std::string& hbt_id) const;
    
    // 检查网络是否已分配HBT
    bool isNetAllocated(const std::string& net_name) const;
    
    // 统计信息
    struct Statistics {
        int total_hbts = 0;
        int used_hbts = 0;
        int available_hbts = 0;
        double utilization_rate = 0.0;  // 使用率
        int allocated_nets = 0;
    };
    
    Statistics getStatistics() const;
    
    // 导出分配结果到JSON
    std::string exportAllocationsToJSON() const;
    
    // 重置所有分配
    void resetAllAllocations();
    
    // 获取HBT阵列引用（只读）
    const HBTArray& getHBTArray() const { return hbt_array_; }
    
    // 检查连接数量是否超过HBT数量
    bool checkConnectionCount(int connection_count) const {
        return connection_count <= hbt_array_.getAvailableCount();
    }
    
    // 获取最大可分配连接数
    int getMaxAllocatableConnections() const {
        return hbt_array_.getAvailableCount();
    }
    
private:
    HBTArray hbt_array_;
    int dbu_per_micron_ = 1000;
    
    // 网络-HBT映射
    std::unordered_map<std::string, std::string> net_to_hbt_;
    std::unordered_map<std::string, std::string> hbt_to_net_;
    
    // 分配记录缓存
    std::vector<NetAllocationResult> allocation_history_;

    struct PendingTimingReallocation {
        bool active = false;
        size_t history_size_before_reassign = 0;
        std::vector<std::string> selected_net_names;
        std::vector<std::pair<std::string, std::string>> previous_assignments;
    };

    PendingTimingReallocation pending_timing_reallocation_;
    
    // 辅助方法
    BatchNetAllocationResult allocateHBTsForNetsWithConfig(
        const std::vector<InterLayerConnection>& net_connections,
        const HBTQuadTreeAllocatorConfig& quad_config);
    BatchNetAllocationResult allocateHBTsForNetsWithConfigAndCandidates(
        const std::vector<InterLayerConnection>& net_connections,
        const HBTQuadTreeAllocatorConfig& quad_config,
        const std::vector<const HBTInfo*>* candidate_hbt_override);
    bool validateNetConnection(const InterLayerConnection& connection) const;
    void clearNetMapping(const std::string& net_name);
    void clearPendingTimingDrivenReallocation();
    bool restoreTimingDrivenReallocationSnapshot(
        const PendingTimingReallocation& snapshot);
};

} // namespace hbta
