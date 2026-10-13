// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// HBT分配器头文件
#include "HBTAllocator.h"
#include "HBTArray.h"
#include "HBTTypes.h"

// 前向声明
namespace lgista {
    class PtGraph;
}

namespace tdp {
    class Timingdp;
    class Network3D;
    class Net3D;
    class Node3D;
    struct Point;
    struct Pointdp3D;
}

namespace hbta {

/**
 * @class HBTDatabaseInterface
 * @brief 3D数据库与HBT分配器之间的接口类
 * 
 * 负责从3D数据库中提取跨层网络信息，并将其转换为HBT分配器可以处理的格式。
 */
class HBTDatabaseInterface {
public:
    HBTDatabaseInterface();
    ~HBTDatabaseInterface() = default;
    
    // 初始化接口
    bool initialize(tdp::Timingdp* timingdp, const std::string& hbt_array_json, int dbu_per_micron = 1000);
    
    // 提取跨层连接
    std::vector<InterLayerConnection> extractInterLayerConnections();
    
    // 获取分配器
    HBTAllocator& getAllocator() { return allocator_; }
    const HBTAllocator& getAllocator() const { return allocator_; }
    
    // 统计信息
    struct ExtractionStats {
        int total_nets = 0;
        int inter_layer_nets = 0;
        int bottom_layer_cells = 0;
        int upper_layer_cells = 0;
        int failed_conversions = 0;
        
        void print() const;
    };
    
    ExtractionStats getExtractionStats() const { return extraction_stats_; }
    void syncHBTAssignmentToNetwork3D(
        const InterLayerConnection& connection);
    void syncHBTAssignmentsToNetwork3D(
        const std::vector<InterLayerConnection>& connections);
    bool runTimingDrivenReallocation(
        const std::vector<InterLayerConnection>& connections);
    tdp::Point convertToTDPPoint(const Point& point) const;
    int64_t convertToTDPDBU(double value) const;
    
    // 重置
    void reset();
    
private:
    // 私有辅助方法
    InterLayerConnection convertNetToConnection(tdp::Net3D* net);
    Point calculateSuggestedHBTLocation(tdp::Net3D* net);
    std::unordered_map<std::string, double> collectNetSlacks(
        const std::vector<InterLayerConnection>& connections) const;
    std::optional<double> collectNetSlack(
        const InterLayerConnection& connection) const;
    std::optional<double> collectNetSlackFromPtGraph(
        const InterLayerConnection& connection,
        lgista::PtGraph* pt_graph) const;
    bool isInterLayerNetwork(tdp::Net3D* net) const;
    std::string getCellLayer(tdp::Node3D* cell) const;
    Point convertToHBTAPoint(const tdp::Pointdp3D& point) const;
    
    // 成员变量
    tdp::Timingdp* timingdp_ = nullptr;
    HBTAllocator allocator_;
    std::string hbt_array_json_path_;
    int dbu_per_micron_ = 1000;
    
    // 统计信息
    ExtractionStats extraction_stats_;
    
    // 缓存
    mutable std::map<std::string, std::string> cell_layer_cache_;
};

} // namespace hbta
