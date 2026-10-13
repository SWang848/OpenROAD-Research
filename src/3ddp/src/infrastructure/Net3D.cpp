// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Network3D.h"
#include "Common3D.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace tdp {

Net3D::Net3D() {
    // 默认构造函数
    db_net = nullptr;
    name = "";
    timing_criticality = 0.0;
    slack = 0.0;
    id_ = -1;
}

Net3D::~Net3D() {
    // 析构函数
}

double Net3D::calculate3DHPWL() const {

    return 0.0;
}

bool Net3D::hasInterLayerConnections() const {
    // 检查是否有跨层连接
    if (connected_cells.empty()) return false;
    
    std::string first_layer = connected_cells[0]->layer;
    for (size_t i = 1; i < connected_cells.size(); i++) {
        if (connected_cells[i]->layer != first_layer) {
            return true;
        }
    }
    return false;
}

std::vector<Pointdp3D> Net3D::getPinLocations3D() const {
    std::vector<Pointdp3D> locations;
    for (auto* cell : connected_cells) {
        locations.push_back(cell->location);
    }
    return locations;
}

double Net3D::calculateInterLayerDistance() const {
    // 计算跨层连接的平均距离
    if (!hasInterLayerConnections()) {
        return 0.0;
    }
    
    double total_distance = 0.0;
    int cross_layer_pairs = 0;
    
    // 计算所有跨层单元对之间的距离
    for (size_t i = 0; i < connected_cells.size(); i++) {
        for (size_t j = i + 1; j < connected_cells.size(); j++) {
            if (connected_cells[i]->layer != connected_cells[j]->layer) {
                // 计算3D曼哈顿距离
                double distance = manhattanDistance3D(
                    connected_cells[i]->location,
                    connected_cells[j]->location
                );
                total_distance += distance;
                cross_layer_pairs++;
            }
        }
    }
    
    if (cross_layer_pairs == 0) {
        return 0.0;
    }
    
    return total_distance / cross_layer_pairs;
}

bool Net3D::checkHBTAlignment() const {
    // 检查所有HBT连接是否对齐
    for (const auto& hbt : hbt_connections) {
        if (!hbt.aligned) {
            // 检查HBT间距是否符合要求
            if (hbt.bottom_cell && hbt.upper_cell) {
                Pointdp3D bottom_loc = hbt.bottom_cell->location;
                Pointdp3D upper_loc = hbt.upper_cell->location;
                
                // 检查XY对齐（Z方向已经通过层确定）
                DBU x_diff = std::abs(bottom_loc.x.v - upper_loc.x.v);
                DBU y_diff = std::abs(bottom_loc.y.v - upper_loc.y.v);
                
                // 需要符合HBT间距要求
                if (hbt.pitch > 0) {
                    if (x_diff % hbt.pitch != 0 || y_diff % hbt.pitch != 0) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// ============================================================================
// 新增方法实现
// ============================================================================

int Net3D::getId() const {
    return id_;
}

void Net3D::setId(int id) {
    id_ = id;
}

int Net3D::getNumPins() const {
    return pins_.size();
}

const std::vector<Pin3D*>& Net3D::getPins() const {
    return pins_;
}

void Net3D::addPin(Pin3D* pin) {
    pins_.push_back(pin);
}

void Net3D::removePin(Pin3D* pin) {
    auto it = std::find(pins_.begin(), pins_.end(), pin);
    if (it != pins_.end()) {
        pins_.erase(it);
    }
}

uint64_t Net3D::hpwl() const {
    // 计算2D HPWL（用于向后兼容）
    if (connected_cells.empty()) {
        return 0;
    }
    
    DBU x_min = std::numeric_limits<DBU>::max();
    DBU x_max = std::numeric_limits<DBU>::lowest();
    DBU y_min = std::numeric_limits<DBU>::max();
    DBU y_max = std::numeric_limits<DBU>::lowest();
    
    for (const auto* cell : connected_cells) {
        x_min = std::min(x_min, cell->location.x.v);
        x_max = std::max(x_max, cell->location.x.v);
        y_min = std::min(y_min, cell->location.y.v);
        y_max = std::max(y_max, cell->location.y.v);
    }
    
    return (x_max - x_min) + (y_max - y_min);
}

double Net3D::calculateTimingCriticality() const {
    // 计算时序关键性
    // 简化实现：基于slack计算
    if (slack >= 0) {
        return 0.0;  // 正slack，非关键
    } else {
        return -slack;  // 负slack越大，关键性越高
    }
}

void Net3D::updateTimingInfo(double new_slack) {
    slack = new_slack;
    timing_criticality = calculateTimingCriticality();
}

std::vector<Node3D*> Net3D::getCellsOnLayer(const std::string& layer_name) const {
    std::vector<Node3D*> layer_cells;
    for (auto* cell : connected_cells) {
        if (cell->layer == layer_name) {
            layer_cells.push_back(cell);
        }
    }
    return layer_cells;
}

bool Net3D::hasHBTConnections() const {
    return !hbt_connections.empty();
}

void Net3D::addHBTConnection(const HBTConnection& connection) {
    hbt_connections.push_back(connection);
}

void Net3D::removeHBTConnection(const HBTConnection& connection) {
    auto it = std::find(hbt_connections.begin(), hbt_connections.end(), connection);
    if (it != hbt_connections.end()) {
        hbt_connections.erase(it);
    }
}

void Net3D::updateHBTAlignment() {
    // 更新所有HBT连接的对齐状态
    for (auto& hbt : hbt_connections) {
        hbt.aligned = checkSingleHBTAlignment(hbt);
    }
}

bool Net3D::checkSingleHBTAlignment(const HBTConnection& hbt) const {
    if (!hbt.bottom_cell || !hbt.upper_cell) {
        return false;
    }
    
    Pointdp3D bottom_loc = hbt.bottom_cell->location;
    Pointdp3D upper_loc = hbt.upper_cell->location;
    
    // 检查XY对齐
    DBU x_diff = std::abs(bottom_loc.x.v - upper_loc.x.v);
    DBU y_diff = std::abs(bottom_loc.y.v - upper_loc.y.v);
    
    // 需要符合HBT间距要求
    if (hbt.pitch > 0) {
        return (x_diff % hbt.pitch == 0) && (y_diff % hbt.pitch == 0);
    }
    
    // 如果没有指定pitch，检查是否完全对齐
    return (x_diff == 0) && (y_diff == 0);
}

double Net3D::calculateAlignmentCost() const {
    // 计算对齐成本（未对齐的HBT连接越多，成本越高）
    double cost = 0.0;
    for (const auto& hbt : hbt_connections) {
        if (!hbt.aligned) {
            cost += 1.0;  // 每个未对齐的HBT连接增加1.0成本
        }
    }
    return cost;
}

std::vector<Point> Net3D::getHBTLocations() const {
    std::vector<Point> locations;
    for (const auto& hbt : hbt_connections) {
        locations.push_back(hbt.location);
    }
    return locations;
}

void Net3D::printInfo() const {
    std::cout << "Net: " << name << std::endl;
    std::cout << "  ID: " << id_ << std::endl;
    std::cout << "  Cells: " << connected_cells.size() << std::endl;
    std::cout << "  Pins: " << pins_.size() << std::endl;
    std::cout << "  Has inter-layer connections: " << (hasInterLayerConnections() ? "yes" : "no") << std::endl;
    std::cout << "  HBT connections: " << hbt_connections.size() << std::endl;
    std::cout << "  HBT aligned: " << (checkHBTAlignment() ? "yes" : "no") << std::endl;
    std::cout << "  Slack: " << slack << std::endl;
    std::cout << "  Timing criticality: " << timing_criticality << std::endl;
}

} // namespace tdp
