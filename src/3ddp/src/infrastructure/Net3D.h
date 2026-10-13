// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <string>
#include <vector>

#include "Common3D.h"
#include "Coordinates.h"
#include "odb/db.h"

// Forward declarations
namespace tdp {
    class Node3D;
    class Pin3D;
}

namespace tdp {

/**
 * @class Net3D
 * @brief 3D net representation with inter-layer connections
 * 
 * Represents a net that can connect cells across different layers.
 */
class Net3D
{
public:
    Net3D();
    ~Net3D();
    
    // 基本信息
    odb::dbNet* db_net = nullptr;
    std::string name;
    
    // 连接的单元
    std::vector<Node3D*> connected_cells;
    
    // 层间连接点 (HBTs)
    struct HBTConnection {
        Point location;                 // HBT位置
        Node3D* bottom_cell = nullptr;  // 下层单元
        Node3D* upper_cell = nullptr;   // 上层单元
        DBU pitch;                      // HBT间距
        bool aligned = false;           // 是否对齐
        
        bool operator==(const HBTConnection& other) const {
            return location.x.v == other.location.x.v &&
                   location.y.v == other.location.y.v &&
                   bottom_cell == other.bottom_cell &&
                   upper_cell == other.upper_cell &&
                   pitch == other.pitch &&
                   aligned == other.aligned;
        }
    };
    
    std::vector<HBTConnection> hbt_connections;
    
    // 时序关键性
    double timing_criticality = 0.0;
    double slack = 0.0;
    
    // 基本方法
    double calculate3DHPWL() const;
    bool hasInterLayerConnections() const;
    std::vector<Pointdp3D> getPinLocations3D() const;
    
    // 3D特定方法
    double calculateInterLayerDistance() const;
    bool checkHBTAlignment() const;
    
    // 新增方法
    int getId() const;
    void setId(int id);
    
    int getNumPins() const;
    const std::vector<Pin3D*>& getPins() const;
    void addPin(Pin3D* pin);
    void removePin(Pin3D* pin);
    
    uint64_t hpwl() const;
    
    double calculateTimingCriticality() const;
    void updateTimingInfo(double new_slack);
    
    std::vector<Node3D*> getCellsOnLayer(const std::string& layer_name) const;
    
    bool hasHBTConnections() const;
    void addHBTConnection(const HBTConnection& connection);
    void removeHBTConnection(const HBTConnection& connection);
    void updateHBTAlignment();
    bool checkSingleHBTAlignment(const HBTConnection& hbt) const;
    
    double calculateAlignmentCost() const;
    
    std::vector<Point> getHBTLocations() const;
    
    void printInfo() const;
    
private:
    int id_ = -1;
    std::vector<Pin3D*> pins_;
    
    // 禁止拷贝
    Net3D(const Net3D&) = delete;
    Net3D& operator=(const Net3D&) = delete;
};

} // namespace tdp