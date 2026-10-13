// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <string>
#include <vector>

#include "Common3D.h"
#include "Coordinates.h"
#include "odb/db.h"

// Forward declarations
namespace odb {
    class dbInst;
    class dbNet;
    class dbSite;
    class dbBTerm;
}

namespace tdp {

// Forward declarations
class Architecture3D;
class Pin3D;  

/**
 * @class Node3D
 * @brief 3D cell representation similar to DPL's Node class
 * 
 * Represents a cell in 3D space with layer information.
 */
class Node3D
{
public:
    // 节点类型枚举（参考DPL的Node::Type）
    enum Type {
        UNKNOWN,
        CELL,       // 标准单元
        TERMINAL,   // IO Port (BTerm)
        MACROCELL,  // 宏单元
        FILLER      // 填充单元
    };
    
    Node3D();
    ~Node3D();
    
    // 基本信息
    odb::dbInst* db_inst = nullptr;
    odb::dbBTerm* db_bterm_ = nullptr;  // IO Port
    std::string name;
    
    // 3D位置信息
    Pointdp3D location;
    
    // 层信息
    std::string layer;  // "bottom" or "upper"
    bool is_fixed = false;
    bool is_movable = true;
    
    // 时序信息
    double arrival_time = 0.0;
    double required_time = 0.0;
    double slack = 0.0;
    
    // 连接的网络
    std::vector<odb::dbNet*> connected_nets;
    void addConnectedNet(odb::dbNet* db_net);
    const std::vector<odb::dbNet*>& getConnectedNets() const { return connected_nets; }
    
    // 引脚管理（新增）
    const std::vector<Pin3D*>& getPins() const { return pins_; }
    void addPin(Pin3D* pin);
    void removePin(Pin3D* pin);
    int getNumPins() const { return pins_.size(); }
    // 基本方法
    bool isPlaced() const;
    bool isMultiHeight() const;
    int getHeightInRows(const Architecture3D* arch) const;
    
    // 3D特定方法
    bool canMoveToLayer(const std::string& target_layer, const Architecture3D* arch) const;
    Pointdp3D getLegalLocation(const Architecture3D* arch) const;
    
    // 组管理
    void setGroupId(int group_id) { group_id_ = group_id; }
    int getGroupId() const { return group_id_; }
    void setGroup(void* group) { group_ = group; }
    void* getGroup() const { return group_; }
    // ID 相关方法
    int getId() const { return id_; }
    void setId(int id) { id_ = id; }
    
    // 新增方法（从DPL的Node类借鉴）
    odb::dbOrientType getOrient() const;
    void setOrient(const odb::dbOrientType& orient);
    
    DbuX getWidth() const;
    DbuY getHeight() const;
    
    int getBottomPower() const;
    int getTopPower() const;
    void setBottomPower(int power);
    void setTopPower(int power);
    
    unsigned getSymmetry() const;
    void setSymmetry(unsigned symmetry);
    
    odb::dbSite* getSite() const;
    void setSite(odb::dbSite* site);
    
    DbuX siteWidth() const;
    
    bool isHybrid() const;
    bool isHybridParent() const;
    
    int64_t area() const;
    
    bool isTerminal() const;
    bool isStdCell() const;
    bool isBlock() const;
    bool isFiller() const;
    
    odb::Rect getBBox() const;
    
    bool adjustCurrOrient(const odb::dbOrientType& newOrient);
    
    // Type相关方法（参考DPL）
    Type getType() const { return type_; }
    void setType(Type type) { type_ = type; }
    
    // BTerm相关方法（参考DPL）
    void setBTerm(odb::dbBTerm* bterm) { db_bterm_ = bterm; }
    odb::dbBTerm* getBTerm() const { return db_bterm_; }
    
    // 宽度和高度的setter（用于Terminal）
    void setWidth(DbuX width) { width_ = width; }
    void setHeight(DbuY height) { height_ = height; }

    const std::string& getLayer() const { return layer; }
    const Pointdp3D& getLocation() const { return location; }
    void setLocation(const Pointdp3D& loc) { location = loc; }

    
private:
    // 引脚列表
    std::vector<Pin3D*> pins_;
    // 组信息
    int group_id_ = -1;
    int id_ = 0;
    void* group_ = nullptr;
    
    // 新增私有成员
    odb::dbOrientType orient_ = odb::dbOrientType::R0;
    unsigned symmetry_ = 0;
    int bottom_power_ = 0;
    int top_power_ = 0;
    odb::dbSite* site_ = nullptr;
    
    // 节点类型
    Type type_ = UNKNOWN;
    
    // 宽度和高度（用于Terminal，因为Terminal没有dbInst）
    DbuX width_{0};
    DbuY height_{0};
    
    // 禁止拷贝
    Node3D(const Node3D&) = delete;
    Node3D& operator=(const Node3D&) = delete;
};

} // namespace tdp
