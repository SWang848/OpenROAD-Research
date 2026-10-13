// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Network3D.h"
#include "Architecture3D.h"
#include "Common3D.h"

#include <algorithm>
#include <cmath>

namespace tdp {

Node3D::Node3D() {
    // 默认构造函数
    db_inst = nullptr;
    db_bterm_ = nullptr;
    name = "";
    layer = "bottom";
    is_fixed = false;
    is_movable = true;
    arrival_time = 0.0;
    required_time = 0.0;
    slack = 0.0;
    group_id_ = -1;
    group_ = nullptr;
    orient_ = odb::dbOrientType::R0;
    symmetry_ = 0;
    bottom_power_ = 0;
    top_power_ = 0;
    site_ = nullptr;
    type_ = UNKNOWN;
    width_ = DbuX{0};
    height_ = DbuY{0};
    id_ = 0;
}

Node3D::~Node3D() {
    // 析构函数
}

// 引脚管理方法实现（新增）
void Node3D::addPin(Pin3D* pin) {
    if (!pin) {
        return;
    }
    
    // 检查引脚是否已存在
    auto it = std::find(pins_.begin(), pins_.end(), pin);
    if (it == pins_.end()) {
        pins_.push_back(pin);
    }
}

void Node3D::removePin(Pin3D* pin) {
    if (!pin) {
        return;
    }
    
    auto it = std::find(pins_.begin(), pins_.end(), pin);
    if (it != pins_.end()) {
        pins_.erase(it);
    }
}

bool Node3D::isPlaced() const {
    // 更合理的实现：检查是否被标记为已放置
    // 同时检查位置是否在合法范围内
    if (is_fixed) {
        return true;  // 固定单元总是被认为是已放置的
    }
    
    // 检查位置是否有效（非零或已在合法范围内）
    // 这里我们假设如果位置被显式设置过，就是已放置的
    return location.x.v != 0 || location.y.v != 0 || location.z.v != 0;
}

bool Node3D::isMultiHeight() const {
    // 检查是否为多高单元
    if (!db_inst || !db_inst->getMaster()) {
        return false;
    }
    
    // 从master获取高度信息
    odb::dbMaster* master = db_inst->getMaster();
    
    // 简化实现：检查高度是否大于标准行高
    // 在实际实现中，需要从Architecture3D获取标准行高
    DBU height = master->getHeight();
    
    // 假设标准行高为100 DBU（需要从Architecture3D获取实际值）
    const DBU STANDARD_ROW_HEIGHT = 100;
    
    return height > STANDARD_ROW_HEIGHT;
}

int Node3D::getHeightInRows(const Architecture3D* arch) const {
    if (!arch || !db_inst || !db_inst->getMaster()) {
        return 1;
    }
    
    // 获取单元高度
    DBU cell_height = db_inst->getMaster()->getHeight();
    
    // 获取标准行高（从架构）
    DBU row_height = arch->getRowHeight();
    if (row_height == 0) {
        return 1;  // 避免除以零
    }
    
    // 计算需要多少行
    int rows_needed = static_cast<int>(std::ceil(static_cast<double>(cell_height) / row_height));
    
    return std::max(1, rows_needed);
}

bool Node3D::canMoveToLayer(const std::string& target_layer, const Architecture3D* arch) const {
    // 检查是否可以移动到目标层
    if (!is_movable || is_fixed) {
        return false;  // 不可移动或固定的单元不能移动
    }
    
    // 检查目标层是否有效
    if (target_layer != "bottom" && target_layer != "upper") {
        return false;
    }
    
    // 检查当前层和目标层是否相同
    if (layer == target_layer) {
        return true;  // 已经在目标层
    }
    
    // 如果有架构信息，可以检查更多约束
    if (arch) {
        // 检查单元类型是否允许跨层移动
        // 例如，某些特殊单元可能只能留在特定层
        if (db_inst && db_inst->getMaster()) {
            std::string master_name = db_inst->getMaster()->getName();
            // 这里可以添加特定master的约束检查
        }
    }
    
    return true;
}

Pointdp3D Node3D::getLegalLocation(const Architecture3D* arch) const {
    if (!arch) {
        return location;  // 没有架构信息，返回当前位置
    }
    
    // 获取当前层的合法位置
    Architecture3D::LayerType layer_type = (layer == "bottom") ? 
                                          Architecture3D::LAYER_BOTTOM : 
                                          Architecture3D::LAYER_UPPER;
    
    // 获取当前层的边界
    DbuX min_x = arch->getMinX(layer_type);
    DbuX max_x = arch->getMaxX(layer_type);
    DbuY min_y = arch->getMinY(layer_type);
    DbuY max_y = arch->getMaxY(layer_type);
    
    // 计算合法位置（确保在边界内）
    DbuX legal_x = std::max(min_x, std::min(location.x, max_x));
    DbuY legal_y = std::max(min_y, std::min(location.y, max_y));
    DbuZ legal_z = location.z;  // Z坐标保持不变
    
    return Pointdp3D(legal_x, legal_y, legal_z);
}

// ============================================================================
// 新增方法实现
// ============================================================================

odb::dbOrientType Node3D::getOrient() const {
    if (db_inst) {
        return db_inst->getOrient();
    }
    return orient_;
}

void Node3D::setOrient(const odb::dbOrientType& orient) {
    orient_ = orient;
    if (db_inst) {
        db_inst->setOrient(orient);
    }
}

DbuX Node3D::getWidth() const {
    // 如果是Terminal或Filler，使用存储的宽度
    if (type_ == TERMINAL || type_ == FILLER) {
        return width_;
    }
    if (db_inst && db_inst->getMaster()) {
        return DbuX{db_inst->getMaster()->getWidth()};
    }
    return DbuX{0};
}

DbuY Node3D::getHeight() const {
    // 如果是Terminal或Filler，使用存储的高度
    if (type_ == TERMINAL || type_ == FILLER) {
        return height_;
    }
    if (db_inst && db_inst->getMaster()) {
        return DbuY{db_inst->getMaster()->getHeight()};
    }
    return DbuY{0};
}

int Node3D::getBottomPower() const {
    return bottom_power_;
}

int Node3D::getTopPower() const {
    return top_power_;
}

void Node3D::setBottomPower(int power) {
    bottom_power_ = power;
}

void Node3D::setTopPower(int power) {
    top_power_ = power;
}

unsigned Node3D::getSymmetry() const {
    return symmetry_;
}

void Node3D::setSymmetry(unsigned symmetry) {
    symmetry_ = symmetry;
}

odb::dbSite* Node3D::getSite() const {
    return site_;
}

void Node3D::setSite(odb::dbSite* site) {
    site_ = site;
    
    // 从site获取对称性信息
    if (site_) {
        symmetry_ = 0;
        if (site_->getSymmetryX()) symmetry_ |= 0x01;
        if (site_->getSymmetryY()) symmetry_ |= 0x02;
        if (site_->getSymmetryR90()) symmetry_ |= 0x04;
    }
}

DbuX Node3D::siteWidth() const {
    if (site_) {
        return DbuX{site_->getWidth()};
    }
    return DbuX{0};
}

bool Node3D::isHybrid() const {
    // 检查是否为混合工艺单元
    if (!db_inst || !db_inst->getMaster()) {
        return false;
    }
    
    // 简化实现：检查master名称是否包含混合工艺标识
    std::string master_name = db_inst->getMaster()->getName();
    std::transform(master_name.begin(), master_name.end(), master_name.begin(), ::tolower);
    
    return master_name.find("hybrid") != std::string::npos ||
           master_name.find("mixed") != std::string::npos;
}

bool Node3D::isHybridParent() const {
    // 检查是否为混合工艺父单元
    // 简化实现：返回false
    return false;
}

int64_t Node3D::area() const {
    DbuX width = getWidth();
    DbuY height = getHeight();
    return static_cast<int64_t>(width.v) * static_cast<int64_t>(height.v);
}

bool Node3D::isTerminal() const {
    // 首先检查类型
    if (type_ == TERMINAL) {
        return true;
    }
    // 检查是否有BTerm
    if (db_bterm_ != nullptr) {
        return true;
    }
    // 检查是否为终端（IO单元）
    if (!db_inst || !db_inst->getMaster()) {
        return false;
    }
    
    odb::dbMaster* master = db_inst->getMaster();
    return master->isPad() || master->isCover();
}

bool Node3D::isStdCell() const {
    // 检查是否为标准单元
    if (!db_inst || !db_inst->getMaster()) {
        return false;
    }
    
    odb::dbMaster* master = db_inst->getMaster();
    return master->isCore() || master->isCoreAutoPlaceable();
}

bool Node3D::isFiller() const {
    return type_ == FILLER;
}

bool Node3D::isBlock() const {
    // 检查是否为宏单元
    if (!db_inst || !db_inst->getMaster()) {
        return false;
    }
    
    odb::dbMaster* master = db_inst->getMaster();
    return master->isBlock();
}

odb::Rect Node3D::getBBox() const {
    // 对于Terminal，从BTerm获取
    if (db_bterm_) {
        return db_bterm_->getBBox();
    }
    // 对于Cell，从dbInst获取
    if (db_inst) {
        odb::dbBox* bbox = db_inst->getBBox();
        if (bbox) {
            return odb::Rect(bbox->xMin(), bbox->yMin(), bbox->xMax(), bbox->yMax());
        }
    }
    
    // 如果没有bbox，根据位置和尺寸计算
    DBU x1 = location.x.v;
    DBU y1 = location.y.v;
    DBU x2 = x1 + getWidth().v;
    DBU y2 = y1 + getHeight().v;
    
    return odb::Rect(x1, y1, x2, y2);
}

bool Node3D::adjustCurrOrient(const odb::dbOrientType& newOrient) {
    // 调整当前方向
    // 简化实现：直接设置新方向
    setOrient(newOrient);
    return true;
}

void Node3D::addConnectedNet(odb::dbNet* db_net)
{
    if (!db_net) return;
    if (std::find(connected_nets.begin(), connected_nets.end(), db_net)
        == connected_nets.end()) {
        connected_nets.push_back(db_net);
    }
}

} // namespace tdp
