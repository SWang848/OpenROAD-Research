// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Master3D.h"
#include "Common3D.h"
#include "Architecture3D.h"

#include <algorithm>

namespace tdp {

// ============================================================================
// MasterEdge3D 实现
// ============================================================================

MasterEdge3D::MasterEdge3D(unsigned int type, const odb::Rect& box)
    : edge_type_idx_(type), bbox_(box) {
}

unsigned int MasterEdge3D::getEdgeType() const {
    return edge_type_idx_;
}

const odb::Rect& MasterEdge3D::getBBox() const {
    return bbox_;
}

// ============================================================================
// Master3D 实现
// ============================================================================

Master3D::Master3D() {
    // 默认构造函数
}

Master3D::~Master3D() {
    // 析构函数
}

bool Master3D::isMultiRow() const {
    return is_multi_row_;
}

bool Master3D::isMultiHeight() const {
    return is_multi_height_;
}

const std::vector<MasterEdge3D>& Master3D::getEdges() const {
    return edges_;
}

odb::Rect Master3D::getBBox() const {
    return boundary_box_;
}

int Master3D::getBottomPowerType() const {
    return bottom_pwr_;
}

int Master3D::getTopPowerType() const {
    return top_pwr_;
}

void Master3D::setMultiRow(bool multi_row) {
    is_multi_row_ = multi_row;
}

void Master3D::setMultiHeight(bool multi_height) {
    is_multi_height_ = multi_height;
}

void Master3D::addEdge(const MasterEdge3D& edge) {
    edges_.push_back(edge);
}

void Master3D::clearEdges() {
    edges_.clear();
}

void Master3D::setBBox(const odb::Rect& box) {
    boundary_box_ = box;
}

void Master3D::setBottomPowerType(int bottom_pwr) {
    bottom_pwr_ = bottom_pwr;
}

void Master3D::setTopPowerType(int top_pwr) {
    top_pwr_ = top_pwr;
}

void Master3D::setDbMaster(odb::dbMaster* db_master) {
    db_master_ = db_master;
    
    if (db_master_) {
        // 从db_master获取边界框
        // 注意：odb::dbMaster的getBBox()方法可能不存在，使用其他方式获取边界框
        // 尝试获取第一个obstruction box作为边界框
        odb::dbSet<odb::dbBox> obstructions = db_master_->getObstructions();
        if (!obstructions.empty()) {
            odb::dbBox* bbox = *obstructions.begin();
            boundary_box_ = odb::Rect(bbox->xMin(), bbox->yMin(), bbox->xMax(), bbox->yMax());
        } else {
            // 如果没有obstruction，使用宽度和高度创建边界框
            boundary_box_ = odb::Rect(0, 0, db_master_->getWidth(), db_master_->getHeight());
        }
        
        // 设置多行属性
        // 检查master类型是否为多行（宏单元通常为多行）
        odb::dbMasterType master_type = db_master_->getType();
        is_multi_row_ = (master_type == odb::dbMasterType::BLOCK) || 
                        (master_type == odb::dbMasterType::PAD) ||
                        (db_master_->getHeight() > 0);  // 简化判断
        
        // 设置多高属性
        // 获取高度并检查是否大于标准行高
        DBU height = db_master_->getHeight();
        const DBU STANDARD_ROW_HEIGHT = 100;  // 需要从配置获取
        is_multi_height_ = height > STANDARD_ROW_HEIGHT;
    }
}

odb::dbMaster* Master3D::getDbMaster() const {
    return db_master_;
}

std::string Master3D::getName() const {
    if (db_master_) {
        return db_master_->getName();
    }
    return "";
}

DbuX Master3D::getWidth() const {
    if (db_master_) {
        return DbuX{db_master_->getWidth()};
    }
    return DbuX{0};
}

DbuY Master3D::getHeight() const {
    if (db_master_) {
        return DbuY{db_master_->getHeight()};
    }
    return DbuY{0};
}

bool Master3D::canBePlacedOnLayer(const std::string& layer) const {
    if (!db_master_) {
        return false;
    }
    // 统一依据标准单元名后缀判定层归属。
    Architecture3D::LayerType lt = Architecture3D::layerOfMaster(db_master_);
    if (layer == "bottom") {
        return lt == Architecture3D::LAYER_BOTTOM;
    } else if (layer == "upper") {
        return lt == Architecture3D::LAYER_UPPER;
    }

    return true;
}

bool Master3D::requiresHBTConnection() const {
    if (!db_master_) {
        return false;
    }
    
    // 检查是否为混合工艺单元
    std::string master_name = db_master_->getName();
    std::transform(master_name.begin(), master_name.end(), master_name.begin(), ::tolower);
    
    return master_name.find("hbt") != std::string::npos ||
           master_name.find("hybrid") != std::string::npos ||
           master_name.find("mixed") != std::string::npos;
}

} // namespace tdp