// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Architecture3D.h"
#include "Network3D.h"
#include "Common3D.h"

#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"

namespace tdp {

// ========================================================================
// Architecture3D实现
// ========================================================================

Architecture3D::Architecture3D()
{
    // 初始化所有层
    for (int i = 0; i < LAYER_COUNT; i++) {
        LayerType layer = static_cast<LayerType>(i);
        layers_[layer] = LayerInfo();
    }
}

Architecture3D::~Architecture3D()
{
    clear();
}

// ========================================================================
// 从数据库创建Row（类似DPL的createArchitecture）
// ========================================================================
/*
Architecture3D::Row3D* Architecture3D::createRowFromDatabase(LayerType layer, odb::dbRow* db_row)
{
    if (!db_row || !db_row->getSite()) {
        return nullptr;
    }
    
    odb::dbSite* site = db_row->getSite();
    odb::Point origin = db_row->getOrigin();
    
    // 创建Row3D对象
    Row3D* row = createAndAddRow(layer);
    if (!row) {
        return nullptr;
    }
    
    // 设置参数（类似DPL的createArchitecture）
    row->setSubRowOrigin(DbuX{origin.x()});
    row->setBottom(DbuY{origin.y()});
    row->setSiteSpacing(DbuX{db_row->getSpacing()});
    row->setNumSites(db_row->getSiteCount());
    row->setSiteWidth(DbuX{site->getWidth()});
    row->setHeight(DbuY{site->getHeight()});
    
    // 设置电源（默认值）
    row->setBottomPower(Row3D::POWER_UNKNOWN);
    row->setTopPower(Row3D::POWER_UNKNOWN);
    
    // 设置对称性（从site获取）
    unsigned symmetry = 0;
    if (site->getSymmetryX()) {
        symmetry |= 0x1;  // X对称
    }
    if (site->getSymmetryY()) {
        symmetry |= 0x2;  // Y对称
    }
    row->setSymmetry(symmetry);
    
    // 设置方向（从row获取）
    row->setOrient(db_row->getOrient());
    
    return row;
}
*/
// ========================================================================
// 其他方法的实现
// ========================================================================

const char* Architecture3D::layerToString(LayerType layer)
{
    switch (layer) {
        case LAYER_BOTTOM: return "bottom";
        case LAYER_UPPER: return "upper";
        default: return "unknown";
    }
}

Architecture3D::LayerType Architecture3D::stringToLayer(const std::string& layer_name)
{
    if (layer_name == "bottom") return LAYER_BOTTOM;
    if (layer_name == "upper") return LAYER_UPPER;
    return LAYER_BOTTOM;  // 默认
}

Architecture3D::LayerType Architecture3D::layerOfMaster(odb::dbMaster* master)
{
    if (!master) {
        return LAYER_BOTTOM;
    }
    // 标准单元后缀约定：上层单元名含 "upper"，下层含 "bottom"
    // （PDK 打包保证，同构/异构一致）。回退查 lib 名以增强鲁棒性。
    std::string name = master->getName();
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    if (name.find("upper") != std::string::npos) return LAYER_UPPER;
    if (name.find("bottom") != std::string::npos) return LAYER_BOTTOM;

    if (master->getLib()) {
        std::string lib = master->getLib()->getName();
        std::transform(lib.begin(), lib.end(), lib.begin(), ::tolower);
        if (lib.find("upper") != std::string::npos) return LAYER_UPPER;
        if (lib.find("bottom") != std::string::npos) return LAYER_BOTTOM;
    }
    return LAYER_BOTTOM;  // 默认归底层
}

Architecture3D::LayerType Architecture3D::layerOfBTerm(odb::dbBTerm* bterm,
                                                       odb::dbTech* tech)
{
    odb::dbTechLayer* hb = tech ? tech->findLayer("hb_layer") : nullptr;
    // 取该 IO 所有 pin shape 中最高的 ROUTING 金属层（绝对层序号）。
    int best = -1;
    if (bterm) {
        for (odb::dbBPin* bpin : bterm->getBPins()) {
            for (odb::dbBox* box : bpin->getBoxes()) {
                odb::dbTechLayer* l = box->getTechLayer();
                if (l && l->getType() == odb::dbTechLayerType::ROUTING) {
                    best = std::max(best, l->getNumber());
                }
            }
        }
    }
    // 兜底：无 hb_layer 或该 IO 无 routing 几何时，保持原行为（上层）。
    if (!hb || best < 0) {
        return LAYER_UPPER;
    }
    // 以 hb_layer 为界：pin 金属层在其下 -> 下层，在其上 -> 上层。
    return (best < hb->getNumber()) ? LAYER_BOTTOM : LAYER_UPPER;
}

bool Architecture3D::addLayer(LayerType layer_type, const std::string& layer_name)
{
    if (hasLayer(layer_type)) {
        return false;
    }
    
    layers_[layer_type].name = layer_name;
    return true;
}

bool Architecture3D::removeLayer(LayerType layer_type)
{
    if (!hasLayer(layer_type)) {
        return false;
    }
    
    clearLayer(layer_type);
    layers_.erase(layer_type);
    return true;
}

bool Architecture3D::hasLayer(LayerType layer_type) const
{
    return layers_.find(layer_type) != layers_.end();
}

const std::string& Architecture3D::getLayerName(LayerType layer_type) const
{
    static const std::string empty_string;
    auto it = layers_.find(layer_type);
    if (it != layers_.end()) {
        return it->second.name;
    }
    return empty_string;
}

const std::vector<Architecture3D::Row3D*>& Architecture3D::getRows(LayerType layer) const
{
    static const std::vector<Row3D*> empty_vector;
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.rows;
    }
    return empty_vector;
}

int Architecture3D::getNumRows(LayerType layer) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.rows.size();
    }
    return 0;
}

Architecture3D::Row3D* Architecture3D::getRow(LayerType layer, int row_idx) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end() && row_idx >= 0 && row_idx < it->second.rows.size()) {
        return it->second.rows[row_idx];
    }
    return nullptr;
}

Architecture3D::Row3D* Architecture3D::createAndAddRow(LayerType layer)
{
    if (!hasLayer(layer)) {
        return nullptr;
    }
    
    Row3D* row = new Row3D(layer);
    layers_[layer].rows.push_back(row);
    return row;
}

DbuX Architecture3D::getMinX(LayerType layer) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.xmin;
    }
    return DbuX{std::numeric_limits<DBU>::max()};
}

DbuX Architecture3D::getMaxX(LayerType layer) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.xmax;
    }
    return DbuX{std::numeric_limits<DBU>::lowest()};
}

DbuY Architecture3D::getMinY(LayerType layer) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.ymin;
    }
    return DbuY{std::numeric_limits<DBU>::max()};
}

DbuY Architecture3D::getMaxY(LayerType layer) const
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        return it->second.ymax;
    }
    return DbuY{std::numeric_limits<DBU>::lowest()};
}

void Architecture3D::setMinX(LayerType layer, DbuX xmin)
{
    if (hasLayer(layer)) {
        layers_[layer].xmin = xmin;
    }
}

void Architecture3D::setMaxX(LayerType layer, DbuX xmax)
{
    if (hasLayer(layer)) {
        layers_[layer].xmax = xmax;
    }
}

void Architecture3D::setMinY(LayerType layer, DbuY ymin)
{
    if (hasLayer(layer)) {
        layers_[layer].ymin = ymin;
    }
}

void Architecture3D::setMaxY(LayerType layer, DbuY ymax)
{
    if (hasLayer(layer)) {
        layers_[layer].ymax = ymax;
    }
}

void Architecture3D::clear()
{
    for (auto& layer_pair : layers_) {
        clearLayer(layer_pair.first);
    }
}

void Architecture3D::clearLayer(LayerType layer)
{
    auto it = layers_.find(layer);
    if (it != layers_.end()) {
        // 删除所有Row
        for (auto* row : it->second.rows) {
            delete row;
        }
        it->second.rows.clear();
        
        // 删除所有Region
        for (auto* region : it->second.regions) {
            delete region;
        }
        it->second.regions.clear();
        
        // 重置边界
        it->second.xmin = DbuX::max();
        it->second.xmax = DbuX::lowest();
        it->second.ymin = DbuY::max();
        it->second.ymax = DbuY::lowest();
        it->second.site_width = 0;
    }
}

// ========================================================================
// Row3D实现
// ========================================================================

Architecture3D::Row3D::Row3D(LayerType layer)
    : layer_(layer)
{
}

DBU Architecture3D::Row3D::getZOffset() const
{
    // 简化：下层Z=0，上层Z=芯片高度
    return (layer_ == LAYER_BOTTOM) ? 0 : 1000;  // 示例值
}

Pointdp3D Architecture3D::Row3D::get3DLocation(DbuX x, DbuY y) const
{
    return Pointdp3D(x.v, y.v, getZOffset());
}

void Architecture3D::Row3D::printInfo() const
{
    // 打印Row信息
}

// ========================================================================
// 
// ========================================================================

void* Architecture3D::createAndAddRegion(LayerType layer)
{
    if (!hasLayer(layer)) {
        return nullptr;
    }
    
    Group3D* region = new Group3D();
    region->setId(layers_[layer].regions.size());
    layers_[layer].regions.push_back(region);
    return region;
}

void Architecture3D::setSiteWidth(LayerType layer, const DbuX& site_width)
{
    if (hasLayer(layer)) {
        layers_[layer].site_width = site_width.v;
    }
}

void Architecture3D::printAllStats() const
{
    std::cout << "=== Architecture3D Statistics ===" << std::endl;
    
    for (int i = 0; i < LAYER_COUNT; i++) {
        LayerType layer = static_cast<LayerType>(i);
        if (hasLayer(layer)) {
            const auto& layer_info = layers_.at(layer);
            std::cout << "Layer " << layerToString(layer) << " (" << layer_info.name << "):" << std::endl;
            std::cout << "  Rows: " << layer_info.rows.size() << std::endl;
            std::cout << "  Regions: " << layer_info.regions.size() << std::endl;
            std::cout << "  Site width: " << layer_info.site_width << " DBU" << std::endl;
            std::cout << "  Bounds: X[" << layer_info.xmin.v << ", " << layer_info.xmax.v 
                      << "], Y[" << layer_info.ymin.v << ", " << layer_info.ymax.v << "]" << std::endl;
        }
    }
    
    std::cout << "HBT pitch: " << hbt_pitch_ << " DBU" << std::endl;
    std::cout << "Inter-layer spacing: " << inter_layer_spacing_ << " DBU" << std::endl;
}

// ========================================================================
// 新增方法实现
// ========================================================================

DBU Architecture3D::getRowHeight() const
{
    // 遍历所有层，找到第一个有行的层，返回其行高
    for (const auto& layer_pair : layers_) {
        const auto& rows = layer_pair.second.rows;
        if (!rows.empty()) {
            // 返回第一个行的高度
            return rows[0]->getHeight().v;
        }
    }
    
    // 如果没有行，返回0表示错误
    return 0;
}

bool Architecture3D::isCellInRows(const Node3D* cell) const
{
    if (!cell) {
        return false;
    }
    
    // 获取单元所在的层
    LayerType layer = stringToLayer(cell->layer);
    if (!hasLayer(layer)) {
        return false;
    }
    
    // 获取单元边界框
    odb::Rect cell_bbox = cell->getBBox();
    if (cell_bbox.isInverted()) {
        return false;
    }
    
    // 检查单元是否在任何一行内
    const auto& rows = getRows(layer);
    for (const auto* row : rows) {
        // 检查单元是否与行重叠
        DbuY row_bottom = row->getBottom();
        DbuY row_top = row->getTop();
        DbuX row_left = row->getLeft();
        DbuX row_right = row->getRight();
        
        // 检查Y方向重叠
        bool y_overlap = (cell_bbox.yMin() < row_top.v) && (cell_bbox.yMax() > row_bottom.v);
        
        // 检查X方向重叠
        bool x_overlap = (cell_bbox.xMin() < row_right.v) && (cell_bbox.xMax() > row_left.v);
        
        if (y_overlap && x_overlap) {
            return true;
        }
    }
    
    return false;
}

Architecture3D::LayerType Architecture3D::getCellLayer(odb::dbInst* inst) const
{
    if (!inst) {
        return LAYER_BOTTOM;
    }
    // 统一依据标准单元名后缀判定层归属。
    return layerOfMaster(inst->getMaster());
}

} // namespace tdp
