// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "Network3D.h"
#include "Architecture3D.h"
#include "Master3D.h"
#include "Pin3D.h"
#include "Common3D.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <limits>

namespace tdp {

Network3D::Network3D() {
    // 构造函数实现
}

Network3D::~Network3D() {
    // 析构函数实现
    clear();
}

void Network3D::addCell(Node3D* cell) {
    cells_.push_back(cell);
    name_to_cell_map_[cell->name] = cell;
    if (cell->db_inst) {
        inst_to_cell_map_[cell->db_inst] = cell;
    }
}

void Network3D::removeCell(Node3D* cell) {
    // 从所有向量中移除
    auto it = std::find(cells_.begin(), cells_.end(), cell);
    if (it != cells_.end()) {
        cells_.erase(it);
    }
    
    // 从映射中移除
    name_to_cell_map_.erase(cell->name);
    if (cell->db_inst) {
        inst_to_cell_map_.erase(cell->db_inst);
    }
    
    // 更新分类
    updateCellCategories();
    
    // 删除对象
    delete cell;
}

Node3D* Network3D::createCellFromDatabase(odb::dbInst* inst, Architecture3D* arch) {
    Node3D* node = new Node3D();
    node->db_inst = inst;
    node->name = inst->getName();
    
    // 设置位置
    int x, y;
    inst->getLocation(x, y);
    node->location = Pointdp3D(DbuX{x}, DbuY{y}, DbuZ{0});
    
    // 设置固定状态
    node->is_fixed = inst->isFixed();
    node->is_movable = !inst->isFixed();
    
    // 设置层信息
    if (arch) {
        // 使用Network3D自己的getCellLayer方法
        Architecture3D::LayerType layer_type = getCellLayer(inst);
        node->layer = (layer_type == Architecture3D::LAYER_BOTTOM) ? "bottom" : "upper";
    }
    // 设置ID 信息
    if (inst) {
        node->setId(inst->getId());  // 使用 dbInst 的 ID
    } else {
        // 如果没有 dbInst（如 Terminal），使用其他方式生成 ID
        static int next_id = 1000000;
        node->setId(next_id++);
    }
    return node;
}

void Network3D::addNode(odb::dbInst* inst) {
    Node3D* node = createCellFromDatabase(inst, nullptr);
    addCell(node);
}

void Network3D::addNode(odb::dbBTerm* term) {
    Node3D* node = new Node3D();
    node->setBTerm(term);
    node->name = term->getName();
    node->setType(Node3D::TERMINAL);  // 设置为Terminal类型
    node->is_fixed = true;  // 终端通常是固定的
    node->is_movable = false;
    // 按 IO pin 几何所在金属层相对 hb_layer 判定其所属层（取代写死 "upper"）。
    auto bt_layer = Architecture3D::layerOfBTerm(term, term->getBlock()->getTech());
    node->layer = (bt_layer == Architecture3D::LAYER_BOTTOM) ? "bottom" : "upper";
    
    // 设置位置和尺寸（从term获取）
    odb::Rect bbox = term->getBBox();
    if (!bbox.isInverted()) {
        node->location = Pointdp3D(DbuX{bbox.xMin()}, DbuY{bbox.yMin()}, DbuZ{0});
        node->setWidth(DbuX{bbox.xMax() - bbox.xMin()});
        node->setHeight(DbuY{bbox.yMax() - bbox.yMin()});
    }
    
    addCell(node);
    term_to_node_map_[term] = node;
    terminals_cnt_++;
}

void Network3D::addFillerNode(DbuX left, DbuY bottom, DbuX width, DbuY height) {
    Node3D* node = new Node3D();
    node->name = "filler_" + std::to_string(cells_.size());
    node->is_fixed = true;
    node->is_movable = false;
    node->location = Pointdp3D(left, bottom, DbuZ{0});
    
    addCell(node);
}

Node3D* Network3D::getCell(odb::dbInst* inst) const {
    auto it = inst_to_cell_map_.find(inst);
    if (it != inst_to_cell_map_.end()) {
        return it->second;
    }
    return nullptr;
}

Node3D* Network3D::getCell(const std::string& name) const {
    auto it = name_to_cell_map_.find(name);
    if (it != name_to_cell_map_.end()) {
        return it->second;
    }
    return nullptr;
}

Node3D* Network3D::getNode(odb::dbBTerm* term) const {
    // 首先使用映射查找
    auto it = term_to_node_map_.find(term);
    if (it != term_to_node_map_.end()) {
        return it->second;
    }
    // 回退：通过名称查找
    std::string term_name = term->getName();
    return getCell(term_name);
}

void Network3D::addNet(Net3D* net) {
    nets_.push_back(net);
    name_to_net_map_[net->name] = net;
    if (net->db_net) {
        net_to_net3d_map_[net->db_net] = net;
    }
}

void Network3D::removeNet(Net3D* net) {
    // 从所有向量中移除
    auto it = std::find(nets_.begin(), nets_.end(), net);
    if (it != nets_.end()) {
        nets_.erase(it);
    }
    
    // 从映射中移除
    name_to_net_map_.erase(net->name);
    if (net->db_net) {
        net_to_net3d_map_.erase(net->db_net);
    }
    
    // 删除对象
    delete net;
}

Net3D* Network3D::createNetFromDatabase(odb::dbNet* net) {
    Net3D* net3d = new Net3D();
    net3d->db_net = net;
    net3d->name = net->getName();
    
    // 设置时序信息
    // 这里可以添加从数据库获取时序信息的逻辑
    
    return net3d;
}

void Network3D::addEdge(odb::dbNet* net) {
    Net3D* net3d = createNetFromDatabase(net);
    addNet(net3d);
}

Master3D* Network3D::addMaster(odb::dbMaster* db_master, Architecture3D* arch) {
    // 检查是否已存在
    auto it = master_to_idx_.find(db_master);
    if (it != master_to_idx_.end()) {
        return masters_[it->second].get();
    }
    
    // 创建新的Master3D
    std::unique_ptr<Master3D> master = std::make_unique<Master3D>();
    master->setDbMaster(db_master);
    
    // 设置属性 - 使用Master3D自己的判断逻辑
    // master->setMultiRow()已经在setDbMaster中设置
    
    // 添加到列表
    int idx = masters_.size();
    masters_.push_back(std::move(master));
    master_to_idx_[db_master] = idx;
    
    return masters_[idx].get();
}

Net3D* Network3D::getNet(odb::dbNet* net) const {
    auto it = net_to_net3d_map_.find(net);
    if (it != net_to_net3d_map_.end()) {
        return it->second;
    }
    return nullptr;
}

Net3D* Network3D::getNet(const std::string& name) const {
    auto it = name_to_net_map_.find(name);
    if (it != name_to_net_map_.end()) {
        return it->second;
    }
    return nullptr;
}

Net3D* Network3D::getEdge(odb::dbNet* net) const {
    // 对于3D，Edge和Net是同一个概念
    return getNet(net);
}

void Network3D::clear() {
    // 删除所有动态分配的对象
    for (auto* cell : cells_) {
        delete cell;
    }
    cells_.clear();
    
    for (auto* net : nets_) {
        delete net;
    }
    nets_.clear();
    
    // 清空智能指针容器
    masters_.clear();
    pins_.clear();
    
    // 清空其他容器
    movable_cells_.clear();
    fixed_cells_.clear();
    blockages_.clear();
    
    // 清空映射
    inst_to_cell_map_.clear();
    net_to_net3d_map_.clear();
    name_to_cell_map_.clear();
    name_to_net_map_.clear();
    master_to_idx_.clear();
    net_to_edge_idx_.clear();
    term_to_node_map_.clear();
    cells_cnt_ = 0;
    terminals_cnt_ = 0;
}

void Network3D::addHBTConnection(const std::string& bottom_cell_name, 
                                const std::string& upper_cell_name,
                                const Point& location) {
    Node3D* bottom_cell = getCell(bottom_cell_name);
    Node3D* upper_cell = getCell(upper_cell_name);
    
    if (bottom_cell && upper_cell) {
        // 创建HBT连接
        Net3D::HBTConnection hbt;
        hbt.location = location;
        hbt.bottom_cell = bottom_cell;
        hbt.upper_cell = upper_cell;
        hbt.pitch = 0;  // 需要从配置获取
        hbt.aligned = false;
        
        // 添加到所有相关网络
        for (auto* net : nets_) {
            // 检查网络是否连接这两个单元
            bool has_bottom = false;
            bool has_upper = false;
            
            for (auto* cell : net->connected_cells) {
                if (cell == bottom_cell) has_bottom = true;
                if (cell == upper_cell) has_upper = true;
            }
            
            if (has_bottom && has_upper) {
                net->addHBTConnection(hbt);
            }
        }
    }
}

std::vector<Net3D::HBTConnection> Network3D::getHBTConnections() const {
    std::vector<Net3D::HBTConnection> all_hbt_connections;
    
    for (const auto* net : nets_) {
        const auto& net_hbts = net->hbt_connections;
        all_hbt_connections.insert(all_hbt_connections.end(), net_hbts.begin(), net_hbts.end());
    }
    
    return all_hbt_connections;
}

double Network3D::calculateTotalWirelength() const {
    double total_wl = 0.0;
    
    for (const auto* net : nets_) {
        // 使用2D HPWL计算（向后兼容）
        total_wl += net->hpwl();
    }
    
    return total_wl;
}

double Network3D::calculateTotalTimingCost() const {
    double total_cost = 0.0;
    
    for (const auto* net : nets_) {
        // 基于时序关键性计算成本
        total_cost += net->timing_criticality;
    }
    
    return total_cost;
}

double Network3D::calculateTotalAlignmentCost() const {
    double total_cost = 0.0;
    
    for (const auto* net : nets_) {
        total_cost += net->calculateAlignmentCost();
    }
    
    return total_cost;
}

bool Network3D::checkPlacementLegality(const Architecture3D* arch) const {
    if (!arch) {
        return false;
    }
    
    // 检查所有单元
    for (const auto* cell : cells_) {
        if (!cell->is_fixed) {
            // 检查位置是否合法
            Pointdp3D legal_loc = cell->getLegalLocation(arch);
            if (cell->location.x.v != legal_loc.x.v || 
                cell->location.y.v != legal_loc.y.v ||
                cell->location.z.v != legal_loc.z.v) {
                return false;  // 位置不合法
            }
            
            // 检查是否在行内
            if (!arch->isCellInRows(cell)) {
                return false;
            }
        }
    }
    
    // 检查重叠
    for (size_t i = 0; i < cells_.size(); i++) {
        for (size_t j = i + 1; j < cells_.size(); j++) {
            if (checkOverlap(cells_[i], cells_[j])) {
                return false;  // 有重叠
            }
        }
    }
    
    // 检查HBT对齐
    for (const auto* net : nets_) {
        if (!net->checkHBTAlignment()) {
            return false;  // HBT未对齐
        }
    }
    
    return true;
}

Architecture3D::LayerType Network3D::getCellLayer(odb::dbInst* inst) const {
    Node3D* cell = getCell(inst);
    if (cell) {
        return (cell->layer == "bottom") ? Architecture3D::LAYER_BOTTOM : Architecture3D::LAYER_UPPER;
    }
    
    // 默认返回底层
    return Architecture3D::LAYER_BOTTOM;
}

void Network3D::printStats() const {
    std::cout << "Network3D Stats:" << std::endl;
    std::cout << "  Total cells: " << cells_.size() << std::endl;
    std::cout << "  Movable cells: " << movable_cells_.size() << std::endl;
    std::cout << "  Fixed cells: " << fixed_cells_.size() << std::endl;
    std::cout << "  Nets: " << nets_.size() << std::endl;
    std::cout << "  Masters: " << masters_.size() << std::endl;
    std::cout << "  Pins: " << pins_.size() << std::endl;
    std::cout << "  Blockages: " << blockages_.size() << std::endl;
    
    // 计算统计信息
    double total_wl = calculateTotalWirelength();
    double total_timing = calculateTotalTimingCost();
    double total_alignment = calculateTotalAlignmentCost();
    
    std::cout << "  Total wirelength: " << total_wl << " DBU" << std::endl;
    std::cout << "  Total timing cost: " << total_timing << std::endl;
    std::cout << "  Total alignment cost: " << total_alignment << std::endl;
    
    // 层统计
    int bottom_cells = 0;
    int upper_cells = 0;
    for (const auto* cell : cells_) {
        if (cell->layer == "bottom") bottom_cells++;
        else if (cell->layer == "upper") upper_cells++;
    }
    std::cout << "  Bottom layer cells: " << bottom_cells << std::endl;
    std::cout << "  Upper layer cells: " << upper_cells << std::endl;
    
    // HBT统计
    int total_hbt = 0;
    int aligned_hbt = 0;
    for (const auto* net : nets_) {
        total_hbt += net->hbt_connections.size();
        for (const auto& hbt : net->hbt_connections) {
            if (hbt.aligned) aligned_hbt++;
        }
    }
    std::cout << "  HBT connections: " << total_hbt << std::endl;
    std::cout << "  Aligned HBTs: " << aligned_hbt << " (" 
              << (total_hbt > 0 ? (100.0 * aligned_hbt / total_hbt) : 0.0) << "%)" << std::endl;
}

void Network3D::updateCellCategories() {
    movable_cells_.clear();
    fixed_cells_.clear();
    
    for (auto* cell : cells_) {
        if (cell->is_fixed) {
            fixed_cells_.push_back(cell);
        } else {
            movable_cells_.push_back(cell);
        }
    }
}

// ============================================================================
// 新增方法实现
// ============================================================================

Pin3D* Network3D::addPin(odb::dbITerm* iterm) {
    Pin3D* pin = new Pin3D();
    pin->setDbITerm(iterm);
    pin->setName(iterm->getName());
    
    // 设置连接
    odb::dbInst* inst = iterm->getInst();
    Node3D* node = getCell(inst);
    if (node) {
        pin->setNode(node);
    }
    
    // 添加到列表
    pins_.push_back(std::unique_ptr<Pin3D>(pin));
    return pin;
}

Pin3D* Network3D::addPin(odb::dbBTerm* bterm) {
    Pin3D* pin = new Pin3D();
    pin->setDbBTerm(bterm);
    pin->setName(bterm->getName());
    
    // 设置连接
    Node3D* node = getNode(bterm);
    if (node) {
        pin->setNode(node);
    }
    
    // 添加到列表
    pins_.push_back(std::unique_ptr<Pin3D>(pin));
    return pin;
}

void Network3D::connect(Pin3D* pin, Node3D* node) {
    pin->setNode(node);
    node->addPin(pin);
}

void Network3D::connect(Pin3D* pin, Net3D* net) {
    pin->setNet(net);
    net->addPin(pin);
}

void Network3D::createAndAddBlockage(const odb::Rect& bounds) {
    blockages_.push_back(bounds);
}

Master3D* Network3D::getMaster(odb::dbMaster* db_master) {
    auto it = master_to_idx_.find(db_master);
    if (it != master_to_idx_.end()) {
        return masters_[it->second].get();
    }
    return nullptr;
}

void Network3D::setCore(const odb::Rect& core) {
    core_ = core;
}

const odb::Rect& Network3D::getCore() const {
    return core_;
}

bool Network3D::checkOverlap(const Node3D* cell1, const Node3D* cell2) const {
    if (!cell1 || !cell2) {
        return false;
    }
    
    // 获取边界框
    odb::Rect bbox1 = cell1->getBBox();
    odb::Rect bbox2 = cell2->getBBox();
    
    // 检查重叠
    return bbox1.xMin() < bbox2.xMax() && 
           bbox1.xMax() > bbox2.xMin() && 
           bbox1.yMin() < bbox2.yMax() && 
           bbox1.yMax() > bbox2.yMin();
}

std::vector<odb::Rect> Network3D::getBlockages() const {
    return blockages_;
}

int Network3D::getNumBlockages() const {
    return blockages_.size();
}

void Network3D::updateAllHBTAlignment() {
    for (auto* net : nets_) {
        net->updateHBTAlignment();
    }
}

double Network3D::calculateTotalCost(double wireweight, double timingweight, double alignmentweight) const {
    double wire_cost = calculateTotalWirelength();
    double timing_cost = calculateTotalTimingCost();
    double alignment_cost = calculateTotalAlignmentCost();
    
    return wireweight * wire_cost + 
           timingweight * timing_cost + 
           alignmentweight * alignment_cost;
}

void Network3D::updateTimingInfo() {
    // 更新所有网络的时序信息
    // 这里可以添加实际的时序分析逻辑
    for (auto* net : nets_) {
        // 简化实现：设置默认值
        net->updateTimingInfo(0.0);
    }
}

} // namespace tdp
