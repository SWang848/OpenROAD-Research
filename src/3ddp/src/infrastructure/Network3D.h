// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include "Common3D.h"
#include "Architecture3D.h"
#include "Node3D.h"
#include "Net3D.h"
#include "Master3D.h"
#include "Pin3D.h"

// Forward declarations
namespace odb {
    class dbInst;
    class dbNet;
    class dbITerm;
    class dbBTerm;
    class dbMaster;
}

namespace tdp {

// Forward declarations
class Edge3D;

/**
 * @class Network3D
 * @brief 3D network representation similar to DPL's Network class
 * 
 * Manages all cells and nets in the 3D design.
 */
class Network3D
{
public:
    Network3D();
    ~Network3D();
    
    // 禁止拷贝
    Network3D(const Network3D&) = delete;
    Network3D& operator=(const Network3D&) = delete;
    
    // Cell管理（类似DPL的Network）
    void addCell(Node3D* cell);
    void removeCell(Node3D* cell);
    
    // 从数据库创建Cell
    Node3D* createCellFromDatabase(odb::dbInst* inst, Architecture3D* arch);
    
    // DPL风格的方法
    void addNode(odb::dbInst* inst);
    void addNode(odb::dbBTerm* term);
    void addFillerNode(DbuX left, DbuY bottom, DbuX width, DbuY height);
    
    Node3D* getCell(odb::dbInst* inst) const;
    Node3D* getCell(const std::string& name) const;
    Node3D* getNode(odb::dbBTerm* term) const;
    
    std::vector<Node3D*>& getCells() { return cells_; }
    const std::vector<Node3D*>& getCells() const { return cells_; }
    
    std::vector<Node3D*>& getMovableCells() { return movable_cells_; }
    std::vector<Node3D*>& getFixedCells() { return fixed_cells_; }

    // 重新计算 movable/fixed 分类（在批量添加 cell 后调用）
    void updateCellCategories();
    
    // Net管理（类似DPL的Network）
    void addNet(Net3D* net);
    void removeNet(Net3D* net);
    
    // 从数据库创建Net（类似DPL的createNetwork）
    Net3D* createNetFromDatabase(odb::dbNet* net);
    
    // DPL风格的方法
    void addEdge(odb::dbNet* net);
    
    // Master管理（类似DPL的Network）
    Master3D* addMaster(odb::dbMaster* db_master, Architecture3D* arch);
    Master3D* getMaster(odb::dbMaster* db_master);
    
    Net3D* getNet(odb::dbNet* net) const;
    Net3D* getNet(const std::string& name) const;
    Net3D* getEdge(odb::dbNet* net) const;
    
    // 清除所有数据（类似DPL的clear）
    void clear();
    
    // 层间连接
    void addHBTConnection(const std::string& bottom_cell_name, 
                         const std::string& upper_cell_name,
                         const Point& location);
    
    std::vector<Net3D::HBTConnection> getHBTConnections() const;
    
    // 统计信息
    int getNumCells() const { return cells_.size(); }
    int getNumNets() const { return nets_.size(); }
    int getNumMovableCells() const { return movable_cells_.size(); }
    int getNumFixedCells() const { return fixed_cells_.size(); }
    uint32_t getNumTerminals() const { return terminals_cnt_; }
    
    // 成本计算
    double calculateTotalWirelength() const;
    double calculateTotalTimingCost() const;
    double calculateTotalAlignmentCost() const;
    double calculateTotalCost(double wireweight, double timingweight, double alignmentweight) const;
    
    // 验证
    bool checkPlacementLegality(const Architecture3D* arch) const;
    
    // 层信息获取
    Architecture3D::LayerType getCellLayer(odb::dbInst* inst) const;
    
    // 调试
    void printStats() const;
    
    // 引脚管理（类似DPL的Network）
    Pin3D* addPin(odb::dbITerm* iterm);
    Pin3D* addPin(odb::dbBTerm* bterm);
    void connect(Pin3D* pin, Node3D* node);
    void connect(Pin3D* pin, Net3D* net);
    
    // 障碍物管理
    void createAndAddBlockage(const odb::Rect& bounds);
    std::vector<odb::Rect> getBlockages() const;
    int getNumBlockages() const;
    
    // 核心区域管理
    void setCore(const odb::Rect& core);
    const odb::Rect& getCore() const;
    
    // 重叠检查
    bool checkOverlap(const Node3D* cell1, const Node3D* cell2) const;
    
    // HBT对齐更新
    void updateAllHBTAlignment();
    
    // 时序信息更新
    void updateTimingInfo();
    
private:
    // Cells和Nets
    std::vector<Node3D*> cells_;
    std::vector<Net3D*> nets_;
    
    // 优化访问的数据结构
    std::vector<Node3D*> movable_cells_;
    std::vector<Node3D*> fixed_cells_;
    
    // 映射表
    std::unordered_map<odb::dbInst*, Node3D*> inst_to_cell_map_;
    std::unordered_map<odb::dbNet*, Net3D*> net_to_net3d_map_;
    std::unordered_map<std::string, Node3D*> name_to_cell_map_;
    std::unordered_map<std::string, Net3D*> name_to_net_map_;
    std::unordered_map<odb::dbBTerm*, Node3D*> term_to_node_map_;  // BTerm到Node映射
    
    // Master管理
    std::vector<std::unique_ptr<Master3D>> masters_;
    std::unordered_map<odb::dbMaster*, int> master_to_idx_;
    
    // 引脚管理
    std::vector<std::unique_ptr<Pin3D>> pins_;
    
    // 障碍物
    std::vector<odb::Rect> blockages_;
    
    // 核心区域
    odb::Rect core_;

    // 计数器（参考DPL）
    uint32_t cells_cnt_{0};
    uint32_t terminals_cnt_{0};
    
    // 网络到边索引
    std::unordered_map<odb::dbNet*, int> net_to_edge_idx_;
    
};

} // namespace tdp
