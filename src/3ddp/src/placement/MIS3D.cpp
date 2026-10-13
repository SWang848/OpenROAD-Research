// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "MIS3D.h"
#include "placement.h"
#include "tdp/Timingdp.h"
#include "color.h"
#include "../eva/eva.h"
#include "../eva/TimingEngine3D.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Pin3D.h"
#include "../infrastructure/Common3D.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <random>

#include <lemon/cost_scaling.h>
#include <lemon/cycle_canceling.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/preflow.h>
#include <lemon/smart_graph.h>

#include "utl/Logger.h"


namespace tdp {

using utl::TDP;
////////////////////////////////////////////////////////////////////////////////
// DetailedMis3D 实现
////////////////////////////////////////////////////////////////////////////////

DetailedMis3D::DetailedMis3D(Timingdp* timingdp, LSP3D* lsp3d)
    : arch_(timingdp ? timingdp->getArchitecture3D() : nullptr),
      network_(timingdp ? timingdp->getNetwork3D() : nullptr),
      logger_(timingdp ? timingdp->getLogger() : nullptr),
      timingdp_(timingdp),
      lsp3d_(lsp3d),
      sta_(timingdp ? timingdp->getSta() : nullptr),
      timing_engine_(timingdp ? timingdp->gettimingengine() : nullptr),
      segment_manager_(timingdp ? timingdp->getSegmentManager() : nullptr)
      
{
    // 构造函数体可以为空，或者做额外检查
    if (timingdp) {
        params_ = DetailedMis3DParams();
    }
}

DetailedMis3D::~DetailedMis3D() {
    clearGrid3D();
}

////////////////////////////////////////////////////////////////////////////////
// 主执行函数
////////////////////////////////////////////////////////////////////////////////

void DetailedMis3D::run() {
    
    if (lsp3d_) {
        lsp3d_->initWeight(timingdp_, timing_engine_);
        lsp3d_->updatePinWeights();
    }
    // TODO: 解析参数并设置params_
    // TODO: 记录开始时间
    
    // 主算法流程
    collectMovableCells();
    if (candidates_.empty()) {
        logger_->info(TDP, 1403, "No movable cells found");
        return;
    }
    
    colorCells3D();
    buildGrid3D();
    
    for (int pass = 0; pass < params_.maxPasses; ++pass) {
        sta_->updateTiming(true);  
        place3D();
        lsp3d_->updatePinWeights(); 
        // TODO: 检查收敛条件
    }
    
    // updateStatistics();
}

////////////////////////////////////////////////////////////////////////////////
// 核心算法函数 
////////////////////////////////////////////////////////////////////////////////

void DetailedMis3D::collectMovableCells() {
    candidates_.clear();
    
    // 从 Network3D 获取所有可移动单元
    std::vector<Node3D*>& movable_cells = network_->getMovableCells();
    candidates_.insert(candidates_.end(),
                       movable_cells.begin(),
                       movable_cells.end());
}


void DetailedMis3D::colorCells3D() {
    // 1. 初始化颜色数组（使用 Node3D）
    colors_.resize(network_->getNumCells());
    std::ranges::fill(colors_, -1);
    
    // 2. 标记可移动单元
    movable_.resize(network_->getNumCells());
    movable_.assign(movable_.size(), false);
    for (const Node3D* ndi : candidates_) {
        movable_[ndi->getId()] = true;
    }
    
    // 3. 构建冲突图（需要考虑3D层间连接）
    ColorGraph gr(network_->getNumCells());
    
    // 遍历所有网络（Net3D)
    auto* network_3d = timingdp_->getNetwork3D();
    const auto& nets = timingdp_->getNets();
    for (auto* db_net : nets) {
        auto* net = network_3d->getNet(db_net);
        if(!net){
        std::cout << "  db_net pointer: " << db_net 
            << ", name: " << db_net->getName() 
            << ", signal type: " << (db_net->getSigType().getString()) << std::endl;
            continue;
        }
        const int numPins = net->getNumPins();
        
        // 跳过小型和大型网络
        if (numPins <= 1 || numPins > params_.skipNetsLargerThanThis) {
            continue;
        }
        
        // 为共享网络的单元添加边
        for (int pi = 0; pi < numPins; pi++) {
            const Pin3D* pini = net->getPins()[pi];
            const Node3D* ndi = pini->getNode();
            if (!movable_[ndi->getId()]) continue;
            
            for (int pj = pi + 1; pj < numPins; pj++) {
                const Pin3D* pinj = net->getPins()[pj];
                const Node3D* ndj = pinj->getNode();
                if (!movable_[ndj->getId()]) continue;
                if (ndj == ndi) continue;
                
                
                gr.addEdge(ndi->getId(), ndj->getId());
            }
        }
    }
    
    // 4. 执行贪心着色
    gr.greedyColoring();
    
    // 5. 分配颜色
    const auto& cells = network_->getCells();
    for (int i = 0; i < cells.size(); i++) {
        const Node3D* ndi = cells[i];
        const int color = gr.getColor(i);
        if (movable_[ndi->getId()]) {
            colors_[ndi->getId()] = color;
        }
    }
}

void DetailedMis3D::buildGrid3D()
{
    // Builds separate coarse grids for bottom and upper layers.
    traversal_ = 0;

    // 清理旧网格
    clearGrid3D();

    // 定义层类型变量
    Architecture3D::LayerType bottom_layer = Architecture3D::stringToLayer("bottom");
    Architecture3D::LayerType upper_layer = Architecture3D::stringToLayer("upper");
    
    // 为每层创建独立的网格
    grid_bottom_.clear();
    grid_upper_.clear();
    
    // ==================== 底层网格 ====================
    const double xmin_bottom = arch_->getMinX(bottom_layer).v;
    const double xmax_bottom = arch_->getMaxX(bottom_layer).v;
    const double ymin_bottom = arch_->getMinY(bottom_layer).v;
    const double ymax_bottom = arch_->getMaxY(bottom_layer).v;
    
    
    // 计算底层单元平均尺寸
    double avgH_bottom = 0., avgW_bottom = 0.;
    int bottom_count = 0;
    for (const Node3D* ndi : candidates_) {
        if (ndi->getLayer() == "bottom") {
            avgH_bottom += ndi->getHeight().v;
            avgW_bottom += ndi->getWidth().v;
            bottom_count++;
        }
    }
    
    if (bottom_count > 0) {
        avgH_bottom /= bottom_count;
        avgW_bottom /= bottom_count;
        
        // 设计底层网格步长
        stepX_bottom_ = avgW_bottom * std::sqrt(maxProblemSize_);
        stepY_bottom_ = avgH_bottom * std::sqrt(maxProblemSize_);
        
        // 计算底层网格维度
        dimX_bottom_ = (int) std::ceil((xmax_bottom - xmin_bottom) / stepX_bottom_);
        dimY_bottom_ = (int) std::ceil((ymax_bottom - ymin_bottom) / stepY_bottom_);
        
        // 创建底层网格
        grid_bottom_.resize(dimX_bottom_);
        for (int i = 0; i < dimX_bottom_; i++) {
            grid_bottom_[i].resize(dimY_bottom_);
            for (int j = 0; j < dimY_bottom_; j++) {
                auto bucket = new Bucket3D;
                bucket->xmin = xmin_bottom + i * stepX_bottom_;
                bucket->xmax = xmin_bottom + (i + 1) * stepX_bottom_;
                bucket->ymin = ymin_bottom + j * stepY_bottom_;
                bucket->ymax = ymin_bottom + (j + 1) * stepY_bottom_;
                bucket->i = i;
                bucket->j = j;
                bucket->k = 0;           // 底层索引
                bucket->travId = traversal_;
                bucket->layer = "bottom";
                grid_bottom_[i][j] = bucket;
            }
        }
    }
    
    // ==================== 顶层网格 ====================
    const double xmin_upper = arch_->getMinX(upper_layer).v;
    const double xmax_upper = arch_->getMaxX(upper_layer).v;
    const double ymin_upper = arch_->getMinY(upper_layer).v;
    const double ymax_upper = arch_->getMaxY(upper_layer).v;
    
    // 计算顶层单元平均尺寸
    double avgH_upper = 0., avgW_upper = 0.;
    int upper_count = 0;
    for (const Node3D* ndi : candidates_) {
        if (ndi->getLayer() == "upper") {
            avgH_upper += ndi->getHeight().v;
            avgW_upper += ndi->getWidth().v;
            upper_count++;
        }
    }
    
    if (upper_count > 0) {
        avgH_upper /= upper_count;
        avgW_upper /= upper_count;
        
        // 设计顶层网格步长
        stepX_upper_ = avgW_upper * std::sqrt(maxProblemSize_);
        stepY_upper_ = avgH_upper * std::sqrt(maxProblemSize_);
        
        // 计算顶层网格维度
        dimX_upper_ = (int) std::ceil((xmax_upper - xmin_upper) / stepX_upper_);
        dimY_upper_ = (int) std::ceil((ymax_upper - ymin_upper) / stepY_upper_);
        
        // 创建顶层网格
        grid_upper_.resize(dimX_upper_);
        for (int i = 0; i < dimX_upper_; i++) {
            grid_upper_[i].resize(dimY_upper_);
            for (int j = 0; j < dimY_upper_; j++) {
                auto bucket = new Bucket3D;
                bucket->xmin = xmin_upper + i * stepX_upper_;
                bucket->xmax = xmin_upper + (i + 1) * stepX_upper_;
                bucket->ymin = ymin_upper + j * stepY_upper_;
                bucket->ymax = ymin_upper + (j + 1) * stepY_upper_;
                bucket->i = i;
                bucket->j = j;
                bucket->k = 1;           // 顶层索引
                bucket->travId = traversal_;
                bucket->layer = "upper";
                grid_upper_[i][j] = bucket;
            }
        }
    }
}

void DetailedMis3D::clearGrid3D()
{   
    // 清理底层网格
    for (auto& row : grid_bottom_) {
        for (auto bucket : row) {
            delete bucket;
        }
        row.clear();
    }
    grid_bottom_.clear();
    
    // 清理顶层网格
    for (auto& row : grid_upper_) {
        for (auto bucket : row) {
            delete bucket;
        }
        row.clear();
    }
    grid_upper_.clear();
}

void DetailedMis3D::populateGrid3D()
{
    
    // 1. 清空两层网格
    for (auto& row : grid_bottom_) {
        for (auto& bucket : row) {
            if (bucket) bucket->clear();
        }
    }
    
    for (auto& row : grid_upper_) {
        for (auto& bucket : row) {
            if (bucket) bucket->clear();
        }
    }
    
    // 2. 清空映射表
    cellToBinMap_.clear();
    
    // 3. 获取两层布局区域原点
    Architecture3D::LayerType bottom_layer = Architecture3D::stringToLayer("bottom");
    Architecture3D::LayerType upper_layer = Architecture3D::stringToLayer("upper");
    
    const double xmin_bottom = arch_->getMinX(bottom_layer).v;
    const double ymin_bottom = arch_->getMinY(bottom_layer).v;
    const double xmin_upper = arch_->getMinX(upper_layer).v;
    const double ymin_upper = arch_->getMinY(upper_layer).v;
    
    // 4. 插入单元到对应的层网格
    for (Node3D* ndi : candidates_) {
        // 获取单元所在层
        std::string cell_layer = ndi->getLayer();  // 假设返回 "bottom" 或 "upper"
        
        if (cell_layer == "bottom") {
            // 插入到底层网格
            populateCellToLayerGrid(ndi, grid_bottom_, xmin_bottom, ymin_bottom,
                                   stepX_bottom_, stepY_bottom_,
                                   dimX_bottom_, dimY_bottom_);
        } else if (cell_layer == "upper") {
            // 插入到顶层网格
            populateCellToLayerGrid(ndi, grid_upper_, xmin_upper, ymin_upper,
                                   stepX_upper_, stepY_upper_,
                                   dimX_upper_, dimY_upper_);
        } else {
            std::cerr << "Error: No find cell layer " << std::endl;
            continue;
        }
    }
}

// 辅助函数：将单个单元插入到指定层网格
void DetailedMis3D::populateCellToLayerGrid(
    Node3D* ndi,
    std::vector<std::vector<Bucket3D*>>& grid,
    double xmin, double ymin,
    double stepX, double stepY,
    int dimX, int dimY)
{
    if (!ndi || dimX <= 0 || dimY <= 0 || stepX <= 0 || stepY <= 0) {
        return;
    }
    
    // 计算单元中心坐标
    const double y = ndi->getLocation().y.v + 0.5 * ndi->getHeight().v;
    const double x = ndi->getLocation().x.v + 0.5 * ndi->getWidth().v;
    
    // 计算网格索引（限制在有效范围内）
    const int j = std::max(std::min((int)((y - ymin) / stepY), dimY - 1), 0);
    const int i = std::max(std::min((int)((x - xmin) / stepX), dimX - 1), 0);
    
    // 获取对应的Bucket
    Bucket3D* bucket = grid[i][j];
    if (!bucket) {
        return;
    }
    
    // 将单元插入到Bucket中
    bucket->nodes.push_back(ndi);
    
    // 建立反向映射
    cellToBinMap_[ndi] = bucket;

}

bool DetailedMis3D::gatherNeighbours3D(Node3D* ndi)
{
    // 获取单行高度（用于计算行跨度）
    Architecture3D::LayerType cell_layer_type = Architecture3D::stringToLayer(ndi->getLayer());
    const double singleRowHeight = arch_->getRow(cell_layer_type, 0)->getHeight().v;

    neighbours_.clear();
    neighbours_.push_back(ndi);

    // 扫描网格结构，收集与当前单元兼容的单元

    auto it = cellToBinMap_.find(ndi);
    if (it == cellToBinMap_.end()) {
        return false;
    }

    const int spanned_i = std::lround(ndi->getHeight().v / singleRowHeight);

    // 根据单元所在层选择对应的网格
    std::vector<std::vector<Bucket3D*>>* grid_ptr = nullptr;
    int dimX = 0, dimY = 0;
    
    if (ndi->getLayer() == "bottom") {
        grid_ptr = &grid_bottom_;
        dimX = dimX_bottom_;
        dimY = dimY_bottom_;
    } else if (ndi->getLayer() == "upper") {
        grid_ptr = &grid_upper_;
        dimX = dimX_upper_;
        dimY = dimY_upper_;
    } else {
        return false;  // 未知层
    }
    
    std::vector<std::vector<Bucket3D*>>& grid = *grid_ptr;

    std::queue<Bucket3D*> Q;
    Q.push(it->second);
    ++traversal_;
    while (!Q.empty()) {
        Bucket3D* currPtr = Q.front();
        Q.pop();

        if (currPtr->travId == traversal_) {
            continue;
        }
        currPtr->travId = traversal_;

        // 扫描此桶中的所有单元。如果它们与原始单元兼容，则添加到邻居列表。
        for (Node3D* ndj : currPtr->nodes) {
            // 检查确保单元不是原始单元，它们具有相同的区域，相同的尺寸（如果适用），
            // 以及相同的颜色（如果适用）。

            // 不同的单元
            if (ndj == ndi) {
                continue;
            }

            // 3D特定：必须位于同一层
            if (ndj->getLayer() != ndi->getLayer()) {

                continue;
                // 如果需要跨层搜索，可以特殊处理
                // if (!params_.allowCrossLayer) {
                //     continue;
                // }
                // // 检查层兼容性
                // if (!checkLayerCompatibility(ndi, ndj)) {
                //     continue;
                // }
            }

            // 必须具有相同的颜色以避免共享网络（仅HPWL目标）
            if (useSameColor_ && colors_[ndi->getId()] != colors_[ndj->getId()]) {
                continue;
            }

            // 必须具有相同的尺寸
            if (useSameSize_
                && (ndi->getWidth() != ndj->getWidth()
                    || ndi->getHeight() != ndj->getHeight())) {
                continue;
            }

            // 必须在同一区域
            if (ndj->getGroupId() != ndi->getGroupId()) {
                continue;
            }

            // 必须跨越相同的行数并且电压兼容
            const int spanned_j = std::lround(ndj->getHeight().v / singleRowHeight);
            if (ndi->getBottomPower() != ndj->getBottomPower()
                || ndi->getTopPower() != ndj->getTopPower()
                || spanned_i != spanned_j) {
                continue;
            }

            // 3D特定：检查3D约束（HBT对齐等）
            // if (!check3DConstraints(ndi, ndj)) {
            //     continue;
            // }

            // 如果兼容，包含此当前单元
            neighbours_.push_back(ndj);
        }

        if (neighbours_.size() >= maxProblemSize_) {
            break;
        }

        // 如果尚未收集足够的单元，将更多桶添加到队列
        if (currPtr->i > 0) {
            Q.push(grid[currPtr->i - 1][currPtr->j]);
        }
        if (currPtr->i + 1 < dimX) {
            Q.push(grid[currPtr->i + 1][currPtr->j]);
        }
        if (currPtr->j > 0) {
            Q.push(grid[currPtr->i][currPtr->j - 1]);
        }
        if (currPtr->j + 1 < dimY) {
            Q.push(grid[currPtr->i][currPtr->j + 1]);
        }
    }
    return true;
}

void DetailedMis3D::solveMatch3D()
{
    if (neighbours_.size() <= 1) {
        return;
    }
    const std::vector<Node3D*>& nodes = neighbours_;

    const int nNodes = (int) nodes.size();
    const int nSpots = (int) nodes.size();

    if (!segment_manager_) {
        logger_->error(TDP, 1404, "Segment manager not available");
        return;
    }

    // 保存每个单元的原始位置和段信息
    std::vector<Pointdp3D> original_positions(nNodes);
    std::vector<std::vector<Segment3D*>> original_segments(nNodes);
    
    for (size_t i = 0; i < nodes.size(); i++) {
        Node3D* ndi = nodes[i];
        original_positions[i] = ndi->getLocation();
        original_segments[i] = segment_manager_->getCellSegments(ndi);
    }

    // ======= STA EVALUATION PIPELINE SETUP =======
    SteinerTreeGenerator steiner(timingdp_);
    auto* global_router = timingdp_->getGlobalRouter();
    auto& global_routes = global_router->getRoutes();

    std::set<odb::dbNet*> affected_nets;
    for (auto* nd : nodes) {
        for (auto* db_net : nd->connected_nets) {
            if (db_net) affected_nets.insert(db_net);
        }
    }

    std::map<odb::dbNet*, grt::GRoute> saved_routes;
    for (auto* db_net : affected_nets) {
        auto it = global_routes.find(db_net);
        if (it != global_routes.end())
            saved_routes[db_net] = it->second;
    }

    // ======= PRE-COMPUTE COST MATRIX (full STA per pair) =======
    const double INF_COST = std::numeric_limits<double>::max();
    std::vector<std::vector<double>> cost_matrix(
        nNodes, std::vector<double>(nSpots, INF_COST));

    for (int i = 0; i < nNodes; i++) {
        for (int j = 0; j < nSpots; j++) {
            // if (i == j) {
            //     logger_ -> info (TDP, 1499, "debug cost mis");
            //     // cost_matrix[i][j] = 0.0;
            //     // continue;
            // }
            const Pointdp3D& target_pos = original_positions[j];
            if (!checkDisplacementConstraint3D(nodes[i], target_pos)) {
                continue;
            }
            cost_matrix[i][j] = evalCostAtPos(nodes[i], target_pos, steiner, saved_routes);
        }
    }

    // ======= RESTORE ALL STATE =======
    for (const auto& [db_net, route] : saved_routes)
        global_routes[db_net] = route;
    for (size_t i = 0; i < nodes.size(); i++){
        nodes[i]->setLocation(original_positions[i]);
        timingdp_->setInstanceLocation3D( 
            nodes[i]->db_inst,
            original_positions[i].x.v,
            original_positions[i].y.v,
            original_positions[i].z.v);
    }

    // ======= BUILD MIN-COST FLOW GRAPH =======
    lemon::ListDigraph g;
    std::vector<lemon::ListDigraph::Node> nodeForCell(nNodes);
    std::vector<lemon::ListDigraph::Node> nodeForSpot(nNodes);
    
    for (int i = 0; i < nNodes; i++) {
        nodeForCell[i] = g.addNode();
        nodeForSpot[i] = g.addNode();
    }
    
    lemon::ListDigraph::Node supplyNode = g.addNode();
    lemon::ListDigraph::Node demandNode = g.addNode();

    lemon::ListDigraph::ArcMap<int> lower_bound(g);
    lemon::ListDigraph::ArcMap<int> upper_bound(g);
    lemon::ListDigraph::ArcMap<int> cost(g);
    
    std::map<lemon::ListDigraph::Arc, std::pair<int, int>> reverse_map;

    for (int i = 0; i < nNodes; i++) {
        // source -> cell
        auto arc_s2c = g.addArc(supplyNode, nodeForCell[i]);
        lower_bound[arc_s2c] = 0;
        upper_bound[arc_s2c] = 1;
        cost[arc_s2c] = 0;
        
        // spot -> sink
        auto arc_s2d = g.addArc(nodeForSpot[i], demandNode);
        lower_bound[arc_s2d] = 0;
        upper_bound[arc_s2d] = 1;
        cost[arc_s2d] = 0;
        
        // cell -> spot
        for (int j = 0; j < nSpots; j++) {
            if (cost_matrix[i][j] == INF_COST) continue;

            int int_cost = (cost_matrix[i][j] > std::numeric_limits<int>::max())
                           ? std::numeric_limits<int>::max()
                           : static_cast<int>(cost_matrix[i][j]);
            
            auto arc_c2s = g.addArc(nodeForCell[i], nodeForSpot[j]);
            lower_bound[arc_c2s] = 0;
            upper_bound[arc_c2s] = 1;
            cost[arc_c2s] = int_cost;
            reverse_map[arc_c2s] = std::make_pair(i, j);
        }
    }
    
    // ======= SOLVE FLOW =======
    lemon::Preflow<lemon::ListDigraph> preflow(g, upper_bound, supplyNode, demandNode);
    preflow.run();
    const int max_flow = preflow.flowValue();
    
    if (max_flow != nNodes) {
        return;
    }
    
    lemon::NetworkSimplex<lemon::ListDigraph> mincost(g);
    mincost.lowerMap(lower_bound);
    mincost.upperMap(upper_bound);
    mincost.costMap(cost);
    mincost.stSupply(supplyNode, demandNode, max_flow);
    
    lemon::NetworkSimplex<lemon::ListDigraph>::ProblemType result = mincost.run();
    if (result != lemon::NetworkSimplex<lemon::ListDigraph>::OPTIMAL) {
        logger_->info(TDP, 1405, "Min-cost flow not optimal");
        return;
    }
    
    // ======= APPLY MATCHING =======
    lemon::ListDigraph::ArcMap<int> flow(g);
    mincost.flowMap(flow);
    
    bool has_violation = false;
    
    for (lemon::ListDigraph::ArcMap<int>::ItemIt it(flow); it != lemon::INVALID; ++it) {
        if (g.target(it) != demandNode && g.source(it) != supplyNode && mincost.flow(it) != 0) {
            auto map_it = reverse_map.find(it);
            if (map_it == reverse_map.end()) {
                logger_->error(TDP, 1406, "Cannot find reverse mapping for flow arc");
                continue;
            }
            
            const int cell_idx = map_it->second.first;
            const int spot_idx = map_it->second.second;
            
            if (cell_idx == spot_idx) {
                continue;
            }
            
            Node3D* moving_cell = nodes[cell_idx];
            
            // Remove from original segments
            for (Segment3D* seg : original_segments[cell_idx]) {
                if (!segment_manager_->removeCellFromSegment(moving_cell, seg->getId())) {
                    logger_->error(TDP, 1407, "Failed to remove cell from segment " + std::to_string(seg->getId()));
                    has_violation = true;
                    break;
                }
            }
            if (has_violation) break;
            
            // Update position
            Pointdp3D old_position = moving_cell->getLocation();
            Pointdp3D new_position = original_positions[spot_idx];
            moving_cell->setLocation(new_position);
            timingdp_->setInstanceLocation3D(moving_cell->db_inst, new_position.x.v, new_position.y.v, new_position.z.v);
            
            // Add to new segments
            std::string cell_layer = moving_cell->getLayer();
            std::vector<Segment3D*> new_segments = 
                segment_manager_->getSegmentsAtPosition(new_position, cell_layer);
            
            bool added_to_segment = false;
            for (Segment3D* seg : new_segments) {
                if (segment_manager_->addCellToSegment(moving_cell, seg->getId())) {
                    added_to_segment = true;
                    break;
                }
            }
            
            if (!added_to_segment) {
                logger_->error(TDP, 1408, "Failed to add cell to any segment at new position");
                moving_cell->setLocation(old_position);
                timingdp_->setInstanceLocation3D(moving_cell->db_inst, old_position.x.v, old_position.y.v, old_position.z.v);
                for (Segment3D* seg : original_segments[cell_idx]) {
                    segment_manager_->addCellToSegment(moving_cell, seg->getId());
                }
                has_violation = true;
                break;
            }
            
            timesUsed_[moving_cell->getId()]++;
        }
    }
    
    // ======= VALIDATE OR ROLLBACK =======
    if (has_violation) {
        logger_->info(TDP, 1401, "Placement violation detected, rolling back match");
        for (size_t i = 0; i < nodes.size(); i++) {
            Node3D* ndi = nodes[i];
            auto current_segments = segment_manager_->getCellSegments(ndi);
            for (Segment3D* seg : current_segments) {
                segment_manager_->removeCellFromSegment(ndi, seg->getId());
            }
            for (Segment3D* seg : original_segments[i]) {
                segment_manager_->addCellToSegment(ndi, seg->getId());
            }
            ndi->setLocation(original_positions[i]);
            timingdp_->setInstanceLocation3D(ndi->db_inst, original_positions[i].x.v, original_positions[i].y.v, original_positions[i].z.v);
        }
    } else {
        // logger_->info(TDP, 1402, "Matching applied successfully for " + std::to_string(nodes.size()) + " cells");
        segment_manager_->validateSegments();
    }
}

//place 函数
void DetailedMis3D::place3D()
{
    // 1. 填充3D网格
    populateGrid3D();
    

    timesUsed_.resize(network_->getNumCells());
    std::ranges::fill(timesUsed_, 0);
    // 2. 初始化统计
    int total_matches = 0;
    
    // 3. 遍历候选单元（随机顺序以避免偏差）
    std::vector<Node3D*> shuffled_candidates = candidates_;
    std::shuffle(shuffled_candidates.begin(), shuffled_candidates.end(), 
                 std::default_random_engine(params_.random_seed));

///////debug
    int total_candidates = shuffled_candidates.size();
    int processed = 0;
    int last_pct = -1;
    
///////debug    
    for (Node3D* ndi : shuffled_candidates) {
///////
        processed++;
        int pct = processed * 100 / total_candidates;
        if (pct != last_pct) {
            std::cout << "MIS3D progress: " << pct << "% ("
                      << processed << "/" << total_candidates << ")"
                      << "  matched: " << total_matches << std::endl;
            last_pct = pct;
            double L_lambda = 0.0;
            for (auto* cell : candidates_) {
                L_lambda += lsp3d_->calculateCellCost(cell, 5.0, 1.0);
            }
            std::cout << "L(λ) = " << L_lambda  << " for " << candidates_.size() << " cells" << std::endl;
        }


////////       
        // 检查单元使用次数
        if (timesUsed_[ndi->getId()] >= params_.maxTimesUsed) {
            continue;
        }
        
        // 4. 收集3D兼容的邻居单元
        if (!gatherNeighbours3D(ndi)) {
            continue;
        }
        
        // 如果邻居太少，跳过
        if (neighbours_.size() < params_.minNeighbours) {
            continue;
        }
        
        // 5. 求解3D匹配问题
        solveMatch3D();
        // lsp3d_->updatePinWeights(); 
        
        total_matches++;
        // if (result.success) {
        //     // 更新单元使用计数
        //     timesUsed_[ndi->getId()]++;
        //     for (Node3D* matched_cell : result.matched_cells) {
        //         timesUsed_[matched_cell->getId()]++;
        //     }
        // }
    }
    
    // 6. 更新统计信息
    // updateStatistics3D(total_matches, successful_matches, 
    //                   total_hpwl_improvement, total_disp_improvement);
    
}

////////////////////////////////////////////////////////////////////////////////
// 成本计算函数 - 需要您实现
////////////////////////////////////////////////////////////////////////////////

uint64_t DetailedMis3D::gettiming3D(const Node3D* ndi) {
    if (!lsp3d_) return 0;
    return static_cast<uint64_t>(
        lsp3d_->calculateCellCost(const_cast<Node3D*>(ndi), 5.0, 1.0));
}

// uint64_t DetailedMis3D::getDisp3D(const Node3D* ndi, DbuX xi, DbuY yi, int zi) {
//     // TODO: 计算3D位移成本
//     return 0;
// }

////////////////////////////////////////////////////////////////////////////////
// 3D特定检查函数 
////////////////////////////////////////////////////////////////////////////////

// bool DetailedMis3D::checkLayerCompatibility(const Node3D* cell1, const Node3D* cell2) const {
//     // TODO: 检查两个单元是否可以在同一层或相邻层
//     return false;
// }

// bool DetailedMis3D::check3DConstraints(const Node3D* cell1, const Node3D* cell2) const {
//     // TODO: 检查3D约束（层差异、HBT对齐等）
//     return false;
// }

// Pointdp3D DetailedMis3D::getLegalPosition3D(const Node3D* cell, const Pointdp3D& target) const {
//     // TODO: 计算合法的3D位置
//     return target;
// }

// bool DetailedMis3D::isValid3DPosition(const Node3D* cell, const Pointdp3D& pos) const {
//     // TODO: 验证3D位置是否有效
//     return false;
// }

// double DetailedMis3D::getHBTAlignmentCost(const Node3D* cell, const Pointdp3D& pos) const {
//     // TODO: 计算HBT对齐成本
//     return 0.0;
// }

////////////////////////////////////////////////////////////////////////////////
// 工具函数
////////////////////////////////////////////////////////////////////////////////
bool DetailedMis3D::checkDisplacementConstraint3D(const Node3D* cell, const Pointdp3D& target) const 
{
    // 获取段管理器（类似DPL的DetailedMgr）
    
    if (!segment_manager_) {
        return false;
    }
    
    // 获取当前单元位置
    Pointdp3D current = cell->getLocation();
    
    // 计算3D位移
    DbuX dx = std::abs((target.x - current.x).v);
    DbuY dy = std::abs((target.y - current.y).v);

    
    // 获取最大位移限制
    Architecture3D::LayerType layer = 
    (cell->getLayer() == "bottom") ? Architecture3D::LAYER_BOTTOM : Architecture3D::LAYER_UPPER;
    int limit = std::max(arch_->getWidth(layer).v, arch_->getHeight(layer).v) << 1;
    int maxDispX = limit;
    int maxDispY = limit;
    // 检查X和Y方向的位移（类似DPL的实现）
    if (dx.v > maxDispX) {
        return false;
    }
    if (dy.v > maxDispY) {
        return false;
    }
    
    return true;
}

// Lambda: evaluate timing cost for cell at target position (full STA pipeline)
double DetailedMis3D::evalCostAtPos(Node3D* cell, const Pointdp3D& target_pos,
                                     SteinerTreeGenerator& steiner,
                                     std::map<odb::dbNet*, grt::GRoute>& saved_routes) {
    if (!lsp3d_) return 0.0;
    Pointdp3D saved_pos = cell->getLocation();
    // Restore global routes to saved state
    auto* global_router = timingdp_->getGlobalRouter();
    auto& global_routes = global_router->getRoutes();
    for (const auto& [db_net, route] : saved_routes)
        global_routes[db_net] = route;

    // Temporarily move cell
    cell->setLocation(target_pos);
    odb::dbInst* inst = cell->db_inst;
    timingdp_->setInstanceLocation3D(
        inst, target_pos.x.v, target_pos.y.v, target_pos.z.v);

    // Update Steiner trees
    steiner.updateSteinerTrees(cell);

    // Incremental STA
    timingdp_->executeTimingDrivenI(cell);

    // Compute cost via LSP3D
    double cost = lsp3d_->calculateCellCost(cell, 5.0, 1.0);

    timingdp_->setInstanceLocation3D(
        inst, saved_pos.x.v, saved_pos.y.v, saved_pos.z.v);
    cell->setLocation(saved_pos);
    // steiner.updateSteinerTrees(cell);
    return cost;
}
// void DetailedMis3D::updateStatistics() {
//     // TODO: 更新算法统计信息
// }

// void DetailedMis3D::logInfo(const std::string& msg) const {
//     // TODO: 日志记录
// }

// void DetailedMis3D::placer_db_->logDebug(const std::string& msg) const {
//     // TODO: 调试日志
// }

} // namespace tdp
