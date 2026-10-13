// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <deque>

#include "../infrastructure/Coordinates.h"
#include "../infrastructure/SegmentManager3D.h"
#include "../eva/eva.h"
#include "utl/Logger.h"

// 前向声明（全局命名空间）
namespace sta { class dbSta; }

namespace tdp {

// 前向声明
class Node3D;
class Architecture3D;
class Network3D;
class Timingdp;
class LSP3D;
class TimingEngine3D;

// 3D MIS 算法参数结构体
struct DetailedMis3DParams
{
public:
    enum Strategy
    {
        KDTree = 0,
        Binning = 1,
        Colour = 2,
    };

    double maxDifferenceInMetric = 0.03;  // 允许重新引入重叠的程度
    unsigned maxNumNodes = 15;           // B&B 考虑的最大节点数
    unsigned maxPasses = 10;              // 最大B&B迭代次数
    double sizeTol = 1.1;                // 相同尺寸的容差
    unsigned skipNetsLargerThanThis = 50; // 跳过大于此值的网络
    Strategy strategy = Binning;          // 使用的策略类型
    bool useSameSize = true;              // 是否只与相同尺寸的单元交换
    bool allowCrossLayer = false;         // 是否允许跨层移动
    int maxDisplacementX = 0;             // X方向最大位移（Dbu）
    int maxDisplacementY = 0;             // Y方向最大位移（Dbu）
    int maxDisplacementZ = 0;             // Z方向最大位移（层数）
    int maxTimesUsed = 20;                 // 单元最大使用次数
    int minNeighbours = 2;                // 最小邻居数
    int random_seed = 42;                 // 随机种子
};

// 3D桶结构体
struct Bucket3D
{
    void clear() { nodes.clear(); }

    std::deque<Node3D*> nodes;  // 桶中的单元
    double xmin = 0.0;          // X最小值
    double xmax = 0.0;          // X最大值
    double ymin = 0.0;          // Y最小值
    double ymax = 0.0;          // Y最大值
    int i = 0;                  // X网格索引
    int j = 0;                  // Y网格索引
    int k = 0;                  // Z网格索引（层索引）
    int travId = 0;             // 遍历ID
    std::string layer;          // 层名称
};

// 3D匹配结果结构体
struct MatchResult3D
{
    bool success = false;
    std::vector<Node3D*> matched_cells;
    double timing_improvement = 0.0;
    // double disp_improvement = 0.0;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class DetailedMis3D
{
public:
    // 目标函数类型
    // enum Objective
    // {
    //     timing,   // 半周长线长
    //     Disp,   // 位移
    //     Timing  // 时序
    // };

    // 构造函数和析构函数
    DetailedMis3D(Timingdp* timingdp, LSP3D* lsp3d = nullptr);
    virtual ~DetailedMis3D();
    void run( );

private:
    // 核心算法函数
    void place3D();
    void collectMovableCells();
    void colorCells3D();
    void buildGrid3D();
    void clearGrid3D();
    void populateGrid3D();
    void populateCellToLayerGrid(Node3D* ndi,
                                 std::vector<std::vector<Bucket3D*>>& grid,
                                 double xmin, double ymin,
                                 double stepX, double stepY,
                                 int dimX, int dimY);
    bool gatherNeighbours3D(Node3D* ndi);
    //bool gatherNeighbours3D(Node3D* ndi, std::vector<Node3D*>& neighbours);
    void solveMatch3D();
    //MatchResult3D solveMatch3D(Node3D* ndi, std::vector<Node3D*>& neighbours);

    // 成本计算函数
    // uint64_t getDisp3D(const Node3D* ndi, DbuX xi, DbuY yi, int zi);
    uint64_t gettiming3D(const Node3D* ndi);

    // 3D特定检查函数
    // bool checkLayerAssignmentCompatibility(const Node3D* cell1, const Node3D* cell2) const;
    bool checkDisplacementConstraint3D(const Node3D* cell, const Pointdp3D& target) const;
    // bool hasPlacementViolation3D(const Node3D* cell) const;
    // double get3DConstraintPenalty(const Node3D* cell, const Pointdp3D& pos) const;
    // double getLayerZCoordinate(const std::string& layer) const;
    // bool moveCell3D(Node3D* cell, const Pointdp3D& new_position);

    // 工具函数
    // void updateStatistics();
    // void updateMatchStatistics(const std::vector<Node3D*>& nodes);
    // void logInfo(const std::string& msg) const;
    // void logDebug(const std::string& msg) const;
    // void logError(const std::string& msg) const;
    double evalCostAtPos(Node3D* cell, const Pointdp3D& target_pos,
                         SteinerTreeGenerator& steiner,
                         std::map<odb::dbNet*, grt::GRoute>& saved_routes);

    // 内部数据结构
    // struct ColorGraph
    // {
    //     ColorGraph(int n) : adj(n), color(n, -1) {}
    //     void addEdge(int u, int v) { adj[u].push_back(v); adj[v].push_back(u); }
    //     void greedyColoring();
    //     int getColor(int i) const { return color[i]; }
        
    //     std::vector<std::vector<int>> adj;
    //     std::vector<int> color;
    // };

    // 参数和状态
    DetailedMis3DParams params_;
    Architecture3D* arch_;
    Network3D* network_;
    utl::Logger* logger_;
    Timingdp* timingdp_ = nullptr;
    LSP3D* lsp3d_ = nullptr;
    sta::dbSta* sta_ = nullptr;
    TimingEngine3D* timing_engine_ = nullptr;
    SegmentManager3D* segment_manager_ = nullptr;

    // 候选单元和邻居
    std::vector<Node3D*> candidates_;
    std::vector<bool> movable_;
    std::vector<int> colors_;
    std::vector<Node3D*> neighbours_;
    std::vector<int> timesUsed_;

    // 3D网格结构（分离的底层和顶层网格）
    std::vector<std::vector<Bucket3D*>> grid_bottom_;  // 底层网格
    std::vector<std::vector<Bucket3D*>> grid_upper_;   // 顶层网格
    std::map<Node3D*, Bucket3D*> cellToBinMap_;        // 单元到桶的映射

    // 网格参数
    int dimX_bottom_ = 0;
    int dimY_bottom_ = 0;
    int dimX_upper_ = 0;
    int dimY_upper_ = 0;
    double stepX_bottom_ = 0.0;
    double stepY_bottom_ = 0.0;
    double stepX_upper_ = 0.0;
    double stepY_upper_ = 0.0;

    // 其他参数
    int skipNetsLargerThanThis_ = 100;
    int maxProblemSize_ = 25;
    int traversal_ = 0;
    bool useSameSize_ = true;
    bool useSameColor_ = true;
    // Objective obj_ = DetailedMis3D::Hpwl;
};

} // namespace tdp
