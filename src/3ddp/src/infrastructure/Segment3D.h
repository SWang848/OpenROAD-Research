// include/Segment3D.h
#pragma once

#include "Common3D.h"
#include "Node3D.h"
#include <vector>
#include <string>

namespace tdp {

class Node3D;

/**
 * @brief 3D布局中的段（Segment）类
 * 
 * 段是行中的连续空间，用于管理单元的位置。
 * 在3D布局中，每层有独立的段系统。
 */
class Segment3D {
public:
    Segment3D(int id, const std::string& layer, int row_id, 
              DbuX min_x, DbuX max_x,
              DbuY row_ymin = DbuY{0}, DbuY row_ymax = DbuY{0});
    
    // 基本访问器
    int getId() const { return id_; }
    std::string getLayer() const { return layer_; }
    int getRowId() const { return row_id_; }
    DbuX getMinX() const { return min_x_; }
    DbuX getMaxX() const { return max_x_; }
    DbuY getRowYMin() const { return row_ymin_; }
    DbuY getRowYMax() const { return row_ymax_; }
    DbuX getWidth() const { return max_x_ - min_x_; }
    
    // 单元管理
    void addCell(Node3D* cell);
    void removeCell(Node3D* cell);
    const std::vector<Node3D*>& getCells() const { return cells_; }
    
    // 容量检查
    bool canAccommodate(const Node3D* cell) const;
    int getUtilization() const { return utilization_; }
    int getAvailableWidth() const { return getWidth().v - utilization_; }
    double getUtilizationRatio() const;
    
    // 位置检查
    bool contains(const Pointdp3D& pos) const;
    bool containsX(DbuX x) const;
    
    // 调试信息
    std::string toString() const;
    
private:
    int id_;                    // 段ID
    std::string layer_;         // 所属层（"bottom" 或 "upper"）
    int row_id_;               // 行ID
    DbuX min_x_, max_x_;       // X范围
    DbuY row_ymin_, row_ymax_; // 所在行的Y范围
    std::vector<Node3D*> cells_; // 段中的单元
    int utilization_;          // 已使用宽度（数据库单位）
};

} // namespace tdp
