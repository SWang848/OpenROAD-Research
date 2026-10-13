// src/Segment3D.cpp
#include "Segment3D.h"
#include "Node3D.h"
#include <algorithm>
#include <sstream>

namespace tdp {

Segment3D::Segment3D(int id, const std::string& layer, int row_id,
                     DbuX min_x, DbuX max_x,
                     DbuY row_ymin, DbuY row_ymax)
    : id_(id), layer_(layer), row_id_(row_id),
      min_x_(min_x), max_x_(max_x),
      row_ymin_(row_ymin), row_ymax_(row_ymax),
      utilization_(0) {
    if (min_x_ >= max_x_) {
        throw std::runtime_error("Segment3D: min_x must be less than max_x");
    }
}

void Segment3D::addCell(Node3D* cell) {
    if (!cell) return;
    
    // 检查单元是否已在段中
    auto it = std::find(cells_.begin(), cells_.end(), cell);
    if (it != cells_.end()) {
        // 单元已在段中，可能是错误或需要更新
        return;
    }

    // 按X坐标排序插入
    DbuX cell_center = cell->location.x + cell->getWidth() / 2;
    auto insert_pos = std::lower_bound(cells_.begin(), cells_.end(), cell,
        [](Node3D* a, Node3D* b) {
            return a->location.x < b->location.x;
        });
    
    cells_.insert(insert_pos, cell);
    utilization_ += cell->getWidth().v;
}

void Segment3D::removeCell(Node3D* cell) {
    if (!cell) return;
    
    auto it = std::find(cells_.begin(), cells_.end(), cell);
    if (it != cells_.end()) {
        cells_.erase(it);
        utilization_ -= cell->getWidth().v;
        if (utilization_ < 0) utilization_ = 0;
    }
}

bool Segment3D::canAccommodate(const Node3D* cell) const {
    if (!cell) return false;
    
    // 检查宽度
    if (cell->getWidth().v > getAvailableWidth()) {
        return false;
    }
    
    // 检查层兼容性
    if (cell->layer != layer_) {
        return false;
    }
    
    // TODO: 检查其他3D约束（HBT对齐等）
    
    return true;
}

bool Segment3D::contains(const Pointdp3D& pos) const {
    if (!containsX(pos.x)) return false;
    // if (pos.y != row_ymin_ || pos.y > row_ymax_) 
    if ( pos.y != row_ymin_ ) 
        return false;
    return true;
}

bool Segment3D::containsX(DbuX x) const {
    return x >= min_x_ && x <= max_x_;
}

double Segment3D::getUtilizationRatio() const {
    DbuX width = getWidth();
    if (width.v == 0) return 0.0;
    return static_cast<double>(utilization_) / width.v;
}

std::string Segment3D::toString() const {
    std::stringstream ss;
    ss << "Segment3D[id=" << id_ 
       << ", layer=" << layer_
       << ", row=" << row_id_
       << ", x=[" << min_x_.v << "," << max_x_.v << "]"
       << ", y=[" << row_ymin_.v << "," << row_ymax_.v << "]"
       << ", width=" << getWidth().v
       << ", util=" << utilization_
       << ", cells=" << cells_.size() << "]";
    return ss.str();
}

} // namespace tdp
