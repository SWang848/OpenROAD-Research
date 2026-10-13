// src/SegmentManager3D.cpp
#include "SegmentManager3D.h"
#include "Architecture3D.h"
#include "Node3D.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace tdp {

// SegmentManager3D::SegmentManager3D(Architecture3D* arch)
//     : arch_(arch), next_segment_id_(0) {
//     if (!arch_) {
//         throw std::runtime_error("SegmentManager3D: Architecture cannot be null");
//     }
// }

SegmentManager3D::~SegmentManager3D() {
    clearSegments();
}

void SegmentManager3D::initializeSegments() {
    clearSegments();
    
    // 为每层创建段
    createSegmentsForLayer("bottom");
    createSegmentsForLayer("upper");
    
    std::cout << "SegmentManager3D: Initialized " << getNumSegments() 
              << " segments" << std::endl;
}

void SegmentManager3D::clearSegments() {
    segments_.clear();
    layer_segments_.clear();
    segments_in_row_.clear();
    cell_to_segments_.clear();
    next_segment_id_ = 0;
}

void SegmentManager3D::createSegmentsForLayer(const std::string& layer) {
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(layer);
    
    // 获取该层的行数
    int num_rows = arch_->getNumRows(layer_type);
    if (num_rows <= 0) {
        std::cerr << "Warning: No rows found for layer " << layer << std::endl;
        return;
    }
    
    // 初始化数据结构
    layer_segments_[layer] = std::vector<Segment3D*>();
    segments_in_row_[layer] = std::map<int, std::vector<Segment3D*>>();
    
    // 为每行创建段
    for (int row_id = 0; row_id < num_rows; row_id++) {
        createSegmentsForRow(layer, row_id);
    }
}

void SegmentManager3D::createSegmentsForRow(const std::string& layer, int row_id) {
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(layer);
    
    // 获取行信息
    const auto* row = arch_->getRow(layer_type, row_id);
    if (!row) {
        std::cerr << "Error: Row " << row_id << " not found in layer " << layer << std::endl;
        return;
    }
    
    DbuX row_min_x = row->getLeft();
    DbuX row_max_x = row->getRight();
    DbuY row_ymin = row->getBottom();
    DbuY row_ymax = row->getTop();
    
    // 获取该行的阻塞（简化版本，实际需要从Architecture3D获取）
    auto blockages = getBlockagesForRow(layer, row_id);
    
    if (blockages.empty()) {
        // 整行空闲，创建一个段
        createSingleSegment(layer, row_id, row_min_x, row_max_x, row_ymin, row_ymax);
    } else {
        // 根据阻塞划分段（类似DPL的逻辑）
        // 对阻塞按X排序
        std::sort(blockages.begin(), blockages.end(),
            [](const Blockage3D& a, const Blockage3D& b) {
                return a.x_min < b.x_min;
            });
        
        // 阻塞前的段
        if (blockages[0].x_min > row_min_x) {
            createSingleSegment(layer, row_id, row_min_x, blockages[0].x_min, row_ymin, row_ymax);
        }
        
        // 阻塞间的段
        for (size_t i = 1; i < blockages.size(); i++) {
            if (blockages[i].x_min > blockages[i-1].x_max) {
                createSingleSegment(layer, row_id, 
                                   blockages[i-1].x_max, blockages[i].x_min, row_ymin, row_ymax);
            }
        }
        
        // 阻塞后的段
        if (blockages.back().x_max < row_max_x) {
            createSingleSegment(layer, row_id, 
                               blockages.back().x_max, row_max_x, row_ymin, row_ymax);
        }
    }
}

void SegmentManager3D::createSingleSegment(const std::string& layer, int row_id,
                                          DbuX min_x, DbuX max_x,
                                          DbuY row_ymin, DbuY row_ymax) {
    if (min_x >= max_x) {
        return;  // 无效范围
    }
    
    // 创建新段
    auto segment = std::make_unique<Segment3D>(
        next_segment_id_++, layer, row_id, min_x, max_x, row_ymin, row_ymax);
    
    Segment3D* segment_ptr = segment.get();
    
    // 添加到数据结构
    segments_.push_back(std::move(segment));
    layer_segments_[layer].push_back(segment_ptr);
    segments_in_row_[layer][row_id].push_back(segment_ptr);
}

std::vector<Segment3D*> SegmentManager3D::getCellSegments(const Node3D* cell) const {
    if (!cell) return {};
    
    auto it = cell_to_segments_.find(cell->getId());
    if (it != cell_to_segments_.end()) {
        return it->second;
    }
    
    return {};
}

std::vector<Segment3D*> SegmentManager3D::getSegmentsAtPosition(
    const Pointdp3D& pos, const std::string& layer) const {
    
    std::vector<Segment3D*> result;
    
    // 查找包含该位置的段
    auto layer_it = layer_segments_.find(layer);
    if (layer_it == layer_segments_.end()) {
        return result;
    }
    
    for (Segment3D* seg : layer_it->second) {
        if (seg->contains(pos)) {
            result.push_back(seg);
        }
    }
    
    return result;
}

std::vector<Segment3D*> SegmentManager3D::getSegmentsInRow(
    const std::string& layer, int row_id) const {
    
    std::vector<Segment3D*> result;
    
    auto layer_it = segments_in_row_.find(layer);
    if (layer_it == segments_in_row_.end()) {
        return result;
    }
    
    auto row_it = layer_it->second.find(row_id);
    if (row_it != layer_it->second.end()) {
        result = row_it->second;
    }
    
    return result;
}

bool SegmentManager3D::addCellToSegment(Node3D* cell, int seg_id) {
    if (!cell) return false;
    
    Segment3D* segment = findSegmentById(seg_id);
    if (!segment) {
        std::cerr << "Error: Segment " << seg_id << " not found" << std::endl;
        return false;
    }
    
    // 检查兼容性
    if (!segment->canAccommodate(cell)) {
        // std::cerr << "Error: Segment " << seg_id << " cannot accommodate cell " 
        //           << cell->getId() << std::endl;
        return false;
    }
    
    // 添加到段
    segment->addCell(cell);
    
    // 更新反向映射
    cell_to_segments_[cell->getId()].push_back(segment);
    
    return true;
}

bool SegmentManager3D::removeCellFromSegment(const Node3D* cell, int seg_id) {
    if (!cell) return false;
    
    Segment3D* segment = findSegmentById(seg_id);
    if (!segment) {
        return false;
    }
    
    // 从段中移除
    segment->removeCell(const_cast<Node3D*>(cell));
    
    // 从反向映射中移除
    auto it = cell_to_segments_.find(cell->getId());
    if (it != cell_to_segments_.end()) {
        auto& seg_list = it->second;
        seg_list.erase(std::remove(seg_list.begin(), seg_list.end(), segment),
                      seg_list.end());
        
        if (seg_list.empty()) {
            cell_to_segments_.erase(it);
        }
    }
    
    return true;
}

bool SegmentManager3D::moveCellToSegment(Node3D* cell, int new_seg_id) {
    if (!cell) return false;
    
    // 获取当前段
    auto current_segments = getCellSegments(cell);
    
    // 从所有当前段中移除
    for (Segment3D* seg : current_segments) {
        removeCellFromSegment(cell, seg->getId());
    }
    
    // 添加到新段
    return addCellToSegment(cell, new_seg_id);
}

Segment3D* SegmentManager3D::findSegmentById(int seg_id) const {
    for (const auto& segment : segments_) {
        if (segment->getId() == seg_id) {
            return segment.get();
        }
    }
    return nullptr;
}

Segment3D* SegmentManager3D::findClosestSegment(const Node3D* cell) const {
    if (!cell) return nullptr;
    
    std::string cell_layer = cell->getLayer();
    Pointdp3D cell_pos = cell->getLocation();
    
    // 查找同一层的段
    auto layer_it = layer_segments_.find(cell_layer);
    if (layer_it == layer_segments_.end()) {
        return nullptr;
    }
    
    Segment3D* closest_segment = nullptr;
    double min_distance = std::numeric_limits<double>::max();
    
    for (Segment3D* seg : layer_it->second) {
        // 检查区域兼容性（简化版本）
        // TODO: 实现完整的区域检查
        
        // 检查段是否能容纳单元
        if (!seg->canAccommodate(cell)) {
            continue;
        }
        
        // 计算距离（使用段中心）
        DbuX seg_center_x = (seg->getMinX() + seg->getMaxX()) / 2;
        double distance = std::abs(seg_center_x.v - cell_pos.x.v);
        
        if (distance < min_distance) {
            min_distance = distance;
            closest_segment = seg;
        }
    }
    
    return closest_segment;
}

int SegmentManager3D::getNumSegments() const {
    return static_cast<int>(segments_.size());
}

int SegmentManager3D::getNumSegmentsInLayer(const std::string& layer) const {
    auto it = layer_segments_.find(layer);
    if (it != layer_segments_.end()) {
        return static_cast<int>(it->second.size());
    }
    return 0;
}

void SegmentManager3D::printStatistics() const {
    std::cout << "\n=== SegmentManager3D Statistics ===" << std::endl;
    std::cout << "Total segments: " << getNumSegments() << std::endl;
    
    for (const auto& layer_pair : layer_segments_) {
        std::string layer = layer_pair.first;
        int num_segments = static_cast<int>(layer_pair.second.size());
        
        std::cout << "  Layer " << layer << ": " << num_segments << " segments" << std::endl;
        
        // 计算利用率
        double total_utilization = 0.0;
        int segments_with_cells = 0;
        
        for (Segment3D* seg : layer_pair.second) {
            double util_ratio = seg->getUtilizationRatio();
            total_utilization += util_ratio;
            if (seg->getCells().size() > 0) {
                segments_with_cells++;
            }
        }
        
        if (num_segments > 0) {
            double avg_utilization = total_utilization / num_segments;
            std::cout << "    Average utilization: " << (avg_utilization * 100) << "%" << std::endl;
            std::cout << "    Segments with cells: " << segments_with_cells 
                      << " / " << num_segments << std::endl;
        }
    }
    
    std::cout << "Cells with segment assignments: " << cell_to_segments_.size() << std::endl;
    std::cout << "==================================\n" << std::endl;
}

void SegmentManager3D::validateSegments() const {
    // std::cout << "Validating segments..." << std::endl;
    
    int errors = 0;
    
    // 检查所有段
    for (const auto& segment : segments_) {
        // 检查段宽度
        // if (segment->getWidth().v <= 0) {
        //     std::cerr << "Error: Segment " << segment->getId() 
        //               << " has invalid width: " << segment->getWidth().v << std::endl;
        //     errors++;
        // }
        
        // 检查利用率不超过宽度
        if (segment->getUtilization() > segment->getWidth().v) {
            std::cerr << "Error: Segment " << segment->getId() 
                      << " utilization (" << segment->getUtilization()
                      << ") exceeds width (" << segment->getWidth().v << ")" << std::endl;
            errors++;
        }
        
        // 检查段中的单元
        const auto& cells = segment->getCells();
        for (size_t i = 1; i < cells.size(); i++) {
            Node3D* prev_cell = cells[i-1];
            Node3D* curr_cell = cells[i];
            
            // 检查单元是否按X坐标排序
            if (prev_cell->getLocation().x > curr_cell->getLocation().x) {
                std::cerr << "Error: Cells not sorted in segment " << segment->getId() << std::endl;
                errors++;
            }
            
            // 检查单元是否重叠（简化检查）
            DbuX prev_right = prev_cell->getLocation().x + prev_cell->getWidth();
            DbuX curr_left = curr_cell->getLocation().x;
            
            if (prev_right > curr_left) {
                std::cerr << "Error: Overlap detected in segment " << segment->getId()
                          << " between cells " << prev_cell->getId() 
                          << " and " << curr_cell->getId() << std::endl;
                errors++;
            }
        }
    }
    
    // 检查反向映射的一致性
    for (const auto& cell_pair : cell_to_segments_) {
        int cell_id = cell_pair.first;
        const auto& seg_list = cell_pair.second;
        
        for (Segment3D* seg : seg_list) {
            // 检查单元是否确实在段中
            const auto& cells = seg->getCells();
            bool found = false;
            for (Node3D* cell : cells) {
                if (cell->getId() == cell_id) {
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                std::cerr << "Error: Cell " << cell_id 
                          << " not found in segment " << seg->getId()
                          << " (inconsistent reverse mapping)" << std::endl;
                errors++;
            }
        }
    }
    
    if (errors != 0) 
        std::cout << "Segment validation failed with " << errors << " errors." << std::endl;
        // std::cout << "Segment validation passed!" << std::endl;
    // } else {
    //     std::cout << "Segment validation failed with " << errors << " errors." << std::endl;
    // }
}

// 辅助函数实现

std::vector<SegmentManager3D::Blockage3D> 
SegmentManager3D::getBlockagesForRow(const std::string& layer, int row_id) const {
    std::vector<Blockage3D> blockages;
    
    // 简化版本：从Architecture3D获取阻塞信息
    // 实际实现需要从Architecture3D获取3D阻塞数据
    
    // 示例：添加一些测试阻塞
    if (layer == "bottom" && row_id == 0) {
        // 在底层第0行添加一个测试阻塞
        blockages.push_back({DbuX{100}, DbuX{200}, DbuY{0}, DbuY{10}, "bottom", "test"});
    }
    
    return blockages;
}

// 辅助函数：获取层的Z坐标
// static double getLayerZCoordinate(const std::string& layer) {
//     if (layer == "bottom") return 0.0;
//     if (layer == "upper") return 1.0;
//     return 0.0;  // 默认
// }

} // namespace tdp
