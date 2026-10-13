// include/SegmentManager3D.h
#pragma once

#include "Segment3D.h"
#include "Architecture3D.h"
#include "Node3D.h"
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace tdp {

/**
 * @brief 3D段管理器
 * 
 * 管理所有层的段，提供段查询和操作接口。
 * 类似DPL的DetailedMgr，但支持3D分层。
 */
class SegmentManager3D {
public:
    SegmentManager3D() : arch_(nullptr), next_segment_id_(0) {}
    SegmentManager3D(Architecture3D* arch) : arch_(arch), next_segment_id_(0) {}
    ~SegmentManager3D();
    
    // 初始化
    void initializeSegments();
    void clearSegments();
    
    // 段查询
    std::vector<Segment3D*> getCellSegments(const Node3D* cell) const;
    std::vector<Segment3D*> getSegmentsAtPosition(const Pointdp3D& pos,
                                                 const std::string& layer) const;
    std::vector<Segment3D*> getSegmentsInRow(const std::string& layer, 
                                            int row_id) const;
    
    // 段操作
    bool addCellToSegment(Node3D* cell, int seg_id);
    bool removeCellFromSegment(const Node3D* cell, int seg_id);
    bool moveCellToSegment(Node3D* cell, int new_seg_id);
    
    // 段查找
    Segment3D* findSegmentById(int seg_id) const;
    Segment3D* findClosestSegment(const Node3D* cell) const;
    
    // 统计信息
    int getNumSegments() const;
    int getNumSegmentsInLayer(const std::string& layer) const;
    void printStatistics() const;
    
    // 调试
    void validateSegments() const;
    
private:
    // 段创建
    void createSegmentsForLayer(const std::string& layer);
    void createSegmentsForRow(const std::string& layer, int row_id);
    void createSingleSegment(const std::string& layer, int row_id, 
                             DbuX min_x, DbuX max_x,
                             DbuY row_ymin = DbuY{0}, DbuY row_ymax = DbuY{0});
    
    // 阻塞检测（3D版本）
    struct Blockage3D {
        DbuX x_min, x_max;
        DbuY y_min, y_max;
        std::string layer;
        std::string type;
    };
    
    std::vector<Blockage3D> getBlockagesForRow(const std::string& layer, 
                                              int row_id) const;
    
    // 数据结构
    Architecture3D* arch_;
    
    // 所有段的列表
    std::vector<std::unique_ptr<Segment3D>> segments_;
    
    // 分层索引：layer -> 段列表
    std::map<std::string, std::vector<Segment3D*>> layer_segments_;
    
    // 行索引：layer -> row_id -> 段列表
    std::map<std::string, std::map<int, std::vector<Segment3D*>>> segments_in_row_;
    
    // 反向映射：cell_id -> 段列表
    std::map<int, std::vector<Segment3D*>> cell_to_segments_;
    
    // ID计数器
    int next_segment_id_;
};

} // namespace tdp
