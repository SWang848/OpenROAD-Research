// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <vector>
#include <string>
#include <memory>

#include "Common3D.h"
#include "Coordinates.h"
#include "odb/db.h"
#include "odb/geom.h"

// Forward declarations
namespace odb {
    class dbMaster;
    class dbTechLayer;
}

namespace tdp {

// Forward declarations
class Architecture3D;
class Grid3D;
class PlacementDRC;

/**
 * @class MasterEdge3D
 * @brief 3D master edge representation
 * 
 * Represents edges of a master cell in 3D space.
 */
class MasterEdge3D
{
public:
    MasterEdge3D(unsigned int type, const odb::Rect& box);
    
    unsigned int getEdgeType() const;
    const odb::Rect& getBBox() const;
    
private:
    unsigned int edge_type_idx_{0};
    odb::Rect bbox_;
};

/**
 * @class Master3D
 * @brief 3D master representation similar to DPL's Master class
 * 
 * Represents a master cell with 3D-specific properties.
 */
class Master3D
{
public:
    Master3D();
    ~Master3D();
    
    // 基本信息
    bool isMultiRow() const;
    bool isMultiHeight() const;
    const std::vector<MasterEdge3D>& getEdges() const;
    odb::Rect getBBox() const;
    
    // 电源信息
    int getBottomPowerType() const;
    int getTopPowerType() const;
    
    // 设置方法
    void setMultiRow(bool multi_row);
    void setMultiHeight(bool multi_height);
    void addEdge(const MasterEdge3D& edge);
    void clearEdges();
    void setBBox(const odb::Rect& box);
    void setBottomPowerType(int bottom_pwr);
    void setTopPowerType(int top_pwr);
    void setDbMaster(odb::dbMaster* db_master);
    
    // 获取方法
    odb::dbMaster* getDbMaster() const;
    std::string getName() const;
    DbuX getWidth() const;
    DbuY getHeight() const;
    
    // 3D特定方法
    bool canBePlacedOnLayer(const std::string& layer) const;
    bool requiresHBTConnection() const;
    
private:
    odb::dbMaster* db_master_{nullptr};
    odb::Rect boundary_box_;
    bool is_multi_row_{false};
    bool is_multi_height_{false};
    std::vector<MasterEdge3D> edges_;
    int bottom_pwr_{0};
    int top_pwr_{0};
    
    // 禁止拷贝
    Master3D(const Master3D&) = delete;
    Master3D& operator=(const Master3D&) = delete;
};

} // namespace tdp