// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

//#include "Coordinates.h"
#include "Common3D.h"

// 包含OpenDB头文件（仿照DPL）
#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"

namespace tdp {

// Forward declarations
class Node3D;
class Padding3D;
class Network3D;

// Group3D class following DPL's Group design pattern
class Group3D {
public:
    // getters (following DPL interface)
    std::string getName() const { return name_; }
    const std::vector<odb::Rect>& getRects() const { return region_boundaries_; }
    std::vector<Node3D*> getCells() const { return cells_; }
    const odb::Rect& getBBox() const { return boundary_; }
    double getUtil() const { return util_; }
    int getId() const { return id_; }
    
    // setters (following DPL interface)
    void setId(int id) { id_ = id; }
    void setName(const std::string& in) { name_ = in; }
    void addRect(const odb::Rect& in) { region_boundaries_.push_back(in); }
    void addCell(Node3D* cell) { cells_.push_back(cell); }
    void setBoundary(const odb::Rect& in) { boundary_ = in; }
    void setUtil(double in) { util_ = in; }
    
    // Additional 3D-specific methods
    void setGroup(Node3D* cell) {
        // Simple implementation for compatibility
    }

private:
    int id_{0};
    std::string name_;
    std::vector<odb::Rect> region_boundaries_;
    std::vector<Node3D*> cells_;
    odb::Rect boundary_;
    double util_{0.0};
};

/**
 * @class Architecture3D
 * @brief 3D layout architecture information for multi-layer 3D IC designs
 * 
 * This class extends the traditional 2D Architecture to support 3D stacking,
 * including layer-specific rows, sites, and inter-layer alignment constraints.
 * It is designed to work with TaiWei's 3D PDK for ASAP7-NanGate45 face-to-face stacking.
 */
class Architecture3D
{
public:
    // Nested Row class for 3D
    class Row3D;

    Architecture3D();
    ~Architecture3D();

    // Layer management
    enum LayerType {
        LAYER_BOTTOM = 0,    // Lower layer (ASAP7 in TaiWei PDK)
        LAYER_UPPER = 1,     // Upper layer (NanGate45 in TaiWei PDK)
        LAYER_COUNT = 2
    };

    static const char* layerToString(LayerType layer);
    static LayerType stringToLayer(const std::string& layer_name);

    // 统一的层归属判定：依据标准单元（master）名后缀 "upper"/"bottom" 判定
    // 该 master 属于上层还是下层。与具体工艺、同构/异构无关。
    static LayerType layerOfMaster(odb::dbMaster* master);

    // 按 BTerm（IO terminal）pin 几何所在金属层相对 hb_layer 判定其所属层：
    // pin 金属层在 hb_layer 之下 -> 下层，之上 -> 上层。
    static LayerType layerOfBTerm(odb::dbBTerm* bterm, odb::dbTech* tech);

    // Layer operations
    bool addLayer(LayerType layer_type, const std::string& layer_name);
    bool removeLayer(LayerType layer_type);
    bool hasLayer(LayerType layer_type) const;
    const std::string& getLayerName(LayerType layer_type) const;
    
    // Row management per layer (类似DPL的Row管理)
    const std::vector<Row3D*>& getRows(LayerType layer) const;
    int getNumRows(LayerType layer) const;
    Row3D* getRow(LayerType layer, int row_idx) const;
    Row3D* createAndAddRow(LayerType layer);
    
    // 从数据库创建Row（类似DPL的createArchitecture）
    Row3D* createRowFromDatabase(LayerType layer, odb::dbRow* db_row);
    
    // Region management per layer
    const std::vector<Group3D*>& getRegions(LayerType layer) const;
    int getNumRegions(LayerType layer) const;
    Group3D* getRegion(LayerType layer, int region_idx) const;
    void* createAndAddRegion(LayerType layer);
    
    // Boundary management
    DbuX getMinX(LayerType layer) const;
    DbuX getMaxX(LayerType layer) const;
    DbuY getMinY(LayerType layer) const;
    DbuY getMaxY(LayerType layer) const;
    
    DbuX getWidth(LayerType layer) const { return getMaxX(layer) - getMinX(layer); }
    DbuY getHeight(LayerType layer) const { return getMaxY(layer) - getMinY(layer); }
    
    void setMinX(LayerType layer, DbuX xmin);
    void setMaxX(LayerType layer, DbuX xmax);
    void setMinY(LayerType layer, DbuY ymin);
    void setMaxY(LayerType layer, DbuY ymax);
    
    // 3D-specific operations
    DBU getLayerZOffset(LayerType layer) const;
    void setLayerZOffset(LayerType layer, DBU z_offset);
    
    DBU getInterLayerSpacing() const { return inter_layer_spacing_; }
    void setInterLayerSpacing(DBU spacing) { inter_layer_spacing_ = spacing; }
    
    DBU getHBTPitch() const { return hbt_pitch_; }
    void setHBTPitch(DBU pitch) { hbt_pitch_ = pitch; }
    
    // Typed coordinate helpers
    static DbuX makeDbuX(DBU value) { return DbuX(value); }
    static DbuY makeDbuY(DBU value) { return DbuY(value); }
    static DbuZ makeDbuZ(DBU value) { return DbuZ(value); }
    
    // Helper to convert DBU to typed coordinates
    static DbuX dbuToDbuX(DBU value) { return DbuX(value); }
    static DbuY dbuToDbuY(DBU value) { return DbuY(value); }
    static DbuZ dbuToDbuZ(DBU value) { return DbuZ(value); }
    
    // Cell compatibility checks
    bool isSingleHeightCell(const Node3D* node) const;
    bool isMultiHeightCell(const Node3D* node) const;
    int getCellHeightInRows(const Node3D* node) const;
    
    // Row height and cell placement checks
    DBU getRowHeight() const;
    bool isCellInRows(const Node3D* cell) const;
    
    // Cell layer determination
    LayerType getCellLayer(odb::dbInst* inst) const;
    
    // Power compatibility with 3D considerations
    bool powerCompatible(const Node3D* node, const Row3D* row, bool& flip) const;
    
    // Inter-layer alignment checks
    bool checkInterLayerAlignment(const Node3D* bottom_node, const Node3D* upper_node) const;
    bool findNearestHBTLocation(const Node3D* bottom_node, Point& hbt_location) const;
    
    // Padding management
    void setPadding(Padding3D* padding);
    void setSiteWidth(LayerType layer, const DbuX& site_width);
    void setUsePadding(bool use_padding = true) { use_padding_ = use_padding; }
    bool getUsePadding() const { return use_padding_; }
    
    bool getCellPadding(const Node3D* node, int& left_padding, int& right_padding) const;
    bool getCellPadding(const Node3D* node, DbuX& left_padding, DbuX& right_padding) const;
    void flipCellPadding(const Node3D* node);
    int getCellSpacing(const Node3D* left_node, const Node3D* right_node) const;
    
    // Clear all data
    void clear();
    void clearLayer(LayerType layer);
    
    // Post-processing after loading
    void postProcess(Network3D* network);
    
    // Find closest row for a given Y coordinate
    int findClosestRow(LayerType layer, DbuY y) const;
    
    // Debug and statistics
    void printLayerStats(LayerType layer) const;
    void printAllStats() const;

private:
    // Layer information
    struct LayerInfo {
        std::string name;
        DBU z_offset = 0;
        DbuX xmin = DbuX::max();
        DbuX xmax = DbuX::lowest();
        DbuY ymin = DbuY::max();
        DbuY ymax = DbuY::lowest();
        DBU site_width = 0;
        std::vector<Row3D*> rows;
        std::vector<Group3D*> regions;
    };
    
    // Layer data
    std::map<LayerType, LayerInfo> layers_;
    
    // 3D-specific parameters
    DBU inter_layer_spacing_ = 0;    // Spacing between layers
    DBU hbt_pitch_ = 0;              // Hybrid Bonding Terminal pitch
    
    // Padding
    bool use_padding_ = false;
    Padding3D* padding_ = nullptr;
    
    // Helper methods
    LayerInfo& getLayerInfo(LayerType layer);
    const LayerInfo& getLayerInfo(LayerType layer) const;
    bool validateLayer(LayerType layer) const;
};

/**
 * @class Architecture3D::Row3D
 * @brief 3D row information with layer identification
 */
class Architecture3D::Row3D
{
public:
    enum PowerType {
        POWER_UNKNOWN = 0,
        POWER_VDD = 1,
        POWER_VSS = 2
    };

    // Constructor
    Row3D(Architecture3D::LayerType layer);
    
    // Getters and setters
    void setId(int id) { id_ = id; }
    int getId() const { return id_; }
    
    LayerType getLayer() const { return layer_; }
    
    void setOrient(const odb::dbOrientType& orient) { site_orient_ = orient; }
    odb::dbOrientType getOrient() const { return site_orient_; }
    
    void setSymmetry(unsigned sym) { site_symmetry_ = sym; }
    unsigned getSymmetry() const { return site_symmetry_; }
    
    void setBottom(DbuY bottom) { row_location_ = bottom; }
    void setHeight(DbuY height) { row_height_ = height; }
    void setSubRowOrigin(DbuX origin) { subrow_origin_ = origin; }
    void setSiteSpacing(DbuX spacing) { site_spacing_ = spacing; }
    void setSiteWidth(DbuX width) { site_width_ = width; }
    void setNumSites(int num_sites) { num_sites_ = num_sites; }
    
    DbuY getBottom() const { return row_location_; }
    DbuY getTop() const { return row_location_ + row_height_; }
    DbuX getLeft() const { return subrow_origin_; }
    DbuX getRight() const { return subrow_origin_ + num_sites_ * site_spacing_.v; }
    DbuY getHeight() const { return row_height_; }
    DbuX getSiteWidth() const { return site_width_; }
    DbuX getSiteSpacing() const { return site_spacing_; }
    int getNumSites() const { return num_sites_; }
    
    void setTopPower(PowerType power) { top_power_ = power; }
    PowerType getTopPower() const { return top_power_; }
    
    void setBottomPower(PowerType power) { bottom_power_ = power; }
    PowerType getBottomPower() const { return bottom_power_; }
    
    // 3D-specific methods
    DBU getZOffset() const;
    Pointdp3D get3DLocation(DbuX x, DbuY y) const;
    
    // Debug
    void printInfo() const;

private:
    int id_ = -1;                       // Row ID
    LayerType layer_;                   // Which layer this row belongs to
    DbuY row_location_ = 0;             // Y-location of the row
    DbuY row_height_ = 0;               // Height of the row
    DbuX subrow_origin_ = 0;            // Starting X location of the row
    DbuX site_spacing_ = 0;             // Spacing between sites
    DbuX site_width_ = 0;               // Width of sites
    int num_sites_ = 0;                 // Number of sites
    odb::dbOrientType site_orient_;     // Orientation of sites
    unsigned site_symmetry_ = 0;        // Symmetry of sites
    PowerType top_power_ = POWER_UNKNOWN;
    PowerType bottom_power_ = POWER_UNKNOWN;
};

} // namespace tdp