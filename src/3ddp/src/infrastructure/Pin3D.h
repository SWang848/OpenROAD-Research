// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <string>

#include "Common3D.h"
#include "Coordinates.h"
#include "odb/db.h"

// Forward declarations
namespace tdp {
    class Node3D;
    class Net3D;
}

namespace tdp {

/**
 * @class Pin3D
 * @brief 3D pin representation similar to DPL's Pin class
 * 
 * Represents a pin in 3D space with layer information.
 */
class Pin3D
{
public:
    enum Direction
    {
        Dir_IN,
        Dir_OUT,
        Dir_INOUT,
        Dir_UNKNOWN
    };

    Pin3D();
    ~Pin3D();
    
    // 方向管理
    void setDirection(int dir);
    int getDirection() const;
    
    // 连接管理
    void setNode(Node3D* node);
    void setNet(Net3D* net);
    Node3D* getNode() const;
    Net3D* getNet() const;
    
    // 偏移管理
    void setOffsetX(DbuX offsetX);
    DbuX getOffsetX() const;
    void setOffsetY(DbuY offsetY);
    DbuY getOffsetY() const;
    void setOffsetZ(DbuZ offsetZ);
    DbuZ getOffsetZ() const;
    
    // 层管理
    void setPinLayer(int layer);
    int getPinLayer() const;
    
    // 尺寸管理
    void setPinWidth(DbuX width);
    DbuX getPinWidth() const;
    void setPinHeight(DbuY height);
    DbuY getPinHeight() const;
    
    // 3D位置计算
    Pointdp3D getLocation3D() const;
    
    // 数据库连接
    void setDbITerm(odb::dbITerm* iterm);
    void setDbBTerm(odb::dbBTerm* bterm);
    odb::dbITerm* getDbITerm() const;
    odb::dbBTerm* getDbBTerm() const;
    
    // 名称
    std::string getName() const;
    void setName(const std::string& name);

    //权重管理
    void setEarlyWeight(double weight);
    double getEarlyWeight() const;
    void setLateWeight(double weight);
    double getLateWeight() const;
    void setEarlyWeightini(double weight);
    double getEarlyWeightini() const;
    void setLateWeightini(double weight);
    double getLateWeightini() const;


    // 时序信息
    void setEarlySlack(double slack);
    double getEarlySlack() const;
    void setLateSlack(double slack);
    double getLateSlack() const;
    void setEarlyArrival(double arrival);
    double getEarlyArrival() const;
    void setLateArrival(double arrival);
    double getLateArrival() const;
    void setEarlyRequired(double required);
    double getEarlyRequired() const;
    void setLateRequired(double required);
    double getLateRequired() const;

    // Arc delay 缓存（库级常数，不随位置变化）
    void setCachedArcDelay(double early, double late) {
        cached_arc_early_ = early;
        cached_arc_late_ = late;
        arc_cached_ = true;
    }
    double getCachedArcEarly() const { return cached_arc_early_; }
    double getCachedArcLate()  const { return cached_arc_late_; }
    bool   isArcCached()       const { return arc_cached_; }

private:
    // 引脚尺寸
    DbuX pinWidth_{0};
    DbuY pinHeight_{0};
    
    // 方向
    int dir_{Dir_INOUT};
    
    // 层
    int pinLayer_{0};
    
    // 连接
    Node3D* node_{nullptr};
    Net3D* net_{nullptr};
    
    // 偏移（相对于单元中心）
    DbuX offsetX_{0};
    DbuY offsetY_{0};
    DbuZ offsetZ_{0};
    
    // 数据库对象
    odb::dbITerm* db_iterm_{nullptr};
    odb::dbBTerm* db_bterm_{nullptr};
    
    // 名称
    std::string name_;
    
    // 禁止拷贝
    Pin3D(const Pin3D&) = delete;
    Pin3D& operator=(const Pin3D&) = delete;


    double early_weight_{1.0};  // 默认值1.0
    double early_weight_ini_{1.0};
    double late_weight_{1.0};
    double late_weight_ini_{1.0};   // 默认值1.0
    
    // 时序信息
    double early_slack_{0.0};
    double late_slack_{0.0};
    double early_arrival_{0.0};
    double late_arrival_{0.0};
    double early_required_{0.0};
    double late_required_{0.0};

    // Arc delay 缓存
    double cached_arc_early_{0.0};
    double cached_arc_late_{0.0};
    bool   arc_cached_{false};
};

} // namespace tdp