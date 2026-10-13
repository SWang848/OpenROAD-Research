// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once
#include "utl/Logger.h"
#include "tdp/Timingdp.h"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <optional>
#include <map>

// 前向声明
namespace sta {
    class dbSta;
    class MinMax;
    class Scene;
    class Path;
    class PathEnd;
    class Parasitic;
    class dbNetwork;
    class Pin;
    class Net;
    class Instance;
    class LibertyLibrary;
    class RiseFall;  
    class Vertex;     
    class Graph;    
}

namespace odb {
    class dbNet;
}

namespace tdp {

// ============================================================================
// 3D时序引擎（合并 TimingAnalyzer3D 和 StaInterface3D 的功能）
// ============================================================================

class TimingEngine3D {
public:
    TimingEngine3D(sta::dbSta* sta, Timingdp* timingdp);
    ~TimingEngine3D();
    
    // 禁止拷贝
    TimingEngine3D(const TimingEngine3D&) = delete;
    TimingEngine3D& operator=(const TimingEngine3D&) = delete;
    
    // ========================================================================
    // 初始化接口（来自 TimingAnalyzer3D）
    // ========================================================================
    
    // 初始化STA
    bool initialize();
    
    // ========================================================================
    // 完整时序分析流程（来自 TimingAnalyzer3D）
    // ========================================================================
    
    // 运行完整的时序分析
    void runFullAnalysis();
    
    // 运行增量时序分析
    void runIncrementalAnalysis();
    
    // ========================================================================
    // 时序驱动布局接口（新增，参考 rsz）
    // ========================================================================
    
    // 时序分析前准备（类似 rsz::resizeSlackPreamble）
    void resizeSlackPreamble();
    
    // 获取网络时序松弛度（类似 rsz::resizeNetSlack）
    std::optional<sta::Slack> getResizeNetSlack(sta::Net* net);
    std::optional<sta::Slack> getResizeNetSlack(odb::dbNet* db_net);
    
    // ========================================================================
    // 时序结果获取（来自 TimingAnalyzer3D）
    // ========================================================================
    
    // 获取最差slack
    sta::Slack getWorstSlack(const sta::MinMax* min_max = nullptr);
    
    // 获取网络slack
    sta::Slack getNetSlack(const std::string& net_name,
                          const sta::MinMax* min_max = nullptr);
    
    // 获取引脚slack
    sta::Slack getPinSlack(const std::string& pin_name,
                          const sta::MinMax* min_max = nullptr);
    
    
    // 获取违反约束的路径
    std::vector<sta::Path*> getViolatingPaths(const sta::MinMax* min_max = nullptr);
    
    // 获取时序统计信息
    struct TimingStats {
        double early_wns = 0.0;
        double early_tns = 0.0;
        double late_wns = 0.0;
        double late_tns = 0.0;
    };
    
    TimingStats getTimingStats();

    struct PinTimingInfo {
        double early_arrival;    // 早期到达时间
        double late_arrival;     // 晚期到达时间
        double early_required;   // 早期要求时间
        double late_required;    // 晚期要求时间
        double early_slack;      // 早期松弛度
        double late_slack;       // 晚期松弛度
        bool is_valid;           // 数据是否有效
    };

     // 获取引脚的完整时序信息
    PinTimingInfo getPinTimingInfo(const std::string& pin_name);
    PinTimingInfo getPinTimingInfo(sta::Pin* pin);
    PinTimingInfo getPinTimingInfoDirect(sta::Pin* pin);
    
    // 获取引脚到达时间
    double getPinArrivalTime(const std::string& pin_name, 
                            bool is_early = false);
    double getPinArrivalTime(sta::Pin* pin, 
                            bool is_early = false);
    
    // 获取引脚要求时间
    double getPinRequiredTime(const std::string& pin_name,
                             bool is_early = false);
    double getPinRequiredTime(sta::Pin* pin,
                             bool is_early = false);
    
    // 获取引脚松弛度
    double getPinSlack(const std::string& pin_name,
                      bool is_early = false);
    double getPinSlack(sta::Pin* pin,
                      bool is_early = false);
    

    // ========================================================================
    // 报告生成（简化版本）
    // ========================================================================
    
    // 生成完整时序报告
    void reportTiming(const std::string& output_file,
                     int max_paths = 10,
                     const sta::MinMax* min_max = nullptr);
    
    
    // ========================================================================
    // 工具函数
    // ========================================================================
    
    // 获取底层STA对象
    sta::dbSta* getSta() const { return sta_; }
    
    // 获取底层网络对象
    sta::dbNetwork* getNetwork() const;
    
    // 获取默认scene
    sta::Scene* getDefaultScene() const;

    sta::Vertex* getStaVertex(class Node3D* cell) const;
    
    // 获取默认MinMax
    const sta::MinMax* getDefaultMinMax() const;
    
    // 获取Timingdp
    Timingdp* getPlacerDB() const { return timingdp_; }
    

   // cell 内弧延迟
    struct ArcDelayInfo {
        bool is_valid;
        double early_delay;  // 早时序弧延迟
        double late_delay;   // 晚时序弧延迟
        std::string from_pin_name;
        std::string to_pin_name;
    };

    ArcDelayInfo getArcDelay(const std::string& from_pin_name, 
                            const std::string& to_pin_name);
    ArcDelayInfo getArcDelay(sta::Pin* from_pin, 
                            sta::Pin* to_pin);

    
private:
    sta::dbSta* sta_ = nullptr;
    Timingdp* timingdp_ = nullptr;
    utl::Logger* logger_;

    
    // 时序驱动布局状态（类似 rsz）
    std::map<const sta::Net*, sta::Slack> net_slack_map_;
    float worst_slack_nets_percent_ = 10.0f;  // 最差时序网络百分比
    
    // 设计统计结构
    struct DesignStats {
        int num_instances = 0;
        int num_nets = 0;
        int num_pins = 0;
        int num_corners = 0;
        double design_area = 0.0;
    };
    
    // 内部辅助函数
    bool setupAnalysis();
    
    // 设计统计
    DesignStats getDesignStats() const;
    
    // 图构建和层级化
    void updateGraph();
    void ensureGraphBuilt();
    void ensureLevelized();
    
    // 结果处理
    void processTimingResults();
    
    // 引脚查找
    sta::Pin* findPinByName(const std::string& pin_name) const;
    
    // 网络查找
    sta::Net* findNetByName(const std::string& net_name) const;
    
    // 调试输出
    void debugPrintSetup() const;
    
    // 私有方法：从odb::dbNet获取sta::Net
    sta::Net* getStaNetFromDbNet(odb::dbNet* db_net) const;

    sta::Vertex* getVertexFromPin(sta::Pin* pin) const;

    const sta::MinMax* getMinMax(bool is_early) const;
    const sta::RiseFall* getDefaultRiseFall() const;
};
} // namespace tdp