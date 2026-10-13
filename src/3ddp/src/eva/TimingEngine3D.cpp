// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "TimingEngine3D.h"
#include "tdp/Timingdp.h"

#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"

#include "sta/Sdc.hh"
#include "sta/PortDirection.hh"
#include "sta/Graph.hh"
#include "sta/Search.hh"
#include "sta/Report.hh"
#include "sta/Units.hh"
#include "sta/PathEnd.hh"
#include "sta/PatternMatch.hh"
#include "sta/MinMax.hh"
#include "sta/Liberty.hh"
#include "sta/Transition.hh"  
#include "sta/Sta.hh"          
#include "sta/Network.hh"   
#include "utl/Logger.h"            

#include "../infrastructure/Node3D.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <tcl.h>

using utl::TDP;

namespace tdp {

// ============================================================================
// TimingEngine3D 实现
// ============================================================================

TimingEngine3D::TimingEngine3D(sta::dbSta* sta, Timingdp* timingdp)
    : sta_(sta), timingdp_(timingdp), 
      logger_(timingdp ? timingdp->getLogger() : nullptr)
{
    if (!sta_) {
        throw std::runtime_error("TimingEngine3D: STA pointer is null");
    }
    if (!timingdp_) {
        throw std::runtime_error("TimingEngine3D: Timingdp pointer is null");
    }
}

TimingEngine3D::~TimingEngine3D()
{
    // 清理资源
}

// ============================================================================
// 初始化接口
// ============================================================================

bool TimingEngine3D::initialize()
{
    if (!sta_) {
        std::cerr << "Error: STA not set" << std::endl;
        return false;
    }
    
    if (!timingdp_) {
        std::cerr << "Error: Timingdp not set" << std::endl;
        return false;
    }
    
    // 设置分析环境
    if (!setupAnalysis()) {
        return false;
    }
    
    std::cout << "TimingEngine3D initialized successfully" << std::endl;
    debugPrintSetup();
    
    return true;
}

// ============================================================================
// 完整时序分析流程
// ============================================================================

void TimingEngine3D::runFullAnalysis()
{

    // 确保图已构建
    ensureGraphBuilt();

    // 验证寄生参数注入


    // 更新时序
    sta_->updateTiming(true);
    sta_->findRequireds();
    
    // // 处理结果
    // processTimingResults();
}

void TimingEngine3D::runIncrementalAnalysis()
{
    // 增量更新时序
    sta_->updateTiming(true);  // incremental update
    
    // 处理结果
    processTimingResults();
}

// ============================================================================
// 时序驱动布局接口（新增，参考 rsz）
// ============================================================================

void TimingEngine3D::resizeSlackPreamble()
{

    // 清除现有的时序松弛度映射
    net_slack_map_.clear();
    
    std::cout << "TimingEngine3D: resize slack preamble completed" << std::endl;
}

std::optional<sta::Slack> TimingEngine3D::getResizeNetSlack(sta::Net* net)
{
    if (!net) {
        return std::nullopt;
    }
    
    auto it = net_slack_map_.find(net);
    if (it != net_slack_map_.end()) {
        return it->second;
    }
    
    // 如果映射中没有，直接计算
    sta::Slack slack = sta_->slack(net, getDefaultMinMax());
    net_slack_map_[net] = slack;
    return slack;
}

std::optional<sta::Slack> TimingEngine3D::getResizeNetSlack(odb::dbNet* db_net)
{
    if (!db_net) {
        return std::nullopt;
    }
    
    // 将odb::dbNet转换为sta::Net
    sta::Net* net = getStaNetFromDbNet(db_net);
    if (!net) {
        return std::nullopt;
    }
    
    return getResizeNetSlack(net);
}

// ============================================================================
// 时序结果获取
// ============================================================================

sta::Slack TimingEngine3D::getWorstSlack(const sta::MinMax* min_max)
{
    
    const sta::MinMax* mm = min_max ? min_max : getDefaultMinMax();
    return sta_->worstSlack(mm);
}

sta::Slack TimingEngine3D::getNetSlack(const std::string& net_name,
                                      const sta::MinMax* min_max)
{
    
    sta::dbNetwork* network = getNetwork();
    if (!network) {
        throw std::runtime_error("Network not available");
    }
    
    // 查找网络
    sta::Net* net = network->findNet(net_name.c_str());
    if (!net) {
        throw std::runtime_error("Net '" + net_name + "' not found");
    }
    
    const sta::MinMax* mm = min_max ? min_max : getDefaultMinMax();
    return sta_->slack(net, mm);
}

sta::Slack TimingEngine3D::getPinSlack(const std::string& pin_name,
                                      const sta::MinMax* min_max)
{
    
    sta::Pin* pin = findPinByName(pin_name);
    if (!pin) {
        throw std::runtime_error("Pin '" + pin_name + "' not found");
    }
    
    const sta::MinMax* mm = min_max ? min_max : getDefaultMinMax();
    return sta_->slack(pin, sta::RiseFallBoth::riseFall(), sta_->scenes(), mm);
}

std::vector<sta::Path*> TimingEngine3D::getViolatingPaths(const sta::MinMax* min_max)
{
    
    const sta::MinMax* mm = min_max ? min_max : getDefaultMinMax();
    
    // 获取违反约束的路径（slack < 0）- 使用正确的OpenSTA API
    // 注意：需要将MinMax转换为MinMaxAll
    const sta::MinMaxAll* min_max_all = (mm == sta::MinMax::max()) ? 
        sta::MinMaxAll::max() : sta::MinMaxAll::min();
    
    // 调用findPathEnds获取违反约束的路径端点（slack < 0）
    sta::StringSeq group_names;
    sta::PathEndSeq path_ends = sta_->findPathEnds(
        nullptr,  // from
        nullptr,  // thrus
        nullptr,  // to
        false,    // unconstrained
        sta_->scenes(), // scenes
        min_max_all, // min_max_all
        2000,     // group_path_count
        1,     // endpoint_path_count
        true,     // unique_pins
        false,    // unique_edges (简化)
        -sta::INF, // slack_min
        0.0,      // slack_max (只获取slack < 0的路径)
        true,     // sort_by_slack
        group_names, // group_names
        true,     // setup
        false,    // hold
        false,    // recovery
        false,    // removal
        false,    // clk_gating_setup
        false     // clk_gating_hold
    );
    
    std::vector<sta::Path*> paths;
    for (sta::PathEnd* path_end : path_ends) {
        if (path_end) {
            sta::Path* path = path_end->path();
            if (path) {
                paths.push_back(path);
            }
        }
    }
    
    return paths;
}

TimingEngine3D::TimingStats TimingEngine3D::getTimingStats()
{
    TimingStats stats;
    stats.early_wns = sta_->worstSlack(sta::MinMax::min());
    stats.early_tns = sta_->totalNegativeSlack(sta::MinMax::min());
    stats.late_wns  = sta_->worstSlack(sta::MinMax::max());
    stats.late_tns  = sta_->totalNegativeSlack(sta::MinMax::max());
    return stats;
}

// ============================================================================
// 报告生成
// ============================================================================

void TimingEngine3D::reportTiming(const std::string& output_file,
                                 int max_paths,
                                 const sta::MinMax* min_max)
{
    
    const sta::MinMax* mm = min_max ? min_max : getDefaultMinMax();
    
    // 打开输出文件
    std::ofstream out(output_file);
    if (!out.is_open()) {
        throw std::runtime_error("TimingEngine3D: Failed to open output file: " + output_file);
    }
    
    // 获取Tcl解释器
    Tcl_Interp* interp = sta_->tclInterp();
    if (!interp) {
        out << "Error: Tcl interpreter not available" << std::endl;
        out.close();
        return;
    }
    
    // 使用Tcl命令生成时序报告
    std::string cmd = "report_checks -path_delay max -max_paths " + std::to_string(max_paths);
    if (mm == sta::MinMax::min()) {
        cmd = "report_checks -path_delay min -max_paths " + std::to_string(max_paths);
    }
    
    // 直接执行命令并捕获输出
    int result = Tcl_Eval(interp, cmd.c_str());
    if (result != TCL_OK) {
        const char* error_msg = Tcl_GetStringResult(interp);
        out << "Error generating timing report: " << error_msg << std::endl;
        out << "Command attempted: " << cmd << std::endl;
    } else {
        // 获取命令输出并写入文件
        const char* output = Tcl_GetStringResult(interp);
        out << output;
    }
    
    out.close();
}

// ============================================================================
// 工具函数
// ============================================================================

sta::dbNetwork* TimingEngine3D::getNetwork() const
{
    if (!sta_) return nullptr;
    return sta_->getDbNetwork();
}

sta::Scene* TimingEngine3D::getDefaultScene() const
{
    if (!sta_) return nullptr;
    return sta_->cmdScene();
}

const sta::MinMax* TimingEngine3D::getDefaultMinMax() const
{
    if (!sta_) return nullptr;
    
    // 默认使用max模式
    return sta::MinMax::max();
}


// ============================================================================
// 私有辅助函数
// ============================================================================

bool TimingEngine3D::setupAnalysis()
{
    if (!sta_) return false;
    
    // 确保图已构建
    ensureGraphBuilt();
    
    // 确保已层级化
    ensureLevelized();
    
    return true;
}

TimingEngine3D::DesignStats TimingEngine3D::getDesignStats() const
{
    DesignStats stats;
    
    if (!sta_) return stats;
    
    sta::dbNetwork* network = getNetwork();
    if (!network) return stats;
    
    // 获取顶层实例
    sta::Instance* top_inst = network->topInstance();
    if (!top_inst) return stats;
    
    // 获取设计统计信息
    sta::PatternMatch all_pattern("*");
    
    // 获取所有实例
    sta::InstanceSeq insts = network->findInstancesMatching(top_inst, &all_pattern);
    stats.num_instances = insts.size();
    
    // 获取所有网络
    sta::NetSeq nets = network->findNetsMatching(top_inst, &all_pattern);
    stats.num_nets = nets.size();
    
    // 获取引脚数量
    stats.num_pins = 0;
    for (const sta::Instance* inst : insts) {
        sta::PinSeq pins = network->findPinsMatching(const_cast<sta::Instance*>(inst), &all_pattern);
        stats.num_pins += pins.size();
    }
    
    // 获取scene数量
    stats.num_corners = sta_->scenes().size();
    
    // 设计面积（简化实现）
    stats.design_area = 0.0;
    
    return stats;
}

void TimingEngine3D::updateGraph()
{
    if (!sta_) return;
    
    sta_->ensureGraph();
}

void TimingEngine3D::ensureGraphBuilt()
{
    if (!sta_) return;
    
    sta_->ensureGraph();
}

void TimingEngine3D::ensureLevelized()
{
    if (!sta_) return;
    
    sta_->ensureLevelized();
}

void TimingEngine3D::processTimingResults()
{
    // 处理时序结果
    // 可以在这里添加自定义处理逻辑
}

sta::Pin* TimingEngine3D::findPinByName(const std::string& pin_name) const
{
    if (!sta_) return nullptr;
    
    sta::dbNetwork* network = getNetwork();
    if (!network) return nullptr;
    
    // 尝试查找引脚
    sta::Pin* pin = network->findPin(pin_name.c_str());
    return pin;
}

sta::Net* TimingEngine3D::findNetByName(const std::string& net_name) const
{
    if (!sta_) return nullptr;
    
    sta::dbNetwork* network = getNetwork();
    if (!network) return nullptr;
    
    // 尝试查找网络
    sta::Net* net = network->findNet(net_name.c_str());
    return net;
}

void TimingEngine3D::debugPrintSetup() const
{
    std::cout << "=== TimingEngine3D Setup ===" << std::endl;
    std::cout << "  STA: " << (sta_ ? "Set" : "Not set") << std::endl;
    std::cout << "  Timingdp: " << (timingdp_ ? "Set" : "Not set") << std::endl;
    
    // 设计统计
    DesignStats design_stats = getDesignStats();
    std::cout << "  Design stats:" << std::endl;
    std::cout << "    Instances: " << design_stats.num_instances << std::endl;
    std::cout << "    Nets: " << design_stats.num_nets << std::endl;
    std::cout << "    Pins: " << design_stats.num_pins << std::endl;
    std::cout << "    Scenes: " << design_stats.num_corners << std::endl;
}

sta::Net* TimingEngine3D::getStaNetFromDbNet(odb::dbNet* db_net) const
{
    if (!db_net || !sta_) return nullptr;
    
    sta::dbNetwork* network = getNetwork();
    if (!network) return nullptr;
    
    // 通过名称查找sta::Net
    const char* net_name = db_net->getConstName();
    return network->findNet(net_name);
}


// 获取引脚的完整时序信息
TimingEngine3D::PinTimingInfo TimingEngine3D::getPinTimingInfo(const std::string& pin_name)
{
    sta::Pin* pin = findPinByName(pin_name);
    if (!pin) {
        throw std::runtime_error("Pin '" + pin_name + "' not found");
    }
    return getPinTimingInfo(pin);
}

TimingEngine3D::PinTimingInfo TimingEngine3D::getPinTimingInfo(sta::Pin* pin)
{
    PinTimingInfo info = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false};
    
    if (!pin || !sta_) {
        return info;
    }
    
    // 获取顶点
    sta::Vertex* vertex = getVertexFromPin(pin);
    if (!vertex) {
        return info;
    }
    
    // 获取默认的上升/下降边沿
    const sta::RiseFall* rf = getDefaultRiseFall();
    const sta::RiseFallBoth* rf_both = rf->asRiseFallBoth();

    // 获取默认scene
    sta::Scene* scene = getDefaultScene();
    if (!scene) {
        return info;
    }
    
    try {
        // 获取到达时间
        info.early_arrival = sta_->arrival(vertex, rf_both, sta_->scenes(), sta::MinMax::min());
        info.late_arrival = sta_->arrival(vertex, rf_both, sta_->scenes(), sta::MinMax::max());
        
        // 获取要求时间
        info.early_required = sta_->required(vertex, rf_both, sta_->scenes(), sta::MinMax::min());
        info.late_required = sta_->required(vertex, rf_both, sta_->scenes(), sta::MinMax::max());
        
        // 计算松弛度
        info.early_slack = info.early_arrival - info.early_required;
        info.late_slack = info.late_required - info.late_arrival;
        
        info.is_valid = true;
    } catch (const std::exception& e) {
        std::cerr << "Error getting timing info for pin: " 
                  << getNetwork()->name(pin) << " - " << e.what() << std::endl;
    }
    
    return info;
}

TimingEngine3D::PinTimingInfo TimingEngine3D::getPinTimingInfoDirect(sta::Pin* pin)
{
    PinTimingInfo info = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false};
    if (!pin || !sta_) return info;

    sta::Vertex* vertex = getVertexFromPin(pin);
    if (!vertex) return info;

    const sta::RiseFall* rf = getDefaultRiseFall();

    sta::Path* early_arr_path = sta_->vertexWorstArrivalPath(vertex, rf, sta::MinMax::min());
    sta::Path* late_arr_path  = sta_->vertexWorstArrivalPath(vertex, rf, sta::MinMax::max());
    sta::Path* early_req_path = sta_->vertexWorstRequiredPath(vertex, rf, sta::MinMax::min());
    sta::Path* late_req_path  = sta_->vertexWorstRequiredPath(vertex, rf, sta::MinMax::max());

    if (early_arr_path) info.early_arrival  = sta::delayAsFloat(early_arr_path->arrival());
    if (late_arr_path)  info.late_arrival   = sta::delayAsFloat(late_arr_path->arrival());
    if (early_req_path) info.early_required = sta::delayAsFloat(early_req_path->required());
    if (late_req_path)  info.late_required  = sta::delayAsFloat(late_req_path->required());

    info.early_slack = info.early_arrival - info.early_required;
    info.late_slack  = info.late_required - info.late_arrival;
    info.is_valid    = true;
    return info;
}

// 获取引脚到达时间
double TimingEngine3D::getPinArrivalTime(const std::string& pin_name, bool is_early)
{
    sta::Pin* pin = findPinByName(pin_name);
    if (!pin) {
        throw std::runtime_error("Pin '" + pin_name + "' not found");
    }
    return getPinArrivalTime(pin, is_early);
}

double TimingEngine3D::getPinArrivalTime(sta::Pin* pin, bool is_early)
{
    if (!pin || !sta_) {
        return 0.0;
    }
    
    sta::Vertex* vertex = getVertexFromPin(pin);
    if (!vertex) {
        return 0.0;
    }
    
    const sta::RiseFall* rf = getDefaultRiseFall();
    return sta_->arrival(vertex, rf->asRiseFallBoth(), sta_->scenes(), getMinMax(is_early));
}

// 获取引脚要求时间
double TimingEngine3D::getPinRequiredTime(const std::string& pin_name, bool is_early)
{
    sta::Pin* pin = findPinByName(pin_name);
    if (!pin) {
        throw std::runtime_error("Pin '" + pin_name + "' not found");
    }
    return getPinRequiredTime(pin, is_early);
}

double TimingEngine3D::getPinRequiredTime(sta::Pin* pin, bool is_early)
{
    if (!pin || !sta_) {
        return 0.0;
    }
    
    sta::Vertex* vertex = getVertexFromPin(pin);
    if (!vertex) {
        return 0.0;
    }
    
    const sta::RiseFall* rf = getDefaultRiseFall();
    const sta::MinMax* minmax = getMinMax(is_early);
    
    return sta_->required(vertex, rf->asRiseFallBoth(), sta_->scenes(), minmax);
}

// 获取引脚松弛度
double TimingEngine3D::getPinSlack(const std::string& pin_name, bool is_early)
{
    sta::Pin* pin = findPinByName(pin_name);
    if (!pin) {
        throw std::runtime_error("Pin '" + pin_name + "' not found");
    }
    return getPinSlack(pin, is_early);
}

double TimingEngine3D::getPinSlack(sta::Pin* pin, bool is_early)
{
    if (!pin || !sta_) {
        return 0.0;
    }
    
    const sta::MinMax* minmax = getMinMax(is_early);
    return sta_->slack(pin, sta::RiseFallBoth::riseFall(), sta_->scenes(), minmax);
}

//获取cell内弧延迟时间
TimingEngine3D::ArcDelayInfo TimingEngine3D::getArcDelay(
    const std::string& from_pin_name, 
    const std::string& to_pin_name) {
    
    ArcDelayInfo info = {false, 0.0, 0.0, from_pin_name, to_pin_name};
    
    sta::Pin* from_pin = findPinByName(from_pin_name);
    sta::Pin* to_pin = findPinByName(to_pin_name);
    
    if (!from_pin || !to_pin) {
        return info;
    }
    
    return getArcDelay(from_pin, to_pin);
}

TimingEngine3D::ArcDelayInfo TimingEngine3D::getArcDelay(
    sta::Pin* from_pin, 
    sta::Pin* to_pin) {
    
    ArcDelayInfo info = {false, 0.0, 0.0, 
                        getNetwork()->name(from_pin), 
                        getNetwork()->name(to_pin)};
    
    if (!from_pin || !to_pin || !sta_) {
        return info;
    }
    
    // 获取图
    sta::Graph* graph = sta_->ensureGraph();
    if (!graph) {
        return info;
    }
    
    // 获取顶点
    sta::Vertex* from_vertex = getVertexFromPin(from_pin);
    sta::Vertex* to_vertex = graph->pinDrvrVertex(to_pin);
    
    // 查找从 from_vertex 到 to_vertex 的边
    sta::VertexOutEdgeIterator edge_iter(from_vertex, graph);
    while (edge_iter.hasNext()) {
        sta::Edge* edge = edge_iter.next();
        
        if (edge->to(graph) == to_vertex) {
            // 找到边，获取时序弧集
            sta::TimingArcSet* arc_set = edge->timingArcSet();
            if (!arc_set) {
                continue;
            }
            
            // 获取默认的scene
            sta::Scene* scene = getDefaultScene();
            if (!scene) {
                continue;
            }
            
            // 获取早时序分析点index（min）
            int early_ap_index = scene->dcalcAnalysisPtIndex(sta::MinMax::min());
            // 获取晚时序分析点index（max）
            int late_ap_index = scene->dcalcAnalysisPtIndex(sta::MinMax::max());
            
            // 遍历所有时序弧，取最坏延迟
            sta::TimingArcSeq& arcs = arc_set->arcs();
            if (arcs.empty()) {
                continue;
            }
            
            // 获取弧延迟
            try {
                for (sta::TimingArc* arc : arcs) {
                    info.early_delay = std::max(info.early_delay,
                        static_cast<double>(delayAsFloat(sta_->arcDelay(edge, arc, early_ap_index))));
                    info.late_delay  = std::max(info.late_delay,
                        static_cast<double>(delayAsFloat(sta_->arcDelay(edge, arc, late_ap_index))));
                }
                info.is_valid = true;
                
                // 调试输出
                // std::cout << "Found arc delay from " << info.from_pin_name 
                //           << " to " << info.to_pin_name
                //           << ": early=" << info.early_delay 
                //           << ", late=" << info.late_delay << std::endl;
                
                return info;
            } catch (const std::exception& e) {
                std::cerr << "Error getting arc delay: " << e.what() << std::endl;
            }
        }
    }
    
    return info;
}


// 私有辅助方法
sta::Vertex* TimingEngine3D::getVertexFromPin(sta::Pin* pin) const
{
    if (!sta_ || !pin) {
        return nullptr;
    }
    
    sta::Graph* graph = sta_->ensureGraph();
    if (!graph) {
        return nullptr;
    }
    
    return graph->pinLoadVertex(pin);
}

sta::Vertex* TimingEngine3D::getStaVertex(Node3D* cell) const
{
    if (!cell || !cell->db_inst) return nullptr;

    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Instance* inst = db_network->dbToSta(cell->db_inst);
    if (!inst) return nullptr;

    sta::InstancePinIterator* pin_iter = db_network->pinIterator(inst);
    sta::Pin* out_pin = nullptr;
    while (pin_iter->hasNext()) {
        sta::Pin* pin = pin_iter->next();
        if (db_network->direction(pin)->isOutput()) {
            out_pin = pin;
            break;
        }
    }
    delete pin_iter;

    if (!out_pin) return nullptr;
    return sta_->graph()->pinDrvrVertex(out_pin);
}

const sta::MinMax* TimingEngine3D::getMinMax(bool is_early) const
{
    // early对应min，late对应max
    return is_early ? sta::MinMax::min() : sta::MinMax::max();
}

const sta::RiseFall* TimingEngine3D::getDefaultRiseFall() const
{
    // 使用上升边沿作为默认
    return sta::RiseFall::rise();
}


} // namespace tdp
