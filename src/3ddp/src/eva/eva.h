// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include <memory>
#include <vector>

#include "stt/SteinerTreeBuilder.h"
#include "tdp/Timingdp.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Pin3D.h"
#include "odb/geom.h"
#include "grt/GRoute.h"
#include "est/SteinerTree.h" 
#include "utl/Logger.h"

#include <utility> 


namespace sta {
class dbSta;
class dbNetwork;
class Network;
class Parasitics;
class Net;
class Pin;
class Scene;
class Parasitic;
class ParasiticNode;
class ArcDelayCalc;
}
namespace est {
class EstimateParasitics;
}

namespace tdp {

using NetRouteMap = std::map<odb::dbNet*, grt::GRoute>;

class Net3D;
class Network3D;
class Pointdp3D;

// 每线程私有的工作上下文（用于 OpenMP 并行 generateSteinerTrees）。
// 定义见 eva.cpp。
struct WorkerCtx;

/**
 * @class SteinerTreeGenerator
 * @brief 为tdp的3D网络生成斯坦纳树
 * 
 * 功能：
 * 1. 单层net的生成逻辑与grt相同（使用FLUTE算法）
 * 2. 多层net投影到同一层
 * 3. 加入事先分配好的HBT位置作为额外斯坦纳点
 */

struct PinCoordinates {
    std::vector<int> x_coords;
    std::vector<int> y_coords;
    int drvr_idx = -1;  // -1 表示未找到驱动引脚
    std::vector<int> layer_coords; // 
    std::vector<const sta::Pin*> pins;
    std::vector<int> conn_layers; // 每个 pin 的连接层（缓存 getPinConnLayer 结果）
    
    bool isValid() const {
        return !x_coords.empty() && 
               x_coords.size() == y_coords.size() &&
               x_coords.size() == layer_coords.size() &&
               drvr_idx >= 0 && 
               drvr_idx < static_cast<int>(x_coords.size());
    }
    
    size_t pinCount() const { return x_coords.size(); }
};

struct Pin3Dlocation{
    int x;
    int y;
    int z;

    Pin3Dlocation() : x(0), y(0), z(0) {}
    Pin3Dlocation(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}
};

 
class SteinerTreeGenerator {
public:
    SteinerTreeGenerator(Timingdp* tdp_db);
    ~SteinerTreeGenerator();
    
    void createSteninerRoute();
    bool updateNetRouteAndParasitics(Net3D* net_3d);
    void updateSteinerTrees(Node3D* cell);
    // 为Net3D生成斯坦纳树
    void generateSteinerTree(Net3D* net_3d, const sta::Pin* drvr_pin, WorkerCtx& ctx);

    Pin3Dlocation getPin3DLocation(const Pin3D* pin) const; 
    int findPinLayer(const est::PinLoc& pinloc, const Net3D* net_3d);
    int getPinConnLayer(const sta::Pin* sta_pin) const;
    void estimateTdpWireParasitic(odb::dbNet* db_net, grt::GRoute& route,
                                  const PinCoordinates& pin_locs,
                                  const std::vector<odb::Point>& all_locs,
                                  WorkerCtx& ctx);
    void writeSpef(const std::string& base_filename);
    void clearAllParasitics();
    sta::ParasiticNode* ensureTdpParasiticNode(
        sta::Parasitics* parasitics, sta::Parasitic* parasitic,
        int x, int y, int layer,
        std::map<grt::RoutePt, sta::ParasiticNode*>& node_map,
        sta::Net* net) const;
    void tdpLayerRC(int wire_length_dbu,
                    int init_x, int init_y, int final_x, int final_y,
                    int layer_id,
                    const sta::Scene* corner,
                     float& res, float& cap) const;
    // 跨层 segment 是否穿过 hb_layer（即 HBT via）。rl_a/rl_b 为布线层序号。
    bool crossesHbtLayer(int rl_a, int rl_b) const;
    // 普通金属 via 电阻：累加该 via 跨越的各 cut 层电阻（值取自 est 的 setRC 表）。
    // rl_a/rl_b 为布线层序号（getPinConnLayer / 代表布线层）。
    double metalViaResistance(int rl_a, int rl_b, const sta::Scene* corner) const;
    double dbuToMeters(int dbu) const;

    //转换route
    grt::GRoute convertSteinerTreeToGRoute(const stt::Tree& tree, int routing_layer);

    // grt::GRoute convertSteinerTreeToGRoute(const stt::Tree& tree, 
    //                                       const std::vector<odb::Point>& pin_locs,
    //                                       int layer = 0) const;
    // // 为Net3D生成GRoute（考虑跨层情况）
    // std::pair<grt::GRoute, grt::GRoute> generateGRouteForNet3D(Net3D* net_3d);
    
    // 批量生成斯坦纳树
    void generateSteinerTrees(std::vector<odb::dbNet*>& nets);
    
    // 获取生成的斯坦纳树统计信息
    struct SteinerTreeStats {
        int num_nets_processed = 0;
        int total_pins = 0;
        int total_hbts = 0;
        int total_steiner_points = 0;
        double total_wirelength = 0.0;

        void merge(const SteinerTreeStats& other) {
            num_nets_processed += other.num_nets_processed;
            total_pins += other.total_pins;
            total_hbts += other.total_hbts;
            total_steiner_points += other.total_steiner_points;
            total_wirelength += other.total_wirelength;
        }
    };
    
    SteinerTreeStats getStats() const { return stats_; }

    grt::GRoute generateLocalRoute(Net3D* net_3d, const Pin3D* moved_pin, int new_x, int new_y);
    // Build the net's wire parasitic network (current geometry, or a trial
    // geometry when moved_pin is given) into an arbitrary Parasitics container
    // instead of the global scene. No reduceParasitic(), no global scene/route
    // writes. moved_pin/new_x/new_y override that pin's position (odb is NOT
    // touched) — used by the per-thread trial parasitics path.
    void buildLocalNetParasiticInto(sta::Parasitics* target, Net3D* net_3d,
                                    const sta::Scene* corner,
                                    const Pin3D* moved_pin = nullptr,
                                    int new_x = 0, int new_y = 0);
    void getSegmentRC(int wire_length_dbu, int init_x, int init_y, int final_x, int final_y, int layer_id, float& res, float& cap) const;
    
private:
   
    
    PinCoordinates estimateWireSteiner(const sta::Net* net, const Net3D* net_3d);
    PinCoordinates estimateWireSteiner(const sta::Pin* drvr_pin, const sta::Net* net, const Net3D* net_3d);
    PinCoordinates collectPinLocations(const sta::Pin* drvr_pin, const Net3D* net_3d);

    // PinCoordinates estimateWireSteiner(const sta::Net* net);
    // PinCoordinates estimateWireSteiner(const sta::Pin* drvr_pin, const sta::Net* net);
    // // 收集引脚位置（考虑多层）
    // PinCoordinates collectPinLocations(const sta::Pin* drvr_pin);
    
    // 投影到单层（处理3D坐标到2D投影）
    std::vector<odb::Point> projectToSingleLayer(const std::vector<Pointdp3D>& pin_locs_3d);
    
    // 添加HBT位置
    std::vector<odb::Point> addHBTLocations(Net3D* net_3d, 
                                          const std::vector<odb::Point>& pin_locs);
    
    // 调用STT生成斯坦纳树
    stt::Tree callFluteAlgorithm(const std::vector<int>& x_coords,
                                const std::vector<int>& y_coords,
                                int drvr_idx,
                                WorkerCtx& ctx);
    
    //生成布线

    // 验证斯坦纳树
    // bool validateSteinerTree(const stt::Tree& tree, Net3D* net_3d);
    bool isInterLayerNetwork(Net3D* net) const ;
    
    void initSceneParasitics();
    
    bool isPadNet(const sta::Net* net) const;
    bool isPadPin(const sta::Pin* pin) const;
    bool isPad(const sta::Instance* inst) const;
    void net2Pins(const sta::Net* net, const sta::Pin*& pin1, const sta::Pin*& pin2) const;
    void makePadParasitic(const sta::Net* net, WorkerCtx& ctx);
    
    Timingdp* tdp_db_;
    std::unique_ptr<stt::SteinerTreeBuilder> stt_builder_;
    utl::Logger* logger_;
    SteinerTreeStats stats_;
    sta::dbSta* sta_ = nullptr;
    grt::GlobalRouter* global_router_ = nullptr;
    sta::Network* network_ = nullptr;
    sta::dbNetwork* db_network_ = nullptr;
    odb::dbBlock* block_ = nullptr;
    odb::dbTech* tech_ = nullptr;
    const sta::MinMax* min_max_ = nullptr;
    est::EstimateParasitics* est_ = nullptr;
};

} // namespace tdp
