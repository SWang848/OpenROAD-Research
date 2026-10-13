// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "placement.h"
#include "MIS3D.h"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "sta/Path.hh"
#include "sta/Sdc.hh"
#include "../eva/eva.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Pin3D.h"
#include "../infrastructure/Common3D.h"

#include <iostream>
#include <cstdlib>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <random>
#include <unordered_map>
#include <limits>
#include <queue>
#include <stack>
#include <list>
#include <vector>
#include <string>

using utl::TDP;

namespace tdp {

namespace {

bool dfsFindCycle(
    int v,
    std::vector<int>& parent,
    std::vector<int>& visited_gen,
    std::vector<int>& on_path_gen,
    int dfs_gen,
    const std::vector<std::vector<int>>& forward_of,
    const std::vector<char>& remaining,
    int& cycle_start,
    int& cycle_end)
{
    visited_gen[v] = dfs_gen;
    on_path_gen[v] = dfs_gen;

    for (int w : forward_of[v]) {
        if (!remaining[w]) continue;
        if (visited_gen[w] != dfs_gen) {
            parent[w] = v;
            if (dfsFindCycle(w, parent, visited_gen, on_path_gen, dfs_gen,
                             forward_of, remaining, cycle_start, cycle_end))
                return true;
        } else if (on_path_gen[w] == dfs_gen) {
            cycle_start = w;
            cycle_end = v;
            return true;
        }
    }

    on_path_gen[v] = 0;
    return false;
}

} // anonymous namespace


Placement3D::Placement3D(Timingdp* timingdp)
    : timingdp_(timingdp),
      sta_(timingdp->getSta()),
      logger_(timingdp ? timingdp->getLogger() : nullptr),
      timing_engine_(timingdp->gettimingengine()),
      clock_period_T_(0.0)
{
    if (sta_) {
        auto& clks = sta_->cmdSdc()->clocks();
        clock_period_T_ = std::numeric_limits<double>::max();
        for (auto* clk : clks)
            clock_period_T_ = std::min(clock_period_T_, static_cast<double>(clk->period()));
        if (clock_period_T_ == std::numeric_limits<double>::max())
            clock_period_T_ = 1.0;
    }
    if (timingdp_) {
        arch_ = timingdp_->getArchitecture3D();
        network_ = timingdp_->getNetwork3D();
    }
}

Placement3D::~Placement3D()
{
    // 清理资源
}

void Placement3D::buildTraverseOrder(Timingdp* timingdp) {
    traverse_list.clear();
    if (!timingdp) {
        logger_->error(TDP, 1300, "Error: Timingdp is null");
        return;
    }

    auto* network_3d = timingdp->getNetwork3D();
    if (!network_3d) {
        logger_->error(TDP, 1301, "Error: Network3D is null");
        return;
    }

    const auto& all_cells = network_3d->getCells();
    const int N = static_cast<int>(all_cells.size());

    std::unordered_map<Node3D*, int> idx_map;
    for (int i = 0; i < N; ++i) {
        idx_map[all_cells[i]] = i;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    // ================================================================
    // Phase 1: out_degree + reverse adjacency (drivers_of)
    // ================================================================
    std::vector<int> out_degree(N, 0);
    std::vector<std::vector<int>> drivers_of(N);
    std::vector<int> downstream_marker(N, -1);

    for (int cell_idx = 0; cell_idx < N; ++cell_idx) {
        auto* cell = all_cells[cell_idx];
        if (!cell) continue;

        int marker = cell_idx + 1;
        int downstream_count = 0;

        if (cell->isTerminal()) {
            bool is_pi = false;
            for (auto* pin : cell->getPins()) {
                if (pin && isPrimaryInput(pin)) { is_pi = true; break; }
            }

            if (is_pi) {
                for (auto* pin : cell->getPins()) {
                    if (!pin) continue;
                    Net3D* net = pin->getNet();
                    if (!net) continue;
                    if (net->db_net && (timingdp->isIdealClockNet(net->db_net) || net->db_net->getSigType().isSupply() || net->db_net->isSpecial())) continue;
                    for (auto* net_pin : net->getPins()) {
                        if (net_pin == pin) continue;
                        Node3D* sink = net_pin->getNode();
                        if (!sink || sink == cell) continue;
                        auto it = idx_map.find(sink);
                        if (it == idx_map.end()) continue;
                        int sink_idx = it->second;
                        if (downstream_marker[sink_idx] != marker) {
                            downstream_marker[sink_idx] = marker;
                            drivers_of[sink_idx].push_back(cell_idx);
                            downstream_count++;
                        }
                    }
                }
            }
        } else if (isSequentialCell(cell)) {
            for (auto* pin : cell->getPins()) {
                if (!pin || pin->getDirection() != Pin3D::Dir_OUT) continue;
                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && (timingdp->isIdealClockNet(net->db_net) || net->db_net->getSigType().isSupply() || net->db_net->isSpecial())) continue;
                for (auto* net_pin : net->getPins()) {
                    if (net_pin == pin) continue;
                    if (net_pin->getDirection() != Pin3D::Dir_IN && !isPrimaryOutput(net_pin)) continue;
                    Node3D* sink = net_pin->getNode();
                    if (!sink || sink == cell) continue;
                    auto it = idx_map.find(sink);
                    if (it == idx_map.end()) continue;
                    int sink_idx = it->second;
                    if (downstream_marker[sink_idx] != marker) {
                        downstream_marker[sink_idx] = marker;
                        drivers_of[sink_idx].push_back(cell_idx);
                    }
                }
            }
        } else {
            for (auto* pin : cell->getPins()) {
                if (!pin || pin->getDirection() != Pin3D::Dir_OUT) continue;
                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && (timingdp->isIdealClockNet(net->db_net) || net->db_net->getSigType().isSupply() || net->db_net->isSpecial())) continue;
                for (auto* net_pin : net->getPins()) {
                    if (net_pin == pin) continue;
                    if (net_pin->getDirection() != Pin3D::Dir_IN && !isPrimaryOutput(net_pin)) continue;
                    Node3D* sink = net_pin->getNode();
                    if (!sink || sink == cell) continue;
                    auto it = idx_map.find(sink);
                    if (it == idx_map.end()) continue;
                    int sink_idx = it->second;
                    if (downstream_marker[sink_idx] != marker) {
                        downstream_marker[sink_idx] = marker;
                        drivers_of[sink_idx].push_back(cell_idx);
                        downstream_count++;
                    }
                }
            }
        }
        out_degree[cell_idx] = downstream_count;
    }

    std::vector<std::vector<int>> forward_of(N);
    for (int sink_idx = 0; sink_idx < N; ++sink_idx) {
        for (int driver_idx : drivers_of[sink_idx]) {
            forward_of[driver_idx].push_back(sink_idx);
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    logger_->info(TDP, 114, "start buildTraverseOrder Phase 2");

    // ================================================================
    // Phase 2: Kahn
    // ================================================================
    std::queue<int> q;
    std::vector<char> visited(N, 0);
    std::vector<char> remaining(N, 0);
    int remaining_count = 0;

    for (int cell_idx = 0; cell_idx < N; ++cell_idx) {
        auto* cell = all_cells[cell_idx];
        if (out_degree[cell_idx] == 0) {
            q.push(cell_idx);
        }
        if (out_degree[cell_idx] > 0 && !isSequentialCell(cell)) {
            if (cell->isTerminal()) continue;
            remaining[cell_idx] = 1;
            remaining_count++;
        }
    }

    int cycle_break_count = 0;

    std::vector<int> parent(N, -1);
    std::vector<int> visited_gen(N, 0);
    std::vector<int> on_path_gen(N, 0);
    int dfs_gen = 0;

    while (!q.empty()) {
        int cell_idx = q.front();
        q.pop();

        traverse_list.push_back(all_cells[cell_idx]);
        visited[cell_idx] = 1;
        if (remaining[cell_idx]) {
            remaining[cell_idx] = 0;
            remaining_count--;
        }

        for (int driver_idx : drivers_of[cell_idx]) {
            if (--out_degree[driver_idx] == 0 && !visited[driver_idx]) {
                if (!isSequentialCell(all_cells[driver_idx])) {
                    q.push(driver_idx);
                }
            }
        }

        if (q.empty() && remaining_count > 0) {
            dfs_gen++;
            int cycle_start = -1, cycle_end = -1;

            for (int i = 0; i < N; ++i) {
                if (remaining[i] && visited_gen[i] != dfs_gen) {
                    if (dfsFindCycle(i, parent, visited_gen, on_path_gen, dfs_gen,
                                     forward_of, remaining, cycle_start, cycle_end))
                        break;
                }
            }

            if (cycle_start >= 0) {
                std::vector<int> cycle_nodes;
                for (int curr = cycle_end; ; curr = parent[curr]) {
                    cycle_nodes.push_back(curr);
                    if (curr == cycle_start) break;
                }

                int cut_idx = -1;
                int min_deg = std::numeric_limits<int>::max();
                for (int c : cycle_nodes) {
                    if (out_degree[c] == 1) { cut_idx = c; break; }
                    if (out_degree[c] < min_deg) min_deg = out_degree[c];
                }
                if (cut_idx == -1) {
                    for (int c : cycle_nodes) {
                        if (out_degree[c] == min_deg) { cut_idx = c; break; }
                    }
                }

                if (cut_idx >= 0) {
                    int target_to_erase = -1;
                    for (int t : forward_of[cut_idx]) {
                        if (remaining[t]) { target_to_erase = t; break; }
                    }
                    if (target_to_erase >= 0) {
                        auto& dvec = drivers_of[target_to_erase];
                        auto dit = std::find(dvec.begin(), dvec.end(), cut_idx);
                        if (dit != dvec.end()) dvec.erase(dit);
                        auto& fvec = forward_of[cut_idx];
                        auto fit = std::find(fvec.begin(), fvec.end(), target_to_erase);
                        if (fit != fvec.end()) fvec.erase(fit);
                    }

                    out_degree[cut_idx] = 0;
                    q.push(cut_idx);
                    remaining[cut_idx] = 0;
                    remaining_count--;
                    cycle_break_count++;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    if (remaining[i]) {
                        out_degree[i] = 0;
                        q.push(i);
                        remaining[i] = 0;
                        remaining_count--;
                    }
                }
            }
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    int total_cells = network_->getCells().size();
    int traversed = traverse_list.size();
    std::cout << "Traverse list: " << traversed 
              << " / Total cells: " << total_cells 
              << " / Remaining: " << remaining_count
              << " / Cycles broken: " << cycle_break_count << std::endl;

    logger_->info(TDP, 116, "buildTraverseOrder Phase1: {}ms Phase2: {}ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(),
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());

    if (traverse_list.empty()) {
        logger_->warn(TDP, 1303, "WARNING: Traverse list is empty!");
    } else {
        logger_->info(TDP, 1304, "Reverse topological order built, {} cells", traverse_list.size());
    }
}


bool Placement3D::isInputPin(const Pin3D* pin) {
    return pin->getDirection() == Pin3D::Dir_IN;
}

// 检查引脚是否是输出引脚
bool Placement3D::isOutputPin(const Pin3D* pin) {
    return pin->getDirection() == Pin3D::Dir_OUT;
}

// 检查引脚是否是主输入（边界输入引脚）
bool Placement3D::isPrimaryInput(const Pin3D* pin) {
    // 边界引脚且方向为输入
    return pin->getDbBTerm() != nullptr && 
           pin->getDirection() == Pin3D::Dir_IN;
}

// 检查引脚是否是主输出（边界输出引脚）
bool Placement3D::isPrimaryOutput(const Pin3D* pin) {
    // 边界引脚且方向为输出
    return pin->getDbBTerm() != nullptr && 
           pin->getDirection() == Pin3D::Dir_OUT;
}

bool Placement3D::isPrimaryInput(const Node3D* cell) const {
    if (!cell || !cell->isTerminal()) return false;
    if (!cell->db_bterm_) return false;
    return cell->db_bterm_->getIoType() == odb::dbIoType::INPUT;
}

bool Placement3D::isPrimaryOutput(const Node3D* cell) const {
    if (!cell || !cell->isTerminal()) return false;
    if (!cell->db_bterm_) return false;
    return cell->db_bterm_->getIoType() == odb::dbIoType::OUTPUT;
}

// 判断一个 cell 是否是时序元件（寄存器/FF/latch）
bool Placement3D::isSequentialCell(const Node3D* cell) const{
    if (!cell || cell->isTerminal()) return false;
    if (cell->db_inst && cell->db_inst->getMaster()) {
        return cell->db_inst->getMaster()->isSequential();
    }
    return false;
}

//find cell row 
Architecture3D::Row3D* Placement3D::findCellRow(Node3D* cell, 
                                                const std::vector<Architecture3D::Row3D*>& rows) {
    
    if (!cell || rows.empty()) {
        return nullptr;
    }
    // 获取cell的位置和尺寸
    Pointdp3D cell_loc = cell->location;
    DBU cell_width = cell->getWidth().v;
    DBU cell_height = cell->getHeight().v;
    
    // Prefer an exact row match for an already legalized cell.
    for (auto* row : rows) {
        if (!row) continue;
        
        // 获取行的边界
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        
        // 检查cell是否在行内（考虑cell的尺寸）
        bool in_y = (cell_loc.y.v >= row_bottom) && 
                    (cell_loc.y.v + cell_height <= row_top);
        
        bool in_x = (cell_loc.x.v >= row_left) && 
                    (cell_loc.x.v + cell_width <= row_right);
        
        if (in_y && in_x) {
            return row;  // 找到cell所在的行
        }
    }

    // The gp2lg handoff contains global-placement coordinates, so a standard
    // cell can straddle two rows before detailed placement.  Use the legal
    // row with the greatest vertical overlap as its search anchor instead of
    // rejecting the cell outright.  Candidate locations remain row-aligned.
    Architecture3D::Row3D* best_row = nullptr;
    DBU best_overlap = -1;
    DBU best_origin_distance = std::numeric_limits<DBU>::max();
    const DBU cell_top = cell_loc.y.v + cell_height;

    for (auto* row : rows) {
        if (!row) continue;

        const DBU row_left = row->getLeft().v;
        const DBU row_right = row->getRight().v;
        if (cell_loc.x.v < row_left || cell_loc.x.v + cell_width > row_right) {
            continue;
        }

        const DBU row_bottom = row->getBottom().v;
        const DBU row_top = row->getTop().v;
        const DBU overlap = std::max<DBU>(
            0, std::min(cell_top, row_top) - std::max(cell_loc.y.v, row_bottom));
        const DBU origin_distance = std::abs(cell_loc.y.v - row_bottom);

        if (overlap > best_overlap
            || (overlap == best_overlap && origin_distance < best_origin_distance)) {
            best_row = row;
            best_overlap = overlap;
            best_origin_distance = origin_distance;
        }
    }

    return best_row;
}


void Placement3D::rundetailedplacement(){


        LRParams params;
        params.lambda_wirelength = 1.0;
        params.lambda_timing = 1.0;
        params.lambda_alignment = 0.5;
        params.max_iterations = 10;
        params.convergence_threshold = 0.001;
        params.window_y = 10;
        params.window_x = 20;
        

        LSP3D tdp_lsp(timingdp_, params);
        DetailedMis3D tdp_mis(timingdp_, &tdp_lsp);
        buildTraverseOrder(timingdp_);
        // verifyTraverseOrder();
        tdp_lsp.traverse_list = traverse_list;

        // 调用 initWeight 函数
        tdp_lsp.place();

        // tdp_mis.run();
        //test file

        
}
// void placeall (){


// }

void Placement3D::verifyTraverseOrder() const
{
    if (!network_) {
        if (logger_)
            logger_->warn(TDP, 9100, "verifyTraverseOrder: network_ is null");
        return;
    }

    // 1. 建位置索引
    std::unordered_map<const Node3D*, int> pos;
    int idx = 0;
    for (auto* cell : traverse_list) {
        pos[cell] = idx++;
    }

    // 2. 统计
    int total_arcs      = 0;
    int violations      = 0;
    int skipped_fixed   = 0;
    int skipped_missing = 0;

    int cur_pos = 0;
    for (auto* cell : traverse_list) {
        if (!cell || cell->is_fixed ) {
            cur_pos++;
            continue;
        }

        for (auto* pin : cell->getPins()) {
            if (!pin) continue;
            if (pin->getDirection() != Pin3D::Dir_IN) continue;

            Net3D* net = pin->getNet();
            if (!net) continue;

            // 在 net 上找 driver pin（输出引脚）
            Pin3D* driver_pin = nullptr;
            for (auto* net_pin : net->getPins()) {
                if (!net_pin || net_pin == pin) continue;
                if (net_pin->getDirection() == Pin3D::Dir_OUT) {
                    driver_pin = net_pin;
                    break;
                }
            }
            if (!driver_pin) continue;

            Node3D* driver_cell = driver_pin->getNode();
            if (!driver_cell || driver_cell == cell) continue;

            if (driver_cell->isTerminal() || driver_cell->is_fixed || isSequentialCell(driver_cell)) {
                skipped_fixed++;
                continue;
            }

            total_arcs++;

            auto it = pos.find(driver_cell);
            if (it == pos.end()) {
                skipped_missing++;
                continue;
            }

            int driver_pos = it->second;

            // 逆拓扑序：消费者(load) 应在 生产者(driver) 之前
            if (cur_pos >= driver_pos) {
                violations++;
                // std::cout << "[REVERSE-TOPO VIOLATION #" << violations << "] "
                //           << cell->name << " (pos=" << cur_pos
                //           << ") >= driver "
                //           << driver_cell->name
                //           << " (pos=" << driver_pos << ")\n";
            }
        }
        cur_pos++;
    }

    std::cout << "\n=== Traverse Order Verification ===\n"
              << "  Total cells:         " << traverse_list.size() << "\n"
              << "  Total driver arcs:   " << total_arcs << "\n"
              << "  Violations:          " << violations << "\n"
              << "  Fixed/terminal skip: " << skipped_fixed << "\n"
              << "  Driver not in list:  " << skipped_missing << "\n"
              << "  Result:              "
              << (violations == 0 ? "PASS" : "FAIL")
              << std::endl;
}

void Placement3D::refreshPinTiming(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;
    auto cell_pins = cell->getPins();
    /////test
    //  auto pctDiff = [](double v1, double v2) -> double {
    //     double denom = std::max(std::abs(v1), std::abs(v2));
    //     const double eps = 1e-15;
    //     return (denom < eps) ? 0.0 : std::abs(v1 - v2) / denom * 100.0;
    // };
    // if (cell->name == "_647_"){
    //     logger_->info(TDP, 1933, " sta im");
    // }
    ////test
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        auto timing_info = timing_engine_->getPinTimingInfo(sta_pin);
        if (timing_info.is_valid) {
            pin->setEarlyArrival(timing_info.early_arrival);
            pin->setLateArrival(timing_info.late_arrival);
            pin->setEarlyRequired(timing_info.early_required);
            pin->setLateRequired(timing_info.late_required);
            pin->setEarlySlack(timing_info.early_slack);
            pin->setLateSlack(timing_info.late_slack);
        }
        //test////////////
        // auto timing_info2 = timing_engine_->getPinTimingInfo(sta_pin);
        // double max_diff = 0.0;
        // // max_diff = std::max(max_diff, pctDiff(timing_info.early_arrival, full.early_arrival));
        // // max_diff = std::max(max_diff, pctDiff(timing_info.late_arrival,  full.late_arrival));
        // // max_diff = std::max(max_diff, pctDiff(timing_info2.early_required, full.early_required));
        // // max_diff = std::max(max_diff, pctDiff(bypass.late_required,  full.late_required));
        // max_diff = std::max(max_diff, pctDiff(timing_info.early_slack,   timing_info2.early_slack));
        // max_diff = std::max(max_diff, pctDiff(timing_info.late_slack,    timing_info2.late_slack));
        // if (max_diff > 3.0) {
        //     logger_->info(TDP, 1943, "all sta im");
        // }
        //test////////////
    }
}

void Placement3D::refreshPinTiming2(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;
    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        auto timing_info = timing_engine_->getPinTimingInfoDirect(sta_pin);
        if (timing_info.is_valid) {
            pin->setEarlyArrival(timing_info.early_arrival);
            pin->setLateArrival(timing_info.late_arrival);
            pin->setEarlyRequired(timing_info.early_required);
            pin->setLateRequired(timing_info.late_required);
            pin->setEarlySlack(timing_info.early_slack);
            pin->setLateSlack(timing_info.late_slack);
        }
    }
}

void Placement3D::compareLocalVsFullSTA(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;

    struct PinTimingPair {
        std::string name;
        TimingEngine3D::PinTimingInfo local;
        TimingEngine3D::PinTimingInfo full;
    };
    std::vector<PinTimingPair> results;

    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        PinTimingPair p;
        p.name = pin->getName();
        p.local = timing_engine_->getPinTimingInfoDirect(sta_pin);
        results.push_back(p);
    }

    for (auto& p : results) {
        sta::Pin* sta_pin = network->findPin(p.name.c_str());
        if (!sta_pin) continue;
        p.full = timing_engine_->getPinTimingInfo(sta_pin);
    }

    auto pctDiff = [](double a, double b) -> double {
        double denom = std::max(std::abs(a), std::abs(b));
        if (denom < 1e-18) return 0.0;
        return std::abs(a - b) / denom * 100.0;
    };

    double max_arr_pct = 0.0, max_req_pct = 0.0, max_slk_pct = 0.0;
    std::string worst_pin;

    for (const auto& p : results) {
        if (!p.local.is_valid || !p.full.is_valid) continue;
        double arr_pct = pctDiff(p.local.late_arrival, p.full.late_arrival);
        double req_pct = pctDiff(p.local.late_required, p.full.late_required);
        double slk_pct = pctDiff(p.local.late_slack, p.full.late_slack);
        double slk_abs = std::abs(p.local.late_slack - p.full.late_slack);

        if (slk_pct > max_slk_pct) {
            max_slk_pct = slk_pct;
            worst_pin = p.name;
        }
        max_arr_pct = std::max(max_arr_pct, arr_pct);
        max_req_pct = std::max(max_req_pct, req_pct);

        // if (req_pct > 1.0) {
        if (arr_pct > 5 || req_pct > 5) {
            logger_->info(TDP, 890,
                "  pin={} late_arr: local={:.4e} full={:.4e} diff={:.1f}%"
                "  late_req: local={:.4e} full={:.4e} diff={:.1f}%"
                "  late_slk: local={:.4e} full={:.4e} diff={:.1f}% abs={:.4e}",
                p.name,
                p.local.late_arrival, p.full.late_arrival, arr_pct,
                p.local.late_required, p.full.late_required, req_pct,
                p.local.late_slack, p.full.late_slack, slk_pct, slk_abs);
        }
    }

    // logger_->info(TDP, 891,
    //     "cell={} max_arr_diff={:.2f}% max_req_diff={:.2f}% max_slk_diff={:.2f}% worst_pin={}",
    //     cell->name, max_arr_pct, max_req_pct, max_slk_pct, worst_pin);
}

// ============================================================================
// sta::Pin → Node3D 映射辅助函数
// ============================================================================
static Node3D* staPinToNode3D(sta::dbNetwork* db_network,
                               sta::Pin* sta_pin,
                               Network3D* network_3d)
{
    odb::dbITerm* iterm = nullptr;
    odb::dbBTerm* bterm = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    db_network->staToDb(sta_pin, iterm, bterm, moditerm);
    
    if (iterm) {
        odb::dbInst* inst = iterm->getInst();
        if (inst) return network_3d->getCell(inst);
    }
    if (bterm) {
        return network_3d->getNode(bterm);
    }
    return nullptr;
}

// ============================================================================
// 构建关键路径优先遍历顺序
//   取 top-k 违反约束路径，收集沿途 cell（排除终端 PI/PO），
//   路径回溯顺序即为逆拓扑序（PO→PI），剩余 cell 保持原顺序追加。
// ============================================================================
std::vector<Node3D*> LSP3D::buildCriticalPathPriorityOrder(
    TimingEngine3D* timing_engine,
    double top_k)
{
    std::vector<Node3D*> priority;
    std::unordered_set<Node3D*> in_priority;

    if (!timing_engine || !network_) return priority;

    sta::dbSta* sta = timing_engine->getSta();
    sta::dbNetwork* db_network = sta->getDbNetwork();

    auto paths = timing_engine->getViolatingPaths();

    int path_count = 0;
    for (sta::Path* p : paths) {
        if (path_count >= top_k) break;

        sta::Path* step = p;
        while (step) {
            sta::Pin* sta_pin = step->pin(sta);
            if (sta_pin) {
                Node3D* cell = staPinToNode3D(db_network, sta_pin, network_);
                if (cell
                    && !cell->isTerminal()
                    && !cell->is_fixed
                    && !in_priority.count(cell)) {
                    in_priority.insert(cell);
                    priority.push_back(cell);
                }
            }
            step = step->prevPath();
        }
        path_count++;
    }

    // for (auto* cell : traverse_list) {
    //     if (cell
    //         && !cell->is_fixed
    //         && !in_priority.count(cell)) {
    //         priority.push_back(cell);
    //     }
    // }

    // S1: how much of total Late TNS is covered by the selected violating paths.
    // low coverage => cell selection misses TNS mass (env-gated, log only).
    if (std::getenv("TDP_DIAG_S1")) {
        double covered_tns = 0.0;
        int nep = 0;
        for (sta::Path* p : paths) {
            sta::Pin* ep = p ? p->pin(sta) : nullptr;
            if (!ep) continue;
            double s = timing_engine->getPinSlack(ep, /*is_early=*/false);
            if (s < 0.0) { covered_tns += s; nep++; }
        }
        double total_tns = timing_engine->getTimingStats().late_tns;
        logger_->info(utl::TDP, 2074,
            "[S1] violating_paths={} covered_ep={} covered_TNS={:.6f}ns "
            "total_TNS={:.6f}ns coverage={:.1f}%",
            (int) paths.size(), nep, covered_tns * 1e9, total_tns * 1e9,
            total_tns != 0.0 ? 100.0 * covered_tns / total_tns : 0.0);
    }

    return priority;
}

} // namespace tdp
