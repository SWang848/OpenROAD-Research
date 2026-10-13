// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "placement.h"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "../eva/eva.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Pin3D.h"
#include "../infrastructure/Common3D.h"

#include "PtGraph.hh"
#include "lgista/LocalGraphIsta.hh"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "sta/Parasitics.hh"
#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/TimingRole.hh"
#include "sta/TimingArc.hh"
#include "sta/Path.hh"
#include "search/Tag.hh"
#include "search/TagGroup.hh"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <cmath>
#include <limits> 
#include <algorithm>
#include <memory>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <vector>
#include <omp.h>

using utl::TDP;

namespace tdp {

// forward decl: defined later in this file
static bool ptVertexLocalArrival(lgista::PtGraph* pg, lgista::PtVertex& ptv,
                                 sta::StaState* sta_state,
                                 double& early_out, double& late_out);

// ---- P0 diagnostic: correlation between local-eval and global-STA costs ----
// Guards: enabled only when TDP_SPIKE_CORR is set in the environment.
// Accumulated per DLS iteration, dumped after each iteration's commit loop,
// and reset to start fresh for the next iteration.
// P0-A enhancement: also tracks baseline_gap (|local_cur − global_cur|) and
// moved_gap (|local_best − global_after|) to decompose whether errors come from
// serial-commit pollution (baseline dominated) or override-vs-odb mismatch
// (moved dominated).  Separate stats for all-accepted and opp_sign only.
struct CorrStats {
    int count = 0;
    double sum_delta_local  = 0.0;
    double sum_delta_global = 0.0;
    double sum_sq_local     = 0.0;
    double sum_sq_global    = 0.0;
    double sum_cross        = 0.0;
    double sum_baseline_gap = 0.0;
    double sum_moved_gap    = 0.0;
    int    same_sign = 0;
    int    opp_sign  = 0;
    int    local_improve_global_worse = 0;
    int    local_worse_global_improve = 0;

    int    opp_count = 0;
    double opp_sum_baseline = 0.0;
    double opp_sum_moved    = 0.0;
    int    opp_printed = 0;

    double sum_current_cost = 0.0;
    double min_current_cost = 1e100;
    double max_current_cost = 0.0;

    void record(double dlocal, double dglobal, double baseline_gap, double moved_gap, double cur_cost) {
        ++count;
        sum_delta_local  += dlocal;
        sum_delta_global += dglobal;
        sum_sq_local     += dlocal  * dlocal;
        sum_sq_global    += dglobal * dglobal;
        sum_cross        += dlocal  * dglobal;
        sum_baseline_gap += baseline_gap;
        sum_moved_gap    += moved_gap;
        if (dlocal > 0 && dglobal > 0) ++same_sign;
        else if (dlocal < 0 && dglobal < 0) ++same_sign;
        else if (dlocal != 0.0 || dglobal != 0.0) ++opp_sign;
        if (dlocal > 0 && dglobal < 0) ++local_improve_global_worse;
        if (dlocal < 0 && dglobal > 0) ++local_worse_global_improve;

        sum_current_cost += cur_cost;
        if (cur_cost < min_current_cost) min_current_cost = cur_cost;
        if (cur_cost > max_current_cost) max_current_cost = cur_cost;
    }

    void record_opp(double baseline_gap, double moved_gap) {
        ++opp_count;
        opp_sum_baseline += baseline_gap;
        opp_sum_moved    += moved_gap;
    }

    double pearsonR() const {
        if (count < 2) return 0.0;
        double num = count * sum_cross - sum_delta_local * sum_delta_global;
        double den1 = count * sum_sq_local  - sum_delta_local  * sum_delta_local;
        double den2 = count * sum_sq_global - sum_delta_global * sum_delta_global;
        double den = std::sqrt(den1 * den2);
        return den > 1e-30 ? num / den : 0.0;
    }
    bool is_active() const { return count > 0; }
};

// ---- TDP_DIAG_CAND candidate distribution diagnostic -----------------------
struct CandStats {
    int    cells_evaluated    = 0;
    int    cells_found_better = 0;
    double sum_improve_ratio  = 0.0;   // Σ (current - best) / current
    double max_improve_ratio  = 0.0;
    double sum_n_better       = 0.0;   // Σ n_better per cell
    double sum_pct_better     = 0.0;   // Σ (n_better / N) * 100
    double sum_cost_range_ratio = 0.0; // Σ (max-min) / current
    double max_cost_range_ratio = 0.0;
    double sum_current_cost   = 0.0;   // for avg reference
    void reset() { *this = CandStats(); }
};
static CandStats g_cand;

static bool diagCandEnabled() {
    const char* v = std::getenv("TDP_DIAG_CAND");
    return v != nullptr;
}

static void dumpCandStats(utl::Logger* logger) {
    if (!diagCandEnabled() || g_cand.cells_evaluated == 0) return;
    static int cand_iter = 0;
    double N = static_cast<double>(g_cand.cells_evaluated);
    double M = static_cast<double>(g_cand.cells_found_better);
    logger->info(utl::TDP, 2068,
        "[CAND] Iter {}: {} cells evaluated, {} found better ({:.1f}%)",
        cand_iter, g_cand.cells_evaluated, g_cand.cells_found_better, 100.0 * M / N);
    logger->info(utl::TDP, 2069,
        "[CAND]   avg_improve_ratio={:.4f}  max_improve_ratio={:.4f}  avg_current_cost={:.6e}",
        M > 0 ? g_cand.sum_improve_ratio / M : 0.0,
        g_cand.max_improve_ratio,
        N > 0 ? g_cand.sum_current_cost / N : 0.0);
    logger->info(utl::TDP, 2080,
        "[CAND]   avg_n_better={:.1f}/100  avg_pct_better={:.1f}%  avg_range_ratio={:.4f}  max_range_ratio={:.4f}",
        N > 0 ? g_cand.sum_n_better / N : 0.0,
        N > 0 ? g_cand.sum_pct_better / N : 0.0,
        N > 0 ? g_cand.sum_cost_range_ratio / N : 0.0,
        g_cand.max_cost_range_ratio);
    g_cand.reset();
    ++cand_iter;
}

// ---- TDP_DIAG_W2 weight distribution diagnostic ----------------------------
struct W2Stats {
    int    n_pins        = 0;
    int    n_crit_late   = 0;    // late_arrival > late_required
    int    n_crit_early  = 0;    // early_arrival < early_required
    double sum_lw_all    = 0.0;  // Σ late_weight
    double sum_lw_crit   = 0.0;  // Σ late_weight on critical pins
    double sum_lw_ncrit  = 0.0;  // Σ late_weight on non-critical pins
    double min_lw        = 1e100;
    double max_lw        = 0.0;
    double sum_ew_all    = 0.0;  // Σ early_weight
    double sum_ew_crit   = 0.0;
    double sum_ew_ncrit  = 0.0;
    double min_ew        = 1e100;
    double max_ew        = 0.0;
    void reset() { *this = W2Stats(); }
};
static W2Stats g_w2;

static bool diagW2Enabled() {
    const char* v = std::getenv("TDP_DIAG_W2");
    return v != nullptr;
}

void LSP3D::collectW2Stats() {
    if (!diagW2Enabled()) return;
    g_w2.reset();
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        for (auto* pin : cell->getPins()) {
            if (!pin) continue;
            double lw = pin->getLateWeight();
            double ew = pin->getEarlyWeight();
            double la = pin->getLateArrival();
            double lr = pin->getLateRequired();
            double ea = pin->getEarlyArrival();
            double er = pin->getEarlyRequired();

            ++g_w2.n_pins;
            g_w2.sum_lw_all += lw;
            g_w2.sum_ew_all += ew;
            if (lw < g_w2.min_lw) g_w2.min_lw = lw;
            if (lw > g_w2.max_lw) g_w2.max_lw = lw;
            if (ew < g_w2.min_ew) g_w2.min_ew = ew;
            if (ew > g_w2.max_ew) g_w2.max_ew = ew;

            bool crit_l = (lr != 0.0 && la > lr);
            if (crit_l) {
                ++g_w2.n_crit_late;
                g_w2.sum_lw_crit += lw;
            } else {
                g_w2.sum_lw_ncrit += lw;
            }
            bool crit_e = (ea != 0.0 && er != 0.0 && ea < er);
            if (crit_e) {
                ++g_w2.n_crit_early;
                g_w2.sum_ew_crit += ew;
            } else {
                g_w2.sum_ew_ncrit += ew;
            }
        }
    }
}

static void dumpW2Stats(utl::Logger* logger) {
    if (!diagW2Enabled() || g_w2.n_pins == 0) return;
    static int w2_iter = 0;
    double crit_mean_lw = g_w2.n_crit_late > 0
        ? g_w2.sum_lw_crit / g_w2.n_crit_late : 0.0;
    double ncrit_mean_lw = (g_w2.n_pins - g_w2.n_crit_late) > 0
        ? g_w2.sum_lw_ncrit / (g_w2.n_pins - g_w2.n_crit_late) : 0.0;
    double ratio_lw = ncrit_mean_lw > 0 ? crit_mean_lw / ncrit_mean_lw : 0.0;
    double crit_mean_ew = g_w2.n_crit_early > 0
        ? g_w2.sum_ew_crit / g_w2.n_crit_early : 0.0;
    double ncrit_mean_ew = (g_w2.n_pins - g_w2.n_crit_early) > 0
        ? g_w2.sum_ew_ncrit / (g_w2.n_pins - g_w2.n_crit_early) : 0.0;
    double ratio_ew = ncrit_mean_ew > 0 ? crit_mean_ew / ncrit_mean_ew : 0.0;
    logger->info(utl::TDP, 2081,
        "[W2] Iter {}: n_pins={}  n_crit_late={} ({:.1f}%)  n_crit_early={} ({:.1f}%)",
        w2_iter, g_w2.n_pins,
        g_w2.n_crit_late, 100.0 * g_w2.n_crit_late / g_w2.n_pins,
        g_w2.n_crit_early, 100.0 * g_w2.n_crit_early / g_w2.n_pins);
    logger->info(utl::TDP, 2082,
        "[W2]   late_weight:  min={:.4f} max={:.4f} mean={:.4f}  crit_mean={:.4f} ncrit_mean={:.4f} ratio={:.4f}",
        g_w2.min_lw, g_w2.max_lw,
        g_w2.n_pins > 0 ? g_w2.sum_lw_all / g_w2.n_pins : 0.0,
        crit_mean_lw, ncrit_mean_lw, ratio_lw);
    logger->info(utl::TDP, 2083,
        "[W2]   early_weight: min={:.4f} max={:.4f} mean={:.4f}  crit_mean={:.4f} ncrit_mean={:.4f} ratio={:.4f}",
        g_w2.min_ew, g_w2.max_ew,
        g_w2.n_pins > 0 ? g_w2.sum_ew_all / g_w2.n_pins : 0.0,
        crit_mean_ew, ncrit_mean_ew, ratio_ew);
    ++w2_iter;
}

// ---- TDP_DIAG_W3 critical path weight update tracking ----------------------
struct W3Stats {
    int n_crit_cells     = 0;
    int n_seq = 0, n_comb = 0;
    int n_updated_seq = 0, n_unchanged_seq = 0;
    int n_updated_comb = 0, n_unchanged_comb = 0;
    double max_lw_delta_seq = 0, max_lw_delta_comb = 0;
    double sum_lw_delta_seq = 0, sum_lw_delta_comb = 0;
    void reset() { *this = W3Stats(); }
};
static W3Stats g_w3;

static bool diagW3Enabled() {
    const char* v = std::getenv("TDP_DIAG_W3");
    return v != nullptr;
}

void LSP3D::collectW3Stats(
    const std::unordered_set<Node3D*>& crit_cells,
    const std::unordered_map<Pin3D*, std::pair<double,double>>& snap)
{
    if (!diagW3Enabled()) return;
    g_w3.reset();
    for (auto* cell : crit_cells) {
        if (!cell) continue;
        ++g_w3.n_crit_cells;
        bool is_seq = isSequentialCell(cell);
        bool any_changed = false;
        double max_delta = 0;
        for (auto* pin : cell->getPins()) {
            if (!pin) continue;
            auto it = snap.find(pin);
            if (it == snap.end()) continue;
            double old_lw = it->second.first;
            double cur_lw = pin->getLateWeight();
            double delta = std::abs(cur_lw - old_lw);
            if (delta > 1e-15) {
                any_changed = true;
                max_delta = std::max(max_delta, delta);
            }
        }
        if (is_seq) {
            ++g_w3.n_seq;
            if (any_changed) {
                ++g_w3.n_updated_seq;
                g_w3.sum_lw_delta_seq += max_delta;
                g_w3.max_lw_delta_seq = std::max(g_w3.max_lw_delta_seq, max_delta);
            } else {
                ++g_w3.n_unchanged_seq;
            }
        } else {
            ++g_w3.n_comb;
            if (any_changed) {
                ++g_w3.n_updated_comb;
                g_w3.sum_lw_delta_comb += max_delta;
                g_w3.max_lw_delta_comb = std::max(g_w3.max_lw_delta_comb, max_delta);
            } else {
                ++g_w3.n_unchanged_comb;
            }
        }
    }
}

static void dumpW3Stats(utl::Logger* logger) {
    if (!diagW3Enabled() || g_w3.n_crit_cells == 0) return;
    static int w3_iter = 0;
    logger->info(utl::TDP, 2084,
        "[W3] Iter {}: critical cells={}  (seq={}, comb={})",
        w3_iter, g_w3.n_crit_cells, g_w3.n_seq, g_w3.n_comb);
    double avg_seq = g_w3.n_updated_seq > 0
        ? g_w3.sum_lw_delta_seq / g_w3.n_updated_seq : 0.0;
    logger->info(utl::TDP, 2085,
        "[W3]   sequential: updated={}/{} unchanged={}  avg_lw_delta={:.4f} max_lw_delta={:.4f}",
        g_w3.n_updated_seq, g_w3.n_seq, g_w3.n_unchanged_seq, avg_seq, g_w3.max_lw_delta_seq);
    double avg_comb = g_w3.n_updated_comb > 0
        ? g_w3.sum_lw_delta_comb / g_w3.n_updated_comb : 0.0;
    logger->info(utl::TDP, 2086,
        "[W3]   combinational: updated={}/{} unchanged={}  avg_lw_delta={:.4f} max_lw_delta={:.4f}",
        g_w3.n_updated_comb, g_w3.n_comb, g_w3.n_unchanged_comb, avg_comb, g_w3.max_lw_delta_comb);
    ++w3_iter;
}

static CorrStats g_corr;        // accumulates during one DLS iteration
static int g_corr_iter = 0;     // iteration counter for the dump header

static bool corr_enabled() {
    static int enabled = -1;
    if (enabled < 0) enabled = (std::getenv("TDP_SPIKE_CORR") != nullptr) ? 1 : 0;
    return enabled == 1;
}

// ---- Phase profiler for par mode (TDP_PROF=1) ----------------------------
// Splits per-cell discreteLocalSearchParallel time into serial phases
// (candgen / prime / baseline / commit) and the parallel candidate-eval phase
// (wall time + per-thread CPU sums of its sub-steps). Reveals whether runtime
// is serial-bound (prime+baseline+commit) or parallel-bound (parasitic/analyze).
struct ProfStats {
    double candgen = 0, prime = 0, baseline = 0, commit = 0, par_wall = 0, misc = 0;
    double par_parasitic = 0, par_analyze = 0, par_maps = 0, par_cost = 0;
    // analyzeOnly sub-step granularity (cpu accumulated over threads)
    double par_ana_makept   = 0, par_ana_recomp  = 0,
           par_ana_delays   = 0, par_ana_arrivals = 0,
           par_ana_requireds = 0;
    // commit sub-steps
    double commit_steiner = 0, commit_exec = 0, commit_refresh = 0;
    long   cells = 0, commits = 0, candidates = 0;
};
static ProfStats g_prof;
static int g_prof_iter = 0;

static bool prof_enabled() {
    static int enabled = -1;
    if (enabled < 0) enabled = (std::getenv("TDP_PROF") != nullptr) ? 1 : 0;
    return enabled == 1;
}

static inline double now_ms() {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---- Lightweight diagnostics: W1(weight dist)/S1(coverage)/P1(gate-vs-wire)/
//      T1(in-sweep TNS). All env-gated (TDP_DIAG_*), default OFF, emit logs only
//      and do NOT change placement decisions. NOTE: T1 calls getTimingStats()
//      mid-sweep which forces an STA update, so enabling T1 perturbs the run.
static void diagLogStats(utl::Logger* logger, const char* label,
                         std::vector<double>& v) {
    if (v.empty()) { logger->info(utl::TDP, 2070, "[{}] empty", label); return; }
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    double sum = 0.0, nzsum = 0.0;
    size_t zeros = 0, nz = 0;
    for (double x : v) {
        sum += x;
        if (std::abs(x) < 1e-12) zeros++; else { nz++; nzsum += x; }
    }
    const double med = v[n / 2];
    const double p99 = v[std::min(n - 1, (size_t)std::floor(0.99 * n))];
    logger->info(utl::TDP, 2071,
        "[{}] n={} zero={}({:.1f}%) min={:.3e} med={:.3e} mean={:.3e} "
        "nz_mean={:.3e} p99={:.3e} max={:.3e}",
        label, n, zeros, 100.0 * zeros / n, v.front(), med, sum / n,
        nz ? nzsum / nz : 0.0, p99, v.back());
}

// P1: decompose the violating paths into gate(cell)-arc vs net(wire)-arc delay,
// classifying each arc canonically via Edge::isWire(). small wire_frac =>
// single-cell placement has little TNS headroom (gate-bound). By default all
// violating paths are summed; env TDP_DIAG_P1_N=<N> limits to the worst N.
static void diagPathDecomp(TimingEngine3D* te, utl::Logger* logger) {
    if (!te) return;
    sta::dbSta* sta = te->getSta();
    auto paths = te->getViolatingPaths();
    int topN = static_cast<int>(paths.size());
    if (const char* v = std::getenv("TDP_DIAG_P1_N")) {
        int n = std::atoi(v);
        if (n > 0) topN = std::min(topN, n);
    }
    double tot_gate = 0.0, tot_wire = 0.0;
    int np = 0;
    for (sta::Path* p : paths) {
        if (np >= topN) break;
        sta::Path* step = p;
        sta::Path* prev = step ? step->prevPath() : nullptr;
        while (step && prev) {
            const double d = sta::delayAsFloat(step->arrival())
                           - sta::delayAsFloat(prev->arrival());
            // Classify the arc feeding 'step' via Edge::isWire() (canonical);
            // fall back to the timing arc's role if the edge is unavailable.
            bool is_wire;
            sta::Edge* e = step->prevEdge(sta);
            if (e) {
                is_wire = e->isWire();
            } else {
                sta::TimingArc* a = step->prevArc(sta);
                is_wire = (a != nullptr) && a->role() && a->role()->isWire();
            }
            if (is_wire) tot_wire += d; else tot_gate += d;
            step = prev;
            prev = step->prevPath();
        }
        np++;
    }
    const double tot = tot_gate + tot_wire;
    logger->info(utl::TDP, 2072,
        "[P1] {}/{} violating paths: gate={:.6f}ns wire={:.6f}ns wire_frac={:.1f}%",
        np, static_cast<int>(paths.size()), tot_gate * 1e9, tot_wire * 1e9,
        tot != 0.0 ? 100.0 * tot_wire / tot : 0.0);
}

static void dumpProfStats(utl::Logger* logger) {
    if (g_prof.cells == 0) return;
    double serial = g_prof.candgen + g_prof.prime + g_prof.baseline + g_prof.commit + g_prof.misc;
    double par_cpu = g_prof.par_parasitic + g_prof.par_analyze + g_prof.par_maps + g_prof.par_cost;
    double wall = serial + g_prof.par_wall;
    auto pct = [&](double v) { return wall > 0 ? 100.0 * v / wall : 0.0; };
    logger->info(utl::TDP, 2060,
        "[PROF] Iter {} | cells={} commits={} candidates={} | wall(sum)={:.0f}ms",
        g_prof_iter, g_prof.cells, g_prof.commits, g_prof.candidates, wall);
    logger->info(utl::TDP, 2061,
        "[PROF]  SERIAL {:.0f}ms ({:.1f}%): candgen={:.0f} prime={:.0f} baseline={:.0f} commit={:.0f} misc={:.0f}",
        serial, pct(serial), g_prof.candgen, g_prof.prime, g_prof.baseline, g_prof.commit, g_prof.misc);
    logger->info(utl::TDP, 2062,
        "[PROF]  PARALLEL wall={:.0f}ms ({:.1f}%)  cpu_sum={:.0f}ms  effective_parallelism={:.1f}x",
        g_prof.par_wall, pct(g_prof.par_wall), par_cpu,
        g_prof.par_wall > 0 ? par_cpu / g_prof.par_wall : 0.0);
    logger->info(utl::TDP, 2063,
        "[PROF]  PARALLEL cpu breakdown: parasitic={:.0f}ms analyze={:.0f}ms maps={:.0f}ms cost={:.0f}ms",
        g_prof.par_parasitic, g_prof.par_analyze, g_prof.par_maps, g_prof.par_cost);
    auto amake = g_prof.par_ana_makept, arecomp = g_prof.par_ana_recomp,
         adelays = g_prof.par_ana_delays, aarr = g_prof.par_ana_arrivals,
         areq  = g_prof.par_ana_requireds;
    logger->info(utl::TDP, 2064,
        "[PROF]    analyzeOnly sub: makePtGraph={:.0f} recomputePara={:.0f} delays={:.0f} arrivals={:.0f} requireds={:.0f}",
        amake, arecomp, adelays, aarr, areq);
    logger->info(utl::TDP, 2065,
        "[PROF]    commit sub: steiner={:.0f} executeTiming={:.0f} refresh={:.0f}",
        g_prof.commit_steiner, g_prof.commit_exec, g_prof.commit_refresh);
    g_prof = ProfStats();
    ++g_prof_iter;
}

static void dumpCorrStats(utl::Logger* logger) {
    if (!g_corr.is_active()) return;
    logger->info(utl::TDP, 2050,
        "[SPIKE_CORR] Iter {} summary: {} accepted moves", g_corr_iter, g_corr.count);
    logger->info(utl::TDP, 2051,
        "[SPIKE_CORR]   Pearson r = {:.4f}", g_corr.pearsonR());
    logger->info(utl::TDP, 2052,
        "[SPIKE_CORR]   same_sign={} opp_sign={}  local-improve-but-global-worse={}  local-worse-but-global-improve={}",
        g_corr.same_sign, g_corr.opp_sign,
        g_corr.local_improve_global_worse, g_corr.local_worse_global_improve);
    logger->info(utl::TDP, 2053,
        "[SPIKE_CORR]   ALL  {} moves | mean|baseline_gap|={:.6e}  mean|moved_gap|={:.6e}",
        g_corr.count,
        g_corr.count ? g_corr.sum_baseline_gap / g_corr.count : 0.0,
        g_corr.count ? g_corr.sum_moved_gap    / g_corr.count : 0.0);
    if (g_corr.count > 0) {
        double avg_cost       = g_corr.sum_current_cost / g_corr.count;
        double avg_delta      = g_corr.sum_delta_local  / g_corr.count;
        double ratio          = avg_cost > 0.0 ? avg_delta / avg_cost : 0.0;
        logger->info(utl::TDP, 2066,
            "[SPIKE_CORR]   avg_current_cost={:.6e}  min_cost={:.6e}  max_cost={:.6e}",
            avg_cost, g_corr.min_current_cost, g_corr.max_current_cost);
        logger->info(utl::TDP, 2067,
            "[SPIKE_CORR]   avg_delta_local={:.6e}  ratio=delta/cost={:.4f}  (signal strength)",
            avg_delta, ratio);
    }
    if (g_corr.opp_count > 0)
        logger->info(utl::TDP, 2054,
            "[SPIKE_CORR]   OPP  {} moves | mean|baseline_gap|={:.6e}  mean|moved_gap|={:.6e}",
            g_corr.opp_count,
            g_corr.opp_sum_baseline / g_corr.opp_count,
            g_corr.opp_sum_moved    / g_corr.opp_count);
    g_corr = CorrStats();
    ++g_corr_iter;
}

LSP3D::LSP3D(Timingdp* timingdp, const LRParams& params)
    : Placement3D(timingdp), params_(params)
{
    if (timingdp_) {
        arch_ = timingdp_->getArchitecture3D();
        network_ = timingdp_->getNetwork3D();
    }
}

LSP3D::~LSP3D()
{
    // 清理资源
}

// ============================================================================
// 主要布局方法
// ============================================================================

void LSP3D::place()
{
    if (!timingdp_ || !arch_ || !network_) {
        logger_->error(TDP, 1305, "Error: Database, architecture, or network not initialized");
        return;
    }

    if (std::getenv("TDP_SPIKE1")) {
        runSearchConcurrencySpike();
        return;
    }
    if (std::getenv("TDP_SPIKE2")) {
        runTrialParasiticSpike();
        return;
    }
    if (std::getenv("TDP_SPIKE3")) {
        runConcurrentTrialSpike();
        return;
    }
    if (std::getenv("TDP_SPIKE2B")) {
        runTrialMovedSpike();
        return;
    }
    if (std::getenv("TDP_SPIKE_EXT")) {
        runExtractorCompareSpike();
        return;
    }
    if (std::getenv("TDP_SPIKE_COST")) {
        runCostAttributionSpike();
        return;
    }

    initWeight(timingdp_,timing_engine_);
    // 获取遍历顺序（反向拓扑顺序）
    logger_ -> info(TDP, 1934, "Starting discrete local search...");
    
    double L_lambda = 0.0;
    for (int iter = 0; iter < params_.max_iterations; ++iter) {
        auto ti0 = std::chrono::steady_clock::now();


        auto stats = timing_engine_->getTimingStats();
        logger_->info(TDP, 1294, "Iter {}: Early WNS={:.6f}ns TNS={:.6f}ns | Late WNS={:.6f}ns TNS={:.6f}ns",
            iter,
            stats.early_wns * 1e9, stats.early_tns * 1e9,
            stats.late_wns * 1e9, stats.late_tns * 1e9);
        // ---- P2 inline: gate vs wire delay decomposition on worst path (TDP_DIAG_P2=N) ----
        if (const char* p2env = std::getenv("TDP_DIAG_P2")) {
            int topN = std::atoi(p2env); if (topN <= 0) topN = 5;
            auto paths = timing_engine_->getViolatingPaths();
            int np = std::min(topN, (int)paths.size());
            sta::dbNetwork* nw = sta_->getDbNetwork();
            for (int pi = 0; pi < np && pi < 1; ++pi) {   // 只打第 1 条（最差）
                sta::Path* pp = paths[pi]; if (!pp) continue;
                std::vector<sta::Path*> chain;
                for (sta::Path* s = pp; s; s = s->prevPath()) chain.push_back(s);
                std::reverse(chain.begin(), chain.end());
                double gsum = 0.0, wsum = 0.0;
                for (size_t k = 1; k < chain.size(); ++k) {
                    double d = sta::delayAsFloat(chain[k]->arrival()) - sta::delayAsFloat(chain[k-1]->arrival());
                    sta::Edge* e = chain[k]->prevEdge(sta_);
                    if (e && e->isWire()) wsum += d; else gsum += d;
                }
                double total = gsum + wsum;
                logger_->info(TDP, 2087,
                    "[P2] Iter {} worst path: slack={:.6f}ns gate={:.6f}ns wire={:.6f}ns wire_frac={:.1f}%",
                    iter, sta::delayAsFloat(pp->slack(sta_))*1e9, gsum*1e9, wsum*1e9,
                    total > 0 ? 100.0*wsum/total : 0.0);
            }
        }
        // ---- TDP_DIAG_W3: snapshot critical path cell weights ----
        std::unordered_set<Node3D*> crit_cells;
        std::unordered_map<Pin3D*, std::pair<double,double>> snap_weights;
        if (diagW3Enabled()) {
            auto prio = buildCriticalPathPriorityOrder(timing_engine_, 1000);
            for (auto* c : prio) if (c) crit_cells.insert(c);
            for (auto* c : crit_cells)
                for (auto* p : c->getPins())
                    if (p) snap_weights[p] = {p->getLateWeight(), p->getEarlyWeight()};
        }
        updatePinWeights();
        logger_->info(TDP, 1292, " update pin weight over" );
        // ---- weight overflow scan ----
        {
            int ovf_cnt = 0;
            double max_lw = 0.0;
            std::string worst_pin;
            for (auto* cell : network_->getCells()) {
                if (!cell) continue;
                for (auto* pin : cell->getPins()) {
                    if (!pin) continue;
                    double lw = pin->getLateWeight();
                    if (std::isinf(lw) || std::abs(lw) > 1e10) {
                        ++ovf_cnt;
                        if (std::abs(lw) > max_lw) {
                            max_lw = std::abs(lw);
                            worst_pin = cell->name + "/" + pin->getName();
                        }
                    }
                }
            }
            if (ovf_cnt > 0)
                logger_->warn(TDP, 2094, "[WGT-OVF] {} pins overflow, max={:.6e} worst_pin={}",
                    ovf_cnt, max_lw, worst_pin);
        }
        collectW2Stats();
        if (diagW2Enabled()) dumpW2Stats(logger_);
        if (diagW3Enabled()) { collectW3Stats(crit_cells, snap_weights); dumpW3Stats(logger_); }

        auto ti1 = std::chrono::steady_clock::now();

        discreteLocalSearchAll(timing_engine_);

        auto ti2 = std::chrono::steady_clock::now();

        SteinerTreeGenerator stn_gen(timingdp_);
        stn_gen.createSteninerRoute();

        auto ti3 = std::chrono::steady_clock::now();

        sta_->updateTiming(false);
        sta_->findRequireds();

        auto ti4 = std::chrono::steady_clock::now();

#define DUR_MS(a,b) std::chrono::duration<double>(b-a).count()
        logger_->info(TDP, 3003,
            "[ITER] {} weights={:.1f}s dls={:.1f}s steiner={:.1f}s sta={:.1f}s total={:.1f}s",
            iter, DUR_MS(ti0,ti1), DUR_MS(ti1,ti2), DUR_MS(ti2,ti3), DUR_MS(ti3,ti4), DUR_MS(ti0,ti4));
#undef DUR_MS

        // ---- per-iteration cost checksum for optimization verification ----
        {
            double total_lambda = 0.0;
            for (auto* c : network_->getCells()) {
                if (c) total_lambda += calculateCellCost(c, 5.0, 1.0);
            }
            logger_->info(TDP, 3004, "[CHKSUM] Iter {} total_cost={:.15e}", iter, total_lambda);
        }

            
    }
    
}

// ============================================================================
// Spike 1: sta::Search 并发只读安全性验证
// 多线程各用一套 per-thread LocalGraphIsta，并发对同一冻结全局图跑 analyzeOnly
// （不 commit / 不移动 cell / 不 updateSteinerTrees），对比并行结果是否与串行
// 位级一致。若一致且无崩溃 -> 局部 STA 计算内核可并发只读。
// 通过 TDP_SPIKE1=1 启用；TDP_SPIKE1_THREADS / _CELLS / _REPEAT / _DEPTH 可调。
// ============================================================================
void LSP3D::runSearchConcurrencySpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int num_threads = env_int("TDP_SPIKE1_THREADS", 8);
    const int max_cells   = env_int("TDP_SPIKE1_CELLS", 128);
    const int repeat      = env_int("TDP_SPIKE1_REPEAT", 5);
    const int depth       = env_int("TDP_SPIKE1_DEPTH", 2);

    logger_->info(TDP, 1960,
        "[SPIKE1] search-concurrency: threads={} cells<={} repeat={} depth={}",
        num_threads, max_cells, repeat, depth);

    // 冻结一份干净的全局时序图
    sta_->updateTiming(true);
    sta_->findRequireds();

    sta::dbNetwork* db_network = sta_->getDbNetwork();

    // 串行解析 cell -> sta::Instance（避免并行区内的 dbToSta 查找）
    std::vector<sta::Instance*> insts;
    insts.reserve(max_cells);
    for (auto* cell : network_->getCells()) {
        if (static_cast<int>(insts.size()) >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst)
            continue;
        sta::Instance* si = db_network->dbToSta(cell->db_inst);
        if (si) insts.push_back(si);
    }
    const int K = static_cast<int>(insts.size());
    if (K == 0) {
        logger_->warn(TDP, 1961, "[SPIKE1] no movable cells collected");
        return;
    }

    // 每线程一套独立引擎
    std::vector<std::unique_ptr<lgista::LocalGraphIsta>> pool;
    pool.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t)
        pool.push_back(std::make_unique<lgista::LocalGraphIsta>(sta_));

    // 一个 PtGraph 的位稳定指纹：所有真实顶点 (early+late) 之和 + 计数
    auto fingerprint = [&](lgista::PtGraph* pg) -> std::pair<double,int> {
        if (!pg) return {0.0, -1};
        double acc = 0.0; int n = 0;
        for (auto& ptv : pg->ptVertices()) {
            if (!ptv.vertex() || !ptv.paths()) continue;
            double e = 0.0, l = 0.0;
            if (ptVertexLocalArrival(pg, ptv, sta_, e, l)) {
                acc += e + l; ++n;
            }
        }
        return {acc, n};
    };

    // 串行预热：ensureInit/copyState 必须串行完成，避免懒初始化竞争
    for (int t = 0; t < num_threads; ++t) {
        auto pg = pool[t]->analyzeOnly(insts[0], depth);
        (void)pg;
    }

    // 串行参考
    std::vector<std::pair<double,int>> ref(K);
    auto ts0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < K; ++i) {
        auto pg = pool[0]->analyzeOnly(insts[i], depth);
        ref[i] = fingerprint(pg.get());
    }
    auto ts1 = std::chrono::high_resolution_clock::now();
    const double serial_ms =
        std::chrono::duration<double, std::milli>(ts1 - ts0).count();

    // 并行多次
    int total_mismatch = 0;
    double worst_delta = 0.0;
    double par_ms_best = 1e300;
    for (int r = 0; r < repeat; ++r) {
        std::vector<std::pair<double,int>> par(K, {0.0, -2});
        auto tp0 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (int i = 0; i < K; ++i) {
            int tid = omp_get_thread_num();
            auto pg = pool[tid]->analyzeOnly(insts[i], depth);
            par[i] = fingerprint(pg.get());
        }
        auto tp1 = std::chrono::high_resolution_clock::now();
        const double par_ms =
            std::chrono::duration<double, std::milli>(tp1 - tp0).count();
        par_ms_best = std::min(par_ms_best, par_ms);

        int mismatch = 0;
        for (int i = 0; i < K; ++i) {
            if (par[i].second != ref[i].second || par[i].first != ref[i].first) {
                ++mismatch;
                worst_delta = std::max(worst_delta,
                                       std::abs(par[i].first - ref[i].first));
            }
        }
        total_mismatch += mismatch;
        logger_->info(TDP, 1962, "[SPIKE1] run {}: mismatch {}/{}  par_ms={:.1f}",
                      r, mismatch, K, par_ms);
    }

    logger_->info(TDP, 1963,
        "[SPIKE1] DONE cells={} threads={} serial_ms={:.1f} par_ms(best)={:.1f} "
        "speedup={:.2f}x total_mismatch={} worst_delta={:.3e}",
        K, num_threads, serial_ms, par_ms_best,
        (par_ms_best > 0 ? serial_ms / par_ms_best : 0.0),
        total_mismatch, worst_delta);

    if (total_mismatch == 0)
        logger_->info(TDP, 1964,
            "[SPIKE1] PASS: parallel == serial (bit-identical) over {} runs", repeat);
    else
        logger_->warn(TDP, 1965,
            "[SPIKE1] FAIL: {} mismatches -> concurrent read NOT safe", total_mismatch);
}

// ============================================================================
// Spike 2: 验证 per-thread（非 scene）trial 寄生能否被局部 STA 一致消费
// 对每个 cell：Path A = updateSteinerTrees(写全局 scene) + analyzeOnly；
// Path B = 把该 cell 的相连 net 寄生构进 per-engine 本地容器（不写 scene/odb）
// + analyzeOnly（findLocalParasiticNetwork 优先读本地）。比较两者到达时间。
// 一致 -> 寄生源可重定向到 per-thread 容器、局部 STA 正确消费、精度吻合。
// 通过 TDP_SPIKE2=1 启用；TDP_SPIKE2_CELLS / _DEPTH 可调。
// ============================================================================
void LSP3D::runTrialParasiticSpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int max_cells = env_int("TDP_SPIKE2_CELLS", 64);
    const int depth     = env_int("TDP_SPIKE2_DEPTH", 2);

    logger_->info(TDP, 1970,
        "[SPIKE2] trial-parasitic: cells<={} depth={}", max_cells, depth);

    sta_->updateTiming(true);
    sta_->findRequireds();

    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) { logger_->warn(TDP, 1971, "[SPIKE2] no scene"); return; }

    SteinerTreeGenerator steiner(timingdp_);

    std::vector<Node3D*> cells;
    std::vector<sta::Instance*> insts;
    for (auto* cell : network_->getCells()) {
        if (static_cast<int>(cells.size()) >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst)
            continue;
        sta::Instance* si = db_network->dbToSta(cell->db_inst);
        if (si) { cells.push_back(cell); insts.push_back(si); }
    }
    if (cells.empty()) { logger_->warn(TDP, 1972, "[SPIKE2] no movable cells"); return; }

    auto engine = std::make_unique<lgista::LocalGraphIsta>(sta_);
    { auto pg = engine->analyzeOnly(insts[0], depth); (void)pg; }  // warm up -> local_sta_
    sta::Parasitics* local_par = engine->localParasiticsContainer();
    if (!local_par) {
        logger_->warn(TDP, 1973, "[SPIKE2] no local parasitics container");
        return;
    }

    auto fingerprint = [&](lgista::PtGraph* pg) -> std::pair<double,int> {
        if (!pg) return {0.0, -1};
        double acc = 0.0; int n = 0;
        for (auto& ptv : pg->ptVertices()) {
            if (!ptv.vertex() || !ptv.paths()) continue;
            double e = 0.0, l = 0.0;
            if (ptVertexLocalArrival(pg, ptv, sta_, e, l)) { acc += e + l; ++n; }
        }
        return {acc, n};
    };

    int mismatch = 0;
    double worst_rel = 0.0;
    const double tol = 1e-3;
    for (size_t k = 0; k < cells.size(); ++k) {
        Node3D* cell = cells[k];

        // 确保 scene 有该 cell 的 net（tdp RC 模型，当前几何）
        steiner.updateSteinerTrees(cell);

        // Path A: 从全局 scene 读寄生
        auto fpA = fingerprint(engine->analyzeOnly(insts[k], depth).get());

        // Path B: 把该 cell 的 net 构进 per-engine 本地容器
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            if (n3) steiner.buildLocalNetParasiticInto(local_par, n3, corner);
        }
        auto fpB = fingerprint(engine->analyzeOnly(insts[k], depth).get());

        // 清理本地容器，供下一个 cell
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (sn) local_par->deleteParasiticNetwork(sn);
        }

        const double a = fpA.first, b = fpB.first;
        const double denom = std::max(std::abs(a), std::abs(b));
        const double rel = denom > 1e-18 ? std::abs(a - b) / denom : 0.0;
        worst_rel = std::max(worst_rel, rel);
        if (fpA.second != fpB.second || rel > tol) {
            ++mismatch;
            logger_->info(TDP, 1974,
                "[SPIKE2] cell {} MISMATCH fpA={:.6e}(n={}) fpB={:.6e}(n={}) rel={:.3e}",
                cell->name, a, fpA.second, b, fpB.second, rel);
        }
    }

    logger_->info(TDP, 1975,
        "[SPIKE2] DONE cells={} mismatch={} worst_rel={:.3e} (tol={:.1e})",
        static_cast<int>(cells.size()), mismatch, worst_rel, tol);
    if (mismatch == 0)
        logger_->info(TDP, 1976,
            "[SPIKE2] PASS: local-container parasitics == scene path within tol");
    else
        logger_->warn(TDP, 1977,
            "[SPIKE2] FAIL: {} cells diverge -> trial-parasitic path NOT faithful",
            mismatch);
}

// ============================================================================
// Spike 3: 验证 Phase-2 真实内循环并发安全
// 每线程一套独立 {LocalGraphIsta 引擎, SteinerTreeGenerator, 寄生容器}，并发对
// 每个 cell 跑 [构 trial 寄生 -> analyzeOnly -> 清理]，与串行做位级对比，并测
// 含 trial-RC 构建的真实并行加速。通过 TDP_SPIKE3=1 启用。
// ============================================================================
void LSP3D::runConcurrentTrialSpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int num_threads = env_int("TDP_SPIKE3_THREADS", 8);
    const int max_cells   = env_int("TDP_SPIKE3_CELLS", 128);
    const int repeat      = env_int("TDP_SPIKE3_REPEAT", 5);
    const int depth       = env_int("TDP_SPIKE3_DEPTH", 2);

    logger_->info(TDP, 1980,
        "[SPIKE3] concurrent trial-RC: threads={} cells<={} repeat={} depth={}",
        num_threads, max_cells, repeat, depth);

    sta_->updateTiming(true);
    sta_->findRequireds();

    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) { logger_->warn(TDP, 1981, "[SPIKE3] no scene"); return; }

    // collect cells + insts + per-cell net lists (serial; no shared lookups in
    // the parallel region beyond what the build itself does)
    std::vector<sta::Instance*> insts;
    std::vector<std::vector<std::pair<Net3D*, sta::Net*>>> cell_nets;
    for (auto* cell : network_->getCells()) {
        if (static_cast<int>(insts.size()) >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst)
            continue;
        sta::Instance* si = db_network->dbToSta(cell->db_inst);
        if (!si) continue;
        std::vector<std::pair<Net3D*, sta::Net*>> nets;
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (n3 && sn) nets.emplace_back(n3, sn);
        }
        insts.push_back(si);
        cell_nets.push_back(std::move(nets));
    }
    const int K = static_cast<int>(insts.size());
    if (K == 0) { logger_->warn(TDP, 1982, "[SPIKE3] no movable cells"); return; }

    // per-thread engines + Steiner generators + parasitic containers
    std::vector<std::unique_ptr<lgista::LocalGraphIsta>> engines;
    std::vector<std::unique_ptr<SteinerTreeGenerator>>   steiners;
    std::vector<sta::Parasitics*>                        local_par(num_threads, nullptr);
    for (int t = 0; t < num_threads; ++t) {
        engines.push_back(std::make_unique<lgista::LocalGraphIsta>(sta_));
        steiners.push_back(std::make_unique<SteinerTreeGenerator>(timingdp_));
    }
    // serial warm-up: create local_sta_ + trigger FLUTE LUT init per thread
    for (int t = 0; t < num_threads; ++t) {
        { auto pg = engines[t]->analyzeOnly(insts[0], depth); (void)pg; }
        local_par[t] = engines[t]->localParasiticsContainer();
        for (auto& pr : cell_nets[0])
            steiners[t]->buildLocalNetParasiticInto(local_par[t], pr.first, corner);
        if (local_par[t])
            for (auto& pr : cell_nets[0])
                local_par[t]->deleteParasiticNetwork(pr.second);
    }

    auto fingerprint = [&](lgista::PtGraph* pg) -> std::pair<double,int> {
        if (!pg) return {0.0, -1};
        double acc = 0.0; int n = 0;
        for (auto& ptv : pg->ptVertices()) {
            if (!ptv.vertex() || !ptv.paths()) continue;
            double e = 0.0, l = 0.0;
            if (ptVertexLocalArrival(pg, ptv, sta_, e, l)) { acc += e + l; ++n; }
        }
        return {acc, n};
    };

    // evaluate cell k using thread tid's private engine/steiner/container
    auto eval = [&](int tid, int k) -> std::pair<double,int> {
        for (auto& pr : cell_nets[k])
            steiners[tid]->buildLocalNetParasiticInto(local_par[tid], pr.first, corner);
        auto fp = fingerprint(engines[tid]->analyzeOnly(insts[k], depth).get());
        if (local_par[tid])
            for (auto& pr : cell_nets[k])
                local_par[tid]->deleteParasiticNetwork(pr.second);
        return fp;
    };

    // serial reference
    std::vector<std::pair<double,int>> ref(K);
    auto ts0 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < K; ++k) ref[k] = eval(0, k);
    auto ts1 = std::chrono::high_resolution_clock::now();
    const double serial_ms =
        std::chrono::duration<double, std::milli>(ts1 - ts0).count();

    int total_mismatch = 0;
    double worst_rel = 0.0;
    double par_ms_best = 1e300;
    for (int r = 0; r < repeat; ++r) {
        std::vector<std::pair<double,int>> par(K, {0.0, -2});
        auto tp0 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (int k = 0; k < K; ++k) {
            int tid = omp_get_thread_num();
            par[k] = eval(tid, k);
        }
        auto tp1 = std::chrono::high_resolution_clock::now();
        const double par_ms =
            std::chrono::duration<double, std::milli>(tp1 - tp0).count();
        par_ms_best = std::min(par_ms_best, par_ms);

        int mismatch = 0;
        for (int k = 0; k < K; ++k) {
            const double a = ref[k].first, b = par[k].first;
            const double denom = std::max(std::abs(a), std::abs(b));
            const double rel = denom > 1e-18 ? std::abs(a - b) / denom : 0.0;
            if (ref[k].second != par[k].second || rel > 1e-9) {
                ++mismatch;
                worst_rel = std::max(worst_rel, rel);
            }
        }
        total_mismatch += mismatch;
        logger_->info(TDP, 1983, "[SPIKE3] run {}: mismatch {}/{}  par_ms={:.1f}",
                      r, mismatch, K, par_ms);
    }

    logger_->info(TDP, 1984,
        "[SPIKE3] DONE cells={} threads={} serial_ms={:.1f} par_ms(best)={:.1f} "
        "speedup={:.2f}x total_mismatch={} worst_rel={:.3e}",
        K, num_threads, serial_ms, par_ms_best,
        (par_ms_best > 0 ? serial_ms / par_ms_best : 0.0),
        total_mismatch, worst_rel);
    if (total_mismatch == 0)
        logger_->info(TDP, 1985,
            "[SPIKE3] PASS: concurrent trial-RC == serial over {} runs", repeat);
    else
        logger_->warn(TDP, 1986,
            "[SPIKE3] FAIL: {} mismatches -> concurrent inner loop NOT safe",
            total_mismatch);
}

// ============================================================================
// Spike 2B (R2): moved-pin override trial geometry vs physically moving odb.
// Reference = setInstanceLocation3D(candidate)+updateSteinerTrees(scene)+
// analyzeOnly; Trial = override the cell's pin on each net to (cur+delta) and
// build into the per-engine local container (odb untouched) + analyzeOnly.
// 通过 TDP_SPIKE2B=1 启用。
// ============================================================================
void LSP3D::runTrialMovedSpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int max_cells = env_int("TDP_SPIKE2B_CELLS", 64);
    const int depth     = env_int("TDP_SPIKE2B_DEPTH", 2);

    logger_->info(TDP, 1998,
        "[SPIKE2B] moved override vs odb-move: cells<={} depth={}", max_cells, depth);

    sta_->updateTiming(true);
    sta_->findRequireds();
    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) { logger_->warn(TDP, 1999, "[SPIKE2B] no scene"); return; }

    SteinerTreeGenerator steiner(timingdp_);

    std::vector<Node3D*> cells;
    std::vector<sta::Instance*> insts;
    for (auto* cell : network_->getCells()) {
        if (static_cast<int>(cells.size()) >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst)
            continue;
        sta::Instance* si = db_network->dbToSta(cell->db_inst);
        if (si) { cells.push_back(cell); insts.push_back(si); }
    }
    if (cells.empty()) { logger_->warn(TDP, 2000, "[SPIKE2B] no cells"); return; }

    auto engine = std::make_unique<lgista::LocalGraphIsta>(sta_);
    { auto pg = engine->analyzeOnly(insts[0], depth); (void)pg; }
    sta::Parasitics* local_par = engine->localParasiticsContainer();
    if (!local_par) { logger_->warn(TDP, 2001, "[SPIKE2B] no local container"); return; }

    auto fingerprint = [&](lgista::PtGraph* pg) -> std::pair<double,int> {
        if (!pg) return {0.0, -1};
        double acc = 0.0; int n = 0;
        for (auto& ptv : pg->ptVertices()) {
            if (!ptv.vertex() || !ptv.paths()) continue;
            double e = 0.0, l = 0.0;
            if (ptVertexLocalArrival(pg, ptv, sta_, e, l)) { acc += e + l; ++n; }
        }
        return {acc, n};
    };

    int tested = 0, mismatch = 0;
    double worst_rel = 0.0;
    const double tol = 1e-3;
    for (size_t k = 0; k < cells.size(); ++k) {
        Node3D* cell = cells[k];
        odb::dbInst* inst = cell->db_inst;
        Pointdp3D cur = cell->location;

        // pick a same-layer candidate that actually moves
        auto cands = timingDrivenWindowSearch(cell, 8, params_.window_x, params_.window_y);
        Pointdp3D cand; bool found = false;
        for (auto& c : cands) {
            if (c.z.v == cur.z.v && (c.x.v != cur.x.v || c.y.v != cur.y.v)) {
                cand = c; found = true; break;
            }
        }
        if (!found) continue;
        ++tested;
        const long dx = static_cast<long>(cand.x.v) - static_cast<long>(cur.x.v);
        const long dy = static_cast<long>(cand.y.v) - static_cast<long>(cur.y.v);

        // Reference: move odb to candidate, build the cell's nets into the
        // LOCAL container (generateLocalRoute reads odb=candidate), with NO
        // updateSteinerTrees / delaysInvalidFromFanin — so the global timing
        // seeds stay frozen, identical to the trial. This isolates the override
        // GEOMETRY from the (separate, intended) frozen-seed approximation.
        timingdp_->setInstanceLocation3D(inst, cand.x.v, cand.y.v, cand.z.v);
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            if (n3) steiner.buildLocalNetParasiticInto(local_par, n3, corner);
        }
        auto fpRef = fingerprint(engine->analyzeOnly(insts[k], depth).get());
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (sn) local_par->deleteParasiticNetwork(sn);
        }
        // restore odb (no scene write needed — we never touched the scene)
        timingdp_->setInstanceLocation3D(inst, cur.x.v, cur.y.v, cur.z.v);

        // Trial: override the cell's pin on each net (odb stays at cur)
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            if (!n3) continue;
            Pin3D* cellpin = nullptr;
            for (auto* p : cell->getPins())
                if (p && p->getNet() == n3) { cellpin = p; break; }
            if (!cellpin) continue;
            Pin3Dlocation pl = steiner.getPin3DLocation(cellpin);
            steiner.buildLocalNetParasiticInto(
                local_par, n3, corner, cellpin,
                static_cast<int>(pl.x + dx), static_cast<int>(pl.y + dy));
        }
        auto fpTrial = fingerprint(engine->analyzeOnly(insts[k], depth).get());
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (sn) local_par->deleteParasiticNetwork(sn);
        }

        const double a = fpRef.first, b = fpTrial.first;
        const double denom = std::max(std::abs(a), std::abs(b));
        const double rel = denom > 1e-18 ? std::abs(a - b) / denom : 0.0;
        worst_rel = std::max(worst_rel, rel);
        if (fpRef.second != fpTrial.second || rel > tol) {
            ++mismatch;
            logger_->info(TDP, 2002,
                "[SPIKE2B] cell {} MISMATCH fpRef={:.6e}(n={}) fpTrial={:.6e}(n={}) "
                "rel={:.3e} d=({},{})",
                cell->name, a, fpRef.second, b, fpTrial.second, rel,
                static_cast<long long>(dx), static_cast<long long>(dy));
        }
    }

    logger_->info(TDP, 2003,
        "[SPIKE2B] DONE tested={} mismatch={} worst_rel={:.3e} (tol={:.1e})",
        tested, mismatch, worst_rel, tol);
    if (mismatch == 0)
        logger_->info(TDP, 2004,
            "[SPIKE2B] PASS: override trial-RC == odb-move within tol");
    else
        logger_->warn(TDP, 2005,
            "[SPIKE2B] FAIL: {} cells diverge -> override geometry NOT equivalent",
            mismatch);
}

// ============================================================================
// Extractor compare: at the CURRENT geometry, compare the local STA result fed
// by the two parasitic extractors:
//   SCENE: updateSteinerTrees(cell) -> estimateTdpWireParasitic (global scene)
//   LOCAL: buildLocalNetParasiticInto(...) -> per-thread container
// Both go through analyzeOnly (-> recomputePtParasitics). Cell nets differ only
// in the extractor used; neighbours read the scene in both. If the fingerprints
// differ, the two extractors are NOT equivalent. TDP_SPIKE_EXT=1.
// ============================================================================
void LSP3D::runExtractorCompareSpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int max_cells = env_int("TDP_SPIKE_EXT_CELLS", 64);
    const int depth     = env_int("TDP_SPIKE_EXT_DEPTH", 2);

    logger_->info(TDP, 2020,
        "[SPIKE_EXT] scene-extractor vs local-extractor (same geometry): cells<={} depth={}",
        max_cells, depth);

    sta_->updateTiming(true);
    sta_->findRequireds();
    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) { logger_->warn(TDP, 2021, "[SPIKE_EXT] no scene"); return; }

    SteinerTreeGenerator steiner(timingdp_);

    std::vector<Node3D*> cells;
    std::vector<sta::Instance*> insts;
    for (auto* cell : network_->getCells()) {
        if (static_cast<int>(cells.size()) >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst)
            continue;
        sta::Instance* si = db_network->dbToSta(cell->db_inst);
        if (si) { cells.push_back(cell); insts.push_back(si); }
    }
    if (cells.empty()) { logger_->warn(TDP, 2022, "[SPIKE_EXT] no cells"); return; }

    auto engine = std::make_unique<lgista::LocalGraphIsta>(sta_);
    { auto pg = engine->analyzeOnly(insts[0], depth); (void)pg; }
    sta::Parasitics* local_par = engine->localParasiticsContainer();
    if (!local_par) { logger_->warn(TDP, 2023, "[SPIKE_EXT] no local container"); return; }

    auto fingerprint = [&](lgista::PtGraph* pg) -> std::pair<double,int> {
        if (!pg) return {0.0, -1};
        double acc = 0.0; int n = 0;
        for (auto& ptv : pg->ptVertices()) {
            if (!ptv.vertex() || !ptv.paths()) continue;
            double e = 0.0, l = 0.0;
            if (ptVertexLocalArrival(pg, ptv, sta_, e, l)) { acc += e + l; ++n; }
        }
        return {acc, n};
    };

    int mismatch = 0, printed = 0;
    double worst_rel = 0.0;
    const double tol = 1e-6;
    for (size_t k = 0; k < cells.size(); ++k) {
        Node3D* cell = cells[k];

        // SCENE extractor
        steiner.updateSteinerTrees(cell);
        auto fpScene = fingerprint(engine->analyzeOnly(insts[k], depth).get());

        // LOCAL extractor (same current geometry)
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            if (n3) steiner.buildLocalNetParasiticInto(local_par, n3, corner);
        }
        auto fpLocal = fingerprint(engine->analyzeOnly(insts[k], depth).get());
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (sn) local_par->deleteParasiticNetwork(sn);
        }

        const double a = fpScene.first, b = fpLocal.first;
        const double denom = std::max(std::abs(a), std::abs(b));
        const double rel = denom > 1e-18 ? std::abs(a - b) / denom : 0.0;
        worst_rel = std::max(worst_rel, rel);
        if (fpScene.second != fpLocal.second || rel > tol) {
            ++mismatch;
            if (printed < 15) {
                logger_->info(TDP, 2024,
                    "[SPIKE_EXT] cell {} DIFF fpScene={:.6e}(n={}) fpLocal={:.6e}(n={}) rel={:.3e}",
                    cell->name, a, fpScene.second, b, fpLocal.second, rel);
                // Per-net total-cap / node-count / resistor-count: scene vs local.
                struct NetPar { double cap; int nodes; int res; double totR; };
                auto capNodes = [](sta::Parasitics* par, sta::Net* net) -> NetPar {
                    if (!par) return {-1.0, -1, -1, -1.0};
                    sta::Parasitic* nw = par->findParasiticNetwork(net);
                    if (!nw) return {-2.0, -2, -2, -2.0};
                    double totR = 0.0;
                    auto rs = par->resistors(nw);
                    for (auto* r : rs) totR += par->value(r);
                    return {(double)par->capacitance(nw), (int)par->nodes(nw).size(),
                            (int)rs.size(), totR};
                };
                sta::Parasitics* pmax = corner->parasitics(sta::MinMax::max());
                for (auto* dbnet : cell->connected_nets) {
                    if (!dbnet) continue;
                    sta::Net* sn = db_network->dbToSta(dbnet);
                    if (!sn) continue;
                    Net3D* n3 = network_->getNet(dbnet);
                    if (n3) steiner.buildLocalNetParasiticInto(local_par, n3, corner);
                    NetPar s = capNodes(pmax, sn);
                    NetPar l = capNodes(local_par, sn);
                    if (n3) local_par->deleteParasiticNetwork(sn);
                    if (s.nodes != l.nodes || s.res != l.res
                        || std::abs(s.totR - l.totR) > 1e-6 * std::max(std::abs(s.totR), std::abs(l.totR)))
                        logger_->info(TDP, 2028,
                            "[SPIKE_EXT]   net {} scene[cap={:.3e} nodes={} res={} R={:.3e}] local[cap={:.3e} nodes={} res={} R={:.3e}]",
                            dbnet->getName(), s.cap, s.nodes, s.res, s.totR,
                            l.cap, l.nodes, l.res, l.totR);
                }
                ++printed;
            }
        }
    }

    logger_->info(TDP, 2025,
        "[SPIKE_EXT] DONE cells={} mismatch={} worst_rel={:.3e} (tol={:.1e})",
        static_cast<int>(cells.size()), mismatch, worst_rel, tol);
    if (mismatch == 0)
        logger_->info(TDP, 2026,
            "[SPIKE_EXT] PASS: scene-extractor == local-extractor (same geometry)");
    else
        logger_->warn(TDP, 2027,
            "[SPIKE_EXT] DIVERGENCE: {} cells differ -> the two parasitic extractors are NOT equivalent",
            mismatch);
}

// ============================================================================
// Cost attribution spike: SPIKE_EXT and SPIKE2B already proved the trial STA /
// geometry / extractor are identical between par and no-commit. So any par-vs-
// no-commit divergence must live in the COST-CONSUMPTION code. For K movable
// cells, at a same-layer moved candidate built with the (proven-faithful) par
// override trial, on ONE shared PtGraph, we compute the candidate cost via:
//   cost_par    = buildLocalArrivalMap (driver INCLUDED) + calculateCellTimingCostLocal
//   cost_nodrv  = same map but drivers REMOVED           + calculateCellTimingCostLocal
//   cost_nc     = refreshAffectedPinTimingFromPtGraph (driver NOT refreshed) + calculateCellCost
// and decompose:
//   driver-frame effect (③) = cost_par  - cost_nodrv
//   arc frozen/live effect (②) = cost_nodrv - cost_nc
// Plus baseline effect (⑥) at dx=0: |cost_par(cur) - calculateCellCost(cur)|.
// TDP_SPIKE_COST=1 ; TDP_SPIKE_COST_CELLS to cap.
// ============================================================================
void LSP3D::runCostAttributionSpike()
{
    auto env_int = [](const char* k, int def) {
        const char* v = std::getenv(k);
        return v ? std::max(1, atoi(v)) : def;
    };
    const int max_cells = env_int("TDP_SPIKE_COST_CELLS", 200);
    const int depth     = 2;
    const double alpha = 5.0, omega = 1.0;

    logger_->info(TDP, 2040,
        "[SPIKE_COST] par-vs-no-commit cost attribution: cells<={} depth={}",
        max_cells, depth);

    // Weights (=1.0) + a clean global timing baseline in the Pin3D caches.
    initWeight(timingdp_, timing_engine_);
    sta_->updateTiming(true);
    sta_->findRequireds();
    for (auto* c : network_->getCells()) refreshPinTiming2(c);

    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) { logger_->warn(TDP, 2041, "[SPIKE_COST] no scene"); return; }
    ensureParallelPool(1);
    lgista::LocalGraphIsta* engine = par_engines_[0].get();
    SteinerTreeGenerator*   steiner = par_steiners_[0].get();
    sta::Parasitics*        lp = par_local_par_[0];
    if (!engine || !steiner || !lp) { logger_->warn(TDP, 2042, "[SPIKE_COST] pool init failed"); return; }

    int tested = 0;
    double sum_gap = 0.0, sum_driver = 0.0, sum_arc = 0.0, sum_base = 0.0;
    double worst_gap = 0.0;
    int printed = 0, driver_reconciles = 0;
    double sum_realgap = 0.0; int real_tested = 0;
    const int real_cap = env_int("TDP_SPIKE_COST_REAL", 40);

    for (auto* cell : network_->getCells()) {
        if (tested >= max_cells) break;
        if (!cell || cell->is_fixed || isSequentialCell(cell) || !cell->db_inst) continue;
        sta::Instance* inst_sta = db_network->dbToSta(cell->db_inst);
        if (!inst_sta) continue;

        // Prime arc caches at the ORIGINAL position (matches par line ~2260).
        calculateCellCost(cell, alpha, omega, timing_engine_);

        // Pick a same-layer candidate that actually moves.
        auto cands = timingDrivenWindowSearch(cell, 12, params_.window_x, params_.window_y);
        Pointdp3D cur = cell->location, cand; bool found = false;
        for (auto& c : cands)
            if (c.z.v == cur.z.v && (c.x.v != cur.x.v || c.y.v != cur.y.v)) { cand = c; found = true; break; }
        if (!found) continue;
        const long dx = (long)cand.x.v - (long)cur.x.v;
        const long dy = (long)cand.y.v - (long)cur.y.v;

        // Collect the cell's pin on each connected net (+ its current coord).
        struct TNet { Net3D* net; sta::Net* sn; Pin3D* pin; int px; int py; };
        std::vector<TNet> tnets;
        for (auto* dbnet : cell->connected_nets) {
            if (!dbnet) continue;
            Net3D* n3 = network_->getNet(dbnet);
            sta::Net* sn = db_network->dbToSta(dbnet);
            if (!n3 || !sn) continue;
            Pin3D* cp = nullptr;
            for (auto* p : cell->getPins()) if (p && p->getNet() == n3) { cp = p; break; }
            if (!cp) continue;
            Pin3Dlocation pl = steiner->getPin3DLocation(cp);
            tnets.push_back({n3, sn, cp, (int)pl.x, (int)pl.y});
        }
        if (tnets.empty()) continue;

        auto build_trial = [&](long ddx, long ddy) {
            for (auto& tn : tnets)
                steiner->buildLocalNetParasiticInto(lp, tn.net, corner, tn.pin,
                                                    tn.px + (int)ddx, tn.py + (int)ddy);
        };
        auto clear_trial = [&]() { for (auto& tn : tnets) lp->deleteParasiticNetwork(tn.sn); };

        // Identify driver pins of the cell's input nets (par includes these).
        std::unordered_set<Pin3D*> drivers;
        for (auto* p : cell->getPins()) {
            if (!p || p->getDirection() != Pin3D::Dir_IN) continue;
            Net3D* net = p->getNet();
            if (!net) continue;
            if (Pin3D* d = findDriverPin(net, p)) drivers.insert(d);
        }

        // ---- baseline (⑥): par local path at cur vs no-commit global cache ----
        build_trial(0, 0);
        {
            auto pg0 = engine->analyzeOnly(inst_sta, depth);
            std::unordered_map<Pin3D*, std::pair<double,double>> arr0, arc0;
            buildLocalArrivalMap(cell, pg0.get(), arr0);
            buildLocalArcMap(cell, pg0.get(), arc0);
            double base_par = calculateCellTimingCostLocal(cell, alpha, omega, arr0, arc0);
            double base_nc  = calculateCellCost(cell, alpha, omega, timing_engine_);
            sum_base += std::abs(base_par - base_nc);
        }
        clear_trial();

        // ---- candidate (moved): three costs on the SAME PtGraph ----
        build_trial(dx, dy);
        auto pg = engine->analyzeOnly(inst_sta, depth);

        std::unordered_map<Pin3D*, std::pair<double,double>> arr;
        buildLocalArrivalMap(cell, pg.get(), arr);
        std::unordered_map<Pin3D*, std::pair<double,double>> arc_map;
        buildLocalArcMap(cell, pg.get(), arc_map);
        double cost_par = calculateCellTimingCostLocal(cell, alpha, omega, arr, arc_map);

        std::unordered_map<Pin3D*, std::pair<double,double>> arr_nodrv = arr;
        for (auto* d : drivers) arr_nodrv.erase(d);
        double cost_nodrv = calculateCellTimingCostLocal(cell, alpha, omega, arr_nodrv, arc_map);

        refreshAffectedPinTimingFromPtGraph(cell, pg.get());
        double cost_nc = calculateCellCost(cell, alpha, omega, timing_engine_);

        clear_trial();

        // Restore the Pin3D caches that refreshAffectedPinTimingFromPtGraph wrote
        // (odb was never moved, so global STA is unchanged -> clean restore).
        refreshPinTiming2(cell);
        for (auto* p : cell->getPins()) {
            if (p->getDirection() != Pin3D::Dir_OUT) continue;
            Net3D* net = p->getNet(); if (!net) continue;
            for (auto* np : net->getPins()) { Node3D* nc = np->getNode(); if (nc && nc != cell) refreshPinTiming2(nc); }
        }

        const double gap    = cost_par - cost_nc;
        const double driver = cost_par - cost_nodrv;   // ③
        const double arc    = cost_nodrv - cost_nc;    // ②
        sum_gap += std::abs(gap); sum_driver += std::abs(driver); sum_arc += std::abs(arc);
        worst_gap = std::max(worst_gap, std::abs(gap));
        if (std::abs(gap) > 1e-18 && std::abs(driver) >= 0.9 * std::abs(gap)) driver_reconciles++;
        ++tested;

        // ---- REAL no-commit path: physically move odb, then compute cost the
        // way no-commit actually does. Compares against par's override cost to
        // expose the arc(②)/moved-geometry effect that the fixed-geometry
        // comparison above cannot see. Heavy -> capped by TDP_SPIKE_COST_REAL.
        if (real_tested < real_cap) {
            odb::dbInst* inst = cell->db_inst;
            timingdp_->setInstanceLocation3D(inst, cand.x.v, cand.y.v, cand.z.v);
            steiner->updateSteinerTrees(cell);
            {
                std::unique_ptr<lgista::PtGraph> pgm(timingdp_->analyzeOnly(cell, depth));
                refreshAffectedPinTimingFromPtGraph(cell, pgm.get());
            }
            double cost_nc_real = calculateCellCost(cell, alpha, omega, timing_engine_);
            // restore odb + scene + timing
            timingdp_->setInstanceLocation3D(inst, cur.x.v, cur.y.v, cur.z.v);
            steiner->updateSteinerTrees(cell);
            sta_->updateTiming(false);
            refreshPinTiming2(cell);
            for (auto* p : cell->getPins()) {
                if (p->getDirection() != Pin3D::Dir_OUT) continue;
                Net3D* net = p->getNet(); if (!net) continue;
                for (auto* np : net->getPins()) { Node3D* nc = np->getNode(); if (nc && nc != cell) refreshPinTiming2(nc); }
            }
            sum_realgap += std::abs(cost_par - cost_nc_real);
            if (printed <= 12)
                logger_->info(TDP, 2048,
                    "[SPIKE_COST]   REAL-MOVE {} par(override)={:.6e} nc(odb-move)={:.6e} gap={:.6e}",
                    cell->name, cost_par, cost_nc_real, cost_par - cost_nc_real);
            ++real_tested;
        }

        if (printed < 12) {
            logger_->info(TDP, 2043,
                "[SPIKE_COST] {} par={:.6e} nc={:.6e} gap={:.6e} | driver(③)={:.6e} arc(②)={:.6e} #drv={}",
                cell->name, cost_par, cost_nc, gap, driver, arc, (int)drivers.size());
            ++printed;
        }
    }

    if (tested == 0) { logger_->warn(TDP, 2044, "[SPIKE_COST] no movable cells tested"); return; }
    logger_->info(TDP, 2045,
        "[SPIKE_COST] DONE tested={} | mean|candidate gap|={:.6e} worst={:.6e}",
        tested, sum_gap / tested, worst_gap);
    logger_->info(TDP, 2046,
        "[SPIKE_COST]   decomposition: mean|driver-frame ③|={:.6e}  mean|arc ②|={:.6e}  mean|baseline ⑥|={:.6e}",
        sum_driver / tested, sum_arc / tested, sum_base / tested);
    logger_->info(TDP, 2047,
        "[SPIKE_COST]   cells where driver-frame ③ explains >=90%% of the gap: {}/{}",
        driver_reconciles, tested);
    if (real_tested > 0)
        logger_->info(TDP, 2049,
            "[SPIKE_COST]   REAL-MOVE (odb-move) mean|par - no-commit gap|={:.6e} over {} cells",
            sum_realgap / real_tested, real_tested);
}

void LSP3D::updatePinWeights()
{
    if (!network_) {
        logger_->error(TDP, 1306, "Error: Network not initialized");
        return;
    }
    
    // 获取所有单元
    // const auto& all_cells = network_->getCells();
    
    int updated_pins = 0;
    // 遍历所有单元
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        refreshPinTiming2(cell);
        
        // 获取单元的所有引脚
        auto cell_pins = cell->getPins();

        if(isSequentialCell(cell)){


            for (auto* pin : cell_pins) {
                if (!pin || pin->getDirection() != Pin3D::Dir_IN) continue;

                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && (timingdp_->isIdealClockNet(net->db_net) || !timingdp_->isValidSigType(net->db_net) || net->db_net->isSpecial())) continue;

                double current_late_weight = pin->getLateWeight();
                double current_early_weight = pin->getEarlyWeight();
                
                // 获取时序信息
                double late_arrival = pin->getLateArrival();
                double late_required = pin->getLateRequired();
                double early_arrival = pin->getEarlyArrival();
                double early_required = pin->getEarlyRequired();
                
                // 避免除零错误
                if (late_required != 0.0 && early_arrival != 0.0) {

                    if(early_required < 0.0){
                            early_required = early_arrival / 10;
                    }

                    double new_late_weight = current_late_weight * (late_arrival / late_required);
                    double new_early_weight = current_early_weight * (early_required / early_arrival); 
                    if(new_late_weight < 0.0  || new_early_weight < 0.0){
                        logger_ -> info(TDP, 1321, "DDX weight flase");
                    }                  
                    pin->setLateWeight(new_late_weight);
                    pin->setEarlyWeight(new_early_weight);
                }
            }
            continue;

        }
        else{

            if (cell->getBTerm() != nullptr){
            for (auto* pin : cell_pins) {
                if (!pin) continue;
                if (pin->getNet() && pin->getNet()->db_net && !timingdp_->isValidSigType(pin->getNet()->db_net)) continue;
                updated_pins++;
                // 只更新输出引脚
                // if (pin->getDbBTerm() != nullptr && pin->getDirection() == Pin3D::Dir_OUT) {
                if (pin->getDbBTerm() != nullptr && pin->getDirection() == Pin3D::Dir_OUT) {
                //     std::cout << "this is outport: " << pin->getName() << std::endl;}
                // if (pin->getDirection() == Pin3D::Dir_OUT) {
                    // 获取当前权重
                    double current_late_weight = pin->getLateWeight();
                    double current_early_weight = pin->getEarlyWeight();
                    
                    // 获取时序信息
                    double late_arrival = pin->getLateArrival();
                    double late_required = pin->getLateRequired();
                    double early_arrival = pin->getEarlyArrival();
                    double early_required = pin->getEarlyRequired();
                    
                    // 避免除零错误
                    if (late_required != 0.0 && early_arrival != 0.0) {

                        if(early_required < 0.0){
                            early_required = early_arrival / 10;
                        }

                        double new_late_weight = current_late_weight * (late_arrival / late_required);
                        double new_early_weight = current_early_weight * (early_required / early_arrival);                   
                        if(new_late_weight < 0.0 || new_early_weight <0.0){
                            logger_ -> info(TDP, 1322, "PO weight flase");
                        }                            
                        pin->setLateWeight(new_late_weight);
                        pin->setEarlyWeight(new_early_weight);
                        

                    }
                }
            }
            }
            

        }
        
 
        balanceCellWeights(cell);      // commented
        distributeCellWeights(cell);   // commented
    }
    
    if (std::getenv("TDP_DIAG_W1")) {
        std::vector<double> all_w, crit_w;
        for (auto* c : network_->getCells()) {
            if (!c) continue;
            for (auto* p : c->getPins())
                if (p && p->getDirection() == Pin3D::Dir_IN)
                    all_w.push_back(p->getLateWeight());
        }
        auto crit = buildCriticalPathPriorityOrder(timing_engine_, 1000);
        for (auto* c : crit) {
            if (!c) continue;
            for (auto* p : c->getPins())
                if (p && p->getDirection() == Pin3D::Dir_IN)
                    crit_w.push_back(p->getLateWeight());
        }
        diagLogStats(logger_, "W1 all-in-late", all_w);
        diagLogStats(logger_, "W1 crit-in-late", crit_w);
    }

    std::cout << "Updated weights for " << updated_pins << " output pins" << std::endl;
}

// ---- LR3 weight scheme: wl (slack/WNS) + wp (propagated) — disabled via #if 0
#if 0
void LSP3D::updatePinWeights__LR3()
{
    if (!network_ || !timing_engine_) return;
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        refreshPinTiming2(cell);
    }
    computeWl();
    propagateWp();
}
void LSP3D::computeWl()
{
    auto stats = timing_engine_->getTimingStats();
    double late_wns = stats.late_wns;
    double early_wns = stats.early_wns;

    for (auto* cell : traverse_list) {
        if (!cell) continue;
        for (auto* pin : cell->getPins()) {
            if (!pin) continue;

            // Late wl = max(0, min(1, slack/WNS))
            pin->setLateWeight(0.0);
            double lr = pin->getLateRequired(), la = pin->getLateArrival();
            if (lr != 0.0 && late_wns < 0.0) {
                double slack = lr - la;
                pin->setLateWeight(std::max(0.0, std::min(1.0, slack / late_wns)));
            }

            // Early wl = max(0, min(1, slack/WNS))
            pin->setEarlyWeight(0.0);
            double er = pin->getEarlyRequired(), ea = pin->getEarlyArrival();
            if (er != 0.0 && ea != 0.0 && early_wns < 0.0) {
                double slack = ea - er;
                pin->setEarlyWeight(std::max(0.0, std::min(1.0, slack / early_wns)));
            }
        }
    }
}

void LSP3D::propagateWp()
{
    // Phase 0: wp = 0 init
    for (auto* cell : traverse_list)
        for (auto* pin : cell->getPins())
            pin->setLateWeightini(0.0), pin->setEarlyWeightini(0.0);

    // Phase 1: endpoint init wp = wl
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        // sequential D pin = timing endpoint
        if (isSequentialCell(cell)) {
            for (auto* pin : cell->getPins())
                if (pin->getDirection() == Pin3D::Dir_IN)
                    pin->setLateWeightini(pin->getLateWeight()),
                    pin->setEarlyWeightini(pin->getEarlyWeight());
        }
        // BTerm output pin = timing endpoint
        if (cell->getBTerm()) {
            for (auto* pin : cell->getPins())
                if (pin->getDbBTerm() && pin->getDirection() == Pin3D::Dir_OUT)
                    pin->setLateWeightini(pin->getLateWeight()),
                    pin->setEarlyWeightini(pin->getEarlyWeight());
        }
    }

    // Phase 2: reverse topological propagation
    for (auto* cell : traverse_list) {
        if (!cell) continue;

        // A: aggregate output pin wp = Σ sink wp
        double sum_L_out = 0, sum_E_out = 0;
        for (auto* out_pin : getCellOutputPins(cell)) {
            if (!out_pin) continue;
            Net3D* net = out_pin->getNet();
            if (!net || !net->db_net || !timingdp_->isValidSigType(net->db_net)) continue;

            double wpL = 0, wpE = 0;
            for (auto* sink : net->getPins()) {
                if (sink == out_pin) continue;
                wpL += sink->getLateWeightini();
                wpE += sink->getEarlyWeightini();
            }
            if (wpL == 0) wpL = out_pin->getLateWeight();
            if (wpE == 0) wpE = out_pin->getEarlyWeight();
            out_pin->setLateWeightini(wpL), out_pin->setEarlyWeightini(wpE);
            sum_L_out += wpL, sum_E_out += wpE;
        }

        // B: distribute to input pins by wl ratio
        if (sum_L_out > 0 || sum_E_out > 0) {
            auto in_pins = getCellInputPins(cell);
            double sum_wl_L = 0, sum_wl_E = 0;
            for (auto* in : in_pins) {
                if (!in) continue;
                Net3D* n = in->getNet();
                if (!n || !n->db_net || !timingdp_->isValidSigType(n->db_net)) continue;
                if (isSequentialCell(cell) && in->getDirection() == Pin3D::Dir_IN) continue;
                sum_wl_L += in->getLateWeight(), sum_wl_E += in->getEarlyWeight();
            }
            for (auto* in : in_pins) {
                if (!in) continue;
                Net3D* n = in->getNet();
                if (!n || !n->db_net || !timingdp_->isValidSigType(n->db_net)) continue;
                if (isSequentialCell(cell) && in->getDirection() == Pin3D::Dir_IN) continue;
                if (sum_wl_L > 0)
                    in->setLateWeightini(sum_L_out * in->getLateWeight() / sum_wl_L);
                if (sum_wl_E > 0)
                    in->setEarlyWeightini(sum_E_out * in->getEarlyWeight() / sum_wl_E);
            }
        }
    }

    // Phase 3: normalize by max
    double max_L = 1e-15, max_E = 1e-15;
    for (auto* cell : traverse_list)
        for (auto* pin : cell->getPins())
            max_L = std::max(max_L, pin->getLateWeightini()),
            max_E = std::max(max_E, pin->getEarlyWeightini());
    for (auto* cell : traverse_list)
        for (auto* pin : cell->getPins()) {
            if (max_L > 1e-15) pin->setLateWeightini(pin->getLateWeightini() / max_L);
            if (max_E > 1e-15) pin->setEarlyWeightini(pin->getEarlyWeightini() / max_E);
        }
}

void LSP3D::combineWeights()
{
    double b = params_.beta;
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        for (auto* pin : cell->getPins()) {
            if (!pin) continue;
            pin->setLateWeight(pin->getLateWeight() + b * pin->getLateWeightini());
            pin->setEarlyWeight(pin->getEarlyWeight() + b * pin->getEarlyWeightini());
        }
    }
}
#endif

/**
 * @brief 平衡单元输入引脚net权重
 * 
 * 对每个输入引脚 i，找到上游驱动引脚 driver，对每个输出引脚 j：
 *   λ_i,j^L = λ_i,j^L × (a_i^L + d_i,j^L) / a_j^L
 *   λ_i,j^E = λ_i,j^E × a_i^E / (a_i^E + d_i,j^E)
 *   λ_i,j^L = λ_i,j^L × ( 1 - ( r_j^L - ( a_i^L + d_i,j^L ) ) / T )^K
 *   λ_i,j^E = λ_i,j^E × ( 1 - ( ( a_i^E + d_i,j^E ) - r_j^E ) / T )^K
 * 对每路输出 j 的结果取算术平均，更新输入引脚权重。
 */
void LSP3D::balanceCellWeights(Node3D* cell) {
    if (!cell || !network_ || !timing_engine_) {
        return;
    }

    auto input_pins  = getCellInputPins(cell);
    auto output_pins = getCellOutputPins(cell);
    if (input_pins.empty() || output_pins.empty()) {
        return;
    }

    for (auto* in_pin : input_pins) {
        if (!in_pin) continue;
        Net3D* net = in_pin->getNet();
        if (!net) continue;
        if (net->db_net && (timingdp_->isIdealClockNet(net->db_net) || !timingdp_->isValidSigType(net->db_net) || net->db_net->isSpecial())) {
            continue;
        }

        // Pin3D* driver = findDriverPin(in_pin->getNet());
        // if (!driver) continue;

        double a_i_L = in_pin->getLateArrival();
        double a_i_E = in_pin->getEarlyArrival();

        double sum_late  = 0.0;
        double sum_early = 0.0;
        int    cnt       = 0;

        for (auto* out_pin : output_pins) {
            if (!out_pin) continue;

            auto arc = timing_engine_->getArcDelay(in_pin->getName(),
                                                    out_pin->getName());
            if (!arc.is_valid) continue;

            // double a_j_L = out_pin->getLateArrival();
            // double a_j_E = out_pin->getEarlyArrival();
            // if (a_j_L != 0.0) {
            //     sum_late += in_pin->getLateWeight()
            //               * (a_i_L + arc.late_delay) / a_j_L;
            // }
            // if (a_i_E + arc.early_delay != 0.0) {
            //     sum_early += in_pin->getEarlyWeight()
            //                * a_j_E / (a_i_E + arc.early_delay);
            // }

            double r_j_L = out_pin->getLateRequired();
            double r_j_E = out_pin->getEarlyRequired();

            double late_slack  = r_j_L - (a_i_L + arc.late_delay);
            double early_slack = (a_i_E + arc.early_delay) - r_j_E;

            constexpr double K = 4.0;
            double late_factor  = pow(std::max(1e-10, 1.0 - late_slack  / clock_period_T_), K);
            double early_factor = pow(std::max(1e-10, 1.0 - early_slack / clock_period_T_), K);

            sum_late  += in_pin->getLateWeight()  * late_factor;
            sum_early += in_pin->getEarlyWeight() * early_factor;
            cnt++;
        }

        if (cnt > 0) {
            if(sum_late < 0.0 || sum_early < 0.0){
                logger_ -> info(TDP, 1323, "weight flase");
            }    
            in_pin->setLateWeight(sum_late / cnt);
            in_pin->setEarlyWeight(sum_early / cnt);
        }
    }
}

/**
 * @brief 在单元内部重新分布权重，确保输入总权重和输出总权重平衡
 * @param cell 要重新分布权重的单元
 */
void LSP3D::distributeCellWeights(Node3D* cell) {
    if (!cell || !network_) {
        return;
    }
    
    // 初始化权重总和
    double in_late_weight = 0.0;
    double in_early_weight = 0.0;
    double out_late_weight = 0.0;
    double out_early_weight = 0.0;
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    // 1. 计算输入引脚和输出引脚的总权重
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_IN) {
            Net3D* net = pin->getNet();
            if (!net) continue;
            if (net->db_net && (timingdp_->isIdealClockNet(net->db_net) || !timingdp_->isValidSigType(net->db_net) || net->db_net->isSpecial())) {
                continue;
            }
            in_late_weight += pin->getLateWeight();
            in_early_weight += pin->getEarlyWeight();
        } else if (pin->getDirection() == Pin3D::Dir_OUT) {
            // 输出引脚
            // out_late_weight += pin->getLateWeight();
            // out_early_weight += pin->getEarlyWeight();
            Net3D* net = pin->getNet();
                if (net) {
                    auto net_pins = net->getPins();
                    for (auto* sink_pin : net_pins) {
                        if (sink_pin == pin) continue;
                        out_late_weight += sink_pin->getLateWeight();
                        out_early_weight += sink_pin->getEarlyWeight();
                    }
                }
        }
    }
    // 2. 重新分布输入引脚权重
    if (in_late_weight > 0.0 && in_early_weight > 0.0) {
        // 计算权重调整比例
        double late_ratio = (out_late_weight > 0.0) ? (out_late_weight / in_late_weight) : 1.0;
        double early_ratio = (out_early_weight > 0.0) ? (out_early_weight / in_early_weight) : 1.0;
        
        //std::cout << "  Adjustment ratios: late=" << late_ratio << ", early=" << early_ratio << std::endl;
        
        // 调整所有输入引脚权重
        for (auto* pin : cell_pins) {
            if (!pin || pin->getDirection() != Pin3D::Dir_IN) {
                continue;
            }
            
            // 保存原始权重用于调试
            double old_late = pin->getLateWeight();
            double old_early = pin->getEarlyWeight();
            
            // 应用调整
            double new_late = old_late * late_ratio;
            double new_early = old_early * early_ratio;
            
            pin->setLateWeight(new_late);
            pin->setEarlyWeight(new_early);

        }
    }
}
// ============================================================================
// 窗口寻找
// ============================================================================
void LSP3D::computeTimingVectors(
    Node3D* cell_j,
    Node3D* connected_cell,
    bool is_fan_in,
    double late_lm,
    double early_lm,
    TimingVector& late_vector,
    TimingVector& early_vector) {
    
    if (!cell_j || !connected_cell) {
        return;
    }
    
    // 获取单元中心位置
    Pointdp3D center_j = getCellCenter(cell_j);
    Pointdp3D center_conn = getCellCenter(connected_cell);
    
    // 计算方向向量（从cell_j指向connected_cell）
    double dx = center_conn.x.v - center_j.x.v;
    double dy = center_conn.y.v - center_j.y.v;
    
    // 计算距离
    double distance = std::sqrt(dx * dx + dy * dy);
    if (distance < 1e-10) {
        return;  // 避免除零
    }
    
    // 单位向量
    double unit_x = dx / distance;
    double unit_y = dy / distance;
    
    if (is_fan_in) {
        // 扇入单元：
        // - 晚时序向量：指向互连单元（吸引力）
        // - 早时序向量：背向互连单元（排斥力）
        late_vector = TimingVector(unit_x * late_lm, unit_y * late_lm, late_lm);
        early_vector = TimingVector(-unit_x * early_lm, -unit_y * early_lm, early_lm);
    } else {
        // 扇出单元：
        // - 晚时序向量：指向互连单元（吸引力）
        // - 早时序向量：背向互连单元（排斥力）
        late_vector = TimingVector(unit_x * late_lm, unit_y * late_lm, late_lm);
        early_vector = TimingVector(-unit_x * early_lm, -unit_y * early_lm, early_lm);
    }
}

/**
 * @brief 获取单元中心位置
 */
Pointdp3D LSP3D::getCellCenter(Node3D* cell) const {
    if (!cell) {
        return Pointdp3D(0, 0, 0);
    }
    
    DBU center_x = cell->location.x.v + cell->getWidth().v / 2;
    DBU center_y = cell->location.y.v + cell->getHeight().v / 2;
    
    return Pointdp3D(DbuX{center_x}, DbuY{center_y}, cell->location.z);
}

/**
 * @brief 计算合成向量和方向分量
 * @param cell 当前单元
 * @param late_components 输出的晚时序方向分量
 * @param early_components 输出的早时序方向分量
 */
bool LSP3D::computeDirectionComponents(
    Node3D* cell,
    DirectionComponents& late_components,
    DirectionComponents& early_components) {
    
    if (!cell || !network_) {
        return false;
    }
    
    // 初始化合成向量
    TimingVector late_synthetic(0.0, 0.0, 0.0);
    TimingVector early_synthetic(0.0, 0.0, 0.0);
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        Net3D* net = pin->getNet();
        if (!net) continue;
        
        // 获取引脚方向
        int dir_int = pin->getDirection();
        Pin3D::Direction dir = static_cast<Pin3D::Direction>(dir_int);
        
        
        if (dir == Pin3D::Dir_IN) {
            // 输入引脚：查找驱动单元
            // 获取LM值（这里需要从引脚或网络获取）
            double late_lm = pin->getLateWeight();  // 使用权重作为LM的代理
            double early_lm = pin->getEarlyWeight();
            Pin3D* driver_pin = findDriverPin(net, pin);
            if (!driver_pin) continue;
            
            Node3D* driver_cell = driver_pin->getNode();
            if (!driver_cell) continue;
            
            // 计算向量
            TimingVector late_vec, early_vec;
            computeTimingVectors(cell, driver_cell, true, late_lm, early_lm, 
                                late_vec, early_vec);
            
            // 累加到合成向量
            late_synthetic.x += late_vec.x;
            late_synthetic.y += late_vec.y;
            early_synthetic.x += early_vec.x;
            early_synthetic.y += early_vec.y;
            
        } else if (dir == Pin3D::Dir_OUT) {
            // 输出引脚：查找所有接收单元
            auto net_pins = net->getPins();
            
            for (auto* sink_pin : net_pins) {
                if (sink_pin == pin) continue;  // 跳过驱动引脚
                // 获取LM值（这里需要从引脚或网络获取）
                double late_lm = sink_pin->getLateWeight();  // 使用权重作为LM的代理
                double early_lm = sink_pin->getEarlyWeight();
                
                Node3D* sink_cell = sink_pin->getNode();
                if (!sink_cell) continue;
                
                // 计算向量
                TimingVector late_vec, early_vec;
                computeTimingVectors(cell, sink_cell, false, late_lm, early_lm,
                                    late_vec, early_vec);
                
                // 累加到合成向量
                late_synthetic.x += late_vec.x;
                late_synthetic.y += late_vec.y;
                early_synthetic.x += early_vec.x;
                early_synthetic.y += early_vec.y;
            }
        }
    }
    
    // 分解为方向分量
    decomposeToComponents(late_synthetic, late_components);
    decomposeToComponents(early_synthetic, early_components);
    
    return true;
}

/**
 * @brief 将向量分解为方向分量
 */
void LSP3D::decomposeToComponents(
    const TimingVector& vector,
    DirectionComponents& components) {
    
    // 水平分量
    if (vector.x > 0) {
        components.right += std::abs(vector.x);
    } else {
        components.left +=  std::abs(vector.x);
    }
    
    // 垂直分量
    if (vector.y > 0) {
        components.up += std::abs(vector.y);
    } else {
        components.down +=  std::abs(vector.y);
    }
}

/**
 * @brief 确定搜索窗口位置
 * @param cell 当前单元
 * @param window_width 搜索窗口宽度
 * @param window_height 搜索窗口高度
 * @return 搜索窗口
 */
LSP3D::SearchWindow LSP3D::determineSearchWindow(
    Node3D* cell,
    DBU window_width,
    DBU window_height) {
    
    if (!cell) {
        return SearchWindow();
    }
    
    // 计算方向分量
    DirectionComponents late_components, early_components;
    if (!computeDirectionComponents(cell, late_components, early_components)) {
        // 如果计算失败，返回以单元为中心的默认窗口
        // std::cerr << "Error: No find cell windows " << std::endl;
        DBU x0 = cell->location.x.v - window_width / 2;
        DBU y0 = cell->location.y.v - window_height / 2;
        return SearchWindow(x0, y0, window_width, window_height);
    }
    
    // 合并晚时序和早时序分量
    DirectionComponents total_components;
    total_components.left = late_components.left + early_components.left;
    total_components.right = late_components.right + early_components.right;
    total_components.up = late_components.up + early_components.up;
    total_components.down = late_components.down + early_components.down;
    
    // 获取单元左下角位置
    DBU xj = cell->location.x.v;
    DBU yj = cell->location.y.v;
    
    // 初始窗口位置（单元位于窗口中心）
    DBU x0 = xj - window_width / 2 ;
    DBU y0 = yj - window_height / 2 ;
    
    // 计算调整后的窗口位置
    DBU x0_adjusted = x0;
    DBU y0_adjusted = y0;
    
    // 计算水平调整
    double total_horizontal = total_components.left - total_components.right;
    if (total_components.left >= 0.0 && total_components.right >= 0.0) {
        double left_ratio = total_horizontal / (std::max(total_components.left, total_components.right) * 2) ;
        x0_adjusted = x0 - left_ratio * window_width;
    }
    else{
        logger_->info(TDP, 1543, "window error");
    }
    
    // 计算垂直调整
    double total_vertical = total_components.down - total_components.up;
    if (total_components.down >= 0.0 &&  total_components.up >= 0.0) {
        double down_ratio = total_vertical / (std::max(total_components.down, total_components.up) * 2);
        y0_adjusted = y0 - down_ratio * window_height;
    }
    else{
        logger_->info(TDP, 1544, "window error");
    }
    return SearchWindow(x0_adjusted, y0_adjusted, window_width, window_height);
}

/**
 * @brief 在搜索窗口覆盖的多个行中生成候选位置（考虑空位）
 * @param cell 当前单元
 * @param window 搜索窗口
 * @param num_candidates 候选位置数量
 * @param step_x X方向步长（站点间距）
 * @return 候选位置列表
 */
std::vector<Pointdp3D> LSP3D::generateCandidatePositions(
    Node3D* cell,
    const SearchWindow& window,
    int num_candidates,
    DBU step_x) {
    
    std::vector<Pointdp3D> candidates;
    
    if (!cell || !arch_ || !network_ || window.width <= 0 || num_candidates <= 0) {
        return candidates;
    }
    
    // 获取单元所在的层
    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    const auto& all_rows = arch_->getRows(layer_type);
    
    if (all_rows.empty()) {
        return candidates;
    }
    
    // 1. 找到搜索窗口覆盖的行
    std::vector<Architecture3D::Row3D*> covered_rows;
    
    for (auto* row : all_rows) {
        if (!row) continue;
        
        // 检查行是否与搜索窗口重叠
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;
        
        // Y方向重叠检查
        bool y_overlap = (row_bottom <= window.y + window.height) && 
                        (row_top >= window.y);
        
        // // X方向重叠检查
        // bool x_overlap = (row_left <= window.x + window.width) && 
        //                 (row_right >= window.x);
        
        if (y_overlap) {
            covered_rows.push_back(row);
        }
    }
    
    if (covered_rows.empty()) {
        return candidates;
    }

    // 方案A：先枚举窗口内所有合法落脚点，再均匀降采样到 num_candidates。
    // 行从下到上排序，保证枚举顺序（行主序：下→上、行内左→右）稳定。
    std::sort(covered_rows.begin(), covered_rows.end(),
              [](Architecture3D::Row3D* a, Architecture3D::Row3D* b) {
                  return a->getBottom().v < b->getBottom().v;
              });

    DBU cell_width = cell->getWidth().v;
    DBU cell_height = cell->getHeight().v;

    // 第一遍：全枚举所有合法（站点对齐、避开占用、落在窗口内）的候选位置
    struct Cand { DBU x; DBU y; DbuZ z; };
    std::vector<Cand> all;

    for (auto* row : covered_rows) {
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;

        // 检查单元高度是否适合行高度
        if (cell_height > (row_top - row_bottom)) {
            continue;  // 单元太高，不适合这个行
        }

        // 计算行与搜索窗口的交集
        DBU intersect_left = std::max(window.x, row_left);
        DBU intersect_right = std::min(window.x + window.width, row_right);
        if (intersect_right - intersect_left <= 0) {
            continue;  // 没有交集
        }

        DBU site_spacing = row->getSiteSpacing().v;
        DBU site_width = row->getSiteWidth().v;

        // 细步长：优先用调用方传入的 step_x，否则用站点间距，做非零兜底
        DBU stride = (step_x != 0) ? step_x : site_spacing;
        if (stride == 0) stride = site_width;
        if (stride == 0) stride = 100;

        // 找到行中的可用空位（站点对齐、避开已占用单元）
        std::vector<std::pair<DBU, DBU>> available_slots =
            findAvailableSlotsInRowAligned(row, cell, intersect_left, intersect_right, site_width, site_spacing);

        for (const auto& slot : available_slots) {
            DBU slot_start = slot.first;
            DBU slot_end = slot.second;

            if (slot_end - slot_start < cell_width) {
                continue;  // 空位太小，放不下单元
            }

            // 在空位内以固定步长从左到右枚举所有合法位置
            for (DBU x = slot_start; x + cell_width <= slot_end; x += stride) {
                // 确保位置在搜索窗口内
                if (x < window.x || x + cell_width > window.x + window.width) {
                    continue;
                }
                all.push_back({x, row_bottom, cell->location.z});
            }
        }
    }

    // 第二遍：均匀降采样。M<=N 全给；M>N 按等距下标抽 N 个（保序，均匀覆盖）。
    const int M = static_cast<int>(all.size());
    if (M == 0) {
        return candidates;
    }

    if (M <= num_candidates) {
        candidates.reserve(M);
        for (const auto& c : all) {
            candidates.emplace_back(DbuX{c.x}, DbuY{c.y}, c.z);
        }
        return candidates;
    }

    candidates.reserve(num_candidates);
    const double stepf = static_cast<double>(M) / num_candidates;
    for (int k = 0; k < num_candidates; ++k) {
        int idx = static_cast<int>(std::floor(k * stepf));
        if (idx >= M) idx = M - 1;
        candidates.emplace_back(DbuX{all[idx].x}, DbuY{all[idx].y}, all[idx].z);
    }

    return candidates;
}

std::vector<Pointdp3D> LSP3D::generateCandidatePositions2(
    Node3D* cell,
    const SearchWindow& window,
    int num_candidates)
{
    std::vector<Pointdp3D> candidates;
    if (num_candidates <= 0 || !cell || !arch_) return candidates;

    DBU W = window.width;
    DBU H = window.height;

    DBU cell_h = cell->getHeight().v;
    int num_y = (cell_h > 0) ? std::max(1, (int)(H / cell_h)) : 1;
    num_y = std::max(1, std::min(num_y, num_candidates));
    int num_x = std::max(1, num_candidates / num_y);

    DBU step_x = W / num_x;
    DBU step_y = cell_h;

    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    const auto& rows = arch_->getRows(layer_type);
    if (rows.empty()) return candidates;

    // 锚定到离窗口底部 (window.y) 最近的行
    Architecture3D::Row3D* base_row = nullptr;
    DBU base_dist = std::numeric_limits<DBU>::max();
    for (auto* row : rows) {
        if (!row) continue;
        DBU d = std::abs(window.y - row->getBottom().v);
        if (d < base_dist) { base_dist = d; base_row = row; }
    }
    if (!base_row) return candidates;
    DBU y0 = base_row->getBottom().v;

    // 从最底行逐行向上（step_y = cell_h），跳过重复行，凑够 num_y 行
    int rows_used = 0;
    bool have_last = false;
    DBU last_snapped_y = 0;
    for (int step = 0;
         rows_used < num_y && (int)candidates.size() < num_candidates;
         step++) {
        DBU y = y0 + step * step_y;
        if (y > window.y + H) break;  // 越过窗口顶部则停止

        Architecture3D::Row3D* nearest_row = nullptr;
        DBU min_dist = std::numeric_limits<DBU>::max();
        for (auto* row : rows) {
            if (!row) continue;
            DBU row_bottom = row->getBottom().v;
            DBU dist = std::abs(y - row_bottom);
            if (dist < min_dist) {
                min_dist = dist;
                nearest_row = row;
            }
        }
        if (!nearest_row) continue;
        DBU snapped_y = nearest_row->getBottom().v;

        // 跳过重复行（缺行/宏遮挡导致 snap 到同一行时继续向上）
        if (have_last && snapped_y == last_snapped_y) continue;
        have_last = true;
        last_snapped_y = snapped_y;
        rows_used++;

        DBU site_spacing = nearest_row->getSiteSpacing().v;
        DBU row_left = nearest_row->getLeft().v;
        DBU row_right = nearest_row->getRight().v;

        for (int i = 0; i < num_x && (int)candidates.size() < num_candidates; i++) {
            DBU x = window.x + i * step_x;

            x = alignToSite(x, site_spacing, /*round_up=*/false);
            x = std::max(row_left, x);
            x = std::min(row_right, x);

            candidates.emplace_back(DbuX{x}, DbuY{snapped_y}, DbuZ{0});
        }
    }
    return candidates;
}

/**
 * @brief 时序驱动的搜索窗口定位（主函数）
 * @param cell 要优化的单元
 * @param num_candidates 候选位置数量
 * @param window_width 搜索窗口宽度（如果为0则使用默认值）
 * @param window_height 搜索窗口高度（如果为0则使用默认值）
 * @return 候选位置列表
 */
std::vector<Pointdp3D> LSP3D::timingDrivenWindowSearch(
    Node3D* cell,
    int num_candidates,
    int window_x,
    int window_y){
    
    std::vector<Pointdp3D> candidates;
    DBU window_width, window_height;

    if (!cell || !arch_ || num_candidates <= 0) {
        return candidates;
    }
    
    // 确定搜索窗口尺寸
    if (window_x == 0 || window_y == 0) {
        // 使用默认窗口尺寸
        window_width = cell->getWidth().v * 5;
        window_height = cell->getHeight().v * 3;
    } else{
        window_width = cell->getWidth().v * window_x;
        window_height = cell->getHeight().v * window_y;
    }
    
    // 确定搜索窗口位置
    SearchWindow window = determineSearchWindow(cell, window_width, window_height);
    
    // 确保窗口在合法范围内
    window = clampSearchWindow(window, cell);
    
    // 生成候选位置
    candidates = generateCandidatePositions(cell, window, num_candidates);
    // candidates = generateCandidatePositions2(cell, window, num_candidates);
    
    // 为每个候选位置设置正确的层信息
    std::string cell_layer = cell->layer;
    DBU z_offset = (cell_layer == "bottom") ? 0 : 1000;  // 简化：底层Z=0，上层Z=1000
    
    for (auto& candidate : candidates) {
        candidate.z = DbuZ{z_offset};
    }
    
    return candidates;
}

/**
 * @brief 确保搜索窗口在合法范围内
 */
LSP3D::SearchWindow LSP3D::clampSearchWindow(
    const SearchWindow& window,
    Node3D* cell) {
    
    if (!cell || !arch_) {
        return window;
    }
    
    SearchWindow clamped = window;
    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    
    // 获取层的边界
    DbuX layer_min_x = arch_->getMinX(layer_type);
    DbuX layer_max_x = arch_->getMaxX(layer_type);
    DbuY layer_min_y = arch_->getMinY(layer_type);
    DbuY layer_max_y = arch_->getMaxY(layer_type);
    
    // 调整X方向
    if (clamped.x < layer_min_x.v) {
        clamped.x = layer_min_x.v;
    }
    
    if (clamped.x + clamped.width > layer_max_x.v) {
        clamped.x = layer_max_x.v - clamped.width;
        if (clamped.x < layer_min_x.v) {
            clamped.x = layer_min_x.v;
            clamped.width = layer_max_x.v - layer_min_x.v;
        }
    }
    
    // 调整Y方向
    if (clamped.y < layer_min_y.v) {
        clamped.y = layer_min_y.v;
    }
    
    if (clamped.y + clamped.height > layer_max_y.v) {
        clamped.y = layer_max_y.v - clamped.height;
        if (clamped.y < layer_min_y.v) {
            clamped.y = layer_min_y.v;
            clamped.height = layer_max_y.v - layer_min_y.v;
        }
    }
    
    return clamped;
}


// ============================================================================
// Row occupancy index — replaces the O(all_cells) scan in findAvailableSlotsInRow
// with O(cells in the queried row). Bucketed by geometric row index computed
// from Y (getRows() vector order is not guaranteed Y-sorted, so we key both the
// cell and the query row through rowIdxForY for consistency).
// ============================================================================
int LSP3D::rowIdxForY(int layer, DBU y_bottom) const {
    if (!arch_ || layer < 0 || layer > 1) return -1;
    auto lt = static_cast<Architecture3D::LayerType>(layer);
    const auto& rows = arch_->getRows(lt);
    if (rows.empty()) return -1;
    DBU min_y = arch_->getMinY(lt).v;
    DBU rh = rows[0]->getHeight().v;
    if (rh <= 0) return -1;
    long idx = (static_cast<long>(y_bottom) - static_cast<long>(min_y)) / static_cast<long>(rh);
    if (idx < 0 || idx >= static_cast<long>(rows.size())) return -1;
    return static_cast<int>(idx);
}

void LSP3D::buildRowIndexInto(RowIndex& dst) const {
    // NOTE: any code path that moves a cell (changes cell->location via
    // setInstanceLocation3D) MUST call updateRowIndexMove(cell, old, new) to keep
    // this index consistent, since the index is built once and then maintained
    // incrementally. Currently hooked: discreteLocalSearch / NoCommit / Parallel
    // net-commit sites. If MIS or any new mover is enabled, hook it too.
    dst.clear();
    if (!arch_ || !network_) return;
    dst.resize(2);
    for (int L = 0; L < 2; ++L) {
        auto lt = static_cast<Architecture3D::LayerType>(L);
        dst[L].assign(arch_->getRows(lt).size(), std::vector<Node3D*>());
    }
    for (auto* cell : network_->getCells()) {
        if (!cell) continue;
        int L = static_cast<int>(Architecture3D::stringToLayer(cell->layer));
        if (L < 0 || L > 1) continue;
        int ri = rowIdxForY(L, cell->location.y.v);
        if (ri < 0) continue;
        dst[L][ri].push_back(cell);
    }
}

void LSP3D::updateRowIndexMove(Node3D* cell,
                               const Pointdp3D& old_loc,
                               const Pointdp3D& new_loc) {
    if (row_index_.empty() || !cell || !arch_) return;
    // DLS candidates stay on the cell's own layer; use its current layer string.
    int L = static_cast<int>(Architecture3D::stringToLayer(cell->layer));
    if (L < 0 || L >= static_cast<int>(row_index_.size())) return;
    int old_ri = rowIdxForY(L, old_loc.y.v);
    int new_ri = rowIdxForY(L, new_loc.y.v);
    if (old_ri == new_ri) return;
    if (old_ri >= 0 && old_ri < static_cast<int>(row_index_[L].size())) {
        auto& b = row_index_[L][old_ri];
        for (size_t i = 0; i < b.size(); ++i)
            if (b[i] == cell) { b[i] = b.back(); b.pop_back(); break; }
    }
    if (new_ri >= 0 && new_ri < static_cast<int>(row_index_[L].size()))
        row_index_[L][new_ri].push_back(cell);
}


/**
 * @brief 在行中找到可用的空位
 * @param row 要检查的行
 * @param exclude_cell 要排除的单元（当前单元）
 * @param range_start 检查范围的起始X坐标
 * @param range_end 检查范围的结束X坐标
 * @return 可用空位列表（每个空位是[start, end]对）
 */
std::vector<std::pair<DBU, DBU>> LSP3D::findAvailableSlotsInRow(
    Architecture3D::Row3D* row,
    Node3D* exclude_cell,
    DBU range_start,
    DBU range_end) {
    
    std::vector<std::pair<DBU, DBU>> available_slots;
    
    if (!row || !network_) {
        return available_slots;
    }
    
    // 获取行参数
    DBU row_left = row->getLeft().v;
    DBU row_right = row->getRight().v;
    DBU row_bottom = row->getBottom().v;
    DBU row_top = row->getTop().v;
    
    // 调整检查范围，确保在行内
    range_start = std::max(range_start, row_left);
    range_end = std::min(range_end, row_right);
    
    if (range_start >= range_end) {
        return available_slots;
    }
    
    // 1. 获取行中所有单元（排除当前单元）
    std::vector<Node3D*> row_cells;

    // 用行占用索引取该行的 cell（O(行内cell)）；索引未建时回退到全设计扫描。
    int _L = static_cast<int>(row->getLayer());
    int _ri = rowIdxForY(_L, row_bottom);
    const std::vector<Node3D*>* bucket = nullptr;
    if (!row_index_.empty() && _L >= 0 && _L < static_cast<int>(row_index_.size())
        && _ri >= 0 && _ri < static_cast<int>(row_index_[_L].size())) {
        bucket = &row_index_[_L][_ri];
    }
    static std::vector<Node3D*> s_all_fallback;
    const std::vector<Node3D*>& scan_cells =
        bucket ? *bucket : network_->getCells();

    // 过滤出在当前行内的单元（保留与原实现完全一致的 Y/X 判定）
    for (auto* cell : scan_cells) {
        if (!cell || cell == exclude_cell) {
            continue;
        }
        
        // 检查单元是否在当前行内
        DBU cell_bottom = cell->location.y.v;
        DBU cell_top = cell_bottom + cell->getHeight().v;
        
        // Y方向检查：单元是否在当前行内
        bool y_in_row = (cell_bottom >= row_bottom) && (cell_top <= row_top);
        
        if (!y_in_row) {
            continue;  // 单元不在当前行内
        }
        
        // 检查单元是否在检查范围内（X方向）
        DBU cell_left = cell->location.x.v;
        DBU cell_right = cell_left + cell->getWidth().v;
        
        // 检查是否有重叠
        bool x_overlap = (cell_right >= range_start) && (cell_left <= range_end);
        
        if (x_overlap) {
            row_cells.push_back(cell);
        }
    }
    
    // 2. 按X坐标排序单元
    std::sort(row_cells.begin(), row_cells.end(),
              [](Node3D* a, Node3D* b) {
                  return a->location.x.v < b->location.x.v;
              });
    
    // 3. 找到占用区域
    std::vector<std::pair<DBU, DBU>> occupied_ranges;
    
    for (auto* cell : row_cells) {
        DBU cell_left = cell->location.x.v;
        DBU cell_right = cell_left + cell->getWidth().v;
        
        // 只考虑在检查范围内的单元
        if (cell_right <= range_start || cell_left >= range_end) {
            continue;
        }
        
        // 调整到检查范围内
        cell_left = std::max(cell_left, range_start);
        cell_right = std::min(cell_right, range_end);
        
        if (cell_left < cell_right) {
            occupied_ranges.emplace_back(cell_left, cell_right);
        }
    }
    
    // 4. 合并重叠的占用区域
    std::vector<std::pair<DBU, DBU>> merged_ranges;
    
    if (!occupied_ranges.empty()) {
        std::sort(occupied_ranges.begin(), occupied_ranges.end());
        
        DBU current_start = occupied_ranges[0].first;
        DBU current_end = occupied_ranges[0].second;
        
        for (size_t i = 1; i < occupied_ranges.size(); ++i) {
            if (occupied_ranges[i].first <= current_end) {
                // 重叠，合并
                current_end = std::max(current_end, occupied_ranges[i].second);
            } else {
                // 不重叠，保存当前区域
                merged_ranges.emplace_back(current_start, current_end);
                current_start = occupied_ranges[i].first;
                current_end = occupied_ranges[i].second;
            }
        }
        
        // 保存最后一个区域
        merged_ranges.emplace_back(current_start, current_end);
    }
    
    // 5. 从占用区域推导出可用空位
    DBU current_pos = range_start;
    
    for (const auto& occupied : merged_ranges) {
        if (current_pos < occupied.first) {
            // 找到可用空位
            available_slots.emplace_back(current_pos, occupied.first);
        }
        current_pos = std::max(current_pos, occupied.second);
    }
    
    // 检查最后一个区域之后的空位
    if (current_pos < range_end) {
        available_slots.emplace_back(current_pos, range_end);
    }
    
    return available_slots;
}

/**
 * @brief 在行中找到可用的空位（考虑站点对齐）
 * @param row 要检查的行
 * @param exclude_cell 要排除的单元（当前单元）
 * @param range_start 检查范围的起始X坐标
 * @param range_end 检查范围的结束X坐标
 * @param site_width 站点宽度
 * @param site_spacing 站点间距
 * @return 可用空位列表（每个空位是[start, end]对，已对齐到站点）
 */
std::vector<std::pair<DBU, DBU>> LSP3D::findAvailableSlotsInRowAligned(
    Architecture3D::Row3D* row,
    Node3D* exclude_cell,
    DBU range_start,
    DBU range_end,
    DBU site_width,
    DBU site_spacing) {
    
    std::vector<std::pair<DBU, DBU>> available_slots;
    
    if (!row || !network_) {
        return available_slots;
    }
    
    // 获取站点参数
    if (site_spacing == 0) {
        site_spacing = row->getSiteSpacing().v;
    }
    if (site_width == 0) {
        site_width = row->getSiteWidth().v;
    }
    
    if (site_spacing == 0) {
        site_spacing = 100;  // 默认值
    }
    
    // 1. 先找到原始的空位
    auto raw_slots = findAvailableSlotsInRow(row, exclude_cell, range_start, range_end);
    
    // 2. 对每个空位进行站点对齐
    for (const auto& raw_slot : raw_slots) {
        DBU slot_start = raw_slot.first;
        DBU slot_end = raw_slot.second;
        // DBU slot_width = slot_end - slot_start;
        
        // 对齐到站点边界
        DBU aligned_start = alignToSite(slot_start, site_spacing, true);  // 向上对齐
        DBU aligned_end = alignToSite(slot_end, site_spacing, false);     // 向下对齐
        
        // 确保对齐后的空位有效
        if (aligned_start < aligned_end) {
            available_slots.emplace_back(aligned_start, aligned_end);
        }
    }
    
    return available_slots;
}

/**
 * @brief 将坐标对齐到站点边界
 * @param coord 要对齐的坐标
 * @param site_spacing 站点间距
 * @param round_up 是否向上对齐（true=向上，false=向下）
 * @return 对齐后的坐标
 */
DBU LSP3D::alignToSite(DBU coord, DBU site_spacing, bool round_up) const {
    if (site_spacing == 0) {
        return coord;
    }
    
    DBU remainder = coord % site_spacing;
    
    if (remainder == 0) {
        return coord;  // 已经在站点边界上
    }
    
    if (round_up) {
        return coord + (site_spacing - remainder);
    } else {
        return coord - remainder;
    }
}

// ============================================================================
// 翻转方法
// ============================================================================

/**
 * @brief 翻转单元（使用MX方向）
 * @param cell 要翻转的单元
 */
void LSP3D::flipCell(Node3D* cell) const {
    if (!cell) {
        return;
    }
    
    // 获取当前方向
    odb::dbOrientType current_orient = cell->getOrient();
    
    // 根据当前方向计算翻转后的方向
    odb::dbOrientType new_orient = current_orient;
    
    // 应用MX翻转（关于X轴镜像）
    switch (current_orient) {
        case odb::dbOrientType::R0:
            new_orient = odb::dbOrientType::MX;
            break;
        case odb::dbOrientType::R90:
            new_orient = odb::dbOrientType::MXR90;
            break;
        case odb::dbOrientType::R180:
            new_orient = odb::dbOrientType::R0;  // MX + R180 = R0
            break;
        case odb::dbOrientType::R270:
            new_orient = odb::dbOrientType::MYR90;  // MX + R270 = MYR90
            break;
        case odb::dbOrientType::MX:
            new_orient = odb::dbOrientType::R0;  // 两次MX = R0
            break;
        case odb::dbOrientType::MY:
            new_orient = odb::dbOrientType::R180;  // MX + MY = R180
            break;
        case odb::dbOrientType::MYR90:
            new_orient = odb::dbOrientType::R270;  // MX + MYR90 = R270
            break;
        case odb::dbOrientType::MXR90:
            new_orient = odb::dbOrientType::R90;  // MX + MXR90 = R90
            break;
        default:
            new_orient = odb::dbOrientType::MX;  // 默认使用MX
            break;
    }
    
    // 设置新方向
    cell->setOrient(new_orient);
    
    // std::cout << "  Cell " << cell->name << " flipped from " 
    //           << orientToString(current_orient) << " to " 
    //           << orientToString(new_orient) << std::endl;
}

/**
 * @brief 将方向转换为字符串
 */
std::string LSP3D::orientToString(const odb::dbOrientType& orient) const {
    switch (orient) {
        case odb::dbOrientType::R0: return "R0";
        case odb::dbOrientType::R90: return "R90";
        case odb::dbOrientType::R180: return "R180";
        case odb::dbOrientType::R270: return "R270";
        case odb::dbOrientType::MY: return "MY";
        case odb::dbOrientType::MYR90: return "MYR90";
        case odb::dbOrientType::MX: return "MX";
        case odb::dbOrientType::MXR90: return "MXR90";
        default: return "UNKNOWN";
    }
}



// ============================================================================
// 成本计算方法
// ============================================================================
double LSP3D::calculateCellCost(Node3D* cell, double alpha, double omega,
                               TimingEngine3D* timing_engine) const {
    if (!cell || !network_) {
        return 0.0;
    }
    
    //计算线长成本
    // double wirelength_cost = calculateCellHPWL(cell);
    
    //计算时序成本
    double timing_cost = calculateCellTimingCost(cell, alpha, omega, timing_engine_);
    
    //组合成本（可以根据需要调整权重）
    //double total_cost = wirelength_cost + timing_cost;
    
    return timing_cost;
}

std::vector<Pin3D*> LSP3D::getCellInputPins(Node3D* cell) const {
    std::vector<Pin3D*> input_pins;
    
    if (!cell || !network_) {
        return input_pins;
    }
    
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_IN) {
            input_pins.push_back(pin);
        }
    }
    
    return input_pins;
}

std::vector<Pin3D*> LSP3D::getCellOutputPins(Node3D* cell) const {
    std::vector<Pin3D*> output_pins;
    
    if (!cell || !network_) {
        return output_pins;
    }
    
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_OUT) {
            output_pins.push_back(pin);
        }
    }
    
    return output_pins;
}

void LSP3D::cellArcDelays(Node3D* cell) const {
    if (!network_ || !timing_engine_) return;
    if (!cell || cell->is_fixed) return;
    auto output_pins = getCellOutputPins(cell);
    auto input_pin = getCellInputPins(cell);
    if (input_pin.empty()) return;

    for (auto* out_pin : output_pins) {
        double min_early = 0.0 , max_late = 0.0;
        for (const auto& in_name : input_pin) {
            auto info = timing_engine_->getArcDelay(in_name->getName(), out_pin->getName());
            if (info.is_valid) {
                min_early = std::min(min_early, info.early_delay);
                max_late  = std::max(max_late,  info.late_delay);
            }
        }
        out_pin->setCachedArcDelay(min_early, max_late);
    }

    for (auto* in_pin : input_pin) {
        double min_early = 0.0 , max_late = 0.0;
        for (const auto& out_pin : output_pins) {
            auto info = timing_engine_->getArcDelay(in_pin->getName(), out_pin->getName());
            if (info.is_valid) {
                min_early = std::min(min_early, info.early_delay);
                max_late  = std::max(max_late,  info.late_delay);
            }
        }
        in_pin->setCachedArcDelay(min_early, max_late);
    }
    // std::cout << "Cached arc delays for " << cached << " output pins" << std::endl;
}
/**
 * @brief 计算单元的时序成本（包含单元内部弧延迟）
 * @param cell 要计算时序成本的单元
 * @param alpha 晚时序权重
 * @param omega 早时序权重
 * @param timing_engine 时序引擎
 * @return 时序成本
 */
double LSP3D::calculateCellTimingCost(Node3D* cell, double alpha, double omega, 
                                     TimingEngine3D* timing_engine) const {
    if (!cell || !network_) {
        return 0.0;
    }
    
    double cost_late = 0.0;
    double cost_early = 0.0;
    
    cellArcDelays(cell);
    // //对每个sink cell的arc进行更新
    std::unordered_set<Node3D*> visited;
    auto out_pins = getCellOutputPins(cell);
    for (auto* out_pin : out_pins) {
        Net3D* net = out_pin->getNet();
        if (!net) continue;
        for (auto* sink_pin : net->getPins()) {
            if (sink_pin == out_pin) continue;
            Node3D* sink_cell = sink_pin->getNode();
            if (sink_cell && visited.insert(sink_cell).second) {
                cellArcDelays(sink_cell);
            }
        }
    }
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        // 获取引脚方向
        int dir_int = pin->getDirection();
        Pin3D::Direction dir = static_cast<Pin3D::Direction>(dir_int);

        Net3D* net = pin->getNet();
        if (net && net->db_net && (timingdp_->isIdealClockNet(net->db_net) || !timingdp_->isValidSigType(net->db_net) || net->db_net->isSpecial())) {
            continue;
        }
        
        if (dir == Pin3D::Dir_IN) {
            // 输入引脚：计算从驱动引脚到当前引脚的延迟
            // Net3D* net = pin->getNet();
            if (!net) continue;
            
            // 找到驱动引脚（source）
            Pin3D* driver_pin = findDriverPin(net, pin);
            if (!driver_pin) continue;

            // 获取单元的内部弧延迟
            double arc_delay_early = 0.0;
            double arc_delay_late = 0.0;
            
            if (pin->isArcCached()) {
                arc_delay_early = pin->getCachedArcEarly();
                arc_delay_late  = pin->getCachedArcLate();
            }
            
            // 计算延迟（线网延迟 + 驱动单元的内部弧延迟）
            double delay_early = pin->getEarlyArrival() - driver_pin->getEarlyArrival() + arc_delay_early;
            double delay_late = pin->getLateArrival() - driver_pin->getLateArrival() + arc_delay_late;
            // double delay_early = pin->getEarlyArrival() - driver_pin->getEarlyArrival();
            // double delay_late = pin->getLateArrival() - driver_pin->getLateArrival();
            // 获取权重
            double early_weight = pin->getEarlyWeight();
            double late_weight = pin->getLateWeight();
            
            cost_early += delay_early * early_weight;
            cost_late += delay_late * late_weight;
            
        } else if (dir == Pin3D::Dir_OUT) {
            // 输出引脚：计算到所有接收引脚的延迟
            // Net3D* net = pin->getNet();
            if (!net) continue;
            
            auto net_pins = net->getPins();
            
            for (auto* sink_pin : net_pins) {
                if (sink_pin == pin) continue;  // 跳过驱动引脚
                
                // 计算线网延迟
                double net_delay_early = sink_pin->getEarlyArrival() - pin->getEarlyArrival();
                double net_delay_late = sink_pin->getLateArrival() - pin->getLateArrival();
                
                // 获取接收单元的内部弧延迟
                double arc_delay_early = 0.0;
                double arc_delay_late = 0.0;

                if(!isSequentialCell(sink_pin->getNode())){
                    arc_delay_early = sink_pin->getCachedArcEarly();
                    arc_delay_late = sink_pin->getCachedArcLate();
                }

            
                // 总延迟 = 线网延迟 + 单元内部弧延迟
                double total_delay_early = net_delay_early + arc_delay_early;
                double total_delay_late = net_delay_late + arc_delay_late;
                
                // 获取权重
                double early_weight = sink_pin->getEarlyWeight();
                double late_weight = sink_pin->getLateWeight();
                
                cost_early += total_delay_early * early_weight;
                cost_late += total_delay_late * late_weight;
                
            }
        }
    }
    
    // 组合时序成本
    double total_timing_cost = alpha * cost_late - omega * cost_early;

    // ---- cost overflow: per-pin delay trace ----
    if (std::abs(total_timing_cost) > 1e15) {
        static int dbg_cnt = 0;
        if (dbg_cnt < 20) {
            logger_->warn(TDP, 2097, "[COST-TRC] cell={} cost_late={:.6e} cost_early={:.6e} total={:.6e}",
                cell->name, cost_late, cost_early, total_timing_cost);
            for (auto* pin : cell_pins) {
                if (!pin) continue;
                Pin3D::Direction d = static_cast<Pin3D::Direction>(pin->getDirection());
                Net3D* n = pin->getNet();
                if (d == Pin3D::Dir_IN) {
                    Pin3D* dp = findDriverPin(n, pin);
                    double ae = pin->getEarlyArrival(), al = pin->getLateArrival();
                    double de = dp ? dp->getEarlyArrival() : 0, dl = dp ? dp->getLateArrival() : 0;
                    double ee = 0, el = 0;
                    if (pin->isArcCached()) { ee = pin->getCachedArcEarly(); el = pin->getCachedArcLate(); }
                    logger_->warn(TDP, 2098, "[COST-TRC]   IN pin={} arrE={:.6e} arrL={:.6e} drvE={:.6e} drvL={:.6e} arcE={:.6e} arcL={:.6e} wE={:.6e} wL={:.6e}",
                        pin->getName(), ae, al, de, dl, ee, el, pin->getEarlyWeight(), pin->getLateWeight());
                } else if (d == Pin3D::Dir_OUT && n) {
                    for (auto* sp : n->getPins()) {
                        if (sp == pin) continue;
                        double ne = sp->getEarlyArrival() - pin->getEarlyArrival();
                        double nl = sp->getLateArrival() - pin->getLateArrival();
                        double ae = sp->getCachedArcEarly(), al = sp->getCachedArcLate();
                        logger_->warn(TDP, 2099, "[COST-TRC]   OUT sink={} netE={:.6e} netL={:.6e} arcE={:.6e} arcL={:.6e} wE={:.6e} wL={:.6e}",
                            sp->getName(), ne, nl, ae, al, sp->getEarlyWeight(), sp->getLateWeight());
                    }
                }
            }
            ++dbg_cnt;
        }
    }

    return total_timing_cost;
}

// ============================================================================
// Phase-2 parallel: per-thread pool, thread-local arrival map + cost
// ============================================================================
void LSP3D::ensureParallelPool(int num_threads)
{
    if (num_threads < 1) num_threads = 1;
    if (par_nthreads_ >= num_threads && !par_engines_.empty()) return;

    par_engines_.clear();
    par_steiners_.clear();
    par_local_par_.assign(num_threads, nullptr);
    par_engines_.reserve(num_threads);
    par_steiners_.reserve(num_threads);

    // a valid instance to warm up ensureInit()/local_sta_ serially per engine
    sta::Instance* warm = nullptr;
    sta::dbNetwork* db_network = sta_->getDbNetwork();
    for (auto* c : network_->getCells()) {
        if (c && !c->is_fixed && c->db_inst) {
            warm = db_network->dbToSta(c->db_inst);
            if (warm) break;
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        par_engines_.push_back(std::make_unique<lgista::LocalGraphIsta>(sta_));
        par_steiners_.push_back(std::make_unique<SteinerTreeGenerator>(timingdp_));
        if (warm) { auto pg = par_engines_[t]->analyzeOnly(warm, 2); (void)pg; }
        par_local_par_[t] = par_engines_[t]->localParasiticsContainer();
    }
    par_nthreads_ = num_threads;
}

void LSP3D::buildLocalArrivalMap(
    Node3D* cell, lgista::PtGraph* pt_graph,
    std::unordered_map<Pin3D*, std::pair<double, double>>& out) const
{
    if (!cell || !pt_graph || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    sta::Graph* graph = sta_ ? sta_->graph() : nullptr;
    if (!network || !graph) return;

    auto add = [&](Pin3D* p) {
        if (!p || out.count(p)) return;
        sta::Pin* sp = network->findPin(p->getName().c_str());
        if (!sp) return;
        sta::Vertex *v = nullptr, *b = nullptr;
        graph->pinVertices(sp, v, b);
        if (!v) return;
        lgista::PtVertex* ptv = pt_graph->ptVertex(v);
        if (!ptv) return;
        double e = 0.0, l = 0.0;
        if (ptVertexLocalArrival(pt_graph, *ptv, sta_, e, l))
            out[p] = {e, l};
    };

    for (auto* p : cell->getPins()) add(p);
    for (auto* p : cell->getPins()) {
        if (!p) continue;
        Net3D* net = p->getNet();
        if (!net) continue;
        if (p->getDirection() == Pin3D::Dir_OUT) {
            for (auto* sp : net->getPins()) if (sp != p) add(sp);
        }
        // NOTE: fanin driver pins are intentionally NOT added here, to match
        // no-commit's refreshAffectedPinTimingFromPtGraph (which refreshes only
        // the cell + its fanout sinks, leaving drivers at the global cache).
        // calculateCellTimingCostLocal falls back to the global Pin3D arrival
        // for pins absent from this map, giving par the same driver reference
        // frame as no-commit.
    }
}

double LSP3D::calculateCellTimingCostLocal(
    Node3D* cell, double alpha, double omega,
    const std::unordered_map<Pin3D*, std::pair<double, double>>& arr) const
{
    static const std::unordered_map<Pin3D*, std::pair<double, double>> kEmptyArc;
    return calculateCellTimingCostLocal(cell, alpha, omega, arr, kEmptyArc);
}

// Per-input-pin internal gate arc delays from the (moved) local PtGraph.
// Mirrors cellArcDelays semantics (init 0.0, min for early / max for late), but
// sourced from the trial PtGraph so par sees the moved-geometry gate delay just
// like no-commit's live getArcDelay. Only cell input pins + fanout sink input
// pins are needed (same pins whose cached arc the cost uses).
void LSP3D::buildLocalArcMap(
    Node3D* cell, lgista::PtGraph* pt_graph,
    std::unordered_map<Pin3D*, std::pair<double, double>>& out) const
{
    if (!cell || !pt_graph || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    sta::Graph* graph = sta_ ? sta_->graph() : nullptr;
    sta::Scene* scene = pt_graph->scene();
    if (!network || !graph || !scene) return;
    const sta::DcalcAPIndex ap_min = scene->dcalcAnalysisPtIndex(sta::MinMax::min());
    const sta::DcalcAPIndex ap_max = scene->dcalcAnalysisPtIndex(sta::MinMax::max());

    auto add_arc = [&](Pin3D* p) {
        if (!p || p->getDirection() != Pin3D::Dir_IN || out.count(p)) return;
        sta::Pin* sp = network->findPin(p->getName().c_str());
        if (!sp) return;
        sta::Vertex *v = nullptr, *b = nullptr;
        graph->pinVertices(sp, v, b);
        if (!v) return;
        lgista::PtVertex* ptv = pt_graph->ptVertex(v);
        if (!ptv) return;
        bool have = false;
        double e_min = 0.0, l_max = 0.0;   // matches cellArcDelays init
        lgista::PtVertexOutEdgeIterator it(ptv->objectIdx(), pt_graph);
        while (it.hasNext()) {
            lgista::PtEdge& pe = it.next();
            if (!pe.hasBase() || pe.isWire()) continue;
            const sta::TimingRole* role = pe.role();
            if (role && role->isTimingCheck()) continue;
            sta::TimingArcSet* arc_set = pe.timingArcSet();
            if (!arc_set) continue;
            for (sta::TimingArc* arc : arc_set->arcs()) {
                double de = sta::delayAsFloat(pt_graph->arcDelay(pe, arc, ap_min));
                double dl = sta::delayAsFloat(pt_graph->arcDelay(pe, arc, ap_max));
                e_min = std::min(e_min, de);
                l_max = std::max(l_max, dl);
                have = true;
            }
        }
        if (have) out[p] = {e_min, l_max};
    };

    for (auto* p : cell->getPins()) {
        if (!p) continue;
        Net3D* net = p->getNet();
        if (!net) continue;
        if (p->getDirection() == Pin3D::Dir_IN) add_arc(p);
        else if (p->getDirection() == Pin3D::Dir_OUT)
            for (auto* sp : net->getPins()) if (sp != p) add_arc(sp);
    }
}

double LSP3D::calculateCellTimingCostLocal(
    Node3D* cell, double alpha, double omega,
    const std::unordered_map<Pin3D*, std::pair<double, double>>& arr,
    const std::unordered_map<Pin3D*, std::pair<double, double>>& arc_map) const
{
    if (!cell || !network_) return 0.0;

    auto early = [&](Pin3D* p) {
        auto it = arr.find(p);
        return it != arr.end() ? it->second.first : p->getEarlyArrival();
    };
    auto late = [&](Pin3D* p) {
        auto it = arr.find(p);
        return it != arr.end() ? it->second.second : p->getLateArrival();
    };
    // Arc delay: prefer the moved-geometry value from the PtGraph arc map;
    // fall back to the frozen cached arc when absent.
    auto arc_early = [&](Pin3D* p, double fb) {
        auto it = arc_map.find(p);
        return it != arc_map.end() ? it->second.first : fb;
    };
    auto arc_late = [&](Pin3D* p, double fb) {
        auto it = arc_map.find(p);
        return it != arc_map.end() ? it->second.second : fb;
    };

    double cost_late = 0.0, cost_early = 0.0;
    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        Pin3D::Direction dir = static_cast<Pin3D::Direction>(pin->getDirection());
        Net3D* net = pin->getNet();
        if (net && net->db_net && (timingdp_->isIdealClockNet(net->db_net) || !timingdp_->isValidSigType(net->db_net) || net->db_net->isSpecial())) continue;

        if (dir == Pin3D::Dir_IN) {
            if (!net) continue;
            Pin3D* driver_pin = findDriverPin(net, pin);
            if (!driver_pin) continue;
            double ae = 0.0, al = 0.0;
            if (pin->isArcCached()) { ae = pin->getCachedArcEarly(); al = pin->getCachedArcLate(); }
            ae = arc_early(pin, ae); al = arc_late(pin, al);
            double de = early(pin) - early(driver_pin) + ae;
            double dl = late(pin)  - late(driver_pin)  + al;
            cost_early += de * pin->getEarlyWeight();
            cost_late  += dl * pin->getLateWeight();
        } else if (dir == Pin3D::Dir_OUT) {
            if (!net) continue;
            for (auto* sink_pin : net->getPins()) {
                if (sink_pin == pin) continue;
                double nde = early(sink_pin) - early(pin);
                double ndl = late(sink_pin)  - late(pin);
                double ae = 0.0, al = 0.0;
                if (!isSequentialCell(sink_pin->getNode())) {
                    ae = sink_pin->getCachedArcEarly();
                    al = sink_pin->getCachedArcLate();
                }
                ae = arc_early(sink_pin, ae); al = arc_late(sink_pin, al);
                cost_early += (nde + ae) * sink_pin->getEarlyWeight();
                cost_late  += (ndl + al) * sink_pin->getLateWeight();
            }
        }
    }
    return alpha * cost_late - omega * cost_early;
}

bool LSP3D::discreteLocalSearchParallel(Node3D* cell, int search_windows,
                                        double alpha, double omega,
                                        TimingEngine3D* timing_engine)
{
    if (!cell || cell->is_fixed || !arch_ || !network_) return false;

    const bool prof = prof_enabled();
    double _t = prof ? now_ms() : 0.0;

    std::vector<Pointdp3D> candidates =
        timingDrivenWindowSearch(cell, search_windows, params_.window_x, params_.window_y);
    if (prof) g_prof.candgen += now_ms() - _t;
    if (candidates.empty()) return false;

    odb::dbInst* inst = cell->db_inst;
    sta::Instance* inst_sta = sta_->getDbNetwork()->dbToSta(inst);
    if (!inst_sta) return false;
    sta::Scene* corner = sta_->scenes().empty() ? nullptr : sta_->scenes()[0];
    if (!corner) return false;
    const int depth = 1;
    const Pointdp3D cur = cell->location;

    // Prime arc caches for cell + sinks (serial; read-only in the region).
    _t = prof ? now_ms() : 0.0;
    calculateCellCost(cell, alpha, omega, timing_engine_);
    if (prof) g_prof.prime += now_ms() - _t;

    const int N = static_cast<int>(candidates.size());
    int nthreads = std::max(1, std::min(omp_get_max_threads(), N));
    ensureParallelPool(nthreads);

    // Per-thread profiler accumulators for the parallel candidate eval.
    std::vector<double> tp(nthreads, 0.0), ta(nthreads, 0.0), tm(nthreads, 0.0), tc(nthreads, 0.0);
    std::vector<double> tmpt(nthreads, 0.0), trp(nthreads, 0.0), tdly(nthreads, 0.0),
                        tarr(nthreads, 0.0), treq(nthreads, 0.0);

    // Precompute, per connected net: the cell's pin on it + current odb coord.
    struct TNet { Net3D* net; sta::Net* sn; Pin3D* pin; int px; int py; };
    std::vector<TNet> tnets;
    SteinerTreeGenerator* st0 = par_steiners_[0].get();
    for (auto* dbnet : cell->connected_nets) {
        if (!dbnet) continue;
        Net3D* n3 = network_->getNet(dbnet);
        if (!n3) continue;
        sta::Net* sn = sta_->getDbNetwork()->dbToSta(dbnet);
        if (!sn) continue;
        Pin3D* cp = nullptr;
        for (auto* p : cell->getPins())
            if (p && p->getNet() == n3) { cp = p; break; }
        if (!cp) continue;
        Pin3Dlocation pl = st0->getPin3DLocation(cp);
        tnets.push_back({n3, sn, cp, static_cast<int>(pl.x), static_cast<int>(pl.y)});
    }

    // Evaluate a position on thread tid via its private engine/steiner/container.
    auto eval = [&](int tid, const Pointdp3D& pos) -> double {
        const long dx = static_cast<long>(pos.x.v) - static_cast<long>(cur.x.v);
        const long dy = static_cast<long>(pos.y.v) - static_cast<long>(cur.y.v);
        sta::Parasitics* lp = par_local_par_[tid];
        double s = prof ? now_ms() : 0.0;
        for (auto& tn : tnets)
            par_steiners_[tid]->buildLocalNetParasiticInto(
                lp, tn.net, corner, tn.pin,
                tn.px + static_cast<int>(dx), tn.py + static_cast<int>(dy));
        if (prof) { tp[tid] += now_ms() - s; s = now_ms(); }
        lgista::LocalGraphIsta::AnalyzeTimings atm;
        lgista::LocalGraphIsta::AnalyzeTimings* patm = prof ? &atm : nullptr;
        auto pg = par_engines_[tid]->analyzeOnly(inst_sta, depth, patm);
        if (prof) { ta[tid] += now_ms() - s; s = now_ms(); }
        if (prof) {
            tmpt[tid] += atm.make_ptgraph_ms;
            trp[tid]  += atm.recompute_para_ms;
            tdly[tid] += atm.delays_ms;
            tarr[tid] += atm.arrivals_ms;
            treq[tid] += atm.requireds_ms;
        }
        std::unordered_map<Pin3D*, std::pair<double, double>> arr;
        buildLocalArrivalMap(cell, pg.get(), arr);
        std::unordered_map<Pin3D*, std::pair<double, double>> arc_map;
        buildLocalArcMap(cell, pg.get(), arc_map);
        if (prof) { tm[tid] += now_ms() - s; s = now_ms(); }
        double c = calculateCellTimingCostLocal(cell, alpha, omega, arr, arc_map);
        // ---- overflow check for ALL cells in PAR eval ----
        if (std::abs(c) > 1e15) {
            static int ovf_cnt = 0;
            if (ovf_cnt < 20) {
                std::string cn = cell->db_inst ? cell->db_inst->getName() : "?";
                logger_->warn(TDP, 2108, "[PAR-OVF] cell={} cost={:.6e}", cn, c);
                for (auto* pin : cell->getPins()) {
                    if (!pin) continue;
                    auto it = arr.find(pin);
                    if (it != arr.end())
                        logger_->warn(TDP, 2109, "[PAR-OVF]   pin={} local_arrE={:.6e} local_arrL={:.6e}",
                            pin->getName(), it->second.first, it->second.second);
                }
                ++ovf_cnt;
            }
        }
        // ---- TDP_DIAG_PAR: check local arrival/arc for matching cells ----
        if (std::getenv("TDP_DIAG_PAR")) {
            static int par_dbg_cnt = 0;
            if (par_dbg_cnt < 30) {
                std::string nm = cell->db_inst ? cell->db_inst->getName() : "";
                if (nm.find("sbe_bp_predict_address") != std::string::npos) {
                    ++par_dbg_cnt;
                    logger_->info(TDP, 2104, "[PAR-DBG] cell={} cost={:.6e}",
                        nm, c);
                    for (auto* pin : cell->getPins()) {
                        if (!pin) continue;
                        auto it = arr.find(pin);
                        if (it != arr.end())
                            logger_->info(TDP, 2105, "[PAR-DBG]   pin={} local_arrE={:.6e} local_arrL={:.6e} glob_arrE={:.6e} glob_arrL={:.6e}",
                                pin->getName(), it->second.first, it->second.second,
                                pin->getEarlyArrival(), pin->getLateArrival());
                        else
                            logger_->info(TDP, 2107, "[PAR-DBG]   pin={} NO_LOCAL glob_arrE={:.6e} glob_arrL={:.6e}",
                                pin->getName(), pin->getEarlyArrival(), pin->getLateArrival());
                    }
                    for (auto* pin : cell->getPins()) {
                        if (!pin || pin->getDirection() != Pin3D::Dir_IN) continue;
                        auto it = arc_map.find(pin);
                        double ae = pin->getCachedArcEarly(), al = pin->getCachedArcLate();
                        if (it != arc_map.end())
                            logger_->info(TDP, 2106, "[PAR-DBG]   arc pin={} local_arcE={:.6e} local_arcL={:.6e} cached_arcE={:.6e} cached_arcL={:.6e}",
                                pin->getName(), it->second.first, it->second.second, ae, al);
                    }
                }
            }
        }
        if (prof) { tc[tid] += now_ms() - s; }
        for (auto& tn : tnets) lp->deleteParasiticNetwork(tn.sn);
        return c;
    };

    // Baseline (current) cost via the real global calculateCellCost, so the
    // accept threshold matches the (possibly commit-polluted) global state the
    // move is actually judged against — same baseline no-commit uses.
    _t = prof ? now_ms() : 0.0;
    const double current_cost = calculateCellCost(cell, alpha, omega, timing_engine_);
    if (prof) g_prof.baseline += now_ms() - _t;

    std::vector<double> costs(N, std::numeric_limits<double>::max());
    double _wall = prof ? now_ms() : 0.0;
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int i = 0; i < N; ++i) {
        int tid = omp_get_thread_num();
        costs[i] = eval(tid, candidates[i]);
    }
    if (prof) {
        g_prof.par_wall += now_ms() - _wall;
        for (int t = 0; t < nthreads; ++t) {
            g_prof.par_parasitic += tp[t];
            g_prof.par_analyze   += ta[t];
            g_prof.par_maps      += tm[t];
            g_prof.par_cost      += tc[t];
            g_prof.par_ana_makept   += tmpt[t];
            g_prof.par_ana_recomp   += trp[t];
            g_prof.par_ana_delays   += tdly[t];
            g_prof.par_ana_arrivals += tarr[t];
            g_prof.par_ana_requireds += treq[t];
        }
        g_prof.cells++;
        g_prof.candidates += N;
    }

    // Deterministic argmin (lowest index on ties).
    int best = -1;
    double best_cost = current_cost;
    for (int i = 0; i < N; ++i)
        if (costs[i] < best_cost) { best_cost = costs[i]; best = i; }

    if (best < 0) return false;  // no improvement; odb untouched

    // ---- candidate distribution diagnostic (TDP_DIAG_CAND) ----
    if (diagCandEnabled()) {
        ++g_cand.cells_evaluated;
        if (best >= 0) {
            ++g_cand.cells_found_better;
            double ratio = (current_cost - best_cost) / current_cost;
            g_cand.sum_improve_ratio += ratio;
            g_cand.max_improve_ratio = std::max(g_cand.max_improve_ratio, ratio);
            g_cand.sum_current_cost += current_cost;

            int n_better = 0;
            double minc = costs[0], maxc = costs[0];
            for (double c : costs) {
                if (c < current_cost) ++n_better;
                if (c < minc) minc = c;
                if (c > maxc) maxc = c;
            }
            g_cand.sum_n_better += static_cast<double>(n_better);
            g_cand.sum_pct_better += 100.0 * n_better / static_cast<double>(N);
            double range_ratio = current_cost > 0.0 ? (maxc - minc) / current_cost : 0.0;
            g_cand.sum_cost_range_ratio += range_ratio;
            g_cand.max_cost_range_ratio = std::max(g_cand.max_cost_range_ratio, range_ratio);
        }
    }

    // ---- P0 correlation diagnostic: capture global cost before commit ----
    double corr_before = 0.0;
    double corr_delta_local = current_cost - best_cost; // >0 (par believes it improves)
    if (corr_enabled())
        corr_before = calculateCellCost(cell, alpha, omega, timing_engine_);

    // Commit the best candidate (serial).
    double _tc = prof ? now_ms() : 0.0;
    const Pointdp3D bestpos = candidates[best];
    updateRowIndexMove(cell, cur, bestpos);   // row index incremental maintenance
    cell->location = bestpos;
    timingdp_->setInstanceLocation3D(inst, bestpos.x.v, bestpos.y.v, bestpos.z.v);
    SteinerTreeGenerator stn(timingdp_);
    double _s = prof ? now_ms() : 0.0;
    stn.updateSteinerTrees(cell);
    if (prof) g_prof.commit_steiner += now_ms() - _s;
    _s = prof ? now_ms() : 0.0;
    timingdp_->executeTimingDrivenI(cell, 2);
    if (prof) g_prof.commit_exec += now_ms() - _s;
    _s = prof ? now_ms() : 0.0;
    refreshAffectedPinTiming2(cell);
    if (prof) g_prof.commit_refresh += now_ms() - _s;
    if (prof) { g_prof.commit += now_ms() - _tc; g_prof.commits++; }

    // ---- P0 correlation diagnostic: capture global cost after commit ----
    if (corr_enabled()) {
        double corr_after = calculateCellCost(cell, alpha, omega, timing_engine_);
        double corr_delta_global = corr_before - corr_after;
        double baseline_gap = std::abs(current_cost - corr_before);   // local_cur vs global_cur
        double moved_gap    = std::abs(best_cost    - corr_after);    // local_best vs global_after
        g_corr.record(corr_delta_local, corr_delta_global, baseline_gap, moved_gap, current_cost);
        bool opp = (corr_delta_local > 0 && corr_delta_global < 0)
                || (corr_delta_local < 0 && corr_delta_global > 0);
        if (opp) {
            g_corr.record_opp(baseline_gap, moved_gap);
            if (g_corr.opp_printed < 10) {
                logger_->info(TDP, 2055,
                    "[SPIKE_CORR] opp {} dloc={:.6e} dglob={:.6e} baseline_gap={:.6e} moved_gap={:.6e}",
                    cell->name, corr_delta_local, corr_delta_global, baseline_gap, moved_gap);
                ++g_corr.opp_printed;
            }
        }
    }

    return true;
}


// ============================================================================
// LR算法特定方法
// ============================================================================
//local search
/**
 * @brief 执行离散局部搜索（针对单个单元）
 * @param cell 要优化的单元
 * @param search_windows 搜索窗口数量（默认33）
 * @param alpha 晚时序权重（默认5.0）
 * @param omega 早时序权重（默认1.0）
 * @return 是否找到了更好的位置
 */
bool LSP3D::discreteLocalSearch(Node3D* cell, int search_windows, double alpha, double omega,
                               TimingEngine3D* timing_engine) {
    if (!cell || cell->is_fixed || !arch_ || !network_) {
        return false;
    }
    double time_sta = 0.0;
    double time_fresh = 0.0;
    // 获取当前单元信息
    Pointdp3D current_location = cell->location;
    std::string current_layer = cell->layer;
    odb::dbOrientType current_orient = cell->getOrient();  // 保存当前方向
    //生成候选位置
    std::vector<Pointdp3D> candidates = timingDrivenWindowSearch(cell, search_windows, params_.window_x, params_.window_y);
    if (candidates.empty()) {
        // std::cerr << "Warning: No candidate positions found for cell " << cell->name << std::endl;
        return false;
    }
    // 2. 计算当前位置的成本
    SteinerTreeGenerator steiner_generator(timingdp_);
    // steiner_generator.writeSpef("tdp_initial");
    double current_cost = calculateCellCost(cell, alpha, omega, timing_engine_);
    // 获取行信息
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(current_layer);
    // 获取当前层的所有行
    const auto& rows = arch_->getRows(layer_type);
    if (rows.empty()) {
        std::cerr << "Error: No rows found for layer " << current_layer << std::endl;
        return false;
    }
    Architecture3D::Row3D* row = findCellRow(cell, rows);
    if (!row) {
        std::cerr << "Error: Cannot find row for cell " << cell->name 
                << " at (" << current_location.x.v << ", " << current_location.y.v 
                << ") on layer " << current_layer << std::endl;
        return false;
    }
    
    // 3. 评估每个候选位置
    double best_cost = current_cost;
    Pointdp3D best_location = current_location;
    bool found_better = false;
    odb::dbOrientType candidate_orient = current_orient;
    odb::dbOrientType best_orient = current_orient;
    odb::dbInst* inst = cell->db_inst;

    // ========= 解析模型模式：不修改路由，不需要 global_router =========
    auto* global_router = timingdp_->getGlobalRouter();
    auto& global_routes = global_router->getRoutes();
    std::vector<odb::dbNet*> affected_nets = cell->connected_nets;
    // ========= 解析模型模式：不修改路由，不需要保存/恢复 =========
    std::map<odb::dbNet*, grt::GRoute> saved_routes;
    for (auto* db_net : affected_nets) {
        if (!db_net) continue;
        auto it = global_routes.find(db_net);
        if (it != global_routes.end())
            saved_routes[db_net] = it->second;
    }

    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& candidate = candidates[i];
        
        // 临时移动单元到候选位置
        cell->location = candidate;
        bool need = (std::abs(((candidate.y.v - row->getBottom().v) /  row->getHeight().v )% 2) == 1);
        
        
        if (!need) {
            // 保存当前方向
            odb::dbOrientType saved_orient = cell->getOrient();
            
            // 执行翻转
            flipCell(cell);
            candidate_orient = cell->getOrient();

            // 恢复方向（用于后续恢复）
            cell->setOrient(saved_orient);

        }
        // 计算候选位置的成本
        timingdp_->setInstanceLocation3D(inst, cell->location.x.v, cell->location.y.v, cell->location.z.v);

        // ========= STA 模式（已注释） =========
        steiner_generator.updateSteinerTrees(cell);
        // steiner_generator.writeSpef("LSP_output");
        auto t1 = std::chrono::high_resolution_clock::now();
        timingdp_->executeTimingDrivenI(cell, 2);
        auto t2 = std::chrono::high_resolution_clock::now();
        refreshAffectedPinTiming(cell);
        // compareLocalVsFullSTA(cell);
        auto t3 = std::chrono::high_resolution_clock::now();
        double sta_ms    = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double refresh_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
        time_sta += sta_ms;
        time_fresh += refresh_ms;
        double candidate_cost = calculateCellCost(cell, alpha, omega, timing_engine_);
        
        // 检查是否找到更好的位置
        if (candidate_cost < best_cost) {
            best_cost = candidate_cost;
            best_location = candidate;
            best_orient = candidate_orient;
            found_better = true;
            
            // std::cout << "  -> Found better position! New cost: " << best_cost 
            //           << " (improvement: " << (current_cost - best_cost) << ")" << std::endl;
        }
        // Restore original routes for next candidate
        // ========= STA 模式（已注释） =========
        for (const auto& [db_net, route] : saved_routes)
            global_routes[db_net] = route;
    }
    // logger_ ->info(TDP, 881, "sta_inter time{}", time_sta);
    // logger_ ->info(TDP, 882, "fresh_inter time{}", time_fresh);
    // 4. 恢复单元位置
    cell->location = current_location;
    // 恢复方向
    cell->setOrient(current_orient);

    // 5. 如果找到更好的位置，应用它
    if (found_better) {
        updateRowIndexMove(cell, current_location, best_location);  // row index maintenance
        cell->location = best_location;
        cell->setOrient(best_orient);
        timingdp_->setInstanceLocation3D(inst, best_location.x.v, best_location.y.v, best_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        // sta_->updateTiming(false);
        timingdp_->executeTimingDrivenI(cell, 3);
        refreshAffectedPinTiming2(cell);
        return true;
    }
    else{
        timingdp_->setInstanceLocation3D(inst, current_location.x.v,
                                  current_location.y.v,
                                  current_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        sta_->updateTiming(false);
        refreshAffectedPinTiming2(cell);
    }
    
    return false;
}

// ============================================================================
// discreteLocalSearchNoCommit — same as discreteLocalSearch but evaluates each
// candidate with a compute-only local STA (no write-back to the global graph).
// All candidates are scored against the SAME frozen global baseline; only the
// chosen best candidate is committed. Kept separate from discreteLocalSearch
// for A/B comparison.
// ============================================================================
bool LSP3D::discreteLocalSearchNoCommit(Node3D* cell, int search_windows,
                                        double alpha, double omega,
                                        TimingEngine3D* timing_engine) {
    if (!cell || cell->is_fixed || !arch_ || !network_) {
        return false;
    }
    double time_sta = 0.0;
    double time_fresh = 0.0;
    Pointdp3D current_location = cell->location;
    std::string current_layer = cell->layer;
    odb::dbOrientType current_orient = cell->getOrient();

    std::vector<Pointdp3D> candidates = timingDrivenWindowSearch(cell, search_windows, params_.window_x, params_.window_y);
    if (candidates.empty()) {
        return false;
    }

    SteinerTreeGenerator steiner_generator(timingdp_);
    double current_cost = calculateCellCost(cell, alpha, omega, timing_engine_);

    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(current_layer);
    const auto& rows = arch_->getRows(layer_type);
    if (rows.empty()) {
        std::cerr << "Error: No rows found for layer " << current_layer << std::endl;
        return false;
    }
    Architecture3D::Row3D* row = findCellRow(cell, rows);
    if (!row) {
        std::cerr << "Error: Cannot find row for cell " << cell->name
                << " at (" << current_location.x.v << ", " << current_location.y.v
                << ") on layer " << current_layer << std::endl;
        return false;
    }

    double best_cost = current_cost;
    Pointdp3D best_location = current_location;
    bool found_better = false;
    odb::dbOrientType candidate_orient = current_orient;
    odb::dbOrientType best_orient = current_orient;
    odb::dbInst* inst = cell->db_inst;

    auto* global_router = timingdp_->getGlobalRouter();
    auto& global_routes = global_router->getRoutes();
    std::vector<odb::dbNet*> affected_nets = cell->connected_nets;
    std::map<odb::dbNet*, grt::GRoute> saved_routes;
    for (auto* db_net : affected_nets) {
        if (!db_net) continue;
        auto it = global_routes.find(db_net);
        if (it != global_routes.end())
            saved_routes[db_net] = it->second;
    }

    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& candidate = candidates[i];

        cell->location = candidate;
        bool need = (std::abs(((candidate.y.v - row->getBottom().v) /  row->getHeight().v )% 2) == 1);

        if (!need) {
            odb::dbOrientType saved_orient = cell->getOrient();
            flipCell(cell);
            candidate_orient = cell->getOrient();
            cell->setOrient(saved_orient);
        }
        timingdp_->setInstanceLocation3D(inst, cell->location.x.v, cell->location.y.v, cell->location.z.v);

        steiner_generator.updateSteinerTrees(cell);

        // ===== compute-only local STA: build local PtGraph, DO NOT commit =====
        auto t1 = std::chrono::high_resolution_clock::now();
        std::unique_ptr<lgista::PtGraph> pg(timingdp_->analyzeOnly(cell, 2));
        auto t2 = std::chrono::high_resolution_clock::now();
        refreshAffectedPinTimingFromPtGraph(cell, pg.get());
        auto t3 = std::chrono::high_resolution_clock::now();
        time_sta   += std::chrono::duration<double, std::milli>(t2 - t1).count();
        time_fresh += std::chrono::duration<double, std::milli>(t3 - t2).count();

        double candidate_cost = calculateCellCost(cell, alpha, omega, timing_engine_);

        if (candidate_cost < best_cost) {
            best_cost = candidate_cost;
            best_location = candidate;
            best_orient = candidate_orient;
            found_better = true;
        }
        // Restore original routes for next candidate (clean baseline)
        for (const auto& [db_net, route] : saved_routes)
            global_routes[db_net] = route;
    }

    // 恢复单元位置/方向
    cell->location = current_location;
    cell->setOrient(current_orient);

    // 如果找到更好的位置，应用并提交一次
    if (found_better) {
        updateRowIndexMove(cell, current_location, best_location);  // row index maintenance
        cell->location = best_location;
        cell->setOrient(best_orient);
        timingdp_->setInstanceLocation3D(inst, best_location.x.v, best_location.y.v, best_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        timingdp_->executeTimingDrivenI(cell, 3);
        refreshAffectedPinTiming2(cell);
        return true;
    }
    else {
        timingdp_->setInstanceLocation3D(inst, current_location.x.v,
                                  current_location.y.v,
                                  current_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        sta_->updateTiming(false);
        refreshAffectedPinTiming2(cell);
    }

    return false;
}

/**
 * @brief 执行全局离散局部搜索（遍历所有单元）
 */
void LSP3D::discreteLocalSearchAll(TimingEngine3D* timing_engine) {
    if (!network_) {
        std::cerr << "Error: Network not initialized" << std::endl;
        return;
    }

    // Row occupancy index: build once (first sweep) and thereafter maintain it
    // incrementally via updateRowIndexMove at every net-committed move. Verified
    // (rebuild-vs-incremental mismatch == 0) that incremental tracking is exact,
    // so no per-sweep rebuild is needed. Between sweeps only createSteninerRoute
    // (no cell moves) and updatePinWeights (no cell moves) run, so the index
    // stays valid across LR iterations.
    if (row_index_.empty())
        buildRowIndexInto(row_index_);

    int moved_cells = 0;
    int total_cells = 0;
    
    auto priority_order = buildCriticalPathPriorityOrder(timing_engine, 2000);

    // P1: gate-vs-wire delay decomposition of the top violating paths (env-gated).
    if (std::getenv("TDP_DIAG_P1"))
        diagPathDecomp(timing_engine, logger_);

    // T1: log true Late TNS every K committed moves (env TDP_DIAG_T1=K, default 200).
    // NOTE: getTimingStats() forces an STA update => enabling T1 perturbs the run.
    int t1_k = 0;
    if (const char* v = std::getenv("TDP_DIAG_T1")) {
        t1_k = std::atoi(v);
        if (t1_k <= 0) t1_k = 200;
    }

    int total_movable = std::count_if(priority_order.begin(), priority_order.end(),
        [](Node3D* c) { return c && !c->is_fixed; });
    
    const int bar_width = 40;
    int report_step = std::max(1, total_movable / 1000);

    double totaltime = 0.0;

    // ---- TDP_SKIP_HIGH_FANOUT=<N>: opt-in skip of cells on ultra-high-fanout nets
    static const int FANOUT_SKIP_THR =
        []() { auto* v = std::getenv("TDP_SKIP_HIGH_FANOUT"); return v ? std::atoi(v) : 0; }();
    
    // 遍历所有单元（关键路径优先）
    for (auto* cell : priority_order) {
        if (!cell || cell->is_fixed || isSequentialCell(cell)) {
            // 固定单元保持原位
            continue;
        }

        if (FANOUT_SKIP_THR > 0) {
            bool skip = false;
            for (auto* p : cell->getPins()) {
                if (!p) continue;
                Net3D* n = p->getNet();
                if (n && (int)n->getPins().size() > FANOUT_SKIP_THR) { skip = true; break; }
            }
            if (skip) continue;
        }
        
        total_cells++;
        auto t4 = std::chrono::high_resolution_clock::now();
        // 执行局部搜索。TDP_PAR=并行；TDP_NOCOMMIT=串行 no-commit；默认=原版。
        static const bool use_par      = (std::getenv("TDP_PAR") != nullptr);
        static const bool use_nocommit = (std::getenv("TDP_NOCOMMIT") != nullptr);
        bool moved;
        if (use_par)
            moved = discreteLocalSearchParallel(cell, 40, 5.0, 1.0, timing_engine_);
        else if (use_nocommit)
            moved = discreteLocalSearchNoCommit(cell, 20, 5.0, 1.0, timing_engine_);
        else
            moved = discreteLocalSearch(cell, 20, 5.0, 1.0, timing_engine_);
        if (moved) {
            moved_cells++;
        }
        if (t1_k > 0 && moved && (moved_cells % t1_k == 0)) {
            auto st = timing_engine_->getTimingStats();
            logger_->info(TDP, 2073,
                "[T1] after {} moves: Late WNS={:.6f}ns TNS={:.6f}ns",
                moved_cells, st.late_wns * 1e9, st.late_tns * 1e9);
        }
        auto t5 = std::chrono::high_resolution_clock::now();

        double dis_ms = std::chrono::duration<double, std::milli>(t5 - t4).count();
        totaltime += dis_ms;

        if (total_cells % report_step == 0 || total_cells == total_movable) {
            double pct = 100.0 * total_cells / total_movable;
            int pos = bar_width * total_cells / total_movable;
            std::cout << "\r进度: [" << std::string(pos, '#')
                      << std::string(bar_width - pos, '.') << "] "
                      << std::fixed << std::setprecision(1)
                      << std::setw(5) << pct << "% ("
                      << total_cells << "/" << total_movable << ")"
                      << std::flush;
        }
    }

    if (corr_enabled())
        dumpCorrStats(logger_);
    if (prof_enabled())
        dumpProfStats(logger_);
    if (diagCandEnabled())
        dumpCandStats(logger_);

    std::cout << std::endl;
    std::cout << "Discrete local search completed. Moved " 
              << moved_cells << " out of " << total_cells << " movable cells" << std::endl;
    logger_ ->info (TDP, 1991, "total time: {}", totaltime);
}

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * 查找网络的驱动引脚
 */
Pin3D* LSP3D::findDriverPin(Net3D* net, Pin3D* exclude_pin) const {
    if (!net) return nullptr;
    
    auto net_pins = net->getPins();
    
    for (auto* pin : net_pins) {
        if (pin == exclude_pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_OUT) {
            return pin;  // 找到输出引脚（驱动）
        }
    }
    
    return nullptr;
}

// ============================================================================
// 时序缓存刷新
// ============================================================================
void LSP3D::refreshAffectedPinTiming(Node3D* cell) {
    refreshPinTiming2(cell);
    for (auto* pin : cell->getPins()) {
        if (pin->getDirection() != Pin3D::Dir_OUT) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;
        for (auto* net_pin : net->getPins()) {
            Node3D* pin_cell = net_pin->getNode();
            if (pin_cell && pin_cell != cell) {
                refreshPinTiming2(pin_cell);
            }
        }
    }
}

void LSP3D::refreshAffectedPinTiming2(Node3D* cell) {
    refreshPinTiming2(cell);
    for (auto* pin : cell->getPins()) {
        if (pin->getDirection() != Pin3D::Dir_OUT) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;
        for (auto* net_pin : net->getPins()) {
            Node3D* pin_cell = net_pin->getNode();
            if (pin_cell && pin_cell != cell) {
                refreshPinTiming2(pin_cell);
            }
        }
    }
}

// Extract worst early(min)/late(max) arrival for a PtVertex from its local
// (non-committed) PtGraph paths. Mirrors getPinTimingInfoDirect's min/max
// worst-arrival selection, but reads the local PtGraph instead of the global
// sta::Graph. Returns false if no usable local path exists.
static bool ptVertexLocalArrival(lgista::PtGraph* pg, lgista::PtVertex& ptv,
                                 sta::StaState* sta_state,
                                 double& early_out, double& late_out) {
    sta::Path* paths = ptv.paths();
    if (!paths) return false;
    sta::TagGroup* tg = pg->tagGroup(ptv);
    if (!tg) return false;
    const size_t count = tg->pathCount();

    bool have_min = false, have_max = false;
    float min_arr = 0.0f, max_arr = 0.0f;
    for (size_t i = 0; i < count; ++i) {
        sta::Tag* tag = paths[i].tag(sta_state);
        if (!tag || tag->isGenClkSrcPath()) continue;
        const float arr = sta::delayAsFloat(paths[i].arrival());
        if (tag->minMax() == sta::MinMax::min()) {
            if (!have_min || arr < min_arr) { min_arr = arr; have_min = true; }
        } else {
            if (!have_max || arr > max_arr) { max_arr = arr; have_max = true; }
        }
    }
    if (!have_min && !have_max) return false;
    early_out = have_min ? min_arr : max_arr;
    late_out  = have_max ? max_arr : min_arr;
    return true;
}

void LSP3D::refreshAffectedPinTimingFromPtGraph(Node3D* cell, lgista::PtGraph* pt_graph) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;
    sta::Graph* graph = sta_ ? sta_->graph() : nullptr;

    auto refresh_one = [&](Node3D* c) {
        if (!c) return;
        for (auto* pin : c->getPins()) {
            if (!pin) continue;
            sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
            if (!sta_pin) continue;

            bool filled = false;
            if (pt_graph && graph) {
                sta::Vertex *v = nullptr, *bidir = nullptr;
                graph->pinVertices(sta_pin, v, bidir);
                if (v) {
                    lgista::PtVertex* ptv = pt_graph->ptVertex(v);
                    if (ptv) {
                        double early = 0.0, late = 0.0;
                        if (ptVertexLocalArrival(pt_graph, *ptv, sta_, early, late)) {
                            pin->setEarlyArrival(early);
                            pin->setLateArrival(late);
                            filled = true;
                        }
                    }
                }
            }
            if (!filled) {
                // Fallback: pin outside the local PtGraph depth → read frozen
                // global baseline value.
                auto ti = timing_engine_->getPinTimingInfoDirect(sta_pin);
                if (ti.is_valid) {
                    pin->setEarlyArrival(ti.early_arrival);
                    pin->setLateArrival(ti.late_arrival);
                }
            }
        }
    };

    refresh_one(cell);
    for (auto* pin : cell->getPins()) {
        if (pin->getDirection() != Pin3D::Dir_OUT) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;
        for (auto* net_pin : net->getPins()) {
            Node3D* pin_cell = net_pin->getNode();
            if (pin_cell && pin_cell != cell) {
                refresh_one(pin_cell);
            }
        }
    }
}

// ============================================================================
// 网络相关辅助方法
// ============================================================================
void LSP3D::initWeight(Timingdp* timingdp, TimingEngine3D* timing_engine)
{

    auto* network_3d = timingdp->getNetwork3D();
        if (!network_3d) {
            std::cout << "Error: Network3D not available" << std::endl;
        }else{

            auto& nets = timingdp->getNets();

            for (auto* db_net : nets) {
                auto* net = network_3d->getNet(db_net);
                if(!net){
                std::cout << "  db_net pointer: " << db_net 
                  << ", name: " << db_net->getName() 
                  << ", signal type: " << (db_net->getSigType().getString()) << std::endl;
                    continue;
                }
                auto pins = net->getPins();
                
                for (auto* pin : pins) {
                    
                    // 设置默认权重为1.0
                    pin->setEarlyWeight(1.0);
                    pin->setLateWeight(1.0);
                    // pin->setEarlyWeightini(1.0);
                    // pin->setLateWeightini(1.0);
                    
                    // 获取引脚名称
                    std::string pin_name = pin->getName();
                    
                    // 尝试从timing_engine获取时序信息
                    try {
                        // 获取引脚时序信息
                        auto timing_info = timing_engine_->getPinTimingInfo(pin_name);
                        
                        if (timing_info.is_valid) {
                            // 设置真实的时序信息
                            pin->setEarlySlack(timing_info.early_slack);
                            pin->setLateSlack(timing_info.late_slack);
                            pin->setEarlyArrival(timing_info.early_arrival);
                            pin->setLateArrival(timing_info.late_arrival);
                            pin->setEarlyRequired(timing_info.early_required);
                            pin->setLateRequired(timing_info.late_required);
                        } else {
                            // 如果无法获取时序信息，使用默认值
                            pin->setEarlySlack(0.0);
                            pin->setLateSlack(0.0);
                            pin->setEarlyArrival(0.0);
                            pin->setLateArrival(0.0);
                            pin->setEarlyRequired(0.0);
                            pin->setLateRequired(0.0);
                            std::cout << "Error timing info "<< std::endl;
                        }
                    } catch (const std::exception& e) {
                        // 如果获取失败，使用默认值
                        // std::cerr << "Warning: Failed to get timing info for pin " 
                                // << pin_name << ": " << e.what() << std::endl;
                        pin->setEarlySlack(0.0);
                        pin->setLateSlack(0.0);
                        pin->setEarlyArrival(0.0);
                        pin->setLateArrival(0.0);
                        pin->setEarlyRequired(0.0);
                        pin->setLateRequired(0.0);
                    }
                }
            }
        }

}

} // namespace tdp
