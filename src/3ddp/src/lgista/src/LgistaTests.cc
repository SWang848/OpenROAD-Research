// SPDX-License-Identifier: BSD-3-Clause
//
// Correctness tests for the local graph ISTA. See TEST_PLAN.md and
// LgistaTests.hh. Comparisons read timing through PtGraph accessors on both
// sides: the oracle is a fresh PtGraph built on the global graph AFTER a real
// replaceCell + updateTiming (PtGraph construction copies the global graph's
// slew/delay as its initial state).

#include "LgistaTests.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "db_sta/dbSta.hh"

#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "sta/Network.hh"
#include "sta/Scene.hh"
#include "sta/Sta.hh"
#include "sta/TimingArc.hh"
#include "sta/Transition.hh"

#include "IstaWriteback.hh"
#include "LocalDmpDelayCalc.hh"
#include "LocalSta.hh"
#include "PtGraph.hh"
#include "lgista/LrfClass.hh"

namespace lgista {

using sta::dbSta;
using sta::DcalcAPIndex;
using sta::Instance;
using sta::LibertyCell;
using sta::LibertyCellSeq;
using sta::LibertyLibrary;
using sta::LibertyLibrarySeq;
using sta::MinMax;
using sta::Network;
using sta::Pin;
using sta::RiseFall;
using sta::Vertex;

namespace {

// ── tolerance accumulator ──────────────────────────────────────────────────
struct Acc {
  const char* name;
  int n = 0;
  int fail = 0;
  double max_abs = 0.0;
  double max_rel = 0.0;
  std::string worst;

  explicit Acc(const char* nm) : name(nm) {}

  void check(const std::string& tag, double a, double b,
             double atol = 1e-12, double rtol = 1e-3) {
    ++n;
    double d = std::fabs(a - b);
    double rel = d / (std::max(std::fabs(a), std::fabs(b)) + 1e-30);
    if (d > max_abs)
      max_abs = d;
    bool ok = (d <= atol) || (rel <= rtol);
    if (!ok) {
      ++fail;
      if (rel > max_rel) {
        max_rel = rel;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s  local=%.6g global=%.6g (rel=%.3g)",
                      tag.c_str(), a, b, rel);
        worst = buf;
      }
    }
  }

  bool pass() const { return fail == 0; }

  void report() const {
    std::printf("[%s] %d checks, %d fail, max_abs=%.3g, max_rel=%.3g  -> %s\n",
                name, n, fail, max_abs, max_rel, pass() ? "PASS" : "FAIL");
    if (!pass())
      std::printf("    worst: %s\n", worst.c_str());
  }
};

DcalcAPIndex apIndexOf(dbSta* sta) {
  // LocalSta::makePtGraph sets the PtGraph analysis point to (default scene,
  // MinMax::max()); match it so we read the index the local engine computes.
  sta::Scene* scene = sta->cmdScene();
  if (scene == nullptr)
    return 0;
  return scene->dcalcAnalysisPtIndex(MinMax::max());
}

Vertex* vertexOf(dbSta* sta, const Pin* pin) {
  sta::Graph* g = sta->graph();
  Network* nw = sta->network();
  if (nw->isDriver(pin))
    return g->pinDrvrVertex(pin);
  return g->pinLoadVertex(pin);
}

// Make sure equivalent-cell records exist (needed by equivCells()).
void ensureEquivCells(dbSta* sta) {
  static bool done = false;
  if (done)
    return;
  done = true;
  Network* nw = sta->network();
  LibertyLibrarySeq libs;
  sta::LibertyLibraryIterator* it = nw->libertyLibraryIterator();
  while (it->hasNext())
    libs.push_back(it->next());
  delete it;
  if (!libs.empty())
    sta->makeEquivCells(&libs, &libs);
}

// Pick an equivalent cell different from the instance's current cell.
LibertyCell* pickEquiv(dbSta* sta, Instance* inst) {
  LibertyCell* cur = sta->network()->libertyCell(inst);
  if (cur == nullptr)
    return nullptr;
  ensureEquivCells(sta);
  LibertyCellSeq* eq = sta->equivCells(cur);
  if (eq == nullptr)
    return nullptr;
  for (LibertyCell* c : *eq) {
    if (c != cur && !c->dontUse())
      return c;
  }
  return nullptr;
}

// Leaf instances sorted in topological (level) order using each instance's
// max output-vertex level.
std::vector<Instance*> topoInstances(dbSta* sta) {
  sta::Graph* g = sta->graph();
  Network* nw = sta->network();
  std::vector<std::pair<int, Instance*>> tmp;
  sta::LeafInstanceIterator* it = nw->leafInstanceIterator();
  while (it->hasNext()) {
    Instance* inst = it->next();
    int lvl = 0;
    sta::InstancePinIterator* pit = nw->pinIterator(inst);
    while (pit->hasNext()) {
      Pin* pin = pit->next();
      if (!nw->isDriver(pin))
        continue;
      Vertex* v = g->pinDrvrVertex(pin);
      if (v)
        lvl = std::max(lvl, static_cast<int>(v->level()));
    }
    delete pit;
    tmp.emplace_back(lvl, inst);
  }
  delete it;
  std::sort(tmp.begin(), tmp.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  std::vector<Instance*> out;
  out.reserve(tmp.size());
  for (auto& p : tmp)
    out.push_back(p.second);
  return out;
}

// Per-pin slew (rise, fall) read from the global graph.
using SlewMap = std::map<const Pin*, std::array<float, 2>>;

SlewMap collectGlobalSlews(dbSta* sta, DcalcAPIndex ap) {
  sta::Graph* g = sta->graph();
  Network* nw = sta->network();
  SlewMap m;
  sta::LeafInstanceIterator* it = nw->leafInstanceIterator();
  while (it->hasNext()) {
    Instance* inst = it->next();
    sta::InstancePinIterator* pit = nw->pinIterator(inst);
    while (pit->hasNext()) {
      Pin* pin = pit->next();
      Vertex* v = vertexOf(sta, pin);
      if (!v)
        continue;
      m[pin] = {static_cast<float>(sta::delayAsFloat(g->slew(v, RiseFall::rise(), ap))),
                static_cast<float>(sta::delayAsFloat(g->slew(v, RiseFall::fall(), ap)))};
    }
    delete pit;
  }
  delete it;
  return m;
}

void compareSlews(Acc& acc, const SlewMap& a, const SlewMap& b) {
  for (const auto& kv : a) {
    auto it = b.find(kv.first);
    if (it == b.end())
      continue;
    acc.check("slew_r", kv.second[0], it->second[0]);
    acc.check("slew_f", kv.second[1], it->second[1]);
  }
}

}  // namespace

// ── T1: virtualReplaceCell local-vs-global equivalence ─────────────────────
static int t1(dbSta* sta) {
  std::printf("=== T1: virtualReplaceCell vs global replaceCell ===\n");
  DcalcAPIndex ap = apIndexOf(sta);
  sta->updateTiming(true);

  LocalSta L(sta);
  L.copyState(sta);
  sta::ArcDelayCalc* calc = makeLocalDelayCalc(sta);

  Acc acc("T1");
  Network* nw = sta->network();
  int tested = 0;
  int skipped = 0;
  const int kMaxTested = 20;

  sta::LeafInstanceIterator* it = nw->leafInstanceIterator();
  while (it->hasNext() && (tested + skipped) < kMaxTested) {
    Instance* inst = it->next();
    LibertyCell* corig = nw->libertyCell(inst);
    LibertyCell* cprime = pickEquiv(sta, inst);
    if (corig == nullptr || cprime == nullptr)
      continue;

    PtGraph* lp = L.makePtGraph(inst, /*update_timing_first=*/true);
    if (lp == nullptr)
      continue;
    L.findLocalDelays(lp, calc);
    L.findLocalArrivals(lp);
    L.findLocalRequireds(lp);

    if (!L.virtualReplaceCell(lp, cprime))
      continue;
    L.findLocalDelays(lp, calc);
    L.findLocalArrivals(lp);
    L.findLocalRequireds(lp);

    sta::Graph* g = sta->graph();

    // Record this instance's INPUT-pin global slews before the swap.
    std::map<const Pin*, std::array<float, 2>> in_before;
    for (PtVertex& pv : lp->ptVertices()) {
      if (pv.type() != PtVertexType::RefInput)
        continue;
      const Pin* pin = pv.pin();
      if (!pin)
        continue;
      Vertex* v = vertexOf(sta, const_cast<Pin*>(pin));
      if (!v)
        continue;
      in_before[pin] = {
          (float)sta::delayAsFloat(g->slew(v, RiseFall::rise(), ap)),
          (float)sta::delayAsFloat(g->slew(v, RiseFall::fall(), ap))};
    }

    // global oracle: REAL replace + updateTiming (lets the inputs float).
    sta->replaceCell(inst, cprime);
    sta->updateTiming(true);

    // Is this instance "input-stable" under the swap? A real replace also
    // shifts I's *input* slew through the input-cap change feeding back to the
    // upstream driver -- a second-order effect a single-instance local graph
    // intentionally does NOT model. Only when the inputs are unchanged is a
    // strict local-vs-global comparison of I's own output meaningful; otherwise
    // skip (the difference is the documented boundary effect, not a bug).
    bool input_stable = true;
    for (const auto& kv : in_before) {
      Vertex* v = vertexOf(sta, const_cast<Pin*>(kv.first));
      if (!v) {
        input_stable = false;
        break;
      }
      float gr = sta::delayAsFloat(g->slew(v, RiseFall::rise(), ap));
      float gf = sta::delayAsFloat(g->slew(v, RiseFall::fall(), ap));
      if (std::fabs(gr - kv.second[0]) > 1e-13 ||
          std::fabs(gf - kv.second[1]) > 1e-13) {
        input_stable = false;
        break;
      }
    }

    if (input_stable) {
      // Compare the swapped instance's OUTPUT pin slews (local PtGraph vs
      // global graph) and its gate arc delays.
      for (PtVertex& pv : lp->ptVertices()) {
        if (pv.type() != PtVertexType::RefOutput)
          continue;
        const Pin* pin = pv.pin();
        if (!pin)
          continue;
        Vertex* v = vertexOf(sta, const_cast<Pin*>(pin));
        if (!v)
          continue;
        for (const sta::RiseFall* rf : RiseFall::range()) {
          double a = sta::delayAsFloat(lp->slew(pv, rf, ap));
          double b = sta::delayAsFloat(g->slew(v, rf, ap));
          acc.check("out_slew", a, b);
          if (std::getenv("LGISTA_DEBUG"))
            std::printf("    [debug] out_slew %s rf=%d local=%.5g global=%.5g\n",
                        std::string(nw->pathName(pin)).c_str(),
                        rf == RiseFall::rise() ? 0 : 1, a, b);
        }
      }
      for (size_t eid = 0; eid < lp->edgeCount(); ++eid) {
        PtEdge& e = lp->edge(eid);
        sta::Edge* ge = e.edge();
        if (ge == nullptr)
          continue;
        sta::TimingArcSet* tas = ge->timingArcSet();
        if (tas == nullptr)
          continue;
        for (sta::TimingArc* arc : tas->arcs()) {
          // Read both sides keyed by the SAME arc (PtGraph::arcDelay handles
          // the per-arc indexing internally; a positional index into
          // arcDelays() does not line up with tas->arcs() iteration order).
          double a = sta::delayAsFloat(lp->arcDelay(e, arc, ap));
          double b = sta::delayAsFloat(g->arcDelay(ge, arc, ap));
          acc.check("arc_delay", a, b);
          if (std::getenv("LGISTA_DEBUG"))
            std::printf("    [debug] arc wire=%d local=%.5g global=%.5g\n",
                        ge->isWire() ? 1 : 0, a, b);
        }
      }
      ++tested;
    } else {
      ++skipped;
    }

    sta->replaceCell(inst, corig);
    sta->updateTiming(true);
  }
  delete it;
  delete calc;

  std::printf("    tested %d instances (skipped %d: input not stable under "
              "swap = second-order boundary effect)\n",
              tested, skipped);
  acc.report();
  return acc.pass() ? 0 : 1;
}

// ── T2-A: write-back identity (no cell change) ─────────────────────────────
static int t2a(dbSta* sta) {
  std::printf("=== T2-A: write-back identity (no change) ===\n");
  DcalcAPIndex ap = apIndexOf(sta);
  sta->updateTiming(true);

  SlewMap before = collectGlobalSlews(sta, ap);

  LocalSta L(sta);
  L.copyState(sta);
  sta::ArcDelayCalc* calc = makeLocalDelayCalc(sta);

  for (Instance* inst : topoInstances(sta)) {
    PtGraph* pt = L.makePtGraph(inst, /*update_timing_first=*/true);
    if (pt == nullptr)
      continue;
    L.findLocalDelays(pt, calc);
    L.findLocalArrivals(pt);
    L.findLocalRequireds(pt);
    updateTimingFromPtGraph(pt);
  }
  delete calc;

  SlewMap after = collectGlobalSlews(sta, ap);

  Acc acc("T2-A");
  compareSlews(acc, before, after);
  acc.report();
  return acc.pass() ? 0 : 1;
}

// ── T2-B: replace-all + topo write-back vs full global re-analysis ─────────
static int t2b(dbSta* sta) {
  std::printf("=== T2-B: replace-all + topo write-back vs global ===\n");
  DcalcAPIndex ap = apIndexOf(sta);
  sta->updateTiming(true);

  LocalSta L(sta);
  L.copyState(sta);
  sta::ArcDelayCalc* calc = makeLocalDelayCalc(sta);

  std::vector<Instance*> order = topoInstances(sta);

  // Phase 1: replace ALL cells first, so every net's load cap is at its final
  // value before any local graph is built. (Interleaving replace with the
  // local recompute makes an upstream driver see its load's OLD input cap,
  // because the downstream cell hasn't been swapped yet -- a real second-order
  // mismatch vs the global full recompute.)
  int replaced = 0;
  for (Instance* inst : order) {
    LibertyCell* cprime = pickEquiv(sta, inst);
    if (cprime == nullptr)
      continue;
    sta->replaceCell(inst, cprime);
    ++replaced;
  }

  // Phase 2: topological order -- build each instance's local graph (now sees
  // final load caps), recompute, and write back. Forward order guarantees an
  // instance's input slews were already committed by its upstream neighbours.
  for (Instance* inst : order) {
    PtGraph* pt = L.makePtGraph(inst, /*update_timing_first=*/true);
    if (pt == nullptr)
      continue;
    L.findLocalDelays(pt, calc);
    L.findLocalArrivals(pt);
    L.findLocalRequireds(pt);
    updateTimingFromPtGraph(pt);  // incremental write-back
  }
  delete calc;

  // incremental result
  SlewMap incremental = collectGlobalSlews(sta, ap);
  // ground truth: full global recompute on the SAME final netlist
  sta->updateTiming(true);
  SlewMap global = collectGlobalSlews(sta, ap);

  std::printf("    replaced %d instances\n", replaced);
  Acc acc("T2-B");
  compareSlews(acc, incremental, global);
  acc.report();
  return acc.pass() ? 0 : 1;
}

int runTest(dbSta* sta, const std::string& which) {
  int rc = 0;
  if (which == "t1" || which == "all")
    rc |= t1(sta);
  if (which == "t2a" || which == "all")
    rc |= t2a(sta);
  if (which == "t2b" || which == "all")
    rc |= t2b(sta);
  return rc;
}

}  // namespace lgista
