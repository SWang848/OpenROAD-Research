// SPDX-License-Identifier: BSD-3-Clause
//
// LocalGraphIsta facade implementation. See LocalGraphIsta.hh.

#include "lgista/LocalGraphIsta.hh"

#include <chrono>
#include <cstdio>
#include <memory>

#include "db_sta/dbSta.hh"
#include "sta/Network.hh"
#include "sta/Sta.hh"

#include "LocalSta.hh"
#include "PtGraph.hh"
#include "LocalDmpDelayCalc.hh"
#include "IstaWriteback.hh"

namespace lgista {

static inline double now_ms() {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

LocalGraphIsta::LocalGraphIsta(sta::dbSta* sta) : sta_(sta)
{
  // Construction is deferred to prepare(): LocalSta::copyState and the local
  // delay calculator need a built/levelized global timing graph, which only
  // exists after updateTiming(). Building them here (before any timing run)
  // dereferences a not-yet-built graph and faults.
}

LocalGraphIsta::~LocalGraphIsta()
{
  delete arc_delay_calc_;
  delete local_sta_;
}

void LocalGraphIsta::prepare()
{
  if (prepared_)
    return;
  // Bring the global timing graph fully up to date (ensures graph build +
  // levelization + arrivals + requireds) so the local graph can read upstream
  // slews / arrivals from it. updateTiming(true) is robust even when the
  // design has no clocks/constraints (findRequireds alone can fault then).
  sta_->updateTiming(/*full=*/true);
  local_sta_ = new LocalSta(sta_);
  local_sta_->copyState(sta_);
  arc_delay_calc_ = lgista::makeLocalDelayCalc(sta_);
  prepared_ = true;
}

PtGraph* LocalGraphIsta::analyzeInstance(sta::Instance* inst, bool print_report)
{
  if (inst == nullptr)
    return nullptr;
  if (!prepared_)
    prepare();

  PtGraph* pt_graph = local_sta_->makePtGraph(inst, /*update_timing_first=*/true);
  if (pt_graph == nullptr)
    return nullptr;

  local_sta_->findLocalDelays(pt_graph, arc_delay_calc_);
  local_sta_->findLocalArrivals(pt_graph);
  local_sta_->findLocalRequireds(pt_graph);

  if (print_report)
    local_sta_->printLocalTiming(pt_graph);

  return pt_graph;
}

PtGraph* LocalGraphIsta::analyzeInstance(const std::string& inst_name,
                                         bool print_report)
{
  sta::Network* network = sta_->network();
  sta::Instance* inst = network->findInstance(inst_name.c_str());
  if (inst == nullptr) {
    std::printf("LocalGraphIsta: instance '%s' not found\n", inst_name.c_str());
    return nullptr;
  }
  return analyzeInstance(inst, print_report);
}

void LocalGraphIsta::commitToGlobal(PtGraph* pt_graph)
{
  if (pt_graph != nullptr)
    updateTimingFromPtGraph(pt_graph);
}

void LocalGraphIsta::ensureInit()
{
  if (prepared_)
    return;
  local_sta_ = new LocalSta(sta_);
  local_sta_->copyState(sta_);
  arc_delay_calc_ = lgista::makeLocalDelayCalc(sta_);
  prepared_ = true;
}

bool LocalGraphIsta::analyzeAndCommit(sta::Instance* inst, int depth)
{
  if (inst == nullptr)
    return false;
  ensureInit();

  PtGraph pt_graph(sta_);
  if (!local_sta_->makePtGraph(&pt_graph, inst, -1, depth))
    return false;

  local_sta_->findLocalDelays(&pt_graph, arc_delay_calc_);
  local_sta_->findLocalArrivals(&pt_graph);
  local_sta_->findLocalRequireds(&pt_graph);

  updateTimingFromPtGraph(&pt_graph);
  return true;
}

std::unique_ptr<PtGraph> LocalGraphIsta::analyzeOnly(sta::Instance* inst, int depth)
{
  return analyzeOnly(inst, depth, nullptr);
}

std::unique_ptr<PtGraph> LocalGraphIsta::analyzeOnly(sta::Instance* inst, int depth,
                                                     AnalyzeTimings* tm)
{
  if (inst == nullptr)
    return nullptr;
  ensureInit();

  double t0 = tm ? now_ms() : 0.0;
  auto pt_graph = std::make_unique<PtGraph>(sta_);
  if (!local_sta_->makePtGraph(pt_graph.get(), inst, -1, depth))
    return nullptr;
  if (tm) { tm->make_ptgraph_ms += now_ms() - t0; t0 = now_ms(); }

  // Reduce parasitics into PtGraph-local PtPiElmore so the delay calc uses
  // them at priority #1 (over the global scene). For nets present in the
  // per-engine local container (trial candidate parasitics) this yields the
  // candidate RC; nets absent there fall back to the scene network. Without
  // this, localParasiticLoad would read the scene PiElmore (current geometry)
  // and ignore any trial parasitics built into the local container.
  local_sta_->localParasitics()->recomputePtParasitics(pt_graph.get());
  if (tm) { tm->recompute_para_ms += now_ms() - t0; t0 = now_ms(); }

  local_sta_->findLocalDelays(pt_graph.get(), arc_delay_calc_);
  if (tm) { tm->delays_ms += now_ms() - t0; t0 = now_ms(); }
  local_sta_->findLocalArrivals(pt_graph.get());
  if (tm) { tm->arrivals_ms += now_ms() - t0; t0 = now_ms(); }
  local_sta_->findLocalRequireds(pt_graph.get());
  if (tm) { tm->requireds_ms += now_ms() - t0; }

  // NOTE: deliberately NO updateTimingFromPtGraph() — global graph stays
  // frozen; caller reads local values from the returned PtGraph.
  return pt_graph;
}

sta::Parasitics* LocalGraphIsta::localParasiticsContainer()
{
  return local_sta_ ? local_sta_->localParasitics() : nullptr;
}

size_t LocalGraphIsta::analyzeAll(bool print_report)
{
  if (!prepared_)
    prepare();

  sta::Network* network = sta_->network();
  size_t built = 0;
  sta::LeafInstanceIterator* iter = network->leafInstanceIterator();
  while (iter->hasNext()) {
    sta::Instance* inst = iter->next();
    if (analyzeInstance(inst, print_report) != nullptr)
      ++built;
  }
  delete iter;
  return built;
}

}  // namespace lgista
