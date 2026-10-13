// SPDX-License-Identifier: BSD-3-Clause
//
// SearchMEEPred implementation. See LrfSearchPred.hh.
// Extracted verbatim from the original lrf TaskArranger.cc.

#include "LrfSearchPred.hh"

#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/TimingRole.hh"

namespace lgista {

SearchMEEPred::SearchMEEPred(sta::StaState* sta) : sta::SearchPred1(sta) {}

bool SearchMEEPred::searchThru(sta::Edge* edge, const sta::Mode* mode) const
{
  const sta::TimingRole* role = edge->role();
  sta::Network* network_ = sta_->network();
  sta::Graph* graph_ = sta_->graph();
  sta::Instance* to_inst
      = network_->instance(graph_->vertex(edge->to())->pin());
  if (to_inst == sta_->network()->topInstance()) {
    return false;
  }
  sta::LibertyCell* to_cell = network_->libertyCell(to_inst);
  if (!to_cell)
    return false;
  return (SearchPred1::searchThru(edge, mode) && role->isWire()
          && !to_cell->hasSequentials());
}

bool SearchMEEPred::searchFrom(const sta::Vertex* from_vertex,
                               const sta::Mode* mode) const
{
  sta::Instance* from_inst = sta_->network()->instance(from_vertex->pin());
  if (from_inst == sta_->network()->topInstance()) {
    return false;
  }
  sta::LibertyCell* from_cell = sta_->network()->libertyCell(from_inst);
  if (!from_cell)
    return false;
  return (!from_cell->hasSequentials()
          && SearchPred0::searchFrom(from_vertex, mode));
}

bool SearchMEEPred::searchTo(const sta::Vertex* to_vertex,
                             const sta::Mode* mode) const
{
  sta::Instance* to_inst = sta_->network()->instance(to_vertex->pin());
  if (to_inst == sta_->network()->topInstance()) {
    return false;
  }
  sta::LibertyCell* to_cell = sta_->network()->libertyCell(to_inst);
  if (!to_cell)
    return false;
  return (!to_cell->hasSequentials()
          && SearchPred0::searchTo(to_vertex, mode));
}

}  // namespace lgista
