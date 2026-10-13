// SPDX-License-Identifier: BSD-3-Clause
//
// Timing search predicate used by the local graph timer. Extracted from
// TaskArranger.hh (part of the sizer-oriented task framework, not included in
// this standalone module). SearchMEEPred restricts the search to combinational
// wire arcs, stopping at sequential cells and the top instance -- the boundary
// used when collecting a local timing graph.

#pragma once

#include "sta/SearchPred.hh"

namespace sta {
class StaState;
class Edge;
class Vertex;
class Mode;
}  // namespace sta

namespace lgista {

class SearchMEEPred : public sta::SearchPred1
{
 public:
  explicit SearchMEEPred(sta::StaState* sta);
  bool searchThru(sta::Edge* edge, const sta::Mode* mode) const override;
  bool searchFrom(const sta::Vertex* from_vertex,
                  const sta::Mode* mode) const override;
  bool searchTo(const sta::Vertex* to_vertex,
                const sta::Mode* mode) const override;
};

}  // namespace lgista
