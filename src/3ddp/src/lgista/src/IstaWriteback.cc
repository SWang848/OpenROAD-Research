// SPDX-License-Identifier: BSD-3-Clause
//
// updateTimingFromPtGraph implementation. See IstaWriteback.hh.
// Extracted verbatim from the original lrf NetlistTransformation.cc.

#include "IstaWriteback.hh"

#include "PtGraph.hh"
#include "lgista/LrfClass.hh"
#include "sta/Graph.hh"

namespace lgista {

void updateTimingFromPtGraph(PtGraph* pt_graph)
{
  for (sta::VertexId vertex_id : pt_graph->sortedVertexIds()) {
    PtVertex& pt_vertex = pt_graph->ptVertex(vertex_id);
    if (!pt_vertex.vertex())
      continue;
    sta::Vertex* sta_vertex = pt_vertex.vertex();
    PtVertexType type = pt_vertex.type();
    if (type == PtVertexType::RefInput || type == PtVertexType::RefOutput
        || type == PtVertexType::SiblingLoad) {
      pt_graph->writeSlewToGraph(pt_vertex, sta_vertex);
      pt_graph->writePathsToGraph(pt_vertex, sta_vertex);
    }
  }
}

}  // namespace lgista
