// SPDX-License-Identifier: BSD-3-Clause
//
// Write API for the local graph ISTA: commit a local timing graph's results
// back into the global sta::Graph. This is the "write" half of the ISTA's
// surface (the "compute" half lives on LocalSta / PtGraph). Traversal order
// and *when* to commit are the optimizer's concern -- the ISTA only provides
// the primitive.

#pragma once

namespace lgista {

class PtGraph;

// Copy slew + arrival/required of the local graph's real (non-virtual)
// vertices back onto their backing sta::Vertex in the global graph. Only
// RefInput / RefOutput / SiblingLoad vertices are written (virtual buffer
// pins have no global vertex and are skipped). Requires that the PtGraph has
// been topologically sorted (sortedVertexIds()).
void updateTimingFromPtGraph(PtGraph* pt_graph);

}  // namespace lgista
