// SPDX-License-Identifier: BSD-3-Clause
//
// Correctness tests for the local graph ISTA. See TEST_PLAN.md.
//   t1   -- virtualReplaceCell local-vs-global equivalence
//   t2a  -- write-back identity (no cell change)
//   t2b  -- replace-all + topo-order write-back vs full global re-analysis
// `which` is one of "t1", "t2a", "t2b", "all". Returns 0 on PASS.

#pragma once

#include <string>

namespace sta {
class dbSta;
}

namespace lgista {

int runTest(sta::dbSta* sta, const std::string& which);

}  // namespace lgista
