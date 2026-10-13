// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "HTreeBuilder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Clock.h"
#include "Clustering.h"
#include "Cts3DDatabase.h"
#include "CtsObserver.h"
#include "SinkClustering.h"
#include "TarpClustering.h"
#include "TechChar.h"
#include "TreeBuilder.h"
#include "Util.h"
#include "odb/db.h"
#include "odb/isotropy.h"
#include "utl/Logger.h"

namespace dcts {

using utl::CTS;

namespace {

bool isNamedHardMacroMaster(const odb::dbMaster* master)
{
  if (master == nullptr) {
    return false;
  }
  std::string name = master->getName();
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return std::tolower(c);
  });
  return name.rfind("fakeram", 0) == 0 || name.rfind("fakeregfile", 0) == 0
         || name.rfind("sram_", 0) == 0;
}

bool isMacroBlockInst(const odb::dbInst* inst)
{
  if (inst == nullptr) {
    return false;
  }
  if (inst->isBlock()) {
    return true;
  }
  const odb::dbMaster* master = inst->getMaster();
  return (master != nullptr && master->isBlock())
         || isNamedHardMacroMaster(master);
}



// V81 KEYFIX (root cause of the dangling clock buffer / CTS-0418 / CTS-0726 class).
//
// mapLocationToSink_ is std::map<Point<double>,...> with an EXACT operator<, but every later
// lookup uses coordinates taken from topLevelSinksClustered_, a std::vector<std::pair<float,float>>
// -- i.e. narrowed to float.  A key computed straight from the double can therefore never be found.
//
// Writing `Point<double>(point.first, point.second)` on a local `const std::pair<float,float>`
// looks equivalent but is not: at -O3 gcc elides the double->float->double round trip, so the key
// went in as the full double (measured bits 0x4037f3dd2b91419d) while the vector held the narrowed
// float (0x4037f3dd20000000).  find() then missed, the leaf buffer was never attached to its parent
// branch net, its input pin stayed on no net, and the FFs under it silently lost their clock.
// `volatile` forces the narrowing to be materialised.
Point<double> narrowedKey(double x, double y)
{
  volatile float fx = static_cast<float>(x);
  volatile float fy = static_cast<float>(y);
  return Point<double>(static_cast<double>(fx), static_cast<double>(fy));
}

}  // namespace

// V73 clock-connectivity fix.
//
// Predicate for the Fix-G macro flat bypass, factored out so it can be evaluated
// BEFORE preSinkClustering() rather than after it.
//
// Root cause it fixes: Fix G skips H-tree construction entirely and wires the macro
// sinks straight to the (shared) root.  But preSinkClustering() ran first and, whenever
// the macro sink count exceeded min_clustering_macro_sinks_ (=10), it had already
// created leaf buffers + leaf subnets that own those macro CLK pins.  Nothing ever
// builds a branch net to drive those leaf buffers, so writeClockNetsToDb() emits them
// with an unconnected input pin.  STA then sees no clock at the macro CLK pins:
//   * no setup/hold checks on those macros, and
//   * vectorless power falls back to default activity, collapsing macro internal power
//     (the macro clk internal_power table is flat in slew, so power == energy x toggle).
// Observed on asap7_nangate45_3D/swerv_wrapper: 5 dangling t1_clkbuf_leaf_*_clk feeding
// 19 of 28 macro clk pins; also 7 dangling clkbuf_leaf_*_clk in the non-per-tier arm.
bool HTreeBuilder::macroFlatBypassApplies(
    int numSinks,
    const std::vector<const ClockInst*>& sinkInsts) const
{
  if (type_ != TreeType::MacroTree) {
    return false;
  }
  int macroFlatMax = 50;
  if (const char* e = std::getenv("CTS_MACRO_FLAT_MAX_SINKS")) {
    macroFlatMax = std::atoi(e);
  }
  if (macroFlatMax <= 0 || numSinks > macroFlatMax) {
    return false;
  }
  for (const ClockInst* sinkInst : sinkInsts) {
    if (sinkInst == nullptr || !isMacroBlockInst(sinkInst->getDbInst())) {
      return false;
    }
  }
  return !sinkInsts.empty();
}

// Removed anonymous namespace containing hasSuffix() and
// mapBufferMasterToTier() — moved to Cts3DDatabase::hasSuffix() and
// Cts3DDatabase::getBufferForTier()

Point<double> HTreeBuilder::legalizeOneBuffer(Point<double> bufferLoc,
                                              const std::string& bufferName)
{
  Point<double> legalLoc
      = TreeBuilder::legalizeOneBuffer(bufferLoc, bufferName);
  return resolveLocationCollision(legalLoc);
}

// Step 3: Skew-aware clustering via skew target penalty
// HB penalty removed — vias (including HB) don't get distance penalty in 2D CTS
// either, and STA already accounts for HB RC in timing. Clustering should remain
// purely spatial + skew-target-driven.
// Removed skewTargetBeta_ dead code block.
// computeDist() is called with candidate branching points (not in
// mapLocationToSink_) → find() always failed → beta penalty never fired.
// Target-aware clustering works via CTS_CKMEANS_TARGET_BETA in
// refineBranchingPointsWithClustering() → Clustering::calcDist().
double HTreeBuilder::computeDist(const Point<double>& x, const Point<double>& y)
{
  return TreeBuilder::computeDist(x, y);
}

// V58 Fix B: Epsilon-tolerant lookup for mapLocationToSink_.
// Falls back to linear scan within eps when exact find() fails.
// Matches the epsilon pattern already used in leaf-sink connection code.
ClockInst* HTreeBuilder::findSinkEps(const Point<double>& loc,
                                      double eps) const
{
  auto it = mapLocationToSink_.find(loc);
  if (it != mapLocationToSink_.end()) {
    return it->second;
  }
  auto hint = mapLocationToSink_.lower_bound(
      Point<double>(loc.getX() - eps, loc.getY() - eps));
  for (auto scan = hint; scan != mapLocationToSink_.end(); ++scan) {
    if (scan->first.getX() > loc.getX() + eps) break;
    if (std::abs(scan->first.getX() - loc.getX()) < eps
        && std::abs(scan->first.getY() - loc.getY()) < eps) {
      return scan->second;
    }
  }
  return nullptr;
}

void HTreeBuilder::preSinkClustering(
    const std::vector<std::pair<float, float>>& sinks,
    const std::vector<const ClockInst*>& sinkInsts,
    const float maxDiameter,
    const unsigned clusterSize,
    const bool secondLevel)
{
  bool maxDiameterSet = (type_ == TreeType::MacroTree)
                            ? options_->isMacroMaxDiameterSet()
                            : options_->isMaxDiameterSet();
  bool clusterSizeSet = (type_ == TreeType::MacroTree)
                            ? options_->isMacroSinkClusteringSizeSet()
                            : options_->isSinkClusteringSizeSet();

  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;

  const std::vector<std::pair<float, float>>& points = sinks;
  if (!secondLevel) {
    clock_.forEachSink([&](ClockInst& inst) {
      // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
      const Point<double> normLocation((double) inst.getX() / wireSegmentUnit_,
                                       (double) inst.getY() / wireSegmentUnit_);
      mapLocationToSink_[normLocation] = &inst;
      if (!fuzzyEqual(inst.getInsertionDelay(), 0.0, 1e-6)) {
        setSinkInsertionDelay(normLocation,
                              inst.getInsertionDelay() / wireSegmentUnit_);
        // clang-format off
	debugPrint(logger_, CTS, "clustering", 1, "sink {} has insDelay {} at {}",
		   inst.getName(), getSinkInsertionDelay(normLocation),
		   normLocation);
        // clang-format on
      }
    });
  }

  if (sinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    topLevelSinksClustered_ = sinks;
    return;
  }

  // V86_SEED: locals for the seeded (off-cluster + refinement) mode.
  std::vector<std::vector<unsigned>> seededSolution;
  bool haveSeeded = false;
  unsigned seedGroupSize = 0;
  float seedDiameter = 0.0f;

  // Hilbert-calibrated TARP.
  // Run Hilbert grid search as lightweight pre-pass to find optimal
  // cluster granularity (bestGroupSize + bestDiameter), then pass
  // to TARP FM. This ensures TARP operates at the right spatial
  // granularity for each design without per-design manual tuning.
  if ((enableTarp_ || tarpSeeded_) && !secondLevel
      && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
      && sinks.size() >= 50u) {

    // Step 1: Lightweight Hilbert grid search (no leaf buffer creation)
    SinkClustering calibration(options_, techChar_, this);
    const unsigned numPts = points.size();
    for (unsigned i = 0; i < numPts; ++i) {
      calibration.addPoint(points[i].first, points[i].second);
      if (sinkInsts[i]->getInputCap() == 0) {
        calibration.addCap(options_->getSinkBufferInputCap());
      } else {
        calibration.addCap(sinkInsts[i]->getInputCap());
      }
    }

    unsigned bestGS = 0;
    float bestDiam = 0.0f;
    // Grid search: groupSize ∈ {10,20,30} × diameter ∈ {50,100,200}μm
    for (unsigned diam : clusterDiameters()) {
      for (unsigned gs : options_->getSinkClusteringSizes()) {
        float maxDiam = diam * static_cast<float>(options_->getDbUnits())
                       / wireSegmentUnit_;
        calibration.run(gs, maxDiam, wireSegmentUnit_, bestGS, bestDiam);
      }
    }

    logger_->info(CTS, 845,
        "TARP Hilbert calibration: bestGroupSize={}, bestDiameter={:.0f}um",
        bestGS, bestDiam / options_->getDbUnits() * wireSegmentUnit_);

    // V86_SEED: seeded mode reuses the calibration grid-search result — which
    // IS the off clustering (same points/caps/grid as the off path's
    // both-unset branch) — as the seed, refines it (critical-pair moves only),
    // and falls through to the OFF consumption path below. This keeps cluster
    // order, tier split, V37 target split and sink registration byte-identical
    // to TARP-off when zero moves fire. runTarpClustering (FM bisection) is
    // NOT used in seeded mode — the bisection itself differs structurally from
    // greedy-scan (V85 CRIT_ONLY evidence) and cannot be made no-regression.
    // Guard: only the both-unset grid branch is replicated by the calibration
    // pre-pass, and macro trees keep the plain off path.
    if (tarpSeeded_ && !clusterSizeSet && !maxDiameterSet
        && type_ != TreeType::MacroTree) {
      seededSolution = calibration.sinkClusteringSolution();  // copy (seed)
      haveSeeded = true;
      seedGroupSize = (bestGS > 0) ? bestGS : 10;
      seedDiameter = bestDiam;
      // fall through to the off consumption path (no return)
    } else if (enableTarp_) {
      // Step 2: Pass calibrated params to TARP (env var overrides take priority)
      runTarpClustering(points, sinkInsts, bestGS, bestDiam);
      return;
    }
    // else: seeded-only requested but ineligible (macro tree / explicit
    // size|diameter flow) -> plain OFF path below (never FM).
  }

  SinkClustering matching(options_, techChar_, this);
  const unsigned numPoints = points.size();

  if (!haveSeeded) {
  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
    const std::pair<float, float>& point = points[pointIdx];
    matching.addPoint(point.first, point.second);
    if (sinkInsts[pointIdx]->getInputCap() == 0) {
      // Comes here in second level since first level buf cap is not set
      matching.addCap(options_->getSinkBufferInputCap());
    } else {
      matching.addCap(sinkInsts[pointIdx]->getInputCap());
    }
  }

  unsigned bestClusterSize = 0;
  float bestDiameter = 0.0;
  if (clusterSizeSet && maxDiameterSet) {
    // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize, maxDiameter, wireSegmentUnit_);
    // clang-format on
    matching.run(clusterSize,
                 maxDiameter,
                 wireSegmentUnit_,
                 bestClusterSize,
                 bestDiameter);
  } else if (!clusterSizeSet && maxDiameterSet) {
    // only diameter is set, try clustering sizes of 10, 20 and 30
    for (unsigned clusterSize2 : options_->getSinkClusteringSizes()) {
      // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize2, maxDiameter, wireSegmentUnit_);
      // clang-format on
      matching.run(clusterSize2,
                   maxDiameter,
                   wireSegmentUnit_,
                   bestClusterSize,
                   bestDiameter);
    }
  } else if (clusterSizeSet && !maxDiameterSet) {
    // only clustering size is set, try diameters of 50, 100 and 200 um
    for (unsigned clusterDiameter2 : options_->getSinkClusteringDiameters()) {
      // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize, clusterDiameter2, wireSegmentUnit_);
      // clang-format on
      float maxDiameter2 = clusterDiameter2 * (float) options_->getDbUnits()
                           / wireSegmentUnit_;
      matching.run(clusterSize,
                   maxDiameter2,
                   wireSegmentUnit_,
                   bestClusterSize,
                   bestDiameter);
    }
  } else {  // neighther clustering size nor diameter is set
    // try diameters of 50, 100 and 200 um
    for (unsigned clusterDiameter2 : clusterDiameters()) {
      // try clustering sizes of 10, 20 and 30
      for (unsigned clusterSize2 : options_->getSinkClusteringSizes()) {
        // clang-format off
        debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                   clusterSize2, clusterDiameter2, wireSegmentUnit_);
        // clang-format on
        float maxDiameter2 = clusterDiameter2 * (float) options_->getDbUnits()
                             / wireSegmentUnit_;
        matching.run(clusterSize2,
                     maxDiameter2,
                     wireSegmentUnit_,
                     bestClusterSize,
                     bestDiameter);
      }
    }
  }

  if (clusterSizeSet || maxDiameterSet) {
    logger_->info(
        CTS,
        204,
        "A clustering solution was found from clustering size of {} "
        "and clustering diameter of {:0.0f}.",
        bestClusterSize,
        std::round(bestDiameter / options_->getDbUnits() * wireSegmentUnit_));
    logger_->info(CTS,
                  205,
                  "Better solution may be possible if either "
                  "-sink_clustering_size, -sink_clustering_max_diameter, or "
                  "both options are omitted to enable automatic clustering.");
  } else {
    logger_->info(
        CTS,
        206,
        "Best clustering solution was found from clustering size of {} "
        "and clustering diameter of {:0.0f}.",
        bestClusterSize,
        std::round(bestDiameter / options_->getDbUnits() * wireSegmentUnit_));
  }
  }  // end if (!haveSeeded) — V86_SEED: matching untouched in seeded mode

  // V86_SEED: refine the off-clustering seed (critical-pair moves only), then
  // consume whichever solution applies through the SAME loop below so the
  // zero-move seeded case is byte-identical to TARP-off.
  if (haveSeeded) {
    refineSeedClusters(
        seededSolution, points, sinkInsts, seedGroupSize, seedDiameter);
  }
  const std::vector<std::vector<unsigned>>& clusteringSolution
      = haveSeeded ? seededSolution : matching.sinkClusteringSolution();

  unsigned clusterCount = 0;

  // Tier-aware clustering: detect cross-tier clusters
  // Calculate HB equivalent distance for cross-tier penalty awareness
  const double wire_res_per_unit = cts3dDb_->getResPerDBU(0) * wireSegmentUnit_;
  const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(0) * wireSegmentUnit_;
  const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
      wire_res_per_unit, wire_cap_per_unit);

  // Note: Current k-means clustering uses geometric distance only.
  // Possible extension: inject HB equivalent distance into SinkClustering/CKMeans
  // as penalty when computing distance between sinks on different tiers.
  // This would actively discourage cross-tier clusters during k-means.

  std::vector<std::pair<float, float>> newSinkLocations;

  // Helper lambda to create one leaf buffer for a group of FFs.
  // This is factored out so that sub-clustering can call it multiple times per cluster.
  int splitCount_2way = 0, splitCount_3way = 0;  // sub-cluster stats
  auto createLeafBuffer = [&](const std::vector<ClockInst*>& subInsts,
                               const std::vector<unsigned>& subPointIdxs,
                               unsigned& leafIdx) {
    if (subInsts.empty()) return;
    float xSum = 0, ySum = 0;
    for (unsigned pi : subPointIdxs) {
      xSum += points[pi].first;
      ySum += points[pi].second;
    }
    const float normCenterX = xSum / subPointIdxs.size();
    const float normCenterY = ySum / subPointIdxs.size();
    Point<double> center((double) normCenterX, (double) normCenterY);

    const int target_tier = cts3dDb_->getDominantTier(subInsts);
    const std::string sink_buffer
        = cts3dDb_->getBufferForTier(options_->getSinkBuffer(), target_tier);

    const std::string baseName = tierPrefix_ + (secondLevel ? "clkbuf_leaf2_" : "clkbuf_leaf_");
    Point<double> rootBufLoc = legalizeOneBuffer(center, sink_buffer);
    commitMoveLoc(center, rootBufLoc);

    ClockInst& rootBuffer = clock_.addClockBuffer(
        baseName + std::to_string(leafIdx),
        sink_buffer,
        rootBufLoc.getX() * wireSegmentUnit_,
        rootBufLoc.getY() * wireSegmentUnit_);
    cts3dDb_->setClockInstTier(rootBuffer, target_tier);

    if (!secondLevel) {
      addFirstLevelSinkDriver(&rootBuffer);
    } else {
      addSecondLevelSinkDriver(&rootBuffer);
    }

    const std::string netBaseName = tierPrefix_ + (secondLevel ? "clknet_leaf2_" : "clknet_leaf_");
    ClockSubNet& clockSubNet
        = clock_.addSubNet(netBaseName + std::to_string(leafIdx));
    clockSubNet.addInst(rootBuffer);
    for (ClockInst* ci : subInsts) {
      clockSubNet.addInst(*ci);
    }
    if (!secondLevel) {
      clockSubNet.setLeafLevel(true);
    }

    const std::pair<float, float> point(rootBufLoc.getX(), rootBufLoc.getY());
    newSinkLocations.emplace_back(point);
    // V81 KEYFIX: key must be the float-narrowed coordinate (see narrowedKey()).
    const Point<double> mapKey = narrowedKey(rootBufLoc.getX(), rootBufLoc.getY());
    mapLocationToSink_[mapKey] = &rootBuffer;

    // Store mean LP skew target for this leaf buffer
    if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
      double clusterSum = 0.0;
      for (ClockInst* ci : subInsts) {
        std::string nm = ci->getName();
        const auto sp = nm.rfind('/');
        if (sp != std::string::npos) nm = nm.substr(0, sp);
        clusterSum += cts3dDb_->getSkewTarget(nm);
      }
      clusterBufTarget_[&rootBuffer] = clusterSum / subInsts.size();
    }
    ++leafIdx;
  };

  for (const std::vector<unsigned>& cluster : clusteringSolution) {
    if (cluster.size() == 1) {
      const std::pair<float, float>& point = points[cluster[0]];
      newSinkLocations.emplace_back(point);
      // V58 Fix A: Register singleton cluster in mapLocationToSink_.
      // Without this, the singleton's coordinate is in topLevelSinksClustered_
      // but missing from mapLocationToSink_ → computeBranchDelayTargets() and
      // getDominantTier() silently skip this sink (exact find() fails).
      // Same fix already exists in TARP path (line ~4819).
      const Point<double> sinkPt((double) point.first, (double) point.second);
      if (mapLocationToSink_.find(sinkPt) == mapLocationToSink_.end()) {
        constexpr double eps = 0.01;
        for (auto& [key, sink] : mapLocationToSink_) {
          if (std::abs(key.getX() - sinkPt.getX()) < eps
              && std::abs(key.getY() - sinkPt.getY()) < eps) {
            mapLocationToSink_[sinkPt] = sink;
            break;
          }
        }
      }
      clusterCount++;
      continue;
    }
    if (cluster.size() > 1) {
      // Collect ClockInst pointers for this cluster
      std::vector<ClockInst*> clusterClockInsts;
      for (auto point_idx : cluster) {
        const std::pair<double, double>& point = points[point_idx];
        const Point<double> mapPoint(point.first, point.second);
        if (mapLocationToSink_.find(mapPoint) == mapLocationToSink_.end()) {
          // FloatFix: epsilon fallback for float→double key mismatch
          bool found = false;
          constexpr double eps = 0.01;
          for (auto& [key, sink] : mapLocationToSink_) {
            if (std::abs(key.getX() - mapPoint.getX()) < eps
             && std::abs(key.getY() - mapPoint.getY()) < eps) {
              mapLocationToSink_[mapPoint] = sink;
              found = true;
              break;
            }
          }
          if (!found) {
            // V58: eps=1e-4 failed — find absolute nearest entry (no distance limit).
            // Skipping would leave this sink unconnected → SIGSEGV in writeDummyLoadsToDb.
            double bestDist = std::numeric_limits<double>::max();
            for (auto& [key, sink] : mapLocationToSink_) {
              double d = std::abs(key.getX() - mapPoint.getX())
                       + std::abs(key.getY() - mapPoint.getY());
              if (d < bestDist) {
                bestDist = d;
                mapLocationToSink_[mapPoint] = sink;
                found = true;
              }
            }
            if (found) {
              logger_->warn(CTS, 79,
                  "Sink at ({:.6f}, {:.6f}) — nearest match dist={:.6f}",
                  mapPoint.getX(), mapPoint.getY(), bestDist);
            }
          }
        }
        clusterClockInsts.push_back(mapLocationToSink_[mapPoint]);
      }

      // Detect cross-tier clustering
      int tier0_count = 0, tier1_count = 0;
      for (const ClockInst* inst : clusterClockInsts) {
        const int inst_tier = cts3dDb_->getInstTier(inst->getDbInst());
        if (inst_tier == 0) tier0_count++;
        else if (inst_tier == 1) tier1_count++;
      }
      if (tier0_count > 0 && tier1_count > 0) {
        logger_->warn(CTS, 571,
                      "Cross-tier cluster detected: cluster={}, tier0={}, tier1={}, "
                      "dominant_tier={}, HB_penalty={:.1f}um",
                      clusterCount, tier0_count, tier1_count,
                      cts3dDb_->getDominantTier(clusterClockInsts),
                      hb_equivalent_dist);
      }

      // Lambda for V37 target-aware split + leaf buffer
      // creation.  Applies target-spread splitting to any (sub-)cluster, then
      // creates leaf buffers.  Used by both tier-split sub-clusters and
      // unsplit clusters so that tier split no longer skips target split.
      auto targetSplitOrCreate = [&](std::vector<ClockInst*>& insts,
                                     const std::vector<unsigned>& idxs) {
        if (enableTargetSplit_ && insts.size() > 2
            && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
          struct FFTarget {
            unsigned pointIdx;
            ClockInst* inst;
            double target;
          };
          std::vector<FFTarget> ffTargets;
          ffTargets.reserve(insts.size());
          double tMin = std::numeric_limits<double>::max();
          double tMax = -std::numeric_limits<double>::max();
          for (size_t i = 0; i < insts.size(); ++i) {
            std::string nm = insts[i]->getName();
            const auto sp = nm.rfind('/');
            if (sp != std::string::npos) nm = nm.substr(0, sp);
            const double t = cts3dDb_->getSkewTarget(nm);
            ffTargets.push_back({idxs[i], insts[i], t});
            tMin = std::min(tMin, t);
            tMax = std::max(tMax, t);
          }
          const double spread = tMax - tMin;

          if (spread > splitThreshold2wayNs_) {
            std::sort(ffTargets.begin(), ffTargets.end(),
                      [](const FFTarget& a, const FFTarget& b) {
                        return a.target < b.target;
                      });

            int numSubs = 1;
            if (spread > splitThreshold3wayNs_ && ffTargets.size() >= 6) {
              numSubs = 3;
              ++splitCount_3way;
            } else if (spread > splitThreshold2wayNs_ && ffTargets.size() >= 4) {
              numSubs = 2;
              ++splitCount_2way;
            }

            if (numSubs > 1) {
              const size_t total = ffTargets.size();
              for (int s = 0; s < numSubs; ++s) {
                const size_t start = s * total / numSubs;
                const size_t end = (s + 1) * total / numSubs;
                if (start >= end) continue;

                std::vector<ClockInst*> splitInsts;
                std::vector<unsigned> splitIdxs;
                for (size_t k = start; k < end; ++k) {
                  splitInsts.push_back(ffTargets[k].inst);
                  splitIdxs.push_back(ffTargets[k].pointIdx);
                }
                if (splitInsts.size() < 2 && s + 1 < numSubs) {
                  continue;
                }

                logger_->info(CTS, 413,
                              "Sub-cluster {}.{}: {} sinks, "
                              "target=[{:.4f}, {:.4f}]ns (spread={:.4f}ns)",
                              clusterCount, s, splitInsts.size(),
                              ffTargets[start].target,
                              ffTargets[end - 1].target,
                              ffTargets[end - 1].target
                                  - ffTargets[start].target);
                createLeafBuffer(splitInsts, splitIdxs, clusterCount);
              }
              return;
            }
          }
        }
        // No target split needed: create leaf buffer directly
        createLeafBuffer(insts, idxs, clusterCount);
      };

      // Tier-aware cluster split.
      // If cluster has both bottom and upper FFs, split into two sub-clusters.
      // Each sub-cluster also gets V37 target-aware splitting if enabled.
      bool wasSplit = false;
      if (enableTierSplit_ && tier0_count > 0 && tier1_count > 0) {
        wasSplit = true;

        // Partition into bottom and upper sub-clusters
        std::vector<ClockInst*> bottomInsts, upperInsts;
        std::vector<unsigned> bottomIdxs, upperIdxs;
        for (size_t i = 0; i < cluster.size(); ++i) {
          const int inst_tier = cts3dDb_->getInstTier(
              clusterClockInsts[i]->getDbInst());
          if (inst_tier == 0) {
            bottomInsts.push_back(clusterClockInsts[i]);
            bottomIdxs.push_back(cluster[i]);
          } else {
            upperInsts.push_back(clusterClockInsts[i]);
            upperIdxs.push_back(cluster[i]);
          }
        }

        logger_->info(CTS, 377,
            "Tier split cluster {}: bottom={}, upper={}",
            clusterCount,
            bottomInsts.size(), upperInsts.size());

        // Apply V37 target split to each tier sub-cluster
        if (!bottomInsts.empty()) targetSplitOrCreate(bottomInsts, bottomIdxs);
        if (!upperInsts.empty()) targetSplitOrCreate(upperInsts, upperIdxs);
      }

      // Non-tier-split: apply V37 target split or direct leaf buffer
      if (!wasSplit) {
        targetSplitOrCreate(clusterClockInsts, cluster);
      }
    }
  }

  // Log sub-clustering statistics
  if (splitCount_2way > 0 || splitCount_3way > 0) {
    logger_->info(CTS, 414,
                  "Target-aware sub-clustering: {} 2-way splits, "
                  "{} 3-way splits, total leaf buffers = {}",
                  splitCount_2way, splitCount_3way, clusterCount);
  }

  topLevelSinksClustered_ = std::move(newSinkLocations);
  if (clusterCount) {
    treeBufLevels_++;
  }

  logger_->info(CTS,
                19,
                " Total number of sinks after clustering: {}.",
                topLevelSinksClustered_.size());

  // Verify every entry in topLevelSinksClustered_ exists in mapLocationToSink_
  int missingCount = 0;
  for (size_t si = 0; si < topLevelSinksClustered_.size(); ++si) {
    const auto& sinkPair = topLevelSinksClustered_[si];
    Point<double> key(sinkPair.first, sinkPair.second);
    if (mapLocationToSink_.find(key) == mapLocationToSink_.end()) {
      logger_->warn(CTS, 418,
                    "Debug: topLevelSinksClustered_[{}] = ({:.10f}, {:.10f}) "
                    "NOT in mapLocationToSink_ (size={})",
                    si, sinkPair.first, sinkPair.second,
                    mapLocationToSink_.size());
      ++missingCount;
    }
  }
  if (missingCount > 0) {
    logger_->warn(CTS, 419,
                  "Debug: {} / {} sinks missing from mapLocationToSink_!",
                  missingCount, topLevelSinksClustered_.size());
  } else {
    logger_->info(CTS, 427,
                  "Debug: All {} sinks verified in mapLocationToSink_ (mapSize={})",
                  topLevelSinksClustered_.size(), mapLocationToSink_.size());
  }
}

Point<double> HTreeBuilder::resolveLocationCollision(
    const Point<double>& legalCenter) const
{
  Point<double> resolvedLocation = legalCenter;

  // Collision check and jittering to ensure unique coordinates
  unsigned jitterCount = 0;
  while (true) {
    // FloatFix: removed (float) cast to match map key precision (double)
    Point<double> checkKey(resolvedLocation.getX(),
                           resolvedLocation.getY());
    if (mapLocationToSink_.find(checkKey) == mapLocationToSink_.end()) {
      break;
    }
    jitterCount++;
    double offset = (double) jitterCount / wireSegmentUnit_;
    resolvedLocation.setX(legalCenter.getX() + offset);
  }

  if (jitterCount > 0) {
    debugPrint(
        logger_,
        utl::CTS,
        "clustering",
        1,
        "Resolved collision via jittering ({} times): ({}, {}) -> ({}, {})",
        jitterCount,
        legalCenter.getX(),
        legalCenter.getY(),
        resolvedLocation.getX(),
        resolvedLocation.getY());
  }
  return resolvedLocation;
}

// Removed getDominantTierFromInsts, getDominantTierFromSinkLocs,
// getDominantTierFromClockSinks — moved to Cts3DDatabase::getDominantTier() and
// Cts3DDatabase::getDominantTierFromClock()

void HTreeBuilder::initSinkRegion()
{
  const unsigned wireSegmentUnitInDbu = techChar_->getLengthUnit();
  const int dbUnits = options_->getDbUnits();
  wireSegmentUnit_ = wireSegmentUnitInDbu;

  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;

  logger_->info(CTS,
                20,
                " Wire segment unit: {}  dbu ({} um).",
                wireSegmentUnit_,
                wireSegmentUnitInDbu / dbUnits);

  if (options_->isSimpleSegmentEnabled()) {
    const int remainingLength
        = options_->getBufferDistance() / (wireSegmentUnitInDbu * 2);
    logger_->info(CTS,
                  21,
                  " Distance between buffers: {} units ({} um).",
                  remainingLength,
                  static_cast<int>(options_->getBufferDistance() / dbUnits));
    if (options_->isVertexBuffersEnabled()) {
      const int vertexBufferLength
          = options_->getVertexBufferDistance() / (wireSegmentUnitInDbu * 2);
      logger_->info(
          CTS,
          22,
          " Branch length for Vertex Buffer: {} units ({} um).",
          vertexBufferLength,
          static_cast<int>(options_->getVertexBufferDistance() / dbUnits));
    }
  }

  std::vector<std::pair<float, float>> topLevelSinks;
  std::vector<const ClockInst*> sinkInsts;
  initTopLevelSinks(topLevelSinks, sinkInsts);

  const float maxDiameter = (clusterDiameter * dbUnits) / wireSegmentUnit_;
  // clang-format off
  debugPrint(logger_, CTS, "clustering", 1, "maxDiameter={:0.3f} = "
             "origMaxDiam={} * dbUnits={} / wireSegmentUnit_={}",
             maxDiameter, clusterDiameter, dbUnits,
             wireSegmentUnit_);
  // clang-format on

  // V73 fix: same guard as runSingleTierTree().  run()'s Fix-G bypass calls
  // createSingleBufferClockNet() and never drives leaf buffers, so clustering here
  // would leave them dangling (7 such buffers observed in the non-per-tier arm).
  const bool flatMacroBypass
      = macroFlatBypassApplies(static_cast<int>(topLevelSinks.size()), sinkInsts);
  if (flatMacroBypass) {
    logger_->info(CTS, 723,
        "V73: macro flat bypass detected before clustering - "
        "skipping preSinkClustering ({} macro sinks connect directly), "
        "sink region from raw sinks",
        topLevelSinks.size());
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
    return;
  }

  preSinkClustering(topLevelSinks, sinkInsts, maxDiameter, clusterSize);
  if (topLevelSinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    logger_->info(CTS, 23, " Original sink region: {}.", sinkRegionDbu);

    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
  } else {
    if (topLevelSinksClustered_.size() > 400
        && options_->getSinkClusteringLevels() > 0) {
      std::vector<std::pair<float, float>> secondLevelLocs;
      std::vector<const ClockInst*> secondLevelInsts;
      initSecondLevelSinks(secondLevelLocs, secondLevelInsts);
      preSinkClustering(secondLevelLocs,
                        secondLevelInsts,
                        maxDiameter * 4,
                        std::ceil(std::sqrt(clusterSize)),
                        true);
    }
    sinkRegion_ = clock_.computeSinkRegionClustered(topLevelSinksClustered_);
  }
  logger_->info(CTS, 24, " Normalized sink region: {}.", sinkRegion_);
  logger_->info(CTS, 25, "    Width:  {:.4f}.", sinkRegion_.getWidth());
  logger_->info(CTS, 26, "    Height: {:.4f}.", sinkRegion_.getHeight());
}

static void plotBlockage(std::ofstream& file,
                         odb::dbDatabase* db_,
                         int scalingFactor)
{
  unsigned i = 0;
  for (odb::dbBlockage* blockage : db_->getChip()->getBlock()->getBlockages()) {
    odb::dbBox* bbox = blockage->getBBox();
    int x = bbox->xMin() / scalingFactor;
    int y = bbox->yMin() / scalingFactor;
    int w = bbox->xMax() / scalingFactor - bbox->xMin() / scalingFactor;
    int h = bbox->yMax() / scalingFactor - bbox->yMin() / scalingFactor;
    file << i++ << " " << x << " " << y << " " << w << " " << h
         << " block  scalingFactor=";
    file << scalingFactor << " " << blockage->getId() << '\n';
  }
}

// distance to move sinks from old loc to new loc
double HTreeBuilder::weightedDistance(const Point<double>& newLoc,
                                      const Point<double>& oldLoc,
                                      const std::vector<Point<double>>& sinks)
{
  double dist = 0;
  for (const Point<double>& sink : sinks) {
    dist += computeDist(newLoc, sink);
    dist += computeDist(newLoc, oldLoc);
  }
  return dist;
}

static void plotSinks(std::ofstream& file,
                      const std::vector<Point<double>>& sinks)
{
  unsigned cnt = 0;
  for (const Point<double>& pt : sinks) {
    double x = pt.getX();
    double y = pt.getY();
    double w = 1;
    double h = 1;
    auto name = "sink_";
    file << cnt++ << " " << x << " " << y << " " << w << " " << h;
    file << " " << name << '\n';
  }
}

unsigned HTreeBuilder::findSibling(LevelTopology& topology,
                                   unsigned i,
                                   unsigned par)
{
  for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
    unsigned k = topology.getBranchingPointParentIdx(idx);
    if (idx != i && k == par) {
      return idx;
    }
  }
  return i;
}

void HTreeBuilder::scalePosition(Point<double>& loc,
                                 const Point<double>& parLoc,
                                 double leng,
                                 double scale)
{
  double px = parLoc.getX();
  double py = parLoc.getY();
  double ax = loc.getX();
  double ay = loc.getY();

  double d = computeDist(loc, parLoc);
  double x, y;
  if (d > 0) {  // yy8
    double delta = d * scale;
    double dx = ax - px;
    double dy = ay - py;
    dx += (dx > 0) ? delta : -delta;
    dy += (dy > 0) ? -delta : delta;
    double scale = leng / d;
    x = px + dx * scale;
    y = py + dy * scale;
  } else {
    x = px + leng / 2;
    y = py + leng / 2;
  }
  loc.setX(x);
  loc.setY(y);
}

static void setSiblingPosition(const Point<double>& a,
                               Point<double>& b,
                               const Point<double>& parLoc)
{
  double px = parLoc.getX();
  double py = parLoc.getY();
  double ax = a.getX();
  double ay = a.getY();
  double bx = 2 * px - ax;
  double by = 2 * py - ay;
  b.setX(bx);
  b.setY(by);
}

// Balance the two branches on the very top level
void HTreeBuilder::adjustToplevelTopology(Point<double>& a,
                                          Point<double>& b,
                                          const Point<double>& parLoc)
{
  double da = computeDist(a, parLoc);
  double db = computeDist(b, parLoc);
  if (da < db) {
    setSiblingPosition(a, b, parLoc);
  } else {
    setSiblingPosition(b, a, parLoc);
  }
}

void HTreeBuilder::findLegalLocations(const Point<double>& parentPoint,
                                      const Point<double>& branchPoint,
                                      double x1,
                                      double y1,
                                      double x2,
                                      double y2,
                                      std::vector<Point<double>>& points)
{
  // add 4 corners of blockage
  addCandidateLoc(x1, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x1, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, y1, parentPoint, x1, y1, x2, y2, points);

  // add straightline neighbors
  double bx = branchPoint.getX();
  double by = branchPoint.getY();
  addCandidateLoc(bx, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(bx, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x1, by, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, by, parentPoint, x1, y1, x2, y2, points);

  double px = parentPoint.getX();
  double py = parentPoint.getY();
  double dx = px - bx;
  double dy = py - by;
  double m = dy / dx;
  // y = m*(x-bx) + by
  addCandidateLoc(x1, m * (x1 - bx) + by, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, m * (x2 - bx) + by, parentPoint, x1, y1, x2, y2, points);
  // x = (y-by)/m + bx
  addCandidateLoc((y1 - by) / m + bx, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc((y2 - by) / m + bx, y2, parentPoint, x1, y1, x2, y2, points);
  // clang-format off
  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("    branchPt:{} is not legal, parentPt:{} blockages:({:0.3f} {:0.3f}) "
        "({:0.3f} {:0.3f})", branchPoint, parentPoint, x1, y1, x2, y2);
    for (Point<double> point : points) {
      logger_->report("      FLL candiate {}", point);
    }
    // clang-format on
  }
}

Point<double> HTreeBuilder::findBestLegalLocation(
    double targetDist,
    const Point<double>& branchPoint,
    const Point<double>& parentPoint,
    const std::vector<Point<double>>& legalLocations,
    const std::vector<Point<double>>& sinks,
    double x1,
    double y1,
    double x2,
    double y2,
    int scalingFactor,
    odb::Direction2D direction)
{
  Point<double> best(0.0, 0.0);
  double minDiff = std::numeric_limits<double>::max();
  for (const Point<double>& loc : legalLocations) {
    double dist = computeDist(loc, parentPoint);
    dist += computeDist(loc, branchPoint);
    double diff = abs(dist - targetDist);
    // clang-format off
    if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
      logger_->report("      Loc {}: curr dist={:0.3f} target dist={:0.3f} sink"
		      " dist={:0.3f}", loc, dist, targetDist,
		      weightedDistance(loc, branchPoint, sinks));
    }
    // clang-format on
    if (diff < minDiff) {
      minDiff = diff;
      best = loc;
    }
  }

  return adjustBestLegalLocation(targetDist,
                                 best,
                                 parentPoint,
                                 sinks,
                                 x1,
                                 y1,
                                 x2,
                                 y2,
                                 scalingFactor,
                                 direction);
}

// Adjust buffer location in two steps:
// step1: try moving buffer to existing blockage boundary (less expensive)
// step2: try moving buffer beyond existing blockage boundary (more expensive)
// In both steps, the first priority is to match target distance from current
// point to parent point.  The second priority is to lower the weighted sink
// distance
Point<double> HTreeBuilder::adjustBestLegalLocation(
    double targetDist,
    const Point<double>& currLoc,
    const Point<double>& parentPoint,
    const std::vector<Point<double>>& sinks,
    double x1,
    double y1,
    double x2,
    double y2,
    int scalingFactor,
    odb::Direction2D direction)
{
  if (fuzzyEqual(targetDist, computeDist(currLoc, parentPoint))) {
    return currLoc;
  }

  // try moving along blockage boundary
  Point<double> bestLoc = currLoc;
  if (adjustAlongBlockage(targetDist,
                          currLoc,
                          parentPoint,
                          sinks,
                          x1,
                          y1,
                          x2,
                          y2,
                          scalingFactor,
                          bestLoc)) {
    return bestLoc;
  }

  // try moving beyond blockage boundary
  return adjustBeyondBlockage(
      currLoc, parentPoint, targetDist, sinks, scalingFactor, direction);
}

bool HTreeBuilder::adjustAlongBlockage(double targetDist,
                                       const Point<double>& currLoc,
                                       const Point<double>& parentPoint,
                                       const std::vector<Point<double>>& sinks,
                                       double x1,
                                       double y1,
                                       double x2,
                                       double y2,
                                       int scalingFactor,
                                       Point<double>& bestLoc)
{
  Point<double> newLoc = currLoc;
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "{} currDist={:0.3f} != "
	     "targetDist={:0.3f}, adjustAlongBlockage...", currLoc,
	     computeDist(currLoc, parentPoint), targetDist);
  // clang-format on
  double x = currLoc.getX();
  double y = currLoc.getY();
  double px = parentPoint.getX();
  double py = parentPoint.getY();

  Point<double> noBest(std::numeric_limits<double>::min(),
                       std::numeric_limits<double>::min());
  bestLoc = noBest;
  double bestSinkDist = std::numeric_limits<double>::max();
  double sinkDist = 0.0;
  std::vector<Point<double>> candidates;

  // move along y axis within blockage
  double newX = x;
  double newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial y#1

  newY = py + targetDist - abs(px - newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial y#2

  // move along x axis within blockage
  newY = y;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial x#1

  newX = px + targetDist - abs(py - newY);
  newLoc.setX(newX);
  candidates.emplace_back(newLoc);  // trial x#2

  for (const Point<double>& candidate : candidates) {
    checkLegalityAndCostSpecial(currLoc,
                                candidate,
                                parentPoint,
                                targetDist,
                                sinks,
                                scalingFactor,
                                x1,
                                y1,
                                x2,
                                y2,
                                bestLoc,
                                sinkDist,
                                bestSinkDist);
  }

  return (bestLoc != noBest);
}

void HTreeBuilder::checkLegalityAndCostSpecial(
    const Point<double>& oldLoc,
    const Point<double>& newLoc,
    const Point<double>& parentPoint,
    double targetDist,
    const std::vector<Point<double>>& sinks,
    int scalingFactor,
    double x1,
    double y1,
    double x2,
    double y2,
    Point<double>& bestLoc,
    double& sinkDist,
    double& bestSinkDist)
{
  if (fuzzyEqual(computeDist(newLoc, parentPoint), targetDist)
      && checkLegalitySpecial(newLoc, x1, y1, x2, y2, scalingFactor)) {
    sinkDist = weightedDistance(newLoc, oldLoc, sinks);
    if (sinkDist < bestSinkDist) {
      bestLoc = newLoc;
      bestSinkDist = sinkDist;
    }
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBestLegalLoc: branchPt "
	       "move:{}=>{} is legal, dist={:0.3f}, sinkDist={:0.3f}",
	       oldLoc, newLoc, targetDist, sinkDist);
    // clang-format on
  }
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "adjustBestLegalLoc: branchPt "
	     "move:{}=>{} is illegal or dist {:0.3f} != {:0.3f}",
	     oldLoc, newLoc, computeDist(newLoc, parentPoint), targetDist);
  // clang-format on
}

// 1) Branch point couldn't be legalized by simply moving it along blockage
//    boundary, or
// 2) Branch point is legal but parent point has moved, so it is necessary to
//    move it to match target topology length
// In both cases, move branch point to match target distance and minimize
// weighted sink distance
Point<double> HTreeBuilder::adjustBeyondBlockage(
    const Point<double>& branchPoint,
    const Point<double>& parentPoint,
    double targetDist,
    const std::vector<Point<double>>& sinks,
    int scalingFactor,
    odb::Direction2D direction)
{
  double px = parentPoint.getX();
  double py = parentPoint.getY();
  std::vector<Point<double>> candidates;
  Point<double> point = branchPoint;

  // try points that are on edges of "Manhattan square"
  // with the parent point in the center
  //              p16 p1 p9
  //            p8         p5
  //          p15            p10
  //        p4        pp       p2
  //         p14             p11
  //           p7           p6
  //              p13 p3 p12
  double leng50 = targetDist * 0.5;
  double leng25 = targetDist * 0.25;
  double leng75 = targetDist * 0.75;
  switch (direction) {
    case odb::Direction2D::North:
      addCandidatePoint(px, py + targetDist, point, candidates);       // p1
      addCandidatePoint(px + leng50, py + leng50, point, candidates);  // p5
      addCandidatePoint(px - leng50, py + leng50, point, candidates);  // p8
      addCandidatePoint(px + leng25, py + leng75, point, candidates);  // p9
      addCandidatePoint(px - leng25, py + leng75, point, candidates);  // p16
      break;
    case odb::Direction2D::East:
      addCandidatePoint(px + targetDist, py, point, candidates);       // p2
      addCandidatePoint(px + leng50, py + leng50, point, candidates);  // p5
      addCandidatePoint(px + leng50, py - leng50, point, candidates);  // p6
      addCandidatePoint(px + leng75, py + leng25, point, candidates);  // p10
      addCandidatePoint(px + leng75, py - leng25, point, candidates);  // p11
      break;
    case odb::Direction2D::South:
      addCandidatePoint(px, py - targetDist, point, candidates);       // p3
      addCandidatePoint(px + leng50, py - leng50, point, candidates);  // p6
      addCandidatePoint(px - leng50, py - leng50, point, candidates);  // p7
      addCandidatePoint(px + leng25, py - leng75, point, candidates);  // p12
      addCandidatePoint(px - leng25, py - leng75, point, candidates);  // p13
      break;
    default:
      addCandidatePoint(px - targetDist, py, point, candidates);       // p4
      addCandidatePoint(px - leng50, py - leng50, point, candidates);  // p7
      addCandidatePoint(px - leng50, py + leng50, point, candidates);  // p8
      addCandidatePoint(px - leng75, py - leng25, point, candidates);  // p14
      addCandidatePoint(px - leng75, py + leng25, point, candidates);  // p15
      // odb::Direction2D::West
      break;
  }

  // check if any corners of Manhanttan square are inside some blockage
  // if so add candidate points that intersect blockage and Manhattan square
  // p1 is top corner

  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    point.setX(px);
    point.setY(py + targetDist);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::North);
  }

  // p2 is right corner
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    point.setX(px + targetDist);
    point.setY(py);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::East);
  }

  // p3 is bottom corner
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    point.setX(px);
    point.setY(py - targetDist);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::South);
  }

  // p4 is left corner
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    point.setX(px - targetDist);
    point.setY(py);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::West);
  }

  // try moving cell along x or y, with some offset
  double bx = branchPoint.getX();
  double by = branchPoint.getY();
  Point<double> noBest(std::numeric_limits<double>::min(),
                       std::numeric_limits<double>::min());
  Point<double> bestLoc = noBest;
  double bestSinkDist = std::numeric_limits<double>::max();
  double sinkDist = 0.0;
  // get small x and y offset to find "channel" through blockages
  double minX = 10.0 * getBufferWidth();
  double minY = 10.0 * getBufferHeight();

  // try small offset in x direction
  Point<double> newLoc(bx, by);
  double newX = bx + minX;
  double newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    candidates.emplace_back(newLoc);  // trial x#1
  }

  newLoc.setY(py + targetDist - abs(px - newX));
  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    candidates.emplace_back(newLoc);  // trial x#2
  }

  newX = bx - minX;
  newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    candidates.emplace_back(newLoc);  // trial x#3
  }

  newLoc.setY(py + targetDist - abs(px - newX));
  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    candidates.emplace_back(newLoc);  // trial x#4
  }

  // try small offset in y direction
  newY = by + minY;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    candidates.emplace_back(newLoc);  // trial y#1
  }

  newLoc.setX(px + targetDist - abs(py - newY));
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    candidates.emplace_back(newLoc);  // trial y#2
  }

  newY = by - minY;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    candidates.emplace_back(newLoc);  // trial y#3
  }

  newLoc.setX(px + targetDist - abs(py - newY));
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    candidates.emplace_back(newLoc);  // trial y#4
  }

  for (const Point<double>& candidate : candidates) {
    checkLegalityAndCost(branchPoint,
                         candidate,
                         parentPoint,
                         targetDist,
                         sinks,
                         scalingFactor,
                         bestLoc,
                         sinkDist,
                         bestSinkDist);
  }

  if (bestLoc != noBest) {
    return bestLoc;
  }
  return branchPoint;
}

void HTreeBuilder::addCandidatePointsAlongBlockage(
    const Point<double>& point,
    const Point<double>& parentPoint,
    double targetDist,
    int scalingFactor,
    std::vector<Point<double>>& candidates,
    odb::Direction2D direction)
{
  double x1, y1, x2, y2;
  if (findBlockage(point, scalingFactor, x1, y1, x2, y2)) {
    Point<double> point2 = point;
    double px = parentPoint.getX();
    double py = parentPoint.getY();
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "  {} corner {} of Manhattan "
	       "sqaure is inside blockage ({:0.3f} {:0.3f}) ({:0.3f} {:0.3f})",
	       direction, point, x1, y1, x2, y2);
    // clang-format on
    switch (direction) {
      case odb::Direction2D::North:
        //        ----------
        //        |        |
        //        |   top  |
        //        |   / \  |
        // (x1, y1)--x---x-(x2, y1)
        //          /  .
        //          (px, py)
        //
        // (px - x) + (y1 - py) = targetDist
        // (x - px) + (y1 - py) = targetDist
        point2.setX(-targetDist - py + px + y1);
        point2.setY(y1);
        candidates.emplace_back(point2);
        point2.setY(targetDist + py + px - y1);
        candidates.emplace_back(point2);
        break;
      case odb::Direction2D::East:
        // (x1, y2)---------
        //       \ |        |
        //         x        |
        //         |\       |
        //   .     | right  |
        // (px,py)| /       |
        //         x        |
        // (x1, y1)---------
        // (x1 - px) + (y - py) = targetDist
        // (x1 - px) + (py - y) = targetDist
        point2.setX(x1);
        point2.setY(targetDist + py + px - x1);
        candidates.emplace_back(point2);
        point2.setY(-targetDist + py - px + x1);
        candidates.emplace_back(point2);
        break;
      case odb::Direction2D::South:
        point2.setX(-targetDist + py + px - y2);
        point2.setY(y2);
        candidates.emplace_back(point2);
        point2.setX(targetDist - py + px + y2);
        candidates.emplace_back(point2);
        break;
      default:
        // odb::Direction2D::West
        point2.setX(x2);
        point2.setY(targetDist + py - px + x2);
        candidates.emplace_back(point2);
        point2.setY(-targetDist + py + px - x2);
        candidates.emplace_back(point2);
        break;
    }
  }
}

void HTreeBuilder::checkLegalityAndCost(const Point<double>& oldLoc,
                                        const Point<double>& newLoc,
                                        const Point<double>& parentPoint,
                                        double targetDist,
                                        const std::vector<Point<double>>& sinks,
                                        int scalingFactor,
                                        Point<double>& bestLoc,
                                        double& sinkDist,
                                        double& bestSinkDist)
{
  if (fuzzyEqual(computeDist(newLoc, parentPoint), targetDist)
      && checkLegalityLoc(newLoc, scalingFactor)) {
    sinkDist = weightedDistance(newLoc, oldLoc, sinks);
    if (sinkDist < bestSinkDist) {
      bestLoc = newLoc;
      bestSinkDist = sinkDist;
    }
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBeyondBlockage: branchPt "
	       "move:{}=>{} is legal, dist={:0.3f}, sinkDist={:0.3f}",
	       oldLoc, newLoc, targetDist, sinkDist);
  } else {
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBeyondBlockage: branchPt "
	       "move:{}=>{} is illegal or dist {:0.3f} != {:0.3f}",
	       oldLoc, newLoc, computeDist(newLoc, parentPoint), targetDist);
    // clang-format on
  }
}

void HTreeBuilder::legalizeDummy()
{
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
      Point<double>& branchPoint = topology.getBranchingPoint(idx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? topLevelBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(idx);

      double leng = topology.getLength();
      Point<double>& sibLoc = findSiblingLoc(topology, idx, parentIdx);

      double d1 = computeDist(branchPoint, sibLoc);
      double d2 = computeDist(branchPoint, parentPoint);
      bool overlap = d1 == 0 || d2 == 0;
      bool dummy = sinks.empty();  // dummy buffers drive no sinks

      // not important, can be removed later?
      if (dummy) {
        setSiblingPosition(sibLoc, branchPoint, parentPoint);
        scalePosition(branchPoint, parentPoint, leng, 0.1);
      } else if (overlap) {
        scalePosition(branchPoint, parentPoint, leng, 0.1);
      } else {
        continue;
      }

      double x1, y1, x2, y2;
      int scalingFactor = wireSegmentUnit_;
      if (!isOccupiedLoc(branchPoint)
          && findBlockage(branchPoint, scalingFactor, x1, y1, x2, y2)) {
        Point<double> legalBranchPoint(branchPoint);
        std::vector<Point<double>> legalLocations;
        findLegalLocations(
            parentPoint, branchPoint, x1, y1, x2, y2, legalLocations);
        legalBranchPoint = findBestLegalLocation(topology.getLength(),
                                                 branchPoint,
                                                 parentPoint,
                                                 legalLocations,
                                                 sinks,
                                                 x1,
                                                 y2,
                                                 x2,
                                                 y2,
                                                 scalingFactor,
                                                 odb::Direction2D::North);
        double d = computeDist(legalBranchPoint, parentPoint);
        // clang-format off
        debugPrint(logger_, CTS, "legalizer", 1,
            "legalizeDummy level index {}: {}->{} d={:0.3f}, leng={:0.3f},"
		   "ratio={:0.3f}", levelIdx, branchPoint, legalBranchPoint,
		   d, leng, d / leng);
        // clang-format on
        commitMoveLoc(branchPoint, legalBranchPoint);
        branchPoint.setX(legalBranchPoint.getX());
        branchPoint.setY(legalBranchPoint.getY());
      }
    }
  }
}

void HTreeBuilder::legalize()
{
  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("HTree before legalization -------");
    printHTree();
  }
  // make sure top level buffer is legal
  Point<double> oldTopBufferLoc = sinkRegion_.getCenter();
  Point<double> newTopBufferLoc
      = legalizeOneBuffer(oldTopBufferLoc, options_->getRootBuffer());
  sinkRegion_.setCenter(newTopBufferLoc);
  commitMoveLoc(oldTopBufferLoc, newTopBufferLoc);
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "legalize: top buf loc:{}->{}",
	     oldTopBufferLoc, newTopBufferLoc);
  // clang-format on
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned bufferIdx = 0; bufferIdx < topology.getBranchingPointSize();
         ++bufferIdx) {
      // bufferIdx is the buffer id at level levelIdx
      Point<double>& branchPoint = topology.getBranchingPoint(bufferIdx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(bufferIdx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? newTopBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      odb::Direction2D::Value branch_point_dir;
      if (isHorizontal(levelIdx + 1)) {
        if (branchPoint.getX() - parentPoint.getX() > 0) {
          branch_point_dir = odb::Direction2D::East;
        } else {
          branch_point_dir = odb::Direction2D::West;
        }
      } else {
        if (branchPoint.getY() - parentPoint.getY() > 0) {
          branch_point_dir = odb::Direction2D::North;
        } else {
          branch_point_dir = odb::Direction2D::South;
        }
      }

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(bufferIdx);

      double leng = computeDist(branchPoint, parentPoint);
      // clang-format off
      if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
        logger_->report("  HTree level*{}* bufId*{}*, parent:{}, branch:{}, "
			"leng:{:0.3f}, sinks:{}", levelIdx, bufferIdx,
			parentPoint, branchPoint, leng, sinks.size());
      }
      // clang-format on
      int scalingFactor = wireSegmentUnit_;
      double x1, y1, x2, y2;
      if (!isOccupiedLoc(branchPoint)
          && findBlockage(branchPoint, scalingFactor, x1, y1, x2, y2)) {
        Point<double> legalBranchPoint(branchPoint);
        std::vector<Point<double>> legalLocations;
        // find all the possible locations off the blockage
        findLegalLocations(
            parentPoint, branchPoint, x1, y1, x2, y2, legalLocations);
        // choose the best new location based on desired topology length
        legalBranchPoint = findBestLegalLocation(topology.getLength(),
                                                 branchPoint,
                                                 parentPoint,
                                                 legalLocations,
                                                 sinks,
                                                 x1,
                                                 y1,
                                                 x2,
                                                 y2,
                                                 scalingFactor,
                                                 branch_point_dir);
        // clang-format off
	debugPrint(logger_, CTS, "legalizer", 1,
		   "findBestLegalLocation branchPt:{}=>{} parentPt:{} new "
		   "branchPt is {} blockage", branchPoint, legalBranchPoint,
		   parentPoint,
		   isInsideBbox(legalBranchPoint.getX(), legalBranchPoint.getY(),
				x1, y1, x2, y2)? "inside" : "outside");
        // clang-format on
        // update branchPoint
        commitMoveLoc(branchPoint, legalBranchPoint);
        branchPoint = legalBranchPoint;
      } else if (isOccupiedLoc(branchPoint)
                 || !fuzzyEqual(leng, topology.getLength(), 0.01)) {
        // legal branch point needs adjustment if parent point moved in previous
        // level
        Point<double> newLocation(branchPoint);
        newLocation = adjustBeyondBlockage(branchPoint,
                                           parentPoint,
                                           topology.getLength(),
                                           sinks,
                                           scalingFactor,
                                           branch_point_dir);
        // clang-format off
	debugPrint(logger_, CTS, "legalizer", 3,
		   "adjustBeyondBlockage applied to legal branchPt:"
		   "{}=>{} parentPt:{} newDist={:0.3f}", branchPoint, newLocation,
		   parentPoint, computeDist(newLocation, parentPoint));
        // clang-format on
        commitMoveLoc(branchPoint, newLocation);
        branchPoint = newLocation;
      } else {
        commitLoc(branchPoint);
      }
    }
  }

  // "further" optimize the location of the "dummy" buffers that drive
  // no sinks (still needed?)
  legalizeDummy();

  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("HTree after legalization -------");
    printHTree();
  }
}

void HTreeBuilder::run()
{
  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  bool useMaxCap = (type_ == TreeType::MacroTree)
                       ? false
                       : options_->getSinkClusteringUseMaxCap();

  logger_->info(
      CTS, 27, "Generating H-Tree topology for net {}.", clock_.getName());
  logger_->info(CTS, 28, " Total number of sinks: {}.", clock_.getNumSinks());
  if (options_->getSinkClustering()) {
    if (useMaxCap) {
      logger_->info(
          CTS, 90, " Sinks will be clustered based on buffer max cap.");
    } else {
      logger_->info(
          CTS,
          29,
          " {} sinks will be clustered in groups of up to {} and with "
          "maximum cluster diameter of {:.1f} um.",
          type_ == TreeType::MacroTree ? "Macro " : "Register",
          clusterSize,
          clusterDiameter);
    }
  }
  logger_->info(
      CTS, 30, " Number of static layers: {}.", options_->getNumStaticLayers());

  clockTreeMaxDepth_ = options_->getClockTreeMaxDepth();
  minInputCap_ = techChar_->getActualMinInputCap();
  numMaxLeafSinks_ = options_->getNumMaxLeafSinks();
  // Autonomous iter2: CTS_NUM_MAX_LEAF_SINKS overrides the H-tree stop criterion
  // (larger => shallower backbone, fewer/bigger leaf regions). Distinct from
  // leaf sink-clustering granularity. Default (env unset) = unchanged. Used to
  // probe whether backbone depth helps the 2 losers (ibex-gated, ariane-macro)
  // without touching winners (they don't set the env).
  if (const char* e = std::getenv("CTS_NUM_MAX_LEAF_SINKS")) {
    int v = std::atoi(e);
    if (v > 0) {
      numMaxLeafSinks_ = (unsigned) v;
      logger_->info(CTS, 407,
                    "CTS_NUM_MAX_LEAF_SINKS override: numMaxLeafSinks_={}", v);
    }
  }
  minLengthSinkRegion_ = techChar_->getMinSegmentLength() * 2;

  // skewTargetBeta_ / CTS_SKEW_TARGET_BETA removed.
  // Was dead code — computeDist() called with candidate branching points
  // not in mapLocationToSink_ → beta penalty never fired.
  // Target-aware clustering: CTS_CKMEANS_TARGET_BETA (read at line ~1725).
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    logger_->info(CTS, 400,
                  "{} skew targets loaded",
                  cts3dDb_->getSkewTargetCount());
  }

  // Initialize useful-skew wire adjustment parameters
  // CTS_WIRE_SKEW_SCALE: 0 = disabled (zero-skew), 1.0 = full LP-guided wire shift
  // CTS_WIRE_DELAY_PS_UM: wire delay per unit length (ps/um)
  //   Used to convert LP target difference (ns) → wire length shift (um)
  //   Reference: Fishburn 1990, LP-SAFETY clock skew optimization
  wireSkewScale_    = 0.0;
  wireDelayPerUnit_ = 0.001;  // default: 1.0 ps/um in ns/um
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    const char* wss_env  = std::getenv("CTS_WIRE_SKEW_SCALE");
    const char* wdpu_env = std::getenv("CTS_WIRE_DELAY_PS_UM");
    if (wss_env != nullptr) {
      wireSkewScale_ = std::atof(wss_env);
    } else {
      wireSkewScale_ = 0.0;  // V58 Fix D: default OFF — hold regression (V51b confirmed)
    }
    if (wdpu_env != nullptr) {
      wireDelayPerUnit_ = std::atof(wdpu_env) * 0.001;  // ps/um → ns/um
    }
    if (wireSkewScale_ > 0.0) {
      logger_->info(CTS, 405,
                    "Useful-skew wire adjustment enabled: scale={:.2f}, "
                    "wireDelay={:.3f}ps/um",
                    wireSkewScale_, wireDelayPerUnit_ * 1000.0);
    }
  }

  // Read target-split env vars BEFORE initSinkRegion(),
  // because preSinkClustering() (called inside initSinkRegion()) uses
  // enableTargetSplit_ at line 328.  Previously read at line 1576 — too late.
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    if (const char* e = std::getenv("CTS_ENABLE_TARGET_SPLIT"))
      enableTargetSplit_ = (std::atoi(e) != 0);
    if (const char* e = std::getenv("CTS_SPLIT_THRESHOLD_2WAY_PS"))
      splitThreshold2wayNs_ = std::atof(e) * 0.001;  // ps → ns
    if (const char* e = std::getenv("CTS_SPLIT_THRESHOLD_3WAY_PS"))
      splitThreshold3wayNs_ = std::atof(e) * 0.001;  // ps → ns
  }

  // Tier-aware cluster split
  if (const char* e = std::getenv("CTS_ENABLE_TIER_SPLIT"))
    enableTierSplit_ = (std::atoi(e) != 0);
  logger_->info(CTS, 376, "Tier split: {}", enableTierSplit_ ? "ON" : "OFF");

  // TARP clustering
  if (const char* e = std::getenv("CTS_ENABLE_TARP"))
    enableTarp_ = (std::atoi(e) != 0);
  if (const char* e = std::getenv("CTS_TARP_ALPHA"))
    tarpAlpha_ = std::stod(e);
  // Seeded refinement (off clustering + critical-pair moves).
  // Clean master switch: CTS_SEED_REFINE=1 enables it STANDALONE —
  // no CTS_ENABLE_TARP needed (seeded is not a TARP sub-mode; it never runs FM).
  // Back-compat: the old pair CTS_ENABLE_TARP=1 + CTS_TARP_SEEDED=1 still works.
  if (const char* e = std::getenv("CTS_TARP_SEEDED"))
    tarpSeeded_ = (std::atoi(e) != 0) && enableTarp_;
  if (const char* e = std::getenv("CTS_SEED_REFINE"))
    tarpSeeded_ = (std::atoi(e) != 0);
  if (const char* e = std::getenv("CTS_ENABLE_PER_TIER_CTS"))
    enablePerTierCts_ = (std::atoi(e) != 0);
  // CTS_ENABLE_UNIFIED_TREE forces a single Pin3D-style H-tree: registers and
  // macros are already merged into one sink set upstream, so per-tier split
  // (shared root + bottom/upper subtrees) must be OFF. Default off → no-op.
  if (const char* e = std::getenv("CTS_ENABLE_UNIFIED_TREE")) {
    if (std::atoi(e) != 0 && enablePerTierCts_) {
      enablePerTierCts_ = false;
      logger_->info(CTS, 866,
          "CTS_ENABLE_UNIFIED_TREE: per-tier CTS forced OFF (unified tree).");
    }
  }
  {
    // V85_TARP: surface surgical mode in the banner (parsed again at tp
    // setup; read here only for logging).
    const char* co = std::getenv("CTS_TARP_CRIT_ONLY");
    const bool crit_only = (co != nullptr && std::atoi(co) != 0);
    logger_->info(CTS, 840,
                  "TARP: {} (alpha={:.2f}, critOnly={}, seeded={}), per-tier: {}",
                  enableTarp_ ? "ON" : "OFF", tarpAlpha_,
                  crit_only ? "ON" : "OFF",
                  tarpSeeded_ ? "ON" : "OFF",
                  enablePerTierCts_ ? "ON" : "OFF");
  }

  // P3 directional asymmetric branch-point embedding (structural useful-skew).
  // Default OFF → the 9 current winners are byte-identical unless opted in.
  if (const char* e = std::getenv("CTS_ENABLE_DIR_SKEW"))
    dirSkewEnable_ = (std::atoi(e) != 0);
  if (const char* e = std::getenv("CTS_DIR_SKEW_ALPHA"))
    dirSkewAlpha_ = std::atof(e);
  if (const char* e = std::getenv("CTS_DIR_MAX_SHIFT"))
    dirSkewMaxShiftUm_ = std::atof(e);
  if (const char* e = std::getenv("CTS_DIR_SKEW_MAX_LEVEL"))
    dirSkewMaxLevel_ = std::atoi(e);
  if (dirSkewEnable_) {
    computeDirSkewVector();  // sets dirSkewVecX/Y_ + dirSkewValid_
    logger_->info(CTS, 863,
        "P3 directional useful-skew: ON (alpha={:.2f}, maxShift={:.1f}um, "
        "upperLevels<={}, proxyValid={})",
        dirSkewAlpha_, dirSkewMaxShiftUm_, dirSkewMaxLevel_, dirSkewValid_);
  }

  // Per-tier independent CTS.
  // MUST be placed BEFORE initSinkRegion() to avoid double preSinkClustering().
  // Env vars for grouped delay, relay, GH-tree etc. are read inside
  // readSkewDeliveryEnvVars() called from runSingleTierTree().
  if (enablePerTierCts_ && cts3dDb_ != nullptr) {
    // Read delivery env vars here (before runSingleTierTree) so they are set
    readSkewDeliveryEnvVars();

    // P6 equalize env (default off → 9 winners unaffected unless opted in).
    if (const char* e = std::getenv("CTS_ENABLE_PER_TIER_EQUALIZE"))
      enablePerTierEqualize_ = (std::atoi(e) != 0);
    if (const char* e = std::getenv("CTS_PER_TIER_EQUALIZE_MAX_BUFS"))
      perTierEqualizeMaxBufs_ = std::atoi(e);

    // P2 dominant-tier graft env (default off → 9 winners byte-identical).
    if (const char* e = std::getenv("CTS_ENABLE_DOMTIER_GRAFT"))
      enableDomTierGraft_ = (std::atoi(e) != 0);
    // P2 and P6 both target the tier seam; keep them mutually exclusive.
    // P2 gives the majority Pin3D-equivalence P6 could not, so P2 wins.
    if (enableDomTierGraft_ && enablePerTierEqualize_) {
      logger_->warn(CTS, 865,
          "P2 dominant-tier graft and P6 per-tier equalize are mutually "
          "exclusive; using P2 (graft) and disabling P6.");
      enablePerTierEqualize_ = false;
    }

    int tier0Count = 0, tier1Count = 0;
    clock_.forEachSink([&](const ClockInst& sink) {
      odb::dbInst* dbInst = sink.getDbInst();
      if (dbInst == nullptr) return;
      int t = cts3dDb_->getInstTier(dbInst);
      if (t == 0) tier0Count++;
      else if (t == 1) tier1Count++;
    });
    logger_->info(CTS, 850, "Per-tier split: tier0={} FFs, tier1={} FFs",
                  tier0Count, tier1Count);

    // ---- P2 dominant-tier graft (default off; skips the shared root) ----
    // Each tier builds its OWN root-direct H-tree via the unified else-branch of
    // createClockSubNets() (sharedRootBuffer_ stays nullptr). The majority root
    // becomes the top clock root; the minority root is grafted onto it below.
    // Skip the small special-cased MACRO tree (Fix-G) — graft only the large
    // register tree where the -2.5% majority-indirection loss lives. Macro tree
    // keeps the shared-root path (crashes in preSinkClustering without it).
    if (enableDomTierGraft_ && tier0Count > 0 && tier1Count > 0
        && type_ != TreeType::MacroTree) {
      // createSharedRoot() (which we skip) normally sets wireSegmentUnit_;
      // set it here or preSinkClustering divides by 0 → SIGSEGV.
      wireSegmentUnit_ = techChar_->getLengthUnit();
      const int majTier = (tier1Count >= tier0Count) ? 1 : 0;
      const int minTier = 1 - majTier;
      graftMajorityTier_ = majTier;               // MUST be set before the builds
      runSingleTierTree(majTier);                 // majority first → sets topBufferName_
      runSingleTierTree(minTier);                 // minority second → root captured
      if (graftMajRootSubNet_ != nullptr && graftMinRootBuf_ != nullptr) {
        graftMajRootSubNet_->addInst(*graftMinRootBuf_);
        logger_->info(CTS, 859,
            "P2 dominant-tier graft: majority tier {} root-direct (top root '{}'), "
            "minority tier {} root '{}' grafted onto majority root subnet '{}'",
            majTier, topBufferName_, minTier,
            graftMinRootBuf_->getName(), graftMajRootSubNet_->getName());
      } else {
        logger_->warn(CTS, 862,
            "P2 dominant-tier graft: capture failed "
            "(majSubNet set={}, minRootBuf set={}) — trees may be disconnected",
            graftMajRootSubNet_ != nullptr, graftMinRootBuf_ != nullptr);
      }
      return;  // Skip shared-root path entirely
    }

    // Create ONE shared root buffer for both tiers
    createSharedRoot();

    // Build per-tier sub-trees (L1+ branches attach to shared root).
    // P6 equalize: build the MAJORITY (deeper) tier first so its level count is
    // known, then pad the MINORITY tier's trunk to match. Ordering swap ONLY
    // under the env gate → default path is byte-identical.
    if (enablePerTierEqualize_ && tier0Count > 0 && tier1Count > 0) {
      const int majTier = (tier1Count >= tier0Count) ? 1 : 0;
      const int minTier = 1 - majTier;
      equalizeActiveTier_ = -1;              // majority: never padded
      runSingleTierTree(majTier);
      equalizeTargetLevels_ = lastTierLevels_;
      equalizeActiveTier_ = minTier;         // minority: pad to majority depth
      runSingleTierTree(minTier);
    } else {
      if (tier0Count > 0) runSingleTierTree(0);
      if (tier1Count > 0) runSingleTierTree(1);
    }
    return;  // Skip unified tree path
  }

  auto prof_t0 = std::chrono::high_resolution_clock::now();
  initSinkRegion();
  auto prof_t1 = std::chrono::high_resolution_clock::now();
  double prof_clustering_s = std::chrono::duration<double>(prof_t1 - prof_t0).count();
  logger_->info(CTS, 911, "Profile: initSinkRegion (TARP+clustering) {:.1f}s", prof_clustering_s);

  // Fix G: Macro flat bypass (same logic as runSingleTierTree path)
  {
    int macroFlatMax = 50;
    if (const char* e = std::getenv("CTS_MACRO_FLAT_MAX_SINKS")) {
      macroFlatMax = std::atoi(e);
    }
    const int nSinks = clock_.getNumSinks();
    if (type_ == TreeType::MacroTree && macroFlatMax > 0 && nSinks <= macroFlatMax) {
      bool allBlock = true;
      clock_.forEachSink([&](const ClockInst& sink) {
        odb::dbInst* dbInst = sink.getDbInst();
        if (!isMacroBlockInst(dbInst)) {
          allBlock = false;
        }
      });
      if (allBlock) {
        logger_->info(CTS, 695,
            "Fix G: Macro flat bypass — {} sinks (all block), threshold={}, "
            "skipping H-tree → flat topology",
            nSinks, macroFlatMax);
        createSingleBufferClockNet();
        treeBufLevels_++;
        return;
      }
    }
  }

  // V58 Fix E: Macro tree depth clamp (fallback for non-flat macro trees)
  int effectiveMaxDepth = clockTreeMaxDepth_;
  if (type_ == TreeType::MacroTree) {
    int macroMaxDepth = 2;
    if (const char* e = std::getenv("CTS_MACRO_TREE_MAX_DEPTH")) {
      macroMaxDepth = std::atoi(e);
    }
    if (macroMaxDepth > 0 && macroMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 692,
          "V58 Fix E: Macro tree depth clamped {} -> {} (sinks={})",
          effectiveMaxDepth, macroMaxDepth, clock_.getNumSinks());
      effectiveMaxDepth = macroMaxDepth;
    }
  }

  auto prof_htree_t0 = std::chrono::high_resolution_clock::now();
  for (int level = 1; level <= effectiveMaxDepth; ++level) {
    const unsigned numSinksPerSubRegion
        = computeNumberOfSinksPerSubRegion(level);
    double regionWidth, regionHeight;
    computeSubRegionSize(level, regionWidth, regionHeight);

    if (isSubRegionTooSmall(regionWidth, regionHeight)) {
      if (options_->isFakeLutEntriesEnabled()) {
        const unsigned minIndex = 1;
        techChar_->createFakeEntries(minLengthSinkRegion_, minIndex);
        minLengthSinkRegion_ = 1;
      } else {
        logger_->info(
            CTS,
            31,
            " Stop criterion found. Min length of sink region is ({}).",
            minLengthSinkRegion_);
        break;
      }
    }

    computeLevelTopology(level, regionWidth, regionHeight);

    if (isNumberOfSinksTooSmall(numSinksPerSubRegion)) {
      logger_->info(CTS,
                    32,
                    " Stop criterion found. Max number of sinks is {}.",
                    options_->getMaxFanout() ? options_->getMaxFanout()
                                             : numMaxLeafSinks_);
      break;
    }
  }

  auto prof_htree_t1 = std::chrono::high_resolution_clock::now();
  double prof_htree_s = std::chrono::duration<double>(prof_htree_t1 - prof_htree_t0).count();
  logger_->info(CTS, 912, "Profile: H-tree level loop ({} levels) {:.1f}s",
                static_cast<int>(topologyForEachLevel_.size()), prof_htree_s);

  if (topologyForEachLevel_.empty()) {
    createSingleBufferClockNet();
    treeBufLevels_++;
    return;
  }

  clock_.setMaxLevel(topologyForEachLevel_.size());

  if (options_->getPlotSolution()
      || logger_->debugCheck(utl::CTS, "HTree", 2)) {
    plotSolution();
  }

  if (CtsObserver* observer = options_->getObserver()) {
    observer->initializeWithClock(this, clock_);
  }

  if (options_->getObstructionAware()) {
    legalize();
  }

  // Extracted to readSkewDeliveryEnvVars().
  // Called from per-tier block (above) or here for unified mode.
  readSkewDeliveryEnvVars();

  auto prof_subnets_t0 = std::chrono::high_resolution_clock::now();
  createClockSubNets();
  auto prof_subnets_t1 = std::chrono::high_resolution_clock::now();
  double prof_subnets_s = std::chrono::duration<double>(prof_subnets_t1 - prof_subnets_t0).count();
  logger_->info(CTS, 913, "Profile: createClockSubNets (leaf+TAP+delay) {:.1f}s", prof_subnets_s);

  // Total run() time
  auto prof_total_t1 = std::chrono::high_resolution_clock::now();
  double prof_total_s = std::chrono::duration<double>(prof_total_t1 - prof_t0).count();
  logger_->info(CTS, 914,
      "Profile: run() total {:.1f}s (clustering={:.1f} htree={:.1f} subnets={:.1f})",
      prof_total_s, prof_clustering_s, prof_htree_s, prof_subnets_s);

  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "Htree file {} has been generated",
             plotHTree());
  debugPrint(logger_, CTS, "legalizer", 3, "Run 'obsAwareCts.py cts.clk.buffer'"
	     "to produce cts.clk.buffer.png");
  // clang-format on
}

// Extracted from run() — reads all skew delivery env vars.
// Safe to call multiple times (idempotent).
void HTreeBuilder::readSkewDeliveryEnvVars()
{
  if (cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) return;

  const char* thr_env = std::getenv("CTS_DELAY_TARGET_THRESHOLD");
  delayTargetThreshold_ = (thr_env != nullptr) ? std::atof(thr_env) : 0.010;

    // Single delay-buffer delay (ns) and max chain length per leaf
    if (const char* e = std::getenv("CTS_LEAF_BUF_DELAY_NS"))
      singleBufDelay_ = std::atof(e);
    if (const char* e = std::getenv("CTS_MAX_LEAF_DELAY_BUFS"))
      maxLeafDelayBufs_ = std::atoi(e);

    // Target-aware sub-clustering parameters
    // NOTE: env var reading moved to before initSinkRegion() (V42-fix).
    // Only logging remains here.
    if (enableTargetSplit_) {
      logger_->info(CTS, 410,
                    "Target-aware sub-clustering enabled: "
                    "2way_threshold={:.1f}ps, 3way_threshold={:.1f}ps",
                    splitThreshold2wayNs_ * 1000.0,
                    splitThreshold3wayNs_ * 1000.0);
    }

    // Per-FF relay buffer parameters
    enablePerFfRelay_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_PER_FF_RELAY"))
      enablePerFfRelay_ = (std::atoi(e) != 0);
    if (enablePerFfRelay_ && cts3dDb_ && cts3dDb_->hasSkewTargets()) {
      if (const char* e = std::getenv("CTS_PER_FF_BUF_DELAY_NS"))
        perFfBufDelay_ = std::atof(e);
      if (const char* e = std::getenv("CTS_PER_FF_MAX_RELAY"))
        perFfMaxRelay_ = std::atoi(e);
      if (const char* e = std::getenv("CTS_PER_FF_HOLD_GUARD_NS"))
        perFfHoldGuard_ = std::atof(e);
      if (const char* e = std::getenv("CTS_PER_FF_MIN_TARGET_NS"))
        perFfMinTarget_ = std::atof(e);
      {
        auto hb_t0 = std::chrono::high_resolution_clock::now();
        if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
          loadPerFfHoldBudgets(csv);
        if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
          loadIoHoldBudgets(ioCsv);
        auto hb_t1 = std::chrono::high_resolution_clock::now();
        double hb_s = std::chrono::duration<double>(hb_t1 - hb_t0).count();
        logger_->info(CTS, 917, "Profile: loadHoldBudgets (CSV parse) {:.1f}s", hb_s);
      }
      // Elmore wire RC parameters for x_useful positioning
      relayRwKOhmPerUm_ = 0.0;
      relayCwFfPerUm_   = 0.0;
      relayDbuPerUm_    = 2000.0;
      if (const char* e = std::getenv("CTS_RELAY_RW_PER_UM"))
        relayRwKOhmPerUm_ = std::atof(e);
      if (const char* e = std::getenv("CTS_RELAY_CW_PER_UM"))
        relayCwFfPerUm_ = std::atof(e);
      if (const char* e = std::getenv("CTS_DBU_PER_UM"))
        relayDbuPerUm_ = std::atof(e);
      const bool useElmore = (relayRwKOhmPerUm_ > 1e-12 && relayCwFfPerUm_ > 1e-12);
      logger_->info(CTS, 420,
                    "Per-FF relay: bufDelay={:.3f}ns maxRelay={} "
                    "holdGuard={:.3f}ns minTarget={:.3f}ns holdBudgets={} "
                    "Elmore={}(rw={:.4f}kOhm/um cw={:.4f}fF/um)",
                    perFfBufDelay_, perFfMaxRelay_, perFfHoldGuard_,
                    perFfMinTarget_,
                    static_cast<int>(perFfHoldBudget_.size()),
                    useElmore ? "ON" : "OFF",
                    relayRwKOhmPerUm_, relayCwFfPerUm_);
    } else if (enablePerFfRelay_) {
      logger_->warn(CTS, 421,
                    "CTS_ENABLE_PER_FF_RELAY=1 but no skew targets; "
                    "falling back to per-branch mode");
      enablePerFfRelay_ = false;
    }

    // Grouped Delay Chain parameters.
    // Groups FFs within each leaf cluster by quantized LP target delta.
    // Shared delay chain per group → fewer nets than per-FF relay → avoids GRT-0183.
    enableGroupedDelay_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_GROUPED_DELAY"))
      enableGroupedDelay_ = (std::atoi(e) != 0);
    if (enableGroupedDelay_ && cts3dDb_ && cts3dDb_->hasSkewTargets()) {
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_MAX_DEPTH"))
        groupedDelayMaxDepth_ = std::atoi(e);
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_MIN_DELTA_NS"))
        groupedDelayMinDelta_ = std::atof(e);
      // Reuse per-FF relay parameters: perFfBufDelay_, perFfHoldGuard_, perFfHoldBudget_
      // Load hold budgets if not already loaded by per-FF relay block
      if (perFfHoldBudget_.empty()) {
        if (const char* e = std::getenv("CTS_PER_FF_BUF_DELAY_NS"))
          perFfBufDelay_ = std::atof(e);
        if (const char* e = std::getenv("CTS_PER_FF_HOLD_GUARD_NS"))
          perFfHoldGuard_ = std::atof(e);
        if (const char* e = std::getenv("CTS_PER_FF_MIN_TARGET_NS"))
          perFfMinTarget_ = std::atof(e);
        {
          auto hb2_t0 = std::chrono::high_resolution_clock::now();
          if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
            loadPerFfHoldBudgets(csv);
          if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
            loadIoHoldBudgets(ioCsv);
          auto hb2_t1 = std::chrono::high_resolution_clock::now();
          double hb2_s = std::chrono::duration<double>(hb2_t1 - hb2_t0).count();
          logger_->info(CTS, 918, "Profile: loadHoldBudgets (per-tier path) {:.1f}s", hb2_s);
        }
      }
      // Cluster-Uniform Depth mode — one chain per cluster, depth = cluster median.
      enableClusterUniform_ = false;
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_CLUSTER_UNIFORM"))
        enableClusterUniform_ = (std::atoi(e) != 0);

      // Liberty-based buffer delay override (from V51-BUF).
      // charBuf_ driving its own input cap = cascaded tap scenario.
      // Overrides env var CTS_PER_FF_BUF_DELAY_NS with Liberty-measured value.
      {
        const double libertyDelay = techChar_->getCharBufDelay();
        const double libInputCap  = techChar_->getCharBufInputCap();
        if (libertyDelay > 0) {
          perFfBufDelay_ = libertyDelay;
          logger_->info(CTS, 860,
              "Liberty buf delay: {:.4f}ns (charBuf input_cap={:.4f}pF)",
              perFfBufDelay_, libInputCap);
        } else {
          logger_->warn(CTS, 861,
              "Liberty delay not available; using default perFfBufDelay={:.3f}ns",
              perFfBufDelay_);
        }
      }

      // GH-Tree Cascaded Delay Tap.
      // When enabled, buildGroupChains uses cascaded tap structure instead of
      // parallel chains. Same depth delivery, ~50% fewer buffers and nets.
      enableGhTree_ = false;
      if (const char* e = std::getenv("CTS_ENABLE_GH_TREE"))
        enableGhTree_ = (std::atoi(e) != 0);
      // Adaptive depth cost ratio (default 1.0 = 1 FF per empty tap to justify)
      if (const char* e = std::getenv("CTS_ADAPTIVE_DEPTH_COST_RATIO"))
        adaptiveDepthCostRatio_ = std::atof(e);

      // Log after all env vars are read so values are accurate
      logger_->info(CTS, 490,
                    "Grouped delay: maxDepth={} minDelta={:.3f}ns "
                    "bufDelay={:.3f}ns holdBudgets={} "
                    "clusterUniform={} cascadedTap={} adaptiveCost={:.1f}",
                    groupedDelayMaxDepth_, groupedDelayMinDelta_,
                    perFfBufDelay_,
                    static_cast<int>(perFfHoldBudget_.size()),
                    enableClusterUniform_ ? 1 : 0,
                    enableGhTree_ ? 1 : 0,
                    adaptiveDepthCostRatio_);
      if (enableGhTree_) {
        logger_->info(CTS, 570,
            "Cascaded tap + adaptive depth: costRatio={:.1f}, bufDelay={:.3f}ns",
            adaptiveDepthCostRatio_, perFfBufDelay_);
      }
    } else if (enableGroupedDelay_) {
      logger_->warn(CTS, 491,
                    "CTS_ENABLE_GROUPED_DELAY=1 but no skew targets; "
                    "falling back to balanced CTS");
      enableGroupedDelay_ = false;
    }

    // Intermediate delay buffer gate (V37 default on, V39 default off)
    enableMidDelayBufs_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_MID_DELAY_BUFS"))
      enableMidDelayBufs_ = (std::atoi(e) != 0);

    // CKMeans target-aware branching weight.
    // When > 0, H-tree 2-way split (CKMeans) groups FFs with similar targets.
    ckmeansTargetBeta_ = 0.0;
    if (const char* e = std::getenv("CTS_CKMEANS_TARGET_BETA"))
      ckmeansTargetBeta_ = std::atof(e);
    if (ckmeansTargetBeta_ > 0.0) {
      logger_->info(CTS, 428,
          "CKMeans target-aware branching: beta={:.1f}",
          ckmeansTargetBeta_);
    }

    computeBranchDelayTargets();
}

bool HTreeBuilder::isNumberOfSinksTooSmall(unsigned numSinksPerSubRegion) const
{
  if (options_->getMaxFanout()) {
    return numSinksPerSubRegion < options_->getMaxFanout();
  }
  return numSinksPerSubRegion < numMaxLeafSinks_;
}

std::string HTreeBuilder::plotHTree()
{
  auto name = std::string("cts.") + clock_.getName() + ".buffer";
  std::ofstream file(name);

  plotBlockage(file, db_, wireSegmentUnit_);

  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();

  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    // clang-format off
    topology.forEachBranchingPoint(
        [&](unsigned idx, Point<double> branchPoint) {
          unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

          Point<double> parentPoint
              = (levelIdx == 0)
                    ? topLevelBufferLoc
                    : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                          parentIdx);

          const std::vector<Point<double>>& sinks
              = topology.getBranchSinksLocations(idx);

          plotSinks(file, sinks);

          double x1 = parentPoint.getX();
          double y1 = parentPoint.getY();
          double x2 = branchPoint.getX();
          double y2 = branchPoint.getY();
          std::string name = "buffer";
          file << levelIdx << " " << x1 << " " << y1 << " " << x2 << " " << y2;
          file << " " << name << '\n';
        });
    // clang-format on
  }

  LevelTopology& leafTopology = topologyForEachLevel_.back();
  unsigned numSinks = 0;
  leafTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        double px = branchPoint.getX();
        double py = branchPoint.getY();

        const std::vector<Point<double>>& sinkLocs
            = leafTopology.getBranchSinksLocations(idx);

        for (const Point<double>& loc : sinkLocs) {
          auto name2 = mapLocationToSink_[loc]->getName();

          file << numSinks << " " << loc.getX() << " " << loc.getY();
          file << " " << px << " " << py << " leafbuffer " << name2;
          file << " z=" << wireSegmentUnit_ << '\n';
          ++numSinks;
        }
      });
  file.close();
  return name;
}

unsigned HTreeBuilder::computeNumberOfSinksPerSubRegion(
    const unsigned level) const
{
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;
  unsigned totalNumSinks = 0;
  if (clock_.getNumSinks() > min_clustering_sinks
      && options_->getSinkClustering()) {
    totalNumSinks = topLevelSinksClustered_.size();
  } else {
    totalNumSinks = clock_.getNumSinks();
  }
  const unsigned numRoots = std::pow(2, level);
  const double numSinksPerRoot = (double) totalNumSinks / numRoots;
  return std::ceil(numSinksPerRoot);
}

void HTreeBuilder::computeSubRegionSize(const unsigned level,
                                        double& width,
                                        double& height) const
{
  unsigned gridSizeX = 0;
  unsigned gridSizeY = 0;
  if (isVertical(1)) {
    gridSizeY = computeGridSizeX(level);
    gridSizeX = computeGridSizeY(level);
  } else {
    gridSizeX = computeGridSizeX(level);
    gridSizeY = computeGridSizeY(level);
  }
  width = sinkRegion_.getWidth() / gridSizeX;
  height = sinkRegion_.getHeight() / gridSizeY;
}

void HTreeBuilder::computeLevelTopology(const unsigned level,
                                        const double width,
                                        const double height)
{
  const unsigned numSinksPerSubRegion = computeNumberOfSinksPerSubRegion(level);
  logger_->report(" Level {}", level);
  logger_->report("    Direction: {}",
                  (isVertical(level)) ? ("Vertical") : ("Horizontal"));
  logger_->report("    Sinks per sub-region: {}", numSinksPerSubRegion);
  logger_->report("    Sub-region size: {:.4f} X {:.4f}", width, height);

  const unsigned minLength = minLengthSinkRegion_;
  const unsigned clampedMinLength = std::max(minLength, 1u);

  unsigned segmentLength
      = std::round(width / (double) clampedMinLength) * minLength / 2;

  if (isVertical(level)) {
    segmentLength
        = std::round(height / (double) clampedMinLength) * minLength / 2;
  }
  segmentLength = std::max<unsigned>(segmentLength, 1);

  LevelTopology topology(segmentLength);

  logger_->info(CTS, 34, "    Segment length (rounded): {}.", segmentLength);

  const int vertexBufferLength
      = options_->getVertexBufferDistance() / (techChar_->getLengthUnit() * 2);
  int remainingLength
      = options_->getBufferDistance() / (techChar_->getLengthUnit());
  int currWl = 0;
  unsigned inputCap = minInputCap_;
  unsigned inputSlew = 1;
  if (level > 1) {
    const LevelTopology& previousLevel = topologyForEachLevel_[level - 2];
    inputCap = previousLevel.getOutputCap();
    inputSlew = previousLevel.getOutputSlew();
    remainingLength = previousLevel.getRemainingLength();
    currWl = previousLevel.getCurrWl();
  }

  const unsigned kSlewThreshold = options_->getMaxSlew();
  const unsigned kInitTolerance = 1;

  int wirelengthThreshold;
  // If max wirelength is 0, set it as slew threshold  * maximum topology
  // wirelength. This will behave as if there was no max wirelength threshold.
  if (!options_->getMaxWl()) {
    wirelengthThreshold = kSlewThreshold * techChar_->getMaxSegmentLength();
  } else {
    wirelengthThreshold = options_->getMaxWl() / options_->getWireSegmentUnit();
  }

  debugPrint(
      logger_, CTS, "tech char", 1, "slew threshold = {}", kSlewThreshold);
  debugPrint(logger_,
             CTS,
             "tech char",
             1,
             "wirelength threshold = {}",
             wirelengthThreshold);
  unsigned length = 0;
  for (int charSegLength = techChar_->getMaxSegmentLength(); charSegLength >= 1;
       --charSegLength) {
    const unsigned numWires = (segmentLength - length) / charSegLength;

    if (numWires >= 1) {
      for (int wireCount = 0; wireCount < numWires; ++wireCount) {
        debugPrint(logger_, CTS, "tech char", 1, "curr wl = {}", currWl);
        unsigned outCap = 0, outSlew = 0;
        unsigned key = 0;
        if (options_->isSimpleSegmentEnabled()) {
          remainingLength -= charSegLength;

          if (segmentLength >= vertexBufferLength && (wireCount + 1 >= numWires)
              && options_->isVertexBuffersEnabled()) {
            remainingLength = 0;
            key = computeMinDelaySegment(charSegLength,
                                         inputSlew,
                                         inputCap,
                                         kSlewThreshold,
                                         kInitTolerance,
                                         outSlew,
                                         outCap,
                                         true,
                                         remainingLength);
            remainingLength
                += options_->getBufferDistance() / (techChar_->getLengthUnit());
          } else {
            if (remainingLength <= 0) {
              key = computeMinDelaySegment(charSegLength,
                                           inputSlew,
                                           inputCap,
                                           kSlewThreshold,
                                           kInitTolerance,
                                           outSlew,
                                           outCap,
                                           true,
                                           remainingLength);
              remainingLength += options_->getBufferDistance()
                                 / (techChar_->getLengthUnit());
            } else {
              key = computeMinDelaySegment(charSegLength,
                                           inputSlew,
                                           inputCap,
                                           kSlewThreshold,
                                           kInitTolerance,
                                           outSlew,
                                           outCap,
                                           false,
                                           remainingLength);
            }
          }
        } else {
          key = computeMinDelaySegment(charSegLength,
                                       inputSlew,
                                       inputCap,
                                       kSlewThreshold,
                                       wirelengthThreshold,
                                       kInitTolerance,
                                       outSlew,
                                       outCap,
                                       currWl);
        }

        if (key == std::numeric_limits<unsigned>::max()) {
          // No tech char entry found.
          continue;
        }

        length += charSegLength;
        techChar_->reportSegment(key);

        inputCap = std::max(outCap, minInputCap_);
        inputSlew = outSlew;
        topology.addWireSegment(key);
        topology.setRemainingLength(remainingLength);
      }

      if (length == segmentLength) {
        break;
      }
    }
  }

  topology.setOutputSlew(inputSlew);
  topology.setOutputCap(inputCap);
  topology.setCurrWl(currWl);

  computeBranchingPoints(level, topology);
  topologyForEachLevel_.push_back(topology);
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();

  techChar_->forEachWireSegment(
      length, 1, 1, [&](unsigned key, const WireSegment& seg) {
        if (!seg.isBuffered()) {
          return;
        }
        if (seg.getDelay() < minDelay) {
          minKey = key;
          minDelay = seg.getDelay();
        }
      });

  return minKey;
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length,
                                              const unsigned inputSlew,
                                              const unsigned inputCap,
                                              const unsigned slewThreshold,
                                              const int wirelengthThreshold,
                                              const unsigned tolerance,
                                              unsigned& outputSlew,
                                              unsigned& outputCap,
                                              int& currWl) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKey = std::numeric_limits<unsigned>::max();
  unsigned minBufDelay = std::numeric_limits<unsigned>::max();

  for (int load = 1; load <= techChar_->getMaxCapacitance(); ++load) {
    for (int outSlew = 1; outSlew <= techChar_->getMaxSlew(); ++outSlew) {
      techChar_->forEachWireSegment(
          length, load, outSlew, [&](unsigned key, const WireSegment& seg) {
            if (std::abs((int) seg.getInputCap() - (int) inputCap) > tolerance
                || std::abs((int) seg.getInputSlew() - (int) inputSlew)
                       > tolerance) {
              return;
            }

            if (seg.isBuffered()
                && (currWl + seg.getWl2FirstBuffer() > wirelengthThreshold)) {
              return;
            }

            if (!seg.isBuffered()
                && (currWl + seg.getWl2FirstBuffer()
                    > wirelengthThreshold - techChar_->getMinSegmentLength())) {
              return;
            }

            if (seg.getDelay() < minDelay) {
              minDelay = seg.getDelay();
              minKey = key;
            }

            if (seg.isBuffered() && seg.getDelay() < minBufDelay) {
              minBufDelay = seg.getDelay();
              minBufKey = key;
            }
          });
    }
  }

  const unsigned MAX_TOLERANCE = 10;
  if (inputSlew >= slewThreshold) {
    if (minBufKey < std::numeric_limits<unsigned>::max()) {
      const WireSegment& bestBufSegment = techChar_->getWireSegment(minBufKey);
      outputSlew = bestBufSegment.getOutputSlew();
      outputCap = bestBufSegment.getLoad();
      currWl = bestBufSegment.getLastWl();
      return minBufKey;
    }
    if (tolerance < MAX_TOLERANCE) {
      // Increasing tolerance
      return computeMinDelaySegment(length,
                                    inputSlew,
                                    inputCap,
                                    slewThreshold,
                                    wirelengthThreshold,
                                    tolerance + 1,
                                    outputSlew,
                                    outputCap,
                                    currWl);
    }
  }

  if (minKey == std::numeric_limits<unsigned>::max()) {
    if (tolerance >= MAX_TOLERANCE) {
      return minKey;
    }
    // Increasing tolerance
    return computeMinDelaySegment(length,
                                  inputSlew,
                                  inputCap,
                                  slewThreshold,
                                  wirelengthThreshold,
                                  tolerance + 1,
                                  outputSlew,
                                  outputCap,
                                  currWl);
  }

  const WireSegment& bestSegment = techChar_->getWireSegment(minKey);
  if (bestSegment.isBuffered()) {
    outputSlew = bestSegment.getOutputSlew();
    currWl = bestSegment.getLastWl();
  } else {
    outputSlew
        = std::max((unsigned) bestSegment.getOutputSlew(), inputSlew + 1);
    currWl += bestSegment.getLastWl();
  }
  outputCap = bestSegment.getLoad();

  return minKey;
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length,
                                              const unsigned inputSlew,
                                              const unsigned inputCap,
                                              const unsigned slewThreshold,
                                              const unsigned tolerance,
                                              unsigned& outputSlew,
                                              unsigned& outputCap,
                                              const bool forceBuffer,
                                              const int expectedLength) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKey = std::numeric_limits<unsigned>::max();
  unsigned minBufDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKeyFallback = std::numeric_limits<unsigned>::max();
  unsigned minDelayFallback = std::numeric_limits<unsigned>::max();

  for (int load = 1; load <= techChar_->getMaxCapacitance(); ++load) {
    for (int outSlew = 1; outSlew <= techChar_->getMaxSlew(); ++outSlew) {
      techChar_->forEachWireSegment(
          length, load, outSlew, [&](unsigned key, const WireSegment& seg) {
            // Same as the other functions, however, forces a segment
            // to have a buffer in a specific location.
            const unsigned normalLength = length;
            if (!seg.isBuffered() && seg.getDelay() < minDelay) {
              minDelay = seg.getDelay();
              minKey = key;
            }
            if (seg.isBuffered() && seg.getDelay() < minBufDelay
                && seg.getNumBuffers() == 1) {
              // If buffer is in the range of 10% of the expected location, save
              // its key.
              if (seg.getBufferLocation(0)
                      > ((((double) normalLength + (double) expectedLength)
                          / (double) normalLength)
                         * 0.9)
                  && seg.getBufferLocation(0)
                         < ((((double) normalLength + (double) expectedLength)
                             / (double) normalLength)
                            * 1.1)) {
                minBufDelay = seg.getDelay();
                minBufKey = key;
              }
              if (seg.getDelay() < minDelayFallback) {
                minDelayFallback = seg.getDelay();
                minBufKeyFallback = key;
              }
            }
          });
    }
  }

  if (forceBuffer) {
    if (minBufKey != std::numeric_limits<unsigned>::max()) {
      return minBufKey;
    }

    if (minBufKeyFallback != std::numeric_limits<unsigned>::max()) {
      return minBufKeyFallback;
    }
  }

  return minKey;
}

void HTreeBuilder::computeBranchingPoints(const unsigned level,
                                          LevelTopology& topology)
{
  if (level == 1) {
    const Point<double> clockRoot(sinkRegion_.getCenter());
    Point<double> low(clockRoot);
    Point<double> high(clockRoot);
    if (isHorizontal(level)) {
      low.setX(low.getX() - topology.getLength());
      high.setX(high.getX() + topology.getLength());
    } else {
      low.setY(low.getY() - topology.getLength());
      high.setY(high.getY() + topology.getLength());
    }
    const unsigned branchPtIdx1
        = topology.addBranchingPoint(low, LevelTopology::NO_PARENT);
    const unsigned branchPtIdx2
        = topology.addBranchingPoint(high, LevelTopology::NO_PARENT);

    refineBranchingPointsWithClustering(topology,
                                        level,
                                        branchPtIdx1,
                                        branchPtIdx2,
                                        clockRoot,
                                        topLevelSinksClustered_);
    return;
  }

  LevelTopology& parentTopology = topologyForEachLevel_[level - 2];
  parentTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> clockRoot) {
        Point<double> low(clockRoot);
        Point<double> high(clockRoot);
        if (isHorizontal(level)) {
          low.setX(low.getX() - topology.getLength());
          high.setX(high.getX() + topology.getLength());
        } else {
          low.setY(low.getY() - topology.getLength());
          high.setY(high.getY() + topology.getLength());
        }
        const unsigned branchPtIdx1 = topology.addBranchingPoint(low, idx);
        const unsigned branchPtIdx2 = topology.addBranchingPoint(high, idx);

        std::vector<std::pair<float, float>> sinks;
        computeBranchSinks(parentTopology, idx, sinks);
        refineBranchingPointsWithClustering(
            topology, level, branchPtIdx1, branchPtIdx2, clockRoot, sinks);
      });
}

// Per-tier independent H-tree.
// Create a single shared root buffer at the center of ALL sinks (both tiers).
// Both per-tier sub-trees attach their L1 branches to this shared root subnet.
// This ensures tier-to-tier clock skew is controlled at the root level.
void HTreeBuilder::createSharedRoot()
{
  wireSegmentUnit_ = techChar_->getLengthUnit();

  // Compute center from ALL sinks (not filtered by tier)
  std::vector<std::pair<float, float>> allSinks;
  std::vector<const ClockInst*> allInsts;
  initTopLevelSinks(allSinks, allInsts);

  float minX = 1e30f, maxX = -1e30f, minY = 1e30f, maxY = -1e30f;
  for (auto& [x, y] : allSinks) {
    minX = std::min(minX, x); maxX = std::max(maxX, x);
    minY = std::min(minY, y); maxY = std::max(maxY, y);
  }
  Point<double> center((minX + maxX) / 2.0, (minY + maxY) / 2.0);

  // Root tier: overall dominant tier from all sinks.
  // V58_CL2: CTS_FORCE_ROOT_TIER overrides dominant-tier selection.
  // Pin3D uses bottom root buffer (same tier as majority of macros).
  // Our getDominantTierFromClock() picks upper (98% upper FFs in swerv),
  // creating a fundamentally different tree from Pin3D even with per-tier OFF.
  int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
  if (const char* e = std::getenv("CTS_FORCE_ROOT_TIER")) {
    int forced = std::atoi(e);
    if (forced >= 0 && forced != root_tier) {
      logger_->info(CTS, 700,
          "V58_CL2: Root tier forced {} -> {} (CTS_FORCE_ROOT_TIER)",
          root_tier, forced);
      root_tier = forced;
    }
  }
  sharedRootTier_ = root_tier;
  const std::string root_buffer
      = cts3dDb_->getBufferForTier(options_->getRootBuffer(), root_tier);

  Point<double> legalCenter = legalizeOneBuffer(center, root_buffer);
  const int centerX = legalCenter.getX() * wireSegmentUnit_;
  const int centerY = legalCenter.getY() * wireSegmentUnit_;

  // Create shared root buffer and subnet (no tier prefix — shared by both tiers)
  sharedRootBuffer_
      = &clock_.addClockBuffer("clkbuf_0", root_buffer, centerX, centerY);
  // Only set topBufferName_ if not already set by forkRegisterClockNetwork().
  // RegisterTree builders have topBufferName_ = "clkbuf_regs_0_..." (separator buffer
  // that connects macro clock net to register clock net). Overwriting it with the
  // shared root name breaks clock connectivity for register trees in multi-clock designs.
  if (topBufferName_.empty()) {
    topBufferName_ = sharedRootBuffer_->getName();
  }
  addTreeLevelBuffer(sharedRootBuffer_);

  sharedRootSubNet_ = &clock_.addSubNet("clknet_0");
  sharedRootSubNet_->setTreeLevel(0);
  sharedRootSubNet_->addInst(*sharedRootBuffer_);
  treeBufLevels_++;

  sharedRootLocation_ = legalCenter;

  logger_->info(CTS, 856,
      "Shared root: tier={}, buffer={}, sinks={}, center=({:.2f}, {:.2f}), "
      "rootBufName='{}', topBufferName='{}', clockName='{}'",
      root_tier, root_buffer, allSinks.size(),
      legalCenter.getX(), legalCenter.getY(),
      sharedRootBuffer_->getName(), topBufferName_, clock_.getName());
}

// Runs the complete H-tree pipeline for a single tier:
//   1. Filter sinks to this tier
//   2. Clear per-tree state (mapLocationToSink_, topLevelSinksClustered_, topology)
//   3. Run clustering → CKMeans → wire segments → createClockSubNets
// When sharedRootBuffer_ is set, L1 branches attach to the shared root
// instead of creating a tier-specific root buffer.
void HTreeBuilder::runSingleTierTree(int tier)
{
  logger_->info(CTS, 851, "Building independent H-tree for tier {}", tier);

  // --- Step 0: Set tier naming prefix to avoid buffer name collisions ---
  // Tier 0 uses no prefix (backward compatible names: clkbuf_0, clkbuf_leaf_0, etc.)
  // Tier 1 uses "t1_" prefix: t1_clkbuf_0, t1_clkbuf_leaf_0, etc.
  tierPrefix_ = (tier == 0) ? "" : "t" + std::to_string(tier) + "_";

  // P6 equalize: track current tier + reset level counter for the sibling tier.
  currentBuildTier_ = tier;
  lastTierLevels_ = 0;

  // --- Step 1: Clear per-tree state ---
  mapLocationToSink_.clear();
  topLevelSinksClustered_.clear();
  topologyForEachLevel_.clear();
  clusterBufTarget_.clear();
  branchDelayTargets_.clear();
  perLevelGlobalMean_.clear();

  // wireSegmentUnit_ comes from TechChar (charBuf height × 10 / 2).
  // charBuf is selected from buf_list (highest max-cap buffer).
  // Do NOT override per-tier: in heterogeneous 3D with upper-only buf_list,
  // using bottom tier buffer height (5x larger for NG45) destroys H-tree grid.
  // Removed per-tier wireSegmentUnit recalculation (V57 hotfix, ).

  // --- Step 2: Filter sinks for this tier and cluster ---
  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  const int dbUnits = options_->getDbUnits();
  const float maxDiameter = (clusterDiameter * dbUnits) / wireSegmentUnit_;

  std::vector<std::pair<float, float>> topLevelSinks;
  std::vector<const ClockInst*> sinkInsts;
  initTopLevelSinksForTier(tier, topLevelSinks, sinkInsts);
  logger_->info(CTS, 852, "Tier {}: {} sinks after filtering", tier,
                topLevelSinks.size());

  if (topLevelSinks.empty()) {
    logger_->info(CTS, 853, "Tier {}: no sinks, skipping tree", tier);
    return;
  }

  // --- Step 4: Pre-sink clustering (Hilbert+Greedy or TARP) ---
  // V73 fix: skip clustering when the Fix-G macro flat bypass (Step 4b) will fire.
  // Clustering would create leaf buffers that the bypass never drives -> dangling
  // clock buffers -> macro CLK pins with no propagated clock.  See
  // macroFlatBypassApplies().
  const bool flatMacroBypass
      = macroFlatBypassApplies(static_cast<int>(topLevelSinks.size()), sinkInsts);
  if (flatMacroBypass) {
    logger_->info(CTS, 724,
        "V73: tier {} macro flat bypass detected before clustering - "
        "skipping preSinkClustering ({} macro sinks connect directly)",
        tier, topLevelSinks.size());
  } else {
    preSinkClustering(topLevelSinks, sinkInsts, maxDiameter, clusterSize);
  }

  // Compute sink region from clustered sinks
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;
  if (flatMacroBypass) {
    // V73: no clustered sinks exist; derive the region from the raw sinks.
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
  } else if (topLevelSinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
  } else {
    if (topLevelSinksClustered_.size() > 400
        && options_->getSinkClusteringLevels() > 0) {
      std::vector<std::pair<float, float>> secondLevelLocs;
      std::vector<const ClockInst*> secondLevelInsts;
      initSecondLevelSinks(secondLevelLocs, secondLevelInsts);
      preSinkClustering(secondLevelLocs, secondLevelInsts,
                        maxDiameter * 4,
                        std::ceil(std::sqrt(clusterSize)),
                        true);
    }
    sinkRegion_ = clock_.computeSinkRegionClustered(topLevelSinksClustered_);
  }
  logger_->info(CTS, 854, "Tier {} sink region: {}", tier, sinkRegion_);

  // --- Step 4b: Macro flat bypass (Fix G) ---
  // Pin3D: bottom-only CTS → macro net has ~7 sinks (<15) → single-buffer flat
  //   topology → macro latency ~563ps, balanced with registers (~442ps).
  // V58 per-tier CTS: all 21 SRAMs on clk net → H-tree builds → 1,095ps latency
  //   → 633ps gap vs registers → catastrophic hold violations.
  // Fix G: when MacroTree AND all sinks are block (SRAM) AND sink_count ≤ threshold,
  //   skip H-tree entirely → createSingleBufferClockNet (flat fan-out from root).
  //   This replicates Pin3D's behavior: root → direct wire → each macro.
  // CTS_MACRO_FLAT_MAX_SINKS env override (0=disable, default 50).
  {
    int macroFlatMax = 50;
    if (const char* e = std::getenv("CTS_MACRO_FLAT_MAX_SINKS")) {
      macroFlatMax = std::atoi(e);
    }
    const int nSinks = static_cast<int>(topLevelSinks.size());
    if (type_ == TreeType::MacroTree && macroFlatMax > 0 && nSinks <= macroFlatMax) {
      // Verify all sinks are block (SRAM/macro) instances
      bool allBlock = true;
      for (const auto* sinkInst : sinkInsts) {
        if (!sinkInst || !isMacroBlockInst(sinkInst->getDbInst())) {
          allBlock = false;
          break;
        }
      }
      if (allBlock) {
        logger_->info(CTS, 693,
            "Fix G: Macro flat bypass — {} sinks (all block), threshold={}, "
            "skipping H-tree → flat topology (Pin3D-equivalent)",
            nSinks, macroFlatMax);

        if (sharedRootSubNet_ != nullptr) {
          // Per-tier mode: shared root already exists (clkbuf_0_clk).
          // V58_CL Fix: Cross-tier macro sub-root buffer.
          // When macro sinks are on a different tier from the shared root,
          // connecting them directly creates long cross-tier wires with huge
          // latency spread (swerv DCCM: 475-1163ps from unbalanced wire RC).
          // Fix: create a per-tier sub-root buffer on the macro's tier,
          // connect it to the shared root, then connect macros to the sub-root.
          // Same-tier macros still connect directly (no change).
          const bool crossTier = (tier != sharedRootTier_
                                  && sharedRootTier_ >= 0);
          ClockSubNet* macroSubNet = sharedRootSubNet_;
          if (crossTier) {
            // Compute centroid of macro sinks on this tier
            float cx = 0, cy = 0;
            int cnt = 0;
            for (const ClockInst* inst : sinkInsts) {
              if (inst != nullptr) {
                cx += (double) inst->getX() / wireSegmentUnit_;
                cy += (double) inst->getY() / wireSegmentUnit_;
                cnt++;
              }
            }
            if (cnt > 0) { cx /= cnt; cy /= cnt; }
            // Create sub-root buffer on the macro's tier
            const std::string subRootMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), tier);
            Point<double> subRootLoc(cx, cy);
            Point<double> legalSubRoot = legalizeOneBuffer(
                subRootLoc, subRootMaster);
            const int sx = legalSubRoot.getX() * wireSegmentUnit_;
            const int sy = legalSubRoot.getY() * wireSegmentUnit_;
            const std::string subRootName
                = tierPrefix_ + "clkbuf_macro_subroot_"
                  + clock_.getName();
            ClockInst& subRootBuf = clock_.addClockBuffer(
                subRootName, subRootMaster, sx, sy);
            if (cts3dDb_ != nullptr) {
              cts3dDb_->setClockInstTier(subRootBuf, tier);
            }
            addTreeLevelBuffer(&subRootBuf);
            // Connect sub-root to shared root subnet
            sharedRootSubNet_->addInst(subRootBuf);
            // Create sub-net for macro sinks
            const std::string subNetName
                = tierPrefix_ + "clknet_macro_subroot_"
                  + clock_.getName();
            ClockSubNet& subNet = clock_.addSubNet(subNetName);
            subNet.addInst(subRootBuf);
            subNet.setLeafLevel(true);
            macroSubNet = &subNet;
            treeBufLevels_++;
            logger_->info(CTS, 696,
                "Fix G: cross-tier macro sub-root '{}' on tier {} "
                "at ({:.2f}, {:.2f}) for {} macro sinks "
                "(shared root on tier {})",
                subRootName, tier, legalSubRoot.getX(),
                legalSubRoot.getY(), nSinks, sharedRootTier_);
          } else {
            // Same tier: connect directly to shared root (original behavior)
            if (!sharedRootSubNet_->isLeafLevel()) {
              sharedRootSubNet_->setLeafLevel(true);
              treeBufLevels_++;
            }
          }
          for (const ClockInst* inst : sinkInsts) {
            if (inst != nullptr) {
              macroSubNet->addInst(*const_cast<ClockInst*>(inst));
            }
          }
          logger_->info(CTS, 694,
              "Fix G: tier {} added {} macro sinks {} '{}'",
              tier, nSinks,
              crossTier ? "via sub-root to" : "directly to shared root",
              crossTier ? (tierPrefix_ + "clkbuf_macro_subroot_"
                           + clock_.getName())
                        : sharedRootBuffer_->getName());
        } else {
          // Non-per-tier mode: use standard single-buffer flat topology.
          createSingleBufferClockNet();
          treeBufLevels_++;
        }
        return;
      }
    }
  }

  // --- Step 5: H-tree level topology (CKMeans bisection) ---
  // V58 Fix E: Macro tree depth clamp (fallback for non-flat macro trees).
  int effectiveMaxDepth = clockTreeMaxDepth_;
  if (type_ == TreeType::MacroTree) {
    int macroMaxDepth = 2;
    if (const char* e = std::getenv("CTS_MACRO_TREE_MAX_DEPTH")) {
      macroMaxDepth = std::atoi(e);
    }
    if (macroMaxDepth > 0 && macroMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 691,
          "V58 Fix E: Macro tree depth clamped {} -> {} (sinks={})",
          effectiveMaxDepth, macroMaxDepth,
          static_cast<int>(topLevelSinksClustered_.size()));
      effectiveMaxDepth = macroMaxDepth;
    }
  }
  // V58_CL: Per-tier register tree depth clamp.
  // Per-tier CTS creates deeper trees than single-tree Pin3D (path 14 vs 8)
  // because each tier builds an independent H-tree from the shared root.
  // Higher base latency → hold violations even for FFs with lp_target=0.
  // CTS_PER_TIER_REG_MAX_DEPTH (default 0=disabled) caps H-tree levels
  // for register trees in per-tier mode. Pin3D typically builds 7-8 levels.
  if (enablePerTierCts_ && type_ != TreeType::MacroTree) {
    int perTierRegMaxDepth = 0;
    if (const char* e = std::getenv("CTS_PER_TIER_REG_MAX_DEPTH")) {
      perTierRegMaxDepth = std::atoi(e);
    }
    if (perTierRegMaxDepth > 0
        && perTierRegMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 697,
          "V58_CL: Per-tier register tree depth clamped {} -> {} "
          "(tier={}, sinks={})",
          effectiveMaxDepth, perTierRegMaxDepth, tier,
          static_cast<int>(topLevelSinksClustered_.size()));
      effectiveMaxDepth = perTierRegMaxDepth;
    }
  }
  for (int level = 1; level <= effectiveMaxDepth; ++level) {
    const unsigned numSinksPerSubRegion
        = computeNumberOfSinksPerSubRegion(level);
    double regionWidth, regionHeight;
    computeSubRegionSize(level, regionWidth, regionHeight);

    if (isSubRegionTooSmall(regionWidth, regionHeight)) {
      if (options_->isFakeLutEntriesEnabled()) {
        const unsigned minIndex = 1;
        techChar_->createFakeEntries(minLengthSinkRegion_, minIndex);
        minLengthSinkRegion_ = 1;
      } else {
        break;
      }
    }
    computeLevelTopology(level, regionWidth, regionHeight);
    if (isNumberOfSinksTooSmall(numSinksPerSubRegion)) {
      break;
    }
  }

  if (topologyForEachLevel_.empty()) {
    // Too few sinks for H-tree — create single buffer for this tier
    logger_->info(CTS, 855, "Tier {}: too few clusters for H-tree, "
                  "creating single buffer", tier);
    createSingleBufferClockNet();
    treeBufLevels_++;
    return;
  }

  clock_.setMaxLevel(topologyForEachLevel_.size());

  if (options_->getObstructionAware()) {
    legalize();
  }

  // --- Step 6: Compute branch delay targets (for grouped delay/TAP) ---
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    computeBranchDelayTargets();
  }

  // P6 equalize: record this tier's H-tree level count for the sibling tier.
  lastTierLevels_ = static_cast<int>(topologyForEachLevel_.size());

  // --- Step 7: Create clock sub-nets (buffers, wires, grouped delay) ---
  createClockSubNets();
}

void HTreeBuilder::initTopLevelSinks(
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  clock_.forEachSink([&](const ClockInst& sink) {
    // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
    sinkLocations.emplace_back((double) sink.getX() / wireSegmentUnit_,
                               (double) sink.getY() / wireSegmentUnit_);
    sinkInsts.emplace_back(&sink);
  });
}

// Filter sinks to a specific tier for per-tier CTS.
// Only includes FFs whose ODB instance is on the requested tier
// (determined by Cts3DDatabase::getInstTier).
void HTreeBuilder::initTopLevelSinksForTier(
    int tier,
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  sinkInsts.clear();
  clock_.forEachSink([&](const ClockInst& sink) {
    if (cts3dDb_ == nullptr) return;
    odb::dbInst* dbInst = sink.getDbInst();
    if (dbInst == nullptr) return;
    int instTier = cts3dDb_->getInstTier(dbInst);
    if (instTier != tier) return;
    // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
    sinkLocations.emplace_back((double) sink.getX() / wireSegmentUnit_,
                               (double) sink.getY() / wireSegmentUnit_);
    sinkInsts.emplace_back(&sink);
  });
}

void HTreeBuilder::initSecondLevelSinks(
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  for (const auto& buf : topLevelSinksClustered_) {
    sinkLocations.emplace_back(buf.first, buf.second);
    const Point<double> bufPos(buf.first, buf.second);
    sinkInsts.emplace_back(mapLocationToSink_[bufPos]);
  }
}

void HTreeBuilder::computeBranchSinks(
    const LevelTopology& topology,
    const unsigned branchIdx,
    std::vector<std::pair<float, float>>& sinkLocations) const
{
  sinkLocations.clear();
  for (const Point<double>& point :
       topology.getBranchSinksLocations(branchIdx)) {
    sinkLocations.emplace_back(point.getX(), point.getY());
  }
}

// P3 dominant data-flow unit vector — PLACEMENT/NETLIST ONLY (no STA, no LP).
// For each clock sink FF, take its first SIGNAL data-input net, walk backward
// through combinational drivers (bounded hops) until a sequential (another sink)
// or macro/block driver is reached, and accumulate (capture_loc - launch_loc).
// The normalized sum is the design's dominant launch->capture direction.
void HTreeBuilder::computeDirSkewVector()
{
  dirSkewValid_ = false;
  dirSkewVecX_ = dirSkewVecY_ = 0.0;
  if (db_ == nullptr) {
    return;
  }
  std::unordered_set<odb::dbInst*> seqSet;
  clock_.forEachSink([&](const ClockInst& sink) {
    if (odb::dbInst* di = sink.getDbInst()) {
      seqSet.insert(di);
    }
  });
  if (seqSet.empty()) {
    return;
  }
  auto instCenter = [](odb::dbInst* inst, double& cx, double& cy) {
    odb::dbBox* bb = inst->getBBox();
    cx = 0.5 * (bb->xMin() + bb->xMax());
    cy = 0.5 * (bb->yMin() + bb->yMax());
  };
  const int maxHops = 6;
  double sumX = 0.0, sumY = 0.0;
  long nCross = 0;
  for (odb::dbInst* capInst : seqSet) {
    odb::dbNet* net = nullptr;
    for (auto* it : capInst->getITerms()) {
      if (it->getIoType() != odb::dbIoType::INPUT) continue;
      if (it->getSigType() != odb::dbSigType::SIGNAL) continue;
      odb::dbNet* n = it->getNet();
      if (n != nullptr && !n->isSpecial()) { net = n; break; }
    }
    if (net == nullptr) continue;
    odb::dbInst* srcInst = nullptr;
    for (int hop = 0; hop < maxHops && net != nullptr; ++hop) {
      odb::dbInst* drv = nullptr;
      for (auto* it : net->getITerms()) {
        if (it->getIoType() == odb::dbIoType::OUTPUT) { drv = it->getInst(); break; }
      }
      if (drv == nullptr) break;
      if (seqSet.count(drv) || isMacroBlockInst(drv)) { srcInst = drv; break; }
      odb::dbNet* nextNet = nullptr;
      for (auto* it : drv->getITerms()) {
        if (it->getIoType() != odb::dbIoType::INPUT) continue;
        if (it->getSigType() != odb::dbSigType::SIGNAL) continue;
        odb::dbNet* n = it->getNet();
        if (n != nullptr && !n->isSpecial() && n != net) { nextNet = n; break; }
      }
      net = nextNet;
    }
    if (srcInst == nullptr || srcInst == capInst) continue;
    double lx, ly, cx, cy;
    instCenter(srcInst, lx, ly);
    instCenter(capInst, cx, cy);
    sumX += (cx - lx);
    sumY += (cy - ly);
    ++nCross;
  }
  if (nCross == 0) return;
  const double mag = std::sqrt(sumX * sumX + sumY * sumY);
  if (mag < 1e-9) return;
  dirSkewVecX_ = sumX / mag;
  dirSkewVecY_ = sumY / mag;
  dirSkewValid_ = true;
  logger_->info(CTS, 864,
      "P3 dir-skew proxy: cross-pairs={}, launch->capture unit=({:+.3f},{:+.3f})",
      nCross, dirSkewVecX_, dirSkewVecY_);
}

void HTreeBuilder::refineBranchingPointsWithClustering(
    LevelTopology& topology,
    const unsigned level,
    const unsigned branchPtIdx1,
    const unsigned branchPtIdx2,
    const Point<double>& rootLocation,
    const std::vector<std::pair<float, float>>& sinks)
{
  CKMeans::Clustering clusteringEngine(
      sinks, rootLocation.getX(), rootLocation.getY(), logger_);

  // Pass per-sink LP targets to CKMeans.
  // When ckmeansTargetBeta_ > 0, CKMeans groups FFs with similar targets
  // into the same branch, enabling larger V31 wire shift differences.
  if (ckmeansTargetBeta_ > 0.0 && cts3dDb_ != nullptr
      && cts3dDb_->hasSkewTargets()) {
    std::vector<float> sinkTargets;
    sinkTargets.reserve(sinks.size());
    for (const auto& s : sinks) {
      const Point<double> sinkLoc(s.first, s.second);
      auto it = mapLocationToSink_.find(sinkLoc);
      if (it != mapLocationToSink_.end() && it->second != nullptr) {
        sinkTargets.push_back(
            static_cast<float>(getClockInstTarget(it->second)));
      } else {
        sinkTargets.push_back(0.0f);
      }
    }
    clusteringEngine.setSinkTargets(
        sinkTargets, static_cast<float>(ckmeansTargetBeta_));
  }

  Point<double>& branchPt1 = topology.getBranchingPoint(branchPtIdx1);
  Point<double>& branchPt2 = topology.getBranchingPoint(branchPtIdx2);

  std::vector<std::pair<float, float>> means;
  means.emplace_back(branchPt1.getX(), branchPt1.getY());
  means.emplace_back(branchPt2.getX(), branchPt2.getY());

  const unsigned cap
      = options_->getMaxFanout()
            ? (unsigned) (sinks.size() * 0.5)
            : (unsigned) (sinks.size() * options_->getClusteringCapacity());
  clusteringEngine.iterKmeans(
      1, means.size(), cap, 5, options_->getClusteringPower(), means);

  if (((int) options_->getNumStaticLayers() - (int) level) < 0) {
    branchPt1 = Point<double>(means[0].first, means[0].second);
    branchPt2 = Point<double>(means[1].first, means[1].second);
  }

  // Retrieve cluster assignments (used both for V31 wire shift and sink assignment)
  std::vector<std::vector<unsigned>> clusters;
  clusteringEngine.getClusters(clusters);

  // Useful-skew wire length adjustment (Fishburn 1990)
  // After k-means sets branch point positions, shift each branch point radially
  // from rootLocation based on its cluster's mean LP skew target.
  //
  // Physical meaning:
  //   shift_um = wireSkewScale_ * (T_cluster - T_global) / wireDelayPerUnit_
  //   T_cluster > T_global → branch farther from root → longer wire → later clock arrival
  //   T_cluster < T_global → branch closer to root  → shorter wire → earlier clock arrival
  //
  // This realizes the Fishburn LP-SAFETY schedule in the physical wire topology,
  // instead of patching with delay buffers after a balanced tree is built.
  if (wireSkewScale_ > 0.0 && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
      && clusters.size() >= 2) {
    double sumT[2] = {0.0, 0.0};
    int    cntT[2] = {0,   0  };
    for (int ci = 0; ci < 2; ++ci) {
      for (unsigned sinkIdx : clusters[ci]) {
        const Point<double> sinkLoc(sinks[sinkIdx].first, sinks[sinkIdx].second);
        auto it = mapLocationToSink_.find(sinkLoc);
        if (it != mapLocationToSink_.end() && it->second != nullptr) {
          // Unified target accessor (cluster buffer or individual FF)
          sumT[ci] += getClockInstTarget(it->second);
          cntT[ci]++;
        }
      }
    }
    const double meanT0     = (cntT[0] > 0) ? sumT[0] / cntT[0] : 0.0;
    const double meanT1     = (cntT[1] > 0) ? sumT[1] / cntT[1] : 0.0;
    const double globalMean = (meanT0 + meanT1) / 2.0;

    // Radially shift each branch point: positive shift = farther from root
    auto shiftBranchPoint = [&](Point<double>& bp, double shiftNs) {
      const double dx   = bp.getX() - rootLocation.getX();
      const double dy   = bp.getY() - rootLocation.getY();
      const double dist = std::sqrt(dx * dx + dy * dy);
      if (dist < 1e-6 || wireDelayPerUnit_ < 1e-12) {
        return;
      }
      const double shiftUm = wireSkewScale_ * shiftNs / wireDelayPerUnit_;
      // Widened clamp from 50% to 80% of branch length.
      // 50% limited wire-skew to ~10ps; LP requests 30-50ps differentials.
      const double maxShift    = 0.8 * dist;
      const double clampedShift = std::max(-maxShift, std::min(maxShift, shiftUm));
      const double scale = (dist + clampedShift) / dist;
      bp.setX(rootLocation.getX() + dx * scale);
      bp.setY(rootLocation.getY() + dy * scale);
    };

    shiftBranchPoint(branchPt1, meanT0 - globalMean);
    shiftBranchPoint(branchPt2, meanT1 - globalMean);

    // Log at debug level to avoid flooding; only log when shift is significant
    const double shift0Ps = wireSkewScale_ * (meanT0 - globalMean) / wireDelayPerUnit_ * 1000.0;
    const double shift1Ps = wireSkewScale_ * (meanT1 - globalMean) / wireDelayPerUnit_ * 1000.0;
    if (std::abs(shift0Ps) > 0.1 || std::abs(shift1Ps) > 0.1) {
      debugPrint(logger_, CTS, "skew_wire", 1,
                 "Level {} wire skew: cluster0 T={:+.1f}ps shift={:+.1f}um, "
                 "cluster1 T={:+.1f}ps shift={:+.1f}um",
                 level,
                 (meanT0 - globalMean) * 1000.0, shift0Ps / 1000.0 * wireDelayPerUnit_ / wireDelayPerUnit_,
                 (meanT1 - globalMean) * 1000.0, shift1Ps / 1000.0 * wireDelayPerUnit_ / wireDelayPerUnit_);
    }
  }
  // End wire-skew adjustment

  // ---- P3: directional asymmetric branch-point embedding (structural useful-skew) ----
  // Rigidly translate BOTH child branch points by +d along the level split axis.
  // Wire-conserving but redistributes: +axis child (downstream/capture) gets a
  // LONGER wire (later clock, borrows time on the dominant path); -axis child
  // (launch) gets SHORTER (earlier). Applied AFTER k-means (which overwrites
  // branchPt1/2 for every level since numStaticLayers default 0) and ONLY at upper
  // levels 1..dirSkewMaxLevel_ so deep/intra-cluster hold pairs are untouched.
  // alpha<0 flips the sign for empirical tuning. Env default off → winners identical.
  if (dirSkewEnable_ && dirSkewValid_
      && level >= 1 && level <= (unsigned) dirSkewMaxLevel_) {
    const bool   horiz = isHorizontal(level);
    const double s     = horiz ? dirSkewVecX_ : dirSkewVecY_;
    double d = dirSkewAlpha_ * topology.getLength() * s;
    const double capLU = dirSkewMaxShiftUm_
        * (double) options_->getDbUnits() / (double) techChar_->getLengthUnit();
    d = std::max(-capLU, std::min(capLU, d));
    if (std::abs(d) > 1e-9) {
      if (horiz) {
        branchPt1.setX(branchPt1.getX() + d);
        branchPt2.setX(branchPt2.getX() + d);
      } else {
        branchPt1.setY(branchPt1.getY() + d);
        branchPt2.setY(branchPt2.getY() + d);
      }
      const double shiftUm = d * (double) techChar_->getLengthUnit()
                             / (double) options_->getDbUnits();
      debugPrint(logger_, CTS, "dirskew", 1,
          "Level {} dir-skew: axis={} s={:+.3f} shift={:+.2f}um (cap {:.1f}um)",
          level, horiz ? "X" : "Y", s, shiftUm, dirSkewMaxShiftUm_);
    }
  }
  // End directional useful-skew

  unsigned movedSinks = 0;
  const double errorFactor = 1.2;
  for (int clusterIdx = 0; clusterIdx < clusters.size(); ++clusterIdx) {
    for (int elementIdx = 0; elementIdx < clusters[clusterIdx].size();
         ++elementIdx) {
      const unsigned sinkIdx = clusters[clusterIdx][elementIdx];
      const Point<double> sinkLoc(sinks[sinkIdx].first, sinks[sinkIdx].second);
      const double dist = clusterIdx == 0 ? computeDist(branchPt1, sinkLoc)
                                          : computeDist(branchPt2, sinkLoc);
      const double distOther = clusterIdx == 0
                                   ? computeDist(branchPt2, sinkLoc)
                                   : computeDist(branchPt1, sinkLoc);
      if (clusterIdx == 0) {
        topology.addSinkToBranch(branchPtIdx1, sinkLoc);
      } else {
        topology.addSinkToBranch(branchPtIdx2, sinkLoc);
      }

      if (dist >= distOther * errorFactor) {
        movedSinks++;
      }
    }
  }

  if (movedSinks > 0) {
    debugPrint(logger_,
               CTS,
               "clustering",
               1,
               " Out of {} sinks, {} sinks closer to other cluster.",
               sinks.size(),
               movedSinks);
  }
}

// Unified LP skew target accessor for any ClockInst.
// After preSinkClustering(), mapLocationToSink_ entries for multi-FF clusters
// point to cluster buffers (clkbuf_leaf_X), NOT individual FFs.
// Cluster buffer names are not in the LP CSV, so direct name lookup returns 0.
// This helper checks clusterBufTarget_ first (set at clustering time from FF members),
// then falls back to stripping "/clk_pin" suffix for individual FF lookup.
double HTreeBuilder::getClockInstTarget(ClockInst* inst) const
{
  if (inst == nullptr || cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) {
    return 0.0;
  }
  // Cluster buffer path: use precomputed mean of member FF targets
  auto cit = clusterBufTarget_.find(inst);
  if (cit != clusterBufTarget_.end()) {
    return cit->second;
  }
  // Individual FF path: ClockInst name = "inst_name/clk_pin" → strip suffix
  std::string nm = inst->getName();
  const auto sp = nm.rfind('/');
  if (sp != std::string::npos) nm = nm.substr(0, sp);
  return cts3dDb_->getSkewTarget(nm);
}

// Step 4: Compute mean arrival target per branch.
// Extended to ALL levels (not just leaf).
// This determines which branches need extra delay buffers at each level.
void HTreeBuilder::computeBranchDelayTargets()
{
  branchDelayTargets_.clear();
  globalMeanTarget_ = 0.0;
  perLevelGlobalMean_.clear();

  if (topologyForEachLevel_.empty() || cts3dDb_ == nullptr
      || !cts3dDb_->hasSkewTargets()) {
    return;
  }

  // Debug: per-branch LP target statistics
  struct BranchDetail {
    int level;
    unsigned idx;
    int count;
    double mean, bmin, bmax, stddev;
  };
  std::vector<BranchDetail> branchDetails;

  // Compute per-branch mean target at ALL levels (not just leaf).
  // At each level, a branch's mean target = average of all sink targets
  // reachable from that branch point downward.
  for (int levelIdx = 0;
       levelIdx < static_cast<int>(topologyForEachLevel_.size());
       ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];
    double levelTargetSum = 0.0;
    int levelSinkCount = 0;

    topology.forEachBranchingPoint([&](unsigned idx, Point<double>) {
      const auto& sinkLocs = topology.getBranchSinksLocations(idx);
      if (sinkLocs.empty()) {
        return;
      }

      double branchSum = 0.0, branchSumSq = 0.0;
      double branchMin = std::numeric_limits<double>::max();
      double branchMax = -std::numeric_limits<double>::max();
      int branchCount = 0;

      for (const auto& loc : sinkLocs) {
        // V58 Fix B: epsilon-tolerant lookup (was exact find → skipped
        // singleton macro clusters whose float→double key didn't match)
        ClockInst* sink = findSinkEps(loc);
        if (sink == nullptr) {
          continue;
        }
        const double tgt = getClockInstTarget(sink);
        branchSum   += tgt;
        branchSumSq += tgt * tgt;
        branchMin    = std::min(branchMin, tgt);
        branchMax    = std::max(branchMax, tgt);
        ++branchCount;
      }

      if (branchCount > 0) {
        const double meanTarget = branchSum / branchCount;
        const double var = (branchSumSq / branchCount) - meanTarget * meanTarget;
        const double stddev = (var > 0.0) ? std::sqrt(var) : 0.0;
        branchDelayTargets_[{levelIdx, idx}] = meanTarget;
        levelTargetSum += branchSum;
        levelSinkCount += branchCount;
        branchDetails.push_back(
            {levelIdx, idx, branchCount, meanTarget,
             branchMin, branchMax, stddev});
      }
    });

    if (levelSinkCount > 0) {
      perLevelGlobalMean_[levelIdx] = levelTargetSum / levelSinkCount;
    }
  }

  // Global mean across all sinks (backward compatible with V32c leaf-only)
  const int leafLevelIdx = static_cast<int>(topologyForEachLevel_.size()) - 1;
  auto it = perLevelGlobalMean_.find(leafLevelIdx);
  globalMeanTarget_ = (it != perLevelGlobalMean_.end()) ? it->second : 0.0;

  // Log per-level statistics
  for (int levelIdx = 0;
       levelIdx < static_cast<int>(topologyForEachLevel_.size());
       ++levelIdx) {
    int posCount = 0, negCount = 0;
    int branchCount = 0;
    const double levelMean
        = perLevelGlobalMean_.count(levelIdx)
            ? perLevelGlobalMean_[levelIdx] : 0.0;
    for (const auto& [key, meanTarget] : branchDelayTargets_) {
      if (key.first != levelIdx) continue;
      ++branchCount;
      const double delta = meanTarget - levelMean;
      if (delta > delayTargetThreshold_) ++posCount;
      else if (delta < -delayTargetThreshold_) ++negCount;
    }
    logger_->info(CTS, 401,
                  "Per-branch delay targets L{}: {} branches, "
                  "levelMean={:.4f}ns, {} need extra delay, {} need less delay",
                  levelIdx + 1, branchCount, levelMean, posCount, negCount);
  }

  // Per-branch LP target distribution detail (leaf level only)
  for (const auto& d : branchDetails) {
    if (d.level != leafLevelIdx) continue;
    // globalMean fix: absolute LP target for planned_N diagnostic
    const double delta = d.mean;
    int plannedN = 0;
    if (delta > 0.0 && singleBufDelay_ > 1e-9) {
      plannedN = std::max(0, std::min(
          static_cast<int>(std::round(delta / singleBufDelay_)),
          maxLeafDelayBufs_));
    }
    debugPrint(logger_, CTS, "HTree", 1,
                  "  Leaf branch {:2d}: sinks={}, "
                  "LP=[min={:.4f} mean={:.4f} max={:.4f} spread={:.4f} "
                  "stddev={:.4f}]ns, delta={:.4f}ns -> planned_N={}",
                  d.idx, d.count,
                  d.bmin, d.mean, d.bmax, d.bmax - d.bmin,
                  d.stddev, delta, plannedN);
  }
}

// Load per-FF hold budgets from timing graph CSV.
// For each capture FF (to_ff column), compute min(slack_min_ns) across all
// edges. This is the FF's hold budget: max additional delay we can add
// to its clock arrival without violating hold.
void HTreeBuilder::loadPerFfHoldBudgets(const std::string& timingGraphCsv)
{
  perFfHoldBudget_.clear();
  std::ifstream file(timingGraphCsv);
  if (!file.is_open()) {
    logger_->warn(CTS, 422,
                  "Cannot open timing graph for hold budgets: {}",
                  timingGraphCsv);
    return;
  }

  std::string line;
  std::getline(file, line);  // skip header

  int edgeCount = 0;
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string from_ff, to_ff, slack_max_str, slack_min_str;
    // CSV format: from_ff,to_ff,slack_max_ns,slack_min_ns,...
    std::getline(ss, from_ff, ',');
    std::getline(ss, to_ff, ',');
    std::getline(ss, slack_max_str, ',');
    std::getline(ss, slack_min_str, ',');

    if (to_ff.empty() || slack_min_str.empty()) continue;

    try {
      const double holdSlack = std::stod(slack_min_str);
      auto it = perFfHoldBudget_.find(to_ff);
      if (it == perFfHoldBudget_.end() || holdSlack < it->second) {
        perFfHoldBudget_[to_ff] = holdSlack;
      }
      ++edgeCount;
    } catch (...) {
      // skip malformed lines
    }
  }

  logger_->info(CTS, 423,
                "Hold budgets loaded: {} capture FFs from {} edges",
                perFfHoldBudget_.size(), edgeCount);
}

// Load IO edge hold budgets from IO timing CSV.
// For PI_TO_FF rows, the hold slack represents the budget for adding delay
// to the capture FF's clock arrival. Merged into perFfHoldBudget_ (take min).
void HTreeBuilder::loadIoHoldBudgets(const std::string& ioCsv)
{
  std::ifstream file(ioCsv);
  if (!file.is_open()) {
    logger_->warn(CTS, 429,
                  "Cannot open IO timing CSV for hold budgets: {}", ioCsv);
    return;
  }

  std::string line;
  std::getline(file, line);  // skip header: edge_type,port_name,ff_name,...

  int ioEdgeCount = 0;
  int mergedCount = 0;
  int newCount = 0;
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string edgeType, portName, ffName, slackSetup, slackHold;
    // CSV: edge_type,port_name,ff_name,slack_setup_ns,slack_hold_ns
    std::getline(ss, edgeType, ',');
    std::getline(ss, portName, ',');
    std::getline(ss, ffName, ',');
    std::getline(ss, slackSetup, ',');
    std::getline(ss, slackHold, ',');

    // Only PI_TO_FF hold matters: adding relay to capture FF hurts hold.
    // FF_TO_PO: adding relay to launch FF helps hold (not a concern).
    if (edgeType != "PI_TO_FF") continue;
    if (ffName.empty() || slackHold.empty()) continue;

    try {
      const double holdSlack = std::stod(slackHold);
      auto it = perFfHoldBudget_.find(ffName);
      if (it == perFfHoldBudget_.end()) {
        perFfHoldBudget_[ffName] = holdSlack;
        ++newCount;
      } else if (holdSlack < it->second) {
        it->second = holdSlack;
        ++mergedCount;
      }
      ++ioEdgeCount;
    } catch (...) {
      // skip malformed lines
    }
  }

  logger_->info(CTS, 430,
                "IO hold budgets: {} PI→FF edges, {} new FFs, "
                "{} tightened (total hold budget FFs: {})",
                ioEdgeCount, newCount, mergedCount,
                perFfHoldBudget_.size());
}

void HTreeBuilder::createClockSubNets()
{
  Point<double> legalCenter{0.0, 0.0};
  ClockSubNet* rootSubNetPtr;

  if (sharedRootBuffer_ != nullptr) {
    // Per-tier mode: reuse the shared root buffer created by createSharedRoot()
    rootSubNetPtr = sharedRootSubNet_;
    legalCenter = sharedRootLocation_;
    logger_->info(utl::CTS, 857,
                  "Per-tier sub-tree (tier prefix='{}') using shared root "
                  "at ({:.2f}, {:.2f})",
                  tierPrefix_, legalCenter.getX(), legalCenter.getY());
  } else {
    // Unified mode: create root buffer here (original behavior)
    Point<double> center = sinkRegion_.getCenter();
    // V58_CL2: CTS_FORCE_ROOT_TIER overrides dominant-tier selection (same as createSharedRoot)
    int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
    if (const char* e = std::getenv("CTS_FORCE_ROOT_TIER")) {
      int forced = std::atoi(e);
      if (forced >= 0 && forced != root_tier) {
        logger_->info(CTS, 701,
            "V58_CL2: Root tier forced {} -> {} (CTS_FORCE_ROOT_TIER, unified)",
            root_tier, forced);
        root_tier = forced;
      }
    }
    const std::string root_buffer
        = cts3dDb_->getBufferForTier(options_->getRootBuffer(), root_tier);
    logger_->info(utl::CTS, 315,
                  "3D-CTS root buffer: tier={}, buffer={}",
                  root_tier, root_buffer);
    legalCenter = legalizeOneBuffer(center, root_buffer);
    sinkRegion_.setCenter(legalCenter);
    commitMoveLoc(center, legalCenter);
    const int centerX = legalCenter.getX() * wireSegmentUnit_;
    const int centerY = legalCenter.getY() * wireSegmentUnit_;

    const std::string rootBufName = tierPrefix_ + "clkbuf_0";
    const std::string rootNetName = tierPrefix_ + "clknet_0";

    ClockInst& rootBuffer
        = clock_.addClockBuffer(rootBufName, root_buffer, centerX, centerY);

    if (topBufferName_.empty()) {
      topBufferName_ = rootBuffer.getName();
    }

    // clang-format off
    if (center != legalCenter) {
      debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
                 "root {}: {} => {}", rootBufName, center, legalCenter);
    } else {
      debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
                 "root {}: {}", rootBufName, center);
    }
    // clang-format on

    addTreeLevelBuffer(&rootBuffer);
    ClockSubNet& rootClockSubNet = clock_.addSubNet(rootNetName);
    rootClockSubNet.setTreeLevel(0);
    rootClockSubNet.addInst(rootBuffer);
    treeBufLevels_++;

    rootSubNetPtr = &rootClockSubNet;

    // P2 dominant-tier graft: capture this tier's own root. Majority tier →
    // remember its root SUBNET (graft point); minority tier → remember its root
    // BUFFER (grafted onto majority). Guard graftMajorityTier_>=0 → only in P2.
    if (enableDomTierGraft_ && graftMajorityTier_ >= 0) {
      if (currentBuildTier_ == graftMajorityTier_) {
        graftMajRootSubNet_ = &rootClockSubNet;
      } else {
        graftMinRootBuf_ = &rootBuffer;
      }
    }
  }

  // P6 equalize: pad the shallow (minority) tier's trunk with series buffers so
  // its root-to-leaf buffer-level depth matches the majority tier, collapsing
  // the tier-to-tier latency seam (ariane ~2.5%). Env-gated; majority tier and
  // balanced designs get deficit<=0 -> no-op. NOT a delay tap: these trunk
  // buffers drive the ENTIRE minority subtree (every bottom FF gains the SAME
  // common-mode delay -> intra-tier hold pairs see delta=0 -> hold-safe).
  if (sharedRootBuffer_ != nullptr && enablePerTierEqualize_
      && currentBuildTier_ == equalizeActiveTier_) {
    // CTS_PER_TIER_EQUALIZE_EXTRA: over-pad beyond level-match to compensate for
    // per-level delay mismatch (NG45 bottom trunk buf delay < ASAP7 upper level).
    int eqExtra = 0;
    if (const char* e = std::getenv("CTS_PER_TIER_EQUALIZE_EXTRA"))
      eqExtra = std::atoi(e);
    int deficit = equalizeTargetLevels_ + eqExtra
                  - static_cast<int>(topologyForEachLevel_.size());
    deficit = std::max(0, std::min(deficit, perTierEqualizeMaxBufs_));
    if (deficit > 0) {
      const std::string eqBufMaster = cts3dDb_->getBufferForTier(
          options_->getRootBuffer(), currentBuildTier_);
      Point<double> eqLoc = legalCenter;
      ClockSubNet* prevSub = rootSubNetPtr;   // starts at shared root subnet
      for (int b = 0; b < deficit; ++b) {
        Point<double> legalEq = legalizeOneBuffer(eqLoc, eqBufMaster);
        const int ex = (int) (legalEq.getX() * wireSegmentUnit_);
        const int ey = (int) (legalEq.getY() * wireSegmentUnit_);
        const std::string eqName
            = tierPrefix_ + "clkbuf_eq_" + std::to_string(b);
        ClockInst& eqBuf = clock_.addClockBuffer(eqName, eqBufMaster, ex, ey);
        if (cts3dDb_ != nullptr)
          cts3dDb_->setClockInstTier(eqBuf, currentBuildTier_);
        addTreeLevelBuffer(&eqBuf);
        prevSub->addInst(eqBuf);              // fanout of upstream trunk net
        const std::string eqNetName
            = tierPrefix_ + "clknet_eq_" + std::to_string(b);
        ClockSubNet& eqSub = clock_.addSubNet(eqNetName);
        eqSub.setTreeLevel(0);
        eqSub.addInst(eqBuf);
        prevSub = &eqSub;
        treeBufLevels_++;
        eqLoc = legalEq;
      }
      rootSubNetPtr = prevSub;                // L1 branches hang off padded trunk
      legalCenter = eqLoc;
      logger_->info(CTS, 858,
          "P6 equalize: tier {} padded with {} trunk buffers "
          "(target={} levels, own={} levels)",
          currentBuildTier_, deficit, equalizeTargetLevels_,
          static_cast<int>(topologyForEachLevel_.size()));
    }
  }

  // First level...
  LevelTopology& topLevelTopology = topologyForEachLevel_[0];
  bool isFirstPoint = true;
  topLevelTopology.forEachBranchingPoint([&](unsigned idx,
                                             Point<double> branchPoint) {
    // If the branch point has no sinks that will be connected to
    // it don't create a clock sub net for it
    if (topLevelTopology.getBranchSinksLocations(idx).empty()) {
      return;
    }
    // Replaced with Cts3DDatabase calls
    const int branch_tier = cts3dDb_->getDominantTier(
        topLevelTopology.getBranchSinksLocations(idx), mapLocationToSink_);
    const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
        options_->getRootBuffer(), branch_tier);
    logger_->info(utl::CTS, 316,
                  "3D-CTS branch L1 idx={}: tier={}, buffer={}",
                  idx, branch_tier, branch_root_buffer);

    // Cross-tier HB delay consideration
    // Elmore 50% delay: t_via = 0.693 * R_HB * C_HB (negligible: ~0 with C=0)
    // In shared-root mode, root_tier comes from the shared root buffer's tier.
    const int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
    if (branch_tier != root_tier) {
      const double hb_res = cts3dDb_->getHbtResistance();
      const double hb_cap = cts3dDb_->getHbtCapacitance();
      const double hb_delay = 0.693 * hb_res * hb_cap;

      const double wire_res_per_unit = cts3dDb_->getResPerDBU(branch_tier) * wireSegmentUnit_;
      const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(branch_tier) * wireSegmentUnit_;
      const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
          wire_res_per_unit, wire_cap_per_unit);

      logger_->warn(CTS, 362,
                    "Cross-tier branch L1: branch_tier={}, root_tier={}, "
                    "HB_delay={:.3e}s, equiv_dist={:.1f}um",
                    branch_tier, root_tier, hb_delay, hb_equivalent_dist);
    }

    Point<double> legalBranchPoint
        = legalizeOneBuffer(branchPoint, branch_root_buffer);
    commitMoveLoc(branchPoint, legalBranchPoint);

    // clang-format off
    if (branchPoint != legalBranchPoint) {
      debugPrint(logger_, CTS, "legalizer", 2,
		 "createClockSubNets level 1 clk_buf_1_{}_ : {} => {}",
		 std::to_string(idx), branchPoint, legalBranchPoint);
    } else {
      debugPrint(logger_, CTS, "legalizer", 2, 
		 "createClockSubNets level 1 clk_buf_1_{}_ : {}",
		 std::to_string(idx), branchPoint);
    }
    // clang-format on

    SegmentBuilder builder(tierPrefix_ + "clkbuf_1_" + std::to_string(idx) + "_",
                           tierPrefix_ + "clknet_1_" + std::to_string(idx) + "_",
                           legalCenter,  // center may have moved, don't use
                                         // sinkRegion_.getCenter()
                           legalBranchPoint,
                           topLevelTopology.getWireSegments(),
                           clock_,
                           *rootSubNetPtr,
                           *techChar_,
                           wireSegmentUnit_,
                           this,
                           db_,
                           branch_tier);
    if (!options_->getTreeBuffer().empty()) {
      // Replaced mapBufferMasterToTier with Cts3DDatabase
      const std::string tree_buffer = cts3dDb_->getBufferForTier(
          options_->getTreeBuffer(), branch_tier);
      builder.build(tree_buffer);
    } else {
      builder.build();
    }
    if (topologyForEachLevel_.size() == 1) {
      builder.forceBufferInSegment(branch_root_buffer);
    }
    if (isFirstPoint) {
      treeBufLevels_ += builder.getNumBufferLevels();
      isFirstPoint = false;
    }

    // Insert intermediate delay buffers at Level 1.
    // If this branch's subtree mean target > level globalMean, add N delay
    // buffers in series to shift arrival for the entire subtree.
    ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
    // gate with enableMidDelayBufs_ (default off in V39)
    if (enableMidDelayBufs_ && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
        && singleBufDelay_ > 1e-9) {
      auto tgtIt = branchDelayTargets_.find({0, idx});  // level 0
      const double levelMean
          = perLevelGlobalMean_.count(0) ? perLevelGlobalMean_[0] : 0.0;
      if (tgtIt != branchDelayTargets_.end()) {
        const double delta = tgtIt->second - levelMean;
        int midN = 0;
        if (delta > 0.0) {
          midN = static_cast<int>(std::round(delta / singleBufDelay_));
          midN = std::max(0, std::min(midN, maxLeafDelayBufs_));
        }
        if (midN > 0) {
          const int bx = legalBranchPoint.getX() * wireSegmentUnit_;
          const int by = legalBranchPoint.getY() * wireSegmentUnit_;
          const std::string delayBufMaster = branch_root_buffer;
          ClockSubNet* curSub = branchDrivingSub;
          for (int b = 0; b < midN; ++b) {
            const std::string bufName = tierPrefix_ + "clkbuf_mid_dly_1_"
                + std::to_string(idx) + "_" + std::to_string(b);
            ClockInst& buf = clock_.addClockBuffer(
                bufName, delayBufMaster, bx, by);
            cts3dDb_->setClockInstTier(buf, branch_tier);
            addTreeLevelBuffer(&buf);
            curSub->addInst(buf);
            ClockSubNet& nextSub = clock_.addSubNet(
                tierPrefix_ + "clknet_mid_dly_1_" + std::to_string(idx)
                + "_" + std::to_string(b));
            nextSub.addInst(buf);
            curSub = &nextSub;
          }
          branchDrivingSub = curSub;
          logger_->info(CTS, 411,
                        "Intermediate delay L1 branch {}: "
                        "target={:.4f}ns, levelMean={:.4f}ns, "
                        "delta={:.4f}ns, N={}",
                        idx, tgtIt->second, levelMean, delta, midN);
        }
      }
    }
    topLevelTopology.setBranchDrivingSubNet(idx, *branchDrivingSub);
  });

  // Others...
  for (int levelIdx = 1; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];
    isFirstPoint = true;
    topology.forEachBranchingPoint([&](unsigned idx,
                                       Point<double> branchPoint) {
      // If the branch point has no sinks that will be connected
      // to it don't create a clock sub net for it
      if (topology.getBranchSinksLocations(idx).empty()) {
        return;
      }
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
      LevelTopology& parentTopology = topologyForEachLevel_[levelIdx - 1];
      Point<double> parentPoint = parentTopology.getBranchingPoint(parentIdx);

      // Replaced with Cts3DDatabase calls
      const int branch_tier = cts3dDb_->getDominantTier(
          topology.getBranchSinksLocations(idx), mapLocationToSink_);
      const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
          options_->getRootBuffer(), branch_tier);
      debugPrint(logger_, CTS, "HTree", 1,
                    "3D-CTS branch L{} idx={}: tier={}, buffer={}",
                    levelIdx+1, idx, branch_tier, branch_root_buffer);

      // Cross-tier HB delay consideration
      // Elmore 50% delay: t_via = 0.693 * R_HB * C_HB (negligible: ~0 with C=0)
      const int parent_tier = cts3dDb_->getDominantTier(
          parentTopology.getBranchSinksLocations(parentIdx), mapLocationToSink_);
      if (branch_tier != parent_tier) {
        const double hb_res = cts3dDb_->getHbtResistance();
        const double hb_cap = cts3dDb_->getHbtCapacitance();
        const double hb_delay = 0.693 * hb_res * hb_cap;

        const double wire_res_per_unit = cts3dDb_->getResPerDBU(branch_tier) * wireSegmentUnit_;
        const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(branch_tier) * wireSegmentUnit_;
        const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
            wire_res_per_unit, wire_cap_per_unit);

        logger_->warn(CTS, 363,
                      "Cross-tier branch L{}: branch_tier={}, parent_tier={}, "
                      "HB_delay={:.3e}s, equiv_dist={:.1f}um",
                      levelIdx+1, branch_tier, parent_tier, hb_delay, hb_equivalent_dist);
      }

      Point<double> legalBranchPoint
          = legalizeOneBuffer(branchPoint, branch_root_buffer);
      commitMoveLoc(branchPoint, legalBranchPoint);

      // clang-format off
      if (branchPoint != legalBranchPoint) {
	debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets level {} "
		   "{} : {} => {}", levelIdx+1, tierPrefix_ + "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint, legalBranchPoint);
      } else {
	debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets level {} "
		   "{} : {}", levelIdx+1, tierPrefix_ + "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint);
      }
      // clang-format on

      SegmentBuilder builder(tierPrefix_ + "clkbuf_" + std::to_string(levelIdx + 1) + "_"
                                 + std::to_string(idx) + "_",
                             tierPrefix_ + "clknet_" + std::to_string(levelIdx + 1) + "_"
                                 + std::to_string(idx) + "_",
                             parentPoint,
                             legalBranchPoint,
                             topology.getWireSegments(),
                             clock_,
                             *parentTopology.getBranchDrivingSubNet(parentIdx),
                             *techChar_,
                             wireSegmentUnit_,
                             this,
                             db_,
                             branch_tier);

      // Set clock tree level the first time only.
      if (builder.getDrivingSubNet()->getTreeLevel() < 0) {
        builder.getDrivingSubNet()->setTreeLevel(levelIdx);
      }

      if (!options_->getTreeBuffer().empty()) {
        // Replaced mapBufferMasterToTier with Cts3DDatabase
        const std::string tree_buffer = cts3dDb_->getBufferForTier(
            options_->getTreeBuffer(), branch_tier);
        builder.build(tree_buffer);
      } else {
        builder.build();
      }
      if (levelIdx == topologyForEachLevel_.size() - 1) {
        builder.forceBufferInSegment(branch_root_buffer);
      }
      if (isFirstPoint) {
        treeBufLevels_ += builder.getNumBufferLevels();
        isFirstPoint = false;
      }

      // Insert intermediate delay buffers at Level 2+.
      // Skip the leaf level (handled by existing N-chain code below).
      const int leafLevel
          = static_cast<int>(topologyForEachLevel_.size()) - 1;
      ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
      // gate with enableMidDelayBufs_ (default off in V39)
      if (enableMidDelayBufs_ && levelIdx < leafLevel
          && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
          && singleBufDelay_ > 1e-9) {
        auto tgtIt = branchDelayTargets_.find({levelIdx, idx});
        const double levelMean
            = perLevelGlobalMean_.count(levelIdx)
                ? perLevelGlobalMean_[levelIdx] : 0.0;
        if (tgtIt != branchDelayTargets_.end()) {
          const double delta = tgtIt->second - levelMean;
          int midN = 0;
          if (delta > 0.0) {
            midN = static_cast<int>(std::round(delta / singleBufDelay_));
            midN = std::max(0, std::min(midN, maxLeafDelayBufs_));
          }
          if (midN > 0) {
            const int bx = legalBranchPoint.getX() * wireSegmentUnit_;
            const int by = legalBranchPoint.getY() * wireSegmentUnit_;
            const std::string delayBufMaster = branch_root_buffer;
            ClockSubNet* curSub = branchDrivingSub;
            for (int b = 0; b < midN; ++b) {
              const std::string bufName = tierPrefix_ + "clkbuf_mid_dly_"
                  + std::to_string(levelIdx + 1) + "_"
                  + std::to_string(idx) + "_" + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(
                  bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);
              curSub->addInst(buf);
              ClockSubNet& nextSub = clock_.addSubNet(
                  tierPrefix_ + "clknet_mid_dly_" + std::to_string(levelIdx + 1)
                  + "_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);
              curSub = &nextSub;
            }
            branchDrivingSub = curSub;
            logger_->info(CTS, 412,
                          "Intermediate delay L{} branch {}: "
                          "target={:.4f}ns, levelMean={:.4f}ns, "
                          "delta={:.4f}ns, N={}",
                          levelIdx + 1, idx, tgtIt->second,
                          levelMean, delta, midN);
          }
        }
      }
      topology.setBranchDrivingSubNet(idx, *branchDrivingSub);
    });
  }

  LevelTopology& leafTopology = topologyForEachLevel_.back();
  // leafLevelIdx removed: no longer needed (branchDelayTargets_ lookup replaced by inline absTgt)
  unsigned numSinks = 0;
  unsigned numDelayBufs = 0;

  // Debug: LP-CTS gap tracking accumulators (reported in CTS-408 after loop).
  // totalLpDelay   = sum(branchMean * sinkCount) over all leaf branches
  // totalDelivered = sum(N * singleBufDelay * sinkCount) for fired branches only
  // Coverage = totalDelivered / totalLpDelay * 100%: how much LP plan was implemented.
  double totalLpDelay_debug   = 0.0;
  double totalDelivered_debug = 0.0;
  int totalBranchSinks_debug  = 0;
  int totalBranches_debug     = 0;
  int firedBranches_debug     = 0;

  leafTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        ClockSubNet* subNet = leafTopology.getBranchDrivingSubNet(idx);
        // If no clock sub net was created for a leaf branch point no sinks
        // connect to it, so just skip.
        if (subNet == nullptr) {
          return;
        }

        // Step 4: Insert delay buffers if branch LP target > globalMean
        // N-buffer chain based on ABSOLUTE LP-TNS target.
        //   Problem (found after V32b run): Absolute targeting causes hold violations.
        //   Branch differentials up to 36-48ps (N=1 vs N=4) violated hold on short paths.
        //   LP guaranteed hold only for per-FF delays; branch-mean averaging breaks this.
        // DELTA-BASED N-buffer chain (relative to globalMean).
        //   N = round(max(0, branchMean - globalMean) / singleBufDelay_)
        //   Limits max inter-branch differential to (maxBranchDelta - 0) ≈ 12-16ps (1-2 bufs).
        //   LP-TNS result (asap7/aes): globalMean=26.6ps, max delta=15.8ps → N≤1 per branch.
        //   Chain: subNet → buf_0 → sub_0 → ... → FFs
        // Debug: branchCnt hoisted to lambda scope for CTS-408 gap tracking
        double absTgt = 0.0;
        int branchCnt = 0;
        if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
            && singleBufDelay_ > 1e-9) {
          const auto& bSinkLocs = leafTopology.getBranchSinksLocations(idx);
          double branchSum = 0.0;
          for (const auto& bLoc : bSinkLocs) {
            // V58 Fix B: epsilon-tolerant lookup for branch target calculation
            ClockInst* bSink = findSinkEps(bLoc);
            if (bSink == nullptr) {
              continue;
            }
            // Unified target accessor (cluster buffer or individual FF)
            branchSum += getClockInstTarget(bSink);
            ++branchCnt;
          }
          if (branchCnt > 0) {
            absTgt = branchSum / branchCnt;
          }
        }

        // Grouped Delay Chain.
        // Instead of per-FF relay (1 chain per FF → ~192 nets → GRT-0183),
        // group FFs within each leaf cluster by quantized LP target delta.
        // Group 0 (delta < min): direct to leaf buffer.
        // Group k>0: shared chain of k delay buffers → all group FFs.
        // ~24 group nets (8 branches × ~3 groups) → no GRT-0183.
        //
        // Structure: branch → leaf_buffer → FFs (group 0, direct)
        //                                 → [grp_buf_0] → FFs (group 1)
        //                                 → [grp_buf_0] → [grp_buf_1] → FFs (group 2)
        if (enableGroupedDelay_ && cts3dDb_ != nullptr
            && cts3dDb_->hasSkewTargets() && perFfBufDelay_ > 1e-9) {
          // static counter for globally unique buffer names.
          // Multiple sub-trees (register tree, macro tree, etc.) each call
          // this code with idx=0,1,... → name collision without global ID.
          static int grpGlobalId = 0;
          const int myGrpId = grpGlobalId++;

          subNet->setLeafLevel(true);
          int branchGroupedFFs = 0, branchGroupBufs = 0;
          int branchDirectFFs = 0, branchHoldSkip = 0;

          const int branch_tier = cts3dDb_->getDominantTier(
              leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
          const std::string bufMaster = cts3dDb_->getBufferForTier(
              options_->getRootBuffer(), branch_tier);

          // Helper: get FF name (strip /CLK suffix)
          auto getFFName = [](ClockInst* ff) -> std::string {
            std::string nm = ff->getName();
            const auto sp = nm.rfind('/');
            if (sp != std::string::npos) nm = nm.substr(0, sp);
            return nm;
          };

          // Helper: get per-FF target delta and hold-safe depth
          struct FfGroupInfo {
            ClockInst* inst;
            double delta;
            int depth;  // quantized, hold-capped
          };

          // Collect all FFs in this branch with their group assignment.
          // When enableClusterUniform_=true, all FFs in the cluster share
          // a single depth k_c computed from the cluster median LP target delta.
          // This creates exactly 1 chain per cluster (vs up to MAX_DEPTH chains
          // in V49), reducing buffer count by ~60% with comparable timing.
          // Hold safety: k_c is capped by the minimum hold budget across all FFs.
          auto collectFFs = [&](ClockInst* driverInst, ClockSubNet* driverSubNet,
                                std::vector<FfGroupInfo>& ffInfos,
                                std::vector<ClockInst*>& memberFFs) {
            if (enableClusterUniform_) {
              // Cluster-Uniform mode — compute median delta, assign to all FFs.
              // Pass 1: collect deltas and hold budgets
              std::vector<double> deltas;
              deltas.reserve(memberFFs.size());
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                // globalMean fix: absolute LP target for TAP delivery
                double ffDelta = cts3dDb_->getSkewTarget(ffName);
                // Debug: trace per-FF target lookup (cluster-uniform path)
                if (ffDelta > 0.01) {
                  logger_->info(CTS, 711,
                      "collectFFs cluster-uniform: '{}' target={:.4f}ns",
                      ffName, ffDelta);
                }
                deltas.push_back(ffDelta);
              }
              // Compute median delta for the cluster
              std::vector<double> sorted = deltas;
              std::sort(sorted.begin(), sorted.end());
              double medianDelta = sorted[sorted.size() / 2];

              // Quantize cluster depth from median
              int clusterDepth = 0;
              if (medianDelta > groupedDelayMinDelta_) {
                clusterDepth = static_cast<int>(
                    std::round(medianDelta / perFfBufDelay_));
                clusterDepth = std::max(
                    0, std::min(clusterDepth, groupedDelayMaxDepth_));
              }

              // Removed redundant C++ holdSkip.
              // LP already enforces hold via hold edges + PI-hold-clip + sigma guard.
              // C++ holdSkip was double-guarding with a more conservative local metric,
              // blocking 60-70% of grouped delay delivery. See CLAUDE.md for details.

              // Pass 2: assign uniform depth to all FFs
              for (size_t i = 0; i < memberFFs.size(); ++i) {
                // V58 Fix C: Macro sinks (SRAM) get depth=0 — tight hold margin.
                int ffDepth = clusterDepth;
                odb::dbInst* dInst = memberFFs[i]->getDbInst();
                if (isMacroBlockInst(dInst)) {
                  ffDepth = 0;
                }
                ffInfos.push_back({memberFFs[i], deltas[i], ffDepth});
              }
            } else {
              // Per-FF depth (original behavior)
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                double ffTarget = cts3dDb_->getSkewTarget(ffName);
                // Debug: trace per-FF target lookup to diagnose delivery gap
                // (3058 non-zero LP targets but only ~29 grouped FFs)
                if (ffTarget > 0.01) {
                  debugPrint(logger_, CTS, "HTree", 1,
                      "collectFFs per-FF: '{}' target={:.4f}ns depth_will_be={}",
                      ffName, ffTarget,
                      static_cast<int>(std::round(ffTarget / perFfBufDelay_)));
                }
                // globalMean fix: absolute LP target
                double ffDelta = ffTarget;

                // Quantize to buffer delay steps
                int depth = 0;
                if (ffDelta > groupedDelayMinDelta_) {
                  depth = static_cast<int>(
                      std::round(ffDelta / perFfBufDelay_));
                  depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
                }

                // Removed redundant holdSkip (same as above).
                // V58 Fix C: Macro sinks (SRAM) get depth=0.
                odb::dbInst* dInst = ff->getDbInst();
                if (isMacroBlockInst(dInst)) {
                  depth = 0;
                }
                ffInfos.push_back({ff, ffDelta, depth});
              }
            }
          };

          // V58_CL: Hold-budget-aware TAP depth cap.
          // After collectFFs assigns depth from LP targets, check each FF's
          // hold budget from perFfHoldBudget_ (loaded from timing graph).
          // If depth * perFfBufDelay_ would exceed the FF's hold budget,
          // cap depth to floor(holdBudget / perFfBufDelay_).
          // This protects FFs where LP missed hold constraints (pruned edges)
          // but the actual timing graph has tight PI→FF hold paths.
          // CTS_ENABLE_HOLD_BUDGET_CAP=1 (default ON when hold budgets loaded).
          bool enableHoldBudgetCap = !perFfHoldBudget_.empty();
          if (const char* e = std::getenv("CTS_ENABLE_HOLD_BUDGET_CAP")) {
            enableHoldBudgetCap = (std::atoi(e) != 0);
          }
          auto applyHoldBudgetCap = [&](std::vector<FfGroupInfo>& ffInfos) {
            if (!enableHoldBudgetCap || perFfHoldBudget_.empty()) return;
            int nCapped = 0;
            for (auto& fi : ffInfos) {
              if (fi.depth <= 0) continue;
              std::string ffName = getFFName(fi.inst);
              auto it = perFfHoldBudget_.find(ffName);
              if (it == perFfHoldBudget_.end()) continue;
              double holdBudgetNs = it->second;
              // Max depth that fits within hold budget
              int maxSafeDepth = static_cast<int>(
                  std::floor(holdBudgetNs / perFfBufDelay_));
              if (maxSafeDepth < 0) maxSafeDepth = 0;
              if (fi.depth > maxSafeDepth) {
                fi.depth = maxSafeDepth;
                ++nCapped;
              }
            }
            if (nCapped > 0) {
              logger_->info(CTS, 698,
                  "V58_CL: Hold budget cap reduced TAP depth for {} FFs "
                  "(bufDelay={:.3f}ns)",
                  nCapped, perFfBufDelay_);
            }
          };

          // Build grouped delay chains for collected FFs
          // No local numSinks here. Uses outer numSinks (line 3002)
          // so buildGroupChains [&] captures the correct per-sink counter.
          // Previously an inner numSinks was declared AFTER this lambda,
          // so the lambda captured the outer numSinks (always 0) → duplicate
          // buffer names across cluster buffers in same branch → CTS-0499.
          auto buildGroupChains = [&](ClockInst* driverInst,
                                      ClockSubNet* driverSubNet,
                                      std::vector<FfGroupInfo>& ffInfos) {
            // Group FFs by quantized depth
            std::map<int, std::vector<ClockInst*>> depthGroups;
            for (auto& fi : ffInfos) {
              depthGroups[fi.depth].push_back(fi.inst);
            }

            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) {
                // Direct connection to driver subnet
                for (auto* ff : ffs) {
                  driverSubNet->addInst(*ff);
                  ++branchDirectFFs;
                }
              } else {
                // Determine group tier from majority vote of member FFs
                int grpBottom = 0, grpUpper = 0;
                for (auto* ff : ffs) {
                  int ft = cts3dDb_->getInstTier(ff->getDbInst());
                  if (ft == 0) ++grpBottom; else ++grpUpper;
                }
                const int grpTier = (grpBottom > grpUpper) ? 0 : 1;

                // Filter out cross-tier minority FFs before building the chain.
                // Minority FFs with a different tier than grpTier are direct-connected
                // to driverSubNet (the leaf net) instead of the group chain output net.
                // Rationale: a small chain-output net spanning two tiers triggers
                // GRT-0183 (heap underflow) because local HB via capacity is exhausted.
                // The leaf net (larger, many sinks spread across the die) handles
                // cross-tier routing flexibly, exactly as balanced CTS already does.
                // No extra chain buffers are added for minority FFs → buffer count
                // stays the same, keeping the LP problem efficient.
                std::vector<ClockInst*> sameTierFFs;
                for (auto* ff : ffs) {
                  int ft = cts3dDb_->getInstTier(ff->getDbInst());
                  if (ft < 0) ft = grpTier;  // unknown tier → assume same as group
                  if (ft != grpTier) {
                    // Cross-tier minority FF: return to leaf net (original placement)
                    driverSubNet->addInst(*ff);
                    ++branchDirectFFs;
                  } else {
                    sameTierFFs.push_back(ff);
                  }
                }

                if (sameTierFFs.empty()) {
                  // All FFs in this depth bucket were cross-tier; nothing to chain.
                  totalBranchSinks_debug += ffs.size();
                  continue;
                }

                // Compute centroid of same-tier FFs only (accurate chain placement)
                long long sumX = 0, sumY = 0;
                for (auto* ff : sameTierFFs) {
                  sumX += ff->getX();
                  sumY += ff->getY();
                }
                const int cx = static_cast<int>(sumX / sameTierFFs.size());
                const int cy = static_cast<int>(sumY / sameTierFFs.size());
                const std::string grpBufMaster = cts3dDb_->getBufferForTier(
                    options_->getRootBuffer(), grpTier);

                // Build shared chain: driver → grp_buf_0 → ... → grp_buf_{depth-1} → FFs
                // All chain buffers and sink FFs are the same tier → no cross-tier
                // chain nets → no GRT-0183.
                // V59: Sequential fixed-spacing placement along driver→centroid vector.
                // minStep = bufWidth → buffer bodies cannot overlap.
                // Replaces equal-fraction interpolation which caused overlap when
                // span < depth × bufWidth (DRT unroutable, GUI-0070).
                int grpBufW = 0;
                {
                  odb::dbMaster* m = db_->findMaster(grpBufMaster.c_str());
                  if (m) grpBufW = m->getWidth();
                }
                int grpSpanX = std::abs(cx - driverInst->getX());
                int grpSpanY = std::abs(cy - driverInst->getY());
                int grpSpan = std::max(grpSpanX, grpSpanY);
                int grpMinStep = std::max(grpBufW, 1);
                double grpFracStep = (grpSpan > 0)
                    ? static_cast<double>(grpMinStep) / grpSpan : 1.0;

                ClockSubNet* cur = driverSubNet;
                int prevGrpBufX = driverInst->getX();
                for (int d = 0; d < depth; ++d) {
                  double frac = grpFracStep * (d + 1);
                  if (frac >= 1.0) break;  // no room → stop chain
                  int bx = driverInst->getX()
                      + static_cast<int>((cx - driverInst->getX()) * frac);
                  int by = driverInst->getY()
                      + static_cast<int>((cy - driverInst->getY()) * frac);
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, grpBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }
                  // Post-snap X-only overlap check: legalizeOneBuffer can snap to
                  // a different ROW, breaking the pre-snap spacing guarantee.
                  if (grpBufW > 0 && std::abs(bx - prevGrpBufX) <= grpBufW) {
                    break;  // overlap → stop chain
                  }
                  prevGrpBufX = bx;

                  // use myGrpId (global counter) instead of idx
                  // to avoid collision across sub-trees.
                  const std::string bName = tierPrefix_ + "clkbuf_grp_"
                      + std::to_string(myGrpId) + "_c"
                      + std::to_string(numSinks) + "_d"
                      + std::to_string(depth) + "_"
                      + std::to_string(d);
                  ClockInst& grpBuf = clock_.addClockBuffer(
                      bName, grpBufMaster, bx, by);
                  cts3dDb_->setClockInstTier(grpBuf, grpTier);
                  addTreeLevelBuffer(&grpBuf);

                  cur->addInst(grpBuf);
                  ClockSubNet& next = clock_.addSubNet(
                      tierPrefix_ + "clknet_grp_" + std::to_string(myGrpId) + "_c"
                      + std::to_string(numSinks) + "_d"
                      + std::to_string(depth) + "_"
                      + std::to_string(d));
                  next.addInst(grpBuf);
                  cur = &next;
                  ++branchGroupBufs;
                }

                // Connect same-tier FFs to final subnet (guaranteed single-tier net)
                cur->setLeafLevel(true);
                for (auto* ff : sameTierFFs) {
                  cur->addInst(*ff);
                  ++branchGroupedFFs;
                }
              }

              // Debug stats
              totalLpDelay_debug += 0;  // grouped: tracked at group level
              totalBranchSinks_debug += ffs.size();
            }
          };

          // GH-Tree Cascaded Delay Tap.
          // Replaces parallel grouped delay chains with cascaded tap structure.
          // Parallel chains (V49/V50):
          //   leaf -> [chain_d1]           -> d1_FFs    (1 buf, 1 net)
          //   leaf -> [buf -> buf]         -> d2_FFs    (2 bufs, 2 nets)
          //   leaf -> [buf -> buf -> buf]  -> d3_FFs    (3 bufs, 3 nets)
          //   Total: 6 bufs, 6 nets for maxDepth=3
          //
          // Cascaded tap (V52 GH-Tree):
          //   leaf -> tap1 -> tap2 -> tap3
          //             |       |       |
          //          d1_FFs  d2_FFs  d3_FFs
          //   Total: 3 bufs, 3 nets for maxDepth=3  (50% reduction)
          //
          // Same depth delivery: depth=k FF passes through exactly k tap buffers.
          // d0_FFs connect directly to driverSubNet (no tap buffer).
          // Same-tier purification (V49a) applied before cascaded chain.
          // Tier-aware delivery (CTS_ENABLE_TIER_DELIVER=1): instead of demoting
          // cross-tier (minority) FFs to depth=0 (no useful-skew delivery), build a
          // SEPARATE same-tier cascaded-TAP chain for them via recursion, delivering
          // skew to BOTH tiers of a mixed cluster. Only the single chain-input net is
          // cross-tier (not per-FF), avoiding the GRT-0183 heap underflow the original
          // per-FF cross-tier relay nets caused (V49a). Fixes the measured
          // heterogeneous-ng45 realization collapse (~50% of targeted FFs were purified
          // to depth=0 -> delivery ratio 0.13 vs LP-achievable 90.7%).
          const bool enableTierDeliver =
              (std::getenv("CTS_ENABLE_TIER_DELIVER")
               && std::atoi(std::getenv("CTS_ENABLE_TIER_DELIVER")) != 0);
          // std::function (not auto) so the lambda can recurse on the minority tier.
          std::function<void(ClockInst*, ClockSubNet*, std::vector<FfGroupInfo>&)>
              buildGHSubTree;
          buildGHSubTree = [&](ClockInst* driverInst,
                                     ClockSubNet* driverSubNet,
                                     std::vector<FfGroupInfo>& ffInfos) {
            // Group FFs by quantized depth
            std::map<int, std::vector<ClockInst*>> depthGroups;
            std::vector<FfGroupInfo> minorityInfos;  // cross-tier FFs -> recurse
            for (auto& fi : ffInfos) {
              depthGroups[fi.depth].push_back(fi.inst);
            }

            // depth=0 FFs: direct connection to driver subnet
            if (depthGroups.count(0)) {
              for (auto* ff : depthGroups[0]) {
                driverSubNet->addInst(*ff);
                ++branchDirectFFs;
              }
            }

            // Find max depth with non-empty FFs
            int maxDepth = 0;
            for (auto& [d, ffs] : depthGroups) {
              if (d > 0 && !ffs.empty()) maxDepth = d;
            }
            if (maxDepth == 0) return;  // All depth=0, no sub-tree needed

            // Determine chain tier from ALL depth>0 FFs (majority vote)
            int chainBottom = 0, chainUpper = 0;
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              for (auto* ff : ffs) {
                int ft = cts3dDb_->getInstTier(ff->getDbInst());
                if (ft == 0) ++chainBottom; else ++chainUpper;
              }
            }
            const int chainTier = (chainBottom > chainUpper) ? 0 : 1;

            // same-tier purification: cross-tier minority FFs.
            // Legacy (enableTierDeliver=false): demote to depth=0 (NO delivery).
            // Tier-deliver: collect into minorityInfos (preserving depth) to build a
            // separate same-tier chain via recursion after the majority chain.
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              auto it = ffs.begin();
              while (it != ffs.end()) {
                int ft = cts3dDb_->getInstTier((*it)->getDbInst());
                if (ft < 0) ft = chainTier;  // unknown -> assume same
                if (ft != chainTier) {
                  if (enableTierDeliver) {
                    minorityInfos.push_back(FfGroupInfo{*it, 0.0, depth});
                  } else {
                    driverSubNet->addInst(**it);
                    ++branchDirectFFs;
                  }
                  it = ffs.erase(it);
                } else {
                  ++it;
                }
              }
            }

            // Tier-aware delivery: build the minority-tier chain NOW (before the
            // majority chain's possible early-return below), so BOTH tiers of a
            // mixed cluster get useful-skew delivery. The recursion sees a
            // single-tier FF set -> no further minority -> terminates after 1 level.
            if (enableTierDeliver && !minorityInfos.empty()) {
              buildGHSubTree(driverInst, driverSubNet, minorityInfos);
            }

            // Recompute maxDepth after purification
            maxDepth = 0;
            for (auto& [d, ffs] : depthGroups) {
              if (d > 0 && !ffs.empty()) maxDepth = d;
            }
            if (maxDepth == 0) return;

            // Cost-aware adaptive depth selection
            // Scan from maxDepth downward. For each depth gap (empty taps),
            // check if the FFs beyond the gap justify the buffer cost.
            // If cost > benefit, reassign those FFs to the lower depth.
            // Proposition: reassign when ff_count * buf_delay * gap < gap * costRatio
            //   i.e., ff_count * buf_delay < costRatio  (per-gap-tap basis)
            {
              // Collect occupied depths in ascending order
              std::vector<int> occupiedDepths;
              for (auto& [d, ffs] : depthGroups) {
                if (d > 0 && !ffs.empty()) {
                  occupiedDepths.push_back(d);
                }
              }
              // Scan from highest depth downward
              bool pruned = true;
              while (pruned && occupiedDepths.size() > 1) {
                pruned = false;
                int highD = occupiedDepths.back();
                int lowD  = occupiedDepths[occupiedDepths.size() - 2];
                int gap   = highD - lowD;  // empty taps between lowD and highD
                int ffCount = static_cast<int>(depthGroups[highD].size());

                // Cost:    gap empty tap buffers (one buf+net per empty depth)
                // Benefit: ffCount FFs get (highD - lowD) * buf_delay more delivery
                // Reassign if ffCount is too small to justify the gap cost.
                // Threshold: ffCount * perFfBufDelay_ (ns) < adaptiveDepthCostRatio_ * gap_cost
                // Simplified: ffCount < costRatio * gap  (1 FF per empty tap)
                bool shouldReassign = false;
                if (gap > 1) {
                  // Gap has (gap-1) empty taps. Reassign if too few FFs.
                  int emptyTaps = gap - 1;
                  shouldReassign = (ffCount <= emptyTaps * adaptiveDepthCostRatio_);
                }
                // Also reassign if only 1 FF at the highest depth and gap >= 2
                if (ffCount == 1 && gap >= 2) {
                  shouldReassign = true;
                }

                if (shouldReassign) {
                  // Move highD FFs to lowD
                  for (auto* ff : depthGroups[highD]) {
                    depthGroups[lowD].push_back(ff);
                  }
                  debugPrint(logger_, CTS, "HTree", 1,
                      "Adaptive depth: reassign {} FFs depth {} -> {} "
                      "(gap={}, saved {} empty taps)",
                      ffCount, highD, lowD, gap, gap - 1);
                  depthGroups.erase(highD);
                  occupiedDepths.pop_back();
                  pruned = true;
                }
              }
              // Recompute maxDepth after adaptive pruning
              maxDepth = 0;
              for (auto& [d, ffs] : depthGroups) {
                if (d > 0 && !ffs.empty()) maxDepth = d;
              }
              if (maxDepth == 0) return;
            }

            // Compute centroid of ALL depth>0 same-tier FFs
            long long sumX = 0, sumY = 0;
            int nDownstream = 0;
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              for (auto* ff : ffs) {
                sumX += ff->getX();
                sumY += ff->getY();
                ++nDownstream;
              }
            }
            if (nDownstream == 0) return;
            const int cx = static_cast<int>(sumX / nDownstream);
            const int cy = static_cast<int>(sumY / nDownstream);
            const std::string chainBufMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), chainTier);

            // Cascaded delay tap chain: driver -> tap1 -> tap2 -> ... -> tap_maxDepth.
            // Each series tap adds one buffer delay; depth=d FFs tap off stage d.
            int tapBufW = 0;
            {
              odb::dbMaster* m = db_->findMaster(chainBufMaster.c_str());
              if (m) tapBufW = m->getWidth();
            }
            // Separation between consecutive nominal tap sites: >= 1.5x buffer width.
            const int tapStep = std::max(tapBufW + tapBufW / 2, wireSegmentUnit_);
            double dirX = static_cast<double>(cx - driverInst->getX());
            double dirY = static_cast<double>(cy - driverInst->getY());
            const double dirLen = std::sqrt(dirX * dirX + dirY * dirY);
            double ux, uy;
            if (dirLen < static_cast<double>(tapStep)) {
              ux = 1.0;  // centroid ~= driver: step into nearby free space (+x)
              uy = 0.0;
            } else {
              ux = dirX / dirLen;
              uy = dirY / dirLen;
            }

            // Decoupled cascaded TAP delivery (V59_Deliver):
            // Useful-skew delay is realized by the COUNT of series tap buffers
            // (each adds ~bufDelay), NOT by physical span. The original code placed
            // taps along the driver->centroid vector and required physical room; for
            // tight clusters the depth>0 centroid coincides with the leaf buffer
            // (span~=0) so no tap fit and the LP skew was never delivered (those FFs
            // were even left disconnected). Here every tap is placed at a free legal
            // site offset from the driver — legalizeOneBuffer()/resolveLocationCollision()
            // guarantee distinct, non-overlapping, on-grid sites — so the series chain
            // always forms regardless of cluster tightness.
            ClockSubNet* curNet = driverSubNet;
            for (int d = 1; d <= maxDepth; ++d) {
              int bx = driverInst->getX()
                  + static_cast<int>(std::lround(ux * tapStep * d));
              int by = driverInst->getY()
                  + static_cast<int>(std::lround(uy * tapStep * d));

              // Snap to a free legal site (collision-resolved -> distinct sites).
              {
                Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                     static_cast<double>(by) / wireSegmentUnit_);
                Point<double> legalLoc = legalizeOneBuffer(rawLoc, chainBufMaster);
                bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
              }

              // Clamp the tap to the core area so it can never land off-die.
              // TreeBuilder::legalizeOneBuffer only snaps around blockages (not the
              // die boundary), and resolveLocationCollision jitters +X for a unique
              // map key WITHOUT any bound check — so a near-edge driver stepping
              // toward the die edge (dirLen<tapStep uses ux=+1) can place a cascaded
              // tap outside the die → GRT-0035 "outside die area" → GRT SIGSEGV in
              // routing. Keep the buffer origin+extent inside the core; DPL then
              // re-legalizes onto the row grid.
              if (db_ != nullptr && db_->getChip() != nullptr) {
                const odb::Rect core
                    = db_->getChip()->getBlock()->getCoreArea();
                odb::dbMaster* cm = db_->findMaster(chainBufMaster.c_str());
                const int bw = (cm != nullptr) ? cm->getWidth() : 0;
                const int bh = (cm != nullptr) ? cm->getHeight() : 0;
                bx = std::max<int>(core.xMin(),
                                   std::min<int>(bx, core.xMax() - bw));
                by = std::max<int>(core.yMin(),
                                   std::min<int>(by, core.yMax() - bh));
              }

              // Create tap buffer. When tier-deliver is on, tag the name with the
              // chain tier so the majority and minority chains (same myGrpId/numSinks)
              // don't collide. Empty tag when off -> byte-identical to baseline.
              const std::string bName = tierPrefix_ + "clkbuf_ghtap_"
                  + std::to_string(myGrpId) + "_c"
                  + std::to_string(numSinks)
                  + (enableTierDeliver ? ("_tr" + std::to_string(chainTier))
                                       : std::string())
                  + "_t"
                  + std::to_string(d);
              ClockInst& tapBuf = clock_.addClockBuffer(
                  bName, chainBufMaster, bx, by);
              cts3dDb_->setClockInstTier(tapBuf, chainTier);
              addTreeLevelBuffer(&tapBuf);
              curNet->addInst(tapBuf);

              // Create tap output net (same tier-tag scheme as the buffer above).
              ClockSubNet& tapNet = clock_.addSubNet(
                  tierPrefix_ + "clknet_ghtap_" + std::to_string(myGrpId) + "_c"
                  + std::to_string(numSinks)
                  + (enableTierDeliver ? ("_tr" + std::to_string(chainTier))
                                       : std::string())
                  + "_t"
                  + std::to_string(d));
              tapNet.addInst(tapBuf);

              // Connect depth=d FFs to this tap net (empty intermediate depths are
              // pass-through delay stages needed to reach deeper FFs).
              if (depthGroups.count(d) && !depthGroups[d].empty()) {
                for (auto* ff : depthGroups[d]) {
                  tapNet.addInst(*ff);
                  ++branchGroupedFFs;
                }
              }

              // Mark final tap as leaf level
              if (d == maxDepth) {
                tapNet.setLeafLevel(true);
              }

              curNet = &tapNet;
              ++branchGroupBufs;
            }
            logger_->info(CTS, 581,
                "Cascaded TAP (decoupled): {} series taps for {} depth>0 FFs "
                "(span-independent)",
                maxDepth, nDownstream);

            totalBranchSinks_debug += nDownstream;
          };

          // ---------------- Shared-TAP backbone (CTS_SHARED_TAP=1) ----------------
          // Per-leaf cascaded chains rebuild the whole delay ladder for every leaf
          // cluster, so the pass-through stages (taps that drive no FF at all) are
          // paid once per leaf.  Measured on the routed DEF: 45-61% of all tap
          // buffers serve ZERO FFs (aes 86/140, ibex 453/765, jpeg 399/886), and
          // TAP buffers outnumber H-tree buffers 1.2-3.0x — they dominate the
          // buffer / rWL / power / HBT overhead vs Pin3D.
          //
          // Shared TAP builds ONE ladder per (branch, tier) and lets every FF tap
          // in at its own rung, so the per-FF delivery granularity that the LP
          // solves for is preserved exactly.  Delay bookkeeping:
          //   per-leaf  : branch -> leafBuf -> tap_1..tap_d -> FF   = d+1 stages
          //   shared    : branch -> tap_1..tap_{d+1}        -> FF   = d+1 stages
          // so the ladder gets ONE EXTRA stage and a depth-d FF attaches at rung
          // d+1.  Realized delay per FF is unchanged; only the number of ladders is.
          //
          // Offline check on the real timing graph (aes, 7854 edges): collapsing
          // delivery to per-leaf granularity instead LOSES the useful skew (setup
          // -30,184 vs balanced -28,345), while shared-ladder + per-FF tap-in keeps
          // 100% of the gain (-25,387 vs per-leaf-chain -25,390) at 54-83% fewer taps.
          //
          // Only the CLUSTERED sink path is shared.  Individual (unclustered) FFs
          // hang straight off the branch net with a d-stage chain (no leaf buffer),
          // so their stage count differs and they keep their own chains.
          const bool sharedTap
              = (std::getenv("CTS_SHARED_TAP")
                 && std::atoi(std::getenv("CTS_SHARED_TAP")) != 0);
          // Sharing lengthens the tap->FF wires (asap7_3D/aes: mean tap-net length
          // 19.1um -> 68.8um, total 3477um -> 7503um).  That extra wire delay lands
          // on the delivered FFs' capture clock, so hold degrades on every design
          // even where setup improves.  CTS_SHARED_TAP_MAX_SPAN_UM > 0 bins the FFs
          // of a (branch, tier) into square cells of that size and builds one ladder
          // per cell, trading some of the buffer saving back for shorter tap nets.
          // 0 = one ladder per (branch, tier), i.e. maximum sharing.
          int sharedTapSpanUm = 0;
          if (const char* e = std::getenv("CTS_SHARED_TAP_MAX_SPAN_UM")) {
            sharedTapSpanUm = std::atoi(e);
          }
          // 1 = a leaf whose FFs are ALL depth>0 keeps a private ladder (V76 behaviour).
          // Empirically design-adaptive: helps ibex (72% such leaves), irrelevant on aes (1.6%).
          const bool sharedTapPrivateIfNoD0
              = (std::getenv("CTS_SHARED_TAP_PRIVATE_IF_NO_D0")
                 && std::atoi(std::getenv("CTS_SHARED_TAP_PRIVATE_IF_NO_D0")) != 0);
          // CROSS-TIER LADDER (V78).
          //
          // xtier=0 bins ladders by (tier, cell): one ladder per tier per branch, so the
          // structure is really "a 2D shared ladder on each die" — nothing about it needs
          // face-to-face stacking.
          //
          // xtier=1 lets a SINGLE ladder feed the FFs of BOTH tiers: bins are keyed by
          // (cell) only, the rung buffers live on the majority tier of the bin, and the
          // minority tier's FFs reach them over a hybrid-bond via.  That is the structure
          // an F2F stack makes possible and a 2D flow cannot.  Two measurements of ours
          // say it should be viable: the inter-tier hop is essentially free in timing
          // (hb via R 63 -> 6300 ohm leaves TNS byte-identical), and cross-tier clock
          // nets do route (asap7_3D/aes leaf nets already carry 73 hb vias).
          //
          // Risks to measure rather than assume away:
          //   * GRT-0183 (heap underflow) previously killed small cross-tier chain nets
          //   * more hb vias -> routing / power / wirelength cost
          //   * delivered delay may be less stable than on a single-tier ladder
          // So: aes + swerv first, routed to completion, before believing anything.
          const bool sharedTapXtier
              = (std::getenv("CTS_SHARED_TAP_XTIER")
                 && std::atoi(std::getenv("CTS_SHARED_TAP_XTIER")) != 0);
          int sharedTapLeafKeepDemotes = 0;
          std::vector<FfGroupInfo> branchInfos;              // depth>0 FFs, this branch
          std::map<ClockInst*, ClockSubNet*> ffHomeNet;      // FF -> its leaf subnet

          // Iterate sinks in this branch (same pattern as per-FF relay)
          const std::vector<Point<double>>& sinkLocs
              = leafTopology.getBranchSinksLocations(idx);
          for (const auto& loc : sinkLocs) {
            // Epsilon-tolerant sink lookup (V37 CTS-0080 fix)
            // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 0.01;
              auto hint = mapLocationToSink_.lower_bound(
                  Point<double>(loc.getX() - eps, loc.getY() - eps));
              for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                if (it->first.getX() > loc.getX() + eps) break;
                if (std::abs(it->first.getX() - loc.getX()) < eps
                    && std::abs(it->first.getY() - loc.getY()) < eps) {
                  sinkIt = it;
                  break;
                }
              }
            }
            if (sinkIt == mapLocationToSink_.end()) {
              // Second-level clustering or legalizeOneBuffer can shift a centroid
              // beyond eps, so the branch-sink location is no longer a map key.
              //
              // The old code skipped the sink, with the comment "connects to the
              // leaf net directly".  That is WRONG when the sink is a leaf-cluster
              // BUFFER: skipping means `subNet->addInst(*sinkInst)` below never
              // runs, so the buffer is emitted with a floating input.  STA then
              // sees no clock on everything it drives -> those FFs silently lose
              // their setup/hold checks (and vectorless power collapses).
              // Observed on asap7_nangate45_3D/jpeg: CTS-0419 reports 1/1193 sinks
              // missing and clkbuf_leaf_168_clk comes out dangling (CTS-0726).
              // Present since long before the invariant existed, i.e. that design's
              // timing has been quietly flattered.
              //
              // Fall back to the nearest map entry (unbounded), the same recovery
              // the CTS-0079 / CTS-0404 lookups already use.
              double bestDist = std::numeric_limits<double>::max();
              auto bestIt = mapLocationToSink_.end();
              for (auto it = mapLocationToSink_.begin();
                   it != mapLocationToSink_.end(); ++it) {
                const double d = std::abs(it->first.getX() - loc.getX())
                               + std::abs(it->first.getY() - loc.getY());
                if (d < bestDist) { bestDist = d; bestIt = it; }
              }
              if (bestIt == mapLocationToSink_.end()) {
                logger_->warn(CTS, 583,
                    "Branch sink at ({:.6f}, {:.6f}) not found and no nearest "
                    "match exists; it will be left unconnected.",
                    loc.getX(), loc.getY());
                continue;
              }
              logger_->warn(CTS, 584,
                  "Branch sink at ({:.6f}, {:.6f}) missing from the location map; "
                  "using nearest match at distance {:.6f} to avoid a dangling "
                  "clock buffer.",
                  loc.getX(), loc.getY(), bestDist);
              sinkIt = bestIt;
            }
            ClockInst* sinkInst = sinkIt->second;

            // Check if cluster buffer or individual FF
            auto cIt = clusterBufTarget_.find(sinkInst);
            if (cIt != clusterBufTarget_.end()) {
              // Cluster buffer: add to branch subnet, then group its member FFs
              subNet->addInst(*sinkInst);

              ClockSubNet* leafSubNet = nullptr;
              clock_.forEachSubNet([&](ClockSubNet& sn) {
                if (sn.getDriver() == sinkInst) {
                  leafSubNet = &sn;
                }
              });

              if (leafSubNet != nullptr) {
                std::vector<ClockInst*> memberFFs;
                leafSubNet->forEachSink([&](ClockInst* ff) {
                  memberFFs.push_back(ff);
                });
                // Remove all sinks (will re-add via group chains)
                std::set<ClockInst*> toRemove(
                    memberFFs.begin(), memberFFs.end());
                leafSubNet->removeSinks(toRemove);

                // Collect FF info and build group chains
                std::vector<FfGroupInfo> ffInfos;
                collectFFs(sinkInst, leafSubNet, ffInfos, memberFFs);
                applyHoldBudgetCap(ffInfos);

                // Shared TAP: hand the depth>0 FFs to the per-branch backbone
                // instead of building a private ladder for this leaf.
                //
                // The leaf subnet must keep at least one sink: a driver with no
                // sinks is an empty clock net (CTS-0885 disconnected-subnet path).
                // Two ways to satisfy that when every member FF wants depth>0:
                //
                //  (A) demote the SHALLOWEST member (smallest LP target -> least skew
                //      lost) to depth 0 and share the rest.  Default.
                //  (B) leave the whole leaf on a private ladder.
                //
                // V76 did (B) unconditionally.  That looked like a bug — on
                // asap7_3D/ibex it left 159 of 177 leaves unshared — but ibex turns out
                // to have 162/225 leaves with NO depth-0 FF (aes: 1/63), and post-route
                // ibex is BEST with (B): route setup -1,063,020 vs -1,086,114 for (A).
                // A leaf whose FFs are all delayed keeps its tap->FF wires short on a
                // private ladder; forcing it into the branch ladder stretches them
                // (mean tap net 37um -> 99um) and the wire delay eats the gain.
                // So (B) is a real, design-adaptive heuristic, not an accident:
                // CTS_SHARED_TAP_PRIVATE_IF_NO_D0=1 restores it.
                if (sharedTap && !ffInfos.empty()) {
                  bool hasDepth0 = false;
                  for (auto& fi : ffInfos) {
                    if (fi.depth <= 0) { hasDepth0 = true; break; }
                  }
                  if (!hasDepth0 && sharedTapPrivateIfNoD0) {
                    // (B): this leaf keeps its own ladder.
                    if (enableGhTree_) {
                      buildGHSubTree(sinkInst, leafSubNet, ffInfos);
                    } else {
                      buildGroupChains(sinkInst, leafSubNet, ffInfos);
                    }
                    ++numSinks;
                    continue;
                  }
                  if (!hasDepth0) {
                    auto shallowest = std::min_element(
                        ffInfos.begin(), ffInfos.end(),
                        [](const FfGroupInfo& a, const FfGroupInfo& b) {
                          return a.depth < b.depth;
                        });
                    shallowest->depth = 0;
                    ++sharedTapLeafKeepDemotes;
                  }
                  for (auto& fi : ffInfos) {
                    if (fi.depth <= 0) {
                      leafSubNet->addInst(*fi.inst);
                      ++branchDirectFFs;
                    } else {
                      branchInfos.push_back(fi);
                      ffHomeNet[fi.inst] = leafSubNet;
                    }
                  }
                } else if (enableGhTree_) {
                  // GH-Tree cascaded tap or legacy parallel chains
                  buildGHSubTree(sinkInst, leafSubNet, ffInfos);
                } else {
                  buildGroupChains(sinkInst, leafSubNet, ffInfos);
                }
              } else {
                ++branchDirectFFs;
              }
            } else {
              // Individual FF (no clustering): treat as single-FF group
              std::string ffName = getFFName(sinkInst);
              double ffTarget = cts3dDb_->getSkewTarget(ffName);
              // Debug: trace individual FF target lookup
              if (ffTarget > 0.01) {
                debugPrint(logger_, CTS, "HTree", 1,
                    "collectFFs individual: '{}' target={:.4f}ns",
                    ffName, ffTarget);
              }
              // globalMean fix: absolute LP target
              double ffDelta = ffTarget;
              int depth = 0;
              if (ffDelta > groupedDelayMinDelta_) {
                depth = static_cast<int>(
                    std::round(ffDelta / perFfBufDelay_));
                depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
              }
              // Removed redundant holdSkip (same as above).
              // V58 Fix C: Macro sinks (SRAM) get depth=0.
              odb::dbInst* dInst = sinkInst->getDbInst();
              if (isMacroBlockInst(dInst)) {
                depth = 0;
              }
              if (depth == 0) {
                subNet->addInst(*sinkInst);
                ++branchDirectFFs;
              } else {
                // Single FF needs its own chain
                int ft = cts3dDb_->getInstTier(sinkInst->getDbInst());
                const int ffTier = (ft >= 0) ? ft : branch_tier;
                const std::string ffBufMaster = cts3dDb_->getBufferForTier(
                    options_->getRootBuffer(), ffTier);
                // V59: Sequential fixed-spacing placement along driver→FF vector.
                // minStep = bufWidth → buffer bodies cannot overlap.
                int bufW = 0;
                {
                  odb::dbMaster* m = db_->findMaster(ffBufMaster.c_str());
                  if (m) bufW = m->getWidth();
                }
                int ffSpanX = std::abs(sinkInst->getX() - subNet->getDriver()->getX());
                int ffSpanY = std::abs(sinkInst->getY() - subNet->getDriver()->getY());
                int ffSpan = std::max(ffSpanX, ffSpanY);
                int ffMinStep = std::max(bufW, 1);
                double ffFracStep = (ffSpan > 0)
                    ? static_cast<double>(ffMinStep) / ffSpan : 1.0;

                if (ffFracStep >= 1.0 || depth == 0) {
                  // No room for any buffer → connect directly
                  subNet->addInst(*sinkInst);
                  ++branchDirectFFs;
                  ++numSinks;
                  continue;
                }
                ClockSubNet* cur = subNet;
                int prevFfBufX = subNet->getDriver()->getX();
                for (int d = 0; d < depth; ++d) {
                  double frac = ffFracStep * (d + 1);
                  if (frac >= 1.0) {
                    // No room → connect FF to current net and stop
                    cur->addInst(*sinkInst);
                    cur->setLeafLevel(true);
                    depth = d;
                    break;
                  }
                  int bx = subNet->getDriver()->getX()
                      + static_cast<int>(
                          (sinkInst->getX() - subNet->getDriver()->getX()) * frac);
                  int by = subNet->getDriver()->getY()
                      + static_cast<int>(
                          (sinkInst->getY() - subNet->getDriver()->getY()) * frac);
                  // Grid-snap per-FF grouped delay buffer
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, ffBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }
                  // Post-snap X-only overlap check
                  if (bufW > 0 && std::abs(bx - prevFfBufX) <= bufW) {
                    cur->addInst(*sinkInst);
                    cur->setLeafLevel(true);
                    depth = d;
                    break;
                  }
                  prevFfBufX = bx;
                  const std::string bName = tierPrefix_ + "clkbuf_grp_"
                      + std::to_string(myGrpId) + "_s"
                      + std::to_string(numSinks) + "_"
                      + std::to_string(d);
                  ClockInst& grpBuf = clock_.addClockBuffer(
                      bName, ffBufMaster, bx, by);
                  cts3dDb_->setClockInstTier(grpBuf, ffTier);
                  addTreeLevelBuffer(&grpBuf);
                  cur->addInst(grpBuf);
                  ClockSubNet& next = clock_.addSubNet(
                      tierPrefix_ + "clknet_grp_" + std::to_string(myGrpId) + "_s"
                      + std::to_string(numSinks) + "_"
                      + std::to_string(d));
                  next.addInst(grpBuf);
                  cur = &next;
                  ++branchGroupBufs;
                }
                cur->setLeafLevel(true);
                cur->addInst(*sinkInst);
                ++branchGroupedFFs;
              }
            }
            ++numSinks;
          }

          // ---- Build the one shared TAP ladder for this branch (CTS_SHARED_TAP) ----
          if (sharedTap && !branchInfos.empty()) {
            if (sharedTapLeafKeepDemotes > 0) {
              debugPrint(logger_, CTS, "HTree", 1,
                  "Shared TAP branch {}: demoted {} shallowest FFs so their leaf "
                  "nets keep a sink",
                  idx, sharedTapLeafKeepDemotes);
            }
            // Bin the FFs by (tier, spatial cell).  One ladder per bin.
            // FFs whose tier is unknown cannot go on a single-tier ladder; leave
            // them on their leaf net (no delivery, hold-safe).
            const int shDbuPerUm
                = (db_ != nullptr && db_->getChip() != nullptr)
                    ? db_->getChip()->getBlock()->getDbUnitsPerMicron()
                    : 1000;
            const long long shCell = (sharedTapSpanUm > 0)
                ? static_cast<long long>(sharedTapSpanUm) * shDbuPerUm
                : 0;
            std::map<std::tuple<int, long long, long long>,
                     std::vector<FfGroupInfo>> shBins;
            for (auto& fi : branchInfos) {
              const int t = cts3dDb_->getInstTier(fi.inst->getDbInst());
              if (t < 0) {
                ffHomeNet[fi.inst]->addInst(*fi.inst);
                ++branchDirectFFs;
                continue;
              }
              const long long gx = shCell
                  ? static_cast<long long>(std::floor(
                        static_cast<double>(fi.inst->getX()) / shCell))
                  : 0;
              const long long gy = shCell
                  ? static_cast<long long>(std::floor(
                        static_cast<double>(fi.inst->getY()) / shCell))
                  : 0;
              // xtier: drop tier from the bin key so one ladder serves both tiers.
              shBins[std::make_tuple(sharedTapXtier ? 0 : t, gx, gy)].push_back(fi);
            }
          int shLadderIdx = 0;
          for (auto& [shKey, shBinInfos] : shBins) {
            const int shTierIter = std::get<0>(shKey);
            const int shLadder = shLadderIdx++;
            std::map<int, std::vector<ClockInst*>> dg;
            for (auto& fi : shBinInfos) {
              dg[fi.depth].push_back(fi.inst);
            }
            int shMaxDepth = 0;
            for (auto& [d, v] : dg) {
              if (d > 0 && !v.empty()) shMaxDepth = d;
            }

            // ONE LADDER PER (TIER, SPATIAL BIN), not one per branch.
            //
            // A chain-output net spanning two tiers exhausts local HB via capacity
            // and trips GRT-0183, so a ladder must stay single-tier.  A first cut
            // used a branch-wide majority vote and demoted the minority FFs back to
            // their leaf net (no useful skew).  Measured on asap7_3D/aes that cost
            // 28 of 123 delivered FFs and ~2% of the setup gain, because "minority"
            // at BRANCH granularity is far more FFs than at LEAF granularity.
            // xtier=0: bins are already single-tier; the filter is a cheap safety net.
            // xtier=1: the bin holds both tiers.  Put the rung buffers on the majority
            //          tier and let the minority tier's FFs cross the hybrid bond —
            //          no FF is dropped, so delivery coverage is unchanged.
            int shTier = shTierIter;
            if (sharedTapXtier) {
              int nb = 0, nu = 0;
              for (auto& [d, v] : dg) {
                if (d == 0) continue;
                for (auto* ff : v) {
                  (cts3dDb_->getInstTier(ff->getDbInst()) == 0) ? ++nb : ++nu;
                }
              }
              shTier = (nb > nu) ? 0 : 1;
            } else {
              for (auto& [d, v] : dg) {
                if (d == 0) continue;
                auto it = v.begin();
                while (it != v.end()) {
                  if (cts3dDb_->getInstTier((*it)->getDbInst()) != shTier) {
                    it = v.erase(it);
                  } else {
                    ++it;
                  }
                }
              }
            }
            int shDownstream = 0;
            for (auto& [d, v] : dg) {
              if (d > 0) shDownstream += static_cast<int>(v.size());
            }
            // Recompute after the tier filter — the deepest rung may have emptied.
            shMaxDepth = 0;
            for (auto& [d, v] : dg) {
              if (d > 0 && !v.empty()) shMaxDepth = d;
            }

            if (shDownstream > 0 && shMaxDepth > 0) {
              const std::string shBufMaster = cts3dDb_->getBufferForTier(
                  options_->getRootBuffer(), shTier);
              int shBufW = 0;
              if (odb::dbMaster* m = db_->findMaster(shBufMaster.c_str())) {
                shBufW = m->getWidth();
              }
              const int shStep = std::max(shBufW + shBufW / 2, wireSegmentUnit_);

              ClockInst* shDrv = subNet->getDriver();
              // Step the ladder from the branch driver toward the centroid of the
              // FFs it will feed, so the tap->FF wires stay short.
              double ccx = 0.0, ccy = 0.0;
              for (auto& [d, v] : dg) {
                if (d == 0) continue;
                for (auto* ff : v) { ccx += ff->getX(); ccy += ff->getY(); }
              }
              ccx /= shDownstream;
              ccy /= shDownstream;

              // --- V89 bridge prototype (CTS_BRIDGE_PROTOTYPE=1) ------------------
              // When this ladder sits on a tier different from the branch driver's
              // tier, the remote-ladder root crossing is otherwise scattered across
              // the branch net (router-placed, per remote sink). Insert a DEDICATED
              // two-pin cross-tier bridge: bridge_src on the driver tier + bridge_dst
              // on the ladder tier, both at the nearest legal site to the remote-FF
              // centroid h, joined by their own 2-pin net. That localizes the ONE F2F
              // crossing; bridge_dst then drives a single-tier ladder. Exact
              // co-location is NOT required — nearest legal sites around h suffice
              // (the 2-pin net still crosses once). Gated + optionally limited to one
              // branch (CTS_BRIDGE_BRANCH) so it can be validated in isolation first.
              const int driverTier = shDrv->getTier();
              const bool bridgeProto
                  = (std::getenv("CTS_BRIDGE_PROTOTYPE")
                     && std::atoi(std::getenv("CTS_BRIDGE_PROTOTYPE")) != 0);
              const int bridgeBranch = std::getenv("CTS_BRIDGE_BRANCH")
                  ? std::atoi(std::getenv("CTS_BRIDGE_BRANCH")) : -1;
              ClockSubNet* shCur = subNet;
              bool bridged = false;
              if (bridgeProto && driverTier >= 0 && shTier != driverTier
                  && (bridgeBranch < 0 || bridgeBranch == static_cast<int>(idx))) {
                const std::string bTag = "_sb" + std::to_string(idx) + "g"
                    + std::to_string(shLadder) + "_tr" + std::to_string(shTier);
                // bridge_src on the driver tier, at the legal site nearest h.
                const std::string bsrcMaster = cts3dDb_->getBufferForTier(
                    options_->getRootBuffer(), driverTier);
                Point<double> hSrc = legalizeOneBuffer(
                    Point<double>(ccx / wireSegmentUnit_, ccy / wireSegmentUnit_),
                    bsrcMaster);
                const int bsx = static_cast<int>(hSrc.getX() * wireSegmentUnit_);
                const int bsy = static_cast<int>(hSrc.getY() * wireSegmentUnit_);
                ClockInst& bsrc = clock_.addClockBuffer(
                    tierPrefix_ + "clkbuf_bridge_src" + bTag, bsrcMaster, bsx, bsy);
                cts3dDb_->setClockInstTier(bsrc, driverTier);
                addTreeLevelBuffer(&bsrc);
                subNet->addInst(bsrc);       // same-tier sink of the branch net
                ClockSubNet& bnet = clock_.addSubNet(
                    tierPrefix_ + "clknet_bridge" + bTag);
                bnet.addInst(bsrc);          // bridge_src drives the dedicated 2-pin net
                // bridge_dst on the ladder tier, at the legal site nearest h.
                Point<double> hDst = legalizeOneBuffer(
                    Point<double>(ccx / wireSegmentUnit_, ccy / wireSegmentUnit_),
                    shBufMaster);
                const int bdx = static_cast<int>(hDst.getX() * wireSegmentUnit_);
                const int bdy = static_cast<int>(hDst.getY() * wireSegmentUnit_);
                ClockInst& bdst = clock_.addClockBuffer(
                    tierPrefix_ + "clkbuf_bridge_dst" + bTag, shBufMaster, bdx, bdy);
                cts3dDb_->setClockInstTier(bdst, shTier);
                addTreeLevelBuffer(&bdst);
                bnet.addInst(bdst);          // sole sink of the 2-pin net = THE crossing
                ClockSubNet& bdstNet = clock_.addSubNet(
                    tierPrefix_ + "clknet_bridgeout" + bTag);
                bdstNet.addInst(bdst);       // bridge_dst drives the single-tier ladder
                // Depth compensation (CTS_BRIDGE_DEPTH_COMP, default ON): the bridge
                // already supplies 2 buffer stages (bridge_src + bridge_dst), so a
                // remote depth-d FF must tap only (d-2) LADDER rungs to stay iso-latency
                // with a source-tier depth-d FF (leafBuf + d taps). Without this the
                // whole remote subtree is shifted +2 stages -> tier-skew -> setup/hold
                // loss (the V89 Gate-3 regression). depth-1 FFs are delivered by the
                // bridge alone (bdst = exactly 2 stages, no rung); depth-0 FFs never
                // enter the ladder (they stay on the leaf net) so there is NO
                // over-delivery floor for depth>=1.
                const bool bridgeDepthComp
                    = !(std::getenv("CTS_BRIDGE_DEPTH_COMP")
                        && std::atoi(std::getenv("CTS_BRIDGE_DEPTH_COMP")) == 0);
                if (bridgeDepthComp) {
                  bridged = true;
                  auto it1 = dg.find(1);          // depth-1 FFs: bridge alone (2 stages)
                  if (it1 != dg.end()) {
                    for (auto* ff : it1->second) {
                      bdstNet.addInst(*ff);
                      ++branchGroupedFFs;
                    }
                  }
                }
                shDrv = &bdst;               // rungs now march from the dest-tier bridge
                shCur = &bdstNet;            // rung 1 attaches here (single-tier)
                logger_->info(CTS, 586,
                    "V89 bridge: branch {} ladder {} tier {}<-{} dedicated 2-pin net, "
                    "src@({},{}) dst@({},{})",
                    idx, shLadder, shTier, driverTier, bsx, bsy, bdx, bdy);
              }

              double sdx = ccx - shDrv->getX(), sdy = ccy - shDrv->getY();
              const double sdLen = std::sqrt(sdx * sdx + sdy * sdy);
              double sux = 1.0, suy = 0.0;
              if (sdLen > 1e-9) { sux = sdx / sdLen; suy = sdy / sdLen; }

              // ONE extra stage: a depth-d FF taps rung d+1, matching the
              // leafBuf + d taps of the per-leaf chain.  With a compensated bridge the
              // bridge's 2 stages replace the leafBuf-equiv + first rung, so a remote
              // depth-d FF taps ladder rung (d-1); the ladder needs (shMaxDepth-1) rungs.
              const int shStages = bridged ? std::max(0, shMaxDepth - 1) : shMaxDepth + 1;
              for (int s = 1; s <= shStages; ++s) {
                int bx = shDrv->getX()
                    + static_cast<int>(std::lround(sux * shStep * s));
                int by = shDrv->getY()
                    + static_cast<int>(std::lround(suy * shStep * s));
                {
                  Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                       static_cast<double>(by) / wireSegmentUnit_);
                  Point<double> legalLoc = legalizeOneBuffer(rawLoc, shBufMaster);
                  bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                  by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                }
                // Same off-die clamp as the per-leaf ladder: legalizeOneBuffer only
                // dodges blockages and resolveLocationCollision jitters +X without a
                // bound check, so a near-edge driver can push a tap outside the die
                // (GRT-0035 -> GRT SIGSEGV).
                if (db_ != nullptr && db_->getChip() != nullptr) {
                  const odb::Rect core = db_->getChip()->getBlock()->getCoreArea();
                  odb::dbMaster* cm = db_->findMaster(shBufMaster.c_str());
                  const int bw = (cm != nullptr) ? cm->getWidth() : 0;
                  const int bh = (cm != nullptr) ? cm->getHeight() : 0;
                  bx = std::max<int>(core.xMin(), std::min<int>(bx, core.xMax() - bw));
                  by = std::max<int>(core.yMin(), std::min<int>(by, core.yMax() - bh));
                }

                // shLadder disambiguates the spatial bins of the same (branch, tier)
                // when CTS_SHARED_TAP_MAX_SPAN_UM > 0.
                const std::string shTag = "_sb" + std::to_string(idx) + "g"
                    + std::to_string(shLadder) + "_c"
                    + std::to_string(numSinks) + "_tr" + std::to_string(shTier)
                    + "_t" + std::to_string(s);
                ClockInst& shBuf = clock_.addClockBuffer(
                    tierPrefix_ + "clkbuf_ghtap" + shTag, shBufMaster, bx, by);
                cts3dDb_->setClockInstTier(shBuf, shTier);
                addTreeLevelBuffer(&shBuf);
                shCur->addInst(shBuf);

                ClockSubNet& shNet = clock_.addSubNet(
                    tierPrefix_ + "clknet_ghtap" + shTag);
                shNet.addInst(shBuf);

                // Rung s carries the depth-(s-1) FFs.  Rung 1 (depth 0) is always
                // empty: those FFs stay on their leaf net.
                // Compensated bridge: the 2 bridge stages already cover leafBuf-equiv +
                // rung 1, so ladder rung s carries depth-(s+1) FFs (depth 1 was already
                // attached to bdstNet). This keeps remote FFs iso-latency with source.
                const int dHere = bridged ? s + 1 : s - 1;
                if (dHere >= 1 && dg.count(dHere)) {
                  for (auto* ff : dg[dHere]) {
                    shNet.addInst(*ff);
                    ++branchGroupedFFs;
                  }
                }
                if (s == shStages) {
                  shNet.setLeafLevel(true);
                }
                shCur = &shNet;
                ++branchGroupBufs;
              }
              // #leaves feeding THIS ladder = distinct leaf subnets of its FFs.
              std::set<ClockSubNet*> shLeaves;
              for (auto& [d, v] : dg) {
                if (d == 0) continue;
                for (auto* ff : v) { shLeaves.insert(ffHomeNet[ff]); }
              }
              logger_->info(CTS, 582,
                  "Shared TAP: branch {} ladder {} tier {} - {} rungs feed {} "
                  "depth>0 FFs from {} leaves",
                  idx, shLadder, shTier, shStages, shDownstream,
                  static_cast<int>(shLeaves.size()));
              totalBranchSinks_debug += shDownstream;
            }
          }  // per-tier ladder loop
          }

          numDelayBufs += branchGroupBufs;
          ++totalBranches_debug;
          if (branchGroupedFFs > 0) ++firedBranches_debug;
          totalDelivered_debug += branchGroupBufs * perFfBufDelay_;
          debugPrint(logger_, CTS, "HTree", 1,
              "{} branch {:2d}: grouped={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns",
              enableGhTree_ ? "ghtree" : "grouped", idx,
              branchGroupedFFs, branchGroupBufs,
              branchDirectFFs, branchHoldSkip, absTgt);
          return;  // skip legacy per-branch N-buffer code
        }

        // Per-FF relay buffer path (rewrote V39).
        // BUG #1 fix: V39 operated at cluster-buffer granularity (mean target),
        //   now iterates member FFs individually for per-FF relay decisions.
        // BUG #2 fix: V39 placed relays at FF position (unplanned wire delay),
        //   now places relays at interpolated positions between leaf buffer and FF.
        // BUG #3 fix: Hold guard now checks member FF hold budgets (min).
        //
        // Structure: branch → leaf_buffer → [relay_0] → [relay_1] → FF
        //            branch → leaf_buffer → FF  (zero target, direct)
        if (enablePerFfRelay_ && cts3dDb_ != nullptr
            && cts3dDb_->hasSkewTargets() && perFfBufDelay_ > 1e-9) {
          subNet->setLeafLevel(true);
          int branchRelayFFs = 0, branchRelayBufs = 0;
          int branchDirectFFs = 0, branchHoldSkip = 0;

          const int branch_tier = cts3dDb_->getDominantTier(
              leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
          const std::string bufMaster = cts3dDb_->getBufferForTier(
              options_->getRootBuffer(), branch_tier);

          // Helper lambda to process a single FF with relay decision
          auto processOneFf = [&](ClockInst* ffInst, ClockInst* driverInst,
                                  ClockSubNet* driverSubNet) {
            // Get per-FF target (NOT cluster mean)
            const double ffTarget = [&]() -> double {
              if (cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) return 0.0;
              std::string nm = ffInst->getName();
              const auto sp = nm.rfind('/');
              if (sp != std::string::npos) nm = nm.substr(0, sp);
              return cts3dDb_->getSkewTarget(nm);
            }();
            // globalMean fix: absolute LP target
            const double ffDelta = ffTarget;
            // Pre-compute dx/dy for Elmore model (moved before nRelay decision)
            const int dx = ffInst->getX() - driverInst->getX();
            const int dy = ffInst->getY() - driverInst->getY();

            // Elmore x_useful relay positioning.
            // When wire RC params are set (CTS_RELAY_RW_PER_UM / CTS_RELAY_CW_PER_UM),
            // compute exact relay position so 1 buffer delivers target_delta:
            //   delta(x) = d_buf - rc * x * (L - x)   [ps]
            //   x_useful = [L - sqrt(L^2 - 4*(d_buf-target)/rc)] / 2  [um]
            // where rc = rw[kOhm/um] * cw[fF/um] = ps/um^2 (kOhm*fF = 1ps).
            // For target > d_buf: fall back to multi-relay equal spacing.
            // For rc == 0: fall back to old integer-quantized V41 method.
            int nRelay = 0;
            double relayFrac = 0.5;  // position fraction for 1-relay Elmore case
            if (ffDelta > perFfMinTarget_) {
              const double d_buf_ps  = perFfBufDelay_ * 1000.0;  // ns -> ps
              const double target_ps = ffDelta * 1000.0;          // ns -> ps
              const double rc = relayRwKOhmPerUm_ * relayCwFfPerUm_;  // ps/um^2
              const bool elmore = (rc > 1e-12
                  && relayDbuPerUm_ > 1e-6
                  && (std::abs(dx) + std::abs(dy)) > 0);
              if (elmore) {
                const double L_um = (std::abs(dx) + std::abs(dy))
                    / relayDbuPerUm_;
                if (target_ps <= d_buf_ps && L_um > 1e-6) {
                  // 1 relay at x_useful covers target exactly
                  const double disc = L_um * L_um
                      - 4.0 * (d_buf_ps - target_ps) / rc;
                  if (disc >= 0.0) {
                    const double x_useful = (L_um - std::sqrt(disc)) / 2.0;
                    relayFrac = std::max(0.05,
                        std::min(0.95, x_useful / L_um));
                  }
                  // disc < 0 means target > max wire contribution → midpoint
                  nRelay = 1;
                } else if (L_um > 1e-6) {
                  // target > d_buf: multiple relays at equal spacing
                  nRelay = std::max(1, std::min(perFfMaxRelay_,
                      static_cast<int>(
                          std::round(target_ps / d_buf_ps))));
                }
              } else {
                // Fallback: integer quantization (V41 method)
                nRelay = static_cast<int>(
                    std::round(ffDelta / perFfBufDelay_));
                nRelay = std::max(0, std::min(nRelay, perFfMaxRelay_));
              }
            }

            // Removed redundant holdSkip (same as above).

            // Debug stats (use delta, not absolute, for accurate coverage)
            totalLpDelay_debug += std::max(0.0, ffDelta);
            ++totalBranchSinks_debug;

            if (nRelay > 0) {
              // dx/dy pre-computed above (before nRelay decision).
              // Single relay uses Elmore relayFrac; multi-relay uses equal spacing.
              // Get FF's tier for relay tier assignment.
              // Relay must be on the SAME tier as its target FF to avoid cross-tier
              // relay→FF nets. Cross-tier relay→FF nets cause GRT-0183 heap underflow
              // in 3D maze routing (observed in V43/V46/V47 for minority-tier FFs in
              // upper-dominant clusters). Using branch_tier (cluster dominant tier)
              // was wrong when FF is on the minority tier of the cluster.
              const int ff_tier_relay = cts3dDb_->getInstTier(ffInst->getDbInst());
              const int relay_tier = (ff_tier_relay >= 0) ? ff_tier_relay : branch_tier;

              ClockSubNet* cur = driverSubNet;
              for (int r = 0; r < nRelay; ++r) {
                // 1-relay → Elmore x_useful fraction; else equal spacing
                const double frac = (nRelay == 1)
                    ? relayFrac
                    : static_cast<double>(r + 1) / (nRelay + 1);
                const int rx = driverInst->getX()
                    + static_cast<int>(dx * frac);
                const int ry = driverInst->getY()
                    + static_cast<int>(dy * frac);

                const std::string rName = tierPrefix_ + "clkbuf_relay_"
                    + std::to_string(idx) + "_"
                    + std::to_string(numSinks) + "_"
                    + std::to_string(r);
                ClockInst& relay = clock_.addClockBuffer(
                    rName, bufMaster, rx, ry);
                cts3dDb_->setClockInstTier(relay, relay_tier);
                addTreeLevelBuffer(&relay);
                cur->addInst(relay);
                ClockSubNet& next = clock_.addSubNet(
                    tierPrefix_ + "clknet_relay_" + std::to_string(idx) + "_"
                    + std::to_string(numSinks) + "_"
                    + std::to_string(r));
                next.addInst(relay);
                cur = &next;
              }
              cur->setLeafLevel(true);
              cur->addInst(*ffInst);
              ++branchRelayFFs;
              branchRelayBufs += nRelay;
              totalDelivered_debug += nRelay * perFfBufDelay_;
            } else {
              driverSubNet->addInst(*ffInst);
              ++branchDirectFFs;
            }
            ++numSinks;
          };

          const std::vector<Point<double>>& sinkLocs
              = leafTopology.getBranchSinksLocations(idx);
          for (const auto& loc : sinkLocs) {
            // Epsilon-tolerant sink lookup (V37 CTS-0080 fix reused)
            // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 0.01;
              auto hint = mapLocationToSink_.lower_bound(
                  Point<double>(loc.getX() - eps, loc.getY() - eps));
              for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                if (it->first.getX() > loc.getX() + eps) break;
                if (std::abs(it->first.getX() - loc.getX()) < eps
                    && std::abs(it->first.getY() - loc.getY()) < eps) {
                  sinkIt = it;
                  break;
                }
              }
            }
            if (sinkIt == mapLocationToSink_.end()) {
              logger_->error(CTS, 425,
                  "Relay: Sink not found at ({:.10f}, {:.10f})",
                  loc.getX(), loc.getY());
            }
            ClockInst* sinkInst = sinkIt->second;

            // Check if this is a cluster buffer
            auto cIt = clusterBufTarget_.find(sinkInst);
            if (cIt != clusterBufTarget_.end()) {
              // Cluster buffer: add to branch subnet, then process
              // each member FF individually with per-FF relay chains.
              subNet->addInst(*sinkInst);

              // Find the cluster's leaf subnet
              ClockSubNet* leafSubNet = nullptr;
              clock_.forEachSubNet([&](ClockSubNet& sn) {
                if (sn.getDriver() == sinkInst) {
                  leafSubNet = &sn;
                }
              });

              if (leafSubNet != nullptr) {
                // Collect member FFs (can't modify while iterating)
                std::vector<ClockInst*> memberFFs;
                leafSubNet->forEachSink([&](ClockInst* ff) {
                  memberFFs.push_back(ff);
                });
                // Remove all sinks from leaf subnet (will re-add or relay)
                std::set<ClockInst*> toRemove(
                    memberFFs.begin(), memberFFs.end());
                leafSubNet->removeSinks(toRemove);

                // Process each member FF individually
                for (ClockInst* memberFF : memberFFs) {
                  processOneFf(memberFF, sinkInst, leafSubNet);
                }
              } else {
                // Leaf subnet not found (shouldn't happen)
                ++branchDirectFFs;
                ++numSinks;
              }
            } else {
              // Individual FF (no clustering): process directly
              processOneFf(sinkInst, sinkInst, subNet);
            }
          }

          numDelayBufs += branchRelayBufs;
          ++totalBranches_debug;
          if (branchRelayFFs > 0) ++firedBranches_debug;
          logger_->info(CTS, 424,
              "Relay branch {:2d}: relay={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns",
              idx, branchRelayFFs, branchRelayBufs,
              branchDirectFFs, branchHoldSkip, absTgt);
          return;  // skip legacy per-branch N-buffer code
        }

        // globalMean fix: absolute LP target
        const double delta = absTgt;
        int num_bufs = 0;
        // V58 Fix F: Skip per-branch delay buffers for MacroTree.
        // Macro tree latency is already too high (1,095ps vs register ~462ps).
        // Adding delay buffers widens the gap → worsens hold on reg→macro edges.
        if (type_ != TreeType::MacroTree
            && delta > 0.0 && singleBufDelay_ > 1e-9) {
          num_bufs = static_cast<int>(std::round(delta / singleBufDelay_));
          num_bufs = std::max(0, std::min(num_bufs, maxLeafDelayBufs_));
        }

        // Debug: accumulate LP-CTS gap stats (CTS-408 summary after loop)
        totalLpDelay_debug   += absTgt * branchCnt;
        totalBranchSinks_debug += branchCnt;
        ++totalBranches_debug;

        if (num_bufs > 0) {
            totalDelivered_debug += num_bufs * singleBufDelay_ * branchCnt;
            ++firedBranches_debug;

            const int branch_tier = cts3dDb_->getDominantTier(
                leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
            const std::string delayBufMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), branch_tier);

            const int bx = branchPoint.getX() * wireSegmentUnit_;
            const int by = branchPoint.getY() * wireSegmentUnit_;

            // Build N-buffer chain in series
            ClockSubNet* curSub = subNet;
            for (int b = 0; b < num_bufs; ++b) {
              const std::string bufName = tierPrefix_ + "clkbuf_delay_"
                                          + std::to_string(idx) + "_"
                                          + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);

              // buf is SINK of curSub (curSub already has a driver)
              curSub->addInst(buf);

              // Create new subnet; first addInst sets buf as its driver
              ClockSubNet& nextSub = clock_.addSubNet(
                  tierPrefix_ + "clknet_delay_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);  // buf becomes driver of nextSub

              curSub = &nextSub;
            }

            // Last subnet is the leaf subnet; connect FFs as sinks
            curSub->setLeafLevel(true);
            const std::vector<Point<double>>& sinkLocs
                = leafTopology.getBranchSinksLocations(idx);
            for (const Point<double>& loc : sinkLocs) {
              auto sinkIt = mapLocationToSink_.find(loc);
              // epsilon fallback for float precision mismatch
              // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
              if (sinkIt == mapLocationToSink_.end()) {
                constexpr double eps = 0.01;
                auto hint = mapLocationToSink_.lower_bound(
                    Point<double>(loc.getX() - eps, loc.getY() - eps));
                for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                  if (it->first.getX() > loc.getX() + eps) break;
                  if (std::abs(it->first.getX() - loc.getX()) < eps
                      && std::abs(it->first.getY() - loc.getY()) < eps) {
                    sinkIt = it;
                    break;
                  }
                }
              }
              if (sinkIt == mapLocationToSink_.end()) {
                // V58: eps failed — find absolute nearest (same pattern as CTS-0079).
                // Cannot skip: sink would be unconnected → SIGSEGV in writeDummyLoadsToDb.
                double bestDist = std::numeric_limits<double>::max();
                for (auto& [key, sink] : mapLocationToSink_) {
                  double d = std::abs(key.getX() - loc.getX())
                           + std::abs(key.getY() - loc.getY());
                  if (d < bestDist) {
                    bestDist = d;
                    sinkIt = mapLocationToSink_.find(key);
                  }
                }
                logger_->warn(CTS, 404,
                    "Delay buf branch: sink at ({:.6f}, {:.6f}) — nearest match dist={:.6f}",
                    loc.getX(), loc.getY(), bestDist);
              }
              curSub->addInst(*sinkIt->second);
              ++numSinks;
            }

            ++numDelayBufs;
            logger_->info(CTS, 402,
                          "Delay buffers inserted at leaf branch {}: "
                          "absTarget={:.4f}ns, globalMean={:.4f}ns, delta={:.4f}ns, N={}, buffer={}",
                          idx, absTgt, globalMeanTarget_, delta, num_bufs, delayBufMaster);
            return;  // skip normal sink connection below
          }

        // Branch has LP target but delta ≤ 0 → skipped by delta trigger.
        // These branches are BELOW globalMean; N=0 even though LP assigned non-zero target.
        // This reveals LP targets that are "wasted" (below average, can't be implemented
        // with delta-based approach). Key gap source if most targets fall here.
        if (absTgt > 0.0) {
          logger_->info(CTS, 406,
                        "  Branch {:2d} skipped: absTgt={:.4f}ns <= globalMean={:.4f}ns "
                        "(delta={:.4f}ns -> N=0, {} sinks; "
                        "LP wanted {:.4f}ns/sink but CTS delivers 0ns)",
                        idx, absTgt, globalMeanTarget_, delta, branchCnt, absTgt);
        }

        subNet->setLeafLevel(true);

        const std::vector<Point<double>>& sinkLocs
            = leafTopology.getBranchSinksLocations(idx);
        for (const Point<double>& loc : sinkLocs) {
          auto sinkIt = mapLocationToSink_.find(loc);
          // float/double precision mismatch from multi-level
          // computeBranchSinks() double→float→double roundtrip.
          // If exact lookup fails, search for nearest entry within epsilon.
          // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
          if (sinkIt == mapLocationToSink_.end()) {
            constexpr double eps = 0.01;
            auto hint = mapLocationToSink_.lower_bound(
                Point<double>(loc.getX() - eps, loc.getY() - eps));
            for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
              if (it->first.getX() > loc.getX() + eps) break;
              if (std::abs(it->first.getX() - loc.getX()) < eps
                  && std::abs(it->first.getY() - loc.getY()) < eps) {
                sinkIt = it;
                debugPrint(logger_, CTS, "clustering", 1,
                           "Leaf branch {}: epsilon match ({:.10f},{:.10f}) -> "
                           "({:.10f},{:.10f}) = {}",
                           idx, loc.getX(), loc.getY(),
                           it->first.getX(), it->first.getY(),
                           it->second->getName());
                break;
              }
            }
          }
          if (sinkIt == mapLocationToSink_.end()) {
            // V58: warn+skip instead of fatal error.
            // Multi-level clustering float precision drift can exceed eps=1e-4
            // for large designs (10K+ sinks, 7+ H-tree levels).
            // Skipped sink stays on parent leaf net (depth=0, no grouped delay).
            logger_->warn(CTS, 80,
                          "Sink not found at ({:.10f}, {:.10f}), mapSize={} — "
                          "skipping (connects to leaf net directly)",
                          loc.getX(), loc.getY(), mapLocationToSink_.size());
            continue;
          }

          subNet->addInst(*sinkIt->second);
          ++numSinks;
        }
      });

  logger_->info(CTS, 35, " Number of sinks covered: {}.", numSinks);
  if (numDelayBufs > 0) {
    logger_->info(CTS, 403,
                  " {} delay buffers inserted for per-branch skew targeting.",
                  numDelayBufs);
  }

  // LP-CTS translation efficiency summary.
  // coverage = totalDelivered / totalLpDelay * 100%
  //   ~0%   → N-buffer chain never fires (all branches below globalMean, or absTgt=0)
  //   ~50%  → half of LP target translated (rounding + below-mean branches)
  //   ~100% → N-buffer chain faithfully implements LP targets
  // This is the primary metric for diagnosing LP→CTS pipeline effectiveness.
  if (totalBranchSinks_debug > 0 && totalLpDelay_debug > 1e-12) {
    const double coverage = 100.0 * totalDelivered_debug / totalLpDelay_debug;
    logger_->info(CTS, 408,
                  "LP-CTS gap: LP_total={:.4f}ns ({:.4f}ns/sink), "
                  "delivered={:.4f}ns ({:.4f}ns/sink), coverage={:.1f}%, "
                  "fired={}/{} leaf branches",
                  totalLpDelay_debug,
                  totalLpDelay_debug / totalBranchSinks_debug,
                  totalDelivered_debug,
                  totalDelivered_debug / totalBranchSinks_debug,
                  coverage,
                  firedBranches_debug, totalBranches_debug);
  }

  // Per-FF relay summary log
  if (enablePerFfRelay_ && totalBranchSinks_debug > 0) {
    const double cov = totalLpDelay_debug > 1e-12
        ? 100.0 * totalDelivered_debug / totalLpDelay_debug : 0.0;
    logger_->info(CTS, 426,
        "Per-FF relay summary: {} relay bufs across {}/{} branches, "
        "{} sinks, coverage={:.1f}%",
        numDelayBufs, firedBranches_debug, totalBranches_debug,
        numSinks, cov);
  }

  // Realizability guard pass — must run after ALL sink-to-net assignment
  // (leaf buffers, TAP chains, delay buffers) so every subnet is final.
  insertRealizabilityGuards();
}

// Per-sink realizability guard.
// MEASURED root cause (swerv post-route): per-sink implementation error scales
// as R_route(driver->sink) x C_pin. Geometry-only clustering put 4 SRAM clk
// pins (25fF each vs ~1.2fF FF median) on one leaf buffer 76um away -> 107um
// leaf net, 68% on M2/M3 -> 3,766 Ohm (5.7x FF-leaf baseline) x 105fF ~ 400ps
// delay + 715ps slew, invisible to the CTS-stage estimate. Splitting macros
// into a separate tree did NOT fix it (post-route -404,785 vs -404,682):
// gated-clock SRAMs bypass top-net classification, and a macro tree reuses the
// same load-blind leaf machinery. The general fix is LOCAL: any sink whose
// input-pin cap >> median gets a dedicated guard buffer AT the sink location,
// so the shared leaf net sees the guard's ~1fF input instead of the heavy pin,
// and the heavy pin is driven over a ~zero-length private stub (R ~ 0).
// No-op for all-FF designs (threshold never met). No per-design config.
void HTreeBuilder::insertRealizabilityGuards()
{
  bool enable = true;
  if (const char* e = std::getenv("CTS_ENABLE_REALIZABILITY_GUARD")) {
    enable = (std::atoi(e) != 0);
  }
  if (!enable) {
    return;
  }

  // 1) Median sink input cap over this clock tree.
  std::vector<float> caps;
  clock_.forEachSink([&](const ClockInst& sink) {
    if (sink.getInputCap() > 0.0f) {
      caps.push_back(sink.getInputCap());
    }
  });
  if (caps.size() < 2) {
    return;
  }
  std::sort(caps.begin(), caps.end());
  const float medianCap = caps[caps.size() / 2];
  // Fire only for genuinely heavy pins: >=10x the median AND >= 5fF absolute.
  // getInputPinCap() -> sta::LibertyPort::capacitance() -> FARADS (OpenSTA
  // normalizes all Liberty units to SI): FF CLK ~1.5e-15, SRAM clk ~2.5e-14.
  // (v1 bug: floor 0.005 assumed pF units -> 0.005 F -> guard never fired.)
  const float threshold = std::max(10.0f * medianCap, 5e-15f);

  // 2) Rewire heavy sinks behind a guard buffer, per subnet.
  //    Collect first (can't mutate subnet list while iterating).
  struct GuardTask
  {
    ClockSubNet* net;
    ClockInst* sink;
  };
  std::vector<GuardTask> tasks;
  clock_.forEachSubNet([&](ClockSubNet& subNet) {
    // Skip guard nets themselves (idempotent under per-tier double calls).
    if (subNet.getName().find("rguard") != std::string::npos) {
      return;
    }
    subNet.forEachSink([&](ClockInst* sink) {
      if (sink == nullptr || sink->isClockBuffer()) {
        return;
      }
      if (sink->getInputCap() > threshold
          && rguardDone_.find(sink) == rguardDone_.end()) {
        tasks.push_back({&subNet, sink});
      }
    });
  });
  // Always log the survey (doubles as the FF-only no-op proof: guards=0).
  const float maxCap = caps.back();
  logger_->info(CTS, 720,
      "Realizability guard survey: {} sinks, median={:.2f}fF max={:.2f}fF "
      "threshold={:.2f}fF -> {} heavy sink(s).",
      static_cast<int>(caps.size()), medianCap * 1e15, maxCap * 1e15,
      threshold * 1e15, static_cast<int>(tasks.size()));
  if (tasks.empty()) {
    return;
  }

  const std::string rootBufMaster = options_->getRootBuffer();
  int made = 0;
  for (auto& t : tasks) {
    ClockInst* sink = t.sink;
    // Tier-correct buffer master; fall back to sink's stored tier then driver.
    int sinkTier = sink->getTier();
    if (sinkTier < 0 && cts3dDb_ != nullptr
        && sink->getDbInputPin() != nullptr) {
      sinkTier = cts3dDb_->getInstTier(sink->getDbInputPin()->getInst());
    }
    std::string master = rootBufMaster;
    if (cts3dDb_ != nullptr && sinkTier >= 0) {
      master = cts3dDb_->getBufferForTier(rootBufMaster, sinkTier);
    }
    // Place the guard AT the sink pin; legalizeOneBuffer snaps off macro
    // blockages and resolves collisions (DPL re-legalizes onto rows later).
    Point<double> rawLoc(
        static_cast<double>(sink->getX()) / wireSegmentUnit_,
        static_cast<double>(sink->getY()) / wireSegmentUnit_);
    Point<double> legalLoc = legalizeOneBuffer(rawLoc, master);
    const int gx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
    const int gy = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);

    // Unique name across builders AND across per-tier createClockSubNets
    // calls: rguardCount_ is a member (v1 bug: local counter restarted at 0
    // on the second per-tier call -> duplicate ODB inst names).
    std::string cn = clock_.getName();
    std::replace(cn.begin(), cn.end(), '/', '_');
    std::replace(cn.begin(), cn.end(), ':', '_');
    const std::string bufName = tierPrefix_ + "clkbuf_rguard_" + cn + "_"
                                + std::to_string(rguardCount_);
    ClockInst& guard = clock_.addClockBuffer(bufName, master, gx, gy);
    if (cts3dDb_ != nullptr && sinkTier >= 0) {
      cts3dDb_->setClockInstTier(guard, sinkTier);
    }
    addTreeLevelBuffer(&guard);

    // Audit distances: how far the ORIGINAL driver was vs the new stub.
    ClockInst* origDriver = t.net->getDriver();
    const long odx = origDriver ? std::labs((long) origDriver->getX() - sink->getX()) : -1;
    const long ody = origDriver ? std::labs((long) origDriver->getY() - sink->getY()) : -1;
    const long gdx = std::labs((long) gx - sink->getX());
    const long gdy = std::labs((long) gy - sink->getY());

    // leaf net now drives the guard instead of the heavy sink...
    t.net->replaceSink(sink, &guard);
    // ...and the guard drives the heavy sink over a private stub.
    ClockSubNet& gNet = clock_.addSubNet(tierPrefix_ + "clknet_rguard_" + cn
                                         + "_" + std::to_string(rguardCount_));
    gNet.addInst(guard);
    gNet.addInst(*sink);
    gNet.setLeafLevel(true);
    rguardDone_.insert(sink);
    rguardCount_++;
    made++;

    // Requirement-10 audit log: heavy sink, caps, tier/loc, guard, distances,
    // and replacement confirmation (guard in the original net, sink out).
    bool replaced = false;
    try {
      t.net->findIndex(&guard);   // throws if guard not in net
      replaced = true;
      t.net->findIndex(sink);     // must throw (sink removed)
      replaced = false;           // sink still present -> replacement FAILED
    } catch (const std::out_of_range&) {
      // expected on the second lookup: guard in, sink out
    }
    logger_->info(CTS, 721,
        "  rguard {} [{}] tier={} at ({},{}): sink {} cap={:.2f}fF "
        "(median {:.2f}fF) origDriver_dist=({}, {}) dbu -> stub_dist=({}, {}) "
        "dbu, net={}, replaced={}",
        bufName, master, sinkTier, gx, gy, sink->getName(),
        sink->getInputCap() * 1e15, medianCap * 1e15,
        odx, ody, gdx, gdy, t.net->getName(), replaced ? 1 : 0);
  }
  logger_->info(CTS, 722,
      "Realizability guard: {} high-cap sinks isolated behind dedicated "
      "guard buffers placed at the sink (threshold {:.2f}fF = max(10x median, "
      "5fF)).",
      made, threshold * 1e15);
}

void HTreeBuilder::createSingleBufferClockNet()
{
  logger_->report(" Building single-buffer clock net.");

  Point<double> center = sinkRegion_.getCenter();
  Point<double> legalCenter
      = legalizeOneBuffer(center, options_->getRootBuffer());
  sinkRegion_.setCenter(legalCenter);
  commitMoveLoc(center, legalCenter);
  const int centerX = legalCenter.getX() * wireSegmentUnit_;
  const int centerY = legalCenter.getY() * wireSegmentUnit_;
  ClockInst& rootBuffer = clock_.addClockBuffer(
      tierPrefix_ + "clkbuf_0", options_->getRootBuffer(), centerX, centerY);

  // clang-format off
  if (center != legalCenter) {
    debugPrint(logger_, CTS, "legalizer", 2, "createSingleBufferClockNet "
	       "legalizeOneBuffer {}clkbuf_0: {} => {}", tierPrefix_, center, legalCenter);
  }
  // clang-format on

  addTreeLevelBuffer(&rootBuffer);
  ClockSubNet& clockSubNet = clock_.addSubNet(tierPrefix_ + "clknet_0");
  clockSubNet.setTreeLevel(0);
  clockSubNet.addInst(rootBuffer);

  clock_.forEachSink([&](ClockInst& inst) { clockSubNet.addInst(inst); });
}

void HTreeBuilder::plotSolution()
{
  auto name = std::string("plot_") + clock_.getName() + ".py";
  std::ofstream file(name);
  file << "import numpy as np\n";
  file << "import matplotlib.pyplot as plt\n";
  file << "import matplotlib.path as mpath\n";
  file << "import matplotlib.lines as mlines\n";
  file << "import matplotlib.patches as mpatches\n";
  file << "from matplotlib.collections import PatchCollection\n\n";

  clock_.forEachSink([&](const ClockInst& sink) {
    file << "plt.scatter(" << (double) sink.getX() / wireSegmentUnit_ << ", "
         << (double) sink.getY() / wireSegmentUnit_ << ", s=1)\n";
  });

  LevelTopology& topLevelTopology = topologyForEachLevel_.front();
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  topLevelTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        if (topLevelBufferLoc.getX() < branchPoint.getX()) {
          file << "plt.plot([" << topLevelBufferLoc.getX() << ", "
               << branchPoint.getX() << "], [" << topLevelBufferLoc.getY()
               << ", " << branchPoint.getY() << "], c = 'r')\n";
        } else {
          file << "plt.plot([" << branchPoint.getX() << ", "
               << topLevelBufferLoc.getX() << "], [" << branchPoint.getY()
               << ", " << topLevelBufferLoc.getY() << "], c = 'r')\n";
        }
      });

  for (int levelIdx = 1; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    const LevelTopology& topology = topologyForEachLevel_[levelIdx];
    topology.forEachBranchingPoint([&](unsigned idx,
                                       Point<double> branchPoint) {
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
      Point<double> parentPoint
          = topologyForEachLevel_[levelIdx - 1].getBranchingPoint(parentIdx);
      std::string color = "orange";
      if (levelIdx % 2 == 0) {
        color = "red";
      }

      if (parentPoint.getX() < branchPoint.getX()) {
        file << "plt.plot([" << parentPoint.getX() << ", " << branchPoint.getX()
             << "], [" << parentPoint.getY() << ", " << branchPoint.getY()
             << "], c = '" << color << "')\n";
      } else {
        file << "plt.plot([" << branchPoint.getX() << ", " << parentPoint.getX()
             << "], [" << branchPoint.getY() << ", " << parentPoint.getY()
             << "], c = '" << color << "')\n";
      }
    });
  }

  file << "plt.show()\n";
  file.close();
}

// print structures of Htree from top level buffer
// incluiding branch point locations, topology length and weighted sink lengths
void HTreeBuilder::printHTree()
{
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  logger_->report("HTree: top buf loc:{}", topLevelBufferLoc);
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
      Point<double>& branchPoint = topology.getBranchingPoint(idx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? topLevelBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(idx);

      double leng = topology.getLength();
      // clang-format off
      logger_->report("HTree: level*{}* bufId*{}*: branchPt:{} topo len:{:0.3f}"
		      " dist to parent:{:0.3f} weighted sink len:{:0.3f} "
		      "parentPt:{}", levelIdx, idx, branchPoint, leng,
		      computeDist(branchPoint, parentPoint),
		      weightedDistance(branchPoint, branchPoint, sinks),
		      parentPoint);
      // clang-format on
    }
    logger_->report("-------------------------------------------------");
  }
}

SegmentBuilder::SegmentBuilder(const std::string& instPrefix,
                               const std::string& netPrefix,
                               const Point<double>& root,
                               const Point<double>& target,
                               const std::vector<unsigned>& techCharWires,
                               Clock& clock,
                               ClockSubNet& drivingSubNet,
                               const TechChar& techChar,
                               const unsigned techCharDistUnit,
                               TreeBuilder* tree,
                               odb::dbDatabase* db,
                               int targetTier)
    : instPrefix_(instPrefix),
      netPrefix_(netPrefix),
      root_(root),
      target_(target),
      techCharWires_(techCharWires),
      techChar_(&techChar),
      techCharDistUnit_(techCharDistUnit),
      clock_(&clock),
      drivingSubNet_(&drivingSubNet),
      tree_(tree),
      db_(db),
      targetTier_(targetTier)
{
}

void SegmentBuilder::build(const std::string& forceBuffer)
{
  const double lengthX = std::abs(root_.getX() - target_.getX());
  const bool isLowToHiX = root_.getX() < target_.getX();
  const bool isLowToHiY = root_.getY() < target_.getY();

  double connectionLength = 0.0;
  for (unsigned techCharWireIdx : techCharWires_) {
    const WireSegment& wireSegment = techChar_->getWireSegment(techCharWireIdx);
    const unsigned wireSegLen = wireSegment.getLength();
    for (int buffer = 0; buffer < wireSegment.getNumBuffers(); ++buffer) {
      const double location
          = wireSegment.getBufferLocation(buffer) * wireSegLen;
      connectionLength += location;

      double x = std::numeric_limits<double>::max();
      double y = std::numeric_limits<double>::max();
      if (connectionLength < lengthX) {
        y = root_.getY();
        x = (isLowToHiX) ? (root_.getX() + connectionLength)
                         : (root_.getX() - connectionLength);
      } else {
        x = target_.getX();
        y = (isLowToHiY) ? (root_.getY() + (connectionLength - lengthX))
                         : (root_.getY() - (connectionLength - lengthX));
      }

      const std::string buffMaster = !forceBuffer.empty()
                                         ? forceBuffer
                                         : wireSegment.getBufferMaster(buffer);
      // Replaced mapBufferMasterToTier with Cts3DDatabase
      const std::string tieredMaster
          = tree_->getCts3DDatabase()->getBufferForTier(buffMaster, targetTier_);
      Point<double> bufferLoc(x, y);
      Point<double> legalBufferLoc
          = tree_->legalizeOneBuffer(bufferLoc, tieredMaster);
      tree_->commitMoveLoc(bufferLoc, legalBufferLoc);
      ClockInst& newBuffer = clock_->addClockBuffer(
          instPrefix_ + std::to_string(numBufferLevels_),
          tieredMaster,
          legalBufferLoc.getX() * techCharDistUnit_,
          legalBufferLoc.getY() * techCharDistUnit_);
      // Tier assignment through SSOT (Cts3DDatabase)
      tree_->getCts3DDatabase()->setClockInstTier(newBuffer, targetTier_);

      // clang-format off
      if (bufferLoc != legalBufferLoc) {
	// adjust for cell movement
	connectionLength -= tree_->computeDist(bufferLoc, legalBufferLoc);
	debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
		   " SegmentBuilder::build {} TCId:{} bufId:{} connLen:{:0.1f}: "
		   "{} => {}", instPrefix_  + std::to_string(numBufferLevels_),
		   techCharWireIdx, buffer, connectionLength,
		   bufferLoc, legalBufferLoc);
      } else {
	debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
		   " SegmentBuilder::build {} TCId:{} bufId:{} connLen:{:0.1f}: "
		   "{}", instPrefix_  + std::to_string(numBufferLevels_),
		   techCharWireIdx, buffer, connectionLength, bufferLoc);
      }
      // clang-format on

      tree_->addTreeLevelBuffer(&newBuffer);

      drivingSubNet_->addInst(newBuffer);
      drivingSubNet_
          = &clock_->addSubNet(netPrefix_ + std::to_string(numBufferLevels_));
      drivingSubNet_->addInst(newBuffer);

      ++numBufferLevels_;
    }
    connectionLength += wireSegLen;
  }
}

void SegmentBuilder::forceBufferInSegment(const std::string& master)
{
  if (numBufferLevels_ != 0) {
    return;
  }

  // Replaced mapBufferMasterToTier with Cts3DDatabase
  const std::string tieredMaster
      = tree_->getCts3DDatabase()->getBufferForTier(master, targetTier_);
  ClockInst& newBuffer
      = clock_->addClockBuffer(instPrefix_ + "_f",
                               tieredMaster,
                               target_.getX() * techCharDistUnit_,
                               target_.getY() * techCharDistUnit_);
  // Tier assignment through SSOT (Cts3DDatabase)
  tree_->getCts3DDatabase()->setClockInstTier(newBuffer, targetTier_);
  tree_->addTreeLevelBuffer(&newBuffer);
  // clang-format off
  debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
	     "  forceBufferInSegment {}: {}", instPrefix_ + "_f", target_);
  // clang-format on

  drivingSubNet_->addInst(newBuffer);
  drivingSubNet_ = &clock_->addSubNet(netPrefix_ + "_leaf");
  drivingSubNet_->addInst(newBuffer);
  numBufferLevels_++;
}

// V86_SEED: off-cluster-seeded refinement.
// Move ONLY FFs incident to violating timing edges (slack_max < crit_slack)
// into the cluster of a violating-edge partner, iff the affinity gain is
// positive and hard constraints hold (dest size < maxClusterSize, dest
// bounding-box diagonal <= maxDiameter, dest tier purity, no new/widened
// V37 target-split). Zero moves leave `solution` untouched -> byte-identical
// to TARP-off downstream. See TARP_ALPHA_SWEEP_RESULTS.md ADDENDUM
// for the evidence chain (V85 CRIT_ONLY bisection-floor) motivating this.
void HTreeBuilder::refineSeedClusters(
    std::vector<std::vector<unsigned>>& solution,
    const std::vector<std::pair<float, float>>& points,
    const std::vector<const ClockInst*>& sinkInsts,
    unsigned maxClusterSize,
    float maxDiameter)
{
  const char* csv_env = std::getenv("CTS_FF_TIMING_GRAPH_CSV");
  if (!csv_env) {
    logger_->warn(CTS, 847,
        "Seeded TARP: CTS_FF_TIMING_GRAPH_CSV not set — 0 moves (seed kept).");
    return;
  }

  const auto refine_t0 = std::chrono::high_resolution_clock::now();
  // Knobs: clean CTS_SEED_* names take precedence; old CTS_TARP_SEED_* kept
  // as back-compat fallbacks.
  double critSlackNs = 0.0;
  if (const char* e = std::getenv("CTS_TARP_CRIT_SLACK_NS"))
    critSlackNs = std::stod(e);
  if (const char* e = std::getenv("CTS_SEED_CRIT_SLACK_NS"))
    critSlackNs = std::stod(e);
  long maxMoves = -1;  // -1 = unlimited
  if (const char* e = std::getenv("CTS_TARP_SEED_MAX_MOVES"))
    maxMoves = std::atol(e);
  if (const char* e = std::getenv("CTS_SEED_MAX_MOVES"))
    maxMoves = std::atol(e);

  // T_period: SDC-derived default (V85), env override.
  double T_period
      = (options_->getSdcPeriodNs() > 0.0) ? options_->getSdcPeriodNs() : 1.0;
  if (const char* e = std::getenv("CTS_CLOCK_PERIOD_NS"))
    T_period = std::stod(e);

  const unsigned N = points.size();

  // FF name (pin suffix stripped) -> point index. Same convention as
  // runTarpClustering (see addSink loop).
  std::unordered_map<std::string, unsigned> nameToIdx;
  nameToIdx.reserve(N * 2);
  for (unsigned i = 0; i < N; ++i) {
    std::string nm = sinkInsts[i]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);
    nameToIdx[nm] = i;
  }

  // point index -> cluster index
  std::vector<int> clusterOf(N, -1);
  for (size_t c = 0; c < solution.size(); ++c) {
    for (unsigned idx : solution[c]) {
      if (idx < N) clusterOf[idx] = static_cast<int>(c);
    }
  }

  // Die diagonal in normalized coords (w_spatial normalization, same model
  // as TarpClustering::computeDieDiagonal).
  double xMin = 1e30, xMax = -1e30, yMin = 1e30, yMax = -1e30;
  for (const auto& pt : points) {
    xMin = std::min(xMin, (double) pt.first);
    xMax = std::max(xMax, (double) pt.first);
    yMin = std::min(yMin, (double) pt.second);
    yMax = std::max(yMax, (double) pt.second);
  }
  double dieDiag = std::sqrt((xMax - xMin) * (xMax - xMin)
                             + (yMax - yMin) * (yMax - yMin));
  if (dieDiag < 1e-9) dieDiag = 1.0;

  const double alpha = tarpAlpha_;
  auto affinity = [&](unsigned a, unsigned b, double slack) -> double {
    const double dx = std::abs((double) points[a].first - points[b].first);
    const double dy = std::abs((double) points[a].second - points[b].second);
    const double w_spatial = std::max(0.0, 1.0 - (dx + dy) / dieDiag);
    const double w_timing
        = std::min(2.0, std::max(0.0, 1.0 - slack / T_period));
    return alpha * w_spatial + (1.0 - alpha) * w_timing;
  };

  // Pass 1: violating edges (same-tier, non-ghost, both endpoints known).
  struct VEdge { double slack; unsigned a, b; };
  std::vector<VEdge> violating;
  std::unordered_set<unsigned> critSet;
  auto parseCsv = [&](auto&& perEdge) {
    std::ifstream csvFile(csv_env);
    if (!csvFile.is_open()) return false;
    std::string line;
    std::getline(csvFile, line);  // header
    while (std::getline(csvFile, line)) {
      std::vector<std::string> tokens;
      std::istringstream ls(line);
      std::string token;
      while (std::getline(ls, token, ',')) tokens.push_back(token);
      if (tokens.size() < 14) continue;
      const int f_tier = std::stoi(tokens[10]);
      const int t_tier = std::stoi(tokens[13]);
      if (f_tier != t_tier) continue;  // same-tier only (edges cross-tier
                                       // cannot be regrouped within a tier)
      const double slack_max = std::stod(tokens[2]);
      const double slack_min = std::stod(tokens[3]);
      const double arr_max = std::stod(tokens[4]);
      const double arr_min = std::stod(tokens[5]);
      if (slack_max == 0 && slack_min == 0 && arr_max == 0 && arr_min == 0)
        continue;  // ghost edge (unmeasured pair)
      auto it_a = nameToIdx.find(tokens[0]);
      auto it_b = nameToIdx.find(tokens[1]);
      if (it_a == nameToIdx.end() || it_b == nameToIdx.end()) continue;
      if (it_a->second == it_b->second) continue;
      perEdge(it_a->second, it_b->second, slack_max);
    }
    return true;
  };

  if (!parseCsv([&](unsigned a, unsigned b, double slack) {
        if (slack < critSlackNs) {
          violating.push_back({slack, a, b});
          critSet.insert(a);
          critSet.insert(b);
        }
      })) {
    logger_->warn(CTS, 849,
        "Seeded TARP: cannot open {} — 0 moves (seed kept).", csv_env);
    return;
  }

  if (violating.empty()) {
    logger_->info(CTS, 848,
        "Seeded TARP: 0 violating edges (crit_slack={:.4f}ns) — seed kept, "
        "byte-identical to TARP-off.", critSlackNs);
    return;
  }

  // Per-point LP target lookup (used by the compatibility factor and H4).
  const bool haveTargets = (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets());
  auto targetOf = [&](unsigned idx) -> double {
    std::string nm = sinkInsts[idx]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);
    return cts3dDb_->getSkewTarget(nm);
  };

  // V86 review fix: criticality != clustering compatibility.
  // A critical launch/capture pair whose LP targets differ should NOT be
  // co-clustered — in this flow delivery is per-FF (rung tap-in), so the harm
  // of co-clustering concentrates where |a*_u - a*_v| is large (within-leaf
  // target variance; see l3graph evidence). Optional compatibility factor:
  //   aff' = aff * exp(-|a*_u - a*_v| / tau)
  // CTS_TARP_SEED_COMPAT_TAU_NS > 0 enables; default 0 = exactly the previous
  // behavior (factor == 1).
  double compatTauNs = 0.0;
  if (const char* e = std::getenv("CTS_TARP_SEED_COMPAT_TAU_NS"))
    compatTauNs = std::stod(e);
  if (const char* e = std::getenv("CTS_SEED_COMPAT_TAU_NS"))
    compatTauNs = std::stod(e);
  std::vector<double> lpTarget;
  if (compatTauNs > 0.0 && haveTargets) {
    lpTarget.resize(N, 0.0);
    for (unsigned i = 0; i < N; ++i) lpTarget[i] = targetOf(i);
  }

  // Pass 2: adjacency restricted to critical FFs (gain evaluation only needs
  // edges incident to move candidates; avoids materializing the full graph).
  // Pairs WITHOUT a timing edge contribute affinity 0 by construction (sparse).
  std::unordered_map<unsigned, std::vector<std::pair<unsigned, double>>> adj;
  adj.reserve(critSet.size() * 2);
  parseCsv([&](unsigned a, unsigned b, double slack) {
    const bool ca = critSet.count(a) > 0;
    const bool cb = critSet.count(b) > 0;
    if (!ca && !cb) return;
    double aff = affinity(a, b, slack);
    if (!lpTarget.empty()) {
      aff *= std::exp(-std::abs(lpTarget[a] - lpTarget[b]) / compatTauNs);
    }
    if (ca) adj[a].emplace_back(b, aff);
    if (cb) adj[b].emplace_back(a, aff);
  });

  // Per-point tier lookup (targetOf/haveTargets hoisted above pass 2).
  auto tierOf = [&](unsigned idx) -> int {
    return cts3dDb_->getInstTier(sinkInsts[idx]->getDbInst());
  };

  // Worst violating slack first; each FF moved at most once (no oscillation).
  std::sort(violating.begin(), violating.end(),
            [](const VEdge& x, const VEdge& y) { return x.slack < y.slack; });

  std::unordered_set<unsigned> movedSet;
  long nMoves = 0, nCand = 0;

  auto tryMove = [&](unsigned f, unsigned partner) -> bool {
    if (movedSet.count(f)) return false;
    const int cs = clusterOf[f];
    const int cd = clusterOf[partner];
    if (cs < 0 || cd < 0 || cs == cd) return false;
    ++nCand;
    std::vector<unsigned>& dst = solution[cd];
    // (1) size cap: dest must stay below the seed's calibrated group size.
    if (dst.size() + 1 > maxClusterSize) return false;
    // (2) tier purity: every dest member must share f's tier, otherwise the
    //     move manufactures a mixed cluster -> downstream tier split fires.
    const int ft = tierOf(f);
    for (unsigned m : dst) {
      if (tierOf(m) != ft) return false;
    }
    // (3) diameter: dest bounding-box diagonal (normalized units, same
    //     convention as the calibration maxDiam) must stay within bound.
    double bxMin = points[f].first, bxMax = points[f].first;
    double byMin = points[f].second, byMax = points[f].second;
    for (unsigned m : dst) {
      bxMin = std::min(bxMin, (double) points[m].first);
      bxMax = std::max(bxMax, (double) points[m].first);
      byMin = std::min(byMin, (double) points[m].second);
      byMax = std::max(byMax, (double) points[m].second);
    }
    const double diag = std::sqrt((bxMax - bxMin) * (bxMax - bxMin)
                                  + (byMax - byMin) * (byMax - byMin));
    if (maxDiameter > 0 && diag > maxDiameter) return false;
    // (4) V37 target-spread guard: the move must not create a new 2/3-way
    //     split nor widen an existing over-threshold spread.
    if (enableTargetSplit_ && haveTargets) {
      double tMin = targetOf(f), tMax = tMin;
      double oMin = 1e30, oMax = -1e30;
      for (unsigned m : dst) {
        const double t = targetOf(m);
        oMin = std::min(oMin, t);
        oMax = std::max(oMax, t);
      }
      tMin = std::min(tMin, oMin);
      tMax = std::max(tMax, oMax);
      const double spreadBefore = (dst.empty() ? 0.0 : oMax - oMin);
      const double spreadAfter = tMax - tMin;
      if (spreadAfter > std::max(splitThreshold2wayNs_, spreadBefore)) {
        return false;
      }
    }
    // (5) affinity gain: sum of f's edge affinities into dest minus those
    //     staying behind in source must be positive.
    double gain = 0.0;
    auto it = adj.find(f);
    if (it == adj.end()) return false;
    for (const auto& [nbr, aff] : it->second) {
      if (clusterOf[nbr] == cd) gain += aff;
      else if (clusterOf[nbr] == cs && nbr != f) gain -= aff;
    }
    if (gain <= 0.0) return false;
    // Apply: erase f from source, append to dest. An emptied source cluster
    // is safe — the consumption loop's size==1 / size>1 branches both skip
    // size 0, so no leaf buffer and no dangling net is created (CTS-0726
    // class avoided by construction: nets do not exist yet at this point).
    std::vector<unsigned>& src = solution[cs];
    auto fit = std::find(src.begin(), src.end(), f);
    if (fit == src.end()) return false;  // defensive: clusterOf out of sync
    src.erase(fit);
    dst.push_back(f);
    clusterOf[f] = cd;
    movedSet.insert(f);
    ++nMoves;
    return true;
  };

  for (const VEdge& ve : violating) {
    if (maxMoves >= 0 && nMoves >= maxMoves) break;
    tryMove(ve.a, ve.b);
    if (maxMoves >= 0 && nMoves >= maxMoves) break;
    tryMove(ve.b, ve.a);
  }

  logger_->info(CTS, 846,
      "Seeded TARP: moved {} of {} candidate FFs ({} violating edges, "
      "{} critical FFs, {} clusters, maxSize={}, maxDiam={:.2f}, compatTau={})",
      nMoves, nCand, violating.size(), critSet.size(), solution.size(),
      maxClusterSize, maxDiameter, compatTauNs);
  const double refine_s = std::chrono::duration<double>(
      std::chrono::high_resolution_clock::now() - refine_t0).count();
  logger_->info(CTS, 867, "Profile: seeded refinement {:.2f}s", refine_s);
}

// TARP: Timing-affinity clustering integration.
void HTreeBuilder::runTarpClustering(
    const std::vector<std::pair<float, float>>& points,
    const std::vector<const ClockInst*>& sinkInsts,
    unsigned hilbertBestGroupSize,
    float hilbertBestDiameter)
{
  const unsigned N = points.size();

  TarpParams tp;
  tp.alpha = tarpAlpha_;
  tp.d_buf = perFfBufDelay_;
  tp.max_depth = groupedDelayMaxDepth_;

  // Hilbert-calibrated cluster size: use grid search result as default.
  // Env var CTS_TARP_MAX_CLUSTER overrides if set.
  tp.max_cluster_size = (hilbertBestGroupSize > 0) ? hilbertBestGroupSize : 10;
  if (const char* e = std::getenv("CTS_TARP_MAX_CLUSTER"))
    tp.max_cluster_size = std::atoi(e);

  if (const char* e = std::getenv("CTS_TARP_TOP_K"))
    tp.top_k_edges = std::atoi(e);
  if (const char* e = std::getenv("CTS_TARP_SLACK_THRESH"))
    tp.slack_threshold_ratio = std::stod(e);
  if (const char* e = std::getenv("CTS_TARP_FM_PASSES"))
    tp.fm_max_passes = std::atoi(e);

  // V85_TARP surgical mode: restrict FM moves to FFs on violating edges.
  if (const char* e = std::getenv("CTS_TARP_CRIT_ONLY"))
    tp.crit_only = (std::atoi(e) != 0);
  if (const char* e = std::getenv("CTS_TARP_CRIT_SLACK_NS"))
    tp.crit_slack_ns = std::stod(e);

  // Hilbert-calibrated diameter: use grid search result as default.
  // Env var CTS_TARP_MAX_DIAMETER_UM overrides if set.
  if (hilbertBestDiameter > 0) {
    tp.max_diameter = hilbertBestDiameter;  // already in wireSegmentUnit
  } else {
    tp.max_diameter = 50.0 * options_->getDbUnits() / wireSegmentUnit_;
  }
  if (const char* e = std::getenv("CTS_TARP_MAX_DIAMETER_UM"))
    tp.max_diameter = std::stod(e) * options_->getDbUnits() / wireSegmentUnit_;

  TarpClustering tarp(logger_, tp);

  // Add sinks with LP targets and hold budgets
  int hb_direct = 0, hb_downstream = 0, hb_fail_closed = 0;
  for (unsigned i = 0; i < N; ++i) {
    std::string nm = sinkInsts[i]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);

    double lp_target = cts3dDb_->hasSkewTargets()
                           ? cts3dDb_->getSkewTarget(nm) : 0.0;
    double hold_budget = 0.0;
    auto it = perFfHoldBudget_.find(nm);
    if (it != perFfHoldBudget_.end()) {
      hold_budget = it->second;
      ++hb_direct;
    } else {
      // Clone/gate sink: trace downstream to find real FFs' hold budgets
      odb::dbInst* dbInst = sinkInsts[i]->getDbInst();
      bool found_downstream = false;
      if (dbInst) {
        for (auto* iterm : dbInst->getITerms()) {
          if (iterm->getIoType() != odb::dbIoType::OUTPUT) continue;
          odb::dbNet* outNet = iterm->getNet();
          if (!outNet) continue;
          for (auto* sink_iterm : outNet->getITerms()) {
            if (sink_iterm == iterm) continue;
            std::string ff_name = sink_iterm->getInst()->getName();
            auto ff_it = perFfHoldBudget_.find(ff_name);
            if (ff_it != perFfHoldBudget_.end()) {
              if (!found_downstream) {
                hold_budget = ff_it->second;
                found_downstream = true;
              } else {
                hold_budget = std::min(hold_budget, ff_it->second);
              }
            }
          }
        }
      }
      if (!found_downstream) {
        ++hb_fail_closed;
      } else {
        ++hb_downstream;
      }
    }

    float cap = (sinkInsts[i]->getInputCap() == 0)
                    ? options_->getSinkBufferInputCap()
                    : sinkInsts[i]->getInputCap();

    tarp.addSink(i, points[i].first, points[i].second, cap,
                 lp_target, hold_budget, nm);
  }

  // Hold budget summary (replaces per-sink CTS-0431 spam)
  logger_->info(CTS, 431,
      "TARP hold budget: {} direct, {} downstream-resolved, {} fail-closed (budget=0)",
      hb_direct, hb_downstream, hb_fail_closed);

  // Read timing graph CSV for setup slacks
  auto prof_csv_t0 = std::chrono::high_resolution_clock::now();
  const char* csv_env = std::getenv("CTS_FF_TIMING_GRAPH_CSV");
  if (csv_env) {
    std::ifstream csvFile(csv_env);
    if (csvFile.is_open()) {
      std::unordered_map<std::string, unsigned> nameToIdx;
      const auto& sinks = tarp.sinks();
      for (unsigned i = 0; i < sinks.size(); ++i)
        nameToIdx[sinks[i].name] = i;

      // V85_TARP: default to the SDC-derived period (set by TritonCTS) instead
      // of the historical 1.0ns constant; env var still overrides.
      double T_period
          = (options_->getSdcPeriodNs() > 0.0) ? options_->getSdcPeriodNs() : 1.0;
      if (const char* e = std::getenv("CTS_CLOCK_PERIOD_NS"))
        T_period = std::stod(e);

      std::string line;
      std::getline(csvFile, line);  // skip header
      int edge_count = 0;
      while (std::getline(csvFile, line)) {
        std::vector<std::string> tokens;
        std::istringstream ls(line);
        std::string token;
        while (std::getline(ls, token, ',')) tokens.push_back(token);
        // CSV column layout (14 or 16 columns):
        //   0:from_ff, 1:to_ff, 2:slack_max, 3:slack_min, 4:arr_max, 5:arr_min,
        //   6:req_max, 7:req_min, 8:from_x, 9:from_y, 10:from_tier,
        //   11:to_x, 12:to_y, 13:to_tier [, 14:from_clock_net, 15:to_clock_net]
        // BUG FIX (V57_CLK2): from_clock_net/to_clock_net columns were appended later,
        // making the CSV 16 columns. The old "16-col legacy" branch (f_tier_idx=5,
        // t_tier_idx=8) incorrectly triggered, reading arrival_min_ns and from_x as tier
        // values — always mismatched — causing TARP to report 0 same-tier edges and pin
        // all FFs to Hilbert positions (timing-aware clustering completely disabled).
        // Fix: always use standard column positions (tier at col 10/13) for >=14 columns.
        if (tokens.size() < 14) continue;

        const int f_tier_idx = 10, t_tier_idx = 13;
        const int slack_max_idx = 2, slack_min_idx = 3;
        const int arr_max_idx = 4, arr_min_idx = 5;

        int f_tier = std::stoi(tokens[f_tier_idx]);
        int t_tier = std::stoi(tokens[t_tier_idx]);
        if (f_tier != t_tier) continue;  // skip cross-tier

        double slack_max = std::stod(tokens[slack_max_idx]);
        double slack_min = std::stod(tokens[slack_min_idx]);
        double arr_max = std::stod(tokens[arr_max_idx]);
        double arr_min = std::stod(tokens[arr_min_idx]);
        if (slack_max == 0 && slack_min == 0 && arr_max == 0 && arr_min == 0)
          continue;  // ghost edge

        auto it_a = nameToIdx.find(tokens[0]);
        auto it_b = nameToIdx.find(tokens[1]);
        if (it_a == nameToIdx.end() || it_b == nameToIdx.end()) continue;

        tarp.addTimingEdge(it_a->second, it_b->second, slack_max, T_period);
        ++edge_count;
      }
      logger_->info(CTS, 841, "TARP: {} same-tier edges from {}",
                    edge_count, csv_env);
    }
  }

  auto prof_csv_t1 = std::chrono::high_resolution_clock::now();
  double prof_csv_s = std::chrono::duration<double>(prof_csv_t1 - prof_csv_t0).count();
  logger_->info(CTS, 915, "Profile: TARP CSV parse {:.1f}s", prof_csv_s);

  auto prof_fm_t0 = std::chrono::high_resolution_clock::now();
  tarp.run();
  auto prof_fm_t1 = std::chrono::high_resolution_clock::now();
  double prof_fm_s = std::chrono::duration<double>(prof_fm_t1 - prof_fm_t0).count();
  logger_->info(CTS, 916, "Profile: TARP FM clustering {:.1f}s", prof_fm_s);

  // Convert clusters to leaf buffers
  std::vector<std::pair<float, float>> newSinkLocations;
  unsigned clusterCount = 0;

  for (const auto& cluster : tarp.clusters()) {
    if (cluster.empty()) continue;
    if (cluster.size() == 1) {
      unsigned orig = tarp.sinks()[cluster[0]].orig_idx;
      newSinkLocations.emplace_back(points[orig]);
      // Fix CTS-0492: register single-FF sink in mapLocationToSink_
      // (multi-FF clusters register at line 4725; single-FF was missing)
      const Point<double> sinkPt((double) points[orig].first,
                                 (double) points[orig].second);
      mapLocationToSink_[sinkPt] = const_cast<ClockInst*>(sinkInsts[orig]);
      ++clusterCount;
      continue;
    }

    std::vector<ClockInst*> clusterInsts;
    float xSum = 0, ySum = 0;
    for (unsigned si : cluster) {
      unsigned orig = tarp.sinks()[si].orig_idx;
      // Use sinkInsts directly — no map lookup needed, avoids float/double mismatch
      clusterInsts.push_back(const_cast<ClockInst*>(sinkInsts[orig]));
      xSum += points[orig].first;
      ySum += points[orig].second;
    }

    Point<double> center(xSum / cluster.size(), ySum / cluster.size());
    const int tier = cts3dDb_->getDominantTier(clusterInsts);
    const std::string buf = cts3dDb_->getBufferForTier(
        options_->getSinkBuffer(), tier);

    Point<double> loc = legalizeOneBuffer(center, buf);
    commitMoveLoc(center, loc);

    ClockInst& rootBuf = clock_.addClockBuffer(
        tierPrefix_ + "clkbuf_leaf_" + std::to_string(clusterCount), buf,
        loc.getX() * wireSegmentUnit_, loc.getY() * wireSegmentUnit_);
    cts3dDb_->setClockInstTier(rootBuf, tier);
    addFirstLevelSinkDriver(&rootBuf);

    ClockSubNet& net = clock_.addSubNet(
        tierPrefix_ + "clknet_leaf_" + std::to_string(clusterCount));
    net.addInst(rootBuf);
    for (ClockInst* ci : clusterInsts) net.addInst(*ci);
    net.setLeafLevel(true);

    const std::pair<float, float> bufPt(loc.getX(), loc.getY());
    newSinkLocations.emplace_back(bufPt);
    // V81 KEYFIX: same as createLeafBuffer.
    mapLocationToSink_[narrowedKey(loc.getX(), loc.getY())] = &rootBuf;

    if (cts3dDb_->hasSkewTargets()) {
      double sum = 0;
      for (unsigned si : cluster) sum += tarp.sinks()[si].lp_target;
      clusterBufTarget_[&rootBuf] = sum / cluster.size();
    }
    ++clusterCount;
  }

  topLevelSinksClustered_ = std::move(newSinkLocations);
  if (clusterCount) treeBufLevels_++;
  logger_->info(CTS, 844, "TARP: {} clusters, {} sinks",
                clusterCount, topLevelSinksClustered_.size());
}

// runTimingRefinement removed.
// FM min-cut with Hilbert initial partition replaces boundary swap.
// See TarpClustering.cpp for the FM implementation.

#if 0  // REMOVED: runTimingRefinement (replaced by TARP FM)
void HTreeBuilder::runTimingRefinement_REMOVED(
    std::vector<std::vector<unsigned>>& clusters,
    const std::vector<std::pair<float, float>>& points,
    const std::vector<const ClockInst*>& sinkInsts)
{
  const unsigned N = points.size();
  if (N == 0 || clusters.empty()) return;

  // Read params
  double slackThreshRatio = 0.3;
  if (const char* e = std::getenv("CTS_TARP_SLACK_THRESH"))
    slackThreshRatio = std::stod(e);

  // Get clock period (ns)
  double T_period = 1.0;
  if (const char* e = std::getenv("CTS_CLOCK_PERIOD_NS"))
    T_period = std::stod(e);
  const double slackThresh = slackThreshRatio * T_period;

  // Build FF name → point index map
  std::unordered_map<std::string, unsigned> nameToIdx;
  for (unsigned i = 0; i < N; ++i) {
    std::string nm = sinkInsts[i]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);
    nameToIdx[nm] = i;
  }

  // Build point index → cluster index map
  std::vector<int> ffCluster(N, -1);
  for (unsigned c = 0; c < clusters.size(); ++c) {
    for (unsigned idx : clusters[c]) {
      if (idx < N) ffCluster[idx] = static_cast<int>(c);
    }
  }

  // Load critical timing edges from CSV
  const char* csv_env = std::getenv("PRE_CTS_FF_TIMING_GRAPH");
  if (!csv_env) {
    logger_->info(CTS, 875, "TARP refinement: no timing graph CSV, skipping");
    return;
  }

  struct CritEdge {
    unsigned ff_a, ff_b;
    double slack;
  };
  std::vector<CritEdge> critEdges;

  std::ifstream csvFile(csv_env);
  if (!csvFile.is_open()) {
    logger_->warn(CTS, 876, "TARP refinement: cannot open {}", csv_env);
    return;
  }

  std::string line;
  std::getline(csvFile, line);  // skip header
  while (std::getline(csvFile, line)) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(ss, token, ',')) tokens.push_back(token);
    if (tokens.size() < 14) continue;

    const std::string& from_ff = tokens[0];
    const std::string& to_ff = tokens[1];
    double slack_max = std::stod(tokens[2]);
    int from_tier = std::stoi(tokens[10]);
    int to_tier = std::stoi(tokens[13]);

    // Same-tier only
    if (from_tier != to_tier) continue;
    // Ghost edge filter
    double slack_min = std::stod(tokens[3]);
    double arr_max = std::stod(tokens[4]);
    double arr_min = std::stod(tokens[5]);
    if (slack_max == 0 && slack_min == 0 && arr_max == 0 && arr_min == 0) continue;
    // Slack threshold filter
    if (slack_max > slackThresh) continue;

    auto it_a = nameToIdx.find(from_ff);
    auto it_b = nameToIdx.find(to_ff);
    if (it_a == nameToIdx.end() || it_b == nameToIdx.end()) continue;

    critEdges.push_back({it_a->second, it_b->second, slack_max});
  }

  // Count cross-cluster critical pairs
  unsigned crossBefore = 0;
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) crossBefore++;
  }

  logger_->info(CTS, 877,
      "TARP refinement: {} critical edges (slack<{:.3f}ns), {} cross-cluster ({:.1f}%)",
      critEdges.size(), slackThresh,
      crossBefore, critEdges.empty() ? 0.0 : 100.0 * crossBefore / critEdges.size());

  if (crossBefore == 0) {
    logger_->info(CTS, 878, "TARP refinement: no cross-cluster critical pairs, skipping");
    return;
  }

  // Compute max diameter from best clustering (approximate from existing clusters)
  double maxDiameter = 0;
  for (const auto& cl : clusters) {
    if (cl.size() <= 1) continue;
    float xn = points[cl[0]].first, xx = xn;
    float yn = points[cl[0]].second, yx = yn;
    for (unsigned idx : cl) {
      xn = std::min(xn, points[idx].first); xx = std::max(xx, points[idx].first);
      yn = std::min(yn, points[idx].second); yx = std::max(yx, points[idx].second);
    }
    double diam = (xx - xn) + (yx - yn);
    maxDiameter = std::max(maxDiameter, diam);
  }
  // Allow 20% diameter expansion for refinement
  maxDiameter *= 1.2;

  // Build per-FF critical edge count (how many critical edges this FF participates in)
  std::vector<unsigned> ffCritCount(N, 0);
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) {
      ffCritCount[e.ff_a]++;
      ffCritCount[e.ff_b]++;
    }
  }

  // Greedy swap: for each cross-cluster critical edge, try moving the FF
  // with fewer critical edges to the other's cluster
  unsigned swaps = 0;
  const unsigned maxSwaps = N / 4;  // cap at 25% of FFs

  // Sort critical edges by slack (tightest first — most important to fix)
  auto sortedEdges = critEdges;
  std::sort(sortedEdges.begin(), sortedEdges.end(),
            [](const CritEdge& a, const CritEdge& b) { return a.slack < b.slack; });

  for (const auto& e : sortedEdges) {
    if (swaps >= maxSwaps) break;

    int ca = ffCluster[e.ff_a];
    int cb = ffCluster[e.ff_b];
    if (ca == cb) continue;  // already co-clustered (maybe by previous swap)

    // Pick the FF to move: the one with fewer critical edges in its current cluster
    // (less disruptive to move)
    unsigned mover, target_cluster;
    if (ffCritCount[e.ff_a] <= ffCritCount[e.ff_b]) {
      mover = e.ff_a;
      target_cluster = cb;
    } else {
      mover = e.ff_b;
      target_cluster = ca;
    }

    int src_cluster = ffCluster[mover];

    // Don't empty a cluster
    if (clusters[src_cluster].size() <= 2) continue;

    // Check diameter constraint after swap
    float xn = points[mover].first, xx = xn;
    float yn = points[mover].second, yx = yn;
    for (unsigned idx : clusters[target_cluster]) {
      xn = std::min(xn, points[idx].first); xx = std::max(xx, points[idx].first);
      yn = std::min(yn, points[idx].second); yx = std::max(yx, points[idx].second);
    }
    double newDiam = (xx - xn) + (yx - yn);
    if (newDiam > maxDiameter) continue;  // would violate spatial constraint

    // Execute swap
    // Remove from source cluster
    auto& srcCluster = clusters[src_cluster];
    srcCluster.erase(std::remove(srcCluster.begin(), srcCluster.end(), mover),
                     srcCluster.end());
    // Add to target cluster
    clusters[target_cluster].push_back(mover);
    ffCluster[mover] = target_cluster;
    swaps++;
  }

  // Count cross-cluster after refinement
  unsigned crossAfter = 0;
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) crossAfter++;
  }

  logger_->info(CTS, 879,
      "TARP refinement: {} swaps, cross-cluster {} -> {} ({:.1f}% reduction)",
      swaps, crossBefore, crossAfter,
      crossBefore > 0 ? 100.0 * (crossBefore - crossAfter) / crossBefore : 0.0);
}
#endif  // REMOVED: runTimingRefinement

}  // namespace dcts
