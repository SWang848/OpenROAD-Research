// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "LatencyBalancer.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <ranges>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "Clock.h"
#include "CtsOptions.h"
#include "TreeBuilder.h"
#include "Util.h"
#include "dcts/TritonCTS.h"
#include "odb/db.h"
#include "odb/dbObject.h"
#include "odb/dbSet.h"
#include "odb/geom.h"
#include "sta/Clock.hh"
#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/Mode.hh"
#include "sta/NetworkClass.hh"
#include "sta/Path.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/Sdc.hh"
#include "sta/TimingArc.hh"
#include "sta/TimingModel.hh"
#include "utl/Logger.h"

namespace dcts {

using utl::CTS;

int LatencyBalancer::run()
{
  logger_->info(CTS,
                33,
                "Balancing latency for clock {}",
                root_->getClock().getSdcName());
  // B' skew-aware balancing (CTS_BALANCER_SKEW_AWARE=1, default off): pad each
  // leaf cluster to worseDelay_ + mean(a_i) instead of the uniform worseDelay_.
  // Requires setSkewTargets() to have been called with the LP arrival map.
  {
    const char* e = std::getenv("CTS_BALANCER_SKEW_AWARE");
    skewAware_ = (e && std::atoi(e) != 0 && skewTargets_ != nullptr
                  && !skewTargets_->empty());
    if (skewAware_) {
      logger_->info(CTS, 895,
                    "B' skew-aware balancing ON: {} FF skew targets loaded",
                    static_cast<int>(skewTargets_->size()));
    }
  }
  initSta();
  findLeafBuilders(root_);
  buildGraph(root_->getTopInputNet());
  bufferDelay_ = computeBufferDelay(0);
  balanceLatencies(0);
  logger_->info(CTS,
                36,
                " inserted {} delay buffers",
                delayBufIndex_,
                root_->getClock().getSdcName());
  return delayBufIndex_;
}

void LatencyBalancer::initSta()
{
  openSta_->ensureGraph();
  for (auto mode : openSta_->modes()) {
    openSta_->ensureClkNetwork(mode);
  }
  openSta_->updateTiming(false);
  timingGraph_ = openSta_->graph();
}

// This SHOULD return float. It would need to call LumpedCapDelayCalc::gateDelay
// instead of accessing the models directly to get POCV parameters.
sta::ArcDelay LatencyBalancer::computeBufferDelay(double extra_out_cap)
{
  odb::dbMaster* bufferMaster
      = db_->findMaster(options_->getRootBuffer().c_str());
  sta::Cell* bufferMasterCell = network_->dbToSta(bufferMaster);
  sta::LibertyCell* buffer_cell = network_->libertyCell(bufferMasterCell);
  sta::ArcDelay max_rise_delay = 0;

  sta::LibertyPort *input, *output;
  buffer_cell->bufferPorts(input, output);
  for (sta::Scene* corner : openSta_->scenes()) {
    const sta::Pvt* pvt
        = openSta_->cmdMode()->sdc()->operatingConditions(sta::MinMax::max());

    for (sta::TimingArcSet* arc_set :
         buffer_cell->sceneCell(corner, sta::MinMax::max())
             ->timingArcSets(input, output)) {
      for (sta::TimingArc* arc : arc_set->arcs()) {
        sta::GateTimingModel* model
            = dynamic_cast<sta::GateTimingModel*>(arc->model());
        const sta::RiseFall* in_rf = arc->fromEdge()->asRiseFall();
        const sta::RiseFall* out_rf = arc->toEdge()->asRiseFall();
        // Only look at rise-rise arcs
        if (model != nullptr && in_rf == sta::RiseFall::rise()
            && out_rf == sta::RiseFall::rise()) {
          double in_cap = input->capacitance(in_rf, sta::MinMax::max());
          double load_cap = in_cap + extra_out_cap;
          float arc_delay, arc_slew;
          model->gateDelay(pvt, 0.0, load_cap, arc_delay, arc_slew);
          // Cycle the arc_slew through the gate delay calculator once more
          model->gateDelay(pvt, arc_slew, load_cap, arc_delay, arc_slew);
          // and once more
          model->gateDelay(pvt, arc_slew, load_cap, arc_delay, arc_slew);

          if (delayGreater(arc_delay, max_rise_delay, openSta_)) {
            max_rise_delay = arc_delay;
          }
        }
      }
    }
  }

  return max_rise_delay;
}

void LatencyBalancer::findLeafBuilders(TreeBuilder* builder)
{
  if (builder->isLeafTree()) {
    std::string topBufferName = builder->getTopBufferName();
    inst2builder_[topBufferName] = builder;
  }
  for (const auto& child : builder->getChildren()) {
    findLeafBuilders(child);
  }
}

void LatencyBalancer::buildGraph(odb::dbNet* clkInputNet)
{
  // Macro-skip fix (CTS_BALANCER_SKIP_MACRO=1): exclude macro (block) sinks and
  // macro-tree builders from worseDelay_ so registers are balanced only to the
  // register worst, not padded up to the macro's large internal clkTreeDelay.
  // Default (env unset) = original behavior (macros raise worseDelay_).
  static int skipMacro = -1;
  if (skipMacro < 0) {
    const char* e = std::getenv("CTS_BALANCER_SKIP_MACRO");
    skipMacro = (e && std::atoi(e) != 0) ? 1 : 0;
  }

  std::string rootSrcName;
  odb::dbITerm* rootOutputITerm = clkInputNet->getFirstOutput();
  if (!rootOutputITerm) {
    rootSrcName = clkInputNet->getName();
  } else {
    rootSrcName = rootOutputITerm->getInst()->getName();
  }
  int builderSrcId = graph_.size();
  GraphNode builderSrcNode
      = GraphNode(builderSrcId, std::move(rootSrcName), rootOutputITerm);
  graph_.push_back(std::move(builderSrcNode));

  std::stack<int> visitNode;
  visitNode.push(builderSrcId);

  while (!visitNode.empty()) {
    int driverId = visitNode.top();
    visitNode.pop();

    odb::dbNet* driverNet;
    if (graph_[driverId].inputTerm != nullptr) {
      odb::dbInst* driverInst = graph_[driverId].inputTerm->getInst();
      odb::dbITerm* firstOutput = driverInst->getFirstOutput();
      driverNet = firstOutput->getNet();
      if (!driverNet) {
        continue;
      }
    } else {
      driverNet = clkInputNet;
    }

    for (odb::dbITerm* sinkIterm : driverNet->getITerms()) {
      if (sinkIterm->getIoType() == odb::dbIoType::INPUT) {
        if (!isSink(sinkIterm) && !propagateClock(sinkIterm)) {
          continue;
        }
        int sinkId = graph_.size();
        odb::dbInst* sinkInst = sinkIterm->getInst();
        std::string sinkName = sinkInst->getName();
        GraphNode sinkNode = GraphNode(sinkId, sinkName, sinkIterm);
        graph_.push_back(std::move(sinkNode));
        graph_[driverId].childrenIds.push_back(sinkId);

        // V60_BalBase: CTS_BALANCER_FULL_TRAVERSE=1 makes the balancer traverse
        // INTO leaf builders (balance every intra-builder sink branch) instead
        // of black-boxing the whole sub-tree with a single AVERAGE arrival.
        // Root cause fixed: averaging hid the per-tier depth spread inside a
        // gated-clock register tree (bottom-gated 75ps vs upper-gated 141ps),
        // so balanceLatencies had no sibling nodes to equalize. Full traversal
        // exposes those branch points → balancer can coordinate them.
        // Default (env unset) = original average black-box behavior.
        if (inst2builder_.find(sinkName) != inst2builder_.end()) {
          static int fullTraverse = -1;
          if (fullTraverse < 0) {
            const char* e = std::getenv("CTS_BALANCER_FULL_TRAVERSE");
            fullTraverse = (e && std::atoi(e) != 0) ? 1 : 0;
          }
          if (!fullTraverse) {
            auto builder = inst2builder_[sinkName];
            float builerAvgArrival = computeAveSinkArrivals(builder);
            const bool macroBuilder
                = (builder->getTreeType() == TreeType::MacroTree);
            graph_[sinkId].isMacro = macroBuilder;
            // Macro-skip: a macro-tree builder must not raise worseDelay_.
            if (!(skipMacro && macroBuilder)) {
              worseDelay_ = std::max(worseDelay_, builerAvgArrival);
            }
            graph_[sinkId].arrival = builerAvgArrival;
            // B' skew-aware: attach this cluster's mean useful-skew target so
            // computeNumberOfDelayBuffers pads it to worseDelay_ + skewTarget.
            if (skewAware_ && !(skipMacro && macroBuilder)) {
              graph_[sinkId].skewTarget = computeAveSkewTarget(builder);
            }
            continue;
          }
          // else: fall through → treat top buffer as a normal internal node →
          // push to visit stack → traverse into the sub-tree.
        }

        if (isSink(sinkIterm)) {
          sta::Pin* pin = network_->dbToSta(sinkIterm);
          if (pin) {
            sta::Vertex* sinkVertex = timingGraph_->pinDrvrVertex(pin);
            float arrival
                = getVertexClkArrival(sinkVertex, clkInputNet, sinkIterm);
            float insDelay = 0.0;
            sta::LibertyCell* libCell
                = network_->libertyCell(network_->dbToSta(sinkInst));
            odb::dbMTerm* mterm = sinkIterm->getMTerm();
            if (libCell && mterm) {
              sta::LibertyPort* libPort
                  = libCell->findLibertyPort(mterm->getConstName());
              if (libPort) {
                const float rise = libPort->clkTreeDelay(
                    0.0, sta::RiseFall::rise(), sta::MinMax::max());
                const float fall = libPort->clkTreeDelay(
                    0.0, sta::RiseFall::fall(), sta::MinMax::max());

                if (rise != 0 || fall != 0) {
                  insDelay = (rise + fall);
                  if (rise != 0 && fall != 0) {
                    insDelay /= 2.0;
                  }
                }
              }
            }
            const bool macroSink = sinkInst->isBlock();
            graph_[sinkId].isMacro = macroSink;
            // Macro-skip: a macro sink must not raise worseDelay_ (its internal
            // clkTreeDelay would drag every register up → reg<->macro hold blowup).
            if (!(skipMacro && macroSink)) {
              worseDelay_ = std::max(worseDelay_, (arrival + insDelay));
            }
            graph_[sinkId].arrival = arrival + insDelay;
            // B' skew-aware: an individual FF sink gets its OWN LP target a_i
            // (full per-FF granularity), converted ns -> seconds. This is where
            // ng45/ibex register sinks land (leaf builders are black-boxed
            // separately). 0 if unmapped -> falls back to uniform worseDelay_.
            if (skewAware_ && skewTargets_ != nullptr
                && !(skipMacro && macroSink)) {
              auto it = skewTargets_->find(sinkInst->getName());
              if (it != skewTargets_->end()) {
                graph_[sinkId].skewTarget
                    = static_cast<float>(it->second * 1e-9);
              }
            }
            debugPrint(logger_,
                       CTS,
                       "insertion delay",
                       2,
                       "Sink {}: average sink arrival is {:0.3e}",
                       sinkIterm->getName(),
                       arrival);
          }
          continue;
        }
        visitNode.push(sinkId);
      }
    }
  }
}

odb::dbITerm* LatencyBalancer::getFirstInput(odb::dbInst* inst) const
{
  odb::dbSet<odb::dbITerm> iterms = inst->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    if (iterm->isInputSignal()) {
      return iterm;
    }
  }

  return nullptr;
}

float LatencyBalancer::getVertexClkArrival(sta::Vertex* sinkVertex,
                                           odb::dbNet* topNet,
                                           odb::dbITerm* iterm)
{
  sta::VertexPathIterator pathIter(sinkVertex, openSta_);
  float clkPathArrival = 0.0;
  while (pathIter.hasNext()) {
    sta::Path* path = pathIter.next();
    const sta::ClockEdge* clockEdge = path->clkEdge(openSta_);
    if (clockEdge == nullptr) {
      continue;
    }

    if (clockEdge->transition() != sta::RiseFall::rise()) {
      // only populate with rising edges
      continue;
    }

    if (path->minMax(openSta_) != sta::MinMax::max()) {
      continue;
      // only populate with max delay
    }

    const sta::Clock* clock = path->clock(openSta_);
    if (clock) {
      sta::PathExpanded expand(path, openSta_);
      const sta::Path* start = expand.startPath();

      odb::dbNet* pathStartNet = nullptr;

      odb::dbITerm* term;
      odb::dbBTerm* port;
      odb::dbModITerm* modIterm;
      network_->staToDb(start->pin(openSta_), term, port, modIterm);
      if (term) {
        pathStartNet = term->getNet();
      }
      if (port) {
        pathStartNet = port->getNet();
      }
      if (pathStartNet == topNet) {
        clkPathArrival = path->arrival();
        return clkPathArrival;
      }
    }
  }
  logger_->warn(CTS, 179, "No paths found for pin {}.", iterm->getName());
  return clkPathArrival;
}

float LatencyBalancer::computeAveSinkArrivals(TreeBuilder* builder)
{
  Clock clock = builder->getClock();
  odb::dbNet* topInputClockNet = builder->getTopInputNet();
  // compute average input arrival at all sinks
  float sumArrivals = 0.0;
  unsigned numSinks = 0;
  clock.forEachSink([&](const ClockInst& sink) {
    odb::dbITerm* iterm = sink.getDbInputPin();
    computeSinkArrivalRecur(topInputClockNet, iterm, sumArrivals, numSinks);
  });
  float aveArrival = 0.0;
  if (numSinks) {
    aveArrival = sumArrivals / (float) numSinks;
  }
  builder->setAveSinkArrival(aveArrival);
  debugPrint(logger_,
             CTS,
             "insertion delay",
             2,
             "{} {}: average sink arrival is {:0.3e}",
             (builder->getTreeType() == TreeType::MacroTree) ? "macro tree"
                                                             : "register tree",
             clock.getName(),
             builder->getAveSinkArrival());

  return aveArrival;
}

// B' skew-aware balancing: mean LP arrival target a_i (converted ns -> seconds)
// over this builder's FF sinks. Returns 0 when disabled or no sink is in the
// skew-target map, so the node falls back to the uniform worseDelay_ target.
float LatencyBalancer::computeAveSkewTarget(TreeBuilder* builder) const
{
  if (!skewAware_ || skewTargets_ == nullptr) {
    return 0.0f;
  }
  Clock clock = builder->getClock();
  double sum = 0.0;
  unsigned n = 0;
  clock.forEachSink([&](const ClockInst& sink) {
    // Strip a trailing "/CLK"-style pin suffix to match the FF instance key.
    std::string name = sink.getName();
    const auto slash = name.rfind('/');
    if (slash != std::string::npos) {
      name = name.substr(0, slash);
    }
    auto it = skewTargets_->find(name);
    if (it != skewTargets_->end()) {
      sum += it->second;  // ns
      ++n;
    }
  });
  if (n == 0) {
    return 0.0f;
  }
  // ns -> seconds to match arrival / worseDelay_ / bufferDelay_ units.
  return static_cast<float>((sum / static_cast<double>(n)) * 1e-9);
}

void LatencyBalancer::computeSinkArrivalRecur(odb::dbNet* topClokcNet,
                                              odb::dbITerm* iterm,
                                              float& sumArrivals,
                                              unsigned& numSinks)
{
  if (iterm) {
    odb::dbInst* inst = iterm->getInst();
    if (inst) {
      if (isSink(iterm)) {
        // either register or macro input pin
        sta::Pin* pin = network_->dbToSta(iterm);
        if (pin) {
          sta::Vertex* sinkVertex = timingGraph_->pinDrvrVertex(pin);
          float arrival = getVertexClkArrival(sinkVertex, topClokcNet, iterm);
          // add insertion delay
          float insDelay = 0.0;
          sta::LibertyCell* libCell
              = network_->libertyCell(network_->dbToSta(inst));
          odb::dbMTerm* mterm = iterm->getMTerm();
          if (libCell && mterm) {
            sta::LibertyPort* libPort
                = libCell->findLibertyPort(mterm->getConstName());
            if (libPort) {
              const float rise = libPort->clkTreeDelay(
                  0.0, sta::RiseFall::rise(), sta::MinMax::max());
              const float fall = libPort->clkTreeDelay(
                  0.0, sta::RiseFall::fall(), sta::MinMax::max());

              if (rise != 0 || fall != 0) {
                insDelay = (rise + fall);
                if (rise != 0 && fall != 0) {
                  insDelay /= 2.0;
                }
              }
            }
          }
          sumArrivals += (arrival + insDelay);
          numSinks++;
        }
      } else {
        // not a sink, but a clock gater
        odb::dbITerm* outTerm = inst->getFirstOutput();
        if (outTerm) {
          odb::dbNet* outNet = outTerm->getNet();
          bool propagate = propagateClock(iterm);
          if (outNet && propagate) {
            odb::dbSet<odb::dbITerm> iterms = outNet->getITerms();
            odb::dbSet<odb::dbITerm>::iterator iter;
            for (iter = iterms.begin(); iter != iterms.end(); ++iter) {
              odb::dbITerm* inTerm = *iter;
              if (inTerm->getIoType() == odb::dbIoType::INPUT) {
                computeSinkArrivalRecur(
                    topClokcNet, inTerm, sumArrivals, numSinks);
              }
            }
          }
        }
      }
    }
  }
}

void LatencyBalancer::computeNumberOfDelayBuffers(int nodeId,
                                                  int srcX,
                                                  int srcY)
{
  GraphNode* node = &graph_[nodeId];
  // Macro-skip: never insert delay buffers to pad a macro node up to
  // worseDelay_ (leave macros at their natural latency).
  static int skipMacro = -1;
  if (skipMacro < 0) {
    const char* e = std::getenv("CTS_BALANCER_SKIP_MACRO");
    skipMacro = (e && std::atoi(e) != 0) ? 1 : 0;
  }
  if (skipMacro && node->isMacro) {
    node->nBuffInsert = 0;
    return;
  }
  if (node->arrival != 0.0) {
    // B' skew-aware: target = worseDelay_ + skewTarget (extra insertion delay =
    // the cluster's useful skew). skewTarget = 0 (default) => uniform worseDelay_.
    const double target = worseDelay_ + node->skewTarget;
    int numBuffers = (int) ((target - node->arrival) / bufferDelay_);

    // adjust buffer delay for wire cap
    int sinkX, sinkY;
    graph_[nodeId].inputTerm->getAvgXY(&sinkX, &sinkY);
    float offsetX = (float) (sinkX - srcX) / (numBuffers + 1);
    float offsetY = (float) (sinkY - srcY) / (numBuffers + 1);
    auto newDelay = computeBufferDelay((std::abs(offsetX) + std::abs(offsetY))
                                       * capPerDBU_);
    numBuffers = (int) ((target - node->arrival) / newDelay);
    if (node->childrenIds.empty()) {
      debugPrint(logger_,
                 CTS,
                 "insertion delay",
                 3,
                 "For node {}, isert {:2f} buffers",
                 node->name,
                 numBuffers);
    }
    node->nBuffInsert = numBuffers;
  }
}

void LatencyBalancer::balanceLatencies(int nodeId)
{
  GraphNode* node = &graph_[nodeId];

  // Compute number of buffer needed for leaf node
  if (node->childrenIds.empty()) {
    return;
  }

  // If it is not a leaf node compute the amount of buffers needed for its
  // children
  std::vector<odb::dbITerm*> sinksInput;
  int previouBufToInsert = 0;
  int srcX, srcY;
  if (node->inputTerm == nullptr) {
    odb::dbNet* rootNet = root_->getTopInputNet();
    odb::dbBTerm* clkInput = rootNet->get1stBTerm();
    odb::Rect clkInputBBox = clkInput->getBBox();
    srcX = clkInputBBox.xCenter();
    srcY = clkInputBBox.yCenter();
  } else {
    node->inputTerm->getAvgXY(&srcX, &srcY);
  }

  double maxArrival = std::numeric_limits<double>::min();
  std::map<int, std::vector<odb::dbITerm*>> buffersNeeded2Childern;
  for (int child : node->childrenIds) {
    balanceLatencies(child);
    computeNumberOfDelayBuffers(child, srcX, srcY);
    maxArrival = std::max(graph_[child].arrival, maxArrival);
    buffersNeeded2Childern[graph_[child].nBuffInsert].push_back(
        graph_[child].inputTerm);
  }

  // If the children need a different amount of buffers insert this difference
  for (auto& [bufToInsert, children] :
       std::ranges::reverse_view(buffersNeeded2Childern)) {
    if (bufToInsert == -1) {
      continue;
    }

    if (!previouBufToInsert) {
      previouBufToInsert = bufToInsert;
      sinksInput.clear();
      sinksInput = std::move(children);
      continue;
    }

    int numBuffers = previouBufToInsert - bufToInsert;
    odb::dbITerm* delauBuffInput
        = insertDelayBuffers(numBuffers, srcX, srcY, sinksInput);

    sinksInput.clear();
    sinksInput = std::move(children);
    sinksInput.push_back(delauBuffInput);

    previouBufToInsert = bufToInsert;
  }

  node->nBuffInsert = previouBufToInsert;
  node->arrival = maxArrival;
}

odb::dbITerm* LatencyBalancer::insertDelayBuffers(
    int numBuffers,
    int srcX,
    int srcY,
    const std::vector<odb::dbITerm*>& sinksInput)
{
  // get bbox of current load pins without driver output pin
  odb::Rect loadPinsBbox = odb::Rect();
  loadPinsBbox.mergeInit();
  odb::dbNet* drivingNet = nullptr;
  debugPrint(logger_,
             CTS,
             "insertion delay",
             3,
             "Inserting {} buffers for sinks:",
             numBuffers);
  for (odb::dbITerm* sinkInput : sinksInput) {
    if (drivingNet == nullptr) {
      drivingNet = sinkInput->getNet();
    }
    debugPrint(
        logger_, CTS, "insertion delay", 3, "  {}", sinkInput->getName());
    sinkInput->disconnect();
    int sinkX, sinkY;
    sinkInput->getAvgXY(&sinkX, &sinkY);
    loadPinsBbox.merge({sinkX, sinkY});
  }

  float offsetX = (float) (loadPinsBbox.xCenter() - srcX) / (numBuffers + 1);
  float offsetY = (float) (loadPinsBbox.yCenter() - srcY) / (numBuffers + 1);

  odb::dbInst* returnBuffer = nullptr;

  for (int i = 0; i < numBuffers; i++) {
    // Set the location
    double locX = (double) (srcX + (offsetX * (i + 1))) / wireSegmentUnit_;
    double locY = (double) (srcY + (offsetY * (i + 1))) / wireSegmentUnit_;
    Point<double> bufferLoc(locX, locY);
    Point<double> legalBufferLoc
        = root_->legalizeOneBuffer(bufferLoc, options_->getRootBuffer());

    odb::Point loc{static_cast<int>(legalBufferLoc.getX() * wireSegmentUnit_),
                   static_cast<int>(legalBufferLoc.getY() * wireSegmentUnit_)};

    // Insert buffer
    std::string clkName = root_->getClock().getSdcName();
    std::string newNetName
        = fmt::format("delaynet_{}_{}", delayBufIndex_, clkName);
    std::string newBufferName
        = fmt::format("delaybuf_{}_{}", delayBufIndex_++, clkName);
    odb::dbMaster* bufferMaster
        = db_->findMaster(options_->getRootBuffer().c_str());

    odb::dbInst* lastBuffer = nullptr;

    // Use load pins buffering at the end
    std::set<odb::dbObject*> load_pins;
    for (odb::dbITerm* sinkInput : sinksInput) {
      load_pins.insert(sinkInput);
    }

    // load_pins are not connected yet. So this option is required.
    bool loads_on_different_nets = true;
    lastBuffer = drivingNet->insertBufferBeforeLoads(
        load_pins,
        bufferMaster,
        &loc,
        newBufferName.c_str(),
        newNetName.c_str(),
        odb::dbNameUniquifyType::IF_NEEDED,
        loads_on_different_nets);

    debugPrint(logger_,
               CTS,
               "insertion delay",
               1,
               "new delay buffer {} is inserted at ({} {})",
               lastBuffer->getName(),
               loc.getX(),
               loc.getY());

    // Update the driving iterm & net to insert a next buffer on it
    odb::dbITerm* drvrPin = lastBuffer->getFirstOutput();
    drivingNet = drvrPin->getNet();

    // Update return buffer (the first buffer inserted)
    if (returnBuffer == nullptr) {
      returnBuffer = lastBuffer;
    }
  }

  return getFirstInput(returnBuffer);
}

bool LatencyBalancer::propagateClock(odb::dbITerm* input)
{
  odb::dbInst* inst = input->getInst();
  sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
  sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);

  if (!libertyCell) {
    return false;
  }
  // Clock tree buffers
  if (libertyCell->isInverter() || libertyCell->isBuffer()) {
    return true;
  }
  // Combinational components
  if (!libertyCell->hasSequentials()) {
    return true;
  }
  sta::LibertyPort* inputPort
      = libertyCell->findLibertyPort(input->getMTerm()->getConstName());

  // Clock Gater / Latch improvised as clock gater
  if (inputPort) {
    return inputPort->isClockGateClock() || inputPort->isLatchData();
  }

  return false;
}

bool LatencyBalancer::isSink(odb::dbITerm* iterm)
{
  odb::dbInst* inst = iterm->getInst();
  sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
  sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);
  if (!libertyCell) {
    return true;
  }

  if (inst->isBlock()) {
    return true;
  }

  sta::LibertyPort* inputPort
      = libertyCell->findLibertyPort(iterm->getMTerm()->getConstName());
  if (inputPort) {
    return inputPort->isRegClk();
  }

  return false;
}

void LatencyBalancer::showGraph()
{
  logger_->report("Graph built:");
  for (const auto& node : graph_) {
    odb::dbITerm* inputTerm = node.inputTerm;
    logger_->report(" Node {}", node.name);
    logger_->report("   id       = {}", node.id);
    logger_->report("   delay    = {}", node.arrival);
    logger_->report("   n buffer = {}", node.nBuffInsert);
    logger_->report("   in Term  = {}",
                    inputTerm == nullptr ? "no dbITerm" : inputTerm->getName());
    logger_->report("   childern [");
    for (int cId : node.childrenIds) {
      logger_->report("              {},", cId);
    }
    logger_->report("            ]");
  }
}

}  // namespace dcts
