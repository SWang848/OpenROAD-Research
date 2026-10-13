// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

// clang-format off
%{

#include "tdp/MakeTimingdp.h"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"

%}

%include "../../Exception.i"

%import <std_vector.i>
%import "dbtypes.i"

%inline %{

namespace tdp {

void
timing_detailed_placement_cmd(int max_displacment_x,
                       int max_displacment_y,
                       double m,
                       const char* report_file_name,
                       const char* tech,
                       const char* hbt){
  tdp::Timingdp *timingdp = tdp::getTimingdp();
  timingdp->detailedPlacement(max_displacment_x, max_displacment_y,
                              std::string(report_file_name), m,
                              std::string(tech), std::string(hbt));
}

void
report_timing_stats_cmd(int max_displacment_x,
                        int max_displacment_y,
                        double m,
                        const char* tech,
                        const char* hbt){
  tdp::Timingdp *timingdp = tdp::getTimingdp();
  timingdp->reportTimingStats(max_displacment_x, max_displacment_y, m,
                              std::string(tech), std::string(hbt));
}

const char*
get_tier_site(const char* tier, const char* tech){
  return tdp::Timingdp::tierSite(std::string(tier), std::string(tech));
}

void
report_legalization_stats()
{
  tdp::Timingdp *timingdp = tdp::getTimingdp();
  timingdp->reportLegalizationStats();
}

void
create_new_row(const char* layer)
{
  tdp::Timingdp *timingdp = tdp::getTimingdp();
  timingdp->createnewrow(std::string(layer));
}

void
design_set()
{
  tdp::Timingdp *timingdp = tdp::getTimingdp();
  timingdp->designset();
}

// void
// check_placement_cmd(bool verbose, const char* report_file_name)
// {
//   dpl::Opendp *opendp = ord::OpenRoad::openRoad()->getOpendp();
//   opendp->checkPlacement(verbose, std::string(report_file_name));
// }

// void
// set_debug_cmd(float min_displacement,
//               const odb::dbInst* debug_instance,
//               int jump_moves,
//               bool iterative_placement)
// {
//   dpl::Opendp* opendp = ord::OpenRoad::openRoad()->getOpendp();
//   opendp->setJumpMoves(jump_moves);
//   opendp->setIterativePlacement(iterative_placement);
//   if (dpl::Graphics::guiActive()) {
//       std::unique_ptr<DplObserver> graphics = std::make_unique<dpl::Graphics>(
//           opendp, min_displacement, debug_instance);
//       opendp->setDebug(graphics);
//   }
// }

} // namespace

%} // inline
