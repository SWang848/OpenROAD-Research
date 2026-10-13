// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

%{
#include <cstdint>

#include "dcts/TritonCTS.h"
#include "CtsOptions.h"
#include "TechChar.h"
#include "CtsGraphics.h"
#include "dcts/Make3DCTS.h"

using namespace dcts;
%}

%include "../../Exception.i"
%include "stdint.i"
%include "std_string.i"

%ignore dcts::CtsOptions::setObserver;
%ignore dcts::CtsOptions::getObserver;

// Enum: CtsOptions::NdrStrategy
%typemap(typecheck) CtsOptions::NdrStrategy {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "NONE") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "ROOT_ONLY") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "HALF") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "FULL") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) CtsOptions::NdrStrategy {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "ROOT_ONLY") == 0) {
    $1 = CtsOptions::NdrStrategy::ROOT_ONLY;
  } else if (strcasecmp(str, "HALF") == 0) {
    $1 = CtsOptions::NdrStrategy::HALF;
  } else if (strcasecmp(str, "FULL") == 0) {
    $1 = CtsOptions::NdrStrategy::FULL;
  } else {
    $1 = CtsOptions::NdrStrategy::NONE;
  };
}

%inline %{

namespace dcts {

void
set_sink_clustering(bool enable)
{
  getTriton3DCts()->getParms()->setSinkClustering(enable);
}

void
set_plot_option(bool plot)
{
  getTriton3DCts()->getParms()->setPlotSolution(plot);
}

void
set_sink_clustering_levels(unsigned levels)
{
  getTriton3DCts()->getParms()->setSinkClusteringLevels(levels);
}

void
set_max_char_cap(double cap)
{
  getTriton3DCts()->getParms()->setMaxCharCap(cap);
}

void
set_max_char_slew(double slew)
{
  getTriton3DCts()->getParms()->setMaxCharSlew(slew);
}

void
set_wire_segment_distance_unit(unsigned unit)
{
  getTriton3DCts()->getParms()->setWireSegmentUnit(unit);
}

void
set_root_buffer(const char* buffer)
{
  getTriton3DCts()->setRootBuffer(buffer);
}

void
set_slew_steps(int steps)
{
  getTriton3DCts()->getParms()->setSlewSteps(steps);
}

void
set_cap_steps(int steps)
{
  getTriton3DCts()->getParms()->setCapSteps(steps);
}

void
set_metric_output(const char* file)
{
  getTriton3DCts()->getParms()->setMetricsFile(file);
}

void
set_debug_cmd()
{
  getTriton3DCts()->getParms()->setObserver(std::make_unique<CtsGraphics>());
}

void
report_cts_metrics()
{
  getTriton3DCts()->reportCtsMetrics();
}

void
set_tree_buf(const char* buffer)
{
  getTriton3DCts()->getParms()->setTreeBuffer(buffer);
}

void
set_distance_between_buffers(int distance)
{
  getTriton3DCts()->getParms()->setSimpleSegmentsEnabled(true);
  getTriton3DCts()->getParms()->setBufferDistance(distance);
}

void
set_branching_point_buffers_distance(int distance)
{
  getTriton3DCts()->getParms()->setVertexBuffersEnabled(true);
  getTriton3DCts()->getParms()->setVertexBufferDistance(distance);
}

void
set_clustering_exponent(unsigned power)
{
  getTriton3DCts()->getParms()->setClusteringPower(power);
}

void
set_clustering_unbalance_ratio(double ratio)
{
  getTriton3DCts()->getParms()->setClusteringCapacity(ratio);
}

void
set_sink_clustering_size(unsigned size)
{
  getTriton3DCts()->getParms()->setSinkClusteringSize(size);
}

void
set_clustering_diameter(double distance)
{
  getTriton3DCts()->getParms()->setMaxDiameter(distance);
}

void
set_macro_clustering_size(unsigned size)
{
  getTriton3DCts()->getParms()->setMacroClusteringSize(size);
}

void
set_macro_clustering_diameter(double distance)
{
  getTriton3DCts()->getParms()->setMacroMaxDiameter(distance);
}

void
set_num_static_layers(unsigned num)
{
  getTriton3DCts()->getParms()->setNumStaticLayers(num);
}

void
set_sink_buffer(const char* buffer)
{
  getTriton3DCts()->setSinkBuffer(buffer);
}

void
report_characterization()
{
  getTriton3DCts()->getCharacterization()->report();
}

void
report_wire_segments(unsigned length, unsigned load, unsigned outputSlew)
{
  getTriton3DCts()->getCharacterization()->reportSegments(length, load, outputSlew);
}

int
set_clock_nets(const char* names)
{
  return getTriton3DCts()->setClockNets(names);
}

void
set_skip_clock_nets(odb::dbNet* net)
{
  getTriton3DCts()->getParms()->setSkipNets(net);
}

void
set_buffer_list(const char* buffers)
{
  getTriton3DCts()->setBufferList(buffers);
}

void
set_obstruction_aware(bool obs)
{
  getTriton3DCts()->getParms()->setObstructionAware(obs);
}

void
set_apply_ndr(CtsOptions::NdrStrategy strategy)
{
  getTriton3DCts()->getParms()->setApplyNDR(strategy);
}

void
set_insertion_delay(bool insDelay)
{
  getTriton3DCts()->getParms()->enableInsertionDelay(insDelay);
}

void
set_sink_buffer_max_cap_derate(double derate)
{
  getTriton3DCts()->getParms()->setSinkBufferMaxCapDerate(derate);
}

void
set_dummy_load(bool dummyLoad)
{
  getTriton3DCts()->getParms()->enableDummyLoad(dummyLoad);
}

void
set_delay_buffer_derate(float derate)
{
  getTriton3DCts()->getParms()->setDelayBufferDerate(derate);
}

void
set_cts_library(const char* name)
{
  getTriton3DCts()->getParms()->setCtsLibrary(name);
}

void
set_repair_clock_nets(bool value)
{
  getTriton3DCts()->getParms()->setRepairClockNets(value);
}
 
void
run_triton_cts()
{
  getTriton3DCts()->runTritonCts();
}

const char*
get_ndr_strategy()
{
  return getTriton3DCts()->getParms()->getApplyNdrName();
}

int
get_branching_buffers_distance()
{
  return getTriton3DCts()->getParms()->getVertexBufferDistance();
}
std::string
get_buffer_list()
{
  return getTriton3DCts()->getParms()->getBufferListToString();
}
unsigned
get_clustering_exponent()
{
  return getTriton3DCts()->getParms()->getClusteringPower();
}

double
get_clustering_unbalance_ratio()
{
  return getTriton3DCts()->getParms()->getClusteringCapacity();
}

float
get_delay_buffer_derate()
{
  return getTriton3DCts()->getParms()->getDelayBufferDerate();
}

int
get_distance_between_buffers()
{
  return getTriton3DCts()->getParms()->getBufferDistance();
}

const char*
get_library()
{
  return getTriton3DCts()->getParms()->getCtsLibrary();
}

double
get_macro_clustering_max_diameter()
{
  return getTriton3DCts()->getParms()->getMacroMaxDiameter();
}

unsigned
get_macro_clustering_size()
{
  return getTriton3DCts()->getParms()->getMacroSinkClusteringSize();
}

unsigned
get_num_static_layers()
{
  return getTriton3DCts()->getParms()->getNumStaticLayers();
}

std::string
get_root_buffer()
{
  return getTriton3DCts()->getRootBufferToString();
}

double
get_sink_buffer_max_cap_derate()
{
  return getTriton3DCts()->getParms()->getSinkBufferMaxCapDerate();
}

unsigned
get_sink_clustering_levels()
{
  return getTriton3DCts()->getParms()->getSinkClusteringLevels();
}

double
get_sink_clustering_max_diameter()
{
  return getTriton3DCts()->getParms()->getMaxDiameter();
}

unsigned
get_sink_clustering_size()
{
  return getTriton3DCts()->getParms()->getSinkClusteringSize();
}

std::string
get_skip_nets()
{
  return getTriton3DCts()->getParms()->getSkipNetsToString();
}
std::string
get_tree_buf()
{
  return getTriton3DCts()->getParms()->getTreeBuffer();
}

unsigned
get_wire_unit()
{
  return getTriton3DCts()->getParms()->getWireSegmentUnit();
}

void
set_db_unit(bool reset)
{
  if (reset) {
    getTriton3DCts()->getParms()->setDbUnits(-1);
  } else {
    odb::dbBlock* block = getTriton3DCts()->getBlock();
    getTriton3DCts()->getParms()->setDbUnits(block->getDbUnitsPerMicron());
  }
}
void
reset_apply_ndr()
{
  getTriton3DCts()->getParms()->resetApplyNDR();
}
void
reset_buffer_list()
{
  getTriton3DCts()->getParms()->resetBufferList();
}
void
reset_branching_point_buffers_distance()
{
  getTriton3DCts()->getParms()->resetVertexBufferDistance();
}
void
reset_clustering_exponent()
{
  getTriton3DCts()->getParms()->resetClusteringPower();
}
void
reset_clustering_unbalance_ratio()
{
  getTriton3DCts()->getParms()->resetClusteringCapacity();
}
void
reset_delay_buffer_derate()
{
  getTriton3DCts()->getParms()->resetDelayBufferDerate();
}
void
reset_distance_between_buffers()
{
  getTriton3DCts()->getParms()->resetBufferDistance();
}
void
reset_cts_library()
{
  getTriton3DCts()->getParms()->resetCtsLibrary();
}
void
reset_macro_clustering_diameter()
{
  getTriton3DCts()->getParms()->resetMacroMaxDiameter();
}
void
reset_macro_clustering_size()
{
  getTriton3DCts()->getParms()->resetMacroClusteringSize();
}
void
reset_num_static_layers()
{
  getTriton3DCts()->getParms()->resetNumStaticLayers();
}
void
reset_root_buffer()
{
  getTriton3DCts()->resetRootBuffer();
}
void
reset_sink_buffer_max_cap_derate()
{
  getTriton3DCts()->getParms()->resetSinkBufferMaxCapDerate();
}
void
reset_sink_clustering_levels()
{
  getTriton3DCts()->getParms()->resetSinkClusteringLevels();
}
void
reset_clustering_diameter()
{
  getTriton3DCts()->getParms()->resetMaxDiameter();
}
void
reset_sink_clustering_size()
{
  getTriton3DCts()->getParms()->resetSinkClusteringSize();
}
void
reset_skip_nets()
{
  getTriton3DCts()->getParms()->resetSkipNets();
}
void
reset_tree_buf()
{
  getTriton3DCts()->getParms()->resetTreeBuffer();
}
void
reset_wire_segment_distance_unit()
{
  getTriton3DCts()->getParms()->resetWireSegmentUnit();
}

void
extract_ff_timing_graph_verilog(const char* verilog_file, const char* output_file)
{
  getTriton3DCts()->extractFFGraphFromVerilog(verilog_file, output_file);
}

// ODB+STA-based FF timing graph extraction.
// No Verilog file needed — registers discovered from ODB Liberty cells,
// edges extracted via STA findPathEnds. Handles macro-included designs.
void
extract_ff_timing_graph_odb(const char* output_file)
{
  getTriton3DCts()->extractFFGraphFromODB(output_file);
}

// Exhaustive IO timing edge extraction via STA per-port.
// Replaces Tcl extract_io_timing_edges_phase1a (worst-N truncation).
void
extract_io_timing_edges(const char* verilog_file, const char* output_file)
{
  getTriton3DCts()->extractIOTimingEdges(verilog_file, output_file);
}

// ODB-based IO timing edge extraction.
// No Verilog file needed — uses ODB register collection + STA per-port.
// Handles macro-included designs where Verilog parsing fails.
void
extract_io_timing_edges_odb(const char* output_file)
{
  getTriton3DCts()->extractIOTimingEdgesFromODB(output_file);
}

// V98_CTRCPP: Phase 4-CTR refinement (re-tap + downsize, batched
// measure-and-revert) — replaces the ctr2_plan.py + ctr_stage.tcl pair.
void
run_ctr_refinement()
{
  getTriton3DCts()->runCtrRefinement();
}

// Pre-CTS skew targets for SG-CTS
void
load_skew_targets(const char* csv_file)
{
  getTriton3DCts()->loadSkewTargets(csv_file);
}

// Estimate per-FF physical achievability bounds.
// Computes HB via delay (t_via = 0.693 * R_HB * C_HB) from tech and writes
// bounds CSV for the LP-SAFETY pre-CTS skew solver (pre_cts_skew_lp.py).
// output_csv:  destination path (ff_name,tier,t_min_ns,t_max_ns)
// max_skew_ns: same-tier CTS delay budget in ns (default 0.1 = 100ps)
void
estimate_leaf_latencies(const char* output_csv, float max_skew_ns)
{
  getTriton3DCts()->estimateLeafLatencies(output_csv, max_skew_ns);
}

// C++ LP solver for skew targeting.
// Extracts timing graph from STA in-memory, solves LP-TNS with OR-Tools GLOP,
// stores results in Cts3DDatabase skewTargetMap_ (no CSV round-trip).
void
solve_skew_lp(const char* verilog_file,
              float sigma_local, float sigma_pi,
              float lambda_reg, float hold_margin,
              float gamma_wns, float weight_io,
              int hard_pi_hold, float max_skew)
{
  getTriton3DCts()->solveSkewLp(verilog_file,
                              sigma_local, sigma_pi,
                              lambda_reg, hold_margin,
                              gamma_wns, weight_io,
                              hard_pi_hold != 0, max_skew);
}

// Old solve_buffer_sizing_lp(float, float) removed — replaced by V53_FM_BUF version below.

// Update cascaded TAP chain connections in ODB.
// Re-reads LP targets CSV and reconnects FFs to the correct TAP depth level
// without rebuilding the tree.
void
update_tap_depths(const char* targets_csv)
{
  getTriton3DCts()->updateTapDepths(targets_csv);
}

// C++ LP-based buffer sizing with Liberty delays.
// Replaces Python buffer_sizing_lp.py. Outputs sizing decisions CSV.
void
solve_buffer_sizing_lp(const char* timing_csv,
                       const char* output_csv,
                       const char* skew_targets_csv,
                       double hold_weight,
                       double reg_weight,
                       double skew_weight,
                       double setup_margin_ps,
                       double hold_margin_ps)
{
  // Set timing CSV path as env var for BufSizingLpSolver to read
  setenv("BUF_SIZING_TIMING_CSV", timing_csv, 1);
  getTriton3DCts()->solveBufferSizingLp(output_csv, skew_targets_csv,
                                       hold_weight, reg_weight, skew_weight,
                                       setup_margin_ps, hold_margin_ps);
}

}  // namespace dcts

%} //inline
