# 3D Clock Tree Synthesis

The OpenROAD 3D clock tree synthesis module (`3dcts`) implements useful-skew
clock delivery for face-to-face 3D ICs. It builds clock structures across two
tiers, models hybrid-bonding crossings, schedules per-sink clock arrivals with
linear programming, and realizes those arrivals with shared cascaded tap
chains.

The primary Tcl command is `clock_tree_synthesis_3d`. On-the-fly wire-segment
characterization is available through `configure_3dcts_characterization`.
Clock routing resistance and capacitance should be initialized with
`set_wire_rc` before synthesis.

## Commands

Parameters in square brackets are optional.

### Configure 3DCTS characterization

Configure the slew/capacitance range and sampling resolution used to
characterize clock buffer and wire segments.

```tcl
configure_3dcts_characterization \
    [-max_slew max_slew] \
    [-max_cap max_cap] \
    [-slew_steps slew_steps] \
    [-cap_steps cap_steps]
```

| Switch | Description |
| --- | --- |
| `-max_slew` | Maximum slew, in the current time unit, used during characterization. If omitted, the value is derived from the selected buffer Liberty data. |
| `-max_cap` | Maximum capacitance, in the current capacitance unit, used during characterization. If omitted, the value is derived from the selected buffer Liberty data. |
| `-slew_steps` | Number of characterization points along the slew dimension. |
| `-cap_steps` | Number of characterization points along the capacitance dimension. |

### Run 3D clock tree synthesis

```tcl
clock_tree_synthesis_3d \
    [-wire_unit wire_unit] \
    [-buf_list buffer_list] \
    [-root_buf root_buffer] \
    [-clk_nets clock_nets] \
    [-tree_buf tree_buffer] \
    [-distance_between_buffers distance] \
    [-branching_point_buffers_distance distance] \
    [-clustering_exponent exponent] \
    [-clustering_unbalance_ratio ratio] \
    [-sink_clustering_size size] \
    [-sink_clustering_max_diameter diameter] \
    [-macro_clustering_size size] \
    [-macro_clustering_max_diameter diameter] \
    [-sink_clustering_enable] \
    [-sink_clustering_levels levels] \
    [-num_static_layers layers] \
    [-sink_clustering_buffer buffer] \
    [-no_obstruction_aware] \
    [-apply_ndr strategy] \
    [-sink_buffer_max_cap_derate derate] \
    [-dont_use_dummy_load] \
    [-delay_buffer_derate derate] \
    [-library liberty_library] \
    [-repair_clock_nets] \
    [-no_insertion_delay]
```

| Switch | Description |
| --- | --- |
| `-buf_list` | Tcl list of clock buffer masters available for characterized wire segments. For 3D construction, the list can contain tier-specific masters. |
| `-root_buf` | Master used for the clock root. If omitted, the first suitable buffer from `-buf_list` is selected. |
| `-tree_buf` | Explicit master used for tree buffers. |
| `-wire_unit` | Minimum characterized wire-segment length between buffers. |
| `-clk_nets` | Clock roots to synthesize. If omitted, clock roots are discovered from the timing constraints. |
| `-distance_between_buffers` | Fixed distance, in microns, between buffers in simplified segment construction. |
| `-branching_point_buffers_distance` | Minimum branch length, in microns, for inserting a buffer at a branch point. |
| `-clustering_exponent` | Distance exponent used by CKMeans sink clustering. |
| `-clustering_unbalance_ratio` | Maximum CKMeans cluster-capacity ratio. |
| `-sink_clustering_enable` | Pre-cluster register sinks before H-tree construction. |
| `-sink_clustering_size` | Maximum number of register sinks in a cluster. |
| `-sink_clustering_max_diameter` | Maximum register-cluster diameter, in microns. |
| `-sink_clustering_levels` | Number of sink-clustering levels. |
| `-sink_clustering_buffer` | Buffer master used to drive sink clusters. |
| `-macro_clustering_size` | Maximum number of macro sinks in a cluster. |
| `-macro_clustering_max_diameter` | Maximum macro-cluster diameter, in microns. |
| `-num_static_layers` | Number of static layers at the root of the tree. |
| `-no_obstruction_aware` | Disable obstruction-aware clock-buffer placement. |
| `-apply_ndr` | Clock non-default-rule strategy: `none`, `root_only`, `half`, or `full`. |
| `-dont_use_dummy_load` | Disable dummy loads at clock-tree leaves. |
| `-sink_buffer_max_cap_derate` | Derate the max-capacitance limit used for automatic sink-buffer selection. |
| `-delay_buffer_derate` | Scale the number of latency-balancing delay buffers. |
| `-library` | Liberty library from which clock buffers are selected. |
| `-repair_clock_nets` | Repair long clock wires before latency balancing. |
| `-no_insertion_delay` | Ignore sink insertion delay during tree construction and balancing. |

The 3D delivery policy is controlled by environment variables loaded by the
flow. Important groups include:

- tier construction: `HETEROGENEOUS_3D`, `CTS_ENABLE_UNIFIED_TREE`,
  `CTS_ENABLE_PER_TIER_CTS`;
- shared delivery: `CTS_SHARED_TAP`, `CTS_ENABLE_GROUPED_DELAY`,
  `CTS_BRIDGE_PROTOTYPE`;
- useful-skew bounds: `PRE_CTS_MAX_SKEW`, `PRE_CTS_HOLD_MARGIN`,
  `PRE_CTS_SIGMA_LOCAL`;
- refinement: `CTS_ENABLE_CTR`, `CTS_CTR_BATCH`, and the CTR TNS/WNS guards.

The TaiWei stage defines these settings in `scripts_3dcts/cts_params.env`.

### Report 3DCTS metrics

Report the number of clock roots, inserted buffers, clock subnets, and sinks
after a successful run.

```tcl
report_3dcts [-out_file file]
```

If `-out_file` is omitted, the report is written to standard output.

### Set 3DCTS configuration

Persistent CTS parameters may be configured before synthesis:

```tcl
set_3dcts_config \
    [-apply_ndr strategy] \
    [-buf_list buffer_list] \
    [-branching_point_buffers_distance distance] \
    [-clock_buffer_footprint footprint] \
    [-clock_buffer_string pattern] \
    [-clustering_exponent exponent] \
    [-clustering_unbalance_ratio ratio] \
    [-delay_buffer_derate derate] \
    [-distance_between_buffers distance] \
    [-library liberty_library] \
    [-macro_clustering_max_diameter diameter] \
    [-macro_clustering_size size] \
    [-num_static_layers layers] \
    [-root_buf root_buffer] \
    [-sink_buffer_max_cap_derate derate] \
    [-sink_clustering_levels levels] \
    [-sink_clustering_max_diameter diameter] \
    [-sink_clustering_size size] \
    [-skip_nets clock_nets] \
    [-tree_buf tree_buffer] \
    [-wire_unit wire_unit]
```

`-clock_buffer_footprint` and `-clock_buffer_string` are mutually exclusive.
Use `report_3dcts_config` to print the active values:

```tcl
report_3dcts_config
```

Use `reset_3dcts_config` without arguments to reset all values, or pass the
corresponding option names to reset selected values:

```tcl
reset_3dcts_config [-buf_list] [-root_buf] [-wire_unit]
```

## Example

The following example assumes that the design, Liberty data, tier-specific
LEFs, timing constraints, and wire RC have already been loaded:

```tcl
configure_3dcts_characterization -slew_steps 12 -cap_steps 34

clock_tree_synthesis_3d \
    -root_buf BUFx4_ASAP7_75t_R_bottom \
    -buf_list {BUFx4_ASAP7_75t_R_bottom BUFx4_ASAP7_75t_R_upper} \
    -sink_clustering_enable

report_3dcts
```

For a complete placement-to-CTS invocation, use the TaiWei `ord-cts-3d` stage
described in `TaiWei-Pin-3D/scripts_3dcts/README.md`.

## Developer commands

The flow uses the following SWIG commands in the `dcts::` namespace:

| Command | Purpose |
| --- | --- |
| `dcts::extract_ff_timing_graph_verilog` | Extract register timing edges using Verilog register discovery. |
| `dcts::extract_ff_timing_graph_odb` | Extract register timing edges directly from ODB and OpenSTA. |
| `dcts::extract_io_timing_edges` | Extract PI-to-register and register-to-PO timing edges. |
| `dcts::extract_io_timing_edges_odb` | Extract I/O timing edges without Verilog parsing. |
| `dcts::estimate_leaf_latencies` | Estimate per-sink physical delivery bounds. |
| `dcts::solve_skew_lp` | Solve the useful-skew LP and store per-sink targets. |
| `dcts::load_skew_targets` | Load previously written clock-arrival targets. |
| `dcts::update_tap_depths` | Reconnect sinks to tap depths described by a target file. |
| `dcts::solve_buffer_sizing_lp` | Solve the Liberty-aware clock-buffer sizing problem. |
| `dcts::run_ctr_refinement` | Run timing-guarded tap reassignment and buffer downsizing. |
| `dcts::clock_tree_synthesis_3d_debug` | Attach the GUI observer used to inspect tree construction. |

These commands are low-level building blocks. The stage driver supplies their
required data files, ordering, and timing-state preparation.

## Architecture and source guide

```text
TritonCTS
  |
  +-- timing graph extraction and useful-skew scheduling
  +-- Clock / TreeBuilder / HTreeBuilder
  |     +-- tier-aware H-tree construction
  |     +-- shared tap and bridge delivery
  +-- ClockLatencyEstimator
  +-- BufSizingLpSolver
  +-- post-CTS tap refinement
```

| Source | Responsibility |
| --- | --- |
| `include/dcts/TritonCTS.h` | Public engine interface used by the Tcl/SWIG layer. |
| `src/TritonCTS.cpp` | Top-level synthesis engine, ODB updates, timing extraction, LP orchestration, sizing, and refinement. |
| `src/HTreeBuilder.cpp` | 3D H-tree topology, sink grouping, tier-aware construction, shared tap chains, and HBT bridge delivery. |
| `src/Clock.cpp` | In-memory clock instances, subnets, and tree objects. |
| `src/Cts3DDatabase.cpp` | Tier metadata, HBT parasitics, per-sink targets, and timing data shared by the 3D algorithms. |
| `src/CtsSkewLpSolver.cpp` | OR-Tools GLOP formulation for bounded useful-skew scheduling. |
| `src/ClockLatencyEstimator.cpp` | Physical arrival bounds used by the scheduler. |
| `src/BufSizingLpSolver.cpp` | Liberty-aware clock-buffer sizing optimization. |
| `src/VerilogFFExtractor.cpp` | Verilog, ODB, and OpenSTA register/timing-edge extraction. |
| `src/TarpClustering.cpp` | Optional timing-affinity sink clustering. |
| `src/LatencyBalancer.cpp` | Register/macro latency balancing. |
| `src/TechChar.cpp` | On-the-fly characterization of buffer and wire segments. |
| `src/CtsOptions.cpp` | CTS configuration and environment-controlled 3D policies. |
| `src/TritonCTS.tcl` | User-facing Tcl commands and argument checking. |
| `src/TritonCTS.i` | SWIG bindings for the `dcts::` developer commands. |
| `src/Make3DCTS.cpp` | Module initialization and engine construction. |

Useful method entry points include:

| Functionality | Function |
| --- | --- |
| 3D sink preparation | `HTreeBuilder::preSinkClustering()` |
| Tree and shared-tap delivery | `HTreeBuilder::createClockSubNets()` |
| Useful-skew optimization | `CtsSkewLpSolver::solveTns()` |
| Delivery-bound estimation | `ClockLatencyEstimator::estimateAndWrite()` |
| ODB timing graph extraction | `TritonCTS::extractFFGraphFromODB()` |
| Buffer-sizing optimization | `BufSizingLpSolver::solve()` |
| Clock tap refinement | `TritonCTS::runCtrRefinement()` |

## Dependencies

In addition to the standard OpenROAD CTS dependencies, the module uses
OR-Tools GLOP for useful-skew and buffer-sizing linear programs. The CMake
target links `ortools::ortools`.
