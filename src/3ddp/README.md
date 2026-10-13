# True-3D Detailed Placement

The `3ddp` module provides timing-driven detailed placement for a two-tier
3D design.  It imports the combined OpenDB design, allocates candidate hybrid
bond terminals (HBTs), builds local timing models, moves cells in three
dimensions, and legalizes the two tiers before writing the updated database.

The flow-facing Tcl command is:

```tcl
3d_detailed_placement \
  -tech asap7_3D \
  -hbt /path/to/hbt_array.json \
  -max_displacement {5 1} \
  -m 0.2 \
  -report_file_name placement.rpt
```

`-tech` selects the upper/bottom technology pairing.  The HBT JSON describes
the candidate inter-tier connection sites for the current run.  The TaiWei
adapter generates this file from the active repaired DEF, technology LEF, and
design configuration; it is not a public downstream handoff.

## Architecture

- `include/tdp/` contains the public TDP engine and initialization interfaces.
- `src/infrastructure/` represents the two-tier placement database.
- `src/hbta/` implements HBT modeling and allocation.
- `src/eva/` builds Steiner routes and evaluates timing and parasitics.
- `src/placement/` implements the 3D placement search and legalization.
- `src/lgista/` is the private local timing-analysis engine used by TDP.
- `src/Timingdp.tcl` and `src/Timingdp.i` define the Tcl command layer.

TDP reuses OpenROAD database, STA, resizer, global-router, and detailed
placement services.  Its engine instance is owned inside this module, so it
does not add state to `ord::OpenRoad` and does not replace the existing DPL or
GRT commands.  The TaiWei stage entry is `ord-detail-place-3d`; see
`TaiWei-Pin-3D/scripts_3ddp/README.md` for the complete handoff contract.
