# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The TDP Authors
sta::define_cmd_args "timing_detailed_placement" { \
                           [-max_displacement disp|{disp_x disp_y}] \
                           [-m margin] \
                           [-report_file_name file_name] \
                           [-tech tech] \
                           [-hbt hbt_file]}

# 解析 timing_detailed_placement / report_timing_stats 的公共参数，
# 完成 -tech/-hbt 校验、displacement 换算、design_set 与建 row。
# 解析结果通过 upvar 写回调用方的局部变量：
#   max_displacement_x max_displacement_y margin file_name tech hbt
proc parse_tdp_args { cmd_name args_name } {
  upvar 1 $args_name args
  upvar 1 max_displacement_x max_displacement_x
  upvar 1 max_displacement_y max_displacement_y
  upvar 1 margin margin
  upvar 1 file_name file_name
  upvar 1 tech tech
  upvar 1 hbt hbt

  sta::parse_key_args $cmd_name args \
    keys {-max_displacement -m -report_file_name -tech -hbt}

  if { [info exists keys(-max_displacement)] } {
    set max_displacement $keys(-max_displacement)
    if { [llength $max_displacement] == 1 } {
      sta::check_positive_integer "-max_displacement" $max_displacement
      set max_displacement_x $max_displacement
      set max_displacement_y $max_displacement
    } elseif { [llength $max_displacement] == 2 } {
      lassign $max_displacement max_displacement_x max_displacement_y
      sta::check_positive_integer "-max_displacement" $max_displacement_x
      sta::check_positive_integer "-max_displacement" $max_displacement_y
    } else {
      sta::error TDP 31 "-max_displacement disp|{disp_x disp_y}"
    }
  } else {
    # use default displacement
    set max_displacement_x 0
    set max_displacement_y 0
  }
  if { [info exists keys(-m)] } {
    set margin $keys(-m)
    sta::check_positive_float "-m" $margin
  } else {
      set margin 2
  }
  set file_name ""
  if { [info exists keys(-report_file_name)] } {
    set file_name $keys(-report_file_name)
  }

  if { ![info exists keys(-tech)] } {
    utl::error "TDP" 40 "-tech is required \
(asap7_3D|nangate45_3D|asap7_nangate45_3D|nangate45_asap7_3D)."
  }
  set tech $keys(-tech)
  if { [lsearch -exact \
        {asap7_3D nangate45_3D asap7_nangate45_3D nangate45_asap7_3D} $tech] < 0 } {
    utl::error "TDP" 41 "invalid -tech '$tech'. Must be one of \
asap7_3D, nangate45_3D, asap7_nangate45_3D, nangate45_asap7_3D."
  }
  if { ![info exists keys(-hbt)] } {
    utl::error "TDP" 42 "-hbt hbt_array_json is required."
  }
  set hbt $keys(-hbt)

  sta::check_argc_eq0 $cmd_name $args
  if { ![ord::db_has_core_rows] } {
    utl::error "TDP" 27 "no rows defined in design. Use initialize_floorplan to add rows."
  }
  set site [tdp::get_row_site]
  # Convert displacement from microns to sites.
  set max_displacement_x [expr [ord::microns_to_dbu $max_displacement_x] \
    / [$site getWidth]]
  set max_displacement_y [expr [ord::microns_to_dbu $max_displacement_y] \
    / [$site getHeight]]
  tdp::design_set
  newrow $margin $tech
}

proc timing_detailed_placement { args } {
  parse_tdp_args "timing_detailed_placement" args
  tdp::timing_detailed_placement_cmd $max_displacement_x $max_displacement_y \
    $margin $file_name $tech $hbt
  tdp::report_legalization_stats
}

# Stable flow-facing name used by TaiWei-Pin-3D.  Keep the existing
# timing_detailed_placement command for compatibility with TDP tests.
sta::define_cmd_args "3d_detailed_placement" { \
                           [-max_displacement disp|{disp_x disp_y}] \
                           [-m margin] \
                           [-report_file_name file_name] \
                           [-tech tech] \
                           [-hbt hbt_file]}

proc 3d_detailed_placement { args } {
  timing_detailed_placement {*}$args
}

sta::define_cmd_args "report_timing_stats" { \
                           [-max_displacement disp|{disp_x disp_y}] \
                           [-m margin] \
                           [-report_file_name file_name] \
                           [-tech tech] \
                           [-hbt hbt_file]}

proc report_timing_stats { args } {
  parse_tdp_args "report_timing_stats" args
  tdp::report_timing_stats_cmd $max_displacement_x $max_displacement_y \
    $margin $tech $hbt
}

proc newrow {margin tech} {

lassign [ord::get_die_area] die_lx die_ly die_ux die_uy

set core_lx [expr {$die_lx + $margin}]
set core_ly [expr {$die_ly + $margin}]
set core_ux [expr {$die_ux - $margin}]
set core_uy [expr {$die_uy - $margin}]

# 上层：先按 -tech 确定的上层 site 建 row，立刻抓进 tdp 数据库并标为 upper。
make_rows -core_area [list $core_lx $core_ly $core_ux $core_uy] \
  -site [tdp::get_tier_site "upper" $tech]
tdp::create_new_row "upper"
# 下层：第二次 make_rows 会清掉上层的 OpenDB row（无妨，已抓进 tdp 库）。
make_rows -core_area [list $core_lx $core_ly $core_ux $core_uy] \
  -site [tdp::get_tier_site "bottom" $tech]
tdp::create_new_row "bottom"

}
# sta::define_cmd_args "check_placement" {[-verbose] \
#                                         [-disallow_one_site_gaps] \
#                                         [-report_file_name file_name]}

# proc check_placement { args } {
#   if { [ord::get_db_block] == "NULL" } {
#     utl::error TDP 103 "No design block found."
#   }

#   sta::parse_key_args "check_placement" args \
#     keys {-report_file_name} flags {-verbose -disallow_one_site_gaps}
#   set verbose [info exists flags(-verbose)]
#   sta::check_argc_eq0 "check_placement" $args
#   set file_name ""
#   if { [info exists keys(-report_file_name)] } {
#     set file_name $keys(-report_file_name)
#   }
#   if { [info exists flags(-disallow_one_site_gaps)] } {
#     utl::warn TDP 4 "-disallow_one_site_gaps is deprecated"
#   }
#   tdp::check_placement_cmd $verbose $file_name
# }

namespace eval tdp {
# min_displacement is the smallest displacement to draw
# measured as a multiple of row_height.
# proc detailed_placement_debug { args } {
#   sta::parse_key_args "detailed_placement_debug" args \
#     keys {-instance -min_displacement -jump_moves} flags {-iterative}


#   if { [info exists keys(-min_displacement)] } {
#     set min_displacement $keys(-min_displacement)
#   } else {
#     set min_displacement 0
#   }

#   set jump_moves 0
#   if { [info exists keys(-jump_moves)] } {
#     set jump_moves $keys(-jump_moves)
#     sta::check_positive_integer "-jump_moves" $jump_moves
#   }

#   if { [info exists keys(-instance)] } {
#     set instance_name $keys(-instance)
#     set block [ord::get_db_block]
#     set debug_instance [$block findInst $instance_name]
#     if { $debug_instance == "NULL" } {
#       utl::error TDP 32 "Debug instance $instance_name not found."
#     }
#   } else {
#     set debug_instance "NULL"
#   }

#   tdp::set_debug_cmd $min_displacement $debug_instance $jump_moves [info exists flags(-iterative)]
# }

proc get_masters_arg { arg_name arg } {
  set masters {}
  # Expand master name regexps
  set db [ord::get_db]
  foreach name $arg {
    set matched 0
    foreach lib [$db getLibs] {
      foreach master [$lib getMasters] {
        set master_name [$master getConstName]
        if { [string match $name $master_name] } {
          lappend masters $master
          set matched 1
        }
      }
    }
    if { !$matched } {
      utl::warn "TDP" 28 "$name did not match any masters."
    }
  }
  if { [llength $arg] > 0 && [llength $masters] == 0 } {
    utl::error "TDP" 39 "\"$arg\" did not match any masters.
This could be due to a change from using regex to glob to search for cell masters. https://github.com/The-OpenROAD-Project/OpenROAD/pull/3210"
  }
  return $masters
}

proc get_inst_bbox { inst_name } {
  set block [ord::get_db_block]
  set inst [$block findInst $inst_name]
  if { $inst != "NULL" } {
    set bbox [$inst getBBox]
    return "[$bbox xMin] [$bbox yMin] [$bbox xMax] [$bbox yMax]"
  } else {
    utl::error "TDP" 29 "cannot find instance $inst_name"
  }
}

proc get_inst_grid_bbox { inst_name } {
  set block [ord::get_db_block]
  set inst [$block findInst $inst_name]
  set rows [$block getRows]
  set site [[lindex $rows 0] getSite]
  set width [$site getWidth]
  set height [$site getHeight]
  if { $inst != "NULL" } {
    set bbox [$inst getBBox]
    return "[format_grid [$bbox xMin] $width] [format_grid [$bbox yMin] $height]\
            [format_grid [$bbox xMax] $width] [format_grid [$bbox yMax] $height]"
  } else {
    utl::error "TDP" 30 "cannot find instance $inst_name"
  }
}

proc format_grid { x w } {
  if { $x % $w == 0 } {
    return [expr $x / $w]
  } else {
    return [format "%.2f" [expr $x / double($w)]]
  }
}

proc get_row_site { } {
  return [[lindex [[ord::get_db_block] getRows] 0] getSite]
}
}
