#-----------------------------------------------------------------------------
# Vivado 2025.2 - build the IK investigation platform.
#
#   Zynq-7000 PS  ->  AXI interconnect  ->  HLS kernels on AXI4-Lite
#
# Run:
#   vivado -mode batch -source scripts/build_vivado.tcl
#   vivado -mode batch -source scripts/build_vivado.tcl -tclargs --no-bit
#   vivado -mode batch -source scripts/build_vivado.tcl -tclargs \
#       --kernels mat_mul_kernel,mat_inv_kernel,ik_analytic_kernel,ik_dls_kernel
#
# Produces build/vivado/ik_platform.xsa for the Vitis application build.
#
# Everything version-sensitive is looked up rather than hard-coded: IP VLNVs
# come from get_ipdefs, and the HLS output directories are found by searching
# for component.xml.  The layout of a Vitis HLS work_dir has changed more than
# once across releases and hard-coding it is the usual reason these scripts
# rot.
#-----------------------------------------------------------------------------

set script_dir [file normalize [file dirname [info script]]]
set root_dir   [file normalize $script_dir/..]

set part      "xc7z020clg400-1"
set proj_name "ik_platform"
set bd_name   "ik_bd"
set build_dir $root_dir/build/vivado
set hls_dir   $root_dir/hls/build
set run_bit   1
set jobs      24
set allow_tns 0

# PL clock, MHz.  HLS schedules the kernels against clock=10 with 12.5%
# uncertainty, but that uncertainty is a guess at routing it has not done yet:
# at 100 MHz this design closed placement and then missed setup by 1.582 ns
# post-route, needing 11.582 ns.  80 MHz (12.5 ns) clears that with ~0.9 ns in
# hand, which is enough margin to survive the paths moving when the clock
# changes.  Raise it with --clk once a build shows positive WNS to spare.
#
# Nothing in the measurement depends on this number: ik_driver.c times with the
# PS global timer and main.c reports nanoseconds, so the PL figure stays
# correct - it just gets proportionally larger.  Report the clock alongside any
# latency number, because a PS-vs-PL ratio is meaningless without it.
set clk_mhz   80

# Which kernels go in the bitstream.
#
# The xc7z020 has 220 DSP48E1 and 53,200 LUT, and both are now binding.
# ik_dls_kernel used to need ~126% of the DSPs and could not be built at all;
# raising II at the five multiply sites listed in hls/include/ik_config.hpp
# brought it to 186 DSP (85%) and 51k LUT (96%).  It fits - on its own.
#
# It does NOT fit alongside ik_analytic_kernel, which is 75% DSP and 86% LUT
# by itself, so the single-bitstream head-to-head is still out of reach on
# this part.  Build them separately and compare at the same clock.  Note the
# LUT figure: at 96% occupancy expect placement to be slow and timing to be
# tight, and be ready to drop --clk below 80 if the gate below trips.
#
# Pass --kernels to build any other combination, e.g.
#   --kernels ik_analytic_kernel
#   --kernels mat_mul_kernel,mat_inv_kernel
# to characterise the matrix primitives in their own bitstream.
#
# sw/src/main.c needs no matching edit for any of these: it resolves each
# kernel's base address from xparameters.h, falls back to 0 when an IP is
# absent, and tests whatever is actually present.
set kernels {ik_dls_kernel}

# ---------------------------------------------------------------- args ----
for {set i 0} {$i < $argc} {incr i} {
    switch -- [lindex $argv $i] {
        "--no-bit"  { set run_bit 0 }
        "--part"    { incr i; set part [lindex $argv $i] }
        "--hls-dir" { incr i; set hls_dir [file normalize [lindex $argv $i]] }
        "--jobs"    { incr i; set jobs [lindex $argv $i] }
        "--kernels" { incr i; set kernels [split [lindex $argv $i] ","] }
        "--clk"     { incr i; set clk_mhz [lindex $argv $i] }
        "--allow-timing-fail" { set allow_tns 1 }
        default     { puts "WARNING: ignoring unknown argument [lindex $argv $i]" }
    }
}

# ------------------------------------------------------- helper procs ----

# Latest VLNV for an IP name, so the script survives a core revision bump.
proc latest_vlnv {name} {
    set defs [get_ipdefs -all "xilinx.com:ip:${name}:*"]
    if {[llength $defs] == 0} {
        error "IP '$name' not found in the catalog for this Vivado install."
    }
    return [lindex [lsort -dictionary $defs] end]
}

# Every directory that directly contains a component.xml is a packaged IP;
# return their parents, which is what ip_repo_paths wants.  Walking beats
# globbing a fixed depth because the nesting inside a Vitis HLS work_dir has
# changed between releases.
proc find_ip_repos {root} {
    set repos {}
    if {![file isdirectory $root]} { return $repos }
    set stack [list $root]
    while {[llength $stack]} {
        set d     [lindex $stack end]
        set stack [lrange $stack 0 end-1]
        if {[file exists [file join $d component.xml]]} {
            lappend repos [file dirname $d]
            continue
        }
        foreach sub [glob -nocomplain -directory $d -types d *] {
            lappend stack $sub
        }
    }
    return [lsort -unique $repos]
}

# ------------------------------------------------------------ project ----
file mkdir $build_dir
create_project $proj_name $build_dir -part $part -force
set_property target_language Verilog [current_project]

# -jobs on launch_runs (below) parallelises across runs/strategies; it does
# not by itself set how many threads place_design/route_design use inside a
# single run.  That is general.maxThreads, capped at 32 by Vivado regardless
# of what is asked for.  Set both from the same $jobs so --jobs actually
# uses the machine.
set_param general.maxThreads $jobs

# PYNQ-Z1/Z2 board files are optional; use them when installed so the PS gets
# the right DDR and clock configuration, otherwise fall back to the raw part.
set board_candidates {
    "tul.com.tw:pynq-z2:part0:*"
    "www.digilentinc.com:pynq-z1:part0:*"
    "digilentinc.com:pynq-z1:part0:*"
}
set board_set 0
foreach bc $board_candidates {
    set hits [get_board_parts -quiet $bc]
    if {[llength $hits]} {
        set_property board_part [lindex [lsort -dictionary $hits] end] [current_project]
        puts "INFO: using board_part [get_property board_part [current_project]]"
        set board_set 1
        break
    }
}
if {!$board_set} {
    puts "INFO: no PYNQ board files found; building against part $part only."
}

# --------------------------------------------------------- HLS IP repo ----
set repos [find_ip_repos $hls_dir]
if {[llength $repos] == 0} {
    error "No packaged HLS IP found under '$hls_dir'.\
           Run 'make -C hls ip' first."
}
puts "INFO: HLS IP repositories:"
foreach r $repos { puts "        $r" }
set_property ip_repo_paths $repos [current_project]
update_ip_catalog -rebuild

# Resolve each kernel's VLNV from the catalog.
puts "INFO: building kernels: $kernels"
set kernel_vlnv {}
foreach k $kernels {
    set hits [get_ipdefs -all "xilinx.com:hls:${k}:*"]
    if {[llength $hits] == 0} {
        error "Packaged IP for '$k' not in the catalog. Run 'make -C hls ip-$k'."
    }
    set v [lindex [lsort -dictionary $hits] end]
    dict set kernel_vlnv $k $v
    puts "INFO: $k -> $v"
}

# ----------------------------------------------------- block design ----
create_bd_design $bd_name

# Zynq PS.
set ps [create_bd_cell -type ip -vlnv [latest_vlnv processing_system7] ps7]
if {$board_set} {
    apply_bd_automation -rule xilinx.com:bd_rule:processing_system7 \
        -config {make_external "FIXED_IO, DDR" apply_board_preset "1" \
                 Master "Disable" Slave "Disable"} $ps
} else {
    apply_bd_automation -rule xilinx.com:bd_rule:processing_system7 \
        -config {make_external "FIXED_IO, DDR" apply_board_preset "0" \
                 Master "Disable" Slave "Disable"} $ps
}

# One general-purpose AXI master; PL clock from $clk_mhz (see the note there).
puts "INFO: PL clock $clk_mhz MHz"
set_property -dict [list \
    CONFIG.PCW_USE_M_AXI_GP0        {1} \
    CONFIG.PCW_FPGA0_PERIPHERAL_FREQMHZ $clk_mhz \
    CONFIG.PCW_EN_CLK0_PORT         {1} \
] $ps

# Instantiate the kernels.
set cells {}
foreach k $kernels {
    lappend cells [create_bd_cell -type ip -vlnv [dict get $kernel_vlnv $k] $k]
}

# Let automation build the interconnect and the reset network.  Each kernel
# exposes a single AXI4-Lite control bundle named s_axi_CTRL.
foreach k $kernels {
    apply_bd_automation -rule xilinx.com:bd_rule:axi4 \
        -config [list Master "/ps7/M_AXI_GP0" Clk "Auto"] \
        [get_bd_intf_pins $k/s_axi_CTRL]
}

regenerate_bd_layout
assign_bd_address
validate_bd_design
save_bd_design

puts "INFO: address map"
foreach seg [get_bd_addr_segs -of_objects [get_bd_addr_spaces ps7/Data]] {
    puts [format "        %-46s %s  %s" $seg \
              [get_property OFFSET $seg] [get_property RANGE $seg]]
}

# ------------------------------------------------------------ build ----
set bd_file [get_files $bd_name.bd]
make_wrapper -files $bd_file -top
add_files -norecurse [file join $build_dir $proj_name.gen sources_1 bd \
                          $bd_name hdl ${bd_name}_wrapper.v]
set_property top ${bd_name}_wrapper [current_fileset]
update_compile_order -fileset sources_1

# Synthesise first and print the per-kernel resource breakdown before
# implementation.  place_design refuses to run at all if the design is over
# budget, and its DRC message only reports the total for the whole device -
# which tells you that you are over but not by whose doing.  This table is the
# thing you actually need, and it costs nothing to always emit it.
launch_runs synth_1 -jobs $jobs
wait_on_run synth_1
if {[get_property PROGRESS [get_runs synth_1]] != "100%"} {
    error "Synthesis failed; see $build_dir/$proj_name.runs/synth_1"
}
open_run synth_1 -name synth_1
puts "INFO: post-synthesis utilisation by kernel"
report_utilization -hierarchical -hierarchical_depth 2

if {$run_bit} {
    launch_runs impl_1 -to_step write_bitstream -jobs $jobs
    wait_on_run impl_1
    if {[get_property PROGRESS [get_runs impl_1]] != "100%"} {
        error "Implementation failed; see $build_dir/$proj_name.runs/impl_1"
    }

    open_run impl_1
    set wns [get_property SLACK [get_timing_paths -delay_type max]]
    set whs [get_property SLACK [get_timing_paths -delay_type min]]
    puts "INFO: timing summary"
    puts "        WNS = $wns"
    puts "        WHS = $whs"

    # A bitstream that misses timing still programs, still runs, and still
    # returns answers - occasionally wrong ones, dependent on temperature and
    # on which path happened to lose the race.  That failure mode is
    # indistinguishable from a bug in the kernel, and this investigation exists
    # to attribute latency differences to the target rather than to chance, so
    # exporting one silently is the worst thing this script could do.  Print
    # the paths that failed and stop.  --allow-timing-fail overrides it when
    # you want the XSA for a flow check rather than for a measurement.
    if {$wns < 0 || $whs < 0} {
        puts "INFO: worst failing paths"
        report_timing_summary -max_paths 1 -report_unconstrained
        report_timing -delay_type max -max_paths 10 -nworst 1 -sort_by slack

        set fmax [format "%.1f" [expr {1000.0 / (1000.0/$clk_mhz - $wns)}]]
        set msg "Implementation missed timing: WNS = $wns, WHS = $whs at\
                 $clk_mhz MHz.\n\
                 The design closes at roughly $fmax MHz; rebuild with\
                 '-tclargs --clk <mhz>' at or below that.\n\
                 No XSA was written - a bitstream with negative slack produces\
                 intermittently wrong results that read as kernel bugs.\n\
                 Pass --allow-timing-fail to export one anyway."
        if {!$allow_tns} { error $msg }
        puts "WARNING: $msg"
        puts "WARNING: exporting anyway; do not trust measurements from this bitstream."
    }

    set xsa [file join $build_dir ${proj_name}.xsa]
    write_hw_platform -fixed -include_bit -force -file $xsa
    puts "INFO: exported $xsa"
} else {
    set xsa [file join $build_dir ${proj_name}_pre.xsa]
    write_hw_platform -fixed -force -file $xsa
    puts "INFO: exported $xsa (no bitstream)"
}

puts "INFO: done."
