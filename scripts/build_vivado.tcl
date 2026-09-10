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
set jobs      8

# ik_dls_kernel does not fit on the xc7z020 yet even after resource-sharing
# work (still ~126% DSP utilisation standalone) - default to the three
# kernels that do, so a bitstream and real on-target numbers are obtainable
# now.  Pass --kernels to override, e.g. once ik_dls fits or to build it in
# isolation: --kernels ik_dls_kernel
set kernels {mat_mul_kernel mat_inv_kernel ik_analytic_kernel}

# ---------------------------------------------------------------- args ----
for {set i 0} {$i < $argc} {incr i} {
    switch -- [lindex $argv $i] {
        "--no-bit"  { set run_bit 0 }
        "--part"    { incr i; set part [lindex $argv $i] }
        "--hls-dir" { incr i; set hls_dir [file normalize [lindex $argv $i]] }
        "--jobs"    { incr i; set jobs [lindex $argv $i] }
        "--kernels" { incr i; set kernels [split [lindex $argv $i] ","] }
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

# One general-purpose AXI master, PL clock at 100 MHz to match the 10 ns
# constraint the kernels were synthesised against.
set_property -dict [list \
    CONFIG.PCW_USE_M_AXI_GP0        {1} \
    CONFIG.PCW_FPGA0_PERIPHERAL_FREQMHZ {100} \
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

if {$run_bit} {
    launch_runs impl_1 -to_step write_bitstream -jobs $jobs
    wait_on_run impl_1
    if {[get_property PROGRESS [get_runs impl_1]] != "100%"} {
        error "Implementation failed; see $build_dir/$proj_name.runs/impl_1"
    }

    open_run impl_1
    puts "INFO: timing summary"
    puts "        WNS = [get_property SLACK [get_timing_paths -delay_type max]]"
    puts "        WHS = [get_property SLACK [get_timing_paths -delay_type min]]"

    set xsa [file join $build_dir ${proj_name}.xsa]
    write_hw_platform -fixed -include_bit -force -file $xsa
    puts "INFO: exported $xsa"
} else {
    set xsa [file join $build_dir ${proj_name}_pre.xsa]
    write_hw_platform -fixed -force -file $xsa
    puts "INFO: exported $xsa (no bitstream)"
}

puts "INFO: done."
