# Vitis HLS 2025.2 build driver.
#
# HLS Classic was removed in 2025.1: open_project / open_solution no longer
# describe a supported flow, and vitis_hls -f is gone.  Everything here goes
# through the config-file interface:
#
#   v++ -c --mode hls --config <cfg> --work_dir <dir>      C synthesis
#   vitis-run --mode hls --csim    --config <cfg> ...      C simulation
#   vitis-run --mode hls --cosim   --config <cfg> ...      C/RTL co-simulation
#   vitis-run --mode hls --package --config <cfg> ...      package as IP
#
# Targets:
#   make csim      run C simulation for every kernel (fixed-point)
#   make syn       run C synthesis, emits the latency reports
#   make cosim     run co-simulation (slow; verifies the RTL, not just the C)
#   make ip        package all four kernels for the Vivado IP catalog
#   make all       syn + ip
#   make reports   summarise latency and utilisation from the syn reports
#   make host      run the host regression (no Vitis needed)
#   make clean
#
# Override the device with:  make PART=xc7z020clg400-1 all
#
# BOARD=vek385 is a convenience for out-of-context resource/latency analysis
# on a Versal part: it sets PART to the VEK385's device and moves the build
# output to build-vek385/ so it never clobbers a PYNQ build's reports. It
# does NOT change anything else about the kernels or the flow_target - this
# is still `flow_target=vivado` Vitis HLS estimation, not a real Versal
# platform (see scripts/build_vivado.tcl, which is Zynq-7000/PYNQ-only: it
# instantiates processing_system7, which VEK385 does not have). Confirm the
# exact part string against your own install
# (`get_parts -filter {board_part=~*vek385*}` in Vivado, or the board_part
# XML under a Versal board_files install) and override with an explicit
# PART= if it differs.
#
#   make BOARD=vek385 syn reports
#   make BOARD=vek385 PART=<exact-part> syn      # override the guessed part

BOARD   ?=
PART    ?=
ifeq ($(BOARD),vek385)
  PART := $(if $(PART),$(PART),xcve2302-sfvc784-2LP-e-S)
endif

KERNELS := matmul matinv ik_analytic ik_dls
BUILD   := $(if $(BOARD),build-$(BOARD),build)
VXX     := v++
VRUN    := vitis-run
PARTOPT := $(if $(PART),--part $(PART))

.PHONY: all csim syn cosim ip reports host clean vectors $(KERNELS)

all: syn ip

# The vectors are inputs to both csim and the host build.
vectors:
	@python3 ../model/gen_vectors.py

host:
	@$(MAKE) -C tb

$(BUILD):
	@mkdir -p $(BUILD)

define hls_rule
csim-$(1): vectors | $(BUILD)
	$(VRUN) --mode hls --csim $(PARTOPT) --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

syn-$(1): | $(BUILD)
	$(VXX) -c --mode hls $(PARTOPT) --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

cosim-$(1): vectors | $(BUILD)
	$(VRUN) --mode hls --cosim $(PARTOPT) --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

ip-$(1): | $(BUILD)
	$(VRUN) --mode hls --package $(PARTOPT) --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)
endef

$(foreach k,$(KERNELS),$(eval $(call hls_rule,$(k))))

csim:  $(addprefix csim-,$(KERNELS))
syn:   $(addprefix syn-,$(KERNELS))
cosim: $(addprefix cosim-,$(KERNELS))
ip:    $(addprefix ip-,$(KERNELS))

# Pull the headline numbers out of the synthesis reports.  The analytic kernel
# should show one latency; the DLS kernel should show a min/max spread.
#
# HLS emits one *_csynth.rpt per pipelined region in the hierarchy, not just
# for the top function - e.g. mat_mul_kernel_Pipeline_MM_LOAD..._csynth.rpt
# alongside mat_mul_kernel_csynth.rpt itself.  `find | head -1` used to grab
# whichever one the filesystem listed first, which was almost always a
# sub-region report with no top-level utilization table.  Match the actual
# top function name instead.
TOP-matmul      := mat_mul_kernel
TOP-matinv      := mat_inv_kernel
TOP-ik_analytic := ik_analytic_kernel
TOP-ik_dls      := ik_dls_kernel

reports:
	@echo "=================================================================="
	@echo " C-SYNTHESIS SUMMARY"
	@echo "=================================================================="
	@for k in $(KERNELS); do \
	  case $$k in \
	    matmul)      top=$(TOP-matmul) ;; \
	    matinv)      top=$(TOP-matinv) ;; \
	    ik_analytic) top=$(TOP-ik_analytic) ;; \
	    ik_dls)      top=$(TOP-ik_dls) ;; \
	  esac; \
	  r=$$(find $(BUILD)/$$k -name "$${top}_csynth.rpt" 2>/dev/null | head -1); \
	  if [ -n "$$r" ]; then \
	    echo ""; echo "--- $$k ($$r)"; \
	    sed -n '/Latency/,/^$$/p' $$r | head -20; \
	    awk '/^== Utilization Estimates/{p=1} p && /^== / && !/Utilization Estimates/{exit} p{print}' $$r; \
	  else \
	    echo ""; echo "--- $$k: no top-level report yet (run 'make syn-$$k')"; \
	  fi; \
	done

clean:
	rm -rf $(BUILD)
	@$(MAKE) -C tb clean
