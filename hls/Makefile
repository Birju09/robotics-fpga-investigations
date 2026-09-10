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

KERNELS := matmul matinv ik_analytic ik_dls
BUILD   := build
VXX     := v++
VRUN    := vitis-run

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
	$(VRUN) --mode hls --csim --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

syn-$(1): | $(BUILD)
	$(VXX) -c --mode hls --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

cosim-$(1): vectors | $(BUILD)
	$(VRUN) --mode hls --cosim --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)

ip-$(1): | $(BUILD)
	$(VRUN) --mode hls --package --config cfg/$(1).cfg --work_dir $(BUILD)/$(1)
endef

$(foreach k,$(KERNELS),$(eval $(call hls_rule,$(k))))

csim:  $(addprefix csim-,$(KERNELS))
syn:   $(addprefix syn-,$(KERNELS))
cosim: $(addprefix cosim-,$(KERNELS))
ip:    $(addprefix ip-,$(KERNELS))

# Pull the headline numbers out of the synthesis reports.  The analytic kernel
# should show one latency; the DLS kernel should show a min/max spread.
reports:
	@echo "=================================================================="
	@echo " C-SYNTHESIS SUMMARY"
	@echo "=================================================================="
	@for k in $(KERNELS); do \
	  r=$$(find $(BUILD)/$$k -name '*_csynth.rpt' 2>/dev/null | head -1); \
	  if [ -n "$$r" ]; then \
	    echo ""; echo "--- $$k ($$r)"; \
	    sed -n '/Latency/,/^$$/p' $$r | head -20; \
	    sed -n '/Utilization Estimates/,/^====/p' $$r | head -18; \
	  else \
	    echo ""; echo "--- $$k: no report yet (run 'make syn-$$k')"; \
	  fi; \
	done

clean:
	rm -rf $(BUILD)
	@$(MAKE) -C tb clean
