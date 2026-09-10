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
