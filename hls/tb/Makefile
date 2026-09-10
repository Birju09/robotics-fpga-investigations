# Host regression build.
#
# Compiles the kernels in their float reference configuration (IK_USE_FLOAT),
# which needs no ap_fixed headers and therefore no Vitis installation.  This
# validates every algorithm - DH conventions, CORDIC, Gauss-Jordan, the DLS
# iteration - independently of the fixed-point question, which is settled
# separately by model/validate.py's quantisation sweep.
#
#   make        build and run everything
#   make clean
#
# For the real fixed-point verification run C-simulation:
#   make -C ../.. csim          (see hls/Makefile)

CXX      ?= g++
INC      := -I../include -I.
CXXFLAGS := -O2 -std=c++14 -Wall -Wno-unknown-pragmas -Wno-unused-label -DIK_USE_FLOAT \
            -DTB_VECTOR_DIR=\"vectors\" $(INC)

SRC   := ../src/matmul.cpp ../src/matinv.cpp ../src/kinematics.cpp \
         ../src/ik_analytic.cpp ../src/ik_dls.cpp
TESTS := tb_matmul tb_matinv tb_fk tb_ik_analytic tb_ik_dls
BIN   := $(addprefix build/,$(TESTS))

.PHONY: all run clean
all: run

build:
	@mkdir -p build

build/%: %.cpp $(SRC) | build
	$(CXX) $(CXXFLAGS) $< $(SRC) -o $@

run: $(BIN)
	@echo "=============================================================="
	@echo " HOST REGRESSION (float reference build)"
	@echo "=============================================================="
	@rc=0; for t in $(BIN); do ./$$t || rc=1; done; \
	 echo "--------------------------------------------------------------"; \
	 if [ $$rc -eq 0 ]; then echo " ALL TESTS PASSED"; else echo " FAILURES"; fi; \
	 exit $$rc

clean:
	rm -rf build
