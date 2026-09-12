#ifndef MATMUL_HPP
#define MATMUL_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Matrix multiply, shared between the standalone mat_mul IP and the DLS
//! solver.
//
//! DLS calls mm::multiply() directly rather than via AXI to the standalone
//! IP: routing every 6x6 product through the PS would add hundreds of
//! cycles per iteration (loop runs 4-64x per solve) and make latency depend
//! on interconnect arbitration - the thing this investigation measures.
//

namespace mm {

    //! C = op(A) * op(B), op(X) = X^T if ta/tb set, else X.
    //!
    //! @param m,k,n  rows of op(A)/C, cols of op(A) & rows of op(B), cols of op(B)/C
    //!
    //! @note Matrices stored row-major in fixed IK_MAT_MAX x IK_MAT_MAX
    //! arrays; loops run unconditionally over the full range regardless of
    //! m/k/n, to keep latency independent of input dimensions.
    void multiply(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX],
                  const ik_real_t B[IK_MAT_MAX][IK_MAT_MAX], int m, int k,
                  int n, bool ta, bool tb, ik_real_t C[IK_MAT_MAX][IK_MAT_MAX]);

}  // namespace mm

//! @defgroup MatMulIP Matrix Multiplication IP
//! Standalone AXI4-Lite mapped kernel for matrix multiplication.
//! @{

//! Standalone AXI4-Lite kernel: C = op(A) * op(B). Matrix elements are
//! flattened Q16.16 words, row-major.
extern "C" void mat_mul_kernel(int m, int k, int n, int ta, int tb,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               const ik_word_t B[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t C[IK_MAT_MAX * IK_MAT_MAX],
                               int* status);
//! @}

#endif  //! MATMUL_HPP
