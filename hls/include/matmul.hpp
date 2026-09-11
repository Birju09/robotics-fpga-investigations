#ifndef MATMUL_HPP
#define MATMUL_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Matrix multiply, shared between the standalone mat_mul IP and the DLS
//! solver.
//
//! The DLS core calls mm::multiply() directly rather than issuing AXI
//! transactions to the standalone IP.  Routing every product of a 6x6 through
//! the PS or an AXI master would add hundreds of cycles per iteration to a loop
//! that runs 4-64 times per solve, and would make the latency depend on
//! interconnect arbitration - which is precisely the thing this investigation
//! is trying to measure cleanly.  Same source, two export targets.
//

namespace mm {

    //! Compute matrix multiplication with optional transpose: C = op(A) *
    //! op(B).
    //!
    //! Performs general matrix multiplication where op(X) = X^T if the
    //! corresponding transpose flag is set, otherwise op(X) = X.
    //!
    //! @param A    Input matrix A (stored in IK_MAT_MAX x IK_MAT_MAX array)
    //! @param B    Input matrix B (stored in IK_MAT_MAX x IK_MAT_MAX array)
    //! @param m    Number of rows in op(A) and output C
    //! @param k    Number of columns in op(A) and rows in op(B)
    //! @param n    Number of columns in op(B) and output C
    //! @param ta   If true, use transpose of A
    //! @param tb   If true, use transpose of B
    //! @param C    Output matrix C (m x n, stored in IK_MAT_MAX x IK_MAT_MAX
    //! array)
    //!
    //! @note All matrices are stored row-major in fixed IK_MAT_MAX x IK_MAT_MAX
    //! arrays;
    //!       runtime dimensions (m, k, n) specify which portion is meaningful.
    //!       The loops execute unconditionally over the full IK_MAT_MAX range
    //!       to maintain constant latency independent of input dimensions.
    void multiply(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX],
                  const ik_real_t B[IK_MAT_MAX][IK_MAT_MAX], int m, int k,
                  int n, bool ta, bool tb, ik_real_t C[IK_MAT_MAX][IK_MAT_MAX]);

}  // namespace mm

//! @defgroup MatMulIP Matrix Multiplication IP
//! Standalone AXI4-Lite mapped kernel for matrix multiplication.
//! @{

//! Standalone matrix multiplication kernel with AXI4-Lite interface.
//!
//! Performs matrix multiplication C = op(A) * op(B) as a standalone IP core.
//! All matrix elements are Q16.16 fixed-point words arranged in row-major order
//! and exposed as addressable register banks through the AXI4-Lite control
//! interface.
//!
//! @param m      Number of rows in op(A) and output C
//! @param k      Number of columns in op(A) and rows in op(B)
//! @param n      Number of columns in op(B) and output C
//! @param ta     If true, use transpose of matrix A
//! @param tb     If true, use transpose of matrix B
//! @param A      Input matrix A (flattened, IK_MAT_MAX * IK_MAT_MAX elements)
//! @param B      Input matrix B (flattened, IK_MAT_MAX * IK_MAT_MAX elements)
//! @param C      Output matrix C (flattened, IK_MAT_MAX * IK_MAT_MAX elements)
//! @param status Output status flag: 0 on success, non-zero on error
extern "C" void mat_mul_kernel(int m, int k, int n, int ta, int tb,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               const ik_word_t B[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t C[IK_MAT_MAX * IK_MAT_MAX],
                               int* status);
//! @}

#endif  //! MATMUL_HPP
