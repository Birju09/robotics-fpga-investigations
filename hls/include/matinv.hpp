#ifndef MATINV_HPP
#define MATINV_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Dense matrix inversion by Gauss-Jordan elimination with partial pivoting.
//
//! Gauss-Jordan (not Cholesky) since this is exported as a general-purpose
//! IP: handles any invertible matrix, needs only a divider, no square root.
//
//! ik_dls no longer calls this (see spd::solve() in spd.hpp, which is
//! cheaper for the SPD case DLS needs). This kernel remains the
//! general-purpose primitive and the reference spd::solve() is
//! cross-checked against in hls/tb/tb_spd.cpp.
//

namespace mi {

    //! Inverts the leading n×n block of A. Sub-6×6 requests are padded to
    //! 6×6 with identity in the unused diagonal, keeping loop trip counts
    //! (and latency) constant regardless of n.
    //!
    //! @param A    IK_MAT_MAX x IK_MAT_MAX
    //! @param n    dimension of leading block to invert (≤ IK_MAT_MAX)
    //! @param Ainv output, IK_MAT_MAX x IK_MAT_MAX
    //!
    //! @return IK_ERR_SINGULAR if singular or pivots fall below IK_PIVOT_EPS
    int invert(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
               ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX]);

}  // namespace mi

//! @defgroup MatInvIP Matrix Inversion IP
//! Standalone AXI4-Lite mapped kernel for matrix inversion.
//! @{

//! Standalone AXI4-Lite kernel: inverts an n×n block. Matrix elements are
//! flattened Q16.16 words, row-major.
//!
//! @param n      dimension (≤ IK_MAT_MAX)
//! @param status 0 on success, non-zero if singular
extern "C" void mat_inv_kernel(int n,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t Ainv[IK_MAT_MAX * IK_MAT_MAX],
                               int* status);
//! @}

#endif  //! MATINV_HPP
