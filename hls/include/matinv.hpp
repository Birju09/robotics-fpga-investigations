#ifndef MATINV_HPP
#define MATINV_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Dense matrix inversion by Gauss-Jordan elimination with partial pivoting.
//
//! Why Gauss-Jordan and not Cholesky: this is exported as a general-purpose
//! IP, and Gauss-Jordan handles any invertible matrix while needing only a
//! divider - no square root.
//
//! ik_dls no longer calls this.  It used to, and the note that used to sit
//! here said a Cholesky-family factorisation would be cheaper and more stable
//! for that caller, since (J J^T + lambda^2 I) is symmetric positive definite.
//! That is now spd::solve() - see spd.hpp - which also skips forming the
//! inverse at all, because DLS only ever wanted the solution vector.  This
//! kernel stays as the general-purpose primitive and as the reference
//! spd::solve() is cross-checked against in hls/tb/tb_spd.cpp.
//

namespace mi {

    //! Compute matrix inverse using Gauss-Jordan elimination with partial
    //! pivoting.
    //!
    //! Inverts the leading n×n block of a matrix stored in an
    //! IK_MAT_MAX×IK_MAT_MAX array. Handles sub-6×6 requests by padding the
    //! working matrix to 6×6 with an identity in the unused diagonal; inverting
    //! the padded matrix yields the desired inverse in the leading block. This
    //! ensures all loop nests maintain constant trip counts, resulting in
    //! deterministic latency regardless of n.
    //!
    //! @param A    Input matrix (stored in IK_MAT_MAX x IK_MAT_MAX array)
    //! @param n    Dimension of the leading block to invert (≤ IK_MAT_MAX)
    //! @param Ainv Output inverse matrix (stored in IK_MAT_MAX x IK_MAT_MAX
    //! array)
    //!
    //! @return Status code:
    //!   - IK_OK: Inversion successful
    //!   - IK_ERR_SINGULAR: Matrix is singular or pivots fall below
    //!   IK_PIVOT_EPS
    int invert(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
               ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX]);

}  // namespace mi

//! @defgroup MatInvIP Matrix Inversion IP
//! Standalone AXI4-Lite mapped kernel for matrix inversion.
//! @{

//! Standalone matrix inversion kernel with AXI4-Lite interface.
//!
//! Computes the inverse of an n×n matrix block using Gauss-Jordan elimination
//! with partial pivoting. All matrix elements are Q16.16 fixed-point words
//! arranged in row-major order and exposed as addressable register banks
//! through the AXI4-Lite control interface.
//!
//! @param n      Dimension of the matrix block to invert (≤ IK_MAT_MAX)
//! @param A      Input matrix (flattened, IK_MAT_MAX * IK_MAT_MAX elements)
//! @param Ainv   Output inverse matrix (flattened, IK_MAT_MAX * IK_MAT_MAX
//! elements)
//! @param status Output status: 0 on success, non-zero if matrix is singular
extern "C" void mat_inv_kernel(int n,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t Ainv[IK_MAT_MAX * IK_MAT_MAX],
                               int* status);
//! @}

#endif  //! MATINV_HPP
