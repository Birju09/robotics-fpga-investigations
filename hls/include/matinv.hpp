#ifndef MATINV_HPP
#define MATINV_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Dense matrix inversion by Gauss-Jordan elimination with partial pivoting.
 *
 * Why Gauss-Jordan and not Cholesky: this is exported as a general-purpose
 * IP, and Gauss-Jordan handles any invertible matrix while needing only a
 * divider - no square root.
 *
 * ik_dls no longer calls this.  It used to, and the note that used to sit
 * here said a Cholesky-family factorisation would be cheaper and more stable
 * for that caller, since (J J^T + lambda^2 I) is symmetric positive definite.
 * That is now spd::solve() - see spd.hpp - which also skips forming the
 * inverse at all, because DLS only ever wanted the solution vector.  This
 * kernel stays as the general-purpose primitive and as the reference
 * spd::solve() is cross-checked against in hls/tb/tb_spd.cpp.
 */

namespace mi {

/*
 * Ainv = A^-1 for the leading n x n block.  Returns IK_OK, or
 * IK_ERR_SINGULAR if any pivot falls below IK_PIVOT_EPS.
 *
 * Sub-6x6 requests are handled by padding the working matrix to 6x6 with an
 * identity in the unused diagonal.  Inverting the padded matrix yields the
 * wanted inverse in the leading block, so all three loop nests keep a constant
 * trip count and the IP has one latency for every dimension.
 */
int invert(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
           ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX]);

} /* namespace mi */

/* ---------------- standalone IP top level ---------------- */
extern "C" void mat_inv_kernel(int n,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t Ainv[IK_MAT_MAX * IK_MAT_MAX],
                               int *status);

#endif /* MATINV_HPP */
