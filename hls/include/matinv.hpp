#ifndef MATINV_HPP
#define MATINV_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Dense matrix inversion by Gauss-Jordan elimination with partial pivoting.
 *
 * Why Gauss-Jordan and not Cholesky: the matrix the DLS solver inverts,
 * (J J^T + lambda^2 I), is symmetric positive definite, so Cholesky would be
 * the cheaper and more stable choice for that caller alone.  But this is also
 * being exported as a general-purpose IP, and Gauss-Jordan handles any
 * invertible matrix while needing only a divider - no square root.  The
 * damping already guarantees the conditioning that Cholesky would have bought.
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
