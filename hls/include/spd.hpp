#ifndef SPD_HPP
#define SPD_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Direct solve for symmetric positive definite systems, by LDL^T.
 *
 * Why this exists alongside mi::invert(): the matrix the DLS solver deals
 * with, (J J^T + lambda^2 I), is symmetric positive definite by construction -
 * that is what the damping buys - and matinv.hpp already notes that a
 * Cholesky-family factorisation would be the cheaper and more stable choice
 * for that caller.  mi::invert() stays general-purpose and keeps its exported
 * IP; this is the specialised path ik_dls actually needs.
 *
 * The saving is not just the factorisation being cheaper.  DLS never wanted
 * an inverse in the first place - it wanted u solving A u = e - so forming
 * A^-1 and then multiplying by e spent a full 6x6 matrix product that
 * substitution does for free.  Against Gauss-Jordan this drops:
 *
 *   - the 6x12 augmented matrix (half the working set, and with it the
 *     muxing that a rolled elimination loop pays for on every access)
 *   - the pivot search and row swap (SPD needs no pivoting)
 *   - one of the three mm::multiply() products per DLS iteration
 *   - roughly two thirds of the multiplies: ~n^3/6 against ~n^3
 *
 * It does NOT drop the divisions - six either way - and at ~35 cycles each
 * those are now the dominant term.  A Newton-Raphson reciprocal seeded from a
 * small table would be the next thing to attack if this is still too slow.
 *
 * No square root, despite the family: LDL^T factors A = L D L^T with L unit
 * lower triangular and D diagonal, where plain Cholesky's A = L L^T would
 * need one per column.  Square roots are more expensive than divides in
 * fixed point and buy nothing here.
 */

namespace spd {

/*
 * Solve A u = b for the leading n x n block of A, A symmetric positive
 * definite.  Only the lower triangle of A is read; the upper is assumed to
 * mirror it.
 *
 * Returns IK_OK, or IK_ERR_SINGULAR if any pivot D[j] falls below
 * IK_PIVOT_EPS - which for a damped system means the caller passed something
 * that was not positive definite, since damping guarantees D[j] >= lambda^2.
 *
 * As in mi::invert(), a flagged pivot is replaced with 1 and the
 * factorisation runs on rather than returning early, and sub-6x6 requests are
 * padded to 6x6 with an identity.  Both keep the trip counts fixed so the
 * solver has one latency for every input - the point being that any latency
 * spread ik_dls shows is its iteration count and nothing else.
 */
int solve(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
          const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]);

} /* namespace spd */

#endif /* SPD_HPP */
