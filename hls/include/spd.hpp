#ifndef SPD_HPP
#define SPD_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Direct solve for symmetric positive definite systems, by LDL^T.
//
//! Why this exists alongside mi::invert(): the matrix the DLS solver deals
//! with, (J J^T + lambda^2 I), is symmetric positive definite by construction -
//! that is what the damping buys - and matinv.hpp already notes that a
//! Cholesky-family factorisation would be the cheaper and more stable choice
//! for that caller.  mi::invert() stays general-purpose and keeps its exported
//! IP; this is the specialised path ik_dls actually needs.
//
//! The saving is not just the factorisation being cheaper.  DLS never wanted
//! an inverse in the first place - it wanted u solving A u = e - so forming
//! A^-1 and then multiplying by e spent a full 6x6 matrix product that
//! substitution does for free.  Against Gauss-Jordan this drops:
//
//! - the 6x12 augmented matrix (half the working set, and with it the
//! muxing that a rolled elimination loop pays for on every access)
//! - the pivot search and row swap (SPD needs no pivoting)
//! - one of the three mm::multiply() products per DLS iteration
//! - roughly two thirds of the multiplies: ~n^3/6 against ~n^3
//
//! It does NOT drop the divisions - six either way - and at ~35 cycles each
//! those are now the dominant term.  A Newton-Raphson reciprocal seeded from a
//! small table would be the next thing to attack if this is still too slow.
//
//! No square root, despite the family: LDL^T factors A = L D L^T with L unit
//! lower triangular and D diagonal, where plain Cholesky's A = L L^T would
//! need one per column.  Square roots are more expensive than divides in
//! fixed point and buy nothing here.
//

namespace spd {

    //! Solve a symmetric positive definite linear system using LDL^T
    //! factorization.
    //!
    //! Computes the solution u to the system A u = b where A is symmetric
    //! positive definite. Uses LDL^T factorization (avoiding the square root
    //! required by Cholesky). Only the lower triangle of A is read; the upper
    //! is assumed to mirror it for symmetry.
    //!
    //! For sub-6×6 systems, the working matrix is padded to 6×6 with an
    //! identity in the unused diagonal, ensuring all loop nests maintain
    //! constant trip counts. This results in deterministic latency independent
    //! of system dimension.
    //!
    //! @param A    Symmetric positive definite matrix (lower triangle used)
    //!             stored in IK_MAT_MAX x IK_MAT_MAX array
    //! @param n    Dimension of the n×n system block (≤ IK_MAT_MAX)
    //! @param b    Right-hand side vector (n elements)
    //! @param u    Output solution vector (n elements)
    //!
    //! @return Status code:
    //!   - IK_OK: Solution computed successfully
    //!   - IK_ERR_SINGULAR: Pivot D[j] fell below IK_PIVOT_EPS, indicating
    //!                      the matrix is not positive definite. For damped
    //!                      systems (A = JJ^T + λ²I), this indicates improper
    //!                      damping or numerical issues.
    int solve(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
              const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]);

}  // namespace spd

#endif  //! SPD_HPP
