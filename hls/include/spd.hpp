#ifndef SPD_HPP
#define SPD_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Direct solve for SPD systems via LDL^T.
//
//! Used instead of mi::invert() for the DLS path: A = J J^T + lambda^2 I is
//! SPD by construction, and DLS only ever needed u solving A u = e, not a
//! full inverse - substitution skips the matrix product Gauss-Jordan would
//! need afterward. Vs. Gauss-Jordan this also drops the augmented matrix,
//! pivoting (unneeded for SPD), and roughly 2/3 of the multiplies. Division
//! count is unchanged (six, ~35 cycles each) and is now the dominant cost;
//! a Newton-Raphson reciprocal would be the next optimization target.
//
//! LDL^T (not Cholesky) to avoid a square root per column - more expensive
//! than a divide in fixed point, for no benefit here.
//

namespace spd {

    //! Solve A u = b for SPD A via LDL^T. Only the lower triangle of A is
    //! read; upper is assumed mirrored.
    //!
    //! Sub-6x6 systems are padded to 6x6 with identity in the unused
    //! diagonal, keeping loop trip counts (and latency) constant.
    //!
    //! @param A    lower triangle used, IK_MAT_MAX x IK_MAT_MAX
    //! @param n    dimension of the n×n block (≤ IK_MAT_MAX)
    //! @param b    right-hand side (n elements)
    //! @param u    output solution (n elements)
    //!
    //! @return IK_ERR_SINGULAR if a pivot D[j] falls below IK_PIVOT_EPS
    //!         (not positive definite - for A = JJ^T + λ²I, indicates
    //!         improper damping or numerical issues).
    int solve(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
              const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]);

}  // namespace spd

#endif  //! SPD_HPP
