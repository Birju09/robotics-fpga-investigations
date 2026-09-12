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
//! count is six, ~35 cycles each without IK_FAST_RECIP, and is the dominant
//! cost - which is what makes the factor/substitute split below worth having.
//
//! LDL^T (not Cholesky) to avoid a square root per column - more expensive
//! than a divide in fixed point, for no benefit here.
//
//! ---------------- why this is split three ways ----------------
//
//! factor() and substitute() are the primitives; solve() and solve_n() are
//! the entry points. The split exists because the nullspace step in
//! iks::dls() needs TWO right-hand sides against the SAME matrix
//! (A^-1 e for the task term, A^-1 (J z) for the projection - see
//! docs/collision_aware_ik.md §5), and factoring twice would pay the six
//! reciprocals twice for nothing.
//
//! Callers should use solve() or solve_n(), not the primitives directly:
//! both are `INLINE off` so the LDL^T datapath is synthesised once, and
//! solve_n() shares ONE substitution instance across its right-hand sides
//! via a loop over a single call site - the same trick ik_dls.cpp's MULT
//! loop uses to share one multiplier bank. Calling substitute() twice from
//! straight-line code would instead instantiate it twice.
//

namespace spd {

    //! Factor A = L D L^T. Only the lower triangle of A is read; upper is
    //! assumed mirrored.
    //!
    //! Sub-6x6 systems are padded to 6x6 with identity in the unused
    //! diagonal, keeping loop trip counts (and latency) constant. That
    //! padding is also what lets the reduced-task DLS modes (task_dim < 6)
    //! cost the solve exactly zero extra cycles.
    //!
    //! L and invD are ik_work_t, not ik_real_t: they never cross an IP
    //! boundary (see ik_types.hpp). Only L's strict lower triangle and its
    //! unit diagonal are meaningful - the upper triangle is left zero.
    //!
    //! ARRAY_PARTITION on L/invD is the CALLER's to declare, since these are
    //! function parameters and both primitives are `INLINE` (on). solve()
    //! and solve_n() do it; a caller reaching past them must too, or the
    //! scalarised indexing below degenerates into memory-port contention.
    //!
    //! @return IK_ERR_SINGULAR if a pivot D[j] falls below IK_PIVOT_EPS
    //!         (not positive definite - for A = JJ^T + λ²I, indicates
    //!         improper damping or numerical issues).
    int factor(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
               ik_work_t L[IK_MAT_MAX][IK_MAT_MAX],
               ik_work_t invD[IK_MAT_MAX]);

    //! Solve L D L^T u = b given factor()'s output. No divisions: invD is
    //! already reciprocated, so this is pure multiply-accumulate and is
    //! roughly a third of a full solve.
    void substitute(const ik_work_t L[IK_MAT_MAX][IK_MAT_MAX],
                    const ik_work_t invD[IK_MAT_MAX], int n,
                    const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]);

    //! Solve A u = b for SPD A. factor() + substitute(); kept as the
    //! single-right-hand-side entry point.
    //!
    //! @param A    lower triangle used, IK_MAT_MAX x IK_MAT_MAX
    //! @param n    dimension of the n×n block (≤ IK_MAT_MAX)
    //! @param b    right-hand side (n elements)
    //! @param u    output solution (n elements)
    int solve(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
              const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]);

    //! Solve A U[k] = B[k] for k = 0 .. IK_SPD_NRHS-1 against one
    //! factorisation.
    //!
    //! Every right-hand side is ALWAYS substituted, even when the caller
    //! only wants the first - latency must not depend on how many are live,
    //! for the same reason iks::dls() always computes and applies its
    //! trust-region shift and always runs the s=1 product after flagging a
    //! singular matrix. Discard what you don't need.
    //!
    //! Bit-identical to calling solve() once per right-hand side; asserted
    //! in hls/tb/tb_spd.cpp.
    int solve_n(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
                const ik_real_t B[IK_SPD_NRHS][IK_MAT_MAX],
                ik_real_t U[IK_SPD_NRHS][IK_MAT_MAX]);

}  // namespace spd

#endif  //! SPD_HPP
