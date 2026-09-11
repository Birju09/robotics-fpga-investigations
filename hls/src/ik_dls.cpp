#include "ik_kernels.hpp"
#include "ik_math.hpp"
#include "kinematics.hpp"
#include "matmul.hpp"
#include "spd.hpp"

//
//! Damped least squares:   q <- q + J^T (J J^T + lambda^2 I)^-1 e
//
//! Reuses mm::multiply() and mi::invert() - the same translation units the
//! standalone matrix IPs export - so the arithmetic here is the arithmetic
//! those IPs were verified against.
//
//! The iteration loop deliberately keeps its data-dependent exit.  Bounding it
//! to a fixed trip count would make the latency constant and the comparison
//! against the analytic core meaningless; the point is to expose the jitter.
//! `max_iter` bounds the worst case so the kernel always terminates.
//
int iks::dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
             const ik_real_t q_seed[IK_DOF], ik_real_t lambda, ik_real_t tol,
             int max_iter, ik_real_t step_max, ik_real_t q[IK_DOF], int* iters,
             ik_real_t* resid) {
#pragma HLS INLINE off
    //! Both products per iteration (A=JJ^T and DQ=J^T u) go through the single
    //! mm::multiply() call site in the MULT loop below - see the comment
    //! there for why, and ik_config.hpp for the DSP budget this is part of.
    //! The linear solve between them is spd::solve(), not an inversion.

DLS_SEED:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = q_seed[i];
    }

    //! Both comparisons live in the Q32.32 accumulator type on purpose.  The
    //! default tolerance of 1e-3 squares to 1e-6, which is two decades below
    //! the Q16.16 LSB - testing it in the storage type would quantise to zero
    //! and report convergence on iteration 1.
    const ik_acc_t tol_sq = (ik_acc_t)tol * (ik_acc_t)tol;
    const ik_acc_t lam_sq = (ik_acc_t)lambda * (ik_acc_t)lambda;

    int cap = max_iter;
    if (cap < 1)
        cap = 1;
    if (cap > IK_DLS_MAX_ITER)
        cap = IK_DLS_MAX_ITER;

    int st = IK_ERR_NO_CONV;
    int used = 0;
    ik_acc_t err_sq = (ik_acc_t)0;

DLS_ITER:
    for (int k = 0; k < IK_DLS_MAX_ITER; k++) {
#pragma HLS LOOP_TRIPCOUNT min = 1 max = 64 avg = 6
        if (k >= cap)
            break;
        used = k + 1;

        ik_real_t Rc[3][3], pcur[3];
        ik_real_t J[IK_MAT_MAX][IK_MAT_MAX];
        ikk::fk_jacobian(q, Rc, pcur, J);

        ik_real_t e[IK_MAT_MAX];
        ikk::pose_error(Rd, pd, Rc, pcur, e);

        err_sq = (ik_acc_t)0;
    DLS_NORM:
        //! Rolled, not unrolled: six concurrent squares cost six multipliers
        //! for a six-term reduction that is nowhere near the critical path.
        //! II=1 over six cycles costs one.  The loop-carried dependency is a
        //! 64-bit add, which closes comfortably at this clock.
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE II = 1
            err_sq += (ik_acc_t)(e[i] * e[i]);
        }
        if (err_sq < tol_sq) {
            st = IK_OK;
            break;
        }

        //! Two products per iteration:
        //
        //! s=0   A  = J J^T      then + lambda^2 I, then solved for u
        //! s=1   DQ = J^T u
        //
        //! A is symmetric positive definite by construction - that is what the
        //! damping buys, and it is what keeps the system well conditioned even
        //! at a singularity.  So u solves A u = e directly by LDL^T, with no
        //! pivoting and no explicit inverse.
        //
        //! This used to form A^-1 by Gauss-Jordan and then multiply it by e,
        //! which is roughly three times the arithmetic and needed a third
        //! product (U = Ainv E) that substitution now does for free.  DLS
        //! never wanted an inverse; it wanted a solve.  See spd.hpp.
        //
        //! Both remaining products go through the one mm::multiply() call site
        //! below.  Straight-line calls to a non-inlined function are not
        //! shared by Vitis HLS the way calls inside a loop are, and #pragma
        //! HLS ALLOCATION instances=... limit=1 does not enforce it in this
        //! release either (tried, verified no effect on the synthesised
        //! instance count) - so two call sites would mean two multiplier
        //! banks.  This is the same pattern that already shares dh_step()
        //! across fk()/rot03()/fk_jacobian().
        //
        //! Vitis HLS does not support arrays of pointers ("pointer to
        //! pointer") for synthesis, so the operands are muxed by value into
        //! fixed staging buffers rather than selected by an array of
        //! pointers.  stageA is J for both stages; only stageB and the
        //! transpose flags differ.
        ik_real_t U[IK_MAT_MAX][IK_MAT_MAX];  //! u staged as a column
        ik_real_t DQ[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageA[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageB[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageC[IK_MAT_MAX][IK_MAT_MAX];

        const bool mulTa[2] = {false, true};
        const bool mulTb[2] = {true, false};
        const int mulN[2] = {IK_DOF, 1};

        bool singular = false;
    MULT:
        for (int s = 0; s < 2; s++) {
#pragma HLS PIPELINE off
        MULT_IN:
            for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
                for (int c = 0; c < IK_MAT_MAX; c++) {
#pragma HLS UNROLL
                    stageA[r][c] = J[r][c];
                    stageB[r][c] = (s == 0) ? J[r][c] : U[r][c];
                }
            }

            mm::multiply(stageA, stageB, IK_DOF, IK_DOF, mulN[s], mulTa[s],
                         mulTb[s], stageC);

            if (s == 0) {
                ik_real_t A[IK_MAT_MAX][IK_MAT_MAX];
            DLS_DAMP:
                for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
                    for (int c = 0; c < IK_MAT_MAX; c++) {
#pragma HLS UNROLL
                        A[r][c] =
                            (r == c && r < IK_DOF)
                                ? (ik_real_t)((ik_acc_t)stageC[r][c] + lam_sq)
                                : stageC[r][c];
                    }
                }

                //! Flagged rather than broken out of, matching spd::solve()'s
                //! own choice not to bail early: s=1 then computes a value
                //! DLS_ITER discards on the way out.  Breaking here would save
                //! one product on an error path that ends the solve anyway, at
                //! the cost of another data-dependent exit in a kernel whose
                //! latency story is the thing under measurement.
                ik_real_t uvec[IK_MAT_MAX];
                if (spd::solve(A, IK_DOF, e, uvec) != IK_OK)
                    singular = true;

            DLS_UVEC:
                for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
                    for (int c = 0; c < IK_MAT_MAX; c++) {
#pragma HLS UNROLL
                        U[r][c] =
                            (c == 0 && r < IK_DOF) ? uvec[r] : (ik_real_t)0;
                    }
                }
            } else {
            MULT_OUT:
                for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
                    for (int c = 0; c < IK_MAT_MAX; c++) {
#pragma HLS UNROLL
                        DQ[r][c] = stageC[r][c];
                    }
                }
            }
        }

        if (singular) {
            st = IK_ERR_SINGULAR;
            break;
        }

        //! ---- trust region ----
        //
        //! Bound |dq|_inf by halving until it fits, then apply the same shift
        //! to all six components.  Scaling rather than clipping per joint:
        //! clipping would change the step's DIRECTION, and the direction is
        //! the one thing the least-squares solve got right.  See the trust
        //! region block in ik_config.hpp for the measurements, the choice of
        //! radius, and why this is a shift and not an exact R/|dq| scale.
        //
        //! Written so the cost does not depend on whether the clamp fires:
        //! the shift count is always computed and always applied (sh = 0 when
        //! the step is already inside the radius).  A conditional step here
        //! would make the per-iteration latency data-dependent, and the whole
        //! value of this kernel's timing model - a constant ns/iteration, so
        //! that latency spread is attributable entirely to iteration count -
        //! rests on it not being.
        int sh = 0;
        if (step_max > (ik_real_t)0) {
            ik_real_t mag = (ik_real_t)0;
        DLS_STEP_MAG:
            for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
                ik_real_t a = (DQ[i][0] < (ik_real_t)0)
                                  ? (ik_real_t)(-DQ[i][0])
                                  : DQ[i][0];
                if (a > mag)
                    mag = a;
            }
        DLS_STEP_SHIFT:
            for (int h = 0; h < IK_DLS_STEP_HALVINGS; h++) {
#pragma HLS UNROLL
                if (mag > step_max) {
                    mag = ikm::halve(mag, 1);
                    sh++;
                }
            }
        }

    //! Wrapping each update keeps the joint state inside the CORDIC range
    //! checked in ik_math.hpp.  FK is 2*pi-periodic so this cannot change
    //! the trajectory, only the representative angle that comes out.
    DLS_UPD:
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
            ik_real_t dq = ikm::halve(DQ[i][0], sh);
            q[i] = ikm::wrap_pi((ik_real_t)((ik_acc_t)q[i] + (ik_acc_t)dq));
        }
    }

    *iters = used;
    *resid = ikm::sqrt_acc(err_sq);
    return st;
}

//! ------------------------------------------------------------------
//! Standalone IP wrapper
//! ------------------------------------------------------------------
extern "C" void ik_dls_kernel(const ik_word_t pose[IK_DOF],
                              const ik_word_t q_seed[IK_DOF], ik_word_t lambda,
                              ik_word_t tol, int max_iter, ik_word_t step_max,
                              ik_word_t q[IK_DOF], int* iters,
                              ik_word_t* resid, int* status) {
#pragma HLS INTERFACE s_axilite port = pose bundle = CTRL
#pragma HLS INTERFACE s_axilite port = q_seed bundle = CTRL
#pragma HLS INTERFACE s_axilite port = lambda bundle = CTRL
#pragma HLS INTERFACE s_axilite port = tol bundle = CTRL
#pragma HLS INTERFACE s_axilite port = max_iter bundle = CTRL
#pragma HLS INTERFACE s_axilite port = step_max bundle = CTRL
#pragma HLS INTERFACE s_axilite port = q bundle = CTRL
#pragma HLS INTERFACE s_axilite port = iters bundle = CTRL
#pragma HLS INTERFACE s_axilite port = resid bundle = CTRL
#pragma HLS INTERFACE s_axilite port = status bundle = CTRL
#pragma HLS INTERFACE s_axilite port = return bundle = CTRL

    ik_real_t pd[3], rpy[3], qs[IK_DOF];
DLS_IN:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        pd[i] = ik_from_word(pose[i]);
        rpy[i] = ik_from_word(pose[3 + i]);
    }
DLS_INQ:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        qs[i] = ik_from_word(q_seed[i]);
    }

    ik_real_t Rd[3][3];
    ikk::rpy_to_rot(rpy[0], rpy[1], rpy[2], Rd);

    ik_real_t qr[IK_DOF], rr;
    int it = 0;
    int st = iks::dls(Rd, pd, qs, ik_from_word(lambda), ik_from_word(tol),
                      max_iter, ik_from_word(step_max), qr, &it, &rr);

DLS_OUT:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = ik_to_word(qr[i]);
    }
    *iters = it;
    *resid = ik_to_word(rr);
    *status = st;
}
