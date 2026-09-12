#include "ik_kernels.hpp"
#include "ik_math.hpp"
#include "kinematics.hpp"
#include "matmul.hpp"
#include "spd.hpp"

//
//! q <- q + J^T (J J^T + lambda^2 I)^-1 e
//
//! Iteration loop keeps its data-dependent exit deliberately: a fixed trip
//! count would make latency constant and the comparison against the
//! analytic core meaningless. `max_iter` bounds the worst case.
//
int iks::dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
             const ik_real_t q_seed[IK_DOF], ik_real_t lambda, ik_real_t tol,
             int max_iter, ik_real_t step_max, int task_dim,
             ik_real_t q[IK_DOF], int* iters, ik_real_t* resid) {
#pragma HLS INLINE off
    //! Both products per iteration (A=JJ^T, DQ=J^T u) share the single
    //! mm::multiply() call site in MULT below - see ik_config.hpp for the
    //! DSP budget. The solve between them is spd::solve(), not an inversion.

    const int td = task_dim;

DLS_SEED:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = q_seed[i];
    }

    //! Rejected, not clamped: a task_dim the solver cannot honour must fail
    //! loudly rather than silently solve a different problem than the caller
    //! asked for. task_dim = 5 is deliberately in this set until the
    //! tool-frame rotation exists - see ik_config.hpp.
    //
    //! After DLS_SEED, not before, so a rejected call returns q = q_seed
    //! rather than whatever was in the caller's buffer. Unlike the
    //! data-dependent exits below, this one is decided by an argument rather
    //! than by the data, so it costs the measurement nothing.
    if (td != IK_TASK_FULL && td != IK_TASK_POS) {
        *iters = 0;
        *resid = (ik_real_t)0;
        return IK_ERR_BADDIM;
    }

    //! Compared in the Q32.32 accumulator, not storage type: tol^2 (1e-6 at
    //! default) is two decades below the Q16.16 LSB and would quantise to
    //! zero, reporting convergence on iteration 1.
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
        //! Rolled, not unrolled: this reduction is nowhere near the critical
        //! path, so II=1 over six cycles (one multiplier) beats six.
        //
        //! Masked to the task rows rather than bounded by them: the trip
        //! count stays six whatever task_dim is, so per-iteration latency
        //! does not depend on the task. Without the mask a reduced task
        //! would never converge - the orientation error it is not solving
        //! for would hold err_sq above tol forever.
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE II = 1
            ik_acc_t t = (ik_acc_t)(e[i] * e[i]);
            err_sq += (i < td) ? t : (ik_acc_t)0;
        }
        if (err_sq < tol_sq) {
            st = IK_OK;
            break;
        }

        //! Two products per iteration, one shared mm::multiply() call site:
        //! s=0  A = J_t J_t^T (+ lambda^2 I), td x td, solved for u
        //! s=1  DQ = J_t^T u
        //
        //! J_t is the TASK Jacobian: rows 0..td-1 of J. At td = IK_TASK_POS
        //! those are the three linear-velocity rows, so the task is position
        //! only and the orientation rows are simply never read - which is
        //! why a reduced task needs no separate Jacobian and no change to
        //! ikk::fk_jacobian(). mm::multiply() takes its dimensions at
        //! runtime and visits only the m x n block asked for, so the
        //! narrowing is free; spd::factor() pads back to 6x6 with identity,
        //! so the solve costs the same either way.
        //! A is SPD by construction (that's what the damping buys), so u
        //! solves A u = e via LDL^T - no pivoting, no explicit inverse. See
        //! spd.hpp for why this replaced a Gauss-Jordan inverse.
        //
        //! Looping over a shared call site (rather than two straight-line
        //! calls) is what makes Vitis HLS reuse one multiplier bank instead
        //! of synthesizing two - same pattern as dh_step() sharing across
        //! fk()/rot03()/fk_jacobian().
        //
        //! No pointer-to-pointer arrays in HLS synthesis, so operands are
        //! muxed by value into fixed staging buffers. stageA is J for both
        //! stages; only stageB and the transpose flags differ.
        ik_real_t U[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t DQ[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageA[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageB[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t stageC[IK_MAT_MAX][IK_MAT_MAX];

        //! m and k are per-stage now rather than both IK_DOF, so the task
        //! Jacobian's shape reaches mm::multiply(). Written as const arrays
        //! indexed by s, matching mulTa/mulTb/mulN - #pragma lines are not
        //! macro-expanded, but these are ordinary data and the table form
        //! keeps the stage definitions in one place.
        const int mulM[2] = {td, IK_DOF};
        const int mulK[2] = {IK_DOF, td};
        const int mulN[2] = {td, 1};
        const bool mulTa[2] = {false, true};
        const bool mulTb[2] = {true, false};

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

            mm::multiply(stageA, stageB, mulM[s], mulK[s], mulN[s], mulTa[s],
                         mulTb[s], stageC);

            if (s == 0) {
                ik_real_t A[IK_MAT_MAX][IK_MAT_MAX];
            DLS_DAMP:
                for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
                    for (int c = 0; c < IK_MAT_MAX; c++) {
#pragma HLS UNROLL
                        A[r][c] =
                            (r == c && r < td)
                                ? (ik_real_t)((ik_acc_t)stageC[r][c] + lam_sq)
                                : stageC[r][c];
                    }
                }

                //! Flagged rather than broken out of: avoids adding another
                //! data-dependent exit to a kernel whose latency is the thing
                //! under measurement. s=1 still runs; its result is discarded.
                ik_real_t uvec[IK_MAT_MAX];
                if (spd::solve(A, td, e, uvec) != IK_OK)
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
        //! Bound |dq|_inf by halving until it fits, applying the same shift
        //! to all six components. Scaling (not per-joint clipping) preserves
        //! the step's direction. See ik_config.hpp for radius rationale.
        //
        //! Shift count is always computed and applied (sh=0 if already
        //! inside radius) so per-iteration latency stays data-independent.
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

    //! Wraps keep joint state inside CORDIC range; FK is 2*pi-periodic so
    //! this only changes the representative angle, not the trajectory.
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
                              int task_dim, ik_word_t q[IK_DOF], int* iters,
                              ik_word_t* resid, int* status) {
#pragma HLS INTERFACE s_axilite port = pose bundle = CTRL
#pragma HLS INTERFACE s_axilite port = q_seed bundle = CTRL
#pragma HLS INTERFACE s_axilite port = lambda bundle = CTRL
#pragma HLS INTERFACE s_axilite port = tol bundle = CTRL
#pragma HLS INTERFACE s_axilite port = max_iter bundle = CTRL
#pragma HLS INTERFACE s_axilite port = step_max bundle = CTRL
#pragma HLS INTERFACE s_axilite port = task_dim bundle = CTRL
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
                      max_iter, ik_from_word(step_max), task_dim, qr, &it, &rr);

DLS_OUT:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = ik_to_word(qr[i]);
    }
    *iters = it;
    *resid = ik_to_word(rr);
    *status = st;
}
