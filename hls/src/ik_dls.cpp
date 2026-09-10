#include "ik_kernels.hpp"
#include "kinematics.hpp"
#include "matmul.hpp"
#include "matinv.hpp"
#include "ik_math.hpp"

/*
 * Damped least squares:   q <- q + J^T (J J^T + lambda^2 I)^-1 e
 *
 * Reuses mm::multiply() and mi::invert() - the same translation units the
 * standalone matrix IPs export - so the arithmetic here is the arithmetic
 * those IPs were verified against.
 *
 * The iteration loop deliberately keeps its data-dependent exit.  Bounding it
 * to a fixed trip count would make the latency constant and the comparison
 * against the analytic core meaningless; the point is to expose the jitter.
 * `max_iter` bounds the worst case so the kernel always terminates.
 */
int iks::dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
             const ik_real_t q_seed[IK_DOF],
             ik_real_t lambda, ik_real_t tol, int max_iter,
             ik_real_t q[IK_DOF], int *iters, ik_real_t *resid)
{
#pragma HLS INLINE off
    /* mm::multiply() is called three times below (A=JJ^T, U, DQ). Straight-
     * line calls to a non-inlined function are not automatically shared by
     * Vitis HLS the way calls inside a loop are - #pragma HLS ALLOCATION
     * instances=... limit=1 does not enforce this in this release either
     * (tried, verified no effect on the synthesised instance count). U and
     * DQ's calls are routed through a shared loop below instead; A=JJ^T
     * can't join them because mi::invert() must run in between. */

DLS_SEED:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = q_seed[i];
    }

    /* Both comparisons live in the Q32.32 accumulator type on purpose.  The
     * default tolerance of 1e-3 squares to 1e-6, which is two decades below
     * the Q16.16 LSB - testing it in the storage type would quantise to zero
     * and report convergence on iteration 1. */
    const ik_acc_t tol_sq = (ik_acc_t)tol * (ik_acc_t)tol;
    const ik_acc_t lam_sq = (ik_acc_t)lambda * (ik_acc_t)lambda;

    int cap = max_iter;
    if (cap < 1) cap = 1;
    if (cap > IK_DLS_MAX_ITER) cap = IK_DLS_MAX_ITER;

    int st = IK_ERR_NO_CONV;
    int used = 0;
    ik_acc_t err_sq = (ik_acc_t)0;

DLS_ITER:
    for (int k = 0; k < IK_DLS_MAX_ITER; k++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=64 avg=6
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
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
            err_sq += (ik_acc_t)(e[i] * e[i]);
        }
        if (err_sq < tol_sq) {
            st = IK_OK;
            break;
        }

        /* A = J J^T + lambda^2 I  (symmetric positive definite by
         * construction, which is what keeps the inversion well conditioned
         * even at a singularity - that is the whole point of the damping). */
        ik_real_t A[IK_MAT_MAX][IK_MAT_MAX];
        mm::multiply(J, J, IK_DOF, IK_DOF, IK_DOF, false, true, A);
DLS_DAMP:
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
            A[i][i] = (ik_real_t)((ik_acc_t)A[i][i] + lam_sq);
        }

        ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX];
        if (mi::invert(A, IK_DOF, Ainv) != IK_OK) {
            st = IK_ERR_SINGULAR;
            break;
        }

        /* u = Ainv e, then dq = J^T u.  Solving for u first keeps both
         * products at 6x6-by-6x1 instead of forming the 6x6 pseudo-inverse. */
        ik_real_t E[IK_MAT_MAX][IK_MAT_MAX];
DLS_EVEC:
        for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS UNROLL
            for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS UNROLL
                E[i][j] = (j == 0 && i < IK_DOF) ? e[i] : (ik_real_t)0;
            }
        }

        ik_real_t U[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t DQ[IK_MAT_MAX][IK_MAT_MAX];

        /* U = Ainv*E, then DQ = J^T*U: two sequential calls to
         * mm::multiply() that were synthesising as two separate 24-DSP
         * instances (three, counting A=JJ^T above) instead of sharing one -
         * an ALLOCATION limit on mm::multiply did not take effect in this
         * Vitis release. Route both through a single call site inside a
         * loop instead, the same pattern that already shares dh_step()
         * across fk()/rot03()/fk_jacobian(). */
        const ik_real_t (*mulA[2])[IK_MAT_MAX] = { Ainv, J };
        const ik_real_t (*mulB[2])[IK_MAT_MAX] = { E,    U };
        ik_real_t       (*mulC[2])[IK_MAT_MAX] = { U,    DQ };
        const bool mulTa[2] = { false, true };
UDQ:
        for (int s = 0; s < 2; s++) {
#pragma HLS PIPELINE off
            mm::multiply(mulA[s], mulB[s], IK_DOF, IK_DOF, 1,
                        mulTa[s], false, mulC[s]);
        }

        /* Wrapping each update keeps the joint state inside the CORDIC range
         * checked in ik_math.hpp.  FK is 2*pi-periodic so this cannot change
         * the trajectory, only the representative angle that comes out. */
DLS_UPD:
        for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
            q[i] = ikm::wrap_pi((ik_real_t)((ik_acc_t)q[i]
                                            + (ik_acc_t)DQ[i][0]));
        }
    }

    *iters = used;
    *resid = ikm::sqrt_acc(err_sq);
    return st;
}

/* ------------------------------------------------------------------ */
/* Standalone IP wrapper                                               */
/* ------------------------------------------------------------------ */
extern "C" void ik_dls_kernel(const ik_word_t pose[IK_DOF],
                              const ik_word_t q_seed[IK_DOF],
                              ik_word_t lambda, ik_word_t tol, int max_iter,
                              ik_word_t q[IK_DOF],
                              int *iters, ik_word_t *resid, int *status)
{
#pragma HLS INTERFACE s_axilite port=pose     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=q_seed   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=lambda   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=tol      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CTRL
#pragma HLS INTERFACE s_axilite port=q        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=iters    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=resid    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=status   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return   bundle=CTRL

    ik_real_t pd[3], rpy[3], qs[IK_DOF];
DLS_IN:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        pd[i]  = ik_from_word(pose[i]);
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
    int st = iks::dls(Rd, pd, qs,
                      ik_from_word(lambda), ik_from_word(tol), max_iter,
                      qr, &it, &rr);

DLS_OUT:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = ik_to_word(qr[i]);
    }
    *iters  = it;
    *resid  = ik_to_word(rr);
    *status = st;
}
