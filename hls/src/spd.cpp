#include "spd.hpp"

#include "ik_math.hpp"

//
//! A = L D L^T, then solve by substitution. See spd.hpp for why this exists
//! alongside mi::invert().
//
//! D[j]     = A[j][j] - sum_{p<j} L[j][p] * LD[j][p]
//! LD[i][j] = A[i][j] - sum_{p<j} L[i][p] * LD[j][p]      (i > j)
//! L[i][j]  = LD[i][j] / D[j]
//
//! Carrying LD[i][j] = L[i][j]*D[j] alongside L keeps every reduction a sum
//! of plain Q16.16 products, exact in the Q32.32 accumulator. The equivalent
//! sum L[i][p]*L[j][p]*D[p] form needs a three-way product and a rounding
//! step inside the reduction - where fixed-point factorisation normally
//! loses digits.
//
int spd::solve(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
               const ik_real_t b[IK_MAT_MAX], ik_real_t u[IK_MAT_MAX]) {
#pragma HLS INLINE off

    //! Fully scalarised: loops below index both dimensions with running
    //! variables (L[i][p] forward, L[p][i] back), needing more memory ports
    //! than a real array has. 36 registers + muxes is cheaper.
    ik_work_t W[IK_MAT_MAX][IK_MAT_MAX];
    ik_work_t L[IK_MAT_MAX][IK_MAT_MAX];
    ik_work_t LD[IK_MAT_MAX][IK_MAT_MAX];
    ik_work_t invD[IK_MAT_MAX];
    ik_work_t x[IK_MAT_MAX];
#pragma HLS ARRAY_PARTITION variable = W complete dim = 0
#pragma HLS ARRAY_PARTITION variable = L complete dim = 0
#pragma HLS ARRAY_PARTITION variable = LD complete dim = 0
#pragma HLS ARRAY_PARTITION variable = invD complete dim = 1
#pragma HLS ARRAY_PARTITION variable = x complete dim = 1

    int st = IK_OK;

//! Pad to 6x6 with identity so every loop runs to IK_MAT_MAX - keeps
//! latency independent of n (cf. mat_inv's MI_INIT).
SPD_PAD:
    for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS PIPELINE II = 1
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS UNROLL
            bool real_cell = (i < n && j < n);
            W[i][j] =
                real_cell ? (ik_work_t)A[i][j] : (ik_work_t)((i == j) ? 1 : 0);
            L[i][j] = (ik_work_t)0;
            LD[i][j] = (ik_work_t)0;
        }
    }

    //! ---- factorise ----
SPD_FACT:
    for (int j = 0; j < IK_MAT_MAX; j++) {
        ik_acc_t d = (ik_acc_t)W[j][j];
    SPD_D:
        //! Bounded by j, not IK_MAT_MAX: avoids ~216 wasted cycles across
        //! the factorisation vs. ~20 useful ones.
        for (int p = 0; p < j; p++) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min = 0 max = 5
            d -= (ik_acc_t)(L[j][p] * LD[j][p]);
        }

        //! Damping guarantees d >= lambda^2 for matrices ik_dls builds; this
        //! only fires for a non-SPD caller. Flag and continue with a unit
        //! pivot - see spd.hpp.
        ik_work_t dj = (ik_work_t)d;
        if (dj < (ik_work_t)IK_PIVOT_EPS) {
            st = IK_ERR_SINGULAR;
            dj = (ik_work_t)1;
        }
        //! ikm::recip(), not '/': divides were the largest remaining cost
        //! once Gauss-Jordan was gone. See ik_math.hpp.
        invD[j] = ikm::recip(dj);

    SPD_COL:
        for (int i = j + 1; i < IK_MAT_MAX; i++) {
#pragma HLS LOOP_TRIPCOUNT min = 0 max = 5
            ik_acc_t s = (ik_acc_t)W[i][j];
        SPD_S:
            for (int p = 0; p < j; p++) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min = 0 max = 5
                s -= (ik_acc_t)(L[i][p] * LD[j][p]);
            }
            ik_work_t sv = (ik_work_t)s;
            LD[i][j] = sv;
            L[i][j] = (ik_work_t)(sv * invD[j]);
        }

        L[j][j] = (ik_work_t)1;
    }

    //! ---- forward substitution:  L y = b  ----
SPD_FWD:
    for (int i = 0; i < IK_MAT_MAX; i++) {
        ik_acc_t s = (i < n) ? (ik_acc_t)b[i] : (ik_acc_t)0;
    SPD_FWD_P:
        for (int p = 0; p < i; p++) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min = 0 max = 5
            s -= (ik_acc_t)(L[i][p] * x[p]);
        }
        x[i] = (ik_work_t)s;
    }

    //! ---- diagonal:  z = y / D  ----
SPD_DIAG:
    for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS PIPELINE II = 1
        x[i] = (ik_work_t)(x[i] * invD[i]);
    }

    //! ---- back substitution:  L^T u = z  ----
SPD_BACK:
    for (int i = IK_MAT_MAX - 1; i >= 0; i--) {
        ik_acc_t s = (ik_acc_t)x[i];
    SPD_BACK_P:
        //! L^T[i][p] is L[p][i] - transpose is free, just the other index.
        for (int p = i + 1; p < IK_MAT_MAX; p++) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min = 0 max = 5
            s -= (ik_acc_t)(L[p][i] * x[p]);
        }
        x[i] = (ik_work_t)s;
    }

SPD_OUT:
    for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS UNROLL
        u[i] = (i < n) ? (ik_real_t)x[i] : (ik_real_t)0;
    }

    return st;
}
