#include "matinv.hpp"

#define AUG (2 * IK_MAT_MAX)

static inline ik_work_t mi_abs(ik_work_t v) {
#pragma HLS INLINE
    return v < (ik_work_t)0 ? (ik_work_t)(-v) : v;
}

int mi::invert(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], int n,
               ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX]) {
#pragma HLS INLINE off

    //! Augmented working matrix [ P | I ].  Rows/cols beyond n carry an
    //! identity so the elimination below can run to a fixed 6 columns.
    //
    //! ik_work_t, not ik_real_t: the elimination stores a value per cycle and
    //! ik_real_t's AP_RND would buy a rounding adder on every one of them.
    //! Saturation is kept - see the type's note in ik_types.hpp.  The IP's
    //! interface is unchanged; only the intermediates truncate.
    ik_work_t M[IK_MAT_MAX][AUG];
#pragma HLS ARRAY_PARTITION variable = M complete dim = 2

MI_INIT:
    for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS PIPELINE II = 1
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS UNROLL
            bool real_cell = (i < n && j < n);
            bool pad_diag = (i >= n && i == j);
            M[i][j] = real_cell ? (ik_work_t)A[i][j]
                                : (pad_diag ? (ik_work_t)1 : (ik_work_t)0);
            M[i][IK_MAT_MAX + j] = (i == j) ? (ik_work_t)1 : (ik_work_t)0;
        }
    }

    int st = IK_OK;

MI_COL:
    for (int col = 0; col < IK_MAT_MAX; col++) {
        //! ---- partial pivot: largest magnitude at or below the diagonal ----
        int prow = col;
        ik_work_t best = mi_abs(M[col][col]);
    MI_PIV:
        for (int r = 0; r < IK_MAT_MAX; r++) {
#pragma HLS UNROLL
            if (r > col) {
                ik_work_t v = mi_abs(M[r][col]);
                if (v > best) {
                    best = v;
                    prow = r;
                }
            }
        }

        if (best < (ik_work_t)IK_PIVOT_EPS) {
            //! Flag and keep going with a unit pivot: bailing out early would
            //! make the IP's latency data-dependent, and the caller is told via
            //! status either way.
            st = IK_ERR_SINGULAR;
            M[col][col] = (ik_work_t)1;
        }

        //! ---- swap rows col <-> prow ----
    MI_SWAP:
        for (int j = 0; j < AUG; j++) {
#pragma HLS UNROLL
            ik_work_t t = M[col][j];
            M[col][j] = M[prow][j];
            M[prow][j] = t;
        }

        //! ---- normalise the pivot row ----
        ik_work_t inv_p = (ik_work_t)((ik_real_t)1 / (ik_real_t)M[col][col]);
    MI_NORM:
        //! Rolled, not unrolled.  Twelve concurrent products cost twelve
        //! multipliers; at II=1 over twelve cycles this costs one, and the
        //! twelve cycles are cheap against the elimination below.  j is now a
        //! variable index into the dim=2 partition, which resolves to a mux
        //! rather than a port conflict because only one column is touched per
        //! cycle.  See the resource/latency note in ik_config.hpp.
        for (int j = 0; j < AUG; j++) {
#pragma HLS PIPELINE II = 1
            M[col][j] = (ik_work_t)((ik_acc_t)(M[col][j] * inv_p));
        }

        //! ---- eliminate the column from every other row ----
    MI_ELIM:
        //! II=6, not 1.  The inner j loop is fully unrolled across AUG=12
        //! columns, so this loop's II decides whether that costs twelve
        //! multipliers or two shared across six cycles.  This was the single
        //! largest DSP consumer in ik_dls_kernel.  Latency goes from six
        //! cycles per column to thirty-six, constant either way.
        for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS PIPELINE II = 6
            if (i != col) {
                ik_work_t f = M[i][col];
                for (int j = 0; j < AUG; j++) {
#pragma HLS UNROLL
                    ik_acc_t d = (ik_acc_t)M[i][j] - (ik_acc_t)(f * M[col][j]);
                    M[i][j] = (ik_work_t)d;
                }
            }
        }
    }

MI_OUT:
    for (int i = 0; i < IK_MAT_MAX; i++) {
#pragma HLS PIPELINE II = 1
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS UNROLL
            Ainv[i][j] = (i < n && j < n) ? (ik_real_t)M[i][IK_MAT_MAX + j]
                                          : (ik_real_t)0;
        }
    }

    return st;
}

//! ------------------------------------------------------------------
//! Standalone IP wrapper
//! ------------------------------------------------------------------
extern "C" void mat_inv_kernel(int n,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t Ainv[IK_MAT_MAX * IK_MAT_MAX],
                               int* status) {
#pragma HLS INTERFACE s_axilite port = n bundle = CTRL
#pragma HLS INTERFACE s_axilite port = A bundle = CTRL
#pragma HLS INTERFACE s_axilite port = Ainv bundle = CTRL
#pragma HLS INTERFACE s_axilite port = status bundle = CTRL
#pragma HLS INTERFACE s_axilite port = return bundle = CTRL

    if (n < 1 || n > IK_MAT_MAX) {
        *status = IK_ERR_BADDIM;
        return;
    }

    ik_real_t Ar[IK_MAT_MAX][IK_MAT_MAX];
    ik_real_t Ir[IK_MAT_MAX][IK_MAT_MAX];

MI_LOAD:
    for (int i = 0; i < IK_MAT_MAX; i++) {
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS PIPELINE II = 1
            Ar[i][j] = ik_from_word(A[i * IK_MAT_MAX + j]);
        }
    }

    int st = mi::invert(Ar, n, Ir);

MI_STORE:
    for (int i = 0; i < IK_MAT_MAX; i++) {
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS PIPELINE II = 1
            Ainv[i * IK_MAT_MAX + j] = ik_to_word(Ir[i][j]);
        }
    }

    *status = st;
}
