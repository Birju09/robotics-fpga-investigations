#include "matmul.hpp"

void mm::multiply(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX],
                  const ik_real_t B[IK_MAT_MAX][IK_MAT_MAX],
                  int m, int k, int n, bool ta, bool tb,
                  ik_real_t C[IK_MAT_MAX][IK_MAT_MAX])
{
#pragma HLS INLINE off
MM_ROW:
    for (int i = 0; i < IK_MAT_MAX; i++) {
MM_COL:
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS PIPELINE II=1
            /* Q32.32 accumulator: the product of two Q16.16 values is exactly
             * Q32.32, and six of them sum without rounding.  One rounding step
             * happens on the store below, which is what makes this bit-
             * comparable to model/ik_model.py. */
            ik_acc_t acc = (ik_acc_t)0;
MM_DOT:
            for (int p = 0; p < IK_MAT_MAX; p++) {
#pragma HLS UNROLL
                if (p < k) {
                    ik_real_t a = ta ? A[p][i] : A[i][p];
                    ik_real_t b = tb ? B[j][p] : B[p][j];
                    acc += (ik_acc_t)(a * b);
                }
            }
            C[i][j] = (i < m && j < n) ? (ik_real_t)acc : (ik_real_t)0;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Standalone IP wrapper                                               */
/* ------------------------------------------------------------------ */
extern "C" void mat_mul_kernel(int m, int k, int n, int ta, int tb,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               const ik_word_t B[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t C[IK_MAT_MAX * IK_MAT_MAX],
                               int *status)
{
#pragma HLS INTERFACE s_axilite port=m      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=k      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=n      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=ta     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=tb     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=A      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=B      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=C      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=status bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL

    if (m < 1 || m > IK_MAT_MAX || k < 1 || k > IK_MAT_MAX ||
        n < 1 || n > IK_MAT_MAX) {
        *status = IK_ERR_BADDIM;
        return;
    }

    ik_real_t Ar[IK_MAT_MAX][IK_MAT_MAX];
    ik_real_t Br[IK_MAT_MAX][IK_MAT_MAX];
    ik_real_t Cr[IK_MAT_MAX][IK_MAT_MAX];
#pragma HLS ARRAY_PARTITION variable=Ar complete dim=2
#pragma HLS ARRAY_PARTITION variable=Br complete dim=1
#pragma HLS ARRAY_PARTITION variable=Br complete dim=2

MM_LOAD:
    for (int i = 0; i < IK_MAT_MAX; i++) {
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS PIPELINE II=1
            Ar[i][j] = ik_from_word(A[i * IK_MAT_MAX + j]);
            Br[i][j] = ik_from_word(B[i * IK_MAT_MAX + j]);
        }
    }

    mm::multiply(Ar, Br, m, k, n, ta != 0, tb != 0, Cr);

MM_STORE:
    for (int i = 0; i < IK_MAT_MAX; i++) {
        for (int j = 0; j < IK_MAT_MAX; j++) {
#pragma HLS PIPELINE II=1
            C[i * IK_MAT_MAX + j] = ik_to_word(Cr[i][j]);
        }
    }

    *status = IK_OK;
}
