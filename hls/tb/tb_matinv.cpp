#include "matinv.hpp"
#include "tb_common.hpp"

int main() {
    std::ifstream f = tb_open("matinv.txt");
    int nvec = 0;
    f >> nvec;

    TbStats abs_st("mat_inv (vs golden inverse)");
    TbStats res_st("mat_inv (residual A*Ainv-I)");

    //! The inverse itself is compared loosely: mi::invert reciprocates the
    //! pivot once and multiplies, where the golden model divides exactly, and
    //! the entries reach ~10^3 so a relative slip shows up as a large absolute
    //! one.  A*Ainv-I is the meaningful accuracy check and is held tight.
    const double TOL_ABS = 5e-2;
    const double TOL_RES = 5e-3;

    for (int v = 0; v < nvec; v++) {
        int n;
        f >> n;

        ik_word_t A[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t Ainv[IK_MAT_MAX * IK_MAT_MAX] = {0};
        std::vector<double> Ad(n * n), gold(n * n);

        for (int i = 0; i < n * n; i++) {
            long w;
            f >> w;
            A[(i / n) * IK_MAT_MAX + (i % n)] = (ik_word_t)w;
            Ad[i] = tb_dbl((ik_word_t)w);
        }
        for (int i = 0; i < n * n; i++) {
            long w;
            f >> w;
            gold[i] = tb_dbl((ik_word_t)w);
        }

        int status = -1;
        mat_inv_kernel(n, A, Ainv, &status);
        if (status != IK_OK) {
            std::printf("  vector %d (n=%d): unexpected status %d\n", v, n,
                        status);
            abs_st.note(1e9, TOL_ABS);
            continue;
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                abs_st.note(std::fabs(tb_dbl(Ainv[i * IK_MAT_MAX + j]) -
                                      gold[i * n + j]),
                            TOL_ABS);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double acc = 0.0;
                for (int k = 0; k < n; k++)
                    acc += Ad[i * n + k] * tb_dbl(Ainv[k * IK_MAT_MAX + j]);
                res_st.note(std::fabs(acc - (i == j ? 1.0 : 0.0)), TOL_RES);
            }
    }

    //! A rank-deficient matrix must be reported, not quietly inverted.
    {
        ik_word_t A[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t Ai[IK_MAT_MAX * IK_MAT_MAX] = {0};
        int status = -1;
        A[0] = 1 << 16;
        A[1] = 2 << 16;  //! [[1,2],[2,4]]
        A[IK_MAT_MAX + 0] = 2 << 16;
        A[IK_MAT_MAX + 1] = 4 << 16;
        mat_inv_kernel(2, A, Ai, &status);
        if (status != IK_ERR_SINGULAR) {
            std::printf("  singular matrix not detected (status=%d)  FAIL\n",
                        status);
            return 1;
        }
    }

    std::printf("[tb_matinv]\n");
    int rc = abs_st.report(TOL_ABS);
    rc |= res_st.report(TOL_RES);
    return rc;
}
