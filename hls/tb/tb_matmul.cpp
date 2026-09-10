#include "tb_common.hpp"
#include "matmul.hpp"

int main()
{
    std::ifstream f = tb_open("matmul.txt");
    int nvec = 0;
    f >> nvec;

    TbStats st("mat_mul_kernel");
    const double TOL = 4e-4;      /* a few Q16.16 LSBs over a 6-term dot product */

    for (int v = 0; v < nvec; v++) {
        int m, k, n, ta, tb;
        f >> m >> k >> n >> ta >> tb;

        int arows = ta ? k : m, acols = ta ? m : k;
        int brows = tb ? n : k, bcols = tb ? k : n;

        ik_word_t A[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t B[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t C[IK_MAT_MAX * IK_MAT_MAX] = {0};
        std::vector<double> gold(m * n);

        for (int i = 0; i < arows; i++)
            for (int j = 0; j < acols; j++) {
                long w; f >> w;
                A[i * IK_MAT_MAX + j] = (ik_word_t)w;
            }
        for (int i = 0; i < brows; i++)
            for (int j = 0; j < bcols; j++) {
                long w; f >> w;
                B[i * IK_MAT_MAX + j] = (ik_word_t)w;
            }
        for (int i = 0; i < m * n; i++) {
            long w; f >> w;
            gold[i] = tb_dbl((ik_word_t)w);
        }

        int status = -1;
        mat_mul_kernel(m, k, n, ta, tb, A, B, C, &status);
        if (status != IK_OK) {
            std::printf("  vector %d: unexpected status %d\n", v, status);
            st.note(1e9, TOL);
            continue;
        }

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                st.note(std::fabs(tb_dbl(C[i * IK_MAT_MAX + j])
                                  - gold[i * n + j]), TOL);
    }

    /* Dimension validation must be rejected, not silently clamped. */
    {
        ik_word_t A[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t B[IK_MAT_MAX * IK_MAT_MAX] = {0};
        ik_word_t C[IK_MAT_MAX * IK_MAT_MAX] = {0};
        int status = -1;
        mat_mul_kernel(IK_MAT_MAX + 1, 1, 1, 0, 0, A, B, C, &status);
        if (status != IK_ERR_BADDIM) {
            std::printf("  oversized dim not rejected (status=%d)  FAIL\n", status);
            return 1;
        }
    }

    std::printf("[tb_matmul]\n");
    return st.report(TOL);
}
