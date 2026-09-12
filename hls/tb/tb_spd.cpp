#include "ik_math.hpp"
#include "matinv.hpp"
#include "spd.hpp"
#include "tb_common.hpp"

//! spd::solve() has no golden vectors: A u = b defines u completely, with no
//! convention to reimplement, so this checks the residual ||A u - b||
//! directly rather than agreeing with one reference implementation.
//! Also cross-checked against mi::invert() to catch the two solvers drifting
//! apart (which would silently change the DLS trajectory).
//
//! Inputs are A = B B^T + lambda^2 I, the shape DLS produces, swept over
//! well scaled, nearly rank deficient (kinematic-singularity-like), and
//! badly row-scaled conditioning.
//

//! dup >= 0 copies row 0 of B onto row `dup`, making B exactly rank
//! deficient, so A's eigenvalue in that direction is exactly lambda^2 -
//! a guaranteed small pivot instead of one that depends on the draw.
static void build_spd(int n, int seed, double lambda, double scale[6], int dup,
                      double A[6][6]) {
    double B[6][6];
    //! Not std::rand: this LCG reproduces identically on any machine.
    unsigned s = (unsigned)seed * 1664525u + 1013904223u;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = s * 1664525u + 1013904223u;
            B[i][j] =
                (((double)((s >> 16) & 0xFFFF) / 32768.0) - 1.0) * scale[i];
        }
    }
    if (dup > 0 && dup < n) {
        for (int j = 0; j < n; j++)
            B[dup][j] = B[0][j];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double acc = 0.0;
            for (int p = 0; p < n; p++)
                acc += B[i][p] * B[j][p];
            A[i][j] = acc + (i == j ? lambda * lambda : 0.0);
        }
    }
}

int main() {
    TbStats res_st("spd::solve (residual Au-b)");
    TbStats agr_st("spd::solve (vs mi::invert)");

    //! TOL_RES: Q16.16 LSB is 1.5e-5; a 6-term reduction rounding once per
    //! store can't beat a few LSBs. TOL_AGR is looser because spd::solve and
    //! mi::invert reciprocate/multiply at different points and disagree in
    //! the last digits by construction.
    const double TOL_RES = 2e-3;
    const double TOL_AGR = 5e-2;

    //! All lambdas here satisfy lambda^2 > IK_PIVOT_EPS (below that floor,
    //! IK_ERR_SINGULAR is the correct response, not a solution - checked
    //! separately below).
    const double lambdas[3] = {0.02, 0.2, 0.05};
    double scales[3][6] = {
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},    //! well scaled
        {1.0, 1.0, 0.02, 0.02, 1.0, 1.0},  //! near rank deficient
        {4.0, 0.25, 1.0, 4.0, 0.25, 1.0},  //! badly scaled rows
    };

    int cases = 0;
    for (int n = 2; n <= IK_MAT_MAX; n++) {
        for (int li = 0; li < 3; li++) {
            for (int si = 0; si < 3; si++) {
                for (int seed = 0; seed < 12; seed++) {
                    double Ad[6][6];
                    build_spd(n, seed + 31 * (si + 3 * li) + 971 * n,
                              lambdas[li], scales[si], -1, Ad);

                    //! Quantised through Q16.16: compare what the kernel
                    //! actually receives, not what the generator produced.
                    ik_real_t A[IK_MAT_MAX][IK_MAT_MAX];
                    ik_real_t b[IK_MAT_MAX], u[IK_MAT_MAX];
                    double Aq[6][6], bq[6];

                    for (int i = 0; i < IK_MAT_MAX; i++) {
                        for (int j = 0; j < IK_MAT_MAX; j++) {
                            double v = (i < n && j < n) ? Ad[i][j] : 0.0;
                            A[i][j] = (ik_real_t)v;
                            Aq[i][j] = (double)A[i][j];
                        }
                        double v = (i < n) ? 0.1 * (double)(i + 1) - 0.35 : 0.0;
                        b[i] = (ik_real_t)v;
                        bq[i] = (double)b[i];
                    }

                    int st = spd::solve(A, n, b, u);
                    if (st != IK_OK) {
                        std::printf(
                            "  spd::solve returned %d for n=%d "
                            "lambda=%g scale=%d seed=%d\n",
                            st, n, lambdas[li], si, seed);
                        return 1;
                    }

                    //! ---- residual ----
                    double worst_res = 0.0;
                    for (int i = 0; i < n; i++) {
                        double acc = 0.0;
                        for (int j = 0; j < n; j++)
                            acc += Aq[i][j] * (double)u[j];
                        double d = acc - bq[i];
                        if (d < 0)
                            d = -d;
                        if (d > worst_res)
                            worst_res = d;
                    }
                    res_st.note(worst_res, TOL_RES);

                    //! ---- same answer as inverting and multiplying ----
                    ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX];
                    if (mi::invert(A, n, Ainv) == IK_OK) {
                        double worst_agr = 0.0;
                        for (int i = 0; i < n; i++) {
                            double acc = 0.0;
                            for (int j = 0; j < n; j++)
                                acc += (double)Ainv[i][j] * bq[j];
                            double d = acc - (double)u[i];
                            if (d < 0)
                                d = -d;
                            if (d > worst_agr)
                                worst_agr = d;
                        }
                        agr_st.note(worst_agr, TOL_AGR);
                    }
                    cases++;
                }
            }
        }
    }

    //! A system whose damping is too weak to hold the smallest pivot above
    //! IK_PIVOT_EPS must be REPORTED (IK_ERR_SINGULAR), not silently solved -
    //! a saturated reciprocal returning IK_OK would put a silently wrong
    //! joint command on the wire.
    //
    //! Rank deficiency is exact (duplicated row) so the dependent
    //! eigenvalue is exactly lambda^2 = 4e-6, well under both the pivot
    //! floor and Q16.16's LSB. (An earlier version merely scaled two rows
    //! down and only 3/12 draws tripped - the test was wrong, not the
    //! solver; hence the exact construction here.)
    int flagged = 0, singular_cases = 0;
    for (int seed = 0; seed < 12; seed++) {
        double Ad[6][6];
        double weak[6] = {1.0, 1.0, 0.02, 0.02, 1.0, 1.0};
        build_spd(6, seed + 5000, 0.002, weak, 1, Ad);

        ik_real_t A[IK_MAT_MAX][IK_MAT_MAX];
        ik_real_t b[IK_MAT_MAX], u[IK_MAT_MAX];
        for (int i = 0; i < IK_MAT_MAX; i++) {
            for (int j = 0; j < IK_MAT_MAX; j++)
                A[i][j] = (ik_real_t)Ad[i][j];
            b[i] = (ik_real_t)(0.1 * (double)(i + 1) - 0.35);
        }

        singular_cases++;
        if (spd::solve(A, IK_MAT_MAX, b, u) == IK_ERR_SINGULAR)
            flagged++;
    }

    //! ikm::recip_q16_raw() needs its own check: spd::solve() only calls it
    //! in the fixed-point build, but this host build runs in double and
    //! takes the exact-division path, so nothing above exercises it.
    //! Swept from just above the IK_PIVOT_EPS floor (raw 7) to 20.0, past
    //! any pivot this arm's J J^T + lambda^2 I can produce.
    long recip_n = 0, recip_bad = 0;
    double recip_worst = 0.0;
    for (long d = 7; d <= 20L * 65536L; d++) {
        double exact = 65536.0 * 65536.0 / (double)d;  //! raw 1/d
        if (exact > 2147483647.0)
            continue;
        double got = (double)ikm::recip_q16_raw((int32_t)d);
        double err = got - exact;
        if (err < 0)
            err = -err;
        recip_n++;
        if (err > recip_worst)
            recip_worst = err;
        if (err > 1.0)
            recip_bad++;  //! > 1 raw LSB
    }

    std::printf("[tb_spd]\n");
    int rc = 0;
    rc |= res_st.report(TOL_RES);
    rc |= agr_st.report(TOL_AGR);
    std::printf("  %-28s %d/%d flagged IK_ERR_SINGULAR  %s\n",
                "spd::solve (underdamped)", flagged, singular_cases,
                flagged == singular_cases ? "pass" : "FAIL");
    if (flagged != singular_cases)
        rc = 1;
    std::printf("  %-28s %ld values, worst=%.2f raw LSB  %s\n",
                "ikm::recip_q16_raw", recip_n, recip_worst,
                recip_bad ? "FAIL" : "pass");
    if (recip_bad) {
        std::printf("      %ld results off by more than 1 LSB\n", recip_bad);
        rc = 1;
    }
    std::printf("  %d systems, n=2..6, 3 damping factors, 3 scalings\n", cases);
    return rc;
}
