#include "tb_common.hpp"
#include "spd.hpp"
#include "matinv.hpp"
#include "ik_math.hpp"

/*
 * spd::solve() has no golden vectors of its own, and deliberately so.
 *
 * The other testbenches compare against model/ik_model.py because those
 * kernels implement a convention (DH order, branch selection, rounding) that
 * a reimplementation could get subtly and silently wrong.  A linear solve has
 * no such convention: A u = b defines u completely.  So this checks the
 * defining property directly - the residual ||A u - b|| - on matrices built
 * the way ik_dls builds them, which is a stronger statement than agreeing
 * with one particular reference implementation.
 *
 * Two independent checks, because they fail differently:
 *
 *   residual   catches a wrong answer
 *   vs mi      catches spd::solve() and mi::invert() having drifted apart,
 *              which would mean the DLS trajectory changed when the solver
 *              was swapped underneath it
 *
 * Inputs are A = B B^T + lambda^2 I, the exact shape DLS produces, over a
 * spread of conditioning: well scaled, nearly rank deficient (which is what a
 * kinematic singularity looks like here), and badly scaled between rows.
 */

/*
 * dup >= 0 copies row 0 of B onto row `dup`, making B exactly rank deficient.
 * A = B B^T + lambda^2 I then has an eigenvalue of exactly lambda^2 in the
 * dependent direction, which is how the singular branch below gets a
 * guaranteed small pivot instead of one that depends on the draw.
 */
static void build_spd(int n, int seed, double lambda, double scale[6],
                      int dup, double A[6][6])
{
    double B[6][6];
    /* Deterministic, portable, and not std::rand: the point is that a failure
     * reproduces on someone else's machine. */
    unsigned s = (unsigned)seed * 1664525u + 1013904223u;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = s * 1664525u + 1013904223u;
            B[i][j] = (((double)((s >> 16) & 0xFFFF) / 32768.0) - 1.0) * scale[i];
        }
    }
    if (dup > 0 && dup < n) {
        for (int j = 0; j < n; j++) B[dup][j] = B[0][j];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double acc = 0.0;
            for (int p = 0; p < n; p++) acc += B[i][p] * B[j][p];
            A[i][j] = acc + (i == j ? lambda * lambda : 0.0);
        }
    }
}

int main()
{
    TbStats res_st("spd::solve (residual Au-b)");
    TbStats agr_st("spd::solve (vs mi::invert)");

    /* Residual is held to Q16.16's own resolution: 1/65536 is 1.5e-5, and a
     * six-term reduction that rounds once per store cannot do better than a
     * few LSBs.  Agreement with mi::invert is looser for the same reason
     * tb_matinv's absolute check is - the two reciprocate and multiply at
     * different points, so they disagree in the last digits by construction. */
    const double TOL_RES = 2e-3;
    const double TOL_AGR = 5e-2;

    /* Every lambda here satisfies lambda^2 > IK_PIVOT_EPS, which is the
     * precondition ik_config.hpp states for that epsilon: it is "below the
     * smallest eigenvalue that damping guarantees".  lambda = 0.002 gives
     * lambda^2 = 4e-6, two decades UNDER the 1e-4 pivot floor, so a solver
     * that reported IK_ERR_SINGULAR there would be correct and a test that
     * demanded a solution would be wrong.  That branch is checked
     * deliberately at the end instead. */
    const double lambdas[3] = { 0.02, 0.2, 0.05 };
    double scales[3][6] = {
        { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },      /* well scaled          */
        { 1.0, 1.0, 0.02, 0.02, 1.0, 1.0 },    /* near rank deficient  */
        { 4.0, 0.25, 1.0, 4.0, 0.25, 1.0 },    /* badly scaled rows    */
    };

    int cases = 0;
    for (int n = 2; n <= IK_MAT_MAX; n++) {
        for (int li = 0; li < 3; li++) {
            for (int si = 0; si < 3; si++) {
                for (int seed = 0; seed < 12; seed++) {
                    double Ad[6][6];
                    build_spd(n, seed + 31 * (si + 3 * li) + 971 * n,
                              lambdas[li], scales[si], -1, Ad);

                    /* Quantise through Q16.16 so the testbench compares what
                     * the kernel actually receives, not what the generator
                     * produced. */
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
                        std::printf("  spd::solve returned %d for n=%d "
                                    "lambda=%g scale=%d seed=%d\n",
                                    st, n, lambdas[li], si, seed);
                        return 1;
                    }

                    /* ---- residual ---- */
                    double worst_res = 0.0;
                    for (int i = 0; i < n; i++) {
                        double acc = 0.0;
                        for (int j = 0; j < n; j++)
                            acc += Aq[i][j] * (double)u[j];
                        double d = acc - bq[i];
                        if (d < 0) d = -d;
                        if (d > worst_res) worst_res = d;
                    }
                    res_st.note(worst_res, TOL_RES);

                    /* ---- same answer as inverting and multiplying ---- */
                    ik_real_t Ainv[IK_MAT_MAX][IK_MAT_MAX];
                    if (mi::invert(A, n, Ainv) == IK_OK) {
                        double worst_agr = 0.0;
                        for (int i = 0; i < n; i++) {
                            double acc = 0.0;
                            for (int j = 0; j < n; j++)
                                acc += (double)Ainv[i][j] * bq[j];
                            double d = acc - (double)u[i];
                            if (d < 0) d = -d;
                            if (d > worst_agr) worst_agr = d;
                        }
                        agr_st.note(worst_agr, TOL_AGR);
                    }
                    cases++;
                }
            }
        }
    }

    /*
     * The other branch: a system whose damping is too weak to hold the
     * smallest pivot above IK_PIVOT_EPS must be REPORTED, not silently
     * solved.  This is what the sweep above deliberately excludes.
     *
     * Reporting is the whole contract here - ik_dls turns it into
     * IK_ERR_SINGULAR and abandons the pose, so a solver that returned
     * IK_OK with a saturated reciprocal would put a silently wrong joint
     * command on the wire.
     *
     * The rank deficiency is exact (a duplicated row) rather than merely
     * ill-scaled, so the dependent direction's eigenvalue is exactly
     * lambda^2 = 4e-6 - two decades under the 1e-4 pivot floor, and below
     * Q16.16's 1.5e-5 LSB besides.  An earlier version of this check just
     * scaled two rows down by 0.02 and asserted all twelve draws would
     * trip; only three did, because four full-scale rows leave plenty of
     * rank.  That was the test being wrong, not the solver.
     */
    int flagged = 0, singular_cases = 0;
    for (int seed = 0; seed < 12; seed++) {
        double Ad[6][6];
        double weak[6] = { 1.0, 1.0, 0.02, 0.02, 1.0, 1.0 };
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

    /*
     * ikm::recip_q16_raw() gets its own check, and needs one.
     *
     * spd::solve() calls it only in the fixed-point build; this build runs in
     * double and takes the exact-division path, so nothing above would notice
     * if the Newton-Raphson were wrong.  It is integer-only and identical in
     * both builds precisely so it can be tested here - a silently wrong
     * reciprocal would corrupt every DLS solve, and there is no other
     * opportunity to catch that before hardware.
     *
     * Swept across the whole supported input range: from just above the
     * IK_PIVOT_EPS floor (raw 7) to 20.0, which is well past any pivot
     * J J^T + lambda^2 I can produce for this arm.
     */
    long recip_n = 0, recip_bad = 0;
    double recip_worst = 0.0;
    for (long d = 7; d <= 20L * 65536L; d++) {
        double exact = 65536.0 * 65536.0 / (double)d;      /* raw 1/d */
        if (exact > 2147483647.0) continue;
        double got = (double)ikm::recip_q16_raw((int32_t)d);
        double err = got - exact;
        if (err < 0) err = -err;
        recip_n++;
        if (err > recip_worst) recip_worst = err;
        if (err > 1.0) recip_bad++;                        /* > 1 raw LSB */
    }

    std::printf("[tb_spd]\n");
    int rc = 0;
    rc |= res_st.report(TOL_RES);
    rc |= agr_st.report(TOL_AGR);
    std::printf("  %-28s %d/%d flagged IK_ERR_SINGULAR  %s\n",
                "spd::solve (underdamped)", flagged, singular_cases,
                flagged == singular_cases ? "pass" : "FAIL");
    if (flagged != singular_cases) rc = 1;
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
