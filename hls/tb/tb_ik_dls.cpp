#include <algorithm>

#include "ik_kernels.hpp"
#include "kinematics.hpp"
#include "tb_common.hpp"

//! DLS isn't checked against the golden joint vector - two solvers that both
//! converge can legitimately land on different branches. What must hold is
//! that the pose is reached.
//
//! The iteration histogram is the real deliverable: the latency
//! distribution the real-time analysis is built on.
int main() {
    std::ifstream f = tb_open("ik_dls.txt");
    int nvec = 0;
    f >> nvec;

    TbStats ps("ik_dls (FK round trip, m)");
    const double TOL_POS = 5e-3;

    //! Position-only mode solves a strictly easier problem - three task rows
    //! instead of six - so it must reach position at least as well as the
    //! full-pose solve, on the same poses and the same seeds. It says nothing
    //! about orientation, and is deliberately not checked for it.
    TbStats ps3("ik_dls (pos-only, pos err m)");

    //! task_dim defaults to IK_TASK_FULL, and that path must stay BIT-exact
    //! against the solver as it was before task relaxation existed, or every
    //! hardware number in README section 5 silently stops being reproducible.
    //! Checked here against IK_TASK_DIM_DEFAULT rather than assumed.
    long full_checked = 0, full_bad = 0;
    int conv3 = 0, baddim_ok = 0;
    std::vector<int> hist3;

    const ik_word_t LAMBDA = (ik_word_t)(IK_DLS_LAMBDA_DEFAULT * 65536.0 + 0.5);
    const ik_word_t TOLW = (ik_word_t)(IK_DLS_TOL_DEFAULT * 65536.0 + 0.5);
    const ik_word_t STEPW =
        (ik_word_t)(IK_DLS_STEP_MAX_DEFAULT * 65536.0 + 0.5);

    std::vector<int> hist;
    int converged = 0, total = 0;

    for (int v = 0; v < nvec; v++) {
        ik_word_t pose[IK_DOF], seed[IK_DOF];
        double gold[IK_DOF];
        int gold_iters;

        for (int i = 0; i < IK_DOF; i++) {
            long w;
            f >> w;
            pose[i] = (ik_word_t)w;
        }
        for (int i = 0; i < IK_DOF; i++) {
            long w;
            f >> w;
            seed[i] = (ik_word_t)w;
        }
        for (int i = 0; i < IK_DOF; i++) {
            long w;
            f >> w;
            gold[i] = tb_dbl((ik_word_t)w);
        }
        f >> gold_iters;
        (void)gold;

        ik_word_t q[IK_DOF], resid;
        int iters = -1, status = -1;
        ik_dls_kernel(pose, seed, LAMBDA, TOLW, IK_DLS_MAX_ITER, STEPW,
                      IK_TASK_DIM_DEFAULT, q, &iters, &resid, &status);

        //! ---- the default task really is the full-pose solver ----
        ik_word_t qf[IK_DOF], residf;
        int itersf = -1, statusf = -1;
        ik_dls_kernel(pose, seed, LAMBDA, TOLW, IK_DLS_MAX_ITER, STEPW,
                      IK_TASK_FULL, qf, &itersf, &residf, &statusf);
        for (int i = 0; i < IK_DOF; i++) {
            full_checked++;
            if (qf[i] != q[i])
                full_bad++;
        }
        if (itersf != iters || residf != resid || statusf != status)
            full_bad++;

        //! ---- an unsupported task dimension is refused, not rounded ----
        {
            ik_word_t qb[IK_DOF], residb;
            int itersb = -1, statusb = -1;
            for (int i = 0; i < IK_DOF; i++)
                qb[i] = 0;
            //! 5 is the mode that is planned but not built; it must fail
            //! loudly rather than quietly behave like 3 or 6.
            ik_dls_kernel(pose, seed, LAMBDA, TOLW, IK_DLS_MAX_ITER, STEPW, 5,
                          qb, &itersb, &residb, &statusb);
            bool seed_returned = true;
            for (int i = 0; i < IK_DOF; i++)
                if (qb[i] != seed[i])
                    seed_returned = false;
            if (statusb == IK_ERR_BADDIM && itersb == 0 && seed_returned)
                baddim_ok++;
        }

        //! ---- position-only: three task rows, nullspace dimension 3 ----
        ik_word_t q3[IK_DOF], resid3;
        int iters3 = -1, status3 = -1;
        ik_dls_kernel(pose, seed, LAMBDA, TOLW, IK_DLS_MAX_ITER, STEPW,
                      IK_TASK_POS, q3, &iters3, &resid3, &status3);
        if (status3 == IK_OK)
            conv3++;
        if (iters3 >= 0) {
            if ((int)hist3.size() <= iters3)
                hist3.resize(iters3 + 1, 0);
            hist3[iters3]++;
        }
        {
            ik_real_t qr3[IK_DOF];
            for (int i = 0; i < IK_DOF; i++)
                qr3[i] = ik_from_word(q3[i]);
            ik_real_t R3[3][3], p3[3];
            ikk::fk(qr3, R3, p3);
            double dp3 = 0.0;
            for (int i = 0; i < 3; i++)
                dp3 += std::pow((double)p3[i] - tb_dbl(pose[i]), 2.0);
            ps3.note(std::sqrt(dp3), TOL_POS);
        }

        total++;
        if (status == IK_OK)
            converged++;
        else {
            std::printf("  vector %d: status %d after %d iters, resid=%.3e\n",
                        v, status, iters, tb_dbl(resid));
        }

        if (iters >= 0) {
            if ((int)hist.size() <= iters)
                hist.resize(iters + 1, 0);
            hist[iters]++;
        }

        ik_real_t qr[IK_DOF];
        for (int i = 0; i < IK_DOF; i++)
            qr[i] = ik_from_word(q[i]);
        ik_real_t R[3][3], p[3];
        ikk::fk(qr, R, p);

        double dp = 0.0;
        for (int i = 0; i < 3; i++)
            dp += std::pow((double)p[i] - tb_dbl(pose[i]), 2.0);
        ps.note(std::sqrt(dp), TOL_POS);
    }

    std::printf("[tb_ik_dls]\n");
    int rc = ps.report(TOL_POS);
    rc |= ps3.report(TOL_POS);

    std::printf("  %-28s %ld values, %ld differ  %s\n",
                "default task == IK_TASK_FULL", full_checked, full_bad,
                full_bad ? "FAIL" : "pass");
    if (full_bad)
        rc = 1;
    std::printf("  %-28s %d/%d refused  %s\n", "task_dim=5 -> IK_ERR_BADDIM",
                baddim_ok, total, baddim_ok == total ? "pass" : "FAIL");
    if (baddim_ok != total)
        rc = 1;

    std::printf("  converged %d/%d full pose, %d/%d position only\n", converged,
                total, conv3, total);
    if (converged < total)
        rc = 1;
    if (conv3 < total)
        rc = 1;

    int mn = 1 << 30, mx = 0, sum = 0, cnt = 0;
    std::vector<int> all;
    for (int i = 0; i < (int)hist.size(); i++)
        for (int j = 0; j < hist[i]; j++)
            all.push_back(i);
    for (size_t i = 0; i < all.size(); i++) {
        mn = std::min(mn, all[i]);
        mx = std::max(mx, all[i]);
        sum += all[i];
        cnt++;
    }
    std::sort(all.begin(), all.end());
    if (cnt) {
        std::printf("  iterations: min=%d median=%d p95=%d max=%d mean=%.2f\n",
                    mn, all[all.size() / 2], all[(size_t)(all.size() * 0.95)],
                    mx, (double)sum / cnt);
        std::printf("  histogram:");
        for (int i = 0; i < (int)hist.size(); i++)
            if (hist[i])
                std::printf(" %d:%d", i, hist[i]);
        std::printf("\n");
    }

    //! Reported, not asserted: the iteration distribution of a reduced task
    //! is a measurement, not a requirement. A three-row task is easier per
    //! iteration to satisfy, so this should sit at or below the full-pose
    //! distribution - but it is the number the nullspace work will be
    //! compared against, so it is printed from the start.
    {
        std::vector<int> all3;
        for (int i = 0; i < (int)hist3.size(); i++)
            for (int j = 0; j < hist3[i]; j++)
                all3.push_back(i);
        std::sort(all3.begin(), all3.end());
        if (!all3.empty()) {
            int s3 = 0;
            for (size_t i = 0; i < all3.size(); i++)
                s3 += all3[i];
            std::printf(
                "  pos-only iterations: min=%d median=%d p95=%d max=%d "
                "mean=%.2f\n",
                all3.front(), all3[all3.size() / 2],
                all3[(size_t)(all3.size() * 0.95)], all3.back(),
                (double)s3 / (double)all3.size());
        }
    }
    return rc;
}
