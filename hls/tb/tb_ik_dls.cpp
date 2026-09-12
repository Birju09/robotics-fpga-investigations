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
        ik_dls_kernel(pose, seed, LAMBDA, TOLW, IK_DLS_MAX_ITER, STEPW, q,
                      &iters, &resid, &status);

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

    std::printf("  converged %d/%d\n", converged, total);
    if (converged < total)
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
    return rc;
}
