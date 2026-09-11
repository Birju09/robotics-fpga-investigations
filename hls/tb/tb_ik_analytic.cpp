#include "ik_kernels.hpp"
#include "kinematics.hpp"
#include "tb_common.hpp"

//
//! Two independent checks per vector:
//
//! 1. Pose round trip - FK(IK(pose)) must land back on the commanded pose,
//! in BOTH position and orientation.  This is the property a robot
//! actually depends on, and it is asserted on every vector without
//! exception.
//
//! 2. Joint agreement with the golden branch.  This is only a meaningful
//! assertion when the pose is well conditioned.  Near the elbow
//! singularity gamma = acos(cos_gamma) has slope 1/|sin(gamma)|, so the
//! +-0.5 LSB of the Q16.16 square root turns into milliradians of joint
//! error while the pose is still hit to ~1e-5 m.  The elbow-up and
//! elbow-down branches are merging there and joint space stops being
//! well determined, so those vectors are reported separately rather than
//! being allowed to fail a check they cannot pass.
//

//! Below this |sin(gamma)| joint-space agreement is reported, not asserted.
static const double COND_MIN = 0.05;

int main() {
    std::ifstream f = tb_open("ik_analytic.txt");
    int nvec = 0;
    f >> nvec;

    TbStats jt("ik_analytic (joints, cond ok)");
    TbStats pp("ik_analytic (round trip, pos m)");
    TbStats pr("ik_analytic (round trip, rot)");

    const double TOL_JOINT = 5e-3;
    const double TOL_POS = 2e-4;
    const double TOL_ROT = 2e-4;

    int n_ill = 0;
    double worst_ill = 0.0;

    for (int v = 0; v < nvec; v++) {
        ik_word_t pose[IK_DOF];
        int cfg;
        double gold[IK_DOF];
        long condw;

        for (int i = 0; i < IK_DOF; i++) {
            long w;
            f >> w;
            pose[i] = (ik_word_t)w;
        }
        f >> cfg;
        for (int i = 0; i < IK_DOF; i++) {
            long w;
            f >> w;
            gold[i] = tb_dbl((ik_word_t)w);
        }
        f >> condw;
        double cond = tb_dbl((ik_word_t)condw);

        ik_word_t q[IK_DOF];
        int status = -1;
        ik_analytic_kernel(pose, cfg, q, &status);

        if (status != IK_OK) {
            std::printf("  vector %d: status %d (expected reachable)\n", v,
                        status);
            jt.note(1e9, TOL_JOINT);
            continue;
        }

        double worst = 0.0;
        for (int i = 0; i < IK_DOF; i++)
            worst =
                std::max(worst, std::fabs(tb_angdiff(tb_dbl(q[i]), gold[i])));

        if (cond >= COND_MIN) {
            for (int i = 0; i < IK_DOF; i++)
                jt.note(std::fabs(tb_angdiff(tb_dbl(q[i]), gold[i])),
                        TOL_JOINT);
        } else {
            n_ill++;
            worst_ill = std::max(worst_ill, worst);
        }

        //! Round trip: rebuild the commanded rotation the same way the kernel
        //! does, then compare against FK of the returned joints.
        ik_real_t qr[IK_DOF];
        for (int i = 0; i < IK_DOF; i++)
            qr[i] = ik_from_word(q[i]);

        ik_real_t R[3][3], p[3];
        ikk::fk(qr, R, p);

        ik_real_t Rd[3][3];
        ikk::rpy_to_rot(ik_from_word(pose[3]), ik_from_word(pose[4]),
                        ik_from_word(pose[5]), Rd);

        double dp = 0.0, dr = 0.0;
        for (int i = 0; i < 3; i++) {
            dp += std::pow((double)p[i] - tb_dbl(pose[i]), 2.0);
            for (int j = 0; j < 3; j++)
                dr =
                    std::max(dr, std::fabs((double)R[i][j] - (double)Rd[i][j]));
        }
        pp.note(std::sqrt(dp), TOL_POS);
        pr.note(dr, TOL_ROT);
    }

    std::printf("[tb_ik_analytic]\n");
    int rc = pp.report(TOL_POS);
    rc |= pr.report(TOL_ROT);
    rc |= jt.report(TOL_JOINT);
    std::printf(
        "  %-28s %d vectors, worst joint delta %.3e (reported, not asserted)\n",
        "near elbow singularity", n_ill, worst_ill);
    return rc;
}
