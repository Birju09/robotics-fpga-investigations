#include "tb_common.hpp"
#include "kinematics.hpp"

/*
 * Forward kinematics and the RPY <-> rotation pair.  These sit underneath both
 * solvers, so a convention slip here would show up as a confusing failure two
 * layers up - worth pinning independently.
 */
int main()
{
    std::ifstream f = tb_open("fk.txt");
    int nvec = 0;
    f >> nvec;

    TbStats rot("fk (rotation)");
    TbStats pos("fk (position, m)");
    TbStats rpy("rpy_to_rot -> rot_to_rpy");

    const double TOL_ROT = 2e-3;
    const double TOL_POS = 5e-4;
    const double TOL_RPY = 2e-3;

    for (int v = 0; v < nvec; v++) {
        ik_real_t q[IK_DOF];
        double gold[12];

        for (int i = 0; i < IK_DOF; i++) { long w; f >> w; q[i] = ik_from_word((ik_word_t)w); }
        for (int i = 0; i < 12; i++)     { long w; f >> w; gold[i] = tb_dbl((ik_word_t)w); }

        ik_real_t R[3][3], p[3];
        ikk::fk(q, R, p);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                rot.note(std::fabs((double)R[i][j] - gold[i * 4 + j]), TOL_ROT);
            pos.note(std::fabs((double)p[i] - gold[i * 4 + 3]), TOL_POS);
        }

        /* Round trip the orientation through the pose representation the PS
         * actually sends. */
        ik_real_t a[3], R2[3][3];
        ikk::rot_to_rpy(R, a);
        ikk::rpy_to_rot(a[0], a[1], a[2], R2);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                rpy.note(std::fabs((double)R2[i][j] - (double)R[i][j]), TOL_RPY);
    }

    std::printf("[tb_fk]\n");
    int rc = rot.report(TOL_ROT);
    rc |= pos.report(TOL_POS);
    rc |= rpy.report(TOL_RPY);
    return rc;
}
