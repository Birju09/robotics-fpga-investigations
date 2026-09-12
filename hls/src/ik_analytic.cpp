#include "ik_kernels.hpp"
#include "ik_math.hpp"
#include "kinematics.hpp"

//
//! Closed-form IK for a 6R arm with a spherical wrist.
//
//! Structure: decouple position from orientation at the wrist centre, solve the
//! arm as a planar 2-link problem, then read the wrist angles out of the
//! residual rotation as a ZYZ Euler triple.
//
//! Note the absence of any runtime division: the law-of-cosines denominator is
//! a geometry constant, pre-inverted in ik_config.hpp.  The only non-trivial
//! operators are five CORDIC evaluations and one square root, all fixed
//! iteration count.  That is what makes this path's latency a constant.
//
int iks::analytic(const ik_real_t Rd[3][3], const ik_real_t pd[3], int cfg,
                  ik_real_t q[IK_DOF]) {
#pragma HLS INLINE off

    const bool sh = (cfg & IK_CFG_SHOULDER) != 0;
    const bool el = (cfg & IK_CFG_ELBOW) != 0;
    const bool wr = (cfg & IK_CFG_WRIST) != 0;

    int st = IK_OK;

    //! ---- 1. wrist centre: step back along the approach axis ----
    ik_real_t pc[3];
AN_WC:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        pc[i] = (ik_real_t)((ik_acc_t)pd[i] -
                            (ik_acc_t)((ik_real_t)IK_D6 * Rd[i][2]));
    }

    //! ---- 2. base rotation, with the lateral shoulder offset ----
    //
    //! d3 translates along z2, which is parallel to z1, so joints 2 and 3 move
    //! the wrist centre in a plane held a fixed distance d3 off the joint-1
    //! axis.  That is an exact invariant, independent of q2..q6:
    //
    //!     pc_x sin(t1) - pc_y cos(t1) = d3
    //
    //! giving t1 = atan2(pc_y, pc_x) + atan2(d3, +-sqrt(rho^2 - d3^2)), the two
    //! signs being the two shoulder configurations.  Note they are NOT pi apart
    //! the way they were on the previous zero-offset arm, where this reduced to
    //! t1 = atan2(pc_y, pc_x) and its negation.
    //
    //! The radicand doubles as a reachability test the zero-offset form did not
    //! need: rho < |d3| is inside a cylinder about the joint-1 axis that the arm
    //! cannot enter.  That cylinder is the reason for the offset - without it
    //! the wrist centre can sit exactly on the axis, where t1 is undefined and
    //! atan2(0, 0) returns zero rather than failing.
    ik_real_t t1b, rho;
    ikm::atan2_hypot(pc[1], pc[0], t1b, rho);

    ik_acc_t disc_a = (ik_acc_t)(rho * rho) - (ik_acc_t)IK_D3SQ;
    if (disc_a < (ik_acc_t)0)
        st = IK_ERR_UNREACH;

    ik_real_t root = ikm::sqrt(ikm::clamp_lo0((ik_real_t)disc_a));
    ik_real_t root_s = sh ? root : (ik_real_t)(-root);

    ik_real_t t1 = ikm::wrap_pi(
        (ik_real_t)((ik_acc_t)t1b +
                    (ik_acc_t)ikm::atan2((ik_real_t)IK_D3, root_s)));

    //! ---- 3. planar 2-link sub-problem ----
    //! u = a2 cos(t2) + L3 cos(beta),  w = a2 sin(t2) + L3 sin(beta),
    //! with beta = t2 + t3 + PHI and gamma = beta - t2.
    //
    //! u carries the shoulder branch implicitly - the lefty solution simply
    //! produces a negative u - so unlike the zero-offset version there is no
    //! separate sign to apply to the radius.
    ik_real_t c1, s1;
    ikm::sincos(t1, s1, c1);
    ik_real_t u = (ik_real_t)((ik_acc_t)((ik_real_t)(pc[0] * c1)) +
                              (ik_acc_t)((ik_real_t)(pc[1] * s1)) -
                              (ik_acc_t)IK_A1);
    ik_real_t s = (ik_real_t)((ik_acc_t)pc[2] - (ik_acc_t)IK_D1);
    ik_real_t r = u;  //! kept as `r` below so the planar solve reads unchanged

    ik_acc_t rs =
        (ik_acc_t)(r * r) + (ik_acc_t)(s * s) - (ik_acc_t)IK_A2SQ_PLUS_L3SQ;
    ik_real_t cg_raw = (ik_real_t)(rs * (ik_acc_t)IK_INV_2A2L3);

    if (cg_raw > (ik_real_t)1 || cg_raw < (ik_real_t)-1)
        st = IK_ERR_UNREACH;

    ik_real_t cg = ikm::clamp1(cg_raw);
    ik_real_t sg_mag =
        ikm::sqrt((ik_real_t)((ik_acc_t)1 - (ik_acc_t)(cg * cg)));
    ik_real_t sg = el ? sg_mag : (ik_real_t)(-sg_mag);

    ik_real_t gamma = ikm::atan2(sg, cg);
    //! t3 = gamma - PHI, not gamma + PHI: alpha3 is negative on this arm, which
    //! flips the sense of the elbow offset.  See ik_geometry.hpp.
    ik_real_t t3 =
        ikm::wrap_pi((ik_real_t)((ik_acc_t)gamma - (ik_acc_t)IK_PHI));

    ik_real_t num = (ik_real_t)((ik_real_t)IK_L3 * sg);
    ik_real_t den =
        (ik_real_t)((ik_acc_t)IK_A2 + (ik_acc_t)((ik_real_t)IK_L3 * cg));
    ik_real_t t2 = ikm::wrap_pi((ik_real_t)((ik_acc_t)ikm::atan2(s, r) -
                                            (ik_acc_t)ikm::atan2(num, den)));

    //! ---- 4. wrist: R_3_6 = R_0_3^T * R_desired is a ZYZ Euler rotation ----
    ik_real_t R03[3][3];
    ikk::rot03(t1, t2, t3, R03);

    ik_real_t R36[3][3];
AN_R36_R:
    //! Rolled - a fully unrolled 3x3x3 costs 27 multipliers with zero
    //! reuse for a one-shot computation; see ik_math.hpp's sincos().
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ik_acc_t acc = (ik_acc_t)0;
            for (int k = 0; k < 3; k++) {
                acc += (ik_acc_t)(R03[k][i] * Rd[k][j]);  //! transpose of R03
            }
            R36[i][j] = (ik_real_t)acc;
        }
    }

    ik_real_t s5 = ikm::hypot(R36[0][2], R36[1][2]);
    ik_real_t t4, t5, t6;

    if (s5 < (ik_real_t)IK_WRIST_EPS) {
        //! Wrist singularity: theta4 and theta6 are no longer independent, only
        //! their sum (t5 = 0) or difference (t5 = pi) is observable.  Pin
        //! theta4 at zero and put the whole rotation on theta6.
        bool up = (R36[2][2] > (ik_real_t)0);
        t4 = (ik_real_t)0;
        t5 = up ? (ik_real_t)0 : (ik_real_t)IK_PI;
        t6 = up ? ikm::atan2((ik_real_t)(-R36[0][1]), R36[0][0])
                : ikm::atan2(R36[0][1], (ik_real_t)(-R36[0][0]));
    } else if (wr) {
        //! Negated arguments relative to the textbook ZYZ extraction: with
        //! alpha4 = +pi/2 and alpha5 = -pi/2 the residual rotation is
        //! Rz(t4) Ry(-t5) Rz(t6), so the middle angle enters negated and t4/t6
        //! each pick up a sign flip.  Getting this wrong is silent - it returns
        //! two joints wrong by pi that still look like plausible angles - which
        //! is why the testbenches round-trip through fk() instead.
        t4 = ikm::atan2((ik_real_t)(-R36[1][2]), (ik_real_t)(-R36[0][2]));
        t5 = ikm::atan2(s5, R36[2][2]);
        t6 = ikm::atan2((ik_real_t)(-R36[2][1]), R36[2][0]);
    } else {
        t4 = ikm::wrap_pi((ik_real_t)((ik_acc_t)ikm::atan2(
                                          (ik_real_t)(-R36[1][2]),
                                          (ik_real_t)(-R36[0][2])) +
                                      (ik_acc_t)IK_PI));
        t5 = ikm::atan2((ik_real_t)(-s5), R36[2][2]);
        t6 = ikm::wrap_pi(ikm::atan2(R36[2][1], (ik_real_t)(-R36[2][0])));
    }

    q[0] = t1;
    q[1] = t2;
    q[2] = t3;
    q[3] = t4;
    q[4] = t5;
    q[5] = t6;
    return st;
}

//! ------------------------------------------------------------------
//! Standalone IP wrapper
//! ------------------------------------------------------------------
extern "C" void ik_analytic_kernel(const ik_word_t pose[IK_DOF], int cfg,
                                   ik_word_t q[IK_DOF], int* status) {
#pragma HLS INTERFACE s_axilite port = pose bundle = CTRL
#pragma HLS INTERFACE s_axilite port = cfg bundle = CTRL
#pragma HLS INTERFACE s_axilite port = q bundle = CTRL
#pragma HLS INTERFACE s_axilite port = status bundle = CTRL
#pragma HLS INTERFACE s_axilite port = return bundle = CTRL

    ik_real_t pd[3], rpy[3];
AN_IN:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        pd[i] = ik_from_word(pose[i]);
        rpy[i] = ik_from_word(pose[3 + i]);
    }

    ik_real_t Rd[3][3];
    ikk::rpy_to_rot(rpy[0], rpy[1], rpy[2], Rd);

    ik_real_t qr[IK_DOF];
    int st = iks::analytic(Rd, pd, cfg, qr);

AN_OUT:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        q[i] = ik_to_word(qr[i]);
    }
    *status = st;
}
