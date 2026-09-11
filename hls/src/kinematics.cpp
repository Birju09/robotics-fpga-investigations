#include "kinematics.hpp"
#include "ik_math.hpp"

/* ------------------------------------------------------------------ */
/* One standard-DH step:  [R|p] <- [R|p] * A_i(theta)                  */
/*                                                                     */
/*        [ ct  -st*ca   st*sa   a*ct ]                                */
/* A_i =  [ st   ct*ca  -ct*sa   a*st ]                                */
/*        [ 0    sa      ca      d    ]                                */
/*                                                                     */
/* cos/sin of alpha_i are tabulated (they are exactly 0 or +-1), so    */
/* only one CORDIC call is needed per joint rather than two.           */
/* ------------------------------------------------------------------ */
static void dh_step(int i, ik_real_t th, ik_real_t R[3][3], ik_real_t p[3])
{
#pragma HLS INLINE
    ik_real_t ct, st;
    ikm::sincos(th, st, ct);

    const ik_real_t ca = (ik_real_t)IK_DH_CA[i];
    const ik_real_t sa = (ik_real_t)IK_DH_SA[i];
    const ik_real_t a  = (ik_real_t)IK_DH_A[i];
    const ik_real_t d  = (ik_real_t)IK_DH_D[i];

    ik_real_t AR[3][3];
    AR[0][0] = ct;  AR[0][1] = (ik_real_t)(-st * ca); AR[0][2] = (ik_real_t)(st * sa);
    AR[1][0] = st;  AR[1][1] = (ik_real_t)(ct * ca);  AR[1][2] = (ik_real_t)(-ct * sa);
    AR[2][0] = (ik_real_t)0; AR[2][1] = sa;           AR[2][2] = ca;

    ik_real_t ap[3];
    ap[0] = (ik_real_t)(a * ct);
    ap[1] = (ik_real_t)(a * st);
    ap[2] = d;

    /* p <- p + R * ap   (must use the pre-update R) */
    ik_real_t np[3];
DH_P:
    /* Rolled, like the CORDIC loops above: one multiplier reused nine
     * times costs far less silicon than nine instantiated in parallel,
     * and this runs once per joint in an already CORDIC-dominated latency
     * budget. */
    for (int r = 0; r < 3; r++) {
        ik_acc_t acc = (ik_acc_t)p[r];
        for (int c = 0; c < 3; c++) {
            acc += (ik_acc_t)(R[r][c] * ap[c]);
        }
        np[r] = (ik_real_t)acc;
    }

    /* R <- R * AR, rolled - see DH_P above. */
    ik_real_t nR[3][3];
DH_R:
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            ik_acc_t acc = (ik_acc_t)0;
            for (int k = 0; k < 3; k++) {
                acc += (ik_acc_t)(R[r][k] * AR[k][c]);
            }
            nR[r][c] = (ik_real_t)acc;
        }
    }

DH_WB:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        p[r] = np[r];
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            R[r][c] = nR[r][c];
        }
    }
}

/* ------------------------------------------------------------------ */
void ikk::rpy_to_rot(ik_real_t roll, ik_real_t pitch, ik_real_t yaw,
                     ik_real_t R[3][3])
{
#pragma HLS INLINE off
    ik_real_t sr, cr, sp, cp, sy, cy;
    ikm::sincos(roll,  sr, cr);
    ikm::sincos(pitch, sp, cp);
    ikm::sincos(yaw,   sy, cy);

    R[0][0] = (ik_real_t)(cy * cp);
    R[0][1] = (ik_real_t)((ik_acc_t)(cy * (ik_real_t)(sp * sr)) - (ik_acc_t)(sy * cr));
    R[0][2] = (ik_real_t)((ik_acc_t)(cy * (ik_real_t)(sp * cr)) + (ik_acc_t)(sy * sr));
    R[1][0] = (ik_real_t)(sy * cp);
    R[1][1] = (ik_real_t)((ik_acc_t)(sy * (ik_real_t)(sp * sr)) + (ik_acc_t)(cy * cr));
    R[1][2] = (ik_real_t)((ik_acc_t)(sy * (ik_real_t)(sp * cr)) - (ik_acc_t)(cy * sr));
    R[2][0] = (ik_real_t)(-sp);
    R[2][1] = (ik_real_t)(cp * sr);
    R[2][2] = (ik_real_t)(cp * cr);
}

void ikk::rot_to_rpy(const ik_real_t R[3][3], ik_real_t rpy[3])
{
#pragma HLS INLINE off
    ik_real_t cyp = ikm::hypot(R[0][0], R[1][0]);
    rpy[0] = ikm::atan2(R[2][1], R[2][2]);          /* roll  */
    rpy[1] = ikm::atan2((ik_real_t)(-R[2][0]), cyp); /* pitch */
    rpy[2] = ikm::atan2(R[1][0], R[0][0]);          /* yaw   */
}

/* ------------------------------------------------------------------ */
void ikk::fk(const ik_real_t q[IK_DOF], ik_real_t R[3][3], ik_real_t p[3])
{
#pragma HLS INLINE off
FK_I:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        p[r] = (ik_real_t)0;
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            R[r][c] = (r == c) ? (ik_real_t)1 : (ik_real_t)0;
        }
    }
FK_CHAIN:
    /* Vitis HLS auto-pipelines small loops with no explicit directive, and
     * pipelining this one would force dh_step()'s already-rolled internal
     * matrix multiply to flatten (fully unroll) to fit the schedule -
     * exactly the resource explosion the rolling fix was for. Keep it
     * sequential explicitly rather than relying on the tool's default. */
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE off
        dh_step(i, q[i], R, p);
    }
}

void ikk::rot03(ik_real_t t1, ik_real_t t2, ik_real_t t3, ik_real_t R[3][3])
{
#pragma HLS INLINE off
    ik_real_t p[3];
    ik_real_t th[3];
    th[0] = t1; th[1] = t2; th[2] = t3;

R3_I:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        p[r] = (ik_real_t)0;
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            R[r][c] = (r == c) ? (ik_real_t)1 : (ik_real_t)0;
        }
    }
R3_CHAIN:
    for (int i = 0; i < 3; i++) {   /* PIPELINE off - see fk()'s FK_CHAIN */
#pragma HLS PIPELINE off
        dh_step(i, th[i], R, p);
    }
}

/* ------------------------------------------------------------------ */
void ikk::fk_jacobian(const ik_real_t q[IK_DOF], ik_real_t R[3][3],
                      ik_real_t p[3], ik_real_t J[IK_MAT_MAX][IK_MAT_MAX])
{
#pragma HLS INLINE off
    ik_real_t Rc[3][3], pc[3];
    ik_real_t zax[IK_DOF][3], org[IK_DOF][3];
#pragma HLS ARRAY_PARTITION variable=zax complete dim=0
#pragma HLS ARRAY_PARTITION variable=org complete dim=0

JC_I:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        pc[r] = (ik_real_t)0;
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            Rc[r][c] = (r == c) ? (ik_real_t)1 : (ik_real_t)0;
        }
    }

    /* Joint i rotates about z_{i-1}, anchored at o_{i-1}: capture the frame
     * BEFORE applying step i. */
JC_CHAIN:
    for (int i = 0; i < IK_DOF; i++) {   /* PIPELINE off - see fk()'s FK_CHAIN */
#pragma HLS PIPELINE off
        zax[i][0] = Rc[0][2];
        zax[i][1] = Rc[1][2];
        zax[i][2] = Rc[2][2];
        org[i][0] = pc[0];
        org[i][1] = pc[1];
        org[i][2] = pc[2];
        dh_step(i, q[i], Rc, pc);
    }

JC_COPY:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        p[r] = pc[r];
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            R[r][c] = Rc[r][c];
        }
    }

    /* Jv_i = z_i x (p_e - o_i),  Jw_i = z_i */
JC_COLS:
    /* II=3, not 1.  Six cross-product multiplies per column at II=1 are six
     * multipliers; at II=3 they are two, shared.  fk_jacobian() is called
     * only from ik_dls, so this does not affect ik_analytic_kernel.  See the
     * resource/latency note in ik_config.hpp. */
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE II=3
        ik_real_t dx = (ik_real_t)((ik_acc_t)p[0] - (ik_acc_t)org[i][0]);
        ik_real_t dy = (ik_real_t)((ik_acc_t)p[1] - (ik_acc_t)org[i][1]);
        ik_real_t dz = (ik_real_t)((ik_acc_t)p[2] - (ik_acc_t)org[i][2]);

        J[0][i] = (ik_real_t)((ik_acc_t)(zax[i][1] * dz) - (ik_acc_t)(zax[i][2] * dy));
        J[1][i] = (ik_real_t)((ik_acc_t)(zax[i][2] * dx) - (ik_acc_t)(zax[i][0] * dz));
        J[2][i] = (ik_real_t)((ik_acc_t)(zax[i][0] * dy) - (ik_acc_t)(zax[i][1] * dx));
        J[3][i] = zax[i][0];
        J[4][i] = zax[i][1];
        J[5][i] = zax[i][2];
    }
}

/* ------------------------------------------------------------------ */
void ikk::pose_error(const ik_real_t Rd[3][3], const ik_real_t pd[3],
                     const ik_real_t Rc[3][3], const ik_real_t pc[3],
                     ik_real_t e[IK_MAT_MAX])
{
#pragma HLS INLINE off
PE_POS:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        e[r] = (ik_real_t)((ik_acc_t)pd[r] - (ik_acc_t)pc[r]);
    }

    /* eo = 0.5 * sum over the three column pairs of (current x desired),
     * rolled - see dh_step()'s DH_P/DH_R in this file. */
PE_ROT:
    for (int r = 0; r < 3; r++) {
        int a = (r + 1) % 3;
        int b = (r + 2) % 3;
        ik_acc_t acc = (ik_acc_t)0;
        for (int k = 0; k < 3; k++) {
            acc += (ik_acc_t)(Rc[a][k] * Rd[b][k]);
            acc -= (ik_acc_t)(Rc[b][k] * Rd[a][k]);
        }
        e[3 + r] = (ik_real_t)(acc * (ik_acc_t)0.5);
    }
}
