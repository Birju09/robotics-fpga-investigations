#ifndef IK_CONFIG_HPP
#define IK_CONFIG_HPP

#include "ik_types.hpp"

/*
 * Robot geometry and solver limits.  These MUST stay in step with
 * model/ik_model.py - the golden vectors are generated from those values.
 *
 * 6R anthropomorphic arm, spherical wrist, standard (distal) DH:
 *
 *   i |  a_i  | alpha_i | d_i
 *   --+-------+---------+-----
 *   1 |  a1   |  +pi/2  | d1
 *   2 |  a2   |    0    | 0
 *   3 |  a3   |  +pi/2  | 0
 *   4 |  0    |  -pi/2  | d4
 *   5 |  0    |  +pi/2  | 0
 *   6 |  0    |    0    | d6
 */

#define IK_DOF     6
#define IK_MAT_MAX 6            /* max dimension handled by matmul / matinv */

/* Link geometry, metres. */
#define IK_D1 0.15
#define IK_A1 0.05
#define IK_A2 0.30
#define IK_A3 0.02
#define IK_D4 0.28
#define IK_D6 0.08

/* Derived: effective forearm length and its offset angle. */
#define IK_L3  0.28071337695236404      /* hypot(a3, d4)   */
#define IK_PHI 1.4994888620096063       /* atan2(d4, a3)   */

/* Law-of-cosines constants for the analytic solver, pre-folded so the closed
 * form needs no runtime division at all:
 *     cos(gamma) = (r^2 + s^2 - (a2^2 + L3^2)) * (1 / (2*a2*L3))          */
#define IK_A2SQ_PLUS_L3SQ 0.16880000000000001
#define IK_INV_2A2L3      5.9372541656591364

/* Pivot magnitude below which mat_inv declares the matrix singular.
 * 1e-4 is chosen so the reciprocal (1e4) still fits Q16.16 with margin, while
 * staying below the smallest eigenvalue that damping guarantees
 * (lambda^2 = 4e-4 at the default lambda = 0.02). */
#define IK_PIVOT_EPS 1.0e-4

/* Below this |sin(theta5)| the wrist is treated as singular and theta4 is
 * pinned, because theta4 and theta6 stop being independent. */
#define IK_WRIST_EPS 1.0e-4

#define IK_PI      3.1415926535897931
#define IK_TWO_PI  6.2831853071795862
#define IK_HALF_PI 1.5707963267948966

/* cos/sin of alpha_i are exactly 0 / +-1, so they are tabulated rather than
 * computed - this removes 12 CORDIC calls from every FK evaluation. */
static const double IK_DH_A[IK_DOF]  = { IK_A1, IK_A2, IK_A3, 0.0,   0.0,   0.0   };
static const double IK_DH_D[IK_DOF]  = { IK_D1, 0.0,   0.0,   IK_D4, 0.0,   IK_D6 };
static const double IK_DH_CA[IK_DOF] = { 0.0,   1.0,   0.0,   0.0,   0.0,   1.0   };
static const double IK_DH_SA[IK_DOF] = { 1.0,   0.0,   1.0,  -1.0,   1.0,   0.0   };

/* ---------------- DLS solver defaults ---------------- */
/* lambda = 0.02 and tol = 1e-3 come from the sweep in model/validate.py:
 * 300/300 convergence, median 4 iterations, p95 6.  Tightening tol to 1e-4
 * pushes p95 to the iteration cap because damping bounds the asymptotic rate -
 * that tail is the phenomenon this investigation is meant to expose, so it is
 * a runtime register, not a hard-coded constant. */
#define IK_DLS_LAMBDA_DEFAULT 0.02
#define IK_DLS_TOL_DEFAULT    0.001
#define IK_DLS_MAX_ITER       64

/* ---------------- PL resource / latency trade ---------------- */
/*
 * ik_dls_kernel synthesised to roughly 126% of the xc7z020's 220 DSP48E1
 * slices, so it could not be built at all - let alone alongside
 * ik_analytic_kernel, which is the bitstream the PL-vs-PL comparison
 * actually needs.
 *
 * The cause is the schedule, not the algorithm.  Every loop that multiplies
 * was either fully unrolled or pipelined at II=1, and both of those oblige
 * Vitis HLS to instantiate one multiplier per concurrent product.  An
 * ap_fixed<32,16> product does not fit DSP48E1's 25x18 multiplier, so each
 * one is decomposed across several slices and the count climbs fast:
 *
 *   matinv  MI_NORM   12 products, fully unrolled
 *   matinv  MI_ELIM   12 products at II=1
 *   matmul  MM_DOT     6 products at II=1, x3 call sites
 *   kinem.  JC_COLS    6 products at II=1
 *   ik_dls  DLS_NORM   6 products, fully unrolled
 *
 * Raising II lets HLS share multipliers across cycles instead of
 * instantiating them side by side: II=3 over six products asks for two
 * multipliers rather than six.  The trade is close to linear in both
 * directions.
 *
 * Crucially it costs nothing in determinism.  Every loop listed above has a
 * fixed trip count, so a higher II lengthens every solve by the same number
 * of cycles.  The jitter this project exists to measure comes from
 * DLS_ITER's data-dependent exit, which none of this touches.  Clock rate is
 * likewise irrelevant to it: the jitter is in cycle count, and max/med is a
 * ratio that survives any uniform rescaling.
 *
 * The II values are written as literals at each site rather than as macros
 * here, because #pragma lines are not macro-expanded by the C preprocessor
 * and a knob that silently failed to apply would be worse than no knob.
 * Tune them against 'make -C hls reports' - how many slices a product costs,
 * and how aggressively HLS shares them, both move between releases.  Lower
 * II means more hardware; raise until the DSP total fits, then check what it
 * cost in latency.
 *
 * Sites, current values:
 *   hls/src/matmul.cpp     MM_COL   II=2   ( 6 products -> 3 multipliers)
 *   hls/src/kinematics.cpp JC_COLS  II=3   ( 6 products -> 2 multipliers)
 *   hls/src/kinematics.cpp DH_P/R   II=1   ( 3 multipliers, shared)
 *   hls/src/ik_dls.cpp     DLS_NORM rolled (II=1, 1 multiplier)
 *   hls/src/spd.cpp        all reductions rolled at II=1
 *   hls/src/matinv.cpp     MI_NORM  rolled (II=1, 1 multiplier)
 *   hls/src/matinv.cpp     MI_ELIM  II=6   (12 products -> 2 multipliers)
 *
 * matinv.cpp is no longer on the ik_dls path - spd::solve() replaced it, see
 * spd.hpp - so its two entries now only affect the standalone
 * mat_inv_kernel.  matmul.cpp is still shared with mat_mul_kernel.  Changing
 * either changes what those IPs characterise, so report the II alongside any
 * latency figure from them.
 *
 * Two things that are not multiplier count at all, both now done:
 *
 *   ikm::sincos_batch()  hoists fk_jacobian's six CORDIC evaluations out of
 *                        the DH chain and pipelines them.  They never
 *                        depended on the chain; six serial CORDIC latencies
 *                        were sitting on the critical path for no reason.
 *
 *   ikm::recip()         replaces six ~35-cycle sequential divides per
 *                        spd::solve() with a table-seeded Newton-Raphson.
 *
 * IK_FAST_RECIP below switches the second one off.  It is on by default, but
 * it is the one change in this file's history that the host regression
 * cannot validate - that build runs in double and never touches the
 * fixed-point path - so keep the switch until csim has confirmed it.
 */

/* Table-seeded Newton-Raphson reciprocal in spd::solve().  Undefine to fall
 * back on the ap_fixed divider; see ikm::recip() in ik_math.hpp, and
 * hls/tb/tb_spd.cpp for the check that the two agree. */
#define IK_FAST_RECIP 1

/*
 * Pipelined CORDIC batch in fk_jacobian.  OFF by default, which is a
 * measurement and not a preference.
 *
 * With all four of the latest optimisations on, ik_dls_kernel synthesised to
 * 53,493 LUT against the 53,200 this part has - over by 293, so it does not
 * place at all.  Of the four, this one costs by far the most area for the
 * least time: ikm::sincos_batch()'s PIPELINE II=1 fully unrolls the 24-stage
 * 40-bit CORDIC datapath, several thousand LUTs, to save ~114 cycles of a
 * ~1,050-cycle iteration - under 8%, and less than a third of what the
 * Newton-Raphson reciprocal buys for a fraction of the area.
 *
 * So it is the first thing to drop when LUTs are short, and the first thing
 * to try again on a larger part.  On Versal, where LUTs are not the
 * constraint, turn it on.
 *
 * Turning this on also turns on ikm::sincos_batch()'s unrolled engine, which
 * is shared with nothing - fk() and rot03() keep the rolled sincos() so
 * ik_analytic is unaffected either way.
 */
#define IK_BATCH_CORDIC 0

/* ---------------- analytic branch selection bits ---------------- */
#define IK_CFG_SHOULDER 0x1     /* 1 = theta1 = atan2(pc_y, pc_x), 0 = +pi   */
#define IK_CFG_ELBOW    0x2     /* 1 = sin(gamma) > 0                        */
#define IK_CFG_WRIST    0x4     /* 1 = theta5 >= 0                           */

/* ---------------- kernel status codes ---------------- */
#define IK_OK            0
#define IK_ERR_UNREACH   1      /* target outside the reachable shell        */
#define IK_ERR_NO_CONV   2      /* DLS hit the iteration cap                 */
#define IK_ERR_SINGULAR  3      /* matrix inversion found a null pivot       */
#define IK_ERR_BADDIM    4      /* requested dimension exceeds IK_MAT_MAX    */

#endif /* IK_CONFIG_HPP */
