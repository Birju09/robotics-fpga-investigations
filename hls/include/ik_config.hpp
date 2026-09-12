#ifndef IK_CONFIG_HPP
#define IK_CONFIG_HPP

#include "ik_geometry.hpp"
#include "ik_types.hpp"

//
//! Solver limits and implementation choices.
//
//! Robot geometry lives in generated ik_geometry.hpp, not here - it used to be
//! a hand-kept copy of model/ik_model.py's constants, and a silent divergence
//! there would look like quantisation error, not a build failure. Everything
//! below is a genuine implementation choice rather than a property of the arm.
//

//! ik_geometry.hpp is generated but IS checked in, so a clone builds.  This
//! catches the case where a consumer here needs something a stale copy does
//! not define.  Regenerate with: python3 model/gen_geometry.py
#if !defined(IK_GEOMETRY_VERSION) || IK_GEOMETRY_VERSION < 1
#error \
    "hls/include/ik_geometry.hpp is missing or out of date. It is generated from the robot definition; run: python3 model/gen_geometry.py"
#endif

#define IK_MAT_MAX 6  //! max dimension handled by matmul / matinv

//! Pivot magnitude below which mat_inv declares the matrix singular.
//! 1e-4 is chosen so the reciprocal (1e4) still fits Q16.16 with margin, while
//! staying below the smallest eigenvalue that damping guarantees
//! (lambda^2 = 4e-4 at the default lambda = 0.02).
#define IK_PIVOT_EPS 1.0e-4

//! Below this |sin(theta5)| the wrist is treated as singular and theta4 is
//! pinned, because theta4 and theta6 stop being independent.
#define IK_WRIST_EPS 1.0e-4

#define IK_PI 3.1415926535897931
#define IK_TWO_PI 6.2831853071795862
#define IK_HALF_PI 1.5707963267948966

//! The DH table itself (IK_DH_A / IK_DH_D / IK_DH_CA / IK_DH_SA) is in
//! ik_geometry.hpp, generated alongside the link constants.

//! ---------------- DLS solver defaults ----------------
//! lambda = 0.02 and tol = 1e-3 come from the sweep in model/validate.py:
//! 300/300 convergence, median 4 iterations, p95 6.  Tightening tol to 1e-4
//! pushes p95 to the iteration cap because damping bounds the asymptotic rate -
//! that tail is the phenomenon this investigation is meant to expose, so it is
//! a runtime register, not a hard-coded constant.
#define IK_DLS_LAMBDA_DEFAULT 0.02
#define IK_DLS_TOL_DEFAULT 0.001
#define IK_DLS_MAX_ITER 64

//! ---------------- DLS trust region ----------------
//
//! Bound on |dq|_inf per iteration, radians. Runtime register like lambda and
//! tol - changes the iteration distribution being measured, so must be
//! sweepable without re-synthesis. 0 disables the clamp.
//
//! Why it exists.  At lambda = 0.02 the damping is nearly nil, so
//! dq = J^T (J J^T + lambda^2 I)^-1 e is essentially the full Newton step.
//! From a seed far from the solution that step is taken well outside the range
//! where the linearisation holds: instrumenting the model over far seeds
//! (model/sweep_dls.py) measured |dq|_inf reaching 12.7 rad on a 6R arm whose
//! joints span +-2.9.  The solver lands somewhere unrelated and spends its
//! iteration budget finding its way back.
//
//! Measured over ten independent 48-pose tables of the harness's own
//! composition - 480 poses, model/sweep_dls.py:
//
//!                      conv     med  p95  max  mean   miss@16
//!   no clamp        480/480       4   18   32  6.38        30
//!   clamp 1.5 rad   477/480       4   12   64  5.82        14
//
//! miss@16 is the count unsolved within a 16-iteration deadline (255 us at the
//! measured 15,952 ns/iteration).  That is the column to read: the clamp
//! roughly halves the fraction of poses that miss a deadline, and takes p95
//! from 18 iterations to 12.  The median does not move - the clamp does not
//! fire on a well-seeded solve at all - and the cost is about 3 poses in 480
//! that stop converging inside the 64-iteration cap.
//
//! On an UNSELECTED far-seeded population (300 poses, seeds +-1.2 rad) the
//! clamp is a clear win with no such cost: convergence rises from 277/300 to
//! 288/300 and the mean falls from 14.82 iterations to 11.10.
//
//! Read the two together, because the first table is biased toward the
//! baseline and the second is not.  gen_vectors.py builds its tail block by
//! keeping poses the UNCLAMPED solver converges on in 8..32 iterations, so on
//! that table the unclamped solver scores 480/480 with max=32 by construction,
//! not by merit.  Any solver change can only lose poses from a set defined by
//! its predecessor succeeding on them.  The selection has to be frozen or the
//! experiment eats itself (see gen_vectors._find_dls_tail), but frozen is not
//! the same as neutral, and the 480-pose table understates the clamp.
//
//! Radius swept over {0.5 ... 4.0} on both populations; 1.5 and 2.0 are joint
//! best and the curve is flat between them.  Tightening below ~0.7 starts
//! suppressing legitimate Newton steps and every column gets worse.  An
//! earlier draft of this comment said 1.0, fitted to a single 48-pose draw;
//! at that sample size the radius choice is inside the noise, which is why the
//! numbers above are over ten tables.
//
//! One thing the clamp cannot fix, from the same experiment: about 40% of the
//! poses that never converge are sitting at a genuine local minimum of
//! ||e||^2, not overshooting.  That needs a restart from a different seed,
//! which is a different experiment.
//
#define IK_DLS_STEP_MAX_DEFAULT 1.5

//
//! The clamp scales by a power of two, not by the exact factor R/|dq|_inf.
//
//! Two reasons, and the second is the one that matters.  It costs no
//! multiplier - a right shift on ap_fixed, at 180 of this part's 220 DSP48E1
//! slices - and it is exact in both builds, so the float reference and the
//! Q16.16 kernel follow the same trajectory through the clamp and the host
//! regression still attributes quantisation error the way it was built to.  An
//! exact R/|dq| scale would need a reciprocal, whose fixed-point and double
//! forms differ by an LSB, and the clamp would then be a second source of
//! divergence between the two builds on top of the arithmetic.
//
//! It costs a little: the effective radius lands in (R/2, R] rather than on R.
//! Measured against an exact step_max/|dq| scale, p95 and max are identical
//! and the mean differs in the second decimal.  Not worth a reciprocal.
//
//! Eight stages bounds |dq|_inf at 256*R, far above the 12.7 rad ever
//! observed; four were already enough on the measured workload.  The extra
//! four are a shallow compare chain on one scalar, not on the six components.
//
#define IK_DLS_STEP_HALVINGS 8

//
//! Adaptive lambda (Levenberg-Marquardt) was tried here and is NOT implemented,
//! which is a measurement rather than an omission.  Both forms, on a 48-pose
//! table (model/sweep_dls.py --lm):
//
//!                             conv    med  p95  max   mean
//!   baseline                 48/48      4   17   32   6.58
//!   trust region             48/48      4   10   17   5.31
//!   LM, raise lambda only    47/48      4   18   64   7.54
//!   LM, with backtracking    43/48      4   64   64  11.17
//
//! Textbook LM rejects a step that increased the residual and retries with more
//! damping.  That is a good trade when re-evaluating the residual is cheap
//! relative to a full iteration.  Here it is not: this kernel's cost is one
//! FK+Jacobian per loop pass, so a rejected step costs a whole 16 us iteration.
//! 4,995 rejections over 300 far-seeded poses, and lambda ratchets up faster
//! than the accept path brings it down, so the solver ends up crawling with
//! heavy damping - convergence fell from 48/48 to 43/48.
//
//! The non-backtracking form (raise lambda when the residual grew, keep the
//! step) avoids the wasted pass and is roughly neutral, and adds nothing at all
//! once the clamp is in.  Neither is worth the registers.
//
//! Note also what the diagnostic found about the poses that never converge:
//! about 40% of them plateau at a local minimum of ||e||^2 - a genuine
//! stationary point, not an overshoot.  No damping policy fixes that; it needs
//! a restart from a different seed, which is a different experiment.
//

//! ---------------- PL resource / latency trade ----------------
//
//! ik_dls_kernel synthesised to roughly 126% of the xc7z020's 220 DSP48E1
//! slices, so it could not be built at all - let alone alongside
//! ik_analytic_kernel, which is the bitstream the PL-vs-PL comparison
//! actually needs.
//
//! The cause is the schedule, not the algorithm.  Every loop that multiplies
//! was either fully unrolled or pipelined at II=1, and both of those oblige
//! Vitis HLS to instantiate one multiplier per concurrent product.  An
//! ap_fixed<32,16> product does not fit DSP48E1's 25x18 multiplier, so each
//! one is decomposed across several slices and the count climbs fast:
//
//! matinv  MI_NORM   12 products, fully unrolled
//! matinv  MI_ELIM   12 products at II=1
//! matmul  MM_DOT     6 products at II=1, x3 call sites
//! kinem.  JC_COLS    6 products at II=1
//! ik_dls  DLS_NORM   6 products, fully unrolled
//
//! Raising II lets HLS share multipliers across cycles instead of
//! instantiating them side by side: II=3 over six products asks for two
//! multipliers rather than six.  The trade is close to linear in both
//! directions.
//
//! Crucially it costs nothing in determinism.  Every loop listed above has a
//! fixed trip count, so a higher II lengthens every solve by the same number
//! of cycles.  The jitter this project exists to measure comes from
//! DLS_ITER's data-dependent exit, which none of this touches.  Clock rate is
//! likewise irrelevant to it: the jitter is in cycle count, and max/med is a
//! ratio that survives any uniform rescaling.
//
//! The II values are written as literals at each site rather than as macros
//! here, because #pragma lines are not macro-expanded by the C preprocessor
//! and a knob that silently failed to apply would be worse than no knob.
//! Tune them against 'make -C hls reports' - how many slices a product costs,
//! and how aggressively HLS shares them, both move between releases.  Lower
//! II means more hardware; raise until the DSP total fits, then check what it
//! cost in latency.
//
//! Sites, current values:
//! hls/src/matmul.cpp     MM_COL   II=2   ( 6 products -> 3 multipliers)
//! hls/src/kinematics.cpp JC_COLS  II=3   ( 6 products -> 2 multipliers)
//! hls/src/kinematics.cpp DH_P/R   II=1   ( 3 multipliers, shared)
//! hls/src/ik_dls.cpp     DLS_NORM rolled (II=1, 1 multiplier)
//! hls/src/spd.cpp        all reductions rolled at II=1
//! hls/src/matinv.cpp     MI_NORM  rolled (II=1, 1 multiplier)
//! hls/src/matinv.cpp     MI_ELIM  II=6   (12 products -> 2 multipliers)
//
//! matinv.cpp is no longer on the ik_dls path - spd::solve() replaced it, see
//! spd.hpp - so its two entries now only affect the standalone
//! mat_inv_kernel.  matmul.cpp is still shared with mat_mul_kernel.  Changing
//! either changes what those IPs characterise, so report the II alongside any
//! latency figure from them.
//
//! Two things that are not multiplier count at all, both now done:
//
//! ikm::sincos_batch()  hoists fk_jacobian's six CORDIC evaluations out of
//! the DH chain and pipelines them.  They never
//! depended on the chain; six serial CORDIC latencies
//! were sitting on the critical path for no reason.
//
//! ikm::recip()         replaces six ~35-cycle sequential divides per
//! spd::solve() with a table-seeded Newton-Raphson.
//
//! IK_FAST_RECIP below switches the second one off.  It is on by default, but
//! it is the one change in this file's history that the host regression
//! cannot validate - that build runs in double and never touches the
//! fixed-point path - so keep the switch until csim has confirmed it.
//

//! Table-seeded Newton-Raphson reciprocal in spd::solve().  Undefine to fall
//! back on the ap_fixed divider; see ikm::recip() in ik_math.hpp, and
//! hls/tb/tb_spd.cpp for the check that the two agree.
#define IK_FAST_RECIP 1

//
//! Pipelined CORDIC batch in fk_jacobian.  OFF by default, which is a
//! measurement and not a preference.
//
//! With all four of the latest optimisations on, ik_dls_kernel synthesised to
//! 53,493 LUT against the 53,200 this part has - over by 293, so it does not
//! place at all.  Of the four, this one costs by far the most area for the
//! least time: ikm::sincos_batch()'s PIPELINE II=1 fully unrolls the 24-stage
//! 40-bit CORDIC datapath, several thousand LUTs, to save ~114 cycles of a
//! ~1,050-cycle iteration - under 8%, and less than a third of what the
//! Newton-Raphson reciprocal buys for a fraction of the area.
//
//! So it is the first thing to drop when LUTs are short, and the first thing
//! to try again on a larger part.  On Versal, where LUTs are not the
//! constraint, turn it on.
//
//! Turning this on also turns on ikm::sincos_batch()'s unrolled engine, which
//! is shared with nothing - fk() and rot03() keep the rolled sincos() so
//! ik_analytic is unaffected either way.
//
#define IK_BATCH_CORDIC 0

//! ---------------- analytic branch selection bits ----------------
#define IK_CFG_SHOULDER 0x1  //! 1 = theta1 = atan2(pc_y, pc_x), 0 = +pi
#define IK_CFG_ELBOW 0x2     //! 1 = sin(gamma) > 0
#define IK_CFG_WRIST 0x4     //! 1 = theta5 >= 0

//! ---------------- kernel status codes ----------------
#define IK_OK 0
#define IK_ERR_UNREACH 1   //! target outside the reachable shell
#define IK_ERR_NO_CONV 2   //! DLS hit the iteration cap
#define IK_ERR_SINGULAR 3  //! matrix inversion found a null pivot
#define IK_ERR_BADDIM 4    //! requested dimension exceeds IK_MAT_MAX

#endif  //! IK_CONFIG_HPP
