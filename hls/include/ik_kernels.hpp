#ifndef IK_KERNELS_HPP
#define IK_KERNELS_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * The two inverse-kinematics solvers under comparison.
 *
 *   analytic() - closed form.  Straight-line code: no loop whose trip count
 *                depends on the data, so C-synthesis reports one latency and
 *                that latency is the worst case.
 *
 *   dls()      - damped least squares.  Iterates until the task-space residual
 *                falls below `tol` or the cap is hit, so its latency is a
 *                function of the pose, the seed and the conditioning of the
 *                Jacobian.  It reports the iteration count it actually used;
 *                that number, and its distribution over a workload, is the
 *                measurement this project exists to take.
 */

namespace iks {

/* Returns IK_OK or IK_ERR_UNREACH.  `cfg` selects one of the eight branches
 * via the IK_CFG_* bits. */
int analytic(const ik_real_t Rd[3][3], const ik_real_t pd[3], int cfg,
             ik_real_t q[IK_DOF]);

/* Returns IK_OK, IK_ERR_NO_CONV or IK_ERR_SINGULAR.
 * `iters` receives the number of iterations executed, `resid` the final
 * task-space residual (2-norm, Q16.16). */
int dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
        const ik_real_t q_seed[IK_DOF],
        ik_real_t lambda, ik_real_t tol, int max_iter,
        ik_real_t q[IK_DOF], int *iters, ik_real_t *resid);

} /* namespace iks */

/* ---------------- standalone IP top levels ---------------- */

/*
 * pose[] = { x, y, z, roll, pitch, yaw } in Q16.16, with
 * R = Rz(yaw) Ry(pitch) Rx(roll).  q[] receives six joint angles in radians,
 * wrapped to [-pi, pi].
 */
extern "C" void ik_analytic_kernel(const ik_word_t pose[IK_DOF], int cfg,
                                   ik_word_t q[IK_DOF], int *status);

extern "C" void ik_dls_kernel(const ik_word_t pose[IK_DOF],
                              const ik_word_t q_seed[IK_DOF],
                              ik_word_t lambda, ik_word_t tol, int max_iter,
                              ik_word_t q[IK_DOF],
                              int *iters, ik_word_t *resid, int *status);

#endif /* IK_KERNELS_HPP */
