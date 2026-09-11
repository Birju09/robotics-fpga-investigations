#ifndef IK_KERNELS_HPP
#define IK_KERNELS_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! The two inverse-kinematics solvers under comparison.
//
//! analytic() - closed form.  Straight-line code: no loop whose trip count
//! depends on the data, so C-synthesis reports one latency and
//! that latency is the worst case.
//
//! dls()      - damped least squares.  Iterates until the task-space residual
//! falls below `tol` or the cap is hit, so its latency is a
//! function of the pose, the seed and the conditioning of the
//! Jacobian.  It reports the iteration count it actually used;
//! that number, and its distribution over a workload, is the
//! measurement this project exists to take.
//

namespace iks {

    //! Solve inverse kinematics using analytic closed-form solution.
    //!
    //! Computes joint angles for a desired end-effector pose using a
    //! closed-form algebraic solution. Execution time is constant regardless of
    //! input data, making this suitable for real-time applications with
    //! deterministic latency requirements.
    //!
    //! @param Rd   Desired end-effector orientation (3x3 rotation matrix)
    //! @param pd   Desired end-effector position [x, y, z] in Cartesian
    //! coordinates
    //! @param cfg  Configuration selector using IK_CFG_* bit flags to choose
    //! among
    //!             the eight possible IK branches
    //! @param q    Output joint angle solution (6 elements, in radians)
    //!
    //! @return Status code:
    //!   - IK_OK: Solution found
    //!   - IK_ERR_UNREACH: Target pose is not reachable by the robot
    int analytic(const ik_real_t Rd[3][3], const ik_real_t pd[3], int cfg,
                 ik_real_t q[IK_DOF]);

    //! Solve inverse kinematics using damped least-squares iterative method.
    //!
    //! Iteratively refines a joint angle seed using the damped least-squares
    //! (DLS) algorithm to minimize task-space pose error. Execution time is
    //! data-dependent, varying with the number of iterations required for
    //! convergence.
    //!
    //! @param Rd       Desired end-effector orientation (3x3 rotation matrix)
    //! @param pd       Desired end-effector position [x, y, z]
    //! @param q_seed   Initial joint angle guess (6 elements, in radians)
    //! @param lambda   Damping coefficient for numerical stability
    //! @param tol      Convergence tolerance for task-space residual norm
    //! @param max_iter Maximum iteration count before forced termination
    //! @param q        Output joint angle solution (6 elements, in radians)
    //! @param iters    Output pointer to receive actual iteration count
    //! executed
    //! @param resid    Output pointer to receive final task-space residual
    //! (2-norm, Q16.16)
    //!
    //! @return Status code:
    //!   - IK_OK: Converged to acceptable solution
    //!   - IK_ERR_NO_CONV: Reached max iterations without convergence
    //!   - IK_ERR_SINGULAR: Jacobian became singular during iteration
    int dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
            const ik_real_t q_seed[IK_DOF], ik_real_t lambda, ik_real_t tol,
            int max_iter, ik_real_t q[IK_DOF], int* iters, ik_real_t* resid);

}  // namespace iks

//! @defgroup IKKernels Inverse Kinematics IP Kernels
//! Standalone IP cores for inverse kinematics computation.
//! @{

//! Analytic inverse kinematics IP kernel with AXI4-Lite interface.
//!
//! Standalone hardware implementation of the closed-form IK solver.
//! All parameters are Q16.16 fixed-point words exposed through AXI4-Lite
//! register banks.
//!
//! @param pose Desired end-effector pose [x, y, z, roll, pitch, yaw] in Q16.16,
//!             where R = Rz(yaw) Ry(pitch) Rx(roll)
//! @param cfg  Configuration selector for IK branch selection
//! @param q    Output joint angles (6 elements) in radians, wrapped to [-π, π]
//! @param status Output status: 0 on success, non-zero on error
extern "C" void ik_analytic_kernel(const ik_word_t pose[IK_DOF], int cfg,
                                   ik_word_t q[IK_DOF], int* status);

//! Damped least-squares inverse kinematics IP kernel with AXI4-Lite interface.
//!
//! Standalone hardware implementation of the iterative DLS IK solver.
//! All parameters are Q16.16 fixed-point words exposed through AXI4-Lite
//! register banks.
//!
//! @param pose     Desired end-effector pose [x, y, z, roll, pitch, yaw] in
//! Q16.16,
//!                 where R = Rz(yaw) Ry(pitch) Rx(roll)
//! @param q_seed   Initial joint angle seed (6 elements) in radians
//! @param lambda    Damping coefficient for numerical stability (Q16.16)
//! @param tol       Convergence tolerance for residual norm (Q16.16)
//! @param max_iter  Maximum iteration count
//! @param q        Output joint angles (6 elements) in radians
//! @param iters    Output pointer to iteration count executed
//! @param resid    Output pointer to final task-space residual norm (Q16.16)
//! @param status   Output status: 0 on convergence, non-zero on error
extern "C" void ik_dls_kernel(const ik_word_t pose[IK_DOF],
                              const ik_word_t q_seed[IK_DOF], ik_word_t lambda,
                              ik_word_t tol, int max_iter, ik_word_t q[IK_DOF],
                              int* iters, ik_word_t* resid, int* status);
//! @}

#endif  //! IK_KERNELS_HPP
