#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Forward kinematics, geometric Jacobian and pose-error helpers.
//
//! Shared by both IK kernels: the analytic solver needs rot03() to decouple the
//! wrist, the DLS solver needs fk_jacobian() and pose_error() once per
//! iteration.  Every routine has a fixed trip count.
//

namespace ikk {

    //! R = Rz(yaw) Ry(pitch) Rx(roll), ZYX convention. Angles in radians.
    void rpy_to_rot(ik_real_t roll, ik_real_t pitch, ik_real_t yaw,
                    ik_real_t R[3][3]);

    //! Inverse of rpy_to_rot().
    void rot_to_rpy(const ik_real_t R[3][3], ik_real_t rpy[3]);

    //! Full 6-DOF forward kinematics, base to tool frame.
    void fk(const ik_real_t q[IK_DOF], ik_real_t R[3][3], ik_real_t p[3]);

    //! Rotation of frame 3 (wrist base) from the first three DH steps; used
    //! by the analytic solver to decouple the wrist from the arm.
    void rot03(ik_real_t t1, ik_real_t t2, ik_real_t t3, ik_real_t R[3][3]);

    //! Tool pose and 6x6 geometric Jacobian (joint velocities -> end-effector
    //! linear/angular velocity in base frame). Rows ordered [vx vy vz wx wy wz].
    void fk_jacobian(const ik_real_t q[IK_DOF], ik_real_t R[3][3],
                     ik_real_t p[3], ik_real_t J[IK_MAT_MAX][IK_MAT_MAX]);

    //! 6-dim task-space error [position; orientation]. Orientation error is
    //! 0.5 * sum_k (c_k x d_k) over column pairs of Rd/Rc - avoids
    //! trigonometry, needed for the DLS inner loop.
    //! @param e output [ex, ey, ez, wx, wy, wz]
    void pose_error(const ik_real_t Rd[3][3], const ik_real_t pd[3],
                    const ik_real_t Rc[3][3], const ik_real_t pc[3],
                    ik_real_t e[IK_MAT_MAX]);

}  // namespace ikk

#endif  //! KINEMATICS_HPP
