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

    //! Convert roll-pitch-yaw Euler angles to a 3x3 rotation matrix.
    //!
    //! Computes R = Rz(yaw) Ry(pitch) Rx(roll) using the ZYX Euler angle
    //! convention.
    //!
    //! @param roll  Roll angle (radians)
    //! @param pitch Pitch angle (radians)
    //! @param yaw   Yaw angle (radians)
    //! @param R     Output 3x3 rotation matrix
    void rpy_to_rot(ik_real_t roll, ik_real_t pitch, ik_real_t yaw,
                    ik_real_t R[3][3]);

    //! Convert a 3x3 rotation matrix to roll-pitch-yaw Euler angles.
    //!
    //! Inverse operation of rpy_to_rot().
    //!
    //! @param R   Input 3x3 rotation matrix
    //! @param rpy Output Euler angles as {roll, pitch, yaw} (radians)
    void rot_to_rpy(const ik_real_t R[3][3], ik_real_t rpy[3]);

    //! Compute full forward kinematics from base to tool frame.
    //!
    //! Applies the complete 6-DOF kinematic chain to compute the tool pose.
    //!
    //! @param q    Six joint angles (radians)
    //! @param R    Output 3x3 rotation matrix of tool frame
    //! @param p    Output tool position vector [x, y, z]
    void fk(const ik_real_t q[IK_DOF], ik_real_t R[3][3], ik_real_t p[3]);

    //! Compute rotation of frame 3 (wrist base) using the first three DH steps.
    //!
    //! Used by the analytic solver to decouple the wrist from the arm.
    //!
    //! @param t1 Joint 1 angle (radians)
    //! @param t2 Joint 2 angle (radians)
    //! @param t3 Joint 3 angle (radians)
    //! @param R  Output 3x3 rotation matrix
    void rot03(ik_real_t t1, ik_real_t t2, ik_real_t t3, ik_real_t R[3][3]);

    //! Compute forward kinematics and geometric Jacobian matrix.
    //!
    //! Calculates the tool pose and the 6x6 geometric Jacobian relating joint
    //! velocities to end-effector linear and angular velocities in the base
    //! frame. Jacobian rows are ordered [vx vy vz wx wy wz].
    //!
    //! @param q    Six joint angles (radians)
    //! @param R    Output 3x3 rotation matrix of tool frame
    //! @param p    Output tool position vector [x, y, z]
    //! @param J    Output 6x6 geometric Jacobian matrix
    void fk_jacobian(const ik_real_t q[IK_DOF], ik_real_t R[3][3],
                     ik_real_t p[3], ik_real_t J[IK_MAT_MAX][IK_MAT_MAX]);

    //! Compute 6-dimensional task-space error between desired and current
    //! poses.
    //!
    //! Returns the stack of position error (Cartesian difference) and
    //! orientation error (computed as 0.5 * sum_k (c_k x d_k) over the three
    //! column pairs of the two rotation matrices). This orientation formulation
    //! is computationally efficient and requires no trigonometry, making it
    //! suitable for the DLS inner loop.
    //!
    //! @param Rd   Desired rotation matrix (3x3)
    //! @param pd   Desired position vector [x, y, z]
    //! @param Rc   Current rotation matrix (3x3)
    //! @param pc   Current position vector [x, y, z]
    //! @param e    Output 6x1 error vector [ex, ey, ez, wx, wy, wz]
    void pose_error(const ik_real_t Rd[3][3], const ik_real_t pd[3],
                    const ik_real_t Rc[3][3], const ik_real_t pc[3],
                    ik_real_t e[IK_MAT_MAX]);

}  // namespace ikk

#endif  //! KINEMATICS_HPP
