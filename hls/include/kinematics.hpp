#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Forward kinematics, geometric Jacobian and pose-error helpers.
 *
 * Shared by both IK kernels: the analytic solver needs rot03() to decouple the
 * wrist, the DLS solver needs fk_jacobian() and pose_error() once per
 * iteration.  Every routine has a fixed trip count.
 */

namespace ikk {

/* R = Rz(yaw) Ry(pitch) Rx(roll) */
void rpy_to_rot(ik_real_t roll, ik_real_t pitch, ik_real_t yaw,
                ik_real_t R[3][3]);

/* Inverse of the above; rpy[] = {roll, pitch, yaw}. */
void rot_to_rpy(const ik_real_t R[3][3], ik_real_t rpy[3]);

/* Full forward kinematics: base -> tool. */
void fk(const ik_real_t q[IK_DOF], ik_real_t R[3][3], ik_real_t p[3]);

/* Rotation of frame 3 only - the first three DH steps. */
void rot03(ik_real_t t1, ik_real_t t2, ik_real_t t3, ik_real_t R[3][3]);

/* Forward kinematics plus the 6x6 geometric Jacobian in the base frame.
 * Rows are ordered [vx vy vz wx wy wz]. */
void fk_jacobian(const ik_real_t q[IK_DOF], ik_real_t R[3][3], ik_real_t p[3],
                 ik_real_t J[IK_MAT_MAX][IK_MAT_MAX]);

/* 6x1 task-space error: position difference stacked on the orientation error
 * 0.5 * sum_k (c_k x d_k) over the three column pairs of the two rotation
 * matrices.  Cheap, and needs no trig - which matters because it lands in the
 * DLS inner loop. */
void pose_error(const ik_real_t Rd[3][3], const ik_real_t pd[3],
                const ik_real_t Rc[3][3], const ik_real_t pc[3],
                ik_real_t e[IK_MAT_MAX]);

} /* namespace ikk */

#endif /* KINEMATICS_HPP */
