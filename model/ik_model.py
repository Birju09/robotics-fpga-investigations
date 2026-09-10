#!/usr/bin/env python3
"""
Floating-point golden model for the 6-DOF IK kernels.

This module is the single source of truth for the kinematic conventions used by
the HLS implementation.  Every equation implemented in `hls/src/*.cpp` is first
validated here, and the C++ testbenches are checked against vectors emitted by
`gen_vectors()`.

Robot: 6R anthropomorphic arm with a spherical wrist (PUMA-style), standard
(distal) Denavit-Hartenberg parameters:

    i |  a_i  | alpha_i | d_i
    --+-------+---------+-----
    1 |  a1   |  +pi/2  | d1
    2 |  a2   |    0    | 0
    3 |  a3   |  +pi/2  | 0
    4 |  0    |  -pi/2  | d4
    5 |  0    |  +pi/2  | 0
    6 |  0    |    0    | d6

Pose convention: [x, y, z, roll, pitch, yaw] with R = Rz(yaw) Ry(pitch) Rx(roll).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Robot geometry (metres).  Must match hls/include/ik_config.hpp.
# ---------------------------------------------------------------------------
D1, A1, A2, A3, D4, D6 = 0.15, 0.05, 0.30, 0.02, 0.28, 0.08

DH_A = np.array([A1, A2, A3, 0.0, 0.0, 0.0])
DH_ALPHA = np.array([np.pi / 2, 0.0, np.pi / 2, -np.pi / 2, np.pi / 2, 0.0])
DH_D = np.array([D1, 0.0, 0.0, D4, 0.0, D6])

# Derived constants used by the closed-form solution.
L3 = float(np.hypot(A3, D4))        # effective forearm length
PHI = float(np.arctan2(D4, A3))     # forearm offset angle


# ---------------------------------------------------------------------------
# Forward kinematics
# ---------------------------------------------------------------------------
def dh_mat(a, alpha, d, th):
    """Standard (distal) DH homogeneous transform."""
    ct, st = np.cos(th), np.sin(th)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0.0, sa,       ca,      d],
        [0.0, 0.0,      0.0,     1.0],
    ])


def fk_all(q):
    """Return [T_0_0, T_0_1, ... T_0_6]."""
    T = np.eye(4)
    out = [T.copy()]
    for i in range(6):
        T = T @ dh_mat(DH_A[i], DH_ALPHA[i], DH_D[i], q[i])
        out.append(T.copy())
    return out


def fk(q):
    return fk_all(q)[-1]


def jacobian(q):
    """Geometric Jacobian in the base frame; rows [vx vy vz wx wy wz]."""
    Ts = fk_all(q)
    pe = Ts[6][:3, 3]
    J = np.zeros((6, 6))
    for i in range(6):
        z = Ts[i][:3, 2]
        p = Ts[i][:3, 3]
        J[:3, i] = np.cross(z, pe - p)
        J[3:, i] = z
    return J


# ---------------------------------------------------------------------------
# Pose <-> matrix helpers
# ---------------------------------------------------------------------------
def rpy_to_rot(roll, pitch, yaw):
    """R = Rz(yaw) Ry(pitch) Rx(roll)."""
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw), np.sin(yaw)
    return np.array([
        [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
        [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
        [-sp,     cp * sr,                cp * cr],
    ])


def rot_to_rpy(R):
    pitch = np.arctan2(-R[2, 0], np.hypot(R[0, 0], R[1, 0]))
    yaw = np.arctan2(R[1, 0], R[0, 0])
    roll = np.arctan2(R[2, 1], R[2, 2])
    return roll, pitch, yaw


def pose_to_T(pose):
    T = np.eye(4)
    T[:3, :3] = rpy_to_rot(pose[3], pose[4], pose[5])
    T[:3, 3] = pose[:3]
    return T


def T_to_pose(T):
    r, p, y = rot_to_rpy(T[:3, :3])
    return np.array([T[0, 3], T[1, 3], T[2, 3], r, p, y])


def wrap_pi(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


# ---------------------------------------------------------------------------
# Closed-form (analytic) inverse kinematics
# ---------------------------------------------------------------------------
def ik_analytic(T, shoulder=+1, elbow=+1, wrist=+1):
    """
    Closed-form IK for one of the eight branches.

    shoulder: +1 -> theta1 = atan2(pc_y, pc_x);  -1 -> that angle + pi
    elbow   : sign of sin(gamma)
    wrist   : +1 -> theta5 >= 0;  -1 -> the flipped wrist solution

    Returns (q, ok).  `ok` is False when the target is out of reach.
    """
    R = T[:3, :3]
    p = T[:3, 3]

    # 1. Wrist centre: frame-5 origin lies d6 back along the approach axis.
    pc = p - D6 * R[:, 2]

    # 2. Base rotation.
    t1 = np.arctan2(pc[1], pc[0])
    if shoulder < 0:
        t1 = wrap_pi(t1 + np.pi)

    # 3. Planar 2-link sub-problem in the frame-1 (x1, y1) plane.
    #    r = a2*c2 + L3*cos(beta),  s = a2*s2 + L3*sin(beta),  beta = t2 + t3 - PHI
    rad = np.hypot(pc[0], pc[1]) * (1.0 if shoulder > 0 else -1.0)
    r = rad - A1
    s = pc[2] - D1

    cos_g = (r * r + s * s - A2 * A2 - L3 * L3) / (2.0 * A2 * L3)
    if cos_g > 1.0 or cos_g < -1.0:
        return np.zeros(6), False          # unreachable
    sin_g = elbow * np.sqrt(max(0.0, 1.0 - cos_g * cos_g))

    gamma = np.arctan2(sin_g, cos_g)
    t3 = wrap_pi(gamma + PHI)
    t2 = wrap_pi(np.arctan2(s, r) - np.arctan2(L3 * sin_g, A2 + L3 * cos_g))

    # 4. Wrist: R_3_6 is a ZYZ Euler rotation in theta4/theta5/theta6.
    R03 = fk_all([t1, t2, t3, 0, 0, 0])[3][:3, :3]
    R36 = R03.T @ R

    sin_t5 = np.hypot(R36[0, 2], R36[1, 2])
    if sin_t5 < 1e-9:
        # Wrist singularity: only theta4 +/- theta6 is observable, so pin
        # theta4 at zero and put the whole rotation on theta6.
        #   t5 = 0  : R36 = Rz(t4 + t6)  -> R36[0,0] =  cos, R36[0,1] = -sin
        #   t5 = pi : R36 = Rz(t4) Ry(pi) Rz(t6)
        #             -> R36[0,0] = -cos(t4 - t6), R36[0,1] = -sin(t4 - t6)
        t4 = 0.0
        if R36[2, 2] > 0:
            t5 = 0.0
            t6 = np.arctan2(-R36[0, 1], R36[0, 0])
        else:
            t5 = np.pi
            t6 = np.arctan2(R36[0, 1], -R36[0, 0])
        t6 = wrap_pi(t6)
    elif wrist > 0:
        t4 = np.arctan2(R36[1, 2], R36[0, 2])
        t5 = np.arctan2(sin_t5, R36[2, 2])
        t6 = np.arctan2(R36[2, 1], -R36[2, 0])
    else:
        t4 = wrap_pi(np.arctan2(R36[1, 2], R36[0, 2]) + np.pi)
        t5 = np.arctan2(-sin_t5, R36[2, 2])
        t6 = wrap_pi(np.arctan2(-R36[2, 1], R36[2, 0]))

    return np.array([t1, t2, t3, t4, t5, t6]), True


def elbow_conditioning(T, shoulder=+1):
    """
    |sin(gamma)| for a target: how far the arm is from the elbow singularity.

    gamma = acos(cos_gamma), so the joint solution's sensitivity to numerical
    error in cos_gamma scales as 1/|sin(gamma)|.  As the arm approaches full
    extension or full fold this diverges and theta2/theta3 stop being well
    determined even though the pose is still reached accurately - the elbow-up
    and elbow-down branches are merging.  Test harnesses use this to decide
    whether joint-space agreement is a meaningful thing to assert.
    """
    R, p = T[:3, :3], T[:3, 3]
    pc = p - D6 * R[:, 2]
    rad = np.hypot(pc[0], pc[1]) * (1.0 if shoulder > 0 else -1.0)
    r = rad - A1
    s = pc[2] - D1
    cos_g = (r * r + s * s - A2 * A2 - L3 * L3) / (2.0 * A2 * L3)
    return float(np.sqrt(max(0.0, 1.0 - min(1.0, abs(cos_g)) ** 2)))


def ik_analytic_all(T):
    """All eight branches; returns list of (config_tuple, q)."""
    out = []
    for sh in (+1, -1):
        for el in (+1, -1):
            for wr in (+1, -1):
                q, ok = ik_analytic(T, sh, el, wr)
                if ok:
                    out.append(((sh, el, wr), q))
    return out


# ---------------------------------------------------------------------------
# Damped least squares inverse kinematics
# ---------------------------------------------------------------------------
def pose_error(T_des, T_cur):
    """6x1 error: position difference + Siciliano orientation error."""
    ep = T_des[:3, 3] - T_cur[:3, 3]
    Rd, Rc = T_des[:3, :3], T_cur[:3, :3]
    eo = 0.5 * (np.cross(Rc[:, 0], Rd[:, 0])
                + np.cross(Rc[:, 1], Rd[:, 1])
                + np.cross(Rc[:, 2], Rd[:, 2]))
    return np.concatenate([ep, eo])


def ik_dls(T_des, q0, lam=0.08, max_iter=64, tol=1e-5, step=1.0):
    """
    q <- q + J^T (J J^T + lam^2 I)^-1 e

    Returns (q, iters, err_norm, converged).  `iters` is the datum the
    real-time study cares about: it is data-dependent and unbounded in general,
    which is exactly what the fixed-latency analytic core is being compared to.
    """
    q = np.array(q0, dtype=float)
    I6 = np.eye(6)
    for k in range(1, max_iter + 1):
        T = fk(q)
        e = pose_error(T_des, T)
        err = float(np.linalg.norm(e))
        if err < tol:
            return q, k, err, True
        J = jacobian(q)
        A = J @ J.T + (lam * lam) * I6
        dq = J.T @ np.linalg.solve(A, e)
        q = q + step * dq
    T = fk(q)
    err = float(np.linalg.norm(pose_error(T_des, T)))
    return q, max_iter, err, err < tol


# ---------------------------------------------------------------------------
# Fixed-point helpers (Q16.16, matching ap_fixed<32,16> in the HLS build)
# ---------------------------------------------------------------------------
FRAC_BITS = 16
SCALE = 1 << FRAC_BITS


def to_fixed(x):
    """Float -> signed 32-bit Q16.16 integer (two's complement, as written to a register)."""
    v = int(np.floor(np.asarray(x, dtype=float) * SCALE + 0.5))
    v = max(-(1 << 31), min((1 << 31) - 1, v))
    return v & 0xFFFFFFFF


def from_fixed(v):
    v = int(v) & 0xFFFFFFFF
    if v >= (1 << 31):
        v -= (1 << 32)
    return v / SCALE


def quantize(x):
    return from_fixed(to_fixed(x))
