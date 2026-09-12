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

from robot import ROBOT, check_analytic_form

# ---------------------------------------------------------------------------
# Robot geometry (metres), taken from the definition layer rather than
# restated.  hls/include/ik_geometry.hpp is generated from the same object by
# model/gen_geometry.py, so the C++ constants cannot drift from these.
# ---------------------------------------------------------------------------
DH_A = ROBOT.a
DH_ALPHA = ROBOT.alpha
DH_D = ROBOT.d

D1, A1, A2, A3, D4, D6 = (ROBOT.d1, ROBOT.a1, ROBOT.a2, ROBOT.a3, ROBOT.d4,
                          ROBOT.d6)
#: Lateral shoulder offset, perpendicular to the arm plane.  Non-zero here
#: (PUMA 560's 149.09 mm), which is what keeps the wrist centre off the
#: joint-1 axis where theta1 would be undefined.
D3 = ROBOT.d3

# Derived constants used by the closed-form solution.
L3 = ROBOT.L3                       # effective forearm length
PHI = ROBOT.phi                     # forearm offset angle

# The closed form below is derived against a particular pattern of zeros, not
# against DH tables in general.  Fail here rather than as wrong joint angles.
check_analytic_form(ROBOT)


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

    shoulder: sign of the sqrt in the theta1 offset term; the two values are
              the lefty/righty pair, which for an offset shoulder are NOT
              simply pi apart the way they are on a zero-offset arm
    elbow   : sign of sin(gamma)
    wrist   : +1 -> theta5 >= 0;  -1 -> the flipped wrist solution

    Returns (q, ok).  `ok` is False when the target is out of reach, which now
    includes being inside the singular cylinder (see step 2).
    """
    R = T[:3, :3]
    p = T[:3, 3]

    # 1. Wrist centre: frame-5 origin lies d6 back along the approach axis.
    pc = p - D6 * R[:, 2]

    # 2. Base rotation, with the lateral shoulder offset.
    #
    # Because d3 translates along z2 (which is parallel to z1), joints 2 and 3
    # move the wrist centre in a plane held a fixed distance d3 off the
    # joint-1 axis.  That gives an exact invariant, independent of q2..q6:
    #
    #     pc_x sin(t1) - pc_y cos(t1) = d3
    #
    # whose solutions are t1 = atan2(pc_y, pc_x) + atan2(d3, +-sqrt(rho^2 -
    # d3^2)).  The two signs are the two shoulder configurations.
    #
    # The radicand is also the reachability test that the zero-offset version
    # of this solver did not need: rho < |d3| puts the target inside a
    # cylinder of radius |d3| about the joint-1 axis that the arm cannot
    # enter.  That cylinder is the whole point of the offset.  Without it the
    # wrist centre can sit exactly on the joint-1 axis, where t1 is undefined
    # and atan2(0, 0) quietly returns zero instead of failing.
    rho_sq = pc[0] * pc[0] + pc[1] * pc[1]
    disc = rho_sq - D3 * D3
    if disc < 0.0:
        return np.zeros(6), False          # inside the singular cylinder
    t1 = wrap_pi(np.arctan2(pc[1], pc[0])
                 + np.arctan2(D3, shoulder * np.sqrt(disc)))

    # 3. Planar 2-link sub-problem in the frame-1 (x1, y1) plane.
    #    u = a2*c2 + L3*cos(beta),  w = a2*s2 + L3*sin(beta),  beta = t2 + t3 + PHI
    #
    # u is the in-plane radial coordinate, and it carries the shoulder branch
    # implicitly: the lefty solution simply produces a negative u.  The old
    # zero-offset code had to apply the branch sign to the radius by hand.
    u = pc[0] * np.cos(t1) + pc[1] * np.sin(t1) - A1
    w = pc[2] - D1

    cos_g = (u * u + w * w - A2 * A2 - L3 * L3) / (2.0 * A2 * L3)
    if cos_g > 1.0 or cos_g < -1.0:
        return np.zeros(6), False          # unreachable
    sin_g = elbow * np.sqrt(max(0.0, 1.0 - cos_g * cos_g))

    gamma = np.arctan2(sin_g, cos_g)
    # theta3 = gamma - PHI, not gamma + PHI: alpha3 is -pi/2 on this arm, which
    # flips the sense of the elbow offset.  Verified against fk() rather than
    # assumed - see the round-trip assertions in validate.py.
    t3 = wrap_pi(gamma - PHI)
    t2 = wrap_pi(np.arctan2(w, u) - np.arctan2(L3 * sin_g, A2 + L3 * cos_g))

    # 4. Wrist.  With alpha4 = +pi/2 and alpha5 = -pi/2,
    #
    #     R_3_6 = Rz(t4) Ry(-t5) Rz(t6)
    #
    # - a ZYZ Euler set in which the middle angle enters NEGATED.  The
    # previous arm had the opposite twist signs and gave Ry(+t5), so the t4
    # and t6 extractions below each pick up a sign flip relative to the
    # textbook ZYZ formulae.  Getting this wrong does not fail loudly: it
    # returns joint angles that are wrong by pi in two joints while still
    # looking plausible, which is why validate.py round-trips every branch
    # through fk() rather than checking joint values.
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
        t4 = np.arctan2(-R36[1, 2], -R36[0, 2])
        t5 = np.arctan2(sin_t5, R36[2, 2])
        t6 = np.arctan2(-R36[2, 1], R36[2, 0])
    else:
        t4 = wrap_pi(np.arctan2(-R36[1, 2], -R36[0, 2]) + np.pi)
        t5 = np.arctan2(-sin_t5, R36[2, 2])
        t6 = wrap_pi(np.arctan2(R36[2, 1], -R36[2, 0]))

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
    disc = pc[0] * pc[0] + pc[1] * pc[1] - D3 * D3
    if disc < 0.0:
        return 0.0                  # inside the singular cylinder
    t1 = np.arctan2(pc[1], pc[0]) + np.arctan2(D3, shoulder * np.sqrt(disc))
    u = pc[0] * np.cos(t1) + pc[1] * np.sin(t1) - A1
    w = pc[2] - D1
    cos_g = (u * u + w * w - A2 * A2 - L3 * L3) / (2.0 * A2 * L3)
    return float(np.sqrt(max(0.0, 1.0 - min(1.0, abs(cos_g)) ** 2)))


def shoulder_conditioning(T):
    """
    sqrt(rho^2 - d3^2) / rho: how far a target is from the singular cylinder.

    The lateral offset d3 removes the singularity that a zero-offset arm has
    where the wrist centre meets the joint-1 axis.  It does not remove it for
    free.  It replaces a singular LINE with a singular CYLINDER of radius |d3|
    that the arm cannot enter, and leaves a boundary layer just outside it in
    which theta1 is finite but ill-conditioned:

        theta1 = atan2(pc_y, pc_x) + atan2(d3, root),  root = sqrt(rho^2 - d3^2)

    so d(theta1)/d(rho) = -d3 / (rho * root), which diverges as root -> 0.  A
    half-LSB of error in the wrist-centre position becomes milliradians of
    theta1, in exactly the way a half-LSB in the square root becomes
    milliradians of theta3 near the elbow singularity.

    Returned quantity is root/rho = cos(atan2(d3, root)), in [0, 1], zero on
    the cylinder - the direct analogue of elbow_conditioning()'s |sin gamma|.

    Measured on the 256-vector fixed-point set: every vector that missed the
    5e-3 rad joint tolerance had root/rho < 0.14, and every vector above 0.30
    agreed to 5.7e-4 or better.  The testbenches assert joint agreement only
    above 0.20 and report the rest, which is the same policy the elbow gets.
    """
    R, p = T[:3, :3], T[:3, 3]
    pc = p - D6 * R[:, 2]
    rho = float(np.hypot(pc[0], pc[1]))
    if rho < 1e-12:
        return 0.0
    disc = rho * rho - D3 * D3
    if disc <= 0.0:
        return 0.0                  # on or inside the cylinder
    return float(np.sqrt(disc) / rho)


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


# Trust region.  Must match IK_DLS_STEP_MAX_DEFAULT / IK_DLS_STEP_HALVINGS in
# hls/include/ik_config.hpp, which is where the measurements behind these
# values are recorded.
DLS_STEP_MAX = 1.5
DLS_STEP_HALVINGS = 8


def clamp_step(dq, step_max, halvings=DLS_STEP_HALVINGS):
    """
    Bound |dq|_inf by halving, exactly as iks::dls() does in hardware.

    Halving rather than scaling by step_max/|dq|_inf: a power of two is exact
    in both the double model and the Q16.16 kernel, so the two follow the same
    trajectory through the clamp.  The effective radius is therefore
    (step_max/2, step_max], not step_max - measured as costing nothing against
    an exact scale (see model/sweep_dls.py).

    Scaling, not per-joint clipping: clipping would change the step's
    direction, which is what the least-squares solve determined.
    """
    if step_max <= 0.0:
        return dq
    m = float(np.max(np.abs(dq)))
    sh = 0
    for _ in range(halvings):
        if m > step_max:
            m /= 2.0
            sh += 1
    return dq / float(1 << sh)


def ik_dls(T_des, q0, lam=0.08, max_iter=64, tol=1e-5, step_max=DLS_STEP_MAX):
    """
    q <- q + clamp( J^T (J J^T + lam^2 I)^-1 e )

    Returns (q, iters, err_norm, converged).  `iters` is the datum the
    real-time study cares about: it is data-dependent and unbounded in general,
    which is exactly what the fixed-latency analytic core is being compared to.

    `step_max` bounds |dq|_inf per iteration.  Pass 0 for the undamped-step
    behaviour this solver had before the trust region was added - that is the
    baseline every measurement in ik_config.hpp's trust region block is quoted
    against, and gen_vectors.py uses it deliberately to keep the pose table
    fixed across solver changes.
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
        dq = clamp_step(J.T @ np.linalg.solve(A, e), step_max)
        q = wrap_pi(q + dq)
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
