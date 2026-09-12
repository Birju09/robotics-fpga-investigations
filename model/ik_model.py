#!/usr/bin/env python3
"""
Floating-point golden model for the 6-DOF IK kernels.

Single source of truth for the kinematic conventions used by the HLS
implementation. Every equation in `hls/src/*.cpp` is validated here first,
and the C++ testbenches are checked against vectors from `gen_vectors()`.

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
# Robot geometry (metres). hls/include/ik_geometry.hpp is generated from the
# same object by model/gen_geometry.py, so the C++ constants can't drift.
# ---------------------------------------------------------------------------
DH_A = ROBOT.a
DH_ALPHA = ROBOT.alpha
DH_D = ROBOT.d

D1, A1, A2, A3, D4, D6 = (ROBOT.d1, ROBOT.a1, ROBOT.a2, ROBOT.a3, ROBOT.d4,
                          ROBOT.d6)
#: Lateral shoulder offset (PUMA 560's 149.09 mm). Keeps the wrist centre off
#: the joint-1 axis where theta1 would be undefined.
D3 = ROBOT.d3

L3 = ROBOT.L3                       # effective forearm length
PHI = ROBOT.phi                     # forearm offset angle

# The closed form below assumes a particular pattern of DH zeros, not DH
# tables in general - fail here rather than as wrong joint angles.
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

    shoulder: sign of the sqrt in the theta1 offset term - lefty/righty, NOT
              simply pi apart here since the shoulder is offset
    elbow   : sign of sin(gamma)
    wrist   : +1 -> theta5 >= 0;  -1 -> flipped wrist solution

    Returns (q, ok); ok is False when unreachable, including inside the
    singular cylinder (step 2).
    """
    R = T[:3, :3]
    p = T[:3, 3]

    # 1. Wrist centre: frame-5 origin lies d6 back along the approach axis.
    pc = p - D6 * R[:, 2]

    # 2. Base rotation, with the lateral shoulder offset.
    #
    # d3 translates along z2 (parallel to z1), so joints 2/3 move the wrist
    # centre in a plane held a fixed distance d3 off the joint-1 axis - an
    # exact invariant independent of q2..q6:
    #
    #     pc_x sin(t1) - pc_y cos(t1) = d3
    #
    # solved by t1 = atan2(pc_y, pc_x) + atan2(d3, +-sqrt(rho^2 - d3^2)), the
    # two signs being the two shoulder configurations.
    #
    # The radicand also gives the reachability test a zero-offset arm doesn't
    # need: rho < |d3| puts the target inside an unreachable cylinder of
    # radius |d3| about the joint-1 axis. Without this check, the wrist
    # centre landing exactly on that axis makes t1 undefined and atan2(0,0)
    # quietly return zero instead of failing.
    rho_sq = pc[0] * pc[0] + pc[1] * pc[1]
    disc = rho_sq - D3 * D3
    if disc < 0.0:
        return np.zeros(6), False          # inside the singular cylinder
    t1 = wrap_pi(np.arctan2(pc[1], pc[0])
                 + np.arctan2(D3, shoulder * np.sqrt(disc)))

    # 3. Planar 2-link sub-problem in the frame-1 (x1, y1) plane.
    #    u = a2*c2 + L3*cos(beta),  w = a2*s2 + L3*sin(beta),  beta = t2 + t3 + PHI
    #
    # u carries the shoulder branch implicitly (lefty solution -> negative u),
    # unlike a zero-offset arm which must apply the branch sign by hand.
    u = pc[0] * np.cos(t1) + pc[1] * np.sin(t1) - A1
    w = pc[2] - D1

    cos_g = (u * u + w * w - A2 * A2 - L3 * L3) / (2.0 * A2 * L3)
    if cos_g > 1.0 or cos_g < -1.0:
        return np.zeros(6), False          # unreachable
    sin_g = elbow * np.sqrt(max(0.0, 1.0 - cos_g * cos_g))

    gamma = np.arctan2(sin_g, cos_g)
    # theta3 = gamma - PHI, not gamma + PHI: alpha3 = -pi/2 here flips the
    # sense of the elbow offset. Verified against fk(), see validate.py.
    t3 = wrap_pi(gamma - PHI)
    t2 = wrap_pi(np.arctan2(w, u) - np.arctan2(L3 * sin_g, A2 + L3 * cos_g))

    # 4. Wrist.  With alpha4 = +pi/2 and alpha5 = -pi/2,
    #
    #     R_3_6 = Rz(t4) Ry(-t5) Rz(t6)
    #
    # - ZYZ Euler with the middle angle NEGATED (opposite-twist arms give
    # Ry(+t5)), so t4/t6 below pick up a sign flip vs. textbook ZYZ. Wrong
    # sign fails silently (angles off by pi but still plausible), which is
    # why validate.py round-trips every branch through fk() instead of
    # checking joint values directly.
    R03 = fk_all([t1, t2, t3, 0, 0, 0])[3][:3, :3]
    R36 = R03.T @ R

    sin_t5 = np.hypot(R36[0, 2], R36[1, 2])
    if sin_t5 < 1e-9:
        # Wrist singularity: only theta4 +/- theta6 observable, so pin theta4
        # at zero and put the whole rotation on theta6.
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

    Sensitivity of the joint solution to numerical error in cos_gamma scales
    as 1/|sin(gamma)|; near full extension/fold this diverges and
    theta2/theta3 stop being well determined even though the pose is still
    reached accurately (elbow-up/down branches merging). Test harnesses use
    this to decide whether joint-space agreement is meaningful to assert.
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

    d3 replaces the zero-offset arm's singular LINE with a singular CYLINDER
    of radius |d3|, leaving a boundary layer just outside it where theta1 is
    finite but ill-conditioned:

        theta1 = atan2(pc_y, pc_x) + atan2(d3, root),  root = sqrt(rho^2 - d3^2)

    so d(theta1)/d(rho) = -d3 / (rho * root), diverging as root -> 0 - a
    half-LSB position error becomes milliradians of theta1, the same way it
    does for theta3 near the elbow singularity.

    Returns root/rho = cos(atan2(d3, root)) in [0, 1], zero on the cylinder -
    the analogue of elbow_conditioning()'s |sin gamma|.

    Measured on the 256-vector fixed-point set: every vector missing the
    5e-3 rad joint tolerance had root/rho < 0.14; every vector above 0.30
    agreed to 5.7e-4 or better. Testbenches assert joint agreement only above
    0.20, same policy as the elbow.
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


def joint_margin_rad(q_i, lo, hi):
    """
    Signed distance from q_i to the allowed interval [lo, hi], correctly
    handling the +-pi wrap for a single-turn joint (hi - lo < 2*pi).

    This is NOT the same as min(q_i - lo, hi - q_i) evaluated on q_i's
    principal value. That naive form picks up whichever of +180 deg/-180 deg
    floating-point noise happens to produce and reports wildly different
    margins for the same physical angle - found while sitting the pentagon
    trajectory (model/gen_vectors.py), where q4 sat exactly at the wrap
    boundary and the naive check alternately reported -10 deg and -70 deg
    for what is the same 10-degree overshoot. Positive means inside by this
    many radians; negative means outside by this many, via the shorter arc.

    A joint whose range already spans a full turn (hi - lo >= 2*pi, e.g. a
    multi-turn wrist roll) cannot be treated this way - there is no unique
    "shorter arc" - so it falls back to the plain linear check.
    """
    if hi - lo >= 2 * np.pi - 1e-9:
        return min(q_i - lo, hi - q_i)
    mid = 0.5 * (lo + hi)
    d = ((q_i - mid + np.pi) % (2 * np.pi)) - np.pi
    qq = mid + d
    if qq < lo:
        return qq - lo
    if qq > hi:
        return hi - qq
    return min(qq - lo, hi - qq)


def joint_limit_margins(q, qlim=None):
    """Per-joint joint_margin_rad(), against ROBOT.qlim unless given."""
    if qlim is None:
        qlim = ROBOT.qlim
    return np.array([joint_margin_rad(q[i], qlim[i, 0], qlim[i, 1])
                      for i in range(6)])


def _seg_dist(p0, p1, q0, q1):
    """Minimum distance between line segments p0-p1 and q0-q1 in R^3."""
    d1, d2, r = p1 - p0, q1 - q0, p0 - q0
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    if a <= 1e-12 and e <= 1e-12:
        return float(np.linalg.norm(p0 - q0))
    if a <= 1e-12:
        s, t = 0.0, np.clip(f / e, 0.0, 1.0)
    else:
        c = d1 @ r
        if e <= 1e-12:
            s, t = np.clip(-c / a, 0.0, 1.0), 0.0
        else:
            b = d1 @ d2
            denom = a * e - b * b
            s = np.clip((b * f - c * e) / denom, 0.0, 1.0) if abs(denom) > 1e-12 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)
    return float(np.linalg.norm((p0 + s * d1) - (q0 + t * d2)))


def self_clearance(q):
    """
    Minimum distance between any two non-adjacent PHYSICAL links of the
    zero-radius kinematic skeleton at configuration q.

    "Physical link" merges consecutive DH frames that coincide (a5 = d5 = 0
    puts frames 4 and 5 at the same point, the spherical wrist centre), so
    the reduced skeleton is 5 segments: pedestal, upper arm, the a3/d3
    shoulder-offset segment, forearm, tool offset. Pairs of segments that
    share an endpoint are adjacent and excluded; pairs separated by exactly
    one physical link are included and are naturally floored at that link's
    own length (its two endpoints are exactly that far apart, which is not a
    collision - just the joint between them), so a genuine close approach
    reads far below that floor rather than at it.

    This is NOT a collision check. No link cross-section, tool, or pedestal
    housing geometry is modeled anywhere in this project - only the
    kinematic point skeleton. A small value here means the physical arm
    WOULD interfere for any realistic link diameter; a large value does not
    prove clearance, since real links have thickness this ignores entirely.
    Found and used to retire a configuration where the pentagon trajectory
    (model/gen_vectors.py) swept the shoulder-offset link across the base
    pedestal's own axis (down to ~4 um) whenever theta2 crossed +-90 deg,
    which a1 = 0 (shoulder joint sits directly above the waist axis, no
    horizontal offset) makes a purely geometric certainty at that angle,
    for any position.
    """
    origins = [T[:3, 3] for T in fk_all(q)]
    nodes = [origins[0]]
    for p in origins[1:]:
        if np.linalg.norm(p - nodes[-1]) > 1e-9:
            nodes.append(p)
    m = len(nodes) - 1
    best = np.inf
    for a in range(m):
        for b in range(a + 2, m):
            best = min(best, _seg_dist(nodes[a], nodes[a + 1], nodes[b], nodes[b + 1]))
    return best


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


# Must match IK_DLS_STEP_MAX_DEFAULT / IK_DLS_STEP_HALVINGS in
# hls/include/ik_config.hpp, where the measurements behind these are recorded.
DLS_STEP_MAX = 1.5
DLS_STEP_HALVINGS = 8


def clamp_step(dq, step_max, halvings=DLS_STEP_HALVINGS):
    """
    Bound |dq|_inf by halving, exactly as iks::dls() does in hardware.

    Halving (not scaling by step_max/|dq|_inf) keeps a power of two exact in
    both the double model and the Q16.16 kernel, so both follow the same
    clamp trajectory. Effective radius is (step_max/2, step_max], not
    step_max - measured as costing nothing vs. an exact scale (sweep_dls.py).

    Scaling, not per-joint clipping, to preserve the step's direction as
    determined by the least-squares solve.
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

    Returns (q, iters, err_norm, converged). `iters` is data-dependent and
    unbounded in general - the point of comparison against the fixed-latency
    analytic core.

    `step_max` bounds |dq|_inf per iteration. Pass 0 for the pre-trust-region
    undamped-step behaviour - the baseline every ik_config.hpp measurement is
    quoted against, and what gen_vectors.py uses to keep the pose table fixed
    across solver changes.
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
    """Float -> signed 32-bit Q16.16 integer (two's complement)."""
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
