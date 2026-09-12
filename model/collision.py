#!/usr/bin/env python3
"""
Capsule collision model and clearance gradient - the golden reference for the
PL collision kernel.

This is the layer that turns `ik_model.self_clearance()` from a proxy into a
check. That function measures the zero-radius kinematic skeleton and says so
explicitly: "This is NOT a collision check. No link cross-section, tool, or
pedestal housing geometry is modeled anywhere in this project." What is added
here is exactly the missing radius, and nothing else - a capsule is a segment
plus a radius, and capsule-capsule distance is segment-segment distance minus
the two radii, exactly, with no approximation to argue about.

Why capsules rigid in a DH frame, not the skeleton polyline
-----------------------------------------------------------
The obvious model - one segment per pair of consecutive DH frame origins -
is NOT a rigid body, so it has no single point Jacobian. o_{i-1} is fixed in
frame i-1 and o_i is fixed in frame i, and the segment between them stretches
and shears as theta_i moves.

Decomposing each DH step into its two limbs fixes this, and the decomposition
is exact. T_{i-1,i} = Rz(theta_i) Tz(d_i) Tx(a_i) Rx(alpha_i), so in frame i:

    o_i        = (0, 0, 0)                               the joint centre
    P_i        = (-a_i, 0, 0)                            after Tz, before Tx
    o_{i-1}    = P_i - d_i * (0, sin alpha_i, cos alpha_i)

All three are CONSTANT in frame i - theta_i only rotates frame i itself. So
the whole polyline o_{i-1} -> P_i -> o_i is rigid in frame i, giving two
capsules ("d limb" along z_{i-1}, "a limb" along x_i), each with one point
Jacobian whose columns are the joints preceding that frame. That is the
property the PL kernel needs: one witness point, one gradient, and the
columns come straight out of the zax[]/org[] registers ikk::fk_jacobian()
already keeps.

It is also how the ecosystem models this. franka_description ships its coarse
self-collision geometry as sphere and cylinder sub-links rigidly attached to
a link frame (the `--with-sc` model), not as a skeleton polyline.

Radii are CHOSEN, not cited - see `LINK_RADIUS` below.
"""

from dataclasses import dataclass

import numpy as np

import ik_model as M
from robot import ROBOT

# ---------------------------------------------------------------------------
# Link radii, metres
# ---------------------------------------------------------------------------
# CHOSEN, NOT CITED. No published cross-section table for the PUMA 560 is
# used anywhere in this project, and none is known to the author. These are
# plausible half-thicknesses for an arm of this scale, ordered thick at the
# base and thin at the wrist, and they exist so the collision layer has a
# non-zero radius to work with at all.
#
# What that means for any result computed from them: a reported clearance is
# only as meaningful as these numbers. Treat a NEGATIVE clearance as a real
# finding (the skeleton is closer than any plausible link thickness allows)
# and a positive one as unproven. This is the same caveat self_clearance()
# carries, moved from "no radius at all" to "an assumed radius".
#
# One capsule per non-degenerate DH limb, indexed by DH step (1-based in the
# comment, 0-based in the array).
LINK_RADIUS = np.array([
    0.090,   # 1  pedestal / waist column
    0.075,   # 2  upper arm
    0.070,   # 3  shoulder offset + elbow
    0.055,   # 4  forearm
    0.045,   # 5  (degenerate on this arm: a5 = d5 = 0, no capsule emitted)
    0.040,   # 6  tool
])

# Adjacency is decided on shared SKELETON NODES, not on frame indices.
#
# A frame-index gap looks like it would work and does not. On this arm
# a5 = d5 = 0, so frame 5 contributes no capsule at all and the forearm
# (frame 4) and the tool (frame 6) are two frames apart yet physically share
# the wrist centre - they are adjacent, and checking them reports a permanent
# -0.095 m "collision" that is just the tool capsule lying inside the
# forearm's. Any arm with a degenerate DH step has this hole.
#
# Node ids are therefore assigned by walking the limb path and merging across
# every zero-length offset, which reproduces exactly the rule
# self_clearance() states: "pairs of segments that share an endpoint are
# adjacent and excluded; pairs separated by exactly one physical link are
# included".


@dataclass(frozen=True)
class Capsule:
    """
    A capsule rigidly attached to DH frame `frame` (an index into
    `ik_model.fk_all()`, so frame 0 is the base and frame 6 the tool).

    `pa`/`pb` are the segment endpoints in that frame's coordinates and are
    constant by construction (see the module docstring). `radius` is from
    LINK_RADIUS and is an assumption, not a measurement.

    `node_a`/`node_b` are ids in the deduplicated skeleton, used only for
    adjacency pruning.
    """
    name: str
    frame: int
    pa: np.ndarray
    pb: np.ndarray
    radius: float
    node_a: int
    node_b: int

    @property
    def length(self):
        return float(np.linalg.norm(self.pb - self.pa))


def capsules(rb=ROBOT, radius=None):
    """
    Build the capsule set from a RobotModel's DH table.

    Two limbs per DH step, either of which is dropped when degenerate
    (zero-length), since a zero-length capsule is a sphere and would only add
    a duplicate of its neighbour's endpoint. On the PUMA 560 that leaves six
    capsules: a5 = d5 = 0 emits nothing for step 5, and a1 = d2 = a4 = a6 = 0
    drop one limb each from steps 1, 2, 4 and 6.
    """
    rad = LINK_RADIUS if radius is None else np.asarray(radius, dtype=float)
    if len(rad) != rb.dof:
        raise ValueError(f"need {rb.dof} radii, got {len(rad)}")

    out = []
    node = 0                           # id of the node the walk is standing on
    for k in range(rb.dof):
        f = k + 1                      # fk_all() index of frame i
        a, d, al = float(rb.a[k]), float(rb.d[k]), float(rb.alpha[k])

        o_i = np.zeros(3)
        p_mid = np.array([-a, 0.0, 0.0])
        z_prev = np.array([0.0, np.sin(al), np.cos(al)])
        o_prev = p_mid - d * z_prev

        # d limb: o_{f-1} -> P_f, along z_{f-1}. A zero offset merges nodes
        # instead of emitting a zero-length capsule.
        if abs(d) > 1e-12:
            out.append(Capsule(f"link{f}_d", f, o_prev, p_mid, float(rad[k]),
                               node, node + 1))
            node += 1
        # a limb: P_f -> o_f, along x_f.
        if abs(a) > 1e-12:
            out.append(Capsule(f"link{f}_a", f, p_mid, o_i, float(rad[k]),
                               node, node + 1))
            node += 1
    return out


def _path_gap(caps, i, j):
    """Skeleton path length strictly between two capsules on the limb path."""
    lo = min(caps[i].node_b, caps[j].node_b)
    hi = max(caps[i].node_a, caps[j].node_a)
    return sum(c.length for c in caps if lo <= c.node_a and c.node_b <= hi)


def capsule_pairs(caps, explain=False):
    """
    Index pairs to check, after pruning pairs that cannot carry information.

    Two rules, and the second is one self_clearance() did not need:

    1. SHARED NODE. The two capsules meet at a joint, so their segment
       distance is zero there by construction and the "collision" reported is
       the joint, not an interference. Excluded - self_clearance()'s rule.

    2. UNSEPARABLE. The skeleton path between them is shorter than the sum of
       their radii, so the two capsule surfaces overlap even with the arm
       fully extended in a straight line: there is NO configuration in which
       the pair is clear. Excluded, because such a pair does not report a
       collision, it reports a modelling artefact - and being permanently the
       minimum, it would mask every real one.

       self_clearance() has no analogue of this because at zero radius the
       floor for a once-separated pair is the intervening limb's length,
       which is small but positive. Adding radii turns that floor negative
       whenever a link is shorter than the surrounding parts are thick. On
       the PUMA 560 exactly one pair trips it: link3_d/link4_d, separated by
       the 20.3 mm a3 elbow offset while the two radii sum to 125 mm. The
       physical arm has no gap there either - the elbow is one casting.

       Which means rule 2 is radius-dependent, so changing LINK_RADIUS can
       change the pair list. Pass explain=True to see what was dropped.

    Returns a list of (i, j), or (pairs, dropped) when explain=True.
    """
    pairs, dropped = [], []
    for i in range(len(caps)):
        for j in range(i + 1, len(caps)):
            nodes = {caps[i].node_a, caps[i].node_b,
                     caps[j].node_a, caps[j].node_b}
            if len(nodes) != 4:
                dropped.append((i, j, "shares a joint"))
                continue
            gap = _path_gap(caps, i, j)
            rsum = caps[i].radius + caps[j].radius
            if gap <= rsum:
                dropped.append((i, j, f"unseparable: path {gap:.4f} m <= "
                                      f"radii {rsum:.4f} m"))
                continue
            pairs.append((i, j))
    return (pairs, dropped) if explain else pairs


# ---------------------------------------------------------------------------
# Segment-segment distance with witness points
# ---------------------------------------------------------------------------
def seg_seg_witness(p0, p1, q0, q1, eps=1e-12):
    """
    Closest points between segments [p0,p1] and [q0,q1].

    Returns (dist, wa, wb, s, t) with wa = p0 + s*(p1-p0) and
    wb = q0 + t*(q1-q0), both clamped to [0,1].

    Same algorithm as `ik_model._seg_dist` (which returns only the distance);
    validate.py asserts the two agree, and cross-checks both against brute
    force sampling. The witness points are what the distance-only form cannot
    give and the clearance gradient cannot do without.

    The branch cascade is written as a cascade here, but it must be
    implemented in HLS as `select` muxes with every branch evaluated - see
    the trust-region and singular-flag precedents in iks::dls(), and
    docs/collision_aware_ik.md section 4.
    """
    p0 = np.asarray(p0, float)
    q0 = np.asarray(q0, float)
    d1 = np.asarray(p1, float) - p0
    d2 = np.asarray(q1, float) - q0
    r = p0 - q0

    a = float(d1 @ d1)
    e = float(d2 @ d2)
    f = float(d2 @ r)

    if a <= eps and e <= eps:
        s = t = 0.0
    elif a <= eps:
        s, t = 0.0, np.clip(f / e, 0.0, 1.0)
    else:
        c = float(d1 @ r)
        if e <= eps:
            s, t = np.clip(-c / a, 0.0, 1.0), 0.0
        else:
            b = float(d1 @ d2)
            denom = a * e - b * b
            s = np.clip((b * f - c * e) / denom, 0.0, 1.0) if abs(denom) > eps \
                else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)

    wa = p0 + s * d1
    wb = q0 + t * d2
    return float(np.linalg.norm(wa - wb)), wa, wb, float(s), float(t)


# ---------------------------------------------------------------------------
# Clearance
# ---------------------------------------------------------------------------
def capsule_endpoints(Ts, cap):
    """Capsule endpoints in the base frame at a configuration."""
    T = Ts[cap.frame]
    R, o = T[:3, :3], T[:3, 3]
    return R @ cap.pa + o, R @ cap.pb + o


def min_clearance(q, caps=None, pairs=None, rb=ROBOT):
    """
    Signed minimum clearance over all checked capsule pairs, and the witness
    data for the pair that achieved it.

    Clearance is segment distance minus both radii, so it goes NEGATIVE on
    interpenetration - the sign carries information and is deliberately not
    clamped at zero. Returns (d, info) where info is None only if there are
    no pairs to check.

    info holds ('pair', 'wa', 'wb', 'normal', 'seg_dist'). `normal` is the
    unit vector from wb toward wa, i.e. the direction capsule a must move to
    increase clearance. It is None when the witness points coincide exactly
    (a genuine crossing), where the direction is undefined.
    """
    caps = capsules(rb) if caps is None else caps
    pairs = capsule_pairs(caps) if pairs is None else pairs
    Ts = M.fk_all(q)

    best = np.inf
    info = None
    for (ia, ib) in pairs:
        a0, a1 = capsule_endpoints(Ts, caps[ia])
        b0, b1 = capsule_endpoints(Ts, caps[ib])
        sd, wa, wb, _s, _t = seg_seg_witness(a0, a1, b0, b1)
        d = sd - caps[ia].radius - caps[ib].radius
        if d < best:
            n = wa - wb
            nn = float(np.linalg.norm(n))
            best = d
            info = {
                "pair": (ia, ib),
                "wa": wa,
                "wb": wb,
                "normal": (n / nn) if nn > 1e-12 else None,
                "seg_dist": sd,
            }
    return (float(best) if info is not None else np.inf), info


def point_jacobian(Ts, p, frame):
    """
    3x6 translational Jacobian of a point rigidly attached to `frame`.

    J[:, m] = z_m x (p - o_m) for m < frame, zero otherwise: joint m+1
    rotates about z_m anchored at o_m, and a point attached to a frame at or
    before m is unaffected by it.

    These are the SAME z_m and o_m that ikk::fk_jacobian() already stores in
    zax[]/org[] with ARRAY_PARTITION complete, which is why the PL cost of a
    clearance gradient is six cross products on data already in registers -
    no extra FK pass and no extra CORDIC evaluation.
    """
    J = np.zeros((3, len(Ts) - 1))
    for m in range(frame):
        z = Ts[m][:3, 2]
        o = Ts[m][:3, 3]
        J[:, m] = np.cross(z, p - o)
    return J


def clearance_gradient(q, caps=None, pairs=None, rb=ROBOT):
    """
    Analytic gradient of the minimum clearance, d(d_min)/dq.

    Returns (d_min, grad, info). grad is zero when there is no pair to check
    or the witness points coincide.

        d = |wa - wb| - ra - rb,   n = (wa - wb)/|wa - wb|
        dd/dq = n^T (J_point(wa) - J_point(wb))

    Only the ACHIEVING pair contributes: d_min is a min over pairs, so its
    gradient is the gradient of whichever pair attains it. Holding the
    witness parameters (s,t) fixed at their optimum is correct by the same
    envelope argument - they are minimisers, so their own variation
    contributes nothing to first order. That is what makes the witness points
    the useful output of a distance query rather than just the distance.

    Not differentiable everywhere, and the non-smooth points are real rather
    than artefacts: where two pairs tie for the minimum, and where a witness
    point crosses a segment endpoint. Both are why the repulsion this feeds
    is gated by a threshold and clamped by the existing trust region rather
    than trusted as an exact descent direction.
    """
    caps = capsules(rb) if caps is None else caps
    pairs = capsule_pairs(caps) if pairs is None else pairs
    d, info = min_clearance(q, caps, pairs, rb)
    grad = np.zeros(rb.dof)
    if info is None or info["normal"] is None:
        return d, grad, info

    Ts = M.fk_all(q)
    ia, ib = info["pair"]
    Ja = point_jacobian(Ts, info["wa"], caps[ia].frame)
    Jb = point_jacobian(Ts, info["wb"], caps[ib].frame)
    grad = info["normal"] @ (Ja - Jb)
    return d, grad, info


# ---------------------------------------------------------------------------
# Nullspace projection (the consumer of the gradient)
# ---------------------------------------------------------------------------
def nullspace_project(q, z, lam=0.02, task_dim=M.TASK_FULL, rb=ROBOT):
    """
    (I - J_t^+ J_t) z with the DAMPED pseudoinverse, computed the way
    iks::dls() will: never forming the projector.

        v = J_t z ;  w = A^-1 v ;  Nz = z - J_t^T w ,   A = J_t J_t^T + lam^2 I

    so the only inverse is the task-dimensional one the solver already
    factors, and the second right-hand side reuses that factorisation - which
    is what spd::solve_n() exists for.

    At task_dim = 6 this is NOT zero, because damping makes J_t^+ an
    approximate inverse: the residual projector has eigenvalues
    lam^2/(sigma_i^2 + lam^2). That is the control case - it predicts
    near-zero avoidance authority away from singularities and growing
    authority as sigma -> 0.
    """
    rows = M.task_rows(task_dim)
    J = M.jacobian(q)[rows]
    A = J @ J.T + (lam * lam) * np.eye(task_dim)
    return np.asarray(z, float) - J.T @ np.linalg.solve(A, J @ np.asarray(z, float))


if __name__ == "__main__":
    caps = capsules()
    pairs, dropped = capsule_pairs(caps, explain=True)
    print(f"{ROBOT.name}: {len(caps)} capsules, {len(pairs)} checked pairs")
    print("   name        frame  node   length      radius")
    for c in caps:
        print(f"   {c.name:<11} {c.frame:^5}  {c.node_a}-{c.node_b}  "
              f"{c.length:8.5f}  {c.radius:8.3f}")
    print("   checked:")
    for i, j in pairs:
        print(f"     {caps[i].name:<11} {caps[j].name}")
    print("   pruned:")
    for i, j, why in dropped:
        print(f"     {caps[i].name:<11} {caps[j].name:<11} {why}")
    print()
    print("   RADII ARE CHOSEN, NOT CITED - see LINK_RADIUS. A negative")
    print("   clearance is a finding; a positive one is not a proof.")
    print()
    qz = np.zeros(ROBOT.dof)
    d, g, info = clearance_gradient(qz)
    print(f"   at q=0:  d_min = {d:+.5f} m  pair = "
          f"{caps[info['pair'][0]].name}/{caps[info['pair'][1]].name}")
    print(f"            grad  = {np.array2string(g, precision=4)}")
