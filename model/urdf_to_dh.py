#!/usr/bin/env python3
"""
Derive a standard-DH RobotModel from a URDF.

URDF describes a robot as a tree of links joined by transforms: each joint
carries a fixed origin (xyz + rpy, parent link frame to joint frame) and a
rotation axis expressed in that joint frame.  It places frames wherever the
CAD happened to put them.  DH places frames where the geometry says they must
go - z on the joint axis, x along the common normal between consecutive axes -
and buys, in exchange for that rigidity, four parameters per joint instead of
six.  This module performs that re-framing.

    urdf -> joint chain -> axis lines at q = 0 -> common-normal construction
         -> (a, alpha, d, theta_offset) + base and tool transforms

Two things are worth saying plainly before relying on the output.

**DH parameters are not unique.**  Axis directions may be flipped, and where
consecutive axes are parallel the common normal is undetermined and a
convention has to be imposed.  Two correct converters can therefore disagree
on the numbers while describing the same arm.  Comparing tables is not a valid
test; comparing forward kinematics is.  `verify()` does the latter and
`from_urdf()` calls it, so a conversion either reproduces the URDF's own FK to
tolerance or raises.

**A general URDF is not expressible in DH at all.**  DH describes a serial
chain of revolute/prismatic joints; a tree, a closed loop, or a joint whose
axis does not admit a common normal with its neighbour has no DH form.  Those
cases raise `UrdfConversionError` rather than returning something plausible.

Prismatic joints are parsed and rejected: the rest of this project - the
analytic solver, the fixed-point format, the Jacobian - assumes six revolute
joints throughout, so accepting them here would produce a model nothing
downstream could consume.

Self-test:  python3 model/urdf_to_dh.py
"""

import os
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from robot import RobotModel  # noqa: E402

#: Two axis directions closer than this in the cross-product norm are treated
#: as parallel, which is a separate construction case rather than a numerical
#: nuisance: the common normal between parallel axes is undetermined and the
#: convention below has to pick one.
PARALLEL_EPS = 1e-9

#: Below this, a common normal has zero length and the axes are taken to
#: intersect (a = 0), so x is chosen from the cross product alone.
LENGTH_EPS = 1e-9


class UrdfConversionError(ValueError):
    """The URDF does not describe something expressible as a DH chain."""


# ---------------------------------------------------------------------------
# URDF parsing
# ---------------------------------------------------------------------------
@dataclass
class UrdfJoint:
    name: str
    jtype: str
    parent: str
    child: str
    origin: np.ndarray      # 4x4, parent link frame -> joint frame
    axis: np.ndarray        # unit, in the joint frame
    lower: float
    upper: float


def rpy_to_rot(r, p, y):
    """URDF's fixed-axis roll-pitch-yaw, R = Rz(y) Ry(p) Rx(r)."""
    cr, sr, cp, sp, cy, sy = (np.cos(r), np.sin(r), np.cos(p), np.sin(p),
                              np.cos(y), np.sin(y))
    return np.array([
        [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
        [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
        [-sp, cp * sr, cp * cr],
    ])


def _origin_of(elem):
    T = np.eye(4)
    o = elem.find("origin")
    if o is None:
        return T
    xyz = [float(v) for v in o.get("xyz", "0 0 0").split()]
    rpy = [float(v) for v in o.get("rpy", "0 0 0").split()]
    T[:3, :3] = rpy_to_rot(*rpy)
    T[:3, 3] = xyz
    return T


def parse_urdf(path):
    """Return {name: UrdfJoint} and the set of link names."""
    root = ET.parse(path).getroot()
    if root.tag != "robot":
        raise UrdfConversionError(f"root element is <{root.tag}>, not <robot>")

    joints = {}
    for j in root.findall("joint"):
        jtype = j.get("type")
        parent = j.find("parent")
        child = j.find("child")
        if parent is None or child is None:
            raise UrdfConversionError(
                f"joint '{j.get('name')}' lacks a parent or child")
        ax = j.find("axis")
        axis = np.array([float(v) for v in ax.get("xyz").split()]) \
            if ax is not None else np.array([1.0, 0.0, 0.0])
        n = np.linalg.norm(axis)
        if jtype != "fixed" and n < LENGTH_EPS:
            raise UrdfConversionError(
                f"joint '{j.get('name')}' has a zero-length axis")
        lim = j.find("limit")
        lo = float(lim.get("lower", -np.pi)) if lim is not None else -np.pi
        hi = float(lim.get("upper", np.pi)) if lim is not None else np.pi
        joints[j.get("name")] = UrdfJoint(
            name=j.get("name"), jtype=jtype,
            parent=parent.get("link"), child=child.get("link"),
            origin=_origin_of(j), axis=axis / max(n, LENGTH_EPS),
            lower=lo, upper=hi)

    links = {l.get("name") for l in root.findall("link")}
    return joints, links, root.get("name", os.path.basename(path))


def joint_chain(joints, base=None, tip=None):
    """
    The ordered joint list from `base` to `tip`.

    Defaults walk the unique serial path: base is the link that is never a
    child, tip the link that is never a parent.  A tree has more than one such
    tip, and there is no DH form for a tree, so that raises here rather than
    silently converting one arbitrary branch.
    """
    children = {j.child for j in joints.values()}
    parents = {j.parent for j in joints.values()}
    if base is None:
        roots = sorted(parents - children)
        if len(roots) != 1:
            raise UrdfConversionError(
                f"expected exactly one base link, found {roots or 'none'}")
        base = roots[0]
    if tip is None:
        tips = sorted(children - parents)
        if len(tips) != 1:
            raise UrdfConversionError(
                f"expected one tip link (a serial chain), found {tips}; "
                "a branching tree has no DH form - pass tip= to choose one")
        tip = tips[0]

    by_parent = {}
    for j in joints.values():
        by_parent.setdefault(j.parent, []).append(j)

    chain, link, guard = [], base, 0
    while link != tip:
        nxt = by_parent.get(link)
        if not nxt:
            raise UrdfConversionError(f"no path from '{base}' to '{tip}'")
        if len(nxt) > 1:
            raise UrdfConversionError(
                f"link '{link}' has {len(nxt)} child joints; the chain from "
                f"'{base}' to '{tip}' is not serial")
        chain.append(nxt[0])
        link = nxt[0].child
        guard += 1
        if guard > len(joints) + 1:
            raise UrdfConversionError("cycle detected in the joint tree")
    return chain


def axis_lines(chain):
    """
    Each movable joint's axis as a line (point, unit direction) at q = 0,
    in the base link frame, with the fixed joints folded into the transforms
    on either side.

    Returns (lines, T_tip): T_tip is the base-to-tip transform at q = 0, which
    the tool frame is measured against.
    """
    T = np.eye(4)
    lines = []
    for j in chain:
        T = T @ j.origin
        if j.jtype == "fixed":
            continue
        if j.jtype in ("prismatic", "planar", "floating"):
            raise UrdfConversionError(
                f"joint '{j.name}' is '{j.jtype}'; this project's solvers, "
                "fixed-point format and Jacobian all assume revolute joints "
                "throughout, so a prismatic chain has no usable output here")
        if j.jtype not in ("revolute", "continuous"):
            raise UrdfConversionError(
                f"joint '{j.name}' has unknown type '{j.jtype}'")
        lines.append((T[:3, 3].copy(), (T[:3, :3] @ j.axis).copy()))
    return lines, T


# ---------------------------------------------------------------------------
# Common-normal construction
# ---------------------------------------------------------------------------
def _unit_perp(z):
    """Any unit vector perpendicular to z, chosen deterministically."""
    seed = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 \
        else np.array([0.0, 1.0, 0.0])
    x = np.cross(z, seed)
    return x / np.linalg.norm(x)


def common_normal(p0, z0, p1, z1):
    """
    The common normal from line (p0, z0) to line (p1, z1).

    Returns (a, x, foot0), where `a` is its length, `x` its unit direction
    pointing from the first line to the second, and `foot0` the point where it
    meets the first line.

    Three cases, and they are cases rather than one formula with a guard:

    skew         the generic one; x is the normalised cross product and the
                 feet come from solving the two-line closest-approach system.
    intersecting the cross product still fixes x, but a = 0 and the foot is
                 the intersection point.
    parallel     the cross product vanishes and the common normal is NOT
                 unique - every perpendicular between the lines has the same
                 length.  The universal convention, and the one taken here, is
                 to choose the normal through p0, which is what makes d_i come
                 out as zero rather than arbitrary for a parallel pair.
    """
    cz = np.cross(z0, z1)
    n = np.linalg.norm(cz)

    if n < PARALLEL_EPS:
        w = p1 - p0
        perp = w - np.dot(w, z0) * z0
        a = float(np.linalg.norm(perp))
        x = perp / a if a > LENGTH_EPS else _unit_perp(z0)
        return (a if a > LENGTH_EPS else 0.0), x, p0.copy()

    x = cz / n
    a = float(np.dot(p1 - p0, x))
    if a < 0.0:                     # orient from line 0 towards line 1
        x, a = -x, -a
    # Closest approach: p0 + t0 z0 + a x = p1 + t1 z1
    A = np.column_stack([z0, -z1])
    t, *_ = np.linalg.lstsq(A, (p1 - p0) - a * x, rcond=None)
    return (a if a > LENGTH_EPS else 0.0), x, p0 + t[0] * z0


def dh_from_axes(lines, T_tip):
    """
    Standard-DH parameters for a chain of joint axes.

    Frame i-1 has z on joint i's axis, so the construction walks pairs of
    consecutive axes.  The last frame has no successor axis to build a common
    normal against; its z is carried over from the previous joint and its
    origin placed at the tip, which is the convention that puts the tool
    offset in d_n where the rest of this project expects it.

    Returns (a, alpha, d, theta_offset, T_base, T_tool).
    """
    n = len(lines)
    if n < 2:
        raise UrdfConversionError(f"{n} revolute joints; need at least 2")

    z = [ln[1] for ln in lines]
    p = [ln[0] for ln in lines]

    # x_i for i = 1..n-1 from consecutive axis pairs.
    xs, origins = [None] * (n + 1), [None] * (n + 1)
    a = np.zeros(n)
    for i in range(1, n):
        a_i, x_i, foot0 = common_normal(p[i - 1], z[i - 1], p[i], z[i])
        a[i - 1] = a_i
        xs[i] = x_i
        origins[i] = foot0 + a_i * x_i        # on axis i, i.e. z_i

    # Frame 0: z0 is joint 1's axis; x0 is free.  Aligning it with x1 makes
    # theta_offset_1 zero, which is what a hand-written table would do, and
    # putting O0 at the foot of x1 on z0 makes d1 the base offset.
    xs[0] = xs[1].copy()
    _, _, origins[0] = common_normal(p[0], z[0], p[1], z[1])

    # Last frame: there is no next axis to build a common normal against, so
    # z_n carries over and the origin goes to the tip.  x_n must then point
    # from z_{n-1} TOWARDS the tip - specifically along the component of the
    # tool offset perpendicular to z_{n-1} - because that perpendicular part
    # is a_n, and a_n is expressed in x_n.  Copying x_{n-1} instead is correct
    # only when the tool sits on the axis; the parallel-axis case in the
    # self-test has a tool offset perpendicular to it and catches the error as
    # a position discrepancy of exactly a_n.
    z.append(z[n - 1].copy())
    origins[n] = T_tip[:3, 3].copy()
    off = origins[n] - origins[n - 1]
    perp = off - np.dot(off, z[n - 1]) * z[n - 1]
    ln = float(np.linalg.norm(perp))
    xs[n] = perp / ln if ln > LENGTH_EPS else xs[n - 1].copy()

    alpha = np.zeros(n)
    d = np.zeros(n)
    theta = np.zeros(n)
    for i in range(1, n + 1):
        z_prev, z_cur, x_prev, x_cur = z[i - 1], z[i], xs[i - 1], xs[i]
        # alpha: z_prev -> z_cur about x_cur
        alpha[i - 1] = np.arctan2(np.dot(np.cross(z_prev, z_cur), x_cur),
                                  np.dot(z_prev, z_cur))
        # d: along z_prev, from O_prev to the foot of x_cur
        d[i - 1] = float(np.dot(origins[i] - origins[i - 1], z_prev))
        # theta offset: x_prev -> x_cur about z_prev
        theta[i - 1] = np.arctan2(np.dot(np.cross(x_prev, x_cur), z_prev),
                                  np.dot(x_prev, x_cur))
        if i < n:
            a[i - 1] = float(np.dot(origins[i] - origins[i - 1]
                                    - d[i - 1] * z_prev, x_cur))

    # a_n: the tool offset perpendicular to z_{n-1}, if any.
    a[n - 1] = float(np.dot(origins[n] - origins[n - 1]
                            - d[n - 1] * z[n - 1], xs[n]))

    T_base = np.eye(4)
    T_base[:3, 0], T_base[:3, 1] = xs[0], np.cross(z[0], xs[0])
    T_base[:3, 2], T_base[:3, 3] = z[0], origins[0]

    T_tool = np.eye(4)
    Rn = np.column_stack([xs[n], np.cross(z[n], xs[n]), z[n]])
    T_tool[:3, :3] = Rn.T @ T_tip[:3, :3]
    return a, alpha, d, theta, T_base, T_tool


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def _dh_mat(a, alpha, d, th):
    ct, st, ca, sa = np.cos(th), np.sin(th), np.cos(alpha), np.sin(alpha)
    return np.array([[ct, -st * ca, st * sa, a * ct],
                     [st, ct * ca, -ct * sa, a * st],
                     [0.0, sa, ca, d],
                     [0.0, 0.0, 0.0, 1.0]])


def dh_fk(rb, q, T_base=None, T_tool=None):
    T = np.eye(4) if T_base is None else T_base.copy()
    for i in range(rb.dof):
        T = T @ _dh_mat(rb.a[i], rb.alpha[i], rb.d[i],
                        q[i] + rb.theta_offset[i])
    return T if T_tool is None else T @ T_tool


def urdf_fk(chain, q):
    T = np.eye(4)
    k = 0
    for j in chain:
        T = T @ j.origin
        if j.jtype == "fixed":
            continue
        c, s = np.cos(q[k]), np.sin(q[k])
        ax = j.axis
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]],
                      [-ax[1], ax[0], 0]])
        R = np.eye(3) + s * K + (1 - c) * (K @ K)      # Rodrigues
        Tj = np.eye(4)
        Tj[:3, :3] = R
        T = T @ Tj
        k += 1
    return T


def verify(rb, chain, T_base, T_tool, n=200, tol=1e-9, seed=0):
    """
    Compare DH forward kinematics against the URDF's own, over random
    configurations.

    This is the only meaningful correctness test for a DH conversion.  DH
    parameters are not unique - axis flips and the parallel-axis convention
    both admit different tables for the same arm - so asserting on the numbers
    would fail correct conversions.  Asserting that the two chains put the tool
    in the same place cannot.

    Returns the worst (position, rotation) deviation; raises past `tol`.
    """
    rng = np.random.default_rng(seed)
    wp = wr = 0.0
    for _ in range(n):
        q = rng.uniform(-np.pi, np.pi, size=rb.dof)
        Ta, Tb = dh_fk(rb, q, T_base, T_tool), urdf_fk(chain, q)
        wp = max(wp, float(np.abs(Ta[:3, 3] - Tb[:3, 3]).max()))
        wr = max(wr, float(np.abs(Ta[:3, :3] - Tb[:3, :3]).max()))
    if max(wp, wr) > tol:
        raise UrdfConversionError(
            f"conversion does not reproduce the URDF's kinematics: worst "
            f"position {wp:.3e} m, rotation {wr:.3e} (tol {tol:.1e}). The "
            "chain is probably not expressible in DH - check for a joint "
            "whose axis has no common normal with its neighbour.")
    return wp, wr


# ---------------------------------------------------------------------------
def from_urdf(path, base=None, tip=None, name=None, verify_tol=1e-9):
    """
    Build a RobotModel from a URDF file.

    Verifies the result against the URDF's own forward kinematics before
    returning it, so the caller receives either a conversion that reproduces
    the source or an exception - never an unchecked table.

    The returned model carries `T_base` and `T_tool` in `provenance`: a DH
    chain generally cannot absorb the URDF's base and tool frames, and
    silently dropping them would move the robot.
    """
    joints, _links, robot_name = parse_urdf(path)
    chain = joint_chain(joints, base, tip)
    lines, T_tip = axis_lines(chain)
    a, alpha, d, theta, T_base, T_tool = dh_from_axes(lines, T_tip)

    movable = [j for j in chain if j.jtype != "fixed"]
    qlim = np.array([[j.lower, j.upper] for j in movable])

    rb = RobotModel(
        name=name or robot_name,
        a=a, alpha=alpha, d=d, theta_offset=theta, qlim=qlim,
        provenance={
            "source": f"converted from {os.path.basename(path)} by "
                      "model/urdf_to_dh.py",
            "chain": " -> ".join(j.name for j in chain),
            "T_base": np.array2string(T_base, precision=6),
            "T_tool": np.array2string(T_tool, precision=6),
        })
    wp, wr = verify(rb, chain, T_base, T_tool, tol=verify_tol)
    rb.provenance["verified"] = (
        f"FK matches the URDF over 200 random configurations to "
        f"{wp:.2e} m / {wr:.2e}")
    return rb, T_base, T_tool


# ---------------------------------------------------------------------------
# URDF emission - used to build the self-test corpus
# ---------------------------------------------------------------------------
def to_urdf(rb, path, T_base=None):
    """
    Write a URDF describing the same kinematics as `rb`.

    Used to build the round-trip corpus below.  The construction relies on
    T_i-1_i(q) = Rz(q + off_i) * C_i with C_i = Tz(d) Tx(a) Rx(alpha): C_i is
    constant, so it becomes the ORIGIN of joint i+1, and each joint itself is a
    bare rotation about z.  C_n goes on a final fixed joint to the tip.
    """
    def xyz_rpy(T):
        R = T[:3, :3]
        p = float(np.arctan2(-R[2, 0], np.hypot(R[0, 0], R[1, 0])))
        y = float(np.arctan2(R[1, 0], R[0, 0]))
        r = float(np.arctan2(R[2, 1], R[2, 2]))
        return (" ".join(f"{v:.17g}" for v in T[:3, 3]),
                f"{r:.17g} {p:.17g} {y:.17g}")

    C = [_dh_mat(rb.a[i], rb.alpha[i], rb.d[i], 0.0) for i in range(rb.dof)]
    off = [np.eye(4)] + C[:-1]

    out = [f'<robot name="{rb.name}">', '  <link name="base_link"/>']
    for i in range(rb.dof):
        pre = off[i].copy()
        if i == 0 and T_base is not None:
            pre = T_base @ pre
        # theta_offset folds into the joint origin as a rotation about z.
        if abs(rb.theta_offset[i]) > 0:
            t = rb.theta_offset[i]
            Rz = np.eye(4)
            Rz[:2, :2] = [[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]]
            pre = pre @ Rz
        xyz, rpy = xyz_rpy(pre)
        parent = "base_link" if i == 0 else f"link{i}"
        out += [f'  <link name="link{i + 1}"/>',
                f'  <joint name="joint{i + 1}" type="revolute">',
                f'    <parent link="{parent}"/>',
                f'    <child link="link{i + 1}"/>',
                f'    <origin xyz="{xyz}" rpy="{rpy}"/>',
                '    <axis xyz="0 0 1"/>',
                f'    <limit lower="{rb.qlim[i, 0]:.17g}" '
                f'upper="{rb.qlim[i, 1]:.17g}" effort="100" velocity="3"/>',
                '  </joint>']
    xyz, rpy = xyz_rpy(C[-1])
    out += ['  <link name="tool"/>',
            '  <joint name="tool_fixed" type="fixed">',
            f'    <parent link="link{rb.dof}"/>',
            '    <child link="tool"/>',
            f'    <origin xyz="{xyz}" rpy="{rpy}"/>',
            '  </joint>', '</robot>', '']
    with open(path, "w") as f:
        f.write("\n".join(out))
    return path


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def _self_test():
    """
    Round-trip corpus: DH table -> URDF -> DH table, checked on FK.

    Building the corpus from known tables rather than collecting real URDFs
    keeps the test hermetic and, more usefully, lets it cover the awkward
    cases on purpose - parallel axes, intersecting axes, zero link lengths -
    which is where a common-normal construction actually breaks.
    """
    import tempfile
    from robot import PUMA560

    cases = [
        ("PUMA 560 (skew, offsets, spherical wrist)", PUMA560),
        ("parallel axes (SCARA-like: alpha all zero)", RobotModel(
            name="parallel", a=np.array([0.3, 0.25, 0.1, 0.0, 0.0, 0.05]),
            alpha=np.zeros(6), d=np.array([0.2, 0.0, 0.0, 0.0, 0.0, 0.0]),
            theta_offset=np.zeros(6),
            qlim=np.tile([[-2.5, 2.5]], (6, 1)))),
        ("intersecting axes (a = 0 throughout)", RobotModel(
            name="intersecting", a=np.zeros(6),
            alpha=np.array([np.pi / 2, -np.pi / 2, np.pi / 2, -np.pi / 2,
                            np.pi / 2, 0.0]),
            d=np.array([0.3, 0.0, 0.25, 0.0, 0.0, 0.1]),
            theta_offset=np.zeros(6),
            qlim=np.tile([[-2.5, 2.5]], (6, 1)))),
        ("mixed twists and a long tool", RobotModel(
            name="mixed", a=np.array([0.1, 0.4, 0.05, 0.0, 0.0, 0.0]),
            alpha=np.array([np.pi / 2, 0.0, np.pi / 2, -np.pi / 2,
                            np.pi / 2, 0.0]),
            d=np.array([0.4, 0.0, 0.1, 0.35, 0.0, 0.2]),
            theta_offset=np.zeros(6),
            qlim=np.tile([[-2.8, 2.8]], (6, 1)))),
    ]

    tmp = tempfile.mkdtemp(prefix="urdf_dh_")
    ok = True
    print("URDF <-> DH round trip  (FK equivalence, not table equality)")
    print(f"{'case':46s} {'pos err':>10s} {'rot err':>10s}   result")
    for label, rb in cases:
        path = to_urdf(rb, os.path.join(tmp, f"{rb.name}.urdf"))
        try:
            rb2, Tb, Tt = from_urdf(path, verify_tol=1e-8)
            joints, _, _ = parse_urdf(path)
            chain = joint_chain(joints)
            wp, wr = verify(rb2, chain, Tb, Tt, n=400, tol=1e-8, seed=5)
            print(f"{label:46s} {wp:10.2e} {wr:10.2e}   pass")
        except UrdfConversionError as e:
            print(f"{label:46s} {'-':>10s} {'-':>10s}   FAIL\n      {e}")
            ok = False

    print("\nRejection cases (these must raise, not return a wrong table)")
    bad = {
        "prismatic joint": """<robot name="p">
 <link name="a"/><link name="b"/><link name="c"/>
 <joint name="j1" type="revolute"><parent link="a"/><child link="b"/>
  <axis xyz="0 0 1"/></joint>
 <joint name="j2" type="prismatic"><parent link="b"/><child link="c"/>
  <axis xyz="0 0 1"/></joint></robot>""",
        "branching tree": """<robot name="t">
 <link name="a"/><link name="b"/><link name="c"/><link name="d"/>
 <joint name="j1" type="revolute"><parent link="a"/><child link="b"/>
  <axis xyz="0 0 1"/></joint>
 <joint name="j2" type="revolute"><parent link="b"/><child link="c"/>
  <axis xyz="0 1 0"/></joint>
 <joint name="j3" type="revolute"><parent link="b"/><child link="d"/>
  <axis xyz="1 0 0"/></joint></robot>""",
        "single joint": """<robot name="s">
 <link name="a"/><link name="b"/>
 <joint name="j1" type="revolute"><parent link="a"/><child link="b"/>
  <axis xyz="0 0 1"/></joint></robot>""",
    }
    for label, xml in bad.items():
        path = os.path.join(tmp, f"bad_{label.replace(' ', '_')}.urdf")
        with open(path, "w") as f:
            f.write(xml)
        try:
            from_urdf(path)
            print(f"  {label:44s} FAIL - returned a table")
            ok = False
        except UrdfConversionError as e:
            print(f"  {label:44s} rejected: {str(e).splitlines()[0][:60]}")

    print()
    print("ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_self_test())
