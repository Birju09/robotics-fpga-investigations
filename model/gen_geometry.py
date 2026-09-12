#!/usr/bin/env python3
"""
Emit hls/include/ik_geometry.hpp from the robot definition.

Generated rather than hand-duplicated: a silent mismatch between this and
ik_model.py would look like quantisation error, the exact quantity this
project measures. ik_config.hpp keeps implementation choices; this file
carries only what the robot is. Constants are printed with repr() precision
so the C++ double is bit-identical to numpy's.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collision as COLL  # noqa: E402
from robot import ROBOT, check_analytic_form  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hls",
                   "include", "ik_geometry.hpp")

#: Bump when a consumer requires something new; ik_config.hpp checks it.
#: 2 added the capsule collision tables (IK_CAP_*, IK_COLL_PAIRS).
VERSION = 2


def c_double(x):
    """Shortest decimal that round-trips to the same double."""
    return repr(float(x))


def carr(name, vals):
    body = ", ".join(c_double(v) for v in vals)
    return f"static const double {name}[IK_DOF] = {{{body}}};"


def capsule_block(rb):
    """
    Emit the capsule set and the checked pair list.

    A static pair list, decided here rather than at run time, is what the
    ecosystem does too: MoveIt's SRDF <disable_collisions> block and
    franka_description's *_sc sub-links are both precomputed. The pruning
    rules (and why the radii are assumptions) are in model/collision.py.
    """
    caps = COLL.capsules(rb)
    pairs, dropped = COLL.capsule_pairs(caps, explain=True)

    def vec3(v):
        #! Snapped like IK_DH_CA/SA above: these coordinates are built from
        #! cos/sin(alpha_i) at multiples of pi/2, so exact zeros arrive as
        #! 1e-17 residues. They would quantise to zero in Q16.16 anyway, but
        #! a literal 0.0 keeps the generated header readable and keeps the
        #! float and fixed builds from differing on a value that is meant to
        #! be exactly zero.
        return "{" + ", ".join(
            c_double(0.0 if abs(x) < 1e-15 else x) for x in v) + "}"

    cap_rows = "\n".join(
        f"//! {i:^3} | {c.name:<11} | {c.frame:^5} | {c.node_a}-{c.node_b} | "
        f"{c.length:8.5f} | {c.radius:6.3f}"
        for i, c in enumerate(caps))
    drop_rows = "\n".join(
        f"//!   {caps[i].name:<11} {caps[j].name:<11} {why}"
        for i, j, why in dropped)

    pa = ",\n    ".join(vec3(c.pa) for c in caps)
    pb = ",\n    ".join(vec3(c.pb) for c in caps)
    fr = ", ".join(str(c.frame) for c in caps)
    rr = ", ".join(c_double(c.radius) for c in caps)
    pr = ", ".join("{" + f"{i}, {j}" + "}" for i, j in pairs)

    return f"""
//! ---------------- capsule collision model ----------------
//
//! One capsule per non-degenerate DH limb, each RIGID in one DH frame, so it
//! has a single point Jacobian:
//
//!     o_i     = (0, 0, 0)                            in frame i
//!     P_i     = (-a_i, 0, 0)                         after Tz, before Tx
//!     o_{{i-1}} = P_i - d_i * (0, sin alpha_i, cos alpha_i)
//
//! all three constant, because theta_i only rotates frame i itself. The
//! obvious alternative - a segment between consecutive frame origins - is
//! NOT a rigid body and has no single point Jacobian. See
//! model/collision.py.
//
//! idx | name        | frame | node | length   | radius
//! ----+-------------+-------+------+----------+-------
{cap_rows}
//
//! RADII ARE CHOSEN, NOT CITED. model/collision.py:LINK_RADIUS explains
//! what that costs: a negative clearance is a real finding, a positive one
//! is not a proof of clearance.
//
//! Pairs pruned from the check, and why:
{drop_rows}
//

#define IK_CAP_COUNT {len(caps)}
#define IK_COLL_PAIR_COUNT {len(pairs)}

//! Capsule endpoints in their own frame's coordinates.
static const double IK_CAP_PA[IK_CAP_COUNT][3] = {{
    {pa}
}};
static const double IK_CAP_PB[IK_CAP_COUNT][3] = {{
    {pb}
}};

//! Frame each capsule is rigid in: an index into the DH chain, so the point
//! Jacobian's non-zero columns are joints 0 .. IK_CAP_FRAME-1.
static const int IK_CAP_FRAME[IK_CAP_COUNT] = {{{fr}}};
static const double IK_CAP_RADIUS[IK_CAP_COUNT] = {{{rr}}};

//! Capsule index pairs to check. Fixed at generation time, so the collision
//! kernel's trip count is a compile-time constant and its latency does not
//! depend on the configuration.
static const int IK_COLL_PAIRS[IK_COLL_PAIR_COUNT][2] = {{{pr}}};
"""


def build(rb):
    check_analytic_form(rb)

    ca = [0.0 if abs(np.cos(x)) < 1e-15 else np.cos(x) for x in rb.alpha]
    sa = [0.0 if abs(np.sin(x)) < 1e-15 else np.sin(x) for x in rb.alpha]

    prov = "\n".join(f"//!   {k}\n//!       {v}"
                     for k, v in rb.provenance.items())

    rows = "\n".join(
        f"//! {i + 1} | {rb.a[i]:9.5f} | {np.degrees(rb.alpha[i]):+7.1f} | "
        f"{rb.d[i]:9.5f}"
        for i in range(rb.dof))

    return f"""#ifndef IK_GEOMETRY_HPP
#define IK_GEOMETRY_HPP

//
//! GENERATED FILE - do not edit.  Produced by model/gen_geometry.py from
//! model/robot.py, the single source of truth for the kinematics.
//
//! Robot: {rb.name}
//
//! Standard (distal) DH.  T_i-1_i = Rz(theta_i) Tz(d_i) Tx(a_i) Rx(alpha_i),
//! so d_i translates along z_i-1, not z_i - that's what makes d3 a lateral
//! offset perpendicular to the arm plane.
//
//! i |     a_i | alpha_i |       d_i
//! --+---------+---------+----------
{rows}
//
//! Provenance:
{prov}
//

#define IK_GEOMETRY_VERSION {VERSION}

#define IK_DOF {rb.dof}

//! ---------------- link geometry, metres ----------------
#define IK_D1 {c_double(rb.d1)}
#define IK_A1 {c_double(rb.a1)}
#define IK_A2 {c_double(rb.a2)}
#define IK_A3 {c_double(rb.a3)}
//! Lateral shoulder offset. Keeps the wrist centre off the joint-1 axis
//! (where theta1 is undefined), turning that singularity into an
//! unreachable cylinder of this radius instead.
#define IK_D3 {c_double(rb.d3)}
#define IK_D4 {c_double(rb.d4)}
#define IK_D6 {c_double(rb.d6)}

//! ---------------- derived ----------------
//! Effective forearm length/offset angle. theta3 = gamma - phi on this arm
//! (alpha3 = {np.degrees(rb.alpha[2]):+.0f} deg flips the sense of the elbow offset).
#define IK_L3 {c_double(rb.L3)}   //! hypot(a3, d4)
#define IK_PHI {c_double(rb.phi)}  //! atan2(d4, a3)

//! Squared lateral offset, for the theta1 discriminant rho^2 - d3^2.
#define IK_D3SQ {c_double(rb.d3 ** 2)}

//! Law-of-cosines constants, pre-folded so the closed form needs no runtime
//! division: cos(gamma) = (u^2 + w^2 - (a2^2 + L3^2)) * (1 / (2*a2*L3))
#define IK_A2SQ_PLUS_L3SQ {c_double(rb.a2 ** 2 + rb.L3 ** 2)}
#define IK_INV_2A2L3 {c_double(1.0 / (2.0 * rb.a2 * rb.L3))}

//! ---------------- DH table ----------------
//! cos/sin of alpha_i are exactly 0 / +-1, tabulated rather than computed -
//! removes 12 CORDIC calls from every FK evaluation.
{carr("IK_DH_A", rb.a)}
{carr("IK_DH_D", rb.d)}
{carr("IK_DH_CA", ca)}
{carr("IK_DH_SA", sa)}

//! ---------------- joint limits, radians ----------------
//! Not enforced by the kernels. Used by gen_vectors.py to sample reachable
//! configurations.
static const double IK_QLIM_LO[IK_DOF] = {{{", ".join(c_double(v) for v in rb.qlim[:, 0])}}};
static const double IK_QLIM_HI[IK_DOF] = {{{", ".join(c_double(v) for v in rb.qlim[:, 1])}}};
{capsule_block(rb)}
#endif  //! IK_GEOMETRY_HPP
"""


def main():
    path = os.path.normpath(OUT)
    with open(path, "w") as f:
        f.write(build(ROBOT))
    print(f"  {'ik_geometry.hpp':20s} {ROBOT.name:12s} -> "
          f"{os.path.relpath(path)}")


if __name__ == "__main__":
    print("Generating HLS geometry header:")
    main()
    print("done.")
