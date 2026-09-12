#!/usr/bin/env python3
"""
Emit hls/include/ik_geometry.hpp from the robot definition.

The geometry used to exist twice: as module constants in ik_model.py and as a
block of #defines in ik_config.hpp, agreeing only because someone kept them in
agreement.  That is a bad arrangement for a project whose entire output is
a comparison between a Python reference and a fixed-point kernel - a silent
divergence there does not look like a build error, it looks like quantisation
error, which is exactly the quantity being measured.

So the C++ side is generated.  ik_config.hpp keeps everything that is a
genuine implementation choice (initiation intervals, epsilons, solver
defaults); this file carries only what the robot is.

The emitted constants are printed with repr() precision so the double the C++
compiler sees is bit-identical to the one numpy used.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from robot import ROBOT, check_analytic_form  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hls",
                   "include", "ik_geometry.hpp")

#: Bump when a consumer might require something this file did not used to
#: define.  ik_config.hpp checks it, the same way main.c checks
#: IK_VECTORS_VERSION.
VERSION = 1


def c_double(x):
    """Shortest decimal that round-trips to the same double."""
    return repr(float(x))


def carr(name, vals):
    body = ", ".join(c_double(v) for v in vals)
    return f"static const double {name}[IK_DOF] = {{{body}}};"


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
//! GENERATED FILE - do not edit.  Produced by model/gen_geometry.py from the
//! robot definition in model/robot.py, which is the single source of truth for
//! the kinematics.  Editing this file by hand puts the kernels out of step
//! with the golden model, and the resulting disagreement is indistinguishable
//! from quantisation error - which is the quantity this project measures.
//
//! Robot: {rb.name}
//
//! Standard (distal) DH.  T_i-1_i = Rz(theta_i) Tz(d_i) Tx(a_i) Rx(alpha_i),
//! so d_i translates along z_i-1, not z_i.  That is what makes d3 a lateral
//! offset perpendicular to the arm plane.
//
//! i |     a_i | alpha_i |       d_i
//! --+---------+---------+----------
{rows}
//
//! Provenance:
{prov}
//

//! Bumped when this file grows something a consumer might require.
#define IK_GEOMETRY_VERSION {VERSION}

#define IK_DOF {rb.dof}

//! ---------------- link geometry, metres ----------------
#define IK_D1 {c_double(rb.d1)}
#define IK_A1 {c_double(rb.a1)}
#define IK_A2 {c_double(rb.a2)}
#define IK_A3 {c_double(rb.a3)}
//! Lateral shoulder offset, perpendicular to the arm plane.  Non-zero is what
//! keeps the wrist centre off the joint-1 axis, where theta1 is undefined; it
//! turns a reachable singularity into an unreachable cylinder of this radius.
#define IK_D3 {c_double(rb.d3)}
#define IK_D4 {c_double(rb.d4)}
#define IK_D6 {c_double(rb.d6)}

//! ---------------- derived ----------------
//! Effective forearm length and its offset angle.  theta3 = gamma - phi on
//! this arm (alpha3 = {np.degrees(rb.alpha[2]):+.0f} deg flips the sense of the elbow offset).
#define IK_L3 {c_double(rb.L3)}   //! hypot(a3, d4)
#define IK_PHI {c_double(rb.phi)}  //! atan2(d4, a3)

//! Squared lateral offset, for the theta1 discriminant rho^2 - d3^2.
#define IK_D3SQ {c_double(rb.d3 ** 2)}

//! Law-of-cosines constants for the analytic solver, pre-folded so the closed
//! form needs no runtime division at all:
//! cos(gamma) = (u^2 + w^2 - (a2^2 + L3^2)) * (1 / (2*a2*L3))
#define IK_A2SQ_PLUS_L3SQ {c_double(rb.a2 ** 2 + rb.L3 ** 2)}
#define IK_INV_2A2L3 {c_double(1.0 / (2.0 * rb.a2 * rb.L3))}

//! ---------------- DH table ----------------
//! cos/sin of alpha_i are exactly 0 / +-1, so they are tabulated rather than
//! computed - this removes 12 CORDIC calls from every FK evaluation.
{carr("IK_DH_A", rb.a)}
{carr("IK_DH_D", rb.d)}
{carr("IK_DH_CA", ca)}
{carr("IK_DH_SA", sa)}

//! ---------------- joint limits, radians ----------------
//! Not enforced by the kernels, which solve whatever pose they are given.
//! Used by model/gen_vectors.py to sample reachable configurations, and worth
//! having here so a consumer can range-check a result.
static const double IK_QLIM_LO[IK_DOF] = {{{", ".join(c_double(v) for v in rb.qlim[:, 0])}}};
static const double IK_QLIM_HI[IK_DOF] = {{{", ".join(c_double(v) for v in rb.qlim[:, 1])}}};

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
