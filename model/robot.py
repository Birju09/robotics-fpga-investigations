#!/usr/bin/env python3
"""
The robot definition layer: one declarative source for the kinematics.

Everything downstream reads its geometry from here - `ik_model` for the golden
model, `gen_geometry` for the C++ header the kernels compile against,
`gen_vectors` for the workloads, `plot_dh` for the figure.  Before this module
existed the DH table was a block of module constants in `ik_model.py` and a
parallel block of `#define`s in `hls/include/ik_config.hpp`, kept in agreement
by hand.

A `RobotModel` carries more than the numbers.  It carries where they came from
(`provenance`) and what structure a closed-form solver may assume about them
(`analytic_form`), because both are load-bearing and neither is recoverable
from a table of floats:

  - Provenance separates a published, citable parameter from a value someone
    picked.  The two are mixed in any real robot description and a report that
    cannot tell them apart is overclaiming.
  - The structural assumptions are what `iks::analytic()` is entitled to rely
    on.  A closed-form solution is not generic over DH tables; it is derived
    for a particular pattern of zeros.  `check_analytic_form()` tests that
    pattern, so a table that silently breaks the solver fails here rather than
    as wrong joint angles later.

`from_urdf()` in `urdf_to_dh.py` produces one of these, which is the point of
the split: a URDF-derived robot is not a special case downstream.

Conventions are standard (distal) DH, matching `ik_model.dh_mat()`:

    T_i-1_i = Rz(theta_i + theta_offset_i) Tz(d_i) Tx(a_i) Rx(alpha_i)

Note that d_i translates along z_{i-1}, not z_i.  That is what makes d3 a
lateral offset perpendicular to the arm plane rather than a length along the
forearm, and it is the whole reason the offset removes the shoulder
singularity.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class RobotModel:
    name: str
    a: np.ndarray             # link lengths, metres
    alpha: np.ndarray         # link twists, radians
    d: np.ndarray             # link offsets, metres
    theta_offset: np.ndarray  # added to the joint variable before Rz
    qlim: np.ndarray          # (n, 2) lower/upper joint limits, radians
    provenance: dict = field(default_factory=dict)

    # ---- derived, all read off the table rather than restated ----
    @property
    def dof(self):
        return len(self.a)

    @property
    def d1(self):
        """Base offset along z0: shoulder height above the mounting face."""
        return float(self.d[0])

    @property
    def a1(self):
        """In-plane shoulder offset, joint-1 axis to joint-2 axis."""
        return float(self.a[0])

    @property
    def a2(self):
        """Upper arm."""
        return float(self.a[1])

    @property
    def a3(self):
        """Elbow offset."""
        return float(self.a[2])

    @property
    def d3(self):
        """Lateral shoulder offset, perpendicular to the arm plane.

        Non-zero is what keeps the wrist centre off the joint-1 axis, where
        theta1 would be undefined; it converts a reachable singularity into an
        unreachable cylinder of radius |d3|.  See check_analytic_form().
        """
        return float(self.d[2])

    @property
    def d4(self):
        """Forearm."""
        return float(self.d[3])

    @property
    def d6(self):
        """Wrist centre to tool point along the approach axis."""
        return float(self.d[5])

    @property
    def L3(self):
        """Effective forearm: the elbow offset folded into one length."""
        return float(np.hypot(self.a3, self.d4))

    @property
    def phi(self):
        """Angle of L3 off the a3 direction; theta3 = gamma + phi."""
        return float(np.arctan2(self.d4, self.a3))

    @property
    def reach_outer(self):
        """Largest wrist-centre distance from the shoulder point."""
        return self.a2 + self.L3

    @property
    def reach_inner(self):
        """Smallest ditto: the radius of the unreachable core."""
        return abs(self.a2 - self.L3)

    def summary(self):
        return "\n".join([
            f"{self.name}  ({self.dof} DOF)",
            "   i        a_i     alpha_i         d_i     qlim (deg)",
            *[f"   {i + 1}  {self.a[i]:9.5f}  {np.degrees(self.alpha[i]):+7.1f}"
              f"  {self.d[i]:10.5f}   "
              f"{np.degrees(self.qlim[i, 0]):+7.1f} .. "
              f"{np.degrees(self.qlim[i, 1]):+7.1f}"
              for i in range(self.dof)],
            f"   L3 = {self.L3:.5f}   phi = {np.degrees(self.phi):.3f} deg",
            f"   wrist-centre reach {self.reach_inner:.4f} .. "
            f"{self.reach_outer:.4f} m",
        ])


class AnalyticFormError(ValueError):
    """The table breaks a structural assumption the closed-form solver makes."""


def check_analytic_form(rb):
    """
    Verify the pattern of zeros `iks::analytic()` is derived against.

    The closed-form solution is not generic over DH tables.  It requires a
    6R arm whose last three axes intersect at a point (a spherical wrist), so
    that position and orientation decouple:

        a4 = a5 = a6 = 0 and d5 = 0     the three wrist axes meet
        alpha4, alpha5 = +-pi/2         and meet mutually perpendicular
        alpha2 = 0                      shoulder and elbow axes parallel,
                                        which makes joints 2-3 planar
        d2 = 0                          the only lateral offset is d3

    Raises AnalyticFormError naming the violated assumption.  DLS needs none
    of this - it differentiates whatever fk_all() computes - so a table that
    fails here is still usable by ik_dls, and the caller may choose to carry
    on with the iterative solver alone.
    """
    if rb.dof != 6:
        raise AnalyticFormError(f"{rb.dof} DOF; the closed form is 6R only")

    for i in (3, 4, 5):
        if abs(rb.a[i]) > 1e-12:
            raise AnalyticFormError(
                f"a{i + 1} = {rb.a[i]} is non-zero: the wrist axes do not "
                "intersect, so position and orientation do not decouple")
    if abs(rb.d[4]) > 1e-12:
        raise AnalyticFormError(
            f"d5 = {rb.d[4]} is non-zero: the wrist is offset, not spherical")
    if abs(rb.alpha[1]) > 1e-12:
        raise AnalyticFormError(
            f"alpha2 = {rb.alpha[1]} is non-zero: joints 2 and 3 are not "
            "parallel, so the position sub-problem is not planar")
    if abs(rb.d[1]) > 1e-12:
        raise AnalyticFormError(
            f"d2 = {rb.d[1]} is non-zero; the solver folds the whole lateral "
            "offset into d3")
    for i in (3, 4):
        if abs(abs(rb.alpha[i]) - np.pi / 2) > 1e-9:
            raise AnalyticFormError(
                f"alpha{i + 1} = {rb.alpha[i]} is not +-pi/2; the wrist is "
                "not a ZYZ Euler set")
    return True


# ---------------------------------------------------------------------------
# PUMA 560
# ---------------------------------------------------------------------------
# Standard (distal) DH, as tabulated in Corke's Robotics Toolbox `mdl_puma560`
# and equivalent to the modified-DH table in Craig, "Introduction to Robotics",
# section 3.6.  The arm was designed in inches and the canonical values are
# exact conversions, which is a useful check that a table has not been
# corrupted in transcription:
#
#     a2 = 431.8 mm  = 17.00 in
#     a3 =  20.32 mm =  0.80 in
#     d3 = 149.09 mm =  5.87 in
#     d4 = 433.07 mm = 17.05 in
#
# d1 and d6 are NOT part of that kinematic table, where both are zero: the
# canonical table places frame 0 at the shoulder and the tool point at the
# wrist centre.  They are added here, and flagged as additions in
# `provenance`, because this project needs a mounting face to measure poses
# against and a non-degenerate approach axis for the orientation half of the
# IK problem to mean anything.  Quoting them as PUMA parameters would be
# wrong; they are this project's, applied to PUMA's arm.
PUMA560 = RobotModel(
    name="PUMA 560",
    a=np.array([0.0, 0.4318, 0.0203, 0.0, 0.0, 0.0]),
    alpha=np.array([np.pi / 2, 0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, 0.0]),
    d=np.array([0.6718, 0.0, 0.15005, 0.4318, 0.0, 0.0565]),
    theta_offset=np.zeros(6),
    # Manufacturer joint limits, degrees, as tabulated by Corke.  Asymmetric,
    # unlike the symmetric +-limits this project used before, and that
    # asymmetry is real: the elbow cannot swing equally both ways.  They also
    # do the job a collision model would otherwise have to - the q3 range is
    # what stops the forearm folding back onto the upper arm, which bare link
    # lengths permit (reach_inner is only 19 mm).
    qlim=np.radians(np.array([
        [-160.0, 160.0],
        [-225.0, 45.0],
        [-45.0, 225.0],
        [-110.0, 170.0],
        [-100.0, 100.0],
        [-266.0, 266.0],
    ])),
    provenance={
        "a, alpha, d3, d4": (
            "PUMA 560 standard-DH table, Corke Robotics Toolbox "
            "mdl_puma560; equivalent to Craig section 3.6. Inch-exact."),
        "d1 = 0.6718": (
            "pedestal height. Not in the kinematic table (which has d1 = 0). "
            "Folded into d1 here so the chain stays pure DH and poses have a "
            "mounting face to be measured against."),
        "d6 = 0.0565": (
            "tool offset along the approach axis. Not in the kinematic table "
            "(which has d6 = 0, putting the tool point at the wrist centre). "
            "Added so the orientation half of the IK problem is "
            "non-degenerate."),
        "qlim": "PUMA 560 manufacturer joint limits, as tabulated by Corke.",
    },
)

# What the rest of the tree uses.  Single assignment, so retargeting the
# project at another arm - including one returned by urdf_to_dh.from_urdf() -
# is this line.
ROBOT = PUMA560


if __name__ == "__main__":
    print(ROBOT.summary())
    print()
    check_analytic_form(ROBOT)
    print("analytic form: ok (spherical wrist, planar shoulder/elbow)")
    print(f"lateral offset d3 = {ROBOT.d3:.5f} m -> the wrist centre cannot "
          f"reach the joint-1 axis")
    print()
    for k, v in ROBOT.provenance.items():
        print(f"  {k}\n      {v}")
