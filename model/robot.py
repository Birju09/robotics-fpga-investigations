#!/usr/bin/env python3
"""
Declarative robot/DH definition, single source of truth for `ik_model`,
`gen_geometry`, `gen_vectors`, and `plot_dh`.

`RobotModel` carries `provenance` (which numbers are cited vs. chosen) and
`analytic_form` assumptions (the zero pattern `iks::analytic()` relies on;
checked by `check_analytic_form()`), since neither is recoverable from a bare
table of floats. `from_urdf()` in `urdf_to_dh.py` also produces one of these.

Standard (distal) DH, matching `ik_model.dh_mat()`:

    T_i-1_i = Rz(theta_i + theta_offset_i) Tz(d_i) Tx(a_i) Rx(alpha_i)

d_i translates along z_{i-1}, not z_i - that's what makes d3 a lateral offset
perpendicular to the arm plane, removing the shoulder singularity.
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
        """Lateral shoulder offset, perpendicular to arm plane. Non-zero keeps
        the wrist centre off the joint-1 axis (theta1 undefined there),
        turning that singularity into an unreachable cylinder of radius |d3|.
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
    Verify the zero pattern `iks::analytic()` requires: a 6R arm with a
    spherical wrist (last three axes meet at a point) so position and
    orientation decouple:

        a4 = a5 = a6 = 0 and d5 = 0     wrist axes meet
        alpha4, alpha5 = +-pi/2         mutually perpendicular
        alpha2 = 0                      shoulder/elbow axes parallel (planar)
        d2 = 0                          only lateral offset is d3

    Raises AnalyticFormError naming the violated assumption. ik_dls doesn't
    need any of this, so a table that fails here can still use it.
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
# Standard DH per Corke's Robotics Toolbox `mdl_puma560` (== Craig sec. 3.6).
# Designed in inches; canonical values are exact conversions, a useful
# transcription check:
#     a2 = 431.8 mm = 17.00 in   a3 = 20.32 mm = 0.80 in
#     d3 = 149.09 mm = 5.87 in   d4 = 433.07 mm = 17.05 in
#
# d1, d6 are NOT in the canonical table (both zero there: frame 0 at shoulder,
# tool point at wrist centre). Added here (flagged in `provenance`) so poses
# have a mounting face and the orientation IK is non-degenerate - they are
# this project's additions, not PUMA parameters.
PUMA560 = RobotModel(
    name="PUMA 560",
    a=np.array([0.0, 0.4318, 0.0203, 0.0, 0.0, 0.0]),
    alpha=np.array([np.pi / 2, 0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, 0.0]),
    d=np.array([0.6718, 0.0, 0.15005, 0.4318, 0.0, 0.0565]),
    theta_offset=np.zeros(6),
    # Manufacturer limits (Corke), asymmetric (real: elbow doesn't swing
    # equally both ways). q3's range also substitutes for a collision model -
    # it stops the forearm folding onto the upper arm (reach_inner = 19 mm).
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

# Single assignment: retargeting the project at another arm is this line.
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
