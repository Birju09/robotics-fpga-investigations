#ifndef IK_GEOMETRY_HPP
#define IK_GEOMETRY_HPP

//
//! GENERATED FILE - do not edit.  Produced by model/gen_geometry.py from
//! model/robot.py, the single source of truth for the kinematics.
//
//! Robot: PUMA 560
//
//! Standard (distal) DH.  T_i-1_i = Rz(theta_i) Tz(d_i) Tx(a_i) Rx(alpha_i),
//! so d_i translates along z_i-1, not z_i - that's what makes d3 a lateral
//! offset perpendicular to the arm plane.
//
//! i |     a_i | alpha_i |       d_i
//! --+---------+---------+----------
//! 1 |   0.00000 |   +90.0 |   0.67180
//! 2 |   0.43180 |    +0.0 |   0.00000
//! 3 |   0.02030 |   -90.0 |   0.15005
//! 4 |   0.00000 |   +90.0 |   0.43180
//! 5 |   0.00000 |   -90.0 |   0.00000
//! 6 |   0.00000 |    +0.0 |   0.05650
//
//! Provenance:
//!   a, alpha, d3, d4
//!       PUMA 560 standard-DH table, Corke Robotics Toolbox mdl_puma560; equivalent to Craig section 3.6. Inch-exact.
//!   d1 = 0.6718
//!       pedestal height. Not in the kinematic table (which has d1 = 0). Folded into d1 here so the chain stays pure DH and poses have a mounting face to be measured against.
//!   d6 = 0.0565
//!       tool offset along the approach axis. Not in the kinematic table (which has d6 = 0, putting the tool point at the wrist centre). Added so the orientation half of the IK problem is non-degenerate.
//!   qlim
//!       PUMA 560 manufacturer joint limits, as tabulated by Corke.
//

#define IK_GEOMETRY_VERSION 2

#define IK_DOF 6

//! ---------------- link geometry, metres ----------------
#define IK_D1 0.6718
#define IK_A1 0.0
#define IK_A2 0.4318
#define IK_A3 0.0203
//! Lateral shoulder offset. Keeps the wrist centre off the joint-1 axis
//! (where theta1 is undefined), turning that singularity into an
//! unreachable cylinder of this radius instead.
#define IK_D3 0.15005
#define IK_D4 0.4318
#define IK_D6 0.0565

//! ---------------- derived ----------------
//! Effective forearm length/offset angle. theta3 = gamma - phi on this arm
//! (alpha3 = -90 deg flips the sense of the elbow offset).
#define IK_L3 0.4322769135635166   //! hypot(a3, d4)
#define IK_PHI 1.5238184104468135  //! atan2(d4, a3)

//! Squared lateral offset, for the theta1 discriminant rho^2 - d3^2.
#define IK_D3SQ 0.022515002499999996

//! Law-of-cosines constants, pre-folded so the closed form needs no runtime
//! division: cos(gamma) = (u^2 + w^2 - (a2^2 + L3^2)) * (1 / (2*a2*L3))
#define IK_A2SQ_PLUS_L3SQ 0.37331457
#define IK_INV_2A2L3 2.6787076895038266

//! ---------------- DH table ----------------
//! cos/sin of alpha_i are exactly 0 / +-1, tabulated rather than computed -
//! removes 12 CORDIC calls from every FK evaluation.
static const double IK_DH_A[IK_DOF] = {0.0, 0.4318, 0.0203, 0.0, 0.0, 0.0};
static const double IK_DH_D[IK_DOF] = {0.6718, 0.0, 0.15005, 0.4318, 0.0, 0.0565};
static const double IK_DH_CA[IK_DOF] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
static const double IK_DH_SA[IK_DOF] = {1.0, 0.0, -1.0, 1.0, -1.0, 0.0};

//! ---------------- joint limits, radians ----------------
//! Not enforced by the kernels. Used by gen_vectors.py to sample reachable
//! configurations.
static const double IK_QLIM_LO[IK_DOF] = {-2.792526803190927, -3.9269908169872414, -0.7853981633974483, -1.9198621771937625, -1.7453292519943295, -4.642575810304916};
static const double IK_QLIM_HI[IK_DOF] = {2.792526803190927, 0.7853981633974483, 3.9269908169872414, 2.9670597283903604, 1.7453292519943295, 4.642575810304916};

//! ---------------- capsule collision model ----------------
//
//! One capsule per non-degenerate DH limb, each RIGID in one DH frame, so it
//! has a single point Jacobian:
//
//!     o_i     = (0, 0, 0)                            in frame i
//!     P_i     = (-a_i, 0, 0)                         after Tz, before Tx
//!     o_{i-1} = P_i - d_i * (0, sin alpha_i, cos alpha_i)
//
//! all three constant, because theta_i only rotates frame i itself. The
//! obvious alternative - a segment between consecutive frame origins - is
//! NOT a rigid body and has no single point Jacobian. See
//! model/collision.py.
//
//! idx | name        | frame | node | length   | radius
//! ----+-------------+-------+------+----------+-------
//!  0  | link1_d     |   1   | 0-1 |  0.67180 |  0.090
//!  1  | link2_a     |   2   | 1-2 |  0.43180 |  0.075
//!  2  | link3_d     |   3   | 2-3 |  0.15005 |  0.070
//!  3  | link3_a     |   3   | 3-4 |  0.02030 |  0.070
//!  4  | link4_d     |   4   | 4-5 |  0.43180 |  0.055
//!  5  | link6_d     |   6   | 5-6 |  0.05650 |  0.040
//
//! RADII ARE CHOSEN, NOT CITED. model/collision.py:LINK_RADIUS explains
//! what that costs: a negative clearance is a real finding, a positive one
//! is not a proof of clearance.
//
//! Pairs pruned from the check, and why:
//!   link1_d     link2_a     shares a joint
//!   link2_a     link3_d     shares a joint
//!   link3_d     link3_a     shares a joint
//!   link3_d     link4_d     unseparable: path 0.0203 m <= radii 0.1250 m
//!   link3_a     link4_d     shares a joint
//!   link4_d     link6_d     shares a joint
//

#define IK_CAP_COUNT 6
#define IK_COLL_PAIR_COUNT 9

//! Capsule endpoints in their own frame's coordinates.
static const double IK_CAP_PA[IK_CAP_COUNT][3] = {
    {0.0, -0.6718, 0.0},
    {-0.4318, 0.0, 0.0},
    {-0.0203, 0.15005, 0.0},
    {-0.0203, 0.0, 0.0},
    {0.0, -0.4318, 0.0},
    {0.0, 0.0, -0.0565}
};
static const double IK_CAP_PB[IK_CAP_COUNT][3] = {
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {-0.0203, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0}
};

//! Frame each capsule is rigid in: an index into the DH chain, so the point
//! Jacobian's non-zero columns are joints 0 .. IK_CAP_FRAME-1.
static const int IK_CAP_FRAME[IK_CAP_COUNT] = {1, 2, 3, 3, 4, 6};
static const double IK_CAP_RADIUS[IK_CAP_COUNT] = {0.09, 0.075, 0.07, 0.07, 0.055, 0.04};

//! Squared segment length and its reciprocal, per capsule.
//
//! These are CONSTANTS because a capsule is rigid: |pb - pa| does not depend
//! on q. That removes two of the three reciprocals a segment-segment query
//! would otherwise need - the a = d1.d1 and e = d2.d2 divisors in the clamp
//! cascade are table lookups, and only the denom = a*e - b*b divisor is
//! computed at run time. On this arm that is 9 runtime reciprocals per
//! query instead of 27.
static const double IK_CAP_LEN_SQ[IK_CAP_COUNT] = {0.45131523999999995, 0.18645124000000002, 0.022515002499999996, 0.00041208999999999994, 0.18645124000000002, 0.00319225};
static const double IK_CAP_INV_LEN_SQ[IK_CAP_COUNT] = {2.2157461378880097, 5.36333252597301, 44.41482962304802, 2426.6543716178508, 5.36333252597301, 313.25867334951835};

//! Capsule index pairs to check. Fixed at generation time, so the collision
//! kernel's trip count is a compile-time constant and its latency does not
//! depend on the configuration.
static const int IK_COLL_PAIRS[IK_COLL_PAIR_COUNT][2] = {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 5}, {3, 5}};

#endif  //! IK_GEOMETRY_HPP
