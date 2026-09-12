#ifndef IK_GEOMETRY_HPP
#define IK_GEOMETRY_HPP

//
//! GENERATED FILE - do not edit. Produced by model/gen_geometry.py from
//! model/robot.py (source of truth). Hand edits would desync the kernels
//! from the golden model, indistinguishable from quantisation error - the
//! quantity this project measures.
//
//! Robot: PUMA 560
//
//! Standard (distal) DH.  T_i-1_i = Rz(theta_i) Tz(d_i) Tx(a_i) Rx(alpha_i),
//! so d_i translates along z_i-1, not z_i.  That is what makes d3 a lateral
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

#define IK_GEOMETRY_VERSION 1

#define IK_DOF 6

//! ---------------- link geometry, metres ----------------
#define IK_D1 0.6718
#define IK_A1 0.0
#define IK_A2 0.4318
#define IK_A3 0.0203
//! Non-zero keeps the wrist centre off the joint-1 axis (theta1 undefined
//! there), turning a reachable singularity into an unreachable cylinder.
#define IK_D3 0.15005
#define IK_D4 0.4318
#define IK_D6 0.0565

//! ---------------- derived ----------------
//! theta3 = gamma - phi (alpha3 = -90 deg flips the sense of the elbow offset).
#define IK_L3 0.4322769135635166   //! hypot(a3, d4)
#define IK_PHI 1.5238184104468135  //! atan2(d4, a3)

//! Squared lateral offset, for the theta1 discriminant rho^2 - d3^2.
#define IK_D3SQ 0.022515002499999996

//! Law-of-cosines constants, pre-folded so the closed form needs no runtime
//! division: cos(gamma) = (u^2 + w^2 - (a2^2 + L3^2)) * (1 / (2*a2*L3))
#define IK_A2SQ_PLUS_L3SQ 0.37331457
#define IK_INV_2A2L3 2.6787076895038266

//! ---------------- DH table ----------------
//! cos/sin of alpha_i are exactly 0/+-1, tabulated to avoid 12 CORDIC calls
//! per FK evaluation.
static const double IK_DH_A[IK_DOF] = {0.0, 0.4318, 0.0203, 0.0, 0.0, 0.0};
static const double IK_DH_D[IK_DOF] = {0.6718, 0.0, 0.15005, 0.4318, 0.0, 0.0565};
static const double IK_DH_CA[IK_DOF] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
static const double IK_DH_SA[IK_DOF] = {1.0, 0.0, -1.0, 1.0, -1.0, 0.0};

//! ---------------- joint limits, radians ----------------
//! Not enforced by the kernels. Used by model/gen_vectors.py to sample
//! reachable configurations; also usable to range-check a result.
static const double IK_QLIM_LO[IK_DOF] = {-2.792526803190927, -3.9269908169872414, -0.7853981633974483, -1.9198621771937625, -1.7453292519943295, -4.642575810304916};
static const double IK_QLIM_HI[IK_DOF] = {2.792526803190927, 0.7853981633974483, 3.9269908169872414, 2.9670597283903604, 1.7453292519943295, 4.642575810304916};

#endif  //! IK_GEOMETRY_HPP
