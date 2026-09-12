#ifndef COLL_HPP
#define COLL_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! Capsule self-collision distance.
//
//! Golden reference is model/collision.py; the capsule set and the checked
//! pair list are generated into ik_geometry.hpp (IK_CAP_*, IK_COLL_PAIRS)
//! from model/robot.py, so nothing here is hand-kept.
//
//! Stage 0 of docs/collision_aware_ik.md: a CHECK, not yet a steering term.
//! Nothing in this header modifies a solve. Its purpose is to establish what
//! a fixed-point capsule distance query costs in area and latency before the
//! nullspace work commits to it.
//
//! ---------------- determinism ----------------
//
//! Every loop here has a compile-time trip count: IK_CAP_COUNT capsules,
//! IK_COLL_PAIR_COUNT pairs, IK_DOF chain steps, all fixed by the generated
//! tables. The segment-segment clamp cascade is written with `select`
//! (ternary) rather than early exits, so every branch is evaluated and the
//! result muxed - the same idiom iks::dls() uses deliberately for its
//! trust-region shift ("shift count is always computed and applied... so
//! per-iteration latency stays data-independent") and its singular flag
//! ("flagged rather than broken out of").
//
//! So this adds a CONSTANT to whatever calls it. That is the property the
//! project is organised around: it keeps latency_ns = c0 + c1 * iterations
//! intact in form, moving only c1.
//
//! ---------------- what a clearance here does and does not mean ----------
//
//! Radii come from model/collision.py:LINK_RADIUS and are CHOSEN, NOT CITED
//! - no published PUMA 560 cross-section table is used anywhere in this
//! project. A NEGATIVE clearance is a real finding: the zero-radius skeleton
//! is closer than any plausible link thickness allows. A positive one is not
//! a proof of clearance.
//

namespace coll {

    //! Capsule endpoints in the base frame, at configuration q.
    //!
    //! Walks the DH chain once and transforms each capsule's (constant)
    //! local endpoints as its frame is reached, rather than retaining all
    //! seven frame rotations - 36 stored values instead of 84.
    //!
    //! Each capsule is rigid in exactly one DH frame, which is what makes
    //! this possible and is not true of the obvious frame-origin-to-
    //! frame-origin skeleton. See model/collision.py's module docstring.
    void endpoints(const ik_real_t q[IK_DOF],
                   ik_real_t A[IK_CAP_COUNT][3],
                   ik_real_t B[IK_CAP_COUNT][3]);

    //! Closest points between segments [p0,p1] and [q0,q1].
    //!
    //! @param a,e          squared segment lengths, |p1-p0|^2 and |q1-q0|^2
    //! @param inv_a,inv_e  their reciprocals
    //! @param wa,wb        witness points, on the first and second segment
    //! @return             the distance between them
    //!
    //! The witness points, not just the distance, are the output that
    //! matters: they are what converts a distance into a joint-space
    //! gradient in Stage 1, and retrofitting them later would mean
    //! rewriting the reduction in min_distance().
    //!
    //! a/e and their reciprocals are PARAMETERS rather than computed here
    //! because a capsule is rigid: |pb - pa| does not depend on q, so both
    //! are generated constants (IK_CAP_LEN_SQ / IK_CAP_INV_LEN_SQ). That
    //! leaves one runtime reciprocal per pair instead of three.
    ik_real_t seg_seg(const ik_real_t p0[3], const ik_real_t p1[3],
                      const ik_real_t q0[3], const ik_real_t q1[3],
                      ik_real_t a, ik_real_t inv_a, ik_real_t e,
                      ik_real_t inv_e, ik_real_t wa[3], ik_real_t wb[3]);

    //! Minimum signed clearance over IK_COLL_PAIRS, and the witness data for
    //! the pair achieving it.
    //!
    //! Clearance is segment distance minus both radii, so it goes NEGATIVE
    //! on interpenetration. The sign carries information and is deliberately
    //! not clamped: "how far in" is what a repulsion term needs.
    //!
    //! @param pair   index into IK_COLL_PAIRS of the achieving pair
    //! @param wa,wb  its witness points, in the base frame
    void min_distance(const ik_real_t A[IK_CAP_COUNT][3],
                      const ik_real_t B[IK_CAP_COUNT][3], ik_real_t* d_min,
                      int* pair, ik_real_t wa[3], ik_real_t wb[3]);

}  // namespace coll

//! AXI4-Lite IP wrapper. Q16.16 fixed-point.
//!
//! d_min is signed: negative means the capsules interpenetrate by that much.
//! status is IK_OK, or IK_ERR_COLLISION when d_min < 0 - a flag for the
//! caller's convenience, not a refusal; d_min is always valid.
extern "C" void coll_dist_kernel(const ik_word_t q[IK_DOF], ik_word_t* d_min,
                                 int* pair, ik_word_t wa[3], ik_word_t wb[3],
                                 int* status);

#endif  //! COLL_HPP
