#include "coll.hpp"

#include "ik_math.hpp"
#include "kinematics.hpp"

//
//! Capsule self-collision distance. Golden reference: model/collision.py.
//

//! Clamp to [0,1]. The segment parameters s and t are barycentric along
//! their segments, so this is where "closest point on the infinite line"
//! becomes "closest point on the segment".
static inline ik_real_t clamp01(ik_real_t v) {
#pragma HLS INLINE
    if (v < (ik_real_t)0)
        return (ik_real_t)0;
    if (v > (ik_real_t)1)
        return (ik_real_t)1;
    return v;
}

//! ------------------------------------------------------------------
void coll::endpoints(const ik_real_t q[IK_DOF], ik_real_t A[IK_CAP_COUNT][3],
                     ik_real_t B[IK_CAP_COUNT][3]) {
#pragma HLS INLINE off
    ik_real_t R[3][3], p[3];
#pragma HLS ARRAY_PARTITION variable = R complete dim = 0
#pragma HLS ARRAY_PARTITION variable = p complete dim = 1

CE_I:
    for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
        p[r] = (ik_real_t)0;
        for (int c = 0; c < 3; c++) {
#pragma HLS UNROLL
            R[r][c] = (r == c) ? (ik_real_t)1 : (ik_real_t)0;
        }
    }

//! One pass down the chain. Each capsule is transformed as its own frame is
//! reached, so only the 2*IK_CAP_COUNT base-frame endpoints are retained
//! rather than all IK_DOF+1 frame rotations - 36 values instead of 84.
//
//! Frame 0 carries no capsule by construction (IK_CAP_FRAME >= 1), so the
//! chain step comes first.
//
//! PIPELINE off - see fk()'s FK_CHAIN: auto-pipelining would force
//! dh_step()'s rolled matrix multiply to unroll, undoing the DSP saving.
CE_CHAIN:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE off
        ikk::dh_step(i, q[i], R, p);
        const int f = i + 1;

    //! Masked, not indexed: every capsule is visited at every step and the
    //! write is predicated, so the trip count is IK_DOF * IK_CAP_COUNT
    //! whatever the frame assignment is. Indexing by frame would need a
    //! variable-bound inner loop and make latency depend on the table.
    CE_CAP:
        for (int c = 0; c < IK_CAP_COUNT; c++) {
#pragma HLS UNROLL
            const bool mine = (IK_CAP_FRAME[c] == f);
        CE_ROW:
            for (int r = 0; r < 3; r++) {
#pragma HLS UNROLL
                ik_acc_t sa = (ik_acc_t)p[r];
                ik_acc_t sb = (ik_acc_t)p[r];
                for (int k = 0; k < 3; k++) {
#pragma HLS UNROLL
                    sa += (ik_acc_t)(R[r][k] * (ik_real_t)IK_CAP_PA[c][k]);
                    sb += (ik_acc_t)(R[r][k] * (ik_real_t)IK_CAP_PB[c][k]);
                }
                if (mine) {
                    A[c][r] = (ik_real_t)sa;
                    B[c][r] = (ik_real_t)sb;
                }
            }
        }
    }
}

//! ------------------------------------------------------------------
//! Segment-segment distance with witness points.
//
//! Mirrors model/collision.py:seg_seg_witness() exactly, but written as a
//! select cascade: every branch is EVALUATED and the result muxed, so the
//! instruction count is the same for every input. Early exits would make
//! the latency of a collision query depend on the configuration, which is
//! the one thing this project cannot afford - cf. iks::dls()'s trust-region
//! shift ("always computed and applied") and its singular flag ("flagged
//! rather than broken out of").
//
//! ---- fixed-point notes, the part that is NOT a transcription ----
//
//! denom = a*e - b*b cancels catastrophically for near-parallel segments.
//! It is formed in ik_acc_t, where each product of two Q16.16 values is
//! EXACT in Q32.32 and so is their difference - so the only error in denom
//! is the rounding already present in a, e and b, not new cancellation
//! noise from the subtraction. Same reasoning as ik_dls's err_sq and
//! matmul's MM_DOT.
//
//! a and e are NOT computed here. A capsule is rigid, so |pb - pa| is
//! constant and both the squares and their reciprocals are generated
//! constants (IK_CAP_LEN_SQ / IK_CAP_INV_LEN_SQ). That removes two of the
//! three reciprocals per pair - only denom's is computed at run time.
//
//! s can transiently exceed the segment before clamping: with denom at the
//! IK_PIVOT_EPS floor the quotient reaches ~1e4, which Q16.16 holds. Should
//! a pathological case exceed 32767, ik_real_t SATURATES (AP_SAT) rather
//! than wrapping, and clamp01 then maps it to 1 - degraded, not sign-
//! flipped. That is why the storage type's saturation is load-bearing here
//! and not only in ik_types.hpp's joint-angle argument.
ik_real_t coll::seg_seg(const ik_real_t p0[3], const ik_real_t p1[3],
                        const ik_real_t q0[3], const ik_real_t q1[3],
                        ik_real_t a, ik_real_t inv_a, ik_real_t e,
                        ik_real_t inv_e, ik_real_t wa[3], ik_real_t wb[3]) {
#pragma HLS INLINE off
    ik_real_t d1[3], d2[3], r[3];
#pragma HLS ARRAY_PARTITION variable = d1 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = d2 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = r complete dim = 1

SS_DIFF:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        d1[i] = (ik_real_t)((ik_acc_t)p1[i] - (ik_acc_t)p0[i]);
        d2[i] = (ik_real_t)((ik_acc_t)q1[i] - (ik_acc_t)q0[i]);
        r[i] = (ik_real_t)((ik_acc_t)p0[i] - (ik_acc_t)q0[i]);
    }

    ik_acc_t bb = (ik_acc_t)0, cc = (ik_acc_t)0, ff = (ik_acc_t)0;
SS_DOT:
    for (int i = 0; i < 3; i++) {
#pragma HLS PIPELINE II = 1
        bb += (ik_acc_t)(d1[i] * d2[i]);
        cc += (ik_acc_t)(d1[i] * r[i]);
        ff += (ik_acc_t)(d2[i] * r[i]);
    }
    const ik_real_t b = (ik_real_t)bb;
    const ik_real_t c = (ik_real_t)cc;
    const ik_real_t f = (ik_real_t)ff;

    //! ---- s = (b f - c e) / (a e - b b), without ever dividing out of
    //! ---- range
    //
    //! Both are formed in ik_acc_t, where each product of two Q16.16 values
    //! is EXACT in Q32.32 and so is the difference - the only error is the
    //! rounding already in a, e, b, not new cancellation noise from the
    //! subtraction. That matters because den = a e sin^2(theta) collapses
    //! for near-parallel segments.
    //
    //! The quotient is NOT computed with ikm::recip(). recip() is built for
    //! spd::solve()'s pivots, which damping bounds below; den has no such
    //! bound, and 1/den overflows Q16.16 as soon as den < 3e-5. An earlier
    //! version guarded that by declaring |den| < IK_PIVOT_EPS "parallel" and
    //! forcing s = 0. That is WRONG, and wrong in the dangerous direction:
    //! den = a e sin^2(theta) puts the 1e-4 threshold at 3 degrees off
    //! parallel, and forcing s = 0 there reports a clearance up to 21 mm
    //! LARGER than the true one - it says clear when it is not. Measured by
    //! sweeping the threshold against brute force in model/validate.py.
    //
    //! Instead: compare |num| against |den| first. The clamped cases need no
    //! division at all, and what remains has |quotient| < 1 by construction,
    //! so the divider is only ever asked for a value it can represent. No
    //! threshold, no tuning constant, and it is exact in both builds - the
    //! float reference and the Q16.16 kernel take the same branch, the same
    //! property the trust region's power-of-two halving was chosen for.
    const ik_acc_t num = (ik_acc_t)(b * f) - (ik_acc_t)(c * e);
    const ik_acc_t den = (ik_acc_t)(a * e) - (ik_acc_t)(b * b);
    const ik_acc_t an = (num < (ik_acc_t)0) ? (ik_acc_t)(-num) : num;
    const ik_acc_t ad = (den < (ik_acc_t)0) ? (ik_acc_t)(-den) : den;
    const bool same_sign = ((num < (ik_acc_t)0) == (den < (ik_acc_t)0));

    ik_real_t s;
    if (ad == (ik_acc_t)0) {
        //! Exactly parallel AND the numerator degenerate with it: any point
        //! serves, and the t fixups below pick a valid one.
        s = (ik_real_t)0;
    } else if (an >= ad) {
        //! |s| >= 1: the closest point on the infinite line lies beyond the
        //! segment, so the clamp is the answer and the division is skipped.
        s = same_sign ? (ik_real_t)1 : (ik_real_t)0;
    } else {
        //! |s| < 1, so this cannot overflow. Division, not recip(): a
        //! reciprocal of an unbounded denominator is the thing being avoided.
        s = clamp01((ik_real_t)(num / den));
    }

    ik_real_t t = (ik_real_t)(
        (ik_real_t)((ik_acc_t)(b * s) + (ik_acc_t)f) * inv_e);

    //! ---- t out of range: pin it and re-solve for s ----
    //! Both fixups are computed unconditionally and selected, so the two
    //! extra products cost the same on every input.
    const ik_real_t s_tlo = clamp01((ik_real_t)((ik_real_t)(-c) * inv_a));
    const ik_real_t s_thi = clamp01((ik_real_t)(
        (ik_real_t)((ik_acc_t)b - (ik_acc_t)c) * inv_a));
    const bool t_lo = (t < (ik_real_t)0);
    const bool t_hi = (t > (ik_real_t)1);
    s = t_lo ? s_tlo : (t_hi ? s_thi : s);
    t = t_lo ? (ik_real_t)0 : (t_hi ? (ik_real_t)1 : t);

    //! ---- degenerate segments ----
    //! A zero-length capsule is not emitted by gen_geometry.py, so these
    //! cannot fire for the generated table. They are kept because seg_seg()
    //! is also called directly by the testbench on degenerate inputs, and
    //! because a future robot's table may contain one.
    //
    //! NOT IK_PIVOT_EPS. That constant is 1e-4 and this arm's shortest
    //! capsule, the a3 elbow offset, has a = 4.12e-4 - within 4x of being
    //! declared degenerate and collapsed to a point. The test wanted here is
    //! "exactly zero", and 1e-9 says that in both builds: it is a genuine
    //! 1e-9 in the float reference, and quantises to exactly 0 in Q16.16,
    //! where a zero-length segment is the only thing that can reach it.
    const bool a_deg = (a <= (ik_real_t)1.0e-9);
    const bool e_deg = (e <= (ik_real_t)1.0e-9);
    if (a_deg) {
        s = (ik_real_t)0;
        t = e_deg ? (ik_real_t)0 : clamp01((ik_real_t)(f * inv_e));
    } else if (e_deg) {
        t = (ik_real_t)0;
        s = clamp01((ik_real_t)((ik_real_t)(-c) * inv_a));
    }

    ik_acc_t dsq = (ik_acc_t)0;
SS_OUT:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        wa[i] = (ik_real_t)((ik_acc_t)p0[i] + (ik_acc_t)(s * d1[i]));
        wb[i] = (ik_real_t)((ik_acc_t)q0[i] + (ik_acc_t)(t * d2[i]));
        ik_real_t dd = (ik_real_t)((ik_acc_t)wa[i] - (ik_acc_t)wb[i]);
        dsq += (ik_acc_t)(dd * dd);
    }
    //! Accumulator, not storage: a 1 mm separation is 65 Q16.16 LSBs but its
    //! SQUARE is 4e-3 of an LSB and quantises to zero. Same trap ik_dls
    //! records for tol_sq. sqrt_acc() takes the Q32.32 value directly.
    return ikm::sqrt_acc(dsq);
}

//! ------------------------------------------------------------------
void coll::min_distance(const ik_real_t A[IK_CAP_COUNT][3],
                        const ik_real_t B[IK_CAP_COUNT][3], ik_real_t* d_min,
                        int* pair, ik_real_t wa[3], ik_real_t wb[3]) {
#pragma HLS INLINE off

    //! Starts above any reachable clearance rather than at a sentinel: the
    //! arm's own extent bounds this, and a plain comparison beats carrying
    //! a validity flag through the reduction.
    ik_real_t best = (ik_real_t)1000;
    int best_pair = 0;
    ik_real_t bwa[3], bwb[3];
#pragma HLS ARRAY_PARTITION variable = bwa complete dim = 1
#pragma HLS ARRAY_PARTITION variable = bwb complete dim = 1

CM_INIT:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        bwa[i] = (ik_real_t)0;
        bwb[i] = (ik_real_t)0;
    }

//! Compile-time trip count from the generated pair table, and one shared
//! seg_seg() call site so HLS instantiates the segment engine once rather
//! than IK_COLL_PAIR_COUNT times - same reason ik_dls's MULT loops over a
//! single mm::multiply().
CM_PAIRS:
    for (int k = 0; k < IK_COLL_PAIR_COUNT; k++) {
#pragma HLS PIPELINE off
        const int ia = IK_COLL_PAIRS[k][0];
        const int ib = IK_COLL_PAIRS[k][1];

        ik_real_t cwa[3], cwb[3];
        ik_real_t sd = coll::seg_seg(
            A[ia], B[ia], A[ib], B[ib], (ik_real_t)IK_CAP_LEN_SQ[ia],
            (ik_real_t)IK_CAP_INV_LEN_SQ[ia], (ik_real_t)IK_CAP_LEN_SQ[ib],
            (ik_real_t)IK_CAP_INV_LEN_SQ[ib], cwa, cwb);

        //! Signed: negative means the capsules interpenetrate by that much,
        //! and the amount is what a repulsion term needs. Not clamped.
        ik_real_t d = (ik_real_t)((ik_acc_t)sd -
                                  (ik_acc_t)((ik_real_t)IK_CAP_RADIUS[ia] +
                                             (ik_real_t)IK_CAP_RADIUS[ib]));

        const bool better = (d < best);
        best = better ? d : best;
        best_pair = better ? k : best_pair;
    CM_KEEP:
        for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
            bwa[i] = better ? cwa[i] : bwa[i];
            bwb[i] = better ? cwb[i] : bwb[i];
        }
    }

    *d_min = best;
    *pair = best_pair;
CM_OUT:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        wa[i] = bwa[i];
        wb[i] = bwb[i];
    }
}

//! ------------------------------------------------------------------
//! Standalone IP wrapper
//! ------------------------------------------------------------------
extern "C" void coll_dist_kernel(const ik_word_t q[IK_DOF], ik_word_t* d_min,
                                 int* pair, ik_word_t wa[3], ik_word_t wb[3],
                                 int* status) {
#pragma HLS INTERFACE s_axilite port = q bundle = CTRL
#pragma HLS INTERFACE s_axilite port = d_min bundle = CTRL
#pragma HLS INTERFACE s_axilite port = pair bundle = CTRL
#pragma HLS INTERFACE s_axilite port = wa bundle = CTRL
#pragma HLS INTERFACE s_axilite port = wb bundle = CTRL
#pragma HLS INTERFACE s_axilite port = status bundle = CTRL
#pragma HLS INTERFACE s_axilite port = return bundle = CTRL

    ik_real_t qr[IK_DOF];
CD_IN:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS UNROLL
        qr[i] = ik_from_word(q[i]);
    }

    ik_real_t A[IK_CAP_COUNT][3], B[IK_CAP_COUNT][3];
    coll::endpoints(qr, A, B);

    ik_real_t d, pa[3], pb[3];
    int k = 0;
    coll::min_distance(A, B, &d, &k, pa, pb);

CD_OUT:
    for (int i = 0; i < 3; i++) {
#pragma HLS UNROLL
        wa[i] = ik_to_word(pa[i]);
        wb[i] = ik_to_word(pb[i]);
    }
    *d_min = ik_to_word(d);
    *pair = k;
    *status = (d < (ik_real_t)0) ? IK_ERR_COLLISION : IK_OK;
}
