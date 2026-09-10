#ifndef IK_MATH_HPP
#define IK_MATH_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Elementary functions, implemented directly rather than pulled from
 * hls_math.h.
 *
 * Two reasons, both specific to this investigation:
 *   1. Determinism.  Every routine below is a fixed trip-count loop, so the
 *      latency reported by C-synthesis is the latency you get for any input.
 *      That is the whole premise of comparing a closed-form solver against an
 *      iterative one - the closed-form path must have no data-dependent cost.
 *   2. ap_fixed support in the vendor math library varies by word length and
 *      release; a CORDIC written here behaves identically in csim, cosim and
 *      the float reference build.
 *
 * All angles are radians.
 */

namespace ikm {

/* Internal CORDIC datapath.  Range +-128 covers every intermediate in this
 * design (max |position| = 0.81 m, |rotation entry| = 1, CORDIC gain 1.65),
 * and 32 fractional bits keep the accumulated shift-truncation error near
 * 1e-9 - about four decades below the Q16.16 output LSB. */
#ifdef IK_USE_FLOAT
typedef double cordic_t;
#else
typedef ap_fixed<40, 8> cordic_t;
#endif

static const int CORDIC_ITER = 24;

/* 1 / prod(sqrt(1 + 4^-i)) for i in [0, CORDIC_ITER) */
static const double CORDIC_K_INV = 0.60725293500888267;

static const double CORDIC_ATAN_TAB[CORDIC_ITER] = {
    0.78539816339744828,    0.46364760900080615,
    0.24497866312686414,    0.12435499454676144,
    0.06241880999595735,    0.031239833430268277,
    0.015623728620476831,   0.0078123410601011111,
    0.0039062301319669718,  0.0019531225164788188,
    0.00097656218955931946, 0.00048828121119489829,
    0.00024414062014936177, 0.00012207031189367021,
    6.1035156174208773e-05, 3.0517578115526096e-05,
    1.5258789061315762e-05, 7.62939453110197e-06,
    3.8146972656064961e-06, 1.907348632810187e-06,
    9.5367431640596084e-07, 4.7683715820308884e-07,
    2.3841857910155797e-07, 1.1920928955078068e-07,
};

/* Exact division by 2^i.  An arithmetic right shift on ap_fixed, a
 * power-of-two scale in the float reference build - both are exact, so the
 * two builds follow bit-identical CORDIC trajectories. */
inline cordic_t cshift(cordic_t v, int i)
{
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    return v / (double)((uint64_t)1 << i);
#else
    return v >> i;
#endif
}

/* ------------------------------------------------------------------ */
/* Wrap an angle into [-pi, pi].                                       */
/* Four conditional folds, unrolled - covers |angle| < 9*pi with no    */
/* data-dependent trip count.                                          */
/* ------------------------------------------------------------------ */
inline ik_real_t wrap_pi(ik_real_t a)
{
#pragma HLS INLINE
    cordic_t z = (cordic_t)a;
WRAP:
    for (int i = 0; i < 4; i++) {
#pragma HLS UNROLL
        if (z > (cordic_t)IK_PI)       z -= (cordic_t)IK_TWO_PI;
        else if (z < -(cordic_t)IK_PI) z += (cordic_t)IK_TWO_PI;
    }
    return (ik_real_t)z;
}

/* ------------------------------------------------------------------ */
/* CORDIC rotation mode -> sin and cos of `angle`.                     */
/* ------------------------------------------------------------------ */
inline void sincos(ik_real_t angle, ik_real_t &sin_o, ik_real_t &cos_o)
{
#pragma HLS INLINE
    cordic_t z = (cordic_t)angle;
    bool negate = false;

    /* Fold to [-pi, pi] ... */
FOLD2PI:
    for (int i = 0; i < 4; i++) {
#pragma HLS UNROLL
        if (z > (cordic_t)IK_PI)       z -= (cordic_t)IK_TWO_PI;
        else if (z < -(cordic_t)IK_PI) z += (cordic_t)IK_TWO_PI;
    }
    /* ... then to [-pi/2, pi/2], which is inside the +-1.7433 rad
     * convergence range of the circular CORDIC. */
    if (z > (cordic_t)IK_HALF_PI) {
        z -= (cordic_t)IK_PI;
        negate = true;
    } else if (z < -(cordic_t)IK_HALF_PI) {
        z += (cordic_t)IK_PI;
        negate = true;
    }

    cordic_t x = (cordic_t)CORDIC_K_INV;   /* pre-scaled so gain cancels */
    cordic_t y = 0;

ROT:
    for (int i = 0; i < CORDIC_ITER; i++) {
#pragma HLS UNROLL
        cordic_t xs = cshift(x, i);
        cordic_t ys = cshift(y, i);
        if (z < 0) {            /* d = -1 */
            x += ys;
            y -= xs;
            z += (cordic_t)CORDIC_ATAN_TAB[i];
        } else {                /* d = +1 */
            x -= ys;
            y += xs;
            z -= (cordic_t)CORDIC_ATAN_TAB[i];
        }
    }

    cos_o = negate ? (ik_real_t)(-x) : (ik_real_t)x;
    sin_o = negate ? (ik_real_t)(-y) : (ik_real_t)y;
}

/* ------------------------------------------------------------------ */
/* CORDIC vectoring mode -> atan2(y, x) over the full circle, and the  */
/* magnitude hypot(x, y) as a free by-product.                         */
/* ------------------------------------------------------------------ */
inline void atan2_hypot(ik_real_t yi, ik_real_t xi,
                        ik_real_t &ang_o, ik_real_t &mag_o)
{
#pragma HLS INLINE
    cordic_t x = (cordic_t)xi;
    cordic_t y = (cordic_t)yi;
    cordic_t z = 0;

    /* Vectoring converges only for x > 0; reflect and correct afterwards. */
    bool reflected = false;
    bool y_nonneg  = (yi >= 0);
    if (x < 0) {
        x = -x;
        y = -y;
        reflected = true;
    }

VEC:
    for (int i = 0; i < CORDIC_ITER; i++) {
#pragma HLS UNROLL
        cordic_t xs = cshift(x, i);
        cordic_t ys = cshift(y, i);
        if (y < 0) {            /* d = +1 */
            x -= ys;
            y += xs;
            z -= (cordic_t)CORDIC_ATAN_TAB[i];
        } else {                /* d = -1 */
            x += ys;
            y -= xs;
            z += (cordic_t)CORDIC_ATAN_TAB[i];
        }
    }

    if (reflected)
        z = y_nonneg ? (z + (cordic_t)IK_PI) : (z - (cordic_t)IK_PI);

    ang_o = (ik_real_t)z;
    mag_o = (ik_real_t)(x * (cordic_t)CORDIC_K_INV);
}

inline ik_real_t atan2(ik_real_t y, ik_real_t x)
{
#pragma HLS INLINE
    ik_real_t a, m;
    atan2_hypot(y, x, a, m);
    return a;
}

inline ik_real_t hypot(ik_real_t x, ik_real_t y)
{
#pragma HLS INLINE
    ik_real_t a, m;
    atan2_hypot(y, x, a, m);
    return m;
}

/* ------------------------------------------------------------------ */
/* Square root by restoring shift-subtract on the raw Q16.16 word.     */
/*                                                                     */
/* sqrt(a) in Q16.16 == isqrt(raw(a) << 16), so the whole thing is     */
/* integer arithmetic and is bit-identical in the float and fixed      */
/* builds.  Fixed 24 iterations, no early exit.                        */
/* ------------------------------------------------------------------ */
inline ik_real_t sqrt(ik_real_t a)
{
#pragma HLS INLINE
    if (a <= (ik_real_t)0)
        return (ik_real_t)0;

    uint64_t v   = ((uint64_t)(uint32_t)ik_to_word(a)) << IK_FRAC_BITS;
    uint64_t res = 0;
    uint64_t bit = (uint64_t)1 << 46;      /* highest even power of 4 needed */

SQRT:
    for (int i = 0; i < 24; i++) {
#pragma HLS UNROLL
        uint64_t t = res + bit;
        res >>= 1;
        if (v >= t) {
            v -= t;
            res += bit;
        }
        bit >>= 2;
    }
    /* `res` is floor(sqrt(.)); the residue test converts it to
     * round-to-nearest so the error is +-0.5 LSB instead of a systematic
     * -1 LSB bias that would accumulate through the DLS iteration. */
    if (v > res)
        res += 1;
    return ik_from_word((ik_word_t)(uint32_t)res);
}

/* ------------------------------------------------------------------ */
/* Square root of an accumulator-typed value.                          */
/*                                                                     */
/* Needed because the DLS convergence test cannot be done in Q16.16:   */
/* a tolerance of 1e-3 squares to 1e-6, which is two decades below the */
/* Q16.16 LSB and would quantise to zero, making the solver report     */
/* convergence on its first iteration.  The residual is therefore      */
/* accumulated in Q32.32 and compared there.                           */
/*                                                                     */
/* Conveniently sqrt(v * 2^-32) = sqrt(v) * 2^-16, so an integer sqrt  */
/* of the raw Q32.32 word lands directly on the raw Q16.16 word.       */
/* ------------------------------------------------------------------ */
inline uint64_t acc_raw(ik_acc_t a)
{
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    double s = (double)a * 4294967296.0;         /* 2^32 */
    if (s < 0.0) s = 0.0;
    if (s > 18446744073709549568.0) s = 18446744073709549568.0;
    return (uint64_t)s;
#else
    return (uint64_t)(ap_uint<64>)a.range(63, 0);
#endif
}

inline ik_real_t sqrt_acc(ik_acc_t a)
{
#pragma HLS INLINE
    if (a <= (ik_acc_t)0)
        return (ik_real_t)0;

    uint64_t v   = acc_raw(a);
    uint64_t res = 0;
    uint64_t bit = (uint64_t)1 << 62;

SQRTA:
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        uint64_t t = res + bit;
        res >>= 1;
        if (v >= t) {
            v -= t;
            res += bit;
        }
        bit >>= 2;
    }
    if (v > res)
        res += 1;
    if (res > 0x7FFFFFFFULL)
        res = 0x7FFFFFFFULL;
    return ik_from_word((ik_word_t)(uint32_t)res);
}

/* Clamp helper used by the analytic solver's law-of-cosines term. */
inline ik_real_t clamp1(ik_real_t v)
{
#pragma HLS INLINE
    if (v > (ik_real_t)1)  return (ik_real_t)1;
    if (v < (ik_real_t)-1) return (ik_real_t)-1;
    return v;
}

} /* namespace ikm */

#endif /* IK_MATH_HPP */
