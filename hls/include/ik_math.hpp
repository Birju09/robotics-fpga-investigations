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
/* Reciprocal by table-seeded Newton-Raphson.                          */
/* ------------------------------------------------------------------ */
/*
 * spd::solve() needs six reciprocals per call, and an ap_fixed divide
 * synthesises a sequential restoring divider of roughly one cycle per result
 * bit - about 35 each, which made division the single largest term in the
 * solver once the Gauss-Jordan elimination was gone.
 *
 * This is integer-only and identical in both builds, so hls/tb/tb_spd.cpp can
 * check it against exact reciprocals even though the host regression runs in
 * double and never exercises the fixed-point path otherwise.  That mattered:
 * a silently wrong reciprocal corrupts every DLS solve, and there would
 * otherwise be no way to catch it before hardware.
 *
 * Swept over every raw input from 7 (just above the IK_PIVOT_EPS floor) to
 * 20.0, no result differs from the exact reciprocal by more than one raw
 * LSB.  It may still differ from the built-in divider by that LSB, which
 * rounds where this truncates.
 *
 * Method: normalise d = u * 2^(e+1) with u in [0.5, 1), seed 1/u from a
 * 128-entry table (7 bits, since the leading bit of u is always set), then
 * two Newton steps r <- r*(2 - u*r).  Quadratic convergence takes 8 bits to
 * 32.  Q16.16 in, Q16.16 out: the result raw value is 2^32/d_raw.
 */
static const uint32_t RECIP_SEED[128] = {
    2139127680u, 2122609320u, 2106344115u, 2090326289u,
    2074550241u, 2059010539u, 2043701910u, 2028619239u,
    2013757560u, 1999112051u, 1984678028u, 1970450946u,
    1956426384u, 1942600049u, 1928967768u, 1915525484u,
    1902269252u, 1889195237u, 1876299706u, 1863579030u,
    1851029676u, 1838648207u, 1826431275u, 1814375623u,
    1802478078u, 1790735550u, 1779145029u, 1767703582u,
    1756408351u, 1745256552u, 1734245470u, 1723372457u,
    1712634934u, 1702030384u, 1691556350u, 1681210440u,
    1670990316u, 1660893698u, 1650918360u, 1641062131u,
    1631322890u, 1621698566u, 1612187138u, 1602786629u,
    1593495113u, 1584310703u, 1575231558u, 1566255880u,
    1557381909u, 1548607926u, 1539932252u, 1531353242u,
    1522869291u, 1514478826u, 1506180312u, 1497972245u,
    1489853154u, 1481821601u, 1473876177u, 1466015504u,
    1458238233u, 1450543045u, 1442928645u, 1435393770u,
    1427937179u, 1420557659u, 1413254020u, 1406025099u,
    1398869755u, 1391786871u, 1384775350u, 1377834120u,
    1370962129u, 1364158347u, 1357421763u, 1350751385u,
    1344146244u, 1337605387u, 1331127879u, 1324712805u,
    1318359266u, 1312066382u, 1305833287u, 1299659134u,
    1293543092u, 1287484342u, 1281482084u, 1275535531u,
    1269643912u, 1263806469u, 1258022457u, 1252291148u,
    1246611823u, 1240983779u, 1235406323u, 1229878778u,
    1224400476u, 1218970763u, 1213588993u, 1208254536u,
    1202966770u, 1197725085u, 1192528880u, 1187377568u,
    1182270568u, 1177207310u, 1172187236u, 1167209796u,
    1162274448u, 1157380661u, 1152527912u, 1147715687u,
    1142943480u, 1138210795u, 1133517142u, 1128862041u,
    1124245018u, 1119665609u, 1115123355u, 1110617806u,
    1106148519u, 1101715058u, 1097316994u, 1092953904u,
    1088625374u, 1084330994u, 1080070361u, 1075843080u,
};

inline int32_t recip_q16_raw(int32_t d_raw)
{
#pragma HLS INLINE off
    bool     neg = d_raw < 0;
    uint32_t d   = neg ? (uint32_t)(-(int64_t)d_raw) : (uint32_t)d_raw;
    if (d == 0)
        return neg ? (int32_t)0x80000000 : (int32_t)0x7FFFFFFF;

    /* Position of the most significant set bit. */
    int e = 0;
RECIP_MSB:
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        if (d >> i) e = i;
    }

    uint32_t m = d << (31 - e);                 /* Q0.32, bit 31 set    */
    uint64_t r = RECIP_SEED[(m >> 24) & 0x7F];  /* Q2.30, 1/u           */

RECIP_NR:
    for (int k = 0; k < 2; k++) {
#pragma HLS UNROLL
        uint64_t ur = ((uint64_t)m * r) >> 32;          /* u*r,  Q2.30  */
        uint64_t t  = ((uint64_t)1 << 31) - ur;         /* 2-u*r, Q2.30 */
        r = (r * t) >> 30;
    }

    /* 1/d = (1/u) * 2^(31-e), and r is (1/u) scaled by 2^30. */
    int      sh  = e - 1;
    uint64_t res = (sh >= 0) ? (r >> sh) : (r << (-sh));
    if (res > 0x7FFFFFFFULL) res = 0x7FFFFFFFULL;
    return neg ? -(int32_t)res : (int32_t)res;
}

/*
 * Reciprocal of a working-type value.  The float reference build keeps exact
 * division on purpose - it is the reference, and quantising it here would
 * remove the very baseline the fixed-point error is measured against.
 * recip_q16_raw() above is still compiled in that build, and still tested.
 */
inline ik_work_t recip(ik_work_t d)
{
#pragma HLS INLINE off
#if defined(IK_USE_FLOAT) || !defined(IK_FAST_RECIP)
    return (ik_work_t)((ik_real_t)1 / (ik_real_t)d);
#else
    ik_work_t out;
    out.range(31, 0) =
        (ap_uint<32>)(uint32_t)recip_q16_raw((int32_t)(uint32_t)d.range(31, 0));
    return out;
#endif
}

/* ------------------------------------------------------------------ */
/* CORDIC rotation mode -> sin and cos of `angle`.                     */
/* ------------------------------------------------------------------ */
/*
 * The engine itself, always inlined, with ROT left unannotated so the two
 * entry points below can schedule it differently.  One implementation, two
 * schedules: a second copy of the CORDIC written out to get a pipelined
 * variant would be an open invitation for the two to drift apart.
 */
inline void sincos_core(ik_real_t angle, ik_real_t &sin_o, ik_real_t &cos_o)
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
        /* Deliberately rolled - fixed trip count already gives constant
         * latency; UNROLL would spatially replicate the 40-bit adder chain
         * 24x per call site instead of reusing one. */
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

/*
 * One angle at a time, ROT left rolled.
 *
 * INLINE off, unlike the trivial helpers above: this holds the 24-cycle
 * CORDIC engine, and iks::analytic() alone has roughly a dozen call sites
 * into this and atan2_hypot() below.  Inlined, each one would get its own
 * copy of the engine even after rolling its loop; off, every call site
 * shares one synthesised instance and just queues for it.
 */
inline void sincos(ik_real_t angle, ik_real_t &sin_o, ik_real_t &cos_o)
{
#pragma HLS INLINE off
    sincos_core(angle, sin_o, cos_o);
}

/*
 * All IK_DOF angles from one pipelined pass.
 *
 * fk_jacobian() used to take its six sin/cos values one per DH step, inside
 * a chain that cannot start step i+1 until step i has updated the frame - so
 * six CORDIC latencies, roughly 24 cycles each, fell on the critical path in
 * series.  But the angles are all known at entry: nothing about sincos(q[i])
 * depends on step i-1.  Hoisting them here turns 6 x 24 sequential cycles
 * into 6 + depth.
 *
 * PIPELINE II=1 forces HLS to unroll ROT, which is exactly the spatial
 * replication sincos() above is careful to avoid - about 24 stages of the
 * 40-bit adder chain.  That is affordable HERE and only here: it buys a
 * result per cycle across six calls, where a single call site would pay the
 * same area for one.  CORDIC is shift-add, so the cost is LUTs and flops,
 * not DSPs - which matters because DSPs are the binding resource.
 *
 * Deliberately NOT used by fk() or rot03(): those are on ik_analytic's path,
 * which is already at 86% LUT and has no room for a second CORDIC instance.
 * fk_jacobian() is reachable only from ik_dls.
 */
inline void sincos_batch(const ik_real_t th[IK_DOF],
                         ik_real_t s[IK_DOF], ik_real_t c[IK_DOF])
{
#pragma HLS INLINE off
SC_BATCH:
    for (int i = 0; i < IK_DOF; i++) {
#pragma HLS PIPELINE II=1
        ik_real_t si, ci;
        sincos_core(th[i], si, ci);
        s[i] = si;
        c[i] = ci;
    }
}

/* ------------------------------------------------------------------ */
/* CORDIC vectoring mode -> atan2(y, x) over the full circle, and the  */
/* magnitude hypot(x, y) as a free by-product.                         */
/* ------------------------------------------------------------------ */
inline void atan2_hypot(ik_real_t yi, ik_real_t xi,
                        ik_real_t &ang_o, ik_real_t &mag_o)
{
#pragma HLS INLINE off   /* shared CORDIC engine - see sincos() above */
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
    for (int i = 0; i < CORDIC_ITER; i++) {   /* rolled - see sincos() above */
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
#pragma HLS INLINE off   /* shared 24-cycle engine - see sincos() above */
    if (a <= (ik_real_t)0)
        return (ik_real_t)0;

    uint64_t v   = ((uint64_t)(uint32_t)ik_to_word(a)) << IK_FRAC_BITS;
    uint64_t res = 0;
    uint64_t bit = (uint64_t)1 << 46;      /* highest even power of 4 needed */

SQRT:
    for (int i = 0; i < 24; i++) {   /* rolled - see sincos() above */
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
#pragma HLS INLINE off   /* shared 32-cycle engine - see sincos() above */
    if (a <= (ik_acc_t)0)
        return (ik_real_t)0;

    uint64_t v   = acc_raw(a);
    uint64_t res = 0;
    uint64_t bit = (uint64_t)1 << 62;

SQRTA:
    for (int i = 0; i < 32; i++) {   /* rolled - see sincos() above */
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
