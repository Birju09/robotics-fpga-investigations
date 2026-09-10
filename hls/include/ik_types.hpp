#ifndef IK_TYPES_HPP
#define IK_TYPES_HPP

/*
 * Numeric types shared by every kernel.
 *
 * Format selection was not arbitrary: model/validate.py sweeps Q16.16 / Q14.18 /
 * Q12.20 / Q14.26 / Q16.32 through the full quantised DLS loop and they all
 * converge identically (300/300, median 4 iterations).  The limit is the
 * damping factor, not the word length, so the cheapest format wins.
 *
 * Build with -DIK_USE_FLOAT to compile the identical algorithms in float; that
 * is how quantisation error is attributed during csim.  Do not synthesise the
 * float build for the PL - it is a reference, not a deliverable.
 */

#include <stdint.h>

#define IK_FRAC_BITS 16
#define IK_FRAC_SCALE 65536.0

#ifdef IK_USE_FLOAT

/*
 * Float reference build.  Deliberately does NOT include <ap_fixed.h>, so this
 * configuration compiles with a stock host compiler and the whole testbench
 * suite can be run without a Vitis installation.  See hls/tb/Makefile.
 */
typedef double ik_real_t;
typedef double ik_acc_t;

#else

#include <ap_fixed.h>
#include <ap_int.h>

/*
 * Storage / IO type: Q16.16.  AP_RND matches the round-to-nearest used by
 * model/ik_model.py:to_fixed(); AP_SAT keeps a pathological pose from wrapping
 * a joint angle to the opposite sign, which is the one failure mode that would
 * be dangerous on real hardware rather than merely inaccurate.
 */
typedef ap_fixed<32, 16, AP_RND, AP_SAT> ik_real_t;

/*
 * MAC accumulator.  A product of two Q16.16 values is exactly Q32.32, and the
 * longest dot product here is 6 terms, so Q32.32 accumulates the whole thing
 * without a single rounding step - by construction, this type never actually
 * needs to round or saturate internally.  Rounding happens once, on the way
 * out, through the (ik_real_t) cast, which does carry AP_RND/AP_SAT.
 *
 * It was nonetheless declared AP_RND/AP_SAT itself, which means every one of
 * the many acc += ... accumulation steps throughout the kernels synthesises a
 * rounding adder and a saturating comparator that can never fire - dead
 * hardware, paid for at every call site, on the widest type in the design.
 * AP_TRN/AP_WRAP cost nothing extra when the value is already exact and in
 * range, which it provably is here; the true precision/safety-critical
 * requirements live entirely in ik_real_t above and are unaffected by this.
 */
typedef ap_fixed<64, 32, AP_TRN, AP_WRAP> ik_acc_t;

#endif

/* Raw 32-bit AXI word carrying a Q16.16 value. */
typedef int32_t ik_word_t;

/* ------------------------------------------------------------------ */
/* Word <-> real conversion.  Bit-exact reinterpretation in the fixed  */
/* build, so a value written by the PS driver arrives unchanged.       */
/* ------------------------------------------------------------------ */
inline ik_real_t ik_from_word(ik_word_t w)
{
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    return (ik_real_t)w / IK_FRAC_SCALE;
#else
    ik_real_t r;
    r.range(31, 0) = (ap_uint<32>)(uint32_t)w;
    return r;
#endif
}

inline ik_word_t ik_to_word(ik_real_t v)
{
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    double f = (double)v * IK_FRAC_SCALE;
    if (f > 2147483647.0)  f = 2147483647.0;
    if (f < -2147483648.0) f = -2147483648.0;
    return (ik_word_t)(int64_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5));
#else
    return (ik_word_t)(int32_t)(ap_int<32>)v.range(31, 0);
#endif
}

#endif /* IK_TYPES_HPP */
