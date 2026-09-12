#ifndef IK_TYPES_HPP
#define IK_TYPES_HPP

//
//! Numeric types shared by every kernel.
//
//! Format choice: model/validate.py sweeps Q16.16 / Q14.18 / Q12.20 / Q14.26
//! / Q16.32 through the quantised DLS loop; all converge identically
//! (300/300, median 4 iters) - the limit is the damping factor, not word
//! length, so the cheapest format wins.
//
//! -DIK_USE_FLOAT compiles identical algorithms in float, for attributing
//! quantisation error during csim. Never synthesise the float build for PL.
//

#include <stdint.h>

#define IK_FRAC_BITS 16
#define IK_FRAC_SCALE 65536.0

#ifdef IK_USE_FLOAT

//
//! Deliberately excludes <ap_fixed.h> so the testbench suite builds with a
//! stock host compiler, no Vitis install needed. See hls/tb/Makefile.
//
typedef double ik_real_t;
typedef double ik_acc_t;
typedef double ik_work_t;

#else

#include <ap_fixed.h>
#include <ap_int.h>

//
//! Storage / IO type: Q16.16.  AP_RND matches round-to-nearest in
//! model/ik_model.py:to_fixed(); AP_SAT stops a pathological pose from
//! wrapping a joint angle to the opposite sign - dangerous, not just
//! inaccurate.
//
typedef ap_fixed<32, 16, AP_RND, AP_SAT> ik_real_t;

//
//! MAC accumulator. A Q16.16 x Q16.16 product is exactly Q32.32, and the
//! longest dot product here is 6 terms, so Q32.32 holds it exactly - never
//! needs to round or saturate internally. Rounding happens once, on exit,
//! via the (ik_real_t) cast (which carries AP_RND/AP_SAT).
//
//! AP_TRN/AP_WRAP here (not AP_RND/AP_SAT) avoids a rounding adder and dead
//! saturating comparator on every acc += step, on the widest type in the
//! design - safe because the value is provably always exact and in range.
//
typedef ap_fixed<64, 32, AP_TRN, AP_WRAP> ik_acc_t;

//
//! Internal working storage: same Q16.16 as ik_real_t, but truncating -
//! saves a rounding adder per store (LUTs are the binding resource here;
//! ik_dls_kernel sits at 96%). Bias cost is <=1 LSB, well below testbench
//! tolerances.
//
//! Saturation is KEPT (unlike ik_acc_t): factorisation intermediates can
//! genuinely grow, and wrapping would sign-flip a result rather than just
//! degrade it.
//
//! Values crossing an IP boundary stay ik_real_t (bit-compared against
//! model/ik_model.py); use this only for values that stay inside a kernel.
//
typedef ap_fixed<32, 16, AP_TRN, AP_SAT> ik_work_t;

#endif

//! Raw 32-bit AXI word carrying a Q16.16 value.
typedef int32_t ik_word_t;

//! Bit-exact reinterpretation in the fixed build, so a value written by
//! the PS driver arrives unchanged.
inline ik_real_t ik_from_word(ik_word_t w) {
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    return (ik_real_t)w / IK_FRAC_SCALE;
#else
    ik_real_t r;
    r.range(31, 0) = (ap_uint<32>)(uint32_t)w;
    return r;
#endif
}

inline ik_word_t ik_to_word(ik_real_t v) {
#pragma HLS INLINE
#ifdef IK_USE_FLOAT
    double f = (double)v * IK_FRAC_SCALE;
    if (f > 2147483647.0)
        f = 2147483647.0;
    if (f < -2147483648.0)
        f = -2147483648.0;
    return (ik_word_t)(int64_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5));
#else
    return (ik_word_t)(int32_t)(ap_int<32>)v.range(31, 0);
#endif
}

#endif  //! IK_TYPES_HPP
