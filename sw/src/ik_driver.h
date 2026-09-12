#ifndef IK_DRIVER_H
#define IK_DRIVER_H

#include <stdint.h>

//
//! Vitis compiles the whole app with g++ (even .c files), so declarations here
//! get C++ name mangling unless guarded - ik_sw_ref.cpp expects extern "C".
//
#ifdef __cplusplus
extern "C" {
#endif

//
//! Bare-metal driver for the four HLS kernels.
//
//! Register offsets are assigned by HLS from argument order/size and must be
//! rechecked against the generated headers after any argument-list change:
//! find hls/build/<kernel> -path '*impl/ip/drivers*' -name 'x<kernel>_hw.h'
//! verify_reg_map() in main.c catches a stale map via write/read-back before
//! it produces plausible-looking wrong answers.
//
//! The AP_CTRL block is common to all of them:
//! bit 0  ap_start   (write 1 to launch; self-clearing)
//! bit 1  ap_done    (read)
//! bit 2  ap_idle    (read)
//! bit 3  ap_ready   (read)
//! bit 7  auto_restart
//
//! Fixed point: every data word is signed Q16.16: value = raw / 65536.0.
//

#define IK_FRAC_BITS 16
#define IK_FRAC_SCALE 65536.0f

#define IK_DOF 6
#define IK_MAT_MAX 6

//! ---- ap_ctrl ----
#define IK_ADDR_AP_CTRL 0x00
#define IK_AP_START 0x1
#define IK_AP_DONE 0x2
#define IK_AP_IDLE 0x4
#define IK_AP_READY 0x8

//! ---- status codes, mirroring hls/include/ik_config.hpp ----
#define IK_OK 0
#define IK_ERR_UNREACH 1
#define IK_ERR_NO_CONV 2
#define IK_ERR_SINGULAR 3
#define IK_ERR_BADDIM 4

//! ---- analytic branch selection ----
#define IK_CFG_SHOULDER 0x1
#define IK_CFG_ELBOW 0x2
#define IK_CFG_WRIST 0x4

typedef int32_t ik_word_t;

typedef struct {
    uintptr_t base;
} ik_dev_t;

//! Q16.16 conversion.
static inline ik_word_t ik_f2q(float v) {
    float s = v * IK_FRAC_SCALE;
    return (ik_word_t)(s >= 0.0f ? (s + 0.5f) : (s - 0.5f));
}

static inline float ik_q2f(ik_word_t w) {
    return (float)w / IK_FRAC_SCALE;
}

//! ---------------- kernel entry points ----------------
//! Each returns the kernel's status word; *cycles receives elapsed PL-clock
//! cycles (PS global timer) when non-NULL.

int ik_analytic_solve(ik_dev_t* dev, const float pose[6], int cfg,
                      float q_out[6], uint32_t* cycles);

//! step_max bounds |dq|_inf per DLS iteration, in radians; 0 disables the
//! clamp.  See the trust region block in hls/include/ik_config.hpp.
int ik_dls_solve(ik_dev_t* dev, const float pose[6], const float q_seed[6],
                 float lambda, float tol, int max_iter, float step_max,
                 float q_out[6], int* iters, float* resid, uint32_t* cycles);

int ik_matmul(ik_dev_t* dev, int m, int k, int n, int ta, int tb,
              const float* A, const float* B, float* C, uint32_t* cycles);

int ik_matinv(ik_dev_t* dev, int n, const float* A, float* Ainv,
              uint32_t* cycles);

//! ---------------- software reference (runs on the A9, float, for PS-vs-PL comparison) ----------------
int ik_analytic_solve_sw(const float pose[6], int cfg, float q_out[6]);
int ik_dls_solve_sw(const float pose[6], const float q_seed[6], float lambda,
                    float tol, int max_iter, float step_max, float q_out[6],
                    int* iters, float* resid);

//! ---------------- cycle counter ----------------
void ik_timer_init(void);
uint64_t ik_timer_read(void);  //! free-running, CPU_3x2x ticks
uint32_t ik_timer_hz(void);

#ifdef __cplusplus
}
#endif

#endif  //! IK_DRIVER_H
