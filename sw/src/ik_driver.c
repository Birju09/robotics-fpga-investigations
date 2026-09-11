#include "ik_driver.h"
#include "ik_regmap.h"

#include "xil_io.h"
#include "xparameters.h"

/*
 * XTime_GetTime/COUNTS_PER_SECOND moved house.  They were in the standalone
 * BSP's xtime_l.h; under the 2025.2 system-device-tree flow the BSP builds
 * libxiltimer instead and the declarations come from xiltimer.h, with
 * xtime_l.h absent entirely.  Both spellings are still in the field, so pick
 * whichever the BSP actually shipped rather than hard-coding one and failing
 * on the other.
 */
#if defined(__has_include)
#  if __has_include("xtime_l.h")
#    include "xtime_l.h"
#    define IK_HAVE_XTIME 1
#  elif __has_include("xiltimer.h")
#    include "xiltimer.h"
#    define IK_HAVE_XTIME 1
#    define IK_HAVE_XILTIMER 1
#  endif
#else
#  include "xtime_l.h"
#  define IK_HAVE_XTIME 1
#endif

#ifndef IK_HAVE_XTIME
#error "Neither xtime_l.h nor xiltimer.h is in the BSP; no time base for the measurements."
#endif

/*
 * Kernel invocation is the same three steps everywhere: stage the arguments,
 * pulse ap_start, poll ap_done.  Polling rather than interrupting is
 * deliberate - an interrupt would add scheduler latency to the very number
 * this project is trying to measure.
 *
 * *cycles, where a caller asks for it, spans argument writes through result
 * reads - the whole PS-side transaction, not just the ap_start/ap_done
 * window - because that AXI4-Lite traffic is part of what a real caller
 * waits on. See ik_wait_done()'s comment for why this moved out of it.
 */

#define RD(d, off)      Xil_In32((d)->base + (off))
#define WR(d, off, val) Xil_Out32((d)->base + (off), (uint32_t)(val))

/* ---------------- cycle counter ---------------- */
/*
 * The Cortex-A9 global timer runs at CPU_3x2x, i.e. half the CPU clock
 * (333 MHz on a -1 speed grade Zynq-7020 at the usual 667 MHz).  XTime_GetTime
 * returns its 64-bit count.  Resolution is ~3 ns, which is fine against
 * kernel latencies in the microsecond range.
 *
 * libxiltimer's xiltimer.c already runs XilSleepTimer_Init(&TimerInst) once
 * from a __attribute__((constructor)) function, which is supposed to fire
 * before main() with no help from the application.  On this BSP it evidently
 * doesn't (or doesn't finish arming the counter hardware in time): without
 * the explicit re-init and the one-tick nudge below, every XTime_GetTime()
 * call read back the same value and every measured latency in this harness
 * came out as exactly zero. The usleep(1) exercises the sleep-timer path
 * once, which is what actually starts the counter.  None of this applies to
 * the older xtime_l.h BSP, which lazily arms the counter on first read.
 */
void ik_timer_init(void)
{
#if defined(IK_HAVE_XILTIMER)
    (void)XilSleepTimer_Init(&TimerInst);
    usleep(1);
#endif
}

uint64_t ik_timer_read(void)
{
    XTime t;
    XTime_GetTime(&t);
    return (uint64_t)t;
}

uint32_t ik_timer_hz(void)
{
    return (uint32_t)COUNTS_PER_SECOND;
}

/* ---------------- low-level helpers ---------------- */

static void ik_write_block(ik_dev_t *dev, uintptr_t base,
                           const float *src, int n)
{
    for (int i = 0; i < n; i++)
        WR(dev, base + 4 * i, ik_f2q(src[i]));
}

static void ik_read_block(ik_dev_t *dev, uintptr_t base, float *dst, int n)
{
    for (int i = 0; i < n; i++)
        dst[i] = ik_q2f((ik_word_t)RD(dev, base + 4 * i));
}

/*
 * Pulse ap_start and poll until ap_done, or until a 1-second bound trips -
 * the kernels finish in microseconds, so that bound only ever fires if a
 * register is mis-mapped or the PL is unprogrammed, not in normal operation.
 *
 * This used to also be where *cycles was measured, timing only ap_start to
 * ap_done.  That excluded the argument writes before it and the result reads
 * after it - 14 of the roughly 20 AXI4-Lite transactions a solve costs - so
 * the reported latency undercounted what a real caller actually waits on.
 * The callers below now time their own argument-write-to-result-read span
 * and this function only launches and waits.
 */
static void ik_wait_done(ik_dev_t *dev, int *timed_out)
{
    const uint64_t limit = (uint64_t)ik_timer_hz();   /* 1 second */
    uint64_t t0 = ik_timer_read();

    if (timed_out) *timed_out = 0;

    WR(dev, IK_ADDR_AP_CTRL, IK_AP_START);

    for (;;) {
        uint32_t c = RD(dev, IK_ADDR_AP_CTRL);
        if (c & IK_AP_DONE)
            break;
        if ((ik_timer_read() - t0) > limit) {
            if (timed_out) *timed_out = 1;
            break;
        }
    }

    /* ap_done is clear-on-read for the AP_CTRL_HS protocol HLS generates
     * here; the read above already acknowledged it. */
}

/* ---------------- analytic IK ---------------- */

int ik_analytic_solve(ik_dev_t *dev, const float pose[6], int cfg,
                      float q_out[6], uint32_t *cycles)
{
    int to = 0;
    uint64_t t0 = ik_timer_read();

    ik_write_block(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_POSE_BASE, pose, IK_DOF);
    WR(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_CFG_DATA, cfg);

    ik_wait_done(dev, &to);
    if (to) {
        if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
        return -1;
    }

    ik_read_block(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_Q_BASE, q_out, IK_DOF);
    int st = (int)RD(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_STATUS_DATA);
    if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
    return st;
}

/* ---------------- DLS IK ---------------- */

int ik_dls_solve(ik_dev_t *dev, const float pose[6], const float q_seed[6],
                 float lambda, float tol, int max_iter,
                 float q_out[6], int *iters, float *resid, uint32_t *cycles)
{
    int to = 0;
    uint64_t t0 = ik_timer_read();

    ik_write_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_POSE_BASE, pose, IK_DOF);
    ik_write_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_Q_SEED_BASE, q_seed, IK_DOF);
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_LAMBDA_DATA,   ik_f2q(lambda));
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_TOL_DATA,      ik_f2q(tol));
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_MAX_ITER_DATA, max_iter);

    ik_wait_done(dev, &to);
    if (to) {
        if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
        return -1;
    }

    ik_read_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_Q_BASE, q_out, IK_DOF);
    if (iters)
        *iters = (int)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_ITERS_DATA);
    if (resid)
        *resid = ik_q2f((ik_word_t)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_RESID_DATA));
    int st = (int)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_STATUS_DATA);
    if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
    return st;
}

/* ---------------- matrix multiply ---------------- */

int ik_matmul(ik_dev_t *dev, int m, int k, int n, int ta, int tb,
              const float *A, const float *B, float *C, uint32_t *cycles)
{
    int to = 0;
    const int NN = IK_MAT_MAX * IK_MAT_MAX;
    uint64_t t0 = ik_timer_read();

    ik_write_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_A_BASE, A, NN);
    ik_write_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_B_BASE, B, NN);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_M_DATA,  m);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_K_DATA,  k);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_N_DATA,  n);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_TA_DATA, ta);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_TB_DATA, tb);

    ik_wait_done(dev, &to);
    if (to) {
        if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
        return -1;
    }

    ik_read_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_C_BASE, C, NN);
    int st = (int)RD(dev, XMAT_MUL_KERNEL_CTRL_ADDR_STATUS_DATA);
    if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
    return st;
}

/* ---------------- matrix inversion ---------------- */

int ik_matinv(ik_dev_t *dev, int n, const float *A, float *Ainv,
              uint32_t *cycles)
{
    int to = 0;
    const int NN = IK_MAT_MAX * IK_MAT_MAX;
    uint64_t t0 = ik_timer_read();

    ik_write_block(dev, XMAT_INV_KERNEL_CTRL_ADDR_A_BASE, A, NN);
    WR(dev, XMAT_INV_KERNEL_CTRL_ADDR_N_DATA, n);

    ik_wait_done(dev, &to);
    if (to) {
        if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
        return -1;
    }

    ik_read_block(dev, XMAT_INV_KERNEL_CTRL_ADDR_AINV_BASE, Ainv, NN);
    int st = (int)RD(dev, XMAT_INV_KERNEL_CTRL_ADDR_STATUS_DATA);
    if (cycles) *cycles = (uint32_t)(ik_timer_read() - t0);
    return st;
}
