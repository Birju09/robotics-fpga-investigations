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
 */

#define RD(d, off)      Xil_In32((d)->base + (off))
#define WR(d, off, val) Xil_Out32((d)->base + (off), (uint32_t)(val))

/* ---------------- cycle counter ---------------- */
/*
 * The Cortex-A9 global timer runs at CPU_3x2x, i.e. half the CPU clock
 * (333 MHz on a -1 speed grade Zynq-7020 at the usual 667 MHz).  XTime_GetTime
 * returns its 64-bit count.  Resolution is ~3 ns, which is fine against
 * kernel latencies in the microsecond range.
 */
void ik_timer_init(void)
{
#if defined(IK_HAVE_XILTIMER)
    /* Unlike the old standalone BSP's xtime_l.h, where XTime_GetTime lazily
     * arms the global timer counter on its first call, libxiltimer's
     * XTime_GetTime is a thin wrapper over a driver instance that must be
     * brought up first.  Skipping this leaves the counter it reads never
     * started, so every XTime_GetTime call returns the same value and every
     * measured delta in this harness comes out as exactly zero - which is
     * what happens without this call. */
    XTimer_Init();
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
 * Launch and wait.  Returns elapsed global-timer ticks.
 *
 * The kernels finish in microseconds, so a bounded spin is both the lowest
 * overhead and the lowest jitter option.  The bound exists so a mis-mapped
 * register or an unprogrammed PL cannot hang the application.
 */
static uint64_t ik_run(ik_dev_t *dev, int *timed_out)
{
    const uint64_t limit = (uint64_t)ik_timer_hz();   /* 1 second */
    uint64_t t0, t1;

    if (timed_out) *timed_out = 0;

    t0 = ik_timer_read();
    WR(dev, IK_ADDR_AP_CTRL, IK_AP_START);

    for (;;) {
        uint32_t c = RD(dev, IK_ADDR_AP_CTRL);
        if (c & IK_AP_DONE)
            break;
        t1 = ik_timer_read();
        if ((t1 - t0) > limit) {
            if (timed_out) *timed_out = 1;
            break;
        }
    }
    t1 = ik_timer_read();

    /* ap_done is clear-on-read for the AP_CTRL_HS protocol HLS generates
     * here; the read above already acknowledged it. */
    return t1 - t0;
}

/* ---------------- analytic IK ---------------- */

int ik_analytic_solve(ik_dev_t *dev, const float pose[6], int cfg,
                      float q_out[6], uint32_t *cycles)
{
    int to = 0;

    ik_write_block(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_POSE_BASE, pose, IK_DOF);
    WR(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_CFG_DATA, cfg);

    uint64_t dt = ik_run(dev, &to);
    if (cycles) *cycles = (uint32_t)dt;
    if (to) return -1;

    ik_read_block(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_Q_BASE, q_out, IK_DOF);
    return (int)RD(dev, XIK_ANALYTIC_KERNEL_CTRL_ADDR_STATUS_DATA);
}

/* ---------------- DLS IK ---------------- */

int ik_dls_solve(ik_dev_t *dev, const float pose[6], const float q_seed[6],
                 float lambda, float tol, int max_iter,
                 float q_out[6], int *iters, float *resid, uint32_t *cycles)
{
    int to = 0;

    ik_write_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_POSE_BASE, pose, IK_DOF);
    ik_write_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_Q_SEED_BASE, q_seed, IK_DOF);
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_LAMBDA_DATA,   ik_f2q(lambda));
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_TOL_DATA,      ik_f2q(tol));
    WR(dev, XIK_DLS_KERNEL_CTRL_ADDR_MAX_ITER_DATA, max_iter);

    uint64_t dt = ik_run(dev, &to);
    if (cycles) *cycles = (uint32_t)dt;
    if (to) return -1;

    ik_read_block(dev, XIK_DLS_KERNEL_CTRL_ADDR_Q_BASE, q_out, IK_DOF);
    if (iters)
        *iters = (int)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_ITERS_DATA);
    if (resid)
        *resid = ik_q2f((ik_word_t)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_RESID_DATA));
    return (int)RD(dev, XIK_DLS_KERNEL_CTRL_ADDR_STATUS_DATA);
}

/* ---------------- matrix multiply ---------------- */

int ik_matmul(ik_dev_t *dev, int m, int k, int n, int ta, int tb,
              const float *A, const float *B, float *C, uint32_t *cycles)
{
    int to = 0;
    const int NN = IK_MAT_MAX * IK_MAT_MAX;

    ik_write_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_A_BASE, A, NN);
    ik_write_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_B_BASE, B, NN);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_M_DATA,  m);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_K_DATA,  k);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_N_DATA,  n);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_TA_DATA, ta);
    WR(dev, XMAT_MUL_KERNEL_CTRL_ADDR_TB_DATA, tb);

    uint64_t dt = ik_run(dev, &to);
    if (cycles) *cycles = (uint32_t)dt;
    if (to) return -1;

    ik_read_block(dev, XMAT_MUL_KERNEL_CTRL_ADDR_C_BASE, C, NN);
    return (int)RD(dev, XMAT_MUL_KERNEL_CTRL_ADDR_STATUS_DATA);
}

/* ---------------- matrix inversion ---------------- */

int ik_matinv(ik_dev_t *dev, int n, const float *A, float *Ainv,
              uint32_t *cycles)
{
    int to = 0;
    const int NN = IK_MAT_MAX * IK_MAT_MAX;

    ik_write_block(dev, XMAT_INV_KERNEL_CTRL_ADDR_A_BASE, A, NN);
    WR(dev, XMAT_INV_KERNEL_CTRL_ADDR_N_DATA, n);

    uint64_t dt = ik_run(dev, &to);
    if (cycles) *cycles = (uint32_t)dt;
    if (to) return -1;

    ik_read_block(dev, XMAT_INV_KERNEL_CTRL_ADDR_AINV_BASE, Ainv, NN);
    return (int)RD(dev, XMAT_INV_KERNEL_CTRL_ADDR_STATUS_DATA);
}
