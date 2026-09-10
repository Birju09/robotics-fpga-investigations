/*
 * Bare-metal timing harness for the IK kernels.
 *
 * What this measures, and what it does not
 * ---------------------------------------
 * The number reported per solve is wall-clock time from the PS, measured with
 * the Cortex-A9 global timer around the whole transaction: argument writes,
 * ap_start, the poll loop, and the result reads.  That is the latency a
 * control loop actually experiences, and for these kernels the AXI4-Lite
 * register traffic is a substantial part of it - roughly 20 single-beat
 * transactions per solve, each costing tens of PS cycles.
 *
 * It is therefore NOT the same as the kernel latency Vitis HLS reports.  The
 * HLS figure is PL cycles between ap_start and ap_done; this one includes the
 * bus overhead the HLS figure excludes.  Both matter, and the gap between them
 * is itself a finding: for a kernel this small, moving data can cost more than
 * computing.  Compare against `make -C hls reports`.
 *
 * Scope: ik_analytic only for now.  ik_dls_kernel still does not fit on the
 * xc7z020 standalone (~126% DSP utilisation after the resource-sharing work
 * in hls/src/ik_dls.cpp and hls/include/ik_math.hpp) - see the README
 * Status section.  scripts/build_vivado.tcl defaults to a bitstream with
 * mat_mul_kernel, mat_inv_kernel and ik_analytic_kernel only, so DLS's AXI
 * base address does not exist in xparameters.h for this build; the DLS
 * driver code (ik_driver.c, ik_sw_ref.cpp) is left in place for when it
 * fits and this harness is extended back to cover it.
 */

#include <stdio.h>
#include <string.h>

#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xil_io.h"

#include "ik_driver.h"
#include "ik_regmap.h"
#include "ik_vectors.h"

static void init_platform_stub(void);

/* ------------------------------------------------------------------ */
/* Base addresses.                                                     */
/*                                                                     */
/* xparameters.h names these from the block design instance names, but */
/* the exact macro spelling varies with the Vivado/Vitis version, so   */
/* resolve them with a small cascade rather than assuming one form.    */
/* ------------------------------------------------------------------ */
#if defined(XPAR_IK_ANALYTIC_KERNEL_0_S_AXI_CTRL_BASEADDR)
#define ANALYTIC_BASE XPAR_IK_ANALYTIC_KERNEL_0_S_AXI_CTRL_BASEADDR
#elif defined(XPAR_XIK_ANALYTIC_KERNEL_0_S_AXI_CTRL_BASEADDR)
#define ANALYTIC_BASE XPAR_XIK_ANALYTIC_KERNEL_0_S_AXI_CTRL_BASEADDR
#elif defined(XPAR_IK_ANALYTIC_KERNEL_0_BASEADDR)
#define ANALYTIC_BASE XPAR_IK_ANALYTIC_KERNEL_0_BASEADDR
#else
#error "Cannot find the analytic kernel base address. Check xparameters.h and update these guards."
#endif

/* ik_dls_kernel is not in the default bitstream (see the file header comment
 * above) - no base address to resolve here until it is added back. */

#if defined(XPAR_MAT_MUL_KERNEL_0_S_AXI_CTRL_BASEADDR)
#define MATMUL_BASE XPAR_MAT_MUL_KERNEL_0_S_AXI_CTRL_BASEADDR
#elif defined(XPAR_MAT_MUL_KERNEL_0_BASEADDR)
#define MATMUL_BASE XPAR_MAT_MUL_KERNEL_0_BASEADDR
#else
#define MATMUL_BASE 0
#endif

#if defined(XPAR_MAT_INV_KERNEL_0_S_AXI_CTRL_BASEADDR)
#define MATINV_BASE XPAR_MAT_INV_KERNEL_0_S_AXI_CTRL_BASEADDR
#elif defined(XPAR_MAT_INV_KERNEL_0_BASEADDR)
#define MATINV_BASE XPAR_MAT_INV_KERNEL_0_BASEADDR
#else
#define MATINV_BASE 0
#endif

/* ------------------------------------------------------------------ */
/* Statistics                                                          */
/* ------------------------------------------------------------------ */
typedef struct {
    uint32_t v[IK_NVEC];
    int      n;
} samples_t;

static void samp_reset(samples_t *s) { s->n = 0; }

static void samp_add(samples_t *s, uint32_t x)
{
    if (s->n < IK_NVEC) s->v[s->n++] = x;
}

static void samp_sort(samples_t *s)
{
    for (int i = 1; i < s->n; i++) {
        uint32_t k = s->v[i];
        int j = i - 1;
        while (j >= 0 && s->v[j] > k) { s->v[j + 1] = s->v[j]; j--; }
        s->v[j + 1] = k;
    }
}

/* Ticks -> nanoseconds.  The global timer runs at COUNTS_PER_SECOND. */
static uint32_t ticks_to_ns(uint32_t ticks)
{
    return (uint32_t)(((uint64_t)ticks * 1000000000ULL) / ik_timer_hz());
}

static void samp_report(const char *label, samples_t *s)
{
    if (s->n == 0) { xil_printf("  %s: no samples\r\n", label); return; }
    samp_sort(s);

    uint64_t sum = 0;
    for (int i = 0; i < s->n; i++) sum += s->v[i];

    uint32_t mn  = s->v[0];
    uint32_t md  = s->v[s->n / 2];
    uint32_t p95 = s->v[(s->n * 95) / 100];
    uint32_t mx  = s->v[s->n - 1];
    uint32_t av  = (uint32_t)(sum / s->n);

    xil_printf("  %-26s n=%3d  min=%6u  med=%6u  p95=%6u  max=%6u  mean=%6u  (ns)\r\n",
               label, s->n,
               ticks_to_ns(mn), ticks_to_ns(md), ticks_to_ns(p95),
               ticks_to_ns(mx), ticks_to_ns(av));
    /* Jitter ratio is the real-time figure of merit: how much worse the worst
     * case is than the typical one. */
    if (md)
        xil_printf("  %-26s max/med = %u.%02u\r\n", "",
                   mx / md, ((mx * 100) / md) % 100);
}

/* ------------------------------------------------------------------ */
/* Register-map self check                                             */
/*                                                                     */
/* The offsets in ik_regmap.h are generated by scripts/gen_regmap.py    */
/* from the HLS output.  If that step was skipped, the checked-in       */
/* provisional map may not match this build - so prove it before        */
/* trusting a single measurement.                                       */
/* ------------------------------------------------------------------ */
static int verify_regmap(ik_dev_t *dev, uintptr_t in_base, uintptr_t out_base,
                         const char *name)
{
    const int32_t pat[6] = { 0x00010000, (int32_t)0xFFFF0000, 0x12345678,
                             0x0000BEEF, (int32_t)0x80000001, 0x7FFFFFFF };
    int ok = 1;

    for (int i = 0; i < 6; i++)
        Xil_Out32(dev->base + in_base + 4 * i, (uint32_t)pat[i]);

    for (int i = 0; i < 6; i++) {
        uint32_t got = Xil_In32(dev->base + in_base + 4 * i);
        if (got != (uint32_t)pat[i]) {
            xil_printf("  %s: input bank readback mismatch at +0x%02x: "
                       "wrote 0x%08x read 0x%08x\r\n",
                       name, (unsigned)(in_base + 4 * i),
                       (unsigned)pat[i], (unsigned)got);
            ok = 0;
        }
    }

    /* An idle kernel must report ap_idle. */
    uint32_t ctrl = Xil_In32(dev->base + IK_ADDR_AP_CTRL);
    if (!(ctrl & IK_AP_IDLE)) {
        xil_printf("  %s: ap_idle not set at reset (ap_ctrl=0x%08x)\r\n",
                   name, (unsigned)ctrl);
        ok = 0;
    }

    (void)out_base;
    return ok;
}

/* ------------------------------------------------------------------ */
static void print_pose_result(int i, const char *tag, int st,
                              const float *q, const float *qg)
{
    float worst = 0.0f;
    for (int j = 0; j < IK_DOF; j++) {
        float d = q[j] - qg[j];
        while (d >  3.14159265f) d -= 6.28318531f;
        while (d < -3.14159265f) d += 6.28318531f;
        if (d < 0) d = -d;
        if (d > worst) worst = d;
    }
    xil_printf("    [%2d] %-16s status=%d  worst joint delta = %d urad\r\n",
               i, tag, st, (int)(worst * 1e6f));
}

int main(void)
{
    ik_dev_t analytic = { ANALYTIC_BASE };
    ik_dev_t matmul   = { MATMUL_BASE };
    ik_dev_t matinv   = { MATINV_BASE };

    samples_t s_an_hw, s_an_sw, s_mm, s_mi;

    init_platform_stub();

    xil_printf("\r\n");
    xil_printf("==================================================================\r\n");
    xil_printf(" 6-DOF IK kernels - PL vs PS latency\r\n");
    xil_printf("==================================================================\r\n");
    xil_printf(" global timer      : %u Hz\r\n", (unsigned)ik_timer_hz());
    xil_printf(" poses             : %d\r\n", IK_NVEC);
#if IK_REGMAP_GENERATED
    xil_printf(" register map      : generated from this build\r\n");
#else
    xil_printf(" register map      : PROVISIONAL - run scripts/gen_regmap.py\r\n");
#endif
    xil_printf("\r\n");

    ik_timer_init();

    /* ---- register map sanity ---- */
    xil_printf("-- register map check --\r\n");
    int ok = 1;
    ok &= verify_regmap(&analytic,
                        XIK_ANALYTIC_KERNEL_CTRL_ADDR_POSE_BASE,
                        XIK_ANALYTIC_KERNEL_CTRL_ADDR_Q_BASE, "analytic");
    if (!ok) {
        xil_printf("\r\n  REGISTER MAP IS WRONG - refusing to report timings.\r\n");
        xil_printf("  Run: make -C hls ip && python3 scripts/gen_regmap.py\r\n");
        xil_printf("  then rebuild this application.\r\n");
        return 1;
    }
    xil_printf("  ok\r\n\r\n");

    /* ---- correctness pass ---- */
    xil_printf("-- correctness (first 8 poses) --\r\n");
    for (int i = 0; i < 8 && i < IK_NVEC; i++) {
        float pose[6], qg[6], q[6];
        for (int j = 0; j < 6; j++) {
            pose[j] = ik_q2f(ik_pose_tbl[i][j]);
            qg[j]   = ik_q2f(ik_qgold_tbl[i][j]);
        }
        int st = ik_analytic_solve(&analytic, pose, ik_cfg_tbl[i], q, NULL);
        print_pose_result(i, ik_tag_tbl[i], st, q, qg);
    }
    xil_printf("\r\n");

    /* ---- analytic: PL ---- */
    samp_reset(&s_an_hw);
    samp_reset(&s_an_sw);
    for (int i = 0; i < IK_NVEC; i++) {
        float pose[6], q[6];
        uint32_t c;
        for (int j = 0; j < 6; j++) pose[j] = ik_q2f(ik_pose_tbl[i][j]);

        ik_analytic_solve(&analytic, pose, ik_cfg_tbl[i], q, &c);
        samp_add(&s_an_hw, c);

        uint64_t t0 = ik_timer_read();
        ik_analytic_solve_sw(pose, ik_cfg_tbl[i], q);
        uint64_t t1 = ik_timer_read();
        samp_add(&s_an_sw, (uint32_t)(t1 - t0));
    }

    /* ---- standalone matrix IPs ---- */
    samp_reset(&s_mm);
    samp_reset(&s_mi);
    if (MATMUL_BASE && MATINV_BASE) {
        float A[36], B[36], C[36];
        for (int i = 0; i < 36; i++) {
            A[i] = 0.1f * (float)((i * 7) % 11) - 0.5f;
            B[i] = 0.1f * (float)((i * 5) % 13) - 0.6f;
        }
        for (int i = 0; i < 6; i++) A[i * 6 + i] += 2.0f;   /* keep it invertible */

        for (int r = 0; r < 16; r++) {
            uint32_t c;
            ik_matmul(&matmul, 6, 6, 6, 0, 0, A, B, C, &c);
            samp_add(&s_mm, c);
            ik_matinv(&matinv, 6, A, C, &c);
            samp_add(&s_mi, c);
        }
    }

    /* ---- report ---- */
    xil_printf("-- latency, PS wall clock around the whole transaction --\r\n");
    samp_report("analytic  (PL)", &s_an_hw);
    samp_report("analytic  (PS, double)", &s_an_sw);
    if (MATMUL_BASE && MATINV_BASE) {
        samp_report("mat_mul 6x6x6 (PL)", &s_mm);
        samp_report("mat_inv 6x6   (PL)", &s_mi);
    }

    xil_printf("\r\n");
    xil_printf("Read this against 'make -C hls reports': the HLS latency counts\r\n");
    xil_printf("PL cycles between ap_start and ap_done, the figures above add\r\n");
    xil_printf("the AXI4-Lite argument traffic.  For kernels this small the\r\n");
    xil_printf("difference is not a rounding error.\r\n");
    xil_printf("==================================================================\r\n");

    return 0;
}

/* Kept separate so the cache policy used for the measurements is explicit and
 * in one place: both caches on, which is the realistic configuration for a
 * control loop and the one that makes the PS reference honest. */
void init_platform_stub(void)
{
    Xil_DCacheEnable();
    Xil_ICacheEnable();
}
