//
//! PROVISIONAL register map - NOT yet generated from your build.
//
//! The offsets below follow the layout Vitis HLS conventionally produces for
//! these four signatures, but HLS assigns AXI4-Lite offsets from argument order
//! and width and the result is not guaranteed to match across releases or after
//! any edit to a kernel's arguments.
//
//! Before running on hardware:
//
//! make -C hls ip
//! python3 scripts/gen_regmap.py     # overwrites this file for real
//
//! ik_verify_regmap() in main.c performs a write/read-back check at startup and
//! refuses to continue if these are wrong, so a stale map fails loudly instead
//! of producing plausible-looking wrong answers.
//
#ifndef IK_REGMAP_H
#define IK_REGMAP_H

#define IK_REGMAP_GENERATED 0

//! Common to every HLS s_axilite control bank.
#define IK_ADDR_GIE 0x04
#define IK_ADDR_IER 0x08
#define IK_ADDR_ISR 0x0c

//! ---- mat_mul_kernel(m, k, n, ta, tb, A[36], B[36], C[36], status) ----
#define XMAT_MUL_KERNEL_CTRL_ADDR_AP_CTRL 0x00
#define XMAT_MUL_KERNEL_CTRL_ADDR_M_DATA 0x10
#define XMAT_MUL_KERNEL_CTRL_ADDR_K_DATA 0x18
#define XMAT_MUL_KERNEL_CTRL_ADDR_N_DATA 0x20
#define XMAT_MUL_KERNEL_CTRL_ADDR_TA_DATA 0x28
#define XMAT_MUL_KERNEL_CTRL_ADDR_TB_DATA 0x30
#define XMAT_MUL_KERNEL_CTRL_ADDR_STATUS_DATA 0x38
#define XMAT_MUL_KERNEL_CTRL_ADDR_A_BASE 0x100
#define XMAT_MUL_KERNEL_CTRL_ADDR_A_HIGH 0x18f
#define XMAT_MUL_KERNEL_CTRL_ADDR_B_BASE 0x200
#define XMAT_MUL_KERNEL_CTRL_ADDR_B_HIGH 0x28f
#define XMAT_MUL_KERNEL_CTRL_ADDR_C_BASE 0x300
#define XMAT_MUL_KERNEL_CTRL_ADDR_C_HIGH 0x38f

//! ---- mat_inv_kernel(n, A[36], Ainv[36], status) ----
#define XMAT_INV_KERNEL_CTRL_ADDR_AP_CTRL 0x00
#define XMAT_INV_KERNEL_CTRL_ADDR_N_DATA 0x10
#define XMAT_INV_KERNEL_CTRL_ADDR_STATUS_DATA 0x18
#define XMAT_INV_KERNEL_CTRL_ADDR_A_BASE 0x100
#define XMAT_INV_KERNEL_CTRL_ADDR_A_HIGH 0x18f
#define XMAT_INV_KERNEL_CTRL_ADDR_AINV_BASE 0x200
#define XMAT_INV_KERNEL_CTRL_ADDR_AINV_HIGH 0x28f

//! ---- ik_analytic_kernel(pose[6], cfg, q[6], status) ----
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_AP_CTRL 0x00
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_CFG_DATA 0x10
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_STATUS_DATA 0x18
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_POSE_BASE 0x40
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_POSE_HIGH 0x5f
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_Q_BASE 0x60
#define XIK_ANALYTIC_KERNEL_CTRL_ADDR_Q_HIGH 0x7f

//! ---- ik_dls_kernel(pose[6], q_seed[6], lambda, tol, max_iter, step_max,
//! q[6], iters, resid, status) ----
//
//! step_max (the DLS trust region) was added after this provisional map was
//! first written, which pushed every scalar after it up by one slot and moved
//! the three arrays.  That is exactly the kind of edit the generated map
//! exists for: run scripts/gen_regmap.py.  main.c's write/read-back check
//! catches it if you do not.
#define XIK_DLS_KERNEL_CTRL_ADDR_AP_CTRL 0x00
#define XIK_DLS_KERNEL_CTRL_ADDR_LAMBDA_DATA 0x10
#define XIK_DLS_KERNEL_CTRL_ADDR_TOL_DATA 0x18
#define XIK_DLS_KERNEL_CTRL_ADDR_MAX_ITER_DATA 0x20
#define XIK_DLS_KERNEL_CTRL_ADDR_STEP_MAX_DATA 0x28
#define XIK_DLS_KERNEL_CTRL_ADDR_ITERS_DATA 0x30
#define XIK_DLS_KERNEL_CTRL_ADDR_RESID_DATA 0x38
#define XIK_DLS_KERNEL_CTRL_ADDR_STATUS_DATA 0x40
#define XIK_DLS_KERNEL_CTRL_ADDR_POSE_BASE 0x80
#define XIK_DLS_KERNEL_CTRL_ADDR_POSE_HIGH 0x9f
#define XIK_DLS_KERNEL_CTRL_ADDR_Q_SEED_BASE 0xa0
#define XIK_DLS_KERNEL_CTRL_ADDR_Q_SEED_HIGH 0xbf
#define XIK_DLS_KERNEL_CTRL_ADDR_Q_BASE 0xc0
#define XIK_DLS_KERNEL_CTRL_ADDR_Q_HIGH 0xdf

#endif  //! IK_REGMAP_H
