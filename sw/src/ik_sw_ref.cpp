/*
 * Software reference running on the Cortex-A9.
 *
 * This deliberately compiles the *same* translation units the PL kernels are
 * built from (the .cpp files under hls/src) with -DIK_USE_FLOAT.  Reimplementing
 * separately would make the PS-vs-PL comparison a comparison of two different
 * programs; sharing the source means any latency difference is attributable to
 * the target, not to the algorithm.
 *
 * Note the numeric type differs by construction: the PL runs Q16.16 fixed
 * point, the PS runs double.  That is the honest comparison - it is what each
 * target is good at - but it means the two will not agree bit for bit, and the
 * PS figure carries the cost of double-precision software floating point on a
 * VFPv3 core.  See docs/timing_methodology.md.
 */

#include "ik_kernels.hpp"
#include "kinematics.hpp"

extern "C" {
#include "ik_driver.h"
}

extern "C" int ik_analytic_solve_sw(const float pose[6], int cfg, float q_out[6])
{
    ik_real_t pd[3], rpy[3], Rd[3][3], q[IK_DOF];

    for (int i = 0; i < 3; i++) {
        pd[i]  = (ik_real_t)pose[i];
        rpy[i] = (ik_real_t)pose[3 + i];
    }
    ikk::rpy_to_rot(rpy[0], rpy[1], rpy[2], Rd);

    int st = iks::analytic(Rd, pd, cfg, q);

    for (int i = 0; i < IK_DOF; i++)
        q_out[i] = (float)q[i];
    return st;
}

extern "C" int ik_dls_solve_sw(const float pose[6], const float q_seed[6],
                               float lambda, float tol, int max_iter,
                               float q_out[6], int *iters, float *resid)
{
    ik_real_t pd[3], rpy[3], Rd[3][3], qs[IK_DOF], q[IK_DOF], rr;

    for (int i = 0; i < 3; i++) {
        pd[i]  = (ik_real_t)pose[i];
        rpy[i] = (ik_real_t)pose[3 + i];
    }
    for (int i = 0; i < IK_DOF; i++)
        qs[i] = (ik_real_t)q_seed[i];

    ikk::rpy_to_rot(rpy[0], rpy[1], rpy[2], Rd);

    int it = 0;
    int st = iks::dls(Rd, pd, qs, (ik_real_t)lambda, (ik_real_t)tol,
                      max_iter, q, &it, &rr);

    for (int i = 0; i < IK_DOF; i++)
        q_out[i] = (float)q[i];
    if (iters) *iters = it;
    if (resid) *resid = (float)rr;
    return st;
}
