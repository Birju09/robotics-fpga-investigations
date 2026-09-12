#ifndef IK_KERNELS_HPP
#define IK_KERNELS_HPP

#include "ik_config.hpp"
#include "ik_types.hpp"

//
//! The two inverse-kinematics solvers under comparison.
//! analytic() - closed form, data-independent latency (worst case == typical case).
//! dls()      - iterates until residual < tol or max_iter; its iteration
//! count (and distribution over a workload) is what this project measures.
//

namespace iks {

    //! Analytic closed-form IK. cfg selects among 8 branches via IK_CFG_* bits.
    //! Returns IK_OK or IK_ERR_UNREACH.
    int analytic(const ik_real_t Rd[3][3], const ik_real_t pd[3], int cfg,
                 ik_real_t q[IK_DOF]);

    //! Damped least-squares IK, iterating from q_seed.
    //! @param step_max Trust region: bound on |dq|_inf per iteration, radians.
    //! 0 disables the clamp. Scaled by a power of two, so the
    //! effective bound is (step_max/2, step_max] - see ik_config.hpp.
    //! @param resid Final task-space residual (2-norm, Q16.16).
    //! Returns IK_OK, IK_ERR_NO_CONV, or IK_ERR_SINGULAR.
    int dls(const ik_real_t Rd[3][3], const ik_real_t pd[3],
            const ik_real_t q_seed[IK_DOF], ik_real_t lambda, ik_real_t tol,
            int max_iter, ik_real_t step_max, ik_real_t q[IK_DOF], int* iters,
            ik_real_t* resid);

}  // namespace iks

//! @defgroup IKKernels Inverse Kinematics IP Kernels
//! Standalone IP cores for inverse kinematics computation.
//! @{

//! AXI4-Lite IP wrapper for analytic(). Q16.16 fixed-point.
//! pose is [x, y, z, roll, pitch, yaw], where R = Rz(yaw) Ry(pitch) Rx(roll).
extern "C" void ik_analytic_kernel(const ik_word_t pose[IK_DOF], int cfg,
                                   ik_word_t q[IK_DOF], int* status);

//! AXI4-Lite IP wrapper for dls(). Q16.16 fixed-point.
//! pose is [x, y, z, roll, pitch, yaw], where R = Rz(yaw) Ry(pitch) Rx(roll).
extern "C" void ik_dls_kernel(const ik_word_t pose[IK_DOF],
                              const ik_word_t q_seed[IK_DOF], ik_word_t lambda,
                              ik_word_t tol, int max_iter, ik_word_t step_max,
                              ik_word_t q[IK_DOF], int* iters,
                              ik_word_t* resid, int* status);
//! @}

#endif  //! IK_KERNELS_HPP
