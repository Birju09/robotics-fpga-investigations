#ifndef MATMUL_HPP
#define MATMUL_HPP

#include "ik_types.hpp"
#include "ik_config.hpp"

/*
 * Matrix multiply, shared between the standalone mat_mul IP and the DLS solver.
 *
 * The DLS core calls mm::multiply() directly rather than issuing AXI
 * transactions to the standalone IP.  Routing every product of a 6x6 through
 * the PS or an AXI master would add hundreds of cycles per iteration to a loop
 * that runs 4-64 times per solve, and would make the latency depend on
 * interconnect arbitration - which is precisely the thing this investigation is
 * trying to measure cleanly.  Same source, two export targets.
 */

namespace mm {

/*
 * C = op(A) * op(B),  op(X) = X^T when the corresponding flag is set.
 *
 *   A is (ta ? k x m : m x k)
 *   B is (tb ? n x k : k x n)
 *   C is m x n
 *
 * Storage is always row-major in a fixed IK_MAT_MAX x IK_MAT_MAX array; the
 * runtime dimensions select which part is meaningful.
 *
 * The loops run to IK_MAT_MAX unconditionally and predicate on the runtime
 * dimensions instead of bounding the trip count.  That costs the full 6x6x6
 * work on every call regardless of the requested size, and buys a latency that
 * is identical for every input - the property the analytic path is being
 * judged on.
 */
void multiply(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX],
              const ik_real_t B[IK_MAT_MAX][IK_MAT_MAX],
              int m, int k, int n, bool ta, bool tb,
              ik_real_t C[IK_MAT_MAX][IK_MAT_MAX]);

} /* namespace mm */

/* ---------------- standalone IP top level ---------------- */
/*
 * AXI4-Lite mapped kernel.  A, B and C are Q16.16 words in row-major order,
 * exposed as addressable register banks in the control interface.
 */
extern "C" void mat_mul_kernel(int m, int k, int n, int ta, int tb,
                               const ik_word_t A[IK_MAT_MAX * IK_MAT_MAX],
                               const ik_word_t B[IK_MAT_MAX * IK_MAT_MAX],
                               ik_word_t C[IK_MAT_MAX * IK_MAT_MAX],
                               int *status);

#endif /* MATMUL_HPP */
