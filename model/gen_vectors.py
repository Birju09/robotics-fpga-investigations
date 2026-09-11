#!/usr/bin/env python3
"""
Emit golden test vectors consumed by the HLS C++ testbenches.

Files written to hls/tb/vectors/ (plain text, one record per line, Q16.16 as
signed decimal so the C++ side can read them with a plain istream):

  matmul.txt        M K N ta tb  A[...]  B[...]  Cgold[...]
  matinv.txt        N  A[...]  Ainv_gold[...]
  fk.txt            q[6]  T[12]            (row-major 3x4)
  ik_analytic.txt   pose[6] cfg  q_gold[6]  |sin(gamma)|
  ik_dls.txt        pose[6] qseed[6]  q_gold[6] iters_gold
"""

import os
import numpy as np
import ik_model as M

OUT = os.path.join(os.path.dirname(__file__), "..", "hls", "tb", "vectors")
rng = np.random.default_rng(0x1153)

QLIM = np.array([
    [-2.9, 2.9], [-1.9, 1.9], [-2.6, 2.6],
    [-2.9, 2.9], [-1.9, 1.9], [-2.9, 2.9],
])

LAMBDA = 0.02
TOL = 1e-3
MAX_ITER = 64


def fx(x):
    """float -> signed decimal Q16.16 (as the C++ tb will parse it)."""
    v = M.to_fixed(x)
    return v - (1 << 32) if v >= (1 << 31) else v


def row(*vals):
    out = []
    for v in vals:
        a = np.asarray(v).reshape(-1)
        out.extend(str(fx(x)) for x in a)
    return " ".join(out)


def q_arr(x):
    """Quantise an array through Q16.16 so the golden value is reproducible."""
    return np.array([M.from_fixed(M.to_fixed(v)) for v in np.asarray(x).reshape(-1)])


def gen_matmul(n=64, path="matmul.txt"):
    lines = []
    for _ in range(n):
        Md, K, N = rng.integers(1, 7, size=3)
        ta, tb = rng.integers(0, 2, size=2)
        # A is (K x M) when transposed, else (M x K); likewise for B.
        A = q_arr(rng.uniform(-2, 2, size=(K, Md) if ta else (Md, K))).reshape(
            (K, Md) if ta else (Md, K))
        B = q_arr(rng.uniform(-2, 2, size=(N, K) if tb else (K, N))).reshape(
            (N, K) if tb else (K, N))
        C = (A.T if ta else A) @ (B.T if tb else B)
        lines.append(f"{Md} {K} {N} {ta} {tb} " + row(A, B, C))
    write(path, lines)


def gen_matinv(n=64, path="matinv.txt"):
    lines = []
    made = 0
    while made < n:
        N = int(rng.integers(1, 7))
        # Build well-conditioned SPD matrices of the same shape DLS produces.
        Jr = rng.uniform(-1, 1, size=(N, N))
        A = q_arr(Jr @ Jr.T + (LAMBDA ** 2 + 0.05) * np.eye(N)).reshape(N, N)
        if np.linalg.cond(A) > 5e3:
            continue
        Ai = np.linalg.inv(A)
        if np.abs(Ai).max() > 20000:      # must fit Q16.16
            continue
        lines.append(f"{N} " + row(A, Ai))
        made += 1
    write(path, lines)


def gen_fk(n=128, path="fk.txt"):
    lines = []
    for q in rng.uniform(QLIM[:, 0], QLIM[:, 1], size=(n, 6)):
        q = q_arr(q)
        T = M.fk(q)
        lines.append(row(q, T[:3, :4]))
    write(path, lines)


def gen_ik_analytic(n=256, path="ik_analytic.txt"):
    lines = []
    made = 0
    while made < n:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        T = M.fk(q)
        pose = q_arr(M.T_to_pose(T))
        # Re-derive T from the quantised pose so tb input == kernel input.
        Tq = np.eye(4)
        Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
        Tq[:3, 3] = pose[:3]
        sh, el, wr = rng.choice([-1, 1], size=3)
        qs, ok = M.ik_analytic(Tq, sh, el, wr)
        if not ok:
            continue
        cfg = (1 if sh > 0 else 0) | (2 if el > 0 else 0) | (4 if wr > 0 else 0)
        # |sin(gamma)| lets the testbench tell an ill-conditioned pose from a
        # wrong answer; see ik_model.elbow_conditioning().
        cond = M.elbow_conditioning(Tq, sh)
        lines.append(row(pose) + f" {cfg} " + row(qs) + " " + row(cond))
        made += 1
    write(path, lines)


def gen_ik_dls(n=256, path="ik_dls.txt"):
    lines = []
    made = 0
    while made < n:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        T = M.fk(q)
        pose = q_arr(M.T_to_pose(T))
        Tq = np.eye(4)
        Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
        Tq[:3, 3] = pose[:3]
        seed = q_arr(q + rng.uniform(-0.3, 0.3, size=6))
        qs, iters, err, ok = M.ik_dls(Tq, seed, lam=LAMBDA,
                                      max_iter=MAX_ITER, tol=TOL)
        if not ok:
            continue
        lines.append(row(pose, seed, qs) + f" {iters}")
        made += 1
    write(path, lines)


def gen_c_header(n=48, path=None):
    """
    Emit sw/src/ik_vectors.h: the pose set the bare-metal timing harness runs.

    Deliberately mixes well-conditioned poses with near-singular ones.  The
    latter are the whole point - a workload of only comfortable poses would
    show the DLS solver converging in 3-4 iterations every time and hide the
    tail that determines whether it can be scheduled at all.
    """
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "..", "sw", "src",
                            "ik_vectors.h")

    recs = []

    # Half ordinary poses.
    while len(recs) < n // 2:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        T = M.fk(q)
        if M.elbow_conditioning(T, +1) < 0.25:
            continue
        recs.append(_mk_rec(q, T, "well-conditioned"))

    # A quarter stressed GEOMETRICALLY: wrist near singularity, or elbow near
    # full extension.  These stress the analytic solver, whose branch
    # selection degenerates there.
    while len(recs) < (3 * n) // 4:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        if len(recs) % 2 == 0:
            q[4] = rng.choice([1e-3, -1e-3, 2e-3])        # wrist singularity
            tag = "wrist-singular"
        else:
            # Push the wrist centre close to full extension.
            q[2] = M.wrap_pi(M.PHI + rng.uniform(-0.02, 0.02))
            tag = "elbow-extended"
        q = q_arr(q)
        T = M.fk(q)
        pose = q_arr(M.T_to_pose(T))
        Tq = np.eye(4)
        Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
        Tq[:3, 3] = pose[:3]
        if not M.ik_analytic(Tq)[1]:
            continue
        recs.append(_mk_rec(q, T, tag))

    # A quarter selected on DLS ITERATION COUNT, which is a different thing
    # entirely and is the reason this block exists.
    #
    # The geometric stress above is a proxy for the analytic solver's
    # difficulty, not the iterative one's.  DLS does not care that a branch
    # has degenerated; it cares how far the seed is from the solution and how
    # flat the Jacobian stays along the way, and those are not the same poses.
    # Combined with _mk_rec()'s +-0.25 rad seed - a warm start, near the
    # answer by construction - every pose in the old table converged in three
    # to five iterations.  Measured on hardware the table topped out at 7,
    # giving max/med = 1.77, while the 256-pose host sweep in gen_ik_dls()
    # reaches 25.  The table was not sampling the tail that decides whether
    # the solver can be scheduled at all.
    #
    # So select on the quantity being measured: run the model's DLS and keep
    # poses by their actual iteration count, spread across the tail rather
    # than clustered at the worst case.
    for _iters, q, T, seed in _find_dls_tail(n - len(recs)):
        recs.append(_mk_rec(q, T, "dls-tail", seed=seed))

    def c_rows(key):
        return "\n".join(
            "    { %s }," % ", ".join("%11d" % fx(v) for v in r[key])
            for r in recs)

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("""/*
 * GENERATED FILE - do not edit.  Produced by model/gen_vectors.py.
 *
 * Pose set for the bare-metal timing harness.  All values are signed Q16.16,
 * the format the kernels' AXI4-Lite registers expect.
 *
 * The set mixes three groups, in halves and quarters:
 *
 *   well-conditioned              ordinary poses, comfortable for both solvers
 *   wrist-singular/elbow-extended geometric stress, which is what degrades the
 *                                 ANALYTIC solver's branch selection
 *   dls-tail                      selected on the model's DLS iteration count,
 *                                 which is a different property entirely
 *
 * The last group exists because the first two do not produce one.  DLS does
 * not care that a branch has degenerated; it cares how far the seed is and
 * how flat the Jacobian stays, so a table stressed only geometrically - and
 * seeded within 0.25 rad of the answer - converges in three to five
 * iterations every time.  Measured on hardware that gave max/med = 1.77,
 * against 25 iterations in the 256-pose host sweep.  The tail that decides
 * whether an iterative solver can be scheduled at all was not being sampled.
 */
#ifndef IK_VECTORS_H
#define IK_VECTORS_H

#include <stdint.h>

/* Bumped whenever this file grows something a consumer might require.  It is
 * gitignored - reproduced by model/gen_vectors.py rather than checked in - so
 * pulling a commit updates the code that reads it but not the file itself,
 * and the resulting error is otherwise just an undeclared identifier.
 *
 *   1  ik_model_iters_tbl, DLS tail poses
 *   2  ik_qdls_tbl
 */
#define IK_VECTORS_VERSION 2

#define IK_NVEC %d

/* { x, y, z, roll, pitch, yaw } */
static const int32_t ik_pose_tbl[IK_NVEC][6] = {
%s
};

/* DLS seed: the true solution perturbed, as a tracking loop would supply. */
static const int32_t ik_seed_tbl[IK_NVEC][6] = {
%s
};

/* Reference joint solution from the floating-point model's ANALYTIC solver,
 * branch (+,+,+) - i.e. what ik_cfg_tbl selects.  Compare ik_analytic_kernel
 * against this. */
static const int32_t ik_qgold_tbl[IK_NVEC][6] = {
%s
};

/* Reference joint solution from the floating-point model's DLS solver, run
 * from ik_seed_tbl.  Compare ik_dls_kernel against THIS, not against
 * ik_qgold_tbl: IK is multi-valued, and the two solvers legitimately land on
 * different branches - often ~pi apart in a wrist joint - while both solve
 * the pose to within tolerance.  A DLS result checked against the analytic
 * branch measures which solution was found, not how accurately. */
static const int32_t ik_qdls_tbl[IK_NVEC][6] = {
%s
};

/* Analytic branch selector (IK_CFG_* bits). */
static const uint8_t ik_cfg_tbl[IK_NVEC] = {
%s
};

static const char *const ik_tag_tbl[IK_NVEC] = {
%s
};

/* Iterations the double-precision model needs for each pose at the same
 * lambda and tol the harness passes.  The Q16.16 kernel should track this
 * closely; a systematic disagreement means quantisation moved the
 * convergence path, which is a finding rather than a test failure. */
static const uint8_t ik_model_iters_tbl[IK_NVEC] = {
%s
};

#endif /* IK_VECTORS_H */
""" % (len(recs),
            c_rows("pose"),
            c_rows("seed"),
            c_rows("qgold"),
            c_rows("qdls"),
            "    " + ", ".join("%d" % r["cfg"] for r in recs),
            "\n".join('    "%s",' % r["tag"] for r in recs),
            "    " + ", ".join("%d" % r["iters"] for r in recs)))

    print(f"  {'ik_vectors.h':20s} {len(recs):5d} poses   -> {os.path.relpath(path)}")


def _pose_to_T(pose):
    Tq = np.eye(4)
    Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
    Tq[:3, 3] = pose[:3]
    return Tq


def _find_dls_tail(count, floor=8, ceiling=MAX_ITER // 2,
                   pool_target=30, tries_max=60000):
    """
    Poses whose DLS iteration count actually lands in the tail.

    Seeds are drawn far (+-1.2 rad) rather than at _mk_rec()'s +-0.25: the
    seed distance is the dominant term in how many iterations DLS needs, so a
    tight perturbation cannot produce a tail no matter which pose it is
    applied to.

    Returns a spread across the tail, not `count` copies of the worst case.
    A table of only extreme poses would overstate the median as badly as the
    old one understated the maximum; what the real-time argument needs is the
    whole distribution, max included.
    """
    pool = []
    tries = 0
    while tries < tries_max and len(pool) < pool_target * count:
        tries += 1
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        T = M.fk(q)
        pose = q_arr(M.T_to_pose(T))
        Tq = _pose_to_T(pose)
        if not M.ik_analytic(Tq)[1]:
            continue
        seed = q_arr(q + rng.uniform(-1.2, 1.2, size=6))
        _, iters, _, conv = M.ik_dls(Tq, seed, lam=LAMBDA,
                                     max_iter=MAX_ITER, tol=TOL)
        # Converged, and with margin below the cap.  A pose that reaches
        # MAX_ITER returns IK_ERR_NO_CONV and would make the correctness
        # column meaningless; one that lands exactly ON it is worse, because
        # the model checks its residual once more after the loop and the
        # kernel does not, so the two disagree about whether it converged at
        # all.  Non-convergence is a real phenomenon and worth its own
        # experiment - it is just not this one.
        if conv and floor <= iters <= ceiling:
            pool.append((iters, q, T, seed))

    if len(pool) <= count:
        return pool

    pool.sort(key=lambda r: r[0])
    # Evenly spaced through the sorted pool, so the slowest pose found is
    # always the last one taken.
    idx = sorted({round(i * (len(pool) - 1) / (count - 1))
                  for i in range(count)})
    return [pool[i] for i in idx]


def _mk_rec(q, T, tag, seed=None):
    pose = q_arr(M.T_to_pose(T))
    Tq = _pose_to_T(pose)
    qs, ok = M.ik_analytic(Tq, +1, +1, +1)
    if not ok:
        qs = q
    if seed is None:
        seed = q_arr(q + rng.uniform(-0.25, 0.25, size=6))
    # The model's own DLS result from this seed, and what it cost.
    #
    # qdls is NOT qgold, and the difference is the point.  IK is multi-valued:
    # qgold is the analytic (+,+,+) branch, while DLS converges to whichever
    # branch its seed is nearest, so the two routinely differ by ~pi in a
    # wrist joint while both solve the pose exactly.  Comparing a DLS result
    # against qgold measures branch choice, not accuracy.
    qd, iters, _, _ = M.ik_dls(Tq, seed, lam=LAMBDA, max_iter=MAX_ITER, tol=TOL)
    return {"pose": pose, "seed": seed, "qgold": qs, "qdls": q_arr(qd),
            "iters": int(iters), "cfg": 1 | 2 | 4, "tag": tag}


def write(name, lines):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    with open(p, "w") as f:
        f.write(f"{len(lines)}\n")
        f.write("\n".join(lines) + "\n")
    print(f"  {name:20s} {len(lines):5d} vectors -> {os.path.relpath(p)}")


if __name__ == "__main__":
    print("Generating golden vectors (Q16.16):")
    gen_matmul()
    gen_matinv()
    gen_fk()
    gen_ik_analytic()
    gen_ik_dls()
    gen_c_header()
    print("done.")
