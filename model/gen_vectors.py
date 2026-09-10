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

    # Two thirds ordinary poses.
    while len(recs) < (2 * n) // 3:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        T = M.fk(q)
        if M.elbow_conditioning(T, +1) < 0.25:
            continue
        recs.append(_mk_rec(q, T, "well-conditioned"))

    # One third stressed: wrist near singularity, or elbow near full extension.
    while len(recs) < n:
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
 * The set intentionally mixes well-conditioned poses with wrist-singular and
 * near-full-extension ones, because the DLS iteration count - and therefore
 * its latency - is governed by the hard cases, not the typical ones.
 */
#ifndef IK_VECTORS_H
#define IK_VECTORS_H

#include <stdint.h>

#define IK_NVEC %d

/* { x, y, z, roll, pitch, yaw } */
static const int32_t ik_pose_tbl[IK_NVEC][6] = {
%s
};

/* DLS seed: the true solution perturbed, as a tracking loop would supply. */
static const int32_t ik_seed_tbl[IK_NVEC][6] = {
%s
};

/* Reference joint solution from the floating-point model. */
static const int32_t ik_qgold_tbl[IK_NVEC][6] = {
%s
};

/* Analytic branch selector (IK_CFG_* bits). */
static const uint8_t ik_cfg_tbl[IK_NVEC] = {
%s
};

static const char *const ik_tag_tbl[IK_NVEC] = {
%s
};

#endif /* IK_VECTORS_H */
""" % (len(recs),
            c_rows("pose"),
            c_rows("seed"),
            c_rows("qgold"),
            "    " + ", ".join("%d" % r["cfg"] for r in recs),
            "\n".join('    "%s",' % r["tag"] for r in recs)))

    print(f"  {'ik_vectors.h':20s} {len(recs):5d} poses   -> {os.path.relpath(path)}")


def _mk_rec(q, T, tag):
    pose = q_arr(M.T_to_pose(T))
    Tq = np.eye(4)
    Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
    Tq[:3, 3] = pose[:3]
    qs, ok = M.ik_analytic(Tq, +1, +1, +1)
    if not ok:
        qs = q
    seed = q_arr(q + rng.uniform(-0.25, 0.25, size=6))
    return {"pose": pose, "seed": seed, "qgold": qs,
            "cfg": 1 | 2 | 4, "tag": tag}


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
