#!/usr/bin/env python3
"""
Emit golden test vectors consumed by the HLS C++ testbenches.

Files written to hls/tb/vectors/ (plain text, one record per line, Q16.16 as
signed decimal so the C++ side can read with a plain istream):

  matmul.txt        M K N ta tb  A[...]  B[...]  Cgold[...]
  matinv.txt        N  A[...]  Ainv_gold[...]
  fk.txt            q[6]  T[12]            (row-major 3x4)
  ik_analytic.txt   pose[6] cfg  q_gold[6]  |sin(gamma)|  root/rho
  ik_dls.txt        pose[6] qseed[6]  q_gold[6] iters_gold
"""

import os
import numpy as np
import ik_model as M

OUT = os.path.join(os.path.dirname(__file__), "..", "hls", "tb", "vectors")
rng = np.random.default_rng(0x1153)

# Asymmetric joint limits do the job a collision model would otherwise have
# to: a2 and L3 differ by <1mm here, so only the q3 range stops the forearm
# folding flat onto the upper arm.
QLIM = M.ROBOT.qlim

LAMBDA = 0.02
TOL = 1e-3
MAX_ITER = 64

# ---------------------------------------------------------------------------
# Pentagon trajectory (see _pentagon_traj)
# ---------------------------------------------------------------------------
# Widest pentagon keeping elbow conditioning |sin gamma| > 0.45 everywhere, so
# the workload measures tracking rather than ill-conditioning. Re-run
# model/plot_dh.py after changing these.
TRAJ_CENTRE = (0.40, 0.0)     # metres, base x-y plane
TRAJ_Z = 0.15                 # tool-tip height above base plane
TRAJ_RADIUS = 0.20            # circumradius of the pentagon
TRAJ_STEPS = 100              # samples per edge -> 500 poses, ~2.35mm/step
# Yaw at each vertex, deg about vertical tool axis; interpolated linearly
# along each edge so orientation never turns discontinuously.
TRAJ_YAW_DEG = (0.0, 90.0, 0.0, 90.0, 0.0)


def fx(x):
    """float -> signed decimal Q16.16."""
    v = M.to_fixed(x)
    return v - (1 << 32) if v >= (1 << 31) else v


def row(*vals):
    out = []
    for v in vals:
        a = np.asarray(v).reshape(-1)
        out.extend(str(fx(x)) for x in a)
    return " ".join(out)


def q_arr(x):
    """Quantise an array through Q16.16."""
    return np.array([M.from_fixed(M.to_fixed(v)) for v in np.asarray(x).reshape(-1)])


def gen_matmul(n=64, path="matmul.txt"):
    lines = []
    for _ in range(n):
        Md, K, N = rng.integers(1, 7, size=3)
        ta, tb = rng.integers(0, 2, size=2)
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
        # well-conditioned SPD, same shape DLS produces
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
        # Elbow and shoulder conditioning, so the testbench can tell an
        # ill-conditioned pose from a wrong answer (see ik_model.py).
        cond = M.elbow_conditioning(Tq, sh)
        cond_sh = M.shoulder_conditioning(Tq)
        lines.append(row(pose) + f" {cfg} " + row(qs) + " " + row(cond)
                     + " " + row(cond_sh))
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
    Emit sw/src/ik_vectors.h, the timing harness's pose set.

    Deliberately mixes well-conditioned poses with near-singular ones: an
    all-comfortable workload would hide the DLS iteration tail that decides
    schedulability.
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
    # full extension. Degrades the analytic solver's branch selection.
    while len(recs) < (3 * n) // 4:
        q = q_arr(rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6))
        if len(recs) % 2 == 0:
            q[4] = rng.choice([1e-3, -1e-3, 2e-3])        # wrist singularity
            tag = "wrist-singular"
        else:
            # sin(gamma)=0 is theta3=-PHI here, not +PHI (alpha3 negative).
            # Verified against elbow_conditioning(): 0.094 at +PHI, 0 at -PHI.
            q[2] = M.wrap_pi(-M.PHI + rng.uniform(-0.02, 0.02))
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

    # A quarter selected on DLS ITERATION COUNT - a different axis than the
    # geometric stress above, which is a proxy for the analytic solver's
    # difficulty, not the iterative one's. DLS cares about seed distance and
    # Jacobian flatness, not branch degeneracy: a table stressed only
    # geometrically (and warm-seeded, per _mk_rec) topped out at 7 iterations
    # on hardware (max/med=1.77), while the 256-pose host sweep reaches 25.
    # So select directly on the model's actual iteration count, spread across
    # the tail rather than clustered at the worst case.
    for _iters, q, T, seed in _find_dls_tail(n - len(recs)):
        recs.append(_mk_rec(q, T, "dls-tail", seed=seed))

    # Deterministic - doesn't touch `rng`, so the random table above is
    # unaffected by adding it.
    traj = _pentagon_traj()

    def c_rows(key, rows=None):
        return "\n".join(
            "    { %s }," % ", ".join("%11d" % fx(v) for v in r[key])
            for r in (recs if rows is None else rows))

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("""/*
 * GENERATED FILE - do not edit.  Produced by model/gen_vectors.py.
 *
 * Pose set for the bare-metal timing harness.  All values are signed Q16.16,
 * the format the kernels' AXI4-Lite registers expect.
 *
 * Three groups, in halves and quarters:
 *
 *   well-conditioned              ordinary poses, comfortable for both solvers
 *   wrist-singular/elbow-extended geometric stress, degrades the ANALYTIC
 *                                 solver's branch selection
 *   dls-tail                      selected on the model's DLS iteration count -
 *                                 a different property; geometric stress alone
 *                                 (seeded within 0.25 rad) converges in 3-5
 *                                 iters every time (hw max/med=1.77 vs. 25
 *                                 iters in the 256-pose host sweep), so this
 *                                 group samples the tail that decides
 *                                 schedulability.
 *
 * ik_traj_* is a SECOND, separate workload - a pentagon tracking path
 * (below), answering a different question than the independent-target pose
 * table above.  The harness runs both.
 */
#ifndef IK_VECTORS_H
#define IK_VECTORS_H

#include <stdint.h>

/* Bumped whenever this file grows something a consumer might require. This
 * header is gitignored (regenerated, not checked in), so a stale copy would
 * otherwise just fail as an undeclared identifier instead of erroring loudly.
 *
 *   1  ik_model_iters_tbl, DLS tail poses
 *   2  ik_qdls_tbl
 *   3  ik_qdls_tbl/ik_model_iters_tbl now from the CLAMPED solver (matches
 *      kernel). Tail POSES stay frozen on the unclamped solver selection.
 *   4  ik_traj_* pentagon trajectory workload. Additive.
 */
#define IK_VECTORS_VERSION 4

#define IK_NVEC %d

/* { x, y, z, roll, pitch, yaw } */
static const int32_t ik_pose_tbl[IK_NVEC][6] = {
%s
};

/* DLS seed: the true solution perturbed, as a tracking loop would supply. */
static const int32_t ik_seed_tbl[IK_NVEC][6] = {
%s
};

/* Analytic (+,+,+) branch solution (matches ik_cfg_tbl). Compare
 * ik_analytic_kernel against this. */
static const int32_t ik_qgold_tbl[IK_NVEC][6] = {
%s
};

/* DLS solution from ik_seed_tbl.  Compare ik_dls_kernel against THIS, never
 * ik_qgold_tbl: IK is multi-valued, so the two solvers legitimately land on
 * different branches (~pi apart in a wrist joint) while both solving the
 * pose to within tolerance. */
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

/* Iterations the double-precision model needs per pose, same lambda/tol as
 * the harness.  A systematic disagreement with the Q16.16 kernel means
 * quantisation moved the convergence path - a finding, not a test failure. */
static const uint8_t ik_model_iters_tbl[IK_NVEC] = {
%s
};

/* ======================================================================
 * Pentagon trajectory - the TRACKING workload.
 *
 * %d samples: closed pentagon, circumradius %.2f m, centred at (%.2f, %.2f)
 * in the base x-y plane, at z = %.2f m, tool pointing straight down, %d
 * samples/edge.  Yaw %s deg at the vertices,
 * interpolated linearly so orientation reaches each vertex smoothly.  ~%.0f mm
 * translation and up to %.0f deg yaw per sample - servo-scale increments.
 *
 * Complements ik_pose_tbl (independent targets, no machine drives an arm
 * that way) with the regime a robot actually runs: the iterative solver
 * seeded from its own solution one control period ago.  ik_traj_seed_tbl is
 * chained accordingly (sample i seeded from sample i-1's DLS solution, sample
 * 0 from the last - closed loop). These are the second lap so sample 0's seed
 * is a real predecessor, not a cold-start outlier.
 *
 * The hardware harness re-chains from the KERNEL's own previous output, not
 * this table's, so fixed-point/double divergence shows as drift against
 * ik_traj_qdls_tbl rather than being masked by resync each sample.  Expect
 * iteration counts nearly flat; see IK_VECTORS_VERSION 4 and main.c.
 * ====================================================================== */
#define IK_NTRAJ %d
#define IK_TRAJ_EDGES 5
#define IK_TRAJ_STEPS %d

/* { x, y, z, roll, pitch, yaw } along the path. */
static const int32_t ik_traj_pose_tbl[IK_NTRAJ][6] = {
%s
};

/* Chained seed: previous sample's DLS solution. NOT a perturbation like
 * ik_seed_tbl. */
static const int32_t ik_traj_seed_tbl[IK_NTRAJ][6] = {
%s
};

/* Analytic (+,+,+) solution per sample. Branch never flips along the path,
 * which is what makes it usable for trajectory following unwrapped. */
static const int32_t ik_traj_qgold_tbl[IK_NTRAJ][6] = {
%s
};

/* DLS solution from the chained seed. Compare ik_dls_kernel against this,
 * never ik_traj_qgold_tbl - same multi-valuedness caveat as ik_qdls_tbl. */
static const int32_t ik_traj_qdls_tbl[IK_NTRAJ][6] = {
%s
};

/* Analytic branch selector (IK_CFG_* bits), constant along the path. */
static const uint8_t ik_traj_cfg_tbl[IK_NTRAJ] = {
%s
};

/* Which pentagon edge each sample lies on, 0..4. */
static const uint8_t ik_traj_edge_tbl[IK_NTRAJ] = {
%s
};

/* 1 at the five vertex samples, 0 along the edges. */
static const uint8_t ik_traj_vtx_tbl[IK_NTRAJ] = {
%s
};

/* Iterations the double-precision model needed per sample, chained. */
static const uint8_t ik_traj_iters_tbl[IK_NTRAJ] = {
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
            "    " + ", ".join("%d" % r["iters"] for r in recs),
            # ---- trajectory block ----
            len(traj), TRAJ_RADIUS, TRAJ_CENTRE[0], TRAJ_CENTRE[1], TRAJ_Z,
            TRAJ_STEPS,
            " - ".join("%g" % y for y in TRAJ_YAW_DEG),
            1000.0 * 2.0 * TRAJ_RADIUS * np.sin(np.pi / 5) / TRAJ_STEPS,
            max(abs(TRAJ_YAW_DEG[(k + 1) % 5] - TRAJ_YAW_DEG[k])
                for k in range(5)) / TRAJ_STEPS,
            len(traj), TRAJ_STEPS,
            c_rows("pose", traj),
            c_rows("seed", traj),
            c_rows("qgold", traj),
            c_rows("qdls", traj),
            "    " + ", ".join("%d" % r["cfg"] for r in traj),
            "    " + ", ".join("%d" % r["edge"] for r in traj),
            "    " + ", ".join("%d" % r["vtx"] for r in traj),
            "    " + ", ".join("%d" % r["iters"] for r in traj)))

    print(f"  {'ik_vectors.h':20s} {len(recs):5d} poses + "
          f"{len(traj)} trajectory samples -> {os.path.relpath(path)}")


def _pentagon_traj():
    """
    Closed pentagon path, sampled as a trajectory with orientation switching
    at the vertices - the regime a real arm runs (seeded from its own
    solution one control period ago), unlike the random stress table above.

    Orientation is (roll, pitch, yaw) = (pi, 0, psi): non-degenerate for
    every psi, unlike a vertical work plane where pitch=+-pi/2 breaks the
    roll/yaw split and the quantised pose stops round-tripping through
    T_to_pose(). psi is TRAJ_YAW_DEG, interpolated linearly per edge so the
    tool never turns discontinuously.

    Expect this table flat (2 iterations, vertices included, confirmed on
    hardware for all 500 samples): DLS behaves as fixed-latency while
    TRACKING; the unbounded tail only shows up when reseeded/disturbed, which
    is what the random table measures instead.

    Seeding is chained (sample i from sample i-1's DLS solution, sample 0
    from the last - closed loop); emitted values are the SECOND lap so
    sample 0 has a genuine predecessor rather than a cold-start outlier.
    """
    cx, cy = TRAJ_CENTRE
    # Vertex 0 at top (+y): symmetric about x axis, pentagon in front of shoulder.
    vtx = [(cx + TRAJ_RADIUS * np.cos(np.pi / 2 + 2 * np.pi * k / 5),
            cy + TRAJ_RADIUS * np.sin(np.pi / 2 + 2 * np.pi * k / 5))
           for k in range(5)]

    poses, edges, is_vtx = [], [], []
    for k in range(5):
        x0, y0 = vtx[k]
        x1, y1 = vtx[(k + 1) % 5]
        yaw0 = np.deg2rad(TRAJ_YAW_DEG[k])
        yaw1 = np.deg2rad(TRAJ_YAW_DEG[(k + 1) % 5])
        for j in range(TRAJ_STEPS):
            t = j / TRAJ_STEPS
            poses.append(q_arr([x0 + t * (x1 - x0), y0 + t * (y1 - y0),
                                TRAJ_Z, np.pi, 0.0,
                                yaw0 + t * (yaw1 - yaw0)]))
            edges.append(k)
            is_vtx.append(1 if j == 0 else 0)

    Ts = [_pose_to_T(p) for p in poses]
    qgold = []
    for i, T in enumerate(Ts):
        qs, ok = M.ik_analytic(T, +1, +1, +1)
        if not ok:
            raise RuntimeError(
                f"pentagon sample {i} at {poses[i][:3]} is out of reach - "
                "adjust TRAJ_CENTRE/TRAJ_RADIUS/TRAJ_Z")
        qgold.append(qs)

    # Two laps: lap 0 warms the chain up, lap 1 is what gets emitted.
    q = qgold[0]
    for _lap in range(2):
        seeds, qdls, iters = [], [], []
        for i, T in enumerate(Ts):
            seeds.append(q)
            qd, it, _, conv = M.ik_dls(T, q, lam=LAMBDA, max_iter=MAX_ITER,
                                       tol=TOL)
            if not conv:
                raise RuntimeError(
                    f"pentagon sample {i} did not converge in {MAX_ITER} "
                    "iterations from its predecessor's solution")
            qdls.append(q_arr(qd))
            iters.append(int(it))
            q = qd

    return [{"pose": poses[i], "seed": q_arr(seeds[i]), "qgold": qgold[i],
             "qdls": qdls[i], "iters": iters[i], "cfg": 1 | 2 | 4,
             "edge": edges[i], "vtx": is_vtx[i]}
            for i in range(len(poses))]


def _pose_to_T(pose):
    Tq = np.eye(4)
    Tq[:3, :3] = M.rpy_to_rot(*pose[3:])
    Tq[:3, 3] = pose[:3]
    return Tq


def _find_dls_tail(count, floor=8, ceiling=MAX_ITER // 2,
                   pool_target=30, tries_max=60000):
    """
    Poses whose DLS iteration count lands in the tail.

    Seeds drawn far (+-1.2 rad), not _mk_rec()'s +-0.25: seed distance is the
    dominant term in iteration count, so a tight perturbation can't produce a
    tail. Selection runs the UNCLAMPED solver, frozen there deliberately (see
    call site below). Returns a spread across the tail, not `count` copies of
    the worst case - the real-time argument needs the whole distribution.
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
        # step_max=0: UNCLAMPED solver, frozen deliberately - if selection
        # tracked solver improvements it would just pick harder poses and
        # report no improvement. Also keeps this comparable to pre-trust-
        # region hardware runs.
        _, iters, _, conv = M.ik_dls(Tq, seed, lam=LAMBDA, max_iter=MAX_ITER,
                                     tol=TOL, step_max=0.0)
        # Margin below the cap: a pose landing exactly on MAX_ITER disagrees
        # with the kernel about convergence (model rechecks residual after
        # the loop, kernel doesn't).
        if conv and floor <= iters <= ceiling:
            pool.append((iters, q, T, seed))

    if len(pool) <= count:
        return pool

    pool.sort(key=lambda r: r[0])
    # Evenly spaced so the slowest pose found is always taken.
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
    # qdls != qgold by design: IK is multi-valued, DLS converges to whichever
    # branch its seed is nearest, routinely ~pi apart in a wrist joint from
    # qgold's fixed (+,+,+) branch while both solve the pose exactly.
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
