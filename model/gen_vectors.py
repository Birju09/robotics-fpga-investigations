#!/usr/bin/env python3
"""
Emit golden test vectors consumed by the HLS C++ testbenches.

Files written to hls/tb/vectors/ (plain text, one record per line, Q16.16 as
signed decimal so the C++ side can read them with a plain istream):

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

# The arm's own joint limits, from the robot definition rather than restated.
# These are asymmetric (the elbow cannot swing equally both ways) and they do
# the job a collision model would otherwise have to: a2 and L3 differ by less
# than a millimetre on this arm, so bare link lengths would let the forearm
# fold flat onto the upper arm, and only the q3 range forbids it.
QLIM = M.ROBOT.qlim

LAMBDA = 0.02
TOL = 1e-3
MAX_ITER = 64

# ---------------------------------------------------------------------------
# Pentagon trajectory (see _pentagon_traj)
# ---------------------------------------------------------------------------
# Sited by search over the reachable set rather than picked: the centre,
# height and radius below are the widest pentagon whose every sample keeps the
# elbow conditioning |sin gamma| above 0.45, so no part of the path is near the
# elbow singularity and the workload measures tracking rather than
# ill-conditioning.  Re-run model/plot_dh.py after changing these.
TRAJ_CENTRE = (0.40, 0.0)     # metres, in the base x-y plane
TRAJ_Z = 0.15                 # tool-tip height above the base plane
TRAJ_RADIUS = 0.20            # circumradius of the pentagon
TRAJ_STEPS = 100              # samples per edge -> 5 * 100 = 500 poses. Roughly 2.35mm per step
# End-effector orientation at each vertex, degrees of yaw about the (vertical)
# tool axis.  Interpolated linearly along each edge, so the tool reaches the
# vertex already at that vertex's orientation and never turns discontinuously.
TRAJ_YAW_DEG = (0.0, 90.0, 0.0, 90.0, 0.0)


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
        # Two conditioning numbers, both letting the testbench tell an
        # ill-conditioned pose from a wrong answer.  The elbow one is
        # |sin(gamma)|; the shoulder one is root/rho, which only exists
        # because the arm has a lateral offset - see the docstrings in
        # ik_model.elbow_conditioning() and shoulder_conditioning().
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
            # Push the wrist centre close to full extension.  sin(gamma) = 0
            # is theta3 = -PHI on this arm, not +PHI: alpha3 is negative, so
            # theta3 = gamma - PHI.  Verified against elbow_conditioning()
            # rather than assumed - at +PHI it returns 0.094, at -PHI exactly 0.
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

    # The trajectory table.  Deterministic - it takes nothing from `rng`, so
    # adding it left every value in the random table above bit-identical.
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
 *
 * ik_traj_* is a SECOND, separate workload - a pentagon path, described at
 * its own tables below.  The pose table above is a stress workload of
 * independent targets, which is not what a robot does; the trajectory table
 * is the tracking regime, and the two answer different questions.  Neither
 * replaces the other and the harness runs both.
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
 *   3  DLS trust region: ik_qdls_tbl and ik_model_iters_tbl now come from the
 *      CLAMPED solver, matching the kernel.  The tail POSES are unchanged -
 *      their selection is deliberately frozen on the unclamped solver so the
 *      workload stays fixed across solver changes.
 *   4  ik_traj_* - the pentagon trajectory workload.  Additive; every value
 *      in the tables above is unchanged.
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

/* ======================================================================
 * Pentagon trajectory - the TRACKING workload.
 *
 * %d samples: a closed pentagon of circumradius %.2f m centred at
 * (%.2f, %.2f) in the base x-y plane, traversed at z = %.2f m with the tool
 * approach axis pointing straight down, %d samples per edge.  Tool yaw is
 * %s deg at the five vertices and is interpolated
 * linearly along each edge, so the orientation reaches each vertex smoothly
 * rather than switching there.  Per sample that is ~%.0f mm of translation
 * and up to %.0f deg of yaw - servo-scale increments.
 *
 * Why this exists next to ik_pose_tbl: that table is independent targets,
 * each seeded from its own perturbation, and no machine drives an arm that
 * way.  A robot follows a path, and the iterative solver is seeded from the
 * solution it produced one control period ago.  Whether DLS is schedulable
 * in a servo loop is a question about THAT regime, and it has to be measured
 * rather than inferred from the random one.
 *
 * ik_traj_seed_tbl is chained: sample i is seeded from sample i-1's DLS
 * solution, and sample 0 from the last sample's, because the path is a closed
 * loop.  These are the second lap, so sample 0's seed is a real predecessor
 * rather than a cold start showing up as a lone outlier.
 *
 * The harness re-chains this on hardware (it feeds the KERNEL's previous
 * output forward, not this table's) - so a divergence between the fixed-point
 * and double chains shows up as drift against ik_traj_qdls_tbl rather than
 * being masked by resynchronising the seed every sample.  Compare per-sample
 * and expect the iteration counts to be nearly flat; see IK_VECTORS_VERSION 4
 * and main.c's trajectory section.
 * ====================================================================== */
#define IK_NTRAJ %d
#define IK_TRAJ_EDGES 5
#define IK_TRAJ_STEPS %d

/* { x, y, z, roll, pitch, yaw } along the path. */
static const int32_t ik_traj_pose_tbl[IK_NTRAJ][6] = {
%s
};

/* Chained seed: the previous sample's DLS solution, as a tracking loop
 * supplies.  NOT a perturbation of the answer like ik_seed_tbl. */
static const int32_t ik_traj_seed_tbl[IK_NTRAJ][6] = {
%s
};

/* Analytic (+,+,+) solution per sample.  Continuous along the whole path -
 * the branch does not flip anywhere on it, which is what makes the analytic
 * solver usable for trajectory following without post-hoc unwrapping. */
static const int32_t ik_traj_qgold_tbl[IK_NTRAJ][6] = {
%s
};

/* The model's DLS solution from the chained seed.  Compare ik_dls_kernel
 * against this, never against ik_traj_qgold_tbl - same multi-valuedness
 * caveat as ik_qdls_tbl above. */
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
    A closed pentagon path, sampled as a trajectory, with the tool orientation
    switching at the vertices.

    This is the table a robotics application would actually generate.  The
    random pose set above is a *stress* workload - independent targets, each
    seeded from its own perturbation - and nothing in a real machine looks
    like that.  A real arm follows a path: consecutive targets are millimetres
    apart, and the iterative solver is seeded from the solution it produced
    one control period ago.  Whether DLS can be scheduled in a servo loop
    depends on that regime, not on the random one, so it has to be measured
    separately rather than inferred.

    Geometry: the tool tip traverses the five edges in the horizontal plane
    z = TRAJ_Z with the approach axis pointing straight down, which puts the
    pose orientation at (roll, pitch, yaw) = (pi, 0, psi) - a non-degenerate
    RPY triple for every psi, unlike a vertical work plane, which lands on
    pitch = +-pi/2 where the roll/yaw split stops being defined and the
    quantised pose no longer round-trips through T_to_pose().

    psi is TRAJ_YAW_DEG at each vertex - 0, 90, 0, 90, 0 degrees - and is
    interpolated linearly along the edge between them, so the tool arrives at
    each vertex already carrying that vertex's orientation and never turns
    discontinuously.  At TRAJ_STEPS = 100 that is 0.9 degrees of yaw and about
    2.4 mm of translation per sample, both of them servo-scale increments.

    Expect this table to be flat: every pose converges in two iterations,
    vertices included, because a smooth path warm-started from its own
    predecessor never presents the solver with a large residual.  Hardware
    agrees at every one of the 500 samples.  That is
    the result, not a defect in the workload - it says the DLS solver behaves
    as a fixed-latency block while it is TRACKING, and that its unbounded
    iteration count only becomes a scheduling problem when the loop is
    disturbed or reseeded, which is what the random table measures.  The two
    tables are answering different questions and both are reported.

    Seeding is chained - each pose is seeded from the previous pose's DLS
    solution, and the first from the last, because the path is a closed loop.
    The tables are the SECOND lap, so the first pose's seed is a genuine
    predecessor solution rather than a cold start that would show up as a
    one-sample outlier in a 50-sample summary.
    """
    cx, cy = TRAJ_CENTRE
    # Vertex 0 at the top (+y) so the path is symmetric about the x axis and
    # the whole pentagon stays in front of the shoulder.
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
    Poses whose DLS iteration count actually lands in the tail.

    Seeds are drawn far (+-1.2 rad) rather than at _mk_rec()'s +-0.25: the
    seed distance is the dominant term in how many iterations DLS needs, so a
    tight perturbation cannot produce a tail no matter which pose it is
    applied to.

    Selection runs the UNCLAMPED solver (step_max=0) and is frozen there on
    purpose - see the comment at the call below.

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
        # step_max=0: the UNCLAMPED solver, deliberately, even though the
        # kernel now clamps.  The selection criterion has to be frozen or the
        # experiment eats itself - improve the solver, and a tail selected with
        # the improved solver simply picks harder poses until the iteration
        # counts look the same as before, reporting no improvement no matter
        # how large one is.  Holding the workload fixed is what makes a solver
        # change measurable.  (It also keeps this table comparable with the
        # hardware runs taken before the trust region existed.)
        _, iters, _, conv = M.ik_dls(Tq, seed, lam=LAMBDA, max_iter=MAX_ITER,
                                     tol=TOL, step_max=0.0)
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
