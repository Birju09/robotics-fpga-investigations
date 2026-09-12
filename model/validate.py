#!/usr/bin/env python3
"""Validate the golden model: FK/IK round-trips, Jacobian, and DLS convergence."""

import numpy as np
import collision as C
import ik_model as M

rng = np.random.default_rng(0xC0FFEE)

# Joint limits used for random sampling (radians).
QLIM = np.array([
    [-2.9, 2.9],
    [-1.9, 1.9],
    [-2.6, 2.6],
    [-2.9, 2.9],
    [-1.9, 1.9],
    [-2.9, 2.9],
])


def rand_q(n=1):
    return rng.uniform(QLIM[:, 0], QLIM[:, 1], size=(n, 6))


def test_analytic_roundtrip(n=4000):
    """FK(IK(FK(q))) must reproduce the original pose for every branch returned."""
    worst = 0.0
    checked = 0
    fails = []
    for q in rand_q(n):
        T = M.fk(q)
        sols = M.ik_analytic_all(T)
        if not sols:
            fails.append(("no solution", q))
            continue
        for cfg, qs in sols:
            T2 = M.fk(qs)
            err = np.abs(T2 - T).max()
            worst = max(worst, err)
            checked += 1
            if err > 1e-8:
                fails.append((cfg, q, err))
    return worst, checked, fails


def test_analytic_recovers_original(n=4000):
    """At least one branch must recover the *joint vector* we started from."""
    worst = 0.0
    misses = 0
    for q in rand_q(n):
        T = M.fk(q)
        best = np.inf
        for _cfg, qs in M.ik_analytic_all(T):
            d = np.abs(M.wrap_pi(qs - q)).max()
            best = min(best, d)
        if best > 1e-7:
            misses += 1
        worst = max(worst, min(best, 1e3))
    return worst, misses


def test_jacobian(n=500, h=1e-7):
    """Compare the analytic Jacobian against a numerical differentiation."""
    worst = 0.0
    for q in rand_q(n):
        J = M.jacobian(q)
        Jn = np.zeros((6, 6))
        T0 = M.fk(q)
        for i in range(6):
            qp = q.copy()
            qp[i] += h
            T1 = M.fk(qp)
            Jn[:3, i] = (T1[:3, 3] - T0[:3, 3]) / h
            dR = (T1[:3, :3] - T0[:3, :3]) / h @ T0[:3, :3].T
            Jn[3:, i] = [dR[2, 1], dR[0, 2], dR[1, 0]]
        worst = max(worst, np.abs(J - Jn).max())
    return worst


def test_dls(n=600, lam=0.08, max_iter=64, tol=1e-5):
    """Seed DLS near the true solution and record convergence statistics."""
    iters = []
    conv = 0
    for q in rand_q(n):
        T = M.fk(q)
        q0 = q + rng.uniform(-0.35, 0.35, size=6)
        qs, k, err, ok = M.ik_dls(T, q0, lam=lam, max_iter=max_iter, tol=tol)
        iters.append(k)
        conv += int(ok)
    return np.array(iters), conv, n


def test_dls_near_singular(n=300, lam=0.08, max_iter=64, tol=1e-5):
    """Same, but seeded at wrist/elbow singularities where DLS is stressed."""
    iters = []
    conv = 0
    for _ in range(n):
        q = rng.uniform(QLIM[:, 0], QLIM[:, 1], size=6)
        q[4] = rng.choice([0.0, 0.0, 1e-3, -1e-3])   # wrist singularity
        T = M.fk(q)
        q0 = q + rng.uniform(-0.3, 0.3, size=6)
        q0[4] = 1e-4
        qs, k, err, ok = M.ik_dls(T, q0, lam=lam, max_iter=max_iter, tol=tol)
        iters.append(k)
        conv += int(ok)
    return np.array(iters), conv, n


def test_seg_dist_bruteforce(n=2000, samples=400):
    """
    Segment-segment distance against brute force, and against the distance-only
    form already in ik_model.

    Brute force samples both segments on a grid, so it is an UPPER bound on
    the true minimum: the closed form must never exceed it, and must come
    within the grid's own resolution of it. Checking only |closed - brute|
    would pass a closed form that returns something too large.

    Degenerate cases are included deliberately - zero-length segments
    (coincident DH frames, which this arm genuinely has at the wrist),
    exactly parallel segments, and exactly touching segments. Every
    segment-distance implementation breaks there, and the fixed-point build
    will break differently again.
    """
    worst_over = 0.0     # closed form above brute force: a real error
    worst_gap = 0.0      # brute force above closed form: grid resolution
    worst_agree = 0.0    # vs ik_model._seg_dist
    worst_witness = 0.0  # |wa - wb| vs the returned distance
    u = np.linspace(0.0, 1.0, samples)

    cases = []
    for _ in range(n):
        cases.append(rng.uniform(-1.0, 1.0, size=(4, 3)))
    # zero-length (point vs segment), and point vs point
    for _ in range(60):
        p = rng.uniform(-1.0, 1.0, size=3)
        q = rng.uniform(-1.0, 1.0, size=(2, 3))
        cases.append(np.array([p, p, q[0], q[1]]))
        r = rng.uniform(-1.0, 1.0, size=3)
        cases.append(np.array([p, p, r, r]))
    # exactly parallel, and collinear
    for _ in range(60):
        p0 = rng.uniform(-1.0, 1.0, size=3)
        d = rng.uniform(-1.0, 1.0, size=3)
        off = rng.uniform(-1.0, 1.0, size=3)
        cases.append(np.array([p0, p0 + d, p0 + off, p0 + off + d]))
        cases.append(np.array([p0, p0 + d, p0 + 2.0 * d, p0 + 3.0 * d]))
    # NEAR-parallel, swept across the whole range of small angles.
    #
    # This block is a regression guard, not padding. The first version of
    # coll::seg_seg() declared |a e - b b| < IK_PIVOT_EPS "parallel" and
    # forced s = 0. Since a e - b b = a e sin^2(theta), that threshold sits
    # around 3 degrees off parallel, and the forced s reported clearances up
    # to 21 mm LARGER than the truth - a solver would be told it was clear
    # when it was not. Uniformly random segment pairs are almost never that
    # close to parallel, so nothing else here caught it.
    #
    # worst_over is what fails: any threshold-based shortcut in the clamp
    # cascade shows up as the closed form exceeding brute force.
    for _ in range(240):
        p0 = rng.uniform(-1.0, 1.0, size=3)
        d = rng.uniform(-1.0, 1.0, size=3)
        ang = 10.0 ** rng.uniform(-6.0, -1.0)      # 1e-6 .. 0.1 rad
        perp = np.cross(d, rng.uniform(-1.0, 1.0, size=3))
        perp /= max(np.linalg.norm(perp), 1e-30)
        d2 = d + ang * np.linalg.norm(d) * perp
        q0 = p0 + rng.uniform(-1.0, 1.0, size=3) * 0.3
        # both overlapping and offset along the shared direction, since the
        # failure depends on where the true closest point falls
        cases.append(np.array([p0, p0 + d, q0, q0 + d2]))
        cases.append(np.array([p0, p0 + d, q0 + 0.9 * d, q0 + 0.9 * d + d2]))
    # exactly touching, at an endpoint and in the interior
    for _ in range(60):
        p0, p1, q1 = rng.uniform(-1.0, 1.0, size=(3, 3))
        cases.append(np.array([p0, p1, p1, q1]))
        mid = 0.5 * (p0 + p1)
        cases.append(np.array([p0, p1, mid, mid + rng.uniform(-1, 1, size=3)]))

    for p0, p1, q0, q1 in cases:
        d, wa, wb, s, t = C.seg_seg_witness(p0, p1, q0, q1)

        A = p0 + u[:, None] * (p1 - p0)
        B = q0 + u[:, None] * (q1 - q0)
        brute = float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2)
                              .sum(-1)).min())

        worst_over = max(worst_over, d - brute)
        worst_gap = max(worst_gap, brute - d)
        worst_agree = max(worst_agree,
                          abs(d - M._seg_dist(p0, p1, q0, q1)))
        worst_witness = max(worst_witness,
                            abs(float(np.linalg.norm(wa - wb)) - d))
        assert -1e-12 <= s <= 1 + 1e-12 and -1e-12 <= t <= 1 + 1e-12

    return worst_over, worst_gap, worst_agree, worst_witness, len(cases)


def test_clearance_gradient(n=400, h=1e-6, tie=1e-3):
    """
    Analytic d(d_min)/dq against numerical differentiation.

    d_min is a min over capsule pairs, so it is not differentiable where two
    pairs tie, nor where a witness point crosses a segment endpoint. Samples
    whose achieving pair changes under the perturbation are skipped rather
    than counted as failures - the non-smoothness is real, and it is why the
    repulsion this feeds is threshold-gated and trust-region clamped instead
    of being trusted as an exact descent direction. The skip count is
    reported so it cannot quietly become the whole sample.
    """
    caps = C.capsules()
    pairs = C.capsule_pairs(caps)
    worst = 0.0
    checked = skipped = 0

    for q in rand_q(n):
        d0, g, info = C.clearance_gradient(q, caps, pairs)
        if info is None or info["normal"] is None:
            skipped += 1
            continue

        gn = np.zeros(6)
        ok = True
        for i in range(6):
            qp, qm = q.copy(), q.copy()
            qp[i] += h
            qm[i] -= h
            dp, ip = C.min_clearance(qp, caps, pairs)
            dm, im = C.min_clearance(qm, caps, pairs)
            # same achieving pair either side, and no near-tie
            if ip["pair"] != info["pair"] or im["pair"] != info["pair"]:
                ok = False
                break
            gn[i] = (dp - dm) / (2.0 * h)
        if not ok:
            skipped += 1
            continue

        # A near-tie makes the one-sided pair test pass but the derivative
        # still kink, so check the runner-up margin explicitly.
        second = np.inf
        Ts = M.fk_all(q)
        for (ia, ib) in pairs:
            if (ia, ib) == info["pair"]:
                continue
            a0, a1 = C.capsule_endpoints(Ts, caps[ia])
            b0, b1 = C.capsule_endpoints(Ts, caps[ib])
            sd, _, _, _, _ = C.seg_seg_witness(a0, a1, b0, b1)
            second = min(second, sd - caps[ia].radius - caps[ib].radius)
        if second - d0 < tie:
            skipped += 1
            continue

        worst = max(worst, float(np.abs(g - gn).max()))
        checked += 1

    return worst, checked, skipped


def test_nullspace_projector(n=300, lam=0.02):
    """
    The damped nullspace projector, checked against what damping guarantees
    rather than against zero.

    N = I - J_t^T (J_t J_t^T + lam^2 I)^-1 J_t is not an exact projector, so
    J_t N is not exactly zero. It is exactly

        J_t N = lam^2 (J_t J_t^T + lam^2 I)^-1 J_t

    which is the bound asserted here, and asserting THAT rather than
    |J_t N z| < eps is the point: it is the identity that fails if the
    projection is built wrong.

    At task_dim = 6 the same algebra says N's eigenvalues are
    lam^2/(sigma_i^2 + lam^2), so a full-pose task has near-zero avoidance
    authority away from singularities. That is the control case - it is
    measured here, not assumed.
    """
    out = {}
    for td in (M.TASK_POS, M.TASK_FULL):
        worst_leak = 0.0
        ratios = []
        for q in rand_q(n):
            J = M.jacobian(q)[M.task_rows(td)]
            A = J @ J.T + lam * lam * np.eye(td)
            z = rng.normal(size=6)
            Nz = C.nullspace_project(q, z, lam=lam, task_dim=td)

            # identity: J N z == lam^2 A^-1 J z
            lhs = J @ Nz
            rhs = lam * lam * np.linalg.solve(A, J @ z)
            worst_leak = max(worst_leak, float(np.abs(lhs - rhs).max()))

            ratios.append(float(np.linalg.norm(Nz) / np.linalg.norm(z)))

        sig_bound = []
        for q in rand_q(60):
            J = M.jacobian(q)[M.task_rows(td)]
            s = np.linalg.svd(J, compute_uv=False)
            sig_bound.append(lam * lam / (s.min() ** 2 + lam * lam))

        out[td] = {
            "identity": worst_leak,
            "ratio_med": float(np.median(ratios)),
            "ratio_max": float(np.max(ratios)),
            "eig_bound_med": float(np.median(sig_bound)),
        }
    return out


def test_fixed_point_headroom(n=3000):
    """Check the dynamic range that Q16.16 has to cover."""
    mx = {"J": 0.0, "JJt": 0.0, "inv": 0.0, "pos": 0.0}
    for q in rand_q(n):
        J = M.jacobian(q)
        A = J @ J.T + 0.08 ** 2 * np.eye(6)
        Ai = np.linalg.inv(A)
        mx["J"] = max(mx["J"], np.abs(J).max())
        mx["JJt"] = max(mx["JJt"], np.abs(A).max())
        mx["inv"] = max(mx["inv"], np.abs(Ai).max())
        mx["pos"] = max(mx["pos"], np.abs(M.fk(q)[:3, 3]).max())
    return mx


if __name__ == "__main__":
    print("=" * 70)
    print("GOLDEN MODEL VALIDATION")
    print("=" * 70)

    w = test_jacobian()
    print(f"\n[Jacobian]  max |analytic - numeric| = {w:.3e}   {'PASS' if w < 1e-5 else 'FAIL'}")

    worst, checked, fails = test_analytic_roundtrip()
    print(f"\n[Analytic IK / pose round-trip]")
    print(f"  branches checked      : {checked}")
    print(f"  max |FK(IK(T)) - T|   : {worst:.3e}")
    print(f"  failures              : {len(fails)}   {'PASS' if not fails else 'FAIL'}")
    for f in fails[:5]:
        print("   ", f)

    worst, misses = test_analytic_recovers_original()
    print(f"\n[Analytic IK / joint recovery]")
    print(f"  max best-branch joint error : {worst:.3e}")
    print(f"  poses with no exact branch  : {misses}   {'PASS' if misses == 0 else 'FAIL'}")

    it, conv, n = test_dls()
    print(f"\n[DLS / well-conditioned seeds]")
    print(f"  converged      : {conv}/{n}")
    print(f"  iterations     : min={it.min()} median={int(np.median(it))} "
          f"p95={int(np.percentile(it,95))} max={it.max()}")

    it, conv, n = test_dls_near_singular()
    print(f"\n[DLS / near-singular seeds]")
    print(f"  converged      : {conv}/{n}")
    print(f"  iterations     : min={it.min()} median={int(np.median(it))} "
          f"p95={int(np.percentile(it,95))} max={it.max()}")

    over, gap, agree, wit, ncase = test_seg_dist_bruteforce()
    print(f"\n[Capsule distance]  {ncase} segment pairs "
          f"(incl. degenerate: zero-length, parallel, touching)")
    # Bound is one Q16.16 LSB (1.53e-5), not zero. The model's own
    # denominator guard (eps = 1e-12) does bite at the extreme end of the
    # near-parallel sweep, around 1e-6 rad, and overestimates by ~5e-7 m -
    # sub-micron, 0.03 of an LSB, invisible to the kernel this validates.
    # The failure this test exists to catch is three orders of magnitude
    # larger: a coarse threshold in the clamp cascade overestimated by 21 mm.
    # Bounding at an LSB catches that and does not flag arithmetic that
    # cannot reach the fixed-point result.
    LSB = 1.0 / 65536.0
    print(f"  closed form above brute force : {over:.3e}   "
          f"{'PASS' if over < LSB else 'FAIL'}  (bound 1 LSB = {LSB:.2e})")
    print(f"  brute force above closed form : {gap:.3e}   "
          f"(grid resolution, not an error)")
    print(f"  vs ik_model._seg_dist         : {agree:.3e}   "
          f"{'PASS' if agree < 1e-12 else 'FAIL'}")
    print(f"  |wa - wb| vs returned dist    : {wit:.3e}   "
          f"{'PASS' if wit < 1e-12 else 'FAIL'}")

    caps = C.capsules()
    pairs, dropped = C.capsule_pairs(caps, explain=True)
    w, checked, skipped = test_clearance_gradient()
    print(f"\n[Clearance gradient]  {len(caps)} capsules, "
          f"{len(pairs)} pairs checked, {len(dropped)} pruned")
    print(f"  max |analytic - numeric|      : {w:.3e}   "
          f"{'PASS' if w < 1e-5 else 'FAIL'}")
    print(f"  samples used / skipped        : {checked} / {skipped}   "
          f"{'PASS' if checked > 0.5 * (checked + skipped) else 'FAIL'}")
    print("  (skips are pair ties and witness-point endpoint crossings -")
    print("   genuinely non-differentiable, not implementation failures)")
    print("  NOTE: LINK_RADIUS is chosen, not cited. A negative clearance is")
    print("        a finding; a positive one is not a proof of clearance.")

    ns = test_nullspace_projector()
    print(f"\n[Nullspace projector]  lam = 0.02")
    print("  task_dim   |J N z - lam^2 A^-1 J z|   |Nz|/|z| med    max     "
          "eig bound med")
    for td, r in sorted(ns.items()):
        print(f"  {td:^8}   {r['identity']:^23.3e}   {r['ratio_med']:10.3e} "
              f"{r['ratio_max']:10.3e}   {r['eig_bound_med']:.3e}")
    ns_ok = all(r["identity"] < 1e-9 for r in ns.values())
    print(f"  projector identity holds      : {'PASS' if ns_ok else 'FAIL'}")
    print("  At task_dim=6 the ratio is the damped leak, not a nullspace:")
    print("  near-zero avoidance authority away from singularities, as")
    print("  predicted. At task_dim=3 the nullspace is genuine (dim 3).")

    mx = test_fixed_point_headroom()
    print(f"\n[Q16.16 headroom]  (integer part must hold these)")
    for k, v in mx.items():
        print(f"  max |{k:4s}| = {v:9.4f}")
    print()
