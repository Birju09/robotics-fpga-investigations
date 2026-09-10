#!/usr/bin/env python3
"""Validate the golden model: FK/IK round-trips, Jacobian, and DLS convergence."""

import numpy as np
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

    mx = test_fixed_point_headroom()
    print(f"\n[Q16.16 headroom]  (integer part must hold these)")
    for k, v in mx.items():
        print(f"  max |{k:4s}| = {v:9.4f}")
    print()
