#!/usr/bin/env python3
"""
The experiment behind the DLS trust region, and the decision NOT to
implement adaptive lambda.  Results recorded in the trust region block of
hls/include/ik_config.hpp; this reproduces them.

    python3 model/sweep_dls.py              # the headline table
    python3 model/sweep_dls.py --radius     # trust radius sweep
    python3 model/sweep_dls.py --lm         # adaptive lambda, both forms
    python3 model/sweep_dls.py --lambda     # lambda itself
    python3 model/sweep_dls.py --all

Three things a naive sweep would get wrong, each of which reversed a
conclusion while this was being built:

1. Counts LOOP PASSES, not accepted steps - what the kernel's `iters`
   register reports and what 15,952 ns/iteration multiplies. A rejected LM
   step costs a full iteration here, which is most of why LM loses.

2. Evaluates on a FROZEN pose table: gen_vectors.py selects a quarter of the
   harness table on the model's DLS iteration count, so an unfrozen
   selection would just hand an improved solver harder poses and report no
   improvement. build_table() freezes on the unclamped solver, matching
   gen_vectors.py.

3. Reports median next to max/med, since max/med alone can be minimized by
   making every pose equally slow (see the lambda sweep below).
"""

import sys

import numpy as np

import ik_model as M

I6 = np.eye(6)
QLIM = np.array([
    [-2.9, 2.9], [-1.9, 1.9], [-2.6, 2.6],
    [-2.9, 2.9], [-1.9, 1.9], [-2.9, 2.9],
])

LAMBDA, TOL, MAX_ITER = 0.02, 1e-3, 64

# Measured on hardware at 80 MHz, 48 poses: latency_ns = 230 + 15952 * iters,
# per-iteration cost constant to 0.4%.  See STATUS.md.
NS_PER_ITER, NS_FIXED = 15952, 230


def dls(T_des, q0, lam=LAMBDA, max_iter=MAX_ITER, tol=TOL, step_max=0.0,
        halve=True, halvings=M.DLS_STEP_HALVINGS, lm=None, lm_up=2.0,
        lm_down=2.0, lam_max=64.0):
    """
    One configurable DLS.  Returns (passes, converged, rejects).

    lm:  None       fixed lambda (what the kernel does)
         'forward'  raise lambda when the residual grew, but keep the step
         'backtrack' textbook LM: undo the step, raise lambda, retry

    halve=True reproduces the kernel's power-of-two clamp; halve=False is an
    exact step_max/|dq| scale, the baseline the clamp is measured against.
    """
    q = np.array(q0, dtype=float)
    lam_sq, lam_sq_min, lam_sq_max = lam * lam, lam * lam, lam_max * lam_max
    tol_sq = tol * tol
    qp = ep = Jp = None
    errp = np.inf
    rejects = 0

    for k in range(1, max_iter + 1):
        e = M.pose_error(T_des, M.fk(q))
        err_sq = float(e @ e)
        if err_sq < tol_sq:
            return k, True, rejects
        J = M.jacobian(q)

        if lm == 'backtrack':
            if qp is not None and err_sq > errp:
                q, e, J, err_sq = qp.copy(), ep.copy(), Jp.copy(), errp
                lam_sq = min(lam_sq * lm_up * lm_up, lam_sq_max)
                rejects += 1
            else:
                lam_sq = max(lam_sq / (lm_down * lm_down), lam_sq_min)
                qp, ep, Jp, errp = q.copy(), e.copy(), J.copy(), err_sq
        elif lm == 'forward':
            if err_sq > errp:
                lam_sq = min(lam_sq * lm_up * lm_up, lam_sq_max)
                rejects += 1
            else:
                lam_sq = max(lam_sq / (lm_down * lm_down), lam_sq_min)
            errp = err_sq

        dq = J.T @ np.linalg.solve(J @ J.T + lam_sq * I6, e)
        if step_max > 0.0:
            if halve:
                dq = M.clamp_step(dq, step_max, halvings)
            else:
                m = float(np.max(np.abs(dq)))
                if m > step_max:
                    dq = dq * (step_max / m)
        q = M.wrap_pi(q + dq)

    e = M.pose_error(T_des, M.fk(q))
    return max_iter, float(e @ e) < tol_sq, rejects


# ---------------------------------------------------------------------------
# Workloads
# ---------------------------------------------------------------------------
def _q(x):
    return np.array([M.quantize(v) for v in np.asarray(x).reshape(-1)])


def make_cases(n, spread, seed):
    """n (T_des, q_seed) pairs, seeded +-spread rad from a true solution."""
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        q = _q(rng.uniform(QLIM[:, 0], QLIM[:, 1]))
        T = M.fk(q)
        pose = _q(M.T_to_pose(T))
        Tq = M.pose_to_T(pose)
        if not M.ik_analytic(Tq)[1]:
            continue
        out.append((Tq, _q(q + rng.uniform(-spread, spread, 6))))
    return out


def build_table(n=48, seed=0x1153):
    """
    The harness pose table, same composition as gen_vectors.gen_c_header():
    half ordinary, a quarter geometrically stressed, a quarter selected on
    DLS iteration count with far seeds. Tail selection runs the UNCLAMPED
    solver and stays there (module docstring, point 2).
    """
    rng = np.random.default_rng(seed)
    recs = []

    def mk(q, T, seed_=None):
        pose = _q(M.T_to_pose(T))
        if seed_ is None:
            seed_ = _q(q + rng.uniform(-0.25, 0.25, 6))
        return (M.pose_to_T(pose), seed_)

    while len(recs) < n // 2:
        q = _q(rng.uniform(QLIM[:, 0], QLIM[:, 1]))
        T = M.fk(q)
        if M.elbow_conditioning(T, +1) < 0.25:
            continue
        recs.append(mk(q, T))

    while len(recs) < (3 * n) // 4:
        q = _q(rng.uniform(QLIM[:, 0], QLIM[:, 1]))
        if len(recs) % 2 == 0:
            q[4] = rng.choice([1e-3, -1e-3, 2e-3])
        else:
            q[2] = M.wrap_pi(M.PHI + rng.uniform(-0.02, 0.02))
        q = _q(q)
        T = M.fk(q)
        if not M.ik_analytic(M.pose_to_T(_q(M.T_to_pose(T))))[1]:
            continue
        recs.append(mk(q, T))

    want = n - len(recs)
    pool, tries = [], 0
    while tries < 60000 and len(pool) < 30 * want:
        tries += 1
        q = _q(rng.uniform(QLIM[:, 0], QLIM[:, 1]))
        T = M.fk(q)
        Tq = M.pose_to_T(_q(M.T_to_pose(T)))
        if not M.ik_analytic(Tq)[1]:
            continue
        s = _q(q + rng.uniform(-1.2, 1.2, 6))
        k, ok, _ = dls(Tq, s, step_max=0.0)
        if ok and 8 <= k <= 32:
            pool.append((k, Tq, s))
    pool.sort(key=lambda r: r[0])
    idx = sorted({round(i * (len(pool) - 1) / (want - 1)) for i in range(want)})
    recs.extend((Tq, s) for _, Tq, s in (pool[i] for i in idx))
    return recs


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def report(label, cases, variants, latency=False):
    print(f"\n=== {label}  (n={len(cases)}) ===")
    head = (f"{'variant':28s} {'conv':>9s} {'min':>4s} {'med':>4s} {'p95':>4s} "
            f"{'max':>4s} {'mean':>6s} {'max/med':>8s}")
    if latency:
        head += f" | {'med us':>7s} {'max us':>7s}"
    print(head)
    for name, kw in variants:
        res = [dls(T, s, **kw) for T, s in cases]
        it = np.sort([r[0] for r in res])
        nconv = sum(r[1] for r in res)
        med, mx = int(it[len(it) // 2]), int(it[-1])
        line = (f"{name:28s} {nconv:4d}/{len(cases):4d} {it[0]:4d} {med:4d} "
                f"{int(it[int(len(it) * .95)]):4d} {mx:4d} {it.mean():6.2f} "
                f"{mx / max(med, 1):8.2f}")
        if latency:
            line += (f" | {(NS_FIXED + NS_PER_ITER * med) / 1000:7.1f} "
                     f"{(NS_FIXED + NS_PER_ITER * mx) / 1000:7.1f}")
        print(line)


R = M.DLS_STEP_MAX

HEADLINE = [
    ("no clamp (baseline)         ", dict(step_max=0.0)),
    (f"clamp {R} rad  (implemented)", dict(step_max=R)),
    (f"clamp {R} rad, exact scale  ", dict(step_max=R, halve=False)),
    ("LM, raise lambda only       ", dict(lm='forward')),
    ("LM, with backtracking       ", dict(lm='backtrack')),
    (f"clamp {R} + LM forward      ", dict(step_max=R, lm='forward')),
]

RADIUS = [("no clamp                    ", dict(step_max=0.0))] + [
    (f"clamp {r:<4g} (halving)        ", dict(step_max=r)) for r in
    (2.0, 1.5, 1.0, 0.7, 0.5, 0.4, 0.3, 0.2, 0.15)]

LM = [
    ("fixed lambda (baseline)     ", dict(step_max=0.0)),
    ("forward  up2 dn2            ", dict(lm='forward')),
    ("forward  up4 dn2            ", dict(lm='forward', lm_up=4.0)),
    ("forward  up2 dn1.4          ", dict(lm='forward', lm_down=1.4)),
    ("backtrack up2 dn2           ", dict(lm='backtrack')),
    ("backtrack up4 dn2           ", dict(lm='backtrack', lm_up=4.0)),
    ("forward + clamp             ", dict(lm='forward', step_max=R)),
    ("backtrack + clamp           ", dict(lm='backtrack', step_max=R)),
]

LAM = [(f"lambda {l:<6g}               ", dict(lam=l, step_max=0.0))
       for l in (0.02, 0.05, 0.1, 0.2, 0.4, 0.8)]


def main():
    args = sys.argv[1:]
    do = lambda k: (not args) or '--all' in args or k in args  # noqa: E731
    table = build_table()
    far = make_cases(300, 1.2, 0x1153)

    if do('--headline'):
        report("frozen 48-pose harness table", table, HEADLINE, latency=True)
        report("300 far-seeded poses (+-1.2 rad)", far, HEADLINE)
    if do('--radius'):
        report("trust radius - frozen table", table, RADIUS, latency=True)
        report("trust radius - far seeds", far, RADIUS)
    if do('--lm'):
        report("adaptive lambda - frozen table", table, LM, latency=True)
        report("adaptive lambda - far seeds", far, LM)
    if do('--lambda'):
        report("lambda itself - far seeds.  Watch max/med improve while every "
               "other column gets worse", far, LAM)


if __name__ == "__main__":
    main()
