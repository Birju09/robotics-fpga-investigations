# Status

Current build and measurement state. `README.md` is the stable description
of the project and its findings; this file is what's built, what's measured,
and what's open.

## Built

| component | state |
|---|---|
| robot definition (`model/robot.py`) | PUMA 560, standard DH |
| URDF → DH converter (`model/urdf_to_dh.py`) | passing, round-trips to 1e-15 |
| golden model (`model/validate.py`) | passing, worst pose error 1.9e-13 |
| host regression (`make -C hls host`) | passing, all four kernels |
| `ik_analytic_kernel`, `ik_dls_kernel` | packaged; not co-resident (§6 of README — DSP budget) |
| `mat_mul_kernel`, `mat_inv_kernel` | packaged, not on the DLS critical path |

## Measured on hardware

| kernel | workload | status |
|---|---|---|
| `ik_dls_kernel` | 48-pose disturbed table | current — PUMA geometry, trust region on, repeated twice, agrees to 0.26% |
| `ik_dls_kernel` | pentagon trajectory | stale — measured on a since-retired pentagon geometry (self-intersected its own skeleton, violated a joint limit); corrected geometry not yet re-run |
| `ik_analytic_kernel` | 48-pose table | stale — predates the current manipulator and trust region |

## Open

- Re-run the trajectory workload on the corrected pentagon (`model/gen_vectors.py` already emits it; no synthesis needed).
- Re-synthesise `ik_analytic_kernel` on current geometry, so both kernels share a matched bitstream and workload.
- One DLS pose converges in 22 iterations where the reference model predicts 12, identically across runs — deterministic, unexplained.
- Per-iteration cost rose 7.7% across the manipulator retarget with no change to the DLS loop body — not yet reconciled against synthesis reports.
- The random 48-pose table is not yet checked against the self-collision proxy used to fix the trajectory (~5% of a uniform `qlim` sample fails it).
- Trust-region (step-clamp) before/after comparison is from the previous manipulator; no unclamped run exists on current geometry.
