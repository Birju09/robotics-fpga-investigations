# Status

Running log of what is built, what has been measured on hardware, and what is
open. `README.md` is the stable description of the project; this file is the
part that changes every build.

Last updated after running the trust-region `ik_dls_kernel` and the pentagon
trajectory workload on hardware (`run.log`). Both were open at the previous
revision; both now have measured numbers, recorded below. The analytic figures
remain those of the earlier analytic-only bitstream.

---

## Where it stands

| | state |
|---|---|
| golden model (`model/validate.py`) | passing — 27,008 branch solutions, worst pose error 1.1e-13 |
| host regression (`make -C hls host`) | passing — all four kernels, float build |
| `ik_analytic_kernel` | **on hardware**, 80 MHz, measured (pre-trust-region bitstream) |
| `ik_dls_kernel` | **on hardware**, 80 MHz, measured with the trust region (`run.log`) |
| pentagon trajectory workload | **on hardware**, 80 MHz, measured (`run.log`) |
| `mat_mul_kernel` / `mat_inv_kernel` | synthesise and package; not in either measured bitstream |
| both IK kernels in one bitstream | does not fit — ~345 DSP against 220 |

One kernel at a time is the operating mode, not a temporary state. Build twice
and compare at the same clock; `scripts/build_vivado.tcl --kernels` selects,
and `sw/src/main.c` compiles out whatever is absent.

---

## Hardware results

### `ik_analytic_kernel`, 80 MHz, 48 poses

| path | min | med | p95 | max | max/med |
|---|---|---|---|---|---|
| PL | 10716 ns | 10744 ns | 10781 ns | 10781 ns | **1.00** |
| PS (double) | 32446 ns | 32612 ns | 32732 ns | 37923 ns | 1.16 |

Constant to the resolution of the timer. That is the control: the spread below
is the solver, not the bus, not the timer, not the PS.

### `ik_dls_kernel`, 80 MHz, trust region enabled, 48-pose table (`run.log`)

| path | min | med | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 48535 ns | 64495 ns | 160692 ns | 513320 ns | 89236 ns | **7.95** |
| PS (double) | 79064 ns | 106864 ns | 298227 ns | 871824 ns | 149556 ns | 8.15 |
| PL iterations | 3 | 4 | 10 | 32 | 5 | 8.00 |
| PS iterations | 3 | 4 | 11 | 32 | 5 | 8.00 |

### `ik_dls_kernel`, 80 MHz, pentagon trajectory, 50 samples (`run.log`)

| path | min | med | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 32390 ns | 32421 ns | 32464 ns | 32473 ns | 32418 ns | **1.00** |
| PS (double) | 51560 ns | 52150 ns | 52440 ns | 52566 ns | 52006 ns | 1.00 |
| PL iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |
| PS iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |

Four observations, in decreasing order of how much they change the argument:

- **The tracking case is flat.** Every one of the 50 samples converges in two
  iterations, and the latency spread across the lap is 83 ns on 32.4 µs — the
  resolution of the measurement, not a property of the solver. Set against 7.95
  on the 48-pose table, this is the quantitative form of the claim the project
  has been making qualitatively: DLS is a fixed-latency block while it is
  following a path, and its data-dependent latency is a property of the
  disturbed case, not of the algorithm as deployed.
- **The chained fixed-point solve is stable.** Cumulative drift of the hardware
  chain from the model's chain peaks at **152 µrad at sample 28** and falls back
  to 76 µrad by sample 40. Q16.16's LSB is 15.3 µrad, so the peak is ten LSB, it
  is not monotone, and it does not accumulate: each sample re-converges to
  `tol`, which re-anchors the chain to the commanded pose every control period.
  On a 0.5 m lever, 152 µrad is ~76 µm of tool position. This was the open
  question the trajectory table was added to answer.
- **Quantisation did not move the trajectory's convergence path.** PL and PS
  iteration counts are identical at every sample, and no `(model needed N)`
  line was printed.
- **PL is 1.61× faster than the A9** on the trajectory and 1.66× on the 48-pose
  table, against 3.0× for the analytic kernel. An iterative solver spends its
  time in a loop the PL runs at 80 MHz and the A9 at 667 MHz; the PL wins on
  parallelism per iteration, and at this size that is a narrow win. The case for
  the PL path here is determinism, not throughput.

### Effect of the trust region, measured

The same 48 poses, before and after, both at 80 MHz. The pose set is frozen by
construction (`_find_dls_tail()` selects on the unclamped solver on purpose), so
this is the same workload twice:

| | med | p95 | max | max/med |
|---|---|---|---|---|
| unclamped (previous `run.log`) | 64036 ns / 4 it | 255470 ns / 16 it | 351175 ns / 22 it | 5.48 |
| trust region 1.5 rad | 64495 ns / 4 it | **160692 ns / 10 it** | **513320 ns / 32 it** | 7.95 |

The prediction from the model was med 4, p95 11, max 32, `max/med` 8.00. The
median (unchanged, +0.7%) and the p95 (11 predicted, 10 measured — the kernel's
`err_sq < tol_sq` test crossing one iteration early, as documented below) match
closely. **The maximum does not.** The model expected 32 iterations on the worst
pose in both configurations; the hardware measured 22 unclamped and 32 clamped.

So `max/med` rising 5.48 → 7.95 is not the clamp failing to help — p95 latency
fell 37% and that is where a scheduler lives — but it is also not the "maximum
unchanged" the model predicted. The honest reading is that this table's extreme
is one pose, that pose is marginal (32 of a 64 cap), and the iteration count of
a marginal solve is not reproducible between double and Q16.16. `max/med` on
n=48 is therefore a statistic with one sample behind it. See Open items.

---

## Diagnosis: where the spread comes from

The trajectory workload settles the pairing question directly. All 50 of its
samples take exactly two iterations, so its latency needs no rank-matching, and
it supplies a clean two-iteration point to anchor against the 48-pose table's
32-iteration extreme. Fitting `a + b·iterations` through those two points alone:

```
latency_ns  =  361  +  16030 * iterations
```

i.e. ~1282 PL cycles per iteration at 80 MHz on a fixed overhead of **361 ns**.
Every other measured row follows from it without further fitting:

| row | iters | predicted | measured | error |
|---|---|---|---|---|
| 48-pose min | 3 | 48451 ns | 48535 ns | +0.17% |
| 48-pose med | 4 | 64481 ns | 64495 ns | +0.02% |
| 48-pose p95 | 10 | 160661 ns | 160692 ns | +0.02% |
| 48-pose max | 32 | 513320 ns | 513320 ns | 0.00% |
| trajectory (all 50) | 2 | 32421 ns | 32390–32473 ns | ±0.16% |

The harness's own `per iter` rows agree to within 3 ns: `(361 + 2b)/2 = 16210`
against 16209 measured on the trajectory, `(361 + 4b)/4 = 16120` against 16123
at the 48-pose median, `(361 + 32b)/32 = 16041` against 16040 at its minimum.

**That last set is worth reading carefully, because it looks like a
contradiction and is not.** The trajectory's per-iteration cost (16209 ns) reads
*higher* than the 48-pose table's (16123 ns), which would suggest its iterations
are more expensive. They are not: the difference is entirely the 361 ns fixed
overhead amortised over two iterations instead of four. The 48-pose table's own
per-iteration *minimum*, 16040 ns, belongs to its 32-iteration pose, where the
overhead nearly vanishes. Nothing inside an iteration is data-dependent, which
is precisely what the per-iteration row was added to establish.

Consequences:

- The AXI4-Lite transaction overhead — a real fraction of the analytic kernel's
  10.7 µs — is **0.6% of a median DLS solve, 1.1% of a trajectory sample and
  0.07% of the worst case.** There is nothing to win on the bus for this kernel.
- Every loop inside an iteration has a fixed trip count, and the measurement
  confirms it across a 16× range of iteration counts: the II choices in
  `ik_config.hpp`, the rolled CORDIC and the padded `spd::solve()` cost the same
  cycles every time.
- **All of the spread is iteration count.** On the 48-pose table that spread
  (3 → 32) is partly by construction: `gen_vectors.py` selects a quarter of the
  table on the model's DLS iteration count in [8, 32], seeded ±1.2 rad away.
  Without that block the table converged in 3–5 every time and reported
  `max/med` = 1.77.

The two workloads therefore bracket the answer rather than competing:
**7.95 is the disturbed case and 1.00 is the tracking case**, both measured on
the same kernel, in the same bitstream, in the same run. The 48-pose table is
built to find the tail on purpose, because the tail decides schedulability; the
trajectory is what the machine does the rest of the time.

---

## Done: the DLS trust region

`iks::dls()` now bounds `|dq|∞` per iteration (`IK_DLS_STEP_MAX_DEFAULT`,
1.5 rad, a runtime AXI register like `lambda` and `tol`). The full design log —
radius sweep, the shift-vs-exact-scale decision, and the adaptive-λ result —
is in the trust region block of `hls/include/ik_config.hpp`;
`model/sweep_dls.py` reproduces every number.

**Why.** At λ = 0.02 the damping is nearly nil, so the step is essentially full
Newton. Instrumenting the model over far seeds measured `|dq|∞` reaching
**12.7 rad** on an arm whose joints span ±2.9 — the solver leaves the region
where the linearisation means anything and spends its budget finding its way
back.

**Predicted effect on the shipped 48-pose table** (from the model; the kernel
is bit-comparable but unbuilt):

| | conv | med | p95 | max | mean | p95 latency |
|---|---|---|---|---|---|---|
| unclamped (what `run.log` ran) | 48/48 | 4 | 18 | 32 | 6.40 | 287 µs |
| trust region 1.5 rad | 48/48 | 4 | **11** | 32 | 5.46 | **176 µs** |

Over ten independently drawn tables (480 poses), p95 goes 18 → 12, mean
6.38 → 5.82, and the count missing a 16-iteration deadline **halves, 30 → 14**,
at a cost of ~3 poses in 480 that stop converging inside the 64 cap.

**What it does not do.** The median does not move — the clamp never fires on a
well-seeded solve — and on this particular table the single worst pose is
unchanged at 32 iterations, so `max/med` stays at 8.00. The gain is in the body
of the distribution and in the deadline-miss rate, which is what a scheduler
actually cares about; it is not a smaller extreme.

### Two things this experiment got wrong first, both worth keeping

**The 48-pose table is too small to choose a parameter on.** A radius of 1.0
looked optimal on one draw and flipped two poses from converging to not on
another. The radius curve is flat between 1.5 and 2.0 and the differences below
that are inside the noise at n=48. Anything tuned against this table needs to
be re-checked over several draws — `sweep_dls.build_table()` takes a seed for
exactly this.

**The tail block's selection is frozen, but frozen is not neutral.**
`gen_vectors._find_dls_tail()` keeps poses the *unclamped* solver converges on
in 8–32 iterations. Freezing it is necessary — if the criterion tracked the
solver, an improved solver would just be handed harder poses and would report
no improvement — but it means the baseline scores 480/480 with max = 32 *by
construction*, and any successor can only lose poses from a set defined by its
predecessor succeeding on them. On an unselected far-seeded population the
clamp raises convergence 277/300 → 288/300 and cuts the mean from 14.82 to
11.10 iterations. The table understates it. Both numbers are in
`ik_config.hpp`; quoting only the table one would be quoting the biased one.

A genuinely neutral tail block would select on seed distance alone rather than
on the baseline's iteration count. That changes which poses are in the table,
so it should be done deliberately and with the old table kept for comparison —
listed under Open items.

### Adaptive λ: measured, and not implemented

Both forms lose on this architecture, so the kernel keeps a fixed λ:

| | conv | med | p95 | max | mean |
|---|---|---|---|---|---|
| baseline | 48/48 | 4 | 17 | 32 | 6.58 |
| trust region | 48/48 | 4 | 10 | 17 | 5.31 |
| LM, raise λ only | 47/48 | 4 | 18 | 64 | 7.54 |
| LM, with backtracking | 43/48 | 4 | 64 | 64 | 11.17 |

Textbook Levenberg–Marquardt rejects a step that increased the residual and
retries with more damping. That trade is good when re-evaluating the residual
is cheap next to a full iteration. Here it is not — this kernel's cost *is* one
FK+Jacobian per loop pass, so a rejected step burns a whole 16 µs iteration.
Over 300 far-seeded poses backtracking produced 4,995 rejections, and λ
ratchets up faster than the accept path brings it down, so the solver ends up
crawling under heavy damping: convergence fell to 43/48. The non-backtracking
form avoids the wasted pass, is roughly neutral on its own, and adds nothing
once the clamp is in.

Raising λ outright is worse still, and produces the clearest warning in this
project about the headline metric:

| λ | conv (of 300) | med | mean | max/med |
|---|---|---|---|---|
| 0.02 | 277 | 9 | 14.82 | 7.11 |
| 0.1 | 239 | 19 | 29.11 | 3.37 |
| 0.4 | 59 | 64 | 59.53 | **1.00** |

`max/med` improves monotonically as the solver gets worse, and reaches a
perfect 1.00 when almost nothing converges at all. **It is a figure of merit,
not an objective.** Quote it with the median and the convergence rate.

The diagnostic also found that ~40% of poses that never converge are sitting at
a genuine local minimum of ‖e‖², not overshooting. No damping policy fixes
that; it needs a restart from a different seed.

---

## Reducing the spread

Two independent axes, and they are worth keeping separate because only one of
them changes `max/med`:

- **Fewer iterations on the hard poses** → shrinks `max/med`. This is the
  real-time problem.
- **Cheaper iterations** → scales every latency by the same factor, `max/med`
  unchanged, but the absolute worst case falls.

### A. Fewer iterations — attacks the ratio

**A1. Step clamping (trust region).** **Done and measured on hardware** — see
above. p95 16 → 10 iterations and p95 latency 255 → 161 µs on the shipped table;
deadline misses halved over ten model tables. The maximum went the other way
(22 → 32 iterations), which the model did not predict and which is discussed
above.

**A2. Adaptive λ (Levenberg–Marquardt).** **Measured, rejected** — see above. A
rejected step costs a full iteration in this architecture, which is fatal to
the usual LM trade.

**A3. Balance the pose-error units.** `e` mixes position in metres (O(0.1))
with an orientation term that is O(1) and dimensionless. The least-squares
weighting between them is therefore accidental. Scaling the rotation block by a
characteristic length (~0.28 m, the forearm) makes the two comparable and
usually improves conditioning. Six multiplies per iteration. Note this
redefines what `tol` means, so the convergence threshold must be re-swept.

**A4. Seed DLS from the analytic solver.** Would essentially eliminate the tail
— the analytic result is exact where it exists, so DLS converges in one
iteration. Blocked two ways: the two kernels do not fit on this part
(~345 DSP against 220), and it would dissolve the comparison the project
exists to make. Note it as the hybrid-architecture answer and leave it.

**A5. Warm start, and measure it.** **Done and measured on hardware** — the
pentagon trajectory table (`ik_vectors.h` version 4) and the trajectory section
of `sw/src/main.c`. It does not improve the solver; it quantifies the gap
between the worst-case characterisation above and the operating case, and the
answer is **7.95× disturbed against 1.00× tracking, at 2 iterations flat**. The
chain is fed the kernel's own previous output rather than a table lookup, so the
drift figure above is a second result from the same run.

**A6. Cap `max_iter` to the deadline.** *Zero engineering cost — it is already a
runtime register.* At the measured 16.03 µs per iteration:

| budget | `max_iter` | notes |
|---|---|---|
| 1.02 ms | 64 (current default) | exceeds a 1 kHz period outright |
| 250 µs | 15 | ~5–10% of this pose table returns `IK_ERR_NO_CONV` |
| 100 µs | 6 | bounded hard; the tail poses return best-effort |

This does not make the solver faster; it converts an unbounded latency into a
bounded latency plus a bounded accuracy loss, which is the trade a scheduler can
actually accept. DLS returns a usable partially converged `q` at any cutoff.
Combined with A1/A2 — which raise the fraction of poses that finish inside the
cap — this is the shape of the real answer.

### B. Cheaper iterations — attacks the absolute numbers

Budget per iteration, ~1282 measured cycles. From the last saved HLS report
(`report2.log`, **stale** — see Open items):

| block | cycles | share |
|---|---|---|
| `fk_jacobian` | 590 | ~46% |
| `spd::solve` | ~250 (est.) | ~20% |
| `mm::multiply` ×2 | 83 | ~7% |
| `pose_error` | 14 | ~1% |
| staging / FSM | remainder | ~26% |

**B1. `CORDIC_ITER` 24 → 18.** `fk_jacobian` is the biggest block and six
serial rolled CORDIC calls are most of it. CORDIC resolves about one bit of
angle per iteration; 18 iterations gives 2⁻¹⁸ ≈ 3.8e-6 rad, still four times
finer than the Q16.16 output LSB of 1.5e-5. 24 is buying six bits that are
thrown away on the store. Estimated saving ~8% of the iteration, plus a
proportional LUT saving. `CORDIC_ITER` is shared by the float and fixed builds,
so the host regression stays meaningful; validate through `make -C hls csim`
against the golden vectors, since this moves the fixed-point trajectory.

**B2. `IK_BATCH_CORDIC`, reconsidered.** Currently off because the kernel came
to 53,493 LUT against 53,200 — over by 293. That measurement was taken *after*
`spd::solve()` replaced `mi::invert()`, so it is current, not stale, and on
this part the switch stays off. But it costs ~12.5k LUT for ~11% of the
iteration, and B1 would shrink the unrolled datapath by a quarter. 41k + ~9.4k
is inside the budget on paper. Worth one synthesis run **after** B1 lands: if
it fits, B1 and B2 together are ~20% off every iteration.

**B3. Retry 100 MHz.** 80 MHz is the default because a 100 MHz build missed
setup by WNS = −1.582 ns — but that was the analytic kernel at 75% DSP
occupancy. `ik_dls_kernel` places at ~180/220 DSP, less congested, and routing
congestion is what ate the margin. A `--clk 100` run is one build and would be
25% off every latency if it closes. Quote any result with its clock.

**B4. Exploit the symmetry of JJᵀ.** `A = JJᵀ` is symmetric; 21 of the 36 dot
products are redundant. Saves ~15 cycles of 1276. Real, but about 1% — do it
for tidiness, not for latency.

### What is not worth doing

Bus-side work. 361 ns of 64,495. The AXI4-Lite commentary in `main.c` and
`docs/timing_methodology.md` is correct for the analytic kernel and does not
transfer to this one.

---

## Recommended order

~~1. Build and run the trust region.~~ **Done** — synthesised, run, numbers
above.

~~6. A5 — warm-start row in the harness.~~ **Done** — the pentagon trajectory,
numbers above.

What remains, in order:

1. **Re-measure `ik_analytic_kernel`.** Its figures in this file are from the
   pre-trust-region bitstream and predate the two-workload harness, so there is
   no analytic row for the trajectory at all. It is one build (`--kernels
   ik_analytic_kernel`) and it is what makes the fixed-latency-vs-data-dependent
   comparison a comparison over the same two workloads rather than one.
2. **Resolve the 22-vs-32 iteration discrepancy** on the worst pose (Open
   items). Until it is understood, `max/med` on this table should be quoted with
   the caveat that its numerator is a single marginal solve.
3. **A6** — cap `max_iter` at the deadline, re-run, report the accuracy cost.
   No build needed, and it composes with the trust region: the clamp halves the
   number of poses that miss a 16-iteration deadline, and the cap bounds what
   the rest cost. Now also worth quoting against the trajectory, where a cap of
   3 would never fire.
4. **B1** — `CORDIC_ITER` 18, validated through csim.
5. **B2** — retry `IK_BATCH_CORDIC=1` on top of B1.
6. **B3** — retry `--clk 100`.

Steps 2–3 change the distribution and need both tables re-run to be meaningful.
Steps 4–6 are uniform scalings and can be verified from `make -C hls reports`
before they ever reach hardware.

---

## Open items

- **`report2.log` is stale.** It shows `grp_invert_fu` at 409 cycles / 96 DSP,
  i.e. a build from before `spd::solve()` replaced `mi::invert()` on the DLS
  path. The `spd` row in the budget table above is therefore an estimate. Run
  `make -C hls reports` and replace it before tuning against it.
- **`run.log` is a partial capture.** It begins mid-way through the trajectory
  vertex lines, so the register-map check, the correctness pass and trajectory
  vertices 0–3 are not in it. The register map was correct — `main.c` returns
  before printing any timing row if `verify_regmap()` fails, so the latency
  sections existing is proof — but the per-pose correctness output for this
  bitstream was not saved. Capture the whole console next run.
- **The worst pose took 22 iterations unclamped and 32 clamped**, where the
  model predicted 32 in both cases. A trust region cannot make a solve take
  *more* iterations in general — it can, on a pose where the unclamped full
  Newton step happened to jump straight into the basin — but this is one pose at
  the edge of the 64 cap, and it alone sets `max/med` for the whole table. Worth
  reproducing in `model/sweep_dls.py` at that pose specifically, comparing the
  clamped and unclamped step sequences, before any weight is put on 7.95 as a
  worst case. It is also the strongest argument yet for the n=48 caution below.
- **`ik_analytic_kernel` has not been run under the two-workload harness.** Its
  numbers here are from the earlier bitstream and cover the 48-pose table only,
  so the trajectory comparison currently has a DLS row and no analytic row.
- **The DLS tail block's selection is biased toward the unclamped solver** (see
  above). Worth replacing with a criterion that does not reference any
  solver's iteration count — seed distance alone — and keeping the current
  table alongside it for continuity. Until then, quote the unselected
  far-seed numbers next to the table numbers.
- **3 of 256 DLS vectors converge one iteration earlier in the kernel than in
  the model**, always by exactly one and always earlier. Verified identical at
  `HEAD`, so it predates the trust region; it is the convergence test being
  `err < tol` in the model and `err_sq < tol_sq` in the kernel, which are
  equivalent in exact arithmetic and differ by a rounding step right at the
  threshold. Benign, and `main.c` already reports it per pose ("model needed
  N"). It is also the check that confirmed the clamp is bit-identical across
  the two builds — the mismatch set did not change when it was added.
- **`validate.py`'s DLS section runs at λ=0.08, tol=1e-5** — settings the
  kernel never uses — and reports 422/600 and 200/300 convergence, which reads
  as a failure and is not one. Confirmed identical at `HEAD`, so it predates
  the trust region. Either run it at the kernel's λ=0.02 / tol=1e-3 or say in
  the output why those settings are deliberately harsher.
- **The harness now reports latency ÷ iterations** (`dls per iter (PL)`), so
  the "all spread is iteration count" claim becomes a measurement rather than
  a fit across two separately sorted arrays. It should come out flat at
  ~15,950 ns; if it does not, something in the iteration has become
  data-dependent. The trust region was written to keep it flat — the shift
  count is always computed and always applied, never conditionally skipped.
- **One lap is not a stability proof.** Chain drift peaked at 152 µrad and came
  back down within 50 samples, which is evidence of a bounded, self-correcting
  chain and is consistent with each solve re-converging to `tol`. It is not the
  same as showing the bound holds over minutes of operation, and a pentagon
  revisits the same five configurations repeatedly. Run several hundred laps,
  and a path that does not return to its own starting neighbourhood, before
  quoting a drift bound rather than a drift measurement.
- The matrix IPs have never been measured on hardware. They are packaged and
  the harness section for them exists; they have simply never been in a
  bitstream with a kernel worth measuring alongside.

---

## Build history

Kept because each entry is a defect that cost real time to find.

- **Over budget by 3×.** The first full `syn`/`ip`/Vivado run failed
  `place_design` needing roughly 3× the part's CARRY4/DSP/LUT. Cause was the
  schedule, not the algorithm: `INLINE` + `UNROLL` on 24–32 iteration
  CORDIC/shift-subtract loops gave every call site its own spatial copy of the
  adder chain. Rolled them, switched to `INLINE off` so call sites share one
  engine. Same pattern in `dh_step()`'s rotation multiply and
  `ik_analytic.cpp`'s wrist rotation.
- **`ik_acc_t` was `AP_RND, AP_SAT`** even though its own header says rounding
  happens once, on the final narrowing cast. Every `acc +=` was paying for
  round/saturate hardware that could never fire. Now `AP_TRN, AP_WRAP`.
  `ik_real_t` is untouched — its `AP_RND`/`AP_SAT` are load-bearing.
- **Auto-pipelining flattened the rolled loops back out.** Vitis HLS pipelines
  small loops with no explicit directive, which fully unrolled `dh_step()`'s
  matrix multiply inside the `fk()`/`rot03()`/`fk_jacobian()` chains. Explicit
  `#pragma HLS PIPELINE off`.
- **Three `mm::multiply()` call sites synthesised as three instances.**
  `#pragma HLS ALLOCATION` did not enforce sharing in this release (tried,
  verified no effect). Routed them through a shared loop with value-muxed
  operand staging — arrays of pointers are not synthesisable in Vitis HLS.
- **`mi::invert()` replaced by `spd::solve()`** on the DLS path. `A = JJᵀ+λ²I`
  is SPD by construction, so an LDLᵀ solve is roughly a third of the arithmetic
  of a full Gauss-Jordan inverse, and it absorbs the third product the old form
  needed. `matinv.cpp` is still built and packaged as `mat_inv_kernel`.
- **Missed timing at 100 MHz, WNS = −1.582 ns** (hold fine, WHS +0.033). HLS
  schedules against a 12.5% clock uncertainty estimated before routing exists,
  and a design at 75% DSP occupancy does not get the routes it assumed. Default
  PL clock is now 80 MHz, ~0.9 ns of margin. `build_vivado.tcl` now refuses to
  write an XSA when timing fails (`--allow-timing-fail` overrides) and prints
  the worst paths first — a bitstream with negative slack still programs and
  still returns answers, occasionally wrong ones, which is indistinguishable
  from a kernel bug.
- **`report_utilization -hierarchical` after synthesis, always.** The `UTLZ-1`
  DRC reports a whole-device total, which says you are over but not which
  kernel spent it. The per-kernel table is what you actually need and is free.
- **Every latency read back as exactly zero.** `libxiltimer` arms its default
  timer from an `__attribute__((constructor))` in `xiltimer.c`, meant to run
  before `main()`. On this BSP it did not finish in time. `ik_timer_init()`
  now calls `XilSleepTimer_Init()` explicitly and does a `usleep(1)`, which is
  what actually starts the counter. Both PL and PS read zero, which is what
  pointed at the timer rather than at either kernel.
- **The DLS pose table was not sampling the tail.** Seeded ±0.25 rad and
  stressed only geometrically, every pose converged in 3–5 iterations and
  hardware reported max/med = 1.77 against 25 iterations in the 256-pose host
  sweep. Geometric stress degrades the *analytic* solver's branch selection; it
  is not what makes DLS iterate. `_find_dls_tail()` now selects a quarter of
  the table on the model's own iteration count, spread across the tail rather
  than clustered at the worst case.
- **DLS was being checked against the analytic branch.** IK is multi-valued;
  DLS converges to whichever branch its seed is nearest and lands ~π from
  `ik_qgold_tbl` in a wrist joint while solving the pose exactly. The column
  read ~3.1 rad of "error" on poses solved to 1e-4 m. Compare against
  `ik_qdls_tbl`.
