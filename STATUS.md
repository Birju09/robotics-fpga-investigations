# Status

Running log of what is built, what has been measured on hardware, and what is
open. `README.md` is the stable description of the project; this file is the
part that changes every build.

Last updated after the **second hardware run on the PUMA 560 geometry**
(`run.log`, `ik_dls_kernel` alone, trust region enabled), which repeats
Workload A unchanged and re-samples the trajectory at 100 points per edge
instead of 10. The DLS numbers below are current. The analytic kernel has
**not** been re-synthesised since the retarget and its figures still come from
the previous 0.71 m arm; every such section says so.

The retarget and its motivation are in README §2.1 and §4.5. In short: the
previous arm had no lateral shoulder offset, so the wrist centre could sit
exactly on the joint-1 axis, where theta1 is undefined and `atan2(0, 0)`
returned zero instead of failing. PUMA's `d3 = 150.05 mm` removes that by
construction.

What the finer sampling changed, and it changed a lot: the tracking case is now
**2 iterations at all 500 samples** with a total latency spread of 0.51%, and
chain drift fell from 1922 to **183 urad** despite ten times as many solves per
lap. The two sampling rates together are more informative than either alone —
they show the tracking case's flatness is a function of step size, and that
drift is bounded by the convergence tolerance rather than accumulated. Both are
reported below.

Two things still want attention before anything is built on them: the worst
pose takes 22 iterations where the model says 12 (in *both* the PL and PS
columns, and identically in both runs, so it is neither quantisation nor
noise), and the per-iteration cost rose 7.7% across the retarget with no change
to the DLS inner loop. Both are under Open items.

---

## Where it stands

| | state |
|---|---|
| robot definition (`model/robot.py`) | PUMA 560, standard DH, provenance recorded per parameter |
| URDF -> DH converter (`model/urdf_to_dh.py`) | passing — 4 round-trip cases to 1e-15, 3 malformed inputs rejected |
| golden model (`model/validate.py`) | passing — 32,000 branch solutions, worst pose error 1.9e-13 |
| host regression (`make -C hls host`) | passing — all four kernels, float build, new geometry |
| `ik_dls_kernel` | **on hardware, PUMA geometry**, trust region on, both workloads, Workload A run twice (`run.log`) |
| pentagon trajectory workload | **on hardware, PUMA geometry**, at 100 and at 10 samples/edge (`run.log`) |
| `ik_analytic_kernel` | **stale** — previous arm, pre-trust-region bitstream, not re-synthesised |
| trust-region before/after comparison | **stale** — no unclamped run exists on this geometry |
| `mat_mul_kernel` / `mat_inv_kernel` | synthesise and package; not in either measured bitstream |
| both IK kernels in one bitstream | does not fit — ~345 DSP against 220 |

One kernel at a time is the operating mode, not a temporary state. Build twice
and compare at the same clock; `scripts/build_vivado.tcl --kernels` selects,
and `sw/src/main.c` compiles out whatever is absent.

---

## Hardware results

### `ik_analytic_kernel`, 80 MHz, 48 poses — **previous arm, stale**

| path | min | med | p95 | max | max/med |
|---|---|---|---|---|---|
| PL | 10716 ns | 10744 ns | 10781 ns | 10781 ns | **1.00** |
| PS (double) | 32446 ns | 32612 ns | 32732 ns | 37923 ns | 1.16 |

Constant to the resolution of the timer. That is the control: the spread below
is the solver, not the bus, not the timer, not the PS. That conclusion is about
the harness and survives the change of arm; the absolute numbers do not.

### `ik_dls_kernel`, trust region enabled, 48-pose table (`run.log`, PUMA)

| path | min | med | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 52144 ns | 69378 ns | 190372 ns | 380086 ns | 89600 ns | **5.47** |
| PS (double) | 78252 ns | 105356 ns | 321166 ns | 591052 ns | 137661 ns | 5.60 |
| PL iterations | 3 | 4 | 11 | 22 | 5 | 5.50 |
| PS iterations | 3 | 4 | 12 | 22 | 5 | 5.50 |

Correctness on the first eight poses: all `status=0`, worst joint delta 30–213
µrad, 3–4 iterations.

**This table has now been measured twice in the same bitstream, in separate
sessions**, which is the first repeatability figure the project has:

| | min | med | p95 | max | mean |
|---|---|---|---|---|---|
| run 1 | 52126 | 69560 | 190298 | 380153 | 89612 |
| run 2 | 52144 | 69378 | 190372 | 380086 | 89600 |
| delta | +0.03% | −0.26% | +0.04% | −0.02% | −0.01% |

Iteration quantiles are identical (3 / 4 / 11 / 22 PL, 3 / 4 / 12 / 22 PS), the
22-iteration outlier included, and the eight correctness lines are byte-for-byte
the same. So the kernel is deterministic over its inputs, the ~0.3% is the
instrument, and the 22-iteration pose is a real property of the kernel rather
than a marginal solve that could land either way. That last point matters: it
was previously written off as "one marginal solve"; it is not marginal, it is
repeatable and unexplained.

### `ik_dls_kernel`, pentagon trajectory, 500 samples at 100/edge (`run.log`)

| path | min | med | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 34913 ns | 35003 ns | 35052 ns | 35092 ns | 35003 ns | **1.00** |
| PS (double) | 51193 ns | 51307 ns | 51427 ns | 51720 ns | 51316 ns | 1.00 |
| PL iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |
| PS iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |

Vertex drift, samples 0/100/200/300/400: 91, 45, 45, 76, 91 µrad, all at 2
iterations. `ik_traj_iters_tbl` is 2 at all 500 entries, so kernel and model
agree at every sample.

### The two sampling rates, compared

The preceding run drove the same pentagon at 10 samples per edge. Same arm,
same path, same bitstream, same kernel — only the step size differs:

| | 10/edge (50 samples) | 100/edge (500 samples) |
|---|---|---|
| travel per sample | 23.5 mm | 2.35 mm |
| at a 1 kHz control period | 23.5 m/s | 2.35 m/s |
| iterations | 2 or 3, ~half each | **2 at every sample** |
| PL latency | 34944–52215 ns | 34913–35092 ns |
| `max/min` | 1.49 | **1.005** |
| peak chain drift | 1922 µrad | **183 µrad** |

Five observations, in decreasing order of how much they change the argument:

- **The tracking case is genuinely flat at servo-realistic step sizes, and the
  qualifier is the finding.** 500 solves, 2 iterations each, 179 ns of spread
  on 35 µs. Against 5.47 and 3–22 iterations on the 48-pose table, that is the
  project's central claim in its strongest measured form. But at 23.5 mm/step
  the same lap needs 2 or 3, so flatness is a property of the *sampling rate*,
  not of DLS. Note also which of the two rates is physical: 23.5 mm per 1 kHz
  period is 23.5 m/s, which no PUMA does. The coarse run is a stress case, not
  a slower-servo case.
- **Chain drift does not accumulate with the number of solves.** Ten times more
  solves per lap produced ten and a half times *less* drift: 1922 → 183 µrad
  (~126 → ~12 Q16.16 LSB, ~1.7 mm → ~160 µm at 0.86 m reach). What sets it is
  per-solve difficulty. Each sample stops as soon as `‖e‖ < tol`, so the two
  chains land somewhere inside the same tolerance ball rather than at the same
  point, and the separation grows with warm-start distance. Coarse stepping
  puts them near opposite sides of the ball; fine stepping keeps them near its
  centre. Drift is bounded by `tol`, not integrated along the path — which is
  what re-converging every period buys, now measured across a 10× sweep instead
  of argued from one lap.
- **This also mostly answers the previous run's open item.** 1922 µrad against
  the old arm's 152 µrad had three confounded causes (arm scale, path scale,
  path conditioning). Holding all three fixed and varying only the sampling
  interval reproduces most of the effect, so sampling interval dominates and
  the change of manipulator does not. Not a controlled cross-arm comparison,
  and none is available without rebuilding for the old geometry — but enough to
  stop treating the retarget as the suspect.
- **Quantisation did not move the convergence path anywhere.** PL and PS
  iteration counts agree at every quantile in both trajectory runs, and no
  `(model needed N)` line was printed.
- **PL is 1.47× faster than the A9** on the trajectory and 1.52× on the 48-pose
  table, against 3.0× for the analytic kernel on the old arm. An iterative
  solver spends its time in a loop the PL runs at 80 MHz and the A9 at 667 MHz;
  the PL wins on parallelism per iteration, and at this size that is a narrow
  win. The case for the PL path here is determinism, not throughput.

### Effect of the trust region — **previous arm, stale**

`run.log` confirms the clamp is live in the current build (`step_max = 1500
milli-rad`), but there is no unclamped run on PUMA geometry, so the before/after
below is the old arm's and is kept only for the reasoning:

| | med | p95 | max | max/med |
|---|---|---|---|---|
| unclamped | 64036 ns / 4 it | 255470 ns / 16 it | 351175 ns / 22 it | 5.48 |
| trust region 1.5 rad | 64495 ns / 4 it | **160692 ns / 10 it** | **513320 ns / 32 it** | 7.95 |

There, p95 latency fell 37% — which is where a scheduler lives — while the
maximum went 22 → 32 against a model that predicted 32 in both configurations.
That discrepancy was never resolved before the arm changed. Re-establishing the
clamp's effect on PUMA costs a re-run, not a re-synthesis: `IK_DLS_STEP_MAX` is
a runtime AXI register, so setting it beyond the joint range disables the clamp
in the existing bitstream.

---

## Diagnosis: where the spread comes from

Run 1 gave a two-point fit through its trajectory (2 iterations) and its
48-pose extreme (22 iterations), with no least-squares step:

```
latency_ns  =  423  +  17260 * iterations
```

i.e. ~1381 PL cycles per iteration at 80 MHz on a fixed overhead of **423 ns**.
**Run 2 was not refitted.** Every row below is that relation predicting a
measurement taken in a different session:

| row | iters | predicted | measured | error |
|---|---|---|---|---|
| trajectory, all 500 | 2 | 34943 ns | 34913–35092 ns | +0.17% at the median |
| 48-pose min | 3 | 52203 ns | 52144 ns | −0.11% |
| 48-pose med | 4 | 69463 ns | 69378 ns | −0.12% |
| 48-pose p95 | 11 | 190283 ns | 190372 ns | +0.05% |
| 48-pose max | 22 | 380143 ns | 380086 ns | −0.01% |

Nothing worse than 0.17%, across an 11× range of iteration counts and two
sessions. An ordinary least-squares refit to run 2's five points returns
`413 + 17260·k` — the slope recovered to better than 0.01%, the intercept to
2.4%. So 17260 ns/iteration is a reproducible constant of this bitstream, not
a curve-fitting artefact.

The harness's own `per iter` rows agree: `(423 + 22b)/22 = 17279` against 17273
measured as the 48-pose per-iteration minimum, and `(423 + 4b)/4 = 17366`
against 17372 at its median.

**The per-iteration rows look like a contradiction and are not.** The
trajectory's per-iteration cost (17501 ns) reads *higher* than the 48-pose
table's median (17372 ns), which would suggest its iterations are more
expensive. They are not: the difference is the 423 ns fixed overhead amortised
over two iterations instead of four. The 48-pose table's per-iteration
*minimum*, 17273 ns, belongs to its 22-iteration pose, where the overhead
nearly vanishes. Nothing inside an iteration is data-dependent, which is what
the per-iteration row was added to establish.

Consequences:

- The AXI4-Lite transaction overhead — a real fraction of the analytic kernel's
  10.7 µs — is **0.6% of a median DLS solve, 1.2% of a trajectory sample and
  0.1% of the worst case.** There is nothing to win on the bus for this kernel.
- Every loop inside an iteration has a fixed trip count, confirmed across an
  11× range of iteration counts and reproduced across two sessions: the II
  choices in `ik_config.hpp`, the rolled CORDIC and the padded `spd::solve()`
  cost the same cycles every time.
- **All of the spread is iteration count.** On the 48-pose table that spread
  (3 → 22) is partly by construction: `gen_vectors.py` selects a quarter of the
  table on the model's DLS iteration count in [8, 32], seeded ±1.2 rad away.
  Without that block the table converged in 3–5 every time and reported
  `max/med` = 1.77.

**The slope moved and the cause is not established.** The same fit on the
previous arm gave `361 + 16030·k`, so per-iteration cost is up 7.7% at
1381 cycles against 1282. The DLS inner loop was not touched by the retarget —
only the values of the DH constants changed, and those are compile-time
constants that should not alter the schedule. A clock difference would have to
be 74.3 MHz to explain it, which `build_vivado.tcl` does not offer. The likely
answer is simply a different synthesis run with a different schedule, but it
has not been checked against `make -C hls reports`. Until it is, absolute
latencies here are comparable among themselves and not against the analytic or
unclamped tables.

---

## Done: the DLS trust region

**All numbers in this section were derived on the previous 0.71 m arm and its
pose table.** The clamp itself is unchanged and enabled in the current build;
the *sizing* of it — the 1.5 rad radius and every figure justifying it — has
not been revisited on PUMA geometry, whose joints span a different range. The
reasoning is what is worth keeping here, and the design log in `ik_config.hpp`
records how it was arrived at.

`iks::dls()` bounds `|dq|∞` per iteration (`IK_DLS_STEP_MAX_DEFAULT`,
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
| unclamped | 48/48 | 4 | 18 | 32 | 6.40 | 287 µs |
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

**A1. Step clamping (trust region).** **Built, shipped, and measured — but only
on the previous arm.** There it took p95 16 → 10 iterations and p95 latency
255 → 161 µs, with deadline misses halved over ten model tables, while the
maximum went the other way (22 → 32 iterations) against a model that predicted
no change. The clamp is enabled in the current PUMA bitstream (`step_max = 1500
milli-rad` in `run.log`) but its *effect* on that geometry is unmeasured. A
re-run with `IK_DLS_STEP_MAX` written beyond the joint range gives the
unclamped arm of the comparison without a re-synthesis; see Recommended order.

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
answer is **52.1–380.1 µs and 3–22 iterations disturbed, against 34.91–35.09 µs
and 2 iterations throughout while tracking**. The chain is fed the kernel's own
previous output rather than a table lookup, so the drift figure above is a
second result from the same run. Running it at two sampling rates turned out to
be worth more than either rate alone — see the comparison table above.

**A6. Cap `max_iter` to the deadline.** *Zero engineering cost — it is already a
runtime register.* At the measured 17.26 µs per iteration:

| budget | `max_iter` | notes |
|---|---|---|
| 1.11 ms | 64 (current default) | exceeds a 1 kHz period outright |
| 250 µs | 14 | no pose in the model's table exceeds 12; the hardware's single 22-iteration pose does |
| 100 µs | 5 | 11 of 48 poses (23%) exceed 5 in the model; the tail returns best-effort |
| **35 µs** | **2** | never fires anywhere on the 100/edge trajectory; 48 of 48 poses truncated |

The last row is the one worth noticing, and the finer sampling is what makes it
interesting. At 100 samples per edge the tracking workload never exceeds two
iterations, so a `max_iter` of 2 caps the *whole system* at 35 µs — a hard,
fixed, sub-40 µs latency — at zero cost while tracking, degrading only when the
loop is disturbed. That is the project's argument reduced to a single register
write. At 10 samples per edge the same cap would fire on roughly half the lap,
which is precisely why the sampling rate has to be quoted with the claim.

This does not make the solver faster; it converts an unbounded latency into a
bounded latency plus a bounded accuracy loss, which is the trade a scheduler can
actually accept. DLS returns a usable partially converged `q` at any cutoff.
Combined with A1/A2 — which raise the fraction of poses that finish inside the
cap — this is the shape of the real answer.

### B. Cheaper iterations — attacks the absolute numbers

Budget per iteration, ~1381 measured cycles at 80 MHz. The block breakdown
below comes from the last saved HLS report (`report2.log`, **stale** — see Open
items) and sums to the *previous* 1282-cycle figure, so the shares are
indicative and the 99-cycle difference is unattributed:

| block | cycles | share |
|---|---|---|
| `fk_jacobian` | 590 | ~43% |
| `spd::solve` | ~250 (est.) | ~18% |
| `mm::multiply` ×2 | 83 | ~6% |
| `pose_error` | 14 | ~1% |
| staging / FSM / unattributed | remainder | ~32% |

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
products are redundant. Saves ~15 cycles of 1381. Real, but about 1% — do it
for tidiness, not for latency.

**B5. Account for the 7.7% slope increase** across the retarget (1282 → 1381
cycles per iteration). This is not a tuning item, it is a "the design changed
and nobody knows how" item, and it is worth more than B4 because it is the same
order as B1's expected saving. `make -C hls reports` on the current tree
against `report2.log` should settle it in one run and with no hardware.

### What is not worth doing

Bus-side work. 423 ns of 69,560. The AXI4-Lite commentary in `main.c` and
`docs/timing_methodology.md` is correct for the analytic kernel and does not
transfer to this one.

---

## Recommended order

~~1. Build and run the trust region.~~ **Done** — synthesised, run, numbers
above.

~~6. A5 — warm-start row in the harness.~~ **Done** — the pentagon trajectory,
numbers above.

~~4. Re-run DLS on the PUMA geometry.~~ **Done** — `run.log`, both workloads,
full console capture, Workload A repeated.

~~5. Decompose the chain drift.~~ **Mostly done** — the 10/edge and 100/edge
runs isolate sampling interval as the dominant term. See the comparison table.

What remains, in order. The first three cost re-runs of the *existing*
bitstream and no synthesis, so they should all be done in one session:

1. **A6 with `max_iter` = 2.** The finer sampling makes this the highest-value
   run in the list: the tracking workload never exceeds two iterations, so a
   cap of 2 pins the whole system at a hard 35 µs and fires only when the loop
   is disturbed. Run it, and report what the 48 disturbed poses cost in
   accuracy at that cap. Also worth capturing at 14 and 5 for the curve.
2. **Unclamped re-run on this geometry.** Write `IK_DLS_STEP_MAX` beyond the
   joint range and re-run both workloads. This is the missing arm of the
   trust-region comparison, which is currently the previous manipulator's.
3. **Instrument the one 22-iteration pose.** The model solves it in 12 and both
   the PL and PS columns say 22, identically in both runs, so this is
   kernel-vs-model and it is deterministic. Dump the per-iteration residual for
   that pose from `main.c` and from `model/sweep_dls.py` and compare the two
   sequences. Until this is understood, `max/med` = 5.47 should be quoted with
   the caveat that its numerator is a solve the reference does not reproduce.

Then, requiring builds:

4. **Re-synthesise and re-measure `ik_analytic_kernel`.** It now carries a
   different θ1 branch and a different wrist extraction, and its figures here
   predate the two-workload harness, so there is no analytic row for the
   trajectory at all. One build (`--kernels ik_analytic_kernel`), and it is what
   makes the fixed-latency-vs-data-dependent comparison a comparison over the
   same arm, the same two workloads, and comparable bitstreams.
5. **B5** — reconcile the 7.7% per-iteration slope increase against
   `make -C hls reports`. No hardware needed.
6. **Bracket the step size at which tracking stops being flat.** The sweep has
   two points an order of magnitude apart; 20, 30 and 50 samples per edge would
   locate the crossover and turn "flat while tracking" into a statement with a
   servo rate attached. Generator change and a re-run, no synthesis.
7. **B1** — `CORDIC_ITER` 18, validated through csim.
8. **B2** — retry `IK_BATCH_CORDIC=1` on top of B1.
9. **B3** — retry `--clk 100`.

Steps 1–3 change the distribution and need both tables re-run to be meaningful.
Steps 7–9 are uniform scalings and can be verified from `make -C hls reports`
before they ever reach hardware.

---

## Open items

- **`report2.log` is stale.** It shows `grp_invert_fu` at 409 cycles / 96 DSP,
  i.e. a build from before `spd::solve()` replaced `mi::invert()` on the DLS
  path. The `spd` row in the budget table above is therefore an estimate. Run
  `make -C hls reports` and replace it before tuning against it.
- **The worst pose takes 22 iterations where the model takes 12.** The model's
  table (`ik_model_iters_tbl`) has min 3, med 4, p95 11, max 12; the hardware
  reproduces the first three exactly and misses the last by a factor of two.
  Crucially, **both** the PL and the PS column read 22, and the PS column is
  this kernel's source compiled in double — so this is not quantisation, it is
  a difference between `iks::dls()` and `M.ik_dls()`. The known candidate is the
  convergence test, `err < tol` in the model against `err_sq < tol_sq` in the
  kernel: equivalent in exact arithmetic, not identical in rounding, and capable
  of compounding on a pose re-approaching the tolerance slowly. That mechanism
  already accounts for the 3-of-256 discrepancy below, but there it moves the
  count by one. Ten is a different claim and needs the residual sequences
  compared directly. Until then `max/med` = 5.47 has one solve behind its
  numerator, and that solve is one the reference cannot reproduce. **It is
  identical in both runs**, so it is not a marginal solve that could land either
  way — it is a repeatable property of the kernel that nothing currently
  explains, which makes it more concerning rather than less.
- **The per-iteration cost rose 7.7% across the retarget**, 1282 → 1381 cycles
  at 80 MHz, with no change to the DLS inner loop — only to the *values* of
  compile-time DH constants, which should not alter the schedule. A clock
  difference would have to be 74.3 MHz, which `build_vivado.tcl` does not offer.
  Most likely a different synthesis run scheduling differently, but that is a
  guess. `make -C hls reports` against `report2.log` settles it with no
  hardware. Until it does, absolute latencies from this run are not comparable
  with the pre-retarget ones.
- ~~**Chain drift is 12.6× worse than on the previous arm.**~~ **Largely
  answered.** It was 1922 µrad against the old arm's 152, with three variables
  confounded (arm scale, path scale, path conditioning). Holding all three fixed
  and taking the sampling interval from 23.5 mm to 2.35 mm per sample brings it
  to 183 µrad, so sampling interval is the dominant term and the retarget is
  not the suspect. What is *not* closed: this is still not a controlled
  cross-arm comparison, and the drift-versus-step-size relation is superlinear
  on two points (10× finer sampling, 10.5× less drift, with 10× more solves) —
  two points do not establish a law.
- **`ik_analytic_kernel` has not been re-synthesised since the retarget.** It
  carries a different θ1 branch and a different wrist extraction, so its
  figures need a rebuild rather than a re-run, and it has never been through
  the two-workload harness at all — the trajectory comparison currently has a
  DLS row and no analytic row.
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
- **The `per iter` rows came out flat and that is the result, not a formality.**
  Across 3–22 iterations on the 48-pose table the per-iteration figure spans
  17273–17427 ns, a 0.89% band, and the variation within it is fully explained
  by the 423 ns fixed overhead amortising. Nothing inside an iteration is
  data-dependent. The trust region was written to keep this true — the shift
  count is always computed and always applied, never conditionally skipped —
  and the measurement now confirms it directly rather than by inference across
  two separately sorted arrays.
- **The tracking case's flatness is conditional on step size, and the crossover
  is not located.** Two iterations at every one of 500 samples is a statement
  about a 2.35 mm step; at 23.5 mm the same lap takes 2 or 3. The sweep has two
  points an order of magnitude apart. 20, 30 and 50 samples per edge would
  bracket it, cost a generator change and a re-run, and would turn "DLS is
  fixed-latency while tracking" into a claim with a servo rate attached. Until
  then, always quote the sampling rate with the ratio.
- **One lap is not a stability proof.** Divergence stays under 100 µrad at every
  vertex, peaks at 183 µrad, and *falls* when the number of solves per lap goes
  up tenfold — strong evidence of a bounded, self-correcting chain consistent
  with each solve re-converging to `tol`. It is still not the same as showing
  the bound holds over minutes of operation: 500 samples is 0.5 s at 1 kHz, and
  a pentagon revisits the same five configurations repeatedly. Run several
  hundred laps, and a path that does not return to its own starting
  neighbourhood, before quoting a drift bound rather than a drift measurement.
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
