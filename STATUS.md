# Status

Running log of what is built, what has been measured on hardware, and what is
open. `README.md` is the stable description of the project; this file is the
part that changes every build.

Last updated after adding the DLS trust region. The hardware numbers below are
still from the pre-trust-region bitstream (`run.log`); the new solver has
passed the host regression but has not been synthesised or run on the board.

---

## Where it stands

| | state |
|---|---|
| golden model (`model/validate.py`) | passing — 27,008 branch solutions, worst pose error 1.1e-13 |
| host regression (`make -C hls host`) | passing — all four kernels, float build |
| `ik_analytic_kernel` | **on hardware**, 80 MHz, measured |
| `ik_dls_kernel` | **on hardware**, 80 MHz, measured (`run.log`) |
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

### `ik_dls_kernel`, 80 MHz, 48 poses (`run.log`)

| path | min | med | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 48089 ns | 64036 ns | 255470 ns | 351175 ns | 97301 ns | **5.48** |
| PS (double) | 77160 ns | 103800 ns | 502458 ns | 741507 ns | 166006 ns | 7.14 |
| PL iterations | 3 | 4 | 16 | 22 | 6 | 7.33 |

Correctness passes on the first 8 poses (`status=0`, 3–4 iterations).

Two things to read off this and not miss:

- **PL is only 1.62× faster than the A9 here**, against 3.0× for the analytic
  kernel. An iterative solver spends its time in a loop the PL runs at 80 MHz
  and the A9 runs at 667 MHz; the PL wins on parallelism per iteration, and at
  this size that is a narrow win.
- **The PS spread (7.14) is worse than the PL's (5.48)** and the PS runs the
  same algorithm. Most of that is the same iteration-count spread; the rest is
  cache and double-precision softfloat. It is not a PL artefact.

---

## Diagnosis: where the 5.48× comes from

Divide each PL latency by the iteration count at the same rank:

| rank | latency | iters | ns per iteration |
|---|---|---|---|
| min | 48089 | 3 | 16030 |
| med | 64036 | 4 | 16009 |
| p95 | 255470 | 16 | 15967 |
| max | 351175 | 22 | 15963 |

**Per-iteration cost is constant to 0.4%.** A least-squares fit over those four
points gives

```
latency_ns  =  230  +  15952 * iterations
```

i.e. ~1276 PL cycles per iteration at 80 MHz, on a fixed overhead of **230 ns**.

So:

- The AXI4-Lite transaction overhead — the thing this harness was rebuilt to
  capture, and which is a real fraction of the analytic kernel's 10.7 µs — is
  **0.5% of a median DLS solve and 0.07% of the worst one.** There is nothing
  to win on the bus for this kernel.
- Every loop inside an iteration has a fixed trip count, and the measurement
  confirms it: the II choices in `ik_config.hpp`, the rolled CORDIC, the
  padded `spd::solve()` all cost the same cycles every time.
- **100% of the spread is iteration count**, and the iteration count spread
  (3 → 22) is itself partly by construction: `gen_vectors.py` selects a
  quarter of the pose table on the model's DLS iteration count, in the range
  [8, 32], seeded ±1.2 rad away. Without that block the table converged in 3–5
  every time and reported max/med = 1.77.

That last point matters for interpreting the number. **5.48× is a designed
worst-case characterisation, not what a tracking loop would see.** A control
loop seeds from the previous solution — a few milliradians away, not 1.2 rad —
and would sit at 3–4 iterations with almost no spread. The table is built to
find the tail on purpose, because the tail is what decides schedulability.

### Caveat on the table above

`run.log`'s latency and iteration summaries are two separately sorted arrays,
so pairing rank-for-rank is an inference, not a measurement. The inference is
strong (the residuals are 0.4%), but the harness should report the pairing
directly — see Open items.

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

**A1. Step clamping (trust region).** **Done** — see above. p95 18 → 11
iterations on the shipped table, deadline misses halved over ten tables.

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

**A5. Warm start, and measure it.** Add a harness mode that walks a trajectory
seeding each solve from the previous solution. This does not improve the
solver; it quantifies the gap between the worst-case characterisation above and
the operating case. Worth having as a second row in the report, because
"5.48× worst case, 1.0x× tracking" is a much more complete answer than either
number alone.

**A6. Cap `max_iter` to the deadline.** *Zero engineering cost — it is already a
runtime register.* At 15.95 µs per iteration:

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

Budget per iteration, ~1276 cycles. From the last saved HLS report
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

Bus-side work. 230 ns of 64,000. The AXI4-Lite commentary in `main.c` and
`docs/timing_methodology.md` is correct for the analytic kernel and does not
transfer to this one.

---

## Recommended order

1. **Build and run the trust region.** It has passed the host regression but
   has never been synthesised. `make -C hls csim` first — the clamp is a
   power-of-two shift and so is bit-identical between the two builds by
   construction, but that is an argument, not a measurement. Then `syn` (watch
   LUT: the clamp adds a compare chain and a variable shift, no DSP), then
   `gen_regmap.py` — **the kernel signature changed**, so the register map has
   moved and the checked-in provisional one is wrong on purpose.
2. **A6** — cap `max_iter` at the deadline, re-run, report the accuracy cost.
   No build needed, and it composes with the trust region: the clamp halves
   the number of poses that miss a 16-iteration deadline, and the cap bounds
   what the rest cost.
3. **B1** — `CORDIC_ITER` 18, validated through csim.
4. **B2** — retry `IK_BATCH_CORDIC=1` on top of B1.
5. **B3** — retry `--clk 100`.
6. **A5** — warm-start row in the harness, to report the operating case
   alongside the worst case.

Steps 1–2 change the distribution and need the pose table re-run to be
meaningful. Steps 3–5 are uniform scalings and can be verified from
`make -C hls reports` before they ever reach hardware.

---

## Open items

- **`report2.log` is stale.** It shows `grp_invert_fu` at 409 cycles / 96 DSP,
  i.e. a build from before `spd::solve()` replaced `mi::invert()` on the DLS
  path. The `spd` row in the budget table above is therefore an estimate. Run
  `make -C hls reports` and replace it before tuning against it.
- **`run.log` predates the current harness.** It has one `dls iterations` row;
  `sw/src/main.c` now prints `dls iters (PL)` and `dls iters (PS, double)`
  separately. The PS iteration counts — which are what would show whether
  quantisation is moving the convergence path — are not in this log. Re-run.
- **`ik_dls_kernel`'s register map has moved.** It gained a `step_max`
  argument, which shifts every scalar after it and relocates the three arrays.
  `sw/src/ik_regmap.h` carries an updated *provisional* guess; run
  `scripts/gen_regmap.py` after `make -C hls ip`. `main.c`'s write/read-back
  check will refuse to report timings if this is skipped.
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
