# Latency Determinism of Closed-Form and Iterative Inverse Kinematics on a Zynq-7000 SoC

An experimental comparison of two inverse-kinematics (IK) solvers for a 6-DOF
anthropomorphic manipulator, implemented as Vitis HLS IP for the programmable
logic of a Xilinx Zynq-7000 (xc7z020clg400-1, PYNQ-Z1/Z2) and invoked from the
Cortex-A9 processing system.

## Abstract

Closed-form and iterative IK solvers are ordinarily compared on accuracy,
generality, and mean execution time. This work holds accuracy and arithmetic
fixed and compares them on a property that matters for their use inside a
real-time control loop: whether execution time is a constant of the
implementation or a function of the input. Both solvers are synthesised from
the same source tree, in the same fixed-point format (Q16.16), and measured on
the same device with the same instrument, so that any difference in latency
distribution is attributable to solver structure rather than to arithmetic,
target, or measurement method.

Two workloads are used. The first comprises 48 independent target poses,
deliberately including geometrically ill-conditioned configurations and a tail
selected on the reference model's iteration count; it characterises the
disturbed case. The second is a 50-sample pentagon trajectory in which the
damped-least-squares (DLS) solver is warm-started from its own previous output,
as a servo loop would drive it; it characterises the tracking case.

On hardware at 80 MHz, the DLS kernel exhibits a ratio of maximum to median
latency of 7.95 over the 48-pose set and 1.00 over the trajectory, with
identical arithmetic and configuration in both cases. Total latency is
described across the full measured range by `361 ns + 16030 ns × iterations`,
with residuals below 0.17%, indicating that the observed variance is entirely
attributable to iteration count rather than to any data dependence within an
iteration. The chained fixed-point solve is shown to be numerically stable over
one lap, with cumulative drift from the double-precision reference peaking at
152 µrad (approximately ten Q16.16 LSB) and subsequently decreasing.

These results are single-run measurements on one device and are qualified
accordingly in [Threats to validity](#7-threats-to-validity).

---

## 1. Motivation

An IK solver placed in a servo loop must satisfy a deadline on every
invocation, not on average. The two solver families differ structurally in this
respect:

| | analytic | DLS |
|---|---|---|
| method | closed form, wrist decoupling | `q ← q + Jᵀ(JJᵀ + λ²I)⁻¹e` |
| loop structure | all trip counts fixed at compile time | data-dependent exit condition |
| latency | a single value, which is also the worst case | a distribution |
| failure mode | returns `IK_ERR_UNREACH` | exhausts its iteration budget |
| generality | this kinematic structure only | any differentiable kinematics |

The closed-form solver is not universally preferable: it exists only for
manipulators with a solvable structure (here, a spherical wrist permitting
position/orientation decoupling), whereas DLS applies to arbitrary kinematics
and degrades gracefully near singularities. The question this work addresses is
narrower and quantitative — what does the iterative solver's latency
distribution actually look like on a specific device under a specific
workload — because that distribution, and not the mean, determines
schedulability.

The reference model, run in double precision over 256 randomly drawn reachable
poses, converges in a median of 4 iterations with a 95th percentile of 7 and a
maximum of 25. The ratio between typical and worst case observed in software
motivated the hardware study reported here.

---

## 2. System under test

### 2.1 Manipulator and conventions

A 6R anthropomorphic arm with a spherical wrist, described by standard (distal)
Denavit–Hartenberg parameters:

| i | aᵢ (m) | αᵢ | dᵢ (m) |
|---|---|---|---|
| 1 | 0.05 | +π/2 | 0.15 |
| 2 | 0.30 | 0 | 0 |
| 3 | 0.02 | +π/2 | 0 |
| 4 | 0 | −π/2 | 0.28 |
| 5 | 0 | +π/2 | 0 |
| 6 | 0 | 0 | 0.08 |

Pose is represented as `{x, y, z, roll, pitch, yaw}` with
`R = Rz(yaw)·Ry(pitch)·Rx(roll)`. All quantities crossing the AXI4-Lite
interface are signed Q16.16 (`value = raw / 65536`).

`model/ik_model.py` is the normative statement of these conventions; every
other artefact in the repository is validated against it.

### 2.2 Synthesised IP

| IP | top function | characteristic |
|---|---|---|
| `mat_mul_kernel` | `hls/src/matmul.cpp` | runtime dimensions ≤6×6, optional transpose on either operand |
| `mat_inv_kernel` | `hls/src/matinv.cpp` | Gauss–Jordan with partial pivoting |
| `ik_analytic_kernel` | `hls/src/ik_analytic.cpp` | compile-time-constant latency |
| `ik_dls_kernel` | `hls/src/ik_dls.cpp` | reports iterations consumed |

`ik_dls_kernel` invokes `mm::multiply()` and `spd::solve()` as C++ functions
rather than issuing AXI transactions to the standalone matrix IPs; the same
sources serve both export targets. Routing each 6×6 product over the
interconnect would add hundreds of cycles per iteration to a loop executing 2–64
times per solve, and would render the experiment a measurement of bus
arbitration. The standalone IPs exist so that the matrix primitives can be
characterised independently.

`spd::solve()` performs an LDLᵀ solve and replaced `mi::invert()` on the DLS
path because `A = JJᵀ + λ²I` is symmetric positive definite by construction, for
which a solve is cheaper than a full inverse. `matinv.cpp` remains built and
packaged but is no longer on the DLS critical path.

### 2.3 Repository layout

```
model/          floating-point reference model and vector generation (Python)
hls/include/    shared headers: numeric types, geometry, CORDIC/sqrt, configuration
hls/src/        kernel sources
hls/tb/         C++ testbenches, executed under csim and as a host regression
hls/cfg/        Vitis HLS 2025.2 configuration files, one per packaged IP
scripts/        Vivado block design, register-map extraction, Vitis application build
sw/src/         bare-metal driver and timing harness
docs/           technical manual and timing methodology
```

---

## 3. Numeric format

The fixed-point format was selected by measurement rather than by convention.
`model/validate.py` executes the complete quantised DLS loop in Q16.16, Q14.18,
Q12.20, Q14.26 and Q16.32 over 300 poses:

```
  fmt       lam     tol     conv  med  p95
  Q16.16   0.02   1e-03  300/300    4    6
  Q14.26   0.02   1e-03  300/300    4    6
  Q16.32   0.02   1e-03  299/300    4    6
```

Convergence in this regime is bounded by the damping factor rather than by word
length, so wider formats yield no measurable benefit and the narrowest adequate
format was adopted. Products accumulate in an exact Q32.32 accumulator
(`ik_acc_t`) and are rounded once on the narrowing store.

Three fixed-point types are used, with distinct rounding and overflow
behaviour justified in `hls/include/ik_types.hpp`:

| type | format | rounding / overflow | rationale |
|---|---|---|---|
| `ik_real_t` | `ap_fixed<32,16>` | `AP_RND, AP_SAT` | crosses IP boundaries and is bit-compared against the model; saturation prevents a joint angle sign-wrapping on overflow |
| `ik_acc_t` | `ap_fixed<64,32>` | `AP_TRN, AP_WRAP` | Q32.32 holds the full six-term dot product of two Q16.16 operands exactly, so rounding and saturation hardware could never fire |
| `ik_work_t` | `ap_fixed<32,16>` | `AP_TRN, AP_SAT` | internal storage; truncates to save an adder per store but retains saturation, since factorisation intermediates in `spd.cpp` can grow |

One consequence merits note. A tolerance of 1e-3 squares to 1e-6, two decades
below the Q16.16 LSB. The DLS convergence test is therefore evaluated in the
accumulator type; evaluated in Q16.16 it would quantise to zero and the solver
would report success on its first iteration.

Compiling with `-DIK_USE_FLOAT` builds the identical kernel algorithms in
double precision with no `<ap_fixed.h>` dependency. This is the mechanism by
which quantisation error is attributed separately from algorithmic error.

---

## 4. Method

### 4.1 Instrument

`sw/src/main.c` reports PS wall-clock time, measured with the Cortex-A9 global
timer, around the complete transaction: AXI4-Lite argument writes, `ap_start`,
the polling loop, and result reads. This is the interval a control loop
experiences.

This quantity is deliberately distinct from the latency Vitis HLS reports,
which counts PL cycles between `ap_start` and `ap_done` and excludes
approximately twenty single-beat AXI4-Lite accesses per solve. For kernels of
this size the difference is not negligible, and the relationship between the two
is itself a subject of the study. See
[docs/timing_methodology.md](docs/timing_methodology.md).

The harness reports a five-number summary (min, median, p95, max, mean) per
path, the ratio of maximum to median as the jitter figure of merit, the
iteration count distribution for DLS on both PL and PS paths, and latency
divided by iteration count. The last of these exists so that the claim that
variance is attributable to iteration count is a direct measurement rather than
an inference drawn across two independently sorted arrays.

### 4.2 Workloads

Both kernels are exercised over two tables emitted by `model/gen_vectors.py`
into `sw/src/ik_vectors.h`, reported separately.

**Workload A — 48 independent poses (`ik_pose_tbl`).** Half well-conditioned;
one quarter geometrically stressed (wrist near singularity, or elbow near full
extension), which degrades the analytic solver's branch selection; one quarter
selected on the reference model's DLS iteration count in the range [8, 32] with
seeds drawn ±1.2 rad from the solution. Each DLS solve is seeded from its own
perturbation of the answer, so that iteration counts are independent of table
order.

The third group exists because the first two do not produce a tail. DLS is
insensitive to branch degeneracy; it is sensitive to seed distance and to
Jacobian conditioning along the path taken. A table stressed only geometrically
and seeded within 0.25 rad converged in three to five iterations on every pose,
yielding a max/median of 1.77 on hardware against 25 iterations in the 256-pose
software sweep — that is, it did not sample the tail that determines
schedulability. Selection is deliberately frozen on the *unclamped* solver so
that the workload remains fixed across solver changes and improvements remain
measurable.

**Workload B — 50-sample pentagon trajectory (`ik_traj_*`).** A closed pentagon
of circumradius 0.10 m centred at (0.35, 0) in the base plane, traversed at
z = 0.10 m with the tool approach axis directed vertically downward, ten samples
per edge. Tool yaw takes the values 0°, 90°, 0°, 90°, 0° at the five vertices
and is interpolated linearly along each edge, so that orientation reaches each
vertex continuously. Successive samples differ by approximately 12 mm of
translation and up to 9° of yaw.

The work plane is horizontal rather than vertical for a representational
reason: a vertical plane places the pose at pitch = ±π/2, where the RPY
parameterisation used on the bus is singular and a quantised pose no longer
round-trips through `T_to_pose()`. The horizontal configuration yields
`(roll, pitch, yaw) = (π, 0, ψ)`, which is well defined for all ψ.

DLS is chained: each sample is seeded from the previous sample's solution, and
the first from the last, the path being closed. The emitted tables are the
second lap, so that sample 0's seed is a genuine predecessor rather than a cold
start appearing as an isolated outlier. On hardware the chain is fed the
*kernel's* own output rather than the table's, so that quantisation error
integrates around the lap rather than being resynchronised to the reference each
sample; the resulting divergence is reported as chain drift.

The two workloads answer different questions and neither substitutes for the
other. Workload A characterises what a scheduler must survive; Workload B
characterises what the machine does the remainder of the time.

### 4.3 Verification

Three layers, each capable of detecting failures the others cannot:

1. **`model/validate.py`** — the reference model against itself: analytic
   Jacobian versus numerical differentiation, and FK(IK(FK(q))) round-trips over
   all eight branches. 27,008 branch solutions, worst pose error 1.1e-13. No
   HLS involved; isolates kinematic and algorithmic error at the source.
2. **`make -C hls host`** — the kernel sources themselves, compiled in their
   double-precision configuration and executed against the reference vectors on
   the development machine. No Vitis installation required.
3. **`make -C hls csim` / `cosim`** — the same testbenches against `ap_fixed`,
   and subsequently against generated RTL. This is the only layer able to
   detect a fixed-point-specific defect.

Current host regression:

```
  mat_mul_kernel               vectors=665   worst=0.000e+00  pass
  mat_inv (residual A*Ainv-I)  vectors=943   worst=3.736e-05  pass
  fk (position, m)             vectors=384   worst=7.649e-06  pass
  ik_analytic (round trip, pos m) vectors=256 worst=1.513e-05  pass
  ik_analytic (round trip, rot)   vectors=256 worst=2.177e-05  pass
  ik_dls  converged 256/256
          iterations: min=3 median=4 p95=7 max=25 mean=4.42
```

Comparison of DLS output against the analytic solver's joint solution is
invalid and is not performed. IK is multi-valued: `ik_qgold_tbl` records the
analytic `(+,+,+)` branch, whereas DLS converges to whichever branch its seed is
nearest, routinely differing by approximately π in a wrist joint while solving
the commanded pose to tolerance. DLS output is compared against `ik_qdls_tbl`,
the reference model's own DLS result from the same seed. Checked against the
analytic branch, this column reported approximately 3.1 rad of apparent error on
poses solved to within 1e-4 m.

### 4.4 A characterised accuracy limitation

Near the elbow singularity the analytic solver's joint output degrades while
its pose output does not. Since γ = acos(cos γ) has slope 1/|sin γ|, as the arm
approaches full extension or full fold the ±0.5 LSB uncertainty of the Q16.16
square root is amplified into milliradians of joint error. At |sin γ| ≈ 2e-3 the
worst joint deviation reaches 3.2e-2 rad, while the commanded pose is still
attained to 1.5e-5 m and 2.2e-5 rad.

This is a conditioning property rather than an implementation defect: the
elbow-up and elbow-down branches merge in that region, so joint space is
genuinely not well determined there. The testbench asserts the pose round-trip
on every vector without exception and asserts joint agreement only where
|sin γ| ≥ 0.05, reporting the remainder separately.
`model/ik_model.elbow_conditioning()` computes the quantity.

---

## 5. Results

All figures below are from the PS wall-clock instrument described in §4.1, at a
PL clock of 80 MHz. The two IK kernels do not fit simultaneously on this device
(§6), so they were measured in separate bitstreams; the analytic figures come
from an earlier bitstream that predates both the trust region and the
two-workload harness, and are reported here for scale rather than as a
matched-run comparison.

### 5.1 Analytic kernel, Workload A

| path | min | median | p95 | max | max/med |
|---|---|---|---|---|---|
| PL | 10716 ns | 10744 ns | 10781 ns | 10781 ns | **1.00** |
| PS (double, same algorithm) | 32446 ns | 32612 ns | 32732 ns | 37923 ns | 1.16 |

The PL distribution is constant to the resolution of the timer. This serves as
the experimental control: it establishes that neither the bus, the timer, nor
the harness contributes measurable spread, and therefore that spread observed
on the DLS kernel is attributable to the solver.

### 5.2 DLS kernel, Workload A (48 poses, trust region enabled)

| path | min | median | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 48535 ns | 64495 ns | 160692 ns | 513320 ns | 89236 ns | **7.95** |
| PS (double) | 79064 ns | 106864 ns | 298227 ns | 871824 ns | 149556 ns | 8.15 |
| PL iterations | 3 | 4 | 10 | 32 | 5 | 8.00 |
| PS iterations | 3 | 4 | 11 | 32 | 5 | 8.00 |

### 5.3 DLS kernel, Workload B (50-sample trajectory, chained warm start)

| path | min | median | p95 | max | mean | max/med |
|---|---|---|---|---|---|---|
| PL | 32390 ns | 32421 ns | 32464 ns | 32473 ns | 32418 ns | **1.00** |
| PS (double) | 51560 ns | 52150 ns | 52440 ns | 52566 ns | 52006 ns | 1.00 |
| PL iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |
| PS iterations | 2 | 2 | 2 | 2 | 2 | 1.00 |

Every sample converges in two iterations. The observed latency spread across
the lap is 83 ns on 32.4 µs, which is at the resolution of the instrument rather
than a property of the solver. Under this workload the same kernel that exhibits
a max/median of 7.95 on Workload A exhibits 1.00.

The PL and PS iteration counts are identical at every sample, so quantisation
did not displace the convergence path anywhere on this trajectory.

### 5.4 Decomposition of latency

Workload B provides 50 samples at a known and identical iteration count, which
removes the need to pair independently sorted latency and iteration arrays.
Fitting `a + b·iterations` through the trajectory median (2 iterations) and the
Workload A maximum (32 iterations) — two points, no least-squares fit — gives

```
latency_ns  =  361  +  16030 × iterations
```

corresponding to approximately 1282 PL cycles per iteration at 80 MHz on a fixed
overhead of 361 ns. Every other measured row follows from this relation without
further fitting:

| row | iterations | predicted | measured | residual |
|---|---|---|---|---|
| Workload A min | 3 | 48451 ns | 48535 ns | +0.17% |
| Workload A median | 4 | 64481 ns | 64495 ns | +0.02% |
| Workload A p95 | 10 | 160661 ns | 160692 ns | +0.02% |
| Workload A max | 32 | 513320 ns | 513320 ns | 0.00% |
| Workload B, all 50 | 2 | 32421 ns | 32390–32473 ns | ±0.16% |

The harness's independently computed per-iteration rows agree to within 3 ns:
`(361 + 2b)/2 = 16210` against 16209 measured on the trajectory,
`(361 + 4b)/4 = 16120` against 16123 at the Workload A median, and
`(361 + 32b)/32 = 16041` against 16040 at its minimum.

The per-iteration figures warrant care in interpretation. The trajectory's
per-iteration cost (16209 ns) reads *higher* than Workload A's (16123 ns), which
might suggest that its iterations are more expensive. They are not: the
difference is the 361 ns fixed overhead amortised over two iterations rather
than four, and Workload A's own per-iteration minimum of 16040 ns belongs to its
32-iteration pose, where that overhead is nearly absent. The consistency of the
relation across a sixteen-fold range of iteration counts supports the conclusion
that no loop within an iteration is data-dependent.

Three consequences follow:

- AXI4-Lite transaction overhead constitutes 0.6% of a median Workload A solve,
  1.1% of a trajectory sample, and 0.07% of the worst case. For this kernel there
  is no material improvement available on the bus, in contrast to the analytic
  kernel where the same overhead is a substantial fraction of 10.7 µs.
- All observed latency variance is iteration count.
- The two workloads bracket rather than contradict one another. 7.95 and 1.00
  were measured on the same kernel, in the same bitstream, in the same run.

### 5.5 Numerical stability of the chained solve

Cumulative drift of the hardware chain from the reference model's chain peaks at
152 µrad at sample 28 and decreases to 76 µrad by sample 40. The Q16.16 LSB is
15.3 µrad, so the peak corresponds to approximately ten LSB. On a 0.5 m lever
this is of the order of 76 µm of tool position, and it is an order of magnitude
inside the solver's own convergence tolerance.

The drift is bounded and non-monotone, which is consistent with the mechanism:
each sample re-converges to `tol`, re-anchoring the chain to the commanded pose
every control period, so error cannot integrate without limit. The same figure
is 7 µrad in the host build, where the chain and the reference share arithmetic;
the twentyfold increase on hardware is the fixed-point contribution.

This result is one lap on one trajectory and should be read as a measurement
rather than as a bound (§7).

### 5.6 Effect of the trust region

`iks::dls()` bounds `|dq|∞` per iteration at 1.5 rad by successive halving, a
power-of-two operation that is exact in both the double model and the Q16.16
kernel. The same 48 poses, before and after, both at 80 MHz:

| | median | p95 | max | max/med |
|---|---|---|---|---|
| unclamped | 64036 ns / 4 it | 255470 ns / 16 it | 351175 ns / 22 it | 5.48 |
| trust region 1.5 rad | 64495 ns / 4 it | **160692 ns / 10 it** | **513320 ns / 32 it** | 7.95 |

The p95 latency fell by 37%, and over ten independently drawn tables in the
model the count of poses missing a 16-iteration deadline halved. The median is
unchanged, the clamp not firing on a well-seeded solve. The maximum, however,
rose from 22 to 32 iterations, where the model had predicted 32 in both
configurations.

The rise in `max/med` from 5.48 to 7.95 should therefore not be read as the
clamp degrading the solver, nor as confirmation of the model's prediction that
the maximum would be unchanged. This table's extreme is a single marginal pose
at 32 iterations against a cap of 64, and the iteration count of a marginal
solve is not reproducible between double and Q16.16 arithmetic. The discrepancy
is unresolved and is recorded in [STATUS.md](STATUS.md).

### 5.7 Processing-system reference

The PL path is 1.61× faster than the A9 executing the same algorithm in double
precision on Workload B, and 1.66× on Workload A, against approximately 3.0× for
the analytic kernel. An iterative solver spends its time in a loop that the PL
runs at 80 MHz and the A9 at 667 MHz; the PL advantage derives from parallelism
within an iteration, and at this problem size that advantage is modest. The
argument for the PL implementation here rests on determinism rather than on
throughput.

It is also notable that the PS spread (8.15) slightly exceeds the PL's (7.95)
while running the same algorithm. Most of this is the shared iteration-count
spread; the remainder is cache behaviour and double-precision software
floating point. It is not an artefact of the PL implementation.

---

## 6. Implementation constraints

An initial full build of all four kernels exceeded the xc7z020's CARRY4,
DSP48E1 and LUT budgets by roughly a factor of three and failed `place_design`.
The following changes account for most of the reduction:

- `ikm::sincos()`, `atan2_hypot()`, `sqrt()` and `sqrt_acc()` combined
  `#pragma HLS INLINE` with `#pragma HLS UNROLL` on 24–32 iteration
  CORDIC/shift-subtract loops, so each call site received a fully spatial copy
  of the adder chain. The loops were rolled (a fixed trip count still yields
  constant latency) and `INLINE` changed to `INLINE off`, so that call sites
  share one synthesised engine.
- The same over-unrolling affected `dh_step()`'s 3×3 rotation product and the
  wrist-rotation computation in `ik_analytic.cpp`.
- `ik_acc_t` was declared `AP_RND, AP_SAT` although rounding is required only
  once, on the narrowing cast to `ik_real_t`; every accumulation step was paying
  for hardware that could not fire. Changed to `AP_TRN, AP_WRAP`.
- Vitis HLS auto-pipelines small loops absent an explicit directive, which
  silently flattened `dh_step()`'s already-rolled product inside the chain loops
  of `fk()`, `rot03()` and `fk_jacobian()`. Explicit `#pragma HLS PIPELINE off`
  restored sequential execution.
- Three sequential calls to `mm::multiply()` in `iks::dls()` synthesised as
  three instances rather than sharing one (an `ALLOCATION` pragma did not
  enforce sharing in this release); two were routed through a shared loop with
  value-muxed operand staging, arrays of pointers not being synthesisable.

`ik_analytic_kernel` subsequently fits standalone at 86% LUT and 75% DSP.
`ik_dls_kernel` fits standalone at approximately 180 DSP and 41k LUT. Together
they require approximately 345 DSP against the device's 220, so the two are
built and measured in separate bitstreams at the same clock. This is the
operating mode rather than a temporary state; `scripts/build_vivado.tcl
--kernels` selects, and `sw/src/main.c` compiles out whatever is absent, so no
source edit is required to switch.

Comparing across bitstreams built at different clocks is a documented hazard: an
earlier round compared the analytic kernel at 80 MHz against DLS at 40 MHz, and
half the apparent difference was the clock.

The PL clock is 80 MHz. A 100 MHz build placed and routed but missed setup
timing at WNS = −1.582 ns (hold was met, WHS +0.033). HLS schedules against a
10 ns clock with 12.5% uncertainty, but that uncertainty is estimated before
routing exists, and a DSP-heavy design at 75% DSP occupancy does not obtain the
routes assumed. The design closes at approximately 86 MHz, leaving roughly 0.9 ns
of margin at 80 MHz. `build_vivado.tcl` refuses to write an XSA when timing
fails, overridable with `--allow-timing-fail`: a bitstream with negative slack
still programs and still returns answers, occasionally incorrect ones varying
with temperature and with which path lost the race, which is indistinguishable
from a kernel defect.

Any latency quoted from this design should be quoted with its clock. Nothing in
the measurement method depends on the clock value — `ik_driver.c` times with the
PS global timer and `main.c` reports nanoseconds — so PL figures remain correct
and simply scale.

---

## 7. Threats to validity

- **Single-run measurements on a single device.** No figure in §5 is repeated
  across devices, temperatures, or bitstream builds. Placement varies between
  builds and the reported PL figures should be expected to move at the level of
  their observed spread.
- **`max/med` on Workload A rests on one pose.** The numerator is a single
  marginal 32-iteration solve whose count is not reproducible between double and
  Q16.16 (§5.6). The statistic is reported because it is the real-time figure of
  merit, but n = 48 is too small to tune a parameter against; an earlier trust
  region radius of 1.0 rad appeared optimal on one draw and flipped two poses
  from converging to not converging on another.
- **The analytic figures are not a matched run.** They come from an earlier
  bitstream that predates both the trust region and the two-workload harness, so
  there is no analytic row for Workload B at all. §5.1 should be read as
  establishing the instrument's noise floor, not as a like-for-like comparison
  against §5.2.
- **One lap is not a stability proof.** §5.5 shows drift bounded and decreasing
  over 50 samples, which is consistent with the re-convergence mechanism, but a
  pentagon revisits the same five configurations repeatedly and 50 samples is
  not minutes of operation.
- **Workload A's tail is constructed.** A quarter of the table is selected on
  the reference model's iteration count with seeds drawn ±1.2 rad away. It is a
  deliberate worst-case characterisation and is not a sample of any operational
  distribution.
- **The workloads are not exhaustive.** Two paths — one adversarial, one
  benign — do not span the space of trajectories a manipulator may be commanded
  to follow. In particular, no workload here drives the arm through a
  singularity while tracking.

---

## 8. Reproduction

Prerequisites: Vivado and Vitis 2025.2, Python 3 with numpy.

Vitis HLS Classic was removed in 2025.1; `open_project` / `open_solution` and
`vitis_hls -f` no longer describe a supported flow. All HLS steps here use the
configuration-file interface (`v++ -c --mode hls`, `vitis-run --mode hls`),
driven by `hls/Makefile`.

```bash
# 0. reference model and vectors (also emits sw/src/ik_vectors.h, gitignored)
python3 model/validate.py
python3 model/gen_vectors.py

# 1. host regression — no Vitis required
make -C hls host

# 2. HLS: C-simulation, synthesis, IP packaging
make -C hls csim
make -C hls syn
make -C hls reports        # latency and utilisation summary per kernel
make -C hls ip

# 3. hardware — one IK kernel at a time (§6)
vivado -mode batch -source scripts/build_vivado.tcl
#   -tclargs --kernels ik_dls_kernel
#   -tclargs --clk 100
#   -tclargs --no-bit / --allow-timing-fail

# 4. register map, then the application
python3 scripts/gen_regmap.py
vitis -s scripts/build_vitis.py    # must run under the Vitis Python interpreter
```

Override the device with `make -C hls PART=<part>` and
`vivado ... -tclargs --part <part>`.

### 8.1 Generated files

`sw/src/ik_vectors.h` and `hls/tb/vectors/` are generated and gitignored, so a
checkout updates the code that consumes them but not the files themselves.
`sw/src/main.c` guards this with `#error` on `IK_VECTORS_VERSION`, naming the
remedy. A stale header implies a stale workload, not merely a stale constant.

### 8.2 The register-map step is not optional

Vitis HLS assigns AXI4-Lite offsets from argument order and width.
`sw/src/ik_regmap.h` ships with a provisional map so that the tree builds;
`scripts/gen_regmap.py` replaces it with the map extracted from the actual
build. The application performs a write/read-back check at startup and refuses
to report any timing if the map is inconsistent, so a stale header fails audibly
rather than producing plausible but incorrect results.

### 8.3 A measurement defect worth recording

`libxiltimer` arms its default timer instance from an
`__attribute__((constructor))` function in `xiltimer.c`, intended to run before
`main()` without application involvement. On this BSP it did not, or did not
complete in time: every `XTime_GetTime()` call returned the same value and every
latency in the first hardware run read as exactly zero, for both PL and PS
paths — which is what identified the timer rather than either kernel as the
cause. `ik_timer_init()` in `sw/src/ik_driver.c` now calls `XilSleepTimer_Init()`
explicitly and performs a `usleep(1)`.

---

## 9. Further documentation

- [STATUS.md](STATUS.md) — running log: what is built, what has been measured,
  what remains open, and the analysis of how the observed spread might be
  reduced.
- [docs/technical_manual.md](docs/technical_manual.md) — derivation from first
  principles of the kinematics, both solvers, the fixed-point format, and the
  HLS realisation, including the build and measurement defects encountered.
- [docs/timing_methodology.md](docs/timing_methodology.md) — why PS wall clock
  rather than HLS-reported cycle latency is the quantity of interest here.
