# Robotic IK on Zynq-7000 — a timing and real-timeness investigation

Two inverse-kinematics solvers for a 6-DOF arm, built as Vitis HLS IP for a
PYNQ-Z1/Z2 (xc7z020clg400-1) and invoked as kernels from the Cortex-A9, plus
the matrix-multiply and matrix-inversion IPs they are built from.

The point is not to compute IK quickly. It is to compare **a solver whose
latency is a constant against one whose latency is a function of its input**,
on the same silicon, with the same arithmetic, and measure the difference.

| | analytic | DLS |
|---|---|---|
| method | closed form, wrist decoupling | `q ← q + Jᵀ(JJᵀ + λ²I)⁻¹e` |
| loops | all trip counts fixed at compile time | data-dependent exit |
| latency | one number; it *is* the worst case | a distribution |
| fails by | returning `IK_ERR_UNREACH` | running out of iterations |

Measured in the host regression over 256 poses, the DLS solver converges in a
median of **4** iterations, p95 **7**, max **25**. That ~6× spread between
typical and worst case is the whole subject of the investigation: it is the
difference between a solver you can put in a 1 kHz control loop and one you
can only put there if you are willing to bound it.

---

## Layout

```
model/          floating-point golden model + vector generation (Python)
hls/include/    shared headers: numeric types, geometry, CORDIC/sqrt
hls/src/        the four kernels
hls/tb/         C++ testbenches (run under csim AND as a host regression)
hls/cfg/        Vitis HLS 2025.2 config files, one per IP
scripts/        Vivado block design, register-map extraction, Vitis app build
sw/src/         bare-metal driver + timing harness
docs/           timing methodology
```

Four separately packaged IPs:

| IP | top function | notes |
|---|---|---|
| `mat_mul_kernel` | `hls/src/matmul.cpp` | runtime dims ≤6×6, optional transpose on either operand |
| `mat_inv_kernel` | `hls/src/matinv.cpp` | Gauss-Jordan, partial pivoting |
| `ik_analytic_kernel` | `hls/src/ik_analytic.cpp` | fixed latency |
| `ik_dls_kernel` | `hls/src/ik_dls.cpp` | reports iterations used |

The DLS kernel **calls** `mm::multiply()` and `mi::invert()` directly rather
than issuing AXI transactions to the standalone matrix IPs. Same source, two
export targets. Routing every 6×6 product over the interconnect would add
hundreds of cycles per iteration to a loop that runs 4–64 times per solve, and
would make the measurement a measurement of arbitration. The standalone IPs
exist so the matrix primitives can be driven and benchmarked on their own.

---

## Robot and conventions

6R anthropomorphic arm with a spherical wrist, standard (distal) DH:

| i | aᵢ | αᵢ | dᵢ |
|---|---|---|---|
| 1 | 0.05 | +π/2 | 0.15 |
| 2 | 0.30 | 0 | 0 |
| 3 | 0.02 | +π/2 | 0 |
| 4 | 0 | −π/2 | 0.28 |
| 5 | 0 | +π/2 | 0 |
| 6 | 0 | 0 | 0.08 |

Pose is `{x, y, z, roll, pitch, yaw}` with `R = Rz(yaw)·Ry(pitch)·Rx(roll)`.
All values on the AXI bus are signed **Q16.16** (`value = raw / 65536`).

`model/ik_model.py` is the single source of truth for these conventions.
Everything else is checked against it.

---

## Numeric format

Q16.16 was chosen by measurement, not by feel. `model/validate.py` runs the
whole quantised DLS loop in Q16.16, Q14.18, Q12.20, Q14.26 and Q16.32:

```
  fmt       lam     tol     conv  med  p95
  Q16.16   0.02   1e-03  300/300    4    6
  Q14.26   0.02   1e-03  300/300    4    6
  Q16.32   0.02   1e-03  299/300    4    6
```

Wider words buy nothing — convergence is bounded by the damping factor, not by
the word length — so the cheapest format wins. Products accumulate in an exact
Q32.32 accumulator (`ik_acc_t`) and round once on the way out.

One consequence worth knowing: a tolerance of 1e-3 squares to 1e-6, two decades
below the Q16.16 LSB. The DLS convergence test therefore runs in the
accumulator type. Testing it in Q16.16 would quantise to zero and the solver
would report success on its first iteration.

---

## Build

Prerequisites: Vivado + Vitis 2025.2, Python 3 with numpy.

Vitis HLS Classic was removed in 2025.1 — `open_project` / `open_solution` and
`vitis_hls -f` no longer describe a supported flow. Everything here uses the
config-file interface (`v++ -c --mode hls`, `vitis-run --mode hls`).

```bash
# 0. golden model + vectors  (also emits sw/src/ik_vectors.h)
python3 model/validate.py
python3 model/gen_vectors.py

# 1. host regression — no Vitis needed, catches algorithm bugs early
make -C hls host

# 2. HLS: C-sim, synthesis, IP packaging
make -C hls csim
make -C hls syn
make -C hls reports        # latency + utilisation summary
make -C hls ip

# 3. hardware — defaults to mat_mul_kernel + mat_inv_kernel + ik_analytic_kernel
#    (ik_dls_kernel does not fit the xc7z020 standalone yet; see Status)
vivado -mode batch -source scripts/build_vivado.tcl

# 4. register map, then the application
python3 scripts/gen_regmap.py
vitis -s scripts/build_vitis.py
```

Override the device with `make -C hls PART=<part>` and
`vivado ... -tclargs --part <part>`.

### The register map step is not optional

Vitis HLS assigns AXI4-Lite offsets from argument order and width. `sw/src/ik_regmap.h`
ships with a **provisional** map so the tree builds; `scripts/gen_regmap.py`
replaces it with the real one extracted from your build. The application
performs a write/read-back check at startup and refuses to report any timing if
the map is wrong, so a stale header fails loudly rather than producing
plausible-looking wrong answers.

---

## Verification

Three layers, each catching something the others cannot:

1. **`model/validate.py`** — the golden model against itself: analytic Jacobian
   vs numerical differentiation, and FK(IK(FK(q))) round-trips over all eight
   branches. 27,008 branch solutions, worst pose error 1.1e-13.

2. **`make -C hls host`** — the actual kernel sources compiled in their float
   configuration and run against the golden vectors on the development machine.
   No Vitis required. This is where algorithm bugs die.

3. **`make -C hls csim` / `cosim`** — the same testbenches against `ap_fixed`,
   then against the generated RTL.

Current host regression:

```
  mat_mul_kernel               vectors=665   worst=0.000e+00  pass
  mat_inv (residual A*Ainv-I)  vectors=943   worst=3.736e-05  pass
  fk (position, m)             vectors=384   worst=7.649e-06  pass
  ik_analytic (round trip, pos m) vectors=256 worst=1.513e-05  pass
  ik_analytic (round trip, rot)   vectors=256 worst=2.177e-05  pass
  ik_dls  converged 256/256
          iterations: min=3 median=4 p95=7 max=25
```

### A known and characterised limitation

Near the **elbow singularity** the analytic solver's joint output degrades while
its pose output does not. γ = acos(cos γ) has slope 1/|sin γ|, so as the arm
approaches full fold the ±0.5 LSB of the Q16.16 square root becomes
milliradians of joint error. At |sin γ| ≈ 2e-3 the joint delta reaches 3.2e-2 rad
— while the commanded pose is still reached to 1.5e-5 m and 2.2e-5 rad.

This is conditioning, not a defect: the elbow-up and elbow-down branches are
merging there, so joint space genuinely stops being well determined. The
testbench asserts the pose round trip on every vector without exception, and
asserts joint agreement only where |sin γ| ≥ 0.05, reporting the rest
separately. `model/ik_model.elbow_conditioning()` computes the number.

---

## What the timing harness measures

`sw/src/main.c` reports PS wall-clock time around the **whole transaction** —
argument writes, `ap_start`, the poll loop, result reads. That is what a control
loop experiences.

Scope is currently `ik_analytic` (plus the standalone `mat_mul`/`mat_inv` IPs)
only — see Status below for why `ik_dls` is not in the harness yet.

It is deliberately *not* the same number Vitis HLS reports. The HLS latency is
PL cycles between `ap_start` and `ap_done`; it excludes roughly 20 single-beat
AXI4-Lite accesses per solve. For kernels this small that gap is not a rounding
error, and comparing the two is part of the exercise. See
[docs/timing_methodology.md](docs/timing_methodology.md).

---

## Status

Verified: the golden model, all four kernels' algorithms (host regression, all
passing), and the numeric format choice.

A first full run of `make -C hls syn/ip` and `scripts/build_vivado.tcl` against
2025.2 got through `opt_design` and failed `place_design`, needing roughly 3x
the xc7z020's CARRY4/DSP48E1/LUT budget with all four kernels in one bitstream.
Several rounds of targeted fixes brought that down substantially:

- `ikm::sincos()`, `atan2_hypot()`, `sqrt()`, `sqrt_acc()` combined
  `#pragma HLS INLINE` with `#pragma HLS UNROLL` on their 24-32 iteration
  CORDIC/shift-subtract loops, so every call site got a fully spatial copy of
  the wide adder chain. Rolled the loops (fixed trip count still gives
  constant latency) and switched `INLINE` to `INLINE off` so call sites share
  one synthesized engine instead of duplicating it.
- The same over-unrolling pattern hit `dh_step()`'s 3x3 rotation-matrix
  multiply and `ik_analytic.cpp`'s wrist-rotation computation.
- `ik_acc_t` (the Q32.32 MAC accumulator) was declared `AP_RND, AP_SAT` even
  though its own header comment says rounding only ever needs to happen once,
  on the final narrowing cast to `ik_real_t` — every `acc += ...` step was
  paying for round/saturate hardware that could never fire. Switched to
  `AP_TRN, AP_WRAP`; `ik_real_t` itself is untouched, since its `AP_RND`/
  `AP_SAT` are load-bearing (bit-exactness with the golden model, and
  stopping a joint angle from wrapping sign on overflow).
- Vitis HLS auto-pipelines small loops with no explicit directive; this
  silently flattened (fully unrolled) `dh_step()`'s already-rolled matrix
  multiply inside `fk()`/`rot03()`/`fk_jacobian()`'s chain loops. Added
  explicit `#pragma HLS PIPELINE off` to keep them sequential.
- `iks::dls()`'s three sequential calls to `mm::multiply()` synthesized as
  three separate instances instead of sharing one (an `ALLOCATION` pragma
  did not enforce this in this release); routed two of the three through a
  shared loop with value-muxed operand staging (arrays of pointers are not
  synthesizable in Vitis HLS).

Net effect: `ik_analytic_kernel` now fits standalone (86% LUT, 75% DSP).
`ik_dls_kernel` improved from ~3x over budget to ~126% DSP utilisation
standalone but still does not fit — `mi::invert()`'s cost alone (96 DSP)
matches the standalone `mat_inv_kernel`, and appears to be close to the
practical floor without a larger architectural change.

Given that, the default `scripts/build_vivado.tcl` build and the
`sw/src/main.c` benchmarking harness are scoped to `mat_mul_kernel`,
`mat_inv_kernel` and `ik_analytic_kernel` for now, so a real bitstream and
on-target PS/PL numbers are obtainable. `ik_dls`'s HLS sources, driver code
(`ik_driver.c`, `ik_sw_ref.cpp`) and register map are all still in place;
`--kernels` on `build_vivado.tcl` can add `ik_dls_kernel` back in (alone or
otherwise) once its resource fit is resolved or if you want to characterise
it in isolation regardless of fit.

Not yet run to completion on real hardware — co-simulation and the on-target
measurements themselves. The code is staged so that
`make -C hls syn && make -C hls reports` and
`vivado -mode batch -source scripts/build_vivado.tcl` produce them.
