# Timing methodology

How to get numbers out of this project that mean something, and which numbers
not to trust.

---

## 1. Three different latencies

There are at least three defensible answers to "how long does the IK take",
and they differ by more than an order of magnitude. Say which one you mean.

### (a) HLS kernel latency — PL cycles, `ap_start` to `ap_done`

From `make -C hls reports`. This is the pure compute cost. It is the right
number for "how much of my PL budget does this consume" and for comparing
pragma variants against each other.

For `ik_analytic_kernel` the report gives a single figure, because every loop
in the design has a compile-time trip count. **That figure is the worst case.**
No input can make it slower.

For `ik_dls_kernel` the report gives a min/max range, because `DLS_ITER` exits
on the residual. The max corresponds to the iteration cap and is pessimistic;
the min corresponds to converging on the first check. Neither is the number you
schedule against — see §3.

### (b) Transaction latency — PS wall clock

What `sw/src/main.c` reports. Includes argument writes, the `ap_start` write,
the `ap_done` poll, and result reads: roughly 20 single-beat AXI4-Lite accesses
per solve, each costing tens of PS cycles through the GP port.

This is the number a control loop actually experiences, and for kernels this
small it can dominate. If (b) is several times (a), the design is
communication-bound and the interesting optimisation is the interface — batching
poses, moving to AXI-Stream, or an `m_axi` master pulling from DDR — not the
arithmetic.

Reporting (a) and calling it the loop latency is the most common way to make an
accelerator look better than it is.

### (c) End-to-end, with the rest of the system running

Not measured here. Add interrupt latency, cache effects from other work, and
DDR contention if anything else is bus-mastering. Always worse than (b), and it
is what actually determines whether a deadline is met.

---

## 2. The measurement setup

**Clock.** The Cortex-A9 global timer, via `XTime_GetTime()`, running at
`CPU_3x2x` — half the CPU clock, ~333 MHz on a -1 grade Zynq-7020. Resolution
~3 ns, against latencies in the microsecond range.

**Polling, not interrupts.** `ik_run()` spins on `ap_done`. An interrupt would
add scheduler latency to the very quantity being measured. The spin is bounded
at one second so a mis-mapped register cannot hang the harness.

**Caches on.** Both I- and D-cache enabled, in `init_platform_stub()`. This is
the realistic configuration and it is what makes the PS software reference
honest — running the comparison with caches off would flatter the PL enormously
and tell you nothing about a real system.

**Same algorithm on both sides.** The PS reference compiles
`hls/src/*.cpp` with `-DIK_USE_FLOAT`. It is not a reimplementation. Any
latency difference is attributable to the target, not to two people having
written two different programs.

The numeric types do differ by construction — PL runs Q16.16, PS runs `double`.
That is the fair comparison, since it is what each target is good at, but it
means the two will not agree bit for bit, and the PS figure carries the cost of
double-precision software floating point on a VFPv3 core. If you want to
separate "fixed vs float" from "PL vs PS", build the PS reference with `float`
and compare all three.

---

## 3. Reading the DLS distribution

The DLS kernel's latency is `fixed_overhead + iterations × per_iteration_cost`.
`iterations` is reported by the kernel in its `iters` register, so you can
reconstruct the per-iteration cost by regressing measured latency against it —
a useful cross-check on the HLS report.

From the host regression over 256 poses (λ=0.02, tol=1e-3):

```
iterations: min=3 median=4 p95=7 max=25
histogram:  3:63  4:128  5:36  6:12  7:5  8:5  9:1  10:2  12:1  18:1  19:1  25:1
```

Three things to take from this:

**The mean is misleading.** Mean 4.42, max 25. Scheduling on the mean
under-provisions by ~6×.

**The tail is thin but real.** 250 of 256 poses finish in ≤8 iterations; four
outliers run 12–25. Those are the near-singular poses. A workload of only
comfortable poses would hide them entirely — which is why
`model/gen_vectors.py` deliberately mixes in wrist-singular and
near-full-extension cases for the on-target set.

**The cap is a design parameter, not a safety net.** `max_iter` converts an
unbounded-latency algorithm into a bounded one by giving up. Whether that is
acceptable depends on what the caller does with `IK_ERR_NO_CONV`. Holding the
previous joint command for one cycle is usually fine; treating a non-converged
result as valid is not.

---

## 4. Sweeps worth running

The parameters are runtime registers precisely so these are cheap:

**λ (damping).** Lower λ converges faster but conditions the inversion worse.
`model/validate.py` shows max|A⁻¹| rising from 44 at λ=0.15 to 8863 at λ=0.01 —
the latter needs 14 integer bits and leaves Q16.16 little headroom. Sweep λ
against both iteration count and residual.

**tol.** The single biggest lever on the tail. At tol=1e-3 everything converges
by iteration 25; at 1e-4 the p95 hits the cap, because damping bounds the
asymptotic rate. If your tail is unacceptable, loosening tol is usually the
cheaper fix than raising the cap.

**Seed quality.** The on-target seeds are the true solution perturbed by
±0.25 rad, which models a tracking loop. Seeding from the previous solution in
a real trajectory will do better; seeding from a fixed home pose will do much
worse. Iteration count is a strong function of this and it is easy to
accidentally measure an unrealistically good case.

**Clock target.** The kernels are constrained at 10 ns. Rebuild at 8 or 7 ns
and watch both Fmax and whether the II of the inner loops degrades.

---

## 5. Things that will bite

**Register map drift.** HLS assigns AXI4-Lite offsets by argument order and
width. Change a signature, rerun `scripts/gen_regmap.py`. The startup check in
`main.c` exists because a stale map produces wrong numbers that look plausible.

**`ap_done` is clear-on-read.** The `AP_CTRL_HS` protocol clears `ap_done` when
you read it. Reading `ap_ctrl` twice and expecting the same value will not work.

**Timer wraparound.** `XTime` is 64-bit at ~333 MHz — no practical wrap. But
`ik_run()` returns the difference as `uint32_t`, which caps a single
measurement at ~12 s. Fine here, worth knowing if you extend it.

**The first call is slower.** Cold I-cache on the driver path. Discard it, or
run a warm-up pass, if you care about the minimum.

**Comparing against the wrong software baseline.** `-O2` double-precision on
one A9 core is a reasonable baseline. `-O0` is not, and neither is a version
that recomputes `sin`/`cos` where the HLS path uses CORDIC. The shared-source
arrangement above avoids the second problem; do not undo it by "optimising" one
side only.
