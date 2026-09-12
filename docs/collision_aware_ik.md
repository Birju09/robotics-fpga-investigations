# Nullspace Collision Avoidance for the DLS Kernel, on the Existing 6-DOF Arm

**Status: proposal. Nothing in this document is built or measured.** Numbers
attributed to this repo come from `README.md` §5, `report.log` /
`report2.log`, and the design log in `hls/include/ik_config.hpp`; everything
else is an estimate and is labelled as one.

**Scope.** Nullspace projection only, inside `ik_dls` only, on the PUMA 560
at 6 DOF. `ik_analytic_kernel` is not touched. No manipulator retarget, no
`IK_DOF` change, no `IK_MAT_MAX` change. (An earlier draft of this document
proposed a 7-DOF Panda port to obtain the nullspace; that is deliberately
abandoned here, and §1 is the replacement argument.)

---

## 1. Where the nullspace comes from on a 6-DOF arm

The nullspace form is

```
dq = J⁺ e  +  (I − J⁺J) z
```

where the first term meets the task and the second moves in directions that
don't disturb it, to first order. `z` will be the gradient of a clearance
potential.

The difficulty is immediate: **for a 6-DOF arm at a full 6-DOF pose task,
`J` is square and `(I − J⁺J)` is zero.** Away from singularities the
solution is locally isolated; there is no motion that improves clearance and
keeps the pose. So a nullspace has to be *created*, and there are exactly
three ways to do it.

### (a) Relax the task — recommended

Take the task Jacobian `J_t` to be `m × 6` with `m < 6`, by dropping task
rows the application does not care about. The nullspace then has dimension
`6 − m`, honestly and exactly:

| m | task | nullspace dim | typical application |
|---|---|---|---|
| 6 | full pose | 0 | current behaviour; the control case (§1b) |
| 5 | position + tool axis, **spin free** | 1 | welding, gluing, dispensing, drilling, most tool work |
| 3 | position only, orientation free | 3 | reaching, pointing, pick approach |

This is not a workaround, it is what the redundancy literature means by
redundancy: a manipulator is redundant *with respect to a task*, not in the
abstract. A 6-DOF arm doing a 5-DOF job is a redundant system, and the free
tool-spin case is overwhelmingly the common one in industrial practice.

**The decisive point for this project: `m` can be a runtime register.** Like
`lambda`, `tol` and `step_max` already are, and for the same reason recorded
in `ik_config.hpp` — it changes the distribution being measured, so it must
be sweepable without re-synthesis. That turns nullspace dimension into a
**swept experimental axis** rather than a property of the hardware. A 7-DOF
arm would have given one fixed nullspace dimension; this gives three, on
silicon that already exists, with `m = 6` as a built-in control. That is a
better experiment, not a compromised one.

### (b) The damped "soft" nullspace at m = 6 — the control case, not the plan

Worth stating because it is a real effect and a falsifiable prediction. With
damping, `J⁺ = Jᵗ(JJᵗ + λ²I)⁻¹` is not a true pseudoinverse, so
`(I − J⁺J)` is *not* zero even for square `J`. Its eigenvalues are

```
λ² / (σᵢ² + λ²)
```

for singular values `σᵢ` of `J`. **Prediction: at `m = 6` the collision term
does essentially nothing, except near singularities where `σ → 0` and the
leak approaches 1.**

**Now measured** (`model/validate.py`, `[Nullspace projector]`, λ = 0.02,
300 random configurations, `‖Nz‖/‖z‖` for random `z`):

| m | median | max |
|---|---|---|
| 3 | 0.726 | 0.986 |
| 6 | 0.021 | 0.859 |

So a full-pose task retains about 2% of the avoidance direction against 73%
for position-only — a 35× difference in authority — and the `m = 6` maximum
of 0.86 is the predicted blow-up at near-singular configurations. The
control behaves as the algebra says it should.

This is a genuinely useful control. If `m = 6` ever shows large avoidance
effects *away* from singularities, the implementation is wrong — probably
the projector. Free to measure, and it falsifies a whole class of bug.
`validate.py` asserts the stronger statement it comes from, `J_t N =
λ²(J_t J_tᵗ + λ²I)⁻¹ J_t`, which holds to 6×10⁻¹⁵.

### (c) Not proposed

- **Exploiting actual kinematic singularities** for a rank-deficient
  nullspace: the nullspace exists only where the arm is degenerate. Useless.
- **Task augmentation / extended Jacobian** (stack the clearance constraint
  as an extra task row, making `J` 7×6 and overdetermined): this solves a
  least-squares compromise between pose and clearance rather than projecting.
  It is architecture (B) from §2 wearing a projector's clothes — the pose is
  no longer met exactly, and the weighting is a tuning surface. Rejected for
  the same reason.

---

## 2. Survey: the ecosystem context, narrowed

(The full five-way survey from the first draft is compressed here to what
bears on the chosen approach. The framing question mentioned "OpenMP"; the
relevant library is **OMPL**, and it is worth stating plainly that **OMPL
does not do IK** — it is a sampling-based planner consuming an opaque
`StateValidityChecker`, with kinematics and collision geometry both supplied
by the embedder, MoveIt in the ROS 2 case. Its contribution here is
architectural: a hard separation between "propose a configuration" and "is
it valid", with a boolean crossing the boundary.)

### What the alternatives are, and why they lose here

- **Reject-and-resample** — what MoveIt 2 actually does. Its
  `kinematics_base` plugins (KDL, TRAC-IK, pick_ik, IKFast) are
  collision-*unaware* by construction; collision-awareness is bolted on via
  `RobotState::setFromIK()`'s `GroupStateValidityCallback`, which calls
  `PlanningScene::checkCollision()`. The solver proposes, the scene vetoes,
  the solver reseeds. **Latency is geometric in the acceptance rate** — the
  tail is unbounded and scene-dependent, which no synthesis report can
  bound. This is the architecture a fixed-latency core must reject, and
  saying why is itself a result.
- **Collision as a weighted cost** — `bio_ik` / `pick_ik` goals,
  `RelaxedIK` / `CollisionIK`. Flexible and graceful, but the pose is no
  longer met exactly and the iteration count is wide-tailed. Note that
  CollisionIK uses a *learned* collision proxy specifically because exact
  distance queries are too slow in the control loop — which is precisely the
  gap a PL distance engine addresses.
- **Hard constraints in a QP** — Drake's `DifferentialInverseKinematics`
  with `MinimumDistanceConstraint`, TSID, newer MoveIt Servo. The only
  formulation that can actually *guarantee* clearance, via a linearised
  inequality `d + (∂d/∂q)ᵀdq ≥ d_safe`. Needs an active-set or ADMM solver
  in the loop; active-set iteration count is data-dependent and
  combinatorially so. Recorded as the ceiling. On a part where the current
  6-DOF DLS kernel needs `IK_BATCH_CORDIC` switched off to place at all,
  this is not a candidate.
- **Batched parallel sampling** — cuRobo. Wrong machine (throughput-bound),
  but the **most useful single reference for the geometry decision**: it
  approximates the robot by *spheres* and the world by *cuboids* precisely
  so that every query has a fixed instruction count, and its FK kernel emits
  sphere centres alongside the end-effector pose. The restriction to
  capsules, cylinders, spheres and boxes is the same instinct, independently
  arrived at.

### Why nullspace projection is the right one for *this* repo

Not just resource cost. The project's central measured result is

```
latency_ns = 423 + 17260 × iterations
```

together with the claim that **all** observed variance is iteration count
and nothing inside an iteration is data-dependent. A nullspace collision
term is closed-form and fixed-cost: it changes the `17260` coefficient and
leaves the *structure* of the claim intact, so it can be re-measured against
the existing methodology unchanged. Reject-and-resample and weighted-cost
approaches would invalidate the methodology instead. That is the strongest
argument, and it is specific to this project rather than general.

### Collision geometry back-ends

FCL (MoveIt 2's default), **Coal** (the 2024 rename of hpp-fcl; used by
Pinocchio, and the choice if analytic distance *derivatives* are wanted),
Bullet, libccd. None is portable to the PL — all are pointer-chasing,
dynamically-allocating, runtime-polymorphic C++. What transfers is the
**closed-form primitive-pair distance formulae**, which for this shortlist
are small and fixed-length, and one of which this repo has already written
(§4).

Worth knowing: `franka_description` ships its coarse self-collision model as
**spheres and cylinders** in `*_sc` sub-links, matching the robot's own
internal check. Not applicable to the PUMA directly, but it establishes that
the primitive set here is a vendor-grade modelling choice rather than a
convenience.

---

## 3. What the existing code already gives you, for free

This is the most important section, because the cost of the plan turns
almost entirely on it. Three properties of the current tree make the 6-DOF
nullspace design substantially cheaper than expected.

### `mm::multiply` already takes runtime dimensions

```cpp
void mm::multiply(const ik_real_t A[IK_MAT_MAX][IK_MAT_MAX], ..., int m, int k,
                  int n, bool ta, bool tb, ...)
```

`m`, `k`, `n` are runtime arguments, clamped, with the destination fully
zeroed first "so the compute loop below can visit only the `m × n` block
actually requested". A `3 × 6` or `5 × 6` task Jacobian needs **no change to
`matmul.cpp` at all** — only the argument values at the call site.

In `ik_dls.cpp` the `MULT` loop currently hardcodes both `m` and `k`:

```cpp
mm::multiply(stageA, stageB, IK_DOF, IK_DOF, mulN[s], mulTa[s], mulTb[s], stageC);
```

so the change is to promote those two literals into the existing
`mulTa`/`mulTb`/`mulN` const-array pattern as `mulM[]` and `mulK[]`. With
`J` holding the task rows in `0 .. m−1`:

- `s = 0`: `A = J_t J_tᵗ` → `(m, 6, m, ta=false, tb=true)`
- `s = 1`: `DQ = J_tᵗ u` → `(6, m, 1, ta=true, tb=false)`

Both read only rows `0 .. m−1` of `J`, which is exactly right.

### `spd::solve` is already latency-independent of `n`

```cpp
//! Pad to 6x6 with identity so every loop runs to IK_MAT_MAX - keeps
//! latency independent of n (cf. mat_inv's MI_INIT).
```

The `m × m` damped system for a reduced task is padded with identity
(pivots = 1, comfortably above `IK_PIVOT_EPS`), `SPD_FWD` zeroes the
right-hand side beyond `n`, and `SPD_OUT` zeroes the result beyond `n`.
Every loop bound depends on the loop index only, never on data. **A reduced
task costs the solve nothing and changes its latency by zero cycles.** This
was written for `mat_mul_kernel`'s narrow products; it happens to be exactly
what task relaxation needs.

### `fk_jacobian` already stores everything a clearance gradient needs

```cpp
ik_real_t zax[IK_DOF][3], org[IK_DOF][3];
#pragma HLS ARRAY_PARTITION variable = zax complete dim = 0
#pragma HLS ARRAY_PARTITION variable = org complete dim = 0
```

Every joint axis `z_i` and frame origin `o_i` is already held fully
partitioned in registers. The geometric Jacobian of **any** point on the arm
— not just the tool — is the two lines already written in `JC_COLS`:

```
J_point[:, i] = z_i × (p_w − o_i)   for joints i preceding the witness link
              = 0                   otherwise
```

So the clearance gradient for a witness point `p_w` on link `k` with
outward unit normal `n̂` is one 3-vector dotted into that:

```
∂d/∂q = n̂ᵗ J_point(p_w, k)
```

Six cross products and six dot products, on data already in registers, with
**no new FK pass and no new CORDIC evaluation**. This is the cheapest
possible place to put a collision gradient, and it is available only because
`fk_jacobian` was written the way it was.

### The one refactor that is genuinely required

`spd::solve()` fuses factorisation and substitution (`spd.hpp:38`). The
nullspace step needs two right-hand sides against the *same* matrix, so
split it:

```
spd::factor(A, n, L, invD)  -> status
spd::substitute(L, invD, n, b, u)
```

The six reciprocals — "the dominant cost" per `spd.hpp`, and the reason
`ikm::recip()` exists — are then paid **once per iteration, not twice**.
The split is behaviour-preserving, independently testable in `tb_spd.cpp`
(extend it to two right-hand sides against one factorisation), and worth
doing on its own merits. **Recommended as the first commit.**

One HLS caveat to expect: `L`, `LD` and `invD` are `ARRAY_PARTITION complete`
today, which is what makes the scalarised indexing affordable ("36 registers
+ muxes is cheaper"). Passing them across a function boundary with
`#pragma HLS INLINE off` risks HLS serialising them through memory ports and
destroying that. Plan on `INLINE` (on) for `factor`/`substitute`, and check
the area report rather than assuming either way.

---

## 4. The collision layer

### Geometry

Approximate each physical link by a **capsule** (segment + radius), and each
static obstacle by a **box (OBB)**, **sphere**, or **capsule**.

The golden model is already half written. `model/ik_model.py:self_clearance()`
already reduces the arm to a segment skeleton with coincident DH frames
merged (`a5 = d5 = 0` puts frames 4 and 5 at the wrist centre) and already
computes segment-segment distance in `_seg_dist`. Capsule distance is
segment distance minus the two radii, exactly, with no approximation to
argue about.

**Correction to an earlier draft of this document, found while building
`model/collision.py`.** It is *not* true that a capsule is "exactly that
skeleton with a radius added", and the difference matters for the gradient.
A segment running from `o_{i-1}` to `o_i` is **not a rigid body**: `o_{i-1}`
is fixed in frame `i-1`, `o_i` is fixed in frame `i`, and the segment
between them shears as `θ_i` moves. It therefore has no single point
Jacobian, and the plan in §5 depends on there being one.

Decomposing each DH step into its two limbs fixes this exactly. Since
`T_{i-1,i} = Rz(θ_i) Tz(d_i) Tx(a_i) Rx(α_i)`, in frame `i`:

```
o_i     = (0, 0, 0)                               the joint centre
P_i     = (-a_i, 0, 0)                            after Tz, before Tx
o_{i-1} = P_i - d_i · (0, sin α_i, cos α_i)
```

All three are **constant** — `θ_i` only rotates frame `i` itself. So the
polyline `o_{i-1} → P_i → o_i` is rigid in frame `i`, giving two capsules
(the `d` limb along `z_{i-1}`, the `a` limb along `x_i`), each with one point
Jacobian whose non-zero columns are the joints preceding that frame. One
witness point, one gradient, straight out of the `zax[]`/`org[]` registers —
which is what §3 claimed and what the naive skeleton would not have
delivered.

This is also how the ecosystem models it: `franka_description`'s coarse
`*_sc` geometry is sphere and cylinder sub-links rigidly attached to a link
frame, not a skeleton polyline.

On the PUMA 560 this yields **6 capsules** (step 5 is fully degenerate and
emits none; `a1 = d2 = a4 = a6 = 0` drop one limb each from steps 1, 2, 4
and 6).

That function has also already earned its keep: it retired the original
pentagon trajectory, which crossed itself to within 4 µm (see the header
comment in `gen_vectors.py`). Adding radii converts its own disclaimer —
"this is NOT a collision check... no link cross-section, tool, or pedestal
housing geometry is modeled anywhere in this project" — into an actual
check.

Note also `gen_vectors.py:22`: the PUMA's asymmetric `qlim` is currently
"doing the job a collision model would otherwise have to". Once a real
clearance term exists, that implicit dependency should be stated in
`robot.py`'s provenance rather than relied upon.

### Primitive kernels

| query | formulation | divisions | trip count |
|---|---|---|---|
| sphere–sphere | `‖c₁−c₂‖ − r₁ − r₂` | 0 | fixed |
| capsule–capsule | segment–segment, clamped `(s,t)` | 1 | **fixed** |
| sphere–capsule | point–segment, clamped `t` | 1 | fixed |
| sphere–box | point clamp in box frame | 0 | fixed |
| capsule–box (OBB) | segment vs OBB, clamp in box frame | 0 | fixed |
| capsule–cylinder | no clean closed form | — | — |

**Avoid capsule–cylinder.** Approximate a cylinder obstacle by a capsule
(conservative if the radius is kept — the spherical caps only add volume) or
by a box. Keeping the primitive-pair set to four cases bounds the mux depth
of the dispatch.

Three notes specific to this codebase:

1. **Keep it branchless and the latency stays a constant.** `_seg_dist`'s
   structure is a cascade of `if` on clamped parameters; in HLS these must
   become `select` muxes with every branch *evaluated* and the result chosen.
   `iks::dls()` already establishes this idiom twice, deliberately — the
   trust-region shift ("shift count is always computed and applied... so
   per-iteration latency stays data-independent") and the singular flag
   ("flagged rather than broken out of: avoids adding another data-dependent
   exit to a kernel whose latency is the thing under measurement"). Follow
   that precedent exactly. Result: **the collision layer adds a compile-time
   constant to the iteration cost**, which is the property the whole project
   is organised around.
2. **The division reuses existing hardware.** `_seg_dist`'s
   `denom = ae − b²` needs one reciprocal; `ikm::recip()` is already built
   and already validated against the `ap_fixed` divider in `tb_spd.cpp`. The
   near-parallel-segments guard has an established idiom in `IK_PIVOT_EPS`.
3. **Squared clearances must live in `ik_acc_t`.** A 1 mm margin is 65 LSBs
   in Q16.16 — fine — but its *square* is 4×10⁻³ LSB, i.e. zero. This is the
   same trap the repo already documented for `tol_sq`: "1e-6 at default is
   two decades below the Q16.16 LSB and would quantise to zero, reporting
   convergence on iteration 1". Get this wrong and the kernel reports
   clearance everywhere. Compare `d² ≥ d_safe²` in Q32.32 and take
   `ikm::sqrt_acc` only where an actual distance is needed.

### Broad phase: a compile-time pair list

No runtime broad phase. Emit the pair list at generation time, exactly as
the ecosystem does — MoveIt's SRDF `<disable_collisions>` is a static
precomputed "never check these" set:

- `model/gen_geometry.py` emits `IK_CAP_*` and `IK_COLL_PAIRS`, derived from
  `RobotModel` — **done**;
- obstacles arrive over AXI4-Lite as a fixed-size array of primitives with a
  type tag and a count, with `IK_MAX_OBSTACLES` setting the area/generality
  trade the way `IK_DLS_MAX_ITER` does.

Pruning needs **two** rules, not one, and the second was not anticipated:

1. **Shared node** — the rule `self_clearance()` states. Two capsules
   meeting at a joint have zero segment distance there by construction.
   Adjacency must be decided on *deduplicated skeleton nodes*, not on frame
   index: a frame-index gap looks sufficient and is not, because `a5 = d5 =
   0` means frames 4 and 6 are two apart yet physically share the wrist
   centre. Any arm with a degenerate DH step has this hole.
2. **Unseparable** — the skeleton path between two capsules is shorter than
   the sum of their radii, so their surfaces overlap *even fully extended*:
   no configuration clears them. `self_clearance()` needs no analogue
   because at zero radius the floor for a once-separated pair is the
   intervening limb's length, small but positive; adding radii turns that
   floor negative wherever a link is shorter than the surrounding parts are
   thick. On the PUMA exactly one pair trips it — `link3_d/link4_d`,
   separated by the 20.3 mm `a3` elbow offset while the radii sum to 125 mm.
   Without this rule that pair is permanently the minimum at −0.105 m and
   masks every real collision. The physical arm has no gap there either; the
   elbow is one casting.

Rule 2 is radius-dependent, so changing `LINK_RADIUS` can change the pair
list — `capsule_pairs(explain=True)` reports what was dropped and why, and
the generated header records it.

Result on the PUMA: **6 capsules, 15 candidate pairs, 6 pruned, 9 checked.**
Obstacles would add `6 × IK_MAX_OBSTACLES`. The self-collision half is
nearly free; `IK_MAX_OBSTACLES` is the knob that decides the cost.

**The radii are chosen, not cited.** No published cross-section table for the
PUMA 560 is used anywhere in this project. `model/collision.py:LINK_RADIUS`
states the consequence and `validate.py` repeats it: a *negative* clearance
is a real finding, a positive one is not a proof of clearance. This is
`self_clearance()`'s own caveat moved from "no radius at all" to "an assumed
radius", which is an improvement but not a measurement.

### Return witness points from the start

The check-only kernel needs only `d_min`. The steering kernel needs the
**witness pair** — the two closest points and the unit vector between them —
because that is what turns a distance into a joint-space gradient. Design
`coll::min_distance()` to return `(d, p_a, p_b, link_index)` even for
Stage 0. Retrofitting it later means rewriting the reduction.

---

## 5. The nullspace DLS step

### Never form the projector

The naive `(I − J⁺J)` is a 6×6 product requiring `J_tᵗ A⁻¹ J_t`. Don't build
it. Use:

```
dq = J_tᵗ A⁻¹ e_t  +  [ z − J_tᵗ A⁻¹ (J_t z) ]        A = J_t J_tᵗ + λ²I_m
```

Per iteration this adds:

- `v = J_t z` — one `m × 6` mat-vec;
- `w = A⁻¹ v` — **one `spd::substitute()` against the factors already
  computed** for the task term. No second factorisation, no extra
  reciprocals;
- `Nz = z − J_tᵗ w` — one `6 × m` mat-vec.

No 6×6 matrix is ever formed, and the two mat-vecs are far cheaper than the
`A = J_t J_tᵗ` product already in the loop.

### It extends the existing shared call site

`ik_dls.cpp`'s `MULT` loop exists specifically so that both products per
iteration share one multiplier bank:

```
//! Looping over a shared call site (rather than two straight-line
//! calls) is what makes Vitis HLS reuse one multiplier bank instead
//! of synthesizing two
```

The two new mat-vecs are two more stages of that same loop — `s = 0..3`
instead of `s = 0..1` — so they inherit the sharing rather than
instantiating new multipliers. This is the single most important
implementation instruction in this document: **extend `MULT`, do not add
straight-line `mm::multiply` calls.**

Suggested stage table, with `j` the task-row count:

| s | product | m | k | n | ta | tb |
|---|---|---|---|---|---|---|
| 0 | `A = J_t J_tᵗ` | j | 6 | j | – | ✓ |
| 1 | `DQ = J_tᵗ u` | 6 | j | 1 | ✓ | – |
| 2 | `v = J_t z` | j | 6 | 1 | – | – |
| 3 | `Nz_partial = J_tᵗ w` | 6 | j | 1 | ✓ | – |

`spd::factor()` runs after `s = 0`; `spd::substitute()` runs twice — once
after `s = 0` for `u`, once after `s = 2` for `w`.

### Task-row selection

- **`m = 3`, position only.** Rows `0..2` of `J` and `e` as they already
  stand. **Zero extra arithmetic.** Nullspace dimension 3 — maximum
  avoidance authority. Recommended for Stage 1.
- **`m = 5`, free tool spin.** Needs the orientation error and the angular
  Jacobian rows rotated into the tool frame before dropping the spin
  component, since the free axis is the tool approach axis `z₆`, not a base
  axis:
  ```
  e_o_tool = Rcᵗ e_o          (keep components x, y)
  J_w_tool = Rcᵗ J_w          (keep rows 0, 1)
  ```
  The `3×3 · 3×6` product is one more `MULT` stage. Task rows become
  `[Jv (3); (Rcᵗ J_w) rows 0,1]` and `e_t = [ep (3); e_o_tool x,y]`.
  Recommended as Stage 2 — it is the practically interesting mode.
- **`m = 6`.** Unchanged from today, plus the §1b soft-nullspace control.

### The potential and its gain

Activate the repulsion only inside a threshold, and make both the threshold
and the gain runtime registers — same argument as `lambda`, `tol` and
`step_max`: they change the distribution being measured, so they must be
sweepable without re-synthesis.

```
z = Σ_pairs  w(d) · (∂d/∂q)ᵗ        w(d) = 0 for d ≥ d_thresh
```

Start with the single *minimum*-clearance pair rather than a sum over all
pairs — one witness point, one gradient, one reduction. It is the cheapest
possible version, it captures most of the benefit, and it makes the
reduction in `coll::min_distance()` do double duty. Summing over all
violating pairs is a later refinement whose cost is linear in the pair
count.

**The existing trust region applies to the whole step and should not be
changed.** `DLS_STEP_MAG` / `DLS_STEP_SHIFT` clamp `|dq|_∞` by power-of-two
halving, and the reason recorded for power-of-two is that it "is exact in
both builds, so the float reference and the Q16.16 kernel follow the same
trajectory through the clamp". Clamping the *combined* task-plus-nullspace
step preserves that, and also bounds the repulsion gain's worst case for
free. Do not add a second clamp on `Nz`.

---

## 6. Measurement design

Four things must be reported for this to be a result rather than a demo. Two
of them are traps this repo has already fallen into and documented.

1. **Refit the per-iteration constant and validate on a held-out run.** The
   claim to defend is that `latency_ns = c₀ + c₁ × iterations` still holds
   with a new `c₁`. The existing standard — fit on one run, predict a second
   to within 0.17% without refitting (README §5.2) — is the bar.
2. **Nullspace dimension is the swept axis.** Report the full
   `m ∈ {3, 5, 6}` sweep. `m = 6` is the control and should show near-zero
   avoidance authority away from singularities (§1b); if it doesn't, the
   projector is wrong. This is the experiment the 6-DOF framing makes
   *possible* rather than merely adequate — a fixed 7-DOF arm would have
   given one data point where this gives a curve.
3. **Do not evaluate the collision-aware solver on a workload selected by
   the collision-unaware one.** `ik_config.hpp`'s trust-region section is an
   extended warning about exactly this failure: `gen_vectors.py` builds its
   tail block by keeping poses the *unclamped* solver converges on, so "any
   solver change can only lose poses from a set defined by its predecessor
   succeeding on them... frozen is not the same as neutral", and the
   480-pose table consequently understates the clamp. The collision workload
   must be sited **independently** — poses whose ordinary DLS solution
   violates clearance — or it will understate the new term the same way.
4. **Separate "avoided" from "never in collision".** Report the fraction of
   poses where the unmodified solver's answer violated clearance and the new
   one's does not, against the fraction where neither did. Without that
   split, a collision term that does nothing is indistinguishable from one
   that works.

Two secondary effects to watch, both plausible and neither obvious in sign:

- **Iteration count.** The repulsion fights the task term, so convergence
  may slow; at `m = 3` the task is easier, so it may speed up. The net is a
  measurement. This also means the max/median ratio — the project's headline
  number — will move, and in an unknown direction.
- **New non-convergence mode.** A clearance potential adds local minima. The
  repo already records that ~40% of non-convergent poses sit at genuine
  stationary points of `‖e‖²` with plain DLS; a repulsion term can only add
  to that. `IK_ERR_NO_CONV` will need to be distinguishable from "converged
  but still in collision", so the new status word needs both bits.

---

## 7. Resource reality

Dropping the 7-DOF port turns this from a blocking problem into a tuning
problem — but the part is still tight, and Stage 0 exists to find out how
tight before committing.

From `report.log` / `report2.log` on the xc7z020 (220 DSP48E1, 53,200 LUT):

| build | DSP | LUT | note |
|---|---|---|---|
| `mat_mul_kernel` | 24 (10%) | 5,773 (10%) | |
| `mat_inv_kernel` | 96 (43%) | 12,374 (23%) | |
| (third block) | 166 (75%) | 46,131 (86%) | |
| `ik_dls_kernel`, early | 434 (197%) | 66,369 (124%) | did not place |
| `ik_dls_kernel`, after II tuning | 278 (126%) | 42,776 (80%) | still did not place |

The current tree fits only after the full `ik_config.hpp` campaign — rolled
multipliers, raised `II`, and `IK_BATCH_CORDIC` off because it overflowed
the part "by 293 LUTs". `ik_types.hpp` records the DLS kernel at 96% LUT.
README §6 records that the two IK kernels cannot be co-resident at all —
which no longer matters here, since `ik_analytic_kernel` is out of scope and
the comparison bitstream is unaffected.

What this plan adds, in descending order of expected cost:

1. the capsule distance engine and witness reduction (**dominant**, and
   scales with `IK_MAX_OBSTACLES`);
2. the point Jacobian and gradient dot product (cheap — registers already
   partitioned, §3);
3. two extra `MULT` stages (cheap — shared multipliers, mat-vecs not
   mat-mats);
4. one extra `spd::substitute()` (cheap — no reciprocals);
5. the `factor`/`substitute` split (should be ~neutral; could regress if
   HLS stops partitioning `L` — watch it, §3).

Items 2–5 are plausibly affordable. Item 1 is the open question, and the
knobs against it are, in order: `IK_MAX_OBSTACLES` (self-collision only at
3 pairs is nearly free — obstacles are what cost), the `II` on the distance
engine's product loops, and the primitive set (sphere-only obstacles are
dramatically cheaper than OBBs).

**This is why Stage 0 is unconditionally first.** A standalone
`coll_dist_kernel`, synthesised and timed on its own, tells you the engine's
area and latency *before* any integration work is committed. If it doesn't
fit alongside `ik_dls`, that number is still a publishable result — "what a
fixed-point capsule distance query costs in PL" — and the fallback (a larger
part, e.g. ZU3EG / Kria K26) is a decision made on evidence rather than
optimism.

---

## 8. Verification

Extend the existing three layers in place; don't parallel them.

1. **`model/validate.py`** — capsule-capsule distance against dense
   brute-force sampling of both segments, and the analytic clearance
   gradient `∂d/∂q` against numerical differentiation of `self_clearance()`.
   This mirrors what the layer already does for the Jacobian (analytic vs.
   numerical differentiation) and is where a sign error in the repulsion
   direction will die. Add one more check specific to this design: **verify
   numerically that `J_t · Nz ≈ 0`** for random `z` at `m = 3` and `m = 5`,
   and that it is *not* ≈ 0 at `m = 6` only by the §1b factor. That single
   assertion validates the projector algebra independently of any geometry.
2. **`make -C hls host`** (`-DIK_USE_FLOAT`) — the collision kernel and the
   nullspace step compiled in double against golden vectors. Per CLAUDE.md,
   "this is where algorithm bugs die"; the branchless `select` cascade, the
   witness-point reduction, and the `MULT` stage table all belong here.
3. **`csim` / `cosim`** — the only layer that can catch the
   Q16.16-specific failures, of which this design has three named
   candidates: the squared-clearance underflow (§4, note 3),
   `ikm::recip()`'s behaviour on a near-zero `denom` for near-parallel
   segments, and whether the `factor`/`substitute` split is bit-exact
   against fused `spd::solve()` (it should be; prove it).

One genuinely new vector set is needed: **degenerate geometry** —
zero-length segments (coincident DH frames, which `self_clearance()` already
merges), exactly parallel segments, and exactly touching capsules (`d = 0`).
These are where every segment-distance implementation breaks, and the
fixed-point build will break differently from the double build.

Finally, the default must stay `m = 6` with the collision gain at zero, so
that **every existing measurement in README §5 remains reproducible
bit-for-bit** on the new build. Anything else and the new work invalidates
the baseline it is supposed to be compared against.

---

## 9. Order of work

Branch: `nullspace-collision-ik`.

1. **DONE — `spd::solve()` → `spd::factor()` + `spd::substitute()`**, plus
   `spd::solve_n()` sharing one substitution instance across right-hand
   sides. `tb_spd.cpp` asserts bit-exactness against the fused version:
   `factor+substitute == solve` over 3240 values, `solve_n == solve` over
   6480, zero differing.
2. **DONE — task relaxation, no collision.** `m` and `k` promoted to the
   `mulM[]`/`mulK[]` table at the `MULT` call site; `task_dim` is a runtime
   register and an AXI4-Lite argument. `IK_TASK_FULL` asserted bit-identical
   to the pre-existing path over 1536 values; `IK_TASK_POS` converges 256/256
   on the vector table; `task_dim = 5` refused with `IK_ERR_BADDIM`.
   Position-only iterations come in at min 2 / median 3 / p95 4 / max 5,
   against 3 / 4 / 6 / 26 for the full pose.
3. **DONE — capsule distance and clearance gradient in the golden model.**
   `model/collision.py`; `LINK_RADIUS` on the robot; `IK_CAP_*` /
   `IK_COLL_PAIRS` emitted by `gen_geometry.py` (geometry header version 2).
   `validate.py` checks: closed-form distance never exceeds brute force over
   2360 segment pairs including degenerate cases (0.0), agreement with
   `ik_model._seg_dist` (0.0), analytic vs numerical `∂d/∂q` (6.9×10⁻¹¹ over
   370 of 400 samples, the rest skipped at genuine non-differentiabilities),
   and the damped-projector identity (6×10⁻¹⁵).
4. **NEXT — Stage 0, `coll_dist_kernel`** as a standalone packaged IP,
   synthesised and timed. *This is the go/no-go for the part.* Also closes
   STATUS.md's open item that the random 48-pose table has never been checked
   against the self-collision proxy (~5% of a uniform `qlim` sample fails
   it).
5. **Stage 1 — nullspace repulsion at `m = 3`**, single minimum-clearance
   pair, and the four-part measurement of §6.
6. **Stage 2 — `m = 5`** free-tool-spin mode, and the `m ∈ {3,5,6}` sweep.
7. Optional refinement: sum over all violating pairs rather than the minimum
   one; joint-limit avoidance as a second `z` term (free once the nullspace
   plumbing exists, and worth having given that `qlim` is currently
   substituting for a collision model).

Steps 1–4 each have standalone value and none of them commits to the rest.
Nothing through step 3 has been near a synthesis tool: all of it is the host
regression and the Python model, so **every area and latency claim in §7
remains unverified.**

---

## References

**Redundancy resolution / nullspace**
- Siciliano & Slotine, *A general framework for managing multiple tasks in
  highly redundant robotic systems*, ICAR 1991 — the task-priority
  formulation.
- Maciejewski & Klein, *Obstacle avoidance for kinematically redundant
  manipulators in dynamically varying environments*, IJRR 1985 — the
  canonical nullspace obstacle-avoidance paper.
- Chiaverini, *Singularity-robust task-priority redundancy resolution*,
  T-RA 1997 — on why the damped projector behaves as §1b describes.

**ROS 2 / MoveIt 2 IK solvers**
- [TRAC-IK in MoveIt 2](https://moveit.picknik.ai/main/doc/how_to_guides/trac_ik/trac_ik_tutorial.html)
- [pick_ik in MoveIt 2](https://moveit.picknik.ai/main/doc/how_to_guides/pick_ik/pick_ik_tutorial.html)
- [`bio_ik`](https://github.com/PickNikRobotics/bio_ik) ([original, TAMS](https://github.com/TAMS-Group/bio_ik))
- [GSoC 2023: MoveIt Servo and IK Benchmarking](https://moveit.ai/moveit/benchmarking/inverse%20kinematics/servo/2023/11/21/GSoC-2023-MoveIt-Servo-and-IK-Benchmarking.html)
  — the `ik_benchmarking` tool
- [Planning Scene tutorial](https://moveit.picknik.ai/main/doc/examples/planning_scene/planning_scene_tutorial.html)
  — `isStateValid`, and the callback mechanism behind collision-aware
  `setFromIK` ([discussion](https://github.com/ros-planning/moveit/issues/639))

**Collision geometry**
- [Coal](https://github.com/coal-library/coal) — FCL extension, the 2024
  rename of hpp-fcl · [docs](https://docs.ros.org/en/jazzy/p/coal/index.html)
- [`franka_description`](https://github.com/frankarobotics/franka_description)
  — the `--with-sc` coarse self-collision model (spheres + cylinders), as
  precedent for the primitive set
- [panda_moveit_config SRDF collision model](https://github.com/moveit/panda_moveit_config/pull/35/files)
  — precedent for a static precomputed pair list

**Parallel / batched collision-aware IK, for the geometry decision**
- cuRobo: [paper](https://curobo.org/reports/curobo_report.pdf) ·
  [arXiv:2310.17274](https://arxiv.org/html/2310.17274v2) ·
  [collision world representation](https://curobo.org/get_started/2c_world_collision.html) ·
  [robot sphere configuration](https://curobo.org/tutorials/1_robot_configuration.html)
