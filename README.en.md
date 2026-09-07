# Krisite — exact, plane-based geometry for point clouds and meshes

*[日本語版はこちら](README.md)*

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/krisite-logo-dark.svg">
    <img src="assets/krisite-logo.svg"
         alt="Krisite — exact, plane-based geometry for point clouds and meshes"
         width="711">
  </picture>
</p>

---

`Krisite` is a header-only C++20 library for 3D data. The long-term goal is to unify
point-cloud compression, meshing, and exact boolean operations.
**Phase 0 (the arithmetic foundation), Phase 1 (minimal validation of output
extraction), Phase 2 (adaptive subdivision and output semantics), Phase 3 (core
redesign) and Phase 4 (parallelism) are complete.** The Phase 1 verdict was:
continue. Work is now in
**Phase 5 (full Thingi10K validation and performance targets)**.
**CP1 — 500 real-data pairs — has completed: 499 succeeded, 1 failed.**
[`docs/ROADMAP.md`](docs/ROADMAP.md) is the single source of truth for where the
project stands (Japanese).

## What exists today

| Layer | Contents | Spec |
|---|---|---|
| `arith/` | Fixed-width exact integers — no dynamic allocation, no exceptions, no global state | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | Plane-based geometric predicates. **Widths live in the type**, so exceeding a derived bound is a compile error | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | Exact booleans ($\cup$ / $\cap$ / $\setminus$, **$n$-ary**), topology checking, **adaptive subdivision + early-out + constructed-point reuse**, **local BSP**, **contact splitting**, **WNV classification**, a **convex split** at the entry | [`SPEC-phase1.md`](docs/SPEC-phase1.md) – [`SPEC-phase3.md`](docs/SPEC-phase3.md) |
| `par/` | **Persistent thread pool.** Parallelism in the core and the exit. **The output is bit-identical regardless of thread count** | [`SPEC-phase4.md`](docs/SPEC-phase4.md) |

**No floating point, no epsilons.** Every decision is the sign of an exact integer.

```cpp
#include <krisite/krisite.hpp>

using namespace krisite::geom;

IPoint a{0, 0, 0}, b{100, 0, 0}, c{0, 100, 0}, d{7, 9, 3};

// Supporting plane of a triangle. Widths of the normal and offset
// are derived from b automatically.
PlaneD pl = plane_from_triangle(a, b, c);

int s  = side(pl, d);          // -1 / 0 / +1, exact
int o  = orient3d(a, b, c, d); // likewise

// The intersection of three planes carries no coordinates: it stays
// a constructed point in homogeneous form.
HPointD v = intersect3(pl, other1, other2);
int s2 = side(pl, v);          // decided exactly, without division
bool lt = lex_less(v, other_v);
```

Booleans sit on the same exactness.

```cpp
#include <krisite/csg/boolean.hpp>       // binary boolean_op
#include <krisite/csg/soup_boolean.hpp>  // n-ary boolean (the PolySoup path)
#include <krisite/csg/to_mesh.hpp>

using namespace krisite;

mesh::TriMesh A = /* closed, oriented triangle mesh on integer coordinates */;
mesh::TriMesh B = /* likewise */, C = /* likewise */;

csg::BoolStats st;
// `depth` is the octree subdivision depth — a runtime parameter that
// must not affect the semantics.
csg::BoolMesh r = csg::boolean_op(A, B, csg::BoolOp::Union, /*depth=*/2, &st);

// Output vertices remain constructed points (intersections of three planes).
auto t = mesh::check_topology(r.triangles);
assert(t.ok());  // edge-manifold, vertex-manifold, consistently oriented, non-degenerate

// Three or more meshes go through the PolySoup path, which **never rounds an
// intermediate result** — the type is closed under CSG.
csg::BoolOptions opt;
opt.depth = 2;
csg::PolySoup s = csg::boolean(csg::from_mesh(A), csg::from_mesh(B), csg::BoolOp::Union, opt);
s = csg::boolean(s, csg::from_mesh(C), csg::BoolOp::Difference, opt);   // (A ∪ B) \ C

assert(s.source_count() == 3);                // all three inputs are still there
assert(s.sources[0].vertices == A.vertices);  // not one bit changed
const csg::SoupMesh out = csg::to_mesh(s);    // stitch, resolve T-vertices, split, triangulate
```

A `PolySoup` carries **the generation-0 input meshes themselves** (`sources`) plus
**an expression tree for the indicator function** (`indicator`). Classification applies
that indicator to each point's winding number vector $\mathbf{w} \in \mathbb{Z}^n$, so a
chain may cross the same surface twice: $(A \cup B) \setminus B = A \setminus B$, which a
single inside/outside bit cannot express. **Bit widths do not grow along a chain** —
CSG introduces no new planes, so a constructed point stays "the intersection of three
planes chosen from the input" (measured: 141 bits at 1, 2 and 3 stages).

Because every predicate reduces to the sign of a fixed-width integer, and because
there is no allocation, no exception and no global state, the whole thing
parallelises as-is.

### Parallelism does not change a single bit of the output

```cpp
#include <krisite/par/thread_pool.hpp>

krisite::par::ThreadPool pool(8);

csg::BoolOptions opt;
opt.threads = 8;
opt.pool = &pool;              // **reuse the pool** — creating one is expensive

csg::ToMeshOptions tm;
tm.threads = 8;
tm.pool = &pool;               // without this the exit runs single-threaded

const csg::SoupMesh out = csg::to_mesh(csg::boolean(X, Y, csg::BoolOp::Union, opt), tm);
// The bytes of `out` are identical for thread counts 1 / 2 / 4 / 8
```

**Determinism is a requirement, not a by-product.** Each task writes only to its own
indexed slot, and results are joined in index order. The sequential implementation
(`threads = 1`) is therefore a *complete* oracle, which turns verification into its
strongest form — **byte equality** rather than "volume and topology agree"
(888 configurations checked).

**The representative of a set of coincident constructed points is chosen
canonically** as well, so that the homogeneous scale factor does not drift with
ordering. What this buys is uniqueness *as long as the set of constructed points is
unchanged*; it is not a function of the geometry alone (that would need GCD
normalisation, hence division).

### Planes and points are both homogeneous 4-vectors

A plane is stored as `[a, b, c, d]` meaning **`N·x + d = 0`** (with `d = -N·p₁`).
A constructed point is `[x, y, z, w]`, its real position being `V/w`. Under this
convention `side` collapses into a **single 4-dimensional inner product**:

```
sign(w) · sign(a·x + b·y + c·z + d·w)
```

Projective duality shows up directly in the types, and one primitive disappears.
See `docs/SPEC-phase0.md` §3.1.

### Bit widths live in the type

This is the load-bearing design decision. Multiplication widens, and the widening
is visible in the type:

```cpp
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;
```

Writing a predicate's expression fixes the required number of limbs at compile
time, and overflow is prevented by the type system rather than by testing.
The only place limb counts appear as numeric literals is
`include/krisite/geom/widths.hpp`.

## Building

C++20. GCC 13+ / Clang 16+ / MSVC 2022+. Header-only, so putting `include/` on the
include path is enough to use it.

```bash
# development default
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKRISITE_CHECKED_ARITH=ON \
  -DKRISITE_BUILD_TESTS_WITH_GMP=ON
cmake --build build
ctest --test-dir build --output-on-failure

# a single test
ctest --test-dir build -R fixed_int --output-on-failure

# benchmarks
cmake -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DKRISITE_CHECKED_ARITH=OFF -DKRISITE_BUILD_BENCH=ON
cmake --build build-rel
./build-rel/bench/pred_bench   # predicate throughput
./build-rel/bench/soup_bench   # entry / core / exit breakdown (SPEC-phase3 §11)
```

### CMake options

| Option | Default | Meaning |
|---|---|---|
| `KRISITE_COORD_BITS` | 21 | Bit width `b` of input coordinates |
| `KRISITE_CHECKED_ARITH` | ON in Debug | Overflow checking on every operation |
| `KRISITE_BUILD_TESTS` | ON | Build the tests |
| `KRISITE_BUILD_TESTS_WITH_GMP` | **OFF** | GMP differential tests (LGPL, tests only) |
| `KRISITE_BUILD_TESTS_WITH_MANIFOLD` | **OFF** | Manifold oracle (Apache-2.0, tests only) |
| `KRISITE_BUILD_MUTANTS` | OFF | Mutation tests (requires checking OFF) |
| `KRISITE_DEFAULT_ADAPTIVE` | OFF | Make adaptive subdivision + early-out + constructed-point reuse the default for boolean operations |
| `KRISITE_DEFAULT_THREADS` | 1 | Default for `BoolOptions::threads` (used by the parallel-mode CI job) |
| `KRISITE_BUILD_BENCH` | OFF | Build the benchmarks |

`KRISITE_CHECKED_ARITH` is independent of `NDEBUG`. Passing
`-DKRISITE_CHECKED_ARITH=ON` enables the checks even in a Release build.

### Platform validation happens in CI

The matrix of Linux (GCC/Clang) / macOS (Apple Silicon) / Windows (MSVC) ×
`b = 21, 26` is defined in [`.github/workflows/ci.yml`](.github/workflows/ci.yml),
and that file is the authority. There is no need to install those toolchains in a
development container.

**"Must be verified" and "runs on every PR" are different things.** The specs require
the former; how often things run is an operational decision. The policy lives in
[`docs/ROADMAP.md`](docs/ROADMAP.md) ("CI の方針") and has four tiers.

| Tier | When | Contents |
|---|---|---|
| PR gate | Every PR | Formatting + the full test suite on one primary configuration |
| Impact-triggered | PRs touching the relevant paths | Mutation tests, platform matrix, Manifold |
| **main gate** | Push to main | **Every job** |
| **Phase completion** | Closing a phase | **Record that every job was green, with the SHA it ran on** |

**The requirement that every job is green before a phase closes has not been
relaxed.** The tiers change frequency only; nothing was removed from the checks.

CI runs:

| Job | Contents |
|---|---|
| **PR gate** | Full test suite on Linux Clang, `b=21`. **Always runs** |
| Changed scope | Decides from the changed paths which of the jobs below to run |
| Fallback paths | Forces the `__int128` and 32-bit schoolbook paths explicitly |
| **UBSan / ASan** | Undefined behaviour on the **shipping configuration (checking OFF)** |
| GMP differential | $10^7$ arithmetic and $10^6$ predicate comparisons, plus volume identities |
| **Adaptive-subdivision mode** | Flips the default to adaptive and confirms **the Phase 1 test suite still passes in full** |
| Manifold oracle | Compares component count and genus of boolean output against an independent implementation |
| **Mutation tests** | Injects deliberate faults and pins down **both** which tests catch them **and** which combinations do not |
| **Mutation tests with GMP** | Runs only the mutations for which volume is meaningful (2 detected + 4 deliberately not) |
| **ThreadSanitizer** | Races that the determinism check cannot see (two mutations are caught by nothing else) |
| **Parallel mode** | Runs **the entire Phase 1–3 test suite** with `KRISITE_DEFAULT_THREADS=8` |
| $n$-ary / soup path | Exercises the convex split and indicator trees of depth 3 or more |
| Benchmarks (baseline) | Predicate throughput, plus the entry/core/exit breakdown |
| clang-format | Formatting |

## Documentation

The design documents are written in Japanese.

| File | Contents |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **Where the project stands. Start here** |
| [`docs/SPEC-phase4.md`](docs/SPEC-phase4.md) | **Phase 4 spec.** The determinism requirement, how to state a prediction, the dispatch floor |
| [`docs/IMPL-phase4.md`](docs/IMPL-phase4.md) | Phase 4 implementation notes (**as of completion**). **Eight things the measurements revealed** |
| [`docs/SPEC-phase3.md`](docs/SPEC-phase3.md) | Phase 3 spec. The $n$-ary contract, WNV, core/post-processing split |
| [`docs/IMPL-phase3.md`](docs/IMPL-phase3.md) | Phase 3 implementation notes (**as of completion**). Decisions, rationale, **and the mistakes that were corrected** |
| [`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) | The discussion log behind the Phase 3 spec |
| [`docs/DECISION-core-contract.md`](docs/DECISION-core-contract.md) | How the core contract ($n$-ary, WNV, core/post-processing split) was decided |
| [`docs/SPEC-phase2.md`](docs/SPEC-phase2.md) | Phase 2 spec: split-plane culling, adaptive subdivision, non-manifold output semantics |
| [`docs/IMPL-phase2.md`](docs/IMPL-phase2.md) | Phase 2 implementation notes (**as of completion**). Decisions, rationale, and how added mechanisms moved the detectors |
| [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) | Phase 1 spec: the stitching question, test corpus, abort conditions |
| [`docs/IMPL-phase1.md`](docs/IMPL-phase1.md) | Phase 1 implementation notes. Decisions, rationale, **and the mistakes that were corrected** |
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 spec: bit-width analysis, predicates, test requirements |
| [`docs/IMPL-phase0.md`](docs/IMPL-phase0.md) | Phase 0 implementation notes. **Why it is built this way**, and how the tests were designed to have detection power |
| [`docs/BENCH.md`](docs/BENCH.md) | Benchmark baseline and the Phase 1 / 2 / 3 / 4 measurements |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | Third-party components and the mechanisms that keep them out of the distributable |
| [`docs/STYLE.md`](docs/STYLE.md) | Coding conventions |
| [`assets/BRAND.md`](assets/BRAND.md) | Logo and theme colours |
| [`tools/README.md`](tools/README.md) | One-off revision scripts for the documents (**not part of the library**) |

## Licence

**MIT.** This constraint takes precedence over everything else.

The library itself (`include/krisite/`) has no external dependencies whatsoever.
GMP (LGPL) and Manifold (Apache-2.0) are used **only as test oracles**, are disabled
by default, and are never part of the distributable.

[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) lists the third-party
components and, more importantly, the mechanisms that enforce those boundaries at
build time rather than in prose.

## Roadmap

| Phase | Contents | Status |
|---|---|---|
| 0 | Fixed-width exact integers + plane-based predicates | Complete (2026-08-26) |
| 1 | Minimal validation of output extraction (fixed-depth subdivision, single-threaded) | Complete (2026-08-27) — **verdict: continue** |
| **2** | **Adaptive subdivision, constructed-point reuse, output semantics** | **Complete (2026-08-28)** |
| 3 | **Core redesign** ($n$-ary, WNV, local BSP, convex split; single-threaded) | Complete (2026-08-29) |
| **4** | **Parallelism** (core + exit; **determinism required**) | **Complete (2026-08-29)** |
| **5** | **Thingi10K full-corpus validation, performance targets** | **In progress (CP1 complete → CP2)** |
| 6+ | Point-cloud codec, GWN, meshing | Not started |

**Phase 1 was the decision point, and the verdict was: continue** (decided
2026-08-27). None of the abort conditions (`SPEC-phase1.md` §11) were met. The most
dangerous one — *needing a mechanism whose bit-width bound cannot be derived* — had
the widest margin of all. **The fixed-width-integer premise held all the way
through**, and that is the reason to keep going.

**Phase 2 passed CP1 through CP5 and is complete** (2026-08-28). Its completion
criteria are about correctness only; no performance target was set.

**Phase 3 rebuilt the core.** Booleans are normally used in chains, so the
requirement is now to evaluate them as an $n$-ary operation and **chain them without
rounding intermediate results**.

```
from_mesh : TriMesh (integers) → PolySoup   entry: quantise, convex-split, build edge planes
boolean   : PolySoup × … → PolySoup         ★ closed under CSG; n-ary
to_mesh   : PolySoup → TriMesh              exit: stitch, resolve T-vertices, split, triangulate
```

Classification moved from sign vectors to **generalized winding number vectors (WNV)**.
Fragment subdivision moved from over-subdivision by every support plane to a **local BSP**
(raw fragments down to 80.3%, canonicalised fragments to 78.6%). The entry now has a
**convex split**, switchable at runtime against plain triangulation (pieces down to 54.5%,
planes to 64.9%).

**Phase 4 parallelised the core and the exit. Its completion criteria are not about
speed but about "do you understand where the time goes"** (`SPEC-phase4.md` §6.3).
**Determinism is a requirement**, and 888 configurations agree byte for byte.
See [`docs/ROADMAP.md`](docs/ROADMAP.md) (Japanese).

### Measurements (Phase 1 → Phase 2)

Across 22 cases × 3 operations × (fixed depths 0–3 + adaptive subdivision).

| Quantity | Phase 1 | Phase 2 | Why it matters |
|---|---|---|---|
| Leaves | 32,256 | **7,434** | Adaptive subdivision; sets the grain for cell parallelism |
| Fragments (after canonicalisation) | 15,624 | **8,217** | Adaptive subdivision → early-out |
| `intersect3` calls | 1,539,819 | **28,836 (1.9%)** | **Constructed-point reuse** |
| `side` : `intersect3` ratio | 1.26 : 1 | **160.57 : 1** | **Phase 1's top-priority problem is resolved.** Mixed-width `det3` has dropped in priority |
| Wall clock, whole corpus | 220.9 ms | **47.6 ms** | Recorded, not a criterion. **Varies ±15%, so a 10% difference means nothing** |
| Non-manifold outputs excluded | 3 configs | **0** | Contact splitting |

The structural numbers Phase 1 established still hold.

| Quantity | Measured | Why it matters |
|---|---|---|
| Max planes meeting at one point | **3** axis-aligned / **9** with slanted faces | **An axis-aligned-only corpus structurally cannot exceed 3** |
| Rate of value-based vertex merging | up to **44%** | Keying on the plane triple alone is not enough |
| Spatial extent of a merge group | **one cell and its neighbours** | No global sort needed; parallelism closes at cell scope |

### Measurements (Phase 3)

The core and the post-processing are timed separately. **Comparing against published
numbers means first checking what they include** — EMBER's 1.6 ms is the time to
produce the soup, nothing after it.

| Stage | Time | Share |
|---|---:|---:|
| `from_mesh` (entry) | 0.8 ms | 0.7% |
| **`boolean` (core)** | **75.9 ms** | **64.6%** |
| `to_mesh` (exit) | 40.9 ms | 34.8% |

| Mechanism | Effect |
|---|---|
| Local BSP (replacing over-subdivision) | raw fragments **80.3%** / canonicalised **78.6%** |
| Convex split (entry, runtime switch) | pieces **54.5%** / planes **64.9%** / fragments **76.8%** |
| Bit width at chain depth 1 / 2 / 3 | **stays at 141** |

Details in [`docs/BENCH.md`](docs/BENCH.md); the reasoning behind them in
[`docs/IMPL-phase1.md`](docs/IMPL-phase1.md),
[`docs/IMPL-phase2.md`](docs/IMPL-phase2.md) and
[`docs/IMPL-phase3.md`](docs/IMPL-phase3.md).

### Measurements (Phase 4)

**The performance corpus is separate from the correctness corpus**
(`SPEC-phase4.md` §8.2). The existing corpus does about 1 ms of work per boolean,
which is far too little to measure scaling.

| Corpus | Threads | $k_c$ (core) | $k_e$ (exit) | Overall |
|---|---|---:|---:|---:|
| Two spheres (uniform) | 8 | 4.69 | 4.41 | **4.62** |
| Two spheres (uniform) | 16 | 5.48 | 5.45 | **5.52** |
| Cube grid × sphere (skewed) | 16 | 3.15 | 4.30 | **3.44** |
| Small corpus (22 cases × 3 ops) | 8 | 1.59 | 1.11 | **1.41** |

**On a skewed distribution the core saturates at 3.15.** The per-leaf workload varies
too much; **load imbalance sets the ceiling**, not the mechanism.

| Check | Result |
|---|---|
| Determinism (3 corpora × threads 1/2/4/8) | **888 configurations, byte-identical** |
| Agreement with Phase 3 (sequential) | **geometry, topology and volume: 264 / 264** |
| ThreadSanitizer | 30 / 30, zero warnings |
| Prediction vs. measurement | **−0.7% … +1.6%** across 13 configurations |

**The prediction composes $k_c$ and $k_e$ taken from two *separate* runs.** Written
from a single run's breakdown,

$$\text{prediction} = \frac{\text{total}_1}{\text{entry} + \text{core}_1/k_c + \text{exit}_1/k_e}$$

has a denominator that is exactly $\text{total}_{th}$, so **agreement is algebraically
forced** (`SPEC-phase4.md` §6.0). Composing a core-only-parallel run with an
exit-only-parallel run makes the prediction fail if the two interfere.

Details in [`docs/BENCH.md`](docs/BENCH.md); the reasoning in
[`docs/IMPL-phase4.md`](docs/IMPL-phase4.md).

### Real data (Phase 5, in progress)

**What a synthetic corpus shows and what real data shows are different things.**
1,000 models were drawn from Thingi10K under the same conditions as EMBER §5.1
(solid, manifold, no self-intersection, 1,000–100,000 faces), paired up, and run
through all three boolean operations.

**Input acceptance** ($b = 21$, 1,000 models): **1,000 accepted, 0 rejected.**
Quantisation removed 1,618 zero-area triangles across 86 models and **left
$\partial S = 0$ intact in every case**. Zero rejections is the first evidence that
the input constraints (PWN, integer coordinates) are not a practical obstacle.

**CP1 (500 pairs) — complete** (2026-09-06):

| | |
|---|---|
| succeeded / **failed** / halted | **499** / **1** / **0** |
| total compute | **3.56 hours** (8 threads, depth 6) |
| per pair (3 ops + all checks) | median **3.0 s** / p95 93.4 s / max 1,418 s |
| largest pair reached | **196,282 triangles** |
| **models that could declare NSI** | **930 / 1,000 (93.0%)** |

All 500 pairs ran under one identical configuration. Every anomaly is concentrated
in the single failing pair: **across 500 pairs there is not one deficient edge, not
one edge of degree 5 or more, and not one unexplained non-manifold vertex.**

> **These are not performance targets.** No performance target has been set for
> Phase 5 yet (that is CP4's job). This is a record of what real data did.

#### The indicator function was wrong: inside is $w > 0$, not $w \ne 0$

**A self-intersecting closed surface creates regions of negative winding number**,
and the old definition classified them as inside. On `78535` the output fell from
725,650 triangles to **58,976** (12.3× fewer); on `110031` the excess edges went from
94 to **0**. **70 of the 1,000 models (7.0%) self-intersect after quantisation**, so
every failure rate measured before this correction contains the error.

#### The one failure belongs to one model, not to the pair

Seven hypotheses were **refuted by measurement**, not by argument: the input being
non-manifold (its excess-edge count is 0), scale (a same-size pair succeeded while
taking 4× longer), self-intersection as such (the controls self-intersect and
succeed), the NSI declaration, octree cell boundaries (**depth 5/6/7 change the
output by 2.5× while the defect does not move at all**), the known
"same face emitted twice" shape (**0 of the 8 degree-4 edges have two incident
triangles with equal vertex sets**), and contact splitting (**disabling it yields the
same 8 edges and the same triangle indices**).

**`135071` alone — self-unioned, without the other operand — reproduces the defect
exactly.** The non-manifold part is a single simple polyline through 9 vertices, 8 of
them constructed points and 1 an original input vertex.

> The failing pair satisfied a *known* signature (no deficient edges, degree 4,
> unresolved splits). **Looking at the actual triangles showed the signature meant
> something different.** Matching a condition and matching its meaning are not the
> same thing.

#### "Self-intersecting" does not explain the failure

The **70 models that self-intersect after quantisation** were selected **by property,
not by listing identifiers**, and self-unioned: **69 succeeded, 1 failed**, in 0.94
hours total. Self-intersection resolution genuinely runs — **50 of the 70 models
produce regions of negative winding number** (median 18, max 3,742), and `135071`
sits at 44, nowhere near the top.

The **topology of the self-intersection curve** was measured too (no new predicate
and no new bit-width derivation were required): **64 of 70 models have an open curve,
and 63 of those succeed.** An open self-intersection curve is not the explanatory
variable either.

> **The predictions were written down before the run, together with what each
> possible outcome would mean.** The spec side predicted further failures; the
> implementation side predicted 2–8 models. **Both were wrong.** Because "exactly one
> failure ⇒ look for something specific to `135071`" had been decided in advance, the
> next step was fixed the moment the result came in — and "zero failures ⇒ suspect the
> procedure" had already been ruled out by reproducing `135071` on its own first.

#### Where the superlinearity actually lives

The original measurement was $t \propto n^{1.86}$. Re-measured over 500 pairs, it
splits into two relations — and the earlier figures, taken from a handful of pairs,
have been withdrawn:

| Relation | Exponent | Reading |
|---|---:|---|
| $n \to P$ (input triangles → output polygons) | **0.91** | **essentially linear — no superlinearity here** |
| **$P \to t$** | **1.23** | **this is where it is** |
| $n \to t$ | **1.19** | |

**Operation counts localise it further** — these are exactly reproducible and do not
depend on timing noise. Over a 10.9× increase in $n$, **per output polygon**:
ray-versus-triangle tests grow **1.98×**, BSP cut candidates **1.50×**, and the
number of regions **1.04× — it does not move.** What grows is how many triangles one
ray has to examine, which is the natural $O(\sqrt{n})$ behaviour. **The place where a
spatial index pays off is now identified by number, not by argument.**

#### Eight improvements, and the line that separates them

| # | Improvement | Effect | Output | How it is guarded |
|---|---|---|---|---|
| 1 | Remove classification memoisation | memory 8,242 → 1,438 MB | unchanged | **hash equality** |
| 2 | Remove duplicated support-plane computation | 1.18× | unchanged | **hash equality** |
| 3 | 2-D hierarchical grid for ray casting | **$P^{1.31} \to P^{1.03}$** | unchanged | **hash equality** |
| 4 | **NSI** (declared absence of self-intersection) | **median 5.5× faster** | **changes** | topology + exact GMP volume |
| 5 | **Key the T-junction index on (leaf, support plane)** | **142× fewer candidates** | unchanged | **hash equality** |
| 6 | **Reject cell assignments whose support plane misses the box** | $\sum P_\ell^2$ **9–25× smaller** | unchanged | **hash equality** |
| 7 | **Projected-AABB pre-test before ray casting** | classification **1.87–2.0×** | unchanged | **hash equality** |
| 8 | **Split single-source leaves when $P_\ell^2$ is large** | **$\ge$ 168×** on the worst pair | **changes** | topology + exact GMP volume |

**1–3 and 5–7 make the same decision faster, so a hash guards them. 4 and 8 change
how the geometry is cut, so a hash cannot.** Confusing the two cost one wasted
attempt: 4 was checked for hash equality and failed on the very first pair, which is
exactly what a 6.93× reduction in $P$ should do.

**5, 6 and 7 all came from the same question:** *what fraction of the candidates this
stage returns does it actually use?* The answers were 1,050×, 21.7× and 85× waste
respectively, and all three were fixed by putting a cheap exact pre-test in front.
**The effects are measured in operation counts, and all three mechanisms are observed
firing across the 500 real pairs** — a separate check makes sure none of them is idle.

For 8, the threshold was **derived from measurement, not invented**: depth 6/7/8 ×
$P \in \{64,\dots,2048\}$, 18 configurations, with $P = 128$ minimal at every depth
and **the exact volume identical in all 18**. The first version made pairs that *can*
declare NSI **1.75× slower**, because splitting a leaf whose local BSP is already
skipped only adds grid cuts. The exclusion is therefore written as a property —
*"skip when the local BSP is skipped"*, not *"skip when NSI was declared"*.

#### Two broken instruments, recorded

**When the instrument is wrong, normal looks abnormal.** `check_topology` returns
early on an empty mesh and leaves `oriented` false, so **116 of 500 pairs (23.2%) read
as "orientation broken" when all of them were simply an empty intersection.** And a
continuation run was launched with a *different* NSI setting from the 351 pairs
already recorded, which inflated the output 2.2× and was briefly misread as another
mechanism's effect — **the evidence to catch it was already in the recorded columns**,
and the driver now prints its full configuration at startup.

#### CP2 is now executable

**CP2 can never declare NSI** — its population is defined by having
self-intersections. Measured on the same 30 pairs with only the NSI declaration
changed, that costs **2.4×** (matching the 2.2–2.5× ratio in operation counts).

| Pairs | Estimated time |
|---:|---:|
| **100** | **3.0 hours** |
| 500 | 15.1 hours |
| 2,200 | **66.2 hours** |

**5.3× better than the previous 350-hour estimate.** The earlier verdict — *not a
volume that can be run on the working machine* — no longer holds.

## References

| Short | Reference |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **FARMA** | Cherchi, Livesu, Scateni, Attene. *Fast and Robust Mesh Arrangements using Floating-point Arithmetic.* ACM TOG 39(6), 2020. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

**The references are a starting point, not a source: the EMBER / OEBSP papers were not
obtained until just before Phase 3.** Until then the design was derived from the
structure of the problem. What was adopted after reading them (the local BSP, half-open
cell assignment, winding-number classification) and what was derived here is recorded in
[`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) (Japanese).

**No GPL/LGPL code (CGAL, Indirect_Predicates, OpenMeshCraft, VCGlib and the like)
has been consulted, quoted, or ported.**
