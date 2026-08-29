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
extraction), Phase 2 (adaptive subdivision and output semantics) and Phase 3 (core
redesign) are complete.** The Phase 1 verdict was: continue. Work now heads toward
**Phase 4 (work-stealing parallelism)**.
[`docs/ROADMAP.md`](docs/ROADMAP.md) is the single source of truth for where the
project stands (Japanese).

## What exists today

| Layer | Contents | Spec |
|---|---|---|
| `arith/` | Fixed-width exact integers — no dynamic allocation, no exceptions, no global state | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | Plane-based geometric predicates. **Widths live in the type**, so exceeding a derived bound is a compile error | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | Exact booleans ($\cup$ / $\cap$ / $\setminus$, **$n$-ary**), topology checking, **adaptive subdivision + early-out + constructed-point reuse**, **local BSP**, **contact splitting**, **WNV classification**, a **convex split** at the entry | [`SPEC-phase1.md`](docs/SPEC-phase1.md) – [`SPEC-phase3.md`](docs/SPEC-phase3.md) |

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
| clang-format | Formatting |

## Documentation

The design documents are written in Japanese.

| File | Contents |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **Where the project stands. Start here** |
| [`docs/SPEC-phase3.md`](docs/SPEC-phase3.md) | **Phase 3 spec.** The $n$-ary contract, WNV, core/post-processing split |
| [`docs/IMPL-phase3.md`](docs/IMPL-phase3.md) | Phase 3 implementation notes (**as of completion**). Decisions, rationale, **and the mistakes that were corrected** |
| [`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) | The discussion log behind the Phase 3 spec |
| [`docs/DECISION-core-contract.md`](docs/DECISION-core-contract.md) | How the core contract ($n$-ary, WNV, core/post-processing split) was decided |
| [`docs/SPEC-phase2.md`](docs/SPEC-phase2.md) | Phase 2 spec: split-plane culling, adaptive subdivision, non-manifold output semantics |
| [`docs/IMPL-phase2.md`](docs/IMPL-phase2.md) | Phase 2 implementation notes (**as of completion**). Decisions, rationale, and how added mechanisms moved the detectors |
| [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) | Phase 1 spec: the stitching question, test corpus, abort conditions |
| [`docs/IMPL-phase1.md`](docs/IMPL-phase1.md) | Phase 1 implementation notes. Decisions, rationale, **and the mistakes that were corrected** |
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 spec: bit-width analysis, predicates, test requirements |
| [`docs/IMPL-phase0.md`](docs/IMPL-phase0.md) | Phase 0 implementation notes. **Why it is built this way**, and how the tests were designed to have detection power |
| [`docs/BENCH.md`](docs/BENCH.md) | Benchmark baseline and the Phase 1 / 2 / 3 measurements |
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
| **3** | **Core redesign** ($n$-ary, WNV, local BSP, convex split; single-threaded) | **Complete (2026-08-29)** |
| **4** | **Work-stealing parallelism** | **Next** |
| 5 | Thingi10K full-corpus validation, performance targets | Not started |
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
**No implementation was ever consulted.**

**No GPL/LGPL code (CGAL, Indirect_Predicates, OpenMeshCraft, VCGlib and the like)
has been consulted, quoted, or ported.**
