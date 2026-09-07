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

**Side-of-plane, orientation and intersection decisions are made in fixed-width exact
integer arithmetic** — no floating point, no tolerances.

The project is currently in **Phase 5 (Thingi10K validation and performance
targets)**. [`docs/ROADMAP.md`](docs/ROADMAP.md) is the single source of truth for
where it stands (Japanese).

## What exists today

| Layer | Contents | Spec |
|---|---|---|
| `arith/` | Fixed-width exact integer arithmetic (no allocation, no exceptions, no global state) | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | Plane-based geometric predicates. **Widths live in the type**, so exceeding a bound is a compile error | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | Exact booleans ($\cup$ / $\cap$ / $\setminus$, **$n$-ary**), topology checks, adaptive subdivision, local BSP, contact splitting, winding-vector classification, convex splitting at the entry | [`SPEC-phase1.md`](docs/SPEC-phase1.md) – [`SPEC-phase3.md`](docs/SPEC-phase3.md) |
| `par/` | Persistent thread pool. **Output is bit-identical regardless of thread count** | [`SPEC-phase4.md`](docs/SPEC-phase4.md) |

## Usage

```cpp
#include <krisite/krisite.hpp>

using namespace krisite::geom;

IPoint a{0, 0, 0}, b{100, 0, 0}, c{0, 100, 0}, d{7, 9, 3};

// Supporting plane of a triangle; the widths of N and d follow from b automatically
PlaneD pl = plane_from_triangle(a, b, c);

int s = side(pl, d);          // -1 / 0 / +1, exact
int o = orient3d(a, b, c, d); // likewise

// The intersection of three planes carries no coordinates: it stays a homogeneous point
HPointD v = intersect3(pl, other1, other2);
int s2 = side(pl, v);         // exact, with no division
```

Booleans sit on the same exactness.

```cpp
#include <krisite/csg/soup_boolean.hpp>  // n-ary boolean (PolySoup path)
#include <krisite/csg/to_mesh.hpp>
#include <krisite/par/thread_pool.hpp>

using namespace krisite;

mesh::TriMesh A = /* closed, oriented triangle mesh on integer coordinates */,
              B = /* ditto */, C = /* ditto */;

par::ThreadPool pool(8);
csg::BoolOptions opt;
opt.depth = 6;                 // octree depth: a runtime parameter, semantics unaffected
opt.threads = 8;
opt.pool = &pool;              // **reuse the pool** — creating one is expensive

// Intermediate results are never rounded (the type is closed under CSG)
csg::PolySoup s = csg::boolean(csg::from_mesh(A), csg::from_mesh(B), csg::BoolOp::Union, opt);
s = csg::boolean(s, csg::from_mesh(C), csg::BoolOp::Difference, opt);   // (A ∪ B) \ C

assert(s.source_count() == 3);                // all three inputs are still there
assert(s.sources[0].vertices == A.vertices);  // not one bit has changed

csg::ToMeshOptions tm;
tm.threads = 8;
tm.pool = &pool;               // the exit needs the pool too, or it runs single-threaded
const csg::SoupMesh out = csg::to_mesh(s, tm);  // stitching, T-junctions, splitting, triangulation
```

A `PolySoup` carries **the generation-zero input meshes themselves** (`sources`) and
**an expression tree for the indicator function** (`indicator`). Classification
applies that indicator to each point's winding vector
$\mathbf{w} \in \mathbb{Z}^n$, so chains that cross the same surface twice — such as
$(A \cup B) \setminus B = A \setminus B$ — are expressible.
**Bit widths do not grow along a chain**: CSG creates no new planes, so a constructed
point is always the intersection of three input planes, however deep the chain.

## Design notes

### Planes and points are both homogeneous 4-vectors

A plane is `[a, b, c, d]` meaning **`N·x + d = 0`** (`d = -N·p₁`); a constructed point
is `[x, y, z, w]` (real coordinates `V/w`). Under this convention `side` becomes

```
sign(w) · sign(a·x + b·y + c·z + d·w)
```

a **single four-dimensional dot product**. Projective duality shows up directly in the
types, and one primitive disappears.

### Bit widths live in the type

```cpp
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;
```

Writing the predicate fixes the required limb count at compile time, and overflow is
prevented by the type system. The only place limb counts appear as numeric literals is
`include/krisite/geom/widths.hpp`.

### Parallelism does not change a single bit of the output

**Determinism is a requirement, not a by-product.** Each task writes only to the slot
at its own index and results are combined in index order. That makes the sequential
implementation (`threads = 1`) a complete oracle, so validation takes the strongest
possible form: **byte equality**, not "the volume and topology agree".

## Building

C++20; GCC 13+ / Clang 16+ / MSVC 2022+. Header-only, so adding `include/` is enough.

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKRISITE_CHECKED_ARITH=ON \
  -DKRISITE_BUILD_TESTS_WITH_GMP=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

### CMake options

| Option | Default | Meaning |
|---|---|---|
| `KRISITE_COORD_BITS` | 21 | Input coordinate bit width `b` |
| `KRISITE_CHECKED_ARITH` | ON in Debug | Overflow checking on every operation (independent of `NDEBUG`) |
| `KRISITE_BUILD_TESTS` | ON | Build the tests |
| `KRISITE_BUILD_TESTS_WITH_GMP` | **OFF** | GMP differential tests (LGPL, tests only) |
| `KRISITE_BUILD_TESTS_WITH_MANIFOLD` | **OFF** | Manifold oracle (Apache-2.0, tests only) |
| `KRISITE_BUILD_MUTANTS` | OFF | Mutation tests |
| `KRISITE_DEFAULT_ADAPTIVE` | OFF | Make adaptive subdivision the default |
| `KRISITE_DEFAULT_THREADS` | 1 | Default for `BoolOptions::threads` |
| `KRISITE_BUILD_BENCH` | OFF | Build the benchmarks |

### Platform validation happens in CI

The Linux (GCC/Clang) / macOS (Apple Silicon) / Windows (MSVC) × `b = 21, 26` matrix
is defined by [`.github/workflows/ci.yml`](.github/workflows/ci.yml), which also runs
the portable fallback paths, UBSan/ASan, GMP differential tests, the Manifold oracle,
mutation tests, ThreadSanitizer, the parallel configuration, the $n$-ary path and the
benchmarks. The operational policy lives in [`docs/ROADMAP.md`](docs/ROADMAP.md).

## Documentation

| File | Contents |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **Where the project stands. Read this first** |
| `docs/SPEC-phase<N>.md` | Per-phase specification (**the spec is authoritative**) |
| `docs/IMPL-phase<N>.md` | Per-phase implementation notes (decisions, and hypotheses that were ruled out) |
| [`docs/DESIGN-phase5-hotspots.md`](docs/DESIGN-phase5-hotspots.md) | Phase 5 performance deliberation log |
| [`docs/DECISION-core-contract.md`](docs/DECISION-core-contract.md) | How the core contract ($n$-ary, WNV, core/post-processing split) was decided |
| [`docs/BENCH.md`](docs/BENCH.md) | **Benchmarks and measurements** (numbers are authoritative here) |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | Third-party components and the machinery that enforces their constraints |
| [`docs/STYLE.md`](docs/STYLE.md) | Coding conventions |

Documentation is written in Japanese.

## Licence

**MIT**, and that constraint takes precedence over everything else.

The library itself (`include/krisite/`) has no external dependencies. GMP (LGPL) and
Manifold (Apache-2.0) are used **only as test oracles**, are disabled by default, and
are never part of a distribution. See
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

**No GPL/LGPL code (CGAL, Indirect_Predicates, OpenMeshCraft, VCGlib, …) has been
consulted, quoted, or ported.**

## Roadmap

| Phase | Contents | Status |
|---|---|---|
| 0 | Fixed-width exact integers + plane-based predicates | Complete (2026-08-26) |
| 1 | Minimal validation of output extraction (fixed-depth subdivision, single-threaded) | Complete (2026-08-27) — **verdict: continue** |
| 2 | Adaptive subdivision, constructed-point reuse, output semantics | Complete (2026-08-28) |
| 3 | Core redesign ($n$-ary, WNV, local BSP, convex split; single-threaded) | Complete (2026-08-29) |
| 4 | Parallelism (core + exit; **determinism required**) | Complete (2026-08-29) |
| **5** | **Thingi10K validation, performance targets** | **In progress** (breakdown below) |
| 6+ | Point-cloud codec, GWN, meshing | Not started |

### Phase 5 breakdown

**Correctness first, performance second** — the order matters.

| Stage | Scope | What it checks | Status |
|---|---|---|---|
| **CP1** | Solid, manifold, **non-self-intersecting** — 1,000 models → **500 pairs** | Correctness on real data; the common ground for comparison with EMBER and FARMA | **Complete** (499 succeeded, **1 failed**) — **the same configuration** as the 12 CP2 failures |
| **CP1.5** | — | **Removing the superlinearity**, a precondition for CP2/CP3 being runnable at all | **Complete** ($P \to t$ from 1.31 to 1.23; the CP2 estimate from 350 to **100 hours**) |
| **CP2** | Models that **do** self-intersect | Whether the operation itself resolves self-intersection into a clean output (self-union) | **84 pairs complete** (72 succeeded, **12 failed**). Every failure is a configuration where **contact splitting is not yet implemented** |
| **CP3** | No constraint on solidity, manifoldness or self-intersection | Whether the entry checks work on **non-PWN and degenerate models** | Not started |
| **CP4 onwards** | — | **Setting and pursuing performance targets** | Not started (**no target has been set yet**) |

Measurements are in [`docs/BENCH.md`](docs/BENCH.md); the reasoning behind each
decision is in [`docs/IMPL-phase5.md`](docs/IMPL-phase5.md).

## References

| Short | Reference |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **FARMA** | Cherchi, Livesu, Scateni, Attene. *Fast and Robust Mesh Arrangements using Floating-point Arithmetic.* ACM TOG 39(6), 2020. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

**The references were a starting point, not a source**: the EMBER and OEBSP papers
were not obtained until just before Phase 3. What was adopted after reading them, and
what was derived independently, is recorded in
[`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md).
