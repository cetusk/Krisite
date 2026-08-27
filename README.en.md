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
**Phase 0 (the arithmetic foundation) and Phase 1 (minimal validation of output
extraction) have both met their completion criteria.**
[`docs/ROADMAP.md`](docs/ROADMAP.md) is the single source of truth for where the
project stands (Japanese).

## What exists today

| Layer | Contents | Spec |
|---|---|---|
| `arith/` | Fixed-width exact integers — no dynamic allocation, no exceptions, no global state | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | Plane-based geometric predicates. **Widths live in the type**, so exceeding a derived bound is a compile error | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | Exact booleans ($\cup$ / $\cap$ / $\setminus$), topology checking, fixed-depth spatial subdivision | [`SPEC-phase1.md`](docs/SPEC-phase1.md) |

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
#include <krisite/csg/boolean.hpp>

using namespace krisite;

mesh::TriMesh A = /* closed, oriented triangle mesh on integer coordinates */;
mesh::TriMesh B = /* likewise */;

csg::BoolStats st;
// `depth` is the octree subdivision depth — a runtime parameter that
// must not affect the semantics.
csg::BoolMesh r = csg::boolean_op(A, B, csg::BoolOp::Union, /*depth=*/2, &st);

// Output vertices remain constructed points (intersections of three planes).
auto t = mesh::check_topology(r.triangles);
assert(t.ok());  // edge-manifold, vertex-manifold, consistently oriented, non-degenerate
```

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
./build-rel/bench/pred_bench
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
| `KRISITE_BUILD_BENCH` | OFF | Build the benchmarks |

`KRISITE_CHECKED_ARITH` is independent of `NDEBUG`. Passing
`-DKRISITE_CHECKED_ARITH=ON` enables the checks even in a Release build.

### Platform validation happens in CI

The matrix of Linux (GCC/Clang) / macOS (Apple Silicon) / Windows (MSVC) ×
`b = 21, 26` is defined in [`.github/workflows/ci.yml`](.github/workflows/ci.yml),
and that file is the authority. There is no need to install those toolchains in a
development container.

CI also runs:

| Job | Contents |
|---|---|
| Fallback paths | Forces the `__int128` and 32-bit schoolbook paths explicitly |
| GMP differential | $10^7$ arithmetic and $10^6$ predicate comparisons, plus volume identities |
| Manifold oracle | Compares component count and genus of boolean output against an independent implementation |
| **Mutation tests** | Injects deliberate faults and pins down **both** which tests catch them **and** which combinations do not |
| clang-format | Formatting |

## Documentation

The design documents are written in Japanese.

| File | Contents |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **Where the project stands. Start here** |
| [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) | Phase 1 spec: the stitching question, test corpus, abort conditions |
| [`docs/IMPL-phase1.md`](docs/IMPL-phase1.md) | Phase 1 implementation notes. Decisions, rationale, **and the mistakes that were corrected** |
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 spec: bit-width analysis, predicates, test requirements |
| [`docs/IMPL-phase0.md`](docs/IMPL-phase0.md) | Phase 0 implementation notes. **Why it is built this way**, and how the tests were designed to have detection power |
| [`docs/BENCH.md`](docs/BENCH.md) | Benchmark baseline and the Phase 1 measurements (per case, per depth) |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | Third-party components and the mechanisms that keep them out of the distributable |
| [`docs/STYLE.md`](docs/STYLE.md) | Coding conventions |
| [`assets/BRAND.md`](assets/BRAND.md) | Logo and theme colours |

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
| **1** | **Minimal validation of output extraction** (fixed-depth subdivision, single-threaded) | **Criteria met (2026-08-27)** |
| 2 | Adaptive recursive subdivision + early-out | Not started |
| 3 | Work-stealing parallelism, seam consistency | Not started |
| 4 | Thingi10K full-corpus validation | Not started |
| 5+ | Point-cloud codec, GWN, meshing | Not started |

**Phase 1 is the decision point.** Whether to continue with this approach is
decided here. If it is abandoned, `csg/` becomes a wrapper around Manifold and the
project skips Phases 2–4 to go straight to the Phase 5 line. See
[`docs/ROADMAP.md`](docs/ROADMAP.md) and [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) §11.

### What Phase 1 measured

These numbers drive the decision and the design of the next phases. Measured across
18 cases × 3 operations × depths 0–3.

| Quantity | Measured | Why it matters |
|---|---|---|
| Max planes meeting at one point | **3** axis-aligned / **9** with slanted faces | **An axis-aligned-only corpus structurally cannot exceed 3** |
| Rate of value-based vertex merging | up to **44%** | Keying on the plane triple alone is not enough |
| Spatial extent of a merge group | **one cell and its neighbours** | No global sort needed; parallelism closes at cell scope |
| `side` : `intersect3` call ratio | **1.26 : 1** | Constructed points are recomputed; caching them would cut this sharply |

Details in [`docs/BENCH.md`](docs/BENCH.md); the reasoning behind them in
[`docs/IMPL-phase1.md`](docs/IMPL-phase1.md).

## References

| Short | Reference |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

Everything is a reimplementation from the papers. No GPL/LGPL code (CGAL,
Indirect_Predicates, OpenMeshCraft, VCGlib and the like) has been consulted,
quoted, or ported.
