<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/krisite-logo-dark.svg">
    <img src="assets/krisite-logo.svg"
         alt="Krisite — exact, plane-based geometry for point clouds and meshes"
         width="711">
  </picture>
</p>

---

`Krisite` は 3D データを扱う C++20 ヘッダオンリーライブラリです。
最終目標は点群圧縮・メッシュ化・厳密ブール演算の統合ですが、
**現在は Phase 0（算術基盤）のみ実装されています。**

名前の由来は **kris**（短剣、かつ crystal の語幹 κρύσταλλος）+ **-ite**（鉱物の接尾辞）。
平面ベース表現を採る本ライブラリにとって、結晶が平面で囲まれた立体であることは
単なる比喩ではありません。

## Phase 0 で提供するもの

固定幅の厳密整数演算と、それに基づく平面ベース幾何述語。仕様は
[`docs/SPEC-phase0.md`](docs/SPEC-phase0.md)。

```cpp
#include <krisite/krisite.hpp>

using namespace krisite::geom;

IPoint a{0, 0, 0}, b{100, 0, 0}, c{0, 100, 0}, d{7, 9, 3};

// 三角形の支持平面（法線と d のビット幅は b から自動導出される）
PlaneD pl = plane_from_triangle(a, b, c);

int s  = side(pl, d);          // -1 / 0 / +1、厳密
int o  = orient3d(a, b, c, d); // 同上

// 3 平面の交点は座標を持たず、同次座標の構成点として表現される
HPointD v = intersect3(pl, other1, other2);
int s2 = side(pl, v);          // 除算なしで厳密に判定
bool lt = lex_less(v, other_v);
```

浮動小数点は一切使いません。すべての述語は固定幅整数の符号だけで決まり、
動的メモリ確保・例外・グローバル状態がないため、そのまま並列化できます。

### 平面も点も同次 4 元ベクトル

平面は `[a, b, c, d]` として **`N·x + d = 0`**（`d = -N·p₁`）の形で保持します。
構成点は `[x, y, z, w]`（実座標は `V/w`）。この流儀により `side` は

```
sign(w) · sign(a·x + b·y + c·z + d·w)
```

という **単一の 4 次元内積** になります。射影双対性が型にそのまま現れ、
プリミティブが一つ減ります。詳しくは `docs/SPEC-phase0.md` §3.1。

### 型でビット幅を表現する

本ライブラリの設計の要です。積は幅が増えることが型に現れます。

```cpp
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;
```

述語の式を書いた時点で必要なリム数がコンパイル時に確定し、オーバーフローが型で防がれます。
リム数を数値リテラルで書く場所は `include/krisite/geom/widths.hpp` の 1 箇所だけです。

## ビルド

C++20。GCC 13+ / Clang 16+ / MSVC 2022+。ヘッダオンリーなので `include/` を通すだけでも使えます。

```bash
# 開発時の既定
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKRISITE_CHECKED_ARITH=ON \
  -DKRISITE_BUILD_TESTS_WITH_GMP=ON
cmake --build build
ctest --test-dir build --output-on-failure

# 単一テスト
ctest --test-dir build -R fixed_int --output-on-failure

# ベンチマーク
cmake -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DKRISITE_CHECKED_ARITH=OFF -DKRISITE_BUILD_BENCH=ON
cmake --build build-rel
./build-rel/bench/pred_bench
```

### CMake オプション

| オプション | 既定 | 意味 |
|---|---|---|
| `KRISITE_COORD_BITS` | 21 | 入力座標のビット幅 `b` |
| `KRISITE_CHECKED_ARITH` | Debug で ON | 全演算のオーバーフロー検査 |
| `KRISITE_BUILD_TESTS` | ON | テストのビルド |
| `KRISITE_BUILD_TESTS_WITH_GMP` | **OFF** | GMP 差分テスト（LGPL、テスト専用） |
| `KRISITE_BUILD_BENCH` | OFF | ベンチマークのビルド |

`KRISITE_CHECKED_ARITH` は `NDEBUG` とは独立に効きます。Release ビルドでも
`-DKRISITE_CHECKED_ARITH=ON` を渡せば検査が入ります。

### プラットフォーム検証は CI で行います

Linux(GCC/Clang) / macOS(Apple Silicon) / Windows(MSVC) のマトリクスは
[`.github/workflows/ci.yml`](.github/workflows/ci.yml) が正です（SPEC §9）。
開発コンテナにこれらのツールチェーンを入れる必要はありません。
CI は可搬経路（`__int128` / 32bit 筆算）と GMP 差分テストも別ジョブで回します。

## ライセンス

**MIT。** この制約は他のすべてに優先します。

ライブラリ本体（`include/krisite/`）は外部依存を一切持ちません。GMP（LGPL）は
`KRISITE_BUILD_TESTS_WITH_GMP=ON` のときにテストバイナリからのみリンクされ、
配布物には含まれません。方針の詳細は [`CLAUDE.md`](CLAUDE.md) を参照してください。

## ロードマップ

| Phase | 内容 | 状態 |
|---|---|---|
| **0** | 固定幅厳密整数 + 平面ベース述語 | **実装済み** |
| 1 | 出力抽出の最小検証（立方体2個、固定深度分割、単スレッド） | 未着手 |
| 2 | 適応的再帰分割 + early-out 判定 | 未着手 |
| 3 | work-stealing 並列化、継ぎ目の整合性検証 | 未着手 |
| 4 | Thingi10K 全件検証 | 未着手 |
| 5+ | 点群コーデック、GWN、メッシュ化 | 未着手 |

## 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

実装はすべて論文からの再実装です。GPL/LGPL のコード（CGAL、Indirect_Predicates、
OpenMeshCraft、VCGlib 等）は参照・引用・移植していません。
