*[English version here](README.en.md)*

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/krisite-logo-dark.svg">
    <img src="assets/krisite-logo.svg"
         alt="Krisite — exact, plane-based geometry for point clouds and meshes"
         width="711">
  </picture>
</p>


# Krisite: Exact, plane-based geometry for point clouds and meshes

`Krisite` は 3D データを扱う C++20 ヘッダオンリーライブラリです。
最終目標は点群圧縮・メッシュ化・厳密ブール演算の統合です。
**Phase 0（算術基盤）と Phase 1（出力抽出の最小検証）は完了し、Phase 1 の判断は「続行」です。**
現在は **Phase 2（適応分割と構成点の保持）** を進めています。現在地は [`docs/ROADMAP.md`](docs/ROADMAP.md)。

## いま提供しているもの

| 層 | 内容 | 仕様 |
|---|---|---|
| `arith/` | 固定幅の厳密整数演算（動的確保なし・例外なし・グローバル状態なし） | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | 平面ベース幾何述語。**幅を型で表す**ので上界超過がコンパイル時に防がれる | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | 厳密ブール演算（$\cup$ / $\cap$ / $\setminus$）、位相検査、固定深度の空間分割 | [`SPEC-phase1.md`](docs/SPEC-phase1.md) |

**浮動小数点も許容誤差も使いません。** 判定はすべて整数の厳密演算です。

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

ブール演算も同じ厳密性の上に乗っています。

```cpp
#include <krisite/csg/boolean.hpp>

using namespace krisite;

mesh::TriMesh A = /* 整数座標の閉じた向き付き三角メッシュ */;
mesh::TriMesh B = /* 同上 */;

csg::BoolStats st;
// depth は八分木の分割深度（実行時パラメータ）。意味論には影響しない
csg::BoolMesh r = csg::boolean_op(A, B, csg::BoolOp::Union, /*depth=*/2, &st);

// 出力の頂点は座標を持たない構成点（3 平面の交点）のまま
auto t = mesh::check_topology(r.triangles);
assert(t.ok());                  // 辺多様体・頂点多様体・向きの整合・退化なし
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
| `KRISITE_BUILD_TESTS_WITH_MANIFOLD` | **OFF** | Manifold 正解器（Apache-2.0、テスト専用） |
| `KRISITE_BUILD_MUTANTS` | OFF | 変異テスト（検査 OFF の構成が必要） |
| `KRISITE_DEFAULT_ADAPTIVE` | OFF | ブール演算の既定を適応分割 + early-out + 構成点の保持にする |
| `KRISITE_BUILD_BENCH` | OFF | ベンチマークのビルド |

`KRISITE_CHECKED_ARITH` は `NDEBUG` とは独立に効きます。Release ビルドでも
`-DKRISITE_CHECKED_ARITH=ON` を渡せば検査が入ります。

### プラットフォーム検証は CI で行います

Linux(GCC/Clang) / macOS(Apple Silicon) / Windows(MSVC) × `b = 21, 26` のマトリクスは
[`.github/workflows/ci.yml`](.github/workflows/ci.yml) が正です（SPEC §9）。
開発コンテナにこれらのツールチェーンを入れる必要はありません。

CI は他に次のジョブを回します。

| ジョブ | 内容 |
|---|---|
| 可搬経路 | `__int128` / 32bit 筆算のフォールバックを明示的に強制 |
| GMP 差分テスト | 算術 $10^7$ 件・述語 $10^6$ 件の突き合わせ、体積の恒等式 |
| **適応分割モード** | 既定を適応分割に反転し、**Phase 1 の検査体系が全通過する**ことを確認 |
| Manifold 正解器 | ブール出力の連結成分数と種数を独立実装と照合 |
| **変異テスト** | 意図的に誤りを埋め、**検出されること**と**検出しない組合せ**の両方を固定 |
| clang-format | 書式 |

## ドキュメント

| ファイル | 内容 |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **現在地とフェーズの全体像。まずここを読む** |
| [`docs/SPEC-phase2.md`](docs/SPEC-phase2.md) | Phase 2 の仕様。分割平面の絞り込み、適応分割、非多様体出力の意味論 |
| [`docs/IMPL-phase2.md`](docs/IMPL-phase2.md) | Phase 2 の実装ノート（**CP1 到達時点**）。判断と根拠、開発環境のつまずき |
| [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) | Phase 1 の仕様。縫合の可否判定、テストコーパス、中止条件 |
| [`docs/IMPL-phase1.md`](docs/IMPL-phase1.md) | Phase 1 の実装ノート（**CP3 到達時点**）。判断と根拠、**つまずいた点と訂正** |
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 の仕様。ビット幅解析、述語一覧、テスト要件 |
| [`docs/IMPL-phase0.md`](docs/IMPL-phase0.md) | Phase 0 の実装ノート。**なぜそう作ったか**、検証の設計と検出力、Phase 1 への申し送り |
| [`docs/BENCH.md`](docs/BENCH.md) | ベンチマークの基準線と、Phase 1 / Phase 2 の計数 |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | 第三者コンポーネントの扱いと、それを機構で保証する仕組み |
| [`docs/STYLE.md`](docs/STYLE.md) | コーディング規約（命名、書式、算術コードの制約） |
| [`assets/BRAND.md`](assets/BRAND.md) | ロゴとテーマカラーの定義 |

## ライセンス

**MIT。** この制約は他のすべてに優先します。

ライブラリ本体（`include/krisite/`）は外部依存を一切持ちません。GMP（LGPL）と
Manifold（Apache-2.0）は**テストの正解器としてのみ**使い、既定では無効です。
配布物には含まれません。

第三者コンポーネントの一覧と、それを機構で保証している仕組みは
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) にあります。

## ロードマップ

| Phase | 内容 | 状態 |
|---|---|---|
| 0 | 固定幅厳密整数 + 平面ベース述語 | 完了（2026-08-26） |
| **1** | **出力抽出の最小検証**（固定深度分割、単スレッド） | **完了条件を充足（2026-08-27）** |
| **2** | **適応分割・構成点の保持・意味論の確定** | **進行中（CP1 通過）** |
| 3 | work-stealing 並列化、継ぎ目の整合性検証 | 未着手 |
| 4 | Thingi10K 全件検証 | 未着手 |
| 5+ | 点群コーデック、GWN、メッシュ化 | 未着手 |

**Phase 1 は分岐点でした。判断は「続行」です。** 中止条件（`SPEC-phase1.md` §11）は
いずれにも該当せず、とくに最も危険だった「ビット幅の上界が定まらない機構が必要になる」に
最も余裕がありました。**固定幅整数という前提が最後まで崩れなかったこと**が根拠です。
詳細は [`docs/ROADMAP.md`](docs/ROADMAP.md)。

### Phase 1 で得られた数値

判断と、次のフェーズの設計に使うものです。全 18 ケース × 3 演算 × 深度 0〜3 の実測。

| 数値 | 実測 | 意味 |
|---|---|---|
| 1 点に集まる平面の最大枚数 | 軸平行 **3** / 斜面あり **9** | **軸平行だけのコーパスでは 3 枚を超えられません** |
| 値ベース併合の発火率 | 最大 **44%** | 平面3つ組をキーにする第1段だけでは足りない |
| 併合グループの空間的な広がり | **1 セルとその隣接まで** | 大域整列は不要。並列化はセル単位で閉じる |
| `side` : `intersect3` の呼び出し比 | **1.26 : 1** | 構成点を作り直しており、保持すれば大きく減る |

詳細は [`docs/BENCH.md`](docs/BENCH.md)、判断の経緯は
[`docs/IMPL-phase1.md`](docs/IMPL-phase1.md)。

## 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

実装はすべて論文からの再実装です。GPL/LGPL のコード（CGAL、Indirect_Predicates、
OpenMeshCraft、VCGlib 等）は参照・引用・移植していません。
