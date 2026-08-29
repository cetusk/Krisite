*[English version here](README.en.md)*

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/krisite-logo-dark.svg">
    <img src="assets/krisite-logo.svg"
         alt="Krisite — exact, plane-based geometry for point clouds and meshes"
         width="711">
  </picture>
</p>


# Krisite — exact, plane-based geometry for point clouds and meshes

`Krisite` は 3D データを扱う C++20 ヘッダオンリーライブラリです。
最終目標は点群圧縮・メッシュ化・厳密ブール演算の統合です。
**Phase 0（算術基盤）・Phase 1（出力抽出の最小検証）・Phase 2（適応分割と意味論の確定）・
Phase 3（中核の再設計）は完了しています。** Phase 1 の判断は「続行」でした。
現在は **Phase 4（work-stealing 並列化）**に向かっています。
現在地は [`docs/ROADMAP.md`](docs/ROADMAP.md)。

## いま提供しているもの

| 層 | 内容 | 仕様 |
|---|---|---|
| `arith/` | 固定幅の厳密整数演算（動的確保なし・例外なし・グローバル状態なし） | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | 平面ベース幾何述語。**幅を型で表す**ので上界超過がコンパイル時に防がれる | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | 厳密ブール演算（$\cup$ / $\cap$ / $\setminus$、**$n$ 項**）、位相検査、**適応分割 + early-out + 構成点の保持**、**局所 BSP**、**接触の分裂**、**WNV による分類**、入口の**凸分割** | [`SPEC-phase1.md`](docs/SPEC-phase1.md) 〜 [`SPEC-phase3.md`](docs/SPEC-phase3.md) |

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
#include <krisite/csg/boolean.hpp>       // 二項の boolean_op
#include <krisite/csg/soup_boolean.hpp>  // n 項の boolean（PolySoup 経路）
#include <krisite/csg/to_mesh.hpp>

using namespace krisite;

mesh::TriMesh A = /* 整数座標の閉じた向き付き三角メッシュ */;
mesh::TriMesh B = /* 同上 */, C = /* 同上 */;

csg::BoolStats st;
// depth は八分木の分割深度（実行時パラメータ）。意味論には影響しない
csg::BoolMesh r = csg::boolean_op(A, B, csg::BoolOp::Union, /*depth=*/2, &st);

// 出力の頂点は座標を持たない構成点（3 平面の交点）のまま
auto t = mesh::check_topology(r.triangles);
assert(t.ok());                  // 辺多様体・頂点多様体・向きの整合・退化なし

// 3 枚以上は PolySoup 経路で。**中間結果を丸めません**（型が CSG について閉じる）
csg::BoolOptions opt;
opt.depth = 2;
csg::PolySoup s = csg::boolean(csg::from_mesh(A), csg::from_mesh(B), csg::BoolOp::Union, opt);
s = csg::boolean(s, csg::from_mesh(C), csg::BoolOp::Difference, opt);   // (A ∪ B) \ C

assert(s.source_count() == 3);                // 入力 3 枚がそのまま残っている
assert(s.sources[0].vertices == A.vertices);  // 1 ビットも変わっていない
const csg::SoupMesh out = csg::to_mesh(s);    // 縫合・T 解決・分裂・三角形化はここだけ
```

`PolySoup` が持つのは**第 0 世代の入力メッシュそのもの**（`sources`）と、
**指示関数の式木**（`indicator`）です。分類は各点の巻き数ベクトル
$\mathbf{w} \in \mathbb{Z}^n$ に指示関数を適用して行うので、
$(A \cup B) \setminus B = A \setminus B$ のように**同じ曲面を 2 度跨ぐ連鎖**も書けます
（内外の 1 ビットでは表せません）。**ビット幅は連鎖で伸びません** —
CSG は新しい平面を作らないので、構成点は何段重ねても
「入力平面から 3 枚を選んだ交点」のままです（実測: 1 段 / 2 段 / 3 段とも 141 ビット）。

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
./build-rel/bench/pred_bench   # 述語のスループット
./build-rel/bench/soup_bench   # 入口・中核・出口の内訳（SPEC-phase3 §11）
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

**「検証されなければならない」と「全 PR で走る」は別です。** 仕様が要求するのは前者で、
実行頻度は運用の裁量です。方針は [`docs/ROADMAP.md`](docs/ROADMAP.md)「CI の方針」が正で、
段は 4 つあります。

| 段 | いつ | 中身 |
|---|---|---|
| PR ゲート | 毎 PR | 書式 + 主要 1 構成の全テスト |
| 影響トリガ | 該当パスが変わった PR | 変異テスト、プラットフォーム行列、Manifold |
| **main ゲート** | main への push | **全ジョブ** |
| **フェーズ完了** | 締めるとき | **全ジョブが green であることを、対象の SHA 付きで記録** |

**フェーズを閉じる前に全ジョブが green であることは緩めていません。** 段分けが変えるのは
頻度だけで、検査の中身は 1 件も減っていません。

CI は次のジョブを回します。

| ジョブ | 内容 |
|---|---|
| **PR ゲート** | 主要 1 構成（Linux Clang, b=21）の全テスト。**毎 PR 必ず走ります** |
| 変更の範囲 | 変更されたパスから、下のジョブを回すかを決める（影響トリガ） |
| 可搬経路 | `__int128` / 32bit 筆算のフォールバックを明示的に強制 |
| **UBSan / ASan** | **検査 OFF（出荷時構成）での未定義動作**を直接見る |
| GMP 差分テスト | 算術 $10^7$ 件・述語 $10^6$ 件の突き合わせ、体積の恒等式 |
| **適応分割モード** | 既定を適応分割に反転し、**Phase 1 の検査体系が全通過する**ことを確認 |
| Manifold 正解器 | ブール出力の連結成分数と種数を独立実装と照合 |
| **変異テスト** | 意図的に誤りを埋め、**検出されること**と**検出しない組合せ**の両方を固定 |
| **GMP 付き変異テスト** | 体積が意味を持つ変異だけを GMP 付きで回す（検出 2 件 + 素通り 4 件） |
| clang-format | 書式 |

## ドキュメント

| ファイル | 内容 |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **現在地とフェーズの全体像。まずここを読む** |
| [`docs/SPEC-phase3.md`](docs/SPEC-phase3.md) | **Phase 3 の仕様。** $n$ 項の契約、WNV、中核と後処理の分離 |
| [`docs/IMPL-phase3.md`](docs/IMPL-phase3.md) | Phase 3 の実装ノート（**完了時点**）。判断と根拠、**踏んだ誤りと訂正** |
| [`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) | Phase 3 の仕様を決めるまでの議論ログ |
| [`docs/DECISION-core-contract.md`](docs/DECISION-core-contract.md) | 中核の契約（$n$ 項・WNV・中核と後処理の分離）を決めた経緯 |
| [`docs/SPEC-phase2.md`](docs/SPEC-phase2.md) | Phase 2 の仕様。分割平面の絞り込み、適応分割、非多様体出力の意味論 |
| [`docs/IMPL-phase2.md`](docs/IMPL-phase2.md) | Phase 2 の実装ノート（**完了時点**）。判断と根拠、機構が検出器を動かした記録 |
| [`docs/SPEC-phase1.md`](docs/SPEC-phase1.md) | Phase 1 の仕様。縫合の可否判定、テストコーパス、中止条件 |
| [`docs/IMPL-phase1.md`](docs/IMPL-phase1.md) | Phase 1 の実装ノート（**完了時点**）。判断と根拠、**つまずいた点と訂正** |
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 の仕様。ビット幅解析、述語一覧、テスト要件 |
| [`docs/IMPL-phase0.md`](docs/IMPL-phase0.md) | Phase 0 の実装ノート。**なぜそう作ったか**、検証の設計と検出力、Phase 1 への申し送り |
| [`docs/BENCH.md`](docs/BENCH.md) | ベンチマークの基準線と、Phase 1 / 2 / 3 の計数 |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | 第三者コンポーネントの扱いと、それを機構で保証する仕組み |
| [`docs/STYLE.md`](docs/STYLE.md) | コーディング規約（命名、書式、算術コードの制約） |
| [`assets/BRAND.md`](assets/BRAND.md) | ロゴとテーマカラーの定義 |
| [`tools/README.md`](tools/README.md) | 文書を直す一度きりの改訂スクリプト（**ライブラリの一部ではありません**） |

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
| 1 | 出力抽出の最小検証（固定深度分割、単スレッド） | 完了（2026-08-27）**続行を決定** |
| **2** | **適応分割・構成点の保持・意味論の確定** | **完了（2026-08-28）** |
| **3** | **中核の再設計**（$n$ 項、WNV、局所 BSP、凸分割。単スレッド） | **完了（2026-08-29）** |
| **4** | **work-stealing 並列化** | **次** |
| 5 | Thingi10K 全件検証・性能目標 | 未着手 |
| 6+ | 点群コーデック、GWN、メッシュ化 | 未着手 |

**Phase 1 は分岐点でした。判断は「続行」です**（2026-08-27 に決定済み）。中止条件
（`SPEC-phase1.md` §11）はいずれにも該当せず、とくに最も危険だった「ビット幅の上界が
定まらない機構が必要になる」に最も余裕がありました。**固定幅整数という前提が最後まで
崩れなかったこと**が根拠です。

**Phase 2 は CP1〜CP5 をすべて通過して完了しました**（2026-08-28）。完了条件は正しさ
だけで、性能目標は置いていません。

**Phase 3 は中核を作り直しました。** ブール演算は連鎖して使うのが常態なので、
$n$ 項演算として一度に評価でき、**中間結果を丸めずに繋げる**ことを要件に加えました。

```
from_mesh : TriMesh(整数) → PolySoup      入口。量子化・凸分割・辺平面の構成
boolean   : PolySoup × … → PolySoup       ★ CSG について閉じる。n 項
to_mesh   : PolySoup → TriMesh            出口。縫合・T 解決・再併合・三角形化
```

分類は符号ベクトルから**一般化巻き数ベクトル（WNV）**に変わりました。
断片の分割も、全支持平面による過剰分割から**局所 BSP**に置き換わっています
（生の断片が 80.3%、正準化後が 78.6% に減りました）。入口には**凸分割**が入り、
実行時に三角形化と切り替えられます（片が 54.5%、平面が 64.9% に減ります）。詳細は
[`docs/ROADMAP.md`](docs/ROADMAP.md)。

### 実測値（Phase 1 → Phase 2）

全 22 ケース × 3 演算 × （固定深度 0〜3 + 適応分割）の実測です。

| 数値 | Phase 1 | Phase 2 | 意味 |
|---|---|---|---|
| 葉 | 32,256 | **7,434** | 適応分割。セル並列の粒度 |
| 断片（正準化後） | 15,624 | **8,217** | 適応分割 → early-out |
| `intersect3` の呼び出し | 1,539,819 | **28,836（1.9%）** | **構成点の保持** |
| `side` : `intersect3` の比 | 1.26 : 1 | **160.57 : 1** | **Phase 1 の最優先課題は解消しました。** 混合幅 `det3` の優先度は下がりました |
| 全コーパスの壁時計 | 220.9 ms | **47.6 ms** | 記録。**±15% 動くので 1 割の差に意味はありません** |
| 多様体でない出力の除外 | 3 構成 | **0 件** | 接触の分裂 |

Phase 1 が明らかにした構造的な数値も残ります。

| 数値 | 実測 | 意味 |
|---|---|---|
| 1 点に集まる平面の最大枚数 | 軸平行 **3** / 斜面あり **9** | **軸平行だけのコーパスでは 3 枚を超えられません** |
| 値ベース併合の発火率 | 最大 **44%** | 平面3つ組をキーにする第1段だけでは足りない |
| 併合グループの空間的な広がり | **1 セルとその隣接まで** | 大域整列は不要。並列化はセル単位で閉じる |

### 実測値（Phase 3）

中核と後処理を分けて測っています。**先行研究の数値がどこまでを含むかを確かめないと
比較になりません** — EMBER の 1.6 ms はスープを出すまでの時間です。

| 区分 | 時間 | 割合 |
|---|---:|---:|
| `from_mesh`（入口） | 0.8 ms | 0.7% |
| **`boolean`（中核）** | **75.9 ms** | **64.6%** |
| `to_mesh`（出口） | 40.9 ms | 34.8% |

| 機構 | 効果 |
|---|---|
| 局所 BSP（過剰分割の置き換え） | 生の断片 **80.3%** / 正準化後 **78.6%** |
| 凸分割（入口、実行時切替） | 片 **54.5%** / 平面 **64.9%** / 断片 **76.8%** |
| 連鎖 1 / 2 / 3 段のビット幅 | **141 のまま伸びない** |

詳細は [`docs/BENCH.md`](docs/BENCH.md)、判断の経緯は
[`docs/IMPL-phase1.md`](docs/IMPL-phase1.md)、[`docs/IMPL-phase2.md`](docs/IMPL-phase2.md)、
[`docs/IMPL-phase3.md`](docs/IMPL-phase3.md)。

## 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **FARMA** | Cherchi, Livesu, Scateni, Attene. *Fast and Robust Mesh Arrangements using Floating-point Arithmetic.* ACM TOG 39(6), 2020. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

**論文は記述だけを読んでいます。実装は 1 行も参照していません。**
設計には論文に述べられた考え方（局所 BSP、半開区間の割り当て、巻き数による分類など）を
取り入れていますが、ビット幅の導出・辺平面の構成・T 頂点の解決・分類の分解などは
自分で導いたものです。

**GPL / LGPL のコード（CGAL、Indirect_Predicates、OpenMeshCraft、VCGlib 等）は
参照・引用・移植していません。**
