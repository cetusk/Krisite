*[English version here](README_en.md)*

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

**固定幅整数の厳密演算による、平面の側・向き・交差の判定です**（浮動小数点や許容誤差を使わない）。

現在は **Phase 5（Thingi10K での実データ検証・性能目標）**を進めています。
現在地は [`docs/ROADMAP.md`](docs/ROADMAP.md) が唯一の情報源です。

## いま提供しているもの

| 層 | 内容 | 仕様 |
|---|---|---|
| `arith/` | 固定幅の厳密整数演算（動的確保なし・例外なし・グローバル状態なし） | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `geom/` | 平面ベース幾何述語。**幅を型で表す**ので上界超過がコンパイル時に防がれる | [`SPEC-phase0.md`](docs/SPEC-phase0.md) |
| `mesh/` `octree/` `csg/` | 厳密ブール演算（$\cup$ / $\cap$ / $\setminus$、**$n$ 項**）、位相検査、適応分割、局所 BSP、**接触の分裂**（辺のまわりの角度順による組の決定を含む）、WNV による分類、**退化三角形を作らない三角形化**、入口の凸分割 | [`SPEC-phase1.md`](docs/SPEC-phase1.md) 〜 [`SPEC-phase3.md`](docs/SPEC-phase3.md) |
| `par/` | 永続スレッドプール。**出力はスレッド数に依らずビット単位で同一** | [`SPEC-phase4.md`](docs/SPEC-phase4.md) |

## 使い方

```cpp
#include <krisite/krisite.hpp>

using namespace krisite::geom;

IPoint a{0, 0, 0}, b{100, 0, 0}, c{0, 100, 0}, d{7, 9, 3};

// 三角形の支持平面（法線と d のビット幅は b から自動導出される）
PlaneD pl = plane_from_triangle(a, b, c);

int s = side(pl, d);          // -1 / 0 / +1、厳密
int o = orient3d(a, b, c, d); // 同上

// 3 平面の交点は座標を持たず、同次座標の構成点として表現される
HPointD v = intersect3(pl, other1, other2);
int s2 = side(pl, v);         // 除算なしで厳密に判定
```

ブール演算も同じ厳密性の上に乗っています。

```cpp
#include <krisite/csg/soup_boolean.hpp>  // n 項の boolean（PolySoup 経路）
#include <krisite/csg/to_mesh.hpp>
#include <krisite/par/thread_pool.hpp>

using namespace krisite;

mesh::TriMesh A = /* 整数座標の閉じた向き付き三角メッシュ */, B = /* 同上 */, C = /* 同上 */;

par::ThreadPool pool(8);
csg::BoolOptions opt;
opt.depth = 6;                 // 八分木の深度。実行時パラメータで、意味論には影響しない
opt.threads = 8;
opt.pool = &pool;              // **プールは持ち回してください**（生成は高い）

// 中間結果を丸めません（型が CSG について閉じる）
csg::PolySoup s = csg::boolean(csg::from_mesh(A), csg::from_mesh(B), csg::BoolOp::Union, opt);
s = csg::boolean(s, csg::from_mesh(C), csg::BoolOp::Difference, opt);   // (A ∪ B) \ C

assert(s.source_count() == 3);                // 入力 3 枚がそのまま残っている
assert(s.sources[0].vertices == A.vertices);  // 1 ビットも変わっていない

csg::ToMeshOptions tm;
tm.threads = 8;
tm.pool = &pool;               // 出口にも渡さないと 1 スレッドで回ります
const csg::SoupMesh out = csg::to_mesh(s, tm);  // 縫合・T 解決・分裂・三角形化はここだけ
```

`PolySoup` が持つのは**第 0 世代の入力メッシュそのもの**（`sources`）と、
**指示関数の式木**（`indicator`）です。分類は各点の巻き数ベクトル
$\mathbf{w} \in \mathbb{Z}^n$ に指示関数を適用して行うので、
$(A \cup B) \setminus B = A \setminus B$ のように**同じ曲面を 2 度跨ぐ連鎖**も書けます。
**ビット幅は連鎖で伸びません** — CSG は新しい平面を作らないので、
構成点は何段重ねても「入力平面から 3 枚を選んだ交点」のままです。

## 設計の要点

### 平面も点も同次 4 元ベクトル

平面は `[a, b, c, d]` として **`N·x + d = 0`**（`d = -N·p₁`）、
構成点は `[x, y, z, w]`（実座標は `V/w`）。この流儀により `side` は

```
sign(w) · sign(a·x + b·y + c·z + d·w)
```

という**単一の 4 次元内積**になります。射影双対性が型にそのまま現れ、プリミティブが一つ減ります。

### 型でビット幅を表現する

```cpp
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;
```

述語の式を書いた時点で必要なリム数がコンパイル時に確定し、オーバーフローが型で防がれます。
リム数を数値リテラルで書く場所は `include/krisite/geom/widths.hpp` の 1 箇所だけです。

### 並列化しても出力は 1 ビットも変わりません

**決定性は要求であって副産物ではありません。** 各タスクは自分の添字のスロットにだけ書き、
結合は添字順に行います。**逐次実装（`threads = 1`）がそのまま完全な正解器**になるので、
検証が「体積と位相が一致」ではなく**「バイト列が一致」**という最も強い形になります。

## ビルド

C++20。GCC 13+ / Clang 16+ / MSVC 2022+。ヘッダオンリーなので `include/` を通すだけでも使えます。

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKRISITE_CHECKED_ARITH=ON \
  -DKRISITE_BUILD_TESTS_WITH_GMP=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

### CMake オプション

| オプション | 既定 | 意味 |
|---|---|---|
| `KRISITE_COORD_BITS` | 21 | 入力座標のビット幅 `b` |
| `KRISITE_CHECKED_ARITH` | Debug で ON | 全演算のオーバーフロー検査（`NDEBUG` とは独立に効きます） |
| `KRISITE_BUILD_TESTS` | ON | テストのビルド |
| `KRISITE_BUILD_TESTS_WITH_GMP` | **OFF** | GMP 差分テスト（LGPL、テスト専用） |
| `KRISITE_BUILD_TESTS_WITH_MANIFOLD` | **OFF** | Manifold 正解器（Apache-2.0、テスト専用） |
| `KRISITE_BUILD_MUTANTS` | OFF | 変異テスト |
| `KRISITE_DEFAULT_ADAPTIVE` | OFF | ブール演算の既定を適応分割にする |
| `KRISITE_DEFAULT_THREADS` | 1 | `BoolOptions::threads` の既定 |
| `KRISITE_BUILD_BENCH` | OFF | ベンチマークのビルド |

### プラットフォーム検証は CI で行います

Linux(GCC/Clang) / macOS(Apple Silicon) / Windows(MSVC) × `b = 21, 26` のマトリクスは
[`.github/workflows/ci.yml`](.github/workflows/ci.yml) が正です。
可搬経路・UBSan/ASan・GMP 差分・Manifold 正解器・変異テスト・ThreadSanitizer・
並列モード・$n$ 項経路・ベンチマークを回します。運用の方針は
[`docs/ROADMAP.md`](docs/ROADMAP.md)「CI の方針」が正です。

## ドキュメント

| ファイル | 内容 |
|---|---|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | **現在地とフェーズの全体像。まずここを読む** |
| [`docs/CONTRACTS.md`](docs/CONTRACTS.md) | **確定した契約の索引。** 5 つの仕様書を辿らずに済む 1 枚 |
| [`docs/HANDOVER.md`](docs/HANDOVER.md) | **引き継ぎ用の現在地。** 記録の読み方と、判断待ちの一覧 |
| `docs/SPEC-phase<N>.md` | 各フェーズの仕様（**仕様が正**） |
| `docs/IMPL-phase<N>.md` | 各フェーズの実装ノート（判断と根拠、外した仮説の記録） |
| [`docs/IMPL-v2.md`](docs/IMPL-v2.md) | **2 人目の実装担当の記録**（2026-09-07 以降） |
| [`docs/DESIGN-phase5-hotspots.md`](docs/DESIGN-phase5-hotspots.md) | Phase 5 の性能の検討ログ（手法の比較） |
| [`docs/RESEARCH-perf.md`](docs/RESEARCH-perf.md) | 文献の性能主張と、その前提が Krisite に当てはまるかの検証 |
| [`docs/PERF.md`](docs/PERF.md) | 費用モデルと、測定の作法 |
| [`docs/DECISION-core-contract.md`](docs/DECISION-core-contract.md) | 中核の契約（$n$ 項・WNV・中核と後処理の分離）を決めた経緯 |
| [`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) | Phase 3 の仕様を決めるまでの議論ログ |
| [`docs/LOG-phase5-checkpoints.md`](docs/LOG-phase5-checkpoints.md) | Phase 5 のチェックポイントの記録 |
| [`docs/BENCH.md`](docs/BENCH.md) | **ベンチマークと実測値**（数値はここが正） |
| [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) | 第三者コンポーネントの扱いと、それを機構で保証する仕組み |
| [`docs/STYLE.md`](docs/STYLE.md) | コーディング規約 |
| [`assets/BRAND.md`](assets/BRAND.md) | ロゴとテーマカラーの定義 |

## ライセンス

**MIT。** この制約は他のすべてに優先します。

ライブラリ本体（`include/krisite/`）は外部依存を一切持ちません。GMP（LGPL）と
Manifold（Apache-2.0）は**テストの正解器としてのみ**使い、既定では無効です。
配布物には含まれません。詳細は
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)。

**GPL / LGPL のコード（CGAL、Indirect_Predicates、OpenMeshCraft、VCGlib 等）は
参照・引用・移植していません。**

## ロードマップ

| Phase | 内容 | 状態 |
|---|---|---|
| 0 | 固定幅厳密整数 + 平面ベース述語 | 完了（2026-08-26） |
| 1 | 出力抽出の最小検証（固定深度分割、単スレッド） | 完了（2026-08-27）**続行を決定** |
| 2 | 適応分割・構成点の保持・意味論の確定 | 完了（2026-08-28） |
| 3 | 中核の再設計（$n$ 項、WNV、局所 BSP、凸分割。単スレッド） | 完了（2026-08-29） |
| 4 | 並列化（中核 + 出口。**決定性を要求**） | 完了（2026-08-29） |
| **5** | **実データ検証（Thingi10K）・性能目標** | **進行中**（内訳は下記） |
| 6+ | 点群コーデック、GWN、メッシュ化 | 未着手 |

### Phase 5 の内訳（2026-09-15）

**正しさを先に、性能を後に**という順序です。
**CP2・CP3 の層化標本の実行と、追加の 5 対の GMP 診断を終え、保存数値の整合を照合しました。**
式 2 の不一致の原因を独立に検証する方法は、今後の設計課題です。
全母集団の検証や Phase 5 の完了を意味しません。

| 段 | 対象 | 何を確かめるか | 状態 |
|---|---|---|---|
| **CP1** | 詰まっていて多様体、**自己交差なし**の 1,000 模型 → **500 対** | 実データ検証の初期基準 | 初期実行は成功 499 / 失敗 1（過去の版の記録） |
| **CP1.5** | — | CP2 / CP3 を実行可能にするための性能改善 | 改善を実施。過去の外挿は現行の投入予算ではありません |
| **CP2** | 自己交差を含む入力（CP2 条件） | 実データでのブール演算と出力検査 | **層化標本 295 対を完走、駆動判定で全件成功** |
| **CP3** | solid・多様体・自己交差を問わない模型（CP3 条件） | 入口の検査と実データでの出力検査 | **層化標本 295 対を完走、駆動判定で全件成功** |
| **GMP 診断** | CP2 の 2 対 + CP3 の 3 対 | 4 出力の厳密体積整合性 | **5 対で式 1 の残差 0 を記録**。独立した正解証明ではありません |
| **CP4 以降** | — | **性能目標の設定と最適化** | 未着手（**目標はまだ置いていません**） |

### 接触の修復と標本検証

曲面が辺で接する配置への角度順の組分けに加え、頂点まわりの接触を扱う修復を実装しました。
2026-09-14〜15 の CP2 / CP3 層化標本は、**各 295 対、合計 590 実行**で
駆動判定の失敗・停止が 0 件でした。この標本では修復が **20 実行・132 辺**で発火し、
修復後の `unresolved` は 0 でした。これは標本の観測であり、全入力への保証ではありません。

### GMP 診断で確認した範囲

追加の 5 対では、合成検定を通した検査器で、4 演算
（$A\cup B$ / $A\cap B$ / $A\setminus B$ / $B\setminus A$）の符号付き体積を厳密に計算し、
**式 1: $V_U = V_{A\setminus B} + V_{B\setminus A} + V_I$** の残差 0 を記録しました。
4 演算の位相検査も通過し、`unresolved` は 0。旧 3 演算の互換ハッシュは、
対応する CP2 / CP3 の記録と **5 / 5 一致**しました。
保存された有理数から両式の残差を独立に再計算し、記録との一致を確認しています。
省略のない結果と指紋は [`docs/evidence/gmp_diag_r1/`](docs/evidence/gmp_diag_r1/README.md) にあります。

**式 1 は独立した正解器ではありません。** 共通の欠落や誤差の相殺は除外できません。
入力体積を使う式 2（$V_U+V_I=V_A+V_B$）は 5 対とも不一致ですが、
領域の体積として使うための巻き数 0/1 の前提は未確認であり、原因は未確定です。

**浮動小数点の篩から得た 261 実行は診断候補であって、厳密検査済み件数ではありません。**
全候補への拡大投入は未承認です。現在地と検証の限定は
[`docs/ROADMAP.md`](docs/ROADMAP.md) と [`docs/SPEC-phase5.md`](docs/SPEC-phase5.md) §5.10.14.74、
実装側の訂正記録は [`docs/DESIGN-phase5-vertex-level.md`](docs/DESIGN-phase5-vertex-level.md) §32 を参照してください。

### 性能の改善（2026-09-07 〜 09-08）

以下は当時の構成での記録です。現行版の所要時間や投入予算を保証するものではありません。

**すべて「出力が変わらないこと」を検査で守っています。**
**効果は演算回数で測っています**（時間は実行ごとに ±15% 動くため）。

| 機構 | 効果 | 出力 |
|---|---|---|
| 分裂の検証を同値な増分計算に置き換え | 出口 **1.8〜2.4 倍**、検証だけなら 11.9〜36.4 倍 | バイト一致 |
| T 解決の照合を多角形あたり 1 回 + 区間の二分探索に | `side` の評価 **7.3〜26 倍**減 | バイト一致 |
| 出力の外接箱を「元の多角形の箱 $\cap$ セル箱」に狭める | 葉の数 **−19〜−50%**、隅のレイキャスト **−74〜−98.6%** | 単発はバイト一致（連鎖は位相と体積で一致） |
| 到達可能性による early-out（EMBER §4.5.2） | 作った断片 **−0.5〜−40.8%**（$\setminus$ と $\cap$ で効き、$\cup$ では効きません） | バイト一致 |
| **代表点の float ヒントの欠陥を修正** | **主経路の成功率 2.5〜25.9% → 70.5〜86.2%** | バイト一致 |

**最後の 1 件は Phase 3 から 3 フェーズにわたって壊れていた欠陥です。**
**固定幅整数を double に落とす補助関数が、【すべての負の値】で $2^{128}$ を
返していました。**

> **厳密性は無傷でした。** ヒントは候補を出すだけで、判定は厳密演算です。
> **壊れていても答えは正しく、遅くなるだけでした。**
> **だから正しさの検査では見つかりません。**
> 経緯は [`docs/DESIGN-phase5-hotspots.md`](docs/DESIGN-phase5-hotspots.md) §20。

測定値は [`docs/BENCH.md`](docs/BENCH.md)、経緯と判断は
[`docs/IMPL-phase5.md`](docs/IMPL-phase5.md) と
[`docs/IMPL-v2.md`](docs/IMPL-v2.md) にあります。
**現在地は [`docs/ROADMAP.md`](docs/ROADMAP.md) が正です。**

## 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **FARMA** | Cherchi, Livesu, Scateni, Attene. *Fast and Robust Mesh Arrangements using Floating-point Arithmetic.* ACM TOG 39(6), 2020. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

**参考文献は設計の出発点であり、EMBER / OEBSP の本文は Phase 3 の直前まで入手して
いません。** 精読したあとに取り入れた考え方と、自分で導いたものとの差分は
[`docs/LOG-phase3-design.md`](docs/LOG-phase3-design.md) に記録しています。
