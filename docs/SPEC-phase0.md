# Krisite — Phase 0 仕様書：固定幅厳密整数演算と平面ベース述語

> **命名について**
> `Krisite` = **kris**（インドネシア・ジャワ由来の短剣。crystal の語幹 κρύσταλλος でもある）
> \+ **-ite**（鉱物名の接尾辞）。「綺麗に削り取る」と「結晶」を一語で表します。
> 平面ベース表現を採る本ライブラリにとって、結晶が平面で囲まれた立体であることは
> 単なる比喩ではありません。
>
> - 名前空間: `krisite::`（短縮別名 `namespace kri = krisite;` を提供してよい）
> - ファイル拡張子: `.kris`
> - マジックバイト: ASCII `"KRIS"` + `uint32` バージョン（Phase 5 以降で使用）
> - CMake オプション接頭辞: `KRISITE_`

---

## 0. この文書の位置づけ

本書は **Phase 0 のみ** を対象とします。Phase 1 以降は `docs/ROADMAP.md` を参照。

Phase 0 の成果物は、以降のすべてのモジュールが依存する **算術基盤** です。ここが不正確だと
上流のすべてが静かに壊れるため、**性能より正しさを優先**します。最適化は正しさの検証が
自動化された後に行ってください。

---

## 1. 背景（なぜこれを作るのか）

最終目標は、点群圧縮・メッシュ化・ブール演算を統合した高速 3D ツールキットです。
ブール演算の中核として **EMBER**（Trettner, Nehring-Wirxel, Kobbelt, SIGGRAPH 2022）の
再現実装を目指します。

EMBER は、頂点を明示的な座標ではなく **3枚の平面の交点** として表現し、すべての幾何述語を
**固定幅の多倍長整数**で厳密に評価します。これにより、

- 動的メモリ確保が発生しない（並列化と相性が良い）
- 浮動小数点フィルタとそのフォールバック経路が不要
- 退化構成（共平面・共線・重複頂点）でも分岐が増えない

という性質が得られ、既存の厳密手法に対して 100 倍規模の高速化が報告されています。

Phase 0 では、この「固定幅厳密整数」と「平面ベース述語」だけを、他から独立した形で作ります。

### 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* Computer-Aided Design 135, 2021. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

**Phase 0 の実装に EMBER 論文本体は必須ではありません。** 本書に必要な数式はすべて
展開済みです。ただし述語の一覧を確定する際（§7）には論文を参照してください。

---

## 2. 座標モデル

入力点はすべて **符号付き整数** です。浮動小数点は境界（ファイル入出力）でのみ扱います。

```
実座標 = origin + scale * integer_coord
origin : double[3]   バウンディングボックス原点
scale  : double      格子幅
coord  : int32[3]    格子座標、各軸 coord ∈ [-2^(b-1), 2^(b-1))
```

**座標は符号付き b ビット**です。符号なし b ビットの Morton 格子 `[0, 2^b)` から
バイアス `2^(b-1)` を引いて中心化したもの、と考えてください。§3.1 の「座標差 = b+1 ビット」
はこの前提のときにのみ成り立ちます。

> **正誤（2026-08 修正）**: 旧版は `|coord| < 2^b` と書いていましたが、これだと座標差が
> b+2 ビットになり §3.1 の表と矛盾します。上記が正です。

**Phase 0 では `b = 21` を既定とします。** 理由は 64bit Morton コードに 3 軸 21bit が
ちょうど収まるため（21 × 3 = 63）。ただしビット幅は **コンパイル時パラメータ** とし、
`b` を変えれば必要な整数幅が自動的に決まる設計にしてください（§3.3）。

---

## 3. ビット幅解析

### 3.1 導出

入力座標は **符号付き $b$ ビット**とする（§2。$-2^{b-1} \le \mathrm{coord} < 2^{b-1}$）。
以下、すべて符号付きのビット数（符号ビット込み）。

**平面** — 三角形 $(p_1, p_2, p_3)$ の支持平面。**同次4元ベクトル $[a,b,c,d]$ として
$N \cdot x + d = 0$ の形で保持します**（$d = -N \cdot p_1$）。

この流儀を採るのは、平面と点がどちらも同次4元ベクトルになり、`side` が単一の
4次元内積で書けるためです。射影双対性が型にそのまま現れ、プリミティブが一つ減ります。

| 量 | 式 | ビット幅 |
|---|---|---|
| 座標差 | $p_2 - p_1$ | $b + 1$ |
| 法線成分 | $N = (p_2-p_1) \times (p_3-p_1)$ | $2b + 3$ |
| オフセット | $d = -N \cdot p_1$ | $3b + 5$ |

**構成点** — 3 平面の交点、同次座標 $V = [x, y, z, w]$（実座標は $V/w$）：

| 量 | 式 | ビット幅 |
|---|---|---|
| $w$ | $\det(N_1, N_2, N_3)$ | $6b + 12$ |
| $x, y, z$ | Cramer（1 列を $-d$ ベクトルに置換） | $7b + 14$ |

> **置換列の符号に注意**: 平面を $N \cdot x + d = 0$ と置いたので、3 平面の交点が
> 満たす連立は $N \cdot X = -d$ です。したがって Cramer で置き換える列は $d$ ではなく
> $-d$ になります。$d$ をそのまま置換すると、符号が反転した点 $-V$ が返ります。
>
> 等価な流儀として $w = -\det(N_1, N_2, N_3)$ と定義して置換列を $d$ のままにする手も
> ありますが、本書は $w = \det(N_1, N_2, N_3)$ を採ります。

**主要述語** — 構成点 $V$ が平面 $(N, d)$ のどちら側か：

$$\mathrm{side}(N, d, V) = \mathrm{sign}(w) \cdot \mathrm{sign}(N \cdot V + d\,w)$$

| 量 | ビット幅 |
|---|---|
| $N \cdot V$ | $(2b+3) + (7b+14) + 2 = 9b + 19$ |
| $d \cdot w$ | $(3b+5) + (6b+12) = 9b + 17$ |
| **総和** | $\mathbf{9b + 20}$ |

**同次点の比較述語** — `cmp_h` の被符号値 $w_2 x_1 - w_1 x_2$：

| 量 | ビット幅 |
|---|---|
| $w_2 x_1$ | $(6b+12) + (7b+14) = 13b + 26$ |
| 差分 | $\mathbf{13b + 27}$ |

**これが Phase 0 で最大の述語です。** `side` ではありません。

### 3.2 検算

$b = 26$ → $9(26) + 20 = 254$ ビット。

平面ベースの厳密ブールで経験的に言われる「入力 26bit・演算 256bit」という規模と
一致します。独立に導いた上界が既知の経験値と合うので、本解析は妥当と判断します。

### 3.3 必要幅の表

| $b$ | `side` $9b{+}20$ | `cmp_h` $13b{+}27$ | 最大リム数（64bit） |
|---|---|---|---|
| 19 | 191 | 274 | 5 |
| **21** | **209** | **300** | **5** |
| 24 | 236 | 339 | 6 |
| 26 | 254 | 365 | 6 |
| 27 | 263 | 378 | 6 |

**$b = 21$ で必要なのは 320bit（5 リム）です。**

> **正誤（2026-08 修正）**: 旧版は「4 リムに 47bit の余裕」としていましたが、これは
> `side` のみの解析でした。`cmp_h` が $13b+27$ を要するため 5 リムが必要です。
> 型で幅を伝播する設計（§5.2）のおかげで、この誤りは実装に波及していません。

### 3.4 実装上の要求

- 各リム数を `constexpr` から自動導出すること。ハードコードした `256` を散らばらせない
- **デバッグビルドではすべての演算でオーバーフロー検査を有効にし、検出したら `assert` で停止**
- リリースビルドでは検査を外してよいが、`KRISITE_CHECKED_ARITH=ON` で有効化できること

> **Phase 0 の述語については解析済み**（`cmp_h` が最大で $13b+27$）。
> ただし **EMBER の radial sort と分類段の述語は未解析**です。Phase 2 以降で述語を
> 追加する際は、必ず同じ手順で導出して本表を更新してください。
>
> **乱択の実測値でリム数を減らさないこと。** 乱択は真の最大に届きません。
> 減らす判断は、再導出した上界が下のリム境界を下回るときだけ行ってください。

---

## 4. モジュール構成（Phase 0 の範囲）

```
include/krisite/
  arith/
    fixed_int.hpp      固定幅符号付き整数（テンプレート、リム数パラメータ）
    ops.hpp            add / sub / mul / neg / cmp / sign / shift
    intrinsics.hpp     x86-64 組み込み関数と可搬フォールバックの切替
  geom/
    point.hpp          IPoint（明示整数点）, HPoint（同次構成点）
    plane.hpp          Plane（平面ベース表現）
    predicates.hpp     幾何述語（§7）
src/
  （ヘッダオンリーが望ましいが、必要なら .cpp を置く）
tests/
  arith/               GMP 差分テスト、性質テスト
  geom/                述語の差分テスト、退化ケーステスト
bench/
  arith_bench.cpp
  pred_bench.cpp
```

---

## 5. `fixed_int` の設計

### 5.1 型

```cpp
namespace krisite::arith {

// N 個の 64bit リムからなる符号付き整数（2 の補数、リトルエンディアン順）
template <std::size_t N>
struct fixed_int {
    std::array<std::uint64_t, N> limb;   // limb[0] が最下位
};

}  // namespace krisite::arith
```

**設計方針**

- POD であること。コンストラクタで初期化コストを払わない
- **動的メモリ確保を一切行わない**
- 例外を投げない（すべて `noexcept`）
- `constexpr` 対応が望ましいが、組み込み関数を使う経路では困難なので必須としない

### 5.2 演算

積は **幅が増える** ことを型で表現してください。これが本設計の要です。

```cpp
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;

template <std::size_t N>
fixed_int<N> add(const fixed_int<N>&, const fixed_int<N>&) noexcept;   // 検査付き

template <std::size_t N>
fixed_int<N + 1> add_widen(const fixed_int<N>&, const fixed_int<N>&) noexcept;

template <std::size_t N>
int sign(const fixed_int<N>&) noexcept;   // -1, 0, +1

template <std::size_t N, std::size_t M>
int cmp(const fixed_int<N>&, const fixed_int<M>&) noexcept;
```

こうすると、述語の式を書いた時点で **必要なリム数がコンパイル時に確定** し、
オーバーフローが型レベルで防げます。§3 の解析を人手で追う必要がなくなります。

### 5.3 実装

- x86-64 (GCC/Clang): `_addcarry_u64` / `_subborrow_u64` / `_mulx_u64`（`<immintrin.h>`）
- 可搬フォールバック: `unsigned __int128`
- **MSVC: `unsigned __int128` が存在しません。** `<intrin.h>` の `_addcarry_u64` /
  `_subborrow_u64` / `_umul128` を使う第3の経路が必要です。「可搬フォールバック」は
  MSVC には可搬ではないので、CI で最初に落ちるのはここだと想定してください。
- 切替は `intrinsics.hpp` に閉じ込め、上位からは見えないようにすること
- ARM64（Apple Silicon）でもビルドが通ること。`__int128` 経路で可

**`sign()` は最上位リムの符号ビットだけで決まります。** 全リムを走査しないでください。
述語は最終的に符号しか使わないので、ここが最頻出のホットパスです。

---

## 6. 幾何型

```cpp
namespace krisite::geom {

// 入力点（格子上の整数座標）
struct IPoint {
    std::int32_t x, y, z;
};

// 平面 N·p + d = 0（同次 4 元ベクトル [a,b,c,d]。d = -N·p1 → §3.1）
template <std::size_t NB, std::size_t DB>
struct Plane {
    arith::fixed_int<NB> a, b, c;   // 法線
    arith::fixed_int<DB> d;         // オフセット d = -N·p1
};

// 同次座標の構成点（実座標は (x/w, y/w, z/w)）
template <std::size_t VB, std::size_t WB>
struct HPoint {
    arith::fixed_int<VB> x, y, z;
    arith::fixed_int<WB> w;         // w != 0 を不変条件とする
};

}  // namespace krisite::geom
```

リム数は §3 の式から `constexpr` で導出してください。手書きの数値を書かないこと。

### 構成

```cpp
// 三角形の支持平面
Plane<...> plane_from_triangle(IPoint p1, IPoint p2, IPoint p3) noexcept;

// 3 平面の交点（Cramer）
// 前提: 3 平面が一点で交わること（呼び出し側が保証）
// w == 0 の場合（平行・共線）は呼び出し側の契約違反 → デバッグ時 assert
HPoint<...> intersect3(const Plane&, const Plane&, const Plane&) noexcept;
```

---

## 7. 述語

### 7.1 Phase 0 で実装するもの

| 名前 | 引数 | 返り値 | 用途 |
|---|---|---|---|
| `orient3d(a,b,c,d)` | `IPoint` × 4 | −1/0/+1 | 入力点の向き |
| `orient2d(a,b,c)` | `IPoint` × 3（投影後） | −1/0/+1 | 2D 向き |
| `side(plane, hp)` | `Plane`, `HPoint` | −1/0/+1 | **最重要**。構成点の平面に対する側 |
| `side(plane, ip)` | `Plane`, `IPoint` | −1/0/+1 | 入力点の平面に対する側 |
| `cmp_h(hp1, hp2, axis)` | `HPoint` × 2 | −1/0/+1 | 同次点の軸別比較 |
| `lex_less(hp1, hp2)` | `HPoint` × 2 | bool | 全順序（頂点テーブル用） |

### 7.2 主要な式

**`side(plane, hpoint)`**

$$\mathrm{sign}(w) \cdot \mathrm{sign}(a x + b y + c z + d\,w)$$

**`cmp_h`** — 除算を使わずに $x_1/w_1$ と $x_2/w_2$ を比較：

$$\mathrm{sign}\left(\frac{x_1}{w_1} - \frac{x_2}{w_2}\right)
= \mathrm{sign}(w_1)\cdot\mathrm{sign}(w_2)\cdot\mathrm{sign}(w_2 x_1 - w_1 x_2)$$

以下の場合は積を計算せずに即答できます（最適化として実装すること）:

- $x_1 = x_2 = 0$ → 0
- $x_1/w_1$ と $x_2/w_2$ の符号が異なる → 符号から即決
- $w_1 = w_2$ → $x_1$ と $x_2$ の比較のみ

**`orient3d(a,b,c,d)`** — 座標差で書くこと（値そのものの行列式にしない）:

$$\det\begin{pmatrix} a_x - d_x & a_y - d_y & a_z - d_z \\ b_x - d_x & \cdots \\ c_x - d_x & \cdots \end{pmatrix}$$

差分形にすると中間値のビット幅が小さくなります。

### 7.3 未確定事項

EMBER の radial sort と分類段が要求する述語は、論文を読んで確定させてください。
**新しい述語を追加する際は必ず §3 と同じ手順でビット幅を導出し、§3.3 の表を更新すること。**

---

## 8. テスト

**このセクションが Phase 0 の実質的な成果物です。** 実装より先にテストの枠組みを作ってください。

### 8.1 GMP 差分テスト

GMP（`mpz_t`）を **正解器** として使います。

> **ライセンス上の重要な制約**: GMP は LGPL です。
> **テストコードからのみ参照し、ライブラリ本体に一切リンクしないでください。**
> CMake オプション `KRISITE_BUILD_TESTS_WITH_GMP=ON`（既定 OFF）でのみ有効化し、
> 配布物には含めないこと。詳細は `THIRD_PARTY_LICENSES.md` を参照。

テスト内容:

1. **算術**: ランダムな `fixed_int` × 10^7 組について `add` / `sub` / `mul` / `cmp` を
   GMP と比較。ビット幅は境界値（全ビット 1、最小値、0、±1）を重点的に。
2. **述語**: ランダムな `IPoint` から平面と構成点を作り、各述語の結果を
   GMP の有理数演算と比較。

### 8.2 退化ケース

**ランダムテストは退化ケースをほとんど生成しません。** 以下を明示的に作ること。

| ケース | 生成方法 |
|---|---|
| 共平面 | 4 点を同一平面上に配置 |
| 共線 | 3 点を同一直線上に |
| 重複頂点 | 同一座標の点を複数 |
| 退化三角形 | 面積 0 |
| 構成点が入力頂点に一致 | 3 平面の交点が既存頂点になるよう逆算 |
| 構成点同士の一致 | 異なる 3 平面組が同一点を与えるよう構成 |
| 格子の端 | 座標が -2^(b-1) と 2^(b-1) − 1（§2 の両端） |
| 極端に細長い三角形 | 1 点だけ大きくずらす |

最後の 2 つは Levy24 が Thingi10K で実際に破綻を報告している構成です。

### 8.3 性質テスト

正解器なしで検証できる不変条件:

- `orient3d(a,b,c,d) == -orient3d(b,a,c,d)`（引数交換で符号反転）
- 偶置換で符号不変、奇置換で反転（全 24 通り）
- `side(plane, p) == 0` ⟺ `p` が `plane` 上
- `cmp_h` が全順序の公理（反対称律・推移律）を満たす
- 同次座標をスカラー倍しても述語の結果が不変

### 8.4 ビット幅の実測

§3 の解析が正しいことを **実測で確認**してください。
デバッグビルドで各演算の実際の使用ビット幅の最大値を記録し、
理論上界を超えないこと、かつ上界が過大でないことを検証します。

---

## 9. 完了条件

以下がすべて満たされたら Phase 0 完了とします。

- [ ] `b = 21` および `b = 26` の両方でビルドが通る
- [ ] Linux(GCC/Clang) / macOS(Apple Silicon) / Windows(MSVC) でビルドが通る
      → **開発コンテナではなく CI（GitHub Actions のマトリクス）で確認すること。**
      ツールチェーンをコンテナに入れるために `settings.json` の deny を緩めないでください。
- [ ] GMP 差分テストが全項目でパス（算術 10^7 件、述語 10^6 件以上）
- [ ] §8.2 の退化ケースがすべてパス
- [ ] §8.3 の性質テストがすべてパス
- [ ] 実測ビット幅が §3 の理論上界内に収まることを確認
- [ ] `side()` のスループットをベンチで計測し、数値を `docs/BENCH.md` に記録
- [ ] ライブラリ本体のバイナリに GMP がリンクされていないことを `ldd` 等で確認

**性能目標は Phase 0 では設定しません。** まず正しさを固め、数値を記録して基準線とします。

---

## 10. Phase 0 の非目標（やらないこと）

以下は Phase 0 の範囲外です。着手しないでください。

- 八分木、Morton コード
- メッシュデータ構造、ファイル入出力
- BSP、arrangement、ブール演算そのもの
- 並列化（`fixed_int` はスレッドセーフであれば十分）
- SIMD 化（まず正しさ。最適化は基準線を取ってから）
- 浮動小数点フィルタ（固定幅で十分速い可能性が高い。実測してから判断）

最後の 2 点は「やらない」というより **「まだ判断材料がない」** という意味です。
Phase 0 のベンチ結果を見てから設計判断します。勝手に入れないでください。
