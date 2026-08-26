# Krisite — ベンチマーク結果

SPEC-phase0.md §9「`side()` のスループットをベンチで計測し、数値を `docs/BENCH.md` に記録」。

> **これは基準線であって目標値ではありません。** Phase 0 では性能目標を設定しません
> （SPEC §9）。SIMD 化・浮動小数点フィルタは、この数値を見てから設計判断します（SPEC §10）。

## 測定環境

| 項目 | 値 |
|---|---|
| 日付 | 2026-08-26（SPEC 2026-08 改訂に追従後） |
| CPU | AMD Ryzen 7 9800X3D (8C/16T) |
| OS | Ubuntu 26.04 (Docker, Linux 6.12.67) |
| コンパイラ | clang 21.1.0（`zig c++` 経由） |
| ビルド | `Release` (`-O3 -DNDEBUG`), `KRISITE_CHECKED_ARITH=OFF` |
| `KRISITE_COORD_BITS` | 21 |
| 測定方法 | 各項目 3 回計測して最小値。入力は L1 に収まる 2048 要素のプールを巡回 |

> **注意**: これは開発コンテナの数値です。CI ランナー（GitHub Actions）の数値は大きく
> 揺れるため合否判定には使いません（`.github/workflows/ci.yml` の `bench` ジョブ）。

再現手順:

```bash
cmake -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DKRISITE_CHECKED_ARITH=OFF -DKRISITE_BUILD_BENCH=ON
cmake --build build-rel
./build-rel/bench/pred_bench
./build-rel/bench/arith_bench
```

## 述語（`bench/pred_bench.cpp`）

`side` は 209 ビット / 4 リム、`cmp_h` は 300 ビット / 5 リムで評価（SPEC §3.3）。

| 演算 | ns/op | Mops/s |
|---|---:|---:|
| **`side(plane, HPoint)`** ★最重要 | **7.80** | **128.3** |
| `side(plane, IPoint)` | 2.35 | 426.0 |
| `orient3d(IPoint × 4)` | 12.87 | 77.7 |
| `orient2d(IPoint × 3, Z)` | 2.08 | 481.6 |
| `cmp_h(HPoint × 2, X)` | 11.25 | 88.9 |
| `lex_less(HPoint × 2)` | 11.42 | 87.6 |
| `plane_from_triangle` | 3.38 | 295.7 |
| `intersect3(Plane × 3)` | 261.76 | 3.8 |

## 固定幅整数（`bench/arith_bench.cpp`）

| 演算 | 1 リム (64bit) | 2 リム (128bit) | 4 リム (256bit) |
|---|---:|---:|---:|
| `add` | 0.23 ns | 0.38 ns | 0.70 ns |
| `mul` | 0.56 ns | 1.90 ns | 11.63 ns |
| `cmp` | 0.43 ns | 0.41 ns | 0.41 ns |
| `sign` | 0.41 ns | 0.40 ns | 0.36 ns |

`cmp` と `sign` がリム数にほぼ依存しないのは、符号ビットと最上位リムで決着する経路が
効いているため（SPEC §5.3）。`mul` はリム数の 2 乗で伸びており、筆算の計算量と一致します。

## 読み取れること

**`side()` は 7.8 ns（128 Mops/s）。** 浮動小数点フィルタを入れるかどうかの判断材料として、
これが基準線になります。EMBER が「フィルタ不要」と主張する根拠と整合する水準です。

**`cmp_h` は `side` より重い（11.3 ns）。** 被符号値が 300 ビット（5 リム）と Phase 0 で
最大幅であること（SPEC §3.3）が素直に出ています。加えて、リム単位の積が生む 7 リムを
§3 の上界である 5 リムへ詰める `resize` のぶん、約 2 ns 上乗せされています。
この詰め直しは検査ビルドで上界超過を捕まえるための投資であり、外すなら
`KRISITE_CHECKED_ARITH=OFF` のときだけ省く形にするのが筋です（Phase 0 では手を入れません）。

**`intersect3()` が 262 ns で突出して重い。** これは実装の選択によるもので、
`intersect3` は 3x3 行列式を 4 回計算しますが、列ごとに幅が違う（法線 1 リム、
オフセット 2 リム）ため、**正しさを優先して全列を 2 リムにそろえてから** 7 リム幅で
計算しています（`geom/plane.hpp`）。列ごとの幅を保った混合幅の行列式にすれば
数倍は縮むはずですが、Phase 0 では手を入れません（SPEC §0「性能より正しさを優先」）。

なお EMBER の想定では構成点は一度作って何度も `side` にかけるので、
`intersect3` の重さが支配的になるとは限りません。**Phase 1 で実際の比率を測ってから**
最適化対象を決めてください。

## 次に測るべきもの

Phase 0 の完了条件には含まれませんが、Phase 1 の設計判断に必要になります。

- `side()` の分岐予測が外れる条件（構成点が平面上に載る割合が高い実データ）
- キャッシュミスを含む現実的なアクセスパターン（プール巡回ではなくランダム参照）
- 混合幅 `det3` にしたときの `intersect3` の改善幅
- `cmp_h` の `resize` を検査ビルド限定にしたときの差
