#!/usr/bin/env python3
"""経路 P（量子化の再実装）と、入力の構造検査・成分体積の集計。

**被検体のブール演算・分類器・述語は 1 つも呼びません。**
幾何の部品は `inspect_geom.py`（合成入力で 73 件の試験を通したもの）を使います。

**段は 1 回の起動につき 1 つ**です（`--stage`）。
外側の駆動が段ごとに期限と資源を掛けるので、**この中では時間を測りません**。

段の順序（`DESIGN-phase5-vertex-level.md` §43.4 の評価順に従います）::

    reconstruct   経路 P の量子化 → 同一性 3・4（C と P のバイト一致）、同一性 1（面数）
    a0            A0-1 / A0-2（T）。**破れたら停止**。A0-3（共線面）は数えるだけ
    volume        同一性 2（符号つき 6 倍体積 と 保存値）、A1〜A4（T と S の両方）
    components    S の成分分解、成分ごとの符号つき 6 倍体積、Σ_j S_j と S(X) の照合

**「参照の妥当性（A0-1 / A0-2）は、体積計算より先」**です。
体積は段 `volume` にあり、A0 は段 `a0` にあります。**順序を入れ替えないでください。**

**結論の限定**（**この道具が示せないこと**）::

  * **A5 を回しません。** 成分が立体の境界であるとは主張しません。
    得られるのは「向き付き三角形の集合についての符号つき 6 倍体積」までです。
  * **巻き数を計算しません。** 内側成分の判定はしません。
  * 経路 C は被検体と同じローダーを使うので、**ローダー自体の誤りはこの比較では出ません。**
  * **経路 P は独立な正解器ではありません。** 手順も浮動小数点の演算順序も
    `loader.hpp` に合わせてあります（合わせないと最下位ビットが正当に変わり、
    バイト一致を要求できません）。**したがって検出できるのは、書き写しの誤り・
    ビルド・エンディアン・逐次化の食い違いだけで、
    C と P が同じ式の誤りを共有する場合は必ず通ります。**
  * 経路 P は**量子化と併合だけ**を対象にします。**変換の生成は検査しません**
    （`std::log` / `std::cos` / `std::sqrt` を通り、libm の実装差でバイト一致を保証できないため）。
  * **同一性 3・4 は、解析したあとの値の比較**です（`cmp` によるバイト比較ではありません）。
  * **踏まなかった経路は検査していません。** 併合・退化面の除去・クランプは、
    発火回数が 0 なら一度も走っていません。**発火回数を出力に書きます。**
  * 保存値との一致は「**いま同じコードを通すと保存値に一致する**」までで、
    ローダーの正しさを示しません。
"""
import argparse, os, struct, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig

# b = 21（`CMakeLists.txt:10`）。`include/krisite/config.hpp:36-37`
COORD_BITS = 21
COORD_MIN = -(1 << (COORD_BITS - 1))
COORD_MAX = (1 << (COORD_BITS - 1)) - 1
FILL = 0.6                      # `loader.hpp` の既定値


# --- 入出力 ---------------------------------------------------------------
def read_kmesh(path):
    """`.kmesh` を読む（`loader.hpp:44-65` と同じ書式。リトルエンディアン固定）。"""
    with open(path, "rb") as f:
        b = f.read()
    if len(b) < 16 or b[:4] != b"KMSH":
        raise ValueError("マジックが違います: %s" % path)
    ver, nv, nf = struct.unpack_from("<III", b, 4)
    if ver != 1:
        raise ValueError("版が違います: %d" % ver)
    off = 16
    need = off + nv * 24 + nf * 12
    if len(b) < need:
        raise ValueError("短すぎます: %d < %d" % (len(b), need))
    V = list(struct.unpack_from("<%dd" % (nv * 3), b, off))
    off += nv * 24
    F = list(struct.unpack_from("<%dI" % (nf * 3), b, off))
    return nv, nf, V, F


def read_transform(path):
    """経路 C が書いた 12 個の `double`（`%a`）を読む。r[0..8] と shift[0..2]。"""
    xs = [float.fromhex(t.strip()) for t in open(path) if t.strip()]
    if len(xs) != 12:
        raise ValueError("変換は 12 個のはずです: %d 個" % len(xs))
    return xs[:9], xs[9:]


def write_bin(path, V, F):
    """`uint32 Nv / uint32 Nf / int32 x 3Nv / uint32 x 3Nf`。リトルエンディアン固定。"""
    with open(path, "wb") as f:
        f.write(struct.pack("<II", len(V), len(F)))
        f.write(struct.pack("<%di" % (len(V) * 3), *[c for v in V for c in v]))
        f.write(struct.pack("<%dI" % (len(F) * 3), *[c for t in F for c in t]))


def read_bin(path):
    with open(path, "rb") as f:
        b = f.read()
    nv, nf = struct.unpack_from("<II", b, 0)
    V = list(struct.unpack_from("<%di" % (nv * 3), b, 8))
    F = list(struct.unpack_from("<%dI" % (nf * 3), b, 8 + nv * 12))
    return ([tuple(V[3 * i:3 * i + 3]) for i in range(nv)],
            [tuple(F[3 * i:3 * i + 3]) for i in range(nf)])


# --- 経路 P（量子化と併合の再実装。`loader.hpp:102-165`）-------------------
def quantize_p(nv, nf, V, F, r, shift, fill=FILL):
    """**手順と、浮動小数点の演算順序まで `loader.hpp` に合わせます。**

    戻り値: `(頂点, 面, 併合した頂点数, 落とした退化面数, 範囲外か, 発火回数)`

    **発火回数**は「その経路が実際に走ったか」の数です。
    **0 なら、その経路はこの入力では検査されていません**（`clamp` / `half`）。
    """
    # 1. 回転（行優先。**左から順に足します**）
    p = [0.0] * (nv * 3)
    for i in range(nv):
        x, y, z = V[3 * i], V[3 * i + 1], V[3 * i + 2]
        for k in range(3):
            p[3 * i + k] = r[3 * k] * x + r[3 * k + 1] * y + r[3 * k + 2] * z

    # 2. AABB
    lo = [p[0], p[1], p[2]]
    hi = [p[0], p[1], p[2]]
    for i in range(1, nv):
        for k in range(3):
            t = p[3 * i + k]
            if t < lo[k]:
                lo[k] = t
            if t > hi[k]:
                hi[k] = t

    # 3. ext / range / scale
    ext = 0.0
    for k in range(3):
        d = hi[k] - lo[k]
        if d > ext:
            ext = d
    if not (ext > 0):
        return [], [], 0, 0, False, {"clamp": 0, "half": 0}
    rng = float(COORD_MAX) * fill
    scale = rng / ext

    # 4. 量子化 + 同一格子点の併合
    #    **代表は最小の元の添字**、**新しい添字は元の添字の昇順**（`loader.hpp:139-148`）
    mid = [0.5 * (lo[k] + hi[k]) for k in range(3)]
    at, remap, VQ = {}, [0] * nv, []
    merged = 0
    fired = {"clamp": 0, "half": 0}
    for i in range(nv):
        c = []
        for k in range(3):
            val = (p[3 * i + k] - mid[k]) * scale + shift[k] * rng
            # `quantize_one`: nearbyint（最近接・偶数優先）→ クランプ
            if abs(val - int(val)) == 0.5:
                fired["half"] += 1      # ちょうど .5（丸めの半数規則が効く点）
            q = _nearbyint(val)
            if q < COORD_MIN:
                q = COORD_MIN
                fired["clamp"] += 1
            elif q > COORD_MAX:
                q = COORD_MAX
                fired["clamp"] += 1
            c.append(q)
        key = (c[0], c[1], c[2])
        j = at.get(key)
        if j is None:
            at[key] = len(VQ)
            remap[i] = len(VQ)
            VQ.append(key)
        else:
            remap[i] = j
            merged += 1
    oor = any(not (COORD_MIN <= x <= COORD_MAX) for v in VQ for x in v)

    # 5. 面。**頂点が潰れたものは落とす**（`loader.hpp:158`。面積ではなく添字で判定）
    FQ, dropped = [], 0
    for i in range(nf):
        a, b, c = remap[F[3 * i]], remap[F[3 * i + 1]], remap[F[3 * i + 2]]
        if a == b or b == c or c == a:
            dropped += 1
            continue
        FQ.append((a, b, c))
    return VQ, FQ, merged, dropped, oor, fired


def _nearbyint(x):
    """C の `nearbyint`（既定の丸めモード = 最近接・偶数優先）。

    **`int(x + 0.5)` や `floor(x + 0.5)` は、ちょうど .5 の点で食い違います。**
    Python の `round()` は偶数優先なので、そちらを使います。
    """
    return int(round(x))


# --- 保存値 ---------------------------------------------------------------
def saved_row(path, key):
    """診断の保存物から 1 行を取り、列を名前で返す。

    **添字は 0 起点**です（`t[5]` は 1 起点で数えた列 6）。
    **列が挿入・並べ替えされると黙って別の列を読む**ので、列数を下から検査します。
    **27〜34 列（1 起点）は有理数になり得ます**（`分子/分母`）。
    この対の列 28 は整数ですが、**他の対では `int()` が例外を投げます**。
    **推測で埋めず、そのまま落とします。**
    """
    need = 34                      # `run_gmp_diag.py:27` の固定列数
    for line in open(path):
        t = line.split()
        if t and t[0] == key:
            if len(t) < need:
                raise ValueError("列が %d 個しかありません（%d 以上のはず）" % (len(t), need))
            return {"faces_A": int(t[5]), "faces_B": int(t[6]),
                    "R": int(t[27]), "S_A": int(t[32]), "S_B": int(t[33]),
                    "ncols": len(t)}
    raise KeyError("保存物に %s がありません" % key)


# --- 幾何（S の作り方）-----------------------------------------------------
def make_S(P, F):
    """T から**外積が零の面（共線面）を除いた**三角形の集合。頂点配列は共有します。"""
    keep = [f for f in F if ig.normal(P, f) != (0, 0, 0)]
    return keep


def volume6(P, F):
    return ig.signed_volume6(P, F)


# --- 段 -------------------------------------------------------------------
def stage_reconstruct(a, out, log):
    ok = True
    ident = {}
    saved = saved_row(a.saved, a.key)
    log("保存物 %s / 鍵 %s / 列数 %d" % (a.saved, a.key, saved["ncols"]))
    log("保存値: 面数 A=%d B=%d / S(A)=%d / S(B)=%d / R=%d"
        % (saved["faces_A"], saved["faces_B"], saved["S_A"], saved["S_B"], saved["R"]))
    for name, mid, src in (("A", a.a_id, a.a_kmesh), ("B", b_id_of(a), a.b_kmesh)):
        r, sh = read_transform(os.path.join(out, "C_%s_transform.hex" % name))
        nv, nf, V, F = read_kmesh(src)
        VQ, FQ, merged, dropped, oor, fired = quantize_p(nv, nf, V, F, r, sh)
        write_bin(os.path.join(out, "quantized_P_%s.bin" % name), VQ, FQ)
        VC, FC = read_bin(os.path.join(out, "C_%s_quantized.bin" % name))
        log("[%s] 模型 %s / 原本 頂点 %d 面 %d" % (name, mid, nv, nf))
        log("[%s] 経路 P: 頂点 %d / 面 %d / 併合 %d / 落とした退化面 %d / 範囲外 %s"
            % (name, len(VQ), len(FQ), merged, dropped, oor))
        log("[%s] 経路 C: 頂点 %d / 面 %d" % (name, len(VC), len(FC)))
        log("[%s] 発火回数: 併合 %d / 退化面の除去 %d / クランプ %d / ちょうど .5 %d"
            % (name, merged, dropped, fired["clamp"], fired["half"]))
        log("[%s]   **発火 0 の経路は、この入力では検査されていません。**" % name)
        i3 = (VQ == VC)
        i4 = (FQ == FC)
        want = saved["faces_A"] if name == "A" else saved["faces_B"]
        i1 = (len(FC) == want)
        ident["%s_3" % name], ident["%s_4" % name], ident["%s_1" % name] = i3, i4, i1
        log("[%s] 同一性 3（整数頂点列が C と P で一致）: %s" % (name, "一致" if i3 else "不一致"))
        log("[%s] 同一性 4（向き付き面列が C と P で一致）: %s" % (name, "一致" if i4 else "不一致"))
        log("[%s] 同一性 1（面数 %d と保存値 %d）: %s"
            % (name, len(FC), want, "一致" if i1 else "不一致"))
        if not i3:
            d = [k for k in range(min(len(VQ), len(VC))) if VQ[k] != VC[k]]
            log("[%s]   頂点の食い違い %d 件。先頭: %s"
                % (name, len(d) + abs(len(VQ) - len(VC)),
                   "" if not d else "添字 %d  P=%s C=%s" % (d[0], VQ[d[0]], VC[d[0]])))
        ok = ok and i3 and i4 and i1
    log("判定: %s" % ("同一性 1・3・4 すべて一致" if ok else "★ 不一致あり。後続を起動しません"))
    return 0 if ok else 1


def b_id_of(a):
    return a.b_id


def stage_a0(a, out, log):
    ok = True
    # ★ 以降の段が読むのは【経路 P の配列】です。
    #   段 reconstruct で同一性 3・4（C と P の一致）が通っているときだけ、
    #   これは経路 C の配列と同じものです。**段を単独で回すと、その保証がありません。**
    log("読む配列: quantized_P_*.bin（経路 P）。"
        "段 reconstruct の同一性 3・4 が通っているときに限り、経路 C と同じ内容です。")
    for name in ("A", "B"):
        P, F = read_bin(os.path.join(out, "quantized_P_%s.bin" % name))
        oor, dup, col = ig.check_a0(P, F)
        log("[%s] A0-1（頂点参照が範囲内）: %s（違反 %d 枚）"
            % (name, "通過" if not oor else "★ 破れ", len(oor)))
        log("[%s] A0-2（3 頂点が相異なる）: %s（違反 %d 枚）"
            % (name, "通過" if not dup else "★ 破れ", len(dup)))
        if oor or dup:
            log("[%s] ★ 参照の妥当性が破れました。体積計算へ進みません。" % name)
            ok = False
            continue
        log("[%s] A0-3（共線面）: %d 枚（**止めません。数えるだけ**）" % (name, len(col)))
        used = {x for f in F for x in f}
        log("[%s] 頂点 %d 個のうち、面から参照されないもの %d 個"
            % (name, len(P), len(P) - len(used)))
    log("判定: %s" % ("A0-1・A0-2 は T の両方で通過" if ok else "★ 参照の妥当性が破れました"))
    return 0 if ok else 1


def stage_volume(a, out, log):
    saved = saved_row(a.saved, a.key)
    ok = True
    for name in ("A", "B"):
        P, F = read_bin(os.path.join(out, "quantized_P_%s.bin" % name))
        S = make_S(P, F)
        # 同一性 2 は T で。**共線面の符号つき体積は恒等的に 0** なので S でも同じ値です。
        vt = volume6(P, F)
        vs = volume6(P, S)
        want = saved["S_A"] if name == "A" else saved["S_B"]
        i2 = (vt == want)
        log("[%s] 符号つき 6 倍体積: T で %d / S で %d（差 %d）" % (name, vt, vs, vt - vs))
        log("[%s] 同一性 2（T の体積 と 保存値 %d）: %s" % (name, want, "一致" if i2 else "★ 不一致"))
        ok = ok and i2
        for tag, G in (("T", F), ("S", S)):
            a1 = ig.check_a1(G)
            a2 = ig.check_a2(G)
            a3 = ig.check_a3(G) if not ig.check_a0(P, G)[1] else [("A0-2 未通過", None)]
            a4 = ig.check_a4(G)
            log("[%s] A1(%s) 各無向辺に面が 2 枚: %s（違反 %d 本）"
                % (name, tag, "通過" if not a1 else "破れ", len(a1)))
            log("[%s] A2(%s) 各有向辺が 1 回: %s（違反 %d 本）"
                % (name, tag, "通過" if not a2 else "破れ", len(a2)))
            log("[%s] A3(%s) 頂点リンクが単一の閉路: %s（違反 %d 個）"
                % (name, tag, "通過" if not a3 else "破れ", len(a3)))
            log("[%s] A4(%s) 重複面なし: %s（重複 %d 組）"
                % (name, tag, "通過" if not a4 else "破れ", len(a4)))
    log("判定: %s" % ("同一性 2 は A・B とも一致" if ok
                    else "★ 同一性 2 が不一致。後続を起動しません"))
    log("**A1〜A4 の破れは記録するだけで止めません。ただし閉立体とは呼びません。**")
    return 0 if ok else 1


def stage_components(a, out, log):
    saved = saved_row(a.saved, a.key)
    ok = True
    log("読む配列: quantized_P_*.bin（経路 P）。"
        "段 reconstruct の同一性 3・4 が通っているときに限り、経路 C と同じ内容です。")
    for name in ("A", "B"):
        P, F = read_bin(os.path.join(out, "quantized_P_%s.bin" % name))
        S = make_S(P, F)
        comps = ig.components(S)
        tot = 0
        log("[%s] S の面数 %d / 連結成分 %d 個" % (name, len(S), len(comps)))
        for j, g in enumerate(comps):
            G = [S[i] for i in g]
            v = volume6(P, G)
            tot += v
            a1 = ig.check_a1(G)
            a2 = ig.check_a2(G)
            log("[%s]   成分 %d: 面 %d 枚 / 符号つき 6 倍体積 %d / A1 違反 %d 本 / A2 違反 %d 本"
                % (name, j, len(G), v, len(a1), len(a2)))
        want = saved["S_A"] if name == "A" else saved["S_B"]
        vs = volume6(P, S)
        log("[%s] Σ_j S_j = %d / S(S) = %d / 保存値 %d" % (name, tot, vs, want))
        log("[%s] 内部整合（Σ_j S_j = S(S)）: %s" % (name, "一致" if tot == vs else "★ 不一致"))
        log("[%s] 保存値との一致: %s" % (name, "一致" if tot == want else "★ 不一致"))
        # ★ 判定を終了値へ届けます。**ログに書くだけでは駆動が成功と読みます。**
        ok = ok and (tot == vs) and (tot == want)
    log("**A5 を回していないので、成分が立体の境界であるとは主張しません。**")
    log("**巻き数を計算していないので、内側成分の判定はしません。**")
    log("判定: %s" % ("成分の和は内部整合・保存値ともに一致" if ok else "★ 不一致あり"))
    return 0 if ok else 1


STAGES = {"reconstruct": stage_reconstruct, "a0": stage_a0,
          "volume": stage_volume, "components": stage_components}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True, choices=sorted(STAGES))
    p.add_argument("--out", required=True, help="出力の置き場所")
    p.add_argument("--a-kmesh", default="data/thingi10k/kmesh/250394.kmesh")
    p.add_argument("--b-kmesh", default="data/thingi10k/kmesh/45413.kmesh")
    p.add_argument("--a-id", default="250394")
    p.add_argument("--b-id", default="45413")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    path = os.path.join(a.out, "stage_%s.txt" % a.stage)
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    log("=== 段 %s ===" % a.stage)
    try:
        rc = STAGES[a.stage](a, a.out, log)
    finally:
        # **成果物は必ず書きます。** 途中で落ちても、どこまで進んだかが残ります。
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
