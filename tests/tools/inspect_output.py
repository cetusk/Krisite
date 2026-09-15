#!/usr/bin/env python3
r"""案 I の経路 P — 出力の読込、$E0$・$E1$・退化面（全数）、限定標本の検査と単価測定。

**被検体を 1 行も呼びません。** 幾何は整数・有理数で厳密に計算します。

段（1 回の起動につき 1 つ）::

    identity   I-1 受理済み配列の SHA256 照合（**被検体へ渡すのと同じファイル**）
    anchor     I-3 錨の照合（4 演算ハッシュと面数）。**同一性の証明ではありません**
    measure    I-4 E0・E1・退化面（**全数**）と、限定標本の検査・単価測定

**結論の限定**（`DESIGN-phase5-vertex-level.md` §44）::

  * **$E2$〜$E5$・C1・C2 の全数判定をせず、$Q$ の認定もしません。**
    **$E0$・$E1$・退化面の計数だけが全数です。**
  * **単価は「その出力と、選んだ標本」についての観測値**です。
    **1 標本の値を全体の保証にしません**（有理数の幅・約分・候補の数・通る分岐に依存）。
  * **錨の不一致は「同じ出力である根拠が無い」**であって、「被検体が壊れた」ではありません。

**記録は 1 行ごとに書き出して flush します**（打ち切られても、そこまでが残るように）。
"""
import argparse, hashlib, os, struct, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EXPECT_IN = {
    "A": "d99d6d3b45132ee22c6ae83e391c47c8e8b2f0215b2c851edbab044921f5cf6b",
    "B": "fd9a03316f3c49e9d5b13a77fcc41ad6879dbdac0c6f0c8cbf69b50a97610596",
}
OPS = ["union", "isect", "diff_ab", "diff_ba"]


# --- 読込 -----------------------------------------------------------------
def read_meta(path):
    m = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith("op"):
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            m[k] = v
    return m


def _signed(limbs):
    """リムの並び（little-endian の `uint64`）を 2 の補数の整数に直します。"""
    v = 0
    for i, w in enumerate(limbs):
        v |= w << (64 * i)
    bits = 64 * len(limbs)
    if v >> (bits - 1):
        v -= 1 << bits
    return v


def read_soup(path, nx, nw):
    """`dump_boolean` が書いた出力を読みます（§44.8 の書式）。

    戻り値: `(頂点の同次座標 [(X,Y,Z,W)] （整数）, 面 [(a,b,c)], sha256)`
    """
    with open(path, "rb") as f:
        b = f.read()
    nv, nf = struct.unpack_from("<II", b, 0)
    per = 3 * nx + nw
    need = 8 + nv * per * 8 + nf * 12
    if len(b) != need:
        raise ValueError("長さが合いません: %d != %d（%s）" % (len(b), need, path))
    V = []
    off = 8
    raw = struct.unpack_from("<%dQ" % (nv * per), b, off)
    for i in range(nv):
        base = i * per
        X = _signed(raw[base:base + nx])
        Y = _signed(raw[base + nx:base + 2 * nx])
        Z = _signed(raw[base + 2 * nx:base + 3 * nx])
        W = _signed(raw[base + 3 * nx:base + per])
        V.append((X, Y, Z, W))
    off = 8 + nv * per * 8
    F = struct.unpack_from("<%dI" % (nf * 3), b, off)
    F = [tuple(F[3 * i:3 * i + 3]) for i in range(nf)]
    return V, F, hashlib.sha256(b).hexdigest()


# --- 幾何（整数・有理数）--------------------------------------------------
def dir_int(p, q):
    """同次座標 p→q の方向ベクトル（整数。共通の分母 $w_p w_q$ を落としたもの）。

    $q/w_q - p/w_p = (q_k w_p - p_k w_q)/(w_p w_q)$ なので、
    **分母は 3 成分で共通**です。**向きの符号は $w_p w_q$ の符号に依ります。**
    """
    return tuple(q[k] * p[3] - p[k] * q[3] for k in range(3))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def degenerate(V, f):
    """外積が零か（面積 0）。**整数だけで判定します。** $w=0$ の頂点があれば `None`。"""
    p, q, r = V[f[0]], V[f[1]], V[f[2]]
    if p[3] == 0 or q[3] == 0 or r[3] == 0:
        return None
    return cross(dir_int(p, q), dir_int(p, r)) == (0, 0, 0)


def pt(v):
    """同次座標を有理点に。"""
    return (Fr(v[0], v[3]), Fr(v[1], v[3]), Fr(v[2], v[3]))


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def classify_e4(P, f1, f2):
    """§44.3 の 4 行に分類。**3 分類に当たらなければ「未判定」**（正規扱いにしません）。"""
    import inspect_geom as ig
    n1, n2 = ig.normal(P, f1), ig.normal(P, f2)
    if n1 == (0, 0, 0) or n2 == (0, 0, 0):
        return "未判定（退化面。案 I では扱いません）"
    cop = ig.cross(n1, n2) == (0, 0, 0) and ig.dot(n1, P[f2[0]]) == ig.dot(n1, P[f1[0]])
    T1 = [P[i] for i in f1]
    T2 = [P[i] for i in f2]
    if cop:
        ax = max(range(3), key=lambda k: abs(n1[k]))
        ij = [k for k in range(3) if k != ax]
        poly = ig._inter2(T1, T2, ij)
        if not poly:
            return "正規（交わり空）"
        if ig._area2(poly, ij) != 0:
            return "不正（正の面積を共有）"
        return _shared_or_undecided(P, f1, f2, poly)
    s, ok, kind = ig.pair_a5(P, f1, f2)
    if ok:
        return "正規（交わり空、または共有する単体と一致）"
    if kind == "交差":
        return "不正（真の交差）"
    return "未判定（3 分類に当たりません）"


def _shared_or_undecided(P, f1, f2, poly):
    """交わりが、**座標で見た**共有頂点・辺と一致するか。**添字では判定しません。**"""
    sh = [P[i] for i in f1 if any(P[i] == P[j] for j in f2)]
    if len(poly) == len(sh) and sorted(poly) == sorted(sh):
        return "正規（共有する単体と一致）"
    return "未判定（3 分類に当たりません）"


def c2_sample(P, f, a):
    """支持平面が入力三角形のどれかと一致し、**単一の**入力三角形に載るか。"""
    import inspect_geom as ig
    import struct as _st
    n = ig.normal(P, f)
    d = ig.dot(n, P[f[0]])
    for name, path in (("A", a.in_a), ("B", a.in_b)):
        with open(path, "rb") as fp:
            b = fp.read()
        nv, nf = _st.unpack_from("<II", b, 0)
        Q = _st.unpack_from("<%di" % (nv * 3), b, 8)
        Q = [tuple(Q[3 * i:3 * i + 3]) for i in range(nv)]
        G = _st.unpack_from("<%dI" % (nf * 3), b, 8 + nv * 12)
        G = [tuple(G[3 * i:3 * i + 3]) for i in range(nf)]
        for gi, g in enumerate(G):
            m = ig.normal(Q, g)
            if m == (0, 0, 0):
                continue
            # 平面の一致（法線が平行、かつ offset も比例）
            if ig.cross(n, m) != (0, 0, 0):
                continue
            e = ig.dot(m, Q[g[0]])
            if n[0] * e != m[0] * d or n[1] * e != m[1] * d or n[2] * e != m[2] * d:
                continue
            ax = max(range(3), key=lambda k: abs(m[k]))
            ij = [k for k in range(3) if k != ax]
            T = [Q[i] for i in g]
            if all(ig._inside2(P[i], T, ij) for i in f):
                same = ig.dot(n, m) > 0
                return ("単一の入力三角形に収まる（%s の面 %d、向きは%s）"
                        % (name, gi, "同じ" if same else "★ 逆"))
    return "未判定（単一の入力三角形に収まりません）"


def c1_sample(P, F, a, log):
    """非退化な直線 1 本で、$[w_X>0]$ と期待値、および $w_X$ の値を比べます。"""
    import inspect_geom as ig
    import struct as _st
    ins = {}
    for name, path in (("A", a.in_a), ("B", a.in_b)):
        with open(path, "rb") as fp:
            b = fp.read()
        nv, nf = _st.unpack_from("<II", b, 0)
        Q = _st.unpack_from("<%di" % (nv * 3), b, 8)
        Q = [tuple(Q[3 * i:3 * i + 3]) for i in range(nv)]
        G = _st.unpack_from("<%dI" % (nf * 3), b, 8 + nv * 12)
        ins[name] = (Q, [tuple(G[3 * i:3 * i + 3]) for i in range(nf)])

    st = [1]

    def nxt():
        x = st[0]
        x ^= (x >> 12) & 0xFFFFFFFFFFFFFFFF
        x ^= (x << 25) & 0xFFFFFFFFFFFFFFFF
        x ^= (x >> 27) & 0xFFFFFFFFFFFFFFFF
        st[0] = x & 0xFFFFFFFFFFFFFFFF
        return ((st[0] * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF) >> 11

    for attempt in range(1, 33):
        o = tuple(Fr(nxt() % 400001 - 200000) for _ in range(3))
        d = tuple(nxt() % 199 - 99 for _ in range(3))
        if d == (0, 0, 0):
            continue
        ev = events(P, F, o, d)
        if ev is None:
            continue
        ea = events(ins["A"][0], ins["A"][1], o, d)
        eb = events(ins["B"][0], ins["B"][1], o, d)
        if ea is None or eb is None:
            continue
        log("[union] C1: 試行 %d 本目で非退化な直線を得ました（交点 X %d / A %d / B %d）"
            % (attempt, len(ev), len(ea), len(eb)))
        ts = sorted({t for t, _ in ev + ea + eb})
        bad_c, bad_b = 0, 0
        for k in range(len(ts) + 1):
            m = (ts[0] - 1 if k == 0 else
                 (ts[-1] + 1 if k == len(ts) else (ts[k - 1] + ts[k]) / 2))
            wx = sum(sg for t, sg in ev if t < m)
            wa = sum(sg for t, sg in ea if t < m)
            wb = sum(sg for t, sg in eb if t < m)
            want = 1 if (wa > 0 or wb > 0) else 0
            if (1 if wx > 0 else 0) != want:
                bad_c += 1
            elif wx != want:
                bad_b += 1
        log("[union] C1: 区間 %d 個 / `C1-c` %d 個 / `C1-b` %d 個" % (len(ts) + 1, bad_c, bad_b))
        log("[union]   **通った分岐**: 非退化の直線 1 本、交点 %d 個。"
            "**この直線が通らない食い違いは見えません。**" % len(ev))
        return "`C1-c` %d 個 / `C1-b` %d 個" % (bad_c, bad_b)
    return "**未測定**（32 回引いても非退化な直線が得られませんでした）"


def events(P, F, o, d):
    """直線 $o+td$ と面の交点を `(t, 符号)` で返す。**退化に当たったら `None`。**"""
    import inspect_geom as ig
    out = []
    for f in F:
        n = ig.normal(P, f)
        if n == (0, 0, 0):
            continue
        den = ig.dot(n, d)
        num = ig.dot(n, tuple(P[f[0]][k] - o[k] for k in range(3)))
        if den == 0:
            if num == 0:
                return None                    # 直線が面の平面に載る
            continue
        t = Fr(num, den)
        x = tuple(o[k] + t * d[k] for k in range(3))
        T = [P[i] for i in f]
        ax = max(range(3), key=lambda k: abs(n[k]))
        ij = [k for k in range(3) if k != ax]
        if ig._strict_inside2(x, T, ij):
            # ★ 符号の向き。**直線の【手前側】の交差を足して巻き数にします。**
            #   外向きの面を $d$ の向きに横切る（`den > 0`）と【外へ出る】ので $-1$。
            #   逆向きに横切る（`den < 0`）と【中へ入る】ので $+1$。
            #   **これを逆にすると、箱の内側で巻き数が $-1$ になります**
            #   （合成対照が検出しました）。
            out.append((t, -1 if den > 0 else 1))
        elif ig._inside2(x, T, ij):
            return None                        # 辺・頂点に当たった
    return out


# --- 段 -------------------------------------------------------------------
def stage_identity(a, log):
    ok = True
    for name, path in (("A", a.in_a), ("B", a.in_b)):
        with open(path, "rb") as f:
            h = hashlib.sha256(f.read()).hexdigest()
        good = (h == EXPECT_IN[name])
        log("[%s] %s" % (name, path))
        log("[%s]   sha256 %s / 期待 %s / %s" % (name, h, EXPECT_IN[name],
                                                "一致" if good else "★ 不一致"))
        ok = ok and good
    log("**この同じファイルを被検体へ渡します。再量子化しません。**")
    log("判定: %s" % ("入力の同一性は成立" if ok else "★ 破れ。後続を起動しません"))
    return 0 if ok else 1


def stage_anchor(a, log):
    meta = read_meta(os.path.join(a.out, "b_run_meta.txt"))
    saved = saved_row(a.saved, a.key)
    log("保存値: 4 演算ハッシュ %s / 面数 %s" % (saved["hash4"], saved["tri"]))
    got_hash = meta.get("hash4", "")
    got_tri = []
    for line in open(os.path.join(a.out, "b_run_meta.txt")):
        if line.startswith("op"):
            for t in line.split():
                if t.startswith("tri="):
                    got_tri.append(int(t[4:]))
    log("再実行: 4 演算ハッシュ %s / 面数 %s" % (got_hash, got_tri))
    h_ok = (got_hash == saved["hash4"])
    t_ok = (got_tri == saved["tri"])
    log("錨（4 演算ハッシュ）: %s" % ("一致" if h_ok else "★ 不一致"))
    log("錨（面数）: %s" % ("一致" if t_ok else "★ 不一致"))
    log("**錨は照合材料であって、同一性の証明ではありません**（列 10 は 64 ビットの合成）。")
    log("**不一致は「同じ出力である根拠が無い」であって、「被検体が壊れた」ではありません。**")
    log("**錨の不一致では後続を止めません**（単価は錨と無関係です）。")
    log("判定: 錨 = %s" % ("一致" if (h_ok and t_ok) else "不一致"))
    return 0


def saved_row(path, key):
    need = 34
    for line in open(path):
        t = line.split()
        if t and t[0] == key:
            if len(t) < need:
                raise ValueError("列が %d 個しかありません" % len(t))
            return {"hash4": t[9], "tri": [int(t[14]), int(t[15]), int(t[16]), int(t[17])],
                    "vol_union": t[28]}
    raise KeyError("保存物に %s がありません" % key)


def stage_measure(a, log):
    meta = read_meta(os.path.join(a.out, "b_run_meta.txt"))
    nx, nw = int(meta["kHomoXyz"]), int(meta["kHomoW"])
    log("リム数: kHomoXyz %d / kHomoW %d" % (nx, nw))
    data = {}
    total_tri = 0
    # ---- E0（全数。停止対象）----
    stop = False
    for op in OPS:
        V, F, sha = read_soup(os.path.join(a.out, "b_out_%s.bin" % op), nx, nw)
        data[op] = (V, F)
        total_tri += len(F)
        oor = [i for i, f in enumerate(F) if max(f) >= len(V) or min(f) < 0]
        dup = [i for i, f in enumerate(F) if len(set(f)) < 3]
        w0 = [i for i, v in enumerate(V) if v[3] == 0]
        log("[%s] 頂点 %d / 面 %d / sha256 %s" % (op, len(V), len(F), sha))
        log("[%s] E0 参照範囲 %s（違反 %d）/ 3 頂点相異 %s（違反 %d）/ 同次分母 0 %s（%d 個）"
            % (op, "通過" if not oor else "★ 破れ", len(oor),
               "通過" if not dup else "★ 破れ", len(dup),
               "無し" if not w0 else "★ あり", len(w0)))
        if oor or dup or w0:
            stop = True
    log("対象の合計 %d 三角形" % total_tri)
    if stop:
        log("★ E0 が破れました。後続の計算が定義できません。ここで止めます。")
        return 1

    # ---- E1（全数。C1 の前提）----
    e1_ok = {}
    for op in OPS:
        V, F = data[op]
        cnt = {}
        for f in F:
            for k in range(3):
                cnt[(f[k], f[(k + 1) % 3])] = cnt.get((f[k], f[(k + 1) % 3]), 0) + 1
        bad = [e for e in cnt if cnt[e] != cnt.get((e[1], e[0]), 0)]
        e1_ok[op] = not bad
        log("[%s] E1 各無向辺で #(u,v)=#(v,u): %s（違反 %d 本）"
            % (op, "通過" if not bad else "★ 破れ", len(bad)))
    log("**E1 が破れた出力では、C1 の巻き数判定を【未評価】にします**（階段関数を巻き数と"
        "読むには閉じた向き付き曲面の前提が要るため）。")

    # ---- 退化面（全数）----
    t0 = time.process_time()
    deg = {}
    for op in OPS:
        V, F = data[op]
        deg[op] = sum(1 for f in F if degenerate(V, f))
        log("[%s] 退化面（外積が零）: %d 枚" % (op, deg[op]))
    dt = time.process_time() - t0
    log("退化面の計数: 合計 %d 三角形で CPU %.3f 秒（1 枚あたり %.3g 秒）"
        % (total_tri, dt, dt / max(total_tri, 1)))
    log("**発火回数**: 退化面の判定を %d 回行いました" % total_tri)

    # ---- $E4$: 両方とも非退化な面対 1 組（∪。AABB が交わる対のうち添字の昇順で最初）----
    V, F = data["union"]
    P = [pt(v) for v in V]
    nondeg = [i for i, f in enumerate(F) if not degenerate(V, f)]
    bb = {}
    for i in nondeg:
        vs = [P[k] for k in F[i]]
        bb[i] = tuple(min(v[k] for v in vs) for k in range(3)) + \
                tuple(max(v[k] for v in vs) for k in range(3))
    t0 = time.process_time()
    pair = None
    tried = 0
    for ii in range(len(nondeg)):
        i = nondeg[ii]
        for jj in range(ii + 1, len(nondeg)):
            j = nondeg[jj]
            x, y = bb[i], bb[j]
            if (x[3] < y[0] or y[3] < x[0] or x[4] < y[1] or y[4] < x[1]
                    or x[5] < y[2] or y[5] < x[2]):
                continue
            tried += 1
            pair = (i, j)
            break
        if pair is not None:
            break
    if pair is None:
        log("[union] E4: AABB が交わる非退化の対がありません（**未測定**）")
    else:
        i, j = pair
        kind = classify_e4(P, F[i], F[j])
        dt = time.process_time() - t0
        log("[union] E4 の標本: 面 %d と %d → **%s**（AABB で外した対を含む探索 CPU %.4f 秒）"
            % (i, j, kind, dt))
        log("[union]   **1 対の値です。全対の保証にしません。**"
            "**「未判定」は正規扱いにしません。**")

    # ---- C2: ∪ の非退化な三角形 1 枚 ----
    t0 = time.process_time()
    tri = nondeg[0] if nondeg else None
    if tri is None:
        log("[union] C2: 非退化な三角形がありません（**未測定**）")
    else:
        res = c2_sample(P, F[tri], a)
        dt = time.process_time() - t0
        log("[union] C2 の標本: 面 %d → %s（CPU %.4f 秒）" % (tri, res, dt))
        log("[union]   **単一の入力三角形に収まらない場合は【未判定】**です"
            "（案 I は複数面の被覆を調べないので、「またがる」と「はみ出す」を区別できません）。")

    # ---- C1: 非退化な直線 1 本（∪ について）----
    t0 = time.process_time()
    if not e1_ok["union"]:
        log("[union] C1: **未評価**（E1 が破れているので、階段関数を巻き数と読めません）")
    else:
        got = c1_sample(P, F, a, log)
        log("[union] C1: %s（CPU %.4f 秒）" % (got, time.process_time() - t0))

    # ---- 単価: 符号つき 6 倍体積（∪ の 1 演算だけ）----
    V, F = data["union"]
    t0 = time.process_time()
    s = Fr(0)
    for f in F:
        s += det3(pt(V[f[0]]), pt(V[f[1]]), pt(V[f[2]]))
    dt = time.process_time() - t0
    saved = saved_row(a.saved, a.key)
    got = "%d/%d" % (s.numerator, s.denominator) if s.denominator != 1 else str(s.numerator)
    log("[union] 符号つき 6 倍体積: 分子 %d 桁 / 分母 %d 桁"
        % (len(str(abs(s.numerator))), len(str(s.denominator))))
    log("[union] 保存値（列 29）との一致: %s" % ("一致" if got == saved["vol_union"] else "★ 不一致"))
    log("[union] 体積の計算: 面 %d 枚で CPU %.3f 秒（1 枚あたり %.3g 秒）"
        % (len(F), dt, dt / max(len(F), 1)))
    log("**この単価は、この出力とこの演算についての観測値です。**"
        "**有理数の幅・約分・通る分岐に依存するので、他の演算や他の対の保証にしません。**")
    return 0


STAGES = {"identity": stage_identity, "anchor": stage_anchor, "measure": stage_measure}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True, choices=sorted(STAGES))
    p.add_argument("--out", required=True)
    p.add_argument("--in-a", default="data/logs/inspect/20260915-205108/C_A_quantized.bin")
    p.add_argument("--in-b", default="data/logs/inspect/20260915-205108/C_B_quantized.bin")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    path = os.path.join(a.out, "stage_%s.txt" % a.stage)
    f = open(path, "w")

    def log(s):
        print(s)
        sys.stdout.flush()
        f.write(s + "\n")
        f.flush()          # ★ 打ち切られても、そこまでが残るように

    log("=== 段 %s ===" % a.stage)
    t0 = time.monotonic()
    try:
        rc = STAGES[a.stage](a, log)
    finally:
        f.write("経過 %.2f 秒\n" % (time.monotonic() - t0))
        f.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
