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

# 保存行（列 29 の照合の基準）の sha256。**基準の側も照合します。**
EXPECT_ROW = {
    "250394x45413":
        "0bc986ec3709acd924d8a8d21dea70fae2681dfaa328330fac10d5bdebc781cd",
}
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


def volume6_rational(V, F):
    """符号つき 6 倍体積（同次座標から。**本体と合成対照が、この同じ関数を呼びます**）。"""
    s = Fr(0)
    for f in F:
        s += det3(pt(V[f[0]]), pt(V[f[1]]), pt(V[f[2]]))
    return s


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def e0_violations(V, F):
    """$E0$ の違反。戻り値: `(参照範囲, 3 頂点相異, 同次分母 0)` の一覧。

    **段 `measure` と合成対照が、この同じ実装を呼びます**
    （**写しを検査して本体が壊れたままになる形**を避けるため）。
    """
    oor = [i for i, f in enumerate(F) if max(f) >= len(V) or min(f) < 0]
    dup = [i for i, f in enumerate(F) if len(set(f)) < 3]
    w0 = [i for i, v in enumerate(V) if v[3] == 0]
    return oor, dup, w0


def e1_violations(F):
    r"""$E1$（各無向辺で $\#(u,v)=\#(v,u)$）の違反した有向辺。"""
    cnt = {}
    for f in F:
        for k in range(3):
            cnt[(f[k], f[(k + 1) % 3])] = cnt.get((f[k], f[(k + 1) % 3]), 0) + 1
    return [e for e in cnt if cnt[e] != cnt.get((e[1], e[0]), 0)]


def shared_by_coord(P, f1, f2):
    """**座標で**共有する点の一覧。**添字では判定しません**（§44.3）。"""
    out = []
    for i in f1:
        if any(P[i] == P[j] for j in f2) and P[i] not in out:
            out.append(P[i])
    return out


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
    # ★ 非共面。**`pair_a5` に丸投げしません** — あれは共有を【添字】で数えるので、
    #   併合漏れで同じ座標の別添字があると、幾何としては辺を共有する 2 枚が
    #   「共有しない」と扱われ、**実在し得る欠陥が「未判定」へ静かに流れます**（§44.3）。
    d1, d2 = ig.dot(n1, P[f1[0]]), ig.dot(n2, P[f2[0]])
    s1, s2 = ig._seg_on_line(P, f1, n2, d2), ig._seg_on_line(P, f2, n1, d1)
    if not s1 or not s2:
        return "正規（交わり空）"
    dv = ig.cross(n1, n2)
    ax = max(range(3), key=lambda k: abs(dv[k]))
    lo = max(min(p[ax] for p in s1), min(p[ax] for p in s2))
    hi = min(max(p[ax] for p in s1), max(p[ax] for p in s2))
    if lo > hi:
        return "正規（交わり空）"
    base = s1[0]

    def at(t):
        lam = Fr(t - base[ax], dv[ax])
        return tuple(base[k] + lam * dv[k] for k in range(3))

    pts = [at(lo)] if lo == hi else [at(lo), at(hi)]
    sh = shared_by_coord(P, f1, f2)
    if len(pts) == len(sh) and sorted(pts) == sorted(tuple(Fr(x) for x in q) for q in sh):
        return "正規（共有する単体と一致）"
    mid = at(Fr(lo + hi, 2)) if lo != hi else pts[0]

    def strict(T, n):
        k0 = max(range(3), key=lambda k: abs(n[k]))
        return ig._strict_inside2(mid, T, [k for k in range(3) if k != k0])

    if strict(T1, n1) and strict(T2, n2):
        return "不正（真の交差）"
    return "未判定（3 分類に当たりません）"


def _shared_or_undecided(P, f1, f2, poly):
    """交わりが、**座標で見た**共有頂点・辺と一致するか。**添字では判定しません。**"""
    sh = shared_by_coord(P, f1, f2)
    if len(poly) == len(sh) and sorted(poly) == sorted(tuple(Fr(x) for x in q) for q in sh):
        return "正規（共有する単体と一致）"
    return "未判定（3 分類に当たりません）"


def c2_sample(P, f, a):
    """支持平面が入力三角形のどれかと一致し、**単一の**入力三角形に載るか。"""
    import inspect_geom as ig
    import struct as _st
    n = ig.normal(P, f)
    d = ig.dot(n, P[f[0]])
    plane_hit = False
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
            plane_hit = True
            if all(ig._inside2(P[i], T, ij) for i in f):
                same = ig.dot(n, m) > 0
                return ("単一の入力三角形に収まる（%s の面 %d、向きは%s）"
                        % (name, gi, "同じ" if same else "★ 逆"))
    # ★ **「支持平面の不一致」と「収まらない」を分けます**（§44.8 の負の対照）。
    #   潰すと、**どの入力平面にも載っていない出力面**が違反 `C2` ではなく
    #   「未判定」として記録され、報告が通ってしまいます。
    if not plane_hit:
        return "**違反 C2**（どの入力三角形の支持平面とも一致しません）"
    return "未判定（支持平面は一致しますが、単一の入力三角形に収まりません）"


def compare_line(ex, ea, eb, op="union"):
    """直線上の階段関数を比べます。戻り値 `(C1-c の区間数, C1-b の区間数, 区間数)`。

    **段 `measure` と合成対照が、この同じ実装を呼びます。**
    """
    ts = sorted({t for t, _ in ex + ea + eb})
    if not ts:
        return None
    c = b = 0
    rows = []
    for k in range(len(ts) + 1):
        lo = "-inf" if k == 0 else _q(ts[k - 1])
        hi = "+inf" if k == len(ts) else _q(ts[k])
        m = (ts[0] - 1 if k == 0 else
             (ts[-1] + 1 if k == len(ts) else (ts[k - 1] + ts[k]) / 2))
        wx = sum(sg for t, sg in ex if t < m)
        wa = sum(sg for t, sg in ea if t < m)
        wb = sum(sg for t, sg in eb if t < m)
        if op == "union":
            want = 1 if (wa > 0 or wb > 0) else 0
        elif op == "isect":
            want = 1 if (wa > 0 and wb > 0) else 0
        elif op == "diff_ab":
            want = 1 if (wa > 0 and not wb > 0) else 0
        else:
            want = 1 if (wb > 0 and not wa > 0) else 0
        if (1 if wx > 0 else 0) != want:
            c += 1
            kind = "C1-c"
        elif wx != want:
            b += 1
            kind = "C1-b"
        else:
            kind = "ok"
        rows.append((lo, hi, _q(m), wx, wa, wb, want, kind))
    return c, b, len(ts) + 1, rows


def _q(x):
    """有理数を `分子/分母` で（分母 1 なら整数）。"""
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else "%d/%d" % (x.numerator, x.denominator)


def c1_sample(P, F, a, log, rows_out=None):
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
    skipped = {}
    c1_rows = []

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
            skipped["退化（出力）"] = skipped.get("退化（出力）", 0) + 1
            continue
        ea = events(ins["A"][0], ins["A"][1], o, d)
        eb = events(ins["B"][0], ins["B"][1], o, d)
        if ea is None or eb is None:
            skipped["退化（入力）"] = skipped.get("退化（入力）", 0) + 1
            continue
        # ★ **A と B の両方に交点を持つことを条件にします。**
        #   **1 回目の実行では、A を 1 度も通らない直線が 1 本目で採られました**
        #   （交点 X 2 / A 0 / B 4）。**最も見たい相互作用の領域を外します。**
        if not ea or not eb:
            skipped["A か B を外した"] = skipped.get("A か B を外した", 0) + 1
            continue
        r = compare_line(ev, ea, eb, "union")
        if r is None:
            skipped["交点なし"] = skipped.get("交点なし", 0) + 1
            continue
        bad_c, bad_b, nseg, rows = r
        c1_rows.extend(rows)
        log("[union] C1: 試行 %d 本目で採用（交点 X %d / A %d / B %d）"
            % (attempt, len(ev), len(ea), len(eb)))
        log("[union] C1: 採った直線 o=%s d=%s" % (tuple(str(x) for x in o), d))
        if rows_out is not None:
            rows_out.append(("line", _q(o[0]), _q(o[1]), _q(o[2]), d[0], d[1], d[2]))
            rows_out.extend(c1_rows)
        log("[union] C1: 区間 %d 個 / `C1-c` %d 個 / `C1-b` %d 個" % (nseg, bad_c, bad_b))
        log("[union] C1: 捨てた直線の内訳 %s" % (skipped if skipped else "無し"))
        log("[union]   **通った分岐**: 非退化の直線 1 本、交点 %d 個。"
            "**この直線が通らない食い違いは見えません。**" % len(ev))
        return "`C1-c` %d 個 / `C1-b` %d 個" % (bad_c, bad_b)
    log("[union] C1: 捨てた直線の内訳 %s" % (skipped if skipped else "無し"))
    return "**未測定**（32 回引いても条件を満たす直線が得られませんでした）"


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


def load_sums(path):
    """`SHA256SUMS` を `{ファイル名: ハッシュ}` に。**無ければ空**（照合できないので偽になります）。"""
    out = {}
    try:
        for line in open(path):
            t = line.split()
            if len(t) == 2:
                out[os.path.basename(t[1])] = t[0]
    except OSError:
        pass
    return out


def verify_inputs(a, log):
    """**幾何処理より前に、読むもの【全部】を照合します。** 1 つでも外れたら偽。"""
    ok = True
    for name, path in (("A", a.in_a), ("B", a.in_b)):
        with open(path, "rb") as f:
            h = hashlib.sha256(f.read()).hexdigest()
        good = (h == EXPECT_IN[name])
        log("照合 入力 %s: %s / %s" % (name, h, "一致" if good else "★ 不一致"))
        ok = ok and good
    src = getattr(a, "read_from", "") or a.out
    expect = load_sums(getattr(a, "sums", "") or os.path.join(src, os.pardir, "SHA256SUMS"))
    for fn in ["b_run_meta.txt"] + ["b_out_%s.bin" % o for o in OPS]:
        q = os.path.join(src, fn)
        with open(q, "rb") as f:
            h = hashlib.sha256(f.read()).hexdigest()
        want = expect.get(fn)
        if want is None:
            log("照合 %s: %s / **期待値が無いので照合できていません**" % (fn, h))
            ok = False
            continue
        good = (h == want)
        log("照合 %s: %s / %s" % (fn, h, "一致" if good else "★ 不一致"))
        ok = ok and good
    # ★ 保存物と保存行は、**計算して書くだけでは照合になりません**。期待値と比べます。
    #   列 29 は G1 の回帰（∪ の体積の一致）の【唯一の基準】なので、
    #   基準の側が無検証だと「一致」に意味がありません。
    with open(a.saved, "rb") as f:
        body = f.read()
    hs = hashlib.sha256(body).hexdigest()
    want_file = load_sums(getattr(a, "saved_sums", "")
                          or os.path.join(os.path.dirname(a.saved), "SHA256SUMS")
                          ).get(os.path.basename(a.saved))
    good = (want_file is not None) and (hs == want_file)
    log("照合 保存物: %s / 期待 %s / %s"
        % (hs, want_file or "**無し**", "一致" if good else "★ 照合できません"))
    ok = ok and good
    row = None
    for line in body.decode("utf-8", "replace").split("\n"):
        if line.split() and line.split()[0] == a.key:
            row = line
            break
    if row is None:
        log("★ 保存物に %s がありません。" % a.key)
        return False
    hr = hashlib.sha256(row.encode()).hexdigest()
    want_row = getattr(a, "expect_row", "") or EXPECT_ROW.get(a.key)
    good = (want_row is not None) and (hr == want_row)
    log("照合 保存行: %s / 期待 %s / %s"
        % (hr, want_row or "**無し**", "一致" if good else "★ 照合できません"))
    return ok and good


def stage_measure(a, log):
    src = getattr(a, "read_from", "") or a.out
    log("読み先 %s / 書き先 %s" % (src, a.out))
    if not verify_inputs(a, log):
        log("★ 入力の照合が通りませんでした。幾何処理へ進みません。")
        return 1
    meta = read_meta(os.path.join(src, "b_run_meta.txt"))
    nx, nw = int(meta["kHomoXyz"]), int(meta["kHomoW"])
    log("リム数: kHomoXyz %d / kHomoW %d" % (nx, nw))
    data = {}
    total_tri = 0
    # ---- E0（全数。停止対象）----
    stop = False
    for op in OPS:
        V, F, sha = read_soup(os.path.join(src, "b_out_%s.bin" % op), nx, nw)
        data[op] = (V, F)
        total_tri += len(F)
        oor, dup, w0 = e0_violations(V, F)
        log("[%s] 頂点 %d / 面 %d / sha256 %s" % (op, len(V), len(F), sha))
        log("[%s] E0 参照範囲 %s（違反 %d）/ 3 頂点相異 %s（違反 %d）/ 同次分母 0 %s（%d 個）"
            % (op, "通過" if not oor else "★ 破れ", len(oor),
               "通過" if not dup else "★ 破れ", len(dup),
               "無し" if not w0 else "★ あり", len(w0)))
        if oor or dup or w0:
            stop = True
    log("対象の合計 %d 三角形" % total_tri)
    log("**発火回数**: E0 の判定を %d 面 + %d 頂点ぶん行いました"
        % (total_tri, sum(len(data[o][0]) for o in OPS)))
    if stop:
        log("★ E0 が破れました。後続の計算が定義できません。ここで止めます。")
        return 1

    # ---- E1（全数。C1 の前提）----
    e1_ok = {}
    for op in OPS:
        V, F = data[op]
        bad = e1_violations(F)
        e1_ok[op] = not bad
        log("[%s] E1 各無向辺で #(u,v)=#(v,u): %s（違反 %d 本、有向辺 %d 本を数えました）"
            % (op, "通過" if not bad else "★ 破れ", len(bad), 3 * len(F)))
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
    c1_rows, c1_agg = [], None
    if not e1_ok["union"]:
        log("[union] C1: **未評価**（E1 が破れているので、階段関数を巻き数と読めません）")
        log("**分岐到達**: E1 違反 → C1 未評価 の分岐に 1 回到達しました")
    else:
        log("**分岐到達**: E1 通過 → C1 評価 の分岐に 1 回到達しました")
        got = c1_sample(P, F, a, log, rows_out=c1_rows)
        c1_agg = (sum(1 for r in c1_rows if r[0] != "line" and r[-1] == "C1-c"),
                  sum(1 for r in c1_rows if r[0] != "line" and r[-1] == "C1-b"))
        log("[union] C1: %s（CPU %.4f 秒）" % (got, time.process_time() - t0))

    # ---- 単価: 符号つき 6 倍体積（∪ の 1 演算だけ）----
    V, F = data["union"]
    t0 = time.process_time()
    s = volume6_rational(V, F)
    dt = time.process_time() - t0
    saved = saved_row(a.saved, a.key)
    got = "%d/%d" % (s.numerator, s.denominator) if s.denominator != 1 else str(s.numerator)
    log("[union] 符号つき 6 倍体積: 分子 %d 桁 / 分母 %d 桁"
        % (len(str(abs(s.numerator))), len(str(s.denominator))))
    vol_ok = (got == saved["vol_union"])
    log("[union] 保存値（列 29）との一致: %s" % ("一致" if vol_ok else "★ 不一致"))
    log("[union] 体積の計算: 面 %d 枚で CPU %.3f 秒（1 枚あたり %.3g 秒）"
        % (len(F), dt, dt / max(len(F), 1)))
    log("**この単価は、この出力とこの演算についての観測値です。**"
        "**有理数の幅・約分・通る分岐に依存するので、他の演算や他の対の保証にしません。**")

    # ---- G4: 区間の一覧を書き、**ファイルを読み直して**数え直す ----
    ipath = os.path.join(a.out, "stage_measure_c1_intervals.txt")
    with open(ipath, "w") as f:
        # ★ 必ず見出しを書きます（未評価・未測定でも空にしません）
        f.write("# line <o の x> <o の y> <o の z> <d の x> <d の y> <d の z>\n")
        f.write("# <下端> <上端> <代表点> <wX> <wA> <wB> <期待> <判定>\n")
        f.write("#   端点と代表点は有理数（分子/分母、分母 1 なら整数）。"
                "±無限は -inf / +inf\n")
        for r in c1_rows:
            f.write(" ".join(str(x) for x in r) + "\n")
    back_c = back_b = 0
    for line in open(ipath):
        t = line.split()
        if len(t) == 8 and t[0] not in ("line",):
            pass
        if len(t) >= 8 and t[-1] == "C1-c":
            back_c += 1
        elif len(t) >= 8 and t[-1] == "C1-b":
            back_b += 1
    log("[union] C1 の区間を %d 行書きました（見出しを除く）: %s" % (len(c1_rows), ipath))
    log("[union] ファイルを読み直して数え直し: C1-c %d / C1-b %d（段の集計 %s）"
        % (back_c, back_b, c1_agg))
    log("[union]   ★ これは**恒等式**です（同じ一覧から数えるため）。"
        "**独立な確認ではありません。確かめているのは書き出しの取りこぼしが無いことだけです。**")
    read_ok = (c1_agg is not None and (back_c, back_b) == c1_agg)
    has_line = any(r[0] == "line" for r in c1_rows)
    log("[union] G4: C1 の採用 %s / 区間の保存 %s / 読戻しの一致 %s"
        % ("あり" if has_line else "**無し**", "あり" if len(c1_rows) > 1 else "**無し**",
           "一致" if read_ok else "★ 不一致"))
    g4_ok = has_line and len(c1_rows) > 1 and read_ok
    if not g4_ok:
        log("[union] ★ **G4 未解消**（C1 が未評価 / 未測定、または書き出しの取りこぼし）")
    if not vol_ok:
        log("[union] ★ **G1 の回帰が破れた**（$\\cup$ の体積が列 29 と一致しません）")
    log("判定: %s" % ("G4 解消・体積の回帰も一致" if (g4_ok and vol_ok) else "★ 未解消あり"))
    return 0 if (g4_ok and vol_ok) else 1


STAGES = {"identity": stage_identity, "anchor": stage_anchor, "measure": stage_measure}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True, choices=sorted(STAGES))
    p.add_argument("--out", required=True)
    p.add_argument("--read-from", default="", help="保存束を【読むだけ】の場所。--out とは分けます")
    p.add_argument("--sums", default="", help="読み先の SHA256SUMS（既定は読み先の親）")
    p.add_argument("--saved-sums", default="", help="保存物の SHA256SUMS（既定は保存物と同じ場所）")
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
