#!/usr/bin/env python3
"""合成入力のための幾何検査の部品（厳密。浮動小数点を使いません）。

**この節の座標・面・方向は、このファイルが原本です。** 文書と二重に持ちません。

**被検体（`include/krisite/`）を 1 行も呼びません。**
**`tests/tools/design_sweep.py` の値を正解として使いません**（同ファイルは凍結）。

扱う対象: 整数座標の三角形メッシュ（頂点配列 P と、頂点の添字 3 つ組の配列 F）。
"""
from fractions import Fraction as Fr
from collections import Counter, defaultdict

# レイの方向。**このファイルが原本です。** 相異なり、零ベクトルを含みません。
#
# ★ 【正負の両方】を入れます。非負の八分区間だけだと、その反対側にある欠損
#   （開いた曲面の穴など）を、どの方向からも見られません。
#   実測: 面を 1 枚除いた立方体の内側で、非負 16 本では全方向が一致して 1 を返し、
#   「測れない」と言えませんでした。正負 32 本にすると食い違いを検出します。
_DIRS_POS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
             (1, 0, 1), (0, 1, 1), (1, 1, 1), (2, 1, 0),
             (1, 2, 0), (2, 0, 1), (1, 0, 2), (0, 2, 1),
             (0, 1, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)]
DIRS = _DIRS_POS + [tuple(-x for x in d) for d in _DIRS_POS]


# --- ベクトル ------------------------------------------------------------
def sub(a, b): return tuple(a[k] - b[k] for k in range(3))
def dot(a, b): return sum(a[k] * b[k] for k in range(3))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def normal(P, f):
    """面の法線（正規化しない整数ベクトル）。共線面では零ベクトル。"""
    return cross(sub(P[f[1]], P[f[0]]), sub(P[f[2]], P[f[0]]))


def signed_volume6(P, F):
    """符号つき 6 倍体積。外向きの閉曲面なら正。"""
    s = 0
    for f in F:
        a, b, c = P[f[0]], P[f[1]], P[f[2]]
        s += (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
              + a[2]*(b[0]*c[1]-b[1]*c[0]))
    return s


# --- 立方体の生成（箱と向きから） ---------------------------------------
def box(lo, hi, sign=+1):
    """軸平行な立方体。sign=+1 で外向き、-1 で内向き。

    頂点は x→y→z の順に lo/hi を回した辞書順（8 個）。
    面は軸 x→y→z、各軸で lo 側→hi 側、各面を 2 枚に割る（12 枚）。
    割り方は「その面の法線の先端から見返して反時計回りに q0..q3、
    q0 は辞書順最小、対角線 q0q2」。
    """
    V = [(x, y, z) for x in (lo[0], hi[0])
         for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    idx = {v: i for i, v in enumerate(V)}
    F = []
    for ax in range(3):
        for side, out in ((lo[ax], -1), (hi[ax], +1)):
            nrm = out * sign
            corners = sorted(v for v in V if v[ax] == side)
            c0 = corners[0]
            rest = [v for v in corners if v != c0]
            # c0 の隣（1 成分だけ違う）2 つと、対角の 1 つ
            adj = [v for v in rest if sum(1 for k in range(3) if sub(v, c0)[k]) == 1]
            far = [v for v in rest if v not in adj][0]
            q1, q2, q3 = adj[0], far, adj[1]
            nvec = tuple(nrm if k == ax else 0 for k in range(3))
            if dot(cross(sub(q1, c0), sub(q2, c0)), nvec) < 0:
                q1, q3 = q3, q1                    # 法線の先端から見て反時計回りに
            F += [(idx[c0], idx[q1], idx[q2]), (idx[c0], idx[q2], idx[q3])]
    return V, F


def join(*meshes):
    """複数のメッシュを 1 本の頂点配列にまとめる（添字は共有しません）。"""
    P, F = [], []
    for v, f in meshes:
        off = len(P)
        P += list(v)
        F += [tuple(x + off for x in t) for t in f]
    return P, F


# --- A0〜A4 --------------------------------------------------------------
def check_a0(P, F):
    """戻り値: (範囲外の面, 3 頂点が相異でない面, 共線面の一覧)"""
    oor = [f for f in F if max(f) >= len(P) or min(f) < 0]
    dup = [f for f in F if len(set(f)) < 3]
    col = [] if (oor or dup) else [f for f in F if normal(P, f) == (0, 0, 0)]
    return oor, dup, col


def check_a1(F):
    """各無向辺に接する面がちょうど 2 枚か。戻り値: 違反した辺の一覧"""
    e = Counter()
    for f in F:
        for k in range(3):
            e[frozenset((f[k], f[(k+1) % 3]))] += 1
    return sorted((tuple(sorted(k)), v) for k, v in e.items() if v != 2)


def check_a2(F):
    """各有向辺がちょうど 1 回か。戻り値: 違反した有向辺の一覧"""
    d = Counter()
    for f in F:
        for k in range(3):
            d[(f[k], f[(k+1) % 3])] += 1
    return sorted((k, v) for k, v in d.items() if v != 1)


def link_cycles(F, v):
    """頂点 v のリンクが何本の閉路か。閉路でなければ None。"""
    es = [tuple(x for x in f if x != v) for f in F if v in f]
    adj = defaultdict(list)
    for a, b in es:
        adj[a].append(b); adj[b].append(a)
    if not adj or any(len(t) != 2 for t in adj.values()):
        return None
    seen, n = set(), 0
    for s in adj:
        if s in seen:
            continue
        n += 1
        st = [s]; seen.add(s)
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
    return n


def check_a3(F):
    """各頂点のリンクが単一の閉路か。戻り値: (頂点, 閉路数 or None) の一覧

    **前提: A0-2 を先に通すこと**（面に同じ頂点が 2 回あると例外になります）。
    """
    bad = []
    for v in sorted({x for f in F for x in f}):
        n = link_cycles(F, v)
        if n != 1:
            bad.append((v, n))
    return bad


def check_a4(F):
    """重複した面（同じ 3 頂点の組）。戻り値: 重複した 3 つ組の一覧"""
    c = Counter(tuple(sorted(f)) for f in F)
    return sorted(k for k, v in c.items() if v > 1)


# --- A5（面どうしの交わり） ---------------------------------------------
def _seg_on_line(P, f, n2, d2):
    """三角形 f を平面 (n2,d2) で切った線分を、[点, 点] で返す（無ければ None）。"""
    vs = [P[i] for i in f]
    ds = [dot(n2, v) - d2 for v in vs]
    if all(x > 0 for x in ds) or all(x < 0 for x in ds):
        return None
    pts = []
    for k in range(3):
        a, b = vs[k], vs[(k+1) % 3]
        da, db = ds[k], ds[(k+1) % 3]
        if da == 0:
            pts.append(tuple(Fr(x) for x in a))
        if (da > 0) != (db > 0) and da != 0 and db != 0:
            t = Fr(-da, db - da)
            pts.append(tuple(Fr(a[i]) + t * (b[i] - a[i]) for i in range(3)))
    uniq = []
    for p in pts:
        if p not in uniq:
            uniq.append(p)
    return uniq if uniq else None


def _inside2(pt, tri, ij):
    """2 次元（軸 ij）で、点が三角形の内部か境界にあるか。"""
    sg = 0
    for k in range(3):
        a, b = tri[k], tri[(k+1) % 3]
        e = (b[ij[0]] - a[ij[0]], b[ij[1]] - a[ij[1]])
        w = (pt[ij[0]] - a[ij[0]], pt[ij[1]] - a[ij[1]])
        cr = e[0]*w[1] - e[1]*w[0]
        if cr != 0:
            if sg and (cr > 0) != (sg > 0):
                return False
            sg = cr
    return True


def _strict_inside2(pt, tri, ij):
    sg = 0
    for k in range(3):
        a, b = tri[k], tri[(k+1) % 3]
        e = (b[ij[0]] - a[ij[0]], b[ij[1]] - a[ij[1]])
        w = (pt[ij[0]] - a[ij[0]], pt[ij[1]] - a[ij[1]])
        cr = e[0]*w[1] - e[1]*w[0]
        if cr == 0:
            return False
        if sg and (cr > 0) != (sg > 0):
            return False
        sg = cr
    return True


def _clip2(poly, a, b, ij):
    """2 次元で、有向辺 a→b の左側に凸多角形を切る（Sutherland-Hodgman）。"""
    def side(p):
        return ((b[ij[0]]-a[ij[0]])*(p[ij[1]]-a[ij[1]])
                - (b[ij[1]]-a[ij[1]])*(p[ij[0]]-a[ij[0]]))
    out = []
    for k in range(len(poly)):
        p, q = poly[k], poly[(k+1) % len(poly)]
        sp, sq = side(p), side(q)
        if sp >= 0:
            out.append(p)
        if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
            t = Fr(sp, sp - sq)
            out.append(tuple(Fr(p[m]) + t * (Fr(q[m]) - Fr(p[m])) for m in range(3)))
    return out


def _inter2(T1, T2, ij):
    """共面な 2 つの三角形の交わり（2 次元）。頂点の一覧を返す（空なら []）。"""
    poly = [tuple(Fr(x) for x in v) for v in T1]
    # T2 の向きを、2 次元で反時計回りに揃えてから切る
    ar = sum((T2[k][ij[0]]*T2[(k+1) % 3][ij[1]] - T2[(k+1) % 3][ij[0]]*T2[k][ij[1]])
             for k in range(3))
    T2o = T2 if ar > 0 else [T2[0], T2[2], T2[1]]
    for k in range(3):
        if not poly:
            break
        poly = _clip2(poly, T2o[k], T2o[(k+1) % 3], ij)
    uniq = []
    for p in poly:
        if p not in uniq:
            uniq.append(p)
    return uniq


def _area2(poly, ij):
    if len(poly) < 3:
        return Fr(0)
    a = Fr(0)
    for k in range(len(poly)):
        p, q = poly[k], poly[(k+1) % len(poly)]
        a += p[ij[0]]*q[ij[1]] - q[ij[0]]*p[ij[1]]
    return a


def pair_a5(P, f1, f2):
    """1 対の判定。戻り値: (共有頂点数 s, 正規か, 分類)

    §43.4 の規則:
      s=0 … 交わりが【空】なら正規
      s=1 … 交わりが【その頂点だけ】なら正規
      s=2 … 交わりが【その辺だけ】なら正規
      s=3 … 同じ 3 頂点の組。A4 が先に数えるが、ここでも不正（接触）として扱う
    分類は "交差"（共面でなく、交わりが両方の相対内部と交わる）/ "接触"（それ以外の不正）。
    """
    s = len(set(f1) & set(f2))
    n1, n2 = normal(P, f1), normal(P, f2)
    if n1 == (0, 0, 0) or n2 == (0, 0, 0):
        return s, True, None                      # 共線面は A0-3 が数える
    d1, d2 = dot(n1, P[f1[0]]), dot(n2, P[f2[0]])
    T1 = [P[i] for i in f1]
    T2 = [P[i] for i in f2]
    shared = sorted(set(f1) & set(f2))
    cop = cross(n1, n2) == (0, 0, 0) and dot(n1, P[f2[0]]) == d1

    def ok_for_shared(pts):
        """交わり pts が、共有する単体（点・辺）とちょうど一致するか。"""
        want = [tuple(Fr(x) for x in P[i]) for i in shared]
        if s == 1:
            return len(pts) == 1 and pts[0] == want[0]
        if s == 2:
            if len(pts) != 2:
                return False
            return sorted(pts) == sorted(want)
        return False

    if cop:
        ax = max(range(3), key=lambda k: abs(n1[k]))
        ij = [k for k in range(3) if k != ax]
        poly = _inter2(T1, T2, ij)
        if not poly:
            return s, True, None                  # 交わらない
        if _area2(poly, ij) != 0:
            return s, False, "接触"               # 共面なので交差ではない
        if ok_for_shared(poly):
            return s, True, None
        return s, False, "接触"                   # 退化した交わりだが空ではない
    # 共面でない: 交線上の区間どうしの重なり
    s1, s2 = _seg_on_line(P, f1, n2, d2), _seg_on_line(P, f2, n1, d1)
    if not s1 or not s2:
        return s, True, None
    dv = cross(n1, n2)
    ax = max(range(3), key=lambda k: abs(dv[k]))
    lo = max(min(p[ax] for p in s1), min(p[ax] for p in s2))
    hi = min(max(p[ax] for p in s1), max(p[ax] for p in s2))
    if lo > hi:
        return s, True, None                      # 交わらない
    base = s1[0]

    def at(t):
        lam = Fr(t - base[ax], dv[ax])
        return tuple(base[k] + lam * dv[k] for k in range(3))
    pts = [at(lo)] if lo == hi else [at(lo), at(hi)]
    if ok_for_shared(pts):
        return s, True, None
    mid = at((lo + hi) / 2)

    def strict(T, n):
        a = max(range(3), key=lambda k: abs(n[k]))
        return _strict_inside2(mid, T, [k for k in range(3) if k != a])
    return s, False, ("交差" if (strict(T1, n1) and strict(T2, n2)) else "接触")


def check_a5(P, F):
    """全対を見る。戻り値: (交差の件数, 接触の件数, 不正な対の一覧)"""
    nx = nt = 0
    bad = []
    for i in range(len(F)):
        for j in range(i + 1, len(F)):
            s, ok, kind = pair_a5(P, F[i], F[j])
            if not ok:
                bad.append((i, j, s, kind))
                if kind == "交差":
                    nx += 1
                else:
                    nt += 1
    return nx, nt, bad


# --- 成分・巻き数・領域 --------------------------------------------------
def components(F):
    """無向辺を共有する面どうしの連結成分。戻り値: 面の添字の集合の一覧"""
    e2f = defaultdict(list)
    for i, f in enumerate(F):
        for k in range(3):
            e2f[frozenset((f[k], f[(k+1) % 3]))].append(i)
    seen, out = set(), []
    for i in range(len(F)):
        if i in seen:
            continue
        st, grp = [i], set()
        seen.add(i)
        while st:
            x = st.pop(); grp.add(x)
            for k in range(3):
                for y in e2f[frozenset((F[x][k], F[x][(k+1) % 3]))]:
                    if y not in seen:
                        seen.add(y); st.append(y)
        out.append(sorted(grp))
    return out


def winding(P, F, q, faces=None):
    """点 q（有理座標）での巻き数。DIRS を順に試し、退化しない方向を使う。

    戻り値: (巻き数, 使った方向の番号, 有効だった方向の本数)
    有効な方向が 2 本未満、または値が食い違うときは (None, None, 本数)。
    """
    idxs = range(len(F)) if faces is None else faces
    vals = []
    for di, d in enumerate(DIRS):
        w, ok = 0, True
        for i in idxs:
            f = F[i]
            n = normal(P, f)
            if n == (0, 0, 0):
                continue
            den = dot(n, d)
            num = dot(n, sub(P[f[0]], tuple(q)))
            T0 = [P[i2] for i2 in f]
            ax0 = max(range(3), key=lambda k: abs(n[k]))
            ij0 = [k for k in range(3) if k != ax0]
            if num == 0 and _inside2(tuple(q), T0, ij0):
                ok = False; break                  # 問い合わせ点が面の上にある
            if den == 0:
                if num == 0:
                    ok = False; break              # 面の平面上を走る
                continue
            t = Fr(num, den)
            if t <= 0:
                continue
            x = tuple(q[k] + t * d[k] for k in range(3))
            T = [P[i2] for i2 in f]
            ax = max(range(3), key=lambda k: abs(n[k]))
            ij = [k for k in range(3) if k != ax]
            if _strict_inside2(x, T, ij):
                w += 1 if den > 0 else -1
            elif _inside2(x, T, ij):
                ok = False; break                  # 辺・頂点に当たった
        if ok:
            vals.append((di, w))
    if len(vals) < 2:
        return None, None, len(vals)
    if len({w for _, w in vals}) != 1:
        return None, None, len(vals)
    return vals[0][1], vals[0][0], len(vals)


def face_centroid(P, f):
    return tuple(Fr(P[f[0]][k] + P[f[1]][k] + P[f[2]][k], 3) for k in range(3))


def analyze(P, F):
    """成分に分解し、包含の森・領域ごとの巻き数と体積・残差を求める。

    A5 を通った（成分どうしが交差も接触もしない）ことを前提にします。
    戻り値は辞書。判定はしません（呼び出し側が行います）。
    """
    comps = components(F)
    S = [signed_volume6(P, [F[i] for i in c]) for c in comps]
    SB = sum(S)
    # 包含: 成分 c の 1 点が、成分 j の内側か
    probe = [face_centroid(P, F[c[0]]) for c in comps]
    inside = [[] for _ in comps]
    for c in range(len(comps)):
        for j in range(len(comps)):
            if c == j:
                continue
            w, _, _ = winding(P, F, probe[c], faces=comps[j])
            if w is None:
                return {"error": "巻き数を測れません（成分 %d を成分 %d から）" % (c, j)}
            if w != 0:
                inside[c].append(j)
    # 親 = c を含むもののうち、最も内側のもの
    #     ＝ 他のすべての包含者に【含まれる】もの。
    #     「j が k に含まれる」は k in inside[j]（inside[x] は x を含む成分の一覧）。
    parent = []
    for c in range(len(comps)):
        cand = inside[c]
        p = None
        for j in cand:
            if all(k == j or k in inside[j] for k in cand):
                p = j; break
        parent.append(p)
    child = defaultdict(list)
    for c, p in enumerate(parent):
        if p is not None:
            child[p].append(c)
    # 領域 j = 成分 j の内側かつ子の外側
    vol = [abs(S[j]) - sum(abs(S[c]) for c in child[j]) for j in range(len(comps))]
    w = []
    for j in range(len(comps)):
        chain, x = [j], parent[j]
        while x is not None:
            chain.append(x); x = parent[x]
        w.append(sum(1 if S[m] > 0 else -1 for m in chain))
    mu = sum(v for v, x in zip(vol, w) if x > 0)
    inner = [j for j in range(len(comps)) if w[j] >= 2]
    return {"components": comps, "S": S, "S_total": SB, "parent": parent,
            "region_vol": vol, "region_w": w, "measure": mu, "R": mu - SB,
            "inner": inner, "S_inner": sum(S[j] for j in inner)}
