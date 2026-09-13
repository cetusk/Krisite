# -*- coding: utf-8 -*-
"""出力の曲面が、問題の頂点の近くで【埋め込まれているか】を厳密に確かめる。

**C++ の実装とは別経路です**（Python の有理数演算。`CLAUDE.md`「正解器は被検体と
別経路で書く。言語ごと変えられればなお良い」）。

リンクが埋め込まれていないなら、交わりは頂点 v のいくらでも近くで起きる。
したがって **v に接する三角形どうし**を調べれば足りる。

判定: 頂点 v を共有する 2 枚 T1, T2 について、
  「T1 の相対内部と T2 が交わるか」を、v から見た【立体角】の重なりで見る。
  v を原点に取り、各三角形を v から出る 2 本の半直線が張る扇形（2 次元の角）とみなす。
  球面上の 2 つの弧が交差すれば、2 枚は v の近くで貫き合っている。
"""
import sys, re
from fractions import Fraction


def parse(path):
    stars = []          # [(v, [tri...], {vid: (x,y,z,w)})]
    cur = None
    for ln in open(path, encoding='utf-8'):
        m = re.match(r'\s*STAR (\d+)', ln)
        if m:
            cur = (int(m.group(1)), [], {})
            stars.append(cur)
            continue
        m = re.match(r'\s*TRI (\d+) (\d+) (\d+) (\d+)', ln)
        if m and cur is not None:
            cur[1].append(tuple(int(x) for x in m.groups()))
            continue
        m = re.match(r'\s*VTX (\d+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+)', ln)
        if m and cur is not None:
            vid = int(m.group(1))
            vals = []
            for g in m.groups()[1:]:
                v = int(g[1:], 16)
                vals.append(-v if g[0] == '-' else v)
            cur[2][vid] = tuple(vals)
    return stars


def euclid(h):
    x, y, z, w = h
    return (Fraction(x, w), Fraction(y, w), Fraction(z, w))


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sgn(x):
    return (x > 0) - (x < 0)


def _inside(p, a, b, n):
    return sgn(dot(cross(a, p), n)) > 0 and sgn(dot(cross(p, b), n)) > 0


def _wedges_overlap(u1, u2, w1, w2, n):
    """同一平面の 2 つの扇が、角度として内部で重なるか。"""
    # n の向きを扇 A に合わせる（扇 B の向きは逆かもしれない）
    nB = cross(w1, w2)
    a1, a2 = (w1, w2) if sgn(dot(nB, n)) > 0 else (w2, w1)
    for p in (u1, u2):
        if _inside(p, a1, a2, n):
            return True
    for p in (a1, a2):
        if _inside(p, u1, u2, n):
            return True
    return False


def arcs_cross(u1, u2, w1, w2):
    """v を原点として、扇 (u1,u2) と扇 (w1,w2) が【内部で】重なるか。

    球面上の弧どうしの交差判定に落とす。弧 A = u1→u2（短いほうではなく、
    扇が張る側）、弧 B = w1→w2。
    2 つの平面 n_A = u1×u2、n_B = w1×w2 の交線 d = n_A×n_B を求め、
    ±d のどちらかが両方の扇の【内部】にあれば交差する。
    """
    nA = cross(u1, u2)
    nB = cross(w1, w2)
    if all(c == 0 for c in nA) or all(c == 0 for c in nB):
        return False           # 退化した扇（v で潰れている）
    d = cross(nA, nB)
    if all(c == 0 for c in d):
        # **同一平面。** 扇が角度として重なっていれば、2 枚は面で重なっている
        # （曲面が埋め込まれていない別の形）。**貫きとは区別して数えます。**
        return "coplanar" if _wedges_overlap(u1, u2, w1, w2, nA) else False
    for s in (1, -1):
        p = (s * d[0], s * d[1], s * d[2])
        # p が扇 (u1,u2) の内部にあるか: u1×p と p×u2 が n_A と同じ向き
        if sgn(dot(cross(u1, p), nA)) <= 0:
            continue
        if sgn(dot(cross(p, u2), nA)) <= 0:
            continue
        if sgn(dot(cross(w1, p), nB)) <= 0:
            continue
        if sgn(dot(cross(p, w2), nB)) <= 0:
            continue
        return True
    return False


def main(path):
    stars = parse(path)
    print(f"読んだ STAR {len(stars)} 個")
    tot_cross = 0
    tot_cop = 0
    per = []
    for v, tris, vtx in stars:
        o = euclid(vtx[v])
        fans = []
        for t in tris:
            ti, a, b, c = t
            others = [x for x in (a, b, c) if x != v]
            if len(others) != 2:
                continue
            u1 = sub(euclid(vtx[others[0]]), o)
            u2 = sub(euclid(vtx[others[1]]), o)
            fans.append((ti, others[0], others[1], u1, u2))
        n = 0
        m = 0
        for i in range(len(fans)):
            for j in range(i + 1, len(fans)):
                ti, p1, p2, u1, u2 = fans[i]
                tj, q1, q2, w1, w2 = fans[j]
                if len({p1, p2} & {q1, q2}) > 0:
                    continue            # 辺を共有する隣どうしは対象外
                r = arcs_cross(u1, u2, w1, w2)
                if r == "coplanar":
                    m += 1
                    if m <= 3:
                        print(f"  v{v}: **t{ti} と t{tj} が同一平面で重なっている**")
                elif r:
                    n += 1
                    if n <= 3:
                        print(f"  v{v}: **t{ti} と t{tj} が v の近くで貫き合っている**")
        per.append((v, len(fans), n, m))
        tot_cross += n
        tot_cop += m
    print()
    print("| 頂点 | 接する三角形 | **貫き合う対** | **同一平面で重なる対** |")
    print("|---|---:|---:|---:|")
    for v, k, n, m in per:
        print(f"| {v} | {k} | **{n}** | **{m}** |")
    print()
    print(f"**貫き 合計 {tot_cross} 対 / 同一平面の重なり 合計 {tot_cop} 対 / STAR {len(stars)} 個**")
    ok = sum(1 for _, _, n, m in per if n > 0)
    ok2 = sum(1 for _, _, n, m in per if m > 0)
    print(f"**貫きがあった頂点 {ok} / {len(per)}**")
    print(f"**同一平面の重なりがあった頂点 {ok2} / {len(per)}**")


main(sys.argv[1])
