# -*- coding: utf-8 -*-
"""**radial_pair が作った組が、本当に「角度順で隣どうし」か**を厳密に確かめる。

**C++ とは別経路**（Python の有理数演算）。
`DESIGN-phase5-vertex-level.md` §3 の前提 (d) を直接検査する。

手順:
  1. 辺 (a,b) の 4 枚について、辺に垂直で三角形の内側を向く方向 m_i を厳密に作る
  2. 辺の方向 d のまわりで m_i を角度順に並べる（符号だけで比較）
  3. 角度順で隣どうしの組（2 通りある）と、radial_pair が作った組を比べる
"""
import sys, re
from fractions import Fraction


def sgn(x):
    return (x > 0) - (x < 0)


def sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(a[i]*b[i] for i in range(3))


def parse(path):
    """(辺, 組, STAR ごとの三角形と座標) を取り出す。"""
    edges = []          # (a, b, [(g0a,g0b),(g1a,g1b)])
    stars = {}          # v -> ([tri], {vid: h})
    cur = None
    for ln in open(path, encoding='utf-8'):
        m = re.match(r'\s*分裂前 (\d+)-(\d+) 次数 (\d+) .*組 (\d+): \(t(\d+),t(\d+)\)\(t(\d+),t(\d+)\)', ln)
        if m:
            g = [int(x) for x in m.groups()]
            edges.append((g[0], g[1], [(g[4], g[5]), (g[6], g[7])]))
            continue
        m = re.match(r'\s*STAR (\d+)', ln)
        if m:
            cur = stars.setdefault(int(m.group(1)), ([], {}))
            continue
        m = re.match(r'\s*TRI (\d+) (\d+) (\d+) (\d+)', ln)
        if m and cur is not None:
            t = tuple(int(x) for x in m.groups())
            if t not in cur[0]:
                cur[0].append(t)
            continue
        m = re.match(r'\s*VTX (\d+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+)', ln)
        if m and cur is not None:
            vals = []
            for g in m.groups()[1:]:
                v = int(g[1:], 16)
                vals.append(-v if g[0] == '-' else v)
            cur[1][int(m.group(1))] = tuple(vals)
    return edges, stars


def eu(h):
    x, y, z, w = h
    return (Fraction(x, w), Fraction(y, w), Fraction(z, w))


def main(path):
    edges, stars = parse(path)
    print(f"辺 {len(edges)} 本 / STAR {len(stars)} 個\n")
    agree = disagree = skipped = 0
    shown = 0
    for a, b, groups in edges:
        if a not in stars:
            skipped += 1
            continue
        tris, vtx = stars[a]
        if b not in vtx:
            skipped += 1
            continue
        A, B = eu(vtx[a]), eu(vtx[b])
        d = sub(B, A)
        four = []
        for (ti, p, q, r) in tris:
            vs = (p, q, r)
            if a not in vs or b not in vs:
                continue
            c = [x for x in vs if x != a and x != b]
            if len(c) != 1 or c[0] not in vtx:
                continue
            C = eu(vtx[c[0]])
            ac = sub(C, A)
            # 辺に垂直な成分（有理数で厳密）
            k = Fraction(dot(ac, d), dot(d, d))
            m = tuple(ac[i] - k * d[i] for i in range(3))
            four.append((ti, m))
        if len(four) != 4:
            skipped += 1
            continue
        # d まわりの角度で並べる（基準 = four[0]）
        base = four[0][1]
        def key(item):
            m = item[1]
            s = sgn(dot(cross(base, m), d))
            if s > 0:
                cls = 1
            elif s < 0:
                cls = 3
            else:
                cls = 0 if sgn(dot(base, m)) > 0 else 2
            return cls
        order = sorted(four, key=key)
        # 同じ類の中は外積の符号で並べ直す（4 枚なので挿入ソートで足りる）
        for i in range(len(order)):
            for j in range(i + 1, len(order)):
                if key(order[i]) == key(order[j]):
                    if sgn(dot(cross(order[i][1], order[j][1]), d)) < 0:
                        order[i], order[j] = order[j], order[i]
        ids = [t for t, _ in order]
        adj1 = {frozenset((ids[0], ids[1])), frozenset((ids[2], ids[3]))}
        adj2 = {frozenset((ids[1], ids[2])), frozenset((ids[3], ids[0]))}
        got = {frozenset(g) for g in groups}
        ok = (got == adj1) or (got == adj2)
        if ok:
            agree += 1
        else:
            disagree += 1
            if shown < 5:
                shown += 1
                print(f"  **辺 {a}-{b}: 組が角度順で隣どうしではありません**")
                print(f"    角度順 {ids}")
                print(f"    組     {[sorted(g) for g in groups]}")
    print()
    print("| 判定 | 辺の本数 |")
    print("|---|---:|")
    print(f"| **組が角度順で隣どうし**（前提 (d) が成立） | **{agree}** |")
    print(f"| **隣どうしでない**（前提 (d) が破れる） | **{disagree}** |")
    print(f"| 調べられなかった（STAR に無い） | {skipped} |")


main(sys.argv[1])
