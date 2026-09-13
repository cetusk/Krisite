# -*- coding: utf-8 -*-
"""**平面の和で作った点が、本当に辺の内部に載るか**を厳密に確かめる。

**C++ とは別経路**（Python の有理数演算。`CLAUDE.md`「正解器は被検体と別経路で書く。
言語ごと変えられればなお良い」）。

案 H は辺 $(v,w)$ を細分します。**点どうしの中点は幅が倍になる**ので
（EMBER §3.2 が「固定精度では一般に計算できない」と書いた操作そのもの）、
**平面による構成**に置き換えられるかを調べます。

  1. 両端に載る平面 2 枚 $P_1, P_2$ で、辺の線 $L$ を張る
  2. $v$ だけに載る $P_3$ と、$w$ だけに載る $P_4$ を取る
  3. $Q = s_3 P_3 + s_4 P_4$ の符号を、$Q(v)$ と $Q(w)$ が逆になるように選ぶ
  4. $m = P_1 \\cap P_2 \\cap Q$ が、辺の内部にあるか

そして**係数のビット幅を実測**し、導出した上界と突き合わせます。
"""
import sys, re
from fractions import Fraction


def rd(tok):
    v = int(tok[1:], 16)
    return -v if tok[0] == '-' else v


def parse(path):
    recs = []
    cur = None
    for ln in open(path, encoding='utf-8'):
        m = re.match(r'\s*EDGEPT (\d+) (\d+)', ln)
        if m:
            cur = {'a': int(m.group(1)), 'b': int(m.group(2)), 'pts': {}, 'planes': []}
            recs.append(cur)
            continue
        m = re.match(r'\s*EP ([ab]) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+)', ln)
        if m and cur is not None:
            cur['pts'][m.group(1)] = tuple(rd(g) for g in m.groups()[1:])
            continue
        m = re.match(r'\s*PL (\d+) (-?\d+) (-?\d+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+) ([+-][0-9a-f]+)', ln)
        if m and cur is not None:
            cur['planes'].append({
                'id': int(m.group(1)), 'sa': int(m.group(2)), 'sb': int(m.group(3)),
                'co': tuple(rd(g) for g in m.groups()[3:]),
            })
    return recs


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def meet3(p, q, r):
    """3 平面の交点（同次座標）。w = det(N), x = -det(d を x 列に入れたもの) …"""
    N = [p[:3], q[:3], r[:3]]
    d = [p[3], q[3], r[3]]
    w = det3(N)
    if w == 0:
        return None
    xs = []
    for col in range(3):
        M = [list(N[i]) for i in range(3)]
        for i in range(3):
            M[i][col] = d[i]
        xs.append(-det3(M))
    return (xs[0], xs[1], xs[2], w)


def ev(pl, h):
    """平面 pl の同次点 h での値（符号は w を掛けて正規化）。"""
    v = pl[0] * h[0] + pl[1] * h[1] + pl[2] * h[2] + pl[3] * h[3]
    return v if h[3] > 0 else -v


def bits(n):
    return n.bit_length() + 1  # 符号ビットを足す


def eu(h):
    return tuple(Fraction(h[i], h[3]) for i in range(3))


def main(path):
    recs = parse(path)
    print(f"辺 {len(recs)} 本\n")
    ok = 0
    fail = []
    maxbits = {'w': 0, 'xyz': 0, 'Qn': 0, 'Qd': 0}
    for r in recs:
        A, B = r['pts']['a'], r['pts']['b']
        onboth = [p['co'] for p in r['planes'] if p['sa'] == 0 and p['sb'] == 0]
        onlya = [p['co'] for p in r['planes'] if p['sa'] == 0 and p['sb'] != 0]
        onlyb = [p['co'] for p in r['planes'] if p['sb'] == 0 and p['sa'] != 0]
        # 平行でない 2 枚を選ぶ
        P1 = P2 = None
        for i in range(len(onboth)):
            for j in range(i + 1, len(onboth)):
                u, v = onboth[i][:3], onboth[j][:3]
                cx = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
                if any(cx):
                    P1, P2 = onboth[i], onboth[j]
                    break
            if P1 is not None:
                break
        if P1 is None or not onlya or not onlyb:
            fail.append((r['a'], r['b'], "平面が足りない"))
            continue
        P3, P4 = onlya[0], onlyb[0]
        # 符号を選ぶ: Q(v) と Q(w) が逆符号になるように
        made = None
        # **符号の選択を外した対照**（`CLAUDE.md`「機構が発火したことを数える」）。
        # 外すと辺の外に出る辺があるはず — 無ければ検査が空回りしています
        cand = (1,) if NOSIGN else (1, -1)
        for s4 in cand:
            Q = tuple(P3[k] + s4 * P4[k] for k in range(4))
            qa, qb = ev(Q, A), ev(Q, B)
            if not NOSIGN and (qa == 0 or qb == 0 or (qa > 0) == (qb > 0)):
                continue
            m = meet3(P1, P2, Q)
            if m is None:
                continue
            made = (Q, m)
            break
        if made is None:
            fail.append((r['a'], r['b'], "符号の選び方が無い"))
            continue
        Q, m = made
        # 厳密に「辺の内部」か: m = A + t (B - A) で 0 < t < 1
        ea, eb, em = eu(A), eu(B), eu(m)
        t = None
        for k in range(3):
            if eb[k] != ea[k]:
                t = (em[k] - ea[k]) / (eb[k] - ea[k])
                break
        on_line = all(em[k] == ea[k] + t * (eb[k] - ea[k]) for k in range(3))
        inside = on_line and 0 < t < 1
        if inside:
            ok += 1
        else:
            fail.append((r['a'], r['b'], f"線上 {on_line} / t = {t}"))
        maxbits['w'] = max(maxbits['w'], bits(abs(m[3])))
        maxbits['xyz'] = max(maxbits['xyz'], max(bits(abs(m[k])) for k in range(3)))
        maxbits['Qn'] = max(maxbits['Qn'], max(bits(abs(Q[k])) for k in range(3)))
        maxbits['Qd'] = max(maxbits['Qd'], bits(abs(Q[3])))
    print("| 判定 | 辺の本数 |")
    print("|---|---:|")
    print(f"| **辺の内部に点ができた** | **{ok}** |")
    print(f"| できなかった | **{len(fail)}** |")
    for f in fail[:5]:
        print(f"    {f}")
    print()
    b = 21
    print("| 量 | 実測の最大 | 導出した上界 | 上界の式 |")
    print("|---|---:|---:|---|")
    print(f"| $Q$ の法線成分 | {maxbits['Qn']} | {2*b+4} | $2b+4$ |")
    print(f"| $Q$ のオフセット | {maxbits['Qd']} | {3*b+6} | $3b+6$ |")
    print(f"| 交点の $w$ | {maxbits['w']} | {6*b+13} | $6b+13$ |")
    print(f"| 交点の $x,y,z$ | {maxbits['xyz']} | {7*b+15} | $7b+15$ |")
    over = (maxbits['Qn'] > 2*b+4 or maxbits['Qd'] > 3*b+6
            or maxbits['w'] > 6*b+13 or maxbits['xyz'] > 7*b+15)
    print()
    print("**上界を超えた量があります**" if over else "**すべて上界の中に収まりました**")


NOSIGN = len(sys.argv) > 2 and sys.argv[2] == 'nosign'
main(sys.argv[1])
