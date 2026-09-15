#!/usr/bin/env python3
r"""保存配列 2 本に限定した A5 確認と、包含・内側成分の判定。

**被検体のブール演算・分類器・述語は 1 つも呼びません。**
面対の判定は既存の `inspect_geom.pair_a5`、巻き数は `inspect_geom.winding` に集約します。
**既存だから正しいとはせず、段 `controls` で必要な正負の経路を先に確認します。**

段（1 回の起動につき 1 つ）::

    controls   合成対照。**この検証で実際に使う経路**の正負を確認する
    a5         A 内部・B 内部それぞれの全面対を整数 AABB で絞り、`pair_a5` へ渡す
               （**B の成分間も含みます。A 対 B の交差を拒む検査ではありません**）
    nesting    A の巻き数の前提、B 両成分の向き・相互包含、内側成分と保存残差の照合

**対象外**: ビルド、再量子化、ブール演算、**格子による高速化**、他の対への拡大。
**整数 AABB による絞り込みは入れます**（**厳密な絞り込み**です。下記）。

**$R = -\sum_j S(\text{内側成分}_j)$ の導出**（**数の一致だけを根拠にしないため、先に書きます**）::

  $B$ の巻き数が 外 0 / 殻 1 / 内 2 の形なら、
  $S_B = \int w_B = \mathrm{vol}(c_0) + \mathrm{vol}(c_1) = S_0 + S_1$、
  集合としての $|B| = \mathrm{vol}(c_0) = S_B - S_1$。
  式 2 の残差は $R = (|A\cup B| + |A\cap B|) - (S_A + S_B)$ で、
  **出力が集合として正しく**、かつ $w_A \in \{0,1\}$ なら
  $|A\cup B| + |A\cap B| = |A| + |B| = S_A + S_B - S_1$。よって $R = -S_1$。

  **したがって一致は、次の 3 つが【同時に】成り立つときに予測される値です。**
  **(i) B の入れ子が測ったとおり、(ii) $w_A \in \{0,1\}$、
  (iii) 出力の $|A\cup B|$・$|A\cap B|$ が集合としての測度に等しい。**
  **一致は、この 3 つを分離しません。** (iii) はこの道具では検証しません。

**$w_A \in \{0,1\}$ の根拠**（**2 点の測定だけでは空間全体について何も言えません**）::

  A0-1 / A0-2（参照の妥当性）・A0-3（退化面 0）・A1（各無向辺に面 2 枚）・
  A2（各有向辺 1 回 ＝ 向きが整合）・A3（頂点リンクが単一閉路 ＝ 頂点多様体）・
  A4（重複面なし）・**A5（自己交差・自己接触なし）**
  ⟹ A は閉じた・向き付き・**埋め込まれた** 2 次元多様体で、連結成分 1 個
  ⟹ **Jordan–Brouwer** により $\mathbb{R}^3 \setminus S_A$ の連結成分はちょうど 2 つで、
     巻き数は各領域で定数
  ⟹ 非有界側 0（実測）・有界側 +1（実測）⟹ **$w_A \in \{0,1\}$ が全空間で成立**

  **この論証は A0〜A5 に依存します。** だから段 `a5` で A0〜A4 も**この実行の中で**測り、
  別実行のログを引かずに済むようにしています。

**整数 AABB の絞り込みが厳密であることの論証**（**controls は補助検査であって根拠ではありません**）::

  三角形は自分の AABB に含まれます。`bb_hit` は厳密不等号で比べるので、
  接する箱は「当たり」に倒します。箱が交わらなければ交わりは空で、
  共有頂点も存在し得ません（共有頂点は両方の箱に入るため）。よって $s=0$ かつ交わり空。
  §43.4 の規則により、この対は**必ず正規**です。

**結論の限定**::

  * **体積の一致から内側を選びません。** 内側かどうかは巻き数で決め、
    そのあとで体積を照合します。
  * 保存残差 R は**被検体の出力体積から作られた量**です。一致しても
    「ブール出力が正しい」ことにはなりません。
  * **この 2 模型についてだけ**の結果です。
"""
import argparse, hashlib, os, struct, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig
from fractions import Fraction as Fr

# 受理済みの保存配列（`data/logs/inspect/20260915-205108/`）。**入力はハッシュで固定します。**
EXPECT = {
    "A": "d99d6d3b45132ee22c6ae83e391c47c8e8b2f0215b2c851edbab044921f5cf6b",
    "B": "fd9a03316f3c49e9d5b13a77fcc41ad6879dbdac0c6f0c8cbf69b50a97610596",
}


def read_bin(path):
    with open(path, "rb") as f:
        b = f.read()
    nv, nf = struct.unpack_from("<II", b, 0)
    need = 8 + nv * 12 + nf * 12
    if len(b) != need:
        raise ValueError("長さが合いません: %d != %d（%s）" % (len(b), need, path))
    V = struct.unpack_from("<%di" % (nv * 3), b, 8)
    F = struct.unpack_from("<%dI" % (nf * 3), b, 8 + nv * 12)
    return ([tuple(V[3 * i:3 * i + 3]) for i in range(nv)],
            [tuple(F[3 * i:3 * i + 3]) for i in range(nf)],
            hashlib.sha256(b).hexdigest())


def load(a, name, log):
    path = os.path.join(a.arrays, "C_%s_quantized.bin" % name)
    P, F, h = read_bin(path)
    want = EXPECT[name]
    log("[%s] 入力 %s" % (name, path))
    log("[%s]   sha256 %s / 期待 %s / %s" % (name, h, want, "一致" if h == want else "★ 不一致"))
    if h != want:
        raise ValueError("入力のハッシュが違います（%s）" % name)
    log("[%s]   頂点 %d / 面 %d" % (name, len(P), len(F)))
    return P, F


def aabb(P, f):
    vs = [P[i] for i in f]
    return (min(v[0] for v in vs), min(v[1] for v in vs), min(v[2] for v in vs),
            max(v[0] for v in vs), max(v[1] for v in vs), max(v[2] for v in vs))


def bb_hit(x, y):
    return not (x[3] < y[0] or y[3] < x[0] or x[4] < y[1] or y[4] < x[1]
                or x[5] < y[2] or y[5] < x[2])


# ===== 段 controls =========================================================
def stage_controls(a, log):
    """**この検証で使う経路**の正負を、合成入力で確認します。"""
    res = []

    def chk(name, got, want):
        ok = (got == want)
        res.append(ok)
        log("  %s %-52s 得た値 %-22s 期待 %s" % ("ok  " if ok else "★NG", name, got, want))

    # --- pair_a5 の正負 ---
    P = [(0, 0, 0), (10, 0, 0), (0, 10, 0), (10, 10, 0),      # 0..3  z=0 平面
         (0, 0, 10), (100, 0, 0), (100, 10, 0), (100, 0, 10),  # 4..7
         (2, 2, -5), (2, 2, 5), (8, 3, 0)]                     # 8,9,10
    f_a = (0, 1, 2)
    chk("s=0 離れた 2 枚は正規", ig.pair_a5(P, f_a, (5, 6, 7))[1:], (True, None))
    chk("s=2 辺を共有し交わりがその辺だけ", ig.pair_a5(P, f_a, (1, 2, 3))[1:], (True, None))
    chk("s=1 頂点だけを共有", ig.pair_a5(P, f_a, (1, 5, 7))[1:], (True, None))
    # 共面でなく、相対内部を横切る → 交差
    chk("非共面・相対内部を横切る → 交差", ig.pair_a5(P, f_a, (8, 9, 10))[1:], (False, "交差"))
    # 共面で面積を共有 → 接触
    P2 = P + [(1, 1, 0), (9, 1, 0), (1, 9, 0)]                 # 11..13
    chk("共面で面積を共有 → 接触", ig.pair_a5(P2, f_a, (11, 12, 13))[1:], (False, "接触"))
    # ★ 共面・面積 0・交わりが空でない（保全コードが読み飛ばしていた経路）
    P3 = [(0, 0, 0), (10, 0, 0), (0, 10, 0),                   # 0..2
          (5, 0, 0), (15, 0, 0), (5, -10, 0)]                  # 3..5  辺 z=0,y=0 上で線分を共有
    chk("共面・面積 0・線分で接する → 接触", ig.pair_a5(P3, (0, 1, 2), (3, 4, 5))[1:],
        (False, "接触"))
    P4 = [(0, 0, 0), (10, 0, 0), (0, 10, 0),
          (10, 0, 0), (20, 0, 0), (10, -10, 0)]                # 頂点 1 点だけで接する（共有添字なし）
    chk("共面・面積 0・1 点で接する → 接触", ig.pair_a5(P4, (0, 1, 2), (3, 4, 5))[1:],
        (False, "接触"))
    # 非共面で、辺上の 1 点だけに触れる（共有添字なし）→ 接触
    P5 = [(0, 0, 0), (10, 0, 0), (0, 10, 0),
          (5, 0, 0), (5, 0, 10), (9, 4, 5)]
    chk("非共面・辺上の 1 点だけ → 接触", ig.pair_a5(P5, (0, 1, 2), (3, 4, 5))[1:],
        (False, "接触"))

    # --- 整数 AABB による絞り込みが【厳密】であること ---
    # AABB が交わらなければ三角形も交わらないので、必ず s=0 で正規。
    # **決定的な擬似乱数**（`loader.hpp:173-178` と同じ xorshift64*、種 1）で作った
    # 三角形で、AABB が交わらない対を全部確かめます。
    st = [1]

    def nxt():
        s = st[0]
        s ^= (s >> 12) & 0xFFFFFFFFFFFFFFFF
        s ^= (s << 25) & 0xFFFFFFFFFFFFFFFF
        s ^= (s >> 27) & 0xFFFFFFFFFFFFFFFF
        st[0] = s & 0xFFFFFFFFFFFFFFFF
        return ((st[0] * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF) >> 11

    QP, QF = [], []
    for _ in range(120):
        base = tuple(nxt() % 60 for _ in range(3))
        tri = []
        for _ in range(3):
            QP.append(tuple(base[k] + nxt() % 9 for k in range(3)))
            tri.append(len(QP) - 1)
        QF.append(tuple(tri))
    BB = [aabb(QP, f) for f in QF]
    miss = skipped = 0
    for i in range(len(QF)):
        for j in range(i + 1, len(QF)):
            if bb_hit(BB[i], BB[j]):
                continue
            skipped += 1
            s, ok, kind = ig.pair_a5(QP, QF[i], QF[j])
            if not (s == 0 and ok):
                miss += 1
    chk("AABB が交わらない対は必ず s=0 かつ正規（%d 対）" % skipped, miss, 0)
    res.append(skipped > 0)
    log("  %s AABB で外した対が 0 件でないこと（対照が空回りしていないか）  得た値 %d"
        % ("ok  " if skipped > 0 else "★NG", skipped))

    # --- winding の正負 ---
    C = ig.box((0, 0, 0), (10, 10, 10), +1)
    chk("立方体の内側 → 1", ig.winding(C[0], C[1], (5, 5, 5))[0], 1)
    chk("立方体の外側 → 0", ig.winding(C[0], C[1], (50, 5, 5))[0], 0)
    chk("面の上の点 → 測れない", ig.winding(C[0], C[1], (0, 5, 5))[0], None)
    # 入れ子（外 +1 / 内 +1）: 内側の内側で総巻き数 2 —— 内側成分の判定に直結
    N = ig.join(ig.box((0, 0, 0), (12, 12, 12), +1), ig.box((3, 3, 3), (9, 9, 9), +1))
    chk("入れ子の内側の内側 → 2", ig.winding(N[0], N[1], (6, 6, 6))[0], 2)
    chk("入れ子の外側の殻の中 → 1", ig.winding(N[0], N[1], (Fr(3, 2), 6, 6))[0], 1)
    chk("入れ子の外 → 0", ig.winding(N[0], N[1], (50, 6, 6))[0], 0)
    # 開いた曲面（面を 1 枚落とす）→ 測れない
    O = (C[0], C[1][:-1])
    chk("開いた曲面の内側 → 測れない", ig.winding(O[0], O[1], (5, 5, 5))[0], None)
    # 成分分解
    D = ig.join(ig.box((0, 0, 0), (5, 5, 5), +1), ig.box((20, 0, 0), (25, 5, 5), +1))
    chk("離れた 2 個 → 成分 2 個", len(ig.components(D[1])), 2)
    chk("入れ子 → 成分 2 個", len(ig.components(N[1])), 2)

    # --- すぐ内側の点の取り方（向きを仮定しないこと）---
    #   **巻き数 0 の側を「内側」として受け取らないこと**も、ここで確かめます。
    def quiet(_):
        pass
    for sign, nm, want in ((+1, "外向き", 1), (-1, "反転（内向き）", -1)):
        BX = ig.box((0, 0, 0), (10, 10, 10), sign)
        q, w, wo, wt = probe_inside(BX[0], BX[1], BX[1], [], quiet, "")
        chk("すぐ内側（%s の箱）の巻き数" % nm, w, want)
        chk("すぐ内側（%s の箱）の他成分は空" % nm, wo, [])
        chk("すぐ内側（%s の箱）の全体の巻き数" % nm, wt, want)
    # 入れ子（外 +1 / 内 +1）で、内側の成分のすぐ内側の総巻き数が 2 になること
    NC = ig.components(N[1])
    NF = [[N[1][i] for i in c] for c in NC]
    got = []
    for j in range(len(NF)):
        _, wj, woj, wtj = probe_inside(N[0], NF[j], N[1],
                                       [NF[k] for k in range(len(NF)) if k != j],
                                       quiet, "")
        got.append((wj, tuple(woj), wtj))
    chk("入れ子: すぐ内側の（自分, 他, 全体）が 1 個は (1,(1,),2)",
        sum(1 for g in got if g[2] == 2), 1)
    chk("入れ子: 残りは総巻き数 1", sorted(g[2] for g in got), [1, 2])

    # --- 符号つき 6 倍体積（結論に直結するのに、対照が無かった）---
    U = ig.box((0, 0, 0), (1, 1, 1), +1)
    chk("単位立方体の符号つき 6 倍体積", ig.signed_volume6(U[0], U[1]), 6)
    UI = ig.box((0, 0, 0), (1, 1, 1), -1)
    chk("向きを反転すると符号が反転", ig.signed_volume6(UI[0], UI[1]), -6)
    chk("辺 10 の立方体は 6*1000", ig.signed_volume6(C[0], C[1]), 6000)
    chk("入れ子は成分の和", ig.signed_volume6(N[0], N[1]),
        ig.signed_volume6(*ig.box((0, 0, 0), (12, 12, 12), +1))
        + ig.signed_volume6(*ig.box((3, 3, 3), (9, 9, 9), +1)))

    ng = res.count(False)
    log("通過 %d / 失敗 %d" % (res.count(True), ng))
    log("判定: %s" % ("合成対照は全件通過" if ng == 0 else "★ 失敗あり。後続を起動しません"))
    return 0 if ng == 0 else 1


# ===== 段 a5 ===============================================================
def stage_a5(a, log):
    """A 内部・B 内部の**全面対**を整数 AABB で絞り、`pair_a5` へ渡します。

    **絞り込みは厳密です**: AABB が交わらなければ三角形は交わらないので、
    §43.4 の s=0 の規則により必ず正規。**段 `controls` で確認済み。**
    **格子は使いません**（候補の漏れを検証する仕事を増やさないため）。
    """
    ok = True
    for name in ("A", "B"):
        P, F = load(a, name, log)
        # ★ A0〜A4 も【この実行の中で】測ります。
        #   $w_A \in \{0,1\}$ の論証がこれらに依存するので、
        #   別実行のログを引くと依存が追跡できなくなります。
        oor, dup, col = ig.check_a0(P, F)
        log("[%s] A0-1 参照が範囲内: %s（違反 %d 枚）"
            % (name, "通過" if not oor else "★ 破れ", len(oor)))
        log("[%s] A0-2 3 頂点が相異: %s（違反 %d 枚）"
            % (name, "通過" if not dup else "★ 破れ", len(dup)))
        if oor or dup:
            log("[%s] ★ 参照の妥当性が破れました。後続を起動しません。" % name)
            return 1
        log("[%s] A0-3 退化面（外積が零）: %d 枚" % (name, len(col)))
        log("[%s]   ★ pair_a5 は退化面を含む対を無条件に正規と返します。"
            "**%d 枚なので、検査していない対は %s。**"
            % (name, len(col), "ありません" if not col else "あります"))
        a1, a2 = ig.check_a1(F), ig.check_a2(F)
        a3, a4 = ig.check_a3(F), ig.check_a4(F)
        for tag, bad in (("A1 各無向辺に面 2 枚", a1), ("A2 各有向辺が 1 回", a2),
                         ("A3 頂点リンクが単一閉路", a3), ("A4 重複面なし", a4)):
            log("[%s] %s: %s（違反 %d）"
                % (name, tag, "通過" if not bad else "★ 破れ", len(bad)))
        if a1 or a2 or a3 or a4:
            log("[%s] ★ A1〜A4 が破れました。埋め込まれた多様体の論証が使えません。" % name)
            return 1
        comps = ig.components(F)
        log("[%s] 成分 %d 個（面数 %s）" % (name, len(comps), [len(c) for c in comps]))
        of = {}
        for ci, c in enumerate(comps):
            for i in c:
                of[i] = ci
        BB = [aabb(P, f) for f in F]
        n = len(F)
        total = n * (n - 1) // 2
        judged = skipped = 0
        cross = {}          # 成分をまたぐ対の件数
        bad = []
        for i in range(n):
            bi = BB[i]
            fi = F[i]
            for j in range(i + 1, n):
                if not bb_hit(bi, BB[j]):
                    skipped += 1
                    continue
                judged += 1
                if of[i] != of[j]:
                    k = tuple(sorted((of[i], of[j])))
                    cross[k] = cross.get(k, 0) + 1
                s, good, kind = ig.pair_a5(P, fi, F[j])
                if not good:
                    bad.append((i, j, s, kind))
        nx = sum(1 for b in bad if b[3] == "交差")
        nt = len(bad) - nx
        log("[%s] 全対 %d / AABB で外した %d / 判定した %d" % (name, total, skipped, judged))
        log("[%s] 成分をまたぐ判定対: %s" % (name, cross if cross else "無し"))
        log("[%s] 交差 %d 件 / 接触 %d 件" % (name, nx, nt))
        if bad:
            for b in bad[:10]:
                log("[%s]   不正 面 %d-%d（共有頂点 %d、%s）" % (name, b[0], b[1], b[2], b[3]))
            if len(bad) > 10:
                log("[%s]   …ほか %d 件" % (name, len(bad) - 10))
        log("[%s] A5: %s" % (name, "通過" if not bad else "★ 破れ"))
        ok = ok and not bad
    log("判定: %s" % ("A5 は A・B とも通過" if ok else "★ A5 が破れました。後続を起動しません"))
    return 0 if ok else 1


# ===== 段 nesting ==========================================================
def probe_inside(P, Fc, F_all, others, log, tag):
    r"""成分 `Fc` の**すぐ内側**の点を 1 つ返す。

    面の重心から法線方向へ eps だけ寄せ、**この点で使うすべての巻き数**
    （自分の成分・他の各成分・全体）が**2 段続けて同じ値**になるまで eps を半分にします。

    **自分の成分の巻き数だけで止めてはいけません。** 他成分の面が重心と探針点の
    あいだを通る配置では、包含の判定が eps に依存して変わり得ます。
    A5 が保証するのは「交わらない」ことだけで、**距離の下界ではありません**。

    **巻き数が 0 の側を「内側」として受け取りません。** 0 を受け取ると、
    外側の点を「すぐ内側」と呼ぶことになります。
    **法線が外向きとは限らないので、向きを仮定せず両側を試します。**

    戻り値: `(点, 自分の巻き数, 他成分の巻き数の一覧, 全体の巻き数)`。
    確定できなければ `(None, None, None, None)`。**推測で埋めません。**
    """
    # 退化面（外積が零）は法線が取れないので避けます。
    f = None
    for g in Fc:
        if ig.normal(P, g) != (0, 0, 0):
            f = g
            break
    if f is None:
        log("  %s ★ 法線の取れる面がありません（全部が退化面）" % tag)
        return None, None, None, None
    c = ig.face_centroid(P, f)
    nv = ig.normal(P, f)
    scale = max(abs(x) for x in nv)

    def measure(q):
        ws, _, _ = ig.winding(P, Fc, q)
        if ws is None or ws == 0:
            return None
        wo = []
        for g in others:
            w, _, _ = ig.winding(P, g, q)
            if w is None:
                return None
            wo.append(w)
        wt, _, _ = ig.winding(P, F_all, q)
        if wt is None:
            return None
        return (ws, tuple(wo), wt)

    for sign, side in ((-1, "法線の逆側"), (+1, "法線の側")):
        prev = None
        for k in range(1, 41):
            eps = Fr(1, (1 << k) * scale)
            q = tuple(c[t] + sign * eps * nv[t] for t in range(3))
            cur = measure(q)
            if cur is None:
                prev = None
                continue
            if prev is not None and prev == cur:
                log("  %s すぐ内側の点を %s・k=%d で確定"
                    "（自分 %+d / 他成分 %s / 全体 %+d。**k-1 と k で全部同値**）"
                    % (tag, side, k, cur[0], list(cur[1]), cur[2]))
                return q, cur[0], list(cur[1]), cur[2]
            prev = cur
    log("  %s ★ すぐ内側の点を確定できませんでした（両側とも、k を 40 まで試しました）" % tag)
    return None, None, None, None


def stage_nesting(a, log):
    saved = saved_row(a.saved, a.key)
    log("保存値: S(A)=%d / S(B)=%d / R=%d" % (saved["S_A"], saved["S_B"], saved["R"]))
    out = {}
    for name in ("A", "B"):
        P, F = load(a, name, log)
        comps = ig.components(F)
        Fs = [[F[i] for i in c] for c in comps]
        S = [ig.signed_volume6(P, g) for g in Fs]
        log("[%s] 成分 %d 個 / 符号つき 6 倍体積 %s" % (name, len(comps), S))

        # 外側の点（AABB の外）で総巻き数が 0 であること
        hi = max(max(v) for v in P) + 1
        w_out, _, nv_out = ig.winding(P, F, (hi, hi, hi))
        log("[%s] 外側の点 (%d,%d,%d) の総巻き数 %s（有効方向 %d 本）"
            % (name, hi, hi, hi, w_out, nv_out))
        if w_out != 0:
            log("[%s] ★ 外側で 0 になりません。後続を起動しません。" % name)
            return 1

        inner = []
        for j, g in enumerate(Fs):
            others = [Fs[k] for k in range(len(Fs)) if k != j]
            oidx = [k for k in range(len(Fs)) if k != j]
            q, ws, wo, wt = probe_inside(P, g, F, others, log,
                                         "[%s] 成分 %d:" % (name, j))
            if q is None:
                log("[%s] ★ 成分 %d のすぐ内側を取れませんでした。" % (name, j))
                return 1
            log("[%s] 成分 %d: 体積の符号 %s / すぐ内側の自分の巻き数 %+d / 一致 %s"
                % (name, j, "+" if S[j] > 0 else "-", ws,
                   "する" if (S[j] > 0) == (ws > 0) else "★ しない"))
            for t, k in enumerate(oidx):
                log("[%s] 成分 %d のすぐ内側で、成分 %d の巻き数 %+d → %s"
                    % (name, j, k, wo[t],
                       "成分 %d は成分 %d の内側" % (j, k) if wo[t] != 0 else "外側"))
            log("[%s] 成分 %d のすぐ内側の総巻き数 %+d（成分ごとの和 %+d）"
                % (name, j, wt, ws + sum(wo)))
            log("[%s]   ★ この 2 つが一致するのは**恒等式**です"
                "（方向の有効性も計数も面ごとの和なので）。**独立な確認ではありません。**"
                % name)
            if wt != ws + sum(wo):
                log("[%s] ★ 恒等式が破れました。実装の誤りです。" % name)
                return 1
            if wt >= 2:
                inner.append(j)
        log("[%s] 内側成分（すぐ内側の総巻き数 ≥ 2）: %s" % (name, inner if inner else "無し"))
        out[name] = {"S": S, "inner": inner,
                     "S_inner": sum(S[j] for j in inner)}

    # --- 保存残差との照合（**体積の一致から内側を選んでいません**）---
    si = out["A"]["S_inner"] + out["B"]["S_inner"]
    log("Σ S(内側成分) = %d（A %d + B %d）" % (si, out["A"]["S_inner"], out["B"]["S_inner"]))
    log("保存残差 R = %d / -Σ S(内側成分) = %d / %s"
        % (saved["R"], -si, "一致" if saved["R"] == -si else "★ 不一致"))
    log("**内側の判定は巻き数で行いました。体積の一致から内側を選んでいません。**")
    log("**R は被検体の出力体積から作られた量です。一致してもブール出力の正しさは示しません。**")
    log("**この 2 模型についてだけの結果です。**")
    log("判定: %s" % ("内側成分の体積と保存残差が一致" if saved["R"] == -si else "★ 不一致"))
    return 0 if saved["R"] == -si else 1


def saved_row(path, key):
    r"""保存物から 1 行。**添字は 0 起点**（`t[27]` は 1 起点の列 28）。

    **列の意味は、出力を書いた側の定義から取っています**（値の一致を根拠にしません）:
    `tests/thingi10k/thingi_cp1.cpp` の列コメントで、
    列 28 = **式 2 の残差**（定義は `tests/gmp_identity.hpp` の `identity2_residual`。
    **$w_A, w_B \in \{0,1\}$ を前提にした式です**）、列 33 / 34 = 入力 A / B の
    符号つき 6 倍体積。
    **列が挿入・並べ替えされると黙って別の列を読む**ので、列数を下から検査します。
    """
    need = 34
    for line in open(path):
        t = line.split()
        if t and t[0] == key:
            if len(t) < need:
                raise ValueError("列が %d 個しかありません（%d 以上のはず）" % (len(t), need))
            return {"S_A": int(t[32]), "S_B": int(t[33]), "R": int(t[27])}
    raise KeyError("保存物に %s がありません" % key)


STAGES = {"controls": stage_controls, "a5": stage_a5, "nesting": stage_nesting}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True, choices=sorted(STAGES))
    p.add_argument("--out", required=True)
    p.add_argument("--arrays", default="data/logs/inspect/20260915-205108",
                   help="受理済みの保存配列の置き場所（**再量子化しません**）")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    path = os.path.join(a.out, "stage_%s.txt" % a.stage)
    lines = []

    def log(s):
        print(s)
        sys.stdout.flush()
        lines.append(s)

    log("=== 段 %s ===" % a.stage)
    t0 = time.monotonic()
    try:
        rc = STAGES[a.stage](a, log)
    finally:
        lines.append("経過 %.2f 秒" % (time.monotonic() - t0))
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
