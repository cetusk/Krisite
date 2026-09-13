# -*- coding: utf-8 -*-
"""`embed_check.py` の判定が【発火すること】を先に確かめる。

**0 件という結果が「埋め込まれている」なのか「検査が壊れている」なのかを
区別できないままでは使えません**（`CLAUDE.md`「機構が空回りしていないことを別に検査」）。
"""
import sys
sys.path.insert(0, sys.argv[1])
from fractions import Fraction
import importlib.util
spec = importlib.util.spec_from_file_location("ec", sys.argv[1] + "/embed_check.py")
ec = importlib.util.module_from_spec(spec)
# main() を走らせないため、ファイルを読み込んで main 呼び出し行だけ落とす
src = open(sys.argv[1] + "/embed_check.py", encoding='utf-8').read()
src = src.replace("main(sys.argv[1])", "")
exec(compile(src, "embed_check.py", "exec"), ec.__dict__)

F = Fraction
def V(*a): return tuple(F(x) for x in a)

cases = [
    # 名前, 扇 A, 扇 B, 期待
    ("直交する 2 枚が貫く", (V(1,0,0), V(0,1,0)), (V(1,1,-1), V(1,1,1)), True),
    ("離れた 2 枚",         (V(1,0,0), V(0,1,0)), (V(-1,0,0), V(0,-1,0)), False),
    ("同一平面で重なる",    (V(1,0,0), V(0,1,0)), (V(1,1,0), V(-1,1,0)), "coplanar"),
    ("同一平面で入れ子",    (V(1,0,0), V(0,1,0)), (V(2,1,0), V(1,2,0)), "coplanar"),
    ("同一平面だが離れている", (V(1,0,0), V(0,1,0)), (V(-1,0,0), V(0,-1,0)), False),
    ("同一平面で辺だけ接する", (V(1,0,0), V(0,1,0)), (V(0,1,0), V(-1,0,0)), False),
]
bad = 0
for name, (u1,u2), (w1,w2), want in cases:
    got = ec.arcs_cross(u1,u2,w1,w2)
    ok = (got == want)
    if not ok: bad += 1
    print(f"  {'OK ' if ok else '**NG**'} {name}: 期待 {want} / 実測 {got}")
print(f"\n**不一致 {bad} 件 / {len(cases)} 件**")
sys.exit(1 if bad else 0)
