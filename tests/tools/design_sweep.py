#!/usr/bin/env python3
"""docs/DESIGN-phase5-vertex-level.md §43 の機械的な突き合わせ。

9 巡のレビューで、確実な誤りの過半が「前の巡の修正が、同じ記号を使う
他の箇所に波及していない」形だった。記号ごとに全出現を並べ、
旧い書き方が残っていたら落ちる。

使い方: python3 tests/tools/design_sweep.py
"""
import io, re, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
p = os.path.join(HERE, "..", "..", "docs", "DESIGN-phase5-vertex-level.md")
L=io.open(p,encoding="utf-8").read().split("\n")
st=next(i for i,l in enumerate(L) if l.startswith("## 43."))
bad=0
# 1) 表の健全性
nt=0;i=st
while i<len(L):
    if L[i].startswith("|"):
        blk=[];j=i
        while j<len(L) and L[j].startswith("|"): blk.append(L[j]);j+=1
        nt+=1
        sep=len(blk)>1 and re.match(r'^\|[\s:\-|]+\|$',blk[1].strip())
        if not sep or len({r.count("|") for r in blk})!=1:
            print("★表が異常 行%d"%(i+1));bad+=1
        i=j
    else: i+=1
print("表 %d 個 / 異常 %d"%(nt,bad))
# 2) 記号ごとの全出現
def sweep(name,pat,forbid=None):
    global bad
    hits=[(k+1,L[k]) for k in range(st,len(L)) if re.search(pat,L[k])]
    bads=[h for h in hits if forbid and re.search(forbid,h[1])]
    print("  %-14s 出現 %2d 件%s"%(name,len(hits)," / ★違反 %d"%len(bads) if bads else ""))
    for b in bads: print("      行%d: %s"%(b[0],b[1][:110]))
    bad+=len(bads)
sweep("A5 の対象",   r"A5",        r"A1〜A5|A5[^|]*\*\*T のみ\*\*")
sweep("A1〜A4",      r"A1〜A4",    None)
sweep("|S|",         r"lvert S",   r"\\lvert S\\rvert")   # 添字なしの |S| は禁止
sweep("内側成分",     r"内側成分",   None)
sweep("H / H′",      r"\bH′|仮説 H|H 成立|H′ 成立", None)
sweep("P-静/P-動",   r"P-静|P-動|前提通過", None)
sweep("2,236,520",   r"2\{,\}236\{,\}520", None)
sweep("k_max",       r"k_\{\\max\}", r"k_\{\\max\} = 40")
sys.exit(1 if bad else 0)
