#!/usr/bin/env python3
"""合成の .kmesh を作ります（駆動そのものの回帰試験のため）。

**実データを「拒否されるはず」と言って試す形を避けるため**に要ります
（`DESIGN-phase5-vertex-level.md` §22.25）。

書式は `tests/thingi10k/loader.hpp` の `load_kmesh` と対:
    "KMSH" + uint32 版(1) + uint32 頂点数 + uint32 面数 + float64[3*nv] + uint32[3*nf]
"""
import struct
import sys


def cube(ox, oy, oz, s):
    v = [(ox, oy, oz), (ox + s, oy, oz), (ox + s, oy + s, oz), (ox, oy + s, oz),
         (ox, oy, oz + s), (ox + s, oy, oz + s), (ox + s, oy + s, oz + s), (ox, oy + s, oz + s)]
    f = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
         (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]
    return v, f


def write(path, v, f):
    with open(path, "wb") as fp:
        fp.write(b"KMSH")
        fp.write(struct.pack("<III", 1, len(v), len(f)))
        for p in v:
            fp.write(struct.pack("<ddd", *[float(x) for x in p]))
        for t in f:
            fp.write(struct.pack("<III", *t))


if __name__ == "__main__":
    # 使い方: mk_synth_kmesh.py <出力先ディレクトリ>
    out = sys.argv[1]
    write(f"{out}/1001.kmesh", *cube(0, 0, 0, 4))
    write(f"{out}/1002.kmesh", *cube(2, 0, 0, 4))
    print("1001 12\n1002 12")
