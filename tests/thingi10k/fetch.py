#!/usr/bin/env python3
"""Thingi10K を取得して、Krisite の検証ハーネスが読む形に変換する（`SPEC-phase5.md` §2）。

**標準ライブラリだけで書いています。** `thingi10k` パッケージは `datasets` /
`polars` / `lagrange-open` を引き込むので使いません。npz は無圧縮の zip に
`.npy` が 2 つ入っているだけで、ヘッダは 1 行の Python 辞書リテラルです。

**データセットは Git に含めません**（`SPEC-phase5.md` §11）。取得先は
`data/thingi10k/`（`.gitignore` 済み）です。

使い方:

    python3 tests/thingi10k/fetch.py meta          # メタデータだけ取る
    python3 tests/thingi10k/fetch.py select cp1    # CP1 の 1,000 件を選んで一覧を書く
    python3 tests/thingi10k/fetch.py get cp1       # 一覧のメッシュを 1 件ずつ取得して変換する
    python3 tests/thingi10k/fetch.py bundle cp3    # **まとめた書庫を 1 回で取り、要るものだけ変換する**

出力する `.kmesh` は次の形式です（**リトルエンディアン固定**）。

    magic  "KMSH"          4 byte
    ver    uint32 = 1
    nv     uint32
    nf     uint32
    verts  float64 * 3nv
    faces  uint32  * 3nf
"""

import argparse
import tarfile
import tempfile
import csv
import hashlib
import os
import random
import struct
import sys
import time
import urllib.request
import zipfile

HF = "https://huggingface.co/datasets/Thingi10K/Thingi10K/resolve/main"
ROOT = os.path.join("data", "thingi10k")
META = os.path.join(ROOT, "metadata")
NPZ = os.path.join(ROOT, "npz")
MESH = os.path.join(ROOT, "kmesh")


def fetch(url, path):
    """未取得なら取る。**キャッシュは消さない**（再実行を安くするため）。"""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    os.replace(tmp, path)
    return path


def read_npy(buf):
    """`.npy` を読む（**numpy を使いません**）。返すのは (dtype, shape, bytes)。"""
    assert buf[:6] == b"\x93NUMPY", "npy ではない"
    major = buf[6]
    if major == 1:
        hlen = struct.unpack("<H", buf[8:10])[0]
        off = 10
    else:
        hlen = struct.unpack("<I", buf[8:12])[0]
        off = 12
    head = eval(buf[off : off + hlen].decode("latin1"))  # noqa: S307 — 自前の固定入力
    assert not head["fortran_order"], "fortran_order は想定外"
    return head["descr"], head["shape"], buf[off + hlen :]


def npz_to_kmesh(npz_path, out_path):
    """npz → `.kmesh`。**dtype は入力に合わせて変換します**（f4/f8、i4/i8 が来ます）。"""
    with zipfile.ZipFile(npz_path) as z:
        names = {n.split(".")[0]: n for n in z.namelist()}
        vd, vs, vb = read_npy(z.read(names["vertices"]))
        fd, fs, fb = read_npy(z.read(names["facets"]))
    nv, nf = vs[0], fs[0]
    assert vs[1] == 3 and fs[1] == 3, f"3 列でない: {vs} {fs}"

    vfmt = {"<f8": "d", "<f4": "f"}[vd]
    ffmt = {"<i4": "i", "<i8": "q", "<u4": "I", "<u8": "Q"}[fd]
    verts = struct.unpack("<" + vfmt * (nv * 3), vb[: nv * 3 * struct.calcsize(vfmt)])
    faces = struct.unpack("<" + ffmt * (nf * 3), fb[: nf * 3 * struct.calcsize(ffmt)])

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".part"
    with open(tmp, "wb") as f:
        f.write(b"KMSH" + struct.pack("<III", 1, nv, nf))
        f.write(struct.pack("<" + "d" * (nv * 3), *verts))
        f.write(struct.pack("<" + "I" * (nf * 3), *faces))
    os.replace(tmp, out_path)


def bundle(rows, which):
    """**まとめた書庫（`Thingi10K_npz.tar.gz`、4.09 GB）を 1 回だけ取り、要る模型だけ `.kmesh` にする。**

    **`get` が遅いのは大きさではなく【往復の回数】です**（実測: 1 件 2.33 s のうち転送は 0.2 s、
    残りは待ち時間。帯域は 2.1 MB/s）。5,535 件なら 3.6 時間ですが、**書庫なら 1 回・約 30〜60 分**です。

    **そして大量の個別取得は、提供側にとっても望ましくありません。**
    Hugging Face はこの書庫を**まとめ取り用に公開している**ので、そちらを使います。

    書庫は `npz/<file_id>.npz` の並び。**gzip は途中から読めない**ので全体を流しますが、
    **書き出すのは要る模型だけ**です（`.npz` は残しません。展開すると 14 GB になるため）。
    中断したら、変換済みを飛ばして流し直します（帯域だけ払い直し）。
    """
    need = {r["file_id"] for r in select(rows, which, 0, 20260830)}
    done = {os.path.splitext(f)[0] for f in os.listdir(MESH)} if os.path.isdir(MESH) else set()
    todo = need - done
    print(f"{which}: 母集団 {len(need)} 模型、変換済み {len(need & done)}、これから {len(todo)} 件を書庫から取り出します")
    if not todo:
        return
    url = f"{HF}/Thingi10K_npz.tar.gz"
    t0 = time.time()
    ok = err = seen = 0
    with urllib.request.urlopen(url, timeout=300) as resp:
        with tarfile.open(fileobj=resp, mode="r|gz") as tar:
            for m in tar:
                seen += 1
                if seen % 1000 == 0:
                    print(f"  書庫の {seen} 件目、取り出し {ok} 件、{time.time()-t0:.0f} s", flush=True)
                if not m.isfile():
                    continue
                fid = os.path.splitext(os.path.basename(m.name))[0]
                if fid not in todo:
                    continue
                try:
                    data = tar.extractfile(m).read()
                    with tempfile.NamedTemporaryFile(suffix=".npz", delete=False) as tf:
                        tf.write(data)
                        tmp = tf.name
                    npz_to_kmesh(tmp, os.path.join(MESH, f"{fid}.kmesh"))
                    os.unlink(tmp)
                    ok += 1
                except Exception as e:  # noqa: BLE001 — 1 件の失敗で全体を止めない
                    err += 1
                    print(f"  変換に失敗 {fid}: {e}", file=sys.stderr)
    el = time.time() - t0
    print(f"{which}: 取り出し {ok} 件 / 失敗 {err} 件 / 書庫の {seen} 件を通過 / {el:.0f} s → {MESH}")


def load_meta():
    fetch(HF + "/metadata/geometry_data.csv", os.path.join(META, "geometry_data.csv"))
    fetch(HF + "/metadata/contextual_data.csv", os.path.join(META, "contextual_data.csv"))
    with open(os.path.join(META, "geometry_data.csv")) as f:
        return list(csv.DictReader(f))


def num(r, k):
    try:
        return float(r.get(k, ""))
    except (TypeError, ValueError):
        return None


def tf(r, k):
    return str(r.get(k, "")).strip().lower() in ("1", "true", "t", "yes")


def select(rows, which, count, seed):
    """**選択は決定的です。** seed と条件が同じなら同じ一覧になります。"""
    if which == "cp1":
        # EMBER §5.1: solid・多様体・自己交差なし、1,000〜100,000 面
        pool = [
            r
            for r in rows
            if tf(r, "solid")
            and tf(r, "vertex_manifold")
            and tf(r, "edge_manifold")
            and (num(r, "num_self_intersections") or 0) == 0
            and 1000 <= (num(r, "num_faces") or 0) <= 100000
        ]
    elif which == "cp2":
        # 交差のあるモデル（自己交差入力。self-union が効くか）
        #
        # **面数の帯は CP1 と同じにします**（`SPEC-phase5.md` §2）。
        # **CP2 で変えるのは「自己交差の有無」だけ**です。
        #
        # 帯を切らないと、自己交差する模型は複雑な形状に偏るぶん規模も大きくなり
        # （実測: 中央 2.0 倍、90% 分位 5.1 倍、最大 19.8 倍）、
        # **「自己交差の効果」と「規模の効果」が交絡します。**
        #
        # **規模を広げるのは別の目的**です。CP3（全 10K）かその先で扱ってください。
        pool = [
            r
            for r in rows
            if (num(r, "num_self_intersections") or 0) > 0
            and 1000 <= (num(r, "num_faces") or 0) <= 100000
        ]
    elif which == "cp3":
        # 全 10K。**solid・多様体・自己交差の条件は課しません**
        # （PWN でない模型や退化した模型を含めて、入口の検査が働くかを試すため）。
        #
        # **ただし面数の帯は CP1 / CP2 と同じにします**（`SPEC-phase5.md` §2）。
        # **CP3 で変えるのは「模型の性質」だけ**です。
        #
        # **規模の上限を試すのは別の工程です**（§7.1）。CP3 に混ぜると、
        # 失敗したときに「入口の検査が捉えられない性質のせい」なのか
        # 「規模が大きすぎたせい」なのかが分かりません。
        pool = [r for r in rows if 1000 <= (num(r, "num_faces") or 0) <= 100000]
    else:
        raise SystemExit(f"未知の集合: {which}")
    pool.sort(key=lambda r: int(r["file_id"]))
    if count and count < len(pool):
        rnd = random.Random(seed)
        pool = sorted(rnd.sample(pool, count), key=lambda r: int(r["file_id"]))
    return pool


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["meta", "select", "get", "bundle"])
    ap.add_argument("which", nargs="?", default="cp1")
    ap.add_argument("--count", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260830)
    a = ap.parse_args()

    rows = load_meta()
    if a.cmd == "meta":
        print(f"メタデータ {len(rows)} 件を {META} に置きました")
        return

    if a.cmd == "bundle":
        bundle(rows, a.which)
        return
    sel = select(rows, a.which, a.count if a.which != "cp3" else 0, a.seed)
    lst = os.path.join(ROOT, f"{a.which}.txt")
    if a.cmd == "select":
        os.makedirs(ROOT, exist_ok=True)
        with open(lst, "w") as f:
            for r in sel:
                f.write(f"{r['file_id']} {int(num(r,'num_faces'))}\n")
        h = hashlib.sha256(open(lst, "rb").read()).hexdigest()[:16]
        print(f"{a.which}: {len(sel)} 件 → {lst}  sha256[:16]={h}  seed={a.seed}")
        return

    ok = err = 0
    for i, r in enumerate(sel):
        fid = r["file_id"]
        try:
            p = fetch(f"{HF}/npz/{fid}.npz", os.path.join(NPZ, f"{fid}.npz"))
            out = os.path.join(MESH, f"{fid}.kmesh")
            if not os.path.exists(out):
                npz_to_kmesh(p, out)
            ok += 1
        except Exception as e:  # noqa: BLE001 — 1 件の失敗で全体を止めない
            err += 1
            print(f"  取得に失敗 {fid}: {e}", file=sys.stderr)
        if (i + 1) % 100 == 0:
            print(f"  {i+1} / {len(sel)}", flush=True)
    print(f"{a.which}: 取得 {ok} 件 / 失敗 {err} 件 → {MESH}")


if __name__ == "__main__":
    main()
