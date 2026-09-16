#!/usr/bin/env python3
"""2 模型限定の入力再構成・構造検査・成分体積を、**段ごとに監督して**回します。

**既存の監督（`run_gmp_diag.py` の `spawn` / `supervise`）を再利用します。**
監督は改修しません。

  * **全段の絶対期限** … `--deadline`（既定 360 秒）。**終了猶予を含みます**
  * 段ごとの枠     … 下の `STAGES`（合計 300 秒）。**終了猶予を含みます**
  * `TERM` → `KILL` … `--grace`（既定 10 秒）。**枠の内側**
  * 回収の枠       … `--reap`（既定 30 秒）。**枠の外**
  * 仮想アドレス空間 … `--as-gib`（既定 4 GiB）。**1 プロセスあたり**

**後続を起動しない条件**（**自動延長も自動再試行もしません**）::

    回収不能 / 打ち切り / 子の終了値が 0 でない / 成果物が無いか空
    絶対期限までに終了猶予を確保できない

**種は一覧から導きます。** `data/thingi10k/cp3.txt` の 0 起点の添字に 1000 を足したもの
（`thingi_cp1.cpp:1690` の `prep[i] = prepare(raw, 1000 + i)`）。
**対の一覧の添字ではありません。**
"""
import argparse, os, subprocess, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_gmp_diag as sup

# (名前, 枠[秒], 命令を作る関数, 段の成果物)
STAGES = [
    ("C_A",         20),
    ("C_B",         20),
    ("reconstruct", 90),
    ("a0",          30),
    ("volume",      90),
    ("components",  50),
]
ARTIFACTS = {
    "C_A": ["C_A_transform.hex", "C_A_quantized.bin"],
    "C_B": ["C_B_transform.hex", "C_B_quantized.bin"],
    "reconstruct": ["stage_reconstruct.txt", "quantized_P_A.bin", "quantized_P_B.bin"],
    "a0": ["stage_a0.txt"],
    "volume": ["stage_volume.txt"],
    "components": ["stage_components.txt"],
}


def seed_of(list_path, model_id):
    """一覧の 0 起点の添字 + 1000。**見つからなければ例外**（推測で埋めません）。"""
    with open(list_path) as f:
        for i, line in enumerate(f):
            t = line.split()
            if t and t[0] == model_id:
                return 1000 + i, i
    raise KeyError("%s に %s がありません" % (list_path, model_id))


def cmd_of(name, a, out):
    if name == "C_A":
        s, _ = seed_of(a.list, a.a_id)
        return [a.dumper, a.a_kmesh, str(s), os.path.join(out, "C_A")]
    if name == "C_B":
        s, _ = seed_of(a.list, a.b_id)
        return [a.dumper, a.b_kmesh, str(s), os.path.join(out, "C_B")]
    return [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "inspect_components.py"),
            "--stage", name, "--out", out,
            "--a-kmesh", a.a_kmesh, "--b-kmesh", a.b_kmesh,
            "--a-id", a.a_id, "--b-id", a.b_id,
            "--saved", a.saved, "--key", a.key]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="")
    p.add_argument("--deadline", type=float, default=360.0)
    p.add_argument("--grace", type=float, default=10.0)
    p.add_argument("--reap", type=float, default=30.0)
    p.add_argument("--as-gib", type=float, default=4.0)
    p.add_argument("--dumper", default="build/tests/dump_quantized")
    p.add_argument("--list", default="data/thingi10k/cp3.txt")
    p.add_argument("--a-kmesh", default="data/thingi10k/kmesh/250394.kmesh")
    p.add_argument("--b-kmesh", default="data/thingi10k/kmesh/45413.kmesh")
    p.add_argument("--a-id", default="250394")
    p.add_argument("--b-id", default="45413")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    # ★ 引数を検査します。**`run_gmp_diag.py:218-224` が同じ検査を持っているのに、
    #   こちらだけ落ちていました。**
    for nm, v in (("--deadline", a.deadline), ("--grace", a.grace),
                  ("--reap", a.reap), ("--as-gib", a.as_gib)):
        if not (v == v and v not in (float("inf"), float("-inf"))) or v <= 0:
            p.error("%s は有限の正の数であること: %r" % (nm, v))
    if a.grace >= a.deadline:
        p.error("--grace (%g) は --deadline (%g) より小さいこと" % (a.grace, a.deadline))

    out = a.out or os.path.join("data", "logs", "inspect",
                                time.strftime("%Y%m%d-%H%M%S"))
    os.makedirs(out, exist_ok=True)
    # ★ 出力は書く前に消します。**前の実行の成果物が残っていると、
    #   「成果物がある」という検査が古いファイルで通ります。**
    stale = 0
    for names in ARTIFACTS.values():
        for f in names:
            q = os.path.join(out, f)
            if os.path.exists(q):
                os.unlink(q)
                stale += 1
    as_bytes = int(a.as_gib * 1024**3)

    t0 = time.monotonic()
    d_abs = t0 + a.deadline
    # ★ 設定と、実際に扱う対象を先に全部出します。**終わってから数えるのでは遅い。**
    sa, ia = seed_of(a.list, a.a_id)
    sb, ib = seed_of(a.list, a.b_id)
    print("これから %d 段を回します: %s" % (len(STAGES), " → ".join(n for n, _ in STAGES)))
    print("  対象 A = %s（一覧の添字 %d → 種 %d） / B = %s（添字 %d → 種 %d）"
          % (a.a_id, ia, sa, a.b_id, ib, sb))
    print("  一覧 %s / 保存物 %s / 鍵 %s" % (a.list, a.saved, a.key))
    print("  段の枠の合計 %g 秒 / 絶対期限 %g 秒 / 猶予 %g 秒（枠の内側）/ "
          "回収 %g 秒（枠の外）/ 1 プロセス %d バイト"
          % (sum(b for _, b in STAGES), a.deadline, a.grace, a.reap, as_bytes))
    print("  出力 %s（前の実行の成果物 %d 個を消しました）" % (out, stale))
    # ★ 試験専用の注入が効いたままだと、黙って挙動が変わります。**設定として出します。**
    inj = {k: os.environ[k] for k in
           ("KRI_DIAG_TEST_UNREAPED_AT", "KRI_DIAG_TEST_SPAWN_DELAY",
            "KRI_DIAG_TEST_CHILD_DELAY", "KRI_DIAG_TEST_SPAWN_LOG")
           if k in os.environ}
    print("  試験専用の注入: %s" % (inj if inj else "無し"))
    sys.stdout.flush()

    for name, budget in STAGES:
        now = time.monotonic()
        if d_abs - now <= a.grace:
            print("[%s] ★ 残り %.1f 秒では終了猶予 %g 秒を確保できません。起動しません。"
                  % (name, d_abs - now, a.grace))
            return 4
        t_kill = min(now + budget, d_abs)
        t_term = max(now, t_kill - a.grace)
        log = os.path.join(out, "child_%s.log" % name)
        cmd = cmd_of(name, a, out)
        print("[%s] 起動: %s" % (name, " ".join(cmd)))
        print("[%s]   TERM は %.1f 秒後 / KILL は %.1f 秒後 / 回収の枠 %g 秒"
              % (name, t_term - now, t_kill - now, a.reap))
        sys.stdout.flush()

        pid = sup.spawn(cmd, dict(os.environ), log, as_bytes)
        code, rss, state = sup.supervise(pid, t_term, t_kill, hard_cap=a.reap)
        dt = time.monotonic() - now
        try:
            sys.stdout.write(open(log, encoding="utf-8", errors="replace").read())
        except OSError:
            pass
        print("[%s] 状態 %s / 終了値 %s / ピーク RSS %s / 経過 %.2f 秒"
              % (name, state, "不明" if code is None else code,
                 "不明" if rss is None else "%d KiB" % rss, dt))
        sys.stdout.flush()

        if state == sup.UNREAPED:
            print("[%s] ★ 回収できていません（pid=%d、群 %d）。後続を起動しません。"
                  % (name, pid, pid))
            return 2
        if state == sup.CUT:
            print("[%s] ★ 打ち切り、または期限後の回収。後続を起動しません。" % name)
            return 3
        if code != 0:
            print("[%s] ★ 子が終了値 %s で終わりました。後続を起動しません。" % (name, code))
            return 1
        miss = [f for f in ARTIFACTS[name]
                if not os.path.exists(os.path.join(out, f))
                or os.path.getsize(os.path.join(out, f)) == 0]
        if miss:
            print("[%s] ★ 成果物が無いか空です: %s。後続を起動しません。" % (name, ", ".join(miss)))
            return 5

    print("全 %d 段を回しました。経過 %.2f 秒（絶対期限 %g 秒）"
          % (len(STAGES), time.monotonic() - t0, a.deadline))
    print("出力: %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
