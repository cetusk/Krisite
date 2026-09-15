#!/usr/bin/env python3
"""GMP 整合性診断の投入（DESIGN-phase5-vertex-level.md §22.26）。

**予算と上限を守るのはこの層の仕事です。** 駆動は 1 対ずつ子プロセスで回します。

  * **全体の期限 D**（既定 1200 秒 = 20 分）。**準備・終了猶予をすべて内側に含みます。**
  * **各対の期限** = ``min(per_pair, 残り − 猶予)``。**0 以下なら起動しません。**
  * **D − 猶予** までに ``TERM``、**D** までに ``KILL``。
  * **1 ワーカープロセスの仮想アドレス空間**を ``RLIMIT_AS`` で制限
    （**RSS の制限ではありません**。RSS は ``wait4`` の ``ru_maxrss`` で**測るだけ**）。

**OS による終了・回収の遅延があれば、超過として報告します。**
**「必ず D 以内に回収できる」とは保証しません。**
"""
import argparse
import os
import resource
import signal
import subprocess
import sys
import time


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="駆動（thingi_cp1_gmp）")
    ap.add_argument("--list", required=True, help="模型の一覧（変換の添字の元）")
    ap.add_argument("--only", required=True, help="対象の対の一覧（1 行 1 対）")
    ap.add_argument("--args", required=True, help="駆動へ渡す引数列（第 2 引数以降）")
    ap.add_argument("--deadline", type=float, default=1200.0, help="全体の期限（秒）")
    ap.add_argument("--grace", type=float, default=10.0, help="終了猶予（秒。期限の内側）")
    ap.add_argument("--per-pair", type=float, default=300.0, help="1 対の上限（秒）")
    ap.add_argument("--as-gib", type=float, default=8.0, help="仮想アドレス空間の上限（GiB）")
    ap.add_argument("--workdir", default=None, help="1 対ぶんの一覧を置く場所")
    ap.add_argument("--logdir", default=None, help="対ごとの出力を残す場所")
    ap.add_argument("--dry-run", action="store_true", help="起動せず、計画だけ出す")
    return ap.parse_args(argv)


def preexec(as_bytes):
    def f():
        resource.setrlimit(resource.RLIMIT_AS, (as_bytes, as_bytes))
        os.setsid()          # **子を独立した群にして、群ごと止められるようにします**
    return f


def run_one(a, key, only_path, deadline_s, logdir):
    """1 対を子プロセスで回します。戻り値: (状態, 終了値, 秒, ピーク RSS の KiB, ログの道)。

    **★ 子の出力は必ずファイルに残します。** 捨てると、**診断の記録そのものが消えます**
    （最初の実装は捨てていました)。
    """
    env = dict(os.environ, KRI_GMP_ONLY=only_path)
    cmd = [a.bin, a.list] + a.args.split()
    log_path = os.path.join(logdir, key.replace("/", "_") + ".log")
    t0 = time.monotonic()
    with open(log_path, "wb") as lf:
        p = subprocess.Popen(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT,
                             preexec_fn=preexec(int(a.as_gib * (1 << 30))))
        status, rc, rss = "ok", None, 0
        try:
            p.wait(timeout=deadline_s)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            status = "打ち切り"
            os.killpg(p.pid, signal.SIGTERM)
            try:
                p.wait(timeout=a.grace)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            rc = p.returncode
    dt = time.monotonic() - t0
    # **ピーク RSS は測るだけ**（制限ではありません）
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if status == "ok" and rc != 0:
        status = "失敗"
    return status, rc, dt, rss, log_path


def main(argv=None):
    a = parse_args(argv)
    for f in (a.bin, a.list, a.only):
        if not os.path.exists(f):
            print(f"**ありません**: {f}")
            return 2
    keys = [ln.split()[0] for ln in open(a.only) if ln.strip()]
    if not keys:
        print(f"**対象が空です**: {a.only}")
        return 2
    wd = a.workdir or os.path.dirname(os.path.abspath(a.only))
    logdir = a.logdir or wd
    os.makedirs(logdir, exist_ok=True)
    t_start = time.monotonic()
    t_end = t_start + a.deadline
    print(f"**これから {len(keys)} 対を回します**（全体の期限 {a.deadline:.0f} 秒、"
          f"各対 {a.per_pair:.0f} 秒、猶予 {a.grace:.0f} 秒、仮想アドレス空間 {a.as_gib} GiB）")
    for k in keys:
        print(f"  対象: {k}")
    if a.dry_run:
        print("**計画だけ出しました。起動していません。**")
        return 0

    done, skipped, bad = 0, 0, 0
    for k in keys:
        remain = t_end - time.monotonic()
        budget = min(a.per_pair, remain - a.grace)
        if budget <= 0:
            print(f"  {k}: **予算切れのため起動しません**（残り {remain:.1f} 秒）")
            skipped += 1
            continue
        only_path = os.path.join(wd, "_gmp_one.txt")
        with open(only_path, "w") as f:
            f.write(k + "\n")
        status, rc, dt, rss, log = run_one(a, k, only_path, budget, logdir)
        print(f"  {k}: {status}（終了値 {rc}、{dt:.1f} 秒、ピーク RSS {rss/1024:.0f} MiB、"
              f"期限 {budget:.1f} 秒）  ログ: {log}")
        if status != "ok":
            # **失敗したときは、その場で末尾を出します**（後で探さずに済むように）
            try:
                tail = open(log, errors="replace").read().splitlines()[-5:]
                for ln in tail:
                    print(f"      | {ln}")
            except OSError:
                pass
        if status == "ok":
            done += 1
        else:
            bad += 1
    over = time.monotonic() - t_end
    print(f"\n**済み {done} / 失敗・打ち切り {bad} / 未起動 {skipped}**")
    if over > 0:
        print(f"**★ 期限を {over:.1f} 秒超過しました**（OS の終了・回収の遅延）")
    return 0 if (bad == 0 and skipped == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
