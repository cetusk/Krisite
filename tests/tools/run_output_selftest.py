#!/usr/bin/env python3
r"""`run_output.py` の停止制御そのものを検定します（§44.8 の運用）。

**被検体は走らせません。** 偽の駆動（`/bin/false` など）と注入で、
**停止条件が実際に効くこと**と**終了値**を確かめます。

  * 子が非零で終わる → 1
  * 回収不能（注入）  → 2
  * 打ち切り          → 3
  * 猶予を確保できない → 4
  * 成果物が無いか空  → 5
  * 起動失敗（127）   → 6

**「テストが緑」と「テストが検証している」は別**なので、**期待する終了値を明示**し、
**外れたら落とします。**
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DRV = os.path.join(HERE, "run_output.py")


def run(args, env=None, cwd=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run([sys.executable, DRV] + args, capture_output=True, text=True,
                       env=e, cwd=cwd, timeout=300)
    return p.returncode, p.stdout + p.stderr


def main():
    root = os.path.abspath(os.path.join(HERE, "..", ".."))
    cases = []
    with tempfile.TemporaryDirectory() as td:
        def out(n):
            d = os.path.join(td, n)
            os.makedirs(d, exist_ok=True)
            return d

        cases.append(("子が非零 → 1", ["--out", out("a"), "--dumper", "/bin/false"], None, 1))
        cases.append(("成果物が無いか空 → 5",
                      ["--out", out("b"), "--dumper", "/bin/true"], None, 5))
        cases.append(("起動失敗（127）→ 6",
                      ["--out", out("c"), "--dumper", "/nonexistent_dumper_xyz"], None, 6))
        # ★★ 【取り下げ】以前ここに「終了値 4 は到達不能」と書きましたが、**誤りです**
        #   （`SPEC-phase5.md` §5.10.14.98）。子の完了確認の時刻と次段の起動判定の時刻の
        #   あいだには、ログの読込み・出力・flush・成果物の確認が入り、親が走らない時間も
        #   入り得ます。**ガードは必要です。**
        #   **そして終了値 4 は、この自己検定では観測できていません（未検定）。**
        #   下の構成で観測されるのは打ち切り（3）です。**「猶予不足の枝を検定済み」とは
        #   扱いません。**
        cases.append(("猶予が枠に近い構成では打ち切りで止まる（**終了値 4 は未検定**）",
                      ["--out", out("d"), "--deadline", "11", "--grace", "10"], None, 3))
        cases.append(("回収不能（注入）→ 2", ["--out", out("e")],
                      {"KRI_DIAG_TEST_UNREAPED_AT": "1"}, 2))
        cases.append(("打ち切り（子を遅らせる）→ 3",
                      ["--out", out("f"), "--deadline", "30"],
                      {"KRI_DIAG_TEST_CHILD_DELAY": "60"}, 3))

        ng = 0
        for name, args, env, want in cases:
            rc, log = run(args, env, cwd=root)
            ok = (rc == want)
            if not ok:
                ng += 1
            print("  %s %-36s 終了値 %s / 期待 %s" % ("ok  " if ok else "★NG", name, rc, want))
            if not ok:
                print("     --- 出力の末尾 ---")
                for line in log.strip().split("\n")[-6:]:
                    print("     " + line)
    print("通過 %d / 失敗 %d" % (len(cases) - ng, ng))
    return 0 if ng == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
