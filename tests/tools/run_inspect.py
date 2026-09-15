#!/usr/bin/env python3
"""合成の幾何検査を、**外部監督のもとで**起動します。

**既存の監督（`run_gmp_diag.py`）を再利用します。** 監督は作り直しません。

  * 計算＋終了猶予 … `--budget`（既定 90 秒）。**終了猶予を【含みます】。**
    この時刻に `KILL`。**`TERM` はその `--grace` 秒【手前】に送ります**
    （`t_kill = t0 + budget`、`t_term = t_kill - grace`）。
    **`budget` に `grace` を足して `KILL` の時刻にしてはいけません。**
    許可された枠は「計算＋終了猶予で 90 秒」なので、二重に足すと 100 秒になります。
  * `TERM` → `KILL` の猶予 … `--grace`（既定 10 秒）。**`budget` の内側**
  * 回収の枠       … `--reap`（既定 30 秒）。`KILL` の後これだけ待って回収。
    **これだけが `budget` の外側**（許可された「回収は追加 30 秒」）
  * 仮想アドレス空間 … `--as-gib`（既定 4 GiB）。子に `RLIMIT_AS` で掛ける
    （**常駐集合（RSS）ではありません。** RSS は測るだけで上限を掛けません）

  状態は `完了` / `打ち切り` / `回収不能` の 3 つに分かれます。
  **`回収不能` を `打ち切り` と同じに扱いません**（子が生きている可能性が残ります）。

親の終了値:

  ==  =================================================================
  0   `完了` かつ子の終了値 0
  1   `完了` だが子の終了値が 0 でない（**被検体の失敗**）
  3   `打ち切り`（**打ち切ったか、期限を過ぎてから回収した**。
      被検体の失敗とは別に立てます。**信号送出の有無は含みません**）
  2   `回収不能`（**子が残っている可能性がある**）
  ==  =================================================================

**この道具が示せないこと**（**主張をここまでに限定してください**）:

  * **正常に終わった実行は、停止・回収の制御を 1 度も発火させません。**
    出力の「停止制御」の行が、発火したかどうかを述べます。
    **枠の値が表示されたことは、枠が働いたことの根拠になりません。**
  * 停止は群（`killpg`）へ送るので、**子が自ら `setsid()` した子孫には届きません。**
    その場合でも直接の子は回収できるので状態は `打ち切り` になります。
    **「打ち切り」は「子孫まで含めて止まった」を意味しません。**
  * 子の終了値 **127** は、`run_gmp_diag.spawn` が子側の例外をすべて畳む値です。
    **起動できなかった（命令が無い／ログを開けない／`RLIMIT_AS` が小さすぎる）場合と、
    被検体が本当に 127 で終わった場合を区別できません。**

使い方::

    python3 tests/tools/run_inspect.py -- python3 tests/tools/inspect_geom_test.py
    python3 tests/tools/run_inspect.py --budget 2 --grace 1 --reap 5 -- <止まらない命令>
"""
import argparse, os, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_gmp_diag as sup          # 監督はここから再利用する（改修しません）


def _fmt(x):
    """設定の表示。**丸めて実値とずらさない**ようにします。"""
    return ("%d" % x) if float(x).is_integer() else ("%r" % x)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--budget", type=float, default=90.0, help="計算＋終了猶予（秒）")
    p.add_argument("--grace", type=float, default=10.0, help="TERM から KILL までの猶予（秒）")
    p.add_argument("--reap", type=float, default=30.0, help="KILL の後に回収を待つ上限（秒）")
    p.add_argument("--as-gib", type=float, default=4.0, help="子の仮想アドレス空間（GiB）")
    p.add_argument("--log", default="",
                   help="子の出力を書く先。**子が O_TRUNC で開くので、既存の内容は消えます。**"
                        "省くと一時ファイルに採り、終わってからこちらへ流します")
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        p.error("起動する命令を -- の後に書いてください")

    # ★ 既定を `/dev/stdout` にしてはいけません。
    #   子は `O_WRONLY|O_CREAT|O_TRUNC` で開くので、親の標準出力がファイルに
    #   向いていると【そのファイルを 0 に切り詰め】、親が終わるときに先頭から
    #   上書きします。**子の出力が黙って消え、消えたことが出力から分かりません。**
    #   さらに取り残した子が親の標準出力を握り続け、下流の読み手が固まります。
    #   → 一時ファイルに採り、監督が返ってからこちらへ流します。
    tmp = None
    if a.log:
        log = a.log
    else:
        fd, tmp = tempfile.mkstemp(prefix="run_inspect_child_", suffix=".log")
        os.close(fd)
        log = tmp

    as_bytes = int(a.as_gib * 1024**3)
    t0 = time.monotonic()
    # ★ `budget` は【終了猶予を含む】枠です。`KILL` の時刻がこれに一致します。
    #   `t_kill = t_term + grace` と書くと、許可された 90 秒が 100 秒になります。
    t_kill = t0 + a.budget
    t_term = max(t0, t_kill - a.grace)     # 猶予が枠より大きければ即座に TERM
    print("起動: %s" % " ".join(cmd))
    print("  計算＋終了猶予 %s 秒（猶予込み。TERM は開始から %.2f 秒、KILL は %s 秒）"
          " / 回収の枠 %s 秒（枠の外）/ 仮想アドレス空間 %s バイト（%s GiB を切り捨て）"
          % (_fmt(a.budget), t_term - t0, _fmt(a.budget),
             _fmt(a.reap), as_bytes, _fmt(a.as_gib)))
    print("  子の出力の行き先: %s%s" % (log, "（一時ファイル。下に流します）" if tmp else ""))
    sys.stdout.flush()

    pid = sup.spawn(cmd, dict(os.environ), log, as_bytes)
    code, rss, state = sup.supervise(pid, t_term, t_kill, hard_cap=a.reap)
    dt = time.monotonic() - t0

    # ★ 回収できていないときは流しません。**子がまだ書き足す可能性がある**ので、
    #   その時点の中身を「子の出力」として見せると、途中までを全部と読ませます。
    if tmp and state != sup.UNREAPED:
        print("  --- 子の出力 ここから ---")
        sys.stdout.flush()
        try:
            with open(tmp, "rb") as f:
                sys.stdout.buffer.write(f.read())
            sys.stdout.buffer.flush()
        except OSError as e:
            print("  （子の出力を読めませんでした: %s）" % e)
        print("  --- 子の出力 ここまで ---")

    print("  状態 %s / 終了値 %s / ピーク RSS %s / 経過 %.2f 秒"
          % (state, "不明" if code is None else code,
             "不明" if rss is None else ("%d KiB" % rss), dt))

    # ★ 「枠を表示した」と「枠が働いた」を分けます。
    #   `supervise` は信号を送ったかを返さないので、**状態から導きます**。
    if state == sup.DONE:
        print("  停止制御: **未発火**（経過 %.2f 秒 < 計算＋終了猶予 %s 秒。"
              "TERM も KILL も送っていません）" % (dt, _fmt(a.budget)))
    elif state == sup.CUT:
        # ★ `CUT` からは信号を送ったことを導けません。
        #   `run_gmp_diag.py:148-149` は、信号を送る【前】でも期限を過ぎてから
        #   回収できれば `CUT` を返します（`late` の枝）。
        #   **実際に信号が届いたことの根拠は、子の側の観測**（ハンドラの出力、
        #   終了値が負であること）**であって、この戻り値ではありません。**
        print("  停止制御: **打ち切り、または期限後の回収**"
              "（信号を送ったかどうかは、この戻り値だけでは分かりません。"
              "終了値 %s が負なら、子はシグナルで終わっています）" % code)
    else:
        # ★ 注入で `UNREAPED` を返させた場合、信号は 1 度も送られていません。
        #   `supervise` は信号を送ったかを返さないので、**発火したとは書けません**。
        print("  停止制御: **回収に至りませんでした**"
              "（信号を送ったかどうかは監督の戻り値から分かりません）")

    if code == 127:
        print("  ★ 終了値 127 は、起動できなかった場合と被検体が 127 で終わった場合を"
              "区別できません（`run_gmp_diag.spawn` が子側の例外を畳む値）。")

    if state == sup.UNREAPED:
        print("  ★ 回収できていません。後続を起動しません。")
        print("     残っている可能性のある子: pid=%d（群 %d）。"
              "後始末は `kill -9 -%d` か、`tests/tools/stopctl.py alive <argv の要素>` で確認してください。"
              % (pid, pid, pid))
        if tmp:
            print("     子の出力は %s に残しました（子がまだ書き足す可能性があるため消しません）。" % tmp)
        return 2

    if tmp:
        try:
            os.unlink(tmp)
        except OSError:
            pass
    if state == sup.CUT:
        return 3
    return 0 if code == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
