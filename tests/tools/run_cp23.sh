#!/bin/bash
# CP2 本体 + CP3 を続けて回す（`SPEC-phase5.md` §5.10.14.73 / §5.10.14.74）。
#
# **前提が揃わなければ計算を始めません。段が失敗したら次段を起動しません。**
#
#   bash tests/tools/run_cp23.sh            新規（結果ファイルが空であること）
#   bash tests/tools/run_cp23.sh --resume   再開（meta が一致し、FAIL が無いこと）
#
# **対応範囲は sbx の中（Linux / bash）だけ**です。Windows 側では動かしません。
#
# 差し替え用（自己検査が使います。既定は実物）:
#   KRI_ROOT / KRI_BIN / KRI_SHA / KRI_ARGS / KRI_BASES
set -u

fail() { printf '**中止**: %s\n' "$1" >&2; exit 2; }

ROOT="${KRI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT" || fail "作業ディレクトリへ移動できません: $ROOT"

BIN="${KRI_BIN:-build/thingi_cp1_o3}"
SHA="${KRI_SHA:-tests/tools/run_cp23.sha256}"
ARGS="${KRI_ARGS:-0 6 8 1 0 1 1 0 0 16 2 1 0 1}"
BASES="${KRI_BASES:-cp2b cp3}"
RESUME=0
[ "${1:-}" = "--resume" ] && RESUME=1

# ---- 起動前の検査（ここで止まれば計算は 1 秒も走りません）--------------------
[ -x "$BIN" ] || fail "実行ファイルがありません: $BIN"
# **指紋は機械ごとの値なので追跡しません。** 投入の直前に作ってください:
#   sha256sum build/thingi_cp1_o3 > tests/tools/run_cp23.sha256
[ -f "$SHA" ] || fail "指紋のファイルがありません: $SHA（sha256sum で作ってください）"
sha256sum -c "$SHA" > /dev/null 2>&1 || fail "指紋が一致しません（$SHA）"

meta_of() { printf '%s_results.meta' "data/thingi10k/$1"; }
res_of()  { printf '%s_results.txt'  "data/thingi10k/$1"; }
only_of() { printf '%s_only.txt'     "data/thingi10k/$1"; }

for base in $BASES; do
    only="$(only_of "$base")"
    [ -f "$only" ] || fail "標本の一覧がありません: $only"
    n=$(grep -c . "$only" || true)
    [ "$n" -gt 0 ] || fail "標本の一覧が空です: $only"
    # **版と予定キー集合**。開始時刻は履歴であって、一致判定には使いません
    want="$(sha256sum "$BIN" | cut -d' ' -f1) args=${ARGS} keys=$(sort "$only" | sha256sum | cut -d' ' -f1) n=${n}"
    res="$(res_of "$base")"
    meta="$(meta_of "$base")"
    if [ -s "$res" ]; then
        [ "$RESUME" = 1 ] || fail "結果が残っています（新規なら空にしてください）: $res"
        [ -f "$meta" ] || fail "再開なのに meta がありません: $meta"
        got="$(grep -v '^started=' "$meta" || true)"
        [ "$got" = "$want" ] || fail "再開の meta が一致しません: $meta"
        awk '{print $1}' "$res" | sort | uniq -d | grep -q . && fail "結果にキーの重複があります: $res"
        awk 'NF<3{exit 1}' "$res" || fail "結果に不正な行があります: $res"
        grep -q ' FAIL ' "$res" && fail "結果に既存の FAIL があります: $res"
    else
        printf '%s\nstarted=%s\n' "$want" "$(date -Is)" > "$meta" || fail "meta を書けません: $meta"
    fi
done

# ---- 実行（計算側と tee 側の状態を直後に両方保存する）-----------------------
for base in $BASES; do
    printf '\n=== %s を回します（%s）===\n' "$base" "$(date -Is)"
    log="data/thingi10k/${base}_run.log"
    set -o pipefail
    "$BIN" "data/thingi10k/${base}.txt" $ARGS 2>&1 | tee "$log"
    # **★ 直後に配列ごと保存します。** 先に `PIPESTATUS[0]` を代入すると、
    # その代入自体が `PIPESTATUS` を書き換えて `[1]` が消えます
    rc=("${PIPESTATUS[@]}")
    set +o pipefail
    rc_run=${rc[0]:-1}
    rc_tee=${rc[1]:-1}
    [ "$rc_run" = 0 ] || fail "${base}: 計算側が異常終了しました（exit ${rc_run}）"
    [ "$rc_tee" = 0 ] || fail "${base}: 記録（tee）が失敗しました（exit ${rc_tee}）"

    # ---- 完了の検査（終了コードだけでは判定しません）----
    res="$(res_of "$base")"; only="$(only_of "$base")"
    [ -s "$res" ] || fail "${base}: 結果が空です"
    awk '{print $1}' "$res" | sort | uniq -d | grep -q . && fail "${base}: 結果にキーの重複"
    if ! diff -q <(sort "$only") <(awk '{print $1}' "$res" | sort) > /dev/null; then
        fail "${base}: 予定キー集合と結果のキー集合が一致しません"
    fi
    nf=$(awk '$2=="FAIL"' "$res" | wc -l)
    [ "$nf" = 0 ] || fail "${base}: FAIL が ${nf} 件あります"
    printf '%s: 予定 %s 件すべてが成功しました\n' "$base" "$(grep -c . "$only")"
done
printf '\n完了（%s）\n' "$(date -Is)"
