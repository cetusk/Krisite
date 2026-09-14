#!/bin/bash
# CP2 本体 + CP3 を続けて回す（`SPEC-phase5.md` §5.10.14.73 / §5.10.14.74）。
#
# **前提が揃わなければ計算を始めません。段が失敗したら次段を起動しません。**
#
#   bash tests/tools/run_cp23.sh            新規（結果ファイルが空であること）
#   bash tests/tools/run_cp23.sh --resume   再開（承認済みの基準と一致すること）
#
# **対応範囲は sbx の中（Linux / bash）だけ**です。Windows 側では動かしません。
#
# ## 承認済みの基準（`KRI_MANIFEST`。**追跡しません**。機械ごとの値）
#
#   bin=<sha256>                             投入を承認したバイナリ
#   <base> keys=<sha256> n=<件数> cols=<列数>  承認した標本と、結果行の【数値の列数】
#
# **`cols` は駆動の版で決まります**（`ps.print` が出す項目の数）。
# **明示を必須にします** — 既定値を置くと、版が変わったときに黙って通ります。
#
# **一覧から作った値どうしを比べても、誤った一覧を弾けません。**
# **独立した基準と照合します。**
#
# 差し替え用（自己検査が使います。既定は実物）:
#   KRI_ROOT / KRI_BIN / KRI_MANIFEST / KRI_ARGS / KRI_BASES
set -u

fail() { printf '**中止**: %s\n' "$1" >&2; exit 2; }

ROOT="${KRI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT" || fail "作業ディレクトリへ移動できません: $ROOT"

BIN="${KRI_BIN:-build/thingi_cp1_o3}"
MANIFEST="${KRI_MANIFEST:-tests/tools/run_cp23.manifest}"
ARGS="${KRI_ARGS:-0 6 8 1 0 1 1 0 0 16 2 1 0 1}"
BASES="${KRI_BASES:-cp2b cp3}"
RESUME=0
[ "${1:-}" = "--resume" ] && RESUME=1

# **前処理は一時ファイルに落とし、【各段の終了状態】を見ます。**
# **`sort "$f" | sha256sum` の形は、`sort` が失敗しても成功を返します**
# （パイプの終了状態は最後のコマンドのもの）。**空入力のハッシュが返ります。**
# **「重複なし」「差分なし」と取り違えないために、段ごとに確かめます。**
TDIR="$(mktemp -d)"
trap 'rm -rf "$TDIR"' EXIT

sorted_of() {  # sorted_of <入力> <出力>
    [ -r "$1" ] || fail "読めません: $1"
    sort "$1" > "$2" || fail "並べ替えに失敗しました: $1"
}
sha_of() {  # sha_of <ファイル> → ハッシュ
    local h
    h="$(sha256sum < "$1")" || fail "ハッシュを作れません: $1"
    printf '%s' "${h%% *}"
}
dups_of() {  # dups_of <整列済み> → 重複行（空なら重複なし）
    local d
    d="$(uniq -d < "$1")" || fail "重複の検査に失敗しました: $1"
    printf '%s' "$d"
}
keys_of_results() {  # keys_of_results <結果> <出力（整列済み）>
    [ -r "$1" ] || fail "読めません: $1"
    awk '{print $1}' "$1" > "$TDIR/k.raw" || fail "キーを取り出せません: $1"
    sort "$TDIR/k.raw" > "$2" || fail "キーの並べ替えに失敗しました: $1"
}

meta_of() { printf 'data/thingi10k/%s_results.meta' "$1"; }
res_of()  { printf 'data/thingi10k/%s_results.txt'  "$1"; }
only_of() { printf 'data/thingi10k/%s_only.txt'     "$1"; }
list_of() { printf 'data/thingi10k/%s.txt'          "$1"; }

# ---- 結果の行の検査（起動前と終了時で【同じもの】を使います）----------------
#
# **「FAIL が無い」と「全件が成功」は別です。**
# 書式が崩れた行・途中で切れた行・知らない状態は、どちらでも拒否します。
check_rows() {  # check_rows <ファイル> <require_ok: 0/1> <cols>
    local f="$1" need_ok="$2" cols="$3"
    awk -v need_ok="$need_ok" -v cols="$cols" '
        function isnum(x) {
            return x ~ /^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$/
        }
        {
            # **途中で切れた行**は、数値の列がそろわないことで分かります
            if (NF < cols)                    { print "列が足りない行(" NF "<" cols "): " NR > "/dev/stderr"; bad=1; next }
            if ($1 !~ /^[0-9]+x[0-9]+$/)      { print "キーの形が違う行: " NR > "/dev/stderr"; bad=1; next }
            if ($2 != "ok" && $2 != "FAIL")   { print "知らない状態の行: " NR > "/dev/stderr"; bad=1; next }
            if ($3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/) { print "三角形数が数でない行: " NR > "/dev/stderr"; bad=1; next }
            # **第 5 列は時間**。負や非数は拒否します
            if (!isnum($5) || $5 + 0 < 0)     { print "時間が不正な行: " NR > "/dev/stderr"; bad=1; next }
            if ($6 !~ /^[0-9a-f]{16}$/)       { print "ハッシュの形が違う行: " NR > "/dev/stderr"; bad=1; next }
            # **第 7 列以降 cols 列までは構造の記録**。すべて数でなければなりません
            for (i = 7; i <= cols; i++) {
                if (!isnum($i)) { print "構造の記録が数でない行(" i "列目): " NR > "/dev/stderr"; bad=1; next }
            }
            if (need_ok == 1 && $2 != "ok")   { print "成功していない行: " NR > "/dev/stderr"; bad=1; next }
        }
        END { exit bad ? 1 : 0 }
    ' "$f"
}

# ---- 起動前の検査（ここで止まれば計算は 1 秒も走りません）--------------------
[ -x "$BIN" ] || fail "実行ファイルがありません: $BIN"
[ -f "$MANIFEST" ] || fail "承認済みの基準がありません: $MANIFEST"

# **実際に起動するバイナリと、承認済みのハッシュを直接照合します**
bin_now="$(sha256sum "$BIN" | cut -d' ' -f1)"
bin_want="$(awk -F= '$1=="bin"{print $2}' "$MANIFEST")"
[ -n "$bin_want" ] || fail "基準に bin= がありません: $MANIFEST"
[ "$bin_now" = "$bin_want" ] || fail "起動するバイナリが承認済みの指紋と違います（$BIN）"

# **出力契約は、手で数えずに【バイナリ自身に報告させます】。**
# **版がずれたときに黙って外れるのを防ぎます**（計算は始めません）
# **取得そのものの成功を先に確かめます。**
# **`|| true` で握り潰すと、値が合っていても異常終了した版で本計算を始めます。**
COLERR="$TDIR/cols.err"
bin_cols="$("$BIN" --cols 2>"$COLERR")"
rc_cols=$?
if [ "$rc_cols" != 0 ]; then
    # **stderr も診断の根拠として残します**
    printf '  --cols の stderr:\n' >&2
    head -5 "$COLERR" >&2
    fail "バイナリの --cols が異常終了しました（exit ${rc_cols}）: $BIN"
fi
case "${bin_cols:-x}" in
    ''|*[!0-9]*) fail "バイナリが列数を報告しません（--cols が必要です）: $BIN" ;;
esac

for base in $BASES; do
    only="$(only_of "$base")"; list="$(list_of "$base")"
    [ -f "$list" ] || fail "入力の一覧がありません: $list"
    [ -f "$only" ] || fail "標本の一覧がありません: $only"
    sorted_of "$only" "$TDIR/only.sorted"
    n_now="$(grep -c . "$TDIR/only.sorted")" || n_now=0
    keys_now="$(sha_of "$TDIR/only.sorted")"
    line="$(awk -v b="$base" '$1==b{print}' "$MANIFEST")"
    [ -n "$line" ] || fail "基準に ${base} の行がありません: $MANIFEST"
    keys_want="$(echo "$line" | sed -n 's/.*keys=\([0-9a-f]*\).*/\1/p')"
    n_want="$(echo "$line" | sed -n 's/.* n=\([0-9]*\).*/\1/p')"
    # **`cols` は【固定部の終端列番号】**です（キー・状態・三角形数 2 つ・時間・
    # ハッシュの 6 列 + 構造の記録）。**値の全体の形と範囲を検査します** —
    # 接頭辞だけを見ると `cols=8junk` を 8 として受理してしまいます
    cols_want="$(echo "$line" | sed -n 's/.*cols=\([^ ]*\).*/\1/p')"
    [ -n "$cols_want" ] || fail "基準に ${base} の cols= がありません: $MANIFEST"
    case "$cols_want" in
        ''|*[!0-9]*) fail "${base}: cols が数ではありません（${cols_want}）" ;;
    esac
    [ "$cols_want" -ge 7 ] && [ "$cols_want" -le 4096 ] \
        || fail "${base}: cols が範囲の外です（${cols_want}。7〜4096）"
    [ "$n_now" = "$n_want" ] || fail "${base}: 標本の件数が基準と違います（${n_now} 対 ${n_want}）"
    [ "$keys_now" = "$keys_want" ] || fail "${base}: 標本のハッシュが基準と違います"
    [ "$bin_cols" = "$cols_want" ] \
        || fail "${base}: バイナリの列数（${bin_cols}）が基準の cols（${cols_want}）と違います"
    [ -z "$(dups_of "$TDIR/only.sorted")" ] || fail "${base}: 標本にキーの重複があります"

    res="$(res_of "$base")"; meta="$(meta_of "$base")"
    # **開始時刻は履歴です。一致判定には使いません**
    want="bin=${bin_now} args=${ARGS} keys=${keys_now} n=${n_now} cols=${cols_want}"
    if [ -s "$res" ]; then
        [ "$RESUME" = 1 ] || fail "結果が残っています（新規なら空にしてください）: $res"
        [ -f "$meta" ] || fail "再開なのに meta がありません: $meta"
        [ "$(grep -v '^started=\|^head=\|^b=' "$meta" || true)" = "$want" ] \
            || fail "再開の meta が一致しません: $meta"
        # **共通の列解析で、既存の結果にも【成功】を要求します。**
        # **文字列 ' FAIL ' の検索では、タブ区切りの FAIL がすり抜けます**
        check_rows "$res" 1 "$cols_want" \
            || fail "${base}: 既存の結果に不正な行か、成功していない行があります"
        keys_of_results "$res" "$TDIR/res.sorted"
        [ -z "$(dups_of "$TDIR/res.sorted")" ] \
            || fail "${base}: 既存の結果にキーの重複があります"
        # **再開の結果は、予定キー集合の【部分集合】でなければなりません**
        uniq < "$TDIR/res.sorted" > "$TDIR/res.uniq" || fail "${base}: 既存の結果を読めません"
        uniq < "$TDIR/only.sorted" > "$TDIR/only.uniq" || fail "${base}: 標本を読めません"
        extra="$(comm -23 "$TDIR/res.uniq" "$TDIR/only.uniq")" \
            || fail "${base}: 集合の比較に失敗しました"
        [ -z "$extra" ] || fail "${base}: 既存の結果に、予定に無いキーがあります"
    else
        # **HEAD と b は記録しますが、バイナリのビルド元の証明ではありません**
        printf '%s\nhead=%s\nb=%s\nstarted=%s\n' "$want" \
            "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)" \
            "${KRISITE_COORD_BITS:-21}" "$(date -Is)" > "$meta" \
            || fail "meta を書けません: $meta"
    fi
done

# ---- 実行（計算側と tee 側の状態を直後に両方保存する）-----------------------
for base in $BASES; do
    printf '\n=== %s を回します（%s）===\n' "$base" "$(date -Is)"
    log="data/thingi10k/${base}_run.log"
    set -o pipefail
    "$BIN" "$(list_of "$base")" $ARGS 2>&1 | tee "$log"
    # **★ 直後に配列ごと保存します。** 先に `PIPESTATUS[0]` を代入すると
    # その代入自体が `PIPESTATUS` を書き換えて `[1]` が消えます
    rc=("${PIPESTATUS[@]}")
    set +o pipefail
    [ "${rc[0]:-1}" = 0 ] || fail "${base}: 計算側が異常終了しました（exit ${rc[0]:-?}）"
    [ "${rc[1]:-1}" = 0 ] || fail "${base}: 記録（tee）が失敗しました（exit ${rc[1]:-?}）"

    # ---- 完了の検査（終了コードだけでは判定しません）----
    res="$(res_of "$base")"; only="$(only_of "$base")"
    cols_of_base="$(awk -v b="$base" '$1==b{print}' "$MANIFEST" | sed -n 's/.*cols=\([0-9]*\).*/\1/p')"
    [ -s "$res" ] || fail "${base}: 結果が空です"
    # **読めないことを「行が不正」と混同しない**（理由を取り違えると診断が遠回りになります）
    [ -r "$res" ] || fail "${base}: 結果を読めません: $res"
    check_rows "$res" 1 "$cols_of_base" \
        || fail "${base}: 結果に不正な行か、成功していない行があります"
    sorted_of "$only" "$TDIR/only.sorted"
    keys_of_results "$res" "$TDIR/res.sorted"
    [ -z "$(dups_of "$TDIR/res.sorted")" ] || fail "${base}: 結果にキーの重複"
    # **`cmp` の終了状態は 0（一致）/ 1（不一致）/ 2 以上（エラー）。**
    # **エラーを「一致」と取り違えないよう、3 つを分けます**
    cmp -s "$TDIR/only.sorted" "$TDIR/res.sorted"
    rc_cmp=$?
    case "$rc_cmp" in
        0) ;;
        1) fail "${base}: 予定キー集合と結果のキー集合が一致しません" ;;
        *) fail "${base}: キー集合の比較に失敗しました（exit ${rc_cmp}）" ;;
    esac
    printf '%s: 予定 %s 件すべてが成功しました\n' "$base" "$(grep -c . "$only")"
done
printf '\n完了（%s）\n' "$(date -Is)"
