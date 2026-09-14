#!/bin/bash
# `run_cp23.sh` の自己検査（`SPEC-phase5.md` §5.10.14.74 の 3 / 4）。
#
# **実データは使いません。** 一時領域に模擬の実行ファイルと一覧を作って回します。
# **待機中の実物（`build/thingi_cp1_o3`、`data/thingi10k/*`）には触れません。**
#
#   bash tests/tools/run_cp23_selftest.sh
#
# **対応範囲は sbx の中（Linux / bash）だけ**です。`ctest` には登録していません。
#
# **CP3 が起動していないことは、結果ファイルの不在ではなく
# 【模擬コマンドの起動記録】で確かめます。**
set -u
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN="$SRC/tests/tools/run_cp23.sh"
TMPROOT="$(mktemp -d)"
trap 'chmod -R u+w "$TMPROOT" 2>/dev/null; rm -rf "$TMPROOT"' EXIT
pass=0; fail=0

setup() {  # setup <名前> [各一覧の件数]
    local d="$TMPROOT/$1"; local n="${2:-2}"; rm -rf "$d"
    mkdir -p "$d/build" "$d/data/thingi10k"
    cat > "$d/build/mock" <<'MOCK'
#!/bin/bash
# 模擬の計算側。**起動を記録**し、`<base>_only.txt` を読んで結果行を書く
set -u
cols="${MOCK_COLS:-8}"
# **正しい列数を出しても非零で終わる**構成を作れるようにします
[ "${1:-}" = "--cols" ] && { echo "$cols"; echo "模擬: --cols の異常" >&2; exit "${MOCK_COLS_RC:-0}"; }
list="$1"; base="${list%.txt}"
only="${base}_only.txt"; res="${base}_results.txt"
printf '%s\n' "$(basename "$base")" >> "${MOCK_CALLS:-/dev/null}"
mode="${MOCK_MODE:-ok}"
[ "$mode" = "crash" ] && { echo "模擬: 異常終了"; exit 3; }
while read -r k; do
    [ -z "$k" ] && continue
    if [ "$mode" = "resume" ] && [ -f "$res" ] && grep -q "^${k} " "$res"; then continue; fi
    out="$k"
    [ "$mode" = "otherkey" ] && out="99${k}"
    st=ok
    [ "$mode" = "fail1" ] && { st=FAIL; mode=ok; }
    case "$mode" in
        badrow)  printf '%s\n' "$out" >> "$res" || exit 4 ;;                       # キーだけ
        badtime) printf '%s %s 10 10 -1 0123456789abcdef 1 2\n' "$out" "$st" >> "$res" || exit 4 ;;
        cutrow)  printf '%s %s 10 10 0.1 0123456789abcdef 1\n' "$out" "$st" >> "$res" || exit 4 ;;
        badcol)  printf '%s %s 10 10 0.1 0123456789abcdef 1 x\n' "$out" "$st" >> "$res" || exit 4 ;;
        *)  # **cols 列ぶんの行**（先頭 6 列 + 構造の記録）
            row="$out $st 10 10 0.1 0123456789abcdef"
            i=7; while [ "$i" -le "$cols" ]; do row="$row $i"; i=$((i+1)); done
            [ "$mode" = "withwhy" ] && row="$row 体積の篩"
            printf '%s\n' "$row" >> "$res" || exit 4 ;;
    esac
done < "$only"
[ "$mode" = "lockres" ] && chmod 000 "$res"
echo "模擬: ${base} 完了"
MOCK
    chmod +x "$d/build/mock"
    for b in cp2b cp3; do
        : > "$d/data/thingi10k/$b.txt"
        : > "$d/data/thingi10k/${b}_only.txt"
        for i in $(seq 1 "$n"); do printf '%d0x%d1\n' "$i" "$i" >> "$d/data/thingi10k/${b}_only.txt"; done
    done
    write_manifest "$d"
    echo "$d"
}

write_manifest() {  # 承認済みの基準（いまの一覧から作る = 正常系の基準）
    local d="$1"
    { printf 'bin=%s\n' "$(sha256sum "$d/build/mock" | cut -d' ' -f1)"
      for b in cp2b cp3; do
          printf '%s keys=%s n=%s cols=%s\n' "$b" \
              "$(sort "$d/data/thingi10k/${b}_only.txt" | sha256sum | cut -d' ' -f1)" \
              "$(grep -c . "$d/data/thingi10k/${b}_only.txt")" "${MOCK_COLS:-8}"
      done
    } > "$d/manifest"
}

write_meta() {  # write_meta <dir> <base>
    local d="$1" b="$2"
    printf 'bin=%s args=0 keys=%s n=%s cols=%s\nhead=x\nb=21\nstarted=x\n' \
        "$(sha256sum "$d/build/mock" | cut -d' ' -f1)" \
        "$(sort "$d/data/thingi10k/${b}_only.txt" | sha256sum | cut -d' ' -f1)" \
        "$(grep -c . "$d/data/thingi10k/${b}_only.txt")" "${MOCK_COLS:-8}" \
        > "$d/data/thingi10k/${b}_results.meta"
}

run() {  # run <dir> [引数...]
    ( cd "$1" && KRI_ROOT="$1" KRI_BIN="$1/build/mock" KRI_MANIFEST="$1/manifest" \
        KRI_ARGS="0" KRI_BASES="cp2b cp3" MOCK_CALLS="$1/calls.log" \
        MOCK_COLS="${MOCK_COLS:-8}" MOCK_COLS_RC="${MOCK_COLS_RC:-0}" \
        bash "$RUN" "${@:2}" > "$1/out.txt" 2>&1 )
    echo $?
}

check() {  # check <名前> <期待コード> <実際> <停止理由の語> <dir> <cp3 が起動してよいか>
    local name="$1" want="$2" got="$3" why="$4" d="$5" cp3ok="$6" ok=1 note=""
    [ "$got" = "$want" ] || { ok=0; note="コード"; }
    if [ -n "$why" ] && ! grep -qF -- "$why" "$d/out.txt"; then ok=0; note="${note} 理由"; fi
    # **起動記録**で見ます（結果ファイルの不在では、書けなかった場合と区別できません）
    if [ "$cp3ok" = no ] && grep -q '^cp3$' "$d/calls.log" 2>/dev/null; then
        ok=0; note="${note} CP3 起動"
    fi
    # **起動前に止まるべき構成では、CP2 も起動していないこと**を見ます
    # （`no_run=yes` を渡した構成だけ。計算側が走ってから失敗する構成とは分けます）
    if [ "${7:-}" = "no_run" ] && [ -s "$d/calls.log" ]; then
        ok=0; note="${note} 本計算が起動"
    fi
    if [ "$ok" = 1 ]; then pass=$((pass+1)); else fail=$((fail+1)); fi
    printf '| %s | %s | %s | %s |\n' "$name" "$want" "$got" \
        "$([ "$ok" = 1 ] && echo OK || echo "**NG**${note}")"
}

printf '## `run_cp23.sh` の自己検査（模擬コマンド。実データ不使用）\n\n'
printf '| # 構成 | 期待 | 実際 | |\n|---|---:|---:|---|\n'

d=$(setup c1);  check "1 正常（新規）" 0 "$(run "$d")" "すべてが成功" "$d" yes
d=$(setup c2);  MOCK_MODE=crash; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "2 CP2 が非零終了" 2 "$g" "計算側が異常終了" "$d" no
d=$(setup c3);  rm -f "$d/build/mock"
check "3 起動失敗（実行ファイル無し）" 2 "$(run "$d")" "実行ファイルがありません" "$d" no no_run
# 4a: **結果ファイルだけ**を書けなくする（meta とログは書ける）
d=$(setup c4);  : > "$d/data/thingi10k/cp2b_results.txt"; chmod a-w "$d/data/thingi10k/cp2b_results.txt"
g=$(run "$d"); chmod u+w "$d/data/thingi10k/cp2b_results.txt"
check "4a 記録失敗（結果ファイルだけ）" 2 "$g" "計算側が異常終了" "$d" no
# 4b: tee のログだけを書けなくする
d=$(setup c5);  mkdir -p "$d/data/thingi10k/cp2b_run.log"
check "4b 記録失敗（tee のログ）" 2 "$(run "$d")" "記録（tee）が失敗" "$d" no
d=$(setup c6);  sed -i 's/^bin=.*/bin=deadbeef/' "$d/manifest"
check "5 指紋不一致" 2 "$(run "$d")" "承認済みの指紋と違います" "$d" no no_run
d=$(setup c7);  : > "$d/data/thingi10k/cp2b_only.txt"
check "6a 件数不一致（空の一覧）" 2 "$(run "$d")" "件数が基準と違います" "$d" no no_run
d=$(setup c8);  printf '99x99\n' >> "$d/data/thingi10k/cp2b_only.txt"
check "6b 件数不一致（非空。基準より 1 多い）" 2 "$(run "$d")" "件数が基準と違います" "$d" no no_run
d=$(setup c9);  sed -i '1s/.*/77x77/' "$d/data/thingi10k/cp2b_only.txt"
check "6c 件数は同じでハッシュが違う" 2 "$(run "$d")" "ハッシュが基準と違います" "$d" no no_run
d=$(setup c10); printf '10x11 FAIL 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
check "7 既存 FAIL（再開。空白区切り）" 2 "$(run "$d" --resume)" "成功していない行" "$d" no no_run
d=$(setup c11); MOCK_MODE=otherkey; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "8 同件数だが別キー" 2 "$g" "予定キー集合と結果のキー集合が一致しません" "$d" no
d=$(setup c12); printf '10x11 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
check "9 新規なのに結果が残っている" 2 "$(run "$d")" "結果が残っています" "$d" no no_run
d=$(setup c13); printf '10x11 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b; sed -i '1s/n=2/n=3/' "$d/data/thingi10k/cp2b_results.meta"
check "10 再開の meta が不一致" 2 "$(run "$d" --resume)" "meta が一致しません" "$d" no no_run
d=$(setup c14); printf '10x11 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b; MOCK_MODE=resume; export MOCK_MODE; g=$(run "$d" --resume); unset MOCK_MODE
check "11 正常な再開" 0 "$g" "すべてが成功" "$d" yes
# 12: 書式が崩れた行（キーだけの行）→ 終了時に拒否
d=$(setup c15); MOCK_MODE=badrow; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "12 結果の行の書式が崩れている" 2 "$g" "不正な行か、成功していない行" "$d" no
# 13: 全件回ったが 1 件 FAIL → 終了時に拒否
d=$(setup c16); MOCK_MODE=fail1; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "13 予定は揃うが 1 件 FAIL" 2 "$g" "不正な行か、成功していない行" "$d" no
# 14: 再開の結果に、予定に無いキーが混ざっている
d=$(setup c17); printf '55x55 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
check "14 再開の結果に予定外のキー" 2 "$(run "$d" --resume)" "予定に無いキー" "$d" no no_run

# 15: 第 5 列（時間）が不正
d=$(setup c18); MOCK_MODE=badtime; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "15 時間の列が不正（負）" 2 "$g" "不正な行か、成功していない行" "$d" no
# 16: 構造の記録が途中で切れている
d=$(setup c19); MOCK_MODE=cutrow; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "16 構造の記録が途中で切れている" 2 "$g" "不正な行か、成功していない行" "$d" no
# 17: 構造の記録に数でない値
d=$(setup c20); MOCK_MODE=badcol; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "17 構造の記録が数でない" 2 "$g" "不正な行か、成功していない行" "$d" no
# 18: **タブ区切りの FAIL** が、再開の検査をすり抜けないこと
d=$(setup c21)
printf '10x11\tFAIL\t10\t10\t0.1\t0123456789abcdef\t1\t2\n' \
    > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
check "18 タブ区切りの FAIL（再開）" 2 "$(run "$d" --resume)" "成功していない行" "$d" no no_run
# 19: 基準に cols が無い
d=$(setup c22); sed -i 's/ cols=8//' "$d/manifest"
check "19 基準に cols が無い" 2 "$(run "$d")" "cols= がありません" "$d" no no_run

# ---- 現行の固定部 184 列での構成（`--cols` が報告する値）----
export MOCK_COLS=184
d=$(setup c23); check "20 正常（184 列）" 0 "$(run "$d")" "すべてが成功" "$d" yes
d=$(setup c24); MOCK_MODE=withwhy; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "21 理由の語つき（184 列 + 1）" 0 "$g" "すべてが成功" "$d" yes
d=$(setup c25); MOCK_MODE=cutrow; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "22 途中欠落（184 列に足りない）" 2 "$g" "不正な行か、成功していない行" "$d" no
unset MOCK_COLS
# 23: バイナリの列数が基準と違う
d=$(setup c26); sed -i 's/cols=8/cols=9/g' "$d/manifest"
check "23 バイナリの列数が基準と違う" 2 "$(run "$d")" "が基準の cols" "$d" no no_run
# 24: cols が数でない（接頭辞だけの受理を防ぐ）
d=$(setup c27); sed -i 's/cols=8/cols=8junk/g' "$d/manifest"
check "24 cols が数でない" 2 "$(run "$d")" "cols が数ではありません" "$d" no no_run
# 25: cols が範囲の外
d=$(setup c28); sed -i 's/cols=8/cols=3/g' "$d/manifest"
check "25 cols が範囲の外" 2 "$(run "$d")" "cols が範囲の外" "$d" no no_run

# 26: **正しい列数を出すが、--cols が非零で終わる** → 起動前に拒否
d=$(setup c29); MOCK_COLS_RC=7; export MOCK_COLS_RC; g=$(run "$d"); unset MOCK_COLS_RC
check "26 --cols が非零終了（値は正しい）" 2 "$g" "--cols が異常終了" "$d" no no_run

# ---- 前処理の失敗を「正常」と取り違えないこと ----
# 27: 標本を読めない → 起動前に止まり、本計算は走らない
d=$(setup c30); chmod 000 "$d/data/thingi10k/cp2b_only.txt"
g=$(run "$d"); chmod 644 "$d/data/thingi10k/cp2b_only.txt"
check "27 標本を読めない（起動前）" 2 "$g" "読めません" "$d" no no_run
# 28: 結果を読めない → 完了の検査で止まる（「差分なし」と取り違えない）
d=$(setup c31); MOCK_MODE=lockres; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
chmod 644 "$d/data/thingi10k/"*_results.txt 2>/dev/null
check "28 結果を読めない（完了の検査）" 2 "$g" "読めません" "$d" no

# ---- 前処理のコマンドを【正常に出力してから非零で終わる】形に差し替える ----
#
# **読み取り可否の検査では止まらない経路**を突きます
# （`SPEC-phase5.md` §5.10.14.74 の指摘 3）。
inject() {  # inject <dir> <コマンド名> [引数に含まれる語。省略すると毎回]
    local d="$1" cmd="$2" pat="${3:-}"
    mkdir -p "$d/binover"
    cat > "$d/binover/$cmd" <<INJ
#!/bin/bash
# **正しい出力を出してから非零で終わります。**
# **引数の語で選べます** — どの呼び出しを壊したかを試験ごとに限定するため
real="\$(PATH=/usr/bin:/bin command -v $cmd)"
"\$real" "\$@"
rc=\$?
pat='$pat'
if [ -z "\$pat" ] || printf '%s' "\$*" | grep -qF -- "\$pat"; then
    echo "模擬: $cmd を非零で終わらせます（引数: \$*）" >&2
    exit 7
fi
exit "\$rc"
INJ
    chmod +x "$d/binover/$cmd"
}
run_inj() {  # run_inj <dir> [引数...]
    ( cd "$1" && PATH="$1/binover:$PATH" KRI_ROOT="$1" KRI_BIN="$1/build/mock" \
        KRI_MANIFEST="$1/manifest" KRI_ARGS="0" KRI_BASES="cp2b cp3" \
        MOCK_CALLS="$1/calls.log" MOCK_COLS="${MOCK_COLS:-8}" \
        bash "$RUN" "${@:2}" > "$1/out.txt" 2>&1 )
    echo $?
}

d=$(setup c32); inject "$d" sort
check "29 sort が正常出力の後に非零" 2 "$(run_inj "$d")" "並べ替えに失敗" "$d" no no_run
d=$(setup c33); inject "$d" uniq
check "30 uniq が正常出力の後に非零" 2 "$(run_inj "$d")" "重複検査に失敗" "$d" no no_run
d=$(setup c34); inject "$d" sha256sum
check "31 sha256sum が正常出力の後に非零" 2 "$(run_inj "$d")" "指紋を作れません" "$d" no no_run
d=$(setup c35); inject "$d" cmp
check "32 cmp が非零（集合の比較）" 2 "$(run_inj "$d")" "キー集合の比較に失敗" "$d" no
# 33: **キーの取り出しだけ**を失敗させる（`check_rows` の awk は壊さない）
d=$(setup c36); inject "$d" awk '{print $1}'
check "33 キーの取り出しだけ非零" 2 "$(run_inj "$d")" "キーを取り出せません" "$d" no
# 34: **基準の解析だけ**を失敗させる（起動前に止まり、CP2 も起動しない）
d=$(setup c37); inject "$d" awk '$1=="bin"'
check "34 基準の解析だけ非零" 2 "$(run_inj "$d")" "基準を読めません" "$d" no no_run
# 35: **CP 別の基準の解析だけ**を失敗させる
d=$(setup c38); inject "$d" awk '$1==b'
check "35 CP 別の基準の解析だけ非零" 2 "$(run_inj "$d")" "基準を読めません" "$d" no no_run

chk2() { if [ "$2" = "$3" ]; then pass=$((pass+1)); r=OK; else fail=$((fail+1)); r='**NG**'; fi
         printf '| %s | %s | %s | %s |\n' "$1" "$2" "$3" "$r"; }

# **両方の CP について、meta・結果・実行ログが作られていないこと**
untouched() {  # untouched <dir> → 作られていたファイル数
    local d="$1" n=0
    for b in cp2b cp3; do
        for f in _results.meta _results.txt _run.log; do
            [ -e "$d/data/thingi10k/${b}${f}" ] && n=$((n+1))
        done
    done
    echo "$n"
}

# 36: **--check-only は照合だけで、計算も meta の作成もしない**
d=$(setup c39); g=$(run "$d" --check-only)
check "36 --check-only（照合のみ）" 0 "$g" "照合だけ行いました" "$d" no
chk2 "36b そのとき本計算は起動していない" 0 "$([ -s "$d/calls.log" ] && echo 1 || echo 0)"
chk2 "36c 両 CP の meta・結果・ログが作られていない" 0 "$(untouched "$d")"

# 37〜41: **誤記・併用・余分な引数・重複**
d=$(setup c40); g=$(run "$d" --check-onyl)
check "37 誤記（--check-onyl）" 2 "$g" "知らない引数" "$d" no no_run
chk2 "37b 何も作られていない" 0 "$(untouched "$d")"
d=$(setup c41); g=$(run "$d" --check-only --resume)
check "38 併用（--check-only --resume）" 0 "$g" "照合だけ行いました" "$d" no
chk2 "38b 何も作られていない" 0 "$(untouched "$d")"
d=$(setup c42); g=$(run "$d" --resume --check-only)
check "39 併用（逆順）" 0 "$g" "照合だけ行いました" "$d" no
chk2 "39b 何も作られていない" 0 "$(untouched "$d")"
d=$(setup c43); g=$(run "$d" --check-only extra)
check "40 余分な引数" 2 "$g" "知らない引数" "$d" no no_run
chk2 "40b 何も作られていない" 0 "$(untouched "$d")"
d=$(setup c44); g=$(run "$d" --resume --resume)
check "41 同じ旗の重複" 2 "$g" "引数が重複" "$d" no no_run
chk2 "41b 何も作られていない" 0 "$(untouched "$d")"

printf '\n**OK %d / NG %d**\n' "$pass" "$fail"
[ "$fail" = 0 ]
