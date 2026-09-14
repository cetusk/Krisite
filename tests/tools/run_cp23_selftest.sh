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
[ "${1:-}" = "--cols" ] && { echo "$cols"; exit 0; }
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
        MOCK_COLS="${MOCK_COLS:-8}" \
        bash "$RUN" "${@:2}" > "$1/out.txt" 2>&1 )
    echo $?
}

check() {  # check <名前> <期待コード> <実際> <停止理由の語> <dir> <cp3 が起動してよいか>
    local name="$1" want="$2" got="$3" why="$4" d="$5" cp3ok="$6" ok=1 note=""
    [ "$got" = "$want" ] || { ok=0; note="コード"; }
    if [ -n "$why" ] && ! grep -q "$why" "$d/out.txt"; then ok=0; note="${note} 理由"; fi
    # **起動記録**で見ます（結果ファイルの不在では、書けなかった場合と区別できません）
    if [ "$cp3ok" = no ] && grep -q '^cp3$' "$d/calls.log" 2>/dev/null; then
        ok=0; note="${note} CP3 起動"
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
check "3 起動失敗（実行ファイル無し）" 2 "$(run "$d")" "実行ファイルがありません" "$d" no
# 4a: **結果ファイルだけ**を書けなくする（meta とログは書ける）
d=$(setup c4);  : > "$d/data/thingi10k/cp2b_results.txt"; chmod a-w "$d/data/thingi10k/cp2b_results.txt"
g=$(run "$d"); chmod u+w "$d/data/thingi10k/cp2b_results.txt"
check "4a 記録失敗（結果ファイルだけ）" 2 "$g" "計算側が異常終了" "$d" no
# 4b: tee のログだけを書けなくする
d=$(setup c5);  mkdir -p "$d/data/thingi10k/cp2b_run.log"
check "4b 記録失敗（tee のログ）" 2 "$(run "$d")" "記録（tee）が失敗" "$d" no
d=$(setup c6);  sed -i 's/^bin=.*/bin=deadbeef/' "$d/manifest"
check "5 指紋不一致" 2 "$(run "$d")" "承認済みの指紋と違います" "$d" no
d=$(setup c7);  : > "$d/data/thingi10k/cp2b_only.txt"
check "6a 件数不一致（空の一覧）" 2 "$(run "$d")" "件数が基準と違います" "$d" no
d=$(setup c8);  printf '99x99\n' >> "$d/data/thingi10k/cp2b_only.txt"
check "6b 件数不一致（非空。基準より 1 多い）" 2 "$(run "$d")" "件数が基準と違います" "$d" no
d=$(setup c9);  sed -i '1s/.*/77x77/' "$d/data/thingi10k/cp2b_only.txt"
check "6c 件数は同じでハッシュが違う" 2 "$(run "$d")" "ハッシュが基準と違います" "$d" no
d=$(setup c10); printf '10x11 FAIL 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
check "7 既存 FAIL（再開。空白区切り）" 2 "$(run "$d" --resume)" "成功していない行" "$d" no
d=$(setup c11); MOCK_MODE=otherkey; export MOCK_MODE; g=$(run "$d"); unset MOCK_MODE
check "8 同件数だが別キー" 2 "$g" "予定キー集合と結果のキー集合が一致しません" "$d" no
d=$(setup c12); printf '10x11 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
check "9 新規なのに結果が残っている" 2 "$(run "$d")" "結果が残っています" "$d" no
d=$(setup c13); printf '10x11 ok 10 10 0.1 0123456789abcdef 1 2\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b; sed -i '1s/n=2/n=3/' "$d/data/thingi10k/cp2b_results.meta"
check "10 再開の meta が不一致" 2 "$(run "$d" --resume)" "meta が一致しません" "$d" no
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
check "14 再開の結果に予定外のキー" 2 "$(run "$d" --resume)" "予定に無いキー" "$d" no

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
check "18 タブ区切りの FAIL（再開）" 2 "$(run "$d" --resume)" "成功していない行" "$d" no
# 19: 基準に cols が無い
d=$(setup c22); sed -i 's/ cols=8//' "$d/manifest"
check "19 基準に cols が無い" 2 "$(run "$d")" "cols= がありません" "$d" no

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
check "23 バイナリの列数が基準と違う" 2 "$(run "$d")" "が基準の cols" "$d" no
# 24: cols が数でない（接頭辞だけの受理を防ぐ）
d=$(setup c27); sed -i 's/cols=8/cols=8junk/g' "$d/manifest"
check "24 cols が数でない" 2 "$(run "$d")" "cols が数ではありません" "$d" no
# 25: cols が範囲の外
d=$(setup c28); sed -i 's/cols=8/cols=3/g' "$d/manifest"
check "25 cols が範囲の外" 2 "$(run "$d")" "cols が範囲の外" "$d" no

printf '\n**OK %d / NG %d**\n' "$pass" "$fail"
[ "$fail" = 0 ]
