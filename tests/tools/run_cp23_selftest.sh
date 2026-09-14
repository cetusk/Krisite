#!/bin/bash
# `run_cp23.sh` の自己検査（`SPEC-phase5.md` §5.10.14.74 の 3 / 4）。
#
# **実データは使いません。** 一時領域に模擬の実行ファイルと一覧を作って回します。
# **待機中の実物（`build/thingi_cp1_o3`、`data/thingi10k/*`）には触れません。**
#
#   bash tests/tools/run_cp23_selftest.sh
#
# **対応範囲は sbx の中（Linux / bash）だけ**です。`ctest` には登録していません。
set -u
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN="$SRC/tests/tools/run_cp23.sh"
TMPROOT="$(mktemp -d)"
trap 'rm -rf "$TMPROOT"' EXIT
pass=0; fail=0

# 一時領域に、実物と同じ形の木を作る
setup() {
    local d="$TMPROOT/$1"; rm -rf "$d"; mkdir -p "$d/build" "$d/data/thingi10k"
    cat > "$d/build/mock" <<'MOCK'
#!/bin/bash
# 模擬の計算側。`<base>_only.txt` を読み、結果行を書く
set -u
list="$1"; base="${list%.txt}"
only="${base}_only.txt"; res="${base}_results.txt"
mode="${MOCK_MODE:-ok}"
[ "$mode" = "crash" ] && { echo "模擬: 異常終了"; exit 3; }
while read -r k; do
    [ -z "$k" ] && continue
    if [ "$mode" = "resume" ] && [ -f "$res" ] && grep -q "^${k} " "$res"; then continue; fi
    if [ "$mode" = "otherkey" ]; then k="XX${k}"; fi
    st=ok
    if [ "$mode" = "fail1" ]; then st=FAIL; mode=ok; fi
    printf '%s %s 10 10 0.1 0000\n' "$k" "$st" >> "$res" || exit 4
done < "$only"
echo "模擬: ${base} 完了"
MOCK
    chmod +x "$d/build/mock"
    for b in cp2b cp3; do
        : > "$d/data/thingi10k/$b.txt"
        printf 'a%sxb%s\nc%sxd%s\n' "$b" "$b" "$b" "$b" > "$d/data/thingi10k/${b}_only.txt"
    done
    sha256sum "$d/build/mock" > "$d/sha256"
    echo "$d"
}

# 再開の試験用に、正しい meta を書く
write_meta() {  # write_meta <dir> <base>
    local d="$1" b="$2"
    printf '%s args=0 keys=%s n=%s\nstarted=x\n' \
        "$(sha256sum "$d/build/mock" | cut -d' ' -f1)" \
        "$(sort "$d/data/thingi10k/${b}_only.txt" | sha256sum | cut -d' ' -f1)" \
        "$(grep -c . "$d/data/thingi10k/${b}_only.txt")" \
        > "$d/data/thingi10k/${b}_results.meta"
}

run() {  # run <dir> [引数...] → 終了コードと出力
    ( cd "$1" && KRI_ROOT="$1" KRI_BIN="$1/build/mock" KRI_SHA="$1/sha256" \
        KRI_ARGS="0" KRI_BASES="cp2b cp3" bash "$RUN" "${@:2}" > "$1/out.txt" 2>&1 )
    echo $?
}

check() {  # check <名前> <期待コード> <実際> <期待する停止理由の語> <dir> <CP3 が走ってよいか>
    local name="$1" want="$2" got="$3" why="$4" d="$5" cp3ok="$6" ok=1
    [ "$got" = "$want" ] || ok=0
    if [ -n "$why" ] && ! grep -q "$why" "$d/out.txt"; then ok=0; fi
    if [ "$cp3ok" = no ] && [ -s "$d/data/thingi10k/cp3_results.txt" ]; then ok=0; fi
    if [ "$ok" = 1 ]; then pass=$((pass+1)); else fail=$((fail+1)); fi
    printf '| %s | %s | %s | %s |\n' "$name" "$want" "$got" \
        "$([ "$ok" = 1 ] && echo OK || echo '**NG**')"
}

printf '## `run_cp23.sh` の自己検査（模擬コマンド。実データ不使用）\n\n'
printf '| # 構成 | 期待コード | 実際 | |\n|---|---:|---:|---|\n'

d=$(setup c1); got=$(run "$d");            check "1 正常（新規）" 0 "$got" "すべてが成功" "$d" yes
d=$(setup c2); MOCK_MODE=crash; export MOCK_MODE
got=$(run "$d"); unset MOCK_MODE;          check "2 CP2 が非零終了" 2 "$got" "計算側が異常終了" "$d" no
d=$(setup c3); rm -f "$d/build/mock"
got=$(run "$d");                           check "3 起動失敗（実行ファイル無し）" 2 "$got" "実行ファイルがありません" "$d" no
d=$(setup c4); chmod a-w "$d/data/thingi10k"
got=$(run "$d"); chmod u+w "$d/data/thingi10k"
                                           check "4a 記録失敗（結果ファイル）" 2 "$got" "" "$d" no
d=$(setup c5); mkdir -p "$d/data/thingi10k/cp2b_run.log"
got=$(run "$d");                           check "4b 記録失敗（tee のログ）" 2 "$got" "記録（tee）が失敗" "$d" no
d=$(setup c6); echo "deadbeef  $d/build/mock" > "$d/sha256"
got=$(run "$d");                           check "5 指紋不一致" 2 "$got" "指紋が一致しません" "$d" no
d=$(setup c7); : > "$d/data/thingi10k/cp2b_only.txt"
got=$(run "$d");                           check "6 件数不一致（一覧が空）" 2 "$got" "一覧が空です" "$d" no
d=$(setup c8); printf 'acp2bxbcp2b FAIL 10 10 0.1 0000\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
got=$(run "$d" --resume);                  check "7 既存 FAIL（再開）" 2 "$got" "既存の FAIL" "$d" no
d=$(setup c9); MOCK_MODE=otherkey; export MOCK_MODE
got=$(run "$d"); unset MOCK_MODE;          check "8 同件数だが別キー" 2 "$got" "予定キー集合と結果のキー集合が一致しません" "$d" no
# 9: 新規なのに結果が残っている → 拒否する
d=$(setup c10); printf 'acp2bxbcp2b ok 10 10 0.1 0000\n' > "$d/data/thingi10k/cp2b_results.txt"
got=$(run "$d");                           check "9 新規なのに結果が残っている" 2 "$got" "結果が残っています" "$d" no
# 10: 再開の meta が一致しない → 拒否する
d=$(setup c11); printf 'acp2bxbcp2b ok 10 10 0.1 0000\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b; sed -i 's/n=2/n=3/' "$d/data/thingi10k/cp2b_results.meta"
got=$(run "$d" --resume);                  check "10 再開の meta が不一致" 2 "$got" "meta が一致しません" "$d" no
# 11: 正常な部分結果からの再開（模擬は済みキーを飛ばす）
d=$(setup c12); printf 'acp2bxbcp2b ok 10 10 0.1 0000\n' > "$d/data/thingi10k/cp2b_results.txt"
write_meta "$d" cp2b
MOCK_MODE=resume; export MOCK_MODE
got=$(run "$d" --resume); unset MOCK_MODE
check "11 正常な再開" 0 "$got" "すべてが成功" "$d" yes

printf '\n**OK %d / NG %d**\n' "$pass" "$fail"
[ "$fail" = 0 ]
