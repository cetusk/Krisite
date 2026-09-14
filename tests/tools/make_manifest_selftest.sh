#!/bin/bash
# `make_manifest.sh` の自己検査（一時領域の模擬環境。**実物に触れません**）。
#
#   bash tests/tools/make_manifest_selftest.sh
set -u
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MK="$SRC/tests/tools/make_manifest.sh"
RUN="$SRC/tests/tools/run_cp23.sh"
T="$(mktemp -d)"; trap 'chmod -R u+w "$T" 2>/dev/null; rm -rf "$T"' EXIT
pass=0; fail=0

setup() {
    local d="$T/$1"; rm -rf "$d"; mkdir -p "$d/build" "$d/data/thingi10k" "$d/tests/tools"
    cat > "$d/build/cand" <<'MOCK'
#!/bin/bash
set -u
[ "${1:-}" = "--cols" ] && { echo "${MOCK_COLS:-184}"; exit "${MOCK_COLS_RC:-0}"; }
list="$1"; base="${list%.txt}"
printf '%s\n' "$(basename "$base")" >> "${MOCK_CALLS:-/dev/null}"
cols="${MOCK_COLS:-184}"
while read -r k; do
    [ -z "$k" ] && continue
    row="$k ok 10 10 0.1 0123456789abcdef"; i=7
    while [ "$i" -le "$cols" ]; do row="$row $i"; i=$((i+1)); done
    printf '%s\n' "$row" >> "${base}_results.txt" || exit 4
done < "${base}_only.txt"
MOCK
    chmod +x "$d/build/cand"
    cp "$d/build/cand" "$d/build/old"        # 旧バイナリ（別物にするため 1 行足す）
    printf '# old\n' >> "$d/build/old"
    for b in cp2b cp3; do
        : > "$d/data/thingi10k/$b.txt"
        printf '10x11\n20x21\n' > "$d/data/thingi10k/${b}_only.txt"
    done
    echo "$d"
}
chk() {  # chk <名前> <期待> <実際>
    if [ "$2" = "$3" ]; then pass=$((pass+1)); r=OK; else fail=$((fail+1)); r='**NG**'; fi
    printf '| %s | %s | %s | %s |\n' "$1" "$2" "$3" "$r"
}

printf '## `make_manifest.sh` の自己検査（模擬環境。実物に触れません）\n\n'
printf '| 構成 | 期待 | 実際 | |\n|---|---:|---:|---|\n'

d=$(setup m1); ( cd "$d" && bash "$MK" build/cand out.manifest > log 2>&1 ); chk "1 正常" 0 $?
d=$(setup m2); ( cd "$d" && bash "$MK" build/cand out.manifest > /dev/null 2>&1
                 bash "$MK" build/cand out.manifest > log 2>&1 ); chk "2 既存を上書きしない" 2 $?
d=$(setup m3); ( cd "$d" && bash "$MK" build/cand out.manifest > /dev/null 2>&1
                 bash "$MK" --force build/cand out.manifest > log 2>&1 ); chk "3 --force なら置換" 0 $?
d=$(setup m4); ( cd "$d" && MOCK_COLS_RC=7 bash "$MK" build/cand out.manifest > log 2>&1 )
chk "4 --cols が非零終了" 2 $?
d=$(setup m5); ( cd "$d" && MOCK_COLS=3 bash "$MK" build/cand out.manifest > log 2>&1 )
chk "5 --cols が範囲の外" 2 $?
d=$(setup m6); ( cd "$d" && bash "$MK" build/none out.manifest > log 2>&1 )
chk "6 実行ファイルが無い" 2 $?
d=$(setup m7); : > "$d/data/thingi10k/cp2b_only.txt"
( cd "$d" && bash "$MK" build/cand out.manifest > log 2>&1 ); chk "7 標本が空" 2 $?
d=$(setup m8); printf '10x11\n10x11\n' > "$d/data/thingi10k/cp2b_only.txt"
( cd "$d" && bash "$MK" build/cand out.manifest > log 2>&1 ); chk "8 標本に重複" 2 $?

# ---- 9: 作った基準で、候補を KRI_BIN で名指しして回せること（移動しない）----
d=$(setup m9)
( cd "$d" && bash "$MK" build/cand m.manifest > /dev/null 2>&1 ) || true
( cd "$d" && KRI_ROOT="$d" KRI_BIN="$d/build/cand" KRI_MANIFEST="$d/m.manifest" \
    KRI_ARGS="0" KRI_BASES="cp2b cp3" MOCK_CALLS="$d/calls.log" MOCK_COLS=184 \
    bash "$RUN" > "$d/run.log" 2>&1 )
chk "9 候補を KRI_BIN で名指しして回す" 0 $?
# ---- 10: 旧バイナリを名指しすると、指紋が合わず起動前に止まる ----
d=$(setup m10)
( cd "$d" && bash "$MK" build/cand m.manifest > /dev/null 2>&1 ) || true
( cd "$d" && KRI_ROOT="$d" KRI_BIN="$d/build/old" KRI_MANIFEST="$d/m.manifest" \
    KRI_ARGS="0" KRI_BASES="cp2b cp3" MOCK_CALLS="$d/calls.log" MOCK_COLS=184 \
    bash "$RUN" > "$d/run.log" 2>&1 )
rc=$?; started=$([ -s "$d/calls.log" ] && echo 1 || echo 0)
chk "10 旧バイナリを名指し → 起動前に拒否" 2 "$rc"
chk "10b そのとき本計算は起動していない" 0 "$started"

# ---- 11: 標本を読めない → 失敗し、既存の基準は残る ----
d=$(setup m11)
( cd "$d" && bash "$MK" build/cand keep.manifest > /dev/null 2>&1 ) || true
before="$(sha256sum "$d/keep.manifest" | cut -d' ' -f1)"
chmod 000 "$d/data/thingi10k/cp2b_only.txt"
( cd "$d" && bash "$MK" --force build/cand keep.manifest > log 2>&1 ); rc=$?
chmod 644 "$d/data/thingi10k/cp2b_only.txt"
after="$(sha256sum "$d/keep.manifest" | cut -d' ' -f1)"
chk "11 標本を読めない" 2 "$rc"
chk "11b そのとき既存の基準は変わらない" "$before" "$after"

# ---- 前処理のコマンドを【正常に出力してから非零で終わる】形に差し替える ----
inject() {  # inject <dir> <コマンド名>
    local d="$1" cmd="$2"
    mkdir -p "$d/binover"
    cat > "$d/binover/$cmd" <<INJ
#!/bin/bash
real="\$(PATH=/usr/bin:/bin command -v $cmd)"
"\$real" "\$@"
echo "模擬: $cmd を非零で終わらせます" >&2
exit 7
INJ
    chmod +x "$d/binover/$cmd"
}
inj_mk() {  # inj_mk <dir> <出力先>
    ( cd "$1" && PATH="$1/binover:$PATH" bash "$MK" --force build/cand "$2" > log 2>&1 )
    echo $?
}
for c in sort uniq sha256sum; do
    d=$(setup "mi_$c")
    ( cd "$d" && bash "$MK" build/cand keep.manifest > /dev/null 2>&1 ) || true
    before="$(sha256sum "$d/keep.manifest" | cut -d' ' -f1)"
    inject "$d" "$c"
    rc=$(inj_mk "$d" keep.manifest)
    after="$(sha256sum "$d/keep.manifest" | cut -d' ' -f1)"
    chk "12 $c が正常出力の後に非零" 2 "$rc"
    chk "12b そのとき既存の基準は変わらない（$c）" "$before" "$after"
done

printf '\n**OK %d / NG %d**\n' "$pass" "$fail"
[ "$fail" = 0 ]
