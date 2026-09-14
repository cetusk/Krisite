#!/bin/bash
# 判定器そのものの検査（模擬ログ）。**判定器が空回りしないことを先に確かめます。**
set -u
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
pass=0; fail=0

mk() { cat > "$T/$1"; }

# (a) 正しい診断: 読み出しのスタックに to_mesh.hpp がある
mk a.log <<'LOG'
=================================================================
==12==ERROR: AddressSanitizer: heap-use-after-free on address 0x603000000010
READ of size 8 at 0x603000000010 thread T0
    #0 0x55 in krisite::csg::detail::repair_unresolved_edges(...) to_mesh.hpp:431
    #1 0x55 in krisite::csg::to_mesh(...) to_mesh.hpp:470
0x603000000010 is located 0 bytes inside of 48-byte region
freed by thread T0 here:
    #0 0x7f in operator delete(void*) asan_new_delete.cpp:160
    #1 0x55 in std::vector<...>::push_back(...) new_allocator.h:168
previously allocated by thread T0 here:
    #0 0x7f in operator new(unsigned long) asan_new_delete.cpp:95
SUMMARY: AddressSanitizer: heap-use-after-free to_mesh.hpp:431
LOG

# (b) 解放側だけに to_mesh.hpp がある（読み出しは別の場所）→ 拒否
mk b.log <<'LOG'
==12==ERROR: AddressSanitizer: heap-use-after-free on address 0x603000000010
READ of size 8 at 0x603000000010 thread T0
    #0 0x55 in somewhere_else(...) other_file.hpp:99
0x603000000010 is located 0 bytes inside of 48-byte region
freed by thread T0 here:
    #0 0x55 in krisite::csg::detail::repair_unresolved_edges(...) to_mesh.hpp:402
SUMMARY: AddressSanitizer: heap-use-after-free other_file.hpp:99
LOG

# (c) 別の診断（stack-buffer-overflow）→ 拒否
mk c.log <<'LOG'
==12==ERROR: AddressSanitizer: stack-buffer-overflow on address 0x7ffd
READ of size 8 at 0x7ffd thread T0
    #0 0x55 in krisite::csg::detail::repair_unresolved_edges(...) to_mesh.hpp:431
SUMMARY: AddressSanitizer: stack-buffer-overflow to_mesh.hpp:431
LOG

# (d) 診断が無い（正常終了のログ）→ 拒否
mk d.log <<'LOG'
## 接触辺の細分（案 H）
NEGCTL realloc=1 cap_before=104 cap_after=208 added=2 site=to_mesh.hpp
**不一致 0 件**
LOG

# (e) 確保側だけに to_mesh.hpp → 拒否
mk e.log <<'LOG'
==12==ERROR: AddressSanitizer: heap-use-after-free on address 0x603000000010
READ of size 8 at 0x603000000010 thread T0
    #0 0x55 in somewhere_else(...) other_file.hpp:99
previously allocated by thread T0 here:
    #0 0x55 in krisite::csg::detail::repair_unresolved_edges(...) to_mesh.hpp:365
LOG

chk() {  # chk <名前> <ログ> <期待 0/1>
    bash "$D/negctl_check.sh" "$T/$2" to_mesh.hpp; got=$?
    if [ "$got" = "$3" ]; then pass=$((pass+1)); r=OK; else fail=$((fail+1)); r='**NG**'; fi
    printf '| %s | %s | %s | %s |\n' "$1" "$3" "$got" "$r"
}
printf '## 負の対照の判定器の検査（模擬ログ）\n\n'
printf '| 構成 | 期待 | 実際 | |\n|---|---:|---:|---|\n'
chk "(a) 読み出しのスタックに一致" a.log 0
chk "(b) 解放側だけに一致" b.log 1
chk "(c) 別の診断（種類が違う）" c.log 1
chk "(d) 診断が無い" d.log 1
chk "(e) 確保側だけに一致" e.log 1
printf '\n**OK %d / NG %d**\n' "$pass" "$fail"
[ "$fail" = 0 ]
