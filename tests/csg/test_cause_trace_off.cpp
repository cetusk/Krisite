// **診断定義 OFF の単独 include。** 通常ビルドで記録が完全に消えることを確かめます。
// **他ヘッダの間接 include に依存しません**（このファイルは 2 本しか include しません）。
#include <cstdio>

#include "krisite/csg/diagnostic_trace.hpp"
#include "thingi10k/cause_trace_io.hpp"

int main() {
    std::printf("=== 診断定義 OFF の単独 include ===\n");
    std::printf("KRI_DIAG_TRACE_ON = %d\n", KRI_DIAG_TRACE_ON);
    int calls = 0;
    // **呼出しが完全に消える**ので、この式は 1 度も評価されません。
    KRI_DIAG_PUSH((void*)nullptr, frags, (++calls, 0));
    const bool ok = (KRI_DIAG_TRACE_ON == 0) && (calls == 0);
    std::printf("  %s OFF では記録の呼出しが消える（評価回数 %d / 期待 0）\n",
                ok ? "ok  " : "**NG**", calls);
    std::printf("判定: %s\n", ok ? "OFF の単独 include は通過" : "**失敗あり**");
    return ok ? 0 : 1;
}
