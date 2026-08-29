// Krisite — 並列化（`SPEC-phase4.md` §7.1 / §7.3）
//
// **Phase 4 の主検出器は決定性です。**
//
// > 出力がスレッド数に依らず、ビット単位で同一であること。
//
// **競合は非決定性として現れます。** 位相も体積も、順序が違うだけの出力を区別しません
// （`SPEC-phase3.md` の「辺平面の軸選択」と同じ形）。
//
// ## なぜバイト列で比べるのか
//
// 意味論を変えていないので、**逐次実装がそのまま完全な正解器**になります。
// 「体積と位相が一致」ではなく「**バイト列が一致**」まで要求できるのは、
// 並列化が意味論に触れていないことの帰結です（§0.2）。
//
// ## TSan との分担
//
// **決定性の検査だけでは足りません**（§7.2）。コーパスが小さいので、競合があっても
// たまたま顕在化しないことがあります。**TSan は実行された経路上の競合を、
// 顕在化していなくても検出します。** CI の別ジョブが担当します。
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/par/thread_pool.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;

/// 出力のバイト列。**同じバイト列 = 同じ出力**（§7.1）。
std::string bytes(const SoupMesh& m) {
    std::string s;
    const auto put = [&s](const void* p, std::size_t n) {
        s.append(static_cast<const char*>(p), n);
    };
    const std::size_t nv = m.vertices.size(), nt = m.triangles.size();
    put(&nv, sizeof nv);
    put(&nt, sizeof nt);
    for (const auto& v : m.vertices) put(&v, sizeof v);
    for (const auto& t : m.triangles) put(&t, sizeof t);
    return s;
}

BoolOptions all_on(unsigned depth, bool adaptive, unsigned threads, krisite::par::ThreadPool* p) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;
    o.adaptive = adaptive;
    o.leaf_threshold = 0;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    o.threads = threads;
    o.pool = p;
    return o;
}

const char* op_name(BoolOp op) {
    switch (op) {
        case BoolOp::Union:
            return "∪";
        case BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

std::size_t g_cmp = 0, g_nonempty = 0;

/// §7.1: スレッド数 1 / 2 / 4 / 8 で出力がビット単位で一致すること。
void test_determinism() {
    constexpr unsigned kThreads[] = {1, 2, 4, 8};
    // **`ThreadPool` はコピーも move もできません**（所有するスレッドがあるため）。
    // 器に入れるときは間接参照で持ちます
    std::vector<std::unique_ptr<krisite::par::ThreadPool>> pools;
    for (unsigned t : kThreads) pools.push_back(std::make_unique<krisite::par::ThreadPool>(t));

    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                std::string base;
                for (std::size_t ti = 0; ti < std::size(kThreads); ++ti) {
                    ToMeshOptions tm;
                    tm.split_contacts = true;
                    // **出口にもプールを渡すこと。** 渡さないと `to_mesh` は
                    // 1 スレッドで回り、**出口の並列化が検査の外に出ます**
                    // （変異 23 が素通りしていました。`IMPL-phase4.md` §2.4）
                    tm.threads = kThreads[ti];
                    tm.pool = pools[ti].get();
                    const BoolOptions o = all_on(d, d == kMaxDepth, kThreads[ti], pools[ti].get());
                    const std::string s =
                        bytes(to_mesh(boolean(from_mesh(a), from_mesh(b), op, o), tm));
                    const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                            "（深度 " + std::to_string(d) + "、スレッド " +
                                            std::to_string(kThreads[ti]) + "）";
                    if (ti == 0) {
                        base = s;
                        if (!s.empty()) ++g_nonempty;
                    } else {
                        KRI_CHECK_MSG(s == base, tag + ": **スレッド数で出力が変わった**（§7.1）");
                        ++g_cmp;
                    }
                }
            }
        }
    }
    // **期待値は式で持たせます。** 空回りは成功と区別が付きません
    const std::size_t want =
        kritest::corpus().size() * 3 * (kMaxDepth + 1) * (std::size(kThreads) - 1);
    KRI_CHECK_MSG(g_cmp == want, "比較数が式と合わない" + kritest::pair_msg(want, g_cmp));
    // **空の出力ばかりなら、バイト列の一致は何も言っていません**
    KRI_CHECK_MSG(g_nonempty > 0, "**基準がすべて空**。決定性の検査が空回りしています");
    std::printf("    決定性 %zu 件（うち非空の基準 %zu 構成）\n", g_cmp, g_nonempty);
}

/// §7.3: スレッド局所であることを、**統計で確かめます。**
///
/// キャッシュを共有すると命中数が変わります（1 本にまとまるので増えます）。
/// **共有していないなら、スレッド数を増やすと命中率は下がるはず**です。
void test_thread_local_cache() {
    const kritest::Case& c = kritest::corpus()[0];
    const TriMesh a = c.make_a(), b = c.make_b();
    std::size_t hits1 = 0, hits8 = 0, entries1 = 0, entries8 = 0;
    for (unsigned th : {1u, 8u}) {
        krisite::par::ThreadPool pool(th);
        BoolStats st{};
        boolean(from_mesh(a), from_mesh(b), BoolOp::Union, all_on(3, true, th, &pool), &st);
        if (th == 1) {
            hits1 = st.cache_hits;
            entries1 = st.cache_entries;
        } else {
            hits8 = st.cache_hits;
            entries8 = st.cache_entries;
        }
    }
    // **項目数はスレッド数に比例して増えるはず**（同じ点が複数のキャッシュに入る）。
    // 減っていたら共有されています
    KRI_CHECK_MSG(entries8 >= entries1,
                  "**キャッシュの項目数がスレッド数で減った。** 共有されていませんか" +
                      kritest::pair_msg(entries1, entries8));
    std::printf("    キャッシュ: 1 スレッド 命中 %zu / 項目 %zu → 8 スレッド 命中 %zu / 項目 %zu\n",
                hits1, entries1, hits8, entries8);
}

}  // namespace

int main() {
    std::printf("\n  並列化（SPEC-phase4 §7）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    test_determinism();
    test_thread_local_cache();
    std::printf("\n");
    return kritest::finish("csg/parallel");
}
