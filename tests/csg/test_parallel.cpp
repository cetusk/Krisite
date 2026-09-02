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
#include <cstring>
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
    for (unsigned t : kThreads) {
        pools.push_back(std::make_unique<krisite::par::ThreadPool>(t));
        // **ディスパッチの下限を外します**（`SPEC-phase4.md` §6.3）。
        //
        // 下限は**性能のための機構**で、項目数が少ない段を逐次に落とします。
        // **そのままだと、このコーパスでは扇が常に逐次になり**（頂点は最大 254、
        // 下限は 256）、**並列の経路が検査されません。** 実際に変異 23 の検出器が
        // 消えます。**「後段で埋める機構は、上流の誤りを覆い隠します」**（`CLAUDE.md`）。
        pools.back()->set_min_items(0);
    }

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

std::size_t g_cmp_soup = 0, g_cmp_nary = 0;

/// §12 の CP3:「**全コーパス**」。主コーパスだけでは通らない経路が 2 つあります。
///
/// | コーパス | ここでしか通らない経路 |
/// |---|---|
/// | `soup_only_corpus` | **凸分割**（`FromMeshOptions::max_merges`）。非凸な面 |
/// | `nary_corpus` | **深さ 3 以上の指示関数**。連鎖した `Indicator` の評価 |
///
/// **主コーパスは二項・三角形化だけです。** 機構を全部有効にすると言うなら、
/// この 2 つを外せません。
void test_determinism_all_corpora() {
    constexpr unsigned kThreads[] = {1, 2, 4, 8};
    std::vector<std::unique_ptr<krisite::par::ThreadPool>> pools;
    for (unsigned t : kThreads) {
        pools.push_back(std::make_unique<krisite::par::ThreadPool>(t));
        // **ディスパッチの下限を外します**（`SPEC-phase4.md` §6.3）。
        //
        // 下限は**性能のための機構**で、項目数が少ない段を逐次に落とします。
        // **そのままだと、このコーパスでは扇が常に逐次になり**（頂点は最大 254、
        // 下限は 256）、**並列の経路が検査されません。** 実際に変異 23 の検出器が
        // 消えます。**「後段で埋める機構は、上流の誤りを覆い隠します」**（`CLAUDE.md`）。
        pools.back()->set_min_items(0);
    }

    // ---- スープ専用（凸分割）----
    FromMeshOptions fm;
    fm.max_merges = kMergeAll;
    for (const kritest::Case& c : kritest::soup_only_corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                std::string base;
                for (std::size_t ti = 0; ti < std::size(kThreads); ++ti) {
                    ToMeshOptions tm;
                    tm.split_contacts = true;
                    tm.threads = kThreads[ti];
                    tm.pool = pools[ti].get();
                    const BoolOptions o = all_on(d, d == kMaxDepth, kThreads[ti], pools[ti].get());
                    const std::string s =
                        bytes(to_mesh(boolean(from_mesh(a, fm), from_mesh(b, fm), op, o), tm));
                    const std::string tag = std::string("スープ ") + c.id + " " + op_name(op) +
                                            "（深度 " + std::to_string(d) + "、スレッド " +
                                            std::to_string(kThreads[ti]) + "）";
                    if (ti == 0) {
                        base = s;
                    } else {
                        KRI_CHECK_MSG(s == base, tag + ": **スレッド数で出力が変わった**（§7.1）");
                        ++g_cmp_soup;
                    }
                }
            }
        }
    }

    // ---- n 項（連鎖した指示関数）----
    //
    // **和で連鎖します。** 各ケースの意味論を再現するのは `test_nary` の仕事で、
    // ここが問うのは「**同じ式がスレッド数に依らず同じ答えを返すか**」だけです。
    for (const kritest::NaryCase& c : kritest::nary_corpus()) {
        const std::vector<TriMesh> m = c.make();
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            std::string base;
            for (std::size_t ti = 0; ti < std::size(kThreads); ++ti) {
                ToMeshOptions tm;
                tm.split_contacts = true;
                tm.threads = kThreads[ti];
                tm.pool = pools[ti].get();
                const BoolOptions o = all_on(d, d == kMaxDepth, kThreads[ti], pools[ti].get());
                PolySoup acc = from_mesh(m[0]);
                for (std::size_t i = 1; i < m.size(); ++i) {
                    acc = boolean(acc, from_mesh(m[i]), BoolOp::Union, o);
                }
                const std::string s = bytes(to_mesh(acc, tm));
                const std::string tag = std::string("n 項 ") + c.id + "（深度 " +
                                        std::to_string(d) + "、スレッド " +
                                        std::to_string(kThreads[ti]) + "）";
                if (ti == 0) {
                    base = s;
                } else {
                    KRI_CHECK_MSG(s == base, tag + ": **スレッド数で出力が変わった**（§7.1）");
                    ++g_cmp_nary;
                }
            }
        }
    }

    const std::size_t want_soup =
        kritest::soup_only_corpus().size() * 3 * (kMaxDepth + 1) * (std::size(kThreads) - 1);
    const std::size_t want_nary =
        kritest::nary_corpus().size() * (kMaxDepth + 1) * (std::size(kThreads) - 1);
    KRI_CHECK_MSG(g_cmp_soup == want_soup,
                  "スープ経路の比較数が式と合わない" + kritest::pair_msg(want_soup, g_cmp_soup));
    KRI_CHECK_MSG(g_cmp_nary == want_nary,
                  "n 項の比較数が式と合わない" + kritest::pair_msg(want_nary, g_cmp_nary));
    std::printf("    決定性（凸分割）%zu 件、（n 項）%zu 件\n", g_cmp_soup, g_cmp_nary);
}

/// §7.3: スレッド局所であることを、**統計で確かめます。**
///
/// キャッシュを共有すると命中数が変わります（1 本にまとまるので増えます）。
/// **共有していないなら、スレッド数を増やすと命中率は下がるはず**です。
void test_thread_local_cache() {
    const kritest::Case& c = kritest::corpus()[0];
    const TriMesh a = c.make_a(), b = c.make_b();
    std::size_t hits1 = 0, hits8 = 0, entries1 = 0, entries8 = 0;
    // **空回りを許さないこと。** `entries8 == entries1` は「共有していない」ではなく
    // 「**8 スレッドで走っていない**」（呼び出し元が全部さらった）ときにも起きます。
    // 等号で通すと、検査が黙って無意味になります（`CLAUDE.md`）。
    //
    // スケジューリングは保証できないので、**複数回試して 1 度も並列に走らなければ
    // 落とします。**
    //
    // **回数はコア数への依存を実測して決めました**（`IMPL-phase4.md` §6.1）。
    //
    //   4 コア   1 回目（GitHub Actions の標準ランナーはここ）
    //   2 コア   1 回目
    //   1 コア   2〜16 回目（呼び出し元が全部さらってしまう）
    //
    // 1 回が約 2 ms なので、32 回でも 70 ms 程度です。
    constexpr int kTries = 32;
    int tries = 0;
    for (; tries < kTries; ++tries) {
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
        if (entries8 > entries1) break;
    }
    KRI_CHECK_MSG(tries < kTries,
                  "**8 スレッドで一度も並列に走りませんでした**（32 回試行）。"
                  "この検査は空回りしています（§7.3）。項目数が増えないなら、"
                  "共有の有無を見ていません。**実行環境のコア数を確認してください**");
    std::printf(
        "    キャッシュ: 1 スレッド 命中 %zu / 項目 %zu → 8 スレッド 命中 %zu / 項目 %zu"
        "（%d 回目で並列を確認）\n",
        hits1, entries1, hits8, entries8, tries + 1);
}

/// §4.4: **同値な構成点の代表が正準であること。**
///
/// `lex_less` は同値な点に順序を付けないので、**何もしないと「どれが残るか」が
/// 入力の並び次第**になります。**多角形の並びを逆にして、頂点の列が
/// バイト単位で変わらないこと**で確かめます。
///
/// > **幾何も位相も体積も、表現の違いを見ません。** §7.1 の決定性検査も、
/// > 同じコードどうしの比較なので素通りします。**この検査だけが捕まえます。**
///
/// **到達点は「構成点の集合が同じ限り一意」です**（§4.4）。幾何だけの関数では
/// ありません。GCD 正規化には除算が要ります。
void test_canonical_representative() {
    std::size_t n = 0, swaps = 0, merged = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                const BoolOptions o = all_on(d, d == kMaxDepth, 1u, nullptr);
                PolySoup r = boolean(from_mesh(a), from_mesh(b), op, o);
                ToMeshOptions tm;
                tm.split_contacts = true;
                ToMeshStats ts{};
                const SoupMesh m1 = to_mesh(r, tm, &ts);
                swaps += ts.canonical_swaps;
                merged += ts.merged_by_value;
                // **並びだけを変えます。** 断片の集合は同じです
                std::reverse(r.polys.begin(), r.polys.end());
                const SoupMesh m2 = to_mesh(r, tm);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                        "（深度 " + std::to_string(d) + "）";
                KRI_CHECK_MSG(m1.vertices.size() == m2.vertices.size(),
                              tag + ": 並べ替えで頂点数が変わった");
                bool eq = true;
                for (std::size_t i = 0; i < m1.vertices.size() && eq; ++i) {
                    eq = std::memcmp(&m1.vertices[i], &m2.vertices[i], sizeof m1.vertices[i]) == 0;
                }
                KRI_CHECK_MSG(eq, tag +
                                      ": **多角形の並びで同次座標の表現が変わった**"
                                      "（§4.4 の代表が正準でない）");
                ++n;
            }
        }
    }
    // **空回りを許さないこと。** 代表が一度も入れ替わらないなら、この検査は
    // 「正準化が効いている」ではなく「正準化する場面が無い」を見ています
    KRI_CHECK_MSG(swaps > 0,
                  "**代表が一度も入れ替わっていません。** §4.4 の機構が空回りしています");
    KRI_CHECK_MSG(merged > 0, "**値で併合した点が 0 件。** 同値な組がありません");
    std::printf("    正準な代表 %zu 構成（代表の入れ替え %zu 回 / 値の併合 %zu 回）\n", n, swaps,
                merged);
}

}  // namespace

/// §7.5 の変異 23（出口の段のバリア）を**確定的に**検出させるための入力。
///
/// **変異 23 は「ID が完了順に配られる」ので、分裂頂点が 2 個以上あって、
/// かつ扇が複数スレッドで動いたときにしか現れません。**
///
/// **コーパスではその機会が 1,008 実行中 72 件（7.1%）しかありませんでした。**
/// そのため検出が確率的になり、**実測で 55〜60%** でした
/// （`IMPL-phase5.md` §24.5 / §27）。CI では実質コイン投げです。
///
/// ここでは**辺で接する 2 箱の組を 16 個並べ**、分裂頂点を **62〜124 個**作ります。
/// 機会が十分にあれば、完了順がスロット順と一致する確率は無視できます。
///
/// **コーパスには入れません。** コーパスに足すと全テストの費用が上がり、
/// **無関係な検査が負荷で落ちます**（`IMPL-phase5.md` §24）。
krisite::mesh::TriMesh contact_chain(int k) {
    krisite::mesh::TriMesh m;
    const int w = 32 / k;
    for (int i = 0; i < k; ++i) {
        const int x0 = -32 + i * w, xm = x0 + w / 2, x1 = x0 + w;
        m = kritest::concat(
            m, kritest::box(kritest::at(x0, 32), kritest::at(-1, 2), kritest::at(-1, 2),
                            kritest::at(xm, 32), 0, kritest::at(1, 2)));
        m = kritest::concat(
            m, kritest::box(kritest::at(xm, 32), 0, kritest::at(-1, 2), kritest::at(x1, 32),
                            kritest::at(1, 2), kritest::at(1, 2)));
    }
    return m;
}

/// §7.1 を、**出口の順序が実際に効く入力**で見ます。
void test_exit_ordering() {
    std::printf("  §7.5: 出口の段の順序（変異 23 の検出器）\n");
    constexpr unsigned kThreads[] = {1, 2, 4, 8};
    std::vector<std::unique_ptr<krisite::par::ThreadPool>> pools;
    for (unsigned t : kThreads) {
        pools.push_back(std::make_unique<krisite::par::ThreadPool>(t));
        pools.back()->set_min_items(0);
    }
    const TriMesh a = contact_chain(16);
    const TriMesh b = kritest::box(kritest::at(-1, 2), kritest::at(-1, 2), kritest::at(1, 4),
                                   kritest::at(1, 2), kritest::at(1, 2), kritest::at(1, 2));

    // **機会を数えます**（`CLAUDE.md`「発火した」と「効いた」を分ける）。
    std::size_t chances = 0, max_split = 0, cmp = 0;
    // **番人自身がスケジューラに依存します。**
    //
    // > `fan_threads_used > 1` は「**複数のワーカーが実際に扇の仕事を取った**」ことで、
    // > **取るかどうかは走行時の競合次第**です。機械が空いていれば 1 本が全部さらい、
    // > 混んでいれば散ります。**単独実行で 23〜24、`ctest -j4` の競合下で 1** まで振れました。
    //
    // **回数を増やして解きます。** 足りなければもう一巡します。
    // **「一度も散らない」なら本物の後退**（機構が壊れている）なので、そこでは落ちます。
    // **通常は 1 巡で足ります**（24 ≫ 12）ので、実行時間はほぼ変わりません。
    int rounds = 0;
    for (; rounds < 3 && chances < 12; ++rounds)
    // **∩ と ＼ が分裂頂点を作ります**（∪ は 0）。実測で選びました
    for (BoolOp op : {BoolOp::Intersection, BoolOp::Difference}) {
        for (unsigned d = 0; d <= 3; ++d) {
            std::string base;
            for (std::size_t ti = 0; ti < std::size(kThreads); ++ti) {
                ToMeshOptions tm;
                tm.split_contacts = true;
                tm.threads = kThreads[ti];
                tm.pool = pools[ti].get();
                ToMeshStats ts;
                const BoolOptions o = all_on(d, d == 3, kThreads[ti], pools[ti].get());
                const std::string sbytes =
                    bytes(to_mesh(boolean(from_mesh(a), from_mesh(b), op, o), tm, &ts));
                max_split = std::max(max_split, ts.split.split_vertices);
                if (kThreads[ti] > 1 && ts.split.split_vertices >= 2 &&
                    ts.split.fan_threads_used > 1) {
                    ++chances;
                }
                const std::string tag = std::string("接触鎖 ") + op_name(op) + "（深度 " +
                                        std::to_string(d) + "、スレッド " +
                                        std::to_string(kThreads[ti]) + "）";
                if (ti == 0) {
                    base = sbytes;
                } else {
                    KRI_CHECK_MSG(sbytes == base, tag + ": **スレッド数で出力が変わった**（§7.1）");
                    ++cmp;
                }
            }
        }
    }
    // **式で持ちます**（実測で書くと、比較が減っても PASS になります）。
    // **巡回数も式に入れます** — 入れないと、もう一巡した瞬間に空回りします
    KRI_CHECK_MSG(cmp == static_cast<std::size_t>(rounds) * 2 * 4 * 3,
                  "比較数が式と合わない" +
                      kritest::pair_msg(static_cast<std::size_t>(rounds) * 2 * 4 * 3, cmp));
    // ★ **機構が動いたことを数えます**（案 A）。
    // **0 なら、この検査は変異 23 について何も言っていません。**
    KRI_CHECK_MSG(chances >= 12,
                  "**出口の順序が効く機会が足りない**（分裂頂点 2 以上 かつ 扇が複数スレッド）"
                  "。変異 23 の検出が確率的に戻ります" +
                      kritest::pair_msg(12, chances));
    std::printf("    分裂頂点 最大 %zu / 順序が効く機会 %zu 件 / 比較 %zu 件 / **%d 巡**\n",
                max_split, chances, cmp, rounds);
}

int main() {
    std::printf("\n  並列化（SPEC-phase4 §7）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    test_determinism();
    test_determinism_all_corpora();
    test_canonical_representative();
    test_thread_local_cache();
    test_exit_ordering();
    std::printf("\n");
    return kritest::finish("csg/parallel");
}
