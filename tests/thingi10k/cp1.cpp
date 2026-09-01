// Krisite — Phase 5 CP1: 実データでの正しさ（`SPEC-phase5.md` §3、§9）
//
// **拒否 / 失敗 / 停止を分けて数えます**（§1）。混ぜると数字が意味を失います。
//
//   拒否   入口の検査で弾いた           正しい動作。**理由を分けて記録**
//   失敗   受け入れたのに壊れた         **1 件でも欠陥**
//   停止   KRISITE_CHECK で止まった     失敗と同じ扱い
//
// **選択はメタデータ、分類は量子化後**（§1.0）。集合は `fetch.py` が決めたものを
// そのまま使い、ここでは選び直しません。
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/par/thread_pool.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

namespace {

/// 拒否の理由（§1 の分類表）。**性質の違うものを混ぜないこと。**
enum class Reject {
    None,
    BoundaryOriginal,    ///< 元から ∂S ≠ 0（開曲面を含む）
    OpenSurface,         ///< 境界辺がある（∂S ≠ 0 の一種。理由を分ける）
    BoundaryAfterQuant,  ///< **量子化 + 面積 0 の除去のあとで** ∂S ≠ 0。b の設計に直結
    OutOfRange,          ///< 座標範囲
    Empty,               ///< 量子化で三角形が消えた
};

const char* reject_name(Reject r) {
    switch (r) {
        case Reject::BoundaryOriginal:
            return "∂S≠0（元から）";
        case Reject::OpenSurface:
            return "開曲面（境界辺あり）";
        case Reject::BoundaryAfterQuant:
            return "∂S≠0（量子化+除去のあと）";
        case Reject::OutOfRange:
            return "座標範囲";
        case Reject::Empty:
            return "空（量子化で消えた）";
        default:
            return "—";
    }
}

/// 1 つの入力を量子化して受け入れ判定まで行う。
struct Prepared {
    mesh::TriMesh mesh;
    Reject reject = Reject::None;
    std::size_t dropped = 0;       ///< 面積 0 で落とした三角形
    std::size_t merged = 0;        ///< 同じ格子点に落ちた頂点
    bool boundary_before = false;  ///< 除去の【前】に ∂S = 0 だったか
};

Prepared prepare(const krithingi::RawMesh& raw, std::uint64_t seed) {
    Prepared p;
    const krithingi::Quantized q = krithingi::quantize(raw, krithingi::make_transform(seed));
    p.dropped = q.dropped_degenerate;
    p.merged = q.merged_vertices;
    if (q.out_of_range) {
        p.reject = Reject::OutOfRange;
        return p;
    }
    if (q.mesh.triangles.empty()) {
        p.reject = Reject::Empty;
        return p;
    }
    p.mesh = q.mesh;
    // **面積 0 の三角形は既に落としてあります**（支持平面が定義できないため）。
    // **除去が ∂S を壊すことがあるので、除去後に検査します**（§1 の分類表）
    p.boundary_before = (p.dropped == 0);
    if (!mesh::boundary_is_zero(p.mesh)) {
        p.reject = (p.dropped > 0) ? Reject::BoundaryAfterQuant : Reject::BoundaryOriginal;
    }
    return p;
}

struct Counts {
    std::size_t pairs = 0, accepted = 0, failed = 0, halted = 0;
    std::size_t reject[6] = {0, 0, 0, 0, 0, 0};
    std::size_t dropped_total = 0, merged_total = 0;
    std::size_t models_with_dropped = 0;
};

/// 出力のバイト単位のハッシュ。**同一版の中の比較に使います**（`CLAUDE.md`）。
///
/// 索引（`ray_index`）は厳密な絞り込みなので、**有無で 1 ビットも変わってはいけません。**
unsigned long long hash_mesh(const csg::SoupMesh& m) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long v) {
        h ^= v;
        h *= 1099511628211ull;
    };
    mix(m.vertices.size());
    mix(m.triangles.size());
    for (const geom::HPointD& v : m.vertices) {
        for (std::size_t l = 0; l < geom::limbs::kHomoXyz; ++l) {
            mix(v.x[l]);
            mix(v.y[l]);
            mix(v.z[l]);
        }
        for (std::size_t l = 0; l < geom::limbs::kHomoW; ++l) mix(v.w[l]);
    }
    for (const mesh::Tri& t : m.triangles) {
        mix(t[0]);
        mix(t[1]);
        mix(t[2]);
    }
    return h;
}

/// §3.1 の検査。**解析的期待値は使えない**ので恒等式と位相で見ます。
bool check_one(const mesh::TriMesh& a, const mesh::TriMesh& b, const csg::BoolOptions& o,
               par::ThreadPool* pool, std::string* why, unsigned long long* hash_out = nullptr) {
    const csg::PolySoup A = csg::from_mesh(a), B = csg::from_mesh(b);
    csg::ToMeshOptions tm;
    tm.split_contacts = true;
    tm.threads = pool->size();
    tm.pool = pool;
    const csg::SoupMesh mu = csg::to_mesh(csg::boolean(A, B, csg::BoolOp::Union, o), tm);
    const csg::SoupMesh mi = csg::to_mesh(csg::boolean(A, B, csg::BoolOp::Intersection, o), tm);
    const csg::SoupMesh md = csg::to_mesh(csg::boolean(A, B, csg::BoolOp::Difference, o), tm);

    for (const auto* pr : {&mu, &mi, &md}) {
        const mesh::TopologyReport t = mesh::check_topology(pr->triangles);
        if (t.empty) continue;
        if (!t.edge_manifold) {
            *why = "辺多様体でない";
            return false;
        }
        if (!t.vertex_manifold) {
            *why = "頂点多様体でない";
            return false;
        }
        if (!t.oriented) {
            *why = "向きが整合しない";
            return false;
        }
        if (!t.no_degenerate) {
            *why = "退化三角形が残った";
            return false;
        }
        if ((t.chi % 2) != 0) {
            *why = "χ が奇数";
            return false;
        }
    }
    if (hash_out != nullptr) {
        *hash_out = hash_mesh(mu) ^ (hash_mesh(mi) * 3) ^ (hash_mesh(md) * 7);
    }
    // **体積の恒等式**（§3.1）。|A∪B| + |A∩B| = |A| + |B| を出力側の 6 倍体積で見る
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = (argc > 1) ? argv[1] : "data/thingi10k/cp1.txt";
    const std::size_t limit = (argc > 2) ? std::strtoul(argv[2], nullptr, 10) : 0;
    const unsigned depth = (argc > 3) ? static_cast<unsigned>(std::atoi(argv[3])) : 6;
    const unsigned nthreads = (argc > 4) ? static_cast<unsigned>(std::atoi(argv[4])) : 16;
    // **索引の ON/OFF でハッシュが一致するかを対ごとに確かめるモード**（CP1.5D）。
    const bool verify_index = (argc > 5) && (std::atoi(argv[5]) != 0);
    // **済みの対をやり直すモード。** 既定は追記（再開）
    const bool redo = (argc > 6) && (std::atoi(argv[6]) != 0);
    std::size_t verified = 0;
    // 資源上限で落ちた対（1 行 1 対）。**手で足すのではなく、監視スクリプトが足します**
    std::vector<std::string> skip;
    {
        std::ifstream f("data/thingi10k/cp1_skip.txt");
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty()) skip.push_back(line.substr(0, line.find(' ')));
        }
    }

    std::vector<std::string> ids;
    {
        std::ifstream f(list);
        std::string id, nf;
        while (f >> id >> nf) ids.push_back(id);
    }
    if (limit != 0 && ids.size() > limit) ids.resize(limit);
    std::printf("CP1: %zu 件、深度 %u、b=%d\n", ids.size(), depth, KRISITE_COORD_BITS);

    // ---- 1. 量子化と受け入れ判定（**モデル単位**）----
    std::vector<Prepared> prep(ids.size());
    Counts c;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const krithingi::RawMesh raw =
            krithingi::load_kmesh("data/thingi10k/kmesh/" + ids[i] + ".kmesh");
        prep[i] = prepare(raw, 1000 + i);
        c.dropped_total += prep[i].dropped;
        c.merged_total += prep[i].merged;
        if (prep[i].dropped > 0) ++c.models_with_dropped;
        ++c.reject[static_cast<int>(prep[i].reject)];
        if ((i + 1) % 200 == 0) std::printf("  量子化 %zu / %zu\n", i + 1, ids.size());
    }
    std::printf("\n## 量子化（b=%d）\n\n", KRISITE_COORD_BITS);
    std::printf("| 事象 | 件数 |\n|---|---:|\n");
    std::printf("| 受け入れ | %zu |\n", c.reject[0]);
    for (int r = 1; r < 6; ++r) {
        if (c.reject[r] != 0)
            std::printf("| 拒否: %s | %zu |\n", reject_name(static_cast<Reject>(r)), c.reject[r]);
    }
    std::printf("| 面積 0 で落とした三角形（延べ） | %zu |\n", c.dropped_total);
    std::printf("| 落ちたモデル数 | %zu |\n", c.models_with_dropped);
    std::printf("| 同じ格子点に落ちた頂点（延べ） | %zu |\n", c.merged_total);

    // ---- 2. 対を作ってブール演算（§2.1）----
    //
    // **面数の小さい順に組みます。** 大きい対に時間を取られて、
    // **coverage が出る前に打ち切られるのを避けるため**です。
    // **並びは決定的**なので、再開しても同じ対になります。
    csg::BoolOptions o;
    o.depth = depth;
    o.adaptive = true;
    o.leaf_threshold = 0;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    par::ThreadPool pool(nthreads);
    o.threads = nthreads;
    o.pool = &pool;

    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (prep[i].reject == Reject::None) order.push_back(i);
    }
    std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
        const std::size_t nx = prep[x].mesh.triangles.size(), ny = prep[y].mesh.triangles.size();
        return (nx != ny) ? (nx < ny) : (ids[x] < ids[y]);
    });

    // **結果は 1 対ごとに追記します。** 途中で止まっても、
    // どこまで通ったかが残り、再開できます（§3.3 の追跡に要る）
    // **やり直しモードは別のファイルに書きます。** 追記すると再開の記録が濁ります
    const std::string done_path =
        redo ? "data/thingi10k/cp1_verify.txt" : "data/thingi10k/cp1_results.txt";
    std::vector<std::string> already;
    {
        std::ifstream f("data/thingi10k/cp1_results.txt");
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty()) already.push_back(line.substr(0, line.find(' ')));
        }
    }
    const auto seen = [&already](const std::string& k) {
        return std::find(already.begin(), already.end(), k) != already.end();
    };
    std::ofstream out(done_path, std::ios::app);

    std::printf("\n## ブール演算（対 %zu、スレッド %u）\n\n", order.size() / 2, nthreads);
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t k = 0; k + 1 < order.size(); k += 2) {
        const std::size_t i = order[k], j = order[k + 1];
        const std::string key = ids[i] + "x" + ids[j];
        if (!redo && seen(key)) continue;
        // **資源上限で落ちた対を飛ばします**（`SPEC-phase5.md` §1 の「停止」）。
        // **落ちた対を記録しないと、再開のたびに同じ対で落ち続けます。**
        if (std::find(skip.begin(), skip.end(), key) != skip.end()) {
            ++c.halted;
            std::printf("  停止済みとして飛ばす %s\n", key.c_str());
            continue;
        }
        ++c.pairs;
        // **開始を先に出します。** OOM で殺されると完了行が出ないので、
        // **どの対で落ちたかが分からなくなります**（実際に分からなくなりました）。
        std::printf("  → 開始 %s（入力 %zu+%zu）\n", key.c_str(), prep[i].mesh.triangles.size(),
                    prep[j].mesh.triangles.size());
        std::fflush(stdout);
        const auto tp = std::chrono::steady_clock::now();
        std::string why;
        unsigned long long h = 0;
        const bool ok = check_one(prep[i].mesh, prep[j].mesh, o, &pool, &why, &h);
        // **索引の有無で出力が変わらないこと**（`SPEC-phase5.md` §6.3 の安全弁）。
        // 索引は厳密な絞り込みなので、**バイトで一致しなければ欠陥**です。
        if (ok && verify_index) {
            csg::BoolOptions o2 = o;
            o2.ray_index = !o.ray_index;
            unsigned long long h2 = 0;
            std::string why2;
            const bool ok2 = check_one(prep[i].mesh, prep[j].mesh, o2, &pool, &why2, &h2);
            if (!ok2 || h2 != h) {
                ++c.failed;
                why = "索引の有無で出力が変わった";
                std::printf("**失敗** %s: %s（%016llx vs %016llx）\n", key.c_str(), why.c_str(), h,
                            h2);
                out << key << " FAIL " << prep[i].mesh.triangles.size() << ' '
                    << prep[j].mesh.triangles.size() << " 0 " << why << '\n';
                out.flush();
                continue;
            }
            ++verified;
        }
        const double dt =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        if (ok) {
            ++c.accepted;
        } else {
            ++c.failed;
            std::printf("**失敗** %s: %s\n", key.c_str(), why.c_str());
        }
        out << key << ' ' << (ok ? "ok" : "FAIL") << ' ' << prep[i].mesh.triangles.size() << ' '
            << prep[j].mesh.triangles.size() << ' ' << dt << ' ';
        std::array<char, 24> hb{};
        std::snprintf(hb.data(), hb.size(), "%016llx", h);
        out << hb.data() << ' ' << why << '\n';
        out.flush();
        const double s =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("  %zu 対目 %s（入力 %zu+%zu、%.1f s、累計 %.0f s）\n", c.pairs, key.c_str(),
                    prep[i].mesh.triangles.size(), prep[j].mesh.triangles.size(), dt, s);
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("\n対 %zu / 成功 %zu / **失敗 %zu** / **停止 %zu** / %.1f s\n", c.pairs, c.accepted,
                c.failed, c.halted, s);
    if (verify_index) {
        std::printf("**索引の有無でハッシュ一致: %zu 対**\n", verified);
    }
    return c.failed == 0 ? 0 : 1;
}
