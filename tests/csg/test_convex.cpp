// Krisite — 凸分割（`SPEC-phase3.md` §4.1）
//
// **三角形化してから、凸性が保たれる限り貪欲に併合する。** 併合の上限を 0 にすれば
// 三角形化になるので、コードは一本化されています。**切り替えは実行時**です（§4.1.2）。
//
// ## 何が一致し、何が一致しないか（§4.1.1）
//
// | 量 | |
// |---|---|
// | 体積、$(C, \chi)$ | **厳密に一致** |
// | 三角形の集合、断片数、平面数 | **一致しない**（凸分割のほうが少ない。それが目的） |
//
// **一致しないものを明示しておかないと、正しい実装が落ちます。**
// ここでは「一致しないこと」も**検査します** — 一致してしまったら凸分割が
// 効いていないので、機構が空回りしています。
//
// ## 決定性は幾何の検査では捕まりません
//
// 貪欲な併合は走査順に依存します。**正準な順序**を決めないと同じ入力で違う結果が
// 出ますが、**位相も体積も「幾何として正しいか」しか見ません**
// （`CLAUDE.md`「幾何が正しくても決定性が壊れることがあります」）。
// そこで**同じ入力を 2 度通して完全一致する**ことを別に検査します。
#include <algorithm>
#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;
constexpr std::size_t kPasses = kMaxDepth + 2;  ///< 固定深度 0〜3 + 適応分割

BoolOptions all_on(unsigned depth, bool adaptive) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;
    o.adaptive = adaptive;
    o.leaf_threshold = 0;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    return o;
}

FromMeshOptions tri_opt() {
    FromMeshOptions o;
    o.max_merges = 0;  // **明示します。** 既定に依存させない
    return o;
}
FromMeshOptions cvx_opt() {
    FromMeshOptions o;
    o.max_merges = kMergeAll;
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

/// 断片を「支持平面 ID + 辺平面 ID の列」で表した正準形（決定性の検査に使う）。
std::vector<std::vector<PlaneId>> soup_key(const PolySoup& s) {
    std::vector<std::vector<PlaneId>> out;
    out.reserve(s.polys.size());
    for (const Poly& q : s.polys) {
        std::vector<PlaneId> v;
        v.push_back(q.frag.support);
        v.push_back(q.frag.flipped ? 1u : 0u);
        v.insert(v.end(), q.frag.edge.begin(), q.frag.edge.end());
        out.push_back(std::move(v));
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::size_t g_cmp = 0, g_differ_tris = 0;

/// 断面が任意の**星型**単純多角形の角柱（頂点 0 からの扇で三角形化）。
///
/// コーパスの面はすべて四角形以下なので、**5 辺以上に併合される配置**と
/// **非凸で併合できない配置**をここで作ります。
TriMesh prism(const std::vector<std::int32_t>& px, const std::vector<std::int32_t>& py,
              std::int32_t zlo, std::int32_t zhi) {
    const std::size_t n = px.size();
    TriMesh m;
    for (std::size_t i = 0; i < n; ++i) m.vertices.push_back({px[i], py[i], zlo});
    for (std::size_t i = 0; i < n; ++i) m.vertices.push_back({px[i], py[i], zhi});
    const auto u = [n](std::size_t i) { return static_cast<std::uint32_t>(i); };
    const auto d = [n](std::size_t i) { return static_cast<std::uint32_t>(i + n); };
    for (std::size_t i = 1; i + 1 < n; ++i) {
        m.triangles.push_back({u(0), u(i + 1), u(i)});  // 底面（法線 -z）
        m.triangles.push_back({d(0), d(i), d(i + 1)});  // 上面（法線 +z）
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        m.triangles.push_back({u(i), u(j), d(j)});
        m.triangles.push_back({u(i), d(j), d(i)});
    }
    return m;
}

/// 正六角形（凸）。**面が 6 辺に併合されるはず。**
TriMesh hex_prism() {
    using kritest::at;
    const std::vector<std::int32_t> px = {at(120, 256),  at(60, 256),  at(-60, 256),
                                          at(-120, 256), at(-60, 256), at(60, 256)};
    const std::vector<std::int32_t> py = {at(7, 256), at(111, 256), at(111, 256),
                                          at(7, 256), at(-97, 256), at(-97, 256)};
    return prism(px, py, at(-128, 256), at(128, 256));
}

/// L 字（**非凸**）。**面は 1 つに併合できない。**
///
/// 併合の凸性判定が壊れていれば、半空間の交わりは**凸包**になるので、
/// 点集合が変わって差が空になりません。
TriMesh l_prism() {
    using kritest::at;
    const std::vector<std::int32_t> px = {at(-121, 256), at(135, 256), at(135, 256),
                                          at(7, 256),    at(7, 256),   at(-121, 256)};
    const std::vector<std::int32_t> py = {at(-121, 256), at(-121, 256), at(7, 256),
                                          at(7, 256),    at(135, 256),  at(135, 256)};
    return prism(px, py, at(-128, 256), at(128, 256));
}

/// **点集合として一致するか**（両方向の差が空。`test_nary.cpp` と同じ形）。
///
/// $(C,\chi)$ は形が変わっても保たれることがあります。**凸性の判定が壊れて
/// 凸包になった**ような誤りは、ここでしか出ません。
void check_same_solid(const TriMesh& m, const std::string& tag) {
    BoolOptions o = all_on(1, false);
    ToMeshOptions tm;
    tm.split_contacts = true;
    const PolySoup a = from_mesh(m, tri_opt()), b = from_mesh(m, cvx_opt());
    KRI_CHECK_MSG(!to_mesh(a, tm).triangles.empty(), tag + ": 三角形化の側が空");
    for (int dir = 0; dir < 2; ++dir) {
        const PolySoup& x = (dir == 0) ? a : b;
        const PolySoup& y = (dir == 0) ? b : a;
        const SoupMesh dm = to_mesh(boolean(x, y, BoolOp::Difference, o), tm);
        KRI_CHECK_MSG(
            dm.triangles.empty(),
            tag +
                (dir == 0 ? ": 三角形化 \\ 凸分割 が空でない" : ": 凸分割 \\ 三角形化 が空でない") +
                "（三角形 " + std::to_string(dm.triangles.size()) + " 枚）");
    }
}

/// 0. 局所メッシュ — 5 辺以上への併合と、非凸で併合しないこと。
void test_shapes() {
    struct Item {
        const char* name;
        TriMesh (*make)();
        std::size_t want_max_edges;
    };
    const Item items[] = {{"六角柱（凸）", hex_prism, 6}, {"L 字柱（非凸）", l_prism, 4}};
    for (const Item& it : items) {
        const TriMesh m = it.make();
        const PolySoup b = from_mesh(m, cvx_opt());
        std::size_t max_edges = 0, on_top = 0;
        PlaneId top = kNoPlane;
        for (const Poly& q : b.polys) {
            max_edges = std::max(max_edges, q.frag.edge.size());
        }
        // 上面（$+z$）の片が何枚残ったか
        for (const Poly& q : b.polys) {
            if (q.frag.edge.size() == max_edges) {
                top = q.frag.support;
                break;
            }
        }
        for (const Poly& q : b.polys) {
            if (q.frag.support == top) ++on_top;
        }
        KRI_CHECK_MSG(max_edges == it.want_max_edges,
                      std::string(it.name) + ": 最大辺数が想定と違う" +
                          kritest::pair_msg(it.want_max_edges, max_edges));
        check_same_solid(m, it.name);
        std::printf("    %s: 三角形 %zu → 片 %zu（最大辺数 %zu、同一平面 %zu 枚）\n", it.name,
                    m.triangles.size(), b.polys.size(), max_edges, on_top);
    }
}

/// 1. 入口だけの比較 — 片数・平面数が減り、往復が壊れないこと。
void test_entry() {
    std::size_t t = 0, p_tri = 0, p_cvx = 0, pl_tri = 0, pl_cvx = 0, n = 0;
    std::size_t max_edges = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        for (int w = 0; w < 2; ++w) {
            const TriMesh m = (w == 0) ? c.make_a() : c.make_b();
            const PolySoup a = from_mesh(m, tri_opt());
            const PolySoup b = from_mesh(m, cvx_opt());
            const std::string tag = std::string("ケース ") + c.id + (w == 0 ? " A" : " B");

            // **決定性**（§4.1.2）。位相も体積も見ないので、ここでしか捕まりません
            KRI_CHECK_MSG(soup_key(b) == soup_key(from_mesh(m, cvx_opt())),
                          tag + ": **凸分割が決定的でない**（同じ入力で結果が変わった）");

            // 往復が壊れないこと（入口と出口だけで閉じる検査）
            const TopologyReport r0 = check_topology(m.triangles);
            const TopologyReport r1 = check_topology(to_mesh(b).triangles);
            KRI_CHECK_MSG(r1.ok(), tag + ": 凸分割の往復で多様体でなくなった");
            KRI_CHECK_MSG(r0.components == r1.components,
                          tag + ": 凸分割の往復で C が変わった" +
                              kritest::pair_msg(r0.components, r1.components));
            KRI_CHECK_MSG(r0.chi == r1.chi, tag + ": 凸分割の往復で χ が変わった" +
                                                kritest::pair_msg(r0.chi, r1.chi));

            KRI_CHECK_MSG(b.polys.size() <= a.polys.size(), tag + ": 凸分割で片が増えた");
            check_same_solid(m, tag);
            for (const Poly& q : b.polys) max_edges = std::max(max_edges, q.frag.edge.size());
            t += m.triangles.size();
            p_tri += a.polys.size();
            p_cvx += b.polys.size();
            pl_tri += a.table.size();
            pl_cvx += b.table.size();
            ++n;
        }
    }
    // **一致しないことも検査する。** 一致したら凸分割が空回りしています
    KRI_CHECK_MSG(p_cvx < p_tri, "**片数が減っていない**（凸分割が効いていない）");
    KRI_CHECK_MSG(pl_cvx < pl_tri, "**平面数が減っていない**（§4.1 の表）");
    KRI_CHECK_MSG(max_edges > 3, "**4 辺以上の片が 1 つも無い**（併合が起きていない）");
    std::printf("    入口 %zu 件: 片 %zu → %zu（%.1f%%）, 平面 %zu → %zu（%.1f%%）, 最大辺数 %zu\n",
                n, p_tri, p_cvx, 100.0 * static_cast<double>(p_cvx) / static_cast<double>(p_tri),
                pl_tri, pl_cvx, 100.0 * static_cast<double>(pl_cvx) / static_cast<double>(pl_tri),
                max_edges);
    (void)t;
}

/// 2. ブール演算を通した比較 — $(C, \chi)$ は一致し、三角形の集合は一致しない。
void test_boolean() {
    std::size_t tri_frags = 0, cvx_frags = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            for (std::size_t pass = 0; pass < kPasses; ++pass) {
                const bool adaptive = (pass == kPasses - 1);
                const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(pass);
                const BoolOptions o = all_on(depth, adaptive);
                ToMeshOptions tm;
                tm.split_contacts = true;
                BoolStats st_t{}, st_c{};
                const SoupMesh mt = to_mesh(
                    boolean(from_mesh(a, tri_opt()), from_mesh(b, tri_opt()), op, o, &st_t), tm);
                const SoupMesh mc = to_mesh(
                    boolean(from_mesh(a, cvx_opt()), from_mesh(b, cvx_opt()), op, o, &st_c), tm);
                const TopologyReport rt = check_topology(mt.triangles);
                const TopologyReport rc = check_topology(mc.triangles);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（" +
                                        (adaptive ? "適応" : "深度 " + std::to_string(depth)) +
                                        "）";
                KRI_CHECK_MSG(rt.components == rc.components,
                              tag + ": 凸分割で C が変わった" +
                                  kritest::pair_msg(rt.components, rc.components));
                KRI_CHECK_MSG(rt.chi == rc.chi,
                              tag + ": 凸分割で χ が変わった" + kritest::pair_msg(rt.chi, rc.chi));
                KRI_CHECK_MSG(rc.empty || rc.ok(), tag + ": 凸分割の出力が多様体でない");
                if (mt.triangles.size() != mc.triangles.size()) ++g_differ_tris;
                tri_frags += st_t.fragments;
                cvx_frags += st_c.fragments;
                ++g_cmp;
            }
        }
    }
    // **一致しないことも検査する**（§4.1.1）
    KRI_CHECK_MSG(g_differ_tris > 0,
                  "**三角形の集合が 1 件も変わらない**（凸分割が下流に届いていない）");
    KRI_CHECK_MSG(cvx_frags < tri_frags, "**断片数が減っていない**（§4.1.1）");
    const std::size_t want = kritest::corpus().size() * 3 * kPasses;
    KRI_CHECK_MSG(g_cmp == want, "比較数が式と合わない" + kritest::pair_msg(want, g_cmp));
    std::printf("    ブール %zu 件: 断片 %zu → %zu（%.1f%%）, 三角形が変わった %zu 件\n", g_cmp,
                tri_frags, cvx_frags,
                100.0 * static_cast<double>(cvx_frags) / static_cast<double>(tri_frags),
                g_differ_tris);
}

}  // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    std::printf("\n  凸分割（SPEC-phase3 §4.1）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    test_shapes();
    test_entry();
    test_boolean();
    std::printf("\n");
    return kritest::finish("csg/convex");
}
