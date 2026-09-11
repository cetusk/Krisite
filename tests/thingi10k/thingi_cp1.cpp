// **★ この実行ファイルは `build/tests/cp1`（コーパスのテスト、`tests/csg/test_cp1.cpp`）
// とは【別物】です。** 以前は両方とも `cp1` という名前になり得ました
// （`CLAUDE.md`「同じ名前の実行ファイルを 2 つ作らないでください」）。
//
//     tests/csg/test_cp1.cpp      コーパスのテスト。`ctest` に入っている
//     tests/thingi10k/thingi_cp1.cpp  **実データのドライバ。`ctest` に入っていない**
//
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
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/par/thread_pool.hpp"

#include "corpus_expect.hpp"
#include "thingi10k/loader.hpp"
#include "volume_fp.hpp"

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

/// **対ごとの構造**（`SPEC-phase5.md` §1.5.0）。
///
/// **測定は、要求されなければ実装されません。** CP1 は対ごとの構造を 1 つも
/// 記録しておらず、**失敗と時間の切り分けが事後にできませんでした。**
///
/// **3 演算ぶんを合算します**（最大の項目は最大を取る）。1 対で 1 行にするためです。
struct PairStruct {
    std::size_t polys = 0, fragments = 0, regions = 0;
    std::size_t raycasts = 0, ray_tri_tests = 0;
    std::size_t leaf_nonempty = 0, leaf_input_total = 0, leaf_input_max = 0;
    std::size_t max_planes_per_cell = 0;
    std::size_t leaf_single_src = 0, leaf_single_max = 0, leaf_both_max = 0;
    std::size_t bsp_slots = 0, bsp_used = 0, bsp_slots_single = 0, bsp_used_single = 0;
    std::size_t bsp_cells_skipped = 0;
    /// **浮動小数点の体積恒等式の相対誤差**（`SPEC-phase5.md` §3.0 の篩）。
    /// **検査ではありません。** GMP に回すものを絞るために記録します。
    double vol_err = 0, diff_err = 0;
    double ms_arrange = 0, ms_classify = 0, ms_stitch = 0;
    std::size_t unresolved = 0, edges_excess = 0, edges_deficient = 0;
    std::size_t max_edge_degree = 0;
    // ---- 除外の判定に要る 4 項目（`SPEC-phase2.md` §9.3.1）----------------------
    //
    // **除外の 7 条件のうち、記録していたのは 3 つだけでした**
    // （`unresolved > 0` / 不足辺 0 / 最大次数 ≤ 4）。**残る 4 つが無いので、
    // 実データの記録からは除外かどうかを判定できませんでした**（`IMPL-phase5.md` §48）。
    //
    // **`check_topology` は既に呼ばれており、`TopologyReport` に全項目が入っています。
    // 記録していないだけでした。**
    std::size_t edges_odd_degree = 0;
    std::size_t nonmanifold_unexplained = 0;
    /// 3 演算すべてで成立したか（**論理積**。1 つでも崩れたら除外できません）。
    ///
    /// **★ 空の出力を論理積に入れないこと**（2026-09-06 に踏みました）。
    /// `check_topology` は**空のメッシュで早期に返し、`oriented` も `no_degenerate` も
    /// 偽のまま**です。そのまま論理積を取ると、
    /// **「向きが崩れた」と「出力が空だった」が区別できません。**
    /// 実際 500 対のうち 116 対（23.2%）が 0 でしたが、**その全部が
    /// 「∩ が空」などの正常な事象**でした（離れた 2 立体の共通部分など）。
    /// **偽が 2 つの意味を持つ真偽値を作らないこと。**
    int all_oriented = 1, all_no_degenerate = 1;
    /// **空の出力を出した演算の数**（0〜3）。上の論理積から外した分をここで数えます
    int empty_ops = 0;
    // ---- 三角形化が残した退化三角形（`SPEC-phase2.md` §2.4.4 (2)）----
    //
    // **★ Phase 2 が「残した枚数を §11 に記録すること」と要求していた項目です。**
    // **Phase 5 の記録項目に引き継いでいませんでした**（`IMPL-phase5.md` §95）。
    //
    // コーパスでは 20 枚でしたが、**実データでは 1 演算あたり最大 59,465 枚
    // （出力の 0.95%）**です。3 桁違います。
    //
    // **判定は組合せです**（幾何の述語を使うと $20b+43$ が要る。§2.4.4 (2) の禁止）。
    std::size_t degenerate_kept = 0;
    std::size_t apex_fallback = 0;  ///< 扇の起点を選べなかった多角形
    /// **radial sort**（`SPEC-phase2.md` §5.1.2.1）の試行と解決
    std::size_t radial_attempted = 0;
    std::size_t radial_resolved = 0;
    /// **★ 接触の分裂が発火した回数**（`SPEC-phase5.md` §1.5.0。2026-09-06 の要求）。
    ///
    /// **仕様が要求していたのに、記録項目に入っていませんでした**（`IMPL-v2.md` §3）。
    /// 同じ改訂で足された `degenerate_kept` と radial は入っており、**これだけ落ちました。**
    ///
    /// > **`unresolved` は「解けなかった数」であって「解いた数」ではありません**
    /// > （§1.5.0.1）。**分裂が空回りしていた可能性を、これがないと否定できません。**
    std::size_t split_vertices = 0;
    /// **`unresolved` のうち、事後の検査（分裂の【後】に非多様体）で数えたもの**。
    ///
    /// もう一方は「対応付けできなかった辺」です（`unsplit_edges`）。
    /// **`IMPL-phase5.md` §98 の残件はこちらです。分けて数えないと区別できません。**
    std::size_t unresolved_post = 0;
    std::size_t unsplit_edges = 0;  ///< 対応付けできず、分裂させずに残した辺
    // ---- ★ 出口の 5 段（`SPEC-phase5.md` §5.11.1 / §6）------------------------
    //
    // **§6 が「出口の内訳」を要求しているのに、記録項目に入っていませんでした。**
    // **`ms_arrange` / `ms_classify` / `ms_stitch`（中核）だけがありました。**
    double ms_construct = 0;  ///< 構成点
    double ms_merge = 0;      ///< 値で併合
    double ms_index = 0;      ///< 平面索引（候補の整列を含む）
    double ms_tri = 0;        ///< T 解決 + 三角形化
    double ms_split = 0;      ///< 接触の分裂
    double ms_tomesh = 0;     ///< 出口の全体（5 段の和と一致するはず）
    /// **中核（`boolean`）の全体**。`ms_arrange` などの和より大きい（前処理・葉の列挙を含む）
    double ms_core = 0;
    /// **入口（`from_mesh`）**。`fm_seconds` は検査込みなので、別に採る
    double ms_inlet = 0;
    /// **★ 駆動側の検査・検算**（`SPEC-phase5.md`
    /// §5.10.14。実務では切る部分。分母から外すために別に採る）
    double ms_topo = 0;  ///< `check_topology`（演算ごとに 2 回: 統計用と合否用）
    double ms_vol = 0;   ///< 体積の検算（`volume6_fp` × 5）
    double ms_hash = 0;  ///< 出力のハッシュ（× 3）
    /// **分類と arrange の内訳**（§5.10.14.7。3 演算の和）
    double cl_prep = 0, cl_rep = 0, cl_ray = 0, cl_out = 0, cl_par_wall = 0;
    double ar_gather = 0, ar_present = 0, ar_prep = 0, ar_frag = 0, ar_coplanar = 0, ar_stitch = 0;
    /// **索引の粒度**（§5.10.14.11。段ごとの候補と項目。3 演算の和）
    std::size_t cand_level[12] = {}, items_level[12] = {};
    std::size_t tri_level[12] = {}, fit1_level[12] = {}, cells0_level[12] = {};
    std::size_t hist_tri[6] = {}, hist_cells[6] = {};
    std::size_t split_attempts = 0, split_actual = 0, uncut = 0, uncut_tris = 0;
    /// **版をまたいだ意味論の比較のため**: 3 演算の出力の 6 倍体積（浮動小数の篩）と、三角形数
    double vol3[3] = {0, 0, 0};
    std::size_t skip_box = 0, skip_exact = 0, skip_boxside = 0;
    double fr_clip = 0, fr_prep = 0, fr_cut = 0, fr_commit = 0;
    std::size_t tri3[3] = {0, 0, 0};
    /// **除外できた演算の数**（0〜3）。`3` なら 3 演算すべてが除外の条件を満たす
    int excluded_ops = 0;
    /// **NSI を宣言できたか**（-1 = 検査していない / 0 = 自己交差あり / 1 = 宣言した）。
    /// **検査は量子化後に走る**ので、模型ごとにキャッシュできません（§33）。
    int nsi_a = -1, nsi_b = -1;
    // ---- Phase 5 の 3 機構の発火（2026-09-05 追加）--------------------------------
    //
    // **実データの分布で、どれだけ効いたかを記録します。**
    // **記録していなければ、後から「ログにない」となります**（`CLAUDE.md`）。
    std::size_t cell_index_groups = 0;   ///< A-3: (葉, 支持平面) の組。**0 なら退避**
    std::size_t assign_rejected = 0;     ///< 10-a: 落とした割り当て
    std::size_t ray_kept = 0;            ///< D: 前判定を通った候補
    std::size_t regions_negative_w = 0;  ///< **負の巻き数を持つ領域**（§2.9 の訂正が効いたか）
    std::size_t regions_w_ge2 = 0;       ///< 巻き数 2 以上の領域
    /// **`from_mesh` 2 回ぶんの秒数（検査込み）。**
    ///
    /// **検査だけの時間ではありません。** 分けて測ろうとすると `from_mesh` を
    /// 2 度走らせることになるので、**含んだまま記録して、そう書きます。**
    /// 検査単独の費用は `quantstat` が別に出しています（実測 0.09 秒/模型）。
    double fm_seconds = 0;

    void add(const csg::BoolStats& b, const csg::ToMeshStats& t, const mesh::TopologyReport& r,
             std::size_t np) {
        polys += np;
        fragments += b.fragments;
        regions += b.regions;
        raycasts += b.raycasts;
        ray_tri_tests += b.ray_tri_tests;
        leaf_nonempty += b.leaf_nonempty;
        leaf_input_total += b.leaf_input_total;
        leaf_input_max = std::max(leaf_input_max, b.leaf_input_max);
        max_planes_per_cell = std::max(max_planes_per_cell, b.max_planes_per_cell);
        leaf_single_src += b.leaf_single_src;
        leaf_single_max = std::max(leaf_single_max, b.leaf_single_src_input_max);
        leaf_both_max = std::max(leaf_both_max, b.leaf_both_input_max);
        bsp_slots += b.bsp_cut_slots;
        bsp_used += b.bsp_cuts_used;
        bsp_slots_single += b.bsp_cut_slots_single;
        bsp_used_single += b.bsp_cuts_used_single;
        bsp_cells_skipped += b.bsp_cells_skipped_nsi;
        // Phase 5 の 3 機構（2026-09-05 追加）
        cell_index_groups += t.cell_index_groups;
        assign_rejected += b.assign_rejected_plane;
        ray_kept += b.ray_tri_kept;
        regions_negative_w += b.regions_negative_w;
        regions_w_ge2 += b.regions_w_ge2;
        degenerate_kept += t.t.degenerate_kept;
        apex_fallback += t.t.apex_fallback;
        radial_attempted += t.split.radial_attempted;
        radial_resolved += t.split.radial_resolved;
        split_vertices += t.split.split_vertices;
        unresolved_post += t.split.unresolved_post;
        unsplit_edges += t.split.unsplit_edges;
        ms_construct += t.ms_construct;
        ms_merge += t.ms_merge;
        ms_index += t.ms_index;
        ms_tri += t.ms_tri;
        ms_split += t.ms_split;
        ms_tomesh += t.ms_total;
        cl_prep += b.ms_cl_prep;
        cl_rep += b.ms_cl_rep;
        cl_ray += b.ms_cl_ray;
        cl_out += b.ms_cl_out;
        cl_par_wall += b.ms_cl_par_wall;
        ar_gather += b.ms_arr_gather;
        ar_present += b.ms_arr_present;
        ar_prep += b.ms_arr_prep;
        ar_frag += b.ms_arr_frag;
        ar_coplanar += b.ms_arr_coplanar;
        ar_stitch += b.ms_arr_stitch;
        for (int l = 0; l < 12; ++l) {
            cand_level[l] += b.ray_cand_level[l];
            items_level[l] += b.ray_items_level[l];
            tri_level[l] += b.ray_tri_level[l];
            fit1_level[l] += b.ray_fit1_level[l];
            cells0_level[l] += b.ray_cells0_level[l];
        }
        for (int k = 0; k < 6; ++k) {
            hist_tri[k] += b.ray_hist_tri[k];
            hist_cells[k] += b.ray_hist_cells[k];
        }
        split_attempts += b.bsp_split_attempts;
        skip_box += b.bsp_skip_box;
        fr_clip += b.ms_fr_clip;
        fr_prep += b.ms_fr_prep;
        fr_cut += b.ms_fr_cut;
        fr_commit += b.ms_fr_commit;
        skip_exact += b.bsp_skip_exact;
        skip_boxside += b.bsp_skip_boxside;
        split_actual += b.bsp_split_actual;
        uncut += b.frags_uncut;
        uncut_tris += b.frags_uncut_tris;
        ms_arrange += b.ms_arrange;
        ms_classify += b.ms_classify;
        ms_stitch += b.ms_stitch;
        unresolved += t.split.unresolved;
        edges_excess += r.edges_excess;
        edges_deficient += r.edges_deficient;
        max_edge_degree = std::max(max_edge_degree, r.max_edge_degree);
        // **除外の判定に要る 4 項目**（§9.3.1）
        edges_odd_degree += r.edges_odd_degree;
        nonmanifold_unexplained += r.nonmanifold_vertices_unexplained;
        if (r.empty) {
            ++empty_ops;  // **空は論理積に入れない**（上の注記）
        } else {
            all_oriented &= r.oriented ? 1 : 0;
            all_no_degenerate &= r.no_degenerate ? 1 : 0;
        }
        {
            const kritest::Exclusion ex = kritest::exclusion_when_split(t.split.unresolved, r);
            std::string why;
            if (ex != kritest::Exclusion::None && kritest::exclusion_conditions_ok(ex, r, &why)) {
                ++excluded_ops;
            }
        }
    }
    void print(std::ostream& o) const {
        o << polys << ' ' << fragments << ' ' << regions << ' ' << raycasts << ' ' << ray_tri_tests
          << ' ' << leaf_nonempty << ' ' << leaf_input_total << ' ' << leaf_input_max << ' '
          << max_planes_per_cell << ' ' << (long long)ms_arrange << ' ' << (long long)ms_classify
          << ' ' << (long long)ms_stitch << ' ' << unresolved << ' ' << edges_excess << ' '
          << edges_deficient << ' ' << max_edge_degree << ' ' << leaf_single_src << ' '
          << leaf_single_max << ' ' << leaf_both_max << ' ' << bsp_slots << ' ' << bsp_used << ' '
          << bsp_slots_single << ' ' << bsp_used_single << ' ' << vol_err << ' ' << diff_err << ' '
          << nsi_a << ' ' << nsi_b << ' ' << fm_seconds << ' ' << edges_odd_degree << ' '
          << nonmanifold_unexplained << ' ' << all_oriented << ' ' << all_no_degenerate << ' '
          << excluded_ops
          // Phase 5 の 3 機構（2026-09-05 追加）
          << ' ' << cell_index_groups << ' ' << assign_rejected << ' ' << ray_kept << ' '
          << regions_negative_w << ' '
          << regions_w_ge2
          // **空の出力の数**（2026-09-06 追加。列を足したので `cp1_results.txt` の
          // 既存 500 行にはありません。**CP2 以降の記録に入ります**）
          << ' '
          << empty_ops
          // **仕様が要求していた記録**（`SPEC-phase2.md` §2.4.4 (2) / §5.1.2.1）
          << ' ' << degenerate_kept << ' ' << apex_fallback << ' ' << radial_attempted << ' '
          << radial_resolved
          // **§1.5.0 が要求していたのに入っていなかった項目**（2026-09-07 追加。`IMPL-v2.md` §3）
          << ' ' << split_vertices << ' ' << unresolved_post << ' '
          << unsplit_edges
          // **出口の 5 段と、入口・中核・出口の合計**（2026-09-08 追加。§5.11.1）
          << ' ' << (long long)ms_construct << ' ' << (long long)ms_merge << ' '
          << (long long)ms_index << ' ' << (long long)ms_tri << ' ' << (long long)ms_split << ' '
          << (long long)ms_tomesh << ' ' << (long long)ms_core << ' ' << (long long)ms_inlet << ' '
          << (long long)ms_topo << ' ' << (long long)ms_vol << ' ' << (long long)ms_hash << ' '
          << (long long)cl_prep << ' ' << (long long)cl_rep << ' ' << (long long)cl_ray << ' '
          << (long long)cl_out << ' ' << (long long)cl_par_wall << ' ' << (long long)ar_gather
          << ' ' << (long long)ar_present << ' ' << (long long)ar_prep << ' ' << (long long)ar_frag
          << ' ' << (long long)ar_coplanar << ' ' << (long long)ar_stitch;
        for (int l = 0; l < 12; ++l) o << ' ' << cand_level[l];
        for (int l = 0; l < 12; ++l) o << ' ' << items_level[l];
        for (int l = 0; l < 12; ++l) o << ' ' << tri_level[l];
        for (int l = 0; l < 12; ++l) o << ' ' << fit1_level[l];
        for (int l = 0; l < 12; ++l) o << ' ' << cells0_level[l];
        for (int k = 0; k < 6; ++k) o << ' ' << hist_tri[k];
        for (int k = 0; k < 6; ++k) o << ' ' << hist_cells[k];
        o << ' ' << split_attempts << ' ' << split_actual << ' ' << uncut << ' ' << uncut_tris;
        for (int k = 0; k < 3; ++k) o << ' ' << std::setprecision(17) << vol3[k];
        for (int k = 0; k < 3; ++k) o << ' ' << tri3[k];
        o << ' ' << skip_box << ' ' << skip_exact << ' ' << (long long)fr_clip << ' '
          << (long long)fr_prep << ' ' << (long long)fr_cut << ' ' << (long long)fr_commit << ' '
          << skip_boxside;
    }
};

/// §3.1 の検査。**解析的期待値は使えない**ので恒等式と位相で見ます。
bool check_one(const mesh::TriMesh& a, const mesh::TriMesh& b, const csg::BoolOptions& o,
               par::ThreadPool* pool, std::string* why, unsigned long long* hash_out = nullptr,
               PairStruct* ps = nullptr, int nsi_decl = 0, bool verify_delta = true) {
    // **NSI は呼び出し側が宣言します**（`SPEC-phase3.md` §5.6、EMBER §4.5.1）。
    // ライブラリは検証しません。**宣言してよいかを確かめるのは呼び出し側の仕事**で、
    // `from_mesh` の `verify_nsi` がその補助です。
    //
    // **CP1 の母集団はメタデータで「自己交差なし」と分類されたモデルですが、
    // 量子化後の性質は保証されません**（§1.0。実測で 1,000 件中 59〜70 件が該当）。
    // だからこそ**既定は宣言しない**で、`nsi_decl == 1` のときだけ検査して通ったものを宣言します。
    //
    //   0 … 宣言しない（従来）
    //   1 … **検査して、通ったものだけ宣言する**（既定の運用）
    //   2 … 検査せずに宣言する（旧挙動。§20 の測定を再現するためだけに残す）
    csg::FromMeshOptions fm;
    fm.verify_nsi = (nsi_decl == 1);
    const auto tsi = std::chrono::steady_clock::now();
    csg::PolySoup A = csg::from_mesh(a, fm), B = csg::from_mesh(b, fm);
    if (ps != nullptr) {
        ps->ms_inlet =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tsi)
                .count();
    }
    if (ps != nullptr && nsi_decl == 1) {
        ps->fm_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tsi).count();
        ps->nsi_a = A.nsi[0] ? 1 : 0;
        ps->nsi_b = B.nsi[0] ? 1 : 0;
    }
    if (nsi_decl == 2) {
        A.nsi.assign(A.sources.size(), 1);
        B.nsi.assign(B.sources.size(), 1);
    }
    csg::ToMeshOptions tm;
    tm.split_contacts = true;
    tm.threads = pool->size();
    tm.pool = pool;
    // **§5.5 の検算は CP1〜CP3 では ON**（`SPEC-phase5.md` §3.2）。
    //
    // **既定は偽です**（純粋な診断を本番の経路に置かないため。`HANDOVER.md` §5.1）。
    // **実データの駆動プログラムは明示的に立てます。** 中身は同値な増分計算なので、
    // 立てても費用はほとんど増えません（`IMPL-v2.md` §2）。
    //
    // > **★ ただし §4.3.2（EMBER と比較可能な数字）は「検算 OFF」を要求します。**
    // > **測る目的が違うので、引数で切り替えます。**
    // > **CP1〜CP3 では既定（真）のまま回してください。**
    tm.verify_split_delta = verify_delta;
    // **対ごとの構造を採ります**（`SPEC-phase5.md` §1.5.0）。3 演算ぶんを合算。
    csg::SoupMesh out3[3];
    int k3 = 0;
    // **★ 段ごとの進捗を出します**（2026-09-05 追加）。
    //
    // **対の中で「いまどの段にいるか」を出していませんでした。**
    // 1 対に 2 時間以上かかったとき、**「進んでいるか」を判断する材料がありませんでした**
    // （プロセスの CPU 時間しか見えず、`gdb` もこの環境にありません）。
    //
    // **`CLAUDE.md`「測定は、要求されなければ実装されません」**に当たります。
    //
    // **費用はほぼゼロです**（対あたり 6 行）。**打ち切りの判断に直接効きます。**
    const char* kOpName[3] = {"∪", "∩", "＼"};
    const auto t_pair = std::chrono::steady_clock::now();
    const auto lap = [&t_pair]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t_pair).count();
    };
    for (csg::BoolOp op :
         {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
        csg::BoolStats bs;
        csg::ToMeshStats ts;
        std::printf("      [%6.1f s] %s 中核…\n", lap(), kOpName[k3]);
        std::fflush(stdout);
        const auto t_core = std::chrono::steady_clock::now();
        const csg::PolySoup soup = csg::boolean(A, B, op, o, &bs);
        if (ps != nullptr) {
            ps->ms_core +=
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_core)
                    .count();
        }
        std::printf("      [%6.1f s] %s 中核 完了（多角形 %zu、断片 %zu）→ 出口…\n", lap(),
                    kOpName[k3], soup.polys.size(), bs.raw_fragments);
        std::fflush(stdout);
        out3[k3] = csg::to_mesh(soup, tm, &ts);
        std::printf("      [%6.1f s] %s 出口 完了（三角形 %zu、頂点 %zu）\n", lap(), kOpName[k3],
                    out3[k3].triangles.size(), out3[k3].vertices.size());
        std::fflush(stdout);
        if (ps != nullptr) {
            const auto t_tp = std::chrono::steady_clock::now();
            const mesh::TopologyReport tr = mesh::check_topology(out3[k3].triangles);
            ps->ms_topo +=
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_tp)
                    .count();
            ps->add(bs, ts, tr, soup.polys.size());
        }
        ++k3;
    }
    const csg::SoupMesh& mu = out3[0];
    const csg::SoupMesh& mi = out3[1];
    const csg::SoupMesh& md = out3[2];

    const auto t_vol = std::chrono::steady_clock::now();
    if (ps != nullptr) {
        // **浮動小数点の体積恒等式**（`SPEC-phase5.md` §3.0）。**篩であって検査ではありません。**
        // 入力側は `signed_volume6` が厳密なので、丸めは出力側だけに乗ります。
        const double va = kritest::volume6_fp(a), vb = kritest::volume6_fp(b);
        const double vu = kritest::volume6_fp(mu), vi = kritest::volume6_fp(mi);
        const double vd = kritest::volume6_fp(md);
        ps->vol_err = kritest::identity_error(vu, vi, va, vb);
        ps->vol3[0] = vu;
        ps->vol3[1] = vi;
        ps->vol3[2] = vd;
        ps->tri3[0] = mu.triangles.size();
        ps->tri3[1] = mi.triangles.size();
        ps->tri3[2] = md.triangles.size();
        ps->diff_err = kritest::difference_error(vd, va, vi);
        ps->ms_vol =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_vol)
                .count();
    }
    const auto t_tp2 = std::chrono::steady_clock::now();

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
    if (ps != nullptr) {
        ps->ms_topo +=
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_tp2)
                .count();
    }
    if (ps != nullptr && ps->vol_err > kritest::kIdentityTol) {
        // **篩に引っかかりました。** 失敗ではありません — **GMP に回す印**です（§3.0）。
        *why = "体積の篩（GMP へ）";
    }
    if (hash_out != nullptr) {
        const auto t_h = std::chrono::steady_clock::now();
        *hash_out = hash_mesh(mu) ^ (hash_mesh(mi) * 3) ^ (hash_mesh(md) * 7);
        if (ps != nullptr) {
            ps->ms_hash =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_h)
                    .count();
        }
    }
    // **体積の恒等式**（§3.1）。|A∪B| + |A∩B| = |A| + |B| を出力側の 6 倍体積で見る
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = (argc > 1) ? argv[1] : "data/thingi10k/cp1.txt";
    // **付随ファイルは一覧の名前から機械的に決めます**（`cp1.txt` → `cp1_*`）。
    //
    // **CP2 / CP3 で書き換えないため**です。手で書くと、CP1 の記録に CP2 の結果を
    // 追記する事故が起きます（§40 の変換の取り違えと同じ形 — **手で書ける口を塞ぐ**）。
    const std::string stem = [&list] {
        const std::size_t sl = list.find_last_of('/');
        std::string b = (sl == std::string::npos) ? list : list.substr(sl + 1);
        const std::size_t dot = b.rfind(".txt");
        return (dot == std::string::npos) ? b : b.substr(0, dot);
    }();
    const std::string base = "data/thingi10k/" + stem;
    const std::size_t limit = (argc > 2) ? std::strtoul(argv[2], nullptr, 10) : 0;
    const unsigned depth = (argc > 3) ? static_cast<unsigned>(std::atoi(argv[3])) : 6;
    const unsigned nthreads = (argc > 4) ? static_cast<unsigned>(std::atoi(argv[4])) : 16;
    // **索引の ON/OFF でハッシュが一致するかを対ごとに確かめるモード**（CP1.5D）。
    const bool verify_index = (argc > 5) && (std::atoi(argv[5]) != 0);
    // **済みの対をやり直すモード。** 既定は追記（再開）
    const bool redo = (argc > 6) && (std::atoi(argv[6]) != 0);
    // **NSI の扱い**:
    //   0 … 宣言しない（従来）
    //   1 … **検査して、通ったものだけ宣言する**（(a) の既定の運用）
    //   2 … 検査つきの宣言と、宣言なしの両方を回して突き合わせる
    //   3 … 検査せずに宣言する（旧 1。§20 の測定を再現するためだけに残す）
    const int nsi_mode = (argc > 7) ? std::atoi(argv[7]) : 0;
    // **§5.5 の検算を切る**（`SPEC-phase5.md` §4.3.2 の「検算 OFF」）。
    //
    // > **既定は真です。CP1〜CP3 は §3.2 が ON を要求しています。**
    // > **切るのは EMBER と比較可能な数字を採るときだけ**で、
    // > **そのときは正しさの判定に使わないでください。**
    const bool verify_delta = (argc > 8) ? (std::atoi(argv[8]) != 0) : true;
    // **レイ索引の細かい割り当ての上限 K**（第 9 引数。既定 0 = 従来。§5.10.14.15）
    const std::size_t fine_k = (argc > 9) ? std::strtoul(argv[9], nullptr, 10) : 0;
    // 宣言の内訳（(a) の空回り検査）。**「検査を入れた」と「検査が効いた」は別**なので、
    // **宣言できた数と落ちた数を必ず出します。**
    std::size_t nsi_declared = 0, nsi_rejected = 0;
    std::size_t verified = 0;
    // **この対だけを回す**（空なら全件）。失敗の分類のように、
    // **少数の対だけ構造を採りたい**場面のためです（`SPEC-phase5.md` §1.5.4）。
    std::vector<std::string> only;
    {
        std::ifstream f(base + "_only.txt");
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty()) only.push_back(line.substr(0, line.find(' ')));
        }
    }
    // 資源上限で落ちた対（1 行 1 対）。**手で足すのではなく、監視スクリプトが足します**
    std::vector<std::string> skip;
    {
        std::ifstream f(base + "_skip.txt");
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
    o.measure_classify = true;   // 分類の内訳（§5.10.14.7。領域ごとに時計 4 回）
    o.measure_frag = true;       // 断片の生成の内訳（§5.10.14.30。多角形ごとに時計 4 回）
    o.record_ray_levels = true;  // 索引の粒度（§5.10.14.11）
    o.ray_index_fine_cells = fine_k;
    // **第 10 引数: 記憶の上限（項目 / 三角形）から K を導く形**
    const std::size_t fine_budget = (argc > 10) ? std::strtoul(argv[10], nullptr, 10) : 0;
    o.ray_index_fine_budget = fine_budget;
    // **第 11 引数: 記憶の上限（絶対量、MB / 軸）。既定 16。0 で使わない**
    const std::size_t fine_mb = (argc > 11) ? std::strtoul(argv[11], nullptr, 10) : 16;
    o.ray_index_fine_bytes = fine_mb << 20;
    // **第 12 引数: O3 の判定（0 = 従来、1 = 整数の箱だけ、2 = 箱 + 厳密）。既定 0**
    // **引数が無ければライブラリの既定（2）を使う。** 0
    // を既定にすると既定の経路を測っていないことになる（実際に踏んだ）
    if (argc > 12) o.bsp_skip_disjoint = std::atoi(argv[12]);
    // **第 13 引数: 箱が片側なら飛ばす（1 = 既定、0 = 切る）**
    if (argc > 13) o.bsp_skip_boxside = std::atoi(argv[13]) != 0;
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
    // **やり直しは b ごとに別ファイル**。混ぜると意味が変わります
    const std::string done_path =
        redo ? (base + "_struct_b" + std::to_string(KRISITE_COORD_BITS) + ".txt")
             : (base + "_results.txt");
    std::vector<std::string> already;
    {
        std::ifstream f(base + "_results.txt");
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty()) already.push_back(line.substr(0, line.find(' ')));
        }
    }
    const auto seen = [&already](const std::string& k) {
        return std::find(already.begin(), already.end(), k) != already.end();
    };
    std::ofstream out(done_path, std::ios::app);

    // **★ 回す前に「これから何対か」を出します**（`CLAUDE.md`「対象の数も設定の一部」）。
    // **終わってから数えるのでは遅い。** 絞り込みの一覧が残っていて 2 対で終わった事例が
    // あります（2026-09-05）。
    {
        std::size_t planned = 0, skipped_seen = 0, skipped_only = 0, skipped_halt = 0;
        for (std::size_t k = 0; k + 1 < order.size(); k += 2) {
            const std::string key = ids[order[k]] + "x" + ids[order[k + 1]];
            if (!only.empty() && std::find(only.begin(), only.end(), key) == only.end()) {
                ++skipped_only;
                continue;
            }
            if (!redo && seen(key)) {
                ++skipped_seen;
                continue;
            }
            if (std::find(skip.begin(), skip.end(), key) != skip.end()) {
                ++skipped_halt;
                continue;
            }
            ++planned;
        }
        std::printf(
            "\n**これから %zu 対を回します**（全 %zu 対。既済 %zu / 一覧で除外 %zu / "
            "停止済み %zu）\n",
            planned, order.size() / 2, skipped_seen, skipped_only, skipped_halt);
        if (!only.empty()) {
            std::printf(
                "**★ `%s_only.txt` が %zu 件に絞っています。**全件回すなら退避してください\n",
                base.c_str(), only.size());
        }
    }
    // **★ 設定を全部出します**（`CLAUDE.md`「測定の出力に、設定と実際に測った対象を
    // 必ず書いてください」）。**2026-09-06 に、NSI の扱いが出力に無かったために、
    // 記録済みの 351 対（検査して宣言する = 1）に対して、続きを「宣言しない = 0」で
    // 回してしまいました。** 記録が版ではなく【設定】で割れており、
    // 出力の多角形数が 3.6 倍ずれて「案 1 の効果」に見えていました。
    static const char* kNsiName[4] = {"宣言しない", "検査して通ったものだけ宣言する",
                                      "宣言あり・なしの両方を回して突き合わせる",
                                      "検査せずに宣言する"};
    std::printf("\n## ブール演算（対 %zu、スレッド %u）\n\n", order.size() / 2, nthreads);
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 一覧 | `%s` |\n", list.c_str());
    std::printf("| 記録先 | `%s` |\n", done_path.c_str());
    std::printf("| 深度 | %u |\n", depth);
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| b（座標ビット） | %d |\n", KRISITE_COORD_BITS);
    std::printf("| **NSI の扱い** | **%d = %s** |\n", nsi_mode,
                (nsi_mode >= 0 && nsi_mode < 4) ? kNsiName[nsi_mode] : "?");
    std::printf(
        "| レイ索引の細かい割り当て K | %zu（0 = 従来）、上限 %zu 項目/三角形、絶対量 %zu MB/軸 "
        "|\n",
        fine_k, fine_budget, fine_mb);
    std::printf(
        "| O3（触れない三角形の平面で切らない） | %d（0 = 従来、1 = 整数の箱、2 = 箱 + 厳密） |\n",
        o.bsp_skip_disjoint);
    std::printf("| 箱が片側なら飛ばす（7 例目） | %d |\n", o.bsp_skip_boxside ? 1 : 0);
    std::printf("| 単一 source を割る閾値 P^2 | %zu |\n", o.single_src_sq);
    std::printf("| 索引の ON/OFF 突き合わせ | %s |\n", verify_index ? "する" : "しない");
    std::printf("| 済みの対 | %s |\n", redo ? "やり直す" : "飛ばす（再開）");
    // **★ 設定はすべて出します**（`HANDOVER.md` §6.2 の第三の条件）。
    // **記録済みと違う設定で回して、差を機能の効果と読み違えた事例があります。**
    std::printf("| **§5.5 の検算** | **%s**%s |\n", verify_delta ? "ON" : "**OFF**",
                verify_delta ? "（`SPEC-phase5.md` §3.2）"
                             : "（**§4.3.2 の EMBER 比較用。正しさの判定に使わないこと**）");
    std::printf("\n");
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t k = 0; k + 1 < order.size(); k += 2) {
        const std::size_t i = order[k], j = order[k + 1];
        const std::string key = ids[i] + "x" + ids[j];
        if (!only.empty() && std::find(only.begin(), only.end(), key) == only.end()) continue;
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
        PairStruct ps;
        const bool ok = check_one(prep[i].mesh, prep[j].mesh, o, &pool, &why, &h, &ps,
                                  nsi_mode == 3 ? 2 : (nsi_mode == 0 ? 0 : 1), verify_delta);
        if (ps.nsi_a >= 0) {
            (ps.nsi_a ? nsi_declared : nsi_rejected) += 1;
            (ps.nsi_b ? nsi_declared : nsi_rejected) += 1;
        }
        const double dt_first =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        // **NSI の宣言で出力が変わらないこと**（`SPEC-phase5.md` (c) の検証）。
        //
        // **新しい正解器は要りません。** 宣言しない側が従来の答えで、
        // **一致すればその入力で NSI が成り立っていた証拠**、
        // **不一致なら量子化がその模型の自己交差を作った**ことになります。
        if (nsi_mode == 2) {
            PairStruct ps0;
            unsigned long long h0 = 0;
            std::string why0;
            const auto t0n = std::chrono::steady_clock::now();
            const bool ok0 = check_one(prep[i].mesh, prep[j].mesh, o, &pool, &why0, &h0, &ps0, 0);
            const double dt0 =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0n).count();
            std::printf(
                "  NSI %s | バイト %s | P %zu -> %zu (%.2f 倍) | 秒 %.1f -> %.1f (%.2f 倍) "
                "| 省いたセル %zu | 宣言 A=%d B=%d(from_mesh %.2f 秒) | 判定 %s -> %s\n",
                key.c_str(), (h0 == h) ? "**一致**" : "**不一致**", ps0.polys, ps.polys,
                ps.polys ? double(ps0.polys) / double(ps.polys) : 0.0, dt0, dt_first,
                dt_first > 0 ? dt0 / dt_first : 0.0, ps.bsp_cells_skipped, ps.nsi_a, ps.nsi_b,
                ps.fm_seconds, ok0 ? "ok" : why0.c_str(), ok ? "ok" : why.c_str());
            std::fflush(stdout);
            // **ハッシュ不一致で止めません**（`IMPL-phase5.md` §20.2）。
            // NSI は切断を減らすので**三角形分割が変わります。**
            // **バイトでは守れない種類の変更**なので、判定（位相）で見ます。
            if (ok0 != ok) {
                std::printf("**判定が変わった** %s: %s -> %s\n", key.c_str(),
                            ok0 ? "ok" : why0.c_str(), ok ? "ok" : why.c_str());
                return 2;
            }
            continue;
        }
        // **索引の有無で出力が変わらないこと**（`SPEC-phase5.md` §6.3 の安全弁）。
        // 索引は厳密な絞り込みなので、**バイトで一致しなければ欠陥**です。
        if (ok && verify_index) {
            csg::BoolOptions o2 = o;
            o2.ray_index = !o.ray_index;
            unsigned long long h2 = 0;
            std::string why2;
            // **NSI の扱いを揃えます。** 揃えないと「索引の比較」のつもりで
            // **NSI の有無を比べる**ことになり、番人が黙って無意味になります。
            const bool ok2 = check_one(prep[i].mesh, prep[j].mesh, o2, &pool, &why2, &h2, nullptr,
                                       nsi_mode == 3 ? 2 : (nsi_mode == 0 ? 0 : 1));
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
        out << hb.data() << ' ';
        ps.print(out);
        out << ' ' << why << '\n';
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
    // **機構が空回りしていないことを別に検査します**（`CLAUDE.md`）。
    // 検査つきの宣言では、**宣言できた数と落ちた数の両方が出ます。**
    // 片方が 0 なら、検査が効いていないか、母集団が偏っています。
    if (nsi_mode == 1 || nsi_mode == 2) {
        const std::size_t tot = nsi_declared + nsi_rejected;
        std::printf("**NSI: 宣言 %zu / 却下 %zu（%zu 模型ぶん、%.1f%%）**\n", nsi_declared,
                    nsi_rejected, tot, tot ? 100.0 * double(nsi_declared) / double(tot) : 0.0);
    }
    return c.failed == 0 ? 0 : 1;
}
