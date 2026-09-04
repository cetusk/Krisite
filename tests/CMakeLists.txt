# Krisite — テスト
#
# 外部テストフレームワークは使わない。各テストは独自 main を持ち、
# 失敗時に非零で終了することで ctest に結果を伝える。

add_library(krisite_test_util INTERFACE)
target_include_directories(krisite_test_util INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(krisite_test_util INTERFACE krisite)

function(krisite_add_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE krisite krisite_test_util)
    add_test(NAME ${name} COMMAND ${name} ${ARGN})
    set_tests_properties(${name} PROPERTIES TIMEOUT 3600)
endfunction()

# ---- 除外ラベル `heavy_geometry`（CI の副構成で外すもの）---------------------
#
# **印を付ける向きに注意してください**（`CLAUDE.md`「ラベルで検査範囲を絞るときは、
# 失敗の向きを安全側にすること」）。
#
#   対象に印   … 付け忘れると**テストが走らない。網が黙って縮む**
#   対象外に印 … 付け忘れると**テストが走る。時間が少し増えるだけ**  ← こちらを採る
#
# ここで付けるのは**外してよい側**です。**新しいテストは、何もしなければすべての
# 構成で走ります。**
#
# 外してよい根拠があるのは 2 つの構成だけです。
#
#   b=26 の副構成   b=26 が固有に見るのはリム数の違いで、プラットフォームには依りません。
#                   **Linux Clang の b=26 では全部走ります**（幾何配置での上界超過の網）
#   可搬経路        変わるのは arith の実装だけ。直接の網は fixed_int と GMP 差分 10^7 件
#
# **`widths.hpp` を参照するテストは入れません**（`topology` / `plane_order` /
# `bitwidth` / `volume_gmp`）。幅に直接依存するので b=26 で見る価値があります。
function(krisite_mark_heavy)
    foreach(t ${ARGN})
        set_tests_properties(${t} PROPERTIES LABELS "heavy_geometry")
    endforeach()
endfunction()

# ---- GMP を必要としないテスト（既定で常に実行）------------------------------

krisite_add_test(fixed_int      arith/test_fixed_int.cpp)
krisite_add_test(bitwidth       arith/test_bitwidth.cpp)
krisite_add_test(predicates     geom/test_predicates.cpp)
krisite_add_test(degenerate     geom/test_degenerate.cpp)
krisite_add_test(plane_order    geom/test_plane_order.cpp)
krisite_add_test(plane_box      geom/test_plane_box.cpp)
krisite_add_test(topology      mesh/test_topology.cpp)
krisite_add_test(uniform_grid  octree/test_uniform_grid.cpp)
krisite_add_test(faces         csg/test_faces.cpp)
krisite_add_test(fragment      csg/test_fragment.cpp)
krisite_add_test(cp1           csg/test_cp1.cpp)
# SPEC-phase2 §9.1 / CP1。絞り込みの有無を同一プロセスで比較する
krisite_add_test(split_strategy csg/test_split_strategy.cpp)
# SPEC-phase2 §2.4.4 (3) の【正の対照】。負の対照（固定深度で挿入 0 件）は
# split_strategy が持つ。**負の対照だけでは何もしないスタブが通る。**
krisite_add_test(tjunction     csg/test_tjunction.cpp)
# SPEC-phase2 §3 / §9.1 / CP2。固定深度（Phase 1 の挙動）を正解器に、適応分割と比較する
krisite_add_test(adaptive      csg/test_adaptive.cpp)
# SPEC-phase2 §3.2 / §3.3 / CP2。**early-out 無効側が正解器**。省いても答えが
# 変わらないことを直接検査する
krisite_add_test(early_out     csg/test_early_out.cpp)
# SPEC-phase2 §4 / CP3。**キャッシュ無効側が正解器**。値まで一致することを要求する
krisite_add_test(cache         csg/test_cache.cpp)
# SPEC-phase2 §5 / CP4。**§5.5.1 の壊れ方 2 通りを両方見る**
#   誤った組で分裂が起きる   → C 不変性
#   誤った組で分裂が起きない → 分裂後の多様体性 + ΔE の予測不一致
krisite_add_test(split_semantics csg/test_split_semantics.cpp)
krisite_add_test(cp2           csg/test_cp2.cpp)
krisite_add_test(cp25          csg/test_cp25.cpp)
krisite_add_test(cp3           csg/test_cp3.cpp)
# **レイキャストの 2 次元索引が出力を変えないこと**（SPEC-phase5 §6.3、CP1.5-2）。
# 索引は厳密な絞り込みなので、出力はバイトで一致しなければならない。
# **索引が実際に検査数を減らした回数も数え、0 なら落とす**（空回りの検査）
krisite_add_test(ray_index     csg/test_ray_index.cpp)
# **NSI の宣言が意味論を変えないこと**（SPEC-phase3 §5.6、SPEC-phase5 (c)）。
# 三角形分割は変わるのでバイトでは守れない。位相で見る（体積は volume_gmp）。
# **ケース 26 が唯一 NSI を発火させる**ので、空回りの番人を置いてある
krisite_add_test(nsi           csg/test_nsi.cpp)
# **自己交差の検査**（SPEC-phase3 §5.6）。NSI を宣言してよいかを機械的に確かめる。
# 健全側に倒した検査なので、**意図的に自己交差させた入力が落ちること**を番人にする
krisite_add_test(self_intersect csg/test_self_intersect.cpp)
# **浮動小数点の体積恒等式の校正**（SPEC-phase5 §3.0）。閾値を実測から決めるための
# 分布を出し、丸めの床が上がっていないことを見る回帰の番人でもある
krisite_add_test(fp_volume_calib csg/fp_volume_calib.cpp)
# SPEC-phase2 §13 の CP5。**機構どうしの相互作用**（適応分割 × 分裂 /
# 構成点の保持 × T 解決 / early-out × 分裂）。3 経路を踏んだ回数を数え、
# 0 なら落とす。§9.4.2 の「分裂の有無で幾何が一致する」もここが常時見る
krisite_add_test(cp5           csg/test_cp5.cpp)
# **Phase 3 の CP5**（n 項 / スープ経路）。Phase 2 の `cp5` は二項経路のままです。
krisite_add_test(cp5n          csg/test_cp5n.cpp)
# **n 項のケース 17 / 18**（SPEC-phase3 §9）。別の演算木で同じ立体を作る恒等式。
krisite_add_test(nary          csg/test_nary.cpp)
# SPEC-phase3 §5.2.0 / §9 のケース 27。**「内側」の定義（w > 0）**。
# 既存の検査（位相・体積・分割戦略・決定性・正解器）がすべて素通りするので、
# **絶対値で検査します**（§9.1）。
krisite_add_test(winding csg/test_winding.cpp)
# **凸分割**（SPEC-phase3 §4.1）。三角形化との比較と、決定性。
krisite_add_test(convex        csg/test_convex.cpp)
# **並列化**（SPEC-phase4 §7）。決定性が主検出器。
krisite_add_test(parallel      csg/test_parallel.cpp)
# **この検査は CPU を要求します**（`SPEC-phase4.md` §7.3）。
#
# 8 スレッドのプールが**実際に並列に走ったか**を統計で確かめる検査があり、
# 走らなければ「共有の有無を見ていない」ので落とします。**空回りを許さない設計**です。
#
# **`ctest --parallel 4` で他のテストと同時に走ると、4 コアのランナーでは
# プールのワーカーが起動する前に呼び出し元がキューをさらいます。**
# 実際に Windows の CI で落ちました（`IMPL-phase5.md` §24）。
#
# **検査を緩めるのは誤りです** — 「並列に走らなかった」という報告自体は正しい。
# **必要な CPU を予約するのが正しい対処**です。ctest はこの値を見て、
# 他のテストと同時に走らせません。
set_tests_properties(parallel PROPERTIES PROCESSORS 4)
# SPEC-phase3 §3.1 / §2.1 の CP1。辺平面の構成（軸方向との外積）と、
# 断片の相対内部の点。**両方の経路が使われることを固定する**（空回りの番人）
krisite_add_test(interior      csg/test_interior.cpp)
# SPEC-phase3 §14 の CP2。入口・中核・出口の分離。**二項正解器との一致**と、
# **連鎖でビット幅が伸びないこと**（§10.3.1。契約の成立条件）
krisite_add_test(soup          csg/test_soup.cpp)
krisite_add_test(idempotence   csg/test_idempotence.cpp)
krisite_add_test(depth_monotonic csg/test_depth_monotonic.cpp)

# §12 の呼び出し比率。**この実行ファイルだけ**計数を有効にする
# （既定ビルドに計数のコストを持ち込まない）。
krisite_add_test(call_ratio    csg/test_call_ratio.cpp)
target_compile_definitions(call_ratio PRIVATE KRISITE_COUNT_PREDICATES=1)

# **重い幾何テストに「外してよい」印を付ける。** 根拠は krisite_mark_heavy の上の注記。
krisite_mark_heavy(
    uniform_grid faces fragment
    cp1 cp2 cp25 cp3 cp5
    split_strategy tjunction adaptive early_out cache split_semantics
    idempotence depth_monotonic call_ratio)

# ---- GMP 差分テスト（LGPL、テスト専用）--------------------------------------

if(KRISITE_GMP_FOUND)
    # SPEC §9: 算術 10^7 件、述語 10^6 件以上
    krisite_add_test(fixed_int_gmp  arith/test_fixed_int_gmp.cpp 10000000)
    target_link_libraries(fixed_int_gmp PRIVATE krisite_gmp)

    krisite_add_test(predicates_gmp geom/test_predicates_gmp.cpp 1000000)
    target_link_libraries(predicates_gmp PRIVATE krisite_gmp)

    krisite_add_test(plane_order_gmp geom/test_plane_order_gmp.cpp 300000)
    target_link_libraries(plane_order_gmp PRIVATE krisite_gmp)

    # SPEC-phase2 §2.3 / §7。正解器は 8 隅すべてを mpz で評価する別経路
    krisite_add_test(plane_box_gmp  geom/test_plane_box_gmp.cpp 200000)
    target_link_libraries(plane_box_gmp PRIVATE krisite_gmp)

    # §10.3 体積の恒等式。出力の頂点は有理数なので GMP が要る
    krisite_add_test(volume_gmp     csg/test_volume_gmp.cpp)
    target_link_libraries(volume_gmp PRIVATE krisite_gmp)

    # 開発中に素早く回すための短縮版
    add_test(NAME fixed_int_gmp_quick COMMAND fixed_int_gmp 40000)
    add_test(NAME predicates_gmp_quick COMMAND predicates_gmp 5000)
    add_test(NAME plane_order_gmp_quick COMMAND plane_order_gmp 3000)
    add_test(NAME plane_box_gmp_quick COMMAND plane_box_gmp 3000)
endif()

# ---- Manifold 正解器（Apache-2.0、テスト専用）--------------------------------
#
# **ライブラリ本体には絶対にリンクしません**（THIRD_PARTY_LICENSES.md）。この実行ファイルだけです。

if(KRISITE_MANIFOLD_FOUND)
    krisite_add_test(manifold_oracle csg/test_manifold.cpp)
    target_link_libraries(manifold_oracle PRIVATE manifold::manifold)
endif()

# ---- 変異テスト（SPEC-phase1.md §10.5）--------------------------------------
#
# **「変異テストがある」ことと「変異テストが効いている」ことは別です**（§10.5）。
# そこで期待値を表として書き下し、ctest に検証させます。
#
#   expect = FAIL … その変異をそのテストが検出すること（`WILL_FAIL`）
#   expect = PASS … **検出しないこと**。これも固定します。素通りする組合せを
#                   記録しておかないと、後から「変異テストがあるから安心」と
#                   誤解されます
#
# 変異 3 は**2 つの独立した機構**が捕まえます。検査ビルドでは `fragment_sign` の
# 不変条件が停止させ、検査なしビルドでは §10.1 の位相検査が落とします。
# ここで測るのは後者です（前者はライブラリ本体の防御であってコーパスの検出力ではない）。

if(KRISITE_BUILD_MUTANTS)
    # **変異テストは検査 OFF で回してください。**
    #
    # 検査 ON だと変異 3 は `KRISITE_CHECK`（`fragment_sign` の不変条件）で abort します。
    # それも立派な検出ですが、**測っているものが違います。** 変異テストが答えるべきは
    # 「テスト体系がこの誤りを検出するか」であって「表明が止めるか」ではありません。
    # 表明はライブラリ本体の防御であり、コーパスの検出力とは独立です。
    #
    # 加えて CTest は abort（シグナル終了）を `WILL_FAIL` で拾えません。
    if(KRISITE_CHECKED_ARITH OR CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(FATAL_ERROR
            "KRISITE_BUILD_MUTANTS には KRISITE_CHECKED_ARITH=OFF と非 Debug 構成が要ります。"
            " 検査 ON では変異 3 が表明で停止し、テスト体系の検出力を測れません"
            "（SPEC-phase1 §10.5）。")
    endif()

    # 5 番目の引数に GMP を渡すと GMP をリンクする（体積検査の変異用）。
    #
    # **ラベルを分けます**（SPEC-phase2 §9.4.1）。GMP 版は `mutation_gmp` だけを持つので、
    # **`mutation` の件数は GMP の有無で動きません。** 件数の番人が構成に依存しなくなります。
    function(krisite_add_mutant tag define test expect)
        set(name "mut_${tag}_${test}")
        add_executable(${name} csg/test_${test}.cpp)
        target_link_libraries(${name} PRIVATE krisite krisite_test_util)
        set(label "mutation")
        if("${ARGV4}" STREQUAL "GMP")
            target_link_libraries(${name} PRIVATE krisite_gmp)
            set(label "mutation_gmp")
        endif()
        target_compile_definitions(${name} PRIVATE ${define})
        add_test(NAME ${name} COMMAND ${name})
        set_tests_properties(${name} PROPERTIES TIMEOUT 3600 LABELS "${label}")
        # **本体と同じ資源を予約します。**
        #
        # `parallel` の §7.3 / §7.1 は「8 スレッドで実際に並列に走ったか」に依存します。
        # 本体にだけ付けて変異側に付け忘れると、**実行条件が揃いません。**
        #
        # **ただし、これは変異 23 の不安定さの対策ではありません。**
        # 変異 23 の検出率は**実測で 55〜60%** で、`92a477a`（CI 緑）でも同じでした
        # （`IMPL-phase5.md` §24.5）。**資源を与えても解決しません。**
        if(test STREQUAL "parallel")
            set_tests_properties(${name} PROPERTIES PROCESSORS 4)
        endif()
        if(expect STREQUAL "FAIL")
            set_tests_properties(${name} PROPERTIES WILL_FAIL TRUE)
        endif()
    endfunction()

    # **主検出器も固定します**（SPEC-phase2 §9.3.1）。
    #
    # 「素通りする組合せを固定する」の自然な拡張です。**通ることだけでなく、どれが
    # 捕まえているかを固定します。** 検出器が黙って移動したら、この表が落ちて可視化されます。
    # 実際、案 D の導入で変異 3 / 4 の主検出器が位相検査から §9.1 に移動しました（§2.4.6）。
    #
    # **表には「CI で検証されているか」も併記します**（§9.3.1）。CI で回らない検出器が
    # 混ざっていると、**表そのものが「緑だが検証していない」状態**になります。
    #
    #   変異                       主検出器                       CI
    #   ------------------------   ----------------------------   --------------------
    #    1 第2段の無効化           位相（cp25 / cp3）             mutation
    #    2 割り当てを開区間に       位相（cp2 / cp25 / cp3）       mutation
    #    3 分割平面を局所に         split_strategy（§9.1）★       mutation
    #    4 到達三角形の取り違え     split_strategy（§9.1）★       mutation
    #    6 early-out を緩める       early_out（有無の比較）        mutation + mutation_gmp
    #    7 キャッシュのキー落ち     cache（値の完全一致）          mutation + mutation_gmp
    #    8 分裂の無効化             split_semantics（多様体性）    mutation + mutation_gmp
    #    9 対応付けを owner に      split_semantics（深度で交替）  mutation + mutation_gmp
    #   10 T 解決の無効化           位相（adaptive。次数 1）       mutation
    #   11 T を 1 個ずつ            位相（adaptive。次数 1）       mutation
    #   12 割り当てを閉領域に       位相（adaptive。C の破裂）     mutation
    #   13 退化三角形を捨てる       位相（adaptive。次数 1）       mutation
    #   14 辺平面の軸選択を非正準に interior（正準性）★           mutation
    #
    # **CI で回らない登録はありません。** `mutation`（GMP OFF）と `mutation_gmp`
    # （§9.4.1）の 2 ジョブが、上のすべてを回します。

    # 変異 1: §5.3 の第2段（値ベースの併合）を無効化
    #   **主検出器: 位相検査（cp25 / cp3）**
    #
    # **軸平行コーパスでは素通りします。** 値併合が 0 件（第1段が全頂点を捕まえる）
    # なので、第2段を消しても何も変わりません。斜面コーパスが必須です。
    krisite_add_mutant(1 KRISITE_MUTATION_NO_STAGE2=1 cp25 FAIL)
    krisite_add_mutant(1 KRISITE_MUTATION_NO_STAGE2=1 cp3  FAIL)
    krisite_add_mutant(1 KRISITE_MUTATION_NO_STAGE2=1 cp2  PASS)
    # **冪等性は 3 変異とも検出しません。** 比較側が頂点を値で同一視するので、
    # 第2段を止めても一致します（構造的な盲点）。∩ と ∪ を互いに比べる形なので、
    # 両方に同じ変異がかかると差が出ません。**素通りを固定して記録します。**
    krisite_add_mutant(1 KRISITE_MUTATION_NO_STAGE2=1 idempotence PASS)

    # 変異 2: §4.2 の割り当てを開領域に変える
    #   **主検出器: 位相検査（cp2 / cp25 / cp3）**
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp2             FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp25            FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp3             FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 depth_monotonic FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 idempotence     PASS)

    # 変異 3: §4.3.1 の分割平面を割り当て集合に変える（T 字接合が生じる）
    # 変異 4: §2.3 の絞り込みの基準を「三角形がセルに届くか」に戻す
    #
    #   **主検出器: split_strategy（§9.1 の分割戦略不変性）** ★ 位相検査から移動しました
    #
    # **案 D（T 頂点の解決）が、この 2 つが作る T 字接合を埋めてしまいます**（§2.4.6）。
    # Phase 1 では位相検査が主検出器でしたが、いまは素通りします。
    #
    # **原因は「T 解決が覆い隠している」ことで確定しています**（§9.3.2 の一度きりの診断）。
    # 変異 10（T 解決の無効化）と組み合わせると cp3 が 218/2085 件で落ちます。
    # **コーパスは突く能力を保持しています。**
    #
    # **`adaptive` も素通りします。** §9.1 を $(C,\chi)$ で比べる側は、
    # 「同じ立体が別の三角形分割で出る」ことを許すためです。**検出するのは
    # $(V,E,F)$ の完全一致を要求する `split_strategy` の側だけです。**
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 split_strategy  FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 depth_monotonic FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 adaptive        PASS)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp1             PASS)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp2             PASS)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp25            PASS)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp3             PASS)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 idempotence     PASS)

    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 split_strategy  FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 depth_monotonic FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 adaptive        PASS)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp1             PASS)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp2             PASS)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp25            PASS)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp3             PASS)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 idempotence     PASS)

    # 変異 6: §3.2 の early-out の判定を緩め、**相手が居るセルでも省く**。
    #   **主検出器: early_out（early-out 有無の比較）**
    #
    # **検出器は 2 つあり、役割が違います**（SPEC-phase1 §10.2 の但し書き）。
    #
    #   early_out   … 位相。省いたことで答えが変わっていないか（268 件）
    #   volume_gmp  … **体積。断片の選択の誤り**（35 件）。位相は「割れているか」しか見ない
    #
    # **early_out 無効のテストは素通りします。** 変異が発火しないので当然です。
    krisite_add_mutant(6 KRISITE_MUTATION_LOOSE_EARLY_OUT=1 early_out       FAIL)
    krisite_add_mutant(6 KRISITE_MUTATION_LOOSE_EARLY_OUT=1 adaptive        PASS)
    krisite_add_mutant(6 KRISITE_MUTATION_LOOSE_EARLY_OUT=1 cp3             PASS)
    krisite_add_mutant(6 KRISITE_MUTATION_LOOSE_EARLY_OUT=1 depth_monotonic PASS)
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(6 KRISITE_MUTATION_LOOSE_EARLY_OUT=1 volume_gmp FAIL GMP)
    endif()

    # 変異 7: §4.2 の構成点キャッシュのキーから 1 平面を落とす。
    #   **主検出器: cache（キャッシュ有無の比較。値の完全一致）**
    #
    # **異なる3つ組が同じキーに落ち、別の点を同一視します。** 静かに壊れる種類の誤りです。
    #
    #   cache       … 出力の**値**が変わる（1101 件）
    #   volume_gmp  … **体積。別の立体を作っている**（97 件）
    #
    # **キャッシュ無効のテストは素通りします。** 変異が発火しないので当然です。
    krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 cache     FAIL)
    krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 nary      FAIL)
    krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 adaptive  PASS)
    krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 early_out PASS)
    krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 cp3       PASS)
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(7 KRISITE_MUTATION_CACHE_KEY_DROP=1 volume_gmp FAIL GMP)
    endif()

    # 変異 8: §5 の接触の分裂を無効化する（辺・頂点の両方）。
    #   **主検出器: split_semantics（位相。多様体性）**
    #
    # **8a / 8b の 2 変種も用意します**（§9.3）。**11b は辺分裂でしか直りません。**
    #
    #   8   両方を無効        FAIL 28 / 2359
    #   8a  辺の分裂だけ無効   FAIL 45 / 2359（= 頂点分裂だけ。11b が直らない）
    #   8b  頂点の分裂だけ無効 FAIL 22 / 2359（= 辺分裂だけ。4T / 4T′ が直らない）
    #
    # **体積は分裂の誤りを検出しません**（§5.1.3 のとおり体積不変）。位相でしか出ません。
    # cp3 / adaptive は分裂 OFF で回るので素通りします。
    krisite_add_mutant(8 KRISITE_MUTATION_NO_SPLIT=1 split_semantics FAIL)
    krisite_add_mutant(8 KRISITE_MUTATION_NO_SPLIT=1 cp3             PASS)
    krisite_add_mutant(8 KRISITE_MUTATION_NO_SPLIT=1 adaptive        PASS)
    krisite_add_mutant(8a KRISITE_MUTATION_NO_EDGE_SPLIT=1 split_semantics FAIL)
    krisite_add_mutant(8a KRISITE_MUTATION_NO_EDGE_SPLIT=1 cp3             PASS)
    krisite_add_mutant(8b KRISITE_MUTATION_NO_VERTEX_SPLIT=1 split_semantics FAIL)
    krisite_add_mutant(8b KRISITE_MUTATION_NO_VERTEX_SPLIT=1 cp3             PASS)

    # 変異 9: §5.1.2 の対応付けを owner に戻す。
    #   **主検出器: split_semantics。ただし深度で入れ替わります**（§5.5.1 の訂正）
    #
    #   深度 0      分裂が起きない  → 多様体性 / ΔE（予測 4・実測 0）/ Δχ / unresolved
    #   深度 1〜3   誤って分裂する  → **C 不変性のみ**（C が 4 → 1）
    #
    # **片方の検査だけでは深度の半分を取り逃がします。** 9 件の失敗の内訳がこれです。
    #
    # **ケース 16（自己接触）でのみ発火します**（§8.1）。13 / 14 / 15 は自己接触を
    # 作らないので、owner でも連結成分でも同じ結果になり**素通りします。**
    #
    # **$\chi$ でも体積でも捕まりません。** cp3 / adaptive は素通り。体積は 0 件
    # （§5.1.3。分裂は体積を変えません）。**volume_gmp は GMP ビルドでのみ登録されます。**
    krisite_add_mutant(9 KRISITE_MUTATION_SPLIT_BY_OWNER=1 split_semantics FAIL)
    krisite_add_mutant(9 KRISITE_MUTATION_SPLIT_BY_OWNER=1 cp3             PASS)
    krisite_add_mutant(9 KRISITE_MUTATION_SPLIT_BY_OWNER=1 adaptive        PASS)

    # **分裂の変異は体積では捕まりません**（§5.1.3 の定理。§9.4.1 で素通りを固定する）。
    # **不変量そのものは volume_gmp の §9.4.2 の検査が直接見ています。**
    # ここで固定するのは「体積は分裂の誤りの検出器ではない」という役割分担です。
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(8  KRISITE_MUTATION_NO_SPLIT=1         volume_gmp PASS GMP)
        krisite_add_mutant(8a KRISITE_MUTATION_NO_EDGE_SPLIT=1    volume_gmp PASS GMP)
        krisite_add_mutant(8b KRISITE_MUTATION_NO_VERTEX_SPLIT=1  volume_gmp PASS GMP)
        krisite_add_mutant(9  KRISITE_MUTATION_SPLIT_BY_OWNER=1   volume_gmp PASS GMP)
    endif()

    # 変異 10: §2.4.3 の T 頂点解決を無効化する（SPEC-phase2 §9.3）。
    #   **主検出器: 位相検査（adaptive。辺の次数 1）**
    #
    # **T 解決は継ぎ目の正しさの単一障害点です。** 検出の経路は 3 つあり、
    # **測っているものが違います。**
    #
    #   adaptive       … 位相で検出（186 件）。**T 字接合がそのまま残る**
    #   tjunction      … 正の対照。幾何で検出（効くべき場所で効かない）
    #   split_strategy … 負の対照の【非空虚の番人】。計数で検出（索引が候補を返さない）
    #
    # **固定深度のコーパス（cp3 / depth_monotonic）は素通りします。** 固定深度では
    # T 解決がそもそも no-op なので、無効化しても出力が変わりません。
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 adaptive        FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 tjunction       FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 split_strategy  FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 cp3             PASS)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 depth_monotonic PASS)

    # 変異 11: T 頂点を 1 本の辺に 1 個しか入れない（一括でなく）。
    #   **主検出器: 位相検査（adaptive。辺の次数 1）**
    #
    # **1 本の辺に T 頂点が 2 個以上載る配置でしか検出できません**（§8 のケース 13 の要件）。
    # 適応分割でケース 13 がその配置を作るので、位相で検出されます（35 件）。
    # **固定深度は T 頂点が 0 個なので構造的に素通りします。**
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 adaptive       FAIL)
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 tjunction      FAIL)
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 split_strategy PASS)
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 cp3            PASS)

    # 変異 12: §3.4 の割り当てを閉領域に戻す。
    #   **主検出器: 位相検査（adaptive。C の破裂）**
    #
    # **適応分割 × 面がセル境界に乗るケースでしか検出できません。** 粗い側の 1 枚と
    # 細かい側の 4 枚が入れ子になり、併合されずに同じ面積が二重に出力されて C が破裂します。
    # **固定深度は素通りします**（両側の断片が一致して `region_key` が併合するため）。
    #
    # `depth_monotonic` は別経路で捕まえます。規律違反の配置で**生の断片数が増える**ように
    # なるためで、§3.4 の「重複が生じなくなった」ことの裏返しです。
    krisite_add_mutant(12 KRISITE_MUTATION_CLOSED_CELLS=1 adaptive        FAIL)
    krisite_add_mutant(12 KRISITE_MUTATION_CLOSED_CELLS=1 depth_monotonic FAIL)
    krisite_add_mutant(12 KRISITE_MUTATION_CLOSED_CELLS=1 split_strategy  PASS)
    krisite_add_mutant(12 KRISITE_MUTATION_CLOSED_CELLS=1 cp3             PASS)

    # 変異 13: 退化三角形を「作らない」のではなく「捨てる」（§2.4.4 (2)）。
    #   **主検出器: 位相検査（adaptive。辺の次数 1）**
    #
    # **捨てると境界辺がメッシュから消え、次数 1 の辺が残ります。**
    # **面積は保たれるので体積検査（volume_gmp）では捕まりません。** 位相でしか出ません。
    # 固定深度は T 頂点が 0 個なので素通りします。
    krisite_add_mutant(13 KRISITE_MUTATION_DROP_DEGENERATE=1 adaptive        FAIL)
    krisite_add_mutant(13 KRISITE_MUTATION_DROP_DEGENERATE=1 tjunction       FAIL)
    krisite_add_mutant(13 KRISITE_MUTATION_DROP_DEGENERATE=1 split_strategy  PASS)
    krisite_add_mutant(13 KRISITE_MUTATION_DROP_DEGENERATE=1 cp3             PASS)

    # 変異 15: 出力の向きを **`f.flipped` を見ずに `in_front` だけ**で決める（CP4 で修正）。
    #   **主検出器: volume_gmp（体積。符号が狂う）**
    #
    # **これは実際に入っていた誤りです**（`IMPL-phase3.md` §6.1）。断片の巻き順は
    # 元の三角形の法線基準なので、`f.flipped` が真の面だけ向きが反転します。
    #
    #   volume_gmp … 体積。**符号が狂う**（`A ∪ A` で −1/3 になりました）
    #   soup       … 位相。`oriented`（各辺で #(u,v) = #(v,u)）が崩れる
    #
    # **(C, χ) は面の向きを反転しても変わりません。** 正解器との (C, χ) 比較は
    # **素通りします** — それがこの誤りが CP2 から CP3 まで残った理由です。
    krisite_add_mutant(15 KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED=1 soup            FAIL)
    krisite_add_mutant(15 KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED=1 nary            FAIL)
    krisite_add_mutant(15 KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED=1 cp3             PASS)
    krisite_add_mutant(15 KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED=1 depth_monotonic PASS)
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(15 KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED=1 volume_gmp FAIL GMP)
    endif()

    # 変異 16: §5.4 の局所 BSP が **1 枚も切らない**。
    #   **主検出器: soup（局所 BSP ↔ 過剰分割の比較）**
    #
    # 断片の相対内部をシートが横切るので、**代表点 1 点で分類できる前提が崩れます**
    # （§7.2 の仮定）。位相と体積の両方に出ます。
    #
    # **二項実装は過剰分割のままなので、`cp3` などは素通りします。**
    krisite_add_mutant(16 KRISITE_MUTATION_BSP_NEVER_CUT=1 soup FAIL)
    krisite_add_mutant(16 KRISITE_MUTATION_BSP_NEVER_CUT=1 nary FAIL)
    krisite_add_mutant(16 KRISITE_MUTATION_BSP_NEVER_CUT=1 cp3  PASS)
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(16 KRISITE_MUTATION_BSP_NEVER_CUT=1 volume_gmp FAIL GMP)
    endif()

    # 変異 17: 「この曲面がこのセルに存在するか」を、**割り当て用の半開区間の述語**で
    # 判定する（`assign_to_cell`）。
    #   **主検出器: cp5n（early-out 有無の一致）**
    #
    # **これは実際に入っていた誤りです**（`IMPL-phase3.md` §7.1）。面がセルの上側境界に
    # ちょうど乗る配置でしか出ません。存在の問いは**閉領域**（`overlaps_cell`）です。
    #
    #   cp5n … early-out 有無の一致が崩れる。n 項 3 source × 適応分割で非多様体になる
    #
    # **二項実装は別経路なので `cp3` などは素通りします。**
    krisite_add_mutant(17 KRISITE_MUTATION_HALF_OPEN_PRESENCE=1 cp5n FAIL)
    krisite_add_mutant(17 KRISITE_MUTATION_HALF_OPEN_PRESENCE=1 cp3  PASS)
    # **`nary` も素通りします。** 面がセルの上側境界に乗る配置がケース 17 / 18 に
    # 無いためで、**`cp5n` が唯一の検出器**であることの記録です。
    # ここが FAIL に転じたら検出器が増えたということなので、表を更新してください。
    krisite_add_mutant(17 KRISITE_MUTATION_HALF_OPEN_PRESENCE=1 nary PASS)

    # 変異 18: **局所 BSP の交差線分を 1 本落とす**（`SPEC-phase3.md` §10.5）。
    #   **主検出器: soup（局所 BSP ↔ 過剰分割の比較）**
    #
    # 落とすと断片が相手の曲面をまたいだまま残り、**代表点 1 点で分類できる前提
    # （定理 7.2）が崩れます。** 位相にも体積にも出ます。
    #
    #   soup       … 局所 BSP と過剰分割で (C, χ) が食い違う
    #   cp5n       … 全機構有効の全域
    #   nary       … 恒等式の両方向の差が空にならない
    #   volume_gmp … 体積。**別の立体を作っている**
    #
    # **二項実装は過剰分割のままなので `cp3` は素通りします。**
    #
    # > 仕様の「保守的分割をやめる」は **retire** しました（§10.5.1）。`Fragment` が
    # > 凸専用なので、交差線分だけで切る defect は**表現できません。** これが代替です。
    krisite_add_mutant(18 KRISITE_MUTATION_BSP_DROP_ONE_CUT=1 soup FAIL)
    krisite_add_mutant(18 KRISITE_MUTATION_BSP_DROP_ONE_CUT=1 cp5n FAIL)
    krisite_add_mutant(18 KRISITE_MUTATION_BSP_DROP_ONE_CUT=1 nary FAIL)
    krisite_add_mutant(18 KRISITE_MUTATION_BSP_DROP_ONE_CUT=1 cp3  PASS)
    if(KRISITE_GMP_FOUND)
        krisite_add_mutant(18 KRISITE_MUTATION_BSP_DROP_ONE_CUT=1 volume_gmp FAIL GMP)
    endif()

    # 変異 19: **凸分割で凸性を確かめずに併合する**（`SPEC-phase3.md` §4.1.2）。
    #   **主検出器: convex（点集合として同じか）**
    #
    # `Fragment` は半空間の交わりなので、非凸なループを渡すと**凸包**になります。
    # **位相は保たれることがあります**（$\chi$ も $C$ も変わらない）。
    # 「三角形化と凸分割で両方向の差が空」でしか出ません。
    #
    # **既定は三角形化なので、`soup` や `cp5n` は素通りします。**
    krisite_add_mutant(19 KRISITE_MUTATION_MERGE_IGNORE_CONVEXITY=1 convex FAIL)
    krisite_add_mutant(19 KRISITE_MUTATION_MERGE_IGNORE_CONVEXITY=1 soup   PASS)
    krisite_add_mutant(19 KRISITE_MUTATION_MERGE_IGNORE_CONVEXITY=1 cp5n   PASS)

    # 変異 20: **結合時の正準な整列を外す**（`SPEC-phase4.md` §7.5）。
    #   **主検出器: parallel（決定性）**
    #
    # スロットに書かず、完了した順に共有の器へ積みます。排他は取るので**競合は
    # ありません**が、並びがスケジューリングに依存します。
    #
    # **幾何は同じなので、位相も体積も区別しません。** `SPEC-phase3.md` の
    # 「辺平面の軸選択」と同じ形で、**決定性の検査でしか捕まりません。**
    krisite_add_mutant(20 KRISITE_MUTATION_JOIN_ORDER=1 parallel FAIL)
    krisite_add_mutant(20 KRISITE_MUTATION_JOIN_ORDER=1 soup     PASS)
    krisite_add_mutant(20 KRISITE_MUTATION_JOIN_ORDER=1 nary     PASS)
    # 変異 21 / 22: **キャッシュと統計をスレッド間で共有する**（§7.5）。
    #   **主検出器: ThreadSanitizer（CI の `tsan` ジョブ）**
    #
    # **通常のビルドでは登録しません。** 21 は `std::map` への同時挿入なので
    # **落ち方が未定義**（実測では abort）で、`WILL_FAIL` はシグナル終了を拾えません。
    # 22 は出力に現れないので、**決定性でも位相でも体積でも捕まりません。**
    #
    # 実測（TSan、8 スレッド、コーパス全件）:
    #
    #   KRISITE_MUTATION_SHARE_CACHE   警告 164 件
    #   KRISITE_MUTATION_SHARE_STATS   警告  45 件
    #
    # **「TSan でしか捕まらない」ことが、TSan を置く理由そのものです。**

    # 変異 23: **出口の段の間のバリアを 1 つ外す**（`SPEC-phase4.md` §7.5）。
    #   **主検出器: parallel（決定性）**
    #
    # 扇の計算（並列）と ID の割り当て（逐次）を融合し、**完了した順に**番号を配ります。
    # 排他は取るので**競合はありません**。番号の付き方だけが変わります。
    #
    # **幾何も位相も体積も同じ**なので、`soup` / `nary` / `cp5n` は素通りします。
    #
    # > **この変異は最初、素通りしました。** 決定性の検査が `ToMeshOptions` に
    # > プールを渡しておらず、**出口を 1 スレッドで回していた**からです。
    # > 「機構を足したら、その機構が空回りしていないことを別に検査する」の実例です
    # > （`IMPL-phase4.md` §2.4）。
    krisite_add_mutant(23 KRISITE_MUTATION_NO_EXIT_BARRIER=1 parallel FAIL)
    krisite_add_mutant(23 KRISITE_MUTATION_NO_EXIT_BARRIER=1 soup     PASS)
    krisite_add_mutant(23 KRISITE_MUTATION_NO_EXIT_BARRIER=1 nary     PASS)
    krisite_add_mutant(23 KRISITE_MUTATION_NO_EXIT_BARRIER=1 cp5n     PASS)

    # 変異 24: **同値な構成点の代表を正準に選ぶのをやめる**（`SPEC-phase4.md` §4.4）。
    #   **主検出器: parallel（並べ替え不変性）**
    #
    # 最初に来た点を代表にします。**幾何も位相も体積も変わりません。**
    # 同じ点の同次座標のスカラー倍が変わるだけです。
    #
    # **§7.1 の決定性でも捕まりません。** あれは同じコードどうしの比較なので、
    # 代表の選び方が非正準でも一致します。**多角形の並びを変える検査だけが捕まえます。**
    krisite_add_mutant(24 KRISITE_MUTATION_NO_CANONICAL_REPR=1 parallel FAIL)
    krisite_add_mutant(24 KRISITE_MUTATION_NO_CANONICAL_REPR=1 soup     PASS)
    krisite_add_mutant(24 KRISITE_MUTATION_NO_CANONICAL_REPR=1 nary     PASS)
    krisite_add_mutant(24 KRISITE_MUTATION_NO_CANONICAL_REPR=1 cp5n     PASS)

    # ---- CP5 の組合せでも検出されること（SPEC-phase2 §13）--------------------
    #
    # **CP5 は機構を全部有効にした状態です。** 機構が増えると検出器は移動し、
    # ときに**互いを覆い隠します**（実例: 案 D の T 解決が変異 3 / 4 を覆いました。
    # §9.3.2）。そこで「この組合せでも捕まるか」を別に固定します。
    #
    # **素通りの 2 件も固定します。** 変異 3 / 4 は T 解決に覆われて cp5 では
    # 素通りします（主検出器は split_strategy に移動済み。§2.4.6）。
    # **ここが FAIL に転じたら覆いが外れたということなので、報告してください。**
    #
    # 実測（`IMPL-phase2.md` §7.16）: 14 の変異のうち 12 が cp5 で落ち、
    # 素通りは変異 3 / 4 の 2 件だけでした。ここに登録したのはその代表です。
    krisite_add_mutant(3  KRISITE_MUTATION_LOCAL_PLANES=1      cp5 PASS)
    krisite_add_mutant(4  KRISITE_MUTATION_REACHING_TRIANGLE=1 cp5 PASS)
    krisite_add_mutant(6  KRISITE_MUTATION_LOOSE_EARLY_OUT=1   cp5 FAIL)
    krisite_add_mutant(7  KRISITE_MUTATION_CACHE_KEY_DROP=1    cp5 FAIL)
    krisite_add_mutant(9  KRISITE_MUTATION_SPLIT_BY_OWNER=1    cp5 FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1      cp5 FAIL)
    krisite_add_mutant(12 KRISITE_MUTATION_CLOSED_CELLS=1      cp5 FAIL)

    # 変異 14: §3.1 の辺平面の軸選択を非正準にする（最小ではなく最大の k を選ぶ）。
    #   **主検出器: interior（§3.1.3 の正準性の検査）**
    #
    # **幾何としては正しい平面が返ります。** 辺の直線を含み支持平面と一致しない、
    # という条件は満たすので、位相も体積も変わりません。**壊れるのは決定性だけ**で、
    # 位相・体積の検査はすべて素通りします。決定性を直接見ている検査が無ければ
    # 誰も捕まえません。
    krisite_add_mutant(14 KRISITE_MUTATION_EDGE_AXIS_LAST=1 interior        FAIL)
    krisite_add_mutant(14 KRISITE_MUTATION_EDGE_AXIS_LAST=1 cp3             PASS)
    krisite_add_mutant(14 KRISITE_MUTATION_EDGE_AXIS_LAST=1 adaptive        PASS)
    krisite_add_mutant(14 KRISITE_MUTATION_EDGE_AXIS_LAST=1 split_semantics PASS)
endif()
