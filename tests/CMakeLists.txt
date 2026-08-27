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
krisite_add_test(cp2           csg/test_cp2.cpp)
krisite_add_test(cp25          csg/test_cp25.cpp)
krisite_add_test(cp3           csg/test_cp3.cpp)
krisite_add_test(idempotence   csg/test_idempotence.cpp)
krisite_add_test(depth_monotonic csg/test_depth_monotonic.cpp)

# §12 の呼び出し比率。**この実行ファイルだけ**計数を有効にする
# （既定ビルドに計数のコストを持ち込まない）。
krisite_add_test(call_ratio    csg/test_call_ratio.cpp)
target_compile_definitions(call_ratio PRIVATE KRISITE_COUNT_PREDICATES=1)

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

    function(krisite_add_mutant tag define test expect)
        set(name "mut_${tag}_${test}")
        add_executable(${name} csg/test_${test}.cpp)
        target_link_libraries(${name} PRIVATE krisite krisite_test_util)
        target_compile_definitions(${name} PRIVATE ${define})
        add_test(NAME ${name} COMMAND ${name})
        set_tests_properties(${name} PROPERTIES TIMEOUT 3600 LABELS "mutation")
        if(expect STREQUAL "FAIL")
            set_tests_properties(${name} PROPERTIES WILL_FAIL TRUE)
        endif()
    endfunction()

    # 変異 1: §5.3 の第2段（値ベースの併合）を無効化
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
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp2             FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp25            FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 cp3             FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 depth_monotonic FAIL)
    krisite_add_mutant(2 KRISITE_MUTATION_OPEN_CELLS=1 idempotence     PASS)

    # 変異 3: §4.3.1 の分割平面を割り当て集合に変える（T 字接合が生じる）
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp1             FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp2            FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp25           FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 cp3             FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 depth_monotonic FAIL)
    krisite_add_mutant(3 KRISITE_MUTATION_LOCAL_PLANES=1 idempotence     PASS)

    # ---- Phase 2 ------------------------------------------------------------
    #
    # 変異 4: §2.3 の絞り込みの**基準**を「平面がセルを横切るか」から
    #         「三角形がセルに届くか」に戻す（SPEC-phase2 §9.3）。
    #
    # Phase 1 の変異 3 と同型ですが、あちらは大域集合を割り当て集合に**置き換える**
    # もので、こちらは絞り込みの**判定基準**を差し替えます。plane(T) は無限に延びる
    # ので、T が届かないセルにも切断点を生みます。
    #
    # **split_strategy が最も鋭い検出器です**（118 件）。絞り込み無効の側を
    # 正解器として直接突き合わせるためです。
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 split_strategy  FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp1             FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp2             FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp25            FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 cp3             FAIL)
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 depth_monotonic FAIL)
    # **冪等性は変異 4 も検出しません。** 比較の両側に同じ変異がかかるので差が出ない、
    # という Phase 1 と同じ構造的な盲点です。**素通りを固定して記録します。**
    krisite_add_mutant(4 KRISITE_MUTATION_REACHING_TRIANGLE=1 idempotence     PASS)

    # 変異 10: §2.4.3 の T 頂点解決を無効化する（SPEC-phase2 §9.3）。
    #
    # **T 解決は継ぎ目の正しさの単一障害点です。** 検出の経路は 2 つあり、
    # **測っているものが違います。**
    #
    #   tjunction      … 正の対照。**幾何で**検出する（効くべき場所で効かない）
    #   split_strategy … 負の対照の【非空虚の番人】。**計数で**検出する
    #                    （索引が候補を返しておらず、対照が空虚になっている）
    #
    # **固定深度のコーパスは素通りします。** 固定深度では T 解決がそもそも no-op なので、
    # 無効化しても出力が変わりません。**ケース 13（適応分割）が入るまでの暫定値です。**
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 tjunction       FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 split_strategy  FAIL)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 cp3             PASS)
    krisite_add_mutant(10 KRISITE_MUTATION_NO_TJUNCTION=1 depth_monotonic PASS)

    # 変異 11: T 頂点を 1 本の辺に 1 個しか入れない（一括でなく）。
    #
    # **1 本の辺に T 頂点が 2 個以上載る配置でしか検出できません**（§8 のケース 13 の要件）。
    # 正の対照はその配置を手で作っているので検出しますが、**固定深度のコーパスは
    # T 頂点が 0 個なので構造的に素通りします。**
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 tjunction      FAIL)
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 split_strategy PASS)
    krisite_add_mutant(11 KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME=1 cp3            PASS)
endif()
