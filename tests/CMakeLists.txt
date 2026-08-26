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
krisite_add_test(topology      mesh/test_topology.cpp)
krisite_add_test(uniform_grid  octree/test_uniform_grid.cpp)
krisite_add_test(faces         csg/test_faces.cpp)
krisite_add_test(fragment      csg/test_fragment.cpp)
krisite_add_test(cp1           csg/test_cp1.cpp)
krisite_add_test(cp2           csg/test_cp2.cpp)
krisite_add_test(depth_monotonic csg/test_depth_monotonic.cpp)

# ---- GMP 差分テスト（LGPL、テスト専用）--------------------------------------

if(KRISITE_GMP_FOUND)
    # SPEC §9: 算術 10^7 件、述語 10^6 件以上
    krisite_add_test(fixed_int_gmp  arith/test_fixed_int_gmp.cpp 10000000)
    target_link_libraries(fixed_int_gmp PRIVATE krisite_gmp)

    krisite_add_test(predicates_gmp geom/test_predicates_gmp.cpp 1000000)
    target_link_libraries(predicates_gmp PRIVATE krisite_gmp)

    krisite_add_test(plane_order_gmp geom/test_plane_order_gmp.cpp 300000)
    target_link_libraries(plane_order_gmp PRIVATE krisite_gmp)

    # 開発中に素早く回すための短縮版
    add_test(NAME fixed_int_gmp_quick COMMAND fixed_int_gmp 40000)
    add_test(NAME predicates_gmp_quick COMMAND predicates_gmp 5000)
    add_test(NAME plane_order_gmp_quick COMMAND plane_order_gmp 3000)
endif()
