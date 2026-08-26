# Krisite — exact, plane-based geometry for point clouds and meshes
#
# Phase 0: 固定幅厳密整数演算と平面ベース述語（docs/SPEC-phase0.md）
cmake_minimum_required(VERSION 3.20)
project(Krisite VERSION 0.0.1 LANGUAGES CXX
        DESCRIPTION "Exact, plane-based geometry for point clouds and meshes")

# ---- オプション（CLAUDE.md「CMake オプション」の表）--------------------------

set(KRISITE_COORD_BITS 21 CACHE STRING "入力座標のビット幅 b")
option(KRISITE_BUILD_TESTS "テストのビルド" ON)
option(KRISITE_BUILD_TESTS_WITH_GMP "GMP 差分テスト（LGPL、テスト専用）" OFF)
option(KRISITE_BUILD_BENCH "ベンチマークのビルド" OFF)
option(KRISITE_CHECKED_ARITH "全演算のオーバーフロー検査（Debug では自動で ON）" OFF)

# ---- ライブラリ本体（ヘッダオンリー）----------------------------------------

add_library(krisite INTERFACE)
add_library(Krisite::krisite ALIAS krisite)

target_include_directories(krisite INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(krisite INTERFACE cxx_std_20)
target_compile_definitions(krisite INTERFACE
    KRISITE_COORD_BITS=${KRISITE_COORD_BITS}
    # SPEC §3.4: デバッグビルドでは常に検査を有効にする
    $<$<OR:$<BOOL:${KRISITE_CHECKED_ARITH}>,$<CONFIG:Debug>>:KRISITE_CHECKED_ARITH=1>)

if(MSVC)
    target_compile_options(krisite INTERFACE /W4 /permissive-)
else()
    target_compile_options(krisite INTERFACE -Wall -Wextra)
endif()

# ---- テスト -----------------------------------------------------------------

if(KRISITE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "Krisite: b = ${KRISITE_COORD_BITS}, "
               "checked = ${KRISITE_CHECKED_ARITH} (Debug では常に ON)")
