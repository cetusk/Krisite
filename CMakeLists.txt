# Krisite — exact, plane-based geometry for point clouds and meshes
#
# Phase 0: 固定幅厳密整数演算と平面ベース述語（docs/SPEC-phase0.md）
cmake_minimum_required(VERSION 3.20)
project(Krisite VERSION 0.0.1 LANGUAGES CXX
        DESCRIPTION "Exact, plane-based geometry for point clouds and meshes")

# ---- オプション（README.md「CMake オプション」の表）--------------------------

set(KRISITE_COORD_BITS 21 CACHE STRING "入力座標のビット幅 b")
option(KRISITE_BUILD_TESTS "テストのビルド" ON)
option(KRISITE_BUILD_TESTS_WITH_GMP "GMP 差分テスト（LGPL、テスト専用）" OFF)
option(KRISITE_BUILD_BENCH "ベンチマークのビルド" OFF)
option(KRISITE_CHECKED_ARITH "全演算のオーバーフロー検査（Debug では自動で ON）" OFF)
option(KRISITE_BUILD_MUTANTS "変異テスト（SPEC-phase1 §10.5）のビルド" OFF)
option(KRISITE_BUILD_TESTS_WITH_MANIFOLD "Manifold 正解器（Apache-2.0、テスト専用）" OFF)
# SPEC-phase2 §9.4: **適応分割モードを CI のマトリクスに加える。**
# `BoolOptions` の既定を「適応分割 + early-out + 構成点の保持」に切り替える。
# **固定深度モードは消えません**（§0.1 の正解器。オプションを明示すれば従来どおり）。
option(KRISITE_DEFAULT_ADAPTIVE "BoolOptions の既定を適応分割にする（SPEC-phase2 §9.4）" OFF)
# SPEC-phase4 §7.4: **並列モードを CI のマトリクスに加える。**
# `BoolOptions::threads` の既定をこの値にする。**0 なら既定は逐次のまま。**
# **逐次経路は消えません**（オプションを明示すれば従来どおり。§7.1 の正解器）。
set(KRISITE_DEFAULT_THREADS 0 CACHE STRING "BoolOptions::threads の既定（SPEC-phase4 §7.4）")
# **タグではなくコミットに固定します**（SPEC-phase1 §10.4）。
# タグは付け替えられるので、再現性の保証になりません。
# 下の SHA は v3.0.1 を CI で解決した結果です。更新するときは CI のログで
# 解決結果を確認してから書き換えてください。
set(KRISITE_MANIFOLD_COMMIT "98b8142519d35c13e0e25cfa9fd6e3a271403be6"
    CACHE STRING "Manifold の固定コミット（v3.0.1）")

# ---- ライブラリ本体（ヘッダオンリー）----------------------------------------

add_library(krisite INTERFACE)
add_library(Krisite::krisite ALIAS krisite)

target_include_directories(krisite INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(krisite INTERFACE cxx_std_20)
# SPEC-phase4 §2: `par::ThreadPool` が std::thread を使うので、利用側に伝播させる。
# **ヘッダオンリーのまま**で、リンクするのは pthread だけです。
find_package(Threads REQUIRED)
target_link_libraries(krisite INTERFACE Threads::Threads)
target_compile_definitions(krisite INTERFACE
    KRISITE_COORD_BITS=${KRISITE_COORD_BITS}
    # SPEC §3.4: デバッグビルドでは常に検査を有効にする
    $<$<OR:$<BOOL:${KRISITE_CHECKED_ARITH}>,$<CONFIG:Debug>>:KRISITE_CHECKED_ARITH=1>
    $<$<BOOL:${KRISITE_DEFAULT_ADAPTIVE}>:KRISITE_DEFAULT_ADAPTIVE=1>
    $<$<BOOL:${KRISITE_DEFAULT_THREADS}>:KRISITE_DEFAULT_THREADS=${KRISITE_DEFAULT_THREADS}>)

if(MSVC)
    # /utf-8: ソースも実行文字集合も UTF-8。日本語のコメントと文字列リテラルが
    #         ランナーのコードページに引きずられて化けるのを防ぐ（C4819 も消える）。
    # /Zc:__cplusplus: __cplusplus を実際の規格値にする。
    target_compile_options(krisite INTERFACE /W4 /permissive- /utf-8 /Zc:__cplusplus)
else()
    target_compile_options(krisite INTERFACE -Wall -Wextra)
endif()

# ---- GMP（テスト専用。ライブラリ本体には絶対にリンクしない）------------------
#
# ライセンス方針（THIRD_PARTY_LICENSES.md）/ SPEC §8.1:
#   GMP は LGPL。テストコードからのみ参照し、配布物には含めない。

set(KRISITE_GMP_FOUND OFF)
if(KRISITE_BUILD_TESTS AND KRISITE_BUILD_TESTS_WITH_GMP)
    find_path(GMP_INCLUDE_DIR gmp.h
              HINTS ENV GMP_ROOT PATH_SUFFIXES include include/x86_64-linux-gnu)
    find_library(GMP_LIBRARY NAMES gmp
                 HINTS ENV GMP_ROOT PATH_SUFFIXES lib lib/x86_64-linux-gnu)
    if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
        set(KRISITE_GMP_FOUND ON)
        add_library(krisite_gmp INTERFACE)
        target_include_directories(krisite_gmp INTERFACE ${GMP_INCLUDE_DIR})
        target_link_libraries(krisite_gmp INTERFACE ${GMP_LIBRARY})
        message(STATUS "Krisite: GMP 差分テストを有効化 (${GMP_LIBRARY})")
    else()
        message(FATAL_ERROR
            "KRISITE_BUILD_TESTS_WITH_GMP=ON ですが GMP が見つかりません。"
            " GMP_ROOT を設定してください。")
    endif()
endif()

# ---- Manifold（テスト専用の正解器。ライブラリ本体には絶対にリンクしない）-----
#
#   SPEC-phase1 §10.4 / THIRD_PARTY_LICENSES.md「Manifold」
#   Manifold は Apache-2.0。**テストの正解器としてのみ**使い、本体から呼ばない。
#
#   vendoring しない（自前コードとの境界を曖昧にしないため）。
#   `git submodule add` は settings.json の deny 対象なので FetchContent を使う。
#
#   推移的依存を増やさないため、機能はすべて切る。
#     MANIFOLD_CROSS_SECTION=OFF … Clipper2（BSL-1.0）を引かない
#     MANIFOLD_PAR=OFF           … TBB を引かない
#     MANIFOLD_DOWNLOADS=OFF     … ビルド中の自動取得を止める
set(KRISITE_MANIFOLD_FOUND OFF)
if(KRISITE_BUILD_TESTS AND KRISITE_BUILD_TESTS_WITH_MANIFOLD)
    include(FetchContent)
    set(MANIFOLD_TEST OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_PAR OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_PYBIND OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_CBIND OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_JSBIND OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_CROSS_SECTION OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_EXPORT OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_DOWNLOADS OFF CACHE BOOL "" FORCE)
    # GIT_SHALLOW は使えません。浅いクローンでは任意のコミットを取得できないためです。
    FetchContent_Declare(manifold
        GIT_REPOSITORY https://github.com/elalish/manifold.git
        GIT_TAG ${KRISITE_MANIFOLD_COMMIT})
    FetchContent_MakeAvailable(manifold)
    set(KRISITE_MANIFOLD_FOUND ON)
    # ライセンスの根拠を機械的に固定する。README ではなく LICENSE を読むこと
    # （根拠は THIRD_PARTY_LICENSES.md）。
    if(EXISTS "${manifold_SOURCE_DIR}/LICENSE")
        file(READ "${manifold_SOURCE_DIR}/LICENSE" _kri_manifold_license)
        string(FIND "${_kri_manifold_license}" "Apache License" _kri_apache)
        if(_kri_apache EQUAL -1)
            message(FATAL_ERROR
                "Manifold の LICENSE が Apache License ではありません。"
                " 本プロジェクトのライセンス方針に反するので停止します。")
        endif()
        message(STATUS "Krisite: Manifold の LICENSE を確認 (Apache License)")
    else()
        message(FATAL_ERROR "Manifold の LICENSE が見つかりません。")
    endif()
endif()

# ---- テスト -----------------------------------------------------------------

if(KRISITE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---- ベンチマーク -----------------------------------------------------------

if(KRISITE_BUILD_BENCH)
    add_subdirectory(bench)
endif()

message(STATUS "Krisite: b = ${KRISITE_COORD_BITS}, "
               "checked = ${KRISITE_CHECKED_ARITH} (Debug では常に ON), "
               "gmp = ${KRISITE_GMP_FOUND}, "
               "mutants = ${KRISITE_BUILD_MUTANTS}, "
               "adaptive_default = ${KRISITE_DEFAULT_ADAPTIVE}, "
               "threads_default = ${KRISITE_DEFAULT_THREADS}, "
               "manifold = ${KRISITE_MANIFOLD_FOUND}")
