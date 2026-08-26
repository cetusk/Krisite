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

krisite_add_test(fixed_int arith/test_fixed_int.cpp)
