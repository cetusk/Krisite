# Krisite — ベンチマーク（SPEC-phase0.md §9「side() のスループットを計測」）

add_executable(arith_bench arith_bench.cpp)
target_link_libraries(arith_bench PRIVATE krisite)
target_include_directories(arith_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tests)

add_executable(pred_bench pred_bench.cpp)
target_link_libraries(pred_bench PRIVATE krisite)
target_include_directories(pred_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tests)
