# Krisite — ベンチマーク（SPEC-phase0.md §9「side() のスループットを計測」）

add_executable(arith_bench arith_bench.cpp)
target_link_libraries(arith_bench PRIVATE krisite)
target_include_directories(arith_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tests)

add_executable(pred_bench pred_bench.cpp)
target_link_libraries(pred_bench PRIVATE krisite)
target_include_directories(pred_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tests)

# Krisite — スープ経路の計数と時間（SPEC-phase3.md §11）。
# **コーパスを使うので tests/ を include します**（ベンチも同じ配置で測るため）。
add_executable(soup_bench soup_bench.cpp)
target_link_libraries(soup_bench PRIVATE krisite)
target_include_directories(soup_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tests)
