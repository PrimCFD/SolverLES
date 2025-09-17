set(SIMPLE_BENCH_HDR ${CMAKE_SOURCE_DIR}/tests/performance/simple_bench.hpp)

# Binaries
add_executable(bench_allocate bench_allocate.cpp ${SIMPLE_BENCH_HDR})
add_executable(bench_index bench_index.cpp ${SIMPLE_BENCH_HDR})
add_executable(bench_copy bench_copy.cpp ${SIMPLE_BENCH_HDR})
add_executable(bench_overhead bench_overhead.cpp ${SIMPLE_BENCH_HDR})
add_executable(bench_io_write bench_io_write.cpp ${SIMPLE_BENCH_HDR})

target_include_directories(bench_allocate
                           PRIVATE ${CMAKE_SOURCE_DIR}/tests/performance)
target_include_directories(bench_index
                           PRIVATE ${CMAKE_SOURCE_DIR}/tests/performance)
target_include_directories(bench_copy
                           PRIVATE ${CMAKE_SOURCE_DIR}/tests/performance)
target_include_directories(bench_overhead
                           PRIVATE ${CMAKE_SOURCE_DIR}/tests/performance)
target_include_directories(bench_io_write
                           PRIVATE ${CMAKE_SOURCE_DIR}/tests/performance)

target_link_libraries(bench_allocate PRIVATE core)
target_link_libraries(bench_index PRIVATE core)
target_link_libraries(bench_copy PRIVATE core)
target_link_libraries(bench_overhead PRIVATE core)
target_link_libraries(bench_io_write PRIVATE core)

if(ENABLE_CUDA)
  target_link_libraries(bench_copy PRIVATE CUDA::cudart)
  target_compile_definitions(bench_copy PRIVATE HAVE_CUDA)
endif()

# Tests (serial, labeled "perf", CUDA copy auto-skips when CUDA is off)
add_perf_test(perf_alloc bench_allocate)
add_perf_test(perf_index bench_index)
add_perf_test(perf_overhead bench_overhead)
add_perf_test(perf_io_write bench_io_write)

# 16 MiB and 64 MiB roundtrips (H->D + D->H). If CUDA is off, these report
# "skipped".
add_perf_test(perf_copy_16MB bench_copy 16777216) # 16 MiB
add_perf_test(perf_copy_64MB bench_copy 67108864) # 64 MiB
