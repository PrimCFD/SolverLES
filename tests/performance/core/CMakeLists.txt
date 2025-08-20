# Binaries
add_executable(bench_allocate bench_allocate.cpp simple_bench.hpp)
add_executable(bench_index bench_index.cpp simple_bench.hpp)
add_executable(bench_copy bench_copy.cpp simple_bench.hpp)

target_link_libraries(bench_allocate PRIVATE core)
target_link_libraries(bench_index PRIVATE core)
target_link_libraries(bench_copy PRIVATE core)

if(ENABLE_CUDA)
  target_link_libraries(bench_copy PRIVATE CUDA::cudart)
  target_compile_definitions(bench_copy PRIVATE HAVE_CUDA)
endif()

# Tests (serial, labeled "perf", CUDA copy auto-skips when CUDA is off)
add_perf_test(perf_alloc bench_allocate)
add_perf_test(perf_index bench_index)

# 16 MiB and 64 MiB roundtrips (H->D + D->H). If CUDA is off, these report
# "skipped".
add_perf_test(perf_copy_16MB bench_copy 16777216) # 16 MiB
add_perf_test(perf_copy_64MB bench_copy 67108864) # 64 MiB
