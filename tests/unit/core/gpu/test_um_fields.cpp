#ifdef HAVE_CUDA
#include "master/Master.hpp"
#include <catch2/catch_test_macros.hpp>
#include <cuda_runtime.h>

__global__ void bump(double* p, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        p[i] += 1.0;
}

TEST_CASE("UM field survives orchestration and is device-touchable", "[cuda][um]")
{
    using namespace core;
    double* a = nullptr;
    int n = 64;
    REQUIRE(cudaMallocManaged(&a, n * sizeof(double)) == cudaSuccess);
    for (int i = 0; i < n; ++i)
        a[i] = 0.0;

    master::RunContext rc{};
    mesh::Layout layout;
    mesh::HaloExchange halos;
    mesh::Boundary bcs;
    master::Master m(rc, layout, halos, bcs);

    m.fields().register_scalar("a", a, sizeof(double), {4, 4, 4}, {1, 4, 16},
                               core::master::Stagger::Cell);
    m.fields().select_for_output("a");
    m.set_writer(std::make_unique<master::io::NullWriter>()); // you have this header uploaded

    m.configure_program("noop", {});
    master::TimeControls tc;
    tc.dt = 0.1;
    tc.t_end = 0.1;
    tc.write_every = 100; // 1 step, no writes
    m.run(tc);

    // now touch on device and verify host sees the change
    bump<<<(n + 127) / 128, 128>>>(a, n);
    REQUIRE(cudaDeviceSynchronize() == cudaSuccess);
    REQUIRE(a[0] == Approx(1.0));
    cudaFree(a);
}
#endif
