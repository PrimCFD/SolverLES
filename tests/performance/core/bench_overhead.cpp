#include "master/Master.hpp"
#include "mesh/Mesh.hpp"
#include "simple_bench.hpp"

using namespace core;

static void run_noop_steps(int steps)
{
    master::RunContext rc{};
    // If you build with MPI and want to pass the communicator:
    // rc.mpi_comm = reinterpret_cast<void*>(MPI_COMM_WORLD);

    // Minimal mesh: interior 2x2x2, no ghosts for this overhead micro-bench
    mesh::Mesh mesh;
    mesh.local = {2, 2, 2};
    mesh.ng = 0;

    master::Master m(rc, mesh);

    // Use builtin "noop" program (installed by PluginHost)
    m.configure_program("noop", {});

    master::TimeControls tc;
    tc.dt = 1.0;
    tc.t_end = static_cast<double>(steps);
    tc.write_every = 1 << 30; // effectively disable writes (NullWriter by default)

    m.run(tc);
}

int main()
{
    // warm-up to avoid first-call noise
    run_noop_steps(32);

    constexpr int STEPS = 2000;
    auto [mean, stddev] = bench::run([&] { run_noop_steps(STEPS); });

    bench::report("orchestrator_noop_" + std::to_string(STEPS) + "steps", mean, stddev);
    return 0;
}
