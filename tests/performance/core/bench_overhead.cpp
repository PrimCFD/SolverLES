#include "mesh_seams.hpp"     // minimal no-op Layout/HaloExchange/Boundary
#include "master/Master.hpp"  
#include "simple_bench.hpp"

using namespace core;

static void run_noop_steps(int steps)
{
  master::RunContext rc{};
  mesh::Layout       layout;
  mesh::HaloExchange halos;
  mesh::Boundary     bcs;

  master::Master m(rc, layout, halos, bcs);

  // Use builtin "noop" program (installed by PluginHost)
  m.configure_program("noop", {});

  master::TimeControls tc;
  tc.dt          = 1.0;
  tc.t_end       = static_cast<double>(steps);
  tc.write_every = 1 << 30;   // effectively disable writes

  m.run(tc);
}

int main()
{
  // warm-up to avoid first-call noise
  run_noop_steps(32);

  constexpr int STEPS = 2000;
  auto [mean, stddev] = bench::run([&]{
    run_noop_steps(STEPS);
  });

  bench::report("orchestrator_noop_" + std::to_string(STEPS) + "steps",
                mean, stddev);
  return 0;
}
