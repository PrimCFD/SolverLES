#include "master/Scheduler.hpp"
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"
#include "master/io/IWriter.hpp"
#include <cmath>


// Keep MPI optional: don't include mesh headers that pull <mpi.h>.
// Declare only what we need to call.
namespace core { namespace mesh {
  class HaloExchange { public: void start(); void finish(); };
  class Boundary     { public: void apply(); };
  // class Layout; // already forward-declared in the header; no need to redeclare
}} // namespace core::mesh

namespace core {
namespace master {

using plugin::Phase;
using plugin::has;

Scheduler::Scheduler(const RunContext& rc,
                     const mesh::Layout& layout,
                     mesh::HaloExchange& halos,
                     mesh::Boundary& bcs,
                     FieldCatalog& fields,
                     io::IWriter& writer)
  : rc_(rc), layout_(layout), halos_(halos), bcs_(bcs), fields_(fields), writer_(writer)
{}

void Scheduler::set_program(std::unique_ptr<plugin::IProgram> program) {
  program_ = std::move(program);
}

void Scheduler::run(const TimeControls& tc) {
  if (!program_) return; // nothing to do
  writer_.open_case("output", layout_);

  const int steps = (tc.dt > 0.0 && tc.t_end > 0.0)
      ? static_cast<int>(std::llround(tc.t_end / tc.dt))
      : 0;
  double t = 0.0;

  for (int s = 0; s < steps; ++s, t += tc.dt) {
    auto plan = program_->plan_step(tc.dt);

    // One coarse tile covering the whole domain for now; tiling can be added later.
    MeshTileView tile;
    tile.box = {{0,0,0}, {0,0,0}}; // unknown extents here; scheduler remains a stub
    tile.stream = rc_.device_stream;

    // 1) PreExchange tiled actions (rare)
    for (auto& a : plan.tiled)
      if (has(a->info().phases, Phase::PreExchange)) a->execute(tile, fields_, tc.dt);

    // 2) Halo exchanges
    halos_.start();

    // 3) Interior work during exchange
    for (auto& a : plan.tiled)
      if (has(a->info().phases, Phase::Interior)) a->execute(tile, fields_, tc.dt);

    halos_.finish();

    // 4) After exchange, before BCs
    for (auto& a : plan.tiled)
      if (has(a->info().phases, Phase::PostExchange)) a->execute(tile, fields_, tc.dt);
    for (auto& g : plan.globals)
      if (has(g->phases(), Phase::PostExchange)) g->run(fields_, tc.dt);

    // 5) Apply BCs
    bcs_.apply();

    // 6) PostBC work
    for (auto& a : plan.tiled)
      if (has(a->info().phases, Phase::PostBC)) a->execute(tile, fields_, tc.dt);
    for (auto& g : plan.globals)
      if (has(g->phases(), Phase::PostBC)) g->run(fields_, tc.dt);

    // 7) EndStep hooks
    for (auto& a : plan.tiled)
      if (has(a->info().phases, Phase::EndStep)) a->execute(tile, fields_, tc.dt);
    for (auto& g : plan.globals)
      if (has(g->phases(), Phase::EndStep)) g->run(fields_, tc.dt);

    // I/O cadence
    if (tc.write_every > 0 && (s % tc.write_every == 0)) {
      io::WriteRequest req;
      req.step = s;
      req.time = t;
      for (auto& v : fields_.selected_for_output()) req.fields.push_back(v);
      writer_.write(req);
    }
  }

  writer_.close();
}

} // namespace master
} // namespace core
