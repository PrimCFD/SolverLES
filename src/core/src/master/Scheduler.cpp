#include "master/Scheduler.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/RunContext.hpp"
#include "master/io/IWriter.hpp"
#include "master/plugin/Program.hpp"

#include "master/Views.hpp" // AnyFieldView, MeshTileView, make_interior_copy
#include "mesh/Field.hpp"   // typed Field<T> view (no MPI)
#include "mesh/Mesh.hpp"    // mesh::Mesh for extents/ghosts

#include <cmath>

namespace core::master
{
using plugin::has;
using plugin::Phase;

Scheduler::Scheduler(const RunContext& rc, FieldCatalog& fields, io::IWriter& writer,
                     const core::mesh::Mesh& mesh)
    : rc_(rc), fields_(fields), writer_(writer), mesh_(mesh)
{
}

void Scheduler::set_program(std::unique_ptr<plugin::IProgram> program)
{
    program_ = std::move(program);
}

void Scheduler::run(const TimeControls& tc)
{
    if (!program_)
        return; // nothing to do

    // Writer API
    writer_.open_case(tc.case_name.c_str());

    const int steps =
        (tc.dt > 0.0 && tc.t_end > 0.0) ? static_cast<int>(std::llround(tc.t_end / tc.dt)) : 0;
    double t = 0.0;

    // One coarse tile covering the whole interior domain for now.
    MeshTileView tile;
    tile.box.lo = {0, 0, 0};
    tile.box.hi = {mesh_.local[0], mesh_.local[1], mesh_.local[2]};
    tile.stream = rc_.device_stream;
    tile.mesh = &mesh_;

    for (int s = 0; s < steps; ++s, t += tc.dt)
    {
        auto plan = program_->plan_step(tc.dt);

        // PreExchange (rare)
        for (auto& a : plan.tiled)
            if (has(a->info().phases, Phase::PreExchange))
                a->execute(tile, fields_, tc.dt);

        // Halo exchange (MPI only; no-op otherwise)
        exchange_all_fields(fields_, mesh_, rc_.mpi_comm);

        // Interior phase
        for (auto& a : plan.tiled)
            if (has(a->info().phases, Phase::Interior))
                a->execute(tile, fields_, tc.dt);

        // PostExchange
        for (auto& a : plan.tiled)
            if (has(a->info().phases, Phase::PostExchange))
                a->execute(tile, fields_, tc.dt);
        for (auto& g : plan.globals)
            if (has(g->phases(), Phase::PostExchange))
                g->run(fields_, tc.dt);

        // PostBC
        for (auto& a : plan.tiled)
            if (has(a->info().phases, Phase::PostBC))
                a->execute(tile, fields_, tc.dt);
        for (auto& g : plan.globals)
            if (has(g->phases(), Phase::PostBC))
                g->run(fields_, tc.dt);

        // EndStep
        for (auto& a : plan.tiled)
            if (has(a->info().phases, Phase::EndStep))
                a->execute(tile, fields_, tc.dt);
        for (auto& g : plan.globals)
            if (has(g->phases(), Phase::EndStep))
                g->run(fields_, tc.dt);

        // I/O cadence
        const int out_step = s + 1;
        const double out_time = (s + 1) * tc.dt; // state *after* this step
        if (tc.write_every > 0 && (out_step % tc.write_every == 0))
        {
            io::WriteRequest req;
            req.step = out_step;
            req.time = out_time;
            for (const auto& v : fields_.selected_for_output())
                req.fields.push_back(make_interior_copy(v, mesh_.ng));
            writer_.write(req);
        }
    }

    writer_.close();
}

} // namespace core::master
