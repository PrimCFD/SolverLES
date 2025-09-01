#include "master/Scheduler.hpp"
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/io/IWriter.hpp"
#include "master/plugin/Program.hpp"

#include "mesh/Field.hpp" // typed Field<T> view (no MPI)
#include "mesh/Mesh.hpp"  // mesh::Mesh for extents/ghosts

#include <cmath>

// Only include the MPI-based halo exchange when MPI is enabled.
#ifdef HAVE_MPI
#include "mesh/HaloExchange.hpp"
#endif

namespace core::master
{
using plugin::has;
using plugin::Phase;

// ------------------------------------------------------------------
// Halo exchange helper:
//  - When HAVE_MPI: exchange ghosts for all registered fields.
//  - Otherwise: no-op (builds without MPI).
// ------------------------------------------------------------------
namespace
{
#ifdef HAVE_MPI
    static void exchange_all_fields(FieldCatalog& cat, const core::mesh::Mesh& m,
                                    void* mpi_comm_void)
    {
        // RunContext.mpi_comm is stored opaquely as void*; cast only in MPI builds.
        MPI_Comm comm = static_cast<MPI_Comm>(mpi_comm_void);

        for (const AnyFieldView& v : cat.all_views())
        {
            const auto e = v.extents;
            if (v.elem_size == sizeof(double))
            {
                core::mesh::Field<double> f(static_cast<double*>(v.host_ptr), e, m.ng);
                core::mesh::exchange_ghosts(f, m, comm);
            }
            else if (v.elem_size == sizeof(float))
            {
                core::mesh::Field<float> f(static_cast<float*>(v.host_ptr), e, m.ng);
                core::mesh::exchange_ghosts(f, m, comm);
            }
            else
            {
                // Extend with more primitive types if needed (e.g., int32).
            }
        }
    }
#else
    static void exchange_all_fields(FieldCatalog&, const core::mesh::Mesh&, void*)
    {
        // No MPI: nothing to do.
    }
#endif
} // namespace

// ------------------------------------------------------------------
// Scheduler
// ------------------------------------------------------------------

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

    // Writer API without deprecated Layout parameter.
    writer_.open_case("output");

    const int steps =
        (tc.dt > 0.0 && tc.t_end > 0.0) ? static_cast<int>(std::llround(tc.t_end / tc.dt)) : 0;
    double t = 0.0;

    // One coarse tile covering the whole interior domain for now.
    MeshTileView tile;
    tile.box.lo = {0, 0, 0};
    tile.box.hi = {mesh_.local[0], mesh_.local[1], mesh_.local[2]};
    tile.stream = rc_.device_stream;

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

        // (Optional) BCs hook: once you provide a config, call mesh::apply_* here.

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
        if (tc.write_every > 0 && (s % tc.write_every == 0))
        {
            io::WriteRequest req;
            req.step = s;
            req.time = t;
            for (const auto& v : fields_.selected_for_output())
                req.fields.push_back(v);
            writer_.write(req);
        }
    }

    writer_.close();
}

} // namespace core::master
