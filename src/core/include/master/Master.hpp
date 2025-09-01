#pragma once
#include "master/FieldCatalog.hpp"
#include "master/PluginHost.hpp"
#include "master/RunContext.hpp"
#include "master/Scheduler.hpp"
#include "master/io/IWriter.hpp"
#include "mesh/Mesh.hpp"
#include <filesystem>
#include <memory>
#include <string>

/**
 * @file Master.hpp
 * @brief Application façade that wires context, program, writer, and scheduler.
 *
 * @details
 * Typical usage:
 *
 * @rst
 * .. code-block:: cpp
 *
 *   core::master::RunContext rc{};
 *   core::mesh::Mesh mesh{ .local={nx,ny,nz}, .ng=ng };
 *   core::master::Master M(rc, mesh);
 *
 *   M.fields().register_scalar("rho", rho_ptr, sizeof(double), {nx,ny,nz},
 *                              {8, nx*8, nx*ny*8});
 *   M.fields().select_for_output("rho");
 *
 *   auto W = std::make_unique<core::master::io::XdmfHdf5Writer>(cfg);
 *   M.set_writer(std::move(W));
 *
 *   M.load_plugin_library("libphysics.so");
 *   M.configure_program("noop", {}); // builtin
 *
 *   core::master::TimeControls tc{.dt=1e-3, .t_end=1.0, .write_every=10};
 *   M.run(tc);
 *
 * @endrst
 *
 */

namespace core::master
{
class Master
{
  public:
    Master(RunContext rc, const core::mesh::Mesh& mesh);

    FieldCatalog& fields() noexcept { return fields_; }

    void set_writer(std::unique_ptr<io::IWriter> w)
    {
        writer_ = std::move(w);
        if (writer_) // rebind to avoid dangling refs
            sched_ = std::make_unique<Scheduler>(rc_, fields_, *writer_, mesh_);
    }

    // Plugins
    void load_plugin_library(const std::filesystem::path& lib) { plugins_.load_library(lib); }
    void configure_program(const std::string& key, const plugin::KV& cfg);

    void run(const TimeControls& tc);

  private:
    RunContext rc_;
    core::mesh::Mesh mesh_;
    FieldCatalog fields_;
    PluginHost plugins_;
    std::unique_ptr<io::IWriter> writer_;
    std::unique_ptr<Scheduler> sched_;
};

} // namespace core::master
