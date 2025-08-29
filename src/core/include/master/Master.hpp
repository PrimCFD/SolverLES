#pragma once
#include <filesystem>
#include <memory>
#include <string>
#include "master/RunContext.hpp"
#include "master/Scheduler.hpp"
#include "master/FieldCatalog.hpp"
#include "master/PluginHost.hpp"
#include "master/io/IWriter.hpp"

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
 *    core::master::TimeControls tc;
 *    tc.dt = 1e-3;
 *    tc.t_end = 1.0;
 *    tc.write_every = 10;
 *    m.run(tc);
 * @endrst
 * 
 */


namespace core {
namespace mesh { class Layout; class HaloExchange; class Boundary; }

namespace master {

class Master {
public:
  Master(RunContext rc,
         const mesh::Layout& layout,
         mesh::HaloExchange& halos,
         mesh::Boundary& bcs);

  FieldCatalog& fields() noexcept { return fields_; }
  void set_writer(std::unique_ptr<io::IWriter> w) {
    writer_ = std::move(w);
    // Rebind scheduler to the new writer to avoid dangling references.
    if (writer_) {
      sched_ = std::make_unique<Scheduler>(rc_, layout_, halos_, bcs_, fields_, *writer_);
    }
  }

  // Plugins
  void load_plugin_library(const std::filesystem::path& lib) { plugins_.load_library(lib); }
  void configure_program(const std::string& key, const plugin::KV& cfg);

  void run(const TimeControls& tc);

private:
  RunContext rc_;
  const mesh::Layout& layout_;
  mesh::HaloExchange& halos_;
  mesh::Boundary& bcs_;

  FieldCatalog fields_;
  PluginHost   plugins_;
  std::unique_ptr<io::IWriter> writer_;
  std::unique_ptr<Scheduler>   sched_;
};

} // namespace master
} // namespace core
