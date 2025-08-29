#include "master/Master.hpp"
#include "master/io/NullWriter.hpp"

namespace core {
namespace master {

Master::Master(RunContext rc,
               const mesh::Layout& layout,
               mesh::HaloExchange& halos,
               mesh::Boundary& bcs)
  : rc_(rc)
  , layout_(layout)
  , halos_(halos)
  , bcs_(bcs)
{
  if (!writer_) writer_ = std::make_unique<io::NullWriter>();
  sched_ = std::make_unique<Scheduler>(rc_, layout_, halos_, bcs_, fields_, *writer_);
}

void Master::configure_program(const std::string& key, const plugin::KV& cfg) {
  auto prog = plugins_.make_program(key, cfg, rc_);
  sched_->set_program(std::move(prog));
}

void Master::run(const TimeControls& tc) {
  sched_->run(tc);
}

} // namespace master
} // namespace core
