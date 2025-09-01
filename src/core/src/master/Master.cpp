#include "master/Master.hpp"
#include "master/io/NullWriter.hpp"

namespace core::master
{

Master::Master(RunContext rc, const core::mesh::Mesh& mesh) : rc_(rc), mesh_(mesh)
{
    if (!writer_)
        writer_ = std::make_unique<io::NullWriter>();
    sched_ = std::make_unique<Scheduler>(rc_, fields_, *writer_, mesh_);
}

void Master::configure_program(const std::string& key, const plugin::KV& cfg)
{
    auto prog = plugins_.make_program(key, cfg, rc_);
    sched_->set_program(std::move(prog));
}

void Master::run(const TimeControls& tc)
{
    sched_->run(tc);
}

} // namespace core::master
