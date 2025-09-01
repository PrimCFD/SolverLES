#include "master/io/Preflight.hpp"
#include <sstream>

namespace core::master::io
{

std::pair<bool, std::string> run_preflight(const WriterConfig& cfg, const WritePlan& plan,
                                           int mpi_world_size, std::size_t available_ram_bytes,
                                           std::size_t available_disk_bytes)
{
    std::size_t total_step_bytes = 0;
    for (auto& f : plan.fields)
        total_step_bytes += f.bytes;

    // Rough RAM bound: one staging buffer per field plus some headroom
    const std::size_t overhead = total_step_bytes / 8 + (64ull << 20); // 12.5% + 64MB
    const bool ram_ok = total_step_bytes + overhead < available_ram_bytes;

    // Disk check: reserve at least one step + metadata (5%)
    const bool disk_ok = (total_step_bytes * 1.05) < available_disk_bytes;

    std::ostringstream msg;
    msg << "Preflight: step_bytes=" << total_step_bytes << ", fields=" << plan.fields.size()
        << ", mpi_world_size=" << mpi_world_size << ", backend="
        << (cfg.backend == WriterConfig::Backend::CGNS
                ? "CGNS"
                : (cfg.backend == WriterConfig::Backend::XDMF ? "XDMF" : "Null"))
        << ".";

    bool ok = ram_ok && disk_ok;
    if (!ram_ok)
        msg << " RAM check failed.";
    if (!disk_ok)
        msg << " Disk check failed.";
    return {ok, msg.str()};
}

} // namespace core::master::io