#pragma once
#include <string>
#include <utility>

#include "WritePlan.hpp"
#include "WriterConfig.hpp"

/**
 * @file Preflight.hpp
 * @brief RAM/disk sanity checks for I/O configurations.
 *
 * @details
 * Estimates total bytes per step from a :cpp:struct:`WritePlan` and verifies simple RAM/disk
 * headroom constraints before the run. The caller can enforce `Warn`/`Strict` policies.
 *
 * @return `{ok, message}` where `message` summarizes byte counts and the backend choice.
 */

namespace core::master::io
{

// Returns {ok, message}. If ok=false and strict preflight is enabled, the app should abort.
std::pair<bool, std::string> run_preflight(const WriterConfig& cfg, const WritePlan& plan,
                                           int mpi_world_size, std::size_t available_ram_bytes,
                                           std::size_t available_disk_bytes);

} // namespace core::master::io