#pragma once
#include "mesh/Mesh.hpp"
#include <cstddef>
#include <string>
#include <vector>

/**
 * @file WriterConfig.hpp
 * @brief Configuration for writers (backend selection, precision, policies).
 *
 * @details
 * `WriterConfig` selects the concrete backend (XDMF/CGNS/Null) and optional policies:
 * precision down-cast (`Float64`→`Float32`), series layout, compression, preflight checks,
 * and MPI/monitoring stubs.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   WriterConfig cfg;
 *   cfg.backend   = WriterConfig::Backend::XDMF;
 *   cfg.precision = WriterConfig::Precision::Float32;
 *   cfg.path      = "out/caseA";
 * @endrst
 */

namespace core::master::io
{

struct WriterConfig
{
    enum class Backend
    {
        CGNS,
        XDMF,
        Null
    };
    enum class Series
    {
        Single,
        PerStep
    };

    Backend backend = Backend::CGNS;
    Series series = Series::Single;

    std::string path = "output";     // case directory
    std::vector<std::string> fields; // names to write (empty = all selected)

    bool include_halos = false; // write rind/halo cells
    bool async = true;          // decorate with AsyncWriter

    // Casting on write (helps with many-field memory; eg double→float)
    enum class Precision
    {
        Float32,
        Float64
    };
    Precision precision = Precision::Float32;

    // Compression only where it doesn’t hurt parallel IO (XDMF side)
    bool compress = false;

    // XDMF XML version selector (affects only the XML index; HDF5 stays the same)
    enum class XdmfVersion
    {
        V2,
        V3
    };
    XdmfVersion xdmf_version = XdmfVersion::V2; // default

    struct MPI
    {
        enum class Mode
        {
            Auto,
            ParallelHDF5,
            Aggregators
        };
        Mode mode = Mode::Auto;
        int aggregators_per_node = 1;
    } mpi;

    struct Monitoring
    {
        int stride = 0; // 0 = off; every N steps produce small snapshot
        int decimate[3] = {1, 1, 1};
    } monitoring;

    enum class Preflight
    {
        Off,
        Warn,
        Strict
    };
    Preflight preflight = Preflight::Strict;

    // ---- Parallel IO context (POD only; no MPI includes in public headers) ----
    // If 'mesh' is non-null and 'mpi_cart_comm' is non-null, writers may use parallel IO.
    // Mesh provides: local, global, global_lo (rank's offset), and ng.
    const core::mesh::Mesh* mesh = nullptr; // decomposition & ghost width
    void* mpi_cart_comm = nullptr;          // boxed MPI_Comm (nullptr => serial)
};

} // namespace core::master::io