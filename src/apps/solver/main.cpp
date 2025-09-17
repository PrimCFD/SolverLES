#include "master/FieldCatalog.hpp"
#include "master/Master.hpp"
#include "master/Views.hpp"
#include "master/io/AsyncWriter.hpp"
#include "master/io/CGNSWriter.hpp"
#include "master/io/ConfigYAML.hpp" // AppConfig + load_config_from_yaml()
#include "master/io/NullWriter.hpp"
#include "master/io/Preflight.hpp"
#include "master/io/WritePlan.hpp"
#include "master/io/WriterConfig.hpp"
#include "master/io/XdmfHdf5Writer.hpp"
#include "memory/MemoryManager.hpp"
#include "mesh/Mesh.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef __linux__
#include <sys/sysinfo.h>
#endif

using core::master::FieldCatalog;
using core::master::Master;
using core::master::RunContext;
using core::master::TimeControls;
using core::master::io::AsyncWriter;
using core::master::io::build_write_plan;
using core::master::io::CGNSWriter;
using core::master::io::IWriter;
using core::master::io::NullWriter;
using core::master::io::run_preflight;
using core::master::io::WritePlan;
using core::master::io::WriterConfig;
using core::master::io::XdmfHdf5Writer;
namespace fs = std::filesystem;

// Row-major byte-stride helper
static inline std::array<std::ptrdiff_t, 3> strides_bytes(int nx, int ny, int /*nz*/,
                                                          std::size_t elem)
{
    return {(std::ptrdiff_t) elem, (std::ptrdiff_t)(elem * nx), (std::ptrdiff_t)(elem * nx * ny)};
}

// Map AppConfig backend → WriterConfig backend
static inline WriterConfig::Backend to_writer_backend(AppConfig::Backend b)
{
    using B = AppConfig::Backend;
    using WB = WriterConfig::Backend;
    switch (b)
    {
    case B::Xdmf:
        return WB::XDMF;
    case B::Cgns:
        return WB::CGNS;
    default:
        return WB::Null;
    }
}

// Build plugin::KV from either unordered_map<string,string>
// or vector<pair<string,string>> (both are iterable pairs)
template <class Assoc> static core::master::plugin::KV make_kv(const Assoc& params)
{
    core::master::plugin::KV kv;
    for (const auto& p : params)
        kv.emplace(p.first, p.second);
    return kv;
}

int main(int argc, char** argv)
{
    // 1) Parse YAML config
    const std::string cfg_path = (argc > 1) ? argv[1] : "case.yaml";
    const AppConfig cfg = load_config_from_yaml(cfg_path);

    // 2) Build a RunContext (MPI/stream/memory)
    RunContext rc{};
#ifdef HAVE_MPI
    {
        int init = 0;
        MPI_Initialized(&init);
        if (!init)
        {
            int prov = 0;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);
        }
        rc.mpi_comm = reinterpret_cast<void*>(MPI_COMM_WORLD); // opaque in the ABI
    }
#endif
    rc.device_stream = nullptr; // set a real stream if/when you have a device
    rc.mem = &core::memory::MemoryManager::instance();

    // 3) Mesh from config
    core::mesh::Mesh mesh;
    mesh.local = cfg.local;
    mesh.ng = cfg.ng;
    mesh.periodic = cfg.periodic;

    // 4) App façade
    Master master(rc, mesh);

    // 5) Build WriterConfig
    WriterConfig wcfg{};
    wcfg.backend = to_writer_backend(cfg.io.backend);
    wcfg.path = cfg.io.path;
    // Only set precision override when requested (avoid relying on a “Native” enum)
    if (cfg.io.precision == AppConfig::Precision::F32)
    {
        wcfg.precision = WriterConfig::Precision::Float32;
    }
    else if (cfg.io.precision == AppConfig::Precision::F64)
    {
        wcfg.precision = WriterConfig::Precision::Float64;
    }
    // (If Precision::Native, leave wcfg.precision at its default)

    // Optional XDMF version (accept strings "v2"/"v3" or booly flag in your YAML)
    if (cfg.io.xdmf_version == "v2")
    {
        wcfg.xdmf_version = WriterConfig::XdmfVersion::V2;
    }
    else if (cfg.io.xdmf_version == "v3")
    {
        wcfg.xdmf_version = WriterConfig::XdmfVersion::V3;
    } // else keep default

    // 6) Construct concrete writer (+optional async) and attach to Master
    std::unique_ptr<IWriter> sink;
    switch (cfg.io.backend)
    {
    case AppConfig::Backend::Xdmf:
        sink = std::make_unique<XdmfHdf5Writer>(wcfg);
        break;
    case AppConfig::Backend::Cgns:
        sink = std::make_unique<CGNSWriter>(wcfg);
        break;
    default:
        sink = std::make_unique<NullWriter>();
        break;
    }
    if (cfg.io.async.enabled)
    {
        AsyncWriter::Options opt{};
        opt.max_queue = cfg.io.async.max_queue;
        opt.drop_to_null_on_overflow = cfg.io.async.drop_on_overflow;
        sink = std::make_unique<AsyncWriter>(std::move(sink), opt);
    }
    master.set_writer(std::move(sink));

    // 7) Load plugins & configure program (accepts map or vector<pair>)
    for (const auto& lib : cfg.plugin_libs)
        master.load_plugin_library(lib);
    master.configure_program(cfg.program_key, make_kv(cfg.program_params));

    // 8) Allocate and register fields owned by the app
    auto& mm = core::memory::MemoryManager::instance();
    const auto ext = mesh.extents();
    const int nx_tot = ext[0], ny_tot = ext[1], nz_tot = ext[2];

    struct OwnedField
    {
        std::string name;
        double* ptr{};
    };
    std::vector<OwnedField> owned;
    owned.reserve(std::max<std::size_t>(1, cfg.fields_output.size()));
    if (cfg.fields_output.empty())
    {
        owned.push_back({"rho", mm.allocate<double>((std::size_t) nx_tot * ny_tot * nz_tot)});
    }
    else
    {
        for (const auto& name : cfg.fields_output)
            owned.push_back({name, mm.allocate<double>((std::size_t) nx_tot * ny_tot * nz_tot)});
    }

    FieldCatalog& fc = master.fields();
    const std::size_t elem = sizeof(double);
    for (auto& f : owned)
    {
        std::fill_n(f.ptr, (std::size_t) nx_tot * ny_tot * nz_tot, 0.0);
        fc.register_scalar(f.name.c_str(), f.ptr, elem, {nx_tot, ny_tot, nz_tot},
                           strides_bytes(nx_tot, ny_tot, nz_tot, elem));
        fc.select_for_output(f.name.c_str());
    }

    // 9) Optional preflight (RAM/disk estimate)
    if (cfg.io.preflight.enabled)
    {
        auto views = fc.selected_for_output();

        // Derive byte override directly from AppConfig precision
        std::size_t override_b = 0;
        if (cfg.io.precision == AppConfig::Precision::F32)
            override_b = 4;
        else if (cfg.io.precision == AppConfig::Precision::F64)
            override_b = 8;

        WritePlan plan = build_write_plan(
            std::span<const core::master::AnyFieldView>(views.data(), views.size()), override_b);

        std::size_t avail_disk = 0, avail_ram = 0;
        try
        {
            avail_disk = fs::space(wcfg.path).available;
        }
        catch (...)
        {
            avail_disk = 0;
        }
#ifdef __linux__
        {
            struct sysinfo si
            {
            };
            if (!sysinfo(&si))
                avail_ram = std::size_t(si.freeram) * std::size_t(si.mem_unit);
        }
#endif
#ifdef HAVE_MPI
        int world = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &world);
#else
        int world = 1;
#endif
        auto [ok, msg] = run_preflight(wcfg, plan, world, avail_ram, avail_disk);
        std::cerr << msg << "\n";
        if (!ok)
            std::cerr << "Preflight warnings: consider reducing outputs or precision.\n";
    }

    // 10) Time controls & run
    TimeControls tc{};
    tc.dt = cfg.dt;
    tc.t_end = cfg.t_end;
    tc.write_every = cfg.write_every_steps;
    tc.case_name = cfg.case_name;
    master.run(tc);

    // 11) Cleanup
    for (auto& f : owned)
        mm.release(f.ptr);
#ifdef HAVE_MPI
    {
        int fin = 0;
        MPI_Finalized(&fin);
        if (!fin)
            MPI_Finalize();
    }
#endif
    return 0;
}
