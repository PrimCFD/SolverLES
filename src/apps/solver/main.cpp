#include "master/FieldCatalog.hpp"
#include "master/Log.hpp"
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
#include <iomanip> // setw, left/right formatting
#include <iostream>
#include <memory>
#include <optional>
#include <sstream> // ostringstream
#include <string>
#include <utility>
#include <vector>

#include <cstdlib>
#include <fstream> // for std::ifstream (CPU topology probe)
#include <mpi.h>
#include <omp.h>
#include <set> // if you use std::set in the topology code
#include <thread>

#include <petscsys.h> // PETSc Initialization

#ifdef __linux__
#include <sched.h>
#include <sys/sysinfo.h>
#endif

#include "memory/MpiBox.hpp"

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

// ---- Robust MPI+PETSc once-only lifetime -----------------------------------
// Initializes MPI (with FUNNELED) and PETSc exactly once; finalizes in reverse order
// only if we were the owner who initialized them. Safe under mpirun or standalone.
struct MpiPetscOnce
{
    bool mpi_owner{false};
    bool petsc_owner{false};

    MpiPetscOnce(int& argc, char**& argv)
    {
        int inited = 0;
        MPI_Initialized(&inited);
        if (!inited)
        {
            int provided = MPI_THREAD_SINGLE;
            MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
            mpi_owner = true;
        }
        PetscBool pinited = PETSC_FALSE;
        PetscInitialized(&pinited);
        if (!pinited)
        {
            if (PetscInitialize(&argc, &argv, nullptr, nullptr))
            {
                LOGE("ERROR: PetscInitialize failed\n");
                std::abort();
            }
            petsc_owner = true;
        }
    }
    ~MpiPetscOnce()
    {
        // Finalize PETSc first
        PetscBool pfin = PETSC_FALSE;
        PetscFinalized(&pfin);
        if (!pfin)
        {
            PetscBool pinited = PETSC_FALSE;
            PetscInitialized(&pinited);
            if (pinited && petsc_owner)
            {
                PetscFinalize();
            }
        }
        // Then MPI
        int mfin = 0;
        MPI_Finalized(&mfin);
        if (!mfin)
        {
            int minit = 0;
            MPI_Initialized(&minit);
            if (minit && mpi_owner)
            {
                MPI_Finalize();
            }
        }
    }
};

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

// Build plugin::KV from either unordered_map<string,string> or vector<pair<string,string>>
template <class Assoc> static core::master::plugin::KV make_kv(const Assoc& params)
{
    core::master::plugin::KV kv;
    for (const auto& p : params)
        kv.emplace(p.first, p.second);
    return kv;
}

// -------- OpenMP/CPU setup & warm-up helpers --------------------------------
static int count_available_logical_cpus()
{
#ifdef __linux__
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
        return CPU_COUNT(&mask);
#endif
    unsigned n = std::thread::hardware_concurrency();
    return n ? static_cast<int>(n) : 1;
}

// Best-effort physical core count on Linux; falls back to logical
static int detect_physical_cores_linux()
{
#ifdef __linux__
    try
    {
        namespace fs = std::filesystem;
        std::set<std::pair<int, int>> pkg_core;
        for (const auto& e : fs::directory_iterator("/sys/devices/system/cpu"))
        {
            if (!e.is_directory())
                continue;
            const std::string name = e.path().filename().string();
            if (name.rfind("cpu", 0) != 0)
                continue;
            const std::string topo = (e.path() / "topology").string();
            std::ifstream f1(topo + "/physical_package_id"), f2(topo + "/core_id");
            if (!f1 || !f2)
                continue;
            int pkg = -1, core = -1;
            f1 >> pkg;
            f2 >> core;
            if (pkg >= 0 && core >= 0)
                pkg_core.emplace(pkg, core);
        }
        if (!pkg_core.empty())
            return static_cast<int>(pkg_core.size());
    }
    catch (...)
    {
    }
#endif
    return count_available_logical_cpus();
}

// --- MAC invariants -----------------------------------------

static inline void assert_mac_extents(const core::master::AnyFieldView& v,
                                      const core::mesh::Mesh& m)
{
    auto nx = m.local[0], ny = m.local[1], nz = m.local[2], ng = m.ng;
    auto want = std::array<int, 3>{nx + 2 * ng, ny + 2 * ng, nz + 2 * ng};
    if (v.stagger == core::master::Stagger::IFace)
        want = {nx + 1 + 2 * ng, ny + 2 * ng, nz + 2 * ng};
    if (v.stagger == core::master::Stagger::JFace)
        want = {nx + 2 * ng, ny + 1 + 2 * ng, nz + 2 * ng};
    if (v.stagger == core::master::Stagger::KFace)
        want = {nx + 2 * ng, ny + 2 * ng, nz + 1 + 2 * ng};
    if (v.extents[0] != want[0] || v.extents[1] != want[1] || v.extents[2] != want[2])
    {
        LOGE("[MAC] Extents for '%s' do not match MAC totals. got=%dx%dx%d want=%dx%dx%d\n",
             v.name.c_str(), v.extents[0], v.extents[1], v.extents[2], want[0], want[1], want[2]);
        std::abort();
    }
}

static void setenv_if_empty(const char* k, const char* v)
{
    if (!std::getenv(k))
        setenv(k, v, 1);
}

static void setup_openmp_defaults()
{
    // Choose threads = physical cores when we can; else logical
    const int n_phys = detect_physical_cores_linux();
    const int n_log = count_available_logical_cpus();
    const int n_thr = std::max(1, n_phys > 0 ? n_phys : n_log);

    // Respect user overrides; otherwise set sane defaults
    {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%d", n_thr);
        setenv_if_empty("OMP_NUM_THREADS", buf);
    }
    setenv_if_empty("OMP_PLACES", "cores");
    setenv_if_empty("OMP_PROC_BIND", "close");
    setenv_if_empty("OMP_DYNAMIC", "FALSE");

    // Apply programmatically too (in case env was empty)
    omp_set_dynamic(0);
    omp_set_num_threads(n_thr);
}

static void print_omp_summary()
{
    const char* nt = std::getenv("OMP_NUM_THREADS");
    const char* pl = std::getenv("OMP_PLACES");
    const char* pb = std::getenv("OMP_PROC_BIND");
    LOGI("OpenMP: threads=%s places=%s bind=%s (max_threads=%d)\n", (nt ? nt : "?"),
         (pl ? pl : "?"), (pb ? pb : "?"), omp_get_max_threads());
}

// --------- Pretty logging helpers (mini/compact/full) ----------
static const char* log_level()
{
    const char* s = std::getenv("SOLVER_LOG");
    if (!s || !*s)
        return "mini"; // default
    return s;
}

static void print_global_line(MPI_Comm comm, const AppConfig& cfg, const core::mesh::Mesh& m)
{
    int size = 1, rank = 0;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int topo = MPI_UNDEFINED;
    MPI_Topo_test(comm, &topo);
    int cd[3] = {1, 1, 1}, cp[3] = {0, 0, 0}, dummy[3] = {0, 0, 0};
    if (topo == MPI_CART)
        MPI_Cart_get(comm, 3, cd, cp, dummy);
    else
        MPI_Dims_create(size, 3, cd);

    long long NX = 0, NY = 0, NZ = 0;
    if (cfg.global.has_value())
    {
        NX = (*cfg.global)[0];
        NY = (*cfg.global)[1];
        NZ = (*cfg.global)[2];
    }
    else
    {
        NX = 1LL * m.local[0] * cd[0];
        NY = 1LL * m.local[1] * cd[1];
        NZ = 1LL * m.local[2] * cd[2];
    }
    const long long cells = NX * NY * NZ, u = (NX + 1) * NY * NZ, v = NX * (NY + 1) * NZ,
                    w = NX * NY * (NZ + 1);
    if (rank == 0)
    {
        LOGI("[run]   mpi=%d | omp=%d | dims=%dx%dx%d | per=%d,%d,%d | log=%s\n", size,
             omp_get_max_threads(), cd[0], cd[1], cd[2], cp[0], cp[1], cp[2], log_level());
        LOGI("[global] %lldx%lldx%lld | cells=%lld u=%lld v=%lld w=%lld\n", NX, NY, NZ, cells, u, v,
             w);
    }
}

// Gather per-rank lines to rank 0 and print as one block (prevents interleaving under CTest).
static void gather_and_print_lines(MPI_Comm comm, const std::string& line)
{
    int size = 1, rank = 0;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    int len = static_cast<int>(line.size()) + 1;
    std::vector<int> lens, displs;
    if (rank == 0)
        lens.resize(size), displs.resize(size);
    MPI_Gather(&len, 1, MPI_INT, rank == 0 ? lens.data() : nullptr, 1, MPI_INT, 0, comm);
    int total = 0;
    if (rank == 0)
    {
        displs.resize(size);
        for (int i = 0; i < size; ++i)
        {
            displs[i] = total;
            total += lens[i];
        }
    }
    std::vector<char> recvbuf(rank == 0 ? total : 0);
    MPI_Gatherv(line.c_str(), len, MPI_CHAR, rank == 0 ? recvbuf.data() : nullptr,
                rank == 0 ? lens.data() : nullptr, rank == 0 ? displs.data() : nullptr, MPI_CHAR, 0,
                comm);
    if (rank == 0)
    {
        for (int i = 0; i < size; ++i)
        {
            LOGD("%s\n", recvbuf.data() + displs[i]);
        }
    }
}

template <class T> static void parallel_fill(T* p, std::size_t n, T val)
{
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i)
        p[i] = val;
}

int main(int argc, char** argv)
{

    // Ensure MPI & PETSc lifetimes are handled exactly once and in the right order.
    MpiPetscOnce runtime(argc, argv);

    // ---------- Logging (rank0-gated INFO/DEBUG by default; SOLVER_LOG controls level) ----------
    core::master::logx::init(
        {core::master::logx::Level::Info, /*color*/ true, /*rank0_only*/ true});
    // Users can override with: SOLVER_LOG=quiet|error|warn|info|debug

    setup_openmp_defaults();
    // Print OMP summary once (rank 0) to avoid interleaved spam
    {
        int r0 = 1, rk = 0;
        MPI_Initialized(&r0);
        if (r0)
            MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        if (!r0 || rk == 0)
        {
            print_omp_summary();
        }
    }

    // 1) Parse YAML config
    const std::string cfg_path = (argc > 1) ? argv[1] : "case.yaml";
    const AppConfig cfg = load_config_from_yaml(cfg_path);

    // 2) Build a RunContext (MPI/stream/memory)
    RunContext rc{};
    {
        // This enables MPI halo exchange (skipped on non-Cartesian comms).
        MPI_Comm world = MPI_COMM_WORLD;

        // Build Cartesian dims with partial pin support:
        // Any positive entry in cfg.proc_grid is treated as "pinned".
        // Zeros mean "choose for me" and will be filled by MPI_Dims_create.
        int size = 1;
        MPI_Comm_size(world, &size);
        int dims[3] = {0, 0, 0};

        // Seed any pinned axes from YAML (e.g., [0,0,1] → dims=[0,0,1]).
        if (cfg.proc_grid[0] > 0)
            dims[0] = cfg.proc_grid[0];
        if (cfg.proc_grid[1] > 0)
            dims[1] = cfg.proc_grid[1];
        if (cfg.proc_grid[2] > 0)
            dims[2] = cfg.proc_grid[2];

        // If user fixed all 3 but product != size, warn and fall back to auto (keep zeros).
        if (dims[0] > 0 && dims[1] > 0 && dims[2] > 0)
        {
            const int prod = dims[0] * dims[1] * dims[2];
            if (prod != size)
            {
                int wrank = -1;
                MPI_Comm_rank(world, &wrank);
                if (wrank == 0)
                {
                    LOGW("[mpi] WARNING: proc_grid=%dx%dx%d (product=%d) != world size %d — "
                         "falling back to MPI_Dims_create().\n",
                         dims[0], dims[1], dims[2], prod, size);
                }
                dims[0] = dims[1] = dims[2] = 0;
            }
        }

        // Fill zeros near-cubically while preserving any pinned entries
        // (e.g., [0,0,1] with -np 8 → [2,4,1] or [4,2,1]).
        MPI_Dims_create(size, 3, dims);

        // Periodicity must match mesh.periodic[] so halos will wrap when requested.
        int periods[3] = {
            cfg.periodic[0] ? 1 : 0,
            cfg.periodic[1] ? 1 : 0,
            cfg.periodic[2] ? 1 : 0,
        };

        // Let MPI reorder for better placement if it wants to.
        int reorder = 1;
        MPI_Comm cart = MPI_COMM_NULL;
        MPI_Cart_create(world, 3, dims, periods, reorder, &cart);

        // Fallback: if Cart create failed (e.g., size==0?!), keep using WORLD.
        if (cart == MPI_COMM_NULL)
            cart = world;

        // Store the communicator we actually want to use everywhere.
        rc.mpi_comm = mpi_box(cart);
    }

    {
        // Unbox the communicator we actually intend to use everywhere (Cartesian)
        MPI_Comm comm = mpi_unbox(rc.mpi_comm);

        // Fallbacks in case something odd happened
        if (comm == MPI_COMM_NULL)
        {
            comm = MPI_COMM_WORLD;
            // Optional: re-box so downstream code sees a valid communicator
            if (rc.mpi_comm)
                mpi_box_free(rc.mpi_comm);
            rc.mpi_comm = mpi_box(comm);
        }

        // Fill in RunContext's rank/size from *that* communicator
        MPI_Comm_rank(comm, &rc.world_rank);
        MPI_Comm_size(comm, &rc.world_size);

        // Make errors return codes instead of aborts during bring-up
        MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    }

    rc.device_stream = nullptr; // set a real stream if/when you have a device
    rc.mem = &core::memory::MemoryManager::instance();

    // ---------------- Mesh from config (global→local derivation if requested) ----------------
    // 3) Mesh from config
    core::mesh::Mesh mesh;

    mesh.local = cfg.local;
    mesh.ng = cfg.ng;
    mesh.periodic = cfg.periodic;
    mesh.proc_grid = cfg.proc_grid;

    // If user supplied a GLOBAL size in YAML, derive rank-local extents via the communicator dims.
    // This keeps solver semantics (mesh.local is per-rank interior) while allowing human-friendly
    // input.
    if (cfg.global.has_value())
    {
        MPI_Comm comm = mpi_unbox(rc.mpi_comm);
        int cd[3] = {0, 0, 0}, cp[3] = {0, 0, 0}, dummy_coords[3] = {0, 0, 0};
        int topo = MPI_UNDEFINED;
        MPI_Topo_test(comm, &topo);
        if (topo == MPI_CART)
        {
            MPI_Cart_get(comm, 3, cd, cp, dummy_coords);
        }
        else
        {
            // Fall back to factorization (shouldn't happen here)
            MPI_Dims_create(rc.world_size, 3, cd);
        }

        const auto G = *cfg.global; // nx_g, ny_g, nz_g
        const int factors[3] = {std::max(cd[0], 1), std::max(cd[1], 1), std::max(cd[2], 1)};

        bool ok_div =
            (G[0] % factors[0] == 0) && (G[1] % factors[1] == 0) && (G[2] % factors[2] == 0);
        if (!ok_div)
        {
            int r = -1;
            MPI_Comm_rank(comm, &r);
            if (r == 0)
            {
                LOGE("[mesh] ERROR: mesh.global=%dx%dx%d not divisible by proc dims=%dx%dx%d\n",
                     G[0], G[1], G[2], factors[0], factors[1], factors[2]);
            }
            std::abort();
        }
        mesh.local = {G[0] / factors[0], G[1] / factors[1], G[2] / factors[2]};
    }

    {
        // Prefer the actual communicator dims (what halos/neighbors will use),
        // only backfilling non-positive entries so YAML pins remain honored.
        int cd[3] = {0, 0, 0}, cp[3] = {0, 0, 0}, cc[3] = {0, 0, 0};
        MPI_Comm comm = mpi_unbox(rc.mpi_comm);
        int topo = MPI_UNDEFINED;
        MPI_Topo_test(comm, &topo);
        if (topo == MPI_CART)
        {
            MPI_Cart_get(comm, 3, cd, cp, cc);
        }
        else
        {
            // Extremely unlikely here, but keep a reasonable fallback
            int sz = 1;
            MPI_Comm_size(comm, &sz);
            int tmp[3] = {0, 0, 0};
            MPI_Dims_create(sz, 3, tmp);
            cd[0] = tmp[0];
            cd[1] = tmp[1];
            cd[2] = tmp[2];
        }
        if (mesh.proc_grid[0] <= 0)
            mesh.proc_grid[0] = std::max(1, cd[0]);
        if (mesh.proc_grid[1] <= 0)
            mesh.proc_grid[1] = std::max(1, cd[1]);
        if (mesh.proc_grid[2] <= 0)
            mesh.proc_grid[2] = std::max(1, cd[2]);
    }

    // Sanity check: communicator dims must match mesh.proc_grid for correct halo neighbors.
    {
        MPI_Comm comm = mpi_unbox(rc.mpi_comm ? rc.mpi_comm : nullptr);
        if (comm != MPI_COMM_NULL)
        {
            int topo = MPI_UNDEFINED;
            MPI_Topo_test(comm, &topo);
            if (topo == MPI_CART)
            {
                int cd[3] = {0, 0, 0}, cp[3] = {0, 0, 0}, cc[3] = {0, 0, 0};
                MPI_Cart_get(comm, 3, cd, cp, cc);
                if (mesh.proc_grid[0] != cd[0] || mesh.proc_grid[1] != cd[1] ||
                    mesh.proc_grid[2] != cd[2])
                {
                    int r = -1;
                    MPI_Comm_rank(comm, &r);
                    if (r == 0)
                    {
                        LOGW("[mpi] WARNING: mesh.proc_grid=%dx%dx%d but communicator "
                             "dims=%dx%dx%d — ghost exchanges may use wrong neighbors.\n",
                             mesh.proc_grid[0], mesh.proc_grid[1], mesh.proc_grid[2], cd[0], cd[1],
                             cd[2]);
                    }
                }
            }
        }
    }

    // Build mesh, create Master, run simulation
    // Wrap PETSc users in a scope so their destructors run BEFORE PetscFinalize().
    {

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

        // Optional XDMF version (accept strings "v2"/"v3" or booly flag in YAML)
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

        // 8) Allocate and register fields (MAC-aware, via FieldCatalog creators)
        FieldCatalog& fc = master.fields();

        // If user didn't request outputs, pick a sensible fluids default and emit warning
        std::vector<std::string> want = cfg.fields_output;
        if (want.empty())
        {
            want = {"p", "u", "v", "w"}; // sensible default
            LOGW("No output fields selected, picked sensible defaults: p, u, v, w\n");
        };

        auto ensure = [&](const std::string& name)
        {
            if (name == "u")
                fc.create_face_scalar("u", /*axis=*/0, mesh);
            else if (name == "v")
                fc.create_face_scalar("v", /*axis=*/1, mesh);
            else if (name == "w")
                fc.create_face_scalar("w", /*axis=*/2, mesh);
            else
                fc.create_center_scalar(name, mesh);
            if (fc.contains(name))
                fc.select_for_output(name);
        };
        for (const auto& n : want)
            ensure(n);

        // ---- Tight, non-clobbered startup summary ----
        {
            MPI_Comm comm = mpi_unbox(rc.mpi_comm);
            int size = 1, rank = 0;
            MPI_Comm_size(comm, &size);
            MPI_Comm_rank(comm, &rank);

            // Single brief header + global line
            print_global_line(comm, cfg, mesh);

            // Build one short line per rank, gather on rank 0, then print the block.
            int topo = MPI_UNDEFINED;
            MPI_Topo_test(comm, &topo);
            int cd[3] = {1, 1, 1}, cp[3] = {0, 0, 0}, cc[3] = {0, 0, 0};
            if (topo == MPI_CART)
            {
                MPI_Cart_get(comm, 3, cd, cp, cc);
                MPI_Cart_coords(comm, rank, 3, cc);
            }
            else
            {
                MPI_Dims_create(size, 3, cd);
            }

            const int nx = mesh.local[0], ny = mesh.local[1], nz = mesh.local[2], ng = mesh.ng;
            const long long cells = 1LL * nx * ny * nz;
            const long long ucnt = 1LL * (nx + 1) * ny * nz;
            const long long vcnt = 1LL * nx * (ny + 1) * nz;
            const long long wcnt = 1LL * nx * ny * (nz + 1);

            std::ostringstream os;
            const std::string level = log_level();
            if (level == "full")
            {
                const int cx = nx + 2 * ng, cy = ny + 2 * ng, cz = nz + 2 * ng;
                const int ux = nx + 1 + 2 * ng, uy = cy, uz = cz;
                const int vx = cx, vy = ny + 1 + 2 * ng, vz = cz;
                const int wx = cx, wy = cy, wz = nz + 1 + 2 * ng;
                os << "[r" << rank << "/" << size << "] c(" << cc[0] << "," << cc[1] << "," << cc[2]
                   << ") d(" << cd[0] << "," << cd[1] << "," << cd[2] << ") " << "L " << nx << "x"
                   << ny << "x" << nz << " | N=" << cells << " U=" << ucnt << " V=" << vcnt
                   << " W=" << wcnt << " | ng=" << ng << " | Lg " << cx << "x" << cy << "x" << cz
                   << " Ug " << ux << "x" << uy << "x" << uz << " Vg " << vx << "x" << vy << "x"
                   << vz << " Wg " << wx << "x" << wy << "x" << wz;
            }
            else if (level == "compact")
            {
                os << "[r" << rank << "/" << size << "] c(" << cc[0] << "," << cc[1] << "," << cc[2]
                   << ") d(" << cd[0] << "," << cd[1] << "," << cd[2] << ") " << "L " << nx << "x"
                   << ny << "x" << nz << " | N=" << cells << " U=" << ucnt << " V=" << vcnt
                   << " W=" << wcnt << " | ng=" << ng;
            }
            else
            { // mini (default)
                os << "[r" << rank << "/" << size << "] " << "c(" << cc[0] << "," << cc[1] << ","
                   << cc[2] << ") " << "d(" << cd[0] << "," << cd[1] << "," << cd[2] << ") " << "L "
                   << nx << "x" << ny << "x" << nz << " | N=" << cells << " | ng=" << ng;
            }
            gather_and_print_lines(comm, os.str());
        }

        // Verify every field’s extents (selected for output is fine; or iterate over all)
        for (const auto& v : fc.selected_for_output())
        {
            assert_mac_extents(v, mesh);
        }

        // Prefetch selected fields to host (UM or mirrored paths do the right thing)
        for (const auto& v : fc.selected_for_output())
        {
            const std::size_t N = static_cast<std::size_t>(v.extents[0]) *
                                  static_cast<std::size_t>(v.extents[1]) *
                                  static_cast<std::size_t>(v.extents[2]) * v.elem_size;
            rc.mem->to_host(v.host_ptr, N, nullptr);
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
                std::span<const core::master::AnyFieldView>(views.data(), views.size()),
                override_b);

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
            int world = rc.world_size;
            auto [ok, msg] = run_preflight(wcfg, plan, world, avail_ram, avail_disk);
            LOGI("%s\n", msg.c_str());
            if (!ok)
                LOGW("Preflight warnings: consider reducing outputs or precision.\n");
        }

        // 10) Time controls & run
        TimeControls tc{};
        tc.dt = cfg.dt;
        tc.t_end = cfg.t_end;
        tc.write_every = cfg.write_every_steps;
        tc.case_name = cfg.case_name;
        master.run(tc);
    }

    return 0;
}
