#include "master/io/CGNSWriter.hpp"
#include "master/Views.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"
#include "memory/MpiBox.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include <cgnslib.h>
#include <mpi.h>
#include <pcgnslib.h>

namespace core::master::io
{
namespace fs = std::filesystem;

struct CGNSWriter::Impl
{
    int file = -1;
    int base = 1;
    int zone = 0;
    bool opened = false;
    bool planned = false;

    int step_index = 0;
    std::vector<double> times;
    std::vector<std::string> sol_names;

    int nx = 0, ny = 0, nz = 0;
    std::string filepath;
    int rank = 0, size = 1;
};

// -- Helper: get communicator dims & this-rank coords with a robust fallback -------------
static inline void cart_dims_coords(MPI_Comm in_comm, int rank, int dims[3], int coords[3])
{
    dims[0] = dims[1] = dims[2] = 0;
    coords[0] = coords[1] = coords[2] = 0;

    MPI_Comm comm = (in_comm == MPI_COMM_NULL) ? MPI_COMM_WORLD : in_comm;
    int topo = MPI_UNDEFINED;
    MPI_Topo_test(comm, &topo);
    if (topo == MPI_CART)
    {
        int periods[3] = {0, 0, 0}, dummy_coords[3] = {0, 0, 0};
        MPI_Cart_get(comm, 3, dims, periods, dummy_coords);
        MPI_Cart_coords(comm, rank, 3, coords);
        return;
    }

    // Fallback if not a Cartesian communicator: derive dims and flatten rank into coords.
    int size = 1;
    MPI_Comm_size(comm, &size);
    int td[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, td);
    dims[0] = td[0]; dims[1] = td[1]; dims[2] = td[2];
    int r = rank;
    coords[0] = (dims[0] > 0) ? (r % dims[0]) : 0; r /= (dims[0] > 0 ? dims[0] : 1);
    coords[1] = (dims[1] > 0) ? (r % dims[1]) : 0; r /= (dims[1] > 0 ? dims[1] : 1);
    coords[2] = (dims[2] > 0) ? (r % dims[2]) : 0;
}

static inline bool ok(int code)
{
    return code == CG_OK;
}

static inline bool is_contig(const AnyFieldView& v, std::size_t elem_size)
{
    const int nx = v.extents[0], ny = v.extents[1], nz = v.extents[2];
    const auto sx = v.strides[0], sy = v.strides[1], sz = v.strides[2];
    return sx == (std::ptrdiff_t) elem_size && sy == (std::ptrdiff_t)(elem_size * nx) &&
           sz == (std::ptrdiff_t)(elem_size * nx * ny);
}

// Small helper to infer CGNS GridLocation from field extents relative to the cell counts
static inline GridLocation_t infer_location(int fnx, int fny, int fnz, int nx, int ny, int nz)
{
    if (fnx == nx && fny == ny && fnz == nz)
        return CellCenter;
    if (fnx == nx + 1 && fny == ny && fnz == nz)
        return IFaceCenter;
    if (fnx == nx && fny == ny + 1 && fnz == nz)
        return JFaceCenter;
    if (fnx == nx && fny == ny && fnz == nz + 1)
        return KFaceCenter;
    // Fallback: treat anything else as Vertex to avoid lying; readers can still consume raw arrays
    return Vertex;
}

template <class SrcT, class DstT>
static void pack_cast_3d(DstT* __restrict d, const SrcT* __restrict s, int nx, int ny, int nz,
                         std::ptrdiff_t sx, std::ptrdiff_t sy, std::ptrdiff_t sz)
{
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            auto* line = reinterpret_cast<const unsigned char*>(s) + std::size_t(k) * sz +
                         std::size_t(j) * sy;
            for (int i = 0; i < nx; ++i)
            {
                const auto* sp = reinterpret_cast<const SrcT*>(line + std::size_t(i) * sx);
                *d++ = static_cast<DstT>(*sp);
            }
        }
    }
}

static void pack_bytes_3d(unsigned char* __restrict d, const unsigned char* __restrict s, int nx,
                          int ny, int nz, std::ptrdiff_t sx, std::ptrdiff_t sy, std::ptrdiff_t sz,
                          std::size_t elem_size)
{
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            auto* line = s + std::size_t(k) * sz + std::size_t(j) * sy;
            for (int i = 0; i < nx; ++i)
            {
                std::memcpy(d, line + std::size_t(i) * sx, elem_size);
                d += elem_size;
            }
        }
    }
}

static WritePlan make_plan_from_request(const std::vector<AnyFieldView>& views,
                                        WriterConfig::Precision prec)
{
    std::size_t override_b = (prec == WriterConfig::Precision::Float32   ? 4
                              : prec == WriterConfig::Precision::Float64 ? 8
                                                                         : 0);
    return build_write_plan(std::span<const AnyFieldView>(views.data(), views.size()), override_b);
}

// ---------- helpers: average face-centered → cell-centered ----------
// NOTE: v.extents[] include ghosts; nx,ny,nz are *cell* counts (interior).
// We must offset the source pointer by the ghost width before averaging.
template <class SrcT, class DstT>
static inline void avg_i_to_cell(DstT* dst, const AnyFieldView& v, int nx, int ny, int nz)
{
    const auto sx = v.strides[0] / sizeof(SrcT), sy = v.strides[1] / sizeof(SrcT),
               sz = v.strides[2] / sizeof(SrcT);
    // IFace has (nx+1) interior in x; ghosts in y,z only
    const int ngx = (v.extents[0] - (nx + 1)) / 2;
    const int ngy = (v.extents[1] - ny) / 2;
    const int ngz = (v.extents[2] - nz) / 2;
    const SrcT* base = static_cast<const SrcT*>(v.host_ptr);
    const SrcT* s = base + ngx * sx + ngy * sy + ngz * sz;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                dst[(k * ny + j) * nx + i] = DstT(
                    (s[k * sz + j * sy + (i) *sx] + s[k * sz + j * sy + (i + 1) * sx]) * SrcT(0.5));
}
template <class SrcT, class DstT>
static inline void avg_j_to_cell(DstT* dst, const AnyFieldView& v, int nx, int ny, int nz)
{
    const auto sx = v.strides[0] / sizeof(SrcT), sy = v.strides[1] / sizeof(SrcT),
               sz = v.strides[2] / sizeof(SrcT);
    // JFace has (ny+1) interior in y; ghosts in x,z only
    const int ngx = (v.extents[0] - nx) / 2;
    const int ngy = (v.extents[1] - (ny + 1)) / 2;
    const int ngz = (v.extents[2] - nz) / 2;
    const SrcT* base = static_cast<const SrcT*>(v.host_ptr);
    const SrcT* s = base + ngx * sx + ngy * sy + ngz * sz;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                dst[(k * ny + j) * nx + i] = DstT(
                    (s[k * sz + j * sy + i * sx] + s[k * sz + (j + 1) * sy + i * sx]) * SrcT(0.5));
}
template <class SrcT, class DstT>
static inline void avg_k_to_cell(DstT* dst, const AnyFieldView& v, int nx, int ny, int nz)
{
    const auto sx = v.strides[0] / sizeof(SrcT), sy = v.strides[1] / sizeof(SrcT),
               sz = v.strides[2] / sizeof(SrcT);
    // KFace has (nz+1) interior in z; ghosts in x,y only
    const int ngx = (v.extents[0] - nx) / 2;
    const int ngy = (v.extents[1] - ny) / 2;
    const int ngz = (v.extents[2] - (nz + 1)) / 2;
    const SrcT* base = static_cast<const SrcT*>(v.host_ptr);
    const SrcT* s = base + ngx * sx + ngy * sy + ngz * sz;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                dst[(k * ny + j) * nx + i] = DstT(
                    (s[k * sz + j * sy + i * sx] + s[(k + 1) * sz + j * sy + i * sx]) * SrcT(0.5));
}

// ---- iterative metadata (finalize at close) ------------------------------------

static void finalize_iter_meta(int fn, int base, int zone, const std::vector<double>& times,
                               const std::vector<std::string>& sol_names)
{
    const int nsteps = (int) times.size();
    if (nsteps <= 0)
        return;

    // strict monotone: keep physical times, only nudge ties
    std::vector<double> tvals(times.begin(), times.end());
    for (int i = 1; i < nsteps; ++i)
    {
        if (!(tvals[i] > tvals[i - 1]))
        {
            const double eps = 1e-12 * (std::abs(tvals[i - 1]) + 1.0);
            tvals[i] = tvals[i - 1] + eps;
        }
    }

    // BaseIterativeData: recreate cleanly with final NumberOfSteps
    if (ok(cg_goto(fn, base, "end")))
    {
        cg_delete_node("BaseIterativeData");
    }
    ok(cg_biter_write(fn, base, "BaseIterativeData", nsteps));
    if (ok(cg_goto(fn, base, "BaseIterativeData_t", 1, "end")))
    {
        cg_delete_node("TimeValues");
        cg_delete_node("IterationValues");
        cgsize_t n = (cgsize_t) nsteps;
        ok(cg_array_write("TimeValues", RealDouble, 1, &n, tvals.data()));
        std::vector<int> iters(nsteps);
        for (int i = 0; i < nsteps; ++i)
            iters[i] = i + 1;
        ok(cg_array_write("IterationValues", Integer, 1, &n, iters.data()));
    }

    // ZoneIterativeData: recreate with final FlowSolutionPointers
    if (ok(cg_goto(fn, base, "Zone_t", zone, "end")))
    {
        cg_delete_node("ZoneIterativeData");
    }
    ok(cg_ziter_write(fn, base, zone, "ZoneIterativeData"));
    if (ok(cg_goto(fn, base, "Zone_t", zone, "ZoneIterativeData_t", 1, "end")))
    {
        cg_delete_node("FlowSolutionPointers");
        // leave GridCoordinatesPointers out for static grid

        constexpr int WIDTH = 32; // fixed-width, Fortran-style, space-padded
        const cgsize_t dims[2] = {(cgsize_t) WIDTH, (cgsize_t) nsteps};

        std::vector<char> solptrs((size_t) WIDTH * nsteps, ' ');
        for (int s = 0; s < nsteps; ++s)
        {
            const std::string& nm = (s < (int) sol_names.size()) ? sol_names[s] : sol_names.back();
            const size_t L = std::min(nm.size(), (size_t) WIDTH);
            std::memcpy(&solptrs[(size_t) s * WIDTH], nm.data(), L);
        }
        ok(cg_array_write("FlowSolutionPointers", Character, 2, dims, solptrs.data()));
    }
}

// --------------------------------------------------------------------------------

CGNSWriter::CGNSWriter(WriterConfig cfg) : cfg_(cfg), impl_(new Impl) {}
CGNSWriter::~CGNSWriter()
{
    close();
}

void CGNSWriter::open_case(const std::string& case_name)
{
    fs::create_directories(cfg_.path);
    impl_->filepath = cfg_.path + "/" + case_name + ".cgns";

    // Establish MPI rank/size for rank-0 responsibilities
    {
        MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm); // may be a Cartesian comm
        if (comm == MPI_COMM_NULL)
            comm = MPI_COMM_WORLD; // fallback
        MPI_Comm_rank(comm, &impl_->rank);
        MPI_Comm_size(comm, &impl_->size);
        cgp_mpi_comm(comm);
        cg_set_file_type(CG_FILE_HDF5);
    }

    
    // PARALLEL open (all ranks)
    {
        MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm);
        if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

        // Open brand-new file collectively in WRITE mode and create metadata in parallel session
        int ok_local = (cgp_open(impl_->filepath.c_str(), CG_MODE_WRITE, &impl_->file) == CG_OK);
        int ok_all = 0;
        MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm);
        if (!ok_all) {
            if (!ok_local) {
                std::fprintf(stderr, "[CGNS] cgp_open(%s, WRITE) failed: %s\n",
                             impl_->filepath.c_str(), cg_get_error());
            }
            // Ensure everyone leaves consistently
            impl_->opened = false;
            return;
        }

        // Create Base/Zone COLLECTIVELY (all ranks enter cg_* with identical args)
        int base_id = 0, zone_id = 0;
        {
            const int cell_dim = 3, phys_dim = 3;
            if (cg_base_write(impl_->file, "Base", cell_dim, phys_dim, &base_id) != CG_OK) {
                std::fprintf(stderr, "[CGNS] cg_base_write failed: %s\n", cg_get_error());
                base_id = 0;
            }
            const int GX_cells = cfg_.mesh->global[0];
            const int GY_cells = cfg_.mesh->global[1];
            const int GZ_cells = cfg_.mesh->global[2];
            cgsize_t size[9] = {
                (cgsize_t)(GX_cells + 1), (cgsize_t)(GY_cells + 1), (cgsize_t)(GZ_cells + 1),
                (cgsize_t) GX_cells,      (cgsize_t) GY_cells,      (cgsize_t) GZ_cells,
                0, 0, 0};
            if (base_id != 0 &&
                cg_zone_write(impl_->file, base_id, "Zone", size, Structured, &zone_id) != CG_OK) {
                std::fprintf(stderr, "[CGNS] cg_zone_write failed: %s\n", cg_get_error());
                zone_id = 0;
            }
            if (base_id != 0) {
                ok(cg_simulation_type_write(impl_->file, base_id, TimeAccurate));
            }
        }
        // If root failed to create metadata, all ranks must back out collectively.
        int meta_ok = (base_id != 0 && zone_id != 0);
        int meta_all = 0;
        MPI_Allreduce(&meta_ok, &meta_all, 1, MPI_INT, MPI_MIN, comm);
        if (!meta_all) {
            if (impl_->file >= 0) cgp_close(impl_->file);
            impl_->opened = false;
            return;
        }
        impl_->base = base_id;
        impl_->zone = zone_id;

        MPI_Barrier(comm);
    }

    impl_->step_index = 0;
    impl_->times.clear();
    impl_->sol_names.clear();

    impl_->opened = true;
    impl_->planned = false;
}

void CGNSWriter::write(const WriteRequest& req)
{
    if (!impl_->opened)
        return;

    if (!impl_->planned)
    {
        plan_ = make_plan_from_request(req.fields, cfg_.precision);
        // Cache *global* cell counts for shapes / aliasing
        const int GX_cells = cfg_.mesh->global[0];
        const int GY_cells = cfg_.mesh->global[1];
        const int GZ_cells = cfg_.mesh->global[2];
        impl_->nx = GX_cells;
        impl_->ny = GY_cells;
        impl_->nz = GZ_cells;

        // Static grid coordinates (unit lattice) -- create node once, then write in parallel
        const std::size_t vx = (std::size_t) GX_cells + 1;
        const std::size_t vy = (std::size_t) GY_cells + 1;
        const std::size_t vz = (std::size_t) GZ_cells + 1;

        // Create coord arrays (collective) via cgp_
        int cx = 0, cy = 0, cz = 0;
        {
            MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm);
            if (comm == MPI_COMM_NULL)
                comm = MPI_COMM_WORLD;
            int ok_local = (cgp_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble,
                                            "CoordinateX", &cx) == CG_OK);
            int ok_all = 0;
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write(CoordinateX) failed: %s\n",
                                 cg_get_error());
                return;
            }
            ok_local = (cgp_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble,
                                        "CoordinateY", &cy) == CG_OK);
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write(CoordinateY) failed: %s\n",
                                 cg_get_error());
                return;
            }
            ok_local = (cgp_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble,
                                        "CoordinateZ", &cz) == CG_OK);
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write(CoordinateZ) failed: %s\n",
                                 cg_get_error());
                return;
            }
        }
        // All ranks: write their vertex hyperslab. On each axis, only the tile that actually
        // touches the global high end contributes the extra "+1" vertex. This detection is based
        // purely on tile extents (ox+lx vs global size), making it robust to any communicator
        // shape/reordering.
        {

            const int ox = cfg_.mesh->global_lo[0];
            const int oy = cfg_.mesh->global_lo[1];
            const int oz = cfg_.mesh->global_lo[2];
            const int lx = cfg_.mesh->local[0];
            const int ly = cfg_.mesh->local[1];
            const int lz = cfg_.mesh->local[2];

            // Global cell counts (cached earlier)
            const int GX = impl_->nx;
            const int GY = impl_->ny;
            const int GZ = impl_->nz;

            // Decide "+1" strictly from geometric coverage of the global high end.
            const bool last_x = (ox + lx == GX);
            const bool last_y = (oy + ly == GY);
            const bool last_z = (oz + lz == GZ);

            const int vxn = lx + (last_x ? 1 : 0);
            const int vyn = ly + (last_y ? 1 : 0);
            const int vzn = lz + (last_z ? 1 : 0);

            std::vector<double> X_local((size_t) vxn * vyn * vzn),
                Y_local((size_t) vxn * vyn * vzn), Z_local((size_t) vxn * vyn * vzn);
            size_t p = 0;
            for (int kz = 0; kz < vzn; ++kz)
                for (int jy = 0; jy < vyn; ++jy)
                    for (int ix = 0; ix < vxn; ++ix, ++p)
                    {
                        X_local[p] = double(ox + ix);
                        Y_local[p] = double(oy + jy);
                        Z_local[p] = double(oz + kz);
                    }
            const cgsize_t rminV[3] = {(cgsize_t) (ox + 1), (cgsize_t) (oy + 1),
                                       (cgsize_t) (oz + 1)};
            const cgsize_t rmaxV[3] = {(cgsize_t) (ox + vxn), (cgsize_t) (oy + vyn),
                                       (cgsize_t) (oz + vzn)};
            // Coordinate slabs must also succeed on all ranks collectively.
            MPI_Comm comm2 = mpi_unbox(cfg_.mpi_cart_comm);
            if (comm2 == MPI_COMM_NULL) comm2 = MPI_COMM_WORLD;
            int ok_local = (cgp_coord_write_data(impl_->file, impl_->base, impl_->zone, cx, rminV,
                                                 rmaxV, X_local.data()) == CG_OK);
            int ok_all = 0;
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm2);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write_data(X) failed: %s\n",
                                 cg_get_error());
                return;
            }
            ok_local = (cgp_coord_write_data(impl_->file, impl_->base, impl_->zone, cy, rminV,
                                             rmaxV, Y_local.data()) == CG_OK);
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm2);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write_data(Y) failed: %s\n",
                                 cg_get_error());
                return;
            }
            ok_local = (cgp_coord_write_data(impl_->file, impl_->base, impl_->zone, cz, rminV,
                                             rmaxV, Z_local.data()) == CG_OK);
            MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, comm2);
            if (!ok_all)
            {
                if (!ok_local)
                    std::fprintf(stderr, "[CGNS] cgp_coord_write_data(Z) failed: %s\n",
                                 cg_get_error());
                return;
            }
        }

        std::size_t maxb = 0;
        for (auto& f : plan_.fields)
            maxb = std::max(maxb, f.bytes);
        pool_.reserve(plan_.fields.size(), maxb);

        impl_->planned = true;
    }

    // Step bookkeeping
    const int s = impl_->step_index++;
    impl_->times.push_back(req.time);

    char solname[64];
    std::snprintf(solname, sizeof(solname), "FlowSolutionAtStep%06d", s + 1);
    impl_->sol_names.emplace_back(solname);

    // We may need multiple FlowSolution nodes per step (one per GridLocation)
    // Name them deterministically and lazily create on first use.
    // Keep a canonical solution name for iterative metadata (prefer CellCenter).
    int sol_cell = 0, sol_i = 0, sol_j = 0, sol_k = 0, sol_other = 0;

    // Helper (COLLECTIVE): find FlowSolution id by name (linear scan).
    auto find_solution_id_root = [&](const char* wanted) -> int
    {
        // All ranks must enter these cg_* calls collectively in a pcgns session.
        int nsol = 0, id = 0;
        if (cg_nsols(impl_->file, impl_->base, impl_->zone, &nsol) != CG_OK) nsol = 0;
        for (int si = 1; si <= nsol; ++si) {
            char nm[256] = {0};
            GridLocation_t l = Vertex;
            if (cg_sol_info(impl_->file, impl_->base, impl_->zone, si, nm, &l) != CG_OK) continue;
            if (std::strcmp(nm, wanted) == 0) { id = si; break; }
        }
        return id;
    };

    auto get_or_make_sol = [&](GridLocation_t loc) -> int
    {
        const char* tag = (loc == CellCenter)    ? "Cell"
                          : (loc == IFaceCenter) ? "IFace"
                          : (loc == JFaceCenter) ? "JFace"
                          : (loc == KFaceCenter) ? "KFace"
                                                 : "Vertex";
        char name[80];
        std::snprintf(name, sizeof(name), "%s_%s", solname, tag);
        int* slot = (loc == CellCenter)    ? &sol_cell
                    : (loc == IFaceCenter) ? &sol_i
                    : (loc == JFaceCenter) ? &sol_j
                    : (loc == KFaceCenter) ? &sol_k
                                           : &sol_other;
        if (*slot != 0)
            return *slot;

        MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm);
        if (comm == MPI_COMM_NULL)
            comm = MPI_COMM_WORLD;

        // COLLECTIVE creation: first check if it exists (all ranks call the same query),
        // then if missing, ALL ranks enter cg_sol_write with identical args.
        int id = 0;
        {
            int nsol = 0;
            if (cg_nsols(impl_->file, impl_->base, impl_->zone, &nsol) != CG_OK)
                nsol = 0;
            for (int si = 1; si <= nsol; ++si)
            {
                char nm[256] = {0};
                GridLocation_t l = Vertex;
                if (cg_sol_info(impl_->file, impl_->base, impl_->zone, si, nm, &l) != CG_OK)
                    continue;
                if (std::strcmp(nm, name) == 0) { id = si; break; }
            }
            int need_create = (id == 0);
            int need_all = 0;
            MPI_Allreduce(&need_create, &need_all, 1, MPI_INT, MPI_MAX, comm);
            if (need_all) {
                if (cg_sol_write(impl_->file, impl_->base, impl_->zone, name, loc, &id) != CG_OK) {
                    std::fprintf(stderr, "[CGNS] cg_sol_write(%s) failed: %s\n", name, cg_get_error());
                    id = 0;
                }
            }
        }
        if (id == 0)
            return 0;
        *slot = id;
        return id;
    };

    // Fields
    bool have_cell_solution = false;
    for (std::size_t i = 0; i < plan_.fields.size(); ++i)
    {
        const auto& fp = plan_.fields[i];
        const auto& view = req.fields[i];
        // Decide GridLocation from LOCAL tile sizes (not globals).
        const int lx_loc = cfg_.mesh->local[0];
        const int ly_loc = cfg_.mesh->local[1];
        const int lz_loc = cfg_.mesh->local[2];
        GridLocation_t loc = CellCenter;
        if (fp.shape.nx == lx_loc + 1 && fp.shape.ny == ly_loc && fp.shape.nz == lz_loc)
            loc = IFaceCenter;
        else if (fp.shape.nx == lx_loc && fp.shape.ny == ly_loc + 1 && fp.shape.nz == lz_loc)
            loc = JFaceCenter;
        else if (fp.shape.nx == lx_loc && fp.shape.ny == ly_loc && fp.shape.nz == lz_loc + 1)
            loc = KFaceCenter;
        else if (!(fp.shape.nx == lx_loc && fp.shape.ny == ly_loc && fp.shape.nz == lz_loc))
            loc = Vertex; // conservative fallback if shapes don't match expectations
        int sol_id = get_or_make_sol(loc);
        if (sol_id == 0)
            return;
        if (loc == CellCenter)
            have_cell_solution = true;
        DataType_t dtype = (fp.shape.elem_size == 8)   ? RealDouble
                           : (fp.shape.elem_size == 4) ? RealSingle
                                                       : DataTypeNull;

        void* staging = pool_.data(i);

        const bool same_size = (fp.shape.elem_size == view.elem_size);
        const bool contig = is_contig(view, view.elem_size);

        if (same_size && contig)
        {
            std::memcpy(staging, view.host_ptr, fp.bytes);
        }
        else if (view.elem_size == 8 && fp.shape.elem_size == 4)
        {
            pack_cast_3d<double, float>(static_cast<float*>(staging),
                                        static_cast<const double*>(view.host_ptr), fp.shape.nx,
                                        fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                                        view.strides[2]);
        }
        else if (view.elem_size == 4 && fp.shape.elem_size == 8)
        {
            pack_cast_3d<float, double>(static_cast<double*>(staging),
                                        static_cast<const float*>(view.host_ptr), fp.shape.nx,
                                        fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                                        view.strides[2]);
        }
        else
        {
            pack_bytes_3d(static_cast<unsigned char*>(staging),
                          static_cast<const unsigned char*>(view.host_ptr), fp.shape.nx,
                          fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                          view.strides[2], fp.shape.elem_size);
        }

        // Parallel CGNS does not support IFace/JFace/KFace FlowSolution locations on structured
        // zones. For face-centered velocities, skip writing the native face field and go straight
        // to the CellCenter alias (u_cell/v_cell/w_cell) below.
        int fld_id{};
        const bool is_face_loc = (loc == IFaceCenter || loc == JFaceCenter || loc == KFaceCenter);
        bool wrote_native = false;
        if (!is_face_loc)
        {
            if (cgp_field_write(impl_->file, impl_->base, impl_->zone, sol_id, dtype,
                                fp.shape.name.c_str(), &fld_id) != CG_OK)
            {
                std::fprintf(stderr, "[CGNS] cgp_field_write(%s) failed: %s\n",
                             fp.shape.name.c_str(), cg_get_error());
                // Don't abort the whole write step; continue to any CellCenter alias below.
            }
            else
            {
                wrote_native = true;
            }
        }
        // Compute rmin/rmax (Fortran 1-based). Use mesh globals and local sizes (+1 on stagger
        // axis).
        const int ox = cfg_.mesh->global_lo[0], oy = cfg_.mesh->global_lo[1],
                  oz = cfg_.mesh->global_lo[2];
        int lx = fp.shape.nx, ly = fp.shape.ny, lz = fp.shape.nz; // already interior (+1 if face)
        cgsize_t rmin[3] = {(cgsize_t) (ox + 1), (cgsize_t) (oy + 1), (cgsize_t) (oz + 1)};
        cgsize_t rmax[3] = {(cgsize_t) (ox + lx), (cgsize_t) (oy + ly), (cgsize_t) (oz + lz)};
        if (wrote_native)
        {
            if (cgp_field_write_data(impl_->file, impl_->base, impl_->zone, sol_id, fld_id, rmin,
                                     rmax, staging) != CG_OK)
            {
                std::fprintf(stderr, "[CGNS] cgp_field_write_data(%s) failed: %s\n",
                             fp.shape.name.c_str(), cg_get_error());
                // Keep going; we'll still try to emit the CellCenter alias for velocities.
            }
        }

        // If this is a face-centered velocity, also emit a CellCenter alias
        // named u_cell/v_cell/w_cell by averaging onto cells (ParaView-friendly).
        if (loc != CellCenter &&
            (fp.shape.name == "u" || fp.shape.name == "v" || fp.shape.name == "w"))
        {
            // Ensure CellCenter FlowSolution for this step exists.
            int sol_cell_id = get_or_make_sol(CellCenter);
            if (!sol_cell_id)
                return;
            have_cell_solution = true;

            // Alias name (avoid duplicate "u" etc. inside the same FlowSolution).
            const std::string alias = (fp.shape.name == "u"   ? "u_cell"
                                       : fp.shape.name == "v" ? "v_cell"
                                                              : "w_cell");

            {
                // Build averaged buffer into the same staging slot
                const bool src_f64 = (view.elem_size == 8);
                const bool dst_f64 = (fp.shape.elem_size == 8); // keep plan dtype for the alias
                void* tmp = pool_.data(i); // face buffer ≥ cell buffer, safe to reuse

                // local interior sizes for this rank's tile
                const int nx_loc = cfg_.mesh->local[0];
                const int ny_loc = cfg_.mesh->local[1];
                const int nz_loc = cfg_.mesh->local[2];

                if (dst_f64)
                {
                    auto* d = static_cast<double*>(tmp);
                    if (loc == IFaceCenter)
                    {
                        src_f64 ? avg_i_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_i_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                    else if (loc == JFaceCenter)
                    {
                        src_f64 ? avg_j_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_j_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                    else
                    { // KFaceCenter
                        src_f64 ? avg_k_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_k_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                }
                else
                {
                    auto* d = static_cast<float*>(tmp);
                    if (loc == IFaceCenter)
                    {
                        src_f64 ? avg_i_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_i_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                    else if (loc == JFaceCenter)
                    {
                        src_f64 ? avg_j_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_j_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                    else
                    { // KFaceCenter
                        src_f64 ? avg_k_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                                : avg_k_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                    }
                }

                // Create alias field and write our cell-centered tile
                int fld_cc{};
                if (cgp_field_write(impl_->file, impl_->base, impl_->zone, sol_cell_id, dtype,
                                    alias.c_str(), &fld_cc) != CG_OK)
                {
                    std::fprintf(stderr, "[CGNS] cgp_field_write(%s) failed: %s\n", alias.c_str(),
                                 cg_get_error());
                    return;
                }
                // Cell-centered alias write ranges are strictly the local (cell) tile:
                cgsize_t rmin_cc[3] = {(cgsize_t) (cfg_.mesh->global_lo[0] + 1),
                                       (cgsize_t) (cfg_.mesh->global_lo[1] + 1),
                                       (cgsize_t) (cfg_.mesh->global_lo[2] + 1)};
                cgsize_t rmax_cc[3] = {(cgsize_t) (cfg_.mesh->global_lo[0] + cfg_.mesh->local[0]),
                                       (cgsize_t) (cfg_.mesh->global_lo[1] + cfg_.mesh->local[1]),
                                       (cgsize_t) (cfg_.mesh->global_lo[2] + cfg_.mesh->local[2])};
                if (cgp_field_write_data(impl_->file, impl_->base, impl_->zone, sol_cell_id, fld_cc,
                                         rmin_cc, rmax_cc, tmp) != CG_OK)
                {
                    std::fprintf(stderr, "[CGNS] cgp_field_write_data(CellCenter %s) failed: %s\n",
                                 alias.c_str(), cg_get_error());
                    return;
                }
            } // alias creation & write (collective)
        }
    }

    // NOTE: we keep only one pointer per step in iterative metadata.
    // Prefer Cell solution; otherwise point to the first created one.
    if (have_cell_solution)
    {
        char solname_cell[80];
        std::snprintf(solname_cell, sizeof(solname_cell), "%s_Cell", solname);
        impl_->sol_names.back() = solname_cell;
    }
    // else impl_->sol_names already contains the base name; finalize_iter_meta will use it.
}

void CGNSWriter::close()
{
    if (!impl_ || !impl_->opened)
        return;

    // Use the same communicator the writer was initialized with.
    MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm);
    if (comm == MPI_COMM_NULL)
        comm = MPI_COMM_WORLD;

    // Everyone arrives before serial metadata update
    MPI_Barrier(comm);

    // FINAL ITERATIVE METADATA MUST BE COLLECTIVE IN A PARALLEL SESSION
    if (impl_->zone != 0 && !impl_->times.empty())
    {
        finalize_iter_meta(impl_->file, impl_->base, impl_->zone, impl_->times, impl_->sol_names);
    }

    // Ensure metadata is visible before collective close
    MPI_Barrier(comm);

    cgp_close(impl_->file);
    impl_->opened = false;
}

} // namespace core::master::io
