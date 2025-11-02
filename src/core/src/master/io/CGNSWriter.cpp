#include "master/io/CGNSWriter.hpp"
#include "master/Views.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"

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
};

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

    // Fresh file each run
    if (cg_open(impl_->filepath.c_str(), CG_MODE_WRITE, &impl_->file))
        return;

    // Base
    int nbases = 0;
    if (ok(cg_nbases(impl_->file, &nbases)) && nbases >= 1)
    {
        impl_->base = 1;
    }
    else
    {
        const int cell_dim = 3, phys_dim = 3;
        if (cg_base_write(impl_->file, case_name.c_str(), cell_dim, phys_dim, &impl_->base))
            return;
    }
    ok(cg_simulation_type_write(impl_->file, impl_->base, TimeAccurate));

    impl_->zone = 0;
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
        // Use MIN across all selected fields to get the true cell counts.
        // (Face fields are +1 along their normal; MIN is the cell size.)
        impl_->nx = std::numeric_limits<int>::max();
        impl_->ny = std::numeric_limits<int>::max();
        impl_->nz = std::numeric_limits<int>::max();
        for (const auto& f : plan_.fields)
        {
            impl_->nx = std::min(impl_->nx, f.shape.nx);
            impl_->ny = std::min(impl_->ny, f.shape.ny);
            impl_->nz = std::min(impl_->nz, f.shape.nz);
        }
        if (impl_->nx == std::numeric_limits<int>::max())
            impl_->nx = 0;
        if (impl_->ny == std::numeric_limits<int>::max())
            impl_->ny = 0;
        if (impl_->nz == std::numeric_limits<int>::max())
            impl_->nz = 0;

        // Zone (create once)
        cgsize_t size[9] = {(cgsize_t) (impl_->nx + 1),
                            (cgsize_t) (impl_->ny + 1),
                            (cgsize_t) (impl_->nz + 1),
                            (cgsize_t) impl_->nx,
                            (cgsize_t) impl_->ny,
                            (cgsize_t) impl_->nz,
                            0,
                            0,
                            0};

        int nzones = 0;
        if (ok(cg_nzones(impl_->file, impl_->base, &nzones)) && nzones >= 1)
        {
            impl_->zone = 1;
        }
        else
        {
            if (cg_zone_write(impl_->file, impl_->base, "Zone", size, Structured, &impl_->zone))
                return;
        }

        // Static grid coordinates (unit lattice)
        const std::size_t vx = (std::size_t) impl_->nx + 1;
        const std::size_t vy = (std::size_t) impl_->ny + 1;
        const std::size_t vz = (std::size_t) impl_->nz + 1;
        const std::size_t nvert = vx * vy * vz;
        std::vector<double> X(nvert), Y(nvert), Z(nvert);
        std::size_t p = 0;
        for (std::size_t k = 0; k < vz; ++k)
            for (std::size_t j = 0; j < vy; ++j)
                for (std::size_t i = 0; i < vx; ++i, ++p)
                {
                    X[p] = double(i);
                    Y[p] = double(j);
                    Z[p] = double(k);
                }
        int gc;
        if (cg_grid_write(impl_->file, impl_->base, impl_->zone, "GridCoordinates", &gc))
            return;
        int cx, cy, cz;
        if (cg_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble, "CoordinateX",
                           X.data(), &cx))
            return;
        if (cg_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble, "CoordinateY",
                           Y.data(), &cy))
            return;
        if (cg_coord_write(impl_->file, impl_->base, impl_->zone, RealDouble, "CoordinateZ",
                           Z.data(), &cz))
            return;

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
        if (*slot == 0)
        {
            if (cg_sol_write(impl_->file, impl_->base, impl_->zone, name, loc, slot))
                return 0;
        }
        return *slot;
    };

    // Fields
    bool have_cell_solution = false;
    for (std::size_t i = 0; i < plan_.fields.size(); ++i)
    {
        const auto& fp = plan_.fields[i];
        const auto& view = req.fields[i];
        // Decide location from extents
        GridLocation_t loc =
            infer_location(fp.shape.nx, fp.shape.ny, fp.shape.nz, impl_->nx, impl_->ny, impl_->nz);
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

        int fld_id{};
        // Write the raw array under the correct GridLocation
        if (cg_field_write(impl_->file, impl_->base, impl_->zone, sol_id, dtype,
                           fp.shape.name.c_str(), staging, &fld_id) != CG_OK)
        {
            std::fprintf(stderr, "[CGNS] cg_field_write(%s) failed: %s\n", fp.shape.name.c_str(),
                         cg_get_error());
            return;
        }

        // If this is a face-centered velocity, also emit a CellCenter alias
        // named u_cell/v_cell/w_cell by averaging onto cells (ParaView-friendly & no collisions).
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

            // If alias already exists in the CellCenter solution for this step, skip it.
            int nfld = 0;
            bool alias_exists = false;
            if (cg_nfields(impl_->file, impl_->base, impl_->zone, sol_cell_id, &nfld) == CG_OK)
            {
                for (int fi = 1; fi <= nfld; ++fi)
                {
                    CGNS_ENUMT(DataType_t) fdtype;
                    char fname[33] = {0};
                    if (cg_field_info(impl_->file, impl_->base, impl_->zone, sol_cell_id, fi,
                                      &fdtype, fname) != CG_OK)
                        continue;
                    if (std::strcmp(fname, alias.c_str()) == 0)
                    {
                        alias_exists = true;
                        break;
                    }
                }
            }
            if (!alias_exists)
            {
                // Build averaged buffer into the same staging slot
                const bool src_f64 = (view.elem_size == 8);
                const bool dst_f64 = (fp.shape.elem_size == 8); // keep plan dtype for the alias
                void* tmp = pool_.data(i); // face buffer ≥ cell buffer, safe to reuse

                if (dst_f64)
                {
                    auto* d = static_cast<double*>(tmp);
                    if (loc == IFaceCenter)
                    {
                        src_f64 ? avg_i_to_cell<double, double>(d, view, impl_->nx, impl_->ny,
                                                                impl_->nz)
                                : avg_i_to_cell<float, double>(d, view, impl_->nx, impl_->ny,
                                                               impl_->nz);
                    }
                    else if (loc == JFaceCenter)
                    {
                        src_f64 ? avg_j_to_cell<double, double>(d, view, impl_->nx, impl_->ny,
                                                                impl_->nz)
                                : avg_j_to_cell<float, double>(d, view, impl_->nx, impl_->ny,
                                                               impl_->nz);
                    }
                    else
                    { // KFaceCenter
                        src_f64 ? avg_k_to_cell<double, double>(d, view, impl_->nx, impl_->ny,
                                                                impl_->nz)
                                : avg_k_to_cell<float, double>(d, view, impl_->nx, impl_->ny,
                                                               impl_->nz);
                    }
                }
                else
                {
                    auto* d = static_cast<float*>(tmp);
                    if (loc == IFaceCenter)
                    {
                        src_f64
                            ? avg_i_to_cell<double, float>(d, view, impl_->nx, impl_->ny, impl_->nz)
                            : avg_i_to_cell<float, float>(d, view, impl_->nx, impl_->ny, impl_->nz);
                    }
                    else if (loc == JFaceCenter)
                    {
                        src_f64
                            ? avg_j_to_cell<double, float>(d, view, impl_->nx, impl_->ny, impl_->nz)
                            : avg_j_to_cell<float, float>(d, view, impl_->nx, impl_->ny, impl_->nz);
                    }
                    else
                    { // KFaceCenter
                        src_f64
                            ? avg_k_to_cell<double, float>(d, view, impl_->nx, impl_->ny, impl_->nz)
                            : avg_k_to_cell<float, float>(d, view, impl_->nx, impl_->ny, impl_->nz);
                    }
                }

                int fld_cc{};
                if (cg_field_write(impl_->file, impl_->base, impl_->zone, sol_cell_id, dtype,
                                   alias.c_str(), tmp, &fld_cc) != CG_OK)
                {
                    std::fprintf(stderr, "[CGNS] cg_field_write(CellCenter alias %s) failed: %s\n",
                                 alias.c_str(), cg_get_error());
                    return;
                }
            }
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

    // Write final iterative metadata (all steps)
    if (impl_->zone != 0 && !impl_->times.empty())
    {
        finalize_iter_meta(impl_->file, impl_->base, impl_->zone, impl_->times, impl_->sol_names);
    }

    cg_close(impl_->file);
    impl_->opened = false;
}

} // namespace core::master::io
