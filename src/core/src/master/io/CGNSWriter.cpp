#include "master/io/CGNSWriter.hpp"
#include "master/Views.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#ifdef SOLVERLES_WITH_CGNS
#include <cgnslib.h>
#endif

namespace core::master::io
{
namespace fs = std::filesystem;

/// \cond DOXYGEN_EXCLUDE

struct CGNSWriter::Impl
{
#ifdef SOLVERLES_WITH_CGNS
    int file = -1;
    int base = 1;
    int zone = 1;
    bool have_coords = false;
#endif
    std::string filepath;
    bool opened = false;
    bool planned = false;
    int step_index = 0;
    std::vector<double> times;
    int nx = 0, ny = 0, nz = 0;
};

static inline bool is_contig(const AnyFieldView& v, std::size_t elem_size)
{
    const int nx = v.extents[0], ny = v.extents[1], nz = v.extents[2];
    const auto sx = v.strides[0], sy = v.strides[1], sz = v.strides[2];
    return sx == (std::ptrdiff_t) elem_size && sy == (std::ptrdiff_t)(elem_size * nx) &&
           sz == (std::ptrdiff_t)(elem_size * nx * ny);
}

template <class SrcT, class DstT>
static void pack_cast_3d(DstT* __restrict d, const SrcT* __restrict s, int nx, int ny, int nz,
                         std::ptrdiff_t sx, std::ptrdiff_t sy, std::ptrdiff_t sz)
{
    for (int k = 0; k < nz; ++k)
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

static void pack_bytes_3d(unsigned char* __restrict d, const unsigned char* __restrict s, int nx,
                          int ny, int nz, std::ptrdiff_t sx, std::ptrdiff_t sy, std::ptrdiff_t sz,
                          std::size_t elem_size)
{
    for (int k = 0; k < nz; ++k)
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

static WritePlan make_plan_from_request(const std::vector<AnyFieldView>& views,
                                        WriterConfig::Precision prec)
{
    std::size_t override_b = (prec == WriterConfig::Precision::Float32   ? 4
                              : prec == WriterConfig::Precision::Float64 ? 8
                                                                         : 0);
    return build_write_plan(std::span<const AnyFieldView>(views.data(), views.size()), override_b);
}

/// \endcond

CGNSWriter::CGNSWriter(WriterConfig cfg) : cfg_(cfg), impl_(new Impl) {}
CGNSWriter::~CGNSWriter()
{
    close();
}

void CGNSWriter::open_case(const std::string& case_name)
{
#ifdef SOLVERLES_WITH_CGNS
    fs::create_directories(cfg_.path);
    impl_->filepath = cfg_.path + "/" + case_name + ".cgns";
    if (cg_open(impl_->filepath.c_str(), CG_MODE_WRITE, &impl_->file))
    {
        return;
    }
    const int cell_dim = 3, phys_dim = 3;
    if (cg_base_write(impl_->file, case_name.c_str(), cell_dim, phys_dim, &impl_->base))
    {
        return;
    }
    impl_->opened = true;
    impl_->step_index = 0;
    impl_->times.clear();
#else
    (void) case_name;
#endif
}

void CGNSWriter::write(const WriteRequest& req)
{
#ifdef SOLVERLES_WITH_CGNS
    if (!impl_->opened)
        return;

    if (!impl_->planned)
    {
        plan_ = make_plan_from_request(req.fields, cfg_.precision);
        impl_->nx = plan_.fields.front().shape.nx;
        impl_->ny = plan_.fields.front().shape.ny;
        impl_->nz = plan_.fields.front().shape.nz;

        cgsize_t size[9] = {(cgsize_t) (impl_->nx + 1),
                            (cgsize_t) (impl_->ny + 1),
                            (cgsize_t) (impl_->nz + 1),
                            (cgsize_t) impl_->nx,
                            (cgsize_t) impl_->ny,
                            (cgsize_t) impl_->nz,
                            0,
                            0,
                            0};
        if (cg_zone_write(impl_->file, impl_->base, "Zone", size, Structured, &impl_->zone))
            return;

        // coordinates once (unit spacing)
        const std::size_t vx = (std::size_t) impl_->nx + 1, vy = (std::size_t) impl_->ny + 1,
                          vz = (std::size_t) impl_->nz + 1;
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
        impl_->have_coords = true;

        cg_biter_write(impl_->file, impl_->base, "BaseIterativeData", 0);
        cg_ziter_write(impl_->file, impl_->base, impl_->zone, "ZoneIterativeData");

        std::size_t maxb = 0;
        for (auto& f : plan_.fields)
            maxb = std::max(maxb, f.bytes);
        pool_.reserve(plan_.fields.size(), maxb);
        impl_->planned = true;
    }

    impl_->times.push_back(req.time);
    char solname[64];
    std::snprintf(solname, sizeof(solname), "FlowSolution_%06d", req.step);
    int sol_id{};
    if (cg_sol_write(impl_->file, impl_->base, impl_->zone, solname, CellCenter, &sol_id))
        return;

    for (std::size_t i = 0; i < plan_.fields.size(); ++i)
    {
        const auto& fp = plan_.fields[i];
        const auto& view = req.fields[i];
        CG_DataType_t dtype = (fp.shape.elem_size == 8) ? RealDouble : RealSingle;

        const void* src = view.host_ptr;
        void* staging = pool_.data(i);

        const bool same_size = (fp.shape.elem_size == view.elem_size);
        const bool contig = is_contig(view, view.elem_size);

        if (same_size && contig)
        {
            // direct
        }
        else if (view.elem_size == 8 && fp.shape.elem_size == 4)
        {
            pack_cast_3d<float, double>(static_cast<float*>(staging),
                                        static_cast<const double*>(view.host_ptr), fp.shape.nx,
                                        fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                                        view.strides[2]);
            src = staging;
        }
        else if (view.elem_size == 4 && fp.shape.elem_size == 8)
        {
            pack_cast_3d<double, float>(static_cast<double*>(staging),
                                        static_cast<const float*>(view.host_ptr), fp.shape.nx,
                                        fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                                        view.strides[2]);
            src = staging;
        }
        else
        {
            pack_bytes_3d(static_cast<unsigned char*>(staging),
                          static_cast<const unsigned char*>(view.host_ptr), fp.shape.nx,
                          fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                          view.strides[2], fp.shape.elem_size);
            src = staging;
        }

        int fld_id{};
        cg_field_write(impl_->file, impl_->base, impl_->zone, sol_id, dtype, fp.shape.name.c_str(),
                       src, &fld_id);
    }

    // Update TimeValues array
    if (cg_goto(impl_->file, impl_->base, "BaseIterativeData_t", 1, "end"))
        return;
    int nsteps = impl_->step_index + 1;
    cg_array_write("TimeValues", RealDouble, 1, (cgsize_t*) &nsteps, impl_->times.data());
    impl_->step_index++;
#else
    (void) req;
#endif
}

void CGNSWriter::close()
{
#ifdef SOLVERLES_WITH_CGNS
    if (impl_ && impl_->opened)
    {
        cg_close(impl_->file);
        impl_->opened = false;
    }
#endif
}

} // namespace core::master::io
