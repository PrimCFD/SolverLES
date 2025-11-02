#include "master/io/XdmfHdf5Writer.hpp"
#include "master/Views.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <hdf5.h>

namespace core::master::io
{
namespace fs = std::filesystem;

/// \cond DOXYGEN_EXCLUDE

struct XdmfHdf5Writer::Impl
{
    hid_t file = -1;
    hid_t xfer_plist = -1;

    std::string h5_path;
    std::string xmf_path;

    bool opened = false;
    bool planned = false;
    int step_index = 0;

    int nx = 0, ny = 0, nz = 0;
    std::vector<double> times;
    std::vector<std::string> field_names;
};

static hid_t h5_dtype(std::size_t bytes)
{
    return bytes == 8 ? H5T_IEEE_F64LE : H5T_IEEE_F32LE;
}

static inline bool is_contig(const AnyFieldView& v, std::size_t elem_size)
{
    const int nx = v.extents[0], ny = v.extents[1], nz = v.extents[2];
    const auto sx = v.strides[0], sy = v.strides[1], sz = v.strides[2];
    return sx == (std::ptrdiff_t) elem_size && sy == (std::ptrdiff_t)(elem_size * nx) &&
           sz == (std::ptrdiff_t)(elem_size * nx * ny);
}

static inline const char* xmf_center_for(int fnx, int fny, int fnz, int nx, int ny, int nz)
{
    if (fnx == nx && fny == ny && fnz == nz)
        return "Cell";
    // XDMF supports "Face" center for staggered face data on structured grids
    if ((fnx == nx + 1 && fny == ny && fnz == nz) || (fnx == nx && fny == ny + 1 && fnz == nz) ||
        (fnx == nx && fny == ny && fnz == nz + 1))
        return "Face";
    // Conservative fallback
    return "Node";
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

XdmfHdf5Writer::XdmfHdf5Writer(WriterConfig cfg) : cfg_(cfg), impl_(new Impl) {}
XdmfHdf5Writer::~XdmfHdf5Writer()
{
    close();
}

void XdmfHdf5Writer::open_case(const std::string& case_name)
{
    fs::create_directories(cfg_.path);
    impl_->h5_path = cfg_.path + "/" + case_name + ".h5";
    impl_->xmf_path = cfg_.path + "/" + case_name + ".xmf";

    impl_->file = H5Fcreate(impl_->h5_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (impl_->file < 0)
        return;
    impl_->xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    impl_->opened = true;
    impl_->planned = false;
    impl_->step_index = 0;
    impl_->times.clear();
}

void XdmfHdf5Writer::write(const WriteRequest& req)
{
    if (!impl_->opened)
        return;

    if (!impl_->planned)
    {
        plan_ = make_plan_from_request(req.fields, cfg_.precision);
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
        impl_->field_names.clear();
        for (auto& f : plan_.fields)
            impl_->field_names.push_back(f.shape.name);
        std::size_t maxb = 0;
        for (auto& f : plan_.fields)
            maxb = std::max(maxb, f.bytes);
        pool_.reserve(plan_.fields.size(), maxb);
        impl_->planned = true;
    }

    impl_->times.push_back(req.time);
    const int step_idx = impl_->step_index++;

    // group /Step_xxxxxx
    char gname[64];
    std::snprintf(gname, sizeof(gname), "Step_%06d", step_idx);
    hid_t g = H5Gcreate2(impl_->file, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0)
        return;

    for (std::size_t i = 0; i < plan_.fields.size(); ++i)
    {
        const auto& fp = plan_.fields[i];
        const auto& view = req.fields[i];
        hid_t dtype = h5_dtype(fp.shape.elem_size);

        const void* src = view.host_ptr;
        void* staging = pool_.data(i);

        const bool same_size = (fp.shape.elem_size == view.elem_size);
        const bool contig = is_contig(view, view.elem_size);

        if (same_size && contig)
        {
            // direct
        }
        else if (view.elem_size == 8 && fp.shape.elem_size == 4)
        { // Src=double -> Dst=float
            pack_cast_3d<double, float>(static_cast<float*>(staging),
                                        static_cast<const double*>(view.host_ptr), fp.shape.nx,
                                        fp.shape.ny, fp.shape.nz, view.strides[0], view.strides[1],
                                        view.strides[2]);
            src = staging;
        }
        else if (view.elem_size == 4 && fp.shape.elem_size == 8)
        { // Src=float  -> Dst=double
            pack_cast_3d<float, double>(static_cast<double*>(staging),
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

        const hsize_t dsz[3] = {(hsize_t) fp.shape.nz, (hsize_t) fp.shape.ny,
                                (hsize_t) fp.shape.nx};
        hid_t space_f = H5Screate_simple(3, dsz, nullptr);
        hid_t dset = H5Dcreate2(g, fp.shape.name.c_str(), dtype, space_f, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);

        if (dset < 0)
        {
            H5Sclose(space_f);
            H5Gclose(g);
            return;
        }
        if (H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, impl_->xfer_plist, src) < 0)
        {
            H5Dclose(dset);
            H5Sclose(space_f);
            H5Gclose(g);
            return;
        }
        H5Dclose(dset);
        H5Sclose(space_f);
    }

    H5Gclose(g);

    // rewrite XDMF index (v2 or v3)
    const int cx = impl_->nx, cy = impl_->ny, cz = impl_->nz;
    auto version_str = (cfg_.xdmf_version == WriterConfig::XdmfVersion::V2) ? "2.0" : "3.0";
    auto topo_token =
        (cfg_.xdmf_version == WriterConfig::XdmfVersion::V2) ? "3DCORECTMESH" : "3DCoRectMesh";
    std::ofstream xmf(impl_->xmf_path, std::ios::trunc);
    xmf << R"(<?xml version="1.0" ?>)" << "\n<Xdmf Version=\"" << version_str << "\">\n  <Domain>\n"
        << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
    for (int s = 0; s < impl_->step_index; ++s)
    {
        char stepname[64];
        std::snprintf(stepname, sizeof(stepname), "Step_%06d", s);
        xmf << "      <Grid Name=\"" << stepname << "\" GridType=\"Uniform\">\n"
            << "        <Time Value=\"" << impl_->times[s] << "\"/>\n"
            << "        <Topology TopologyType=\"" << topo_token << "\" Dimensions=\"" << (cz + 1)
            << " " << (cy + 1) << " " << (cx + 1) << "\"/>\n"
            << "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
            << "          <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" "
               "Format=\"XML\">0 0 0</DataItem>\n"
            << "          <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" "
               "Format=\"XML\">1 1 1</DataItem>\n"
            << "        </Geometry>\n";
        for (std::size_t i = 0; i < plan_.fields.size(); ++i)
        {
            const auto& fp = plan_.fields[i];
            const auto& fname = fp.shape.name;
            const char* center = xmf_center_for(fp.shape.nx, fp.shape.ny, fp.shape.nz, impl_->nx,
                                                impl_->ny, impl_->nz);
            xmf << " <Attribute Name=\"" << fname << "\" AttributeType=\"Scalar\" Center=\""
                << center << "\">\n"
                << " <DataItem Dimensions=\"" << fp.shape.nz << " " << fp.shape.ny << " "
                << fp.shape.nx << "\" NumberType=\"Float\" Precision=\"" << fp.shape.elem_size
                << "\" Format=\"HDF\">" << fs::path(impl_->h5_path).filename().string() << ":/"
                << stepname << "/" << fname << "</DataItem>\n"
                << " </Attribute>\n";
        }
        xmf << "      </Grid>\n";
    }
    xmf << "    </Grid>\n  </Domain>\n</Xdmf>\n";
    xmf.close();
}

void XdmfHdf5Writer::close()
{
    if (impl_)
    {
        if (impl_->xfer_plist >= 0)
        {
            H5Pclose(impl_->xfer_plist);
            impl_->xfer_plist = -1;
        }
        if (impl_->file >= 0)
        {
            H5Fclose(impl_->file);
            impl_->file = -1;
        }
        impl_->opened = false;
    }
}

} // namespace core::master::io
