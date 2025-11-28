#include "master/io/XdmfHdf5Writer.hpp"
#include "master/Views.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"
#include "memory/MpiBox.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

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
    int rank = 0, size = 1;
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

    // Establish MPI rank/size and open HDF5 collectively
    MPI_Comm comm = mpi_unbox(cfg_.mpi_cart_comm);
    if (comm == MPI_COMM_NULL)
        comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &impl_->rank);
    MPI_Comm_size(comm, &impl_->size);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
    impl_->file = H5Fcreate(impl_->h5_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (impl_->file < 0)
        return;
    impl_->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(impl_->xfer_plist, H5FD_MPIO_COLLECTIVE);

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
        // Use Mesh *global* cell counts
        impl_->nx = cfg_.mesh->global[0];
        impl_->ny = cfg_.mesh->global[1];
        impl_->nz = cfg_.mesh->global[2];
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

        // Global dataset dims (account for face staggering)
        // Detect face alignment from LOCAL sizes (cfg_.mesh->local), not global.
        const bool iFace =
            (fp.shape.nx == cfg_.mesh->local[0] + 1 && fp.shape.ny == cfg_.mesh->local[1] &&
             fp.shape.nz == cfg_.mesh->local[2]);
        const bool jFace =
            (fp.shape.nx == cfg_.mesh->local[0] && fp.shape.ny == cfg_.mesh->local[1] + 1 &&
             fp.shape.nz == cfg_.mesh->local[2]);
        const bool kFace =
            (fp.shape.nx == cfg_.mesh->local[0] && fp.shape.ny == cfg_.mesh->local[1] &&
             fp.shape.nz == cfg_.mesh->local[2] + 1);
        const bool isFace = (iFace || jFace || kFace);
        const bool isVel = (fp.shape.name == "u" || fp.shape.name == "v" || fp.shape.name == "w");
        const hsize_t GX = (hsize_t) (impl_->nx + (iFace ? 1 : 0));
        const hsize_t GY = (hsize_t) (impl_->ny + (jFace ? 1 : 0));
        const hsize_t GZ = (hsize_t) (impl_->nz + (kFace ? 1 : 0));

        // Write native dataset unless it is a face-centered velocity (skip like CGNS)
        if (!(isVel && isFace))
        {
            // Global file space (Z,Y,X)
            const hsize_t file_dims[3] = {GZ, GY, GX};
            hid_t space_f = H5Screate_simple(3, file_dims, nullptr);
            hid_t dset = H5Dcreate2(g, fp.shape.name.c_str(), dtype, space_f, H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);
            if (dset < 0)
            {
                H5Sclose(space_f);
                H5Gclose(g);
                return;
            }
            // Select this rank's hyperslab and write collectively
            const hsize_t ox = (hsize_t) cfg_.mesh->global_lo[0];
            const hsize_t oy = (hsize_t) cfg_.mesh->global_lo[1];
            const hsize_t oz = (hsize_t) cfg_.mesh->global_lo[2];
            const hsize_t LX = (hsize_t) fp.shape.nx;
            const hsize_t LY = (hsize_t) fp.shape.ny;
            const hsize_t LZ = (hsize_t) fp.shape.nz;
            hsize_t start[3] = {oz, oy, ox}; // Z,Y,X
            hsize_t count[3] = {LZ, LY, LX};
            H5Sselect_hyperslab(space_f, H5S_SELECT_SET, start, nullptr, count, nullptr);
            const hsize_t mem_dims[3] = {LZ, LY, LX};
            hid_t space_m = H5Screate_simple(3, mem_dims, nullptr);
            if (H5Dwrite(dset, dtype, space_m, space_f, impl_->xfer_plist, src) < 0)
            {
                H5Dclose(dset);
                H5Sclose(space_f);
                H5Sclose(space_m);
                H5Gclose(g);
                return;
            }
            H5Dclose(dset);
            H5Sclose(space_m);
            H5Sclose(space_f);
        }

        // ---- Cell-centered aliases for face-centered velocity (ParaView-friendly) ----
        const bool isU = (fp.shape.name == "u");
        const bool isV = (fp.shape.name == "v");
        const bool isW = (fp.shape.name == "w");
        if ((isU || isV || isW) && isFace)
        {
            // Average onto cell centers using the same helpers as CGNS writer
            void* tmp = pool_.data(i); // reuse staging slot
            const bool src_f64 = (view.elem_size == 8);
            const bool dst_f64 = (fp.shape.elem_size == 8);
            const int nx_loc = cfg_.mesh->local[0];
            const int ny_loc = cfg_.mesh->local[1];
            const int nz_loc = cfg_.mesh->local[2];

            if (dst_f64)
            {
                auto* d = static_cast<double*>(tmp);
                if (iFace)
                {
                    src_f64 ? avg_i_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_i_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                }
                else if (jFace)
                {
                    src_f64 ? avg_j_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_j_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                }
                else
                {
                    src_f64 ? avg_k_to_cell<double, double>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_k_to_cell<float, double>(d, view, nx_loc, ny_loc, nz_loc);
                }
            }
            else
            {
                auto* d = static_cast<float*>(tmp);
                if (iFace)
                {
                    src_f64 ? avg_i_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_i_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                }
                else if (jFace)
                {
                    src_f64 ? avg_j_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_j_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                }
                else
                {
                    src_f64 ? avg_k_to_cell<double, float>(d, view, nx_loc, ny_loc, nz_loc)
                            : avg_k_to_cell<float, float>(d, view, nx_loc, ny_loc, nz_loc);
                }
            }

            // Create /Step_xxxxxx/u_cell (global dims = cells)
            const char* alias = isU ? "u_cell" : isV ? "v_cell" : "w_cell";
            const hsize_t gXc = (hsize_t) impl_->nx, gYc = (hsize_t) impl_->ny,
                          gZc = (hsize_t) impl_->nz;
            const hsize_t file_dims_cc[3] = {gZc, gYc, gXc};
            hid_t space_f_cc = H5Screate_simple(3, file_dims_cc, nullptr);
            hid_t dset_cc =
                H5Dcreate2(g, alias, dtype, space_f_cc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (dset_cc >= 0)
            {
                // Select this rank's *cell* tile (no +1 anywhere)
                const hsize_t ozc = (hsize_t) cfg_.mesh->global_lo[2]; // Z
                const hsize_t oyc = (hsize_t) cfg_.mesh->global_lo[1]; // Y
                const hsize_t oxc = (hsize_t) cfg_.mesh->global_lo[0]; // X
                const hsize_t Lzc = (hsize_t) cfg_.mesh->local[2];     // Z
                const hsize_t Lyc = (hsize_t) cfg_.mesh->local[1];     // Y
                const hsize_t Lxc = (hsize_t) cfg_.mesh->local[0];     // X
                hsize_t start_cc[3] = {ozc, oyc, oxc};
                hsize_t count_cc[3] = {Lzc, Lyc, Lxc};
                H5Sselect_hyperslab(space_f_cc, H5S_SELECT_SET, start_cc, nullptr, count_cc,
                                    nullptr);
                const hsize_t mem_dims_cc[3] = {Lzc, Lyc, Lxc};
                hid_t space_m_cc = H5Screate_simple(3, mem_dims_cc, nullptr);
                H5Dwrite(dset_cc, dtype, space_m_cc, space_f_cc, impl_->xfer_plist, tmp);
                H5Sclose(space_m_cc);
                H5Dclose(dset_cc);
            }
            H5Sclose(space_f_cc);
        }
    } // fields loop

    H5Gclose(g);

    // rewrite XDMF index (v2 or v3) -- rank 0 only
    if (impl_->rank == 0)
    {
        const int cx = impl_->nx, cy = impl_->ny, cz = impl_->nz;
        auto version_str = (cfg_.xdmf_version == WriterConfig::XdmfVersion::V2) ? "2.0" : "3.0";
        auto topo_token =
            (cfg_.xdmf_version == WriterConfig::XdmfVersion::V2) ? "3DCORECTMESH" : "3DCoRectMesh";
        std::ofstream xmf(impl_->xmf_path, std::ios::trunc);
        xmf << R"(<?xml version="1.0" ?>)" << "\n<Xdmf Version=\"" << version_str
            << "\">\n  <Domain>\n"
            << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" "
               "CollectionType=\"Temporal\">\n";
        for (int s = 0; s < impl_->step_index; ++s)
        {
            char stepname[64];
            std::snprintf(stepname, sizeof(stepname), "Step_%06d", s);
            xmf << "      <Grid Name=\"" << stepname << "\" GridType=\"Uniform\">\n"
                << "        <Time Value=\"" << impl_->times[s] << "\"/>\n"
                << "        <Topology TopologyType=\"" << topo_token << "\" Dimensions=\""
                << (cz + 1) << " " << (cy + 1) << " " << (cx + 1) << "\"/>\n"
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
                const char* center = xmf_center_for(fp.shape.nx, fp.shape.ny, fp.shape.nz,
                                                    impl_->nx, impl_->ny, impl_->nz);
                const bool isVel = (fname == "u" || fname == "v" || fname == "w");
                // Global attribute dimensions (account for face staggering)
                const bool iFace = (fp.shape.nx == impl_->nx + 1 && fp.shape.ny == impl_->ny &&
                                    fp.shape.nz == impl_->nz);
                const bool jFace = (fp.shape.nx == impl_->nx && fp.shape.ny == impl_->ny + 1 &&
                                    fp.shape.nz == impl_->nz);
                const bool kFace = (fp.shape.nx == impl_->nx && fp.shape.ny == impl_->ny &&
                                    fp.shape.nz == impl_->nz + 1);
                const bool isFace = (iFace || jFace || kFace);
                // Skip advertising native face velocity attributes
                if (!(isVel && isFace))
                {
                    const int gnx = impl_->nx + (iFace ? 1 : 0);
                    const int gny = impl_->ny + (jFace ? 1 : 0);
                    const int gnz = impl_->nz + (kFace ? 1 : 0);
                    xmf << " <Attribute Name=\"" << fname << "\" AttributeType=\"Scalar\" Center=\""
                        << center << "\">\n"
                        << " <DataItem Dimensions=\"" << gnz << " " << gny << " " << gnx
                        << "\" NumberType=\"Float\" Precision=\"" << fp.shape.elem_size
                        << "\" Format=\"HDF\">" << fs::path(impl_->h5_path).filename().string()
                        << ":/" << stepname << "/" << fname << "</DataItem>\n"
                        << " </Attribute>\n";
                }

                // If velocity was face-centered, also advertise the cell-centered alias
                if (fname == "u" || fname == "v" || fname == "w")
                {
                    const std::string alias = (fname == "u")   ? "u_cell"
                                              : (fname == "v") ? "v_cell"
                                                               : "w_cell";
                    xmf << " <Attribute Name=\"" << alias
                        << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                        << " <DataItem Dimensions=\"" << impl_->nz << " " << impl_->ny << " "
                        << impl_->nx << "\" NumberType=\"Float\" Precision=\"" << fp.shape.elem_size
                        << "\" Format=\"HDF\">" << fs::path(impl_->h5_path).filename().string()
                        << ":/" << stepname << "/" << alias << "</DataItem>\n"
                        << " </Attribute>\n";
                }
            }
            xmf << "      </Grid>\n";
        }
        xmf << "    </Grid>\n  </Domain>\n</Xdmf>\n";
        xmf.close();
    } // rank 0
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
