#include "InitTG.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include "mesh/Mesh.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <mpi.h>
#include <iostream>

#include "memory/MpiBox.hpp"
#include "mesh/Field.hpp"
#include "master/Views.hpp"

using namespace core::master;

namespace fluids
{

InitTG::InitTG(double Lx, double Ly, double Lz, double U0, void* mpi_comm)
    : Lx_(Lx), Ly_(Ly), Lz_(Lz), U0_(U0), mpi_comm_(mpi_comm)
{
    info_.name = "init_taylor_green";
    info_.phases = plugin::Phase::PreExchange; // run before any halos (initial fill)
}

std::shared_ptr<plugin::IAction> make_init_tg(const plugin::KV& kv,
                                              const core::master::RunContext& rc)
{
    // U0 default 1.0; Lx/Ly/Lz default 1.0 if not provided
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    auto to_d = [](const std::string& s, double d)
    {
        char* e = nullptr;
        double v = std::strtod(s.c_str(), &e);
        return (e && *e == 0) ? v : d;
    };
    const double Lx = to_d(get("Lx", "1.0"), 1.0);
    const double Ly = to_d(get("Ly", "1.0"), 1.0);
    const double Lz = to_d(get("Lz", "1.0"), 1.0);
    const double U0 = to_d(get("U0", "1.0"), 1.0);
    return std::make_shared<InitTG>(Lx, Ly, Lz, U0, rc.mpi_comm);
}

void InitTG::execute(const MeshTileView& tile, FieldCatalog& fields, double)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.init_tg] fields u/v/w must be registered by the app.");

    // Views (u: IFace, v: JFace, w: KFace on a MAC grid)
    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");

    // ---- Single Source of Truth (mesh) for halos and interior cell counts ----
    if (!tile.mesh)
        throw std::runtime_error("[fluids.init_tg] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng; // authoritative ghost/halo width
    // Local interior (centers)
    const int nx_c = mesh.local[0], ny_c = mesh.local[1], nz_c = mesh.local[2];
    // Global interior (centers)
    const int NX = mesh.global[0], NY = mesh.global[1], NZ = mesh.global[2];
    // Global starting cell-index (center-based) of this slab (provided by Mesh)
    const int i0 = mesh.global_lo[0];
    const int j0 = mesh.global_lo[1];
    const int k0 = mesh.global_lo[2];

    // Totals implied by MAC staggering (+1 along face-normal)
    const int nxc_tot = nx_c + 2 * ng, nyc_tot = ny_c + 2 * ng, nzc_tot = nz_c + 2 * ng;
    const int nxu_tot = (nx_c + 1) + 2 * ng, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = (ny_c + 1) + 2 * ng, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = (nz_c + 1) + 2 * ng;

    // Optional sanity: allocated views must match mesh totals
    auto check = [](const char* name, const std::array<int, 3>& e, int ex, int ey, int ez)
    {
        if (e[0] != ex || e[1] != ey || e[2] != ez)
            throw std::runtime_error(std::string("[fluids.init_tg] view '") + name +
                                     "' extents do not match mesh totals.");
    };
    check("u", vu.extents, nxu_tot, nyu_tot, nzu_tot);
    check("v", vv.extents, nxv_tot, nyv_tot, nzv_tot);
    check("w", vw.extents, nxw_tot, nyw_tot, nzw_tot);

    // Linear indexers per field (row-major: x fastest)
    auto idxU = [nxu_tot, nyu_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxu_tot) * (std::size_t(j) + std::size_t(nyu_tot) * std::size_t(k));
    };
    auto idxV = [nxv_tot, nyv_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxv_tot) * (std::size_t(j) + std::size_t(nyv_tot) * std::size_t(k));
    };
    auto idxW = [nxw_tot, nyw_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxw_tot) * (std::size_t(j) + std::size_t(nyw_tot) * std::size_t(k));
    };

    double* u = static_cast<double*>(vu.host_ptr);
    double* v = static_cast<double*>(vv.host_ptr);
    double* w = static_cast<double*>(vw.host_ptr);

    // Physical spacings: use GLOBAL counts so all ranks sample the same function
    const double dx = Lx_ / static_cast<double>(NX);
    const double dy = Ly_ / static_cast<double>(NY);
    const double dz = Lz_ / static_cast<double>(NZ);

    // ---- U on IFaces: x at faces, y/z at centers ----
    const int nx_u = nx_c + 1; // +1 along normal
    for (int k = 0; k < nz_c; ++k)
    {
        const double zc = (k0 + k + 0.5) * dz; // GLOBAL z center
        for (int j = 0; j < ny_c; ++j)
        {
            const double yc = (j0 + j + 0.5) * dy; // GLOBAL y center
            for (int i = 0; i < nx_u; ++i)
            {                                    // faces in x
                const double xf = (i0 + i) * dx; // GLOBAL x face
                const int I = i + ng, J = j + ng, K = k + ng;
                const double sX = std::sin(2 * M_PI * xf / Lx_);
                const double cY = std::cos(2 * M_PI * yc / Ly_);
                const double cZ = std::cos(2 * M_PI * zc / Lz_);
                u[idxU(I, J, K)] = U0_ * sX * cY * cZ;
            }
        }
    }

    // ---- V on JFaces: y at faces, x/z at centers ----
    const int ny_v = ny_c + 1; // +1 along normal
    for (int k = 0; k < nz_c; ++k)
    {
        const double zc = (k0 + k + 0.5) * dz; // GLOBAL z center
        for (int j = 0; j < ny_v; ++j)
        {                                    // faces in y
            const double yf = (j0 + j) * dy; // GLOBAL y face
            for (int i = 0; i < nx_c; ++i)
            {
                const double xc = (i0 + i + 0.5) * dx; // GLOBAL x center
                const int I = i + ng, J = j + ng, K = k + ng;
                const double cX = std::cos(2 * M_PI * xc / Lx_);
                const double sY = std::sin(2 * M_PI * yf / Ly_);
                const double cZ = std::cos(2 * M_PI * zc / Lz_);
                v[idxV(I, J, K)] = -U0_ * cX * sY * cZ;
            }
        }
    }

    // ---- W on KFaces: z at faces, x/y at centers ----
    const int nz_w = nz_c + 1; // +1 along normal
    for (int k = 0; k < nz_w; ++k)
    {                                    // faces in z
        const double zf = (k0 + k) * dz; // GLOBAL z face
        for (int j = 0; j < ny_c; ++j)
        {
            const double yc = (j0 + j + 0.5) * dy; // GLOBAL y center
            for (int i = 0; i < nx_c; ++i)
            {
                const double xc = (i0 + i + 0.5) * dx; // GLOBAL x center
                const int I = i + ng, J = j + ng, K = k + ng;
                // Standard 3D TG often has w=0; keep it zero unless you use a symmetric variant.
                (void) zf;
                (void) xc;
                (void) yc;
                w[idxW(I, J, K)] = 0.0;
            }
        }
    }

    // Pressure = 0 by default (app initializes arrays to zero).

    // One-time halo fill so ghosts are valid at t=0 on multi-rank runs.
    if (tile.mesh)
    {
        core::master::exchange_named_fields(fields, *tile.mesh, mpi_comm_, {"u", "v", "w"});
    }

}

} // namespace fluids
