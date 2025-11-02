#include "VelocityCorrector.hpp"
#include "master/FieldCatalog.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <stdexcept>
#include <vector>
#include "kernels_fluids.h"

using namespace core::master;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* e = nullptr;
    double v = std::strtod(s.c_str(), &e);
    return (e && *e == 0) ? v : dflt;
}

Corrector::Corrector(double rho, double dx, double dy, double dz)
    : rho_(rho), dx_(dx), dy_(dy), dz_(dz)
{
    info_.name = "velocity_corrector";
    info_.phases = plugin::Phase::PostBC; // keep: ghosts valid before finalize
}

std::shared_ptr<plugin::IAction> make_corrector(const plugin::KV& kv)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    return std::make_shared<Corrector>(to_d(get("rho", "1.0"), 1.0), to_d(get("dx", "1.0"), 1.0),
                                       to_d(get("dy", "1.0"), 1.0), to_d(get("dz", "1.0"), 1.0));
}

void Corrector::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{
    (void) tile;
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w") ||
        !fields.contains("p"))
        throw std::runtime_error("[fluids.correct] u/v/w/p required.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    if (!tile.mesh)
        throw std::runtime_error("[fluids.correct] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;
    const int nx_i = mesh.local[0], ny_i = mesh.local[1], nz_i = mesh.local[2];
    const int nxc_tot = nx_i + 2 * ng, nyc_tot = ny_i + 2 * ng, nzc_tot = nz_i + 2 * ng;
    const int nxu_tot = (nx_i + 1) + 2 * ng, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = (ny_i + 1) + 2 * ng, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = (nz_i + 1) + 2 * ng;

    auto check = [](const char* name, const std::array<int, 3>& e, int ex, int ey, int ez)
    {
        if (e[0] != ex || e[1] != ey || e[2] != ez)
            throw std::runtime_error(std::string("[fluids.correct] view '") + name +
                                     "' extents do not match mesh totals.");
    };
    check("p", vp.extents, nxc_tot, nyc_tot, nzc_tot);
    check("u", vu.extents, nxu_tot, nyu_tot, nzu_tot);
    check("v", vv.extents, nxv_tot, nyv_tot, nzv_tot);
    check("w", vw.extents, nxw_tot, nyw_tot, nzw_tot);

    // Face-gradients sized to face arrays
    std::vector<double> dpx_u((std::size_t) vu.extents[0] * vu.extents[1] * vu.extents[2], 0.0);
    std::vector<double> dpy_v((std::size_t) vv.extents[0] * vv.extents[1] * vv.extents[2], 0.0);
    std::vector<double> dpz_w((std::size_t) vw.extents[0] * vw.extents[1] * vw.extents[2], 0.0);

    // Compute pressure gradients on faces (MAC)
    gradp_faces_c(static_cast<const double*>(vp.host_ptr), nxc_tot, nyc_tot, nzc_tot, ng, dx_, dy_,
                  dz_, dpx_u.data(), nxu_tot, nyu_tot, nzu_tot, dpy_v.data(), nxv_tot, nyv_tot,
                  nzv_tot, dpz_w.data(), nxw_tot, nyw_tot, nzw_tot);

    // Velocity correction on faces (variable-ρ if present; else constant ρ_)
    if (fields.contains("rho"))
    {
        auto vr = fields.view("rho");
        correct_velocity_varrho_mac_c(
            static_cast<double*>(vu.host_ptr), static_cast<double*>(vv.host_ptr),
            static_cast<double*>(vw.host_ptr), dpx_u.data(), dpy_v.data(), dpz_w.data(), nxu_tot,
            nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot, nxw_tot, nyw_tot, nzw_tot, nxc_tot,
            nyc_tot, nzc_tot, ng, static_cast<const double*>(vr.host_ptr), dt);
    }
    else
    {
        correct_velocity_mac_c(static_cast<double*>(vu.host_ptr), static_cast<double*>(vv.host_ptr),
                               static_cast<double*>(vw.host_ptr), dpx_u.data(), dpy_v.data(),
                               dpz_w.data(), nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot,
                               nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng, rho_, dt);
    }
}

} // namespace fluids
