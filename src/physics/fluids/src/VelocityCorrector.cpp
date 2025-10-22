#include "VelocityCorrector.hpp"
#include "master/FieldCatalog.hpp"
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include "kernels_fluids.h"
#include <algorithm>
#include <iostream>
#include <limits>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

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

    // Use p (centers) to define ng; face arrays have their own totals
    const int nx_i = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (vp.extents[0] - nx_i) / 2;

    // Face-gradients sized to face arrays
    std::vector<double> dpx_u((std::size_t) vu.extents[0] * vu.extents[1] * vu.extents[2], 0.0);
    std::vector<double> dpy_v((std::size_t) vv.extents[0] * vv.extents[1] * vv.extents[2], 0.0);
    std::vector<double> dpz_w((std::size_t) vw.extents[0] * vw.extents[1] * vw.extents[2], 0.0);

    // Compute pressure gradients on faces (MAC)
    gradp_faces_c(static_cast<const double*>(vp.host_ptr), vp.extents[0], vp.extents[1],
                  vp.extents[2], ng, dx_, dy_, dz_, dpx_u.data(), vu.extents[0], vu.extents[1],
                  vu.extents[2], dpy_v.data(), vv.extents[0], vv.extents[1], vv.extents[2],
                  dpz_w.data(), vw.extents[0], vw.extents[1], vw.extents[2]);

    // Velocity correction on faces (variable-ρ if present; else constant ρ_)
    if (fields.contains("rho"))
    {
        auto vr = fields.view("rho");
        correct_velocity_varrho_mac_c(
            static_cast<double*>(vu.host_ptr), static_cast<double*>(vv.host_ptr),
            static_cast<double*>(vw.host_ptr), dpx_u.data(), dpy_v.data(), dpz_w.data(),
            vu.extents[0], vu.extents[1], vu.extents[2], vv.extents[0], vv.extents[1],
            vv.extents[2], vw.extents[0], vw.extents[1], vw.extents[2], vp.extents[0],
            vp.extents[1], vp.extents[2], ng, static_cast<const double*>(vr.host_ptr), dt);
    }
    else
    {
        correct_velocity_mac_c(static_cast<double*>(vu.host_ptr), static_cast<double*>(vv.host_ptr),
                               static_cast<double*>(vw.host_ptr), dpx_u.data(), dpy_v.data(),
                               dpz_w.data(), vu.extents[0], vu.extents[1], vu.extents[2],
                               vv.extents[0], vv.extents[1], vv.extents[2], vw.extents[0],
                               vw.extents[1], vw.extents[2], vp.extents[0], vp.extents[1],
                               vp.extents[2], ng, rho_, dt);
    }
}

} // namespace fluids
