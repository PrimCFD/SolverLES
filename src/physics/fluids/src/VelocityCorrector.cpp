#include "VelocityCorrector.hpp"
#include "master/FieldCatalog.hpp"
#include <cstdlib>
#include <stdexcept>
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
    info_.phases = plugin::Phase::PostBC; // after BCs so ghosts are valid before we finalize
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

    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = (vp.extents[0] + vu.extents[0]) / 2;
    (void) nx; // silence warnings
    const int ng = (vu.extents[0] - (tile.box.hi[0] - tile.box.lo[0])) / 2;

    std::vector<double> dpx((std::size_t) nx_tot * ny_tot * nz_tot, 0.0),
        dpy((std::size_t) nx_tot * ny_tot * nz_tot, 0.0),
        dpz((std::size_t) nx_tot * ny_tot * nz_tot, 0.0);

    gradp_c(static_cast<const double*>(vp.host_ptr), nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_,
            dpx.data(), dpy.data(), dpz.data());

    if (fields.contains("rho"))
    {
        auto vr = fields.view("rho");
        const double* rho = static_cast<const double*>(vr.host_ptr);
        correct_velocity_varrho_c(static_cast<double*>(vu.host_ptr),
                                  static_cast<double*>(vv.host_ptr),
                                  static_cast<double*>(vw.host_ptr), dpx.data(), dpy.data(),
                                  dpz.data(), nx_tot, ny_tot, nz_tot, ng, rho, dt);
    }
    else
    {
        correct_velocity_c(static_cast<double*>(vu.host_ptr), static_cast<double*>(vv.host_ptr),
                           static_cast<double*>(vw.host_ptr), dpx.data(), dpy.data(), dpz.data(),
                           nx_tot, ny_tot, nz_tot, ng, rho_, dt);
    }
}

} // namespace fluids
