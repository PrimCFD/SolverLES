#include "InitTG.hpp"
#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include <cmath>
#include <stdexcept>

using namespace core::master;

namespace fluids
{

InitTG::InitTG(double Lx, double Ly, double Lz, double U0) : Lx_(Lx), Ly_(Ly), Lz_(Lz), U0_(U0)
{
    info_.name = "init_taylor_green";
    info_.phases = plugin::Phase::PreExchange; // run before any halos (initial fill)
}

std::shared_ptr<plugin::IAction> make_init_tg(const plugin::KV& kv)
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
    return std::make_shared<InitTG>(Lx, Ly, Lz, U0);
}

void InitTG::execute(const MeshTileView& tile, FieldCatalog& fields, double)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.init_tg] fields u/v/w must be registered by the app.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");

    auto nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    auto nx = tile.box.hi[0] - tile.box.lo[0];
    auto ny = tile.box.hi[1] - tile.box.lo[1];
    auto nz = tile.box.hi[2] - tile.box.lo[2];
    int ng = (nx_tot - nx) / 2;

    double* u = static_cast<double*>(vu.host_ptr);
    double* v = static_cast<double*>(vv.host_ptr);
    double* w = static_cast<double*>(vw.host_ptr);

    auto idx = [nx_tot, ny_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nx_tot) * (std::size_t(j) + std::size_t(ny_tot) * std::size_t(k));
    };

    // Fill interior with Taylor–Green vortex at t=0
    for (int k = 0; k < nz; ++k)
    {
        const double z = (k + 0.5) * (Lz_ / nz);
        for (int j = 0; j < ny; ++j)
        {
            const double y = (j + 0.5) * (Ly_ / ny);
            for (int i = 0; i < nx; ++i)
            {
                const double x = (i + 0.5) * (Lx_ / nx);
                const int I = i + ng, J = j + ng, K = k + ng;
                const double s = std::sin(2 * M_PI * x / Lx_);
                const double c = std::cos(2 * M_PI * x / Lx_);
                const double sy = std::sin(2 * M_PI * y / Ly_);
                const double cy = std::cos(2 * M_PI * y / Ly_);
                const double sz = std::sin(2 * M_PI * z / Lz_);
                const double cz = std::cos(2 * M_PI * z / Lz_);
                u[idx(I, J, K)] = U0_ * s * cy * cz;
                v[idx(I, J, K)] = -U0_ * c * sy * cz;
                w[idx(I, J, K)] = 0.0;
            }
        }
    }

    // Pressure = 0 by default (app initializes arrays to zero).
}

} // namespace fluids
