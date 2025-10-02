#include "SGS.hpp"
#include "master/FieldCatalog.hpp"
#include "master/Views.hpp"
#include <cstdlib>
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

SGS::SGS(double dx, double dy, double dz, double Cs) : dx_(dx), dy_(dy), dz_(dz), Cs_(Cs)
{
    info_.name = "sgs_smagorinsky";
    info_.phases = plugin::Phase::Interior;
}

std::shared_ptr<plugin::IAction> make_sgs(const plugin::KV& kv)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    return std::make_shared<SGS>(to_d(get("dx", "1.0"), 1.0), to_d(get("dy", "1.0"), 1.0),
                                 to_d(get("dz", "1.0"), 1.0), to_d(get("Cs", "0.17"), 0.17));
}

void SGS::execute(const MeshTileView& tile, FieldCatalog& fields, double)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.sgs] fields u/v/w must be registered.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");

    const int nx_tot = vu.extents[0], ny_tot = vu.extents[1], nz_tot = vu.extents[2];
    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nx_tot - nx) / 2;

    double* u = static_cast<double*>(vu.host_ptr);
    double* v = static_cast<double*>(vv.host_ptr);
    double* w = static_cast<double*>(vw.host_ptr);

    double* nu_t = nullptr;
    wrote_to_field_ = false;
    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t");
        nu_t = static_cast<double*>(vt.host_ptr);
        wrote_to_field_ = true;
    }
    else
    {
        if (scratch_.size() != static_cast<std::size_t>(nx_tot * ny_tot * nz_tot))
            scratch_.assign(static_cast<std::size_t>(nx_tot * ny_tot * nz_tot), 0.0);
        nu_t = scratch_.data();
    }

    sgs_smagorinsky_mac_c(u, v, w, nx_tot, ny_tot, nz_tot, ng, dx_, dy_, dz_, Cs_, nu_t);
}

} // namespace fluids
