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

    auto vu = fields.view("u"); // face-x
    auto vv = fields.view("v"); // face-y
    auto vw = fields.view("w"); // face-z

    // Get ng and center totals from a center field (prefer p, or nu_t if present)
    int ng = 0, nxc_tot = 0, nyc_tot = 0, nzc_tot = 0;

    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t"); // center field
        nxc_tot = vt.extents[0];
        nyc_tot = vt.extents[1];
        nzc_tot = vt.extents[2];
        // infer ng by comparing with tile box like you do elsewhere
        const int nx_c = tile.box.hi[0] - tile.box.lo[0];
        ng = (nxc_tot - nx_c) / 2;
    }
    else if (fields.contains("p"))
    {
        auto vp = fields.view("p"); // center field
        nxc_tot = vp.extents[0];
        nyc_tot = vp.extents[1];
        nzc_tot = vp.extents[2];
        const int nx_c = tile.box.hi[0] - tile.box.lo[0];
        ng = (nxc_tot - nx_c) / 2;
    }
    else
    {
        throw std::runtime_error(
            "[fluids.sgs] need a center-located field (p or nu_t) to infer center extents.");
    }

    double* nu_t = nullptr;
    wrote_to_field_ = false;
    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t"); // center
        nu_t = static_cast<double*>(vt.host_ptr);
        wrote_to_field_ = true;
    }
    else
    {
        const std::size_t ncent = (std::size_t) nxc_tot * nyc_tot * nzc_tot;
        if (scratch_.size() != ncent)
            scratch_.assign(ncent, 0.0);
        nu_t = scratch_.data();
    }

    sgs_smagorinsky_mac_c(static_cast<const double*>(vu.host_ptr),
                          static_cast<const double*>(vv.host_ptr),
                          static_cast<const double*>(vw.host_ptr),
                          /* u faces */ vu.extents[0], vu.extents[1], vu.extents[2],
                          /* v faces */ vv.extents[0], vv.extents[1], vv.extents[2],
                          /* w faces */ vw.extents[0], vw.extents[1], vw.extents[2],
                          /* centers  */ nxc_tot, nyc_tot, nzc_tot, ng, dx_, dy_, dz_, Cs_, nu_t);
}

} // namespace fluids
