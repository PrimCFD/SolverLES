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

    if (!tile.mesh)
        throw std::runtime_error("[fluids.sgs] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;
    const int nx_i = mesh.local[0], ny_i = mesh.local[1], nz_i = mesh.local[2];
    const int nxc_tot = nx_i + 2 * ng;
    const int nyc_tot = ny_i + 2 * ng;
    const int nzc_tot = nz_i + 2 * ng;
    const int nxu_tot = (nx_i + 1) + 2 * ng, nyu_tot = nyc_tot, nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot, nyv_tot = (ny_i + 1) + 2 * ng, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot, nyw_tot = nyc_tot, nzw_tot = (nz_i + 1) + 2 * ng;

    // Sanity: views should match mesh totals
    auto check = [](const char* name, const std::array<int, 3>& e, int ex, int ey, int ez)
    {
        if (e[0] != ex || e[1] != ey || e[2] != ez)
            throw std::runtime_error(std::string("[fluids.sgs] view '") + name +
                                     "' extents do not match mesh totals.");
    };
    check("u", vu.extents, nxu_tot, nyu_tot, nzu_tot);
    check("v", vv.extents, nxv_tot, nyv_tot, nzv_tot);
    check("w", vw.extents, nxw_tot, nyw_tot, nzw_tot);
    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t");
        check("nu_t", vt.extents, nxc_tot, nyc_tot, nzc_tot);
    }
    if (fields.contains("p"))
    {
        auto vp = fields.view("p");
        check("p", vp.extents, nxc_tot, nyc_tot, nzc_tot);
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
                          /* u faces */ nxu_tot, nyu_tot, nzu_tot,
                          /* v faces */ nxv_tot, nyv_tot, nzv_tot,
                          /* w faces */ nxw_tot, nyw_tot, nzw_tot,
                          /* centers  */ nxc_tot, nyc_tot, nzc_tot, ng, dx_, dy_, dz_, Cs_, nu_t);
}

} // namespace fluids
