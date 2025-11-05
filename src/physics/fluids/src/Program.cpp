#include "Program.hpp"
#include "ApplyBCs.hpp"
#include "ProjectionLoop.hpp"
#include "master/FieldCatalog.hpp"
#include "master/plugin/Action.hpp"
#include <algorithm>
#include <cstdlib>
#include "kernels_fluids.h"

using namespace core::master;
using namespace core::master::plugin;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return (end && *end == 0) ? v : dflt;
}
static inline int to_i(const std::string& s, int dflt)
{
    char* end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    return (end && *end == 0) ? int(v) : dflt;
}

static std::string lower(std::string x)
{
    std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c) { return std::tolower(c); });
    return x;
}

static void maybe_parse_bc(const std::string& key, const std::string& val, Params& p)
{
    if (key.rfind("bc.", 0) != 0)
        return;

    // bc.<field>.<face>   (e.g., "bc.u.north")
    const auto dot1 = key.find('.', 3);
    if (dot1 == std::string::npos)
        return;
    const std::string field = key.substr(3, dot1 - 3); // u|v|w|p
    const std::string face = key.substr(dot1 + 1);     // west|east|south|north|bottom|top|*

    auto colon = val.find(':');
    if (colon == std::string::npos)
        return;
    const std::string type_s = lower(val.substr(0, colon));
    const std::string v_s = val.substr(colon + 1);

    BcSpec::Type t;
    if (type_s == "dirichlet")
        t = BcSpec::Type::dirichlet;
    else if (type_s == "neumann")
        t = BcSpec::Type::neumann;
    else if (type_s == "mirror" || type_s == "symmetry")
        t = BcSpec::Type::mirror;
    else if (type_s == "extrap" || type_s == "extrapolate")
        t = BcSpec::Type::extrap;
    else
        return;

    const double v = to_d(v_s, 0.0);

    auto apply = [&](const char* f) { p.bcs[field + std::string(".") + f] = BcSpec{t, v}; };

    if (face == "*")
    {
        apply("west");
        apply("east");
        apply("south");
        apply("north");
        apply("bottom");
        apply("top");
    }
    else
    {
        apply(face.c_str());
    }
}

Params parse_params(const KV& kv)
{
    Params p;
    auto get = [&](const char* k) -> const std::string*
    {
        auto it = kv.find(k);
        return (it == kv.end() ? nullptr : &it->second);
    };

    if (auto s = get("dx"))
        p.dx = to_d(*s, p.dx);
    if (auto s = get("dy"))
        p.dy = to_d(*s, p.dy);
    if (auto s = get("dz"))
        p.dz = to_d(*s, p.dz);

    if (auto s = get("rho"))
        p.rho = to_d(*s, p.rho);
    if (auto s = get("nu"))
        p.nu = to_d(*s, p.nu);
    if (auto s = get("Cs"))
        p.Cs = to_d(*s, p.Cs);

    if (auto s = get("adv_blend"))
        p.adv_blend = to_d(*s, p.adv_blend);
    if (auto s = get("cfl"))
        p.cfl = to_d(*s, p.cfl);

    for (const auto& [k, v] : kv)
        maybe_parse_bc(k, v, p);
    return p;
}

extern std::shared_ptr<IAction> make_apply_bcs(const KV&);
extern std::shared_ptr<IAction> make_init_tg(const KV&, const RunContext&);
extern std::shared_ptr<IAction> make_sgs(const KV&);
extern std::shared_ptr<IAction> make_predictor(const KV&, const RunContext&);
extern std::shared_ptr<IAction> make_poisson(const KV&, const RunContext&);
extern std::shared_ptr<IAction> make_corrector(const KV&);

LESProgram::LESProgram(const KV& kv, const RunContext& rc) : p_(parse_params(kv))
{
    auto it = kv.find("init");
    std::string init =
        (it == kv.end() ? "none" : lower(it->second)); // 'lower' already defined above

    if (init == "taylor_green" || init == "taylor-green" || init == "tg" || init == "tgv")
    {
        init_ = make_init_tg(kv, rc);
    }
    else
    {
        init_.reset(); // no initialization action
    }

    sgs_ = make_sgs(kv);
    pred_ = make_predictor(kv, rc);
    psolve_ = make_poisson(kv, rc);
    corr_ = make_corrector(kv);
    bc_ = make_apply_bcs(kv);

    // Build one loop action that can do FE or BE based on kv.
    loop_ = make_projection_loop(kv, rc, sgs_, bc_, pred_, psolve_, corr_);
}

StepPlan LESProgram::plan_step(double)
{
    StepPlan plan;
    if (!did_init_ && init_)
    { // run only once
        plan.tiled.push_back(init_);
        did_init_ = true;
    }

    plan.tiled.push_back(loop_);

    return plan;
}

LESProgram::~LESProgram()
{
    // Make sure no actions are executing anymore.
    fluids_kernels_free_scratch();
}

} // namespace fluids
