#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fluids
{

enum class Face
{
    west,
    east,
    south,
    north,
    bottom,
    top
};

struct BcSpec
{
    enum class Type
    {
        dirichlet,
        neumann,
        mirror,
        extrap
    };
    Type type{Type::dirichlet};
    double value{0.0};
};

using BcTable = std::unordered_map<std::string, BcSpec>;

struct Params
{
    double dx{1.0}, dy{1.0}, dz{1.0};
    double rho{1.0}, nu{1e-3}, Cs{0.16};
    double alpha_u{0.7}, alpha_p{0.3};
    double adv_blend{0.0};
    double cfl{0.7};
    BcTable bcs;
};

Params parse_params(const core::master::plugin::KV& kv);

std::shared_ptr<core::master::plugin::IAction> make_apply_bcs(const core::master::plugin::KV& kv);
std::shared_ptr<core::master::plugin::IAction> make_init_tg(const core::master::plugin::KV&);
std::shared_ptr<core::master::plugin::IAction> make_sgs(const core::master::plugin::KV&);
std::shared_ptr<core::master::plugin::IAction> make_predictor(const core::master::plugin::KV&,
                                                              const core::master::RunContext&);
std::shared_ptr<core::master::plugin::IAction> make_poisson(const core::master::plugin::KV&,
                                                            const core::master::RunContext&);
std::shared_ptr<core::master::plugin::IAction> make_corrector(const core::master::plugin::KV&);

struct LESProgram final : core::master::plugin::IProgram
{
    explicit LESProgram(const core::master::plugin::KV& kv, const core::master::RunContext& rc);
    core::master::plugin::StepPlan plan_step(double dt) override;
    bool did_init_ = false;
    ~LESProgram();

  private:
    Params p_;
    bool first_{true};
    std::shared_ptr<core::master::plugin::IAction> init_, sgs_, pred_, psolve_, corr_, bc_;
    std::shared_ptr<core::master::plugin::IAction> loop_;
};

} // namespace fluids
