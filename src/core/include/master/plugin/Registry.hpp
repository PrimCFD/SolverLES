#pragma once
#include "master/RunContext.hpp"
#include "master/plugin/Program.hpp"
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

/**
 * @file Registry.hpp
 * @brief Runtime factory registry for programs, actions, and globals.
 *
 * @details
 * Shared libraries register factories under string keys using the exported function
 * :cpp:func:`physics_register_v1`. The application resolves keys (e.g., from YAML)
 * and constructs the selected components.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   extern "C" bool physics_register_v1(Registry* R) {
 *     R->add_program("explicit_rk3", [](const KV&, const RunContext&){ ... });
 *     R->add_action ("compute_residual", [](const KV&, const RunContext&){ ... });
 *     return true;
 *   }
 * @endrst
 */

namespace core::master::plugin
{

class Registry
{
  public:
    using CreateProgram =
        std::function<std::unique_ptr<IProgram>(const KV&, const core::master::RunContext&)>;
    using CreateAction =
        std::function<std::shared_ptr<IAction>(const KV&, const core::master::RunContext&)>;
    using CreateGlobal =
        std::function<std::shared_ptr<IGlobal>(const KV&, const core::master::RunContext&)>;

    void add_program(std::string key, CreateProgram f) { programs_[std::move(key)] = std::move(f); }
    void add_action(std::string key, CreateAction f) { actions_[std::move(key)] = std::move(f); }
    void add_global(std::string key, CreateGlobal f) { globals_[std::move(key)] = std::move(f); }

    std::unique_ptr<IProgram> make_program(const std::string& key, const KV&,
                                           const core::master::RunContext&) const;
    std::shared_ptr<IAction> make_action(const std::string& key, const KV&,
                                         const core::master::RunContext&) const;
    std::shared_ptr<IGlobal> make_global(const std::string& key, const KV&,
                                         const core::master::RunContext&) const;

  private:
    std::unordered_map<std::string, CreateProgram> programs_;
    std::unordered_map<std::string, CreateAction> actions_;
    std::unordered_map<std::string, CreateGlobal> globals_;
};

// Plugin entry point
using RegisterFn = bool (*)(Registry*);
inline constexpr const char* kRegisterSymbol = "physics_register_v1";

} // namespace core::master::plugin
