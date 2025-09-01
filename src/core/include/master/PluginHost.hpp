#pragma once
#include "master/plugin/Registry.hpp"
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

/**
 * @file PluginHost.hpp
 * @brief Loader for shared libraries and registry owner.
 *
 * @details
 * The host owns OS handles (``dlopen``/``LoadLibrary``) and exposes helpers to build
 * programs/actions/globals from the registered factories. A builtin **noop** program is
 * installed so the solver can run without external plugins.
 *
 * @warning On Linux the application must link with ``dl`` to resolve ELF loader calls.
 */

namespace core::master
{

struct RunContext;

class PluginHost
{
  public:
    PluginHost();
    ~PluginHost();
    PluginHost(const PluginHost&) = delete;
    PluginHost& operator=(const PluginHost&) = delete;
    PluginHost(PluginHost&&) noexcept;
    PluginHost& operator=(PluginHost&&) noexcept;

    void load_library(const std::filesystem::path& lib);
    std::unique_ptr<plugin::IProgram> make_program(const std::string& key, const plugin::KV& cfg,
                                                   const RunContext& rc) const;

    // Convenience: make individual actions/globals if the app wants to assemble its own program
    std::shared_ptr<plugin::IAction> make_action(const std::string& key, const plugin::KV& cfg,
                                                 const RunContext& rc) const;
    std::shared_ptr<plugin::IGlobal> make_global(const std::string& key, const plugin::KV& cfg,
                                                 const RunContext& rc) const;

    const plugin::Registry& registry() const noexcept { return reg_; }

  private:
    plugin::Registry reg_;
    std::vector<void*> handles_;
};

} // namespace core::master
