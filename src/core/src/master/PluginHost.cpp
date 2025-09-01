#include "master/PluginHost.hpp"
#include "master/RunContext.hpp"
#include <stdexcept>
#include <utility>

/// \cond DOXYGEN_EXCLUDE

#if defined(_WIN32)
#include <windows.h>

static void* load_so(const std::string& p)
{
    return (void*) ::LoadLibraryA(p.c_str());
}
static void close_so(void* h)
{
    if (h)
        ::FreeLibrary((HMODULE) h);
}
static void* load_sym(void* h, const char* s)
{
    return (void*) ::GetProcAddress((HMODULE) h, s);
}

#else
#include <dlfcn.h>

static void* load_so(const std::string& p)
{
    return ::dlopen(p.c_str(), RTLD_NOW);
}
static void close_so(void* h)
{
    if (h)
        ::dlclose(h);
}
static void* load_sym(void* h, const char* s)
{
    return ::dlsym(h, s);
}

#endif

/// \endcond

using namespace core::master;

// ----- builtin "noop" program so tests/app can run with no plugins -----------
namespace
{

using namespace core::master::plugin;

struct NoopProgram final : IProgram
{
    StepPlan plan_step(double) override { return {}; } // no tiled actions, no globals
};

inline void register_builtin_noop(plugin::Registry& r)
{
    r.add_program("noop",
                  [](const KV&, const RunContext&) { return std::make_unique<NoopProgram>(); });
}

} // namespace
// -----------------------------------------------------------------------------

PluginHost::PluginHost()
{
    register_builtin_noop(reg_);
}

PluginHost::~PluginHost()
{
    for (void* h : handles_)
        close_so(h);
}

PluginHost::PluginHost(PluginHost&& o) noexcept
    : reg_(std::move(o.reg_)), handles_(std::move(o.handles_))
{
    o.handles_.clear();
}

PluginHost& PluginHost::operator=(PluginHost&& o) noexcept
{
    if (this != &o)
    {
        for (void* h : handles_)
            close_so(h);
        reg_ = std::move(o.reg_);
        handles_ = std::move(o.handles_);
        o.handles_.clear();
    }
    return *this;
}

void PluginHost::load_library(const std::filesystem::path& lib)
{
    auto* h = load_so(lib.string());
    if (!h)
        throw std::runtime_error("Failed to load plugin library: " + lib.string());
    handles_.push_back(h);

    auto* sym = load_sym(h, plugin::kRegisterSymbol);
    if (!sym)
        throw std::runtime_error("Missing symbol in plugin: " +
                                 std::string(plugin::kRegisterSymbol));
    auto reg_fn = reinterpret_cast<plugin::RegisterFn>(sym);
    if (!reg_fn(&reg_))
        throw std::runtime_error("Plugin registration returned failure");
}

std::unique_ptr<plugin::IProgram>
PluginHost::make_program(const std::string& key, const plugin::KV& cfg, const RunContext& rc) const
{
    return reg_.make_program(key, cfg, rc);
}

std::shared_ptr<plugin::IAction>
PluginHost::make_action(const std::string& key, const plugin::KV& cfg, const RunContext& rc) const
{
    return reg_.make_action(key, cfg, rc);
}

std::shared_ptr<plugin::IGlobal>
PluginHost::make_global(const std::string& key, const plugin::KV& cfg, const RunContext& rc) const
{
    return reg_.make_global(key, cfg, rc);
}
