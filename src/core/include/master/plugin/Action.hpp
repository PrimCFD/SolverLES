#pragma once
#include "master/Views.hpp"
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/**
 * @file Action.hpp
 * @brief Domain-agnostic plugin actions and phases.
 *
 * @details
 * An :cpp:class:`core::master::plugin::IAction` runs per tile and declares
 * a phase mask and data-access needs. A :cpp:class:`core::master::plugin::IGlobal`
 * runs once per step (reductions, diagnostics).
 *
 * @rst
 * **Phases**
 *
 * - ``PreExchange``: before halo posts (rare)
 * - ``Interior``: while halos are in flight (interior-only work)
 * - ``PostExchange``: after halo completion, before BCs
 * - ``PostBC``: after BCs
 * - ``EndStep``: final hooks
 * @endrst
 */

namespace core::master
{
class FieldCatalog;
struct RunContext;
} // namespace core::master

namespace core::master::plugin
{

// Where in the step an action can run.
enum class Phase : unsigned
{
    PreExchange = 1 << 0,  // before halo start (rare)
    Interior = 1 << 1,     // safe while halos are in flight
    PostExchange = 1 << 2, // requires fresh ghost cells
    PostBC = 1 << 3,       // after boundary conditions applied
    EndStep = 1 << 4       // final book-keeping / diagnostics
};

// Bitmask helpers
inline Phase operator|(Phase a, Phase b)
{
    return Phase(unsigned(a) | unsigned(b));
}
inline bool has(Phase mask, Phase bit)
{
    return (unsigned(mask) & unsigned(bit)) != 0;
}

// Declarative data-access (used by the scheduler for halos/I/O decisions)
struct Access
{
    std::vector<std::string> reads;
    std::vector<std::string> writes;
    // Per-field halo depth requirement (0 = interior only)
    std::unordered_map<std::string, int> halos;
};

struct ActionInfo
{
    std::string name; // stable id
    Phase phases;     // where it may run
    Access access;
};

// Tile action: coarse-grained compute on one sub-box.
struct IAction
{
    virtual ~IAction() = default;
    virtual const ActionInfo& info() const = 0;
    virtual void execute(const MeshTileView& tile, FieldCatalog& fields, double dt) = 0;
};

// Global action: reductions, solvers, etc. (no tiling)
struct IGlobal
{
    virtual ~IGlobal() = default;
    virtual std::string_view name() const = 0;
    virtual Phase phases() const = 0;
    virtual void run(FieldCatalog& fields, double dt) = 0;
};

using KV = std::unordered_map<std::string, std::string>;

} // namespace core::master::plugin
