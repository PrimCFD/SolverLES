#pragma once
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @file   ConfigYAML.hpp
 * @brief  YAML → AppConfig loader and schema for the solver app.
 *
 * @details
 * @rst
 * The **YAML configuration** drives the solver app and maps cleanly onto the core's
 * orchestration and I/O. This file defines:
 *
 * - :cpp:struct:`AppConfig` — the strongly-typed config object (mesh, time, io, program…)
 * - :cpp:func:`load_config_from_yaml` — loader that parses a YAML file into :cpp:struct:`AppConfig`
 *
 * **Schema (v0)**
 *
 * .. code-block:: yaml
 *
 *    case: <string>                 # case name (used for writer artifacts)
 *
 *    mesh:
 *      local: [nx, ny, nz]          # interior sizes (ints), per-rank
 *      ng: 2                        # uniform ghost width
 *      periodic: [false,false,false]# reserved
 *
 *    time:
 *      dt: 1.0e-3                   # seconds
 *      t_end: 1.0                   # seconds
 *      write_every:                 # choose one; if both set, 'steps' wins
 *        steps: 10                  # every N steps
 *        # time: 0.05               # OR every T seconds (ceil(T/dt))
 *
 *    io:
 *      backend: xdmf                # xdmf | cgns | null
 *      path: out                    # output directory
 *      precision: native            # native | float64 | float32 (packing override)
 *      async:
 *        enabled: true
 *        max_queue: 8               # 0 = unbounded
 *        drop_on_overflow: true     # protect cadence under backpressure
 *      preflight:
 *        enabled: true
 *        ram_bytes: auto            # auto | integer bytes
 *        disk_bytes: auto           # auto | integer bytes
 *
 *    plugins:
 *      - lib: libphysics_fluids.so  # physics DSOs to load (in order)
 *
 *    program:
 *      key: rk3
 *      params:                      # free-form KV (string → string)
 *        cfl: "0.7"
 *
 *    fields:
 *      output: [rho]                # names to register & write
 *
 * **Semantics**
 *
 * - ``write_every.steps`` takes precedence; else ``write_every.time`` derives
 *   :cpp:member:`AppConfig::write_every_steps` via ``ceil(time / dt)``.
 * - ``precision: native`` keeps element sizes; ``float32``/``float64`` request a packing cast.
 * - ``async`` wraps whichever backend is chosen using a bounded worker queue.
 * - ``preflight`` estimates bytes-per-step (via WritePlan) and checks RAM/disk headroom.
 *
 * **Core integration**
 *
 * - :cpp:func:`core::master::Master::set_writer` consumes a writer built from ``io``
 * (XDMF/CGNS/Null).
 * - :cpp:func:`core::master::Master::configure_program` consumes ``program.key`` and
 * ``program.params``.
 * - :cpp:class:`core::master::FieldCatalog` selections come from ``fields.output``.
 *
 * See the *Core Orchestration & I/O — Developer Guide* for how this flows through Master/Scheduler
 * and writers.  The formatting of this comment block follows our docs style guide
 * (Doxygen + embedded reST).  For authoring tips, see *Documentation — build, workflow, and style
 * guide*.
 * @endrst
 */

struct AppConfig
{
    // mesh
    std::array<int, 3> local{32, 32, 32};
    int ng = 2;
    std::array<bool, 3> periodic{false, false, false};

    // time
    double dt = 1e-3;
    double t_end = 1e-2;
    int write_every_steps = 5;              // derived
    std::optional<double> write_every_time; // optional

    // io
    enum class Backend
    {
        Null,
        Xdmf,
        Cgns
    };
    enum class Precision
    {
        Native,
        F32,
        F64
    };
    struct IO
    {
        Backend backend = Backend::Null;
        std::string path = "out";
        Precision precision = Precision::Native;
        // XDMF XML index version selector for the XDMF backend ("v2" or "v3")
        std::string xdmf_version = "v2";

        struct Async
        {
            bool enabled = false;
            int max_queue = 0;
            bool drop_on_overflow = false;
        } async;

        struct Preflight
        {
            bool enabled = false;
            std::optional<std::size_t> ram_bytes;
            std::optional<std::size_t> disk_bytes;
        } preflight;

    } io;

    // plumbing
    std::string case_name = "case";
    std::vector<std::string> plugin_libs{};
    std::string program_key = "noop";
    std::unordered_map<std::string, std::string> program_params{}; // <— map, not vector
    std::vector<std::string> fields_output{"rho"};
};

static inline std::string to_lower(std::string s)
{
    for (auto& c : s)
        c = (char) std::tolower((unsigned char) c);
    return s;
}
static inline AppConfig::Backend parse_backend(const std::string& s)
{
    auto v = to_lower(s);
    if (v == "xdmf")
        return AppConfig::Backend::Xdmf;
    if (v == "cgns")
        return AppConfig::Backend::Cgns;
    return AppConfig::Backend::Null;
}
static inline AppConfig::Precision parse_precision(const std::string& s)
{
    auto v = to_lower(s);
    if (v == "float32" || v == "f32")
        return AppConfig::Precision::F32;
    if (v == "float64" || v == "f64")
        return AppConfig::Precision::F64;
    return AppConfig::Precision::Native;
}

inline AppConfig load_config_from_yaml(const std::string& path)
{
    AppConfig cfg;
    YAML::Node root = YAML::LoadFile(path);

    if (auto n = root["case"])
        cfg.case_name = n.as<std::string>();

    if (auto m = root["mesh"])
    {
        if (auto L = m["local"])
        {
            auto v = L.as<std::vector<int>>();
            if (v.size() == 3)
                cfg.local = {v[0], v[1], v[2]};
        }
        if (auto N = m["ng"])
            cfg.ng = N.as<int>();
        if (auto P = m["periodic"])
        {
            auto v = P.as<std::vector<bool>>();
            if (v.size() == 3)
                cfg.periodic = {v[0], v[1], v[2]};
        }
    }

    if (auto t = root["time"])
    {
        if (auto n = t["dt"])
            cfg.dt = n.as<double>();
        if (auto n = t["t_end"])
            cfg.t_end = n.as<double>();
        if (auto we = t["write_every"])
        {
            bool set = false;
            if (auto s = we["steps"])
            {
                cfg.write_every_steps = std::max(1, s.as<int>());
                set = true;
            }
            if (!set)
            { // fall back to time→steps
                if (auto ti = we["time"])
                    cfg.write_every_time = ti.as<double>();
            }
        }
    }

    if (auto I = root["io"])
    {
        if (auto n = I["backend"])
            cfg.io.backend = parse_backend(n.as<std::string>());
        if (auto n = I["path"])
            cfg.io.path = n.as<std::string>();
        if (auto n = I["precision"])
            cfg.io.precision = parse_precision(n.as<std::string>());
        if (auto n = I["xdmf_version"])
        {
            auto v = to_lower(n.as<std::string>());
            // accept only v2 or v3; keep default otherwise
            if (v == "v2" || v == "2")
                cfg.io.xdmf_version = "v2";
            else if (v == "v3" || v == "3")
                cfg.io.xdmf_version = "v3";
        }

        if (auto A = I["async"])
        {
            if (auto n = A["enabled"])
                cfg.io.async.enabled = n.as<bool>();
            if (auto n = A["max_queue"])
                cfg.io.async.max_queue = n.as<int>();
            if (auto n = A["drop_on_overflow"])
                cfg.io.async.drop_on_overflow = n.as<bool>();
        }
        if (auto P = I["preflight"])
        {
            if (auto n = P["enabled"])
                cfg.io.preflight.enabled = n.as<bool>();
            if (auto n = P["ram_bytes"])
            {
                const auto s = to_lower(n.as<std::string>());
                if (s != "auto")
                    cfg.io.preflight.ram_bytes = (std::size_t) std::stoull(s);
            }
            if (auto n = P["disk_bytes"])
            {
                const auto s = to_lower(n.as<std::string>());
                if (s != "auto")
                    cfg.io.preflight.disk_bytes = (std::size_t) std::stoull(s);
            }
        }
    }

    if (auto P = root["plugins"])
    {
        for (const auto& item : P)
        {
            if (auto n = item["lib"])
                cfg.plugin_libs.push_back(n.as<std::string>());
        }
    }

    if (auto pr = root["program"])
    {
        if (auto n = pr["key"])
            cfg.program_key = n.as<std::string>();
        if (auto K = pr["params"])
        {
            for (auto it = K.begin(); it != K.end(); ++it)
                cfg.program_params.emplace(it->first.as<std::string>(),
                                           it->second.as<std::string>());
        }
    }

    if (auto F = root["fields"])
    {
        if (auto out = F["output"])
            cfg.fields_output = out.as<std::vector<std::string>>();
    }

    // derive write_every_steps from time if provided
    if (cfg.write_every_time.has_value() && cfg.dt > 0.0)
        cfg.write_every_steps = std::max(1, (int) std::ceil(*cfg.write_every_time / cfg.dt));

    return cfg;
}
