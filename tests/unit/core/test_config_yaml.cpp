#include <catch2/catch_all.hpp>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "master/io/ConfigYAML.hpp" // <-- same header main.cpp uses

namespace fs = std::filesystem;

static std::string write_temp_yaml(const std::string& stem, const std::string& content)
{
    const auto tmp =
        fs::temp_directory_path() / (stem + "_" + std::to_string(std::rand()) + ".yml");
    std::ofstream out(tmp);
    out << content;
    out.close();
    return tmp.string();
}

TEST_CASE("config.hpp parses hello_mesh.yml correctly", "[config][yaml]")
{
    // NOTE: proper indentation matters in YAML.
    const char* YAML_TXT = R"YAML(
#=======================
# KolmoPlas config (v0)
#=======================

# Short identifier for outputs (writer will create path/<case>.h5, path/<case>.xmf, etc.)
case: cavity_64  # default: "case"

# Mesh (local box on each rank, ghosts are uniform)
mesh:
  local: [64, 64, 64]            # required (ints)
  ng: 2                          # required (uniform ghost width)
  periodic: [false, false, false]# optional (reserved)

# Time integration and writer cadence
time:
  dt: 1.0e-3          # required (seconds)
  t_end: 1.0          # required (seconds)
  write_every:        # choose one; if both set, 'steps' wins
    steps: 10         # write every N steps
    # time: 0.05

# I/O policy & backend
io:
  backend: xdmf       # xdmf | cgns | null
  path: out           # output directory
  precision: native   # native | float64 | float32
  async:
    enabled: true
    max_queue: 8
    drop_on_overflow: true
  preflight:
    enabled: true
    ram_bytes: auto
    disk_bytes: auto

# Runtime-loaded physics libraries (in this order)
plugins:
  - lib: libphysics_fluids.so

# Program selector + free-form KV passed to plugin (strings on both sides)
program:
  key: rk3
  params:
    cfl: "0.7"

# Field output selection by name
fields:
  output: [rho]
)YAML";

    const std::string path = write_temp_yaml("hello_mesh", YAML_TXT);
    const AppConfig cfg = load_config_from_yaml(path);

    // Capture everything for easier failure diagnostics
    CAPTURE(cfg.case_name);
    CAPTURE(cfg.local[0], cfg.local[1], cfg.local[2]);
    CAPTURE(cfg.ng, cfg.periodic[0], cfg.periodic[1], cfg.periodic[2]);
    CAPTURE(cfg.dt, cfg.t_end, cfg.write_every_steps, cfg.write_every_time.has_value());
    CAPTURE((int) cfg.io.backend, cfg.io.path, (int) cfg.io.precision);
    CAPTURE(cfg.io.async.enabled, cfg.io.async.max_queue, cfg.io.async.drop_on_overflow);
    CAPTURE(cfg.io.preflight.enabled, cfg.io.preflight.ram_bytes.has_value(),
            cfg.io.preflight.disk_bytes.has_value());
    CAPTURE(cfg.plugin_libs.size(), cfg.program_key);
    CAPTURE(cfg.fields_output);

    // Case
    CHECK(cfg.case_name == std::string("cavity_64"));

    // Mesh
    CHECK(cfg.local[0] == 64);
    CHECK(cfg.local[1] == 64);
    CHECK(cfg.local[2] == 64);
    CHECK(cfg.ng == 2);
    CHECK(cfg.periodic[0] == false);
    CHECK(cfg.periodic[1] == false);
    CHECK(cfg.periodic[2] == false);

    // Time
    CHECK(cfg.dt == Catch::Approx(1.0e-3));
    CHECK(cfg.t_end == Catch::Approx(1.0));
    CHECK(cfg.write_every_steps == 10); // 'steps' wins
    CHECK_FALSE(cfg.write_every_time.has_value());

    // IO
    CHECK(cfg.io.backend == AppConfig::Backend::Xdmf);
    CHECK(cfg.io.path == std::string("out"));
    CHECK(cfg.io.precision == AppConfig::Precision::Native);
    CHECK(cfg.io.async.enabled == true);
    CHECK(cfg.io.async.max_queue == 8);
    CHECK(cfg.io.async.drop_on_overflow == true);
    CHECK(cfg.io.preflight.enabled == true);
    CHECK_FALSE(cfg.io.preflight.ram_bytes.has_value()); // 'auto' ⇒ std::nullopt
    CHECK_FALSE(cfg.io.preflight.disk_bytes.has_value());

    // Plugins
    REQUIRE(cfg.plugin_libs.size() == 1);
    CHECK(cfg.plugin_libs[0] == std::string("libphysics_fluids.so"));

    // Program
    CHECK(cfg.program_key == std::string("rk3"));
    // AppConfig::program_params can be unordered_map<string,string> OR vector<pair<...>> depending
    // on your header. Make this robust to either:
    bool found_cfl = false;
    if constexpr (std::is_same_v<decltype(cfg.program_params),
                                 std::unordered_map<std::string, std::string>>)
    {
        auto it = cfg.program_params.find("cfl");
        found_cfl = (it != cfg.program_params.end() && it->second == "0.7");
    }
    else
    {
        for (const auto& kv : cfg.program_params)
            if (kv.first == "cfl" && kv.second == "0.7")
            {
                found_cfl = true;
                break;
            }
    }
    CHECK(found_cfl);

    // Fields
    CHECK(cfg.fields_output == std::vector<std::string>({"rho"}));
}

TEST_CASE("config.hpp derives write_every from time when provided", "[config][yaml]")
{
    const char* YAML_TXT = R"YAML(
case: tcase
mesh: { local: [8,8,8], ng: 1 }
time: { dt: 0.01, t_end: 1.0, write_every: { time: 0.05 } }
io: { backend: null, path: out, precision: native }
plugins: []
program: { key: noop, params: {} }
fields: { output: [] }
)YAML";
    const std::string path = write_temp_yaml("hello_time", YAML_TXT);
    const AppConfig cfg = load_config_from_yaml(path);
    CHECK(cfg.write_every_steps == 5); // ceil(0.05/0.01)=5
}
