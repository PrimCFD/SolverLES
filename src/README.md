# `src/` — Core, Plugins & Apps

This README explains how the **domain‑agnostic core** orchestrates a run, how **physics plugins** implement actions/programs, and how **apps** wire everything together. It also reflects the actual shape of the `src/` tree so newcomers can find things quickly.

---

## Directory layout (source tree)

```text
src/
├─ apps/                       # Executable entry points (CLIs)
│  └─ solver/                  # Reference solver app
├─ core/                       # Orchestration, memory, mesh, I/O, plugin host
│  ├─ include/
│  │  ├─ master/               # façade API (Master, Scheduler, RunContext, FieldCatalog, Views)
│  │  ├─ io/                   # writer contracts & configs
│  │  ├─ memory/               # MemoryManager & allocation policies
│  │  └─ mesh/                 # structured grid, halos, boundary ops
│  └─ src/master/              # implementations (no heavy deps in public headers)
├─ physics/                    # Runtime‑loaded shared libraries per domain
│  ├─ fluids/                  # example: actions (flux/sgs/time) + programs
│  └─ …                        # other domains (heat, EM, multiphysics)
├─ kernels/                    # Shared math kernels (Fortran/C/CUDA), reused by plugins
├─ bindings/                   # C/Fortran/CUDA interop shims & headers
├─ ipc/                        # App↔GUI IPC (optional)
└─ gui/                        # Qt/VTK front‑end (optional build)
```

> **Notes**
>
> * `core/` is **domain‑agnostic** and exposes only ABI‑stable POD types to plugins/writers.
> * Numerics (fluxes, SGS, time schemes, etc.) live **inside** each physics plugin.
> * Apps are thin: they parse a case, allocate fields, load plugins, pick a program, and run.

---

## Big picture

* The **core** owns: run loop & phases, halo/BC call sites, field registry/views, memory residency, and writer cadence.
* **Plugins** contribute:

  * **Actions** — coarse‑grained kernels that run per **tile** and declare their *Access* & required halo depths.
  * **Programs** — assemble actions into a per‑step **StepPlan** (e.g., explicit RK, PISO).
* **Apps** (e.g., `apps/solver`) load one or more physics DSOs, select a program by key, and turn the crank.

---

## Orchestration — step phases

At each time step the scheduler executes phases in order (MPI halos overlap with interior work):

1. Build or fetch the **StepPlan** from the active program.
2. `halos.start()` — post asynchronous exchanges (MPI builds).
3. Run tiled actions tagged **Interior** (overlaps with halo traffic).
4. `halos.finish()` — ensure ghost cells are valid.
5. Apply **boundary conditions**.
6. Run remaining tiled/global actions for **PostExchange**, **PostBC**, **EndStep**.
7. Trigger **I/O** at the configured cadence.

---

## Public API (what plugins & apps see)

### Minimal POD views across the ABI

```cpp
struct FieldView {
  void* host;     // core‑owned host pointer (UM or host mirror)
  void* device;   // device mirror or same as host in UM mode
  int   nx, ny, nz;
  ptrdiff_t sx, sy, sz; // byte strides
  void* stream;   // device stream/queue (opaque to the ABI)
};

struct MeshTileView { int i0,i1, j0,j1, k0,k1; /* + metadata */ };
```

*Plugins receive only trivial structs and never own storage.*

### FieldCatalog

A registry mapping names → `FieldView`s. Apps register fields once, select which ones should be written, and pass the catalog to the scheduler/plugins.

### RunContext

Opaque handles the core threads through the system: MPI communicator, device stream/queue, and a `MemoryManager*`. Public headers keep these as `void*` to avoid heavy includes.

### MemoryManager (residency)

Centralizes allocation and host↔device movement for two primary strategies:

* **Unified Memory** — `cudaMallocManaged` + `cudaMemPrefetchAsync`.
* **Mirrored** — explicit device mirrors + async H2D/D2H copies (host buffers may be pinned).

Plugins request movement via `to_device(host_ptr, bytes, stream)` / `to_host(...)`; the *same stream* flows through all subsystems to keep ordering correct.

### Plugin interfaces

* **Registration symbol** in each DSO: `extern "C" bool physics_register_v1(Registry*);`
* **Actions** implement `execute(tile, fields, dt)` and declare `ActionInfo { phases, access, halo }`.
* **Programs** implement `plan_step(dt)` and return a `StepPlan` describing tiled and global work.

### Writer interfaces

Implement `IWriter` (`open_case / write / close`) and optionally wrap with `AsyncWriter` for threaded output. `WritePlan` synthesizes layout/packing (incl. precision casts) from selected `FieldView`s.

---

## How apps wire it together

1. Build a `Mesh` and allocate ghosted arrays for all fields.
2. Register those arrays in `FieldCatalog` and select outputs.
3. Create a writer (e.g., XDMF/HDF5 or CGNS) and hand it to `Master`.
4. Load physics DSOs via `PluginHost`.
5. Choose a **program** & parameters by registration key (usually from YAML).
6. Run the loop with `TimeControls { dt, t_end, write_every }`.

```cpp
core::master::Master M(rc, mesh);
M.fields().register_scalar("rho", rho_ptr, sizeof(double), {nx,ny,nz}, {8,8*nx,8*nx*ny});
M.set_writer(std::make_unique<io::XdmfHdf5Writer>(cfg));
M.load_plugin_library("libphysics_fluids.so");
M.configure_program("explicit_rk3", {{"cfl","0.7"}});
M.run({ .dt=1e-3, .t_end=1.0, .write_every=10 });
```

---

## YAML — selecting programs & actions

The solver app reads a YAML file and passes **string keys** to the registry. Example :

```yaml
#=======================
# SolverLES config (v0)
#=======================


# Short identifier for outputs (writer will create path/<case>.h5, path/<case>.xmf, etc.)
case: cavity_64 # default: "case"


# Mesh (local box on each rank, ghosts are uniform)
mesh:
  local: [64, 64, 64] # required (ints)
  ng: 2 # required (uniform ghost width)
  periodic: [false, false, false] # optional (reserved)


# Time integration and writer cadence
time:
  dt: 1.0e-3 # required (seconds)
  t_end: 1.0 # required (seconds)
  write_every: # choose one; if both set, 'steps' wins
    steps: 10 # write every N steps
    # time: 0.05 # OR every T seconds (converted with ceil(T/dt))


# I/O policy & backend
io:
  backend: xdmf # xdmf | cgns | null
  path: out # output directory
  precision: native # native | float64 | float32 (packing override)
  async: # optional background writer
    enabled: true
    max_queue: 8 # 0 = unbounded
    drop_on_overflow: true # drop writes when queue is full (protects cadence)
  preflight: # optional resource checks (RAM/disk)
    enabled: true
    ram_bytes: auto # auto | integer bytes (e.g., "17179869184")
    disk_bytes: auto # auto | integer bytes


# Runtime-loaded physics libraries (in this order)
plugins:
- lib: libphysics_fluids.so
# - lib: libphysics_heat.so


# Program selector + free-form KV passed to plugin (strings on both sides)
program:
  key: rk3
  params:
    cfl: "0.7"
    # any other plugin-defined keys → strings


# Field output selection by name (the app allocates/registers these as doubles for now)
fields:
  output: [rho]
```

---

## Relationship to `kernels/`, `bindings/`, `ipc/`, and `gui/`

* **kernels/** holds shared numerical kernels (Fortran/C/CUDA). Plugins call into these to keep hot loops **inside** the DSO and avoid chatty ABI calls.
* **bindings/** provides the C/Fortran interop shims used by both core and plugins (e.g., `bind(C)`, ISO C pointers).
* **ipc/** is used by apps/GUI to exchange progress/controls; the core itself doesn’t depend on it.
* **gui/** is an optional front‑end that drives the same `Master` façade as the CLI app.

---

## CMake targets (at a glance)

* `core` — static library (orchestration, memory, mesh, I/O contracts).
* `physics_<domain>` — shared libraries (plugins: actions + programs).
* `solver` — the reference CLI app (parses YAML → sets up `Master`).
* Writers (optional): `XdmfHdf5Writer`, `CGNSWriter`; decorator: `AsyncWriter`.

---

## Extending the system

**New action**

1. In a plugin, implement `IAction` and its `ActionInfo` (phases, access, halo depths).
2. Use the `FieldView` strides/sizes and request residency via `MemoryManager` before kernels.
3. Register it under a string key during `physics_register_v1`.

**New program**

1. Implement `IProgram::plan_step(dt)` to return a `StepPlan` (tiled actions + globals per phase).
2. Register it under a program key.

**New writer**

1. Subclass `IWriter`, build a `WritePlan` on first call, use `StagingPool` if packing is needed.
2. Optionally wrap with `AsyncWriter` for background I/O with backpressure.

**New domain plugin**

1. Create `physics/<domain>/` with actions/programs and any domain‑specific kernels.
2. Export `physics_register_v1` and add a shared library target.

---

## Conventions & pitfalls

* Keep ABI surfaces trivial: POD structs, byte strides, no STL containers across DSO boundaries.
* Use **coarse‑grained** calls (per tile/sub‑stage). Don’t call back into the core per cell.
* Field storage is owned by the app/core; plugins **must not** allocate/free those arrays.
* Keep `<mpi.h>`, CUDA, HDF5, CGNS includes out of public headers; guard usage in `.cpp` only.
* Ghost width is **uniform** in all directions; BC/halo helpers assume this.
* Always use the **provided stream** for transfers and kernels to maintain ordering.

---

## Glossary

* **Tile** — rectangular sub‑box of the local mesh the scheduler passes to actions.
* **View** — non‑owning storage description (name, pointer, extents, byte strides).
* **Plan** — synthesized description of work (StepPlan) or output (WritePlan).

---

## See also

* `src/core/` Developer Guide for deep dives into MemoryManager, mesh layout, writers, and feature flags.
* Repository‑level README for quick‑start, CI, and developer workflow.
