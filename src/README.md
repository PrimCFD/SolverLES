# `src/` Core, Plugins & Apps

This document describes the source tree under `src/`, focusing on the **core orchestration**, the **plugin model** (actions/programs), and how apps wire everything together. It replaces older notes that referred to hard‑coded domain terms inside the core.

---

## Directory layout (relevant parts)

```text
src/
├─ apps/
│  └─ solver/                 # Executables / entry points
├─ core/                      # Orchestration: run loop, halos, memory, registry, I/O
│  ├─ include/
│  │  ├─ master/              # Public API façade for apps
│  │  │  ├─ Master.hpp        # top-level façade
│  │  │  ├─ Scheduler.hpp     # step loop + phase runner
│  │  │  ├─ RunContext.hpp    # MPI/stream/MemoryManager handles
│  │  │  ├─ FieldCatalog.hpp  # field registry + POD views
│  │  │  ├─ Views.hpp         # POD types passed across ABI
│  │  │  └─ plugin/
│  │  │     ├─ Action.hpp     # IAction / IGlobal + Access/Phase
│  │  │     ├─ Program.hpp    # IProgram → StepPlan
│  │  │     └─ Registry.hpp   # factories + register symbol
│  │  └─ io/
│  │     ├─ IWriter.hpp       # I/O abstraction
│  │     └─ NullWriter.hpp    # baseline no-op
│  └─ src/master/             # implementations (Master, Scheduler, PluginHost, Registry)
└─ physics/…                  # Optional shared libraries loaded at runtime (per domain)
```

> Notes
>
> * The **core** is *domain‑agnostic*. It doesn’t know about “flux/LES/time”; those are just
>   examples of **actions** implemented inside physics plugins.
> * Each physics plugin (e.g., fluids, heat, EM) owns its numerics and hot kernels and is built as a
>   **shared library**. The app loads one or more plugins at runtime.

---

## Big picture

* The **core** owns the run loop, halo/BC call sites, field registry/views, memory residency hooks, and I/O cadence.
* **Plugins** provide two things:

  * **Actions**: coarse‑grained operations that run **per tile** (e.g., compute residuals, update SGS, apply sources). Each declares its **data access** and **phase**.
  * **Programs**: assemble actions into a **StepPlan** for the scheduler (e.g., explicit RK, PISO, relaxation). Swapping the program changes the cadence without touching the core.
* **FieldCatalog** exposes POD *views* (pointers + extents/strides + stream) to plugins; ownership stays in the app/core.
* **RunContext** carries low‑level handles (MPI communicator, device stream, MemoryManager\*).

---

## Runtime flow (phases)

At each time step, the scheduler executes phases in order and tiles work when possible:

1. Build/get a `StepPlan` from the selected program.
2. `halos.start()` → post asynchronous exchanges.
3. Run tiled actions whose phase includes **Interior** (overlaps with halo traffic).
4. `halos.finish()` → ensure ghost cells are valid.
5. Apply **BCs**.
6. Run remaining tiled/global actions for **PostExchange**, **PostBC**, **EndStep**.
7. Trigger **I/O** at the configured cadence.

---

## Plugin API at a glance

* A plugin shared library must export: `extern "C" bool physics_register_v1(Registry*);`
* During registration, the plugin adds factories keyed by **strings** (these are what YAML selects).
* The app uses a `PluginHost` to load libraries and create a program; the **core** calls into the program each step to obtain the `StepPlan` and then executes the listed actions.
* **Actions**:

  * Receive a `MeshTileView` and a set of field views from the `FieldCatalog`.
  * Declare their **Access** (read/write) and **halo depths**; the scheduler places them in suitable phases.
  * Run entire kernel sequences **inside the DSO** (no per‑cell callbacks across the ABI).

### Data views passed across the ABI

Minimal POD types (illustrative):

```cpp
struct FieldView {
  void*    host;     // base host pointer (core‑owned)
  void*    device;   // device mirror or UM pointer
  int      nx, ny, nz;
  ptrdiff_t sx, sy, sz; // element strides
  void*    stream;   // device stream/queue (nullptr = default)
};

struct MeshTileView { int i0,i1, j0,j1, k0,k1; /* + metadata as needed */ };
```

---

## Memory residency & streams

* Field arrays are **owned by the core**; plugins never allocate/free them.
* `core::MemoryManager` maintains host↔device residency and supports two build modes:

  * **Unified Memory**: host pointers are UM; transfers are **prefetches**.
  * **Mirrored**: explicit device mirrors; transfers are **async H2D/D2H copies**.
* Typical calls from an action before/after kernels:

  * `MemoryManager::to_device(host_ptr, bytes, stream)`
  * `MemoryManager::to_host  (host_ptr, bytes, stream)`
* Always use the **stream** carried in the field view/run context for transfers and kernels.

**Do**

* Pass POD views (pointers, extents, strides, stream) into plugins.
* Keep hot loops **inside** the plugin (select once per tile/stage).
* Use core‑provided streams/queues.

**Don’t**

* Allocate/free core‑owned field arrays inside plugins.
* Make per‑cell calls back into the core.
* Pass STL/Eigen types across the ABI.

---

## YAML configuration (selector keys)

The app’s YAML chooses **programs** and **actions/models** by their **registration keys** (strings). For example, a fluids plugin might register `flux_central2`, `les_smag`, `time_piso`, and a program `rk3`. The core only sees opaque keys and runs the `StepPlan` provided by the selected program.

```yaml
mesh:
  nx: 64; ny: 64; nz: 64
  Lx: 1;  Ly: 1;  Lz: 1
  periodic: [false, false, false]

time:
  dt: 1e-3
  t_end: 1.0
  write_every: 0.05

linear solver:
  - name: poisson
    params:
      tol: 1e-10
      max_it: 200

physics:
  - physics: fluids # Matrix if multiphysics
    scheme: flux_central2
    les: [gaussian, smag] # [Filter, SGS model]
    time: time_piso_incomp

    bc: fluids # Matrix if multiphysics
      - { patch: 0, type: wall }
      - { patch: 1, type: inflow, U: [1, 0, 0] }
```

---

## I/O

* Implement `IWriter` to provide output; the default `NullWriter` is a no‑op.
* Writers receive a `WriteRequest` describing selected fields and step/time metadata.
* An async writer can be layered later without changing the public interface.

---

## Building blocks (CMake targets)

* **core** — static library implementing orchestration, registry, memory, I/O abstractions.
* **physics\_* plugins*\* — shared libraries; each registers actions and at least one program.
* **solver app** — executable; parses YAML, instantiates fields, loads plugins, selects a program, runs the loop.

---

## Tests

* Unit tests: field catalog behavior, scheduler cadence, halo/BC seams.
* Integration tests: sample plugin exercising loader/registry and a trivial program.
* (Planned) Micro‑benchmarks for ABI call cost and overlap efficiency.

---

## Conventions

* Use POD views across the ABI.
* Keep calls **coarse‑grained** (per tile/sub‑stage); avoid fine virtual calls across DSOs.
* Prefer forward declarations in headers; include heavy headers (`cuda.h`, MPI, etc.) only in `.cpp`.

---

## Migration notes (deprecations)

* Older references to `src/plugins/{flux,sgs,time}` as top‑level directories are obsolete. Those numerics now live **inside** each physics plugin’s own source tree.
* The core no longer names domain concepts (flux/LES/time); it orchestrates **actions** arranged by a **program**.
* Memory ownership and residency are centralized in `core::MemoryManager`; plugins operate on views and request residency via the provided stream.

---

## See also

* `docs/master.rst` for a deeper narrative about the orchestration and design.
* The repository‑level README for quick‑start, CI, and developer workflow.
