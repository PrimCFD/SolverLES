# `src/core` Core architecture

This document expands the high-level view of the **core orchestration** that lives under `src/core/`
and how it connects to apps, plugins, and I/O. It is intentionally concise and practical.

---

## Directory layout (relevant parts)

```text
src/
├─ apps/
│ └─ solver/ # entry points (executables)
├─ core/
│ ├─ include/
│ │ ├─ master/ # orchestration public API
│ │ │ ├─ Master.hpp # façade for apps
│ │ │ ├─ Scheduler.hpp # time loop + phase runner
│ │ │ ├─ RunContext.hpp # MPI/stream/memory handles
│ │ │ ├─ FieldCatalog.hpp # field registry + POD views
│ │ │ ├─ Views.hpp # POD types passed across ABI
│ │ │ ├─ plugin/ # domain-agnostic plugin API
│ │ │ │ ├─ Action.hpp # IAction / IGlobal + Access/Phase
│ │ │ │ ├─ Program.hpp # IProgram → StepPlan
│ │ │ │ └─ Registry.hpp # factories + register symbol
│ │ │ └─ io/
│ │ │ ├─ IWriter.hpp # I/O abstraction
│ │ │ └─ NullWriter.hpp # baseline no-op
│ │ └─ (mesh/, memory/ … existing)
│ └─ src/
│ └─ master/
│ ├─ FieldCatalog.cpp
│ ├─ Master.cpp
│ ├─ PluginHost.cpp
│ ├─ Scheduler.cpp
│ └─ plugin/Registry.cpp
└─ physics/… (optional shared libraries loaded at runtime)
```

---

## Big picture

- **Core is domain-agnostic.** It owns the run loop, halo/BC call sites, memory residency hooks,
  and I/O cadence. It never references “flux/LES/time” or other domain terms.
- **Plugins** provide *actions* and *programs*:
  - An **action** is a coarse-grained operation that runs **per tile** (e.g. “compute residuals”,
    “update SGS”, “apply source”). Actions declare their **data access** and **phase**.
  - A **program** assembles actions into a **StepPlan** for the scheduler (e.g. explicit RK3, PISO,
    Jacobi relaxation). The program can be replaced without touching the core.
- **FieldCatalog** exposes POD *views* to plugins and writers; ownership remains in the app/core.
- **RunContext** carries low-level handles (MPI communicator, device stream, MemoryManager*).

---

## Run loop (phases)

At each step, the scheduler:

1. Builds/gets a `StepPlan` from the selected `IProgram`.
2. `halos.start()` → post asynchronous exchanges.
3. Run tiled actions whose `Phase` includes **Interior** (can overlap with halo traffic).
4. `halos.finish()` → complete exchanges.
5. Apply **BCs**.
6. Run remaining tiled/global actions for **PostExchange**, **PostBC**, **EndStep** phases.
7. Trigger **I/O** at the configured cadence.

```cpp
for (s : steps) {
  auto plan = program.plan_step(dt);

  halos.start();
  for (a : plan.tiled if a.phase&Interior) a.execute(tile, fields, dt);
  halos.finish();

  bcs.apply();

  for (a : plan.tiled if a.phase&PostBC) a.execute(tile, fields, dt);
  for (g : plan.globals if g.phase&PostBC) g.run(fields, dt);

  if (s % write_every == 0) writer.write(request(fields));
}
```
---

## Plugins in one minute

- A plugin shared library exports: `extern "C" bool physics_register_v1(Registry*)`.

- During registration, it adds factories keyed by strings (what YAML selects).

- The app loads one or more libraries via `PluginHost::load_library(...)` and builds a program:
`auto P = host.make_program("noop", cfg, rc)`.

- The core calls into the program to obtain the per-step plan; each action gets passed a `MeshTileView` + `FieldCatalog` and runs the entire kernel sequence inside the DSO (no per-cell callbacks).

**Access and halos**: Each action declares the fields it reads/writes and halo depths it needs. The scheduler uses this info to place the action in the right phase and to ensure ghost cells are valid before execution.

---

## Memory and residency

- Plugins operate on host pointers and (if needed) request residency through `MemoryManager::{to_device,to_host}` using the stream in the view.

- With Unified Memory builds, those become prefetches; with mirrored builds, they are async H2D/D2H copies.

- Plugins do not allocate/free core-owned field arrays.

---

## I/O

- Writers implement IWriter (open_case, write, close). The default is NullWriter.

- Writers receive a WriteRequest containing a list of AnyFieldView and step/time metadata.

- An AsyncWriter can be added later without changing the interface (wrapper pattern).

---

## Tests

- Unit tests cover: field catalog, scheduler cadence, and seam stubs for halo/BC.

- Integration tests add a sample plugin to exercise the loader and registry.

- Performance micro-benchmarks (future work) measure ABI call cost and overlap efficiency.

---

## Extending

- Add an action: implement IAction, fill out ActionInfo (phase + access), register it.

- Add a program: implement IProgram::plan_step, returning a `StepPlan` with your actions
arranged per sub-stage; register under a unique key.

- Add I/O: subclass `IWriter`, pass it to `Master::set_writer(...)`.

---

## Conventions

- POD views across the ABI (no STL/Eigen in signatures).

- Coarse-grained calls per tile/step (no per-cell or small virtual calls across DSOs).

- Use the core stream for all transfers and kernels.

- Keep headers light (forward declarations); include heavy headers only in `.cpp` to keep optional MPI/CUDA clean.