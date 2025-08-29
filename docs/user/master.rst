.. _master-architecture:

Master / Orchestration Architecture
===================================

This page explains the **domain-agnostic orchestration** implemented under :code:`src/core`.
It focuses on *how* the scheduler, plugins, and I/O cooperate. For memory, fields, and mesh
details see :ref:`memory`.

Overview
--------

The core owns the **time loop**, **halo/BC call sites**, **data residency**, and **I/O cadence**.
It never references domain-specific constructs (e.g., “flux/LES/time”). Instead, domain logic
arrives through **plugins** providing *actions* (tile-wise computations) and *programs*
(step planners).

Key components
--------------

- :cpp:class:`core::master::Master` — façade used by apps; wires a writer, a program, and runs.
- :cpp:class:`core::master::Scheduler` — advances time; overlaps halo exchanges with interior work;
  enforces :ref:`phases <phases>`.
- :cpp:class:`core::master::FieldCatalog` — field registry; exposes POD views to plugins and writers.
- :cpp:struct:`core::master::RunContext` — runtime handles (MPI communicator, device stream,
  memory manager).
- :cpp:class:`core::master::PluginHost` — loads shared libraries; holds the
  :cpp:class:`core::master::plugin::Registry` of factories.
- :cpp:class:`core::master::io::IWriter` — output abstraction (rank-local VTK/HDF5, etc.).

.. _phases:

Phases and step plan
--------------------

Plugins register **actions** with **phase** tags and **access** metadata.
A **program** assembles a per-step :cpp:struct:`core::master::plugin::StepPlan` that the scheduler executes.

.. code-block:: none

   PreExchange   : runs before halo posts (rare)
   Interior      : runs while halos are in flight (safe on interior cells)
   PostExchange  : runs after halo completion, before BCs
   PostBC        : runs after boundary conditions
   EndStep       : final hooks (diagnostics, CFL checks, etc.)

Execution flow
--------------

.. graphviz::

   digraph G {
     rankdir=LR;
     node [shape=box];
     plan   [label="program.plan_step(dt)"];
     halos1 [label="halos.start()"];
     int    [label="run tiled actions (Interior)"];
     halos2 [label="halos.finish()"];
     bcs    [label="bcs.apply()"];
     post   [label="run tiled/globals (PostExchange/PostBC/EndStep)"];
     io     [label="writer.write() (cadence)"];

     plan -> halos1 -> int -> halos2 -> bcs -> post -> io;
   }

Data access and halos
---------------------

Each action declares its **reads**, **writes**, and required **halo depths** per field. The scheduler
uses this metadata to order actions and ensure ghost zones are valid before execution. Interior-only
actions can overlap with halo traffic.

Memory residency
----------------

Plugins operate on **views** exposing host pointers and strides. When they need the data on the
device (or back on the host), they call :cpp:func:`core::memory::MemoryManager::to_device()` /
:cpp:func:`core::memory::MemoryManager::to_host()` with the provided stream. In *Unified Memory*
builds, these are **prefetch** operations; in mirrored builds, they are **async copies**.

I/O
---

Writers implement :cpp:class:`core::master::io::IWriter` and receive a
:cpp:struct:`core::master::io::WriteRequest`. The default :cpp:class:`core::master::io::NullWriter`
is useful for performance baselines and smoke tests.

Extending
---------

- **New action:** implement :cpp:class:`core::master::plugin::IAction`, describe your phase + access,
  and register a factory under a string key.
- **New program:** implement :cpp:class:`core::master::plugin::IProgram` returning a
  :cpp:struct:`core::master::plugin::StepPlan` per step.
- **New writer:** implement :cpp:class:`core::master::io::IWriter` and pass it to
  :cpp:func:`core::master::Master::set_writer`.

Testing
-------

Unit tests verify `FieldCatalog` semantics and scheduler cadence on a **builtin noop program**.
Halo exchange and boundary application are exercised through **link-time stubs** in tests so
MPI remains optional in default builds. Integration tests can load a sample plugin to exercise
runtime registry & loading.