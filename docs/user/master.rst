.. _master:

Orchestration & Master
======================

This page describes the project’s orchestration layer in depth: a **single “Master” façade**
that wires **PLUGINS** (physics), a **Scheduler** (phases/time loop with halo/BC ordering),
a shared **RunContext** (MPI communicator + device stream + memory manager), a **FieldCatalog**
for non-owning views, and **I/O backends** (XDMF-ORCH/HDF5, CGNS-ORCH) with staging and optional
precision downcast. The design emphasizes **separation of concerns** (physics ↔ orchestration ↔ I/O),
**overlap of communication and compute** (non-blocking MPI), and **portable GPU residency**
(UM or mirrored buffers).  [BPG-ORCH]_, [MPI-NB-ORCH]_, [UM-API-ORCH]_, [XDMF-ORCH]_, [CGNS-ORCH]_.

High-level goals
----------------

- Keep kernels simple and fast; the runtime guarantees phase/halo/BC ordering.  [MPI-NB-ORCH]_
- Centralize ownership (RAII) and expose **non-owning views** for zero-overhead access. (Analogy: ``std::SPAN-ORCH`` / unmanaged views).  [SPAN-ORCH]_, [KOKKOS-VIEW-ORCH]_
- Overlap communication and compute with **non-blocking MPI**; avoid incidental global sync.  [MPI-NB-ORCH]_
- Offer **policy-agnostic residency**: explicit mirrors by default; **Unified Memory** as an opt-in with prefetch/advice.  [UM-API-ORCH]_
- Standard, tool-friendly output: **XDMF-ORCH (Light XML + Heavy HDF5)** or **CGNS-ORCH**.  [XDMF-ORCH]_, [CGNS-ORCH]_.

Core concepts
-------------

RunContext
^^^^^^^^^^
Opaque, immutable handles passed across the orchestration: MPI communicator, device
stream/queue, pointer to the MemoryManager, flags (CUDA-aware MPI, rank/size). Treat as
read-only during a run; rebuild the scheduler if handles change.

FieldCatalog
^^^^^^^^^^^^
Registry of named **non-owning** views (host pointer, element size **in bytes**,
logical extents, **byte** strides). The catalog also maintains an **output selection set**.
Ownership of storage remains with your application. (This keeps writer ABIs trivial
and decouples them from compile-time element types.)

PluginHost & Registry
^^^^^^^^^^^^^^^^^^^^^
Loads physics DSOs and fills a **Registry** with factories for programs/actions/globals.
A user-selected key + KV options produce a **Program** that assembles a per-step plan.
This “registry + dlopen” pattern is standard in performant, extensible C/C++ codes.  [PLUGINS-ORCH]_.

Program & Actions
^^^^^^^^^^^^^^^^^
A **Program** returns a **StepPlan** each step:

- **Tiled actions**: run on sub-boxes in specific **phases** (``PreExchange``, ``Interior``,
  ``PostExchange``, ``PostBC``, ``EndStep``).
- **Global actions**: run once per step (reductions/diagnostics).

Each action declares **Access**: reads/writes and **required halo depth** per field. The
scheduler uses this to order exchange/BCs and to decide what can run during overlap.

Scheduler & phases
------------------

Step order (per rank):

1. ``PreExchange`` hooks (rare)
2. Post non-blocking halo **exchange** for fields that need ghosts
3. ``Interior`` actions **while exchange progresses** (overlap)
4. Complete exchange → apply **boundary conditions** (scalar/vector)
5. ``PostBC`` then ``EndStep`` hooks
6. On cadence: issue a **write** of the selected fields

This is the canonical non-blocking stencil workflow; overlap is **implementation-dependent**,
but this pattern exposes the opportunity and preserves correctness.  [MPI-NB-ORCH]_.

Parallel & overlap (MPI)
------------------------

- Use **MPI_Isend/Irecv** + completion (``MPI_Wait*``) to progress halos; keep interior
  kernels independent of in-flight ghosts.  [MPI-NB-ORCH]_
- Two-wave exchange: faces first; the completion handles edges/corners (wide halos extend
  the band accordingly).  [HALO-PATTERN-ORCH]_
- Avoid unconditional barriers inside actions. If only rank 0 writes, **broadcast a status**
  before any barrier to keep ranks in lock-step.  [MPI-NB-ORCH]_.

Accelerators & residency
------------------------

Unified Memory (opt-in)
^^^^^^^^^^^^^^^^^^^^^^^
UM simplifies ownership (single pointer), but **prefetch** to GPU for compute, and **prefetch**
to CPU before heavy I/O to avoid on-demand migrations. ``cudaMemAdvise`` hints (preferred
location / read-mostly / accessed-by) help reduce thrash.  [UM-API-ORCH]_.

Explicit mirrors (default)
^^^^^^^^^^^^^^^^^^^^^^^^^^
Host owns canonical SoA arrays; device mirrors are persistent and updated with async copies on
the shared stream. For MPI, **CUDA-aware** stacks can send/recv device or UM buffers directly;
otherwise exchange stages through pinned host buffers.

Pinned & streams
^^^^^^^^^^^^^^^^
Use **pinned** host buffers with **``cudaMemcpyAsync`` in streams** to overlap copies with kernels
and to reach high PCIe/NVLink bandwidth. Pageable memory degrades to sync paths.  [BPG-ORCH]_.

Boundary conditions
-------------------
Face-wise operators for scalars/vectors:

- **Dirichlet**, **NeumannZero** (copy), **Extrapolate1** (linear), and **Mirror** (vector with
  component-wise parity). Operators write **faces**; edges/corners follow the last face write.
- Ensure halo depth ≥ widest stencil consumed in ``PostExchange``/``PostBC``.

Output & post-processing
------------------------

Backends
^^^^^^^^
- **XDMF-ORCH/HDF5** — Light metadata in XML; Heavy arrays in HDF5. First-class in ParaView/VisIt; ideal
  for long series.  [XDMF-ORCH]_.
- **CGNS-ORCH** — Standard CFD hierarchy (structured zone; GridCoordinates_t; FlowSolution arrays) for
  maximum tool interoperability.  [CGNS-ORCH]_.

WritePlan & staging
^^^^^^^^^^^^^^^^^^^
The writer synthesizes a **WritePlan** from selected views: logical extents, **byte** strides,
contiguity flags, element size after precision policy, and total byte counts. It chooses between
**direct** writes and **staging** (page-aligned host buffers) and can apply **precision downcast**
(e.g., ``float64 → float32``) during packing. For HDF5 at scale, enable **collective metadata I/O** and
tune chunking/alignment to mitigate metadata hot-spots.  [HDF5-TUNE-ORCH]_.

Error handling & invariants
---------------------------

- Views never free; **single owner** (MemoryManager/driver code) controls lifetime.
- Host allocations meet alignment invariants (≥64 B on x86 is typical).  [INTEL-OPT-ORCH]_
- ``std::aligned_alloc`` sizes are **integral multiples of alignment**.  [ALIGNED-ALLOC-ORCH]_
- Asynchronous transfers originate from **pinned** buffers; MPI requests complete before re-use.  [BPG-ORCH]_, [MPI-NB-ORCH]_

References
----------

.. [BPG-ORCH] **NVIDIA**, *CUDA C++ Best Practices Guide* — coalesced access, pinned+async copies, overlap. https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/

.. [UM-API-ORCH] **NVIDIA**, *CUDA C++ Programming Guide — Unified Memory (prefetch/advise)*. https://docs.nvidia.com/cuda/cuda-c-programming-guide/  :contentReference

.. [MPI-NB-ORCH] **ENCCS**, *Non-blocking MPI for stencil workflows* — overlap pattern & caveats. https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/  

.. [HALO-PATTERN-ORCH] **Pearson (IWAPT’20)**, *Node-Aware Stencil Communication* — pack/send/unpack and face/edge/corner treatment. https://www.carlpearson.net/pdf/20200522_pearson_iwapt.pdf 

.. [XDMF-ORCH] **XDMF-ORCH Project**, *Model & Format* — Light (XML) / Heavy (HDF5). https://www.XDMF-ORCH.org/index.php/XDMF-ORCH_Model_and_Format  

.. [CGNS-ORCH] **CGNS-ORCH SIDS / NASA User Guide**, structured conventions, FlowSolution arrays. https://CGNS-ORCH.org/standard/SIDS/CGNS-ORCH_SIDS.html ; https://ntrs.nasa.gov/api/citations/20010110762/downloads/20010110762.pdf 

.. [HDF5-TUNE-ORCH] **HDF Group**, *Parallel HDF5 performance & collective metadata*. https://support.hdfgroup.org/documentation/hdf5-docs/hdf5_topics/Parallel-HDF5-Performance-Tuning.pdf ; https://support.hdfgroup.org/documentation/hdf5/latest/_par_compr.html 

.. [SPAN-ORCH] **cppreference**, ``std::SPAN-ORCH`` — contiguous, non-owning view. https://en.cppreference.com/w/cpp/container/SPAN-ORCH.html 

.. [KOKKOS-VIEW-ORCH] **Kokkos**, *View — unmanaged/wrapping existing allocations*. https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html  

.. [INTEL-OPT-ORCH] **Intel**, *Optimization Reference Manual* — 64-byte cache lines & alignment guidance. https://cdrdv2-public.intel.com/814198/248966-Optimization-Reference-Manual-V1-049.pdf  

.. [ALIGNED-ALLOC-ORCH] **cppreference**, ``std::aligned_alloc`` — size multiple of alignment. https://en.cppreference.com/w/cpp/memory/c/aligned_alloc  

.. [PLUGINS-ORCH] Examples of C/C++ “registry + dlopen” plugin patterns and discussion. https://ecs.syr.edu/faculty/fawcett/Handouts/CppShortCourse/PLUGINS.htm 
