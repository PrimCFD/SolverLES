.. _memory:

Memory & Data-Layout
====================

This page describes the project’s memory subsystem in depth: **Structure-of-Arrays (SoA) fields with linear indexing**, a **single owner** with non-owning **views**, explicit **alignment** and **pools** on CPU/GPU, efficient **host↔device transfers** (pinned + streams), optional **Unified Memory** with prefetch/advice, and scalable **halo exchange** (faces/edges/corners) via non-blocking MPI.

High-level goals
----------------

- Keep inner loops pointer-fast and vectorizable on CPU and coalesced on GPU.  [BPG]_
- Centralize ownership (RAII) and expose zero-overhead views for kernels.  [SPAN]_
- Minimize (de)allocation cost via pooling on device and aligned pools on host.  [CUDA-POOLS]_
- Overlap data movement and compute where practical (streams + non-blocking MPI).  [BPG]_, [MPI-NB]_
- Keep the API policy-agnostic: explicit copies by default; UM as an opt-in.  [UM-API]_

Core concepts
-------------

Structure-of-Arrays (SoA)
^^^^^^^^^^^^^^^^^^^^^^^^^
Each primitive (``rho``, ``ux``, ``uy``, ``uz``, ``p``, …) lives in a separate, contiguous 1-D array. Neighboring threads/lanes touch neighboring elements, enabling:
(1) **CPU SIMD** unit-stride vector loads/stores; (2) **GPU** coalesced global memory transactions.  [INTEL-OPT]_, [BPG]_

Layout & linear indexing
^^^^^^^^^^^^^^^^^^^^^^^^
Data are stored row-major and indexed linearly:

.. code-block:: c++

   // (i,j,k) → idx for dimensions Nx, Ny, Nz (row-major: x varies fastest)
   inline size_t idx(size_t i, size_t j, size_t k,
                     size_t Nx, size_t Ny) noexcept {
     return (k*Ny + j)*Nx + i;
   }

Row-major formulas and tradeoffs are standard; use unit-stride in the innermost loop.  [ROWMAJOR]_

Ownership & views
^^^^^^^^^^^^^^^^^
``MemoryManager`` owns raw buffers and returns non-owning views that carry shape/strides (host) or device pointers for kernels. On CPU, views model ``std::span`` semantics (contiguous, non-owning). For portability, Kokkos “unmanaged views” are an analogous concept.  [SPAN]_, [KOKKOS-VIEW]_

Alignment & allocators
----------------------

Host (CPU)
^^^^^^^^^^
- **Cache lines are 64 B** on modern x86; aligning arrays to ≥64 B avoids line splits and false sharing.  [INTEL-OPT]_
- Use ``std::aligned_alloc`` for explicitly aligned pools; **size must be a multiple of alignment** (C++17 rule).  [ALIGNED-ALLOC]_

Device (GPU)
^^^^^^^^^^^^
- CUDA runtime/driver allocations are **≥256-byte aligned**; some I/O paths (e.g., GDS) require larger (4 KiB) alignment.  [RMM-ALIGN]_
- Prefer **stream-ordered memory pools** (``cudaMallocAsync``/mem pools) to amortize allocation overhead and reduce sync.  [CUDA-POOLS]_

Host↔Device data movement
-------------------------

Pinned (page-locked) memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use **pinned** host buffers with **``cudaMemcpyAsync`` in streams** to overlap copies with kernels and to reach higher PCIe/NVLink bandwidth. Pageable memory falls back to synchronous behavior.  [BPG]_

Policy: explicit mirrors (default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Host owns canonical SoA arrays; device mirrors are created once and reused.
- Transfers: pack halo faces (if needed), enqueue H2D/D2H on dedicated streams, record events, and overlap with compute.  [BPG]_

Policy: Unified Memory (opt-in)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
UM simplifies ownership (single pointer) but still benefits from **prefetch** and **advice** for performance-critical paths:

- ``cudaMemPrefetchAsync(ptr, nbytes, device, stream)`` to stage pages near the next kernel.  
- ``cudaMemAdvise`` (``SetPreferredLocation``, ``SetAccessedBy``, ``SetReadMostly``) to reduce page thrash.  
  [UM-API]_, [UM-BLOG]_, [UM-ORNL]_, [UM-NASA]_

Halo exchange
-------------

Ghost layers & neighborhoods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We maintain a one-cell (configurable) ghost layer around local subdomains and exchange **faces, edges, and corners** (26-neighbor in 3-D) each step for stencil updates.  [MPI-HALO]_

Non-blocking progression & overlap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The exchange uses ``MPI_Irecv/Isend`` + ``MPI_Waitall``; interior compute proceeds while messages progress. Overlap is **implementation-dependent**, but the non-blocking pattern is the standard route to expose concurrency.  [MPI-NB]_

Datatype option (packing-free faces)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Where convenient, use **``MPI_Type_create_subarray``** (or vector/contiguous types) to describe faces/edges directly in memory and avoid manual pack/unpack.  [MPI-SUBARRAY]_

Threading notes
^^^^^^^^^^^^^^^
On hybrid nodes, OpenMP tasks/threads can dedicate a team to halo progress while others compute local cells.  [MPI-OMP]_

Error handling & invariants
---------------------------

- All allocations come from a **single owner**; views never free.  
- Host allocations meet **alignment** invariants (≥64 B); device meets **≥256 B** alignment.  
- Transfers that claim asynchrony **must** originate from **pinned** buffers.  
- MPI requests are completed before buffer reuse.  
- UM mode must prefetch before first-touch kernels in tight loops.

References
----------

.. [BPG] NVIDIA, *CUDA C++ Best Practices Guide*. Coalesced access, pinned memory & async copies with streams; guidance on overlapping copy/compute. https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/  (accessed Aug 25 2025)

.. [CUDA-POOLS] NVIDIA, *CUDA Runtime API — Memory Pools / Stream-Ordered Allocator* (``cudaMallocAsync``, ``cudaMemPool*``). https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY__POOLS.html

.. [RMM-ALIGN] RAPIDS RMM Docs, *Memory Resources* — CUDA allocations are aligned to **at least 256 bytes**; some paths (e.g., GDS) need larger alignment. https://docs.rapids.ai/api/rmm/nightly/librmm_docs/memory_resources/

.. [INTEL-OPT] Intel, *Intel® 64 and IA-32 Architectures Optimization Reference Manual* — cache line is 64 B; unit-stride & alignment guidance. https://cdrdv2-public.intel.com/814198/248966-Optimization-Reference-Manual-V1-049.pdf

.. [ALIGNED-ALLOC] cppreference, ``std::aligned_alloc`` (C++17) — **size must be an integral multiple of alignment**. https://en.cppreference.com/w/cpp/memory/c/aligned_alloc

.. [ROWMAJOR] Wikipedia, *Row- and column-major order* — linear index formulas & row-major background. https://en.wikipedia.org/wiki/Row-_and_column-major_order

.. [SPAN] cppreference, ``std::span`` — non-owning view over a contiguous sequence (analogy for host views). https://en.cppreference.com/w/cpp/container/span.html

.. [KOKKOS-VIEW] Kokkos, *View — Multidimensional array* — unmanaged/wrapping existing allocations. https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html

.. [UM-API] NVIDIA Docs, *CUDA C++ Programming Guide / Runtime API — Unified Memory* (``cudaMemPrefetchAsync``, ``cudaMemAdvise``). https://docs.nvidia.com/cuda/cuda-c-programming-guide/

.. [UM-BLOG] NVIDIA Developer Blog, *Maximizing Unified Memory Performance in CUDA* — when/why to prefetch & advise. https://developer.nvidia.com/blog/maximizing-unified-memory-performance-cuda/

.. [UM-ORNL] ORNL OLCF Training, *CUDA Unified Memory slides* — concise overview & best practices. https://www.olcf.ornl.gov/wp-content/uploads/2019/06/06_Managed_Memory.pdf

.. [UM-NASA] NASA HECC (2025), *Simplifying GPU Programming with Unified Memory*. https://www.nas.nasa.gov/hecc/support/kb/simplifying-gpu-programming-with-unified-memory_703.html

.. [MPI-HALO] SC’24 Poster / arXiv (2025), *Persistent and Partitioned MPI for Stencil Communication* — defines halo exchange (3-D faces/edges/corners). https://arxiv.org/html/2508.13370v1

.. [MPI-NB] ENCCS, *Non-blocking point-to-point — performant stencil workflow* — overlap is implementation-dependent, pattern for correctness. https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/

.. [MPI-SUBARRAY] RookieHPC, *MPI_Type_create_subarray* — using subarray datatypes for strided faces. https://rookiehpc.org/mpi/docs/mpi_type_create_subarray/index.html

.. [MPI-OMP] ENCCS, *MPI and threads in practice* — OpenMP tasking with halo exchange. https://enccs.github.io/intermediate-mpi/mpi-and-threads-pt2/

