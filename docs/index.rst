FVM–LES–Plasma Solver
=====================

.. rst-class:: lead

Open, modular **finite‑volume LES solver for thermal–plasma flows** focused on fast iteration on new SGS models,
with a clean docs site built for the `Furo`_ theme.

|docs| |license| |linux| |docsci| |style|

.. admonition:: What is this project?
   :class: tip

   A research‑oriented solver to explore **large‑eddy simulation (LES)** models for thermal plasma jets, designed to scale on modern CPU/GPU machines while keeping a modular plugin architecture.

----

Quick links
-----------

- :ref:`User guide <user-guide>` — how to run examples, inputs, and workflow.
- :ref:`Developer guide <developer-guide>` — coding standards, build scripts, CI.
- :ref:`API reference <api-reference>` — C++/Fortran symbols bridged via Breathe/Exhale.
- :ref:`Documentation guide <documentation-guide>` — Guide about documentation structure, tools and build.

Get started
-----------

**Fast path (CPU, Release)**

.. code-block:: bash

   # Configure + build (Release by default) into ./build
   ./scripts/build.sh

   # Run the hello-mesh example
   ./build/bin/solver examples/hello_mesh.yaml

**Docs build & local preview**

.. code-block:: bash

   # Build Doxygen XML + Sphinx HTML, then serve locally
   ./scripts/build_docs.sh --serve --open

.. note::
   CI builds the docs and publishes to GitHub Pages from ``main``.
   Local builds write to ``build-docs/docs/html``. See the :ref:`Developer guide <developer-guide>` for details.

Highlights
----------

- **Modular architecture:** hot‑swappable plugins (flux, SGS, time integration).
- **Mixed‑language core:** modern C++ orchestration with shared Fortran kernels.
- **Scalable:** MPI‑aware and GPU‑ready math kernels.
- **Reproducible builds:** helper scripts for offline vendor cache & deterministic CI.
- **First‑class docs:** Furo theme, MyST, Doxygen → Breathe → Exhale pipeline.

.. note::
   See the :ref:`Developer guide <developer-guide>` for details.

Project navigation
------------------

.. _user-guide:

User guide
~~~~~~~~~~

.. toctree::
   :caption: Guide for users / SGS developers
   :maxdepth: 2

   user/index

.. _developer-guide:

Developer guide
~~~~~~~~~~~~~~~

.. toctree::
   :caption: Guide for core developers
   :maxdepth: 2

   developer/index

.. _api-reference:

API reference
~~~~~~~~~~~~~

.. toctree::
   :caption: Full code API reference
   :maxdepth: 2

   api/library_root
   api/unabridged_orphan

.. _documentation-guide:

Documentation guide
~~~~~~~~~~~~~~~~~~~

.. toctree::
   :caption: Guide about documentation structure, tools and build
   :maxdepth: 1

   README

License
-------

Distributed under the **MIT License**. See :file:`LICENSE`.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`

.. |docs| image:: https://img.shields.io/badge/docs-online-blue
   :target: https://primcfd.github.io/SolverLES/
   :alt: Documentation

.. |license| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: MIT License

.. |linux| image:: https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml/badge.svg?branch=main
   :target: https://github.com/PrimCFD/SolverLES/actions/workflows/linux.yml
   :alt: Linux CI

.. |docsci| image:: https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml/badge.svg?branch=main
   :target: https://github.com/PrimCFD/SolverLES/actions/workflows/docs.yml
   :alt: Docs CI

.. |style| image:: https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml/badge.svg?branch=main
   :target: https://github.com/PrimCFD/SolverLES/actions/workflows/style.yml
   :alt: Style CI

.. _Furo: https://pradyunsg.me/furo/
