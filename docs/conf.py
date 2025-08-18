import os

extensions = ["breathe", "exhale"]

breathe_projects = {
    "SolverLES": os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../build-docs/docs/_build/doxygen/xml")
    )
}
breathe_default_project = "SolverLES"

exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": "../src",
    "createTreeView": True,
}

project = "SolverLES"
author = "RICHARD Anthony"
release = "0.0.0"
