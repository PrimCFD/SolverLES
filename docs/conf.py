import os
import pathlib
from pathlib import Path as _Path
import re

# Project info 
project = "SolverLES"
author = "RICHARD Anthony"
release = "0.1.0"

# Extensions
extensions = [
    "breathe",
    "exhale",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.graphviz",
    "sphinx.ext.todo",
    "sphinx.ext.napoleon",
    "myst_parser",            # for Markdown in docs/
    "sphinxcontrib.mermaid",  # for Mermaid in docs/ 
]

# MyST: treat ```mermaid as a directive (no Pygments lexer warning)
myst_fence_as_directive = ["mermaid"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Theme 
html_theme = os.environ.get("SPHINX_THEME", "furo")
if (_Path(__file__).parent / '_static').exists():
    html_static_path = ['_static']
else:
    html_static_path = []

# Breathe / Exhale 
breathe_projects = {}
xml_env = os.environ.get("DOXYGEN_XML_DIR")
if xml_env:
    breathe_projects["SolverLES"] = xml_env
breathe_default_project = "SolverLES"

breathe_use_project_refids = True

# Exhale: generate a tidy API index under docs/api/
repo_root = pathlib.Path(__file__).resolve().parents[1]
exhale_args = {
    "containmentFolder": "api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": str(repo_root),
    "createTreeView": True,
}

# Build strictness & nitpicky mode
# Make warnings fail in CI: SPHINX_STRICT=1 ./scripts/build_docs.sh
if os.environ.get("SPHINX_STRICT") == "1":
    nitpicky = True
    nitpick_ignore = [
        ("cpp:identifier", "core"),
        ("cpp:identifier", "core::detail"),
        ("cpp:identifier", "uint8_t"),
        ("cpp:identifier", "MPI_Comm"),
        ("cpp:identifier", "subroutine"),
    ]
    # Future-proof a bit without being too broad.
    nitpick_ignore_regex = [
        (r"cpp:identifier", r"^(std::)?u?int(8|16|32|64)_t$"),
        (r"cpp:identifier", r"^MPI_.*$"),
        # Exhale/Breathe sometimes emit unqualified identifiers for top-level namespaces
        (r"cpp:identifier", r"^(mesh|io|plugin)(::.*)?$"),
        (r"cpp:identifier", r"^core(::.*)?$"),
        (r"cpp:identifier", r"^core::memory(::.*)?$"),
    ]


# Autosectionlabel so :ref: works across pages without collisions
autosectionlabel_prefix_document = True

# Enable TODOs when desired:
todo_include_todos = bool(int(os.environ.get("SPHINX_TODOS", "0")))

# Graphviz defaults
graphviz_output_format = "svg"

# -- Exhale + Furo: silence 'contents' warning on generated pages ----------
def _tag_contents_for_furo(app, docname, source):
    # Only touch Exhale-generated pages (anything under api/)
    if not docname.startswith("api/"):
        return
    text = source[0]
    # Skip if already present (idempotent)
    if "this-will-duplicate-information-and-it-is-still-useful-here" in text:
        source[0] = text
        return
    # Add the special class Furo looks for
    text = re.sub(
        r"(^\.\. contents::[^\n]*\n(?:\s*:\w+:.*\n)*)",
        r"\1   :class: this-will-duplicate-information-and-it-is-still-useful-here\n",
        text,
        flags=re.M,
    )
    source[0] = text

def setup(app):
    app.connect("source-read", _tag_contents_for_furo)
