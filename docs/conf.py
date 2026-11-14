import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "miPyRNA"
copyright = "2026, Naveen Duhan"
author = "Naveen Duhan"
release = "1.0.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
autodoc_mock_imports = ["tensorflow", "keras", "RNA", "forgi", "rpy2", "pyfaidx", "HTSeq", "waiting"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
