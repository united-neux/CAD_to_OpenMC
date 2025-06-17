# Configuration file for the Sphinx documentation builder.
#
# -- Project information -----------------------------------------------------
import os, sys

project = "CAD_to_OpenMC Project"
copyright = "2025, CAD_to_OpenMC-team"
author = "ebknudsen the Docs core team"

sys.path.insert(0, os.path.abspath('../../src/'))

master_doc = 'index'

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "myst_parser",
]

bibtex_bibfiles = ["bibliography.bib"]

myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
]
intersphinx_mapping = {
    "rtd": ("https://docs.readthedocs.io/en/stable/", None),
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

source_suffix = {
        ".rst" : "restructuredtext",
        ".md" : "markdown"
}


numfig = True

# -- Options for EPUB output
epub_show_urls = "footnote"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ["_static"]
