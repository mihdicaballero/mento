# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import datetime
from importlib.metadata import version

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "mento"
author = "Mehdí Caballero & Juan Pablo Romaris"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

try:  # pragma: no cover
    version = version(project)  # type: ignore
except Exception:  # pragma: no cover
    # we seem to have a local copy not installed without setuptools
    # so the reported version will be unknown
    version = "unknown"  # type: ignore

release = version
this_year = datetime.date.today().year
copyright = f"2025-{this_year}, mento Developers"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "sphinx_copybutton",
    "sphinx.ext.mathjax",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_design",
]

# Hide input and output prompts
nbsphinx_input_prompt = "%.0s"
nbsphinx_output_prompt = "%.0s"

master_doc = "index"

templates_path = ["_templates"]

source_suffix = {
    ".rst": "restructuredtext",
}

exclude_patterns = ["build"]

autodoc_inherit_docstrings = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = [
    "custom.css",  # Add a custom CSS file
]
html_logo = "_static/logo/mento_isotipo_transparente.png"
html_favicon = "_static/logo/mento_isotipo_transparente.png"
html_theme_options = {
    "repository_url": "https://github.com/mihdicaballero/mento",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_issues_button": True,
    "logo_only": True,
    "display_version": False,
}

# Output file base name for HTML help builder.
htmlhelp_basename = "mentodoc"
