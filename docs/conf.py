from __future__ import annotations

import importlib.metadata

project = "coincident"
copyright = "2025, UW TACO Lab"
author = "Scott Henderson"
version = release = importlib.metadata.version("coincident")

# NOTE: it seems order of extensions matters if one depends on another
extensions = [
    # "myst_parser",
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_design",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".md": "myst-nb",
}
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Other API docs using this theme for inspiration:
# https://virtualizarr.readthedocs.io

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "announcement": "Welcome! Coincident is in early development. APIs may change significantly",
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/uw-cryo/coincident",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
    "show_version_warning_banner": True,
    "footer_center": ["last-updated"],
}
html_title = "Coincident"
html_context = {
    "github_user": "uw-cryo",
    "github_repo": "coincident",
    "github_version": "main",
    "doc_path": "docs",
}
html_show_sourcelink = False

# remove sidebar, see GH issue #82
html_css_files = [
    "custom.css",
]

html_static_path = ["_static"]

# NOTE: consider adding back in once for distinct sections (user guide, examples, API reference)
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#primary-sidebar-left
# html_sidebars = {"**": []}
# Remove primary sidebar from just API autodoc page
html_sidebars = {"api": []}

myst_enable_extensions = [
    "colon_fence",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "geopandas": ("https://geopandas.org/en/latest", None),
    "pystac_client": ("https://pystac-client.readthedocs.io/en/stable", None),
    "pystac": ("https://pystac.readthedocs.io/en/stable", None),
    "stac-asset": ("https://stac-asset.readthedocs.io/en/latest", None),
    "xarray": ("https://docs.xarray.dev/en/stable", None),
}

nitpick_ignore = [
    ("py:class", "_io.StringIO"),
    ("py:class", "_io.BytesIO"),
]

always_document_param_types = True
# autodoc_typehints = "none"
nb_execution_mode = "auto"  # off, on
nb_execution_show_tb = True
nb_execution_timeout = 90
nb_execution_allow_errors = False
# nb_execution_excludepatterns = ["elevation_plotting.ipynb"]
