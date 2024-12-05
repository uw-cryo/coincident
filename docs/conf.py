from __future__ import annotations

import importlib.metadata

project = "coincident"
copyright = "2024, Scott Henderson"
author = "Scott Henderson"
version = release = importlib.metadata.version("coincident")

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

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/uw-cryo/coincident",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}
html_title = "Coincident"
html_context = {
    "github_user": "uw-cryo",
    "github_repo": "coincident",
    "github_version": "main",
    "doc_path": "docs",
}

# remove sidebar, see GH issue #82
html_css_files = [
    "custom.css",
]

html_static_path = ["_static"]

# NOTE: consider adding back in once for distinct sections (user guide, examples, API reference)
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#primary-sidebar-left
# html_sidebars = {"**": []}

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
nb_execution_excludepatterns = ["sliderule.ipynb", "contextual_data.ipynb"]
