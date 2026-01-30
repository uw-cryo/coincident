from __future__ import annotations

import importlib.metadata
import warnings

from sphinx.deprecation import RemovedInSphinx10Warning

# Ignore specific RemovedInSphinx warnings
warnings.filterwarnings("ignore", category=RemovedInSphinx10Warning)

project = "coincident"
copyright = "2026, UW TACO Lab"
author = "Scott Henderson"
version = release = importlib.metadata.version("coincident")

# NOTE: it seems order of extensions matters if one depends on another
extensions = [
    # "myst_parser", # myst-nb includes myst-parser
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    # "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_design",
]


# Code block syntax highlighting
pygments_style = "sphinx"

# Show simplified type names
typehints_fully_qualified = False

autodoc_type_aliases = {
    "gpd.GeoDataFrame": "geopandas.GeoDataFrame",
    "gpd.GeoSeries": "geopandas.GeoSeries",
    "xr.DataArray": "xarray.DataArray",
    "xr.Dataset": "xarray.Dataset",
    "pd.Series": "pandas.Series",
    "pd.Timestamp": "pandas.Timestamp",
    "pd.DataFrame": "pandas.DataFrame",
    "pystac.item_collection.ItemCollection": "pystac.ItemCollection",
    "pystac_client.client.Client": "pystac_client.Client",
    "pystac_client.item_search.ItemSearch": "pystac_client.ItemSearch",
    "arro3.core.Table": "arro3.core.Table",  # don't abbreviate this one
}

autosummary_generate = True
autosummary_generate_overwrite = False  # Don't regenerate existing files

# Saves RAM by not importing and inspecting all dependencies!
autodoc_mock_imports = [
    # Heavy geospatial libraries with C/C++ extensions
    "rasterio",  # GDAL wrapper, very memory-intensive
    "osgeo",  # GDAL Python bindings
    "gdal",
    "gdal_array",
    "pyogrio",  # Fast geospatial I/O
    # Heavy data processing libraries
    "odc.stac",  # Open Data Cube STAC loader
    "xarray",  # Large array operations
    "rioxarray",  # Rasterio + xarray
    # Cloud/AWS libraries
    "boto3",
    "botocore",
    "aiobotocore",
    "s3fs",
    # Heavy geospatial analysis
    # "geopandas",  # mocking breaks type hint hyperlinks
    # "shapely",  # mocking breaks type hint hyperlinks
    "pyproj",
    # Specialized data access
    "sliderule",
    "maxar_platform",
    # STAC libraries (less critical but can add up)
    # NOTE: pystac and pystac_client are needed for intersphinx type linking
    # "pystac",
    # "pystac_client",
    "rustac",
    "stac_asset",
    # Arrow ecosystem
    # "arro3.core", # mocking breaks type hint hyperlinks
    # "pyarrow", # mocking breaks type hint hyperlinks
    # Visualization (probably safe to mock for docs)
    "matplotlib",
    "contextily",
    "matplotlib_scalebar",
    "xyzservices",
    # Misc heavy dependencies
    "scipy",
    "tqdm",
    "requests",
]


autodoc_preserve_defaults = True
autodoc_typehints = "description"  # Show type hints in docstring only, not in signature
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "undoc-members": False,  # Skip undocumented members
    "private-members": False,
    "special-members": False,
    "inherited-members": False,  # Don't document inherited members
    "show-inheritance": False,
}

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".md": "myst-nb",
}
exclude_patterns = [
    "_build",
    # "datasets", # UNCOMMENT to temporarily disable executing notebooks
    # "examples",
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
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "geopandas": ("https://geopandas.org/en/latest", None),
    "pystac_client": ("https://pystac-client.readthedocs.io/en/stable", None),
    "pystac": ("https://pystac.readthedocs.io/en/stable", None),
    "stac-asset": ("https://stac-asset.readthedocs.io/en/latest", None),
    "xarray": ("https://docs.xarray.dev/en/stable", None),
    "arro3": ("https://kylebarron.dev/arro3/latest", None),
}

nitpick_ignore = [
    ("py:class", "_io.StringIO"),
    ("py:class", "_io.BytesIO"),
]

# always_document_param_types = True
# autodoc_typehints = "none"
nb_execution_mode = "auto"  # off, on
nb_execution_show_tb = True
nb_execution_timeout = 90
nb_execution_allow_errors = True
# NOTE: skip execution of specific notebooks that require large downloads or long processing times
nb_execution_excludepatterns = ["noaa-coastal-lidar.ipynb"]
