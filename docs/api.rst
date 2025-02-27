.. _reference:

API reference
=============

The API reference provides an overview of public objects, functions, and methods implemented in ``coincident``.

Search
--------

.. currentmodule:: coincident.search

.. autosummary::
   :nosignatures:
   :toctree: generated/

   search
   cascading_search
   to_pystac_items
   to_geopandas
   wesm.read_wesm_csv
   wesm.load_by_fid
   wesm.get_swath_polygons


Overlaps
--------

.. currentmodule:: coincident.overlaps

.. autosummary::
   :nosignatures:
   :toctree: generated/

   geographic_area
   subset_by_minimum_area
   subset_by_temporal_overlap
   subset_by_maximum_duration

IO
--

.. currentmodule:: coincident.io

.. autosummary::
   :nosignatures:
   :toctree: generated/

   sliderule.subset_gedi02a
   sliderule.subset_atl06
   sliderule.process_atl06sr
   sliderule.sample_3dep
   xarray.to_dataset
   xarray.open_maxar_browse
   download.download_item

Plot
----

.. currentmodule:: coincident.plot

.. autosummary::
   :nosignatures:
   :toctree: generated/

   plot_maxar_browse
   plot_esa_worldcover
   plot_dem
   compare_dems
   plot_altimeter_points
   boxplot_terrain_diff
   hist_esa

Datasets
--------

.. currentmodule:: coincident.datasets

.. autosummary::
   :nosignatures:
   :toctree: generated/

   general.Dataset
