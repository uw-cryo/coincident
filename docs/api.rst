.. _reference:

API reference
=============

The API reference provides an overview of public objects, functions, and methods implemented in ``coincident``.

Search
--------

.. currentmodule:: coincident.search

.. autosummary::
   :toctree: generated/

   search
   to_pystac_items
   wesm.load_by_fid
   wesm.get_swath_polygons


Overlaps
--------

.. currentmodule:: coincident.overlaps

.. autosummary::
   :toctree: generated/

   geographic_area
   subset_by_minimum_area
   subset_by_temporal_overlap


Datasets
--------

.. currentmodule:: coincident.datasets

.. autosummary::
   :toctree: generated/

   general.Dataset
   maxar.open_browse
   maxar.plot_browse


IO
--

.. currentmodule:: coincident.io

.. autosummary::
   :toctree: generated/

   sliderule.subset_gedi02a
   sliderule.subset_atl06
   sliderule.process_atl06sr
   sliderule.sample_3dep
   xarray.to_dataset
   xarray.plot_esa_worldcover
