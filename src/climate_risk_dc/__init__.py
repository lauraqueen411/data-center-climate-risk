"""climate_risk_dc: physical climate risk assessment of data centers (western USA).

This package provides reusable, tested logic for building an *exposure asset
library* -- a unified, geospatially-consistent representation of data-center
assets and the supporting infrastructure / resources they depend on
(fiber networks, public water supply, watersheds).

Design principles
-----------------
* Reusable logic lives here in ``src/``; notebooks are for exploration only.
* Vector data is handled with :mod:`geopandas`. Gridded climate *hazard*
  layers (added in later project stages) should be handled with :mod:`xarray`.
* Public functions carry type hints and docstrings.

Sub-modules
-----------
* :mod:`climate_risk_dc.config` -- dataset paths and project constants.
* :mod:`climate_risk_dc.geo` -- CRS handling and spatial-join helpers.
* :mod:`climate_risk_dc.exposure` -- per-dataset loaders + library assembler.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
