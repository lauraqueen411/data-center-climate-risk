"""Project configuration: dataset locations and shared constants.

Paths are resolved relative to the repository root so the package works
regardless of the current working directory (e.g. when imported from a
notebook in ``notebooks/``). If you relocate the data, update
``config/datasets.yml`` rather than editing code.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

# ---------------------------------------------------------------------------
# Repository layout
# ---------------------------------------------------------------------------
# This file lives at <repo>/src/climate_risk_dc/config.py, so the repo root is
# three parents up.
REPO_ROOT: Path = Path(__file__).resolve().parents[2]
DATA_DIR: Path = REPO_ROOT / "data"
OUTPUTS_DIR: Path = REPO_ROOT / "outputs"
CONFIG_DIR: Path = REPO_ROOT / "config"

# ---------------------------------------------------------------------------
# Coordinate reference systems
# ---------------------------------------------------------------------------
#: Geographic CRS used for input lon/lat coordinates and web mapping (WGS84).
CRS_GEOGRAPHIC: str = "EPSG:4326"

#: Equal-area CRS for distance / area computations in the CONUS.
#: NAD83 / CONUS Albers -- also the native CRS of the water-service-area layer.
CRS_EQUAL_AREA: str = "EPSG:5070"

#: Region of interest for the first project stage (western USA).
#: State abbreviations used to subset national datasets.
WESTERN_STATES: tuple[str, ...] = (
    "WA", "OR", "CA", "ID", "NV", "MT", "WY", "UT", "AZ", "CO", "NM",
)


@dataclass(frozen=True)
class DatasetPaths:
    """Resolved on-disk locations of the raw exposure datasets.

    Attributes
    ----------
    data_centers_gpkg:
        GeoPackage of the IM3/PNNL open-source data-center atlas. Contains
        ``point``, ``building`` and ``campus`` layers.
    data_centers_csv:
        Tabular version of the same atlas (lon/lat columns).
    fiber_csv:
        FCC Broadband Data Collection "fiber to the premises" location file.
    water_service_areas_shp:
        USEPA/USGS water-service-area polygons (``WSA_v1``).
    watershed_gdb_zip:
        Zipped USGS Watershed Boundary Dataset (national GDB, HUC polygons).
    """

    data_centers_gpkg: Path
    data_centers_csv: Path
    fiber_csv: Path
    water_service_areas_shp: Path
    watershed_gdb_zip: Path


def _default_paths() -> DatasetPaths:
    """Return the default dataset paths based on the current repo contents."""
    return DatasetPaths(
        data_centers_gpkg=DATA_DIR
        / "data_centers_im3_pnnl"
        / "im3_open_source_data_center_atlas_v2026.02.09.gpkg",
        data_centers_csv=DATA_DIR
        / "data_centers_im3_pnnl"
        / "im3_open_source_data_center_atlas_v2026.02.09.csv",
        fiber_csv=DATA_DIR
        / "fiber"
        / "bdc_41_FibertothePremises_fixed_broadband_D25_08jul2026.csv",
        water_service_areas_shp=DATA_DIR
        / "public_water_sources"
        / "WSA_v1"
        / "WSA_v1.shp",
        watershed_gdb_zip=DATA_DIR
        / "watershed_boundaries"
        / "WBD_National_GDB.zip",
    )


def load_dataset_paths(config_file: Path | None = None) -> DatasetPaths:
    """Load dataset paths, optionally overriding defaults from a YAML file.

    Parameters
    ----------
    config_file:
        Path to a YAML file with a top-level ``datasets`` mapping whose keys
        match :class:`DatasetPaths` fields and whose values are paths relative
        to the repository root. If ``None`` (default) and
        ``config/datasets.yml`` exists, that file is used; otherwise the
        built-in defaults are returned.

    Returns
    -------
    DatasetPaths
        Fully resolved (absolute) dataset paths.
    """
    if config_file is None:
        candidate = CONFIG_DIR / "datasets.yml"
        config_file = candidate if candidate.exists() else None

    if config_file is None:
        return _default_paths()

    with open(config_file, "r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    datasets: dict[str, str] = raw.get("datasets", {})
    defaults = _default_paths()

    def _resolve(field: str) -> Path:
        if field in datasets:
            return (REPO_ROOT / datasets[field]).resolve()
        return getattr(defaults, field)

    return DatasetPaths(
        data_centers_gpkg=_resolve("data_centers_gpkg"),
        data_centers_csv=_resolve("data_centers_csv"),
        fiber_csv=_resolve("fiber_csv"),
        water_service_areas_shp=_resolve("water_service_areas_shp"),
        watershed_gdb_zip=_resolve("watershed_gdb_zip"),
    )
