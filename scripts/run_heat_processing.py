"""Deliverable 1: Heat index processing for Oregon data-center risk workflow.

Assumptions log
---------------
- Heat-only scope for this run; flood and wildfire are deferred.
- LOCA2 tasmin/tasmax are loaded from NetCDF and converted explicitly from
  Kelvin to Celsius before any index calculations.
- Outputs are period-mean annual indices for two windows:
  historical (1950-1979) and SSP585 future (2070-2099).
- Common CRS for this phase is EPSG:4326.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

from climate_risk_dc.climate.heat_indices import (
    HeatIndexConfig,
    clip_to_bbox,
    compute_period_heat_indices,
    mask_to_polygon,
    threshold_frequency_summary,
)


@dataclass(frozen=True)
class RunConfig:
    """Configuration for heat processing run."""

    model: str = "ACCESS-CM2"
    member: str = "r1i1p1f1"
    scenario_hist: str = "historical"
    scenario_fut: str = "ssp585"
    hist_start: int = 1950
    hist_end: int = 1979
    fut_start: int = 2070
    fut_end: int = 2099
    threshold_c: float = 35.0
    cdd_base_c: float = 18.3
    min_duration_days: int = 3
    lon_min: float = -125.0
    lon_max: float = -116.0
    lat_min: float = 42.0
    lat_max: float = 46.5


def _loca2_file(
    loca2_root: Path,
    model: str,
    member: str,
    scenario: str,
    variable: str,
) -> Path:
    """Resolve the expected LOCA2 NetCDF path for the n_west regional split."""
    year_part = "1950-2014" if scenario == "historical" else "2015-2100"
    # Some future files are split by shorter windows; we first try canonical full file.
    candidate = (
        loca2_root
        / model
        / "n_west"
        / "0p0625deg"
        / member
        / scenario
        / variable
        / f"{variable}.{model}.{scenario}.{member}.{year_part}.LOCA_16thdeg_v20220413.n_west.nc"
    )
    if candidate.exists():
        return candidate

    # Fallback: select a daily file in the variable directory.
    var_dir = loca2_root / model / "n_west" / "0p0625deg" / member / scenario / variable
    all_matches = sorted(var_dir.glob(f"{variable}.{model}.{scenario}.{member}*.n_west.nc"))
    daily_matches = [
        p for p in all_matches if ".monthly." not in p.name and ".yearly." not in p.name
    ]
    if not daily_matches:
        raise FileNotFoundError(f"No LOCA2 file found in {var_dir}")

    def _year_span(path: Path) -> tuple[int, int] | None:
        m = re.search(r"\.(\d{4})-(\d{4})\.LOCA", path.name)
        if not m:
            return None
        return int(m.group(1)), int(m.group(2))

    # Prefer files that contain the later years we need for SSP585 and full span for historical.
    if scenario != "historical":
        candidates = []
        for p in daily_matches:
            span = _year_span(p)
            if span is None:
                continue
            start, end = span
            if end >= 2099:
                candidates.append((start, end, p))
        if candidates:
            candidates.sort(key=lambda x: (x[0], x[1]))
            return candidates[-1][2]

    return daily_matches[0]


def _oregon_boundary(watershed_gdb: Path, layer: str = "WBDHU8") -> gpd.GeoDataFrame:
    """Build Oregon boundary from WBD layer dissolve."""
    wbd = gpd.read_file(watershed_gdb, layer=layer)
    if wbd.crs is None:
        wbd = wbd.set_crs("EPSG:4326")
    wbd = wbd.to_crs("EPSG:4326")
    if "states" in wbd.columns:
        wbd = wbd[wbd["states"].fillna("").str.contains("OR")].copy()
    if wbd.empty:
        raise ValueError("Could not build Oregon boundary from WBD data.")
    boundary = gpd.GeoDataFrame(geometry=[wbd.union_all()], crs="EPSG:4326")
    return boundary


def _da_to_geotiff(da_2d, out_path: Path) -> None:
    """Write 2D DataArray with dims (lat, lon) to GeoTIFF in EPSG:4326."""
    if tuple(da_2d.dims) != ("lat", "lon"):
        raise ValueError(f"Expected dims ('lat', 'lon'), got {da_2d.dims}")

    lat = da_2d["lat"].values
    lon = da_2d["lon"].values
    data = da_2d.values.astype("float32")

    if lat[0] < lat[-1]:
        lat = lat[::-1]
        data = data[::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        data = data[:, ::-1]

    yres = float(np.abs(lat[1] - lat[0]))
    xres = float(np.abs(lon[1] - lon[0]))
    transform = from_origin(float(lon.min()) - xres / 2, float(lat.max()) + yres / 2, xres, yres)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        out_path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=np.nan,
    ) as dst:
        dst.write(data, 1)


def _plot_quicklook(index_da, boundary_gdf, dc_gdf, title: str, out_path: Path) -> None:
    """Create quicklook map for an index with Oregon boundary and data centers."""
    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
    index_da.plot(ax=ax, cmap="inferno", cbar_kwargs={"shrink": 0.8})
    boundary_gdf.boundary.plot(ax=ax, color="black", linewidth=1.2)
    dc_gdf.plot(ax=ax, color="cyan", edgecolor="black", markersize=45, alpha=0.9)
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def _load_oregon_data_centers(csv_path: Path) -> gpd.GeoDataFrame:
    df = pd.read_csv(csv_path)
    df = df[df["state_abb"] == "OR"].copy()
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326")
    return gdf


def _select_temp_var(ds, base_name: str):
    if base_name in ds.data_vars:
        return ds[base_name]
    # Some LOCA2 products use suffixed names such as tasmax_tavg.
    for name in ds.data_vars:
        if name.lower().startswith(base_name.lower()):
            return ds[name]
    if len(ds.data_vars) == 1:
        return next(iter(ds.data_vars.values()))
    raise KeyError(f"Could not resolve variable '{base_name}' from {list(ds.data_vars)}")


def _first_finite_value(da) -> float:
    arr = np.asarray(da.values, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("No finite values found for temperature sanity check.")
    return float(finite[0])


def _log_unit_sanity(tasmax_file: Path, tasmin_file: Path) -> None:
    import xarray as xr

    dsx = xr.open_dataset(tasmax_file)
    dsn = xr.open_dataset(tasmin_file)
    tasmax_k = _select_temp_var(dsx, "tasmax")
    tasmin_k = _select_temp_var(dsn, "tasmin")

    sample_max = _first_finite_value(tasmax_k)
    sample_min = _first_finite_value(tasmin_k)
    print(f"Kelvin sample tasmax={sample_max:.3f}, tasmin={sample_min:.3f}")
    print(f"Celsius sample tasmax={sample_max - 273.15:.3f}, tasmin={sample_min - 273.15:.3f}")

    for v in (sample_max, sample_min):
        if v < 180.0 or v > 330.0:
            raise ValueError(
                f"Temperature sample {v:.3f} is outside plausible Kelvin bounds [180, 330]."
            )


def run(args: argparse.Namespace) -> None:
    cfg = RunConfig()
    hi_cfg = HeatIndexConfig(
        cdd_base_c=cfg.cdd_base_c,
        threshold_c=cfg.threshold_c,
        min_duration_days=cfg.min_duration_days,
    )

    root = Path(args.repo_root).resolve()
    outputs = root / "outputs" / "month2_heat"
    raster_dir = outputs / "rasters"
    map_dir = outputs / "maps"
    summary_dir = outputs / "summaries"

    watershed_gdb = root / "data" / "watershed_boundaries" / "WBD_National_GDB.gdb"
    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    loca2_root = Path(args.loca2_root)

    boundary = _oregon_boundary(watershed_gdb)
    boundary_poly = boundary.geometry.iloc[0]
    dc_gdf = _load_oregon_data_centers(dc_csv)

    periods = [
        ("historical_1950_1979", cfg.scenario_hist, cfg.hist_start, cfg.hist_end),
        ("ssp585_2070_2099", cfg.scenario_fut, cfg.fut_start, cfg.fut_end),
    ]

    for label, scenario, start_year, end_year in periods:
        tasmax_file = _loca2_file(loca2_root, cfg.model, cfg.member, scenario, "tasmax")
        tasmin_file = _loca2_file(loca2_root, cfg.model, cfg.member, scenario, "tasmin")
        print(f"Using {label}:\n  tasmax={tasmax_file}\n  tasmin={tasmin_file}")
        _log_unit_sanity(tasmax_file, tasmin_file)

        idx = compute_period_heat_indices(
            tasmax_path=tasmax_file,
            tasmin_path=tasmin_file,
            start_year=start_year,
            end_year=end_year,
            config=hi_cfg,
        )

        tmax_c_clipped = clip_to_bbox(idx["tmax_c"], cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        freq = threshold_frequency_summary(tmax_c_clipped, cfg.threshold_c)
        summary_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([freq]).to_csv(summary_dir / f"threshold_frequency_{label}.csv", index=False)

        print(
            "Threshold-days sanity summary "
            f"(p05={freq['p05']:.2f}, p50={freq['p50']:.2f}, p95={freq['p95']:.2f}, max={freq['max']:.2f})"
        )

        for key in ("cdd", "tmax_annual_max", "tmax_threshold_days", "heatwave_events"):
            da = clip_to_bbox(idx[key], cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
            da = mask_to_polygon(da, boundary_poly)

            tif_path = raster_dir / f"{key}_{label}.tif"
            png_path = map_dir / f"{key}_{label}.png"

            _da_to_geotiff(da, tif_path)
            _plot_quicklook(
                da,
                boundary,
                dc_gdf,
                title=f"{key} ({label.replace('_', ' ')})",
                out_path=png_path,
            )

    manifest = pd.DataFrame(
        [
            {
                "model": cfg.model,
                "member": cfg.member,
                "historical_window": f"{cfg.hist_start}-{cfg.hist_end}",
                "future_window": f"{cfg.fut_start}-{cfg.fut_end}",
                "threshold_c": cfg.threshold_c,
                "cdd_base_c": cfg.cdd_base_c,
                "min_duration_days": cfg.min_duration_days,
                "crs": "EPSG:4326",
            }
        ]
    )
    manifest.to_csv(outputs / "run_manifest.csv", index=False)
    print(f"Wrote outputs to {outputs}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LOCA2 heat processing for Oregon.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root path (default: current directory)",
    )
    parser.add_argument(
        "--loca2-root",
        default="/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2/CONUS_regions_split",
        help="Root path for LOCA2 CONUS_regions_split files",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
