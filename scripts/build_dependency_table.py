"""Deliverable 2: Build Oregon data-center dependency table.

Assumptions log
---------------
- Power linkage uses a moderate proxy: facility -> nearby substation(s) ->
  connected transmission endpoint substations -> generators within a configured
  radius of those connected substations.
- This is a spatial proxy, not contractual electrical topology.
- Water source type is set to 'unknown' in this phase because the local WSA
  schema does not expose a direct surface/groundwater field.
- Flood and wildfire are deferred; placeholder columns are emitted as null.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import h3
import networkx as nx
import numpy as np
import pandas as pd
import rasterio

from climate_risk_dc.geo import load_oregon_data_centers


@dataclass(frozen=True)
class TableConfig:
    substation_service_radius_miles: float = 15.0
    generation_proxy_radius_miles: float = 50.0
    transmission_min_kv: float = 100.0
    heat_threshold_c: float = 35.0
    watershed_level: str = "WBDHU8"
    max_candidate_substations: int = 3
    include_transmission_one_hop: bool = False

    @property
    def substation_service_radius_m(self) -> float:
        return self.substation_service_radius_miles * 1609.344

    @property
    def generation_proxy_radius_m(self) -> float:
        return self.generation_proxy_radius_miles * 1609.344


def _col_pick(df: pd.DataFrame, candidates: list[str]) -> str | None:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def _load_generators(path: Path) -> gpd.GeoDataFrame:
    xls = pd.ExcelFile(path)
    frames: list[pd.DataFrame] = []

    for sheet in xls.sheet_names:
        raw = pd.read_excel(path, sheet_name=sheet)
        if raw.empty:
            continue

        # EIA workbook has an extra row containing display labels.
        if len(raw) > 2:
            raw.columns = raw.iloc[1].fillna("").astype(str)
            df = raw.iloc[2:].copy()
        else:
            df = raw.copy()

        frames.append(df)

    if not frames:
        raise ValueError("No generator sheets with tabular content were loaded.")

    df = pd.concat(frames, ignore_index=True)

    lat_col = _col_pick(df, ["Latitude", "Plant Latitude"])
    lon_col = _col_pick(df, ["Longitude", "Plant Longitude"])
    if lat_col is None or lon_col is None:
        raise KeyError("Could not find Latitude/Longitude columns in generator workbook.")

    df[lat_col] = pd.to_numeric(df[lat_col], errors="coerce")
    df[lon_col] = pd.to_numeric(df[lon_col], errors="coerce")
    df = df.dropna(subset=[lat_col, lon_col]).copy()

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[lon_col], df[lat_col]), crs="EPSG:4326")

    # Standardized optional fields for downstream summaries.
    fuel_col = _col_pick(gdf, ["Energy Source 1", "Fuel1", "Fuel", "Technology"])
    prime_col = _col_pick(gdf, ["Prime Mover", "PrimeMover", "Prime Mover Code"])
    capacity_col = _col_pick(gdf, ["Nameplate Capacity (MW)", "Nameplate Capacity", "Summer Capacity (MW)"])
    plant_col = _col_pick(gdf, ["Plant Name", "Plant Name and Unit", "Plant"])

    gdf["fuel_type"] = gdf[fuel_col].astype(str) if fuel_col else "unknown"
    gdf["prime_mover"] = gdf[prime_col].astype(str) if prime_col else "unknown"
    gdf["capacity_mw"] = pd.to_numeric(gdf[capacity_col], errors="coerce") if capacity_col else np.nan
    gdf["generator_name"] = gdf[plant_col].astype(str) if plant_col else "unknown"

    # Collapse unit-level rows to a plant-level proxy to avoid over-counting
    # independent generation sources at a single physical site.
    gdf["gen_lat"] = gdf.geometry.y.round(5)
    gdf["gen_lon"] = gdf.geometry.x.round(5)
    gdf["plant_key"] = (
        gdf["generator_name"].fillna("unknown").astype(str)
        + "|"
        + gdf["gen_lon"].astype(str)
        + "|"
        + gdf["gen_lat"].astype(str)
    )

    def _pick_first_non_unknown(values: pd.Series) -> str:
        vals = [str(v) for v in values.dropna().tolist() if str(v).strip() and str(v).lower() != "unknown"]
        return vals[0] if vals else "unknown"

    agg = (
        gdf.groupby("plant_key", as_index=False)
        .agg(
            generator_name=("generator_name", "first"),
            fuel_type=("fuel_type", _pick_first_non_unknown),
            prime_mover=("prime_mover", _pick_first_non_unknown),
            capacity_mw=("capacity_mw", "sum"),
            gen_lon=("gen_lon", "first"),
            gen_lat=("gen_lat", "first"),
            generator_unit_count=("plant_key", "size"),
        )
        .copy()
    )
    agg["generator_id"] = np.arange(len(agg), dtype=int)
    plant_gdf = gpd.GeoDataFrame(agg, geometry=gpd.points_from_xy(agg["gen_lon"], agg["gen_lat"]), crs="EPSG:4326")
    return plant_gdf


def _load_power_layers(substations_path: Path, transmission_path: Path, min_kv: float) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    subst = gpd.read_file(substations_path)
    lines = gpd.read_file(transmission_path)
    if subst.crs is None:
        subst = subst.set_crs("EPSG:4326")
    if lines.crs is None:
        lines = lines.set_crs("EPSG:4326")
    subst = subst.to_crs("EPSG:4326")
    lines = lines.to_crs("EPSG:4326")

    kv_col = _col_pick(lines, ["VOLTAGE", "VOLT_CLASS", "VOLTCLASS"])
    if kv_col is not None:
        kv = pd.to_numeric(lines[kv_col], errors="coerce")
        lines = lines[(kv.isna()) | (kv >= min_kv)].copy()

    return subst, lines


def _build_transmission_graph(lines: gpd.GeoDataFrame) -> nx.Graph:
    sub1_col = _col_pick(lines, ["SUB_1", "sub_1"])
    sub2_col = _col_pick(lines, ["SUB_2", "sub_2"])
    if sub1_col is None or sub2_col is None:
        raise KeyError("Transmission layer missing SUB_1/SUB_2 fields for topology proxy.")

    graph = nx.Graph()
    valid = lines[[sub1_col, sub2_col]].dropna()
    for _, row in valid.iterrows():
        a = str(row[sub1_col]).strip()
        b = str(row[sub2_col]).strip()
        if a and b:
            graph.add_edge(a, b)
    return graph


def _substation_id_column(subst: gpd.GeoDataFrame) -> str:
    col = _col_pick(subst, ["OBJECTID", "OBJECTID_1", "ID", "SUB_ID", "NAME"])
    if col is None:
        raise KeyError("Could not identify substation ID column.")
    return col


def _sample_raster_values(raster_path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    with rasterio.open(raster_path) as src:
        points = list(zip(lon, lat, strict=False))
        values = np.array([v[0] for v in src.sample(points)], dtype=float)
    return values


def _nanmean_safe(values: np.ndarray) -> float:
    """Return nanmean without emitting warnings for all-NaN arrays."""
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.all(np.isnan(arr)):
        return float("nan")
    return float(np.nanmean(arr))


def run(args: argparse.Namespace) -> None:
    cfg = TableConfig()
    root = Path(args.repo_root).resolve()
    out_dir = root / "outputs" / "month2_heat"
    out_dir.mkdir(parents=True, exist_ok=True)

    dc_path = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    sub_path = root / "data" / "HIFLD" / "oregon_hifld_substations.gpkg"
    line_path = root / "data" / "HIFLD" / "oregon_hifld_transmission_lines.gpkg"
    gen_path = root / "data" / "EIA-generators" / "june_generator2026.xlsx"
    water_path = root / "data" / "public_water_sources" / "WSA_v1" / "WSA_v1.shp"
    wbd_path = root / "data" / "watershed_boundaries" / "WBD_National_GDB.gdb"
    fiber_path = root / "data" / "fiber" / "bdc_41_FibertothePremises_fixed_broadband_D25_08jul2026.csv"

    rasters = {
        "cdd_hist": out_dir / "rasters" / "cdd_historical_1950_1979.tif",
        "tmax_days_hist": out_dir / "rasters" / "tmax_threshold_days_historical_1950_1979.tif",
        "heatwave_hist": out_dir / "rasters" / "heatwave_events_historical_1950_1979.tif",
        "cdd_fut": out_dir / "rasters" / "cdd_ssp585_2070_2099.tif",
        "tmax_days_fut": out_dir / "rasters" / "tmax_threshold_days_ssp585_2070_2099.tif",
        "heatwave_fut": out_dir / "rasters" / "heatwave_events_ssp585_2070_2099.tif",
    }

    dc = load_oregon_data_centers(dc_path)
    subst, lines = _load_power_layers(sub_path, line_path, cfg.transmission_min_kv)
    gens = _load_generators(gen_path)

    graph = _build_transmission_graph(lines)
    sub_id_col = _substation_id_column(subst)
    subst = subst.copy()
    subst["substation_id"] = subst[sub_id_col].astype(str)

    dc_ea = dc.to_crs("EPSG:5070")
    subst_ea = subst.to_crs("EPSG:5070")
    gens_ea = gens.to_crs("EPSG:5070")

    # Nearest and in-radius substation candidates.
    nearest = gpd.sjoin_nearest(
        dc_ea,
        subst_ea[["substation_id", "geometry"]],
        how="left",
        distance_col="nearest_substation_m",
    )
    nearest = nearest[~nearest.index.duplicated(keep="first")]

    sindex_sub = subst_ea.sindex
    sindex_gen = gens_ea.sindex
    sub_geom_lookup = subst_ea.set_index("substation_id").geometry

    water = gpd.read_file(water_path)
    if water.crs is None:
        water = water.set_crs("EPSG:4326")
    water = water.to_crs("EPSG:4326")
    water_id_col = _col_pick(water, ["WSA_AGIDF", "PWSID", "ID", "SYSTEM", "NAME"])
    water_name_col = _col_pick(water, ["SYSTEM", "SYSTEM_NAME", "NAME", "PWS_NAME"])

    wbd = gpd.read_file(wbd_path, layer=cfg.watershed_level)
    if wbd.crs is None:
        wbd = wbd.set_crs("EPSG:4326")
    wbd = wbd.to_crs("EPSG:4326")
    if "states" in wbd.columns:
        wbd = wbd[wbd["states"].fillna("").str.contains("OR")].copy()
    huc_col = _col_pick(wbd, ["HUC8", "huc8", "HUC_8", "huc"])

    dc_water = gpd.sjoin(
        dc[["id", "name", "geometry"]],
        water[[c for c in [water_id_col, water_name_col] if c] + ["geometry"]],
        how="left",
        predicate="within",
    )
    dc_water = dc_water[~dc_water.index.duplicated(keep="first")]

    dc_huc = gpd.sjoin(
        dc[["id", "name", "geometry"]],
        wbd[[c for c in [huc_col] if c] + ["geometry"]],
        how="left",
        predicate="within",
    )
    dc_huc = dc_huc[~dc_huc.index.duplicated(keep="first")]

    # Fiber provider counts by H3 cell at facility location.
    fiber = pd.read_csv(fiber_path)
    fiber = fiber[(fiber["state_usps"] == "OR") & (fiber["max_advertised_download_speed"] >= 1000)].copy()
    fiber_counts = fiber.groupby("h3_res8_id")["provider_id"].nunique().to_dict()

    edge_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for idx, dc_row in dc.iterrows():
        dc_id = str(dc_row["id"])
        dc_name = dc_row["name"]
        lon = float(dc_row["lon"])
        lat = float(dc_row["lat"])

        nearest_sub = nearest.loc[idx]
        nearest_sub_id = str(nearest_sub["substation_id"]) if pd.notna(nearest_sub["substation_id"]) else None
        nearest_dist_m = float(nearest_sub["nearest_substation_m"]) if pd.notna(nearest_sub["nearest_substation_m"]) else np.nan

        center = dc_ea.loc[idx, "geometry"]
        cand_idx = list(sindex_sub.query(center.buffer(cfg.substation_service_radius_m), predicate="intersects"))
        if cand_idx:
            cand_sub = subst_ea.iloc[cand_idx].copy()
            cand_sub["distance_m"] = cand_sub.geometry.distance(center)
            cand_sub = cand_sub[cand_sub["distance_m"] <= cfg.substation_service_radius_m].copy()
        else:
            cand_sub = subst_ea.iloc[[]].copy()

        candidate_sub_ids_uncapped = set(cand_sub["substation_id"].astype(str))
        cand_sub_for_generation = cand_sub.copy()
        if not cand_sub_for_generation.empty and cfg.max_candidate_substations > 0:
            cand_sub_for_generation = (
                cand_sub_for_generation.sort_values("distance_m")
                .head(cfg.max_candidate_substations)
                .copy()
            )

        candidate_sub_ids_for_generation = set(cand_sub_for_generation["substation_id"].astype(str))
        if nearest_sub_id:
            candidate_sub_ids_uncapped.add(nearest_sub_id)
            candidate_sub_ids_for_generation.add(nearest_sub_id)

        connected_sub_ids_uncapped = set(candidate_sub_ids_uncapped)
        connected_sub_ids_for_generation = set(candidate_sub_ids_for_generation)
        if cfg.include_transmission_one_hop:
            for sid in list(candidate_sub_ids_uncapped):
                if sid in graph:
                    connected_sub_ids_uncapped.update(str(n) for n in graph.neighbors(sid))
            for sid in list(candidate_sub_ids_for_generation):
                if sid in graph:
                    connected_sub_ids_for_generation.update(str(n) for n in graph.neighbors(sid))

        connected_geoms = [
            sub_geom_lookup[s]
            for s in connected_sub_ids_for_generation
            if s in sub_geom_lookup.index
        ]
        reachable_gen_idx: set[int] = set()
        for geom in connected_geoms:
            nearby = sindex_gen.query(geom.buffer(cfg.generation_proxy_radius_m), predicate="intersects")
            for gi in nearby:
                if gens_ea.geometry.iloc[gi].distance(geom) <= cfg.generation_proxy_radius_m:
                    reachable_gen_idx.add(int(gi))

        if reachable_gen_idx:
            gen_sel = gens.iloc[sorted(reachable_gen_idx)].copy()
            heat_vals = {}
            for key, raster in rasters.items():
                if raster.exists():
                    heat_vals[key] = _sample_raster_values(
                        raster,
                        lon=gen_sel["gen_lon"].to_numpy(),
                        lat=gen_sel["gen_lat"].to_numpy(),
                    )
                else:
                    heat_vals[key] = np.full(len(gen_sel), np.nan)

            for i, gen_row in gen_sel.iterrows():
                i_pos = gen_sel.index.get_loc(i)
                edge_rows.append(
                    {
                        "facility_id": dc_id,
                        "facility_name": dc_name,
                        "facility_lon": lon,
                        "facility_lat": lat,
                        "nearest_substation_id": nearest_sub_id,
                        "nearest_substation_distance_m": nearest_dist_m,
                        "generator_id": int(gen_row["generator_id"]),
                        "generator_name": gen_row["generator_name"],
                        "generator_unit_count": int(gen_row.get("generator_unit_count", 1)),
                        "fuel_type": gen_row["fuel_type"],
                        "prime_mover": gen_row["prime_mover"],
                        "capacity_mw": gen_row["capacity_mw"],
                        "generator_lon": float(gen_row["gen_lon"]),
                        "generator_lat": float(gen_row["gen_lat"]),
                        "cdd_historical_1950_1979": float(heat_vals["cdd_hist"][i_pos]),
                        "tmax_threshold_days_historical_1950_1979": float(heat_vals["tmax_days_hist"][i_pos]),
                        "heatwave_events_historical_1950_1979": float(heat_vals["heatwave_hist"][i_pos]),
                        "cdd_ssp585_2070_2099": float(heat_vals["cdd_fut"][i_pos]),
                        "tmax_threshold_days_ssp585_2070_2099": float(heat_vals["tmax_days_fut"][i_pos]),
                        "heatwave_events_ssp585_2070_2099": float(heat_vals["heatwave_fut"][i_pos]),
                    }
                )

            # Aggregated facility-level heat stats from reachable generators.
            cdd_hist_mean = _nanmean_safe(heat_vals["cdd_hist"])
            tmax_days_hist_mean = _nanmean_safe(heat_vals["tmax_days_hist"])
            hw_hist_mean = _nanmean_safe(heat_vals["heatwave_hist"])
            cdd_fut_mean = _nanmean_safe(heat_vals["cdd_fut"])
            tmax_days_fut_mean = _nanmean_safe(heat_vals["tmax_days_fut"])
            hw_fut_mean = _nanmean_safe(heat_vals["heatwave_fut"])

            generation_source_count = int(len(gen_sel))
            generation_unit_count = int(gen_sel.get("generator_unit_count", pd.Series(np.ones(len(gen_sel)))).sum())
            generation_by_fuel = int(pd.Series(gen_sel["fuel_type"]).nunique())
        else:
            cdd_hist_mean = np.nan
            tmax_days_hist_mean = np.nan
            hw_hist_mean = np.nan
            cdd_fut_mean = np.nan
            tmax_days_fut_mean = np.nan
            hw_fut_mean = np.nan
            generation_source_count = 0
            generation_unit_count = 0
            generation_by_fuel = 0

        h3_cell = h3.latlng_to_cell(lat, lon, 8)
        fiber_provider_count = int(fiber_counts.get(h3_cell, 0))

        wrow = dc_water.loc[idx] if idx in dc_water.index else None
        hrow = dc_huc.loc[idx] if idx in dc_huc.index else None
        water_system_id = wrow[water_id_col] if (wrow is not None and water_id_col) else np.nan
        water_system_name = wrow[water_name_col] if (wrow is not None and water_name_col) else np.nan
        watershed_huc = hrow[huc_col] if (hrow is not None and huc_col) else np.nan

        summary_rows.append(
            {
                "facility_id": dc_id,
                "facility_name": dc_name,
                "facility_lon": lon,
                "facility_lat": lat,
                "nearest_substation_id": nearest_sub_id,
                "nearest_substation_distance_m": nearest_dist_m,
                "substations_within_service_radius_count": int(len(candidate_sub_ids_uncapped)),
                "substations_used_for_generation_count": int(len(candidate_sub_ids_for_generation)),
                "connected_substations_count": int(len(connected_sub_ids_uncapped)),
                "connected_substations_used_for_generation_count": int(len(connected_sub_ids_for_generation)),
                "generation_source_count": generation_source_count,
                "generation_plant_count": generation_source_count,
                "generation_unit_count_proxy": generation_unit_count,
                "generation_fuel_type_count": generation_by_fuel,
                "water_system_id": water_system_id,
                "water_system_name": water_system_name,
                "water_source_type": "unknown",
                "watershed_huc8": watershed_huc,
                "water_source_count": 1 if pd.notna(water_system_id) else 0,
                "fiber_provider_count_gte_1gbps": fiber_provider_count,
                "cdd_historical_1950_1979_gen_mean": cdd_hist_mean,
                "tmax_threshold_days_historical_1950_1979_gen_mean": tmax_days_hist_mean,
                "heatwave_events_historical_1950_1979_gen_mean": hw_hist_mean,
                "cdd_ssp585_2070_2099_gen_mean": cdd_fut_mean,
                "tmax_threshold_days_ssp585_2070_2099_gen_mean": tmax_days_fut_mean,
                "heatwave_events_ssp585_2070_2099_gen_mean": hw_fut_mean,
                # Deferred hazards placeholders for schema stability.
                "flood_exposure_score": np.nan,
                "wildfire_exposure_score": np.nan,
                "flood_data_status": "deferred",
                "wildfire_data_status": "deferred",
                "assumption_power_proxy": (
                    "Generation linked at plant level to nearest candidate substations "
                    f"(max {cfg.max_candidate_substations}, one-hop={cfg.include_transmission_one_hop}) "
                    f"within {cfg.generation_proxy_radius_miles} miles. "
                    "Reported substation redundancy columns are uncapped."
                ),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    edges_df = pd.DataFrame(edge_rows)

    summary_df.to_csv(out_dir / "dependency_table_facility_summary.csv", index=False)
    edges_df.to_csv(out_dir / "dependency_table_generation_edges.csv", index=False)
    print(f"Wrote {len(summary_df)} facility rows and {len(edges_df)} generator-edge rows to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Oregon data-center dependency table.")
    parser.add_argument("--repo-root", default=".", help="Repository root path")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
