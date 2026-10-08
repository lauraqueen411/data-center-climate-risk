"""A4 + A5: grids, time axes and source files, per model and period.

A4: the runner pairs temperature and humidity by array position, not by
date, so the two time axes must be identical. For every model and period,
compares the clipped temperature time axis (as the runner builds it) with
the cached humidity time axis, and the lat/lon of the cells the facilities
are sampled from in each. Spot-checks one cached humidity value against the
raw humidity file, which confirms the cache file belongs to that model,
scenario and period.

A5: lists the exact temperature and humidity files resolved for every model
and period, with scenario and member checks.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from audit_common import (HUMIDITY_ROOT, LOCA2_ROOT, PERIODS, ROOT, facilities, gcms, lon_query, member, scenario,
                          write)
from climate_risk_dc.climate.loca2_io import humidity_path_config, resolve_loca2_file, temperature_path_config
from run_pue_wue_gridded import GriddedRunConfig

CACHE = ROOT / "outputs/ensemble_pue_wue/humidity_cache"


def nearest_idx(values: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.abs(values[:, None] - q[None, :]).argmin(axis=0)


def main() -> None:
    cfg = GriddedRunConfig()
    fac = facilities()
    lon, lat = fac["lon"].to_numpy(), fac["lat"].to_numpy()
    files, grids = [], []
    for gcm in gcms():
        mem = member(gcm)
        for period, (start, end) in PERIODS.items():
            scen = scenario(period)
            tcfg, hcfg = temperature_path_config(LOCA2_ROOT, gcm, mem), humidity_path_config(HUMIDITY_ROOT, gcm, mem)
            p = {v: resolve_loca2_file(tcfg, scen, v, start, end) for v in ("tasmax", "tasmin")}
            p.update({v: resolve_loca2_file(hcfg, scen, v, start, end) for v in ("hursmax", "hursmin")})
            for v, path in p.items():
                files.append({"gcm": gcm, "period": period, "variable": v, "file": str(path),
                              "scenario_in_path": f"/{scen}/" in str(path), "member_in_path": f"/{mem}/" in str(path),
                              "manifest_member": mem})

            # Coordinates only (no data load): the runner's clip_to_bbox + select_year_range, reproduced on coords.
            traw = xr.open_dataset(p["tasmax"])
            t_lon = traw["lon"].values
            t_lon = np.where(t_lon > 180.0, t_lon - 360.0, t_lon)
            t_lat = traw["lat"].values
            t_years = traw["time"].dt.year.values
            t_sel = traw["time"].values[(t_years >= start) & (t_years <= end)]
            in_lon = (t_lon >= cfg.lon_min) & (t_lon <= cfg.lon_max)
            in_lat = (t_lat >= cfg.lat_min) & (t_lat <= cfg.lat_max)
            t_lon, t_lat = np.sort(t_lon[in_lon]), t_lat[in_lat]
            cache_path = CACHE / f"hursmax_hursmin_{gcm}_{scen}_{period}.nc"
            row = {"gcm": gcm, "period": period, "cache_exists": cache_path.exists()}
            if cache_path.exists():
                hc = xr.open_dataset(cache_path)
                t_time = pd.Index([str(t)[:10] for t in t_sel])
                h_time = pd.Index([str(t)[:10] for t in hc["time"].values])
                row.update({
                    "n_time_temperature": len(t_time), "n_time_humidity": len(h_time),
                    "time_axes_identical": t_time.equals(h_time),
                    "first_time_mismatch": next((f"{a} vs {b}" for a, b in zip(t_time, h_time) if a != b), ""),
                    "temperature_calendar": str(traw["time"].encoding.get("calendar", "")),
                    "humidity_calendar": str(hc["time"].encoding.get("calendar", "")),
                })
                ti_lat, ti_lon = nearest_idx(t_lat, lat), nearest_idx(t_lon, lon)
                hi_lat, hi_lon = nearest_idx(hc["lat"].values, lat), nearest_idx(hc["lon"].values, lon)
                row["max_cell_lat_diff_deg"] = float(np.abs(t_lat[ti_lat] - hc["lat"].values[hi_lat]).max())
                row["max_cell_lon_diff_deg"] = float(np.abs(t_lon[ti_lon] - hc["lon"].values[hi_lon]).max())
                # Spot check: cache vs raw humidity file, first facility, middle of the period.
                k = len(h_time) // 2
                cached_val = float(hc["hursmin"].isel(time=k, lat=hi_lat[0], lon=hi_lon[0]).values)
                raw = xr.open_dataset(p["hursmin"])
                raw_da = raw[list(raw.data_vars)[0]]
                q = lon_query(lon[0], raw_da["lon"].values)
                raw_val = float(raw_da.sel(time=h_time[k]).sel(lat=float(hc["lat"].values[hi_lat[0]]), lon=q,
                                                               method="nearest").squeeze().values)
                row.update({"spot_date": h_time[k], "spot_cache_hursmin": cached_val, "spot_raw_hursmin": raw_val,
                            "spot_abs_diff": abs(cached_val - raw_val)})
            grids.append(row)
        print(f"  {gcm} done")
    files = pd.DataFrame(files)
    grids = pd.DataFrame(grids)
    write(files, "a5_source_files.csv")
    write(grids, "a4_grids_time_axes.csv")
    print(f"  A5: {int((~files['scenario_in_path']).sum())} files without the expected scenario, "
          f"{int((~files['member_in_path']).sum())} without the manifest member")
    g = grids[grids["cache_exists"]]
    print(f"  A4: {int((~g['time_axes_identical']).sum())} of {len(g)} model-periods with differing time axes; "
          f"max cell offset lat {g['max_cell_lat_diff_deg'].max():.2e}, lon {g['max_cell_lon_diff_deg'].max():.2e} deg; "
          f"max cache-vs-raw spot diff {g['spot_abs_diff'].max():.2e}; missing caches {int((~grids['cache_exists']).sum())}")


if __name__ == "__main__":
    main()
