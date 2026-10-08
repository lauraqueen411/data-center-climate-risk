"""A2: raw LOCA2 / humidity values vs. the daily CSVs, end to end.

For one Columbia Plateau and one Willamette Valley facility, two models and
a handful of dates, read tasmax, tasmin, hursmax and hursmin straight from
the source NetCDF files at the nearest cell (an independent lookup, using
xarray ``sel(method="nearest")`` rather than the pipeline's sampler) and
compare with the daily CSV. Catches unit conversion, cell choice, longitude
convention, date alignment and a hursmax/hursmin swap.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from audit_common import ROOT, HUMIDITY_ROOT, LOCA2_ROOT, daily_path, facilities, lon_query, member, scenario, write
from climate_risk_dc.climate.loca2_io import humidity_path_config, resolve_loca2_file, temperature_path_config

GCMS = ["ACCESS-CM2", "MPI-ESM1-2-HR"]
DATES = {"historical": ["1990-01-15", "1990-07-15", "1988-02-29", "2014-12-31"], "midcentury": ["2060-07-15"]}
WINDOWS = {"historical": (1985, 2014), "midcentury": (2045, 2074)}


def raw_value(path, lon, lat, date):
    ds = xr.open_dataset(path)
    da = ds[list(ds.data_vars)[0]]
    q = lon_query(lon, da["lon"].values)
    try:
        pt = da.sel(time=date).sel(lat=lat, lon=q, method="nearest")
    except KeyError:
        return np.nan, np.nan, np.nan, da.attrs.get("units", ""), ""
    pt = pt.squeeze()
    return (float(pt.values), float(pt["lat"]), float(pt["lon"]), da.attrs.get("units", ""), list(ds.data_vars)[0])


def main() -> None:
    fac = facilities()
    picks = []
    for eco in ("Columbia Plateau", "Willamette Valley"):
        sub = fac[fac["ecoregion"] == eco].sort_values("facility_id")
        picks.append(sub.iloc[0])
    rows = []
    for gcm in GCMS:
        mem = member(gcm)
        for period, dates in DATES.items():
            start, end = WINDOWS[period]
            scen = scenario(period)
            paths = {
                "tasmax_c": resolve_loca2_file(temperature_path_config(LOCA2_ROOT, gcm, mem), scen, "tasmax", start, end),
                "tasmin_c": resolve_loca2_file(temperature_path_config(LOCA2_ROOT, gcm, mem), scen, "tasmin", start, end),
                "hursmax_pct": resolve_loca2_file(humidity_path_config(HUMIDITY_ROOT, gcm, mem), scen, "hursmax", start, end),
                "hursmin_pct": resolve_loca2_file(humidity_path_config(HUMIDITY_ROOT, gcm, mem), scen, "hursmin", start, end),
            }
            csv = pd.read_csv(daily_path(gcm, "ae-chiller", period),
                              usecols=["facility_id", "date", "tasmax_c", "tasmin_c", "hursmax_pct", "hursmin_pct"])
            csv["day"] = csv["date"].astype(str).str.slice(0, 10)
            for f in picks:
                for date in dates:
                    c = csv[(csv["facility_id"] == f["facility_id"]) & (csv["day"] == date)]
                    for col, path in paths.items():
                        val, rlat, rlon, units, var = raw_value(path, f["lon"], f["lat"], date)
                        if col.endswith("_c") and np.isfinite(val):
                            val = val - 273.15 if val > 150 else val
                        csv_val = float(c[col].iloc[0]) if len(c) else np.nan
                        rows.append({
                            "gcm": gcm, "period": period, "ecoregion": f["ecoregion"], "facility_id": f["facility_id"],
                            "date": date, "variable": col, "raw_value": val, "csv_value": csv_val,
                            "abs_diff": abs(val - csv_val) if np.isfinite(val) and np.isfinite(csv_val) else np.nan,
                            "raw_units": units, "raw_var_name": var, "cell_lat": rlat, "cell_lon": rlon,
                            "pipeline_grid_lat": f.get("grid_lat", np.nan), "pipeline_grid_lon": f.get("grid_lon", np.nan),
                            "date_in_csv": len(c) > 0, "file": path.name,
                        })
    out = pd.DataFrame(rows)
    write(out, "a2_raw_inputs.csv")
    ok = out.dropna(subset=["abs_diff"])
    print(f"  compared {len(ok)} values; max |diff| {ok['abs_diff'].max():.2e}; "
          f"{(ok['abs_diff'] > 1e-3).sum()} differ by > 0.001; "
          f"{(~out['date_in_csv']).sum()} requested dates absent from the CSV (e.g. 29 Feb in no-leap calendars)")


if __name__ == "__main__":
    main()
