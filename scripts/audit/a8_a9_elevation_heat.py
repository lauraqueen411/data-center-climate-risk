"""A8 + A9: elevation/pressure per facility, and heat-index ensemble membership.

A8: elevation and pressure exactly as the runner computes them, the elevation
file's units and grid, and each facility's distance to its LOCA2 cell centre.
A9: the models in the heat-index facility outputs vs. the PUE/WUE manifest,
and the threshold/CDD settings the heat run used.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import xarray as xr

from audit_common import ELEVATION, ROOT, facilities, gcms, manifest, pressure_by_facility, write
from climate_risk_dc.climate.heat_indices import HeatIndexConfig


def main() -> None:
    fac = facilities()
    p = pressure_by_facility(fac)
    ds = xr.open_dataset(ELEVATION)
    el = ds["Elevation"]
    from run_pue_wue_facilities import sample_elevation_at_points

    fac["elevation_m"] = sample_elevation_at_points(ELEVATION, fac["lon"].to_numpy(), fac["lat"].to_numpy())
    fac["pressure_pa"] = fac["facility_id"].map(p)
    if "grid_lat" in fac:
        km_lat = 111.2
        fac["dist_to_cell_centre_km"] = np.hypot((fac["lat"] - fac["grid_lat"]) * km_lat,
                                                 (fac["lon"] - fac["grid_lon"]) * km_lat * np.cos(np.radians(fac["lat"])))
    write(fac, "a8_elevation_pressure.csv")
    lon = el["lon"].values
    print(f"  elevation units '{el.attrs.get('units', '')}', grid lon {lon.min():.4f}..{lon.max():.4f} "
          f"step {np.diff(lon).mean():.5f}, lat step {np.diff(el['lat'].values).mean():.5f}")
    print(f"  facility elevation {fac['elevation_m'].min():.0f}-{fac['elevation_m'].max():.0f} m, "
          f"pressure {fac['pressure_pa'].min():.0f}-{fac['pressure_pa'].max():.0f} Pa")
    if "dist_to_cell_centre_km" in fac:
        print(f"  distance to LOCA2 cell centre: max {fac['dist_to_cell_centre_km'].max():.2f} km "
              f"(half a 1/16 deg cell is about 3.5 km N-S)")

    heat_dir = ROOT / "outputs/month2_heat_ensemble/facilities"
    heat = sorted({re.match(r"facility_heat_indices_(.+)_ssp370_\w+\.csv", f.name).group(1)
                   for f in heat_dir.glob("facility_heat_indices_*_ssp370_*.csv")})
    pue = gcms()
    cfg = HeatIndexConfig()
    out = pd.DataFrame([{"heat_models": len(heat), "pue_models": len(pue),
                         "in_heat_not_pue": "; ".join(sorted(set(heat) - set(pue))),
                         "in_pue_not_heat": "; ".join(sorted(set(pue) - set(heat))),
                         "threshold_c": cfg.threshold_c, "cdd_base_c": cfg.cdd_base_c,
                         "members": "; ".join(f"{r.model}={r.member}" for r in manifest().itertuples())}])
    write(out, "a9_heat_ensemble.csv")
    print(out.drop(columns="members").to_string(index=False))


if __name__ == "__main__":
    main()
