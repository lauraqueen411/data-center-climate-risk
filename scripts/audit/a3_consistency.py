"""A3: physical consistency and calendars across all daily weather.

Reads the weather columns of every Case 1 daily file (the weather is the
same for all three cases; that is checked on one model per period). Per
model and period, reports: facilities, days, years, days per year (which
reveals the calendar), 29 Feb / 30 Feb / 31st-of-30-day-month presence,
tasmin > tasmax, hursmin > hursmax, RH outside 0-100%, NaNs, and whether the
years match the period.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from audit_common import ARCHETYPES, PERIODS, daily_path, gcms, write

WEATHER = ["tasmax_c", "tasmin_c", "hursmax_pct", "hursmin_pct"]


def main() -> None:
    rows, same = [], []
    for gcm in gcms():
        for period, (start, end) in PERIODS.items():
            d = pd.read_csv(daily_path(gcm, "ae-chiller", period), usecols=["facility_id", "date"] + WEATHER)
            day = d["date"].astype(str).str.slice(0, 10)
            y = day.str.slice(0, 4).astype(int)
            md = day.str.slice(5, 10)
            one = d["facility_id"] == d["facility_id"].iloc[0]
            per_year = y[one].value_counts()
            rows.append({
                "gcm": gcm, "period": period, "n_facilities": d["facility_id"].nunique(), "n_rows": len(d),
                "year_min": y.min(), "year_max": y.max(), "years_match_period": (y.min(), y.max()) == (start, end),
                "n_years": y.nunique(), "days_per_year_min": per_year.min(), "days_per_year_max": per_year.max(),
                "has_feb29": (md[one] == "02-29").any(), "has_feb30": (md[one] == "02-30").any(),
                "has_31st_of_30day_month": md[one].isin(["04-31", "06-31", "09-31", "11-31"]).any(),
                "date_example": d["date"].iloc[0],
                "n_tasmin_gt_tasmax": int((d["tasmin_c"] > d["tasmax_c"]).sum()),
                "n_hursmin_gt_hursmax": int((d["hursmin_pct"] > d["hursmax_pct"]).sum()),
                "n_rh_below_0": int(((d["hursmin_pct"] < 0) | (d["hursmax_pct"] < 0)).sum()),
                "n_rh_above_100": int(((d["hursmin_pct"] > 100) | (d["hursmax_pct"] > 100)).sum()),
                "rh_max_seen": float(d[["hursmin_pct", "hursmax_pct"]].max().max()),
                "n_nan": int(d[WEATHER].isna().sum().sum()),
                "tasmax_min": float(d["tasmax_c"].min()), "tasmax_max": float(d["tasmax_c"].max()),
                "tasmin_min": float(d["tasmin_c"].min()),
            })
            if gcm == gcms()[0]:
                for arch in ARCHETYPES[1:]:
                    o = pd.read_csv(daily_path(gcm, arch, period), usecols=["facility_id", "date"] + WEATHER)
                    same.append({"gcm": gcm, "period": period, "archetype": arch, "same_rows": len(o) == len(d),
                                 "max_abs_weather_diff": float(np.nanmax(np.abs(o[WEATHER].to_numpy() - d[WEATHER].to_numpy())))
                                 if len(o) == len(d) else np.nan})
        print(f"  {gcm} done")
    out = pd.DataFrame(rows)
    write(out, "a3_consistency.csv")
    write(pd.DataFrame(same), "a3_weather_same_across_cases.csv")
    print("  calendars (days per year min-max) by model:")
    cal = out.groupby("gcm").agg(dmin=("days_per_year_min", "min"), dmax=("days_per_year_max", "max"),
                                 feb29=("has_feb29", "any"), feb30=("has_feb30", "any"))
    print(cal.to_string())
    flags = ["n_tasmin_gt_tasmax", "n_hursmin_gt_hursmax", "n_rh_below_0", "n_rh_above_100", "n_nan"]
    print("  totals:", {k: int(out[k].sum()) for k in flags}, "| periods with wrong years:",
          int((~out["years_match_period"]).sum()))


if __name__ == "__main__":
    main()
