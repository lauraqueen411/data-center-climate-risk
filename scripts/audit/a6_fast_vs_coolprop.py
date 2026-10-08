"""A6: fast backend vs. CoolProp reference on real Oregon weather.

Compares every ``*_coolprop.csv`` in the facility output folder with the fast
output for the same model, case and period. The CoolProp files predate the
2026-10-08 parameter correction, so each is compared with both the old fast
file (superseded folder) and the current one; the closer match shows which
parameter set the CoolProp run used, and that pair is the like-for-like
backend comparison.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from audit_common import FAC_DIR, OLD_DIR, write

KEYS = ["facility_id", "date"]
COLS = ["pue_daily", "wue_daily", "pue_hot_pairing", "pue_cool_pairing"]


def compare(ref: pd.DataFrame, fast: pd.DataFrame) -> dict:
    m = ref.merge(fast, on=KEYS, suffixes=("_cp", "_fast"))
    out = {"n_matched_rows": len(m)}
    for c in COLS:
        if f"{c}_cp" not in m:
            continue
        d = m[f"{c}_fast"] - m[f"{c}_cp"]
        rel = d / m[f"{c}_cp"].where(m[f"{c}_cp"].abs() > 1e-9)
        out.update({f"{c}_mean_diff": d.mean(), f"{c}_p95_absdiff": d.abs().quantile(0.95),
                    f"{c}_max_absdiff": d.abs().max(), f"{c}_max_reldiff": rel.abs().max()})
    ann = m.assign(year=m["date"].astype(str).str.slice(0, 4)).groupby(["facility_id", "year"])[
        ["pue_daily_cp", "pue_daily_fast", "wue_daily_cp", "wue_daily_fast"]].mean().groupby(level="facility_id").mean()
    out["facility_period_pue_max_absdiff"] = (ann["pue_daily_fast"] - ann["pue_daily_cp"]).abs().max()
    out["facility_period_wue_max_absdiff"] = (ann["wue_daily_fast"] - ann["wue_daily_cp"]).abs().max()
    return out


def main() -> None:
    rows = []
    cp_files = sorted(FAC_DIR.glob("facility_pue_wue_*_coolprop.csv"))
    if not cp_files:
        print("  no *_coolprop.csv files found")
    for f in cp_files:
        meta = re.match(r"facility_pue_wue_(.+)_ssp370_(ae-chiller|we-chiller|chiller-only)_(\w+?)_coolprop\.csv", f.name)
        gcm, arch, period = meta.groups()
        ref = pd.read_csv(f)
        cols = KEYS + [c for c in COLS if c in ref]
        ref = ref[cols]
        for label, folder in (("old_fast", OLD_DIR), ("current_fast", FAC_DIR)):
            fp = folder / f.name.replace("_coolprop.csv", "_fast.csv")
            if not fp.exists():
                continue
            r = {"coolprop_file": f.name, "gcm": gcm, "archetype": arch, "period": period, "compared_with": label,
                 "coolprop_n_rows": len(ref), "coolprop_n_facilities": ref["facility_id"].nunique()}
            r.update(compare(ref, pd.read_csv(fp, usecols=cols)))
            rows.append(r)
    out = pd.DataFrame(rows)
    write(out, "a6_fast_vs_coolprop.csv")
    if len(out):
        show = ["archetype", "period", "compared_with", "n_matched_rows", "pue_daily_max_reldiff",
                "facility_period_pue_max_absdiff", "facility_period_wue_max_absdiff"]
        print(out[show].to_string(index=False))


if __name__ == "__main__":
    main()
