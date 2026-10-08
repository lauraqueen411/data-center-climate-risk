"""Count the Case 1 weather pairings whose economizer decision depends on the
upper supply-air RH setpoint (60% as run vs. the Table B.1 midpoint 77.5%).

In ``cooling_archetypes_fast._air_side_economizer_fast`` the setpoint only
changes the outcome in branch A (``t_lw <= T_oa < t_up``, 14-31 C): the
outdoor humidity ratio exceeds ``min(W(T_oa, RH_up), W(Tdp = dp_up))``.
Below ``t_lw`` the RH limit carries no energy penalty, and at or above
``t_up`` the 21 C dewpoint limit is always the tighter one. Comparing the two
setpoints at the same temperature and pressure, a pairing is affected when

    t_lw <= T < t_up  and  rh_up_run < RH <= rh_up_mid  and  Tdp <= dp_up

(dehumidification with the chiller at 60%, free cooling at 77.5%). The test
is pressure-independent, so no elevation is needed. Each daily PUE averages a
hot pairing (tasmax, hursmin) and a cool pairing (tasmin, hursmax); both are
counted.

Reads the Case 1 (``ae-chiller``) historical and mid-century daily files in
``outputs/ensemble_pue_wue_facilities/`` (no model is run) and writes
``outputs/ensemble_pue_wue_facilities/analysis/rh_up_sensitivity_counts.csv``:
one row per facility, model, period, month and pairing with the number of
days and the number of affected days.

Assumptions log
----------------
- Dewpoint uses the Magnus formula (Alduchov and Eskridge 1996 constants);
  the error near 21 C is well under 0.1 C, negligible for a count.
- Thresholds come from ``Case1Config`` (t_lw 14 C, t_up 31 C, dp_up 21 C);
  ``RH_UP_MID`` is the Table B.1 midpoint of 60-95%.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.cooling_archetypes import Case1Config

RH_UP_MID = 77.5
PERIODS = ("historical", "midcentury")
FNAME_RE = re.compile(r"facility_pue_wue_(?P<gcm>.+)_ssp370_ae-chiller_(?P<period>historical|midcentury)_fast\.csv")


def dewpoint_c(t_c: np.ndarray, rh_pct: np.ndarray) -> np.ndarray:
    a, b = 17.625, 243.04
    g = np.log(np.clip(rh_pct, 1e-6, None) / 100.0) + a * t_c / (b + t_c)
    return b * g / (a - g)


def affected(t_c: np.ndarray, rh_pct: np.ndarray, cfg: Case1Config) -> np.ndarray:
    return ((t_c >= cfg.t_lw_c) & (t_c < cfg.t_up_c) & (rh_pct > cfg.rh_up_pct) & (rh_pct <= RH_UP_MID)
            & (dewpoint_c(t_c, rh_pct) <= cfg.dp_up_c))


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    fac_dir = root / "outputs/ensemble_pue_wue_facilities"
    cfg = Case1Config()
    files = sorted(f for f in fac_dir.glob("facility_pue_wue_*_ssp370_ae-chiller_*_fast.csv") if FNAME_RE.match(f.name))
    if len(files) != 40:
        raise FileNotFoundError(f"Expected 40 Case 1 historical/mid-century files, found {len(files)}")
    parts = []
    for f in files:
        meta = FNAME_RE.match(f.name).groupdict()
        d = pd.read_csv(f, usecols=["facility_id", "date", "tasmax_c", "tasmin_c", "hursmax_pct", "hursmin_pct"])
        d["month"] = d["date"].str.slice(5, 7).astype(int)
        d["hot"] = affected(d["tasmax_c"].to_numpy(), d["hursmin_pct"].to_numpy(), cfg)
        d["cool"] = affected(d["tasmin_c"].to_numpy(), d["hursmax_pct"].to_numpy(), cfg)
        g = d.groupby(["facility_id", "month"]).agg(n_days=("hot", "size"), hot=("hot", "sum"), cool=("cool", "sum"))
        g = g.reset_index().melt(id_vars=["facility_id", "month", "n_days"], var_name="pairing", value_name="n_affected")
        g["gcm"], g["period"] = meta["gcm"], meta["period"]
        parts.append(g)
        print(f"{f.name}: hot {d['hot'].mean():.3%}, cool {d['cool'].mean():.3%} of days affected")
    out = pd.concat(parts, ignore_index=True)
    path = fac_dir / "analysis" / "rh_up_sensitivity_counts.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    print(f"\nWrote {path} ({len(out)} rows)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Count Case 1 pairings sensitive to the upper RH setpoint.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
