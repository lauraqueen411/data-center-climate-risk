"""A7: chiller COP extrapolation, plus a check that saved outputs reproduce.

1. The training range of ``COP_2.pkl``, read from the GP itself
   (``X_train_``), and the GP's COP curve over -15 to 40 C wet-bulb for each
   case's chiller load (``a7_cop_curve.csv``).
2. For every model and period (Case 1 weather files), the share of hot and
   cool pairings whose ``Twb + AT_CT`` lies outside the training range, by
   ecoregion and month (``a7_out_of_range_share.csv``).
3. Sensitivity: for a sample of models, periods and cases, recompute PUE/WUE
   from the saved weather with (a) the GP as used and (b) the GP inputs
   clamped to the training range. (a) is compared with the saved
   ``pue_daily``/``wue_daily`` (does the current code reproduce the files?);
   (b) minus (a) is the effect of extrapolation on the period mean and on the
   mid-century change, per facility (``a7_clamp_effect.csv``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from audit_common import PERIODS, ROOT, daily_path, facilities, gcms, pressure_by_facility, write
from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case2Config, Case5Config, load_cop_water_cooled_gp
from climate_risk_dc.climate.cooling_archetypes_fast import (_wetbulb_from_rh, pue_wue_case1, pue_wue_case2,
                                                             pue_wue_case5)

CASES = {"ae-chiller": (pue_wue_case1, Case1Config()), "we-chiller": (pue_wue_case2, Case2Config()),
         "chiller-only": (pue_wue_case5, Case5Config())}
SAMPLE = [("ACCESS-CM2", "ae-chiller"), ("MPI-ESM1-2-HR", "ae-chiller"), ("ACCESS-CM2", "chiller-only"),
          ("ACCESS-CM2", "we-chiller")]
SAMPLE_PERIODS = ("historical", "midcentury")


class ClampedGP:
    def __init__(self, gp, lo, hi):
        self.gp, self.lo, self.hi = gp, lo, hi

    def predict(self, x):
        return self.gp.predict(np.clip(x, self.lo, self.hi))


def main() -> None:
    gp = load_cop_water_cooled_gp(ROOT / "data/lei_masanet_2022/COP_2.pkl")
    xt = np.asarray(gp.X_train_)
    lo, hi = xt.min(axis=0), xt.max(axis=0)
    print(f"  GP training range: Twb+AT {lo[0]:.2f}-{hi[0]:.2f} C, load {lo[1]:.2f}-{hi[1]:.2f} ({len(xt)} points)")
    at = Case1Config().at_ct_c
    twb = np.linspace(-15, 40, 221)
    curve = []
    for name, (_, cfg) in CASES.items():
        x = np.column_stack([twb + cfg.at_ct_c, np.full_like(twb, cfg.chiller_load)])
        curve.append(pd.DataFrame({"archetype": name, "twb_c": twb, "twb_plus_at_c": twb + cfg.at_ct_c,
                                   "chiller_load": cfg.chiller_load, "cop_gp": gp.predict(x),
                                   "cop_clamped": gp.predict(np.clip(x, lo, hi)),
                                   "in_training_range": (x[:, 0] >= lo[0]) & (x[:, 0] <= hi[0])}))
    write(pd.concat(curve), "a7_cop_curve.csv")
    write(pd.DataFrame({"feature": ["twb_plus_at_c", "chiller_load"], "train_min": lo, "train_max": hi}), "a7_training_range.csv")

    fac = facilities()
    pressure = pressure_by_facility(fac)
    eco = dict(zip(fac["facility_id"], fac.get("ecoregion", pd.Series(dtype=str))))

    share = []
    for gcm in gcms():
        for period in PERIODS:
            d = pd.read_csv(daily_path(gcm, "ae-chiller", period),
                            usecols=["facility_id", "date", "tasmax_c", "tasmin_c", "hursmax_pct", "hursmin_pct"])
            p = d["facility_id"].map(pressure).to_numpy()
            month = d["date"].astype(str).str.slice(5, 7).astype(int)
            for pairing, t, rh in (("hot", "tasmax_c", "hursmin_pct"), ("cool", "tasmin_c", "hursmax_pct")):
                x = _wetbulb_from_rh(d[t].to_numpy(), d[rh].to_numpy() / 100.0, p) + at
                g = pd.DataFrame({"ecoregion": d["facility_id"].map(eco), "month": month,
                                  "below": x < lo[0], "above": x > hi[0]})
                s = g.groupby(["ecoregion", "month"]).agg(n=("below", "size"), n_below=("below", "sum"),
                                                          n_above=("above", "sum")).reset_index()
                share.append(s.assign(gcm=gcm, period=period, pairing=pairing))
        print(f"  share: {gcm} done")
    share = pd.concat(share, ignore_index=True)
    write(share, "a7_out_of_range_share.csv")
    tot = share.groupby(["period", "pairing"])[["n", "n_below", "n_above"]].sum()
    tot["pct_below"], tot["pct_above"] = 100 * tot["n_below"] / tot["n"], 100 * tot["n_above"] / tot["n"]
    print(tot.round(2).to_string())

    effect = []
    clamped = ClampedGP(gp, lo, hi)
    for gcm, arch in SAMPLE:
        fn, cfg = CASES[arch]
        per = {}
        for period in SAMPLE_PERIODS:
            d = pd.read_csv(daily_path(gcm, arch, period))
            p = d["facility_id"].map(pressure).to_numpy()
            res = {}
            for label, model in (("gp", gp), ("clamped", clamped)):
                ph, wh = fn(d["tasmax_c"].to_numpy(), d["hursmin_pct"].to_numpy(), p, model, cfg)
                pc, wc = fn(d["tasmin_c"].to_numpy(), d["hursmax_pct"].to_numpy(), p, model, cfg)
                res[f"pue_{label}"], res[f"wue_{label}"] = (ph + pc) / 2, (wh + wc) / 2
            frame = pd.DataFrame({"facility_id": d["facility_id"], "year": d["date"].astype(str).str.slice(0, 4),
                                  "pue_saved": d["pue_daily"], "wue_saved": d["wue_daily"], **res})
            repro_pue = float(np.nanmax(np.abs(frame["pue_gp"] - frame["pue_saved"])))
            repro_wue = float(np.nanmax(np.abs(frame["wue_gp"] - frame["wue_saved"])))
            print(f"  {gcm} {arch} {period}: reproduce saved files, max |diff| PUE {repro_pue:.2e}, WUE {repro_wue:.2e}")
            ann = frame.groupby(["facility_id", "year"]).mean(numeric_only=True).groupby(level="facility_id").mean()
            ann["repro_max_abs_pue"], ann["repro_max_abs_wue"] = repro_pue, repro_wue
            per[period] = ann
        for period, ann in per.items():
            for fid, r in ann.iterrows():
                row = {"gcm": gcm, "archetype": arch, "period": period, "facility_id": fid, "ecoregion": eco.get(fid),
                       "pue_gp": r["pue_gp"], "pue_clamped": r["pue_clamped"], "wue_gp": r["wue_gp"],
                       "wue_clamped": r["wue_clamped"], "pue_saved": r["pue_saved"],
                       "repro_max_abs_pue": r["repro_max_abs_pue"], "repro_max_abs_wue": r["repro_max_abs_wue"]}
                if period != "historical":
                    h = per["historical"].loc[fid]
                    row["d_pue_gp"] = r["pue_gp"] - h["pue_gp"]
                    row["d_pue_clamped"] = r["pue_clamped"] - h["pue_clamped"]
                    row["d_wue_gp"] = r["wue_gp"] - h["wue_gp"]
                    row["d_wue_clamped"] = r["wue_clamped"] - h["wue_clamped"]
                effect.append(row)
    effect = pd.DataFrame(effect)
    write(effect, "a7_clamp_effect.csv")
    e = effect.dropna(subset=["d_pue_gp"])
    e = e.assign(change_effect=e["d_pue_clamped"] - e["d_pue_gp"], level_effect=e["pue_clamped"] - e["pue_gp"])
    print(e.groupby(["gcm", "archetype", "ecoregion"])[["d_pue_gp", "change_effect", "level_effect"]].mean()
          .round(5).to_string())


if __name__ == "__main__":
    main()
