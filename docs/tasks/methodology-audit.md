# Task: final methodology audit of the heat chapter (code vs. methods vs. results)

## Why

On 2026-10-08 a line-by-line check against Lei and Masanet (2022) Table B.1 found three equipment parameters that were not the midpoints the methods claim (Case 1 upper supply-air RH 60% instead of 77.5%, and two Case 5 values). Correcting them halved the Willamette Valley signal. The methods text, the code comments and an earlier "cross-checked" note all said the values were right. This audit looks for any other mismatch between what the chapter says, what the code does, and what the source papers specify, before server access ends.

## Constraints

- **Server access ends 2026-10-11.** After that only text and figure changes are possible. Anything that needs the server, whether a check or a rerun, goes first. Finish Part A by the end of 2026-10-09 so that one rerun day remains if needed.
- **Audit first, fix second.** Don't change code, outputs or text while auditing. Record every finding, then stop and report. Laura decides what gets fixed.
- **Verify, don't trust.** Docstrings, comments, Assumptions logs and earlier task reports are claims to check, not evidence. Where possible, confirm each item from the source: the paper's PDF in `papers/`, the raw data file, or a computed value. Extract the PDFs with `pypdf`, which can be installed to the scratchpad with `pip install --target`.

## Sources of truth

| What | Where |
|---|---|
| **Chapter methods text (current)** | `docs/results-figures/methodology.md`: shared methods plus the method and caption for each of Figures 2–5 |
| Captions | `docs/results-figures/*_caption.txt` (must match `methodology.md`) |

All other text drafts (`docs/heat-methodology-draft.md`, `docs/methods-document.md`, `docs/results-draft.md`, the outline and earlier task reports) are out of date. Don't read them or audit against them.
| Source papers | `papers/Lei and Masanet - 2020 ...pdf`, `papers/Lei and Masanet - 2022 ...pdf` |
| Model code | `src/climate_risk_dc/climate/cooling_archetypes.py`, `cooling_archetypes_fast.py`, `heat_indices.py`, `loca2_io.py`, `geo.py` |
| Run scripts | `scripts/run_pue_wue_facilities.py`, `run_pue_wue_gridded.py` (shared helpers), `run_heat_indices_ensemble.py`, `aggregate_heat_ensemble.py`, `build_pue_wue_ensemble_summary.py` |
| Figure and table scripts | `scripts/fig_show_everything.py`, `fig3_maps.py`, `fig4_monthly.py`, `fig5_mechanism.py`, `build_results_table.py` |

## Output

Write `docs/audit/methodology-audit.md`, plus `docs/audit/findings.csv` with one row per finding and these columns:

| Column | Content |
|---|---|
| id | A1, B3, … |
| area | data / model / aggregation / text / figure |
| finding | One sentence |
| evidence | File and line, paper page, or the computed value |
| type | **bug** (code does something wrong) / **mismatch** (text and code disagree) / **undocumented** (a choice the text doesn't state) / **ok** |
| severity | **high** (changes a reported number or conclusion) / **medium** (changes a number slightly, or a method statement is false) / **low** (wording) |
| needs_server | yes / no |
| fix_cost | e.g. "text only", "rebuild figures locally", "rerun Case 1, ~1 h" |

Report the high-severity, server-dependent findings to Laura as soon as they're found. Don't wait for the full report.

---

## Part A: checks that need the server (do first)

Run these from the server repo root, with `PYTHONPATH=$PWD/src`. Keep each check to a small script under `scripts/audit/` that writes a CSV, and sync the CSVs back.

**A1. Runner sanity warnings.** `run_pue_wue_facilities.py` prints `Sanity check WARNING` and carries on when `check_pue_wue_plausible` fails. Grep every per-model log from the rerun (`outputs/ensemble_pue_wue_facilities/logs/rerun_2026-10-08/`), and any logs from the original 2026-09-10 run, for warnings. Report which models, cases and periods triggered them, and why.

**A2. Raw inputs end to end.** For two facilities (one CP, one WV), two models and three dates (a winter day, a summer day and 29 February where it exists), read `tasmax`, `tasmin`, `hursmax` and `hursmin` directly from the LOCA2 and humidity NetCDF files at the nearest cell. Compare with the daily CSV. This checks the unit conversion (K to °C), which grid cell is picked, the longitude convention (0–360 vs ±180), date alignment, and that `hursmax` and `hursmin` are not swapped.

**A3. Physical consistency across all daily files.** Report counts of:
- days with `tasmin > tasmax`
- days with `hursmin > hursmax`
- RH outside 0–100%
- missing values
- the number of days per year for each model

That last count shows each model's calendar (standard, no-leap or 360-day). Then check that the month grouping in `facility_month_summary.csv` and the period year ranges are right under every calendar. A 360-day calendar, for example, has no 31st of any month.

**A4. Temperature and humidity on the same grid.** Confirm that the humidity cache's lat/lon arrays match the LOCA2 temperature grid exactly, so the same index means the same cell. Also confirm that the humidity cache was built for the right model, scenario and period: it is keyed by file name, so check each key against its content.

**A5. Scenario stitching.** From the run logs or by re-resolving, list the exact LOCA2 and humidity files used for each model and period. Historical must come from the `historical` experiment and future periods from `ssp370`. Confirm the member ID matches `gcm_manifest.csv`, and that it is the same member for temperature and humidity.

**A6. Fast vs CoolProp on real Oregon weather.** The CoolProp reference outputs for ACCESS-CM2, historical period (`*_historical_coolprop.csv`), were computed with the old parameters, so compare them with the old fast outputs in `superseded_rh_up60_2026-10-08/`, which used the same parameters. Report the PUE and WUE differences per facility (mean, 95th percentile, maximum), the difference in the facility-level annual means, and whether any economizer branch flips. If the differences are larger than the 0.1–0.35% claimed in the docstring, flag them.

**A7. Chiller COP extrapolation.** `COP_2.pkl` was trained on `Twb_oa + AT_CT` from 15.6 to 29.4 °C and load fractions 0.25–1.0. For every Case 1 and Case 5 day, report the share that falls outside that range, by ecoregion and period. Also plot the predicted COP against `Twb + AT_CT` across the full Oregon range, to see whether the GP drifts back toward its prior mean outside the training range. The Assumptions log says this "should be checked per-GCM" and there is no record that it was. If COP behaves oddly in winter or on extreme days, that feeds straight into PUE.

**A8. Elevation and pressure.** List the elevation and pressure at each facility. Check the elevation units (m), that the nearest-cell sampling matches the temperature grid, and the pressure range (Prineville about 91 kPa, Hillsboro about 100 kPa).

**A9. Heat-index rasters.** For the exceedance count, confirm the threshold (`> 35` vs `>= 35`), that the count is per year and then averaged over the period, and that the CDD formula uses base 18.3 °C on `(tasmax + tasmin)/2`. Confirm the same 20 models and members as the PUE/WUE ensemble.

## Part B: model code vs. the papers (local)

**B1. Every parameter, with units.** Re-verify all 92 fields of `Case1Config`, `Case2Config` and `Case5Config` against Table B.1 (pages 8–9 of the 2022 PDF). Each value must be the midpoint. Each conversion must also be right: kPa to Pa, % to fraction (`ups_e`, `pd_lr`, `l_percentage`, `shr`, `windage_p`), RH kept in %, and `pcop` as a fraction. Also check that every Table B.1 row that applies to Cases 1, 2 and 5 has a field, and that no field is missing from the table.

**B2. Equations vs the paper.** Check each energy and water term in `pue_wue_case1`, `case2` and `case5` (in both modules) against the equations in the 2022 paper and its Appendix A: UPS and distribution losses, lighting, fans, pumps, chiller power from COP, cooling-tower evaporation, blowdown from cycles of concentration, drift from windage, humidification and adiabatic water, and the PUE and WUE denominators (IT energy). Note any term in the paper that the code omits, and any term in the code the paper doesn't have.

**B3. The Case 1 indexing fix.** The module says the published `PUE_WUE_AE_Chiller` reads the wrong indices from `Air_side_economizer` and that the port corrects it. Verify the corrected mapping independently, against the documented return order and the physics. Then check that the economizer branch logic (three temperature regimes and their humidity sub-cases) matches the paper's description of Case 1 control.

**B4. Psychrometrics.** Using CoolProp locally, test the fast module's humidity ratio, enthalpy, dewpoint, wet-bulb and specific volume on a grid covering Oregon conditions: −15 to 42 °C, 5–100% RH, 88–102 kPa. Report the maximum errors. Separately, check `elevation_adjusted_pressure_pa` against the standard barometric formula.

**B5. Daily pairing.** Each day's PUE is the mean of a hot pairing (tasmax with hursmin) and a cool pairing (tasmin with hursmax). Find out what weather resolution the paper uses (hourly TMY3?). Check that `methodology.md` describes the pairing correctly. It is a methodological choice, and it changes how often economizer thresholds are crossed.

**B6. Wet-bulb for the figures.** `fig_show_everything.py` recomputes afternoon wet-bulb with `_wetbulb_from_rh` at facility pressure. Confirm it's the same function and pressure the model uses, and that "afternoon" means the hot pairing.

**B7. 2020 model.** Determine whether anything in the chapter (text, a figure, the validation claim) depends on `pue.py` or the 2020 paper's model. If it does, run the same Table A.1 parameter check on `pue.py`.

## Part C: analysis and aggregation code vs. the stated methods (local)

**C1.** Re-derive by hand one Table 1 row and one Figure 4 month: CP, mid-century, Case 1 PUE, from the per-model cache and `facility_month_summary.csv`. Confirm the order of operations: annual mean, then period mean, then the change per facility and model, then the ecoregion mean within each model, then the ensemble mean and the IQR.

**C2.** Check `facility_month_summary.csv`, which pools all days, against an annual-then-monthly calculation for one model. Make sure the difference is negligible.

**C3.** Check the facility inventory:
- the `state_abb == "OR"` filter
- duplicate IDs or coordinates
- the 109 count
- the 21 facilities with no name
- ecoregion assignment for facilities near polygon edges (distance to boundary < 1 km)

**C4.** Check that `gcm_manifest.csv`, the heat ensemble and the PUE/WUE ensemble use the same 20 models and members, and that no model is counted twice. CESM2-LENS in particular could easily carry more than one member.

**C5.** Check every script that reads `model_values_cache.csv`, `facility_month_summary.csv` or `ensemble_delta.csv` for stale-cache risk. After the rerun, all of them must come from the new outputs. Compare file timestamps and the `run_manifest.csv` contents.

## Part D: text vs. code and outputs (local)

**D1. Methods text.** Go through `docs/results-figures/methodology.md` sentence by sentence and record every factual claim with the code or output that supports or contradicts it. Its shared-methods section covers the facilities, climate data, periods, daily PUE/WUE, exceedance days, wet-bulb, aggregation, model spread and ecoregion order. Each figure section covers the inputs, steps and caption. Also record what it **leaves out** that a reader needs, and check each item against the code:
- that equipment parameters are held at the Table B.1 midpoints, and the 2026-10-08 correction
- the chiller COP model and its training range
- the elevation-adjusted pressure formula
- the exceedance threshold's source (ASHRAE A2) and whether it is `>` or `>=`
- which ensemble members are used
- that facilities sharing a grid cell have identical values
- that uncertainty reflects the model ensemble only, not equipment parameters

**D2. Captions.** Check each `*_caption.txt` against `methodology.md` and the current scripts, including the n values, the periods, the band definition and the colors.

---

## Stopping point

After Part A (by the end of 2026-10-09), report the server-side findings, and say for each whether it needs a rerun and how long that would take. Then do Parts B–D and report the full findings table. Don't fix anything until Laura has reviewed the findings.
