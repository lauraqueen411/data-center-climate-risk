# Task: Table 1 and a number check for the heat risk chapter Results (run locally)

Use the analysis outputs behind Figures 2–5 (`outputs/ensemble_pue_wue_facilities/analysis/`). Follow the conventions in `docs/results-figures/methodology.md`:
- annual means
- change = future minus historical, per facility and model
- ecoregion mean of facilities within each model
- ensemble mean, with the interquartile range across the 20 models

Don't change any figures.

## 1. Table 1 → `docs/results-figures/table1.csv`

**Rows:** for each ecoregion in wet-bulb order (ECSF, BM, CP, KM, WV), one row for 2045–2074 and one for 2075–2100.

**Columns**

| Column | Content | Decimals |
|---|---|---|
| ecoregion | name | – |
| n_facilities | count | – |
| hist_exceedance_days | historical annual exceedance days | 1 |
| hist_wetbulb_C | historical annual mean afternoon wet-bulb (°C) | 1 |
| period | 2045–2074 or 2075–2100 | – |
| d_exceedance | change in exceedance days per year | 1 |
| d_exceedance_p25 | 25th percentile across models | 1 |
| d_exceedance_p75 | 75th percentile across models | 1 |
| d_pue | change in Case 1 annual mean PUE | 3 |
| d_pue_p25 | 25th percentile across models | 3 |
| d_pue_p75 | 75th percentile across models | 3 |
| d_wue | change in Case 1 annual mean WUE (L/kWh) | 2 |
| d_wue_p25 | 25th percentile across models | 2 |
| d_wue_p75 | 75th percentile across models | 2 |

Also write `table1_cases.csv` with the same layout for Cases 2 and 5, for reference only.

## 2. Number check → `docs/results-figures/results_number_check.csv`

Compute each value below exactly. Use the columns: item, draft value, computed value, match (yes/no at the stated precision). The draft values were measured from figure images, so small differences are expected. Report every mismatch.

**3.1**
- Historical exceedance days per year: CP 19.2, KM 17.6, BM 7.8, WV 4.1, ECSF 1.3.
- Change in exceedance days, mid-century: CP 30.4, WV 10.9. Confirm WV has the smallest increase of the five ecoregions.
- Change in exceedance days, end-of-century: CP 48.9, WV 22.6.
- Case 1 change in PUE, mid-century: WV 0.021, KM 0.019, CP 0.017, BM 0.011, ECSF 0.011.
- Case 1 change in WUE, mid-century: WV 0.27, KM 0.24, CP 0.21, BM 0.13, ECSF 0.13.
- Case 1 change in PUE, end-of-century: WV 0.037, CP 0.025.
- **Model agreement:** the number of models (out of 20) in which the WV change in Case 1 PUE exceeds the CP change, at mid- and end-century. Also report this for WUE, and for Cases 2 and 5.
- Range of facility-level Case 1 change in PUE within WV, mid-century: 0.014–0.023.
- Case 2:
  - change in PUE below 0.006 in every ecoregion at end-of-century
  - same ecoregion order as Case 1 (yes/no, both periods)
  - change in WUE below 0.01 L/kWh everywhere
- Case 5:
  - ecoregion change in PUE of 0.006–0.007 at mid-century and 0.010–0.012 at end-of-century
  - change in WUE below 0.03 L/kWh everywhere

**3.2**
- Facility range of historical Case 1 PUE: 1.125–1.139.
- Facility range of the mid-century change in Case 1 PUE: 0.011–0.023.

**3.3** (monthly change at mid-century, ecoregion mean)
- Case 1 PUE:
  - WV July 0.060, August 0.066, October 0.013, November 0.015
  - CP June 0.043, July 0.032, September 0.049
- Case 1 WUE peaks: WV August 0.82 L/kWh; CP September 0.63 L/kWh.
- Case 5 PUE peak month and value for CP and WV (draft: July–August, about 0.011).
- Case 2 PUE August value: WV 0.011, CP 0.006.

**3.4**
- Mid-century change in annual mean daily Tmax, as facility ranges:
  - ECSF, BM and CP: 2.7–2.8 °C
  - KM and WV: 2.3–2.5 °C
- Historical annual afternoon wet-bulb ranges:
  - WV facilities: 10.7–11.7 °C
  - BM 8.1 °C and ECSF 6.8 °C
- Confirm no facility falls between 8.1 and 9.8 °C.
- Confirm the WV facility with the smallest Case 1 change in PUE (draft 0.014) has the lowest wet-bulb in WV. Report its ID and location.

**Methods**
- Confirm there are 22 distinct grid cells.
- Report the exact number of equipment parameters held at the Lei and Masanet (2022) Table B.1 midpoints for each of Cases 1, 2 and 5.

## 3. Saturation check (report only)

For CP and WV, report historical monthly Case 1 PUE (ensemble mean, 1985–2014) for each month, and say which month is highest. If the model exposes cooling modes, also report the monthly share of time steps in supplemental-chiller mode, both historically and at mid-century.

The question is whether CP's July dip in the mid-century change coincides with July already being near its maximum historically.

## Output

Print a short summary listing:
- every mismatch
- the model-agreement counts
- the parameter counts
- the saturation result