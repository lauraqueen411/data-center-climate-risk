# Task: export the values needed to write revised draft 3 (run locally, final outputs)

Write `docs/results-figures/draft3_values.md`, plus CSVs for the larger tables. Use the final, audited outputs and the conventions in `methodology.md`:
- annual means
- change = future minus historical, per facility and model
- ecoregion mean of facilities within each model
- ensemble mean, with the interquartile range across the 20 models

Use these decimals throughout: exceedance days 1, PUE 3 (4 for differences under 0.01), WUE 2 (3 for differences under 0.1), temperature 1.

## A. Table 1 → `table1.csv` (Case 1)

**Rows:** for each ecoregion (ECSF, BM, CP, KM, WV), one row for 2045–2074 and one for 2075–2100.

**Columns:**
- n facilities
- historical annual exceedance days
- historical annual mean afternoon wet-bulb (°C)
- historical annual mean Case 1 PUE and WUE
- change in exceedance days (mean, p25, p75)
- change in PUE (mean, p25, p75)
- change in WUE (mean, p25, p75)

## B. Cases 2 and 5 → `table_cases25.csv`

The same layout as A for PUE and WUE (means plus p25 and p75), with historical values.

## C. Facility-level facts (Case 1)

- **Historical PUE:** the facility minimum and maximum, and which ecoregions they're in.
- **Mid-century change in PUE:**
  - the facility minimum and maximum overall
  - the minimum and maximum within each ecoregion
- **The facility with the smallest mid-century change:** ID, ecoregion, grid cell, historical wet-bulb.
- **The number of distinct grid cells:** in total, and per ecoregion.

## D. Model agreement (Case 1, PUE and WUE, mid- and end-century)

Count the models (out of 20) in which:
- the WV change exceeds the CP change
- the WV change exceeds the KM change
- the CP change exceeds the KM change

## E. Monthly (Figure 4) → `monthly.csv`

**For CP and WV:**
- the mid-century change in monthly mean PUE for Cases 1, 2 and 5, all 12 months (mean, p25, p75)
- the mid-century change in monthly mean Case 1 WUE, all 12 months
- historical monthly mean Case 1 PUE, all 12 months

In the markdown, state:
- the peak month and value for each case and region
- whether CP's July change is a local minimum
- CP's hottest month by historical mean Tmax
- CP's highest month of historical Case 1 PUE

## F. Humidity setpoint sensitivity (Figure 5) → `rh_sensitivity.csv`

**By ecoregion, at `rh_up` 60% and 77.5%:**
- historical Case 1 PUE and WUE
- mid- and end-century change in PUE and WUE (mean, p25, p75)
- the difference between the two setpoints for each of these

**Also include:**
- confirmation that `rh_up` is the only parameter that differs between the two runs
- the model agreement (WV vs CP) at 60%
- cooling-mode shares, if the model exposes them: the share of June–September time steps in supplemental-chiller mode for CP and WV, historical and mid-century, at both setpoints

## G. Warming

By ecoregion, the mid- and end-century change in annual mean daily Tmax: the ecoregion mean, plus the facility minimum and maximum.

## H. Methods facts

- The number of equipment parameters held at Table B.1 midpoints, for each case, with a list.
- The wet-bulb calculation used: the function name and method.
- The 20 models and the member used for each.
- Confirmation that every parameter is at its Table B.1 midpoint except `rh_up` in the 60% run.

## I. Files to hand over

Copy these into one folder:
- the final PNGs of Figures 2–5
- `methodology.md` with final captions
- the CSVs above
- `draft3_values.md`