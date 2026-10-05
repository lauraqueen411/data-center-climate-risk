# Task: build the "show everything" results figure for the heat risk chapter

## Context

The heat risk chapter is dropping the relative-change figure (statewide change factors). It is being replaced by one figure that shows every cooling case, every future window and every ecoregion, so the regional pattern is visible. The layout is modelled on Lei and Masanet (2022), Fig. 5. A mockup with placeholder values is in `show-everything-figure-mockup.png`. Match its layout, but use real data.

**Scope:** this figure only. Don't change other figures, the pipeline outputs or the chapter text. If something below can't be done with the data as stored, stop and report it. Don't approximate.

## Step 1: inspect, then report before plotting

Find where the following are stored. Report the file paths, variable names, dimensions and the averaging already applied.

1. **Daily PUE and WUE** for each facility, each of the 20 models, Cases 1, 2 and 5, and all three periods.
2. **Annual threshold-exceedance day counts** for each facility, each model and each year.
3. **Historical afternoon wet-bulb** for each facility (from Tmax, RHmin and elevation-adjusted pressure).
4. **Ecoregion assignment** for each facility.
5. **LOCA2 grid cell** for each facility.

If only ensemble means are stored, say so. The shaded band needs per-model values. Annual means per facility, per model and per year are enough for this figure.

## Step 2: compute (conventions are fixed; don't change them)

**Periods**
- historical: 1985–2014
- mid-century: 2045–2074
- end-of-century: 2075–2100

**Metrics** (per facility, per model, per period):
- Exceedance days: the annual count of days with Tmax > 35°C, averaged over the years in the period.
- PUE and WUE: the annual mean of daily values (all days), averaged over the years in the period. Cases 1, 2 and 5 only. Don't use May–October averages anywhere in this figure.
- Wet-bulb: the annual mean of daily afternoon wet-bulb, 1985–2014, then the mean across models. Set this with a flag, `WETBULB_WINDOW = "annual"`.

**Change:** future period minus historical, per facility and per model.

**What gets plotted**
- **Facility point:** the mean across the 20 models of that facility's change.
- **Ecoregion mean line:** the mean of the facility points in that ecoregion, with equal weight per facility.
- **Shaded bar:** for each model, compute the ecoregion mean of facility changes. The bar shows the spread of those 20 values, chosen with a flag at the top of the script: `BAND = "iqr"`.
  - `"iqr"`: 25th–75th percentile. This is the band the Seventh Oregon Climate Assessment precipitation chapter shades (Pierce and Cayan 2025).
  - `"p05_p95"`: 5th–95th percentile. This is what the Sixth Oregon Climate Assessment reported.
  - `"minmax"`: the full range across models.
  
  Use `np.percentile` with default linear interpolation. Write all three in the CSV whatever is plotted, and take the legend and caption wording from the flag.

**Ecoregion order:** ascending by ecoregion-mean historical annual afternoon wet-bulb. Compute the order; don't hard-code it. The draft's May–October values gave the order Eastern Cascades Slopes and Foothills, Blue Mountains, Columbia Plateau, Klamath Mountains, Willamette Valley. Report whether annual wet-bulb gives a different order. Include only the five ecoregions with facilities.

## Step 3: plot

**Layout** (matplotlib, gridspec, about 7.5 × 5.5 in, minimum font size 7 pt):
- Panel (a): change in exceedance days. A tall panel on the left spanning both rows; the same for every case.
- Panels (b–d): change in PUE for Cases 1, 2 and 5. Top row.
- Panels (e–g): change in WUE (L/kWh) for Cases 1, 2 and 5. Bottom row.
- Share the y-axis within each row (`SHARE_Y_ROWS = True`, as a flag at the top of the script). Hide y tick labels on the inner panels of each row.
- No x tick labels on the top row.
- Panel titles, top-left: "(a) Exceedance days (all designs)", "(b) PUE, Case 1", and so on.

**Each ecoregion** sits at integer x, with two strips:
- mid-century at x − 0.19, colour `#6da7ec`
- end-of-century at x + 0.19, colour `#184f95`

**Each strip has three layers**, back to front:
1. **Shaded bar:** the `BAND` range, width 0.26, the strip colour at alpha 0.22, no edge, `zorder=1`.
2. **Facility points:** jittered uniformly within ±0.075 of the strip centre with a fixed random seed. About 12 pt² markers in the strip colour, with a white edge of linewidth 0.6, `zorder=3`.
3. **Ecoregion mean:** a horizontal line spanning ±0.12, colour `#1f1f1d`, linewidth 1.6, round caps, `zorder=4`.

**Axes and labels**
- Horizontal zero line in `#6b6a66`, linewidth 0.8.
- Light y-grid only (`#e4e3df`). Remove the top and right spines; draw the other spines and the ticks in `#6b6a66`.
- X tick labels, three lines each: short ecoregion name ("Cascades", "Blue Mtns", "Columbia", "Klamath", "Willamette"), then wet-bulb to one decimal ("15.6 °C"), then "n = 33".
- Y labels:
  - (a) "Change in annual threshold-exceedance days"
  - top row "Change in annual mean PUE"
  - bottom row "Change in annual mean WUE (L/kWh)"

**Legend**, one row below the panels, no frame, four entries:
- Mid-century (2045–2074)
- End-of-century (2075–2100)
- Ecoregion mean of facilities
- [band name from `BAND`] of ecoregion mean across 20 models

**Don't** add a figure title or the mockup banner, and don't draw error bars per facility.

**Save** PNG at 300 dpi and PDF, wherever the repo keeps figures. Name them `fig_show_everything.*`.

## Step 4: outputs and checks to report back

1. **The script**, following repo conventions (for example `scripts/fig_show_everything.py`), runnable end to end.
2. **Figure files**: PNG and PDF.
3. **Tidy CSVs** of everything plotted:
   - facility level: facility, ecoregion, case, metric, period, change, grid cell
   - ecoregion level: ecoregion, case, metric, period, mean, p25, p75, p05, p95, min, max, n facilities, n distinct grid cells, wet-bulb
4. **A short printed summary:**
   - the ecoregion order with wet-bulb, n facilities and n distinct grid cells
   - the ecoregion mean and the plotted band for each panel and period
5. **Reconciliation check.** The draft's Table 1 has Case 1 mid-century changes of 0.017 PUE / 0.21 L/kWh WUE for the Columbia Plateau and 0.02 / 0.27 for the Willamette Valley. These look like annual means. Confirm they match the new annual values, and flag any that don't.
6. **Model agreement.** Count the models, out of 20, in which the Willamette Valley change exceeds the Columbia Plateau change, for each case, metric and period, using annual means. The draft's agreement counts may have used May–October values.
7. **A look at the rendered PNG.** Check for label collisions or clipped markers, and list anything that looks off: overplotting, a strip with a single value, a range that looks implausible.

## Caption (reference; the figure should match it)

Projected change in (a) annual number of threshold-exceedance days, (b–d) power usage effectiveness (PUE), and (e–g) water usage effectiveness (WUE, L/kWh) for cooling-system Cases 1, 2, and 5. Changes are between the historical period (1985–2014) and mid-century (2045–2074; light blue) or end-of-century (2075–2100; dark blue). PUE and WUE are annual means of daily values. Ecoregions are ordered by historical annual mean afternoon wet-bulb temperature, shown below each name with the number of data centers (n). Points are facilities (mean of 20 climate models), jittered horizontally; facilities within the same LOCA2 grid cell have identical values. Black lines are ecoregion means, and shaded bars are the interquartile range (25th–75th percentile) of the ecoregion mean across the 20 models. [Update this phrase if `BAND` changes.]
