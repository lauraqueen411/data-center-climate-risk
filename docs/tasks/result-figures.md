# Task: build the four Results figures for the heat risk chapter (run locally)

This runs on the local clone of the repo. The four figures tell one story in order:

- **Figure 2:** where risk changes and for which cooling design ("show everything")
- **Figure 3:** acute risk grows where it already is, while chronic risk shifts to the Willamette Valley (maps)
- **Figure 4:** when the change happens (monthly)
- **Figure 5:** why (mechanism)

Work through one figure at a time. Pause after each one so Laura can review it.

## Step 0: check the data first (don't plot yet)

Check whether the local clone has what each figure needs. Report a table with one row per figure: the files found, the variables, the dimensions, and what's missing. If anything is missing, stop and list exactly what to export from the server. Don't approximate.

| Figure | Needs (per facility × model, all 20 models) |
|---|---|
| 2 | annual PUE and WUE (Cases 1, 2, 5) by year; annual exceedance-day counts by year; historical annual afternoon wet-bulb; ecoregion; grid cell |
| 3 | as for Figure 2, plus Level III ecoregion polygons and the Oregon boundary (the same ones the facility map uses) |
| 4 | monthly mean of daily PUE and WUE (Cases 1, 2, 5) for each period |
| 5 | annual mean daily Tmax by year (historical and mid-century); historical annual afternoon wet-bulb; annual PUE for Cases 1 and 5 |

## Shared conventions (all figures)

**Style:** use the shared style module from the Figure 2 restyle (`figstyle.py`, or whatever was created). If the restyle in `task-show-everything-style.md` hasn't been applied yet, apply it first. Figure 2 is then the style reference for the other three.
- Figures are 6.5 in wide, with large text matching the other chapter figures.
- Black ink, a plain frame, no gridlines, no rounded marks.

**Periods**
- historical: 1985–2014
- mid-century: 2045–2074
- end-of-century: 2075–2100

**Change:** future minus historical, per facility and model.

**Averaging:** annual means of daily values only. Use no May–October averages.

**Plotted values:** the ensemble mean across the 20 models.

**Bands:** the interquartile range (25th–75th percentile) across models of the plotted aggregate (ecoregion or statewide mean), at alpha 0.3.

**Period colors** (Okabe–Ito): mid-century `#0072B2`, end-of-century `#D55E00`.

**Ecoregions:** order by historical annual wet-bulb, and label with ECSF, BM, CP, KM, WV.

## Figure 2: "show everything" (`fig2_show_everything`)

This figure already exists. Rename the output files only, and check it matches the conventions above.

## Figure 3: maps (`fig3_maps`)

**Layout:** 2 × 2 panels.
- (a) historical annual exceedance days
- (b) historical annual mean PUE, Case 1
- (c) change in exceedance days, mid-century
- (d) change in annual mean PUE, Case 1, mid-century

Label the row of change panels "Change, 2045–2074 minus 1985–2014".

**Fill and markers**
- Fill each ecoregion that has facilities with its ecoregion mean. Show the other ecoregions as outlines with a light grey fill.
- Plot facilities as small dots colored by their own value, with a thin black edge, so the spread within the Willamette Valley stays visible.
- Label the five ecoregions with their abbreviations.

**Colors**
- Exceedance days: one sequential palette, such as `Oranges`.
- PUE: a different sequential palette, not blue, to avoid a clash with the mid-century color. Use `Purples` or `cividis`.
- Give each panel its own horizontal colorbar with units.

Reuse the projection and boundary code from the existing facility map.

## Figure 4: monthly change (`fig4_monthly`)

**Layout:** (a) PUE and (b) WUE.
- x: January to December.
- y: monthly change in the statewide mean of all 109 facilities. Don't write "fleet" anywhere.

**Lines**
- Color shows the period (mid-century and end-of-century).
- Line style shows the case: Case 1 solid, Case 2 dashed, Case 5 dotted.
- Each line has its interquartile-range band.

**Regional panels:** add (c) PUE and (d) WUE for Case 1 only.
- These show CP and WV separately: the mean across that ecoregion's facilities, CP solid and WV dashed.
- Color by period, and give each line its band.
- Use a shared y-axis with (a) and (b) respectively, so the regional peaks are directly comparable to the statewide ones.

## Figure 5: mechanism (`fig5_mechanism`)

Plot one point per facility, mid-century only. Show ecoregion by marker shape, all in black: ECSF diamond, BM triangle, CP circle, KM inverted triangle, WV square. Put a legend with the abbreviations below the panels.

- **(a)** x: historical annual afternoon wet-bulb (°C). y: change in annual mean daily Tmax (°C). The expected result is a nearly flat cloud.
- **(b)** Same x. y: change in annual mean PUE.
  - Case 1 uses filled markers and Case 5 open markers.
  - Draw an ordinary least-squares line for each case.
  - Annotate each line with its slope (PUE per °C) and R².

## Outputs

For each figure:
- PNG at 300 dpi and PDF
- a CSV of the plotted values
- a short printed summary of the key numbers for the text:
  - Figure 3: historical value and mid-century change for each ecoregion
  - Figure 4: the peak month and its value for each case and period statewide, and for Case 1 in CP and WV
  - Figure 5: the range of Tmax change across facilities, and the slope, R² and p-value for each case

After all four figures, also report: the number of models (out of 20) in which the Willamette Valley change exceeds the Columbia Plateau change, for Cases 1, 2 and 5, PUE and WUE, at mid- and end-century.

Look at each render at 6.5 in wide, check for label collisions and legibility, and report anything that looks off.