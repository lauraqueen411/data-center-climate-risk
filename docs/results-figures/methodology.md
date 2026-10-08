# Results figures: methodology

Methods behind Figures 2–5 of the heat risk chapter. Each figure is built by one script in `scripts/`, writes a PNG (300 dpi) and PDF to `docs/results-figures/`, and writes its plotted values to `outputs/ensemble_pue_wue_facilities/analysis/<figure>/`.

## Shared methods

**Facilities and ecoregions.** The 109 Oregon facilities come from the IM3/PNNL Open Source Data Center Atlas (v2026.02.09). Each facility is assigned to the EPA Level III ecoregion polygon that contains it. Five ecoregions contain facilities: Eastern Cascades Slopes and Foothills (ECSF, n = 1), Blue Mountains (BM, n = 13), Columbia Plateau (CP, n = 61), Klamath Mountains (KM, n = 1) and Willamette Valley (WV, n = 33).

**Climate data.** Daily maximum and minimum temperature come from LOCA2, and daily maximum and minimum relative humidity from the LOCA CMIP6 humidity companion on the same 1/16° grid. The ensemble is 20 CMIP6 models under SSP3-7.0, one member per model (`outputs/ensemble_pue_wue/gcm_manifest.csv`). Each facility takes the values of its nearest LOCA2 grid cell. The 109 facilities fall in 22 grid cells, and facilities sharing a cell have identical values.

**Periods.** Historical 1985–2014, mid-century 2045–2074, end-of-century 2075–2100.

**Daily PUE and WUE.** The Lei and Masanet (2022) cooling archetypes compute daily PUE and WUE for Case 1 (airside economizer with adiabatic cooling and chiller), Case 2 (waterside economizer with water-cooled chiller) and Case 5 (chiller only). Air pressure at each facility is adjusted for elevation with the standard barometric formula, using the LOCA elevation field. Each day's value is the mean of two runs: a hot pairing of maximum temperature with minimum relative humidity, and a cool pairing of minimum temperature with maximum relative humidity.

**Exceedance days.** An exceedance day is a day with maximum temperature above 35 °C. Annual counts are computed on the LOCA2 grid for each model and period, and the period mean of the annual counts is sampled at each facility's nearest grid cell.

**Afternoon wet-bulb temperature.** Daily afternoon wet-bulb temperature is computed from daily maximum temperature and minimum relative humidity at the facility's elevation-adjusted pressure, using the psychrometric functions of the PUE/WUE model. Each facility's value is the mean of all days in each historical year, averaged over 1985–2014 and then over the 20 models.

**Aggregation.** All values are annual means of daily values, averaged over the years in each period. No seasonal (May–October) averages are used.
1. Per facility and model, the period value is the mean of the annual means.
2. Change is the future period value minus the historical value, per facility and model.
3. The plotted facility value is the mean across the 20 models.
4. An ecoregion or statewide value is the mean of its facilities, with each facility weighted equally.

**Model spread.** Bands show the interquartile range, the 25th to 75th percentile, across the 20 models. The percentiles are taken of the aggregate that is plotted: each model's ecoregion mean is computed first, then the percentiles across models. Bands are drawn at alpha 0.3.

**Ecoregion order.** Where ecoregions are ordered, they run from lowest to highest historical annual afternoon wet-bulb temperature, which is the mean of each ecoregion's facilities: ECSF, BM, CP, KM, WV.

**Style.** Figures are 6.5 in wide and share the settings in `src/climate_risk_dc/figstyle.py`. Mid-century is Okabe–Ito blue (`#0072B2`) and end-of-century vermillion (`#D55E00`). Ecoregion colors are ECSF orange, BM sky blue, CP bluish green, KM black and WV reddish purple.

## Figure 2: change by cooling case and ecoregion

Script: `scripts/fig_show_everything.py`. Output: `fig2_show_everything.{png,pdf}`.

1. For each facility and model, compute the change in annual exceedance days and in annual mean PUE and WUE for Cases 1, 2 and 5, at mid-century and end-of-century.
2. Plot each facility's ensemble-mean change as a point, jittered horizontally within its ecoregion and colored by period.
3. Plot the ecoregion mean of the facility points as a black line.
4. Plot the interquartile range across models of the ecoregion mean as a shaded bar.
5. Order ecoregions by historical annual wet-bulb temperature. Panels in the same row share a y-axis.

**Caption.** Projected change in (a) annual number of threshold-exceedance days, (b–d) annual mean power usage effectiveness (PUE), and (e–g) annual mean water usage effectiveness (WUE, L/kWh) for cooling-system Cases 1, 2, and 5. Changes are between the historical period (1985–2014) and mid-century (2045–2074; blue) or end-of-century (2075–2100; vermillion). Ecoregions are ordered from left to right by historical annual mean afternoon wet-bulb temperature: Eastern Cascades Slopes and Foothills (ECSF; n = 1), Blue Mountains (BM; n = 13), Columbia Plateau (CP; n = 61), Klamath Mountains (KM; n = 1), and Willamette Valley (WV; n = 33). Points are facilities (mean of 20 climate models), jittered horizontally; facilities within the same LOCA2 grid cell have identical values. Black lines are ecoregion means, and shaded bars are the interquartile range (25th–75th percentile) of the ecoregion mean across the 20 models.

## Figure 3: maps

Script: `scripts/fig3_maps.py`. Output: `fig3_maps.{png,pdf}`.

1. For each facility, compute the historical value and the mid-century change in annual exceedance days and in annual mean Case 1 PUE, each as the mean across models.
2. Fill each ecoregion that contains facilities with the mean of its facility values. Fill all other Level III ecoregions light grey.
3. Draw each facility as a dot with a thin black edge, colored by its own value on the same scale as the fill.
4. Give each panel its own color scale, running from the lowest to the highest facility value. Exceedance days use Oranges and PUE uses Purples. The palest 15% of each colormap is left out so that the lowest-valued ecoregion is not confused with the grey fill.
5. Draw the maps in an Albers equal-area projection (central meridian 120.5° W, standard parallels 43° and 45.5° N) with Natural Earth 50m state lines and coastline, the same projection as the facility orientation map.

**Caption.** Historical (1985–2014) values and mid-century change (2045–2074 minus 1985–2014) in (a, c) annual number of threshold-exceedance days (daily maximum temperature above 35 °C) and (b, d) annual mean power usage effectiveness (PUE) for cooling-system Case 1. Values are the mean of 20 climate models (SSP3-7.0). Ecoregions containing facilities are filled with the mean of their facilities: Eastern Cascades Slopes and Foothills (ECSF; n = 1), Blue Mountains (BM; n = 13), Columbia Plateau (CP; n = 61), Klamath Mountains (KM; n = 1) and Willamette Valley (WV; n = 33). Dots are facilities, colored on the same scale; other Level III ecoregions are shown in grey.

## Figure 4: monthly change

Script: `scripts/fig4_monthly.py`. Output: `fig4_monthly.{png,pdf}`.

1. For each facility, model, case and calendar month, take the monthly mean PUE: the mean of all daily values in that month over the period's years, from `facility_month_summary.csv`. Every year contributes the same number of days to a month except February in leap years, so this equals the mean of the yearly monthly means to within that one day.
2. Compute the mid-century change per facility, model, case and month.
3. Average the changes over the CP facilities and over the WV facilities for each model. Plot the mean across models as a line and the interquartile range across models as a band.
4. Draw Case 1 solid, Case 2 dashed and Case 5 dotted, colored by region. Both panels share the y-axis.

**Caption.** Monthly change from 1985–2014 to mid-century (2045–2074) in mean power usage effectiveness (PUE) for facilities in (a) the Columbia Plateau (n = 61) and (b) the Willamette Valley (n = 33), for cooling-system Cases 1 (solid), 2 (dashed) and 5 (dotted). Lines are the mean of the ecoregion's facilities, averaged over 20 climate models (SSP3-7.0); shading is the interquartile range (25th–75th percentile) across models. Both panels share the y-axis.

## Figure 5: warming and PUE response against wet-bulb

Script: `scripts/fig5_mechanism.py`. Output: `fig5_mechanism.{png,pdf}`.

1. For each facility and model, compute the annual mean of daily maximum temperature in each year and average over each period. The annual means come from `facility_annual_tmax.csv`, which was produced on the server from the Case 1 daily files.
2. Compute the mid-century change in annual mean daily maximum temperature and in annual mean Case 1 PUE, per facility and model, then average across models.
3. Plot each facility against its historical annual afternoon wet-bulb temperature, with marker shape and color giving the ecoregion. Both y-axes start at zero.
4. Draw no fit lines. The 109 facilities fall in 22 grid cells and five ecoregion clusters, so they are not independent observations. Least-squares fits over the facilities and over the 22 grid cells are written to `fits.csv` for reference only.

**Caption.** Mid-century (2045–2074 minus 1985–2014) change in (a) annual mean daily maximum temperature and (b) annual mean power usage effectiveness (PUE) for cooling-system Case 1, against historical annual mean afternoon wet-bulb temperature. Each marker is one facility (mean of 20 climate models, SSP3-7.0); marker shape and color give the ecoregion: Eastern Cascades Slopes and Foothills (ECSF; n = 1), Blue Mountains (BM; n = 13), Columbia Plateau (CP; n = 61), Klamath Mountains (KM; n = 1) and Willamette Valley (WV; n = 33). Facilities in the same LOCA2 grid cell share identical values.
