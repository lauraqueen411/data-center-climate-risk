# Code Development Plan: Heat & Flood Risk Analysis
### Oregon Data Center Physical Climate Risk — `data-center-climate-risk` repo

## Current status (from reading the actual repo, not guesswork)

**Shared foundation — mostly built.** `scripts/plot_im3_data_centers.py` loads and maps the IM3/PNNL data center atlas cleanly. `scripts/build_dependency_table.py` ("Deliverable 2") already does real work: it links each Oregon data center to nearby HIFLD substations within a service radius, walks the HIFLD transmission graph (`networkx`) to find connected substations, and uses that to identify EIA-860 generators within a proxy radius — a genuine spatial proxy for the Rinaldi/Peerenboom/Kelly dependency chain the methods doc describes. It's a solid base for both 1.3 and 2.5. Two caveats: it also carries fiber-provider and public-water-source columns that belong to an earlier, broader project scope — the current methods doc explicitly defers water and fiber, so that code should be set aside rather than extended for now. And it samples heat rasters only at *generator* locations, not at substations/transmission lines or at the facilities themselves.

**1.1 Acute heat risk — substantially built.** `src/climate_risk_dc/climate/heat_indices.py` + `scripts/run_heat_processing.py` load LOCA2 tasmax/tasmin, convert Kelvin→Celsius with an explicit sanity check, and compute annual CDD, annual Tmax max, threshold-exceedance days (35°C, matching the ASHRAE Class A2 upper bound), and heatwave event counts, for historical (1950–1979) vs. SSP585 (2070–2099). Outputs are GeoTIFFs, quicklook PNGs, and an Oregon-wide sanity summary. Gaps: no facility-level table (the six raster values still need to be sampled at each DC's own lon/lat — right now only generators get sampled); no elevation-based threshold adjustment at high-elevation sites; locked to a single GCM (ACCESS-CM2, r1i1p1f1) with no documented rationale for not using an ensemble.

**1.2 Chronic heat risk (PUE) — not started.** No MACA humidity loading, no Lei & Masanet model code, anywhere in the repo.

**1.3 Infrastructure heat dependency — partially started.** The dependency graph exists (see above) and generators get heat values attached. The actual derating formulas (IEEE C57.91 for transformers, IEEE 738 for transmission lines, ISO 8528-1 for generators) don't exist yet, and substations/lines aren't sampled for heat at all yet.

**2.1–2.5 Flood risk — essentially unstarted.** Two thin exploratory notebooks: `JRI_flood_explore.ipynb` opens a single file (`glofas_rp100_oregon.tif`) and plots it. `ISIMIP_flood_explore.ipynb` has only boilerplate imports, no logic. No change-factor scaling, no facility extraction, no 3DEP elevation, no Hazus, no AAD, no infrastructure flood exposure.

*Resolved: the GloFAS-vs-JRC naming question.* Confirmed via the source portal (`jeodpp.jrc.ec.europa.eu/.../CEMS-GLOFAS/flood_hazard`) — this is the **JRC Global River Flood Hazard Maps** (v2.1, Baugh, Colonese, D'Angelo et al. 2024), distributed under Copernicus Emergency Management Service (CEMS) branding that folds in "GLOFAS" because it builds on the GloFAS river-flow framework (LISFLOOD/LISFLOOD-FP). It's the same dataset the methods doc describes — 90m resolution, seven return periods (RP10/20/50/75/100/200/500), water depth in meters — not a different product. `glofas_rp100_oregon.tif` is correctly named after its `RP100` source folder; no rename needed. Worth doing: cite Baugh et al. (2024) alongside ref [16] in the methods doc's flood section, since it's the more precise source for the JRC depth grids than the current generic description.

**Environment.** Conda env (`environment.yml`) already has the full stack needed: geopandas, rasterio, xarray, dask, h3, networkx, cartopy, mapclassify. Nothing extra to install for the work below.

## Decisions locked in for this plan

- **1.3 scope:** build the lightweight version now (heat-sampling extended to substations/transmission lines, simple "extreme heat day coincides with a dependency" flag — no capacity-loss percentages), structured so the full IEEE/ISO derating math can be swapped in later.
- **Sequencing:** flood work starts now, in parallel with closing out heat's remaining gaps — not gated behind finishing heat first. Phases 1 and 2 can run concurrently; Phase 0 is a short prerequisite for both.

---

## Phase 0 — Foundation cleanup (do first, ~1–2 days, unblocks everything else)

- [x] Confirm the real source of `glofas_rp100_oregon.tif` — verified as JRC Global River Flood Hazard Maps v2.1 (CEMS-GLOFAS distribution); no rename needed, naming is already correct
- [ ] Leave `exposure/fiber.py`, `exposure/water.py`, `exposure/watersheds.py`, `DC_fiber_exposure.ipynb`, `DC_water_exposure.ipynb`, and `test_public_water_sources.py` untouched — don't extend them
- [ ] Strip or clearly mark the fiber/water columns in `build_dependency_table.py`'s active output path as legacy, so new heat/flood work doesn't tangle with deferred scope
- [ ] Add `config/datasets.yml` covering HIFLD, EIA-860, JRC, ISIMIP, and (new) 3DEP paths
- [ ] Route `run_heat_processing.py` and `build_dependency_table.py` through `config.load_dataset_paths()` instead of hardcoded `root/data/...` paths
- [ ] Extend the dependency table to emit substation and transmission-line identifiers/geometries as first-class per-asset output (not just generator edges) — needed by both 1.3 and 2.5

## Phase 1 — Heat risk close-out (runs alongside Phase 2)

**1.1 close-out**
- [ ] Sample the six existing heat rasters directly at each Oregon DC's own coordinates (reuse `_sample_raster_values`), producing the real per-facility exceedance-day/heatwave table
- [ ] Add elevation-based threshold adjustment at high-elevation sites (share the 3DEP loader being built in Phase 2.2 rather than duplicating it)
- [ ] Document the single-GCM-vs-ensemble decision explicitly (commit to ACCESS-CM2 only with a stated limitation, or extend to a small multi-model set)

**1.2 build (from scratch)**
- [ ] Create `src/climate_risk_dc/climate/pue.py`: Lei & Masanet PUE model, outside-air/evaporative variant, as pure functions over temperature + humidity arrays
- [ ] Add a MACA humidity loader, mirroring the LOCA2 loader pattern in `heat_indices.py`
- [ ] Create `scripts/run_pue_processing.py` (structured like `run_heat_processing.py`): historical vs. future PUE per facility
- [ ] Validate model output against reported Oregon utility PUE figures (ref [12]) before trusting the future projection

**1.3 lightweight build**
- [ ] Extend the Phase 0 per-asset dependency table so substations and transmission-line midpoints get the same heat-raster sampling generators already get
- [ ] Produce a simple coincidence flag per facility: extreme-heat day overlaps with a dependency (substation, line, or generator) — no capacity-loss math
- [ ] Leave a marked extension point (e.g. a `derating.py` stub with IEEE C57.91/738/ISO 8528-1 references in the docstring) for a possible full version later

## Phase 2 — Flood risk build-out (runs alongside Phase 1, starts immediately)

**2.1 Flood hazard data and future scaling**
- [x] Confirm the JRC flood hazard source — verified as JRC Global River Flood Hazard Maps v2.1 (Baugh et al. 2024), 90m, RP10–RP500, distributed via CEMS-GLOFAS; `glofas_rp100_oregon.tif` naming is correct, no rename needed
- [ ] Download/acquire the remaining 6 return periods for Oregon (RP10, RP20, RP50, RP75, RP200, RP500 — only RP100 pulled so far)
- [ ] Confirm/acquire ISIMIP historical + future files
- [ ] Add Baugh et al. (2024) as a citation alongside ref [16] in the methods doc's flood section
- [ ] Build `src/climate_risk_dc/flood/hazard.py`: nearest-cell change-factor computation (ISIMIP future/historical ratio per return period) applied to the JRC baseline

**2.2 Facility-level flood exposure**
- [ ] Build `src/climate_risk_dc/flood/elevation.py`: USGS 3DEP extraction at facility coordinates (shared with Phase 1's elevation adjustment)
- [ ] Extract scaled flood depth at each DC for each return period × time horizon, relative to facility ground elevation
- [ ] Flag facilities exceeding a minimum damage threshold

**2.3 Depth-damage translation**
- [ ] Build `src/climate_risk_dc/flood/hazus.py`: apply the chosen Hazus building-category curve
- [ ] Test against at least one alternate building category for comparison
- [ ] Apply the first-floor-height offset assumption

**2.4 Average annual damage**
- [ ] Implement trapezoidal-rule integration of damage over return-period probability
- [ ] Document the edge-case assumptions (no damage at ~1-year return period; damage held constant beyond the largest available return period)

**2.5 Infrastructure flood dependency**
- [ ] Apply the Phase 2.1/2.2 flood-depth extraction to the Phase 0 per-asset dependency table
- [ ] Flag facilities with any exposed substation, line, or generator (same pattern as the lightweight 1.3 flag, using flood depth instead of heat)

## Phase 3 — Integration, outputs, and QA (after Phases 1 & 2 substantially land)

- [ ] Combine the heat (1.1–1.3) and flood (2.1–2.5) facility-level tables into one master Oregon DC risk table, one row per facility — following the existing `dependency_table_facility_summary.csv` pattern
- [ ] Add bounds/sanity checks for the new outputs (plausible flood-depth ranges, PUE ratios, damage percentages), matching the style of `_log_unit_sanity`
- [ ] Spot-check a couple of well-known facilities by hand (e.g., Google's The Dalles campus) against pipeline output
- [ ] Add unit tests for the new pure-math modules — PUE, Hazus depth-damage, AAD integration, and the lightweight derating flag are all testable on small synthetic inputs
- [ ] Fill in `README.md` (currently a placeholder) with setup instructions and how to run each script

## Open items to keep an eye on

- [x] Resolve the GloFAS-vs-JRC filename question — confirmed same dataset (JRC Global River Flood Hazard Maps v2.1 / CEMS-GLOFAS), no action needed
- [ ] Decide (later, not now) whether the water/fiber code should come back into scope as a future extension
- [ ] Document the single-GCM choice for LOCA2 (ACCESS-CM2) in the methods doc once finalized
