# Code Development Plan: Heat & Flood Risk Analysis
### Oregon Data Center Physical Climate Risk — `data-center-climate-risk` repo

## Current status (from reading the actual repo, not guesswork)

**Shared foundation — mostly built.** `scripts/plot_im3_data_centers.py` loads and maps the IM3/PNNL data center atlas cleanly. `scripts/build_dependency_table.py` ("Deliverable 2") already does real work: it links each Oregon data center to nearby HIFLD substations within a service radius, walks the HIFLD transmission graph (`networkx`) to find connected substations, and uses that to identify EIA-860 generators within a proxy radius — a genuine spatial proxy for the Rinaldi/Peerenboom/Kelly dependency chain the methods doc describes. It's a solid base for both 1.3 and 2.5. Two caveats: it also carries fiber-provider and public-water-source columns that belong to an earlier, broader project scope — the current methods doc explicitly defers water and fiber, so that code should be set aside rather than extended for now. And it samples heat rasters only at *generator* locations, not at substations/transmission lines or at the facilities themselves.

**1.1 Acute heat risk — CLOSED OUT (2026-08-21): window-resolver bug fixed, outputs rerun and verified.** `src/climate_risk_dc/climate/heat_indices.py` + `scripts/run_heat_processing.py` load LOCA2 tasmax/tasmin, convert Kelvin→Celsius with an explicit sanity check, and compute annual CDD, annual Tmax max, threshold-exceedance days (35°C, matching the ASHRAE Class A2 upper bound), and heatwave event counts. Outputs are GeoTIFFs, quicklook PNGs, and an Oregon-wide sanity summary, now for **three periods**: historical (1985-2014), mid-century (2045-2074), end-of-century (2075-2100). `scripts/extract_facility_heat_indices.py` samples all twelve rasters at each Oregon DC's own coordinates (`outputs/month2_heat/facility_heat_indices.csv`, 109/109 facilities, both future-period sanity checks show 0 negative deltas).

**Bug fixed (2026-08-21), see `docs/tasks/0.2-fix-loca2-window-resolver-and-time-horizons.md`:** `run_heat_processing.py`'s old `_loca2_file` future-file fallback picked a single file whose window ended `>= 2099`, which silently returned only the `2075-2100` file for the old `2070-2099` configured period — actually running on **2075-2099 only (25 years)**, silently dropping 2070-2074. Fixed by moving file resolution into a shared `src/climate_risk_dc/climate/loca2_io.py` module (`resolve_loca2_file`) that requires the requested period to be fully contained in exactly one native-window file and raises a clear error otherwise — no more silent "closest file" fallback. Re-verified against the HPC (`ls`, 2026-08-21): ssp585 tasmax/tasmin are split into `2015-2044`/`2045-2074`/`2075-2100`, confirming the original diagnosis. Spot-check for Google's The Dalles facility: corrected end-of-century (2075-2100) `tmax_threshold_days` = 84.85 vs. old truncated (2075-2099) = 84.56 — a small, expected +1-year difference; mid-century (2045-2074) = 59.33, a genuinely new number that didn't exist before (old pipeline never loaded that file). All `outputs/month2_heat/` artifacts (rasters, PNGs, summaries, `facility_heat_indices.csv`) regenerated with the fix; old buggy-window files (`*_1950_1979`, `*_2070_2099`) removed. Unit tests added in `tests/test_loca2_io.py`, including a regression guard asserting the old `2070-2099` period now raises rather than silently resolving.

`geo.get_facility_elevation_m()` + `heat_indices.elevation_adjusted_threshold()` exist and are unit-tested but are **not wired into real facility data** — no DEM exists anywhere in the repo yet (real DEM acquisition is a follow-up, not a silent gap; the placeholder ISA lapse rate of 6.5°C/km used in `elevation_adjusted_threshold` also still needs ASHRAE confirmation). The single-GCM (ACCESS-CM2, r1i1p1f1) limitation is documented in `run_heat_processing.py`'s docstring; whether to extend to a multi-model ensemble remains an open scope/cost decision for a human.

**1.2 Chronic heat risk (PUE) — BUILT (2026-08-21).** `src/climate_risk_dc/climate/pue.py` ports the Lei & Masanet (2020) AE-Chiller (airside economizer + adiabatic cooling) model from the paper's own reference implementation, with a `PUEModelConfig` dataclass defaulted to Table A.1 midpoints (corrected against the paper directly, not just the GitHub sensitivity notebook, for `l_to_g`). Regression-tested against the source repo's own worked example (`tests/test_pue.py`): reproduces `PUE ≈ 1.0746922969797756` exactly. Humidity comes from the LOCA CMIP6 companion dataset (Pierce & Cayan), not MACA — resolved/loaded via `src/climate_risk_dc/climate/loca2_io.py` (`humidity_path_config`, `resolve_loca2_file`, `load_hursmax_hursmin_pct`), which required a fix to avoid materializing the humidity archive's full CONUS-wide grids (~40GB per file, unlike the pre-clipped ~3GB LOCA2 temperature files) — see that module and `scripts/run_pue_processing.py`'s Assumptions logs. `scripts/run_pue_processing.py` computes three-period (historical 1985-2014, mid-century 2045-2074, end-of-century 2075-2100) facility-level PUE for all 109 Oregon facilities using the documented Tmax+RHmin / Tmin+RHmax daily pairing (Section 1.2's stated approximation for the model's hourly-vs-daily resolution gap — see that script's docstring), writing both a period-mean summary (`outputs/month2_pue/facility_pue.csv`) and the full daily per-facility time series (`outputs/month2_pue/facility_pue_timeseries.csv`, ~3.4M rows) so the expensive per-day model evaluation doesn't need to be rerun. Result: Oregon-wide mean PUE 1.1199 (historical) → 1.1209 (mid-century) → 1.1225 (end-of-century), all facilities within the plausible [1.0, 3.0] sanity range and 0 negative deltas under SSP585. **Open items:** the hourly-vs-daily approximation is a genuine, documented gap against how the model was validated (not just a minor detail); the Prineville/Dalles comparison in the source paper is a plausibility anchor for facility-calibrated parameters, not an apples-to-apples validation of this repo's generic Table A.1 defaults; `p_atm_pa` is fixed at sea-level standard atmosphere for every facility (same no-DEM limitation as 1.1's elevation adjustment).

**1.3 Infrastructure heat dependency — partially started.** The dependency graph exists (see above) and generators get heat values attached. The actual derating formulas (IEEE C57.91 for transformers, IEEE 738 for transmission lines, ISO 8528-1 for generators) don't exist yet, and substations/lines aren't sampled for heat at all yet.

**2.1–2.5 Flood risk — essentially unstarted.** Two thin exploratory notebooks: `JRI_flood_explore.ipynb` opens a single file (`glofas_rp100_oregon.tif`) and plots it. `ISIMIP_flood_explore.ipynb` has only boilerplate imports, no logic. No change-factor scaling, no facility extraction, no 3DEP elevation, no Hazus, no AAD, no infrastructure flood exposure.

*Resolved: the GloFAS-vs-JRC naming question.* Confirmed via the source portal (`jeodpp.jrc.ec.europa.eu/.../CEMS-GLOFAS/flood_hazard`) — this is the **JRC Global River Flood Hazard Maps** (v2.1, Baugh, Colonese, D'Angelo et al. 2024), distributed under Copernicus Emergency Management Service (CEMS) branding that folds in "GLOFAS" because it builds on the GloFAS river-flow framework (LISFLOOD/LISFLOOD-FP). It's the same dataset the methods doc describes — 90m resolution, seven return periods (RP10/20/50/75/100/200/500), water depth in meters — not a different product. `glofas_rp100_oregon.tif` is correctly named after its `RP100` source folder; no rename needed. Worth doing: cite Baugh et al. (2024) alongside ref [16] in the methods doc's flood section, since it's the more precise source for the JRC depth grids than the current generic description.

**Environment.** Conda env (`environment.yml`) already has the full stack needed: geopandas, rasterio, xarray, dask, h3, networkx, cartopy, mapclassify. Nothing extra to install for the work below.

## Decisions locked in for this plan

- **1.3 scope:** build the lightweight version now (heat-sampling extended to substations/transmission lines, simple "extreme heat day coincides with a dependency" flag — no capacity-loss percentages), structured so the full IEEE/ISO derating math can be swapped in later.
- **Sequencing:** flood work starts now, in parallel with closing out heat's remaining gaps — not gated behind finishing heat first. Phases 1 and 2 can run concurrently; Phase 0 is a short prerequisite for both.
- **Time horizons (2026-08-20, revised same day):** switching from a single future window (2070–2099, previously mislabeled "mid-century" — it isn't) to two future windows, snapped to the LOCA2 archive's own native ssp585 file windows rather than round decades (avoids needing any cross-file concatenation): **mid-century (2045–2074, 30 years)** as the primary analysis period, since a data center's own operational/mechanical lifespan is a few decades, and **end-of-century (2075–2100, 26 years)** kept as a conservative bookend, since the power infrastructure a facility depends on (substations, transmission, generation) is planned and operated on much longer horizons than the facility itself. Applies to 1.1, 1.2, and (once flood work starts) 2.x. See `docs/tasks/0.2-fix-loca2-window-resolver-and-time-horizons.md`.

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
- [x] Sample the heat rasters (three periods × CDD/threshold-days/heatwave-events) directly at each Oregon DC's own coordinates (reuse `sample_raster_at_points`), producing the real per-facility exceedance-day/heatwave table — `scripts/extract_facility_heat_indices.py` → `outputs/month2_heat/facility_heat_indices.csv`, 109/109 Oregon facilities, all deltas non-negative for both mid-century and end-of-century
- [x] Add elevation-based threshold adjustment at high-elevation sites (share the 3DEP loader being built in Phase 2.2 rather than duplicating it) — `geo.get_facility_elevation_m()` + `heat_indices.elevation_adjusted_threshold()` built and unit-tested; **not yet wired into real facility data** — no DEM exists in the repo yet (user decision: build+test now, source a real DEM later)
- [x] Document the single-GCM-vs-ensemble decision explicitly (commit to ACCESS-CM2 only with a stated limitation, or extend to a small multi-model set) — limitation documented in `run_heat_processing.py`'s module docstring; whether to extend to an ensemble is still an open scope/cost decision for a human, not resolved here
- [x] **(2026-08-21)** Fix confirmed LOCA2 future-window resolver bug and switch to the three-period structure (historical 1985-2014, mid-century 2045-2074, end-of-century 2075-2100) snapped to native ssp585 file windows — `src/climate_risk_dc/climate/loca2_io.py` (`resolve_loca2_file`, unit-tested in `tests/test_loca2_io.py`), all `outputs/month2_heat/` artifacts regenerated and reverified. See `docs/tasks/0.2-fix-loca2-window-resolver-and-time-horizons.md`.

**1.2 build (from scratch) — DONE (2026-08-21)**
- [x] Create `src/climate_risk_dc/climate/pue.py`: Lei & Masanet AE-Chiller (airside economizer + adiabatic cooling) PUE model, ported from the paper's own reference implementation as documented pure functions over temperature + humidity arrays; `PUEModelConfig` defaults to Table A.1 midpoints
- [x] Add a LOCA CMIP6 humidity companion loader (not MACA — superseded per methods doc, since MACA is CMIP5-era and would mismatch the CMIP6 temperature pipeline's GCM lineage), in `src/climate_risk_dc/climate/loca2_io.py` alongside the shared LOCA2 file resolver
- [x] Create `scripts/run_pue_processing.py` (structured like `run_heat_processing.py`): three-period (historical/mid-century/end-of-century) PUE per facility, writing both a period-mean summary and the full daily per-facility time series
- [x] Validate: regression test against the source repo's own worked example (exact match, `tests/test_pue.py`); Prineville/Dalles compared against the paper's reported figures as a plausibility anchor (not a strict pass/fail bar — see the 1.2 status note above). Ref [12] (Guidi et al.) was not pursued further per the task doc's guidance, since its Oregon utility-level PUE breakdown wasn't retrievable and the source-repo regression + paper's own validation set were judged sufficient.

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