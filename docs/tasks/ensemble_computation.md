# Task for Claude Code: Gridded, Full-Ensemble PUE/WUE Projections for Oregon

## Context (for you, before handing this off)

This replaces Section 1.2 of the Oregon data center physical climate risk study. The current pipeline computes PUE at each of ~109 facility point-locations, using one GCM (ACCESS-CM2) and one archetype (Lei & Masanet's Case 1: airside economizer + adiabatic humidification + water-cooled chiller). This task extends that to the full LOCA2 CMIP6 ensemble, every grid cell in Oregon, and multiple cooling archetypes — closer in spirit to Lei & Masanet's own national study, but future-facing and at grid-cell rather than climate-sub-region resolution.

Two things to reconcile before this lands in the paper, not blocking the build itself: Section 1.1 (heat risk) runs SSP5-8.5, while this runs SSP3-7.0 — that inconsistency needs a sentence in Limitations. And this is a substantially bigger lift than what's currently scheduled for the Sept 29 draft — treat the compute run as the fast part; validation, figures, and a Methods rewrite follow after.

Run Phase 0 first and review its output before starting Phase 2 — it determines how much of the rest needs to change.

---

## Phase 0: Inventory the archetype reference implementation (do this first, report back)

Clone/read `github.com/nuoaleon/Data-Center-Water-footprint` in full — not just the README. Confirmed so far (from `simulation_funs_DC.py`): three pickled `GaussianProcessRegressor` models back three distinct cooling archetypes —

- `COP_2.pkl` — water-cooled chiller, `predict([Twb_oa+AT_CT, Chiller_load])` (wet-bulb + cooling-tower approach temp, chiller load). This is the one already feeding the ported Case 1 AE archetype.
- `COP_DX.pkl` — direct expansion, `predict([T_oa])` — dry-bulb only, univariate, no load or humidity input.
- `COP_AC.pkl` — air-cooled chiller, `predict([T_oa, Chiller_load])` — dry-bulb and load, no humidity input.

Note DX and air-cooled take no humidity input at all in the COP call — they may turn out to have little or no modeled water use, which is a legitimate result to report for the archetype comparison, not a bug to chase. Still to confirm in Phase 0, not yet verified:

- Whether/how each archetype's water balance (cooling tower evaporation/windage/blowdown, adiabatic/space humidification) is computed elsewhere in the file, since DX and air-cooled may skip this entirely — check rather than assume.
- The equipment parameters each needs (the Case-specific column in Lei & Masanet 2022's Table A.1: `pcop` and others) and water-balance parameters (Table B.1) where applicable.
- Cross-reference each of the three against the named Cases in Lei & Masanet (2020, *Energy* 201, 117556) and (2022, *Resources, Conservation and Recycling* 182, 106323) to confirm which published Case each corresponds to. These papers are downloaded in /papers/
- **The training range each GP was fit on.** Each `.predict()` call embeds the training data through its kernel — extrapolating `Twb_oa+AT_CT` or `T_oa` into SSP3-7.0 end-century territory the model never saw behaves however that kernel does outside its envelope (commonly reverting toward the prior mean with inflated variance, but this is kernel-dependent and model-specific — check it, don't assume it). Pull the calibration/validation range from the papers (Table A.1 bounds, Figure 4/6 ranges) and compare against this study's actual SSP3-7.0 mid- and end-century range for `Twb_oa+AT_CT` and `T_oa` over the Oregon grid, per the approach already laid out in `climate-extension-justification.md`.

Report this back before writing the full archetype library. Don't guess at parameter tables from function code alone if the papers' tables are the authoritative source — pull the actual published values.

---

## Phase 1: Data and domain

**Climate inputs** (already on this machine, read-only shared data mount):
- LOCA2 tasmax/tasmin: `/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2/CONUS_regions_split/` — **confirmed pre-split to a Northwest region**, not full CONUS. Read directly from the Northwest split; no CONUS-wide read or subsetting pass needed for temperature (Phase 4's I/O optimization applies to humidity, not this).
- LOCA2-consistent CMIP6 humidity fields: `/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2_humidity/` — **confirmed full CONUS**, no regional pre-split. This is the actual target of Phase 4's Oregon-clipping/subsetting pass — the original reason parallelization/I/O was flagged as needing real attention for this task.
- LOCA elevation/orography field: `/home/server/pi/homes/queenl/projects/climate-risk-dc/data/elevation.LOCA_2016-04-02.nc` — **confirmed**: filename is `LOCA`, not `LOCA2`. Static field, load once and reuse across all GCMs/periods rather than re-reading per combination.

**Ensemble**: LOCA2 has up to 27 CMIP6 GCMs, but SSP3-7.0 availability varies by model — don't assume all 27 have it downscaled. Use however many of the 27 actually have SSP3-7.0 present in the downloaded data, one member per model. Auto-detect which member is present per model from the downloaded files rather than hardcoding an ID — models weren't all downloaded with the same realization/initialization/physics/forcing index. Log the exact member ID used per model, and the final list of which GCMs were actually included vs. skipped for lacking SSP3-7.0, to a manifest file; both need to be reported in the paper's Methods/Table for reproducibility.

**Time periods**: historical 1985–2014 (30yr), mid-century 2045–2074 (30yr), end-century 2075–2100 (26yr — note in output metadata that this window is 4 years shorter than the other two, since it's a real asymmetry in climatological robustness worth flagging in Limitations).

**Domain**: clip to the Oregon state boundary polygon (not a bounding box) — reuse whatever boundary source `scripts/plot_im3_data_centers.py` already uses, don't introduce a second one. Mask out any ocean/no-data cells outside LOCA2's land mask.

---

## Phase 2: Core computation

- Extend the existing `pue.py` (already vectorized for facility-level, single-timeseries Case 1) to operate on full gridded arrays — dims `(lat, lon, time)` — instead of per-facility 1D series. Don't fork a separate module per archetype if the underlying physics functions can share structure; but don't force a shared abstraction Phase 0 shows doesn't fit.
- The chiller COP model is a `GaussianProcessRegressor` per archetype (`COP_2.pkl` / `COP_DX.pkl` / `COP_AC.pkl`), each currently called as `.predict(np.array([...]).reshape(1, n_features))` on one facility-timestep at a time — that pattern is the reference implementation's, written for point-based use. At grid scale, build one `(N, n_features)` array spanning every grid cell × timestep for a given (GCM, period, archetype) and call `.predict()` once — this is the main vectorization opportunity, and it's a real order-of-magnitude change, not just tidying: GP prediction cost scales with kernel evaluations against the training set, so batching amortizes that instead of paying Python-loop overhead per query. Pin the scikit-learn version the pickled GP models were fit under — check the reference repo's own `environment.yml`/`requirements.txt` (or any pinned version in its README) against what `pip show scikit-learn` reports on this machine before running anything, rather than assuming they match; a version mismatch can silently change predictions, not just error out.
- Add elevation-adjusted atmospheric pressure to the chiller/cooling-tower physics using the LOCA elevation field (Phase 1), replacing the standard-atmosphere assumption the facility-level pipeline still uses. This matters much more here than it did for a handful of valley-floor sites — Oregon's grid spans sea level to >3000m in the Cascades.
- Total output files: N_gcms (however many of the 27 actually have SSP3-7.0, per Phase 1) × 1 scenario × N archetypes (N from Phase 0) × 3 time windows — each file spans the full Oregon grid × daily time series for that one window (see Phase 5 for the output spec).

**CLI interface — this is the actual deliverable, not an implementation detail.** Build one script, e.g.:

```
python run_pue_wue_gridded.py --gcms ACCESS-CM2 [CESM2 ...]
```

`--gcms` is the required, primary argument. For every GCM passed, the script loops through every archetype (from Phase 0) and every one of the 3 time windows, writing one output file per (GCM, archetype, period) — see Phase 5 — under the fixed SSP3-7.0 scenario. Archetypes, scenario, and periods are not separately selectable on the command line; the GCM list is the only thing that varies per invocation. Keep the script itself simple: a plain loop over the given GCM(s), then nested loops over archetypes and periods — the concurrency happens outside this script, not inside it (see Phase 3). Within one invocation, process the given GCM(s) **sequentially**, not in parallel internally — this is deliberate, so a run with 1 or 2 GCMs gives a clean, uncontended wall-clock benchmark before committing to the full set (Phase 1 — however many of the 27 actually have SSP3-7.0).

Log wall-clock time for every (GCM, archetype, period) output file as it completes, plus a per-GCM total when that GCM's full set finishes, to a structured log (CSV or JSON lines) rather than only stdout — this is the actual benchmarking data the concurrency decision in Phase 3 depends on.

---

## Phase 3: Parallelization and performance (shared 64-core lab server, no job scheduler)

Two separate levels of parallelism here — don't conflate them:

1. **Within a single GCM's run** (inside one script invocation, processing one GCM at a time per the Phase 2 CLI design): vectorize fully over space, time, and archetypes with numpy/xarray broadcasting — no Python-level loops over grid cells or days — and batch the GP `.predict()` calls per Phase 2. This is what makes one GCM's run fast; it is not multiprocessing.
2. **Across GCMs**: this is on you, by design — you're starting several processes yourself (e.g. one `--gcms <model>` invocation per process, or a small batch of models per process), not asking the script to fan out internally. That keeps the script itself a plain loop (Phase 2) and keeps each process's resource use and timing easy to reason about. The benchmark-then-scale workflow this enables: run `--gcms ACCESS-CM2`, then `--gcms ACCESS-CM2 CESM2`, use the logged per-GCM wall-clock and peak memory from those to decide how many concurrent processes the box can actually absorb, then launch the rest (the full N_gcms from Phase 1, not necessarily 27) in batches.

This is a shared multi-user machine, not a dedicated box: 64 logical cores, 376GB RAM, but ~22 users regularly logged in and load average sitting around 2–3.5 even at idle-ish moments, with other people's jobs (MATLAB, Python) actively running. There is no SLURM or cgroup fair-share enforcing limits — self-impose them instead.

- **Thread oversubscription**: numpy/xarray operations often use multi-threaded BLAS under the hood by default (can silently try to grab all 64 cores per process). If several GCM processes run concurrently, each defaulting to all-core BLAS threading will massively oversubscribe the box. Explicitly cap per-process thread counts (`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, e.g. set to `cores / concurrent_processes`) before launching concurrent invocations.
- **Concurrency level**: use the benchmark runs' observed wall-clock and memory to decide concurrent process count — don't default to `cores - 2` on a shared box like this; start conservative (e.g. a handful of concurrent GCM processes, not dozens), watch `uptime`/load average once running, and scale up only if the machine stays responsive for other users. Check whether there's an informal big-job etiquette on this server (lab Slack, shared calendar, an admin to loop in) before committing a multi-day batch at higher concurrency.
- **Memory**: cap memory per process explicitly rather than letting each one expand to whatever's free — 339GB is currently available, but nothing stops another user's job from needing it mid-run, and nothing here would stop this job from taking it first. Budget well under the full amount in aggregate across however many concurrent processes you run.
- **Checkpointing (important for an unattended multi-day run, and for the incremental `--gcms` workflow itself)**: the (GCM, archetype, period) output file (Phase 5) is the checkpoint unit — write it as soon as that combination finishes, not held in memory alongside others for a single combined write. On startup, check for already-existing (GCM, archetype, period) output files and skip them — this way both a crash and a deliberate follow-up invocation with an overlapping or extended `--gcms` list resume/extend cleanly instead of recomputing finished work.
- **Logging**: per-file and per-GCM timing to a structured log file (Phase 2), not just stdout, so status across many concurrent processes is checkable without attaching to each one's terminal.

---

## Phase 4: I/O

- **Disk placement — confirmed**: `/home/redmond/data` and `/home/redmond/data2` are readable (the `pub/` climate data lives there) but not writable for this account, so neither is a valid output target. `<OUTPUT_DIR>` = the repo's own `outputs/` directory at `/home/server/pi/homes/queenl/projects/climate-risk-dc/outputs/`, on the `homes-new:/export/pi` mount — confirmed writable (the elevation file lives under the same tree), with ~56TB free pool-wide and current usage under this home directory at only ~21GB total. No per-user quota tool is installed to check a hard personal cap, but the ~100GB+ full-pipeline-output estimate (Phase 5) is nowhere near either number, so this isn't worth blocking on further. `/` and `/tmp` are still too small (~30GB and ~24GB respectively) regardless — don't fall back to either for anything beyond trivial scratch.
- **This is specifically a humidity problem, not a temperature problem**: tasmax/tasmin already comes pre-split to a Northwest region (Phase 1) — just clip that down to the Oregon boundary directly, no separate CONUS-scale subsetting pass needed. Humidity is full CONUS, so do a one-time up-front subset pass per GCM to an Oregon-plus-buffer extract (written to `<OUTPUT_DIR>`), and reuse that smaller extract across all archetype/period runs rather than re-reading the full CONUS humidity domain each time.
- If this pipeline will be re-run or iterated on, consider converting the Oregon-clipped extracts to Zarr for faster repeated chunked access over plain NetCDF — again, on `<OUTPUT_DIR>`, not the root or home filesystem.
- If anything in the pipeline ends up using dask (not required by the current plain-loop design, but possible for chunked reads of the humidity fields), point its local/spill directory explicitly at `<OUTPUT_DIR>` rather than the default, which would otherwise land on the small `/tmp`.

---

## Phase 5: Outputs

- **Primary output — one NetCDF per (GCM, archetype, period)**, e.g. `pue_wue_ACCESS-CM2_ssp370_ae-chiller_historical.nc`, `pue_wue_ACCESS-CM2_ssp370_ae-chiller_midcentury.nc`, `pue_wue_ACCESS-CM2_ssp370_ae-chiller_endcentury.nc` — three separate files per (GCM, archetype), not one file spanning all three windows. Each contains **daily** PUE and WUE for every Oregon grid cell across that window's real calendar dates. This is also the checkpoint/resume unit referenced in Phase 3.
- Use NetCDF compression (`zlib`, a reasonable `complevel`) and `float32` rather than `float64` — daily, full-grid output is a lot bigger than an annual-aggregate product would be. Rough order of magnitude: a few thousand cells × ~31,000 total days across the three windows × 2 variables × N_gcms (up to 27) × N archetypes lands somewhere around 100GB+ uncompressed (less if fewer than 27 GCMs have SSP3-7.0), same total regardless of the three-way file split. `<OUTPUT_DIR>`'s headroom is already confirmed comfortable for this (Phase 4) — worth a real size check once N archetypes is known, but not a blocking concern.
- Everything downstream of these files — ensemble median/spread across the 27 GCMs, monthly climatology, annual means, GeoTIFF/quicklook figures matching the existing 1.1 heat-processing output pattern — is a separate analysis step that reads them. It is not this script's job. Keep the compute script's scope exactly as stated in Phase 2: GCMs in, daily PUE/WUE NetCDFs out, nothing else.

---

## Phase 6: Validation

- Confirm the historical-period ensemble (whatever GCMs actually made it in, per Phase 1) brackets the previously validated single-GCM (ACCESS-CM2) facility-level result and the reported Oregon utility PUE figures (ref [12] in the methods doc). Note: MACA was never actually used anywhere in this project, despite the methods doc/`pue-explainer.md`/`four-week-timeline.md` describing a MACA pairing — those docs are stale on that point (worth a correction pass, separate from this task). That rules out a humidity-data-source mismatch as a confound here, but don't assume the facility-level benchmark used the exact same humidity product as this pipeline without checking — confirm what it actually used before treating this as a clean comparison.
- Re-run the Prineville/The Dalles spot-check against the full ensemble spread, not just one model — pull it directly from the (GCM, archetype, period) daily NetCDFs at those two grid cells, no separate validation-only output needed now that daily is the standard everywhere.
- Sanity-check output ranges against Lei & Masanet (2022) Figure 4/6 published ranges by climate zone, separately for each archetype.

---

## Limitations to carry into the paper (document regardless of how the run goes)

- SSP3-7.0 here vs. SSP5-8.5 in Section 1.1 — different forcing pathway for heat vs. cooling-efficiency risk within the same study.
- End-century window is 26 years (2075–2100), not the full 30-year window used for the other two periods.
- One ensemble member per model — this captures structural (inter-model) uncertainty, not a given model's internal variability.
- Elevation-adjusted pressure is new in this version; the facility-level pipeline it replaces did not have it — note the two are not directly comparable on that basis if both appear anywhere in the paper.
- [Archetype-specific limitations, once Phase 0 reports back what's actually implemented and what its calibration range is.]