# Task for Claude Code: Analyze Site-Based PUE/WUE Ensemble Output for the OCA7 Chapter

## Context (for you, before starting)

This is the analysis step downstream of the site-based (109 facilities) LOCA2-ensemble PUE/WUE run — the grid-based version of this (`pue-wue-gridded-ensemble-task.md`) was shelved as computationally infeasible; this run stayed site-based instead. The run already exists: your job is to analyze its output and produce the figures and tables listed in `data-center-chapter-outline.md`, not to rerun the PUE/WUE computation itself.

Confirmed by Laura (Sept 11): 109 Oregon data center facilities, full LOCA2 CMIP6 ensemble (however many GCMs actually have the run scenario downscaled), SSP3-7.0, three time windows (historical / mid-century / end-century), three cooling archetypes — Lei & Masanet (2022) Case 1, Case 3, and Case 5. Their exact definitions (confirmed from the source paper, Table 2):

| Case | Size class | Configuration | Direct evaporation | Space humidification | Cooling tower |
|---|---|---|---|---|---|
| 1 | Large-scale | Airside economizer + adiabatic cooling, water-cooled chiller supplemental | ✓ | adiabatic | ✓ |
| 3 | Midsize | Airside economizer only (no adiabatic), water-cooled chiller supplemental | n/a | adiabatic | ✓ |
| 5 | Midsize | Water-cooled chiller only, no economizer | n/a | adiabatic | ✓ |

Two published, falsifiable expectations from that same paper, worth explicitly checking against this Oregon output rather than assuming they transfer: (1) absolute PUE is markedly higher for Case 5 than Case 1 (paper's reported medians: large-scale cases ~1.12–1.25, midsize cases ~1.39–1.98); (2) Case 1's WUE is *lower* than Case 3's despite adiabatic cooling adding direct water use, because it reduces how often the more water-intensive chiller+tower system runs.

Run Phase 0 first and report back before building any figures — it determines file paths, schema, and whether the assumptions above (SSP, time windows, archetype-to-parameter mapping) actually match what was run, which prior project docs have gotten wrong before (MACA was named as the humidity source in three docs but was never actually used, per `pue-wue-gridded-ensemble-task.md`'s validation note — don't repeat that kind of drift here).

---

## Phase 0: Inventory the actual output (do this first, report back)

- Locate the site-based ensemble run's actual output files (ask Laura for the path if it's not obvious from the repo/outputs directory) and report: file format (NetCDF/CSV/Parquet/other), one row/record's actual dimensions (site × GCM × archetype × period × variable × time-resolution — is it daily or already aggregated to annual/monthly?), and the exact list of the 109 facility identifiers and their lat/lon.
- Confirm the exact list of GCMs actually included in the ensemble (don't assume all LOCA2 models — cross-reference against which ones have SSP3-7.0 downscaled, same caveat the shelved grid task flagged) and how many realizations per model (one per model, per that task doc's convention, but confirm here rather than re-assume).
- Confirm the exact year ranges used for historical / mid-century / end-century in this run — the shelved grid task used 1985–2014 / 2045–2074 / 2075–2100 (with the end-century window 4 years short), but don't assume the site-based run matches without checking the actual run config/metadata.
- Confirm which humidity data source actually fed this run. Check directly rather than trusting `methods-document.md`, `pue-explainer.md`, or `four-week-timeline.md`, all of which name MACA — already flagged elsewhere in this project as likely stale.
- Confirm the equipment parameters actually used for Case 1, 3, and 5 in this pipeline (pull from whatever config/parameter file the run used) against Lei & Masanet (2022) Table A.1/B.1 published values — report any deviation from the paper's published midpoints, since the existing single-archetype pipeline (`pue-explainer.md`) used "all 26 equipment parameters at the midpoint of Table A.1's published bounds," and that convention should either carry over for Cases 3/5 too or be flagged if it doesn't.
- Confirm whether Section 1.1's extreme-heat exceedance-day output (for Figure 7 in the outline) already exists per-facility for this same set of 109 sites, same time windows — if it's SSP5-8.5 rather than this run's SSP3-7.0, that scenario mismatch needs to be visible in Figure 7's caption, not silently glossed over.
- Report back: file locations, schema, confirmed GCM list and count, confirmed time windows, confirmed humidity source, confirmed archetype parameters, and whether 1.1's heat output is available and on what scenario. Don't proceed to Phase 1 until this is reviewed.

---

## Phase 1: Ensemble summary statistics

Build a reusable summary module (not one-off per figure) that, for any (metric ∈ {PUE, WUE}, archetype, period) combination, computes across the confirmed GCM ensemble:
- mean, median
- 5th and 95th percentile (matching OCA7's own reporting convention — see `oca7-structure-and-style-reference.md`)
- both at statewide-fleet level (mean/percentile across all 109 sites × all GCMs) and at individual-site level (mean/percentile across GCMs only, per site)

Also compute, per site and per GCM: Δ(metric) = future period value − historical value, for both mid-century and end-century, for each of PUE and WUE, for each archetype. This is the core quantity most of the outline's figures need.

---

## Phase 2: Figures

Produce each of the following (see `data-center-chapter-outline.md` for full rationale per figure — this phase is the concrete build spec):

1. **Site map** — 109 facility locations, colored/sized by historical baseline fleet-mean PUE (Case 1, or all three archetypes as small multiples if that reads better).
2. **Fleet-mean PUE by period × archetype** — grouped bars or box plots, historical/mid/end-century × Case 1/3/5, ensemble mean + 5th–95th percentile whiskers.
3. **Fleet-mean WUE by period × archetype** — same structure as #2. Explicitly annotate or caption whether Case 1 < Case 3 (the paper's finding) holds here, for each period.
4. **PUE-vs-WUE scatter** — all sites × archetypes × periods, colored by period, shaped by archetype.
5. **Monthly climatology small multiples** — one panel per archetype, historical/mid/end-century overlaid, for PUE (and WUE if it shows a clear seasonal signal). Check whether the previously observed "taller and later peak" pattern (July→September peak shift, Case 1 single-GCM result) holds across the full ensemble and across Case 3/5.
6. **Site-level ΔPUE map** (end-century minus historical) — point map across the 109 facilities, colored by magnitude, one map per archetype or a faceted set.
7. **Extreme heat index + PUE trend** — two-panel or dual-axis, exceedance-day counts alongside fleet-mean PUE, historical→end-century. Only build this once Phase 0 confirms 1.1's heat output is available for the same sites/windows; caption the SSP mismatch explicitly if it's SSP5-8.5 vs. this run's SSP3-7.0.
8. **Ensemble spread / spaghetti plot** — statewide-mean PUE trajectory per individual GCM, historical→end-century, one line per model, for Case 1 (or faceted by archetype if it's not too busy).
9. **(Only if Laura confirms the sub-region analysis is happening)** ΔPUE by Oregon climate sub-region — hold this until she decides; don't build speculatively.

Save each as both a standalone PNG/vector file and, where practical, the underlying summary table (CSV) it was built from — reviewers and co-authors will want to check numbers without re-deriving them from raw ensemble output.

---

## Phase 3: Tables

1. **Data sources** — adapt `methods-document.md`'s table, filled in with Phase 0's confirmed specifics (GCM count/list, humidity source, archetype parameter source).
2. **Archetype definitions** — Case 1/3/5, reproduced from the box above plus Phase 0's confirmed equipment parameters.
3. **Historical validation** — modeled vs. published/utility PUE at Prineville and The Dalles specifically (not fleet-wide), ensemble spread vs. the previously validated single-GCM (ACCESS-CM2) point estimate, by archetype where the facility's actual archetype is known/assumed.
4. **Summary statistics** — ensemble mean ± 5th–95th percentile range, PUE and WUE, by period × archetype, statewide. This is the number set the chapter's executive-summary paragraph would quote directly.

---

## Phase 4: Validation / sanity checks

- Re-run the Prineville/The Dalles spot-check against the full ensemble spread (not just ACCESS-CM2), for whichever archetype most plausibly matches each facility's real cooling system (Case 1/AE is the one the original single-GCM pipeline validated against both facilities — confirm Case 3/5 aren't expected to match these two real facilities as closely, since they're a different size class/technology).
- Confirm the Case 5 > Case 3 > Case 1 absolute-PUE ordering from the source paper holds in the Oregon output; report clearly if it doesn't; don't force the comparison to match by adjusting parameters after the fact.
- Confirm or refute the Case 1 WUE < Case 3 WUE ordering; either result is reportable, but state which one was found.
- Sanity-check output ranges against Lei & Masanet (2022) Figure 4/6 published ranges by climate zone, separately for each of Case 1/3/5.
- Flag anything in the actual output that contradicts an assumption baked into `data-center-chapter-outline.md` (e.g., if the humidity source, time windows, or GCM count turn out different from what's assumed there) so the outline gets corrected before the chapter draft is written from it.

---

## Deliverable

A short written report (not just files) covering: Phase 0's findings, the built summary module and where it lives, all figures/tables with file paths, and Phase 4's validation results — specifically calling out anywhere this analysis diverged from what `data-center-chapter-outline.md` assumed, so that doc can be corrected before it's used to draft chapter prose.