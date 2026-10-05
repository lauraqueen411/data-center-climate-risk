# Chapter Outline: Climate Risk to Data Center Cooling in Oregon — PUE, WUE, and Extreme Heat

Draft outline for the OCA7-style chapter on this project's ultimate goal. Scoped to heat only (flood dropped per `four-week-timeline.md`); built around the site-based, full-LOCA2-ensemble, 3-archetype PUE/WUE run that replaced the shelved grid-based approach (`pue-wue-gridded-ensemble-task.md`). Structured against the conventions documented in `oca7-structure-and-style-reference.md` — single-author(s) IMRaD, no key-messages box, calibrated language, ensemble mean + 5th–95th percentile range, DJF/MAM/JJA/SON seasonal convention.

**Confirmed run parameters:** 109 site-based facilities, 20-model LOCA2 CMIP6 ensemble (of 27 candidates — the ones with SSP3-7.0 downscaled for both temperature and its humidity companion; see `outputs/ensemble_pue_wue/gcm_manifest.csv`), SSP3-7.0, three time windows (historical / mid-century / end-century), three cooling archetypes per Lei & Masanet (2022) Table 2:

| Case | DC size class | Cooling system configuration | Direct evaporation | Space humidification | Cooling tower |
|---|---|---|---|---|---|
| 1 | Large-scale | Airside economizer + adiabatic cooling, **with water-cooled chiller as supplemental** (used only when the economizer + adiabatic system can't meet the cooling load) | ✓ | ✓ (adiabatic) | ✓ |
| 2 | Large-scale | Waterside economizer + water-cooled chiller | n/a | n/a | ✓ |
| 5 | Midsize | Water-cooled chiller only — **no economizer at all** | n/a | n/a | ✓ |

**Correction (2026-09-14):** this box previously named the second archetype "Case 3" (midsize, airside economizer, no adiabatic) — confirmed with the user to be a typo. **Case 2** (large-scale, waterside economizer) is what was actually implemented and run; it appears nowhere as "Case 3" in the codebase. Case 2 is not a stripped-down Case 1 (it uses a different economizer *type* — waterside, not airside — not "Case 1 minus adiabatic cooling"), so the comparison across the three archetypes is really: two different large-scale economizer technologies (Case 1, Case 2) against no economizer at all (Case 5), not a single graded axis of "how much free cooling is layered on."

Two things worth carrying into Results and Discussion, both stated directly in the paper and **now confirmed against the actual Oregon ensemble output** (`outputs/ensemble_pue_wue_facilities/tables/table4_summary_statistics.csv`): (1) reported PUE medians are markedly higher for the midsize case than the large-scale ones (large-scale Cases 1-2: median 1.12–1.25; midsize cases including 5: median 1.39–1.98) — **confirmed**: historical fleet-mean PUE is 1.134 (Case 1), 1.136 (Case 2, essentially tied with Case 1 as the paper's own size-class grouping predicts), and 1.648 (Case 5); (2) counterintuitively, the paper found Case 1's *water* use is observably **lower** than Case 2's, even though adiabatic cooling adds direct evaporative water use — because it also cuts how often the more water-intensive chiller+cooling-tower system has to run. **Confirmed, and dramatically so**: historical fleet-mean WUE is 0.25 L/kWh (Case 1) vs. 2.08 L/kWh (Case 2), roughly an 8x gap — much larger than "counterintuitively lower" might suggest on its own.

**SSP note (resolved 2026-09-14):** Section 1.1 (acute heat/exceedance-day risk) has since been extended to the same 20-GCM SSP3-7.0 ensemble used here (`scripts/run_heat_indices_ensemble.py`), replacing its original single-GCM/SSP5-8.5 figures for this chapter. Figure 7 pairs both metrics from the same ensemble/scenario — no mismatch caveat needed there. The original single-GCM/SSP5-8.5 heat output remains in the repo for uses outside this chapter.

---

## Proposed placement in OCA7's taxonomy

Built Environment and Infrastructure, alongside the offshore wind chapter — see the structure reference doc for why (it's the closest existing analog, and it's already the chapter that name-checks data centers as an Oregon electricity-demand driver, giving this chapter a ready in-text citation to open with).

## Working title

"Climate Risk to Data Center Cooling Efficiency in Oregon: Projected Changes in Power and Water Usage Effectiveness Under Extreme Heat" (placeholder — tighten once scope is final)

---

## 1. Introduction

- Data centers as physical infrastructure (servers, networking, power, cooling); Oregon's early role (Google's The Dalles facility, 2006, >1.3M sq ft) — reuse methods-document.md's opening framing.
- National context: electricity demand from data centers projected 4.4%→~12% of U.S. consumption by 2028; direct water consumption for cooling roughly doubling to quadrupling over the same window (Shehabi et al. 2024).
- **Direct link into OCA7's existing narrative**: the offshore wind chapter (p.171) already identifies data centers (AI + cryptocurrency) as a driver of Oregon's accelerating electricity-demand growth (5-yr growth projections nearly doubled, 2.6%→4.7%, 2022–2023; Wilson and Zimmerman 2023). Open by building on that sentence rather than reintroducing data centers as a new topic — this is the seam that ties the new chapter into the existing report.
- IPCC hazard/exposure/vulnerability framing (hazard = heat; exposure = 109 Oregon facilities; vulnerability = cooling-system dependence on ambient temperature and humidity).
- Define PUE and WUE briefly, in plain language, on first use (matches OCA7's convention of defining technical terms inline — see how "baseload"/"dispatchable" are handled in the wind chapter).
- Literature gap: Kollar & Grady (2025) and Esparza et al. (2025) are both historical-hazard, facility-location-based, national-scale studies; neither uses future climate projections or models cooling-system response. This chapter is projection-based and Oregon-specific, extending Lei & Masanet's (2020, 2022) validated PUE/WUE model with future climate forcing instead of historical weather (see `climate-extension-justification.md` for the "morphing" precedent).
- **Fills an actual gap in OCA7 itself**: no existing chapter treats data center energy/water infrastructure, and extreme heat is currently scattered (June 2021 heat-dome case study buried in the precipitation chapter, plus mentions in health/legal/agriculture chapters) rather than centralized — this chapter is one of the more systematic heat treatments in the assessment.
- Close with an explicit scope statement: heat only (flood named as future work, consistent with the rest of the paper); site-based rather than full-grid, because full-grid was computationally infeasible (worth one sentence on why, for transparency); SSP3-7.0; 3 archetypes.

## 2. Data and Methods

- **Site inventory**: 109 Oregon data center facilities (IM3/PNNL Data Center Inventory) — table or map (Figure 1) showing count, size distribution if available, and geographic spread (Willamette Valley cluster vs. Columbia Gorge / Prineville vs. elsewhere).
- **Climate inputs**: LOCA2 CMIP6 ensemble (tasmax/tasmin) paired with the LOCA CMIP6 humidity companion dataset (Pierce & Cayan, Scripps/UCSD) — **confirmed (2026-09-14), not MACA**. MACA was named in early drafts of this methodology but never actually used anywhere in this project; see `docs/methods-document.md`'s Time Horizons section.
- **Scenario and time windows**: SSP3-7.0; historical baseline (1985-2014), mid-century (2045-2074), end-century (2075-2100, 26 years not 30) — **confirmed** to match the shelved grid task's windows exactly.
- **Ensemble**: 20 of 27 LOCA2 GCMs have SSP3-7.0 downscaled for both temperature and humidity and are included, one realization per model — full list and member IDs in `outputs/ensemble_pue_wue/gcm_manifest.csv`.
- **Extreme heat index**: ASHRAE Class A2 threshold exceedance-day counts (Section 1.1's method, now run across the same 20-GCM SSP3-7.0 ensemble as PUE/WUE rather than the original single-GCM/SSP5-8.5 run) — daily high/low tracked together since cooling systems need cool nights to recover during multi-day events.
- **PUE/WUE model**: Lei & Masanet (2020, 2022) hybrid physics/statistical model; three archetypes — Case 1 (large-scale, airside economizer + adiabatic cooling + supplemental water-cooled chiller), **Case 2** (large-scale, waterside economizer + water-cooled chiller — corrected from "Case 3," a typo, 2026-09-14), Case 5 (midsize, water-cooled chiller only, no economizer) — see Table 2 in the box above for the full definitions.
- **Validation**: historical-period ensemble bracketing the previously validated single-GCM result and Oregon utility PUE figures (Guidi et al. 2026); Prineville and The Dalles spot-check against the full ensemble spread — **done**: ensemble mean matches the single-GCM (ACCESS-CM2) point estimate almost exactly at both sites (`outputs/ensemble_pue_wue_facilities/tables/table3_historical_validation.csv`).
- **Uncertainty framing**: ensemble mean plus 5th–95th percentile range across GCMs — matches OCA7's own convention (see Table 2's caption style in the precipitation chapter) and is used consistently in every results figure/table here.

## 3. Results

**Superseded 2026-09-15**: the numbered plan and the "Candidate figures and tables" catalog below reflect the original nine-figure, raw-value-led design. That design has been replaced by a four-figure, delta-first structure, arrived at by requiring each figure to state one verifiable claim rather than hold as many dimensions as possible. `docs/results-draft.md` is now the authoritative Results outline; its own header explains what was tried and rejected along the way (a percent-change chart needing three separate axis reads, and a facility-level correlation that turned out to be a correlation-of-ratios artifact). The list below is kept for historical reference only.

Proposed order, each tied to a figure/table (full list below):

1. Historical (baseline) fleet characterization — statewide mean PUE and WUE by archetype, with ensemble spread. Establishes where Oregon's fleet sits today (recall the earlier single-GCM/Case-1-only result: fleet PUE ~1.12, low end of industry range — now check whether that holds up across the full ensemble and the other two archetypes).
2. Projected change — mid-century and end-century, by archetype, ensemble mean ± range. This is the chapter's central result.
3. Extreme heat index trends (exceedance-days) alongside PUE trends — connects the acute-risk metric (Section 1.1) to the chronic-risk metric (Section 1.2) directly, which is exactly the "compound risk" framing the methods doc already sets up.
4. Archetype comparison — Case 1 and Case 2 are both large-scale designs with an economizer (airside+adiabatic vs. waterside) while Case 5 has none; does that split (rather than a single graded "amount of free cooling" axis, see the corrected Table 2 above) change climate sensitivity, on top of the baseline absolute-PUE gap the source paper already documents? **Confirmed**: Case 1 shows real climate sensitivity (WUE nearly doubles from historical to end-of-century, seasonal peak shifts July→August) while Case 2 and Case 5 are nearly flat across periods — Case 1's WUE < Case 2's WUE ordering (a specific published finding) holds strongly under Oregon's climate (0.25 vs. 2.08 L/kWh historically), and persists under warming.
5. Site-level / spatial variation across the 109 facilities — does location within Oregon matter (Willamette Valley vs. high desert vs. Gorge)?
6. **Climate sub-region breakdown (resolved 2026-09-15)** — EPA Level III ecoregions (Omernik 1987), matching OCA7's own convention of presenting results averaged across U.S. Level III ecoregions in Oregon (Fig. 1). Rather than averaging temperature/humidity across each ecoregion's facilities before running the PUE/WUE model — which would run into Jensen's-inequality problems given the model's thresholded economizer/chiller switching — one representative facility per ecoregion was selected (nearest to that ecoregion's own facility-cluster centroid) and its already-computed daily PUE/WUE pulled directly from the 109-site run, following Lei & Masanet (2022)'s own convention of representing a climate zone with one representative point (their per-IECC/ASHRAE-zone representative city, run on a full year of TMY data) rather than an average. Of Oregon's 9 Level III ecoregions, 5 contain at least one facility: Columbia Plateau (61), Willamette Valley (33), Blue Mountains (13), Eastern Cascades Slopes and Foothills (1), Klamath Mountains/California High North Coast Range (1).

## 4. Discussion / Implications for Adaptation

- What rising PUE/WUE means for Oregon's electricity and water systems, explicitly tied back to the offshore wind chapter's demand-growth narrative — this chapter quantifies one piece of *why* that demand keeps rising (cooling gets less efficient, not just more numerous facilities).
- Siting and cooling-technology choice as adaptation levers, if the archetype comparison shows a meaningful spread — a genuinely actionable takeaway for the "Implications" framing OCA7 chapters close with.
- Positions this chapter as the infrastructure-dependency piece Kollar & Grady and Esparza et al. both flagged as missing from the literature, but note plainly that this chapter, per the four-week scope cut, does the *lightweight* version of infrastructure dependency (Section 1.3's coincidence-flag, if included) rather than full derating math — don't overclaim it closes that gap.

## 5. Limitations

- SSP3-7.0 (this analysis) vs. SSP5-8.5 — **resolved (2026-09-14)**: Section 1.1 heat risk now also runs the same 20-GCM SSP3-7.0 ensemble for this chapter, so PUE/WUE and extreme-heat figures share one forcing pathway. SSP5-8.5 remains the scenario for the original single-GCM heat output used outside this chapter (e.g. Section 1.3).
- Site-based (109 discrete points), not full statewide grid — grid-based was attempted and found computationally infeasible; state this plainly rather than silently narrowing scope.
- Generic Table A.1 midpoint equipment parameters for every facility, not facility-specific specs (proprietary, unavailable).
- GP-model extrapolation beyond training range for the chiller-COP component all three cases share (water-cooled chiller) — carry over the same due-diligence framing already written in `climate-extension-justification.md`; this applies equally to all three cases here since they all fall back on the same chiller type, unlike the DX/air-cooled cases considered (and shelved) in the grid-based task doc.
- WUE has a weaker validation basis than PUE — no facility-level ground truth for Prineville or The Dalles; validated only against Lei & Masanet's (2022) published climate-zone ranges.
- One ensemble member per model — captures inter-model (structural) uncertainty, not any single model's internal variability.
- End-century window may be shorter than the other two (confirm; was 26 vs. 30 years in the shelved grid plan) — flag if it carries over.
- MACA-vs-actual-humidity-source correction (see Methods note above) — **resolved (2026-09-14)**: confirmed the LOCA CMIP6 humidity companion was used throughout, MACA was never actually used anywhere in this project. `docs/methods-document.md` was already correct on this point before this correction pass.

## 6. References

Carry over methods-document.md's reference list wholesale (Lei & Masanet 2020/2022, Shehabi et al. 2024, Kollar & Grady 2025, Esparza et al. 2025, IPCC 2022, etc.) plus Wilson and Zimmerman (2023) and Guidi et al. (2026), already cited elsewhere in this project. Add whatever the site-based ensemble pipeline's own documentation cites for GCM provenance (LOCA2 core papers: Pierce et al. 2014, 2023).

---

## Candidate figures and tables

**Superseded 2026-09-15** — see the note under Section 3 above. This nine-figure/five-table catalog is kept for historical reference; the current figure set (four figures, built by `scripts/build_results_delta_figures.py`) is documented in `docs/results-draft.md`.

Ordered roughly as they'd appear in Results. Each entry notes the OCA7 convention it's modeled on and what the CLI analysis task needs to produce it.

**Figure 1 — Site map.** All 109 Oregon data center locations, sized/colored by historical baseline PUE (or by IT capacity if available). Modeled directly on OCA7's Figure 1 (offshore wind resource map) — a context-setting map before any results. Good opener. Built: `outputs/ensemble_pue_wue_facilities/figures/figure1_site_map.png`.

**Figure 2 — Fleet-mean PUE by period and archetype.** Grouped bars or box plots: historical / mid-century / end-century × Case 1 / 2 / 5, each showing ensemble mean and 5th–95th percentile whiskers. Mirrors OCA7's Table 2 percentile convention rendered as a figure. This is likely the single most important figure in the chapter. Built: `outputs/ensemble_pue_wue_facilities/figures/figure2_pue_by_period_archetype.png`.

**Figure 3 — Fleet-mean WUE by period and archetype.** Same structure as Figure 2, for WUE. Since all three archetypes here use a water-cooled chiller + cooling tower, none should show near-zero WUE — the interesting test is whether Case 1 < Case 2 (the source paper's counterintuitive finding: adiabatic cooling adds direct water use but cuts chiller/tower runtime enough to net lower overall) holds for Oregon's climate and strengthens or weakens under future warming. **Confirmed**, and by a wider margin than the paper's framing might suggest (0.25 vs. 2.08 L/kWh historically, ~8x). Built: `outputs/ensemble_pue_wue_facilities/figures/figure3_wue_by_period_archetype.png`.

**Figure 4 — PUE-vs-WUE relationship.** Scatter across all sites × archetypes × periods, colored by period, shaped by archetype — answers the "are PUE and WUE correlated, and does that depend on cooling technology" question the four-week timeline already flagged as a planned figure, now enriched with the ensemble/archetype dimensions that weren't available when that plan was written. Built: `outputs/ensemble_pue_wue_facilities/figures/figure4_pue_vs_wue_scatter.png`.

**Figure 5 — Monthly climatology, historical vs. mid- vs. end-century.** Small multiples (one panel per archetype) showing the "taller and later peak season" signal already observed in the Case-1-only single-GCM run (peak shifts later, not just intensifies) — check whether that pattern holds across the full ensemble and other archetypes, or was partly a single-model artifact. **Checked**: holds, but the shift is smaller than the single-GCM run suggested — July→August (one month), not July→September, across the ensemble for Case 1 and Case 5; Case 2 stays essentially flat year-round (consistent with its much weaker climate sensitivity overall). Built: `outputs/ensemble_pue_wue_facilities/figures/figure5_monthly_climatology.png`.

**Figure 6 — Site-level map of projected change (ΔPUE, end-century minus historical).** Choropleth/point map across the 109 facilities, colored by magnitude of change — the natural "does location within Oregon matter" figure, and a direct visual complement to Figure 1. Faceted by all three archetypes; Case 1 shows the largest and most spatially varied delta, Case 2 the smallest and flattest. Built: `outputs/ensemble_pue_wue_facilities/figures/figure6_delta_map.png`.

**Figure 7 — Extreme heat index and PUE trend together.** Two-panel figure: exceedance-day counts (Section 1.1's metric, now on the same 20-GCM SSP3-7.0 ensemble as PUE — no scenario mismatch to caption) alongside fleet-mean PUE (Case 1), historical through end-century. Exceedance days roughly quadruple (13→37→53 days/yr) while PUE degrades only modestly (1.13→1.15→1.16) — a legitimate, reportable finding that raw heat exposure and cooling-system efficiency respond very differently to the same warming. Built: `outputs/ensemble_pue_wue_facilities/figures/figure7_heat_pue_trend.png`.

**Figure 8 — Ensemble spread / structural uncertainty.** Spaghetti plot or fan chart of statewide-mean PUE trajectory across individual GCMs, historical through end-century — makes the inter-model spread visible in a way Figure 2's whiskers only summarize; useful if reviewers will ask "how much does the choice of GCM matter here." Built: inter-model spread turns out to be small relative to the multi-decadal warming trend for all three archetypes. `outputs/ensemble_pue_wue_facilities/figures/figure8_ensemble_spread.png`.

**Figure 9 — PUE by EPA Level III ecoregion.** Resolved (2026-09-15): grouped bars by ecoregion x period, Case 1 only, one representative facility per ecoregion (see Results item 6 above) — ensemble mean + 5th-95th pct. whiskers, same convention as Figure 2. Built: `outputs/ensemble_pue_wue_facilities/figures/figure9_ecoregion_representative_sites.png`.

**Table 1 — Data sources.** Adapt directly from methods-document.md's existing table, updated for the site-based/ensemble/multi-archetype specifics (GCM count actually used, humidity source once confirmed, archetype definitions once confirmed). Built: `outputs/ensemble_pue_wue_facilities/tables/table1_data_sources.csv`.

**Table 2 — Archetype definitions.** Reproduce the Case 1/2/5 definitions confirmed above (adapted from Lei & Masanet 2022 Table 2) — all three use the same Table B.1 midpoint equipment-parameter convention already established for the single-GCM 1.2 pipeline. Built: `outputs/ensemble_pue_wue_facilities/tables/table2_archetype_definitions.csv`.

**Table 3 — Historical validation.** Modeled vs. published/utility PUE at Prineville and The Dalles, ensemble spread vs. single-GCM point estimate, for Case 1 (the archetype both real facilities' documented cooling systems and the paper's own validation set actually match — Case 2/5 weren't checked against these two specific facilities, since they're a different economizer type/size class). Built: `outputs/ensemble_pue_wue_facilities/tables/table3_historical_validation.csv` — ensemble mean matches the single-GCM point estimate almost exactly at both sites.

**Table 4 — Summary statistics.** Ensemble mean ± 5th–95th percentile range, PUE and WUE, by period × archetype, statewide — the numeric backbone the Executive-Summary paragraph (bolded lead-in, OCA7 style) would be drawn from. Built: `outputs/ensemble_pue_wue_facilities/tables/table4_summary_statistics.csv`.

**Table 5 — Ecoregion representative sites.** Ensemble mean ± 5th–95th percentile range, PUE and WUE, by EPA Level III ecoregion × period × archetype, one representative facility per ecoregion (see Results item 6). Built: `outputs/ensemble_pue_wue_facilities/tables/table5_ecoregion_representative_sites.csv`.