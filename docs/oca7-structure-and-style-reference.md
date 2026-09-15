# OCA7 (Seventh Oregon Climate Assessment) — structure & style reference

Source: `OCA7_full.pdf`, 314 pp., published January 2025 by the Oregon Climate Change Research Institute (OCCRI) at Oregon State University. Editor: Erica Fleishman. Recommended citation format:

> Fleishman, E., editor. 2025. Seventh Oregon climate assessment. Oregon Climate Change Research Institute, Oregon State University, Corvallis, Oregon. https://doi.org/10.5399/osu/1181.

Individual chapters are cited the same way but with chapter authors and title as a chapter-in-report citation (standard practice, confirmed by the multi-author structure — see below). License: CC BY-NC-SA 4.0 (text/figures not individually credited to outside photographers/sources).

## Report structure (sector taxonomy)

OCA7 is organized into five parts, each containing one or more standalone chapters:

1. **State of Climate Science** — Trends in Climate and Advances in Climate Science; USDA Plant Hardiness Map changes; ENSO impacts
2. **Climate-Related Natural Hazards** — Projected Precipitation Changes; Freezing Rain/Ice Accretion (Willamette Basin); Drought History and Projections
3. **Adaptation Sectors**, subdivided into:
   - *Economy* — Wildfire/Timberland Value; Wildfire Smoke Economic Impacts; Business and Climate Change
   - *Natural Systems* — Afforestation/Reforestation Carbon Sequestration; Oregon Season Tracker
   - *Built Environment and Infrastructure* — Floating Offshore Wind Energy Infrastructure; Urban/Suburban Tree Planting Trade-offs
   - *Public Health* — Zoonotic Disease Transmission; Wildfire Smoke Health/Cost Scenarios; Drought and Health
   - *Social Systems* — Climate Litigation; Wildfire Local Solutions; Farmer Perspectives in Climate Modeling; Wine Industry Responses
4. Online Appendices (data/methods supplements for 3 of the chapters)

**There is no standalone "energy infrastructure" or "extreme heat" hazard chapter.** This is the key structural gap a data-center PUE/WUE/heat chapter would fill — it doesn't yet exist in OCA7's sector taxonomy. Given its subject matter (data center cooling infrastructure, electricity/water demand), the closest home is **Adaptation Sectors → Built Environment and Infrastructure**, alongside the offshore wind chapter. It could also straddle Economy (electricity demand/growth framing) depending on emphasis.

## Existing data-center content already in OCA7 (important anchor point)

The Floating Offshore Wind Energy Infrastructure chapter (p. 171) already frames data centers as a driver of Oregon electricity demand growth, in the context of motivating more renewable generation:

> "An emerging and substantive concern is that the demand for energy is increasing faster than previously projected. Five-year growth projections almost doubled, from 2.6 to 4.7 percent, between 2022 and 2023 (Wilson and Zimmerman 2023). New data centers (including cryptocurrency and artificial intelligence) and industrial facilities... are two of the main drivers of this sudden growth in energy load (Wilson and Zimmerman 2023)."

There's also a passing mention (p. ~150 area) of industry lobbying against emissions regulation for energy-intensive data centers (citing O'Donovan et al.). A new data-center chapter should explicitly reference/build on this existing framing rather than reintroduce the topic cold — it signals continuity with the rest of the assessment and gives a ready-made in-text citation (Wilson and Zimmerman 2023) for the electricity-demand-growth motivation.

**No mentions of PUE, WUE, cooling, or water use efficiency anywhere in OCA7** — this chapter would be introducing that vocabulary fresh to the assessment; plan to define PUE/WUE briefly for a non-specialist policy audience on first use, consistent with how other chapters define technical terms (e.g., "baseload," "dispatchable," "penetration" are defined inline in the wind chapter when first used).

## Extreme heat: how it's currently treated (fragmented, not centralized)

Heat is not a standalone hazard chapter; it's addressed piecemeal:
- A "**Extreme Heat Events**" subsection sits inside the *Projected Changes in Oregon Precipitation* chapter, opening with a detailed case study of the June 2021 Pacific Northwest heat dome (station records, AMS heat-dome definition, day-by-day account of record-breaking temperatures).
- Heat recurs in Public Health (smoke+heat compound effects), Social Systems (Multnomah County's 2021 heat-wave litigation), and the wine industry / agriculture chapter (heat-wave damage to vineyards).
- No standardized "heat index" or wet-bulb-temperature metric appears anywhere in the report.

This means a data-center chapter's extreme-heat-index analysis would be one of the more rigorous, systematic heat treatments in the assessment — worth flagging in the intro as complementary to (not overlapping with) the existing heat-dome case study, and possibly worth citing that case study as the reference historical event when contextualizing projected changes.

## Chapter template (structural pattern to follow)

Most chapters (freezing rain, drought, wildfire economics, zoonotic disease) follow a standard scientific-report IMRaD pattern with **no abstract or "key messages" box at the chapter level**:

```
[Chapter Title, centered]
[Author names, centered, no affiliations — affiliations only listed once on the report-wide Authors page]

Introduction
Methods  (sometimes "Approach")
Results
Discussion  (sometimes folded into Results, or a closing "Implications" framing)
[Figures/Tables inline, numbered sequentially within chapter]
References
```

The **Built Environment and Infrastructure** chapter (offshore wind) is the one exception: it's a multi-author compilation with named "Editors," organized as thematic subsections each with its own sub-byline (e.g., "Kaustubha Raghukumar" under Underwater Sound, "Sarah Henkel, Kyle Newton, and Taylor Chapple" under Electromagnetic Fields) rather than strict IMRaD. Given the data-center chapter is coming from a single, tightly-scoped research contract (one methodology, one model), the **standard IMRaD single-author(s) format is the better fit**, not the compilation format.

## Executive Summary style (for whatever gets rolled up there)

The Executive Summary is organized by the same sector taxonomy, and each item leads with a **bolded topical phrase** followed by 2-5 tightly quantified sentences, e.g.:

> "**Temperature and precipitation.** Oregon's annual average temperature increased by 2.2°F per century since 1895. Without considerable reductions in greenhouse gas emissions, annual temperature in Oregon is projected to increase by at least 5°F by 2074 and 7.6°F by 2100..."

> "**Built environment.** Offshore winds on the U.S. West Coast represent one of the most energetic and consistent renewable energy resources in the nation..."

Pattern to reuse for a data-center executive-summary paragraph: bolded lead-in (e.g., "**Data center cooling.**"), then 2-4 sentences with specific quantified projections (percent change in PUE/WUE, degree-days, number of hours above a threshold), no hedging beyond the calibrated-language conventions below.

## Climate methodology conventions already used in OCA7 (directly reusable for this project's methods section)

This project already uses LOCA2 downscaling (temperature, paired with the LOCA CMIP6 humidity companion dataset — not MACA, which was named in early drafts of this methodology but never actually used, see `docs/methods-document.md`'s Time Horizons section) — OCA7 uses closely related conventions, so the new chapter's methods section can match reporting style exactly. Scenario-wise, the data-center chapter itself now runs a 20-GCM SSP3-7.0 ensemble (confirmed 2026-09-14, see `docs/code-development-plan.md`); the project's original single-GCM baseline work used SSP5-8.5 and remains that way outside the chapter:
- **Downscaling**: LOCA (Localized Constructed Analogs; Pierce et al. 2014, 2023) is OCA7's own precipitation/temperature downscaling method — directly citable as precedent.
- **Scenario set**: OCA7 typically reports **three SSPs together** — SSP2-4.5, SSP3-7.0, SSP5-8.5 — rather than SSP5-8.5 alone. Table captions specify ensemble size per SSP (e.g., "18 ensemble members from 8 global climate models (SSP2-4.5), 52 ensemble members from 11 global climate models (SSP3-7.0), and 18 ensemble members from 7 global climate models (SSP5-8.5)"). Worth deciding explicitly whether the data-center chapter stays SSP5-8.5-only (as the current methods doc does) or adds a lower-scenario comparison to match OCA7 convention — flagging this as a design decision rather than assuming.
- **Baseline period**: 1950–2014 historical baseline is OCA7's standard reference period for percent-change framing.
- **Uncertainty reporting**: ensemble mean plus **5th–95th percentile range across models**, reported in parentheses — this is the exact format to mirror for PUE/WUE and heat-index projection tables.
- **Seasonal breakdown**: winter (DJF) / spring (MAM) / summer (JJA) / autumn (SON), consistently used for every seasonal table.

## Citation and reference style

In-text: author-date, e.g. "(Rupp et al. 2017, Zhang et al. 2023)" or narrative "Pierce et al. (2014)... ". Multiple citations comma-separated, chronological or as-relevant order, no "et al." period-comma inconsistency to worry about — matches standard ecology/climate-journal convention.

Reference list: alphabetical by first author surname, journal names spelled out in full (not abbreviated), volume:pages, DOI given as a full https://doi.org/... URL. Example:

> Pierce, D.W., D.R. Cayan, and B.L. Thrasher. 2014. Statistical downscaling using localized constructed analogs (LOCA). Journal of Hydrometeorology 15:2558–2585.

## Voice and register

- Calibrated/IPCC-style confidence language throughout: "likely," "very likely," "projected to," "may," "is uncertain whether" — never bare unqualified claims about the future.
- Third person, no first person ("we"/"our") in the scientific chapters (the Introduction and Acknowledgments do use "we" as OCCRI's institutional voice, but chapter bodies don't).
- Definitions given inline at first use for any technical term a non-specialist policy reader wouldn't know (modeled directly on how "baseload," "dispatchable," "penetration" are defined in the wind energy chapter — do the same for PUE, WUE, wet-bulb temperature, etc.).
- Numbers are precise and consistently unit-doubled (metric first or imperial first depending on chapter, always with the other in parentheses — e.g. "100 m (328 ft.)").
- Chapters close with adaptation-relevant framing even without a formal "Implications" header — the takeaway for decision-makers is woven into the Discussion rather than bullet-pointed.

## Open decision for the new chapter

Whether to write it as a **standalone Built Environment and Infrastructure sub-chapter** (IMRaD, single-topic, matching freezing-rain/drought style) or propose it as a **new addition to the Economy section** (framed around electricity/water demand and cost, extending the "Business and Climate Change" and offshore-wind chapters' demand-growth narrative). IMRaD-under-Built-Environment looks like the stronger structural fit given the project's methods (physics-based PUE/WUE model + climate projections is squarely a "built environment / infrastructure" analysis), but this is worth confirming with Laura/Erica rather than assuming.