# Methodology audit: findings

Audit of the heat chapter's code, data and methods text (`docs/results-figures/methodology.md`) against the source papers and the authors' reference code, run 2026-10-08 per `docs/tasks/methodology-audit.md`. One row per finding is in `findings.csv`. The server checks' raw outputs are in `outputs/audit_part_a/`.

## Decisions needed (server access ends 2026-10-11)

| Finding | What it affects | Options |
|---|---|---|
| **A7b: Case 5 change is roughly half chiller-COP extrapolation** | All Case 5 numbers: Figure 2(d, g), the Case 5 lines in Figure 4 (including the nonzero winter change), and Table 1 Cases 2/5 | 1. Keep it and state it as an inherited limitation of the Lei and Masanet model. 2. Clamp the COP inputs to the training range and rerun Case 5 (~40 min); Cases 1 and 2 change by < 0.00003 |
| **B2b: Case 2 cooling-tower pump bug** | Case 2 PUE level (up to −0.004 on hot days) and change (about −1.4%) | 1. Document it. 2. Fix the one line and rerun Case 2 (~40 min) |
| **C3a: atlas records counted as facilities** | The n values in captions and Table 1, and the word "facility" | 1. Define a facility as an atlas record in the text. 2. Merge campus and building records at the same site and rebuild locally (no server needed) |

If both server options are chosen, Cases 2 and 5 can be rerun together in about an hour.

## Findings by severity

**High**
- **A7b.** Below about 10.8 °C wet-bulb, the chiller COP model is outside its training data, and it keeps raising COP linearly (about 10 at −10 °C wet-bulb). Most cold-weather pairings fall there. Case 1 and Case 2 don't run the chiller in cold weather, so they are unaffected. Case 5 runs its chiller all year, so warmer winters lower its extrapolated COP. Clamping the inputs to the training range cuts Case 5's mid-century PUE change by 46–69% (CP 0.0084 → 0.0046, WV 0.0090 → 0.0049; ACCESS-CM2). The source papers use the same unbounded model, so this is an inherited limitation, not a porting error.

**Medium**
- **C3a.** The 109 "facilities" are atlas records (100 buildings, 6 campuses, 3 points). The Meta and Apple Prineville campuses are listed alongside their own buildings, and the Google campuses at The Dalles each have a building listed beside them. Values barely change, because duplicates share grid cells, but n overstates the number of distinct sites, especially in BM.
- **B5.** The papers run the model on hourly typical-year weather. This chapter uses two daily pairings, because LOCA2 is daily. `methodology.md` describes the pairing but doesn't say it approximates the papers' hourly method.
- **A7a, D1a–c.** `methodology.md` doesn't say that:
  - equipment parameters are held at the Table B.1 midpoints (three were corrected on 2026-10-08)
  - chiller COP comes from a model fit to manufacturer data, extrapolated outside its training range
  - the uncertainty bands cover climate models only, not equipment parameters

**Low**
- **B2b.** Case 2 cooling-tower pump bug (see above).
- **B2c.** Inherited from the reference code: Case 2's cooling-tower water use is computed from the IT heat load, not the heat actually rejected.
- **B2d.** The supply-air humidity is set to the midpoint of its range, where the reference draws it at random.
- **A3b.** Minimum RH exceeds maximum RH on 0.19% of facility-days. This is in the source data.
- **C3b.** Two sites at The Dalles lie within 0.4 km of an ecoregion boundary.
- **D1d, D1e.** The exceedance threshold's source (ASHRAE A2) isn't stated, and Table 1 isn't described.

## Confirmed correct

- **Raw inputs, end to end (A2):** units, grid cell, longitude convention, dates, and the order of maximum and minimum RH.
- **Data consistency (A3–A5):** standard calendars, correct period years, temperature and humidity on identical time axes and grids, correct humidity caches, correct scenario and member for all 240 source files.
- **Model code (B1–B4):**
  - All 92 parameters are Table B.1 midpoints, in the units the reference code expects.
  - The Case 1 and Case 5 ports reproduce the authors' reference code to 1e-14.
  - The Case 1 indexing correction is verified.
  - The fast backend is within 0.1% of CoolProp on real Oregon weather.
- **Saved outputs (A7c):** the current code reproduces the saved daily files to 1e-8.
- **Analysis code (C1, C2, C4, C5, A8, A9):** Table 1 re-derived independently, no stale caches, the same 20 models throughout, and correct elevation, pressure and grid-cell sampling.
