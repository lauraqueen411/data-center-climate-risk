# climate-risk-dc — project instructions for Claude Code

Oregon data-center physical climate risk analysis (heat + flood). `docs/methods-document.md` (full methods) and `docs/code-development-plan.md` (live task checklist) are the source of truth — read both before starting substantive work. The plan tracks what's done vs. open as checkboxes; keep it current as you go.

## Scope boundaries — read before editing

- **In current scope:** heat risk (1.1–1.3) and flood risk (2.1–2.5), per the methods doc.
- **Out of current scope — do not extend:** `exposure/fiber.py`, `exposure/water.py`, `exposure/watersheds.py`, `DC_fiber_exposure.ipynb`, `DC_water_exposure.ipynb`, `test_public_water_sources.py`, and the fiber/water columns emitted by `build_dependency_table.py`. This is legacy work from an earlier, broader project phase. Leave it alone unless explicitly asked to revisit it.
- **1.3 (infrastructure heat dependency) is intentionally lightweight for now**: a heat/dependency coincidence flag, not full IEEE C57.91 / IEEE 738 / ISO 8528-1 derating math. Don't build the full derating formulas unless asked — there's a marked extension point (`derating.py` stub) for that.

## Repo layout

- `src/climate_risk_dc/` — the installable package. `config.py` resolves dataset paths (via optional `config/datasets.yml`, currently not yet created — most scripts still hardcode `root/data/...` paths, which is a known cleanup item). `geo.py` has shared CRS/distance/join helpers — reuse these rather than writing new spatial logic inline. `climate/heat_indices.py` is the most mature module (LOCA2 loading, CDD/exceedance/heatwave indices).
- `scripts/` — argparse-driven entry points: `run_heat_processing.py` (heat rasters + maps), `build_dependency_table.py` (facility → substation → transmission → generator dependency graph, via HIFLD + EIA-860), `plot_im3_data_centers.py` (atlas map).
- `notebooks/` — exploratory work, uneven maturity. Flood work (`JRI_flood_explore.ipynb`, `ISIMIP_flood_explore.ipynb`) is still at the "opens one file and plots it" stage — essentially unstarted.
- `tests/` — currently just one test (`test_public_water_sources.py`). New pure-math modules (heat indices, PUE, Hazus, AAD, derating flag) are good candidates for unit tests on synthetic inputs.
- `data/` and `outputs/` are gitignored — present on disk, not in the repo.

## Conventions to follow

- Frozen `@dataclass` for run/index configuration (see `HeatIndexConfig`, `RunConfig`, `TableConfig`).
- Module docstrings include an "Assumptions log" section documenting simplifications and proxies — follow this pattern for any new module, especially anything involving an engineering approximation (derating rates, damage-function category choice, elevation adjustment rate, etc.). Don't silently pick a number the methods doc doesn't specify — document it as an assumption.
- Scripts are `argparse`-driven, take `--repo-root`, write to `outputs/<phase>/` with rasters/maps/summaries in separate subfolders plus a manifest CSV of run parameters.
- CRS handling is explicit: `EPSG:4326` (geographic) for I/O and web mapping, `EPSG:5070` (`CRS_EQUAL_AREA` in `config.py`) for any distance/area computation — always reproject via `geo.py`'s helpers rather than computing distances in degrees.
- Raster point-sampling: reuse `_sample_raster_values` (currently in `build_dependency_table.py`) rather than writing a new sampler — this belongs in a shared module (`geo.py`) if it isn't already.
- Unit sanity checks on load (see `_log_unit_sanity` for Kelvin bounds) — apply the same pattern to new data sources (flood depth bounds, humidity 0–100%, PUE plausible range, etc.).

## Data source notes

- **JRC flood hazard data** = JRC Global River Flood Hazard Maps v2.1 (Baugh et al. 2024), distributed as `CEMS-GLOFAS` — this is why files are named `glofas_rp*.tif`, not a different/wrong dataset. 90m resolution, return periods RP10/20/50/75/100/200/500. Only RP100 has been pulled for Oregon so far; the other six are still needed.
- LOCA2 heat pipeline is currently locked to a single GCM/member (ACCESS-CM2, r1i1p1f1) — a documented limitation, not an oversight. Don't silently switch to or add an ensemble without flagging it.

## Workflow

- After completing checklist items, update `docs/code-development-plan.md`'s checkboxes to reflect what's actually done.
- Flag open engineering/methodological decisions (assumption rates, scope calls) back to the user rather than resolving them silently — several are intentionally left open in the plan.
