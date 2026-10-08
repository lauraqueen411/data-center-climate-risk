#!/usr/bin/env bash
# Rerun the facility PUE/WUE ensemble for Cases 1 and 5 after the Table B.1
# midpoint corrections (Case 1 rh_up_pct 60 -> 77.5; Case 5 delta_t_water_c
# 5 -> 7.5 and fan_pressure_ct_pa 250 -> 300), then rebuild the summaries the
# figures read. Case 2 is unchanged and is not rerun.
#
# Run on the server from the repo root, with the climate-risk env active:
#   nohup bash scripts/rerun_case1_case5.sh > rerun_case1_case5.log 2>&1 &
#   tail -f rerun_case1_case5.log
#
# Steps:
#   0. Preflight: corrected defaults are the ones Python imports; 180 old files present.
#   1. RH-setpoint diagnostic on the existing Case 1 weather (no model run).
#   2. Move the old Case 1 and Case 5 daily files, and copy the old analysis
#      CSVs, to outputs/ensemble_pue_wue_facilities/superseded_rh_up60_2026-10-08/.
#   3. Rerun Cases 1 and 5 for the 20 ensemble models, JOBS models at a time.
#   4. Check that all 180 daily files exist again.
#   5. Rebuild analysis/ (build_pue_wue_ensemble_summary.py) and the Figure 2
#      per-model cache (fig_show_everything.py --refresh-cache).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
PY="${PY:-python}"
JOBS="${JOBS:-4}"
FAC="outputs/ensemble_pue_wue_facilities"
OLD="$FAC/superseded_rh_up60_2026-10-08"
LOGS="$FAC/logs/rerun_2026-10-08"
T0=$(date +%s)
stamp() { echo "[$(date '+%F %T')] $*"; }

# ---- 0. Preflight -----------------------------------------------------------
stamp "Preflight"
"$PY" - <<EOF
import climate_risk_dc, sys
from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case5Config
assert climate_risk_dc.__file__.startswith("$ROOT/src"), f"climate_risk_dc imported from {climate_risk_dc.__file__}, not $ROOT/src"
c1, c5 = Case1Config(), Case5Config()
assert (c1.rh_up_pct, c5.delta_t_water_c, c5.fan_pressure_ct_pa) == (77.5, 7.5, 300.0), "corrected defaults not loaded"
print(f"OK: rh_up_pct={c1.rh_up_pct}, Case 5 delta_t_water_c={c5.delta_t_water_c}, fan_pressure_ct_pa={c5.fan_pressure_ct_pa}")
EOF
if [ -d "$OLD" ]; then
  echo "ERROR: $OLD already exists -- this script has run before. Inspect it before rerunning." >&2
  exit 1
fi
n_old=$(ls "$FAC"/facility_pue_wue_*_fast.csv | wc -l)
[ "$n_old" -eq 180 ] || { echo "ERROR: expected 180 daily files, found $n_old" >&2; exit 1; }
GCMS=$("$PY" -c "import pandas as pd; m = pd.read_csv('outputs/ensemble_pue_wue/gcm_manifest.csv'); print(' '.join(m.loc[m['included_in_ensemble'] == 'yes', 'model']))")
n_gcm=$(echo "$GCMS" | wc -w)
[ "$n_gcm" -eq 20 ] || { echo "ERROR: expected 20 ensemble models in the manifest, found $n_gcm" >&2; exit 1; }
stamp "Models: $GCMS"

# ---- 1. Diagnostic (weather only, unaffected by the rerun) -----------------
stamp "RH-setpoint diagnostic"
"$PY" scripts/diagnose_rh_up_sensitivity.py --repo-root . | tail -3

# ---- 2. Move old outputs aside ---------------------------------------------
stamp "Moving old Case 1 and Case 5 files to $OLD"
mkdir -p "$OLD/analysis" "$LOGS"
mv "$FAC"/facility_pue_wue_*_ssp370_ae-chiller_*_fast.csv "$FAC"/facility_pue_wue_*_ssp370_chiller-only_*_fast.csv "$OLD"/
cp -p "$FAC"/analysis/*.csv "$OLD/analysis/"
[ -d "$FAC/analysis/fig_show_everything" ] && cp -rp "$FAC/analysis/fig_show_everything" "$OLD/analysis/"
stamp "Moved $(ls "$OLD"/*.csv | wc -l) daily files (expected 120)"

# ---- 3. Rerun ---------------------------------------------------------------
stamp "Running Cases 1 and 5, $JOBS models at a time; per-model logs in $LOGS"
export PY LOGS
set +e
echo "$GCMS" | tr ' ' '\n' | xargs -P "$JOBS" -I{} bash -c \
  '"$PY" scripts/run_pue_wue_facilities.py --gcms {} --archetypes ae-chiller chiller-only --repo-root . > "$LOGS/{}.log" 2>&1 \
   && echo "[$(date "+%F %T")] done {}" || echo "[$(date "+%F %T")] FAILED {} (see $LOGS/{}.log)"'
set -e

# ---- 4. Check ----------------------------------------------------------------
n_new=$(ls "$FAC"/facility_pue_wue_*_fast.csv | wc -l)
stamp "Daily files present: $n_new (expected 180)"
if [ "$n_new" -ne 180 ]; then
  echo "ERROR: incomplete rerun. Rerunning this script is blocked by $OLD; rerun the missing models by hand:" >&2
  echo "  $PY scripts/run_pue_wue_facilities.py --gcms <model> --archetypes ae-chiller chiller-only --repo-root ." >&2
  exit 1
fi

# ---- 5. Rebuild summaries ---------------------------------------------------
stamp "Rebuilding analysis/ summaries"
"$PY" scripts/build_pue_wue_ensemble_summary.py --repo-root . > "$LOGS/build_summary.log" 2>&1
stamp "Rebuilding Figure 2 per-model cache"
"$PY" scripts/fig_show_everything.py --repo-root . --refresh-cache > "$LOGS/fig_show_everything.log" 2>&1

stamp "All done in $(( ($(date +%s) - T0) / 60 )) min. Sync $FAC/analysis/ back to the laptop."
