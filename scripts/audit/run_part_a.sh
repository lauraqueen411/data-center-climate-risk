#!/usr/bin/env bash
# Part A of docs/tasks/methodology-audit.md: the server-side checks.
# Read-only: nothing here changes model outputs. Results go to outputs/audit_part_a/.
#
# Run on the server from the repo root, with the climate-risk env active:
#   nohup bash scripts/audit/run_part_a.sh > outputs/audit_part_a.log 2>&1 &
#   tail -f outputs/audit_part_a.log
#
# Each check runs even if an earlier one fails; failures are listed at the end.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src:$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}"
PY="${PY:-python}"
OUT="outputs/audit_part_a"
mkdir -p "$OUT"
T0=$(date +%s)
FAILED=()
stamp() { echo "[$(date '+%F %T')] $*"; }

run() {  # name, command...
  local name="$1"; shift
  stamp "=== $name"
  local t=$(date +%s)
  if "$@"; then stamp "    $name ok ($(( $(date +%s) - t ))s)"; else stamp "    $name FAILED"; FAILED+=("$name"); fi
}

a1() {
  {
    echo "file,line"
    grep -rn "Sanity check WARNING" outputs/ensemble_pue_wue_facilities/logs rerun_case1_case5.log nohup.out \
      --include='*.log' --include='nohup.out' 2>/dev/null | sed 's/,/;/g; s/:/,/'
  } > "$OUT/a1_sanity_warnings.csv"
  echo "  $(( $(wc -l < "$OUT/a1_sanity_warnings.csv") - 1 )) sanity warnings in run logs"
  echo "  other logs that mention warnings or errors:"
  grep -rlE "WARNING|Traceback|Error" outputs/ensemble_pue_wue_facilities/logs 2>/dev/null | head -20 | sed 's/^/    /'
  ls -la outputs/ensemble_pue_wue_facilities/*coolprop* 2>/dev/null | sed 's/^/  coolprop file: /'
}

run "A1 sanity warnings" a1
run "A2 raw inputs end to end" "$PY" scripts/audit/a2_raw_inputs.py
run "A3 consistency and calendars" "$PY" scripts/audit/a3_consistency.py
run "A4/A5 grids, time axes, source files" "$PY" scripts/audit/a4_a5_grids_files.py
run "A6 fast vs CoolProp" "$PY" scripts/audit/a6_fast_vs_coolprop.py
run "A7 chiller COP extrapolation" "$PY" scripts/audit/a7_cop_extrapolation.py
run "A8/A9 elevation, heat ensemble" "$PY" scripts/audit/a8_a9_elevation_heat.py

stamp "Part A finished in $(( ($(date +%s) - T0) / 60 )) min. Failed: ${FAILED[*]:-none}"
stamp "Sync $OUT/ and outputs/audit_part_a.log back to the laptop."
