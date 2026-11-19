#!/bin/bash
# run_cal.sh <work dir> [<exe>]: the CroptimizR calibration of the package decks (calibrate_face.R), run on copies in
# <work dir>/decks; the package decks are not changed. Steps: initial values (initial_values.csv) written to the copies;
# phase phen (9 control decks); phase cult (WDF25 Ruta/Zahir to kernel number, MARICOPA_906 leaf only, its fruiting
# efficiency is measured); all 38 decks run; grain-mineral calibration step (recentre.py).
# Results: <work dir>/cal/{phen,cult}/calibrated_*.csv and the calibrated decks in <work dir>/decks.
# Needs R with CroptimizR, python3 and the model built with the package Makefile (default <package>/mmco2crop).
set -e
H=$(cd "$(dirname "$0")" && pwd); PKG=$(cd "$H/../.." && pwd)
W=${1:?usage: run_cal.sh <work dir> [<exe>]}; mkdir -p "$W"; W=$(cd "$W" && pwd)
EXE=${2:-$PKG/mmco2crop}; EXE=$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")
V=$W/decks; rm -rf "$V"; mkdir -p "$V"
for d in $(cd "$PKG/decks" && ls -d */ | tr -d /); do cp -R "$PKG/decks/$d" "$V/"; done
export FCAL_HOME=$H FCAL_DECKS=$V FCAL_EXE=$EXE PYTHONDONTWRITEBYTECODE=1
cd "$H"
python3 -I init_values.py "$V" initial_values.csv > "$W/init.log"
Rscript calibrate_face.R phen targets.csv "$W/cal/phen" > "$W/phen.log" 2>&1
python3 -I apply.py "$V" "$W/cal/phen/calibrated_phen.csv" > "$W/apply_phen.log"
KN_DECKS=WDF25_Ruta_aCO2_WW,WDF25_Zahir_aCO2_WW NO_FE=MARICOPA_906 ONLY_DECKS=WDF25_Ruta_aCO2_WW,WDF25_Zahir_aCO2_WW,MARICOPA_906 \
  Rscript calibrate_face.R cult targets.csv "$W/cal/cult" > "$W/cult.log" 2>&1
python3 -I apply.py "$V" "$W/cal/cult/calibrated_cult.csv" > "$W/apply_cult.log"
bash "$PKG/paper_scripts/table4_and_effects/run_decks.sh" "$EXE" "$V" "$W/runs" 8 > "$W/runs.log" 2>&1
echo "runs: $(grep -c 'exit=0' "$W/runs.log") of 38 ok"
python3 -I recentre.py "$W/runs" "$V" > "$W/recentre.log"
echo "calibration done: $W"
