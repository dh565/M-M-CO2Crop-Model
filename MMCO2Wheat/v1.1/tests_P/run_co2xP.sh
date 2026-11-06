#!/bin/bash
# CO2 x P test (design section 9, step 4). Runs the reference case
# (WDF25_Ruta_aCO2_WW) with the P module on, at a range of soil-test P and at
# two CO2 levels, plus one P-off run. Usage: [EPS=0.3] [METHOD=2] [ECASE=0] run_co2xP.sh <run_dir>
# EPS sets a measured eps_test (default -1 = use the default eps of METHOD: 1 Olsen,
# 2 Bray-1 (Rodriguez soils, default), 3 Mehlich-3, 4 AL; ECASE -1/0/+1 = low/central/high).
# The emergent prediction to test: the eCO2 gain in biomass and grain shrinks
# as soil P falls, because leaf cytosolic Pi limits TPU and P limits construction.
set -e
B=$(cd "$(dirname "$0")/.." && pwd)
R=${1:?run dir}
mkdir -p "$R"
run() {  # name co2 soiltest(-1 = P off)
  d="$R/$1"; rm -rf "$d"; mkdir -p "$d/output"; cp -R "$B/testcase/input" "$d/"
  sed -i '' "22s/^420.0 /$2 /" "$d/input/driver_main.txt"
  if [ "$3" != "-1" ]; then
    sed -e "s/^soiltest_mgkg,[^,]*,/soiltest_mgkg,$3,/" -e "s/^eps_test,[^,]*,/eps_test,${EPS:--1},/" -e "s/^test_method,[^,]*,/test_method,${METHOD:-2},/" -e "s/^eps_case,[^,]*,/eps_case,${ECASE:-0},/" "$B/input_templates/driver_P.csv" > "$d/input/driver_P.csv"
  fi
  (cd "$d" && "$B/mmco2crop" > run.log 2>&1; echo "exit=$?" >> run.log)
  echo "$1 $(tail -1 $d/run.log)"
}
for co2 in 420 550; do
  run "off_${co2}" "$co2.0" -1
  for st in 40 15 8 5 3; do run "P${st}_${co2}" "$co2.0" "$st"; done
done
