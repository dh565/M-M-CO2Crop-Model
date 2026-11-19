<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Reference simulation of Part 1 — `WDF25_Ruta_aCO2_WW`

This directory is **the reference simulation of M&M-CO2Wheat v1.0** (release of
18 Nov 2026): the WheatDryFACE 2024/25 experiment, cultivar **Ruta**, **ambient CO2**,
**well-watered** treatment. It is the case to reproduce when checking an installation. It
was generated on 18 Nov 2026 from the reference deck `decks/WDF25_Ruta_aCO2_WW`
with the build of this release (gfortran 16.2.0, Apple silicon, `-O3 -Wall -ffree-line-length-none`; executable md5 ae912216adbb5164fc3b11767694b7d1), after the one calibration of the decks on this model (18 Nov 2026).

| item | value |
|---|---|
| deck name | `WDF25_Ruta_aCO2_WW` |
| site | WheatDryFACE, Rehovot, Israel (31.8953 N, 34.8106 E) |
| season | 2024/25 |
| cultivar | Ruta |
| CO2 | ambient, 437.3 umol mol-1 (driver record 22; measured daytime mean of the ambient rings, WDF ICASA delivery) |
| water | well-watered |
| simulated days | 160 |
| `input/` | 21 files, byte-copied from `decks/WDF25_Ruta_aCO2_WW/input` at this commit |
| `output/` | the 90 result files written by the v1.0 build on this input |
| `run.log` | the console log of the same run, with the runner's `exit=0` line appended |
| `MD5SUMS` | 112 digests: 21 input files, 90 output files, `run.log`. It covers the reference DATA only; this `README.md`, `compare.py` and `MD5SUMS` itself are excluded, so editing the text cannot make the integrity check fail. |
| `compare.py` | the acceptance check: compares a run with `output/` and gives IDENTICAL, PASS, DIFFERENT-DIAGNOSTIC or FAIL per file (section "Compare it") |
| headline values (last day of `output/daily_growth.csv`) | grain yield 426.760 g m-2, above-ground biomass 1335.27 g m-2, LAI 0.0648 |

## Run it

From the directory that holds the `Makefile` (one level above this one):

```
make
mkdir -p run && cp -R testcase/input run/ && mkdir -p run/output \
  && (cd run && ../mmco2crop > run.log 2>&1; echo "exit=$?" >> run.log) \
  && tail -1 run/run.log
```

The model reads its inputs from the relative path `input/...` and writes its results to
the relative path `output/...`, so it must be run from inside a case directory and
`output/` must exist before it starts; the command above does both. The last line printed
should be `exit=0`. The run takes well under a minute.

`testcase/` itself is never written to: the command works on a copy in `run/`. Keep this
directory unchanged — it is the reference, and `MD5SUMS` lets you prove at any time that
it still is:

```
cd testcase && md5 -r $(find . -type f ! -name MD5SUMS ! -name README.md ! -name compare.py | sed 's|^\./||' | LC_ALL=C sort) \
  | sed 's|^\([0-9a-f]*\) |\1  |' | diff - MD5SUMS && echo "reference intact"
```

The digests are in the standard `<md5>  <path>` format, so on Linux the same check is
simply:

```
cd testcase && md5sum -c MD5SUMS
```

## Compare it

```
python3 testcase/compare.py run/output testcase/output
```

This is the acceptance check. It compares all 90 files your run produced with the 90
shipped expected files. Requirements: `python3` 3.6 or newer, standard library only.

For each file, `compare.py` prints one of:

- **IDENTICAL**: the file is the same byte for byte.
- **PASS**: the file differs, but every number is within tolerance (see below).
- **DIFFERENT-DIAGNOSTIC**: the only numbers outside tolerance are in the diagnostic columns listed below. Each such cell is listed by file, column and row. This is not a failure.
- **FAIL**: a number outside a diagnostic column is outside tolerance, or the rows, columns or text differ.
- **MISSING**: the file was not produced.

Each line also gives the largest absolute and the largest relative difference in that file. The last line is `OVERALL: PASS` or `OVERALL: FAIL`. OVERALL is PASS only if no file is FAIL or MISSING.

A number is within tolerance if any of these holds:

1. The relative difference is at most 1e-6: |a − b| ≤ 1e-9 + 1e-6 · max(|a|, |b|).
2. The difference is at most half a unit of the last printed digit. For example, this is 0.000005 for `2.26563`, and 0.00000005 for `0.585938E-01`.
3. The number is in a near-zero residual column listed below, and the difference is within that column's absolute tolerance.

Near-zero residual columns:

| file | column | tolerance | reason |
|---|---|---|---|
| `daily_soil_heat.csv`, `hourly_soil_heat.csv` | `heat_residual_J_m2` | 1e-3 J m-2 | The column is the soil heat closure: the surface heat flux integral minus the change in heat storage. Each term is typically 1e5 to 1e6 J m-2 (median 4e5 J m-2 per day and 1.6e5 J m-2 per hour). Rounding leaves 1e-8 to 1e-6 J m-2, so the relative test cannot be used on a value that is zero up to rounding. The tolerance of 1e-3 J m-2 is 1e-8 of the median daily term and still shows any real closure error. |

Diagnostic columns. We checked in the source code that each is written to output only, and that no state variable reads it.

| file | column | reason |
|---|---|---|
| `daily_n_construction.csv` | `c_actual_L`, `c_actual_S`, `c_actual_R` | The nitrogen concentration of the tissue built that day: the N paid divided by the dry matter built (`growth_module.f90`). On a day when the increment is a rounding remainder (about 1e-14 g m-2) instead of exactly 0, the ratio is the target concentration instead of 0. Growth is not affected. |
| `hourly_energy_balance.csv` | `Tac_minus_Tair_demand`, `Tcan_minus_Tair_demand` | The difference of two values that `compare.py` checks in their own columns (`Tac_demand_C` or `Tcan_demand_C`, and `Tair_C`). Some of these differences sit exactly on a rounding tie of the printed digits, such as 2.265625 and 0.05859375, so a last-bit difference changes the last printed digit. |

The tolerances are fixed. Do not change them; report the full table instead.

For a detailed report per numeric column (largest absolute and relative difference and the
row of each), use `python3 tools/compare_outputs.py run/output testcase/output`.

To compare the console log as well:

```
cmp run/run.log testcase/run.log && echo "run.log identical"
```

## What a pass looks like

With the build described above, all 90 files are IDENTICAL and the check ends with
`OVERALL: PASS`. With another compiler, compiler version or platform, byte identity is not
expected (repository `README.md`, section "Verifying a run"); the check then reports PASS
or DIFFERENT-DIAGNOSTIC for the files that differ, and `OVERALL: PASS` if no file is FAIL
or MISSING.

## About `run.log`

The `run.log` shipped here was written by the same run that produced `output/`. The only
difference between a log you produce and the one shipped is the `exit=0` line, which the
run command above appends — it is not written by the model.

`run.log` is reported but is **not part of the identity verdict** of the acceptance
harness: it carries the configuration banner and the season summary, which are model
output, but a runner is free to add wall-clock lines to it.
