<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# START HERE — M&M-CO2Wheat v1.0

Three steps, in this order: **compile**, **run the test case**, **compare the results**.
Paste each command from the top directory of this package — the directory that holds
`Makefile`.

Requirements: `gfortran` and `make` for the model, `python3` (3.6 or newer, standard
library only) for the comparison. Nothing else: no libraries, no netCDF, no Python for
the model itself.

## 1. Compile

```
make
```

## 2. Run the test case

```
mkdir -p run && cp -R testcase/input run/ && mkdir -p run/output && (cd run && ../mmco2crop > run.log 2>&1; echo "exit=$?" >> run.log) && tail -1 run/run.log
```

## 3. Compare the results

```
python3 tools/compare_outputs.py run/output testcase/output
```

That is the whole procedure. Notes on each step:

- **Step 1** builds the executable `mmco2crop` in the current directory from the 19
  Fortran source files beside the `Makefile` — 18 compiled units and one `include` body,
  `csv_file_1d.f90` — with the release flags `-O3 -Wall`. It takes
  well under a minute. Compiler **warnings are normal and expected**; what matters is that
  `make` finishes with no error and that `mmco2crop` exists.
- **Step 2** runs the Part 1 reference simulation, `WDF25_Ruta_aCO2_WW` (WheatDryFACE
  2024/25, cultivar Ruta, ambient CO2, well watered; 160 simulated days, about half a
  minute). The model reads its inputs from the relative path `input/...` and writes its
  results to the relative path `output/...`, so it must be run from inside a case
  directory and `output/` must exist before it starts — the command above does both. The
  last line printed should be `exit=0`. `testcase/` itself is never written to: the work
  happens on a copy in `run/`. Keep `testcase/` unchanged — it is the reference, and
  `testcase/MD5SUMS` lets you prove at any time that it still is.
- **Step 3** compares all 88 result files your run produced against the 88 shipped
  expected files. For every file it reports byte identity or, where the file is not
  identical, the maximum absolute and maximum relative difference of each numeric column
  with the data row at which each maximum occurs, then a summary table ordered worst
  first and a count per verdict class. Add `--csv differences.csv` to get the per-file
  table as a file. The model's console log can be compared as well, with
  `cmp run/run.log testcase/run.log`. In a matched environment that is byte-identical
  too: the shipped log already carries the `exit=0` line that the run command appends, and
  this was measured for this package.

In a matched environment step 3 ends with:

```
88 of 88 expected file(s) byte-identical.
VERDICT: byte-identical reproduction of the reference case.
```

---

# The two tolerance statements — read before interpreting any difference

These are two different claims and they must never be blurred into one. They are stated
in full in the repository `README.md`, section "Verifying a run"; in short:

**1. In a matched environment, reproduction is byte-for-byte, and that was measured.**
A matched environment means the same compiler, the same compiler version, the same
optimization flags and the same platform. In that environment all 88 output files of the
test case are reproduced byte for byte. This is not an expectation: it is the result of
the release acceptance test, which required all 42 release decks to exit 0 and all 3 696
output files to be byte-identical to the frozen reference package, and which was re-run on
every release commit. The environment in which it was measured for this package is
recorded in `README.md`, section "Dependencies, and the versions this release was tested
with", and in `PROVENANCE.md`.

**2. On a different compiler, compiler version, optimization level or platform, byte
identity is NOT expected.** One further platform has been measured and is reported in
Part 1, Sect. 13.3 (gfortran 13.3 on Linux x86-64: 55 of 88 output tables byte-identical);
a single second platform does not define a tolerance. Different
compilers and versions legitimately reorder floating-point operations, contract
multiply-add pairs differently and evaluate intrinsics to different last bits; in a
season-long iterative simulation those last-bit differences grow. No measured
cross-platform tolerance can be quoted for this model, because only one platform was
available when it was released, so no second platform could be measured and any tolerance
stated here would be an invention. `tools/compare_outputs.py` exists precisely to
quantify whatever difference you actually get: run it and read the numbers it prints. A
difference it reports on a different compiler is information about compiler sensitivity,
not evidence that your installation is wrong, and the right response is to report the
numbers rather than to adjust anything.

---

## Where to look next

| you want | read |
|---|---|
| what the model is and how it is built | `README.md` |
| every input file, its format and units, and every driver record | `INPUT_GUIDE.md` |
| what each module does | `docs/MODULES.md` |
| the module call and data-flow diagram | `docs/ARCHITECTURE.png` (source: `.dot`) |
| every file that ships, every input and all 88 outputs | `docs/FILE_INVENTORY.md` |
| optional diagnostics that v1.0 cannot produce | `docs/V1_KNOWN_LIMITATIONS.md` |
| how the shipped architecture diagram is regenerated | `scripts/render_architecture.py` |
| the v1.0 driver configuration, record by record | `docs/V1_CONFIGURATION_driver_main_v1.0.txt` |
| the test case itself | `testcase/README.md` |
| licence and how to cite | `LICENSE`, `CITATION.cff` |

If any of the three commands fails, send the full console output of the command that
failed, together with the version lines of your compiler, `make` and `python3`, to the
contact address in `README.md`.
