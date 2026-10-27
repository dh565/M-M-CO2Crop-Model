<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Reference simulation of Part 1 — `WDF25_Ruta_aCO2_WW`

This directory is **the reference simulation of M&M-CO2Wheat v1.0 Part 1**: the
WheatDryFACE 2024/25 experiment, cultivar **Ruta**, **ambient CO2**, **well-watered**
treatment. It is the case to reproduce when checking an installation, and the case
referred to in Part 1 of the companion papers.

| item | value |
|---|---|
| deck name | `WDF25_Ruta_aCO2_WW` |
| site | WheatDryFACE, Rehovot, Israel (31.8953 N, 34.8106 E) |
| season | 2024/25 |
| cultivar | Ruta |
| CO2 | ambient, 420 umol mol-1 (driver record 22) |
| water | well-watered |
| simulated days | 160 |
| `input/` | 31 files, byte-copied from the frozen v1.0 reference run |
| `output/` | the 88 expected result files, byte-copied from the same arm |
| `run.log` | **produced by the current build of this release** — see "About `run.log`" below |
| `MD5SUMS` | 120 digests: 31 input files, 88 output files, `run.log`. It covers the reference DATA only; this `README.md` and `MD5SUMS` itself are excluded, so editing the text cannot make the integrity check fail. |

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
cd testcase && md5 -r $(find . -type f ! -name MD5SUMS ! -name README.md | sed 's|^\./||' | LC_ALL=C sort) \
  | sed 's|^\([0-9a-f]*\) |\1  |' | diff - MD5SUMS && echo "reference intact"
```

The digests are in the standard `<md5>  <path>` format, so on Linux the same check is
simply:

```
cd testcase && md5sum -c MD5SUMS
```

## Compare it

```
python3 tools/compare_outputs.py run/output testcase/output
```

This compares all 88 files your run produced against the 88 shipped expected files. For
each file it prints whether it is byte-for-byte identical and, where it is not, the
maximum absolute and maximum relative difference per numeric column with the row at which
each maximum occurs, then a summary table ordered worst first and a final count.
Requirements: `python3` 3.6 or newer, standard library only.

To compare the console log as well:

```
cmp run/run.log testcase/run.log && echo "run.log identical"
```

## What a pass looks like

In a matched environment the comparison ends with:

```
88 of 88 expected file(s) byte-identical.
VERDICT: byte-identical reproduction of the reference case.
```

The two tolerance statements — what byte identity does and does not mean on another
compiler or platform — are in the repository `README.md`, section "Verifying a run". Read
them before interpreting any difference.

## About `run.log`

The `run.log` shipped here was produced by **the current build of this release**, not by
the frozen reference run that supplied the 88 expected output files. The frozen log has
321 lines, this one has 295, and the whole difference is **console text, with no numeric
value changed**:

| difference | lines | what it is |
|---|---|---|
| removed | 24 | `CANOPY_DBG hr/j:` lines from a hard-coded day-15 debug write that was removed from `CANOPY_MODEL.f90` during the release clean-up |
| removed | 2 | a `DEV DEFAULT EB mode = 1` annotation that contradicted the run's actual setting; the energy-balance mode is still reported, four ways, in the configuration banner |
| reworded | 83 | configuration-banner lines: internal cross-reference tags removed, literature shorthand replaced by author-year citations, internal development labels and document file names dropped, and the source file name corrected from its pre-release name to `main_MMCO2Crop.f90` |

The reworded lines carry the same numbers: comparing the `= value` fields of each
reworded line before and after gives no difference, and all 88 output files are
byte-identical between the two runs.

Shipping the older log would present a line this release cannot produce as this release's
own output, so the current log is shipped instead and the difference is stated here. The
only other difference between a log you produce and the one shipped is the `exit=0` line,
which the run command above appends — it is not written by the model.

`run.log` is reported but is **not part of the identity verdict** of the acceptance
harness: it carries the configuration banner and the season summary, which are model
output, but a runner is free to add wall-clock lines to it.
