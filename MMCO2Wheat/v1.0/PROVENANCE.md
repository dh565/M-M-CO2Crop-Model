<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->

# Provenance of this package

## Identity

| item | value |
|---|---|
| model | M&M-CO2Wheat v1.0 (family: M&M-CO2Crop) |
| source | branch `release-v1.0-20261118` of the development repository |
| phosphorus module | not part of the code (it was never evaluated at a site) |
| decks | calibrated once on this model (CroptimizR, control treatments only); grain-mineral availability set on the control treatments |
| executable | md5 `ae912216adbb5164fc3b11767694b7d1`, `-O3 -Wall`, GNU Fortran (Homebrew GCC 16.2.0) 16.2.0, Apple silicon (arm64), macOS 27.2 |
| date | 2026-11-18 |

## What this archive contains, and what it does not

| promised by the papers' code-and-data-availability statement | here? |
|---|---|
| the model code | yes — 20 Fortran files (19 compiled units plus one `include` body) and the `Makefile` |
| the input decks of the evaluated experiments | yes — all 38 FACE decks, 544 input files, under `decks/` |
| the scripts that produce every figure and table | yes — `paper_scripts/`, 59 files: 20 Python files (`make_*.py` drivers, `run_sensitivity.py`, `config.py` and the modules under `lib/`), 9 files under `table4_and_effects/` (7 Python scripts, `run_decks.sh` and its `README.md`: the Table 4 rows and text numbers that the `make_*` scripts do not write), 25 files under `observations/`, `README.md`, `FIGURE_TABLE_MAP.md`, `requirements.txt`, `MANIFEST.md5`, `LICENSE` |

The figure and table scripts change no model file: they are read-only consumers of the companion reference-outputs archive. They
reproduce every figure and table of both papers and of Supplement S2, and `FIGURE_TABLE_MAP.md`
names, for each one, the script and function that draws it and the model and observation files it
reads. Three qualifications are stated rather than left to be discovered:

1. **The WheatDryFACE observations are not in this archive.** They are deposited openly at Zenodo (ICASA-format workbooks, Helman and Mulero 2026a,
   https://doi.org/10.5281/zenodo.24352317;
   the tables in the `OBS_DIR` layout https://doi.org/10.5281/zenodo.23688951) and the reader supplies them
   under `OBS_DIR` in the layout
   `paper_scripts/README.md` documents. Without them, `make_part2_figures.py`,
   `make_part2_tables.py` (Table 4) and the figures and Tables S1-S4 of `make_supplement_s2.py`
   stop with a `FileNotFoundError` that names the missing file. What runs from this archive and
   the tables in `paper_scripts/observations/` alone: Part 2 Table 5, Supplement Table S5 and
   all of Part 1.
2. **The Braunschweig and Maricopa observation tables in `paper_scripts/observations/` are
   derived from two openly published data sets** (Manderscheid et al. 2020,
   doi 10.18174/odjar.v6i0.16397; Kimball et al. 2017, doi 10.18174/odjar.v3i1.15826, data
   doi 10.7910/DVN/YRTQ3Z) and are redistributed here with the source of every file named in
   `paper_scripts/observations/README_observations.md`. Every one of them is at treatment level
   -- means with standard deviations and replicate counts; none carries a plot, ring or
   replicate identifier.
3. **Part 1 Figs. 2-5 and Table 5 additionally need the executable**, because their 125
   perturbation runs are not part of the reference archive; `run_sensitivity.py` performs them
   (about 1.5 h on two cores).

The reference outputs of all 38 decks are a **companion archive**,
`MMCO2Wheat_v1.0_FACE_reference_outputs` (https://doi.org/10.5281/zenodo.24373926), laid out exactly as
`scripts/release_acceptance.py` expects for `--ref-arm`. One deck's outputs ship here, in
`testcase/`, which is what an installation check needs.

## What was verified for this package

On 2026-11-18 the committed tree was extracted to an empty directory, compiled with the
shipped `Makefile`, and the reference case `WDF25_Ruta_aCO2_WW` was run from a copy of
`testcase/input`:

```
exit=0
90 of 90 expected file(s) byte-identical.
VERDICT: byte-identical reproduction of the reference case.
```

The executable of that build has the md5 given above, which is also the one recorded in
`testcase/README.md`. The headline values of the run (last day of `output/daily_growth.csv`)
are grain yield 426.760 g m-2, above-ground biomass 1335.27 g m-2 and LAI 0.0648.
`testcase/MD5SUMS` (112 digests) and `MD5SUMS` were checked after the last edit.

Not re-measured for this package: runs of the other 37 decks, and the runtime-check
build (`-O0 -fcheck=all`).

## Notes on the shipped input files

`decks/*/input/driver_main.txt` carries 198 positional records; `testcase/input/driver_main.txt`,
`decks/WDF25_Ruta_aCO2_WW/input/driver_main.txt` are byte-identical, and
`docs/V1_CONFIGURATION_driver_main_v1.0.txt` has the same 198 record values (only the comments of ten
records differ, and its last line states that the record values define v1.0). Some of the comment text
inside the input files and the `.f90` comments refers to development steps; no comment is
part of any value the model parses.

## Integrity

`MD5SUMS` lists every file shipped except itself and `testcase/MD5SUMS`. Check it with
`md5sum -c MD5SUMS` (Linux) or the equivalent on macOS. `testcase/MD5SUMS` covers the
reference data only, so editing `testcase/README.md` cannot make that check fail.
