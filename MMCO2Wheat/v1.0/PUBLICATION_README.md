<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# M&M-CO2Wheat v1.0 — publication package

This is the archive that accompanies the two companion papers. It holds the model source,
the input decks of every experiment the papers evaluate, one fully reproducible reference
simulation, and the acceptance harness for deck-by-deck comparison. Its provenance is in
`PROVENANCE.md`; the phosphorus module is not part of it.

| what | where | count |
|---|---|---|
| model source, Fortran 2008 | `*.f90`, `Makefile` | 20 source files (19 compiled units + one `include` body), 1 build file |
| documentation | `README.md`, `INPUT_GUIDE.md`, `docs/` | 2 + 7 |
| input decks of the evaluated experiments | `decks/<deck>/input/` | 38 decks, 544 input files |
| reference simulation for verification | `testcase/` | 21 inputs, 90 expected outputs, `run.log` |
| acceptance harness and comparison tools | `scripts/` (3 files), `tools/` (22 files) | 25 |
| figure and table scripts of both papers | `paper_scripts/` | 50 files (20 Python, 25 under `observations/`, 5 docs/manifest) |
| integrity | `MD5SUMS`, `testcase/MD5SUMS` | — |
| provenance | `PROVENANCE.md` | — |

## Verify the installation in three commands

Exactly the procedure in `docs/START_HERE.md`, which also explains what a difference on
another compiler does and does not mean:

```
make
mkdir -p run && cp -R testcase/input run/ && mkdir -p run/output && (cd run && ../mmco2crop > run.log 2>&1; echo "exit=$?" >> run.log) && tail -1 run/run.log
python3 tools/compare_outputs.py run/output testcase/output
```

In a matched environment the last command ends with `90 of 90 expected file(s)
byte-identical.` (measured on 18 Nov 2026, `PROVENANCE.md`).

## Run any of the evaluated experiments

Each deck directory is self-contained: it holds the exact input bytes that produced the
results in the papers. The model reads `input/...` and writes `output/...` relative to the
directory it is started in, so a deck is run by starting the executable inside its own
directory:

```
mkdir -p decks/<DECK>/output && (cd decks/<DECK> && /full/path/to/mmco2crop > run.log 2>&1; echo "exit=$?" >> run.log)
```

`decks/DECK_MANIFEST.csv` lists all 38 decks with their site. The sites are the three FACE
experiments the papers evaluate: Braunschweig (16 decks), Maricopa (10), WheatDryFACE
2024/25 (8) and WheatDryFACE 2022/23 (4).

`scripts/release_acceptance.py` is the acceptance harness: it runs a whole reference arm
deck by deck and compares every output file byte for byte. It requires a reference arm —
a directory of deck folders each holding `input/`, `output/` and `run.log`. This package
ships the reference **outputs** for one deck only (`testcase/`); the other 37 are in the
companion archive `MMCO2Wheat_v1.0_FACE_reference_outputs` (https://doi.org/10.5281/zenodo.24373926), which is
laid out exactly as the harness expects, so it can be passed straight to `--ref-arm`.

## Observation data

Some deck inputs are observations, not model settings. They are a subset of the **processed**
WheatDryFACE observations — treatment means and standard deviations, campaign-date means and
daily ring-mean series — and **no plot-level observation is in this archive**; plot-level and raw data
are open in the WheatDryFACE records of Helman and Mulero (2026a, https://doi.org/10.5281/zenodo.24352317;
2026b, https://doi.org/10.5281/zenodo.24352486). Those files carry the licence of their own source
archive (CC BY 4.0 for the WheatDryFACE observations, described by Mulero et al. 2026a,b; the
Maricopa and Braunschweig observations are the open data sets of Kimball et al. 2017 and
Manderscheid et al. 2020), not the software licence below. `README.md`, section "Observation
data in this package", lists each file and its level.

## Licence and citation

`LICENSE` and `CITATION.cff`. The model is M&M-CO2Wheat v1.0, part of the M&M-CO2Crop
family. The two companion papers are submitted to Geoscientific Model Development; the DOI of this
package is https://doi.org/10.5281/zenodo.23687586.
