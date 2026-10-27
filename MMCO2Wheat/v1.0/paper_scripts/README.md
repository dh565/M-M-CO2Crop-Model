<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# M&M-CO2Wheat v1.0 — figure and table scripts of the companion papers

These scripts produce every figure and table of

* **Part 1** — model description (Figs. 1–5, Table 5),
* **Part 2** — evaluation at three wheat FACE experiments (Figs. 1–11, Tables 4 and 5), and
* **Supplement S2** of Part 2 (Figs. S1–S8, Tables S1–S5)

from the published reference outputs of the model. No model is run by the Part 2 and Supplement
scripts: every simulated number is read from the reference-outputs archive. Only the Part 1
sensitivity analysis (Figs. 2–5, Table 5) requires the model executable, because its 125
perturbation runs are not part of the archive.

`FIGURE_TABLE_MAP.md` lists, for each figure and table, the script and function that draws it and
the model output files (by deck) and observation files it reads.

## Contents

| path | what it is |
|---|---|
| `config.py` | every path, overridable by environment variables (see below) |
| `make_part2_figures.py` | Part 2, Figs. 1–11 |
| `make_part2_tables.py` | Part 2, Tables 4 and 5 and the numbers quoted in the text |
| `make_supplement_s2.py` | Supplement S2, Figs. S1–S8 and Tables S1–S5 |
| `make_part1_figures.py` | Part 1, Figs. 1–5 |
| `make_part1_tables.py` | Part 1, Table 5 and the numbers of Sects. 13.3 and 14 |
| `run_sensitivity.py` | the 125 sensitivity runs behind Part 1 Figs. 2–5 and Table 5 (needs the model executable) |
| `lib/` | the modules the scripts share: deck registry and model readers (`face_results.py`), observation–model pairing of the energy-balance, soil, gas-exchange and canopy-temperature records (`tleaf_eval.py`, `ebwb_gx.py`), the WheatDryFACE file reader (`wdf_ingest.py`), the Maricopa archive readers (`maricopa_archive.py`, `maricopa_eb.py`, `maricopa_gx.py`, `maricopa_tsoil.py`, `extract_maricopa_eb.py`), the porometer figure (`wdf_poro.py`), the figure code (`paperfigs.py`, `paperfigs2.py`, `p1_figures.py`) |
| `observations/` | observation tables derived from the open Braunschweig and Maricopa data sets, redistributed here (see `observations/README_observations.md`) |
| `requirements.txt` | Python packages, pinned |
| `MANIFEST.md5` | md5 of every file in this folder |

## Software

The figures of the papers were made with Python 3.13.16, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1,
Matplotlib 3.11.2, odfpy 1.4.1 and cmcrameri 1.10 on Linux, with the fonts Liberation Sans /
DejaVu Sans. Any Python 3.9+ with the packages of `requirements.txt` reproduces the content of
every figure; a pixel-identical PNG additionally needs Matplotlib 3.11.2 and the same fonts.

## The three inputs a reader must provide

### 1. The reference-outputs archive (`REF_OUTPUTS`)

`MMCO2Wheat_v1.0_FACE_reference_outputs` — 38 deck folders, each with `input/`, `output/`
(88 files) and `run.log`, as published with the code package `MMCO2Wheat_v1.0` (source commit
`1cab3d55` with the 27 Oct 2026 re-seal, executable md5 `65eb0273a9a543ab807d9f0920072bf1`; see `PROVENANCE.md`). Every Part 2 and Supplement
number comes from these files.

### 2. The model executable (`MODEL_EXE`) — Part 1 Figs. 2–5 and Table 5 only

Built from the code package with its `Makefile`. `run_sensitivity.py` copies the archive deck
`WDF25_Ruta_aCO2_WW/input` (the reference simulation of Part 1) 125 times, applies one input
change per copy and runs the executable in each (about 45 s per run on one core), then collects
the harvest and season metrics of every run in `SENS_DIR/all_metrics.json`. A build in another
environment may differ from the published numbers at the last printed digit (see the code
package's `docs/START_HERE.md`).

### 3. The observation data (`OBS_DIR`)

Observations from the two open data sets are redistributed in `observations/` and need no
action. The WheatDryFACE files the scripts read are the **processed** per-treatment files of
the experiment (treatment means and standard deviations over the replicate rings, daily
sensor series, campaign means, porometer and thermal-image means); these are deposited openly
with the WheatDryFACE data set (the plot-level raw data stay under embargo until 2029-01-01).
They must be placed under `OBS_DIR` with exactly this layout and these names:

```
OBS_DIR/
  WheatDryFACE/
    WDF25/                   S3_<genotype>_<CO2>_<water>.csv   8 files   (season 2024/25; season record,
                                                               thermal-image and porometer columns)
    WDF23/                   WDF23_<genotype>_<CO2>.csv        4 files   (season 2022/23)
    obs_registry_wheatdryface.csv                              per-treatment harvest and canopy scalars
                                                               (obs, obs_sd), a table derived from the files
                                                               above and deposited with them; if it is absent
                                                               the scripts print a warning and read the same
                                                               scalars from the treatment files
  Maricopa_FACE_archive/                                       the workbooks of Kimball et al. (2017),
                                                               data doi 10.7910/DVN/YRTQ3Z, as deposited:
    Weather Canopy & Soil Temperatures Energy Balance 1994 hourly.ods
    Weather Canopy & Soil Temperatures Energy Balance CO2 1996 15-min.ods
    Weather Canopy & Soil Temperatures Energy Balance CO2 O3 1997 15-min.ods
    Physiology Gas Exchange Water Relations.ods
    Quality Wheat Grain Baking 1993-1997.ods
```

`<genotype>` is `Ruta` or `Zahir`, `<CO2>` is `aCO2` or `eCO2`, `<water>` is `WW` or `WL`.
The deposit `WheatDryFACE_figure_data` has this layout: place (or link) it as `OBS_DIR/WheatDryFACE/`.
Its files carry only the columns the scripts read, with the two-row header of the processed
treatment files (row 1 block labels Night/Day and Grain/Stem/Spike (head)/Leaf/Chaff/Straw, row 2
variable names, each value column followed by its `Std` column) and one row per date; the
deposit's `MD5SUMS` lists every file. The Maricopa
workbooks are read directly (the first read of each takes one to two minutes and is cached under
`OUT_DIR/_cache/`); the five workbooks above are the only ones needed.

**What cannot be reproduced without the WheatDryFACE files.** Every Part 2 figure has WheatDryFACE
panels, and the observation registry is loaded before any Part 2 figure or Table 4 is drawn, so
`make_part2_figures.py`, `make_part2_tables.py` (Table 4) and the figures and Tables S1–S4 of
`make_supplement_s2.py` stop with a `FileNotFoundError` naming the missing file. What runs
from the archive and the redistributed tables alone: Part 2 Table 5
(`make_part2_tables.py --only T5`), Supplement Table S5 (`make_supplement_s2.py --only T5`)
and all of Part 1.

## Running

```
export REF_OUTPUTS=/path/to/MMCO2Wheat_v1.0_FACE_reference_outputs
export OBS_DIR=/path/to/observation_data          # layout above
export OUT_DIR=/path/to/generated_results         # default: ../generated_results next to this folder
python3 make_part2_figures.py                     # Part 2 Figs. 1-11      -> OUT_DIR/part2_figures/ (+ gmd/f01..f11)
python3 make_part2_tables.py                      # Part 2 Tables 4, 5     -> OUT_DIR/part2_tables/
python3 make_supplement_s2.py                     # Supplement S2          -> OUT_DIR/supplement_s2/
python3 make_part1_figures.py --only 1            # Part 1 Fig. 1 (no data needed)
export MODEL_EXE=/path/to/mmco2crop               # built from the code package
python3 run_sensitivity.py                        # 125 runs               -> SENS_DIR/all_metrics.json (about 1.5 h on two cores)
python3 make_part1_figures.py                     # Part 1 Figs. 1-5       -> OUT_DIR/part1_figures/ (+ gmd/f01..f05)
python3 make_part1_tables.py                      # Part 1 Table 5         -> OUT_DIR/part1_tables/
```

`SENS_DIR` defaults to `OUT_DIR/sensitivity_runs`. Each `make_*` script accepts `--only` to
draw a subset (see its docstring). Every figure is written as PDF and PNG (300 dpi; Part 1
Fig. 1 at 400 dpi) under its descriptive name, and the main-paper figures also as `f01.png`…
in a `gmd/` subfolder. Nothing is written into the archive or into this folder.

## Notes on the derived observation tables

* `observations/obs_registry_braunschweig_maricopa.csv` and `OBS_DIR/WheatDryFACE/obs_registry_wheatdryface.csv`
  are the per-treatment harvest and canopy scalars (observed mean and sd) that the evaluation
  registers before the raw treatment files are read; the loaders of `lib/face_results.py` fill
  in only what the registry does not carry.
* `observations/consolidated_obs_braunschweig_maricopa.csv` carries the few observations that
  neither the registry nor the loaders provide (Braunschweig grain Fe and Zn and the 2014/15
  harvest index, Maricopa peak LAI and harvested yield); each row names its source file and column.
* The kernel mass (TKW) of the WheatDryFACE 2024/25 treatments shown in Part 2 Figs. 6, 7 and 11
  and Tables S2 and 4 is the least-square mean by genotype x CO2 x water of the 2024/25 season
  study (Mulero et al., submitted; Supplementary Table S4B), held in `lib/face_results.py`; a
  `TKW` column in the 2024/25 treatment files would take precedence over it (the deposit carries none).
* `observations/maricopa/maricopa_eb_hourly.csv` can be rebuilt from the archive workbooks with
  `python3 lib/extract_maricopa_eb.py` (it writes to `OUT_DIR/observations_rebuilt/`).
