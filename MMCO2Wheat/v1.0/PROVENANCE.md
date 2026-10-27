<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->

# Provenance of this package

## Identity

| item | value |
|---|---|
| model | M&M-CO2Wheat v1.0 (family: M&M-CO2Crop) |
| source commit | `1cab3d557391058e784e8932168aa3dfb9c7d3ea` (25 Oct seal); the re-sealed source of this package is the content tagged `v1.0.0` in the public repository |
| frozen original reference | D40, commit `a4810b49`, tag `v1.0.0-rc3` |
| executable, release build | md5 `65eb0273a9a543ab807d9f0920072bf1`, `-O3 -Wall` (re-seal of 2026-10-27; the 2026-10-25 seal was `f4f6f352708609f3521c6347e06dbbb4`, see below) |
| executable, runtime-check build | md5 `a8b96c16c71aacc7bcbdbe0cdafeabcd`, `-O0 -fcheck=all -Wall` (2026-10-25 seal; not re-measured after the re-seal) |
| build environment of the measured results | GNU Fortran (Homebrew GCC 9.3.0) 9.3.0 emitting x86_64 under Rosetta 2, on Darwin 25.6.0, `uname -m` = `arm64` (Apple silicon) |
| sealed | 2026-10-25; re-sealed 2026-10-27 (see "Re-seal of 27 October 2026") |

## What this archive contains, and what it does not

| promised by the papers' code-and-data-availability statement | here? |
|---|---|
| the model code | yes — 19 Fortran files, 18 compiled units plus one `include` body, and the `Makefile` |
| the input decks of the evaluated experiments | yes — all 38 FACE decks, 598 input files, under `decks/` |
| the scripts that produce every figure and table | yes — `paper_scripts/`, 50 files: 5 `make_*.py` drivers, `run_sensitivity.py`, 14 modules under `lib/`, 23 derived observation tables under `observations/`, `config.py`, `README.md`, `FIGURE_TABLE_MAP.md`, `requirements.txt`, `MANIFEST.md5`, `LICENSE` |

The figure and table scripts were added on 2026-10-26, after the source seal, and changed no
model file: they are read-only consumers of the companion reference-outputs archive. They
reproduce every figure and table of both papers and of Supplement S2, and `FIGURE_TABLE_MAP.md`
names, for each one, the script and function that draws it and the model and observation files it
reads. Three qualifications are stated rather than left to be discovered:

1. **The WheatDryFACE observations are not in this archive.** They are deposited openly with
   the WheatDryFACE data set and the reader supplies them under `OBS_DIR` in the layout
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

The full reference outputs of all 38 decks are a **companion archive**,
`MMCO2Wheat_v1.0_FACE_reference_outputs` (3 344 output files, 1.3 GB uncompressed, with its
own `MD5SUMS`), laid out exactly as `scripts/release_acceptance.py` expects for `--ref-arm`.
One deck's outputs ship here, in `testcase/`, which is what an installation check needs.

The 4 decks of a separate harmonized-input exercise are part of the 42-deck
release-acceptance arm. No result from them is reported in either paper; Part 1 names
them in Sect. 13.3 and Table 3 as "the 4 test decks" where it reports the switch-identity
and solver checks. Their inputs and outputs are deliberately not here.

## The two gates, both PASS on all 42 decks

**Release comparison**, from a fresh detached worktree of the release commit with a clean
build, against the sealed reference arm:

```
VERDICT=PASS  exit0=42/42 live=0 bad_input=0
identical=3696 differing=0 missing=0 extra=0
run.log raw-identical 42/42
```

**Runtime checks**, built `-O0 -fcheck=all -Wall`, every deck and not a sample:

```
decks run 42 ; exit 0 42/42 ; 88 output files 42/42
decks with a Fortran runtime error 0
output files byte-identical 3696/3696
VERDICT: PASS
```

## What else was verified, and how

1. **The decks here are the gated bytes.** All 598 input files under `decks/` compared
   against the reference arm the gates measured: 598 compared, 0 differing.
2. **The source here is the gated source.** All 19 `.f90` files and the `Makefile`
   compared against commit `1cab3d557391058e784e8932168aa3dfb9c7d3ea`: 20 compared, 0 differing
   (25 Oct seal; after it, comment text and the 14 printed strings listed under "Re-seal of
   27 October 2026" changed, and nothing else).
3. **Fresh-extraction check.** Extracted to an empty directory, checked against its own
   `MD5SUMS` and `testcase/MD5SUMS`, compiled with the shipped `Makefile` (the binary
   reproduces the release build's md5 exactly), run on the reference case, compared: 88 of
   88 byte-identical, `run.log` byte-identical, `testcase/` unchanged by the run.
4. **Decks run from the extraction**, one of each site class, against the companion
   archive: 88 of 88 output files byte-identical each.

## Public-input comments

The comments of the shipped input files carry no internal development wording. Reaching that
state took four passes, and the first three did not achieve it: they covered only
`driver_main.txt`, so the claim made after them was **wrong** while `driver_minerals.csv`,
`driver_N.csv`, `tsoil_init.csv`, `cultivar_yield.csv` and `phenology_cultivar.txt` still
carried it. The first three passes ran over all 40 shipped `driver_main.txt` copies — the 38
deck drivers, `testcase/` and the `docs/` configuration copy — and made the three uniform in
form:

- records 61-183, removing the internal cross-reference tags and the development commentary;
- **records 1-60**, the original hand-written header block, which the first pass had not
  reached. Record 39 still announced `DEV DEFAULT = 1`, `0=DEPRECATED legacy SW (regression
  reference only)` and a run line for an executable that is not in this release; records 20
  and 21 said "hack for now"; records 52-60 carried internal version tags;
- **record 86**, whose comment described it as reserved and not implemented when it is the
  active grain-nitrogen sink law.

Three accuracy corrections were made in the same passes, each on the evidence of the code:

1. Nine records (23, 24, 26, 27, 31, 35, 36, 39, 50) presented a value as selectable, and
   usually as the default, that the model refuses at startup. A user following the file
   would have had the run stop. Each now states the accepted set the code enforces.
2. Eleven records (2, 6, 9, 14, 15, 21, 28, 34, 41, 42, 48) are read and then used by no
   statement of the release, and seven more (13, 18, 20, 29, 45, 49, 57) are read but have
   their effect blocked by another record or overwritten from an input file. Each now says
   so, and names what blocks it. Inertness was established by experiment, not by reading:
   each record was perturbed, the reference case rerun, and all 88 output files compared —
   byte-identical in every case, while the same test detects an effect for the controls
   (plant density and CO2).
3. Record 86's comment, as above.

Every record position and every value the model parses is unchanged: 183 records in all 40
copies, no value field altered, and record 1 still carries no comment, as it must, because
it is read whole-line without comment stripping. `testcase/input/driver_main.txt` is
byte-identical to the `WDF25_Ruta_aCO2_WW` deck driver, and
`docs/V1_CONFIGURATION_driver_main_v1.0.txt` is byte-identical to that, so the published
configuration, the installation check and the deck agree exactly.

The proof that none of this changed the model is the clean PASS above, measured on every
deck and every output file, and repeated after each of the three passes. Physical notation,
units, numeric ranges and literature citations in those comments are untouched: every DOI
present before the clean-up is still present, and none was invented.

### The fourth pass, 2026-10-26

The fourth pass covered the file types the first three had never opened. It changed **425
files** in this archive and the companion outputs archive, in **2 528** substitutions, and it
replaced five classes of internal development wording with public descriptions:

| what was there | what it says now |
|---|---|
| `THE DEFECT`, and its variants, as the heading of the reason a mode-0 form is not used (1 413 occurrences over 79 files) | `THE LIMITATION ADDRESSED`, `THE LIMITATION OF THAT BRANCH`, and so on — the same statement of the same flaw, without the project's internal word for it |
| `owner ruling <date>`, `ruling D<nn>`, `INPUT FIX`, `ENVELOPE RULE (<internal doc>)`, `OWNER RULING OWED`, lane labels (452 occurrences over 220 files) | the fact the decision recorded — the rule applied, the measurement used, the status of the value — with the decision-process label removed |
| `both are diagnostic attribution modes, not alternative formulations` | `both are diagnostic attribution modes, not alternative formulations` |
| `it moved by more than 10 % between two binaries one day apart` | `it is not transferable between model versions` |
| ~700 pointers to internal design notes, calibration logs and run directories that are **not** part of this release | a named description ending `(internal; not distributed)`, so a reader is told the record exists and is not sent to a path that is not here |

What the pass did **not** touch, by rule: every input value, every record position, every
field count, every citation and page number, every measured quantity, every `PLACEHOLDER`
row and its explanation, every pre-registered prediction, and every `NOT MODELLED,
DECLARED` and `LIMITATION` disclosure. A scientific limitation is a disclosure, not internal
commentary, and none was removed.

In the `.f90` files only comment text was edited — never a string literal — so the compiled
code is unchanged. This was verified, not assumed: the cleaned source compiles to an
executable whose md5 is `f4f6f352708609f3521c6347e06dbbb4`, the sealed release build, bit for
bit. Every deck was then rerun from the cleaned inputs and every output file compared against
the companion reference archive.

**Re-seal of 27 October 2026.** After the sealed build above, two further revisions were made
and the package was re-sealed:

1. The author revised the comments of the 16 `.f90` files that carry them (comment text only;
   verified by a comment-stripped comparison: 0 differing lines of code).
2. Fourteen printed string literals that carried internal wording or pointed to an internal
   document were reworded (`main_MMCO2Crop.f90` lines 3628, 3807, 4031, 4055, 4096, 4128, 5943,
   5995, 8354, 8405, 13769 and 13828; `phenology_module.f90` lines 241 and 321). Thirteen of them
   reach `run.log` only, or sit on error-stop paths that do not fire in the shipped
   configuration; the one at line 13769 is a row label in `output/canopy_sink_summary.csv`.

Because printed text changed, the executable's md5 changed: the re-sealed release build is
md5 `65eb0273a9a543ab807d9f0920072bf1` (same compiler and flags as above). Before and after the string rewording the
two builds were run on the test case and on all 38 decks from identical inputs: in every run
87 of the 88 output files are byte-identical and the 88th (`canopy_sink_summary.csv`) differs
in that one label only; `run.log` differs only in the reworded lines. The reference test case
(`testcase/run.log`, `testcase/output/canopy_sink_summary.csv`, `testcase/MD5SUMS`) and the
companion reference archive (the same two files of every deck, and its `MD5SUMS`) were
regenerated from the re-sealed build on the release machine; all other reference files are
unchanged from the 25 October seal. The WheatDryFACE measurement tables that the frozen run
directories held (`input/WDF25_Data/`, `input/swc_WheatDryFACE_2024_5.csv`; read by no
compiled unit) were removed from the shipped decks and test case; they are the separate
observation deposit cited in `README.md`.

`INPUT_GUIDE.md`, the record-by-record user manual, was corrected at the same time: it had
recorded record 76 as unrestricted, where the code accepts only the value 1. All 124
records that the model guards at startup now agree with the code.

## Integrity

`MD5SUMS` lists every file shipped except itself and `testcase/MD5SUMS`. Check it with
`md5sum -c MD5SUMS` (Linux) or the equivalent on macOS. `testcase/MD5SUMS` covers the
reference data only, so editing `testcase/README.md` cannot make that check fail.
