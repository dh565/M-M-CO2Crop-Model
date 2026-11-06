<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Figure and table map

For every figure and table of Part 1, Part 2 and Supplement S2: the script and function that
produces it, the model files it reads (per deck of the reference-outputs archive, `REF_OUTPUTS/<deck>/...`)
and the observation files it reads. Paths in `observations/` are redistributed with the scripts;
paths beginning with `OBS_DIR/` must be supplied by the reader (the WheatDryFACE processed files
from the open WheatDryFACE deposit, the Maricopa workbooks from the archive of Kimball et al., 2017).

## Model files read by every Part 2 figure and table (loaded once by `lib/paperfigs.init`)

For each of the 38 decks: `output/daily_growth.csv`, `daily_yield_components.csv`, `daily_nitrogen.csv`,
`daily_phenology.csv`, `daily_swc_layers.csv`, `daily_swc_layers_mean.csv`, `daily_stem_wsc.csv`,
`daily_carbon_reserve.csv`, `flag_leaf_diag.csv`, `daily_minerals.csv`, `stem_wsc_law.csv`,
`grain_filling_law.csv`, `grain_number_chain.csv`; `input/driver_swtr.txt` (soil layer boundaries);
`input/phenology_cultivar.txt` (Maricopa: sowing date and observed stage DAP).

Observation files read by `init` (observation registry, phenology, soil water series, grain minerals,
Maricopa grain protein): `observations/obs_registry_braunschweig_maricopa.csv`,
`observations/braunschweig/BRAUNSCHWEIG_<trno>_obs.csv` (16), `observations/maricopa/MARICOPA_YECORA_ROJO_MEASURED.csv`,
`observations/maricopa/MARICOPA_YECORA_ROJO_REPS_SD.csv`, `observations/maricopa/maricopa_archive_summary.csv`,
`observations/consolidated_obs_braunschweig_maricopa.csv`;
**WheatDryFACE deposit:** `OBS_DIR/WheatDryFACE/obs_registry_wheatdryface.csv`, `OBS_DIR/WheatDryFACE/WDF25/S3_<genotype>_<CO2>_<water>.csv` (8),
`OBS_DIR/WheatDryFACE/WDF23/WDF23_<genotype>_<CO2>.csv` (4); Maricopa archive:
`OBS_DIR/Maricopa_FACE_archive/Quality Wheat Grain Baking 1993-1997.ods`.

## Part 2 (`make_part2_figures.py`, `make_part2_tables.py`)

| item | file stem | function | model files beyond `init` | observation files beyond `init` |
|---|---|---|---|---|
| Fig. 1 energy balance, both sites | `FigN01_energy_balance_all_sites` | `lib/paperfigs2.fig_eb` | `output/hourly_radiation.csv`, `hourly_soil_heat.csv`, `input/input_daily_meteo.csv` (WDF decks); `output/hourly_energy_balance.csv`, `daily_phenology.csv` (Maricopa 906, 908, 910, 912–916) | WDF treatment files (radiation Day/Night blocks, soil heat flux); `observations/maricopa/maricopa_eb_hourly.csv` |
| Fig. 2 soil water and temperature | `Fig02_soil_water_temperature` | `lib/paperfigs.fig02` | `output/daily_swc_layers_mean.csv`, `daily_tsoil_profile.csv` | WDF treatment files (SWC 20/55 cm, Tsoil 20/55 cm); `observations/braunschweig/*_obs.csv` (d0_20, d20_40); `observations/maricopa/MARICOPA_SOIL_MOISTURE_AVG_OVER_REPS.csv`; Maricopa archive: the three `Weather Canopy & Soil Temperatures Energy Balance ... .ods` (soil temperature profiles) |
| Fig. 3 flag-leaf gas exchange | `Fig03_flag_leaf_gas_exchange_all_sites` | `lib/paperfigs2.fig_gx` | `output/flag_leaf_diag.csv`, `input/driver_main.txt` (longitude, record 4) | WDF treatment files (An, gs, Tr); Maricopa archive: `Physiology Gas Exchange Water Relations.ods` |
| Fig. 4 leaf area | `FigN04_leaf_area_all_sites` | `lib/paperfigs.figN_canopy` | `output/daily_growth.csv` | WDF treatment files (LAI); Braunschweig `LGrAI`; `observations/maricopa/maricopa_archive_daily.csv` (LAID) |
| Fig. 5 biomass | `FigN05_biomass_all_sites` | `lib/paperfigs.figN_biomass` | `output/daily_growth.csv` | WDF harvest biomass; Braunschweig `SUMDW`; `maricopa_archive_daily.csv` (CWAD) |
| Fig. 6 harvest variables | `FigN06_harvest_all_sites` | `lib/paperfigs.figN_harvest` | `init` scalars | registry, consolidated table, `maricopa_archive_summary.csv` (HWAH, HIAM) |
| Fig. 7 CO2 responses, agronomy | `FigN07_eCO2_response_agronomy_all_sites` | `lib/paperfigs.figN_eco2_agro` | as Figs. 3 and 6 | as Figs. 3 and 6 |
| Fig. 8 CO2 responses, quality | `FigN08_eCO2_response_quality_all_sites` | `lib/paperfigs.figN_eco2_quality` | as Fig. 6 | as Fig. 6 (Braunschweig grain Fe/Zn from the consolidated table) |
| Fig. 9 water and N responses | `FigN09_water_and_N_responses_all_sites` | `lib/paperfigs.figN_water_N` | as Figs. 3 and 6 | as Figs. 3 and 6 |
| Fig. 10 canopy temperature | `FigN10_canopy_temperature_all_sites` | `lib/paperfigs2.fig_tcan` | `output/hourly_energy_balance.csv`, `daily_phenology.csv` | **WheatDryFACE deposit:** `OBS_DIR/WheatDryFACE/WDF25/S3_*.csv` (Tleaf_flir, DTleaf_air_flir); `observations/maricopa/maricopa_eb_hourly.csv` |
| Fig. 11 sign agreement | `FigN11_sign_agreement_site_two_bars` | `lib/paperfigs.figN_signs_site2` | `init` scalars and minerals | as Fig. 6 |
| Table 4 agreement by site and variable; text numbers | `table4_agreement_by_site_variable.csv`, `table4_series_statistics.csv`, `part2_numbers.json` | `make_part2_tables.numbers` | all of the above | all of the above |
| Table 5 numerical checks | `table5_numerical_checks.csv/.json` | `make_part2_tables.table5` | `output/daily_eb_actual.csv`, `daily_nitrogen.csv`, `daily_minerals.csv`, `run.log` (all 38 decks) | none |

## Supplement S2 (`make_supplement_s2.py`)

| item | file stem | function | model files beyond `init` | observation files beyond `init` |
|---|---|---|---|---|
| Figs. S1–S4 soil water, every treatment | `FigS1_soil_water_WDF23` … `FigS4_soil_water_MARICOPA` | `make_supplement_s2.fig_soil_water` (via `lib/paperfigs._swc_course`) | `output/daily_swc_layers_mean.csv` | as Fig. 2 (soil water) |
| Fig. S5 harvest variables by treatment | `FigS5_harvest_by_treatment` | `make_supplement_s2.fig_s5` | `init` scalars | as Fig. 6 |
| Fig. S6 genotype contrast | `FigS6_genotype_contrast_WDF` | `lib/paperfigs.figN_genotype` | as Figs. 3 and 6 | as Figs. 3 and 6 (WheatDryFACE only) |
| Fig. S7 Maricopa soil temperature | `FigS7_soil_temperature_MARICOPA` | `lib/maricopa_tsoil.figS7` | `output/daily_tsoil_profile.csv`, `input/driver_swtr.txt` | Maricopa archive: the three `Weather Canopy & Soil Temperatures Energy Balance ... .ods` |
| Fig. S8 porometer | `FigS8_porometer_WDF25` (+ `FigS8_porometer_pairs.csv`) | `lib/wdf_poro.figS8` | `output/flag_leaf_diag.csv`, `hourly_energy_balance.csv`, `input/driver_main.txt` | **WheatDryFACE deposit:** `OBS_DIR/WheatDryFACE/WDF25/S3_*.csv` (gs_poro, Tr_poro, Tleaf_poro, DTleaf_air_poro; Tleaf_flir; LI-6800 An, gs, Tr) |
| Tables S1–S4 simulated / observed by treatment | `TableS1_WDF23.csv` … `TableS4_MARICOPA.csv` | `make_supplement_s2.tables_s1_s4` | `init` scalars | as Fig. 6 |
| Table S5 closures per simulation | `TableS5_closures_by_simulation.csv` | `make_supplement_s2.table_s5` | `run.log` of every deck | none |

## Part 1 (`make_part1_figures.py`, `make_part1_tables.py`, `run_sensitivity.py`)

| item | file stem | function | inputs |
|---|---|---|---|
| Fig. 1 model structure | `P1_Fig01_model_structure` | `lib/p1_figures.fig1` | none (schematic) |
| Fig. 2 driver responses | `P1_Fig02_driver_responses` | `lib/p1_figures.fig2` | `SENS_DIR/all_metrics.json` |
| Fig. 3 CO2 response and elasticities | `P1_Fig03_alt_cascade_elasticity` | `lib/p1_figures.fig3` | `SENS_DIR/all_metrics.json` |
| Fig. 4 CO2 x water, N, temperature | `P1_Fig04_CO2_interactions` | `lib/p1_figures.fig4` | `SENS_DIR/all_metrics.json` |
| Fig. 5 parameter sensitivity | `P1_Fig05_alt_matrix_tornado` | `lib/p1_figures.fig5` | `SENS_DIR/all_metrics.json` |
| Table 5 sensitivity indices; Sect. 13.3/14 numbers | `table5_parameter_sensitivity.csv`, `part1_numbers.txt` | `make_part1_tables.report` | `SENS_DIR/all_metrics.json` |
| `all_metrics.json` | — | `run_sensitivity.py` | `REF_OUTPUTS/WDF25_Ruta_aCO2_WW/input/` (control deck: `driver_main.txt` records 7, 8, 16, 22, 28; `driver_N.csv`; `cultivar_yield.csv`; `input_daily_meteo.csv`; `N_fertilizer.csv`; `phenology_cultivar.txt`), the model executable `MODEL_EXE`; per run `output/daily_growth.csv`, `daily_yield_components.csv`, `daily_fluxes.csv`, `daily_gas_exchange.csv`, `daily_nitrogen.csv`, `daily_eb_actual.csv`, `daily_phenology.csv`, `hourly_energy_balance.csv`, `daily_minerals.csv` |

## WheatDryFACE files the scripts read (not shipped here; for the open deposit)

All of them are processed per-treatment files; none is a plot-level raw record. "Treatment" =
genotype x CO2 (x water in 2024/25); each value is the mean over the replicate rings of that
treatment with its standard deviation in the adjacent `Std` column.

| expected path under `OBS_DIR` | granularity of what is read | columns read | used by |
|---|---|---|---|
| `WheatDryFACE/WDF25/S3_{Ruta,Zahir}_{aCO2,eCO2}_{WW,WL}.csv` | per treatment: dated 'Phenological stage' column; daily ring-mean sensor series (SWC 20/55 cm, Tsoil 20/55 cm); daily Day/Night block means of the radiation components; campaign means (An, gs, Tr, LAI, Height; thermal-image canopy temperature and canopy-minus-air difference; porometer gs, Tr, leaf temperature and leaf-minus-air difference); harvest means (organ Biomass, [N], grain Fe, Mg, Ca, Cu, Mn, Zn). **Fig. 1 reads the daily energy-balance values (Rn, LWin, LWout, SWin, SWout, albedo, soil heat flux, Day and Night blocks) of every treatment file: the radiometers of a treatment are one or two rings (Zahir treatments one ring), so these daily values are ring-level (per-plot) daily means for those treatments, and the paper counts them as plot-days.** | `Phenological stage`; `An`, `gs`, `Tr`, `LAI`, `Height`; `SWC 20cm`, `SWC 55cm`, `Tsoil 20cm`, `Tsoil55 cm`; Night: `Rn`, `LWin`, `LWout`, `Soil heat flux`; Day: `albedo`, `Rn`, `LWin`, `LWout`, `SWin`, `SWout`, `Soil heat flux` (soil heat flux in the Ruta files only); `Biomass` (Grain, Stem, Spike (head), Leaf, Chaff); `[N]` (Grain, Stem, Leaf, Chaff); `Fe`, `Mg`, `Ca`, `Cu`, `Mn`, `Zn` (Grain); each with `Std`. `Tleaf_flir`, `DTleaf_air_flir`; `gs_poro`, `Tr_poro`, `Tleaf_poro`, `DTleaf_air_poro`; each with `std` | every Part 2 figure and table; Supplement Figs. S2, S5, S6, S8, Table S2; thermal-image columns: Part 2 Fig. 10, Table 4 (canopy-temperature rows), Supplement Fig. S8; porometer columns: Supplement Fig. S8 |
| `WheatDryFACE/WDF23/WDF23_{Ruta,Zahir}_{aCO2,eCO2}.csv` | per treatment, as the 2024/25 files for the 2022/23 season (no soil heat flux, thermal-image or porometer columns; harvest block: Straw, Grain, Spike (head), Chaff) | as the 2024/25 files, except: `Biomass (leaf + stem) (t/ha)` and `[N]` (Straw); `Biomass`, `[N]`, `Fe`, `Mg`, `Ca`, `Cu`, `Mn`, `Zn` (Grain); `Biomass` (Spike (head)); `[N]` (Chaff) | every Part 2 figure and table; Supplement Figs. S1, S5, S6, Table S1 |
| `WheatDryFACE/obs_registry_wheatdryface.csv` | per treatment: harvest and canopy scalars (obs, obs_sd) derived from the files above | `deck, site, quantity, obs, obs_sd, basis` | every Part 2 figure and table; Supplement Figs. S5, S6, Tables S1-S2 |

## Maricopa archive workbooks the scripts read (open data, Kimball et al. 2017; not shipped)

| expected path under `OBS_DIR` | content read | used by |
|---|---|---|
| `Maricopa_FACE_archive/Weather Canopy & Soil Temperatures Energy Balance 1994 hourly.ods`, `... CO2 1996 15-min.ods`, `... CO2 O3 1997 15-min.ods` | hourly / 15-min soil temperature profiles per instrumented plot (Fig. 2, Fig. S7); the same workbooks are the source of the shipped `observations/maricopa/maricopa_eb_hourly.csv` | Part 2 Fig. 2; Supplement Fig. S7 |
| `Maricopa_FACE_archive/Physiology Gas Exchange Water Relations.ods` | sheet `FACE_93-97_data`: leaf gas exchange per plot, date and hour (averaged per replicate plot, then over plots) | Part 2 Figs. 3, 7, 9, Table 4 |
| `Maricopa_FACE_archive/Quality Wheat Grain Baking 1993-1997.ods` | whole-grain protein and moisture per plot (averaged over plots) | Part 2 Figs. 6, 8, 11, Table 4; Supplement Fig. S5, Table S4 |
