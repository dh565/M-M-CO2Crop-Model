<!-- This file (its text) (c) 2018-2026 David Helman and The Hebrew University of Jerusalem, under the M&M-CO2Crop Research Licence v1.0 (see LICENSE). -->
<!-- The DATA TABLES in this folder are NOT under that licence: they are third-party open data redistributed under the source licences stated below (CC BY 4.0 for Braunschweig; CC0 1.0 for Maricopa). See also the licence scope in the package README.md. -->

# Observation tables distributed with the scripts

Every file here is derived from an openly published data set and is redistributed so that the
Braunschweig and Maricopa parts of every figure and table run without any further download.
The WheatDryFACE observations are **not** here; they come from the open WheatDryFACE deposit (see ../README.md).

## Licence of these tables, and what was changed

**These tables are not covered by the M&M-CO2Crop Research Licence.** Each carries the licence of
its own source, and in both cases that licence permits the redistribution done here:

* **Braunschweig** (the 16 `braunschweig/BRAUNSCHWEIG_*_obs.csv` files, and the Braunschweig rows
  of `obs_registry_braunschweig_maricopa.csv` and `consolidated_obs_braunschweig_maricopa.csv`)
  --- (c) Manderscheid, R., Dier, M., Erbs, M., Luig, A., Oldenburg, E., Sickora, J., and
  Weigel, H.-J. (2020), used and redistributed under the **Creative Commons Attribution 4.0
  International Licence (CC BY 4.0)**, https://creativecommons.org/licenses/by/4.0/.
  **Modified:** these files are not the published workbooks. They are aggregations and
  reshapings of them --- treatment means with the sample standard deviation and the number of
  replicate rings, in long format by sheet, date and variable. No value was altered other than by
  that aggregation, and no plot-level record is reproduced.
* **Maricopa** (the 6 `maricopa/*.csv` files, and the Maricopa rows of the two tables above) ---
  derived from the Harvard Dataverse deposit doi 10.7910/DVN/YRTQ3Z of Kimball, B. A.,
  Pinter Jr., P. J., LaMorte, R. L., Leavitt, S. W., Hunsaker, D. J., Wall, G. W., Wechsung, F.,
  Wechsung, G., Bloom, A. J., and White, J. W., which is released under the **CC0 1.0 Universal
  Public Domain Dedication**, https://creativecommons.org/publicdomain/zero/1.0/. CC0 imposes no
  conditions on reuse; the citation below is given as a scholarly obligation. **Modified** in the
  same sense: aggregated and reshaped, as the table states per file.

The two `*_braunschweig_maricopa.csv` tables mix both sources and therefore follow the stricter
of the two, CC BY 4.0.

Anyone redistributing the Braunschweig tables further must keep this notice, the attribution and
the statement that the data were modified.

| file | content | derived from |
|---|---|---|
| `braunschweig/BRAUNSCHWEIG_<trno>_obs.csv` (16 files) | long-format treatment means and sample sd over the ring replicates per sheet, date and variable (organ dry matter, LAI/LGrAI, N concentrations, soil moisture 0-20 and 20-40 cm, stem WSC, dated phenology rows) | Manderscheid, R., Dier, M., Erbs, M., Luig, A., Oldenburg, E., Sickora, J., and Weigel, H.-J. (2020). Open Data Journal for Agricultural Research 6, doi 10.18174/odjar.v6i0.16397 (treatment numbers 101-122 of that data set) |
| `maricopa/MARICOPA_YECORA_ROJO_MEASURED.csv` | harvest end points per treatment (HWUM, H#UM, P#AM, CWAM), sheet `Obs_summary_Avg_over_Reps` | Kimball, B. A., et al. (2017). Open Data Journal for Agricultural Research 3, 29-38, doi 10.18174/odjar.v3i1.15826; data doi 10.7910/DVN/YRTQ3Z, workbook 'Biomass Yield Area Phenology Management Weather Soil Moisture.ods' |
| `maricopa/MARICOPA_YECORA_ROJO_REPS_SD.csv` | replicate means and sd of HWAH, CWAM, LAIX, sheet `Obs_summary_by_Reps` | same workbook |
| `maricopa/MARICOPA_SOIL_MOISTURE_AVG_OVER_REPS.csv` | volumetric soil moisture by depth span, treatment means, sheet `Soil_moisture_Avg_over_Reps` | same workbook |
| `maricopa/maricopa_archive_summary.csv` | per-treatment summary (anthesis and maturity dates, HWAH, CWAM, HWUM, H#UM, P#AM, HIAM, LAIX with replicate sd) | same workbook |
| `maricopa/maricopa_archive_daily.csv` | dated LAI (LAID), above-ground (CWAD), grain (GWAD), leaf (LWAD) and stem (SWAD) dry matter: replicate mean, sd and count, sheet `Obs_daily` | same workbook |
| `maricopa/maricopa_eb_hourly.csv` | hourly canopy temperature, air temperature, net radiation, sensible, latent and soil heat flux per modelled treatment (906, 908, 910, 912-916), replicate plots averaged with their sd; 15-min records averaged to the hour ending | the three 'Weather Canopy & Soil Temperatures Energy Balance ...' workbooks of the same archive (Kimball et al. 1999, Water Resour. Res. 35:1179; 2017); rebuilt by `lib/extract_maricopa_eb.py` |
| `obs_registry_braunschweig_maricopa.csv` | the per-treatment observed harvest and canopy scalars (obs, obs_sd, basis) registered for the Part 2 evaluation | the Braunschweig and Maricopa files above |
| `consolidated_obs_braunschweig_maricopa.csv` | every observation of the consolidated observation-model table of the Part 2 evaluation for the Braunschweig and Maricopa decks (obs, obs_sd, obs_source); the figures read the rows the loaders do not carry (Braunschweig grain Fe and Zn and the 2014/15 harvest index, Maricopa peak LAI and harvested yield) | the files above |

Units are stated in the header comment of each file (lines beginning with `#`).
