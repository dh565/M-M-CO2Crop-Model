<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# File inventory — M&M-CO2Wheat v1.0

Written 2026-10-21. Three parts: **what ships** (sources, build, documentation, tools,
test case), **what a run reads** (the input files), and **what a run writes** (the 88
output files, one line each).

The output list was derived from the reference simulation itself — the 88 files in
`testcase/output/`, with their own header rows — not from memory. Units are given as the
file states them; where a column name already carries the unit (`_mm_d`, `_gN_m2`,
`_W_m2`, `_degCd`) the line names the quantity and the basis instead of repeating every
column.

---

## 1. What ships

### 1.1 Sources — 19 files, 18 compiled units

| file | one line |
|---|---|
| `TYPES.f90` | working-precision kind parameters `wp`, `ik`; declaration only |
| `SHARED.f90` | the global store: constants, measured parameters, the 183 driver switches, every pool and ledger, six helper routines |
| `csv_file.f90` | minimal CSV writer (reduced form of Arjen Markus's public-domain module) |
| `csv_file_1d.f90` | `include` body of `csv_file.f90`'s array write routine; not a separate program unit |
| `count_lines.f90` | counts the records of an input file so the daily arrays can be sized |
| `th.f90` | solar declination, sunrise and daylength for the day |
| `solar.f90` | two-leaf PAR partition: extinction coefficients, sunlit and shaded LAI, absorbed PAR per class |
| `meteo.f90` | hourly weather from one daily record: irradiance, air temperature, humidity, VPD, sky emissivity |
| `radiation_module.f90` | component shortwave split and gray-body longwave radiosity; net radiation of canopy and soil |
| `soil_heat_module.f90` | 1-D soil heat conduction on a grid extended to about 5 m; backward Euler, Thomas algorithm |
| `soiln_module.f90` | residue and microbial C-N turnover; mineralization and immobilization |
| `assim.f90` | leaf gas exchange: Farquhar-Busch photosynthesis, Leuning conductance, optional leaf water potential |
| `energy_balance_module.f90` | three-node canopy / soil / canopy-airspace energy balance, damped Newton on a 3x3 Jacobian |
| `hydro.f90` | mixed-form Richards soil water in adaptive implicit sub-steps; van Genuchten-Mualem, Feddes uptake |
| `mineral_module.f90` | six-element grain mineral economy (Fe, Zn, Ca, Cu, Mg, Mn) |
| `phenology_module.f90` | thermal time, vernalization, photoperiod, the eight stage dates and the 0-2 development stage |
| `CANOPY_MODEL.f90` | the hourly canopy integration and the canopy-temperature solve |
| `growth_module.f90` | the daily carbon and nitrogen economy, organ-level canopy, grain number and filling |
| `main_MMCO2Crop.f90` | `PROGRAM mmco2crop`: reads the decks, refuses retired options, runs the day loop, writes the 88 outputs |

Per-module detail, including what calls what and which module consumes each output:
`docs/MODULES.md`.

### 1.2 Build

| file | one line |
|---|---|
| `Makefile` | the only build file; `make` compiles the 19 objects in dependency order and links `mmco2crop` with `-O3 -Wall` |

### 1.3 Documentation

| file | one line |
|---|---|
| `README.md` | the release README: what the model is, dependencies, build, run, verify, the two tolerance statements, source table, citation, licence, contact |
| `INPUT_GUIDE.md` | the input manual: every input file, its format and units, mandatory versus defaulted, and the allowed value of every driver record |
| `LICENSE` | M&M-CO2Crop Research Licence v1.0 |
| `CITATION.cff` | machine-readable citation; version 1.0.0, both companion papers, DOIs marked pending |
| `docs/MODULES.md` | per-module reference: role, calls, called by, outputs and consumers, time step, spatial unit |
| `docs/ARCHITECTURE.dot` | the module call and data-flow graph, the source of truth for the diagram |
| `docs/ARCHITECTURE.png` | the rendered diagram: modules as boxes grouped by time step, arrows labelled with the main variables |
| `docs/FILE_INVENTORY.md` | this file |
| `docs/V1_KNOWN_LIMITATIONS.md` | the optional diagnostics that v1.0 cannot produce, with the record and the file name of each |
| `testcase/README.md` | how to run and compare the Part 1 reference simulation |

### 1.4 Tools

| file | one line |
|---|---|
| `tools/compare_outputs.py` | compares a produced `output/` tree against the shipped expected one; per file, byte identity or the maximum absolute and relative difference per numeric column with the row at which it occurs. Standard library only |
| `scripts/release_acceptance.py` | the release acceptance harness: rebuilds nothing, runs all 42 release decks against a sealed reference arm and reports identical / differing / missing / extra. Standard library only |
| `scripts/render_architecture.py` | renders `docs/ARCHITECTURE.png` from `docs/ARCHITECTURE.dot` without Graphviz; fails if the two disagree. Standard library plus the macOS system thumbnailer |

`scripts/release_acceptance.py` stays in `scripts/`, where every existing evidence file
and the release handoff already point; it was not moved into `tools/`.

### 1.5 Test case

| path | one line |
|---|---|
| `testcase/input/` | the 31 input files of the Part 1 reference simulation `WDF25_Ruta_aCO2_WW`, byte-copied from the frozen v1.0 reference run |
| `testcase/output/` | the 88 expected result files of that simulation, byte-copied from the same arm |
| `testcase/run.log` | the console log of that simulation, produced by the CURRENT build (see `testcase/README.md`) |
| `testcase/MD5SUMS` | 120 digests covering `input/`, `output/` and `run.log` |
| `testcase/README.md` | what the case is, the exact run command, the exact compare command, and the `run.log` note |

### 1.6 Not part of the release set

The working tree also carries development material: the remaining `docs/` and `scripts/`
files, the analysis and figure scripts at the top level, the local literature library, the
input decks other than the reference case, and the historical plan and status documents.
None of it is needed to build, run or verify the model, and none of it is covered by this
inventory.

---

## 2. What a run reads

The model opens its inputs by **relative path**, so it must be run from a directory that
contains an `input/` subdirectory. Formats, units, mandatory-versus-default status and
the driver record table are in `INPUT_GUIDE.md`.

### 2.1 Read by the model

| file | one line |
|---|---|
| `input/driver_main.txt` | the 183 positional run-control records: weather file path, site, canopy, soil, CO2, and every process switch. Record 1 names the weather file |
| `input/input_daily_meteo.csv` | the daily weather forcing named by record 1: year, month, day, sunrise hour, Tmax, Tmin, RHmax, RHmin, wind, precipitation, irrigation, radiation |
| `input/driver_swtr.txt` | soil-profile setup: hydraulic CSV path, Richards sub-step count, initial root depth, water-table flag, layer count, initial volumetric water content per layer, layer depth boundaries |
| `input/SoilLayer_vgm.csv` | per-layer van Genuchten-Mualem parameters; header `Layer,theta_r,theta_s,alpha_m-1,n,l,ksat_m_d` |
| `input/tsoil_init.csv` | initial soil temperature profile (C) |
| `input/phenology_cultivar.txt` | cultivar phenology and early-growth file: thermal-time method, cardinal temperatures, sensitivities, stage thresholds, sowing date, observed stage anchors, bootstrap organ masses and N fractions, and the optional canopy-architecture traits |
| `input/cultivar_yield.csv` | cultivar yield coefficients: fruiting efficiency, potential kernel mass, chaff share, end of ear construction, and the optional kernel-mass envelope and grain-filling clock |
| `input/driver_N.csv` | soil and plant nitrogen parameters, read positionally: partition, reductase, residue and microbial kinetics, hydrolysis and nitrification constants, initial mineral N |
| `input/driver_minerals.csv` | one row per element (Fe, Zn, Ca, Cu, Mg, Mn): buffer power, uptake kinetics, ceilings, phloem mobility, remobilization rate, initial labile pool |
| `input/N_fertilizer.csv` | fertilizer applications: `DAS, N_kg_per_ha[, form]`; form defaults to urea |
| `input/N_organic_treatments.csv` | organic residue and amendment events: count line, then type, DAS, dry mass, N%, C%, fast fraction, age before sowing |

### 2.2 Carried in the reference deck but not read by the model

Byte-copied with the deck so that `testcase/input/` is exactly what the frozen run was
fed. They change nothing.

| file | one line |
|---|---|
| `input/driver_main.template` | the driver template the deck assembler used; not opened by the model |
| `input/driver_main_silledBottom.txt` | an alternative driver with a silled bottom boundary; not opened unless a user renames it |
| `input/SoilLayer_vgm_bottomOnly.csv`, `_restrictL7.csv`, `_restrictL8.csv`, `_retentL8.csv` | alternative soil-hydraulic profiles kept for sensitivity work; only the file named in `driver_swtr.txt` is read |
| `input/driver_resp_mainten.txt` | an orphan input: read by no compiled unit of v1.0, kept in the decks |
| `input/site_observed_phenology.csv` | observed stage dates for the site; reference data for the user, not read by the model |
| `input/site_observed_targets.csv` | observed yield and quality targets; reference data, not read by the model |

---

## 3. What a run writes — the 88 output files

All are comma-separated with a header row, written to `output/`. "Daily" means one row
per simulated day (160 rows in the reference case), "per layer" means one row per day and
soil layer (10 layers), "per element" one row per day and mineral element (6 elements),
"hourly" one row per day and hour (24 per day). Season tables carry one row per axis
class, rank, or named quantity.

### 3.1 Water, energy and radiation

| file | content and units |
|---|---|
| `daily_fluxes.csv` | daily net radiation, latent, sensible and soil heat flux (MJ m-2 d-1) with transpiration, soil evaporation, interception, drainage and runoff (mm d-1) |
| `daily_fluxes_full.csv` | the same with the canopy / soil split of LE and H, precipitation and irrigation, and the season cumulatives (MJ m-2 d-1, mm d-1, mm) |
| `daily_eb_demand.csv` | the day's energy-balance DEMAND terms before the soil-water repartition: Rn, G, LEc, LEs, Hc, Hs (MJ m-2 d-1) |
| `daily_eb_actual.csv` | the committed daily repartition: potential and actual canopy transpiration, interception evaporation and soil evaporation (mm), pond evaporation, and potential versus actual canopy, airspace and soil temperature (C) |
| `daily_resistances.csv` | daily and daylight-mean stomatal conductance (mol m-2 s-1) and the resistances rcs, rca, rsa, rss (s m-1) with the LE and flux terms they produced |
| `daily_radiation_aero.csv` | daily mixed albedo, LAI, canopy height, net radiation and its canopy / soil split (MJ m-2 d-1), frontal area index, displacement height and roughness length (m), and daylight-mean rca and rsa (s m-1) |
| `daily_tleaf_tair.csv` | daily leaf temperature, air temperature, their difference and the committed canopy temperature (C) |
| `daily_tsoil.csv` | daily mean soil surface temperature (C) |
| `daily_tsoil_profile.csv` | daily soil temperature of the ten hydrology layers (C) |
| `daily_soil_heat.csv` | daily integrated surface soil heat flux (J m-2), storage change, heat residual, and surface and top-layer temperature statistics (C) |
| `daily_residue_evap.csv` | surface residue load (t ha-1) and the residue evaporation suppression factor (-) |
| `hourly_energy_balance.csv` | per hour: the solved three-node state (Tair, Tac, Tcan, Tsoil surface, C), vapour pressures (Pa), assimilation, conductance and every resistance (s m-1) |
| `hourly_radiation.csv` | per hour: shortwave down, absorbed and reflected by canopy and soil, the five longwave terms, net radiation per node (W m-2) and the albedos |
| `hourly_soil_heat.csv` | per hour: surface and bottom soil heat flux (W m-2), storage change and residual (J m-2), node temperatures (C), and the soil component fluxes |
| `hourly_esoil_pairing.csv` | per hour: soil evaporation and transpiration sent versus delivered (mm) with the committed latent-heat fluxes (W m-2) — the demand-to-supply pairing check |
| `hourly_node_pairing.csv` | per hour: the canopy temperature armed by the solver versus the one committed after the repartition (C) |
| `diurnal_energy_balance.csv` | per hour: the solver's own diurnal record — leaf temperature, VPD, the flux terms, the energy-balance residual (W m-2), conductance, resistances and assimilation |
| `test_energy_bal_2025.csv` | the compact daily energy-balance diagnostic written through `csv_file`: five unlabelled columns (air and leaf temperature C, and three flux terms) |
| `soil_heat_fallback.csv` | one row: counts of soil-heat fallback events and days, no-root and flux-feedback events, seed resets and hydrology boundary switches |

### 3.2 Soil water

| file | content and units |
|---|---|
| `daily_swc_layers.csv` | daily water storage of each of the ten layers and the profile total (mm) |
| `daily_swc_layers_raw.csv` | the same before any reporting adjustment (mm) |
| `daily_swc_layers_mean.csv` | the day-mean rather than end-of-day storage per layer and total (mm) |
| `daily_swc_water.csv` | daily water supply, total and rooted-zone storage, and plant-available water (mm) |
| `daily_soil_profile_diag.csv` | per day and layer: water content, saturation and residual content (m3 m-3), matric head (m), conductivity and saturated conductivity (m d-1), effective saturation, and the fluxes in and out (m d-1) |
| `daily_hydro_ledger.csv` | the day's water ledger: requested and actual top flux, infiltration, soil and pond evaporation, runoff, potential and actual root uptake, bottom flux, storage (mm d-1) and the season cumulatives (mm) |
| `daily_hydro_diag.csv` | solver diagnostics: failure flag, time-step halvings, minimum-step reached, boundary-condition type, representative fluxes and conductivities (m d-1), ponding (mm), mass-balance terms (m) |
| `wb_components.csv` | the closed daily water balance: effective rain and irrigation, runoff, ponding and its change, storage and its change, actual evaporation and transpiration, ET and the top and bottom fluxes (mm d-1) |
| `substep_diag_all.csv` | one row per Richards sub-step: step length (d), phase, the surface fluxes and caps (m d-1), ponding and infiltration (mm), top-layer state |
| `hydro_wetfront_trace.csv` | one row per sub-step: the wetting-front limiter trace — boundary type, ponding, potential and actual evaporation, requested top flux, top-layer water content, heads and conductivity, the risk flag, the step cap and whether it was applied |
| `daily_root_diagnostics.csv` | daily root depth (m), the water-stress and pressure-response factors, and the water content of the ten layers (m3 m-3) |
| `daily_feddes_node_diag.csv` | the Feddes response: the p2 node in force and the alternative rule, the node classes, the largest difference between the rules, and the per-layer alpha actually used |
| `daily_root_comp_diag.csv` | the Jarvis compensatory-uptake diagnostic: omega, the raw and normalized compensation factors, root uptake and transpiration demand (mm) |
| `daily_root_comp_layers.csv` | per day and layer: root distribution, Feddes alpha, uncompensated and compensated weights and the compensated uptake (mm) |
| `daily_ft_outer_diag.csv` | one row per outer coupling iteration: the water-stress factor in and out, actual and the several potential transpiration terms (mm), the unstressed-sum residual and the convergence flag |
| `daily_spac_diag.csv` | the soil-plant-atmosphere state: leaf water potential (MPa), plant conductance (mol m-2 s-1 MPa-1), the candidate and used soil water potentials (MPa), solver call and iteration counts and the residual |
| `daily_spac_sink_coupling.csv` | the outer SPAC-to-sink loop: iteration count and convergence, demand sent and delivered (mm), the residual, the repartition factor and the relaxation scale |
| `vgm_params.csv` | the van Genuchten-Mualem parameters actually in force per layer: theta_r, theta_s (m3 m-3), alpha (m-1), n, m, l, saturated conductivity (m d-1) |

### 3.3 Phenology, carbon and growth

| file | content and units |
|---|---|
| `daily_phenology.csv` | calendar date and day of year, days after planting, stage name and index, cumulative thermal development units (degCd), the vernalization and photoperiod factors, photoperiod (h) and the pre- and post-anthesis development temperatures (C) |
| `daily_growth.csv` | daily LAI (m2 m-2), plant height (m), the organ dry-matter pools, above-ground biomass and grain yield (g m-2) and the realized specific leaf area (m2 kg-1) |
| `daily_carbon_growth.csv` | the daily carbon budget: net and gross assimilation (g CO2 m-2 d-1), net assimilate as CH2O, maintenance and growth respiration (g CH2O m-2 d-1), the four partition fractions, the organ pools (g m-2), LAI and height |
| `daily_gas_exchange.csv` | daily and daylight-mean canopy gas exchange: net and gross assimilation (g CO2 m-2 d-1), nitrate reduced (g N m-2 d-1), assimilation per ground and per leaf area (umol m-2 s-1), conductance (mol m-2 s-1), transpiration (mm d-1 and mmol m-2 s-1), Ci (umol mol-1), VPD (Pa), leaf and air temperature (C) |
| `flag_leaf_diag.csv` | the flag-leaf diagnostic, per hour: leaf nitrogen per area (g N m-2), assimilation and day respiration (umol m-2 s-1), conductance (mol m-2 s-1), Ci/Ca, leaf temperature (C), leaf-surface VPD (Pa). Diagnostic only; it enters no canopy total |
| `daily_carbon_reserve.csv` | the reserve ledger: pool size, carbon to the grain from current assimilate and from the reserve, grain sink demand and the unfilled part, realized thousand-kernel weight, reserve respired and the ledger residual (g m-2) |
| `daily_stem_wsc.csv` | the stem water-soluble-carbohydrate reserve: pool and stem mass (g m-2), the reserve fraction of the stem, deposition, the amounts refused by the cap and by the window, the draw, the net gain and the closure residual (g m-2 d-1) |
| `daily_wsc_mode2.csv` | the record-105 reserve law: thermal time and time remaining (degCd), the selected window member, deposition, release, release to grain, conversion loss and respiration (g m-2 d-1) |
| `daily_wsc_fate.csv` | the record-148 fate law: reserve before and after (g m-2), the refused and the released terms, the deficit-driven draw and the respired part (g m-2 d-1) |
| `daily_carbon_overflow.csv` | surplus carbon when the sink is full: the day's surplus, what the reserve stored, what was lost, their cumulatives (g m-2) and the sink-full flag |
| `daily_repro_surplus_fate.csv` | the reproductive-phase surplus: cumulative surplus and labelled reserve, carbon refused to the residual pool per day and cumulatively (g m-2), and the number of refused days |
| `daily_root_front.csv` | root front descent: the daily rate (mm d-1), the legacy and the current root depth and the maximum (m) |

### 3.4 Canopy architecture (organ-level)

| file | content and units |
|---|---|
| `daily_canopy_ranks.csv` | daily leaf-rank chronology: thermal age (degCd), Haun stage, ranks appeared and ligulated on both clocks, the local phyllochron (degCd) and the final leaf number |
| `daily_canopy_blade_area.csv` | daily blade area: ranks appearing, ligulated and expanding, potential and realized shoot blade area (cm2 and m2 per shoot) |
| `daily_canopy_tillering.csv` | daily axis population: live and newly born axes (m-2 and per plant), the density of the main stem and five tiller classes, potential and realized blade area index, and why tillering ceased |
| `daily_canopy_senescence.csv` | daily senescence state: Haun stage, green and senescing leaf number per axis class, the phase, green area index, standing dead index and pool, and the total blade index (m2 m-2) |
| `daily_canopy_sheath_ear.csv` | daily sheath, exposed internode and ear-plus-peduncle area indices with their alternative lags and bounds (m2 m-2) |
| `daily_canopy_sink.csv` | the lamina sink ledger: the green-lamina LAI, the LI-2200C-basis plant area index, the nitrogen-capacity-only LAI (m2 m-2), area built and its daily increment, lamina demand, allocation and backlog (g DM m-2), the binding limit and the realized specific leaf area |
| `daily_canopy_n_profile.csv` | daily canopy nitrogen profile: LAI, flag-leaf green area fraction, leaf nitrogen (g N m-2), mean and top-of-canopy nitrogen per area, the top-to-mean ratio, and the absorbed-PAR consistency check |
| `daily_d2_organs.csv` | the per-axis supply and demand ledger of the organ-level canopy: previous-day surplus carbon, supply and demand per tiller class (g DM), their ratios and the class densities (m-2) |
| `canopy_rank_chronology.csv` | season table, one row per leaf rank: appearance and ligulation thermal age and calendar day on both clocks, flag-leaf flag, final leaf number and the reference phyllochron (degCd) |
| `canopy_blade_area_ranks.csv` | season table, one row per rank: blade length and width (cm), form factor, potential and realized area (cm2), appearance and ligulation ages (degCd), and the potential area on each tiller class |
| `canopy_sheath_ear_ranks.csv` | season table, one row per rank: sheath length and diameter (cm), sheath silhouette area (cm2), internode length, sheath-covered and exposed internode length (cm) |
| `canopy_tiller_classes.csv` | season table, one row per axis class: measured emergence Haun stage, birth thermal age (degCd), the interval window, occurrence probability, born and final density (m-2), leaf number and the potential, realized and sterile axis areas (cm2) |
| `canopy_senescence_axes.csv` | season table, one row per axis class: leaf number, Haun stage, green and senescing leaf counts, phase, the three senescence onset ages (degCd), green and dead area (cm2) and class density (m-2) |
| `canopy_sink_summary.csv` | season table of named quantities with their units: final leaf number cut, cohorts aborted, days of suspended birth, the partition-sum error, the realized specific leaf area and the three end-of-season leaf area indices, plus the declared inputs and assumptions for an instrument forward model |
| `d2_cohort_table.csv` | season table, one row per cohort: birth age (degCd), window length, occurrence probability, window surplus carbon and cohort demand (g DM), the supply ratio, born and final density (m-2) and leaf number |

### 3.5 Yield formation

| file | content and units |
|---|---|
| `daily_yield_components.csv` | daily grain, chaff and whole-spike dry matter (g m-2), kernel number (m-2), thousand-kernel weight (g), grain yield (g m-2) and harvest index |
| `daily_kernel_potential.csv` | the kernel-potential law: kernel number (m-2), grain mass and grain nitrogen (g m-2), nitrogen per kernel (mg N), the setting thermal time and cumulative supply, the potential kernel mass (mg) and the nitrogen-ceiling terms (g N m-2) |
| `grain_number_chain.csv` | season table, one row per axis class: class density (m-2), spike dry weight per axis (g), competent florets and kernels per spike, the spike dry-matter pool (g m-2), floret survival and the ear fraction |
| `daily_grain_n_sink.csv` | the Martre grain-nitrogen source-sink, one row per filling day: post-anthesis thermal time (degCd), the per-kernel potential, cumulative potential, daily sink and supply, what bound it, and the import split between current assimilate, the mobile pool, lamina, stem and chaff (g N m-2) |
| `daily_grain_n_floor.csv` | the minimum grain-nitrogen concentration law: grain mass and nitrogen (g m-2), grain N concentration (%), the dry-matter room each ceiling allows, the transfer before and after the floor, the refused dry matter per day and cumulatively (g m-2) and the binding flag |
| `grain_filling_law.csv` | season table of the grain-filling law's own constants and the value each took, with units or source |
| `grain_n_sink_law.csv` | season table of the grain-nitrogen sink law's constants and values, with units or source |
| `stem_wsc_law.csv` | season table of the stem reserve law's constants and values, with units or source |

### 3.6 Nitrogen

| file | content and units |
|---|---|
| `daily_nitrogen.csv` | the daily soil and plant nitrogen balance: mineral soil nitrogen, fertilizer, mineralization, leaching, uptake, gross and retained biochemical assimilation, demand, the transpiration and CO2 factors and the remobilization terms (g N m-2, g N m-2 d-1) |
| `daily_soilN_layers.csv` | per day: nitrate and ammonium in each of the ten layers (g N m-2), bottom leaching (g N m-2 d-1), bottom water flux (mm d-1) and the advection sub-step count |
| `daily_organic_N.csv` | the residue and microbial pools (g C m-2, g N m-2) with the day's decomposition, respiration, net mineralization, immobilization from ammonium and nitrate, microbial turnover and the soil-nitrogen closure |
| `daily_n_uptake_terms.csv` | per day and layer: mass-flow, diffusive and active uptake, the cap and the realized uptake (g N m-2), the uptake efficiency, mobile-N concentration and the demand before and after |
| `daily_n_payment.csv` | the construction-nitrogen payment: the three organ concentrations (g N g-1 DM), the mobile and assimilate pools, luxury nitrogen available and used, the organ demands and what each was paid (g N m-2) |
| `daily_n_construction.csv` | the construction-nitrogen status law: crop mass (t ha-1), the critical concentration and the status factor, and per organ the floor, target, effective target, actual concentration and the two payment passes (g N m-2) |
| `daily_leafN_capacity.csv` | the lamina nitrogen capacity route: assimilating LAI, green mass, leaf nitrogen before and after, the exported and surplus nitrogen, nitrogen per assimilating area (g N m-2) and leaf nitrogen as percent of dry matter |
| `daily_lai_n_arrest.csv` | the nitrogen arrest of leaf expansion: the mode and whether it is live, the priority fraction, the area increment, area removed and shed with their cumulatives (m2 m-2), the dry matter removed and the specific leaf nitrogen before and after |
| `daily_organ_n_turnover.csv` | the NEMA organ-nitrogen turnover, one row per post-anthesis day: the mobile pool before, the uptake credit, the demand and supply of lamina, stem, chaff and root, the grain import and what bound it (g N m-2) |
| `daily_itemF_nclose.csv` | the nitrogen closure check: assimilating LAI, leaf nitrogen, the specific-leaf-nitrogen cap and the mean, the synthesis and mobile clipping per day and cumulatively, and the post-anthesis grain demand (g N m-2) |
| `daily_tiller_n_afford.csv` | the nitrogen term of the tiller-survival ledger: the priority fraction, the lamina concentration, the days it bound, and per axis class the density, nitrogen demand and supply and the three affordability factors |
| `daily_tiller_afford4.csv` | the carbon tiller-survival ledger (record 99 mode 4): the previous day's supply, leaf and floret totals, the identity check, cumulative deaths, the target densities and, per axis class, density, Haun stage, survival and death terms |

### 3.7 Grain minerals

| file | content and units |
|---|---|
| `daily_minerals.csv` | per day and element: mass-flow and active uptake, remobilization to the grain, the labile soil pool and the root, leaf, stem, chaff and grain pools (mg m-2), the balance residual, the grain concentration (mg kg-1) and the cap flag |
| `daily_mineral_supply.csv` | per day and element: development stage, transpiration (mm d-1), soil-solution concentration (mg L-1), the uptake terms and total (mg m-2), the labile pool, shoot content (mg m-2) and shoot concentration (mg kg-1) |
| `daily_mineral_uptake_terms.csv` | per day and element: the mass-flow, diffusive and high-affinity transport terms, the cap and the realized uptake (mg m-2) |
| `mineral_sink_cap.csv` | season table, one row per element: the grain ceiling (mg kg-1), the cumulative amount refused by it (mg m-2), the days it bound, and the final grain pool (mg m-2) and concentration (mg kg-1) |
