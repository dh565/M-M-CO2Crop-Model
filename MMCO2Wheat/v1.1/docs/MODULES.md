<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Module reference — M&M-CO2Wheat v1.0

Written 2026-10-21. One section per compiled source file: what it does, what it calls,
what calls it, what it produces and which unit consumes each product, its time step and
its spatial unit. The content is taken from the release header of each file and was
verified against the code of this commit: every `USE` dependency below was re-extracted
from the sources, and the external-subroutine calls were checked at their call sites.

Nineteen source files are compiled. `TYPES.f90` and `SHARED.f90` carry declarations only;
`csv_file_1d.f90` is an `include` body, not a program unit. The build order is the `OBJS`
line of the `Makefile` and is fixed by these dependencies: `TYPES` first, `SHARED` second,
everything else above them, `main_MMCO2Crop` last.

## Dependency summary

| file | `USE`s | called by | time step | spatial unit |
|---|---|---|---|---|
| `TYPES.f90` | — | all 19 others (kind parameters) | none | none |
| `SHARED.f90` | `TYPES` | 15 of the 20 units | all (store) | all (store) |
| `csv_file.f90` | — | `main_MMCO2Crop` | per write call | none |
| `csv_file_1d.f90` | (include body of `csv_file`) | `csv_file` (textual include) | per write call | none |
| `count_lines.f90` | `SHARED` | `main_MMCO2Crop` | once per file, at init | none |
| `th.f90` | `SHARED`, `TYPES` | `CANOPY_MODEL` | daily | point |
| `solar.f90` | `SHARED`, `TYPES` | `CANOPY_MODEL` | hourly | point, two leaf classes |
| `meteo.f90` | `SHARED`, `TYPES` | `CANOPY_MODEL` | hourly, from a daily record | point, reference height |
| `radiation_module.f90` | `SHARED`, `TYPES` | `energy_balance_module` | hourly, per energy-balance iterate | canopy layer, soil surface |
| `soil_heat_module.f90` | `SHARED`, `TYPES` | `main_MMCO2Crop`, `energy_balance_module` | energy-balance sub-step (hourly) | soil layer, heat grid |
| `soiln_module.f90` | `SHARED`, `TYPES` | `main_MMCO2Crop` | daily | bulk soil stocks |
| `assim.f90` | `SHARED`, `TYPES` | `CANOPY_MODEL` | hourly, per leaf class | one leaf |
| `energy_balance_module.f90` | `SHARED`, `TYPES`, `radiation_module`, `soil_heat_module` | `CANOPY_MODEL`, `main_MMCO2Crop` | hourly, per iterate | canopy, soil surface, canopy airspace |
| `hydro.f90` | `SHARED`, `TYPES`, `IEEE_ARITHMETIC` | `main_MMCO2Crop` | adaptive Richards sub-step within a day | 10 soil layers to 1.20 m + surface store |
| `mineral_module.f90` | `SHARED`, `TYPES` | `main_MMCO2Crop` | daily, per element | one labile soil pool (0-0.6 m) + organ pools |
| `phenology_module.f90` | `TYPES` | `main_MMCO2Crop` | daily | whole plant (the apex) |
| `CANOPY_MODEL.f90` | `SHARED`, `TYPES`, `energy_balance_module` (+ external `assim`, `meteo`, `solar`, `th`, `et`) | `main_MMCO2Crop` | hourly, 24 per day, with a nested fixed point | canopy as two leaf classes |
| `growth_module.f90` | `SHARED`, `TYPES` | `main_MMCO2Crop` | daily | plant organ, tiller axis, leaf rank |
| `main_MMCO2Crop.f90` | `SHARED`, `TYPES`, `csv_file`, `CANOPY_MODEL`, `energy_balance_module`, `soil_heat_module`, `growth_module`, `mineral_module`, `soiln_module`, `phenology_module` (+ external `hydro`, `count_lines`) | nothing — it is the main program | all four | all |

`CANOPY_MODEL` reaches `assim`, `meteo`, `solar`, `th` and `et` as external subroutines
through explicit interfaces rather than through `USE`; `main_MMCO2Crop` reaches `hydro`
the same way. That is why those five names do not appear in the `USE` column of their
callers above.

---

## 1. `TYPES.f90` (32 lines)

**Role.** Declares the working-precision kind parameters: `wp` for double-precision reals
(`SELECTED_REAL_KIND(15, 307)`) and `ik` for integers (`SELECTED_INT_KIND(9)`).
Declaration-only; no state, no executable statement.

**Calls.** Nothing (the intrinsics `SELECTED_REAL_KIND` / `SELECTED_INT_KIND` only).

**Called by.** Every other compiled unit, for the kind parameters only.

**Outputs and their consumers.** `wp`, `ik` — every model variable in `SHARED`,
`main_MMCO2Crop`, `assim`, `CANOPY_MODEL`, `growth_module`, `phenology_module`, `hydro`,
`et`, `meteo`, `radiation_module`, `soil_heat_module`, `soiln_module`, `mineral_module`,
`energy_balance_module`, `solar` and `th` is declared with them.

**Time step / spatial unit.** None.

## 2. `SHARED.f90` (4 388 lines)

**Role.** The model's single global store of parameters and state: the physical constants,
the measured plant, soil, canopy and soil-nitrogen parameters with their sources, the
driver-record switches of `input/driver_main.txt`, every prognostic pool, every diagnostic
ledger, and six small helper routines (state synchronization, hourly mass binning,
Gauss-weighted hour-to-day summation, driver-comment stripping, an event-reachability
check). No process law is solved here.

**Calls.** Nothing outside itself.

**Called by.** `USE SHARED` appears in 14 of the 18 compiled units. The four that do not
use it are `TYPES.f90`, `csv_file.f90` and `phenology_module.f90` — and `SHARED` itself.
(`csv_file_1d.f90` is the `include` body of `csv_file.f90`, not a compiled unit of its
own, so it is not one of the 18.) The six contained routines are called from `main_MMCO2Crop.f90`
(`sync_hydro_step_state`, `strip_driver_comment`), `CANOPY_MODEL.f90` (the hourly binning
and the Gauss summation) and `hydro.f90`.

**Outputs and their consumers.** The module *is* the state. It publishes the constants,
the per-record configuration switches, the soil hydraulic and thermal fields, the canopy,
organ, nitrogen and mineral pools, the energy-balance and Richards ledgers, and the output
unit numbers. Consumers: `hydro`, `et`, `th` (soil water, evaporation, daylength);
`energy_balance_module`, `radiation_module`, `soil_heat_module` (surface energy balance,
soil heat); `solar`, `meteo`, `CANOPY_MODEL` (radiation, weather disaggregation, hourly
canopy); `assim` (leaf gas exchange); `growth_module`, `soiln_module`, `mineral_module`
(crop C and N, residue N, grain minerals); `count_lines` (input sizing); and
`main_MMCO2Crop`, which fills it from the decks and writes it out.

**Time step / spatial unit.** All of them: the Richards adaptive sub-step, the hourly
canopy and energy-balance step, the daily crop and soil step, and the season ledgers; the
leaf and the sunlit/shaded class, the leaf rank and tiller axis, the organ, the soil layer
and the point.

## 3. `csv_file.f90` (87 lines)

**Role.** Minimal CSV writer: one double-precision value, or a one-dimensional array of
them, to the current record of an already-open unit, with or without advancing. Reduced
form of the public-domain `csv_file` module of Arjen Markus (v1.2, 2006), keeping only the
double-precision routines the model calls. CSV conventions after Kernighan and Pike, *The
Practice of Programming* (1999).

**Calls.** Nothing; the body of the array routine is the `include` file `csv_file_1d.f90`.

**Called by.** `main_MMCO2Crop.f90` (`USE csv_file`; `CALL csv_write`).

**Outputs and their consumers.** Formatted CSV records on the unit the caller passes in.
In v1.0 that is the daily energy-balance diagnostic `output/test_energy_bal_2025.csv`,
written by `main_MMCO2Crop` through the generic `csv_write`.

**Time step / spatial unit.** Per write call; neither applies.

## 4. `csv_file_1d.f90` (33 lines)

**Role.** An `include` **body**, not a program unit: these statements form the body of
`csv_write_dble_1d` inside `csv_file.f90`. They write the elements of a 1-D array to the
current CSV record, advancing only after the last element.

**Calls.** The single-value `csv_write` of `csv_file`.

**Called by.** Nothing directly; textually included by `csv_file.f90` inside
`csv_write_dble_1d`. It has no entry in the `Makefile` `OBJS` list, correctly, because it
is not separately compiled.

**Outputs and their consumers.** One CSV record on the caller's unit.

**Time step / spatial unit.** Per write call; neither applies.

## 5. `count_lines.f90` (43 lines)

**Role.** Counts the records of the input file currently named in `file_name` (`SHARED`)
so the main program can size its daily arrays before reading them.

**Calls.** Nothing.

**Called by.** `main_MMCO2Crop.f90` (`CALL count_lines(NDays)`).

**Outputs and their consumers.** `nlines`, the number of records read including the header
line. `main_MMCO2Crop` subtracts the header to obtain `NDays` and allocates every daily
array with it.

**Time step / spatial unit.** Once per file at initialization, before the simulation
starts; no spatial unit.

## 6. `th.f90` (62 lines)

**Role.** Solar geometry of the day at the site: solar declination from day of year,
sunset and sunrise solar time at the site latitude, and daylength as their difference.

**Calls.** Nothing.

**Called by.** `CANOPY_MODEL.f90` — `CALL th (dTsr, Daylength, Solar_Decl)`.

**Outputs and their consumers.** `SolarDecl` (solar declination, rad), `Sunrise` (sunrise
solar time, h) and `dIval` (daylength, h), returned to `CANOPY_MODEL` as `Solar_Decl`,
`dTsr` and `Daylength`; `meteo` then takes `Solar_Decl` for the hourly solar elevation and
`dTsr` for the sunrise and sunset times of the day.

**Time step / spatial unit.** Once per simulated day, before the 24 hourly steps; the
point (one site latitude).

## 7. `solar.f90` (187 lines)

**Role.** Canopy radiation geometry for the two-big-leaf photosynthesis scheme. From solar
elevation and the direct/diffuse split of global radiation it computes the beam and
diffuse extinction coefficients, the penetration of direct and diffuse PAR through the
canopy, the PAR absorbed per unit sunlit and per unit shaded leaf area, and the sunlit and
shaded leaf area indices. Spherical leaf-angle distribution (Goudriaan 1977; de Pury and
Farquhar 1997).

**Calls.** Nothing.

**Called by.** `CANOPY_MODEL.f90`, once per hourly step, before the canopy-temperature and
assimilation solution.

**Outputs and their consumers.**
`Kdr` (beam extinction coefficient) — `radiation_module` for the shortwave split,
`energy_balance_module` and `et` for the net-radiation partition, `CANOPY_MODEL` for the
sunlit-weighted canopy nitrogen profile.
`Absor_IPAR_sunlit_h`, `Absor_IPAR_shaded_h` (absorbed PAR per unit sunlit / shaded leaf
area, W m-2 leaf) — converted to photon flux and passed to `assim` by `CANOPY_MODEL`.
`sunlit_lai`, `shaded_lai` — scale leaf fluxes to the canopy in `CANOPY_MODEL`.
`I_drrad_h`, `I_dfrad_h` (intercepted direct and diffuse solar radiation, W m-2 ground) —
diagnostics with no downstream consumer in v1.0.

**Time step / spatial unit.** Hourly; the point, canopy as two leaf classes.

## 8. `meteo.f90` (410 lines)

**Role.** Builds the hourly weather forcing from one daily weather record: solar geometry;
hourly global, direct and diffuse irradiance; hourly air temperature and relative humidity
(including the sunset anchor of the night branch and the offset that holds the 24-hour
air-temperature mean at the daily mean); vapor pressure and vapor pressure deficit; sky
emissivity; bulk net radiation.

**Calls.** No other routine.

**Called by.** `CANOPY_MODEL` (`canopy_model`): once per hour over 24 hours to form the
air-temperature offset, then once per committed time step.

**Outputs and their consumers.**
`SolarElev`, `dfrad_h`, `drrad_h` — `solar` (two-leaf PAR partition) and `CANOPY_MODEL`.
`HourTemp`, `Svp_s`, `HourVP_s`, `VpdRef_met` — `CANOPY_MODEL` gas exchange and energy
balance. `sky_emiss_comp` — `CANOPY_MODEL` component longwave (`ebh_LW_down` /
`LW_down_s`). `dTotrad` (daily global irradiance, J m-2 d-1) — `growth_module` and
`main_MMCO2Crop`. `alpha_mix_s`, `HourSW_net_s`, `HourLW_net_s`, `HourNetRad` — the
hourly radiation and aerodynamic diagnostics.

**Time step / spatial unit.** The daily record expanded to hourly, one call per hour; the
point, at the above-canopy reference height.

## 9. `radiation_module.f90` (174 lines)

**Role.** Component radiation for the three-node energy balance. Splits incoming shortwave
between canopy and soil by Beer's law with separate canopy and moisture-dependent soil
albedo, solves the sky-canopy-soil longwave exchange as a two-layer gray-body radiosity
problem to infinite order of reflection, and returns the net radiation of canopy and soil
as functions of the solved temperatures plus a shortwave-driven soil heat flux.

**Calls.** Nothing (it uses `SHARED` state only).

**Called by.** `energy_balance_module.f90` — `eb_residuals`, the Newton solve and the
actual-flux repartition, on every iterate of every hourly step.

**Outputs and their consumers.** `Rnc`, `Rns` (net radiation of canopy and soil, W m-2) and
`G` (soil heat flux, W m-2) — `energy_balance_module` (three-node temperature solution and
flux repartition). The `SHARED` `rad_*` ledger (shortwave absorbed and reflected by canopy
and soil, longwave down, emitted, exchanged and net per node, outgoing longwave at the top,
soil and mixed albedo, shortwave and longwave transmissivity) — read by `CANOPY_MODEL` and
`energy_balance_module` for the hourly energy-balance diagnostics.

**Time step / spatial unit.** Hourly, once per energy-balance iterate; the canopy layer and
the soil surface.

## 10. `soil_heat_module.f90` (368 lines)

**Role.** One-dimensional soil heat conduction. Builds the heat grid (the hydrology layers
extended geometrically to about 5 m), initializes the prognostic temperature profile,
updates the de Vries volumetric heat capacity and the Johansen / Cote-Konrad thermal
conductivity from the current water content, supplies the linear surface-conductance
coefficients for the energy balance, and advances the profile with a backward-Euler
finite-volume step solved by the Thomas algorithm.

**Calls.** Nothing outside itself.

**Called by.** `main_MMCO2Crop.f90` (`initialize_soil_heat`, once before the day loop);
`energy_balance_module.f90` (`update_soil_thermal_properties`,
`soil_heat_surface_coeffs`, `soil_heat_surface_flux`, `advance_soil_heat`, every
energy-balance sub-step).

**Outputs and their consumers.** `Tsoil_prof(:)` (C), `Cv_soil(:)` (J m-3 K-1),
`lambda_soil(:)` (W m-1 K-1), `tk_heat(:)` / `n_heat_tot` (heat-grid geometry), the surface
conductance `k_face` and node-1 temperature, the conductive surface flux `G` (W m-2) and
the storage residual `heat_resid` (J m-2). The fluxes and coefficients are consumed by
`energy_balance_module`; the profile and the residual are read by `main_MMCO2Crop` for the
soil-temperature output and the heat-closure diagnostics.

**Time step / spatial unit.** One energy-balance sub-step (hourly); the soil layer on the
heat grid.

## 11. `soiln_module.f90` (149 lines)

**Role.** Crop-residue decomposition and the microbial carbon-nitrogen turnover that feeds
the soil nitrogen chain. Decays a fast (metabolic) and a slow (structural) residue pool by
first-order kinetics with temperature and water factors, and resolves from the carbon-use
efficiency and the microbial C:N ratio whether each pool's decay mineralizes nitrogen to
ammonium or immobilizes it from the mineral pools; microbial biomass then turns over and
re-mineralizes its nitrogen.

**Calls.** Nothing outside itself.

**Called by.** `main_MMCO2Crop.f90` — once per day in the soil-nitrogen sequence, and once
per pre-sowing day of the residue spin-up.

**Outputs and their consumers.** The residue (`res_fast_C/N`, `res_slow_C/N`), microbial
(`soil_mic_C/N`) and mineral (`N_NH4_soil`, `N_soil_min`) pools are updated in place, and
the day's fluxes are returned: decomposed C and N, respired C, net mineralization to
ammonium, immobilization from ammonium and from nitrate, microbial turnover N (all
g m-2 d-1). `main_MMCO2Crop` consumes them: the ammonium gained enters the NH4 to NO3
nitrification pathway and, in layered mode, is distributed over the layered NH4/NO3
stocks; the respired carbon drives denitrification; the fluxes are written to
`output/daily_organic_N.csv` and its closure check.

**Time step / spatial unit.** Daily; the residue pools and the soil mineral nitrogen as
bulk stocks (g m-2).

## 12. `assim.f90` (672 lines)

**Role.** Leaf gas exchange for one leaf class in one hour. Solves the coupled
photosynthesis-conductance system: the Farquhar et al. (1980) minimum of the Rubisco-,
electron-transport- and triose-phosphate-limited carboxylation rates with the Busch et al.
(2018) photorespiratory nitrogen-assimilation extension, the Leuning (1995) stomatal
conductance, and where selected the soil-plant-atmosphere leaf water potential.

**Calls.** Nothing outside itself.

**Called by.** `CANOPY_MODEL.f90` — `eval_canopy_assim` (sunlit and shaded classes) and
`tleaf_classes_diag` (the flag-leaf diagnostic call).

**Outputs and their consumers.** `Rd` (day respiration), `A_assim` (gross assimilation),
`An_assim` (net assimilation), `gsw_1` (stomatal conductance to water vapour), `Ci`
(intercellular CO2) and `N_total_assim` (nitrate reduced in the chloroplast), all per unit
leaf area — consumed by `CANOPY_MODEL`, which integrates them to the canopy and hands the
canopy fluxes to the energy balance and to `growth_module`. The soil-plant-atmosphere
diagnostics and the leaf-water-potential solve counters are written back into `SHARED`.

**Time step / spatial unit.** Hourly, per leaf class; one leaf (sunlit, shaded or flag
leaf).

## 13. `energy_balance_module.f90` (1 001 lines)

**Role.** The three-node canopy / soil / canopy-airspace energy balance, solved once per
hour and once per iterate of the canopy-temperature solve. The canopy-airspace temperature
`Tac` is a solved node, not a post-hoc diagnostic. Series-resistance sensible-heat network
after Norman et al. (1995) and Kustas and Norman (1999); latent heat uses an explicit
canopy-airspace vapour node `e_ac`. Given `Tcan` and `rcs` it solves the three unknowns
(`Tsoil`, `Tac`, `e_ac`) by a damped Newton iteration on a numerical 3x3 Jacobian; the
outer `Tcan` root-find and the assimilation to `gsw` to `rcs` fixed point live in
`CANOPY_MODEL::solve_biophysics`.

**Calls.** `radiation_module` (`component_radiation`); `soil_heat_module`
(`soil_heat_surface_flux`, `soil_heat_surface_coeffs`, `advance_soil_heat`,
`update_soil_thermal_properties`).

**Called by.** `CANOPY_MODEL.f90` — `solve_biophysics` and `compute_daily_canopy_fluxes`;
`main_MMCO2Crop.f90` — the daily repartition (`repartition_day_actual`) and the hourly
writers.

**Outputs and their consumers.** The resistances `raa`, `rca`, `rsa`, `rss`; the solved
soil-surface temperature, canopy-airspace temperature `Tac` and vapour pressure `e_ac`;
the committed component fluxes `Hc`, `Hs`, `LEc`, `LEs`, `G` and the actual daily
repartition. Consumed by `CANOPY_MODEL` (canopy temperature and the leaf-to-air vapour
pressure deficit that drive `assim`), by `hydro` through the soil-evaporation and
transpiration demand, and by `main_MMCO2Crop` for the hourly energy-balance and radiation
diagnostic files.

**Time step / spatial unit.** Hourly, once per iterate; the canopy as one big leaf, the
soil surface, and the canopy airspace between them, with the soil heat flux taken from the
1-D soil temperature profile layer by layer.

## 14. `hydro.f90` (1 661 lines)

**Role.** Soil water movement by the one-dimensional mixed-form Richards equation, solved
over one day in adaptive sub-steps by implicit backward Euler with Picard (and Newton)
iteration. Hydraulics are van Genuchten (1980) retention with Mualem (1976) conductivity;
root water uptake is distributed over the rooted layers by the Jackson et al. (1996)
profile and reduced by the Feddes et al. (1978) response.

**Calls.** `IEEE_ARITHMETIC` finiteness tests on the iterates; nothing else outside itself.

**Called by.** `main_MMCO2Crop.f90` — the daily loop, once per pass of the
evapotranspiration / Richards coupling iteration, through an explicit interface.

**Outputs and their consumers.** The updated layer water contents `theta_ly` and matric
heads, the root-uptake sink per layer, actual transpiration and soil evaporation,
infiltration, runoff, bottom flux and the surface store, and the per-day water-balance
ledger. Consumed by `main_MMCO2Crop` (the water balance and its diagnostic files), by
`energy_balance_module` (the committed surface head and the soil latent-heat supply limit),
by `CANOPY_MODEL` through the soil water potential that enters the stomatal response, and
by `growth_module` through the water-stress state.

**Time step / spatial unit.** The adaptive Richards sub-step inside one day; the ten
cell-centered soil layers to 1.20 m plus the surface store.

## 15. `mineral_module.f90` (1 086 lines)

**Role.** Grain minerals, six elements (Fe, Zn, Ca, Cu, Mg, Mn). One daily step per
element: soil-solution concentration, uptake by mass flow, diffusion to the root surface
and element-specific active terms, allocation to the growing organs, remobilization from
senescing tissue to the grain, and a per-element mass balance.

**Calls.** Nothing outside itself.

**Called by.** `main_MMCO2Crop.f90` — `read_mineral_params` and `init_mineral_pools` once
per run, `update_minerals_daily` once per day, after grain filling and before the
senescence update.

**Outputs and their consumers.** The labile soil pool `M_soil_labile` and the organ pools
`M_root`, `M_leaf`, `M_stem`, `M_chaff`, `M_grain`, `M_dead`; the daily diagnostics
`M_uptake_mf`, `M_uptake_diff`, `M_uptake_ac`, `M_uptake_hats`, `M_remob_grain`, `M_soln`,
`M_cap_binding`, `M_balance_resid`; and the season ledgers of the mineral records. All are
consumed by `main_MMCO2Crop`, which converts the grain pool to a concentration (mg kg-1)
for output and writes `output/daily_minerals.csv`, `daily_mineral_supply.csv`,
`daily_mineral_uptake_terms.csv` and `mineral_sink_cap.csv`. No carbon, nitrogen, water or
energy state is written here.

**Time step / spatial unit.** Daily, per element; one labile soil pool over the whole
0-0.6 m profile plus the crop's organ pools, every pool in mg element m-2 of ground.

## 16. `phenology_module.f90` (828 lines)

**Role.** Wheat development on one thermal clock: the daily thermal-time increment and its
vernalization and photoperiod factors, the cumulative thermal development units, the dates
of the eight development stages, and the continuous 0-2 development stage that the growth,
grain-filling and senescence laws run on.

**Calls.** Its own cultivar-file reader; no routine outside itself. The daily driving
temperatures are supplied by the caller.

**Called by.** `main_MMCO2Crop` (`PROGRAM mmco2crop`), which calls
`initialize_phenology_params`, `read_phenology_cultivar`, `initialize_phenology_state`,
`bootstrap_phenology_state` and, once per day, `update_phenology_day`.

**Outputs and their consumers.** The derived types `PhenologyParams` (genotype thresholds,
cardinal temperatures, sensitivities, sowing date, observed stage dates, the initial organ
mass and nitrogen fractions, the optional cultivar canopy-architecture traits) and
`PhenologyState` (`cum_tdu`, `cum_vd`, `current_dtt`, `current_dtdu`, `current_vf`,
`current_df`, `photoperiod`, `dap`, `istage`, `stage_dap`, `stage_reached`,
`tdu_at_emerge`, `ts_done`); the functions `calc_dev_stage_0_2` and
`calc_dev_stage_rate_0_2`, whose 0-2 development stage and daily increment are read by
`growth_module` (partition fractions, specific leaf area, root depth, senescence) and by
the grain-filling and canopy clocks in `main_MMCO2Crop`; and `calc_dtt` and
`TS_TDU_MCMASTER`, exported for the optional soil-temperature canopy clock.

**Time step / spatial unit.** Daily, one call per simulated day; the whole plant, the
stages being properties of the apex.

## 17. `CANOPY_MODEL.f90` (887 lines)

**Role.** The canopy integration of leaf gas exchange and the canopy-temperature solve:
24 hourly steps per day, each a bracketed root-find for the canopy temperature wrapped
around a nested fixed point that makes assimilation, stomatal conductance and the
leaf-surface vapour pressure deficit mutually consistent, and each returning the hourly
demand profile the soil-water repartition later replays.

**Calls.** `energy_balance_module` (`compute_eb_resistances`, `solve_eb_3node_inner`,
`esat_pa`, `h_soil_surface`, `raa_heat_stability`) and the external routines `assim`
(through the explicit interface in the file), `meteo`, `solar` and `th`.

**Called by.** `main_MMCO2Crop.f90` — `compute_daily_canopy_fluxes` once per pass of the
daily coupling loop; `canopy_n0_top` is also called there, for the canopy-nitrogen profile
output and the lamina-nitrogen ceiling.

**Outputs and their consumers.** Through the argument list: transpiration and
soil-evaporation demand, gross and net assimilation, nitrate reduced, canopy, air and
soil-surface temperature, vapour pressure deficit, net radiation, soil heat flux, the
latent and sensible heat fluxes of canopy and soil, stomatal conductance, intercellular
CO2 and the five resistances, each also as a daylight mean. They are consumed by
`main_MMCO2Crop`, which drives the soil-water and soil-heat solvers, `growth_module` and
`phenology_module` with them. Written into `SHARED`: the hourly demand profiles the
post-Richards repartition replays (`ebh_*`), the buffered hourly energy-balance and
radiation rows, the solved three-node state, the sunlit and shaded leaf temperatures and
the solver counters.

**Time step / spatial unit.** Hourly with a sub-hourly iteration; the canopy as two leaf
classes, sunlit and shaded, with depth entering only through the vertical leaf-nitrogen
profile. The day's means and totals are returned to the caller.

## 18. `growth_module.f90` (5 549 lines)

**Role.** The crop carbon and nitrogen economy on a daily step. It partitions the day's
assimilate between lamina, true stem, root and spike, charges maintenance and growth
respiration, builds and senesces leaf area, deposits and draws the stem
water-soluble-carbohydrate reserve, pays every structural increment in reduced nitrogen,
turns organ nitrogen over after flowering and loads the grain, sets grain number and
kernel mass, and descends the root front.

**Calls.** Nothing outside itself (`USE SHARED, ONLY:` an explicit list, and `USE TYPES,
ONLY: wp`).

**Called by.** `main_MMCO2Crop.f90` only (`USE GROWTH_MODULE, ONLY: ...`), once per day, in
the order `initialize_growth_state` and the two bootstrap routines at the start of the run,
then per day `allocate_daily_biomass`, `update_structure`, the organ-canopy routines
(`update_canopy_leaf_ranks`, `update_canopy_blade_area`, `update_canopy_tillering`,
`update_leaf_expansion_water`, `update_canopy_leaf_senesc`, `update_canopy_sheath_ear`,
`update_canopy_sink_state`), `apply_leaf_senescence`,
`allocate_yield_components_daily`, `apply_organN_turnover_daily`,
`allocate_grain_N_source_sink_daily` and `update_root_depth`.

**Outputs and their consumers.** The organ dry-matter pools `W_green`, `W_dead`, `W_stem`,
`W_root`, `W_storage`, `W_grain`, `W_chaff` and the reserve `W_reserve`; leaf area index
and, under the organ law, the per-rank and per-axis green and dead areas; canopy height;
root depth; the organ nitrogen pools `N_leaf`, `N_stem`, `N_storage`, `N_chaff`,
`N_dead_leaf` and the mobile vegetative pool `N_mobile_veg`; grain number and kernel mass;
and the daily carbon and nitrogen ledgers and closure residuals. Leaf area and leaf
nitrogen are consumed by `CANOPY_MODEL` and `assim` (the canopy nitrogen profile and gas
exchange), the root depth and the transpiration demand by `hydro`, the organ dry-matter
increments by `mineral_module` (the mineral allocation shares), and every pool and ledger
by `main_MMCO2Crop`, which writes them out.

**Time step / spatial unit.** Daily; the plant organ and, where the organ-level canopy is
active, the tiller axis, the leaf rank and the individual blade or sheath. All pools are
per square metre of ground.

## 19. `main_MMCO2Crop.f90` (13 975 lines)

**Role.** The program. It reads the 183 positional records of `input/driver_main.txt` and
the deck input files (weather, soil water and van Genuchten-Mualem layers, soil nitrogen,
fertilizer and organic-residue schedules, cultivar phenology, yield and mineral
parameters), refuses every retired driver option, initializes the soil water, heat,
nitrogen and mineral states, then runs the daily loop for the whole season: phenology, the
organ-level canopy chronology, the hourly canopy energy balance and gas exchange iterated
against the Richards soil-water solver, interception, infiltration and ponding, soil heat,
the soil and plant nitrogen economy, carbon allocation and yield formation, the mineral
economy, and the daily water, energy, carbon and nitrogen closures.

**Calls.** `SHARED`; `TYPES, ONLY: wp, ik`; `csv_file`; `CANOPY_MODEL`
(`compute_daily_canopy_fluxes`, `canopy_n0_top`); `energy_balance_module`
(`repartition_day_actual` and the hourly writers); `soil_heat_module`
(`initialize_soil_heat`); `GROWTH_MODULE` (allocation, structure, senescence, the
organ-level canopy chain, grain number and filling, organ nitrogen turnover, the internode
height law); `MINERAL_MODULE`; `SOILN_MODULE` (`decompose_residue_day`);
`PHENOLOGY_MODULE`; and the external subroutines `hydro` and `count_lines` through explicit
interfaces.

**Called by.** Nothing — this is the executable's main program (`PROGRAM mmco2crop`).

**Outputs and their consumers.** The 88 files of `output/` — the daily and hourly
energy-balance, radiation, gas-exchange and leaf-temperature files, the soil profile,
water-balance and Feddes diagnostics, the soil- and plant-nitrogen, mineral,
canopy-architecture, reserve, grain-number, grain-filling and organ-nitrogen files, the
season tables — and the season summary on standard output, which the documented run command
captures as `run.log`. The consumer is the user; nothing inside the model reads them. The
per-file list is in `FILE_INVENTORY.md`.

**Time step / spatial unit.** All four time steps (the Richards adaptive sub-step, the
hourly canopy and energy-balance step, the daily crop and soil step, the season) and all
spatial units (the sunlit and shaded leaf class, the leaf rank and tiller axis, the plant
organ, the soil layer, the point).
