<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Input manual — M&M-CO2Wheat v1.0

Updated 2026-10-21 for the v1.0 release. This file states what a user must provide, in
which file, in which order and in which units; which entries are mandatory and which have
a documented default; and what value each of the 183 driver records may take in this
release.

A complete working example of everything below is `testcase/input/`, the Part 1 reference
simulation. The quickest way to build a deck of your own is to copy that directory and
edit it.

## 0. How the model finds its inputs

The executable opens every input by the **relative** path `input/...` and writes to the
relative path `output/...`, so it must be run from a directory holding both, and
`output/` must exist before the run starts (see `README.md`, "Running the test case").
Only one input is found by a name taken from another file: `input/driver_main.txt`
record 1 names the daily weather file. Every other name below is fixed in the code.

---

## 1. The minimum input set

| file | status | format |
|---|---|---|
| `input/driver_main.txt` | **mandatory** | 183 positional records, one value per line, optional `!` comment after the value (section 2) |
| the weather file named by driver record 1 | **mandatory** | CSV, header row then one row per day (section 3) |
| `input/driver_swtr.txt` | **mandatory** | positional; a missing file stops the run (section 4) |
| `input/SoilLayer_vgm.csv` | **mandatory** | CSV with a fixed header; fixed file name (section 5) |
| `input/phenology_cultivar.txt` | **mandatory** | positional; rows 1-35 required, 36+ optional with compiled defaults (section 6) |
| `input/cultivar_yield.csv` | **mandatory** | positional, value first; 6, 8 or 10 values (section 7) |
| `input/driver_N.csv` | **mandatory** | positional, 70 values (section 8) |
| `input/driver_minerals.csv` | **mandatory in practice** | CSV, header plus one row per element (section 9) |
| `input/N_fertilizer.csv` | **mandatory** | CSV, count line then one row per application (section 10) |
| `input/N_organic_treatments.csv` | optional | CSV, count line then one row per event; absent, or a count of 0, is the no-residue baseline (section 11) |
| `input/tsoil_init.csv` | optional | CSV `z_m,T_C` anchors; absent means the driver record-50 initialization is used (section 12) |
| `input/cultivar_partition.csv` | optional | CSV; absent means the WOFOST baseline leaf-partition window, and the run prints a NOTE saying so (section 13) |

**Mandatory versus defaulted, in one rule.** A file marked mandatory must exist: the run
stops with a message naming it. Inside a file, an entry is mandatory unless this manual
names a default for it; where a default exists it is an **end-of-file default** — the
reader detects the missing trailing record, keeps the compiled value, and the run is
byte-identical to one in which the record is absent. A record that is present but
unreadable always stops the run with an `ERROR STOP` naming it; the model never
substitutes a value silently.

`input/driver_minerals.csv` is the one partial exception. A missing file produces the
warning `input/driver_minerals.csv not read fully; mineral pools may be zero` and the run
continues with zero mineral pools, which is never a result you want; and with driver
record 131 = 1 — the shipped value — a file carrying fewer than 18 numeric columns is
refused outright rather than run with silently zeroed availabilities.

---

## 2. `input/driver_main.txt` — the run-control file

> **Which copy is the released configuration.** A deck's own `input/driver_main.txt` is
> the file the model reads. The configuration v1.0 was frozen, gated and evaluated on is
> `docs/V1_CONFIGURATION_driver_main_v1.0.txt`, which is byte-identical to
> `testcase/input/driver_main.txt` and to the deck `decks/WDF25_Ruta_aCO2_WW/input/driver_main.txt`.
> Cite and compare against those. The development repository also carries a template
> `driver_main.txt` that differs from the released configuration at record 38 (`-0.35`
> against the released `-0.3854`); that template is **not distributed in this package**, and
> no file in this archive is a template.

**The file is POSITIONAL.** 183 records, in order, one value per line. Anything after a
`!` on a line is a comment and is stripped before the value is read. Record 1 is a file
path; records 2-183 are numbers. Blank-looking records (record 73 in the reference deck)
are permitted where this manual says so.

> **Never delete or renumber a record.** Retired records keep their position and are
> marked RESERVED, because removing one line silently reassigns every switch below it.

### 2.1 What v1.0 accepts

The release contains the code of the process options it uses and of the alternative
formulations of Part 1 Table 6, which are kept and documented as built-and-inactive. Code
for every other option was removed, so **a retired value is refused at startup**: the run
prints

```
REFUSING TO RUN: driver_main.txt record <N> (<variable>) = <value>
  M&M-CO2Wheat v1.0 implements only these values of record <N>: <set>
```

and exits non-zero. It never falls through to a different mechanism or to an undefined
value. 123 mode and switch records are restricted this way; four more (119, 120, 170,
182) are RESERVED rows that keep their positions. The remaining 56 records — the plain
site, soil, weather, canopy, management and cultivar parameters — are **not restricted at
all**, so a new site can use its own values freely.

The refusal behaviour was tested case by case before release, in 32 cases covering every
class of restricted record, all of which refused the retired value and exited non-zero.
The allowed sets tabulated below are taken from the `reject_retired_option` calls in
`main_MMCO2Crop.f90` itself, and this table is the authoritative list for a user of the
release.

### 2.2 The 183 records

"accepted in v1.0" is the set the validator admits. "not restricted" means the record is
an ordinary parameter and any value its own physical meaning allows is accepted. The
deck-value column is the Part 1 reference simulation `testcase/input/driver_main.txt`.

| rec | variable | v1.0 deck value | accepted in v1.0 | what it sets |
|---:|---|---|---|---|
| 1 | -- | `input/input_daily_meteo.csv` | not restricted | path of the daily weather file |
| 2 | -- | `2` | not restricted | number of developmental stages |
| 3 | -- | `31.8953` | not restricted | site latitude (deg) |
| 4 | -- | `34.8106` | not restricted | site longitude (deg) |
| 5 | -- | `0.000` | not restricted | initial leaf area index (m2 m-2) |
| 6 | -- | `7.00` | not restricted | maximum leaf area index (m2 m-2); read but reaches no state variable in v1.0 |
| 7 | -- | `0.022` | not restricted | specific leaf area (m2 g-1); inert while record 69 = 1 (see note below the table) |
| 8 | -- | `200.0` | not restricted | plant density (plants m-2) |
| 9 | -- | `0.040` | not restricted | initial leaf nitrogen fraction (g N g-1 DM) |
| 10 | -- | `0.0` | not restricted | initial plant height (m) |
| 11 | -- | `2.0` | not restricted | weather reference height (m) |
| 12 | -- | `2.5` | not restricted | canopy wind attenuation coefficient (-) |
| 13 | -- | `2.5` | not restricted | canopy eddy attenuation coefficient (-) |
| 14 | -- | `2.8927` | not restricted | Jarvis (1976) stomatal coefficient a1; inert, no consumer in v1.0 |
| 15 | -- | `0.0040` | not restricted | Jarvis (1976) stomatal coefficient a2; inert, no consumer in v1.0 |
| 16 | -- | `0.027` | not restricted | mean leaf width (m) |
| 17 | -- | `0.02` | not restricted | dry surface layer thickness (m) |
| 18 | -- | `0.42` | not restricted | soil porosity (m3 m-3) |
| 19 | -- | `0.18` | not restricted | pore size distribution index (-) |
| 20 | -- | `0.39` | not restricted | saturation water content (m3 m-3) |
| 21 | -- | `0.20` | not restricted | initial water content of the top layer (m3 m-3) |
| 22 | -- | `420.0` | not restricted | atmospheric CO2 (umol mol-1) |
| 23 | `water_stress_mode` | `1` | 1 | locus of the water-stress factor: 1 = on stomatal conductance, 2 = on assimilation |
| 24 | `use_soil_water_stress` | `1` | 1 | soil-water stress on or off |
| 25 | -- | `1` | not restricted | allow drainage beyond the last layer |
| 26 | `use_explicit_infiltration` | `0` | 0 | explicit infiltration scheme on or off |
| 27 | `integration_mode` | `2` | 2 | time integration: 1 = Gauss daylight quadrature, 2 = full 24 hours |
| 28 | -- | `0.001` | not restricted | root growth rate (m d-1) |
| 29 | -- | `1.21` | not restricted | shoot nitrate-reduction ceiling N_max (umol m-2 s-1) |
| 30 | -- | `1` | not restricted | dynamic eddy attenuation on or off |
| 31 | `meteo_debug` | `0` | 0 | meteorology debug prints |
| 32 | -- | `2` | not restricted | irrigation placement: 1 = above canopy, 2 = below canopy |
| 33 | -- | `0.0` | not restricted | pre-start precipitation (mm) |
| 34 | -- | `1.65` | not restricted | soil bulk density (g cm-3) |
| 35 | `k_avg_mode` | `0` | 0 | internodal conductivity averaging: 0 = arithmetic, 1 = geometric, 2 = upstream |
| 36 | `hyd_model` | `0` | 0 | soil hydraulic model: 0 = van Genuchten-Mualem, 1 = Brooks-Corey, 2 = Campbell |
| 37 | `init_profile_mode` | `2` | 0, 2 | initial water profile: 0 = as read, 1 = hydrostatic, 2 = monotonic top to bottom |
| 38 | -- | `-0.3854` | not restricted | initial bottom pressure head (m, negative = suction); used when record 37 > 0 |
| 39 | `energy_balance_mode` | `2` | 2 | energy balance: 2 = the three-node canopy-soil-airspace network, the only accepted value |
| 40 | `eb_freeconv_fallback` | `0` | 0 | free-convection fallback in the energy balance |
| 41 | -- | `1000.0` | not restricted | boundary-layer height z_i (m) |
| 42 | -- | `1.2` | not restricted | gustiness coefficient beta (-) |
| 43 | -- | `0.60` | not restricted | soil quartz fraction (-) |
| 44 | -- | `2.00e6` | not restricted | volumetric heat capacity of soil minerals (J m-3 K-1) |
| 45 | -- | `2.50e6` | not restricted | volumetric heat capacity of soil organic matter (J m-3 K-1) |
| 46 | -- | `4.18e6` | not restricted | volumetric heat capacity of water (J m-3 K-1) |
| 47 | -- | `1.25e3` | not restricted | volumetric heat capacity of air (J m-3 K-1) |
| 48 | -- | `0` | not restricted | soil-heat bottom boundary type |
| 49 | -- | `15.0` | not restricted | deep soil temperature (C) |
| 50 | `soil_T_init_method` | `1` | 1 | soil temperature initialization: 0 = uniform, 1 = damped |
| 51 | -- | `1.0` | not restricted | Jarvis compensatory root uptake omega (-); 1.0 = off |
| 52 | `residue_evap_mode` | `1` | 1 | residue evaporation suppression mode |
| 53 | `residue_Am_ha_t` | `0.45` | not restricted | residue area per unit mass (ha t-1) |
| 54 | `residue_fsupp` | `0.75` | not restricted | residue evaporation suppression fraction (-) |
| 55 | `carbon_reserve_mode` | `1` | 1 | carbon reserve mode |
| 56 | `stem_reserve_frac` | `0.20` | not restricted | stem reserve fraction of stem mass (-) |
| 57 | `reserve_remob_rate` | `0.05` | not restricted | reserve remobilization rate (d-1) |
| 58 | `reserve_start_dvs2` | `1.00` | not restricted | development stage at which reserve use starts (dvs 0-2) |
| 59 | `reserve_end_dvs2` | `2.00` | not restricted | development stage at which reserve use ends (dvs 0-2) |
| 60 | `reserve_conv_eff` | `0.95` | not restricted | reserve conversion efficiency (-) |
| 61 | `stem_spike_realloc_delta` | `0.0` | not restricted | stem-to-spike reallocation increment (-) |
| 62 | `leafN_drawdown` | `0` | 0 | leaf nitrogen drawdown law |
| 63 | `grainN_partition` | `0` | 0 | grain nitrogen partition law |
| 64 | `canopy_leaf_ranks` | `1` | 1 | organ-level canopy: leaf rank chronology |
| 65 | `canopy_blade_area` | `1` | 1 | organ-level canopy: blade area |
| 66 | `canopy_tillering` | `1` | 1 | organ-level canopy: tillering |
| 67 | `canopy_leaf_senesc` | `1` | 1 | organ-level canopy: leaf senescence |
| 68 | `canopy_sheath_ear` | `1` | 1 | organ-level canopy: sheath and ear area |
| 69 | `canopy_sink_lai` | `1` | 1 | organ-level canopy: sink-limited LAI |
| 70 | `rca_wind_attenuation` | `1` | 1 | canopy aerodynamic resistance from wind attenuation |
| 71 | `wind_hourly_mode` | `0` | 0 | hourly wind mode; 0 = the daily mean held flat |
| 72 | `rwa_lai_blend` | `1.0` | not restricted | LAI blend for the record-70 resistance (m2 m-2) |
| 73 | `wind_hourly_file` | `(blank)` | not restricted | path of a measured sub-daily wind file; used only when record 71 >= 1 |
| 74 | `wind_meas_hgt_m` | `-1.0` | not restricted | measurement height of that wind series (m); -1 = not supplied |
| 75 | `rss_theta_s_consistency` | `1` | 1 | consistency of the soil-surface resistance with theta_s |
| 76 | `canopy_n_base_profile` | `1` | 1 | vertical leaf-nitrogen profile: 1 = the measured wheat form (uniform over the flag-leaf zone, exponential below it), the only accepted value |
| 77 | `ci_solve_mode` | `1` | 1 | intercellular CO2 solution mode |
| 78 | `ft_pot_unstressed` | `1` | 1 | unstressed potential transpiration basis |
| 79 | `ws_two_loci` | `1` | not restricted | two-loci soil-water stress |
| 80 | `feddes_p2_hydrus` | `1` | not restricted | Feddes p2-node unit convention |
| 81 | `spac_mode` | `3 -1.5` | not restricted | soil-plant-atmosphere mode and its leaf water potential threshold (MPa) |
| 82 | `tiller_carbon_gate` | `1` | 1 | carbon gate on tiller birth |
| 83 | `leaf_expansion_water` | `1` | 1 | water limit on leaf expansion |
| 84 | `lai_from_organs` | `1` | 1 | leaf area index taken from the organ-level canopy |
| 85 | `rca_teh_full` | `0` | 0 | full Teh canopy aerodynamic coupling |
| 86 | `grainN_sink_mech` | `1` | 1 | grain nitrogen sink law |
| 87 | `leafN_cap_route` | `1` | 1 | route of the leaf nitrogen capacity limit |
| 88 | `sowing_model_day` | `0` | not restricted | sowing day as a model day; 0 = taken from the cultivar file |
| 89 | `grain_number_mech` | `1` | 1 | grain number mechanism |
| 90 | `floret_survival_mech` | `4` | 4 | floret survival mechanism |
| 91 | `grain_fill_mech` | `1` | 1 | grain filling mechanism |
| 92 | `stem_wsc_mech` | `1` | 1 | stem reserve mechanism |
| 93 | `soil_n_mech` | `0` | 0 | soil nitrogen mechanism |
| 94 | `canopy_n_flag_zone` | `3` | 3 | flag-leaf zone of the canopy nitrogen profile |
| 95 | `tleaf_class_diag` | `0` | 0 | per-class leaf temperature diagnostic |
| 96 | `sw_radiosity_mode` | `0` | 0 | shortwave radiosity mode |
| 97 | `itemF_Nclose_mode` | `1` | 1 | nitrogen closure mode |
| 98 | `spike_sink_mech` | `0` | 0 | spike sink mechanism |
| 99 | `tiller_survival_mech` | `4` | 4 | tiller survival and fertile ear population mechanism |
| 100 | `stem_demand_mech` | `0` | 0 | stem demand mechanism |
| 101 | `surplus_fate_mech` | `0` | 0 | fate of refused (surplus) carbon |
| 102 | `stem_sheath_sep_mode` | `0` | 0 | separation of true stem from leaf sheath in the stem pool |
| 103 | `sdm4_dur_bracket` | `0` | 0 | duration bracket of record 100 mode 4 |
| 104 | `sdm4_rate_bracket` | `0` | 0 | rate bracket of record 100 mode 4 |
| 105 | `wsc_reserve_law` | `2` | 2 | stem water-soluble-carbohydrate reserve law |
| 106 | `wsc_law_member` | `0` | 0 | ceiling-clock member of record 105 |
| 107 | `wsc_rel_window_member` | `2` | 2 | release-window member of record 105 |
| 108 | `height_reserve_basis` | `0` | 0 | plant-height allometry basis |
| 109 | `height_internode_law` | `0` | 0, 1 | canopy height from internode elongation |
| 110 | `hgt_int_window_member` | `0` | 0, 1, 2, 3, 4, 5, 6 | onset and duration bracket of record 109 |
| 111 | `hgt_int_length_member` | `0` | 0, 1, 2, 3, 4, 5 | length-potential bracket of record 109 |
| 112 | `hgt_int_lineic_member` | `0` | 0, 1, 2, 3, 4 | lineic-mass bracket of record 109 |
| 113 | `hgt_int_ear_member` | `0` | 0, 1, 2, 3 | ear and peduncle bracket of record 109 |
| 114 | `mineral_share_basis` | `0` | 0 | basis of the mineral allocation shares |
| 115 | `gnf_ng0_basis` | `0` | 0 | per-kernel grain-nitrogen sink basis |
| 116 | `limb_census_mode` | `0` | 0 | limiting-limb census |
| 117 | `chaff_growth_mode` | `0` | 0 | chaff growth mode |
| 118 | `kernel_mass_ceiling` | `1` | 0, 1 | measured kernel-mass ceiling on grain filling |
| 119 | `RESERVED` | `0` | **0 only** (retired, unavailable in v1.0) | RESERVED (rejected; no active code) |
| 120 | `RESERVED` | `0` | **0 only** (retired, unavailable in v1.0) | RESERVED (rejected; no active code) |
| 121 | `n_payment_gate` | `2` | 2 | nitrogen payment gate |
| 122 | `lai_n_arrest` | `2` | 2 | nitrogen arrest of leaf-area expansion |
| 123 | `tiller_n_afford` | `1` | 1 | nitrogen term in the tiller survival ledger |
| 124 | `root_front_law` | `2` | 2 | root front descent law (1/2/3 = 7.9/9.6/11.3 mm d-1) |
| 125 | `kernel_potential_law` | `1` | 1, 2, 3 | kernel potential law |
| 126 | `mineral_grain_sink_cap` | `1` | 1 | cap on the grain mineral sink |
| 127 | `repro_surplus_fate` | `1` | 1 | fate of the reproductive-phase surplus |
| 128 | `n_construction_status` | `2` | 2 | status-dependent construction nitrogen |
| 129 | `sink_full_carbon_fate` | `1` | 1 | fate of the carbon surplus when the grain sink is full |
| 130 | `grain_n_min_conc` | `3` | 1, 2, 3 | minimum grain nitrogen concentration law |
| 131 | `mineral_supply_law` | `1` | 1 | site-level mineral availability law |
| 132 | `leafn_hold_law` | `0` | 0 | leaf nitrogen holding capacity law |
| 133 | `leafn_turnover_law` | `0` | 0 | leaf nitrogen turnover law |
| 134 | `n_afford_priority_law` | `0` | 0 | nitrogen affordability and mortality-priority law |
| 135 | `n_payment_luxury_pass2` | `0` | 0 | second-pass stem luxury nitrogen draw |
| 136 | `n_afford_standing_basis` | `0` | 0 | re-point the record-134 signal onto standing nitrogen status |
| 137 | `floret_n_law` | `0` | 0 | nitrogen term on floret survival per spikelet |
| 138 | `spike_n_priority` | `0` | 0 | spike and chaff as their own structural-nitrogen tier |
| 139 | `n_uptake_diffusion` | `1` | 1 | diffusive term of soil nitrogen uptake |
| 140 | `mineral_supply_mech` | `1` | 1 | mineral supply mechanism |
| 141 | `lamina_n_remob_pre` | `1` | 1 | pre-anthesis lamina nitrogen remobilization |
| 142 | `lamina_n_expansion` | `0` | 0 | nitrogen-limited lamina expansion |
| 143 | `n_status_construction_price` | `0` | 0, 1 | whole-crop nitrogen-status-dependent construction price |
| 144 | `carbon_coupled_leaf_expansion` | `0` | 0, 1 | carbon-coupled leaf-area expansion |
| 145 | `leaf_maintenance_basis` | `1` | 1 | leaf maintenance respiration basis |
| 146 | `par_reflection_basis` | `1` | 1 | PAR reflection basis |
| 147 | `lamina_n_area_retention` | `0` | 0 | nitrogen limit on retained green lamina area |
| 148 | `wsc_fate_law` | `1` | 1 | fate of the stem reserve |
| 149 | `wb_surface_limit` | `3` | 3 | surface limit of the water balance |
| 150 | `pheno_dev_factor_window` | `1` | 1 | window of the phenology development factor |
| 151 | `fe_trait_mode` | `1` | 1 | fruiting-efficiency trait mode |
| 152 | `census_files` | `1` | 1 | census diagnostic files on or off |
| 153 | `n_init_placement` | `1` | 1 | initial soil nitrogen placement |
| 154 | `organ_n_turnover_nema` | `1` | 1 | NEMA organ nitrogen turnover |
| 155 | `n_uptake_regulation` | `1` | 1 | regulation of nitrogen uptake |
| 156 | `n_construction_ceiling_curve` | `1` | 1 | construction-nitrogen ceiling curve |
| 157 | `lamina_n_storage` | `3` | 3 | lamina nitrogen storage |
| 158 | `nitrate_drain_rate_basis` | `1` | 1 | basis of the nitrate electron-drain term |
| 159 | `mobile_n_nema_estimate` | `0` | 0 | NEMA estimate of the mobile nitrogen pool |
| 160 | `justes_nmax_floor_basis` | `1` | 1 | floor basis of the Justes N_max curve |
| 161 | `canopy_clock_driver` | `0` | 0, 1 | driver of the canopy clock |
| 162 | `canopy_clock_basis` | `1` | 1 | basis of the canopy clock |
| 163 | `postanth_uptake_regulation` | `0` | 0, 1 | post-anthesis uptake regulation |
| 164 | `lamina_capacity_profile_basis` | `0` | 0, 1 | lamina capacity profile basis |
| 165 | `wsc_release_delay_basis` | `1` | 1 | delay basis of reserve release |
| 166 | `grain_ca_xylem_throughput` | `1` | 1 | grain calcium xylem throughput |
| 167 | `lamina_synthesis_sheath_basis` | `2` | 2 | sheath basis of lamina synthesis |
| 168 | `mineral_impedance_basis` | `1` | 1 | mineral impedance basis |
| 169 | `zn_uptake_kinetic_basis` | `1` | 1 | zinc uptake kinetic basis |
| 170 | `RESERVED` | `0` | 0 (retired, unavailable in v1.0; read, no effect) | RESERVED (diffusion temperature; not built) |
| 171 | `mineral_n_coupling_basis` | `1` | 1 | mineral-nitrogen coupling basis |
| 172 | `stem_reserve_deposition_n_basis` | `0` | 0 | nitrogen basis of stem reserve deposition |
| 173 | `spac_sink_coupling` | `1` | 1 | SPAC-to-sink coupling |
| 174 | `canopy_fln_prediction` | `0` | 0 | final-leaf-number prediction |
| 175 | `grain_n_partition_basis` | `4` | 4 | grain nitrogen partition basis |
| 176 | `vegetative_n_demand_basis` | `0` | 0 | vegetative nitrogen demand basis |
| 177 | `grain_n_potential_trait` | `0` | 0 | grain nitrogen potential as a cultivar trait |
| 178 | `chaff_structural_n_basis` | `0` | 0, 1, 2 | measured structural-nitrogen floor of the chaff |
| 179 | `lamina_expansion_n_limit` | `1` | 1 | nitrogen limit on the lamina expansion increment |
| 180 | `grain_ca_loading_basis` | `1` | 1 | grain calcium loading basis |
| 181 | `mg_solution_basis` | `0` | 0 | magnesium solution basis |
| 182 | `RESERVED` | `0` | 0 (retired, unavailable in v1.0; read, no effect) | RESERVED (vegetative-N demand basis; not built) |
| 183 | `mineral_buffer_power_basis` | `1` | 1 | soil buffer power of the mineral diffusion term |

### 2.3 Notes on individual records

* **Record 6** (`maximum LAI`) is read into `gParam_lai_max` and a full-source search finds
  no other reference to that name, so it reaches no state variable in v1.0. It is kept for
  positional stability. It is not a RESERVED row: it is a live read with no consumer.
* **Record 7** (`specific leaf area`) is inert while record 69 = 1, the shipped value: the
  area-to-mass conversion then runs on the organ-level canopy's own leaf mass per area and
  the realized specific leaf area becomes an *output*, `SLA_out_m2_kg` in
  `output/daily_growth.csv`. It is live only with record 69 = 0, which v1.0 does not accept.
* **Records 14 and 15** (Jarvis `a1`, `a2`) are read and declared and referenced nowhere
  else; the Jarvis stomatal branch is not part of the v1.0 stack. Do not list them as
  sampled parameters in a sensitivity analysis.
* **Records 73 and 74** (measured sub-daily wind file and its measurement height) are read
  only when record 71 >= 1. Record 71 is restricted to 0 in v1.0, so record 73 may be blank
  and record 74 may be `-1.0`, as in the reference deck.
* **Record 81** carries **two** values on one line: the SPAC mode and its leaf water
  potential threshold in MPa (`3 -1.5` in the reference deck).
* **Records 119 and 120** are RESERVED and are checked: a value other than 0 stops the run
  with `ERROR rec 119/120: RESERVED rows must be 0`.
* **Records 170 and 182** are RESERVED rows read as integers with an end-of-file default of
  0 and no range check; they have no effect anywhere and exist to hold the positions of
  records 171-172 and 183.
* **Record 39** is also settable from the command line. The validation block runs *after*
  that override, so a retired energy-balance mode cannot enter that way either.

### 2.4 Two settings of `input/phenology_cultivar.txt` are restricted in the same way

| setting | accepted in v1.0 | note |
|---|---|---|
| `GDD_Method` | 1 | `calc_dtt` implements method 1 only; the code of every other method was removed |
| `Use_Beta_Emerge` | 0, 1 | 0 is the v1.0 setting; 1 selects the beta-function emergence model of Part 1 Table 6, built and inactive. Any other integer is refused rather than read as "not 1" |

### 2.5 Optional diagnostics that v1.0 cannot produce

Seven driver settings would have selected an optional diagnostic file whose writer the
clean-up removed, and all seven are now refused: **record 99 = 2, records 134, 135, 136 and
137 at 1 or above, and record 139 modes 2 and 3**. The diagnostic files they would have
written — `daily_tiller_afford.csv`, `daily_n_afford_priority.csv`,
`daily_n_payment_luxury_pass2.csv`, `daily_n_afford_standing_basis.csv`,
`daily_floret_n.csv`, `daily_n_uptake_bound.csv` and `daily_n_uptake_mode3.csv` — are
**unavailable in v1.0**. Eight further conditional diagnostics are likewise unreachable
because the only accepted value of their record leaves them off.

The full list, with the record, the variable name and the file name of each, verified
against the code, is `docs/V1_KNOWN_LIMITATIONS.md`. None of this affects the 88
output files a v1.0 run writes.

---

## 3. The daily weather file (named by driver record 1)

CSV, one header row then one row per simulated day, in calendar order. The number of data
rows sets the length of the simulation. Twelve columns, read positionally:

| col | quantity | unit |
|---:|---|---|
| 1 | year | — |
| 2 | month | — |
| 3 | day of month | — |
| 4 | sunshine hours | h; may be empty, as in the reference deck |
| 5 | maximum air temperature | C |
| 6 | minimum air temperature | C |
| 7 | maximum relative humidity | % |
| 8 | minimum relative humidity | % |
| 9 | wind speed at the reference height (driver record 11) | m s-1 |
| 10 | rainfall | mm d-1 |
| 11 | irrigation | mm d-1 |
| 12 | global radiation | MJ m-2 d-1 |

The same file drives phenology. Do not create a separate weather file for it.

---

## 4. `input/driver_swtr.txt` — soil profile setup

Positional, one value per line, `!` comments after the value.

| # | entry | unit / form |
|---:|---|---|
| 1 | a file-path field, read as text. The shipped decks repeat the weather file path here; the soil hydraulic file is **not** taken from this line but from the fixed name `input/SoilLayer_vgm.csv` | path |
| 2 | `NTsteps`, the number of Richards internal time steps | count |
| 3 | initial root depth | m |
| 4 | `wt_pres`, water-table presence | Fortran logical, `.TRUE.` / `.FALSE.` |
| 5 | `NSoil_layers`, the number of soil layers | count; 10 in the reference deck |
| 6..5+N | initial volumetric water content, one line per layer, top first | m3 m-3 |
| then | layer depth boundaries, one line each | m, negative downward |

---

## 5. `input/SoilLayer_vgm.csv` — layer hydraulics

CSV with this exact header, then one row per soil layer, top first:

```
Layer,theta_r,theta_s,alpha_m-1,n,l,ksat_m_d
```

| column | quantity | unit |
|---|---|---|
| `Layer` | layer index, 1 = top | — |
| `theta_r` | residual water content | m3 m-3 |
| `theta_s` | saturated water content | m3 m-3 |
| `alpha_m-1` | van Genuchten alpha | m-1 |
| `n` | van Genuchten n | — |
| `l` | Mualem pore-connectivity exponent | — |
| `ksat_m_d` | saturated hydraulic conductivity | m d-1 |

This is the single supported source of soil hydraulics; v1.0 has no pedotransfer path
from soil texture. The parameters actually in force are echoed to
`output/vgm_params.csv`, including the derived `m`, so a run always documents the
hydraulics it used.

---

## 6. `input/phenology_cultivar.txt` — cultivar phenology and early growth

Positional: one numeric value per line, followed by whitespace and a label that the reader
ignores. Rows 1-35 are mandatory. Rows 36 onward are optional and have compiled defaults;
they are read in order, so a file cannot supply a later optional row without the earlier
ones.

| # | entry | unit |
|---:|---|---|
| 1 | `GDD_Method` | must be 1 (section 2.4) |
| 2-5 | `T_Base`, `T_Opt_Lower`, `T_Opt_Upper`, `T_Max` | C |
| 6 | `Photo_Sens`, photoperiod sensitivity | — |
| 7 | `Vern_Factor`, vernalization sensitivity | — |
| 8-15 | the eight stage thresholds in cumulative thermal development units: Emergence, 3_leaves, Tillering, Booting, Heading, Anthesis, Grain_filling, Maturity | degCd |
| 16-18 | sowing year, month, day | — |
| 19-27 | the nine observed stage anchors in days after planting: emergence, 2 leaves, 3 leaves, tillering, booting, heading, anthesis, grain filling, maturity | d |
| 28-31 | bootstrap organ masses per plant: green, stem, root, storage | g plant-1 |
| 32-35 | initial organ nitrogen concentrations: leaf, stem, root, storage | g N g-1 DM |
| 36 | `Use_Beta_Emerge` | 0 or 1 (section 2.4) |
| 37+ | the Wang et al. beta-emergence parameters, read only when `Use_Beta_Emerge` = 1 | — |

### Optional cultivar trait rows, appended at the end of the file

Each has an end-of-file default: the row absent means the compiled value is kept and the
run is byte-identical. They are read in the order given, so a deck needing only the last
one must still carry the earlier ones at their defaults.

| row | trait | unit | default | source of the default |
|---|---|---|---|---|
| `Phyllochron_degCd` | cultivar phyllochron, air basis | degCd | **73.9** | McMaster et al. 2003, Ann Bot 91:697-705, Table 2 p.701 (DOI 10.1093/aob/mcg074): mean of the two season means 73.0 and 74.8. Measured |
| `Blade_L_max_cm` | cap on the blade-length-by-rank profile | cm | **20.8** | the profile's own maximum; Abichou 2016 Ch.3 Fig. 8a p.104, NI series, median of five cultivars. Literature, figure-level, provisional |
| `Blade_W_max_cm` | cap on the blade-width-by-rank profile | cm | **1.45** | the profile's own maximum; Abichou 2016 Ch.2 Fig. 3e p.67, central of the 15-treatment Soissons panel. Literature, figure-level, declared gap |
| `Sheath_L_max_cm` | cap on the sheath-length-by-rank profile | cm | **17.0** | the profile's own maximum; Abichou 2016 Ch.2 Fig. 3b p.67, medians of the 15-treatment Soissons panel; range 14-20 cm. Literature |
| `Sheath_D_cm` | mean diameter of the upper sheaths | cm | **0.45** | Abichou 2016 Ch.2 Table S2 p.87, median of eight printed values; measured bracket 0.36-0.64 cm. Literature |
| `Axes_max_per_plant` | cultivar tillering capacity, main stem plus tillers carried alive at once | axes plant-1 | **0.0** = no cap | the crop-scale bound is this value times driver record 8 (plant density). Admissible band from Abichou 2016 Ch.3 Fig. 3 p.100, where live-axis density plateaus at 600-800 axes m-2, divided by the user's own sowing density |

At a default value the reference profile array is copied verbatim rather than multiplied
by 1.0 — an exact-equality branch, because IEEE arithmetic does not guarantee that a
multiply by one reproduces the original bit for bit. A non-default value scales every
element of the profile proportionally.

---

## 7. `input/cultivar_yield.csv` — cultivar yield coefficients

Positional, **value first**, then name, unit and source on the same comma-separated line;
`!` lines are skipped. **Accepted lengths: 6, 8 or 10 values** — each optional pair is
both or neither, and a file with 7 or 9 values is refused with an `ERROR STOP` naming the
first missing field.

| # | value | unit | note |
|---:|---|---|---|
| 1 | `fruiting_efficiency` | grains g-1 ear DW | grain number = this times ear dry weight at ear set |
| 2 | `TKW_potential` | g per 1000 kernels | the potential kernel mass. Under driver record 118 = 1, the shipped value, this must be the **measured** cultivar final kernel mass: it is a measured input, not a calibration parameter |
| 3 | `chaff_share` | fraction | share of ear-window reproductive dry matter going to the ear and chaff structure |
| 4 | `ear_end_dvs2` | dvs 0-2 | end of ear construction and grain-number set; 1.0 = anthesis in every shipped deck |
| 5 | `storage_onset_dvs2` | dvs 0-2 | start of the spike-growth partition window |
| 6 | `storage_peak_frac` | fraction | peak of that window |
| 7-8 | `TKW_env_min`, `TKW_env_max` | g per 1000 kernels | optional pair: the record-125 cultivar kernel-mass envelope. When present, a `TKW_potential` outside `[min, max]` is refused |
| 9-10 | `gfb_tbase_C`, `gfb_dfill_Cd` | C, and degCd **above that base** | optional pair: the grain-filling clock as a cultivar trait. Compiled defaults **6.618733 C** and **405.58109 degCd**, so a deck that omits the pair is byte-identical |

`gfb_dfill_Cd` is a thermal duration **above `gfb_tbase_C`**, not a base-0 thermal time.
The reader checks only that each value is a readable, finite, strictly positive number, so
entering a base-0 literature duration here without converting it is a user error the model
will not catch. Every run prints the clock in force and where each of its two numbers came
from, so check that line.

---

## 8. `input/driver_N.csv` — soil and plant nitrogen parameters

Positional: 70 values, one per non-comment line, each line carrying its own citation as a
comment. The groups are the root/shoot nitrate partition (`f_shoot_min`, `f_shoot_max`,
`K_partition_NO3`, after Andrews 1986 and Andrews et al. 1992), the root nitrate-reductase
capacity `r_root_NR` (g N g-1 root DM d-1), the grain structural nitrogen `n_stru_grain`
(g N g-1 DM) and the end of endosperm cell division `dvs_grain_stru_end` (after Martre et
al. 2003), the shoot reduction ceiling `N_max_ref` (umol m-2 s-1, Busch 2017), the
initial mineral nitrogen `N_soil_init` (nitrate) and `N_soil_init_NH4`, the urea
hydrolysis and nitrification kinetics (`k_hyd_ref`, `Q10_hyd`, `k_nit_ref`, `Q10_nit`,
`T_ref_hyd`, `T_ref_nit`), and the residue and microbial kinetics (`residue_C_frac`,
`k_fast_ref`, `k_slow_ref`, `Q10_dec`, `T_ref_dec`, `theta_ref_dec`, `e_mic`, `CN_mic`,
`k_mic`, and the `presow_*` pre-decomposition conditions).

Some positions are retained for **index stability only** and reach no state variable:
`k_N_canopy` (position 9) is not the canopy leaf-nitrogen extinction coefficient — that is
a compiled value pinned on Moreau et al. 2012 and never read from input, and the two
coefficients have near-identical names; `Ci_ref_N` and the `f_postanth_*` fractions belong
to mechanisms the model no longer uses. Do not list any of them as sampled parameters.
Do not delete them either: the file is positional.

---

## 9. `input/driver_minerals.csv` — the six elements

CSV; the header row is skipped by position but names the columns, then **one row per
element in the order Fe, Zn, Ca, Cu, Mg, Mn**. 21 columns:

```
name,modelled,pool_mg_kg,kd_L_kg,ceiling_mg_L,vmax_hats_nmol_gFW_h,km_hats_mg_L,
vmax_lats_nmol_gFW_h,km_lats_mg_L,ps_rate_nmol_gDW_d,ps_rate_def_nmol_gDW_d,
fe_def_switch,export_frac,xylem_route,grain_h2o_g_gDM,field_rate_mg_gDW_d,
grain_max_mg_kg,avail_A_mgL_CAL,xylem_Ca_mgL_CAL,buffer_power_b,grain_eN_mg_gN
```

Internal pools are mg element m-2 of ground; grain concentration is output as mg kg-1.
`grain_max_mg_kg` is a true guard on total grain mineral, **not a target**. `export_frac`
is the phloem mobility: calcium's is 0, which is why calcium reaches the grain only by the
xylem route. With driver record 131 = 1, the shipped value, a file carrying fewer than 18
numeric columns is refused rather than run with silently zeroed availabilities.

Several columns are inert in v1.0 and are kept for positional stability: `modelled`, the
HATS and LATS kinetic columns, the `km_*` columns, `fe_def_switch` and `ps_rate_def`.

---

## 10. `input/N_fertilizer.csv` — the fertilizer schedule

`!` comment lines, then the number of applications on the first non-comment line, then one
row per application:

```
DAS, N_kg_per_ha[, form]
```

`DAS` is days after sowing; `N_kg_per_ha` is the application rate. The **optional** third
field is a form keyword and defaults to `urea`, so existing two-column files keep working.

| keyword | (NO3, NH4, urea) fractions |
|---|---|
| `urea` | (0, 0, 1) |
| `AN` (ammonium nitrate) | (0.5, 0.5, 0) |
| `CAN` | (0.5, 0.5, 0) |
| `AS` (ammonium sulfate) | (0, 1, 0) |
| `KNO3` or `nitrate` | (1, 0, 0) |

Fertilizer enters the urea, ammonium and nitrate pools by those fractions; urea hydrolyzes
to ammonium and ammonium nitrifies to nitrate, both first-order and temperature-dependent
with `k(T) = k_ref Q10^((T_soil - T_ref)/10)`. **Only nitrate is leachable and taken up**;
ammonium is retained on the exchange complex. The kinetic constants live in
`input/driver_N.csv`.

---

## 11. `input/N_organic_treatments.csv` — residues and amendments (optional)

First non-comment line: the number of **active** treatments. **0 means no treatment and
reproduces the no-residue baseline**, as does an absent file. Then one row per event:

```
treatment_type, das, dry_mass_kg_ha, N_percent, C_percent, fast_fraction, age_days_before_sowing, notes
```

Units and precedence: `dry_mass_kg_ha * 0.1` = g DM m-2; `C_percent` and `N_percent` are
percent. `N_percent < 0` derives nitrogen from carbon and the type's C:N; `C_percent < 0`
uses `residue_C_frac` from `driver_N.csv`; `fast_fraction < 0` uses the type default.
`age_days_before_sowing > 0` makes the event **pre-sowing**: it is pre-decomposed at
initialization under `presow_Tsoil` and `presow_fwater`, with decomposition,
mineralization, immobilization, hydrolysis and nitrification but **no uptake and no
leaching**.

| `treatment_type` | default C:N | default fast fraction |
|---|---:|---:|
| `cereal_residue` | 80 | 0.20 |
| `immature_cereal` | 40 | 0.35 |
| `green_manure` | 20 | 0.55 |
| `legume_residue` | 20 | 0.55 |
| `manure` | 12 | 0.40 |
| `compost` | 18 | 0.10 |
| `custom` | — | — |

Diagnostics go to `output/daily_organic_N.csv`, including the soil-nitrogen closure.

---

## 12. `input/tsoil_init.csv` — measured initial soil temperature (optional)

CSV of `z_m,T_C` anchors, one per line, `!` comments and blank lines allowed. Depth is
positive downward in metres, temperature in C. The profile is initialized by piecewise
interpolation between the anchors and held uniform below the deepest one; with the file
absent, the driver record-50 initialization is used (0 = uniform, 1 = damped, with driver
records 49 `Tdeep` and 50).

The file is optional, but **once present it must be usable**: the run stops with a message
naming the offending line if a non-comment line is not a pair of numbers, if there are
fewer than two or more than eight anchors, if the depths do not increase strictly down the
file, or if a depth falls outside 0-50 m or a temperature outside -60 to 80 C (which is
what a `-999` fill value or a depth written in cm looks like). A malformed file is never
partly read and never quietly replaced by the computed profile of record 50, because the
point of supplying it is that the initial condition is measured. To use record 50 instead,
remove the file.

---

## 13. `input/cultivar_partition.csv` — leaf partition window (optional)

Read if present; if absent the run prints
`NOTE: input/cultivar_partition.csv absent; leaf window = WOFOST baseline defaults.` and
uses those compiled defaults. No shipped deck supplies it.

---

## 14. Reference data that the decks carry but the model never reads

`testcase/input/` also holds measurement and provenance files that are byte-copied with
the deck so that the directory is exactly what the frozen run was fed. They change
nothing: `driver_main.template`, `driver_main_silledBottom.txt`, the four alternative
`SoilLayer_vgm_*.csv` profiles, the orphan `driver_resp_mainten.txt` (read by no compiled
unit of v1.0), `site_observed_phenology.csv` and `site_observed_targets.csv`. See
`docs/FILE_INVENTORY.md` section 2.2. (The WheatDryFACE measurement tables that the frozen
run directories also held are not read by the model and are not shipped here; they are the
separate observation deposit cited in `README.md`.)

---

## 15. Where to look next

| question | file |
|---|---|
| how do I build and run the model? | `README.md` |
| what do the 88 output files contain? | `docs/FILE_INVENTORY.md` section 3 |
| why was a driver value refused? | this file, section 2.2 |
| which optional diagnostics are gone? | `docs/V1_KNOWN_LIMITATIONS.md` |
| what does each source file do? | `docs/MODULES.md` |
