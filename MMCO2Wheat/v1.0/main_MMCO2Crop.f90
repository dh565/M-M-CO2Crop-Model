!===============================================================================!
! M&M-CO2Wheat v1.0  |  main_MMCO2Crop.f90
!
! Purpose: The program. It reads the 183 positional records of
!          input/driver_main.txt and the deck input files (weather and driver
!          records, soil water and van Genuchten-Mualem layers, soil nitrogen,
!          fertilizer and organic-residue schedules, cultivar phenology, yield
!          and mineral parameters), initializes the soil water, heat, nitrogen
!          and mineral states, then runs the DAILY loop for the whole season:
!          phenology, the organ-level canopy chronology, the hourly canopy
!          energy balance and gas exchange iterated against the Richards soil-
!          water solver, interception, infiltration and ponding, soil heat, the
!          soil and plant nitrogen economy, carbon allocation and yield
!          formation, the mineral economy, and the daily water, energy, carbon
!          and nitrogen closures. It writes the 88 output files of output/ and
!          the season summary to run.log. Time steps: the Richards adaptive
!          sub-step, the hourly canopy and energy-balance step, the daily crop
!          and soil step, and the season. Spatial units: the sunlit and shaded
!          leaf class, the leaf rank and tiller axis, the plant organ, the soil
!          layer, and the point as a whole.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: SHARED (all parameters, switches and state); TYPES, ONLY: wp, ik;
!       csv_file; CANOPY_MODEL (compute_daily_canopy_fluxes, canopy_n0_top);
!       energy_balance_module (repartition_day_actual and the hourly writers);
!       soil_heat_module (initialize_soil_heat); GROWTH_MODULE (allocation,
!       structure, senescence, the organ-level canopy chain, grain number and
!       filling, organ nitrogen turnover, the internode height law);
!       MINERAL_MODULE; SOILN_MODULE (decompose_residue_day); PHENOLOGY_MODULE;
!       and the external subroutine hydro through an explicit interface.
! Output: the 88 files of output/ - the daily and hourly energy-balance,
!       radiation, gas-exchange and leaf-temperature files, the soil profile,
!       water-balance and Feddes diagnostics, the soil- and plant-nitrogen,
!       mineral, canopy-architecture, reserve, grain-number, grain-filling and
!       organ-nitrogen files, the season tables, and run.log.
! Called by: nothing - this is the executable's main program (mmco2crop).
!===============================================================================!

PROGRAM mmco2crop



!------------------------------------------------------------------------------!
USE SHARED   ! Shared variables.
USE csv_file ! Enable reading & writing csv files.
USE TYPES,   ONLY: wp, ik
USE CANOPY_MODEL, ONLY: compute_daily_canopy_fluxes, canopy_n0_top   ! rec 164: canopy_n0_top is the SAME PURE profile function that produces the n_top_g_m2 column of output/daily_canopy_n_profile.csv; it is used READ-ONLY and is not altered
USE energy_balance_module, ONLY: repartition_day_actual, soil_heat_resid_max, &
    write_hourly_eb_actual, write_hourly_radiation_committed
USE energy_balance_module, ONLY: write_hourly_tleaf_classes   ! (rec 95)
USE soil_heat_module, ONLY: initialize_soil_heat
USE GROWTH_MODULE, ONLY: initialize_growth_state, bootstrap_emergence_state, bootstrap_nitrogen_state, &
  allocate_daily_biomass, update_root_depth, update_structure, &
  apply_leaf_senescence, allocate_grain_N_source_sink_daily, &
  allocate_yield_components_daily, sla_from_dvs, partition_fractions, &
  apply_organN_turnover_daily, &  ! flag grainN_partition (rec 63)
  update_canopy_leaf_ranks, &     ! flag canopy_leaf_ranks (rec 64)
  update_canopy_blade_area, &     ! flag canopy_blade_area (rec 65)
  cba_L_at, &
  update_canopy_tillering, &      ! flag canopy_tillering (rec 66)
  ct_nl_tiller, ct_L_tiller, ct_W_tiller, ct_hs_axis, &
  update_canopy_leaf_senesc, &
  update_canopy_sheath_ear, &     ! flag canopy_sheath_ear (rec 68)
  cse_Ls_at, cse_Li_at, cse_Ds_ms, cse_ko_tiller, &
  update_canopy_sink_state, &    ! flag canopy_sink_lai (rec 69)
  update_leaf_expansion_water, & ! flag leaf_expansion_water (rec 83)
  update_leaf_carbon_bonus, &    ! rec 144 carbon_coupled_leaf_expansion
  update_lamina_n_arrest_prov, & ! rec 122 lai_n_arrest
  update_spike_axes, update_spike_window, & ! flag grain_number_mech (rec 89)
  bgn_set_spike_demand, &                  ! record 98, flag spike_sink_mech (rec 98)
  sdm4_set_demand, &                       ! rec 100 MODE 4, SHARED (recs 103/104)
  inl_advance                              ! rec 109 internode height law, SHARED (recs 110-113)
USE MINERAL_MODULE, ONLY: NMIN, mineral_names, read_mineral_params, init_mineral_pools, &
  update_minerals_daily, M_soil_labile, M_root, M_leaf, M_stem, M_chaff, M_grain, &
  M_uptake_mf, M_uptake_ac, M_remob_grain, M_balance_resid, M_cap_binding, &
  mm_grain_max, M_cap_refused, M_cap_days, &   ! rec 126
  M_soln, mm_xylem_ca_cal, &       ! rec 131
  M_uptake_diff, M_uptake_hats, msm_s_mf_cum, msm_s_diff_cum, msm_s_hats_cum, &   ! rec 140
  mm_grain_eN, mnc_refused, mnc_days, &   ! rec 171 reporting
  mm_kd, mm_ceiling, mm_pool_mgkg, mm_modelled, mm_b   ! rec 183 reporting ONLY: the column-3 Kd, column-4 ceiling, column-2 pool, modeled flag and the OFF column-19 constant b_e replaces; READ-ONLY in this file
USE SOILN_MODULE, ONLY: decompose_residue_day
USE PHENOLOGY_MODULE, ONLY: PhenologyParams, PhenologyState, NUM_PHENO_STAGES, &
  initialize_phenology_params, read_phenology_cultivar, initialize_phenology_state, &
  bootstrap_phenology_state, update_phenology_day, calendar_doy, calc_dev_stage_0_2, calc_dev_stage_rate_0_2, &
  calc_dtt, TS_TDU_MCMASTER   ! rec 161 canopy_clock_driver (SHARED.f90): the SAME degree-day law and the SAME terminal-spikelet threshold, not re-implemented here
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
IMPLICIT NONE
INTERFACE
  SUBROUTINE hydro (theta_ly, day_id, t_end)
    USE SHARED, ONLY: NSoil_layers
    USE TYPES, ONLY: wp
    REAL(wp), INTENT(INOUT) :: theta_ly(NSoil_layers)
    INTEGER, INTENT(IN) :: day_id
    REAL(wp), INTENT(IN), OPTIONAL :: t_end
  END SUBROUTINE hydro
END INTERFACE
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
REAL(wp) :: avail_water_mm
REAL(wp) :: avail_water
REAL(wp) :: runoff_mm_d
REAL(wp) :: pond_evap_mm
REAL(wp) :: infil_from_supply_mm_d, infil_from_pond_mm_d
REAL(wp) :: ponding_start_day, delta_ponding_mm, total_pond_evap
REAL(wp) :: ponding_min_day, ponding_max_day, ponding_mean_day, ponding_sum_day
REAL(wp) :: shortfall_mm
REAL(wp) :: irr
REAL(wp) :: irr_eff
REAL(wp) :: water_supply_mm_d
REAL(wp) :: soil_mb_resid_mm
REAL(wp) :: prestart_precip_mm
REAL(wp) :: bulk_density_g_cm3
REAL(wp) :: soil_storage_mm
REAL(wp) :: hydro_mb_resid_mm
REAL(wp) :: hydro_mb_correction_mm  ! rec 149 wb_surface_limit=2 ON-only (mirrors hydro_mb_resid_mm's plumbing)
REAL(wp) :: t_i_day
REAL(wp) :: e_i_day
REAL(wp) :: trans_act_pot
! Driver record 78: the potential of the canopy pass whose hourly stores the
! actual-flux repartition scales. Kept separate from trans_act_pot so that the
! latter can hold the day UNSTRESSED potential for the water-stress factor and
! the water ledger while the energy repartition still pairs each pass actual
! flux with its own potential.
REAL(wp) :: trans_pot_pass_mm
REAL(wp) :: trans_pot_demandpass_mm   ! DIAGNOSTIC: SUM of this iteration's demand-pass series
REAL(wp) :: unstr_sum_resid_mm        ! DIAGNOSTIC: SUM(unstressed series) - trans_pot_unstressed_mm
REAL(wp) :: f_t_iter_in               ! DIAGNOSTIC: f_t the demand pass of this outer iteration used
! rec 173: this pass's RAW daily transpiration potential, captured before the outer-loop
! seed relaxation so the hourly demand series can be rescaled by the same factor.
REAL(wp) :: t_i_day_raw_ssc
INTEGER  :: unit_ft_outer             ! output/daily_ft_outer_diag.csv (both modes; read by nothing)
INTEGER  :: unit_feddes_node          ! output/daily_feddes_node_diag.csv (both modes; read by nothing)
REAL(wp) :: evapo_act_pot
REAL(wp) :: trans_act_used
REAL(wp) :: Ci_day
REAL(wp) :: LEc_day
REAL(wp) :: LEs_day
REAL(wp) :: Hc_day
REAL(wp) :: Hs_day
REAL(wp) :: rca_day, rcs_day, rsa_day, rss_day
REAL(wp) :: gsw_day_mean, rca_daylight, rcs_daylight, rsa_daylight, rss_daylight
REAL(wp) :: tr_daylight_mm, es_daylight_mm
REAL(wp) :: Ir_day
REAL(wp) :: Ir_mm_d
REAL(wp) :: drain_mm_d
REAL(wp) :: storage_init, storage_final, total_precip, total_irrig, total_et_out, total_drain, total_runoff
REAL(wp) :: total_supply_eff, total_intercept
REAL(wp) :: total_shortfall_mm
REAL(wp) :: wb_resid_mm_d, total_wb_resid
REAL(wp) :: total_evap_head_exfil, total_evap_atm_flux
REAL(wp) :: total_hydro_top_flux_m, total_hydro_root_sink_m, total_hydro_bottom_flux_m
REAL(wp) :: total_balance_soil_mm, total_balance_surface_mm, total_balance_mm
REAL(wp) :: total_soil_mb_resid_mm
REAL(wp) :: total_soil_evap, total_trans_act
REAL(wp) :: hydro_top_flux_m_day, hydro_top_infil_m_day, hydro_root_sink_m_day, hydro_bottom_flux_m_day
REAL(wp) :: hydro_top_req_m_day
REAL(wp) :: infil_mm_d
REAL(wp) :: pond_evap_mm_d
REAL(wp) :: soil_evap_mm_d
INTEGER, ALLOCATABLE :: hydro_fail_days(:)
INTEGER :: hydro_fail_count
! Richards failure counting. hydro_fail_count is incremented only on the
! ACCEPTED outer pass, while the Richards FAIL line of hydro.f90 prints on
! every outer pass, so the two numbers differ. This counter records the printed
! class itself, over all passes, and is reported beside the accepted-pass count.
! It reads the counter hydro.f90 already maintains (hydro_fail_events_call,
! reset at every entry to hydro) and changes nothing there.
INTEGER :: hydro_fail_events_all = 0
INTEGER :: hydro_fail_calls_all  = 0
REAL(wp) :: total_hydro_mb_resid
REAL(wp) :: total_hydro_mb_correction  ! rec 149 wb_surface_limit=2 ON-only (0 under modes 0/1)
REAL(wp) :: Ci_mean_total
REAL(wp) :: cum_tr, cum_e, cum_ir, cum_et, cum_et_ir, cum_precip
REAL(wp) :: storage_mm, delta_storage, delta_storage_total, balance_mm, balance_no_ir, prev_soil_storage
REAL(wp) :: balance_soil_mm, balance_surface_mm
REAL(wp) :: flux_top_mass_m_day, flux_top_mismatch_mm_d
REAL(wp) :: top_flux_net_mm_d, top_flux_up_mm_d
REAL(wp) :: storage_day1_end
REAL(wp) :: raw_soil_storage_mm
REAL(wp) :: hydro_mb_resid_day_m
REAL(wp) :: hydro_mb_correction_day_m  ! rec 149 wb_surface_limit=2 ON-only; same reset/accumulate pattern as hydro_mb_resid_day_m
INTEGER :: init_profile_mode
REAL(wp) :: init_bottom_head_m
REAL(wp), ALLOCATABLE :: vg_in(:,:)
REAL(wp) :: Epot_soil_mm
REAL(wp) :: pr_eff_base   ! base effective precip (post FAO-interception) before D6 drip add-back
REAL(wp) :: lai_growth_diag, height_growth_diag
REAL(wp) :: W_stem_hgt   ! rec 108: the argument of the plant-height allometry (= W_stem when rec 108 = 0)
REAL(wp) :: dvs2_growth, dvs2_rate_growth, rm_growth_day, rg_growth_day
REAL(wp) :: leafN_conc, stemN_conc, rootN_conc, storageN_conc
REAL(wp) :: agb_biomass, grain_yield_biomass, grain_protein_pct
REAL(wp) :: grainN_conc_b3, chaffN_conc_b3   ! B-3 whole-organ [N]
! Soil nitrogen module local variables.
! The sowing model day is a runtime value set from driver record 88
! (sowing_model_day); it anchors the fertilizer event, the residue event and
! the residue-evaporation mass.
REAL(wp) :: N_leach_d         ! daily N leaching loss (g N m-2 d-1)
REAL(wp) :: N_hydrolysis_day, N_nitrification_day  ! urea->NH4, NH4->NO3 daily fluxes (g N m-2 d-1)
REAL(wp) :: soilN_before_day, soilN_balance_resid  ! soil-N closure check (g N m-2; incl organic pools)
! Organic-residue decomposition daily diagnostics (g m-2 d-1) and pre-sowing closure.
REAL(wp) :: decomp_C_day, decomp_N_day, respired_C_day, net_min_NH4_day, immob_NH4_day, immob_NO3_day
REAL(wp) :: turnover_N_day, residue_N_in_day, pre_sow_N_resid
INTEGER  :: unit_organic, i_org
INTEGER  :: unit_humusN            ! output/daily_humus_N.csv
LOGICAL  :: humusN_open = .FALSE.  ! REPAIR: the file is ON-only (rec 93 >= 1)
INTEGER  :: unit_fnc               ! rec 97: output/daily_itemF_nclose.csv (ON-only)
LOGICAL  :: fnc_open = .FALSE.     ! rec 97: the file is ON-only (rec 97 >= 1), file-set clause
INTEGER  :: unit_ont               ! rec 154: output/daily_organ_n_turnover.csv (ON-only)
LOGICAL  :: ont_open = .FALSE.     ! rec 154: the file is ON-only (rec 154 >= 1), file-set clause
INTEGER  :: unit_npg               ! rec 121: output/daily_n_payment.csv (ON-only)
LOGICAL  :: npg_open = .FALSE.     ! rec 121: ON-only file (file-set clause)
INTEGER  :: unit_lna               ! rec 122: output/daily_lai_n_arrest.csv (ON-only)
LOGICAL  :: lna_open = .FALSE.     ! rec 122: ON-only file (file-set clause)
INTEGER  :: unit_wfl               ! rec 148: output/daily_wsc_fate.csv (ON-only)
LOGICAL  :: wfl_open = .FALSE.     ! rec 148: ON-only file (file-set clause); OFF writes nothing
INTEGER  :: unit_tna               ! rec 123: output/daily_tiller_n_afford.csv (ON-only)
LOGICAL  :: tna_open = .FALSE.     ! rec 123: ON-only file (file-set clause)
INTEGER  :: unit_rfl               ! rec 124: output/daily_root_front.csv (ON-only)
LOGICAL  :: rfl_open = .FALSE.     ! rec 124: ON-only file (file-set clause)
INTEGER  :: unit_kpt               ! rec 125: output/daily_kernel_potential.csv (ON-only)
LOGICAL  :: kpt_open = .FALSE.     ! rec 125: ON-only file (file-set clause)
INTEGER  :: unit_rsf               ! rec 127: output/daily_repro_surplus_fate.csv (ON-only)
LOGICAL  :: rsf_open = .FALSE.     ! rec 127: ON-only file
INTEGER  :: unit_ncs               ! rec 128: output/daily_n_construction.csv (ON-only)
CHARACTER(LEN=512) :: ncs_hdr_base ! rec 156: the incumbent daily_n_construction.csv header, so the two ON-only columns can be appended without restating it
LOGICAL  :: ncs_open = .FALSE.     ! rec 128: ON-only file
INTEGER  :: unit_cfo               ! rec 129: output/daily_carbon_overflow.csv (ON-only)
LOGICAL  :: cfo_open = .FALSE.     ! rec 129: ON-only file
INTEGER  :: unit_gnm               ! rec 130: output/daily_grain_n_floor.csv (ON-only)
LOGICAL  :: gnm_open = .FALSE.     ! rec 130: ON-only file
INTEGER  :: unit_nap               ! rec 134: output/daily_n_afford_priority.csv (ON-only)
LOGICAL  :: nap_open = .FALSE.     ! rec 134: ON-only file
INTEGER  :: unit_npl2              ! rec 135: output/daily_n_payment_luxury_pass2.csv (ON-only)
LOGICAL  :: npl2_open = .FALSE.    ! rec 135: ON-only file
INTEGER  :: unit_nsb               ! rec 136: output/daily_n_afford_standing_basis.csv (ON-only)
LOGICAL  :: nsb_open = .FALSE.     ! rec 136: ON-only file
INTEGER  :: unit_fnl               ! rec 137: output/daily_floret_n.csv (ON-only)
LOGICAL  :: fnl_open = .FALSE.     ! rec 137: ON-only file
REAL(wp) :: kpt_mgn_k = 0.0_wp     ! rec 125: grain N per kernel today, mg N kernel-1 (diagnostic)
INTEGER  :: rec119_reserved = 0, rec120_reserved = 0   ! RESERVED rows, must read 0 (INPUT_GUIDE 119-120)
INTEGER  :: unit_nud               ! rec 139: output/daily_n_uptake_terms.csv (ON-only)
LOGICAL  :: nud_open = .FALSE.     ! rec 139: ON-only file (file-set clause)
INTEGER  :: unit_nub               ! rec 139 mode 2: output/daily_n_uptake_bound.csv (ON-only, mode>=2)
LOGICAL  :: nub_open = .FALSE.     ! rec 139 mode 2: ON-only file, mode>=2 only
INTEGER  :: unit_nuc               ! rec 139 mode 3: output/daily_n_uptake_mode3.csv (ON-only, mode>=3)
LOGICAL  :: nuc_open = .FALSE.     ! rec 139 mode 3: ON-only file, mode>=3 only
INTEGER  :: unit_msm               ! rec 140: output/daily_mineral_uptake_terms.csv (ON-only)
LOGICAL  :: msm_open = .FALSE.     ! rec 140: ON-only file
INTEGER  :: i_msm                  ! rec 140: per-element loop index for output/daily_mineral_uptake_terms.csv
INTEGER  :: i_nud                  ! rec 139: per-layer loop index for output/daily_n_uptake_terms.csv
INTEGER  :: unit_resevap, i_resv
INTEGER  :: unit_cres
INTEGER  :: unit_wsc          ! (rec 92): output/daily_stem_wsc.csv, ON-only
INTEGER  :: unit_wsl2         ! rec 105 = 2 (THE REPAIR): output/daily_wsc_mode2.csv, ON-only
INTEGER  :: unit_canE         ! (rec 94): output/daily_canopy_n_profile.csv, ON-only
INTEGER  :: unit_inl          ! rec 109: output/internode_ranks.csv, ON-only (rec 109 = 1)
INTEGER  :: unit_ins          ! rec 109: output/internode_shortfall.csv, ON-only (rec 109 = 1)
INTEGER  :: i_inl             ! rec 109: inversed-rank loop index for output/internode_ranks.csv
LOGICAL  :: inl_diag_open = .FALSE.   ! rec 109: header written? the OFF arm writes NO such file, so the byte gate is untouched
CHARACTER(LEN=6) :: inl_bind_txt      ! rec 109: the binding-limb enum as printed text (CLOSED/POT/MASS/CAP)
INTEGER  :: unit_gnf          ! (rec 86): output/daily_grain_n_sink.csv, ON-only
LOGICAL  :: gnf_diag_open = .FALSE.   ! header written? (the file starts at anthesis, not day 1)
REAL(wp) :: N_soil_uptake_day ! daily plant N uptake from soil (g N m-2 d-1)
REAL(wp) :: N_gain_total      ! total daily N gain = de-novo + soil (g N m-2 d-1)
REAL(wp) :: N_assim_biochem_retained_day ! biochemical N retained in plant pools (g N m-2 d-1)
REAL(wp) :: N_NO3_root, N_NO3_shoot      ! inorganic nitrate pools (g N m-2)
REAL(wp) :: N_assim_root_day, N_assim_shoot_day ! daily nitrate reduction by route (g N m-2 d-1)
REAL(wp) :: f_photoresp_diag             ! emergent shoot/total assim. fraction (-)
REAL(wp) :: N_plant_total     ! current total plant N (g N m-2)
REAL(wp) :: N_plant_before, N_plant_after, N_balance_resid
REAL(wp) :: N_fert_day, N_mineral_day
REAL(wp) :: N_remob_leaf_day, N_remob_stem_day
REAL(wp) :: AGB_kg_m2         ! above-ground biomass (t ha-1) for dilution curve
REAL(wp) :: Nc_pct            ! critical N concentration (%) Greenwood 1990
! Driver record 86: the shoot mass of the critical-nitrogen dilution curve, at
! the six places that read it - the chain AGB -> Nc -> g_nmax_shoot ->
! N_demand_day and the chain w_shoot_tha -> NNI -> f_Nstatus -> uptake. The sum
! is W_green + W_dead + W_dead_stem + W_stem + W_storage, and W_storage IS the
! reproductive pool, so the ears enter the sum exactly once; replacing
! W_storage by W_grain + W_chaff keeps them in exactly once on the
! grain-nitrogen-sink arms. These quantities feed STATE, not reporting. The
! dilution curve for wheat is defined on total shoot dry matter including the
! ears (Justes et al. 1994 Ann Bot 74:397-407).
LOGICAL  :: gnf_critN_basis = .FALSE.   ! rec 86 = 1 or 3 -> critical-N shoot mass on W_grain + W_chaff
REAL(wp) :: N_demand_day      ! daily plant N demand (g N m-2 d-1)
REAL(wp) :: f_Tr_N            ! transpiration-based N uptake scaling (-)
REAL(wp) :: f_CO2_N           ! CO2-dependent N uptake scaling (Gojon 2022) — DEPRECATED, Process 4
REAL(wp) :: N_uptake_massflow_day, N_uptake_active_day  ! Gojon B mass-flow vs active uptake (g N m-2 d-1)
REAL(wp) :: N_uptake_diffusion_day = 0.0_wp   ! rec 139: realized diffusion-term uptake, today (g N m-2 d-1)
! daily diagnostics (pre-competition potentials + uptake-law factors)
REAL(wp) :: S_active_pot_day = 0.0_wp, S_massflow_pot_day = 0.0_wp, U_pot_day = 0.0_wp
REAL(wp) :: S_diffusion_pot_day = 0.0_wp      ! rec 139: potential (pre-cap) diffusion-term supply, today (g N m-2 d-1)
REAL(wp) :: f_Nstatus_diag = 1.0_wp, a_root_diag = 1.0_wp
INTEGER  :: daa_root = 0               ! completed days after anthesis (dvs2>=1) for a_root (Zheng 2022 DAA basis)
LOGICAL  :: anth_crossed_L2c = .FALSE. ! anthesis-day latch: the crossing day itself is DAA = 0
REAL(wp) :: W_storage_prev_day, dW_storage_day, N_assim_avail_pa  ! grain growth + post-anthesis N source (Process 5)
REAL(wp) :: dW_grain_day  ! daily grain (kernel) growth, sink-limited (Task 7)
REAL(wp) :: W_green_prev_day, W_stem_prev_day, W_root_prev_day, W_chaff_prev_day  ! organ growth capture (Task 8 minerals)
INTEGER :: unit_minerals, i_min
INTEGER :: unit_msl                ! rec 131: output/daily_mineral_supply.csv (ON-only)
LOGICAL :: msl_diag_open = .FALSE. ! rec 131: ON-only file
REAL(wp) :: grain_import_day, N_from_assim_day, N_surplus_veg_day ! Martre grain N source-sink diagnostics
INTEGER  :: I_fert             ! loop index for fertilizer check
REAL(wp) :: N_alloc_leaf_day, N_alloc_stem_day, N_alloc_root_day, N_alloc_storage_day
REAL(wp) :: An_growth_day
REAL(wp) :: T_leaf_growth_day
REAL(wp) :: t_leaf_growth_mean
INTEGER :: day_dt_halvings
REAL(wp) :: Epot_soil_demand_mm, trans_pot_demand_mm  ! P2-e: demand actually sent to Richards
INTEGER :: day_fail_events   ! P2-b: solver-failure events this day (all iterations)
LOGICAL :: day_fail_evap     ! P2-b: any failure under evaporative demand this day
LOGICAL :: day_min_dt
LOGICAL :: day_failed
CHARACTER(LEN=512)  :: swc_header
CHARACTER(LEN=2048) :: swc_line
CHARACTER(LEN=64)   :: swc_tmp
CHARACTER(LEN=1024) :: input_line
CHARACTER(LEN=512)  :: diurnal_header
CHARACTER(LEN=512)  :: root_header
INTEGER :: count_neg_hc
INTEGER :: count_neg_hs
INTEGER :: count_le_gt_rn
INTEGER :: count_et_gt_avail
INTEGER :: count_rn_neg
INTEGER :: count_le_neg
INTEGER :: unit_flux
INTEGER :: unit_swc
INTEGER :: unit_swc_raw
INTEGER :: unit_swc_mean
INTEGER :: unit_wbdiag
INTEGER :: unit_diurnal
INTEGER :: unit_rootdiag
INTEGER :: unit_rootcomp, unit_rootcomp_lyr   ! compensatory-uptake diagnostics
INTEGER :: unit_tsoilprof = -1                ! daily per-layer 1-D soil-heat profile (mode 2)
LOGICAL :: tsoilprof_open = .FALSE.
INTEGER :: unit_hydrodiag
INTEGER :: unit_hydroledger
INTEGER :: unit_soilprof
INTEGER :: unit_tleaf
INTEGER :: unit_resist
INTEGER :: unit_growth
INTEGER :: unit_yield
INTEGER :: unit_pheno
INTEGER :: unit_nitrogen
INTEGER :: unit_gas_exch
INTEGER :: unit_rad_aero
INTEGER :: unit_carbon_growth
! Driver record 87: an on-only diagnostic file. It records the leaf area index
! and the leaf-nitrogen pool that THIS day gas-exchange call actually read, so
! that the two cannot be joined from different days by a reader.
INTEGER :: unit_leafn_cap = -1
LOGICAL :: leafn_cap_open = .FALSE.
REAL(wp) :: N_leaf_before_lnc = 0.0_wp
! Driver record 64: the two diagnostic output files of the per-culm leaf-rank
! chronology. Both are opened only when the record is on.
INTEGER :: unit_clr
LOGICAL :: clr_diag_open = .FALSE.
! (flag canopy_blade_area, rec 65): same ON-only pattern.
INTEGER :: unit_cba
LOGICAL :: cba_diag_open = .FALSE.
! (flag canopy_tillering, rec 66): same ON-only pattern.
INTEGER :: unit_ct
LOGICAL :: ct_diag_open = .FALSE.
! (flag canopy_leaf_senesc, rec 67): same ON-only pattern.
INTEGER :: unit_cls
LOGICAL :: cls_diag_open = .FALSE.
! (flag canopy_sheath_ear, rec 68): same ON-only pattern.
INTEGER :: unit_cse
LOGICAL :: cse_diag_open = .FALSE.
! Driver record 69: the organ-law closure and limiter diagnostic. This step
! moves state, so its file is not an inert emission: it carries the lamina area
! demand, the carbon actually debited for it, the hard non-negative backlog,
! the demographic limiter counters, the realized specific leaf area, and the
! two leaf area indices kept explicitly apart.
INTEGER :: unit_csl
LOGICAL :: csl_diag_open = .FALSE.
! D2 (recs 82/83/84): one ON-only daily diagnostic + one end-of-run cohort table.
INTEGER :: unit_d2
LOGICAL :: d2_diag_open = .FALSE.
INTEGER :: unit_tca                      ! rec 99 = 2: output/daily_tiller_afford.csv
LOGICAL :: tca_diag_open = .FALSE.
INTEGER :: tca_days_written = 0          ! rows written to daily_tiller_afford.csv (one per stepped window day)
INTEGER :: i_tc                          ! implied-do index of the daily_tiller_afford.csv row
INTEGER :: unit_tca4                     ! rec 99 = 4: output/daily_tiller_afford4.csv
LOGICAL :: tca4_diag_open = .FALSE.
INTEGER :: tca4_days_written = 0         ! rows written to daily_tiller_afford4.csv (one per stepped window day)
LOGICAL :: tca4_endpoint_written = .FALSE.   ! the endpoint row (F = 1) written once
INTEGER :: phen_io
TYPE(PhenologyParams) :: phen_params
REAL(wp) :: cum_tdu_prev_gnx = 0.0_wp   ! (rec 89): previous day's cum_tdu, to form the day's thermal increment
TYPE(PhenologyState) :: phen_state
CHARACTER(LEN=256) :: phenology_cultivar_file
CHARACTER(LEN=20) :: phen_stage_name
REAL(wp) :: t_dev_pre, t_dev_post
! Driver records 161 and 162: clk_dtt_canopy is the increment handed to the
! CANOPY clock. At both records at their released values it is
! phen_state%current_dtt and nothing else.
REAL(wp) :: clk_dtt_canopy, clk_gate_canopy
REAL(wp) :: prev_height_state
INTEGER :: unit_substep
INTEGER :: ios
INTEGER :: imin_loop = 0   ! rec 171: element index for the ON-only season summary
!------------------------------------------------------------------------------!
! Output parameters (1D for now).
!------------------------------------------------------------------------------!
INTEGER, ALLOCATABLE :: pheno_stage (:)! phenological stage (0-N)
!------------------------------------------------------------------------------!
REAL(wp), ALLOCATABLE :: days_sim     (:)! days elapsed since simulation (days)
REAL(wp), ALLOCATABLE :: lai_sim      (:)! simulated LAI (m2 m-2)
REAL(wp), ALLOCATABLE :: wg_gleaves   (:)! weight of green leaves (g m-2)
REAL(wp), ALLOCATABLE :: wg_dleaves   (:)! weight of dead leaves (g m-2)
REAL(wp), ALLOCATABLE :: wg_Tleaves   (:)! weight of total leaves (g m-2)
REAL(wp), ALLOCATABLE :: wg_stem      (:)! weight of stem (g m-2)
REAL(wp), ALLOCATABLE :: wg_roots     (:)! weight of roots (g m-2)
REAL(wp), ALLOCATABLE :: wg_storage   (:)! weight of storage organs (g m-2)
REAL(wp), ALLOCATABLE :: plant_height (:)! plant height (m)
REAL(wp), ALLOCATABLE :: A_canopy     (:)! gross assim. by canopy (gCH2O m-2 day-1)
REAL(wp), ALLOCATABLE :: main_res     (:)! maintenance respiration (gCH2O m-2 day-1)
REAL(wp), ALLOCATABLE :: grow_res     (:)! growth respiration (gCH2O m-2 day-1)
REAL(wp), ALLOCATABLE :: soil_le      (:)! soil latent heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: plant_le     (:)! plant latent heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: tot_le       (:)! total latent heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: soil_h       (:)! soil sensible heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: plant_h      (:)! plant sensible heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: tot_h        (:)! total sensible heat flux (MJ m-2 day-1)
REAL(wp), ALLOCATABLE :: E_soil       (:)! actual soil evaporation (mm day-1)
REAL(wp), ALLOCATABLE :: E_plant      (:)! actual plant transpiration (mm day-1)
REAL(wp), ALLOCATABLE :: ET_a         (:)! actual ET (mm day-1)
REAL(wp), ALLOCATABLE :: swc_layers   (:)! soil water content (mm3 mm-3)
REAL(wp), ALLOCATABLE :: theta_raw    (:)! raw Richards theta before balance correction
!------------------------------------------------------------------------------!
! Output array (for saving results).
!------------------------------------------------------------------------------!
REAL(wp), ALLOCATABLE :: A_data     (:,:)
REAL(wp), ALLOCATABLE :: flux_daily (:,:)
REAL(wp), ALLOCATABLE :: flux_daily_full (:,:)
REAL(wp), ALLOCATABLE :: swc_daily  (:,:)
REAL(wp), ALLOCATABLE :: swc_layers_daily (:,:)
REAL(wp), ALLOCATABLE :: swc_layers_raw_daily (:,:)
REAL(wp), ALLOCATABLE :: swc_layers_mean_daily (:,:)
REAL(wp), ALLOCATABLE :: theta_start(:)
REAL(wp), ALLOCATABLE :: theta_init(:)
REAL(wp), ALLOCATABLE :: theta_day_sum(:)
REAL(wp), ALLOCATABLE :: theta_day0(:)
REAL(wp) :: t_day_sum
INTEGER :: outer_iter, outer_max
REAL(wp) :: tleaf_prev_out, tleaf_mean_out, tol_tleaf_out
REAL(wp) :: t_i_day_prev, e_i_day_prev
LOGICAL :: outer_converged, force_final, accum_diag, do_write
!------------------------------------------------------------------------------!
! Season counters for the canopy / Richards OUTER fixed point: the days on
! which the loop exhausted outer_max without reaching tol_tleaf_out, and the
! days on which it ran at all. Reported in one run.log line at season end; the
! per-day outer_converged and outer_resid_K columns are in
! output/daily_eb_actual.csv.
!------------------------------------------------------------------------------!
INTEGER :: outer_days_total = 0, outer_days_unconverged = 0
REAL(wp) :: outer_resid_max_season = -1.0_wp
INTEGER :: unit_eb_demand   ! A3b demand-state daily diagnostics
LOGICAL :: pass2_runs   ! post-Richards canopy pass is the single hourly-CSV writer
REAL(wp) :: Rn_demand_day, G_demand_day, LEc_demand_day, LEs_demand_day, Hc_demand_day, Hs_demand_day  ! A3b diagnostics
LOGICAL :: eb_mode_from_cli = .FALSE.   ! .TRUE. if energy_balance_mode came from a CLI arg
INTEGER :: unit_eb_act = -1             ! daily_eb_actual.csv; -1 = closed
LOGICAL :: eb_act_open = .FALSE.
INTEGER :: unit_tsoil
REAL(wp), PARAMETER :: ponding_head_min_mm = 1.0e-1_wp
!------------------------------------------------------------------------------!

! **************************************************************************** !
! 1. Read basic information from driver.txt file.
! **************************************************************************** !
! Open driver file:
OPEN (1,FILE='input/driver_main.txt',STATUS='OLD',IOSTAT=ios)
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not open input/driver_main.txt'
  STOP
END IF

! Read input file name.
READ (1,'(A)',IOSTAT=ios) file_name
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read file_name from driver_main.txt'
  STOP
END IF

! Detect no. of variables in input file from the first data line.
NinVar = detect_ninvar(file_name)
IF (NinVar < 1) THEN
  WRITE (*,*) 'ERROR: Could not detect number of input variables from: ', TRIM(file_name)
  STOP
END IF

! Read maxuimum no. of developmental stages of the plant.
READ (1,*,IOSTAT=ios) DVS_max
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read DVS_max from driver_main.txt'
  STOP
END IF

! Read latitude (1D for now).
READ (1,*,IOSTAT=ios) lat_s
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read lat_s from driver_main.txt'
  STOP
END IF
! Convert longitude to radians.
lat_s_rad = lat_s * (Param_pi / 180.0_wp)

READ (1,*,IOSTAT=ios) lon_s
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read lon_s from driver_main.txt'
  STOP
END IF
! Convert longitude to radians.
lon_s_rad = lon_s * (Param_pi / 180.0_wp)

! Read initial leaf area index.
READ (1,*,IOSTAT=ios) gParam_lai
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_lai from driver_main.txt'
  STOP
END IF

! Read maximum leaf area index.
READ (1,*,IOSTAT=ios) gParam_lai_max
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_lai_max from driver_main.txt'
  STOP
END IF

! Read specific leaf area (SLA; m2 g-1).
READ (1,*,IOSTAT=ios) gParam_sla
IF (ios /= 0) THEN
  gParam_sla = 0.020_wp
END IF

! Read plant density (plants m-2).
READ (1,*,IOSTAT=ios) plant_density
IF (ios /= 0) THEN
  plant_density = 200.0_wp
END IF

! Read leaf N fraction (g N g-1 DM).
READ (1,*,IOSTAT=ios) gParam_leafN_frac
IF (ios /= 0) THEN
  gParam_leafN_frac = 0.040_wp
END IF

! Read initial plant height.
READ (1,*,IOSTAT=ios) gParam_hgt
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_hgt from driver_main.txt'
  STOP
END IF

! Read reference height.
READ (1,*,IOSTAT=ios) gParam_refhgt
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_refhgt from driver_main.txt'
  STOP
END IF

! Read wind attenuation coefficient (alpha).
READ (1,*,IOSTAT=ios) a_wind
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read a_wind from driver_main.txt'
  STOP
END IF

! Read eddy attenuation coefficient (n).
READ (1,*,IOSTAT=ios) n_eddy
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read n_eddy from driver_main.txt'
  STOP
END IF

! Read Jarvis (1976) empirical coefficient a1.
READ (1,*,IOSTAT=ios) gParam_rstA1
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_rstA1 from driver_main.txt'
  STOP
END IF

! Read Jarvis (1976) empirical coefficient a2.
READ (1,*,IOSTAT=ios) gParam_rstA2
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_rstA2 from driver_main.txt'
  STOP
END IF

! Read mean leaf width (m).
READ (1,*,IOSTAT=ios) gParam_leafwidth
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_leafwidth from driver_main.txt'
  STOP
END IF

! Read upper layer thickness (m).
READ (1,*,IOSTAT=ios) dry_layer_thick
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read dry_layer_thick from driver_main.txt'
  STOP
END IF

! Read Porosity
READ (1,*,IOSTAT=ios) gParam_porosity
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_porosity from driver_main.txt'
  STOP
END IF

! Read pore size distance
READ (1,*,IOSTAT=ios) gParam_poredist
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_poredist from driver_main.txt'
  STOP
END IF

! Read theta_s
READ (1,*,IOSTAT=ios) gParam_vwcsat
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_vwcsat from driver_main.txt'
  STOP
END IF

! Read initail theta
READ (1,*,IOSTAT=ios) gParam_vwc01
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read gParam_vwc01 from driver_main.txt'
  STOP
END IF

! Read atmospheric CO2 (umol mol-1).
READ (1,*,IOSTAT=ios) Ca_s
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read Ca_s from driver_main.txt'
  STOP
END IF

! Optional: water stress mode (1=reduce gsw, 2=reduce An); default to 1 if absent.
READ (1,*,IOSTAT=ios) water_stress_mode
IF (ios /= 0) THEN
  water_stress_mode = 1
END IF
!------------------------------------------------------------------------------!
! Optional: soil water stress on/off (1=on, 0=off); default to on if absent.
! Re-using ios and dummy variable for simplicity
!------------------------------------------------------------------------------!
BLOCK
  INTEGER :: dummy_int
  READ (1,*,IOSTAT=ios) dummy_int
  IF (ios /= 0) THEN
    use_soil_water_stress = .TRUE.
  ELSE
    use_soil_water_stress = (dummy_int /= 0)
  END IF
  
  READ (1,*,IOSTAT=ios) dummy_int
  IF (ios /= 0) THEN
    allow_drain = .TRUE.
  ELSE
    allow_drain = (dummy_int /= 0)
  END IF
  drainage_type = MERGE(1, 0, allow_drain)

  READ (1,*,IOSTAT=ios) dummy_int
  IF (ios /= 0) THEN
    use_explicit_infiltration = .TRUE.
  ELSE
    use_explicit_infiltration = (dummy_int /= 0)
  END IF

  READ (1,*,IOSTAT=ios) dummy_int
  IF (ios /= 0) THEN
    integration_mode = 1
  ELSE
    integration_mode = dummy_int
  END IF

END BLOCK

! Optional: root growth rate (m day-1); default to 0.001 if absent.
READ (1,*,IOSTAT=ios) root_growth_rate
IF (ios /= 0) THEN
  root_growth_rate = 0.001_wp
END IF
! Optional: N_max (umol m-2 s-1); default to 1.21 if absent.
READ (1,*,IOSTAT=ios) N_max_param
IF (ios /= 0) THEN
  N_max_param = 1.21_wp
END IF
! Optional: enable dynamic n_eddy (0=off, 1=on); default off.
READ (1,*,IOSTAT=ios) use_dynamic_n_eddy
IF (ios /= 0) THEN
  use_dynamic_n_eddy = 0
END IF

! Optional: meteo debug prints (0=off, 1=on); default off.
BLOCK
  INTEGER :: dummy_int
  READ (1,*,IOSTAT=ios) dummy_int
  IF (ios /= 0) THEN
    meteo_debug = .FALSE.
  ELSE
    meteo_debug = (dummy_int /= 0)
  END IF
END BLOCK

! Optional: irrigation mode (1=above canopy, 2=below canopy); default to 2 if absent.
READ (1,*,IOSTAT=ios) irrig_mode
IF (ios /= 0) THEN
  irrig_mode = 2
END IF
! Optional: prestart precipitation (mm) applied once before day 1; default 0.
READ (1,*,IOSTAT=ios) prestart_precip_mm
IF (ios /= 0) THEN
  prestart_precip_mm = 0.0_wp
END IF
! Optional: bulk density (g/cm3) for porosity-based theta_s; default 1.65 for very sandy soils.
READ (1,*,IOSTAT=ios) bulk_density_g_cm3
IF (ios /= 0) THEN
  bulk_density_g_cm3 = 1.65_wp
END IF
! Optional: K averaging mode (0=arithmetic, 1=geometric, 2=upstream); default 0 if absent.
READ (1,*,IOSTAT=ios) k_avg_mode
IF (ios /= 0) THEN
  k_avg_mode = 1
END IF
! Optional: hydraulic model selector
! 0=VGM (van Genuchten-Mualem), 1=Brooks-Corey, 2=Campbell.
READ (1,*,IOSTAT=ios) hyd_model
IF (ios /= 0) THEN
  hyd_model = 0
END IF
WRITE (*,*) 'Hydraulic model:', hyd_model

! Optional: initialize theta profile mode
! 0=use driver_swtr theta as-is, 1=hydrostatic from bottom head,
! 2=monotonic profile from top theta to bottom theta (bottom head-derived).
READ (1,*,IOSTAT=ios) init_profile_mode
IF (ios /= 0) THEN
  init_profile_mode = 0
END IF
! Optional: bottom pressure head (m, negative for suction); used when init_profile_mode>0.
READ (1,*,IOSTAT=ios) init_bottom_head_m
IF (ios /= 0) THEN
  init_bottom_head_m = -0.35_wp
END IF
!------------------------------------------------------------------------------!
! Driver record 39, the energy-balance architecture. v1.0 implements 2 only:
! the three-node canopy / soil / airspace energy balance with the prognostic
! one-dimensional soil-heat profile supplying the ground heat flux. Every other
! value is refused below. The read is positional, so a driver that lacks the
! record hits end of file and keeps 0, which is then refused.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) energy_balance_mode
IF (ios /= 0) THEN
  energy_balance_mode = 0
END IF
!------------------------------------------------------------------------------!
! Driver record 40, the optional free-convection fallback of the three-node
! balance. End-of-file default: absent from the driver means the fallback is
! off and the surface layer is pure Monin-Obukhov. eb_zi_m and eb_gust_beta
! are explicit validated inputs, not calibration knobs.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) eb_freeconv_fallback
IF (ios /= 0) eb_freeconv_fallback = 0
READ (1,*,IOSTAT=ios) eb_zi_m
IF (ios /= 0) eb_zi_m = 1000.0_wp
READ (1,*,IOSTAT=ios) eb_gust_beta
IF (ios /= 0) eb_gust_beta = 1.2_wp
!------------------------------------------------------------------------------!
! Soil-heat parameters of the record-39 mode 2 profile, read positionally at
! the tail of the driver so that a driver lacking them hits end of file and
! keeps the SHARED defaults. The order is fixed and is documented in
! INPUT_GUIDE.md; the declarations, the defaults and the literature sources are
! in SHARED.f90.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) quartz_fraction
IF (ios /= 0) quartz_fraction = 0.60_wp
READ (1,*,IOSTAT=ios) Cv_mineral
IF (ios /= 0) Cv_mineral = 2.00e6_wp
READ (1,*,IOSTAT=ios) Cv_org
IF (ios /= 0) Cv_org = 2.50e6_wp
READ (1,*,IOSTAT=ios) Cv_water
IF (ios /= 0) Cv_water = 4.18e6_wp
READ (1,*,IOSTAT=ios) Cv_air
IF (ios /= 0) Cv_air = 1.25e3_wp
READ (1,*,IOSTAT=ios) bottom_bc_type
IF (ios /= 0) bottom_bc_type = 0
READ (1,*,IOSTAT=ios) Tdeep_C
IF (ios /= 0) Tdeep_C = 15.0_wp
READ (1,*,IOSTAT=ios) soil_T_init_method
IF (ios /= 0) soil_T_init_method = 1
! Compensatory root water uptake (Jarvis omega); EOF-defaulted to 1.0 (OFF) for back-compat.
READ (1,*,IOSTAT=ios) root_comp_omega
IF (ios /= 0) root_comp_omega = 1.0_wp
IF (root_comp_omega <= 0.0_wp .OR. root_comp_omega > 1.0_wp) root_comp_omega = 1.0_wp
! Residue-cover evaporation suppression (records 52-54); EOF-defaulted OFF.
READ (1,*,IOSTAT=ios) residue_evap_mode
IF (ios /= 0) residue_evap_mode = 0
IF (residue_evap_mode /= 0 .AND. residue_evap_mode /= 1) residue_evap_mode = 0
READ (1,*,IOSTAT=ios) residue_Am_ha_t
IF (ios /= 0) residue_Am_ha_t = 0.45_wp
IF (residue_Am_ha_t < 0.2_wp) residue_Am_ha_t = 0.2_wp
IF (residue_Am_ha_t > 0.8_wp) residue_Am_ha_t = 0.8_wp
READ (1,*,IOSTAT=ios) residue_fsupp
IF (ios /= 0) residue_fsupp = 0.75_wp
IF (residue_fsupp < 0.0_wp) residue_fsupp = 0.0_wp
IF (residue_fsupp > 0.9_wp) residue_fsupp = 0.9_wp
! Carbon-reserve remobilization (records 55-60); EOF-defaulted OFF.
READ (1,*,IOSTAT=ios) carbon_reserve_mode
IF (ios /= 0) carbon_reserve_mode = 0
IF (carbon_reserve_mode /= 0 .AND. carbon_reserve_mode /= 1) carbon_reserve_mode = 0
READ (1,*,IOSTAT=ios) stem_reserve_frac
IF (ios /= 0) stem_reserve_frac = 0.20_wp
IF (stem_reserve_frac < 0.0_wp) stem_reserve_frac = 0.0_wp
IF (stem_reserve_frac > 0.35_wp) stem_reserve_frac = 0.35_wp
READ (1,*,IOSTAT=ios) reserve_remob_rate
IF (ios /= 0) reserve_remob_rate = 0.05_wp
IF (reserve_remob_rate < 0.0_wp) reserve_remob_rate = 0.0_wp
IF (reserve_remob_rate > 0.15_wp) reserve_remob_rate = 0.15_wp
READ (1,*,IOSTAT=ios) reserve_start_dvs2
IF (ios /= 0) reserve_start_dvs2 = 1.00_wp
IF (reserve_start_dvs2 < 0.9_wp) reserve_start_dvs2 = 0.9_wp
IF (reserve_start_dvs2 > 1.3_wp) reserve_start_dvs2 = 1.3_wp
READ (1,*,IOSTAT=ios) reserve_end_dvs2
IF (ios /= 0) reserve_end_dvs2 = 2.00_wp
IF (reserve_end_dvs2 < 1.5_wp) reserve_end_dvs2 = 1.5_wp
IF (reserve_end_dvs2 > 2.0_wp) reserve_end_dvs2 = 2.0_wp
IF (reserve_end_dvs2 <= reserve_start_dvs2) THEN
  reserve_start_dvs2 = 1.00_wp; reserve_end_dvs2 = 2.00_wp
END IF
READ (1,*,IOSTAT=ios) reserve_conv_eff
IF (ios /= 0) reserve_conv_eff = 0.95_wp
IF (reserve_conv_eff < 0.70_wp) reserve_conv_eff = 0.70_wp
IF (reserve_conv_eff > 1.00_wp) reserve_conv_eff = 1.00_wp
!------------------------------------------------------------------------------!
! Driver record 61 stem_spike_realloc_delta: the fraction of the partition
! table's own stem share redirected to the spike from the existing spike-onset
! development stage. The at-anthesis above-ground stem share it is set against
! is 0.53 (Tyagi et al. 2004 Indian J Agron 49(3):163-167 Tables 2-3; Zhiipao
! et al. 2023 Front Plant Sci, doi:10.3389/fpls.2023.1087343 Table 2; Kunwar
! et al. 2026, doi:10.3389/fpls.2026.1740337). The released value is 0.0, which
! is the partition table exactly; the derivation is documented at
! growth_module.f90 partition_fractions.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) stem_spike_realloc_delta
IF (ios /= 0) stem_spike_realloc_delta = 0.0_wp
IF (stem_spike_realloc_delta < 0.0_wp) stem_spike_realloc_delta = 0.0_wp
IF (stem_spike_realloc_delta > 0.5_wp) stem_spike_realloc_delta = 0.5_wp
!------------------------------------------------------------------------------!
! Driver record 62 leafN_drawdown: the living-leaf nitrogen drawdown during
! grain filling together with the leaf-nitrogen capacity bound. Bounds fluxes
! only; no concentration output is capped. Constants and their sources:
! Bertheloot, Martre & Andrieu (2008) Plant Physiol 148:1707-1720, the
! first-order relative rate 0.0035 (degree-Cd)-1 on thermal time base 0 C and
! the stem residual 0.0035 g N g-1 DM p.1711; Silva-Perez et al. (2020) J Exp
! Bot 71(7):2299-2311 Fig. 4B p.2305, the holding capacity 3.8 g N m-2. The
! released value is 0: record 154 carries the post-anthesis organ-nitrogen
! turnover instead.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) leafN_drawdown
IF (ios /= 0) leafN_drawdown = 0
IF (leafN_drawdown /= 1) leafN_drawdown = 0
!------------------------------------------------------------------------------!
! Driver record 63 grainN_partition: the post-anthesis grain-nitrogen partition
! package - (a) stem, sheath and chaff senescence on the tissue-death rule
! d_tp = 0.4 of Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot
! 108:1097-1109 Table 2 and the mineral export it unlocks, (b) the organ-
! nitrogen turnover of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
! 35:781-796 with a uniform degradation rate 0.008 (degree-Cd)-1 on
! photosynthetic nitrogen and organ-type synthesis rates 1.5e-4 (lamina) and
! 3.0e-5 (internode and chaff), (c) the chaff nitrogen ceiling
! 0.017 g N g-1 DM, and (d) the chaff as an exporting source rather than a
! terminal sink. The released value is 0.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) grainN_partition
IF (ios /= 0) grainN_partition = 0
IF (grainN_partition /= 1) grainN_partition = 0
!------------------------------------------------------------------------------!
! Driver record 64 canopy_leaf_ranks: the per-culm leaf-rank chronology on the
! degree-Cd base-0 clock. The Haun stage follows the phyllochron of McMaster,
! Wilhelm, Palic, Porter & Jamieson (2003) Ann Bot 91:697-705 carrying the
! appearance-rate decline of Streck, Weiss & Baenziger (2003) Ann Bot
! 92:181-190; tip appearance and ligulation thermal times are resolved per rank
! (ligulation = appearance + the measured leaf-expansion duration of Bos 1999),
! and the emergent final leaf number is read off the anchored booting stage on
! the Zadoks ordering, with no final-leaf-number input anywhere. The alternative
! split phyllochron of Evers (2006) is computed side by side and never stacked
! with the rate decline; it is cited through secondary literature, not from that
! thesis directly. The released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_leaf_ranks
IF (ios /= 0) canopy_leaf_ranks = 0
IF (canopy_leaf_ranks /= 1) canopy_leaf_ranks = 0
!------------------------------------------------------------------------------!
! Driver record 65 canopy_blade_area: potential and realized one-sided blade
! area per main-stem rank. Blade LENGTH from Abichou (2016) PhD thesis
! (NNT tel-02316154) Ch.3 Fig. 8a p.104, median of the five cultivars; the
! order-4 polynomial of Ch.2 Table 2 Eqn T2.1 with the Table S2 coefficients
! is not usable, because at the printed precision its four terms of order 1e4
! cancel to a result of order 20 cm, so rounding alone spreads the length over
! -55.9 to +99.8 cm at rank 9. Blade WIDTH from Ch.2 Fig. 3e p.67, read at
! figure level, since the Eqn T2.4 coefficients are not printed in the thesis.
! The form factor is rank-resolved from Dornbusch, Watt, Baccar, Fournier &
! Andrieu (2011) Ann Bot 107:865-873, doi:10.1093/aob/mcq181, Fig. 3A for cv.
! Soissons, and not the flat 0.8, which is its adult-phase value only.
! Expansion runs over the leaf-elongation window between tip appearance and
! ligulation; the kinetic SHAPE is a linear ramp, which no source supplies.
! Tiller constants are carried from Abichou Ch.1 Table 4 p.71 and multiplied by
! no density - cohort densities arrive with record 66. The released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_blade_area
IF (ios /= 0) canopy_blade_area = 0
IF (canopy_blade_area /= 1) canopy_blade_area = 0
!------------------------------------------------------------------------------!
! Driver record 66 canopy_tillering: the axis (tiller) population, and with it
! the first ground-area quantity the canopy chain can produce, since record 65
! emits blade area per axis only. Cohorts T1 to T5 are born at the measured
! main-stem Haun stages of Abichou (2016) Ch.1 Table 6 p.31; each cohort
! occurrence probability comes from Evers et al. (2006) Ann Bot 97:649-658,
! doi:10.1093/aob/mcl020, Fig. 5A - 1 below a canopy-intercepted fraction of
! 0.40, falling to 0 at 0.55 - evaluated over that cohort's own window of
! opportunity, which the record-64 ligulation chronology supplies; the
! intercepted fraction is Beer on the model's own green leaf area index with the
! model's own effective PAR extinction. Cessation is the first of that light law
! and the empirical axis cap 700 m-2 (Abichou Ch.3 Fig. 3 p.100, envelope
! 600-800), and every run reports which one fired. Mortality follows Abichou
! Ch.1 Table 4: start at a main-stem Haun stage linear in the final leaf number
! (the coefficients are not printed and are reconstructed from the printed
! 5.5-8.0 and 10-14 ranges), end at flag-leaf ligulation, pattern Eqn T4.4 with
! theta = 3.32 and delta = 1.77, youngest cohorts culled first, which is that
! thesis's own rule. The final ear density is not taken from observation: it
! follows from the measured tiller die fraction. The released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_tillering
IF (ios /= 0) canopy_tillering = 0
IF (canopy_tillering /= 1) canopy_tillering = 0
!------------------------------------------------------------------------------!
! Driver record 67 canopy_leaf_senesc: per-rank blade senescence and the
! standing-dead pool. The green leaf number follows the four-phase law of
! Abichou (2016) Ch.2 Table 3 p.65, with phase 2 ending structurally at the
! ligulation of leaf Nl-4 and phase 3 at the ligulation of the flag leaf, both
! supplied by the record-64 chronology, so no absolute degree-Cd is
! transferred; the tiller equations T3.5 to T3.8 run on each axis's own Haun ramp
! with the maximum scaled by delta = 0.85 (Ch.2 Table 6 p.71, median of the
! measured values 0.79/0.85/0.96). The fitted g1, g2, g3 and the phase-4
! polynomial are from Ch.2 Table S1 p.86, median of the four main-stem rows at
! normal density with nitrogen, selected on site properties alone; g1 is
! transferred as a HAUN STAGE and not as the printed thermal time, which is
! phyllochron-specific. Blade area is retired rank by rank from the base using
! Ch.2 Eqn 3 p.59, senesced leaf number = Haun stage minus green leaf number,
! and dead area is retained STANDING, since no measured shedding law exists.
! Rank-progressive lifespan is structural rather than constant: Robert et al.
! (2018) Ann Bot 121:975-989, doi:10.1093/aob/mcx192, falsifies a constant
! per-rank lifespan. The onset conflicts with the Fournier (2003) onset quoted in an internal note
! value, which puts it later than Haun stage 6 while the nine measured Abichou
! treatments put it at 2.9-5.4; the read-at-source measurement is used and the
! disagreement is reported. The released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_leaf_senesc
IF (ios /= 0) canopy_leaf_senesc = 0
IF (canopy_leaf_senesc /= 1) canopy_leaf_senesc = 0
!------------------------------------------------------------------------------!
! Driver record 68 canopy_sheath_ear: the NON-LAMINA standing area, i.e. the
! green sheaths, the exposed stem and the ear plus exposed peduncle. This is
! the step that makes the model comparable like for like with a gap-fraction
! plant-area-index observation, so it emits two separated bases: the
! destructive GREEN basis of Abichou (2016) Ch.3 p.95, where the plant area
! index was measured at flowering by dissecting 45 plants per treatment and
! scanning all phytoelements, and the LI-2200C gap-fraction basis, which adds
! the standing dead pools and puts the cylindrical organs on the
! half-total-surface convention a gap-fraction inversion returns.
! Sheath LENGTH by inversed rank from Abichou Ch.2 Fig. 3b p.67, the same
! 15-treatment panel as the blade profile, with rank 1 anchored at the midpoint
! of the printed 14-20 cm largest-sheath range of p.66; the printed Eqn T2.2
! and Table S2 polynomial is not usable at scan resolution, exactly as for the
! blades. Sheath DIAMETER from Eqns T2.5a and T2.5b p.63 with the measured
! constant 0.124 cm (Ch.1 Table 4 p.71) and 0.45 cm, the median of the eight
! values printed in Table S2 p.87 (bracket 0.36-0.64); the break rank Nl-4 is
! the model's own first elongated internode, itself the median of the 17
! treatments in that table. Tillers follow Eqns T2.8 and T2.11 on the
! record-65 shifts with the glossary coefficients p.81.
! PERSISTENCE is the sequential n+2 rule of Abichou Annex 2 p.170, where the
! senescence of sheath n is simulated as synchronous with that of blade n+2;
! it is a model assumption, since that annex is a software option and not a
! field measurement, and the n+1 and n+3 arms are computed every day beside it.
! EXTENSION follows the same page sequential rule (blade first, then sheath;
! leaf extension starts 0.4 phyllochron after tip appearance and lasts 2
! phyllochrons), also a model assumption, with a no-lag upper-bound arm emitted
! beside it. The EXPOSED STEM is a result and not an input: every elongated
! internode is sheathed by the leaf below it, so only the uncovered length
! intercepts, with internode length from Ch.2 Fig. 3c p.67 and diameter
! 0.385 cm, the median of the ten printed Table S2 values; the printed Eqn T2.6
! shape returns a negative diameter at the flag rank and is not used.
! Two quantities have no measured source: no field-measured source for sheath-
! versus-blade senescence timing was located, so the n+2 rule stands with its
! sensitivity arms; and the thesis prints the ear schema at p.167 including the
! projected area of the awnless ear but no numeric value anywhere, so 15 cm2
! per fertile culm is carried as a bounded residual with its 10-20 bracket
! swept every day. The ear WINDOW is measured: 1.8 phyllochron from that axis
! own flag ligulation (Ch.3 p.96, 1333 degree-Cd, which p.100 places close to
! flowering). The released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_sheath_ear
IF (ios /= 0) canopy_sheath_ear = 0
IF (canopy_sheath_ear /= 1) canopy_sheath_ear = 0
!------------------------------------------------------------------------------!
! Driver record 69 canopy_sink_lai: the switch. This is the only step of the
! organ-canopy chain that moves state. The leaf area index stops being the
! bookkeeping product of green lamina mass and specific leaf area and becomes
! the SINK green lamina area index that the carbon-blind organ law of records
! 64 to 68 built; lamina carbon stops being a partition FRACTION and becomes an
! AREA DEMAND, the area increment times the leaf mass per area, with the
! remainder partitioned to stem, root and storage on the retained fractions
! renormalized over those three so that they still sum to exactly 1. An unmet
! demand is carried as a hard non-negative BACKLOG and the DEMOGRAPHIC limiter
! acts on future demand only: it suspends tiller initiation, aborts cohorts
! that have not yet produced a ligulated leaf, and cuts the emergent final leaf
! number by at most 2 (Evers 2006 PhD thesis, doi:10.18174/24302, Ch.3 Fig. 4,
! measured 9 -> 8). NO TIPPED LEAF IS EVER SCALED: in Bos (1999) Ch.2 the leaf
! elongation rate and the final width do not respond to a 2.6-fold range of
! photosynthetic photon flux density, so a carbon scaling of expansion is
! falsified and the response is demographic.
! The carbon and nitrogen balances still close at machine zero, the water
! balance and Richards are unaffected, the non-leaf fractions renormalize to
! exactly 1, and the backlog is hard non-negative and reported. Requires
! records 64 to 68. The leaf mass per area carries a source conflict: the
! bracket implied by the tillering thresholds of Evers Ch.4 (specific leaf area
! 30.8-50 m2 kg-1) does not overlap the measured realized specific leaf area of
! Bos et al. (2000) Ann Bot 85:243-249 Table 1 (20-28) and Araus, Reynolds &
! Acevedo (1993) Crop Sci 33:1273-1279 (18-24 m2 kg-1). The specific-leaf-area
! sources govern, because they measure the quantity the constant is. The
! released value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_sink_lai
IF (ios /= 0) canopy_sink_lai = 0
IF (canopy_sink_lai /= 1) canopy_sink_lai = 0
! ---------------------------------------------------------------------------- !
! Driver records 70 to 74, the canopy boundary-layer resistance and the hourly-
! wind pair. Both flags are end-of-file defaulted, so a short deck is
! unaffected and the daily path is untouched. Sources, derivation and the
! vanishing-attenuation assertion are in SHARED.f90.
! ---------------------------------------------------------------------------- !
READ (1,*,IOSTAT=ios) rca_wind_attenuation
IF (ios /= 0) rca_wind_attenuation = 0
IF (rca_wind_attenuation /= 1) rca_wind_attenuation = 0
READ (1,*,IOSTAT=ios) wind_hourly_mode
IF (ios /= 0) wind_hourly_mode = 0
IF (wind_hourly_mode /= 1) wind_hourly_mode = 0
READ (1,*,IOSTAT=ios) rwa_lai_blend
IF (ios /= 0) rwa_lai_blend = 1.0_wp
IF (rwa_lai_blend <= 0.0_wp) rwa_lai_blend = 1.0_wp
READ (1,'(A)',IOSTAT=ios) wind_hourly_file
IF (ios /= 0) wind_hourly_file = ''
CALL strip_driver_comment(wind_hourly_file)
READ (1,*,IOSTAT=ios) wind_meas_hgt_m
IF (ios /= 0) wind_meas_hgt_m = -1.0_wp
!------------------------------------------------------------------------------!
! Driver record 75 rss_theta_s_consistency: the soil-surface-resistance formula
! reads the layer-1 saturated water content that Richards and the van
! Genuchten-Mualem retention curve already use, instead of two mutually
! inconsistent driver constants. No new constant is introduced. The released
! value is 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) rss_theta_s_consistency
IF (ios /= 0) rss_theta_s_consistency = 0
IF (rss_theta_s_consistency /= 1) rss_theta_s_consistency = 0
!------------------------------------------------------------------------------!
! Driver record 76, the canopy-nitrogen base (structural) term. The released
! value is 1: only the non-structural pool is distributed exponentially and the
! base term is added back uniformly, i.e. Eq. (2) of Bertheloot, Martre &
! Andrieu (2008) Plant Physiol 148:1707-1720 with a base of 0.4 g N m-2 from
! that same page. Value 0 distributed the whole leaf-nitrogen pool
! exponentially with no base term; its code was removed from v1.0 and the value
! is refused below, so the record must be present and must read 1.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_n_base_profile
IF (ios /= 0) canopy_n_base_profile = 0
! Driver record 77 ci_solve_mode: mode 1 solves the intercellular CO2
! concentration from the flux identity by bisection, with no literal bounds.
! The released value is 1. See SHARED.f90 and assim.f90.
READ (1,*,IOSTAT=ios) ci_solve_mode
IF (ios /= 0) ci_solve_mode = 0
IF (ci_solve_mode /= 0) ci_solve_mode = 1
! Driver record 78 ft_pot_unstressed: mode 1 captures the unstressed potential
! transpiration (the hourly series and its daily sum) once per day on outer
! iteration 1 and uses it for BOTH the Richards demand and the water-stress
! factor denominator. The released value is 1.
READ (1,*,IOSTAT=ios) ft_pot_unstressed
IF (ios /= 0) ft_pot_unstressed = 0
IF (ft_pot_unstressed /= 0) ft_pot_unstressed = 1
! Driver record 79 ws_two_loci: mode 1 applies the soil-water factor to the
! stomatal SLOPE only, leaving the residual conductance and the photosynthetic
! capacity untouched. Leuning (1995) Plant Cell Environ 18:339-355, p.350,
! states the alternatives as either the slope or the capacity, and the slope is
! the one taken.
READ (1,*,IOSTAT=ios) ws_two_loci
IF (ios /= 0) ws_two_loci = 0
IF (ws_two_loci /= 0) ws_two_loci = 1
!------------------------------------------------------------------------------!
! Driver record 80 feddes_p2_hydrus: mode 1 compares the daily UNSTRESSED
! potential transpiration in cm d-1 against r2H = 0.5 and r2L = 0.1 cm d-1 and
! interpolates the Feddes p2 node linearly between p2h and p2l, which are the
! units and the rule HYDRUS-1D 4.17 Table 12.7 p.207 defines; the alternative
! node selection compares a demand in m d-1 against the same 0.5 threshold. The
! diagnostic output/daily_feddes_node_diag.csv is written in both modes.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) feddes_p2_hydrus
IF (ios /= 0) feddes_p2_hydrus = 0
IF (feddes_p2_hydrus /= 0) feddes_p2_hydrus = 1
!------------------------------------------------------------------------------!
! Driver record 81 spac_mode: the soil-plant-atmosphere branch of the stomatal
! and leaf-water solution. The syntax is "<mode> [psi_f [K_plant]]": modes 1, 2
! and 3 select the root-weighted soil water potential weighting (root fraction;
! root distribution times alpha; root distribution times the conductivity
! ratio); the optional second token is the reference leaf water potential in
! MPa, default -1.5; the optional third is the plant hydraulic conductance in
! mol m-2 s-1 MPa-1, default 0.0222. Fewer tokens keep the defaults. Refused
! unless records 77, 78, 23 and 24 are all 1. See SHARED.f90 and the leaf-water
! solve in assim.f90.
!------------------------------------------------------------------------------!
BLOCK
  CHARACTER(LEN=512) :: rec81_line
  INTEGER  :: ios2
  REAL(wp) :: psi_f_in, kpl_in
  READ (1,'(A)',IOSTAT=ios) rec81_line
  IF (ios /= 0) THEN
    spac_mode = 0
  ELSE
    CALL strip_driver_comment(rec81_line)
    READ (rec81_line,*,IOSTAT=ios2) spac_mode, psi_f_in, kpl_in
    IF (ios2 == 0) THEN
      spac_psi_f   = psi_f_in
      spac_K_plant = kpl_in
    ELSE
      READ (rec81_line,*,IOSTAT=ios2) spac_mode, psi_f_in
      IF (ios2 == 0) THEN
        spac_psi_f = psi_f_in
      ELSE
        READ (rec81_line,*,IOSTAT=ios2) spac_mode
        IF (ios2 /= 0) spac_mode = 0
      END IF
    END IF
  END IF
  IF (spac_mode < 0 .OR. spac_mode > 3) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 81: spac_mode must be 0, 1, 2 or 3; read ', spac_mode
    ERROR STOP 'rec 81: invalid spac_mode'
  END IF
  IF (spac_mode >= 1) THEN
    IF (ci_solve_mode /= 1) THEN
      WRITE (*,'(A)') ' ERROR rec 81: SPAC-Leuning (rec 81 >= 1) requires rec 77 = 1 (D1 flux-identity Ci solve);'
      WRITE (*,'(A)') '   the OFF Ci path carries the two literal Ci bounds the SPAC no-clamp claim excludes.'
      ERROR STOP 'rec 81 requires rec 77 = 1'
    END IF
    IF (ft_pot_unstressed /= 1) THEN
      WRITE (*,'(A)') ' ERROR rec 81: SPAC-Leuning requires rec 78 = 1 (Richards forced with the UNSTRESSED series);'
      WRITE (*,'(A)') '   with rec 78 = 0 the psi_L-limited E would drive the root sink, a coupling that v1.0 does not adopt.'
      ERROR STOP 'rec 81 requires rec 78 = 1'
    END IF
    IF (water_stress_mode /= 1 .OR. .NOT. use_soil_water_stress) THEN
      WRITE (*,'(A)') ' ERROR rec 81: SPAC-Leuning replaces f_t at the stomata and needs rec 23 = 1 and rec 24 = 1.'
      ERROR STOP 'rec 81 requires rec 23 = 1 and rec 24 = 1'
    END IF
    IF (spac_psi_f >= 0.0_wp) THEN
      WRITE (*,'(A,G14.6)') ' ERROR rec 81: psi_f must be negative (MPa); read ', spac_psi_f
      ERROR STOP 'rec 81: psi_f >= 0'
    END IF
    IF (spac_K_plant <= 0.0_wp) THEN
      WRITE (*,'(A,G14.6)') ' ERROR rec 81: K_plant must be > 0 (mol m-2 s-1 MPa-1); read ', spac_K_plant
      ERROR STOP 'rec 81: K_plant <= 0'
    END IF
    WRITE (*,'(A,I0,A,F7.3,A,G12.5,A)') ' rec 81 SPAC-Leuning ON: psi_s weighting ', spac_mode, &
      ', psi_f = ', spac_psi_f, ' MPa, K_plant = ', spac_K_plant, ' mol m-2 s-1 MPa-1'
  END IF
END BLOCK
! Driver record 82 tiller_carbon_gate: mode 1 sets cohort birth density by the
! constant-free SUPPLY RATIO (the window assimilate surplus over the organ-law
! cost of the cohort first leaves), keeps the Evers light law as TIMING only,
! and retires the axis cap and the die fraction to reported envelopes. Requires
! records 64, 65, 66 and 69. The released value is 1.
READ (1,*,IOSTAT=ios) tiller_carbon_gate
IF (ios /= 0) tiller_carbon_gate = 0
IF (tiller_carbon_gate /= 1) tiller_carbon_gate = 0
! Driver record 83 leaf_expansion_water: mode 1 scales the blade expansion
! increments by the wheat leaf-elongation response of Acker, Zhang, Wang, Yu &
! Zhou (2024) Ann Bot 134:787, doi:10.1093/aob/mcae114, Fig. 4B, linear in soil
! relative water content with no breakpoint, evaluated on the root-fraction-
! weighted water content of the existing soil state. The hydrology is untouched.
! Requires records 64 and 65. The released value is 1.
READ (1,*,IOSTAT=ios) leaf_expansion_water
IF (ios /= 0) leaf_expansion_water = 0
IF (leaf_expansion_water /= 1) leaf_expansion_water = 0
! Driver record 84 lai_from_organs: mode 1 makes the leaf area index the organ
! sum with the three integer demographic limiters of record 69 removed, and the
! specific leaf area becomes an OUTPUT column of output/daily_growth.csv.
! Requires records 64 to 69. The released value is 1.
READ (1,*,IOSTAT=ios) lai_from_organs
IF (ios /= 0) lai_from_organs = 0
IF (lai_from_organs /= 1) lai_from_organs = 0
IF (tiller_carbon_gate == 1) THEN
  IF (canopy_leaf_ranks /= 1 .OR. canopy_blade_area /= 1 .OR. canopy_tillering /= 1 .OR. canopy_sink_lai /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 82: tiller_carbon_gate requires rec 64, 65, 66 and 69 = 1' // &
      ' (S1-S3 chronology/geometry and the S6 lamina ledger).'
    ERROR STOP 'rec 82 requires rec 64, 65, 66, 69 = 1'
  END IF
END IF
IF (leaf_expansion_water == 1) THEN
  IF (canopy_leaf_ranks /= 1 .OR. canopy_blade_area /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 83: leaf_expansion_water requires rec 64 and 65 = 1 (rank chronology and blade geometry).'
    ERROR STOP 'rec 83 requires rec 64, 65 = 1'
  END IF
END IF
IF (lai_from_organs == 1) THEN
  IF (canopy_leaf_ranks /= 1 .OR. canopy_blade_area /= 1 .OR. canopy_tillering /= 1 .OR. &
      canopy_leaf_senesc /= 1 .OR. canopy_sheath_ear /= 1 .OR. canopy_sink_lai /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 84: lai_from_organs requires rec 64, 65, 66, 67, 68 and 69 = 1.'
    ERROR STOP 'rec 84 requires rec 64-69 = 1'
  END IF
END IF
IF (tiller_carbon_gate == 1 .OR. leaf_expansion_water == 1 .OR. lai_from_organs == 1) &
  WRITE (*,'(A,3I2)') ' D2 organ-level records ON (82 83 84): ', tiller_carbon_gate, leaf_expansion_water, lai_from_organs
! Driver record 85 rca_teh_full: mode 1 evaluates the canopy boundary-layer
! resistance as Teh (2006) Eq. 4.78 p.97 in full at every leaf area index, with
! the attenuation coefficient of record 12 and the leaf width of record 16, no
! leaf-area blending and no floor on the attenuation. The released value is 0.
! See SHARED.f90 and energy_balance_module.f90.
READ (1,*,IOSTAT=ios) rca_teh_full
IF (ios /= 0) rca_teh_full = 0
IF (rca_teh_full /= 1) rca_teh_full = 0
! Driver record 86 grainN_sink_mech: the grain-nitrogen sink of Bertheloot,
! Martre & Andrieu (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064,
! Eqs. (11) and (12) p.785, on the per-ear basis, with the measured initial
! grain nitrogen 2.4 mg culm-1 of that paper Table 1 p.784. Requires record 89
! at or above 1. The released value is 1; the full statement and the derived
! quantities are in SHARED.f90.
READ (1,*,IOSTAT=ios) grainN_sink_mech
IF (ios /= 0) grainN_sink_mech = 0
! Mode 5 is mode 1 with the Eq. (6) synthesis limb and the mobile-pool basis
! both off, so the difference between modes 1 and 5 isolates the Eq. (6) leg by
! a single subtraction.
IF (grainN_sink_mech < 0 .OR. grainN_sink_mech > 5) grainN_sink_mech = 0
! Driver record 87 leafN_cap_route: mode 1 applies the model's own leaf-nitrogen
! capacity law, leaf nitrogen at most sln_max_N times the leaf area index -
! already enforced at every other leaf-nitrogen inflow - to the post-anthesis
! surplus route as well, on the leaf area index gas exchange actually read that
! day, with the nitrogen the lamina cannot hold going to the existing mobile
! vegetative pool. The released value is 1.
READ (1,*,IOSTAT=ios) leafN_cap_route
IF (ios /= 0) leafN_cap_route = 0
IF (leafN_cap_route /= 1) leafN_cap_route = 0
! Driver record 88 sowing_model_day: this deck's own sowing model day, so that
! day 0 after sowing lands on model day n at all three anchored places (the
! fertilizer event, the residue event and the residue-evaporation mass). Model
! days are 1-based, so 0 cannot be a legitimate declaration and means the
! compiled default.
READ (1,*,IOSTAT=ios) sowing_model_day
IF (ios /= 0) sowing_model_day = 0
IF (sowing_model_day < 0) sowing_model_day = 0
! Driver record 89 grain_number_mech: mode 1 predicts grain number from the
! spike carbon each axis accrued, through Fischer (1985) J Agric Sci Camb
! 105:447-461 - the spike growth period from penultimate-leaf appearance p.449,
! the partitioning coefficient 0.44 of the above-ground increment p.452, one
! competent floret per 10 mg spike dry weight p.449, and 78 +/- 2 kernels per
! gram of spike dry weight at anthesis p.447 - resolved per axis, so that ears
! per square meter and grains per ear are both predicted. Modes 2 and 3 are the
! two halves taken separately. The released value is 1.
READ (1,*,IOSTAT=ios) grain_number_mech
IF (ios /= 0) grain_number_mech = 0
IF (grain_number_mech < 0 .OR. grain_number_mech > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 89: grain_number_mech must be 0, 1, 2 or 3; read ', grain_number_mech
  ERROR STOP 'rec 89: invalid grain_number_mech'
END IF
!------------------------------------------------------------------------------!
! Driver record 90 floret_survival_mech: the floret-site chain. 194 floret
! primordia per spike (19.4 spikelets, Fischer 1985 J Agric Sci Camb
! 105:447-461 p.450, times 10 per spikelet, Guo & Schnurbusch 2015 J Exp Bot
! 66:5945), a survival fraction that saturates with the spike assimilate each
! axis accrued per primordium (Guo, Slafer & Schnurbusch 2016 J Exp Bot
! 67:4241: the number of primordia that become fertile florets is adjusted to
! the actual assimilate availability), and a fertile-tiller test. The
! axis-level rule is Abichou (2016) Ch.1 sec.1.3 pp.20-21, where a proportion of
! tillers die while the others complete their cycle and produce an ear, the one
! that dies is the less developed, and the number surviving is strongly
! influenced by the availability of light, water and nutrients, with Fischer
! (1985) p.452 on higher-order tillers that bear no spikes being suppressed
! sooner in shaded crops. Culled youngest first (Abichou Ch.1 sec.2.5), with
! the spike carbon of the culled axes returned to the survivors.
! Mode 4, the released value, replaces the class-level bar by the continuous
! ear-bearing fraction MIN(1, spike dry matter per axis / (the record-98 spike
! rate times the window degree-Cd)) - the one affordability law record 99
! carries, applied to the spike sink in the Fischer window; the rest of the
! class stands spikeless (Fischer 1985 p.452) and its pool returns to the
! ear-bearers, and the main stem is never earless. No bar and no new constant.
! Modes 1 to 3 are the bar form, the survival-only arm and a rectangular-
! hyperbola sensitivity arm.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) floret_survival_mech
IF (ios /= 0) floret_survival_mech = 0
IF (floret_survival_mech < 0 .OR. floret_survival_mech > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 90: floret_survival_mech must be 0, 1, 2, 3 or 4; read ', floret_survival_mech
  ERROR STOP 'rec 90: invalid floret_survival_mech'
END IF
!------------------------------------------------------------------------------!
! Driver record 91 grain_fill_mech: the potential kernel growth rate of Yin,
! Goudriaan, Lantinga, Vos & Spiertz (2003) Ann Bot 91(3):361-371,
! doi:10.1093/aob/mcg029, Table 2, Eq. (9) with the Eq. (A1) rate shape,
! integrated over a thermal filling duration whose base temperature 6.619 degC
! and span 405.581 degree-Cd are derived from that same table by the
! reciprocal-of-duration method of Fischer (1985) J Agric Sci Camb 105:447-461
! Table 6 footnote p.455, with the day increment limited by the source. The
! potential kernel mass is not read on this path. Modes 2 and 3 are attribution
! arms (the rate shape removed; the capacity cap kept with the thermal duration
! stop added). The released value is 1. Without a rate law and a duration, 20
! of 22 decks ended with a kernel mass exactly equal to the potential, grain
! grew on 17 to 26 days only, and 58.5 to 451.8 g m-2 of carbon was stranded in
! the reserve.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) grain_fill_mech
IF (ios /= 0) grain_fill_mech = 0
IF (grain_fill_mech < 0 .OR. grain_fill_mech > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 91: grain_fill_mech must be 0, 1, 2 or 3; read ', grain_fill_mech
  ERROR STOP 'rec 91: invalid grain_fill_mech'
END IF
!------------------------------------------------------------------------------!
! Driver record 92 stem_wsc_mech: the stem water-soluble-carbohydrate
! (fructan) reserve. Mode 1 adds (a) a capacity, the reserve at most 0.50 of
! the stem dry mass, (b) an end to net deposition 15 days after anthesis, and
! (c) a DEMAND-DRIVEN draw equal to the day unmet grain demand, so that the
! fixed remobilization rate is not read on this path. Sources read at first
! hand: Schnyder (1993) New Phytol 123:233-245, pp.233, 238 and 239, and
! Gebbing & Schnyder (1999) Plant Physiol 121:871-878, pp.875 and 877, whose
! measured mobilization efficiency 0.97 +/- 0.07 is why the conversion
! efficiency 0.95 is kept unchanged. Modes 2, 3 and 4 are attribution and
! sensitivity arms (capacity and close only; demand-driven draw only;
! mode 1 with deposition restricted to a later development stage). The released
! value is 1.
! The pool is NOT filled only from post-anthesis overflow: record 55 is 1 on
! every deck and growth_module.f90 already banks 0.20 of every pre-anthesis
! stem increment, so at anthesis the pool holds 17.7-19.3 % of stem dry matter.
! The measured shortfalls the mode repairs are that the pool had no capacity and 8
! of 22 decks exceeded the measured 50 % ceiling; that net deposition never
! ended, the pool peaking 3 to 51 days after anthesis against a measured 10 to
! 20 and 17 of 22 decks finishing at their maximum; that the peak-to-anthesis
! ratio was 1.13 to 8.34 against a measured 2 to 4; and that remobilization was
! capped by a rate of 0.05 d-1 that bound on every draw day while 629.6 g m-2
! of grain demand went unmet.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) stem_wsc_mech
IF (ios /= 0) stem_wsc_mech = 0
IF (stem_wsc_mech < 0 .OR. stem_wsc_mech > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 92: stem_wsc_mech must be 0, 1, 2, 3 or 4; read ', stem_wsc_mech
  ERROR STOP 'rec 92: invalid stem_wsc_mech'
END IF
!------------------------------------------------------------------------------!
! Driver record 93 soil_n_mech: humus mineralization. Mode 1 replaces the flat
! constant background mineralization rate - which has no temperature response,
! no moisture response and no pool to deplete, and which creates nitrogen from
! nothing - by the humus mineralization of Brisson et al. (1998) Agronomie
! 18(5-6):311-346, doi:10.1051/agro:19980501, Eq. (31) p.341 and Eqs. (32),
! (33), (42) and (44) p.342 with Table VIII p.330, read at first hand from
! 200-dpi page renders. This is a construction ported from another model from
! its published description. The humus pool becomes a real depleting state and
! the constant rate is not read. Mode 2 is the attribution arm with both
! response functions set to 1, so that the difference from mode 1 isolates the
! temperature and moisture response plus pool depletion. The three soil-analysis
! values (clay 7.0 %, calcium carbonate 0.0 %, initial humus nitrogen
! 252 g N m-2) come from Dier et al. (2019) Field Crops Res 234:107-118 p.108,
! doi:10.1016/j.fcr.2019.02.013, read at first hand. The released value is 0.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) soil_n_mech
IF (ios /= 0) soil_n_mech = 0
IF (soil_n_mech < 0 .OR. soil_n_mech > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 93: soil_n_mech must be 0, 1 or 2; read ', soil_n_mech
  ERROR STOP 'rec 93: invalid soil_n_mech'
END IF
!------------------------------------------------------------------------------!
! Driver record 94 canopy_n_flag_zone: where the canopy leaf-nitrogen
! exponential applies. Bertheloot, Martre & Andrieu (2008) Plant Physiol
! 148:1707-1720 measured NO gradient above the flag-leaf collar - p.1710, no
! specific leaf nitrogen gradient was observed in the upper part of the canopy,
! the height of the flag leaf collar delimiting this upper zone with constant
! specific leaf nitrogen - and p.1713 gives the reason: the nitrogen was
! constant over the length of the flag leaf laminae although the difference in
! photon flux between their top and bottom was 40 %, because in Poaceae the
! light signal is integrated at the whole lamina level. That paper fitted the
! exponential only for the canopy layers below the flag leaf collar (p.1710).
! Applying it above the collar drives the top-of-canopy nitrogen per unit area
! to 4.6-6.0 g N m-2 against the 3.85 maximum measured in Silva-Perez et al.
! (2020) J Exp Bot 71(7):2299-2311 Fig. 4B p.2305, which supplies this model
! own canopy capacity.
! THE LAW: flat over the flag zone, exponential below, with the top value fixed
! by CONSERVATION of leaf nitrogen. The zone depth is the model's OWN flag-rank
! green lamina area index, an emergent state of the organ chain, and neither a
! fitted depth nor a constant. The canopy extinction coefficient 0.41 and the
! Vcmax-nitrogen slope 71 are not refitted; the own value 0.42 x 0.92 = 0.386
! of Bertheloot et al. (2008) p.1710 corroborates 0.41.
!   1 = the flat-top profile, composing with the deck's own record 76.
!   2 = attribution: flat top with the record-76 base term forced to zero.
!   3 = the released value: mode 1 plus the layered fraction of absorbed PAR,
!       so that the Eq. (6) synthesis limb of record 86 reads a per-lamina PAR
!       after extinction, as Eq. (17) p.785 of Bertheloot, Martre & Andrieu
!       (2008) Funct Plant Biol 35:781-796 defines it, instead of canopy-top
!       PAR on the whole green mass.
!   4 = attribution: the layered fraction only, no flat top.
!   5 = diagnostic only: physics identical to mode 0, writing
!       output/daily_canopy_n_profile.csv so that the profile can be measured
!       on the baseline trajectory.
!------------------------------------------------------------------------------!
READ (1,*,IOSTAT=ios) canopy_n_flag_zone
IF (ios /= 0) canopy_n_flag_zone = 0
IF (canopy_n_flag_zone < 0 .OR. canopy_n_flag_zone > 5) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 94: canopy_n_flag_zone must be 0, 1, 2, 3, 4 or 5; read ', &
    canopy_n_flag_zone
  ERROR STOP 'rec 94: invalid canopy_n_flag_zone'
END IF
! Driver record 95 tleaf_class_diag: mode 1 adds the diagnostic sunlit and
! shaded leaf-temperature solve at the converged demand state and the file
! output/hourly_tleaf_classes.csv. Nothing computed under this record is read
! by carbon, water or the energy balance. The record exists because the
! WheatDryFACE infrared camera measures EXPOSED LEAVES and the Maricopa
! radiometers measure the CANOPY, while the model has one node. The equations
! are Teh (2006) pp.68-69. The released value is 0.
READ (1,*,IOSTAT=ios) tleaf_class_diag
IF (ios /= 0) tleaf_class_diag = 0
IF (tleaf_class_diag < 0 .OR. tleaf_class_diag > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 95: tleaf_class_diag must be 0 or 1; read ', tleaf_class_diag
  ERROR STOP 'rec 95: invalid tleaf_class_diag'
END IF
! Driver record 96 sw_radiosity_mode: mode 1 replaces the transitional
! single-scatter shortwave by the closed-form two-layer shortwave radiosity,
! solved to infinite reflection order by the same 1 - rho_canopy rho_soil
! closure the longwave block in the same routine already uses, with a
! machine-precision conservation assertion. No new constant. The released value
! is 0.
READ (1,*,IOSTAT=ios) sw_radiosity_mode
IF (ios /= 0) sw_radiosity_mode = 0
IF (sw_radiosity_mode < 0 .OR. sw_radiosity_mode > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 96: sw_radiosity_mode must be 0 or 1; read ', sw_radiosity_mode
  ERROR STOP 'rec 96: invalid sw_radiosity_mode'
END IF
! Driver record 97 itemF_Nclose_mode: two corrections to the grain-nitrogen
! sink of record 86.
!   1 = the released value: the same-day Eq. (6) lamina synthesis respects the
!       leaf-nitrogen holding capacity IN THE SAME STEP, inside the same
!       routine and on the same day, instead of the capacity being enforced
!       before the inflow that violates it.
!   2 = the post-anthesis grain-nitrogen UPTAKE demand becomes the potential
!       RATE of Eq. (12) p.785 net of the standing mobile pool, instead of a
!       cumulative STOCK deficit that kept drawing nitrogen after the fill had
!       ended, plus a terminal close of the mobile vegetative pool into stem
!       nitrogen at maturity.
!   3 = both. No new constant. The Eq. (11) and (12) statements were read at
!       first hand from 200-dpi page renders.
READ (1,*,IOSTAT=ios) itemF_Nclose_mode
IF (ios /= 0) itemF_Nclose_mode = 0
IF (itemF_Nclose_mode < 0 .OR. itemF_Nclose_mode > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 97: itemF_Nclose_mode must be 0, 1, 2 or 3; read ', itemF_Nclose_mode
  ERROR STOP 'rec 97: invalid itemF_Nclose_mode'
END IF
! Driver record 98 spike_sink_mech: mode 2 makes the spike a SINK WITH PRIORITY
! OVER THE STEM RESERVE, its demand being the potential spike growth per axis
! per degree-day, 1.0767e-3 g axis-1 (degree-Cd)-1 from Fischer (1985) J Agric
! Sci Camb 105:447-461 Table 4 p.451 with p.449 (32.3 florets per spike, 10 mg
! spike dry weight per competent floret, a 300 degree-Cd spike-growth period),
! so no new constant enters. Mode 3 co-limits florets by the spike's own FLORAL
! dry matter, the structural pool (rachis, awns, glumes) being 0.635 and the
! floral pool (lemma, palea) 0.364 of non-grain spike dry matter (Rivera-Amado
! et al. 2019 Field Crops Res 240:154-167 p.160 sec. 3.4 and Table 7), so that
! fruiting efficiency becomes an outcome of composition. Mode 1 is both. The
! released value is 0.
READ (1,*,IOSTAT=ios) spike_sink_mech
IF (ios /= 0) spike_sink_mech = 0
IF (spike_sink_mech < 0 .OR. spike_sink_mech > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 98: spike_sink_mech must be 0, 1, 2 or 3; read ', spike_sink_mech
  ERROR STOP 'rec 98: invalid spike_sink_mech'
END IF
! Driver record 99 tiller_survival_mech: per-axis mortality over the measured
! mortality phase of Abichou (2016) PhD thesis (HAL tel-02316154) Ch.1 Table 4 -
! start Eqn T4.2, end Eqn T4.3 at flag-leaf ligulation of the main stem,
! youngest axes culled first (Ch.1 sec. 2.5 p.29), scaled by that thesis's own
! normalized mortality clock (p.45). An axis that cannot pay, out of its OWN
! light-shared share of the crop assimilate, for its OWN potential lamina
! increment does not complete its cycle. No new constant, and no spike, chaff,
! floret or kernel quantity is read, so the ear count is not circular with the
! spike mass record 90 culls on. The thesis's own ear number is an INPUT variable
! of its tillering model (p.45), so its Eqn T4.4 cannot supply the number
! without calibrating on an observed ear count.
READ (1,*,IOSTAT=ios) tiller_survival_mech
IF (ios /= 0) tiller_survival_mech = 0
!   2 = tiller carbon affordability: the same window and death order, but a
!       per-cohort DEMAND in g CH2O (the cohort's own lamina increment times the
!       leaf mass per area times the glucose requirement) against a SUPPLY that
!       is the cohort exact-layer light share of the previous day growth
!       supply, resolved as the minimum of supply and demand (the pattern of
!       Yin & Schapendonk 2004).
!   4 = the released value: the same ledgers, order and endpoint as mode 2, but
!       the realized cumulative deaths are capped by the fitted
!       mortality-phase shape of Abichou (2016) Ch.1 Eqn T4.4,
!       F = 1 - exp(-theta tau / (delta - tau)) with theta = 3.32 and
!       delta = 1.77 (p.36, sign read consistently with Fig. 11), applied to
!       the running estimate of the eventual deaths; nothing dies on day 1 and
!       F = 1 at flag-leaf ligulation.
IF (tiller_survival_mech == 3) THEN
  WRITE (*,'(A)') ' ERROR rec 99: tiller_survival_mech = 3 is RESERVED (mode 2 + a named eCO2 residual term; not enacted).'
  ERROR STOP 'rec 99 = 3: reserved'
END IF
IF (tiller_survival_mech == 5) THEN
  WRITE (*,'(A)') ' ERROR rec 99: tiller_survival_mech = 5 is RESERVED (endpoint-only cull at tau_em, sensitivity arm;' // &
    ' not enacted).'
  ERROR STOP 'rec 99 = 5: reserved'
END IF
IF (tiller_survival_mech < 0 .OR. tiller_survival_mech > 5) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 99: tiller_survival_mech must be 0, 1, 2 or 4 (3 and 5 reserved); read ', tiller_survival_mech
  ERROR STOP 'rec 99: invalid tiller_survival_mech'
END IF
! Driver record 100 stem_demand_mech: mode 1 makes the TRUE STEM state a DEMAND
! instead of being the algebraic residual of the partition table - its own
! maintenance respiration (Teh 2006 Table 7.1 p.146) plus a thermal-time
! structural growth demand built from the internode and sheath profiles of
! Abichou (2016) PhD thesis (HAL tel-02316154) Ch.2 Fig. 3(b) and 3(c) p.67
! through Eq. (6) p.419 of Bakhshandeh et al. (2012) Cereal Res Commun
! 40(3):413-422 used INVERTED - and resolves that demand against supply. Carbon
! the stem does not demand is REFUSED: it enters no organ and is not
! renormalized onto leaf, root or storage, so crop supply stops equaling crop
! allocation. No new constant.
!   2 = diagnostic only: the incumbent allocation is executed unchanged and
!       only the identity columns are written.
!   3 = mode 1 with the Eq. (6) height basis corrected. Modes 1 and 2 feed the
!       Eq. (6) inverse the telescoped COLLAR height, while the dependent
!       variable of Eq. (6) is the plant height measured from the soil surface
!       to the uppermost part of the plant with a ruler (Bakhshandeh et al.
!       2012, folio 415) - a basis mismatch that caps the demand at a 56.5 cm
!       collar on every deck and ends it before anthesis. Mode 3 simulates the
!       omitted organ instead of adding an observed number to the height: the
!       potential peduncle length is the model's OWN top internode length scaled
!       by the measured dimensionless ratio 36.7/20.5 = 1.79024 (Rivera-Amado
!       et al. 2019 Field Crops Res 240:154-167, doi:10.1016/j.fcr.2019.04.016,
!       Table 6 folio 161, the rank pair that table itself defines), ramped
!       over the deck's own window from flag-leaf ligule to ANTHESIS where the
!       stated validity of Eq. (6) ends; the height the model carries advances
!       only by the realized mass fraction of that potential, so height is
!       carbon-limited and emergent rather than prescribed. No observed length
!       is added to any height anywhere, and the inversion itself is unchanged.
!       Record 101 is not available with mode 3.
! The released value is 0.
READ (1,*,IOSTAT=ios) stem_demand_mech
IF (ios /= 0) stem_demand_mech = 0
IF (stem_demand_mech < 0 .OR. stem_demand_mech > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 100: stem_demand_mech must be 0, 1, 2, 3 or 4; read ', stem_demand_mech
  ERROR STOP 'rec 100: invalid stem_demand_mech'
END IF
! Driver record 101 surplus_fate_mech: mode 1 offers the carbon the true stem
! refuses under record 100 to the stem water-soluble-carbohydrate (fructan)
! pool up to its MEASURED capacity. Schnyder (1993) New Phytol 123:233-245,
! doi:10.1111/j.1469-8137.1993.tb03731.x, p.238, read at first hand from
! 200-dpi page images: stem storage accommodates excess photosynthate when the
! demand by the ears and other sinks, e.g. respiration and extension growth of
! the peduncle, is less than the current net rate of photosynthesis, with a
! capacity of up to 50 % of the dry weight of internodes. The conversion uses
! the glucose requirement for carbohydrates, 1.242 g CH2O g-1 DM (Teh 2006
! Table 7.2 p.147). No assimilation-feedback branch is built, because the same
! page measures its absence in wheat (the absence of a negative feedback on
! leaf photosynthesis following ear removal), so carbon the capacity refuses is
! an explicitly reported residual and no coefficient is invented. Mode 2 is the
! attribution arm with deposition closed at anthesis. Requires record 100 = 1.
! The released value is 0.
READ (1,*,IOSTAT=ios) surplus_fate_mech
IF (ios /= 0) surplus_fate_mech = 0
IF (surplus_fate_mech < 0 .OR. surplus_fate_mech > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 101: surplus_fate_mech must be 0, 1 or 2; read ', surplus_fate_mech
  ERROR STOP 'rec 101: invalid surplus_fate_mech'
END IF
! Driver record 102 stem_sheath_sep_mode: a derived, read-only split of the
! lumped stem pool into the structural true stem, a sheath estimate and the
! already labeled reserve, using the printed leaf-sheath dry-matter
! partitioning index of Rivera-Amado et al. (2019) Field Crops Res 240:154-167,
! doi:10.1016/j.fcr.2019.04.016, folio 158. The printed RANGE 0.16-0.20 is
! carried at both ends as two modes (1 = 0.16, 2 = 0.20) because no mean is
! printed. The true-stem mass is read in exactly two places, the on-only
! diagnostic output/daily_stem_sep.csv and the SHARED variable exposed for the
! height numerator; no allocation, nitrogen, reserve, respiration, mineral or
! structure statement changes. The released value is 0.
READ (1,*,IOSTAT=ios) stem_sheath_sep_mode
IF (ios /= 0) stem_sheath_sep_mode = 0
IF (stem_sheath_sep_mode < 0 .OR. stem_sheath_sep_mode > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 102: stem_sheath_sep_mode must be 0, 1 or 2; read ', stem_sheath_sep_mode
  ERROR STOP 'rec 102: invalid stem_sheath_sep_mode'
END IF
! Driver records 103 sdm4_dur_bracket and 104 sdm4_rate_bracket: the bracket
! members of the per-axis, per-degree-day, bounded structural stem mass demand
! of record 100 mode 4. Both records are read and stored on every run and are
! read by nothing unless record 100 = 4; with record 100 = 4 a value of 0 on
! either record is an error stop, because no member of the sourced bracket is
! privileged and no hidden midpoint exists.
!   record 103: 1, 2 and 3 = the collar onset rule of Fournier (2003), field,
!               with internode durations 277.0, 350.8 and 476.0 degree-Cd
!               (Kirby 1988 folio 135; Kang 2008 folio 1103-1104; Kirby folio
!               135); 4, 5 and 6 = the Malvoisin (1984) chamber onset rule with
!               the same three durations.
!   record 104: 1 = 0.7 x 1.0, 2 = 0.7 x 1.7, 3 = 1.7 mg DM (degree-Cd)-1 per
!               internode, capped at 0.7 x 3.8, 0.7 x 3.8 and 3.8 mg DM
!               (degree-Cd)-1 per axis (Kirby 1988 folio 139-140; the
!               structural factor 0.7 is from folio 137).
! Both released values are 0.
READ (1,*,IOSTAT=ios) sdm4_dur_bracket
IF (ios /= 0) sdm4_dur_bracket = 0
IF (sdm4_dur_bracket < 0 .OR. sdm4_dur_bracket > 6) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 103: sdm4_dur_bracket must be 0 (not selected) or 1-6; read ', sdm4_dur_bracket
  ERROR STOP 'rec 103: invalid sdm4_dur_bracket'
END IF
READ (1,*,IOSTAT=ios) sdm4_rate_bracket
IF (ios /= 0) sdm4_rate_bracket = 0
IF (sdm4_rate_bracket < 0 .OR. sdm4_rate_bracket > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 104: sdm4_rate_bracket must be 0 (not selected) or 1-3; read ', sdm4_rate_bracket
  ERROR STOP 'rec 104: invalid sdm4_rate_bracket'
END IF
! Driver records 105 wsc_reserve_law and 106 wsc_law_member: the bounded,
! deficit-driven draw on the stem reserve.
!   record 105 = 1: the record-92 deficit draw gets one ceiling, the pool times
!               the day thermal increment over the remaining filling thermal
!               time, i.e. the pool spread no faster than uniformly over the
!               remaining filling time. The timing is sourced from Schnyder
!               (1993) pp.233, 239 and 240 and Gebbing & Schnyder (1999) p.875;
!               the uniform spread itself is a model construction. No new
!               constant: the clock constants are those of record 91 and the
!               conversion efficiency is that of record 60.
!   record 106: 1 = the record-91 filling clock, 2 = the record-55 window
!               clock, 3 = no ceiling. Read by nothing unless record 105 = 1.
! Requires records 92 = 1 and 91 at or above 1. The released values are 2 and 0.
READ (1,*,IOSTAT=ios) wsc_reserve_law
IF (ios /= 0) wsc_reserve_law = 0
IF (wsc_reserve_law < 0 .OR. wsc_reserve_law > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 105: wsc_reserve_law must be 0, 1 or 2; read ', wsc_reserve_law
  ERROR STOP 'rec 105: invalid wsc_reserve_law'
END IF
READ (1,*,IOSTAT=ios) wsc_law_member
IF (ios /= 0) wsc_law_member = 0
IF (wsc_law_member < 0 .OR. wsc_law_member > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 106: wsc_law_member must be 0 (not selected) or 1-3; read ', wsc_law_member
  ERROR STOP 'rec 106: invalid wsc_law_member'
END IF
! Driver record 107 wsc_rel_window_member: the release-window-end bracket of
! record 105 mode 2. Inert unless record 105 = 2; with record 105 = 2 a value of
! 0 is an error stop, since no member is privileged. 1 = the record-91 filling
! clock, 2 = the record-59 clock. The released value is 2.
READ (1,*,IOSTAT=ios) wsc_rel_window_member
IF (ios /= 0) wsc_rel_window_member = 0
IF (wsc_rel_window_member < 0 .OR. wsc_rel_window_member > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 107: wsc_rel_window_member must be 0 (not selected) or 1-2; read ', &
    wsc_rel_window_member
  ERROR STOP 'rec 107: invalid wsc_rel_window_member'
END IF
! Driver record 108 height_reserve_basis: mode 1 makes the argument of
! Eq. (6) of Bakhshandeh et al. (2012) the stem mass less the POST-anthesis
! reserve EXCURSION, i.e. the stem mass minus the difference between the
! current reserve and the reserve latched at anthesis, which neutralizes only
! that excursion. It is independent of records 105, 106 and 107: it is applied
! whatever the reserve law is. The released value is 0.
READ (1,*,IOSTAT=ios) height_reserve_basis
IF (ios /= 0) height_reserve_basis = 0
IF (height_reserve_basis < 0 .OR. height_reserve_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 108: height_reserve_basis must be 0 (OFF) or 1; read ', height_reserve_basis
  ERROR STOP 'rec 108: invalid height_reserve_basis'
END IF
! ==================================================================== !
! Driver records 109 to 113, CANOPY HEIGHT FROM INTERNODE ELONGATION.
!   record 109 height_internode_law  0 = off / 1 = the law
!   record 110 hgt_int_window_member 1-6 (onset arm x duration bracket)
!   record 111 hgt_int_length_member 1-5 (the length potential bracket)
!   record 112 hgt_int_lineic_member 1-4 (the lineic-mass bracket)
!   record 113 hgt_int_ear_member    1-3 (the ear and peduncle-end bracket)
! Interlocks: record 109 = 1 requires records 108 = 0 and 100 = 0, since two
! height laws or two stem laws may not run together, and requires records
! 64 = 1 and 66 = 1, which build the leaf-rank chronology and the axis
! population the law reads; otherwise an error stop. A member value of 0 with
! record 109 = 1 is an error stop: no member is privileged and none may be
! defaulted. The released value of record 109 is 0 and of the four members 0.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) height_internode_law
IF (ios /= 0) height_internode_law = 0
IF (height_internode_law < 0 .OR. height_internode_law > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 109: height_internode_law must be 0 (OFF) or 1; read ', height_internode_law
  ERROR STOP 'rec 109: invalid height_internode_law'
END IF
READ (1,*,IOSTAT=ios) hgt_int_window_member
IF (ios /= 0) hgt_int_window_member = 0
IF (hgt_int_window_member < 0 .OR. hgt_int_window_member > 6) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 110: hgt_int_window_member must be 0 (not selected) or 1-6; read ', &
    hgt_int_window_member
  ERROR STOP 'rec 110: invalid hgt_int_window_member'
END IF
READ (1,*,IOSTAT=ios) hgt_int_length_member
IF (ios /= 0) hgt_int_length_member = 0
IF (hgt_int_length_member < 0 .OR. hgt_int_length_member > 5) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 111: hgt_int_length_member must be 0 (not selected) or 1-5; read ', &
    hgt_int_length_member
  ERROR STOP 'rec 111: invalid hgt_int_length_member'
END IF
READ (1,*,IOSTAT=ios) hgt_int_lineic_member
IF (ios /= 0) hgt_int_lineic_member = 0
IF (hgt_int_lineic_member < 0 .OR. hgt_int_lineic_member > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 112: hgt_int_lineic_member must be 0 (not selected) or 1-4; read ', &
    hgt_int_lineic_member
  ERROR STOP 'rec 112: invalid hgt_int_lineic_member'
END IF
READ (1,*,IOSTAT=ios) hgt_int_ear_member
IF (ios /= 0) hgt_int_ear_member = 0
IF (hgt_int_ear_member < 0 .OR. hgt_int_ear_member > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 113: hgt_int_ear_member must be 0 (not selected) or 1-3; read ', &
    hgt_int_ear_member
  ERROR STOP 'rec 113: invalid hgt_int_ear_member'
END IF
! ==================================================================== !
! Driver record 114 mineral_share_basis: 0 = the mineral allocation share is
! the positive part of the day organ dry-mass increment; 1 = the
! phloem-delivery basis, the share of element i being the delivery to organ i
! over the total delivery. input/driver_main.txt is positional, so nothing is
! renumbered. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_share_basis
IF (ios /= 0) mineral_share_basis = 0
IF (mineral_share_basis < 0 .OR. mineral_share_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 114: mineral_share_basis must be 0 (OFF) or 1; read ', mineral_share_basis
  ERROR STOP 'rec 114: invalid mineral_share_basis'
END IF
! ==================================================================== !
! Driver record 115 gnf_ng0_basis: the basis of the initial grain-nitrogen
! sink. 0 = per EAR, the incumbent; 1 = per KERNEL, the per-ear value scaled by
! the kernel number of the cohort. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) gnf_ng0_basis
IF (ios /= 0) gnf_ng0_basis = 0
IF (gnf_ng0_basis < 0 .OR. gnf_ng0_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 115: gnf_ng0_basis must be 0 (PER EAR, OFF) or 1 (PER KERNEL); read ', &
    gnf_ng0_basis
  ERROR STOP 'rec 115: invalid gnf_ng0_basis'
END IF
! ==================================================================== !
! Driver record 116 limb_census_mode: 0 = no file is opened; 1 = write
! output/limb_census.csv, which changes no arithmetic and no flux. The record
! adds no mechanism and has no prerequisite flag, so it is valid in every
! configuration. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) limb_census_mode
IF (ios /= 0) limb_census_mode = 0
IF (limb_census_mode < 0 .OR. limb_census_mode > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 116: limb_census_mode must be 0 (OFF) or 1 (ON); read ', &
    limb_census_mode
  ERROR STOP 'rec 116: invalid limb_census_mode'
END IF
! ==================================================================== !
! Driver record 117 chaff_growth_mode: the post-anthesis spike-non-grain
! (chaff) growth sink, as a share of the post-anthesis above-ground increment.
!   1 = 0.15  Fischer (1985) J Agric Sci Camb 105:447-461 folio 451
!   2 = 0.32  Slafer et al. (2015) Food Energy Secur 4:92-109 folio 97
!   3 = 0.43  Slafer et al. (2015) folio 97
!   4 = 0.20  Slafer et al. (2015) folio 97, low end
!   5 = 0.50  Slafer et al. (2015) folio 97, high end
! The member is a SELECTOR over published values, not a free fraction, so no
! member can be chosen by fit. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) chaff_growth_mode
IF (ios /= 0) chaff_growth_mode = 0
IF (chaff_growth_mode < 0 .OR. chaff_growth_mode > 5) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 117: chaff_growth_mode must be 0 (OFF) or 1-5 (bracket member); read ', &
    chaff_growth_mode
  ERROR STOP 'rec 117: invalid chaff_growth_mode'
END IF
! ==================================================================== !
! Driver record 118 kernel_mass_ceiling: 1 bounds the per-kernel mass by the
! MEASURED cultivar potential thousand-kernel weight of
! input/cultivar_yield.csv, inside the rate x duration law of record 91.
! Requires record 91 = 1 or 2; under record 91 = 0 or 3 it is an error stop,
! because the stock cap there already reads the potential and a second cap
! would be a silent no-op. Read before the driver file is closed. The released
! value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) kernel_mass_ceiling
IF (ios /= 0) kernel_mass_ceiling = 0
IF (kernel_mass_ceiling < 0 .OR. kernel_mass_ceiling > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 118: kernel_mass_ceiling must be 0 (OFF) or 1 (ON); read ', kernel_mass_ceiling
  ERROR STOP 'rec 118: invalid kernel_mass_ceiling'
END IF
IF (kernel_mass_ceiling == 1 .AND. (grain_fill_mech == 0 .OR. grain_fill_mech == 3)) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 118: kernel_mass_ceiling = 1 requires the rec 91 rate law (grain_fill_mech 1 or 2);' // &
    ' under rec 91 = 0 or 3 the legacy stock cap already reads TKW_potential. Read rec 91 = ', grain_fill_mech
  ERROR STOP 'rec 118 requires rec 91 = 1 or 2'
END IF
! ==================================================================== !
! Driver records 119 and 120 are RESERVED rows. They are read so that the
! positional numbering reaches record 121 and are required to be 0; a non-zero
! value is refused rather than defaulted. End of file gives 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) rec119_reserved
IF (ios /= 0) rec119_reserved = 0
READ (1,*,IOSTAT=ios) rec120_reserved
IF (ios /= 0) rec120_reserved = 0
IF (rec119_reserved /= 0 .OR. rec120_reserved /= 0) THEN
  WRITE (*,'(A,I0,A,I0)') ' ERROR rec 119/120: RESERVED rows must be 0 (rejected diagnostics, no active code); read ', &
    rec119_reserved, ' and ', rec120_reserved
  ERROR STOP 'rec 119/120: reserved rows must be 0'
END IF
! ==================================================================== !
! Driver record 121 n_payment_gate: ORGAN NITROGEN PAYMENT FOR STRUCTURAL
! GROWTH. 1 = pay at the structural floors (diagnostic); 2 = pay at the
! construction concentrations, the released value. At or above 1 the
! nitrate-reduction block runs BEFORE the carbon allocation, so that today
! carbon is paid with today's reduced nitrogen. Read before the driver file is
! closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_payment_gate
IF (ios /= 0) n_payment_gate = 0
IF (n_payment_gate < 0 .OR. n_payment_gate > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 121: n_payment_gate must be 0 (OFF), 1 (floors) or 2 (construction); read ', n_payment_gate
  ERROR STOP 'rec 121: invalid n_payment_gate'
END IF
! ==================================================================== !
! Driver record 122 lai_n_arrest: LAMINA AREA CANNOT STAND AHEAD OF ITS
! NITROGEN. 1 = arrest the day lamina-area increment at the paid fraction of
! record 121; 2 = arrest and shed the green area below the structural specific-
! leaf-nitrogen floor, oldest rank first. Requires record 121 at or above 1 and
! record 84 = 1; refused otherwise. The released value is 2.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lai_n_arrest
IF (ios /= 0) lai_n_arrest = 0
IF (lai_n_arrest < 0 .OR. lai_n_arrest > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 122: lai_n_arrest must be 0 (OFF), 1 (arrest) or 2 (arrest + shed); read ', lai_n_arrest
  ERROR STOP 'rec 122: invalid lai_n_arrest'
END IF
IF (lai_n_arrest >= 1 .AND. n_payment_gate == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 122: lai_n_arrest >= 1 requires rec 121 n_payment_gate >= 1 (the paid fraction p_L)'
  ERROR STOP 'rec 122: requires rec 121'
END IF
IF (lai_n_arrest >= 1 .AND. lai_from_organs /= 1) THEN
  WRITE (*,'(A)') ' ERROR rec 122: lai_n_arrest >= 1 requires rec 84 lai_from_organs = 1 (LAI from the organ chain)'
  ERROR STOP 'rec 122: requires rec 84'
END IF
! ==================================================================== !
! Driver record 123 tiller_n_afford: the NITROGEN TERM OF THE TILLER
! AFFORDABILITY TEST. 1 = the per-cohort ledger of record 99 mode 4 also
! carries the cohort lamina construction nitrogen (demand at the record-121
! construction concentration, supply at the last gate paid fraction), and the
! affordability is the minimum of the carbon and nitrogen fractions. Requires
! record 99 = 4. With record 121 = 0 it is identical to record 99 mode 4. The
! released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) tiller_n_afford
IF (ios /= 0) tiller_n_afford = 0
IF (tiller_n_afford < 0 .OR. tiller_n_afford > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 123: tiller_n_afford must be 0 (OFF) or 1; read ', tiller_n_afford
  ERROR STOP 'rec 123: invalid tiller_n_afford'
END IF
IF (tiller_n_afford >= 1 .AND. tiller_survival_mech /= 4) THEN
  WRITE (*,'(A)') ' ERROR rec 123: tiller_n_afford = 1 requires rec 99 tiller_survival_mech = 4' // &
    ' (the tiller carbon-affordability ledger it extends)'
  ERROR STOP 'rec 123: requires rec 99 = 4'
END IF
! ==================================================================== !
! Driver record 124 root_front_law: 0 = the sinusoidal descent of Borg &
! Grimes (1986); 1, 2 and 3 = linear descent from emergence at 7.9, 9.6 and
! 11.3 mm d-1 (the abstract range of Kirkegaard & Lilley 2007: lower bound,
! midpoint, upper bound) from the 0.10 m seedling floor to the deck's own maximum
! root depth. No other quantity changes. The released value is 2.
! ==================================================================== !
READ (1,*,IOSTAT=ios) root_front_law
IF (ios /= 0) root_front_law = 0
IF (root_front_law < 0 .OR. root_front_law > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 124: root_front_law must be 0 (Borg & Grimes), 1, 2 or 3 (7.9/9.6/11.3 mm d-1); read ', &
    root_front_law
  ERROR STOP 'rec 124: invalid root_front_law'
END IF
! ==================================================================== !
! Driver record 125 kernel_potential_law: the per-kernel grain-nitrogen ceiling
! and the kernel-mass potential. 1 = the ceiling, 1.38 mg N kernel-1 (Martre,
! Porter, Jamieson & Triboi 2003 Plant Physiol 133:1959-1967 Table I p.1961);
! 2 = ceiling and potential; 3 = potential only. Modes 2 and 3 require
! record 118 = 0, record 91 = 1 or 2, and the cultivar rows 7-8. Read before
! the driver file is closed. The released value is 1; modes 2 and 3 are an
! inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) kernel_potential_law
IF (ios /= 0) kernel_potential_law = 0
IF (kernel_potential_law < 0 .OR. kernel_potential_law > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 125: kernel_potential_law must be 0, 1, 2 or 3; read ', kernel_potential_law
  ERROR STOP 'rec 125: invalid kernel_potential_law'
END IF
IF (kernel_potential_law >= 2 .AND. kernel_mass_ceiling == 1) THEN
  WRITE (*,'(A)') ' ERROR rec 125: modes 2/3 (kernel potential) require rec 118 = 0 -- two ceilings on one kernel' // &
    ' would make one a silent no-op'
  ERROR STOP 'rec 125 >= 2 requires rec 118 = 0'
END IF
IF (kernel_potential_law >= 2 .AND. (grain_fill_mech == 0 .OR. grain_fill_mech == 3)) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 125: modes 2/3 require the rec 91 rate law (1 or 2); read rec 91 = ', grain_fill_mech
  ERROR STOP 'rec 125 >= 2 requires rec 91 = 1 or 2'
END IF
! ==================================================================== !
! Driver record 126 mineral_grain_sink_cap: 1 bounds each element daily grain
! loading by its maximum grain concentration (optional column 16 of
! input/driver_minerals.csv) times the day grain dry-matter increment; refused
! amounts stay in the leaf (uptake share) or in the dead tissue (export). Read
! before the driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_grain_sink_cap
IF (ios /= 0) mineral_grain_sink_cap = 0
IF (mineral_grain_sink_cap < 0 .OR. mineral_grain_sink_cap > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 126: mineral_grain_sink_cap must be 0 (OFF) or 1; read ', mineral_grain_sink_cap
  ERROR STOP 'rec 126: invalid mineral_grain_sink_cap'
END IF
! ==================================================================== !
! Driver record 127 repro_surplus_fate: 1 books the part of the
! grain-sink-full reproductive surplus that record 92 refuses to label as stem
! reserve to the explicit residual, instead of letting it become unpaid stem
! structure. Read before the driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) repro_surplus_fate
IF (ios /= 0) repro_surplus_fate = 0
IF (repro_surplus_fate < 0 .OR. repro_surplus_fate > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 127: repro_surplus_fate must be 0 (OFF) or 1; read ', repro_surplus_fate
  ERROR STOP 'rec 127: invalid repro_surplus_fate'
END IF
! ==================================================================== !
! Driver record 128 n_construction_status: TWO-POOL STATUS-DEPENDENT
! CONSTRUCTION NITROGEN. 1 replaces the single-tier payment of record 121 by a
! two-pass floor and target payment (pass 1 at the floor with luxury allowed,
! pass 2 up to a status-capped target with no luxury). Modes 2 and 3 are
! diagnostic variants of mode 1 that change only the lamina floor read in
! pass 1: 2 = the model's own minimum functional specific leaf nitrogen of record
! 122; 3 = the minimum specific leaf nitrogen 0.8 g N m-2 attributed to
! Jamieson & Semenov (2000) and Martre et al. (2006), which could not be
! verified in either paper. Requires record 121 at or above 1 in every mode.
! The released value is 2.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_construction_status
IF (ios /= 0) n_construction_status = 0
IF (n_construction_status < 0 .OR. n_construction_status > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 128: n_construction_status must be 0, 1, 2 or 3 (0 OFF, 1 CANDIDATE, 2/3 DIAGNOSTIC); read ', &
    n_construction_status
  ERROR STOP 'rec 128: invalid n_construction_status'
END IF
IF (n_construction_status >= 1 .AND. n_payment_gate == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 128: n_construction_status >= 1 requires rec 121 n_payment_gate >= 1'
  ERROR STOP 'rec 128 requires rec 121 >= 1'
END IF
! ==================================================================== !
! Driver record 129 sink_full_carbon_fate: 1 lifts the record-92 deposition
! WINDOW for the day surplus while the grain sink is full, so that storage
! continues to the existing capacity (Schnyder 1993 New Phytol 123:233-245
! pp.239-240); whatever the capacity still refuses is booked to an explicit,
! reported unallocated-assimilate carbon loss, which is neither respiration nor
! stem structure, superseding the record-127 residual for this surplus.
! Requires record 127 = 1. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) sink_full_carbon_fate
IF (ios /= 0) sink_full_carbon_fate = 0
IF (sink_full_carbon_fate < 0 .OR. sink_full_carbon_fate > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 129: sink_full_carbon_fate must be 0 (OFF) or 1; read ', sink_full_carbon_fate
  ERROR STOP 'rec 129: invalid sink_full_carbon_fate'
END IF
IF (sink_full_carbon_fate == 1 .AND. repro_surplus_fate /= 1) THEN
  WRITE (*,'(A)') ' ERROR rec 129: sink_full_carbon_fate = 1 requires rec 127 repro_surplus_fate = 1' // &
    ' (there is no explicit residual mechanism for rec 129 to supersede).'
  ERROR STOP 'rec 129 requires rec 127 = 1'
END IF
! ==================================================================== !
! Driver record 130 grain_n_min_conc: the GRAIN NITROGEN MINIMUM
! CONCENTRATION. 1 = cumulative mode, the room being the grain nitrogen over
! the minimum concentration less the grain mass; 2 = guarded daily-increment
! mode, the day grain increment bounded by yesterday's grain-nitrogen import over
! the minimum concentration, not evaluated below 10 % of the potential kernel
! mass; 3 = the soft structural floor of Martre, Porter, Jamieson & Triboi
! (2003) Plant Physiol 133:1959-1967 Eq. (3) p.1964, which replaces the fixed
! 1.5 %DM ceiling by a time-varying per-grain structural and metabolic nitrogen
! demand. Requires record 127 = 1, because the refused dry matter of the floor
! needs the defined fate the record-127 and record-129 chain supplies. Read
! before the driver file is closed. The released value is 3; modes 1 and 2 are
! an inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) grain_n_min_conc
IF (ios /= 0) grain_n_min_conc = 0
IF (grain_n_min_conc < 0 .OR. grain_n_min_conc > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 130: grain_n_min_conc must be 0 (OFF), 1 (cumulative), 2 (guarded daily increment)' // &
    ' or 3 (Martre 2003 soft structural floor); read ', grain_n_min_conc
  ERROR STOP 'rec 130: invalid grain_n_min_conc'
END IF
IF (grain_n_min_conc >= 1 .AND. repro_surplus_fate /= 1) THEN
  WRITE (*,'(A)') ' ERROR rec 130: grain_n_min_conc >= 1 requires rec 127 repro_surplus_fate = 1' // &
    ' (the floor''s own refused dry matter needs a defined fate; the refused DM must never be silently lost).'
  ERROR STOP 'rec 130 requires rec 127 = 1'
END IF
! Record 130 mode 3 divides by the kernel-mass potential of records 118 and 125,
! computed in the same routine; under record 118 = 0 with record 125 below 2
! that quantity is zero and the floor refuses ALL grain filling. The dependency
! is stated here rather than met by accident. Reader-only.
IF (grain_n_min_conc == 3 .AND. kernel_mass_ceiling /= 1 .AND. kernel_potential_law < 2) THEN
  WRITE (*,'(A)') ' ERROR rec 130: grain_n_min_conc = 3 requires a kernel-mass potential, rec 118 kernel_mass_ceiling = 1' // &
    ' or rec 125 kernel_potential_law >= 2 (the Martre floor is scaled by grain_cap, which is 0 without one).'
  ERROR STOP 'rec 130 = 3 requires rec 118 = 1 or rec 125 >= 2'
END IF
! ==================================================================== !
! Driver record 131 mineral_supply_law: 1 = the solution concentration and
! generalized active uptake of iron, zinc, copper and manganese, and the
! calcium xylem concentration, driven by the calibrated parameters in columns
! 17 and 18 of input/driver_minerals.csv. Read before the driver file is
! closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_supply_law
IF (ios /= 0) mineral_supply_law = 0
IF (mineral_supply_law < 0 .OR. mineral_supply_law > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 131: mineral_supply_law must be 0 (OFF) or 1; read ', mineral_supply_law
  ERROR STOP 'rec 131: invalid mineral_supply_law'
END IF
! ==================================================================== !
! Driver record 132 leafn_hold_law: leaf-nitrogen holding capacity and
! post-anthesis surplus routing. 1 = retire the post-anthesis Eq. (6) lamina
! re-synthesis, so that the stem and chaff nitrogen drawn by the ladder stays
! in the mobile pool, available to the grain through the existing import law,
! instead of being credited back into the leaf pool. 2 = mode 1 plus the
! re-sourced holding-capacity ceiling 2.5 g N m-2 (Pask 2009 PhD thesis,
! University of Nottingham, p.29, printed) in place of 3.8 in the general
! post-anthesis capacity ladder. Read before the driver file is closed. The
! released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) leafn_hold_law
IF (ios /= 0) leafn_hold_law = 0
IF (leafn_hold_law < 0 .OR. leafn_hold_law > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 132: leafn_hold_law must be 0 (OFF), 1 or 2; read ', leafn_hold_law
  ERROR STOP 'rec 132: invalid leafn_hold_law'
END IF
! ==================================================================== !
! Driver record 133 leafn_turnover_law: the intended production mechanism for
! the same process as record 132, which stays in the tree as the diagnostic
! that measured late lamina nitrogen as too high and is required to be 0
! whenever this record is at or above 1 - the two are alternatives and are
! never stacked. 1 = the state-driven two-stage order-of-claims rule: stage A
! resolves the existing grain sink first against current assimilation plus the
! mobile pool, drawing stem and chaff for the grain only; stage B lets the
! lamina re-synthesis draw ONLY the leftover of assimilation and the mobile
! pool, with no stem or chaff of its own, the unmet remainder being refused and
! reported. 2 = mode 1 plus the same re-sourced 2.5 g N m-2 ceiling record 132
! already defines. Read before the driver file is closed. The released value
! is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) leafn_turnover_law
IF (ios /= 0) leafn_turnover_law = 0
IF (leafn_turnover_law < 0 .OR. leafn_turnover_law > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 133: leafn_turnover_law must be 0 (OFF), 1 or 2; read ', leafn_turnover_law
  ERROR STOP 'rec 133: invalid leafn_turnover_law'
END IF
IF (leafn_turnover_law >= 1 .AND. leafn_hold_law >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A)') ' ERROR rec 133 requires rec 132 (leafn_hold_law) = 0: the two are ALTERNATIVES,' // &
    ' never stacked (rec 132 is the diagnostic that proved late lamina N is too high; rec 133 is the real' // &
    ' mechanism). Read rec 132 = ', leafn_hold_law, ', rec 133 = ', leafn_turnover_law, '.'
  ERROR STOP 'rec 133 requires rec 132 = 0'
END IF
! ==================================================================== !
! Driver record 134 n_afford_priority_law: NITROGEN AFFORDABILITY AND
! MORTALITY-PRIORITY LAW. 1 = the combined-trigger eligibility and rate law
! gates tiers 1 and 2 (leaf senescence and tiller-axis mortality) and a harder
! gate adds a nitrogen term to tier 3 (floret and spike survival); the tiller
! mortality TARGET becomes recoverable. Requires record 123 at or above 1, so
! that this record hooks into the record-123 per-cohort ledger rather than
! building a parallel one. Read before the driver file is closed. The released
! value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_afford_priority_law
IF (ios /= 0) n_afford_priority_law = 0
IF (n_afford_priority_law < 0 .OR. n_afford_priority_law > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 134: n_afford_priority_law must be 0 (OFF) or 1; read ', n_afford_priority_law
  ERROR STOP 'rec 134: invalid n_afford_priority_law'
END IF
IF (n_afford_priority_law >= 1 .AND. tiller_n_afford == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 134: n_afford_priority_law = 1 requires rec 123 tiller_n_afford >= 1' // &
    ' (rec 134 hooks into rec 123''s own per-cohort ledger)'
  ERROR STOP 'rec 134 requires rec 123 >= 1'
END IF
! ==================================================================== !
! Driver record 135 n_payment_luxury_pass2: 1 lets pass 2 of record 128 (the
! lamina tier only) draw a SECOND, bounded luxury tranche from the SAME
! stem-luxury pool pass 1 already accesses, capped at the same release rate
! applied again to whatever of that pool pass 1 left unused, so no new
! parameter enters. Requires record 128 at or above 1. Read before the driver
! file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_payment_luxury_pass2
IF (ios /= 0) n_payment_luxury_pass2 = 0
IF (n_payment_luxury_pass2 < 0 .OR. n_payment_luxury_pass2 > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 135: n_payment_luxury_pass2 must be 0 (OFF) or 1; read ', n_payment_luxury_pass2
  ERROR STOP 'rec 135: invalid n_payment_luxury_pass2'
END IF
IF (n_payment_luxury_pass2 >= 1 .AND. n_construction_status == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 135: n_payment_luxury_pass2 = 1 requires rec 128 n_construction_status >= 1' // &
    ' (rec 135 draws into rec 128''s own pass-2 lamina demand)'
  ERROR STOP 'rec 135 requires rec 128 >= 1'
END IF
! ==================================================================== !
! Driver record 136 n_afford_standing_basis: 1 re-points tiers 1 and 2 of
! record 134 (leaf senescence and tiller-axis mortality, both eligibility and
! per-cohort severity) to a whole-canopy STANDING nitrogen ratio - the leaf
! nitrogen stock against the construction target times the green leaf area
! index, cumulative since the mortality window opens - instead of the
! record-134 own increment ratio, leaving tier 3 unchanged. 2 = mode 1 with
! tier 3 also re-pointed. Requires record 134 at or above 1, since this record
! only changes WHICH basis the record-134 trigger reads. Read before the driver
! file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_afford_standing_basis
IF (ios /= 0) n_afford_standing_basis = 0
IF (n_afford_standing_basis < 0 .OR. n_afford_standing_basis > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 136: n_afford_standing_basis must be 0 (OFF), 1 or 2; read ', n_afford_standing_basis
  ERROR STOP 'rec 136: invalid n_afford_standing_basis'
END IF
IF (n_afford_standing_basis >= 1 .AND. n_afford_priority_law == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 136: n_afford_standing_basis >= 1 requires rec 134 n_afford_priority_law >= 1' // &
    ' (rec 136 only re-points rec 134''s own combined trigger)'
  ERROR STOP 'rec 136 requires rec 134 >= 1'
END IF
! ==================================================================== !
! Driver record 137 floret_n_law: a NITROGEN TERM ON FLORET SURVIVAL PER
! SPIKELET. 1 = a window-cumulative standing nitrogen ratio, on the same
! ingredients and convention as record 136 but with its own accumulators and
! its own window, scales the constant-hazard floret-survival saturation of the
! Fischer chain. One nitrogen term on grain number and never two: when this
! record is on, the tier-3 ear-fraction nitrogen cap of record 134 is skipped.
! Refused if record 136 = 2, which would re-point the now-skipped cap. Read
! before the driver file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) floret_n_law
IF (ios /= 0) floret_n_law = 0
IF (floret_n_law < 0 .OR. floret_n_law > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 137: floret_n_law must be 0 (OFF) or 1; read ', floret_n_law
  ERROR STOP 'rec 137: invalid floret_n_law'
END IF
IF (floret_n_law >= 1 .AND. n_afford_standing_basis == 2) THEN
  WRITE (*,'(A)') ' ERROR rec 137: floret_n_law = 1 is incompatible with rec 136 n_afford_standing_basis' // &
    ' = 2 (that mode re-points record 134''s tier-3 ear-fraction N cap, which rec 137 skips -- one' // &
    ' nitrogen locus on grain number, never two)'
  ERROR STOP 'rec 137 incompatible with rec 136 = 2'
END IF
! ==================================================================== !
! Driver record 138 spike_n_priority: the SPIKE AND CHAFF AS ITS OWN
! STRUCTURAL-NITROGEN TIER, PAID BEFORE THE TRUE STEM. 1 = the spike and chaff
! become their own structural-nitrogen tier, paid before the true stem, at the
! SAME floor the combined tier already used, so this is a priority and order
! change with no new constant. 2 = mode 1 with the chaff tier floor replaced by
! the measured chaff, internode and peduncle structural-nitrogen minimum
! 3.5e-3 g N g-1 DM of Bertheloot, Martre & Andrieu (2008) Plant Physiol
! 148:1707-1720 p.1711; diagnostic. Requires record 121 at or above 1, since
! there is otherwise no tier structure to split. Read before the driver file is
! closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) spike_n_priority
IF (ios /= 0) spike_n_priority = 0
IF (spike_n_priority < 0 .OR. spike_n_priority > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 138: spike_n_priority must be 0 (OFF), 1 or 2; read ', spike_n_priority
  ERROR STOP 'rec 138: invalid spike_n_priority'
END IF
IF (spike_n_priority >= 1 .AND. n_payment_gate == 0) THEN
  WRITE (*,'(A)') ' ERROR rec 138: spike_n_priority >= 1 requires rec 121 n_payment_gate >= 1' // &
    ' (rec 138 splits rec 121''s own combined stem+chaff tier)'
  ERROR STOP 'rec 138 requires rec 121 >= 1'
END IF
! ==================================================================== !
! Driver record 139 n_uptake_diffusion: DIFFUSION TO THE ROOT SURFACE, the
! third term of the nitrogen-uptake supply law.
!   1 = the released value: per layer and ion, the Barber-Cushman form
!       S_diff = MAX(0, c_ion - C_min_ion) RLD De tk 86400 x 14.0, with the
!       effective diffusion coefficient De = Dl theta f_l(theta) / b and
!       f_l the same impedance pair mineral_module.f90 already uses, ADDED to
!       the per-layer per-ion supply cap alongside the mass-flow and
!       transporter terms: p_ion = MIN(pool, active + mass flow + diffusion).
!   2 = the kinetic (active) term is instead BOUNDED by transport delivery
!       (mass flow plus diffusion): p_ion = MIN(pool, MIN(active, mass flow +
!       diffusion)). Roose & Fowler (2004) J Theor Biol 228:173-184 Eq. (4.1)
!       is the root-surface Michaelis-Menten boundary condition, recovered
!       exactly in the transport-replete limit, and this is its
!       transport-limited asymptote, in which uptake is throttled toward the
!       delivery rate as diffusive resupply slows. The three terms themselves
!       are unchanged from mode 1; only their combination changes, and no new
!       constant enters.
!   3 = mode 2 with the diffusion term re-sourced. (a) The impedance pair
!       f_l = MAX(0, 1.60 theta - 0.17), which has no source and is zero at the
!       rooted-layer water content of the sandy decks (0.085-0.109), is
!       replaced by f_l = 1.1 MAX(0, theta - theta_th), Olesen, Moldrup,
!       Yamaguchi & Rolston (2001) Soil Sci 166:89-96 p.94 Eqs. (7a) and (8a),
!       the constant-slope impedance-factor model fitted across 22 soils from 0
!       to 53 % clay including true sands, with the threshold water content
!       from that paper's own p.93 Eq. (5) on the deck's own clay fraction, silt
!       fraction and bulk density (input/driver_N.csv rows 68-70). (b) The
!       lumped diffusion term gains the steady-state radial-diffusion
!       prefactor 2 pi / ln(r1 / (1.65 r0)) of Baldwin, Nye & Tinker (1973)
!       Plant Soil 38:621-635 pp.624-625 Eqs. [ii] to [iv], with
!       r1 = (pi RLD)^-0.5 and a non-wheat-specific root radius 0.05 cm, the
!       worked-example value shared by Baldwin (1973) p.626 and Cushman (1979)
!       Part II Fig. 1. The mode-1 and mode-2 subtraction of the
!       half-saturation constant as a zero-flux threshold is REMOVED, since no
!       source supports that use; the diffusion term now uses the full
!       concentration as the driving gradient, and the smooth approach to zero
!       as concentration falls lives entirely in the unchanged
!       Michaelis-Menten form of the active term. No nitrate minimum
!       concentration was found in any source on this machine.
! Requires the layered soil-nitrogen path of input/driver_N.csv, since the
! per-layer loop this record hooks does not execute on the single-pool path;
! mode 3 additionally requires rows 68-70. Both checks run once driver_N.csv
! has been read, because that file is read after this one. Read before the
! driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_uptake_diffusion
IF (ios /= 0) n_uptake_diffusion = 0
IF (n_uptake_diffusion < 0 .OR. n_uptake_diffusion > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 139: n_uptake_diffusion must be 0 (OFF), 1 (mode 1), 2 (mode 2)' // &
    ' or 3 (mode 3, sourced impedance+geometry); read ', n_uptake_diffusion
  ERROR STOP 'rec 139: invalid n_uptake_diffusion'
END IF
! ==================================================================== !
! Driver record 140 mineral_supply_mech: the SAME three-component supply law as
! record 139, per element, replacing the calibrated generalized active term of
! record 131 mode 1. 1 = per element (iron, zinc, copper, magnesium,
! manganese; calcium excluded, its xylem route unchanged) mass flow plus
! diffusion to the root surface plus, for zinc only, a bounded high-affinity
! transporter ceiling. Requires record 131 at or above 1, checked here, since
! record 131 is read earlier in this file. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_supply_mech
IF (ios /= 0) mineral_supply_mech = 0
IF (mineral_supply_mech < 0 .OR. mineral_supply_mech > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 140: mineral_supply_mech must be 0 (OFF) or 1; read ', mineral_supply_mech
  ERROR STOP 'rec 140: invalid mineral_supply_mech'
END IF
IF (mineral_supply_mech >= 1 .AND. mineral_supply_law < 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 140: mineral_supply_mech >= 1 requires rec 131 mineral_supply_law >= 1' // &
    ' (the calibrated solution concentration C_soln,e this record''s mass-flow and diffusion terms both read' // &
    ' is record 131''s own state). Read mineral_supply_law = ', mineral_supply_law
  ERROR STOP 'rec 140 requires mineral_supply_law >= 1'
END IF
! ==================================================================== !
! Driver record 141 lamina_n_remob_pre: PRE-ANTHESIS LAMINA NITROGEN
! REMOBILIZATION TO THE STEM AND CHAFF STRUCTURAL TIER. 1 = the lamina's own
! nitrogen above its functional minimum specific leaf nitrogen is the LAST
! resort for the structural tiers, released at the model's own live Eq. (5) rate,
! with no new constant. Requires record 121 at or above 1 (there is otherwise
! no structural tier to pay into) and record 122 = 2 (the shed rule that
! retires nitrogen-emptied area must be live, or the lamina would keep
! nitrogen-empty standing area). The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_n_remob_pre
IF (ios /= 0) lamina_n_remob_pre = 0
IF (lamina_n_remob_pre < 0 .OR. lamina_n_remob_pre > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 141: lamina_n_remob_pre must be 0 (OFF) or 1; read ', lamina_n_remob_pre
  ERROR STOP 'rec 141: invalid lamina_n_remob_pre'
END IF
IF (lamina_n_remob_pre >= 1 .AND. (n_payment_gate == 0 .OR. lai_n_arrest /= 2)) THEN
  WRITE (*,'(A,I0,A,I0)') ' ERROR rec 141: lamina_n_remob_pre >= 1 requires rec 121 n_payment_gate >= 1' // &
    ' AND rec 122 lai_n_arrest = 2 (the shed rule that retires N-emptied area must be live). Read rec 121 = ', &
    n_payment_gate, ', rec 122 = ', lai_n_arrest
  ERROR STOP 'rec 141 requires rec 121 >= 1 and rec 122 = 2'
END IF
! ==================================================================== !
! Driver record 142 lamina_n_expansion. The record is RESERVED and its only
! accepted value is 0: v1.0 does not limit lamina expansion at a functional
! specific leaf nitrogen. The record is read so that the file stays positional,
! and any other value is refused below rather than ignored. It is retained here,
! and not deleted, because deleting a line of driver_main.txt would silently
! reassign every record below it.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_n_expansion
IF (ios /= 0) lamina_n_expansion = 0
! ==================================================================== !
! Driver record 143 n_status_construction_price: a WHOLE-CROP
! NITROGEN-STATUS-DEPENDENT CONSTRUCTION PRICE FOR EVERY STRUCTURAL ORGAN.
! 1 = each structural organ pass-1 price is interpolated LINEARLY, by the
! whole-shoot vegetative nitrogen nutrition index that is already computed
! unconditionally, between its existing floor and its existing record-128
! status-capped target. No new constant. Requires records 121 and 128 both at
! or above 1, since there is otherwise no floor and target pair to interpolate
! between. Read before the driver file is closed. The released value is 0.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_status_construction_price
IF (ios /= 0) n_status_construction_price = 0
IF (n_status_construction_price < 0 .OR. n_status_construction_price > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 143: n_status_construction_price must be 0 (OFF) or 1; read ', &
    n_status_construction_price
  ERROR STOP 'rec 143: invalid n_status_construction_price'
END IF
IF (n_status_construction_price >= 1 .AND. (n_payment_gate == 0 .OR. n_construction_status == 0)) THEN
  WRITE (*,'(A,I0,A,I0)') ' ERROR rec 143: n_status_construction_price >= 1 requires rec 121' // &
    ' n_payment_gate >= 1 AND rec 128 n_construction_status >= 1. Read rec 121 = ', n_payment_gate, &
    ', rec 128 = ', n_construction_status
  ERROR STOP 'rec 143 requires rec 121 >= 1 and rec 128 >= 1'
END IF
! ==================================================================== !
! Driver record 144 carbon_coupled_leaf_expansion: CARBON-COUPLED LEAF-AREA
! EXPANSION, the positive source-sink coupling. 1 = a per-axis per-rank bonus
! accumulator on the potential blade area, capped at 1.20 (Masle 2000
! Table II), accrues on days when the lamina sink is carbon-saturated with
! surplus AND water is non-limiting; the nitrogen-payment gate of record 121
! and the shed test of record 122 then price and arrest the extra area exactly
! as they do any other lamina area. Requires records 121 and 122 both at or
! above 1. Read before the driver file is closed. The released value is 0.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) carbon_coupled_leaf_expansion
IF (ios /= 0) carbon_coupled_leaf_expansion = 0
IF (carbon_coupled_leaf_expansion < 0 .OR. carbon_coupled_leaf_expansion > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 144: carbon_coupled_leaf_expansion must be 0 (OFF) or 1; read ', &
    carbon_coupled_leaf_expansion
  ERROR STOP 'rec 144: invalid carbon_coupled_leaf_expansion'
END IF
IF (carbon_coupled_leaf_expansion >= 1 .AND. (n_payment_gate == 0 .OR. lai_n_arrest == 0)) THEN
  WRITE (*,'(A,I0,A,I0)') ' ERROR rec 144: carbon_coupled_leaf_expansion >= 1 requires rec 121' // &
    ' n_payment_gate >= 1 AND rec 122 lai_n_arrest >= 1. Read rec 121 = ', n_payment_gate, &
    ', rec 122 = ', lai_n_arrest
  ERROR STOP 'rec 144 requires rec 121 >= 1 and rec 122 >= 1'
END IF
! ==================================================================== !
! Driver record 145 leaf_maintenance_basis: 1 removes the leaf
! maintenance-respiration term from the organ maintenance sum. No dependency on
! any other record. Read before the driver file is closed. The released value
! is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) leaf_maintenance_basis
IF (ios /= 0) leaf_maintenance_basis = 0
IF (leaf_maintenance_basis < 0 .OR. leaf_maintenance_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 145: leaf_maintenance_basis must be 0 (OFF) or 1; read ', &
    leaf_maintenance_basis
  ERROR STOP 'rec 145: invalid leaf_maintenance_basis'
END IF
! ==================================================================== !
! Driver record 146 par_reflection_basis: 1 removes the extra total-solar
! reflection factor from the two incident-PAR lines, so that the PAR reflection
! coefficient is applied exactly once. No dependency on any other record. Read
! before the driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) par_reflection_basis
IF (ios /= 0) par_reflection_basis = 0
IF (par_reflection_basis < 0 .OR. par_reflection_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 146: par_reflection_basis must be 0 (OFF) or 1; read ', &
    par_reflection_basis
  ERROR STOP 'rec 146: invalid par_reflection_basis'
END IF
! ==================================================================== !
! Driver record 147 lamina_n_area_retention: a NITROGEN LIMIT ON RETAINED GREEN
! LAMINA AREA. 1 = the retention criterion of growth_module.f90 replaces the
! minimum-specific-leaf-nitrogen criterion by the delivery path of record 122.
! Requires record 122 = 2 (the delivery path), record 121 at or above 1 (which
! supplies the leaf nitrogen pool) and record 84 = 1 (the single senescence
! clock). Read before the driver file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_n_area_retention
IF (ios /= 0) lamina_n_area_retention = 0
IF (lamina_n_area_retention < 0 .OR. lamina_n_area_retention > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 147: lamina_n_area_retention must be 0 (OFF), 1, 2, 3 or 4; read ', &
    lamina_n_area_retention
  ERROR STOP 'rec 147: invalid lamina_n_area_retention'
END IF
IF (lamina_n_area_retention >= 1 .AND. (lai_n_arrest /= 2 .OR. n_payment_gate < 1 .OR. lai_from_organs /= 1)) THEN
  WRITE (*,'(A,I0,A,I0,A,I0)') ' ERROR rec 147: lamina_n_area_retention >= 1 requires rec 122' // &
    ' lai_n_arrest = 2 AND rec 121 n_payment_gate >= 1 AND rec 84 lai_from_organs = 1. Read rec 122 = ', &
    lai_n_arrest, ', rec 121 = ', n_payment_gate, ', rec 84 = ', lai_from_organs
  ERROR STOP 'rec 147 requires rec 122 = 2, rec 121 >= 1 and rec 84 = 1'
END IF
! ==================================================================== !
! Driver record 148 wsc_fate_law: THE FATE OF THE STEM RESERVE. Term A is the
! refusal path of the deposition capacity and term B the release limb with its
! maturity clock. Requires record 92 at or above 1 (the pool itself), records
! 127 = 1 and 129 = 1 (the term-A refusal path), records 105 = 2 and 107 = 2
! (the term-B release limb and its maturity clock) and record 91 = 1 or 2 (the
! filling end bound that gates term B exists only under the rate x duration
! law); otherwise an error stop. Read before the driver file is closed. The
! released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) wsc_fate_law
IF (ios /= 0) wsc_fate_law = 0
IF (wsc_fate_law < 0 .OR. wsc_fate_law > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 148: wsc_fate_law must be 0 (OFF), 1 (the law), 2 (term A only)' // &
    ' or 3 (term B only); read ', wsc_fate_law
  ERROR STOP 'rec 148: invalid wsc_fate_law'
END IF
IF (wsc_fate_law >= 1) THEN
  IF (stem_wsc_mech < 1 .OR. repro_surplus_fate /= 1 .OR. sink_full_carbon_fate /= 1 .OR. &
      wsc_reserve_law /= 2 .OR. wsc_rel_window_member /= 2 .OR. &
      .NOT. (grain_fill_mech == 1 .OR. grain_fill_mech == 2)) THEN
    WRITE (*,'(A,6(A,I0))') ' ERROR rec 148: wsc_fate_law >= 1 requires rec 92 stem_wsc_mech >= 1' // &
      ' AND rec 127 repro_surplus_fate = 1 AND rec 129 sink_full_carbon_fate = 1 AND rec 105' // &
      ' wsc_reserve_law = 2 AND rec 107 wsc_rel_window_member = 2 AND rec 91 grain_fill_mech = 1 or 2.', &
      ' Read rec 92 = ', stem_wsc_mech, ', rec 127 = ', repro_surplus_fate, &
      ', rec 129 = ', sink_full_carbon_fate, ', rec 105 = ', wsc_reserve_law, &
      ', rec 107 = ', wsc_rel_window_member, ', rec 91 = ', grain_fill_mech
    ERROR STOP 'rec 148 requires rec 92 >= 1, rec 127 = 1, rec 129 = 1, rec 105 = 2, rec 107 = 2, rec 91 in {1,2}'
  END IF
END IF
! ==================================================================== !
! Driver record 149 wb_surface_limit: the water-balance conservation treatment
! of the Richards solver. The quantity it addresses is the per-sub-step
! Newton-Picard mass-balance error of hydro.f90, whose tolerances are SHARED
! module variables.
!   1 = tightened tolerances only (the step tolerance from 1e-6 to 1e-9 m and
!       the head tolerance at the minimum time step from 1.5e-5 to 1e-8 m),
!       with the relaxed zeroing of the error at the minimum time step skipped
!       so that the true error accumulates.
!   2 = mode 1 plus an exact post-step conservative closure.
!   3 = the released value, the CONSERVATIVE SCHEME: the mode-1 true-residual
!       ledger plus (i) one consistent top control volume (the first layer
!       thickness) in both the flux-boundary balance row and the storage
!       ledger, and (ii) no post-solve modification of the ledgered top flux,
!       so that the evaporation-demand cap and the infiltration-cap
!       re-application are reported rather than applied, the Picard-loop
!       boundary-condition shifts owning both limits. Numerics only, no new
!       physics, tolerances unchanged, no redistribution.
! Under any mode at or above 1 the reported surface balance uses the
! interception-net supply. No co-record dependency: this is water-balance
! conservation accounting and not a crop-mechanism switch. Read before the
! driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) wb_surface_limit
IF (ios /= 0) wb_surface_limit = 0
IF (wb_surface_limit < 0 .OR. wb_surface_limit > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 149: wb_surface_limit must be 0 (OFF), 1 (true-residual' // &
    ' ledger), 2 (+ post-step closure) or 3 (conservative scheme); read ', wb_surface_limit
  ERROR STOP 'rec 149: invalid wb_surface_limit'
END IF
! ==================================================================== !
! Driver record 150 pheno_dev_factor_window: 1 = Eq. (1) p.563 of McMaster et
! al. (2008) as cited, the day thermal time times the minimum of the
! vernalization and photoperiod factors, the two factors acting only from
! seedling emergence until terminal spikelet (400 degree-Cd of thermal
! development units after emergence, p.563). The numerics of the stage engine
! are otherwise unchanged, and the mode is identical to mode 0 on any deck with
! no photoperiod or vernalization sensitivity. Read before the driver file is
! closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) pheno_dev_factor_window
IF (ios /= 0) pheno_dev_factor_window = 0
IF (pheno_dev_factor_window < 0 .OR. pheno_dev_factor_window > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 150: pheno_dev_factor_window must be 0 (legacy) or 1' // &
    ' (McMaster 2008 window, min(VF,DF)); read ', pheno_dev_factor_window
  ERROR STOP 'rec 150: invalid pheno_dev_factor_window'
END IF
! ==================================================================== !
! Driver record 151 fe_trait_mode: 1 = the cultivar row fruiting_efficiency
! (grains per gram of spike dry weight at anthesis) sets the floret mass scale
! as the kernels per floret over that efficiency; the scale is applied after
! input/cultivar_yield.csv is read. The Slafer et al. (2015) sources are in
! SHARED.f90. Read before the driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) fe_trait_mode
IF (ios /= 0) fe_trait_mode = 0
IF (fe_trait_mode < 0 .OR. fe_trait_mode > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 151: fe_trait_mode must be 0 (Fischer constants) or 1' // &
    ' (cultivar fruiting-efficiency row); read ', fe_trait_mode
  ERROR STOP 'rec 151: invalid fe_trait_mode'
END IF
! Driver record 152 census_files: an output switch only. 0 = write the two
! hourly census files; 1 = neither file is opened. No process and no gate reads
! either file, and every other output is untouched in both modes. The released
! value is 1.
READ (1,*,IOSTAT=ios) census_files
IF (ios /= 0) census_files = 0
IF (census_files < 0 .OR. census_files > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 152: census_files must be 0 (write) or 1 (suppress); read ', census_files
  ERROR STOP 'rec 152: invalid census_files'
END IF
IF (census_files == 1) WRITE (*,'(A)') ' rec 152 census_files = 1: ci_clamp_census.csv and' // &
  ' spac_census.csv NOT written (output switch only)'
! ==================================================================== !
! Driver record 153 n_init_placement: INITIAL MINERAL NITROGEN ON ITS OWN
! SAMPLING DEPTH. 0 = the stated pool spread at uniform solution
! concentration over the WHOLE Richards profile; 1 = spread over the layers
! whose top lies above the deck's own reported sampling depth (optional row 71 of
! input/driver_N.csv, in m), proportional to the water content times the layer
! thickness within that depth. A deck that does not report row 71 is
! arithmetically identical to mode 0 even at mode 1: the record places a
! measurement on its measured support and invents no depth. Read before the
! driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_init_placement
IF (ios /= 0) n_init_placement = 0
IF (n_init_placement < 0 .OR. n_init_placement > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 153: n_init_placement must be 0 (whole profile) or 1' // &
    ' (declared sampling depth); read ', n_init_placement
  ERROR STOP 'rec 153: invalid n_init_placement'
END IF
! ==================================================================== !
! Driver record 154 organ_n_turnover_nema: POST-ANTHESIS ORGAN NITROGEN
! TURNOVER. 1 = the organ nitrogen turnover of Bertheloot, Cournede & Andrieu
! (2011) Ann Bot 108:1085-1096 and Bertheloot, Wu, Cournede & Andrieu (2011)
! Ann Bot 108:1097-1109 after anthesis: unconditional degradation of the
! remobilizable nitrogen of lamina, stem, chaff AND ROOT into the culm-common
! mobile pool, the day reduced nitrogen credited to the same pool, the grain
! taking the minimum of its sink and the pool, and concentration- and
! irradiance-driven protein synthesis in every photosynthetic entity from what
! remains. 2 = mode 1 without the root limb, so that the root contribution is
! separable. The alternative is the single-pool law, in which the Eq. (11)
! grain claim is pooled with the Eq. (6) lamina claim and the stem and chaff are
! degraded on demand to serve both. Read before the driver file is closed. The
! released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) organ_n_turnover_nema
IF (ios /= 0) organ_n_turnover_nema = 0
IF (organ_n_turnover_nema < 0 .OR. organ_n_turnover_nema > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 154: organ_n_turnover_nema must be 0 (OFF), 1 (Bertheloot 2011a turnover' // &
    ' incl. roots) or 2 (Bertheloot 2011a turnover, no root limb); read ', organ_n_turnover_nema
  ERROR STOP 'rec 154: invalid organ_n_turnover_nema'
END IF
! ==================================================================== !
! Driver record 155 n_uptake_regulation: REGULATION OF ROOT NITROGEN UPTAKE
! DEMAND.
!   1 = (a) + (b) + (c);  2 = (a) only;  3 = (b) + (c) only, where
!   (a) = the pre-anthesis demand driven by the MAXIMUM dilution curve of
!         Justes et al. (1994) in its marginal form, i.e. luxury uptake between
!         the critical and the maximum concentration,
!   (b) = the down-regulation of the whole-plant uptake demand by the culm
!         mobile-nitrogen concentration, E_N = exp(-130 [N]mob), Eq. (4) of
!         Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096,
!   (c) = the post-anthesis grain demand NET OF the mobile pool, i.e. the
!         record-97 mode-2 statement made reachable, applied only when record
!         154 is at or above 1 and holds the culm-common pool.
! Read before the driver file is closed. The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_uptake_regulation
IF (ios /= 0) n_uptake_regulation = 0
IF (n_uptake_regulation < 0 .OR. n_uptake_regulation > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 155: n_uptake_regulation must be 0 (OFF), 1 ((a)+(b)+(c)),' // &
    ' 2 ((a) only) or 3 ((b)+(c) only); read ', n_uptake_regulation
  ERROR STOP 'rec 155: invalid n_uptake_regulation'
END IF
! ==================================================================== !
! Driver record 156 n_construction_ceiling_curve: WHICH NITROGEN DILUTION CURVE
! CAPS THE ORGAN CONSTRUCTION-NITROGEN TARGETS. 0 = the pass-2 ceiling of
! record 128 keeps the CRITICAL dilution pair a = 5.35, b = 0.442; 1 = the same
! statement, the same marginal form a(1-b) W^-b and the same 1.55 t ha-1
! validity floor with the MAXIMUM dilution pair a = 8.3, b = 0.44, which are
! the model's own existing maximum-curve coefficients and are the STICS
! maximum-curve demand (Brisson et al. 1998 Agronomie 18:311-346 p.342
! Eq. (47) and p.329, with Table VIII p.330 attributing the pair to Justes et
! al. 1994). Inert when record 128 = 0. Read before the driver file is closed.
! The released value is 1.
! ==================================================================== !
READ (1,*,IOSTAT=ios) n_construction_ceiling_curve
IF (ios /= 0) n_construction_ceiling_curve = 0
IF (n_construction_ceiling_curve < 0 .OR. n_construction_ceiling_curve > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 156: n_construction_ceiling_curve must be 0 (critical pair' // &
    ' 5.35/0.442) or 1 (maximum pair 8.3/0.44); read ', n_construction_ceiling_curve
  ERROR STOP 'rec 156: invalid n_construction_ceiling_curve'
END IF
! ==================================================================== !
! Driver record 157 lamina_n_storage: LUXURY (STORAGE) NITROGEN IN THE GREEN
! LAMINA. 0 = the lamina tier of the pre-anthesis ladder keeps the construction
! target 1.5 g N m-2 as its ceiling.
!   1 = the lamina tier ceiling becomes the printed holding capacity
!       2.5 g N m-2 of leaf area (Pask 2009 PhD thesis, University of
!       Nottingham, p.29), which at this model leaf mass per area of
!       48.5 g DM m-2 is 5.155 %DM, on the SAME area basis and with the
!       incumbent order (lamina storage before stem and chaff luxury);
!   2 = the same ceiling with the stem and chaff LUXURY tier paid BEFORE the
!       lamina storage tier, the alternative order;
!   3 = the released value: mode 1 plus the record-155 limb-(b) mobile-nitrogen
!       NUMERATOR widened from the mobile vegetative pool alone to that pool
!       plus the shoot nitrate store plus the stem and chaff nitrogen above the
!       structural floor. The STORED LAMINA nitrogen is deliberately EXCLUDED:
!       Bertheloot, Cournede & Andrieu (2011) define the numerator as the
!       MOBILE common pool of amino acids and nitrate, and define the lamina
!       remobilizable nitrogen as only the nitrogen associated with the
!       photosynthetic apparatus, which is what the lamina store is in this
!       model, since it drives the specific leaf nitrogen and so Vcmax and
!       Jmax. No new constant.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_n_storage
IF (ios /= 0) lamina_n_storage = 0
IF (lamina_n_storage < 0 .OR. lamina_n_storage > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 157: lamina_n_storage must be 0 (OFF, ceiling = n4_sln_target' // &
    ' 1.5 g N m-2), 1 (ceiling = sln_max_N_lnh 2.5 g N m-2, lamina storage first),' // &
    ' 2 (same ceiling, stem luxury first) or 3 (mode 1 + stem/chaff luxury and the shoot' // &
    ' nitrate store in the rec-155 limb-(b) [N]mob numerator); read ', lamina_n_storage
  ERROR STOP 'rec 157: invalid lamina_n_storage'
END IF
! ==================================================================== !
! Driver record 158 nitrate_drain_rate_basis: WHAT DRIVES THE
! PHOTORESPIRATORY NITROGEN-ASSIMILATION ELECTRON DRAIN. 0 = the supply term is
! the standing shoot nitrate store divided by the photoperiod and the leaf area
! index, i.e. a stock asserted to be a flux.
!   1 = the released value: the supply term becomes the model's OWN daily de novo
!       nitrate-reduction DEMAND in the shoot - yesterday's nitrogen paid into
!       lamina, stem and chaff protein plus the grain draw on newly assimilated
!       nitrogen - through the SAME unit chain and under the SAME enzymatic
!       ceiling 1.21 umol m-2 s-1. Busch, Sage & Farquhar (2018) Nat Plants
!       4:46-54, p.49: the rate will be limited by the maximum rate of de novo
!       nitrogen supply to the chloroplast (Nmax), and Nmax can be limited by
!       the nitrogen availability, the demand for nitrogen in the form of
!       glycine and serine, or the activities of the enzymes involved in the
!       reduction of nitrate to ammonium. No new constant.
! There is no mode 2: the only fitted printed cap, 1.21, is already applied by
! the unchanged minimum in both modes, and the only other printed value (1.5,
! their Fig. 4) is an illustrative model input. Read before the driver file is
! closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) nitrate_drain_rate_basis
IF (ios /= 0) nitrate_drain_rate_basis = 0
IF (nitrate_drain_rate_basis < 0 .OR. nitrate_drain_rate_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 158: nitrate_drain_rate_basis must be 0 (OFF, Busch N_max' // &
    ' supply term = the standing shoot nitrate store / photoperiod / LAI) or 1 (supply term' // &
    ' = the day-before shoot de novo nitrate-reduction demand, same unit chain, same' // &
    ' N_max_ref ceiling).  There is no mode 2 (SHARED.f90); read ', nitrate_drain_rate_basis
  ERROR STOP 'rec 158: invalid nitrate_drain_rate_basis'
END IF
! ==================================================================== !
! Driver record 159 mobile_n_nema_estimate: the published ESTIMATE of the
! MOBILE NITROGEN POOL in the record-155 limb-(b) numerator.
!   1 = numerator = the mobile fraction times the remobilizable pools of
!       lamina, stem and chaff of record 154, with
!         N_rem,x = MAX(0, N_x - c_struct,x W_x) / (1 + mobile fraction),
!       and the structural concentrations 0.0079 (lamina) and 0.0035 (stem and
!       chaff) g N g-1 DM. Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot
!       108(6):1097-1109, doi:10.1093/aob/mcr125: remobilizable plus mobile
!       nitrogen mass allocated to a module was estimated as the difference
!       between total nitrogen mass and structural nitrogen mass, and mobile
!       nitrogen mass was estimated as in Bertheloot et al. (2008a), assuming
!       that it represented 10 % of the remobilizable nitrogen mass. The
!       estimate REPLACES the numerator, because it is an estimate OF that
!       single mobile pool and not an addition to it.
!   2 = mode 1 plus the model explicit transport pools, the alternative which
!       double-counts that same pool.
! No new numeric constant: the mobile fraction is a parameter alias of the
! existing 0.10 of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781
! Eq. (1) p.783, value p.786, which is the constant the 2011 paper cites. Read
! before the driver file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mobile_n_nema_estimate
IF (ios /= 0) mobile_n_nema_estimate = 0
IF (mobile_n_nema_estimate < 0 .OR. mobile_n_nema_estimate > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 159: mobile_n_nema_estimate must be 0 (OFF, the rec-155 limb-(b)' // &
    ' [N]mob numerator is whatever rec 157 formed), 1 (numerator = n2_mobile_frac * the rec-154' // &
    ' remobilizable N_ph pools of lamina + stem + chaff -- Bertheloot 2011b''s printed 10 % estimate, the' // &
    ' faithful arm) or 2 (mode 1 + N_mobile_veg + N_NO3_shoot, the disclosed double-counting' // &
    ' alternative); read ', mobile_n_nema_estimate
  ERROR STOP 'rec 159: invalid mobile_n_nema_estimate'
END IF
! ==================================================================== !
! Driver record 160 justes_nmax_floor_basis: WHAT THE MAXIMUM NITROGEN DILUTION
! CURVE IS WORTH BELOW ITS PRINTED VALIDITY FLOOR.
!   0 = the outer cap is held at a constant 6.6 % while the shoot is at or
!       below 1.55 t ha-1. That 6.6 % is NOT in Justes et al. (1994) Ann Bot
!       74:397-407: the paper prints Eq. (17) p.405, Nmax = 8.3 DM^-0.44,
!       valid from 1.5 to 14 t ha-1 (Fig. 5, p.404 col. 2 to p.405 col. 1), and
!       no maximum-concentration statement at all below 1.5 t ha-1.
!   1 = the released value: the outer cap below the floor is held at the
!       printed curve's own value at the floor, 8.3 x 1.55**(-0.44) = 6.844 %,
!       computed from the coefficients this file already carries. This is the
!       same construction the paper applies to the CRITICAL curve (a horizontal
!       line below the break point, p.402 col. 2) and the same construction
!       records 128, 156 and 155 limb (a) already use; the level is set by
!       continuity because the paper states no level to quote.
! There is no mode 2: the only other defensible floor is the printed validity
! floor 1.5 t ha-1 itself (6.943 %, 1.4 % higher), and the build keeps one
! floor, 1.55, the paper's own break point. Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) justes_nmax_floor_basis
IF (ios /= 0) justes_nmax_floor_basis = 0
IF (justes_nmax_floor_basis < 0 .OR. justes_nmax_floor_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 160: justes_nmax_floor_basis must be 0 (OFF, legacy 6.6 %' // &
    ' Nmax plateau below 1.55 t ha-1, a constant that is NOT in Justes et al. 1994) or 1' // &
    ' (the Nmax outer cap below the floor held at the printed curve''s own value there,' // &
    ' 8.3 * 1.55**(-0.44); SHARED.f90).  There is no mode 2; read ', &
    justes_nmax_floor_basis
  ERROR STOP 'rec 160: invalid justes_nmax_floor_basis'
END IF
! ==================================================================== !
! Driver record 161 canopy_clock_driver.
!   0 = the released value: the canopy clock - and with it leaf appearance, the
!       final leaf number, blade expansion and the tiller-birth windows - is
!       advanced by the day SCREEN-AIR degree-days.
!   1 = from emergence until the model's own terminal-spikelet marker the canopy
!       clock is advanced by the degree-days of the TOP SOIL LAYER, the layer
!       that contains the measured 3 cm on the deck grids, and by screen air
!       thereafter. Source, read at first hand: Abichou (2016) PhD thesis HAL
!       tel-02316154, p.23 (soil temperature at 3 cm depth; the apex
!       temperature was estimated based on soil temperature before stem
!       extension and then based on air temperature) with the justification on
!       p.38. Zero new constants.
! THE PRIMARY FIELD FALSIFIER of this record, read at first hand: McMaster,
! Wilhelm et al. (2003) "Spring wheat leaf appearance and temperature:
! extending the paradigm?", Ann Bot 91(6):697-705, doi:10.1093/aob/mcg074.
! Heating the soil at crown depth by 3 C did NOT increase the leaf appearance
! rate. p.700: a mean phyllochron of 74 air growing degree days (their Table 2)
! would suggest about 2.7, 2.3 and 1.8 more leaves appearing in the +3 C
! treatments, or alternatively, for a constant final leaf number, the flag leaf
! would have appeared about 65, 55 or 45 d earlier in the three successive
! planting-date treatments; this was not observed. And p.702, opening sentence
! of the Discussion: in contrast to the prediction of the paradigm, heating the
! soil near crown depth did not increase the rate of leaf appearance of spring
! wheat under field conditions. That is a direct field falsification of a
! soil-driven canopy clock over exactly this record window, since the +3 C
! treatment raised the soil temperature at crown depth, confirmed by the
! earlier seedling emergence it produced, and leaf appearance did not follow.
! Nagelmueller et al. (2016) concurs. This is why record 161 is 0 in the
! released configuration: the record exists so that the soil-temperature
! hypothesis is testable and its falsifier is on the record.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) canopy_clock_driver
IF (ios /= 0) canopy_clock_driver = 0
IF (canopy_clock_driver < 0 .OR. canopy_clock_driver > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 161: canopy_clock_driver must be 0 (OFF, the canopy clock runs on' // &
    ' screen-air degree-days as before) or 1 (the canopy clock runs on Tsoil_prof(1) degree-days from' // &
    ' emergence until terminal spikelet, then on screen air; SHARED.f90).  read ', &
    canopy_clock_driver
  ERROR STOP 'rec 161: invalid canopy_clock_driver'
END IF
! ==================================================================== !
! Driver record 162 canopy_clock_basis.
!   0 = the canopy clock is advanced by the RAW degree-day increment.
!   1 = the released value: the canopy clock is advanced by the
!       vernalization- and photoperiod-gated development increment the stage
!       engine already uses, so that development and canopy run on ONE time
!       basis instead of two.
! The law is McMaster et al. (2008) Ann Bot 102:561-569,
! doi:10.1093/aob/mcn115, Eq. (1) p.563 and Eqs. (2) to (4) pp.564-565. It is
! an INFERENCE and not printed architecture: the printed scope of that paper is
! stage progression from emergence to terminal spikelet and it does not state
! that leaf appearance is gated by those two factors. The phyllochron is held
! at 73.9 degree-Cd and no constant moves.
! ==================================================================== !
READ (1,*,IOSTAT=ios) canopy_clock_basis
IF (ios /= 0) canopy_clock_basis = 0
IF (canopy_clock_basis < 0 .OR. canopy_clock_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 162: canopy_clock_basis must be 0 (OFF, the canopy clock is advanced by' // &
    ' the raw increment current_dtt as before) or 1 (advanced by the VF/DF-gated development increment' // &
    ' current_dtdu the stage engine already uses; SHARED.f90).  read ', &
    canopy_clock_basis
  ERROR STOP 'rec 162: invalid canopy_clock_basis'
END IF
! ==================================================================== !
! Driver record 163 postanth_uptake_regulation: AFTER ANTHESIS, DOES THE GRAIN
! NITROGEN REQUIREMENT REACH THE ROOT?
!   0 = the released value: the post-anthesis uptake potential is the minimum
!       of the construction demand and the supply, the demand being the
!       vegetative demand plus the post-anthesis grain demand. That second term
!       is 0.0000 on EVERY post-anthesis day at the Braunschweig high-nitrogen
!       decks while 44 to 57 kg N ha-1 of mineral nitrogen stands in the soil,
!       so the grain is paid for by stripping the lamina.
!   1 = after anthesis the DEMAND GATE IS REMOVED and the uptake potential is
!       regulated as Bertheloot, Cournede & Andrieu (2011) Ann Bot
!       108:1085-1096, doi:10.1093/aob/mcr119, Eqs. (3) to (5) print it:
!       potential = supply x E_N, with the supply the unchanged existing
!       soil-nitrogen supply law and E_N the unchanged record-155 limb-(b)
!       value. No new numeric constant: the record deletes a gate and adds
!       nothing.
! The carbon term E_C of Eq. (5) is NOT implemented. Its form is no longer
! inferred - the typeset paper is on disk and Eq. (5) was read at first hand on
! p.1090, E_C = 0 if the root remobilizable mass is at or below its minimum and
! 1 - exp(-beta_C x that mass) otherwise - but its argument is a per-culm root
! dry-mass influx in g d-1 with a 0.0013 g d-1 threshold, while this model
! carries only a per-ground-area root increment. Since E_C is at most 1, mode 1
! is an UPPER BOUND on the uptake of that law.
! The mode-1 run.log banner further below still prints the superseded wording
! about the typeset form. It is an executable WRITE statement, so it was not
! edited in this documentation pass; record 163 is 0 in the released
! configuration, so the line is never printed on any shipped deck.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) postanth_uptake_regulation
IF (ios /= 0) postanth_uptake_regulation = 0
IF (postanth_uptake_regulation < 0 .OR. postanth_uptake_regulation > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 163: postanth_uptake_regulation must be 0 (OFF, post-anthesis' // &
    ' U_pot = MIN(construction demand, supply)) or 1 (after anthesis the demand gate is removed and' // &
    ' U_pot = the existing soil-N supply law x the record-155 E_N, as Bertheloot 2011a Eqn 3-5 prints it with' // &
    ' no demand term; SHARED.f90); read ', postanth_uptake_regulation
  ERROR STOP 'rec 163: invalid postanth_uptake_regulation'
END IF
! ==================================================================== !
! Driver record 164 lamina_capacity_profile_basis: AT WHICH LEVEL OF THE CANOPY
! DOES THE 2.5 g N m-2 LAMINA HOLDING CAPACITY APPLY?
!   0 = the released value: the record-157 ceiling is applied as the
!       CANOPY-MEAN lamina capacity, leaf nitrogen over leaf area index. The
!       model's own output/daily_canopy_n_profile.csv then prints a canopy mean of
!       exactly 2.5000 before anthesis with a top value of 3.20-3.24, against a
!       MEASURED canopy-mean specific leaf nitrogen at anthesis of 1.65-1.90
!       (Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707-1720 Fig. 2A
!       p.1711) and 1.73-2.11 at the Braunschweig high-nitrogen decks.
!   1 = the SAME printed 2.5 g N m-2 applied at the level the source sentence
!       names - especially in the UPPER LEAVES of the canopy (Pask 2009 PhD
!       thesis, University of Nottingham, p.29) - i.e. at the TOP of the model
!       OWN existing exponential canopy nitrogen profile, so that the
!       canopy-mean capacity each day becomes 2.5 divided by the ratio of the
!       top value to the mean, taken from that profile's own daily state. No new
!       constant: the top value is the same pure function that writes the
!       n_top_g_m2 column and the mean is leaf nitrogen over leaf area index.
!       The factor is at least 1, so the capacity can only FALL, from 2.50
!       toward about 1.93-1.95 at anthesis, inside the measured band. On a day
!       with a degenerate profile (a vanishing leaf area index, leaf nitrogen
!       or extinction coefficient, or a computed ratio below 1) the capacity
!       falls back to the mode-0 value and the day is counted.
! Inert when record 157 = 0, in which case the lamina storage tier runs on the
! construction target; reported in the banner rather than refused.
! Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_capacity_profile_basis
IF (ios /= 0) lamina_capacity_profile_basis = 0
IF (lamina_capacity_profile_basis < 0 .OR. lamina_capacity_profile_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 164: lamina_capacity_profile_basis must be 0 (OFF, the record-157' // &
    ' ceiling 2.5 g N m-2 applied as the canopy MEAN lamina capacity) or 1 (the same printed' // &
    ' 2.5 g N m-2 applied at the TOP of the model''s own exponential canopy N profile, so the' // &
    ' canopy-mean capacity = 2.5 / (n_top/n_mean); SHARED.f90); read ', &
    lamina_capacity_profile_basis
  ERROR STOP 'rec 164: invalid lamina_capacity_profile_basis'
END IF
! ==================================================================== !
! Driver record 165 wsc_release_delay_basis: MAY THE STEM RESERVE BE MOBILIZED
! BEFORE ITS OWN ACCUMULATION WINDOW HAS CLOSED?
!   0 = the record-92 deficit draw limb carries only the anthesis condition, so
!       it opens on the FIRST post-anthesis day with any unmet grain demand -
!       3 to 8 d after anthesis on 16 of 16 Braunschweig decks, with the pool
!       at 0.00 g m-2 by 18 to 23 days post-anthesis, against a measurement in
!       which the reserve is still RISING at the +18/20 d sample on 15 of 16
!       decks (Manderscheid et al. 2020 ODJAR 6,
!       doi:10.18174/odjar.v6i0.16397).
!   1 = the released value: the deficit limb additionally requires the EXISTING
!       conjunct that the day be past the end of the accumulation window - the
!       identical conjunct the record-105 mode-2 release limb already carries -
!       so that mobilization cannot begin before the printed accumulation
!       window ends. Schnyder (1993) New Phytol 123:233-245,
!       doi:10.1111/j.1469-8137.1993.tb03731.x, p.238: water-soluble
!       carbohydrates accumulate in the stem until 10-20 d after anthesis and,
!       at that stage, the maximum amount is 2-4 times greater than at
!       anthesis; Summary p.233, accumulation usually terminates within three
!       weeks following anthesis, after which fructan is gradually lost until
!       grain maturity. Second, independent first-hand witness, measured
!       internode by internode: Ehdaie, Alloush, Madore & Waines (2006) Crop
!       Sci 46:735-746, doi:10.2135/cropsci2005.04-0033, Table 2 p.738 with the
!       text on p.737 - lower internodes peak at 10 d post-anthesis (822 mg),
!       peduncle (389 mg) and penultimate (417 mg) at 20 d, then progressively
!       decreased over time until maturity - and p.744, maximum stem weight in
!       all genotypes except one was attained 20 d post-anthesis, with Rawson
!       and Evans (1971) reporting 7 to 21 d. Third: Ruuska et al. (2006) Funct
!       Plant Biol 33:799-809, doi:10.1071/FP06062, p.803, the store increases
!       from before anthesis until about midway through grain filling and then
!       decreases again. Fourth: Gebbing, Schnyder & Kuhbauch (1998) Plant Cell
!       Environ 21:301-313, p.302, storage and mobilization mainly occur in the
!       internodes and leaf sheaths and accumulation may continue until 10-20 d
!       after anthesis, and p.310, the loss of stem weight from the time of
!       maximum content, at about 10-20 d after anthesis, to maturity was
!       closely related to that mobilization.
!       ZERO NEW CONSTANTS: the 15-day accumulation end is the EXISTING
!       record-92 parameter with its own printed provenance, the midpoint of
!       that measured 10-20 d. The draw inside the limb, its destination and
!       the conversion efficiency 0.95 (Gebbing & Schnyder 1999 Plant Physiol
!       121:871-878, doi:10.1104/pp.121.3.871, p.877, a very high efficiency of
!       0.97 +/- 0.07 of mobilized reserve utilization in grain filling) are
!       unchanged.
! Inert where the reserve mode is off or the deficit limb is already
! unreachable; reported in the banner rather than refused.
! ==================================================================== !
READ (1,*,IOSTAT=ios) wsc_release_delay_basis
IF (ios /= 0) wsc_release_delay_basis = 0
IF (wsc_release_delay_basis < 0 .OR. wsc_release_delay_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 165: wsc_release_delay_basis must be 0 (OFF, the rec-92 deficit' // &
    ' draw limb opens at anthesis with no post-anthesis delay) or 1 (the limb additionally requires' // &
    ' the EXISTING conjunct (iday_in - wsc_day_anth) > wsc_accum_end_dpa, so mobilization cannot' // &
    ' start before the printed 10-20 d accumulation window ends; SHARED.f90); read ', &
    wsc_release_delay_basis
  ERROR STOP 'rec 165: invalid wsc_release_delay_basis'
END IF
! ==================================================================== !
! Driver record 166 grain_ca_xylem_throughput: THE GRAIN CALCIUM ROUTE.
!   0 = grain calcium is the product of two constants, the grain water content
!       per unit dry matter and the calibrated xylem calcium concentration, so
!       no CO2, water, nitrogen, cultivar or season response is possible.
!   1 = the released value, the ear xylem-throughput route: calcium reaches
!       grain and chaff with the xylem water the EAR receives (the ear share of
!       the model's own standing green area times transpiration), split at the
!       grain-base xylem discontinuity by the measured partition of Hocking
!       (1994) Table 1(c) p.1294, with the campaign-calibrated xylem
!       concentration NOT READ: calcium falls back to the sourced
!       calcite-equilibrium 54.9 mg L-1.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) grain_ca_xylem_throughput
IF (ios /= 0) grain_ca_xylem_throughput = 0
IF (grain_ca_xylem_throughput < 0 .OR. grain_ca_xylem_throughput > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 166: grain_ca_xylem_throughput must be 0 (OFF, grain Ca =' // &
    ' grain_h2o_g_gDM x xylem_Ca_mgL_CAL, a product of two constants) or 1 (the ear' // &
    ' xylem-throughput route, SHARED.f90); read ', grain_ca_xylem_throughput
  ERROR STOP 'rec 166: invalid grain_ca_xylem_throughput'
END IF
! ==================================================================== !
! Driver record 167 lamina_synthesis_sheath_basis: ON WHICH MASS DOES THE
! LAMINA SYNTHESIS RATE ACT?
!   0 = the lamina rate 1.5e-4 multiplies the GREEN LAMINA MASS ALONE and the
!       stem rate 3.0e-5 multiplies the whole lumped stem mass, which CONTAINS
!       the leaf sheaths. Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot
!       108:1097-1109, doi:10.1093/aob/mcr125, Table 2 p.1101 prints 0.00015*
!       and 0.00003+ for that rate with its own footnotes on the same page:
!       "* Laminae and sheaths." and "+ Internodes, peduncle and chaff." Mode 0
!       therefore has the sheath under the wrong footnote. The cost is an
!       analytic fixed point in green-lamina nitrogen percentage at 2.22-2.26 %
!       on the +16-18 d Braunschweig sample (model 2.23-2.27, measured
!       3.42-3.65) whose ceiling at both saturating factors equal to 1 is
!       2.85 %, i.e. the measurement is unreachable on this basis at any
!       factor.
!   1 = the basis becomes the green lamina mass plus the model's OWN green sheath
!       area index (record 68) priced at its OWN green-lamina mass per leaf
!       area. This is a LOWER BOUND: the areas printed in Table 1c p.1100 make
!       the measured sheath 4.1 times the lamina mass per area.
!   2 = the released value: the basis becomes the green lamina mass times
!       (1 + 272/420), the PRINTED MEASURED dry masses at flowering of Table 1c
!       p.1100 (sheath structural 103+69+43+9 = 224 plus remobilizable
!       22+15+9+2 = 48, i.e. 272 mg culm-1; lamina 320 + 100 = 420 mg).
! In BOTH modes 1 and 2 the same sheath mass is REMOVED from the stem basis, so
! no gram synthesizes photosynthetic nitrogen twice and the lamina plus stem
! dry mass is still partitioned exactly once. No constant is changed.
! Pre-registered for mode 2 before the run: canopy-mean specific leaf nitrogen
! at +16-18 d of 1.62/1.64/1.66/1.65 g N m-2 at the four Braunschweig decks
! (then 1.122/1.132/1.140/1.119; measured 2.11/1.81/1.99/1.80), i.e. about half
! the gap closed and still below the observation.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_synthesis_sheath_basis
IF (ios /= 0) lamina_synthesis_sheath_basis = 0
IF (lamina_synthesis_sheath_basis < 0 .OR. lamina_synthesis_sheath_basis > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 167: lamina_synthesis_sheath_basis must be 0 (OFF, the Bertheloot 2011a sigma_lam' // &
    ' acts on the green lamina mass alone and sigma_stem on the whole lumped stem), 1 (sheath mass' // &
    ' from the model''s own green sheath area index cse_sai x W_green/lai) or 2 (sheath mass from' // &
    ' the printed measured 272/420 of Bertheloot 2011b Table 1c p.1100); SHARED.f90; read ', &
    lamina_synthesis_sheath_basis
  ERROR STOP 'rec 167: invalid lamina_synthesis_sheath_basis'
END IF
! ==================================================================== !
! Driver record 168 mineral_impedance_basis: the solute diffusion impedance
! factor inside mineral_module.f90, and nowhere else.
!   0 = f_l(theta) = MAX(0, 1.60 theta - 0.17), the liquid impedance factor
!       of Barraclough & Tinker (1981) J Soil Sci 32:225-236,
!       doi:10.1111/j.1365-2389.1981.tb01702.x, Fig. 1, whose printed best-fit
!       line for the sieved Woburn soil is 1.58 theta - 0.17 (the intercept is
!       exact, the slope 1.3 % high) and for the Rothamsted soil
!       0.99 theta - 0.17. The pair is
!       also SITE-BLIND: the same slope and intercept on the loamy sand, the
!       sand and the clay loam, so the whole between-site contrast in diffusive
!       mineral supply is carried by soil water alone.
!   1 = the released value: Olesen, Moldrup, Yamaguchi & Rolston (2001) Soil
!       Science 166:89-96, doi:10.1097/00010694-200102000-00002, p.94
!       Eqs. (7a) and (7b), f_l = 1.1 (theta - theta_th) for theta above the
!       threshold and 0 below, with the threshold from that paper's own p.93
!       Eq. (5) and the deck's OWN measured clay fraction, silt fraction and
!       topsoil bulk density (input/driver_N.csv positions 68, 69 and 70). No
!       new numerals: the slope and the six Eq. (5) coefficients are the
!       constants record 139 mode 3 already carries. The nitrogen path of
!       record 139 is not touched.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_impedance_basis
IF (ios /= 0) mineral_impedance_basis = 0
IF (mineral_impedance_basis < 0 .OR. mineral_impedance_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 168: mineral_impedance_basis must be 0 (OFF, the DECLARED AND' // &
    ' pair f_l = MAX(0, 1.60*vwc_root - 0.17), re-labeled 2026-10-18 to Barraclough & Tinker' // &
    ' 1981 J Soil Sci 32:225 p.234 Fig. 1 best fit 1.58*theta - 0.17, intercept exact / slope' // &
    ' 1.3 % off; record-168 conflict disclosed) or 1 (Olesen et al. 2001 Soil Sci' // &
    ' 166:89-96 p.94 Eq.(7a)/(7b) CSIF, f_l = 1.1*MAX(0, vwc_root - theta_th), SHARED.f90' // &
    '); read ', mineral_impedance_basis
  ERROR STOP 'rec 168: invalid mineral_impedance_basis'
END IF
! ==================================================================== !
! Driver record 169 zn_uptake_kinetic_basis: ON WHAT CONCENTRATION DOES THE
! ZINC TRANSPORTER SATURATE?
!   0 = the zinc limb keeps the transporter supply at the minimum of its raw
!       value and the sum of mass flow and diffusion. That minimum binds on
!       7499 of 7499 active deck-days across all 42 decks, so it is not a bound
!       but an identity, and the day raw zinc supply is exactly twice the sum
!       of mass flow and diffusion.
!   1 = the released value, the surface-concentration kinetic: the transporter
!       supply becomes the field capacity rate times the root mass times the
!       root activity times c_surf / (Km + c_surf), with c_surf the bulk
!       concentration minus the depletion the diffusion limb has already
!       charged, and no minimum against mass flow plus diffusion. The pairing
!       of the FIELD capacity of Liu et al. (2019) with the SEEDLING
!       half-saturation constant of Hacisalihoglu et al. (2001) is a model
!       construction and not a measurement.
!       DESCRIBED AS IT BEHAVES: the diffusion gradient IS the concentration
!       minus that same half-saturation constant, so on every day with
!       diffusion active the root-surface concentration is PINNED at the
!       half-saturation constant and the saturation factor is exactly 0.5; it
!       is 1.0 on the days the impedance factor is zero, and all of the term
!       variation comes from those days. The law is therefore the field
!       capacity, halved where diffusion is active and whole where it is not -
!       a model construction rather than a Michaelis-Menten kinetic.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) zn_uptake_kinetic_basis
IF (ios /= 0) zn_uptake_kinetic_basis = 0
IF (zn_uptake_kinetic_basis < 0 .OR. zn_uptake_kinetic_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 169: zn_uptake_kinetic_basis must be 0 (OFF, S_hats =' // &
    ' MIN(S_hats_raw, S_mf + S_diff), a MIN that binds on 7499/7499 active deck-days) or 1' // &
    ' (the root-surface-concentration kinetic, SHARED.f90); read ', zn_uptake_kinetic_basis
  ERROR STOP 'rec 169: invalid zn_uptake_kinetic_basis'
END IF

! ==================================================================== !
! Driver record 170 is a RESERVED row. input/driver_main.txt is POSITIONAL, so
! the line is read as an integer and has no effect anywhere, which keeps the
! numbers of records 171 and 172 fixed. End-of-file default 0, and no range
! check: the legal range of a reserved row belongs to whatever builds it. The
! process it is held for is a temperature correction of the mineral diffusion
! coefficients, which are quoted at about 25 C and carry none.
! ==================================================================== !
READ (1,*,IOSTAT=ios) reserved_rec_170
IF (ios /= 0) reserved_rec_170 = 0
! ==================================================================== !
! Driver record 171 mineral_n_coupling_basis: THE GRAIN PROTEIN SINK FOR
! PHLOEM-MOBILE MINERALS.
!   0 = mineral_module.f90 reads no nitrogen state on this path, so grain iron,
!       zinc, copper and magnesium can respond to nitrogen supply only through
!       the biomass denominator, i.e. only by DILUTION, and 13 of the 16
!       Braunschweig nitrogen-dose concentration contrasts have the wrong sign.
!   1 = the released value, FLUX (daily) basis: the day grain loading of
!       element e is bounded by its stoichiometry times the positive part of
!       the day grain-nitrogen increment.
!   2 = POOL (cumulative) basis: the accumulated grain content of element e is
!       bounded by its stoichiometry times the grain nitrogen of the day.
! The stoichiometry, in mg element per g grain nitrogen, is column 20 of
! input/driver_minerals.csv, read only when this record is at or above 1; it
! ships at 0.0, i.e. NO BOUND, for all six elements, because no first-hand
! wheat-grain primary source on this machine prints the stoichiometry. Xia et
! al. (2020) Food Energy Secur 9:e243, doi:10.1002/fes3.243, pp.6-7 and p.10.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_n_coupling_basis
IF (ios /= 0) mineral_n_coupling_basis = 0
IF (mineral_n_coupling_basis < 0 .OR. mineral_n_coupling_basis > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 171: mineral_n_coupling_basis must be 0 (OFF, the mineral module' // &
    ' reads no nitrogen state and the grain N-dose mineral response is dilution only), 1 (FLUX' // &
    ' basis: the day''s grain mineral loading <= grain_eN_mg_gN x the day''s grain N increment) or' // &
    ' 2 (POOL basis: the grain''s accumulated mineral <= grain_eN_mg_gN x grain N);' // &
    ' SHARED.f90; read ', mineral_n_coupling_basis
  ERROR STOP 'rec 171: invalid mineral_n_coupling_basis'
END IF

! ==================================================================== !
! Driver record 172 stem_reserve_deposition_n_basis: IS THE PRE-ANTHESIS STEM
! RESERVE DEPOSITION SHARE NITROGEN-BLIND?
!   0 = the released value: the pre-anthesis bank keeps the constant record-56
!       deposition share 0.20, so the realized at-anthesis stem reserve
!       concentration is 0.1914-0.1990 g g-1 on 16 of 16 Braunschweig decks and
!       its slope on nitrogen is exactly zero, against a measured
!       -0.68 mg g-1 per kg N ha-1 in two independent seasons.
!   1 = the deposition share becomes (a + b NNI) x 0.001, clamped to the
!       EXISTING record-56 bounds [0.0, 0.35], with the nitrogen nutrition
!       index on the critical curve of Justes et al. (1994) the model already
!       computes. The two constants are the unweighted means of the three site
!       regressions of van Herwaarden (1997) PhD thesis Ch.4, Tables 4.7
!       (p.108), 4.9 (p.110) and 4.11 (p.113); no numeral is fitted to any
!       campaign observation.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) stem_reserve_deposition_n_basis
IF (ios /= 0) stem_reserve_deposition_n_basis = 0
IF (stem_reserve_deposition_n_basis < 0 .OR. stem_reserve_deposition_n_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 172: stem_reserve_deposition_n_basis must be 0 (OFF, the' // &
    ' pre-anthesis WSC bank keeps the constant record-56 stem_reserve_frac) or 1 (the van' // &
    ' Herwaarden 1997 PhD Ch.4 NNI law, SHARED.f90); read ', &
    stem_reserve_deposition_n_basis
  ERROR STOP 'rec 172: invalid stem_reserve_deposition_n_basis'
END IF

! Record 172 dependency, refused rather than left a silent no-op: the statement
! the law replaces is live only when record 55 = 1. With record 55 = 0 the
! pre-anthesis bank does not exist, so mode 1 would change nothing and the run
! would report a repair it never made.
IF (stem_reserve_deposition_n_basis >= 1 .AND. carbon_reserve_mode /= 1) THEN
  WRITE (*,'(A,I0,A)') ' ERROR rec 172: stem_reserve_deposition_n_basis = 1 replaces the' // &
    ' pre-anthesis WSC deposition share, which exists only at rec 55 carbon_reserve_mode = 1;' // &
    ' read rec 55 = ', carbon_reserve_mode, '. Set rec 55 = 1 or rec 172 = 0.'
  ERROR STOP 'rec 172 requires rec 55 carbon_reserve_mode = 1'
END IF

! ==================================================================== !
! Driver record 173 spac_sink_coupling: THE RICHARDS SINK IS FORCED BY THE
! COMMITTED SOIL-PLANT-ATMOSPHERE LEAF.
!   0 = the Richards root-sink demand on every outer pass at or above 2 is the
!       record-78 UNSTRESSED hourly series from the outer-iteration-1
!       plain-Leuning pass, while the leaf the model commits is the
!       soil-plant-atmosphere leaf, whose transpiration is larger; the
!       difference is written as a shortfall and never fed back.
!   1 = the released value: on outer passes at or above 2 the sink demand IS
!       this pass's own converged hourly canopy transpiration and the daily
!       ledger potential is that series' own sum. Nothing is scaled, capped or
!       tuned, and the record-78 unstressed series stays written as the
!       Leuning-branch potential diagnostic it is.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) spac_sink_coupling
IF (ios /= 0) spac_sink_coupling = 0
IF (spac_sink_coupling < 0 .OR. spac_sink_coupling > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 173: spac_sink_coupling must be 0 (OFF, Richards is forced' // &
    ' with rec 78''s unstressed series) or 1 (forced with the committed SPAC pass''s own' // &
    ' hourly transpiration, SHARED.f90); read ', spac_sink_coupling
  ERROR STOP 'rec 173: invalid spac_sink_coupling'
END IF
! ==================================================================== !
! Driver record 174 canopy_fln_prediction: A BLADE FINAL DIMENSIONS ARE SET
! ONCE, FROM THE CULM'S OWN FINAL LEAF NUMBER.
!   0 = the released value: the blade-dimension profiles are indexed at a
!       preliminary top of 11 until the final leaf number is frozen at booting,
!       and on that day the FINAL length and width of every already-ligulated
!       blade are rewritten. At the Braunschweig deck with a final leaf number
!       of 7, seven finished blades are promoted four profile positions in one
!       day: the potential shoot area from 36.2047 to 72.1211 cm2, the leaf
!       area index from 1.166 to 3.088 on one day, 54.8 % of the season lamina
!       area, a specific leaf area of 38.7 m2 kg-1 against the model's own
!       18.5-20.6 and a specific leaf nitrogen of 0.865 g N m-2 against every
!       measured 1.96-2.75.
!   1 = the final leaf number is predicted on the EMERGENCE day as the integer
!       part of the closed-form Haun stage at the booting clock reading, and
!       used as the profile top from emergence; the booting freeze only
!       CONFIRMS it and records any difference. Zero new constants: the closed
!       form is the continuous solution of the SAME phyllochron law the
!       chronology integrates.
! Pre-registered before the run: the prediction returns 7 on all 16
! Braunschweig decks, 8-10 at Maricopa, 9-10 on the AgMIP decks, 10 and 12 on
! the two 2022/23 WheatDryFACE decks, 9 on the 2024/25 Zahir decks and 11 on the
! 2024/25 Ruta decks, so exactly 4 of 42 decks are unchanged by construction.
! ==================================================================== !
READ (1,*,IOSTAT=ios) canopy_fln_prediction
IF (ios /= 0) canopy_fln_prediction = 0
IF (canopy_fln_prediction < 0 .OR. canopy_fln_prediction > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 174: canopy_fln_prediction must be 0 (OFF, the provisional' // &
    ' panel-median top cba_Nl_fit = 11 stands until the Booting freeze, which then rewrites the final' // &
    ' dimensions of every finished blade) or 1 (the final leaf number is predicted on the emergence day' // &
    ' and used as n_top from emergence); SHARED.f90; read ', canopy_fln_prediction
  ERROR STOP 'rec 174: invalid canopy_fln_prediction'
END IF
! ==================================================================== !
! Driver record 175 grain_n_partition_basis: WHERE THE GRAIN MISSING NITROGEN
! SITS. The modes are CUMULATIVE.
!   1 = restore the post-anthesis continuity the header of record 158 declares
!       and the record-154 limb breaks (the newly assimilated nitrogen reaching
!       the grain being zero by construction), so that the previous-day
!       reduction rate is no longer identically zero after anthesis and the
!       shoot can still reduce nitrate. Zero new constants; all five terms are
!       already computed and already printed.
!   2 = mode 1 plus the MEASURED initial condition of 2.4 mg grain nitrogen per
!       ear (Bertheloot, Martre & Andrieu 2008 Funct Plant Biol 35:781-796
!       p.786), charged to the chaff pool inside the ear on the first filling
!       day, since the model applies the rate of that law and otherwise
!       discards its initial condition. Zero new constants.
!   3 = mode 2 plus the grain unmet sink as one more claimant in the same
!       pro-rata resolution as the synthesis limbs.
! Read before the driver file is closed. The released value is 4.
! ==================================================================== !
READ (1,*,IOSTAT=ios) grain_n_partition_basis
IF (ios /= 0) grain_n_partition_basis = 0
IF (grain_n_partition_basis < 0 .OR. grain_n_partition_basis > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 175: grain_n_partition_basis must be 0 (OFF), 1 (restore the' // &
    ' post-anthesis continuity of ndb_rate_prev that rec 158 declares and rec 154 breaks),' // &
    ' 2 (1 + FPB 2008 p.786 Ngtot(0) = 2.4 mg/ear charged to the chaff on the first fill day)' // &
    ' 3 (2 + the grain''s unmet sink in the synthesis pro-rata claim set)' // &
    ' or 4 (3, but the Ngtot(0) charge is paid against the grain''s OWN dry mass at the' // &
    ' MEASURED flowering grain [N] of 2.0 % DM -- Bertheloot, Wu, Cournede & Andrieu 2011' // &
    ' Ann Bot 108:1097-1109 Table 1(a) Grains, 2.4 mg N on 120 mg DM, DOI 10.1093/aob/mcr125' // &
    ' -- instead of as a lump on a grain of zero dry mass); SHARED.f90;' // &
    ' read ', grain_n_partition_basis
  ERROR STOP 'rec 175: invalid grain_n_partition_basis'
END IF
! ==================================================================== !
! Driver record 176 vegetative_n_demand_basis: ON WHICH QUANTITY DOES THE
! PRE-ANTHESIS NITROGEN DEMAND ACT? The modes are CUMULATIVE.
!   1 = definition correction: the record-155 limb-(b) numerator reads the
!       CIRCULATING pool (the mobile vegetative pool plus the shoot nitrate
!       store) and not the stem and chaff luxury ORGAN nitrogen that record 157
!       mode 3 adds. Bertheloot, Cournede & Andrieu (2011) p.1090 define that
!       quantity as circulating amino acids and nitrate. Today 94-100 % of the
!       pre-anthesis signal is organ luxury nitrogen, and the regulation factor
!       discards 58.3 % of the demand while being blind to nitrogen status.
!   2 = mode 1 plus the critical-curve STOCK deficit of Justes et al. (1994)
!       Eq. (11) p.401 in the pre-anthesis demand, since a pure flux demand
!       compared with a stock ceiling that never binds cannot recover arrears;
!       still capped by the unchanged Eq. (17) p.405 maximum-curve stock gap.
! Read before the driver file is closed. The released value is 0.
! ==================================================================== !
READ (1,*,IOSTAT=ios) vegetative_n_demand_basis
IF (ios /= 0) vegetative_n_demand_basis = 0
IF (vegetative_n_demand_basis < 0 .OR. vegetative_n_demand_basis > 2) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 176: vegetative_n_demand_basis must be 0 (OFF, the E_N numerator' // &
    ' keeps the rec-157 stem and chaff luxury organ N), 1 (the numerator is the CIRCULATING pool' // &
    ' only, Bertheloot 2011a p.1090) or 2 (1 + the Justes 1994 Eqn 11 p.401 critical-curve STOCK deficit in' // &
    ' the pre-anthesis demand); SHARED.f90; read ', vegetative_n_demand_basis
  ERROR STOP 'rec 176: invalid vegetative_n_demand_basis'
END IF
! ==================================================================== !
! Driver record 177 grain_n_potential_trait: THE GRAIN-NITROGEN POTENTIAL AS A
! CULTIVAR TRAIT.
!   0 = the released value: one compiled per-kernel ceiling, 1.38 mg N
!       kernel-1, for all eight cultivars, while the kernel-MASS potential is
!       already a per-cultivar trait (row 2 of input/cultivar_yield.csv).
!   1 = the optional, name-keyed row grainN_pot_mgN_kernel of
!       input/cultivar_yield.csv overrides the ceiling for that cultivar; a
!       cultivar whose file carries no such row keeps the compiled default. The
!       value is ENFORCED inside the CONTROL (intact-ear) range 0.75-1.11 mg N
!       kernel-1 of Martre, Porter, Jamieson & Triboi (2003) Plant Physiol
!       133:1959-1967 Table I p.1961 with a 10 % margin: outside
!       [0.675, 1.221] is an error stop, inside the margin but outside the
!       range runs and warns.
! One code path for every case: no cultivar name and no site name appears in it.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) grain_n_potential_trait
IF (ios /= 0) grain_n_potential_trait = 0
IF (grain_n_potential_trait < 0 .OR. grain_n_potential_trait > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 177: grain_n_potential_trait must be 0 (OFF, the compiled' // &
    ' per-kernel grain-N ceiling gnf_mgN_kernel_hi = 1.38 mg N kernel-1 applies to every' // &
    ' cultivar) or 1 (the optional name-keyed cultivar_yield.csv row grainN_pot_mgN_kernel' // &
    ' overrides it per cultivar, enforced inside the Martre 2003 Table I control range' // &
    ' 0.75-1.11 mg N kernel-1); SHARED.f90; read ', grain_n_potential_trait
  ERROR STOP 'rec 177: invalid grain_n_potential_trait'
END IF
! ==================================================================== !
! Driver record 178 chaff_structural_n_basis: THE CHAFF'S OWN MEASURED
! STRUCTURAL-NITROGEN FLOOR.
!   0 = the released value: the chaff keeps the sheath and internode floor
!       3.5e-3 g N g-1 DM (Bertheloot, Martre & Andrieu 2008 Plant Physiol
!       148:1707-1720 p.1711, doi:10.1104/pp.108.124156), which is a STEM
!       number.
!   1 = the chaff uses its OWN measured floor 5.10e-3 g N g-1 DM =
!       1.07 mg structural N / 210 mg structural DM (Bertheloot, Wu, Cournede &
!       Andrieu 2011 Ann Bot 108(6):1097-1109 Table 1 section (c),
!       "Photosynthetic entities", column "Chaff", doi:10.1093/aob/mcr125) at
!       the ten chaff places that read a floor and nowhere else.
!   2 = the remobilization-efficiency law: the chaff post-anthesis release is
!       bounded by its MEASURED remobilization efficiency, the chaff nitrogen
!       staying at or above (1 - 0.56) times its anthesis value, with the
!       measured chaff structural minimum 3.80e-3 g N g-1 DM as the hard floor.
!       The efficiency 0.56 is from Pask, Sylvester-Bradley, Jamieson & Foulkes
!       (2012) Field Crops Res 126:104-118 Table 4 p.110 column CH (three-site
!       mean 55.7 %, abstract "chaff (56%)"),
!       doi:10.1016/j.fcr.2011.09.021, corroborated at 0.563 and 0.565 by
!       Barraclough, Lopez-Bellido & Hawkesford (2014) Field Crops Res
!       156:242-248 p.247. The 3.80e-3 floor is Pask et al. (2012) Methods
!       2.9(i) pp.107-108 and Results p.111, chaff 0.38 %, the lowest harvest
!       chaff nitrogen percentage in any plot at zero nitrogen. The nitrogen
!       dependence of the chaff residual is EMERGENT: the floor is a constant
!       fraction of the anthesis pool the model already computes.
!   3 = the same law with the 5.10e-3 floor instead, so that the law and the
!       floor can be ruled on separately.
! One floor for every chaff on every deck; no cultivar name and no site name in
! any statement this record controls. Modes 1 and 2 are retained as an inactive
! alternative, see Part 1 Table 6.
! ==================================================================== !
READ (1,*,IOSTAT=ios) chaff_structural_n_basis
IF (ios /= 0) chaff_structural_n_basis = 0
IF (chaff_structural_n_basis < 0 .OR. chaff_structural_n_basis > 3) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 178: chaff_structural_n_basis must be 0 (OFF, the chaff keeps' // &
    ' the sheath/internode structural-N floor n2_struct_stem = 3.5e-3 g N g-1 DM) or 1 (the' // &
    ' chaff uses its OWN measured floor n2_struct_chaff = 5.10e-3 g N g-1 DM = 1.07 mg N /' // &
    ' 210 mg DM, Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109 Table 1(c)' // &
    ' column Chaff, DOI 10.1093/aob/mcr125) or 2 (the NRE LAW: N_chaff >= (1 - 0.56) x' // &
    ' N_chaff(anthesis), Pask et al. 2012 Field Crops Res 126:104 Table 4 p.110 column CH, with' // &
    ' Pask''s measured chaff structural minimum 3.80e-3 g N g-1 DM as the hard floor) or 3 (the' // &
    ' same law on Bertheloot 2011b''s 5.10e-3 floor); SHARED.f90; read ', &
    chaff_structural_n_basis
  ERROR STOP 'rec 178: invalid chaff_structural_n_basis'
END IF
! ==================================================================== !
! Driver record 179 lamina_expansion_n_limit: A NITROGEN LIMIT ON THE
! LAMINA-AREA EXPANSION INCREMENT.
!   1 = the released value: the minimum photosynthetic nitrogen per unit area
!       is a MIN(running maximum, b) + c, Bertheloot, Martre & Andrieu (2008)
!       Funct Plant Biol 35:781-796 Eq. (8) with Table 1 p.784 (a = 0.46,
!       b = 1.9, c = 0.24 g N m-2), doi:10.1071/FP08064, and the supportable
!       area is the lamina photosynthetic nitrogen over that minimum.
!   2 = the Eq. (8) plateau a b + c = 1.114 g N m-2, with no running maximum.
!   3 = the total lamina nitrogen at 1.5 g N m-2 of Jamieson & Semenov (2000)
!       Field Crops Res 68:21-29 p.23, where their own word is "assumed".
!   4 = attribution arm: the logarithmic form of Yin, Schapendonk, Kropff, van
!       Oijen & Bindraban (2004) NJAS 51(4):407-427 Appendix 2 p.424 with the
!       model's own canopy extinction coefficient 0.41 and base leaf nitrogen
!       0.40.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) lamina_expansion_n_limit
IF (ios /= 0) lamina_expansion_n_limit = 0
IF (lamina_expansion_n_limit < 0 .OR. lamina_expansion_n_limit > 4) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 179: lamina_expansion_n_limit must be 0 (OFF), 1' // &
    ' (PRIMARY, Bertheloot 2008 FPB Eqn 8 minN^ph,a with the running maximum), 2 (the' // &
    ' Eqn 8 plateau 1.114 g N m-2), 3 (Jamieson & Semenov 1.5 g N m-2 on total lamina N)' // &
    ' or 4 (the Yin l_N log form); SHARED.f90; read ', lamina_expansion_n_limit
  ERROR STOP 'rec 179: invalid lamina_expansion_n_limit'
END IF
! ==================================================================== !
! Driver record 180 grain_ca_loading_basis: THE GRAIN CALCIUM LOADING LAW.
!   0 = the record-166 calcium route decides grain calcium and nothing in this
!       record is read.
!   1 = the released value, the XYLEM-STREAM LOADING LAW. Calcium travels in
!       the transpiration stream at one concentration and is deposited
!       irreversibly where that water leaves the xylem; from anthesis to
!       maturity the measured share arriving at the head is divided at the
!       grain-base xylem discontinuity (Zee & O'Brien 1970 Aust J Biol Sci
!       23:783, pp.783-784) between grain and chaff:
!         dCa_grain/dt = f_grain      c_xyl Tr_act(t)
!         dCa_chaff/dt = (f_head - f_grain) c_xyl Tr_act(t)
!         remainder    -> the vegetative terminus, never remobilized.
!       f_grain = 0.30549 and f_head = 0.5827 are arithmetic on Hocking (1994)
!       J Plant Nutr 17:1289 Tables 1(b) and 1(c) p.1294, Table 2 and Table 3
!       p.1295, doi:10.1080/01904169409364807. The xylem concentration is the
!       EXISTING calcium-row column 19 of input/driver_minerals.csv, sourced
!       default 15.0 mg L-1, ENFORCED inside the measured xylem-sap envelope
!       12.0-661 mg L-1 (White & Broadley 2003 Ann Bot 92:487 p.490; Atkinson,
!       Ruiz & Mansfield 1992 J Exp Bot 43:1315 Table 3 p.1319) with a 10 %
!       margin. The record-166 ear-area share and the grain-water route are
!       retired for calcium under this mode, and the record-166-era constructed
!       concentrations 842.0246 and 654.339 mg L-1 are withdrawn: they lie
!       outside the measured envelope and cannot be run under mode 1.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) grain_ca_loading_basis
IF (ios /= 0) grain_ca_loading_basis = 0
IF (grain_ca_loading_basis < 0 .OR. grain_ca_loading_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 180: grain_ca_loading_basis must be 0 (OFF, the record-166' // &
    ' calcium route stands) or 1 (the xylem-stream loading law: Ca = c_xyl_Ca x Tr_act,' // &
    ' partitioned from anthesis by the Hocking 1994 Table 1(c) p.1294 measured fractions;' // &
    ' SHARED.f90); read ', grain_ca_loading_basis
  ERROR STOP 'rec 180: invalid grain_ca_loading_basis'
END IF
! ==================================================================== !
! Driver record 181 mg_solution_basis: MAGNESIUM JOINS THE SOLUTION BRANCH.
!   0 = the released value: magnesium keeps the linear partition, its solution
!       concentration being the pool over the distribution coefficient.
!   1 = column 17 of the magnesium row of input/driver_minerals.csv IS the
!       soil-solution magnesium concentration, read in exactly the place iron,
!       zinc, copper and manganese already read theirs, feeding the unchanged
!       mass-flow and Barber-Cushman diffusion terms. The source of the shipped
!       11.4 mg L-1 and of the enforced window is Barakat et al. (2025) Sci Rep
!       15:44670 Table 2 p.3, doi:10.1038/s41598-025-29041-y.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mg_solution_basis
IF (ios /= 0) mg_solution_basis = 0
IF (mg_solution_basis < 0 .OR. mg_solution_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 181: mg_solution_basis must be 0 (OFF, magnesium keeps the' // &
    ' linear Kd partition (pool_mg_kg)/Kd for its solution concentration) or 1 (magnesium' // &
    ' joins the record-131 solution branch: the Mg row avail_A_mgL_CAL of' // &
    ' input/driver_minerals.csv is the soil-solution Mg concentration, mg L-1);' // &
    ' SHARED.f90; read ', mg_solution_basis
  ERROR STOP 'rec 181: invalid mg_solution_basis'
END IF
! ==================================================================== !
! Driver record 182 is a RESERVED row. input/driver_main.txt is POSITIONAL, so
! the line is read as an integer and has no effect anywhere, which keeps the
! number of record 183 fixed. The process it is held for is an alternative
! vegetative-nitrogen demand basis, which is not built. No range check: the
! legal range of a reserved row belongs to whatever builds it. Read before the
! driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) reserved_rec_182
IF (ios /= 0) reserved_rec_182 = 0
! ==================================================================== !
! Driver record 183 mineral_buffer_power_basis: THE SOIL BUFFER POWER OF THE
! MINERAL DIFFUSION TERM.
!   0 = the Barber-Cushman effective diffusion coefficient divides by the
!       column-19 CONSTANT buffer power of input/driver_minerals.csv, which is
!       1.0 on every shipped element row and which record 140 declares to have
!       no source.
!   1 = the released value: the buffer power b_e = theta + rho_b Kd_e is
!       recomputed DAILY per element from the model's own soil water, the deck
!       own MEASURED topsoil bulk density (input/driver_N.csv position 70,
!       already read for records 139 and 168) and the element's own column-3
!       distribution coefficient - the linear-isotherm buffer power of
!       Baldwin, Nye & Tinker (1973) Plant and Soil 38:621-635 (p.622, b is the
!       buffer power of the soil, dC/dCl; p.625, Db = Dl W fl), decomposed as
!       in Roose & Fowler (2004) J Theor Biol 228:173-184 Eq. (2.8) p.175.
! Read before the driver file is closed.
! ==================================================================== !
READ (1,*,IOSTAT=ios) mineral_buffer_power_basis
IF (ios /= 0) mineral_buffer_power_basis = 0
IF (mineral_buffer_power_basis < 0 .OR. mineral_buffer_power_basis > 1) THEN
  WRITE (*,'(A,I0)') ' ERROR rec 183: mineral_buffer_power_basis must be 0 (OFF, the mineral' // &
    ' diffusion coefficient divides by the driver_minerals.csv column-19 constant' // &
    ' buffer_power_b) or 1 (b_e = theta + rho_b*Kd_e recomputed daily per element from the' // &
    ' model''s own water content, the deck''s measured bulk density and the element''s own' // &
    ' column-3 Kd); SHARED.f90; read ', mineral_buffer_power_basis
  ERROR STOP 'rec 183: invalid mineral_buffer_power_basis'
END IF
CLOSE (1)
! Record 179 dependencies, refused rather than left a silent no-op: without
! record 84 there is no organ-level green lamina area to constrain, and without
! records 121 and 122 there is no increment ledger to commit at a fraction.
IF (lamina_expansion_n_limit >= 1 .AND. &
    (lai_from_organs /= 1 .OR. n_payment_gate < 1 .OR. lai_n_arrest < 1)) THEN
  WRITE (*,'(A,I0,A,I0,A,I0)') ' ERROR rec 179: lamina_expansion_n_limit >= 1 requires' // &
    ' rec 84 lai_from_organs = 1, rec 121 n_payment_gate >= 1 and rec 122 lai_n_arrest >= 1;' // &
    ' read rec 84 = ', lai_from_organs, ', rec 121 = ', n_payment_gate, &
    ', rec 122 = ', lai_n_arrest
  ERROR STOP 'rec 179 requires rec 84 = 1 and rec 121 >= 1 and rec 122 >= 1'
END IF
! rec 183: DEPENDENCY REFUSED, NEVER A SILENT NO-OP. There is no
! diffusion term to divide unless record 140 mineral_supply_mech >= 1.
IF (mineral_buffer_power_basis >= 1 .AND. mineral_supply_mech < 1) THEN
  WRITE (*,'(A,I0,A)') ' ERROR rec 183 = 1: the buffer power b_e only enters the record-140' // &
    ' diffusion-to-the-root-surface term De_e = Dl_e*theta*f_l/b_e, and there is no such term' // &
    ' unless rec 140 mineral_supply_mech >= 1; read rec 140 = ', mineral_supply_mech, &
    '. Set rec 140 >= 1 or rec 183 = 0.'
  ERROR STOP 'rec 183 requires rec 140 mineral_supply_mech >= 1'
END IF
! rec 173: DEPENDENCIES REFUSED, NEVER A SILENT NO-OP.
IF (spac_sink_coupling >= 1) THEN
  IF (spac_mode < 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 173: spac_sink_coupling = 1 couples the SPAC leaf''s own' // &
      ' transpiration to the Richards sink, and there is no SPAC leaf unless rec 81' // &
      ' spac_mode >= 1; read rec 81 = ', spac_mode, '. Set rec 81 >= 1 or rec 173 = 0.'
    ERROR STOP 'rec 173 requires rec 81 spac_mode >= 1'
  END IF
  IF (ft_pot_unstressed /= 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 173: the statement this record removes (ebh_etr_dem_mm =' // &
      ' ebh_etr_dem_unstr_mm, and the trans_act / trans_act_pot override) exists ONLY on the' // &
      ' rec 78 = 1 branch; read rec 78 = ', ft_pot_unstressed, '. Set rec 78 = 1 or rec 173 = 0.'
    ERROR STOP 'rec 173 requires rec 78 ft_pot_unstressed = 1'
  END IF
END IF
! Record 166 dependency, refused rather than left a silent zero. The ear-area
! state the route integrates is built only by record 68, which itself requires
! records 64 to 67. Without it the ear share would be identically zero and the
! grain would receive only the net-import term.
IF (grain_ca_xylem_throughput == 1 .AND. canopy_sheath_ear /= 1) THEN
  WRITE (*,'(A)') ' ERROR rec 166: grain_ca_xylem_throughput = 1 integrates the EAR GREEN AREA' // &
    ' (cse_eai_meas against cls_gai + cse_sai + cse_iai), which exists only when rec 68' // &
    ' canopy_sheath_ear = 1 (itself requiring recs 64-67). Set rec 68 = 1 or rec 166 = 0.'
  ERROR STOP 'rec 166 requires rec 68 canopy_sheath_ear = 1'
END IF
! Record 167 refusals. The statements this record rebases exist only inside the
! record-154 limb, and mode 1 reads a state variable only record 68 builds. A
! flag that cannot act does not ship as a silent no-op.
IF (lamina_synthesis_sheath_basis >= 1) THEN
  IF (organ_n_turnover_nema < 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 167: the lamina/stem synthesis statements this record rebases' // &
      ' exist only inside the rec-154 Bertheloot 2011a organ-N turnover limb; read rec 154 = ', &
      organ_n_turnover_nema, '. A flag that cannot act is not shipped as a silent no-op.'
    ERROR STOP 'rec 167 requires rec 154 >= 1'
  END IF
  IF (lamina_synthesis_sheath_basis == 1 .AND. canopy_sheath_ear /= 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 167 mode 1: the area-derived sheath mass reads cse_sai, which is' // &
      ' built only when rec 68 canopy_sheath_ear = 1; read rec 68 = ', canopy_sheath_ear, &
      '. Use rec 167 = 2 (the printed Bertheloot 2011b Table 1c ratio) or turn rec 68 on.'
    ERROR STOP 'rec 167 mode 1 requires rec 68 = 1'
  END IF
END IF
! Record 174 refusals. Mode 1 needs the rank chronology (record 64) and the
! blade-area step (record 65) to exist at all, and it needs the canopy clock to
! BE the gated development clock, because the validity of the prediction rests
! on the identity between the canopy clock reading and the development units
! since emergence, which holds only at record 162 = 1 with record 161 = 0. With
! record 161 at or above 1 the canopy clock is a SOIL-temperature clock and the
! booting threshold is no longer the clock reading at booting, so the
! prediction would answer a different question.
IF (canopy_fln_prediction >= 1) THEN
  IF (canopy_leaf_ranks /= 1 .OR. canopy_blade_area /= 1) THEN
    WRITE (*,'(A,I0,A,I0,A)') ' ERROR rec 174: the predicted final leaf number is the index of the' // &
      ' blade-dimension profiles, which exist only when rec 64 canopy_leaf_ranks = 1 and rec 65' // &
      ' canopy_blade_area = 1; read rec 64 = ', canopy_leaf_ranks, ', rec 65 = ', canopy_blade_area, &
      '. A flag that cannot act is not shipped as a silent no-op.'
    ERROR STOP 'rec 174 requires rec 64 = 1 and rec 65 = 1'
  END IF
  IF (canopy_clock_basis /= 1 .OR. canopy_clock_driver /= 0) THEN
    WRITE (*,'(A,I0,A,I0,A)') ' ERROR rec 174: the prediction FLN_pred = FLOOR(HS(stage_thresh(4) -' // &
      ' tdu_at_emerge)) rests on the identity clr_tau_Cd == cum_tdu - tdu_at_emerge, which holds only' // &
      ' when the canopy clock IS the gated development increment: rec 162 = 1 AND rec 161 = 0. Read' // &
      ' rec 162 = ', canopy_clock_basis, ', rec 161 = ', canopy_clock_driver, &
      '. SHARED.f90: a prediction on a wrong premise is refused, not approximated.'
    ERROR STOP 'rec 174 requires rec 162 = 1 and rec 161 = 0'
  END IF
END IF
! Record 175 refusals. Mode 1 repairs a continuity that exists only inside the
! record-158 substitution and is broken only by the record-154 limb; modes 2
! and 3 edit statements that live only inside that limb.
IF (grain_n_partition_basis >= 1) THEN
  IF (organ_n_turnover_nema < 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 175: every mode of this record edits statements that exist' // &
      ' only inside the rec-154 Bertheloot 2011a organ-N turnover limb (mode 1 restores the continuity term' // &
      ' that limb zeroes, modes 2-3 edit its grain statements); read rec 154 = ', &
      organ_n_turnover_nema, '. A flag that cannot act is not shipped as a silent no-op.'
    ERROR STOP 'rec 175 requires rec 154 >= 1'
  END IF
  IF (nitrate_drain_rate_basis < 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 175: mode 1 repairs the numerator ndb_rate_prev, which is' // &
      ' READ only when rec 158 nitrate_drain_rate_basis >= 1 (at rec 158 = 0 the Busch N_max' // &
      ' supply term is the shoot nitrate STORE and the continuity term is never read); read' // &
      ' rec 158 = ', nitrate_drain_rate_basis, '. Set rec 158 = 1 or rec 175 = 0.'
    ERROR STOP 'rec 175 requires rec 158 >= 1'
  END IF
END IF
! Record 176 refusals. Mode 1 narrows a numerator that is FORMED by record 157
! mode 3 and READ only by record 155 limb (b); mode 2 enlarges the pre-anthesis
! demand that limb then multiplies. Without either, the record is a no-op.
IF (vegetative_n_demand_basis >= 1) THEN
  IF (.NOT. (n_uptake_regulation == 1 .OR. n_uptake_regulation == 3)) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 176: the E_N numerator this record redefines is READ only by' // &
      ' rec 155 n_uptake_regulation limb (b), i.e. rec 155 = 1 or 3; read rec 155 = ', &
      n_uptake_regulation, '. A flag that cannot act is not shipped as a silent no-op.'
    ERROR STOP 'rec 176 requires rec 155 = 1 or 3'
  END IF
  IF (lamina_n_storage /= 3) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 176 mode 1: the stem and chaff luxury terms it drops are ADDED' // &
      ' to the numerator only by rec 157 lamina_n_storage = 3; read rec 157 = ', lamina_n_storage, &
      '. At rec 157 /= 3 the numerator is already the circulating pool and the record is a no-op.'
    ERROR STOP 'rec 176 requires rec 157 = 3'
  END IF
  IF (mobile_n_nema_estimate >= 1) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 176 mode 1: rec 159 mobile_n_nema_estimate REPLACES the whole' // &
      ' E_N numerator downstream of this record, so the two are alternative resolutions of the' // &
      ' same quantity; read rec 159 = ', mobile_n_nema_estimate, '. Set rec 159 = 0 or rec 176 = 0.'
    ERROR STOP 'rec 176 is incompatible with rec 159 >= 1'
  END IF
END IF
! Record 154 refusals. This record replaces statements that exist ONLY inside
! the record-86 Eq. (6) branch, and three other records are alternative
! resolutions of the very claim it replaces. A run that asks for both is
! refused.
IF (organ_n_turnover_nema >= 1) THEN
  IF (.NOT. (grainN_sink_mech == 1 .OR. grainN_sink_mech == 2 .OR. grainN_sink_mech == 4)) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 154: the Bertheloot 2011a organ-N turnover replaces the lamina Eqn-6' // &
      ' synthesis and the pooled grain+lamina claim, which exist only under rec 86' // &
      ' grainN_sink_mech = 1, 2 or 4; read rec 86 = ', grainN_sink_mech, '.'
    ERROR STOP 'rec 154 requires rec 86 = 1, 2 or 4'
  END IF
  IF (grainN_partition == 1) THEN
    WRITE (*,'(A)') ' ERROR rec 154: rec 63 grainN_partition = 1 calls apply_organN_turnover_daily,' // &
      ' which would apply the SAME Bertheloot 2011a turnover a second time on the same day. Mutually exclusive.'
    ERROR STOP 'rec 154 is mutually exclusive with rec 63'
  END IF
  IF (leafn_hold_law >= 1 .OR. leafn_turnover_law >= 1) THEN
    WRITE (*,'(A,I0,A,I0)') ' ERROR rec 154: records 132 (leafn_hold_law) and 133 (leafn_turnover_law)' // &
      ' are alternative resolutions of the lamina/grain claim this record replaces; read rec 132 = ', &
      leafn_hold_law, ', rec 133 = ', leafn_turnover_law
    ERROR STOP 'rec 154 is mutually exclusive with rec 132 and rec 133'
  END IF
  WRITE (*,'(A,I0,A)') ' rec 154 organ_n_turnover_nema = ', organ_n_turnover_nema, &
    ': Bertheloot 2011a post-anthesis organ N turnover ACTIVE (delta 0.008 degCd-1; sigma lam 1.5e-4,' // &
    ' stem/chaff 3.0e-5 g g-1 degCd-1; root limb ' // &
    MERGE('ON ', 'OFF', organ_n_turnover_nema == 1) // ')'
END IF
! Record 155 banner. Nothing is printed at mode 0. The banner names which limbs
! are live, including the record-154 condition on limb (c): a mode 1 or 3 run
! with record 154 = 0 is not refused, it runs limb (b) only and says so.
IF (n_uptake_regulation >= 1) THEN
  WRITE (*,'(A,I0,A,A,A,A,A,A)') ' rec 155 n_uptake_regulation = ', n_uptake_regulation, &
    ': root N uptake DEMAND regulation ACTIVE -- (a) Justes 1994 maximum-dilution marginal' // &
    ' demand ', MERGE('ON ', 'OFF', n_uptake_regulation == 1 .OR. n_uptake_regulation == 2), &
    '; (b) Bertheloot 2011a Eqn 4 E_N = exp(-130 [N]mob) ', &
    MERGE('ON ', 'OFF', n_uptake_regulation == 1 .OR. n_uptake_regulation == 3), &
    '; (c) grain demand net of the mobile pool ', &
    MERGE('ON ', 'OFF', (n_uptake_regulation == 1 .OR. n_uptake_regulation == 3) .AND. &
                        organ_n_turnover_nema >= 1)
  IF ((n_uptake_regulation == 1 .OR. n_uptake_regulation == 3) .AND. organ_n_turnover_nema == 0) &
    WRITE (*,'(A)') ' rec 155 NOTE: limb (c) is NOT applied because rec 154 = 0 holds no' // &
      ' culm-common mobile pool (SHARED.f90 rec 155 block (c)); limb (b) still applies.'
END IF
! Record 117 dependency, refused rather than defaulted. The on path must RETIRE
! the anchored storage-peak construction of input/cultivar_yield.csv rather
! than sit beside it. The model already holds the non-anchored replacement:
! under record 89 = 1 or 2 the spike partition replaces both nodes with the
! measured 0.44 of Fischer (1985) folio 452 on his own developmental onset, and
! the anchored node is then not read at all. Records 89 = 0 and 3 leave it live
! and are refused.
IF (chaff_growth_mode >= 1) THEN
  IF (.NOT. (grain_number_mech == 1 .OR. grain_number_mech == 2)) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 117: the ON path must retire the ANCHORED storage_peak_frac' // &
      ' construction, which only rec 89 grain_number_mech = 1 or 2 does (ITEM D replaces both' // &
      ' spike-partition nodes with Fischer 1985 p.452); read rec 89 = ', grain_number_mech, &
      '. A chaff sink beside an observation-anchored anthesis target is not shipped.'
    ERROR STOP 'rec 117 requires rec 89 = 1 or 2'
  END IF
END IF
! Record 115 dependencies, refused rather than defaulted. The per-kernel
! divisor is written only inside the grain-number routine, which is called only
! under record 89 = 1 or 3. Under record 89 = 0 or 2 the model holds no
! per-cohort kernel state at all and no divisor may be invented. Record 86 must
! also reach the line being changed, or the flag would be a silent no-op.
IF (gnf_ng0_basis == 1) THEN
  IF (.NOT. (grain_number_mech == 1 .OR. grain_number_mech == 3)) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 115: the per-kernel basis needs the per-cohort kernel state' // &
      ' gnx_kn, which is written ONLY under rec 89 grain_number_mech = 1 or 3; read rec 89 = ', &
      grain_number_mech, '. No divisor is invented where the model holds no kernel state.'
    ERROR STOP 'rec 115 requires rec 89 = 1 or 3'
  END IF
  IF (.NOT. (grainN_sink_mech >= 1 .AND. grainN_sink_mech /= 3)) THEN
    WRITE (*,'(A,I0,A)') ' ERROR rec 115: the basis it changes lives inside the rec 86 item-F grain-N' // &
      ' sink, which is reached only when rec 86 grainN_sink_mech >= 1 and /= 3; read rec 86 = ', &
      grainN_sink_mech, '. A flag that cannot act is not shipped as a silent no-op.'
    ERROR STOP 'rec 115 requires rec 86 >= 1 and /= 3'
  END IF
END IF
IF (sowing_model_day > 0) THEN
  iD_sowing = sowing_model_day
  WRITE (*,'(A,I0,A)') ' rec 88 sowing_model_day = ', iD_sowing, &
    ' : DAS 0 = this model day (fertiliser, residue event, residue-evap mass).'
ELSE
  WRITE (*,'(A,I0,A)') ' rec 88 sowing_model_day = 0 (not declared): legacy anchor iD_sowing = ', &
    iD_sowing, ' (WheatDryFACE).'
END IF
IF (leafN_cap_route == 1) WRITE (*,'(A)') &
  ' rec 87 leafN_cap_route ON: leaf N bounded by sln_max_N * LAI(assim); overflow to N_mobile_veg.'
IF (itemF_Nclose_mode >= 1) WRITE (*,'(A,I0,A)') &
  ' rec 97 itemF_Nclose_mode = ', itemF_Nclose_mode, &
  ' : 1/3 = same-day Eqn-6 synthesis respects the leaf-N capacity; 2/3 = grain-N uptake demand' // &
  ' = FPB Eqn 12 potential RATE net of the mobile pool, + terminal close of N_mobile_veg into N_stem.'
IF (n_payment_gate >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 121 n_payment_gate = ', n_payment_gate, &
    ' : organ structural growth is PAID in reduced N before it is built (1 = floors, 2 = construction concentrations);' // &
    ' unpaid carbon -> stem WSC up to wsc_cap_frac, remainder = explicit residual (never structure). SHARED [NPG-0..6].'
  WRITE (*,'(A)') ' rec 121 UNITS: tier concentrations g N g-1 DM of the increment; pool g N m-2 d-1;' // &
    ' unpaid carbon g CH2O m-2 d-1; WSC deposit g DM m-2 d-1'
  WRITE (*,'(A,F8.5,A,F8.5,A,F8.5,A)') '   mode 1 floors        c_L = ', npg_floor_L, ' (Bertheloot 2008a 0.79 % DM)  c_S = ', &
    n5_stem_nstruct, ' (Jamieson & Semenov 2000 0.5 %)  c_R = ', npg_floor_R, ' (DECLARED)'
  WRITE (*,'(A,F8.5,A,F8.5,A,F8.5,A)') '   mode 2 construction  c_L = ', n4_sln_target / csl_LMA, &
    ' (= n4_sln_target 1.5 g N m-2 / csl_LMA 48.5 g m-2)  c_S = ', n4_stem_nmax, ' (Jamieson & Semenov 2000 1.5 %)  c_R = ', &
    npg_build_R, ' (Hicks 1928 1.29 %)'
  WRITE (*,'(A,F6.3,A)') '   stem luxury release  k_lux = ', npg_k_lux, &
    ' d-1 (DECLARED); order: root-reduced N -> root tier first, surplus exported; shoot pool -> lamina -> stem/chaff' // &
    ' (DECLARED) -> root remainder (phloem return, DECLARED); stem luxury N last (option b)'
END IF
IF (lai_n_arrest >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 122 lai_n_arrest = ', lai_n_arrest, &
    ' : the day''s lamina-area increment is committed at the rec-121 paid fraction p_L (unpaid area removed the same' // &
    ' evening, never used for photosynthesis, no carbon charged); 2 also sheds green area below SLN_min_N (0.50 g N m-2),' // &
    ' oldest rank first, through the rec-84 leaf-death coupling (PRECEDENT rule). SHARED [LNA-1..4]. No new constant.'
END IF
IF (tiller_n_afford >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 123 tiller_n_afford = ', tiller_n_afford, &
    ' : the rec 99 = 4 cohort ledger also carries lamina construction N (demand at the rec-121 c_L, supply = last' // &
    ' gate''s p_L x demand); an axis class survives only if carbon AND nitrogen shares are paid, phi = MIN(phi_C, phi_N).' // &
    ' SHARED [TNA-1..3]. No new constant; byte-identical to rec 99 = 4 when rec 121 = 0.'
END IF
IF (root_front_law >= 1) THEN
  WRITE (*,'(A,I0,A,F5.1,A)') ' rec 124 root_front_law = ', root_front_law, &
    ' : root front descends linearly from emergence at ', rfl_rz_mode(root_front_law), &
    ' mm d-1 (Kirkegaard & Lilley 2007 abstract range 7.9-11.3; 1 lower bound, 2 midpoint DECLARED primary,' // &
    ' 3 upper bound) from the 0.10 m seedling floor to root_max_depth; replaces the Borg & Grimes 1986' // &
    ' season envelope. SHARED [RFL-1..3].'
END IF
IF (kernel_potential_law >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 125 kernel_potential_law = ', kernel_potential_law, &
    ' : 1 = per-kernel grain-N ceiling 1.38 mg N kernel-1 (Martre 2003 Table I); 2 = + kernel-mass potential' // &
    ' from the cell-division-window supply within the cultivar envelope (cultivar_yield.csv rows 7-8);' // &
    ' 3 = potential only. SHARED [KPT-1..4].'
END IF
IF (mineral_grain_sink_cap == 1) THEN
  WRITE (*,'(A)') ' rec 126 mineral_grain_sink_cap = 1 : daily grain loading of each element bounded by grain_max_mg_kg' // &
    ' (driver_minerals.csv column 16; 0 = no bound) x the grain DM increment; refused share stays in leaf / dead tissue.' // &
    ' SHARED [MGS-1..3].'
END IF
IF (repro_surplus_fate == 1) THEN
  WRITE (*,'(A)') ' rec 127 repro_surplus_fate = 1 : the grain-sink-full reproductive surplus that rec 92 does not label as' // &
    ' stem WSC (deposition window closed or capacity full) is booked to the explicit residual, never to stem structure;' // &
    ' the labeled part stays in W_stem / W_reserve as before. SHARED [RSF-1..3].'
END IF
IF (n_construction_status >= 1) THEN
  ! The record-156 banner names the curve it is actually using: the mode-1 arm
  ! says maximum-nitrogen curve (a = 8.3, b = 0.44).
  IF (n_construction_ceiling_curve >= 1) THEN
    WRITE (*,'(A,I0,A)') ' rec 128 n_construction_status = ', n_construction_status, &
      ' : rec 121''s single-tier payment replaced by a two-pass' // &
      ' floor/target law -- pass 1 pays the organ structural floor (luxury allowed, fixes p_L/p_S/p_R); pass 2 pays,' // &
      ' from what pass 1 left in the shoot pool (no luxury), the gap up to a target capped by the Justes/Brisson shoot' // &
      ' MAXIMUM-N curve (a=8.3, b=0.44, W = above-ground DM t/ha, floored at 1.55) [rec 156 = 1; the critical pair' // &
      ' 5.35/0.442 is NOT in use]. Mode 1 lamina floor = Bertheloot' // &
      ' 2008 structural N; mode 2 = SLN_min_N/csl_LMA; mode 3 = 0.8 g N m-2 after Martre et al. (2006)' // &
      ' = ncs_floor3_SLN_martre/csl_LMA; modes 2 and 3 are not used in v1.0. SHARED.'
  END IF
END IF
! Justes et al. (1994) Ann Bot 74:397-407 was read at first hand. It prints
! Eq. (17) p.405, the maximum nitrogen dilution curve Nmax = 8.3 DM^-0.44
! (Fig. 5 p.402, valid from 1.5 to 14 t ha-1, p.404 col. 2), so the model's own
! 8.3 IS the value of that paper and the 8.47 of the STICS table is that model
! own number and not a transcription. It prints the critical pair as Eq. (11)
! p.401 (Fig. 4 p.401) with the 4.4 % plateau at p.402 col. 2, and - p.405, the
! primary source for the marginal law this file uses - the maximum rate of
! nitrogen accumulation as k times the maximum growth rate with
! k = 46.5 DM^-0.44, of which Eq. (47) of Brisson et al. (1998) is the STICS
! implementation. The paper states NOTHING about the maximum concentration
! below 1.5 t ha-1, so the 6.6 % plateau of mode 0 has no source; record 160
! mode 1 is the continuity construction for it.
! Record 156 banner: nothing is printed at mode 0, and the inert case
! (record 128 = 0) is reported rather than silently passed.
IF (n_construction_ceiling_curve >= 1) THEN
  IF (n_construction_status >= 1) THEN
    WRITE (*,'(A)') ' rec 156 n_construction_ceiling_curve = 1: the rec-128 pass-2 target ceiling reads the' // &
      ' Justes/Brisson MAXIMUM dilution pair (a = 8.3, b = 0.44) instead of the critical pair (5.35, 0.442);' // &
      ' same marginal form a(1-b)W^-b, same 1.55 t ha-1 validity floor, zero new constants. Source: Brisson et al.' // &
      ' 1998 Agronomie 18:311-346 p.342 Eqn (47) + p.329, Table VIII p.330 (pair attributed to Justes et al. 1994).' // &
      ' CAVEAT CLOSED 2026-10-16 (record 161/162 batch; the record-156 banner text itself is the only' // &
      ' thing that changed and it moves run.log only): Justes 1994 verified first hand 2026-10-16' // &
      ' (Justes 1994 Ann Bot 74:397; Eqn 17 p.405, k = 46.5 DM^-0.44);' // &
      ' STICS ADILMAX 8.47 is STICS''s own value.  The 6.6 % plateau item alone still stands and is carried' // &
      ' by record 160; it is NEVER READ here.'
  ELSE
    WRITE (*,'(A)') ' rec 156 n_construction_ceiling_curve = 1 but rec 128 n_construction_status = 0: the' // &
      ' ceiling block is not entered, so this record is INERT on this deck (reported, not refused).'
  END IF
END IF
! rec 157 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (lamina_n_storage >= 1) THEN
  WRITE (*,'(A,I0,A,F5.3,A,F6.3,A)') ' rec 157 lamina_n_storage = ', lamina_n_storage, &
    ' : pre-anthesis N above the organ construction demand is STORED in the green lamina up to the lamina' // &
    ' HOLDING CAPACITY sln_max_N_lnh = ', sln_max_N_lnh, ' g N m-2 of leaf area (= ', &
    100.0_wp * sln_max_N_lnh / csl_LMA, &
    ' %DM at this model''s csl_LMA 48.5), replacing the n4_sln_target = 1.5 g N m-2 construction TARGET as the' // &
    ' ladder ceiling only -- the target itself is UNCHANGED. Source: Pask (2009) Nottingham PhD p.29 printed' // &
    ' sentence, the SAME PARAMETER records 132/133 mode 2 use; ZERO new constants. Stored lamina N is' // &
    ' remobilizable (it is inside N_leaf, which rec 154 Bertheloot 2011a turnover degrades; the Bertheloot 2008 p.1711' // &
    ' structural floor stays) and feeds SLN -> Vcmax/Jmax through the existing assimilation path unchanged.'
  WRITE (*,'(A)') ' rec 157 mode 1: payment order = organ construction targets (rec 128) -> lamina storage' // &
    ' to capacity -> stem/chaff structural + root -> stem/chaff luxury (n4_stem_nmax) -> overflow to' // &
    ' N_mobile_veg / N_NO3_shoot, i.e. the INCUMBENT order; only the ceiling changed.'
  IF (lamina_n_storage == 3) THEN
    WRITE (*,'(A)') ' rec 157 mode 3: the rec-155 limb-(b) [N]mob NUMERATOR is widened from N_mobile_veg' // &
      ' alone to N_mobile_veg + N_NO3_shoot + MAX(0, N_stem - 0.005*W_stem) + MAX(0, N_chaff - 0.005*W_chaff);' // &
      ' denominator W_green + W_stem UNCHANGED. The record-157 LAMINA STORE is EXCLUDED, on Bertheloot 2011a''s own' // &
      ' sentence that a photosynthetic entity''s remobilizable N is "only N associated with the photosynthetic' // &
      ' apparatus" while the E_N numerator is the MOBILE common pool, "amino acids and nitrate" (Bertheloot' // &
      ' et al. 2011 Ann Bot 108:1085 eqn 4 and the Plant-structure paragraph). Stem/chaff luxury = Jamieson &' // &
      ' Semenov 2000 Field Crops Res 68:21-29 p.23 "may have a concentration of up to 1.5% N to store excess N".' // &
      ' ZERO new constants; Bertheloot 2011b''s own 10%-of-remobilizable estimate is NOT applied, so this is an UPPER' // &
      ' bound on [N]mob. Falsification target: Bertheloot 2011b prints E_N 0.20-0.27 on its high-N treatments.'
    IF (n_uptake_regulation /= 1 .AND. n_uptake_regulation /= 3) THEN
      WRITE (*,'(A,I0,A)') ' rec 157 mode 3 but rec 155 n_uptake_regulation = ', n_uptake_regulation, &
        ': limb (b) is not evaluated, so the widened numerator is INERT on this deck (reported, not refused).'
    END IF
  END IF
END IF
! rec 158 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (nitrate_drain_rate_basis >= 1) THEN
  WRITE (*,'(A,I0,A,F6.3,A)') ' rec 158 nitrate_drain_rate_basis = ', nitrate_drain_rate_basis, &
    ' : the Busch N_max SUPPLY term is no longer the STANDING vacuolar shoot nitrate store N_NO3_shoot but the' // &
    ' model''s own daily de novo nitrate-reduction DEMAND in the shoot (the previous day''s N_alloc_leaf_day +' // &
    ' N_alloc_stem_day + N_alloc_storage_day + MAX(0, N_from_assim_day), g N m-2 d-1), through the SAME unit chain' // &
    ' (/14, x1e6, /photoperiod, /LAI) and under the SAME unchanged Busch enzymatic ceiling N_max_ref = ', &
    N_max_ref, ' umol m-2 s-1.  Root N is excluded (root nitrate is reduced on the photorespiration-independent' // &
    ' Andrews limb, not in the chloroplast).  Source, Busch et al. (2018) Nature Plants 4:46-54 p.49, verbatim:' // &
    ' "will be limited by the maximum rate of de novo nitrogen supply to the chloroplast (Nmax). Nmax can be' // &
    ' limited by the nitrogen availability, the demand for nitrogen in the form of glycine and serine or by the' // &
    ' activities of the enzymes involved in the reduction of NO3- to NH4+."  ZERO new constants.'
  WRITE (*,'(A)') ' rec 158 DISCLOSED: a construction-N demand is a GROWTH flux (order 0.05 g N m-2 d-1),' // &
    ' one to two orders of magnitude below both the incumbent store-derived N_max and Busch''s own FITTED' // &
    ' N_max = 1.21 (Fig. 3a, an ENZYME-CAPACITY number).  alpha_G and alpha_S therefore fall toward 0 on most' // &
    ' days and A RISES on EVERY deck: this is a BASE-STATE move, gated on the 42-deck scorecard, not an N-axis' // &
    ' repair adopted on the NL sign alone.  The per-day size is in the five ON-only columns of' // &
    ' output/daily_nitrogen.csv and in the rec 158 season-summary line.'
END IF
! Record 159 banner: nothing is printed at mode 0.
! Record 160 banner: nothing is printed at mode 0.
IF (justes_nmax_floor_basis >= 1) THEN
  jnf_cap_pct = 8.3_wp * 1.55_wp**(-0.44_wp)
  WRITE (*,'(A,F6.3,A)') ' rec 160 justes_nmax_floor_basis = 1 : below the validity floor the Nmax OUTER' // &
    ' CAP on daily crop N demand is no longer the undocumented 6.6 % plateau but the printed curve''s own' // &
    ' value at the floor, 8.3 * 1.55**(-0.44) = ', jnf_cap_pct, ' %DM.  Justes, Mary, Meynard, Machet &' // &
    ' Thelier-Huche (1994) Ann Bot 74:397-407, DOI 10.1006/anbo.1994.1133, read at first hand: Eqn (17)' // &
    ' p.405 Nmax = 8.3 DM^-0.44, "valid from 1.5 to 14 t ha-1 (Fig. 5)" (p.404 col.2 - p.405 col.1).' // &
    '  THE PAPER PRINTS NO Nmax STATEMENT BELOW 1.5 t ha-1: the 6.6 % the mode-0 arm reads is an' // &
    ' undocumented constant whose in-code attribution to "p.399-400" is false.  ZERO new constants' // &
    ' (8.3, 0.44 and 1.55 are already in this file).'
  WRITE (*,'(A)') ' rec 160 DECLARED AS A MODEL CONSTRUCTION, NOT A QUOTATION: holding a curve at its' // &
    ' floor value below the floor is the SAME construction Justes applies to the CRITICAL curve (p.402' // &
    ' col.2, the 4.4 % horizontal line meeting Eqn (11) at DM = 1.55 t ha-1) and the SAME one records' // &
    ' 128/156 (growth_module.f90 MAX(W, 1.55)) and record 155 limb (a) already use; for the MAXIMUM curve' // &
    ' the paper states no level to quote, so the level is set by continuity.  The cap can only RISE' // &
    ' (+3.70 %, 6.6 -> 6.844), only on days with shoot DM <= 1.55 t ha-1, and only as the OUTER bound.'
END IF
! rec 161 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (canopy_clock_driver >= 1) THEN
  WRITE (*,'(A)') ' rec 161 canopy_clock_driver = 1 : the CANOPY clock (clr_tau_Cd, and with it leaf' // &
    ' appearance clr_HS, the final-leaf-number freeze, the blade-expansion ramp and the tiller-birth' // &
    ' windows) is advanced from emergence until the model''s own TERMINAL-SPIKELET marker by the' // &
    ' degree-days of the TOP SOIL LAYER, Tsoil_prof(1), and by screen air thereafter.  The degree-day' // &
    ' law, the base temperature and the cultivar cardinal rows are UNCHANGED (the same calc_dtt).' // &
    '  Source read at first hand: Abichou M (2016) these Univ. Paris-Saclay / AgroParisTech, NNT' // &
    ' 2016SACLA022, HAL tel-02316154, printed p.23: "soil temperature at 3 cm depth was monitored in' // &
    ' the experimental plots ... The apex temperature was estimated based on soil temperature before' // &
    ' stem extension and then based on air temperature"; justification printed p.38.  Corroborated by' // &
    ' Evers et al. (2006) Ann Bot 97:650 (soil 6 cm while the apex is below the surface), Gauthier et' // &
    ' al. (2020) JXB 71:5454 SI Table S1.6 p.28, and Jamieson et al. (1995) Field Crops Res 41:35-44' // &
    ' (the primary source all three attribute the rule to; NOT on disk, quoted SECOND HAND).'
  WRITE (*,'(A)') ' rec 161 DISCLOSURES: (i) the soil driver is a CHAPTER 1 + CHAPTER 2' // &
    ' protocol in Abichou, NOT thesis-wide -- her Chapter 4 drops it; (ii) the DEPTH IS A MEASUREMENT' // &
    ' CHOICE IN THE SOURCES, NOT A CONSTANT THIS MODEL FITS (Abichou 3 cm, Evers 6 cm, STICS sowing' // &
    ' depth): layer 1 is taken AS PRINTED, with no interpolation, weighting or smoothing.  MEASURED ON' // &
    ' THE 42 DECKS AND DISCLOSED: each deck prints its own soil grid, so layer 1''s bottom is -0.0333 m' // &
    ' at AGMIP_NL (where the measured 3 cm falls inside it), -0.1000 m on the 16 Braunschweig decks (a' // &
    ' 0-10 cm bulk), -0.0500 m at Maricopa and -0.0200 m on AGMIP_AR / WDF25 / WDF23 (where 3 cm is' // &
    ' BELOW layer 1).  What is implemented is THE MODEL''S TOP SOIL LAYER, whose effective depth each' // &
    ' deck''s own grid sets; selecting "the layer containing 3 cm" instead would make 3 cm a compiled' // &
    ' depth constant, i.e. the fitted quantity the sources say it is not; (iii) the' // &
    ' model declares no stage named "stem extension", so its nearest own declared marker, TERMINAL' // &
    ' SPIKELET (TDU = 400 degCd after emergence, McMaster et al. 2008 p.563), is SUBSTITUTED -- a' // &
    ' DECLARED DEPARTURE from the printed wording.  ZERO NEW CONSTANTS.'
  WRITE (*,'(A)') ' rec 161 FIRST-HAND COUNTER-DATUM, REPORTED NOT GLOSSED: Nagelmueller et al.' // &
    ' (2016) J Exp Bot 67(6):1897-1906 DOI 10.1093/jxb/erw003 compares the three candidate thermometers' // &
    ' directly in its WHEAT FIELD experiment and ranks soil LAST: "air temperature at plant height (5cm' // &
    ' above ground, in field experiments) and temperature in the meristem ... explained leaf elongation' // &
    ' better (higher R2-values) THAN AIR TEMPERATURE 2 M ABOVE GROUND OR SOIL TEMPERATURE (5CM BELOW' // &
    ' GROUND)."  This CORROBORATES that screen air at 2 m is the wrong thermometer -- this record''s' // &
    ' premise -- but names AIR AT PLANT HEIGHT, not soil, as the best measured replacement.  No direct' // &
    ' contradiction (Abichou''s and Evers''s rule is for the period while the APEX IS BELOW GROUND, and' // &
    ' Nagelmueller tracked leaves already above ground), BUT THIS RECORD''S WINDOW IS LONGER THAN THE' // &
    ' EVIDENCE: terminal spikelet is weeks after the apex rises above the surface.'
END IF
! rec 162 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (canopy_clock_basis >= 1) THEN
  WRITE (*,'(A)') ' rec 162 canopy_clock_basis = 1 : the CANOPY clock is advanced by the VF/DF-gated' // &
    ' development increment current_dtdu the stage engine already uses, not by the raw increment' // &
    ' current_dtt, so development and canopy run on ONE time basis instead of two.  The law is already' // &
    ' coded and already cited: McMaster GS, White JW, Hunt LA, Jamieson PD, Dhillon SS,' // &
    ' Ortiz-Monasterio JI (2008) Ann Bot 102(4):561-569, DOI 10.1093/aob/mcn115, Eq. (1) p.563' // &
    ' DTDU = DTT x [minimum(VF,DF)], Eqns (2)-(3) p.564, Eq. (4) p.565.  ZERO new equations and ZERO' // &
    ' new constants: only which EXISTING variable the canopy module is handed.'
  WRITE (*,'(A)') ' rec 162 LABELLED AN INFERENCE, NOT PRINTED ARCHITECTURE: McMaster et al.' // &
    ' (2008)''s printed scope is STAGE PROGRESSION from emergence to terminal spikelet (pp.564-565) and' // &
    ' THE PAPER DOES NOT PRINT THAT LEAF APPEARANCE IS VF/DF-GATED.  What it does print is a phase' // &
    ' duration expressed in phyllochrons OF THE TDU CLOCK (p.563) and, p.567, that "development rates' // &
    ' for many processes (both phenology and others such as leaf appearance)" respond to temperature.' // &
    '  The WINDOW (emergence to terminal spikelet, 400 degCd TDU) IS verbatim faithful.  The' // &
    ' PHYLLOCHRON IS HELD at the compiled 73.9 degCd (A3b): the A3a rescale to' // &
    ' 54.6 degCd is REJECTED because it is 20 % outside the 68.4-77.4 degCd printed bracket.'
  IF (canopy_clock_driver >= 1) THEN
    WRITE (*,'(A)') ' rec 161 + rec 162 BOTH ON: the two records compose multiplicatively (a' // &
      ' thermometer and a gate).  The gate factor phenology applied today is recovered as' // &
      ' current_dtdu / current_dtt and multiplied onto the SOIL-driven increment; on a day with' // &
      ' current_dtt = 0 the ratio is undefined and the gate is rebuilt from current_vf / current_df' // &
      ' under the record-150 window.'
  END IF
END IF
! rec 163 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (postanth_uptake_regulation >= 1) THEN
  WRITE (*,'(A)') ' rec 163 postanth_uptake_regulation = 1 : AFTER ANTHESIS the construction-DEMAND gate on' // &
    ' the root-N uptake potential is REMOVED and the potential is regulated solely as Bertheloot 2011a prints it' // &
    ' -- Bertheloot, Cournede & Andrieu (2011) Ann Bot 108(6):1085-1096, DOI 10.1093/aob/mcr119, Eqn 3-5,' // &
    ' f(soil [N]) x E_C x E_N with NO demand term of any kind.  U_pot = supply x E_N on days with' // &
    ' dvs2_growth >= 1, where supply is the UNCHANGED existing soil-N supply law (record 139 mass flow +' // &
    ' diffusion + active Michaelis-Menten, each capped by its layer pool) and E_N is the UNCHANGED' // &
    ' record-155 limb-(b) exp(-beta_N [N]mob), beta_N = 130 (Bertheloot 2011b Table 2).  ZERO new constants.'
  WRITE (*,'(A)') ' rec 163: post-anthesis root uptake supplies' // &
    ' 3.7-22.7 % of grain N in the model (16-deck mean 12.1 %) where the SAME Braunschweig decks measured' // &
    ' 26.7-48.1 % (mean 35.6 %) and Bertheloot 2008 p.1709 measured 39 % and 45 % independently, with' // &
    ' 44-57 kg N ha-1 of mineral N standing unused in the soil throughout the fill.'
  WRITE (*,'(A)') ' rec 163 DISCLOSED: E_C (Bertheloot 2011a Eqn 5) is NOT implemented -- its typeset form is' // &
    ' an image in the only full text on disk and its argument is a per-culm root dry-mass influx with a' // &
    ' 0.0013 g d-1 threshold that this model does not carry on that basis.  Since E_C <= 1, mode 1 is an' // &
    ' UPPER BOUND on Bertheloot 2011a uptake, which is the falsifiable direction.'
  IF (n_uptake_regulation /= 1 .AND. n_uptake_regulation /= 3) THEN
    WRITE (*,'(A,I0,A)') ' rec 163 but rec 155 n_uptake_regulation = ', n_uptake_regulation, &
      ': limb (b) is not evaluated, so E_N is identically 1.0 and the statement is the PURE supply law' // &
      ' (reported, not refused).'
  END IF
END IF
! rec 164 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (lamina_capacity_profile_basis >= 1) THEN
  WRITE (*,'(A,F5.3,A)') ' rec 164 lamina_capacity_profile_basis = 1 : the record-157 lamina holding capacity' // &
    ' sln_max_N_lnh = ', sln_max_N_lnh, ' g N m-2 is applied at the TOP of the model''s own existing' // &
    ' exponential canopy N profile instead of as the canopy MEAN, so the canopy-mean capacity each day' // &
    ' = sln_max_N_lnh / (n_top/n_mean) with n_top = canopy_n0_top(N_leaf, LAI, kN_canopy) and n_mean =' // &
    ' N_leaf/LAI -- the two quantities output/daily_canopy_n_profile.csv already prints.  ZERO new' // &
    ' constants: the same printed value, the same profile, the same extinction coefficient.'
  WRITE (*,'(A)') ' rec 164 WHY: Pask (2009) PhD thesis, University of Nottingham, p.29 prints' // &
    ' "Wheat leaves commonly accumulate contents in excess of 2.5 g m-2 (about 4% N in dry matter)' // &
    ' especially in the upper leaves of the canopy" -- an UPPER-LEAF observation; the MEASURED canopy' // &
    ' mean at anthesis under non-limiting N is 1.65 (Apache) - 1.90 (Isengrain) g N m-2 (Bertheloot,' // &
    ' Martre & Andrieu 2008 Plant Physiol 148:1707-1720 Fig. 2A p.1711) and 1.73-2.11 at the Braunschweig' // &
    ' high-N decks.  The factor is >= 1, so the capacity can only FALL.'
  WRITE (*,'(A)') ' rec 164 WELL-DEFINED BEFORE THE PROFILE EXISTS: on a day with LAI <= 1.0e-6,' // &
    ' N_leaf <= 1.0e-9, kN_canopy <= 1.0e-6 or a computed ratio below 1.0 the capacity FALLS BACK to the' // &
    ' mode-0 value sln_max_N_lnh itself and the fallback day is counted in the season-summary line.'
  IF (lamina_n_storage < 1) THEN
    WRITE (*,'(A,I0,A)') ' rec 164 but rec 157 lamina_n_storage = ', lamina_n_storage, &
      ': the lamina storage tier runs on n4_sln_target, so this record is INERT on this deck' // &
      ' (reported, not refused).'
  END IF
END IF
! rec 165 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (wsc_release_delay_basis >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 165 wsc_release_delay_basis = 1 : the rec-92 DEFICIT draw limb now also requires' // &
    ' (iday_in - wsc_day_anth) > wsc_accum_end_dpa = ', wsc_accum_end_dpa, &
    ' d -- the IDENTICAL conjunct the rec-105 = 2 release limb already carries -- so stem WSC' // &
    ' mobilization cannot begin before the accumulation window closes.  ZERO new constants: the same' // &
    ' EXISTING PARAMETER (SHARED.f90:7101), the same anthesis latch, the same draw expression.'
  WRITE (*,'(A)') ' rec 165 WHY: Schnyder 1993 New Phytol 123:233-245 p.238 prints "WSC' // &
    ' accumulate in the stem until 10-20 d after anthesis and - at that stage - the maximum amount of' // &
    ' WSC is 2-4 times greater than at anthesis", and p.238 col.2 that storage "is not competitive with' // &
    ' grain filling"; the pool is a SINK inside that window and a SOURCE only after it.  MEASURED: the' // &
    ' Braunschweig stem reserve is still RISING at the +18/20 d sample on 15 of 16 decks while this' // &
    ' model peaks at 3-8 dpa and is empty by 18-23 dpa (Manderscheid et al. 2020 ODJAR 6).'
  WRITE (*,'(A)') ' rec 165 SECOND AND THIRD FIRST-HAND WITNESSES FOR THE SAME WINDOW [WRD-1b]: Ehdaie,' // &
    ' Alloush, Madore & Waines 2006 Crop Sci 46:735-746 Table 2 p.738 / text p.737 col.2 -- lower' // &
    ' internodes peak 10 d postanthesis (822 mg), peduncle (389 mg) and penultimate (417 mg) at 20 d,' // &
    ' then decrease progressively to maturity -- and p.744 col.1, "Maximum stem weight in all genotypes,' // &
    ' except Ramona 50, was attained 20 d postanthesis"; Ruuska et al. 2006 Funct Plant Biol 33:799 p.803' // &
    ' col.2, "WSC increases from before anthesis until about midway through grain filling and then' // &
    ' decreases again"; and Gebbing, Schnyder & Kuhbauch 1998 Plant Cell Environ 21:301-313 p.302 col.1,' // &
    ' "accumulation may continue until 10-20 d after anthesis", with p.310 col.2 tying the stem-weight' // &
    ' loss "from the time of maximum WSC content, at about 10-20 d after anthesis, to maturity" to WSC' // &
    ' mobilization.  THREE independent first-hand sources print the same window; no constant is taken' // &
    ' from any of them.'
  WRITE (*,'(A)') ' rec 165 NOT CHANGED: the draw is still the grain''s own unmet demand flux with NO' // &
    ' rate ceiling invented and none removed (Schnyder p.239: mobilization is "induced in response to a' // &
    ' deficiency in current photosynthate supply to grains", and water stress/shading at maximum WSC' // &
    ' "failed to accelerate the loss"); reserve_conv_eff, the deposition side, the capacity cap, the' // &
    ' grain demand law and every N, water, energy and phenology statement are untouched.'
  IF (carbon_reserve_mode /= 1) THEN
    WRITE (*,'(A,I0,A)') ' rec 165 but rec 55 carbon_reserve_mode = ', carbon_reserve_mode, &
      ': there is no stem WSC pool, so this record is INERT on this deck (reported, not refused).'
  ELSE IF (wsc_reserve_law == 2 .AND. wsc_fate_law /= 1 .AND. wsc_fate_law /= 3) THEN
    WRITE (*,'(A,I0,A,I0,A)') ' rec 165 but rec 105 wsc_reserve_law = ', wsc_reserve_law, &
      ' with rec 148 wsc_fate_law = ', wsc_fate_law, &
      ': the rec-92 deficit limb is already unreachable (its HEAD guard), so this record is INERT' // &
      ' on this deck (reported, not refused).'
  END IF
END IF
! rec 172 stem_reserve_deposition_n_basis -- ON-only echo (the OFF run.log is unchanged)
! rec 173 spac_sink_coupling -- ON-only echo (the OFF run.log is unchanged)
IF (spac_sink_coupling >= 1) THEN
  WRITE (*,'(A)') ' rec 173 spac_sink_coupling = 1 : the RICHARDS ROOT-SINK DEMAND is now the' // &
    ' COMMITTED SPAC leaf''s own converged hourly canopy transpiration on every outer pass >= 2,' // &
    ' and the daily ledger potential (trans_act / trans_act_pot, hence the f_t denominator) is' // &
    ' that series'' own sum.  One assignment is removed; nothing is scaled, capped or fitted.'
  WRITE (*,'(A)') ' rec 173 WHY: at mode 0 the series that forces Richards is rec 78''s' // &
    ' UNSTRESSED outer-iteration-1 series (plain Leuning, g1_leu = 10.58 WITH h(Ds), f_psi = 1),' // &
    ' while every later pass runs the SPAC leaf (g1_tuzet = 5.776, NO h(Ds)) whose transpiration is' // &
    ' LARGER -- measured on the 42-deck verification run as demand / pass-2 SPAC E /' // &
    ' delivered sink (mm): WDF25_Ruta_aCO2_WW 174.8 / 260.5 / 154.0; WDF25_Ruta_aCO2_WL' // &
    ' 176.0 / 265.4 / 151.9; WDF23_Ruta_aCO2 107.8 / 159.6 / 107.8; MARICOPA_906 274.6 / 356.5 /' // &
    ' 265.1; MARICOPA_905 246.0 / 386.1 / 192.3; MARICOPA_913 171.0 / 259.7 / 164.1;' // &
    ' BRAUNSCHWEIG_103 173.9 / 195.0 / 173.8.  Pass 2 then scaled its latent store to the' // &
    ' delivered sink (f_trans ~ 0.59 at WDF25 Ruta aCO2 WW), so the latent deficit became' // &
    ' sensible heat and the canopy warmed for a supply the soil was never asked to meet.'
  WRITE (*,'(A)') ' rec 173 THE ITERATION IS THE EXISTING ONE: the demand <-> delivered-sink' // &
    ' fixed point runs inside the EXISTING outer loop (cap 5 passes, convergence' // &
    ' |dTcan_act| < 0.05 K); NO new iteration and NO new tolerance.  spac_psi_s_MPa, the only soil' // &
    ' term the SPAC leaf reads, is a ONCE-PER-DAY quantity built from the PREVIOUS day''s committed' // &
    ' per-layer heads (rec 81 review C1), so an extra within-day iteration has nothing to converge:' // &
    ' the restored feedback closes at the DAILY step (more water asked -> drier committed profile' // &
    ' -> lower psi_s -> lower f_psi -> smaller E).  The per-day water residual |E_pass2 - sink| and' // &
    ' the pass count are written to output/daily_spac_sink_coupling.csv, not asserted.'
  WRITE (*,'(A)') ' rec 173 WHAT STAYS: f_psi(psi_L) remains the ONLY stomatal water-stress' // &
    ' locus (no f_t anywhere in assim under rec 23 = 1); rec 78''s unstressed series and' // &
    ' trans_pot_unstressed_mm stay written as the Leuning-branch potential diagnostics they are' // &
    '; and the pass-2 repartition code path is KEPT, so any residual' // &
    ' between pass 2''s potential and the delivered sink is still conserved into sensible heat and' // &
    ' still written -- unity by construction when the two agree.'
END IF
! rec 166 grain_ca_xylem_throughput -- ON-only echo (the OFF run.log is unchanged)
IF (grain_ca_xylem_throughput >= 1) THEN
  WRITE (*,'(A)') ' rec 166 grain_ca_xylem_throughput = 1 : grain and chaff calcium now arrive with the XYLEM' // &
    ' WATER THE EAR RECEIVES -- ear share of the model''s own standing green area (cse_eai_meas against' // &
    ' cls_gai + cse_sai + cse_iai, ONE basis) times actual transpiration -- split at the grain-base xylem' // &
    ' discontinuity (Zee & O''Brien 1970 Aust J Biol Sci 23:783) by the Hocking 1994 J Plant Nutr 17:1289' // &
    ' Table 1(c) p.1294 measured partition 10.8/(10.8+9.8) = 0.524272.  The grain''s own net water import' // &
    ' (Lizana 2010 J Exp Bot 61:1147 Tables 2+3, grain_h2o_g_gDM = 0.7625, UNCHANGED) is retained in full' // &
    ' and only the REMAINDER of the ear''s xylem water is partitioned, so nothing is double counted.'
  WRITE (*,'(A)') ' rec 166 WHY: at mode 0 grain Ca is the product of TWO CONSTANTS,' // &
    ' grain_h2o_g_gDM x xylem_Ca_mgL_CAL = 0.7625 x 842.0246 = 642.04 mg kg-1 on every WDF25 deck' // &
    ' (0.7625 x 654.339 = 498.93 at WDF23), so no CO2, water, nitrogen, cultivar or season response is' // &
    ' possible; and 842.0246 mg L-1 is not a measurement but the linear rescale of the sourced calcite' // &
    ' bracket start 54.9 mg L-1 onto deck WDF25_Ruta_aCO2_WW''s OWN OBSERVED grain Ca of 642.043 mg kg-1' // &
    ' (BATCH11 rec-131 result note Sect. 3).  21.0 mmol L-1 of xylem-sap Ca is 4.6x the highest monocot' // &
    ' value measured by Atkinson, Ruiz & Mansfield 1992 J Exp Bot 43:1315 Table 3 p.1319 (Zea mays' // &
    ' 4.6 +- 0.6 mol m-3; no wheat in the table), and it implies a season crop Ca uptake of 129.7 g m-2' // &
    ' against the 2.32 g m-2 whole plant measured at harvest by Gregory, Crawford & McGowan 1979' // &
    ' J Agric Sci 93:485 Table 1 p.491 -- a factor of 56.'
  WRITE (*,'(A)') ' rec 166 THE CAMPAIGN NUMBER IS NOT READ: calcium''s solution/xylem concentration' // &
    ' falls through to each deck''s OWN UNCHANGED driver law -- the calcite ceiling_mg_L = 54.9 mg L-1 (1.37 mmol L-1,' // &
    ' Sci Rep 2025 Table 2, PROVISIONAL) at WDF and Maricopa, the linear Kd partition at Braunschweig; the' // &
    ' calibrated column xylem_Ca_mgL_CAL is not read for Ca on this mode.  This record REMOVES a' // &
    ' calibrated constant and introduces ONE measured partition, no multiplier and no fitted number.'
  WRITE (*,'(A)') ' rec 166 FALSIFIERS: Gregory 1979 Table 1 p.491 ear Ca 0.08 (anthesis) -> 0.18' // &
    ' -> 0.27 g m-2 (harvest), i.e. 70 % of final ear Ca arrives after anthesis by xylem; Hocking 1994' // &
    ' Table 6 p.1297 grain Ca 0.52 and chaff Ca 1.51 mg g-1 DM at maturity; Hocking Table 4 p.1296, Ca' // &
    ' redistribution below the 10 % reporting threshold while N is 73.0 %; White & Broadley 2003 Ann Bot' // &
    ' 92:487 for phloem immobility; Gilliham 2011 J Exp Bot 62:2233 Abstract p.2233 for the' // &
    ' transpiration coupling and the irreversibility of deposition.'
  WRITE (*,'(A)') ' rec 166 NOT REPAIRED, DISCLOSED: whole-crop Ca uptake is still transpiration x' // &
    ' solution with no root-exclusion term (8.5 against Gregory''s 2.32 g m-2 -- 3.6x, down from 56x);' // &
    ' Atkinson 1992 Abstract p.1315 itself reports that xylem Ca transport "was not primarily by mass' // &
    ' flow", a disagreement declared and not resolved; the ear enters on a PROJECTED area basis against' // &
    ' one-sided leaf area with no ear stomatal conductance; and no grain-Ca observation exists at' // &
    ' Braunschweig or Maricopa, so 26 of the 42 decks carry no falsifier for this route.'
END IF
! rec 167 ON banner. Nothing is
! printed at mode 0, so a mode-0 run.log is unchanged.
IF (lamina_synthesis_sheath_basis >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 167 lamina_synthesis_sheath_basis = ', lamina_synthesis_sheath_basis, &
    ' : the Bertheloot 2011a sigma_lam now multiplies the LAMINA PLUS SHEATH green mass and the SAME sheath mass is' // &
    ' removed from the stem''s sigma basis, so no gram of dry mass synthesizes photosynthetic N twice.' // &
    '  ZERO constants changed -- only the mass each PRINTED sigma multiplies.'
  WRITE (*,'(A)') ' rec 167 WHY: Bertheloot, Cournede & Andrieu 2011 Ann Bot 108:1097-1109 (DOI' // &
    ' 10.1093/aob/mcr125) Table 2 p.1101 prints sigma^Nph_tp,i = "0.00015*, 0.00003+" with its OWN' // &
    ' footnotes on the same page, verbatim: "* Laminae and sheaths." and (dagger) "Internodes, peduncle' // &
    ' and chaff."  The model''s W_green is the LAMINA ALONE and its lumped W_stem CONTAINS the sheaths,' // &
    ' so mode 0 has the sheath under the WRONG FOOTNOTE.'
  WRITE (*,'(A)') ' rec 167 WHAT IT COSTS: the rec-154 lamina limb has an analytic fixed point' // &
    ' %N* = 100*[(1+r)*(sigma/delta)*f_N*f_PAR*(M_basis/W_green) + n2_struct_lam] which on a lamina-only' // &
    ' basis returns 2.22-2.26 % at BRAUNSCHWEIG_103/106/115/118 on the +16-18 d sample (model 2.23-2.27,' // &
    ' MEASURED 3.42-3.65) and whose ceiling at f_N = f_PAR = 1 is 2.85 %, i.e. the measurement is' // &
    ' UNREACHABLE on that basis at any f.'
  IF (lamina_synthesis_sheath_basis == 2) THEN
    WRITE (*,'(A,F7.4,A)') ' rec 167 MODE 2: M_sheath = ', lsb_a_sheath, &
      ' x W_green, the PRINTED MEASURED dry masses of Bertheloot 2011b Table 1c p.1100 at flowering -- sheath' // &
      ' M^struct 103+69+43+9 = 224 plus M^rem 22+15+9+2 = 48 = 272 mg culm-1, lamina 320 + 100 = 420 mg' // &
      ' culm-1.  A measured allometry, held FIXED, never fitted.'
  ELSE
    WRITE (*,'(A)') ' rec 167 MODE 1: M_sheath = cse_sai x (W_green/lai), the model''s OWN' // &
      ' green sheath area index (rec 68) priced at its OWN green-lamina mass per leaf area.  DECLARED A' // &
      ' LOWER BOUND: Bertheloot 2011b Table 1c p.1100''s own areas make the measured sheath 16.0/3.91 = 4.1x the' // &
      ' lamina''s mass per unit area.  Days with lai <= 0 or cse_sai <= 0 fall back to the mode-0 basis' // &
      ' and are COUNTED (lsb_days_no_area).'
  END IF
  WRITE (*,'(A)') ' rec 167 NOT CHANGED: sigma_lam = 1.5e-4 and sigma_stem = 3.0e-5 keep their' // &
    ' printed Bertheloot 2011b Table 2 p.1101 values; delta = 0.008, k1 = 0.0018, k2 = 864000 and r = 0.10 are' // &
    ' untouched; the DEGRADATION side keeps n2_struct_lam = 0.0079 g N g-1 on the LAMINA mass, which is' // &
    ' the basis Bertheloot et al. 2008 Plant Physiol 148:1707 p.1711 measured it on and the basis the' // &
    ' pre-registered prediction was evaluated with; the rec-87/97 lamina capacity clip STAYS and is' // &
    ' REPORTED (ont_days_cap_clip), never assumed non-binding.'
  WRITE (*,'(A)') ' rec 167 PRE-REGISTERED, MODE 2: canopy-mean SLN at +16-18 d 1.62/1.64/1.66/1.65' // &
    ' g N m-2 at BRAUNSCHWEIG_103/106/115/118 (today 1.122/1.132/1.140/1.119; MEASURED 2.11/1.81/1.99/' // &
    '1.80) -- about half the gap closed and still BELOW the observation.  Rows named at risk in advance:' // &
    ' grain_N_pct/grain_protein_pct, stem_N_pct at harvest, AGB at harvest.'
  IF (organ_n_turnover_nema < 1) THEN
    WRITE (*,'(A)') ' rec 167 INTERNAL ERROR: reached with rec 154 = 0; the refusal above should have' // &
      ' stopped this run.'
  END IF
END IF
! rec 168 mineral_impedance_basis -- ON-only echo (the OFF run.log is unchanged)
IF (mineral_impedance_basis >= 1) THEN
  WRITE (*,'(A)') ' rec 168 mineral_impedance_basis = 1 : the solute DIFFUSION IMPEDANCE FACTOR inside' // &
    ' mineral_module.f90 is now Olesen, Moldrup, Yamaguchi & Rolston 2001 Soil Science 166:89-96' // &
    ' DOI 10.1097/00010694-200102000-00002, printed p.94 Eq.(7a)/(7b): f_l = 1.1*(theta - theta_th) for' // &
    ' theta > theta_th and 0 otherwise, with theta = vwc_root and theta_th from the same paper''s printed' // &
    ' p.93 Eq.(5) and THIS DECK''S OWN measured clay fraction, silt fraction and topsoil bulk density' // &
    ' (driver_N.csv positions 68-70).  TWO LOCI, BOTH IN mineral_module.f90: the live rec-140' // &
    ' diffusion impedance fl_theta_m, and the already-retired rec-131 active term.'
  WRITE (*,'(A)') ' rec 168 WHY: at mode 0 BOTH loci compute f_l = MAX(0, 1.60*vwc_root - 0.17),' // &
    ' a pair re-labeled 2026-10-18 to Barraclough & Tinker 1981 J Soil Sci 32:225 p.234 Fig. 1 best fit' // &
    ' 1.58*theta - 0.17 (intercept exact, slope 1.3 % off); the Oates & Barber 1987 attribution was REMOVED on 2026-10-07 as' // &
    ' unobtainable on every legitimate route (see the development record),' // &
    ' and 0.17 is in fact the printed STANDARD DEVIATION of Olesen''s slope, not an intercept.  The pair' // &
    ' is also SITE-BLIND: identical slope and intercept on the WDF loamy sand (clay 0.04, silt 0.03,' // &
    ' rho_b 1.63), the Braunschweig sand (0.07/0.24/1.55) and the Maricopa clay loam (0.34/0.31/1.28),' // &
    ' so the entire between-site contrast in diffusive mineral supply was carried by soil WATER and none' // &
    ' of it by soil TEXTURE.'
  WRITE (*,'(A)') ' rec 168 ZERO NEW NUMERALS: nud_csif_slope = 1.1 (p.94 Eq.7a) and' // &
    ' nud_th_c1..c6 = 0.81, -0.90, -0.07, -0.60, 0.22, 0.42 (p.93 Eq.5) are the SAME constants record' // &
    ' 139 mode 3 already carries (SHARED.f90 [NUD-9a]); theta_th is computed once at init from the' // &
    ' deck''s own measured texture.  No multiplier, no fitted coefficient, nothing refitted.  The' // &
    ' module already computed De_e = Dl*theta*f_l/b, which IS Olesen Eq.(8a) over the buffer power, so' // &
    ' this is a two-constant swap inside an already-correct functional form.'
  WRITE (*,'(A)') ' rec 168 NOT TOUCHED: record 139''s NITROGEN impedance in main_MMCO2Crop.f90 is' // &
    ' UNCHANGED (option 2 for the N chain is a separate later side copy), as are mass' // &
    ' flow, the Fe phytosiderophore term, the Zn HATS bound, mm_b, mm_Dl_e, nud_SRL, mm_depth_m, the' // &
    ' remobilization law and every export fraction.  mineral_module.f90 is MG0-isolated, so no carbon,' // &
    ' nitrogen, water, energy or phenology state is reachable from either edited statement.'
  WRITE (*,'(A)') ' rec 168 THE TRADE, STATED NOT SMOOTHED: the slope falls 1.60 -> 1.1 (x0.6875),' // &
    ' so the sensitivity of diffusive supply to the eCO2 soil-water signal is DAMPED by 1.45x and this' // &
    ' record slightly WORSENS the eCO2 grain-mineral CONTENT shortfall while removing the only unsourced' // &
    ' constant that carried site contrast.  Grain mineral LEVELS will need the avail_A re-centering' // &
    ' (an input change; see INPUT_GUIDE.md), which this record deliberately' // &
    ' does NOT carry.'
END IF
! rec 169 zn_uptake_kinetic_basis -- ON-only echo (the OFF run.log is unchanged)
IF (zn_uptake_kinetic_basis >= 1) THEN
  WRITE (*,'(A)') ' rec 169 zn_uptake_kinetic_basis = 1 : the ZINC active term is now a SATURATING FUNCTION OF' // &
    ' THE ROOT-SURFACE CONCENTRATION, S_hats = mm_field_rate(Zn) x MAX(0,W_root) x a_root x c_surf /' // &
    ' (km_zn_mgL + c_surf), with c_surf = soln - S_diff / (RLD_e x De_e x mm_depth_m x 86400 x 1e3), i.e. the' // &
    ' bulk solution concentration MINUS the Barber-Cushman depletion the diffusion limb has already' // &
    ' charged, and with NO MIN against S_mf + S_diff.'
  WRITE (*,'(A)') ' rec 169 WHY: at mode 0 the "ceiling" S_hats = MIN(S_hats_raw, S_mf + S_diff) BINDS ON' // &
    ' 7499 OF 7499 ACTIVE DECK-DAYS across all 42 decks, so it is not a bound but the identity S_hats =' // &
    ' S_mf + S_diff, and the day''s raw zinc supply is EXACTLY 2 x (mass flow + diffusion) -- a factor of two' // &
    ' dressed as a transporter. It doubles the level (model +22.73 mg kg-1 above Carducci et al. 2025''s' // &
    ' printed wheat grain-Zn regression where the observations sit +10.59, 20 decks) and it makes zinc''s' // &
    ' eCO2, water and nitrogen responses a scaled copy of the mass-flow response, so zinc cannot differ from' // &
    ' iron -- and it does not (at most 5 pp apart in all six Braunschweig eCO2 cells), which Gojon et al.' // &
    ' 2022 Trends Plant Sci 27:882 DOI 10.1016/j.tplants.2022.04.007 p.2 rules out as a mechanism.'
  WRITE (*,'(A)') ' rec 169 THIS IS A DECLARED CONSTRUCTION, NOT A MEASUREMENT: no Michaelis-Menten zinc' // &
    ' influx kinetic for FIELD-grown wheat per unit root mass or root surface area exists on disk or, within' // &
    ' the information available to this module, in print.' // &
    ' The SHAPE (a high-affinity transporter half-saturated at 2.3 nM) is taken' // &
    ' from Hacisalihoglu, Hart & Kochian 2001 Plant Physiol 125:456 DOI 10.1104/pp.125.1.456 Table II p.459' // &
    ' (high-affinity row, apparent Km 2.3 +- 1.0 nM for (-)Zn-grown cv BDME-10; the paper''s printed range is' // &
    ' 0.6-2.3 nM), which is a 21-d-old HYDROPONIC SEEDLING parameter and whose companion Vmax was' // &
    ' found inconsistent with the model''s own earlier test runs (grain Zn 50-80x too high). The scale is' // &
    ' mm_field_rate(Zn) = 2.36e-3 mg Zn gDW-1 d-1, the EXISTING driver_minerals.csv column-15 value derived' // &
    ' from Liu, Liu, Zhang, Chen & Zou 2019 Front Plant Sci 10:426 DOI 10.3389/fpls.2019.00426 -- Table 6' // &
    ' p.8 shoot Zn at maturity 259 g ha-1 + Table 4 p.6 root Zn 16.1 g ha-1 = 27.51 mg m-2, / Table 3 p.6' // &
    ' RDW 483 kg ha-1 = 48.3 g m-2, / a 241-d season (p.3, sown Oct 10/11, harvested Jun 8/10). Neither' // &
    ' number is fitted to any campaign observation and NO NEW NUMERAL IS INTRODUCED.'
  WRITE (*,'(A)') ' rec 169 GUARD: IF c_surf <= 0 THE ACTIVE TERM IS 0.0 EXACTLY -- a root surface at or' // &
    ' below zero available concentration supports no influx. Where the water-equivalent diffusive flux' // &
    ' density q_diff = RLD_e x De_e x mm_depth_m x 86400 x 1e3 is <= 0 (no root, zero effective diffusion' // &
    ' coefficient, or a soil water content at which the impedance factor is zero) NO diffusive depletion has' // &
    ' been charged, so c_surf falls back to the bulk solution concentration; that is a well-posedness guard' // &
    ' against a division by zero, not a second law. Both day counts are printed in the season summary.'
  WRITE (*,'(A)') ' rec 169 NOT CHANGED: Fe, Cu, Mn, Mg and Ca are untouched; the mass-flow, diffusion,' // &
    ' phytosiderophore and retired limbs, the pool cap, the organ allocation, the remobilization' // &
    ' route and the export fractions are untouched; records 131 and 140 are untouched; the nitrogen chain''s' // &
    ' own record-139 supply law is NOT touched (this record lives inside mineral_module.f90 only, which' // &
    ' keeps its MG0 isolation and writes no C, N, water or energy state). INERT unless rec 131' // &
    ' mineral_supply_law = 1 and rec 140 mineral_supply_mech >= 1 put the zinc limb on the path.'
  WRITE (*,'(A)') ' rec 169 NOT REPAIRED, DISCLOSED: avail_A(Zn) is NOT re-centered by this patch, so' // &
    ' the Braunschweig zinc LEVEL on BR102/BR114 is expected to move off its calibration until the lead' // &
    ' re-centers it as an INPUT change; the zinc remobilization share (model 26-34 % against Waters, Uauy,' // &
    ' Dubcovsky & Grusak 2009 J Exp Bot 60:4263 p.4267 Table 2''s 42.6 %) is owned by repair R5b, not by' // &
    ' this record; and Mn and Cu still have no active uptake term of any kind.'
END IF
! rec 171 mineral_n_coupling_basis -- ON-only echo (the OFF run.log is unchanged)
IF (mineral_n_coupling_basis >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 171 mineral_n_coupling_basis = ', mineral_n_coupling_basis, &
    ' : the grain''s PHLOEM SINK CAPACITY for a phloem-mobile mineral is set by the grain''s PROTEIN' // &
    ' (nitrogen), at the per-element grain element:N stoichiometry r_e = driver_minerals.csv column 20' // &
    ' grain_eN_mg_gN [mg element per g grain N].  1 = FLUX (daily) basis, the design note''s R5 leg A' // &
    ' wording: the day''s grain loading <= r_e x MAX(0, dN_grain today).  2 = POOL (cumulative) basis,' // &
    ' Xia''s own wording: the grain''s accumulated content <= r_e x N_grain today.  Refused loading is' // &
    ' RETURNED TO ITS SOURCE (uptake share -> leaf terminus, senescence share -> dead tissue), exactly' // &
    ' as record 126 already does, so the per-element daily mass balance is unchanged.'
  WRITE (*,'(A)') ' rec 171 SOURCE: Xia, Wang, Qiao, Kong, Xue, Wang, Kong, Xue & Sizmur 2020,' // &
    ' Food and Energy Security 9:e243, DOI 10.1002/fes3.243, printed p.7 of 22: "Analysis of the' // &
    ' speciation and localization of Zn in wheat grains indicates that Zn interacts with proteins, and' // &
    ' therefore, grain protein content (as a physiological sink for Zn) may determine the sink capacity";' // &
    ' printed pp.6-7 for the SHAPE: "zero or relatively low N supply results in low wheat grain Zn' // &
    ' concentration, whereas optimal N supply achieves high grain Zn concentration, but excess N supply' // &
    ' does not further increase grain Zn, grain N, or protein concentrations"; printed p.7 for the three' // &
    ' steps N acts on.  Corroborated first hand by Cakmak 2008 Plant Soil 302:1 printed p.6 and by' // &
    ' Waters et al. 2009 J Exp Bot 60:4263 pp.4267-4268, already this module''s senescence-export basis.'
  WRITE (*,'(A)') ' rec 171 THE CONSTANT IS DECLARED ABSENT: grain_eN_mg_gN ships 0.0 (= no bound)' // &
    ' for all six elements, the same convention record 126''s grain_max_mg_kg already uses, because NO' // &
    ' first-hand wheat-grain primary on disk prints a grain element:N or element:protein stoichiometry' // &
    ' (null-safe sweep of both textbook trees; Kutman 2010, Peleg 2008, Fan 2008 and Zhao 2009 are absent).' // &
    ' Taking r_e from the campaign''s own measured grain would be fitting to the falsification target and' // &
    ' is forbidden.  With all six at 0.0 this record is a PROVABLE NO-OP at every mode.'
  WRITE (*,'(A)') ' rec 171 NOT CLAIMED: nitrogen enters as a CEILING only, so the N-driven' // &
    ' INCREASE in total uptake that Xia p.7 also reports is NOT modeled and the N-dose CONTENT magnitude' // &
    ' shortfall is not addressed; root uptake, root-to-shoot transport, the senescence export fractions,' // &
    ' the soil impedance pair and the Zn HATS bound are untouched (other repairs, other records); and' // &
    ' because grain protein FALLS under eCO2 a binding ceiling makes the eCO2 CONTENT response slightly' // &
    ' MORE negative -- this record and the eCO2 magnitude shortfall pull in opposite directions.'
  WRITE (*,'(A)') ' rec 171 ISOLATION: N_storage and N_chaff enter mineral_module.f90 READ-ONLY' // &
    ' (N_grain = MAX(0, N_storage - N_chaff), the SAME expression main_MMCO2Crop.f90 already uses for' // &
    ' grainN_conc_b3 and grain_protein_pct) and appear on no left-hand side there; nothing is written to' // &
    ' the nitrogen chain or the carbon chain, and the module''s MG0 isolation is unchanged.'
END IF
! rec 174 ON banner. Nothing is printed at mode
! 0, so a mode-0 run.log is unchanged.
IF (sink_full_carbon_fate == 1) THEN
  WRITE (*,'(A)') ' rec 129 sink_full_carbon_fate = 1 : while the grain sink is full, the rec-92 deposition window is' // &
    ' lifted for that surplus (storage continues to the existing wsc_cap_frac capacity); the remainder is booked to' // &
    ' an explicit, reported unallocated-assimilate carbon loss (cfo_loss_day/cum), not to respiration and not to stem' // &
    ' structure, superseding the rec-127 residual for this surplus. No wheat overflow-respiration rate law exists in' // &
    ' the library on disk, and photosynthetic down-regulation by a full sink is not supported in wheat (Schnyder 1993' // &
    ' p.238); neither is modeled. SHARED.'
END IF
! rec 175 ON banner. Nothing is printed at mode 0.
IF (grain_n_partition_basis >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 175 grain_n_partition_basis = ', grain_n_partition_basis, &
    ' : mode 1 corrects an earlier formulation.  Record 158''s own header declares that ndb_rate_prev is' // &
    ' "continuous across anthesis"; record 154 sets N_from_assim_day = 0 by construction, the three' // &
    ' N_alloc_* terms are already zeroed post-anthesis, so the sum is IDENTICALLY ZERO after anthesis,' // &
    ' N_max = MIN(1.21, 0) = 0 and the shoot stops reducing nitrate the day after anthesis.'
  WRITE (*,'(A)') ' rec 175 MODE 1: the sum gains the record-154 limb''s OWN de novo reduction' // &
    ' demand on the day -- grain import + lamina + stem + chaff synthesis, every term already computed' // &
    ' and already printed in output/daily_organ_n_turnover.csv.  Busch, Sage & Farquhar 2018 Nature' // &
    ' Plants 4:46-54 (DOI 10.1038/s41477-017-0065-x) p.49 limiter (b), the DEMAND, which is the limiter' // &
    ' record 158 mode 1 elected.  ZERO NEW CONSTANTS; the Busch ceiling N_max = 1.21 umol m-2 s-1' // &
    ' (printed in the Fig. 3a parameter block p.52) stays outermost.'
  IF (grain_n_partition_basis >= 2) THEN
    WRITE (*,'(A)') ' rec 175 MODE 2: the model applies FPB Eqn 12''s (p.785) RATE and discards' // &
      ' its INITIAL CONDITION.  Bertheloot, Andrieu, Fournier & Martre 2008 Funct Plant Biol 35:781-796' // &
      ' (DOI 10.1071/FP08064) p.786 verbatim: "grain N mass at anthesis, N_g^tot(0), was set at 2.4 mg' // &
      ' (Girard 1997)", Table 1 p.784 "N_g^tot(0) ... kg ... 2.4 x 10^-6" -- the number the model ALREADY carries as' // &
      ' gnf_ng0_per_ear and ALREADY uses to seed the potential.  Mode 2 charges the deficit to the' // &
      ' spike''s own N inside the ear (N_chaff down, N_storage unchanged), down to the chaff structural' // &
      ' floor; a refused draw is COUNTED and PRINTED, never silent.  ZERO NEW CONSTANTS.'
  END IF
  IF (grain_n_partition_basis >= 3) THEN
    WRITE (*,'(A)') ' rec 175 MODE 3: the grain''s UNMET sink joins the synthesis' // &
      ' limbs'' pro-rata claim set (Bertheloot 2011a p.1090; Bertheloot 2011a Eqn 12' // &
      ' p.1091 = FPB Eqn 6 p.784; sigma 1.5e-4/3.0e-5,' // &
      ' p_r 0.1, delta 0.008, k1 0.0018, k2 864000, all Bertheloot 2011b Table 2 p.1101 and all already in the' // &
      ' code).  PRE-REGISTERED AS INERT: the grain is served FIRST from the pool, so on a short day the' // &
      ' pool is exactly 0 when the limbs resolve and on a served day the residual is 0.  A mode-3 run' // &
      ' must be BYTE-IDENTICAL to a mode-2 run; if it is not, the prediction is falsified and is reported.'
  END IF
  ! MODE 4 banner. Printed only at mode 4, so modes 0-3 keep their run.log bytes.
  IF (grain_n_partition_basis >= 4) THEN
    WRITE (*,'(A)') ' rec 175 MODE 4: the Ngtot(0) HANDOVER is paid against the GRAIN''S OWN' // &
      ' DRY MASS, not as a lump on a grain of zero mass.  Mode 2 relabels gnf_ng0 (1.550 g N m-2 at' // &
      ' WDF25_Ruta_aCO2_WW) from chaff to grain on the first fill day, on which W_grain = 0.00 g m-2,' // &
      ' so the grain''s N concentration on that day is unbounded.  FPB 2008 p.786/Table 1 p.784 print' // &
      ' Ngtot(0) = 2.4 mg with NO dry-mass clause, so mode 2 is NOT repaired in place.  Bertheloot, Wu,' // &
      ' Cournede & Andrieu 2011 Ann Bot 108:1097-1109 Table 1(a) ''Grains'' (DOI 10.1093/aob/mcr125)' // &
      ' publishes the SAME initial condition PAIRED with a dry mass: ''Total N mass N_grain^tot mg 2.4''' // &
      ' and ''Total dry mass M_grain^tot mg 120'', i.e. 0.020 g N g-1 DM = 2.0 % DM.  Mode 4 therefore' // &
      ' charges MIN(gnf_ng0, 0.020 x W_grain) cumulatively, the chaff structural floor still refusing' // &
      ' what it cannot pay.  ONE constant, a measured quotient of two cells of one table.'
  END IF
END IF
! ==================================================================== !
! Record 178 banner and THE ONLY ASSIGNMENT OF THE LIVE CHAFF FLOOR. Nothing is
! printed and nothing is assigned at mode 0, where the live floor already holds
! the stem value by its own initializer.
! ==================================================================== !
IF (chaff_structural_n_basis == 1) THEN
  n2_struct_chaff_use = n2_struct_chaff
  WRITE (*,'(A,ES12.4,A,ES12.4,A)') ' rec 178 chaff_structural_n_basis = 1: the CHAFF''S' // &
    ' structural-N floor is now its OWN measured value, ', n2_struct_chaff, ' g N g-1 DM, instead of the' // &
    ' sheath/internode value ', n2_struct_stem, ' g N g-1 DM this model applied to it.'
  WRITE (*,'(A)') ' rec 178 SOURCE, first hand: Bertheloot J, Wu Q, Cournede P-H, Andrieu B (2011)' // &
    ' ''NEMA, a functional-structural model of nitrogen economy within wheat culms after flowering. II.' // &
    ' Evaluation and sensitivity analysis'', Annals of Botany 108(6):1097-1109, DOI 10.1093/aob/mcr125,' // &
    ' TABLE 1 section (c) ''Photosynthetic entities'', column ''Chaff'': ''Structural N mass' // &
    ' N_tp,i^struct mg 1.07'' and ''Structural dry mass M_tp,i^struct mg 210'' -> 1.07/210 = 5.10e-3' // &
    ' g N g-1 DM.'
  WRITE (*,'(A)') ' rec 178 THE SAME TABLE CORROBORATES THE VALUE IT REPLACES: the mass-weighted' // &
    ' internode + peduncle structural N of Table 1(c) is (0.85+0.33+0.14+0.05+1.30)/(188+180+154+43+257)' // &
    ' = 3.25e-3 g N g-1 DM, so 3.5e-3 is RIGHT for the stem and WRONG for the chaff.  It is a MEASURED' // &
    ' floor, same paper, Methods, verbatim: ''For each lamina and chaff, structural dry mass and N mass' // &
    ' were given by the lowest values observed during the grain-filling period over all treatments.'''
  WRITE (*,'(A)') ' rec 178 SCOPE: the ten CHAFF loci only (flowering snapshot, FPB Eqn 5' // &
    ' degradation basis, chaff tissue death nph and n_to_dead, the rec-175 Ngtot(0) draw floor, the' // &
    ' rec-154 synthesis potential, the rec-86/154 limb, the legacy nJ_R0 chaff limb, and both rec-159' // &
    ' mne_nrem_chaff branches).  n2_struct_lam, the SIX stem/sheath/internode/peduncle loci of' // &
    ' n2_struct_stem, and the root floors are UNTOUCHED.  It is a THRESHOLD, not a pool: no N is' // &
    ' created or destroyed and the Gate-D organ-state N closure must NOT move.'
  WRITE (*,'(A)') ' rec 178 RESIDUAL DISCLOSED: with the chaff on its own measured floor the' // &
    ' model still holds about 1.33 g N m-2 at WDF25_Ruta_aCO2_WW against 2.806 measured.  The measured' // &
    ' WDF25 chaff retains 1.88 % DM at harvest and NO paper on this disk has a wheat chaff ending grain' // &
    ' filling above Bertheloot 2011b''s own measured 0.51 % DM floor, so the remaining gap has no sourced' // &
    ' mechanism and none is invented (section C.2).'
END IF
! ==================================================================== !
! Record 178 modes 2 and 3, the remobilization-efficiency law: the only
! assignment of the live chaff concentration floor on these modes, and the only
! banner. Modes 0 and 1 do not enter. Inactive alternative, see Part 1 Table 6.
! ==================================================================== !
IF (chaff_structural_n_basis >= 2) THEN
  IF (chaff_structural_n_basis == 2) THEN
    n2_struct_chaff_use = n2_struct_chaff_pask
  ELSE
    n2_struct_chaff_use = n2_struct_chaff
  END IF
  WRITE (*,'(A,I0,A,F6.3,A)') ' rec 178 chaff_structural_n_basis = ', chaff_structural_n_basis, &
    ': THE CHAFF''S POST-ANTHESIS RELEASE IS NOW BOUNDED BY ITS MEASURED' // &
    ' N-REMOBILIZATION EFFICIENCY.  N_chaff(maturity) >= (1 - NRE) x N_chaff(anthesis), NRE =', &
    cnr_nre_chaff, ' , so the N-dependence of the chaff residual is EMERGENT from the anthesis' // &
    ' ear-N pool this model already computes -- not a parameterised concentration.'
  WRITE (*,'(A)') ' rec 178 NRE SOURCE, first hand: Pask AJD, Sylvester-Bradley R, Jamieson PD,' // &
    ' Foulkes MJ (2012) ''Quantifying how winter wheat crops accumulate and use nitrogen reserves' // &
    ' during growth'', Field Crops Research 126:104-118, DOI 10.1016/j.fcr.2011.09.021, TABLE 4' // &
    ' p.110 column ''CH'' of the N-remobilisation-efficiency block, whose footnote a DEFINES NRE as' // &
    ' ''percentage of N in crop or crop component at anthesis not recovered in straw or straw' // &
    ' component at harvest'' -- this law''s own form.  Chaff means: Nil 63.3, Rec 55.7, Max 51.0 %' // &
    ' (SED 3.7); abstract p.104 ''chaff (56%)''.'
  WRITE (*,'(A)') ' rec 178 NRE CORROBORATED INDEPENDENTLY: Barraclough PB, Lopez-Bellido R,' // &
    ' Hawkesford MJ (2014) Field Crops Research 156:242-248 p.247, 20 varieties, two N rates:' // &
    ' ''Vegetative-N and chaff-N averaged 6 kg/ha at N0, and 23 and 19 kg/ha, respectively, at' // &
    ' N200'' against anthesis ear N of 0.23 x 189 = 43.5 (N200) and 0.30 x 46 = 13.8 (N0) kg/ha,' // &
    ' i.e. NRE = 0.563 at N200 and 0.565 at N0 -- the SAME efficiency at a fourfold N supply range.'
  WRITE (*,'(A,ES12.4,A)') ' rec 178 HARD FLOOR UNDER THE LAW: the concentration floor in' // &
    ' force is ', n2_struct_chaff_use, ' g N g-1 DM.  Mode 2 = Pask 2012''s own measured chaff' // &
    ' structural minimum (Methods 2.9(i) p.107-108 and Results p.111, ''chaff (0.38%; TT06)'', the' // &
    ' lowest harvest chaff N% in ANY plot at N0 across three site-seasons); mode 3 = Bertheloot 2011b''s' // &
    ' 5.10e-3 (Bertheloot et al. 2011 Table 1(c) column Chaff), which the Braunschweig N1 decks' // &
    ' falsify as a universal minimum (they MEASURE 0.334-0.379 % DM).'
  WRITE (*,'(A)') ' rec 178 CONSTANCY ACROSS N LEVELS TESTED, NOT ASSUMED: Pask''s N effect' // &
    ' on chaff NRE is significant (P < 0.001) with the means DECREASING in N supply (63.3 -> 55.7' // &
    ' -> 51.0 %), while Barraclough''s two N rates differ by 0.002.  A single constant 0.56 holds' // &
    ' to about +-0.07 over the measured N range and is CONSERVATIVE at low N against Pask''s own' // &
    ' trend.  No N-dependent NRE is implemented: no source on this disk gives its form.'
  WRITE (*,'(A)') ' rec 178 SCOPE: NINE of the ten chaff loci read' // &
    ' MAX(concentration floor, (1 - NRE) x N_chaff(anthesis)).  The tenth, the chaff tissue-death' // &
    ' n_to_dead, does NOT: that is structural N leaving WITH dead dry mass, not a remobilization,' // &
    ' so the ratio law has no claim on it.  One snapshot per season, taken on the first live fill' // &
    ' day BEFORE the rec-175 handover and BEFORE that day''s Eqn-5 degradation.  It is a THRESHOLD,' // &
    ' not a pool: no N is created or moved and the Gate-D organ-state closure must NOT move.'
  WRITE (*,'(A)') ' rec 178 NOT USED, DECLARED: Kichey et al. 2007 Field Crops Res 102:22 reports' // &
    ' a chaff NRE of 0.73 by 15N labelling (as CITED by Gaju et al. 2014 Field Crops Res 155:213' // &
    ' p.214); that paper is not on this disk, so it is recorded as the upper bound of the' // &
    ' published range only.  A HIGHER NRE gives a LOWER floor, so 0.56 is the choice that CAN bind.'
END IF
! rec 176 ON banner. Nothing is printed at mode 0.
! rec 148 wsc_fate_law -- ON-only echo (the OFF run.log is unchanged)
IF (wsc_fate_law >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 148 wsc_fate_law = ', wsc_fate_law, &
    ' : 1 = the law (term A + term B), 2 = term A only, 3 = term B only. Term A restores the rec-92' // &
    ' DEPOSITION CLOSE for the grain-sink-full surplus (Schnyder 1993 New Phytol 123:233 p.238, WSC' // &
    ' accumulate "until 10-20 d after anthesis"; p.233 Summary, "terminates within three weeks following' // &
    ' anthesis"), the refused carbon travelling the existing rec-127/129 overflow path to cfo_loss_day/cum;' // &
    ' the rec-129 window LIFT and its withdrawn justification are superseded (misquote corrected at' // &
    ' SHARED). Term B opens the rec-105 = 2 release at TERMINAL GRAIN-SINK CLOSURE' // &
    ' (gfb_end_bound /= 0; Schnyder 1993 p.240, grain sink capacity "fixed at 2-3 wk after anthesis and' // &
    ' ... sink-limited thereafter") and spreads the standing pool to maturity entirely to cres_respired at' // &
    ' unity (p.240, whole-plant respiration rises while the stem WSC is lost), the grain never fed (room = 0' // &
    ' at and after that gate, so gain2 = 0 by the existing arithmetic). Zero new numeric constants.' // &
    ' SHARED.'
END IF
IF (grain_n_min_conc >= 1) THEN
  WRITE (*,'(A,I0,A)') ' rec 130 grain_n_min_conc = ', grain_n_min_conc, &
    ' : grain DM filling is capped so cumulative grain N never implies a concentration below 1.5% (Sinclair & Amir 1992,' // &
    ' via Jamieson & Semenov 2000 FCR 68:21 p.23-24, "grain growth can be limited directly by lack of N");' // &
    ' 1 = cumulative room = N_grain/0.015 - W_grain; 2 = guarded daily increment dW_grain <= dN_grain(t-1)/0.015, not' // &
    ' evaluated below 10% of the rec-118 potential kernel mass. A THIRD ceiling alongside rec 118 / rec 125 (MIN of all' // &
    ' active ceilings'' room); the refused DM flows into the existing rec 127/129 fate chain unchanged. SHARED.'
END IF
IF (mineral_supply_law == 1) THEN
  WRITE (*,'(A)') ' rec 131 mineral_supply_law = 1 : Fe/Zn/Cu/Mn soil-solution concentration + a generalized' // &
    ' per-root-mass active term are driven by ONE calibrated site-level availability parameter each (driver_minerals.csv' // &
    ' column 17 avail_A_mgL_CAL); Ca xylem concentration driven by ONE calibrated parameter (column 18' // &
    ' xylem_Ca_mgL_CAL), replacing the borrowed 54.9 mg/L analogue-soil ceiling; the Lizana 2010 grain-water import ratio' // &
    ' (0.7625 g H2O g-1 final grain DM) stays measured, unchanged. Mg untouched. SHARED/mineral_module.f90.'
END IF
IF (mineral_supply_mech == 1) THEN
  WRITE (*,'(A)') ' rec 140 mineral_supply_mech = 1 : rec 131''s calibrated generalized active term is RETIRED for' // &
    ' Fe/Zn/Cu/Mn; mass flow (unchanged) + diffusion to the root surface (Barber-Cushman, reusing rec 139''s own' // &
    ' RLD/SRL and the SAME declared f_l(theta) pair) + (Zn only) a bounded HATS ceiling (Hacisalihoglu 2001 Table II,' // &
    ' capped at S_mf+S_diff, never a field rate) now carry Fe/Zn/Cu/Mg/Mn uptake; Ca xylem route unchanged.' // &
    ' SHARED/mineral_module.f90.'
END IF
IF (tiller_survival_mech >= 1) THEN
  ! Record 99 dependencies: the cohort densities (record 66), the per-axis blade
  ! areas the demand is built from (record 65), the rank chronology those stand
  ! on (record 64) and the daily gain the supply is shared out of (record 69,
  ! the only writer of that quantity).
  IF (canopy_leaf_ranks /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 99: requires rec 64 = 1 (rank chronology).'
    ERROR STOP 'rec 99 requires rec 64 = 1'
  END IF
  IF (canopy_blade_area /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 99: requires rec 65 = 1 (per-axis blade areas: the demand basis).'
    ERROR STOP 'rec 99 requires rec 65 = 1'
  END IF
  IF (canopy_tillering /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 99: requires rec 66 = 1 (cohort densities).'
    ERROR STOP 'rec 99 requires rec 66 = 1'
  END IF
  IF (canopy_sink_lai /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 99: requires rec 69 = 1 (csl_dm_gain_day, the supply shared out).'
    ERROR STOP 'rec 99 requires rec 69 = 1'
  END IF
  IF (tiller_survival_mech == 4) THEN
    WRITE (*,'(A)') ' rec 99 tiller_survival_mech = 4 : tiller carbon affordability, deaths capped day by day by Abichou' // &
      ' Eqn T4.4 F(tau_n) (theta 3.32, delta 1.77) x U(t), the cull target F(tau_n)*U locked by the monotone MAX;' // &
      ' endpoint at flag-leaf ligulation. SHARED.'
  END IF
END IF
! Record 100 dependencies: the Haun-stage chronology the elongation trajectory
! is written in (record 64) and the preliminary final leaf number the rank
! profiles are indexed on (record 65). Both are diagnostic-only records.
IF (stem_demand_mech >= 1) THEN
  IF (canopy_leaf_ranks /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 100: requires rec 64 = 1 (Haun-stage chronology clr_HS / clr_FLN).'
    ERROR STOP 'rec 100 requires rec 64 = 1'
  END IF
  IF (canopy_blade_area /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 100: requires rec 65 = 1 (cba_Nl_fit, the provisional final leaf number).'
    ERROR STOP 'rec 100 requires rec 65 = 1'
  END IF
  sdm_dens = plant_density
  ! Record 100 mode 4 dependencies and member selection. The demand is per AXIS
  ! CLASS and its floral-transition proxy is the axis spikelet-marker clock
  ! (record 66); the bracket members are named by records 103 and 104, and 0 on
  ! either is refused, since no member is privileged.
  IF (stem_demand_mech == 4) THEN
    IF (canopy_tillering /= 1) THEN
      WRITE (*,'(A)') ' ERROR rec 100 = 4: requires rec 66 = 1 (ct_dens, ct_tau_birth, ct_tau_sm).'
      ERROR STOP 'rec 100 = 4 requires rec 66 = 1'
    END IF
    IF (sdm4_dur_bracket < 1) THEN
      WRITE (*,'(A)') ' ERROR rec 103: must name a bracket member when rec 100 = 4' // &
        ' (1/2/3 = onset F with Delta_I lo/mid/hi, 4/5/6 = onset M with lo/mid/hi); 0 = not selected.'
      ERROR STOP 'rec 103 must name a bracket member when rec 100 = 4'
    END IF
    IF (sdm4_rate_bracket < 1) THEN
      WRITE (*,'(A)') ' ERROR rec 104: must name a bracket member when rec 100 = 4' // &
        ' (1 = lo 0.7 x 1.0, 2 = mid 0.7 x 1.7, 3 = hi 1.7 mg DM degCd-1 internode-1); 0 = not selected.'
      ERROR STOP 'rec 104 must name a bracket member when rec 100 = 4'
    END IF
  END IF
END IF
! rec 101 dependency: with rec 100 = 0 or 2 the stem refuses nothing, so
! the record would be inert -- exactly as the released configuration requires rec 55 = 1.
IF (surplus_fate_mech >= 1) THEN
  ! The widening covers mode 4, which refuses carbon by construction; mode 3 is
  ! deliberately not widened, because its refusal is an artefact of the unbounded
  ! height recharge and giving that a sourced destination would launder it.
  IF (stem_demand_mech /= 1 .AND. stem_demand_mech /= 4) THEN
    WRITE (*,'(A)') ' ERROR rec 101: requires rec 100 = 1 or 4 (there is no refused carbon otherwise).'
    ERROR STOP 'rec 101 requires rec 100 = 1 or 4'
  END IF
END IF
! Record 89 dependencies: the window is the penultimate-leaf event of the
! record-64 rank chronology; the per-axis share needs the record-65 blade areas
! and the record-66 cohort densities; and the partition is written into the
! record-69 sink split.
IF (grain_number_mech >= 1) THEN
  IF (canopy_leaf_ranks /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 89: requires rec 64 = 1 (rank chronology: the penultimate-leaf window anchor).'
    ERROR STOP 'rec 89 requires rec 64 = 1'
  END IF
  IF (canopy_blade_area /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 89: requires rec 65 = 1 (per-axis blade areas: the per-axis spike share).'
    ERROR STOP 'rec 89 requires rec 65 = 1'
  END IF
  IF (canopy_tillering /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 89: requires rec 66 = 1 (axis population: there are no axes to bear ears without it).'
    ERROR STOP 'rec 89 requires rec 66 = 1'
  END IF
  IF (canopy_sink_lai /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 89: requires rec 69 = 1 (the sink partition the Fischer spike share is held through).'
    ERROR STOP 'rec 89 requires rec 69 = 1'
  END IF
  SELECT CASE (grain_number_mech)
  CASE (1)
    WRITE (*,'(A)') ' rec 89 grain_number_mech = 1 (A+D): Fischer 1985 spike window + 0.44 partition;' // &
      ' grain number PREDICTED per axis (100 florets g-1, 0.78 kernels per floret). FE, storage_onset_dvs2' // &
      ' and storage_peak_frac are NOT READ.'
  END SELECT
  ! Record 89 at or above 1 prints that the anchored spike-onset node is retired.
  ! growth_module.f90 gates that node under record 89 = 0, so a deck that asks
  ! for both the mechanism and the anchored node is refused outright, as records
  ! 64, 65, 66 and 69 already are.
  IF (stem_spike_realloc_delta > 0.0_wp) THEN
    WRITE (*,'(A,G14.6)') ' ERROR rec 89: storage_onset_dvs2 / rec 61 stem_spike_realloc_delta is RETIRED' // &
      ' under rec 89 >= 1, but rec 61 was declared > 0: ', stem_spike_realloc_delta
    ERROR STOP 'rec 89 with rec 61 > 0: the retired pre-anthesis stem->spike node cannot also be live'
  END IF
END IF
! rec 90 dependency and banner.
IF (floret_survival_mech >= 1) THEN
  IF (grain_number_mech /= 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 90: the floret-site chain requires rec 89 = 1; rec 89 read ', grain_number_mech
    ERROR STOP 'rec 90 requires rec 89 = 1'
  END IF
  SELECT CASE (floret_survival_mech)
  CASE (4)
    WRITE (*,'(A)') ' rec 90 floret_survival_mech = 4: floret-site chain; ear-bearing by SPIKE AFFORDABILITY' // &
      ' phi_sp = MIN(1, W_axis / (bgn_r_spike x window degCd)) per class, continuous, no bar; earless axes standing' // &
      ' (Fischer 1985 p.452); pool returned to ear-bearers. SHARED.'
  END SELECT
END IF
! rec 91 dependency and banner.
IF (grain_fill_mech >= 1) THEN
  IF (grain_number_mech /= 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 91: the grain-filling law requires rec 89 = 1 (there is no' // &
      ' PREDICTED kernel number to grow without it); rec 89 read ', grain_number_mech
    ERROR STOP 'rec 91 requires rec 89 = 1'
  END IF
  SELECT CASE (grain_fill_mech)
  CASE (1)
    WRITE (*,'(A)') ' rec 91 grain_fill_mech = 1: kernel growth RATE x filling DURATION.' // &
      ' Yin 2003 Ann Bot 91:361 Tab 2 beta rate (c_m 1.9283 mg kernel-1 d-1 at 17.5 C,' // &
      ' +0.11067 per K, t_m/t_e 0.4331); duration 405.581 degCd above 6.6187 C, DERIVED from' // &
      ' that table by Fischer 1985 Tab 6 p.455. TKW_potential is NOT READ.'
  END SELECT
END IF
! rec 92 dependency and banner.
IF (stem_wsc_mech >= 1) THEN
  IF (carbon_reserve_mode /= 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 92: the stem-WSC law requires rec 55 = 1 (there is no reserve' // &
      ' pool to cap, close or mobilise without it); rec 55 read ', carbon_reserve_mode
    ERROR STOP 'rec 92 requires rec 55 = 1'
  END IF
  SELECT CASE (stem_wsc_mech)
  CASE (1)
    WRITE (*,'(A)') ' rec 92 stem_wsc_mech = 1: stem WSC gets a CAPACITY (0.50 of stem DM),' // &
      ' net deposition ENDS 15 d after anthesis, and remobilisation is DEMAND-DRIVEN' // &
      ' (Schnyder 1993 New Phytol 123:233 pp.233/238/239). reserve_remob_rate is NOT READ.'
  END SELECT
  WRITE (*,'(A)') '   NOTE: overflow_to_reserve_total now reports the surplus LABELLED as WSC,' // &
    ' not the surplus moved to the stem. No carbon moves either way.'
END IF
! rec 105 dependencies and banner (POST-ANTHESIS-ONLY ARM, note sec.10;
! SHARED). rec 55 = 1 is implied by the rec 92 gate just above.
IF (wsc_reserve_law == 1) THEN
  IF (stem_wsc_mech /= 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 105: the bounded WSC draw requires rec 92 = 1 (the deficit cue, capacity' // &
      ' and deposition close it builds on); rec 92 read ', stem_wsc_mech
    ERROR STOP 'rec 105 requires rec 92 = 1'
  END IF
  IF (grain_fill_mech < 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 105: requires rec 91 >= 1 (item B''s filling clock gfb_tt); rec 91 read ', &
      grain_fill_mech
    ERROR STOP 'rec 105 requires rec 91 >= 1'
  END IF
  IF (wsc_law_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 106: must name a member (1 item-B clock / 2 legacy clock / 3 no ceiling)' // &
      ' when rec 105 = 1.'
    ERROR STOP 'rec 106 must name a member when rec 105 = 1'
  END IF
END IF
! Record 105 mode 2 dependencies and banner: the deficit-independent window
! release with a deficit-decided split. Record 55 = 1 is implied by the
! record-92 gate above. Four consistency stops, all pre-registered before the
! code was written.
IF (wsc_reserve_law == 2) THEN
  IF (stem_wsc_mech /= 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 105: mode 2 requires rec 92 = 1 (the deficit cue room, the 0.50 capacity' // &
      ' cap, the deposition close and the stem-WSC closure all live under it); rec 92 read ', stem_wsc_mech
    ERROR STOP 'rec 105 = 2 requires rec 92 = 1'
  END IF
  IF (grain_fill_mech < 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 105: mode 2 requires rec 91 >= 1 (the item-B filling clock gfb_tt is the' // &
      ' member-1 release clock and is reported on BOTH members); rec 91 read ', grain_fill_mech
    ERROR STOP 'rec 105 = 2 requires rec 91 >= 1'
  END IF
  IF (wsc_rel_window_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 107: must name a release-window member (1 = item-B filling clock /' // &
      ' 2 = legacy rec-59 clock) when rec 105 = 2. No member is privileged.'
    ERROR STOP 'rec 107 must name a member when rec 105 = 2'
  END IF
  IF (wsc_law_member /= 0) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 106: is the MODE-1 ceiling clock and is NOT READ under rec 105 = 2;' // &
      ' it must be 0 so no mode-1 member code is silently reused; rec 106 read ', wsc_law_member
    ERROR STOP 'rec 106 must be 0 when rec 105 = 2'
  END IF
  WRITE (*,'(A,I0,A)') ' rec 105 wsc_reserve_law = 2 (rec 107 release-window member ', wsc_rel_window_member, &
    ') : THE REPAIR -- the stem WSC pool is RELEASED over the window REGARDLESS of the grain deficit' // &
    ' (rel = MIN(W_reserve, W_stem, cap), cap = W_reserve x dtau/tau_rem), and the deficit decides ONLY' // &
    ' the grain share: CASE 1 rel_g = MIN(rel, room/eff) -> grain at the rec 60 Gebbing efficiency;' // &
    ' CASE 2 rel_f = rel - rel_g -> the NAMED ledger cres_respired AT UNITY (no efficiency).' // &
    ' Members: 1 = item-B filling clock (gfb_dfill - gfb_tt) / 2 = legacy rec-59 clock' // &
    ' (reserve_end_dvs2 - dvs2); SWEPT, never collapsed. rec 106 NOT READ. rec 100 = 0 legal.' // &
    ' SHARED.' // &
    ' Off by default.'
END IF
! Records 109 to 113 interlocks and the one echo line. The echo is emitted only
! when the law is on.
IF (height_internode_law == 1) THEN
  IF (height_reserve_basis /= 0) THEN
    WRITE (*,'(A)') ' ERROR rec 109: the internode height law requires rec 108 height_reserve_basis = 0' // &
      ' -- two height laws may not run together (design note sec.5).'
    ERROR STOP 'rec 109 requires rec 108 = 0'
  END IF
  IF (stem_demand_mech /= 0) THEN
    WRITE (*,'(A)') ' ERROR rec 109: the internode height law requires rec 100 stem_demand_mech = 0' // &
      ' -- two stem laws may not run together (design note sec.5).'
    ERROR STOP 'rec 109 requires rec 100 = 0'
  END IF
  IF (canopy_leaf_ranks /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 109: the internode height law requires rec 64 canopy_leaf_ranks = 1' // &
      ' -- E1/E2 read that per-culm ligulation chronology (design note sec.6/sec.10).'
    ERROR STOP 'rec 109 requires rec 64 = 1'
  END IF
  IF (canopy_tillering /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 109: the internode height law requires rec 66 canopy_tillering = 1' // &
      ' -- E5'' converts the ground-area structural increment to a per-axis one through ct_axes_live.'
    ERROR STOP 'rec 109 requires rec 66 = 1'
  END IF
  IF (hgt_int_window_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 110: must name a window member (1-3 = Fournier collar onset x Delta_I' // &
      ' 277/360/476 degCd; 4-6 = Malvoisin i-2 onset x the same three). No member is privileged.'
    ERROR STOP 'rec 110 must name a member when rec 109 = 1'
  END IF
  IF (hgt_int_length_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 111: must name a length-potential member (1 = Kirby folio 133 glasshouse;' // &
      ' 2/3/4 = Rivera-Amado folio 161 mean/min/max; 5 = cse_li_irank Abichou, no peduncle).'
    ERROR STOP 'rec 111 must name a member when rec 109 = 1'
  END IF
  IF (hgt_int_lineic_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 112: must name a lineic member (1/2/3 = Rivera folio 161 specific weight' // &
      ' mean/min/max, total DM; 4 = Gauthier Table S1.4 SLM_IN, structural).'
    ERROR STOP 'rec 112 must name a member when rec 109 = 1'
  END IF
  IF (hgt_int_ear_member == 0) THEN
    WRITE (*,'(A)') ' ERROR rec 113: must name an ear member (1 = L_ear 0 stem-only; 2 = L_ear 10.1 cm with' // &
      ' peduncle end Delta_I; 3 = L_ear 10.1 cm with peduncle end cse_ep_win*P).'
    ERROR STOP 'rec 113 must name a member when rec 109 = 1'
  END IF
  WRITE (*,'(A,4(A,I0),A)') &
    ' rec 109 height_internode_law = 1 : canopy height from INTERNODE ELONGATION (E1-E8/E5'',', &
    ' window member ', hgt_int_window_member, ', length member ', hgt_int_length_member, &
    ', lineic member ', hgt_int_lineic_member, ', ear member ', hgt_int_ear_member, &
    '). The model''s OWN structural stem increment constrains realised length;' // &
    ' the carbon shortfall is a MODEL-INTERNAL diagnostic that licenses NO parameter change.' // &
    ' Height is reported on THREE bases (B1 collar / B2 peduncle-included / B3 spike+awns)' // &
    ' because the WDF25 height basis is UNKNOWN. NOT accepted, NOT audited, NOT the default.'
END IF
! Records 114, 115, 116, 117 and 118: the one echo line (or block) each,
! emitted only when the record is on.
IF (kernel_mass_ceiling == 1) THEN
  WRITE (*,'(A,G14.6,A)') ' rec 118 kernel_mass_ceiling = 1 : per-kernel mass bounded by the MEASURED cultivar' // &
    ' TKW_potential = ', TKW_potential, ' g per 1000 kernels (input/cultivar_yield.csv row 2), inside the rec-91' // &
    ' rate x duration law (Yin et al. 2003 Ann Bot 91:361 Table 2 w_max; Sofield et al. 1977 Aust J Plant' // &
    ' Physiol 4:785).  SHARED.f90.  v1.0.0 fix A, approved 2026-10-05.'
END IF
! rec 94 dependency and banner.
IF (canopy_n_flag_zone >= 1) THEN
  ! Refused rather than defaulted: the zone depth of record 94 is the model's own
  ! flag-rank green area and does not exist without the organ chain. No fallback
  ! depth is invented.
  IF (canopy_leaf_ranks /= 1 .OR. canopy_blade_area /= 1 .OR. canopy_tillering /= 1 .OR. &
      canopy_leaf_senesc /= 1 .OR. canopy_sheath_ear /= 1 .OR. canopy_sink_lai /= 1 .OR. &
      lai_from_organs /= 1) THEN
    WRITE (*,'(A)') ' ERROR rec 94: the flag-zone N profile needs the organ chain that DEFINES' // &
      ' its depth; rec 64, 65, 66, 67, 68, 69 and 84 must ALL be 1. No fallback depth exists.'
    ERROR STOP 'rec 94 requires rec 64-69 and rec 84 = 1'
  END IF
  SELECT CASE (canopy_n_flag_zone)
  CASE (3)
    WRITE (*,'(A)') ' rec 94 canopy_n_flag_zone = 3: THE ITEM. Mode 1 plus the layered f_PAR' // &
      ' of FPB Eqn 17 p.785 inside item F''s Eqn-6 limb (audit note (iii)).'
  END SELECT
END IF
! rec 86 dependency and banner. SHARED.f90
gnf_critN_basis = (grainN_sink_mech == 1 .OR. grainN_sink_mech == 3 .OR. grainN_sink_mech == 5)
! FPB Eqn 6 (and the Eqn-1 (1+r) pool basis) are ON for every live mode EXCEPT the
! mode-5 attribution arm, which reproduces the rejected arm. SHARED.f90.
gnf_eqn6 = (grainN_sink_mech == 1 .OR. grainN_sink_mech == 2 .OR. grainN_sink_mech == 4)
IF (grainN_sink_mech >= 1) THEN
  IF (grain_number_mech < 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 86: the grain-N sink law is on a PER-EAR basis and requires' // &
      ' rec 89 >= 1 (without item A chain gnx_ears_m2 is 0, the basis does not exist, and the' // &
      ' grain would silently receive no N); rec 89 read ', grain_number_mech
    ERROR STOP 'rec 86 requires rec 89 >= 1'
  END IF
  SELECT CASE (grainN_sink_mech)
  CASE (1)
    WRITE (*,'(A)') ' rec 86 grainN_sink_mech = 1: THE ITEM. Grain N = min(potential per-ear sink,' // &
      ' available N) -- Bertheloot, Martre & Andrieu 2008 FPB 35:781 Eqn 11-12 p.785, Table 1 p.784' // &
      ' (N_g(0) = 2.4 mg ear-1, g = 0.005 degCd-1 DECLARED DERIVED NOT MEASURED, tau = 250 degCd,' // &
      ' thermal time base 0 C). Supply = assimilation + mobile pool + delta-limited organ draw above' // &
      ' the p.1711 structural floors. Post-anthesis UPTAKE demand read off the SAME potential' // &
      ' (the Jamieson & Semenov NHI trajectory is retired). Critical-N shoot mass on W_grain+W_chaff.'
    WRITE (*,'(A)') '   EQN-6 REPAIR 2026-09-21: lamina photosynthetic N now runs FPB Eqn 4 in FULL --' // &
      ' Eqn 5 degradation (delta = 0.008 [degCd]-1) AND Eqn 6 synthesis (sigma = 0.0015 [degCd]-1,' // &
      ' k1 = 0.0018 kg kg-1, k2 = 10 J m-2 s-1, all FPB Table 1 p.784), on the Eqn-1 p.783 pool basis' // &
      ' N^ph = (N^tot - N^struct)/(1+r) with r = 0.10 (p.786). FPB p.786 fitted sigma, delta, k1 and k2' // &
      ' JOINTLY across Eqns 5 AND 6, so delta is only sourced when sigma is present.'
  END SELECT
END IF
! rec 85 refusals: one law per resistance, and Eq. 4.78 is 0/0 at alpha = 0.
IF (rca_teh_full == 1 .AND. rca_wind_attenuation == 1) THEN
  WRITE (*,*) 'REFUSING TO RUN: rec 85 (rca_teh_full) = 1 AND rec 70 (rca_wind_attenuation) = 1.'
  WRITE (*,*) '  Both select a canopy boundary-layer law; set rec 70 = 0 to run Teh Eq. 4.78 in full.'
  ERROR STOP 'rec 85 and rec 70 both ON'
END IF
IF (rca_teh_full == 1 .AND. a_wind <= 0.0_wp) THEN
  WRITE (*,*) 'REFUSING TO RUN: rec 85 (rca_teh_full) = 1 with alpha (rec 12) <= 0:', a_wind
  WRITE (*,*) '  Teh Eq. 4.78 p.97 is 0/0 at alpha = 0; Table 4.2 p.93 gives wheat alpha = 2.5.'
  ERROR STOP 'rec 85 with alpha <= 0'
END IF


! Optional command-line override of the energy-balance mode, without mutating
! the driver. A bare run keeps the driver value; one integer argument forces
! that mode reproducibly. The override is validated by the same range guards as
! the driver value.
IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
  BLOCK
    CHARACTER(LEN=32) :: eb_arg
    INTEGER :: eb_cli_val, ios_cli
    CALL GET_COMMAND_ARGUMENT(1, eb_arg)
    READ (eb_arg, *, IOSTAT=ios_cli) eb_cli_val
    IF (ios_cli == 0) THEN
      energy_balance_mode = eb_cli_val
      eb_mode_from_cli = .TRUE.
    ELSE
      WRITE (*,*) 'WARNING: could not parse energy_balance_mode CLI arg "', &
        TRIM(eb_arg), '"; using driver_main.txt value.'
    END IF
  END BLOCK
END IF

! Mode implemented: 2 = the three-node canopy, soil and airspace energy balance
! with the prognostic one-dimensional soil-heat profile supplying the ground heat
! flux, through energy_balance_module.f90. The retired Shuttleworth-Wallace mode
! and its source file were removed from v1.0; record 39 accepts 2 only.
IF (energy_balance_mode < 0 .OR. energy_balance_mode > 3) THEN
  WRITE (*,*) 'ERROR: energy_balance_mode must be 0-3. Got:', energy_balance_mode
  STOP 1
END IF
IF (energy_balance_mode >= 3) THEN
  WRITE (*,*) '------------------------------------------------------------'
  WRITE (*,*) 'ERROR: energy_balance_mode =', energy_balance_mode, 'is RESERVED and NOT implemented.'
  WRITE (*,*) 'Implemented: 2 (three-node EB + 1-D soil heat).'
  WRITE (*,*) 'Mode 3 (full CO2ST attribution) arrives in a later stage.'
  WRITE (*,*) '------------------------------------------------------------'
  STOP 1
END IF
! Mode 2 marches a PROGNOSTIC soil-temperature profile, which requires
! sub-steps that tile real clock time, so it requires the full 24-hour
! integration (record 27 = 2); the daylight-only quadrature of record 27 = 1
! does not tile the day.
IF (energy_balance_mode == 2 .AND. integration_mode /= 2) THEN
  WRITE (*,*) 'ERROR: energy_balance_mode=2 (1-D soil heat) requires integration_mode=2 (full 24-hr).'
  WRITE (*,*) '       Set record 27 (integration mode) = 2 in input/driver_main.txt.'
  STOP 1
END IF

! **************************************************************************** !
! RELEASED-OPTION VALIDATION
!
! The release contains the code of the process options it uses and of the
! alternative formulations of Part 1 Table 6, which are kept and documented as
! built and inactive. A value whose code is not in the release must therefore
! be REFUSED here: if it were accepted, the absent branch would let the run
! fall through to a different mechanism, or to an undefined value.
!
! The allowed set, per restricted record, is the value or values used by the 42
! released decks plus that record's own Part 1 Table 6 retained alternatives.
!
! Restricted here: the 123 mode and switch records. NOT restricted: the 51
! plain-parameter records, every site, soil, weather, management and cultivar
! input, every physical parameter, record 88 (the sowing day), the reserved
! records 6, 119, 120, 170 and 182, whose own must-be-zero guards stay, and the
! mode records with a single implemented branch (25, 30, 32, 76, 79, 80).
!
! This block runs after all 183 reads, after the energy-balance command-line
! override and after the record-combination guards, so that a refused mode
! cannot enter through the command line either.
! **************************************************************************** !
IF (water_stress_mode /= 1) &
  CALL reject_retired_option (23, 'water_stress_mode', water_stress_mode, '1')
IF (.NOT. use_soil_water_stress) CALL reject_retired_option (24, 'use_soil_water_stress', 0, '1')
IF (use_explicit_infiltration) CALL reject_retired_option (26, 'use_explicit_infiltration', 1, '0')
IF (integration_mode /= 2) &
  CALL reject_retired_option (27, 'integration_mode', integration_mode, '2')
IF (meteo_debug) CALL reject_retired_option (31, 'meteo_debug', 1, '0')
IF (k_avg_mode /= 0) &
  CALL reject_retired_option (35, 'k_avg_mode', k_avg_mode, '0')
IF (hyd_model /= 0) &
  CALL reject_retired_option (36, 'hyd_model', hyd_model, '0')
IF (init_profile_mode /= 0 .AND. init_profile_mode /= 2) &
  CALL reject_retired_option (37, 'init_profile_mode', init_profile_mode, '0, 2')
IF (energy_balance_mode /= 2) &
  CALL reject_retired_option (39, 'energy_balance_mode', energy_balance_mode, '2')
IF (eb_freeconv_fallback /= 0) &
  CALL reject_retired_option (40, 'eb_freeconv_fallback', eb_freeconv_fallback, '0')
IF (soil_T_init_method /= 1) &
  CALL reject_retired_option (50, 'soil_T_init_method', soil_T_init_method, '1')
IF (residue_evap_mode /= 1) &
  CALL reject_retired_option (52, 'residue_evap_mode', residue_evap_mode, '1')
IF (carbon_reserve_mode /= 1) &
  CALL reject_retired_option (55, 'carbon_reserve_mode', carbon_reserve_mode, '1')
IF (leafN_drawdown /= 0) &
  CALL reject_retired_option (62, 'leafN_drawdown', leafN_drawdown, '0')
IF (grainN_partition /= 0) &
  CALL reject_retired_option (63, 'grainN_partition', grainN_partition, '0')
IF (canopy_leaf_ranks /= 1) &
  CALL reject_retired_option (64, 'canopy_leaf_ranks', canopy_leaf_ranks, '1')
IF (canopy_blade_area /= 1) &
  CALL reject_retired_option (65, 'canopy_blade_area', canopy_blade_area, '1')
IF (canopy_tillering /= 1) &
  CALL reject_retired_option (66, 'canopy_tillering', canopy_tillering, '1')
IF (canopy_leaf_senesc /= 1) &
  CALL reject_retired_option (67, 'canopy_leaf_senesc', canopy_leaf_senesc, '1')
IF (canopy_sheath_ear /= 1) &
  CALL reject_retired_option (68, 'canopy_sheath_ear', canopy_sheath_ear, '1')
IF (canopy_sink_lai /= 1) &
  CALL reject_retired_option (69, 'canopy_sink_lai', canopy_sink_lai, '1')
IF (rca_wind_attenuation /= 1) &
  CALL reject_retired_option (70, 'rca_wind_attenuation', rca_wind_attenuation, '1')
IF (wind_hourly_mode /= 0) &
  CALL reject_retired_option (71, 'wind_hourly_mode', wind_hourly_mode, '0')
IF (rss_theta_s_consistency /= 1) &
  CALL reject_retired_option (75, 'rss_theta_s_consistency', rss_theta_s_consistency, '1')
IF (ci_solve_mode /= 1) &
  CALL reject_retired_option (77, 'ci_solve_mode', ci_solve_mode, '1')
IF (ft_pot_unstressed /= 1) &
  CALL reject_retired_option (78, 'ft_pot_unstressed', ft_pot_unstressed, '1')
IF (tiller_carbon_gate /= 1) &
  CALL reject_retired_option (82, 'tiller_carbon_gate', tiller_carbon_gate, '1')
IF (leaf_expansion_water /= 1) &
  CALL reject_retired_option (83, 'leaf_expansion_water', leaf_expansion_water, '1')
IF (lai_from_organs /= 1) &
  CALL reject_retired_option (84, 'lai_from_organs', lai_from_organs, '1')
IF (rca_teh_full /= 0) &
  CALL reject_retired_option (85, 'rca_teh_full', rca_teh_full, '0')
IF (grainN_sink_mech /= 1) &
  CALL reject_retired_option (86, 'grainN_sink_mech', grainN_sink_mech, '1')
IF (leafN_cap_route /= 1) &
  CALL reject_retired_option (87, 'leafN_cap_route', leafN_cap_route, '1')
IF (grain_number_mech /= 1) &
  CALL reject_retired_option (89, 'grain_number_mech', grain_number_mech, '1')
IF (canopy_n_base_profile /= 1) &
  CALL reject_retired_option (76, 'canopy_n_base_profile', canopy_n_base_profile, '1')
IF (floret_survival_mech /= 4) &
  CALL reject_retired_option (90, 'floret_survival_mech', floret_survival_mech, '4')
IF (grain_fill_mech /= 1) &
  CALL reject_retired_option (91, 'grain_fill_mech', grain_fill_mech, '1')
IF (stem_wsc_mech /= 1) &
  CALL reject_retired_option (92, 'stem_wsc_mech', stem_wsc_mech, '1')
IF (soil_n_mech /= 0) &
  CALL reject_retired_option (93, 'soil_n_mech', soil_n_mech, '0')
IF (canopy_n_flag_zone /= 3) &
  CALL reject_retired_option (94, 'canopy_n_flag_zone', canopy_n_flag_zone, '3')
IF (tleaf_class_diag /= 0) &
  CALL reject_retired_option (95, 'tleaf_class_diag', tleaf_class_diag, '0')
IF (sw_radiosity_mode /= 0) &
  CALL reject_retired_option (96, 'sw_radiosity_mode', sw_radiosity_mode, '0')
IF (itemF_Nclose_mode /= 1) &
  CALL reject_retired_option (97, 'itemF_Nclose_mode', itemF_Nclose_mode, '1')
IF (spike_sink_mech /= 0) &
  CALL reject_retired_option (98, 'spike_sink_mech', spike_sink_mech, '0')
IF (tiller_survival_mech /= 4) &
  CALL reject_retired_option (99, 'tiller_survival_mech', tiller_survival_mech, '4')
IF (stem_demand_mech /= 0) &
  CALL reject_retired_option (100, 'stem_demand_mech', stem_demand_mech, '0')
IF (surplus_fate_mech /= 0) &
  CALL reject_retired_option (101, 'surplus_fate_mech', surplus_fate_mech, '0')
IF (stem_sheath_sep_mode /= 0) &
  CALL reject_retired_option (102, 'stem_sheath_sep_mode', stem_sheath_sep_mode, '0')
IF (sdm4_dur_bracket /= 0) &
  CALL reject_retired_option (103, 'sdm4_dur_bracket', sdm4_dur_bracket, '0')
IF (sdm4_rate_bracket /= 0) &
  CALL reject_retired_option (104, 'sdm4_rate_bracket', sdm4_rate_bracket, '0')
IF (wsc_reserve_law /= 2) &
  CALL reject_retired_option (105, 'wsc_reserve_law', wsc_reserve_law, '2')
IF (wsc_law_member /= 0) &
  CALL reject_retired_option (106, 'wsc_law_member', wsc_law_member, '0')
IF (wsc_rel_window_member /= 2) &
  CALL reject_retired_option (107, 'wsc_rel_window_member', wsc_rel_window_member, '2')
IF (height_reserve_basis /= 0) &
  CALL reject_retired_option (108, 'height_reserve_basis', height_reserve_basis, '0')
IF (height_internode_law /= 0 .AND. height_internode_law /= 1) &
  CALL reject_retired_option (109, 'height_internode_law', height_internode_law, '0, 1')
IF (hgt_int_window_member /= 0 .AND. hgt_int_window_member /= 1 .AND. hgt_int_window_member /= 2 .AND. &
  hgt_int_window_member /= 3 .AND. hgt_int_window_member /= 4 .AND. hgt_int_window_member /= 5 .AND. &
  hgt_int_window_member /= 6) &
  CALL reject_retired_option (110, 'hgt_int_window_member', hgt_int_window_member, '0, 1, 2, 3, 4, 5, 6')
IF (hgt_int_length_member /= 0 .AND. hgt_int_length_member /= 1 .AND. hgt_int_length_member /= 2 .AND. &
  hgt_int_length_member /= 3 .AND. hgt_int_length_member /= 4 .AND. hgt_int_length_member /= 5) &
  CALL reject_retired_option (111, 'hgt_int_length_member', hgt_int_length_member, '0, 1, 2, 3, 4, 5')
IF (hgt_int_lineic_member /= 0 .AND. hgt_int_lineic_member /= 1 .AND. hgt_int_lineic_member /= 2 .AND. &
  hgt_int_lineic_member /= 3 .AND. hgt_int_lineic_member /= 4) &
  CALL reject_retired_option (112, 'hgt_int_lineic_member', hgt_int_lineic_member, '0, 1, 2, 3, 4')
IF (hgt_int_ear_member /= 0 .AND. hgt_int_ear_member /= 1 .AND. hgt_int_ear_member /= 2 .AND. &
  hgt_int_ear_member /= 3) &
  CALL reject_retired_option (113, 'hgt_int_ear_member', hgt_int_ear_member, '0, 1, 2, 3')
IF (mineral_share_basis /= 0) &
  CALL reject_retired_option (114, 'mineral_share_basis', mineral_share_basis, '0')
IF (gnf_ng0_basis /= 0) &
  CALL reject_retired_option (115, 'gnf_ng0_basis', gnf_ng0_basis, '0')
IF (limb_census_mode /= 0) &
  CALL reject_retired_option (116, 'limb_census_mode', limb_census_mode, '0')
IF (chaff_growth_mode /= 0) &
  CALL reject_retired_option (117, 'chaff_growth_mode', chaff_growth_mode, '0')
! Record 118 keeps 0 in its allowed set as well: the guard above refuses record
! 125 modes 2 and 3 unless record 118 = 0, and those two modes are retained
! Table 6 alternatives, so forbidding 0 would make them unreachable.
IF (kernel_mass_ceiling /= 0 .AND. kernel_mass_ceiling /= 1) &
  CALL reject_retired_option (118, 'kernel_mass_ceiling', kernel_mass_ceiling, '0, 1')
IF (n_payment_gate /= 2) &
  CALL reject_retired_option (121, 'n_payment_gate', n_payment_gate, '2')
IF (lai_n_arrest /= 2) &
  CALL reject_retired_option (122, 'lai_n_arrest', lai_n_arrest, '2')
IF (tiller_n_afford /= 1) &
  CALL reject_retired_option (123, 'tiller_n_afford', tiller_n_afford, '1')
IF (root_front_law /= 2) &
  CALL reject_retired_option (124, 'root_front_law', root_front_law, '2')
IF (kernel_potential_law /= 1 .AND. kernel_potential_law /= 2 .AND. kernel_potential_law /= 3) &
  CALL reject_retired_option (125, 'kernel_potential_law', kernel_potential_law, '1, 2, 3')
IF (mineral_grain_sink_cap /= 1) &
  CALL reject_retired_option (126, 'mineral_grain_sink_cap', mineral_grain_sink_cap, '1')
IF (repro_surplus_fate /= 1) &
  CALL reject_retired_option (127, 'repro_surplus_fate', repro_surplus_fate, '1')
IF (n_construction_status /= 2) &
  CALL reject_retired_option (128, 'n_construction_status', n_construction_status, '2')
IF (sink_full_carbon_fate /= 1) &
  CALL reject_retired_option (129, 'sink_full_carbon_fate', sink_full_carbon_fate, '1')
IF (grain_n_min_conc /= 1 .AND. grain_n_min_conc /= 2 .AND. grain_n_min_conc /= 3) &
  CALL reject_retired_option (130, 'grain_n_min_conc', grain_n_min_conc, '1, 2, 3')
IF (mineral_supply_law /= 1) &
  CALL reject_retired_option (131, 'mineral_supply_law', mineral_supply_law, '1')
IF (leafn_hold_law /= 0) &
  CALL reject_retired_option (132, 'leafn_hold_law', leafn_hold_law, '0')
IF (leafn_turnover_law /= 0) &
  CALL reject_retired_option (133, 'leafn_turnover_law', leafn_turnover_law, '0')
IF (n_afford_priority_law /= 0) &
  CALL reject_retired_option (134, 'n_afford_priority_law', n_afford_priority_law, '0')
IF (n_payment_luxury_pass2 /= 0) &
  CALL reject_retired_option (135, 'n_payment_luxury_pass2', n_payment_luxury_pass2, '0')
IF (n_afford_standing_basis /= 0) &
  CALL reject_retired_option (136, 'n_afford_standing_basis', n_afford_standing_basis, '0')
IF (floret_n_law /= 0) &
  CALL reject_retired_option (137, 'floret_n_law', floret_n_law, '0')
IF (spike_n_priority /= 0) &
  CALL reject_retired_option (138, 'spike_n_priority', spike_n_priority, '0')
IF (n_uptake_diffusion /= 1) &
  CALL reject_retired_option (139, 'n_uptake_diffusion', n_uptake_diffusion, '1')
IF (mineral_supply_mech /= 1) &
  CALL reject_retired_option (140, 'mineral_supply_mech', mineral_supply_mech, '1')
IF (lamina_n_remob_pre /= 1) &
  CALL reject_retired_option (141, 'lamina_n_remob_pre', lamina_n_remob_pre, '1')
IF (lamina_n_expansion /= 0) &
  CALL reject_retired_option (142, 'lamina_n_expansion', lamina_n_expansion, '0')
IF (n_status_construction_price /= 0 .AND. n_status_construction_price /= 1) &
  CALL reject_retired_option (143, 'n_status_construction_price', n_status_construction_price, '0, 1')
IF (carbon_coupled_leaf_expansion /= 0 .AND. carbon_coupled_leaf_expansion /= 1) &
  CALL reject_retired_option (144, 'carbon_coupled_leaf_expansion', carbon_coupled_leaf_expansion, '0, 1')
IF (leaf_maintenance_basis /= 1) &
  CALL reject_retired_option (145, 'leaf_maintenance_basis', leaf_maintenance_basis, '1')
IF (par_reflection_basis /= 1) &
  CALL reject_retired_option (146, 'par_reflection_basis', par_reflection_basis, '1')
IF (lamina_n_area_retention /= 0) &
  CALL reject_retired_option (147, 'lamina_n_area_retention', lamina_n_area_retention, '0')
IF (wsc_fate_law /= 1) &
  CALL reject_retired_option (148, 'wsc_fate_law', wsc_fate_law, '1')
IF (wb_surface_limit /= 3) &
  CALL reject_retired_option (149, 'wb_surface_limit', wb_surface_limit, '3')
IF (pheno_dev_factor_window /= 1) &
  CALL reject_retired_option (150, 'pheno_dev_factor_window', pheno_dev_factor_window, '1')
IF (fe_trait_mode /= 1) &
  CALL reject_retired_option (151, 'fe_trait_mode', fe_trait_mode, '1')
IF (census_files /= 1) &
  CALL reject_retired_option (152, 'census_files', census_files, '1')
IF (n_init_placement /= 1) &
  CALL reject_retired_option (153, 'n_init_placement', n_init_placement, '1')
IF (organ_n_turnover_nema /= 1) &
  CALL reject_retired_option (154, 'organ_n_turnover_nema', organ_n_turnover_nema, '1')
IF (n_uptake_regulation /= 1) &
  CALL reject_retired_option (155, 'n_uptake_regulation', n_uptake_regulation, '1')
IF (n_construction_ceiling_curve /= 1) &
  CALL reject_retired_option (156, 'n_construction_ceiling_curve', n_construction_ceiling_curve, '1')
IF (lamina_n_storage /= 3) &
  CALL reject_retired_option (157, 'lamina_n_storage', lamina_n_storage, '3')
IF (nitrate_drain_rate_basis /= 1) &
  CALL reject_retired_option (158, 'nitrate_drain_rate_basis', nitrate_drain_rate_basis, '1')
IF (mobile_n_nema_estimate /= 0) &
  CALL reject_retired_option (159, 'mobile_n_nema_estimate', mobile_n_nema_estimate, '0')
IF (justes_nmax_floor_basis /= 1) &
  CALL reject_retired_option (160, 'justes_nmax_floor_basis', justes_nmax_floor_basis, '1')
IF (canopy_clock_driver /= 0 .AND. canopy_clock_driver /= 1) &
  CALL reject_retired_option (161, 'canopy_clock_driver', canopy_clock_driver, '0, 1')
IF (canopy_clock_basis /= 1) &
  CALL reject_retired_option (162, 'canopy_clock_basis', canopy_clock_basis, '1')
IF (postanth_uptake_regulation /= 0 .AND. postanth_uptake_regulation /= 1) &
  CALL reject_retired_option (163, 'postanth_uptake_regulation', postanth_uptake_regulation, '0, 1')
IF (lamina_capacity_profile_basis /= 0 .AND. lamina_capacity_profile_basis /= 1) &
  CALL reject_retired_option (164, 'lamina_capacity_profile_basis', lamina_capacity_profile_basis, '0, 1')
IF (wsc_release_delay_basis /= 1) &
  CALL reject_retired_option (165, 'wsc_release_delay_basis', wsc_release_delay_basis, '1')
IF (grain_ca_xylem_throughput /= 1) &
  CALL reject_retired_option (166, 'grain_ca_xylem_throughput', grain_ca_xylem_throughput, '1')
IF (lamina_synthesis_sheath_basis /= 2) &
  CALL reject_retired_option (167, 'lamina_synthesis_sheath_basis', lamina_synthesis_sheath_basis, '2')
IF (mineral_impedance_basis /= 1) &
  CALL reject_retired_option (168, 'mineral_impedance_basis', mineral_impedance_basis, '1')
IF (zn_uptake_kinetic_basis /= 1) &
  CALL reject_retired_option (169, 'zn_uptake_kinetic_basis', zn_uptake_kinetic_basis, '1')
IF (mineral_n_coupling_basis /= 1) &
  CALL reject_retired_option (171, 'mineral_n_coupling_basis', mineral_n_coupling_basis, '1')
IF (stem_reserve_deposition_n_basis /= 0) &
  CALL reject_retired_option (172, 'stem_reserve_deposition_n_basis', stem_reserve_deposition_n_basis, '0')
IF (spac_sink_coupling /= 1) &
  CALL reject_retired_option (173, 'spac_sink_coupling', spac_sink_coupling, '1')
IF (canopy_fln_prediction /= 0) &
  CALL reject_retired_option (174, 'canopy_fln_prediction', canopy_fln_prediction, '0')
IF (grain_n_partition_basis /= 4) &
  CALL reject_retired_option (175, 'grain_n_partition_basis', grain_n_partition_basis, '4')
IF (vegetative_n_demand_basis /= 0) &
  CALL reject_retired_option (176, 'vegetative_n_demand_basis', vegetative_n_demand_basis, '0')
IF (grain_n_potential_trait /= 0) &
  CALL reject_retired_option (177, 'grain_n_potential_trait', grain_n_potential_trait, '0')
IF (chaff_structural_n_basis /= 0 .AND. chaff_structural_n_basis /= 1 .AND. chaff_structural_n_basis /= 2) &
  CALL reject_retired_option (178, 'chaff_structural_n_basis', chaff_structural_n_basis, '0, 1, 2')
IF (lamina_expansion_n_limit /= 1) &
  CALL reject_retired_option (179, 'lamina_expansion_n_limit', lamina_expansion_n_limit, '1')
IF (grain_ca_loading_basis /= 1) &
  CALL reject_retired_option (180, 'grain_ca_loading_basis', grain_ca_loading_basis, '1')
IF (mg_solution_basis /= 0) &
  CALL reject_retired_option (181, 'mg_solution_basis', mg_solution_basis, '0')
IF (mineral_buffer_power_basis /= 1) &
  CALL reject_retired_option (183, 'mineral_buffer_power_basis', mineral_buffer_power_basis, '1')

! physics_mode is an internal incumbent alias, frozen at 0. It is NEVER set from
! energy_balance_mode: mode 1 uses energy_balance_module.f90, not the old TSEB.
physics_mode = 0

! Report the actual active configuration so the hidden-mode error cannot recur.
WRITE (*,*) '============================================================'
WRITE (*,*) ' Energy balance / hydrology-init configuration'
WRITE (*,*) '------------------------------------------------------------'
WRITE (*,*) ' init_profile_mode    =', init_profile_mode
WRITE (*,*) ' init_bottom_head_m   =', init_bottom_head_m
WRITE (*,*) ' energy_balance_mode  =', energy_balance_mode
IF (eb_mode_from_cli) THEN
  WRITE (*,*) ' EB mode source        = COMMAND-LINE override (driver_main.txt unchanged)'
ELSE
  WRITE (*,*) ' EB mode source        = driver_main.txt'
END IF
WRITE (*,*) ' active EB path        = three-node canopy-soil-airspace EB (scientific path)'
WRITE (*,*) ' three-node EB         = ACTIVE (component radiation, Stage 4)'
WRITE (*,*) ' soil heat mode        = 1-D prognostic profile (mode 2, WIRED: conductive G'
WRITE (*,*) '                         in the EB solve; profile advances on the accepted pass)'
WRITE (*,*) ' raa stability         = Monin-Obukhov (Paulson 1970) on above-canopy heat'
WRITE (*,*) ' free-conv fallback    = off (MO only)'
WRITE (*,*) ' physics_mode (alias)  =', physics_mode
WRITE (*,*) '============================================================'

! **************************************************************************** !
! 2. Maximum number of simulated days, from the days present in the input data.
! **************************************************************************** !
CALL count_lines(NDays)
NDays     = NDays - 1
WRITE (*,*) 'Days in input data ', NDays
! **************************************************************************** !
! 3. Allocate the daily arrays to that number of days.
! **************************************************************************** !
ALLOCATE (pheno_stage  (NDays), days_sim(NDays), lai_sim(NDays))
ALLOCATE (wg_gleaves(NDays), wg_dleaves(NDays), wg_Tleaves(NDays))
ALLOCATE (wg_stem(NDays), wg_roots(NDays), wg_storage(NDays))
ALLOCATE (plant_height(NDays), A_canopy(NDays), main_res(NDays), grow_res(NDays))
ALLOCATE (soil_le(NDays), plant_le(NDays), tot_le(NDays))
ALLOCATE (soil_h(NDays), plant_h(NDays), tot_h(NDays))
ALLOCATE (E_soil(NDays), E_plant(NDays), ET_a(NDays))

ALLOCATE (array_input_data (NDays, NinVar))
ALLOCATE (Year(NDays), Month(NDays), Day(NDays), s_hr(NDays))
ALLOCATE (Tmax_d(NDays), Tmin_d(NDays), RHmax_d(NDays), RHmin_d(NDays))
ALLOCATE (w_d(NDays), WS_d(NDays), Irr_d(NDays), Rad_d(NDays))

  ALLOCATE (A_data(NDays,5), flux_daily(NDays,9), flux_daily_full(NDays,59), swc_daily(NDays,4))
ALLOCATE (hydro_fail_days(NDays))
hydro_fail_days = 0
hydro_fail_count = 0

! **************************************************************************** !
! 4. Read the daily weather and management input data.
! **************************************************************************** !
OPEN (5,FILE = file_name, STATUS='OLD', IOSTAT=ios)
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not open weather input file: ', TRIM(file_name)
  STOP
END IF
! Skip header line (non-numeric first token), if present.
READ(5,'(A)',IOSTAT=ios) input_line
IF (ios == 0) THEN
  READ(input_line,*,IOSTAT=ios) prestart_precip_mm
  IF (ios == 0) THEN
    BACKSPACE(5)
  END IF
END IF
DO iD = 1, NDays
  current_day = iD
  READ(5,'(A)',IOSTAT=ios) input_line
  IF (ios /= 0) THEN
    WRITE (*,*) 'ERROR: Could not read input line for day ', iD
    STOP
  END IF
  CALL parse_csv_reals(input_line, NinVar, array_input_data(iD,1:NinVar))
  Year    (iD) = int(array_input_data(iD,1)); Month(iD) = int(array_input_data(iD,2)); Day(iD) = int(array_input_data(iD,3))
  s_hr    (iD) = array_input_data(iD,4); Tmax_d(iD) = array_input_data(iD,5); Tmin_d(iD) = array_input_data(iD,6)
  RHmax_d (iD) = array_input_data(iD,7); RHmin_d(iD) = array_input_data(iD,8)
  w_d     (iD) = array_input_data(iD,9)
  WS_d(iD) = array_input_data(iD,10)
  Irr_d(iD) = array_input_data(iD,11)
  Rad_d(iD) = array_input_data(iD,12)
END DO
CLOSE (5)

CALL initialize_phenology_params(phen_params)
phenology_cultivar_file = 'input/phenology_cultivar.txt'
CALL read_phenology_cultivar(phenology_cultivar_file, phen_params, phen_io)
IF (phen_io /= 0) THEN
  WRITE (*,*) 'ERROR: Could not read phenology cultivar file: ', TRIM(phenology_cultivar_file)
  STOP
END IF
CALL initialize_phenology_state(phen_state, phen_params%sowing_year, phen_params%sowing_month, phen_params%sowing_day, lat_s)
CALL bootstrap_phenology_state(phen_params, phen_state, Year(1), Month(1), Day(1))

! Transport the three optional cultivar canopy-architecture rows of
! input/phenology_cultivar.txt into the SHARED variables growth_module.f90
! reads, before its first use of any of them. The reader already leaves the
! phenology parameters at the compiled literals (73.9, 20.8, 1.45) when a row
! is absent, so this copy is unconditional and is a plain assignment and never
! arithmetic.
clr_phyllochron_Cd = phen_params%phyllochron_degCd
Blade_L_max_cm      = phen_params%blade_l_max_cm
Blade_W_max_cm       = phen_params%blade_w_max_cm
! Exact-equality branch: when the cultivar scalar equals the compiled default,
! the REFERENCE array is copied with no arithmetic, because IEEE floating point
! does not guarantee that multiplying and dividing by the same number returns
! the operand, so the branch skips the multiply entirely rather than relying on
! it to cancel.
IF (Blade_L_max_cm == 20.8_wp) THEN
  cba_l_irank = cba_l_irank_ref
ELSE
  cba_l_irank = cba_l_irank_ref * (Blade_L_max_cm / 20.8_wp)
END IF
IF (Blade_W_max_cm == 1.45_wp) THEN
  cba_w_irank = cba_w_irank_ref
ELSE
  cba_w_irank = cba_w_irank_ref * (Blade_W_max_cm / 1.45_wp)
END IF

! Transport the two optional cultivar sheath-geometry rows, same plumbing as
! the block above. The sheath diameter is a plain assignment, since the
! quantity IS the trait and there is no profile shape to preserve; the maximum
! sheath length uses the same exact-equality branch as the blade dimensions,
! for the same floating-point reason. The internode profile is not touched
! here; the reason is at its SHARED.f90 declaration.
Sheath_L_max_cm = phen_params%sheath_l_max_cm
Sheath_D_cm     = phen_params%sheath_d_cm
IF (Sheath_L_max_cm == 17.0_wp) THEN
  cse_ls_irank = cse_ls_irank_ref
ELSE
  cse_ls_irank = cse_ls_irank_ref * (Sheath_L_max_cm / 17.0_wp)
END IF
cse_ds2 = Sheath_D_cm
! tillering-capacity row (plain assignment; 0 = no cap, so the default path changes
! nothing and prints nothing -- the banner below is ON-only, file-set and run.log identical OFF).
Axes_max_per_plant = phen_params%axes_max_per_plant
IF (Axes_max_per_plant > 0.0_wp) THEN
  WRITE (*,'(A,F8.3,A,F10.2,A)') ' Axes_max_per_plant = ', Axes_max_per_plant, &
    ' axes plant-1 (cultivar tillering capacity; crop bound = ', Axes_max_per_plant * plant_density, ' axes m-2)'
END IF
! rec 150: hand the driver switch to the phenology engine (plain assignment; 0 = incumbent
! path untouched, prints nothing -- the banner is ON-only).
phen_params%dev_factor_window = pheno_dev_factor_window
IF (pheno_dev_factor_window == 1) THEN
  WRITE (*,'(A,F6.1,A,F6.1,A)') ' rec 150 = 1: VF/DF act from emergence to terminal spikelet (400 degCd TDU),' // &
    ' DTDU = DTT*min(VF,DF) (McMaster et al. 2008 Eq. 1); Photo_Sens = ', phen_params%photo_sens, &
    ', Vern_Factor = ', phen_params%vern_factor, ' (both 0 -> identical to legacy)'
END IF

! ----------------------------------------------------------------------!
! Open the soil-water driver file and read its data.
! ----------------------------------------------------------------------!
OPEN (30,FILE='input/driver_swtr.txt',STATUS='OLD',IOSTAT=ios)
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not open input/driver_swtr.txt'; STOP
END IF
READ (30,'(A)') file_name_swtr; READ (30,*) NTsteps; READ (30,*) root_dep; READ (30,*) wt_pres; READ (30,*) NSoil_layers
WRITE (*,*) 'Number of soil layers = ', NSoil_layers

! Allocate soil arrays
ALLOCATE (hydro_root_dist_diag(NSoil_layers), hydro_alpha_diag(NSoil_layers))
hydro_root_dist_diag = 0.0_wp; hydro_alpha_diag = 0.0_wp
ALLOCATE (hydro_alpha_alt_diag(NSoil_layers)); hydro_alpha_alt_diag = 0.0_wp   ! rec 80 diagnostics
! rec 81: committed per-layer head / K twins of h1_committed_m (filled in both modes at the commit sites).
ALLOCATE (h_committed_m(NSoil_layers), k_committed_m_d(NSoil_layers))
h_committed_m = 0.0_wp; k_committed_m_d = 0.0_wp
! Layered soil-nitrogen states and the read-only hydrology aggregates.
! Allocated unconditionally; the single-pool path never reads them, and the
! hydrology accumulators never feed back into the water physics.
ALLOCATE (NO3_lay(NSoil_layers), NH4_lay(NSoil_layers))
ALLOCATE (q_iface_day_mm(NSoil_layers), rootup_layer_day_mm(NSoil_layers))
ALLOCATE (hydro_q_iface_m(NSoil_layers), hydro_sink_layer_m(NSoil_layers))
NO3_lay = 0.0_wp; NH4_lay = 0.0_wp
q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp
hydro_q_iface_m = 0.0_wp; hydro_sink_layer_m = 0.0_wp
! Record 139: the per-layer diagnostic arrays of
! output/daily_n_uptake_terms.csv, allocated unconditionally; mode 0 never
! opens the file that reads them.
ALLOCATE (nud_mf_lay(NSoil_layers), nud_diff_lay(NSoil_layers), nud_act_lay(NSoil_layers))
ALLOCATE (nud_cap_lay(NSoil_layers), nud_real_lay(NSoil_layers))
nud_mf_lay = 0.0_wp; nud_diff_lay = 0.0_wp; nud_act_lay = 0.0_wp
nud_cap_lay = 0.0_wp; nud_real_lay = 0.0_wp
! rec 139 mode 2: per-layer active-bound flag array, same allocate-
! unconditionally/cheap/read-only-at-OFF discipline as the mode-1 arrays above.
ALLOCATE (nud_bound_lay(NSoil_layers))
nud_bound_lay = 0
! rec 139 mode 3: per-layer diagnostic arrays for
! output/daily_n_uptake_mode3.csv, same allocate-unconditionally discipline.
ALLOCATE (nud_fl3_lay(NSoil_layers), nud_rld3_lay(NSoil_layers), nud_r13_lay(NSoil_layers), nud_geom3_lay(NSoil_layers))
ALLOCATE (nud_cno33_lay(NSoil_layers), nud_cnh43_lay(NSoil_layers))
nud_fl3_lay = 0.0_wp; nud_rld3_lay = 0.0_wp; nud_r13_lay = 0.0_wp; nud_geom3_lay = 0.0_wp
nud_cno33_lay = 0.0_wp; nud_cnh43_lay = 0.0_wp
ALLOCATE (zb(NSoil_layers+1), dz(NSoil_layers), tk(NSoil_layers))
ALLOCATE (theta_1500(NSoil_layers), theta_33(NSoil_layers), theta_s33(NSoil_layers), theta_0(NSoil_layers))
ALLOCATE (vwc_sum(NSoil_layers), theta(NSoil_layers), wc(NSoil_layers), root_frac(NSoil_layers))
ALLOCATE (inflx(NSoil_layers), outflx(NSoil_layers), net_flx(NSoil_layers))
ALLOCATE (h_m(NSoil_layers), hm(NSoil_layers), hg(NSoil_layers), k_s(NSoil_layers), k_x(NSoil_layers))
ALLOCATE (theta_r(NSoil_layers), vg_alpha(NSoil_layers), vg_n(NSoil_layers), vg_m(NSoil_layers), vg_l(NSoil_layers))
ALLOCATE (swc_layers(NSoil_layers), swc_layers_daily(NDays,NSoil_layers))
ALLOCATE (theta_raw(NSoil_layers))
ALLOCATE (swc_layers_raw_daily(NDays,NSoil_layers))
ALLOCATE (swc_layers_mean_daily(NDays,NSoil_layers))
ALLOCATE (theta_start(NSoil_layers))
ALLOCATE (theta_init(NSoil_layers))
ALLOCATE (theta_day_sum(NSoil_layers))
ALLOCATE (theta_day0(NSoil_layers))

! Allocate Richards Solver Arrays (Shared Module)
ALLOCATE (r_a(NSoil_layers), r_b(NSoil_layers), r_c(NSoil_layers), r_r(NSoil_layers))
ALLOCATE (r_h(NSoil_layers), r_k(NSoil_layers), r_cap(NSoil_layers))
ALLOCATE (r_theta(NSoil_layers), sink_source(NSoil_layers))

OPEN (newunit=unit_rootcomp,file='output/daily_root_comp_diag.csv',status='replace',action='write')
WRITE (unit_rootcomp,'(A)') 'day,omega,F_raw,denom,f_comp,root_uptake_mm,trans_demand_mm'
OPEN (newunit=unit_rootcomp_lyr,file='output/daily_root_comp_layers.csv',status='replace',action='write')
WRITE (unit_rootcomp_lyr,'(A)') 'day,layer,root_dist,alpha,uncomp_weight,comp_weight,comp_uptake_mm'
! Record 78: the per-outer-iteration census of the water-stress factor, written
! in both modes into its own file and read by no gate and no process. One row
! per day and outer iteration: the factor the demand pass used, the factor the
! hydrology produced, the Feddes deliverable fraction at the last hydrology
! call entry, the actual transpiration and the four potentials (sent to the
! solver, unstressed, demand-pass and pass-2 stressed).
OPEN (newunit=unit_ft_outer,file='output/daily_ft_outer_diag.csv',status='replace',action='write')
WRITE (unit_ft_outer,'(A)') 'day,outer_iter,ft_pot_unstressed,f_t_in,f_t_out,F_raw_lastcall,' // &
  'trans_act_mm,trans_pot_sent_mm,trans_pot_unstressed_mm,trans_pot_demandpass_mm,' // &
  'trans_pot_pass2_stressed_mm,unstressed_sum_resid_mm,accum_diag,converged'
! Record 81 diagnostics: two separate files, with headers written in both
! modes, read by no gate and no process. daily_spac_diag.csv carries one row
! per day with the three root-weighted soil-potential conventions (computed in
! both modes from the committed profile), the one used, the leaf-water solve
! census counters, and the water ledger of the coupling. spac_census.csv
! carries one row per gas-exchange call on the writer pass when the record is
! on.
OPEN (newunit=unit_spac_daily,file='output/daily_spac_diag.csv',status='replace',action='write')
WRITE (unit_spac_daily,'(A)') 'day,spac_mode,psi_f_MPa,K_plant_mol_m2_s_MPa,psi_s_valid,psi_s_fallback,' // &
  'psi_s_rootfrac_MPa,psi_s_alpha_MPa,psi_s_kwt_MPa,psi_s_used_MPa,spac_calls,spac_iters_max,' // &
  'spac_resid_max_MPa,trans_pot_unstressed_mm,E_spac_pass2_mm,trans_act_mm,supply_shortfall_mm'
! Record 173 on-only diagnostic: one row per day with the demand actually sent
! to the Richards sink, what the Feddes sink delivered, the committed pass-2
! canopy transpiration, their residual and ratio, the outer-pass count and the
! relaxation factor.
IF (spac_sink_coupling >= 1) THEN
  OPEN (newunit=unit_ssc_daily,file='output/daily_spac_sink_coupling.csv',status='replace',action='write')
  WRITE (unit_ssc_daily,'(A)') 'day,outer_iters,outer_converged,demand_sent_to_sink_mm,' // &
    'sink_delivered_mm,E_spac_pass2_mm,resid_E_minus_sink_mm,f_trans_repartition,relax_scale,' // &
    'psiL_fallback_accepts_cum'
  ssc_daily_open = .TRUE.
END IF

! Record 80 census, written in both modes and read by nothing: one row per day
! at the entry of the last Richards call of the day, carrying the demand handed
! to the solver and the daily unstressed potential in cm d-1, the p2 node and
! node class the active rule used, the node and class the other rule would have
! used at the same pressure heads, the largest absolute difference between the
! two root-uptake reduction profiles, and both profiles.
OPEN (newunit=unit_feddes_node,file='output/daily_feddes_node_diag.csv',status='replace',action='write')
swc_line = 'day,feddes_p2_hydrus,ft_pot_unstressed,tp_sent_cm_d,tp_pot_cm_d,p2_used_m,node_class_used,' // &
           'p2_other_rule_m,node_class_other_rule,alpha_maxabsdiff_rules'
DO I = 1, NSoil_layers; WRITE (swc_tmp,'(A,I0)') ',alpha_used_L', I; swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
DO I = 1, NSoil_layers; WRITE (swc_tmp,'(A,I0)') ',alpha_other_rule_L', I; swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
WRITE (unit_feddes_node,'(A)') TRIM(swc_line)

OPEN (newunit=unit_rootdiag,file='output/daily_root_diagnostics.csv',status='replace',action='write')
! header was missing f_pr, shifting every name-based
! parse of this file by one column (the data row writes day,root_dep,f_t,f_pr,theta.).
root_header = 'day,root_dep_m,f_t,f_pr'
DO I = 1, NSoil_layers
  WRITE (swc_tmp,'(",theta",I0)') I
  root_header = TRIM (root_header) // TRIM (swc_tmp)
END DO
WRITE (unit_rootdiag,'(A)') TRIM (root_header)

DO I = 1, NSoil_layers; READ (30,*) theta (I); END DO
theta_init = theta
! Initialize Richards state from the input theta on day 1.
swc_layers = theta
DO I = 1, NSoil_layers + 1; READ (30,*) zb (I); zb (I) = -1.0_wp * zb (I); END DO
! Optional trailing field: the user-controlled maximum root depth (m). The
! positional read means a driver_swtr.txt file that ends at the last layer
! interface hits end of file here and the depth defaults to 1.5 m.
READ (30,*,IOSTAT=ios) root_max_depth
IF (ios /= 0) root_max_depth = 1.5_wp
CLOSE (30)
depth_lim  = ABS(zb (NSoil_layers+1))

DO I = 1, NSoil_layers
  tk (I)  = ABS(zb (I+1) - zb (I))
  dz (I)  = 0.5_wp * (zb (I) + zb (I+1))
END DO

! -------------------------------------------------------------------------- !
! Soil hydraulic parameters: the van Genuchten-Mualem parameters per layer from
! input/SoilLayer_vgm.csv. This is the single supported source; no soil-texture
! file and no pedotransfer-function selector is used.
! File format (CSV header): Layer,theta_r,theta_s,alpha_m-1,n,l,ksat_m_d
! -------------------------------------------------------------------------- !
ALLOCATE (vg_in(NSoil_layers,7))
OPEN (41,FILE = 'input/SoilLayer_vgm.csv', STATUS='OLD', IOSTAT=ios)
IF (ios /= 0) THEN
  WRITE (*,*) 'ERROR: Could not open input/SoilLayer_vgm.csv'; STOP
END IF
READ(41,*)
DO I = 1, NSoil_layers
  BLOCK
    CHARACTER(LEN=256) :: line
    REAL(wp) :: vals(10)
    INTEGER :: nvals
    READ(41,'(A)',IOSTAT=ios) line
    IF (ios /= 0) THEN
      WRITE (*,*) 'ERROR: Could not read SoilLayer_vgm.csv line for layer ', I
      STOP
    END IF
    nvals = count_csv_vals(line)
    IF (nvals < 7) THEN
      WRITE (*,*) 'ERROR: SoilLayer_vgm.csv needs 7 columns (Layer,theta_r,theta_s,alpha,n,l,ksat_m_d).'
      STOP
    END IF
    CALL parse_csv_reals(line, nvals, vals(1:nvals))
    vg_in(I,1) = vals(2) ! theta_r
    vg_in(I,2) = vals(3) ! theta_s
    vg_in(I,3) = vals(4) ! alpha (m-1)
    vg_in(I,4) = vals(5) ! n
    vg_in(I,5) = vals(6) ! l
    vg_in(I,6) = vals(7) ! ksat (m/d)
    ! Publish the saturated water content of layer 1 so that the
    ! soil-surface-resistance formula reads THE SAME number Richards uses
    ! (record 75). Assigned unconditionally and inert unless that record is on.
    IF (I == 1) vgm_theta_s_L1 = vals(3)
  END BLOCK
END DO
CLOSE (41)

! Assign the van Genuchten-Mualem parameters and derive field capacity and
! wilting point from the retention curve
! theta(h) = theta_r + (theta_s - theta_r) / (1 + (alpha |h|)^n)^m, the same
! form as hydro.f90: field capacity at h = -3.3 m and wilting point at
! h = -150 m.
DO I = 1, NSoil_layers
  theta_r(I)  = vg_in(I,1)
  theta_0(I)  = vg_in(I,2)
  vg_alpha(I) = vg_in(I,3)
  vg_n(I)     = vg_in(I,4)
  vg_l(I)     = vg_in(I,5)
  vg_m(I)     = 1.0_wp - 1.0_wp / vg_n(I)
  k_s(I)      = vg_in(I,6)
  IF (k_s(I) <= 0.0_wp) THEN
    WRITE (*,*) 'ERROR: SoilLayer_vgm.csv ksat_m_d must be > 0 for layer ', I
    STOP
  END IF
  theta_33(I)   = theta_r(I) + (theta_0(I) - theta_r(I)) / &
                  (1.0_wp + (vg_alpha(I) * 3.3_wp)**vg_n(I))**vg_m(I)
  theta_1500(I) = theta_r(I) + (theta_0(I) - theta_r(I)) / &
                  (1.0_wp + (vg_alpha(I) * 150.0_wp)**vg_n(I))**vg_m(I)
  theta_s33(I)  = theta_0(I) - theta_33(I)
END DO

! Write VGM parameters for diagnostics (per run).
OPEN (newunit=unit_soilprof, file='output/vgm_params.csv', status='replace', action='write')
WRITE (unit_soilprof,'(A)') 'layer,theta_r,theta_s,alpha_m-1,n,m,l,ksat_m_d'
DO I = 1, NSoil_layers
  WRITE (unit_soilprof,'(I0,7(",",G14.6))') I, theta_r(I), theta_0(I), vg_alpha(I), &
    vg_n(I), vg_m(I), vg_l(I), k_s(I)
END DO
CLOSE (unit_soilprof)
WRITE (*,*) 'Wrote output/vgm_params.csv (VGM parameters from input/SoilLayer_vgm.csv)'

! (Surface crust adjustment removed; revert to base Ksat profile.)

! -------------------------------------------------------------------------- !
! Initialize the water-content profile from the bottom head if requested.
! Mode 1: a hydrostatic profile from the bottom head (h + z constant).
! Mode 2: a monotonic profile from the top water content to the bottom one.
! -------------------------------------------------------------------------- !
IF (init_profile_mode > 0) THEN
  BLOCK
    REAL(wp) :: z_top, z_bot, z_mid, frac
    REAL(wp) :: h_bot
    REAL(wp) :: th_top, th_bot, th_s, th_r, alpha_vg, n_vg, m_vg, Se, H_suc
    INTEGER :: i

    z_top = 0.5_wp * (zb(1) + zb(2))
    z_bot = 0.5_wp * (zb(NSoil_layers) + zb(NSoil_layers+1))
    h_bot = init_bottom_head_m
    th_top = theta(1)

    ! Compute theta at bottom from specified bottom head.
    th_s = MAX(1.0e-8_wp, theta_0(NSoil_layers))
    th_r = MAX(0.0_wp, theta_r(NSoil_layers))
    alpha_vg = MAX(1.0e-9_wp, vg_alpha(NSoil_layers))
    n_vg = MAX(1.01_wp, vg_n(NSoil_layers))
    m_vg = MAX(1.0e-6_wp, vg_m(NSoil_layers))
    H_suc = MAX(1.0e-6_wp, -h_bot)
    Se = 1.0_wp / (1.0_wp + (alpha_vg * H_suc)**n_vg)**m_vg
    th_bot = th_r + Se * (th_s - th_r)
    th_bot = MIN(th_s, MAX(th_r, th_bot))

    DO i = 1, NSoil_layers
      z_mid = 0.5_wp * (zb(i) + zb(i+1))
      ! Monotonic: increase theta from top to bottom.
      frac = (z_mid - z_top) / MAX(1.0e-9_wp, (z_bot - z_top))
      frac = MIN(1.0_wp, MAX(0.0_wp, frac))
      theta(i) = th_top + frac * (th_bot - th_top)
      theta(i) = MIN(theta_0(i), MAX(theta_r(i), theta(i)))
    END DO
    theta_init = theta
    swc_layers = theta
  END BLOCK
END IF

! Initialize soil evaporation parameters for stage-1/2 (top 10 cm).
! Extend the evaporation layer to 0.10 m (covers layer 1 fully + top 5 cm of layer 2).
BLOCK
  REAL(wp) :: evap_depth_m, w1, w2, theta0_eff, theta1500_eff
  evap_depth_m = 0.10_wp
  w1 = MIN(tk(1), evap_depth_m)
  w2 = MAX(0.0_wp, evap_depth_m - w1)
  IF (w2 > 0.0_wp) w2 = MIN(w2, tk(2))
  theta0_eff = (theta_0(1) * w1 + theta_0(2) * w2) / MAX(1.0e-9_wp, (w1 + w2))
  theta1500_eff = (theta_1500(1) * w1 + theta_1500(2) * w2) / MAX(1.0e-9_wp, (w1 + w2))
  TEW_top_mm = 1000.0_wp * MAX(0.0_wp, (theta0_eff - theta1500_eff) * evap_depth_m)
END BLOCK
REW_top_mm = MIN(8.0_wp, TEW_top_mm) ! default readily evaporable water (mm)
De_top_mm  = 0.0_wp

! **************************************************************************** !
! 6. Start of the daily simulation loop.
! **************************************************************************** !
IF (prestart_precip_mm > 0.0_wp) THEN
  ! Apply a one-time wetting event before day 1 (e.g., wet day prior to simulation).
  pr_eff = prestart_precip_mm / 1000.0_wp ! mm -> m
  trans_act = 0.0_wp
  evapo_act = 0.0_wp
  Ir_mm_d = 0.0_wp
  ponding_mm = 0.0_wp
  top_bc_type = 0
  h_top = 0.0_wp
  hydro_failed = .FALSE.
  CALL hydro (swc_layers, 0)
  hydro_fail_events_all = hydro_fail_events_all + hydro_fail_events_call
  IF (hydro_fail_events_call > 0) hydro_fail_calls_all = hydro_fail_calls_all + 1
  swc_layers = theta
END IF

! -------------------------------------------------------------------------- !
! Hydrostatic equilibration: relax the initial water-content profile to a
! no-flux hydrostatic state before day 1, so that day 1 does not open with a
! large redistribution jump.
! -------------------------------------------------------------------------- !
! The save and restore block below reads the effective precipitation and
! irrigation before the first simulated day, and irrigation is not assigned
! until the daily water-supply step, so both are defined here as zero supply.
irr     = 0.0_wp
irr_eff = 0.0_wp
IF (use_soil_water_stress .AND. init_profile_mode /= 2) THEN
  BLOCK
    REAL(wp) :: pr_eff_prev, irr_eff_prev, pr_prev, irr_prev
    REAL(wp) :: trans_prev, evap_prev, ir_prev
    INTEGER :: top_bc_prev
    REAL(wp) :: h_top_prev
    pr_eff_prev = pr_eff; irr_eff_prev = irr_eff; pr_prev = pr; irr_prev = irr
    trans_prev = trans_act; evap_prev = evapo_act; ir_prev = Ir_mm_d
    top_bc_prev = top_bc_type; h_top_prev = h_top
    pr_eff = 0.0_wp; irr_eff = 0.0_wp; pr = 0.0_wp; irr = 0.0_wp
    trans_act = 0.0_wp; evapo_act = 0.0_wp; Ir_mm_d = 0.0_wp
    ponding_mm = 0.0_wp; top_bc_type = 0; h_top = 0.0_wp
    hydro_failed = .FALSE.
    CALL hydro (swc_layers, 0)
    hydro_fail_events_all = hydro_fail_events_all + hydro_fail_events_call
    IF (hydro_fail_events_call > 0) hydro_fail_calls_all = hydro_fail_calls_all + 1
    swc_layers = theta
    pr_eff = pr_eff_prev; irr_eff = irr_eff_prev; pr = pr_prev; irr = irr_prev
    trans_act = trans_prev; evapo_act = evap_prev; Ir_mm_d = ir_prev
    top_bc_type = top_bc_prev; h_top = h_top_prev
  END BLOCK
END IF

pheno_stage  (:) = 0; plant_height (:) = gParam_hgt; lai_sim (:) = zero
Total_trans = 0.0_wp; Total_assim = 0.0_wp; Total_gsw_1 = 0.0_wp; Total_N_ass = 0.0_wp; Ci_last = 0.0_wp
T_leaf_mean_total = 0.0_wp; T_air_mean_total = 0.0_wp; Ci_mean_total = 0.0_wp
count_neg_hc = 0; count_neg_hs = 0; count_le_gt_rn = 0; count_et_gt_avail = 0; count_rn_neg = 0; count_le_neg = 0
total_precip = 0.0_wp; total_irrig = 0.0_wp; total_et_out = 0.0_wp; total_drain = 0.0_wp; total_runoff = 0.0_wp
total_supply_eff = 0.0_wp; total_intercept = 0.0_wp
total_shortfall_mm = 0.0_wp
total_wb_resid = 0.0_wp
total_evap_head_exfil = 0.0_wp; total_evap_atm_flux = 0.0_wp
total_hydro_top_flux_m = 0.0_wp; total_hydro_root_sink_m = 0.0_wp; total_hydro_bottom_flux_m = 0.0_wp
total_balance_soil_mm = 0.0_wp; total_balance_surface_mm = 0.0_wp; total_balance_mm = 0.0_wp
total_soil_mb_resid_mm = 0.0_wp
total_soil_evap = 0.0_wp
total_trans_act = 0.0_wp
total_pond_evap = 0.0_wp
total_hydro_mb_resid = 0.0_wp
total_hydro_mb_correction = 0.0_wp
cum_tr = 0.0_wp; cum_e = 0.0_wp; cum_ir = 0.0_wp; cum_et = 0.0_wp; cum_et_ir = 0.0_wp; cum_precip = 0.0_wp
hydro_ledger = HydroCumLedger()

! Initialize biomass pools; start from zero before emergence.
CALL initialize_growth_state(W_green, W_dead, W_stem, W_root, W_storage, gParam_lai)
! Task 7: yield-component sub-pools start empty; grain number unset until ear_end.
W_grain = 0.0_wp; W_chaff = 0.0_wp; kernel_number_m2 = 0.0_wp; grain_number_set = .FALSE.
N_leaf = 0.0_wp
N_dead_leaf = 0.0_wp
N_stem = 0.0_wp
N_root = 0.0_wp
N_storage = 0.0_wp
N_NO3_root = 0.0_wp
N_NO3_shoot = 0.0_wp
N_assim_root_day = 0.0_wp
N_assim_shoot_day = 0.0_wp
f_photoresp_diag = 0.0_wp
N_assim_biochem_day = 0.0_wp
! -------------------------------------------------------------------- !
! Read the soil-nitrogen parameters from input/driver_N.csv.
! Format: value,label, one parameter per non-comment line; lines beginning
! with an exclamation mark are skipped. The values are read positionally in
! SHARED order.
! -------------------------------------------------------------------- !
BLOCK
  INTEGER :: unit_N_drv, ios_N, ip_N
  CHARACTER(LEN=512) :: line_N
  ! The mandatory record count stays 64, so a 64-row driver_N.csv reads exactly
  ! as before. Positions 65, 66 and 67 carry the MEASURED per-deck soil analysis
  ! (clay, calcium carbonate and the initial humus nitrogen) in place of a single
  ! transplanted soil; they are read only if present, and record 93 at or above 1
  ! refuses to run without them.
  INTEGER, PARAMETER :: N_PARAMS_SOIL = 71  ! last optional position (rec 153: N_soil_init_depth_m, row 71, same present-only convention as rows 65-67 and 68-70)
  INTEGER, PARAMETER :: N_PARAMS = 64   ! +4 (P2) Andrews; +2 (P5) Martre; +6 urea/NH4/NO3; +12 initial-N + residue decomposition; +5 (pH_surf, Kd_NH4, BD_surf, z_surf_N, Km_den); +13 root-scaled active uptake: 5 kinetics/root rows + 8 f_Nstatus(NNI) node values; +5 (layered-N) u_max_nh4_fw, Km_nh4, n_min_depth_m, cfl_frac, soilN_layered (appended AFTER existing rows so old positions keep meaning) — N_PARAMS itself is a STRUCTURAL row count, not a physical constant
  REAL(wp) :: N_pv(N_PARAMS_SOIL)
  N_pv = 0.0_wp; ip_N = 0
  OPEN(newunit=unit_N_drv, FILE='input/driver_N.csv', STATUS='OLD', IOSTAT=ios_N)
  IF (ios_N /= 0) THEN; WRITE(*,*) 'ERROR: Cannot open input/driver_N.csv'; STOP; END IF
  ! Sequential read: skip comment/blank lines, fill N_pv in order.
  DO WHILE (ip_N < N_PARAMS_SOIL)
    READ(unit_N_drv, '(A)', IOSTAT=ios_N) line_N
    IF (ios_N /= 0) EXIT
    line_N = ADJUSTL(line_N)
    IF (LEN_TRIM(line_N) == 0 .OR. line_N(1:1) == '!') CYCLE
    ip_N = ip_N + 1
    READ(line_N, *, IOSTAT=ios_N) N_pv(ip_N)
  END DO
  CLOSE(unit_N_drv)
  IF (ip_N < N_PARAMS) THEN
    WRITE(*,'(A,I0,A,I0,A)') 'ERROR: driver_N.csv has only ', ip_N, &
      ' values, expected ', N_PARAMS, '.'; STOP
  END IF
  ! Assign in exact order matching driver_N.csv (see file for sources).
  N_soil_init      = N_pv(1);  V_Nmax          = N_pv(2)
  Km_N             = N_pv(3);  trans_ref_N     = N_pv(4)
  SLN_min_N        = N_pv(5);  k_remob_leaf_N  = N_pv(6)
  k_remob_stem_N   = N_pv(7);  n_struct_stem_N = N_pv(8)
  k_N_can          = N_pv(9);  n_b_can         = N_pv(10)
  N_max_ref        = N_pv(11); Km_N_busch      = N_pv(12)
  N_mineral_rate   = N_pv(13); Ci_ref_N        = N_pv(14)
  f_postanth_stor  = N_pv(15); f_postanth_lf   = N_pv(16)
  f_postanth_st    = N_pv(17)
  f_shoot_min      = N_pv(18); f_shoot_max     = N_pv(19)
  K_partition_NO3  = N_pv(20); r_root_NR       = N_pv(21)
  n_stru_grain     = N_pv(22); dvs_grain_stru_end = N_pv(23)
  k_hyd_ref        = N_pv(24); Q10_hyd         = N_pv(25)
  k_nit_ref        = N_pv(26); Q10_nit         = N_pv(27)
  T_ref_hyd        = N_pv(28); T_ref_nit       = N_pv(29)
  N_soil_init_NH4  = N_pv(30); residue_C_frac  = N_pv(31)
  k_fast_ref       = N_pv(32); k_slow_ref      = N_pv(33)
  Q10_dec          = N_pv(34); T_ref_dec       = N_pv(35)
  theta_ref_dec    = N_pv(36); e_mic           = N_pv(37)
  CN_mic           = N_pv(38); k_mic           = N_pv(39)
  presow_Tsoil     = N_pv(40); presow_fwater   = N_pv(41)
  pH_surf          = N_pv(42); Kd_NH4          = N_pv(43)
  BD_surf          = N_pv(44); z_surf_N        = N_pv(45)
  Km_den           = N_pv(46)
  ! root-scaled active uptake
  ! rows appended after Km_den — V_Nmax (pos 2) and Km_N (pos 3) are RETIRED (read, unused).
  u_max_fw         = N_pv(47); fw_dw_root      = N_pv(48)
  Km_c_no3         = N_pv(49); aroot_floor     = N_pv(50)
  aroot_decline_d  = N_pv(51)
  ! The nodes of the nitrogen-status function of the nitrogen nutrition index:
  ! four pairs, interleaved positionally (Taulemesse et al. 2015 PLoS ONE
  ! 10(3):e0120291, Table 1 p.7 with the flux sentence p.8 normalized to the
  ! lowest nitrogen treatment; the per-row sources are in driver_N.csv).
  fNstat_nni(1) = N_pv(52); fNstat_f(1) = N_pv(53)
  fNstat_nni(2) = N_pv(54); fNstat_f(2) = N_pv(55)
  fNstat_nni(3) = N_pv(56); fNstat_f(3) = N_pv(57)
  fNstat_nni(4) = N_pv(58); fNstat_f(4) = N_pv(59)
  ! Layered soil-nitrogen and ammonium-uptake rows, with their sources in
  ! driver_N.csv: the ammonium transporter maximum rate and half-saturation
  ! constant (Goyal & Huffaker 1986 pp.1051 and 1053), the depth receiving
  ! background mineralization, the advection sub-step guard, and the
  ! layered-or-single-pool mode flag.
  u_max_nh4_fw     = N_pv(60); Km_nh4          = N_pv(61)
  n_min_depth_m    = N_pv(62); cfl_frac        = N_pv(63)
  soilN_layered    = NINT(N_pv(64))
  ! The MEASURED per-deck soil analysis for the humus mineralization of record
  ! 93, in the same class as the surface pH (position 42) and bulk density
  ! (position 44): a measured site soil analysis is an INPUT and not a fitted
  ! parameter. Absent rows leave the negative sentinel from SHARED.f90, and
  ! record 93 at or above 1 then stops with an error; there is no silent
  ! fallback.
  IF (ip_N >= 65) stics_argi_pct    = N_pv(65)
  IF (ip_N >= 66) stics_calc_pct    = N_pv(66)
  IF (ip_N >= 67) stics_nhum0_gN_m2 = N_pv(67)
  ! Record 139 mode 3: the measured, or class-centroid-inferred, per-deck texture
  ! for the sourced sandy-soil impedance law, in the same class as the soil
  ! analysis above - an input, not a fitted parameter. Absent rows leave the
  ! negative sentinel from SHARED.f90, and mode 3 then stops with an error.
  IF (ip_N >= 68) nud_clay_frac    = N_pv(68)
  IF (ip_N >= 69) nud_silt_frac    = N_pv(69)
  IF (ip_N >= 70) nud_bulk_density = N_pv(70)
  ! Record 153: the SAMPLING DEPTH the initial mineral-nitrogen measurement of
  ! row 1 was derived from, in m below the surface. An OPTIONAL row at the end of
  ! the file, read with the same present-only and explicit-reset idiom as rows 65
  ! to 70 and as the optional cultivar rows: an absent row is not an error and
  ! the field is explicitly reset to the negative sentinel rather than left at
  ! whatever the array holds. The sentinel, a non-positive value, or a value at
  ! or beyond the profile depth, all mean the whole profile, so record 153 mode 1
  ! is then arithmetically identical to mode 0 on that deck. This is a
  ! measurement basis and not a fitted parameter, in the same class as the
  ! surface pH (row 42) and bulk density (row 44).
  IF (ip_N >= 71) THEN
    n_init_depth_m = N_pv(71)
  ELSE
    n_init_depth_m = -1.0_wp
  END IF
  WRITE(*,'(A,F5.1,A,F4.2,A,F5.3)') ' N model: N_soil_init=', N_soil_init, &
    '  u_max_fw=', u_max_fw, '  k_remob_leaf=', k_remob_leaf_N
  WRITE(*,'(A,F4.2,A,F4.2,A,F4.1,A,F6.4)') ' N partition (Andrews): f_shoot_min=', f_shoot_min, &
    '  f_shoot_max=', f_shoot_max, '  K_partition_NO3=', K_partition_NO3, '  r_root_NR=', r_root_NR
  ! Layered-mode banner ONLY when the flag is on (incumbent stdout unchanged, G9).
  IF (soilN_layered == 1) THEN
    WRITE(*,'(A,F4.1,A,F5.3,A,F5.2,A,F4.2)') ' Soil-N LAYERED mode: u_max_nh4_fw=', u_max_nh4_fw, &
      '  Km_nh4=', Km_nh4, '  n_min_depth_m=', n_min_depth_m, '  cfl_frac=', cfl_frac
  END IF
  ! Record 139 hooks the layered per-layer uptake loop only; the single-pool path
  ! never reaches it, so being on there would be a silent no-op and is an error
  ! stop instead.
  IF (n_uptake_diffusion >= 1 .AND. soilN_layered < 1) THEN
    WRITE (*,'(A,I0)') ' ERROR rec 139: n_uptake_diffusion >= 1 requires driver_N.csv soilN_layered >= 1' // &
      ' (the layered per-layer uptake loop this record extends does not execute otherwise). Read soilN_layered = ', &
      soilN_layered
    ERROR STOP 'rec 139 requires soilN_layered >= 1'
  END IF
  ! Record 139 mode 3: a hard refusal rather than a silent fallback, the same
  ! discipline as the soil-analysis gate of record 93. The sourced impedance law
  ! of mode 3 needs the deck's own clay, silt and bulk density to compute the
  ! threshold water content (Olesen et al. 2001 p.93 Eq. (5)); a deck lacking all
  ! three rows cannot be run at mode 3.
  ! Record 168: the threshold water content IS NOW COMPUTED ON EVERY DECK, and
  ! the assignment is moved out of the mode-3 branch below. Inside that branch it
  ! was identically zero on every released deck, which record 168 mode 1 cannot
  ! use, since it would degrade Eq. (7a) to a threshold-free form, a third law
  ! with no source. Moving it changes no output at record 168 = 0: every other
  ! reader and printer of the threshold is itself gated on record 139 at or above
  ! 3, and the mode-3 refusal below is unchanged. The 42-deck identity gate is
  ! what proves that, not this comment. Olesen et al. (2001) Soil Sci 166:89-96
  ! p.93 Eq. (5), the same single law on every deck.
  IF (nud_clay_frac >= 0.0_wp .AND. nud_silt_frac >= 0.0_wp .AND. nud_bulk_density >= 0.0_wp) THEN
    nud_theta_th = nud_th_c1*nud_clay_frac + nud_th_c2*nud_clay_frac**2 + nud_th_c3*nud_silt_frac &
                 + nud_th_c4*nud_bulk_density + nud_th_c5*nud_bulk_density**2 + nud_th_c6
  END IF
  ! Record 168: refused rather than left a silent zero. Mode 1 needs all three
  ! texture rows; without them the threshold would stay zero and the impedance
  ! factor would become a threshold-free form that no page of Olesen prints.
  IF (mineral_impedance_basis >= 1) THEN
    IF (nud_clay_frac < 0.0_wp .OR. nud_silt_frac < 0.0_wp .OR. nud_bulk_density < 0.0_wp) THEN
      WRITE (*,'(A,3(F8.4,1X))') ' ERROR rec 168 mineral_impedance_basis = 1: requires driver_N.csv' // &
        ' rows 68-70 (clay_fraction, silt_fraction, bulk_density_topsoil_Mg_m3) all present, because' // &
        ' Olesen et al. 2001 p.93 Eq.(5) computes theta_th from them. Read (clay,silt,rho_b) = ', &
        nud_clay_frac, nud_silt_frac, nud_bulk_density
      ERROR STOP 'rec 168 mode 1 requires driver_N.csv rows 68-70'
    END IF
    ! ON-only echo (the OFF run.log carries no trace of this record).
    WRITE(*,'(A,F6.4,A,F6.4,A,F6.4,A,F7.4,A,F5.2)') ' rec 168 theta_th from this deck: clay_frac=', &
      nud_clay_frac, '  silt_frac=', nud_silt_frac, '  bulk_density=', nud_bulk_density, &
      '  theta_th=', nud_theta_th, '  CSIF slope=', nud_csif_slope
  END IF
END BLOCK

! -------------------------------------------------------------------- !
! Read the cultivar yield-component coefficients from
! input/cultivar_yield.csv. A positional read (value first, comment lines
! skipped), like driver_N.csv.
! -------------------------------------------------------------------- !
! A missing, unreadable or PARTIAL cultivar file is a HARD ERROR naming the
! first missing field and the file it was expected in. It is not a warning with
! a fall-through to compiled defaults: those defaults belong to one cultivar at
! one site, so a silent fallback would make any other site run that cultivar
! traits.
! The IOSTAT test alone is not sufficient, because Fortran list-directed input
! has two success-with-no-assignment cases - a null item, and a record starting
! with the list terminator - and because a value test is needed as well. This
! block therefore enforces THREE conditions per value: the IOSTAT is zero, the
! slot actually moved away from its unset sentinel, and the number is finite
! and strictly positive. Without them a sentinel row could leave a slot at a
! seed of 0.0, the all-present branch would fire, and the run would exit 0 with
! a fruiting efficiency and a potential kernel mass of zero, i.e. a silent
! zero-grain deck. Behavior is unchanged on every well-formed deck.
BLOCK
  INTEGER :: unit_cy, ios_cy, ip_cy
  CHARACTER(LEN=512) :: line_cy
  INTEGER, PARAMETER :: CY_PARAMS = 10  ! 6 required + 2 optional (rec 125 envelope, rows 7-8)
                                        ! Plus two optional rows (the grain-filling clock, rows 9 and 10), which
                                        ! promote the two compiled grain-filling constants, the base temperature and
                                        ! the filling duration, to CULTIVAR-FILE TRAITS. The accepted record lengths
                                        ! are 6, 8 and 10 - both or neither in each optional pair - so a six-row or
                                        ! eight-row file is unaffected.
  REAL(wp), PARAMETER :: CY_UNSET = -HUGE(1.0_wp)   ! unreachable by any successful read
  CHARACTER(LEN=20), PARAMETER :: CY_NAME(CY_PARAMS) = &
    [ CHARACTER(LEN=20) :: 'fruiting_efficiency', 'TKW_potential', 'chaff_share', &
                           'ear_end_dvs2', 'storage_onset_dvs2', 'storage_peak_frac', &
                           'TKW_env_min', 'TKW_env_max', &
                           'gfb_tbase_C', 'gfb_dfill_Cd' ]
  REAL(wp) :: cy_pv(CY_PARAMS)
  ! Record 177: the NAME-KEYED optional-trait scan. The cultivar file's own format
  ! is value,parameter_name,unit,source, as its own header line says, so the
  ! parameter NAME is already field 2 of every row. A line whose token is a known
  ! name-keyed trait is consumed WITHOUT advancing the positional counter, so
  ! such a row may sit anywhere in the file and can never trip the partial-record
  ! error. Rows 1 to 10 stay strictly positional, since making them name-keyed
  ! would change how every existing file parses.
  CHARACTER(LEN=64) :: cy_tok
  INTEGER :: ic1_cy, ic2_cy
  CHARACTER(LEN=*), PARAMETER :: CY_TRAIT_GNP = 'grainN_pot_mgN_kernel'
  ! The slots are seeded with a value that CANNOT be produced by a successful
  ! read, and not with 0.0. See the block above for why: a list-directed read has
  ! two success-with-no-assignment cases, and with a 0.0 seed both of them
  ! install a silent and physically meaningful zero.
  cy_pv = CY_UNSET; ip_cy = 0
  OPEN(newunit=unit_cy, FILE='input/cultivar_yield.csv', STATUS='OLD', IOSTAT=ios_cy)
  IF (ios_cy /= 0) THEN
    WRITE(*,*) ' ERROR: cannot open input/cultivar_yield.csv (relative to the case directory).'
    WRITE(*,*) '        This deck supplies NO cultivar record. All ', CY_PARAMS, ' values are'
    WRITE(*,*) '        required, in this order, one numeric value per non-comment line:'
    DO ip_cy = 1, CY_PARAMS
      WRITE(*,*) '          ', ip_cy, ' ', TRIM(CY_NAME(ip_cy))
    END DO
    WRITE(*,*) '        Refusing to fall back to the compiled-in cv. Ruta defaults, which'
    WRITE(*,*) '        belong to the Rehovot WheatDryFACE site and to no other.'
    ERROR STOP 'cultivar_yield.csv missing: input/cultivar_yield.csv'
  ELSE
    DO
      READ(unit_cy, '(A)', IOSTAT=ios_cy) line_cy
      IF (ios_cy /= 0) EXIT
      line_cy = ADJUSTL(line_cy)
      IF (LEN_TRIM(line_cy) == 0 .OR. line_cy(1:1) == '!') CYCLE
      ! ---- Record 177: NAME-KEYED optional trait rows, not positional ----
      ! Field 2 of the line is the parameter name the file itself writes. A line
      ! whose name token matches a known name-keyed trait is consumed here and does
      ! not advance the positional counter.
      cy_tok = ''
      ic1_cy = INDEX(line_cy, ',')
      IF (ic1_cy > 0) THEN
        ic2_cy = INDEX(line_cy(ic1_cy+1:), ',')
        IF (ic2_cy > 1) THEN
          cy_tok = ADJUSTL(line_cy(ic1_cy+1:ic1_cy+ic2_cy-1))
        ELSE IF (ic2_cy == 0) THEN
          cy_tok = ADJUSTL(line_cy(ic1_cy+1:))
        END IF
      END IF
      IF (TRIM(cy_tok) == CY_TRAIT_GNP) THEN
        READ(line_cy, *, IOSTAT=ios_cy) gnp_row_value
        IF (ios_cy /= 0 .OR. .NOT. (gnp_row_value > 0.0_wp) .OR. &
            gnp_row_value > HUGE(1.0_wp)/2.0_wp) THEN
          WRITE(*,'(A)') ' ERROR rec 177: input/cultivar_yield.csv row ' // CY_TRAIT_GNP // &
            ' did not supply a finite positive number (mg N kernel-1).'
          WRITE(*,*) '        Offending line: ', TRIM(line_cy)
          ERROR STOP 'rec 177: unreadable grainN_pot_mgN_kernel row in cultivar_yield.csv'
        END IF
        CYCLE
      END IF
      ! A positional row beyond the last defined slot is ignored, exactly as before (the old
      ! DO WHILE simply stopped reading at CY_PARAMS); no shipping file has one.
      IF (ip_cy >= CY_PARAMS) CYCLE
      ip_cy = ip_cy + 1
      READ(line_cy, *, IOSTAT=ios_cy) cy_pv(ip_cy)
      IF (ios_cy /= 0) THEN
        WRITE(*,*) ' ERROR: input/cultivar_yield.csv value ', ip_cy, ' (', &
                   TRIM(CY_NAME(ip_cy)), ') is not a readable number.'
        WRITE(*,*) '        Offending line: ', TRIM(line_cy)
        ERROR STOP 'cultivar_yield.csv unparseable value: input/cultivar_yield.csv'
      END IF
      ! An IOSTAT of zero does NOT mean a value was assigned. A record that is a NULL
      ! ITEM (it begins with the separator) or that begins with the list terminator
      ! completes with IOSTAT zero and leaves the variable exactly as it was. With a
      ! 0.0 seed that is a silent zero; with the unset sentinel it is detectable, and
      ! it is refused here.
      IF (cy_pv(ip_cy) <= CY_UNSET) THEN
        WRITE(*,*) ' ERROR: input/cultivar_yield.csv value ', ip_cy, ' (', &
                   TRIM(CY_NAME(ip_cy)), ') supplied NO number (null item or "/").'
        WRITE(*,*) '        Offending line: ', TRIM(line_cy)
        ERROR STOP 'cultivar_yield.csv empty value: input/cultivar_yield.csv'
      END IF
      ! Every one of the six values is strictly positive by construction (two rates,
      ! a mass, two development-stage nodes and two fractions), so a zero, a negative
      ! or a non-finite value expresses no cultivar. This is the check that catches a
      ! fruiting efficiency of zero, which would make the kernel number identically
      ! zero and the deck yield no grain at all.
      IF (.NOT. (cy_pv(ip_cy) > 0.0_wp) .OR. cy_pv(ip_cy) > HUGE(1.0_wp)/2.0_wp) THEN
        WRITE(*,*) ' ERROR: input/cultivar_yield.csv value ', ip_cy, ' (', &
                   TRIM(CY_NAME(ip_cy)), ') = ', cy_pv(ip_cy), &
                   ' is not a finite positive number.'
        WRITE(*,*) '        Offending line: ', TRIM(line_cy)
        ERROR STOP 'cultivar_yield.csv non-positive value: input/cultivar_yield.csv'
      END IF
    END DO
    CLOSE(unit_cy)
    IF (ip_cy == 6 .OR. ip_cy == 8 .OR. ip_cy == CY_PARAMS) THEN
      fruiting_efficiency = cy_pv(1); TKW_potential = cy_pv(2)
      chaff_share         = cy_pv(3); ear_end_dvs2  = cy_pv(4)
      ! Spike-window trait nodes (runtime cultivar inputs since;
      ! matrix finding: traits must never require a rebuild)
      storage_onset_dvs2  = cy_pv(5); storage_peak_frac = cy_pv(6)
      ! rec 125 envelope (rows 7-8, OPTIONAL; SHARED.f90): absent -> -1 and
      ! modes 2/3 refuse below; present -> read as g per 1000 kernels = mg kernel-1.
      IF (ip_cy >= 8) THEN
        kpt_env_min = cy_pv(7); kpt_env_max = cy_pv(8)
        IF (kpt_env_max < kpt_env_min) THEN
          WRITE(*,*) ' ERROR: input/cultivar_yield.csv TKW_env_max < TKW_env_min: ', kpt_env_max, kpt_env_min
          ERROR STOP 'cultivar_yield.csv envelope inverted (rows 7-8)'
        END IF
        ! Thousand-kernel-weight policy: when rows 7 and 8 are present the cultivar
        ! potential MUST lie inside the stated envelope; a value outside it is
        ! refused, never silently clipped or defaulted.
        IF (kpt_env_min > 0.0_wp .AND. (TKW_potential < kpt_env_min .OR. TKW_potential > kpt_env_max)) THEN
          WRITE(*,'(A,F7.2,A,F6.2,A,F6.2,A)') ' ERROR: input/cultivar_yield.csv TKW_potential ', TKW_potential, &
            ' g is outside the declared envelope rows 7-8 [', kpt_env_min, ', ', kpt_env_max, ']'
          ERROR STOP 'cultivar_yield.csv TKW_potential outside TKW_env_min..TKW_env_max (rows 7-8)'
        END IF
      END IF
      ! The grain-filling clock as a cultivar trait (rows 9 and 10, optional). Absent
      ! rows leave the base temperature and the filling duration at their stated
      ! values, which are the constants compiled in this release (6.618733 degC and
      ! 405.58109 degree-Cd above that base). Present rows replace them, and the
      ! provenance flags make the run.log say which deck did that. The duration is
      ! degree-Cd ABOVE the base temperature and not a base-0 thermal time, so a
      ! measured base-0 duration must be converted before it is entered. The three
      ! per-value conditions above already apply, and no further numeric bound is
      ! compiled in.
      IF (ip_cy == CY_PARAMS) THEN
        gfb_tbase = cy_pv(9);  gfb_tbase_from_file = .TRUE.
        gfb_dfill = cy_pv(10); gfb_dfill_from_file = .TRUE.
      END IF
    ELSE
      WRITE(*,*) ' ERROR: input/cultivar_yield.csv is PARTIAL: ', ip_cy, ' of ', CY_PARAMS, &
                 ' values were found.'
      WRITE(*,*) '        First missing field: value ', ip_cy + 1, ' = ', TRIM(CY_NAME(ip_cy + 1))
      WRITE(*,*) '        Expected order: fruiting_efficiency, TKW_potential, chaff_share,'
      WRITE(*,*) '                        ear_end_dvs2, storage_onset_dvs2, storage_peak_frac'
      WRITE(*,*) '                        [optional rows 7-8: TKW_env_min, TKW_env_max -- BOTH or neither]'
      WRITE(*,*) '                        [optional rows 9-10: gfb_tbase_C, gfb_dfill_Cd -- BOTH or neither;'
      WRITE(*,*) '                         rows 9-10 require rows 7-8 to be present, i.e. 6, 8 or 10 values]'
      WRITE(*,*) '        Refusing to fall back to the compiled-in cv. Ruta defaults, which'
      WRITE(*,*) '        belong to the Rehovot WheatDryFACE site and to no other.'
      ERROR STOP 'cultivar_yield.csv partial: input/cultivar_yield.csv'
    END IF
  END IF
  WRITE(*,'(A,F5.1,A,F5.1,A,F4.2,A,F4.2,A,F5.3,A,F5.3)') ' Cultivar yield: FE=', fruiting_efficiency, &
    '  TKW_pot=', TKW_potential, '  chaff_share=', chaff_share, '  ear_end_dvs2=', ear_end_dvs2, &
    '  storage_onset=', storage_onset_dvs2, '  storage_peak=', storage_peak_frac
  IF (kpt_env_min > 0.0_wp) WRITE(*,'(A,F6.2,A,F6.2,A)') ' Cultivar kernel-mass envelope (rec 125 rows 7-8): ', &
    kpt_env_min, ' - ', kpt_env_max, ' g per 1000 kernels'
  ! The grain-filling clock in force on this deck, and where each of its two
  ! numbers came from. Printed on every deck, so the provenance travels with the
  ! run rather than with a document.
  WRITE(*,'(A,F10.6,A,A,A,F11.5,A,A,A)') ' Cultivar grain-fill clock: T_base=', gfb_tbase, &
    ' degC (', TRIM(MERGE('cultivar_yield.csv row 9 ', 'compiled default         ', gfb_tbase_from_file)), &
    ')  D_fill=', gfb_dfill, ' degCd above that base (', &
    TRIM(MERGE('cultivar_yield.csv row 10', 'compiled default         ', gfb_dfill_from_file)), ')'
  ! ================================================================== !
  ! Record 177: APPLY AND ENFORCE THE CULTIVAR GRAIN-NITROGEN POTENTIAL. The same
  ! three lines run for every genotype - row supplied means used, row absent
  ! means the compiled default - and no cultivar name and no site name appears
  ! anywhere in this block.
  !
  ! THE ENFORCED WINDOW is the CONTROL (intact-ear) range of Martre, Porter,
  ! Jamieson & Triboi (2003) Plant Physiol 133:1959-1967 Table I p.1961, column
  ! mg N grain-1, rows Control: Arche 0.80 +/- 0.02, Recital 0.75 +/- 0.02, Renan
  ! 1.11 +/- 0.02, Tamaro 1.02 +/- 0.01 (field, non-limiting nitrogen supply),
  ! read at first hand. A user value outside that range by MORE than the 10 %
  ! margin is refused; a value inside the margin but outside the range runs and
  ! warns. The compiled default is deliberately not range-checked: it is the
  ! widest value of the same table, 1.38 mg N kernel-1 from an ABLATED ear, i.e.
  ! a reported ceiling and not a cultivar trait, and checking it would refuse
  ! every deck that supplies no trait row.
  ! ================================================================== !
  ! Under record 177 = 1 the compiled default for a cultivar that supplies NO
  ! trait row becomes the control MAXIMUM 1.11 mg N kernel-1 of the same table,
  ! i.e. a value inside the enforced range, so that a default lies in its own
  ! documented window. The 1.38 of the ablated ear stays only on the record
  ! 177 = 0 path, which is the released configuration.
  IF (grain_n_potential_trait == 1 .AND. .NOT. (gnp_row_value > 0.0_wp)) THEN
    gnp_mgN_kernel = gnp_range_hi
  END IF
  IF (gnp_row_value > 0.0_wp) THEN
    WRITE(*,'(A,F8.4,A)') ' WARNING rec 177 (grain_n_potential_trait = 0):' // &
      ' input/cultivar_yield.csv supplies grainN_pot_mgN_kernel = ', gnp_row_value, &
      ' mg N kernel-1 and it is IGNORED -- the record is OFF, so the compiled ceiling' // &
      ' gnf_mgN_kernel_hi = 1.38 mg N kernel-1 stands.  Set driver_main.txt rec 177 = 1' // &
      ' to use the cultivar value.'
  END IF
  ! ON-ONLY, so a rec 177 = 0 run.log is unchanged.
END BLOCK
! rec 151: fruiting efficiency as a cultivar trait.
gnx_floret_scale = 1.0_wp
IF (fe_trait_mode == 1) THEN
  IF (fruiting_efficiency < 50.0_wp .OR. fruiting_efficiency > 200.0_wp) THEN
    WRITE (*,'(A,F8.2,A)') ' ERROR rec 151: cultivar fruiting_efficiency = ', fruiting_efficiency, &
      ' grains per g spike DW is outside the measured range 50-200 (Slafer et al. 2015 Fig. 4 p.97).'
    ERROR STOP 'rec 151: fruiting_efficiency outside [50, 200]'
  END IF
  gnx_floret_scale = (gnx_kern_per_fl / fruiting_efficiency) / gnx_m_floret
  WRITE (*,'(A,F7.2,A,F7.3,A,F6.3,A)') ' rec 151 FRUITING-EFFICIENCY TRAIT ON: FE = ', fruiting_efficiency, &
    ' grains g-1 spike DW -> m_floret,cv = ', 1000.0_wp * gnx_kern_per_fl / fruiting_efficiency, &
    ' mg per competent floret (Fischer 10.0); floret-mass scale = ', gnx_floret_scale, &
    ' applied to the floret-survival hazard and the rec-98 floral bound'
END IF
! Record 177 dependency, refused rather than left a silent no-op. The cultivar
! ceiling is read by exactly one expression, the per-kernel grain-nitrogen
! ceiling of record 125 inside growth_module.f90, and that expression is
! entered only when record 125 is 1 or 2.
IF (grain_n_potential_trait == 1 .AND. &
    .NOT. (kernel_potential_law == 1 .OR. kernel_potential_law == 2)) THEN
  WRITE (*,'(A,I0,A)') ' ERROR rec 177: the cultivar grain-N potential is applied ONLY by the' // &
    ' record-125 per-kernel grain-N ceiling, which is entered only at rec 125' // &
    ' kernel_potential_law = 1 or 2; read rec 125 = ', kernel_potential_law, &
    '. Set rec 125 = 1 or 2, or rec 177 = 0.'
  ERROR STOP 'rec 177 requires rec 125 kernel_potential_law = 1 or 2'
END IF
IF (kernel_potential_law >= 2 .AND. kpt_env_min <= 0.0_wp) THEN
  WRITE (*,'(A)') ' ERROR rec 125: modes 2/3 need the cultivar kernel-mass envelope (input/cultivar_yield.csv rows 7-8' // &
    ' TKW_env_min, TKW_env_max); the file has only the six required values'
  ERROR STOP 'rec 125 >= 2 requires cultivar_yield.csv rows 7-8'
END IF

! Read the T2 leaf-partition-window parameters from input/cultivar_partition.csv (optional).
! File absent or partial -> SHARED WOFOST defaults (0.2741, 0.0, 0.95) = incumbent table.
BLOCK
  INTEGER :: unit_cp, ios_cp, ip_cp
  CHARACTER(LEN=512) :: line_cp
  INTEGER, PARAMETER :: CP_PARAMS = 3
  REAL(wp) :: cp_pv(CP_PARAMS)
  cp_pv = 0.0_wp; ip_cp = 0
  OPEN(newunit=unit_cp, FILE='input/cultivar_partition.csv', STATUS='OLD', IOSTAT=ios_cp)
  IF (ios_cp /= 0) THEN
    WRITE(*,*) ' NOTE: input/cultivar_partition.csv absent; leaf window = WOFOST baseline defaults.'
  ELSE
    DO WHILE (ip_cp < CP_PARAMS)
      READ(unit_cp, '(A)', IOSTAT=ios_cp) line_cp
      IF (ios_cp /= 0) EXIT
      line_cp = ADJUSTL(line_cp)
      IF (LEN_TRIM(line_cp) == 0 .OR. line_cp(1:1) == '!') CYCLE
      ip_cp = ip_cp + 1
      READ(line_cp, *, IOSTAT=ios_cp) cp_pv(ip_cp)
    END DO
    CLOSE(unit_cp)
    IF (ip_cp == CP_PARAMS) THEN
      fleaf_elong = cp_pv(1); fleaf_heading = cp_pv(2); fleaf_end_dvs = cp_pv(3)
    ELSE
      WRITE(*,*) ' WARNING: cultivar_partition.csv had ', ip_cp, ' of ', CP_PARAMS, ' values; using WOFOST defaults.'
    END IF
  END IF
  ! Fail-fast validation (applies to read values AND defaults):
  IF (fleaf_elong < 0.0_wp .OR. fleaf_elong > 1.0_wp) THEN
    WRITE(*,*) ' ERROR: fleaf_elong out of [0,1]: ', fleaf_elong; STOP 1
  END IF
  IF (fleaf_heading < 0.0_wp .OR. fleaf_heading > 1.0_wp) THEN
    WRITE(*,*) ' ERROR: fleaf_heading out of [0,1]: ', fleaf_heading; STOP 1
  END IF
  IF (fleaf_end_dvs < 0.95_wp .OR. fleaf_end_dvs > 1.50_wp) THEN
    WRITE(*,*) ' ERROR: fleaf_end_dvs out of [0.95,1.50]: ', fleaf_end_dvs; STOP 1
  END IF
  IF (fleaf_heading > 1.0e-9_wp .AND. fleaf_end_dvs <= 0.95_wp + 1.0e-9_wp) THEN
    WRITE(*,*) ' ERROR: fleaf_heading > 0 requires fleaf_end_dvs > 0.95 (a closure node);'
    WRITE(*,*) '        otherwise leaf allocation would continue through grain fill. Got:', &
      fleaf_heading, fleaf_end_dvs
    STOP 1
  END IF
  WRITE(*,'(A,F6.4,A,F6.4,A,F5.3)') ' Leaf window: fleaf_elong=', fleaf_elong, &
    '  fleaf_heading=', fleaf_heading, '  fleaf_end_dvs=', fleaf_end_dvs
END BLOCK

! Read mineral parameters (Task 8) and initialize labile soil pools.
BLOCK
  INTEGER :: ios_min, ncols_min
  INTEGER :: i_min_g          ! rec 180: Ca-row search counter
  CALL read_mineral_params('input/driver_minerals.csv', ios_min, ncols_min)
  IF (ios_min /= 0) THEN
    WRITE(*,*) ' WARNING: input/driver_minerals.csv not read fully; mineral pools may be zero.'
  ELSE
    WRITE(*,*) ' Mineral module (Task 8): Fe Zn Ca Cu Mg Mn from calcareous soil; loaded.'
  END IF
  ! Record 131 guard: the cascading reader above defaults columns 17 and 18 to
  ! zero when a deck driver_minerals.csv predates the record and carries 16
  ! columns or fewer. With record 131 off that is harmless, since the values are
  ! never read downstream; with it on, the run would complete with every element
  ! calibrated availability silently zero - a positional failure rather than a
  ! refusal. Hard refusal, not a silent fallback.
  IF (mineral_supply_law == 1 .AND. ncols_min < 18) THEN
    WRITE (*,'(A)') ' ERROR rec 131 = 1: input/driver_minerals.csv does not carry the columns'
    WRITE (*,'(A)') '   this mechanism requires.'
    WRITE (*,'(A,I0,A)') '   Found ', ncols_min, ' trailing numeric columns on at least one mineral row;'
    WRITE (*,'(A)') '   rec 131 needs 18: column 17 (avail_A_mgL_CAL, calibrated site-level'
    WRITE (*,'(A)') '   availability, Fe/Zn/Cu/Mn) and column 18 (xylem_Ca_mgL_CAL, calibrated'
    WRITE (*,'(A)') '   Ca xylem concentration).'
    WRITE (*,'(A)') '   Fix: either supply an 18-column input/driver_minerals.csv for this deck'
    WRITE (*,'(A)') '   (see INPUT_GUIDE.md record 131), or set'
    WRITE (*,'(A)') '   rec 131 = 0 in driver_main.txt for this deck.'
    ERROR STOP 'rec 131: driver_minerals.csv missing columns 17/18 (avail_A_mgL_CAL, xylem_Ca_mgL_CAL)'
  END IF
  ! Record 140 guard: column 19, the buffer power, is OPTIONAL by design, and the
  ! mineral reader's own cascading fallback defaults a missing or non-positive
  ! value to the stated 1.0 and never to 0.0. A missing column 19 is therefore
  ! not a silent failure, and this is a warning rather than a refusal.
  IF (mineral_supply_mech >= 1 .AND. ncols_min < 19) THEN
    WRITE (*,'(A,I0,A)') ' NOTE rec 140 = 1: input/driver_minerals.csv carries only ', ncols_min, &
      ' trailing numeric columns (no column 19, buffer_power_b); every element''s buffer power b_e falls' // &
      ' back to the DECLARED default 1.0 (read_mineral_params'' own non-positive/absent fallback).'
  END IF
  ! ================================================================== !
  ! Record 180: the xylem calcium-concentration guard and the on-only echo. Three
  ! things happen here and nothing else:
  !   (i)   a REFUSAL if the deck driver_minerals.csv does not carry the column
  !         the law reads (column 19, the calcium xylem concentration). The
  !         cascading reader defaults a missing column to zero, which under mode 1
  !         would silently deliver zero calcium - a positional failure. Same
  !         shape as the record-131 guard above.
  !   (ii)  ENFORCEMENT of the measured xylem-sap envelope 12.0-661 mg L-1 (White
  !         & Broadley 2003 Ann Bot 92:487 p.490, a xylem calcium concentration
  !         between 300 umol L-1 and 16.5 mmol L-1 has been reported; Atkinson,
  !         Ruiz & Mansfield 1992 J Exp Bot 43:1315 Table 3 p.1319, 0.4 to
  !         15.3 mol m-3) with a 10 % margin, the same rule record 177 applies:
  !         outside the envelope by more than the margin is an error stop, inside
  !         the margin but outside the envelope runs and warns. The value is
  !         never clipped and never replaced by the sourced default.
  !   (iii) the on-only echo.
  ! ================================================================== !
  IF (grain_ca_loading_basis == 1) THEN
    BLOCK
      INTEGER  :: i_ca_g
      REAL(wp) :: cxl_lo_hard, cxl_hi_hard
      i_ca_g = 0
      DO i_min_g = 1, NMIN
        IF (mineral_names(i_min_g) == 'Ca') i_ca_g = i_min_g
      END DO
      IF (i_ca_g == 0) THEN
        WRITE (*,'(A)') ' ERROR rec 180: no Ca row in the mineral element list; the grain' // &
          ' calcium loading law has no element to act on.'
        ERROR STOP 'rec 180: no Ca element row'
      END IF
      IF (ncols_min < 18) THEN
        WRITE (*,'(A)') ' ERROR rec 180 = 1: input/driver_minerals.csv does not carry the column'
        WRITE (*,'(A)') '   this mechanism reads.'
        WRITE (*,'(A,I0,A)') '   Found ', ncols_min, ' trailing numeric columns on at least one mineral row;'
        WRITE (*,'(A)') '   rec 180 needs 18, because c_xyl_Ca IS the existing Ca-row column 19'
        WRITE (*,'(A)') '   (xylem_Ca_mgL_CAL).  The reader defaults a missing column to 0.0 mg L-1,'
        WRITE (*,'(A)') '   which would deliver ZERO calcium on every day -- a silent positional'
        WRITE (*,'(A)') '   failure, not a refusal.  Fix: supply an 18-column driver_minerals.csv'
        WRITE (*,'(A)') '   for this deck (INPUT_GUIDE.md record 180),'
        WRITE (*,'(A)') '   or set rec 180 = 0.'
        ERROR STOP 'rec 180: driver_minerals.csv missing column 19 (xylem_Ca_mgL_CAL = c_xyl_Ca)'
      END IF
      gcl_cxyl_used = mm_xylem_ca_cal(i_ca_g)
      cxl_lo_hard = gcl_cxyl_lo * (1.0_wp - gcl_cxyl_margin)
      cxl_hi_hard = gcl_cxyl_hi * (1.0_wp + gcl_cxyl_margin)
      IF (gcl_cxyl_used < cxl_lo_hard .OR. gcl_cxyl_used > cxl_hi_hard) THEN
        WRITE (*,'(A,F12.4,A,F7.1,A,F7.1,A,F7.1,A,F7.1,A)') &
          ' ERROR rec 180 (grain_ca_loading_basis): the site calibration variable c_xyl_Ca' // &
          ' (driver_minerals.csv Ca row, column 19 xylem_Ca_mgL_CAL) = ', gcl_cxyl_used, &
          ' mg L-1 is outside the MEASURED xylem-sap envelope [', gcl_cxyl_lo, ', ', gcl_cxyl_hi, &
          '] mg L-1 by more than the 10 % margin [', cxl_lo_hard, ', ', cxl_hi_hard, ']'
        WRITE (*,'(A)') '        ENVELOPE SOURCE: White PJ & Broadley MR (2003) Ann Bot' // &
          ' 92:487-511, DOI 10.1093/aob/mcg164, p.490: "[Ca]xylem between 300 uM and 16.5 mM' // &
          ' has been reported"; independently Atkinson, Ruiz & Mansfield (1992) J Exp Bot' // &
          ' 43:1315-1324, DOI 10.1093/jxb/43.10.1315, Table 3 p.1319, 0.4-15.3 mol m-3.'
        WRITE (*,'(A)') '        The record-166-era values 842.0246 and 654.339 mg L-1 are' // &
          ' WITHDRAWN by record 180: both are linear rescales of the calcite bracket start' // &
          ' 54.9 mg L-1 onto one deck''s OWN OBSERVED grain Ca, i.e. the campaign observation' // &
          ' written back in as an input, and 842.0246 mg L-1 = 21.0 mmol L-1 is 1.27x the' // &
          ' measured upper bound.  The sourced default is 15.0 mg L-1 (SHARED.f90).'
        WRITE (*,'(A)') '        The value is NOT clipped and NOT replaced by the default:' // &
          ' a refused calibration record does not run (SHARED.f90).'
        ERROR STOP 'rec 180: c_xyl_Ca outside the measured xylem-sap envelope'
      END IF
      IF (gcl_cxyl_used < gcl_cxyl_lo .OR. gcl_cxyl_used > gcl_cxyl_hi) THEN
        WRITE (*,'(A,F12.4,A,F7.1,A,F7.1,A)') &
          ' WARNING rec 180 (grain_ca_loading_basis): c_xyl_Ca = ', gcl_cxyl_used, &
          ' mg L-1 is OUTSIDE the measured xylem-sap envelope [', gcl_cxyl_lo, ', ', gcl_cxyl_hi, &
          '] mg L-1 but inside the 10 % margin, so the run PROCEEDS.  White & Broadley (2003)' // &
          ' Ann Bot 92:487 p.490; Atkinson (1992) J Exp Bot 43:1315 Table 3 p.1319.' // &
          '  This is a calibration value outside its own documented range.'
      END IF
      WRITE (*,'(A)') ' rec 180 grain_ca_loading_basis = 1 : GRAIN CALCIUM NOW ARRIVES IN THE' // &
        ' TRANSPIRATION STREAM.  dCa_grain/dt = f_grain_post x c_xyl_Ca x Tr_act and' // &
        ' dCa_chaff/dt = (f_head_post - f_grain_post) x c_xyl_Ca x Tr_act from anthesis to' // &
        ' maturity, remainder to the vegetative terminus, nothing ever remobilized' // &
        ' (SHARED.f90).'
      WRITE (*,'(A,F8.5,A,F8.5,A)') '        f_head_post = ', gcl_f_head_post, &
        ' and f_grain_post = ', gcl_f_grain_post, &
        ' : Hocking PJ (1994) J Plant Nutr 17:1289-1308, DOI 10.1080/01904169409364807,' // &
        ' Tables 1(b)+1(c) p.1294 with Tables 2+3 p.1295 (SHARED.f90).'
      WRITE (*,'(A,F12.4,A,F7.1,A,F7.1,A,F6.1,A)') '        c_xyl_Ca in force = ', gcl_cxyl_used, &
        ' mg L-1, inside the measured envelope [', gcl_cxyl_lo, ', ', gcl_cxyl_hi, &
        '] mg L-1; sourced default ', gcl_cxyl_ref, ' mg L-1 (SHARED.f90).'
      WRITE (*,'(A,F6.1,A,F6.1,A)') '        PRE-REGISTERED eCO2 SIGN for grain Ca:' // &
        ' NEGATIVE, magnitude ', gcl_eco2_lo_pct, ' to ', gcl_eco2_hi_pct, &
        ' % -- Loladze I (2014) eLife 3:e02245, DOI 10.7554/eLife.02245, p.6 lines 153-155,' // &
        ' Fig. 2, pooled over m = 117 means.  CARRIED BESIDE IT, SAME SIGN, NOT SIGNIFICANT:' // &
        ' Hogy P et al. (2009) Plant Biol 11(Suppl 1):60-69, DOI' // &
        ' 10.1111/j.1438-8677.2009.00230.x, the HOHENHEIM mini-FACE (cv. TRISO, 2004-2006,' // &
        ' ambient + 150 umol mol-1), Fig. 2 p.63 and its Results text p.63: "several' // &
        ' macro-elements such as Ca, P and S ... were not significantly lowered in the' // &
        ' high-CO2 treatment", Discussion p.67 "Although not significant, concentrations of' // &
        ' other minerals, except for Na and Cr, were slightly decreased".  A negative model' // &
        ' response smaller than Loladze''s band is consistent with Hogy and is not by itself' // &
        ' a failure; the SIGN plus the WDF25 campaign is the gate (SHARED.f90).'
      IF (grain_ca_xylem_throughput == 1) THEN
        WRITE (*,'(A)') ' NOTE rec 180 = 1 SUPERSEDES rec 166 FOR CALCIUM: the record-166' // &
          ' ear-area share and the Lizana et al. (2010) grain net-water import are not' // &
          ' evaluated on the calcium route while rec 180 = 1 (the grain''s own water import' // &
          ' is phloem and cannot carry a phloem-immobile element -- Zee & O''Brien 1970 Aust J' // &
          ' Biol Sci 23:783 p.783-784; SHARED.f90).  Record 166 is left at its own' // &
          ' value and governs nothing else.'
      END IF
    END BLOCK
  END IF
  ! ================================================================== !
  ! Record 181: the magnesium solution-basis guard and its enforced window. Mode 1
  ! promotes column 17 of the magnesium row to the soil-solution magnesium
  ! concentration, so (a) the column must EXIST, the same positional refusal
  ! record 131 makes above, never a silent zero concentration, and (b) the value
  ! must lie inside the measured range of SHARED.f90 under the same 10 % margin
  ! rule record 177 applies to the grain-nitrogen trait: outside by more than the
  ! margin is an error stop, inside the margin but outside the range is a
  ! warning. A pre-registered gate is never widened to admit a model value.
  ! At mode 0 this block does not execute and prints nothing.
  ! ================================================================== !
  ! ================================================================== !
  ! Record 183: the buffer-power input guard and the per-element table. Mode 1
  ! builds the buffer power as the water content plus the bulk density times the
  ! distribution coefficient, so (a) column 3 must EXIST and (b) every modeled
  ! non-calcium element must have SOME linear partition to invert - a positive
  ! column-3 coefficient, or, for iron, a positive column-4 ceiling against a
  ! positive column-2 pool - otherwise an error stop naming the element, and
  ! never a silent fallback to the unsourced 1.0 this record exists to replace.
  ! The table printed here is the per-element buffer power at this deck's own
  ! measured bulk density and at a reference water content. At mode 0 nothing
  ! executes and nothing prints.
  ! ================================================================== !
  IF (mineral_buffer_power_basis == 1) THEN
    IF (ncols_min < 3) THEN
      WRITE (*,'(A,I0,A)') ' ERROR rec 183 = 1: input/driver_minerals.csv carries only ', ncols_min, &
        ' trailing numeric columns; column 3 (kd_L_kg) is the solid-solution partition the buffer' // &
        ' power b_e = theta + rho_b*Kd_e is built from. Supply a >=3-column file (INPUT_GUIDE.md' // &
        ' record 183) or set rec 183 = 0.'
      ERROR STOP 'rec 183: driver_minerals.csv missing column 3 (kd_L_kg)'
    END IF
    IF (nud_bulk_density > 0.0_wp) THEN
      bpb_rho_b_used = nud_bulk_density
    ELSE
      bpb_rho_b_used = 978.0_wp / 600.0_wp     ! = soil_mass_kg_m2/(mm_depth_m*1000) = 1.63 kg L-1, the module's OWN stated basis
      bpb_rho_fallback = .TRUE.
      WRITE (*,'(A)') ' NOTE rec 183 = 1: this deck supplies no driver_N.csv position 70' // &
        ' (bulk_density_topsoil_Mg_m3), so the buffer power uses the mineral module''s OWN' // &
        ' soil-mass basis rho_b = 978/600 = 1.6300 kg L-1 (BD 1.63 g cm-3 x 0.6 m, the SAME' // &
        ' conversion the labile pool is built on). DECLARED, never silent; SHARED.f90.'
    END IF
    DO imin_loop = 1, NMIN
      IF (imin_loop == 3) CYCLE                ! 3 = idx_Ca: excluded from the diffusion term by its own gate
      IF (mm_modelled(imin_loop) == 1) THEN
        bpb_kd_used(imin_loop) = mm_kd(imin_loop)
        IF (bpb_kd_used(imin_loop) <= 0.0_wp .AND. mm_ceiling(imin_loop) > 0.0_wp &
            .AND. mm_pool_mgkg(imin_loop) > 0.0_wp) &
          bpb_kd_used(imin_loop) = mm_pool_mgkg(imin_loop) / mm_ceiling(imin_loop)
        IF (bpb_kd_used(imin_loop) <= 0.0_wp) THEN
          WRITE (*,'(A,A,A)') ' ERROR rec 183 = 1: element ', TRIM(mineral_names(imin_loop)), &
            ' has neither a positive driver_minerals.csv column-3 Kd nor a positive column-4' // &
            ' ceiling against a positive column-2 pool, so no linear solid-solution partition' // &
            ' can be formed and b_e = theta + rho_b*Kd_e is undefined. Supply a Kd for that row' // &
            ' or set rec 183 = 0 (SHARED.f90).'
          ERROR STOP 'rec 183: no Kd available for a modelled non-Ca element'
        END IF
        WRITE (*,'(A,A,A,A,F12.4,A,F8.4,A,F12.2,A,F12.2,A,A)') ' rec 183 mineral_buffer_power_basis = 1:', &
          ' element ', TRIM(mineral_names(imin_loop)), '  Kd = ', bpb_kd_used(imin_loop), &
          ' L kg-1   rho_b = ', bpb_rho_b_used, ' kg L-1   b_e(theta=0.15) = ', &
          0.15_wp + bpb_rho_b_used * bpb_kd_used(imin_loop), '   b_e(theta=0.25) = ', &
          0.25_wp + bpb_rho_b_used * bpb_kd_used(imin_loop), &
          '   [b_e = theta + rho_b*Kd_e; Baldwin/Nye/Tinker 1973 p.622 + p.625;' // &
          ' Roose/Fowler 2004 Eq.(2.8) p.175; SHARED.f90]'
      END IF
    END DO
    WRITE (*,'(A)') ' rec 183 mineral_buffer_power_basis = 1: the driver_minerals.csv column-19' // &
      ' constant buffer_power_b (1.0 on every shipped row, DECLARED UNSOURCED by record 140' // &
      ') is NO LONGER READ by the diffusion coefficient; b_e is computed daily per' // &
      ' element from the model''s own vwc_root, this deck''s measured rho_b and the element''s' // &
      ' own Kd. Iron''s Kd is the OPERATIONAL pool/ceiling inversion 6.7/0.119 = 56.3025 L kg-1' // &
      ' (SHARED.f90); calcium is excluded from the diffusion term and is untouched.'
  END IF
  CALL init_mineral_pools()
END BLOCK

! -------------------------------------------------------------------- !
! Read the fertilizer schedule from input/N_fertilizer.csv.
! Format: the number of applications on the first non-comment line, then one
! row per application, days after sowing and nitrogen in kg ha-1. Nitrogen is
! converted from kg ha-1 to g m-2.
! -------------------------------------------------------------------- !
BLOCK
  INTEGER :: unit_N_fert, ios_f, ia, n_apps_f, das_f
  REAL(wp) :: nkg_f
  CHARACTER(LEN=256) :: line_f
  OPEN(newunit=unit_N_fert, FILE='input/N_fertilizer.csv', STATUS='OLD', IOSTAT=ios_f)
  IF (ios_f /= 0) THEN
    WRITE(*,*) 'ERROR: Could not open input/N_fertilizer.csv'; STOP
  END IF
  ! Skip comment lines; read n_apps from first data line.
  n_apps_f = 0
  DO
    READ(unit_N_fert, '(A)', IOSTAT=ios_f) line_f
    IF (ios_f /= 0) EXIT
    line_f = ADJUSTL(line_f)
    IF (LEN_TRIM(line_f) == 0 .OR. line_f(1:1) == '!') CYCLE
    READ(line_f, *, IOSTAT=ios_f) n_apps_f
    EXIT
  END DO
  IF (n_apps_f <= 0) THEN
    WRITE(*,*) 'ERROR: N_fertilizer.csv has no applications (n_apps <= 0)'; STOP
  END IF
  N_fert_n_apps = n_apps_f
  ALLOCATE(N_fert_DAS_arr(N_fert_n_apps), N_fert_N_arr(N_fert_n_apps))
  ALLOCATE(N_fert_fNO3_arr(N_fert_n_apps), N_fert_fNH4_arr(N_fert_n_apps), N_fert_furea_arr(N_fert_n_apps))
  BLOCK
    INTEGER :: ios_form, das_tmp
    REAL(wp) :: nkg_tmp
    CHARACTER(LEN=16) :: form_str
    DO ia = 1, N_fert_n_apps
      ! Skip comment/blank lines within the application list.
      DO
        READ(unit_N_fert, '(A)', IOSTAT=ios_f) line_f
        IF (ios_f /= 0) EXIT
        line_f = ADJUSTL(line_f)
        IF (LEN_TRIM(line_f) == 0 .OR. line_f(1:1) == '!') CYCLE
        EXIT
      END DO
      IF (ios_f /= 0) THEN
        WRITE(*,*) 'ERROR: Could not read N_fertilizer.csv row ', ia; STOP
      END IF
      READ(line_f, *, IOSTAT=ios_f) das_f, nkg_f          ! always: DAS, N_kg_per_ha
      IF (ios_f /= 0) THEN
        WRITE(*,*) 'ERROR: Could not parse N_fertilizer.csv row ', ia; STOP
      END IF
      form_str = 'urea'                                   ! backward-compatible default
      READ(line_f, *, IOSTAT=ios_form) das_tmp, nkg_tmp, form_str   ! optional 3rd field: form keyword
      IF (ios_form /= 0) form_str = 'urea'
      N_fert_DAS_arr(ia) = das_f
      N_fert_N_arr(ia)   = nkg_f / 10.0_wp                 ! kg/ha -> g/m2 (1 kg/ha = 0.1 g/m2)
      ! Map fertilizer form -> (NO3, NH4, urea) fractions.
      SELECT CASE (TRIM(ADJUSTL(form_str)))
        CASE ('AN', 'an', 'AmmoniumNitrate'); N_fert_fNO3_arr(ia)=0.5_wp; N_fert_fNH4_arr(ia)=0.5_wp; N_fert_furea_arr(ia)=0.0_wp
        CASE ('CAN', 'can');                  N_fert_fNO3_arr(ia)=0.5_wp; N_fert_fNH4_arr(ia)=0.5_wp; N_fert_furea_arr(ia)=0.0_wp
        CASE ('AS', 'as', 'AmmoniumSulfate'); N_fert_fNO3_arr(ia)=0.0_wp; N_fert_fNH4_arr(ia)=1.0_wp; N_fert_furea_arr(ia)=0.0_wp
        CASE ('KNO3', 'kno3', 'nitrate', 'NO3'); N_fert_fNO3_arr(ia)=1.0_wp; N_fert_fNH4_arr(ia)=0.0_wp; N_fert_furea_arr(ia)=0.0_wp
        CASE DEFAULT  ! 'urea' and unrecognized -> urea
          N_fert_fNO3_arr(ia)=0.0_wp; N_fert_fNH4_arr(ia)=0.0_wp; N_fert_furea_arr(ia)=1.0_wp
      END SELECT
    END DO
  END BLOCK
  CLOSE(unit_N_fert)
  WRITE(*,'(A,I0,A,F6.1,A)') ' N fertilizer: ', N_fert_n_apps, &
    ' applications, total ', SUM(N_fert_N_arr) * 10.0_wp, ' kg N/ha (default form: urea)'
END BLOCK

! Soil mineral N pool initialization. NO3 = N_soil_init; NH4 = N_soil_init_NH4 (default 0).
N_soil_min = N_soil_init
N_NH4_soil = N_soil_init_NH4
N_urea_soil = 0.0_wp
N_soil_uptake_total = 0.0_wp
N_soil_leach_total  = 0.0_wp
N_fert_applied = 0

! Record 93 at or above 1: initialize the humus organic-nitrogen pool and
! evaluate its potential rate constant once from the soil analysis. Brisson et
! al. (1998) Eq. (33) p.342:
!   K2POT = 0.21 / ((12.5 + clay)(54.5 + calcium carbonate)).
! The clay, the calcium carbonate and the initial humus nitrogen are MEASURED
! PER-DECK INPUTS at driver_N.csv positions 65, 66 and 67, cited by file and
! field in SHARED.f90, and not one of them may be adjusted to improve an
! output. Under record 93 = 0 every one of these stays zero and nothing
! downstream reads them.
IF (soil_n_mech >= 1) THEN
  ! Hard refusal, not a silent fallback. A deck that does not carry its own
  ! measured soil analysis must not mineralize on somebody else's soil.
  IF (stics_argi_pct < 0.0_wp .OR. stics_calc_pct < 0.0_wp .OR. &
      stics_nhum0_gN_m2 < 0.0_wp) THEN
    WRITE (*,'(A)') ' ERROR rec 93 >= 1: this deck carries no MEASURED soil analysis.'
    WRITE (*,'(A)') '   input/driver_N.csv must supply positions 65/66/67 ='
    WRITE (*,'(A)') '   soil_clay_pct (ARGI), soil_calc_pct (CALC), soil_humusN_gN_m2 (NHUM(0)).'
    WRITE (*,'(A)') '   Provenance per deck: input/soil_analysis_measured.csv, SHARED.f90.'
    WRITE (*,'(A,3ES13.5)') '   read: ', stics_argi_pct, stics_calc_pct, stics_nhum0_gN_m2
    ERROR STOP 'rec 93: measured soil analysis missing (driver_N.csv pos 65-67)'
  END IF
  NHUM_soil      = stics_nhum0_gN_m2
  NHUM_soil_init = stics_nhum0_gN_m2
  humus_k2pot    = stics_k2pot_num / &
                   ((stics_k2pot_argi0 + stics_argi_pct) * (stics_k2pot_calc0 + stics_calc_pct))
  WRITE (*,'(A,I0,A,ES13.6,A,F8.3,A)') ' rec 93 soil_n_mech = ', soil_n_mech, &
    ' : humus mineralisation ON (Brisson et al. 1998 formulation). K2POT = ', humus_k2pot, &
    ' d-1, NHUM(0) = ', NHUM_soil, ' g N m-2 (N_mineral_rate NOT read).'
  WRITE (*,'(A,F7.3,A,F7.3,A)') '   MEASURED soil analysis: ARGI = ', stics_argi_pct, &
    ' % clay, CALC = ', stics_calc_pct, ' % CaCO3 (driver_N pos 65/66/67).'
END IF

! -------------------------------------------------------------------- !
! Read the organic-nitrogen and residue treatments from
! input/N_organic_treatments.csv. An absent file, or a count of zero, leaves
! the treatments inactive and reproduces the no-residue baseline.
! -------------------------------------------------------------------- !
BLOCK
  INTEGER :: unit_org, ios_org, n_org, io_i
  CHARACTER(LEN=512) :: line_org
  CHARACTER(LEN=24)  :: ttype
  INTEGER  :: das_o, age_o
  REAL(wp) :: mass_o, npct_o, cpct_o, fast_o, cn_def, fast_def
  org_treatments_active = .FALSE.; N_org_n = 0
  OPEN(newunit=unit_org, FILE='input/N_organic_treatments.csv', STATUS='OLD', IOSTAT=ios_org)
  IF (ios_org == 0) THEN
    n_org = -1
    DO   ! first non-comment line = count
      READ(unit_org,'(A)',IOSTAT=ios_org) line_org
      IF (ios_org /= 0) EXIT
      line_org = ADJUSTL(line_org)
      IF (LEN_TRIM(line_org) == 0 .OR. line_org(1:1) == '!') CYCLE
      READ(line_org,*,IOSTAT=ios_org) n_org
      EXIT
    END DO
    IF (n_org > 0) THEN
      N_org_n = n_org
      ALLOCATE(N_org_das(n_org), N_org_age(n_org), N_org_mass(n_org), N_org_Npct(n_org), &
               N_org_Cpct(n_org), N_org_fast(n_org), N_org_CN(n_org))
      DO io_i = 1, n_org
        DO
          READ(unit_org,'(A)',IOSTAT=ios_org) line_org
          IF (ios_org /= 0) EXIT
          line_org = ADJUSTL(line_org)
          IF (LEN_TRIM(line_org) == 0 .OR. line_org(1:1) == '!') CYCLE
          EXIT
        END DO
        IF (ios_org /= 0) THEN
          WRITE(*,*) ' ERROR: N_organic_treatments.csv truncated at row ', io_i; STOP
        END IF
        READ(line_org,*,IOSTAT=ios_org) ttype, das_o, mass_o, npct_o, cpct_o, fast_o, age_o
        IF (ios_org /= 0) THEN
          WRITE(*,*) ' ERROR: could not parse N_organic_treatments.csv row ', io_i; STOP
        END IF
        ! Default C:N and fast_fraction by treatment_type.
        SELECT CASE (TRIM(ADJUSTL(ttype)))
          CASE ('cereal_residue');  cn_def=80.0_wp; fast_def=0.20_wp
          CASE ('immature_cereal'); cn_def=40.0_wp; fast_def=0.35_wp
          CASE ('green_manure');    cn_def=20.0_wp; fast_def=0.55_wp
          CASE ('legume_residue');  cn_def=20.0_wp; fast_def=0.55_wp
          CASE ('manure');          cn_def=12.0_wp; fast_def=0.40_wp
          CASE ('compost');         cn_def=18.0_wp; fast_def=0.10_wp
          CASE DEFAULT;             cn_def=30.0_wp; fast_def=0.30_wp   ! custom/unrecognized
        END SELECT
        N_org_das(io_i)  = das_o
        N_org_age(io_i)  = age_o
        N_org_mass(io_i) = mass_o
        N_org_Npct(io_i) = npct_o
        N_org_Cpct(io_i) = cpct_o
        N_org_fast(io_i) = MERGE(fast_def, fast_o, fast_o < 0.0_wp)
        N_org_CN(io_i)   = cn_def
      END DO
      org_treatments_active = .TRUE.
    END IF
    CLOSE(unit_org)
  END IF
  IF (org_treatments_active) THEN
    WRITE(*,'(A,I0,A)') ' Organic treatments: ', N_org_n, ' event(s) loaded.'
  ELSE
    WRITE(*,*) ' Organic treatments: none (no-residue baseline).'
  END IF
END BLOCK

! Pre-sowing residue decomposition, for events with an age before sowing: the
! residue is pre-decomposed at initialization at the configured soil
! temperature and water factor (decomposition, mineralization, immobilization,
! hydrolysis and nitrification, with no uptake and no leaching), which seeds
! the sowing-day ammonium and nitrate pools.
CALL check_event_reachability(); pre_sow_N_resid = 0.0_wp   ! refusal, see SHARED.f90
IF (org_treatments_active) THEN
  BLOCK
    INTEGER :: io_i, d_ps
    REAL(wp) :: nin, soilN0, dC, dN, rspC, nmin, imNH4, imNO3, tovN, khy, kni
    soilN0 = N_urea_soil + N_NH4_soil + N_soil_min + res_fast_N + res_slow_N + soil_mic_N
    DO io_i = 1, N_org_n
      IF (N_org_age(io_i) > 0) THEN
        CALL apply_residue_event(io_i, nin)
        soilN0 = soilN0 + nin                     ! residue N added pre-sowing
        DO d_ps = 1, N_org_age(io_i)
          CALL decompose_residue_day(presow_Tsoil, presow_fwater, dC, dN, rspC, nmin, imNH4, imNO3, tovN)
          ! hydrolysis + nitrification at pre-season temperature (no urea pre-sowing normally, but general)
          khy = k_hyd_ref * Q10_hyd ** ((presow_Tsoil - T_ref_hyd) / 10.0_wp)
          kni = k_nit_ref * Q10_nit ** ((presow_Tsoil - T_ref_nit) / 10.0_wp)
          BLOCK
            REAL(wp) :: dU, dNit
            dU   = MIN(N_urea_soil, MAX(0.0_wp, khy) * N_urea_soil)
            N_urea_soil = N_urea_soil - dU; N_NH4_soil = N_NH4_soil + dU
            dNit = MIN(N_NH4_soil, MAX(0.0_wp, kni) * N_NH4_soil)
            N_NH4_soil = N_NH4_soil - dNit; N_soil_min = N_soil_min + dNit
          END BLOCK
        END DO
      END IF
    END DO
    ! Pre-sowing closure: no uptake/leaching, all transfers internal -> residual ~0.
    pre_sow_N_resid = (N_urea_soil + N_NH4_soil + N_soil_min + res_fast_N + res_slow_N + soil_mic_N) - soilN0
    WRITE(*,'(A,G12.4,A)') ' Pre-sowing soil-N closure residual: ', pre_sow_N_resid, ' g N m-2'
  END BLOCK
END IF

! -------------------------------------------------------------------- !
! LAYERED soil-nitrogen initialization: distribute the post-spin-up initial
! mineral nitrogen over the Richards layers at UNIFORM SOLUTION CONCENTRATION,
! i.e. each layer share proportional to its water content times its thickness.
! The ammonium initial condition, default zero, is distributed by the same
! rule for generality. Urea stays the layer-1 scalar and the surface
! ammoniacal pool is unchanged. The scalar pools are retained as EXACT
! AGGREGATE MIRRORS, updated after every layered operation, so that every
! downstream diagnostic keeps its meaning.
! -------------------------------------------------------------------- !
IF (soilN_layered == 1) THEN
  BLOCK
    REAL(wp) :: wsum
    ! ================================================================ !
    ! Record 153 n_init_placement. Mode 0 is the whole-profile share below. Mode 1
    ! restricts the support to the deck's own reported sampling depth (optional row
    ! 71 of driver_N.csv) with the SAME uniform-solution-concentration rule,
    ! weighted by water content times thickness WITHIN that depth, so a layer
    ! straddling the sampling depth contributes only its sampled fraction of
    ! thickness. Mass is conserved by construction: the weights are normalized and
    ! the aggregate mirrors are re-formed from the layer sums below. A sentinel,
    ! non-positive or over-deep reported depth leaves the weights equal to the
    ! whole-profile weights, so mode 1 is then arithmetically identical to mode 0 -
    ! the record never invents a depth.
    ! ================================================================ !
    IF (n_init_placement == 1) THEN
      BLOCK
        REAL(wp) :: w_nip(NSoil_layers), dsamp, frac_l, z_top, z_bot
        INTEGER  :: l_nip, n_lay_used
        dsamp = n_init_depth_m
        IF (dsamp <= 0.0_wp .OR. dsamp >= depth_lim) dsamp = depth_lim
        w_nip = 0.0_wp; n_lay_used = 0
        DO l_nip = 1, NSoil_layers
          z_top = ABS(zb(l_nip))
          z_bot = ABS(zb(l_nip+1))
          IF (z_top >= dsamp - 1.0e-12_wp) CYCLE      ! layer top at or below the sampling depth
          frac_l = (MIN(z_bot, dsamp) - z_top) / MAX(1.0e-12_wp, z_bot - z_top)
          frac_l = MIN(1.0_wp, MAX(0.0_wp, frac_l))
          w_nip(l_nip) = theta(l_nip) * tk(l_nip) * frac_l
          IF (w_nip(l_nip) > 0.0_wp) n_lay_used = n_lay_used + 1
        END DO
        wsum = SUM(w_nip)
        IF (wsum > 1.0e-9_wp) THEN
          NO3_lay = N_soil_min * w_nip / wsum
          NH4_lay = N_NH4_soil * w_nip / wsum
        ELSE
          ! No sampled water volume at all: fall back to the whole-profile rule rather
          ! than lose the reported mass. A refusal is not available here, because the
          ! pool is already reported and must be placed somewhere.
          wsum = SUM(theta * tk)
          IF (wsum > 1.0e-9_wp) THEN
            NO3_lay = N_soil_min * theta * tk / wsum
            NH4_lay = N_NH4_soil * theta * tk / wsum
          ELSE
            NO3_lay = N_soil_min / REAL(NSoil_layers, wp)
            NH4_lay = N_NH4_soil / REAL(NSoil_layers, wp)
          END IF
          WRITE (*,'(A)') ' rec 153 n_init_placement = 1: sampled water volume is zero;' // &
            ' whole-profile placement used (mode-0 rule)'
        END IF
        ! ON-only banner (mode 0 stdout is unchanged, line for line).
        WRITE (*,'(A,F6.3,A,I0,A,I0,A,F9.4,A)') ' rec 153 n_init_placement = 1:' // &
          ' initial mineral N placed on the declared sampling depth ', dsamp, &
          ' m -> layers 1-', n_lay_used, ' of ', NSoil_layers, &
          ' (NO3 = ', SUM(NO3_lay), ' g N m-2)'
        IF (n_init_depth_m <= 0.0_wp) WRITE (*,'(A)') ' rec 153: driver_N.csv row 71' // &
          ' N_soil_init_depth_m ABSENT or <= 0 -> whole profile (identical to mode 0)'
      END BLOCK
    ELSE
    wsum = SUM(theta * tk)
    IF (wsum > 1.0e-9_wp) THEN
      NO3_lay = N_soil_min * theta * tk / wsum
      NH4_lay = N_NH4_soil * theta * tk / wsum
    ELSE
      NO3_lay = N_soil_min / REAL(NSoil_layers, wp)
      NH4_lay = N_NH4_soil / REAL(NSoil_layers, wp)
    END IF
    END IF
    ! Mirrors: force exact aggregate identity from the start.
    N_soil_min = SUM(NO3_lay)
    N_NH4_soil = SUM(NH4_lay)
  END BLOCK
END IF

IF (energy_balance_mode == 2) THEN
  OPEN (newunit=unit_flag_diag, file='output/flag_leaf_diag.csv', status='replace', action='write')
  WRITE (unit_flag_diag,'(A)') 'day,j,hour,n0_Narea_gNm2,An_umol_m2_s,gsw_mol_m2_s,CiCa,Tleaf_C,Ds_Pa,Rd_umol_m2_s'
  flag_diag_open = .TRUE.
END IF
! rec 87 ON-only diagnostic. Opened only when rec 87 = 1, so
! the OFF arm's output directory is unchanged, file for file.
IF (leafN_cap_route == 1) THEN
  OPEN (newunit=unit_leafn_cap, file='output/daily_leafN_capacity.csv', status='replace', action='write')
  WRITE (unit_leafn_cap,'(A)') 'day,lai_assim_m2_m2,W_green_g_m2,N_leaf_before_gN_m2,' // &
    'N_excess_exported_gN_m2,N_surplus_to_mobile_gN_m2,N_leaf_after_gN_m2,' // &
    'Narea_assim_basis_gN_m2,leafN_pct_DM,N_mobile_veg_gN_m2'
  leafn_cap_open = .TRUE.
END IF
! A census file for the intercellular-CO2 diagnostic, written to its own file
! and read by no gate.
! ==================================================================== !
! Record 116, the limiting-limb census: an ON-ONLY file. With the record at 0
! or absent the file is not opened and does not exist, exactly as the
! record-87 and record-115 diagnostics behave, so the output directory of the
! released configuration is unchanged file for file. It is a model internal:
! there is no measured limiting limb in any deck and this file is never placed
! beside an observation. One row per day, time step and leaf class, flushed
! once per committed step by CANOPY_MODEL and joinable to
! output/diurnal_energy_balance.csv on day, step and hour.
! ==================================================================== !
OPEN (newunit=unit_diurnal,file='output/diurnal_energy_balance.csv',status='replace',action='write')
diurnal_header = 'day,gauss_idx,iter,solar_hour,Tair_C,Tleaf_C,VPD_Pa,Rn_W_m2,' // &
  'LE_W_m2,LEc_W_m2,LEs_W_m2,EB_err_W_m2,Hc_W_m2,Hs_W_m2,H_W_m2,' // &
  'gsw_An_mol_m2_s,' // &
  'rca_s_m_s,rcs_s_m_s,rsa_s_m_s,rss_s_m_s,ac_W_m2,as_W_m2,VpdMcf_Pa,' // &
  'An_ground_umol_m2_s,A_ground_umol_m2_s,N_assim_umol_m2_s,' // &
  'Ci_umol_mol,sunlit_LAI,shaded_LAI,LAI_m2_m2'
WRITE (unit_diurnal,'(A)') TRIM(diurnal_header)
! Hourly three-node energy-balance closure diagnostics (energy_balance_mode>=1).
IF (energy_balance_mode >= 1) THEN
  OPEN (newunit=unit_eb_diag,file='output/hourly_energy_balance.csv',status='replace',action='write')
  OPEN (newunit=unit_node_pair,file='output/hourly_node_pairing.csv',status='replace',action='write')
  node_pair_open = .TRUE.
  WRITE (unit_node_pair,'(A)') 'day,hour_idx,node_armed_C,node_committed_C'
  OPEN (newunit=unit_esoil_pair,file='output/hourly_esoil_pairing.csv',status='replace',action='write')
  esoil_pair_open = .TRUE.
  WRITE (unit_esoil_pair,'(A)') 'day,hour_idx,sent_mm,delivered_mm,LEs_committed_W_m2,LEs_pot_endstate_W_m2,' // &
    'etr_sent_mm,etr_delivered_mm,LEc_tr_committed_W_m2'
  OPEN (newunit=unit_eb_demand,file='output/daily_eb_demand.csv',status='replace',action='write')
  WRITE (unit_eb_demand,'(A)') 'day,Rn_demand_MJ,G_demand_MJ,LEc_demand_MJ,LEs_demand_MJ,Hc_demand_MJ,Hs_demand_MJ'
  eb_diag_open = .TRUE.
  WRITE (unit_eb_diag,'(A)') 'day,gauss_idx,energy_balance_mode,solar_hour,Tair_C,Tac_C,Tcan_C,' // &
    'Tsoil_surf_C,VPD_ref_Pa,VPD_can_air_Pa,Ds_leaf_assim_Pa,e_ac_Pa,esat_Tcan_Pa,' // &
    'An_ground_umol_m2_s,gsw_An_mol_m2_s,rcs_s_m,raa_h_s_m,raa_neutral_s_m,rca_s_m,rsa_s_m,rss_s_m,' // &
    'L_obukhov_m,zeta_zoverL,psi_h,psi_m,ustar_m_s,stab_regime,' // &
    'Rnc_W_m2,Rns_W_m2,Rn_total_W_m2,LEc_W_m2,LEs_W_m2,Hc_W_m2,Hs_W_m2,H_total_W_m2,G_W_m2,' // &
    'Tac_minus_Tair_C,Tcan_minus_Tair_C,' // &
    'F_can_W_m2,F_soil_W_m2,F_air_W_m2,F_vap_W_m2,h_soil,inner_iters,outer_iters,conv_flag,' // &
    'Tac_demand_C,Tcan_demand_C,Tsoil_demand_C,Rnc_demand,Rns_demand,Rn_demand,LEc_demand,' // &
    'LEs_demand,Hc_demand,Hs_demand,H_total_demand,G_demand,Tac_minus_Tair_demand,Tcan_minus_Tair_demand'
  ! daily actual-flux (post-hydrology repartition) diagnostics.
  OPEN (newunit=unit_eb_act,file='output/daily_eb_actual.csv',status='replace',action='write')
  eb_act_open = .TRUE.
  WRITE (unit_eb_act,'(A)') 'day,Ir_day_mm,Eintercept_act_mm,interception_drip_mm,' // &
    'LEc_trans_pot_mm,LEc_trans_act_mm,LEc_int_pot_mm,LEc_int_act_mm,LEs_pot_mm,LEs_act_mm,' // &
    'Tr_act_mm,Esoil_act_mm,ET_total_act_mm,pond_evap_mm,' // &
    'Tcan_pot_C,Tcan_act_C,Tac_act_C,Tsoil_act_C,Hc_act_W_m2,Hs_act_W_m2,Tair_C,' // &
    'canopy_resid_W_m2,outer_iters,outer_converged,outer_resid_K,' // &
    'LEs_demand_sent_mm,Tr_demand_sent_mm'
  ! The hourly component-radiation ledger. The named longwave and net-radiation
  ! columns report the COMMITTED (accepted) scene; the demand-state values are
  ! appended as their own columns.
  OPEN (newunit=unit_rad,file='output/hourly_radiation.csv',status='replace',action='write')
  rad_open = .TRUE.
  WRITE (unit_rad,'(A)') 'day,gauss_idx,solar_hour,SW_down,SW_abs_can,SW_abs_soil,' // &
    'SW_ref_can,SW_ref_soil,LW_down,LW_emit_can,LW_emit_soil,LW_can_to_soil,LW_soil_to_can,' // &
    'LW_can_net,LW_soil_net,LW_up_sky,Rnc,Rns,Rn_total,alpha_soil,alpha_canopy,alpha_mix_diag,tau_sw,tau_lw,' // &
    'LW_emit_can_demand,LW_emit_soil_demand,LW_can_to_soil_demand,LW_soil_to_can_demand,' // &
    'LW_can_net_demand,LW_soil_net_demand,LW_up_sky_demand,Rnc_demand,Rns_demand,Rn_total_demand'
  ! (rec 95): ON-ONLY. Gate A compares the FILE SET, not only contents, so
  ! this file must NOT exist when the record is 0.
  IF (tleaf_class_diag >= 1) THEN
    OPEN (newunit=unit_tlc,file='output/hourly_tleaf_classes.csv',status='replace',action='write')
    tlc_open = .TRUE.
    WRITE (unit_tlc,'(A)') 'day,gauss_idx,solar_hour,LAI_m2_m2,f_sun,Tair_C,Tac_C,e_ac_Pa,' // &
      'Tcan_demand_C,Tcan_committed_C,T_sun_demand_C,T_shade_demand_C,' // &
      'T_area_mean_demand_C,T_radiometric_mean_demand_C,' // &
      'T_sun_committed_C,T_shade_committed_C,' // &
      'sun_minus_shade_K,sun_minus_canopy_K,sun_minus_air_K,canopy_minus_air_K,' // &
      'Rn_leaf_sun_W_m2,Rn_leaf_shade_W_m2,rb_leaf_s_m,rs_sun_s_m,rs_shade_s_m,' // &
      'bulk_check_minus_Tcan_K,LEc_cap_gap_W_m2,LEc_cap_factor'
  END IF
END IF
OPEN (newunit=unit_tleaf,file='output/daily_tleaf_tair.csv',status='replace',action='write')
! The Tleaf_C column is the DEMAND-PASS (Leuning-branch potential) canopy
! temperature, bit-identical to the potential canopy temperature of
! daily_eb_actual.csv and up to +4.25 K above the committed state on a day.
! Demand-pass quantities are labeled as potential diagnostics everywhere, so
! the COMMITTED canopy temperature is appended here as Tcan_act_C, the quantity
! an infrared measurement compares with; the existing columns are untouched and
! their basis is stated in INPUT_GUIDE.md.
WRITE (unit_tleaf,'(A)') 'day,Tleaf_C,Tair_C,Tleaf_minus_Tair_C,Tcan_act_C'
OPEN (newunit=unit_tsoil,file='output/daily_tsoil.csv',status='replace',action='write')
WRITE (unit_tsoil,'(A)') 'day,Tsoil_surf_C'
OPEN (newunit=unit_resist,file='output/daily_resistances.csv',status='replace',action='write')
WRITE (unit_resist,'(A)') 'day,gsw_mean_mol_m2_s,rcs_s_m,rca_s_m,rsa_s_m,rss_s_m,LEc_MJ_m2_d,LEs_MJ_m2_d,Tr_mm_d,Es_mm_d,' // &
  'gsw_day_mol_m2_s,rcs_day_s_m,rca_day_s_m,rsa_day_s_m,rss_day_s_m,Tr_day_mm,Es_day_mm'
OPEN (newunit=unit_growth,file='output/daily_growth.csv',status='replace',action='write')
IF (lai_from_organs == 1) THEN
  ! (rec 84 = 1 ONLY): APPEND SLA_out_m2_kg; the OFF header is the incumbent line below.
  WRITE (unit_growth,'(A)') 'day,LAI_m2_m2,plant_height_m,W_green_g_m2,W_dead_g_m2,' // &
    'W_stem_g_m2,W_root_g_m2,W_storage_g_m2,AGB_g_m2,GrainYield_g_m2,SLA_out_m2_kg'
ELSE
  WRITE (unit_growth,'(A)') 'day,LAI_m2_m2,plant_height_m,W_green_g_m2,W_dead_g_m2,' // &
    'W_stem_g_m2,W_root_g_m2,W_storage_g_m2,AGB_g_m2,GrainYield_g_m2'
END IF
OPEN (newunit=unit_yield,file='output/daily_yield_components.csv',status='replace',action='write')
WRITE (unit_yield,'(A)') 'day,dvs2,W_grain_g_m2,W_chaff_g_m2,W_spike_g_m2,' // &
  'kernel_number_m2,TKW_g,grain_yield_g_m2,harvest_index'
! Mineral diagnostics (Task 8): long format, one row per (day, mineral).
OPEN (newunit=unit_minerals,file='output/daily_minerals.csv',status='replace',action='write')
WRITE (unit_minerals,'(A)') 'day,mineral,uptake_massflow_mg_m2,uptake_active_mg_m2,' // &
  'remob_to_grain_mg_m2,soil_labile_mg_m2,root_mg_m2,leaf_mg_m2,stem_mg_m2,chaff_mg_m2,' // &
  'grain_mg_m2,balance_resid_mg_m2,grain_mg_kg,cap_binding'
OPEN (newunit=unit_organic,file='output/daily_organic_N.csv',status='replace',action='write')
WRITE (unit_organic,'(A)') 'day,res_fast_C_g_m2,res_fast_N_g_m2,res_slow_C_g_m2,res_slow_N_g_m2,' // &
  'soil_mic_C_g_m2,soil_mic_N_g_m2,decomp_C_g_m2,decomp_N_g_m2,respired_C_g_m2,' // &
  'net_mineralization_NH4_g_m2,immobilization_NH4_g_m2,immobilization_NO3_g_m2,' // &
  'microbial_turnover_N_g_m2,residue_N_in_g_m2,soilN_balance_resid_g_m2'
! The humus diagnostic of record 93, in its own file, so that no existing CSV
! gains a column. It is opened, written and closed ONLY when the record is at
! or above 1, the same on-only pattern as the record-87 diagnostic, so that the
! configuration with the record off does not gain a file the baseline does not
! have.
humusN_open = (soil_n_mech >= 1)
! rec 97: ON-only diagnostic, same pattern as rec 93's above,
! so the OFF arm's output FILE SET is unchanged (Gate A file-set clause).
fnc_open = (itemF_Nclose_mode >= 1)
IF (fnc_open) THEN
  OPEN (newunit=unit_fnc,file='output/daily_itemF_nclose.csv',status='replace',action='write')
  WRITE (unit_fnc,'(A)') 'day,itemF_Nclose_mode,dvs2,LAI_assim,N_leaf_gN_m2,' // &
    'sln_cap_gN_m2,Narea_mean_gN_m2,syn_clip_day_gN_m2,syn_clip_cum_gN_m2,' // &
    'N_mobile_veg_gN_m2,mob_close_day_gN_m2,mob_close_cum_gN_m2,' // &
    'grain_N_demand_postanth_day_gN_m2,gnf_sink_day_gN_m2'
END IF
! rec 121: ON-only diagnostic, same file-set clause as rec 97 above.
npg_open = (n_payment_gate >= 1)
IF (npg_open) THEN
  OPEN (newunit=unit_npg,file='output/daily_n_payment.csv',status='replace',action='write')
  WRITE (unit_npg,'(A)',ADVANCE='NO') 'day,dvs2,mode,live,c_L_gN_gDM,c_S_gN_gDM,c_R_gN_gDM,' // &
    'pool_mobile_in_gN_m2,pool_assim_in_gN_m2,lux_avail_gN_m2,lux_used_gN_m2,' // &
    'D_L_gN_m2,D_S_gN_m2,D_R_gN_m2,p_L,p_S,p_R,paid_L_gN_m2,paid_S_stem_gN_m2,paid_S_chaff_gN_m2,paid_R_gN_m2,' // &
    'pool_out_gN_m2,mobile_left_gN_m2,dw_leaf_pot_g_m2,dw_leaf_g_m2,dw_stem_pot_g_m2,dw_stem_g_m2,' // &
    'dw_stor_pot_g_m2,dw_stor_g_m2,dw_root_pot_g_m2,dw_root_g_m2,unpaid_CH2O_g_m2,wsc_deposit_gDM_m2,' // &
    'resid_CH2O_g_m2,lamina_unpaid_gDM_m2,lamina_unpaid_cum_gDM_m2,N_leaf_gN_m2,N_stem_gN_m2,N_root_gN_m2,W_reserve_g_m2,' // &
    'pool_root_in_gN_m2,root_export_gN_m2,paid_R_root_gN_m2,paid_R_shoot_gN_m2'
  ! Record 138: four on-only extra columns, appended to the SAME header line by a
  ! non-advancing write, so that the file is unchanged when the record is off
  ! (the blank advancing write below supplies the record terminator either way).
  ! Record 141: three on-only extra columns by the same technique.
  IF (lamina_n_remob_pre >= 1) WRITE (unit_npg,'(A)',ADVANCE='NO') &
    ',rem_lam_avail_gN_m2,rem_lam_used_gN_m2,rem_lam_used_cum_gN_m2'
  ! rec 143: eight ON-only extra columns, SAME non-advancing-WRITE
  ! technique, OFF (n_status_construction_price = 0) file.
  IF (n_status_construction_price >= 1) WRITE (unit_npg,'(A)',ADVANCE='NO') &
    ',nsc_nni_veg,nsc_price_L,nsc_price_S,nsc_price_chaff,nsc_price_R,nsc_cut_L_gDM,nsc_cut_S_gDM,nsc_cut_R_gDM'
  WRITE (unit_npg,'(A)') ''
END IF
lna_open = (lai_n_arrest >= 1)
IF (lna_open) THEN
  OPEN (newunit=unit_lna,file='output/daily_lai_n_arrest.csv',status='replace',action='write')
  WRITE (unit_lna,'(A)',ADVANCE='NO') 'day,dvs2,mode,live,p_L,inc_area_m2_m2,removed_m2_m2,removed_cum_m2_m2,' // &
    'removed_dm_gDM_m2,gai_before_m2_m2,gai_after_m2_m2,N_leaf_gN_m2,SLN_green_in_gN_m2,shed_need_m2_m2,' // &
    'shed_m2_m2,shed_cum_m2_m2,sl_bump_ms'
  ! Record 144: three on-only extra columns on the same row, by the same
  ! non-advancing-write technique records 141, 142 and 143 use on
  ! daily_n_payment.csv.
  IF (carbon_coupled_leaf_expansion >= 1) WRITE (unit_lna,'(A)',ADVANCE='NO') &
    ',cbc_surplus_frac,cbc_f_mean_expanding,cbc_area_bonus_m2_m2_day'
  ! Record 179: seven on-only extra columns on the same row, by the same
  ! non-advancing-write technique records 141 to 144 use.
  IF (lamina_expansion_n_limit >= 1) WRITE (unit_lna,'(A)',ADVANCE='NO') &
    ',lnl_Nph_gN_m2,lnl_slnph_gN_m2,lnl_maxNpha_gN_m2,lnl_minNpha_gN_m2,' // &
    'lnl_a_sup_m2_m2,lnl_fN,lnl_refused_m2_m2'
  WRITE (unit_lna,'(A)') ''
END IF
! Record 147: the on-only daily retention diagnostic; with the record at 0
! nothing is opened. The released nitrogen column is the nitrogen resorbed on
! the days the criterion fired and is not an attributed flux.
! Record 148: an on-only daily file; nothing is created at 0.
wfl_open = (wsc_fate_law >= 1)
IF (wfl_open) THEN
  OPEN (newunit=unit_wfl,file='output/daily_wsc_fate.csv',status='replace',action='write')
  WRITE (unit_wfl,'(A)') 'day,dap_anth,dvs2,mode,gfb_end_bound,W_reserve_pre_g_m2,' // &
    'term_a_refused_g_m2_d,term_b_released_g_m2_d,draw_deficit_g_m2_d,cres_respired_g_m2_d,' // &
    'W_reserve_post_g_m2,W_stem_g_m2,wsc_frac_of_stem'
END IF
tna_open = (tiller_n_afford >= 1)
IF (tna_open) THEN
  OPEN (newunit=unit_tna,file='output/daily_tiller_n_afford.csv',status='replace',action='write')
  WRITE (unit_tna,'(A)') 'day,taun,p_L_used,c_L_gN_gDM,days_bind,' // &
    'dens_0,dN_0,DN_0,SN_0,phiN_0,phiC_0,phi_0,dens_1,dN_1,DN_1,SN_1,phiN_1,phiC_1,phi_1,' // &
    'dens_2,dN_2,DN_2,SN_2,phiN_2,phiC_2,phi_2,dens_3,dN_3,DN_3,SN_3,phiN_3,phiC_3,phi_3,' // &
    'dens_4,dN_4,DN_4,SN_4,phiN_4,phiC_4,phi_4,dens_5,dN_5,DN_5,SN_5,phiN_5,phiC_5,phi_5'
END IF
rfl_open = (root_front_law >= 1)
IF (rfl_open) THEN
  OPEN (newunit=unit_rfl,file='output/daily_root_front.csv',status='replace',action='write')
  WRITE (unit_rfl,'(A)') 'day,dvs2,days_emerged,r_z_mm_d,z_legacy_m,z_root_m,z_max_m'
END IF
rsf_open = (repro_surplus_fate == 1)
IF (rsf_open) THEN
  OPEN (newunit=unit_rsf,file='output/daily_repro_surplus_fate.csv',status='replace',action='write')
  WRITE (unit_rsf,'(A)') 'day,dvs2,W_grain_g_m2,W_stem_g_m2,W_reserve_g_m2,surplus_cum_g_m2,wsc_labeled_cum_g_m2,' // &
    'refused_to_residual_day_g_m2,refused_to_residual_cum_g_m2,days_refused'
END IF
ncs_open = (n_construction_status >= 1)
IF (ncs_open) THEN
  OPEN (newunit=unit_ncs,file='output/daily_n_construction.csv',status='replace',action='write')
  ! Record 156: two on-only columns naming the dilution-curve coefficient pair
  ! that produced the critical concentration of that day.
  ncs_hdr_base = 'day,dvs2,mode,live,W_t_ha,c_crit_gN_gDM,f_crit,' // &
    'c_floor_L,c_target_L,c_target_eff_L,c_actual_L,paid_pass1_L_gN_m2,paid_pass2_L_gN_m2,' // &
    'c_floor_S,c_target_S,c_target_eff_S,c_actual_S,paid_pass1_S_gN_m2,paid_pass2_S_gN_m2,' // &
    'c_floor_R,c_target_R,c_target_eff_R,c_actual_R,paid_pass1_R_gN_m2,paid_pass2_R_gN_m2,' // &
    'p_L,p_S,p_R,lamina_below_target_cum_gDM_m2,days_fcrit,days_pL_lt1'
  IF (n_construction_ceiling_curve >= 1) THEN
    ncs_hdr_base = TRIM(ncs_hdr_base) // ',fcrit_curve_a,fcrit_curve_b'
  END IF
  ! Record 157: three on-only columns - the lamina storage headroom offered today,
  ! the nitrogen actually stored, and the capacity-binding flag.
  IF (lamina_n_storage >= 1) THEN
    ncs_hdr_base = TRIM(ncs_hdr_base) // ',lamina_store_cap_gN_m2,lamina_stored_gN_m2,lamina_cap_bind'
  END IF
  ! rec 157 mode 3: the four [N]mob numerator terms, each as its own column, so
  ! the reader can see exactly what each contributes and what dropping one would do.
  IF (lamina_n_storage == 3) THEN
    ncs_hdr_base = TRIM(ncs_hdr_base) // ',nmob_num_gN_m2,nmob_t_mobile_gN_m2,' // &
      'nmob_t_NO3shoot_gN_m2,nmob_t_stemlux_gN_m2,nmob_t_chafflux_gN_m2'
  END IF
  WRITE (unit_ncs,'(A)') TRIM(ncs_hdr_base)
END IF
cfo_open = (sink_full_carbon_fate == 1)
IF (cfo_open) THEN
  OPEN (newunit=unit_cfo,file='output/daily_carbon_overflow.csv',status='replace',action='write')
  WRITE (unit_cfo,'(A)') 'day,dvs2,W_grain_g_m2,W_stem_g_m2,W_reserve_g_m2,surplus_day_g_m2,wsc_stored_day_g_m2,' // &
    'loss_day_g_m2,loss_cum_g_m2,wsc_stored_cum_g_m2,sink_full_flag'
END IF
gnm_open = (grain_n_min_conc >= 1)
IF (gnm_open) THEN
  OPEN (newunit=unit_gnm,file='output/daily_grain_n_floor.csv',status='replace',action='write')
  WRITE (unit_gnm,'(A)',ADVANCE='NO') 'day,dvs2,W_grain_g_m2,N_grain_gN_m2,grain_N_conc_pct,' // &
    'room_dm_118_g_m2_d,room_dm_125_g_m2_d,room_dm_floor_g_m2_d,' // &
    'to_grain_before_g_m2_d,to_grain_after_g_m2_d,refused_dm_day_g_m2,refused_dm_cum_g_m2,floor_binding_flag'
  ! Record 130 mode 3: seven on-only extra columns, appended by the same
  ! non-advancing write records 138 and 141 use on daily_n_payment.csv.
  IF (grain_n_min_conc == 3) WRITE (unit_gnm,'(A)',ADVANCE='NO') &
    ',gnm_Tt_degCd,gnm_Nstru_pergrain_ugN,gnm_c_stru_pct,gnm_room_floor_gDM,' // &
    'gnm_Nstru_paid_gN_m2,gnm_Nsto_gN_m2,gnm_phase_flag'
  WRITE (unit_gnm,'(A)') ''
END IF
! rec 139: ON-only, NEW file (day, layer, mf/diff/act/cap/realized, g N m-2 d-1).
nud_open = (n_uptake_diffusion >= 1)
IF (nud_open) THEN
  OPEN (newunit=unit_nud,file='output/daily_n_uptake_terms.csv',status='replace',action='write')
  ! Record 155: six WHOLE-PLANT columns appended on-only. They are whole-plant
  ! quantities written on every layer row of the day, because the row key of this
  ! file is day and layer and the record adds no per-layer quantity, so the value
  ! is repeated down the day ten rows by design.
  IF (n_uptake_regulation >= 1) THEN
    WRITE (unit_nud,'(A)') 'day,layer,mf_gN_m2,diff_gN_m2,act_gN_m2,cap_gN_m2,realized_gN_m2,' // &
      'E_N,Nmob_g_g,cap_Nmax_gN_m2,cap_Nc_gN_m2,demand_before_gN_m2,demand_after_gN_m2'
  ELSE
    WRITE (unit_nud,'(A)') 'day,layer,mf_gN_m2,diff_gN_m2,act_gN_m2,cap_gN_m2,realized_gN_m2'
  END IF
END IF
! rec 140: ON-only, NEW file (day, mineral, mf/diff/hats/cap/realized, mg element m-2 d-1).
msm_open = (mineral_supply_mech >= 1)
IF (msm_open) THEN
  OPEN (newunit=unit_msm,file='output/daily_mineral_uptake_terms.csv',status='replace',action='write')
  WRITE (unit_msm,'(A)') 'day,mineral,mf_mg_m2,diff_mg_m2,hats_mg_m2,cap_mg_m2,realized_mg_m2'
END IF
kpt_open = (kernel_potential_law >= 1)
IF (kpt_open) THEN
  OPEN (newunit=unit_kpt,file='output/daily_kernel_potential.csv',status='replace',action='write')
  WRITE (unit_kpt,'(A)') 'day,dvs2,kernel_number_m2,W_grain_g_m2,N_grain_gN_m2,mgN_per_kernel,kpt_tt_degCd,kpt_sup_cum_g_m2,' // &
    'kpt_set,W_pot_mg,gnc_cap_gN_m2,gnc_room_gN_m2,gnc_cut_day_gN_m2,gnc_cut_cum_gN_m2,gnc_days_bound,gfb_end_bound'
END IF
IF (humusN_open) THEN
  OPEN (newunit=unit_humusN,file='output/daily_humus_N.csv',status='replace',action='write')
  WRITE (unit_humusN,'(A)') 'day,soil_n_mech,NHUM_gN_m2,VMINH_day_gN_m2,QMINH_gN_m2,' // &
    'FTH,FH_mean,K2POT_d-1,T_soil_day_c,NHUM_drawdown_frac'
END IF
OPEN (newunit=unit_pheno,file='output/daily_phenology.csv',status='replace',action='write')
WRITE (unit_pheno,'(A)') 'day,year,month,dom,doy,dap,stage_name,stage_idx,cum_tdu,vf,df,photoperiod_h,t_dev_pre_c,t_dev_post_c'
! Record 64: the daily rank-chronology diagnostic, opened only when the record
! is on, so that with it off no file is created and the output file set is
! exactly what it was.
IF (canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_clr,file='output/daily_canopy_ranks.csv',status='replace',action='write')
  WRITE (unit_clr,'(A)') 'day,dap,stage_idx,tau_Cd,HS,HS_evers,ranks_appeared,ranks_ligulated,' // &
    'ranks_appeared_evers,ranks_ligulated_evers,phyllochron_local_Cd,FLN_emergent,FLN_tip,' // &
    'FLN_evers,FLN_evers_tip,tau_app_flag_Cd,tau_lig_flag_Cd'
  clr_diag_open = .TRUE.
END IF
! Record 65: the daily blade-area diagnostic, opened only when the record is on.
IF (canopy_blade_area == 1 .AND. canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_cba,file='output/daily_canopy_blade_area.csv',status='replace',action='write')
  WRITE (unit_cba,'(A)') 'day,dap,stage_idx,tau_Cd,HS,ranks_appeared,ranks_ligulated,ranks_expanding,' // &
    'A_shoot_pot_cm2,A_shoot_cm2,A_shoot_pot_m2_per_shoot,A_shoot_m2_per_shoot,FLN_emergent'
  cba_diag_open = .TRUE.
END IF
! Record 66: the daily axis-population diagnostic, same on-only pattern. The
! blade-area indices are ZERO unless record 65 is also on, because without it
! there is no per-axis area to weight.
IF (canopy_tillering == 1 .AND. canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_ct,file='output/daily_canopy_tillering.csv',status='replace',action='write')
  WRITE (unit_ct,'(A)') 'day,dap,stage_idx,tau_Cd,HS,fPAR_model,axes_live_m2,axes_born_m2,' // &
    'axes_live_per_plant,dens_MS,dens_T1,dens_T2,dens_T3,dens_T4,dens_T5,' // &
    'blade_area_index_pot,blade_area_index_real,cessation_reason,site_filling_implied'
  ct_diag_open = .TRUE.
END IF
! Record 67: the daily senescence diagnostic, same on-only pattern. The green
! area index and the standing-dead index are GROUND-AREA quantities (m2 m-2),
! and the standing-dead total adds the blade area of axes culled by the
! mortality phase, which stays standing.
IF (canopy_leaf_senesc == 1 .AND. canopy_tillering == 1 .AND. &
    canopy_blade_area == 1 .AND. canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_cls,file='output/daily_canopy_senescence.csv',status='replace',action='write')
  WRITE (unit_cls,'(A)') 'day,dap,stage_idx,tau_Cd,HS_MS,Gl_MS,Sl_MS,phase_MS,' // &
    'Gl_T1,Gl_T2,Gl_T3,Gl_T4,Gl_T5,green_area_index,standing_dead_index,' // &
    'standing_dead_pool,standing_dead_culled,standing_dead_total,total_blade_index'
  cls_diag_open = .TRUE.
END IF
! Record 68: the daily non-lamina diagnostic, same on-only pattern. TWO BASES
! are written side by side and are never mixed: the green destructive basis is
! GREEN organs only on projected (silhouette) areas, the quantity measured by
! dissection and scanning at flowering; the gap-fraction basis adds the
! standing dead blade and sheath pools and converts the cylindrical organs to
! the half-total surface a gap-fraction inversion returns. The n+1 and n+3
! sheath-persistence arms and the 10-20 cm2 ear-area bracket are emitted as
! their own columns, so the sensitivities are visible in every run.
IF (canopy_sheath_ear == 1 .AND. canopy_leaf_senesc == 1 .AND. &
    canopy_tillering == 1 .AND. canopy_blade_area == 1 .AND. &
    canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_cse,file='output/daily_canopy_sheath_ear.csv',status='replace',action='write')
  WRITE (unit_cse,'(A)') 'day,dap,stage_idx,tau_Cd,HS_MS,Sl_MS,' // &
    'sheath_green_index,sheath_dead_instant,sheath_dead_pool,' // &
    'sheath_green_n1_arm,sheath_green_n3_arm,sheath_green_nolag_arm,' // &
    'stem_exposed_index,ear_ped_index,ear_ped_index_lo,ear_ped_index_hi,' // &
    'ear_only_index_measured,' // &
    'green_blade_index,dead_blade_total,nonlamina_green,nonlamina_over_lamina,' // &
    'PAIg_green_Abichou_basis,PAI_LI2200C_basis,PAIg_cm2_per_live_axis'
  cse_diag_open = .TRUE.
END IF
! Record 69: the closure and limiter diagnostic of the switch. Unlike records
! 64 to 68 this step MOVES STATE, so these columns are how it is audited: the
! lamina area demand and the carbon actually debited for it, the hard
! non-negative backlog, the three demographic-limiter counters, the realized
! specific leaf area against the measured ranges, and the two leaf area indices
! kept explicitly apart.
IF (canopy_sink_lai == 1 .AND. canopy_sheath_ear == 1 .AND. &
    canopy_leaf_senesc == 1 .AND. canopy_tillering == 1 .AND. &
    canopy_blade_area == 1 .AND. canopy_leaf_ranks == 1) THEN
  OPEN (newunit=unit_csl,file='output/daily_canopy_sink.csv',status='replace',action='write')
  WRITE (unit_csl,'(A)') 'day,dap,stage_idx,LAI_sink_green_lamina,PAI_diag_LI2200C_basis,' // &
    'lai_basis_Ncapacity_only,area_built_cum,dA_day,lamina_demand_gDM,lamina_alloc_gDM,' // &
    'backlog_gDM,limit_active,fln_cut,cohorts_aborted,births_suspended,' // &
    'realized_SLA_m2_per_kg,W_green_gDM,partition_sum_err'
  csl_diag_open = .TRUE.
END IF
! Records 82, 83 and 84: one on-only daily diagnostic carrying the supply-ratio
! gate (per-cohort window supply, demand, ratio and births), the root-zone soil
! relative water content and the expansion factor with its two printed
! brackets, the arrested against the clock blade index, the axis envelope flags
! and the ear envelope (reported, never applied), the organ leaf area index and
! the output specific leaf area. Opened only when one of the three is on.
IF (tiller_carbon_gate == 1 .OR. leaf_expansion_water == 1 .OR. lai_from_organs == 1 .OR. &
    tiller_survival_mech >= 1) THEN
  OPEN (newunit=unit_d2,file='output/daily_d2_organs.csv',status='replace',action='write')
  WRITE (unit_d2,'(A)') 'day,dap,stage_idx,tau_Cd,HS,C_sur_prevday_gDM,' // &
    'sup_T1,sup_T2,sup_T3,sup_T4,sup_T5,dem_T1,dem_T2,dem_T3,dem_T4,dem_T5,' // &
    'ratio_T1,ratio_T2,ratio_T3,ratio_T4,ratio_T5,dens_MS,dens_T1,dens_T2,dens_T3,dens_T4,dens_T5,' // &
    'axes_live_m2,axes_born_m2,env600_days,env800_days,Ne_env_lo_m2,Ne_env_hi_m2,' // &
    'SRWC_rw,Se_rw_alt,theta_rw,theta_fc_rw,f_w,f_w_200line,f_w_800line,' // &
    'blade_index_arrested,blade_index_clock,LAI_organs_green,LAI_used,SLA_out_m2_kg,W_green_gDM,' // &
    'lamina_demand_gDM,lamina_alloc_gDM,backlog_gDM,' // &
    'gai_died_day,f_died_day,sen_dm_gDM,sen_N_out_gN,sen_N_dead_gN,sen_N_resorb_gN,N_leaf_gN,Narea_state_gNm2'
  d2_diag_open = .TRUE.
END IF
! Record 99 mode 2: the daily affordability row, one per stepped window day,
! opened only when mode 2 is selected.
! Record 99 mode 4: the mode-4 daily row - the mode-2 columns plus the
! mortality-shape fraction, the running estimate, the target, the target deaths
! and the endpoint flag - opened only when mode 4 is selected.
IF (tiller_survival_mech == 4) THEN
  OPEN (newunit=unit_tca4,file='output/daily_tiller_afford4.csv',status='replace',action='write')
  WRITE (unit_tca4,'(A)') 'day,dap,tau_Cd,taun,G_prev_gCH2O,L_tot,F_tot,sum_s,ident_gCH2O,dead_cum_m2,days_skip,' // &
    'F_taun,U_m2,target_m2,dead_target_m2,endpoint_flag,' // &
    'dens_0,hs_0,L_0,s_0,S_ax_0,D_ax_0,cumS_0,cumD_0,phi_0,' // &
    'dens_1,hs_1,L_1,s_1,S_ax_1,D_ax_1,cumS_1,cumD_1,phi_1,' // &
    'dens_2,hs_2,L_2,s_2,S_ax_2,D_ax_2,cumS_2,cumD_2,phi_2,' // &
    'dens_3,hs_3,L_3,s_3,S_ax_3,D_ax_3,cumS_3,cumD_3,phi_3,' // &
    'dens_4,hs_4,L_4,s_4,S_ax_4,D_ax_4,cumS_4,cumD_4,phi_4,' // &
    'dens_5,hs_5,L_5,s_5,S_ax_5,D_ax_5,cumS_5,cumD_5,phi_5'
  tca4_diag_open = .TRUE.
END IF
OPEN (newunit=unit_nitrogen,file='output/daily_nitrogen.csv',status='replace',action='write')
IF (soilN_layered == 1) THEN
  ! Layered mode: same header + the two new ion-resolved uptake columns
  ! appended at row END (incumbent header when the flag is 0).
  WRITE (unit_nitrogen,'(A)',ADVANCE='NO') 'day,N_soil_min_gN_m2,N_fert_day_gN_m2,N_mineral_day_gN_m2,' // &
  'N_leach_day_gN_m2,N_soil_uptake_day_gN_m2,N_assim_biochem_gross_day_gN_m2,' // &
  'N_assim_biochem_retained_day_gN_m2,N_demand_day_gN_m2,' // &
  'f_Tr_N,f_CO2_N,N_remob_leaf_day_gN_m2,N_remob_stem_day_gN_m2,' // &
  'N_alloc_leaf_day_gN_m2,N_alloc_stem_day_gN_m2,N_alloc_root_day_gN_m2,N_alloc_storage_day_gN_m2,' // &
  'N_leaf_gN_m2,N_dead_leaf_gN_m2,N_stem_gN_m2,N_root_gN_m2,N_storage_gN_m2,N_balance_resid_gN_m2,' // &
  'leafN_conc_gN_gDM,stemN_conc_gN_gDM,rootN_conc_gN_gDM,storageN_conc_gN_gDM,grain_protein_pct,' // &
  'N_NO3_root_gN_m2,N_NO3_shoot_gN_m2,N_assim_root_day_gN_m2,N_assim_shoot_day_gN_m2,f_photoresp_diag,' // &
  'N_uptake_massflow_day_gN_m2,N_uptake_active_day_gN_m2,' // &
  'grain_N_import_day_gN_m2,grain_N_from_assim_day_gN_m2,N_surplus_veg_day_gN_m2,' // &
  'N_NO3_soil_gN_m2,N_NH4_soil_gN_m2,N_urea_soil_gN_m2,N_hydrolysis_day_gN_m2,' // &
  'N_nitrification_day_gN_m2,soilN_balance_resid_gN_m2,' // &
  'N_volat_day_gN_m2,N_denit_day_gN_m2,TAN_surf_gN_m2,' // &
  'veg_N_demand_day_gN_m2,grain_N_demand_postanth_day_gN_m2,' // &
  'N_uptake_pre_cum_gN_m2,N_uptake_post_cum_gN_m2,NHI_diag,' // &
  'N_mobile_veg_gN_m2,N_from_mobile_day_gN_m2,wholeleafN_conc_gN_gDM,' // &
  'N_grain_gN_m2,N_chaff_gN_m2,grainN_conc_gN_gDM,chaffN_conc_gN_gDM,' // &
  'S_active_pot_gN_m2,S_massflow_pot_gN_m2,U_pot_gN_m2,f_Nstatus,a_root,' // &
  'U_act_NO3_gN_m2,U_act_NH4_gN_m2'
  ! Record 158: five on-only columns appended at the end of the row by a
  ! non-advancing write on the data rows, the same technique records 138 and 141
  ! to 143 use on daily_n_payment.csv.
  IF (nitrate_drain_rate_basis >= 1) WRITE (unit_nitrogen,'(A)',ADVANCE='NO') &
    ',ndb_nmax_store_umol_m2leaf_s,ndb_nmax_rate_umol_m2leaf_s,' // &
    'ndb_rate_prev_gN_m2_d,ndb_alpha_g,ndb_alpha_s'
  ! rec 159: five ON-only columns appended at row END, same technique.
  WRITE (unit_nitrogen,'(A)') ''
  ! Per-layer soil mineral-N state + bottom leach flux (layered mode only).
  BLOCK
    CHARACTER(LEN=2048) :: hdr
    CHARACTER(LEN=16)   :: tag
    INTEGER :: il
    OPEN (newunit=unit_soilN_lay,file='output/daily_soilN_layers.csv',status='replace',action='write')
    hdr = 'day'
    DO il = 1, NSoil_layers
      WRITE(tag,'(A,I0)') ',NO3_l', il
      hdr = TRIM(hdr) // TRIM(tag)
    END DO
    DO il = 1, NSoil_layers
      WRITE(tag,'(A,I0)') ',NH4_l', il
      hdr = TRIM(hdr) // TRIM(tag)
    END DO
    hdr = TRIM(hdr) // ',leach_bot_gN_m2_d,q_bot_mm_d,adv_substeps'
    WRITE (unit_soilN_lay,'(A)') TRIM(hdr)
  END BLOCK
ELSE
WRITE (unit_nitrogen,'(A)',ADVANCE='NO') 'day,N_soil_min_gN_m2,N_fert_day_gN_m2,N_mineral_day_gN_m2,' // &
  'N_leach_day_gN_m2,N_soil_uptake_day_gN_m2,N_assim_biochem_gross_day_gN_m2,' // &
  'N_assim_biochem_retained_day_gN_m2,N_demand_day_gN_m2,' // &
  'f_Tr_N,f_CO2_N,N_remob_leaf_day_gN_m2,N_remob_stem_day_gN_m2,' // &
  'N_alloc_leaf_day_gN_m2,N_alloc_stem_day_gN_m2,N_alloc_root_day_gN_m2,N_alloc_storage_day_gN_m2,' // &
  'N_leaf_gN_m2,N_dead_leaf_gN_m2,N_stem_gN_m2,N_root_gN_m2,N_storage_gN_m2,N_balance_resid_gN_m2,' // &
  'leafN_conc_gN_gDM,stemN_conc_gN_gDM,rootN_conc_gN_gDM,storageN_conc_gN_gDM,grain_protein_pct,' // &
  'N_NO3_root_gN_m2,N_NO3_shoot_gN_m2,N_assim_root_day_gN_m2,N_assim_shoot_day_gN_m2,f_photoresp_diag,' // &
  'N_uptake_massflow_day_gN_m2,N_uptake_active_day_gN_m2,' // &
  'grain_N_import_day_gN_m2,grain_N_from_assim_day_gN_m2,N_surplus_veg_day_gN_m2,' // &
  'N_NO3_soil_gN_m2,N_NH4_soil_gN_m2,N_urea_soil_gN_m2,N_hydrolysis_day_gN_m2,' // &
  'N_nitrification_day_gN_m2,soilN_balance_resid_gN_m2,' // &
  'N_volat_day_gN_m2,N_denit_day_gN_m2,TAN_surf_gN_m2,' // &
  'veg_N_demand_day_gN_m2,grain_N_demand_postanth_day_gN_m2,' // &
  'N_uptake_pre_cum_gN_m2,N_uptake_post_cum_gN_m2,NHI_diag,' // &
  'N_mobile_veg_gN_m2,N_from_mobile_day_gN_m2,wholeleafN_conc_gN_gDM,' // &
  'N_grain_gN_m2,N_chaff_gN_m2,grainN_conc_gN_gDM,chaffN_conc_gN_gDM,' // &
  'S_active_pot_gN_m2,S_massflow_pot_gN_m2,U_pot_gN_m2,f_Nstatus,a_root'
  ! rec 158: five ON-only columns appended at row END, same technique.
  IF (nitrate_drain_rate_basis >= 1) WRITE (unit_nitrogen,'(A)',ADVANCE='NO') &
    ',ndb_nmax_store_umol_m2leaf_s,ndb_nmax_rate_umol_m2leaf_s,' // &
    'ndb_rate_prev_gN_m2_d,ndb_alpha_g,ndb_alpha_s'
  ! rec 159: five ON-only columns appended at row END, same technique.
  WRITE (unit_nitrogen,'(A)') ''
END IF
OPEN (newunit=unit_carbon_growth,file='output/daily_carbon_growth.csv',status='replace',action='write')
WRITE (unit_carbon_growth,'(A)') 'day,dvs2,stage_name,An_gCO2_m2_d,A_gCO2_m2_d,' // &
  'net_assim_ch2o_g_m2_d,Rm_gCH2O_m2_d,Rg_gCH2O_m2_d,dm_gain_g_m2_d,' // &
  'f_leaf,f_stem,f_root,f_storage,' // &
  'W_green_g_m2,W_dead_g_m2,W_stem_g_m2,W_root_g_m2,W_storage_g_m2,LAI_m2_m2,plant_height_m'
OPEN (newunit=unit_rad_aero,file='output/daily_radiation_aero.csv',status='replace',action='write')
WRITE (unit_rad_aero,'(A)') 'day,alpha_mix,LAI_m2_m2,plant_height_m,Rn_MJ_m2_d,Ac_MJ_m2_d,As_MJ_m2_d,' // &
  'lambda_f,d_m,z0_m,rca_daylight_s_m,rca_all_s_m,rsa_daylight_s_m'
! The Tleaf_C column below is the DEMAND-PASS (Leuning-branch potential) canopy
! temperature and not the committed state. The committed state is appended as
! the last column Tcan_act_C; the existing columns are untouched and
! INPUT_GUIDE.md states the basis of each.
OPEN (newunit=unit_gas_exch,file='output/daily_gas_exchange.csv',status='replace',action='write')
IF (ft_pot_unstressed == 1) THEN
  ! Record 78 only: append the daylight-masked mean intercellular CO2
  ! concentration, on the same daylight mask as the daily conductance and
  ! transpiration. The 24-hour column is kept for compatibility and is not an
  ! evaluation quantity, since the night hours, where the flux identity puts the
  ! intercellular concentration above ambient, are averaged in.
  WRITE (unit_gas_exch,'(A)') 'day,An_gCO2_m2_d,A_gCO2_m2_d,N_assim_gN_m2_d,' // &
    'An_ground_umol_m2_s_daylight,An_leaf_umol_m2_s_daylight,' // &
    'gsw_day_mol_m2_s,Tr_daylight_mm_d,Tr_daylight_mmol_m2_s,' // &
    'Ci_umol_mol,VPD_Pa,Tleaf_C,Tair_C,LAI_m2_m2,plant_height_m,Ci_daylight_umol_mol,Tcan_act_C'
ELSE
  WRITE (unit_gas_exch,'(A)') 'day,An_gCO2_m2_d,A_gCO2_m2_d,N_assim_gN_m2_d,' // &
    'An_ground_umol_m2_s_daylight,An_leaf_umol_m2_s_daylight,' // &
    'gsw_day_mol_m2_s,Tr_daylight_mm_d,Tr_daylight_mmol_m2_s,' // &
    'Ci_umol_mol,VPD_Pa,Tleaf_C,Tair_C,LAI_m2_m2,plant_height_m,Tcan_act_C'
END IF
OPEN (newunit=unit_hydrodiag,file='output/daily_hydro_diag.csv',status='replace',action='write')
WRITE (unit_hydrodiag,'(A)') 'day,hydro_failed,dt_halvings,min_dt_reached,explicit_infiltration,flux_top_m_d,' // &
  'ponding_mm,top_bc_type,' // &
  'flux_3to4_m_d,K_3to4_m_d,hydro_top_m,hydro_top_infil_m,hydro_root_m,hydro_bot_m,hydro_mb_m,' // &
  'cap_active_frac,flux_top_req_m_d,flux_top_cap_m_d,flux_top_applied_m_d,k1_m_d,h1_m,theta1,h_top_m,' // &
  'h_crit_atm_m,qmax_evap_m_d,fail_events'
! wet-front Richards trace (diagnostics only; rows emitted by hydro.f90
! for days 97, 113-121, 131-146). Stable header includes Task-2 limiter columns.
IF (energy_balance_mode >= 1) THEN
  OPEN (newunit=wetfront_trace_unit,file='output/hydro_wetfront_trace.csv', &
        status='replace',action='write')
  WRITE (wetfront_trace_unit,'(A)') &
    'day,substep,t_day,dt_d,top_bc_type,ponding_mm,pr_eff_mm_d,evap_pot_mm_d,' // &
    'evap_act_mm_d,flux_top_req_m_d,flux_top_mass_m_d,top_infil_mm,theta1,Se1,' // &
    'h1_m,K1_m_d,h_top_m,h_crit_atm_m,wetfront_risk,dt_limit_d,limiter_applied,' // &
    'newton_iters,picard_iters,dt_halvings,min_dt_reached,mb_err_mm'
  wetfront_trace_enabled = .TRUE.
END IF
OPEN (newunit=unit_wbdiag,file='output/wb_components.csv',status='replace',action='write')
! Record 149 mode 2 appends an on-only 26th field, the per-day sum of the
! solver mass-conservation correction. Modes 0 and 1 write the original
! 25-field header.
  WRITE (unit_wbdiag,'(A)') 'day,pr_eff_mm_d,irr_eff_mm_d,runoff_mm_d,ponding_mm,delta_ponding_mm,' // &
    'ponding_start_mm,ponding_min_mm,ponding_max_mm,ponding_mean_mm,' // &
    'soil_storage_mm,delta_storage_mm,evapo_act_mm_d,trans_act_mm_d,ET_act_mm_d,' // &
    'hydro_top_flux_mm_d,hydro_top_flux_mass_mm_d,hydro_top_flux_mismatch_mm_d,' // &
    'hydro_root_sink_mm_d,hydro_bottom_flux_mm_d,' // &
    'soil_mb_resid_mm,balance_soil_mm,balance_surface_mm,balance_mm,wb_resid_mm_d'
OPEN (newunit=unit_hydroledger,file='output/daily_hydro_ledger.csv',status='replace',action='write')
WRITE (unit_hydroledger,'(A)') 'day,top_flux_req_mm_d,top_flux_act_mm_d,infiltration_mm_d,soil_evap_mm_d,pond_evap_mm_d,' // &
  'runoff_mm_d,root_uptake_pot_mm_d,root_uptake_act_mm_d,bottom_flux_mm_d,storage_mm,' // &
  'cum_top_flux_req_mm,cum_top_flux_act_mm,cum_infiltration_mm,cum_soil_evap_mm,cum_pond_evap_mm,' // &
  'cum_runoff_mm,cum_root_uptake_pot_mm,cum_root_uptake_act_mm,cum_bottom_flux_mm'
OPEN (newunit=unit_soilprof,file='output/daily_soil_profile_diag.csv',status='replace',action='write')
WRITE (unit_soilprof,'(A)') 'day,layer,theta_vv,theta_s,theta_r,psi_m,K_m_d,Ksat_m_d,Se,flux_in_m_d,flux_out_m_d'

! Initialize the prognostic one-dimensional soil-heat profile ONCE for the
! record-39 mode-2 path, before the day loop, where the water contents, the
! layer thicknesses and the pre-sowing soil temperature are all set. The
! preflight guard makes a failed or undersized allocation fail loudly rather
! than silently become a zero ground heat flux in the energy-balance solve.
IF (energy_balance_mode == 2) THEN
  ! An optional sensor-anchored initial soil-temperature profile. The pre-sowing
  ! soil temperature keeps only its decomposition role in the pre-season loop;
  ! the thermal initial condition comes from input/tsoil_init.csv when the deck
  ! provides it.
  BLOCK
    INTEGER :: u_tic, ios_tic, n_usable, ln_tic
    REAL(wp) :: zv, tv
    CHARACTER(len=256) :: line_tic
    OPEN (newunit=u_tic, file='input/tsoil_init.csv', status='old', action='read', iostat=ios_tic)
    IF (ios_tic == 0) THEN
      ! The file is OPTIONAL but, once present, it is REQUIRED TO BE USABLE. A
      ! malformed file stops the run with a message naming the line, instead of
      ! being partly skipped or silently abandoned for the synthetic profile:
      ! the whole point of the file is that the initial condition is measured,
      ! so quietly substituting a computed one would mislabel the run. The
      ! absent-file fallback (driver record 50) is unchanged.
      tsoil_init_n = 0; n_usable = 0; ln_tic = 0
      DO
        READ (u_tic,'(A)',iostat=ios_tic) line_tic
        IF (IS_IOSTAT_END(ios_tic)) EXIT
        IF (ios_tic /= 0) THEN
          WRITE (*,'(A,I0,A)') ' ERROR: input/tsoil_init.csv is present but could not be read' // &
            ' past line ', ln_tic, '.'
          WRITE (*,'(A)')      '        A truncated or corrupt file is not read as a short one.'
          STOP 1
        END IF
        ln_tic = ln_tic + 1
        line_tic = ADJUSTL(line_tic)
        IF (LEN_TRIM(line_tic) == 0 .OR. line_tic(1:1) == '!') CYCLE
        READ (line_tic,*,iostat=ios_tic) zv, tv
        IF (ios_tic /= 0) THEN
          WRITE (*,'(A,I0,A)') ' ERROR: input/tsoil_init.csv line ', ln_tic, &
            ' is not a "z_m,T_C" pair of numbers.'
          WRITE (*,'(A,A)')    '        line: ', TRIM(line_tic)
          WRITE (*,'(A)')      '        Expected two reals; comment lines must start with "!".'
          STOP 1
        END IF
        n_usable = n_usable + 1
        IF (n_usable > SIZE(tsoil_init_z)) THEN
          WRITE (*,'(A,I0,A,I0,A)') ' ERROR: input/tsoil_init.csv has more than ', &
            SIZE(tsoil_init_z), ' anchors (line ', ln_tic, ' is one too many).'
          WRITE (*,'(A)')      '        Dropping the deeper anchors would change the initial' // &
            ' condition silently.'
          STOP 1
        END IF
        IF (zv < 0.0_wp .OR. zv > 50.0_wp) THEN
          WRITE (*,'(A,I0,A,ES12.4)') ' ERROR: input/tsoil_init.csv line ', ln_tic, &
            ': depth outside 0 to 50 m: ', zv
          WRITE (*,'(A)')      '        Depth is positive downward in METRES (not cm), and no' // &
            ' fill value is accepted.'
          STOP 1
        END IF
        IF (tv < -60.0_wp .OR. tv > 80.0_wp) THEN
          WRITE (*,'(A,I0,A,ES12.4)') ' ERROR: input/tsoil_init.csv line ', ln_tic, &
            ': temperature outside -60 to 80 C: ', tv
          WRITE (*,'(A)')      '        Soil temperature is in C, and no fill value (-999) is' // &
            ' accepted.'
          STOP 1
        END IF
        IF (n_usable >= 2) THEN
          IF (zv <= tsoil_init_z(n_usable-1)) THEN
            WRITE (*,'(A,I0,A)') ' ERROR: input/tsoil_init.csv line ', ln_tic, &
              ': depths must increase strictly down the file.'
            WRITE (*,'(A,ES12.4,A,ES12.4,A)') '        this depth ', zv, ' m does not exceed the' // &
              ' previous ', tsoil_init_z(n_usable-1), ' m.'
            WRITE (*,'(A)')    '        The interpolation reads the anchors in file order; an' // &
              ' unsorted file would be interpolated wrongly.'
            STOP 1
          END IF
        END IF
        tsoil_init_n = n_usable
        tsoil_init_z(tsoil_init_n) = zv
        tsoil_init_T(tsoil_init_n) = tv
      END DO
      CLOSE (u_tic)
      IF (tsoil_init_n < 2) THEN
        WRITE (*,'(A,I0,A)') ' ERROR: input/tsoil_init.csv holds ', tsoil_init_n, &
          ' usable anchor(s); interpolation needs at least 2.'
        WRITE (*,'(A)')      '        Remove the file to use the driver record-50' // &
          ' initialization instead.'
        STOP 1
      END IF
      tsoil_init_from_file = .TRUE.
    END IF
  END BLOCK
  CALL initialize_soil_heat(NSoil_layers, presow_Tsoil)
  IF (.NOT. (ALLOCATED(Tsoil_prof) .AND. ALLOCATED(Cv_soil) .AND. ALLOCATED(lambda_soil) .AND. &
             ALLOCATED(tk))) THEN
    WRITE (*,*) 'ERROR: mode 2 soil-heat profile failed to allocate during initialization.'
    STOP 1
  END IF
  IF (SIZE(Tsoil_prof) < NSoil_layers .OR. SIZE(Cv_soil) < NSoil_layers .OR. &
      SIZE(lambda_soil) < NSoil_layers .OR. SIZE(tk) < NSoil_layers) THEN
    WRITE (*,*) 'ERROR: mode 2 soil-heat arrays smaller than NSoil_layers after initialization.'
    STOP 1
  END IF
  WRITE (*,'(A,I0,A,F6.2,A)') ' Stage 6 mode 2: 1-D soil-heat profile initialized (', &
    NSoil_layers, ' layers, surface ', presow_Tsoil, ' C).'
  ! Open the mode-2 soil-heat diagnostics, written per sub-step and day from the
  ! actual-flux repartition. The node temperatures are POST-ADVANCE values for
  ! that sub-step; the bottom flux is zero (zero-flux boundary); the
  ! one-dimensional ground heat flux is the surface flux of the production
  ! solver; and the transitional placeholder flux is a comparator only, never fed
  ! to the solver.
  OPEN (newunit=unit_soilheat_h,file='output/hourly_soil_heat.csv',status='replace',action='write')
  soilheat_h_open = .TRUE.
  WRITE (unit_soilheat_h,'(A)') 'day,hour,G_surface_W_m2,G_bottom_W_m2,' // &
    'heat_storage_change_J_m2,heat_residual_J_m2,Tsoil_surf_C,Tsoil_node1_C,Tsoil_node2_C,' // &
    'Tsoil_node3_C,Tsoil_nodeN_C,Rns_W_m2,Hs_W_m2,LEs_W_m2,G_1D_heat_W_m2,G_placeholder_W_m2'
  OPEN (newunit=unit_soilheat_d,file='output/daily_soil_heat.csv',status='replace',action='write')
  soilheat_d_open = .TRUE.
  WRITE (unit_soilheat_d,'(A)') 'day,G_surface_integral_J_m2,storage_change_J_m2,' // &
    'heat_residual_J_m2,Tsoil_surf_mean_C,Tsoil_surf_min_C,Tsoil_surf_max_C,Tsoil_top_mean_C'
  ! Per-layer 1-D soil-temperature profile (end-of-day Tsoil_prof), for exact-depth Tsoil (e.g. 20/55 cm).
  OPEN (newunit=unit_tsoilprof,file='output/daily_tsoil_profile.csv',status='replace',action='write')
  tsoilprof_open = .TRUE.
  WRITE (swc_line,'(A)') 'day'
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",Tsoil_L",I0)') I; swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  WRITE (unit_tsoilprof,'(A)') TRIM(swc_line)
END IF

DO iD = 1, NDays
  IF (MOD(iD, 10) == 1 .OR. iD == NDays) THEN
    WRITE (*,'(A,I0,A,I0)') 'DAY TRACE day=', iD, ' of ', NDays
  END IF
  DVS = MAX(0, phen_state%istage - 1)
  theta_start = theta
  theta_day0 = theta
  theta_day_sum = 0.0_wp
  t_day_sum = 0.0_wp
  prev_height_state = 0.0_wp
  IF (iD > 1) prev_height_state = plant_height(iD-1)
  ! The development-stage-dependent specific leaf area, computed here because the
  ! growth development stage is set later in the day.
  ! Record 108: the reserve mass is LATCHED while the development stage is below
  ! anthesis, so after that date it holds the anthesis-time reserve, the reserve
  ! that WAS inside the fitted predictor of Bakhshandeh et al. (2012). With
  ! record 108 = 1 the argument of Eq. (6) becomes the stem mass minus the
  ! post-anthesis reserve excursion. With record 108 = 0 the height predictor is
  ! an exact copy of the stem mass and no arithmetic is performed. The
  ! record-108 = 0 argument is an EMPIRICAL PLACEHOLDER and not mechanistic
  ! height physics: its 14.43 cm mean absolute error against the measured heights
  ! is a compensating error, and record 108 = 1 worsens it to 23.62 cm while
  ! cutting the reserve-to-height-to-roughness coupling. The mechanistic
  ! successor is internode elongation and structural stem growth, records 109 to
  ! 113.
  IF (calc_dev_stage_0_2(phen_params, phen_state) < 1.0_wp) hgt_wres_anth = W_reserve
  W_stem_hgt = W_stem
  IF (height_reserve_basis == 1) W_stem_hgt = W_stem - (W_reserve - hgt_wres_anth)
  CALL update_structure(sla_from_dvs(calc_dev_stage_0_2(phen_params, phen_state)), &
    W_green, W_stem_hgt, plant_density, DVS, &
    prev_height_state, lai_growth_diag, height_growth_diag)
  ! Record 109: with the internode law on, the structure routine is NOT the
  ! height author - the freeze and the running maximum inside growth_module.f90
  ! are bypassed entirely and the height is the three-basis assembly on the
  ! REALIZED per-rank lengths. This is the MORNING read: it publishes yesterday
  ! realized state, which is exactly what the record-109 = 0 path publishes here
  ! too, since the stem mass has not moved since the previous evening. The law is
  ! ADVANCED once per day, at the evening call.
  IF (height_internode_law == 1) height_growth_diag = inl_h_b3 / 100.0_wp
    gParam_lai = lai_growth_diag
    gParam_hgt_s = height_growth_diag
    lai_sim(iD) = gParam_lai
  DOY = calendar_doy(Year(iD), Month(iD), Day(iD))
  
  s_hr_s = s_hr(iD); Tmax_d_s = Tmax_d(iD); Tmin_d_s = Tmin_d(iD); temp_s = 0.5_wp * (Tmax_d_s + Tmin_d_s)
  RH_max_s = RHmax_d(iD); RH_min_s = RHmin_d(iD); RH_s = 0.5_wp * (RH_max_s + RH_min_s); w_d_s = w_d(iD)
  WS_d_s = WS_d(iD); Irr_d_s = Irr_d(iD); dTotrad_s = Rad_d(iD); lai_s = lai_sim(iD)
  T_soil_ref = temp_s
  ! Pre-emergence development runs on the SEED-DEPTH SOIL temperature from the
  ! committed profile, node 2 at 0.04 m, and not on air temperature: the
  ! germinating seed and coleoptile sit in the soil, and McMaster et al. (2008)
  ! Ann Bot 102:561-569, p.567, state that air temperature differs from the
  ! temperature of the developing tissue.
  IF (ALLOCATED(Tsoil_prof)) THEN
    IF (SIZE(Tsoil_prof) >= 2) THEN
      t_dev_pre = Tsoil_prof(2)
    ELSE
      t_dev_pre = temp_s
    END IF
  ELSE
    t_dev_pre = temp_s
  END IF

  const_ABS = [0.0469101_wp, 0.2307534_wp, 0.5000000_wp, 0.7692465_wp, 0.9530899_wp]
  const_WGT = [0.1184635_wp, 0.2393144_wp, 0.2844444_wp, 0.2393144_wp, 0.1184635_wp]

  pr = WS_d_s / 1000.0_wp
  irr = Irr_d_s / 1000.0_wp
  IF (lai_s <= 0.0_wp) THEN
    f_pr = 1.0_wp
  ELSE
    f_pr = MAX(1.0_wp - (0.2_wp/3.0_wp) * lai_s, 0.8_wp)
  END IF
  ! Water supply and interception per FAO-56 water balance (Allen et al. 1998, Eq. 5).
  IF (irrig_mode == 1) THEN
    irr_eff = f_pr * irr
  ELSE
    irr_eff = irr
  END IF
  pr_eff = f_pr * pr + irr_eff
  pr_eff_base = pr_eff   ! D6: drip-to-soil from un-evaporated interception is added back per outer iter (mode 1)
  Ir_day = MAX(zero, (pr + irr) - pr_eff)
  Ir_mm_d = Ir_day * 1000.0_wp
  water_supply_mm_d = (pr + irr) * 1000.0_wp
  runoff_mm_d = 0.0_wp

  ! For Richards mode, allow full surface flux; no explicit cap.

  ! -------------------------------------------------------------------------- !
  ! 1. Root development, Borg & Grimes (1986), driven by the thermal-time
  !    development stage:
  !      D(dvs) = D_max [0.5 + 0.5 sin(3.03 (dvs2/2) - 1.47)]
  !    The development stage used is the current (previous-day) phenology state,
  !    so the descent runs on the development clock and not on the calendar day.
  ! -------------------------------------------------------------------------- !
  CALL update_root_depth(calc_dev_stage_0_2(phen_params, phen_state), root_max_depth, depth_lim, root_dep)

  ! -------------------------------------------------------------------------- !
  ! 2. Soil-water stress factor, from the available water in the ROOT ZONE.
  ! -------------------------------------------------------------------------- !
  wc_root = zero; wc_sat = zero; wc_fc = zero; wc_pwp = zero
  DO I = 1, NSoil_layers
    BLOCK
      REAL(wp) :: layer_thick
      REAL(wp) :: d_top, d_bot, d_root, overlap
      
      wc(I) = theta(I) * tk(I)
      layer_thick = tk(I)
      
      ! Overlap of layer with root zone
      d_top = ABS(zb(I))
      d_bot = ABS(zb(I+1))
      d_root = root_dep
      
      overlap = MAX(zero, MIN(d_bot, d_root) - d_top)
      
      IF (overlap > 0.0_wp) THEN
        wc_root = wc_root + overlap * theta(I)
        wc_sat = wc_sat + overlap * theta_0(I)
        wc_fc = wc_fc + overlap * theta_33(I)
        wc_pwp = wc_pwp + overlap * theta_1500(I)
      END IF
  END BLOCK
END DO
  
  vwc_root = wc_root / MAX(1.0e-6_wp, root_dep)
  vwc_sat = wc_sat / MAX(1.0e-6_wp, root_dep)
  vwc_fc = wc_fc / MAX(1.0e-6_wp, root_dep)
  vwc_pwp = wc_pwp / MAX(1.0e-6_wp, root_dep)
  ! Critical water content for transpiration stress (FAO-56 style).
  ! Use field capacity, not saturation.
  vwc_cr = vwc_fc - 0.55_wp * (vwc_fc - vwc_pwp)
  
  IF (.NOT. use_soil_water_stress) THEN
    f_t = 1.0_wp
  ELSE
    ! f_t (Transpiration Reduction Factor)
    ! Linear ramp from PWP to Critical Point (Feddes-type)
    IF (vwc_root >= vwc_cr) THEN
        f_t = 1.0_wp
    ELSE IF (vwc_root <= vwc_pwp) THEN
        f_t = 0.01_wp
    ELSE
        f_t = (vwc_root - vwc_pwp) / (vwc_cr - vwc_pwp)
    END IF

  END IF
  ! The root-zone-mean ramp above is retained as a DIAGNOSTIC only. Seeding the
  ! day with it deadlocks the supply loop when yesterday's mean water content sits
  ! below the ramp: a near-zero factor gives a near-zero demand and the
  ! achieved-fraction update cannot recover. The day therefore starts neutral and
  ! the outer evapotranspiration / Richards iterations converge the factor to the
  ! Feddes-consistent deliverable fraction of the CURRENT soil state.
  IF (use_soil_water_stress) f_t = 1.0_wp

  ! Record 81: the root-weighted soil water potential of the day, from the
  ! COMMITTED per-layer profile, i.e. the accepted Richards state of the previous
  ! day. All three weighting conventions are computed in both modes as
  ! diagnostics and the selected one is installed only when the record is on. The
  ! weights are those of the last hydrology call entry and the committed
  ! conductivity. With no committed profile yet, on day 1, the potential is zero
  ! and the leaf-water factor is exactly 1, flagged as not valid.
  spac_psi_s_conv = 0.0_wp; spac_psi_s_fallback = 0; spac_psi_s_valid = 0; spac_psi_s_MPa = 0.0_wp
  spac_calls_day = 0; spac_iters_max_day = 0; spac_resid_max_day = 0.0_wp
  spac_pass_on = .TRUE.
  IF (h_committed_vec_valid) THEN
    BLOCK
      REAL(wp) :: w_spac(NSoil_layers), ws_spac
      LOGICAL  :: conv_ok(3)
      INTEGER  :: kconv, jl
      conv_ok = .FALSE.
      DO kconv = 1, 3
        DO jl = 1, NSoil_layers
          SELECT CASE (kconv)
          CASE (1); w_spac(jl) = hydro_root_dist_diag(jl)
          CASE (2); w_spac(jl) = hydro_root_dist_diag(jl) * hydro_alpha_diag(jl)
          CASE (3); w_spac(jl) = hydro_root_dist_diag(jl) * k_committed_m_d(jl) / MAX(1.0e-30_wp, k_s(jl))
          END SELECT
        END DO
        ws_spac = SUM(w_spac)
        IF (ws_spac > 0.0_wp) THEN
          spac_psi_s_conv(kconv) = spac_rho_g * SUM(w_spac * h_committed_m(1:NSoil_layers)) / ws_spac
          conv_ok(kconv) = .TRUE.
        END IF
      END DO
      IF (conv_ok(1)) spac_psi_s_valid = 1
      IF (spac_mode >= 1 .AND. spac_psi_s_valid == 1) THEN
        IF (conv_ok(spac_mode)) THEN
          spac_psi_s_MPa = spac_psi_s_conv(spac_mode)
        ELSE
          spac_psi_s_MPa = spac_psi_s_conv(1); spac_psi_s_fallback = 1
        END IF
      END IF
    END BLOCK
  END IF

  ! Store initial storage for mass balance
  IF (iD == 1) THEN
    storage_init = SUM(theta * tk) * 1000.0_wp + ponding_mm
    prev_soil_storage = SUM(theta * tk) * 1000.0_wp
    storage_day1_end = storage_init
  END IF

  ! -------------------------------------------------------------------------- !
  ! 3. Canopy model: energy balance and assimilation. The soil-water stress
  !    factor reaches the gas-exchange routine through the SHARED module.
  ! -------------------------------------------------------------------------- !
  An_growth_day = 0.0_wp
  T_leaf_growth_day = 0.0_wp
  t_leaf_growth_mean = temp_s
  outer_max = 5
  tol_tleaf_out = 0.05_wp
  tleaf_prev_out = 1.0e9_wp
  ! Damp changes in potential ET between outer iterations to stabilize Richards.
  ! outer_relax in SHARED (0-1), default 0.5.
  outer_converged = .FALSE.
  eb3_outer_converged = .FALSE.
  force_final = .FALSE.
  ! Iteration 1 of each day conducts the demand pass against the frozen
  ! end-of-previous-day profile, since no committed trajectory exists yet; the
  ! actual-flux repartition validates the surface node for iterations 2 and above.
  eba_node1_valid = .FALSE.
  eb3_outer_resid_day = -1.0_wp

  theta_start = theta
  ponding_start_day = ponding_mm
  theta_day0 = theta_start

  t_i_day_prev = 0.0_wp
  e_i_day_prev = 0.0_wp
  outer_iter_loop: DO outer_iter = 1, outer_max
    f_t_iter_in = f_t   ! A2 diagnostic: the f_t this iteration's demand pass will use
    theta = theta_start
    swc_layers = theta
    ponding_mm = ponding_start_day

    ! Reset daily ledgers each outer iteration; only final iteration will accumulate.
    runoff_mm_d = 0.0_wp
    pond_evap_mm = 0.0_wp
    pond_evap_mm_d = 0.0_wp
    soil_evap_mm_d = 0.0_wp
    infil_mm_d = 0.0_wp
    infil_from_supply_mm_d = 0.0_wp
    infil_from_pond_mm_d = 0.0_wp
    top_flux_net_mm_d = 0.0_wp
    top_flux_up_mm_d = 0.0_wp
    hydro_top_flux_m_day = 0.0_wp
    hydro_top_infil_m_day = 0.0_wp
    hydro_root_sink_m_day = 0.0_wp
    hydro_bottom_flux_m_day = 0.0_wp
    q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp   ! layered-N daily aggregates
    hydro_mb_resid_day_m = 0.0_wp
    hydro_mb_correction_day_m = 0.0_wp
    ponding_min_day = ponding_mm
    ponding_max_day = ponding_mm
    ponding_sum_day = 0.0_wp
    day_dt_halvings = 0
    day_min_dt = .FALSE.
    day_failed = .FALSE.
    day_fail_events = 0
    day_fail_evap = .FALSE.

    accum_diag = force_final .OR. (outer_iter == outer_max)
    do_write = accum_diag

    ! When the post-Richards canopy pass (pass 2) will run, it is the single writer
    ! of the hourly diagnostics; the demand pass writing as well would produce two
    ! rows per day and hour in every hourly CSV.
    pass2_runs = use_soil_water_stress .AND. (.NOT. use_explicit_infiltration)
    ! Record 81: the outer-iteration-1 demand pass is the UNSTRESSED potential
    ! record 78 captures (the stress factor is 1 there by construction); with the
    ! soil-plant-atmosphere leaf it runs with a leaf-water factor of exactly 1, so
    ! the series handed to Richards stays exactly the record-78 series. Every later
    ! pass runs that leaf.
    spac_pass_on = .NOT. (spac_mode >= 1 .AND. ft_pot_unstressed == 1 .AND. outer_iter == 1)
    CALL compute_daily_canopy_fluxes(iD, unit_diurnal, (do_write .AND. .NOT. pass2_runs), &
                                     t_i_day, e_i_day, An_day, A_day, N_day, &
                                     T_leaf_day, T_air_day, T_soil_day, VPD_day, Rn_day, G_day, LE_day, H_day, &
                                     LEc_day, LEs_day, Hc_day, Hs_day, gsw_mm_day_3, Ci_day, &
                                     Ac_day, As_day, rca_day, rcs_day, rsa_day, rss_day, &
                                     gsw_day_mean, rca_daylight, rcs_daylight, rsa_daylight, rss_daylight, &
                                     tr_daylight_mm, es_daylight_mm)
    spac_pass_on = .TRUE.
    ! Record 173: the relaxation below damps the DAILY potentials but never the
    ! hourly series. This pass RAW daily transpiration potential is captured so
    ! that the series can be rescaled by the same factor and the sum of the series
    ! keeps equaling the daily value, the identity record 78 asserts on iteration
    ! 1. The scale is 1 whenever the relaxation does not fire.
    t_i_day_raw_ssc = t_i_day
    ssc_relax_scale = 1.0_wp
    IF (outer_iter > 1) THEN
      ! Relax the potential evapotranspiration between outer iterations. The previous
      ! state is used ONLY if it is finite and physically plausible (0 to 0.05 m d-1);
      ! otherwise the relaxation is skipped this iteration, the fresh value stands,
      ! and the reset is counted. No output is clamped.
      IF (t_i_day_prev >= 0.0_wp .AND. t_i_day_prev < 0.05_wp .AND. &
          e_i_day_prev >= 0.0_wp .AND. e_i_day_prev < 0.05_wp) THEN
        t_i_day = outer_relax * t_i_day + (1.0_wp - outer_relax) * t_i_day_prev
        e_i_day = outer_relax * e_i_day + (1.0_wp - outer_relax) * e_i_day_prev
      ELSE
        eb_seed_resets = eb_seed_resets + 1
      END IF
    END IF
    IF (t_i_day_raw_ssc > 1.0e-12_wp) ssc_relax_scale = t_i_day / t_i_day_raw_ssc
    t_i_day_prev = t_i_day
    e_i_day_prev = e_i_day

    ! ---------------------------------------------------------------------- !
    ! Record 78: pair the Richards demand with the water-stress factor
    ! denominator. Iteration 1 runs the canopy with a stress factor of 1, set at
    ! the top of the day, so its transpiration demand is the UNSTRESSED potential;
    ! the hourly series and its daily sum are stored, and the two are asserted
    ! equal when the record is on. On every later iteration, when the record is on,
    ! the Richards solver is forced with the unstressed series instead of the
    ! stressed one this iteration canopy pass just produced, and the sum of the
    ! stressed series is kept as a diagnostic. The capture itself runs in both
    ! modes, as a diagnostic when the record is off, where no model quantity reads
    ! it.
    ! ---------------------------------------------------------------------- !
    trans_pot_demandpass_mm = SUM(ebh_etr_dem_mm)
    IF (outer_iter == 1) THEN
      ebh_etr_dem_unstr_mm = ebh_etr_dem_mm
      trans_pot_unstressed_mm = t_i_day * 1000.0_wp
      unstr_sum_resid_mm = SUM(ebh_etr_dem_unstr_mm) - trans_pot_unstressed_mm
      IF (ft_pot_unstressed == 1) THEN
        IF (use_soil_water_stress .AND. f_t /= 1.0_wp) THEN
          WRITE (*,'(A,I0,A,G14.6)') ' ERROR A2: outer iteration 1 of day ', iD, &
            ' did not run the canopy with f_t = 1; f_t = ', f_t
          ERROR STOP 'A2: unstressed potential requires f_t = 1 on outer iteration 1'
        END IF
        ! The series and t_i_day are built by the same conversion in CANOPY_MODEL, so
        ! they agree to round-off; 1e-9 relative is a numerical tolerance, not a constant.
        IF (ABS(unstr_sum_resid_mm) > 1.0e-9_wp * MAX(1.0_wp, trans_pot_unstressed_mm)) THEN
          WRITE (*,'(A,I0,2(A,G16.8))') ' ERROR A2: day ', iD, ' SUM(unstressed series) = ', &
            SUM(ebh_etr_dem_unstr_mm), ' /= daily unstressed potential = ', trans_pot_unstressed_mm
          ERROR STOP 'A2: unstressed transpiration series does not sum to its daily potential'
        END IF
      END IF
    ELSE IF (ft_pot_unstressed == 1) THEN
      IF (spac_sink_coupling >= 1) THEN
        ! Record 173 = 1: this pass's OWN converged hourly transpiration series stays in
        ! force as the Richards sink demand - the overwrite below is not executed -
        ! rescaled by the daily relaxation factor so that hourly and daily agree.
        ebh_etr_dem_mm = ebh_etr_dem_mm * ssc_relax_scale
      ELSE
        ebh_etr_dem_mm = ebh_etr_dem_unstr_mm   ! force hydro with the UNSTRESSED hourly series
      END IF
    END IF

    ! ------------------------------------------------------------------------ !
    ! Interception: the wet canopy evaporates intercepted water at the wet-canopy
    ! potential (a vanishing canopy surface resistance), water-limited by the
    ! intercepted amount of the day. Un-evaporated water drips to the soil and is
    ! added back to infiltration, so nothing is lost silently, and the energy of
    ! the interception evaporation enters the canopy balance in the actual-flux
    ! repartition.
    ! ------------------------------------------------------------------------ !
    Ir_day_mm = Ir_mm_d
    Eintercept_act_mm    = MIN(Ir_day_mm, MAX(0.0_wp, LEc_int_pot_day))
    interception_drip_mm = MAX(0.0_wp, Ir_day_mm - Eintercept_act_mm)
    pr_eff = pr_eff_base + interception_drip_mm / 1000.0_wp   ! drip returns to infiltration

  ! -------------------------------------------------------------------------- !
  ! 3b. Surface water balance: ponding and evaporation. Ponding evaporation is
  !     handled inside the Richards sub-steps.
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  ! 4. Hydrology: the Richards equation solver.
  ! -------------------------------------------------------------------------- !
  IF (residue_evap_mode == 1) THEN
    M_res_surf_t_ha = 0.0_wp
    DO i_resv = 1, N_org_n
      ! Same DAS-from-sowing basis as fertilizer. NOTE: this
      ! locus ignores the age field tested at :3498 (no N_org_age term here).
      IF (N_org_das(i_resv) + iD_sowing <= iD) &
        M_res_surf_t_ha = M_res_surf_t_ha + N_org_mass(i_resv)
    END DO
    M_res_surf_t_ha = M_res_surf_t_ha / 1000.0_wp   ! kg/ha -> t/ha
    f_res_evap = 1.0_wp - residue_fsupp * (1.0_wp - EXP(-residue_Am_ha_t * M_res_surf_t_ha))
  ELSE
    f_res_evap = 1.0_wp
  END IF
  t_i = t_i_day; e_i = e_i_day * f_res_evap
  
  trans_act = t_i_day * 1000.0_wp ! mm/day
  evapo_act = e_i_day * 1000.0_wp * f_res_evap ! mm/day (potential, energy-balance driven; residue factor)
  trans_act_pot = trans_act
  ! A2: the repartition ratio of THIS pass's stores uses this pass's own potential
  ! (bit-identical to trans_act_pot when rec 78 = 0).
  trans_pot_pass_mm = trans_act_pot
  IF (ft_pot_unstressed == 1) THEN
    IF (spac_sink_coupling >= 1 .AND. outer_iter > 1) THEN
      ! Record 173 = 1: the demand handed to the solver and the ledger potential are
      ! THIS pass's own converged canopy transpiration, already relaxed and equal to
      ! the sum of the hourly demand series to round-off by the rescale above; the
      ! record-78 override is not applied. Outer iteration 1 is untouched, where the
      ! daily value IS the unstressed potential bit for bit.
      CONTINUE
    ELSE
      ! Record 78 = 1: the demand handed to the solver and the water-stress factor
      ! denominator are the day UNSTRESSED potential, equal to the actual
      ! transpiration on iteration 1 by construction, and the stressed value is not
      ! allowed to overwrite it.
      trans_act = trans_pot_unstressed_mm
      trans_act_pot = trans_pot_unstressed_mm
    END IF
  END IF
  evapo_act_pot = evapo_act
  Epot_soil_mm = evapo_act
  ! Snapshot the demand actually sent to the solver on this outer iteration; pass
  ! 2 later reassigns the potentials to its own end-state values, and the ledger
  ! reports both, so that the actual-and-demand pairing stays auditable.
  Epot_soil_demand_mm = Epot_soil_mm
  trans_pot_demand_mm = trans_act_pot
  ET_act = trans_act + evapo_act + Ir_mm_d
  
  avail_water = wc_root + pr_eff; avail_water_mm = avail_water * 1000.0_wp
  drain_day = 0.0_wp
  pond_evap_mm_d = 0.0_wp
  soil_evap_mm_d = 0.0_wp
  infil_mm_d = 0.0_wp
  hydro_top_flux_m_day = 0.0_wp
  hydro_top_infil_m_day = 0.0_wp
  hydro_root_sink_m_day = 0.0_wp
  hydro_bottom_flux_m_day = 0.0_wp
  q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp   ! layered-N daily aggregates
  hydro_top_req_m_day = 0.0_wp
  hydro_mb_resid_day_m = 0.0_wp
  hydro_mb_correction_day_m = 0.0_wp

  IF (use_soil_water_stress .AND. (.NOT. use_explicit_infiltration)) THEN
    BLOCK
      REAL(wp) :: pr_eff_day, evapo_pot_mm, dt_sub, actual_es_day_mm, drain_day_accum
      REAL(wp) :: ponding_start, supply_step_mm, infil_step_mm
      REAL(wp) :: es_act_mm, ec_act_mm
      REAL(wp) :: pond_evap_mm_day
      REAL(wp) :: soil_evap_mm_day
      REAL(wp) :: rain_flux_day
      REAL(wp) :: cap_active_time, flux_req_sum, flux_cap_sum, flux_applied_sum
      REAL(wp) :: k1_sum, h1_sum, theta1_sum, htop_sum, hcrit_sum
      INTEGER :: sub_i
      LOGICAL :: debug_substeps
      pr_eff_day = pr_eff
      evapo_pot_mm = Epot_soil_mm
      actual_es_day_mm = 0.0_wp
      pond_evap_mm_day = 0.0_wp
      soil_evap_mm_day = 0.0_wp
      infil_mm_d = 0.0_wp
      infil_from_supply_mm_d = 0.0_wp
      infil_from_pond_mm_d = 0.0_wp
      top_flux_net_mm_d = 0.0_wp
      top_flux_up_mm_d = 0.0_wp
      drain_day_accum = 0.0_wp
      ponding_start = ponding_mm
      theta_day_sum = 0.0_wp
      t_day_sum = 0.0_wp
      day_dt_halvings = 0
      day_min_dt = .FALSE.
      day_fail_events = 0
      day_fail_evap = .FALSE.
      day_failed = .FALSE.
      debug_substeps = accum_diag
      IF (debug_substeps) THEN
        OPEN (newunit=unit_substep, file='output/substep_diag_all.csv', &
              status='unknown', position='append', action='write')
        IF (iD == 1) THEN
          WRITE (unit_substep,'(A)') &
            'day,substep,phase,dt_sub_d,pr_eff_m_d,evapo_act_m_d,flux_top_actual_m_d,' // &
            'ponding_mm,infil_step_mm,sum_avail_mm,evap_cap_m_d,evap_pot_m_d,evap_day_m_d,' // &
            'flux_top_m_d,flux_req_m_d,flux_cap_m_d,flux_applied_m_d,cap_active,' // &
            'k1_m_d,h1_m,theta1,theta2,theta3,theta_r1,theta_r2,theta_r3,top_bc_type,h_top_m,h_crit_atm_m'
        END IF
      END IF

      ! HYDRUS-style: adaptive sub-stepping with rejection if convergence fails.
      BLOCK
        INTEGER :: n_sub, attempt, max_attempts, base_sub, max_sub
        REAL(wp) :: evap_step_mm, pond_evap_step, rem_mm
        REAL(wp) :: theta_day0(NSoil_layers), ponding_day0
        REAL(wp) :: theta_substep_start(NSoil_layers)
        LOGICAL :: day_converged
        IF ((pr_eff_day > 1.0e-9_wp) .OR. (irr_eff > 1.0e-9_wp)) THEN
          base_sub = 240
        ELSE
          base_sub = 48
        END IF
        ! Let hydro.f90 own the HYDRUS-style dt reduction locally; do not
        ! replay the full day from theta_day0 with progressively finer grids.
        max_attempts = 1
        max_sub = 4000
        theta_day0 = theta
        ponding_day0 = ponding_mm
        day_converged = .FALSE.
        DO attempt = 1, max_attempts
          ! Increase resolution (≈ HYDRUS dt/3) if previous attempt failed.
          n_sub = MIN(max_sub, base_sub * (3 ** (attempt-1)))
          dt_sub = 1.0_wp / REAL(n_sub, wp)
          rain_flux_day = pr_eff_day
          ! P2-e2: zero the delivered-evaporation bins (per attempt).
          esoil_bin_mm = 0.0_wp; pondevap_bin_mm = 0.0_wp; esoil_bin_valid = .FALSE.
          esoil_sent_mm = ebh_esoil_dem_mm * f_res_evap   ! sent-vintage snapshot
          etr_bin_mm = 0.0_wp; etr_sent_mm = ebh_etr_dem_mm   ! T2
          ! Reset state and accumulators for this attempt.
          theta = theta_day0
          swc_layers = theta
          ponding_mm = ponding_day0
          top_bc_type = 0
          h_top = 0.0_wp
          IF (ponding_mm > ponding_head_min_mm) THEN
            top_bc_type = 1
            h_top = 0.0_wp
          END IF
          actual_es_day_mm = 0.0_wp
          pond_evap_mm_day = 0.0_wp
          soil_evap_mm_day = 0.0_wp
          infil_mm_d = 0.0_wp
          infil_from_supply_mm_d = 0.0_wp
          infil_from_pond_mm_d = 0.0_wp
          top_flux_net_mm_d = 0.0_wp
          top_flux_up_mm_d = 0.0_wp
          drain_day_accum = 0.0_wp
          theta_day_sum = 0.0_wp
          t_day_sum = 0.0_wp
          day_dt_halvings = 0
          day_min_dt = .FALSE.
          day_failed = .FALSE.
          hydro_top_flux_m_day = 0.0_wp
          hydro_top_infil_m_day = 0.0_wp
          hydro_root_sink_m_day = 0.0_wp
          hydro_bottom_flux_m_day = 0.0_wp
          q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp   ! layered-N daily aggregates
          hydro_top_req_m_day = 0.0_wp
          hydro_mb_resid_day_m = 0.0_wp
          hydro_mb_correction_day_m = 0.0_wp
          cap_active_time = 0.0_wp
          flux_req_sum = 0.0_wp
          flux_cap_sum = 0.0_wp
          flux_applied_sum = 0.0_wp
          k1_sum = 0.0_wp
          h1_sum = 0.0_wp
          theta1_sum = 0.0_wp
          htop_sum = 0.0_wp
          hcrit_sum = 0.0_wp

        DO sub_i = 1, n_sub
          theta_substep_start = theta
          ! Simultaneous rain + evaporation in each sub-step.
          pr_eff = rain_flux_day

          ! Evaporate ponding first, energy-limited: if ponding exists, soil evaporation
          ! is not also demanded in the same sub-step. The window demand comes from the
          ! hourly demand series, so it carries the diurnal shape, and its daily total is
          ! the same to round-off.
          evap_step_mm = bin_series_get(ebh_esoil_dem_mm, &
            REAL(sub_i-1,wp)*dt_sub, REAL(sub_i,wp)*dt_sub) * f_res_evap
          pond_evap_step = MIN(evap_step_mm, ponding_mm)
          ponding_mm = ponding_mm - pond_evap_step
          rem_mm = evap_step_mm - pond_evap_step
          pond_evap_mm_day = pond_evap_mm_day + pond_evap_step
          IF (ponding_mm > 0.0_wp) THEN
            evapo_act = 0.0_wp
          ELSE
            IF (rem_mm > 0.0_wp) THEN
              evapo_act = rem_mm / dt_sub
            ELSE
              evapo_act = 0.0_wp
            END IF
          END IF

          ! Preserve the solver-updated atmospheric-limit head state across sub-steps. A
          ! ponded head is forced only when explicit surface water exists, and a stale
          ! ponded head is released once that surface storage is gone.
          IF (ponding_mm > ponding_head_min_mm) THEN
            top_bc_type = 1
            h_top = 0.0_wp
          ELSEIF (top_bc_type == 1 .AND. h_top >= -1.0e-9_wp) THEN
            top_bc_type = 0
            h_top = 0.0_wp
          END IF

          ! T1: window transpiration demand for the Feddes sink (mm/d over
          ! this window); hydro reads trans_act at entry (hydro.f90:107).
          IF (energy_balance_mode == 2) THEN
            trans_act = bin_series_get(ebh_etr_dem_mm, &
              REAL(sub_i-1,wp)*dt_sub, REAL(sub_i,wp)*dt_sub) / dt_sub
          END IF
          hydro_failed = .FALSE.
          CALL hydro (swc_layers, iD, dt_sub)
          hydro_fail_events_all = hydro_fail_events_all + hydro_fail_events_call
          IF (hydro_fail_events_call > 0) hydro_fail_calls_all = hydro_fail_calls_all + 1
          ! P2-e2: bin this call's delivered surface-evaporation mass.
          CALL bin_series_add(esoil_bin_mm, REAL(sub_i-1,wp)*dt_sub, &
               REAL(sub_i,wp)*dt_sub, &
               MAX(0.0_wp,(hydro_step%top_infil_m - hydro_step%top_flux_m)*1000.0_wp))
          CALL bin_series_add(pondevap_bin_mm, REAL(sub_i-1,wp)*dt_sub, &
               REAL(sub_i,wp)*dt_sub, hydro_step%pond_evap_step_mm)
          CALL bin_series_add(etr_bin_mm, REAL(sub_i-1,wp)*dt_sub, &
               REAL(sub_i,wp)*dt_sub, MAX(0.0_wp, hydro_step%root_sink_m*1000.0_wp))
          hydro_top_flux_m_day = hydro_top_flux_m_day + hydro_step%top_flux_m
          hydro_top_infil_m_day = hydro_top_infil_m_day + hydro_step%top_infil_m
          hydro_root_sink_m_day = hydro_root_sink_m_day + hydro_step%root_sink_m
          hydro_bottom_flux_m_day = hydro_bottom_flux_m_day + hydro_step%bottom_flux_m
          ! Layered-N daily aggregates: signed inter-layer water fluxes (mm) and
          ! per-layer root extraction (mm), read-only from the accepted solver state.
          q_iface_day_mm = q_iface_day_mm + hydro_q_iface_m * 1000.0_wp
          rootup_layer_day_mm = rootup_layer_day_mm + hydro_sink_layer_m * 1000.0_wp
          hydro_top_req_m_day = hydro_top_req_m_day + hydro_step%flux_top_req_m_d * dt_sub
          hydro_mb_resid_day_m = hydro_mb_resid_day_m + hydro_step%mb_resid_m
          hydro_mb_correction_day_m = hydro_mb_correction_day_m + hydro_step%mb_correction_m
          cap_active_time = cap_active_time + MERGE(dt_sub, 0.0_wp, hydro_step%cap_active == 1)
          flux_req_sum = flux_req_sum + hydro_step%flux_top_req_m_d * dt_sub
          flux_cap_sum = flux_cap_sum + hydro_step%flux_top_cap_m_d * dt_sub
          flux_applied_sum = flux_applied_sum + hydro_step%flux_top_applied_m_d * dt_sub
          k1_sum = k1_sum + hydro_step%k1_m_d * dt_sub
          h1_sum = h1_sum + hydro_step%h1_m * dt_sub
          theta1_sum = theta1_sum + hydro_step%theta1 * dt_sub
          htop_sum = htop_sum + hydro_step%h_top_m * dt_sub
          hcrit_sum = hcrit_sum + hydro_step%h_crit_atm_m * dt_sub
          drain_day_accum = drain_day_accum + drain_day
          drain_day = 0.0_wp
          day_dt_halvings = day_dt_halvings + hydro_dt_halvings
          day_fail_events = day_fail_events + hydro_fail_events_call
          day_fail_evap = day_fail_evap .OR. hydro_fail_evap
          day_min_dt = day_min_dt .OR. hydro_min_dt_reached
          IF (hydro_failed) THEN
            day_failed = .TRUE.
            EXIT
          END IF

          supply_step_mm = hydro_step%supply_step_mm
          infil_step_mm = hydro_step%infil_step_mm
          pond_evap_mm_day = pond_evap_mm_day + hydro_step%pond_evap_step_mm
          infil_from_supply_mm_d = infil_from_supply_mm_d + hydro_step%infil_from_supply_step_mm
          infil_from_pond_mm_d = infil_from_pond_mm_d + hydro_step%infil_from_pond_step_mm
          runoff_mm_d = runoff_mm_d + hydro_step%runoff_step_mm
          ponding_mm = hydro_step%ponding_end_step_mm
          infil_mm_d = infil_mm_d + hydro_step%infil_step_mm

          ponding_min_day = MIN(ponding_min_day, ponding_mm)
          ponding_max_day = MAX(ponding_max_day, ponding_mm)
          ponding_sum_day = ponding_sum_day + ponding_mm * dt_sub

          actual_es_day_mm = actual_es_day_mm + hydro_step%soil_evap_step_mm
          IF (accum_diag) THEN
            total_evap_atm_flux = total_evap_atm_flux + hydro_step%soil_evap_step_mm
          END IF
          soil_evap_mm_day = actual_es_day_mm

          IF (pr_eff > 0.0_wp .OR. ponding_mm > 0.0_wp) THEN
            top_flux_net_mm_d = top_flux_net_mm_d + hydro_step%flux_top_actual_m_d * 1000.0_wp * dt_sub
          ELSE
            top_flux_net_mm_d = top_flux_net_mm_d + MIN(0.0_wp, hydro_step%flux_top_actual_m_d) * 1000.0_wp * dt_sub
          END IF
          IF (top_bc_type == 0 .OR. (top_bc_type == 1 .AND. h_top < -1.0e-9_wp)) THEN
            top_flux_up_mm_d = top_flux_up_mm_d + MAX(0.0_wp, -hydro_step%flux_top_actual_m_d * 1000.0_wp * dt_sub)
          END IF

          IF (debug_substeps) THEN
            WRITE (unit_substep,'(I0,",",I0,",",A,26(",",G14.6))') &
              iD, sub_i, 'mix', dt_sub, pr_eff, evapo_act, flux_top_actual, ponding_mm, infil_step_mm, &
              hydro_step%sum_avail_m * 1000.0_wp, hydro_step%evap_cap_m_d, hydro_step%evap_rate_pot_m_d, &
              hydro_step%evap_rate_day_m_d, hydro_flux_top_m_d, hydro_step%flux_top_req_m_d, hydro_step%flux_top_cap_m_d, &
              hydro_step%flux_top_applied_m_d, REAL(hydro_step%cap_active,wp), hydro_step%k1_m_d, hydro_step%h1_m, &
              theta(1), theta(2), theta(3), theta_r(1), theta_r(2), theta_r(3), &
              REAL(top_bc_type,wp), hydro_step%h_top_m, hydro_step%h_crit_atm_m
          END IF

          theta_day_sum = theta_day_sum + 0.5_wp * (theta_substep_start + theta) * dt_sub
          t_day_sum = t_day_sum + dt_sub
          END DO
          IF (.NOT. day_failed) THEN
            day_converged = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. day_converged) THEN
          ! Reject day: restore state and mark as failed (HYDRUS-style).
          theta = theta_day0
          swc_layers = theta
          ponding_mm = ponding_day0
          actual_es_day_mm = 0.0_wp
          pond_evap_mm_day = 0.0_wp
          soil_evap_mm_day = 0.0_wp
          infil_mm_d = 0.0_wp
          infil_from_supply_mm_d = 0.0_wp
          infil_from_pond_mm_d = 0.0_wp
          top_flux_net_mm_d = 0.0_wp
          top_flux_up_mm_d = 0.0_wp
          drain_day_accum = 0.0_wp
          ! Layered-N: a rejected (held-state) day moved no water -> no N advection.
          q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp
          theta_day_sum = theta
          t_day_sum = 1.0_wp
          IF (.NOT. day_failed) day_failed = .TRUE.
          IF (accum_diag) THEN
            hydro_fail_count = hydro_fail_count + 1
            hydro_fail_days(hydro_fail_count) = iD
            IF (day_fail_evap) hydro_fail_evap_days = hydro_fail_evap_days + 1
          END IF
        END IF
        IF (debug_substeps) THEN
          CLOSE(unit_substep)
        END IF
      END BLOCK

      pr_eff = pr_eff_day
      hydro_top_req_m_day = flux_req_sum
      drain_mm_d = drain_day_accum * 1000.0_wp
      pond_evap_mm = pond_evap_mm_day
      pond_evap_mm_d = pond_evap_mm_day
      IF (t_day_sum > 0.0_wp) THEN
        ponding_mean_day = ponding_sum_day / t_day_sum
      ELSE
        ponding_mean_day = ponding_mm
      END IF
      ! Actual soil evaporation from the solver's own mass integrals: the gross
      ! infiltration minus the net top flux IS the total upward surface extraction,
      ! by definition, on every boundary-condition path, rain falling at the same
      ! time included. This form is exact and needs no gating; a daily net-upward
      ! proxy reports about zero on a day with both rain and evaporation, and a gated
      ! per-sub-step accumulator under-counts under the head boundary condition.
      actual_es_day_mm = MAX(0.0_wp, (hydro_top_infil_m_day - hydro_top_flux_m_day) * 1000.0_wp)
      soil_evap_mm_d = actual_es_day_mm
      es_act_mm = actual_es_day_mm
      ec_act_mm = hydro_root_sink_m_day * 1000.0_wp
      evapo_act = es_act_mm
      trans_act = ec_act_mm
      ! Reconnect stomatal water stress to the ACTUAL supply the solver just
      ! delivered, instead of a once-per-day pre-hydrology value. The Richards sink
      ! is proportional to the demand, S = D F with F the deliverable fraction of the
      ! current soil state, so the achieved ratio of actual to potential
      ! transpiration IS F, independently of the stress factor, and the
      ! supply-consistent stomatal factor is therefore F itself - a multiplicative
      ! update would accumulate powers of F and crash to the floor. The next outer
      ! iteration, and pass 2, then run the canopy with supply-consistent stomata,
      ! closing assimilation, conductance, transpiration and canopy temperature
      ! together. No new constants: the floor 0.01 and the cap 1.0 are the existing
      ! bounds of the factor, and Feddes-limited uptake remains the sole
      ! authoritative supply.
      IF (use_soil_water_stress .AND. trans_act_pot > 1.0e-9_wp) THEN
        f_t = MAX(0.01_wp, MIN(1.0_wp, trans_act / trans_act_pot))
      END IF
      IF (energy_balance_mode >= 1 .AND. .NOT. pass2_runs) THEN
        ! When pass 2 runs, which is the released configuration, the repartition moves
        ! AFTER it and consumes the pass-2 stores; this block remains for
        ! configurations without the post-Richards canopy pass.
        ! The actual-flux repartition redistributes the daily actual fluxes over the
        ! hours by potential weight and recomputes the three-node temperatures with the
        ! ACTUAL latent heat, so that the canopy is warmer when transpiration is
        ! water-limited and the latent deficit becomes sensible heat.
        BLOCK
          REAL(wp) :: f_trans, f_int, f_soil
          f_trans = 0.0_wp; f_int = 0.0_wp; f_soil = 0.0_wp
          ! P2-e2: the committed soil latent uses the BINNED Richards mass, not
          ! f_soil x potential; f_soil remains a diagnostic ratio only.
          esoil_bin_valid = (energy_balance_mode == 2)
          IF (trans_pot_pass_mm > 1.0e-9_wp) f_trans = MIN(1.0_wp, trans_act / trans_pot_pass_mm)
          IF (Epot_soil_mm  > 1.0e-9_wp) f_soil  = MIN(1.0_wp, evapo_act / Epot_soil_mm)
          IF (LEc_int_pot_day > 1.0e-9_wp) f_int  = MIN(1.0_wp, Eintercept_act_mm / LEc_int_pot_day)
          ! The accumulate-and-write flag is true only on the accepted final outer
          ! iteration, so the soil-heat profile advances, and its diagnostics are
          ! written, exactly once per day, the profile being held during trial
          ! iterations.
          CALL repartition_day_actual(f_trans, f_int, f_soil, accum_diag)
        ! A3b: the REPORTED daily EB ledger is the COMMITTED (actual) state; the
        ! demand-state values are retained for the *_demand diagnostics.
        Rn_demand_day = Rn_day; G_demand_day = G_day
        LEc_demand_day = LEc_day; LEs_demand_day = LEs_day
        Hc_demand_day = Hc_day; Hs_demand_day = Hs_day
        Rn_day = Rn_day_act; G_day = G_day_act
        LEc_day = LEc_day_act; LEs_day = LEs_day_act
        Hc_day = Hc_day_act; Hs_day = Hs_day_act
        LE_day = LEc_day + LEs_day; H_day = Hc_day + Hs_day
        IF (accum_diag) CALL write_hourly_eb_actual(unit_eb_diag, iD)
        IF (accum_diag) CALL write_hourly_radiation_committed(iD)
        IF (accum_diag) CALL write_hourly_tleaf_classes(iD)      ! (rec 95), ON-only
        IF (accum_diag .AND. node_pair_open .AND. energy_balance_mode == 2) THEN
          BLOCK
            INTEGER :: hh
            DO hh = 1, eb_nsteps_s
              WRITE (unit_node_pair,'(I0,",",I0,2(",",G14.6))') iD, hh, node_armed_row(hh), eba_node1(hh)
            END DO
          END BLOCK
        END IF
        IF (accum_diag .AND. esoil_pair_open .AND. energy_balance_mode == 2) THEN
          BLOCK
            INTEGER :: hh
            REAL(wp) :: L_hh, dtj
            DO hh = 1, eb_nsteps_s
              L_hh = (2500.8_wp - 2.36_wp*ebh_Tair(hh) + 0.0016_wp*ebh_Tair(hh)**2) * 1000.0_wp
              dtj = ebh_weight(hh) * 86400.0_wp
              WRITE (unit_esoil_pair,'(I0,",",I0,7(",",G14.6))') iD, hh, &
                esoil_sent_mm(hh), esoil_bin_mm(hh), &
                esoil_bin_mm(hh) * L_hh / dtj, ebh_LEs_pot(hh), &
                etr_sent_mm(hh), etr_bin_mm(hh), etr_bin_mm(hh) * L_hh / dtj
            END DO
          END BLOCK
        END IF
        IF (accum_diag) WRITE (unit_eb_demand,'(I0,6(",",G14.6))') iD, &
          Rn_demand_day, G_demand_day, LEc_demand_day, LEs_demand_day, Hc_demand_day, Hs_demand_day
          ! actual latent daily totals (mm) for ET reconciliation / diagnostics.
          LEc_trans_act_day = trans_act
          LEc_int_act_day   = Eintercept_act_mm
          LEs_act_day       = evapo_act
        END BLOCK
      END IF
      ! Soil-sourced evapotranspiration only (transpiration plus soil evaporation
      ! plus pond evaporation) for the SOIL water ledger. Wet-canopy interception
      ! evaporation leaves the CANOPY store and not the soil column, so it must not
      ! enter the soil mass balance, where it would be counted twice, its water
      ! having already left through reduced infiltration. Total evapotranspiration
      ! including interception is reported separately.
      ET_act = trans_act + evapo_act + pond_evap_mm_d

      BLOCK
        REAL(wp) :: cap_active_frac, flux_req_mean, flux_cap_mean, flux_applied_mean
        REAL(wp) :: k1_mean, h1_mean, theta1_mean, htop_mean, hcrit_mean
        IF (t_day_sum > 0.0_wp) THEN
          cap_active_frac = cap_active_time / t_day_sum
          flux_req_mean = flux_req_sum / t_day_sum
          flux_cap_mean = flux_cap_sum / t_day_sum
          flux_applied_mean = flux_applied_sum / t_day_sum
          k1_mean = k1_sum / t_day_sum
          h1_mean = h1_sum / t_day_sum
          theta1_mean = theta1_sum / t_day_sum
          htop_mean = htop_sum / t_day_sum
          hcrit_mean = hcrit_sum / t_day_sum
        ELSE
          cap_active_frac = 0.0_wp
          flux_req_mean = 0.0_wp
          flux_cap_mean = 0.0_wp
          flux_applied_mean = 0.0_wp
          k1_mean = 0.0_wp
          h1_mean = 0.0_wp
          theta1_mean = 0.0_wp
          htop_mean = 0.0_wp
          hcrit_mean = 0.0_wp
        END IF
        IF (accum_diag) THEN
          WRITE (unit_hydrodiag,'(I0,",",L1,",",I0,",",L1,",",I0,10(",",G14.6),11(",",G14.6))') iD, &
            day_failed, day_dt_halvings, day_min_dt, MERGE(1,0,use_explicit_infiltration), &
            flux_top_actual, ponding_mm, REAL(top_bc_type,wp), &
            hydro_q34_m_d, hydro_k34_m_d, hydro_top_flux_m_day, hydro_top_infil_m_day, &
            hydro_root_sink_m_day, hydro_bottom_flux_m_day, hydro_mb_resid_day_m, &
            cap_active_frac, flux_req_mean, flux_cap_mean, flux_applied_mean, k1_mean, h1_mean, theta1_mean, htop_mean, &
            hcrit_mean, hydro_qmax_evap_m_d, REAL(day_fail_events, wp)
          ! P2-d: commit the accepted surface head for the next day's demand pass.
          h1_committed_m = hydro_h1_m
          h1_committed_valid = .TRUE.
          ! rec 81 (review C1): committed per-layer head and K, same vintage as h1_committed_m.
          h_committed_m(1:NSoil_layers) = r_h(1:NSoil_layers)
          k_committed_m_d(1:NSoil_layers) = r_k(1:NSoil_layers)
          h_committed_vec_valid = .TRUE.
        END IF
      END BLOCK
      ! HYDRUS-style surface ledger:
      ! Use sub-step accumulated infiltration (infil_mm_d); do not overwrite with per-call hydro_top_infil_m.
      delta_ponding_mm = ponding_mm - ponding_start_day
      shortfall_mm = MAX(0.0_wp, (pr_eff_day * 1000.0_wp + irr_eff * 1000.0_wp) - infil_mm_d - delta_ponding_mm)
      IF (accum_diag) THEN
        total_shortfall_mm = total_shortfall_mm + shortfall_mm
      END IF
    END BLOCK
    IF (day_failed) THEN
      ! Reject hydrology day after Richards failure: restore state and zero fluxes.
      theta = theta_day0
      swc_layers = theta
      ! Preserve the failed day supply in surface storage first, then spill only
      ! the excess above the allowed ponding store to runoff.
      ponding_mm = MIN(ponding_max_mm, ponding_start_day + water_supply_mm_d)
      top_bc_type = 0
      h_top = 0.0_wp
      IF (ponding_start_day > ponding_head_min_mm) THEN
        top_bc_type = 1
        h_top = 0.0_wp
      END IF
      IF (ponding_mm > ponding_head_min_mm) THEN
        top_bc_type = 1
        h_top = 0.0_wp
      END IF
      runoff_mm_d = MAX(0.0_wp, ponding_start_day + water_supply_mm_d - ponding_max_mm)
      drain_mm_d = 0.0_wp
      infil_mm_d = 0.0_wp
      infil_from_supply_mm_d = 0.0_wp
      infil_from_pond_mm_d = 0.0_wp
      top_flux_net_mm_d = 0.0_wp
      top_flux_up_mm_d = 0.0_wp
      hydro_top_flux_m_day = 0.0_wp
      hydro_top_infil_m_day = 0.0_wp
      hydro_root_sink_m_day = 0.0_wp
      hydro_bottom_flux_m_day = 0.0_wp
      q_iface_day_mm = 0.0_wp; rootup_layer_day_mm = 0.0_wp   ! layered-N daily aggregates
      hydro_top_req_m_day = 0.0_wp
      hydro_mb_resid_day_m = 0.0_wp
      hydro_mb_correction_day_m = 0.0_wp
      delta_ponding_mm = 0.0_wp
      pond_evap_mm = 0.0_wp
      pond_evap_mm_d = 0.0_wp
      soil_evap_mm_d = 0.0_wp
      evapo_act = 0.0_wp
      trans_act = 0.0_wp
      ET_act = 0.0_wp
      theta_day_sum = theta
      t_day_sum = 1.0_wp
      ponding_min_day = ponding_mm
      ponding_max_day = ponding_mm
      ponding_sum_day = 0.0_wp
    END IF
  ELSE
    drain_mm_d = 0.0_wp
    soil_evap_mm_d = evapo_act
    pond_evap_mm_d = 0.0_wp
    ET_act = trans_act + evapo_act + pond_evap_mm_d
  END IF
  delta_ponding_mm = ponding_mm - ponding_start_day
  IF (t_day_sum <= 0.0_wp) THEN
    theta_day_sum = theta
    t_day_sum = 1.0_wp
  END IF

  ! -------------------------------------------------------------------------- !
  ! 4a. Recompute the canopy energy balance with the actual soil evaporation, by
  !     actual-to-potential flux scaling of the committed energy ledger. This
  !     updates the leaf temperature, the assimilation and the energy-balance
  !     partitioning without altering the hydrologic fluxes already applied to
  !     the soil.
  ! -------------------------------------------------------------------------- !
  IF (use_soil_water_stress .AND. (.NOT. use_explicit_infiltration)) THEN
    trans_act_used = trans_act
    CALL compute_daily_canopy_fluxes(iD, unit_diurnal, do_write, t_i_day, e_i_day, An_day, A_day, N_day, &
                                     T_leaf_day, T_air_day, T_soil_day, VPD_day, Rn_day, G_day, LE_day, H_day, &
                                     LEc_day, LEs_day, Hc_day, Hs_day, gsw_mm_day_3, Ci_day, &
                                     Ac_day, As_day, rca_day, rcs_day, rsa_day, rss_day, &
                                     gsw_day_mean, rca_daylight, rcs_daylight, rsa_daylight, rss_daylight, &
                                     tr_daylight_mm, es_daylight_mm)
    ! Pass 2 is the authoritative pass. Its potential totals, on the post-hydrology
    ! water content and the supply-consistent stress factor, are captured before
    ! the restore, and the actual-flux repartition then runs on ITS hourly stores,
    ! so that the reported committed canopy and soil temperatures, the soil-heat
    ! advance and the outer-loop convergence variable are all computed from one
    ! consistent hourly state with the correct per-sub-step solar geometry.
    ! Pass 2 ran with the stressed stress factor; its potential scales its own
    ! stores in the repartition below and is reported as the stressed-potential
    ! diagnostic. Under record 78 = 1 it must not overwrite the unstressed
    ! potential, which is the stress-factor denominator and the water ledger.
    trans_pot_pass_mm = t_i_day * 1000.0_wp
    trans_pot_stressed_mm = trans_pot_pass_mm
    IF (ft_pot_unstressed == 0) trans_act_pot = trans_pot_pass_mm
    Epot_soil_mm  = e_i_day * 1000.0_wp * f_res_evap
    trans_act = trans_act_used
    t_i_day = trans_act / 1000.0_wp
    e_i_day = evapo_act / 1000.0_wp
    IF (energy_balance_mode >= 1) THEN
      BLOCK
        REAL(wp) :: f_trans, f_int, f_soil
        f_trans = 0.0_wp; f_int = 0.0_wp; f_soil = 0.0_wp
        ! P2-e2: committed soil latent from the BINNED Richards mass (see above).
        esoil_bin_valid = (energy_balance_mode == 2)
        IF (trans_pot_pass_mm > 1.0e-9_wp) f_trans = MIN(1.0_wp, trans_act / trans_pot_pass_mm)
        IF (Epot_soil_mm  > 1.0e-9_wp) f_soil  = MIN(1.0_wp, evapo_act / Epot_soil_mm)
        IF (LEc_int_pot_day > 1.0e-9_wp) f_int  = MIN(1.0_wp, Eintercept_act_mm / LEc_int_pot_day)
        CALL repartition_day_actual(f_trans, f_int, f_soil, accum_diag)
        ! A3b: the REPORTED daily EB ledger is the COMMITTED (actual) state; the
        ! demand-state values are retained for the *_demand diagnostics.
        Rn_demand_day = Rn_day; G_demand_day = G_day
        LEc_demand_day = LEc_day; LEs_demand_day = LEs_day
        Hc_demand_day = Hc_day; Hs_demand_day = Hs_day
        Rn_day = Rn_day_act; G_day = G_day_act
        LEc_day = LEc_day_act; LEs_day = LEs_day_act
        Hc_day = Hc_day_act; Hs_day = Hs_day_act
        LE_day = LEc_day + LEs_day; H_day = Hc_day + Hs_day
        IF (accum_diag) CALL write_hourly_eb_actual(unit_eb_diag, iD)
        IF (accum_diag) CALL write_hourly_radiation_committed(iD)
        IF (accum_diag) CALL write_hourly_tleaf_classes(iD)      ! (rec 95), ON-only
        IF (accum_diag .AND. node_pair_open .AND. energy_balance_mode == 2) THEN
          BLOCK
            INTEGER :: hh
            DO hh = 1, eb_nsteps_s
              WRITE (unit_node_pair,'(I0,",",I0,2(",",G14.6))') iD, hh, node_armed_row(hh), eba_node1(hh)
            END DO
          END BLOCK
        END IF
        IF (accum_diag .AND. esoil_pair_open .AND. energy_balance_mode == 2) THEN
          BLOCK
            INTEGER :: hh
            REAL(wp) :: L_hh, dtj
            DO hh = 1, eb_nsteps_s
              L_hh = (2500.8_wp - 2.36_wp*ebh_Tair(hh) + 0.0016_wp*ebh_Tair(hh)**2) * 1000.0_wp
              dtj = ebh_weight(hh) * 86400.0_wp
              WRITE (unit_esoil_pair,'(I0,",",I0,7(",",G14.6))') iD, hh, &
                esoil_sent_mm(hh), esoil_bin_mm(hh), &
                esoil_bin_mm(hh) * L_hh / dtj, ebh_LEs_pot(hh), &
                etr_sent_mm(hh), etr_bin_mm(hh), etr_bin_mm(hh) * L_hh / dtj
            END DO
          END BLOCK
        END IF
        IF (accum_diag) WRITE (unit_eb_demand,'(I0,6(",",G14.6))') iD, &
          Rn_demand_day, G_demand_day, LEc_demand_day, LEs_demand_day, Hc_demand_day, Hs_demand_day
        LEc_trans_act_day = trans_act
        LEc_int_act_day   = Eintercept_act_mm
        LEs_act_day       = evapo_act
      END BLOCK
    END IF
  END IF

  ! Add the pond evaporation to the reported soil and surface latent-heat ledger.
  ! The final canopy recompute scales the soil latent heat to soil evaporation
  ! only, and the pond evaporation is accounted for separately in the water
  ! balance, so it is added here to keep the reported latent and sensible heat
  ! consistent with the final daily evapotranspiration terms.
  BLOCK
    REAL(wp) :: lambda_pond_mj_mm, pond_le_day, tair_mean_c
    ! The daily mean air temperature is still the 24-hour SUM at this point: the
    ! canopy routine returns the Gauss-weighted sum and the division by 24 happens
    ! further below. The latent-heat-of-vaporization polynomial is therefore
    ! normalized locally rather than by moving that division, so that no other
    ! consumer of the accumulator changes. Every other call site of the same
    ! polynomial passes an instantaneous temperature.
    tair_mean_c = T_air_day / 24.0_wp
    lambda_pond_mj_mm = (2500.8_wp - 2.36_wp * tair_mean_c + 0.0016_wp * tair_mean_c**2) / 1000.0_wp
    pond_le_day = pond_evap_mm_d * lambda_pond_mj_mm
    LEs_day = LEs_day + pond_le_day
    LE_day = LEc_day + LEs_day
    Hs_day = Rn_day - G_day - LE_day - Hc_day
    H_day = Hc_day + Hs_day
  END BLOCK

  ! mode 1 converges the outer loop on ACTUAL Tcan (post-repartition);
  ! mode 0 (incumbent) keeps its potential daily-mean Tleaf criterion.
  tleaf_mean_out = Tcan_act_day
  ! A2 diagnostic census row (both modes): this outer iteration's f_t pairing.
  WRITE (unit_ft_outer,'(I0,",",I0,",",I0,9(",",G14.6),",",L1,",",L1)') iD, outer_iter, &
    ft_pot_unstressed, f_t_iter_in, f_t, hydro_Fraw_diag, trans_act, trans_pot_demand_mm, &
    trans_pot_unstressed_mm, trans_pot_demandpass_mm, trans_pot_stressed_mm, unstr_sum_resid_mm, &
    accum_diag, (ABS(tleaf_mean_out - tleaf_prev_out) < tol_tleaf_out)
  eb3_outer_iters_day = outer_iter
  ! COUPLING: record the outer fixed-point residual (|dTcan_act| between outer
  ! iterations, K) for the convergence census; -1 = no previous iterate yet.
  IF (tleaf_prev_out < 1.0e8_wp) eb3_outer_resid_day = ABS(tleaf_mean_out - tleaf_prev_out)
  IF (ABS(tleaf_mean_out - tleaf_prev_out) < tol_tleaf_out) THEN
    IF (.NOT. accum_diag) THEN
      force_final = .TRUE.
      tleaf_prev_out = tleaf_mean_out
      CYCLE outer_iter_loop
    END IF
    outer_converged = .TRUE.
    eb3_outer_converged = .TRUE.
    EXIT outer_iter_loop
  END IF
  tleaf_prev_out = tleaf_mean_out
  IF (accum_diag) EXIT outer_iter_loop
  END DO outer_iter_loop
  ! R9: season convergence census of the outer (canopy<->Richards) fixed point.
  outer_days_total = outer_days_total + 1
  IF (.NOT. outer_converged) outer_days_unconverged = outer_days_unconverged + 1
  IF (eb3_outer_resid_day > outer_resid_max_season) outer_resid_max_season = eb3_outer_resid_day
  ! Carbon-reserve diagnostics: one row per day, written after the outer loop.
  IF (carbon_reserve_mode == 1) THEN
    IF (iD == 1) THEN
      OPEN (NEWUNIT=unit_cres, FILE='output/daily_carbon_reserve.csv', STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit_cres,'(A)') 'day,reserve_C_pool,daily_C_to_grain_from_current_assim,' // &
        'daily_C_to_grain_from_reserve,grain_sink_demand,grain_sink_unfilled,realized_TKW,' // &
        'reserve_respired,ledger_residual'
    END IF
    WRITE (unit_cres,'(I0,8(",",G14.6))') iD, W_reserve, cres_to_grain_assim, cres_to_grain_res, &
      cres_sink_demand, cres_sink_unfilled, &
      1000.0_wp*W_grain/MAX(kernel_number_m2, 1.0_wp), &
      cres_respired, cres_ledger_resid
  END IF
  ! Record 92 only: everything needed to test the reserve law against Schnyder
  ! (1993) p.238 without re-deriving anything - the labeled pool, its share of
  ! the stem, what was deposited and what the capacity and the deposition close
  ! refused, the draw, and the day organ-state closure residual.
  IF (stem_wsc_mech >= 1) THEN
    IF (iD == 1) THEN
      OPEN (NEWUNIT=unit_wsc, FILE='output/daily_stem_wsc.csv', STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit_wsc,'(A)') 'day,dap_anth,W_reserve_g_m2,W_stem_g_m2,wsc_frac_of_stem,' // &
        'deposited_today_g_m2,refused_cap_cum_g_m2,refused_window_cum_g_m2,' // &
        'draw_today_g_m2,gain_today_g_m2,closure_resid_g_m2_d'
    END IF
    WRITE (unit_wsc,'(I0,",",I0,9(",",G14.6))') iD, &
      MERGE(iD - wsc_day_anth, -999, wsc_day_anth > 0), &
      W_reserve, W_stem, W_reserve/MAX(W_stem, 1.0e-9_wp), &
      wsc_dep_day, wsc_refused_cap, wsc_refused_win, &
      cres_to_grain_res/MAX(reserve_conv_eff, 1.0e-9_wp), cres_to_grain_res, &
      cres_ledger_resid
  END IF
  ! Record 109 only: everything needed to recompute the windows, the length
  ! profile, the modeled lineic mass and the height assembly from the file
  ! itself. THE SHORTFALL PRINTED HERE IS A MODEL-INTERNAL quantity - the
  ! difference between a bracketed literature length potential and what this
  ! model's own carbon paid for. It is not a measured discrepancy, it is not
  ! evidence of a defect in the axis count, tillering, photosynthesis or
  ! allocation, and it licenses no parameter change.
  IF (height_internode_law == 1) THEN
    IF (.NOT. inl_diag_open) THEN
      OPEN (NEWUNIT=unit_inl, FILE='output/internode_ranks.csv', STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit_inl,'(A)') 'day,tau_Cd,HS,nl_ms,axes_live_m2,window_member,length_member,' // &
        'lineic_member,ear_member,onset_arm,dI_Cd,ir,abs_rank,tau_on_Cd,tau_end_Cd,w_open,' // &
        'L_f_cm,L_cm,dL_pot_cm,dL_act_cm,dL_short_cm,lambda_mg_cm,m_struct_mg_axis,' // &
        'lambda_mod_mg_cm,bind_flag,dM_struct_g_m2,dM_offer_g_m2,dM_len_g_m2,e8_resid_g_m2,' // &
        'h_b1_cm,h_b2_cm,h_b3_cm,L_ear_cm'
      inl_diag_open = .TRUE.
    END IF
    DO i_inl = 0, INL_MAXIR
      IF (inl_Lf_use(i_inl) <= 0.0_wp) CYCLE
      SELECT CASE (inl_bind(i_inl))
      CASE (1) ; inl_bind_txt = 'POT'
      CASE (2) ; inl_bind_txt = 'MASS'
      CASE (3) ; inl_bind_txt = 'CAP'
      CASE DEFAULT ; inl_bind_txt = 'CLOSED'
      END SELECT
      WRITE (unit_inl,'(I0,2(",",G14.6),",",I0,",",G14.6,5(",",I0),",",G14.6,",",I0,",",I0,' // &
        '2(",",G14.6),",",I0,8(",",G14.6),",",A,8(",",G14.6))') &
        iD, clr_tau_Cd, clr_HS, MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1), &
        inl_axes, hgt_int_window_member, hgt_int_length_member, hgt_int_lineic_member, &
        hgt_int_ear_member, inl_onset_arm, inl_dI_use, i_inl, inl_absrank(i_inl), &
        inl_tau_on(i_inl), inl_tau_end(i_inl), MERGE(1, 0, inl_bind(i_inl) > 0), &
        inl_Lf_use(i_inl), inl_L(i_inl), inl_dLpot(i_inl), inl_dLact(i_inl), &
        inl_dLpot(i_inl) - inl_dLact(i_inl), inl_lam_use(i_inl), inl_m(i_inl), &
        inl_m(i_inl) / MAX(inl_L(i_inl), 1.0e-12_wp), inl_bind_txt, &
        inl_dM_struct, inl_dM_offer, inl_dM_len, inl_e8_resid, &
        inl_h_b1, inl_h_b2, inl_h_b3, inl_Lear_cm
    END DO
  END IF
  ! Record 94 only: everything needed to falsify the law without re-deriving
  ! anything - the zone depth the model built for itself, the flat top-of-canopy
  ! nitrogen per unit area against the measured maximum 3.85 g N m-2
  ! (Silva-Perez et al. 2020 J Exp Bot 71(7):2299-2311 Fig. 4B p.2305), the
  ! canopy mean it must always exceed, and BOTH fractions of absorbed PAR.
  IF (canopy_n_flag_zone >= 1) THEN
    IF (iD == 1) THEN
      OPEN (NEWUNIT=unit_canE, FILE='output/daily_canopy_n_profile.csv', STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit_canE,'(A)') 'day,dvs,LAI,Fc_flag_gai,Fc_over_LAI,N_leaf_g_m2,' // &
        'n_mean_g_m2,n_top_g_m2,n_top_over_mean,over_3p85,fPAR_top,fPAR_eff,fPAR_rel_err'
    END IF
    WRITE (unit_canE,'(I0,12(",",G14.6))') iD, dvs, gParam_lai, cls_gai_flag, &
      cls_gai_flag/MAX(gParam_lai, 1.0e-9_wp), N_leaf, &
      N_leaf/MAX(gParam_lai, 1.0e-9_wp), canE_n_top, &
      canE_n_top/MAX(N_leaf/MAX(gParam_lai, 1.0e-9_wp), 1.0e-9_wp), &
      MERGE(1.0_wp, 0.0_wp, canE_n_top > 3.85_wp .AND. gParam_lai > 0.5_wp), &
      canE_fPAR_top, canE_fPAR_eff, &
      (canE_fPAR_top - canE_fPAR_eff)/MAX(canE_fPAR_eff, 1.0e-9_wp)
  END IF
  ! Residue-evap diagnostics: one row per day, written after the outer loop.
  IF (residue_evap_mode == 1) THEN
    IF (iD == 1) THEN
      OPEN (NEWUNIT=unit_resevap, FILE='output/daily_residue_evap.csv', STATUS='REPLACE', ACTION='WRITE')
      WRITE (unit_resevap,'(A)') 'day,M_res_t_ha,f_res_evap'
    END IF
    WRITE (unit_resevap,'(I0,2(",",G14.6))') iD, M_res_surf_t_ha, f_res_evap
  END IF

  ! write daily actual-flux EB diagnostics (mode 1). All latent components
  ! reported in mm; energy state from the actual-flux repartition.
  IF (energy_balance_mode >= 1 .AND. eb_act_open) THEN
    WRITE (unit_eb_act,'(I0,26(",",G14.6))') &
      iD, Ir_day_mm, Eintercept_act_mm, interception_drip_mm, &
      trans_act_pot, LEc_trans_act_day, LEc_int_pot_day, LEc_int_act_day, Epot_soil_mm, LEs_act_day, &
      trans_act, evapo_act, (trans_act + evapo_act + pond_evap_mm_d + Eintercept_act_mm), pond_evap_mm_d, &
      ! The air-temperature column: this write precedes the division of the 24-hour
      ! accumulator by 24, so the value is normalized inline here, as the leaf
      ! temperature on the same line already was.
      (T_leaf_day / 24.0_wp), Tcan_act_day, Tac_act_day, Tsoil_act_day, Hc_act_day, Hs_act_day, &
      (T_air_day / 24.0_wp), &
      eb_canopy_resid_day, REAL(eb3_outer_iters_day, wp), MERGE(1.0_wp, 0.0_wp, eb3_outer_converged), &
      eb3_outer_resid_day, Epot_soil_demand_mm, trans_pot_demand_mm
  END IF

  ! NOTE: WB-override ET disabled to keep reported ET consistent with Richards.

  ! -------------------------------------------------------------------------- !
  ! 4b. Full-season soil profile diagnostics, per layer and per day.
  ! -------------------------------------------------------------------------- !
  ! ======================================================================== !
  ! DIAGNOSTIC-ONLY: this block re-derives, for printing, the interface fluxes
  ! hydro.f90 assembled, and no solver statement is changed. Both terms are the
  ! solver's own expressions: the gravity term carries the sign the solver uses
  ! with the depth coordinate positive downward, and the inter-node conductivity
  ! is the geometric or arithmetic mean the active averaging mode selects. The
  ! published flux columns of output/daily_soil_profile_diag.csv are therefore
  ! the model's own fluxes. The water balance, the water contents, the ledgers and
  ! every gate are unaffected: they never read these two columns.
  ! ======================================================================== !
  BLOCK
    REAL(wp) :: th_s, th_r, Se, q_in, q_out, K_mid, dz_mid, eps_len
    INTEGER :: k
    eps_len = 1.0e-8_wp
    DO k = 1, NSoil_layers
      th_s = MAX(1.0e-8_wp, theta_0(k))
      th_r = MAX(0.0_wp, theta_r(k))
      Se = (theta(k) - th_r) / MAX(1.0e-8_wp, (th_s - th_r))
      Se = MIN(1.0_wp, MAX(0.0_wp, Se))

      ! Flux into layer k (positive downward). R8: solver expression,
      ! hydro.f90:710-724 (k_avg_mode switch + `- 1.0_wp` gravity term).
      IF (k == 1) THEN
        q_in = flux_top_actual
      ELSE
        K_mid = 0.5_wp * (r_k(k) + r_k(k-1))
        dz_mid = MAX(eps_len, ABS(dz(k) - dz(k-1)))
        q_in = -1.0_wp * K_mid * ((r_h(k) - r_h(k-1)) / dz_mid - 1.0_wp)
      END IF

      ! Flux out of layer k (positive downward). R8: solver expression,
      ! hydro.f90:728-743 (k_avg_mode switch + `- 1.0_wp` gravity term).
      IF (k < NSoil_layers) THEN
        K_mid = 0.5_wp * (r_k(k) + r_k(k+1))
        dz_mid = MAX(eps_len, ABS(dz(k+1) - dz(k)))
        q_out = -1.0_wp * K_mid * ((r_h(k+1) - r_h(k)) / dz_mid - 1.0_wp)
      ELSE
        IF (drainage_type == 1) THEN
          q_out = r_k(NSoil_layers)
        ELSE
          q_out = 0.0_wp
        END IF
      END IF

      WRITE (unit_soilprof,'(I0,",",I0,9(",",G14.6))') iD, k, theta(k), th_s, th_r, r_h(k), r_k(k), k_s(k), Se, q_in, q_out
    END DO
  END BLOCK

  ! -------------------------------------------------------------------------- !
  ! 5. Mass balance and reporting.
  ! -------------------------------------------------------------------------- !
  raw_soil_storage_mm = SUM(theta * tk) * 1000.0_wp
  hydro_mb_resid_mm = hydro_mb_resid_day_m * 1000.0_wp
  hydro_mb_correction_mm = hydro_mb_correction_day_m * 1000.0_wp  ! rec 149 wb_surface_limit=2 ON-only (0 under modes 0/1)
  ! Preserve raw Richards solution before any balance correction.
  theta_raw = theta

  ! Carry forward the raw Richards solution directly.
  swc_layers = theta

  soil_storage_mm = SUM(theta * tk) * 1000.0_wp
  storage_mm = soil_storage_mm + ponding_mm
  delta_storage = soil_storage_mm - prev_soil_storage
  delta_storage_total = storage_mm - (prev_soil_storage + ponding_start_day)
  ! Storage-implied top flux after conservative correction.
  flux_top_mass_m_day = (soil_storage_mm - prev_soil_storage) / 1000.0_wp + &
                        hydro_root_sink_m_day + hydro_bottom_flux_m_day
  flux_top_mismatch_mm_d = (hydro_top_flux_m_day - flux_top_mass_m_day) * 1000.0_wp
  IF (day_failed .AND. ABS(delta_storage) > 1.0e-6_wp) THEN
    WRITE (*,'(A,I0,2(A,G14.6))') 'Water-balance check: storage change not closed on day ', iD, &
      ' delta_storage_mm=', delta_storage, ' soil_storage_mm=', soil_storage_mm
  END IF
  soil_mb_resid_mm = (pr_eff * 1000.0_wp - trans_act - evapo_act - drain_mm_d - runoff_mm_d) - delta_storage
  ! The effective precipitation already includes irrigation; it is not counted
  ! twice. DRAINAGE belongs in this ledger: bottom drainage leaves the soil
  ! column exactly as runoff and evaporation do, and the sibling ledger one line
  ! above subtracts it for the same control volume.
  wb_resid_mm_d = (pr_eff * 1000.0_wp) - runoff_mm_d - drain_mm_d &
                  - (soil_evap_mm_d + pond_evap_mm_d + trans_act) - delta_storage_total
  total_wb_resid = total_wb_resid + wb_resid_mm_d
  ! Mass balance with explicit soil and surface storage (HYDRUS-style ledger).
  total_pond_evap = total_pond_evap + pond_evap_mm
  ! Use hydro-integrated infiltration for soil balance to avoid ledger drift.
  ! Soil balance directly from Richards-integrated fluxes (HYDRUS-style).
  balance_soil_mm = (soil_storage_mm - prev_soil_storage) - &
                    (flux_top_mass_m_day * 1000.0_wp - hydro_root_sink_m_day * 1000.0_wp - hydro_bottom_flux_m_day * 1000.0_wp)
  ! The surface balance is driven strictly by the net top flux of the Richards
  ! solver rather than by a ledgered infiltration.
  ! Record 149 at or above 1, a reported quantity only: the gross supply is
  ! precipitation plus irrigation, while the real inflow to the surface and
  ! ponding state is the interception-net supply. Using the gross term makes this
  ! print equal the negative cumulative interception evaporation every season, a
  ! print defect and not a real surface loss. Mode 0 keeps the gross-supply
  ! print.
  IF (wb_surface_limit >= 1) THEN
    balance_surface_mm = (ponding_mm - ponding_start_day) - &
                         (pr_eff * 1000.0_wp - pond_evap_mm - runoff_mm_d - hydro_top_infil_m_day * 1000.0_wp)
  ELSE
    balance_surface_mm = (ponding_mm - ponding_start_day) - &
                         (water_supply_mm_d - pond_evap_mm - runoff_mm_d - hydro_top_infil_m_day * 1000.0_wp)
  END IF
  balance_mm = balance_soil_mm + balance_surface_mm
  balance_no_ir = balance_mm
  prev_soil_storage = soil_storage_mm
  IF (iD == 1) THEN
    storage_day1_end = storage_mm
  END IF

  ! Record 149 mode 2: append the on-only mass-correction field (the 26th),
  ! matching the mode-2 header above. Modes 0 and 1 write the original 25-field
  ! row.
    WRITE (unit_wbdiag,'(I0,24(",",G14.6))') iD, pr_eff * 1000.0_wp, irr_eff * 1000.0_wp, &
      runoff_mm_d, ponding_mm, delta_ponding_mm, ponding_start_day, ponding_min_day, ponding_max_day, ponding_mean_day, &
      soil_storage_mm, delta_storage, evapo_act, trans_act, ET_act, &
      hydro_top_flux_m_day * 1000.0_wp, flux_top_mass_m_day * 1000.0_wp, flux_top_mismatch_mm_d, &
      hydro_root_sink_m_day * 1000.0_wp, hydro_bottom_flux_m_day * 1000.0_wp, &
      soil_mb_resid_mm, balance_soil_mm, balance_surface_mm, balance_mm, wb_resid_mm_d

  ! One-day debug print for hydrology ledger (day 1).
  IF (iD == 1) THEN
    WRITE (*,'(A)') 'Mass-balance check, day 1:'
    WRITE (*,'(A,1X,G14.6)') ' delta_storage_mm=', delta_storage
    WRITE (*,'(A,1X,G14.6)') ' hydro_top_flux_mm=', hydro_top_flux_m_day * 1000.0_wp
    WRITE (*,'(A,1X,G14.6)') ' hydro_root_sink_mm=', hydro_root_sink_m_day * 1000.0_wp
    WRITE (*,'(A,1X,G14.6)') ' hydro_bottom_flux_mm=', hydro_bottom_flux_m_day * 1000.0_wp
    WRITE (*,'(A,1X,G14.6)') ' balance_soil_mm=', balance_soil_mm
    WRITE (*,'(A,1X,G14.6)') ' balance_surface_mm=', balance_surface_mm
    WRITE (*,'(A,1X,G14.6)') ' soil_mb_resid_mm=', soil_mb_resid_mm
    WRITE (*,'(A,1X,G14.6)') ' flux_top_actual_mm_d=', flux_top_actual * 1000.0_wp
    WRITE (*,'(A,1X,G14.6)') ' ponding_mm=', ponding_mm
    WRITE (*,'(A,1X,G14.6)') ' runoff_mm_d=', runoff_mm_d
  END IF

  VPD_day = VPD_day / 24.0_wp; T_leaf_day = T_leaf_day / 24.0_wp; T_air_day = T_air_day / 24.0_wp
  T_soil_day = T_soil_day / 24.0_wp
  ! Post-emergence development runs on the COMMITTED CANOPY temperature daily
  ! mean - the quantity validated against the infrared measurements, with a
  ! six-date mean bias of 0.03-0.10 K - and not on air temperature. The basis is
  ! McMaster et al. (2008) p.567: air temperature usually differs from shoot apex
  ! or canopy temperature. The apex-specific arbitration could not be obtained
  ! from the two sources that carry it. The day falls back to air temperature on
  ! any day the canopy state is not yet defined. The phenology was calibrated
  ! against air temperature, and the demand-pass canopy temperature, which runs
  ! several K warmer at low leaf area index, accelerates the stages.
  IF (Tcan_act_day > -90.0_wp .AND. lai_growth_diag > 0.0_wp) THEN
    t_dev_post = Tcan_act_day
  ELSE
    t_dev_post = T_air_day
  END IF
  CALL update_phenology_day(phen_params, phen_state, Year(iD), Month(iD), Day(iD), Tmax_d_s, Tmin_d_s, t_dev_pre, t_dev_post)
  ! Record 89: open the spike growth period of Fischer (1985) on his own measured
  ! duration, 300 degree-Cd above 4.5 degC before anthesis (p.449). The
  ! remaining development is EXACT - the anthesis threshold of the phenology
  ! module minus the accumulated development units, both a cultivar input and a
  ! past state - and only the conversion of that remainder into degree-Cd is
  ! projected, from the season-to-date weather. It must precede the allocation,
  ! whose partition reads the window flag, and it returns on its first line when
  ! the record is off.
  BLOCK
    REAL(wp) :: tdu_rem_gnx, dtdu_gnx
    tdu_rem_gnx = MAX(0.0_wp, phen_params%stage_thresh(6) - phen_state%cum_tdu)  ! stage 6 = Anthesis
    dtdu_gnx    = MAX(0.0_wp, phen_state%cum_tdu - cum_tdu_prev_gnx)
    cum_tdu_prev_gnx = phen_state%cum_tdu
    CALL update_spike_window(iD, T_air_day, tdu_rem_gnx, dtdu_gnx)
  END BLOCK
  pheno_stage(iD) = MAX(0, phen_state%istage - 1)
  DVS = pheno_stage(iD)
  ! Expose the day thermal time to the organ-nitrogen laws. The increment is the
  ! positive part of the development increment from the existing phenology clock
  ! (base temperature 0 C), i.e. the SAME thermal-time definition as the
  ! first-order rate fit of Bertheloot, Martre & Andrieu (2008) Plant Physiol
  ! 148:1707-1720 (p.1718: daily degree-days above 0 C from the mean of the
  ! maximum and minimum air temperature), as the degradation and synthesis rates
  ! of the organ turnover package, and as the grain-nitrogen sink of record 86
  ! and its rate-limited supply. The assignment is therefore made whenever any of
  ! those laws is live, and skipped when none is.
  IF (leafN_drawdown == 1 .OR. grainN_partition == 1 .OR. grainN_sink_mech >= 1) &
    dTT_day_N = MAX(0.0_wp, phen_state%current_dtt)
  ! Record 64: advance the per-culm leaf-rank chronology on the SAME day anchored
  ! phenology state. Placed immediately after the phenology update because it
  ! consumes that call outputs and nothing else: the stage index (1 = emergence
  ! starts the clock and 4 = booting freezes the emergent final leaf number) and
  ! the day degree-Cd on base 0 C. It returns on its first line when the record
  ! is off.
  ! ================================================================== !
  ! Records 161 and 162: THE CANOPY CLOCK DAILY INCREMENT, and nothing else in
  ! the model, is formed here.
  !
  ! Why here and only here: update_canopy_leaf_ranks is the ONLY consumer that
  ! takes a clock argument. update_canopy_blade_area and update_canopy_tillering
  ! take none - they read the module state this call sets - so leaf appearance,
  ! the final-leaf-number freeze, the blade-expansion ramp and the tiller-birth
  ! windows ALL follow from this one argument. The phenology stage clock, the
  ! accumulated development units, the growth development stage and the
  ! partition are untouched.
  !
  ! At both records 0 the increment is the screen-air development increment and
  ! nothing else.
  !
  ! Record 161 = 1: the day degree-days are recomputed from the TOP SOIL LAYER
  ! temperature through the SAME degree-day law and the SAME cultivar cardinal
  ! temperatures, for the window Abichou (2016) printed - soil temperature before
  ! stem extension and air temperature thereafter (p.23, justification p.38). The
  ! window is closed on the model's own nearest stated marker, TERMINAL
  ! SPIKELET, which is a substitution for the start of stem extension; the
  ! thermal test is written out as well as the marker flag so that the window
  ! also closes under record 150 = 0, where that flag is never set. The
  ! allocation and size guards fall back to the air increment rather than to an
  ! undefined value.
  !
  ! Record 162 = 1: the increment becomes the vernalization- and
  ! photoperiod-gated development increment (McMaster et al. 2008 Eq. (1) p.563,
  ! already coded in phenology_module.f90). With record 161 at 0 the EXACT field
  ! is passed and nothing is reconstructed. With record 161 on, the gate factor
  ! is recovered as the ratio of the gated to the raw increment and applied to
  ! the soil-driven increment; on a day with a zero raw increment that ratio is
  ! undefined, so the gate is rebuilt from the two factors under the record-150
  ! window, which is safe because the development units cannot advance on such a
  ! day and so the marker cannot have flipped during it.
  ! ================================================================== !
  clk_dtt_canopy = phen_state%current_dtt
  IF (canopy_clock_driver >= 1) THEN
    IF (phen_state%emerge_done .AND. .NOT. phen_state%ts_done .AND. &
        (phen_state%cum_tdu - phen_state%tdu_at_emerge) < TS_TDU_MCMASTER) THEN
      IF (ALLOCATED(Tsoil_prof)) THEN
        IF (SIZE(Tsoil_prof) >= 1) THEN
          clk_dtt_canopy = calc_dtt(Tsoil_prof(1), phen_params%t_base, phen_params%t_opt_lower, &
                                    phen_params%t_opt_upper, phen_params%t_max, phen_params%gdd_method)
        END IF
      END IF
    END IF
  END IF
  IF (canopy_clock_basis >= 1) THEN
    IF (canopy_clock_driver >= 1) THEN
      IF (phen_state%current_dtt > 0.0_wp) THEN
        clk_gate_canopy = phen_state%current_dtdu / phen_state%current_dtt
      ELSE IF (phen_params%dev_factor_window == 1) THEN
        IF (phen_state%emerge_done .AND. .NOT. phen_state%ts_done) THEN
          clk_gate_canopy = MIN(phen_state%current_vf, phen_state%current_df)
        ELSE
          clk_gate_canopy = 1.0_wp
        END IF
      ELSE
        clk_gate_canopy = phen_state%current_vf * phen_state%current_df
      END IF
      clk_dtt_canopy = clk_dtt_canopy * clk_gate_canopy
    ELSE
      clk_dtt_canopy = phen_state%current_dtdu
    END IF
  END IF
  ! Record 174: the canopy-clock reading at which the anchored booting stage will
  ! be stated, formed here and only here, from the same two phenology numbers
  ! the stage engine itself uses - the cultivar file booting threshold and the
  ! accumulated development units recorded on the emergence day. It is a pure
  ! diagnostic hand-off, consumed by the chronology on the emergence day alone,
  ! and it is assigned on every day regardless of the record, which costs nothing
  ! and leaves the record-174 = 0 path unchanged, since the prediction stays zero
  ! and nothing reads the target.
  IF (phen_state%emerge_done) THEN
    clr_tau_boot_target = phen_params%stage_thresh(4) - phen_state%tdu_at_emerge
  END IF
  CALL update_canopy_leaf_ranks(pheno_stage(iD), clk_dtt_canopy, iD)
  ! Record 65: the blade dimensions and the expansion read the SAME day rank
  ! chronology, so this must follow the record-64 call. It returns on its first
  ! lines when either record is off.
  CALL update_canopy_blade_area()
  ! Record 66: the axis population reads the SAME day chronology (the windows of
  ! opportunity, the cohort birth stages, the flag ligulation) and the SAME day
  ! per-axis areas, so it must follow both calls above. It is passed the model
  ! own green leaf area index, the light state the bud-outgrowth law is evaluated
  ! on, and the plant density.
  CALL update_canopy_tillering(gParam_lai, plant_density)
  ! Record 98: form the spike DEMAND for today, here, because the cohort
  ! densities are now today's ones and the allocation that spends the demand is
  ! still below. It returns with a zero demand when the record is off or the
  ! Fischer window is shut.
  CALL bgn_set_spike_demand(T_air_day)
  ! Record 83: the water-arrested expansion accumulates on the SAME day rank
  ! geometry and cohort births and must precede the senescence step, which reads
  ! it.
  CALL update_leaf_expansion_water()
  ! Record 144: the carbon-coupled bonus accumulator reads TODAY'S water factor,
  ! just set above, and YESTERDAY'S lamina-sink carbon surplus, set at the end of
  ! the lamina split inside the allocation far below, which is the one-day lag;
  ! its update is read by TOMORROW'S record-65 and record-66 calls, the potential
  ! areas being already fixed for today.
  CALL update_leaf_carbon_bonus()
  ! Record 122: the pre-commit nitrogen-arrest increment of every rank, on the
  ! SAME day geometry, must precede the senescence step, which reads it, and is
  ! committed at the paid fraction after the allocation.
  CALL update_lamina_n_arrest_prov()
  ! Record 67: senescence reads the SAME day chronology, per-rank areas and axis
  ! densities, so it must follow all three calls above.
  CALL update_canopy_leaf_senesc()
  ! (flag canopy_sheath_ear, rec 68): the non-lamina standing
  ! area reads S4's Sl and S3's densities and writes only the cse_* block.
  CALL update_canopy_sheath_ear(pheno_stage(iD))
  ! Record 100 mode 4: form today's bounded per-axis structural stem mass demand
  ! here, after today's chronology and axis population are final and before the
  ! allocation spends it. It returns on its first line when record 100 is not 4.
  CALL sdm4_set_demand(phen_state%current_dtt)
  ! Record 69: publish the sink state the rest of the day reads. It must run
  ! AFTER the whole record-64-to-68 chain, because the area it publishes has to
  ! be carbon-blind.
  CALL update_canopy_sink_state()
  CALL bootstrap_emergence_state(DVS, plant_density, phen_params%init_green_g_plant, phen_params%init_stem_g_plant, &
    phen_params%init_root_g_plant, phen_params%init_storage_g_plant, W_green, W_dead, W_stem, W_root, W_storage)
  ! The last argument passes the bootstrap-day specific leaf area (m2 g-1) so
  ! that the nitrogen bootstrap can bound the initial leaf nitrogen by the
  ! capacity times the initial leaf area index, with that index being the
  ! specific leaf area times the green lamina mass, which is the model's own leaf-
  ! area identity. The specific leaf area comes from the existing development-
  ! stage profile at the current stage.
  CALL bootstrap_nitrogen_state(DVS, phen_params%init_leafN_frac, phen_params%init_stemN_frac, &
    phen_params%init_rootN_frac, phen_params%init_storageN_frac, W_green, W_dead, W_stem, W_root, W_storage, &
    N_leaf, N_dead_leaf, N_stem, N_root, N_storage, sla_from_dvs(calc_dev_stage_0_2(phen_params, phen_state)))
    IF (An_day < 0.0_wp) An_day = 0.0_wp
    ! Maintenance respiration is driven by the mean AIR temperature (Teh 2006
    ! Eqs. 7.3-7.4). The demand-pass canopy temperature can reach 29 degC at low
    ! leaf area index, because limited transpiration and a small boundary-layer
    ! resistance prevent canopy cooling, and using it inflates the maintenance cost
    ! about 2.5-fold and traps the model in a state where maintenance exceeds gross
    ! assimilation.
    t_leaf_growth_mean = T_air_day
    An_growth_day = An_day
    N_assim_biochem_day = MAX(0.0_wp, N_day)
  dvs2_growth = calc_dev_stage_0_2(phen_params, phen_state)
  dvs2_rate_growth = calc_dev_stage_rate_0_2(phen_params, phen_state)
  ! The days-after-anthesis counter for the living-root capacity decline (Zheng
  ! et al. 2022 Front Plant Sci 12:814658 Tables 2-3; the development stage 1
  ! marks anthesis). The crossing day itself is day 0, where the capacity is
  ! exactly 1, and the counter advances on each subsequent day. The basis is the
  ! calendar day after anthesis because that paper reports no site temperatures
  ! or degree-days, so no thermal-time form can be written.
  IF (dvs2_growth >= 1.0_wp) THEN
    IF (anth_crossed_L2c) daa_root = daa_root + 1
    anth_crossed_L2c = .TRUE.
  END IF
  N_fert_day = 0.0_wp
  N_mineral_day = 0.0_wp
  N_leach_d = 0.0_wp
  N_soil_uptake_day = 0.0_wp
  N_uptake_massflow_day = 0.0_wp
  N_uptake_active_day = 0.0_wp
  S_active_pot_day = 0.0_wp; S_massflow_pot_day = 0.0_wp; U_pot_day = 0.0_wp  ! diagnostics reset
  N_assim_biochem_retained_day = 0.0_wp
  N_demand_day = 0.0_wp
  ! rec 155: day-reset of the ON-only uptake-regulation diagnostics.
  ! New variables only -- nothing printed or read at mode 0, so this is inert there.
  nur_EN_day = 1.0_wp ; nur_Nmob_day = 0.0_wp ; nur_cap_nmax_day = 0.0_wp
  nur_cap_nc_day = 0.0_wp ; nur_dem_before = 0.0_wp ; nur_dem_after = 0.0_wp
  nur_D_max_day = 0.0_wp
  f_Tr_N = 0.0_wp
  f_CO2_N = 0.0_wp
  N_remob_leaf_day = 0.0_wp
  N_remob_stem_day = 0.0_wp
  N_balance_resid = 0.0_wp
  ! The dead stem and chaff nitrogen joins the daily nitrogen closure so that the
  ! stem and chaff death flux is mass-conserving. It is identically zero when
  ! record 63 is off.
  N_plant_before = N_leaf + N_dead_leaf + N_dead_stem + N_stem + N_root + N_storage + N_NO3_root + N_NO3_shoot &
                   + N_mobile_veg   ! P4: mobile vegetative pool in the closure

  ! -------------------------------------------------------------------- !
  ! The maximum rate of de novo nitrogen supply to the chloroplast. Busch, Sage &
  ! Farquhar (2018) Nat Plants 4:46-54, p.49: it can be limited by the nitrogen
  ! availability, by the demand for nitrogen in the form of glycine and serine,
  ! or by the activities of the enzymes involved in the reduction of nitrate to
  ! ammonium (their Fig. 3b and Eqs. 16-21). With an explicit shoot nitrate pool
  ! it is the lesser of the enzymatic ceiling 1.21 umol m-2 s-1, their Fig. 3b
  ! value at non-limiting nitrate, and the rate at which the shoot pool can
  ! supply nitrate over the photoperiod.
  ! The unit chain from a ground-area pool in g N m-2 to a leaf-area rate in
  ! umol m-2 s-1 divides by the molar mass 14 g mol-1, multiplies by 1e6, and
  ! divides by the daylight seconds and by the leaf area index, because the
  ! quantity is leaf-area based in assim.f90 and the canopy re-multiplies by the
  ! leaf area index. It feeds the NEXT day canopy, since the canopy is computed
  ! earlier in the day: an explicit one-day lag. The CO2 sensitivity of the shoot
  ! reduction enters through the photorespiratory ratio.
  ! -------------------------------------------------------------------- !
  ! ==================================================================== !
  ! Record 158: WHICH QUANTITY IS THE SUPPLY TERM. Mode 0 is the STANDING
  ! vacuolar store divided by the photoperiod and the leaf area index, i.e. a
  ! stock asserted to be a flux. Mode 1 replaces the NUMERATOR ONLY by the model
  ! own daily de novo nitrate-reduction DEMAND in the shoot, set at the end of
  ! the PREVIOUS day, the same one-day lag. The unit chain, the operand order and
  ! the enzymatic ceiling are unchanged, so mode 0 performs the identical
  ! arithmetic on the identical operand.
  ! ==================================================================== !
  BLOCK
    REAL(wp) :: daylight_s_nmax, lai_nmax, N_max_supply
    daylight_s_nmax = MAX(1.0_wp, Daylength * 3600.0_wp)
    lai_nmax        = MAX(0.1_wp, gParam_lai)
    IF (nitrate_drain_rate_basis >= 1) THEN
      N_max_supply  = (ndb_rate_prev / 14.0_wp) * 1.0e6_wp / daylight_s_nmax / lai_nmax  ! umol m-2leaf s-1
    ELSE
    N_max_supply    = (N_NO3_shoot / 14.0_wp) * 1.0e6_wp / daylight_s_nmax / lai_nmax  ! umol m-2leaf s-1
    END IF
    N_max_param     = MIN(N_max_ref, N_max_supply)   ! Busch ceiling vs shoot-pool supply
    ! Record 158 ledgers: both bases side by side on the SAME day, so that the size
    ! of the change is readable per deck. Nothing reads or prints them at mode 0.
    IF (nitrate_drain_rate_basis >= 1) THEN
      ndb_nmax_store_day = (N_NO3_shoot  / 14.0_wp) * 1.0e6_wp / daylight_s_nmax / lai_nmax
      ndb_nmax_rate_day  = (ndb_rate_prev / 14.0_wp) * 1.0e6_wp / daylight_s_nmax / lai_nmax
      ndb_days_live      = ndb_days_live + 1
      IF (N_max_ref <= N_max_supply) ndb_days_cap = ndb_days_cap + 1
      ndb_nmax_used_cum  = ndb_nmax_used_cum  + N_max_param
      ndb_nmax_store_cum = ndb_nmax_store_cum + ndb_nmax_store_day
      ndb_rate_cum       = ndb_rate_cum       + ndb_rate_prev
    END IF
  END BLOCK

  ! Fertilizer applications from the input/N_fertilizer.csv schedule: the
  ! nitrogen enters the soil mineral pool on the specified day after sowing, and
  ! the availability timing follows the fertilizer type and the soil conditions
  ! the deck specifies in that file.
  ! -------------------------------------------------------------------- !
  ! Soil nitrogen (urea, ammonium, nitrate, residue and microbial) before today
  ! external inputs, for the closure.
  soilN_before_day = N_urea_soil + N_NH4_soil + TAN_surf + N_soil_min + res_fast_N + res_slow_N + soil_mic_N
  ! Record 93 at or above 1: the humus pool is INSIDE the ledger on both sides
  ! and humus mineralization stops counting as an external input, because it is
  ! an internal transfer from humus to ammonium. A dropped decrement of the pool,
  ! or a flux that does not match the decrement, then surfaces in the balance
  ! residual instead of hiding in it. Under record 93 = 0 the sum is the same as
  ! the line above.
  IF (soil_n_mech >= 1) soilN_before_day = soilN_before_day + NHUM_soil
  N_hydrolysis_day = 0.0_wp; N_nitrification_day = 0.0_wp
  decomp_C_day = 0.0_wp; decomp_N_day = 0.0_wp; respired_C_day = 0.0_wp
  net_min_NH4_day = 0.0_wp; immob_NH4_day = 0.0_wp; immob_NO3_day = 0.0_wp
  turnover_N_day = 0.0_wp; residue_N_in_day = 0.0_wp
  ! (1) The fertilizer is split into urea, ammonium and nitrate pools by the form
  !     fractions of the deck (urea, ammonium nitrate, calcium ammonium nitrate,
  !     ammonium sulfate, potassium nitrate). input/N_fertilizer.csv documents
  !     the dates as days after SOWING, so the application day is the day after
  !     sowing plus the sowing model day of record 88.
  DO I_fert = 1, N_fert_n_apps
    IF (iD == N_fert_DAS_arr(I_fert) + iD_sowing .AND. N_fert_applied < I_fert) THEN
      N_urea_soil = N_urea_soil + N_fert_N_arr(I_fert) * N_fert_furea_arr(I_fert)
      ! Surface-applied ammonium-form fertilizer enters the SURFACE ammoniacal pool,
      ! the emission layer, and not the bulk pool, because ammonia volatilization is
      ! a surface process (Sommer et al. 2004 Adv Agron 82:557, sec. III.B).
      ! Nitrate-form nitrogen goes to the bulk, since nitrate does not volatilize.
      TAN_surf    = TAN_surf    + N_fert_N_arr(I_fert) * N_fert_fNH4_arr(I_fert)
      IF (soilN_layered == 1) THEN
        ! Layered: surface-applied NO3-form fertilizer enters LAYER 1
        ! (design A: TAN/urea surface semantics preserved; NO3 -> NO3_1).
        NO3_lay(1) = NO3_lay(1) + N_fert_N_arr(I_fert) * N_fert_fNO3_arr(I_fert)
        N_soil_min = N_soil_min + N_fert_N_arr(I_fert) * N_fert_fNO3_arr(I_fert)  ! aggregate mirror
      ELSE
      N_soil_min  = N_soil_min  + N_fert_N_arr(I_fert) * N_fert_fNO3_arr(I_fert)
      END IF
      N_fert_day  = N_fert_day  + N_fert_N_arr(I_fert)
      N_fert_applied = I_fert
    END IF
  END DO

  ! -------------------------------------------------------------------- !
  ! Nitrogen leaching from the root zone with the drainage water. Burns (1975)
  ! J Soil Sci 26:255-270, Eq. 1, piston-flow nitrate transport:
  !   N_leach = N_soil drain / (theta_rz z_root 1000)
  ! with the bottom drainage from the Richards solver, the mean volumetric water
  ! content of the root zone, and the current root depth. Ramos et al. (2002)
  ! Agric Water Manag 55:97-118, Eq. 2 and Table 3, report piston-flow leaching
  ! agreeing with lysimeter data to within 10-15 % for coarse-textured soils
  ! under similar irrigation regimes.
  ! -------------------------------------------------------------------- !
  ! (2) Mineralization (ammonification): background soil organic nitrogen to
  !     ammonium, and not directly to nitrate. The background rate is read from
  ! driver_N.csv and is site-specific; the typical range for agricultural soils
  ! is 0.005-0.03 g N m-2 d-1 (Stanford & Smith 1972 SSSAJ 36:465-472).
  ! -------------------------------------------------------------------- !
  ! Record 93 at or above 1: that flat, site-specific rate is REPLACED by humus
  ! mineralization drawn out of a real depleting pool. Brisson et al. (1998)
  ! Agronomie 18(5-6):311-346, doi:10.1051/agro:19980501, read at first hand from
  ! 200-dpi page renders:
  !   (31) p.341  V(z)    = Vpot(z) FTH(z) FH(z)
  !   (32) p.342  Vpot(z) = K2POT NHUM
  !   (42) p.342  FTH(z)  = exp(0.115 (Tsoil(z) - 15))
  !   (44) p.342  FH(z)   = 0.2 + 0.8 (HUR - HUMIN)/(HUCC - HUMIN)
  ! This is a construction ported from another model. The humus pool is
  ! distributed uniformly over the existing mineralization depth, so summing
  ! Eq. (31) over the profile gives the pool times the rate constant times the
  ! temperature response times the thickness-weighted mean water response. The
  ! mineralized nitrogen is then routed to ammonium by the existing
  ! depth-weighting code below: the record changes the RATE and not where it
  ! lands.
  ! Two guards, neither of them a new constant: (a) Eq. (42) is unbounded, so the
  ! 45 degC cap on the temperature entering the response is reused from the
  ! residue temperature cap in soiln_module.f90; (b) Eq. (44) is capped at 1,
  ! because in the source model its water term is the content of a tipping bucket
  ! that never exceeds field capacity, while a Richards water content CAN exceed
  ! it transiently after rain or irrigation, and the cap keeps the equation
  ! inside the domain it was published for.
  ! -------------------------------------------------------------------- !
  humus_min_day = 0.0_wp; humus_fth_day = 0.0_wp; humus_fh_day = 0.0_wp
  IF (soil_n_mech >= 1) THEN
    BLOCK
      REAL(wp) :: zt_h, ov_h, ovsum_h, fth_h, fh_h, fhw_h, tdec_h
      INTEGER  :: iLh
      tdec_h = MIN(T_soil_day, soiln_Tdec_cap_c)
      fth_h  = EXP(stics_fth_slope * (tdec_h - stics_fth_tref))
      zt_h = 0.0_wp; ovsum_h = 0.0_wp; fhw_h = 0.0_wp
      DO iLh = 1, NSoil_layers
        ov_h = MAX(0.0_wp, MIN(zt_h + tk(iLh), n_min_depth_m) - zt_h)
        IF (ov_h > 0.0_wp) THEN
          fh_h = stics_fh_min + stics_fh_span * (theta(iLh) - theta_1500(iLh)) / &
                 MAX(1.0e-6_wp, theta_33(iLh) - theta_1500(iLh))
          fh_h = MIN(1.0_wp, MAX(stics_fh_min, fh_h))
          fhw_h   = fhw_h   + fh_h * ov_h
          ovsum_h = ovsum_h + ov_h
        END IF
        zt_h = zt_h + tk(iLh)
        IF (zt_h >= n_min_depth_m) EXIT
      END DO
      IF (ovsum_h > 1.0e-9_wp) THEN
        humus_fh_day = fhw_h / ovsum_h
      ELSE
        humus_fh_day = 1.0_wp
      END IF
      ! rec 93 = 2: ATTRIBUTION ARM -- same pool and same K2POT, responses off.
      IF (soil_n_mech == 2) THEN
        fth_h = 1.0_wp; humus_fh_day = 1.0_wp
      END IF
      humus_fth_day = fth_h
      humus_min_day = MIN(NHUM_soil, &
                          MAX(0.0_wp, humus_k2pot * NHUM_soil * fth_h * humus_fh_day))
      NHUM_soil       = NHUM_soil       - humus_min_day
      humus_min_total = humus_min_total + humus_min_day
    END BLOCK
    N_mineral_day = humus_min_day
  ELSE
  N_mineral_day = N_mineral_rate
  END IF
  IF (soilN_layered == 1) THEN
    ! Layered: background mineralization to the ammonium of the layers spanning the
    ! mineralization depth (a driver row, on the basis that the measured organic
    ! matter is a topsoil value), weighted by each layer overlap with that span.
    ! The constant itself is unchanged.
    BLOCK
      REAL(wp) :: zt, ov, ovsum, ovl(NSoil_layers)
      INTEGER  :: iLm
      zt = 0.0_wp; ovsum = 0.0_wp; ovl = 0.0_wp
      DO iLm = 1, NSoil_layers
        ov = MAX(0.0_wp, MIN(zt + tk(iLm), n_min_depth_m) - zt)
        ovl(iLm) = ov; ovsum = ovsum + ov
        zt = zt + tk(iLm)
        IF (zt >= n_min_depth_m) EXIT
      END DO
      IF (ovsum > 1.0e-9_wp) THEN
        DO iLm = 1, NSoil_layers
          IF (ovl(iLm) > 0.0_wp) NH4_lay(iLm) = NH4_lay(iLm) + N_mineral_day * ovl(iLm) / ovsum
        END DO
      ELSE
        NH4_lay(1) = NH4_lay(1) + N_mineral_day
      END IF
      N_NH4_soil = N_NH4_soil + N_mineral_day   ! aggregate mirror
    END BLOCK
  ELSE
  N_NH4_soil = N_NH4_soil + N_mineral_day
  END IF

  ! (2b) Organic residue treatments. In-season events are added on their day
  !      after sowing, then decomposition, immobilization and microbial turnover
  !      run daily. Ordered BEFORE hydrolysis and nitrification so that
  !      high-carbon residues can immobilize fresh ammonium before uptake and
  !      leaching.
  IF (org_treatments_active) THEN
    DO i_org = 1, N_org_n
      ! An in-season residue event needs an already-applied guard, like the
      ! fertilizer block above: without one, an event dated at sowing re-fires every
      ! day while the days after planting are still zero, and the deck receives many
      ! times its residue nitrogen, silently inflating the soil nitrogen economy of
      ! every run with organic treatments. The guard marks the event applied by
      ! setting its age to one day; the pre-sowing branch runs only at
      ! initialization. The basis is the same day-after-sowing one the fertilizer
      ! uses.
      IF (N_org_age(i_org) <= 0 .AND. iD == N_org_das(i_org) + iD_sowing) THEN
        BLOCK
          REAL(wp) :: nin
          CALL apply_residue_event(i_org, nin)
          residue_N_in_day = residue_N_in_day + nin
          N_org_age(i_org) = 1
        END BLOCK
      END IF
    END DO
    IF (soilN_layered == 1) THEN
      ! Layered bridging for the scalar residue module, which is itself untouched:
      ! the module reads and writes the scalar mirrors, and its net ammonium and
      ! nitrate changes are routed into the layer states here - additions to layer 1,
      ! whose thickness is deck-specific (2.0-5.0 cm from the soil-water driver), and
      ! removals (immobilization) in proportion to the layer pools. This bridging
      ! rule is a model construction.
      BLOCK
        REAL(wp) :: nh4_b, no3_b, d_nh4, d_no3, ssum
        nh4_b = N_NH4_soil; no3_b = N_soil_min
        CALL decompose_residue_day(T_soil_day, MIN(1.0_wp, vwc_root / MAX(1.0e-6_wp, theta_ref_dec)), &
          decomp_C_day, decomp_N_day, respired_C_day, net_min_NH4_day, immob_NH4_day, immob_NO3_day, turnover_N_day)
        d_nh4 = N_NH4_soil - nh4_b; d_no3 = N_soil_min - no3_b
        IF (d_nh4 >= 0.0_wp) THEN
          NH4_lay(1) = NH4_lay(1) + d_nh4
        ELSE
          ssum = SUM(NH4_lay)
          IF (ssum > 1.0e-12_wp) NH4_lay = NH4_lay * (1.0_wp + d_nh4 / ssum)
        END IF
        IF (d_no3 >= 0.0_wp) THEN
          NO3_lay(1) = NO3_lay(1) + d_no3
        ELSE
          ssum = SUM(NO3_lay)
          IF (ssum > 1.0e-12_wp) NO3_lay = NO3_lay * (1.0_wp + d_no3 / ssum)
        END IF
      END BLOCK
    ELSE
    CALL decompose_residue_day(T_soil_day, MIN(1.0_wp, vwc_root / MAX(1.0e-6_wp, theta_ref_dec)), &
      decomp_C_day, decomp_N_day, respired_C_day, net_min_NH4_day, immob_NH4_day, immob_NO3_day, turnover_N_day)
    END IF
  END IF

  ! (3) Urea hydrolysis (urea to ammonium) and (4) nitrification (ammonium to
  !     nitrate), both first-order and temperature-dependent,
  !     k(T) = k_ref Q10^((T_soil - T_ref)/10) d-1. The hydrolysis constants are
  !     0.07, 0.12 and 0.20 d-1 at 5, 15 and 25 C with a Q10 of about 2 (Lasisi &
  !     Akinremi 2020 Can J Soil Sci, doi:10.1139/cjss-2020-0044, their Table 3);
  !     the nitrification reference 0.15 d-1 at 25 C is from Taylor (2017) and
  !     Stanford & Smith (1972).
  BLOCK
    REAL(wp) :: k_hyd, k_nit
    k_hyd = k_hyd_ref * Q10_hyd ** ((T_soil_day - T_ref_hyd) / 10.0_wp)
    k_nit = k_nit_ref * Q10_nit ** ((T_soil_day - T_ref_nit) / 10.0_wp)
    N_hydrolysis_day    = MIN(N_urea_soil, MAX(0.0_wp, k_hyd) * N_urea_soil)
    N_urea_soil = N_urea_soil - N_hydrolysis_day
    ! Hydrolysis of surface-applied urea feeds the SURFACE ammoniacal pool, the
    ! emission layer, and not the bulk pool. In layered mode the routing is
    ! identical, urea being layer-1 only by design, which preserves the
    ! volatilization semantics.
    TAN_surf    = TAN_surf    + N_hydrolysis_day
    ! Nitrification acts on both compartments at the same first-order,
    ! temperature-dependent rate (same microbial process).
    IF (soilN_layered == 1) THEN
      ! Layered: PER-LAYER nitrification with the existing law and constants
      ! unchanged, driven by a single daily soil temperature, since no soil-
      ! temperature profile state exists on this path. The nitrified surface
      ! ammoniacal nitrogen goes to the layer-1 nitrate, a surface process. Ammonium
      ! is NOT advected, being held by cation exchange on the sandy-soil archetype;
      ! that is a model assumption.
      BLOCK
        REAL(wp) :: nit_surf, dnit_l
        INTEGER  :: iLn
        nit_surf = MIN(TAN_surf, MAX(0.0_wp, k_nit) * TAN_surf)
        TAN_surf = TAN_surf - nit_surf
        NO3_lay(1) = NO3_lay(1) + nit_surf
        N_nitrification_day = nit_surf
        DO iLn = 1, NSoil_layers
          dnit_l = MIN(NH4_lay(iLn), MAX(0.0_wp, k_nit) * NH4_lay(iLn))
          NH4_lay(iLn) = NH4_lay(iLn) - dnit_l
          NO3_lay(iLn) = NO3_lay(iLn) + dnit_l
          N_nitrification_day = N_nitrification_day + dnit_l
        END DO
        ! Aggregate mirrors (exact sums; conserve to round-off)
        N_NH4_soil = SUM(NH4_lay)
        N_soil_min = SUM(NO3_lay)
      END BLOCK
    ELSE
    BLOCK
      REAL(wp) :: nit_surf
      nit_surf = MIN(TAN_surf, MAX(0.0_wp, k_nit) * TAN_surf)
      TAN_surf = TAN_surf - nit_surf
      N_nitrification_day = MIN(N_NH4_soil, MAX(0.0_wp, k_nit) * N_NH4_soil) + nit_surf
      N_NH4_soil  = N_NH4_soil  - (N_nitrification_day - nit_surf)
      N_soil_min  = N_soil_min  + N_nitrification_day    ! NO3 pool now set for uptake/leaching competition
    END BLOCK
    END IF
  END BLOCK

  ! -------------------------------------------------------------------- !
  ! Mechanistic ammonia volatilization from the surface ammoniacal pool: a chain
  ! of measured chemistry and the model's own resistance network, with no fitted
  ! constants.
  !
  ! (1) Dissolved concentration in the emission layer, with the ammonium
  !     solid / solution partition on the exchange complex:
  !       C_aq [g N m-3 water] = (TAN / z_surf) / (theta_1 + BD Kd)
  !     the linear-partition form, with the distribution coefficient measured at
  !     1-4 L kg-1 (Venterea et al. 2015).
  ! (2) The un-ionized fraction from the measured dissociation equilibrium:
  !       f = 1 / (1 + 10^(pKa - pH)),  pKa(T) = 0.09018 + 2729.92/T_K
  !     (Emerson et al. 1975 J Fish Res Board Can 32:2379, Table 1).
  ! (3) The gas-liquid partition from the measured Henry constant:
  !       KH(T) = 10^(1477.7/T_K - 1.6937)  (dimensionless, aqueous over gas)
  !     (Hales & Drewes 1979, as compiled in Genermont & Cellier 1997 Agric For
  !     Meteorol 88:145, Eq. A6).
  ! (4) The flux through the model's own Shuttleworth-Wallace resistances, the
  !     soil-to-canopy-air and canopy-air-to-reference resistances:
  !       F = chi_s / (r_sa + r_ca) x 86400
  !     (the resistance analogy of Genermont & Cellier 1997).
  ! (5) Infiltration advects the surface ammoniacal nitrogen into the bulk pool
  !     by piston flow, the same Burns (1975) form as the leaching, which
  !     reproduces the measured result that 10 mm or more of rain stops
  !     volatilization with no new parameter (Sommer et al. 2004, sec. IV.C).
  ! -------------------------------------------------------------------- !
  N_volat_day = 0.0_wp
  BLOCK
    REAL(wp) :: T_K, pKa, f_nh3, KH, C_aq, chi_s, r_path, F_vol, F_down, w_surf
    IF (TAN_surf > 1.0e-9_wp) THEN
      T_K   = T_soil_day + 273.15_wp
      pKa   = 0.09018_wp + 2729.92_wp / T_K
      f_nh3 = 1.0_wp / (1.0_wp + 10.0_wp**(pKa - pH_surf))
      KH    = 10.0_wp**(1477.7_wp / T_K - 1.6937_wp)
      C_aq  = (TAN_surf / z_surf_N) / MAX(1.0e-3_wp, theta(1) + BD_surf * Kd_NH4)
      chi_s = C_aq * f_nh3 / MAX(1.0e-6_wp, KH)          ! g N m-3 air at surface
      r_path = MAX(10.0_wp, rsa_day + rca_day)           ! s m-1 (model's own SW network)
      F_vol  = chi_s / r_path * 86400.0_wp               ! g N m-2 d-1
      N_volat_day = MIN(TAN_surf, MAX(0.0_wp, F_vol))
      TAN_surf    = TAN_surf - N_volat_day
      ! (5) infiltration advection surface -> bulk (piston flow)
      w_surf = MAX(1.0_wp, theta(1) * z_surf_N * 1000.0_wp)   ! mm water in layer
      F_down = TAN_surf * MIN(1.0_wp, MAX(0.0_wp, infil_mm_d) / w_surf)
      TAN_surf   = TAN_surf   - F_down
      IF (soilN_layered == 1) THEN
        ! Layered: TAN incorporation enters LAYER-1 NH4 (design A — the
        ! emission layer 0-z_surf_N sits inside Richards layer 1).
        NH4_lay(1) = NH4_lay(1) + F_down
        N_NH4_soil = N_NH4_soil + F_down   ! aggregate mirror
      ELSE
      N_NH4_soil = N_NH4_soil + F_down
      END IF
    END IF
  END BLOCK
  N_volat_total = N_volat_total + N_volat_day

  ! -------------------------------------------------------------------- !
  ! Mechanistic denitrification:
  !   D = 0.933 respired_C f_W(WFPS) NO3/(NO3 + Km)
  ! - 0.933 g N per g C is the ELECTRON-BALANCE stoichiometry of complete
  !   heterotrophic denitrification, 5 CH2O + 4 NO3- + 4 H+ -> 5 CO2 + 2 N2 +
  !   7 H2O, i.e. 4 x 14 / (5 x 12) g N per g C oxidized (Firestone 1982, in
  !   Nitrogen in Agricultural Soils, Agron Monogr 22, the equation on p.289), so
  !   the electrons come from the same decomposition carbon flux the model
  !   already simulates, carrying its own temperature and moisture response with
  !   no double counting.
  ! - f_W is the anaerobic activity fraction against the water-filled pore
  !   space: zero below 0.6, rising to 1 at saturation, fitted to the MEASURED
  !   denitrifying-activity points of Linn & Doran (1984) SSSAJ 48:1267, Fig. 4,
  !   on the anaerobic-microsite basis of Smith (1980) J Soil Sci 31:263,
  !   Eqs. 1-4.
  ! - The Michaelis nitrate term uses a half-saturation constant measured at
  !   22-53 mg N L-1 (Murray et al. 1989), with the nitrate solution
  !   concentration from the same root-zone water the uptake and the leaching
  !   use, and the water-filled pore space from the model's own Richards state.
  ! -------------------------------------------------------------------- !
  N_denit_day = 0.0_wp
  BLOCK
    REAL(wp) :: wfps, f_w, no3_L, f_n, D_stoich
    IF (N_soil_min > 1.0e-9_wp .AND. respired_C_day > 0.0_wp) THEN
      ! The water-filled pore space is taken where the ELECTRON DONOR is: heterotrophic
      ! carbon oxidation is concentrated in the residue- and organic-matter-rich
      ! topsoil, so the anaerobic exposure that matters is that of the top 0.3 m,
      ! denitrifier activity being co-located with decomposable carbon (Burford &
      ! Bremner 1975; Parkin 1987 SSSAJ 51:1194, hot spots at buried residue). The
      ! root-zone value over-states the anaerobiosis early in the season and the
      ! whole-column value erases the post-irrigation topsoil pulses.
      BLOCK
        REAL(wp) :: sw, ss, zt, ov
        INTEGER  :: iL
        sw = 0.0_wp; ss = 0.0_wp; zt = 0.0_wp
        DO iL = 1, NSoil_layers
          ov = MAX(0.0_wp, MIN(zt + tk(iL), 0.30_wp) - zt)
          sw = sw + ov * theta(iL)
          ss = ss + ov * theta_0(iL)
          zt = zt + tk(iL)
          IF (zt >= 0.30_wp) EXIT
        END DO
        wfps = MIN(1.0_wp, sw / MAX(1.0e-6_wp, ss))
      END BLOCK
      IF (wfps > 0.6_wp) THEN
        ! A quartic fit to the MEASURED denitrifying-activity points of Linn & Doran
        ! (1984) SSSAJ 48:1267, Fig. 4: near zero at a water-filled pore space of 0.75
        ! or below (about 0.004 at 0.70, 0.06 at 0.80, 0.32 at 0.90, 1.0 at
        ! saturation). Their curve is much steeper than quadratic, and a quadratic
        ! over-predicts seasonal denitrification several-fold above the measured field
        ! ranges of Aulakh et al. (1992) Adv Soil Sci 18:1, so the exponent follows the
        ! SOURCE data.
        f_w   = ((wfps - 0.6_wp) / 0.4_wp)**4
        IF (soilN_layered == 1) THEN
          ! Layered, with the law untouched: the nitrate concentration argument becomes
          ! the ROOT-ZONE AGGREGATE of the layered nitrate, the layer pools weighted by
          ! their overlap with the root depth over the root-zone water, and the removal
          ! is taken from the root-zone layers in proportion to their contribution.
          BLOCK
            REAL(wp) :: zt2, ov2, no3_rz_sum, wrz(NSoil_layers)
            INTEGER  :: iLd
            zt2 = 0.0_wp; no3_rz_sum = 0.0_wp; wrz = 0.0_wp
            DO iLd = 1, NSoil_layers
              ov2 = MAX(0.0_wp, MIN(zt2 + tk(iLd), MAX(0.1_wp, root_dep)) - zt2)
              IF (tk(iLd) > 1.0e-9_wp) wrz(iLd) = (ov2 / tk(iLd)) * NO3_lay(iLd)
              no3_rz_sum = no3_rz_sum + wrz(iLd)
              zt2 = zt2 + tk(iLd)
            END DO
            no3_L = no3_rz_sum / MAX(1.0_wp, vwc_root * MAX(0.1_wp, root_dep) * 1000.0_wp)  ! g N L-1
            f_n   = no3_L / (no3_L + Km_den)
            D_stoich = 0.933_wp * respired_C_day * f_w * f_n
            N_denit_day = MIN(no3_rz_sum, MAX(0.0_wp, D_stoich))
            IF (N_denit_day > 0.0_wp .AND. no3_rz_sum > 1.0e-12_wp) THEN
              DO iLd = 1, NSoil_layers
                NO3_lay(iLd) = NO3_lay(iLd) - N_denit_day * wrz(iLd) / no3_rz_sum
              END DO
              N_soil_min = SUM(NO3_lay)   ! aggregate mirror
            END IF
          END BLOCK
        ELSE
        no3_L = N_soil_min / MAX(1.0_wp, vwc_root * MAX(0.1_wp, root_dep) * 1000.0_wp)  ! g N L-1
        f_n   = no3_L / (no3_L + Km_den)
        D_stoich = 0.933_wp * respired_C_day * f_w * f_n
        N_denit_day = MIN(N_soil_min, MAX(0.0_wp, D_stoich))
        N_soil_min  = N_soil_min - N_denit_day
        END IF
      END IF
    END IF
  END BLOCK
  N_denit_total = N_denit_total + N_denit_day

  ! -------------------------------------------------------------------- !
  ! Plant nitrogen uptake, demand-limited Michaelis-Menten.
  !
  ! The demand comes from the critical dilution curve
  !   Nc(%) = 5.35 W^(-0.442), W in t ha-1 of shoot dry matter,
  ! Greenwood et al. (1990) Ann Bot 66:585-597, Eq. 3 and Table 1, and for winter
  ! wheat Justes et al. (1994) Ann Bot 74:397-407, p.402, Table 2.
  !
  ! The uptake kinetics are root-scaled active uptake,
  !   S_active = u_max_dw f_Nstatus a_root W_root NO3_rz/(Km + NO3_rz),
  ! which replaces a flat whole-crop Michaelis-Menten law. The constant that law
  ! used descends from the abstract of Barraclough (1986) J Agric Sci 106:45-52,
  ! 2.75 kg N ha-1 d-1 of whole-crop realized uptake per unit ground area, which
  ! is not a root kinetic; both of its parameters are retired and their driver
  ! rows are kept positionally.
  !
  ! The transpiration scaling of nitrogen uptake is the ratio of the actual
  ! transpiration to its reference, after Lemaire & Gastal (2009) Ann Bot
  ! 106:35-49, p.38, Eq. 2: nitrogen uptake is proportional to the water flux
  ! through the root system by mass flow when nitrogen is not limiting, scaling
  ! with the actual transpiration.
  ! -------------------------------------------------------------------- !
  N_soil_uptake_day = 0.0_wp
  N_leach_d = 0.0_wp
  IF (soilN_layered == 1) THEN
  ! ==================================================================== !
  ! The LAYERED soil-nitrogen branch; the single-pool path is the branch below.
  ! The constants live in input/driver_N.csv, one source and DOI per row:
  !   the ammonium transporter maximum rate 3.5 umol g-1 fresh weight h-1 and
  !   half-saturation constant 0.05 mol m-3 (Goyal & Huffaker 1986 Plant Physiol
  !   82:1051-1056, doi:10.1104/pp.82.4.1051, p.1053);
  !   the mineralization depth 0.20 m, on the measured topsoil organic-matter
  !   basis;
  !   the advection sub-step guard 0.5, pure numerics; and the layered mode flag.
  ! The per-layer root shares reuse the EXISTING beta root-mass profile of
  ! Jackson et al. (1996) from hydro.f90, read-only.
  !
  ! (i) NITRATE ADVECTION: explicit upwind convection on the day SIGNED
  !     inter-layer water fluxes from the Richards accumulators, read-only; this
  !     block writes no water state. The interface nitrogen flux is the water
  !     flux times the UPSTREAM layer solution concentration, and the
  !     bottom-interface export is the leaching. A Courant guard sub-steps the
  !     day so that no layer advects out more than the guard fraction of its
  !     nitrate per sub-step. No dispersive or diffusive term is included. An
  !     upward bottom flux imports nothing, since no nitrogen state exists below
  !     the domain.
  ! (ii) AMMONIUM IS NOT ADVECTED: it is held by cation exchange on the
  !     sandy-soil archetype, a model assumption.
  ! (iii) UPTAKE PER LAYER, the nitrate structure kept and an ammonium system
  !     added:
  !       C_NO3 = u_max_dw     f_Nstatus a_root W_root
  !       C_NH4 = u_max_NH4_dw f_Nstatus a_root W_root
  !     with u_max_NH4_dw = 3.5 umol g-1 FW h-1 x 18.6 g FW per g DW x 14e-6 x 24
  !     = 2.188e-2 g N g-1 DW d-1 (Goyal & Huffaker 1986: a half-saturation
  !     constant of 0.05 mM and a maximum rate of 3.5 umol per g fresh weight per
  !     hour were calculated for the ammonium transporter, p.1053; the same
  !     intact-wheat-seedling system as their nitrate numbers, 2.6 and 0.027 mM,
  !     pp.1052 and 1055). The saturation their Fig. 3d shows is inside the
  !     Michaelis-Menten form, so no extra term is needed. The nitrogen-status
  !     function and the root activity, both derived on nitrate, are applied to
  !     the ammonium system unchanged.
  !     The per-layer shares are the existing Jackson beta root-mass profile,
  !     reused read-only. The active term per layer and ion is the capacity times
  !     the root share times the Michaelis-Menten factor in the layer solution
  !     water, and the mass-flow term is the layer's own root water uptake times
  !     its total solution concentration. The potential uptake is the lesser of
  !     the daily demand and the summed supply; the realized uptake scales every
  !     term by the common ratio, and the removal per layer and ion is
  !     proportional to its term. The solution volumes use the end-of-day water
  !     content with the daily integrated fluxes, an operator-splitting
  !     approximation.
  ! ==================================================================== !
  BLOCK
    REAL(wp) :: w_lay(NSoil_layers), fr_lay(NSoil_layers)
    REAL(wp) :: act_no3(NSoil_layers), act_nh4(NSoil_layers)
    REAL(wp) :: mf_no3(NSoil_layers), mf_nh4(NSoil_layers)
    REAL(wp) :: p_no3(NSoil_layers), p_nh4(NSoil_layers)
    REAL(wp) :: diff_no3(NSoil_layers), diff_nh4(NSoil_layers)   ! rec 139/: per-layer diffusion-term supply (g N m-2 d-1), 0.0 when OFF
    REAL(wp) :: u_max_dw, u_max_nh4_dw, C_NO3_cap, C_NH4_cap
    REAL(wp) :: f_Nstatus, a_root, c_no3, c_nh4, mf_tot, denomN, tot_ion
    REAL(wp) :: w_shoot_tha, n_shoot_up, nc_crit_pct, nni_up
    REAL(wp) :: U_pot_sum, U_pot, ratio, frac_out, frac_max, F_if, c_up, dt_adv
    REAL(wp) :: fl_theta, De_lay, RLD_lay   ! rec 139: per-layer diffusion-law working scalars
    REAL(wp) :: act_no3_eff(NSoil_layers), act_nh4_eff(NSoil_layers)  ! rec 139 mode 2: active term, EFFECTIVE for the cap + the tot_ion split (= act_ion unless mode 2 bounds it to mf_ion+diff_ion); PER LAYER (consumed by the second, ratio-scaling DO loop below)
    INTEGER  :: k_nni, iLu, n_adv, is_adv
    ! ---- (a) the daily crop nitrogen demand. The expression is written out here
    !      rather than shared, so that the single-pool branch below keeps its own
    !      behavior; here the demand is computed whenever the crop is present,
    !      because ammonium may supply it.
    N_demand_day = 0.0_wp
    IF (DVS > 0) THEN
      N_plant_total = N_leaf + N_stem + N_root + N_storage + N_dead_leaf
      ! The dead stem, sheath and chaff dry matter stays ABOVE GROUND - it is
      ! harvested straw and the observed above-ground biomass includes it - so it
      ! enters every above-ground dry-matter total. It is identically zero when
      ! record 63 is off.
      AGB_kg_m2 = MAX(0.01_wp, (W_green + W_dead + W_dead_stem + W_stem + &
                MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) / 100.0_wp) ! now t/ha
      ! The maximum-dilution demand cap, the same expression as the single-pool
      ! branch below: Justes et al. (1994) Ann Bot 74:397-407,
      ! doi:10.1006/anbo.1994.1133, maximum dilution curve Nmax = 8.3 W^-0.44, their
      ! Eq. (17) p.405 and Fig. 5 p.402, printed valid from 1.5 to 14 t ha-1 (p.404
      ! col. 2), with W in t ha-1. The 6.6 % value read below 1.55 t ha-1 is not in
      ! that paper; see record 160.
      ! ================================================================== !
      ! Record 160: WHAT THE OUTER MAXIMUM CAP IS WORTH BELOW THE VALIDITY FLOOR.
      ! Mode 0 is the 6.6 % plateau on its own line below. Mode 1 holds the PRINTED
      ! curve at its own value at the floor, 8.3 x 1.55**(-0.44) = 6.844 %DM, the
      ! same construction the floor in growth_module.f90 and record 155 limb (a)
      ! already use, and the same one the paper applies to the CRITICAL curve
      ! (p.402 col. 2). Zero new constants; the level is a continuity construction,
      ! because the paper prints no statement below 1.5 t ha-1.
      ! ================================================================== !
      IF (AGB_kg_m2 <= 1.55_wp) THEN
        IF (justes_nmax_floor_basis >= 1) THEN
          Nc_pct = 8.3_wp * 1.55_wp**(-0.44_wp)
          jnf_days_floor = jnf_days_floor + 1
        ELSE
        Nc_pct = 6.6_wp
        END IF
      ELSE
        Nc_pct = 8.3_wp * AGB_kg_m2**(-0.44_wp)
      END IF
      BLOCK
        REAL(wp) :: g_nmax_shoot, n_shoot, d_growth, d_reservoir, d_root, d_leaf_maint
        ! N buried in dead stem/chaff is still shoot N.
        n_shoot = N_leaf + N_dead_leaf + N_dead_stem + N_stem + N_storage
        ! same above-ground basis as AGB_kg_m2 above.
        g_nmax_shoot = MAX(0.0_wp, (Nc_pct / 100.0_wp) * (W_green + W_dead + W_dead_stem + W_stem + &
                  MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) - n_shoot)
        IF (dvs2_growth < 1.0_wp) THEN
          ! The leaf share of the prospective demand is bounded by the remaining
          ! leaf-nitrogen holding capacity, as coded:
          !   d_leaf = MIN(n4_sln_target dLAI_prev, MAX(0, sln_max_N LAI - N_leaf))
          ! with the construction target 1.5 g N m-2 (Jamieson & Semenov 2000 Field Crops
          ! Res 68:21-29 p.23) and the holding capacity 3.8 g N m-2 (Silva-Perez et al.
          ! 2020 J Exp Bot 71(7):2299-2311 Fig. 4B p.2305, the maximum measured nitrogen
          ! per unit leaf area); the full sources are at their SHARED.f90 declarations.
          IF (leafN_drawdown == 1) THEN
            d_growth = MIN(n4_sln_target * n5_dlai_prev, &
                           MAX(0.0_wp, sln_max_N * lai_growth_diag - N_leaf)) &
                     + n5_stem_nstruct * (n5_dwstem_prev + n5_dwchaff_prev)
          ELSE
            d_growth = n4_sln_target * n5_dlai_prev &
                     + n5_stem_nstruct * (n5_dwstem_prev + n5_dwchaff_prev)
          END IF
          ! The luxury-reservoir DEMAND is read from the same two organ pools the
          ! allocation fills, or the chaff share of the reservoir would be demanded twice
          ! - once as the gap in the stem pool against a stem-plus-chaff basis and again
          ! as the gap in the chaff pool - and the post-anthesis uptake would be
          ! inflated. As coded:
          !   d_reservoir = MAX(0, n4_stem_nmax W_stem  - N_stem)
          !               + MAX(0, n4_stem_nmax W_chaff - N_chaff)
          ! with the stem luxury ceiling 0.015 g N g-1 DM (Jamieson & Semenov 2000 Field
          ! Crops Res 68:21-29 p.23, a stem concentration of up to 1.5 % N to store
          ! excess nitrogen; the locus is at its SHARED.f90 declaration). The zeros are
          ! positivity guards.
          ! Under record 63 the chaff DEMAND is read against the same ceiling the chaff
          ! ALLOCATION fills, or the reservoir would be under-demanded: the chaff ceiling
          ! is 0.017 g N g-1 DM, derived from Bertheloot, Martre & Andrieu (2008) Plant
          ! Physiol Fig. 1D and Table I at its SHARED.f90 declaration, and not from any
          ! observation of this campaign. With record 63 off the inherited stem ceiling
          ! is selected.
          d_reservoir = MAX(0.0_wp, n4_stem_nmax * W_stem  - N_stem) &
                      + MAX(0.0_wp, MERGE(n2_chaff_nmax, n4_stem_nmax, grainN_partition == 1) &
                                    * W_chaff - N_chaff)
          d_root = phen_params%init_rootN_frac * n5_dwroot_prev
          veg_N_demand_day = d_growth + d_reservoir + d_root
          ! ================================================================ !
          ! Record 176 at or above 2: THE PRE-ANTHESIS DEMAND GAINS ITS STOCK-DEFICIT
          ! TERM. The vegetative demand above is a pure FLUX - a construction price on
          ! yesterday's dry-matter increment plus a luxury-reservoir gap - and it is
          ! compared with a ceiling that is a STOCK gap: a stock gap cannot bind a daily
          ! flux (zero binding days on 42 of 42 decks) and a pure flux can never recover
          ! arrears once the crop has fallen below the critical curve. The mirror defect
          ! on the post-anthesis side is named at the record-97 block: the post-anthesis
          ! uptake demand is a cumulative STOCK deficit while the grain sink is a RATE.
          !   Nc% = 4.4 for W at or below 1.55 t ha-1, else 5.35 W^-0.442 - Justes et al.
          ! (1994) Ann Bot 74:397-407, doi:10.1006/anbo.1994.1133, Eq. (11) p.401,
          ! printed valid from 1.55 to 12 t ha-1, with the 4.4 % plateau and the
          ! 1.55 t ha-1 break at p.402 col. 2, verified at first hand. The same two
          ! numbers are already coded at the nitrogen-nutrition-index diagnostic below,
          ! and the shoot dry-matter and shoot nitrogen bases are the IDENTICAL
          ! expressions the stock cap uses four lines above, so the deficit and its cap
          ! are on ONE basis by construction. The MAXIMUM curve (Eq. (17) p.405) stays as
          ! the outer cap, unchanged. The term is SELF-EXTINGUISHING and is a DEMAND
          ! only: it is still served through the unchanged single minimum against the
          ! layer supply, so a pool-limited low-nitrogen deck refuses it. Zero new
          ! constants.
          ! ================================================================ !
          grain_N_demand_postanth_day = 0.0_wp
          N_demand_day = MIN(veg_N_demand_day, g_nmax_shoot)
          ! ================================================================ !
          ! Record 155 LIMB (a): the PRE-ANTHESIS LUXURY DEMAND ON THE MAXIMUM DILUTION
          ! CURVE. The outer stock cap above already carries the MAXIMUM curve
          ! (8.3 W^-0.44, plateau 6.6 %) and it binds on ZERO pre-anthesis days on 42 of
          ! 42 released decks, so swapping the curve in the cap is a provable no-op. What
          ! binds is the vegetative demand itself, the sum of organ-MINIMUM construction
          ! targets, which cannot reach the nitrogen nutrition index of 1.3 the
          ! Braunschweig crop measured. Limb (a) therefore puts the SAME curve on the
          ! demand DRIVER, in the marginal form this code already uses for the CRITICAL
          ! curve in growth_module.f90, c_marg = a(1-b) W^-b with W floored at the
          ! validity floor 1.55 t ha-1 (Pask 2009 PhD thesis, University of Nottingham,
          ! p.151). That marginal form IS the STICS maximum-curve DEMAND, held at first
          ! hand: Brisson et al. (1998) Agronomie 18:311-346,
          ! doi:10.1051/agro:19980501, p.342 Eq. (47), DEMANDE = 10 ADILMAX DLTAMS
          ! (1 - BDILMAX) MASEC^-BDILMAX, and p.329, the product of the daily growth rate
          ! and the derivative of the maximum curve of nitrogen dilution. Its Table VIII
          ! p.330 prints a coefficient of 8.47 for wheat where this model codes 8.3; the
          ! model's OWN 8.3 is used here, so the record adds no coefficient. Justes et al.
          ! (1994) was read at first hand and prints Eq. (17) p.405,
          ! Nmax = 8.3 DM^-0.44, so 8.3 IS the value of that paper and 8.47 is the STICS
          ! own number. That paper prints NO maximum value below 1.5 t ha-1, so the 6.6 %
          ! plateau has no source; record 160 mode 1 replaces it with the curve's own floor
          ! value, and limb (a) never reads it either way. Also at first hand, p.405: the
          ! maximum rate of nitrogen accumulation is k times the maximum growth rate with
          ! k = 46.5 DM^-0.44, which is the primary source for this marginal law, of
          ! which Eq. (47) of Brisson et al. (1998) is the STICS form.
          ! The crop therefore demands at least enough nitrogen to build TODAY'S OWN new
          ! above-ground dry matter at the MAXIMUM dilution concentration, and is still
          ! never allowed past the maximum-curve stock gap. The dry-matter increments are
          ! the existing one-day-lagged realized increments the demand law already reads.
          ! No new constant: 8.3, 0.44, 1.55 and 0.01 are all already in this file or in
          ! growth_module.f90.
          ! ================================================================ !
          IF (n_uptake_regulation == 1) THEN
            BLOCK
              REAL(wp) :: nur_W_eff, nur_c_max_marg, nur_dem_inc
              nur_W_eff      = MAX(AGB_kg_m2, 1.55_wp)
              nur_c_max_marg = (8.3_wp * (1.0_wp - 0.44_wp) * nur_W_eff**(-0.44_wp)) * 0.01_wp
              nur_D_max_day  = nur_c_max_marg * (nur_dwgreen_prev + n5_dwstem_prev + n5_dwchaff_prev) &
                             + d_root
              nur_dem_inc    = MIN(MAX(veg_N_demand_day, nur_D_max_day), g_nmax_shoot)
              IF (nur_dem_inc > N_demand_day + 1.0e-12_wp) THEN
                nur_days_a_bind = nur_days_a_bind + 1
                nur_cum_a_gain  = nur_cum_a_gain + (nur_dem_inc - N_demand_day)
              END IF
              N_demand_day = nur_dem_inc
            END BLOCK
          END IF
        ELSE
          ! The maintenance demand is bounded identically: the target is the lesser of
          ! the construction target and the holding capacity, which keeps the demand
          ! within the leaf capacity. At the released values (1.5 below 3.8) it is a
          ! structural guard only.
          IF (leafN_drawdown == 1) THEN
            d_leaf_maint = MAX(0.0_wp, MIN(n4_sln_target, sln_max_N) * lai_growth_diag - N_leaf)
          ELSE
            d_leaf_maint = MAX(0.0_wp, n4_sln_target * lai_growth_diag - N_leaf)
          END IF
          veg_N_demand_day = d_leaf_maint
          grain_N_demand_postanth_day = n4_grain_deficit_prev
          N_demand_day = MIN(veg_N_demand_day + grain_N_demand_postanth_day, g_nmax_shoot)
        END IF
        ! ================================================================ !
        ! Record 155 LIMB (b): the DOWN-REGULATION OF THE WHOLE-PLANT UPTAKE DEMAND BY
        ! THE CULM MOBILE-NITROGEN CONCENTRATION, with its two ceiling diagnostics.
        !   [N]mob = N_mobile_veg / (W_green + W_stem)        g g-1
        !   E_N    = exp(-beta_N [N]mob),  beta_N = 130
        !   N_demand_day <- N_demand_day E_N
        ! Eq. (4) of Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096 with
        ! the coefficient of Bertheloot, Wu, Cournede & Andrieu (2011) Table 2. It is
        ! applied on EVERY live day because that model applies it on every day;
        ! pre-anthesis the pool is near zero in this model, so the factor is near 1
        ! there and the realized value is REPORTED per day rather than asserted. The
        ! denominator is the lamina-plus-stem basis the record-154 nitrogen factor
        ! uses. The critical-curve stock gap is computed beside it as a DIAGNOSTIC
        ! read by no statement, so that both ceilings are visible side by side.
        ! ================================================================ !
        nur_dem_before   = N_demand_day
        nur_cap_nmax_day = g_nmax_shoot
        ! Record 160: on a sub-floor day (shoot dry matter at or below 1.55 t ha-1),
        ! did the raised maximum-curve STOCK GAP bind today's demand? Diagnostic only: no
        ! statement reads the counter and it is printed in the season summary.
        IF (justes_nmax_floor_basis >= 1 .AND. AGB_kg_m2 <= 1.55_wp) THEN
          IF (N_demand_day >= g_nmax_shoot - 1.0e-12_wp) jnf_days_bind = jnf_days_bind + 1
        END IF
        BLOCK
          REAL(wp) :: nur_nc_pct
          IF (AGB_kg_m2 <= 1.55_wp) THEN
            nur_nc_pct = 4.4_wp
          ELSE
            nur_nc_pct = 5.35_wp * AGB_kg_m2**(-0.442_wp)
          END IF
          nur_cap_nc_day = MAX(0.0_wp, (nur_nc_pct / 100.0_wp) * (W_green + W_dead + W_dead_stem + W_stem + &
                    MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) - n_shoot)
        END BLOCK
        ! Record 157 mode 3: WHAT ENTERS the limb-(b) mobile-nitrogen NUMERATOR. At
        ! modes 0, 1 and 2 the numerator is the positive part of the mobile vegetative
        ! pool. At mode 3 it also carries the shoot nitrate store and the stem and
        ! chaff nitrogen above the structural floor of Jamieson & Semenov (2000) p.23,
        ! both of them circulating amino acids and nitrate in the definition of
        ! Bertheloot, Cournede & Andrieu (2011). The record-157 LAMINA STORE is NOT
        ! counted: that paper calls the remobilizable nitrogen of a photosynthetic
        ! entity only the nitrogen associated with the photosynthetic apparatus, which
        ! is what the lamina store is here. The DENOMINATOR is unchanged.
        lns_nmob_t_mobile  = MAX(0.0_wp, N_mobile_veg)
        lns_nmob_t_no3     = 0.0_wp ; lns_nmob_t_stemlux = 0.0_wp ; lns_nmob_t_chaflux = 0.0_wp
        IF (lamina_n_storage == 3) THEN
          lns_nmob_t_no3     = MAX(0.0_wp, N_NO3_shoot)
          lns_nmob_t_stemlux = MAX(0.0_wp, N_stem  - n5_stem_nstruct * W_stem)
          lns_nmob_t_chaflux = MAX(0.0_wp, N_chaff - n5_stem_nstruct * W_chaff)
        END IF
        lns_nmob_num_day = lns_nmob_t_mobile + lns_nmob_t_no3 + lns_nmob_t_stemlux + lns_nmob_t_chaflux
        ! Record 176 at or above 1: Bertheloot, Cournede & Andrieu (2011) p.1090 define
        ! the regulating quantity as the circulating amino acids and nitrate. Stem and
        ! chaff nitrogen held at the luxury ceiling is STORED and not circulating, and
        ! it carries 94-100 % of the pre-anthesis signal, which is why the regulation
        ! factor discards 58.3 % of the pre-anthesis demand and reads the SAME value at
        ! 40 and at 240 kg N ha-1. The numerator returns to the two circulating terms
        ! the model has. The coefficient 130 and the DENOMINATOR are unchanged, and a
        ! term is dropped rather than added.
        ! Record 159: the published ESTIMATE of the mobile pool REPLACES the numerator
        ! formed above. Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot 108(6):1097,
        ! doi:10.1093/aob/mcr125: remobilizable plus mobile nitrogen mass allocated to
        ! a module was estimated as the difference between total nitrogen mass and
        ! structural nitrogen mass, and mobile nitrogen mass was estimated as in
        ! Bertheloot et al. (2008a), assuming that it represented 10 % of the
        ! remobilizable nitrogen mass. Those two sentences solve exactly for
        ! N_rem = (N_tot - N_struct)/(1+r) and N_mob = r N_rem, which is already this
        ! model's own photosynthetic-nitrogen definition under record 154, so the three
        ! terms below ARE the record-154 remobilizable pools, organ for organ and floor
        ! for floor. The mobile fraction is a parameter alias of the existing 0.10, so
        ! no new number enters the build. Mode 1 is faithful, the estimate being an
        ! estimate OF that single mobile pool and so a replacement; mode 2 adds the
        ! explicit transport pools and is the double count. The DENOMINATOR below is
        ! unchanged.
        IF (n_uptake_regulation == 1) THEN
          nur_Nmob_day = lns_nmob_num_day / MAX(1.0e-9_wp, W_green + W_stem)
          nur_EN_day   = EXP(-n1_beta_N * nur_Nmob_day)
          N_demand_day = N_demand_day * nur_EN_day
          nur_days_live = nur_days_live + 1
          nur_EN_sum    = nur_EN_sum + nur_EN_day
          nur_EN_min    = MIN(nur_EN_min, nur_EN_day)
        END IF
        nur_dem_after = N_demand_day
      END BLOCK
    END IF
    ! ---- (b) NO3 advection: explicit upwind, CFL sub-stepped ---- !
    DO iLu = 1, NSoil_layers
      w_lay(iLu) = MAX(1.0e-6_wp, theta(iLu) * tk(iLu) * 1000.0_wp)   ! mm solution water
    END DO
    frac_max = 0.0_wp
    DO iLu = 1, NSoil_layers
      frac_out = 0.0_wp
      IF (q_iface_day_mm(iLu) > 0.0_wp) frac_out = frac_out + q_iface_day_mm(iLu) / w_lay(iLu)
      IF (iLu > 1) THEN
        IF (q_iface_day_mm(iLu-1) < 0.0_wp) frac_out = frac_out - q_iface_day_mm(iLu-1) / w_lay(iLu)
      END IF
      frac_max = MAX(frac_max, frac_out)
    END DO
    n_adv = 1
    IF (frac_max > cfl_frac) n_adv = MIN(2000, CEILING(frac_max / MAX(1.0e-6_wp, cfl_frac)))
    N_adv_substeps_day = REAL(n_adv, wp)
    dt_adv = 1.0_wp / REAL(n_adv, wp)
    DO is_adv = 1, n_adv
      DO iLu = 1, NSoil_layers - 1
        IF (q_iface_day_mm(iLu) >= 0.0_wp) THEN
          c_up = NO3_lay(iLu)   / w_lay(iLu)      ! g N per mm water, upstream above
        ELSE
          c_up = NO3_lay(iLu+1) / w_lay(iLu+1)    ! upstream below (upward flow)
        END IF
        F_if = q_iface_day_mm(iLu) * dt_adv * c_up
        NO3_lay(iLu)   = NO3_lay(iLu)   - F_if
        NO3_lay(iLu+1) = NO3_lay(iLu+1) + F_if
      END DO
      IF (q_iface_day_mm(NSoil_layers) > 0.0_wp) THEN
        F_if = q_iface_day_mm(NSoil_layers) * dt_adv * NO3_lay(NSoil_layers) / w_lay(NSoil_layers)
        F_if = MIN(F_if, NO3_lay(NSoil_layers))
        NO3_lay(NSoil_layers) = NO3_lay(NSoil_layers) - F_if
        N_leach_d = N_leach_d + F_if                 ! bottom export = leaching
      END IF
    END DO
    ! ---- (c) per-layer uptake ---- !
    U_act_NO3_day = 0.0_wp; U_act_NH4_day = 0.0_wp
    N_uptake_massflow_day = 0.0_wp; N_uptake_active_day = 0.0_wp
    N_uptake_diffusion_day = 0.0_wp   ! rec 139
    S_active_pot_day = 0.0_wp; S_massflow_pot_day = 0.0_wp; U_pot_day = 0.0_wp
    S_diffusion_pot_day = 0.0_wp      ! rec 139
    ! Unit chain, identical to the single-pool branch below: the transporter
    ! maximum rates in umol N g-1 fresh weight h-1 (Goyal & Huffaker 1986 Plant
    ! Physiol 82:1051-1056, doi:10.1104/pp.82.4.1051, nitrate 2.6 p.1055 and
    ! ammonium 3.5 p.1053) times the root fresh-to-dry ratio 18.6 (Khan et al. 2024
    ! Heliyon 10:e29042 Table 3) times 14.0e-6 g N umol-1 times 24 h d-1, giving
    ! g N g-1 DW d-1.
    u_max_dw     = u_max_fw     * fw_dw_root * 14.0e-6_wp * 24.0_wp
    u_max_nh4_dw = u_max_nh4_fw * fw_dw_root * 14.0e-6_wp * 24.0_wp   ! 2.188e-2 g N gDW-1 d-1 (G&H86 p.1053; same FW:DW as NO3)
    ! The nitrogen-status function of the nitrogen nutrition index, the same law as
    ! the single-pool branch below. The critical dilution curve is
    ! Nc% = 5.35 W^-0.442 for shoot dry matter above 1.55 t ha-1 and a plateau of
    ! 4.4 % below it (Justes et al. 1994 Ann Bot 74:397-407,
    ! doi:10.1006/anbo.1994.1133). The node mapping is from the four measured wheat
    ! points of Taulemesse et al. (2015) PLoS ONE 10(3):e0120291,
    ! doi:10.1371/journal.pone.0120291: the abscissae are the nitrogen nutrition
    ! indices at flowering of Table 1 p.7, and the ordinates come from the printed
    ! sentence p.8, that relative to the highest nitrogen treatment the daily
    ! nitrate uptake fluxes were 5.7, 5 and 2.4 times higher at the three lower
    ! treatments, normalized to the lowest. The above-ground basis is the same one
    ! the demand uses.
    w_shoot_tha = (W_green + W_dead + W_dead_stem + W_stem + &
              MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) / 100.0_wp
    n_shoot_up  = N_leaf + N_dead_leaf + N_stem + N_storage
    nni_up = 0.0_wp
    IF (w_shoot_tha > 1.0e-6_wp) THEN
      IF (w_shoot_tha <= 1.55_wp) THEN
        nc_crit_pct = 4.4_wp
      ELSE
        nc_crit_pct = 5.35_wp * w_shoot_tha**(-0.442_wp)
      END IF
      nni_up = (n_shoot_up / w_shoot_tha) / nc_crit_pct
    END IF
    IF (nni_up <= fNstat_nni(1)) THEN
      f_Nstatus = fNstat_f(1)
    ELSE IF (nni_up >= fNstat_nni(4)) THEN
      f_Nstatus = fNstat_f(4)
    ELSE
      f_Nstatus = fNstat_f(4)
      DO k_nni = 1, 3
        IF (nni_up <= fNstat_nni(k_nni+1)) THEN
          f_Nstatus = fNstat_f(k_nni) + (fNstat_f(k_nni+1) - fNstat_f(k_nni)) * &
                      (nni_up - fNstat_nni(k_nni)) / (fNstat_nni(k_nni+1) - fNstat_nni(k_nni))
          EXIT
        END IF
      END DO
    END IF
    ! a_root: the incumbent law, unchanged (Zheng 2022; DAA basis, quantity with no published source).
    a_root = 1.0_wp
    IF (daa_root > 0) a_root = MAX(aroot_floor, 1.0_wp - aroot_decline_d * REAL(daa_root, wp))
    f_Nstatus_diag = f_Nstatus
    a_root_diag    = a_root
    ! Record 134: the SAME critical curve restricted to the organs the validity
    ! statement of Pask (2009) describes, i.e. leaf and stem and not the
    ! reproductive pool, which is the vegetative-basis nitrogen nutrition index. It
    ! is computed unconditionally as read-only context and mirrored into SHARED for
    ! the cross-module diagnostic snapshot of growth_module.f90.
    BLOCK
      REAL(wp) :: w_veg_tha, n_veg_up, nc_crit_veg
      w_veg_tha = (W_green + W_dead + W_dead_stem + W_stem) / 100.0_wp
      n_veg_up  = N_leaf + N_dead_leaf + N_stem
      nap_nni_veg = 1.0_wp
      IF (w_veg_tha > 1.0e-6_wp) THEN
        IF (w_veg_tha <= 1.55_wp) THEN
          nc_crit_veg = 4.4_wp
        ELSE
          nc_crit_veg = 5.35_wp * w_veg_tha**(-0.442_wp)
        END IF
        nap_nni_veg = (n_veg_up / w_veg_tha) / nc_crit_veg
      END IF
    END BLOCK
    nap_f_Nstatus_day = f_Nstatus
    C_NO3_cap = u_max_dw     * f_Nstatus * a_root * W_root
    C_NH4_cap = u_max_nh4_dw * f_Nstatus * a_root * W_root
    fr_lay = hydro_root_dist_diag   ! existing Jackson profile (read-only reuse)
    diff_no3 = 0.0_wp; diff_nh4 = 0.0_wp   ! rec 139: reset, filled below ONLY when ON
    U_pot_sum = 0.0_wp
    DO iLu = 1, NSoil_layers
      c_no3 = NO3_lay(iLu) * 1000.0_wp / (w_lay(iLu) * 14.0_wp)   ! mol m-3 solution
      c_nh4 = NH4_lay(iLu) * 1000.0_wp / (w_lay(iLu) * 14.0_wp)
      act_no3(iLu) = C_NO3_cap * fr_lay(iLu) * c_no3 / (Km_c_no3 + c_no3)
      act_nh4(iLu) = C_NH4_cap * fr_lay(iLu) * c_nh4 / (Km_nh4   + c_nh4)
      mf_tot = rootup_layer_day_mm(iLu) * (NO3_lay(iLu) + NH4_lay(iLu)) / w_lay(iLu)
      denomN = NO3_lay(iLu) + NH4_lay(iLu)
      IF (denomN > 1.0e-12_wp) THEN
        mf_no3(iLu) = mf_tot * NO3_lay(iLu) / denomN
        mf_nh4(iLu) = mf_tot * NH4_lay(iLu) / denomN
      ELSE
        mf_no3(iLu) = 0.0_wp; mf_nh4(iLu) = 0.0_wp
      END IF
      ! Record 139: diffusion to the root surface, the third supply term. At mode 0
      ! this block does not execute and the diffusion terms stay at the zero reset
      ! above, so the extended cap below is an exact addition of zero.
      IF (n_uptake_diffusion >= 1 .AND. n_uptake_diffusion < 3) THEN
        ! f_l(theta): the SAME pair already coded at
        ! mineral_module.f90 -- reused verbatim, NOT a second pair (Barraclough & Tinker 1981 J Soil Sci 32:225-236 Fig. 1, printed 1.58*theta - 0.17).
        fl_theta = MAX(0.0_wp, 1.60_wp * theta(iLu) - 0.17_wp)
        ! De: Barber-Cushman effective diffusion coefficient in bulk soil
        ! (m2 s-1); nud_Dl_no3/nud_b_no3 stated, SHARED.f90.
        De_lay  = nud_Dl_no3 * theta(iLu) * fl_theta / nud_b_no3
        ! RLD: root length density (m root m-3 soil); nud_SRL stated,
        ! SHARED.f90.
        RLD_lay = W_root * fr_lay(iLu) * nud_SRL / MAX(tk(iLu), 1.0e-6_wp)
        ! Unit-chain closure: the literal product of root length density, the effective
        ! diffusion coefficient and the concentration difference is a volumetric soil
        ! sink (mol N m-3 soil s-1); multiplying by the layer thickness integrates it
        ! over the layer's own bulk thickness to a per-ground-area flux (mol m-2 s-1),
        ! and the 86400 s d-1 and 14.0 g N mol-1 are unit conversions.
        diff_no3(iLu) = MAX(0.0_wp, c_no3 - Km_c_no3) * RLD_lay * De_lay * tk(iLu) * 86400.0_wp * 14.0_wp
        diff_nh4(iLu) = MAX(0.0_wp, c_nh4 - Km_nh4)   * RLD_lay * De_lay * tk(iLu) * 86400.0_wp * 14.0_wp
        S_diffusion_pot_day = S_diffusion_pot_day + diff_no3(iLu) + diff_nh4(iLu)
      END IF
      ! The per-layer per-ion availability cap, for mass conservation; the demand
      ! competition below scales all terms by the common ratio. Record 139 mode 1
      ! extends the supply sum with the diffusion term, which is zero when off.
      ! Record 139 mode 2 instead BOUNDS the active (kinetic) term by the transport
      ! delivery, the sum of mass flow and diffusion, so the cap becomes the minimum
      ! of the pool and the minimum of the active term and the delivery - which is
      ! algebraically the minimum of the pool and the bounded active term, since the
      ! delivery is at least that bounded term by construction. The bounded term
      ! reduces to the raw active term, i.e. mode-1 behavior, whenever the mode is
      ! below 2 or the delivery is not the binding constraint.
        act_no3_eff(iLu) = act_no3(iLu)
        act_nh4_eff(iLu) = act_nh4(iLu)
        p_no3(iLu) = MIN(MAX(0.0_wp, NO3_lay(iLu)), act_no3_eff(iLu) + mf_no3(iLu) + diff_no3(iLu))
        p_nh4(iLu) = MIN(MAX(0.0_wp, NH4_lay(iLu)), act_nh4_eff(iLu) + mf_nh4(iLu) + diff_nh4(iLu))
        nud_bound_lay(iLu) = 0
      U_pot_sum = U_pot_sum + p_no3(iLu) + p_nh4(iLu)
      S_active_pot_day   = S_active_pot_day   + act_no3(iLu) + act_nh4(iLu)
      S_massflow_pot_day = S_massflow_pot_day + mf_tot
      ! rec 139: per-layer diagnostics for output/daily_n_uptake_terms.csv
      ! (NEW file, ON-only; computed unconditionally -- cheap, read-only here).
      nud_mf_lay(iLu)   = mf_no3(iLu) + mf_nh4(iLu)
      nud_diff_lay(iLu) = diff_no3(iLu) + diff_nh4(iLu)
      nud_act_lay(iLu)  = act_no3(iLu) + act_nh4(iLu)
      nud_cap_lay(iLu)  = p_no3(iLu) + p_nh4(iLu)
    END DO
    U_pot = 0.0_wp
    ! ================================================================== !
    ! Record 163, the first of its two places (the layered branch; the single-pool
    ! branch below carries the same statement, and the two are mutually exclusive,
    ! so exactly one runs per day and the counters cannot double-count).
    ! Mode 0 is the demand gate, on its own line below. Mode 1, on POST-ANTHESIS
    ! days only, removes that gate and lets the potential be the existing soil-
    ! nitrogen supply law - the record-139 per-layer mass flow, diffusion and
    ! active Michaelis-Menten terms, each already capped by its own layer pool -
    ! times the existing record-155 limb-(b) regulation factor, which is Eqs. (3)
    ! to (5) of Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085,
    ! doi:10.1093/aob/mcr119, with no demand term. That factor is 1.0 unless
    ! record 155 mode 1 or 3 wrote it today, so this is the pure supply law when
    ! that limb is off. The demand itself is NOT altered and is still reported, so
    ! the demand the model would have presented stays visible beside the uptake it
    ! no longer gates. Zero new constants.
    ! Inactive alternative, see Part 1 Table 6.
    ! ================================================================== !
    IF (postanth_uptake_regulation >= 1 .AND. DVS > 0 .AND. dvs2_growth >= 1.0_wp) THEN
      U_pot = U_pot_sum * nur_EN_day
      pur_gain_day  = U_pot - MIN(N_demand_day, U_pot_sum)
      pur_cum_gain  = pur_cum_gain + MAX(0.0_wp, pur_gain_day)
      pur_EN_sum    = pur_EN_sum + nur_EN_day
      pur_EN_min    = MIN(pur_EN_min, nur_EN_day)
      pur_days_post = pur_days_post + 1
      IF (pur_gain_day > 1.0e-12_wp) pur_days_bind = pur_days_bind + 1
    ELSE
    IF (DVS > 0) U_pot = MIN(N_demand_day, U_pot_sum)
    END IF
    U_pot_day = U_pot
    ratio = 0.0_wp
    IF (U_pot_sum > 1.0e-12_wp) ratio = U_pot / U_pot_sum
    DO iLu = 1, NSoil_layers
      NO3_lay(iLu) = NO3_lay(iLu) - ratio * p_no3(iLu)
      NH4_lay(iLu) = NH4_lay(iLu) - ratio * p_nh4(iLu)
      U_act_NO3_day = U_act_NO3_day + ratio * p_no3(iLu)
      U_act_NH4_day = U_act_NH4_day + ratio * p_nh4(iLu)
      nud_real_lay(iLu) = ratio * (p_no3(iLu) + p_nh4(iLu))   ! rec 139
      ! The mass-flow, active and diffusion diagnostic split: within a layer and ion
      ! the cap scales all terms proportionally. Record 139 mode 1 extends the total
      ! with the diffusion term, zero when off. Mode 2 uses the delivery-bounded
      ! active terms instead of the raw ones, so that the split stays consistent with
      ! the cap that actually bounded the availability above.
      tot_ion = act_no3_eff(iLu) + mf_no3(iLu) + diff_no3(iLu)
      IF (tot_ion > 1.0e-15_wp) THEN
        N_uptake_massflow_day  = N_uptake_massflow_day  + ratio * p_no3(iLu) * mf_no3(iLu)      / tot_ion
        N_uptake_active_day    = N_uptake_active_day    + ratio * p_no3(iLu) * act_no3_eff(iLu) / tot_ion
        N_uptake_diffusion_day = N_uptake_diffusion_day + ratio * p_no3(iLu) * diff_no3(iLu)    / tot_ion
      END IF
      tot_ion = act_nh4_eff(iLu) + mf_nh4(iLu) + diff_nh4(iLu)
      IF (tot_ion > 1.0e-15_wp) THEN
        N_uptake_massflow_day  = N_uptake_massflow_day  + ratio * p_nh4(iLu) * mf_nh4(iLu)      / tot_ion
        N_uptake_active_day    = N_uptake_active_day    + ratio * p_nh4(iLu) * act_nh4_eff(iLu) / tot_ion
        N_uptake_diffusion_day = N_uptake_diffusion_day + ratio * p_nh4(iLu) * diff_nh4(iLu)    / tot_ion
      END IF
    END DO
    N_soil_uptake_day = U_act_NO3_day + U_act_NH4_day
    ! rec 139: season-cumulative REALIZED totals + the layer-pool-cap-
    ! binding day counter, gated so OFF touches none of these.
    IF (n_uptake_diffusion >= 1) THEN
      nud_s_mf_cum   = nud_s_mf_cum   + N_uptake_massflow_day
      nud_s_act_cum  = nud_s_act_cum  + N_uptake_active_day
      nud_s_diff_cum = nud_s_diff_cum + N_uptake_diffusion_day
      IF (ANY(p_no3 + 1.0e-12_wp < act_no3 + mf_no3 + diff_no3) .OR. &
          ANY(p_nh4 + 1.0e-12_wp < act_nh4 + mf_nh4 + diff_nh4)) THEN
        nud_days_cap_bound = nud_days_cap_bound + 1
      END IF
    END IF
    ! Record 139 mode 2: the season day-count on which any layer active term was
    ! bounded by delivery, gated so that a lower mode touches nothing.
    ! Record 139 mode 3: the same day-count and the season-cumulative realized
    ! share of the mode-3 diffusion term, gated so that a lower mode touches
    ! nothing.
    IF (dvs2_growth < 1.0_wp) THEN
      N_uptake_pre_cum  = N_uptake_pre_cum  + N_soil_uptake_day
    ELSE
      N_uptake_post_cum = N_uptake_post_cum + N_soil_uptake_day
    END IF
    ! Exact aggregate mirrors for all downstream readers + G8 closure.
    N_soil_min = SUM(NO3_lay)
    N_NH4_soil = SUM(NH4_lay)
    N_soil_uptake_total = N_soil_uptake_total + N_soil_uptake_day
    N_soil_leach_total  = N_soil_leach_total  + N_leach_d
  END BLOCK
  f_Tr_N  = MIN(1.0_wp, trans_act / MAX(0.01_wp, trans_ref_N))   ! diagnostic only
  f_CO2_N = 1.0_wp                                               ! DEPRECATED (Gojon D deferred), as incumbent
  ELSE
  IF (N_soil_min > 0.0_wp) THEN
    ! Critical-N demand (Greenwood 1990, Eq. 3) only when the crop is present (DVS>0);
    ! leaching of NO3 can occur even pre-emergence.
    N_demand_day = 0.0_wp
    IF (DVS > 0) THEN
      ! W must be in t/ha: 1 g/m² = 0.01 t/ha → divide by 100.
      N_plant_total = N_leaf + N_stem + N_root + N_storage + N_dead_leaf
      ! The dead stem, sheath and chaff dry matter stays ABOVE GROUND - it is
      ! harvested straw and the observed above-ground biomass includes it - so it
      ! enters every above-ground dry-matter total. It is identically zero when
      ! record 63 is off.
      AGB_kg_m2 = MAX(0.01_wp, (W_green + W_dead + W_dead_stem + W_stem + &
                MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) / 100.0_wp) ! now t/ha
      ! The critical-nitrogen dilution curve with its published LOW-BIOMASS PLATEAU.
      !   Justes et al. (1994) Ann Bot 74:397-407, read at first hand, Eq. (11) p.401
      !   and Fig. 4 p.401: Nc% = 5.35 W^-0.442, printed valid for
      !   1.55 < W < 12 t ha-1, with the PRINTED constant plateau Nc = 4.4 % below
      !   1.55 t ha-1 (p.402 col. 2: the mean value, 4.4 %, can be retained, the two
      !   curves meeting at an abscissa of 1.55 t ha-1); and the maximum curve
      !   Nmax% = 8.3 W^-0.44, Eq. (17) p.405 and Fig. 5 p.402, printed valid from
      !   1.5 to 14 t ha-1 (p.404 col. 2 bottom). THE PAPER PRINTS NO MAXIMUM VALUE
      !   BELOW 1.5 t ha-1: the 6.6 % read below is an undocumented constant and not
      !   a value of that paper (record 160).
      ! Extrapolating the power law below the plateau would give a small crop an
      ! unbounded critical concentration - 41 % at 0.01 t ha-1 - and hence an
      ! effectively uncapped early-season demand and luxury uptake while the soil
      ! pool is fullest. Uptake stays demand-bounded exactly as before; only the
      ! curve is applied within its published validity range. This is the
      ! nitrogen-status feedback of Devienne-Barret et al. (2000) Ann Bot
      ! 86:995-1005, in which the uptake capacity is limited by the crop position
      ! relative to the critical dilution curve, with no site-fitted constant.
      ! The uptake CAPACITY is bounded by the MAXIMUM dilution curve, the
      ! luxury-storage ceiling, and not by the critical curve: the critical curve is
      ! the minimum concentration for unrestricted growth and defines the nitrogen
      ! nutrition index, while the maximum curve is the upper envelope of observed
      ! crop nitrogen. Capping the demand at the critical curve leaves the crop no
      ! physiological buffer and starves later rapid growth.
      ! ================================================================== !
      ! Record 160: mode 0 is the 6.6 % plateau, on its own line below; mode 1 holds
      ! the printed curve at its own value at the floor,
      ! 8.3 x 1.55**(-0.44) = 6.844 %DM, a continuity construction, since the paper
      ! prints no value there. Zero new constants.
      ! ================================================================== !
      IF (AGB_kg_m2 <= 1.55_wp) THEN
        IF (justes_nmax_floor_basis >= 1) THEN
          Nc_pct = 8.3_wp * 1.55_wp**(-0.44_wp)
          jnf_days_floor = jnf_days_floor + 1
        ELSE
        Nc_pct = 6.6_wp
        END IF
      ELSE
        Nc_pct = 8.3_wp * AGB_kg_m2**(-0.44_wp)
      END IF
      ! The maximum-curve gap is the OUTER CAP ONLY; the demand DRIVER is
      ! component-bounded before anthesis and grain-deficit-driven after it. One
      ! global process, with no site factors.
      BLOCK
        REAL(wp) :: g_nmax_shoot, n_shoot, d_growth, d_reservoir, d_root, d_leaf_maint
        ! The dilution curves are SHOOT-based: the target is on the above-ground
        ! biomass against the SHOOT nitrogen, and the root nitrogen is not subtracted
        ! against a shoot target. Nitrogen buried in dead stem and chaff is still shoot
        ! nitrogen.
        n_shoot = N_leaf + N_dead_leaf + N_dead_stem + N_stem + N_storage
        ! same above-ground basis as AGB_kg_m2 above.
        g_nmax_shoot = MAX(0.0_wp, (Nc_pct / 100.0_wp) * (W_green + W_dead + W_dead_stem + W_stem + &
                  MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) - n_shoot)
        IF (dvs2_growth < 1.0_wp) THEN
          ! The demand is PROSPECTIVE: the previous day realized growth increments
          ! (Jamieson & Semenov 2000 p.24, the first priority is given to the new
          ! increment of green area index) plus the bounded stem luxury reservoir (at
          ! most 1.5 %, p.23) plus a separate root structural demand, with the
          ! shoot-consistent maximum-curve gap as the outer cap.
          ! The leaf share of the prospective demand is bounded by the remaining
          ! leaf-nitrogen holding capacity, as coded:
          !   d_leaf = MIN(n4_sln_target dLAI_prev, MAX(0, sln_max_N LAI - N_leaf))
          ! with the construction target 1.5 g N m-2 (Jamieson & Semenov 2000 Field Crops
          ! Res 68:21-29 p.23) and the holding capacity 3.8 g N m-2 (Silva-Perez et al.
          ! 2020 J Exp Bot 71(7):2299-2311 Fig. 4B p.2305); the full sources are at their
          ! SHARED.f90 declarations.
          IF (leafN_drawdown == 1) THEN
            d_growth = MIN(n4_sln_target * n5_dlai_prev, &
                           MAX(0.0_wp, sln_max_N * lai_growth_diag - N_leaf)) &
                     + n5_stem_nstruct * (n5_dwstem_prev + n5_dwchaff_prev)
          ELSE
            d_growth = n4_sln_target * n5_dlai_prev &
                     + n5_stem_nstruct * (n5_dwstem_prev + n5_dwchaff_prev)
          END IF
          ! The luxury-reservoir DEMAND is read from the same two organ pools the
          ! allocation fills, or the chaff share would be demanded twice and the
          ! post-anthesis uptake inflated. As coded:
          !   d_reservoir = MAX(0, n4_stem_nmax W_stem  - N_stem)
          !               + MAX(0, n4_stem_nmax W_chaff - N_chaff)
          ! with the stem luxury ceiling 0.015 g N g-1 DM (Jamieson & Semenov 2000 Field
          ! Crops Res 68:21-29 p.23). The zeros are positivity guards.
          ! Under record 63 the chaff DEMAND is read against the same ceiling the chaff
          ! ALLOCATION fills, or the reservoir would be under-demanded: 0.017 g N g-1 DM,
          ! derived from Bertheloot, Martre & Andrieu (2008) Plant Physiol Fig. 1D and
          ! Table I at its SHARED.f90 declaration, and not from any observation of this
          ! campaign. With record 63 off the inherited stem ceiling is selected.
          d_reservoir = MAX(0.0_wp, n4_stem_nmax * W_stem  - N_stem) &
                      + MAX(0.0_wp, MERGE(n2_chaff_nmax, n4_stem_nmax, grainN_partition == 1) &
                                    * W_chaff - N_chaff)
          d_root = phen_params%init_rootN_frac * n5_dwroot_prev
          veg_N_demand_day = d_growth + d_reservoir + d_root
          ! ================================================================ !
          ! Record 176 at or above 2: THE PRE-ANTHESIS DEMAND GAINS ITS STOCK-DEFICIT
          ! TERM, the single-pool copy of the statement above. The vegetative demand is a
          ! pure FLUX and is compared with a ceiling that is a STOCK gap, which cannot
          ! bind a daily flux (zero binding days on 42 of 42 decks), while a pure flux
          ! can never recover arrears once the crop has fallen below the critical curve.
          !   Nc% = 4.4 for W at or below 1.55 t ha-1, else 5.35 W^-0.442 - Justes et al.
          ! (1994) Ann Bot 74:397-407, doi:10.1006/anbo.1994.1133, Eq. (11) p.401,
          ! printed valid from 1.55 to 12 t ha-1, with the 4.4 % plateau and the
          ! 1.55 t ha-1 break at p.402 col. 2. The shoot dry-matter and shoot nitrogen
          ! bases are the IDENTICAL expressions the stock cap uses above, so the deficit
          ! and its cap are on ONE basis. The MAXIMUM curve stays as the outer cap. The
          ! term is self-extinguishing and is a DEMAND only, still served through the
          ! unchanged single minimum against supply. Zero new constants.
          ! ================================================================ !
          grain_N_demand_postanth_day = 0.0_wp
          N_demand_day = MIN(veg_N_demand_day, g_nmax_shoot)
          ! ================================================================ !
          ! Record 155 LIMB (a), the single-pool copy: the PRE-ANTHESIS LUXURY DEMAND ON
          ! THE MAXIMUM DILUTION CURVE. The outer stock cap above already carries that
          ! curve and binds on ZERO pre-anthesis days on 42 of 42 released decks, so
          ! swapping the curve in the cap is a provable no-op; what binds is the
          ! vegetative demand itself, the sum of organ-MINIMUM construction targets.
          ! Limb (a) therefore puts the SAME curve on the demand DRIVER, in the marginal
          ! form this code already uses for the CRITICAL curve, c_marg = a(1-b) W^-b with
          ! W floored at 1.55 t ha-1 (Pask 2009 PhD thesis, University of Nottingham,
          ! p.151). That marginal form IS the STICS maximum-curve DEMAND: Brisson et al.
          ! (1998) Agronomie 18:311-346, doi:10.1051/agro:19980501, p.342 Eq. (47) and
          ! p.329, the product of the daily growth rate and the derivative of the maximum
          ! dilution curve. Its Table VIII p.330 prints 8.47 for wheat where this model
          ! codes 8.3, and the model's OWN 8.3 is used, so the record adds no coefficient;
          ! Justes et al. (1994) prints Eq. (17) p.405, Nmax = 8.3 DM^-0.44, so 8.3 IS
          ! the value of that paper. That paper prints no maximum value below 1.5 t ha-1,
          ! so the 6.6 % plateau has no source and record 160 mode 1 replaces it with the
          ! curve's own floor value; limb (a) never reads it either way. Also at first
          ! hand, p.405: the maximum rate of nitrogen accumulation is k times the maximum
          ! growth rate with k = 46.5 DM^-0.44, the primary source for this marginal law.
          ! The crop therefore demands at least enough nitrogen to build TODAY'S OWN new
          ! above-ground dry matter at the MAXIMUM dilution concentration, and is still
          ! never allowed past the maximum-curve stock gap. The increments are the
          ! existing one-day-lagged realized increments. No new constant.
          ! ================================================================ !
          IF (n_uptake_regulation == 1) THEN
            BLOCK
              REAL(wp) :: nur_W_eff, nur_c_max_marg, nur_dem_inc
              nur_W_eff      = MAX(AGB_kg_m2, 1.55_wp)
              nur_c_max_marg = (8.3_wp * (1.0_wp - 0.44_wp) * nur_W_eff**(-0.44_wp)) * 0.01_wp
              nur_D_max_day  = nur_c_max_marg * (nur_dwgreen_prev + n5_dwstem_prev + n5_dwchaff_prev) &
                             + d_root
              nur_dem_inc    = MIN(MAX(veg_N_demand_day, nur_D_max_day), g_nmax_shoot)
              IF (nur_dem_inc > N_demand_day + 1.0e-12_wp) THEN
                nur_days_a_bind = nur_days_a_bind + 1
                nur_cum_a_gain  = nur_cum_a_gain + (nur_dem_inc - N_demand_day)
              END IF
              N_demand_day = nur_dem_inc
            END BLOCK
          END IF
        ELSE
          ! Post-anthesis: the leaf maintenance demand, which delays senescence (de
          ! Ruiter & Brooking 1996 through Jamieson & Semenov 2000 p.24), plus yesterday
          ! unmet grain-nitrogen demand on the harvest-index trajectory. The maintenance
          ! demand is bounded identically, its target being the lesser of the
          ! construction target and the holding capacity, which at the released values
          ! (1.5 below 3.8) is a structural guard only.
          IF (leafN_drawdown == 1) THEN
            d_leaf_maint = MAX(0.0_wp, MIN(n4_sln_target, sln_max_N) * lai_growth_diag - N_leaf)
          ELSE
            d_leaf_maint = MAX(0.0_wp, n4_sln_target * lai_growth_diag - N_leaf)
          END IF
          veg_N_demand_day = d_leaf_maint
          grain_N_demand_postanth_day = n4_grain_deficit_prev
          N_demand_day = MIN(veg_N_demand_day + grain_N_demand_postanth_day, g_nmax_shoot)
        END IF
        ! ================================================================ !
        ! Record 155 LIMB (b), the single-pool copy: the DOWN-REGULATION OF THE
        ! WHOLE-PLANT UPTAKE DEMAND BY THE CULM MOBILE-NITROGEN CONCENTRATION, with its
        ! two ceiling diagnostics.
        !   [N]mob = N_mobile_veg / (W_green + W_stem)        g g-1
        !   E_N    = exp(-beta_N [N]mob),  beta_N = 130
        !   N_demand_day <- N_demand_day E_N
        ! Eq. (4) of Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096 with
        ! the coefficient of Bertheloot, Wu, Cournede & Andrieu (2011) Table 2, applied
        ! on EVERY live day because that model applies it on every day. Pre-anthesis
        ! the pool is near zero here, so the factor is near 1 and its realized value is
        ! REPORTED per day. The denominator is the lamina-plus-stem basis the
        ! record-154 nitrogen factor uses. The critical-curve stock gap beside it is a
        ! DIAGNOSTIC read by no statement.
        ! ================================================================ !
        nur_dem_before   = N_demand_day
        nur_cap_nmax_day = g_nmax_shoot
        ! Record 160: on a sub-floor day, did the raised maximum-curve STOCK GAP bind
        ! today's demand? Diagnostic only: no statement reads the counter and it is
        ! printed in the season summary.
        IF (justes_nmax_floor_basis >= 1 .AND. AGB_kg_m2 <= 1.55_wp) THEN
          IF (N_demand_day >= g_nmax_shoot - 1.0e-12_wp) jnf_days_bind = jnf_days_bind + 1
        END IF
        BLOCK
          REAL(wp) :: nur_nc_pct
          IF (AGB_kg_m2 <= 1.55_wp) THEN
            nur_nc_pct = 4.4_wp
          ELSE
            nur_nc_pct = 5.35_wp * AGB_kg_m2**(-0.442_wp)
          END IF
          nur_cap_nc_day = MAX(0.0_wp, (nur_nc_pct / 100.0_wp) * (W_green + W_dead + W_dead_stem + W_stem + &
                    MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) - n_shoot)
        END BLOCK
        ! Record 157 mode 3: WHAT ENTERS the limb-(b) mobile-nitrogen NUMERATOR. At
        ! modes 0, 1 and 2 it is the positive part of the mobile vegetative pool; at
        ! mode 3 it also carries the shoot nitrate store and the stem and chaff
        ! nitrogen above the structural floor of Jamieson & Semenov (2000) p.23, both
        ! of them circulating amino acids and nitrate in the definition of Bertheloot,
        ! Cournede & Andrieu (2011). The record-157 LAMINA STORE is not counted: that
        ! paper calls the remobilizable nitrogen of a photosynthetic entity only the
        ! nitrogen associated with the photosynthetic apparatus. The DENOMINATOR is
        ! unchanged.
        lns_nmob_t_mobile  = MAX(0.0_wp, N_mobile_veg)
        lns_nmob_t_no3     = 0.0_wp ; lns_nmob_t_stemlux = 0.0_wp ; lns_nmob_t_chaflux = 0.0_wp
        IF (lamina_n_storage == 3) THEN
          lns_nmob_t_no3     = MAX(0.0_wp, N_NO3_shoot)
          lns_nmob_t_stemlux = MAX(0.0_wp, N_stem  - n5_stem_nstruct * W_stem)
          lns_nmob_t_chaflux = MAX(0.0_wp, N_chaff - n5_stem_nstruct * W_chaff)
        END IF
        lns_nmob_num_day = lns_nmob_t_mobile + lns_nmob_t_no3 + lns_nmob_t_stemlux + lns_nmob_t_chaflux
        ! Record 176 at or above 1: Bertheloot, Cournede & Andrieu (2011) p.1090 define
        ! the regulating quantity as the circulating amino acids and nitrate. Stem and
        ! chaff nitrogen held at the luxury ceiling is STORED and not circulating, and
        ! it carries 94-100 % of the pre-anthesis signal, which is why the regulation
        ! factor discards 58.3 % of the pre-anthesis demand and reads the SAME value at
        ! 40 and at 240 kg N ha-1. The numerator returns to the two circulating terms
        ! the model has; the coefficient 130 and the denominator are unchanged.
        ! Record 159: the published ESTIMATE of the mobile pool REPLACES the numerator
        ! formed above. Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot 108(6):1097,
        ! doi:10.1093/aob/mcr125: remobilizable plus mobile nitrogen mass allocated to
        ! a module was estimated as the difference between total nitrogen mass and
        ! structural nitrogen mass, and mobile nitrogen mass was estimated as in
        ! Bertheloot et al. (2008a), assuming that it represented 10 % of the
        ! remobilizable nitrogen mass. Those two sentences solve exactly for
        ! N_rem = (N_tot - N_struct)/(1+r) and N_mob = r N_rem, which is already this
        ! model's own photosynthetic-nitrogen definition under record 154, so the three
        ! terms below ARE the record-154 remobilizable pools, organ for organ and floor
        ! for floor. The mobile fraction is a parameter alias of the existing 0.10, so
        ! no new number enters. Mode 1 is faithful, the estimate being an estimate OF
        ! that single mobile pool; mode 2 adds the explicit transport pools and is the
        ! double count. The denominator below is unchanged.
        IF (n_uptake_regulation == 1) THEN
          nur_Nmob_day = lns_nmob_num_day / MAX(1.0e-9_wp, W_green + W_stem)
          nur_EN_day   = EXP(-n1_beta_N * nur_Nmob_day)
          N_demand_day = N_demand_day * nur_EN_day
          nur_days_live = nur_days_live + 1
          nur_EN_sum    = nur_EN_sum + nur_EN_day
          nur_EN_min    = MIN(nur_EN_min, nur_EN_day)
        END IF
        nur_dem_after = N_demand_day
      END BLOCK
    END IF
    ! Michaelis-Menten uptake scaled by transpiration (Lemaire & Gastal 2009 Ann
    ! Bot 106:35, p.38 Eq. 2).
    !
    ! The CO2 route works through the model's own iteratively computed daily
    ! intercellular CO2 concentration. The mechanism is Gojon (2022) Trends Plant
    ! Sci, pp.4-5, Fig. 1C and Fig. 2B: elevated CO2 raises the intercellular
    ! concentration, which lowers the photorespiratory ratio, which exports less
    ! 2-oxoglutarate through the malate-oxaloacetate shuttle to the cytosol, which
    ! leaves less reductant for cytosolic nitrate reductase and down-regulates the
    ! transporters; that paper reports elevated CO2 significantly reducing root
    ! nitrate uptake by almost 50 %. Zhao et al. (2021) Plant Physiol 187:1812,
    ! Table 1, measure a 5-9 % decrease in nitrogen assimilation at 800 against
    ! 400 ppm through the same ratio.
    ! The reference intercellular concentration is read from driver_N.csv, the
    ! default 300 umol mol-1 being typical of wheat near 400 ppm ambient.
    ! The transpiration mass-flow delivery of nitrate (Lemaire & Gastal 2009 Eq. 2;
    ! Gojon 2022 Fig. 1B) is the transpired water volume times the soil-solution
    ! nitrate concentration, and it is the ONLY route by which CO2 affects uptake
    ! here, through stomatal closure. Beside it stands the root-scaled active
    ! transporter term, which is energy-driven and transpiration-independent. The
    ! total is bounded by the critical-nitrogen demand AND by the available soil
    ! nitrogen, for mass conservation. The empirical direct-CO2 multiplier is NOT
    ! applied: it would stand in for the transporter down-regulation of that
    ! review, whose molecular mechanism is unknown.
    BLOCK
      REAL(wp) :: soil_water_mm_up, no3_soln_gN_L, S_massflow, S_active, U_pot, L_pot, dsum, sfac
      REAL(wp) :: u_max_dw, NO3_rz_up, f_Nstatus, a_root
      REAL(wp) :: w_shoot_tha, n_shoot_up, nc_crit_pct, nni_up
      INTEGER  :: k_nni
      soil_water_mm_up = MAX(1.0_wp, vwc_root * MAX(0.1_wp, root_dep) * 1000.0_wp)  ! L water m-2
      no3_soln_gN_L    = N_soil_min / soil_water_mm_up                              ! g N L-1 (NO3)
      S_massflow = trans_act * no3_soln_gN_L                                        ! Gojon B (mass flow)
      ! ---------------------------------------------------------------- !
      ! ROOT-SCALED ACTIVE UPTAKE.
      !
      ! Unit chain:
      !   u_max_fw   2.6 umol NO3- g-1 FW h-1  [Goyal & Huffaker 1986 Plant Physiol
      !              82:1051, p.1055 from their Fig. 1c and 1d, WHEAT; the 25 C basis
      !              is carried as stated, with no temperature term, because the only
      !              alternative source has a cultivar-dependent sign]
      !   x fw_dw_root 18.6 g FW per g DW     [Khan et al. 2024 Heliyon Table 3]
      !              = 48.36 umol NO3- g-1 DW h-1
      !   x 14e-6 g N umol-1                  = 6.770e-4 g N g-1 DW h-1
      !   x 24 h d-1                          = 1.625e-2 g N g-1 DW d-1
      !   x the root dry mass (g DW m-2) and the dimensionless factors, giving the
      !     active uptake in g N m-2 d-1.
      ! ---------------------------------------------------------------- !
      u_max_dw = u_max_fw * fw_dw_root * 14.0e-6_wp * 24.0_wp
      ! The root-zone nitrate concentration (mol m-3), the same expression as the
      ! nitrate-partition block, recomputed here on the PRE-uptake pool, while that
      ! block runs after uptake on the depleted pool: one formula, each on its own
      ! point in the daily sequence.
      NO3_rz_up = N_soil_min / MAX(1.0e-6_wp, vwc_root * MAX(0.1_wp, root_dep)) / 14.0_wp
      ! The nitrogen-status function of the nitrogen nutrition index: the whole-plant
      ! nitrogen-status feedback on root high-affinity uptake. A root-concentration
      ! basis is not usable, because the root pins at its structural cap while the
      ! shoot is nitrogen-deficient, which is an allocation artefact and not plant
      ! repletion. The mapping is taken directly from the four MEASURED wheat points,
      ! piecewise linear and clamped:
      !   (NNI 0.50, 1.000) (1.12, 0.877) (1.53, 0.421) (1.95, 0.175)
      ! Taulemesse et al. (2015) PLoS ONE 10(3):e0120291,
      ! doi:10.1371/journal.pone.0120291: the abscissae are Table 1 p.7, the nitrogen
      ! nutrition index at flowering of the four treatments with their standard
      ! errors; the ordinates are the printed sentence p.8, that relative to the
      ! highest treatment the daily nitrate uptake fluxes were 5.7, 5 and 2.4 times
      ! higher at the three lower ones, normalized to the lowest.
      ! Three statements are assumptions and not measurements: (a) at-flowering
      ! measurements are applied over the whole season; (b) that paper's own caveat,
      ! that the index was developed for field crops and should be used carefully
      ! under their semi-hydroponic conditions; and (c) the erratum to Justes et al.
      ! (1996) Ann Bot 78:267 could not be retrieved.
      !
      ! The index is the shoot nitrogen percentage over the critical concentration,
      ! with the CRITICAL dilution curve
      !   Nc% = 5.35 W^(-0.442) for shoot dry matter from 1.55 to 12 t ha-1,
      !   Nc% = 4.4 (the plateau) below 1.55 t ha-1
      ! (Justes et al. 1994 Ann Bot 74:397, doi:10.1006/anbo.1994.1133). This is the
      ! CRITICAL curve, which defines the index, and it stands alongside the maximum
      ! cap curve of the demand block. The shoot sums are the demand block's own
      ! conventions: shoot dry matter is the green lamina, dead lamina, stem and
      ! reproductive pools (g m-2, divided by 100 for t ha-1) and shoot nitrogen the
      ! matching four nitrogen pools, both evaluated on the beginning-of-day state, a
      ! stated one-day lag. The above-ground basis is the same one the demand uses.
      w_shoot_tha = (W_green + W_dead + W_dead_stem + W_stem + &
                MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) / 100.0_wp
      n_shoot_up  = N_leaf + N_dead_leaf + N_stem + N_storage
      nni_up = 0.0_wp
      IF (w_shoot_tha > 1.0e-6_wp) THEN
        IF (w_shoot_tha <= 1.55_wp) THEN
          nc_crit_pct = 4.4_wp
        ELSE
          nc_crit_pct = 5.35_wp * w_shoot_tha**(-0.442_wp)
        END IF
        nni_up = (n_shoot_up / w_shoot_tha) / nc_crit_pct
      END IF
      IF (nni_up <= fNstat_nni(1)) THEN
        f_Nstatus = fNstat_f(1)              ! = 1.000 clamped for NNI <= 0.50
      ELSE IF (nni_up >= fNstat_nni(4)) THEN
        f_Nstatus = fNstat_f(4)              ! = 0.175 clamped for NNI >= 1.95
      ELSE
        f_Nstatus = fNstat_f(4)
        DO k_nni = 1, 3
          IF (nni_up <= fNstat_nni(k_nni+1)) THEN
            f_Nstatus = fNstat_f(k_nni) + (fNstat_f(k_nni+1) - fNstat_f(k_nni)) * &
                        (nni_up - fNstat_nni(k_nni)) / (fNstat_nni(k_nni+1) - fNstat_nni(k_nni))
            EXIT
          END IF
        END DO
      END IF
      ! The living-root uptake capacity fraction (Zheng et al. 2022 Front Plant Sci
      ! 12:814658, Tables 2-3: the product of the root-length-density fraction and
      ! the activity fraction is 0.28-0.48, midpoint 0.38, at 20 days after
      ! anthesis). It is 1 up to anthesis, since no readable primary source shows a
      ! pre-anthesis decline, and after anthesis it falls linearly in DAYS AFTER
      ! ANTHESIS at (1 - 0.38)/20 = 0.031 d-1, floored at 0.28, the measured lower
      ! bound held beyond the measured window. The calendar-day basis is forced by
      ! that paper reporting no site temperatures or degree-days.
      a_root = 1.0_wp
      IF (daa_root > 0) a_root = MAX(aroot_floor, 1.0_wp - aroot_decline_d * REAL(daa_root, wp))
      ! The law (no cap constant anywhere; soil pool still bounds realized
      ! uptake via the competition step below):
      S_active = u_max_dw * f_Nstatus * a_root * W_root * NO3_rz_up / (Km_c_no3 + NO3_rz_up)
      f_Nstatus_diag = f_Nstatus
      a_root_diag    = a_root
      ! rec 134/: same vegetative-basis NNI and f_Nstatus mirror
      ! as the layered branch above (incumbent soilN_layered = 0 branch).
      BLOCK
        REAL(wp) :: w_veg_tha2, n_veg_up2, nc_crit_veg2
        w_veg_tha2 = (W_green + W_dead + W_dead_stem + W_stem) / 100.0_wp
        n_veg_up2  = N_leaf + N_dead_leaf + N_stem
        nap_nni_veg = 1.0_wp
        IF (w_veg_tha2 > 1.0e-6_wp) THEN
          IF (w_veg_tha2 <= 1.55_wp) THEN
            nc_crit_veg2 = 4.4_wp
          ELSE
            nc_crit_veg2 = 5.35_wp * w_veg_tha2**(-0.442_wp)
          END IF
          nap_nni_veg = (n_veg_up2 / w_veg_tha2) / nc_crit_veg2
        END IF
      END BLOCK
      nap_f_Nstatus_day = f_Nstatus
      ! (6) Potential NO3 uptake: demand-bounded only (NOT soil-bounded; competition does that).
      U_pot = 0.0_wp
      ! ================================================================ !
      ! Record 163, the second of its two places (the single-pool branch, mutually
      ! exclusive with the layered one above). The statement is identical, with this
      ! branch's own supply quantity, the positive part of the mass-flow plus active
      ! terms, in place of the summed layer potential. Mode 0 is the demand gate on
      ! its own line. Zero new constants. Inactive alternative, see Part 1 Table 6.
      ! ================================================================ !
      IF (postanth_uptake_regulation >= 1 .AND. DVS > 0 .AND. dvs2_growth >= 1.0_wp) THEN
        U_pot = MAX(0.0_wp, S_massflow + S_active) * nur_EN_day
        pur_gain_day  = U_pot - MIN(N_demand_day, MAX(0.0_wp, S_massflow + S_active))
        pur_cum_gain  = pur_cum_gain + MAX(0.0_wp, pur_gain_day)
        pur_EN_sum    = pur_EN_sum + nur_EN_day
        pur_EN_min    = MIN(pur_EN_min, nur_EN_day)
        pur_days_post = pur_days_post + 1
        IF (pur_gain_day > 1.0e-12_wp) pur_days_bind = pur_days_bind + 1
      ELSE
      IF (DVS > 0) U_pot = MIN(N_demand_day, MAX(0.0_wp, S_massflow + S_active))
      END IF
      ! (5) Potential NO3 leaching (Burns 1975 piston flow) from the SAME pre-removal NO3 pool.
      L_pot = 0.0_wp
      IF (drain_mm_d > 0.0_wp) L_pot = N_soil_min * MIN(1.0_wp, drain_mm_d / soil_water_mm_up)
      ! (7) Simultaneous competition: scale BOTH if uptake + leaching exceed available NO3.
      dsum = U_pot + L_pot
      IF (dsum > N_soil_min .AND. dsum > 0.0_wp) THEN
        sfac = N_soil_min / dsum
        N_soil_uptake_day = U_pot * sfac
        N_leach_d         = L_pot * sfac
      ELSE
        N_soil_uptake_day = U_pot
        N_leach_d         = L_pot
      END IF
      ! Report mass-flow/active components scaled to the actual (competed) uptake.
      sfac = 0.0_wp
      IF (S_massflow + S_active > 0.0_wp) sfac = N_soil_uptake_day / (S_massflow + S_active)
      N_uptake_massflow_day = S_massflow * sfac
      N_uptake_active_day   = S_active   * sfac
      ! daily diagnostics: PRE-competition potentials (ceiling-finding
      ! does S_active ever bind U_pot?) -> daily_nitrogen.csv.
      S_active_pot_day   = S_active
      S_massflow_pot_day = S_massflow
      U_pot_day          = U_pot
    END BLOCK
    f_Tr_N  = MIN(1.0_wp, trans_act / MAX(0.01_wp, trans_ref_N))   ! diagnostic only (mass flow now explicit)
    f_CO2_N = 1.0_wp                                               ! DEPRECATED: removed from core (Gojon D deferred)
    ! diagnostics (§5.4): pre/post-anthesis uptake split
    IF (dvs2_growth < 1.0_wp) THEN
      N_uptake_pre_cum  = N_uptake_pre_cum  + N_soil_uptake_day
    ELSE
      N_uptake_post_cum = N_uptake_post_cum + N_soil_uptake_day
    END IF
    N_soil_min = MAX(0.0_wp, N_soil_min - N_soil_uptake_day - N_leach_d)
    N_soil_uptake_total = N_soil_uptake_total + N_soil_uptake_day
    N_soil_leach_total  = N_soil_leach_total  + N_leach_d
  END IF
  END IF   ! soilN_layered branch (G9)
  ! (8) The in-season soil-nitrogen closure (urea, ammonium, nitrate, residue and
  !     microbial): the pool change equals the external inputs (fertilizer,
  !     background mineralization and in-season residue nitrogen) minus the
  !     removals (leaching, uptake, volatilization and denitrification, the last
  !     two being gaseous exports). Hydrolysis, nitrification, decomposition,
  !     immobilization and microbial turnover are internal transfers, and
  !     pre-sowing residue is already in the pools. The surface ammoniacal pool
  !     is part of the pool sum. The residual should be about zero.
  IF (soil_n_mech >= 1) THEN
    ! humus N inside the pool sum; humus mineralisation is an
    ! INTERNAL transfer and is NOT subtracted as an external input.
    soilN_balance_resid = (N_urea_soil + N_NH4_soil + TAN_surf + N_soil_min + res_fast_N + res_slow_N &
                           + soil_mic_N + NHUM_soil) &
                          - soilN_before_day - N_fert_day - residue_N_in_day &
                          + N_leach_d + N_soil_uptake_day + N_volat_day + N_denit_day
  ELSE
  soilN_balance_resid = (N_urea_soil + N_NH4_soil + TAN_surf + N_soil_min + res_fast_N + res_slow_N + soil_mic_N) &
                        - soilN_before_day - N_fert_day - N_mineral_day - residue_N_in_day &
                        + N_leach_d + N_soil_uptake_day + N_volat_day + N_denit_day
  END IF

  ! -------------------------------------------------------------------- !
  ! Nitrogen allocation follows wheat physiology:
  !
  ! PRE-ANTHESIS (development stage below 1): both soil uptake and de novo
  ! reduction build the vegetative pools (leaf, stem, root), and the tabled
  ! storage fraction is zero before anthesis. Martre et al. (2003) Plant Cell
  ! Environ 26:561-571, Fig. 1: pre-anthesis nitrogen accumulation in vegetative
  ! organs is the primary source of final grain nitrogen in wheat.
  !
  ! POST-ANTHESIS (development stage 1 and above): grain nitrogen comes primarily
  ! from remobilization. Post-anthesis soil uptake maintains leaf nitrogen and
  ! delays senescence; it is not channelled through the tabled storage fraction,
  ! which would send almost all of it to the reproductive pool and give
  ! unrealistic protein concentrations. The remaining 10-30 % from post-anthesis
  ! uptake maintains the vegetative nitrogen and reaches the grain by
  ! remobilization (Martre et al. 2003). The 70-90 % remobilization share has a
  ! citation that could not be verified: neither of the two papers of that author
  ! and year held here carries it and no such file is on this machine, and the
  ! rate it supported is not used on any live path.
  ! -------------------------------------------------------------------- !
  W_storage_prev_day = W_storage   ! capture grain biomass before growth (Process 5 sink)
  W_green_prev_day = W_green; W_stem_prev_day = W_stem; W_root_prev_day = W_root; W_chaff_prev_day = W_chaff
  ! Record 92: today's labeled reserve deposition is accumulated by BOTH routes -
  ! the pre-anthesis stem banking in the allocation and the sink-full surplus in
  ! the yield-component routine - so it is zeroed here, ahead of both. Nothing is
  ! written at record 92 = 0.
  IF (stem_wsc_mech >= 1) wsc_dep_day = 0.0_wp
  ! Record 105: today's terms are zeroed ahead of the draw, which may not run
  ! today, and the alternative-clock increment is handed in. Nothing here
  ! executes at record 105 = 0 or 1.
  ! Record 105 mode 2: today's mode-2 terms are zeroed ahead of the release, which
  ! may not run today, and the same clock increment is handed in. On a day the
  ! release branch does not run, the pre-release pool and the flux columns are
  ! zero, so the ledger identity is tested from the post-release column across
  ! rows, which is written every day.
  IF (wsc_reserve_law == 2) THEN
    wsl2_rel_day = 0.0_wp ; wsl2_relg_day = 0.0_wp ; wsl2_relf_day = 0.0_wp
    wsl2_gain_day = 0.0_wp ; wsl2_lconv_day = 0.0_wp ; wsl2_cres_day = 0.0_wp
    wsl2_wres_pre = 0.0_wp ; wsl2_tau_day = 0.0_wp ; wsl2_taurem_day = 0.0_wp
    wsl2_bound_day = 0
    wsl_ddvs2 = dvs2_rate_growth
  END IF
  ! Record 148: today's term-A, term-B and deficit-draw terms are zeroed ahead of
  ! the allocation call, since any of them may not run today. Nothing here
  ! executes at record 148 = 0.
  IF (wsc_fate_law >= 1) THEN
    wfl_a_refused_day = 0.0_wp ; wfl_b_released_day = 0.0_wp ; wfl_draw_def_day = 0.0_wp
  END IF
  ! Record 121: under the payment gate today's nitrate reduction is resolved BEFORE
  ! the carbon is allocated. Its inputs, the day soil uptake from the block above
  ! and the day biochemical reduction from the canopy pass, are final at this
  ! point. Nothing here executes at record 121 = 0, where the block runs at its
  ! later position.
  IF (n_payment_gate >= 1) THEN
    CALL andrews_nitrate_partition_day()
    npg_n_assim_root_in  = N_assim_root_day      ! pays the root tier first step (1)
    npg_n_assim_shoot_in = N_assim_shoot_day
  END IF
  ! Record 128: the shoot critical-nitrogen curve reads TODAY'S above-ground dry
  ! matter, set here just before the allocation, on the IDENTICAL basis the live
  ! nitrogen-nutrition-index term uses. No new above-ground definition is
  ! introduced. Modes 1 to 3 all need it.
  IF (n_construction_status >= 1) THEN
    ncs_W_t_ha = (W_green + W_dead + W_dead_stem + W_stem + &
                  MERGE(W_grain + W_chaff, W_storage, gnf_critN_basis)) / 100.0_wp
  END IF
  ! Record 172: publish today's above-ground dry matter and shoot nitrogen on the
  ! IDENTICAL basis the live nitrogen-nutrition-index term uses, for the
  ! deposition law inside the allocation. Nothing executes and nothing is read at
  ! record 172 = 0.
  CALL allocate_daily_biomass(dvs2_growth, t_leaf_growth_mean, An_growth_day, W_green, W_dead, W_stem, W_root, W_storage, &
    rm_growth_day, rg_growth_day)
  dW_storage_day = MAX(0.0_wp, W_storage - W_storage_prev_day)   ! daily reproductive DM gain (g m-2)
  ! Record 89: split TODAY'S spike dry-matter increment among the live axis
  ! classes, on each axis's OWN spike growth period (Fischer 1985 p.449), and
  ! accumulate the window diagnostics. It must run BEFORE the grain number is set
  ! below and AFTER the allocation, which publishes the day above-ground and
  ! spike increments on the Fischer basis. It returns on its first line when
  ! record 89 = 0.
  ! Then the reproductive increment is split into grain and ear or chaff, the
  ! grain number is set at the end of the ear-construction window, and the grain
  ! is sink-limited by the kernel number times the potential kernel mass, with
  ! the grain and chaff masses summing to the reproductive pool.
  ! The air temperature is the SCREEN AIR daily mean at this point, normalized
  ! upstream. The degree-days of Fischer (1985) are on screen air temperature:
  ! p.449, Materials and Methods, temperature was measured in standard
  ! meteorological screens.
  CALL update_spike_axes(iD, T_air_day, dvs2_growth)
  ! Record 91: the same air temperature and the model day are handed to the
  ! grain-filling law. Both are ignored when record 91 = 0.
  ! Record 117: today's development-stage increment is handed in only when the
  ! record is on, exactly as record 105 hands in its own. THE WITHIN-WINDOW
  ! DISTRIBUTION IS AN IMPLEMENTATION ASSUMPTION AND NOT A PUBLISHED LAW - the
  ! sources give a TOTAL and not a time course - and it is never tuned.
  ! The reproductive pool is handed in so that the routine can book the
  ! sink-full surplus (pool to stem) and the reserve draw (stem to grain) on the
  ! reproductive pool itself.
  ! Record 130: the grain nitrogen entering today's dry-matter allocation is the
  ! reproductive nitrogen minus the chaff nitrogen AS OF YESTERDAY, because the
  ! grain-nitrogen source-sink routine for today has not run yet, being called
  ! later in this day loop. It is set unconditionally and read only when
  ! record 130 is at or above 1.
  gnm_N_grain_in = MAX(0.0_wp, N_storage - N_chaff)
  CALL allocate_yield_components_daily(dvs2_growth, dW_storage_day, W_stem, W_grain, W_chaff, &
    kernel_number_m2, grain_number_set, dW_grain_day, T_air_day, iD, W_storage)
  ! Record 105, on only: the day terms and the RESERVE LEDGER identity, the pool
  ! today equaling the pool yesterday plus the deposition minus the draw, from
  ! the five and only five writers of the pool - the record-101 store, the
  ! record-55 and 56 bank, the record-92 post-anthesis label, the record-121
  ! reserve payment and the draw - all of which ran above. The bound code is
  ! 0 none, 1 deficit, 2 pool, 3 stem, 4 ceiling.
  ! Record 117, on only: one row per model day, the whole sink on one line, so
  ! that supply limitation and carbon conservation can be tested FROM THIS FILE
  ! ALONE. It is a model internal and is never placed beside an observation. The
  ! chaff mass is printed so that the daily identity between its difference and
  ! the amount taken is checkable across rows on every day the loss-only
  ! record-63 turnover path is not running.
  ! Record 105 mode 2, on only: the day release and split terms and the reserve
  ! ledger, with the columns needed to test the four conservation identities from
  ! this file alone; the run.log summary carries their maxima at full double
  ! precision, because this print format cannot resolve 1e-12. The deposition is
  ! that of the three labeled routes plus the record-101 store when record 101
  ! is on, the same five-writer set as mode 1. The bound code is 0 none, 2 pool,
  ! 3 stem, 4 window ceiling.
  IF (wsc_reserve_law == 2) THEN
    BLOCK
      REAL(wp) :: dep_today2, resid2
      dep_today2 = wsc_dep_day
      IF (surplus_fate_mech >= 1) dep_today2 = dep_today2 + sfl_dw_wsc_day
      ! Record 148: under modes 1 and 3 the DEFICIT limb is reachable again at
      ! record 105 = 2, the limbs being sequential rather than exclusive, so the
      ! ledger must see its draw as well as the release. The term is an exact zero at
      ! modes 0 and 2.
      resid2 = W_reserve - (wsl2_wres_prev + dep_today2 - wsl2_rel_day - wfl_draw_def_day)
      IF (ABS(resid2) > wsl2_i4_max) wsl2_i4_max = ABS(resid2)
      IF (iD == 1) THEN
        OPEN (NEWUNIT=unit_wsl2, FILE='output/daily_wsc_mode2.csv', STATUS='REPLACE', ACTION='WRITE')
        WRITE (unit_wsl2,'(A)') 'day,tau,tau_rem,window_member,W_reserve_pre_g_m2,dep_today_g_m2_d,' // &
          'rel_g_m2_d,rel_grain_g_m2_d,gain_g_m2_d,loss_conv_g_m2_d,rel_fate_g_m2_d,' // &
          'cres_respired_g_m2_d,W_reserve_post_g_m2,W_stem_g_m2,bound'
      END IF
      WRITE (unit_wsl2,'(I0,2(",",G14.6),",",I0,10(",",G14.6),",",I0)') iD, wsl2_tau_day, wsl2_taurem_day, &
        wsc_rel_window_member, wsl2_wres_pre, dep_today2, wsl2_rel_day, wsl2_relg_day, wsl2_gain_day, &
        wsl2_lconv_day, wsl2_relf_day, wsl2_cres_day, W_reserve, W_stem, wsl2_bound_day
      wsl2_wres_prev = W_reserve
    END BLOCK
  END IF
  ! Record 148, on only: one row per model day, the whole fate of the reserve on
  ! one line, so that term A (the refused deposition, which leaves the stem to
  ! the record-129 loss accumulator), term B (released at unity) and the
  ! pre-gate DEFICIT draw are separable from this file alone, and the pool
  ! trajectory - peak, gate, harvest - is readable without the mode-2 ledger. It
  ! is a model internal and is never placed beside an observation.
  IF (wfl_open) THEN
    WRITE (unit_wfl,'(I0,",",I0,",",G14.6,",",I0,",",I0,8(",",G14.6))') iD, &
      MERGE(iD - wsc_day_anth, -999, wsc_day_anth > 0), dvs2_growth, wsc_fate_law, gfb_end_bound, &
      wsl2_wres_pre, wfl_a_refused_day, wfl_b_released_day, wfl_draw_def_day, cres_respired, &
      W_reserve, W_stem, W_reserve / MAX(1.0e-9_wp, W_stem)
  END IF
  ! Mechanistic mineral uptake (mass flow plus the effective active term) and
  ! grain loading, keyed to the organ growth increments.
  ! Record 114: the five arguments are the model's OWN SAME-DAY DELIVERY LEDGER
  ! (root, leaf, stem with the reserve excluded, chaff and grain), each
  ! non-negative by construction, so the caller positivity clips are not needed
  ! on this path. The grain dry-matter increment is not passed here, because it
  ! carries the reserve draw the remobilization blocks already move.
    CALL update_minerals_daily(trans_act, vwc_root, root_max_depth, W_root, W_grain, dvs2_growth, &
      MAX(0.0_wp, W_root - W_root_prev_day), MAX(0.0_wp, W_green - W_green_prev_day), &
      MAX(0.0_wp, W_stem - W_stem_prev_day), MAX(0.0_wp, W_chaff - W_chaff_prev_day), dW_grain_day, &
      dW_grain_day)
  ! Record 131, on only: the daily mineral-supply trace. Nothing here is read by
  ! any physics path.
  IF (mineral_supply_law == 1) THEN
    BLOCK
      REAL(wp) :: msl_shoot_dm, msl_shoot_content, msl_shoot_conc
      INTEGER :: i_msl
      IF (.NOT. msl_diag_open) THEN
        OPEN (NEWUNIT=unit_msl, FILE='output/daily_mineral_supply.csv', STATUS='REPLACE', ACTION='WRITE')
        WRITE (unit_msl,'(A)') 'day,mineral,dvs2,transpiration_mm_d,soln_mg_L,uptake_massflow_mg_m2,' // &
          'uptake_active_mg_m2,uptake_total_mg_m2,soil_labile_mg_m2,shoot_content_mg_m2,shoot_conc_mg_kg'
        msl_diag_open = .TRUE.
      END IF
      msl_shoot_dm = W_green + W_stem + W_chaff
      DO i_msl = 1, NMIN
        msl_shoot_content = M_leaf(i_msl) + M_stem(i_msl) + M_chaff(i_msl)
        msl_shoot_conc = 0.0_wp
        IF (msl_shoot_dm > 1.0e-9_wp) msl_shoot_conc = 1000.0_wp * msl_shoot_content / msl_shoot_dm
        WRITE (unit_msl,'(I0,",",A,10(",",G14.6))') iD, mineral_names(i_msl), dvs2_growth, trans_act, &
          M_soln(i_msl), M_uptake_mf(i_msl), M_uptake_ac(i_msl), M_uptake_mf(i_msl) + M_uptake_ac(i_msl), &
          M_soil_labile(i_msl), msl_shoot_content, msl_shoot_conc
      END DO
    END BLOCK
  END IF
  ! Record 140, on only: the daily mineral-uptake-terms trace. The mass-flow,
  ! diffusion and transporter columns are REALIZED, i.e. post-pool-cap, per-term
  ! uptake, and are zero whenever the record is off or the element is calcium or
  ! not zinc. The cap column is the labile pool REMAINING after today's uptake,
  ! i.e. context for whether the pool is binding, and the realized column is the
  ! day TOTAL realized uptake across every term - the single
  ! conservation-consistent total this module's own uptake variable represents.
  ! There is no separate whole-plant demand-competition stage downstream, unlike
  ! the nitrogen-uptake law of record 139, so the cap and the realized total are
  ! not required to match.
  IF (msm_open) THEN
    DO i_msm = 1, NMIN
      WRITE (unit_msm,'(I0,",",A,5(",",G14.6))') iD, mineral_names(i_msm), M_uptake_mf(i_msm), &
        M_uptake_diff(i_msm), M_uptake_hats(i_msm), M_soil_labile(i_msm), &
        M_uptake_mf(i_msm) + M_uptake_ac(i_msm) + M_uptake_diff(i_msm)
    END DO
  END IF
  ! Record 132, on only: the daily leaf-nitrogen holding-capacity trace. Nothing
  ! here is read by any physics path.
  ! Record 133, on only: the daily lamina-nitrogen-turnover trace. Nothing here
  ! is read by any physics path.

  ! -------------------------------------------------------------------- !
  ! The root and shoot nitrate-reduction partition is the internal procedure
  ! andrews_nitrate_partition_day in the CONTAINS section. At record 121 = 0 it
  ! is executed HERE, after the allocation; at record 121 at or above 1 it is
  ! executed BEFORE the allocation, so that today's reduced nitrogen pays today
  ! carbon.
  ! -------------------------------------------------------------------- !
  IF (n_payment_gate == 0) CALL andrews_nitrate_partition_day()

  IF (dvs2_growth < 1.0_wp) THEN
    ! PRIORITY allocation, in the ordering of Jamieson & Semenov (2000) p.24: (1)
    ! the leaf to its specific-leaf-nitrogen target, including the day expansion
    ! through the previous-day increment; (2) the stem and chaff structural 0.5 %;
    ! (3) the root structural demand; (4) the stem LUXURY up to 1.5 %; (5) the
    ! overflow REMAINS unassimilated in the shoot nitrate pool, i.e. vacuolar
    ! storage, which is mass-conservative and represents demand-limited
    ! assimilation. Fixed tabled fractions instead force the stem to 4.4-5.1 % N,
    ! three times the maximum of 1.5 % that paper states on p.23.
    BLOCK
      REAL(wp) :: pool, cap, take, mobile_in
      ! The mobile vegetative nitrogen, resorbed before anthesis, is available for
      ! new construction ahead of the day assimilation overflow; whatever is not used
      ! returns to the mobile pool and never to the shoot nitrate store, since it is
      ! reduced organic nitrogen and not nitrate.
      IF (n_payment_gate >= 1 .AND. npg_live_today) THEN
        ! rec 121: the gate already drew the mobile pool and today's reduced N and
        ! booked the structural payments; the ladder receives what it left.
        mobile_in = npg_mobile_left
        N_mobile_veg = 0.0_wp
        pool = npg_pool_out
      ELSE
      mobile_in = N_mobile_veg
      N_mobile_veg = 0.0_wp
      pool = N_gain_total + mobile_in
      END IF
      N_alloc_leaf_day = 0.0_wp; N_alloc_stem_day = 0.0_wp
      N_alloc_root_day = 0.0_wp; N_alloc_storage_day = 0.0_wp
      ! The prospective leaf allocation cap is bounded by the leaf nitrogen-holding
      ! capacity, as coded:
      !   cap = MAX(0, MIN(n4_sln_target (LAI + dLAI_prev), sln_max_N LAI) - N_leaf)
      ! The leaf area index here is the last REALIZED one, today being computed after
      ! the allocation, which is the conservative choice for a capacity bound. The
      ! constants are the construction target 1.5 g N m-2 (Jamieson & Semenov 2000
      ! Field Crops Res 68:21-29 p.23) and the holding capacity 3.8 g N m-2
      ! (Silva-Perez et al. 2020 J Exp Bot 71(7):2299-2311 Fig. 4B p.2305, the
      ! maximum measured nitrogen per unit leaf area); the full sources are at their
      ! SHARED.f90 declarations. This bounds every pre-anthesis nitrogen INFLOW to
      ! the leaf; no concentration output is capped.
      ! ================================================================ !
      ! Record 157: the PER-AREA CEILING of this lamina tier. Mode 0 is the
      ! construction TARGET 1.5 g N m-2 (Jamieson & Semenov 2000 p.23), which the
      ! record-121 and record-128 payment gate has ALREADY paid the lamina to, so the
      ! tier headroom is near zero and the whole pre-anthesis surplus runs past it
      ! into the stem luxury pool and then into the shoot nitrate store. Modes 1 and
      ! 2 use instead the lamina HOLDING CAPACITY 2.5 g N m-2 of leaf area, the
      ! printed sentence of Pask (2009) PhD thesis, University of Nottingham, p.29,
      ! already carried by records 132 and 133 mode 2 - no new constant, the same
      ! area basis, the same statement, and the record-62 branch minimum untouched.
      ! At mode 0 the ceiling variable is set to the construction target itself, so
      ! the arithmetic is on the same operand.
      ! ================================================================ !
      lns_sln_ceil_used = n4_sln_target
      IF (lamina_n_storage >= 1) lns_sln_ceil_used = sln_max_N_lnh
      ! ================================================================ !
      ! Record 164: AT WHICH LEVEL OF THE CANOPY DOES THE 2.5 g N m-2 CEILING APPLY?
      ! Mode 0 leaves the line above untouched and this block does not execute.
      ! Mode 1 applies THE SAME PRINTED CONSTANT at the level the source's own sentence
      ! names - especially in the upper leaves of the canopy (Pask 2009 PhD thesis,
      ! University of Nottingham, p.29) - through the MODEL'S OWN existing exponential
      ! canopy nitrogen profile, so that the canopy-mean capacity the tier below uses
      ! becomes
      !     sln_max_N_lnh / (n_top / n_mean)
      ! with the top value the same pure profile function that writes the n_top_g_m2
      ! column of output/daily_canopy_n_profile.csv and the mean the leaf nitrogen
      ! over the leaf area index, evaluated on the operands this statement itself
      ! multiplies. The leaf area index used is the last REALIZED one, which is the
      ! value the gas-exchange leaf area is set equal to earlier in this same day
      ! loop, so this is the printed diagnostic's own ratio; using the ceiling's own area
      ! basis is the consistent choice in every case.
      ! Zero new constants: the printed 2.5, the model's own profile and its own
      ! extinction coefficient, the divisor being a model STATE and not a
      ! coefficient.
      ! WELL-DEFINED BEFORE THE PROFILE EXISTS: with a vanishing leaf area index,
      ! leaf nitrogen or extinction coefficient - the same three guards the canopy
      ! profile branch itself uses - or a ratio below 1, which a monotonically
      ! decreasing profile cannot produce and which would therefore signal a
      ! degenerate state, the capacity FALLS BACK to the mode-0 value and the day is
      ! COUNTED rather than passed over.
      ! ================================================================ !
      IF (lamina_capacity_profile_basis >= 1 .AND. lamina_n_storage >= 1) THEN
        lcp_factor_day = 0.0_wp
        IF (lai_growth_diag > 1.0e-6_wp .AND. N_leaf > 1.0e-9_wp .AND. kN_canopy > 1.0e-6_wp) &
          lcp_factor_day = canopy_n0_top(N_leaf, lai_growth_diag, kN_canopy) &
                           * lai_growth_diag / N_leaf
        IF (lcp_factor_day >= 1.0_wp) THEN
          lns_sln_ceil_used = sln_max_N_lnh / lcp_factor_day
          lcp_factor_sum    = lcp_factor_sum + lcp_factor_day
          lcp_days_profile  = lcp_days_profile + 1
        ELSE
          lcp_factor_day    = 0.0_wp
          lcp_days_fallback = lcp_days_fallback + 1
        END IF
        lcp_cap_day = lns_sln_ceil_used
      END IF
      ! Record 157: mode 2 defers this whole tier to AFTER the stem and chaff luxury
      ! tier, the alternative order; modes 0 and 1 run it here, in the incumbent
      ! position.
      IF (lamina_n_storage /= 2) THEN
      IF (leafN_drawdown == 1) THEN
        cap = MAX(0.0_wp, MIN(lns_sln_ceil_used * (lai_growth_diag + n5_dlai_prev), &
                              sln_max_N * lai_growth_diag) - N_leaf)
      ELSE
        cap = MAX(0.0_wp, lns_sln_ceil_used * (lai_growth_diag + n5_dlai_prev) - N_leaf)
      END IF
      take = MIN(pool, cap); N_leaf = N_leaf + take
      N_alloc_leaf_day = take; pool = pool - take
      ! rec 157 ON-only ledger (no flux, no state): what the tier was
      ! offered, what it stored, and whether the capacity bound it today.
      IF (lamina_n_storage >= 1) THEN
        lns_cap_day    = cap
        lns_stored_day = take
        lns_stored_cum = lns_stored_cum + take
        lns_bind_day   = (cap > 0.0_wp .AND. take >= cap - 1.0e-12_wp)
        lns_days_live  = lns_days_live + 1
        IF (lns_bind_day) lns_days_bind = lns_days_bind + 1
        IF (take > 0.0_wp) lns_days_store = lns_days_store + 1
      END IF
      END IF
      ! ------------------------------------------------------------------ !
      ! CHAFF STRUCTURAL NITROGEN.
      ! The chaff - glumes, lemmas, paleas and rachis, i.e. spike structural tissue -
      ! must receive the nitrogen its structural DEMAND is charged for. That demand
      ! is charged here, against the pool, at the SAME concentrations as the stem and
      ! on the stem-plus-chaff basis. Depositing all of it into the stem pool instead
      ! leaves the chaff with exactly zero nitrogen, inflates the stem concentration
      ! above the 1.5 %DM ceiling it is supposed to obey, because the chaff mass sits
      ! in the cap basis, and makes the chaff nitrogen fully remobilizable to the
      ! grain through the stem pool.
      ! EQUATION AS CODED, per tier, organ-resolved, pro rata by unmet demand so that
      ! no organ is privileged by statement order:
      !   d_stem  = MAX(0, c W_stem  - N_stem)
      !   d_chaff = MAX(0, c W_chaff - N_chaff)
      !   take    = MIN(pool, d_stem + d_chaff)
      !   N_stem    += take d_stem /(d_stem + d_chaff)
      !   N_chaff   += take d_chaff/(d_stem + d_chaff)
      !   N_storage += take d_chaff/(d_stem + d_chaff)
      ! applied twice, with c the structural concentration (tier 2) and then the
      ! luxury ceiling (tier 4). The TIER TOTAL is unchanged whenever the pool is not
      ! limiting; only the receiving pool changes. The chaff nitrogen is by
      ! declaration the chaff-tagged SUBSET of the reproductive nitrogen, so the same
      ! amount is added to BOTH, which preserves the grain nitrogen as the
      ! reproductive pool minus the chaff share and keeps the daily nitrogen closure
      ! at machine zero.
      ! CONSTANTS: none new. The stem structural minimum 0.005 and luxury ceiling
      ! 0.015 g N g-1 DM are the existing values of Jamieson & Semenov (2000) Field
      ! Crops Res 68:21-29 p.23: a stem minimum requirement of 0.5 % N while the
      ! plant is actively growing, and a concentration of up to 1.5 % N to store
      ! excess nitrogen; the full sources are at their SHARED.f90 declarations. The
      ! chaff IS spike structural tissue and therefore obeys the stem structural-
      ! nitrogen law unchanged; no chaff-specific concentration is introduced and no
      ! observed chaff concentration is used as a constant. The zeros are positivity
      ! guards and the positive-total test is a division guard.
      ! ------------------------------------------------------------------ !
      BLOCK
        REAL(wp) :: d_stem, d_chaff, d_tot, to_stem, to_chaff
        d_stem  = MAX(0.0_wp, n5_stem_nstruct * W_stem  - N_stem)
        d_chaff = MAX(0.0_wp, n5_stem_nstruct * W_chaff - N_chaff)
        d_tot   = d_stem + d_chaff
        take    = MIN(pool, d_tot)
        IF (d_tot > 0.0_wp) THEN
          to_stem  = take * d_stem  / d_tot
          to_chaff = take - to_stem
          N_stem    = N_stem    + to_stem
          N_chaff   = N_chaff   + to_chaff
          N_storage = N_storage + to_chaff
          N_alloc_stem_day    = N_alloc_stem_day    + to_stem
          N_alloc_storage_day = N_alloc_storage_day + to_chaff
          pool = pool - take
        END IF
      END BLOCK
      cap = MAX(0.0_wp, phen_params%init_rootN_frac * W_root - N_root)
      take = MIN(pool, cap); N_root = N_root + take
      N_alloc_root_day = take; pool = pool - take
      ! luxury tier — same law, c = n4_stem_nmax (header above).
      BLOCK
        REAL(wp) :: d_stem, d_chaff, d_tot, to_stem, to_chaff
        d_stem  = MAX(0.0_wp, n4_stem_nmax * W_stem  - N_stem)
        ! Record 63: the CHAFF-SPECIFIC CONCENTRATION CEILING. The chaff otherwise
        ! inherits the true-stem luxury ceiling 0.015 g N g-1 DM, which is a TRUE STEM
        ! statement of Jamieson & Semenov (2000) and sits BELOW the chaff's own measured
        ! concentration, so no allocation could ever reach it. Under this record the
        ! chaff gets its own ceiling 0.017 g N g-1 DM, derived from Bertheloot, Martre
        ! & Andrieu (2008) Plant Physiol Fig. 1D, the chaff nitrogen per unit ground at
        ! anthesis, over Table I, the chaff dry matter at anthesis, for both cultivars;
        ! the full derivation and its status are at the SHARED.f90 declaration. It is
        ! not from any observation of this campaign. The stem keeps its own ceiling
        ! unchanged, and with the record off the inherited stem ceiling is used.
        IF (grainN_partition == 1) THEN
          d_chaff = MAX(0.0_wp, n2_chaff_nmax * W_chaff - N_chaff)
        ELSE
          d_chaff = MAX(0.0_wp, n4_stem_nmax * W_chaff - N_chaff)
        END IF
        d_tot   = d_stem + d_chaff
        take    = MIN(pool, d_tot)
        IF (d_tot > 0.0_wp) THEN
          to_stem  = take * d_stem  / d_tot
          to_chaff = take - to_stem
          N_stem    = N_stem    + to_stem
          N_chaff   = N_chaff   + to_chaff
          N_storage = N_storage + to_chaff
          N_alloc_stem_day    = N_alloc_stem_day    + to_stem
          N_alloc_storage_day = N_alloc_storage_day + to_chaff
          pool = pool - take
        END IF
      END BLOCK
      ! rec 121: the structural payments are allocation too (diagnostic columns)
      IF (n_payment_gate >= 1 .AND. npg_live_today) THEN
        N_alloc_leaf_day    = N_alloc_leaf_day    + npg_paid_L
        N_alloc_stem_day    = N_alloc_stem_day    + npg_paid_S
        N_alloc_storage_day = N_alloc_storage_day + npg_paid_chaff
        N_alloc_root_day    = N_alloc_root_day    + npg_paid_R
      END IF
      ! P4: overflow split - the mobile-origin share returns to the mobile
      ! pool; only assimilation-origin overflow goes to the shoot NO3 store.
      BLOCK
        REAL(wp) :: back_to_mobile
        back_to_mobile = MIN(pool, mobile_in)
        N_mobile_veg = N_mobile_veg + back_to_mobile
        N_NO3_shoot = N_NO3_shoot + (pool - back_to_mobile)
      END BLOCK
    END BLOCK
    N_assim_avail_pa = 0.0_wp
  ELSE
    ! Post-anthesis (Martre et al. 2003): today's assimilated nitrogen is NOT pushed
    ! to the pools by fixed fractions. It becomes the post-anthesis SOURCE for the
    ! grain-nitrogen source-sink below, and whatever the grain does not import
    ! sustains the vegetative nitrogen there.
    N_assim_avail_pa    = N_gain_total
    N_alloc_storage_day = 0.0_wp
    N_alloc_leaf_day    = 0.0_wp
    N_alloc_stem_day    = 0.0_wp
    N_alloc_root_day    = 0.0_wp
    ! rec 157: the pre-anthesis ladder did not run today, so its
    ! day ledger is zeroed rather than left stale in the ON-only CSV.
    IF (lamina_n_storage >= 1) THEN
      lns_cap_day = 0.0_wp ; lns_stored_day = 0.0_wp ; lns_bind_day = .FALSE.
    END IF
  END IF
  IF (npg_open) THEN
    WRITE (unit_npg,'(I0,",",G14.6,",",I0,",",I0,40(",",G14.6))',ADVANCE='NO') iD, dvs2_growth, n_payment_gate, &
      MERGE(1, 0, npg_live_today), npg_c_L, npg_c_S, npg_c_R, &
      npg_pool_mobile_in, npg_pool_assim_in, npg_lux_avail, npg_lux_used, &
      npg_D_L, npg_D_S, npg_D_R, npg_p_L, npg_p_S, npg_p_R, npg_paid_L, npg_paid_S, npg_paid_chaff, npg_paid_R, &
      npg_pool_out, npg_mobile_left, npg_dw_leaf_pot, npg_dw_leaf, npg_dw_stem_pot, npg_dw_stem, &
      npg_dw_stor_pot, npg_dw_stor, npg_dw_root_pot, npg_dw_root, npg_unpaid_ch2o_day, npg_dw_wsc_day, &
      npg_resid_day, npg_lamina_unpaid_day, npg_lamina_unpaid_cum, N_leaf, N_stem, N_root, W_reserve, &
      npg_pool_root_in, npg_root_export, npg_paid_R_root, npg_paid_R_shoot
    ! rec 141: three ON-only extra columns on the SAME row, SAME
    ! non-advancing-WRITE technique, OFF.
    IF (lamina_n_remob_pre >= 1) WRITE (unit_npg,'(3(",",G14.6))',ADVANCE='NO') &
      lnr_rem_avail_day, lnr_rem_used_day, lnr_rem_used_cum
    ! rec 143: eight ON-only extra columns on the SAME row, SAME
    ! non-advancing-WRITE technique, OFF.
    IF (n_status_construction_price >= 1) WRITE (unit_npg,'(8(",",G14.6))',ADVANCE='NO') &
      nap_nni_veg, nsc_price_L, nsc_price_S, nsc_price_chaff, nsc_price_R, &
      nsc_cut_L_day, nsc_cut_S_day, nsc_cut_R_day
    WRITE (unit_npg,'(A)') ''
  END IF
  IF (lna_open) THEN
    WRITE (unit_lna,'(I0,",",G14.6,",",I0,",",I0,13(",",G14.6))',ADVANCE='NO') iD, dvs2_growth, lai_n_arrest, &
      MERGE(1, 0, lna_live_today), lna_pL_day, lna_inc_area_day, lna_removed_day, lna_removed_cum, lna_removed_dm_day, &
      lna_gai_before, lna_gai_after, N_leaf, lna_sln_green_in, lna_shed_need_day, lna_shed_day, lna_shed_cum, lna_sl_bump(0)
    ! rec 144: three ON-only extra columns on the SAME row, SAME
    ! non-advancing-WRITE technique, OFF.
    IF (carbon_coupled_leaf_expansion >= 1) WRITE (unit_lna,'(3(",",G14.6))',ADVANCE='NO') &
      cbc_surplus_day, cbc_f_mean_expanding, cbc_area_bonus_day
    ! rec 179: seven ON-only extra columns, OFF.
    IF (lamina_expansion_n_limit >= 1) WRITE (unit_lna,'(7(",",G14.6))',ADVANCE='NO') &
      lnl_Nph_day, lnl_slnph_day, lnl_maxNpha, lnl_minNpha_day, &
      lnl_a_sup_day, lnl_fN_day, lnl_refused_day
    WRITE (unit_lna,'(A)') ''
  END IF
  IF (rfl_open) THEN
    IF (rfl_day_zmax == 0 .AND. rfl_z_new >= rfl_z_max - 1.0e-9_wp .AND. rfl_z_max > 0.0_wp) rfl_day_zmax = iD
    WRITE (unit_rfl,'(I0,",",G14.6,",",I0,4(",",G14.6))') iD, dvs2_growth, rfl_days_emerged, rfl_rz, rfl_z_legacy, &
      rfl_z_new, rfl_z_max
  END IF
  ! rec 139: one row per soil layer, ON-only.
  IF (nud_open) THEN
    DO i_nud = 1, NSoil_layers
      ! rec 155: ON-only, the six appended whole-plant columns.
      IF (n_uptake_regulation >= 1) THEN
        WRITE (unit_nud,'(I0,",",I0,11(",",G14.6))') iD, i_nud, nud_mf_lay(i_nud), nud_diff_lay(i_nud), &
          nud_act_lay(i_nud), nud_cap_lay(i_nud), nud_real_lay(i_nud), &
          nur_EN_day, nur_Nmob_day, nur_cap_nmax_day, nur_cap_nc_day, nur_dem_before, nur_dem_after
      ELSE
        WRITE (unit_nud,'(I0,",",I0,5(",",G14.6))') iD, i_nud, nud_mf_lay(i_nud), nud_diff_lay(i_nud), &
          nud_act_lay(i_nud), nud_cap_lay(i_nud), nud_real_lay(i_nud)
      END IF
    END DO
  END IF
  ! Record 139 mode 2 and above, on only: one row per soil layer. The delivery
  ! column is the mass-flow plus diffusion terms, both already summed over ions;
  ! the delivery potential is the season-day transport-delivery potential summed
  ! over layers; the active potential is the existing UNBOUNDED kinetic-capacity
  ! diagnostic, unchanged by mode 2; and the bound fraction is the realized
  ! active uptake today over the unbounded active capacity today, which is 1.0 on
  ! a day the delivery never binds and below 1.0 when it does.
  IF (nub_open) THEN
    DO i_nud = 1, NSoil_layers
      WRITE (unit_nub,'(I0,",",I0,",",G14.6,",",I0,3(",",G14.6))') iD, i_nud, &
        nud_mf_lay(i_nud) + nud_diff_lay(i_nud), nud_bound_lay(i_nud), &
        S_massflow_pot_day + S_diffusion_pot_day, S_active_pot_day, &
        MERGE(N_uptake_active_day / S_active_pot_day, 1.0_wp, S_active_pot_day > 1.0e-12_wp)
    END DO
  END IF
  ! rec 139 mode 3: one row per soil layer, ON-only, mode>=3 ONLY.
  IF (nuc_open) THEN
    DO i_nud = 1, NSoil_layers
      WRITE (unit_nuc,'(I0,",",I0,13(",",G14.6),",",I0)') iD, i_nud, &
        theta(i_nud), nud_theta_th, nud_fl3_lay(i_nud), nud_rld3_lay(i_nud), &
        nud_r13_lay(i_nud), nud_geom3_lay(i_nud), nud_cno33_lay(i_nud), nud_cnh43_lay(i_nud), &
        nud_mf_lay(i_nud), nud_diff_lay(i_nud), nud_act_lay(i_nud), &
        nud_mf_lay(i_nud) + nud_diff_lay(i_nud), nud_real_lay(i_nud), nud_bound_lay(i_nud)
    END DO
  END IF
  IF (rsf_open) THEN
    WRITE (unit_rsf,'(I0,8(",",G14.6),",",I0)') iD, dvs2_growth, W_grain, W_stem, W_reserve, rsf_surplus_cum, &
      rsf_label_cum, rsf_resid_day, rsf_resid_cum, rsf_days
  END IF
  IF (ncs_open) THEN
    ! Record 156, on only: the two appended curve-basis columns.
    ! Record 157: the base row is written NON-ADVANCING so that the three
    ! record-157 columns can be appended; at record 157 = 0 the record is closed by
    ! an empty advancing write, which appends nothing.
    IF (n_construction_ceiling_curve >= 1) THEN
      WRITE (unit_ncs,'(I0,",",G14.6,",",I0,",",I0,25(",",G14.6),2(",",I0),2(",",G14.6))',ADVANCE='NO') &
      iD, dvs2_growth, n_construction_status, &
      MERGE(1, 0, npg_live_today), ncs_W_t_ha, ncs_c_crit, ncs_f_crit, &
      ncs_c_floor_L_eff, ncs_c_target_L, ncs_c_target_eff_L, ncs_c_actual_L, ncs_paid1_L, ncs_paid2_L, &
      n5_stem_nstruct, ncs_c_target_S, ncs_c_target_eff_S, ncs_c_actual_S, &
      ncs_paid1_S + ncs_paid1_chaff, ncs_paid2_S + ncs_paid2_chaff, &
      npg_floor_R, ncs_c_target_R, ncs_c_target_eff_R, ncs_c_actual_R, ncs_paid1_R, ncs_paid2_R, &
      npg_p_L, npg_p_S, npg_p_R, ncs_lamina_below_target_cum, ncs_days_fcrit, npg_days_bind_L, &
        ncc_a_used, ncc_b_used
    END IF
    IF (lamina_n_storage == 3) THEN
      WRITE (unit_ncs,'(2(",",G14.6),",",I0,5(",",G14.6))') lns_cap_day, lns_stored_day, &
        MERGE(1, 0, lns_bind_day), lns_nmob_num_day, lns_nmob_t_mobile, lns_nmob_t_no3, &
        lns_nmob_t_stemlux, lns_nmob_t_chaflux
    END IF
  END IF
  IF (cfo_open) THEN
    WRITE (unit_cfo,'(I0,9(",",G14.6),",",I0)') iD, dvs2_growth, W_grain, W_stem, W_reserve, &
      cfo_surplus_day, cfo_wsc_day, cfo_loss_day, cfo_loss_cum, cfo_wsc_cum, MERGE(1, 0, cfo_surplus_day > 0.0_wp)
  END IF
  IF (gnm_open) THEN
    WRITE (unit_gnm,'(I0,11(",",G14.6),",",I0)',ADVANCE='NO') iD, dvs2_growth, W_grain, gnm_N_grain_in, &
      100.0_wp * gnm_N_grain_in / MAX(1.0e-9_wp, W_grain), &
      gnm_room_dm_118, gnm_room_dm_125, MIN(gnm_room_dm_floor, 1.0e30_wp), &
      gnm_to_grain_before, gnm_to_grain_after, gnm_refused_dm_day, gnm_refused_dm_cum, MERGE(1, 0, gnm_floor_bind_day)
    ! rec 130 mode 3: 7 extra ON-only columns, SAME non-advancing-
    ! WRITE technique; modes 0/1/2 (nothing below executes).
    IF (grain_n_min_conc == 3) WRITE (unit_gnm,'(6(",",G14.6),",",I0)',ADVANCE='NO') &
      gnm_Tt_degCd, 1.0e6_wp * gnm_Nstru_pergrain, 100.0_wp * gnm_c_stru_day, &
      MIN(gnm_room_dm_floor, 1.0e30_wp), gnm_Nstru_paid, gnm_Nsto, gnm_phase_day
    WRITE (unit_gnm,'(A)') ''
  END IF
  IF (kpt_open) THEN
    kpt_mgn_k = 0.0_wp
    IF (kernel_number_m2 > 0.0_wp) kpt_mgn_k = 1000.0_wp * MAX(0.0_wp, N_storage - N_chaff) / kernel_number_m2
    WRITE (unit_kpt,'(I0,7(",",G14.6),",",I0,5(",",G14.6),2(",",I0))') iD, dvs2_growth, kernel_number_m2, W_grain, &
      MAX(0.0_wp, N_storage - N_chaff), kpt_mgn_k, kpt_tt, kpt_sup_cum, MERGE(1, 0, kpt_set), kpt_w_pot, gnc_cap, &
      gnc_room_day, gnc_cut_day, gnc_cut_cum, gnc_days_bound, gfb_end_bound
  END IF
  IF (tna_open .AND. tna_step_today) THEN
    WRITE (unit_tna,'(I0,",",ES23.15,2(",",G14.6),",",I0,42(",",G14.6))') iD, tca_taun, tna_pL_used, npg_c_L, tna_days_bind, &
      (ct_dens(i_tc), tna_dN_day(i_tc), tna_DN(i_tc), tna_SN(i_tc), tna_phiN(i_tc), tna_phiC(i_tc), tca_phi(i_tc), &
       i_tc = 0, CT_MAXTIL)
  END IF
  IF (nap_open .AND. tna_step_today) THEN
    WRITE (unit_nap,'(I0,3(",",G14.6),3(",",I0),7(",",G14.6))') iD, dvs2_growth, nap_phi_realized, &
      nap_nni_veg, MERGE(1,0,nap_gate12_open), MERGE(1,0,nap_gate12_open), MERGE(1,0,nap_gate3_open), nap_m, &
      lna_shed_cum, tca_dead, nap_tier3_cap, nap_f_Nstatus_at_gate, nap_N_mobile_veg_at_gate, nap_identity_resid
  END IF
  IF (npl2_open) THEN
    WRITE (unit_npl2,'(I0,",",G14.6,4(",",G14.6),",",I0)') iD, dvs2_growth, ncs_D2_L, ncs_lux2_cap_day, &
      ncs_lux2_used_day, ncs_lux2_used_cum, ncs_lux2_days
  END IF
  IF (nsb_open .AND. tna_step_today) THEN
    WRITE (unit_nsb,'(I0,",",G14.6,5(",",G14.6),",",I0)') iD, dvs2_growth, nap_phi_realized, nap_phi_standing, &
      nap_phi_elig12, nap_phi_elig3, nap_m, n_afford_standing_basis
  END IF
  IF (fnl_open .AND. fnl_step_today) THEN
    WRITE (unit_fnl,'(I0,",",G14.6,",",I0,27(",",G14.6))') iD, dvs2_growth, MERGE(1,0,gnx_win_open), &
      fnl_num, fnl_den, fnl_f_N, &
      (fnl_r(i_tc), fnl_surv_c_only(i_tc), gnx_surv(i_tc), gnx_nfl(i_tc), i_tc = 0, CT_MAXTIL)
  END IF

  ! The grain-nitrogen source-sink (Martre et al. 2003). The grain nitrogen
  ! import is the lesser of the source and the sink: the source is the
  ! post-anthesis assimilation plus the remobilizable vegetative nitrogen, and
  ! the sink is the structural demand during cell division, i.e. sink-limited
  ! then, and unlimited in the storage phase, i.e. source-limited there. No
  ! observed-protein cap is applied anywhere.
  ! The green lamina mass is passed so that the routine can form the same-day
  ! leaf-area identity basis for the capacity bound.
  ! ORDER OF OPERATIONS: the organ-nitrogen turnover runs BEFORE the grain draw
  ! on the same day, because it is what fills the culm-common mobile pool the
  ! grain then draws from - nitrogen distribution within the culm precedes grain
  ! accumulation in the source model's own ordering. The full equations, constants
  ! and sources are at the routine header and at the SHARED.f90 block. With the
  ! record off the routine is not called at all.
  IF (grainN_partition == 1) THEN
    CALL apply_organN_turnover_daily(dvs2_growth, dTT_day_N, W_green, W_stem, W_chaff, &
      W_storage, N_leaf, N_stem, N_chaff, N_storage)
  END IF
  N_leaf_before_lnc = N_leaf                                   ! rec 87 diagnostic only
  ! Record 154: the two root quantities the root limb needs, handed in through
  ! SHARED, since growth_module.f90 has no root mass in scope. The previous-day
  ! root mass was captured before the allocation ran, so this difference is TODAY
  ! root dry-matter increment - the model's OWN existing expression, used
  ! identically at the record-114 mineral call. Set every day: these are
  ! variables nothing else reads.
  ont_W_root_in   = W_root
  ont_dw_root_day = MAX(0.0_wp, W_root - W_root_prev_day)
  CALL allocate_grain_N_source_sink_daily(dvs2_growth, dW_grain_day, N_assim_avail_pa, &
    lai_growth_diag, W_stem, W_green, W_chaff, N_leaf, N_stem, N_storage, &
    grain_import_day, N_from_assim_day, N_remob_leaf_day, N_remob_stem_day, N_surplus_veg_day)
  gnm_dN_grain_prev = grain_import_day   ! rec 130: one-day lag read by TOMORROW's guarded mode-2 ceiling
  ! ==================================================================== !
  ! Record 158: TODAY'S shoot de novo nitrate-reduction DEMAND, captured here -
  ! after the construction ladder and the payment gate have booked the organ
  ! nitrogen payments and after the grain source-sink has resolved how much of
  ! today's grain nitrogen came from NEWLY ASSIMILATED nitrogen - and read by
  ! TOMORROW'S supply block. Pre-anthesis the three organ allocation terms carry it
  ! and the newly assimilated grain term is zero; post-anthesis the allocation
  ! terms are zeroed and that term carries it, so the sum is continuous across
  ! anthesis and never double-counts. The root allocation is EXCLUDED, because
  ! root nitrate is reduced on the photorespiration-independent root limb and not
  ! in the chloroplast. Nothing reads it at mode 0, and it is nevertheless set
  ! unconditionally so that the two arms differ in one place only.
  ! ==================================================================== !
  ! Record 175 at or above 1: THE STATED CONTINUITY, RESTORED. Under record 154
  ! at or above 1 that limb sets the newly assimilated grain nitrogen to zero by
  ! construction and the three organ allocation terms are already zeroed
  ! post-anthesis, so the sum above is IDENTICALLY ZERO after anthesis, the
  ! supply ceiling becomes zero, and the shoot stops reducing nitrate the day
  ! after anthesis - machine-verified on every post-anthesis day at the
  ! Braunschweig deck, with 3.395 g N m-2 stranded in the shoot nitrate store.
  ! The added term is that limb's OWN de novo reduction demand on this day, the
  ! nitrogen it commits to grain protein and to photosynthetic-apparatus
  ! re-synthesis, read from quantities the limb already computes and already
  ! prints in output/daily_organ_n_turnover.csv. The root synthesis is excluded
  ! for the same reason the root allocation is. The term is read ONLY on a day
  ! that limb actually ran, so no stale value can enter, and the enzymatic
  ! ceiling 1.21 umol m-2 s-1 stays outermost and unchanged. Zero new constants.
  ! ==================================================================== !
  gnp_ndb_add_day = 0.0_wp
  IF (grain_n_partition_basis >= 1 .AND. ont_live_today) THEN
    gnp_ndb_add_day = MAX(0.0_wp, ont_grain_import + ont_S_lam + ont_S_stem + ont_S_chaff)
    gnp_ndb_add_cum = gnp_ndb_add_cum + gnp_ndb_add_day
    IF (gnp_ndb_add_day > 0.0_wp) gnp_days_ndb_add = gnp_days_ndb_add + 1
  END IF
  ndb_rate_prev = MAX(0.0_wp, N_alloc_leaf_day) + MAX(0.0_wp, N_alloc_stem_day) &
                + MAX(0.0_wp, N_alloc_storage_day) + MAX(0.0_wp, N_from_assim_day) &
                + gnp_ndb_add_day
  ! Record 87: the leaf area index here is the TOP-OF-DAY value, i.e. exactly the
  ! one this day gas-exchange call was handed, and the leaf nitrogen is the pool
  ! the NEXT day call will read, the organ-clock shed that follows scaling the
  ! leaf nitrogen and the published area by the same fraction. Written only when
  ! record 87 = 1.
  IF (leafn_cap_open) THEN
    WRITE (unit_leafn_cap,'(I0,9(",",G14.6))') iD, lai_growth_diag, W_green, N_leaf_before_lnc, &
      lnc_excess_day, lnc_to_mobile_day, N_leaf, &
      N_leaf / MAX(1.0e-9_wp, lai_growth_diag), &
      100.0_wp * N_leaf / MAX(1.0e-9_wp, W_green), N_mobile_veg
  END IF
  ! Record 154, on only: the daily organ-nitrogen-turnover trace, one row per day
  ! on which the turnover block actually ran (post-anthesis, filling live).
  ! Nothing here is read by any physics path. The closure residual is this
  ! routine's own organ-state nitrogen residual, which carries the root term in
  ! mode 1, printed at full precision so that the conservation statement is
  ! testable from this file alone.
  IF (organ_n_turnover_nema >= 1 .AND. ont_live_today) THEN
    IF (.NOT. ont_open) THEN
      OPEN (NEWUNIT=unit_ont, FILE='output/daily_organ_n_turnover.csv', STATUS='REPLACE', ACTION='WRITE')
      ! Record 167: the base header is written NON-ADVANCING and the six record-167
      ! columns are appended to the SAME record only when the record is at or above 1,
      ! so that the column list has one ordered source for both arms and the base
      ! string is never duplicated.
      WRITE (unit_ont,'(A)',ADVANCE='NO') 'day,dvs2,dTT_degCd,pool_before_gN_m2,uptake_credit_gN_m2,' // &
        'D_lam_gN_m2,D_stem_gN_m2,D_chaff_gN_m2,D_root_gN_m2,sink_pot_gN_m2,grain_import_gN_m2,' // &
        'bound_by,S_lam_pot_gN_m2,S_lam_gN_m2,S_stem_gN_m2,S_chaff_gN_m2,S_root_gN_m2,s_scale,' // &
        'f_N,f_PAR,N_leaf_gN_m2,N_stem_gN_m2,N_chaff_gN_m2,N_grain_gN_m2,N_root_gN_m2,' // &
        'N_mobile_after_gN_m2,organ_closure_resid_gN_m2'
      IF (lamina_synthesis_sheath_basis >= 1) THEN
        WRITE (unit_ont,'(A)',ADVANCE='NO') ',lsb_mode,lsb_M_sheath_gDM_m2,lsb_basis_ratio,' // &
          'lsb_days_applied_cum,lsb_syn_lam_gain_cum_gN_m2,lsb_syn_stem_loss_cum_gN_m2'
      END IF
      WRITE (unit_ont,'(A)') ''
      ont_open = .TRUE.
    END IF
    WRITE (unit_ont,'(I0,10(",",G14.6),",",I0,14(",",G14.6),",",ES23.15)',ADVANCE='NO') iD, dvs2_growth, ont_dTT, &
      ont_pool_before, ont_uptake_credit, ont_D_lam, ont_D_stem, ont_D_chaff, ont_D_root, &
      ont_sink_pot, ont_grain_import, ont_bound_by, ont_S_lam_pot, ont_S_lam, ont_S_stem, &
      ont_S_chaff, ont_S_root, ont_s_scale, ont_f_N, ont_f_PAR, N_leaf, N_stem, N_chaff, &
      MAX(0.0_wp, N_storage - N_chaff), N_root, ont_mob_after, ont_resid
    ! rec 167: TODAY's realized sheath mass and realized basis ratio,
    ! plus the running season ledgers -- ON-only, diagnostic, never read back.
    IF (lamina_synthesis_sheath_basis >= 1) THEN
      WRITE (unit_ont,'(",",I0,2(",",G14.6),",",I0,2(",",G14.6))',ADVANCE='NO') &
        lamina_synthesis_sheath_basis, lsb_m_sheath_day, lsb_ratio_day, &
        lsb_days_applied, lsb_syn_lam_gain_cum, lsb_syn_stem_loss_cum
    END IF
    WRITE (unit_ont,'(A)') ''
  END IF
  ! The grain-nitrogen accrual target (Jamieson & Semenov 2000 p.24, hypothesis
  ! H1: a source-determined thermal-time trajectory to a nitrogen harvest index
  ! of 0.80, linear in the development stage). The unmet part becomes TOMORROW
  ! root uptake demand, an explicit one-day signaling lag, served by the
  ! unchanged supply machinery; the grain draws stay through the source ordering.
  IF (dvs2_growth >= 1.0_wp) THEN
    IF (.NOT. n4_anth_snap) THEN
      n4_ag_anth = N_leaf + N_dead_leaf + N_stem + N_storage
      n4_anth_snap = .TRUE.
    END IF
    ! The accrual trajectory is a GRAIN nitrogen harvest index (Jamieson & Semenov
    ! 2000 p.24, hypothesis H1: grain nitrogen over above-ground nitrogen), so the
    ! achieved side must be GRAIN nitrogen, the reproductive pool minus the chaff
    ! share, and not the whole reproductive pool: with the chaff charged, using the
    ! whole pool would let chaff structural nitrogen, which never enters a kernel,
    ! retire the grain's own nitrogen deficit and silently switch off post-anthesis
    ! uptake. The same substitution is applied at the concentration and protein
    ! places. The anthesis above-ground nitrogen is unchanged and legitimately
    ! includes the chaff. The target 0.80 is Jamieson & Semenov (2000) Field Crops
    ! Res 68:21-29 p.24, hypothesis H1, with its locus at the SHARED.f90
    ! declaration; the 0.0 and 1.0 are a positivity guard and the completion bound
    ! of the trajectory.
    ! Record 86 modes 1 and 4: that trajectory is retired AS A PUMP. It is a TARGET
    ! the crop chases with root uptake - a curve, not a process - and this is the
    ! only place it is still live on the released decks. Under record 86 the demand
    ! is read off the SAME potential the grain's own sink uses, so that the
    ! post-anthesis grain-nitrogen demand is defined ONCE for the whole model:
    !   deficit = MAX(0, potential grain nitrogen(t) - N_grain)
    ! with the potential from Eq. (12) integrated and the grain nitrogen the
    ! reproductive pool minus the chaff share, the same achieved side. The one-day
    ! signaling lag and the supply machinery are unchanged and no new constant
    ! enters. The source is not being misused: on their own data H1 beat the
    ! per-grain sink H2 (coefficient of determination 0.95 against 0.88 and root
    ! mean square deviation 0.18 against 0.19 %N, pp.27-28). What is retired is its
    ! use as a pump; inside a minimum of sink and source the sink is only a
    ! CEILING, which is not what they tested.
    IF (grainN_sink_mech == 1 .OR. grainN_sink_mech == 4 .OR. grainN_sink_mech == 5) THEN
      ! ================================================================ !
      ! THE POST-ANTHESIS GRAIN-NITROGEN UPTAKE DEMAND. Record 97 mode 2 or 3.
      ! The demand below, as a cumulative STOCK deficit against an import that is a
      ! potential RATE, is one step behind the potential curve and therefore creates
      ! a PERMANENT uptake demand the grain is structurally incapable of absorbing,
      ! the difference landing in the mobile vegetative pool every day for the rest
      ! of the season. It also keeps running AFTER grain filling has ENDED: at one
      ! Maricopa deck the cumulative potential and the grain nitrogen both freeze and
      ! the deficit is frozen at a constant 0.6629 g N m-2 d-1 for about 30 d, i.e.
      ! about 60 kg N ha-1 of post-maturity uptake with no sink of any kind.
      ! THE STATEMENT: the demand is the potential RATE of Eq. (12) p.785 - the SAME
      ! quantity the grain's own sink uses, so that the post-anthesis grain-nitrogen
      ! demand is defined once for the whole model - NET OF the standing mobile pool,
      ! because Eq. (11) p.785 is explicit that that pool IS the grain source (the
      ! only source of nitrogen for the grains is the common pool of mobile
      ! nitrogen, read from a 200-dpi page render). The rate is identically zero once
      ! the filling has ended, so the post-maturity uptake closes outright.
      ! Subtracting what the plant already holds is this model's OWN demand idiom: the
      ! leaf maintenance and the luxury reservoir are both written that way, and the
      ! grain term was the only one that ignored the plant's own stock.
      ! The source model has NO root-nitrogen-uptake demand law - its whole-plant
      ! uptake is a forced variable, p.785 - so this is THIS MODEL CONSTRUCTION built
      ! out of that paper's own rate and own pool. No new constant.
      ! Ruled out by measurement rather than by argument: the maximum-dilution cap
      ! NEVER BINDS on those days (the daily demand equals the post-anthesis grain
      ! demand on every post-anthesis day), so adding the mobile pool to the shoot
      ! nitrogen would not fix this. That basis omission is left open and stated.
      ! ================================================================ !
      ! Record 155 LIMB (c): under record 155 mode 1 or 3, AND ONLY when record 154
      ! is at or above 1 and holds the culm-common mobile pool, the post-anthesis
      ! uptake demand takes the SAME statement record 97 modes 2 and 3 code on the
      ! line below - today's potential RATE net of the standing pool - instead of the
      ! cumulative STOCK deficit. No new statement and no new constant: the record
      ! only makes the existing one reachable without moving the record-97 mode on
      ! the deck, the released decks running record 97 = 1.
      IF (itemF_Nclose_mode == 2 .OR. itemF_Nclose_mode == 3 .OR. &
          (n_uptake_regulation == 1 .AND. organ_n_turnover_nema >= 1)) THEN
        IF (itemF_Nclose_mode /= 2 .AND. itemF_Nclose_mode /= 3) nur_days_c = nur_days_c + 1
        n4_grain_deficit_prev = MAX(0.0_wp, gnf_sink_day - MAX(0.0_wp, N_mobile_veg))
      ELSE
      n4_grain_deficit_prev = MAX(0.0_wp, gnf_pot_cum - MAX(0.0_wp, N_storage - N_chaff))
      END IF
    ELSE
    n4_grain_deficit_prev = MAX(0.0_wp, n4_nhi_target * n4_ag_anth * &
        MIN(1.0_wp, dvs2_growth - 1.0_wp) - MAX(0.0_wp, N_storage - N_chaff))
    END IF
  ELSE
    n4_grain_deficit_prev = 0.0_wp
  END IF
  ! Record 86 at or above 1 and not 3, on only: everything needed to test the
  ! grain-nitrogen sink against Eqs. (11) and (12) and against the measured
  ! sink-and-source timing of Martre et al. (2003) p.1963 without re-deriving
  ! anything - the thermal clock, the potential the law carries, WHICH of the two
  ! terms of the minimum bound the day, every source in the ladder, and the day
  ! organ-state nitrogen residual.
  ! THE PLACEMENT IS DELIBERATE and is the one thing a reader must check: this
  ! write sits AFTER the grain-nitrogen source-sink routine and after the
  ! uptake-demand statement, so every column is THIS day value. Placed earlier in
  ! the loop it would carry a ONE-DAY LAG.
  ! The row is written only from the first post-anthesis day, the only days the
  ! law covers. Before that, the development stage and the two grain terms are
  ! not yet defined on day 1 and printing them yields denormals that are
  ! unparseable as floating-point numbers. No information is lost: the law's own
  ! state is zero by initialization until anthesis.
  IF (grainN_sink_mech >= 1 .AND. grainN_sink_mech /= 3 .AND. gnf_snap) THEN
    IF (.NOT. gnf_diag_open) THEN
      gnf_diag_open = .TRUE.
      OPEN (NEWUNIT=unit_gnf, FILE='output/daily_grain_n_sink.csv', STATUS='REPLACE', ACTION='WRITE')
      ! A COLUMN SPLIT in the grain-SOURCE block. The old single column carried the
      ! lamina draw, which is the degradation minus the synthesis, i.e. a SIGNED NET
      ! exchange, and it goes negative on 227 day-rows across 28 of 42 decks -
      ! verified on all 227 to equal the Eq. (5) degradation minus the Eq. (6)
      ! synthesis, i.e. the days the lamina re-synthesis out of the common pool
      ! exceeded its own degradation, so the lamina was a net nitrogen IMPORTER. The
      ! physics is right and is unchanged here; only the reporting is split.
      ! Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781-796,
      ! doi:10.1071/FP08064, Eq. (4) p.784 defines the lamina NET nitrogen change as
      ! the Eq. (6) synthesis minus the Eq. (5) degradation, an intrinsically signed
      ! quantity, while the ONE-WAY remobilization flux is Eq. (5) alone, which is
      ! already printed unchanged. The signed column is therefore named as a net
      ! exchange, positive meaning export to the common pool, and two strictly
      ! one-way columns are added beside it. No state and no published row moves.
      WRITE (unit_gnf,'(A)') 'day,dvs2,tt_post_anth_degCd,n_g0_gN_m2,pot_cum_gN_m2,' // &
        'sink_day_gN_m2,supply_day_gN_m2,bound_by,import_day_gN_m2,' // &
        'from_assim_gN_m2,from_mobile_gN_m2,net_lamina_exchange_gN_m2,' // &
        'lamina_export_gN_m2,lamina_import_gN_m2,from_stem_gN_m2,' // &
        'from_chaff_gN_m2,N_grain_gN_m2,N_leaf_gN_m2,N_stem_gN_m2,N_chaff_gN_m2,' // &
        'N_mobile_gN_m2,organ_state_N_resid_gN_m2_d,' // &
        'eqn5_deg_lam_gN_m2,eqn6_syn_lam_gN_m2,eqn6_syn_pot_gN_m2,f_N,f_PAR,pool_scale'
    END IF
    WRITE (unit_gnf,'(I0,6(",",G14.6),",",A,20(",",G14.6))') iD, dvs2_growth, &
      gnf_tt, gnf_ng0, gnf_pot_cum, gnf_sink_day, gnf_supply_day, &
      MERGE('FLOOR ', MERGE('ENDED ', MERGE('SINK  ', 'SOURCE', &
            MERGE(gnf_pool_scale >= 1.0_wp, gnf_sink_day <= gnf_supply_day, gnf_eqn6)), &
            dvs2_growth >= 2.0_wp .OR. canopy_senesced), &
            grain_n_min_conc >= 1 .AND. gnm_floor_bind_day), &
      grain_import_day, N_from_assim_day, N_from_mobile_day, &
      gnf_draw_lam, MAX(0.0_wp, gnf_draw_lam), MAX(0.0_wp, -gnf_draw_lam), &
      gnf_draw_stem, gnf_draw_chaff, &
      MAX(0.0_wp, N_storage - N_chaff), N_leaf, N_stem, N_chaff, N_mobile_veg, &
      gnf_nres_day, gnf_deg_lam_day, gnf_syn_lam_day, gnf_syn_pot_day, &
      gnf_fN_day, gnf_fPAR_day, gnf_pool_scale
  END IF
  ! N_dead_stem in the closure (see N_plant_before).
  N_plant_after = N_leaf + N_dead_leaf + N_dead_stem + N_stem + N_root + N_storage + N_NO3_root + N_NO3_shoot &
                  + N_mobile_veg   ! P4: mobile vegetative pool in the closure
  N_balance_resid = N_plant_after - N_plant_before - N_soil_uptake_day
  IF (pheno_stage(iD) <= 0) THEN
    phen_stage_name = 'PreEmergence'
  ELSE
    phen_stage_name = phen_params%stage_name(MIN(pheno_stage(iD), NUM_PHENO_STAGES))
  END IF
  prev_height_state = 0.0_wp
  IF (iD > 1) prev_height_state = plant_height(iD-1)
  ! rec 108 (SHARED.f90): the SECOND height call of the day takes the same
  ! reserve-excursion-corrected argument; hgt_wres_anth was latched this morning. The released setting rec 108 = 0 keeps the empirical height law.
  W_stem_hgt = W_stem
  IF (height_reserve_basis == 1) W_stem_hgt = W_stem - (W_reserve - hgt_wres_anth)
  CALL update_structure(sla_from_dvs(dvs2_growth), W_green, W_stem_hgt, plant_density, DVS, &
    prev_height_state, lai_growth_diag, height_growth_diag)
  ! Record 109: ADVANCE the internode law by one day, HERE, because this is the
  ! point at which today's leaf-rank chronology (record 64), today's axis population
  ! (record 66) and today's allocation are all final, so the stem mass less the
  ! reserve is TODAY'S OWN structural stem mass. The law reads that increment and
  ! nothing else, and the reserve is NOT an operand of any length, which is why
  ! record 105 cannot move the height by construction. The law, and not the
  ! structure routine, is then the height author for the day.
  CALL inl_advance(W_stem - W_reserve, ct_axes_live, phen_state%current_dtt)
  IF (height_internode_law == 1) height_growth_diag = inl_h_b3 / 100.0_wp
  plant_height(iD) = height_growth_diag
  ! capture today's realized growth increments as
  ! TOMORROW's prospective N demand (explicit 1-day lag).
  n5_dlai_prev    = MAX(0.0_wp, lai_growth_diag - n5_lai_prev)
  n5_lai_prev     = lai_growth_diag
  n5_dwstem_prev  = MAX(0.0_wp, W_stem  - W_stem_prev_day)
  n5_dwchaff_prev = MAX(0.0_wp, W_chaff - W_chaff_prev_day)
  n5_dwroot_prev  = MAX(0.0_wp, W_root  - W_root_prev_day)
  ! The green-lamina sibling of the three increments above, captured by the SAME
  ! rule at the SAME place. Written unconditionally as a variable nothing reads
  ! at record 155 = 0, and read only by limb (a).
  nur_dwgreen_prev = MAX(0.0_wp, W_green - W_green_prev_day)
  ! R6: + Tcan_act_C (committed Tuzet/SPAC canopy temperature; energy_balance_mode
  ! >= 1 only -- 0.0 under mode 0, where the model holds no committed canopy node).
  WRITE (unit_tleaf,'(I0,4(",",G14.6))') iD, T_leaf_day, T_air_day, (T_leaf_day - T_air_day), &
    Tcan_act_day
  WRITE (unit_tsoil,'(I0,1(",",G14.6))') iD, T_soil_day
  WRITE (unit_resist,'(I0,16(",",G14.6))') iD, gsw_mm_day_3, rcs_day, rca_day, rsa_day, rss_day, &
    LEc_day, LEs_day, trans_act, evapo_act, gsw_day_mean, rcs_daylight, rca_daylight, rsa_daylight, &
    rss_daylight, tr_daylight_mm, es_daylight_mm
  ! The reported above-ground biomass. The sum includes the dead stem, sheath and
  ! chaff dry matter, which is harvested straw and is part of the observed
  ! above-ground biomass, and it uses the grain and chaff masses rather than the
  ! reproductive pool.
  ! The reason the basis matters: the allocation adds the storage share to the
  ! reproductive pool, and the yield-component routine then moves the sink-full
  ! SURPLUS into the stem reserve pool and later credits the reserve draw times
  ! the conversion efficiency to the grain. If the surplus did not also decrement
  ! the reproductive pool, the difference between that pool and the grain plus
  ! chaff would be a phantom, and summing the stem pool AND the reproductive pool
  ! would count the surplus twice; the post-emergence carbon closure residual was
  ! then a median 22.2 % of the summed daily gain on the reproductive basis
  ! against a median 1.70 % on this one. That phantom is closed in
  ! growth_module.f90, where the surplus and the reserve draw are booked on the
  ! reproductive pool, so the reproductive pool equals the grain plus the chaff
  ! on every day and the two above-ground bases coincide. The reported basis here
  ! is kept as it is, because it was already correct.
  ! The reported above-ground biomass is WRITE-ONLY: it feeds no state update and
  ! no flux, so this is a reporting basis and the physics is untouched.
  ! NOT CHANGED HERE: the same stem-plus-reproductive sum is the
  ! shoot mass of the critical-nitrogen dilution curve at SIX places - the chain
  ! above-ground biomass -> critical concentration -> maximum-curve stock gap ->
  ! daily demand, and the chain shoot dry matter -> nitrogen nutrition index ->
  ! nitrogen-status factor -> uptake. Those FEED STATE and are not reporting. The
  ! sum is wrong in TWO ways: it ADDS the reproductive pool, which is already
  ! inside the stem pool on the phantom basis, and it OMITS the grain and the
  ! chaff, although the dilution curve for wheat is defined on TOTAL shoot dry
  ! matter INCLUDING the ears. After anthesis both errors are large and they do
  ! not cancel. Correcting it is a physics change and belongs to the
  ! grain-nitrogen sink item, which cannot be built on a shoot mass that omits
  ! the grain it is filling; it is not done here.
  agb_biomass = W_green + W_dead + W_dead_stem + W_stem + W_grain + W_chaff
  ! Record 102, on only: the DERIVED true-stem and sheath split, evaluated AFTER
  ! the day last write to the stem pool (the allocation, the yield components,
  ! the organ-nitrogen turnover and the structure update all precede this line)
  ! and written to its own file. It reads the stem and reserve pools, the
  ! reported above-ground sum, the model plant height and the live axis count,
  ! writes ONLY its own diagnostics, and feeds no flux.
  grain_yield_biomass = W_grain                          ! Task 7: grain only (not grain+chaff)
  IF (lai_from_organs == 1) THEN
    WRITE (unit_growth,'(I0,10(",",G14.6))') iD, lai_growth_diag, height_growth_diag, W_green, W_dead, W_stem, W_root, W_storage, &
      agb_biomass, grain_yield_biomass, lfo_sla_out
  ELSE
    WRITE (unit_growth,'(I0,9(",",G14.6))') iD, lai_growth_diag, height_growth_diag, W_green, W_dead, W_stem, W_root, W_storage, &
      agb_biomass, grain_yield_biomass
  END IF
  BLOCK
    REAL(wp) :: tkw_diag, hi_diag, spike_diag
    tkw_diag = 0.0_wp
    IF (kernel_number_m2 > 1.0e-9_wp) tkw_diag = 1000.0_wp * W_grain / kernel_number_m2  ! g per 1000 kernels
    hi_diag = 0.0_wp
    IF (agb_biomass > 1.0e-9_wp) hi_diag = W_grain / agb_biomass
    spike_diag = W_grain + W_chaff
    WRITE (unit_yield,'(I0,8(",",G14.6))') iD, dvs2_growth, W_grain, W_chaff, spike_diag, &
      kernel_number_m2, tkw_diag, W_grain, hi_diag
  END BLOCK
  ! Task 8: per-mineral diagnostics (mg m-2; grain concentration mg kg-1 at output only).
  DO i_min = 1, NMIN
    BLOCK
      REAL(wp) :: grain_mgkg
      grain_mgkg = 0.0_wp
      IF (W_grain > 1.0e-9_wp) grain_mgkg = 1000.0_wp * M_grain(i_min) / W_grain   ! mg/m2 / (g/m2) * 1000 = mg/kg
      WRITE (unit_minerals,'(I0,",",A,12(",",G14.6))') iD, mineral_names(i_min), &
        M_uptake_mf(i_min), M_uptake_ac(i_min), M_remob_grain(i_min), M_soil_labile(i_min), &
        M_root(i_min), M_leaf(i_min), M_stem(i_min), M_chaff(i_min), M_grain(i_min), &
        M_balance_resid(i_min), grain_mgkg, M_cap_binding(i_min)
    END BLOCK
  END DO
  IF (fnc_open) &
    WRITE (unit_fnc,'(I0,",",I0,12(",",G14.6))') iD, itemF_Nclose_mode, dvs2_growth, &
      lai_growth_diag, N_leaf, sln_max_N * MAX(0.0_wp, lai_growth_diag), &
      N_leaf / MAX(1.0e-9_wp, lai_growth_diag), fnc_syn_clip_day, fnc_syn_clip_cum, &
      N_mobile_veg, fnc_mob_close_day, fnc_mob_close_cum, &
      n4_grain_deficit_prev, gnf_sink_day
  IF (humusN_open) &
    WRITE (unit_humusN,'(I0,",",I0,8(",",G14.6))') iD, soil_n_mech, NHUM_soil, humus_min_day, &
      humus_min_total, humus_fth_day, humus_fh_day, humus_k2pot, T_soil_day, &
      humus_min_total / MAX(1.0e-9_wp, NHUM_soil_init)
  WRITE (unit_organic,'(I0,15(",",G14.6))') iD, res_fast_C, res_fast_N, res_slow_C, res_slow_N, &
    soil_mic_C, soil_mic_N, decomp_C_day, decomp_N_day, respired_C_day, &
    net_min_NH4_day, immob_NH4_day, immob_NO3_day, turnover_N_day, residue_N_in_day, soilN_balance_resid
  WRITE (unit_pheno,'(I0,5(",",I0),",",A,",",I0,6(",",G14.6))') iD, Year(iD), Month(iD), Day(iD), DOY, &
    phen_state%dap, TRIM(phen_stage_name), pheno_stage(iD), phen_state%cum_tdu, phen_state%current_vf, &
    phen_state%current_df, phen_state%current_photoperiod_h, t_dev_pre, t_dev_post
  ! Record 64: the daily rank-chronology diagnostic row. The local phyllochron
  ! column is the chronology's OWN emergent degree-Cd per leaf at today's Haun stage,
  ! the reference value times (HS/2)**0.3 above Haun stage 2 and times the
  ! seed-reserve ratio below it, i.e. the reciprocal of the rate actually
  ! integrated, printed so that the reference scale 73.9 degree-Cd and the
  ! appearance-rate decline are both auditable from the output alone. It is a
  ! print of state and not a state itself.
  IF (clr_diag_open) THEN
    BLOCK
      REAL(wp) :: phyll_local
      IF (clr_HS < clr_fC_HS_ref) THEN
        phyll_local = clr_phyllochron_Cd * clr_seed_ratio
      ELSE
        phyll_local = clr_phyllochron_Cd * (clr_HS / clr_fC_HS_ref) ** (-clr_fC_exp)
      END IF
      WRITE (unit_clr,'(I0,2(",",I0),3(",",G14.6),4(",",I0),",",G14.6,4(",",I0),2(",",G14.6))') &
        iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, clr_HS, clr_HS_ev, &
        clr_n_app, clr_n_lig, clr_n_app_ev, clr_n_lig_ev, phyll_local, &
        clr_FLN, clr_FLN_tip, clr_FLN_ev, clr_FLN_ev_tip, clr_tau_app_flag, clr_tau_lig_flag
    END BLOCK
  END IF
  ! Record 65: the daily per-shoot blade-area diagnostic. The areas are printed
  ! both in cm2 per main stem, the unit every source measurement is in, so that
  ! the comparison can be read directly against the source, and in m2 per shoot,
  ! a pure unit conversion and the basis record 66 will multiply by an axis
  ! density. NO density is applied here.
  IF (cba_diag_open) THEN
    WRITE (unit_cba,'(I0,2(",",I0),2(",",G14.6),3(",",I0),4(",",G14.6),",",I0)') &
      iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, clr_HS, &
      clr_n_app, clr_n_lig, cba_n_expanding, &
      cba_A_shoot_pot, cba_A_shoot, cba_A_shoot_pot * 1.0e-4_wp, cba_A_shoot * 1.0e-4_wp, &
      clr_FLN
  END IF
  ! Record 66: the daily axis-population row. The densities are axes m-2 by class
  ! (index 0 is the main stem) and the two area indices are m2 m-2 of GROUND, the
  ! first ground-area quantity the canopy chain produces.
  IF (ct_diag_open) THEN
    WRITE (unit_ct,'(I0,2(",",I0),14(",",G14.6),",",I0,",",G14.6)') &
      iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, clr_HS, &
      ct_fpar, ct_axes_live, ct_axes_born, &
      MERGE(ct_axes_live / plant_density, 0.0_wp, plant_density > 0.0_wp), &
      ct_dens(0), ct_dens(1), ct_dens(2), ct_dens(3), ct_dens(4), ct_dens(5), &
      ct_bai_pot, ct_bai, ct_cess_why, ct_fs_implied
  END IF
  ! the daily senescence row.
  IF (cls_diag_open) THEN
    WRITE (unit_cls,'(I0,2(",",I0),4(",",G14.6),",",I0,11(",",G14.6))') &
      iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, ct_hs_axis(0, &
      MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)), &
      cls_gl(0), cls_sl(0), cls_phase(0), &
      cls_gl(1), cls_gl(2), cls_gl(3), cls_gl(4), cls_gl(5), &
      cls_gai, cls_dai, cls_dai_pool, cls_dai_culled, &
      cls_dai_pool + cls_dai_culled, cls_gai + cls_dai_pool + cls_dai_culled
  END IF
  ! the daily non-lamina row (two bases, kept apart).
  IF (cse_diag_open) THEN
    WRITE (unit_cse,'(I0,2(",",I0),21(",",G14.6))') &
      iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, ct_hs_axis(0, &
      MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)), cls_sl(0), &
      cse_sai, cse_sdi, cse_sdi_pool, cse_sai_n1, cse_sai_n3, cse_sai_full, &
      cse_iai, cse_eai, cse_eai_lo, cse_eai_hi, cse_eai_meas, &
      cls_gai, cls_dai_pool + cls_dai_culled, cse_nonlam, cse_ratio, &
      cse_paig, cse_pai_li, cse_paig_axis
  END IF
  ! the daily closure and limiter row.
  IF (csl_diag_open) THEN
    WRITE (unit_csl,'(I0,2(",",I0),8(",",G14.6),",",I0,3(",",I0),3(",",G14.6))') &
      iD, phen_state%dap, pheno_stage(iD), &
      csl_lai_sink, csl_pai_diag, csl_lai_Ncap, csl_area_built, csl_dA_day, &
      csl_demand_day, csl_alloc_day, csl_backlog, &
      MERGE(1, 0, csl_limit_active), csl_fln_cut, csl_cohorts_aborted, &
      csl_births_suspended, csl_sla_real, W_green, csl_fsum_err
  END IF
  ! D2 (recs 82/83/84): the daily organ-level row.
  IF (d2_diag_open) THEN
    ! (rec 84 leaf-death coupling): 8 columns appended — the organ law's retired green area and fraction,
    ! the lamina DM / N it carried out (split dead / resorbed), N_leaf and the leaf-N-per-area STATE N_leaf/LAI.
    WRITE (unit_d2,'(I0,2(",",I0),3(",",G14.6),23(",",G14.6),",",I0,",",I0,2(",",G14.6),7(",",G14.6),9(",",G14.6),8(",",G14.6))') &
      iD, phen_state%dap, pheno_stage(iD), clr_tau_Cd, clr_HS, tcg_csur_day, &
      tcg_sup(1:5), tcg_dem(1:5), tcg_ratio(1:5), ct_dens(0:5), &
      ct_axes_live, ct_axes_born, tcg_env600_days, tcg_env800_days, tcg_ne_env_lo, tcg_ne_env_hi, &
      lew_srwc, lew_se_rw, lew_theta_rw, lew_thetafc_rw, lew_fw, lew_fw_200, lew_fw_800, &
      lew_bai, lew_bai_clock, cls_gai, lai_growth_diag, lfo_sla_out, W_green, &
      csl_demand_day, csl_alloc_day, csl_backlog, &
      lfo_gai_died_day, lfo_fdied_day, lfo_sen_dm_day, lfo_sen_N_out_day, lfo_sen_N_dead_day, lfo_sen_N_resorb_day, &
      N_leaf, MERGE(N_leaf / lai_growth_diag, 0.0_wp, lai_growth_diag > 1.0e-9_wp)
  END IF
  ! Record 99 mode 2: one row per stepped window day. The cohort sum and the
  ! identity residual are printed at full double precision, the identities being
  ! gated below the print floor of the other columns.
  IF (tca_diag_open) THEN
    IF (tca_days > tca_days_written) THEN
      tca_days_written = tca_days
      WRITE (unit_tca,'(I0,",",I0,5(",",G14.6),2(",",ES23.15),",",G14.6,",",I0,54(",",G14.6))') &
        iD, phen_state%dap, clr_tau_Cd, tca_taun, tca_gsup_prev, tca_ltot, tca_ftot, tca_sum_s, tca_ident, &
        tca_dead, tca_days_skip, &
        (ct_dens(i_tc), tca_hs(i_tc), tca_L(i_tc), tca_share(i_tc), tca_sax(i_tc), tca_dax(i_tc), &
         tca_S(i_tc), tca_D(i_tc), tca_phi(i_tc), i_tc = 0, CT_MAXTIL)
    END IF
  END IF
  ! Record 99 mode 4: one row per stepped window day, plus the endpoint row
  ! (mortality-shape fraction 1 and the endpoint flag set); its per-cohort share
  ! and ledger columns are those of the last stepped day, and the densities and
  ! the deaths are after the endpoint cull.
  IF (tca4_diag_open) THEN
    IF (tca_days > tca4_days_written .OR. (tca4_endpoint_day .AND. .NOT. tca4_endpoint_written)) THEN
      tca4_days_written = tca_days
      IF (tca4_endpoint_day) tca4_endpoint_written = .TRUE.
      ! taun is printed at full double precision here (the mode-2 file prints it G14.6): F_taun is
      ! recomputed from it by gate I6 to 1e-12, which a 6-digit taun could not support.
      WRITE (unit_tca4,'(I0,",",I0,",",G14.6,",",ES23.15,3(",",G14.6),2(",",ES23.15),",",G14.6,",",I0,' // &
                        '",",ES23.15,3(",",G14.6),",",I0,54(",",G14.6))') &
        iD, phen_state%dap, clr_tau_Cd, tca_taun, tca_gsup_prev, tca_ltot, tca_ftot, tca_sum_s, tca_ident, &
        tca_dead, tca_days_skip, tca4_F, tca4_U, tca4_target, tca_unaff, MERGE(1, 0, tca4_endpoint_day), &
        (ct_dens(i_tc), tca_hs(i_tc), tca_L(i_tc), tca_share(i_tc), tca_sax(i_tc), tca_dax(i_tc), &
         tca_S(i_tc), tca_D(i_tc), tca_phi(i_tc), i_tc = 0, CT_MAXTIL)
    END IF
  END IF
  leafN_conc = 0.0_wp
  stemN_conc = 0.0_wp
  rootN_conc = 0.0_wp
  storageN_conc = 0.0_wp
  IF (W_green > 1.0e-12_wp) leafN_conc = N_leaf / W_green
  IF (W_stem > 1.0e-12_wp) stemN_conc = N_stem / W_stem
  IF (W_root > 1.0e-12_wp) rootN_conc = N_root / W_root
  ! Concentrations on WHOLE-ORGAN oven-dry dry matter, the observation basis: the
  ! grain concentration is the grain nitrogen over the grain dry matter and the
  ! chaff concentration the chaff structural nitrogen over the chaff dry matter.
  ! The spike-pool concentration is retained as a diagnostic only.
  IF (W_storage > 1.0e-12_wp) storageN_conc = N_storage / W_storage
  grainN_conc_b3 = 0.0_wp; chaffN_conc_b3 = 0.0_wp
  IF (W_grain > 1.0e-12_wp) grainN_conc_b3 = MAX(0.0_wp, N_storage - N_chaff) / W_grain
  IF (W_chaff > 1.0e-12_wp) chaffN_conc_b3 = N_chaff / W_chaff
  grain_protein_pct = 0.0_wp
  ! Wheat grain protein from the GRAIN nitrogen concentration, the denominator
  ! being the grain dry matter and not the reproductive pool, and the chaff
  ! structural nitrogen excluded. The factor 5.7 is the standard
  ! nitrogen-to-protein conversion for wheat (Jones 1931 USDA Circular 183, which
  ! is not held on this machine).
  IF (W_grain > 1.0e-12_wp) grain_protein_pct = 100.0_wp * 5.7_wp * MAX(0.0_wp, N_storage - N_chaff) / W_grain
  BLOCK
    REAL(wp) :: nhi_diag_v, n_ag_now
    ! N buried in dead stem/chaff is still above-ground N.
    n_ag_now = N_leaf + N_dead_leaf + N_dead_stem + N_stem + N_storage
    nhi_diag_v = 0.0_wp
    ! The nitrogen harvest index is GRAIN nitrogen over above-ground nitrogen
    ! (Jamieson & Semenov 2000 p.24), so the numerator excludes the chaff
    ! structural share now that it is charged; the same substitution as the accrual
    ! target and the concentration and protein places. The denominator, the
    ! above-ground nitrogen, legitimately keeps it.
    IF (n_ag_now > 1.0e-9_wp) nhi_diag_v = MAX(0.0_wp, N_storage - N_chaff) / n_ag_now
    IF (soilN_layered == 1) THEN
    ! Layered: incumbent columns + U_act_NO3/U_act_NH4 appended at row END.
    WRITE (unit_nitrogen,'(I0,65(",",G14.6))',ADVANCE='NO') iD, N_soil_min, N_fert_day, N_mineral_day, &
      N_leach_d, N_soil_uptake_day, N_assim_biochem_day, N_assim_biochem_retained_day, &
      N_demand_day, f_Tr_N, f_CO2_N, &
      N_remob_leaf_day, N_remob_stem_day, &
      N_alloc_leaf_day, N_alloc_stem_day, N_alloc_root_day, N_alloc_storage_day, &
      N_leaf, N_dead_leaf, N_stem, N_root, N_storage, N_balance_resid, &
      leafN_conc, stemN_conc, rootN_conc, storageN_conc, grain_protein_pct, &
      N_NO3_root, N_NO3_shoot, N_assim_root_day, N_assim_shoot_day, f_photoresp_diag, &
      N_uptake_massflow_day, N_uptake_active_day, &
      grain_import_day, N_from_assim_day, N_surplus_veg_day, &
      N_soil_min, N_NH4_soil, N_urea_soil, N_hydrolysis_day, N_nitrification_day, soilN_balance_resid, &
      N_volat_day, N_denit_day, TAN_surf, &
      veg_N_demand_day, grain_N_demand_postanth_day, &
      N_uptake_pre_cum, N_uptake_post_cum, nhi_diag_v, &
      N_mobile_veg, N_from_mobile_day, &
      (N_leaf + N_dead_leaf) / MAX(1.0e-9_wp, W_green + W_dead), &
      MAX(0.0_wp, N_storage - N_chaff), N_chaff, grainN_conc_b3, chaffN_conc_b3, &
      S_active_pot_day, S_massflow_pot_day, U_pot_day, f_Nstatus_diag, a_root_diag, &
      U_act_NO3_day, U_act_NH4_day
    ! rec 158: five ON-only extra columns on the SAME row, same
    ! non-advancing-WRITE technique, OFF.
    IF (nitrate_drain_rate_basis >= 1) WRITE (unit_nitrogen,'(5(",",G14.6))',ADVANCE='NO') &
      ndb_nmax_store_day, ndb_nmax_rate_day, ndb_rate_prev, ndb_ag_last, ndb_as_last
    ! rec 159: five ON-only columns, same technique.
    WRITE (unit_nitrogen,'(A)') ''
    ! Per-layer state + bottom leach flux (daily_soilN_layers.csv).
    WRITE (unit_soilN_lay,'(I0,*(",",G14.6))') iD, NO3_lay(1:NSoil_layers), NH4_lay(1:NSoil_layers), &
      N_leach_d, q_iface_day_mm(NSoil_layers), N_adv_substeps_day
    ELSE
    WRITE (unit_nitrogen,'(I0,63(",",G14.6))',ADVANCE='NO') iD, N_soil_min, N_fert_day, N_mineral_day, &
      N_leach_d, N_soil_uptake_day, N_assim_biochem_day, N_assim_biochem_retained_day, &
      N_demand_day, f_Tr_N, f_CO2_N, &
      N_remob_leaf_day, N_remob_stem_day, &
      N_alloc_leaf_day, N_alloc_stem_day, N_alloc_root_day, N_alloc_storage_day, &
      N_leaf, N_dead_leaf, N_stem, N_root, N_storage, N_balance_resid, &
      leafN_conc, stemN_conc, rootN_conc, storageN_conc, grain_protein_pct, &
      N_NO3_root, N_NO3_shoot, N_assim_root_day, N_assim_shoot_day, f_photoresp_diag, &
      N_uptake_massflow_day, N_uptake_active_day, &
      grain_import_day, N_from_assim_day, N_surplus_veg_day, &
      N_soil_min, N_NH4_soil, N_urea_soil, N_hydrolysis_day, N_nitrification_day, soilN_balance_resid, &
      N_volat_day, N_denit_day, TAN_surf, &
      veg_N_demand_day, grain_N_demand_postanth_day, &
      N_uptake_pre_cum, N_uptake_post_cum, nhi_diag_v, &
      N_mobile_veg, N_from_mobile_day, &
      (N_leaf + N_dead_leaf) / MAX(1.0e-9_wp, W_green + W_dead), &
      MAX(0.0_wp, N_storage - N_chaff), N_chaff, grainN_conc_b3, chaffN_conc_b3, &
      S_active_pot_day, S_massflow_pot_day, U_pot_day, f_Nstatus_diag, a_root_diag
    ! rec 158: five ON-only extra columns on the SAME row, same technique.
    IF (nitrate_drain_rate_basis >= 1) WRITE (unit_nitrogen,'(5(",",G14.6))',ADVANCE='NO') &
      ndb_nmax_store_day, ndb_nmax_rate_day, ndb_rate_prev, ndb_ag_last, ndb_as_last
    ! rec 159: five ON-only columns, same technique.
    WRITE (unit_nitrogen,'(A)') ''
    END IF
  END BLOCK
  ! Daily carbon/growth diagnostics.
  !
  ! WHAT THE f_leaf/f_stem/f_root/f_storage COLUMNS OF
  ! output/daily_carbon_growth.csv ACTUALLY ARE: the TABLED partition, obtained
  ! from a fresh standalone call to partition_fractions below, purely for this
  ! file. They are NOT the realized partition that moved the carbon. The realized
  ! fractions are formed afterwards inside growth_module's allocate_daily_biomass
  ! and differ by three steps this file does not see:
  !   (a) the stem closure f_stem = MAX(0, 1 - f_leaf - f_root - f_storage),
  !   (b) the renormalization over the non-leaf organs in sink_lamina_split
  !       (with the spike share held under record 89 >= 1), and
  !   (c) the one corrective glucose iteration that re-solves the increment on
  !       the fractions (a) and (b) produced.
  ! The carbon actually allocated uses those realized values (dw_root =
  ! f_root * dm_gain in growth_module), so the four columns here will not in
  ! general reproduce the organ increments in this same file, and the printed
  ! quartet need not sum to one.
  !
  ! Likewise dm_gain_g_m2_d below is an INDEPENDENT diagnostic recomputation
  ! from the daily assimilation and respiration totals, not growth_module's own
  ! dm_gain. The two agree only where the corrective iteration of (c) is inert.
  ! Both are reported as written; nothing here feeds back into the model state.
  BLOCK
    REAL(wp) :: net_ch2o, dm_gain_diag, f_lf, f_st, f_rt, f_sr
    REAL(wp), PARAMETER :: co2_to_ch2o_diag = 30.0_wp / 44.0_wp
    net_ch2o = MAX(0.0_wp, An_growth_day) * co2_to_ch2o_diag
    ! dm_gain = net_ch2o - Rm - Rg (conservation: supply = maintenance + growth_resp + DM)
    dm_gain_diag = MAX(0.0_wp, net_ch2o - rm_growth_day - rg_growth_day)
    ! rec 121(b): rg_growth_day = supply - dm_gain - explicit N-unpaid residual, so the
    ! subtraction above would report the residual as dry matter; remove it (ON only).
    IF (n_payment_gate >= 1) dm_gain_diag = MAX(0.0_wp, dm_gain_diag - npg_resid_day)
    CALL partition_fractions(dvs2_growth, f_lf, f_st, f_rt, f_sr)
    WRITE (unit_carbon_growth,'(I0,",",G14.6,",",A,17(",",G14.6))') &
      iD, dvs2_growth, TRIM(phen_stage_name), An_growth_day, A_day, &
      net_ch2o, rm_growth_day, rg_growth_day, dm_gain_diag, &
      f_lf, f_st, f_rt, f_sr, &
      W_green, W_dead, W_stem, W_root, W_storage, lai_growth_diag, height_growth_diag
  END BLOCK
  ! Daily radiation/albedo/aerodynamics diagnostics.
  WRITE (unit_rad_aero,'(I0,13(",",G14.6))') iD, alpha_mix_s, lai_growth_diag, height_growth_diag, &
    Rn_day, Ac_day, As_day, &
    frontal_area_index_s, disp_height_s, rough_len_s, &
    rca_day, rca_daylight, rsa_daylight
  ! Daily gas-exchange diagnostics. The conversions are those of the output
  ! specification; the measurement window from flag leaf to maturity is applied
  ! in the analysis scripts and not here.
  BLOCK
    REAL(wp) :: daylight_s_ge, An_ground_umol, An_leaf_umol, Tr_mmol_s
    daylight_s_ge = MAX(1.0_wp, Daylength * 3600.0_wp)
    ! An_day [g CO2 m-2 d-1] -> umol CO2 m-2 s-1 averaged over daylight
    An_ground_umol = (An_day / 44.0e-6_wp) / daylight_s_ge
    An_leaf_umol   = An_ground_umol / MAX(1.0e-6_wp, lai_growth_diag)
    ! tr_daylight_mm [kg H2O m-2 daylight-1] -> mmol H2O m-2 s-1.
    Tr_mmol_s = (tr_daylight_mm / daylight_s_ge) / 0.01801528_wp * 1000.0_wp
    IF (ft_pot_unstressed == 1) THEN
      WRITE (unit_gas_exch,'(I0,16(",",G14.6))') iD, An_day, A_day, N_day, &
        An_ground_umol, An_leaf_umol, &
        gsw_day_mean, tr_daylight_mm, Tr_mmol_s, &
        Ci_day, VPD_day, T_leaf_day, T_air_day, lai_growth_diag, height_growth_diag, &
        ci_daylight_mean_diag, Tcan_act_day
    ELSE
      WRITE (unit_gas_exch,'(I0,15(",",G14.6))') iD, An_day, A_day, N_day, &
        An_ground_umol, An_leaf_umol, &
        gsw_day_mean, tr_daylight_mm, Tr_mmol_s, &
        Ci_day, VPD_day, T_leaf_day, T_air_day, lai_growth_diag, height_growth_diag, &
        Tcan_act_day
    END IF
  END BLOCK
  T_leaf_mean_total = T_leaf_mean_total + T_leaf_day; T_air_mean_total = T_air_mean_total + T_air_day
  Ci_mean_total = Ci_mean_total + Ci_day
  total_precip = total_precip + WS_d_s
  total_irrig = total_irrig + Irr_d_s
  total_supply_eff = total_supply_eff + pr_eff * 1000.0_wp
  total_intercept = total_intercept + MAX(0.0_wp, (pr + irr - pr_eff)) * 1000.0_wp
  total_et_out = total_et_out + ET_act
  total_soil_evap = total_soil_evap + evapo_act
  total_trans_act = total_trans_act + trans_act
  total_drain  = total_drain + drain_mm_d
  total_runoff = total_runoff + runoff_mm_d
  total_hydro_mb_resid = total_hydro_mb_resid + hydro_mb_resid_mm
  total_hydro_mb_correction = total_hydro_mb_correction + hydro_mb_correction_mm
  total_hydro_top_flux_m = total_hydro_top_flux_m + hydro_top_flux_m_day
  total_hydro_root_sink_m = total_hydro_root_sink_m + hydro_root_sink_m_day
  total_hydro_bottom_flux_m = total_hydro_bottom_flux_m + hydro_bottom_flux_m_day
  hydro_ledger%top_flux_req_mm = hydro_ledger%top_flux_req_mm + hydro_top_req_m_day * 1000.0_wp
  hydro_ledger%top_flux_act_mm = hydro_ledger%top_flux_act_mm + hydro_top_flux_m_day * 1000.0_wp
  hydro_ledger%infiltration_mm = hydro_ledger%infiltration_mm + hydro_top_infil_m_day * 1000.0_wp
  hydro_ledger%evaporation_mm = hydro_ledger%evaporation_mm + soil_evap_mm_d
  hydro_ledger%pond_evap_mm = hydro_ledger%pond_evap_mm + pond_evap_mm_d
  hydro_ledger%runoff_mm = hydro_ledger%runoff_mm + runoff_mm_d
  hydro_ledger%root_uptake_pot_mm = hydro_ledger%root_uptake_pot_mm + trans_act_pot
  hydro_ledger%root_uptake_act_mm = hydro_ledger%root_uptake_act_mm + hydro_root_sink_m_day * 1000.0_wp
  hydro_ledger%bottom_flux_mm = hydro_ledger%bottom_flux_mm + hydro_bottom_flux_m_day * 1000.0_wp
  hydro_ledger%storage_mm = storage_mm
  WRITE (unit_hydroledger,'(I0,19(",",G14.6))') iD, &
    hydro_top_req_m_day * 1000.0_wp, hydro_top_flux_m_day * 1000.0_wp, hydro_top_infil_m_day * 1000.0_wp, &
    soil_evap_mm_d, pond_evap_mm_d, runoff_mm_d, trans_act_pot, hydro_root_sink_m_day * 1000.0_wp, &
    hydro_bottom_flux_m_day * 1000.0_wp, storage_mm, hydro_ledger%top_flux_req_mm, hydro_ledger%top_flux_act_mm, &
    hydro_ledger%infiltration_mm, hydro_ledger%evaporation_mm, hydro_ledger%pond_evap_mm, hydro_ledger%runoff_mm, &
    hydro_ledger%root_uptake_pot_mm, hydro_ledger%root_uptake_act_mm, hydro_ledger%bottom_flux_mm
  total_balance_soil_mm = total_balance_soil_mm + balance_soil_mm
  total_balance_surface_mm = total_balance_surface_mm + balance_surface_mm
  total_balance_mm = total_balance_mm + balance_mm
  total_soil_mb_resid_mm = total_soil_mb_resid_mm + soil_mb_resid_mm

  Total_trans = Total_trans + trans_act; Total_assim = Total_assim + An_day
  Total_N_ass = Total_N_ass + N_assim_biochem_day; Total_gsw_1 = Total_gsw_1 + gsw_mm_day_3
  ! The total biomass includes the dead stem, sheath and chaff dry matter, so that
  ! the total is conserved, and the roots. It feeds the reported biomass and
  ! nitrogen-to-biomass columns only.
  Biomass = W_green + W_dead + W_dead_stem + W_stem + W_root + W_grain + W_chaff

  A_data (iD,1)     = T_leaf_day
  A_data (iD,2)     = T_air_day
  A_data (iD,3)     = H_day
  A_data (iD,4)     = LE_day
  A_data (iD,5)     = Rn_day

  flux_daily (iD,1) = Rn_day
  flux_daily (iD,2) = LE_day
  flux_daily (iD,3) = H_day
  flux_daily (iD,4) = G_day
  flux_daily (iD,5) = trans_act
  flux_daily (iD,6) = evapo_act
  flux_daily (iD,7) = Ir_mm_d
  flux_daily (iD,8) = drain_mm_d
  flux_daily (iD,9) = runoff_mm_d

  cum_tr = cum_tr + trans_act
  cum_e = cum_e + evapo_act
  cum_ir = cum_ir + Ir_mm_d
  cum_et = cum_et + (trans_act + evapo_act)
  cum_et_ir = cum_et_ir + (trans_act + evapo_act + Ir_mm_d)
  cum_precip = cum_precip + WS_d_s

  flux_daily_full (iD,1)  = Rn_day
  flux_daily_full (iD,2)  = G_day
  flux_daily_full (iD,3)  = LE_day
  flux_daily_full (iD,4)  = LEc_day
  flux_daily_full (iD,5)  = LEs_day
  flux_daily_full (iD,6)  = H_day
  flux_daily_full (iD,7)  = Hc_day
  flux_daily_full (iD,8)  = Hs_day
  flux_daily_full (iD,9)  = Rn_day
  flux_daily_full (iD,10) = trans_act
  flux_daily_full (iD,11) = evapo_act
  flux_daily_full (iD,12) = Ir_mm_d
  flux_daily_full (iD,13) = trans_act + evapo_act
  flux_daily_full (iD,14) = ET_act
  flux_daily_full (iD,15) = cum_tr
  flux_daily_full (iD,16) = cum_e
  flux_daily_full (iD,17) = cum_ir
  flux_daily_full (iD,18) = cum_et
  flux_daily_full (iD,19) = cum_et_ir
  flux_daily_full (iD,20) = WS_d_s
  flux_daily_full (iD,21) = cum_precip
  flux_daily_full (iD,22) = Irr_d_s
  flux_daily_full (iD,23) = total_irrig
  flux_daily_full (iD,24) = T_leaf_day
  flux_daily_full (iD,25) = T_air_day
  flux_daily_full (iD,26) = T_leaf_day - T_air_day
  flux_daily_full (iD,27) = dTotrad_s
  flux_daily_full (iD,28) = Tmax_d_s
  flux_daily_full (iD,29) = Tmin_d_s
  flux_daily_full (iD,30) = temp_s
  flux_daily_full (iD,31) = Ac_day
  flux_daily_full (iD,32) = As_day
  flux_daily_full (iD,33) = storage_mm
  flux_daily_full (iD,34) = delta_storage
  flux_daily_full (iD,35) = balance_mm
  flux_daily_full (iD,36) = balance_no_ir
  flux_daily_full (iD,37) = Biomass
  flux_daily_full (iD,38) = Total_N_ass
  IF (Biomass > 0.0_wp) THEN
    flux_daily_full (iD,39) = Total_N_ass / Biomass
  ELSE
    flux_daily_full (iD,39) = 0.0_wp
  END IF
  flux_daily_full (iD,40) = runoff_mm_d
  flux_daily_full (iD,41) = soil_mb_resid_mm
  flux_daily_full (iD,42) = ponding_mm
  flux_daily_full (iD,43) = hydro_mb_resid_mm
  flux_daily_full (iD,44) = trans_act_pot + evapo_act_pot
  flux_daily_full (iD,45) = pond_evap_mm_d
  flux_daily_full (iD,46) = soil_evap_mm_d
  flux_daily_full (iD,47) = infil_mm_d
  flux_daily_full (iD,48) = hydro_top_infil_m_day * 1000.0_wp
  flux_daily_full (iD,49) = infil_from_pond_mm_d
  flux_daily_full (iD,50) = balance_soil_mm
  flux_daily_full (iD,51) = balance_surface_mm
  flux_daily_full (iD,52) = top_flux_net_mm_d
  flux_daily_full (iD,53) = top_flux_up_mm_d
  flux_daily_full (iD,54) = hydro_top_flux_m_day * 1000.0_wp
  flux_daily_full (iD,55) = hydro_top_infil_m_day * 1000.0_wp
  flux_daily_full (iD,56) = hydro_root_sink_m_day * 1000.0_wp
  flux_daily_full (iD,57) = hydro_bottom_flux_m_day * 1000.0_wp
  flux_daily_full (iD,58) = hydro_mb_resid_day_m * 1000.0_wp
  flux_daily_full (iD,59) = delta_storage - &
                            (hydro_top_flux_m_day * 1000.0_wp - hydro_root_sink_m_day * 1000.0_wp - &
                             hydro_bottom_flux_m_day * 1000.0_wp)

  swc_daily (iD,1)  = WS_d_s
  swc_daily (iD,2)  = SUM(theta*tk)*1000.0_wp
  swc_daily (iD,3)  = vwc_root*root_dep*1000.0_wp
  swc_daily (iD,4)  = avail_water_mm

  DO I = 1, NSoil_layers
    swc_layers_raw_daily(iD,I) = theta_raw(I)*tk(I)*1000.0_wp
    swc_layers_daily(iD,I) = theta(I)*tk(I)*1000.0_wp
    swc_layers_mean_daily(iD,I) = (theta_day_sum(I) / t_day_sum) * tk(I) * 1000.0_wp
  END DO
  WRITE (swc_line,'(I0,",",G14.6,",",G14.6,",",G14.6)') iD, root_dep, f_t, f_pr
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') theta (I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  WRITE (unit_rootdiag,'(A)') TRIM(swc_line)

  ! Compensatory-uptake diagnostics (day snapshot). f_comp = F/denom = fraction of demand taken.
  WRITE (unit_rootcomp,'(I0,6(",",G14.6))') iD, root_comp_omega, hydro_Fraw_diag, hydro_comp_denom_diag, &
    hydro_Fraw_diag / MAX(root_comp_eps, hydro_comp_denom_diag), hydro_root_sink_m_day * 1000.0_wp, trans_act
  DO I = 1, NSoil_layers
    WRITE (unit_rootcomp_lyr,'(I0,",",I0,5(",",G14.6))') iD, I, hydro_root_dist_diag(I), hydro_alpha_diag(I), &
      hydro_root_dist_diag(I) * hydro_alpha_diag(I), &
      hydro_root_dist_diag(I) * hydro_alpha_diag(I) / MAX(root_comp_eps, hydro_comp_denom_diag), &
      trans_act * hydro_root_dist_diag(I) * hydro_alpha_diag(I) / MAX(root_comp_eps, hydro_comp_denom_diag)
  END DO
  ! rec 80 Feddes-node census row (both modes; see the OPEN for the columns).
  WRITE (swc_line,'(I0,",",I0,",",I0,3(",",G14.6),",",I0,",",G14.6,",",I0,",",G14.6)') iD, feddes_p2_hydrus, &
    ft_pot_unstressed, hydro_tp_sent_cm_d_diag, hydro_tp_pot_cm_d_diag, hydro_p2_used_m_diag, hydro_p2_class_diag, &
    hydro_p2_alt_m_diag, hydro_p2_class_alt_diag, MAXVAL(ABS(hydro_alpha_diag - hydro_alpha_alt_diag))
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') hydro_alpha_diag(I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') hydro_alpha_alt_diag(I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  WRITE (unit_feddes_node,'(A)') TRIM(swc_line)

  ! Record 81 daily row, written in both modes. The pass-2 canopy transpiration
  ! is the accepted value with the soil-plant-atmosphere leaf; the supply
  ! shortfall is the unstressed demand minus that value, neither capped nor fed
  ! back, and zero when the record is off; and the actual transpiration is what
  ! the Feddes sink delivered against the unstressed demand.
  BLOCK
    REAL(wp) :: e_spac_mm, shortfall_spac_mm
    e_spac_mm = 0.0_wp; shortfall_spac_mm = 0.0_wp
    IF (spac_mode >= 1) THEN
      e_spac_mm = trans_pot_stressed_mm
      shortfall_spac_mm = trans_pot_unstressed_mm - e_spac_mm
    END IF
    WRITE (unit_spac_daily,'(I0,",",I0,2(",",G14.6),2(",",I0),4(",",G14.6),2(",",I0),5(",",G14.6))') iD, spac_mode, &
      spac_psi_f, spac_K_plant, spac_psi_s_valid, spac_psi_s_fallback, spac_psi_s_conv(1), spac_psi_s_conv(2), &
      spac_psi_s_conv(3), spac_psi_s_MPa, spac_calls_day, spac_iters_max_day, spac_resid_max_day, &
      trans_pot_unstressed_mm, e_spac_mm, trans_act, shortfall_spac_mm
    ! Record 173, on only: what was SENT to the sink this day, what the sink
    ! DELIVERED, the committed pass-2 canopy transpiration, their residual -
    ! conserved into sensible heat by the repartition path - and the repartition
    ! factor that carries it.
    IF (ssc_daily_open) THEN
      BLOCK
        REAL(wp) :: resid_mm, f_tr_rep
        resid_mm = e_spac_mm - trans_act
        f_tr_rep = 1.0_wp
        IF (e_spac_mm > 1.0e-9_wp) f_tr_rep = MIN(1.0_wp, trans_act / e_spac_mm)
        ssc_resid_abs_max_mm = MAX(ssc_resid_abs_max_mm, ABS(resid_mm))
        WRITE (unit_ssc_daily,'(I0,",",I0,",",I0,6(",",G14.6),",",I0)') iD, eb3_outer_iters_day, &
          MERGE(1,0,eb3_outer_converged), trans_pot_demand_mm, trans_act, e_spac_mm, resid_mm, &
          f_tr_rep, ssc_relax_scale, spac_fallback_accepts_run
      END BLOCK
    END IF
  END BLOCK

  ! Per-layer 1-D soil-temperature profile (mode 2 only; end-of-day Tsoil_prof).
  IF (tsoilprof_open .AND. ALLOCATED(Tsoil_prof)) THEN
    WRITE (swc_line,'(I0)') iD
    DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') Tsoil_prof(I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
    WRITE (unit_tsoilprof,'(A)') TRIM(swc_line)
  END IF

  ! ------------------------------------------------------------------------ !
  ! 6. Leaf senescence. The age and self-shading death rate follows the leaf
  !    death-rate structure of Teh (2006) Listing 8.7, and the nitrogen-driven
  !    component is the relative senescence rate of Yin et al. (2004) NJAS
  !    51:407-427, Appendix 2.
  ! ------------------------------------------------------------------------ !
  ! Record 69: the senescence law itself is unchanged, but it converts AREA to
  ! MASS with the specific leaf area it is handed. With the record on, the leaf
  ! area index it receives is the SINK area, so it must be handed the
  ! SINK-CONSISTENT specific leaf area, the reciprocal of the leaf mass per area,
  ! rather than the development-stage profile, or the area-to-mass conversion
  ! inside it would be on a different basis from the one that built the area.
  ! With the record off this is exactly the development-stage profile value as
  ! before.
  !
  ! The sink-consistent specific leaf area is the RECIPROCAL of the leaf mass per
  ! area, in m2 g-1, which is the stated unit of this argument: it is consumed
  ! in m2 g-1 at both places inside the senescence routine - the nitrogen-driven
  ! shed, which divides the excess leaf area index by it, and the senesced lamina
  ! area, which multiplies by it, against which the structural residual specific
  ! leaf nitrogen (g N m-2 of LEAF) is charged - and the other arm of the same
  ! selection, the development-stage profile, returns 0.018-0.026 m2 g-1. The
  ! leaf mass per area is in g DM m-2, so a factor of 1000 larger value would be
  ! in m2 kg-1 and would be wrong by three orders of magnitude. The consequences
  ! of that wrong unit are both verifiable: the senesced area would be 1000 times
  ! too large, which makes the nitrogen sent to the dead pool collapse to the
  ! whole carried amount, so that all senescing lamina nitrogen is buried and the
  ! mobile vegetative pool is identically zero; and the nitrogen-driven shed
  ! would be 1000 times too small, i.e. silently disabled. A wrong-unit value
  ! expresses no law, so there is no alternative behavior to select between: the
  ! reciprocal form is used on every path. The selection reaches this branch only
  ! when record 69 is on, and otherwise the development-stage profile is used
  ! unchanged.
  CALL apply_leaf_senescence(dvs2_growth, dvs2_rate_growth, &
    MERGE(1.0_wp / MAX(1.0e-9_wp, csl_LMA), sla_from_dvs(dvs2_growth), &
          canopy_sink_lai == 1 .AND. csl_active), lai_growth_diag, &
    W_green, W_dead, N_leaf, N_dead_leaf)
  ! B-2: termination-state bookkeeping (Martre 2003 conditions; days logged).
  IF (grain_fill_ended .AND. day_fill_ended == 0) day_fill_ended = iD
  IF (.NOT. canopy_senesced .AND. W_green <= 1.0e-6_wp .AND. dvs2_growth > 1.0_wp) THEN
    canopy_senesced = .TRUE.
  END IF
  IF (canopy_senesced .AND. day_canopy_dead == 0) day_canopy_dead = iD
END DO

WRITE (*,*) 'Simulation complete. Biomass: ', Biomass, ' gr m-2 year-1'
! validation hook: report the worst per-step heat-storage closure
! residual over the run (mode 2 only; should be ~round-off if the advance is exact).
IF (energy_balance_mode == 2) THEN
  WRITE (*,'(A,ES12.4,A)') ' Mode 2 soil-heat: max |heat-storage closure residual| = ', &
    soil_heat_resid_max, ' J m-2 per substep.'
END IF
! end-of-run AGB, same correction as the daily locus above.
agb_biomass = W_green + W_dead + W_dead_stem + W_stem + W_grain + W_chaff   ! O2(b): + dead stem/chaff
grain_yield_biomass = W_grain                          ! Task 7: grain only (not grain+chaff)
grain_protein_pct = 0.0_wp
! This end-of-run line uses GRAIN nitrogen over the grain dry matter, the same
! expression as the per-day CSV column, so that the chaff nitrogen is not
! attributed to the grain. The factor 5.7 is the standard
! nitrogen-to-protein conversion for wheat, carrying the same status as at the
! CSV place (Jones 1931 USDA Circular 183, not held on this machine); no new
! constant.
IF (W_grain > 1.0e-12_wp) grain_protein_pct = 100.0_wp * 5.7_wp * &
  (MAX(0.0_wp, N_storage - N_chaff) / W_grain)
WRITE (*,*) 'Final AGB:       ', agb_biomass, ' gr m-2'
WRITE (*,*) 'Final GrainYield:', grain_yield_biomass, ' gr m-2'
WRITE (*,*) 'Final GrainProt: ', grain_protein_pct, ' %'

storage_final = SUM(theta * tk) * 1000.0_wp + ponding_mm
WRITE (*,*) 'Mass Balance Report (mm):'
WRITE (*,*) 'Initial Storage (day1 start): ', storage_init
WRITE (*,*) 'Initial Storage (day1 end): ', storage_day1_end
WRITE (*,*) 'Final Storage:   ', storage_final
WRITE (*,*) 'Total Precip:    ', total_precip
WRITE (*,*) 'Total Irrig:     ', total_irrig
WRITE (*,*) 'Total Eff Supply:', total_supply_eff
WRITE (*,*) 'Total Intercept:', total_intercept
WRITE (*,*) 'Total ET (actual, soil-sourced):', total_et_out
IF (energy_balance_mode >= 1) THEN
  ! Wet-canopy interception evaporation (= total_intercept in mode 1) leaves the
  ! canopy store; ET_total = soil-sourced ET + interception evaporation.
  WRITE (*,*) 'Total ET incl interception (mode 1):', total_et_out + total_intercept
END IF
WRITE (*,*) 'Total Soil Evap: ', total_soil_evap
WRITE (*,*) 'Total Trans:     ', total_trans_act
WRITE (*,*) 'Total Drain:     ', total_drain
WRITE (*,*) 'Total Runoff:    ', total_runoff
WRITE (*,*) 'Total Pond Evap: ', total_pond_evap
WRITE (*,*) 'Total Evap (atm flux, mm): ', total_evap_atm_flux
WRITE (*,*) 'Total Exfil (head->pond, mm): ', total_evap_head_exfil
! This season total is the sum of the daily surface-ledger shortfall,
!   shortfall = MAX(0, (pr_eff + irr_eff) 1000 - infiltration - d(ponding)),
! i.e. the part of the day surface supply the surface ledger did not see enter
! either the profile or the pond. Irrigation is applied OUTSIDE that
! infiltration accounting, so the ledger books irrigation as shortfall, and the
! quantity can reach hundreds of mm without any water being lost. It is NOT a
! water-balance residual: conservation is the hydrology mass-balance line of
! the same log. It is also a DIFFERENT quantity from the supply shortfall of
! record 81, which is the unstressed demand minus the actual canopy
! transpiration and is a potential diagnostic.
WRITE (*,*) 'Cumulative surface-ledger unaccounted supply (mm) [was mislabelled Total Shortfall; ' // &
            'NOT a water-balance residual, NOT lost water -- see Total Hydro MB]: ', total_shortfall_mm
! One season-end line reporting the canopy / Richards OUTER fixed-point
! convergence census. The criterion is a change in the committed canopy
! temperature below the outer tolerance of 0.05 K between outer iterations,
! with a cap of 5 iterations; the state of the last iteration is committed
! whether or not the criterion was met, so this line says on how many days that
! happened. It is a report and not a gate.
WRITE (*,'(A,I0,A,I0,A,F8.4,A)') ' OUTER LOOP (canopy-Richards) unconverged days = ', &
  outer_days_unconverged, ' of ', outer_days_total, ' (criterion 0.05 K); max resid = ', &
  outer_resid_max_season, ' K'
! rec 154: three ON-only season-summary lines. NOTHING is written at
! rec 154 = 0, so a mode-0 run.log is unchanged.
IF (organ_n_turnover_nema >= 1) THEN
  WRITE (*,'(A,I0,A,4(A,F10.4))') ' rec 154 Bertheloot 2011a turnover (mode ', organ_n_turnover_nema, &
    ') season DEGRADATION, g N m-2:', '  lamina ', ont_cum_D_lam, '  stem ', ont_cum_D_stem, &
    '  chaff ', ont_cum_D_chaff, '  root ', ont_cum_D_root
  WRITE (*,'(A,4(A,F10.4),A,F10.4)') ' rec 154 Bertheloot 2011a turnover season SYNTHESIS, g N m-2:', &
    '  lamina ', ont_cum_S_lam, '  stem ', ont_cum_S_stem, '  chaff ', ont_cum_S_chaff, &
    '  root ', ont_cum_S_root, '  rec-97 cap clip ', ont_cum_syn_clip
  WRITE (*,'(A,F10.4,A,F10.4,A,I0,A,I0,A,I0,A,ES12.4)') ' rec 154 grain N from the pool, g N m-2: ', &
    ont_cum_grain, '  (uptake credited to the pool ', ont_cum_uptake, '); turnover days = ', ont_days, &
    ', grain POOL-BOUND days = ', ont_days_pool_bound, ', synthesis pro-rata days = ', ont_days_scaled, &
    '; max |organ-state N closure| = ', ont_resid_max
END IF
! rec 155: one ON-only season-summary line. NOTHING is written at
! rec 155 = 0, so a mode-0 run.log is unchanged.
IF (n_uptake_regulation >= 1) THEN
  WRITE (*,'(A,I0,A,F10.4,A,F10.4,A,F8.4,A,F8.4,A,I0,A,I0,A,F10.4)') &
    ' rec 155 uptake regulation (mode ', n_uptake_regulation, &
    ') season realized soil N uptake, g N m-2: PRE-anthesis ', N_uptake_pre_cum, &
    '  POST-anthesis ', N_uptake_post_cum, '; E_N mean ', &
    MERGE(nur_EN_sum / REAL(MAX(nur_days_live,1), wp), 1.0_wp, nur_days_live > 0), &
    ', min ', nur_EN_min, ' over ', nur_days_live, ' live days; limb (a) raised the demand on ', &
    nur_days_a_bind, ' days by a season total of ', nur_cum_a_gain
  WRITE (*,'(A,I0,A)') ' rec 155 limb (c) (grain demand net of the mobile pool) applied on ', &
    nur_days_c, ' post-anthesis days beyond what rec 97 already applies'
END IF
WRITE (*,*) 'Total Hydro MB:  ', total_hydro_mb_resid
! rec 149 wb_surface_limit >= 1 (ON-only echo; mode 0 run.log is, no new line).
IF (wb_surface_limit >= 1) THEN
  WRITE (*,*) 'Total Hydro MB Correction (mm, rec149 mode>=1 only): ', total_hydro_mb_correction
END IF
WRITE (*,*) 'Total Hydro Top Flux Net:', total_hydro_top_flux_m * 1000.0_wp
WRITE (*,*) 'Total Hydro Root Sink:  ', total_hydro_root_sink_m * 1000.0_wp
WRITE (*,*) 'Total Hydro Bottom Flux:', total_hydro_bottom_flux_m * 1000.0_wp
WRITE (*,*) 'Sum Balance Soil (mm):  ', total_balance_soil_mm
WRITE (*,*) 'Sum Balance Surface(mm):', total_balance_surface_mm
WRITE (*,*) 'Sum Balance Total (mm): ', total_balance_mm
WRITE (*,*) 'Sum Soil MB Resid (mm): ', total_soil_mb_resid_mm
WRITE (*,*) 'Full-System Closure (raw P+I, mm): ', &
  (storage_init + (total_precip + total_irrig) - &
   (total_soil_evap + total_pond_evap + total_trans_act + total_runoff + total_drain) - &
   storage_final)
WRITE (*,*) 'Full-System Closure (eff supply, mm): ', &
  (storage_init + total_supply_eff - &
   (total_soil_evap + total_pond_evap + total_trans_act + total_runoff + total_drain) - &
   storage_final)
WRITE (*,*) 'Evap Required (from storage+P+I, mm): ', &
  (storage_init + total_precip + total_irrig - total_runoff - total_drain - storage_final)
WRITE (*,*) 'Evap Reported (soil+pond+trans, mm): ', &
  (total_soil_evap + total_pond_evap + total_trans_act)
WRITE (*,*) 'Evap Missing (required - reported, mm): ', &
  (storage_init + total_precip + total_irrig - total_runoff - total_drain - storage_final - &
   (total_soil_evap + total_pond_evap + total_trans_act))
WRITE (*,*) 'WB Residual (daily sum, mm): ', total_wb_resid
WRITE (*,*) 'MB Error (eff):  ', (storage_final - storage_day1_end) - &
  (total_supply_eff - total_et_out - total_drain - total_runoff) - &
  total_hydro_mb_resid
WRITE (*,*) 'MB Error (raw):  ', (storage_final - storage_day1_end) - &
  ((total_precip + total_irrig) - total_et_out - total_drain - total_runoff) - &
  total_hydro_mb_resid
WRITE (*,*) 'MB Error (eff, start):  ', (storage_final - storage_init) - &
  (total_supply_eff - total_et_out - total_drain - total_runoff) - &
  total_hydro_mb_resid
WRITE (*,*) 'MB Error (raw, start):  ', (storage_final - storage_init) - &
  ((total_precip + total_irrig) - total_et_out - total_drain - total_runoff) - &
  total_hydro_mb_resid
WRITE (*,*) 'MB Error (hydro, start):  ', (storage_final - storage_init) - &
  ((total_hydro_top_flux_m - total_hydro_root_sink_m - total_hydro_bottom_flux_m) * 1000.0_wp)
Ci_last = Ci_mean_total / REAL(NDays,wp)
WRITE (*,*) 'Summary Ca/Ci/Tleaf/Tcanopy/Biomass/N/Tr/gsw', &
  Ca_s, Ci_last, (T_leaf_mean_total/REAL(NDays,wp)), (T_air_mean_total/REAL(NDays,wp)), &
  Biomass, Total_N_ass, Total_trans, Total_gsw_1
IF (hydro_fail_count > 0) THEN
  WRITE (*,'(A,I0,A)') 'Richards failed on ', hydro_fail_count, ' day(s):'
  WRITE (*,'(*(I0,1X))') (hydro_fail_days(i), i=1, hydro_fail_count)
ELSE
  WRITE (*,*) 'Richards had no failed days.'
END IF
! the same class the solver PRINTS, counted over every outer pass.
WRITE (*,'(A,I0,A,I0,A)') 'Richards non-convergence events over ALL outer passes: ', &
  hydro_fail_events_all, ' event(s) in ', hydro_fail_calls_all, ' solver call(s).'
IF (hydro_fail_events_all > 0 .AND. hydro_fail_count == 0) THEN
  WRITE (*,'(A)') '  NOTE: those events occurred on DISCARDED outer iterations only; the accepted' // &
    ' pass converged, which is why the accepted-pass count above is zero.'
END IF
WRITE (*,'(A,I0)') 'Richards evap-failure days (accepted pass, uptake-denied class): ', hydro_fail_evap_days
WRITE (*,'(A,I0)') 'Richards evaporative head-BC fallback retries taken: ', hydro_evap_fallback_events

! ---------------------------------------------------------------------------- !
! The canopy boundary-layer resistance and hourly-wind pair: the end-of-run
! report, emitted only when one of the records is on.
! ---------------------------------------------------------------------------- !
IF (rca_wind_attenuation == 1) THEN
  WRITE (*,'(A)')          'PART A rca_wind_attenuation (rec 70) = ON — Teh (2006) Eq. 4.78 p.97 in full'
  WRITE (*,'(A,F10.6)')    '  alpha (driver rec 12; Teh Table 4.2 p.93, Cionco 1972)        : ', a_wind
  WRITE (*,'(A,F10.6)')    '  attenuation factor alpha/(1-exp(-alpha/2)) the code restores  : ', rwa_att_factor
  WRITE (*,'(A,F10.6)')    '  sparse-canopy blend LAI (rec 72; DECLARED ASSUMPTION)         : ', rwa_lai_blend
  WRITE (*,'(A,F10.6)')    '  last-substep realised rca(ON)/rca(free-air)                   : ', rwa_rca_ratio
  ! The vanishing-attenuation assertion: half the attenuation integral must equal
  ! 1 + alpha/4 + alpha^2/48 + O(alpha^3), i.e. the resistance must collapse onto
  ! Teh Eq. 4.73, which is the unattenuated form, as the attenuation vanishes. It
  ! is reported as the worst second-order residual expressed as a FRACTION of its
  ! analytic bound, so that a value below 1 is a pass and the margin is legible
  ! rather than asserted.
  IF (rwa_lim_failed) THEN
    WRITE (*,'(A,F10.6)')  '  ALPHA->0 LIMIT ASSERTION: **FAILED**, worst residual/bound     : ', rwa_lim_worst
  ELSE
    WRITE (*,'(A,F10.6)')  '  alpha->0 limit assertion (Teh Eq. 4.73): PASS, worst resid/bound: ', rwa_lim_worst
  END IF
END IF
WRITE (*,'(A,I0,A,G12.5)') 'P2-e2 binned-LEs substeps above substep potential: ', les_bin_exceed_count, &
  '  max ratio: ', les_bin_exceed_max
WRITE (*,'(A,I0,A,I0)') 'B-2 grain-N termination: fill-ended day=', day_fill_ended, &
  '  canopy-senesced day=', day_canopy_dead
WRITE (*,'(A,G12.5)') 'B-2 assim-N kept vegetative after termination (gN/m2): ', assimN_surplus_postfill
WRITE (*,'(A,G12.5,A,G12.5)') 'B-1 post-maturity unallocated C (g CH2O/m2): ', unalloc_C_total, &
  '  sink-full surplus to stem reserve (g DM/m2): ', overflow_to_reserve_total

! Save results
OPEN (newunit=unit_flux,file='output/daily_fluxes.csv',status='replace',action='write')
WRITE (unit_flux,'(A)') 'day,Rn_MJ_m2_d,LE_MJ_m2_d,H_MJ_m2_d,G_MJ_m2_d,Tr_mm_d,E_mm_d,Ir_mm_d,Drain_mm_d,Runoff_mm_d'
DO iD = 1, NDays
  WRITE (unit_flux,'(I0,9(",",G14.6))') iD, (flux_daily(iD,J), J=1,9)
END DO
CLOSE(unit_flux)

OPEN (newunit=unit_flux,file='output/daily_fluxes_full.csv',status='replace',action='write')
WRITE (unit_flux,'(A)') 'day,Rn_MJ_m2_d,G_MJ_m2_d,LE_MJ_m2_d,LEc_MJ_m2_d,LEs_MJ_m2_d,' // &
  'H_MJ_m2_d,Hc_MJ_m2_d,Hs_MJ_m2_d,dRn_MJ_m2_d,Tr_mm_d,E_mm_d,Ir_mm_d,ET_mm_d,ET_Ir_mm_d,' // &
  'Cum_Tr_mm,Cum_E_mm,Cum_Ir_mm,Cum_ET_mm,Cum_ET_Ir_mm,Precip_mm_d,Cum_Precip_mm,' // &
  'Irrig_mm_d,Cum_Irrig_mm,Tleaf_C,Tair_C,Tleaf_minus_Tair_C,Rad_input_MJ_m2_d,' // &
  'Tmax_input_C,Tmin_input_C,Tmean_input_C,Ac_MJ_m2_d,As_MJ_m2_d,Storage_mm,Delta_Storage_mm,' // &
  'Daily_Balance_mm,Daily_Balance_NoIr_mm,Biomass_g_m2,Nitrogen_g_m2,N_to_B,' // &
  'Runoff_mm_d,Soil_MB_Residual_mm_d,Ponding_mm,Hydro_MB_Residual_mm_d,ET_pot_mm_d,' // &
  'PondEvap_mm_d,SoilEvap_mm_d,Infil_mm_d,Infil_Supply_mm_d,Infil_Pond_mm_d,' // &
  'Balance_Soil_mm_d,Balance_Surface_mm_d,TopFluxNet_mm_d,TopFluxUp_mm_d,' // &
  'HydroTopFluxNet_mm_d,HydroTopInfil_mm_d,HydroRootSink_mm_d,HydroBottom_mm_d,HydroMB_mm_d,HydroBalance_mm_d'
DO iD = 1, NDays
  WRITE (unit_flux,'(I0,59(",",G14.6))') iD, (flux_daily_full(iD,J), J=1,59)
END DO
CLOSE(unit_flux)
CLOSE(unit_pheno)

OPEN (newunit=unit_swc,file='output/daily_swc_layers.csv',status='replace',action='write')
OPEN (newunit=unit_swc_raw,file='output/daily_swc_layers_raw.csv',status='replace',action='write')
swc_header = 'day'
DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",layer",I0,"_mm")') I; swc_header = TRIM(swc_header)//TRIM(swc_tmp); END DO
WRITE (unit_swc,'(A)') TRIM(swc_header)//',total_mm'
WRITE (unit_swc_raw,'(A)') TRIM(swc_header)//',total_mm'
! Write day 0 (initial theta from driver_swtr.txt) so plots can show true initial state.
WRITE (swc_line,'(I0)') 0
DO I = 1, NSoil_layers
  WRITE (swc_tmp,'(",",G14.6)') theta_init(I)*tk(I)*1000.0_wp
  swc_line = TRIM(swc_line)//TRIM(swc_tmp)
END DO
WRITE (swc_tmp,'(",",G14.6)') SUM(theta_init * tk * 1000.0_wp)
WRITE (unit_swc,'(A)') TRIM(swc_line)//TRIM(swc_tmp)
WRITE (unit_swc_raw,'(A)') TRIM(swc_line)//TRIM(swc_tmp)
DO iD = 1, NDays
  WRITE (swc_line,'(I0)') iD
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') swc_layers_daily(iD,I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  WRITE (swc_tmp,'(",",G14.6)') swc_daily(iD,2); WRITE (unit_swc,'(A)') TRIM(swc_line)//TRIM(swc_tmp)
  WRITE (swc_line,'(I0)') iD
  DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",",G14.6)') swc_layers_raw_daily(iD,I); swc_line = TRIM(swc_line)//TRIM(swc_tmp); END DO
  WRITE (swc_tmp,'(",",G14.6)') SUM(swc_layers_raw_daily(iD,1:NSoil_layers))
  WRITE (unit_swc_raw,'(A)') TRIM(swc_line)//TRIM(swc_tmp)
END DO
CLOSE(unit_swc)
CLOSE(unit_swc_raw)

OPEN (newunit=unit_swc_mean,file='output/daily_swc_layers_mean.csv',status='replace',action='write')
swc_header = 'day'
DO I = 1, NSoil_layers; WRITE (swc_tmp,'(",layer",I0,"_mm")') I; swc_header = TRIM(swc_header)//TRIM(swc_tmp); END DO
WRITE (unit_swc_mean,'(A)') TRIM(swc_header)//',total_mm'
DO iD = 1, NDays
  WRITE (swc_line,'(I0)') iD
  DO I = 1, NSoil_layers
    WRITE (swc_tmp,'(",",G14.6)') swc_layers_mean_daily(iD,I)
    swc_line = TRIM(swc_line)//TRIM(swc_tmp)
  END DO
  WRITE (swc_tmp,'(",",G14.6)') SUM(swc_layers_mean_daily(iD,1:NSoil_layers))
  WRITE (unit_swc_mean,'(A)') TRIM(swc_line)//TRIM(swc_tmp)
END DO
CLOSE(unit_swc_mean)
CLOSE(unit_hydrodiag)
IF (wetfront_trace_enabled) THEN
  CLOSE(wetfront_trace_unit)
  wetfront_trace_enabled = .FALSE.
  wetfront_trace_unit    = -1
END IF
CLOSE(unit_wbdiag)
CLOSE(unit_hydroledger)
CLOSE(unit_soilprof)

OPEN (newunit=unit_swc,file='output/daily_swc_water.csv',status='replace',action='write')
WRITE (unit_swc,'(A)') 'day,water_supply_mm_d,swc_total_mm,swc_root_mm,avail_water_mm'
DO iD = 1, NDays
  WRITE (unit_swc,'(I0,4(",",G14.6))') iD, swc_daily(iD,1), swc_daily(iD,2), swc_daily(iD,3), swc_daily(iD,4)
END DO
CLOSE(unit_swc)

OPEN (unit=1,file='output/test_energy_bal_2025.csv',status='unknown')
DO iD = 1, NDays; CALL csv_write(1, A_data (iD,1:5),.true.); END DO
CLOSE(1)
CLOSE(unit_diurnal)
IF (flag_diag_open) THEN ; CLOSE(unit_flag_diag) ; flag_diag_open = .FALSE. ; END IF
CLOSE(unit_spac_daily)
IF (ssc_daily_open) THEN ; CLOSE(unit_ssc_daily) ; ssc_daily_open = .FALSE. ; END IF
! psi_L solver robustification: report the bounded fallback ONLY if it ever fired, so a run that
! never needs it writes nothing (the OFF run.log stays). SHARED.f90 (a)/(b).
IF (spac_fallback_accepts_run > 0) THEN
  WRITE (*,'(A,I0,A)') ' SPAC PSI_L: bounded fallback accepted ', spac_fallback_accepts_run, &
    ' solve(s) with 1e-6 MPa <= |g| < 1e-5 MPa (HEAD would have ERROR STOPped the run).'
END IF
IF (eb_diag_open) CLOSE(unit_eb_diag)
IF (eb_act_open) CLOSE(unit_eb_act)
IF (rad_open) CLOSE(unit_rad)
IF (tlc_open) CLOSE(unit_tlc)    ! (rec 95) ON-only diagnostic
IF (esoil_pair_open) CLOSE(unit_esoil_pair)
IF (node_pair_open) CLOSE(unit_node_pair)
IF (soilheat_h_open) CLOSE(unit_soilheat_h)
IF (soilheat_d_open) CLOSE(unit_soilheat_d)
IF (tsoilprof_open) CLOSE(unit_tsoilprof)
IF (energy_balance_mode >= 1 .AND. eb3_implausible_count > 0) THEN
  WRITE (*,'(A,I0,A)') ' WARNING (three-node EB): ', eb3_implausible_count, &
    ' daytime substeps converged with Tcan outside [Tair-15, Tair+25] C.'
  WRITE (*,*) '   These are energy-conserving solutions (not clamped); inspect'
  WRITE (*,*) '   output/hourly_energy_balance.csv. Likely driver: large raa at'
  WRITE (*,*) '   sparse canopy / low gsw routing available energy into sensible heat.'
END IF
CLOSE(unit_rootdiag)
CLOSE(unit_tleaf)
CLOSE(unit_resist)
CLOSE(unit_growth)
IF (leafn_cap_open) CLOSE(unit_leafn_cap)   ! rec 87 ON-only diagnostic
CLOSE(unit_yield)
CLOSE(unit_minerals)
CLOSE(unit_organic)
IF (humusN_open) CLOSE(unit_humusN)   ! rec 93 ON-only diagnostic
IF (fnc_open) CLOSE(unit_fnc)         ! rec 97 ON-only diagnostic
! rec 121 end-of-run summary (ON only; unconditional position, outside every other record's block)
IF (n_payment_gate >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,I0,A,I0)') ' rec 121 N PAYMENT: live days = ', npg_days_live, &
    ' ; binding days L / S / R = ', npg_days_bind_L, ' / ', npg_days_bind_S, ' / ', npg_days_bind_R
  WRITE (*,'(A,G14.6,A,G14.6,A,G14.6)') ' rec 121 N paid cum L / S+chaff / R (g N m-2) = ', npg_paid_cum_L, ' / ', &
    npg_paid_cum_S, ' / ', npg_paid_cum_R
  WRITE (*,'(A,G14.6,A,G14.6)') ' rec 121 root tier paid from root-reduced N / from the shoot pool (g N m-2) = ', &
    npg_paid_cum_R_root, ' / ', npg_paid_cum_R_shoot
  WRITE (*,'(A,G14.6,A,G14.6,A,G14.6)') ' rec 121 unpaid C cum (g CH2O m-2) = ', npg_unpaid_ch2o_cum, &
    ' ; to stem WSC (g DM m-2) = ', npg_wsc_cum, ' ; explicit residual (g CH2O m-2) = ', npg_resid_cum
  WRITE (*,'(A,G14.6,A)') ' rec 121 gate N identity max|resid| = ', npg_ident_N_max, ' g N m-2 (machine zero expected)'
  WRITE (*,'(A,G14.6,A)') ' rec 121 lamina DM built ahead of its nitrogen, cum (g DM m-2) = ', npg_lamina_unpaid_cum, &
    '  (Batch 2 closes this)'
  WRITE (*,'(A,G14.6,A,G14.6,A,G14.6)') ' rec 121 CARBON CLOSURE max|organ-state - dm_gain| = ', bgn_dm_close_max, &
    ' ; max|[SFL] supply ledger| = ', sfl_ledger_max, ' ; max|rec-100 identity - refused - N-unpaid| = ', sdm_close_max
END IF
IF (lai_n_arrest >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G14.6,A,G14.6,A)') ' rec 122 LAMINA N ARREST: arrest days = ', lna_days_arrest, &
    ' ; shed days = ', lna_days_shed, ' ; unpaid area removed cum (m2 m-2) = ', lna_removed_cum, &
    ' (= ', lna_removed_dm_cum, ' g DM m-2 at csl_LMA)'
  WRITE (*,'(A,G14.6,A,I0,A,I0)') ' rec 122 shed area cum (m2 m-2) = ', lna_shed_cum, ' ; shed-unmet days = ', &
    lna_shed_unmet_days, ' ; expanding-while-senesced rank-days = ', lna_expdead_rank_days
  WRITE (*,'(A)') ' rec 122 standing lamina area built ahead of its nitrogen at any day''s end = 0 by construction (i)'
END IF
IF (tiller_n_afford >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G14.6,A,G14.6)') ' rec 123 TILLER N AFFORDABILITY: window days stepped = ', tna_days_step, &
    ' ; days with phi_N < phi_C = ', tna_days_bind, ' ; final phi_N main stem = ', tna_phiN(0), ' ; T1 = ', tna_phiN(1)
  WRITE (*,'(A,6(G12.4,1X))') ' rec 123 final phi (MIN of C and N) per cohort 0..5 = ', (tca_phi(i_tc), i_tc = 0, CT_MAXTIL)
END IF
IF (root_front_law >= 1) THEN
  WRITE (*,'(A,F5.1,A,I0,A,G14.6,A,G14.6,A,I0)') ' rec 124 ROOT FRONT: r_z = ', rfl_rz, ' mm d-1 ; days since emergence = ', &
    rfl_days_emerged, ' ; final z_root (m) = ', rfl_z_new, ' (legacy Borg & Grimes would be ', rfl_z_legacy, &
    ') ; day z_max reached = ', rfl_day_zmax
END IF
IF (kernel_potential_law >= 1) THEN
  WRITE (*,'(A,I0,A,G14.6,A,G14.6,A,I0)') ' rec 125 GRAIN-N CEILING: mode = ', kernel_potential_law, &
    ' ; cap (g N m-2) = ', gnc_cap, ' ; cut cumulative (g N m-2) = ', gnc_cut_cum, ' ; days bound = ', gnc_days_bound
  IF (kernel_potential_law >= 2) WRITE (*,'(A,I0,A,G14.6,A,G14.6,A,G14.6,A,I0,A,G14.6,A,G14.6)') &
    ' rec 125 KERNEL POTENTIAL: set on day ', kpt_day_set, ' ; s_k mean (mg kernel-1 d-1) = ', kpt_sk_mean, &
    ' ; W_raw (mg) = ', kpt_w_raw, ' ; W_pot (mg) = ', kpt_w_pot, ' ; bound (0 interior/1 min/2 max) = ', kpt_bound, &
    ' ; envelope ', kpt_env_min, ' - ', kpt_env_max
END IF
IF (mineral_grain_sink_cap == 1) THEN
  BLOCK
    INTEGER :: u_mgs, i_mgs
    REAL(wp) :: gm
    OPEN (newunit=u_mgs, file='output/mineral_sink_cap.csv', status='replace', action='write')
    WRITE (u_mgs,'(A)') 'mineral,grain_max_mg_kg,refused_cum_mg_m2,days_bound,grain_mg_m2_final,grain_mg_kg_final'
    DO i_mgs = 1, NMIN
      gm = 0.0_wp
      IF (W_grain > 1.0e-9_wp) gm = 1000.0_wp * M_grain(i_mgs) / W_grain
      WRITE (u_mgs,'(A,",",G14.6,",",G14.6,",",I0,",",G14.6,",",G14.6)') mineral_names(i_mgs), mm_grain_max(i_mgs), &
        M_cap_refused(i_mgs), M_cap_days(i_mgs), M_grain(i_mgs), gm
      WRITE (*,'(A,A,A,G12.4,A,G12.4,A,I0,A,G12.4)') ' rec 126 SINK BOUND ', mineral_names(i_mgs), ': m_max (mg kg-1) = ', &
        mm_grain_max(i_mgs), ' ; refused cumulative (mg m-2) = ', M_cap_refused(i_mgs), ' ; days bound = ', M_cap_days(i_mgs), &
        ' ; final grain (mg kg-1) = ', gm
    END DO
    CLOSE (u_mgs)
  END BLOCK
END IF
IF (repro_surplus_fate == 1) THEN
  WRITE (*,'(A,G12.4,A,G12.4,A,G12.4,A,I0)') ' rec 127 REPRO SURPLUS FATE: sink-full surplus cum (g DM m-2) = ', rsf_surplus_cum, &
    ' ; labeled stem WSC = ', rsf_label_cum, ' ; refused to explicit residual (never structural) = ', rsf_resid_cum, &
    ' ; days refused = ', rsf_days
END IF
IF (n_construction_status >= 1) THEN
  WRITE (*,'(A,I0,A,G12.4,A,I0)') ' rec 128 CONSTRUCTION N STATUS: days f_crit < 1 = ', ncs_days_fcrit, &
    ' ; cumulative lamina DM built below target (g DM m-2) = ', ncs_lamina_below_target_cum, &
    ' ; days p_L < 1 = ', npg_days_bind_L
END IF
! rec 156: one ON-only season-summary line. Nothing is written at
! rec 156 = 0, so a mode-0 run.log is unchanged.
IF (n_construction_ceiling_curve >= 1 .AND. n_construction_status >= 1) THEN
  WRITE (*,'(A,F6.3,A,F6.4,A,I0,A)') ' rec 156 CONSTRUCTION CEILING CURVE: pass-2 target capped on the MAXIMUM' // &
    ' dilution curve, a = ', ncc_a_used, ', b = ', ncc_b_used, ' , evaluated on ', ncc_days_live, &
    ' days (the days f_crit < 1 count above is on that curve).'
END IF
! rec 157: one ON-only season-summary line. Nothing is written at
! rec 157 = 0, so a mode-0 run.log is unchanged.
IF (lamina_n_storage >= 1) THEN
  WRITE (*,'(A,I0,A,F5.3,A,G12.4,A,I0,A,I0,A,I0)') ' rec 157 LAMINA N STORAGE (mode ', lamina_n_storage, &
    '): ceiling = ', sln_max_N_lnh, ' g N m-2 leaf ; season N stored in the lamina above the construction' // &
    ' target (g N m-2) = ', lns_stored_cum, ' ; ladder days = ', lns_days_live, &
    ' ; days storing > 0 = ', lns_days_store, ' ; days the holding capacity BOUND = ', lns_days_bind
  IF (lamina_n_storage == 3) THEN
    WRITE (*,'(A,G12.4,A,G12.4,A,G12.4,A,G12.4,A,G12.4)') ' rec 157 mode 3 [N]mob NUMERATOR at the last' // &
      ' evaluated day (g N m-2): total = ', lns_nmob_num_day, ' = N_mobile_veg ', lns_nmob_t_mobile, &
      ' + N_NO3_shoot ', lns_nmob_t_no3, ' + stem luxury ', lns_nmob_t_stemlux, &
      ' + chaff luxury ', lns_nmob_t_chaflux
  END IF
END IF
! rec 158 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (nitrate_drain_rate_basis >= 1) THEN
  WRITE (*,'(A,I0,A,G12.4,A,G12.4,A,G12.4,A,I0,A,I0)') &
    ' rec 158 NITRATE DRAIN RATE BASIS (mode ', nitrate_drain_rate_basis, &
    '): season mean N_max ACTUALLY used (umol m-2leaf s-1) = ', &
    ndb_nmax_used_cum / MAX(1.0_wp, REAL(ndb_days_live, wp)), &
    ' ; season mean N_max the incumbent STORE basis would have used = ', &
    ndb_nmax_store_cum / MAX(1.0_wp, REAL(ndb_days_live, wp)), &
    ' ; season sum of the shoot de novo N-reduction demand (g N m-2) = ', ndb_rate_cum, &
    ' ; N_max days = ', ndb_days_live, &
    ' ; days the Busch enzymatic ceiling N_max_ref BOUND = ', ndb_days_cap
END IF
! rec 159 season summary -- ON-only, so a mode-0 run.log is unchanged.
! rec 160 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (justes_nmax_floor_basis >= 1) THEN
  WRITE (*,'(A,F6.3,A,I0,A,I0)') &
    ' rec 160 JUSTES Nmax FLOOR BASIS (mode 1): sub-floor Nmax cap used (%DM) = ', jnf_cap_pct, &
    ' (the printed Eqn (17) curve at its own floor, 8.3 x 1.55^-0.44; the legacy mode-0 value is the undocumented 6.6 %)' // &
    ' ; days with shoot DM <= 1.55 t ha-1 = ', jnf_days_floor, &
    ' ; of those, days the Nmax stock gap BOUND the day''s N demand = ', jnf_days_bind
END IF
! rec 163 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (postanth_uptake_regulation >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G12.4,A,G12.4,A,G12.4)') &
    ' rec 163 POST-ANTHESIS UPTAKE REGULATION (mode 1, Bertheloot 2011a Eqn 3-5, no demand term): post-anthesis' // &
    ' days the demand gate was removed = ', pur_days_post, &
    ' ; of those, days the removal actually RAISED U_pot = ', pur_days_bind, &
    ' ; season sum of the extra uptake POTENTIAL granted (g N m-2) = ', pur_cum_gain, &
    ' ; mean E_N over those days = ', pur_EN_sum / MAX(1.0_wp, REAL(pur_days_post, wp)), &
    ' ; minimum E_N = ', pur_EN_min
END IF
! rec 164 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (lamina_capacity_profile_basis >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G12.4,A,G12.4,A,F5.3)') &
    ' rec 164 LAMINA CAPACITY PROFILE BASIS (mode 1, Pask 2009 p.29 2.5 g N m-2 at the TOP of the' // &
    ' model''s own canopy N profile): days the profile factor was used = ', lcp_days_profile, &
    ' ; days the mode-0 value was used because the profile was not defined = ', lcp_days_fallback, &
    ' ; mean n_top/n_mean over the profile days = ', &
    lcp_factor_sum / MAX(1.0_wp, REAL(lcp_days_profile, wp)), &
    ' ; canopy-mean capacity AT that mean factor, 2.5 / mean(n_top/n_mean) (g N m-2) = ', &
    sln_max_N_lnh * REAL(lcp_days_profile, wp) / MAX(1.0e-9_wp, lcp_factor_sum), &
    ' ; the top of the canopy is held at the printed ', sln_max_N_lnh
END IF
! rec 165 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (wsc_release_delay_basis >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G12.4,A,G12.4,A,I0)') &
    ' rec 165 WSC RELEASE DELAY (mode 1, Schnyder 1993 p.238 10-20 d accumulation window,' // &
    ' wsc_accum_end_dpa = ', wsc_accum_end_dpa, &
    ' d): days the rec-92 deficit limb WOULD have drawn and was HELD = ', wrd_days_delayed, &
    ' ; cumulative unmet grain demand the pool was not allowed to serve on those days (g DM m-2) = ', &
    wrd_room_held_cum, &
    ' ; largest pool held (g DM m-2) = ', wrd_pool_held_max, &
    ' ; first draw at (d after anthesis; -999 = never) = ', &
    MERGE(wsc_day_mob0 - wsc_day_anth, -999, wsc_day_mob0 > 0 .AND. wsc_day_anth > 0)
END IF
! rec 166 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (grain_ca_xylem_throughput >= 1) THEN
  WRITE (*,'(A,G12.4,A,G12.4,A,G12.4,A,G12.4,A,G12.4,A,I0)') &
    ' rec 166 GRAIN CALCIUM XYLEM THROUGHPUT: cumulative ear xylem water V_ear (L m-2) = ', gcx_v_ear_cum, &
    ' ; cumulative Lizana net grain water import V_net (L m-2) = ', gcx_v_net_cum, &
    ' ; Ca into grain by this route (mg m-2) = ', gcx_ca_grain_cum, &
    ' ; Ca into chaff (mg m-2) = ', gcx_ca_chaff_cum, &
    ' ; largest daily ear share of green area (-) = ', gcx_f_ear_max, &
    ' ; days the route moved Ca = ', gcx_days_route
END IF
! ==================================================================== !
! Record 180 season summary, on only. Ledgers, never fluxes.
! ==================================================================== !
IF (grain_ca_loading_basis >= 1) THEN
  WRITE (*,'(A,F12.4,A,G12.4,A,G12.4,A,G12.4,A,G12.4,A,G12.4,A,I0,A,I0)') &
    ' rec 180 GRAIN CALCIUM LOADING BASIS (xylem stream, SHARED.f90): c_xyl_Ca (mg L-1) = ', &
    gcl_cxyl_used, &
    ' ; cumulative Tr_act seen by the law (L m-2) = ', gcl_tr_season_cum, &
    ' ; of which post-anthesis (L m-2) = ', gcl_tr_post_cum, &
    ' ; Ca to the GRAIN terminus (mg m-2) = ', gcl_ca_grain_cum, &
    ' ; Ca to the CHAFF terminus (mg m-2) = ', gcl_ca_chaff_cum, &
    ' ; Ca to the VEGETATIVE terminus (mg m-2) = ', gcl_ca_veg_cum, &
    ' ; days the post-anthesis partition ran = ', gcl_days_post, &
    ' ; days the labile Ca pool bounded the law = ', gcl_days_pool
  WRITE (*,'(A,F10.4,A)') '        INDEPENDENT CHECK, never a tuning target: whole-crop Ca' // &
    ' delivered this season (g m-2) = ', &
    (gcl_ca_grain_cum + gcl_ca_chaff_cum + gcl_ca_veg_cum) * 1.0e-3_wp, &
    ' against 2.32 g m-2 measured at harvest by Gregory PJ, Crawford DV & McGowan M (1979)' // &
    ' J Agric Sci 93:485-494, DOI 10.1017/S0021859600038268, Table 1 p.491 (leaf 0.89 +' // &
    ' stem 0.74 + ear 0.27 + root 0.42); above ground 1.90 g m-2.'
END IF
! ==================================================================== !
! Record 167 season summary, on only.
! ==================================================================== !
IF (lamina_synthesis_sheath_basis >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,I0,A,G12.4,A,G12.4,A,G12.4,A,G12.4)') &
    ' rec 167 LAMINA SYNTHESIS SHEATH BASIS (mode ', lamina_synthesis_sheath_basis, &
    '): days the rebasing applied = ', lsb_days_applied, &
    ' ; mode-1 days with no area, mode-0 basis stood = ', lsb_days_no_area, &
    ' ; mean M_sheath over the applied days (g DM m-2) = ', &
    lsb_msheath_cum / MAX(1.0_wp, REAL(lsb_days_applied, wp)), &
    ' ; realized M_sheath/W_green on the last applied day (-) = ', lsb_ratio_last, &
    ' ; cumulative EXTRA lamina synthesis asked for before s_scale and the clip (g N m-2) = ', &
    lsb_syn_lam_gain_cum, &
    ' ; cumulative stem synthesis removed (g N m-2) = ', lsb_syn_stem_loss_cum
END IF
! rec 168 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (mineral_impedance_basis >= 1) THEN
  WRITE (*,'(A,G12.4,A,G12.4,A,G12.4,A,I0,A,I0,A,G12.4,A,G12.4)') &
    ' rec 168 MINERAL IMPEDANCE BASIS (Olesen 2001 p.94 Eq.7a CSIF): theta_th (m3 m-3) = ', nud_theta_th, &
    ' ; season-mean f_l under Eq.(7a) = ', mib_fl_sum / MAX(1.0_wp, REAL(mib_fl_days, wp)), &
    ' ; season-mean f_l under the retired MAX(0, 1.60*theta - 0.17) pair on the SAME days = ', &
    mib_fl_off_sum / MAX(1.0_wp, REAL(mib_fl_days, wp)), &
    ' ; ELEMENT-days the live diffusion locus ran (5 non-Ca elements x days) = ', mib_fl_days, &
    ' ; element-days theta <= theta_th so Eq.(7b) returned exactly 0 = ', mib_days_at_zero, &
    ' ; vwc_root range on those days (m3 m-3) = ', mib_vwc_min, ' to ', mib_vwc_max
END IF
! Record 183 season summary, on only. ONE LINE PER ELEMENT: the distribution
! coefficient the buffer power was built from, the bulk density, and the season
! minimum, mean and maximum of the daily buffer power the diffusion coefficient
! actually divided by. The per-element supply decomposition it explains is in
! output/daily_mineral_uptake_terms.csv (record 140, the mass-flow and
! diffusion columns).
IF (mineral_buffer_power_basis >= 1) THEN
  WRITE (*,'(A,F8.4,A,L1,A)') ' rec 183 MINERAL BUFFER POWER BASIS (b_e = theta + rho_b*Kd_e;' // &
    ' Baldwin, Nye & Tinker 1973 Plant Soil 38:621 p.622 "b = buffer power of soil = dC/dCl"' // &
    ' and p.625 "Db = Dl.W.fl"; Roose & Fowler 2004 J Theor Biol 228:173 Eq.(2.8) p.175):' // &
    ' rho_b used (kg L-1) = ', bpb_rho_b_used, ' ; module-basis fallback used = ', bpb_rho_fallback, &
    ' ; the driver_minerals.csv column-19 constant buffer_power_b was NOT read by De_e.'
  DO imin_loop = 1, NMIN
    IF (bpb_b_days(imin_loop) > 0) THEN
      WRITE (*,'(A,A,A,F12.4,A,F12.2,A,F12.2,A,F12.2,A,I0,A,F12.4,A,I0)') &
        ' rec 183   element ', TRIM(mineral_names(imin_loop)), ': Kd (L kg-1) = ', &
        bpb_kd_used(imin_loop), ' ; b_e min = ', bpb_b_min(imin_loop), &
        ' ; b_e season mean = ', bpb_b_sum(imin_loop) / MAX(1.0_wp, REAL(bpb_b_days(imin_loop), wp)), &
        ' ; b_e max = ', bpb_b_max(imin_loop), &
        ' ; days = ', bpb_b_days(imin_loop), &
        ' ; the OFF constant it replaced (column 19) = ', mm_b(imin_loop), &
        ' ; days b_e < 1 (WARNING) = ', bpb_warn_lo(imin_loop)
      IF (bpb_warn_lo(imin_loop) > 0) THEN
        WRITE (*,'(A,A,A,I0,A)') ' WARNING rec 183: element ', TRIM(mineral_names(imin_loop)), &
          ' computed a buffer power below 1.0 on ', bpb_warn_lo(imin_loop), &
          ' day(s), i.e. this record INCREASED its diffusive supply instead of damping it.' // &
          ' b = 1 is the NON-ADSORBED limit (So & Nye 1989 J Soil Sci 40:743 p.744 Eq.(3)' // &
          ' "for a non-adsorbed ion, D = Dl fl"); a value below it is not expected for a' // &
          ' sorbing cation. SHARED.f90.'
      END IF
    END IF
  END DO
END IF
! Record 171 season summary, on only. One line per element with a non-zero
! stoichiometry; a run in which every stoichiometry is 0.0, which is the
! delivered input, prints the inert line and nothing else, which is itself the
! proof that the bound did not act.
IF (mineral_n_coupling_basis >= 1) THEN
  IF (MAXVAL(mm_grain_eN) <= 0.0_wp) THEN
    WRITE (*,'(A,I0,A)') ' rec 171 GRAIN PROTEIN MINERAL SINK (mode ', mineral_n_coupling_basis, &
      ') : INERT on this deck -- driver_minerals.csv column 20 grain_eN_mg_gN is 0.0 for all six' // &
      ' elements (the stoichiometry is DECLARED ABSENT, SHARED.f90), so no day was bounded' // &
      ' and this run is numerically identical to mode 0.'
  ELSE
    DO imin_loop = 1, NMIN
      IF (mm_grain_eN(imin_loop) > 0.0_wp) THEN
        WRITE (*,'(A,I0,A,A,A,G12.4,A,G12.4,A,I0)') &
          ' rec 171 GRAIN PROTEIN MINERAL SINK (mode ', mineral_n_coupling_basis, ') ', &
          TRIM(mineral_names(imin_loop)), &
          ' : stoichiometry r_e (mg per g grain N) = ', mm_grain_eN(imin_loop), &
          ' ; grain loading REFUSED by the protein sink (cum, mg m-2) = ', mnc_refused(imin_loop), &
          ' ; days the sink bound = ', mnc_days(imin_loop)
      END IF
    END DO
  END IF
END IF
! rec 174 season summary -- ON-only, so a mode-0 run.log is unchanged.
! rec 175 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (grain_n_partition_basis >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,G12.4,A,G12.4)') &
    ' rec 175 GRAIN N PARTITION BASIS (mode ', grain_n_partition_basis, &
    '): MODE 1, days the restored continuity term was > 0 = ', gnp_days_ndb_add, &
    ' ; cumulative continuity term entering ndb_rate_prev (g N m-2) = ', gnp_ndb_add_cum, &
    ' ; shoot nitrate standing at the last day (g N m-2) = ', N_NO3_shoot
  IF (grain_n_partition_basis >= 4) THEN
    WRITE (*,'(A,G12.4,A,G12.4,A,F8.4,A,G12.4,A,G12.4)') &
      ' rec 175 MODE 4 (Ngtot(0) paid against the GRAIN''S OWN dry mass at the measured' // &
      ' flowering grain [N] 0.020 g N g-1 DM, Bertheloot 2011b Table 1(a)): Ngtot(0) target (g N m-2) = ', &
      gnp_g0_target, ' ; cumulative paid by the chaff pool = ', gnp_g0_take, &
      ' ; fraction of the law''s own initial condition the grain''s dry mass earned = ', gnp_g0_frac, &
      ' ; STANDING arrears the chaff structural floor refused (DISCLOSED, a state not a sum) = ', &
      gnp_g0_refused, ' ; still outstanding against Ngtot(0) = ', gnp_g0_deficit
  END IF
  IF (grain_n_partition_basis >= 3) THEN
    WRITE (*,'(A,I0,A,G12.4,A)') &
      ' rec 175 MODE 3 -- the PRE-REGISTERED prediction is INERT: days the pro-rata solve' // &
      ' released N to the grain = ', gnp_m3_days, ' ; cumulative release (g N m-2) = ', &
      gnp_m3_release_cum, '.  A non-zero value here FALSIFIES the prediction and must be reported, not absorbed.'
  END IF
END IF
! rec 178 season summary -- modes 2/3 ONLY, so a mode-0 and a mode-1
! run.log keep their bytes. Every number here is STATE at the last day.
IF (chaff_structural_n_basis >= 2) THEN
  WRITE (*,'(A,I0,A,F6.3,A,G12.4,A,G12.4)') &
    ' rec 178 CHAFF NRE LAW (mode ', chaff_structural_n_basis, &
    '): NRE_chaff (Pask 2012 FCR 126:104 Table 4 p.110 col CH; Barraclough 2014 FCR 156:242 p.247) = ', &
    cnr_nre_chaff, ' ; chaff N at anthesis, snapshot on the first live fill day (g N m-2) = ', &
    cnr_N_chaff_anth, ' ; the floor the law sets, (1 - NRE) x that (g N m-2) = ', cnr_floor_gN
  WRITE (*,'(A,G12.4,A,G12.4,A,I0,A,G12.4)') &
    ' rec 178 WHICH LOWER BOUND BINDS: chaff N at the last day (g N m-2) = ', N_chaff, &
    ' ; the concentration floor in force x W_chaff (g N m-2) = ', &
    n2_struct_chaff_use * MAX(0.0_wp, W_chaff), &
    ' ; fill days on which the NRE floor was the binding one = ', cnr_days_bound, &
    ' ; its excess over the concentration floor at the last day (g N m-2) = ', cnr_held_gN
  WRITE (*,'(A,G12.4)') &
    ' rec 178 REALIZED CHAFF NRE, model (1 - N_chaff(last) / N_chaff(anthesis)) -- a FALSIFIABLE' // &
    ' number, compare it with the 0.510-0.633 band of Pask Table 4 N-level means = ', &
    MERGE(1.0_wp - N_chaff / MAX(cnr_N_chaff_anth, 1.0e-12_wp), -999.0_wp, cnr_N_chaff_anth > 1.0e-12_wp)
END IF
! rec 176 season summary -- ON-only, so a mode-0 run.log is unchanged.
IF (sink_full_carbon_fate == 1) THEN
  WRITE (*,'(A,G12.4,A,G12.4,A,G12.4,A,I0,A,I0)') &
    ' rec 129 SINK-FULL CARBON FATE: surplus stored as stem WSC under the lifted window (cum, g DM m-2) = ', cfo_wsc_cum, &
    ' ; unallocated-assimilate carbon loss (cum, g DM m-2, NOT respiration, NOT structure) = ', cfo_loss_cum, &
    ' ; W_reserve at harvest (g m-2) = ', W_reserve, &
    ' ; days storing = ', cfo_days_store, ' ; days with overflow loss = ', cfo_days_loss
END IF
IF (grain_n_min_conc >= 1) THEN
  WRITE (*,'(A,I0,A,G12.4,A,G12.4,A,G12.4)') &
    ' rec 130 GRAIN N MINIMUM: mode = ', grain_n_min_conc, &
    ' ; days the floor bound grain DM import = ', gnm_days_bind, &
    ' ; refused grain DM (cum, g DM m-2, routed to the rec 127/129 fate chain) = ', gnm_refused_dm_cum, &
    ' ; final grain N concentration (%) = ', 100.0_wp * MAX(0.0_wp, N_storage - N_chaff) / MAX(1.0e-9_wp, W_grain)
  IF (grain_n_min_conc == 3) THEN
    WRITE (*,'(A,G12.4,A,G12.4,A,G12.4)') &
      ' rec 130 mode 3 (Martre 2003 soft structural floor): final gnm_c_stru_day (%DM) = ', 100.0_wp * gnm_c_stru_day, &
      ' ; final Nstru(Tt) (g N m-2, Eq. 3 cumulative demand) = ', gnm_Nstru_total, &
      ' ; cumulative structural N paid (g N m-2, diagnostic) = ', gnm_Nstru_paid
  END IF
END IF
IF (mineral_supply_law == 1) THEN
  BLOCK
    INTEGER :: i_msl2
    REAL(wp) :: gm2
    WRITE (*,'(A)') ' rec 131 MINERAL SUPPLY LAW (calibrated availability), final grain concentration (mg/kg):'
    DO i_msl2 = 1, NMIN
      gm2 = 0.0_wp
      IF (W_grain > 1.0e-9_wp) gm2 = 1000.0_wp * M_grain(i_msl2) / W_grain
      WRITE (*,'(A,A,A,G14.6,A,G14.6,A,G14.6)') '   ', mineral_names(i_msl2), ' : soln/xylem (mg/L) = ', M_soln(i_msl2), &
        ' ; soil labile remaining (mg/m2) = ', M_soil_labile(i_msl2), ' ; grain (mg/kg) = ', gm2
    END DO
  END BLOCK
END IF
IF (mineral_supply_mech == 1) THEN
  BLOCK
    INTEGER :: i_msm2
    WRITE (*,'(A)') ' rec 140 MINERAL SUPPLY MECH (diffusion/HATS), season REALIZED uptake (mg/m2) and final grain (mg/kg):'
    DO i_msm2 = 1, NMIN
      WRITE (*,'(A,A,A,G14.6,A,G14.6,A,G14.6,A,G14.6)') '   ', mineral_names(i_msm2), ' : S_mf = ', msm_s_mf_cum(i_msm2), &
        ' ; S_diff = ', msm_s_diff_cum(i_msm2), ' ; S_hats = ', msm_s_hats_cum(i_msm2), &
        ' ; grain (mg/kg) = ', MERGE(1000.0_wp * M_grain(i_msm2) / MAX(W_grain, 1.0e-9_wp), 0.0_wp, W_grain > 1.0e-9_wp)
    END DO
  END BLOCK
END IF
! rec 169 zn_uptake_kinetic_basis -- ON-only season summary (one line; the OFF run.log is unchanged)
IF (zn_uptake_kinetic_basis >= 1) THEN
  WRITE (*,'(A,I0,A,I0,A,I0,A,F8.5,A,G14.6)') &
    ' rec 169 ZN SURFACE-CONCENTRATION KINETIC: days on the zinc limb = ', znk_days_active, &
    ' ; days q_diff <= 0 (c_surf = bulk soln) = ', znk_days_qdiff0, &
    ' ; days c_surf <= 0 (active term 0.0 exactly, guard) = ', znk_days_csurf0, &
    ' ; mean saturation c_surf/(Km+c_surf) = ', &
    MERGE(znk_sat_frac_sum / REAL(MAX(znk_days_active - znk_days_csurf0, 1), wp), 0.0_wp, znk_days_active > znk_days_csurf0), &
    ' ; cumulative S_hats before the pool cap (mg/m2) = ', znk_s_hats_cum
END IF
IF (npg_open) CLOSE(unit_npg)         ! rec 121 ON-only diagnostic
IF (lna_open) CLOSE(unit_lna)         ! rec 122 ON-only diagnostic
IF (wfl_open) CLOSE(unit_wfl)         ! rec 148 ON-only diagnostic
IF (tna_open) CLOSE(unit_tna)         ! rec 123 ON-only diagnostic
IF (rfl_open) CLOSE(unit_rfl)         ! rec 124 ON-only diagnostic
IF (kpt_open) CLOSE(unit_kpt)         ! rec 125 ON-only diagnostic
IF (rsf_open) CLOSE(unit_rsf)         ! rec 127 ON-only diagnostic
IF (ncs_open) CLOSE(unit_ncs)         ! rec 128 ON-only diagnostic
IF (cfo_open) CLOSE(unit_cfo)         ! rec 129 ON-only diagnostic
IF (gnm_open) CLOSE(unit_gnm)         ! rec 130 ON-only diagnostic
IF (nud_open) CLOSE(unit_nud)         ! rec 139 ON-only diagnostic
IF (msm_open) CLOSE(unit_msm)         ! rec 140 ON-only diagnostic
IF (msl_diag_open) CLOSE(unit_msl)    ! rec 131 ON-only diagnostic
IF (ont_open) CLOSE(unit_ont)         ! rec 154 ON-only diagnostic
CLOSE(unit_nitrogen)
IF (soilN_layered == 1 .AND. unit_soilN_lay >= 0) CLOSE(unit_soilN_lay)
CLOSE(unit_gas_exch)
CLOSE(unit_ft_outer)
CLOSE(unit_feddes_node)
CLOSE(unit_rad_aero)
CLOSE(unit_carbon_growth)
! Record 64: the per-rank chronology summary - one row per main-stem phytomer
! rank with its tip-appearance and ligulation thermal times on the degree-Cd
! base-0 clock and the model days they fall on, plus the same appearance times
! on the alternative side-by-side chronology. Written once at the end of the
! run; a pure state print.
IF (clr_diag_open) THEN
  BLOCK
    INTEGER :: i_rank
    CLOSE(unit_clr)
    OPEN (newunit=unit_clr,file='output/canopy_rank_chronology.csv',status='replace',action='write')
    ! Record 174: the prediction, its closed-form Haun stage, the target clock
    ! reading and the confirmation are appended as FOUR EXTRA COLUMNS only when the
    ! record is on, so that a mode-0 file keeps its exact column list.
      WRITE (unit_clr,'(A)') 'rank,tau_app_Cd,tau_lig_Cd,day_app,day_lig,tau_app_evers_Cd,' // &
        'tau_lig_evers_Cd,is_flag_leaf,FLN_emergent,FLN_tip,FLN_evers,FLN_evers_tip,' // &
        'tau_boot_Cd,LED_Cd,phyllochron_ref_Cd'
    DO i_rank = 1, CLR_MAXRANK
      IF (clr_tau_lig(i_rank) < 0.0_wp .AND. clr_tau_lig_ev(i_rank) < 0.0_wp) CYCLE
        WRITE (unit_clr,'(I0,2(",",G14.6),2(",",I0),2(",",G14.6),5(",",I0),3(",",G14.6))') &
          i_rank, clr_tau_app(i_rank), clr_tau_lig(i_rank), &
          clr_day_app(i_rank), clr_day_lig(i_rank), &
          clr_tau_app_ev(i_rank), clr_tau_lig_ev(i_rank), &
          MERGE(1, 0, i_rank == clr_FLN), &
          clr_FLN, clr_FLN_tip, clr_FLN_ev, clr_FLN_ev_tip, &
          clr_tau_boot, clr_LED_Cd, clr_phyllochron_Cd
    END DO
    CLOSE(unit_clr)
    clr_diag_open = .FALSE.
  END BLOCK
END IF
! Record 65: the per-rank blade table - one row per main-stem phytomer rank
! with its final length, maximum width and form factor, the potential and
! realized one-sided blade area, the inversed rank the measured panels are
! indexed by, and the potential area the SAME rank would carry on each of the
! six primary tiller classes through the phytomer shifts. The tiller columns
! are POTENTIAL AREAS PER AXIS and are multiplied by NO cohort density -
! densities arrive with record 66.
IF (cba_diag_open) THEN
  BLOCK
    INTEGER  :: i_rank, i_t, nmax_used
    REAL(wp) :: xs, lt, wt
    CLOSE(unit_cba)
    OPEN (newunit=unit_cba,file='output/canopy_blade_area_ranks.csv',status='replace',action='write')
    WRITE (unit_cba,'(A)') 'rank,inversed_rank,L_blade_cm,W_blade_cm,form_factor,A_pot_cm2,A_realized_cm2,' // &
      'tau_app_Cd,tau_lig_Cd,A_pot_T1_cm2,A_pot_T2_cm2,A_pot_T3_cm2,A_pot_T4_cm2,A_pot_T5_cm2,A_pot_T6_cm2,' // &
      'A_pot_T1_sterile_cm2,FLN_emergent'
    nmax_used = MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)
    ! rec 174: the diagnostic table must show the SAME top the model
    ! actually used; clr_FLN_pred is 0 on a mode-0 run, so this is inert there.
    DO i_rank = 1, MIN(MIN(MAX(clr_n_app, clr_n_lig), nmax_used), CLR_MAXRANK)
      WRITE (unit_cba,'(I0,",",I0,7(",",G14.6))',ADVANCE='NO') &
        i_rank, nmax_used - i_rank + 1, cba_L_pot(i_rank), cba_W_pot(i_rank), cba_ft_rank(i_rank), &
        cba_A_pot(i_rank), cba_A(i_rank), clr_tau_app(i_rank), clr_tau_lig(i_rank)
      DO i_t = 1, 6
        ! Abichou (2016) Ch.2 Table 2 Eqn T2.7b: a tiller phytomer of rank x carries the
        ! MAIN STEM dimension at the SHIFTED rank x - s, with the shifts from Ch.1
        ! Table 4 p.71. Below the shift, Eqn T2.7a ramps from the measured constant
        ! first-blade length 7.023 cm; the ramp needs the rank of the first elongated
        ! internode, which the median panel used here does not carry, so that length is
        ! applied flat there. The width follows Eqn T2.10 with a coefficient of 0.877.
        xs = REAL(i_rank, wp) - cba_sTi(i_t)
        IF (xs >= 1.0_wp) THEN
          lt = cba_L_at(xs, nmax_used)
        ELSE
          lt = cba_lb1
        END IF
        wt = cba_kw * cba_W_pot(i_rank)
        WRITE (unit_cba,'(",",G14.6)',ADVANCE='NO') lt * wt * cba_ft_rank(i_rank)
      END DO
      ! Sterile T1 (Eqns T2.13 / T2.16): k_l_st on length, k_w_st on width.
      xs = REAL(i_rank, wp) - cba_sTi(1)
      IF (xs >= 1.0_wp) THEN
        lt = cba_L_at(xs, nmax_used)
      ELSE
        lt = cba_lb1
      END IF
      WRITE (unit_cba,'(",",G14.6,",",I0)') &
        (cba_kl_st * lt) * (cba_kw_st * cba_kw * cba_W_pot(i_rank)) * cba_ft_rank(i_rank), clr_FLN
    END DO
    CLOSE(unit_cba)
    cba_diag_open = .FALSE.
  END BLOCK
END IF

! Record 66: the per-axis-class table - one row per axis class (main stem and
! T1 to T5) with its measured birth stage, the light state in its own window of
! opportunity, the resulting occurrence probability, its birth and final
! densities, its final leaf number (Eqn T3.5) and its per-axis blade area on
! both the fertile and the sterile reduction - then the population summary. No
! carbon or nitrogen ledger reads any of it.
IF (ct_diag_open) THEN
  BLOCK
    INTEGER  :: i_c, nl_used
    REAL(wp) :: a_p, a_r, a_s
    CLOSE(unit_ct)
    OPEN (newunit=unit_ct,file='output/canopy_tiller_classes.csv',status='replace',action='write')
    WRITE (unit_ct,'(A)') 'axis_class,HS_emergence_measured,tau_birth_Cd,Pint_window,p_occurrence,' // &
      'density_born_m2,density_final_m2,Nl_axis,area_axis_pot_cm2,area_axis_real_cm2,area_axis_sterile_cm2'
    nl_used = MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)
    WRITE (unit_ct,'(A,10(",",G14.6))') 'MS', 0.0_wp, ct_tau_birth(0), 0.0_wp, 1.0_wp, &
      ct_dens_born(0), ct_dens(0), ct_nl_axis(0), ct_area_axis(0), ct_area_axis_real(0), 0.0_wp
    DO i_c = 1, CT_MAXTIL
      a_p = ct_area_axis(i_c) ; a_r = ct_area_axis_real(i_c) ; a_s = ct_area_axis_st(i_c)
      WRITE (unit_ct,'("T",I0,10(",",G14.6))') i_c, ct_hs_emerg(i_c), ct_tau_birth(i_c), &
        ct_pint_win(i_c), ct_p_occ(i_c), ct_dens_born(i_c), ct_dens(i_c), &
        ct_nl_tiller(i_c, nl_used), a_p, a_r, a_s
    END DO
    WRITE (unit_ct,'(A)') ''
    WRITE (unit_ct,'(A)') 'summary,value,units'
    WRITE (unit_ct,'("Namax_axes_m2,",G14.6,",axes m-2")')        ct_axes_max
    WRITE (unit_ct,'("Namax_axes_per_plant,",G14.6,",axes")')     MERGE(ct_axes_max/plant_density, 0.0_wp, plant_density > 0.0_wp)
    WRITE (unit_ct,'("axes_live_end_m2,",G14.6,",axes m-2")')     ct_axes_live
    WRITE (unit_ct,'("Ne_predicted_m2,",G14.6,",ears m-2")')      ct_ears
    WRITE (unit_ct,'("ears_per_plant,",G14.6,",ears")')           MERGE(ct_ears/plant_density, 0.0_wp, plant_density > 0.0_wp)
    WRITE (unit_ct,'("tau_start_mortality_Cd,",G14.6,",degCd")')  ct_tau_sm
    WRITE (unit_ct,'("tau_end_mortality_Cd,",G14.6,",degCd")')    ct_tau_em
    WRITE (unit_ct,'("HS_start_mortality,",G14.6,",leaves")')     ct_hs_sm
    WRITE (unit_ct,'("cessation_reason,",I0,",1=Evers light law 2=axis cap 0=never")') ct_cess_why
    WRITE (unit_ct,'("blade_area_index_pot,",G14.6,",m2 m-2")')   ct_bai_pot
    WRITE (unit_ct,'("blade_area_index_real,",G14.6,",m2 m-2")')  ct_bai
    WRITE (unit_ct,'("site_filling_implied,",G14.6,",tillers tiller-1 phyllochron-1")') ct_fs_implied
    CLOSE(unit_ct)
    ct_diag_open = .FALSE.
  END BLOCK
END IF

! Record 67: the per-axis senescence table - one row per axis class with the
! green-leaf-number law parameters actually used, the three phase-transition
! times taken from the chronology, the axis final green leaf number and
! senesced leaf number, and its blade area split between green and STANDING
! DEAD - then the population summary and the two comparison quantities (the
! Haun stage at which senescence actually started, and the green leaf number at
! flag ligulation). No carbon or nitrogen ledger reads any of it.
IF (cls_diag_open) THEN
  BLOCK
    INTEGER  :: i_c, nl_used
    CLOSE(unit_cls)
    OPEN (newunit=unit_cls,file='output/canopy_senescence_axes.csv',status='replace',action='write')
    WRITE (unit_cls,'(A)') 'axis_class,Nl_axis,HS_axis,Gl_axis,Sl_axis,phase,' // &
      'tau1_onset_Cd,tau2_Nl_minus_4_Cd,tau3_flag_lig_Cd,' // &
      'area_green_cm2,area_dead_cm2,density_m2'
    nl_used = MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)
    WRITE (unit_cls,'(A,4(",",G14.6),",",I0,6(",",G14.6))') 'MS', REAL(nl_used, wp), &
      ct_hs_axis(0, nl_used), cls_gl(0), cls_sl(0), cls_phase(0), &
      cls_tau1(0), cls_tau2(0), cls_tau3(0), &
      cls_area_green(0), cls_area_dead(0), ct_dens(0)
    DO i_c = 1, CT_MAXTIL
      WRITE (unit_cls,'("T",I0,4(",",G14.6),",",I0,6(",",G14.6))') i_c, &
        ct_nl_tiller(i_c, nl_used), ct_hs_axis(i_c, nl_used), cls_gl(i_c), cls_sl(i_c), &
        cls_phase(i_c), cls_tau1(i_c), cls_tau2(i_c), cls_tau3(i_c), &
        cls_area_green(i_c), cls_area_dead(i_c), ct_dens(i_c)
    END DO
    WRITE (unit_cls,'(A)') ''
    WRITE (unit_cls,'(A)') 'summary,value,units'
    WRITE (unit_cls,'("green_area_index_end,",G14.6,",m2 m-2")')        cls_gai
    WRITE (unit_cls,'("standing_dead_instantaneous_end,",G14.6,",m2 m-2")') cls_dai
    WRITE (unit_cls,'("standing_dead_pool_end,",G14.6,",m2 m-2")')      cls_dai_pool
    WRITE (unit_cls,'("standing_dead_culled_end,",G14.6,",m2 m-2")')    cls_dai_culled
    WRITE (unit_cls,'("standing_dead_total_end,",G14.6,",m2 m-2")')     cls_dai_pool + cls_dai_culled
    WRITE (unit_cls,'("standing_dead_max,",G14.6,",m2 m-2")')           cls_dai_max
    WRITE (unit_cls,'("instantaneous_dip_days,",I0,",days")')           cls_mono_viol
    WRITE (unit_cls,'("instantaneous_dip_max,",G14.6,",m2 m-2")')       cls_dai_dip_max
    WRITE (unit_cls,'("HS_at_senescence_onset,",G14.6,",leaves")')      cls_hs_onset
    WRITE (unit_cls,'("Gl_MS_at_flag_ligulation,",G14.6,",leaves")')    cls_gl_at_flag
    WRITE (unit_cls,'("g1_g2_g3_used,",G14.6,",",G14.6,",",G14.6)')     cls_g1, cls_g2, cls_g3
    WRITE (unit_cls,'("delta_tiller_used,",G14.6,",-")')                cls_delta
    CLOSE(unit_cls)
    cls_diag_open = .FALSE.
  END BLOCK
END IF

! Record 68: the end-of-run per-rank and per-axis non-lamina dump, plus the
! two-basis summary. Written only when the record is on.
IF (cse_diag_open) THEN
  BLOCK
    INTEGER :: i_c, i_x, nl_used
    REAL(wp) :: ls_x, ds_x, li_x, lsb_x
    CLOSE(unit_cse)
    OPEN (newunit=unit_cse,file='output/canopy_sheath_ear_ranks.csv',status='replace',action='write')
    nl_used = MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)
    WRITE (unit_cse,'(A)') 'rank_abs,rank_inversed,sheath_length_cm,sheath_diam_cm,' // &
      'sheath_silhouette_cm2,internode_length_cm,sheath_below_cm,internode_exposed_cm'
    DO i_x = 1, MIN(nl_used, CLR_MAXRANK)
      ls_x  = cse_Ls_at(REAL(i_x, wp), nl_used)
      ds_x  = cse_Ds_ms(REAL(i_x, wp), nl_used)
      li_x  = cse_Li_at(REAL(i_x, wp), nl_used)
      lsb_x = cse_Ls_at(REAL(i_x - 1, wp), nl_used)
      WRITE (unit_cse,'(I0,",",I0,6(",",G14.6))') i_x, nl_used - i_x + 1, &
        ls_x, ds_x, ls_x * ds_x, li_x, lsb_x, MAX(0.0_wp, li_x - lsb_x)
    END DO
    WRITE (unit_cse,'(A)') ''
    WRITE (unit_cse,'(A)') 'axis_class,Sl_axis,sheath_green_cm2,sheath_dead_cm2,' // &
      'stem_exposed_cm2,ear_ped_cm2,ko_diameter_factor,density_m2'
    WRITE (unit_cse,'(A,7(",",G14.6))') 'MS', cls_sl(0), cse_sh_green(0), &
      cse_sh_dead(0), cse_st_axis(0), cse_ep_axis(0), 1.0_wp, ct_dens(0)
    DO i_c = 1, CT_MAXTIL
      WRITE (unit_cse,'("T",I0,7(",",G14.6))') i_c, cls_sl(i_c), cse_sh_green(i_c), &
        cse_sh_dead(i_c), cse_st_axis(i_c), cse_ep_axis(i_c), cse_ko_tiller(i_c), ct_dens(i_c)
    END DO
    WRITE (unit_cse,'(A)') ''
    WRITE (unit_cse,'(A)') 'summary,value,units'
    WRITE (unit_cse,'("sheath_green_index_end,",G14.6,",m2 m-2")')      cse_sai
    WRITE (unit_cse,'("sheath_dead_pool_end,",G14.6,",m2 m-2")')        cse_sdi_pool
    WRITE (unit_cse,'("sheath_dead_instant_end,",G14.6,",m2 m-2")')     cse_sdi
    WRITE (unit_cse,'("sheath_dip_days,",I0,",days")')                  cse_mono_viol
    WRITE (unit_cse,'("sheath_dip_max,",G14.6,",m2 m-2")')              cse_sdi_dip_max
    WRITE (unit_cse,'("sheath_green_n1_arm_end,",G14.6,",m2 m-2")')     cse_sai_n1
    WRITE (unit_cse,'("sheath_green_n3_arm_end,",G14.6,",m2 m-2")')     cse_sai_n3
    WRITE (unit_cse,'("sheath_green_nolag_arm_end,",G14.6,",m2 m-2")')  cse_sai_full
    WRITE (unit_cse,'("stem_exposed_index_end,",G14.6,",m2 m-2")')      cse_iai
    WRITE (unit_cse,'("ear_ped_index_end,",G14.6,",m2 m-2")')           cse_eai
    WRITE (unit_cse,'("ear_ped_index_bracket,",G14.6,",",G14.6)')       cse_eai_lo, cse_eai_hi
    WRITE (unit_cse,'("ear_only_index_measured_Zhang2023,",G14.6,",m2 m-2")') cse_eai_meas
    WRITE (unit_cse,'("nonlamina_over_lamina_end,",G14.6,",-")')        cse_ratio
    WRITE (unit_cse,'("nonlamina_over_lamina_at_flowering,",G14.6,",-")') cse_ratio_flower
    WRITE (unit_cse,'("PAIg_green_at_flowering,",G14.6,",m2 m-2")')     cse_paig_flower
    WRITE (unit_cse,'("live_axes_at_flowering,",G14.6,",axes m-2")')    cse_axis_flower
    WRITE (unit_cse,'("PAIg_cm2_per_axis_at_flowering,",G14.6,",cm2 axis-1")') &
      MERGE(cse_paig_flower / MAX(cse_axis_flower, 1.0e-9_wp) * 1.0e4_wp, 0.0_wp, &
            cse_axis_flower > 0.0_wp)
    WRITE (unit_cse,'("PAI_LI2200C_basis_max,",G14.6,",m2 m-2")')       cse_pai_li_max
    WRITE (unit_cse,'("A_ep_used,",G14.6,",cm2 per fertile culm")')     cse_a_ep
    WRITE (unit_cse,'("ds2_used,",G14.6,",cm")')                        cse_ds2
    WRITE (unit_cse,'("persistence_rule_n_plus,",I0,",ranks")')         cse_persist_n
    CLOSE(unit_cse)
    cse_diag_open = .FALSE.
  END BLOCK
END IF

! Records 82, 83 and 84: the end-of-run cohort table - the supply-ratio gate
! per cohort (window days, summed surplus, organ-law demand, ranks counted,
! ratio, births, final density), the envelope counters and the ear envelope.
IF (d2_diag_open) THEN
  BLOCK
    INTEGER :: i_c
    CLOSE(unit_d2)
    OPEN (newunit=unit_d2,file='output/d2_cohort_table.csv',status='replace',action='write')
    WRITE (unit_d2,'(A)') 'axis_class,HS_emergence_measured,tau_birth_Cd,window_days,Pint_window,p_occurrence_TIMING_ONLY,' // &
      'sum_Csur_window_gDM,D_cohort_gDM,ranks_in_D_cohort,supply_ratio,density_born_m2,density_final_m2,Nl_axis'
    WRITE (unit_d2,'(A,3(",",G14.6),",",A,",",A,2(",",G14.6),",",I0,4(",",G14.6))') 'MS', 0.0_wp, ct_tau_birth(0), 0.0_wp, &
      'NA', 'NA', 0.0_wp, 0.0_wp, 0, 1.0_wp, ct_dens_born(0), ct_dens(0), ct_nl_axis(0)
    DO i_c = 1, CT_MAXTIL
      WRITE (unit_d2,'("T",I0,5(",",G14.6),2(",",G14.6),",",I0,4(",",G14.6))') i_c, ct_hs_emerg(i_c), ct_tau_birth(i_c), &
        tcg_nwin(i_c), ct_pint_win(i_c), ct_p_occ(i_c), tcg_sup(i_c), tcg_dem(i_c), tcg_nrank(i_c), &
        tcg_ratio(i_c), ct_dens_born(i_c), ct_dens(i_c), ct_nl_axis(i_c)
    END DO
    WRITE (unit_d2,'(A)') ''
    WRITE (unit_d2,'(A)') 'summary,value,units'
    WRITE (unit_d2,'("rec82_tiller_carbon_gate,",I0,",-")')     tiller_carbon_gate
    WRITE (unit_d2,'("rec83_leaf_expansion_water,",I0,",-")')   leaf_expansion_water
    WRITE (unit_d2,'("rec84_lai_from_organs,",I0,",-")')        lai_from_organs
    WRITE (unit_d2,'("Namax_axes_m2,",G14.6,",axes m-2")')       ct_axes_max
    WRITE (unit_d2,'("axes_live_end_m2,",G14.6,",axes m-2")')    ct_axes_live
    WRITE (unit_d2,'("axes_born_m2,",G14.6,",axes m-2")')        ct_axes_born
    WRITE (unit_d2,'("env600_days,",I0,",days live axes > 600 (Abichou Ch.3 Fig.3 lower bound) REPORTED")') tcg_env600_days
    WRITE (unit_d2,'("env800_days,",I0,",days live axes > 800 (Abichou Ch.3 Fig.3 upper bound) REPORTED")') tcg_env800_days
    WRITE (unit_d2,'("Ne_envelope_lo_m2,",G14.6,",ears m-2 IF Evers die fraction 0.67 on Namax - NOT APPLIED")') &
      tcg_ne_env_lo
    WRITE (unit_d2,'("Ne_envelope_hi_m2,",G14.6,",ears m-2 IF Evers die fraction 0.57 on Namax - NOT APPLIED")') &
      tcg_ne_env_hi
    WRITE (unit_d2,'("productive_fraction,",A)') 'NOT DETERMINED (D2-4 blocked: per-axis structural demand unmeasured)'
    WRITE (unit_d2,'("cessation_reason_end,",I0,",0 open / 1 Evers light law / 2 cap (never under rec 82)")') ct_cess_why
    WRITE (unit_d2,'("LMA_used,",G14.6,",g DM m-2 (Sales 2022 Table 1 anthesis mean)")') csl_LMA
    WRITE (unit_d2,'("Acker_wheat_200ppm_line,",A)') 'LER = 1.6 + 3.4*SRWC (Fig. 4B) - bracket'
    WRITE (unit_d2,'("Acker_wheat_800ppm_line,",A)') 'LER = 0.8 + 4.7*SRWC (Fig. 4B) - bracket'
    WRITE (unit_d2,'("f_w_applied,",A)') 'mean of the two lines normalised at SRWC = 1: (1.2 + 4.05*SRWC)/5.25'
    WRITE (unit_d2,'("SRWC_basis,",A)') 'root-fraction-weighted theta/theta_33 (Acker Eq. 1: water / field capacity);' // &
      ' Se alternative reported'
    WRITE (unit_d2,'("SLA_out_end_m2_kg,",G14.6,",m2 kg-1")')    lfo_sla_out
    ! (rec 84 leaf-death coupling) season totals; all zero when rec 84 is OFF.
    IF (lai_from_organs == 1) THEN
      WRITE (unit_d2,'("leaf_death_clock,",A)') 'organ law (Abichou M-law) carries W_green and N_leaf out with the dead area (v14)'
    ELSE
      WRITE (unit_d2,'("leaf_death_clock,",A)') 'legacy Yin/Teh biomass law (rec 84 OFF)'
    END IF
    WRITE (unit_d2,'("sen_dm_cum_gDM,",G14.6,",g DM m-2 lamina shed by the organ clock")')            lfo_sen_dm_cum
    WRITE (unit_d2,'("sen_N_out_cum_gN,",G14.6,",g N m-2 out of N_leaf with the dead area")')         lfo_sen_N_out_cum
    WRITE (unit_d2,'("sen_N_dead_cum_gN,",G14.6,",g N m-2 structural (SLN_min x area) to the dead leaf")') lfo_sen_N_dead_cum
    WRITE (unit_d2,'("sen_N_resorb_cum_gN,",G14.6,",g N m-2 resorbed to N_mobile_veg")')               lfo_sen_N_resorb_cum
    ! (rec 99): the fertile ear population and the evidence it rests on.
    ! Written ONLY when the record is ON, so the OFF file is unchanged.
    IF (tiller_survival_mech >= 1) THEN
      WRITE (unit_d2,'(A)') ''
      WRITE (unit_d2,'(A)') 'rec99_tiller_survival,value,units'
      WRITE (unit_d2,'("rec99_tiller_survival_mech,",I0,",-")')                 tiller_survival_mech
      WRITE (unit_d2,'("tau_sm_Cd,",G14.6,",degCd start of mortality phase (Abichou Eqn T4.2)")')   ct_tau_sm
      WRITE (unit_d2,'("tau_em_Cd,",G14.6,",degCd end = flag-leaf ligulation of the main stem (Eqn T4.3)")') ct_tau_em
      WRITE (unit_d2,'("mortality_window_days,",I0,",days")')                   MERGE(tsv_days, tca_days, tiller_survival_mech == 1)
      WRITE (unit_d2,'("ears_predicted_m2,",G14.6,",ears m-2 PREDICTED, spike-mass-free (surviving axes at tau_em)")') ct_ears
      WRITE (unit_d2,'("axes_max_m2,",G14.6,",axes m-2 before mortality (Namax)")')                  ct_axes_max
      WRITE (unit_d2,'("plants_m2,",G14.6,",plants m-2 = the Ne floor, main stem never culled")')    ct_dens_born(0)
      WRITE (unit_d2,'("die_fraction_realised,",G14.6,",- (Namax - Ne)/(Namax - Np), EMERGENT not prescribed")') &
        MERGE((ct_axes_max - ct_ears) / MAX(1.0e-12_wp, ct_axes_max - ct_dens_born(0)), 0.0_wp, &
              ct_axes_max > ct_dens_born(0))
      IF (tiller_survival_mech == 1) THEN
      DO i_c = 0, CT_MAXTIL
        WRITE (unit_d2,'("class,",I0,4(",",G14.6))') i_c, tsv_kill(i_c), tsv_ilight(i_c), tsv_sup(i_c), tsv_dem(i_c)
      END DO
      WRITE (unit_d2,'(A)') 'class_columns,axis_class,kill_fraction,light_fraction_Beer,' // &
        'cum_supply_gDM_per_axis,cum_potential_lamina_demand_gDM_per_axis'
      ELSE
      ! rec 99 = 2: mode, the unaffordable density, the skipped days and
      ! the per-class ledgers in g CH2O (NOT g DM: the growth currency).
      WRITE (unit_d2,'("mode,",I0,",- (2 = tiller carbon affordability, SHARED)")') tiller_survival_mech
      WRITE (unit_d2,'("dead_axes_m2,",G14.6,",axes m-2 dead = SUM_c dead_c (each cohort monotone; youngest first)")') &
        tca_dead
      WRITE (unit_d2,'("unaff_target_m2,",G14.6,",axes m-2 monotone population target MAX_t SUM dens_pre (1 - phi)")') &
        tca_unaff
      WRITE (unit_d2,'("days_skip,",I0,",window days skipped because L_tot = 0")') tca_days_skip
      WRITE (unit_d2,'("frozen,",I0,",1 = tau_em reached and the population final")') MERGE(1, 0, tca_frozen)
      DO i_c = 0, CT_MAXTIL
        WRITE (unit_d2,'("class,",I0,7(",",G14.6))') i_c, tca_phi(i_c), tca_share(i_c), tca_S(i_c), tca_D(i_c), &
          tca_dens_pre(i_c), ct_dens(i_c), tca_dead_c(i_c)
      END DO
      WRITE (unit_d2,'(A)') 'class_columns,axis_class,affordability_phi,light_share_last_day,' // &
        'cum_supply_gCH2O_per_axis,cum_lamina_demand_gCH2O_per_axis,dens_entering_m2,dens_final_m2,dead_m2'
      ! rec 99 = 4: the timing repair's endpoint quantities (mode-2 lines above reused).
      IF (tiller_survival_mech == 4) THEN
        WRITE (unit_d2,'("mode4_label,",A)') 'survival timed by T4.4 cap on the matured ledger (SHARED);' // &
          ' the mode line above carries the mode-2 label text'
        WRITE (unit_d2,'("U_endpoint_m2,",G14.6,",axes m-2 U = SUM dens_pre (1 - phi) on the last stepped window day")') tca4_U
        WRITE (unit_d2,'("dead_target_m2,",G14.6,",axes m-2 MAX_t F(tau_n) U(t) incl. the endpoint (F = 1)")') tca_unaff
        WRITE (unit_d2,'("premature_m2,",G14.6,",axes m-2 dead - U_endpoint (deaths taken on an immature U; PK gate)")') &
          tca_dead - tca4_U
        WRITE (unit_d2,'("F_endpoint,",G14.6,",- F on the endpoint day (:= 1 by Eqn T4.3 / Fig. 2 phi4 plateau)")') tca4_F
        WRITE (unit_d2,'("Ne_source_mode4,",A)') 'Eqn T4.4 used as TIMING SHAPE only (theta 3.32, delta 1.77, doc p.36);' // &
          ' N_e EMERGENT; the Ne_source line below is the mode-2 wording'
      END IF
      END IF
      WRITE (unit_d2,'("Ne_source,",A)') 'PREDICTED by rec 99; Abichou 2016 doc p.45 lists Ne as an INPUT' // &
        ' VARIABLE of his tillering model, so Eqn T4.4 is NOT used and no observed ear count enters'
      WRITE (unit_d2,'("cullgreen_seen_m2_m2,",G14.6,",m2 m-2 green blade area the rec-99 cull removed")') &
        tsv_cullgreen_seen
      WRITE (unit_d2,'("cullgreen_charged_m2_m2,",G14.6,",m2 m-2 the same area charged to the leaf-death coupling")') &
        tsv_cullgreen_charged
      WRITE (unit_d2,'("cullgreen_closure_resid_m2_m2,",G14.6,",m2 m-2 |seen - charged|; GATE D closure for the rec-99 path")') &
        ABS(tsv_cullgreen_seen - tsv_cullgreen_charged)
      WRITE (unit_d2,'("new_constants_introduced,",I0,",count")') 0
    END IF
    CLOSE(unit_d2)
    d2_diag_open = .FALSE.
  END BLOCK
END IF
IF (tca_diag_open) THEN
  CLOSE(unit_tca)
  tca_diag_open = .FALSE.
END IF
IF (tca4_diag_open) THEN
  CLOSE(unit_tca4)
  tca4_diag_open = .FALSE.
END IF

! ==================================================================== !
! Record 89: the grain-number chain at the end of the run, written only when
! the record is on. Every row states its source or says that the quantity is
! not determined; nothing here is fitted to any observation.
! ==================================================================== !
IF (grain_number_mech >= 1) THEN
  BLOCK
    INTEGER :: unit_gnx, i_c
    OPEN (newunit=unit_gnx,file='output/grain_number_chain.csv',status='replace',action='write')
    ! rec 90 = 0 writes EXACTLY the item-A+D file, so that arm stays;
    ! the two chain columns appear only when the chain is on.
    IF (floret_survival_mech == 4) THEN
      ! rec 90 = 4: the mode-1 columns + the ear-bearing fraction of the class (ears of the class = density x phi_sp)
      WRITE (unit_gnx,'(A)') 'axis_class,density_m2,spike_DW_per_axis_g,competent_florets_per_spike,' // &
        'kernels_per_spike,bears_ear,spike_DM_pool_g_m2,floret_survival,spike_DM_pool_pre_mortality_g_m2,ear_fraction_phi_sp'
    ELSE IF (floret_survival_mech >= 1) THEN
      WRITE (unit_gnx,'(A)') 'axis_class,density_m2,spike_DW_per_axis_g,competent_florets_per_spike,' // &
        'kernels_per_spike,bears_ear,spike_DM_pool_g_m2,floret_survival,spike_DM_pool_pre_mortality_g_m2'
    ELSE
      WRITE (unit_gnx,'(A)') 'axis_class,density_m2,spike_DW_per_axis_g,competent_florets_per_spike,' // &
        'kernels_per_spike,bears_ear,spike_DM_pool_g_m2'
    END IF
    DO i_c = 0, CT_MAXTIL
      IF (floret_survival_mech == 4) THEN
        IF (i_c == 0) THEN
          WRITE (unit_gnx,'(A,4(",",G14.6),",",A,4(",",G14.6))') 'MS', gnx_dens(0), gnx_wsp(0), gnx_nfl(0), gnx_kn(0), &
            MERGE('YES','NO ',gnx_is_ear(0)), gnx_pool(0), gnx_surv(0), gnx_pool_pre(0), gnx_phi_sp(0)
        ELSE
          WRITE (unit_gnx,'("T",I0,4(",",G14.6),",",A,4(",",G14.6))') i_c, gnx_dens(i_c), gnx_wsp(i_c), gnx_nfl(i_c), &
            gnx_kn(i_c), MERGE('YES','NO ',gnx_is_ear(i_c)), gnx_pool(i_c), gnx_surv(i_c), gnx_pool_pre(i_c), gnx_phi_sp(i_c)
        END IF
      ELSE IF (floret_survival_mech >= 1) THEN
        IF (i_c == 0) THEN
          WRITE (unit_gnx,'(A,4(",",G14.6),",",A,3(",",G14.6))') 'MS', gnx_dens(0), gnx_wsp(0), gnx_nfl(0), gnx_kn(0), &
            MERGE('YES','NO ',gnx_is_ear(0)), gnx_pool(0), gnx_surv(0), gnx_pool_pre(0)
        ELSE
          WRITE (unit_gnx,'("T",I0,4(",",G14.6),",",A,3(",",G14.6))') i_c, gnx_dens(i_c), gnx_wsp(i_c), gnx_nfl(i_c), &
            gnx_kn(i_c), MERGE('YES','NO ',gnx_is_ear(i_c)), gnx_pool(i_c), gnx_surv(i_c), gnx_pool_pre(i_c)
        END IF
      ELSE
        IF (i_c == 0) THEN
          WRITE (unit_gnx,'(A,4(",",G14.6),",",A,",",G14.6)') 'MS', gnx_dens(0), gnx_wsp(0), gnx_nfl(0), gnx_kn(0), &
            MERGE('YES','NO ',gnx_is_ear(0)), gnx_pool(0)
        ELSE
          WRITE (unit_gnx,'("T",I0,4(",",G14.6),",",A,",",G14.6)') i_c, gnx_dens(i_c), gnx_wsp(i_c), gnx_nfl(i_c), &
            gnx_kn(i_c), MERGE('YES','NO ',gnx_is_ear(i_c)), gnx_pool(i_c)
        END IF
      END IF
    END DO
    WRITE (unit_gnx,'(A)') ''
    WRITE (unit_gnx,'(A)') 'quantity,value,units_or_source'
    WRITE (unit_gnx,'("rec89_grain_number_mech,",I0,",0 OFF / 1 A+D / 2 D only / 3 A only")') grain_number_mech
    WRITE (unit_gnx,'("chain_fired,",A,",-")') MERGE('YES','NO ', gnx_set)
    WRITE (unit_gnx,'(A)') '# --- the spike growth period (Fischer 1985 p.449: penultimate leaf appearance -> anthesis) ---'
    WRITE (unit_gnx,'("window_open_model_day,",I0,",day the MAIN STEM reached Haun FLN-1")') gnx_win_open_day
    WRITE (unit_gnx,'("window_days,",I0,",d (Fischer: the 30 days or so preceding anthesis, p.447)")') gnx_win_days
    WRITE (unit_gnx,'("window_degCd_above_4.5,",G14.6,",screen air; FALSIFIER vs Fischer p.449 ~300")') &
      gnx_win_tt45
    WRITE (unit_gnx,'("window_open_rule,",A)') &
      'projected remaining degCd above 4.5 C <= 300 (Fischer p.449 thermal form)'
    WRITE (unit_gnx,'("tdu_remaining_at_open,",G14.6,",phenology units still to anthesis")') gnx_tdu_rem_open
    WRITE (unit_gnx,'("tt_ratio_degCd_per_tdu,",G14.6,",season-to-date unit conversion; DECLARED")') gnx_tt_ratio
    WRITE (unit_gnx,'("penultimate_leaf_model_day,",I0,",Fischer DEVELOPMENTAL form; CROSS-CHECK ONLY, read by nothing")') &
      gnx_penult_day
    WRITE (unit_gnx,'(A)') '# --- the partitioning coefficient (Fischer 1985 p.452: 0.44 of the ABOVE-GROUND increment) ---'
    WRITE (unit_gnx,'("window_AGB_increment_g_m2,",G14.6,",SUM dm_gain*(1-f_root) in window")') gnx_win_dagb
    WRITE (unit_gnx,'("window_spike_increment_g_m2,",G14.6,",SUM dm_gain*f_storage in window")') gnx_win_dspike
    WRITE (unit_gnx,'("window_spike_share_realised,",G14.6,",TARGET Fischer 0.44 p.452 = ITEM-D GATE")') gnx_share_real
    WRITE (unit_gnx,'("spike_share_clamp_days,",I0,",days lamina left < 0.44*(1-f_root) room")') gnx_clamp_days
    WRITE (unit_gnx,'("spike_share_clamp_max,",G14.6,",share units")') gnx_clamp_max
    WRITE (unit_gnx,'(A)') '# --- the chain (Fischer 1985 p.449 florets, p.447 kernels) ---'
    WRITE (unit_gnx,'("spike_DW_at_anthesis_g_m2,",G14.6,",= W_chaff when the number was set")') gnx_spike_dw_anth
    WRITE (unit_gnx,'("florets_per_g_spike_DW,",G14.6,",Fischer p.449 1 floret / 10 mg; Tab 4 102.3")') &
      1.0_wp/gnx_m_floret
    WRITE (unit_gnx,'("kernels_per_competent_floret,",G14.6,",78 g-1 p.447 / 100 g-1 p.449")') &
      gnx_kern_per_fl
    WRITE (unit_gnx,'("competent_florets_m2,",G14.6,",florets m-2 PREDICTED")') gnx_florets_m2
    IF (floret_survival_mech == 4) THEN
      WRITE (unit_gnx,'("ears_m2,",G14.6,",PREDICTED: SUM dens x phi_sp (spike affordability, no bar)")') &
        gnx_ears_m2
    ELSE IF (floret_survival_mech >= 1) THEN
      WRITE (unit_gnx,'("ears_m2,",G14.6,",PREDICTED: axes carrying >= 1 fertile floret PER SPIKELET")') &
        gnx_ears_m2
    ELSE
      WRITE (unit_gnx,'("ears_m2,",G14.6,",PREDICTED: axis class with >= 1 kernel per spike")') &
        gnx_ears_m2
    END IF
    WRITE (unit_gnx,'("grains_per_ear,",G14.6,",kernels ear-1 PREDICTED")') gnx_grains_per_ear
    WRITE (unit_gnx,'("grain_number_m2,",G14.6,",kernels m-2 PREDICTED")') kernel_number_m2
    IF (floret_survival_mech >= 1) THEN
      WRITE (unit_gnx,'("kernels_per_g_spike_DW_realised,",G14.6,",THE HEADLINE TEST: must VARY across decks;' // &
        ' measured 52.4-82.9 Rivera-Amado 2019 FCR 240:154, 73-94 Fischer 1985 p.454, 78+-2 p.447")') gnx_fe_realised
      WRITE (unit_gnx,'("kernels_per_spikelet_realised,",G14.6,",FALSIFIER vs Fischer 1985 p.450 measured 1.98")') &
        MERGE(gnx_grains_per_ear/gnx_spikelets, 0.0_wp, gnx_spikelets > 0.0_wp)
      WRITE (unit_gnx,'("ears_per_plant,",G14.6,",FALSIFIER vs Abichou 2016 Ch.1 sec.3.6 measured 1.4-5.0")') &
        MERGE(gnx_ears_m2/plant_density, 0.0_wp, plant_density > 0.0_wp)
      WRITE (unit_gnx,'("spike_DW_per_ear_g,",G14.6,",FALSIFIER vs Fischer Tab 4 0.29-0.32 g, Rivera-Amado 0.56-1.33 g")') &
        MERGE(gnx_spike_dw_anth/gnx_ears_m2, 0.0_wp, gnx_ears_m2 > 0.0_wp)
    END IF
    ! ---- Record 98: the spike sink and the composition split, printed only when
    ! the record is on.
    WRITE (unit_gnx,'(A)') '# --- carbon bookkeeping (the per-axis split creates and destroys nothing) ---'
    WRITE (unit_gnx,'("pool_rescale_to_W_chaff,",G14.6,",1.0 = pools already summed to W_chaff")') &
      gnx_pool_scale
    WRITE (unit_gnx,'("spike_DM_unshared_g_m2,",G14.6,",to main stem on days no axis window open")') &
      gnx_unshared
    WRITE (unit_gnx,'("spike_DM_in_sterile_classes_g_m2,",G14.6,",no-ear classes; STAYS in W_chaff")') &
      gnx_sterile_C
    WRITE (unit_gnx,'(A)') '# --- DECLARED: not a measurement ---'
    WRITE (unit_gnx,'("per_axis_share_rule,",A)') &
      'green lamina area of the axis class (ct_dens*ct_area_axis_real), from the CROP' // &
      ' window onward - DECLARED CONSTRUCTION; no measured per-axis sharing law is on' // &
      ' this machine. A per-axis sub-window was built first and REMOVED: see growth_module'
    IF (floret_survival_mech >= 1) THEN
      IF (floret_survival_mech == 4) THEN
        WRITE (unit_gnx,'("rec90_floret_survival_mech,",I0,",0 OFF / 1 chain / 2 survival only / 3 hyperbola arm' // &
          ' / 4 spike-affordability ear test")') floret_survival_mech
      ELSE
        WRITE (unit_gnx,'("rec90_floret_survival_mech,",I0,",0 OFF / 1 chain / 2 survival only / 3 hyperbola arm")') &
          floret_survival_mech
      END IF
      WRITE (unit_gnx,'("spikelets_per_spike,",G14.6,",Fischer 1985 p.450 verbatim, cv. Yecora 70 CIANO' // &
        ' (436 spikes/m2, 19.4 spikelets/spike, 1.98 kernels/spikelet)")') gnx_spikelets
      WRITE (unit_gnx,'("floret_primordia_per_spikelet,",G14.6,",Guo & Schnurbusch 2015 JXB 66:5945 DOI' // &
        ' 10.1093/jxb/erv303: around 10, range 9-12; Guo 2016 8-11; Backhaus 2023 mean 9.13")') gnx_mfs
      WRITE (unit_gnx,'("floret_primordia_per_spike,",G14.6,",= spikelets x primordia per spikelet")') gnx_ninit
      WRITE (unit_gnx,'("floret_survival_law,",A)') 'S = 1 - exp(-W_spike_per_axis / (194 x 10 mg)).' // &
        ' DECLARED CONSTRUCTION with NO constant of its own: pinned at W->0 on Fischer''s measured' // &
        ' proportionality (1 floret per 10 mg) and at W->inf on the measured primordia ceiling' // &
        ' (Guo & Schnurbusch 2015). Backhaus 2023 states survival is a saturation function; every' // &
        ' named alternative (cumulative normal / Hill / logistic) carries a shape parameter, this does not.'
      IF (floret_survival_mech == 4) THEN
      WRITE (unit_gnx,'("fertile_tiller_test,",A)') 'rec 90 = 4: a class bears ears on the FRACTION' // &
        ' phi_sp = MIN(1, W_axis / W_pot) of its axes, W_axis = spike DM per axis actually accrued (gnx_pool_pre /' // &
        ' ct_dens, SHARED gnx_pool_pre) and W_pot = bgn_r_spike x gnx_win_tt45 (SHARED bgn_r_spike = 0.323 g / 300' // &
        ' degCd, Fischer 1985 Tab 4 p.451 x p.449; gnx_win_tt45 = the realised window degCd above 4.5 C); the rest' // &
        ' STAND without a spike (Fischer 1985 p.452 verbatim "higher order tillers which do not bear spikes in the' // &
        ' control crop yet amount to 5-10% of total anthesis dry weight") and their pool returns to the ear-bearers' // &
        ' pro rata (carbon-neutral); main stem never earless (Abichou 2016 N_e floor); NO bar, NO new constant; the' // &
        ' ONE affordability law of rec 99 (phi = MIN(1, supply/demand), SHARED) applied to the spike sink' // &
        ' in its own sourced window. Abichou 2016 Ch.1 sec.1.3 doc pp.20-21: one dichotomy (die / complete the' // &
        ' cycle and produce an ear), one cause (competition), one order (last emerged), NO criterion at the ear.' // &
        ' The ROLE of Fischer''s control spike as the per-axis potential is DECLARED (the rec 98 [B-GN3]' // &
        ' declaration). Youngest-first is EMERGENT. Code: growth_module.f90 gnx_set_grain_number, rec 90 = 4 branch.'
      ELSE
      WRITE (unit_gnx,'("fertile_tiller_test,",A)') 'an axis bears an ear iff its spike carries >= 1' // &
        ' fertile floret PER SPIKELET. THE AXIS-LEVEL RULE IS CITED, THE BAR IS A DECLARED MODEL' // &
        ' CONSTRUCTION. Rule: Abichou 2016 (PhD thesis, NNT 2016SACLA022, HAL tel-02316154) Ch.1' // &
        ' sec.1.3 p.20 verbatim "a phase of mortality during which a proportion of tillers die while' // &
        ' the others complete their cycle and produce an ear", pp.20-21 the tiller that dies is the' // &
        ' "last emerged" or "less developed" and the number completing the cycle "is highly influenced' // &
        ' by availability of light, water, and nutriments"; field corroboration Fischer 1985 p.452' // &
        ' "higher order tillers which do not bear spikes in the control crop yet amount to 5-10% of' // &
        ' total anthesis dry weight ... likely to be suppressed sooner in shaded crops". Bar: one' // &
        ' fertile floret per spikelet is THIS MODEL''S carbon criterion for an axis with no' // &
        ' reproductive sink, NOT a cited law; the choice of one rather than two is the implementer''s' // &
        ' and rec 90 = 2 prices it. Backhaus et al. 2023 JXB 74:5088 is an ANALOGY ONLY: its unit is' // &
        ' the rudimentary basal SPIKELET of an otherwise fertile spike, it removes no spike and' // &
        ' reports no axis-level analysis (independent audit 59a4951 sec.3). Youngest culled first' // &
        ' (Abichou 2016 Ch.1 sec.2.5 p.29, verbatim); the main stem is never culled.'
      END IF
      IF (floret_survival_mech == 4) THEN
        WRITE (unit_gnx,'("axis_classes_culled,",I0,",classes with ear-bearing fraction phi_sp < 1 (part reduction,' // &
          ' no class removed)")') gnx_axes_culled
      ELSE
        WRITE (unit_gnx,'("axis_classes_culled,",I0,",classes that failed the fertile-tiller test")') gnx_axes_culled
      END IF
      WRITE (unit_gnx,'("cull_passes,",I0,",monotone loop passes")') gnx_cull_passes
      WRITE (unit_gnx,'("spike_DM_returned_to_survivors_g_m2,",G14.6,",carbon-NEUTRAL: SUM_c pool(c) unchanged")') &
        gnx_pool_moved
      IF (floret_survival_mech == 4) THEN
        ! rec 90 = 4: the three fates and the two identities, in double precision from STATE (falsifiers F3, G9d)
        BLOCK
          REAL(wp) :: axes_std, psum4
          INTEGER  :: j_c
          axes_std = 0.0_wp ; psum4 = 0.0_wp
          DO j_c = 0, CT_MAXTIL
            axes_std = axes_std + gnx_dens(j_c)
            psum4    = psum4    + gnx_pool(j_c)
          END DO
          WRITE (unit_gnx,'("W_pot_axis_g,",G14.6,",g per axis = bgn_r_spike x window_degCd_above_4.5 (Fischer''s unshaded' // &
            ' spike over the REALISED window; DECLARED role, rec 98 [B-GN3])")') gnx_wpot_axis
          WRITE (unit_gnx,'("axes_standing_m2,",G14.6,",axes m-2 standing at anthesis = SUM_c density_m2 (rec 99 deaths' // &
            ' already removed; d2_cohort_table.csv)")') axes_std
          WRITE (unit_gnx,'("earless_axes_m2,",G14.6,",axes m-2 standing WITHOUT a spike = SUM_c density_m2 (1 - phi_sp);' // &
            ' Fischer 1985 p.452 class; Jiang 2022 Tab 7 np-T# band 4-32 % (harvest, pot) REPORTED")') gnx_earless_m2
          WRITE (unit_gnx,'("standing_spikeless_fraction,",G14.6,",earless_axes_m2 / axes_standing_m2")') &
            MERGE(gnx_earless_m2/axes_std, 0.0_wp, axes_std > 0.0_wp)
          WRITE (unit_gnx,'("fate_identity_resid,",ES14.6,",axes_standing_m2 - ears_m2 - earless_axes_m2, from STATE;' // &
            ' must be 0 to 1e-9 (F3)")') axes_std - gnx_ears_m2 - gnx_earless_m2
          WRITE (unit_gnx,'("pool_return_resid_g_m2,",ES14.6,",SUM_c spike_DM_pool_g_m2 - spike_DW_at_anthesis_g_m2, from' // &
            ' STATE; carbon-neutral return, must be 0 to 1e-12 (F3)")') psum4 - gnx_spike_dw_anth
        END BLOCK
      END IF
      ! ---- Record 134: the tier-3 nitrogen cap, the window-end value and the
      ! test-day recompute side by side, because any floret or spike nitrogen cap
      ! used at anthesis must be computed AT the anthesis or spike-test day from the
      ! CURRENT standing nitrogen status and never inherited from the tiller window.
      ! Written only when the cap could apply (record 90 = 4) AND record 134 is on.
      ! ---- Record 137: the floret-survival nitrogen term at the test day, written
      ! only when the record is on. The carbon-only comparator - the same survival
      ! with the nitrogen factor forced to 1 - is printed beside the nitrogen-scaled
      ! survival actually used, per cohort, so that the effect of the nitrogen term
      ! is readable without a second run.
      WRITE (unit_gnx,'("env_MFS_range,",A)') '9 - 12 primordia per spikelet (Guo & Schnurbusch 2015) REPORTED'
      ! The like-for-like pair: the published floret survival is a per-SPIKELET
      ! fraction on MAIN SHOOTS, so the model quantity that matches it is the
      ! MAIN-STEM class and not the crop mean.
      WRITE (unit_gnx,'("floret_survival_MAIN_STEM,",G14.6,",REPORTED NOT GATED (the source gives no' // &
        ' spike DW, so no equal-mass comparison is possible) vs Guo & Schnurbusch 2015' // &
        ' control (tillered) MAIN-SHOOT means 0.337-0.441: greenhouse 3.11/9.23, 3.94/10.20, 3.78/10.09;' // &
        ' field 3.56/9.72, 4.69/10.64, 4.64/10.67")') gnx_surv(0)
      WRITE (unit_gnx,'("fertile_florets_per_spikelet_MAIN_STEM,",G14.6,",REPORTED NOT GATED vs Guo &' // &
        ' Schnurbusch 2015 control MAIN-SHOOT measured 3.11-4.69 florets per spikelet")') &
        MERGE(gnx_nfl(0)/gnx_spikelets, 0.0_wp, gnx_spikelets > 0.0_wp)
      WRITE (unit_gnx,'("floret_survival_CROP_MEAN,",G14.6,",THE GATE (like-for-like: a CROP quantity' // &
        ' against a FIELD-CROP measurement) = Fischer-derived 0.216-0.232; NOT the Guo main-shoot value")') &
        MERGE(gnx_florets_m2/(gnx_ears_m2*gnx_ninit), 0.0_wp, gnx_ears_m2*gnx_ninit > 0.0_wp)
      WRITE (unit_gnx,'("fertile_florets_per_spikelet_CROP,",G14.6,",THE GATE, LIKE-FOR-LIKE vs Fischer 1985 p.450' // &
        ' FIELD CROP 1.98 kernels/spikelet / 0.853-0.916 kernels per competent floret (Tab 4 p.451)' // &
        ' = 2.16-2.32 fertile florets per spikelet")') &
        MERGE(gnx_florets_m2/(gnx_ears_m2*gnx_spikelets), 0.0_wp, gnx_ears_m2*gnx_spikelets > 0.0_wp)
      WRITE (unit_gnx,'("env_floret_survival_range,",A)') 'CORRECTED 2026-09-17, bound AND basis.' // &
        ' WAS "0.40 - 0.60 across 30 modern cultivars (Guo, Slafer & Schnurbusch 2016 JXB 67:4241)",' // &
        ' printed beside the model CROP MEAN. The source sentence is "the relative variation in floret' // &
        ' survival (~50% difference FROM <0.4 TO ~0.55-0.6)" -- "<0.4" is not a measured floor of 0.40 --' // &
        ' and it defines survival as "the proportion of the maximum number of floret primordia that' // &
        ' reached the stage of fertile florets (FFS MFS-1)", measured on MAIN SHOOTS (Guo & Schnurbusch' // &
        ' 2015 Materials and methods: "main shoots of three plants for each cultivar were randomly' // &
        ' selected to measure floret primordia and fertile floret number per spikelet") at apical,' // &
        ' central and basal spikelet positions. IS NOW: MAIN-SHOOT survival 0.337-0.441, the six control' // &
        ' (normally tillered) means of Guo & Schnurbusch 2015 -- NARROWER at the top than the old gate,' // &
        ' not widened to admit the model -- with the 30-cultivar genotypic spread "<0.4 to ~0.55-0.6"' // &
        ' reported alongside but not used as the gate. WHICH COMPARISON IS THE GATE: the model grain' // &
        ' number is a CROP quantity and Fischer 1985 p.450 measured a FIELD CROP (1.98 kernels/spikelet,' // &
        ' 436 spikes/m2), while Guo measured MAIN SHOOTS of selected plants and reports NO spike dry' // &
        ' weight, so Guo cannot be compared at equal spike mass - and survival here is a function of' // &
        ' spike assimilate per axis. Under the like-for-like rule the GATE is the CROP comparator' // &
        ' (fertile_florets_per_spikelet_CROP vs 2.16-2.32) and the main-shoot value is REPORTED beside' // &
        ' it, not gated. Guo et al. 2016 also states at first hand that "the rate of floret survival is' // &
        ' generally low" and that "increased plant density reduces both the number of floret primordia' // &
        ' initiated and floret survival", which is the direction of the remaining main-shoot gap.'
      WRITE (unit_gnx,'("env_kernels_per_fertile_floret,",A)') '0.853 - 0.916 (Fischer 1985 Table 4 p.451,' // &
        ' kernels at maturity / competent florets at anthesis) REPORTED; the model applies 0.78 (p.447/p.449)'
      WRITE (unit_gnx,'("NOT_VERIFIED_Fischer1983,",A)') 'the "1 floret per 10 mg" SLOPE at Fischer 1985' // &
        ' p.449 is a SECONDARY citation of Fischer (1983), which is not on this machine. The VALUE is' // &
        ' rested on first-hand measurement instead: Fischer Table 4 p.451 102.3/103.3 florets g-1 and' // &
        ' Prieto et al. 2020 JXB 71:2703 100.9-123.8 florets g-1.'
      WRITE (unit_gnx,'("NOT_MODELLED,",A)') 'genotypic variation in spikelet number; post-anthesis grain' // &
        ' abortion as a separate process (folded into kernels per fertile floret); Abichou Eqn T4.4 ramp' // &
        ' (its Ne is an experimental INPUT, and the shipped die fraction 0.632 rests on Evers 2006, whose' // &
        ' file here is a 1816-byte anti-bot HTML page = NOT VERIFIED).'
    ELSE
      ! Record 90 = 0 reproduces the record-89 file including these three
      ! statements. They are true of the record-90 = 0 path and false of the floret
      ! chain, which is why they are branched.
      WRITE (unit_gnx,'("spikelets_per_spike,",A)') 'NOT DETERMINED - not modelled; no primary' // &
        ' source on this machine (Kirby 1985 = anti-bot HTML page; Evers 2006 and Moreau 2012 =' // &
        ' 1816-byte stubs; Ferrante / Gonzalez / Serrago / Miralles absent)'
      WRITE (unit_gnx,'("florets_INITIATED_per_spikelet,",A)') 'NOT DETERMINED - not modelled; same sweep'
      WRITE (unit_gnx,'("tiller_mortality_law,",A)') 'NOT DETERMINED - no axis dies under rec 82;' // &
        ' ears m-2 above is the axis population that carries >= 1 kernel'
    END IF
    WRITE (unit_gnx,'(A)') '# --- Fischer 1985 measured components: REPORTED ENVELOPE, read by nothing ---'
    WRITE (unit_gnx,'("env_Fischer_T3_spikes_m2,",G14.6,",Table 3 control p.450, Yecora 70 CIANO")') gnx_env_spikes_t3
    WRITE (unit_gnx,'("env_Fischer_T3_spikelets_per_spike,",G14.6,",Tab 3 p.450")') gnx_env_spklt_t3
    WRITE (unit_gnx,'("env_Fischer_T3_kernels_per_spikelet,",G14.6,",Table 3 control p.450")') gnx_env_kspk_t3
    WRITE (unit_gnx,'("env_Fischer_T4_spikes_m2,",G14.6,",Table 4 p.451 at anthesis, cv. Egret")') gnx_env_spikes_t4
    WRITE (unit_gnx,'("env_Fischer_T4_florets_per_spike,",G14.6,",Tab 4 p.451 anth, Egret")') gnx_env_flspk_t4
    WRITE (unit_gnx,'("env_Fischer_T4_kernels_per_spike,",G14.6,",Table 4 p.451 at maturity, cv. Egret")') gnx_env_kspk_t4
    WRITE (unit_gnx,'("env_Fischer_spike_growth_period_degCd,",G14.6,",p.449 Yecora 70")') gnx_tt_ref_Cd
    WRITE (unit_gnx,'(A)') '# --- inputs this record RETIRES (read by nothing on this path) ---'
    IF (grain_number_mech == 1) THEN
      IF (fe_trait_mode == 1) THEN
        WRITE (unit_gnx,'("fruiting_efficiency_trait,",G14.6,",READ under rec 151 = 1: m_floret,cv = 0.78/FE; scale ",G14.6)') &
          fruiting_efficiency, gnx_floret_scale
      END IF
      WRITE (unit_gnx,'("retired_storage_onset_dvs2,",G14.6,",NOT READ under rec 89 = 1")') storage_onset_dvs2
      WRITE (unit_gnx,'("retired_storage_peak_frac,",G14.6,",NOT READ under rec 89 = 1")') storage_peak_frac
    END IF
    ! Record 138: on-only summary rows, written only when the record is at or above
    ! 1.
    ! Record 139: on-only summary rows, written only when the record is at or above
    ! 1.
    IF (n_uptake_diffusion >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 139 n_uptake_diffusion: diffusion to the root surface, N uptake supply law ---'
      WRITE (unit_gnx,'("rec139_n_uptake_diffusion,",I0,",0 OFF / 1 ON")') n_uptake_diffusion
      WRITE (unit_gnx,'("rec139_S_diff_season_gN_m2,",G14.6,",season cumulative REALIZED diffusion-term N uptake' // &
        ' (post-cap, post-ratio)")') nud_s_diff_cum
      WRITE (unit_gnx,'("rec139_S_mf_season_gN_m2,",G14.6,",season cumulative REALIZED mass-flow-term N uptake' // &
        ' (post-cap, post-ratio)")') nud_s_mf_cum
      WRITE (unit_gnx,'("rec139_S_act_season_gN_m2,",G14.6,",season cumulative REALIZED active-transporter-term' // &
        ' N uptake (post-cap, post-ratio)")') nud_s_act_cum
      WRITE (unit_gnx,'("rec139_days_cap_bound,",I0,",days ANY layer''s standing pool (not the rate terms' // &
        ' mf+act+diff) was the binding cap")') nud_days_cap_bound
    END IF
    ! Record 139 mode 2: an on-only summary row, modes 2 and above only.
    ! Record 139 mode 3: on-only summary rows, modes 3 and above only.
    ! Record 140: on-only summary rows, written only when the record is at or above
    ! 1.
    IF (mineral_supply_mech >= 1) THEN
      BLOCK
        INTEGER :: i_gmsm
        REAL(wp) :: gmsm_conc
        WRITE (unit_gnx,'(A)') '# --- rec 140 mineral_supply_mech: diffusion/HATS mineral supply mechanics ---'
        WRITE (unit_gnx,'("rec140_mineral_supply_mech,",I0,",0 OFF / 1 ON")') mineral_supply_mech
        DO i_gmsm = 1, NMIN
          gmsm_conc = 0.0_wp
          IF (W_grain > 1.0e-9_wp) gmsm_conc = 1000.0_wp * M_grain(i_gmsm) / W_grain
          WRITE (unit_gnx,'("rec140_",A,"_S_mf_season_mg_m2,",G14.6,",season cumulative REALIZED mass-flow uptake")') &
            TRIM(mineral_names(i_gmsm)), msm_s_mf_cum(i_gmsm)
          WRITE (unit_gnx,'("rec140_",A,"_S_diff_season_mg_m2,",G14.6,",season cumulative REALIZED diffusion' // &
            '-term uptake")') TRIM(mineral_names(i_gmsm)), msm_s_diff_cum(i_gmsm)
          WRITE (unit_gnx,'("rec140_",A,"_S_hats_season_mg_m2,",G14.6,",season cumulative REALIZED Zn-HATS' // &
            '-term uptake (0 for every element but Zn)")') TRIM(mineral_names(i_gmsm)), msm_s_hats_cum(i_gmsm)
          WRITE (unit_gnx,'("rec140_",A,"_grain_conc_final,",G14.6,",mg kg-1, final grain concentration")') &
            TRIM(mineral_names(i_gmsm)), gmsm_conc
        END DO
      END BLOCK
    END IF
    ! rec 141: ON-only summary rows, ONLY when lamina_n_remob_pre >= 1
    ! (OFF file above is: nothing in this IF executes).
    IF (lamina_n_remob_pre >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 141 lamina_n_remob_pre: pre-anthesis lamina N remobilization to stem/chaff ---'
      WRITE (unit_gnx,'("rec141_lamina_n_remob_pre,",I0,",0 OFF / 1 ON")') lamina_n_remob_pre
      WRITE (unit_gnx,'("rec141_rem_used_cum_gN_m2,",G14.6,",season cumulative lamina N remobilized to stem/chaff")') &
        lnr_rem_used_cum
      WRITE (unit_gnx,'("rec141_days_rem_used,",I0,",-")') lnr_days_rem_used
      WRITE (unit_gnx,'("rec141_rem_used_window_gN_m2,",G14.6,",SUM rem_used over the Fischer' // &
        ' pre-anthesis window (gnx_win_open .AND. dvs2 < 1)")') lnr_win_rem_used_sum
    END IF
    ! Record 142: on-only summary rows, written only when the record is at or above
    ! 1.
    ! Record 143: on-only summary rows, written only when the record is at or above
    ! 1.
    IF (n_status_construction_price >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 143 n_status_construction_price: whole-crop N-status construction price ---'
      WRITE (unit_gnx,'("rec143_n_status_construction_price,",I0,",0 OFF / 1 ON")') n_status_construction_price
      WRITE (unit_gnx,'("rec143_cut_L_cum_gDM_m2,",G14.6,",season cumulative lamina DM refused relative to a' // &
        ' floor-only counterfactual, DIAGNOSTIC-ONLY")') nsc_cut_L_cum
      WRITE (unit_gnx,'("rec143_cut_S_cum_gDM_m2,",G14.6,",season cumulative stem+chaff DM refused relative to a' // &
        ' floor-only counterfactual, DIAGNOSTIC-ONLY")') nsc_cut_S_cum
      WRITE (unit_gnx,'("rec143_cut_R_cum_gDM_m2,",G14.6,",season cumulative root DM refused relative to a' // &
        ' floor-only counterfactual, DIAGNOSTIC-ONLY")') nsc_cut_R_cum
      IF (nsc_nni_veg_days > 0) THEN
        WRITE (unit_gnx,'("rec143_nni_veg_mean_veg,",G14.6,",mean nap_nni_veg over this record''s own live' // &
          ' (dvs2 < 1) days")') nsc_nni_veg_sum / REAL(nsc_nni_veg_days, wp)
      ELSE
        WRITE (unit_gnx,'("rec143_nni_veg_mean_veg,",A,",no live days counted")') 'NA'
      END IF
    END IF
    ! Record 130 mode 3: on-only summary rows, written only at mode 3, the soft
    ! structural floor.
    IF (grain_n_min_conc == 3) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 130 grain_n_min_conc = 3: Martre (2003) soft structural grain-N floor ---'
      WRITE (unit_gnx,'("rec130_mode,",I0,",0 OFF / 1 cumulative hard / 2 guarded daily hard / 3 Martre soft structural")') &
        grain_n_min_conc
      WRITE (unit_gnx,'("rec130_c_stru_final_pct,",G14.6,",final-day gnm_c_stru_day x 100, Option C reference %DM")') &
        100.0_wp * gnm_c_stru_day
      WRITE (unit_gnx,'("rec130_days_floor_bound,",I0,",days the mode-3 ceiling bound room (reuses the shared' // &
        ' gnm_days_bind counter)")') gnm_days_bind
      WRITE (unit_gnx,'("rec130_Nstru_final_gN_m2,",G14.6,",final-day kernel_number x Nstru(Tt), Eq. 3 cumulative demand")') &
        gnm_Nstru_total
    END IF
    ! Record 144: five on-only summary rows, written only when the record is at or
    ! above 1.
    IF (carbon_coupled_leaf_expansion >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 144 carbon_coupled_leaf_expansion: carbon-coupled leaf-area expansion ---'
      WRITE (unit_gnx,'("rec144_carbon_coupled_leaf_expansion,",I0,",0 OFF / 1 ON")') carbon_coupled_leaf_expansion
      WRITE (unit_gnx,'("rec144_days_surplus,",I0,",days cbc_surplus_day > 0 AND lew_fw >= 1 (the accumulator' // &
        ' was actually incremented)")') cbc_days_surplus
      WRITE (unit_gnx,'("rec144_f_max_reached,",G14.6,",running MAX of cbc_f(c,x) over every axis/rank cell;' // &
        ' ceiling cbc_fmax = 1.20 (Masle 2000 Table II)")') cbc_f_max_reached
      WRITE (unit_gnx,'("rec144_area_bonus_cum_m2_m2,",G14.6,",running MAX of cbc_area_bonus_day over the' // &
        ' season -- the peak standing bonus ground-area reached (a LEVEL, not a flux-sum)")') cbc_area_bonus_cum
      WRITE (unit_gnx,'("rec144_LAImax,",G14.6,",running MAX of cls_gai (the model''s own green area index)' // &
        ' over the season")') cbc_laimax
    END IF
    ! Records 145 and 146: a per-run summary line each, on only. Both are static
    ! configuration switches with no per-day accumulator, so one row each records
    ! the setting actually used this run.
    IF (leaf_maintenance_basis >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 145 leaf_maintenance_basis: leaf term removed from maint_raw, Rd only ---'
      WRITE (unit_gnx,'("rec145_leaf_maintenance_basis,",I0,",0 OFF (Teh leaf term + Rd) / 1 ON' // &
        ' (leaf term removed, Rd only)")') leaf_maintenance_basis
    END IF
    IF (par_reflection_basis >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 146 par_reflection_basis: single Teh per-waveband PAR reflection coefficient ---'
      WRITE (unit_gnx,'("rec146_par_reflection_basis,",I0,",0 OFF (double reflection, const_P_sr then' // &
        ' const_P_par) / 1 ON (const_P_par only, solar.f90:88,92)")') par_reflection_basis
    END IF
    ! Record 147: on-only summary rows.
    ! Record 179: on-only summary rows.
    IF (lamina_expansion_n_limit >= 1) THEN
      WRITE (unit_gnx,'(A)') '# --- rec 179 lamina_expansion_n_limit: N limit on the lamina-area EXPANSION increment ---'
      WRITE (unit_gnx,'("rec179_lamina_expansion_n_limit,",I0,",0 OFF / 1 PRIMARY (Bertheloot 2008' // &
        ' FPB Eqn 8 minN^ph,a, a 0.46 b 1.9 c 0.24, Table 1 p.784) / 2 arm (Eqn 8 plateau 1.114)' // &
        ' / 3 arm (Jamieson & Semenov 1.5 on total lamina N) / 4 attribution arm (Yin l_N log' // &
        ' form, kN 0.41, nb 0.40)")') lamina_expansion_n_limit
      WRITE (unit_gnx,'("rec179_minNpha_end_gN_m2,",G14.6,",the Eqn-8 minimum photosynthetic N per' // &
        ' unit lamina area in force on the last evaluated day")') lnl_minNpha_day
      WRITE (unit_gnx,'("rec179_maxNpha_gN_m2,",G14.6,",the running maximum photosynthetic SLN the' // &
        ' lamina reached (Eqn 8 maxN^ph,a); modes 2-4 do not use it")') lnl_maxNpha
      WRITE (unit_gnx,'("rec179_refused_area_cum_m2_m2,",G14.6,",season sum of the lamina area THIS' // &
        ' record refused BEYOND what records 121/122 already refused")') lnl_refused_cum
      WRITE (unit_gnx,'("rec179_days_refuse,",I0,",days this record refused area beyond rec 121/122")') lnl_days_refuse
      ! rec 179: the green-area guard on Eqn 8's maxN^ph,a.
      WRITE (unit_gnx,'("rec179_gai_max_min_m2_m2,",G14.6,",GREEN-AREA GUARD on the running maximum' // &
        ' maxN^ph,a: a DECLARED DOMAIN guard (SHARED.f90 [LNL-G-1]), not a literature constant;' // &
        ' below it N^ph/cls_gai is not the SLN of any lamina (Bertheloot 2008 FPB p.784, maxN^ph,a' // &
        ' is reached BY THE LAMINA); mode 1 only")') lnl_gai_max_min
      WRITE (unit_gnx,'("rec179_days_maxNpha_updated,",I0,",days the guard admitted a new running' // &
        ' maximum (mode 1 only; 0 at modes 2-4, which never read maxN^ph,a)")') lnl_days_maxupd
      WRITE (unit_gnx,'("rec179_gai_at_last_maxupd_m2_m2,",G14.6,",cls_gai on the last day the' // &
        ' running maximum was updated")') lnl_gai_at_maxupd
    END IF
    CLOSE(unit_gnx)
  END BLOCK
END IF

! ==================================================================== !
! Record 91: the grain-filling law at the end of the run, written only when the
! record is on. Every row carries its source or says that the quantity is not
! modeled; nothing here is fitted to any observation.
! ==================================================================== !
IF (grain_fill_mech >= 1) THEN
  BLOCK
    INTEGER :: unit_gfb
    REAL(wp) :: tmean_fill, tkw_real, rate_mean
    tmean_fill = 0.0_wp
    IF (gfb_days > 0) tmean_fill = gfb_tsum / REAL(gfb_days, wp)
    tkw_real = 0.0_wp
    IF (kernel_number_m2 > 1.0e-9_wp) tkw_real = 1000.0_wp * W_grain / kernel_number_m2
    rate_mean = 0.0_wp
    IF (gfb_days > 0 .AND. kernel_number_m2 > 1.0e-9_wp) &
      rate_mean = 1000.0_wp * gfb_supplied_tot / kernel_number_m2 / REAL(gfb_days, wp)
    OPEN (newunit=unit_gfb,file='output/grain_filling_law.csv',status='replace',action='write')
    WRITE (unit_gfb,'(A)') 'quantity,value,units_or_source'
    WRITE (unit_gfb,'("rec91_grain_fill_mech,",I0,",0 OFF / 1 rate x duration / 2 shape-off arm / 3 cap-kept arm")') &
      grain_fill_mech
    WRITE (unit_gfb,'(A)') '# --- the constants, all from ONE table (SHARED.f90) ---'
    WRITE (unit_gfb,'("c_m_at_T_ref,",G14.6,",mg kernel-1 d-1; mean c_m of Yin et al. 2003 Ann Bot 91:361' // &
      ' Table 2 regime T1 (20/15 C), six CIMMYT genotypes; DOI 10.1093/aob/mcg029, PMC4244967")') gfb_cm_ref
    WRITE (unit_gfb,'("c_m_slope,",G14.6,",mg kernel-1 d-1 K-1 = (2.481667 - 1.928333)/5 from the same table")') &
      gfb_cm_slope
    WRITE (unit_gfb,'("T_ref,",G14.6,",degC; daily mean of Yin''s 20/15 regime at a 12 h photoperiod.' // &
      ' DECLARED ASSUMPTION: the paper states no photoperiod. D_fill is INVARIANT to it (405.581 for' // &
      ' 12/14/16 h); only T_base moves, 6.6187 -> 7.4521")') gfb_t_ref
    WRITE (unit_gfb,'("T_base,",G14.6,",degC; DERIVED, not asserted: reciprocal-of-duration method of' // &
      ' Fischer 1985 Table 6 footnote p.455 (printed; pdf p.9) applied to Yin''s two mean t_e (37.2733 / 25.5383 d)")') gfb_tbase
    WRITE (unit_gfb,'("D_fill,",G14.6,",degCd above T_base; same derivation; identical from BOTH Yin regimes")') &
      gfb_dfill
    WRITE (unit_gfb,'("t_m_over_t_e,",G14.6,",beta shape; mean of ALL TWELVE rows of Yin Table 2")') gfb_u
    WRITE (unit_gfb,'("beta_shape_integral,",G14.6,",mean rate / c_m at that shape (analytic 0.685711)")') &
      gfb_shape_int
    WRITE (unit_gfb,'(A)') '# --- what this deck actually did ---'
    WRITE (unit_gfb,'("filling_start_model_day,",I0,",day the grain number was set (= anthesis)")') gfb_day0
    WRITE (unit_gfb,'("filling_end_model_day,",I0,",0 = still filling at the last simulated day")') gfb_day_end
    WRITE (unit_gfb,'("filling_days_realised,",I0,",d; THE DURATION -- an OUTCOME, not an input")') gfb_days
    WRITE (unit_gfb,'("filling_degCd_realised,",G14.6,",degCd above T_base accumulated")') gfb_tt
    BLOCK
      CHARACTER(LEN=17) :: eb_txt
      SELECT CASE (gfb_end_bound)
      CASE (1); eb_txt = 'THERMAL_DURATION '
      CASE (2); eb_txt = 'PHENOLOGY_dvs2_2 '
      CASE (3); eb_txt = 'KERNEL_MASS_CAP  '
      CASE (4); eb_txt = 'KERNEL_POTENTIAL '
      CASE DEFAULT; eb_txt = 'STILL_FILLING    '
      END SELECT
      IF (kernel_mass_ceiling == 1) THEN
        WRITE (unit_gfb,'("rec118_kernel_mass_ceiling,1,ON: per-kernel mass bounded by TKW_potential = ",G14.6,' // &
          '" g per 1000 kernels (SHARED.f90); a third stop KERNEL_MASS_CAP can fire")') TKW_potential
      END IF
      IF (kernel_potential_law >= 2) THEN
        WRITE (unit_gfb,'("rec125_kernel_potential_mg,",G14.6,",ON: W_pot = MIN(env_max MAX(env_min s_k_mean x D_fill))' // &
          ' set on day ",I0," (SHARED.f90); bound 0 interior / 1 env_min / 2 env_max = ",I0)') &
          kpt_w_pot, kpt_day_set, kpt_bound
        WRITE (unit_gfb,'("rec125_window_supply_per_kernel,",G14.6,",mg kernel-1 d-1 mean over the 250 degCd window;' // &
          ' W_raw = ",G14.6," mg; Yin potential at the window mean T = ",G14.6," mg; envelope ",G14.6," - ",G14.6,' // &
          '" mg (cultivar_yield.csv rows 7-8)")') kpt_sk_mean, kpt_w_raw, kpt_w_yin, kpt_env_min, kpt_env_max
      END IF
      IF (kernel_potential_law == 1 .OR. kernel_potential_law == 2) THEN
        ! This line prints THE VALUE the per-kernel grain-nitrogen ceiling actually
        ! used, which is the kernel number times the live per-kernel potential, and
        ! names where that value came from: the compiled 1.38 mg N kernel-1 on the
        ! record-177 = 0 path, or, under record 177 = 1, the cultivar
        ! grainN_pot_mgN_kernel row of input/cultivar_yield.csv, or the control maximum
        ! 1.11 mg N kernel-1 of Martre, Porter, Jamieson & Triboi (2003) Plant Physiol
        ! 133:1959 Table I p.1961 when the cultivar has no row.
        WRITE (unit_gfb,'("rec125_grainN_ceiling_gN_m2,",G14.6,",kernel_number x ",F6.4," mg N kernel-1' // &
          ' IN FORCE (rec 177 = 1: the cultivar_yield.csv grainN_pot_mgN_kernel row if present, else the' // &
          ' Martre et al. 2003 Plant Physiol 133:1959 Table I p.1961 control maximum 1.11; rec 177 = 0:' // &
          ' the compiled gnf_mgN_kernel_hi = 1.38 from the ABLATED-ear treatment of the same table);' // &
          ' cut cumulative ",G14.6," g N m-2 on ",I0," days (SHARED.f90)")') &
          gnc_cap, gnp_mgN_kernel, gnc_cut_cum, gnc_days_bound
      END IF
      ! Record 177, on only: say in the run.log WHETHER the per-kernel
      ! grain-nitrogen ceiling BOUND this deck, and with which value and provenance.
      ! It can bind exactly - on four decks the grain nitrogen concentration of
      ! 3.2093 % is 1.38 mg N kernel-1 over a 43 mg kernel - and nothing else in the
      ! run.log says so.
      IF (kernel_potential_law >= 2) THEN
        WRITE (unit_gfb,'("filling_end_bound,",A)') eb_txt // ',which of the four stops fired -- REPORTED, never assumed'
      ELSE IF (kernel_mass_ceiling == 1) THEN
        WRITE (unit_gfb,'("filling_end_bound,",A)') eb_txt // ',which of the three stops fired -- REPORTED, never assumed'
      ELSE
        WRITE (unit_gfb,'("filling_end_bound,",A)') eb_txt // ',which of the two stops fired -- REPORTED, never assumed'
      END IF
    END BLOCK
    WRITE (unit_gfb,'("T_mean_filling_C,",G14.6,",degC, screen air, mean over the filling days.' // &
      ' Yin''s measured bracket is 17.5-22.5: values inside it are INTERPOLATION")') tmean_fill
    WRITE (unit_gfb,'("c_m_min_over_window,",G14.6,",mg kernel-1 d-1")') MERGE(gfb_cm_min, 0.0_wp, gfb_days > 0)
    WRITE (unit_gfb,'("c_m_max_over_window,",G14.6,",mg kernel-1 d-1")') gfb_cm_max
    WRITE (unit_gfb,'("c_m_zero_clamped,",A,",TRUE would mean a filling day below 0.08 degC")') &
      MERGE('YES', 'NO ', gfb_cm_clamped)
    WRITE (unit_gfb,'("realised_rate_max,",G14.6,",mg kernel-1 d-1, largest POTENTIAL daily rate reached")') &
      gfb_rate_max
    WRITE (unit_gfb,'("realised_rate_mean,",G14.6,",mg kernel-1 d-1, carbon ACTUALLY delivered / kernels / days.' // &
      ' FALSIFIER vs Yang & Zhang 2006 New Phytol 169:223 Table 3 p.228 wheat MEASURED 0.95-1.94")') rate_mean
    IF (grain_fill_mech == 3) THEN
      WRITE (unit_gfb,'("sink_demand_total_g_m2,",A)') 'NOT DEFINED on this arm: rec 91 = 3 keeps the' // &
        ' legacy STOCK capacity cap, so a cumulative FLUX demand has no meaning here'
      WRITE (unit_gfb,'("unmet_demand_total_g_m2,",A)') 'NOT DEFINED on this arm, same reason'
    ELSE
      WRITE (unit_gfb,'("sink_demand_total_g_m2,",G14.6,",g m-2, integral of the potential demand")') gfb_demand_tot
      WRITE (unit_gfb,'("unmet_demand_total_g_m2,",G14.6,",g m-2 CUMULATIVE. This REPLACES the OFF path''s' // &
        ' stock ''sink room left'': under a flux law that stock has no meaning. > 0 = SOURCE-limited")') &
        gfb_unmet_tot
    END IF
    WRITE (unit_gfb,'("supplied_to_grain_g_m2,",G14.6,",g m-2, assimilate + reserve actually delivered")') &
      gfb_supplied_tot
    WRITE (unit_gfb,'("reserve_C_pool_end_g_m2,",G14.6,",g m-2 W_reserve after supply on the LAST day the law' // &
      ' was live (report fix 2026-10-15: the harvest-day value printed before this was zeroed by the rec 105' // &
      ' release). FALSIFIER: it was 58.5-451.8 across the 20 ceiling-pinned decks under the capacity cap, and' // &
      ' must FALL where the real sink is larger")') gfb_reserve_end
    WRITE (unit_gfb,'("N_floor_refused_dm_total_g_m2,",G14.6,",g m-2 CUMULATIVE grain DM refused by the rec 130' // &
      ' structural-N floor AFTER supply (same quantity as run.log refused grain DM)")') gnm_refused_dm_cum
    WRITE (unit_gfb,'("demand_cut_before_supply_g_m2,",G14.6,",g m-2 = sink_demand_total - supplied_to_grain' // &
      ' - unmet_demand_total: the part of the potential demand removed by the N floor and the other caps' // &
      ' BEFORE the carbon supply saw it (report fix 2026-10-15; 0 where no cap bound)")') &
      gfb_demand_tot - gfb_supplied_tot - gfb_unmet_tot
    WRITE (unit_gfb,'("grain_number_m2,",G14.6,",from the grain-number chain (rec 89/90);' // &
      ' the grain-filling law does not touch it")') &
      kernel_number_m2
    WRITE (unit_gfb,'("TKW_g_realised,",G14.6,",g per 1000 kernels; THE OUTCOME = integral of the realised' // &
      ' rate / kernel number. NOT bounded by any deck constant")') tkw_real
    WRITE (unit_gfb,'("TKW_potential_input_g,",G14.6,",the deck constant. NOT READ under rec 91 = 1 or 2;' // &
      ' still read under rec 91 = 3 (the attribution arm that keeps the cap)")') TKW_potential
    WRITE (unit_gfb,'(A)') '# --- measured envelopes: falsification targets, never tuning targets ---'
    WRITE (unit_gfb,'("env_filling_days,",A)') '21.60 - 42.95 d (t_e, Yin 2003 Table 2, twelve rows);' // &
      ' 19 - 34 d (active grain-filling period, Yang & Zhang 2006 Table 3 p.228, wheat cv. Yangmai 158)'
    WRITE (unit_gfb,'("env_mean_rate,",A)') '0.95 - 1.94 mg kernel-1 d-1 (Yang & Zhang 2006 Table 3 p.228,' // &
      ' wheat, four N x soil-moisture treatments); the law gives 1.32 at 17.5 C and 1.70 at 22.5 C'
    WRITE (unit_gfb,'("env_kernel_weight,",A)') '36 - 45 mg (Yang & Zhang 2006 Table 3 p.228, wheat);' // &
      ' 39.73 - 52.20 mg (w_max, Yin 2003 Table 2)'
    WRITE (unit_gfb,'("env_duration_base0_degCd,",A)') '652.3 / 574.6 (Yin Table 2 means on a base-0 basis);' // &
      ' 498 (Apache) and 567 (Renan), Baillot et al. 2018 PLoS ONE 13:e0209597 Results, base 0;' // &
      ' 614 - 711, Girousse et al. 2023 PLoS ONE 18:e0285218 (RANGE ONLY, transcription discrepancy flagged)'
    WRITE (unit_gfb,'("env_duration_response_bracket,",A)') 'THE OPEN TENSION. Shortening factor per +5 degC:' // &
      ' 1.24 Wardlaw 2002 Ann Bot 90:469 (18/13 -> 27/22 C, -32% over 9 C); 1.46 CODED (Yin 2003 Table 2);' // &
      ' 1.73 Wang et al. 2017 Nature Plants 3:17102 Fig. 4b post-anthesis cardinals (0, 33, 44 C), derived' // &
      ' from Hot Serial Cereal at MARICOPA cv. Yecora Rojo + INRA over 12.6-24.7 C. Three MEASURED' // &
      ' estimates that DISAGREE; the coded value sits between the brackets, not at an extreme.'
    WRITE (unit_gfb,'("NOT_MODELLED,",A)') 'a cultivar kernel-size trait (so a small-grained cultivar such as' // &
      ' WDF25 cv. Ruta, measured TKW 30.0 g, is OVER-PREDICTED -- the fix by scaling on TKW_potential is' // &
      ' REFUSED because at AR and AU that constant IS the observation and elsewhere it is a 38.0' // &
      ' placeholder); any direct water- or N-stress term on the rate (Fang et al. 2024 JXB' // &
      ' doi 10.1093/jxb/erae310 measured that both suppress it -- here stress reaches the grain ONLY' // &
      ' through the source); any decline of c_m above the optimum (c_m is LINEAR in T and only' // &
      ' interpolated inside Yin''s 17.5-22.5 C bracket).'
    WRITE (unit_gfb,'("NOT_VERIFIED,",A)') 'Sofield et al. 1977 (10.1071/PP9770785), Wardlaw & Moncur 1995' // &
      ' (10.1071/PP9950391), Chowdhury & Wardlaw 1978 (Aust J Agric Res 29:205, 10.1071/AR9780205 -- the' // &
      ' "Aust J Plant Physiol 5:195" form does not exist), Panozzo & Eagles 1999 (10.1071/AR98146) and' // &
      ' Pararajasingham & Hunt 1991 (10.4141/cjps91-007, the one paper that DERIVES a grain-filling base' // &
      ' temperature) are all CLOSED ACCESS with no route from this machine. No number from any of them is used.'
    CLOSE(unit_gfb)
  END BLOCK
END IF

! ==================================================================== !
! Record 92: the stem-reserve law at the end of the run, written only when the
! record is on. Every envelope row is a FALSIFICATION target from a source read
! at first hand, never a tuning target.
! ==================================================================== !
IF (stem_wsc_mech >= 1) THEN
  BLOCK
    INTEGER :: unit_wsl
    REAL(wp) :: f_anth, r_peak, f_after16
    f_anth = -1.0_wp
    IF (wsc_stem_anth > 1.0e-9_wp) f_anth = wsc_pool_anth / wsc_stem_anth
    r_peak = -1.0_wp
    IF (wsc_pool_anth > 1.0e-9_wp) r_peak = wsc_pool_peak / wsc_pool_anth
    f_after16 = -1.0_wp
    IF (wsc_draw_tot > 1.0e-9_wp) f_after16 = wsc_draw_after16 / wsc_draw_tot
    OPEN (newunit=unit_wsl,file='output/stem_wsc_law.csv',status='replace',action='write')
    WRITE (unit_wsl,'(A)') 'quantity,value,units_or_source'
    WRITE (unit_wsl,'("rec92_stem_wsc_mech,",I0,",0 OFF / 1 THE ITEM / 2 cap+close arm / 3 draw arm' // &
      ' / 4 elongation-start sensitivity arm")') stem_wsc_mech
    WRITE (unit_wsl,'(A)') '# --- the constants (SHARED.f90) ---'
    WRITE (unit_wsl,'("wsc_cap_frac,",G14.6,",g WSC per g stem DM; Schnyder 1993 New Phytol 123:233-245' // &
      ' p.238: total WSC concentration in the stem can reach UP TO 50 % of the dry weight of internodes")') &
      wsc_cap_frac
    WRITE (unit_wsl,'("wsc_accum_end_dpa,",I0,",d after anthesis; MIDPOINT of Schnyder 1993 p.238' // &
      ' measured 10-20 d, and p.233 three weeks. Corroborated by Gebbing & Schnyder 1999 Plant Physiol' // &
      ' 121:871 p.875 (mobilisation starts about mid-grain filling)")') wsc_accum_end_dpa
    WRITE (unit_wsl,'("reserve_conv_eff,",G14.6,",g grain DM per g reserve DM. UNCHANGED by item C.' // &
      ' Gebbing & Schnyder 1999 p.877: MEWSC = 0.97 +- 0.07 se; 0.95 is inside one se")') reserve_conv_eff
    WRITE (unit_wsl,'("stem_reserve_frac,",G14.6,",share of the pre-anthesis stem increment LABELLED as' // &
      ' WSC. NOT a measured law and NOT touched by item C: Ruuska et al. 2006 Funct Plant Biol 33:799 is' // &
      ' NOT ON THIS MACHINE, so it is NOT VERIFIED. The OUTCOME is what is tested (next rows)")') &
      stem_reserve_frac
    WRITE (unit_wsl,'("reserve_remob_rate_input,",G14.6,",d-1, the legacy first-order throttle. NOT READ' // &
      ' under rec 92 = 1, 3 or 4; still read under rec 92 = 2 (the attribution arm that keeps it)")') &
      reserve_remob_rate
    WRITE (unit_wsl,'(A)') '# --- what the crop actually did ---'
    WRITE (unit_wsl,'("anthesis_model_day,",I0,",model day the grain number was set (dvs2 = ear_end_dvs2,' // &
      ' = 1.0 anthesis on every shipped deck)")') wsc_day_anth
    WRITE (unit_wsl,'("W_reserve_at_anthesis_g_m2,",G14.6,",g DM m-2")') wsc_pool_anth
    WRITE (unit_wsl,'("W_stem_at_anthesis_g_m2,",G14.6,",g DM m-2 (before that day''s surplus)")') wsc_stem_anth
    WRITE (unit_wsl,'("wsc_frac_of_stem_at_anthesis,",G14.6,",-  THE headline biology number")') f_anth
    WRITE (unit_wsl,'("W_reserve_peak_g_m2,",G14.6,",g DM m-2")') wsc_pool_peak
    WRITE (unit_wsl,'("peak_days_after_anthesis,",I0,",d")') MERGE(wsc_day_peak - wsc_day_anth, -999, &
      wsc_day_anth > 0)
    WRITE (unit_wsl,'("peak_over_anthesis_ratio,",G14.6,",-")') r_peak
    WRITE (unit_wsl,'("W_reserve_end_g_m2,",G14.6,",g DM m-2 stranded at the last day")') W_reserve
    WRITE (unit_wsl,'("deposited_pre_anthesis_g_m2,",G14.6,",g DM m-2 LABELLED by the stem-increment route")') &
      wsc_dep_pre
    WRITE (unit_wsl,'("deposited_post_anthesis_g_m2,",G14.6,",g DM m-2 LABELLED by the sink-full surplus route")') &
      wsc_dep_post
    WRITE (unit_wsl,'("refused_by_capacity_g_m2,",G14.6,",g DM m-2 NOT labelled because the pool was at' // &
      ' 0.50 of stem DM. THE CARBON IS NOT LOST -- it is in W_stem, exactly where it always was")') &
      wsc_refused_cap
    WRITE (unit_wsl,'("refused_by_deposition_close_g_m2,",G14.6,",g DM m-2 NOT labelled because the' // &
      ' deposition window had closed. Same: the carbon is in W_stem")') wsc_refused_win
    WRITE (unit_wsl,'("days_capacity_bound,",I0,",days")') wsc_days_capbnd
    WRITE (unit_wsl,'("days_deposition_closed,",I0,",days")') wsc_days_winbnd
    WRITE (unit_wsl,'("remobilisation_onset_dap_anth,",I0,",d after anthesis of the FIRST draw; -999 = never")') &
      MERGE(wsc_day_mob0 - wsc_day_anth, -999, wsc_day_mob0 > 0 .AND. wsc_day_anth > 0)
    WRITE (unit_wsl,'("drawn_total_g_m2,",G14.6,",g DM m-2 out of the pool")') wsc_draw_tot
    WRITE (unit_wsl,'("delivered_to_grain_g_m2,",G14.6,",g DM m-2 = draw x conversion efficiency")') wsc_gain_tot
    WRITE (unit_wsl,'("conversion_loss_g_m2,",G14.6,",g DM m-2 respired")') wsc_resp_tot
    WRITE (unit_wsl,'("reserve_share_of_grain,",G14.6,",-  delivered_to_grain / final W_grain. FALSIFICATION' // &
      ' BAND: Gebbing & Schnyder 1999 p.877 Conclusions measured pre-anthesis reserves supplying 8-27 % of' // &
      ' total carbohydrate-C in grains")') wsc_gain_tot/MAX(W_grain, 1.0e-9_wp)
    WRITE (unit_wsl,'("fraction_of_draw_after_16dpa,",G14.6,",-  FALSIFICATION: Gebbing & Schnyder 1999' // &
      ' p.875 measured 83 % of stem pre-anthesis C mobilisation occurring after their 16-d sample")') f_after16
    WRITE (unit_wsl,'("carbon_closure_max_g_m2_d,",G14.6,",g DM m-2 d-1, max |dW_grain + dW_stem + dW_chaff' // &
      ' + respired - assimilate_taken| over the run. THIS REPLACES THE RETIRED TAUTOLOGY (item-B audit' // &
      ' note N4): the old cres_ledger_resid was (-draw)+gain+(draw-gain), identically 0 for any inputs")') &
      wsc_closure_max
    WRITE (unit_wsl,'(A)') '# --- measured envelopes: falsification targets, never tuning targets ---'
    WRITE (unit_wsl,'("env_wsc_frac_of_stem_at_anthesis,",A)') '0.125 - 0.250, DERIVED from Schnyder 1993' // &
      ' p.238''s own pair of measurements: the maximum is "up to 50 % of the dry weight of internodes"' // &
      ' and that maximum is "2-4 times greater than at anthesis"'
    WRITE (unit_wsl,'("env_peak_over_anthesis,",A)') '2 - 4 (Schnyder 1993 p.238, verbatim)'
    WRITE (unit_wsl,'("env_peak_days_after_anthesis,",A)') '10 - 20 d (Schnyder 1993 p.238: "WSC accumulate' // &
      ' in the stem until 10-20 d after anthesis"); "within three weeks following anthesis" (Summary p.233)'
    WRITE (unit_wsl,'("env_reserve_share_of_grain,",A)') '0.08 - 0.27 of grain carbohydrate-C' // &
      ' (Gebbing & Schnyder 1999 Plant Physiol 121:871 Conclusions p.877, 13C/12C steady-state labelling,' // &
      ' two spring wheat cultivars x two N levels x two years)'
    WRITE (unit_wsl,'("NOT_MODELLED,",A)') 'the ONSET of deposition (it still starts at emergence, not at' // &
      ' stem elongation: Schnyder p.238 states the onset internode by internode and no source here gives a' // &
      ' CROP-scale node in a unit this model has -- rec 92 = 4 MEASURES how much it matters); any rate' // &
      ' ceiling on the draw (the only rate the sources name is the deficit itself); any WSC in leaf' // &
      ' sheaths, glumes or the peduncle as separate pools.'
    WRITE (unit_wsl,'("DECLARED_TENSIONS,",A)') 'the 15 d close is CALENDAR, as both sources state it, on' // &
      ' temperate crops -- at a hot short-season site it is a larger fraction of the fill; Schnyder''s 50 %' // &
      ' is of INTERNODES while this W_stem includes sheaths, so the cap is GENEROUS; and at NL/AU the cap' // &
      ' inherits item M''s early stem over-build, so their capacity is over-large.'
    ! rec 148 -- ON-ONLY rows (nothing below is written at
    ! wsc_fate_law = 0, so the OFF file is unchanged).
    IF (wsc_fate_law >= 1) THEN
      WRITE (unit_wsl,'(A)') '# --- rec 148 wsc_fate_law: the fate of the stem WSC reserve ---'
      WRITE (unit_wsl,'("rec148_wsc_fate_law,",I0,",0 OFF / 1 the law (term A + term B) / 2 term A' // &
        ' only / 3 term B only")') wsc_fate_law
      WRITE (unit_wsl,'("wfl_term_a_refused_g_m2,",G14.6,",g DM m-2 refused by the RESTORED deposition' // &
        ' close for the grain-sink-full surplus and booked to the existing rec-129 overflow loss' // &
        ' cfo_loss_cum. Schnyder 1993 New Phytol 123:233 p.238 (10-20 d), p.233 Summary (three weeks)")') &
        wfl_a_refused_cum
      WRITE (unit_wsl,'("wfl_term_a_days,",I0,",days term A refused a deposition")') wfl_a_days
      WRITE (unit_wsl,'("wfl_gate_day,",I0,",model day of TERMINAL GRAIN-SINK CLOSURE as term B saw it' // &
        ' (gfb_end_bound /= 0); 0 = never reached")') wfl_gate_day
      WRITE (unit_wsl,'("wfl_gate_bound,",I0,",which bound closed the grain sink: 1 thermal duration /' // &
        ' 2 dvs2 >= 2 / 3 cultivar kernel-mass ceiling / 4 rec 125 kernel potential")') wfl_gate_bound
      WRITE (unit_wsl,'("wfl_pool_at_gate_g_m2,",G14.6,",g DM m-2 W_reserve on the gate day before term' // &
        ' B released anything; -1 = the gate was never reached")') wfl_pool_at_gate
      WRITE (unit_wsl,'("wfl_term_b_released_g_m2,",G14.6,",g DM m-2 released from the pool after the' // &
        ' gate, ALL of it to cres_respired at unity (Schnyder 1993 p.240: whole-plant respiration rises' // &
        ' during the WSC loss with the ear removed). ENVELOPE, read-out only: Spiertz & Ellen 1978 as' // &
        ' quoted on p.240, 220 and 290 g DM m-2 lost between maximum mass and maturity at high and low N")') &
        wfl_b_released_cum
      WRITE (unit_wsl,'("wfl_term_b_days,",I0,",days term B released")') wfl_b_days
      WRITE (unit_wsl,'("wfl_pool_end_g_m2,",G14.6,",g DM m-2 W_reserve at the last day. FALSIFICATION' // &
        ' read-out: the only stem-WSC measurement in the deck set is Braunschweig 11.6-20.3 mg g-1 of' // &
        ' stem DM at maturity (Manderscheid et al. 2020 ODjAR 6:28-33). NO FLOOR IS APPLIED (missing' // &
        ' quantity M-A1): the law empties to ~0 and the undershoot is reported, not calibrated away")') &
        W_reserve
      WRITE (unit_wsl,'("wfl_draw_deficit_g_m2,",G14.6,",g DM m-2 drawn by the PRE-GATE deficit limb,' // &
        ' which this record makes reachable again under rec 105 = 2 (the E1 blackout removed)")') &
        wfl_draw_def_cum
      WRITE (unit_wsl,'("wfl_draw_deficit_days,",I0,",days the pre-gate deficit limb drew")') wfl_draw_def_days
      WRITE (unit_wsl,'("wfl_days_grain_fed,",I0,",FALSIFIER: days term B delivered any carbon to the' // &
        ' grain. MUST be 0 -- room is 0 at and after the gate, so gain2 = 0 by the existing arithmetic")') &
        wfl_days_grain_fed
      WRITE (unit_wsl,'("wfl_grain_fed_cum_g_m2,",G14.6,",g DM m-2 of that. MUST be 0.0")') wfl_grain_fed_cum
    END IF
    CLOSE(unit_wsl)
  END BLOCK
END IF

! Record 105: one run.log summary line, written only when the record is on. The
! late-draw fraction reuses the record-92 accumulator against the 83 % after
! 16 days post-anthesis of Gebbing & Schnyder (1999) p.875, a falsification
! read and never enforced; -1 means no draw.

! Record 105 mode 2: one run.log summary line, written only when the mode is
! on. The two streams are reported SEPARATELY. The stranding statistic is the
! reserve at maturity over the reserve at anthesis plus the summed
! post-anthesis deposition, pre-registered and REPORTED here and enforced by
! nothing in the code. The four conservation residual maxima are printed at
! full precision, because the daily CSV format cannot resolve 1e-12.
IF (wsc_reserve_law == 2) THEN
  BLOCK
    REAL(wp) :: fres2
    fres2 = W_reserve / MAX(1.0e-9_wp, wsc_pool_anth + wsc_dep_post)
    WRITE (*,'(A,I0,6(A,G14.6),A,G14.6,A,I0,4(A,G14.6))') &
      ' rec 105 WSC MODE2 SUMMARY member=', wsc_rel_window_member, &
      ' sum_rel_g_m2=', wsl2_rel_tot, ' sum_rel_grain_g_m2=', wsl2_relg_tot, &
      ' sum_gain_g_m2=', wsl2_gain_tot, ' sum_rel_fate_g_m2=', wsl2_relf_tot, &
      ' sum_loss_conv_g_m2=', wsl2_loss_conv_tot, ' sum_loss_fate_g_m2=', wsl2_loss_fate_tot, &
      ' f_res_maturity=', fres2, ' days_release_zero_deficit=', wsl2_days_zerodef, &
      ' I1_max=', wsl2_i1_max, ' I2_max=', wsl2_i2_max, ' I3_max=', wsl2_i3_max, &
      ' I4_max=', wsl2_i4_max
  END BLOCK
END IF

! Record 117: ONE run.log SUMMARY BLOCK, written only when the sink is on.
! These are the per-deck quantities the sink is judged on: the anthesis chaff
! dry matter; the post-anthesis increment; the fraction of the sourced increment
! REALIZED, which can fall BELOW the sourced total on a supply-limited deck and
! is then reported rather than hidden; the carbon taken from the reproductive
! stream, i.e. from what the grain would otherwise have had, and from the
! reserve, which is identically zero by construction and is printed so that it
! can be checked; and ANY UNMET CHAFF DEMAND, in grams and in days. The
! nitrogen and mineral handling status is printed as it stands.

! Record 109 only: the per-deck SEASON SUMMARY of the carbon-shortfall
! diagnostic. The row with an inversed rank of -1 is the TOTAL row, and the
! required-mass fraction is taken against THE DECK'S OWN season-maximum
! structural stem mass - no literature mass, density or lineic value appears as
! a denominator, a target or a bound anywhere here.
! A shortfall reported here is a MODEL-INTERNAL quantity: the difference
! between a bracketed literature length potential and what this model's own
! carbon paid for. It is not a measured discrepancy, it is not evidence of a
! defect in the axis count, tillering, photosynthesis or allocation, and it
! licenses no parameter change.
IF (height_internode_law == 1) THEN
  BLOCK
    REAL(wp) :: lsh, mreq_ax, mreq_m2, lsh_tot, mreq_ax_tot, mreq_m2_tot, mden
    INTEGER  :: ir_s
    OPEN (NEWUNIT=unit_ins, FILE='output/internode_shortfall.csv', STATUS='REPLACE', ACTION='WRITE')
    WRITE (unit_ins,'(A)') 'window_member,length_member,lineic_member,ear_member,onset_arm,dI_Cd,' // &
      'axes_live_m2,Mstruct_max_g_m2,ir,abs_rank,L_pot_final_cm,L_final_cm,L_short_cm,lambda_mg_cm,' // &
      'M_req_mg_axis,M_req_g_m2,M_req_frac,lambda_mod_mg_cm,m_struct_mg_axis,bind_final,' // &
      'h_b1_cm,h_b2_cm,h_b3_cm,L_ear_cm,e8_resid_max_g_m2,e8_over_max_g_m2,' // &
      'days_mass,days_pot,days_thicken,open_days'
    lsh_tot = 0.0_wp ; mreq_ax_tot = 0.0_wp ; mreq_m2_tot = 0.0_wp
    mden = MAX(inl_Mstruct_max, 1.0e-12_wp)
    DO ir_s = 0, INL_MAXIR
      IF (inl_Lf_use(ir_s) <= 0.0_wp) CYCLE
      lsh     = MAX(0.0_wp, inl_Lf_use(ir_s) - inl_L(ir_s))
      mreq_ax = inl_lam_use(ir_s) * lsh
      mreq_m2 = mreq_ax * inl_axes / 1000.0_wp
      lsh_tot = lsh_tot + lsh ; mreq_ax_tot = mreq_ax_tot + mreq_ax ; mreq_m2_tot = mreq_m2_tot + mreq_m2
      SELECT CASE (inl_bind(ir_s))
      CASE (1) ; inl_bind_txt = 'POT'
      CASE (2) ; inl_bind_txt = 'MASS'
      CASE (3) ; inl_bind_txt = 'CAP'
      CASE DEFAULT ; inl_bind_txt = 'CLOSED'
      END SELECT
      WRITE (unit_ins,'(5(I0,","),G14.6,",",G14.6,",",G14.6,",",I0,",",I0,9(",",G14.6),",",A,' // &
        '6(",",G14.6),4(",",I0))') &
        hgt_int_window_member, hgt_int_length_member, hgt_int_lineic_member, hgt_int_ear_member, &
        inl_onset_arm, inl_dI_use, inl_axes, inl_Mstruct_max, ir_s, inl_absrank(ir_s), &
        inl_Lf_use(ir_s), inl_L(ir_s), lsh, inl_lam_use(ir_s), mreq_ax, mreq_m2, mreq_m2 / mden, &
        inl_m(ir_s) / MAX(inl_L(ir_s), 1.0e-12_wp), inl_m(ir_s), inl_bind_txt, &
        inl_h_b1, inl_h_b2, inl_h_b3, inl_Lear_cm, inl_e8_resid_max, inl_e8_over_max, &
        inl_days_mass, inl_days_pot, inl_days_thicken, inl_open_days
    END DO
    WRITE (unit_ins,'(5(I0,","),G14.6,",",G14.6,",",G14.6,",",I0,",",I0,9(",",G14.6),",",A,' // &
      '6(",",G14.6),4(",",I0))') &
      hgt_int_window_member, hgt_int_length_member, hgt_int_lineic_member, hgt_int_ear_member, &
      inl_onset_arm, inl_dI_use, inl_axes, inl_Mstruct_max, -1, -2, &
      SUM(inl_Lf_use), SUM(inl_L), lsh_tot, 0.0_wp, mreq_ax_tot, mreq_m2_tot, mreq_m2_tot / mden, &
      SUM(inl_m) / MAX(SUM(inl_L), 1.0e-12_wp), SUM(inl_m), 'TOTAL', &
      inl_h_b1, inl_h_b2, inl_h_b3, inl_Lear_cm, inl_e8_resid_max, inl_e8_over_max, &
      inl_days_mass, inl_days_pot, inl_days_thicken, inl_open_days
    CLOSE (unit_ins)
    WRITE (*,'(A,4(A,I0),3(A,G14.6),A,G14.6,A,G14.6,4(A,I0))') &
      ' rec 109 INTERNODE SUMMARY', ' window=', hgt_int_window_member, ' length=', hgt_int_length_member, &
      ' lineic=', hgt_int_lineic_member, ' ear=', hgt_int_ear_member, &
      ' h_b1_cm=', inl_h_b1, ' h_b2_cm=', inl_h_b2, ' h_b3_cm=', inl_h_b3, &
      ' E8_resid_max_massday_g_m2=', inl_e8_resid_max, ' E8_overpay_max_g_m2=', inl_e8_over_max, &
      ' days_mass=', inl_days_mass, ' days_pot=', inl_days_pot, &
      ' days_thicken=', inl_days_thicken, ' open_days=', inl_open_days
  END BLOCK
END IF

! ==================================================================== !
! Record 86: the grain-nitrogen sink law at the end of the run, written only
! when the record is on and not on the critical-nitrogen-basis-only arm. Every
! envelope row is a FALSIFICATION target from a source read at first hand on
! this machine, never a tuning target.
! ==================================================================== !
IF (grainN_sink_mech >= 1 .AND. grainN_sink_mech /= 3) THEN
  BLOCK
    INTEGER :: unit_gnl
    REAL(wp) :: kn_f, ngrain_f, mg_kern, mg_ear, pot_mg_kern, pot_mg_ear, tot_src
    kn_f     = MAX(0.0_wp, kernel_number_m2)
    ngrain_f = MAX(0.0_wp, N_storage - N_chaff)
    mg_kern  = 0.0_wp ; IF (kn_f > 0.0_wp) mg_kern = 1000.0_wp * ngrain_f / kn_f
    mg_ear   = 0.0_wp ; IF (gnx_ears_m2 > 0.0_wp) mg_ear = 1000.0_wp * ngrain_f / gnx_ears_m2
    pot_mg_kern = 0.0_wp ; IF (kn_f > 0.0_wp) pot_mg_kern = 1000.0_wp * gnf_pot_cum / kn_f
    pot_mg_ear  = 0.0_wp ; IF (gnx_ears_m2 > 0.0_wp) pot_mg_ear = 1000.0_wp * gnf_pot_cum / gnx_ears_m2
    ! The GRAIN'S own sources. Under the released law (Eq. (11) p.785: the only
    ! source of nitrogen for the grains is the common pool of mobile nitrogen) the
    ! grain draws exactly two things, today's assimilation and the culm-common mobile
    ! pool. The lamina, stem and chaff rows below are POOL INFLOWS and not grain
    ! sources - they reach the grain only through the pool, and adding them here
    ! would double-count. Under mode 5 the organs ARE direct grain sources and the
    ! five-term sum is the right one.
    IF (gnf_eqn6) THEN
      tot_src = gnf_cum_assim + gnf_cum_mob
    ELSE
      tot_src = gnf_cum_assim + gnf_cum_mob + gnf_cum_lam + gnf_cum_stem + gnf_cum_chaff
    END IF
    OPEN (NEWUNIT=unit_gnl, FILE='output/grain_n_sink_law.csv', STATUS='REPLACE', ACTION='WRITE')
    WRITE (unit_gnl,'(A)') 'quantity,value,units_or_source'
    WRITE (unit_gnl,'("rec86_grainN_sink_mech,",I0,",0 OFF / 1 THE ITEM / 2 law-only arm /' // &
      ' 3 critical-N-basis-only arm / 4 law+uptake-demand arm / 5 mode-1-without-Eqn-6 arm = the' // &
      ' REJECTED v24 build, for attribution only")') grainN_sink_mech
    WRITE (unit_gnl,'("eqn6_synthesis_active,",L1,",FPB Eqn 6 lamina photosynthetic-N synthesis ON' // &
      ' (.T. on modes 1, 2 and 4; .F. on mode 5 = the v24 build)")') gnf_eqn6
    WRITE (unit_gnl,'(A)') '# --- the constants, ONE table, ONE paper (SHARED.f90) ---'
    WRITE (unit_gnl,'("n_g0_per_ear_mgN,",G14.6,",mg N ear-1 at anthesis; Bertheloot, Martre & Andrieu' // &
      ' 2008 Funct Plant Biol 35:781-796 doi 10.1071/FP08064 Table 1 p.784 (2.4e-6 kg) and p.786' // &
      ' ""grain N mass at anthesis, Ngtot(0), was set at 2.4 mg (Girard 1997)"". PER CULM = PER EAR:' // &
      ' FPB is a single-culm model.")') 1000.0_wp * gnf_ng0_per_ear
    WRITE (unit_gnl,'("g_rate_per_degCd,",G14.6,",[degCd]-1; FPB Table 1 p.784. *** DECLARED: NOT' // &
      ' MEASURED ***. FPB p.786: potNg(t=tau) was estimated as the slope of a linear regression on' // &
      ' treatment M15, then ""parameter g was calculated analytically from the values of Ngtot(t=0),' // &
      ' tau and potNgtot(t=250 degCd)"". g exists ONLY as a pair with 2.4 mg and 250 degCd and is NOT' // &
      ' refit against any observation here")') gnf_g_rate
    WRITE (unit_gnl,'("tau_cell_division_degCd,",G14.6,",degCd base 0 C; FPB Table 1 p.784, taken from' // &
      ' Gleadow et al. 1982 and Singh & Jenner 1982 (p.786) -- an EXTERNAL MEASUREMENT, unlike g")') &
      gnf_tau_degCd
    WRITE (unit_gnl,'("delta_organ_turnover_per_degCd,",G14.6,",[degCd]-1 relative rate of' // &
      ' photosynthetic-N degradation, FPB Eqn 5 + Table 1 p.784. UNCONDITIONAL on the LAMINA, paired' // &
      ' with Eqn 6 below exactly as FPB p.786 fitted it; a demand-gated POTENTIAL on stem and chaff,' // &
      ' for which FPB supplies no synthesis term -- see DECLARED_DEPARTURE")') n2_delta_TT
    WRITE (unit_gnl,'("sigma_lamina_synthesis_per_degCd,",G14.6,",[degCd]-1 relative rate of lamina' // &
      ' photosynthetic-N SYNTHESIS, FPB Eqn 6 + Table 1 p.784 ""sigma ... 0.0015"", read at first' // &
      ' hand. RESTORED 2026-09-21 after the independent audit REJECTED the v24 build for deleting' // &
      ' Eqn 6: FPB p.786 says sigma, delta, k1 and k2 (Eqns 5 AND 6) were estimated by minimising the' // &
      ' RMSE, so delta is one of FOUR JOINTLY FITTED parameters and is un-sourced without sigma. NOT' // &
      ' this tree n2_sigma_lam = 1.5e-4: that value''s Bertheloot 2011b source is NOT on this machine and its' // &
      ' misprint argument is refuted by FPB p.791 -- see SHARED.f90")') gnf_sigma_lam
    WRITE (unit_gnl,'("k1_Nmob_MM_gN_gDM,",G14.6,",kg kg-1; FPB Eqn 6 + Table 1 p.784")') n2_k1_Nmob
    WRITE (unit_gnl,'("k2_PAR_MM_J_m2_d,",G14.6,",J m-2 d-1 = FPB Table 1 p.784 k2 = 10 J m-2 s-1 x' // &
      ' 86400 s d-1, an EXACT unit conversion; the ratio PAR/(PAR+k2) is unchanged by it")') n2_k2_PAR
    WRITE (unit_gnl,'("r_mobile_to_photosynthetic_N,",G14.6,",dimensionless; FPB Eqn 1 p.783 gives' // &
      ' N^tot = (1+r)N^ph + N^struct, so N^ph = (N^tot - N^struct)/(1+r), and p.786 sets r = 0.10,' // &
      ' reached after the first time steps irrespective of the initial value. The v24 build omitted' // &
      ' the (1+r) and so degraded about 10 per cent MORE POOL than Eqn 5 addresses")') gnf_r_mob_ph
    WRITE (unit_gnl,'("struct_N_lamina_gN_gDM,",G14.6,",Bertheloot et al. 2008 Plant Physiol 148:1707' // &
      ' p.1711; the 0.0035 [degCd]-1 RATE with the same digits is on p.1712. PAGE RE-CHECKED BY PDF' // &
      ' PAGE INDEX 2026-09-21 -- item F 2026-09-20 changed this to p.1710 and the audit verified that' // &
      ' wrongly; both misread a running head in an interleaved two-column dump")') n2_struct_lam
    WRITE (unit_gnl,'("struct_N_sheath_internode_chaff_gN_gDM,",G14.6,",same page, same sentence")') &
      n2_struct_stem
    WRITE (unit_gnl,'("thermal_time_basis,base 0 C after anthesis,",A)') 'FPB p.782 ""Time is expressed' // &
      ' as thermal time after anthesis with a base temperature of 0 C"" == the model dTT_day_N' // &
      ' (cultivar T_Base 0.0). NO conversion, NO new constant'
    WRITE (unit_gnl,'(A)') '# --- what this deck actually did ---'
    WRITE (unit_gnl,'("ears_m2_from_item_A,",G14.6,",ears m-2 PREDICTED by rec 89/90; the per-ear' // &
      ' basis stands or falls with it and it is UNOBSERVED except at Maricopa")') gnx_ears_m2
      WRITE (unit_gnl,'("n_g0_gN_m2,",G14.6,",g N m-2 = 2.4e-3 * ears m-2")') gnf_ng0
    WRITE (unit_gnl,'("tt_post_anthesis_degCd,",G14.6,",degCd base 0 C at the last simulated day")') gnf_tt
    WRITE (unit_gnl,'("potential_grain_N_gN_m2,",G14.6,",g N m-2, potNg(t) = FPB Eqn 12 integrated")') &
      gnf_pot_cum
    WRITE (unit_gnl,'("realised_grain_N_gN_m2,",G14.6,",g N m-2 = N_storage - N_chaff")') ngrain_f
    WRITE (unit_gnl,'("fill_days_with_live_sink,",I0,",d")') gnf_days_fill
    WRITE (unit_gnl,'("days_SINK_bound,",I0,",d on which min() took the SINK term")') gnf_days_bind
    WRITE (unit_gnl,'("first_SINK_bound_dpa,",I0,",d after anthesis; -999 = never")') gnf_first_bind_dpa
    WRITE (unit_gnl,'("last_SINK_bound_dpa,",I0,",d after anthesis; -999 = never")') gnf_last_bind_dpa
    WRITE (unit_gnl,'("SINK_bound_days_after_15dpa,",I0,",d. FALSIFIER: Martre et al. 2003 Plant Physiol' // &
      ' 133:1959 p.1963 measured grain N ""sink limited or colimited ... for the first 10 to 15 d after' // &
      ' anthesis"" and ""always source limited during the grain filling"" thereafter. A deck sink-bound' // &
      ' through the linear phase FALSIFIES this parameterisation")') gnf_bind_after_15d
    WRITE (unit_gnl,'(A)') '# --- the source ladder, season totals (g N m-2) ---'
    WRITE (unit_gnl,'("from_assimilation,",G14.6,",g N m-2 the GRAIN took straight from today assimilation' // &
      ' (Jamieson & Semenov 2000 FCR 68:21 p.24 source order). Under the repaired law this EXCLUDES the' // &
      ' share FPB Eqn 6 routes to lamina photosynthetic-N synthesis, which is not a grain source")') gnf_cum_assim
    WRITE (unit_gnl,'("from_common_mobile_pool,",G14.6,",g N m-2 the GRAIN drew from the culm-common mobile' // &
      ' pool -- FPB Eqn 11 p.785: the ONLY other source of N for the grains. Filled by senescence' // &
      ' resorption, the rec-87 leaf-capacity export, and the organ inflows listed next")') gnf_cum_mob
    WRITE (unit_gnl,'("ROW_BASIS_organ_rows,",A)') 'Under the repaired law (rec 86 = 1, 2 or 4) the three' // &
      ' organ rows below are POOL INFLOWS, not grain sources: the organs reach the grain only through' // &
      ' the common pool, so source_total is assimilation + pool and adding the organs would DOUBLE-COUNT.' // &
      ' Under rec 86 = 5 (the rejected v24 build) the organs ARE direct grain sources and source_total' // &
      ' is the five-term sum.'
    WRITE (unit_gnl,'("from_lamina,",G14.6,",NET lamina export = Eqn 5 degradation MINUS Eqn 6' // &
      ' synthesis, above 0.0079 gN gDM-1. NEGATIVE means the lamina was a net N SINK over the fill,' // &
      ' which is a DECLARED FAILURE, not a result: FPB p.789 verbatim ""laminae N mass decreased in' // &
      ' all treatments, including Apache ... in all cases"". The p.791 ""or even increases"" defence' // &
      ' is WITHDRAWN (item E 2026-09-22): it is doubly conditioned on JUST AFTER ANTHESIS and on' // &
      ' significant post-anthesis root uptake, so it licenses a TRANSIENT, never a season total.")') gnf_cum_lam
    WRITE (unit_gnl,'("eqn5_lamina_degradation_total,",G14.6,",g N m-2, delta*N^ph*dTT summed' // &
      ' UNCONDITIONALLY over the fill (FPB Eqn 5 p.784)")') gnf_cum_deg_lam
    WRITE (unit_gnl,'("eqn6_lamina_synthesis_total,",G14.6,",g N m-2 returned to the lamina from the' // &
      ' common pool (FPB Eqn 6 p.784). ZERO on mode 5 and in the REJECTED v24 build")') gnf_cum_syn_lam
    WRITE (unit_gnl,'("from_stem,",G14.6,",delta-limited above 0.0035 gN gDM-1; see ROW_BASIS_organ_rows")') gnf_cum_stem
    ! This row prints the chaff structural-nitrogen floor IN FORCE, so that the
    ! published diagnostic can never disagree with the executable: record 178 = 1
    ! makes the chaff limb apply its OWN measured floor 5.10e-3 g N g-1 DM
    ! (Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109 Table 1(c)
    ! column Chaff, doi:10.1093/aob/mcr125) instead of the stem value 0.0035.
    WRITE (unit_gnl,'("from_chaff,",G14.6,",delta-limited above ",ES10.3," gN gDM-1 IN FORCE' // &
      ' (rec 178 = 0: the sheath/internode n2_struct_stem = 3.5e-3; rec 178 = 1: the chaff own' // &
      ' measured n2_struct_chaff = 5.10e-3, Bertheloot 2011b Table 1(c) column Chaff; rec 178 = 2: Pask 2012' // &
      ' chaff structural minimum 3.80e-3 PLUS the NRE floor (1 - 0.56) x N_chaff(anthesis); rec 178' // &
      ' = 3: the same NRE law on the 5.10e-3 floor). THIS ROUTE DID NOT' // &
      ' EXIST BEFORE THE GRAIN-N SINK LAW: N_chaff was FROZEN from anthesis to maturity on 22 of 22 baseline decks,' // &
      ' stranding 0.6-3.6 g N m-2 that PP148 p.1711 measures as mobile; see ROW_BASIS_organ_rows")') &
      gnf_cum_chaff, n2_struct_chaff_use
    WRITE (unit_gnl,'("source_total,",G14.6,",g N m-2, must equal realised grain N to machine precision")') tot_src
    WRITE (unit_gnl,'(A)') '# --- measured envelopes: falsification targets, never tuning targets ---'
    WRITE (unit_gnl,'("N_per_kernel_mg,",G14.6,",mg N kernel-1 REALISED. ENVELOPE 0.75-1.38 mg,' // &
      ' Martre et al. 2003 Plant Physiol 133:1959 Table I p.1961 (0.75-1.11 intact ears; 1.38 with' // &
      ' ablation). A PLAUSIBILITY band, not a target")') mg_kern
    WRITE (unit_gnl,'("N_per_ear_mg,",G14.6,",mg N ear-1 REALISED. ENVELOPE 18.9-42.6 mg, Bertheloot' // &
      ' 2008 FPB p.787 observed final grain N over 11 treatments")') mg_ear
    WRITE (unit_gnl,'("potential_N_per_kernel_mg,",G14.6,",mg N kernel-1 the LAW allows")') pot_mg_kern
    WRITE (unit_gnl,'("potential_N_per_ear_mg,",G14.6,",mg N ear-1 the LAW allows; MUST sit inside' // &
      ' 18.9-42.6 on a healthy deck or the (2.4 mg, g, tau) triple is falsified")') pot_mg_ear
    WRITE (unit_gnl,'("organ_state_N_closure_max_gN_m2_d,",G14.6,",max |dN_leaf + dN_stem + dN_storage' // &
      ' + dN_mobile - N_assimilated| over the run, read from ORGAN STATE captured on entry -- NOT the' // &
      ' flux identity re-arranged. A dropped decrement or a double-counted import surfaces at its own' // &
      ' size. Written because the item-C auditor showed the old carbon check was a tautology")') gnf_nres_max
    WRITE (unit_gnl,'("DECLARED_DEPARTURE,",A)') 'CONFINED TO STEM AND CHAFF as of the 2026-09-21' // &
      ' Eqn-6 repair. The LAMINA now runs FPB Eqn 4 in full (Eqn 5 degradation, unconditional, PLUS' // &
      ' Eqn 6 synthesis), so delta acts inside the calibration that produced it. Stem and chaff keep' // &
      ' the demand-gated delta potential because FPB models the stem as a FORCED variable (Eqn 14' // &
      ' p.784) and gives it no sigma; a sigma for them exists only in Bertheloot 2011b, which is NOT on this' // &
      ' machine. Making their degradation unconditional with no synthesis limb would repeat the very' // &
      ' error being repaired, one organ over. CONSEQUENCE: where the sink binds, STEM and CHAFF N is' // &
      ' retained rather than degraded; the LAMINA is no longer stripped bare.'
    WRITE (unit_gnl,'("NOT_MODELLED,",A)') 'any CO2 term in the grain-N law -- FPB Eqn 11-12 has none,' // &
      ' so this model reproduces ONLY the growth-DILUTION half of the eCO2 protein decline. Pleijel &' // &
      ' Uddling 2012 Glob Change Biol 18:596-605 p.599 measured, over 57 experiments, y = 0.75x - 0.072' // &
      ' with intercept 0.072 +- 0.038 P < 0.001: at ZERO CO2-induced yield change grain protein yield' // &
      ' still falls about 7 %. THE YIELD-INDEPENDENT HALF IS ABSENT BY CONSTRUCTION and is NOT absorbed' // &
      ' into a tuned sink. Also not modelled: FPB Eqn 6 for the STEM and CHAFF (no sigma exists for' // &
      ' them in any source on this machine); FPB p.786 [N]mob = 1e-4 kg kg-1 grain cut-off, which that' // &
      ' page introduces expressly to avoid TRANSITORY NEGATIVE N_c^mob -- a solver guard, not a' // &
      ' physiological floor, and the daily pool here cannot go negative; any grain protein' // &
      ' FRACTION (albumin/gliadin/glutenin, Martre 2003); nitrate-assimilation inhibition (Bloom 2014).'
    WRITE (unit_gnl,'("RETIRED_BY_THIS_ITEM,",A)') 'the grain taking 100 % of N_assim_avail and 100 % of' // &
      ' N_mobile_veg unconditionally; the organ draw paced by n3_eff*N_organ(0)*prog with prog linear in' // &
      ' dvs2 and the Pask 2012 efficiency caps 0.76/0.48; and (on modes 1 and 4) the Jamieson & Semenov' // &
      ' 2000 NHI trajectory AS THE POST-ANTHESIS UPTAKE DEMAND. J&S is NOT called wrong: on their own' // &
      ' data H1 BEAT the per-grain sink H2 (R2 0.95 / RMSD 0.18 N% vs 0.19 / 0.88, pp.27-28). What is' // &
      ' retired is its use as a PUMP; inside min(sink, source) the sink is only a CEILING.'
    CLOSE(unit_gnl)
  END BLOCK
END IF

! ==================================================================== !
! Record 115: the per-kernel grain-nitrogen sink basis at the end of the run,
! ON ONLY, so that this file ABSENCE is itself part of the comparison. Every
! cohort row is printed so that the three-branch identity and the
! ear-fraction cancellation of record 90 mode 4 are re-derivable FROM THIS FILE
! ALONE.
! ==================================================================== !

! (flag canopy_sink_lai, rec 69): the end-of-run INVARIANT 6
! audit. Written only when the flag is on, so the OFF run is unchanged.
IF (csl_diag_open) THEN
  CLOSE(unit_csl)
  OPEN (newunit=unit_csl,file='output/canopy_sink_summary.csv',status='replace',action='write')
  WRITE (unit_csl,'(A)') 'quantity,value,units'
  WRITE (unit_csl,'("LMA_used,",G14.6,",g DM m-2")')                   csl_LMA
  WRITE (unit_csl,'("LMA_implied_SLA,",G14.6,",m2 kg-1")')             MERGE(1000.0_wp/csl_LMA, 0.0_wp, csl_LMA > 0.0_wp)
  WRITE (unit_csl,'("lamina_area_built_total,",G14.6,",m2 m-2")')      csl_area_built
  WRITE (unit_csl,'("lamina_C_demanded_total,",G14.6,",g DM m-2")')    csl_area_built * csl_LMA
  WRITE (unit_csl,'("backlog_end,",G14.6,",g DM m-2")')                csl_backlog
  WRITE (unit_csl,'("backlog_max,",G14.6,",g DM m-2")')                csl_backlog_max
  WRITE (unit_csl,'("limited_days,",I0,",days")')                      csl_limited_days
  WRITE (unit_csl,'("fln_cut,",I0,",leaves")')                         csl_fln_cut
  WRITE (unit_csl,'("cohorts_aborted,",I0,",count")')                  csl_cohorts_aborted
  WRITE (unit_csl,'("births_suspended,",I0,",days")')                  csl_births_suspended
  WRITE (unit_csl,'("partition_sum_err_max,",G14.6,",-")')             csl_fsum_err
  WRITE (unit_csl,'("realized_SLA_end,",G14.6,",m2 kg-1")')            csl_sla_real
  WRITE (unit_csl,'("LAI_sink_end,",G14.6,",m2 m-2")')                 csl_lai_sink
  WRITE (unit_csl,'("PAI_diag_end,",G14.6,",m2 m-2")')                 csl_pai_diag
  WRITE (unit_csl,'("lai_basis_Ncapacity_end,",G14.6,",m2 m-2")')      csl_lai_Ncap
  ! ---- INPUTS FOR THE INSTRUMENT FORWARD MODEL AND ITS INVERSION -------- !
  ! The per-element-class standing areas at DAILY resolution are already emitted
  ! in output/daily_canopy_sheath_ear.csv: green blade, standing-dead blade,
  ! green sheath, standing-dead sheath, bare internode, and ear plus peduncle,
  ! all in m2 m-2 of ground. What that file cannot carry is the ANGULAR
  ! information, because THE MODEL HAS NONE; it is stated here per class rather
  ! than left to the reader to assume. This block records the geometry and the
  ! inclination ASSUMPTIONS, so that an inversion has the complete input and
  ! knows exactly what is measured, what is assumed and what is absent.
  WRITE (unit_csl,'(A)') ''
  WRITE (unit_csl,'(A)') 'inversion_input,value,units_or_status'
  WRITE (unit_csl,'(A)') 'row_spacing,NOT A MODEL VARIABLE,' // &
    'the model is 1-D and carries no row geometry; the SITE value ' // &
    '(0.175 m rows) is a site input (see INPUT_GUIDE.md)'
  WRITE (unit_csl,'(A)') 'inclination_blade_green,NONE IN MODEL,' // &
    'areas are one-sided projected areas and no angular distribution is ' // &
    'carried. MEASURED SOURCE AVAILABLE BUT NOT IMPLEMENTED: Abichou 2016 ' // &
    'Ch.4 Figs 5-7 (leaf inclination and insertion angle vs phyllochronic age)'
  WRITE (unit_csl,'(A)') 'inclination_blade_dead,NONE IN MODEL,' // &
    'as green; no post-senescence reorientation is carried'
  WRITE (unit_csl,'(A)') 'inclination_sheath,VERTICAL BY CONSTRUCTION,' // &
    'a sheath wraps the culm; entered as a cylinder and converted to half ' // &
    'total surface (pi/2 x L x D) on the LI-2200C basis'
  WRITE (unit_csl,'(A)') 'inclination_internode_bare,VERTICAL BY CONSTRUCTION,' // &
    'same cylinder treatment as the sheath'
  WRITE (unit_csl,'(A)') 'inclination_ear_peduncle,VERTICAL BY CONSTRUCTION,' // &
    'A_ep is a projected area on the same convention; the AWN STATUS of the ' // &
    'model ear is UNDEFINED because A_ep is a declared gap'
  WRITE (unit_csl,'(A)') 'clumping,NOT MODELLED,' // &
    'the model reports a GEOMETRIC total area; no clumping index is applied'
  CLOSE(unit_csl)
  csl_diag_open = .FALSE.
END IF

! Soil-heat guard QA summary: always written under mode 2
! so QA can assert an explicit zero.
IF (energy_balance_mode == 2) THEN
  BLOCK
    INTEGER :: unit_shfb
    OPEN (NEWUNIT=unit_shfb, FILE='output/soil_heat_fallback.csv', STATUS='REPLACE', ACTION='WRITE')
    WRITE (unit_shfb,'(A)') 'soil_heat_fallback_events,soil_heat_fallback_days,' // &
      'eb_noroot_events,eb_noroot_days,eb_seed_resets,eb_fluxfb_events,eb_fluxfb_days,' // &
      'hydro_bcswitch_events,hydro_bcswitch_days'
    WRITE (unit_shfb,'(I0,8(",",I0))') soil_heat_fb_events, soil_heat_fb_days, &
      eb_noroot_events, eb_noroot_days, eb_seed_resets, eb_fluxfb_events, eb_fluxfb_days, &
      hydro_bcswitch_events, hydro_bcswitch_days
    CLOSE (unit_shfb)
  END BLOCK
END IF

!==============================================================================!
CONTAINS

! ---------------------------------------------------------------------------- !
! Refuse a driver-record value whose code the release does not contain. Prints
! the record number, the record variable name, the offending value and the
! implemented set, then stops: a value whose branch was removed never reaches a
! mechanism.
! ---------------------------------------------------------------------------- !
SUBROUTINE reject_retired_option (rec_no, rec_var, got, allowed_txt)
  INTEGER,          INTENT(IN) :: rec_no
  CHARACTER(LEN=*), INTENT(IN) :: rec_var
  INTEGER,          INTENT(IN) :: got
  CHARACTER(LEN=*), INTENT(IN) :: allowed_txt
  WRITE (*,'(A)') ' ------------------------------------------------------------'
  WRITE (*,'(A,I0,A,A,A,I0)') ' REFUSING TO RUN: driver_main.txt record ', rec_no, &
    ' (', TRIM(rec_var), ') = ', got
  WRITE (*,'(A,I0,A,A)') '   M&M-CO2Wheat v1.0 implements only these values of record ', &
    rec_no, ': ', TRIM(allowed_txt)
  WRITE (*,'(A)') '   The code of every other value of this record was removed from the'
  WRITE (*,'(A)') '   v1.0 release. It is refused here so that it can never fall through'
  WRITE (*,'(A)') '   to a different mechanism or to an undefined value.'
  WRITE (*,'(A)') '   This option is not available in v1.0 (see Part 1, Table 6).'
  WRITE (*,'(A)') ' ------------------------------------------------------------'
  ERROR STOP 'retired driver option: see the record number printed above'
END SUBROUTINE reject_retired_option

  ! The root and shoot nitrate-reduction block, as an internal procedure so that
  ! it can run at either of two points of the day; by host association every name
  ! is the program's own. Called once per day.
  SUBROUTINE andrews_nitrate_partition_day()
  ! -------------------------------------------------------------------- !
  ! THE ROOT AND SHOOT NITRATE PARTITION. The mass-balance foundation is
  ! retained: soil uptake is the ONLY entry of nitrogen atoms; it feeds the
  ! inorganic root pool, and all reduction draws from the two pools, so the
  ! assimilation can never exceed the uptake and no nitrogen comes from nowhere.
  ! The two pools are the substrate for the transpiration mass-flow uptake and
  ! for the photorespiratory shoot reduction.
  !
  ! Andrews (1986) Plant Cell Environ 9:511-519: the root and shoot split of
  ! nitrate assimilation is set by the external nitrate concentration (p.511
  ! abstract and p.515, a constant partitioning is achieved by balancing nitrate
  ! reductase levels in the root with the rates of nitrate uptake); when the
  ! uptake outruns the root reductase, the unassimilated overflow is passed into
  ! the xylem to the shoot. Their Table 1b, Triticum aestivum: the shoot share is
  ! about 50 % on in-vivo reductase activity, rising to 65-83 % on xylem sap.
  ! Andrews et al. (1992) Ann Bot 70:271-282: the cereal transition is at 1 to
  ! 5 mol m-3, so the half-saturation of the partition is about 3.
  !   f_shoot = f_min + (f_max - f_min) NO3_rz / (K + NO3_rz)
  !
  ! Root reduction is the photorespiration-INDEPENDENT baseline, capped by the
  ! root reductase ceiling, the specific rate times the root mass (the reductase
  ! level balanced with the uptake of Andrews; the specific rate is derived in
  ! driver_N.csv). Shoot reduction is capped by the photorespiratory capacity.
  ! Ordering: the xylem overflow is taken first, the shoot share of the root
  ! pool, then the root reduces the RETAINED share up to its ceiling and the
  ! remainder carries over in the root pool - which is mass-conservative and
  ! avoids the double draw a literal reading would produce, where both routes
  ! draw the full pool. The photorespiratory flux itself stays a DIAGNOSTIC.
  ! -------------------------------------------------------------------- !
  N_NO3_root = N_NO3_root + N_soil_uptake_day          ! uptake S -> root pool (Gojon B substrate)
  BLOCK
    REAL(wp) :: NO3_rz, f_shoot, R_root, N_xylem
    ! Root-zone nitrate concentration (mol m-3): g N m-2 / (water depth, m) / 14 g mol-1.
    ! vwc_root (line ~1146) is the current root-zone mean theta; root_dep the current depth.
    IF (soilN_layered == 1) THEN
      ! Layered: the concentration the root system senses is the ROOT-ZONE AGGREGATE
      ! of the layered nitrate, the layer pools weighted by their overlap with the
      ! root depth, the same convention as the denitrification argument. The
      ! partition LAW is unchanged.
      BLOCK
        REAL(wp) :: zt3, ov3, no3_rz_sum
        INTEGER  :: iLa
        zt3 = 0.0_wp; no3_rz_sum = 0.0_wp
        DO iLa = 1, NSoil_layers
          ov3 = MAX(0.0_wp, MIN(zt3 + tk(iLa), MAX(0.1_wp, root_dep)) - zt3)
          IF (tk(iLa) > 1.0e-9_wp) no3_rz_sum = no3_rz_sum + (ov3 / tk(iLa)) * NO3_lay(iLa)
          zt3 = zt3 + tk(iLa)
        END DO
        NO3_rz = no3_rz_sum / MAX(1.0e-6_wp, vwc_root * MAX(0.1_wp, root_dep)) / 14.0_wp
      END BLOCK
    ELSE
    NO3_rz  = N_soil_min / MAX(1.0e-6_wp, vwc_root * MAX(0.1_wp, root_dep)) / 14.0_wp
    END IF
    f_shoot = f_shoot_min + (f_shoot_max - f_shoot_min) * NO3_rz / (K_partition_NO3 + NO3_rz)  ! Andrews 1986/1992
    N_xylem     = f_shoot * N_NO3_root                 ! overflow -> xylem -> shoot (Andrews p.515)
    N_NO3_root  = N_NO3_root  - N_xylem
    N_NO3_shoot = N_NO3_shoot + N_xylem
    R_root = r_root_NR * W_root                        ! root-NR ceiling ~ root DM (Andrews p.511/515)
    N_assim_root_day  = MIN(N_NO3_root, R_root)        ! photoresp-INDEPENDENT root reduction
    N_NO3_root  = N_NO3_root  - N_assim_root_day        ! unreduced nitrate carries over in pool
    ! The shoot reduction is capped by the photorespiratory capacity, the glycine
    ! and serine flux, which is CO2-sensitive through the photorespiratory ratio
    ! (Busch, Sage & Farquhar 2018 Eqs. 16-21 and Fig. 3b). Unassimilated nitrate
    ! stays in the shoot pool, which is the measured behavior: shoot nitrate
    ! assimilation is inhibited at elevated CO2 and nitrate accumulates
    ! (Rachmilevitch et al. 2004; Bloom et al. 2010). The minimum also guarantees
    ! mass conservation.
    N_assim_shoot_day = MIN(N_NO3_shoot, N_assim_biochem_day)  ! Process 3 (Gojon C)
    N_NO3_shoot = N_NO3_shoot - N_assim_shoot_day
  END BLOCK
  N_gain_total = N_assim_root_day + N_assim_shoot_day  ! <= N_soil_uptake_day by construction
  IF (N_gain_total > 0.0_wp) THEN
    f_photoresp_diag = N_assim_shoot_day / N_gain_total ! EMERGENT shoot/total (Andrews range ~0.5-0.75)
  ELSE
    f_photoresp_diag = 0.0_wp
  END IF
  N_assim_biochem_retained_day = 0.0_wp                ! Busch flux retained as phantom N: REMOVED (P1)
  END SUBROUTINE andrews_nitrate_partition_day

  ! Add one organic-treatment event carbon and nitrogen into the fast and slow
  ! residue pools. The precedence is: dry mass in kg ha-1 times 0.1 gives g m-2;
  ! the carbon and nitrogen percentages are divided by 100; a negative nitrogen
  ! percentage means derive the nitrogen from the carbon input and the default
  ! carbon-to-nitrogen ratio, and a negative carbon percentage means use the
  ! default residue carbon fraction.
  SUBROUTINE apply_residue_event(i, n_added)
    INTEGER, INTENT(IN)   :: i
    REAL(wp), INTENT(OUT) :: n_added
    REAL(wp) :: mass_g, Cfrac, C_in, N_in, ff
    mass_g = N_org_mass(i) * 0.1_wp
    IF (N_org_Cpct(i) >= 0.0_wp) THEN
      Cfrac = N_org_Cpct(i) / 100.0_wp
    ELSE
      Cfrac = residue_C_frac
    END IF
    C_in = mass_g * Cfrac
    IF (N_org_Npct(i) >= 0.0_wp) THEN
      N_in = mass_g * N_org_Npct(i) / 100.0_wp
    ELSE
      N_in = C_in / MAX(1.0_wp, N_org_CN(i))
    END IF
    ff = N_org_fast(i)
    res_fast_C = res_fast_C + C_in * ff
    res_fast_N = res_fast_N + N_in * ff
    res_slow_C = res_slow_C + C_in * (1.0_wp - ff)
    res_slow_N = res_slow_N + N_in * (1.0_wp - ff)
    n_added = N_in
  END SUBROUTINE apply_residue_event

  INTEGER FUNCTION count_fields(line)
    CHARACTER(LEN=*), INTENT(IN) :: line
    INTEGER :: i, ncommas, llen
    llen = LEN_TRIM(line)
    IF (llen == 0) THEN
      count_fields = 0
      RETURN
    END IF
    ncommas = 0
    DO i = 1, llen
      IF (line(i:i) == ',') ncommas = ncommas + 1
    END DO
    count_fields = ncommas + 1
  END FUNCTION count_fields

  INTEGER FUNCTION detect_ninvar(fname)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    CHARACTER(LEN=1024) :: line
    INTEGER :: ios_local, unit_local
    REAL(wp) :: tmp
    detect_ninvar = -1
    OPEN (newunit=unit_local, FILE=TRIM(fname), STATUS='OLD', IOSTAT=ios_local)
    IF (ios_local /= 0) RETURN
    READ (unit_local,'(A)',IOSTAT=ios_local) line
    IF (ios_local /= 0) THEN
      CLOSE(unit_local)
      RETURN
    END IF
    READ (line,*,IOSTAT=ios_local) tmp
    IF (ios_local /= 0) THEN
      READ (unit_local,'(A)',IOSTAT=ios_local) line
      IF (ios_local /= 0) THEN
        CLOSE(unit_local)
        RETURN
      END IF
    END IF
    detect_ninvar = count_fields(line)
    CLOSE(unit_local)
  END FUNCTION detect_ninvar

  SUBROUTINE parse_csv_reals(line, nvals, vals)
    CHARACTER(LEN=*), INTENT(IN) :: line
    INTEGER, INTENT(IN) :: nvals
    REAL(wp), INTENT(OUT) :: vals(nvals)
    INTEGER :: pos, start, lenl, idx, ios
    CHARACTER(LEN=128) :: token
    vals = 0.0_wp
    lenl = LEN_TRIM(line)
    start = 1
    idx = 1
    DO pos = 1, lenl + 1
      IF (pos == lenl + 1 .OR. line(pos:pos) == ',') THEN
        IF (idx <= nvals) THEN
          IF (pos > start) THEN
            token = ADJUSTL(line(start:pos-1))
          ELSE
            token = ''
          END IF
          IF (LEN_TRIM(token) == 0) THEN
            vals(idx) = 0.0_wp
          ELSE
            READ (token,*,IOSTAT=ios) vals(idx)
            IF (ios /= 0) vals(idx) = 0.0_wp
          END IF
        END IF
        idx = idx + 1
        start = pos + 1
      END IF
    END DO
  END SUBROUTINE parse_csv_reals

  FUNCTION count_csv_vals(line) RESULT(nvals)
    CHARACTER(LEN=*), INTENT(IN) :: line
    INTEGER :: nvals
    INTEGER :: i, lenl
    LOGICAL :: in_field
    nvals = 0
    lenl = LEN_TRIM(line)
    in_field = .FALSE.
    DO i = 1, lenl
      IF (line(i:i) == ',') THEN
        IF (in_field) nvals = nvals + 1
        in_field = .FALSE.
      ELSEIF (line(i:i) /= ' ') THEN
        in_field = .TRUE.
      END IF
    END DO
    IF (in_field) nvals = nvals + 1
  END FUNCTION count_csv_vals

END PROGRAM mmco2crop
!==============================================================================!