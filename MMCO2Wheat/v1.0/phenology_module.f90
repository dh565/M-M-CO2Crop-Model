!===============================================================================!
! M&M-CO2Wheat v1.0  |  phenology_module.f90
!
! Purpose: Wheat development on one thermal clock: the daily thermal-time
!          increment and its vernalization and photoperiod factors, the
!          cumulative thermal development units, the dates of the eight
!          development stages, and the continuous 0-2 development stage that
!          the growth, grain-filling and senescence laws run on. Time step:
!          daily (one call per simulated day). Spatial unit: the whole plant,
!          the stages being properties of the apex.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision) only. The module reads the cultivar parameter
!       file phenology_cultivar.txt through its own reader and calls no routine
!       outside itself; the daily driving temperatures are supplied by the
!       caller.
! Output: the derived types PhenologyParams (genotype thresholds, cardinal
!         temperatures, sensitivities, sowing date, observed stage dates, the
!         initial organ mass and nitrogen fractions and the optional cultivar
!         canopy-architecture traits) and PhenologyState (cum_tdu, cum_vd,
!         current_dtt, current_dtdu, current_vf, current_df, photoperiod, dap,
!         istage, stage_dap, stage_reached, tdu_at_emerge, ts_done); the
!         functions calc_dev_stage_0_2 and calc_dev_stage_rate_0_2, whose 0-2
!         development stage and its daily increment are read by growth_module
!         (partition fractions, specific leaf area, root depth, senescence) and
!         by the grain-filling and canopy clocks in main_MMCO2Crop; calc_dtt
!         and TS_TDU_MCMASTER, exported for the optional soil-temperature
!         canopy clock.
! Called by: main_MMCO2Crop (PROGRAM mmco2crop), which calls
!            initialize_phenology_params, read_phenology_cultivar,
!            initialize_phenology_state, bootstrap_phenology_state and, once
!            per day, update_phenology_day.
!===============================================================================!
MODULE PHENOLOGY_MODULE
  USE TYPES, ONLY: wp
  IMPLICIT NONE
  PRIVATE
  ! --------------------------------------------------------------------------- !
  ! Development follows one thermal clock (McMaster et al., 2008, Annals of
  ! Botany 102(4):561-569). The daily thermal development unit is
  ! DTDU = DTT * min(VF, DF) inside the emergence-to-terminal-spikelet window
  ! and DTT outside it (Eq. (1), p. 563; p. 565 "DF affects development from
  ! emergence to terminal spikelet"; p. 564, vernalization days "begin at
  ! germination and proceed until ... terminal spikelet initiation"). The
  ! temperature response, the vernalization factor VF and the photoperiod
  ! factor DF are given with the functions that compute them below. Eight
  ! stages are declared when the cumulative thermal development units pass the
  ! genotype thresholds of the cultivar file.
  ! Pre-emergence and post-emergence thermal forcing are passed in separately
  ! by the caller; this module owns no weather input.
  ! --------------------------------------------------------------------------- !
  INTEGER, PARAMETER, PUBLIC :: NUM_PHENO_STAGES = 8
  REAL(wp), PARAMETER :: PI_PHENO = 3.14159265358979323846_wp
  ! Sowing depth (cm) above which no subcrown internode has to be elongated
  ! before the coleoptile reaches the surface (Wang et al., 2009).
  REAL(wp), PARAMETER :: DEPTH_THRESH_CM = 2.5_wp
  ! McMaster et al. (2008) p. 563: "Terminal spikelet is reached when TDU reaches 400 degCd",
  ! counted from seedling emergence, which is also the start of the VF/DF window.
  REAL(wp), PARAMETER :: TS_TDU_MCMASTER = 400.0_wp

  TYPE, PUBLIC :: PhenologyParams
    INTEGER :: gdd_method = 1
    REAL(wp) :: t_base = 0.0_wp
    REAL(wp) :: t_opt_lower = 25.0_wp
    REAL(wp) :: t_opt_upper = 50.0_wp
    REAL(wp) :: t_max = 60.0_wp
    REAL(wp) :: photo_sens = 0.0_wp
    REAL(wp) :: vern_factor = 0.0_wp
    REAL(wp) :: stage_thresh(NUM_PHENO_STAGES) = 0.0_wp
    CHARACTER(LEN=20) :: stage_name(NUM_PHENO_STAGES)
    INTEGER :: sowing_year = 1900
    INTEGER :: sowing_month = 1
    INTEGER :: sowing_day = 1
    INTEGER :: obs_stage_dap(NUM_PHENO_STAGES) = 0
    INTEGER :: obs_two_leaves_dap = 0
    REAL(wp) :: init_green_g_plant = 0.020_wp
    REAL(wp) :: init_stem_g_plant = 0.012_wp
    REAL(wp) :: init_root_g_plant = 0.015_wp
    REAL(wp) :: init_storage_g_plant = 0.0_wp
    REAL(wp) :: init_leafN_frac = 0.050_wp
    REAL(wp) :: init_stemN_frac = 0.020_wp
    REAL(wp) :: init_rootN_frac = 0.020_wp
    REAL(wp) :: init_storageN_frac = 0.015_wp
    LOGICAL :: use_beta_emerge = .FALSE.
    REAL(wp) :: sowing_depth_cm = 5.0_wp
    REAL(wp) :: beta_tb = 0.0_wp
    REAL(wp) :: beta_tm = 39.0_wp
    REAL(wp) :: mu_germ = -5.32_wp
    REAL(wp) :: alpha_germ = 1.42_wp
    REAL(wp) :: beta_germ = 0.178_wp
    REAL(wp) :: mu_cole = -1.663_wp
    REAL(wp) :: alpha_cole = 0.9526_wp
    REAL(wp) :: beta_cole = 0.3335_wp
    REAL(wp) :: a_sub = 0.0217_wp
    REAL(wp) :: b_sub = -0.0015_wp
    REAL(wp) :: alpha_sub = 1.8075_wp
    REAL(wp) :: beta_sub = 0.7899_wp
    ! ------------------------------------------------------------------------- !
    ! Three OPTIONAL cultivar canopy-architecture rows, read from the end of
    ! phenology_cultivar.txt after Use_Beta_Emerge and its own optional
    ! beta-model block, when that is present. A row absent at end of file
    ! leaves the default declared here, which is the same literal as the
    ! compiled default in SHARED.f90 (clr_phyllochron_Cd, cba_l_irank_ref(2),
    ! cba_w_irank_ref(1:2)). All three are cultivar-scoped: the same value at
    ! every site and season for a given cultivar, never a site- or
    ! treatment-specific multiplier.
    ! ------------------------------------------------------------------------- !
    REAL(wp) :: phyllochron_degCd = 73.9_wp  ! degree-Cd; default = McMaster et al. 2003 Ann Bot 91:697-705 Table 2 p.701 (DOI 10.1093/aob/mcg074) -- see SHARED.f90 clr_phyllochron_Cd
    REAL(wp) :: blade_l_max_cm = 20.8_wp     ! cm; default = Abichou 2016 Ch.3 Fig. 8a doc p.104 profile maximum -- see SHARED.f90 cba_l_irank_ref(2)
    REAL(wp) :: blade_w_max_cm = 1.45_wp     ! cm; default = Abichou 2016 Ch.2 Fig. 3e doc p.67 profile maximum -- see SHARED.f90 cba_w_irank_ref(1:2)
    ! ------------------------------------------------------------------------- !
    ! Two FURTHER optional rows, read immediately after Blade_W_max_cm with the
    ! same end-of-file default: a row absent leaves the literal declared here,
    ! identical to the compiled SHARED.f90 default. Sheath and internode
    ! geometry is a cultivar canopy-architecture trait, not site- or
    ! treatment-specific. Sheath_D_cm is read only by update_canopy_sheath_ear
    ! (growth_module.f90), which is diagnostic: it writes
    ! daily_canopy_sheath_ear.csv and enters no carbon or nitrogen ledger, no
    ! partition fraction, no flux and no leaf-area identity. Sheath_L_max_cm is
    ! also read by the stem-demand and internode-elongation-height mechanisms,
    ! which are switched off in every deck of this release; see the SHARED.f90
    ! cse_ls_irank_ref note. The internode profile cse_li_irank is deliberately
    ! not exposed here, for the reason given in the SHARED.f90 cse_li_irank
    ! note.
    ! ------------------------------------------------------------------------- !
    REAL(wp) :: sheath_l_max_cm = 17.0_wp  ! cm; default = Abichou 2016 Ch.2 Fig. 3b doc p.67 profile maximum -- see SHARED.f90 cse_ls_irank_ref(1)
    REAL(wp) :: sheath_d_cm     = 0.45_wp  ! cm; default = Abichou 2016 Ch.2 Table S2 doc p.87 median upper-sheath diameter -- see SHARED.f90 cse_ds2_ref
    ! ------------------------------------------------------------------------- !
    ! A SIXTH optional row: tillering capacity as a genotype trait, the maximum
    ! number of live axes per plant (axes plant-1). 0 = no cap and is the
    ! end-of-file default. The crop-scale bound is this value times the deck's
    ! own plant density, so the 600-800 axes m-2 at 150-250 plants m-2 of
    ! Abichou (2016, Ch. 3, Fig. 3) maps onto 2.4-5.3 axes plant-1. See
    ! SHARED.f90 Axes_max_per_plant.
    ! ------------------------------------------------------------------------- !
    REAL(wp) :: axes_max_per_plant = 0.0_wp  ! axes plant-1; 0 = no cap (identity)
    ! ------------------------------------------------------------------------- !
    ! Development-factor window, driver record 150, set by main_MMCO2Crop after
    ! the cultivar file is read. 1 selects McMaster et al. (2008) as cited,
    ! DTT*min(VF,DF) acting only from emergence to terminal spikelet, and is
    ! the value shipped with v1.0; with 0 no vernalization, photoperiod or
    ! thermal-time accumulation takes place in update_phenology_day at all.
    ! ------------------------------------------------------------------------- !
    INTEGER :: dev_factor_window = 0
  END TYPE PhenologyParams

  TYPE, PUBLIC :: PhenologyState
    INTEGER  :: planting_doy = 1
    INTEGER  :: planting_year = 1900
    REAL(wp) :: latitude_deg = 0.0_wp
    REAL(wp) :: cum_tdu = 0.0_wp
    REAL(wp) :: cum_vd = 0.0_wp
    REAL(wp) :: current_dtt = 0.0_wp
    REAL(wp) :: current_dtdu = 0.0_wp   ! today's thermal development-unit increment
    REAL(wp) :: current_vf = 1.0_wp
    REAL(wp) :: current_df = 1.0_wp
    REAL(wp) :: current_photoperiod_h = 0.0_wp
    INTEGER  :: dap = 0
    INTEGER  :: istage = 1
    INTEGER  :: stage_dap(NUM_PHENO_STAGES) = 0
    LOGICAL  :: stage_reached(NUM_PHENO_STAGES) = .FALSE.
    REAL(wp) :: cum_germ = 0.0_wp
    REAL(wp) :: cum_sub = 0.0_wp
    REAL(wp) :: cum_cole = 0.0_wp
    LOGICAL  :: germ_done = .FALSE.
    LOGICAL  :: sub_done = .FALSE.
    LOGICAL  :: cole_done = .FALSE.
    LOGICAL  :: emerge_done = .FALSE.
    INTEGER  :: dap_germ = 0
    INTEGER  :: dap_sub = 0
    INTEGER  :: dap_cole = 0
    REAL(wp) :: tdu_at_emerge = 0.0_wp  ! cum_tdu the day emergence is declared: start of the VF/DF window
    LOGICAL  :: ts_done = .FALSE.       ! terminal spikelet reached: end of the VF/DF window
  END TYPE PhenologyState

  PUBLIC :: initialize_phenology_params
  PUBLIC :: read_phenology_cultivar
  PUBLIC :: initialize_phenology_state
  PUBLIC :: bootstrap_phenology_state
  PUBLIC :: update_phenology_day
  PUBLIC :: calendar_doy
  PUBLIC :: calc_dev_stage_0_2
  PUBLIC :: calc_dev_stage_rate_0_2
  ! --------------------------------------------------------------------------- !
  ! The optional canopy clock driven by top-soil-layer temperature, which lives
  ! in main_MMCO2Crop and is an alternative formulation tested and not adopted
  ! (Part 1, Table 6), must pass the soil temperature through the SAME
  ! degree-day law and the SAME cultivar cardinal-temperature rows as the
  ! air-temperature clock, and must close its window on the model's own
  ! terminal-spikelet threshold. Both are therefore exported here rather than
  ! re-implemented in the caller.
  ! --------------------------------------------------------------------------- !
  PUBLIC :: calc_dtt
  PUBLIC :: TS_TDU_MCMASTER

CONTAINS

  SUBROUTINE initialize_phenology_params(params)
    TYPE(PhenologyParams), INTENT(OUT) :: params

    ! Stage labels. The thermal thresholds that declare each stage are supplied
    ! through the cultivar parameter file.
    params%stage_name(1) = 'Emergence'
    params%stage_name(2) = '3_leaves'
    params%stage_name(3) = 'Tillering'
    params%stage_name(4) = 'Booting'
    params%stage_name(5) = 'Heading'
    params%stage_name(6) = 'Anthesis'
    params%stage_name(7) = 'Grain_filling'
    params%stage_name(8) = 'Maturity'
  END SUBROUTINE initialize_phenology_params

  SUBROUTINE read_phenology_cultivar(file_name, params, ios_out)
    CHARACTER(LEN=*), INTENT(IN) :: file_name
    TYPE(PhenologyParams), INTENT(INOUT) :: params
    INTEGER, INTENT(OUT) :: ios_out

    INTEGER :: iunit, i, use_beta_int
    CHARACTER(LEN=20) :: dummy_label
    ! ------------------------------------------------------------------------- !
    ! Reads the cultivar parameter file: one value per line, value first and
    ! label second, blank lines and lines beginning with '!' skipped. The
    ! records are positional, so the order of the reads below is the file
    ! format.
    ! ------------------------------------------------------------------------- !
    CALL initialize_phenology_params(params)
    ios_out = 0

    iunit = 21
    OPEN(unit=iunit, file=TRIM(file_name), status='old', iostat=ios_out)
    IF (ios_out /= 0) RETURN

    CALL read_int(iunit, params%gdd_method, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    ! ------------------------------------------------------------------------- !
    ! RETIRED-OPTION REJECTION (v1.0). calc_dtt below implements GDD_Method 1 only
    ! (simple growing-degree day, max(0, T - Tb)); the code of every other method was
    ! removed by the v1.0 clean-up, so any other value is refused here rather than
    ! left to return an undefined thermal-time increment.
    ! ------------------------------------------------------------------------- !
    IF (params%gdd_method /= 1) THEN
      WRITE (*,'(A)') ' ------------------------------------------------------------'
      WRITE (*,'(A,A)') ' REFUSING TO RUN: GDD_Method in ', TRIM(file_name)
      WRITE (*,'(A,I0)') '   read ', params%gdd_method
      WRITE (*,'(A)') '   M&M-CO2Wheat v1.0 implements only GDD_Method = 1 (simple'
      WRITE (*,'(A)') '   growing-degree day, max(0, T_mean - T_base)). The code of every'
      WRITE (*,'(A)') '   other method was removed from the release; it is refused here so'
      WRITE (*,'(A)') '   that calc_dtt can never return an undefined value.'
      WRITE (*,'(A)') '   This option is not available in v1.0 (see Part 1, Table 6).'
      WRITE (*,'(A)') ' ------------------------------------------------------------'
      ERROR STOP 'phenology_cultivar.txt: retired GDD_Method'
    END IF
    CALL read_real(iunit, params%t_base, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%t_opt_lower, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%t_opt_upper, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%t_max, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%photo_sens, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%vern_factor, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900

    DO i = 1, NUM_PHENO_STAGES
      CALL read_real(iunit, params%stage_thresh(i), dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
    END DO

    CALL read_int(iunit, params%sowing_year, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%sowing_month, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%sowing_day, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(1), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_two_leaves_dap, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(2), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(3), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(4), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(5), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(6), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(7), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_int(iunit, params%obs_stage_dap(8), dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_green_g_plant, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_stem_g_plant, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_root_g_plant, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_storage_g_plant, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_leafN_frac, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_stemN_frac, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_rootN_frac, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    CALL read_real(iunit, params%init_storageN_frac, dummy_label, ios_out)
    IF (ios_out /= 0) GOTO 900
    ! ------------------------------------------------------------------------- !
    ! inactive alternative, see Part 1 Table 6
    ! Selector and the eleven parameters of the beta-function emergence model.
    ! ------------------------------------------------------------------------- !
    CALL read_int(iunit, use_beta_int, dummy_label, ios_out)
    IF (ios_out /= 0) THEN
      use_beta_int = 0
      ios_out = 0
    END IF
    ! ------------------------------------------------------------------------- !
    ! Use_Beta_Emerge keeps both of its values in v1.0: 0 is the release setting and
    ! 1 selects the beta-function emergence model of Part 1 Table 6, which is built and
    ! off. No other value has code, so it is refused rather than read as 'not 1'.
    ! ------------------------------------------------------------------------- !
    IF (use_beta_int /= 0 .AND. use_beta_int /= 1) THEN
      WRITE (*,'(A)') ' ------------------------------------------------------------'
      WRITE (*,'(A,A)') ' REFUSING TO RUN: Use_Beta_Emerge in ', TRIM(file_name)
      WRITE (*,'(A,I0)') '   read ', use_beta_int
      WRITE (*,'(A)') '   Implemented values: 0 (observed emergence dates, the v1.0'
      WRITE (*,'(A)') '   setting) and 1 (beta-function germination and coleoptile model,'
      WRITE (*,'(A)') '   Part 1 Table 6, built and off).'
      WRITE (*,'(A)') '   This option is not available in v1.0 (see Part 1, Table 6).'
      WRITE (*,'(A)') ' ------------------------------------------------------------'
      ERROR STOP 'phenology_cultivar.txt: invalid Use_Beta_Emerge'
    END IF
    params%use_beta_emerge = (use_beta_int == 1)

    IF (params%use_beta_emerge) THEN
      CALL read_real(iunit, params%sowing_depth_cm, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%beta_tb, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%beta_tm, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%mu_germ, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%alpha_germ, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%beta_germ, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%mu_cole, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%alpha_cole, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%beta_cole, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%a_sub, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%b_sub, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%alpha_sub, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
      CALL read_real(iunit, params%beta_sub, dummy_label, ios_out)
      IF (ios_out /= 0) GOTO 900
    END IF
    ! ------------------------------------------------------------------------- !
    ! The optional cultivar canopy-architecture rows at the END of the file,
    ! after Use_Beta_Emerge and its own optional block, present or not. Each
    ! read is attempted in order; a row absent at end of file is not an error:
    ! the field is explicitly reassigned its declared default (read_real leaves
    ! its INTENT(OUT) argument undefined after a failed READ, so the default is
    ! written back rather than relied upon to survive) and ios_out is reset to
    ! 0 so the file still closes cleanly below. Reading stops at the first
    ! missing row, a file not being able to carry row N+1 without row N, so
    ! once one row is absent every row after it keeps its default too.
    ! ------------------------------------------------------------------------- !
    CALL read_real(iunit, params%phyllochron_degCd, dummy_label, ios_out)
    IF (ios_out /= 0) THEN
      params%phyllochron_degCd = 73.9_wp
      ios_out = 0
    ELSE
      CALL read_real(iunit, params%blade_l_max_cm, dummy_label, ios_out)
      IF (ios_out /= 0) THEN
        params%blade_l_max_cm = 20.8_wp
        ios_out = 0
      ELSE
        CALL read_real(iunit, params%blade_w_max_cm, dummy_label, ios_out)
        IF (ios_out /= 0) THEN
          params%blade_w_max_cm = 1.45_wp
          ios_out = 0
        ELSE
          ! Two further optional rows, continuing the same sequential chain
          ! immediately after Blade_W_max_cm, with the same end-of-file default
          ! and explicit reset as above.
          CALL read_real(iunit, params%sheath_l_max_cm, dummy_label, ios_out)
          IF (ios_out /= 0) THEN
            params%sheath_l_max_cm = 17.0_wp
            ios_out = 0
          ELSE
            CALL read_real(iunit, params%sheath_d_cm, dummy_label, ios_out)
            IF (ios_out /= 0) THEN
              params%sheath_d_cm = 0.45_wp
              ios_out = 0
            ELSE
              ! Sixth optional row, Axes_max_per_plant (tillering capacity),
              ! same end-of-file default and explicit reset; absent -> 0, no
              ! cap.
              CALL read_real(iunit, params%axes_max_per_plant, dummy_label, ios_out)
              IF (ios_out /= 0) THEN
                params%axes_max_per_plant = 0.0_wp
                ios_out = 0
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF

900 CONTINUE
    CLOSE(iunit)
  END SUBROUTINE read_phenology_cultivar

  SUBROUTINE initialize_phenology_state(state, planting_year, planting_month, planting_day, latitude_deg)
    TYPE(PhenologyState), INTENT(OUT) :: state
    INTEGER, INTENT(IN) :: planting_year, planting_month, planting_day
    REAL(wp), INTENT(IN) :: latitude_deg

    ! Stores planting calendar and site latitude for the daily phenology update.
    state%planting_doy = calendar_doy(planting_year, planting_month, planting_day)
    state%planting_year = planting_year
    state%latitude_deg = latitude_deg
  END SUBROUTINE initialize_phenology_state

  SUBROUTINE bootstrap_phenology_state(params, state, first_year, first_month, first_day)
    TYPE(PhenologyParams), INTENT(IN) :: params
    TYPE(PhenologyState), INTENT(INOUT) :: state
    INTEGER, INTENT(IN) :: first_year, first_month, first_day

    INTEGER :: elapsed_before_day, highest_stage, i
    ! ------------------------------------------------------------------------- !
    ! Aligns the state with the start of the weather record:
    ! - forcing starting before sowing, or on the sowing day, leaves the state
    !   inactive, so the crop starts at sowing;
    ! - forcing starting after sowing initializes the state from the observed
    !   stage dates (days after planting) already passed, so that the crop is
    !   not restarted from sowing on the first available weather day.
    ! ------------------------------------------------------------------------- !
    elapsed_before_day = serial_day(first_year, first_month, first_day) - &
                         serial_day(params%sowing_year, params%sowing_month, params%sowing_day)
    IF (elapsed_before_day <= 0) RETURN

    state%dap = elapsed_before_day
    highest_stage = 0
    DO i = 1, NUM_PHENO_STAGES
      IF (params%obs_stage_dap(i) > 0 .AND. params%obs_stage_dap(i) <= elapsed_before_day) THEN
        highest_stage = i
      END IF
    END DO

    IF (highest_stage <= 0) RETURN

    DO i = 1, highest_stage
      state%stage_dap(i) = params%obs_stage_dap(i)
      state%stage_reached(i) = .TRUE.
    END DO
    state%emerge_done = state%stage_reached(1)
    state%istage = MIN(highest_stage + 1, NUM_PHENO_STAGES + 1)
    state%cum_tdu = params%stage_thresh(highest_stage)
  END SUBROUTINE bootstrap_phenology_state

  SUBROUTINE update_phenology_day(params, state, year_in, month_in, day_in, tmax_air, tmin_air, t_pre_emerge, t_post_emerge)
    TYPE(PhenologyParams), INTENT(IN) :: params
    TYPE(PhenologyState), INTENT(INOUT) :: state
    INTEGER, INTENT(IN) :: year_in, month_in, day_in
    REAL(wp), INTENT(IN) :: tmax_air, tmin_air, t_pre_emerge, t_post_emerge

    REAL(wp) :: t_driver, vd, dtdu, rate_germ, rate_sub, rate_cole
    INTEGER :: doy_in
    ! ------------------------------------------------------------------------- !
    ! Daily phenology update, in order: the day's thermal time DTT, the
    ! vernalization and photoperiod factors and the cumulative thermal
    ! development units (McMaster et al., 2008: Eq. (1), p. 563; Fig. 1,
    ! p. 564; Eqs. (2)-(3), p. 564; Eq. (4), p. 565), then the stage loop that
    ! declares every stage whose threshold the cumulative units have passed.
    ! ------------------------------------------------------------------------- !
    IF (year_in < state%planting_year) RETURN
    doy_in = calendar_doy(year_in, month_in, day_in)
    IF (serial_day(year_in, month_in, day_in) < serial_day(params%sowing_year, params%sowing_month, params%sowing_day)) RETURN
    ! ------------------------------------------------------------------------- !
    ! Days after planting are counted before the terminal-stage RETURN below,
    ! so that dap and every quantity derived from it keep advancing after
    ! thermal maturity is declared.
    ! ------------------------------------------------------------------------- !
    state%dap = serial_day(year_in, month_in, day_in) - serial_day(params%sowing_year, params%sowing_month, params%sowing_day) + 1
    IF (state%istage > NUM_PHENO_STAGES) RETURN

    ! Use separate thermal drivers before and after emergence.
    IF (state%emerge_done) THEN
      t_driver = t_post_emerge
    ELSE
      t_driver = t_pre_emerge
    END IF

    state%current_dtt = calc_dtt(t_driver, params%t_base, params%t_opt_lower, params%t_opt_upper, params%t_max, params%gdd_method)

    IF (params%dev_factor_window == 1) THEN
      ! ----------------------------------------------------------------------- !
      ! Development-factor window 1, the configuration shipped with v1.0.
      ! Vernalization days accumulate from the first phenology day (sowing, the
      ! germination sub-model being optional) until terminal spikelet
      ! (McMaster et al., 2008, p. 564), and are reduced again on a hot day
      ! before ten of them have accumulated (devernalization).
      ! ----------------------------------------------------------------------- !
      IF (.NOT. state%ts_done) THEN
        vd = calc_vern_days(t_driver)
        state%cum_vd = state%cum_vd + vd
        IF (state%cum_vd < 10.0_wp .AND. tmax_air > 30.0_wp) THEN
          state%cum_vd = state%cum_vd - 0.05_wp * (tmax_air - 30.0_wp)
          state%cum_vd = MAX(0.0_wp, state%cum_vd)
        END IF
      END IF
      state%current_vf = calc_vf(state%cum_vd, params%vern_factor)
      state%current_photoperiod_h = calc_daylength(doy_in, state%latitude_deg)
      state%current_df = calc_df(state%current_photoperiod_h, params%photo_sens)
      ! Window: the factors act only after seedling emergence (p. 565) and until
      ! terminal spikelet, 400 degCd TDU after emergence (p. 563); Eq. (1) p. 563
      ! takes the MINIMUM of the two factors, not their product.
      IF (state%emerge_done .AND. .NOT. state%ts_done) THEN
        dtdu = state%current_dtt * MIN(state%current_vf, state%current_df)
      ELSE
        dtdu = state%current_dtt
      END IF
      state%current_dtdu = dtdu
      state%cum_tdu = state%cum_tdu + dtdu
      IF (state%emerge_done .AND. .NOT. state%ts_done) THEN
        IF (state%cum_tdu - state%tdu_at_emerge >= TS_TDU_MCMASTER) state%ts_done = .TRUE.
      END IF
    END IF
    ! ------------------------------------------------------------------------- !
    ! inactive alternative, see Part 1 Table 6
    ! Beta-function germination, subcrown internode and coleoptile model
    ! (Wang et al., 2009, NJAS 57:65-71): emergence is declared when all three
    ! sub-stages are complete, in place of the thermal threshold
    ! stage_thresh(1) used by the stage loop below.
    ! ------------------------------------------------------------------------- !
    IF (params%use_beta_emerge .AND. .NOT. state%emerge_done) THEN
      IF (.NOT. state%germ_done) THEN
        rate_germ = calc_beta_rate(t_pre_emerge, params%beta_tb, params%beta_tm, &
                                   params%mu_germ, params%alpha_germ, params%beta_germ)
        state%cum_germ = state%cum_germ + rate_germ
        IF (state%cum_germ >= 1.0_wp) THEN
          state%germ_done = .TRUE.
          state%dap_germ = state%dap
        END IF
      END IF

      IF (state%germ_done .AND. .NOT. state%sub_done) THEN
        IF (params%sowing_depth_cm > DEPTH_THRESH_CM) THEN
          rate_sub = calc_beta_rate_subcrown(t_pre_emerge, params%beta_tb, params%beta_tm, params%sowing_depth_cm, &
                                             params%a_sub, params%b_sub, params%alpha_sub, params%beta_sub)
          state%cum_sub = state%cum_sub + rate_sub
          IF (state%cum_sub >= 1.0_wp) THEN
            state%sub_done = .TRUE.
            state%dap_sub = state%dap
          END IF
        ELSE
          state%sub_done = .TRUE.
          state%dap_sub = state%dap
        END IF
      END IF

      IF (state%germ_done .AND. state%sub_done .AND. .NOT. state%cole_done) THEN
        rate_cole = calc_beta_rate(t_pre_emerge, params%beta_tb, params%beta_tm, &
                                   params%mu_cole, params%alpha_cole, params%beta_cole)
        state%cum_cole = state%cum_cole + rate_cole
        IF (state%cum_cole >= 1.0_wp) THEN
          state%cole_done = .TRUE.
          state%dap_cole = state%dap
        END IF
      END IF

      IF (state%germ_done .AND. state%sub_done .AND. state%cole_done) THEN
        state%emerge_done = .TRUE.
        state%tdu_at_emerge = state%cum_tdu  ! start of the VF/DF window
        state%stage_dap(1) = state%dap
        state%stage_reached(1) = .TRUE.
        IF (state%istage == 1) state%istage = 2
      END IF
    END IF

    ! Stage loop: declare every stage whose genotype threshold the cumulative
    ! thermal development units have reached. Emergence is skipped here when the
    ! beta-function emergence model above owns it.
    DO WHILE (state%istage <= NUM_PHENO_STAGES)
      IF (state%istage == 1 .AND. params%use_beta_emerge) EXIT
      IF (state%cum_tdu >= params%stage_thresh(state%istage)) THEN
        state%stage_dap(state%istage) = state%dap
        state%stage_reached(state%istage) = .TRUE.
        IF (state%istage == 1) THEN
          state%emerge_done = .TRUE.
          state%tdu_at_emerge = state%cum_tdu  ! start of the VF/DF window
        END IF
        state%istage = state%istage + 1
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE update_phenology_day

  INTEGER FUNCTION calendar_doy(year_in, month_in, day_in)
    INTEGER, INTENT(IN) :: year_in, month_in, day_in
    INTEGER :: mdays(12)

    ! Day of year from a calendar date, leap years included.
    mdays = [31,28,31,30,31,30,31,31,30,31,30,31]
    IF (is_leap_year(year_in)) mdays(2) = 29
    calendar_doy = day_in
    IF (month_in > 1) calendar_doy = calendar_doy + SUM(mdays(1:month_in-1))
  END FUNCTION calendar_doy

  REAL(wp) FUNCTION calc_dev_stage_0_2(params, state)
    TYPE(PhenologyParams), INTENT(IN) :: params
    TYPE(PhenologyState), INTENT(IN) :: state
    REAL(wp) :: th_e, th_a, th_m
    ! ------------------------------------------------------------------------- !
    ! The single development scale the growth, yield and nitrogen equations run
    ! on:
    !   0 = Emergence reached, 1 = Anthesis reached, 2 = Maturity reached,
    !   linear in the CUMULATIVE THERMAL DEVELOPMENT UNITS within each phase
    !   and clamped to [0, 2].
    !
    ! It is formed from the same cum_tdu that drives the stage engine,
    ! cum_tdu = sum(DTT * min(VF, DF)) inside the emergence-to-terminal-spikelet
    ! window and sum(DTT) outside it (McMaster et al., 2008, Annals of Botany
    ! 102:561-569, Eq. (1) p. 563; DTT: Fig. 1 p. 564; VF: Eqs. (2)-(3) p. 564;
    ! DF: Eq. (4) p. 565). Temperature, vernalization and photoperiod therefore
    ! act on growth development exactly as they act on stage development: the
    ! day the engine declares Anthesis (stage_thresh(6)) this function returns
    ! exactly 1, and the day it declares Maturity (stage_thresh(8)) exactly 2.
    ! There is no calendar anchor in it, so it responds to a changed
    ! temperature in the same way as the stage dates do.
    ! ------------------------------------------------------------------------- !
    th_e = params%stage_thresh(1)
    th_a = MAX(params%stage_thresh(6), th_e + 1.0_wp)
    th_m = MAX(params%stage_thresh(8), th_a + 1.0_wp)

    IF (state%cum_tdu <= th_e) THEN
      calc_dev_stage_0_2 = 0.0_wp
    ELSE IF (state%cum_tdu < th_a) THEN
      calc_dev_stage_0_2 = (state%cum_tdu - th_e) / (th_a - th_e)
    ELSE IF (state%cum_tdu < th_m) THEN
      calc_dev_stage_0_2 = 1.0_wp + (state%cum_tdu - th_a) / (th_m - th_a)
    ELSE
      calc_dev_stage_0_2 = 2.0_wp
    END IF
  END FUNCTION calc_dev_stage_0_2

  REAL(wp) FUNCTION calc_dev_stage_rate_0_2(params, state)
    TYPE(PhenologyParams), INTENT(IN) :: params
    TYPE(PhenologyState), INTENT(IN) :: state
    REAL(wp) :: th_e, th_a, th_m
    ! ------------------------------------------------------------------------- !
    ! Today's increment of the same development scale: current_dtdu, the day's
    ! thermal development units stored by update_phenology_day, divided by the
    ! span of the phase the crop is in. It is read by the senescence ageing
    ! law. It is zero before emergence and after thermal maturity, development
    ! then being complete.
    ! ------------------------------------------------------------------------- !
    th_e = params%stage_thresh(1)
    th_a = MAX(params%stage_thresh(6), th_e + 1.0_wp)
    th_m = MAX(params%stage_thresh(8), th_a + 1.0_wp)

    IF (state%cum_tdu <= th_e) THEN
      calc_dev_stage_rate_0_2 = 0.0_wp
    ELSE IF (state%cum_tdu < th_a) THEN
      calc_dev_stage_rate_0_2 = state%current_dtdu / (th_a - th_e)
    ELSE IF (state%cum_tdu < th_m) THEN
      calc_dev_stage_rate_0_2 = state%current_dtdu / (th_m - th_a)
    ELSE
      calc_dev_stage_rate_0_2 = 0.0_wp
    END IF
  END FUNCTION calc_dev_stage_rate_0_2

  SUBROUTINE read_real(iunit_in, val, label, ios_out)
    INTEGER, INTENT(IN) :: iunit_in
    REAL(wp), INTENT(OUT) :: val
    CHARACTER(LEN=*), INTENT(OUT) :: label
    INTEGER, INTENT(OUT) :: ios_out
    CHARACTER(LEN=256) :: buf

    ios_out = 0
    DO
      READ(iunit_in, '(A)', iostat=ios_out) buf
      IF (ios_out /= 0) RETURN
      buf = ADJUSTL(buf)
      IF (LEN_TRIM(buf) == 0) CYCLE
      IF (buf(1:1) == '!') CYCLE
      READ(buf, *, iostat=ios_out) val, label
      RETURN
    END DO
  END SUBROUTINE read_real

  SUBROUTINE read_int(iunit_in, val, label, ios_out)
    INTEGER, INTENT(IN) :: iunit_in
    INTEGER, INTENT(OUT) :: val
    CHARACTER(LEN=*), INTENT(OUT) :: label
    INTEGER, INTENT(OUT) :: ios_out
    CHARACTER(LEN=256) :: buf

    ios_out = 0
    DO
      READ(iunit_in, '(A)', iostat=ios_out) buf
      IF (ios_out /= 0) RETURN
      buf = ADJUSTL(buf)
      IF (LEN_TRIM(buf) == 0) CYCLE
      IF (buf(1:1) == '!') CYCLE
      READ(buf, *, iostat=ios_out) val, label
      RETURN
    END DO
  END SUBROUTINE read_int

  LOGICAL FUNCTION is_leap_year(year_in)
    INTEGER, INTENT(IN) :: year_in

    is_leap_year = (MOD(year_in,400) == 0) .OR. ((MOD(year_in,4) == 0) .AND. (MOD(year_in,100) /= 0))
  END FUNCTION is_leap_year

  INTEGER FUNCTION serial_day(year_in, month_in, day_in)
    INTEGER, INTENT(IN) :: year_in, month_in, day_in
    INTEGER :: y

    y = year_in - 1
    serial_day = 365 * y + y / 4 - y / 100 + y / 400 + calendar_doy(year_in, month_in, day_in)
  END FUNCTION serial_day

  REAL(wp) FUNCTION calc_dtt(tavg, tbase, toptl, toptu, tmax, method)
    REAL(wp), INTENT(IN) :: tavg, tbase, toptl, toptu, tmax
    INTEGER, INTENT(IN) :: method
    ! ------------------------------------------------------------------------- !
    ! Daily thermal time. Method 1, the simple growing-degree-day form
    ! max(0, T - Tb) with the genotype base temperature, is the only method
    ! implemented: the code of every other method was removed by the v1.0
    ! clean-up. read_phenology_cultivar refuses any other GDD_Method, so the
    ! ELSE below is unreachable in a run that started; it exists so that this
    ! function can never return an undefined value.
    ! ------------------------------------------------------------------------- !
    IF (method == 1) THEN
      calc_dtt = MAX(0.0_wp, tavg - tbase)
    ELSE
      calc_dtt = 0.0_wp
      WRITE (*,'(A,I0)') ' REFUSING TO RUN: calc_dtt reached with GDD_Method = ', method
      WRITE (*,'(A)') '   Only GDD_Method = 1 is implemented in M&M-CO2Wheat v1.0.'
      ERROR STOP 'calc_dtt: retired GDD_Method'
    END IF
  END FUNCTION calc_dtt

  REAL(wp) FUNCTION calc_vern_days(tavg)
    REAL(wp), INTENT(IN) :: tavg

    ! Daily vernalization increment (vernalization days) from the temperature
    ! response of McMaster et al. (2008), Fig. 2, p. 565: full effect from 0 to
    ! 7 deg C, falling linearly to zero at -5 and at 15 deg C.
    IF (tavg < -5.0_wp) THEN
      calc_vern_days = 0.0_wp
    ELSE IF (tavg < 0.0_wp) THEN
      calc_vern_days = (tavg + 5.0_wp) / 5.0_wp
    ELSE IF (tavg <= 7.0_wp) THEN
      calc_vern_days = 1.0_wp
    ELSE IF (tavg <= 15.0_wp) THEN
      calc_vern_days = 1.0_wp - (tavg - 7.0_wp) / 8.0_wp
    ELSE
      calc_vern_days = 0.0_wp
    END IF
  END FUNCTION calc_vern_days

  REAL(wp) FUNCTION calc_vf(cumvd_in, p1v)
    REAL(wp), INTENT(IN) :: cumvd_in, p1v
    REAL(wp) :: vf0

    ! Vernalization factor VF from McMaster et al. (2008), Eqs. (2)-(3),
    ! p. 564, clamped to [0, 1]. A genotype with zero vernalization
    ! requirement (p1v <= 0) has VF = 1 throughout, as a spring wheat does.
    IF (p1v <= 0.0_wp) THEN
      calc_vf = 1.0_wp
      RETURN
    END IF
    vf0 = 1.0_wp - (p1v / 50.0_wp)
    calc_vf = vf0 + cumvd_in / p1v
    calc_vf = MAX(0.0_wp, MIN(1.0_wp, calc_vf))
  END FUNCTION calc_vf

  REAL(wp) FUNCTION calc_daylength(doy_in, lat)
    INTEGER, INTENT(IN) :: doy_in
    REAL(wp), INTENT(IN) :: lat
    REAL(wp) :: lat_rad, dec, depression, cos_ha, ha

    ! Daylength (h) including civil twilight, a solar depression of 6 degrees:
    ! the photoperiod the daylength factor DF below is defined on.
    lat_rad = lat * PI_PHENO / 180.0_wp
    dec = 0.4093_wp * SIN(2.0_wp * PI_PHENO * REAL(doy_in - 82, wp) / 365.0_wp)
    depression = 6.0_wp * PI_PHENO / 180.0_wp
    cos_ha = (-SIN(depression) - SIN(lat_rad) * SIN(dec)) / (COS(lat_rad) * COS(dec))
    cos_ha = MAX(-1.0_wp, MIN(1.0_wp, cos_ha))
    ha = ACOS(cos_ha)
    calc_daylength = 2.0_wp * ha * 12.0_wp / PI_PHENO
  END FUNCTION calc_daylength

  REAL(wp) FUNCTION calc_df(pp_in, p1d)
    REAL(wp), INTENT(IN) :: pp_in, p1d

    ! Daylength factor DF from McMaster et al. (2008), Eq. (4), p. 565, clamped
    ! to [0, 1]. A photoperiod-insensitive genotype (p1d = 0) has DF = 1
    ! throughout.
    calc_df = 1.0_wp - (p1d / 10000.0_wp) * (20.0_wp - pp_in)**2
    calc_df = MAX(0.0_wp, MIN(1.0_wp, calc_df))
  END FUNCTION calc_df

  REAL(wp) FUNCTION calc_beta_rate(t, tb, tm, mu, alpha_p, beta_p)
    REAL(wp), INTENT(IN) :: t, tb, tm, mu, alpha_p, beta_p

    ! inactive alternative, see Part 1 Table 6
    ! Wang et al. (2009), Eq. (7) with the coefficients of their Table 2: the
    ! beta-form temperature response of the daily rate, used for germination
    ! and for coleoptile growth.
    IF (t <= tb .OR. t >= tm) THEN
      calc_beta_rate = 0.0_wp
    ELSE
      calc_beta_rate = EXP(mu) * (t - tb)**alpha_p * (tm - t)**beta_p
    END IF
  END FUNCTION calc_beta_rate

  REAL(wp) FUNCTION calc_beta_rate_subcrown(t, tb, tm, sdepth, a_p, b_p, alpha_p, beta_p)
    REAL(wp), INTENT(IN) :: t, tb, tm, sdepth, a_p, b_p, alpha_p, beta_p
    REAL(wp) :: scale

    ! inactive alternative, see Part 1 Table 6
    ! Wang et al. (2009), Eq. (8) with the coefficients of their Table 2: the
    ! subcrown internode elongation rate, the same beta-form temperature
    ! response scaled by the sowing depth below the 2.5 cm reference.
    IF (t <= tb .OR. t >= tm) THEN
      calc_beta_rate_subcrown = 0.0_wp
    ELSE
      scale = a_p + b_p * (sdepth - 2.5_wp)
      scale = MAX(0.0_wp, scale)
      calc_beta_rate_subcrown = scale * (t - tb)**alpha_p * (tm - t)**beta_p
    END IF
  END FUNCTION calc_beta_rate_subcrown

END MODULE PHENOLOGY_MODULE