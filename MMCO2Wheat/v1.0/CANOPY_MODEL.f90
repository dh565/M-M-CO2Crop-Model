!===============================================================================!
! M&M-CO2Wheat v1.0  |  CANOPY_MODEL.f90
!
! Purpose: The canopy integration of leaf gas exchange and the canopy-temperature
!          solve. Time step HOURLY, with a sub-hourly iteration: 24 hourly steps
!          per day, each a bracketed root-find for the canopy temperature wrapped
!          around a nested fixed point that makes assimilation, stomatal
!          conductance and the leaf-surface vapor pressure deficit mutually
!          consistent, and each returning the hourly demand profile the
!          soil-water repartition later replays. Spatial unit: the canopy as two
!          leaf classes, sunlit and shaded, with depth entering only through the
!          vertical leaf-nitrogen profile; the day's means and totals are
!          returned to the caller.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision); SHARED (the hourly forcing, the canopy, soil
!       and leaf-nitrogen state, and every process switch); energy_balance_module
!       (compute_eb_resistances, solve_eb_3node_inner, esat_pa, h_soil_surface,
!       raa_heat_stability); and the external routines assim (leaf gas exchange,
!       through the explicit interface below), meteo (hourly forcing), solar
!       (beam extinction, sunlit and shaded leaf area, absorbed PAR per class)
!       and th (daylength and solar declination).
! Output: the day's canopy quantities, returned through the argument list -
!       transpiration and soil-evaporation demand, gross and net assimilation,
!       nitrate reduced, canopy, air and soil-surface temperature, vapor pressure
!       deficit, net radiation, soil heat flux, the latent and sensible heat
!       fluxes of canopy and soil, stomatal conductance, intercellular CO2 and
!       the five resistances, each also as a daylight mean. They are consumed by
!       main_MMCO2Crop.f90, which drives the soil-water and soil-heat solvers,
!       growth_module and phenology with them. Written into SHARED: the hourly
!       demand profiles the post-Richards repartition replays (ebh_*), the
!       buffered hourly energy-balance and radiation rows, the solved three-node
!       state, the sunlit and shaded leaf temperatures and the solver counters.
! Called by: main_MMCO2Crop.f90 - compute_daily_canopy_fluxes once per pass of
!       the daily coupling loop; canopy_n0_top is also called there, for the
!       canopy-nitrogen profile output and the lamina-nitrogen ceiling.

!===============================================================================!
MODULE CANOPY_MODEL
  USE TYPES, ONLY: wp
  USE SHARED
  USE energy_balance_module, ONLY: compute_eb_resistances, solve_eb_3node_inner, &
                                   esat_pa, h_soil_surface, raa_heat_stability
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE assim (Rd, A_assim, An_assim, gsw_1, Ca, Ci, N_total_assim, Ds, Narea_in)
      USE TYPES, ONLY: wp
      REAL(wp), INTENT(OUT)   :: Rd, A_assim, An_assim, gsw_1, Ci, N_total_assim
      REAL(wp), INTENT(IN)    :: Ca, Ds
      REAL(wp), INTENT(IN), OPTIONAL :: Narea_in
    END SUBROUTINE assim
  END INTERFACE

CONTAINS

  ! ------------------------------------------------------------------ !
  ! CANOPY LEAF-NITROGEN VERTICAL PROFILE. One definition, used by BOTH call
  ! sites - the flag-leaf diagnostic and the sunlit/shaded class means - so the
  ! two cannot drift apart. Two driver records select the form, and n(l) is the
  ! leaf nitrogen per unit area at the cumulative lamina area index l counted
  ! from the canopy top:
  !   Record 76 = 0: n(l) = n0 exp(-kN l), n0 = N_leaf kN / (1 - exp(-kN LAI)).
  !        The whole lamina-nitrogen pool is distributed exponentially with no
  !        base term inside the profile, the structural minimum being subtracted
  !        uniformly afterwards in assim.f90. The method is the single-term
  !        exponential with no base term inside the profile, as published by
  !        Yin, Struik & Goudriaan 2021, Field Crops Res 271:108254,
  !        Appendix A Eq. (A3), "ni = kn N e^-kn Li /(1 - e^-kn L)".
  !        PROVENANCE: that equation is the canopy-nitrogen profile of the
  !        GECROS model; it is a model formulation, not a measurement.
  !   Record 76 = 1: n(l) = nb_eff + (n0 - nb_eff) exp(-kN l), the measured wheat
  !        form of Bertheloot, Martre & Andrieu 2008, Plant Physiol
  !        148:1707-1720, DOI 10.1104/pp.108.124156, Materials and Methods Eq. 2,
  !        with the structural base nb = 0.40 g N m-2 from the same section and
  !        nb_eff = MIN(nb, N_leaf/LAI) as the degenerate-case guard. n0 follows
  !        from CONSERVATION of N_leaf over 0..LAI:
  !        nb_eff LAI + (n0 - nb_eff)(1 - exp(-kN LAI))/kN = N_leaf.
  !   Record 94, modes 1 to 3: the profile is UNIFORM over the green area of the
  !        flag-leaf rank and declines exponentially below it. The zone depth Fc
  !        is cls_gai_flag, the model's own flag-rank green lamina area index, an
  !        emergent state of the organ-level canopy and neither a constant nor a
  !        fitted depth; it is clamped into [0, LAI], and Fc = 0 returns the
  !        profile above exactly. Mode 2 forces nb_eff = 0, which separates the
  !        zone from the base term; the modes above 3 leave the profile alone and
  !        act on the absorbed-radiation profile instead.
  ! Both forms return a TOTAL Narea, structural plus photosynthetic; assim.f90
  ! subtracts the structural base once. The extinction coefficient kN = 0.41 is
  ! the product of the nitrogen-to-light extinction ratio 0.96 and the light
  ! extinction coefficient 0.425 reported for wheat by Moreau et al. 2012,
  ! Plant Physiol 160:1479-1490, DOI 10.1104/pp.112.199935. Full argument and
  ! every source at items N-1..N-6 and E-1..E-7 of the revision record (internal; not distributed).
  ! ------------------------------------------------------------------ !
  PURE SUBROUTINE canopy_n_zone(N_leaf_in, lai_in, Fc_out, nb_out)
    REAL(wp), INTENT(IN)  :: N_leaf_in, lai_in
    REAL(wp), INTENT(OUT) :: Fc_out, nb_out
    Fc_out = 0.0_wp
    IF (canopy_n_flag_zone >= 1 .AND. canopy_n_flag_zone <= 3) &
      Fc_out = MIN(MAX(0.0_wp, cls_gai_flag), lai_in)
    nb_out = 0.0_wp
    IF (canopy_n_flag_zone /= 2) &
      nb_out = MIN(nb_canopy, N_leaf_in / lai_in)
  END SUBROUTINE canopy_n_zone

  PURE FUNCTION canopy_n0_top(N_leaf_in, lai_in, kN_in) RESULT(n0)
    REAL(wp), INTENT(IN) :: N_leaf_in, lai_in, kN_in
    REAL(wp) :: n0, denom, nb_eff, Fc, Lb, Eint
    IF (canopy_n_flag_zone >= 1 .AND. canopy_n_flag_zone <= 3) THEN
      ! Flat over 0..Fc, exponential below it; item E-4 of the revision record (internal; not distributed).
      CALL canopy_n_zone(N_leaf_in, lai_in, Fc, nb_eff)
      Lb   = MAX(0.0_wp, lai_in - Fc)
      Eint = (1.0_wp - EXP(-kN_in * Lb)) / kN_in
      n0   = (N_leaf_in - nb_eff * (Lb - Eint)) / (Fc + Eint)
    ELSE
      denom = 1.0_wp - EXP(-kN_in * lai_in)
      nb_eff = MIN(nb_canopy, N_leaf_in / lai_in)
      n0 = nb_eff + (N_leaf_in - nb_eff * lai_in) * kN_in / denom
    END IF
  END FUNCTION canopy_n0_top
  ! ------------------------------------------------------------------ !
  ! Cumulative SUNLIT leaf nitrogen per unit ground area: the integral 
  ! over 0..LAI of n(l) exp(-kb l) dl on the same profile as 
  ! canopy_n0_top, with kb the beam extinction coefficient. 
  ! The first term is nb_eff x sunlit_lai exactly, since 
  ! solar.f90 defines sunlit_lai = (1 - exp(-kb LAI))/kb with kb = Kdr.
  ! ------------------------------------------------------------------ !
  PURE FUNCTION canopy_n_sunlit(N_leaf_in, lai_in, kN_in, kb_in) RESULT(Nsun)
    REAL(wp), INTENT(IN) :: N_leaf_in, lai_in, kN_in, kb_in
    REAL(wp) :: Nsun, n0, nb_eff, Fc, Lb
    n0 = canopy_n0_top(N_leaf_in, lai_in, kN_in)
    IF (canopy_n_flag_zone >= 1 .AND. canopy_n_flag_zone <= 3) THEN
      ! -------------------------------------------------------------- !
      ! The same profile under the beam weight; item E-4 of the revision 
      ! record (internal; not distributed).
      ! -------------------------------------------------------------- !
      CALL canopy_n_zone(N_leaf_in, lai_in, Fc, nb_eff)
      Lb   = MAX(0.0_wp, lai_in - Fc)
      Nsun = n0     * (1.0_wp - EXP(-kb_in * Fc)) / kb_in                &
           + nb_eff * (EXP(-kb_in * Fc) - EXP(-kb_in * lai_in)) / kb_in  &
           + (n0 - nb_eff) * EXP(-kb_in * Fc)                            &
           * (1.0_wp - EXP(-(kN_in + kb_in) * Lb)) / (kN_in + kb_in)
    ELSE
      nb_eff = MIN(nb_canopy, N_leaf_in / lai_in)
      Nsun = nb_eff * (1.0_wp - EXP(-kb_in * lai_in)) / kb_in &
           + (n0 - nb_eff) * (1.0_wp - EXP(-(kN_in + kb_in) * lai_in)) / (kN_in + kb_in)
    END IF
  END FUNCTION canopy_n_sunlit

  SUBROUTINE compute_daily_canopy_fluxes(day_idx, unit_diurnal, do_write, &
                                         t_i_day, e_i_day, &
                                         An_day_out, A_day_out, N_day_out, &
                                         T_leaf_day_out, T_air_day_out, T_soil_day_out, &
                                         VPD_day_out, Rn_day_out, G_day_out, &
                                         LE_day_out, H_day_out, LEc_day_out, &
                                         LEs_day_out, Hc_day_out, Hs_day_out, &
                                         gsw_mm_day_out, Ci_day_out, &
                                         Ac_day_out, As_day_out, rca_day_out, &
                                         rcs_day_out, rsa_day_out, &
                                         rss_day_out, gsw_day_out, &
                                         rca_daylight_out, rcs_daylight_out, &
                                         rsa_daylight_out, rss_daylight_out, &
                                         tr_daylight_out, es_daylight_out)
    INTEGER, INTENT(IN)   :: day_idx
    INTEGER, INTENT(IN)   :: unit_diurnal
    LOGICAL, INTENT(IN)   :: do_write
    REAL(wp), INTENT(OUT) :: t_i_day, e_i_day
    REAL(wp), INTENT(OUT) :: An_day_out, A_day_out, N_day_out
    REAL(wp), INTENT(OUT) :: T_leaf_day_out, T_air_day_out, T_soil_day_out, VPD_day_out
    REAL(wp), INTENT(OUT) :: Rn_day_out, G_day_out, LE_day_out, H_day_out
    REAL(wp), INTENT(OUT) :: LEc_day_out, LEs_day_out, Hc_day_out, Hs_day_out
    REAL(wp), INTENT(OUT) :: gsw_mm_day_out, Ci_day_out
    REAL(wp), INTENT(OUT) :: Ac_day_out, As_day_out
    REAL(wp), INTENT(OUT) :: rca_day_out, rcs_day_out, rsa_day_out, rss_day_out
    REAL(wp), INTENT(OUT) :: gsw_day_out, rca_daylight_out, rcs_daylight_out
    REAL(wp), INTENT(OUT) :: rsa_daylight_out, rss_daylight_out, tr_daylight_out, es_daylight_out

    ! Local variables
    INTEGER :: J, iter_eb, nsteps
    REAL(wp) :: weight, hour_val, temp_sum
    REAL(wp) :: T2, f2
    REAL(wp) :: gsw_sun, gsw_shade, lai_sum
    REAL(wp) :: Ds_s, Ds_new, Ci_step
    REAL(wp) :: EB_err, H_total_s
    REAL(wp) :: L_Temp_h, mass_conv
    REAL(wp) :: VpdMcf_s
    REAL(wp) :: par_save, n0_flag, Rd_f, A_f, An_f, gsw_f, Ci_f, N_f
    REAL(wp) :: day_weight_sum
    REAL(wp) :: gsw_day_sum, rca_day_sum, rcs_day_sum, rsa_day_sum, rss_day_sum
    REAL(wp) :: t_i_daylight, e_i_daylight
    REAL(wp) :: ci_day_sum          ! daylight-masked Ci accumulator
    LOGICAL  :: is_day

    ! Initializations
    LEc_int_pot_day  = zero   ! daily wet-canopy (interception) evaporation demand (mm/d)
    t_i_day          = zero; e_i_day = zero
    An_day_out       = zero; A_day_out = zero; N_day_out = zero
    T_leaf_day_out   = zero; T_air_day_out = zero; T_soil_day_out = zero; VPD_day_out = zero
    Rn_day_out       = zero; G_day_out = zero; LE_day_out = zero; H_day_out = zero
    LEc_day_out      = zero; LEs_day_out = zero; Hc_day_out = zero; Hs_day_out = zero
    gsw_mm_day_out   = zero; Ci_day_out = zero
    Ac_day_out       = zero; As_day_out = zero
    rca_day_out      = zero; rcs_day_out = zero; rsa_day_out = zero; rss_day_out = zero
    gsw_day_out      = zero; rca_daylight_out = zero; rcs_daylight_out = zero
    rsa_daylight_out = zero; rss_daylight_out = zero
    day_weight_sum   = zero; gsw_day_sum = zero
    rca_day_sum      = zero; rcs_day_sum = zero; rsa_day_sum = zero; rss_day_sum = zero
    t_i_daylight     = zero; e_i_daylight = zero
    ci_day_sum       = zero; ci_daylight_mean_diag = zero

    CALL th (dTsr, Daylength, Solar_Decl)
    nsteps = 24

    temp_offset = 0.0_wp
    temp_sum = 0.0_wp
    DO J = 1, 24
      hour_val = REAL(J - 1, wp)
      th_s = hour_val
      CALL meteo (SolarElev_s, dfrad_s, drrad_s, VpdRef_s, HourTemp_s, HourNetRad_s)
      temp_sum = temp_sum + HourTemp_s
    END DO
    temp_offset = temp_s - (temp_sum / 24.0_wp)

    DO J = 1, nsteps
      hour_val = REAL(J - 1, wp)
      weight = 1.0_wp / 24.0_wp
      th_s = hour_val

      w_eb_s = w_d_s

      CALL meteo (SolarElev_s, dfrad_s, drrad_s, VpdRef_s, HourTemp_s, HourNetRad_s)
      Rn_s = HourNetRad_s
      CALL solar (Kdr_s, I_drrad_h_s, I_dfrad_h_s, Absor_IPAR_sunlit_h_s, Absor_IPAR_shaded_h_s, sunlit_lai_s, shaded_lai_s)
      ! Effective LAI for canopy resistance scaling (no max cap).
      dLe_s = gParam_lai

      Ds_s = VpdRef_s
      ! ----------------------------------------------------------------- !
      ! Substep-consistent soil-heat state for the demand solve. 
      ! From the second pass of the day's coupling loop onward, the soil 
      ! heat flux is conducted against the COMMITTED node-1 trajectory of 
      ! the previous pass, captured in the repartition, instead of the 
      ! end-of-previous-day frozen profile. The first pass of each day 
      ! keeps the frozen-profile seed; the outer fixed point then converges 
      ! the demand and the committed soil heat flux onto one trajectory.
      ! ----------------------------------------------------------------- !
      IF (energy_balance_mode == 2 .AND. eba_node1_valid) THEN
        shf_node1_use = .TRUE.
        shf_node1_val = eba_node1(J)
        node_armed_row(J) = shf_node1_val   ! armed-node identity record
        ! --------------------------------------------------------------- !
        ! Supply cap on the demand-state soil latent heat: the previous pass's
        ! committed value for this substep, on the same trajectory as node 1.
        ! --------------------------------------------------------------- !
        eb_les_cap_use = .TRUE.
        eb_les_cap_val = eba_LEs(J)
        ! --------------------------------------------------------------- !
        ! Supply cap on the demand-state canopy latent heat: the previous 
        ! pass's committed value, delivered transpiration plus interception.
        ! --------------------------------------------------------------- !
        eb_lec_cap_use = .TRUE.
        eb_lec_cap_val = eba_LEc(J)
      ELSE
        shf_node1_use = .FALSE.
        eb_les_cap_use = .FALSE.
        eb_lec_cap_use = .FALSE.
        node_armed_row(J) = -999.0_wp       ! unarmed (first pass)
      END IF

      IF (energy_balance_mode >= 1 .AND. (gParam_lai <= 0.0_wp .OR. gParam_hgt_s <= 0.0_wp)) THEN
        ! --------------------------------------------------------------- !
        ! No canopy, before emergence or on bare soil: the canopy 
        ! temperature is undefined and the canopy residual is about zero 
        ! for any value of it, so a bracketed search would lock onto a 
        ! meaningless endpoint. Evaluate once at air temperature; 
        ! the inner solve still returns the soil and canopy-air nodes.
        ! --------------------------------------------------------------- !
        T2 = HourTemp_s; iter_eb = 0
        CALL solve_biophysics(T2, Ds_s, f2, H_total_s, LE_s, LEc_s, LEs_s, gsw_An, Ds_new, Ci_step, VpdMcf_s)
        Ds_s = Ds_new
      ELSE IF (energy_balance_mode >= 1) THEN
        ! --------------------------------------------------------------- !
        ! Bracketed root-find for the canopy temperature. The canopy 
        ! residual is monotone decreasing in it, sensible and latent heat 
        ! both rising with canopy temperature, so bisection cannot miss 
        ! the root. The bracket [Tair - 15, Tair + 70] is deliberately 
        ! WIDE because a dry, sparsely coupled canopy, with a large 
        ! above-canopy resistance, can sit far above air temperature. 
        ! A final evaluation at the accepted temperature guarantees that 
        ! the state left in SHARED is the reported one.
        ! --------------------------------------------------------------- !
        BLOCK
          REAL(wp) :: Tlo, Thi, flo, fhi, Tm, fm
          INTEGER  :: kit
          Tlo = HourTemp_s - 15.0_wp
          Thi = HourTemp_s + 70.0_wp
          CALL solve_biophysics(Tlo, Ds_s, flo, H_total_s, LE_s, LEc_s, LEs_s, gsw_An, Ds_new, Ci_step, VpdMcf_s)
          Ds_s = Ds_new
          CALL solve_biophysics(Thi, Ds_s, fhi, H_total_s, LE_s, LEc_s, LEs_s, gsw_An, Ds_new, Ci_step, VpdMcf_s)
          Ds_s = Ds_new
          T2 = HourTemp_s; f2 = flo; iter_eb = 0
          IF ((.NOT. (ABS(flo) < 1.0e30_wp)) .OR. (.NOT. (ABS(fhi) < 1.0e30_wp)) .OR. &
              flo * fhi > 0.0_wp) THEN
            ! ----------------------------------------------------------- !
            ! No sign change in the bracket, or a non-finite residual: the hour
            ! is evaluated at air temperature instead of silently accepting an
            ! endpoint, and is counted per event and per day for the run log. No
            ! output is clamped.
            ! ----------------------------------------------------------- !
            T2 = HourTemp_s
            eb_noroot_events = eb_noroot_events + 1
            IF (day_idx /= eb_noroot_lastday) THEN
              eb_noroot_days = eb_noroot_days + 1
              eb_noroot_lastday = day_idx
              WRITE (*,'(A,I0)') ' EB NOROOT day=', day_idx
            END IF
          ELSE
            DO kit = 1, 60
              iter_eb = kit
              Tm = 0.5_wp * (Tlo + Thi)
              CALL solve_biophysics(Tm, Ds_s, fm, H_total_s, LE_s, LEc_s, LEs_s, gsw_An, Ds_new, Ci_step, VpdMcf_s)
              Ds_s = Ds_new
              T2 = Tm; f2 = fm
              IF (.NOT. (ABS(fm) < 1.0e30_wp)) THEN
                ! ------------------------------------------------------- !
                ! Non-finite residual mid-bisection: flag and fall back 
                ! to air temperature.
                ! ------------------------------------------------------- !
                T2 = HourTemp_s
                eb_noroot_events = eb_noroot_events + 1
                IF (day_idx /= eb_noroot_lastday) THEN
                  eb_noroot_days = eb_noroot_days + 1
                  eb_noroot_lastday = day_idx
                  WRITE (*,'(A,I0)') ' EB NOROOT day=', day_idx
                END IF
                EXIT
              END IF
              IF (ABS(fm) < 0.2_wp .OR. (Thi - Tlo) < 1.0e-3_wp) EXIT
              IF (flo * fm < 0.0_wp) THEN; Thi = Tm; fhi = fm; ELSE; Tlo = Tm; flo = fm; END IF
            END DO
          END IF
          ! ------------------------------------------------------------- !
          ! Final evaluation at the accepted canopy temperature, so that 
          ! all the state left in SHARED is consistent with it.
          ! ------------------------------------------------------------- !
          CALL solve_biophysics(T2, Ds_s, f2, H_total_s, LE_s, LEc_s, LEs_s, gsw_An, Ds_new, Ci_step, VpdMcf_s)
          Ds_s = Ds_new
        END BLOCK
        ! --------------------------------------------------------------- !
        ! Plausibility gate on the component available energy: warn and 
        ! count, do NOT clamp, so the solver never silently masks a 
        ! pathological state. The count is summarized once at the end 
        ! of the run.
        ! --------------------------------------------------------------- !
        IF ((Rnc_s + Rns_s) > 0.0_wp .AND. &
            (T2 < HourTemp_s - 15.0_wp .OR. T2 > HourTemp_s + 25.0_wp)) THEN
          eb3_implausible_count = eb3_implausible_count + 1
        END IF
      END IF

      T_leaf_s = T2; EB_err = f2
      ! ----------------------------------------------------------------- !
      ! Flux-side well-posedness guard: the temperature root can exist 
      ! while the latent heat demand that goes with it is impossible. 
      ! Detect that and revert THIS HOUR's demand to equilibrium evaporation, 
      ! partitioned between canopy and soil by the canopy cover fraction; 
      ! count it, and never modify a sane hour.
      ! ----------------------------------------------------------------- !
      IF (energy_balance_mode >= 1) THEN
        IF ((.NOT. (ABS(LEc_s) < 1.0e30_wp)) .OR. (.NOT. (ABS(LEs_s) < 1.0e30_wp)) .OR. &
            LEc_s < -50.0_wp .OR. LEs_s < -50.0_wp .OR. &
            (LEc_s + LEs_s) > MAX(0.0_wp, Rnc_s + Rns_s) + 600.0_wp) THEN
          BLOCK
            REAL(wp) :: delta_svp, le_eq, f_canfrac, tt
            tt = HourTemp_s
            delta_svp = 4098.0_wp * (0.6108_wp * EXP(17.27_wp * tt / (tt + 237.3_wp))) &
                        / ((tt + 237.3_wp) ** 2)                    ! kPa/C (FAO-56 Eq.13)
            le_eq = MAX(0.0_wp, (delta_svp / (delta_svp + 0.066_wp)) * &
                    MAX(0.0_wp, Rnc_s + Rns_s))
            f_canfrac = 1.0_wp - EXP(-0.5_wp * MAX(0.0_wp, gParam_lai))
            LEc_s = f_canfrac * le_eq
            LEs_s = (1.0_wp - f_canfrac) * le_eq
            eb_fluxfb_events = eb_fluxfb_events + 1
            IF (day_idx /= eb_fluxfb_lastday) THEN
              eb_fluxfb_days = eb_fluxfb_days + 1
              eb_fluxfb_lastday = day_idx
              WRITE (*,'(A,I0)') ' EB FLUXFB day=', day_idx
            END IF
          END BLOCK
        END IF
      END IF
      L_Temp_h = (2500.8_wp - 2.36_wp * T_leaf_s + 0.0016_wp * T_leaf_s**2) * 1000.0_wp
      mass_conv = 86400.0_wp / 1000.0_wp

      t_i_day = t_i_day + (LEc_s / L_Temp_h) * mass_conv * weight
      e_i_day = e_i_day + (LEs_s / L_Temp_h) * mass_conv * weight
      ! ----------------------------------------------------------------- !
      ! Per-substep soil-evaporation demand (mm in this window), by the same
      ! conversion as e_i_day, so the series sums to it exactly. Refreshed on
      ! every pass of the coupling loop.
      ! ----------------------------------------------------------------- !
      IF (J <= EB_MAXSTEP) ebh_esoil_dem_mm(J) = (LEs_s / L_Temp_h) * mass_conv * weight * 1000.0_wp
      ! ----------------------------------------------------------------- !
      ! Per-substep transpiration demand, by the same conversion as t_i_day.
      ! ----------------------------------------------------------------- !
      IF (J <= EB_MAXSTEP) ebh_etr_dem_mm(J) = (LEc_s / L_Temp_h) * mass_conv * weight * 1000.0_wp
      ! ----------------------------------------------------------------- !
      An_day_out = An_day_out + (Assim_net_total * 44.0_wp * 1.0E-6_wp) * weight * 86400.0_wp
      A_day_out  = A_day_out  + (Assim_total * 44.0_wp * 1.0E-6_wp) * weight * 86400.0_wp
      ! ----------------------------------------------------------------- !
      ! N_total is returned from assim() as umol N m-2 s-1.
      ! Convert to g N m-2 d-1 using 14 g mol-1 nitrogen.
      ! ----------------------------------------------------------------- !
      N_day_out  = N_day_out  + (N_total * 14.0_wp * 1.0E-6_wp) * weight * 86400.0_wp
      ! ----------------------------------------------------------------- !
      ! These five are DAILY MEANS: the weighted average SUM(x*weight) with
      ! SUM(weight) = 1, the weight being 1/24 over the 24 hourly steps, which is
      ! also why t_i_day above applies the 86400 s d-1 conversion exactly once.
      ! The temperature and vapor-pressure-deficit accumulators below carry a
      ! factor 24 and are divided by 24 again in main_MMCO2Crop.f90; these five
      ! never were, and the daylight means a few lines below divide by the
      ! daylight weight sum instead.
      ! ----------------------------------------------------------------- !
      gsw_mm_day_out = gsw_mm_day_out + gsw_An * weight
      rca_day_out = rca_day_out + rca_s * weight
      rcs_day_out = rcs_day_out + rcs_p * weight
      rsa_day_out = rsa_day_out + rsa_s * weight
      rss_day_out = rss_day_out + rss_s * weight
      ! ----------------------------------------------------------------- !
      ! Daylight gate on the component available energy.
      ! ----------------------------------------------------------------- !
      is_day = ((Rnc_s + Rns_s) > 0.0_wp)
      IF (is_day) THEN
        day_weight_sum = day_weight_sum + weight
        gsw_day_sum = gsw_day_sum + gsw_An * weight
        rca_day_sum = rca_day_sum + rca_s * weight
        rcs_day_sum = rcs_day_sum + rcs_p * weight
        rsa_day_sum = rsa_day_sum + rsa_s * weight
        rss_day_sum = rss_day_sum + rss_s * weight
        t_i_daylight = t_i_daylight + (LEc_s / L_Temp_h) * mass_conv * weight
        e_i_daylight = e_i_daylight + (LEs_s / L_Temp_h) * mass_conv * weight
        ! --------------------------------------------------------------- !
        ! Intercellular CO2 under the same daylight mask as the 
        ! conductance and transpiration means.
        ! --------------------------------------------------------------- !
        ci_day_sum = ci_day_sum + Ci_step * weight
      END IF
      Ci_day_out = Ci_day_out + Ci_step * weight

      G_day_out   = G_day_out  + dHour_to_dDay (G_s , weight)
      LE_day_out  = LE_day_out  + dHour_to_dDay (LE_s, weight)
      LEc_day_out = LEc_day_out + dHour_to_dDay (LEc_s, weight)
      LEs_day_out = LEs_day_out + dHour_to_dDay (LEs_s, weight)
      H_day_out   = H_day_out   + dHour_to_dDay (H_total_s , weight)
      Hc_day_out  = Hc_day_out  + dHour_to_dDay (Hc_s, weight)
      Hs_day_out  = Hs_day_out  + dHour_to_dDay (Hs_s, weight)
      ! ----------------------------------------------------------------- !
      ! Reported daily net radiation from the COMPONENT system.
      ! ----------------------------------------------------------------- !
      Rn_day_out  = Rn_day_out + dHour_to_dDay (Rnc_s + Rns_s, weight)
      Ac_day_out  = Ac_day_out + dHour_to_dDay (ac_s, weight)
      As_day_out  = As_day_out + dHour_to_dDay (as_s, weight)

      T_leaf_day_out = T_leaf_day_out + T_leaf_s * weight * 24.0_wp
      T_air_day_out  = T_air_day_out  + HourTemp_s * weight * 24.0_wp
      T_soil_day_out = T_soil_day_out + Tsoil_surf_s * weight * 24.0_wp
      VPD_day_out    = VPD_day_out    + Ds_s * weight * 24.0_wp

      IF (do_write) THEN
        ! --------------------------------------------------------------- !
        ! Hourly diurnal energy-balance row: 3 index columns and 27 payload
        ! columns, the last seven being ground-area net and gross assimilation,
        ! nitrate reduced, intercellular CO2, sunlit and shaded leaf area index
        ! and leaf area index. The header is written in main_MMCO2Crop.f90 and
        ! must match this list.
        ! --------------------------------------------------------------- !
        WRITE (unit_diurnal,'(I0,",",I0,",",I0,27(",",G14.6))') &
          day_idx, J, iter_eb, hour_val, HourTemp_s, T_leaf_s, Ds_s, Rn_s, &
          LE_s, LEc_s, LEs_s, EB_err, Hc_s, Hs_s, H_total_s, &
          gsw_An, rca_s, rcs_p, rsa_s, rss_s, ac_s, as_s, VpdMcf_s, &
          Assim_net_total, Assim_total, N_total, Ci_step, &
          sunlit_lai_s, shaded_lai_s, gParam_lai
      END IF
      ! ----------------------------------------------------------------- !
      ! Top-of-canopy gas exchange for leaf-level comparison: the canopy-top
      ! Narea n0 from the profile at the head of this module, the sunlit absorbed
      ! PAR, and the converged canopy temperature and leaf-surface vapor pressure
      ! deficit, through the same stomatal model. DIAGNOSTIC ONLY - it is never
      ! added to the canopy totals, and par is saved and restored so the call has
      ! no side effect. The class tag is set around the call so the limiting-limb
      ! census inside assim can exclude it: it is taken once per committed time
      ! step, inside the same writer block that writes the hourly row above and
      ! after the temperature solve has settled, so the census records committed
      ! states and not solver trials.
      ! ----------------------------------------------------------------- !
      IF (do_write .AND. energy_balance_mode == 2 .AND. flag_diag_open .AND. &
          kN_canopy > 1.0e-6_wp .AND. gParam_lai > 1.0e-6_wp .AND. N_leaf > 1.0e-9_wp) THEN
        n0_flag  = canopy_n0_top(N_leaf, gParam_lai, kN_canopy)   ! records 76 and 94 select the profile
        par_save = par
        par      = Absor_IPAR_sunlit_h_s * 4.57_wp
        ci_census_tag = 3; ci_census_day = day_idx; ci_census_j = J; ci_census_hour = hour_val
        CALL assim (Rd_f, A_f, An_f, gsw_f, Ca_s, Ci_f, N_f, Ds_s, Narea_in=n0_flag)
        ci_census_tag = 0
        par      = par_save
        ! --------------------------------------------------------------- !
        ! Flag-leaf day respiration rides along from the same assim call
        ! (umol CO2 m-2 s-1), diagnostic only.
        ! --------------------------------------------------------------- !
        WRITE (unit_flag_diag,'(I0,",",I0,8(",",G14.6))') day_idx, J, hour_val, &
          n0_flag, An_f, gsw_f, Ci_f/Ca_s, T_leaf_s, Ds_s, Rd_f
      END IF
      ! ----------------------------------------------------------------- !
      ! Hourly three-node energy-balance closure diagnostics. These use the
      ! SOLVED, pre-repartition component fluxes, so the residuals reflect the
      ! solver itself rather than the committed latent heat.
      ! ----------------------------------------------------------------- !
      IF (do_write .AND. energy_balance_mode >= 1 .AND. eb_diag_open) THEN
        ! --------------------------------------------------------------- !
        ! Buffer the row: the file is written after the repartition, from the
        ! ACTUAL state, with these demand values kept as the _demand columns.
        ! --------------------------------------------------------------- !
        ebh_row(:, J) = (/ hour_val, HourTemp_s, Tac_s, Tcan_s, Tsoil_surf_s, &
          VpdRef_s, VPD_can_s, Ds_leaf_s, e_ac_s, esat_pa(Tcan_s), &
          Assim_net_total, gsw_An, rcs_p, raa_s, raa_neutral_s, rca_s, rsa_s, rss_s, &
          L_obukhov_s, zeta_s, psih_s, psim_s, ustar_s, REAL(stab_regime_s, wp), &
          Rnc_s, Rns_s, (Rnc_s + Rns_s), LEc_solved_s, LEs_solved_s, Hc_solved_s, Hs_solved_s, &
          (Hc_solved_s + Hs_solved_s), G_s, (Tac_s - HourTemp_s), (Tcan_s - HourTemp_s), &
          F_can_s, F_soil_s, F_air_s, F_vap_s, h_soil_s, &
          REAL(eb3_inner_iters, wp), REAL(eb3_outer_iters, wp), REAL(eb3_conv_flag, wp) /)
      END IF
      ! ----------------------------------------------------------------- !
      ! Buffer the demand-state radiation row; the radiation file is written
      ! after the repartition, from the committed scene, with these values
      ! appended as the _demand columns.
      ! ----------------------------------------------------------------- !
      IF (do_write .AND. energy_balance_mode >= 1 .AND. rad_open) THEN
        rad_row(:, J) = (/ hour_val, rad_SW_down, rad_SW_abs_can, rad_SW_abs_soil, &
          rad_SW_ref_can, rad_SW_ref_soil, rad_LW_down, rad_LW_emit_can, rad_LW_emit_soil, &
          rad_LW_can_to_soil, rad_LW_soil_to_can, rad_LW_can_net, rad_LW_soil_net, &
          rad_LW_up_sky, &
          Rnc_s, Rns_s, (Rnc_s + Rns_s), rad_alpha_soil, alpha_canopy_wheat, rad_alpha_mix, &
          rad_tau_sw, rad_tau_lw /)
      END IF
      ! ----------------------------------------------------------------- !
      ! Store this substep's POTENTIAL profile, so that the post-Richards
      ! repartition can redistribute the day's actual fluxes over the hours and
      ! recompute the actual three-node temperatures. The resistances are
      ! daily-constant, so storing them per substep is exact.
      ! ----------------------------------------------------------------- !
      IF (energy_balance_mode >= 1 .AND. J <= EB_MAXSTEP) THEN
        ebh_weight(J)        = weight
        ebh_LEc_trans_pot(J) = LEc_pot_unc_s   ! potential, not the capped state
        ! --------------------------------------------------------------- !
        ! Wet-canopy (interception) potential evaporation: the stomata removed,
        ! latent heat passing through the leaf boundary layer only. Zero when
        ! there is no canopy, the boundary-layer resistance then being huge.
        ! --------------------------------------------------------------- !
        ebh_LEc_int_pot(J)   = MAX(0.0_wp, (rho * Cp / gamma_const) * &
                               (esat_pa(Tcan_s) - e_ac_s) / MAX(1.0_wp, rca_s))
        ebh_LEs_pot(J)       = LEs_pot_unc_s   ! potential, not the capped state
        ebh_Rnc(J)           = Rnc_s
        ebh_Rns(J)           = Rns_s
        ebh_G(J)             = G_s
        ebh_rca(J)           = rca_s
        ebh_rsa(J)           = rsa_s
        ebh_raa_neutral(J)   = raa_neutral_s
        ebh_Tair(J)          = HourTemp_s
        ebh_wind(J)          = w_eb_s            ! the day's mean wind speed
        ebh_hour(J)          = hour_val
        ebh_Tcan_pot(J)      = Tcan_s
        ebh_lai(J)           = gParam_lai
        ebh_hgt(J)           = gParam_hgt_s
        ebh_SW_in(J)         = MAX(0.0_wp, dfrad_s + drrad_s)
        ! --------------------------------------------------------------- !
        ! The component path uses the Brutsaert clear-sky emissivity with the
        ! Crawford and Duchon cloud modification, set in meteo.
        ! --------------------------------------------------------------- !
        ebh_LW_down(J)       = sky_emiss_comp * const_SB * ((HourTemp_s + 273.15_wp)**4)
        ebh_theta_top(J)     = theta(1)
        ! --------------------------------------------------------------- !
        ! This substep's solar geometry, stored for the actual-flux replay.
        ! --------------------------------------------------------------- !
        ebh_Kdr(J)           = Kdr_s
        ebh_SolarElev(J)     = SolarElev_s
        ! --------------------------------------------------------------- !
        ! Driver record 95: diagnostic sunlit and shaded leaf temperature 
        ! at the converged demand state. DIAGNOSTIC ONLY - it writes 
        ! nothing that any flux, state or solver reads. 
        ! Full provenance at SHARED.f90 H-1..H-6.
        ! --------------------------------------------------------------- !
        IF (tleaf_class_diag >= 1) CALL tleaf_classes_diag(J)
        ! --------------------------------------------------------------- !
        ! Daily interception-evaporation demand (mm/d), 
        ! same conversion as t_i_day.
        ! --------------------------------------------------------------- !
        LEc_int_pot_day = LEc_int_pot_day + (ebh_LEc_int_pot(J) / L_Temp_h) * mass_conv * weight * 1000.0_wp
      END IF
    END DO
    IF (energy_balance_mode >= 1) eb_nsteps_s = MIN(nsteps, EB_MAXSTEP)

    IF (day_weight_sum > 0.0_wp) THEN
      ci_daylight_mean_diag = ci_day_sum / day_weight_sum
      gsw_day_out      = gsw_day_sum / day_weight_sum
      rca_daylight_out = rca_day_sum / day_weight_sum
      rcs_daylight_out = rcs_day_sum / day_weight_sum
      rsa_daylight_out = rsa_day_sum / day_weight_sum
      rss_daylight_out = rss_day_sum / day_weight_sum
    END IF
    tr_daylight_out    = t_i_daylight * 1000.0_wp
    es_daylight_out    = e_i_daylight * 1000.0_wp

  CONTAINS

    SUBROUTINE solve_biophysics(T_guess, Ds_in, EB_residual, H_tot, LE_tot, LEc_val, LEs_val, gsw_canopy, Ds_out, &
                                Ci_out, VpdMcf_out)
      REAL(wp), INTENT(IN) :: T_guess, Ds_in
      REAL(wp), INTENT(OUT) :: EB_residual, H_tot, LE_tot, LEc_val, LEs_val, gsw_canopy, Ds_out, Ci_out
      REAL(wp), INTENT(OUT) :: VpdMcf_out

      lai_sum = sunlit_lai_s + shaded_lai_s
      T_leaf_s = T_guess ! Update shared state for assim

      ! ---------------------------------------------------------------------- !
      ! energy_balance_mode >= 1: the three-node canopy-soil-airspace energy
      ! balance. The outer root-find on the canopy temperature is the caller's;
      ! here, for THIS trial temperature, a nested fixed point is converged so
      ! that assimilation, stomatal conductance and the bulk stomatal resistance
      ! are mutually consistent with the solved canopy-air vapor pressure BEFORE
      ! the canopy residual is returned. Without it the energy balance would be
      ! another post-hoc coupling, the conductance seeing a stale vapor pressure
      ! deficit.
      ! ---------------------------------------------------------------------- !
      IF (energy_balance_mode >= 1) THEN
        BLOCK
          REAL(wp) :: raa, raa0, raa_prev, rsa, rca, rss, rcs, gsw_eff, rst
          REAL(wp) :: SW_in_s, LW_down_s, theta_top_s
          REAL(wp) :: Ds_leaf, Ds_leaf_new, e_ac, Tac, Tsoil, e_ref, H_est, Tair_K
          REAL(wp) :: raa_neut, Lob, zeta, psih, psim, ustar
          INTEGER  :: stabreg
          REAL(wp) :: LEc_l, LEs_l, Hc_l, Hs_l, G_l, Rnc_l, Rns_l
          REAL(wp) :: Fsoil, Fair, Fvap, Fcan
          INTEGER  :: kfp, inner_it
          LOGICAL  :: inner_conv, fp_conv
          REAL(wp), PARAMETER :: fp_tol = 1.0_wp   ! Pa tolerance on leaf-surface VPD
          INTEGER,  PARAMETER :: fp_max = 40
          ! ---------------------------------------------------------------- !
          ! Neutral resistances (also sets disp_height_s/rough_len_s used by the
          ! MO correction). rsa/rca/rss stay neutral; only the above-canopy heat
          ! resistance raa receives the Monin-Obukhov stability correction below.
          ! ---------------------------------------------------------------- !
          CALL compute_eb_resistances(raa0, rsa, rca, rss)
          e_ref  = HourVP_s
          Tair_K = HourTemp_s + 273.15_wp
          ! ---------------------------------------------------------------- !
          ! Component-radiation inputs: incoming shortwave, and sky downwelling
          ! longwave from the Brutsaert clear-sky emissivity with the Crawford
          ! and Duchon cloud modification, set in meteo from the wired humidity.
          ! The top-layer water content drives the soil albedo.
          ! ---------------------------------------------------------------- !
          SW_in_s     = MAX(0.0_wp, dfrad_s + drrad_s)
          LW_down_s   = sky_emiss_comp * const_SB * (Tair_K**4)
          theta_top_s = theta(1)
          Tac    = HourTemp_s; Tsoil = HourTemp_s; e_ac = e_ref   ! warm start
          Ds_leaf = MAX(0.0_wp, Ds_in)
          H_est   = 0.0_wp                                        ! neutral on first pass
          raa = raa0; raa_prev = raa0; raa_neut = raa0
          Lob = 1.0e6_wp; zeta = zero; psih = zero; psim = zero; ustar = zero
          fp_conv = .FALSE.; inner_conv = .FALSE.; inner_it = 0
          LEc_l = zero; LEs_l = zero; Hc_l = zero; Hs_l = zero
          G_l = zero; Rnc_l = zero; Rns_l = zero
          Fsoil = zero; Fair = zero; Fvap = zero

          DO kfp = 1, fp_max
            ! -------------------------------------------------------------- !
            ! Stability-corrected above-canopy heat resistance from the current H
            ! (neutral on the first pass; couples raa<->H within the fixed point).
            ! Bare soil (no canopy): MO canopy formula is N/A -> keep neutral raa0.
            ! -------------------------------------------------------------- !
            IF (gParam_lai > 0.0_wp .AND. gParam_hgt_s > 0.0_wp) THEN
              raa = raa_heat_stability(w_eb_s, H_est, Tair_K, raa_neut, Lob, zeta, psih, psim, ustar, stabreg)   ! the day's mean wind speed
            ELSE
              raa = raa0; raa_neut = raa0; stabreg = 0
            END IF
            ! -------------------------------------------------------------- !
            ! assim at the current leaf-surface VPD -> canopy gsw -> rcs.
            ! -------------------------------------------------------------- !
            CALL eval_canopy_assim(Ds_leaf, gsw_canopy, Ci_out)
            gsw_An  = gsw_canopy
            gsw_eff = MAX(1.0e-6_wp, gsw_An)
            ! -------------------------------------------------------------- !
            ! gsw (mol m-2 s-1) -> rst (s/m); rcs = rst/effLAI.
            ! -------------------------------------------------------------- !
            rst = 1.0_wp / (gsw_eff * (gasc * (273.15_wp + HourTemp_s)) / 101325.0_wp)
            rcs = rst / MAX(0.01_wp, dLe_s)

            CALL solve_eb_3node_inner(T_guess, rcs, raa, rsa, rca, rss, SW_in_s, LW_down_s, &
                   gParam_lai, theta_top_s, HourTemp_s, e_ref, &
                   Tsoil, Tac, e_ac, LEc_l, LEs_l, Hc_l, Hs_l, G_l, Rnc_l, Rns_l, &
                   Fsoil, Fair, Fvap, inner_it, inner_conv)
            ! -------------------------------------------------------------- !
            ! Relax the H that drives the stability correction. The raa<->H feedback
            ! can otherwise oscillate on near-closed-stomata hot full-canopy hours.
            ! -------------------------------------------------------------- !
            H_est = 0.5_wp * (Hc_l + Hs_l) + 0.5_wp * H_est

            Ds_leaf_new = MAX(0.0_wp, esat_pa(T_guess) - e_ac)
            IF (ABS(Ds_leaf_new - Ds_leaf) < fp_tol .AND. ABS(raa - raa_prev) < 1.0_wp &
                .AND. kfp > 1) THEN
              Ds_leaf = Ds_leaf_new
              fp_conv = inner_conv
              EXIT
            END IF
            raa_prev = raa
            Ds_leaf = 0.5_wp * Ds_leaf_new + 0.5_wp * Ds_leaf   ! relaxation
          END DO

          Fcan = Rnc_l - Hc_l - LEc_l
          ! ---------------------------------------------------------------- !
          ! Solved (pre-hydrology-scaling) state -> SHARED for closure diagnostics.
          ! ---------------------------------------------------------------- !
          rca_s = rca; rsa_s = rsa; rss_s = rss; raa_s = raa; rcs_p = rcs
          raa_neutral_s = raa_neut; L_obukhov_s = Lob; zeta_s = zeta
          psih_s = psih; psim_s = psim; ustar_s = ustar; stab_regime_s = stabreg
          ac_s = Rnc_l; as_s = Rns_l - G_l; G_s = G_l
          Hc_s = Hc_l; Hs_s = Hs_l
          Tsoil_surf_s = Tsoil; T_leaf_s = T_guess; Tcan_s = T_guess; Tac_s = Tac
          e_ac_s = e_ac; VPD_can_s = MAX(0.0_wp, esat_pa(T_guess) - e_ac); Ds_leaf_s = Ds_leaf
          h_soil_s = h_soil_surface(Tsoil)
          LEc_solved_s = LEc_l; LEs_solved_s = LEs_l; Hc_solved_s = Hc_l; Hs_solved_s = Hs_l
          ! ---------------------------------------------------------------- !
          ! The uncapped network soil latent heat at the solved state is the
          ! POTENTIAL.
          ! ---------------------------------------------------------------- !
          LEs_pot_unc_s = LEs_unc_last
          LEc_pot_unc_s = LEc_unc_last   ! uncapped canopy potential
          Rnc_s = Rnc_l; Rns_s = Rns_l
          F_can_s = Fcan; F_soil_s = Fsoil; F_air_s = Fair; F_vap_s = Fvap
          eb3_inner_iters = inner_it; eb3_outer_iters = kfp
          eb3_conv_flag = 0
          IF (fp_conv) eb3_conv_flag = 1

          Ds_out = VPD_can_s; VpdMcf_out = VPD_can_s
          EB_residual = Fcan          ! objective of the outer canopy-temperature solve (W m-2)
          ! ---------------------------------------------------------------- !
          ! The demand pass returns POTENTIAL latent heat, with no post-hoc
          ! scaling. Water limitation is applied afterwards, by the hydrology
          ! actual-flux repartition, not here. The returned values are the
          ! UNCAPPED network ones; the capped state stays in LEc_solved_s and
          ! LEs_solved_s for the closure diagnostics.
          ! ---------------------------------------------------------------- !
          LEc_val = LEc_pot_unc_s
          LEs_val = LEs_pot_unc_s
          LE_tot  = LEc_val + LEs_val
          H_tot   = Hc_l + Hs_l
        END BLOCK
      END IF

    END SUBROUTINE solve_biophysics
    ! ---------------------------------------------------------------------- !
    ! Sunlit plus shaded assimilation with the two-leaf canopy integration
    ! (Goudriaan 1977; de Pury and Farquhar 1997). Factored out so the three-node
    ! fixed point can re-evaluate it at an updated leaf-surface vapor pressure
    ! deficit. Leaf temperature and absorbed PAR come from SHARED.
    ! ---------------------------------------------------------------------- !
    SUBROUTINE eval_canopy_assim(Ds_in, gsw_canopy, Ci_out)
      REAL(wp), INTENT(IN)  :: Ds_in
      REAL(wp), INTENT(OUT) :: gsw_canopy, Ci_out
      REAL(wp) :: A_s_local, An_s_local, gsw_s_local, Rd_local, N_local, Ci_local
      REAL(wp) :: Ci_sun, Ci_shade
      REAL(wp) :: n0_canN, Nsun_canN, Narea_sun, Narea_shd, kb_canN, canN_resid
      LOGICAL  :: use_canN

      lai_sum = sunlit_lai_s + shaded_lai_s
      IF (lai_sum > 1.0e-6_wp) THEN
        ! ----------------------------------------------------------------- !
        ! Per-class mean Narea from the canopy leaf-nitrogen profile, with the
        ! total lamina nitrogen conserved: the sunlit class takes the
        ! beam-weighted integral of n(l) divided by the sunlit leaf area index,
        ! the shaded class the remainder divided by the shaded leaf area index.
        ! The profile itself is canopy_n0_top and canopy_n_sunlit at the head of
        ! this module and has two branches, selected by driver records 76 and 94;
        ! their forms and sources are stated there and at SHARED.f90 N-1..N-6 and
        ! items E-1..E-7 of the revision record (internal; not distributed). 
        ! The beam weight is exp(-kb l) with kb = Kdr from solar.f90,
        ! which also defines sunlit_lai = (1 - exp(-kb LAI))/kb, the quantity the
        ! first term of the sunlit integral reduces to. The conservation residual
        ! is formed and warned on, never corrected.
        ! ----------------------------------------------------------------- !
        use_canN=.FALSE.; Narea_sun=-1.0_wp; Narea_shd=-1.0_wp; canN_resid=0.0_wp
        IF (energy_balance_mode == 2 .AND. kN_canopy > 1.0e-6_wp .AND. gParam_lai > 1.0e-6_wp &
            .AND. sunlit_lai_s > 1.0e-6_wp .AND. shaded_lai_s > 1.0e-6_wp &
            .AND. N_leaf > 1.0e-9_wp .AND. Kdr_s > 1.0e-6_wp) THEN
          kb_canN   = Kdr_s
          n0_canN   = canopy_n0_top  (N_leaf, gParam_lai, kN_canopy)            ! records 76 and 94 select the profile
          canE_n_top = n0_canN   ! diagnostic: under record 94 modes 1-3 this is the flat top-zone Narea
          Nsun_canN = canopy_n_sunlit(N_leaf, gParam_lai, kN_canopy, kb_canN)
          Narea_sun = Nsun_canN / sunlit_lai_s
          Narea_shd = (N_leaf - Nsun_canN) / shaded_lai_s
          use_canN  = (Narea_sun > 0.0_wp .AND. Narea_shd > 0.0_wp)
          canN_resid = Narea_sun*sunlit_lai_s + Narea_shd*shaded_lai_s - N_leaf
          IF (use_canN .AND. do_write .AND. ABS(canN_resid) > 1.0e-6_wp) &
            WRITE (*,'(A,ES12.4)') ' WARN canopy-N conservation resid (g N m-2): ', canN_resid
        END IF
        par = Absor_IPAR_sunlit_h_s * 4.57_wp
        ! ----------------------------------------------------------------- !
        ! Class tag for the limiting-limb census inside assim; diagnostic only,
        ! and the rows are written on the writer pass only.
        ! ----------------------------------------------------------------- !
        IF (do_write) THEN
          ci_census_tag = 1; ci_census_day = day_idx; ci_census_j = J; ci_census_hour = hour_val
        END IF
        CALL assim (Rd_local, A_s_local, An_s_local, gsw_s_local, Ca_s, Ci_local, N_local, Ds_in, Narea_in=Narea_sun)
        A_can_sun = A_s_local * sunlit_lai_s; An_can_sun = An_s_local * sunlit_lai_s
        gsw_sun = gsw_s_local; Rd_s_can_sun = Rd_local * sunlit_lai_s
        N_can_sun = N_local * sunlit_lai_s; Ci_sun = Ci_local

        par = Absor_IPAR_shaded_h_s * 4.57_wp
        IF (do_write) ci_census_tag = 2
        CALL assim (Rd_local, A_s_local, An_s_local, gsw_s_local, Ca_s, Ci_local, N_local, Ds_in, Narea_in=Narea_shd)
        ci_census_tag = 0
        A_can_shade = A_s_local * shaded_lai_s; An_can_shade = An_s_local * shaded_lai_s
        gsw_shade = gsw_s_local; Rd_s_can_shade = Rd_local * shaded_lai_s
        N_can_shade = N_local * shaded_lai_s; Ci_shade = Ci_local

        gsw_canopy = (gsw_sun * sunlit_lai_s + gsw_shade * shaded_lai_s) / lai_sum
        Ci_out = (Ci_sun * sunlit_lai_s + Ci_shade * shaded_lai_s) / lai_sum
      ELSE
        A_can_sun = zero; An_can_sun = zero; A_can_shade = zero; An_can_shade = zero
        gsw_canopy = zero; Rd_s_can_sun = zero; Rd_s_can_shade = zero; N_can_sun = zero; N_can_shade = zero
        Ci_out = zero
      END IF

      Assim_total = A_can_sun + A_can_shade; Assim_net_total = An_can_sun + An_can_shade
      Rd_total = Rd_s_can_sun + Rd_s_can_shade; N_total = N_can_sun + N_can_shade
      gsw_An = gsw_canopy
    END SUBROUTINE eval_canopy_assim

    ! --------------------------------------------------------------------- !
    ! Driver record 95: DIAGNOSTIC sunlit and shaded leaf temperature. It solves
    ! the model's OWN canopy-node energy balance divided by the leaf area index,
    ! that is the per-leaf form, separately for the sunlit and the shaded class,
    ! at the SAME converged canopy-air temperature and vapor pressure and the
    ! same per-leaf boundary layer, each class with its own absorbed net
    ! radiation and its own stomatal resistance. NO NEW CONSTANT and NO FEEDBACK:
    ! every result goes to the hourly diagnostic row and the two class
    ! temperatures in SHARED and nowhere else. Full provenance SHARED.f90
    ! H-1..H-6.
    ! --------------------------------------------------------------------- !
    SUBROUTINE tleaf_classes_diag(Jslot)
      INTEGER, INTENT(IN) :: Jslot
      REAL(wp) :: L_tot, conv_h, rb_h, rst_sun_h, rst_shd_h, fcap_h
      REAL(wp) :: w_sun_h, w_shd_h, wsum_h, sw_sun_g, sw_shd_g, lw_sun_g, lw_shd_g
      REAL(wp) :: rn_sun_l, rn_shd_l, rn_bulk_l, t_sun_h, t_shd_h, t_chk_h

      L_tot = gParam_lai
      T_leaf_sun_s = Tcan_s; T_leaf_shd_s = Tcan_s
      ebh_tlc(:, Jslot) = 0.0_wp
      ebh_tlc(10, Jslot) = Tac_s; ebh_tlc(11, Jslot) = e_ac_s
      ebh_tlc(12, Jslot) = Tcan_s; ebh_tlc(13, Jslot) = L_tot
      ebh_tlc(1, Jslot) = Tcan_s; ebh_tlc(2, Jslot) = Tcan_s; ebh_tlc(9, Jslot) = Tcan_s
      IF (L_tot <= 1.0e-6_wp .OR. sunlit_lai_s <= 1.0e-9_wp .OR. shaded_lai_s <= 1.0e-9_wp &
          .OR. rca_s <= 0.0_wp .OR. rca_s > 1.0e5_wp) RETURN
      ! ------------------------------------------------------------------- !
      ! Per-leaf resistances: rca = rb/L and rcs = rst/L, so rb = rca*L. The
      ! class stomatal resistances come from the two-leaf assimilation's OWN
      ! gsw_sun / gsw_shade, converted exactly as the bulk rst is converted.
      ! ------------------------------------------------------------------- !
      conv_h    = (gasc * (273.15_wp + HourTemp_s)) / 101325.0_wp
      rb_h      = rca_s * L_tot
      ! ------------------------------------------------------------------- !
      ! LATENT BASIS. The bulk canopy node is solved against the CAPPED canopy
      ! latent heat wherever the supply cap binds (energy_balance_module takes
      ! the minimum of the network value and the cap), so by definition
      !     LEc_solved = fcap x LEc_potential,  fcap = LEc_solved/LEc_potential,
      ! and the bulk node equation, divided by the leaf area index, is the
      ! per-leaf form WITH fcap on the latent term. The class solves therefore
      ! carry the SAME fcap: the cap is a whole-canopy water-supply limit and the
      ! model has no per-class supply, so the limitation is shared in proportion
      ! to each class's own potential demand. That is a DECLARED CHOICE, and it
      ! is the choice that makes the sunlit-minus-canopy difference a difference
      ! within ONE state rather than a comparison of a supply-limited node with
      ! an unlimited leaf.
      ! ------------------------------------------------------------------- !
      fcap_h = 1.0_wp
      IF (LEc_pot_unc_s > 1.0e-9_wp) &
        fcap_h = MAX(0.0_wp, MIN(1.0_wp, LEc_solved_s / LEc_pot_unc_s))
      rst_sun_h = 1.0_wp / (MAX(1.0e-6_wp, gsw_sun)   * conv_h)
      rst_shd_h = 1.0_wp / (MAX(1.0e-6_wp, gsw_shade) * conv_h)
      ! ------------------------------------------------------------------- !
      ! Class net radiation, conserving Rnc = sw_can_abs + lw_can_net EXACTLY.
      ! ------------------------------------------------------------------- !
      w_sun_h = Absor_IPAR_sunlit_h_s * sunlit_lai_s
      w_shd_h = Absor_IPAR_shaded_h_s * shaded_lai_s
      wsum_h  = w_sun_h + w_shd_h
      IF (wsum_h > 1.0e-12_wp) THEN
        sw_sun_g = rad_SW_abs_can * (w_sun_h / wsum_h)
      ELSE
        sw_sun_g = rad_SW_abs_can * (sunlit_lai_s / L_tot)   ! night / no beam: area split
      END IF
      sw_shd_g = rad_SW_abs_can - sw_sun_g
      lw_sun_g = rad_LW_can_net * (sunlit_lai_s / L_tot)
      lw_shd_g = rad_LW_can_net - lw_sun_g
      rn_sun_l  = (sw_sun_g + lw_sun_g) / sunlit_lai_s
      rn_shd_l  = (sw_shd_g + lw_shd_g) / shaded_lai_s
      rn_bulk_l = (rad_SW_abs_can + rad_LW_can_net) / L_tot

      t_sun_h = leaf_T_newton(rn_sun_l, rb_h, rst_sun_h, fcap_h)
      t_shd_h = leaf_T_newton(rn_shd_l, rb_h, rst_shd_h, fcap_h)
      ! ------------------------------------------------------------------- !
      ! Falsifiability check: the SAME solver on the BULK inputs must return the
      ! canopy temperature (SHARED.f90 H-5).
      ! ------------------------------------------------------------------- !
      t_chk_h = leaf_T_newton(rn_bulk_l, rb_h, rcs_p * L_tot, fcap_h)

      T_leaf_sun_s = t_sun_h; T_leaf_shd_s = t_shd_h
      ebh_tlc(1, Jslot) = t_sun_h
      ebh_tlc(2, Jslot) = t_shd_h
      ebh_tlc(3, Jslot) = sunlit_lai_s / L_tot
      ebh_tlc(4, Jslot) = rn_sun_l
      ebh_tlc(5, Jslot) = rn_shd_l
      ebh_tlc(6, Jslot) = rb_h
      ebh_tlc(7, Jslot) = rst_sun_h
      ebh_tlc(8, Jslot) = rst_shd_h
      ebh_tlc(9, Jslot) = t_chk_h
      ebh_tlc(14, Jslot) = 1.0_wp
      ! ------------------------------------------------------------------- !
      ! The bulk canopy node is solved against the capped canopy latent heat when
      ! the supply cap binds, while this per-leaf solve uses the uncapped network
      ! form. The gap is written out so that the bulk check can be ATTRIBUTED and
      ! is never offered as a clean identity (SHARED.f90 H-5).
      ! ------------------------------------------------------------------- !
      ebh_tlc(15, Jslot) = LEc_pot_unc_s - LEc_solved_s
      ebh_tlc(16, Jslot) = fcap_h

    END SUBROUTINE tleaf_classes_diag
    ! --------------------------------------------------------------------- !
    ! Helper of the class diagnostic. Solves, for T,
    !   Rn_leaf = rho cp (T - Tac)/rb + fcap (rho cp/gamma)(esat(T) - e_ac)/(rb + rst)
    ! by damped Newton with a numerical d(esat)/dT, at the converged demand-state
    ! canopy-air temperature and vapor pressure. This is the model's OWN
    ! canopy-node equation divided by the leaf area index (SHARED.f90 H-2),
    ! called once per class and once on the bulk inputs, where it must return the
    ! canopy temperature (SHARED.f90 H-5).
    ! --------------------------------------------------------------------- !
    FUNCTION leaf_T_newton(rn_leaf, rb_in, rst_in, fcap_in) RESULT(t_out)
      REAL(wp), INTENT(IN) :: rn_leaf, rb_in, rst_in, fcap_in
      REAL(wp) :: t_out, fv, dfv, dsdt, step, rsum, rcp_l, cle_l, rb_l
      INTEGER  :: itl
      rcp_l = rho * Cp
      cle_l = rcp_l / gamma_const
      rb_l  = MAX(1.0_wp, rb_in)
      rsum  = MAX(1.0_wp, rb_in + rst_in)
      t_out = Tcan_s
      DO itl = 1, 40
        dsdt = (esat_pa(t_out + 0.01_wp) - esat_pa(t_out - 0.01_wp)) / 0.02_wp
        fv   = rn_leaf - rcp_l * (t_out - Tac_s) / rb_l &
                       - fcap_in * cle_l * (esat_pa(t_out) - e_ac_s) / rsum
        dfv  = - rcp_l / rb_l - fcap_in * cle_l * dsdt / rsum
        IF (ABS(dfv) < 1.0e-12_wp) EXIT
        step = -fv / dfv
        step = MAX(-10.0_wp, MIN(10.0_wp, step))
        t_out = t_out + step
        t_out = MAX(Tac_s - 40.0_wp, MIN(Tac_s + 50.0_wp, t_out))
        IF (ABS(step) < 1.0e-5_wp) EXIT
      END DO
    END FUNCTION leaf_T_newton


  END SUBROUTINE compute_daily_canopy_fluxes

END MODULE CANOPY_MODEL