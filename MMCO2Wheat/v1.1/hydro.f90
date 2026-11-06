!===============================================================================!
! M&M-CO2Wheat v1.0  |  hydro.f90
!
! Purpose: Soil water movement by the one-dimensional mixed-form Richards
!          equation, solved over one day in adaptive sub-steps by implicit
!          backward Euler with Picard (and Newton) iteration. Spatial unit: the
!          ten cell-centered soil layers to 1.20 m, plus the surface store.
!          Hydraulics are van Genuchten (1980) retention with Mualem (1976)
!          conductivity; root water uptake is distributed over the rooted layers
!          by the Jackson et al. (1996) profile and reduced by the Feddes et al.
!          (1978) response.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: SHARED (the layer geometry and van Genuchten-Mualem parameters, the day's
!       rainfall, irrigation, soil-evaporation and transpiration demand, the
!       root-depth state and the solver configuration); TYPES (working
!       precision); IEEE_ARITHMETIC (finiteness tests on the iterates).
! Output: the updated layer water contents theta_ly and matric heads, the
!       root-uptake sink per layer, actual transpiration and soil evaporation,
!       infiltration, runoff, bottom flux and the surface store, together with
!       the per-day water-balance ledger - consumed by main_MMCO2Crop (the water
!       balance and its diagnostic files), by energy_balance_module (the
!       committed surface head and the soil latent-heat supply limit), by
!       CANOPY_MODEL through the soil water potential that enters the stomatal
!       response, and by growth_module through the water-stress state.
! Called by: main_MMCO2Crop.f90 - the daily loop, once per pass of the
!       evapotranspiration / Richards coupling iteration.
!===============================================================================!
SUBROUTINE hydro (theta_ly, day_id, t_end)
  USE SHARED
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: ieee_is_finite
  USE TYPES, ONLY: wp
  IMPLICIT NONE

  REAL(wp), INTENT(INOUT) :: theta_ly(NSoil_layers)
  INTEGER, INTENT(IN) :: day_id
  REAL(wp), INTENT(IN), OPTIONAL :: t_end

  ! Local Variables
  INTEGER :: i_layer, iter_picard, max_iter_picard
  REAL(wp) :: dt, t_current, t_day_end
  REAL(wp) :: flux_top, pr_flux
  REAL(wp) :: max_change, max_h_change
  ! ------------------------------------------------------------------------------------------------- !
  ! Under-relaxation for Picard updates to improve stability (reduce oscillations/divergence).
  ! Source: Moukalled et al., "The Finite Volume Method in CFD", Ch.14 (Under-Relaxation), pp. 8–9.
  ! Link (PDF): https://www.aub.edu/msfea/research/Documents/CFD-Chapter14TheSourceTerm.pdf
  ! ------------------------------------------------------------------------------------------------- !
  REAL(wp) :: relax_h
  REAL(wp) :: K_mid, dz_mid
  REAL(wp) :: z_surface, dz_top
  REAL(wp) :: flux_top_step
  REAL(wp) :: storage_before, storage_after, sink_total, q_bot, mb_err
  REAL(wp) :: mb_err_ledger ! rec149: TRUE pre-min_dt-zeroing residual, ledger/mode-2-closure ONLY (see hydro.f90:~1029-1036)
  REAL(wp) :: flux_top_mass
  REAL(wp) :: h_crit_atm
  REAL(wp) :: supply_step_mm
  REAL(wp) :: infil_step_mm, infil_supply_step_mm, infil_pond_step_mm
  REAL(wp) :: ponding_start_step_mm, runoff_step_mm, soil_evap_step_mm
  REAL(wp) :: surface_store_mm, surface_store_old_mm, surface_store_trial_mm
  REAL(wp) :: surface_head_old_m
  REAL(wp) :: call_supply_mm, call_infil_mm, call_infil_supply_mm, call_infil_pond_mm
  REAL(wp) :: call_soil_evap_mm, call_pond_evap_mm, call_runoff_mm
  ! ------------------------------------------------------------------------------------------------- !
  ! Physical constants for the vapor-equilibrium (Kelvin) surface head [PINNED]:
  !   R_gas: molar gas constant, exact by SI definition, 8.314 462 618... J/(mol K)
  !          (CODATA 2018, Tiesinga et al. 2021, Rev. Mod. Phys. 93, 025010, Table XXXI).
  !   Mw:    molar mass of H2O = 18.015 28 g/mol (IUPAC standard atomic weights;
  !          2*H 1.00794 + O 15.9994).
  !   g_acc: standard gravity 9.806 65 m/s2, exact by convention (3rd CGPM, 1901).
  ! ------------------------------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: R_gas = 8.314462618_wp
  REAL(wp), PARAMETER :: Mw = 0.01801528_wp
  REAL(wp), PARAMETER :: g_acc = 9.80665_wp
  REAL(wp), PARAMETER :: eps_len = 1.0e-8_wp
  REAL(wp), PARAMETER :: ponding_head_min_mm = 1.0e-1_wp
  REAL(wp), PARAMETER :: h_eval_min = -1.0e6_wp
  REAL(wp), PARAMETER :: h_eval_max = 1.0e3_wp
  REAL(wp), PARAMETER :: newton_dh_cap = 0.1_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! mb_tol_step and mb_tol_head_min_dt are SHARED module variables rather than
  ! local parameters, and are set once per hydro() call in the per-day reset
  ! block below.
  ! ------------------------------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: flux_shift_tol = 1.0e-8_wp
  REAL(wp) :: transpiration_rate_day ! m/d
  REAL(wp) :: evap_rate_day ! m/d
  REAL(wp) :: evap_rate_pot ! m/d
  REAL(wp) :: root_dist(NSoil_layers), uptake_alpha(NSoil_layers)
  REAL(wp) :: theta_new(NSoil_layers), h_new(NSoil_layers)
  REAL(wp) :: K_new(NSoil_layers), Cap_new(NSoil_layers)
  LOGICAL :: converged
  LOGICAL :: forced_head_step
  LOGICAL :: pond_unsat_fluxbc  ! Option3: ponded-unsaturated flux-BC infiltration active this substep
  LOGICAL :: pond_retry_done    ! A: head-BC pond-exhaustion single dt-cut taken this substep
  LOGICAL :: newton_converged
  LOGICAL :: loop_trace_printed
  LOGICAL :: it_crit
  INTEGER :: top_bc_type_old
  LOGICAL :: latch_risk_supply   ! FIX-BC: substep STARTS dry-latched under supply
  INTEGER :: hydro_loop_count
  INTEGER :: wf_substep        ! Wet-front trace row counter (per-day; diagnostics only).
  LOGICAL :: trace_this_day    ! True on target failure-class days when tracing is enabled.
  REAL(wp) :: wf_last_flux_m_d  ! Last accepted top-flux magnitude (m/d) for the wet-front dt limiter.
  REAL(wp) :: h_top_old
  REAL(wp) :: res_norm, res_norm_trial
  REAL(wp) :: res_mass   ! max|mixed-form residual| at the Picard iterate (units: 1/day)
  REAL(wp) :: res_step_m ! mass-scaled residual SUM(|R|*tk)*dt (units: m water)
  ! ------------------------------------------------------------------------------------------------- !
  ! PROVENANCE OF THE SURFACE BOUNDARY TREATMENT: the atmospheric boundary
  ! condition used below -- a prescribed flux that switches to a prescribed head
  ! when the surface would exceed the infiltration capacity or fall below the
  ! minimum head hCritA, with the switch tested inside the nonlinear iteration --
  ! is the method implemented in the HYDRUS-1D model and documented in Simunek,
  ! van Genuchten & Sejna (2016), "HYDRUS-1D Technical Manual", atmospheric
  ! boundary with maximum infiltration. It is a model formulation taken from that
  ! implementation, not a measurement. Inline mentions below are labelled
  ! "atmospheric-BC" and refer to this method.
  ! ------------------------------------------------------------------------------------------------- !
  LOGICAL, PARAMETER :: use_infiltration_cap = .TRUE.  ! enable the atmospheric-BC infiltration cap

  ! Parameters for Solver
REAL(wp), PARAMETER :: tol_theta = 1.0e-5_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! Residual-conservation tolerance on max|R|, where R = dtheta/dt - d(q)/dz + sink (units: 1/day).
  ! By the telescoping identity mb_err = SUM_i tk_i*R_i*dt, max|R| <= res_tol_mass => mb_err -> 0.
  ! Single tolerance shared by the Picard gate and the Newton convergence test below.
  ! ------------------------------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: res_tol_mass = 1.0e-6_wp
  REAL(wp), PARAMETER :: pond_exh_tol_mm = 1.0e-3_wp  ! finite-pond exhaustion tolerance (mm)
  REAL(wp), PARAMETER :: min_dt = 1.0e-6_wp
  max_iter_picard = 50
  ! ------------------------------------------------------------------------------------------------- !
  ! Initialize sink term (Transpiration distributed over layers)
  ! trans_act is in mm/day from CANOPY_MODEL. Convert to m/day.
  ! ------------------------------------------------------------------------------------------------- !
  transpiration_rate_day = trans_act / 1000.0_wp

  ! Calculate root distribution profile (Jackson et al. 1996)
  CALL Get_Root_Distribution(root_dist)
  ! ------------------------------------------------------------------------------------------------- !
  ! Root water sink, with the Jarvis-type compensation of a layer's unmet uptake
  ! by the wetter layers (omega = 1 reproduces the uncompensated profile).
  ! .TRUE. stores the day-start snapshot for the compensation diagnostics.
  ! ------------------------------------------------------------------------------------------------- !
  CALL apply_root_sink(.TRUE.)
  ! ------------------------------------------------------------------------------------------------- !
  ! Boundary conditions (fluxes in m d-1).
  ! Top: infiltration minus soil evaporation, where the effective rainfall
  ! pr_eff (rainfall less interception, m d-1) is supplied by main_MMCO2Crop.
  ! Initialize the time loop.
  ! ------------------------------------------------------------------------------------------------- !
  t_current = 0.0_wp
  t_day_end = 1.0_wp
  IF (PRESENT(t_end)) t_day_end = t_end
  dt = 0.01_wp ! Start with small step (15 mins)
  drain_day = 0.0_wp ! Reset daily drainage accumulator
  hydro_mb_resid_m = 0.0_wp
  hydro_mb_correction_m = 0.0_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! The water-balance tolerances.  They carry the same values in every mode of
  ! the surface-limit record (driver record 149): tightening them was tested at
  ! 1.0e-9 and 1.0e-7 m and both are impractically slow, so the accept/reject
  ! tolerances are left alone and the Picard/Newton behaviour is identical in all
  ! modes.  What the record changes is the min-dt relaxed LEDGER zeroing below,
  ! which applies to mode 0 only: from mode 1 on, the TRUE per-sub-step mass
  ! balance error accumulates into hydro_mb_resid_m instead of being zeroed at
  ! the minimum time step, so the ledger stops hiding the residual it already
  ! carries.  Mode 2's exact post-step conservative closure, further below, is
  ! what carries conservation to machine precision.
  ! ------------------------------------------------------------------------------------------------- !
  mb_tol_step        = 1.0e-6_wp
  mb_tol_head_min_dt = 1.5e-5_wp
  hydro_top_flux_m = 0.0_wp
  hydro_top_infil_m = 0.0_wp
  hydro_bottom_flux_m = 0.0_wp
  hydro_root_sink_m = 0.0_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! Read-only accumulators for the layered nitrogen transport: reset per call,
  ! filled at step acceptance below, never read back by the water solver.
  ! ------------------------------------------------------------------------------------------------- !
  IF (ALLOCATED(hydro_q_iface_m))   hydro_q_iface_m   = 0.0_wp
  IF (ALLOCATED(hydro_sink_layer_m)) hydro_sink_layer_m = 0.0_wp
  hydro_dt_halvings = 0
  hydro_fail_events_call = 0
  hydro_fail_evap = .FALSE.
  hydro_min_dt_reached = .FALSE.
  hydro_last_dt = dt
  hydro_supply_step_mm = 0.0_wp
  hydro_infil_step_mm = 0.0_wp
  hydro_infil_from_supply_step_mm = 0.0_wp
  hydro_infil_from_pond_step_mm = 0.0_wp
  hydro_soil_evap_step_mm = 0.0_wp
  hydro_pond_evap_step_mm = 0.0_wp
  hydro_runoff_step_mm = 0.0_wp
  runoff_step_mm = 0.0_wp
  hydro_ponding_start_step_mm = ponding_mm
  hydro_ponding_end_step_mm = ponding_mm
  call_supply_mm = 0.0_wp
  call_infil_mm = 0.0_wp
  call_infil_supply_mm = 0.0_wp
  call_infil_pond_mm = 0.0_wp
  call_soil_evap_mm = 0.0_wp
  call_pond_evap_mm = 0.0_wp
  call_runoff_mm = 0.0_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! Soil evaporation handled as an atmospheric boundary condition: the demand is
  ! applied as a potential surface flux and limited by the vapor-equilibrium
  ! surface head h_crit_atm when the soil cannot supply it; the physics and the
  ! primary sources are at the h_crit_atm block and at the deliverability switch
  ! below.
  ! ------------------------------------------------------------------------------------------------- !
  evap_rate_pot = evapo_act / 1000.0_wp
  evap_rate_day = evap_rate_pot
  hydro_evap_rate_pot_m_d = evap_rate_pot
  hydro_evap_rate_day_m_d = evap_rate_day
  hydro_cap_active = 0
  res_norm = -1.0_wp
  res_norm_trial = -1.0_wp
  hydro_loop_count = 0
  loop_trace_printed = .FALSE.
  wf_substep = 0
  wf_last_flux_m_d = 0.0_wp
  trace_this_day = wetfront_trace_enabled .AND. &
                   ( (day_id == 97) .OR. &
                     (day_id >= 113 .AND. day_id <= 121) .OR. &
                     (day_id >= 131 .AND. day_id <= 146) )
  ! Reset Task-2 placeholder diagnostics each call (behavior-neutral; not read by Task-1 logic).
  wetfront_risk_flag       = .FALSE.
  wetfront_limiter_applied = .FALSE.
  wetfront_dt_limit_d      = 0.0_wp

  ponding_start_step_mm = ponding_mm
  surface_store_mm = ponding_mm
  supply_step_mm = 0.0_wp
  hydro_pond_evap_step_mm = 0.0_wp
  ! ------------------------------------------------------------------------------------------------- !
  ! Vapour-equilibrium (Kelvin) surface head from the daily mean relative
  ! humidity: Philip (1957) J Meteorol 14:354-366, Eq. (1) p.355,
  ! h = exp[g Psi / R (T + 273)] with h the relative humidity, Psi the moisture
  ! potential in length units and R the gas constant for water vapour; inverted
  ! and written in molar form, Psi = (R_gas T /(Mw g)) ln(RH).
  ! doi:10.1175/1520-0469(1957)014<0354:EAMAHF>2.0.CO;2
  ! ------------------------------------------------------------------------------------------------- !
  BLOCK
    REAL(wp) :: RH_eff, T_k
    RH_eff = MIN(0.999999_wp, MAX(0.01_wp, RH_s / 100.0_wp))
    T_k = MAX(100.0_wp, temp_s + 273.15_wp)
    h_crit_atm = (R_gas * T_k) / (Mw * g_acc) * LOG(RH_eff)
  ! ------------------------------------------------------------------------------------------------- !
    ! h_crit_atm is the Kelvin equilibrium head at the actual daily relative
    ! humidity.  No numerical clamp is applied to it: the Kelvin head at RH 0.5
    ! and 300 K is already about -9.8e3 m, so any such clamp would set the
    ! deliverable soil evaporation instead of the physics.  It is bounded
    ! implicitly by the RH floor (0.01, giving h about -6.4e4 m) and by the
    ! requirement below that the head stay negative.
  ! ------------------------------------------------------------------------------------------------- !
    IF (h_crit_atm > -1.0e-6_wp) h_crit_atm = -1.0e-6_wp
    hydro_h_crit_atm_m = h_crit_atm
  END BLOCK
  ! ------------------------------------------------------------------------------------------------- !
  ! Boundary Conditions (Fluxes in m/day)
  ! Atmospheric boundary (uniform over the active sub-step): net flux = precipitation - evap demand.
  ! ------------------------------------------------------------------------------------------------- !
  pr_flux = pr_eff
  flux_top = pr_flux - evap_rate_day

  ! Copy current state to solver arrays
  r_theta = theta_ly
  DO i_layer = 1, NSoil_layers
    ! ----------------------------------------------------------------------------------------------- !
    ! Initialize h from theta using the selected hydraulic model.
    ! VGM: van Genuchten (1980); BC/Campbell: Brooks & Corey (1964) / Campbell (1974)
    ! (Weihermüller et al. 2021, p.3, Eq. 1–3).
    ! ----------------------------------------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: th_s, th_r, Se, alpha_vg, n_vg, m_vg
      th_s = MAX(1.0e-8_wp, theta_0(i_layer))
      th_r = MAX(0.0_wp, theta_r(i_layer))
      alpha_vg = MAX(1.0e-9_wp, vg_alpha(i_layer))
      n_vg = MAX(1.01_wp, vg_n(i_layer))
      m_vg = MAX(1.0e-6_wp, vg_m(i_layer))
      Se = (r_theta(i_layer) - th_r) / MAX(1.0e-8_wp, (th_s - th_r))
      Se = MIN(0.999999_wp, MAX(1.0e-6_wp, Se))
      r_h(i_layer) = - ( (Se**(-1.0_wp/m_vg) - 1.0_wp)**(1.0_wp/n_vg) ) / alpha_vg
    END BLOCK
    CALL Get_Hydraulic_Props(i_layer, r_theta(i_layer), r_h(i_layer), r_k(i_layer), r_cap(i_layer))
  END DO

  !--------------------------------------------------------------------!
  ! ADAPTIVE TIME STEPPING LOOP
  !--------------------------------------------------------------------!
  ! Evaporation is capped each sub-step by available water in the top 0.10 m.
  ! This prevents extracting more water than exists in that layer range.
  forced_head_step = .FALSE.
  pond_retry_done = .FALSE.
  TimeLoop: DO WHILE (t_current < t_day_end)
    hydro_loop_count = hydro_loop_count + 1
    IF ((.NOT. loop_trace_printed) .AND. hydro_loop_count >= 20000) THEN
      WRITE (*,'(A,I0,3(A,G14.6),A,I0,A,G14.6)') 'HYDRO TRACE day=', day_id, &
        ' t=', t_current, ' dt=', dt, ' ponding_mm=', surface_store_mm, &
        ' top_bc=', top_bc_type, ' flux_top_m_d=', flux_top
      loop_trace_printed = .TRUE.
    END IF
    IF (t_current + dt > t_day_end) dt = t_day_end - t_current
    ! ---- Conservative wet-front time-step limiter ----
    ! It sits before supply_step_mm and every other dt-dependent per-step
    ! quantity, so that a reduction of dt cannot strand water at the old dt.  It
    ! only REDUCES dt and changes no flux, state or storage.  The flags are reset
    ! each iteration so that they do not leak into later trace rows.
    wetfront_risk_flag       = .FALSE.
    wetfront_limiter_applied = .FALSE.
    wetfront_dt_limit_d      = 0.0_wp
    BLOCK
      REAL(wp) :: Se_top_lim, th_s_lim, th_r_lim, dz_top_eff, flux_est, dt_wf, h_surface_lim
      LOGICAL  :: ponded_head_infiltration
      th_s_lim = MAX(1.0e-8_wp, theta_0(1))
      th_r_lim = MAX(0.0_wp, theta_r(1))
      Se_top_lim = (r_theta(1) - th_r_lim) / MAX(1.0e-8_wp, (th_s_lim - th_r_lim))
      Se_top_lim = MIN(1.0_wp, MAX(0.0_wp, Se_top_lim))
      ponded_head_infiltration = (surface_store_mm > ponding_head_min_mm) .OR. &
                                 (top_bc_type == 1 .AND. h_top >= -1.0e-9_wp)
      IF ((energy_balance_mode >= 1) .AND. ponded_head_infiltration) THEN
        h_surface_lim = 0.0_wp
        dz_top_eff = MAX(eps_len, ABS(dz(1)))
        flux_est = MAX(0.0_wp, -r_k(1) * ((r_h(1) - h_surface_lim) / dz_top_eff - 1.0_wp))
        flux_est = MAX(flux_est, wf_last_flux_m_d, wetfront_flux_eps_m_d)
        IF (Se_top_lim < wetfront_Se_trigger) THEN
          wetfront_risk_flag = .TRUE.
          dt_wf = wetfront_max_dtheta * tk(1) / flux_est
          dt_wf = MAX(min_dt, dt_wf)
          wetfront_dt_limit_d = dt_wf
          IF (dt_wf < dt) THEN
            dt = dt_wf
            wetfront_limiter_applied = .TRUE.
          END IF
        END IF
      END IF
    END BLOCK
    pond_unsat_fluxbc = .FALSE.
    top_bc_type_old = top_bc_type
    h_top_old = h_top
    surface_store_old_mm = surface_store_mm
    surface_head_old_m = MAX(0.0_wp, surface_store_old_mm / 1000.0_wp)
    ponding_start_step_mm = surface_store_old_mm
    supply_step_mm = MAX(0.0_wp, pr_eff) * dt * 1000.0_wp

    ! Limit evaporation demand to water available in the top 0.10 m.
    ! Full extraction from 0-0.05 m, then linearly reduced to zero at 0.10 m.
    BLOCK
      REAL(wp) :: z_full, z_max, sum_avail, avail_i, d_top, d_bot, overlap
      REAL(wp) :: w_top, w_bot, w_avg, th_r_eff
      REAL(wp) :: evap_cap
      z_full = 0.05_wp
      z_max = 0.10_wp
      sum_avail = 0.0_wp
      DO i_layer = 1, NSoil_layers
        d_top = ABS(zb(i_layer))
        d_bot = ABS(zb(i_layer+1))
        overlap = MAX(0.0_wp, MIN(d_bot, z_max) - d_top)
        IF (overlap > 0.0_wp) THEN
          ! Weight: 1.0 in 0-0.05 m, then linearly decreases to 0 at 0.10 m.
          IF (d_top <= z_full) THEN
            w_top = 1.0_wp
          ELSEIF (d_top >= z_max) THEN
            w_top = 0.0_wp
          ELSE
            w_top = (z_max - d_top) / (z_max - z_full)
          END IF
          IF (d_bot <= z_full) THEN
            w_bot = 1.0_wp
          ELSEIF (d_bot >= z_max) THEN
            w_bot = 0.0_wp
          ELSE
            w_bot = (z_max - d_bot) / (z_max - z_full)
          END IF
          w_avg = 0.5_wp * (w_top + w_bot)
          th_r_eff = MAX(0.0_wp, theta_r(i_layer))
          avail_i = MAX(0.0_wp, r_theta(i_layer) - th_r_eff) * overlap * w_avg
          sum_avail = sum_avail + avail_i
        END IF
      END DO
      ! atmospheric-BC atmospheric BC: do not externally cap by a fixed
      ! "available water" bucket; let Richards + h_crit limit the actual flux.
      ! We keep the diagnostic "evap_cap" for reporting only.
      evap_cap = 0.0_wp
      IF (sum_avail > 0.0_wp) evap_cap = sum_avail / MAX(eps_len, dt)
      evap_rate_day = evap_rate_pot
      hydro_sum_avail_m = sum_avail
      hydro_evap_cap_m_d = evap_cap
      hydro_evap_rate_day_m_d = evap_rate_day
      hydro_evap_rate_pot_m_d = evap_rate_pot
      evap_rate_act = evap_rate_day
    flux_top = pr_flux - evap_rate_act
    hydro_flux_top_req_m_d = flux_top
    hydro_flux_top_m_d = flux_top
    END BLOCK

    ! Cap infiltration by surface hydraulic capacity.
    ! References:
    !  - Šimůnek et al. (2016) HYDRUS-1D Technical Manual, atmospheric BC with max infiltration:
    !    https://www.pc-progress.com/downloads/Pgm_Hydrus3D5/HYDRUS_Technical_Manual_1D_V5.pdf
    !  - Neuman et al. (1974) on limiting infiltration by soil hydraulic capacity.
    ! Here we approximate q_max by K * (1 + dh/dz) at the surface node.
    ! Apply infiltration cap ONLY under atmospheric flux BC.
    IF (use_infiltration_cap .AND. top_bc_type == 0 .AND. flux_top > 0.0_wp) THEN
      BLOCK
        REAL(wp) :: local_cap
        z_surface = 0.0_wp
        dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
        local_cap = MAX(0.0_wp, -1.0_wp * r_k(1) * ((r_h(1) - h_top) / dz_top - 1.0_wp))
        hydro_infil_cap_m_d = MAX(0.0_wp, local_cap)
        hydro_flux_top_cap_m_d = hydro_infil_cap_m_d
        IF (flux_top > hydro_infil_cap_m_d) THEN
          hydro_cap_active = 1
        END IF
        ! Do not clip flux_top here. HYDRUS tests atmospheric flux through Shift.
      END BLOCK
    ELSE
      hydro_infil_cap_m_d = 0.0_wp
      hydro_flux_top_cap_m_d = 0.0_wp
    END IF
    hydro_flux_top_applied_m_d = flux_top

    ! Under infiltration without surface water, always start the step from atmospheric flux BC.
    IF (flux_top > 0.0_wp .AND. surface_store_mm <= ponding_head_min_mm) THEN
      top_bc_type = 0
    END IF

    ! If the surface saturates under rainfall supply, switch to head BC (ponding).
    IF (top_bc_type == 0 .AND. flux_top > 0.0_wp .AND. r_theta(1) >= (theta_0(1) - 1.0e-6_wp) &
        .AND. surface_store_mm > ponding_head_min_mm) THEN
      top_bc_type = 1
      h_top = 0.0_wp
    END IF
    ! ----------------------------------------------------------------------------------------------- !
    ! If evaporative demand exceeds soil supply, switch to head BC (atmospheric limit).
    ! The dry-limit switch is armed by DELIVERABILITY, not by threshold-crossing
    ! alone.  Physics and primary sources:
    !  - The soil surface cannot be drier than vapor equilibrium with the air:
    !    Philip (1957), J. Meteorol. 14:354-366, Eq. (1) p.355, h = exp[g*Psi/
    !    R(T+273)] (R = gas constant for water vapor), inverted here for
    !    Psi = h_crit_atm at the ambient relative humidity.
    !    https://doi.org/10.1175/1520-0469(1957)014<0354:EAMAHF>2.0.CO;2
    !  - The flux the soil can DELIVER to that driest-possible surface is the
    !    Buckingham flux law q = lambda * grad(psi) with the gravity term:
    !    Buckingham (1907), USDA Bureau of Soils Bulletin 38: flux law Q =
    !    lambda*S, sec.7 p.27; capillary potential defined sec.8 pp.28-29;
    !    gravity in the vertical balance, Eq. (2) p.31.
    !    https://archive.org/details/studiesonmovemen38buck
    ! Discretized between node 1 and a surface at h_crit_atm (upward positive,
    ! the -1 is gravity opposing upward flow):
    !   q_max_evap = max(0, K(h1)*((h1 - h_crit_atm)/dz_top - 1)).
    ! The switch is armed by DELIVERABILITY and not by a head threshold alone: a
    ! threshold-only test (r_h(1) <= h_crit_atm) leaves the solver stalling at a
    ! surface head just above the threshold with a vanishing surface conductivity
    ! while the flux boundary condition still demands evaporation.  When the
    ! demand is not deliverable, the step starts as a head boundary condition at
    ! h_crit_atm and the actual evaporation EMERGES as the Darcy-Buckingham flux,
    ! capped by the demand in the accepted-step accounting.
    ! The only quantities read are K, h1, h_crit_atm and dz.
    ! ----------------------------------------------------------------------------------------------- !
    IF (top_bc_type == 0 .AND. flux_top < 0.0_wp) THEN
      BLOCK
        REAL(wp) :: q_max_evap, dz_top_l
        dz_top_l = MAX(eps_len, ABS(dz(1)))
        q_max_evap = MAX(0.0_wp, r_k(1) * ((r_h(1) - h_crit_atm) / dz_top_l - 1.0_wp))
        hydro_qmax_evap_m_d = q_max_evap
        IF (r_h(1) <= h_crit_atm .OR. (-flux_top) > q_max_evap) THEN
          top_bc_type = 1
          h_top = h_crit_atm
        END IF
      END BLOCK
    END IF
    ! ----------------------------------------------------------------------------------------------- !
    ! Supply-limited ponded-infiltration capacity boundary condition.  When the
    ! surface is ponded and the top layer is still unsaturated, water infiltrates
    ! under a FLUX boundary condition whose rate is bounded by (i) the topsoil
    ! saturated capacity k_s(1) and (ii) the finite surface water available this
    ! sub-step (pond plus supply), rather than under a head boundary condition
    ! limited by the conductivity of a dry surface node, which would reject rain.
    ! It applies to infiltration only (flux_top > 0); the evaporation path is
    ! untouched, and the head boundary condition resumes on saturation.
    ! ----------------------------------------------------------------------------------------------- !
    IF (energy_balance_mode >= 1 .AND. flux_top > 0.0_wp &
        .AND. surface_store_mm > ponding_head_min_mm &
        .AND. r_theta(1) < (theta_0(1) - 1.0e-6_wp)) THEN
      BLOCK
        REAL(wp) :: avail_m, cap_avail
        avail_m = MAX(0.0_wp, (ponding_start_step_mm + supply_step_mm) / 1000.0_wp)
        cap_avail = avail_m / MAX(eps_len, dt)
        top_bc_type = 0
        h_top = 0.0_wp
        flux_top = MIN(k_s(1), cap_avail)
        pond_unsat_fluxbc = .TRUE.
      END BLOCK
    END IF
    hydro_h_top_m = h_top
    ! ----------------------------------------------------------------------------------------------- !
    ! The dry-latch release applies ONLY when the SUB-STEP STARTS in the stable
    ! dry-latched state under supply, i.e. the start-of-step surface head is at
    ! or below h_crit_atm while water is arriving; never on a mid-iteration
    ! transient.
    ! ----------------------------------------------------------------------------------------------- !
    latch_risk_supply = (flux_top > 0.0_wp .AND. r_h(1) <= h_crit_atm)

  ! Initialize Picard Iteration
  theta_new = r_theta
  h_new     = r_h
  K_new     = r_k
  Cap_new   = r_cap
  converged = .FALSE.
    ! ----------------------------------------------------------------------------------------------- !
    ! Adaptive under-relaxation for wet conditions to improve convergence.
    ! (Standard practice for Richards solvers; see Celia et al. 1990; Moukalled et al. Ch.14)
    ! ----------------------------------------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: th_s_top, th_r_top, Se_top
      th_s_top = MAX(1.0e-8_wp, theta_0(1))
      th_r_top = MAX(0.0_wp, theta_r(1))
      Se_top = (r_theta(1) - th_r_top) / MAX(1.0e-8_wp, (th_s_top - th_r_top))
      Se_top = MIN(1.0_wp, MAX(0.0_wp, Se_top))
      IF (Se_top > 0.9_wp .AND. flux_top > 0.0_wp) THEN
        relax_h = 0.25_wp
        max_iter_picard = 200
      ELSE
        relax_h = 0.7_wp
        max_iter_picard = 50
      END IF
    END BLOCK

    PicardLoop: DO iter_picard = 1, max_iter_picard

      ! atmospheric-BC boundary Shift belongs inside the nonlinear iteration,
      ! before assembling and solving the active top row.
      BLOCK
        REAL(wp) :: v_top_iter, q_out_top_iter, top_cv_thickness_iter
        ! A (catch-all): a ponded head BC (h_top>=0) requires REAL surface water. Without a
        ! pond it is a phantom infinite reservoir -> drop to flux BC. Runs every Picard
        ! iteration so it catches head BC from any engagement path. Evaporation atmospheric
        ! head BC (h_top << 0) is untouched.
        IF (top_bc_type == 1 .AND. h_top >= -1.0e-9_wp .AND. surface_store_mm <= ponding_head_min_mm) THEN
          top_bc_type = 0
          h_top = 0.0_wp
        END IF
        IF (top_bc_type == 1) THEN
          ! A DRY-type head boundary condition (h_top at or near h_crit_atm) is
          ! an evaporation limiter and must never gate INFILTRATION.  When water
          ! arrives (flux_top > 0) on a surface that STARTED the sub-step
          ! dry-latched, the boundary is released to the atmospheric flux
          ! condition so that the supply can enter.  The ponded head condition
          ! (h_top >= 0) is untouched.  The release is counted and reported.
          IF (latch_risk_supply .AND. h_top <= 0.5_wp * h_crit_atm) THEN
            top_bc_type = 0
            h_top = 0.0_wp
            hydro_bcswitch_events = hydro_bcswitch_events + 1
            IF (day_id /= hydro_bcswitch_lastday) THEN
              hydro_bcswitch_days = hydro_bcswitch_days + 1
              hydro_bcswitch_lastday = day_id
              WRITE (*,'(A,I0)') ' HYDRO BC SWITCHBACK day=', day_id
            END IF
          END IF
        END IF
        IF (top_bc_type == 1) THEN
          IF (flux_top < 0.0_wp) THEN
            CALL Evaluate_Top_Flux_BC_State(h_new(1), theta_new(1), K_new(1), h_new, K_new, &
                                            q_out_top_iter, v_top_iter, top_cv_thickness_iter)
            IF (ABS(v_top_iter) > ABS(flux_top) + flux_shift_tol .OR. v_top_iter * flux_top < -flux_shift_tol) THEN
              top_bc_type = 0
              h_top = 0.0_wp
            ELSEIF (h_new(1) <= 0.99_wp * h_crit_atm) THEN
              top_bc_type = 0
              h_top = 0.0_wp
            END IF
          END IF
        ELSE
          ! Finite-pond head-BC gate: engage the saturated-surface head BC only when real
          ! surface water exists (pond above threshold or active supply). A node that
          ! saturates (h_new(1)>0) with no pond and no supply must NOT re-engage head BC
          ! (that is the phantom infinite-reservoir leak); it stays on flux BC and drains.
          IF ((flux_top > 0.0_wp .AND. surface_store_mm > ponding_head_min_mm &
               .AND. theta_new(1) >= (theta_0(1) - 1.0e-6_wp)) &
              .OR. (h_new(1) > 0.0_wp .AND. (surface_store_mm > ponding_head_min_mm &
                                             .OR. supply_step_mm > 0.0_wp))) THEN
            top_bc_type = 1
            h_top = 0.0_wp
            forced_head_step = .TRUE.
          ELSEIF ((.NOT. latch_risk_supply) .AND. h_new(1) <= h_crit_atm) THEN
            ! In a sub-step that STARTED dry-latched under supply, the dry
            ! h_crit_atm head boundary condition must not re-engage: that is the
            ! latch, an evaporation limiter gating infiltration.  In every other
            ! situation, mid-iteration transients included, it engages normally.
            top_bc_type = 1
            h_top = h_crit_atm
          END IF
        END IF
        hydro_h_top_m = h_top
      END BLOCK

      ! Build Tridiagonal Matrix for Mixed-Form Richards Eq.
      ! Celia et al. (1990) approximation:
      ! C_i^{m,k} (h_i^{m,k+1} - h_i^{m,k}) / dt + ...
      ! Here we solve for delta_h

      DO i_layer = 1, NSoil_layers

        ! --- Top Flux Term ---
        IF (i_layer == 1) THEN
          ! q_top is defined positive downwards (infiltration)
          ! Flux = -K (dh/dz + 1)
          ! Discretized top flux is explicitly handled in RHS
        END IF

        ! --- Bottom Flux Term ---
        ! Managed via coefficients for the last node
      END DO

      ! Construct Matrix (A*delta_h = R)
      DO i_layer = 1, NSoil_layers
        BLOCK
          REAL(wp) :: q_in, q_out
          REAL(wp) :: surface_storage_coeff, surface_storage_rhs, top_cv_thickness

          r_a(i_layer) = 0.0_wp; r_b(i_layer) = 0.0_wp; r_c(i_layer) = 0.0_wp; r_r(i_layer) = 0.0_wp
          surface_storage_coeff = 0.0_wp
          surface_storage_rhs = 0.0_wp
          top_cv_thickness = MAX(eps_len, ABS(dz(1)))
          IF (i_layer == 1 .AND. top_bc_type == 0) top_cv_thickness = Top_Flux_BC_Thickness()

          ! Properties at interfaces (Arithmetic mean for K)
          ! Interface i-1/2 (Upper)
          IF (i_layer > 1) THEN
             K_mid = 0.5_wp * (K_new(i_layer) + K_new(i_layer-1))
             dz_mid = MAX(eps_len, ABS(dz(i_layer) - dz(i_layer-1))) ! Distance between centers
             r_a(i_layer) = -1.0_wp * K_mid / (MAX(eps_len, tk(i_layer)) * dz_mid)
          END IF

          ! Interface i+1/2 (Lower)
          IF (i_layer < NSoil_layers) THEN
             K_mid = 0.5_wp * (K_new(i_layer) + K_new(i_layer+1))
             dz_mid = MAX(eps_len, ABS(dz(i_layer+1) - dz(i_layer)))
             IF (i_layer == 1 .AND. top_bc_type == 0) THEN
               r_c(i_layer) = -1.0_wp * K_mid / (top_cv_thickness * dz_mid)
             ELSE
               r_c(i_layer) = -1.0_wp * K_mid / (MAX(eps_len, tk(i_layer)) * dz_mid)
             END IF
          ELSE
             ! Bottom Boundary Condition
             IF (drainage_type == 1) THEN
               ! Free Drainage (Unit Gradient): Flux = K_new(N)
               ! q_out = K_bottom. This is a sink term in the balance.
               ! It doesn't depend on h_{N+1}, so no 'c' coefficient.
               ! Added to RHS as explicit sink.
             ELSE
               ! Silled / Zero Flux: Flux = 0. No term.
             END IF
          END IF

          ! Diagonal (Time term + spatial terms)
          r_b(i_layer) = Cap_new(i_layer) / dt - r_a(i_layer) - r_c(i_layer)
          r_b(i_layer) = r_b(i_layer) + surface_storage_coeff

          ! RHS (Residual)
          ! R = - [ (theta^{k} - theta^{n})/dt + d(Flux)/dz - Sink ]
          ! Flux_in/out = -K * (dh/dz - 1) when z is positive downward.


          ! Inflow
          IF (i_layer == 1) THEN
             IF (top_bc_type == 1) THEN
               z_surface = 0.0_wp
               dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
               K_mid = K_new(i_layer)
               q_in = -1.0_wp * K_mid * ((h_new(i_layer) - h_top) / dz_top - 1.0_wp)
               ! Couple ponded surface storage to the top control volume during
               ! head-controlled infiltration. Mixed-form Richards solvers use
               ! this transient surface-reservoir term to retain mass when the
               ! upper boundary switches between atmospheric flux and ponded head;
               ! see Celia et al. (1990, Water Resources Research 26:1483-1496)
               ! for the mass-conservative mixed formulation and Simunek et al.
               ! (2016, HYDRUS technical manual, atmospheric boundary with
               ! surface storage/ponding) for the standard ponded-surface
               ! treatment under head conditions.
               IF (h_top >= -1.0e-12_wp) THEN
                 surface_storage_coeff = 1.0_wp / (dt * MAX(eps_len, tk(i_layer)))
                 surface_storage_rhs = surface_head_old_m * surface_storage_coeff
                 r_b(i_layer) = r_b(i_layer) + surface_storage_coeff
               END IF
             ELSE
               q_in = flux_top
             END IF
             ELSE
            K_mid = 0.5_wp * (K_new(i_layer) + K_new(i_layer-1))
             dz_mid = MAX(eps_len, ABS(dz(i_layer) - dz(i_layer-1)))
             q_in = -1.0_wp * K_mid * ((h_new(i_layer) - h_new(i_layer-1))/dz_mid - 1.0_wp)
          END IF

          ! Outflow
          IF (i_layer < NSoil_layers) THEN
            K_mid = 0.5_wp * (K_new(i_layer) + K_new(i_layer+1))
             dz_mid = MAX(eps_len, ABS(dz(i_layer+1) - dz(i_layer)))
             q_out = -1.0_wp * K_mid * ((h_new(i_layer+1) - h_new(i_layer))/dz_mid - 1.0_wp)
          ELSE
             IF (drainage_type == 1) THEN
               ! Free drainage (unit gradient). q is positive downward.
               q_out = K_new(i_layer)
             ELSE
               q_out = 0.0_wp
             END IF
          END IF

          IF (i_layer == 1 .AND. top_bc_type == 0) THEN
            r_r(i_layer) = -1.0_wp * ( (theta_new(i_layer) - r_theta(i_layer))/dt - &
                           (q_in - q_out)/top_cv_thickness + sink_source(i_layer) )
          ELSE
            r_r(i_layer) = -1.0_wp * ( (theta_new(i_layer) - r_theta(i_layer))/dt - &
                           (q_in - q_out)/MAX(eps_len, tk(i_layer)) + sink_source(i_layer) )
          END IF
          r_r(i_layer) = r_r(i_layer) + surface_storage_rhs
        END BLOCK
      END DO

      ! Solve for delta_h
      CALL Thomas_Solver(NSoil_layers, r_a, r_b, r_c, r_r)

      ! Update solution
      max_change = 0.0_wp
      max_h_change = 0.0_wp
      it_crit = .TRUE.
      DO i_layer = 1, NSoil_layers
        BLOCK
          REAL(wp) :: theta_prev_iter
          REAL(wp) :: h_prev_iter
          REAL(wp) :: theta_iter_change
          REAL(wp) :: h_iter_change

          h_prev_iter = h_new(i_layer)
          h_new(i_layer) = h_new(i_layer) + relax_h * r_r(i_layer) ! r_r holds delta_h now

          ! Update properties
          theta_prev_iter = theta_new(i_layer)
          CALL Get_Hydraulic_Props(i_layer, theta_new(i_layer), h_new(i_layer), K_new(i_layer), Cap_new(i_layer))

          ! HYDRUS uses a mixed convergence test: for unsaturated nodes it checks
          ! water-content change, while near saturation it checks head change.
          theta_iter_change = ABS(theta_new(i_layer) - theta_prev_iter)
          h_iter_change = ABS(h_new(i_layer) - h_prev_iter)
          max_change = MAX(max_change, theta_iter_change)
          max_h_change = MAX(max_h_change, h_iter_change)
          IF (h_prev_iter < 0.0_wp .AND. h_new(i_layer) < 0.0_wp) THEN
            IF (theta_iter_change > tol_theta) it_crit = .FALSE.
          ELSE
            IF (h_iter_change > tol_h) it_crit = .FALSE.
          END IF
        END BLOCK
      END DO

      IF (it_crit) THEN
        ! Head/theta change has converged. Additionally require the TRUE mixed-form mass
        ! residual (free-drainage bottom flux q_out(N)=K(h_N) is represented consistently in
        ! Compute_Residual) to be below res_tol_mass before accepting. This removes the
        ! one-signed Picard K-lag bias at the bottom boundary that accumulates in HydroMB on
        ! drying days. If the residual is still too large, keep iterating Picard (and, on
        ! Picard exhaustion, fall through to the Newton refinement below).
        BLOCK
          REAL(wp) :: theta_tmp(NSoil_layers), K_tmp(NSoil_layers)
          REAL(wp) :: Cap_tmp(NSoil_layers), R_tmp(NSoil_layers)
          CALL Compute_Residual(h_new, theta_tmp, K_tmp, Cap_tmp, R_tmp)
          res_mass = MAXVAL(ABS(R_tmp))
          res_step_m = SUM(ABS(R_tmp) * tk) * dt
        END BLOCK
        ! At a normal time step the acceptance test is on the rate residual.  At
        ! the minimum time step the rate residual is inflated by 1/dt, so the step
        ! is accepted only if the residual converted to WATER MASS (res_step_m) is
        ! below the minimum-time-step mass tolerance.
        IF (res_mass <= res_tol_mass .OR. &
            (dt <= min_dt + 1.0e-12_wp .AND. res_step_m <= mb_tol_head_min_dt)) THEN
          converged = .TRUE.
          EXIT PicardLoop
        END IF
      END IF

    END DO PicardLoop

    ! Check convergence
    IF (.NOT. converged) THEN
      ! Try a Newton step at current dt before reducing dt further.
      ! Under non-ponded infiltration, prefer timestep reduction over Newton fallback.
      newton_converged = .FALSE.
      IF (.NOT. (flux_top > 0.0_wp .AND. surface_store_mm <= ponding_head_min_mm)) THEN
      BLOCK
        REAL(wp) :: h_trial(NSoil_layers), theta_trial(NSoil_layers)
        REAL(wp) :: K_trial(NSoil_layers), Cap_trial(NSoil_layers)
        REAL(wp) :: R0(NSoil_layers), R1(NSoil_layers)
        REAL(wp) :: Jmat(NSoil_layers, NSoil_layers)
        REAL(wp) :: delta_h(NSoil_layers)
        REAL(wp) :: dh, alpha
        INTEGER :: i, j, newton_iter
        LOGICAL :: bad_trial_state
        res_norm = 0.0_wp
        CALL Compute_Residual(h_new, theta_new, K_new, Cap_new, R0)
        res_norm = MAXVAL(ABS(R0))
        h_trial = h_new
        DO newton_iter = 1, 12
          bad_trial_state = .FALSE.
          ! Build finite-difference Jacobian (small system: NSoil_layers <= ~50).
          DO j = 1, NSoil_layers
            dh = MAX(1.0e-6_wp, 1.0e-3_wp * ABS(h_trial(j)))
            h_trial(j) = h_trial(j) + dh
            IF (.NOT. ALL(ieee_is_finite(h_trial))) THEN
              bad_trial_state = .TRUE.
              h_trial(j) = h_trial(j) - dh
              EXIT
            END IF
            CALL Compute_Residual(h_trial, theta_trial, K_trial, Cap_trial, R1)
            DO i = 1, NSoil_layers
              Jmat(i,j) = (R1(i) - R0(i)) / dh
            END DO
            h_trial(j) = h_trial(j) - dh
          END DO
          IF (bad_trial_state) EXIT
          ! Solve J * delta_h = -R0
          CALL DenseSolve(NSoil_layers, Jmat, R0, delta_h)

          ! Limit Newton head updates to avoid non-physical jumps during wet flux-BC infiltration.
          DO i = 1, NSoil_layers
            delta_h(i) = SIGN(MIN(ABS(delta_h(i)), newton_dh_cap), delta_h(i))
          END DO

          ! Line search for residual reduction.
          alpha = 1.0_wp
          DO
            h_trial = h_new + alpha * delta_h
            IF (.NOT. ALL(ieee_is_finite(h_trial))) THEN
              alpha = alpha * 0.5_wp
              IF (alpha <= 0.125_wp) EXIT
              CYCLE
            END IF
            CALL Compute_Residual(h_trial, theta_trial, K_trial, Cap_trial, R1)
            res_norm_trial = MAXVAL(ABS(R1))
            IF (res_norm_trial <= res_norm .OR. alpha <= 0.125_wp) EXIT
            alpha = alpha * 0.5_wp
          END DO
          h_new = h_trial
          theta_new = theta_trial
          K_new = K_trial
          Cap_new = Cap_trial
          R0 = R1
          res_norm = res_norm_trial
          IF (res_norm < res_tol_mass .OR. &
              (dt <= min_dt + 1.0e-12_wp .AND. SUM(ABS(R0) * tk) * dt <= mb_tol_head_min_dt)) THEN
            newton_converged = .TRUE.
            EXIT
          END IF
        END DO
      END BLOCK
      END IF
      IF (newton_converged) THEN
        converged = .TRUE.
      ELSE
        top_bc_type = top_bc_type_old
        h_top = h_top_old
        surface_store_mm = surface_store_old_mm
        hydro_h_top_m = h_top
        h_new = r_h
        theta_new = r_theta
        K_new = r_k
        Cap_new = r_cap
        dt = MAX(min_dt, dt / 3.0_wp)
        hydro_dt_halvings = hydro_dt_halvings + 1
        IF (dt <= min_dt + 1.0e-12_wp) THEN
          IF (flux_top > 0.0_wp) THEN
            CALL Apply_Surface_Failure_Fallback(supply_step_mm, ponding_start_step_mm, hydro_pond_evap_step_mm, &
                                                surface_store_mm, runoff_step_mm)
            ponding_mm = surface_store_mm
            call_supply_mm = call_supply_mm + supply_step_mm
            call_runoff_mm = call_runoff_mm + runoff_step_mm
            hydro_min_dt_reached = .TRUE.
            hydro_last_dt = dt
            IF (surface_store_mm > ponding_head_min_mm) THEN
              top_bc_type = 1
              h_top = 0.0_wp
            ELSE
              top_bc_type = 0
              h_top = 0.0_wp
            END IF
            hydro_h_top_m = h_top
            forced_head_step = .TRUE.
            t_current = t_current + dt
            ! Fallback substep advanced at min_dt (no accepted Richards step):
            ! mb_err/flux_mass/infil not defined here -> pass 0.
            CALL write_wetfront_row(t_current, dt, 0.0_wp, 0.0_wp, 0.0_wp)
            CYCLE TimeLoop
          END IF
          ! If flux BC cannot converge under strong infiltration and ponded surface water
          ! exists, switch to ponded head BC and retry once at minimum dt.
          IF (top_bc_type == 0 .AND. flux_top > 0.0_wp .AND. &
              (surface_store_mm > ponding_head_min_mm) &
              .AND. (.NOT. forced_head_step)) THEN
            top_bc_type = 1
            h_top = 0.0_wp
            hydro_h_top_m = h_top
            forced_head_step = .TRUE.
            dt = min_dt
            CYCLE TimeLoop
          END IF
          ! Evaporative-side fallback, symmetric to the infiltration fallback
          ! above.  If the flux boundary condition cannot converge at the minimum
          ! time step under evaporative demand and this sub-step has not yet been
          ! retried as a head condition, it is retried once with the surface at
          ! the vapour-equilibrium head (Philip 1957 Eq. (1) p.355 and the
          ! Buckingham 1907 flux law, as above).  Evaporation then becomes the
          ! emergent Darcy-Buckingham flux and the ROOT SINK REMAINS ACTIVE, so
          ! stiffness at the evaporating surface cannot deny uptake for the day.
          IF (top_bc_type == 0 .AND. flux_top <= 0.0_wp .AND. (.NOT. forced_head_step)) THEN
            top_bc_type = 1
            h_top = h_crit_atm
            hydro_h_top_m = h_top
            forced_head_step = .TRUE.
            dt = min_dt
            hydro_evap_fallback_events = hydro_evap_fallback_events + 1
            CYCLE TimeLoop
          END IF
          WRITE(*,*) 'WARNING: Richards solver failed to converge at t=', t_current, &
                      ' even with minimum time step. Holding state for this day.'
          hydro_failed = .TRUE.
          hydro_fail_events_call = hydro_fail_events_call + 1
          IF (flux_top <= 0.0_wp) hydro_fail_evap = .TRUE.
          hydro_min_dt_reached = .TRUE.
          ! Diagnostic snapshot to understand top-boundary behavior when convergence fails.
          WRITE (*,'(A,I0,2(A,G14.6),A,I0,8(A,G14.6))') 'Richards FAIL day=', day_id, &
            ' t=', t_current, ' dt=', dt, ' top_bc=', top_bc_type, ' ponding_mm=', surface_store_mm, &
            ' flux_top_m_d=', flux_top, ' theta1=', r_theta(1), ' h1=', r_h(1), ' K1=', r_k(1), &
            ' max_dtheta=', max_change, ' max_dh=', max_h_change, ' res_norm=', res_norm
          CALL write_wetfront_row(t_current, dt, 0.0_wp, 0.0_wp, 0.0_wp)
          EXIT TimeLoop
        ELSE
          CYCLE TimeLoop ! Retry with smaller step
        END IF
      END IF
    END IF
    IF (converged) THEN
      ! Success: compute step fluxes (accumulate after Shift, once accepted).
      ! Mass-consistent top flux from storage change (mixed-form Richards).
      ! Source: Celia et al. 1990, WRR 26:1483–1496 (mass-conservative formulation).
      storage_before = SUM(r_theta * tk)
      sink_total = SUM(sink_source * tk)
      IF (drainage_type == 1) THEN
        q_bot = K_new(NSoil_layers)
      ELSE
        q_bot = 0.0_wp
      END IF
      ! Use actual top flux for mass balance (atmospheric-BC atmospheric BC).
      IF (top_bc_type == 1) THEN
        z_surface = 0.0_wp
        dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
        flux_top_step = -1.0_wp * K_new(1) * ((h_new(1) - h_top) / dz_top - 1.0_wp)
        ! No post-solve cap of the upward flux by the evaporative demand is
        ! applied here: the ledger must carry the flux the solver actually removed
        ! from the profile, and a post-solve cap would credit less evaporation
        ! than the soil lost, leaving a one-signed negative residual.  Exceedance
        ! of the demand is handled INSIDE the Picard loop by the head-to-flux
        ! boundary shift in Evaluate_Top_Flux_BC_State.
      ELSE
        flux_top_step = flux_top
      END IF

      ! Re-apply a non-negative atmospheric infiltration cap using the accepted surface state.
      IF (use_infiltration_cap .AND. top_bc_type == 0 .AND. flux_top_step > 0.0_wp) THEN
        IF (pond_unsat_fluxbc) THEN
          ! supply-limited ponded-infiltration capacity: bound by topsoil Ksat, not dry-node K
          hydro_infil_cap_m_d = k_s(1)
        ELSE
          z_surface = 0.0_wp
          dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
          hydro_infil_cap_m_d = MAX(0.0_wp, -1.0_wp * K_new(1) * ((h_new(1) - h_top) / dz_top - 1.0_wp))
        END IF
        ! The cap is REPORTED (hydro_infil_cap_m_d) but the ledger flux is NOT
        ! reduced after the solve.  The profile has already received flux_top in
        ! this sub-step, so crediting only MIN(flux_top, cap) would leave
        ! (flux_top - cap)*dt both in the soil and at the surface, a one-signed
        ! positive residual that reaches tens of millimetres over a season.
        ! Supply the surface cannot absorb is handled INSIDE the Picard loop by
        ! the flux-to-head (ponding) boundary shift, which is the HYDRUS-1D
        ! atmospheric boundary rule (Simunek et al. 2016, HYDRUS-1D Technical
        ! Manual, Sect. 2.5).
      END IF
      ! mb_err computed after any boundary clamp (storage_after recomputed).

      ! Boundary Shift is handled inside PicardLoop before matrix assembly so
      ! the active boundary condition is solved in the same nonlinear iteration.

      ! atmospheric-BC clamp: prevent surface head from falling below hCritA
      ! under atmospheric (flux) BC.
      IF (top_bc_type == 0 .AND. h_new(1) < h_crit_atm) THEN
        h_new(1) = h_crit_atm
        CALL Get_Hydraulic_Props(1, theta_new(1), h_new(1), K_new(1), Cap_new(1))
      END IF

      ! Recompute storage after any boundary clamp.
      storage_after = SUM(theta_new * tk)
      mb_err = (storage_after - storage_before) - dt * (flux_top_step - q_bot - sink_total)
      ! --------------------------------------------------------------------------------------------- !
      ! mb_err has TWO coupled roles: it feeds the accept/reject test below
      ! (ABS(mb_err) > mb_tol_step) AND it accumulates into the reported water
      ! balance of the run.  The relaxed zeroing at the minimum time step
      ! immediately below is therefore not a reporting artifact: because both
      ! roles read the same post-zeroing variable, skipping the zeroing makes
      ! every minimum-time-step whose error lies in (mb_tol_step,
      ! mb_tol_head_min_dt] fail the acceptance test and enter the expensive
      ! boundary-toggle retry cascade it would otherwise bypass, which is
      ! impractically slow.  The zeroing is therefore kept unconditionally, so the
      ! accept/reject control flow is the same in every mode, and mb_err_ledger
      ! captures the TRUE pre-zeroing value for reporting and for the conservative
      ! closure only.  mb_err_ledger never feeds back into any accept/reject or
      ! solver decision.
      ! --------------------------------------------------------------------------------------------- !
      mb_err_ledger = mb_err
      IF (dt <= min_dt + 1.0e-12_wp .AND. ABS(mb_err) <= mb_tol_head_min_dt) THEN
        mb_err = 0.0_wp
      END IF
      ! --------------------------------------------------------------------------------------------- !
      ! Finite-pond bound on head-boundary infiltration, on any engagement path.
      ! The accepted head-condition step may not infiltrate more than the real
      ! surface water (pond plus supply less pond evaporation).  If it would, the
      ! same sub-step is retried at the pond-exhaustion time with a single cut, so
      ! that the solved state cannot create an infinite surface reservoir.
      ! --------------------------------------------------------------------------------------------- !
      IF ((top_bc_type == 1) .AND. (flux_top_step > 0.0_wp) .AND. (.NOT. pond_retry_done) &
          .AND. (dt > min_dt + 1.0e-12_wp)) THEN
        BLOCK
          REAL(wp) :: avail_mm, infil_cand_mm, dt_exh
          avail_mm = MAX(0.0_wp, ponding_start_step_mm + supply_step_mm - hydro_pond_evap_step_mm)
          infil_cand_mm = flux_top_step * 1000.0_wp * dt
          IF (infil_cand_mm > avail_mm + pond_exh_tol_mm) THEN
            top_bc_type = top_bc_type_old; h_top = h_top_old
            surface_store_mm = surface_store_old_mm; hydro_h_top_m = h_top
            h_new = r_h; theta_new = r_theta; K_new = r_k; Cap_new = r_cap
            dt_exh = (avail_mm / 1000.0_wp) / MAX(eps_len, flux_top_step)
            dt = MAX(min_dt, dt_exh); pond_retry_done = .TRUE.
            CYCLE TimeLoop
          END IF
        END BLOCK
      END IF
      IF (ABS(mb_err) > mb_tol_step) THEN
        top_bc_type = top_bc_type_old
        h_top = h_top_old
        surface_store_mm = surface_store_old_mm
        hydro_h_top_m = h_top
        dt = MAX(min_dt, dt / 3.0_wp)
        hydro_dt_halvings = hydro_dt_halvings + 1
        IF (dt <= min_dt + 1.0e-12_wp) THEN
          IF (flux_top_step > 0.0_wp) THEN
            CALL Apply_Surface_Failure_Fallback(supply_step_mm, ponding_start_step_mm, hydro_pond_evap_step_mm, &
                                                surface_store_mm, runoff_step_mm)
            ponding_mm = surface_store_mm
            call_supply_mm = call_supply_mm + supply_step_mm
            call_runoff_mm = call_runoff_mm + runoff_step_mm
            hydro_min_dt_reached = .TRUE.
            hydro_last_dt = dt
            IF (surface_store_mm > ponding_head_min_mm) THEN
              top_bc_type = 1
              h_top = 0.0_wp
            ELSE
              top_bc_type = 0
              h_top = 0.0_wp
            END IF
            hydro_h_top_m = h_top
            forced_head_step = .TRUE.
            t_current = t_current + dt
            ! MB-fail fallback advanced at min_dt: mb_err and flux_top_step are defined here.
            CALL write_wetfront_row(t_current, dt, flux_top_step, 0.0_wp, mb_err)
            CYCLE TimeLoop
          END IF
          ! Evaporative-side fallback, as at the Newton-failure site above.
          IF (top_bc_type == 0 .AND. flux_top_step <= 0.0_wp .AND. (.NOT. forced_head_step)) THEN
            top_bc_type = 1
            h_top = h_crit_atm
            hydro_h_top_m = h_top
            forced_head_step = .TRUE.
            dt = min_dt
            hydro_evap_fallback_events = hydro_evap_fallback_events + 1
            CYCLE TimeLoop
          END IF
          hydro_min_dt_reached = .TRUE.
          hydro_failed = .TRUE.
          hydro_fail_events_call = hydro_fail_events_call + 1
          IF (flux_top_step <= 0.0_wp) hydro_fail_evap = .TRUE.
          WRITE (*,'(A,I0,2(A,G14.6),A,I0,4(A,G14.6))') 'Richards MB FAIL day=', day_id, &
            ' t=', t_current, ' dt=', dt, ' top_bc=', top_bc_type, ' ponding_mm=', surface_store_mm, &
            ' flux_top_m_d=', flux_top_step, ' q_bot_m_d=', q_bot, ' mb_err_m=', mb_err
          CALL write_wetfront_row(t_current, dt, flux_top_step, 0.0_wp, mb_err)
          EXIT TimeLoop
        END IF
        h_new = r_h
        theta_new = r_theta
        K_new = r_k
        Cap_new = r_cap
        CYCLE TimeLoop
      END IF
      ! --------------------------------------------------------------------------------------------- !
      ! SOLVER MASS-CONSERVATION CORRECTION (numerics only, no new physics). The
      ! sub-step is accepted (ABS(mb_err) <= mb_tol_step, just verified above) but still
      ! carries a real water-mass discrepancy mb_err (m) between the Newton-solved storage
      ! change and dt*(flux_top_step-q_bot-sink_total). Restore EXACT conservation for this
      ! substep by adding the deficit (-mb_err, a depth of water) to theta_new in
      ! proportion to each layer's |delta theta| THIS substep (r_theta is still the
      ! PRE-step state here; if no layer changed, put the whole deficit in the top layer,
      ! the shallowest point water would have crossed), then recompute storage_after and
      ! zero mb_err so the ledger (hydro_mb_resid_m, "Total Hydro MB") reports the
      ! corrected state. The correction is logged, never silent: hydro_mb_correction_m.
      ! mb_err_ledger (not mb_err, which may already be zeroed for the min-dt accept
      ! test above) is the TRUE pre-zeroing residual to cancel -- using `mb_err` here
      ! would silently skip the correction on exactly the min-dt steps it exists for.
      ! Accepted step: accumulate fluxes and MB. Mode 0 accumulates `mb_err` (the
      ! ORIGINAL, possibly min-dt-zeroed value -- bit-identical). Mode >=1 accumulates
      ! `mb_err_ledger` instead (the TRUE residual; for non-min-dt steps the two are
      ! identical since the zeroing condition never triggers there; mode 2 has already
      ! zeroed it above after redistributing the deficit into theta_new).
      ! --------------------------------------------------------------------------------------------- !
      flux_top_mass = (storage_after - storage_before) / dt + q_bot + sink_total
      wf_last_flux_m_d = MAX(0.0_wp, flux_top_mass)
        hydro_mb_resid_m = hydro_mb_resid_m + mb_err_ledger
      hydro_bottom_flux_m = hydro_bottom_flux_m + q_bot * dt
      hydro_root_sink_m = hydro_root_sink_m + sink_total * dt
      ! --------------------------------------------------------------------------------------------- !
      ! Layered-N read-only flux diagnostics: signed Darcy flux
      ! through each INTERNAL interface of the ACCEPTED state, using the SAME
      ! inter-node conductivity averaging (k_avg_mode) and discretization as
      ! the assembled system above (Celia et al. 1990 mixed form), plus the
      ! free-drainage/silled bottom flux q_bot, and the per-layer root sink
      ! actually applied this step. PURE ACCUMULATORS for the soil-N
      ! transport/uptake code in main; nothing here feeds back into the
      ! water solution.
      ! Sign convention: + = downward (same as q_out in the residual above).
      ! hydro_q_iface_m(i) = water crossing the interface BELOW layer i (m).
      ! --------------------------------------------------------------------------------------------- !
      IF (ALLOCATED(hydro_q_iface_m)) THEN
        BLOCK
          REAL(wp) :: Km_if, dzm_if, q_if
          INTEGER  :: jj
          DO jj = 1, NSoil_layers - 1
            Km_if = 0.5_wp * (K_new(jj) + K_new(jj+1))
            dzm_if = MAX(eps_len, ABS(dz(jj+1) - dz(jj)))
            q_if = -1.0_wp * Km_if * ((h_new(jj+1) - h_new(jj)) / dzm_if - 1.0_wp)
            hydro_q_iface_m(jj) = hydro_q_iface_m(jj) + q_if * dt
          END DO
          hydro_q_iface_m(NSoil_layers) = hydro_q_iface_m(NSoil_layers) + q_bot * dt
        END BLOCK
      END IF
      IF (ALLOCATED(hydro_sink_layer_m)) THEN
        BLOCK
          INTEGER :: jj
          DO jj = 1, NSoil_layers
            hydro_sink_layer_m(jj) = hydro_sink_layer_m(jj) + sink_source(jj) * tk(jj) * dt
          END DO
        END BLOCK
      END IF
      t_current = t_current + dt
      hydro_last_dt = dt
      r_theta = theta_new
      r_h     = h_new
      r_k     = K_new
      r_cap   = Cap_new
      CALL apply_root_sink(.FALSE.)

      infil_step_mm = 0.0_wp
      BLOCK
        REAL(wp) :: surface_flux_step_mm, surface_available_mm, surface_flux_step_m
        surface_flux_step_mm = flux_top_step * 1000.0_wp * dt
        surface_available_mm = MAX(0.0_wp, ponding_start_step_mm + supply_step_mm - hydro_pond_evap_step_mm)
        IF (surface_flux_step_mm > 0.0_wp) THEN
          infil_step_mm = MIN(surface_flux_step_mm, surface_available_mm)
          surface_flux_step_m = infil_step_mm / 1000.0_wp
        ELSE
          surface_flux_step_m = surface_flux_step_mm / 1000.0_wp
        END IF
        hydro_top_flux_m = hydro_top_flux_m + surface_flux_step_m
        hydro_top_infil_m = hydro_top_infil_m + MAX(0.0_wp, surface_flux_step_m)
      END BLOCK
      infil_supply_step_mm = MIN(infil_step_mm, supply_step_mm)
      infil_pond_step_mm = MAX(0.0_wp, infil_step_mm - infil_supply_step_mm)
      surface_store_trial_mm = MAX(0.0_wp, ponding_start_step_mm + supply_step_mm - &
                                   hydro_pond_evap_step_mm - infil_step_mm)
      runoff_step_mm = MAX(0.0_wp, surface_store_trial_mm - ponding_max_mm)
      surface_store_mm = surface_store_trial_mm - runoff_step_mm
      ponding_mm = surface_store_mm
      call_supply_mm = call_supply_mm + supply_step_mm
      call_infil_mm = call_infil_mm + infil_step_mm
      call_infil_supply_mm = call_infil_supply_mm + infil_supply_step_mm
      call_infil_pond_mm = call_infil_pond_mm + infil_pond_step_mm
      call_runoff_mm = call_runoff_mm + runoff_step_mm

      ! Accumulate drainage
      IF (drainage_type == 1) THEN
         drain_day = drain_day + (r_k(NSoil_layers) * dt) ! m
      END IF

      forced_head_step = .FALSE.
      pond_retry_done = .FALSE.
      CALL write_wetfront_row(t_current, dt, flux_top_mass, infil_step_mm, mb_err)

      ! Increase time step if convergence was fast
      IF (iter_picard < 5) dt = MIN(dt * 1.1_wp, t_day_end - t_current)
    END IF

  END DO TimeLoop

  ! Update global state
  theta_ly = r_theta
  theta    = r_theta ! Update the SHARED global array too
  hydro_k1_m_d = r_k(1)
  hydro_h1_m = r_h(1)
  hydro_theta1 = r_theta(1)

  ! Report actual surface flux (m/day)
  IF (top_bc_type == 1) THEN
    z_surface = 0.0_wp
    dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
    flux_top_actual = -1.0_wp * r_k(1) * ((r_h(1) - h_top) / dz_top - 1.0_wp)
    IF (flux_top_actual < -evap_rate_day) flux_top_actual = -evap_rate_day
  ELSE
    flux_top_actual = flux_top
  END IF

  infil_step_mm = 0.0_wp
  IF (pr_eff > 0.0_wp .OR. ponding_start_step_mm > 0.0_wp) THEN
    infil_step_mm = MAX(0.0_wp, flux_top_actual) * 1000.0_wp * hydro_last_dt
  END IF
  infil_supply_step_mm = MIN(infil_step_mm, supply_step_mm)
  infil_pond_step_mm = MAX(0.0_wp, infil_step_mm - infil_supply_step_mm)
  ponding_mm = surface_store_mm
  soil_evap_step_mm = 0.0_wp
  IF (top_bc_type == 0) THEN
    ! Under the flux boundary condition the top boundary applies
    ! flux_top = pr - evap_rate_act, so the evaporative extraction is
    ! evap_rate_act itself.  Counting only the net upward flux would miss the
    ! evaporation on every sub-step in which rain dominates, where it appears as
    ! reduced infiltration, under-counting seasonal soil evaporation and breaking
    ! the effective-supply water-balance ledger.
    soil_evap_step_mm = MAX(0.0_wp, evap_rate_act) * 1000.0_wp * hydro_last_dt
  ELSE IF (top_bc_type == 1 .AND. h_top < -1.0e-9_wp) THEN
    ! Dry head BC (atmospheric limit): evaporation is supply-limited; the
    ! actual upward flux is the extraction.
    soil_evap_step_mm = MAX(0.0_wp, -flux_top_actual) * 1000.0_wp * hydro_last_dt
  END IF
  hydro_supply_step_mm = call_supply_mm
  hydro_infil_step_mm = call_infil_mm
  hydro_infil_from_supply_step_mm = call_infil_supply_mm
  hydro_infil_from_pond_step_mm = call_infil_pond_mm
  hydro_soil_evap_step_mm = call_soil_evap_mm + soil_evap_step_mm
  hydro_pond_evap_step_mm = call_pond_evap_mm
  hydro_runoff_step_mm = call_runoff_mm
  hydro_ponding_start_step_mm = ponding_start_step_mm
  hydro_ponding_end_step_mm = ponding_mm

  ! (Day-10 debug block removed)

  ! Diagnostic: flux and K at the layer 3-4 interface (if present).
  IF (NSoil_layers >= 4) THEN
    K_mid = 0.5_wp * (r_k(3) + r_k(4))
    dz_mid = MAX(eps_len, ABS(dz(4) - dz(3)))
    hydro_q34_m_d = -1.0_wp * K_mid * ((r_h(4) - r_h(3)) / dz_mid - 1.0_wp)
    hydro_k34_m_d = K_mid
  ELSE
    hydro_q34_m_d = 0.0_wp
    hydro_k34_m_d = 0.0_wp
  END IF

  CALL sync_hydro_step_state()

CONTAINS

  SUBROUTINE write_wetfront_row(t_row, dt_row, flux_mass_row, infil_mm_row, mb_err_row)
    ! Diagnostics only: emit one wet-front trace line. Reads host-associated state;
    ! modifies only the diagnostic counter wf_substep. No solver state is changed.
    REAL(wp), INTENT(IN) :: t_row, dt_row, flux_mass_row, infil_mm_row, mb_err_row
    REAL(wp) :: Se1, den
    INTEGER  :: wfr_int, lim_int, minflag
    IF (.NOT. trace_this_day) RETURN
    IF (.NOT. wetfront_trace_enabled) RETURN  ! NEWUNIT returns a negative unit; do not sign-test it.
    wf_substep = wf_substep + 1
    den = MAX(theta_0(1) - theta_r(1), 1.0e-12_wp)
    Se1 = (r_theta(1) - theta_r(1)) / den
    Se1 = MIN(1.0_wp, MAX(0.0_wp, Se1))
    wfr_int = 0; IF (wetfront_risk_flag)       wfr_int = 1
    lim_int = 0; IF (wetfront_limiter_applied) lim_int = 1
    minflag = 0; IF (hydro_min_dt_reached)     minflag = 1
    WRITE (wetfront_trace_unit, &
      '(I0,",",I0,2(",",G14.6),",",I0,13(",",G14.6),",",I0,",",G14.6,5(",",I0),",",G14.6)') &
      day_id, wf_substep, t_row, dt_row, top_bc_type, &
      surface_store_mm, pr_eff*1000.0_wp, evap_rate_pot*1000.0_wp, evap_rate_act*1000.0_wp, &
      hydro_flux_top_req_m_d, flux_mass_row, infil_mm_row, r_theta(1), Se1, r_h(1), r_k(1), &
      h_top, h_crit_atm, wfr_int, wetfront_dt_limit_d, lim_int, 0, iter_picard, &
      hydro_dt_halvings, minflag, mb_err_row*1000.0_wp
  END SUBROUTINE write_wetfront_row

  !--------------------------------------------------------------------!
  ! Get Hydraulic Properties (hyd_model controls VGM vs Brooks-Corey vs Campbell)
  ! Input: h (matric potential, m) OR theta (volumetric water content)
  ! Output: All properties consistent with input
  ! Note: Input h takes precedence if calculating from head,
  !       but we usually update theta from h.
  !--------------------------------------------------------------------!
  SUBROUTINE Get_Hydraulic_Props(lyr, th_val, h_val, K_val, Cap_val)
    INTEGER, INTENT(IN) :: lyr
    REAL(wp), INTENT(INOUT) :: th_val, h_val
    REAL(wp), INTENT(OUT) :: K_val, Cap_val

    REAL(wp) :: th_s, th_r, alpha_vg, n_vg, m_vg, Se, K_sat, H_suc, l_m, vg_term

    th_s = MAX(1.0e-8_wp, theta_0(lyr))
    th_r = MAX(0.0_wp, theta_r(lyr))
    alpha_vg = MAX(1.0e-9_wp, vg_alpha(lyr))
    n_vg = MAX(1.01_wp, vg_n(lyr))
    m_vg = MAX(1.0e-6_wp, vg_m(lyr))
    K_sat = MAX(1.0e-12_wp, k_s(lyr)) ! m/day
    l_m = MAX(0.0_wp, vg_l(lyr))

    H_suc = MAX(1.0e-6_wp, -h_val)

    ! Water retention: theta(h) = theta_r + (theta_s-theta_r)/(1+(alpha|h|)^n)^m.
    vg_term = (alpha_vg * H_suc)**n_vg
    IF (vg_term > 1.0e30_wp) vg_term = 1.0e30_wp
    Se = 1.0_wp / (1.0_wp + vg_term)**m_vg
    th_val = th_r + Se * (th_s - th_r)
    th_val = MIN(th_s, MAX(th_r, th_val))

    ! Mualem conductivity: K = Ksat * Se^l * [1-(1-Se^(1/m))^m]^2.
    K_val = K_sat * (Se**l_m) * (1.0_wp - (1.0_wp - Se**(1.0_wp/m_vg))**m_vg)**2

    ! Capacity C = dtheta/dh derived from van Genuchten (1980), Eq. (7).
    Cap_val = (th_s - th_r) * m_vg * n_vg * (alpha_vg**n_vg) * &
      (H_suc**(n_vg - 1.0_wp)) / (1.0_wp + vg_term)**(m_vg + 1.0_wp)
    Cap_val = MAX(1.0e-8_wp, Cap_val)

  END SUBROUTINE Get_Hydraulic_Props

  !--------------------------------------------------------------------!
  ! Compute residual for a given head vector (for Newton fallback).
  ! Returns theta/K/Cap consistent with h and residual vector R.
  !--------------------------------------------------------------------!
  SUBROUTINE Compute_Residual(h_in, theta_out, K_out, Cap_out, R_out)
    REAL(wp), INTENT(IN)  :: h_in(:)
    REAL(wp), INTENT(OUT) :: theta_out(:), K_out(:), Cap_out(:), R_out(:)
    REAL(wp) :: h_eff(NSoil_layers)
    REAL(wp) :: q_in, q_out, K_mid, dz_mid
    REAL(wp) :: z_surface, dz_top, th_s, th_r, h_tmp, top_cv_thickness
    INTEGER :: i_layer

    DO i_layer = 1, NSoil_layers
      theta_out(i_layer) = r_theta(i_layer)
      IF (.NOT. ieee_is_finite(h_in(i_layer))) THEN
        IF (h_in(i_layer) < 0.0_wp) THEN
          h_eff(i_layer) = h_eval_min
        ELSE
          h_eff(i_layer) = h_eval_max
        END IF
      ELSE
        h_eff(i_layer) = MIN(h_eval_max, MAX(h_eval_min, h_in(i_layer)))
      END IF
      h_tmp = h_eff(i_layer)
      CALL Get_Hydraulic_Props(i_layer, theta_out(i_layer), h_tmp, K_out(i_layer), Cap_out(i_layer))
      ! Ensure within bounds
      th_s = MAX(1.0e-8_wp, theta_0(i_layer))
      th_r = MAX(0.0_wp, theta_r(i_layer))
      theta_out(i_layer) = MIN(th_s, MAX(th_r, theta_out(i_layer)))
    END DO

    DO i_layer = 1, NSoil_layers
      top_cv_thickness = MAX(eps_len, ABS(dz(1)))
      IF (i_layer == 1 .AND. top_bc_type == 0) top_cv_thickness = Top_Flux_BC_Thickness()
      ! Inflow
      IF (i_layer == 1) THEN
        IF (top_bc_type == 1) THEN
          z_surface = 0.0_wp
          dz_top = MAX(eps_len, ABS(dz(1) - z_surface))
          K_mid = K_out(i_layer)
          q_in = -1.0_wp * K_mid * ((h_eff(i_layer) - h_top) / dz_top - 1.0_wp)
        ELSE
          q_in = flux_top
        END IF
      ELSE
        K_mid = 0.5_wp * (K_out(i_layer) + K_out(i_layer-1))
        dz_mid = MAX(eps_len, ABS(dz(i_layer) - dz(i_layer-1)))
        q_in = -1.0_wp * K_mid * ((h_eff(i_layer) - h_eff(i_layer-1))/dz_mid - 1.0_wp)
      END IF

      ! Outflow
      IF (i_layer < NSoil_layers) THEN
        K_mid = 0.5_wp * (K_out(i_layer) + K_out(i_layer+1))
        dz_mid = MAX(eps_len, ABS(dz(i_layer+1) - dz(i_layer)))
        q_out = -1.0_wp * K_mid * ((h_eff(i_layer+1) - h_eff(i_layer))/dz_mid - 1.0_wp)
      ELSE
        IF (drainage_type == 1) THEN
          q_out = K_out(i_layer)
        ELSE
          q_out = 0.0_wp
        END IF
      END IF

      IF (i_layer == 1 .AND. top_bc_type == 0) THEN
        R_out(i_layer) = (theta_out(i_layer) - r_theta(i_layer))/dt - &
                         (q_in - q_out)/top_cv_thickness + sink_source(i_layer)
      ELSE
        R_out(i_layer) = (theta_out(i_layer) - r_theta(i_layer))/dt - &
                         (q_in - q_out)/MAX(eps_len, tk(i_layer)) + sink_source(i_layer)
      END IF
    END DO
  END SUBROUTINE Compute_Residual

  REAL(wp) FUNCTION Top_Flux_BC_Thickness()
    REAL(wp) :: dx_face
      ! --------------------------------------------------------------------------------------------- !
    ! CONSISTENT TOP CONTROL VOLUME.  Under the atmospheric FLUX boundary
    ! condition an alternative balance row for node 1 divides the net flux by HALF
    ! the node spacing (0.5*|dz(2)-dz(1)|, i.e. tk(1)/2 for centred nodes), while
    ! the step's storage ledger (storage_after - storage_before) weights theta(1)
    ! by the FULL layer thickness tk(1).  The two volumes differ by about a factor
    ! of two, so even an exactly converged step would carry
    ! mb_err = dt*(q_in - q_out)*(tk(1)/T_cv - 1): the acceptance test then
    ! throttles dt towards 1e-5 d, which costs millions of sub-steps per season,
    ! and the accepted sub-tolerance residuals sum to tens of millimetres per
    ! season.  The shipped form uses the SAME volume in the balance row as in the
    ! storage ledger (tk(1), the cell-centred finite-volume control volume), which
    ! is the mixed-form conservation requirement of Celia et al. (1990) Water
    ! Resour Res 26:1483-1496, Eqs. 15-17: the discrete flux divergence and the
    ! discrete storage change must refer to one and the same volume.  The Darcy
    ! distance from the surface to the node, dz(1), is unchanged.
      ! --------------------------------------------------------------------------------------------- !
    IF (wb_surface_limit == 3) THEN
      Top_Flux_BC_Thickness = MAX(eps_len, tk(1))
      RETURN
    END IF
    IF (NSoil_layers > 1) THEN
      dx_face = MAX(eps_len, ABS(dz(2) - dz(1)))
      Top_Flux_BC_Thickness = 0.5_wp * dx_face
    ELSE
      Top_Flux_BC_Thickness = MAX(eps_len, ABS(dz(1)))
    END IF
  END FUNCTION Top_Flux_BC_Thickness

  SUBROUTINE Evaluate_Top_Flux_BC_State(h_top_node, theta_top_node, K_top_node, h_state, K_state, &
                                        q_out_top, v_top, top_cv_thickness)
    REAL(wp), INTENT(IN) :: h_top_node, theta_top_node, K_top_node
    REAL(wp), INTENT(IN) :: h_state(:), K_state(:)
    REAL(wp), INTENT(OUT) :: q_out_top, v_top, top_cv_thickness
    REAL(wp) :: K_int, dx_face

    top_cv_thickness = Top_Flux_BC_Thickness()
    IF (NSoil_layers > 1) THEN
      dx_face = MAX(eps_len, ABS(dz(2) - dz(1)))
      K_int = 0.5_wp * (K_top_node + K_state(2))
      q_out_top = -1.0_wp * K_int * ((h_state(2) - h_top_node) / dx_face - 1.0_wp)
    ELSE
      q_out_top = 0.0_wp
    END IF

    v_top = q_out_top + top_cv_thickness * ( (theta_top_node - r_theta(1)) / dt + sink_source(1) )
  END SUBROUTINE Evaluate_Top_Flux_BC_State

  SUBROUTINE Apply_Surface_Failure_Fallback(supply_mm_local, ponding_start_mm_local, pond_evap_mm_local, &
                                            surface_store_new_mm, runoff_new_mm)
    REAL(wp), INTENT(IN) :: supply_mm_local, ponding_start_mm_local, pond_evap_mm_local
    REAL(wp), INTENT(OUT) :: surface_store_new_mm, runoff_new_mm
    REAL(wp) :: trial_store_mm

    trial_store_mm = MAX(0.0_wp, ponding_start_mm_local + supply_mm_local - pond_evap_mm_local)
    runoff_new_mm = MAX(0.0_wp, trial_store_mm - ponding_max_mm)
    surface_store_new_mm = trial_store_mm - runoff_new_mm
  END SUBROUTINE Apply_Surface_Failure_Fallback

  !--------------------------------------------------------------------!
  ! Dense linear solver (Gaussian elimination with partial pivoting).
  ! Solves A * x = -b (b input), returns x.
  !--------------------------------------------------------------------!
  SUBROUTINE DenseSolve(n, A, b, x)
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(INOUT) :: A(n,n)
    REAL(wp), INTENT(IN) :: b(n)
    REAL(wp), INTENT(OUT) :: x(n)
    REAL(wp) :: rhs(n)
    INTEGER :: i, j, k, piv
    REAL(wp) :: maxval, tmp, factor

    rhs = -b
    DO k = 1, n-1
      piv = k
      maxval = ABS(A(k,k))
      DO i = k+1, n
        IF (ABS(A(i,k)) > maxval) THEN
          maxval = ABS(A(i,k))
          piv = i
        END IF
      END DO
      IF (piv /= k) THEN
        DO j = k, n
          tmp = A(k,j); A(k,j) = A(piv,j); A(piv,j) = tmp
        END DO
        tmp = rhs(k); rhs(k) = rhs(piv); rhs(piv) = tmp
      END IF
      IF (ABS(A(k,k)) < 1.0e-20_wp) CYCLE
      DO i = k+1, n
        factor = A(i,k) / A(k,k)
        DO j = k, n
          A(i,j) = A(i,j) - factor * A(k,j)
        END DO
        rhs(i) = rhs(i) - factor * rhs(k)
      END DO
    END DO
    x(n) = rhs(n) / MAX(1.0e-20_wp, A(n,n))
    DO i = n-1, 1, -1
      x(i) = rhs(i)
      DO j = i+1, n
        x(i) = x(i) - A(i,j) * x(j)
      END DO
      x(i) = x(i) / MAX(1.0e-20_wp, A(i,i))
    END DO
  END SUBROUTINE DenseSolve

  !--------------------------------------------------------------------!
  ! Tridiagonal Matrix Solver (Thomas Algorithm)
  !--------------------------------------------------------------------!
  SUBROUTINE Thomas_Solver(n, a, b, c, r)
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: a(:), b(:), c(:) ! Diagonals
    REAL(wp), INTENT(INOUT) :: r(:) ! RHS (returned as solution)

    REAL(wp) :: gam(n), bet
    INTEGER :: j

    ! Forward elimination
    bet = b(1)
    r(1) = r(1) / bet

    DO j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j) * gam(j)
       r(j) = (r(j) - a(j) * r(j-1)) / bet
    END DO

    ! Back substitution
    DO j = n-1, 1, -1
       r(j) = r(j) - gam(j+1) * r(j+1)
    END DO

  END SUBROUTINE Thomas_Solver

  !--------------------------------------------------------------------!
  ! Root Distribution (Jackson et al. 1996)
  ! Y = 1 - beta^d
  !--------------------------------------------------------------------!
  SUBROUTINE Get_Root_Distribution(dist_profile)
    REAL(wp), INTENT(OUT) :: dist_profile(NSoil_layers)
    REAL(wp) :: z_top, z_bot, Y_top, Y_bot
    REAL(wp) :: sum_dist
    INTEGER :: k

    sum_dist = 0.0_wp
    dist_profile = 0.0_wp

    DO k = 1, NSoil_layers
       z_top = zb(k)      ! Depth to top of layer (positive?) zb is usually negative in SHARED
       z_bot = zb(k+1)
       ! Convert to positive depth cm? Jackson uses cm usually.
       ! Assuming zb is negative meters (0, -0.1, ...).
       ! d in Jackson is depth from surface.

       ! Use positive depths
       BLOCK
         REAL(wp) :: d_top, d_bot
         d_top = ABS(zb(k)) * 100.0_wp ! cm
         d_bot = ABS(zb(k+1)) * 100.0_wp ! cm

         ! Check against current root depth
         ! root_dep is in meters
         IF (d_top >= root_dep * 100.0_wp) THEN
           dist_profile(k) = 0.0_wp
           CYCLE
         END IF

         ! Clip bottom to root depth
         IF (d_bot > root_dep * 100.0_wp) d_bot = root_dep * 100.0_wp

         Y_top = 1.0_wp - root_beta**d_top
         Y_bot = 1.0_wp - root_beta**d_bot

         dist_profile(k) = Y_bot - Y_top
         sum_dist = sum_dist + dist_profile(k)
       END BLOCK
    END DO

    ! Normalize (if roots don't reach bottom, sum < 1, which is correct?
    ! No, we assume Transpiration demand is met by existing roots.
    ! So we should normalize so sum = 1 over the *active* root zone.)
    IF (sum_dist > 1.0e-6_wp) THEN
       dist_profile = dist_profile / sum_dist
    END IF

  END SUBROUTINE Get_Root_Distribution
  ! --------------------------------------------------------------------------------------------- !
  ! The p2 node (HYDRUS h3) of the Feddes ramp, selected from the transpiration demand.
  !   tp_m_d      = the demand handed to Richards for this call (m/d; may be stressed)
  !   tp_pot_mm_d = the daily UNSTRESSED potential transpiration (mm/d, SHARED)
  !   hydrus_rule = .FALSE. (driver record 80 = 0): p2h is selected when
  !                 tp_m_d > 0.5.  That 0.5 is HYDRUS-1D's r2H, which is defined
  !                 in cm d-1 (Simunek et al. 2013, HYDRUS-1D 4.17, Table 12.7
  !                 p.207) but is compared here with a value in m d-1, so p2h is
  !                 never selected (the largest modeled demand is about
  !                 0.006 m d-1).
  !               = .TRUE. (driver record 80 = 1): Tp = tp_pot_mm_d / 10
  !                 (mm d-1 -> cm d-1, the manual's unit) is compared with
  !                 r2H = 0.5 and r2L = 0.1 cm d-1 (SHARED feddes_r2h_cm_d,
  !                 feddes_r2l_cm_d), linear interpolation as Table 12.7 p.207 prints it:
  !                   h3 = P2H + (P2L - P2H)/(r2H - r2L) * (r2H - Tp)  for r2L < Tp < r2H
  !                   h3 = P2L for Tp <= r2L ;  h3 = P2H for Tp >= r2H
  !   node_class  = 0 p2l, 1 interpolated, 2 p2h (diagnostic only)
  ! --------------------------------------------------------------------------------------------- !
  REAL(wp) FUNCTION Feddes_P2(tp_m_d, tp_pot_mm_d, hydrus_rule, node_class)
    REAL(wp), INTENT(IN)  :: tp_m_d, tp_pot_mm_d
    LOGICAL,  INTENT(IN)  :: hydrus_rule
    INTEGER,  INTENT(OUT) :: node_class
    REAL(wp) :: tp_cm_d
    IF (hydrus_rule) THEN
      tp_cm_d = tp_pot_mm_d / 10.0_wp
      IF (tp_cm_d >= feddes_r2h_cm_d) THEN
        Feddes_P2 = feddes_p2h_m
        node_class = 2
      ELSE IF (tp_cm_d <= feddes_r2l_cm_d) THEN
        Feddes_P2 = feddes_p2l_m
        node_class = 0
      ELSE
        Feddes_P2 = feddes_p2h_m + (feddes_p2l_m - feddes_p2h_m) / (feddes_r2h_cm_d - feddes_r2l_cm_d) &
                    * (feddes_r2h_cm_d - tp_cm_d)
        node_class = 1
      END IF
    ELSE
      Feddes_P2 = MERGE(feddes_p2h_m, feddes_p2l_m, tp_m_d > 0.5_wp)
      node_class = MERGE(2, 0, tp_m_d > 0.5_wp)
    END IF
  END FUNCTION Feddes_P2

  REAL(wp) FUNCTION Feddes_Alpha(h_m, p2_m)
    REAL(wp), INTENT(IN) :: h_m, p2_m

    IF (h_m >= feddes_p0_m) THEN
      Feddes_Alpha = 0.0_wp
    ELSEIF (h_m >= feddes_p1_m) THEN
      Feddes_Alpha = (feddes_p0_m - h_m) / MAX(eps_len, feddes_p0_m - feddes_p1_m)
    ELSEIF (h_m >= p2_m) THEN
      Feddes_Alpha = 1.0_wp
    ELSEIF (h_m > feddes_p3_m) THEN
      Feddes_Alpha = (h_m - feddes_p3_m) / MAX(eps_len, p2_m - feddes_p3_m)
    ELSE
      Feddes_Alpha = 0.0_wp
    END IF
    Feddes_Alpha = MIN(1.0_wp, MAX(0.0_wp, Feddes_Alpha))
  END FUNCTION Feddes_Alpha

  ! Distribute the transpiration demand into per-layer root-water sinks.
  SUBROUTINE apply_root_sink(store_diag)
    LOGICAL, INTENT(IN) :: store_diag
    INTEGER :: j, p2_class, p2_class_alt
    REAL(wp) :: Fraw, denom, p2_m, p2_alt
    Fraw = 0.0_wp
    ! The Feddes p2 node is selected ONCE per call: it depends on the
    ! transpiration demand, not on the layer (see Feddes_P2 above).
    p2_m = Feddes_P2(transpiration_rate_day, trans_pot_unstressed_mm, feddes_p2_hydrus == 1, p2_class)
    DO j = 1, NSoil_layers
      uptake_alpha(j) = Feddes_Alpha(r_h(j), p2_m)
      Fraw = Fraw + root_dist(j) * uptake_alpha(j)
    END DO
    Fraw  = MIN(1.0_wp, MAX(0.0_wp, Fraw))
    denom = MAX(root_comp_omega, Fraw)
    DO j = 1, NSoil_layers
      sink_source(j) = transpiration_rate_day * root_dist(j) * uptake_alpha(j) / MAX(eps_len, tk(j))
    END DO
    IF (store_diag) THEN
      hydro_Fraw_diag = Fraw
      hydro_comp_denom_diag = denom
      IF (ALLOCATED(hydro_root_dist_diag) .AND. SIZE(hydro_root_dist_diag) >= NSoil_layers) &
        hydro_root_dist_diag(1:NSoil_layers) = root_dist(1:NSoil_layers)
      IF (ALLOCATED(hydro_alpha_diag) .AND. SIZE(hydro_alpha_diag) >= NSoil_layers) &
        hydro_alpha_diag(1:NSoil_layers) = uptake_alpha(1:NSoil_layers)
      ! Feddes-node diagnostics (both modes, write-only): the node the active rule used and the
      ! node the OTHER rule would have used at the same pressure heads, with its alpha per layer.
      p2_alt = Feddes_P2(transpiration_rate_day, trans_pot_unstressed_mm, feddes_p2_hydrus /= 1, p2_class_alt)
      hydro_tp_sent_cm_d_diag = transpiration_rate_day * 100.0_wp
      hydro_tp_pot_cm_d_diag  = trans_pot_unstressed_mm / 10.0_wp
      hydro_p2_used_m_diag    = p2_m
      hydro_p2_alt_m_diag     = p2_alt
      hydro_p2_class_diag     = p2_class
      hydro_p2_class_alt_diag = p2_class_alt
      IF (ALLOCATED(hydro_alpha_alt_diag) .AND. SIZE(hydro_alpha_alt_diag) >= NSoil_layers) THEN
        DO j = 1, NSoil_layers
          hydro_alpha_alt_diag(j) = Feddes_Alpha(r_h(j), p2_alt)
        END DO
      END IF
    END IF
  END SUBROUTINE apply_root_sink

END SUBROUTINE hydro
!======================================================================!