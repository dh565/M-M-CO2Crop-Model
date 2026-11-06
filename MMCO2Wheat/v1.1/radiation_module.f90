!===============================================================================!
! M&M-CO2Wheat v1.0  |  radiation_module.f90
!
! Purpose: Component radiation for the three-node energy balance. Splits incoming
!          shortwave between canopy and soil by Beer's law with separate canopy
!          and (moisture-dependent) soil albedo, solves the sky-canopy-soil
!          longwave exchange as a two-layer gray-body radiosity problem to
!          infinite order of reflection, and returns the net radiation of the
!          canopy and of the soil as functions of the solved temperatures, plus a
!          shortwave-driven soil heat flux. Time step: hourly, once per
!          energy-balance iterate. Spatial unit: canopy layer and soil surface.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision wp); SHARED (beam extinction Kdr_s, longwave
!       extinction k_lw, canopy albedo alpha_canopy_wheat, wet/dry soil albedo,
!       canopy and soil emissivities, Stefan-Boltzmann constant, top-layer
!       theta_r/theta_0, solar elevation SolarElev_s).
! Output: Rnc and Rns, the net radiation of the canopy and of the soil (W m-2),
!       the soil heat flux G (W m-2), and the SHARED rad_* ledger (shortwave
!       absorbed and reflected by canopy and soil, longwave down, emitted,
!       exchanged and net per node, outgoing longwave at the top, soil and mixed
!       albedo, shortwave and longwave transmissivity). Rnc/Rns/G are consumed by
!       energy_balance_module.f90 (the three-node temperature solution and its
!       flux repartition); the rad_* ledger is read by CANOPY_MODEL.f90 and
!       energy_balance_module.f90 for the hourly energy-balance diagnostics.
! Called by: energy_balance_module.f90 (eb_residuals, the Newton solve and the
!       actual-flux repartition), every iterate of every hourly step.
!===============================================================================!
MODULE radiation_module
  USE TYPES, ONLY: wp
  USE SHARED
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: component_radiation

CONTAINS
  ! --------------------------------------------------------------------------- !
  ! Soil-surface albedo from top-layer water content: linear between the wet and
  ! the dry value with the relative saturation of layer 1, bounded to [0.05,0.5].
  !   alpha_soil = alpha_wet + (alpha_dry - alpha_wet)*(1 - wet_frac),
  !   wet_frac = clip((theta_top - theta_r)/(theta_s - theta_r), 0, 1).
  ! The 0.12 (wet) to 0.28 (dry) end-members lie inside the bare-soil range
  ! tabulated by Monteith & Unsworth (2008) (values in SHARED).
  ! --------------------------------------------------------------------------- !
  PURE FUNCTION soil_albedo(theta_top) RESULT(a)
    REAL(wp), INTENT(IN) :: theta_top
    REAL(wp) :: a, wet_frac, thr, ths
    thr = 0.0_wp; ths = 0.40_wp
    IF (ALLOCATED(theta_r)) thr = theta_r(1)
    IF (ALLOCATED(theta_0)) ths = theta_0(1)
    wet_frac = (theta_top - thr) / MAX(1.0e-6_wp, ths - thr)
    wet_frac = MAX(0.0_wp, MIN(1.0_wp, wet_frac))
    a = alpha_soil_wet + (alpha_soil_dry_sandy - alpha_soil_wet) * (1.0_wp - wet_frac)
    a = MAX(0.05_wp, MIN(0.5_wp, a))
  END FUNCTION soil_albedo

  ! --------------------------------------------------------------------------- !
  ! Component net radiation for the canopy and soil given the solved Tcan, Tsoil.
  ! Returns Rnc, Rns, G and (via SHARED rad_* diagnostics) the full SW/LW ledger.
  ! Because the outgoing longwave is evaluated at the solved temperatures, Rnc and
  ! Rns are functions of those temperatures and the sigma*T^4 cooling feedback is
  ! carried inside the energy-balance iteration.
  ! Shortwave and longwave framework: Campbell & Norman (1998), chapters on
  ! radiation; two-source canopy/soil partition after Norman et al. (1995); the
  ! shortwave extinction coefficient after Goudriaan (1977).
  ! --------------------------------------------------------------------------- !
  SUBROUTINE component_radiation(Tcan_C, Tsoil_C, SW_in, LW_down, lai, theta_top, Rnc, Rns, G)
    REAL(wp), INTENT(IN)  :: Tcan_C, Tsoil_C, SW_in, LW_down, lai, theta_top
    REAL(wp), INTENT(OUT) :: Rnc, Rns, G
    REAL(wp) :: k_sw, tau_sw, tau_lw, a_soil, a_can, fcov
    REAL(wp) :: sw_can_abs, sw_soil_abs, sw_ref_can, sw_ref_soil
    REAL(wp) :: Tcan_k, Tsoil_k, ec, es, lw_can_net, lw_soil_net
    REAL(wp) :: lw_emit_can, lw_emit_soil, lw_can_to_soil, lw_soil_to_can
    REAL(wp) :: dSolarInc
    ! ------------------------------------------------------------------------- !
    ! Extinction for total shortwave, k_sw = Kdr*sqrt(1 - sigma) with the beam
    ! extinction coefficient Kdr and the leaf scattering coefficient sigma = 0.5
    ! (Goudriaan, 1977); longwave transmissivity exp(-k_lw*LAI) with k_lw = 0.7
    ! (Campbell & Norman, 1998).
    ! ------------------------------------------------------------------------- !
    a_soil = soil_albedo(theta_top)
    a_can  = alpha_canopy_wheat
    k_sw   = Kdr_s * SQRT(0.50_wp)                 ! sqrt(1 - sigma), sigma = 0.5
    tau_sw = EXP(-k_sw * MAX(0.0_wp, lai))         ! SW transmissivity to soil
    tau_lw = EXP(-k_lw * MAX(0.0_wp, lai))         ! LW transmissivity (k_lw SHARED)
    fcov   = 1.0_wp - tau_sw                       ! SW fractional interception by canopy
    ! ------------------------------------------------------------------------- !
    ! Shortwave (Beer's law). Shortwave reflected by the soil leaves the system:
    ! it is not re-absorbed by the canopy.
    ! ------------------------------------------------------------------------- !
    sw_can_abs  = SW_in * fcov     * (1.0_wp - a_can)
    sw_ref_can  = SW_in * fcov     * a_can
    sw_soil_abs = SW_in * tau_sw   * (1.0_wp - a_soil)
    sw_ref_soil = SW_in * tau_sw   * a_soil
    ! ------------------------------------------------------------------------- !
    ! Closed-form two-layer gray-body LONGWAVE RADIOSITY. One partially
    ! transparent canopy layer (transmittance tau_lw; the intercepted fraction has
    ! absorptivity = emissivity = ec and reflectivity 1-ec) over an opaque soil
    ! (absorptivity = emissivity = es, reflectivity 1-es). The inter-layer
    ! exchange is solved EXACTLY to infinite reflection order via the radiosity
    ! pair
    !   D = tau*S + Ec + rho_c*B      (downwelling onto the soil)
    !   B = Es + rho_s*D              (upwelling from the soil)
    ! => D = (tau*S + Ec + rho_c*Es)/(1 - rho_c*rho_s),  B = Es + rho_s*D.
    ! The conservation identity S - LW_up == lw_can_net + lw_soil_net then holds
    ! ALGEBRAICALLY, and the assertion below enforces it to floating-point
    ! precision at every call. Gray-body basis: Kirchhoff's law (absorptivity =
    ! emissivity); emissivities 0.98 (canopy) and 0.95 (soil), after
    ! Campbell & Norman (1998), longwave chapters.
    ! ------------------------------------------------------------------------- !
    Tcan_k  = Tcan_C  + 273.15_wp
    Tsoil_k = Tsoil_C + 273.15_wp
    ec = emiss_canopy; es = emiss_soil
    BLOCK
      REAL(wp) :: a_c, rho_c, rho_s, Ec_face, Es_up, denom, D_dn, B_up
      REAL(wp) :: lw_up_sky, resid_lw
      a_c     = (1.0_wp - tau_lw) * ec              ! layer absorptance per pass
      rho_c   = (1.0_wp - tau_lw) * (1.0_wp - ec)   ! layer reflectance per pass
      rho_s   = 1.0_wp - es                          ! soil reflectance
      Ec_face = a_c * const_SB * Tcan_k**4           ! canopy emission per face
      Es_up   = es  * const_SB * Tsoil_k**4          ! soil emission
      denom   = 1.0_wp - rho_c * rho_s               ! infinite-reflection closure
      D_dn    = (tau_lw * LW_down + Ec_face + rho_c * Es_up) / denom
      B_up    = Es_up + rho_s * D_dn
      lw_soil_net = es * D_dn - Es_up
      lw_can_net  = a_c * (LW_down + B_up) - 2.0_wp * Ec_face
      lw_up_sky   = rho_c * LW_down + tau_lw * B_up + Ec_face
    ! ------------------------------------------------------------------------- !
      ! Machine-precision LW conservation assertion (hard stop on violation).
    ! ------------------------------------------------------------------------- !
      resid_lw = LW_down - lw_up_sky - lw_can_net - lw_soil_net
      IF (ABS(resid_lw) > 1.0e-8_wp * MAX(100.0_wp, LW_down)) THEN
        WRITE(*,*) ' FATAL A2: LW radiosity budget violated, resid W/m2 =', resid_lw
        STOP 1
      END IF
      lw_emit_can    = Ec_face                       ! UPWARD-face canopy emission
      lw_emit_soil   = Es_up
      lw_can_to_soil = Ec_face + rho_c * B_up        ! total canopy-layer downwelling add
      lw_soil_to_can = a_c * B_up                    ! canopy absorption from below (incl. reflections)
      rad_LW_up_sky  = lw_up_sky                     ! model LWout at the top (diagnostic)
    END BLOCK

    Rnc = sw_can_abs  + lw_can_net
    Rns = sw_soil_abs + lw_soil_net
    ! ------------------------------------------------------------------------- !
    ! Soil heat flux of this routine: a fixed fraction 0.35 of the shortwave
    ! absorbed by the soil, scaled by the cosine of the solar zenith angle, and
    ! non-negative. Under energy_balance_mode == 2 the energy balance discards
    ! this value and uses the conductive flux of the one-dimensional soil-heat
    ! solution instead (soil_heat_module.f90), keeping this one as a diagnostic;
    ! it is the operative G only under energy_balance_mode == 1.
    ! ------------------------------------------------------------------------- !
    dSolarInc = (Param_pi / 2.0_wp) - SolarElev_s
    G = 0.35_wp * COS(dSolarInc) * MAX(0.0_wp, sw_soil_abs)

    ! Diagnostic ledger (SHARED).
    rad_SW_down       = SW_in
    rad_SW_abs_can    = sw_can_abs
    rad_SW_abs_soil   = sw_soil_abs
    rad_SW_ref_can    = sw_ref_can
    rad_SW_ref_soil   = sw_ref_soil
    rad_LW_down       = LW_down
    rad_LW_emit_can   = lw_emit_can                ! UPWARD-face canopy emission
    rad_LW_emit_soil  = lw_emit_soil
    rad_LW_can_to_soil = lw_can_to_soil
    rad_LW_soil_to_can = lw_soil_to_can
    rad_LW_can_net    = lw_can_net
    rad_LW_soil_net   = lw_soil_net
    rad_alpha_soil    = a_soil
    rad_alpha_mix     = MERGE((sw_ref_can + sw_ref_soil) / MAX(1.0e-6_wp, SW_in), &
                              alpha_canopy_wheat, SW_in > 1.0e-6_wp)
    rad_tau_sw        = tau_sw
    rad_tau_lw        = tau_lw
  END SUBROUTINE component_radiation

END MODULE radiation_module