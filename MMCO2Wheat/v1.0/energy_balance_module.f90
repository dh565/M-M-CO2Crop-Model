!===============================================================================!
! M&M-CO2Wheat v1.0  |  energy_balance_module.f90
!
! Purpose: The three-node canopy / soil / canopy-airspace energy balance, solved
!          once per hour and once per iterate of the canopy-temperature solve.
!          Spatial units are the canopy as one big leaf, the soil surface, and
!          the canopy airspace between them; the soil heat flux is taken from
!          the 1-D soil temperature profile, layer by layer.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision); SHARED (hourly weather at the reference
!       height, canopy structure, the resistance and stability parameters and the
!       configuration switches); radiation_module (component_radiation, the
!       canopy/soil shortwave and longwave split); soil_heat_module
!       (soil_heat_surface_flux, soil_heat_surface_coeffs, advance_soil_heat,
!       update_soil_thermal_properties).
! Output: the resistances raa, rca, rsa and rss; the solved soil-surface
!       temperature, canopy-airspace temperature Tac and vapour pressure e_ac;
!       the committed component fluxes Hc, Hs, LEc, LEs and G and the actual
!       daily repartition - consumed by CANOPY_MODEL (canopy temperature and the
!       leaf-to-air vapour pressure deficit that drive assim), by hydro through
!       the soil evaporation and transpiration demand, and by main_MMCO2Crop for
!       the hourly energy-balance and radiation diagnostic files.
! Called by: CANOPY_MODEL.f90 - solve_biophysics and compute_daily_canopy_fluxes;
!       main_MMCO2Crop.f90 - the daily repartition and the hourly writers.
!
! The canopy-airspace temperature Tac is a SOLVED node, not a post-hoc
! diagnostic, which keeps the two-source resistance-network concept without
! reconstructing a leaf temperature in a sparse canopy.
!
! Series-resistance sensible heat network (Norman et al. 1995 Agric For Meteorol
! 77:263-293; Kustas & Norman 1999 Agric For Meteorol 94:13-29):
!   Hc = rho*Cp*(Tcan -Tac )/rca
!   Hs = rho*Cp*(Tsoil-Tac )/rsa
!   Hc+Hs = rho*Cp*(Tac-Tair_ref)/raa
! Component energy balances:
!   Rnc = Hc + LEc
!   Rns = Hs + LEs + G
!
! Latent heat uses an explicit canopy-airspace vapour node e_ac:
!   LEc = (rho*Cp/gamma)*(esat(Tcan ) - e_ac)/(rcs+rca)
!   LEs = (rho*Cp/gamma)*(esat(Tsoil)*h_soil - e_ac)/(rss+rsa)
!   airspace latent closure: LEc + LEs = (rho*Cp/gamma)*(e_ac - e_ref)/raa
!
! Given Tcan and rcs, this module solves the three unknowns (Tsoil, Tac, e_ac) by
! a damped Newton iteration on a numerical 3x3 Jacobian. The outer Tcan root-find
! and the assimilation -> gsw -> rcs fixed point live in
! CANOPY_MODEL::solve_biophysics, which calls solve_eb_3node_inner once per
! fixed-point iteration so that e_ac, the canopy vapour pressure deficit and rcs
! are mutually consistent before the canopy residual is returned.

!===============================================================================!
MODULE energy_balance_module
  USE TYPES, ONLY: wp
  USE SHARED
  USE radiation_module, ONLY: component_radiation
  USE soil_heat_module, ONLY: soil_heat_surface_flux, soil_heat_surface_coeffs, &
                              advance_soil_heat, update_soil_thermal_properties
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: compute_eb_resistances, solve_eb_3node_inner, esat_pa, h_soil_surface
  PUBLIC :: raa_heat_stability, repartition_eb_actual, repartition_day_actual, write_hourly_eb_actual
  PUBLIC :: write_hourly_radiation_committed
  PUBLIC :: write_hourly_tleaf_classes

  ! Largest absolute heat-storage closure residual (J m-2) over all soil-heat
  ! advances of the run, printed once at the end by main_MMCO2Crop.
  REAL(wp), PUBLIC :: soil_heat_resid_max = 0.0_wp

CONTAINS

  ! --------------------------------------------------------------------------- !
  ! Paulson (1970) integrated stability functions for the UNSTABLE surface layer
  ! (Businger-Dyer phi forms), and the standard stable-side linear form.
  !   zeta = (z-d)/L. Unstable: zeta<0 ; stable: zeta>0.
  ! Source: Paulson, C.A. (1970), J. Appl. Meteorol. 9:857-861.
  ! --------------------------------------------------------------------------- !
  PURE SUBROUTINE mo_psi(zeta, psim, psih)
    REAL(wp), INTENT(IN)  :: zeta
    REAL(wp), INTENT(OUT) :: psim, psih
    REAL(wp) :: x
    IF (zeta < 0.0_wp) THEN
      x = (1.0_wp - 16.0_wp * zeta) ** 0.25_wp
      psim = 2.0_wp*LOG(0.5_wp*(1.0_wp + x)) + LOG(0.5_wp*(1.0_wp + x*x)) &
             - 2.0_wp*ATAN(x) + 0.5_wp*Param_pi
      psih = 2.0_wp*LOG(0.5_wp*(1.0_wp + x*x))
    ELSE
      ! Stable: Businger-Dyer linear (capped via zeta clamp in caller).
      psim = -5.0_wp * zeta
      psih = -5.0_wp * zeta
    END IF
  END SUBROUTINE mo_psi

  ! --------------------------------------------------------------------------- !
  ! Stability-corrected above-canopy aerodynamic HEAT resistance raa, controlling
  !   H = rho*Cp*(Tac - Tair_ref)/raa.
  ! The neutral SW two-layer raa (above-canopy log term + within-canopy eddy term,
  ! Shuttleworth & Wallace 1985) is corrected for atmospheric stability via
  ! Monin-Obukhov similarity: friction velocity and the above-canopy log term carry
  ! Paulson psi_m/psi_h; the within-canopy term is rescaled by the corrected u*.
  ! Iterated to consistency between u*, the Obukhov length L (from the supplied H),
  ! and raa.
  ! Optional Beljaars (1995) free-convection fallback (eb_freeconv_fallback=1) for
  ! the near-calm / non-convergent limit, using EXPLICIT eb_zi_m, eb_gust_beta.
  !
  ! Sources: Paulson 1970 (J Appl Meteorol 9:857-861); Businger et al. 1971;
  !          Beljaars 1995 (QJRMS 121:255-270) [fallback only].
  ! --------------------------------------------------------------------------- !
  FUNCTION raa_heat_stability(u_in, H, Tair_K, raa_neutral, L, zeta, psih, psim, ustar, stab_regime) RESULT(raa_h)
    REAL(wp), INTENT(IN)  :: u_in, H, Tair_K
    REAL(wp), INTENT(OUT) :: raa_neutral, L, zeta, psih, psim, ustar
    INTEGER,  INTENT(OUT) :: stab_regime   ! 1 = unstable (MO-corrected); 0 = neutral/stable (raa = neutral)
    REAL(wp) :: raa_h
    REAL(wp) :: u, hgt_eff, d, z0, n_eddy_eff, lnzd, lnzh, term2_fac
    REAL(wp) :: ust, raa_prev, pm, ph
    REAL(wp), PARAMETER :: grav = 9.81_wp, zeta_min = -5.0_wp, zeta_max = 1.0_wp
    REAL(wp), PARAMETER :: u_floor = 0.10_wp
    INTEGER :: it
    REAL(wp) :: rho_cp

    rho_cp = rho * Cp
    u = MAX(u_floor, u_in)
    hgt_eff = MAX(0.05_wp, gParam_hgt_s)
    d  = disp_height_s            ! from compute_eb_resistances (Raupach 1994)
    z0 = rough_len_s
    n_eddy_eff = n_eddy
    IF (use_dynamic_n_eddy == 1) n_eddy_eff = MIN(4.0_wp, MAX(2.0_wp, gParam_lai))

    lnzd = LOG(MAX(1.0001_wp, (gParam_refhgt - d) / z0))         ! ln((z-d)/z0)
    lnzh = LOG(MAX(1.0001_wp, (gParam_refhgt - d) / (hgt_eff - d))) ! above-canopy span
    ! --------------------------------------------------------------------- !
    ! Neutral reference (matches compute_eb_resistances raa exactly).
    ! --------------------------------------------------------------------- !
    ust = k_cons * u / lnzd
    term2_fac = (1.0_wp / (n_eddy_eff * k_cons)) * &
                (EXP(n_eddy_eff * (1.0_wp - (z0 + d)/hgt_eff)) - 1.0_wp)  ! * (1/ust)
    raa_neutral = lnzh / (k_cons * ust) + term2_fac / ust
    psim = 0.0_wp; psih = 0.0_wp; raa_h = raa_neutral; ustar = ust; stab_regime = 0
    ! --------------------------------------------------------------------- !
    ! Stability parameter for DIAGNOSTICS in both regimes (sign(H) kept).
    ! --------------------------------------------------------------------- !
    L    = -(ust**3) * rho_cp * Tair_K / (k_cons * grav * SIGN(MAX(ABS(H), 0.5_wp), H))
    zeta = MAX(zeta_min, MIN(zeta_max, (gParam_refhgt - d) / L))
    ! --------------------------------------------------------------------- !
    ! UNSTABLE-ONLY correction (H > 0, surface heating, z/L < 0): the MO correction
    ! may only REDUCE raa (enhanced daytime convective mixing). For stable/neutral
    ! (H <= 0, nights) we keep neutral raa: applying the stable-side psi here drives
    ! a nocturnal decoupling runaway (more stable -> larger raa -> colder -> more
    ! stable), which the canopy heat storage and dew/frost physics needed to bound
    ! it are not present to prevent.
    ! --------------------------------------------------------------------- !
    IF (H > 0.0_wp) THEN
      stab_regime = 1
      raa_prev = raa_neutral
      DO it = 1, 20
        ! Tiny H floor: prevents the 1/H singularity while keeping raa continuous at
        ! the H=0 gate (as H->0+, L->-inf, psi->0, raa->neutral = the stable-side
        ! value). A larger floor would make raa jump at H=0 and trap the Tcan root-find.
        L = -(ust**3) * rho_cp * Tair_K / (k_cons * grav * MAX(H, 0.01_wp))  ! L < 0
        zeta = MAX(zeta_min, MIN(0.0_wp, (gParam_refhgt - d) / L))
        CALL mo_psi(zeta, pm, ph)
        psim = pm; psih = ph
        ust = k_cons * u / MAX(0.1_wp, (lnzd - psim))
        ! Above-canopy heat term carries psi_h; within-canopy term rescaled by u*.
        raa_h = (lnzh - psih) / (k_cons * ust) + term2_fac / ust
        raa_h = MAX(5.0_wp, MIN(raa_h, raa_neutral))   ! never amplify beyond neutral
        IF (ABS(raa_h - raa_prev) < 0.5_wp) EXIT
        raa_prev = raa_h
      END DO
      ustar = ust
    END IF

  END FUNCTION raa_heat_stability

  ! --------------------------------------------------------------------------- !
  ! Saturation vapour pressure (Tetens 1930), Pa.
  ! --------------------------------------------------------------------------- !
  PURE FUNCTION esat_pa(T_C) RESULT(es)
    REAL(wp), INTENT(IN) :: T_C
    REAL(wp) :: es
    es = 610.78_wp * EXP(17.269_wp * T_C / (T_C + 237.3_wp))
  END FUNCTION esat_pa

  ! --------------------------------------------------------------------------- !
  ! Soil-surface relative humidity (0-1) from top-layer matric head via Kelvin.
  !   h_soil = exp( psi * M_w * g / (R * Tsoil_K) ), psi in meteres of water head.
  ! Falls back to a bounded saturation-fraction form if the head is unavailable.
  ! The top-layer water content is supplied by the hydrology of hydro.f90.
  ! --------------------------------------------------------------------------- !
  FUNCTION h_soil_surface(Tsoil_C) RESULT(hr)
    REAL(wp), INTENT(IN) :: Tsoil_C
    REAL(wp) :: hr
    REAL(wp) :: psi_m, Tk, arg
    REAL(wp), PARAMETER :: Mw = 0.018015_wp, grav = 9.81_wp, Rgas = 8.314_wp
    LOGICAL :: use_kelvin

    Tk = Tsoil_C + 273.15_wp
    use_kelvin = .FALSE.
    psi_m = 0.0_wp
    ! --------------------------------------------------------------------- !
    ! The head read here is the ACCEPTED (committed) Richards surface head
    ! h1_committed_m, captured after the accepted hydrology pass, and not the
    ! solver working array r_h(1), whose value is whatever the last Richards
    ! invocation left behind (a rain-substep re-initialization or a rejected
    ! trial iteration).  Any negative head is physical here, because the
    ! exponential underflow guard below bounds the arithmetic; the |psi| < 1e4 m
    ! sanity gate applies only to the r_h fallback path.
    ! --------------------------------------------------------------------- !
    IF (h1_committed_valid) THEN
      psi_m = h1_committed_m
      IF (psi_m <= 0.0_wp) use_kelvin = .TRUE.
    ELSE IF (ALLOCATED(r_h)) THEN
      IF (SIZE(r_h) >= 1) THEN
        psi_m = r_h(1)
        ! Only trust physically sane suctions (negative head, not extreme).
        IF (psi_m <= 0.0_wp .AND. psi_m > -1.0e4_wp) use_kelvin = .TRUE.
      END IF
    END IF

    IF (use_kelvin) THEN
      ! Kelvin vapour equilibrium: Philip (1957) J Meteorol 14:354-366, Eq. (1)
      ! p.355, h = exp[g Psi / R (T + 273)] with R the gas constant for water
      ! vapour; written here in the molar form exp(psi Mw g /(R_univ T)).
      ! doi:10.1175/1520-0469(1957)014<0354:EAMAHF>2.0.CO;2
      arg = psi_m * Mw * grav / (Rgas * Tk)
      arg = MAX(-50.0_wp, MIN(0.0_wp, arg))   ! guard exp underflow / wet>1
      hr = EXP(arg)
    ELSE
      ! Bounded fallback from top-layer saturation fraction (documented transitional).
      IF (ALLOCATED(theta) .AND. ALLOCATED(theta_r) .AND. ALLOCATED(theta_0)) THEN
        hr = (theta(1) - theta_r(1)) / MAX(1.0e-6_wp, theta_0(1) - theta_r(1))
      ELSE
        hr = 1.0_wp
      END IF
    END IF
    hr = MAX(0.05_wp, MIN(1.0_wp, hr))
  END FUNCTION h_soil_surface

  ! --------------------------------------------------------------------------- !
  ! Aerodynamic + surface resistances (raa, rsa, rca, rss).
  ! These depend on wind, canopy structure (Raupach 1994) and soil state only -
  ! NOT on the solved temperatures - so they are computed once per substep.
  ! The canopy stomatal resistance rcs is computed
  ! by the caller from gsw and passed to the inner solver, because it changes at
  ! every fixed-point iteration.
  ! --------------------------------------------------------------------------- !
  SUBROUTINE compute_eb_resistances(raa, rsa, rca, rss)
    REAL(wp), INTENT(OUT) :: raa, rsa, rca, rss
    REAL(wp) :: Zeroplane, RoughLen, FrictionVelocity, n_eddy_eff

    ! Bare soil / no-canopy guard.
    IF (gParam_hgt_s <= 0.0_wp .OR. gParam_lai <= 0.0_wp) THEN
      Zeroplane = 0.0_wp
      RoughLen  = z0_soil
      ! ----------------------------------------------------------------------- !
      ! w_eb_s is the wind speed the energy balance uses: it is assigned from the
      ! daily wind w_d_s at the top of every sub-step in CANOPY_MODEL and is
      ! replaced by the hourly series only when the hourly wind option is on.
      ! ----------------------------------------------------------------------- !
      FrictionVelocity = k_cons * w_eb_s / LOG(gParam_refhgt / RoughLen)
      raa = LOG(gParam_refhgt / RoughLen) / (k_cons * FrictionVelocity)
      rsa = raa
      rca = 1.0e6_wp
      rss = rss_soil_surface(theta(1))
      ! ----------------------------------------------------------------------- !
      ! Publish geometry so raa_heat_stability never sees stale/zero d,z0.
      ! ----------------------------------------------------------------------- !
      frontal_area_index_s = 0.0_wp
      disp_height_s        = Zeroplane
      rough_len_s          = RoughLen
      RETURN
    END IF
    ! ------------------------------------------------------------------------- !
    ! Raupach (1994) LAI-dependent d and z0.
    ! ------------------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: hgt_eff, lambda_f, sqrt_cdl, d_over_h, z0_over_h
      hgt_eff = MAX(0.05_wp, gParam_hgt_s)
      lambda_f = frontal_lai_factor * MAX(0.0_wp, gParam_lai)
      IF (lambda_f > 1.0e-6_wp) THEN
        sqrt_cdl = SQRT(raupach_cd1 * lambda_f)
        d_over_h = 1.0_wp - (1.0_wp - EXP(-sqrt_cdl)) / MAX(1.0e-9_wp, sqrt_cdl)
        z0_over_h = (1.0_wp - d_over_h) * &
                    EXP(-k_cons / SQRT(MAX(1.0e-9_wp, raupach_cd * lambda_f)))
        Zeroplane = d_over_h * hgt_eff
        RoughLen  = MAX(z0_soil, z0_over_h * hgt_eff)
      ELSE
        Zeroplane = 0.0_wp
        RoughLen  = z0_soil
      END IF
      frontal_area_index_s = lambda_f
      disp_height_s        = Zeroplane
      rough_len_s          = RoughLen
      FrictionVelocity = k_cons * w_eb_s / &            ! PART B: = w_d_s when rec 71 = 0
        LOG(MAX(1.0001_wp, (gParam_refhgt - Zeroplane) / RoughLen))
    END BLOCK

    n_eddy_eff = n_eddy
    IF (use_dynamic_n_eddy == 1) n_eddy_eff = MIN(4.0_wp, MAX(2.0_wp, gParam_lai))
    ! ------------------------------------------------------------------------- !
    ! rsa: soil <-> canopy air space.
    ! ------------------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: hgt_eff
      hgt_eff = MAX(0.05_wp, gParam_hgt_s)
      rsa = (hgt_eff * EXP(n_eddy_eff) / (n_eddy_eff * k_cons * FrictionVelocity * hgt_eff)) * &
            (EXP(-n_eddy_eff * z0_soil / hgt_eff) - EXP(-n_eddy_eff * (RoughLen + Zeroplane) / hgt_eff))
    END BLOCK

    ! raa: canopy air space <-> reference height.
    BLOCK
      REAL(wp) :: hgt_eff
      hgt_eff = MAX(0.05_wp, gParam_hgt_s)
      raa = (LOG((gParam_refhgt - Zeroplane)/(hgt_eff - Zeroplane)) / (k_cons * FrictionVelocity)) + &
            (hgt_eff / (n_eddy_eff * k_cons * FrictionVelocity * hgt_eff)) * &
            (EXP(n_eddy_eff * (1.0_wp - (RoughLen + Zeroplane)/hgt_eff)) - 1.0_wp)
    END BLOCK
    ! ------------------------------------------------------------------------- !
    ! rca: leaf boundary-layer resistance (Monteith 1973).
    ! The prefactor 0.012 is 2 x 0.006 (Teh 2006 Eq. 4.72 p.95 = 1.65 x 0.003459,
    ! after Mitchell 1976 Biophys J 16:561-569 and Campbell & Norman 1998): it is
    ! the 2/alpha integration factor, not a two-sided-leaf factor.  Reference
    ! implementation: Teh (2006) Code 4.3 BoundLayerResistance(), p.186.
    !
    ! Driver record 70 selects the wind-attenuation form.  Off (the shipped
    ! value), the expression above stands.  On, it becomes Teh (2006) Eq. 4.78
    ! p.97 in full, restoring the attenuation coefficient alpha in the numerator
    ! and [1 - exp(-alpha/2)] in the denominator and so reading the wind
    ! attenuation coefficient a_wind of driver record 12; alpha = 2.5 for wheat
    ! is Teh Table 4.2 p.93, "Wind attenuation coefficients alpha for some
    ! vegetation types (Cionco, 1972)", after Cionco (1972) Boundary-Layer
    ! Meteorol 3:255-263, doi:10.1007/BF02033923.
    ! ------------------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: u_h_local, hgt_eff, alpha_w, f_att, rca_free, rca_teh, wblend
      hgt_eff = MAX(0.05_wp, gParam_hgt_s)
      u_h_local = MAX(1.0e-4_wp, (FrictionVelocity/k_cons) * &
                  LOG(MAX(1.0e-9_wp, (hgt_eff - Zeroplane) / MAX(RoughLen, 1.0e-6_wp))))
      rca_free = 1.0_wp / (0.012_wp * MAX(gParam_lai, 1.0e-3_wp) * &
                 SQRT(u_h_local / MAX(gParam_leafwidth, 1.0e-4_wp)))
      IF (rca_wind_attenuation == 1) THEN
        alpha_w = MAX(RWA_ALPHA_MIN, a_wind)              ! driver record 12
        f_att   = 1.0_wp - EXP(-0.5_wp * alpha_w)         ! Teh Eq. 4.78 attenuation integral
        rca_teh = alpha_w * rca_free / f_att              ! Eq. 4.78 in full
        rwa_att_factor = alpha_w / f_att                  ! = 3.504 at alpha = 2.5
        ! Sparse-canopy blend: a declared assumption with no literature value.
        ! The uniform leaf-area density of Teh Eq. 4.75 is a closed-canopy
        ! premise, so below rwa_lai_blend the free-air limit is retained instead.
        wblend  = MIN(1.0_wp, MAX(0.0_wp, gParam_lai / MAX(rwa_lai_blend, 1.0e-6_wp)))
        rca     = wblend * rca_teh + (1.0_wp - wblend) * rca_free
        rwa_rca_ratio = rca / MAX(rca_free, 1.0e-12_wp)
        CALL rwa_alpha0_limit_check(alpha_w)
      ELSE
        rca = rca_free
      END IF
    END BLOCK
    ! ------------------------------------------------------------------------- !
    ! rss: soil-surface resistance.
    ! ------------------------------------------------------------------------- !
    rss = rss_soil_surface(theta(1))
  END SUBROUTINE compute_eb_resistances

  ! --------------------------------------------------------------------------- !
  ! Soil-surface (dry-layer vapour-diffusion) resistance.  One definition, used
  ! by both the bare-soil branch and the canopy branch, so that the two cannot
  ! drift apart.
  !
  ! Driver record 75 selects which saturated water content the expression reads.
  ! Off (the shipped value) it reads the two driver constants gParam_vwcsat and
  ! gParam_porosity; on, it reads the layer-1 saturated water content that the
  ! Richards / van Genuchten-Mualem retention already uses, for both.  The
  ! substitution of the porosity prefactor is a substitution of a different
  ! quantity and is declared as such at the SHARED.f90 declarations.
  ! --------------------------------------------------------------------------- !
  FUNCTION rss_soil_surface(theta_top) RESULT(rss_out)
    REAL(wp), INTENT(IN) :: theta_top
    REAL(wp) :: rss_out, phi_use, ths_use
    IF (rss_theta_s_consistency == 1) THEN
      ths_use = MAX(1.0e-6_wp, vgm_theta_s_L1) ! Richards' OWN layer-1 theta_s, as read from SoilLayer_vgm.csv
      phi_use = ths_use                        ! DECLARED substitution, R-4
      rss_theta_s_used = ths_use
    ELSE
      ths_use = gParam_vwcsat                  ! driver record 20
      phi_use = gParam_porosity                ! driver record 18
    END IF
    rss_out = ((2.0_wp * dry_layer_thick) / (phi_use * Dmv)) * &
              EXP(-(1.0_wp / gParam_poredist) * (theta_top / ths_use))
  END FUNCTION rss_soil_surface

  SUBROUTINE rwa_alpha0_limit_check(alpha_in)
    REAL(wp), INTENT(IN) :: alpha_in
    REAL(wp) :: a, f, ratio, tol, resid
    INTEGER  :: i
    ! Sweep the limit on a fixed ladder every time it is called, so the assertion
    ! is exercised on every run that has the flag on, independent of the alpha the
    ! deck happens to carry.
    DO i = 1, 5
      a = 1.0_wp / (10.0_wp ** REAL(i, wp))          ! 1e-1, 1e-2, ... 1e-5
      f = a / (1.0_wp - EXP(-0.5_wp * a))
      ratio = 0.5_wp * f                              ! -> 1 + alpha/4 + alpha^2/48
      resid = ABS(ratio - 1.0_wp - 0.25_wp * a)       ! second-order residual
      tol   = a * a / 24.0_wp + 1.0e-12_wp            ! bound = 2 x (alpha^2/48)
      IF (resid > tol) rwa_lim_failed = .TRUE.
      rwa_lim_worst = MAX(rwa_lim_worst, resid / tol) ! worst residual as a FRACTION of its bound
    END DO
    ! And record the factor at the deck's own alpha for the report.
    IF (alpha_in > 0.0_wp) rwa_att_factor = alpha_in / (1.0_wp - EXP(-0.5_wp * alpha_in))
  END SUBROUTINE rwa_alpha0_limit_check

  ! --------------------------------------------------------------------------- !
  ! Inner solve: given Tcan (C) and rcs (s/m), solve (Tsoil, Tac, e_ac) so that
  !   F_soil = Rns - Hs - LEs - G              = 0   (soil energy balance)
  !   F_air  = (Hc + Hs) - H                   = 0   (airspace sensible closure)
  !   F_vap  = (LEc + LEs) - LE_to_ref         = 0   (airspace latent closure)
  ! Rnc, Rns and G come from component_radiation(Tcan, Tsoil, ...) and from the
  ! soil heat solver, so they DEPEND on the solved temperatures and the numerical
  ! 3x3 Jacobian picks up dRns/dTsoil.  Tcan is the outer variable, fixed here;
  ! the outer caller forms the canopy residual F_can = Rnc - Hc - LEc.
  ! --------------------------------------------------------------------------- !
  SUBROUTINE solve_eb_3node_inner(Tcan_C, rcs, raa, rsa, rca, rss, SW_in, LW_down, lai, theta_top, &
                                  Tair_C, e_ref, &
                                  Tsoil_C, Tac_C, e_ac, &
                                  LEc, LEs, Hc, Hs, G, Rnc, Rns, &
                                  Fsoil, Fair, Fvap, iters, conv)
    REAL(wp), INTENT(IN)    :: Tcan_C, rcs, raa, rsa, rca, rss, SW_in, LW_down, lai, theta_top
    REAL(wp), INTENT(IN)    :: Tair_C, e_ref
    REAL(wp), INTENT(INOUT) :: Tsoil_C, Tac_C, e_ac   ! warm-started by caller
    REAL(wp), INTENT(OUT)   :: LEc, LEs, Hc, Hs, G, Rnc, Rns
    REAL(wp), INTENT(OUT)   :: Fsoil, Fair, Fvap
    INTEGER,  INTENT(OUT)   :: iters
    LOGICAL,  INTENT(OUT)   :: conv

    REAL(wp) :: rho_cp, coef_le
    REAL(wp) :: f(3), x(3), J(3,3), dx(3)
    REAL(wp) :: fp(3), xp(3), h_k
    REAL(wp), PARAMETER :: dT = 0.05_wp, dE = 1.0_wp, tolW = 0.5_wp
    INTEGER,  PARAMETER :: itmax = 30
    INTEGER  :: it, kk
    REAL(wp) :: step_cap

    rho_cp  = rho * Cp
    coef_le = rho_cp / gamma_const

    x = [Tsoil_C, Tac_C, e_ac]
    conv = .FALSE.
    iters = 0

    DO it = 1, itmax
      iters = it
      CALL eb_residuals(x, f)
      IF (MAXVAL(ABS(f)) < tolW) THEN
        conv = .TRUE.
        EXIT
      END IF
      ! Numerical Jacobian (columns = d f / d x_j).
      DO kk = 1, 3
        xp = x
        IF (kk == 3) THEN; h_k = dE; ELSE; h_k = dT; END IF
        xp(kk) = xp(kk) + h_k
        CALL eb_residuals(xp, fp)
        J(:,kk) = (fp - f) / h_k
      END DO
      CALL solve3x3(J, -f, dx, conv)
      IF (.NOT. conv) EXIT          ! singular Jacobian; bail with last state
      conv = .FALSE.
      ! Damp the temperature steps to avoid overshoot (e_ac step left free).
      step_cap = 10.0_wp
      dx(1) = MAX(-step_cap, MIN(step_cap, dx(1)))
      dx(2) = MAX(-step_cap, MIN(step_cap, dx(2)))
      x = x + dx
    END DO

    ! Final fluxes + component radiation at the converged (or last) state.
    CALL eb_fluxes(x, LEc, LEs, Hc, Hs)
    CALL component_radiation(Tcan_C, x(1), SW_in, LW_down, lai, theta_top, Rnc, Rns, G)
    ! Mode 2: conductive 1-D soil-heat flux from the held profile (see eb_residuals).
    IF (energy_balance_mode == 2) G = soil_heat_surface_flux(x(1))
    Fsoil = Rns - Hs - LEs - G
    Fair  = (Hc + Hs) - rho_cp * (x(2) - Tair_C) / raa
    Fvap  = (LEc + LEs) - coef_le * (x(3) - e_ref) / raa
    Tsoil_C = x(1); Tac_C = x(2); e_ac = x(3)

  CONTAINS

    SUBROUTINE eb_fluxes(xx, LEc_o, LEs_o, Hc_o, Hs_o)
      REAL(wp), INTENT(IN)  :: xx(3)
      REAL(wp), INTENT(OUT) :: LEc_o, LEs_o, Hc_o, Hs_o
      REAL(wp) :: es_can, es_soil, hsoil
      es_can  = esat_pa(Tcan_C)
      es_soil = esat_pa(xx(1))
      hsoil   = h_soil_surface(xx(1))
      Hc_o  = rho_cp * (Tcan_C - xx(2)) / rca
      Hs_o  = rho_cp * (xx(1)  - xx(2)) / rsa
      LEc_o = coef_le * (es_can  - xx(3)) / (rcs + rca)
      ! ----------------------------------------------------------------------- !
      ! Supply-limited canopy latent heat in the DEMAND energy state.  The cap is
      ! the previous iteration's COMMITTED canopy latent heat for this sub-step;
      ! the uncapped network value is kept as the potential.
      ! ----------------------------------------------------------------------- !
      LEc_unc_last = LEc_o
      IF (eb_lec_cap_use) LEc_o = MIN(LEc_o, MAX(0.0_wp, eb_lec_cap_val))
      LEs_o = coef_le * (es_soil * hsoil - xx(3)) / (rss + rsa)
      ! ----------------------------------------------------------------------- !
      ! Supply-limited soil latent heat in the DEMAND energy state.  From outer
      ! iteration 2 the cap is the previous iteration's COMMITTED LEs for this
      ! sub-step, so the demand temperatures and G see the soil evaporation the
      ! committed pass can actually sustain; the uncapped network value is
      ! retained as the potential (LEs_unc_last).  At the outer fixed point the
      ! demand and the committed LEs coincide.
      ! ----------------------------------------------------------------------- !
      LEs_unc_last = LEs_o
      IF (eb_les_cap_use) LEs_o = MIN(LEs_o, MAX(0.0_wp, eb_les_cap_val))
    END SUBROUTINE eb_fluxes

    SUBROUTINE eb_residuals(xx, ff)
      REAL(wp), INTENT(IN)  :: xx(3)
      REAL(wp), INTENT(OUT) :: ff(3)
      REAL(wp) :: LEc_l, LEs_l, Hc_l, Hs_l, H_l, LEref_l, Rnc_l, Rns_l, G_l
      CALL eb_fluxes(xx, LEc_l, LEs_l, Hc_l, Hs_l)
      ! Rns(Tsoil, Tcan) and G from component radiation (Tcan fixed in the inner solve).
      CALL component_radiation(Tcan_C, xx(1), SW_in, LW_down, lai, theta_top, Rnc_l, Rns_l, G_l)
      ! The ground heat flux is the conductive 1-D soil-heat flux of the HELD
      ! profile; the profile itself is advanced only after the repartition.
      IF (energy_balance_mode == 2) G_l = soil_heat_surface_flux(xx(1))
      H_l     = rho_cp * (xx(2) - Tair_C) / raa
      LEref_l = coef_le * (xx(3) - e_ref) / raa
      ff(1) = Rns_l - Hs_l - LEs_l - G_l      ! F_soil
      ff(2) = (Hc_l + Hs_l) - H_l             ! F_air  (sensible airspace closure)
      ff(3) = (LEc_l + LEs_l) - LEref_l       ! F_vap  (latent airspace closure)
    END SUBROUTINE eb_residuals

  END SUBROUTINE solve_eb_3node_inner

  ! --------------------------------------------------------------------------- !
  ! Solve a 3x3 linear system A x = b by Gaussian elimination with partial pivot.
  ! ok = .FALSE. if the matrix is (near-)singular.
  ! --------------------------------------------------------------------------- !
  SUBROUTINE solve3x3(A_in, b_in, x, ok)
    REAL(wp), INTENT(IN)  :: A_in(3,3), b_in(3)
    REAL(wp), INTENT(OUT) :: x(3)
    LOGICAL,  INTENT(OUT) :: ok
    REAL(wp) :: A(3,3), b(3), factor, pmax, tmpr
    INTEGER  :: i, k, p

    A = A_in; b = b_in; ok = .TRUE.; x = 0.0_wp
    DO k = 1, 2
      ! partial pivot
      p = k; pmax = ABS(A(k,k))
      DO i = k+1, 3
        IF (ABS(A(i,k)) > pmax) THEN; pmax = ABS(A(i,k)); p = i; END IF
      END DO
      IF (pmax < 1.0e-12_wp) THEN; ok = .FALSE.; RETURN; END IF
      IF (p /= k) THEN
        DO i = 1, 3; tmpr = A(k,i); A(k,i) = A(p,i); A(p,i) = tmpr; END DO
        tmpr = b(k); b(k) = b(p); b(p) = tmpr
      END IF
      DO i = k+1, 3
        factor = A(i,k) / A(k,k)
        A(i,:) = A(i,:) - factor * A(k,:)
        b(i)   = b(i)   - factor * b(k)
      END DO
    END DO
    IF (ABS(A(3,3)) < 1.0e-12_wp) THEN; ok = .FALSE.; RETURN; END IF
    x(3) = b(3) / A(3,3)
    x(2) = (b(2) - A(2,3)*x(3)) / A(2,2)
    x(1) = (b(1) - A(1,2)*x(2) - A(1,3)*x(3)) / A(1,1)
  END SUBROUTINE solve3x3

  ! --------------------------------------------------------------------------- !
  ! ACTUAL-flux repartition (closed form). Given a sub-step's net radiation
  ! split, resistances and the ACTUAL (hydrology-limited) latent components, solve
  ! the three-node temperatures. With latent PRESCRIBED, the sensible heats are the
  ! energy-balance residuals (Hc may be NEGATIVE when total canopy latent exceeds
  ! Rnc - e.g. strong wet-canopy interception evaporation supplied by sensible
  ! heat), and the temperatures follow directly:
  !   Hc  = Rnc - LEc_trans_act - LEc_int_act           (signed)
  !   Hs  = Rns - LEs_act - G                           (signed)
  !   Tac = Tair + (Hc+Hs)*raa/(rho Cp)                 (raa via unstable-only MO)
  !   Tcan= Tac  + Hc*rca/(rho Cp)
  !   Tsoil=Tac  + Hs*rsa/(rho Cp)
  ! Canopy/whole-system energy closes by construction. raa uses the same MO
  ! correction as the demand solve (d/z0 are daily-constant, set by the demand pass).
  ! --------------------------------------------------------------------------- !
  SUBROUTINE repartition_eb_actual(SW_in, LW_down, lai, theta_top, &
                                   LEc_trans_act, LEc_int_act, LEs_act, &
                                   rca, rsa, raa_neutral_in, wind, Tair_C, hgt, &
                                   Tcan, Tsoil, Tac, Hc, Hs, Rnc, Rns, G, G_placeholder, &
                                   raa_used, rp_iters, rp_resid)
    REAL(wp), INTENT(IN)    :: SW_in, LW_down, lai, theta_top
    REAL(wp), INTENT(IN)    :: LEc_trans_act, LEc_int_act, LEs_act
    REAL(wp), INTENT(IN)    :: rca, rsa, raa_neutral_in, wind, Tair_C, hgt
    REAL(wp), INTENT(INOUT) :: Tcan, Tsoil           ! warm-started by caller
    REAL(wp), INTENT(OUT)   :: Tac, Hc, Hs, Rnc, Rns, G, raa_used
    REAL(wp), INTENT(OUT)   :: G_placeholder         ! Stage-4 placeholder G (DIAGNOSTIC only)
    INTEGER,  INTENT(OUT)   :: rp_iters
    REAL(wp), INTENT(OUT)   :: rp_resid
    REAL(wp) :: rho_cp, Htot, rn, L, zeta, ph, pm, ust, Tcan_new, Tsoil_new
    REAL(wp) :: k_face, T1
    REAL(wp), PARAMETER :: rp_tol = 0.02_wp, omega = 0.4_wp   ! under-relaxation
    INTEGER  :: sreg, it
    INTEGER,  PARAMETER :: itmax = 40

    rho_cp = rho * Cp
    raa_used = MAX(5.0_wp, raa_neutral_in)
    rp_resid = 0.0_wp
    G_placeholder = 0.0_wp
    ! Mode 2: held surface-conductance coeffs (constant during this substep's solve;
    ! lambda_soil/Tsoil_prof are refreshed before the loop and advanced only after).
    k_face = 0.0_wp; T1 = 0.0_wp
    IF (energy_balance_mode == 2) CALL soil_heat_surface_coeffs(k_face, T1)
    ! ------------------------------------------------------------------------- !
    ! Rnc, Rns and G depend on (Tcan, Tsoil) through the component radiation, so
    ! the actual-flux repartition is implicit and is re-closed by fixed point.
    ! The sigma T^4 longwave feedback makes plain Picard iteration
    ! non-contractive, so the temperature update is under-relaxed by omega and
    ! clamped while iterating.
    ! ------------------------------------------------------------------------- !
    DO it = 1, itmax
      rp_iters = it
      CALL component_radiation(Tcan, Tsoil, SW_in, LW_down, lai, theta_top, Rnc, Rns, G)
      G_placeholder = G                            ! the radiation-module G, kept as a diagnostic
                                                   ! only; never fed back into the solve.
      Hc   = Rnc - LEc_trans_act - LEc_int_act     ! signed residual
      ! ----------------------------------------------------------------------- !
      ! The conductive G = k_face*(Tsoil - T1) is LINEAR in the soil-surface
      ! temperature, so it is solved IMPLICITLY rather than by explicit Picard, which
      ! is unstable for a stiff k_face. Substituting G and Tac = Tair+(Hc+Hs)*raa/rcp
      ! into the soil balance Rns - LEs = rho_cp/rsa*(Tsoil-Tac) + G gives a closed
      ! form for Hs (uses the lagged raa_used; the outer iteration re-closes raa):
      ! ----------------------------------------------------------------------- !
      Hs = (Rns - LEs_act - k_face*(Tair_C - T1) - k_face*Hc*raa_used/rho_cp) &
           / (1.0_wp + k_face*(raa_used + rsa)/rho_cp)
      Htot = Hc + Hs
      IF (lai > 0.0_wp .AND. hgt > 0.0_wp) THEN
        raa_used = raa_heat_stability(wind, Htot, Tair_C + 273.15_wp, rn, L, zeta, ph, pm, ust, sreg)
      ELSE
        raa_used = MAX(5.0_wp, raa_neutral_in)     ! bare soil: neutral (no canopy MO)
      END IF
      Tac       = Tair_C + Htot * raa_used / rho_cp
      Tcan_new  = Tac + Hc * rca / rho_cp
      Tsoil_new = Tac + Hs * rsa / rho_cp
      ! Under-relax + clamp to damp the radiative (sigma*T^4) feedback.
      Tcan_new  = omega * Tcan_new  + (1.0_wp - omega) * Tcan
      Tsoil_new = omega * Tsoil_new + (1.0_wp - omega) * Tsoil
      Tcan_new  = MAX(Tair_C - 40.0_wp, MIN(Tair_C + 60.0_wp, Tcan_new))
      Tsoil_new = MAX(Tair_C - 40.0_wp, MIN(Tair_C + 70.0_wp, Tsoil_new))
      rp_resid  = MAX(ABS(Tcan_new - Tcan), ABS(Tsoil_new - Tsoil))
      Tcan = Tcan_new; Tsoil = Tsoil_new
      IF (rp_resid < rp_tol) EXIT
    END DO
    ! ------------------------------------------------------------------------- !
    ! G is returned from the CONVERGED soil energy balance, not by re-evaluating
    ! k_face*(Tsoil - T1) on the under-relaxed and clamped Tsoil.  The soil
    ! balance the implicit Hs solves is
    !     Rns - LEs = rho_cp/rsa*(Tsoil - Tac) + G   ==>   G = Rns - LEs - Hs,
    ! and the closed form for Hs already carries the series linkage
    !     dHs/dT1 = k_face / (1 + k_face*(raa+rsa)/rho_cp),
    ! so G = Rns - LEs_act - Hs has dG/dT1 = -k_face/(1 + k_face*(raa+rsa)/rho_cp),
    ! i.e. the series conductance of the surface face with the (raa+rsa) air path.
    ! Re-deriving G from the relaxed and clamped Tsoil would sever that linkage
    ! (giving the bare, stiff dG/dT1 = -k_face) and would leave G inconsistent
    ! with (Rns, LEs, Hs) whenever the relaxation or the clamp moved Tsoil.  It is
    ! the same identity the single-node soil path uses, Hs = Rns - LEs - G.
    ! ------------------------------------------------------------------------- !
    IF (energy_balance_mode == 2) G = Rns - LEs_act - Hs
  END SUBROUTINE repartition_eb_actual
  ! --------------------------------------------------------------------------- !
  ! Write the hourly energy-balance file from the COMMITTED (actual) state.
  ! Column layout: the first 43 columns keep their names, and their physical
  ! temperature and flux slots carry the ACTUAL (repartition) values, while the
  ! solver-diagnostic slots (vapour pressure deficit, e_ac, resistances,
  ! Monin-Obukhov terms, iteration counts and residuals) stay in the demand
  ! state, because they describe the solver.  Fourteen appended *_demand columns
  ! carry the demand-state physical values explicitly.
  ! --------------------------------------------------------------------------- !
  SUBROUTINE write_hourly_eb_actual(unit_eb, day_idx)
    INTEGER, INTENT(IN) :: unit_eb, day_idx
    INTEGER :: h
    REAL(wp) :: r(43)
    DO h = 1, eb_nsteps_s
      r = ebh_row(:, h)
      ! swap physical slots to ACTUAL (indices per the buffered write list)
      r(3)  = eba_Tac(h);  r(4)  = eba_Tcan(h); r(5)  = eba_Tsoil(h)
      r(25) = eba_Rnc(h);  r(26) = eba_Rns(h);  r(27) = eba_Rnc(h) + eba_Rns(h)
      r(28) = eba_LEc(h);  r(29) = eba_LEs(h)
      r(30) = eba_Hc(h);   r(31) = eba_Hs(h);   r(32) = eba_Hc(h) + eba_Hs(h)
      r(33) = eba_G(h)
      r(34) = eba_Tac(h)  - ebh_Tair(h)
      r(35) = eba_Tcan(h) - ebh_Tair(h)
      WRITE (unit_eb,'(I0,",",I0,",",I0,57(",",G14.6))') day_idx, h, energy_balance_mode, &
        r(1:43), &
        ebh_row(3,h), ebh_row(4,h), ebh_row(5,h), ebh_row(25,h), ebh_row(26,h), &
        ebh_row(27,h), ebh_row(28,h), ebh_row(29,h), ebh_row(30,h), ebh_row(31,h), &
        ebh_row(32,h), ebh_row(33,h), ebh_row(34,h), ebh_row(35,h)
    END DO
  END SUBROUTINE write_hourly_eb_actual

  ! --------------------------------------------------------------------------- !
  ! Write hourly_radiation.csv from the COMMITTED scene.  The named longwave and
  ! net-radiation columns take the accepted-repartition radiosity state (eba_rad,
  ! captured per sub-step in repartition_day_actual), and the demand-pass values
  ! (rad_row, buffered in the demand loop) are appended as *_demand.  The
  ! incoming shortwave and longwave forcing are scene-independent and come from
  ! the buffered row.
  ! --------------------------------------------------------------------------- !
  SUBROUTINE write_hourly_radiation_committed(day_idx)
    INTEGER, INTENT(IN) :: day_idx
    INTEGER :: h
    REAL(wp) :: r(22)
    IF (.NOT. rad_open) RETURN
    DO h = 1, eb_nsteps_s
      r = rad_row(:, h)
      ! ----------------------------------------------------------------------- !
      ! swap LW/Rn slots to the COMMITTED scene (indices per the buffered list)
      ! ----------------------------------------------------------------------- !
      r(8)  = eba_rad(1,h);  r(9)  = eba_rad(2,h)
      r(10) = eba_rad(3,h);  r(11) = eba_rad(4,h)
      r(12) = eba_rad(5,h);  r(13) = eba_rad(6,h)
      r(14) = eba_rad(7,h)
      r(15) = eba_rad(8,h);  r(16) = eba_rad(9,h);  r(17) = eba_rad(10,h)
      WRITE (unit_rad,'(I0,",",I0,32(",",G14.6))') day_idx, h, &
        r(1:22), &
        rad_row(8,h), rad_row(9,h), rad_row(10,h), rad_row(11,h), rad_row(12,h), &
        rad_row(13,h), rad_row(14,h), rad_row(15,h), rad_row(16,h), rad_row(17,h)
    END DO
  END SUBROUTINE write_hourly_radiation_committed

  ! --------------------------------------------------------------------------- !
  ! Write output/hourly_tleaf_classes.csv (driver record 95). The class
  ! temperatures are DEMAND-state, the only state that carries a vapour node.
  ! The committed bulk Tcan is written on the SAME row so that the reader can
  ! form either basis, and the committed-basis class columns are the declared
  ! linear transfer Tcan_committed + (T_class_demand - Tcan_demand).
  ! --------------------------------------------------------------------------- !
  SUBROUTINE write_hourly_tleaf_classes(day_idx)
    INTEGER, INTENT(IN) :: day_idx
    INTEGER  :: h
    REAL(wp) :: tsun, tshd, fsun, tcd, tca, tam, trm, dsun, dshd
    IF (.NOT. tlc_open) RETURN
    DO h = 1, eb_nsteps_s
      IF (ebh_tlc(14, h) < 0.5_wp) CYCLE            ! no canopy at this substep
      tsun = ebh_tlc(1, h); tshd = ebh_tlc(2, h); fsun = ebh_tlc(3, h)
      tcd  = ebh_tlc(12, h); tca = eba_Tcan(h)
      tam  = fsun * tsun + (1.0_wp - fsun) * tshd                     ! area mean
      trm  = (fsun * (tsun + 273.15_wp)**4 &
              + (1.0_wp - fsun) * (tshd + 273.15_wp)**4) ** 0.25_wp - 273.15_wp
      dsun = tca + (tsun - tcd)                                        ! DECLARED transfer
      dshd = tca + (tshd - tcd)
      WRITE (unit_tlc,'(I0,",",I0,26(",",G14.6))') day_idx, h, &
        ebh_hour(h), ebh_lai(h), fsun, ebh_Tair(h), ebh_tlc(10,h), ebh_tlc(11,h), &
        tcd, tca, tsun, tshd, tam, trm, dsun, dshd, &
        (tsun - tshd), (tsun - tcd), (tsun - ebh_Tair(h)), (tca - ebh_Tair(h)), &
        ebh_tlc(4,h), ebh_tlc(5,h), ebh_tlc(6,h), ebh_tlc(7,h), ebh_tlc(8,h), &
        (ebh_tlc(9,h) - tcd), ebh_tlc(15,h), ebh_tlc(16,h)
    END DO
  END SUBROUTINE write_hourly_tleaf_classes

  ! --------------------------------------------------------------------------- !
  ! Redistribute the DAILY actual fluxes over the stored sub-steps by potential
  ! latent-heat weight (uniform actual/potential scaling per component) and
  ! repartition each sub-step's temperatures with the actual latent heat.  The
  ! weighted daily means of Tcan, Tac, Tsoil, Hc and Hs and the canopy
  ! energy-closure residual are aggregated into the SHARED *_act_day
  ! diagnostics.  The f_* arguments are the actual/potential ratios (at most 1,
  ! except f_int, which is the intercepted evaporation over its demand).
  ! --------------------------------------------------------------------------- !
  SUBROUTINE repartition_day_actual(f_trans, f_int, f_soil, is_accepted)
    REAL(wp), INTENT(IN) :: f_trans, f_int, f_soil
    ! ------------------------------------------------------------------------- !
    ! is_accepted = .TRUE. only on the accepted, final outer evapotranspiration
    ! and Richards iteration.  The prognostic soil-heat profile is advanced, its
    ! properties refreshed and its diagnostics written ONCE per day, on that
    ! accepted pass only.  Trial outer iterations compute the energy balance and
    ! the repartition against the HELD profile and must not advance or alter
    ! Tsoil_prof, Cv_soil or lambda_soil.
    ! ------------------------------------------------------------------------- !
    LOGICAL,  INTENT(IN) :: is_accepted
    INTEGER  :: h, rp_it
    REAL(wp) :: w, wsum, LEc_tr, LEc_in, LEs_a, Tcan, Tac, Tsoil, Hc, Hs
    REAL(wp) :: Rnc_h, Rns_h, G_h, raa_u, resid, rp_res
    REAL(wp) :: dt_s, heat_resid_h, G_placeholder_h, storage_change_h
    REAL(wp), ALLOCATABLE :: tsp_save(:), tsp_day_save(:)
    ! Daily soil-heat aggregates (J m-2, C).
    REAL(wp) :: shd_Gint, shd_store, shd_Tsurf_sum, shd_Ttop_sum, shd_w
    REAL(wp) :: shd_Tsurf_min, shd_Tsurf_max

    Tcan_act_day = 0.0_wp; Tac_act_day = 0.0_wp; Tsoil_act_day = 0.0_wp
    Hc_act_day = 0.0_wp; Hs_act_day = 0.0_wp; eb_canopy_resid_day = 0.0_wp
    eb_repartition_iters = 0; eb_repartition_resid = 0.0_wp
    ! ------------------------------------------------------------------------- !
    ! Refresh Cv_soil and lambda_soil from this iteration's post-hydrology theta
    ! on EVERY repartition entry, so that the next demand pass conducts with the
    ! same k_face vintage as the committed trajectory it is paired with
    ! (eba_node1).  Refreshing only on the accepted pass breaks the per-sub-step
    ! G identity on wet days, where irrigation or rain steps the conductivity, by
    ! order 100 W m-2 through the conductivity vintage alone.  theta converges
    ! over the outer loop, so the persistent accepted state is unchanged.
    ! ------------------------------------------------------------------------- !
    IF (energy_balance_mode == 2) CALL update_soil_thermal_properties(NSoil_layers)
    ! ------------------------------------------------------------------------- !
    ! The committed pass conducts against the LIVE advancing profile, never
    ! against the previous-iteration override armed for the demand pass; it is
    ! disarmed here and the demand loop re-arms it per sub-step next iteration.
    ! ------------------------------------------------------------------------- !
    shf_node1_use = .FALSE.
    eb_les_cap_use = .FALSE.
    eb_lec_cap_use = .FALSE.
    ! ------------------------------------------------------------------------- !
    ! The profile advances per sub-step on EVERY outer iteration, so that
    ! eba_node1(:) records the committed node-1 trajectory the NEXT demand pass
    ! conducts against.  Non-accepted iterations restore the day-start profile on
    ! exit, so the prognostic state still evolves exactly once per day.
    ! ------------------------------------------------------------------------- !
    IF (energy_balance_mode == 2) tsp_day_save = Tsoil_prof
    ! Zero the daily soil-heat aggregates.
    shd_Gint = 0.0_wp; shd_store = 0.0_wp; shd_Tsurf_sum = 0.0_wp; shd_Ttop_sum = 0.0_wp
    shd_w = 0.0_wp; shd_Tsurf_min = HUGE(1.0_wp); shd_Tsurf_max = -HUGE(1.0_wp)
    wsum = 0.0_wp
    ! Zero the daily ACTUAL energy sums, accumulated below per sub-step.
    Rn_day_act = 0.0_wp; G_day_act = 0.0_wp
    LEc_day_act = 0.0_wp; LEs_day_act = 0.0_wp
    Hc_day_act = 0.0_wp; Hs_day_act = 0.0_wp
    DO h = 1, eb_nsteps_s
      w = ebh_weight(h)
      ! ----------------------------------------------------------------------- !
      ! The committed transpiration latent heat is the ROOT-SINK MASS that
      ! Richards actually delivered in this window, converted with the same
      ! L(Tair) as the soil term.  f_trans is a diagnostic here; interception
      ! keeps its own f_int scaling.
      ! ----------------------------------------------------------------------- !
      IF (esoil_bin_valid) THEN
        BLOCK
          REAL(wp) :: L_c, dtc
          L_c = (2500.8_wp - 2.36_wp*ebh_Tair(h) + 0.0016_wp*ebh_Tair(h)**2) * 1000.0_wp
          dtc = ebh_weight(h) * 86400.0_wp
          LEc_tr = etr_bin_mm(h) * L_c / dtc
        END BLOCK
      ELSE
        LEc_tr = ebh_LEc_trans_pot(h) * f_trans
      END IF
      LEc_in = ebh_LEc_int_pot(h)   * f_int
      ! ----------------------------------------------------------------------- !
      ! The committed soil latent heat is the MASS that Richards delivered in
      ! this sub-step's window (binned, conserving the daily total to round-off),
      ! converted with L(Tair_h), the same latent-heat polynomial as the demand
      ! mass conversion, evaluated at a pre-solve temperature so that the
      ! conversion is not circular.  It is a flux the soil produced, not the
      ! potential shape scaled by f_soil.  Because the demand forcing is hourly
      ! the delivered series is diurnal, so an exceedance of the sub-step
      ! potential (counted below) can only come from window misalignment.
      ! ----------------------------------------------------------------------- !
      IF (esoil_bin_valid) THEN
        BLOCK
          REAL(wp) :: L_h, dtj
          L_h = (2500.8_wp - 2.36_wp*ebh_Tair(h) + 0.0016_wp*ebh_Tair(h)**2) * 1000.0_wp
          dtj = ebh_weight(h) * 86400.0_wp
          LEs_a = esoil_bin_mm(h) * L_h / dtj
          IF (LEs_a > ebh_LEs_pot(h) + 1.0e-9_wp .AND. is_accepted) THEN
            les_bin_exceed_count = les_bin_exceed_count + 1
            IF (ebh_LEs_pot(h) > 1.0e-9_wp) &
              les_bin_exceed_max = MAX(les_bin_exceed_max, LEs_a / ebh_LEs_pot(h))
          END IF
        END BLOCK
      ELSE
        LEs_a  = ebh_LEs_pot(h)       * f_soil
      END IF
      Tcan  = ebh_Tcan_pot(h)        ! warm start from demand temperatures
      Tsoil = ebh_Tair(h)
      ! ----------------------------------------------------------------------- !
      ! Restore THIS sub-step's solar geometry before the implicit
      ! component_radiation call: otherwise the SHARED globals still hold the
      ! last demand-loop sub-step's night values, Kdr_s = 0, and the canopy
      ! absorbs zero shortwave in every actual-flux computation.
      ! ----------------------------------------------------------------------- !
      Kdr_s       = ebh_Kdr(h)
      SolarElev_s = ebh_SolarElev(h)
      ! ----------------------------------------------------------------------- !
      ! Rnc, Rns and G are recomputed from the solved temperatures through the
      ! component radiation, so they are implicit; repartition_eb_actual
      ! re-closes them by fixed point.
      ! ----------------------------------------------------------------------- !
      CALL repartition_eb_actual(ebh_SW_in(h), ebh_LW_down(h), ebh_lai(h), ebh_theta_top(h), &
             LEc_tr, LEc_in, LEs_a, ebh_rca(h), ebh_rsa(h), ebh_raa_neutral(h), &
             ebh_wind(h), ebh_Tair(h), ebh_hgt(h), &
             Tcan, Tsoil, Tac, Hc, Hs, Rnc_h, Rns_h, G_h, G_placeholder_h, raa_u, rp_it, rp_res)
      ! Capture the per-sub-step ACTUAL state and accumulate the daily ACTUAL
      ! energy ledger (W m-2 -> MJ m-2 d-1 through weight*86400*1e-6).
      eba_Tac(h) = Tac; eba_Tcan(h) = Tcan; eba_Tsoil(h) = Tsoil
      eba_Rnc(h) = Rnc_h; eba_Rns(h) = Rns_h
      ! Capture the COMMITTED-scene longwave ledger: the rad_* globals hold the
      ! converged repartition_eb_actual radiosity state, so hourly_radiation.csv
      ! reports outgoing longwave and net radiation on the committed scene.
      eba_rad(1,h) = rad_LW_emit_can;   eba_rad(2,h) = rad_LW_emit_soil
      eba_rad(3,h) = rad_LW_can_to_soil; eba_rad(4,h) = rad_LW_soil_to_can
      eba_rad(5,h) = rad_LW_can_net;    eba_rad(6,h) = rad_LW_soil_net
      eba_rad(7,h) = rad_LW_up_sky
      eba_rad(8,h) = Rnc_h; eba_rad(9,h) = Rns_h; eba_rad(10,h) = Rnc_h + Rns_h
      eba_LEc(h) = LEc_tr + LEc_in; eba_LEs(h) = LEs_a
      eba_Hc(h) = Hc; eba_Hs(h) = Hs; eba_G(h) = G_h
      BLOCK
        REAL(wp) :: wmj
        wmj = ebh_weight(h) * 86400.0_wp * 1.0e-6_wp
        Rn_day_act  = Rn_day_act  + (Rnc_h + Rns_h) * wmj
        G_day_act   = G_day_act   + G_h * wmj
        LEc_day_act = LEc_day_act + (LEc_tr + LEc_in) * wmj
        LEs_day_act = LEs_day_act + LEs_a * wmj
        Hc_day_act  = Hc_day_act  + Hc * wmj
        Hs_day_act  = Hs_day_act  + Hs * wmj
      END BLOCK
      ! ----------------------------------------------------------------------- !
      ! Advance the prognostic soil-temperature profile with this sub-step's
      ! COMMITTED ACTUAL surface flux G_h on EVERY outer iteration, so that
      ! eba_node1(:) records the committed node-1 trajectory the next demand pass
      ! conducts against.  dt_s is the sub-step clock duration (3600 s on the
      ! hourly integration).  Non-accepted iterations restore the day-start
      ! profile on exit below, so the prognostic state still evolves exactly once
      ! per day, and the counters and diagnostics stay once per day.
      ! ----------------------------------------------------------------------- !
      IF (energy_balance_mode == 2) THEN
        ! Capture the START-of-sub-step node value, the state THIS sub-step's
        ! implicit solve conducted against (soil_heat_surface_coeffs reads it
        ! before the advance), for the per-sub-step demand/committed G identity.
        eba_node1(h) = Tsoil_prof(1)
        dt_s = ebh_weight(h) * 86400.0_wp
        tsp_save = Tsoil_prof   ! auto-(re)allocates; guard restore point
        CALL advance_soil_heat(NSoil_layers, dt_s, G_h, heat_resid_h)
        IF (.NOT. (ABS(heat_resid_h) < 1.0e30_wp)) THEN
          ! Non-finite advance: hold the profile for this sub-step instead of
          ! aborting, and count and report the event.
          Tsoil_prof = tsp_save
          ! ------------------------------------------------------------------- !
          ! The residual is NOT zeroed.  heat_resid_h is the channel that reports
          ! the failure (soil_heat_resid_max, storage_change_h and
          ! soil_heat_hourly.csv), so overwriting it with 0.0 would make a
          ! non-finite advance indistinguishable from a perfectly closed one in
          ! every published residual.
          ! ------------------------------------------------------------------- !
          IF (is_accepted) THEN
            soil_heat_fb_events = soil_heat_fb_events + 1
            IF (iD /= soil_heat_fb_lastday) THEN
              soil_heat_fb_days = soil_heat_fb_days + 1
              soil_heat_fb_lastday = iD
              WRITE (*,'(A,I0)') ' SOILHEAT FALLBACK day=', iD
            END IF
          END IF
        END IF
        IF (is_accepted) THEN
          soil_heat_resid_max = MAX(soil_heat_resid_max, ABS(heat_resid_h))
          ! ------------------------------------------------------------------- !
          ! Soil-heat diagnostics.  The storage change follows from the advance
          ! residual: heat_resid = storage_change - dt*G_surface with G_bottom = 0,
          ! so storage_change = heat_resid + dt*G_surface.  The node temperatures
          ! written here are post-advance.
          ! ------------------------------------------------------------------- !
          storage_change_h = heat_resid_h + dt_s * G_h
          IF (soilheat_h_open) THEN
            WRITE (unit_soilheat_h,'(I0,",",G0,14(",",G14.6))') &
              iD, ebh_hour(h), G_h, 0.0_wp, storage_change_h, heat_resid_h, &
              Tsoil, Tsoil_prof(1), Tsoil_prof(2), Tsoil_prof(3), Tsoil_prof(NSoil_layers), &
              Rns_h, Hs, LEs_a, G_h, G_placeholder_h
          END IF
          shd_Gint      = shd_Gint  + G_h * dt_s
          shd_store     = shd_store + storage_change_h
          shd_Tsurf_sum = shd_Tsurf_sum + Tsoil * w
          shd_Ttop_sum  = shd_Ttop_sum  + Tsoil_prof(1) * w
          shd_Tsurf_min = MIN(shd_Tsurf_min, Tsoil)
          shd_Tsurf_max = MAX(shd_Tsurf_max, Tsoil)
          shd_w         = shd_w + w
        END IF
      END IF
      resid = ABS(Rnc_h - Hc - LEc_tr - LEc_in)        ! ~0 by construction
      Tcan_act_day        = Tcan_act_day        + Tcan  * w
      Tac_act_day         = Tac_act_day         + Tac   * w
      Tsoil_act_day       = Tsoil_act_day       + Tsoil * w
      Hc_act_day          = Hc_act_day          + Hc    * w
      Hs_act_day          = Hs_act_day          + Hs    * w
      eb_canopy_resid_day = eb_canopy_resid_day + resid * w
      eb_repartition_iters = MAX(eb_repartition_iters, rp_it)
      eb_repartition_resid = MAX(eb_repartition_resid, rp_res)
      wsum = wsum + w
    END DO
    IF (wsum > 0.0_wp) THEN
      Tcan_act_day = Tcan_act_day / wsum; Tac_act_day = Tac_act_day / wsum
      Tsoil_act_day = Tsoil_act_day / wsum; Hc_act_day = Hc_act_day / wsum
      Hs_act_day = Hs_act_day / wsum; eb_canopy_resid_day = eb_canopy_resid_day / wsum
    END IF
    ! ------------------------------------------------------------------------ !
    ! Non-accepted iterations restore the day-start profile, so the prognostic
    ! state advances exactly once per day, on the accepted pass.  eba_node1(:)
    ! keeps this iteration's committed trajectory for the next demand pass
    ! whether or not the profile is restored.
    ! ------------------------------------------------------------------------ !
    IF (energy_balance_mode == 2) THEN
      IF (.NOT. is_accepted) Tsoil_prof = tsp_day_save
      eba_node1_valid = .TRUE.
    END IF
    ! Daily soil-heat roll-up.  G_bottom = 0, so the daily heat residual is the
    ! storage change minus the surface-flux integral.
    IF (energy_balance_mode == 2 .AND. is_accepted .AND. soilheat_d_open .AND. shd_w > 0.0_wp) THEN
      WRITE (unit_soilheat_d,'(I0,7(",",G14.6))') &
        iD, shd_Gint, shd_store, (shd_store - shd_Gint), &
        shd_Tsurf_sum / shd_w, shd_Tsurf_min, shd_Tsurf_max, shd_Ttop_sum / shd_w
    END IF
  END SUBROUTINE repartition_day_actual

END MODULE energy_balance_module