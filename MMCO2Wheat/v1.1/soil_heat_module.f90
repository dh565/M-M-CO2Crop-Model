!===============================================================================!
! M&M-CO2Wheat v1.0  |  soil_heat_module.f90
!
! Purpose: One-dimensional soil heat conduction. Builds the heat grid (the
!          hydrology layers extended geometrically to about 5 m), initializes the
!          prognostic temperature profile, updates the de Vries volumetric heat
!          capacity and the Johansen/Cote-Konrad thermal conductivity from the
!          current water content, supplies the linear surface-conductance
!          coefficients for the energy balance, and advances the profile with a
!          backward-Euler finite-volume step solved by the Thomas algorithm.
!          Time step: one energy-balance sub-step (hourly). Spatial unit: soil
!          layer on the heat grid.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision wp); SHARED (water content theta, porosity
!       theta_0, residual theta_r, layer thicknesses tk/dz, the heat-grid state
!       Tsoil_prof/Cv_soil/lambda_soil/tk_heat/n_heat_tot, the soil thermal
!       constants, the initialization controls and the surface-node override).
! Output: Tsoil_prof(:) soil temperature profile (C), Cv_soil(:) volumetric heat
!       capacity (J m-3 K-1), lambda_soil(:) thermal conductivity (W m-1 K-1),
!       tk_heat(:)/n_heat_tot heat-grid geometry, the surface conductance k_face
!       and node-1 temperature, the conductive surface flux G (W m-2) and the
!       storage residual heat_resid (J m-2). The fluxes and coefficients are
!       consumed by energy_balance_module.f90; the profile and the residual are
!       read by main_MMCO2Crop.f90 for the soil-temperature output and the heat
!       closure diagnostics.
! Called by: main_MMCO2Crop.f90 (initialize_soil_heat, once before the day loop);
!       energy_balance_module.f90 (update_soil_thermal_properties,
!       soil_heat_surface_coeffs, soil_heat_surface_flux, advance_soil_heat,
!       every energy-balance sub-step).

!===============================================================================!
MODULE soil_heat_module
  USE TYPES,  ONLY: wp
  USE SHARED
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_soil_heat, update_soil_thermal_properties, advance_soil_heat
  PUBLIC :: soil_heat_surface_flux, soil_heat_surface_coeffs

  ! Cote-Konrad (2005) empirical constants for natural mineral soils. Module
  ! constants, not driver inputs; ck_kappa is the unfrozen coarse-soil value of
  ! their Table 4.
  REAL(wp), PARAMETER :: ck_chi   = 0.75_wp ! lambda_dry pre-factor          (Cote-Konrad 2005)
  REAL(wp), PARAMETER :: ck_eta   = 1.2_wp  ! lambda_dry porosity exponent   (Cote-Konrad 2005)
  REAL(wp), PARAMETER :: ck_kappa = 4.60_wp ! Kersten kappa, unfrozen coarse (Cote-Konrad 2005, Table 4)
  ! Damping depth used only by the synthetic initial profile; the prognostic
  ! solver then carries the profile in time, so it is a starting guess.
  REAL(wp), PARAMETER :: z_damp_init = 1.0_wp ! m
  ! No bulk organic volume fraction is supplied, so the de Vries organic term is zero.
  REAL(wp), PARAMETER :: f_org = 0.0_wp        ! bulk organic volume fraction

CONTAINS

  ! --------------------------------------------------------------------------- !
  ! Bounded clamp helper: y = min(max(x,lo),hi).
  PURE FUNCTION clamp(x, lo, hi) RESULT(y)
    REAL(wp), INTENT(IN) :: x, lo, hi
    REAL(wp) :: y
    y = MIN(MAX(x, lo), hi)
  END FUNCTION clamp

  ! --------------------------------------------------------------------------- !
  ! Volumetric heat capacity of one layer (J m-3 K-1) by volumetric mixing of the
  ! constituents, Cv = sum_j x_j (rho c)_j:
  !   solids = 1-phi;  x_org = f_org (=0 here);  x_mineral = max(0,(1-phi)-f_org);
  !   theta_air = max(0, phi - theta).
  ! de Vries, D.A. (1963), "Thermal properties of soils", in van Wijk (ed.),
  ! Physics of Plant Environment, North-Holland, pp.210-235, Table 7.1 (the
  ! constituent capacities Cv_mineral/Cv_org/Cv_water/Cv_air live in SHARED).
  PURE FUNCTION devries_cv(theta_l, phi) RESULT(cv)
    REAL(wp), INTENT(IN) :: theta_l, phi
    REAL(wp) :: cv, th, x_air, x_min
    th    = clamp(theta_l, 0.0_wp, phi)
    x_air = MAX(0.0_wp, phi - th)
    x_min = MAX(0.0_wp, (1.0_wp - phi) - f_org)
    cv = x_min*Cv_mineral + f_org*Cv_org + th*Cv_water + x_air*Cv_air
  END FUNCTION devries_cv

  ! --------------------------------------------------------------------------- !
  ! Bounded thermal conductivity of one layer (W m-1 K-1), the normalized
  ! (Kersten-number) form of Johansen (1975) with the generalized parameters of
  ! Cote, J. & Konrad, J.-M. (2005), Can. Geotech. J. 42(2):443-458,
  ! doi:10.1139/t04-106:
  !   lambda       = lambda_dry + (lambda_sat - lambda_dry)*Ke
  !   lambda_sat   = lambda_solid^(1-phi) * lambda_water^phi   (geometric mean)
  !   lambda_solid = lambda_quartz^q * lambda_other^(1-q)      (Johansen 1975)
  !   lambda_dry   = ck_chi * 10^(-ck_eta*phi)   (natural mineral soils)
  !   Ke           = ck_kappa*Sr / (1 + (ck_kappa-1)*Sr)   (Sr = eff. saturation)
  ! The end-member conductivities lambda_quartz/lambda_other/lambda_water are
  ! Johansen, O. (1975), "Thermal conductivity of soils", PhD thesis, Trondheim
  ! (values in SHARED).
  PURE FUNCTION cote_konrad_lambda(theta_l, phi, thr) RESULT(lam)
    REAL(wp), INTENT(IN) :: theta_l, phi, thr
    REAL(wp) :: lam, sr, ke, lam_solid, lam_sat, lam_dry, q, ph
    ph        = clamp(phi, 0.05_wp, 0.95_wp)                 ! keep porosity physical
    q         = clamp(quartz_fraction, 0.0_wp, 1.0_wp)       ! quartz fraction in [0,1]
    sr        = clamp((theta_l - thr) / MAX(1.0e-6_wp, phi - thr), 0.0_wp, 1.0_wp)
    lam_solid = lambda_quartz ** q * lambda_other ** (1.0_wp - q)
    lam_sat   = lam_solid ** (1.0_wp - ph) * lambda_water ** ph
    lam_dry   = ck_chi * 10.0_wp ** (-ck_eta * ph)
    ke        = clamp(ck_kappa * sr / (1.0_wp + (ck_kappa - 1.0_wp) * sr), 0.0_wp, 1.0_wp)
    lam       = lam_dry + (lam_sat - lam_dry) * ke
    lam       = clamp(lam, lam_dry, lam_sat)                 ! bounded by construction
  END FUNCTION cote_konrad_lambda

  ! --------------------------------------------------------------------------- !
  ! Recompute Cv_soil(:) and lambda_soil(:) from the current SHARED theta(:),
  ! porosity theta_0(:), and residual theta_r(:). Allocation/size guarded so the
  ! routine never writes out of bounds.
  SUBROUTINE update_soil_thermal_properties(nlayer)
    INTEGER, INTENT(IN) :: nlayer
    INTEGER :: i, n
    IF (nlayer < 1) RETURN
    IF (.NOT. (ALLOCATED(Cv_soil) .AND. ALLOCATED(lambda_soil))) RETURN
    IF (.NOT. (ALLOCATED(theta) .AND. ALLOCATED(theta_0) .AND. ALLOCATED(theta_r))) RETURN
    n = MIN(nlayer, SIZE(theta), SIZE(theta_0), SIZE(theta_r))
    DO i = 1, MIN(SIZE(Cv_soil), SIZE(lambda_soil), MAX(n, n_heat_tot))
      BLOCK
        INTEGER :: ih
        ih = MIN(i, n)   ! below the hydrology column the heat grid reuses the
                         ! bottom hydrology layer's water content and porosity
        Cv_soil(i)     = devries_cv(theta(ih), theta_0(ih))
        lambda_soil(i) = cote_konrad_lambda(theta(ih), theta_0(ih), theta_r(ih))
      END BLOCK
    END DO
  END SUBROUTINE update_soil_thermal_properties

  ! --------------------------------------------------------------------------- !
  ! Allocate (idempotent) and initialize the prognostic profile Tsoil_prof(:),
  ! then populate Cv_soil/lambda_soil. Reads soil_T_init_method, Tdeep_C (SHARED)
  ! and the node-center depths ABS(dz(:)). Tsurf_C is the surface temperature
  ! supplied by the caller.
  SUBROUTINE initialize_soil_heat(nlayer, Tsurf_C)
    INTEGER,  INTENT(IN) :: nlayer
    REAL(wp), INTENT(IN) :: Tsurf_C
    INTEGER  :: i
    REAL(wp) :: z_i
    IF (nlayer < 1) RETURN
    ! Heat domain: the heat grid is the hydrology layers plus a geometric
    ! extension (first extra layer 0.25 m, ratio 1.4) down to
    ! heat_domain_target_m, so that the zero-flux bottom sits well below the
    ! seasonal wave (about 5 m with the default target).
    BLOCK
      REAL(wp) :: zbot, thick_next
      INTEGER  :: nx
      zbot = 0.0_wp
      DO i = 1, nlayer
        IF (ALLOCATED(tk) .AND. i <= SIZE(tk)) zbot = zbot + ABS(tk(i))
      END DO
      nx = 0; thick_next = 0.25_wp
      IF (heat_domain_target_m > zbot) THEN
        BLOCK
          REAL(wp) :: z_acc
          z_acc = zbot
          DO WHILE (z_acc < heat_domain_target_m .AND. nx < 30)
            nx = nx + 1; z_acc = z_acc + thick_next; thick_next = thick_next * 1.4_wp
          END DO
        END BLOCK
      END IF
      n_heat_tot = nlayer + nx
      IF (ALLOCATED(tk_heat)) DEALLOCATE(tk_heat)
      ALLOCATE (tk_heat(n_heat_tot))
      thick_next = 0.25_wp
      DO i = 1, n_heat_tot
        IF (i <= nlayer) THEN
          tk_heat(i) = MAX(1.0e-9_wp, ABS(tk(i)))
        ELSE
          tk_heat(i) = thick_next; thick_next = thick_next * 1.4_wp
        END IF
      END DO
      IF (nx > 0) WRITE (*,'(A,I0,A,F6.2,A)') ' DEEP-HEAT: ', nx, ' extra layers to ', &
        SUM(tk_heat), ' m (zero-flux bottom at the extended depth)'
    END BLOCK
    IF (ALLOCATED(Tsoil_prof)) DEALLOCATE(Tsoil_prof)
    IF (ALLOCATED(Cv_soil)) DEALLOCATE(Cv_soil)
    IF (ALLOCATED(lambda_soil)) DEALLOCATE(lambda_soil)
    ALLOCATE (Tsoil_prof(n_heat_tot), Cv_soil(n_heat_tot), lambda_soil(n_heat_tot))
    ! Measured initial profile: where input/tsoil_init.csv supplies at least two
    ! depth-temperature anchors, the profile is interpolated linearly in depth
    ! between them instead of taking the synthetic form below.
    IF (tsoil_init_from_file .AND. tsoil_init_n >= 2) THEN
      DO i = 1, n_heat_tot
        BLOCK
          INTEGER :: kk
          z_i = 0.0_wp
          DO kk = 1, i
            z_i = z_i + tk_heat(kk)
          END DO
          z_i = z_i - 0.5_wp * tk_heat(i)   ! node-center depth on the heat grid
        END BLOCK
        BLOCK
          INTEGER :: k
          IF (z_i <= tsoil_init_z(1)) THEN
            Tsoil_prof(i) = tsoil_init_T(1)
          ELSE IF (z_i >= tsoil_init_z(tsoil_init_n)) THEN
            Tsoil_prof(i) = tsoil_init_T(tsoil_init_n)   ! hold the deepest anchor below it
          ELSE
            DO k = 1, tsoil_init_n - 1
              IF (z_i <= tsoil_init_z(k+1)) THEN
                Tsoil_prof(i) = tsoil_init_T(k) + (tsoil_init_T(k+1) - tsoil_init_T(k)) * &
                  (z_i - tsoil_init_z(k)) / MAX(1.0e-9_wp, tsoil_init_z(k+1) - tsoil_init_z(k))
                EXIT
              END IF
            END DO
          END IF
        END BLOCK
      END DO
      WRITE (*,'(A,I0,A)') ' TSOIL-IC: sensor-anchored profile from tsoil_init.csv (', tsoil_init_n, ' anchors)'
      RETURN
    END IF
    SELECT CASE (soil_T_init_method)
    CASE (0)               ! uniform deep temperature
      DO i = 1, n_heat_tot
        Tsoil_prof(i) = Tdeep_C
      END DO
    CASE DEFAULT           ! 1 (and any other): damped exponential surface -> deep
      DO i = 1, n_heat_tot
        z_i = 0.0_wp
        IF (ALLOCATED(dz) .AND. i <= SIZE(dz)) z_i = ABS(dz(i))
        IF (i > nlayer) z_i = SUM(tk_heat(1:i)) - 0.5_wp*tk_heat(i)
        Tsoil_prof(i) = Tdeep_C + (Tsurf_C - Tdeep_C) * EXP(-z_i / z_damp_init)
      END DO
    END SELECT
    CALL update_soil_thermal_properties(nlayer)
  END SUBROUTINE initialize_soil_heat

  ! --------------------------------------------------------------------------- !
  ! Surface-face conductance and node-1 temperature of the held profile, so the
  ! linear conductive surface flux G = k_face*(Tsurf - T_node1) can be solved
  ! IMPLICITLY in the energy-balance repartition (an explicit Picard update of
  ! this term is unstable for a stiff k_face):
  !   k_face  = lambda_soil(1) / (0.5*ABS(tk(1)))   [W m-2 K-1]
  !   T_node1 = Tsoil_prof(1)                        [C]
  ! Only ever called under energy_balance_mode == 2 after initialize_soil_heat; a
  ! missing/undersized profile is a programming error, so it fails LOUDLY.
  SUBROUTINE soil_heat_surface_coeffs(k_face, T_node1)
    REAL(wp), INTENT(OUT) :: k_face, T_node1
    LOGICAL :: ok
    ok = ALLOCATED(Tsoil_prof) .AND. ALLOCATED(lambda_soil) .AND. ALLOCATED(tk)
    IF (ok) ok = (SIZE(Tsoil_prof) >= 1) .AND. (SIZE(lambda_soil) >= 1) .AND. (SIZE(tk) >= 1)
    IF (.NOT. ok) THEN
      WRITE (*,*) 'ERROR: soil_heat_surface_coeffs called before the 1-D soil-heat profile ' // &
                  'was initialized (mode 2). This is a bug in the EB call order.'
      STOP 1
    END IF
    k_face  = lambda_soil(1) / MAX(1.0e-6_wp, 0.5_wp * ABS(tk(1)))
    T_node1 = Tsoil_prof(1)
    ! When armed by the caller, the demand solve conducts against the committed
    ! node-1 trajectory of the previous outer iteration (override set in SHARED).
    IF (shf_node1_use) T_node1 = shf_node1_val
  END SUBROUTINE soil_heat_surface_coeffs

  ! --------------------------------------------------------------------------- !
  ! Conductive surface heat flux from the held profile (W m-2, downward-positive):
  !   G = k_face*(Tsurf_C - T_node1)  (single-sourced via soil_heat_surface_coeffs).
  ! G > 0 when Tsurf_C > T_node1 -> heat flows downward into the soil (same sign as
  ! advance_soil_heat). Used by the demand Newton solve; the repartition uses the
  ! coeffs directly for the implicit solve.
  FUNCTION soil_heat_surface_flux(Tsurf_C) RESULT(g)
    REAL(wp), INTENT(IN) :: Tsurf_C
    REAL(wp) :: g, k_face, T_node1
    CALL soil_heat_surface_coeffs(k_face, T_node1)
    g = k_face * (Tsurf_C - T_node1)
  END FUNCTION soil_heat_surface_flux

  ! --------------------------------------------------------------------------- !
  ! Backward-Euler (implicit) advance of the 1-D soil-temperature profile over
  ! dt_s seconds, finite-volume on the heat grid. z is positive DOWNWARD; all
  ! face fluxes are downward-positive (G > 0 = heat into the soil). The surface
  ! face takes the prescribed accepted flux G_surface_actual [W m-2]; the bottom
  ! face is ZERO-FLUX (G_bottom = 0; bottom_bc_type/Tdeep_C are NOT used here).
  ! Cv_soil/lambda_soil are refreshed from current SHARED theta at entry and held
  ! fixed across the solve, so the storage closure is exact to round-off:
  !   heat_resid = sum_i Cv_i*thick_i*(T_i^{n+1}-T_i^n) - dt_s*G_surface_actual.
  ! PRECONDITIONS: nlayer>=1, dt_s>0, and Tsoil_prof/Cv_soil/lambda_soil/theta/
  ! theta_0/theta_r/tk all ALLOCATED with SIZE >= nlayer. If any fails the routine
  ! does NOT modify Tsoil_prof and returns heat_resid = -HUGE (clearly invalid);
  ! it never silently solves a partial column.
  SUBROUTINE advance_soil_heat(nlayer, dt_s, G_surface_actual, heat_resid)
    INTEGER,  INTENT(IN)  :: nlayer
    REAL(wp), INTENT(IN)  :: dt_s, G_surface_actual
    REAL(wp), INTENT(OUT) :: heat_resid
    INTEGER  :: i, n
    LOGICAL  :: ok
    REAL(wp), ALLOCATABLE :: sub(:), dia(:), sup(:), rhs(:), Told(:), Tnew(:)
    REAL(wp), ALLOCATABLE :: gf(:), ccap(:), thick(:)
    heat_resid = 0.0_wp
    ok = (nlayer >= 1) .AND. (dt_s > 0.0_wp)
    ok = ok .AND. ALLOCATED(Tsoil_prof) .AND. ALLOCATED(Cv_soil) .AND. ALLOCATED(lambda_soil)
    ok = ok .AND. ALLOCATED(theta) .AND. ALLOCATED(theta_0) .AND. ALLOCATED(theta_r) &
              .AND. ALLOCATED(tk)
    IF (ok) THEN
      ok = (SIZE(Tsoil_prof) >= nlayer) .AND. (SIZE(Cv_soil)  >= nlayer) .AND. &
           (SIZE(lambda_soil) >= nlayer) .AND. (SIZE(theta)   >= nlayer) .AND. &
           (SIZE(theta_0)     >= nlayer) .AND. (SIZE(theta_r) >= nlayer) .AND. &
           (SIZE(tk)          >= nlayer)
    END IF
    IF (.NOT. ok) THEN
      heat_resid = -HUGE(1.0_wp)        ! preconditions unmet: no partial solve
      RETURN
    END IF
    n = MAX(nlayer, n_heat_tot)                      ! solve on the full heat grid
    CALL update_soil_thermal_properties(nlayer)      ! refresh from theta, then freeze
    ALLOCATE (sub(n), dia(n), sup(n), rhs(n), Told(n), Tnew(n), gf(n), ccap(n), thick(n))
    gf = 0.0_wp
    DO i = 1, n
      thick(i) = 1.0e-9_wp                          ! one positive effective thickness
      IF (i <= SIZE(tk)) thick(i) = MAX(1.0e-9_wp, ABS(tk(i)))
      IF (ALLOCATED(tk_heat)) THEN
        IF (i <= SIZE(tk_heat)) thick(i) = MAX(1.0e-9_wp, tk_heat(i))
      END IF
      Told(i)  = Tsoil_prof(i)
      ccap(i)  = Cv_soil(i) * thick(i) / dt_s       ! W m-2 K-1
    END DO
    ! Interface conductance: the harmonic (series) mean of the two half-layer
    ! resistances, so a conductivity contrast between layers is resolved.
    DO i = 1, n - 1
      gf(i) = 1.0_wp / ( 0.5_wp*thick(i)  /MAX(1.0e-9_wp, lambda_soil(i))   &
                       + 0.5_wp*thick(i+1)/MAX(1.0e-9_wp, lambda_soil(i+1)) )
    END DO
    DO i = 1, n
      sub(i) = 0.0_wp; sup(i) = 0.0_wp; dia(i) = ccap(i)
      IF (i > 1) THEN; sub(i) = -gf(i-1); dia(i) = dia(i) + gf(i-1); END IF
      IF (i < n) THEN; sup(i) = -gf(i);   dia(i) = dia(i) + gf(i);   END IF
      rhs(i) = ccap(i) * Told(i)
    END DO
    rhs(1) = rhs(1) + G_surface_actual              ! surface Neumann flux (downward +)
    CALL thomas_solve(n, sub, dia, sup, rhs, Tnew)
    heat_resid = -dt_s * G_surface_actual           ! G_bottom = 0
    DO i = 1, n
      Tsoil_prof(i) = Tnew(i)
      heat_resid = heat_resid + Cv_soil(i) * thick(i) * (Tnew(i) - Told(i))
    END DO
    DEALLOCATE (sub, dia, sup, rhs, Told, Tnew, gf, ccap, thick)
  END SUBROUTINE advance_soil_heat

  ! --------------------------------------------------------------------------- !
  ! Thomas algorithm for a tridiagonal system (sub/dia/sup = the three diagonals,
  ! rhs the right-hand side, x the solution). PRIVATE; local work arrays only, so
  ! it never touches the Richards SHARED scratch (r_a/r_b/r_c/r_r).
  SUBROUTINE thomas_solve(n, sub, dia, sup, rhs, x)
    INTEGER,  INTENT(IN)  :: n
    REAL(wp), INTENT(IN)  :: sub(n), dia(n), sup(n), rhs(n)
    REAL(wp), INTENT(OUT) :: x(n)
    INTEGER  :: i
    REAL(wp) :: m
    REAL(wp), ALLOCATABLE :: cp(:), dp(:)
    ALLOCATE (cp(n), dp(n))
    cp(1) = sup(1) / dia(1)
    dp(1) = rhs(1) / dia(1)
    DO i = 2, n
      m     = dia(i) - sub(i) * cp(i-1)
      cp(i) = sup(i) / m
      dp(i) = (rhs(i) - sub(i) * dp(i-1)) / m
    END DO
    x(n) = dp(n)
    DO i = n - 1, 1, -1
      x(i) = dp(i) - cp(i) * x(i+1)
    END DO
    DEALLOCATE (cp, dp)
  END SUBROUTINE thomas_solve

END MODULE soil_heat_module
