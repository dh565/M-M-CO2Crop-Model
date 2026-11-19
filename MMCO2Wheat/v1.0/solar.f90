!===============================================================================!
! M&M-CO2Wheat v1.0  |  solar.f90
!
! Purpose: Canopy radiation geometry for the two-big-leaf photosynthesis scheme.
!          From the solar elevation and the direct/diffuse split of global
!          radiation it computes the beam and diffuse extinction coefficients,
!          the penetration of direct and diffuse photosynthetically active
!          radiation through the canopy, the PAR absorbed per unit sunlit and per
!          unit shaded leaf area, and the sunlit and shaded leaf area indices.
!          Time step: hourly. Spatial unit: point (whole canopy, two leaf
!          classes).
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: SHARED (solar elevation SolarElev_s, direct and diffuse global radiation
!       drrad_s/dfrad_s, green lamina area index gParam_lai and lai_s, the
!       reflection coefficients const_P_sr/const_P_par, the scattering and
!       absorptance coefficients const_alpha_sr/const_alpha_par, the PAR
!       reflection-basis switch par_reflection_basis, zero); TYPES (working
!       precision wp).
! Output: Kdr, the beam extinction coefficient (consumed by radiation_module.f90
!       for the shortwave split, by energy_balance_module.f90 for the
!       net-radiation partition, and by CANOPY_MODEL.f90 for the sunlit-weighted
!       canopy nitrogen profile); Absor_IPAR_sunlit_h and Absor_IPAR_shaded_h,
!       the absorbed PAR per unit sunlit and shaded leaf area (W m-2 leaf),
!       converted to photon flux and passed to assim.f90 by CANOPY_MODEL.f90;
!       sunlit_lai and shaded_lai, the sunlit and shaded leaf area indices, which
!       scale leaf fluxes to the canopy in CANOPY_MODEL.f90; I_drrad_h and
!       I_dfrad_h, the canopy-intercepted direct and diffuse total solar
!       radiation (W m-2 ground), returned as diagnostics with no downstream
!       consumer.
! Called by: CANOPY_MODEL.f90, once per hourly step before the canopy-temperature
!       and assimilation solution.
!===============================================================================!
SUBROUTINE solar (Kdr,I_drrad_h, I_dfrad_h, Absor_IPAR_sunlit_h, Absor_IPAR_shaded_h,sunlit_lai,shaded_lai)
!------------------------------------------------------------------------------!
! The canopy is treated as two big leaves, sunlit and shaded, with a spherical
! leaf-angle distribution (Goudriaan, 1977; de Pury and Farquhar, 1997); the
! scattered-beam term follows Goudriaan (1977) and Spitters (1986). The
! implemented algebra is that of Teh (2006), Introduction to Mathematical
! Modeling of Crop Growth; the equation numbers quoted below are its own.
!------------------------------------------------------------------------------!
USE SHARED   ! Shared variables.
USE TYPES, ONLY: wp
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
IMPLICIT NONE
!------------------------------------------------------------------------------!
! Parameters.
!------------------------------------------------------------------------------!
REAL(wp)    :: Kdr
REAL(wp)    :: Kdf
REAL(wp)    :: dElev
REAL(wp)    :: dL
REAL(wp)    :: dN
REAL(wp)    :: tao_drrad
REAL(wp)    :: tao_drrad_alpha
REAL(wp)    :: IPAR_pdr
REAL(wp)    :: IPAR_pdrdr

REAL(wp)    :: I_drrad_h
REAL(wp)    :: I_dfrad_h
REAL(wp)    :: IPAR_dr
REAL(wp)    :: IPAR_df
REAL(wp)    :: IPAR_pdf
REAL(wp)    :: IPAR_pdra
REAL(wp)    :: Absor_IPAR_sunlit_h
REAL(wp)    :: Absor_IPAR_shaded_h
REAL(wp)    :: sunlit_lai
REAL(wp)    :: shaded_lai
!------------------------------------------------------------------------------!

! *****************************************************************************!
! 1. Intercepted radiation.
! *****************************************************************************!
!------------------------------------------------------------------------------!
! Extinction coefficient for direct (beam) fluxes, Kdr = 0.5/sin(beta) for a
! spherical leaf-angle distribution (Teh 2006, Eq. 3.3). With the sun below the
! horizon there is no beam, the whole canopy is shaded, and Kdr = 0.
!------------------------------------------------------------------------------!
sunlit_lai = zero
shaded_lai = gParam_lai
dElev = SolarElev_s
IF (dElev < 1.0e-8_wp) THEN    !    sun below the horizon
  Kdr = 0.0_wp
ELSE
  Kdr = 0.5_wp / sin (dElev)
END IF
! ---------------------------------------------------------------------------- !
! Extinction coefficient for diffuse fluxes, a function of leaf area index
! (Teh 2006, Eq. 3.29).
!------------------------------------------------------------------------------!
dL     = sqrt (lai_s)
Kdf    = (1. + 0.1174 * dL) / (1. + 0.3732 * dL)
! ---------------------------------------------------------------------------- !
! Hourly total solar radiation intercepted by the canopy (W m-2 ground), direct
! and diffuse, each as incident minus canopy reflection times one minus the Beer
! transmission with the scattering-corrected extinction (Teh 2006, Eqs. 3.8 and
! 3.32); const_alpha_sr = 0.5 is the leaf scattering coefficient for total solar
! radiation and const_P_sr = 0.11 its canopy reflection coefficient.
!------------------------------------------------------------------------------!
! direct:
I_drrad_h = (1. - const_P_sr) * drrad_s &
        * (1. - exp (-Kdr * sqrt(const_alpha_sr) * gParam_lai))
! diffuse:
I_dfrad_h = (1. - const_P_sr) * dfrad_s &
        * (1. - exp (-Kdf * sqrt(const_alpha_sr) * gParam_lai))
! ---------------------------------------------------------------------------- !
! Penetration functions for the direct beam within the canopy.
!------------------------------------------------------------------------------!
! transmission of beam plus its scattered component:
tao_drrad       = exp ( -Kdr * gParam_lai )                         ! Teh 2006, Eq. 3.1
! transmission of the unintercepted beam only (no scattering):
tao_drrad_alpha = exp ( -Kdr * gParam_lai * sqrt(const_alpha_par) )  ! Teh 2006, Eq. 3.6
!------------------------------------------------------------------------------!
! Flux components within the canopy.
!------------------------------------------------------------------------------!
! PAR from direct and diffuse radiation (W m-2). The PAR fraction of global
! radiation is 0.45.
! --------------------------------------------------- !
! A. Direct-beam PAR incident at canopy top (W m-2 leaf area; Teh 2006, Eq. 3.38).
! This is the incident flux density (drrad_s), not the canopy-intercepted total
! I_drrad_h: the latter carries a (1-exp(-k*LAI)) factor that makes it
! per-ground-area, whereas the formula requires the per-leaf incident flux so
! that Kdr * IPAR_dr gives absorbed PAR per unit sunlit leaf area.
! One reflection coefficient is applied per waveband (Teh 2006, p.67-69,
! Eqs. 3.33-3.39, with the worked example on p.139), so for the PAR stream the
! only reflection subtraction is the PAR coefficient const_P_par = 0.04, applied
! below at IPAR_pdf/IPAR_pdr/IPAR_pdrdr. The alternative branch below additionally
! subtracts the total-solar coefficient const_P_sr = 0.11 from the incident PAR
! flux; the 0.45 PAR fraction is unaffected either way.
IF (par_reflection_basis >= 1) THEN
  IPAR_dr   = 0.45 * drrad_s
ELSE
  IPAR_dr   = 0.45 * (1. - const_P_sr) * drrad_s
END IF
! --------------------------------------------------- !
! B. Diffuse PAR incident at canopy top (W m-2 leaf area; Teh 2006, Eqs. 3.34 and
! 3.39), and its mean over the canopy depth. Same reasoning as A: the incident
! flux, not the canopy-intercepted total, and one reflection coefficient per
! waveband.
IF (par_reflection_basis >= 1) THEN
  IPAR_df   = 0.45 * dfrad_s
ELSE
  IPAR_df   = 0.45 * (1. - const_P_sr) * dfrad_s
END IF
dN          = Kdf * sqrt(const_alpha_par) * gParam_lai
IF (ABS(dN) > 1.0e-9) THEN
  IPAR_pdf    = ( (1. - const_P_par) * IPAR_df * (1. - exp (-dN)) ) / dN
ELSE
  IPAR_pdf    = zero
END IF
! --------------------------------------------------- !
! C. Direct-beam PAR with and without its scattered component, and the scattered
! beam itself as half their difference (Teh 2006, Eqs. 3.33, 3.34 and 3.36).
! --------------------------------------------------- !
IPAR_pdr    = (1. - const_P_par) * tao_drrad_alpha * IPAR_dr
IPAR_pdrdr  = (1. - const_P_par) * tao_drrad * IPAR_dr
IPAR_pdra   = (IPAR_pdr - IPAR_pdrdr) / 2.
!------------------------------------------------------------------------------!
! Hourly absorbed PAR per unit sunlit and per unit shaded leaf area (W m-2 leaf);
! const_alpha_par = 0.8 is the leaf PAR absorptance.
!------------------------------------------------------------------------------!
! sunlit leaves receive beam, diffuse and scattered beam (Teh 2006, Eq. 3.38):
Absor_IPAR_sunlit_h &
              = const_alpha_par * (Kdr * IPAR_dr + IPAR_pdf + IPAR_pdra)
! shaded leaves receive diffuse and scattered beam only (Teh 2006, Eq. 3.39):
Absor_IPAR_shaded_h &
              = const_alpha_par * (IPAR_pdf + IPAR_pdra)

! *****************************************************************************!
! 2. Leaf area index of sunlit and shaded leaves (Teh 2006, Eqs. 3.41 and 3.42).
!    The sunlit fraction is the Beer integral (1-exp(-Kdr*LAI))/Kdr; the rest of
!    the green lamina area is shaded. With no beam both keep the values set above.
! *****************************************************************************!
IF (Kdr > zero) THEN
  sunlit_lai = (1. - tao_drrad) / Kdr
  shaded_lai = gParam_lai - sunlit_lai
END IF
!------------------------------------------------------------------------------!
! 3. Driver record 185 = 1: absorbed PAR that conserves the intercepted PAR.
!    Per unit leaf area at cumulative leaf area l, the absorbed flux is the depth
!    derivative of the attenuated fluxes of Eqs. 3.31 and 3.33-3.34 (same kdr,
!    kdf, sqrt(alpha) and p):
!      diffuse          (1-p) kd' Q_df exp(-kd' l),             kd' = sqrt(alpha) kdf
!      beam + scatter   (1-p) kb' Q_dr exp(-kb' l),             kb' = sqrt(alpha) kdr
!      direct beam      alpha kdr Q_dr exp(-kdr l)  (sunlit leaves only)
!    The sunlit class integrates these weighted by the sunlit probability
!    exp(-kdr l); the shaded class is the canopy total minus the sunlit part.
!    The canopy totals are then exactly the intercepted fluxes of Eqs. 3.8 and
!    3.32, which the default branch (Eqs. 3.36-3.39) does not conserve.
!    Sources (read in full): absorbed irradiance per
!    unit leaf area is the depth derivative of the flux, Spitters (1986) Agric
!    For Meteorol 38:231, Eqs. 9-14 pp. 234-235, and de Pury & Farquhar (1997)
!    Plant Cell Environ 20:537, Table A1 Eqs. A1-A8 and A11 p. 553; the sunlit
!    integrals below are the complements of their shaded-leaf integrals,
!    Eqs. A26b (diffuse) and A26c (scattered beam), p. 555, term for term,
!    with this model's own coefficients (alpha = 1 - sigma, p for rho).
!------------------------------------------------------------------------------!
IF (par_absorb_conserving == 1 .AND. gParam_lai > 1.0e-6_wp) THEN
  BLOCK
    REAL(wp) :: kdp, kbp, sa, L, Ic_df, Ic_bt, Ic_bd, Isun
    sa  = sqrt(const_alpha_par)
    L   = gParam_lai
    kdp = Kdf * sa
    Ic_df = (1. - const_P_par) * IPAR_df * (1. - exp(-kdp * L))
    IF (Kdr > zero) THEN
      kbp   = Kdr * sa
      Ic_bt = (1. - const_P_par) * IPAR_dr * (1. - exp(-kbp * L))
      Ic_bd = const_alpha_par * IPAR_dr * (1. - tao_drrad)
      Isun  = Ic_bd &
            + (1. - const_P_par) * IPAR_df * kdp / (kdp + Kdr) * (1. - exp(-(kdp + Kdr) * L)) &
            + IPAR_dr * ( (1. - const_P_par) * kbp / (kbp + Kdr) * (1. - exp(-(kbp + Kdr) * L)) &
                        - const_alpha_par * (1. - exp(-2. * Kdr * L)) / 2. )
      Absor_IPAR_sunlit_h = Isun / MAX(sunlit_lai, 1.0e-9_wp)
      Absor_IPAR_shaded_h = MAX(zero, Ic_df + Ic_bt - Isun) / MAX(shaded_lai, 1.0e-9_wp)
    ELSE
      Absor_IPAR_sunlit_h = zero
      Absor_IPAR_shaded_h = Ic_df / L
    END IF
  END BLOCK
END IF
!------------------------------------------------------------------------------!
END SUBROUTINE solar
!==============================================================================!