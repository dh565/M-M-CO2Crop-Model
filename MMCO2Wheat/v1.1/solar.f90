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
!       consumer in v1.0.
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
END SUBROUTINE solar
!==============================================================================!