!===============================================================================!
! M&M-CO2Wheat v1.0  |  th.f90
!
! Purpose: Solar geometry of the day at the site: the solar declination from the
!          day of year, the sunset and sunrise solar time (hours) at the site
!          latitude, and the daylength as their difference. Called once per
!          simulated day, before the 24 hourly steps; spatial unit is the point
!          (one site latitude).
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: SHARED (DOY, lat_s_rad, Param_pi, zero, and the work variables dA, dB);
!       TYPES (wp)
! Output: SolarDecl (solar declination, rad), Sunrise (sunrise solar time, h)
!         and dIval (daylength, h), returned to CANOPY_MODEL.f90 as Solar_Decl,
!         dTsr and Daylength; meteo.f90 then takes Solar_Decl for the hourly
!         solar elevation and dTsr for the sunrise/sunset times of the day
! Called by: CANOPY_MODEL.f90 (CALL th (dTsr, Daylength, Solar_Decl))
!===============================================================================!
SUBROUTINE th (Sunrise, dIval, SolarDecl)
!------------------------------------------------------------------------------!
USE SHARED   ! Shared variables.
USE TYPES, ONLY: wp
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
IMPLICIT NONE
!------------------------------------------------------------------------------!
! Parameters.
!------------------------------------------------------------------------------!
REAL(wp)    :: SolarDecl
REAL(wp)    :: Sunset
REAL(wp)    :: Sunrise
REAL(wp)    :: dIval
! *****************************************************************************!
! Solar declination (radians)                                             (2.2).
! *****************************************************************************!
SolarDecl = -0.4093_wp * cos ( 2.0_wp * Param_pi * ((REAL(DOY, wp) + 10.0_wp)/365.0_wp))
! *****************************************************************************!
! Sunset solar time (hours), from the hour angle at which the solar elevation
! crosses zero at this declination and latitude                          (2.11).
! *****************************************************************************!
dA     = zero
dB     = zero
dA     = sin (SolarDecl) * sin (lat_s_rad)
dB     = cos (SolarDecl) * cos (lat_s_rad)
Sunset = 12.0_wp + (12.0_wp / Param_pi) * acos (-dA / dB)
! *****************************************************************************!
! Sunrise solar time (hour), symmetric to sunset about solar noon.
! *****************************************************************************!
Sunrise = (24.0_wp - Sunset)
! *****************************************************************************!
! Daylength (hours).
! *****************************************************************************!
dIval   = Sunset - Sunrise
!------------------------------------------------------------------------------!
END SUBROUTINE th
!===============================================================================!