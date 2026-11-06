!===============================================================================!
! M&M-CO2Wheat v1.0  |  TYPES.f90
!
! Purpose: Declares the working-precision kind parameters of the model: wp for
!          double-precision reals and ik for integers. Declaration-only module;
!          it carries no state, no time step and no spatial unit.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: nothing (intrinsic SELECTED_REAL_KIND / SELECTED_INT_KIND only)
! Output: wp, ik - the kind parameters with which every model variable is
!         declared (SHARED, main_MMCO2Crop, assim, CANOPY_MODEL, growth,
!         phenology, hydro, et, meteo, radiation, soil_heat, soiln, mineral,
!         energy_balance, solar, th)
! Called by: every module (kind parameters only; no executable statements)
!===============================================================================!
MODULE TYPES

  IMPLICIT NONE

  ! Standard floating point precision (double precision).
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(15, 307)

  ! Standard integer kind.
  INTEGER, PARAMETER :: ik = SELECTED_INT_KIND(9)

END MODULE TYPES
!===============================================================================!
