!===============================================================================!
! M&M-CO2Wheat v1.0  |  count_lines.f90
!
! Purpose: Counts the records of the input file currently named in file_name
!          (SHARED) so that the main program can size its daily arrays before
!          reading them. Called once per file at initialization, before the
!          simulation starts; no time step and no spatial unit.
!
! Author: David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: SHARED (file_name - the path of the file to be counted)
! Output: nlines - the number of records read, including the header line;
!         main_MMCO2Crop.f90 subtracts the header to obtain NDays and then
!         allocates all daily arrays with it
! Called by: main_MMCO2Crop.f90 (CALL count_lines(NDays))
!===============================================================================!
SUBROUTINE count_lines(nlines)

!------------------------------------------------------------------------------!
USE SHARED ! Shared variables.
!------------------------------------------------------------------------------!

IMPLICIT NONE
!------------------------------------------------------------------------------!
INTEGER :: nlines
!------------------------------------------------------------------------------!
nlines = 0
OPEN(1,file = file_name)
!------------------------------------------------------------------------------!
! Count records until the end-of-file condition is reached.
!------------------------------------------------------------------------------!
DO ! Loop over lines and count
  READ(1,*, END=10)
  nlines = nlines + 1
END DO
10 CLOSE (1)
!------------------------------------------------------------------------------!
End SUBROUTINE count_lines
!===============================================================================!
