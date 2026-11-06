!===============================================================================!
! M&M-CO2Wheat v1.0  |  csv_file_1d.f90
!
! Purpose: Include BODY, not a program unit: these statements are pulled into
!          csv_file.f90 to form the body of csv_write_dble_1d. They
!          write the elements of a one-dimensional array to the current CSV
!          record, advancing only after the last element. Per write call; no
!          time step and no spatial unit.
!
! Author: David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: the generic csv_write of csv_file (single-value routine)
! Output: one CSV record on the caller's logical unit
! Called by: nothing directly; included by csv_file.f90 inside
!            subroutine csv_write_dble_1d
!===============================================================================!
    integer, intent(in)                 :: lun
    logical, intent(in), optional       :: advance

    logical                             :: adv
    integer                             :: i

    adv = .true.
    if ( present(advance) ) adv = advance

    do i = 1,size(array)-1
        call csv_write( lun, array(i), .false. )
    enddo
    call csv_write( lun, array(size(array)), adv )
