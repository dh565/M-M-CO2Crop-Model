!===============================================================================!
! M&M-CO2Wheat v1.0  |  csv_file.f90
!
! Purpose: Minimal CSV writer used by the output routines: writes one
!          double-precision value, or a one-dimensional array of them, to the
!          current record of an already-open file, with or without advancing to
!          the next record. Reduced form of the public-domain csv_file module of
!          Arjen Markus (csv_file.f90,v 1.2, 2006), keeping only the
!          double-precision routines the model calls. Per write call; no time
!          step and no spatial unit.
!
! Author: David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: nothing; the body of the array routine is the include file
!       csv_file_1d.f90
! Output: formatted CSV records on the logical unit passed in by the caller.
!         main_MMCO2Crop.f90 writes the daily energy-balance diagnostic file
!         output/test_energy_bal_2025.csv through the generic csv_write
! Called by: main_MMCO2Crop.f90 (USE csv_file; CALL csv_write)
!
! File layout written by this module:
! - single items are written to the end of the current record;
! - one-dimensional items are also written to the end of the current record;
! - for single items the caller must state whether to advance to the next
!   record, for arrays the argument is optional and the default is to advance.
! - items are separated by a single comma; trailing blanks are dropped.
!   CSV conventions after Kernighan and Pike, The Practice of Programming,
!   Addison-Wesley, 1999.
! - the file must already be opened by the caller: only a unit number is passed.
!===============================================================================!
module csv_file
    implicit none
    interface csv_write
        module procedure csv_write_dble
        module procedure csv_write_dble_1d
    end interface
contains

! csv_write_dble --
!     Write a single double precision real to the CSV-file
! Arguments:
!     lun        LU-number of the CSV-file
!     value      Value to write
!     advance    Advance (.true.) or not, so that more items can be
!                written to the same record
! Result:
!     The value is written to the current record of the CSV-file
!
subroutine csv_write_dble( lun, value, advance )
    integer, intent(in)                    :: lun
    real(kind=kind(1.0d0)), intent(in)     :: value
    logical, intent(in)                    :: advance

    character(len=40)     :: buffer
    write( buffer, '(G20.12)' ) value
    buffer = adjustl(buffer)
    if ( advance ) then
       write(lun,'(a)') trim(buffer)
    else
       write(lun,'(a,a)',advance='no') trim(buffer), ','
    endif
end subroutine csv_write_dble

! csv_write_dble_1d --
!     Write a one-dimensional array of items to the CSV-file
! Arguments:
!     lun        LU-number of the CSV-file
!     array      Array to write
!     advance    Advance (.true.) or not, so that more items can be
!                written to the same record
! Result:
!     The array is written to the current record of the CSV-file
! Note:
!     The body is held in the include file csv_file_1d.f90.
!
subroutine csv_write_dble_1d( lun, array, advance )
    real(kind=kind(1.0d0)), dimension(:), intent(in)   :: array

    include 'csv_file_1d.f90'

end subroutine csv_write_dble_1d

end module csv_file
