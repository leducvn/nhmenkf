! $Id: ropp_io_vlist_print.f90 3551 2013-02-25 09:51:28Z idculv $

!****h|f|s|c|m|v|d|*|i* Module/Function *
!
! NAME
!    X - save the world
!
! SYNOPSIS
!
!
! DESCRIPTION
!    This ...
!
! ARGUMENTS
!    Use this, or the following two...
!
! INPUTS
!
!
! OUTPUT
!
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!
!
! CALLS
!
!
! REFERENCES
!
!
! TODO
!
!
! BUGS
!
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Print the contents of a vlist
!-------------------------------------------------------------------------------

! 1.1 Generic
! -----------

SUBROUTINE print_vlist(vlist)

! USE ropp_io,       not_this => print_vlist
  USE ropp_io_types, ONLY: Vlisttype

  IMPLICIT NONE

  TYPE(Vlisttype) :: vlist

  CALL print_vlistD0d(vlist%VlistD0d)
  CALL print_vlistD1d(vlist%VlistD1d)
  CALL print_vlistD2d(vlist%VlistD2d)

END SUBROUTINE print_vlist

! 1.2 D0d
! --------

RECURSIVE SUBROUTINE print_vlistD0d(vlist)

! USE ropp_io,       not_this => print_vlistD0d
  USE ropp_io_types, ONLY: VlisttypeD0d

  IMPLICIT NONE

  TYPE(VlisttypeD0d), POINTER :: vlist

  IF (ASSOCIATED(vlist)) THEN
     PRINT *, ""
     PRINT *, "Name:     ", TRIM(vlist%name)
     PRINT *, "Long name:", TRIM(vlist%long_name)
     PRINT *, "Units:    ", TRIM(vlist%units)
     PRINT *, "Value:    ", vlist%data
     IF (ASSOCIATED(vlist%next)) THEN
        CALL print_vlistD0d(vlist%next)
     ENDIF
  ENDIF

END SUBROUTINE print_vlistD0d

! 1.3 D1d
! --------

RECURSIVE SUBROUTINE print_vlistD1d(vlist)

! USE ropp_io, not_this => print_vlistD1d
  USE ropp_io_types, ONLY: VlisttypeD1d

  IMPLICIT NONE

  TYPE(VlisttypeD1d), POINTER :: vlist

  IF (ASSOCIATED(vlist)) THEN
     PRINT *, ""
     PRINT *, "Name:     ", TRIM(vlist%name)
     PRINT *, "Long name:", TRIM(vlist%long_name)
     PRINT *, "Units:    ", TRIM(vlist%units)
     PRINT *, "Value:    ", vlist%data(1:3), '...'
     IF (ASSOCIATED(vlist%next)) THEN
        CALL print_vlistD1d(vlist%next)
     ENDIF
  ENDIF

END SUBROUTINE print_vlistD1d


! 1.4 D2d
! --------

RECURSIVE SUBROUTINE print_vlistD2d(vlist)

! USE ropp_io, not_this => print_vlistD2d
  USE ropp_io_types, ONLY: VlisttypeD2d

  IMPLICIT NONE

  TYPE(VlisttypeD2d), POINTER :: vlist

  IF (ASSOCIATED(vlist)) THEN
     PRINT *, ""
     PRINT *, "Name:     ", TRIM(vlist%name)
     PRINT *, "Long name:", TRIM(vlist%long_name)
     PRINT *, "Units:    ", TRIM(vlist%units)
     PRINT *, "Value:    ", vlist%data(1:3,1), '...'
     IF (ASSOCIATED(vlist%next)) THEN
        CALL print_vlistD2d(vlist%next)
     ENDIF
  ENDIF

END SUBROUTINE print_vlistD2d
