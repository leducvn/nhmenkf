! $Id: ropp_io_vlist_size.f90 3551 2013-02-25 09:51:28Z idculv $

!****f* Variable-lists/size *
!
! NAME
!    size - Number of variables in a variable list.
!
! SYNOPSIS
!    n = size(vlist)
!
! DESCRIPTION
!    This function returns the number of variables in a variable list, i.e.
!    the size of the variable list.
!
! INPUTS
!    type(Vlisttype)    :: vlist       A variable list.
!
!    or one of its sub-types:
!
!    type(VlisttypeD0d) :: vlist_0d
!    type(VlisttypeD1d) :: vlist_1d
!    type(VlisttypeD2d) :: vlist_2d
!
! OUTPUT
!    int                :: int
!
! NOTES
!    If a vlist is empty (i.e. has no variables stored in it), size returns 0.
!
!    In the present implementation, variable lists are limited to double
!    precision scalar, one- or two-dimensional floating point numbers. Internally,
!    these are implemented as linked lists, where seperate lists are kept for
!    the various data types and array shapes. The individual vlists are combined
!    into a single Vlisttype structure.
!
!    This function is overloading the intrinsic size() function for the various
!    (sub-) types of vlists.
!
! SEE ALSO
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
! 1. Size of a vlist
!-------------------------------------------------------------------------------

! 1.1 Generic
! -----------

FUNCTION size_vlist(vlist) RESULT(n)

! USE ropp_io,       not_this => size_vlist
  USE ropp_io_types, ONLY: Vlisttype

  IMPLICIT NONE

  TYPE(Vlisttype)       :: vlist
  INTEGER               :: n

  INTEGER, DIMENSION(3) :: nn

  INTEGER               :: size_vlistD0d
  INTEGER               :: size_vlistD1d
  INTEGER               :: size_vlistD2d

  nn(1) = size_vlistD0d(vlist%VlistD0d)
  nn(2) = size_vlistD1d(vlist%VlistD1d)
  nn(2) = size_vlistD2d(vlist%VlistD2d)

  n = SUM(nn)

END FUNCTION size_vlist

! 1.2 D0d
! --------

RECURSIVE FUNCTION size_vlistD0d(vlist) RESULT(n)

! USE ropp_io,       not_this => size_vlistD0d
  USE ropp_io_types, ONLY: VlisttypeD0d

  IMPLICIT NONE

  TYPE(VlisttypeD0d), POINTER :: vlist
  INTEGER                     :: n

  IF (ASSOCIATED(vlist)) THEN
     IF (ASSOCIATED(vlist%next)) THEN
        n = size_vlistD0d(vlist%next) + 1
     ELSE
        n = 1
     ENDIF
  ELSE
     n = 0
  ENDIF

END FUNCTION size_vlistD0d

! 1.3 D1d
! --------

RECURSIVE FUNCTION size_vlistD1d(vlist) RESULT(n)

! USE ropp_io, not_this => size_vlistD1d
  USE ropp_io_types, ONLY: VlisttypeD1d

  IMPLICIT NONE

  TYPE(VlisttypeD1d), POINTER :: vlist
  INTEGER                     :: n

  IF (ASSOCIATED(vlist)) THEN
     IF (ASSOCIATED(vlist%next)) THEN
        n = size_vlistD1d(vlist%next) + 1
     ELSE
        n = 1
     ENDIF
  ELSE
     n = 0
  ENDIF

END FUNCTION size_vlistD1d

! 1.4 D2d
! --------

RECURSIVE FUNCTION size_vlistD2d(vlist) RESULT(n)

! USE ropp_io, not_this => size_vlistD2d
  USE ropp_io_types, ONLY: VlisttypeD2d

  IMPLICIT NONE

  TYPE(VlisttypeD2d), POINTER :: vlist
  INTEGER                     :: n

  IF (ASSOCIATED(vlist)) THEN
     IF (ASSOCIATED(vlist%next)) THEN
        n = size_vlistD2d(vlist%next) + 1
     ELSE
        n = 1
     ENDIF
  ELSE
     n = 0
  ENDIF

END FUNCTION size_vlistD2d
