! $Id: ropp_io_addvar.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* ExtendingDatatypes/ropp_io_addvar *
!
! NAME
!    ropp_io_addvar - Add an additional variable to an RO profile structure.
!
! SYNOPSIS
!    call ropp_io_addvar(ro_data, name, long_name, units, range, data)
!
! DESCRIPTION
!    This subroutine extends a given RO profile structure with a variable as
!    characterised by the arguments. Upon writing the variable to a netCDF
!    file, additional (compared to the standard set of netCDF) variables will
!    be created, using name, long_name, units and range as standard attributes
!    in the netCDF data structure.
!
! INPUTS
!    type(ROprof)     :: ro_data
!    char(len = *)    :: name
!    char(len = *)    :: long_name
!    char(len = *)    :: units
!    real(wp), dim(2) :: range
!    real(wp), dim(n) :: data
!
! OUTPUT
!    None; however, an additional variable with the attributes and data as
!    specified via the subroutine's arguments is attached to ro_data.
!
! NOTES
!    All arguments are mandatory.
!
! EXAMPLE
!    Adding the value of a cost function to a retrieved profile:
!
!       call ropp_io_addvar(ro_data, &
!                           name      = "J", &
!                           long_name = "Cost function value at convergence", &
!                           units     = "", &
!                           range     = (/ 0.0_wp, 9999.0_wp/), &
!                           data      = J)
!
!    A variable holding relative humidity as additional retrieved quantity
!    can be created with
!
!       call ropp_io_addvar(ro_data, &
!                           name      = "rhum", &
!                           long_name = "Relative humidity", &
!                           units     = "percent", &
!                           range     = (/ 0.0_wp, 9999.0_wp/), &
!                           data      = rhum)
!
!    where rhum is a one-dimensional double precision array.
!
! SEE ALSO
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
! 1. Add variables to an ROdata object
!-------------------------------------------------------------------------------

! 1.1 D0d
! -------

SUBROUTINE ropp_io_addvar_rodataD0d(ROdata, name, long_name, units, range, DATA)

! 1.1.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_rodataD0d
  USE ropp_io_types, ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof),           INTENT(inout) :: ROdata
  CHARACTER(len = *),     INTENT(in)    :: name
  CHARACTER(len = *),     INTENT(in)    :: long_name
  CHARACTER(len = *),     INTENT(in)    :: units
  REAL(wp), DIMENSION(2), INTENT(in)    :: range
  REAL(wp),               INTENT(in)    :: DATA

! 1.1.2 Call the appropriate recursive vlist routine

  CALL ropp_io_addvar_vlistD0d(ROdata%vlist%VlistD0d, name, long_name, units, range, DATA)

END SUBROUTINE ropp_io_addvar_rodataD0d

! 1.2 D1d
! --------

SUBROUTINE ropp_io_addvar_rodataD1d(ROdata, name, long_name, units, range, DATA)

! 1.2.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_rodataD1d
  USE ropp_io_types, ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof),           INTENT(inout) :: ROdata
  CHARACTER(len = *),     INTENT(in)    :: name
  CHARACTER(len = *),     INTENT(in)    :: long_name
  CHARACTER(len = *),     INTENT(in)    :: units
  REAL(wp), DIMENSION(2), INTENT(in)    :: range
  REAL(wp), DIMENSION(:), INTENT(in)    :: DATA

! 1.2.2 Call the appropriate recursive vlist routine

  CALL ropp_io_addvar_vlistD1d(ROdata%vlist%VlistD1d, name, long_name, units, range, DATA)

END SUBROUTINE ropp_io_addvar_rodataD1d


! 1.3 D2d
! --------

SUBROUTINE ropp_io_addvar_rodataD2d(ROdata, name, long_name, units, range, DATA)

! 1.3.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_rodataD2d
  USE ropp_io_types, ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof),             INTENT(inout) :: ROdata
  CHARACTER(len = *),       INTENT(in)    :: name
  CHARACTER(len = *),       INTENT(in)    :: long_name
  CHARACTER(len = *),       INTENT(in)    :: units
  REAL(wp), DIMENSION(2),   INTENT(in)    :: range
  REAL(wp), DIMENSION(:,:), INTENT(in)    :: DATA

! 1.2.2 Call the appropriate recursive vlist routine

  CALL ropp_io_addvar_vlistD2d(ROdata%vlist%VlistD2d, name, long_name, units, range, DATA)

END SUBROUTINE ropp_io_addvar_rodataD2d


!-------------------------------------------------------------------------------
! 2. Add variables to a Vlist
!-------------------------------------------------------------------------------

! 2.1 D0d
! -------

RECURSIVE SUBROUTINE ropp_io_addvar_vlistD0d(vlist, name, long_name, units, range, DATA)

! 2.1.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_vlistD0d
  USE ropp_io_types, ONLY: VlisttypeD0d

  IMPLICIT NONE

  TYPE(VlisttypeD0d),     POINTER    :: vlist
  CHARACTER(len = *),     INTENT(in) :: name
  CHARACTER(len = *),     INTENT(in) :: long_name
  CHARACTER(len = *),     INTENT(in) :: units
  REAL(wp), DIMENSION(2), INTENT(in) :: range
  REAL(wp),               INTENT(in) :: DATA

! 2.1.2 Step down to the next list element...

  IF (ASSOCIATED(vlist)) THEN

     CALL ropp_io_addvar_vlistD0d(vlist%next, name, long_name, units, range, DATA)

  ELSE

! 2.1.3 ...or fill the variable

     ALLOCATE(vlist)

     vlist%name      = name
     vlist%long_name = long_name
     vlist%units     = units
     vlist%range     = range
     vlist%data      = DATA

  ENDIF

END SUBROUTINE ropp_io_addvar_vlistD0d


! 2.2 D1d
! -------

RECURSIVE SUBROUTINE ropp_io_addvar_vlistD1d(vlist, name, long_name, units, range, DATA)

! 2.2.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_vlistD1d
  USE ropp_io_types, ONLY: VlisttypeD1d

  IMPLICIT NONE

  TYPE(VlisttypeD1d),     POINTER    :: vlist
  CHARACTER(len = *),     INTENT(in) :: name
  CHARACTER(len = *),     INTENT(in) :: long_name
  CHARACTER(len = *),     INTENT(in) :: units
  REAL(wp), DIMENSION(2), INTENT(in) :: range
  REAL(wp), DIMENSION(:), INTENT(in) :: DATA

! 2.2.2 Step down to the next list element...

  IF (ASSOCIATED(vlist)) THEN

     CALL ropp_io_addvar_vlistD1d(vlist%next, name, long_name, units, range, DATA)

  ELSE

! 2.2.3 ...or fill the variable

     ALLOCATE(vlist)
     ALLOCATE(vlist%data(SIZE(DATA)))

     vlist%name      = name
     vlist%long_name = long_name
     vlist%units     = units
     vlist%range     = range
     vlist%data      = DATA

  ENDIF

END SUBROUTINE ropp_io_addvar_vlistD1d


! 2.3 D2d
! -------

RECURSIVE SUBROUTINE ropp_io_addvar_vlistD2d(vlist, name, long_name, units, range, DATA)

! 2.3.1 Declarations

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_addvar_vlistD2d
  USE ropp_io_types, ONLY: VlisttypeD2d

  IMPLICIT NONE

  TYPE(VlisttypeD2d),       POINTER    :: vlist
  CHARACTER(len = *),       INTENT(in) :: name
  CHARACTER(len = *),       INTENT(in) :: long_name
  CHARACTER(len = *),       INTENT(in) :: units
  REAL(wp), DIMENSION(2),   INTENT(in) :: range
  REAL(wp), DIMENSION(:,:), INTENT(in) :: DATA

! 2.3.2 Step down to the next list element...

  IF (ASSOCIATED(vlist)) THEN

     CALL ropp_io_addvar_vlistD2d(vlist%next, name, long_name, units, range, DATA)

  ELSE

! 2.3.3 ...or fill the variable

     ALLOCATE(vlist)
     ALLOCATE(vlist%data(SIZE(DATA,1), SIZE(DATA,2)))

     vlist%name      = name
     vlist%long_name = long_name
     vlist%units     = units
     vlist%range     = range
     vlist%data      = DATA

  ENDIF

END SUBROUTINE ropp_io_addvar_vlistD2d
