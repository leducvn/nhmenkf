! $Id: ropp_io.f90 3551 2013-02-25 09:51:28Z idculv $

MODULE ropp_io

!****m* Modules/ropp_io *
!
! NAME
!    ropp_io - The ROPP input / output library.
!
! SYNOPSIS
!    use ropp_io
!
! DESCRIPTION
!    This Fortran module provides interfaces and data types required for the
!    use of the ROPP input / output library.
!
! NOTES
!    For reading data from ROPP data files, the ropp_io_read routine is the
!    primary interface to the library. This routine will return a Fortran
!    derived type - or structure - "ROprof" which, on return, holds the
!    contents of the data file. Writing data to a ROPP data file requires
!    filling an instance of a ROprof data structure with the data to written.
!    The arrays contained in this structure might be initialised by a call to
!    ropp_io_init. The actual writing to a data file is performed using
!    the subroutine ropp_io_write.
!
!    The data format is based on the netCDF scientific data format.
!
! EXAMPLE
!    Assume that a data file shall be read into a Fortran program; the
!    required steps are
!
!       use ropp_io
!         ...
!       type(ROprof)       :: ro_data
!         ...
!       character(len = *) :: file = 'my_data_file.nc'
!         ...
!       call ropp_io_read(ro_data, file)
!
!    Similarly, writing a RO profile contains of the following steps:
!
!       use ropp_io
!         ...
!       type(ROprof)       :: ro_data
!         ...
!       character(len = *) :: file = 'my_data_file.nc'
!         ...
!
!    !  Fill the structure
!    !  ------------------
!
!       call ropp_io_init(ro_data, NLev1a, NLev1b, NLev2a, NLev2b, NLev2c, NLev2d)
!         ...
!
!    !  Fill the structure
!    !  ------------------
!
!       ro_data%Lev1b%impact_L1 = ...
!       ro_data%Lev1b%bangle_L1 = ...
!
!    !  Write the data
!    !  --------------
!
!       call ropp_io_write(ro_data, file)
!
! SEE ALSO
!    ropp_io_types
!    ropp_io_read
!    ropp_io_write
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
! 1. Include derived data types
!-------------------------------------------------------------------------------

  USE ropp_io_types

!-------------------------------------------------------------------------------
! 2. Public interfaces
!-------------------------------------------------------------------------------

! 2.1 Initialisation
! ------------------

  INTERFACE ropp_io_init
     SUBROUTINE ropp_io_init_l1atype(var, n)
       USE ropp_io_types
       TYPE(L1atype), INTENT(inout) :: var
       INTEGER                      :: n
     END SUBROUTINE ropp_io_init_l1atype
     SUBROUTINE ropp_io_init_l1btype(var, n)
       USE ropp_io_types
       TYPE(L1btype), INTENT(inout) :: var
       INTEGER                      :: n
     END SUBROUTINE ropp_io_init_l1btype
     SUBROUTINE ropp_io_init_l2atype(var, n)
       USE ropp_io_types
       TYPE(L2atype), INTENT(inout) :: var
       INTEGER                     :: n
     END SUBROUTINE ropp_io_init_l2atype
     SUBROUTINE ropp_io_init_l2btype(var, n)
       USE ropp_io_types
       TYPE(L2btype), INTENT(inout) :: var
       INTEGER                     :: n
     END SUBROUTINE ropp_io_init_l2btype
     SUBROUTINE ropp_io_init_l2btype_2d(var, n)
       USE ropp_io_types
       TYPE(L2btype_2d), INTENT(inout) :: var
       INTEGER, DIMENSION(2)           :: n
     END SUBROUTINE ropp_io_init_l2btype_2d
     SUBROUTINE ropp_io_init_l2ctype(var, n)
       USE ropp_io_types
       TYPE(L2ctype), INTENT(inout) :: var
       INTEGER                      :: n
     END SUBROUTINE ropp_io_init_l2ctype
     SUBROUTINE ropp_io_init_l2ctype_2d(var, n)
       USE ropp_io_types
       TYPE(L2ctype_2d), INTENT(inout) :: var
       INTEGER, DIMENSION(2)           :: n
     END SUBROUTINE ropp_io_init_l2ctype_2d
     SUBROUTINE ropp_io_init_l2dtype(var, n)
       USE ropp_io_types
       TYPE(L2dtype), INTENT(inout) :: var
       INTEGER                      :: n
     END SUBROUTINE ropp_io_init_l2dtype
     SUBROUTINE ropp_io_init_vlist(var)
       USE ropp_io_types
       TYPE(Vlisttype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_init_vlist
     SUBROUTINE ropp_io_init_rotype(ROdata, NLev1a, NLev1b, &
                                            NLev2a, NLev2b, NLev2c, NLev2d)
       USE ropp_io_types
       TYPE(ROprof), INTENT(inout) :: ROdata
       INTEGER,      INTENT(in)    :: NLev1a
       INTEGER,      INTENT(in)    :: NLev1b
       INTEGER,      INTENT(in)    :: NLev2a
       INTEGER,      INTENT(in)    :: NLev2b
       INTEGER,      INTENT(in)    :: NLev2c
       INTEGER,      INTENT(in)    :: NLev2d
     END SUBROUTINE ropp_io_init_rotype
     SUBROUTINE ropp_io_init_rotype_2d(ROdata, NLev1a, NLev1b, NLev2a, &
                                       NLev2b, NLev2c, NLev2d, NHoriz)
       USE ropp_io_types
       TYPE(ROprof2d), INTENT(inout) :: ROdata
       INTEGER,        INTENT(in)    :: NLev1a
       INTEGER,        INTENT(in)    :: NLev1b
       INTEGER,        INTENT(in)    :: NLev2a
       INTEGER,        INTENT(in)    :: NLev2b
       INTEGER,        INTENT(in)    :: NLev2c
       INTEGER,        INTENT(in)    :: NLev2d
       INTEGER,        INTENT(in)    :: NHoriz
     END SUBROUTINE ropp_io_init_rotype_2d
  END INTERFACE

! 2.2a Structure copy ROprof -> ROprof
! ------------------------------------

  INTERFACE ropp_io_roprof2roprof
     SUBROUTINE ropp_io_l1atype_l1atype(from_var, to_var)
       USE ropp_io_types
       TYPE(L1atype), INTENT(in)    :: from_var
       TYPE(L1atype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l1atype_l1atype
     SUBROUTINE ropp_io_l1btype_l1btype(from_var, to_var)
       USE ropp_io_types
       TYPE(L1btype), INTENT(in)    :: from_var
       TYPE(L1btype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l1btype_l1btype
     SUBROUTINE ropp_io_l2atype_l2atype(from_var, to_var)
       USE ropp_io_types
       TYPE(L2atype), INTENT(in)    :: from_var
       TYPE(L2atype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l2atype_l2atype
     SUBROUTINE ropp_io_l2btype_l2btype(from_var, to_var)
       USE ropp_io_types
       TYPE(L2btype), INTENT(in)    :: from_var
       TYPE(L2btype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l2btype_l2btype
     SUBROUTINE ropp_io_l2ctype_l2ctype(from_var, to_var)
       USE ropp_io_types
       TYPE(L2ctype), INTENT(in)    :: from_var
       TYPE(L2ctype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l2ctype_l2ctype
     SUBROUTINE ropp_io_l2dtype_l2dtype(from_var, to_var)
       USE ropp_io_types
       TYPE(L2dtype), INTENT(in)    :: from_var
       TYPE(L2dtype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_l2dtype_l2dtype
     SUBROUTINE ropp_io_rotype_rotype(from_var, to_var)
       USE ropp_io_types
       TYPE(ROprof), INTENT(in)    :: from_var
       TYPE(ROprof), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_rotype_rotype
  END INTERFACE

! 2.2b Structure copy ROprof -> ROprof (assignment)
! -------------------------------------------------

  INTERFACE ASSIGNMENT(=)
     SUBROUTINE ropp_io_assign_l1atype(to_var, from_var)
       USE ropp_io_types
       TYPE(L1atype), INTENT(in)    :: from_var
       TYPE(L1atype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l1atype
     SUBROUTINE ropp_io_assign_l1btype(to_var, from_var)
       USE ropp_io_types
       TYPE(L1btype), INTENT(in)    :: from_var
       TYPE(L1btype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l1btype
     SUBROUTINE ropp_io_assign_l2atype(to_var, from_var)
       USE ropp_io_types
       TYPE(L2atype), INTENT(in)    :: from_var
       TYPE(L2atype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l2atype
     SUBROUTINE ropp_io_assign_l2btype(to_var, from_var)
       USE ropp_io_types
       TYPE(L2btype), INTENT(in)    :: from_var
       TYPE(L2btype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l2btype
     SUBROUTINE ropp_io_assign_l2ctype(to_var, from_var)
       USE ropp_io_types
       TYPE(L2ctype), INTENT(in)    :: from_var
       TYPE(L2ctype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l2ctype
     SUBROUTINE ropp_io_assign_l2dtype(to_var, from_var)
       USE ropp_io_types
       TYPE(L2dtype), INTENT(in)    :: from_var
       TYPE(L2dtype), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_l2dtype
     SUBROUTINE ropp_io_assign_rotype(to_var, from_var)
       USE ropp_io_types
       TYPE(ROprof), INTENT(in)    :: from_var
       TYPE(ROprof), INTENT(inout) :: to_var
     END SUBROUTINE ropp_io_assign_rotype
  END INTERFACE

! 2.3 Checking data
! ------------------

  INTERFACE ropp_io_rangecheck
     SUBROUTINE ropp_io_rangecheck_1d(ROdata)
       USE ropp_io_types
       TYPE (ROprof), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_rangecheck_1d
     SUBROUTINE ropp_io_rangecheck_2d(ROdata)
       USE ropp_io_types
       TYPE (ROprof2d), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_rangecheck_2d
     SUBROUTINE ropp_io_rangecheck_l1atype(Lev1a)
       USE ropp_io_types
       TYPE (L1atype), INTENT(inout) :: Lev1a
     END SUBROUTINE ropp_io_rangecheck_l1atype
     SUBROUTINE ropp_io_rangecheck_l1btype(Lev1b)
       USE ropp_io_types
       TYPE (L1btype), INTENT(inout) :: Lev1b
     END SUBROUTINE ropp_io_rangecheck_l1btype
     SUBROUTINE ropp_io_rangecheck_l2atype(Lev2a, Lat)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE (L2atype), INTENT(inout) :: Lev2a
       REAL(wp), INTENT(in)          :: Lat
     END SUBROUTINE ropp_io_rangecheck_l2atype
     SUBROUTINE ropp_io_rangecheck_l2btype(Lev2b)
       USE ropp_io_types
       TYPE (L2btype), INTENT(inout) :: Lev2b
     END SUBROUTINE ropp_io_rangecheck_l2btype
     SUBROUTINE ropp_io_rangecheck_l2btype_2d(Lev2b)
       USE ropp_io_types
       TYPE (L2btype_2d), INTENT(inout) :: Lev2b
     END SUBROUTINE ropp_io_rangecheck_l2btype_2d
     SUBROUTINE ropp_io_rangecheck_l2ctype(Lev2c)
       USE ropp_io_types
       TYPE (L2ctype), INTENT(inout) :: Lev2c
     END SUBROUTINE ropp_io_rangecheck_l2ctype
     SUBROUTINE ropp_io_rangecheck_l2ctype_2d(Lev2c)
       USE ropp_io_types
       TYPE (L2ctype_2d), INTENT(inout) :: Lev2c
     END SUBROUTINE ropp_io_rangecheck_l2ctype_2d
     SUBROUTINE ropp_io_rangecheck_l2dtype(Lev2d)
       USE ropp_io_types
       TYPE (L2dtype), INTENT(inout) :: Lev2d
     END SUBROUTINE ropp_io_rangecheck_l2dtype
  END INTERFACE

! 2.4 Freeing memory
! ------------------

  INTERFACE ropp_io_free
     SUBROUTINE ropp_io_free_l1atype(var)
       USE ropp_io_types
       TYPE(L1atype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l1atype
     SUBROUTINE ropp_io_free_l1btype(var)
       USE ropp_io_types
       TYPE(L1btype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l1btype
     SUBROUTINE ropp_io_free_l2atype(var)
       USE ropp_io_types
       TYPE(L2atype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2atype
     SUBROUTINE ropp_io_free_l2btype(var)
       USE ropp_io_types
       TYPE(L2btype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2btype
     SUBROUTINE ropp_io_free_l2btype_2d(var)
       USE ropp_io_types
       TYPE(L2btype_2d), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2btype_2d
     SUBROUTINE ropp_io_free_l2ctype(var)
       USE ropp_io_types
       TYPE(L2ctype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2ctype
     SUBROUTINE ropp_io_free_l2ctype_2d(var)
       USE ropp_io_types
       TYPE(L2ctype_2d), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2ctype_2d
     SUBROUTINE ropp_io_free_l2dtype(var)
       USE ropp_io_types
       TYPE(L2dtype), INTENT(inout) :: var
     END SUBROUTINE ropp_io_free_l2dtype
     SUBROUTINE ropp_io_free_rotype(ROdata)
       USE ropp_io_types
       TYPE(ROprof), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_free_rotype
     SUBROUTINE ropp_io_free_rotype_2d(ROdata)
       USE ropp_io_types
       TYPE(ROprof2d), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_free_rotype_2d
     SUBROUTINE ropp_io_free_vlist(vlist)
       USE ropp_io_types
       TYPE(Vlisttype), INTENT(inout) :: vlist
     END SUBROUTINE ropp_io_free_vlist
     SUBROUTINE ropp_io_free_vlistD0d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD0d), POINTER :: vlist
     END SUBROUTINE ropp_io_free_vlistD0d
     SUBROUTINE ropp_io_free_vlistD1d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD1d), POINTER :: vlist
     END SUBROUTINE ropp_io_free_vlistD1d
     SUBROUTINE ropp_io_free_vlistD2d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD2d), POINTER :: vlist
     END SUBROUTINE ropp_io_free_vlistD2d
     SUBROUTINE ropp_io_free_corcov_sca(corcov)
       USE ropp_io_types
       TYPE(ROcorcov), INTENT(inout) :: corcov
     END SUBROUTINE ropp_io_free_corcov_sca
     SUBROUTINE ropp_io_free_corcov_arr(corcov)
       USE ropp_io_types
       TYPE(ROcorcov), DIMENSION(:), POINTER :: corcov
     END SUBROUTINE ropp_io_free_corcov_arr
  END INTERFACE

! 2.4 Shrinking memory
! --------------------

  INTERFACE ropp_io_shrink
     SUBROUTINE ropp_io_shrink_l1atype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L1atype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l1atype
     SUBROUTINE ropp_io_shrink_l1btype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L1btype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l1btype
     SUBROUTINE ropp_io_shrink_l2atype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L2atype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l2atype
     SUBROUTINE ropp_io_shrink_l2btype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L2btype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l2btype
     SUBROUTINE ropp_io_shrink_l2ctype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L2ctype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l2ctype
     SUBROUTINE ropp_io_shrink_l2dtype(var, imin, imax, stride)
       USE ropp_io_types
       TYPE(L2dtype), INTENT(inout) :: var
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_l2dtype
     SUBROUTINE ropp_io_shrink_rotype(ROdata, imin, imax, stride)
       USE ropp_io_types
       TYPE(ROprof), INTENT(inout) :: ROdata
       INTEGER,       INTENT(in)    :: imin
       INTEGER,       INTENT(in)    :: imax
       INTEGER,       INTENT(in)    :: stride
     END SUBROUTINE ropp_io_shrink_rotype
  END INTERFACE

! 2.5 Adding variables
! --------------------

  INTERFACE ropp_io_addvar
     SUBROUTINE ropp_io_addvar_rodataD0d(ROdata, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(ROprof),           INTENT(inout) :: ROdata
       CHARACTER(len = *),     INTENT(in)    :: name
       CHARACTER(len = *),     INTENT(in)    :: long_name
       CHARACTER(len = *),     INTENT(in)    :: units
       REAL(wp), DIMENSION(2), INTENT(in)    :: range
       REAL(wp),               INTENT(in)    :: DATA
     END SUBROUTINE ropp_io_addvar_rodataD0d
     SUBROUTINE ropp_io_addvar_rodataD1d(ROdata, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(ROprof),           INTENT(inout) :: ROdata
       CHARACTER(len = *),     INTENT(in)    :: name
       CHARACTER(len = *),     INTENT(in)    :: long_name
       CHARACTER(len = *),     INTENT(in)    :: units
       REAL(wp), DIMENSION(2), INTENT(in)    :: range
       REAL(wp), DIMENSION(:), INTENT(in)    :: DATA
     END SUBROUTINE ropp_io_addvar_rodataD1d
     SUBROUTINE ropp_io_addvar_rodataD2d(ROdata, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(ROprof),             INTENT(inout) :: ROdata
       CHARACTER(len = *),       INTENT(in)    :: name
       CHARACTER(len = *),       INTENT(in)    :: long_name
       CHARACTER(len = *),       INTENT(in)    :: units
       REAL(wp), DIMENSION(2),   INTENT(in)    :: range
       REAL(wp), DIMENSION(:,:), INTENT(in)    :: DATA
     END SUBROUTINE ropp_io_addvar_rodataD2d
     RECURSIVE SUBROUTINE ropp_io_addvar_vlistD0d(vlist, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(VlisttypeD0d),     POINTER    :: vlist
       CHARACTER(len = *),     INTENT(in) :: name
       CHARACTER(len = *),     INTENT(in) :: long_name
       CHARACTER(len = *),     INTENT(in) :: units
       REAL(wp), DIMENSION(2), INTENT(in) :: range
       REAL(wp),               INTENT(in) :: DATA
     END SUBROUTINE ropp_io_addvar_vlistD0d
     RECURSIVE SUBROUTINE ropp_io_addvar_vlistD1d(vlist, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(VlisttypeD1d),     POINTER    :: vlist
       CHARACTER(len = *),     INTENT(in) :: name
       CHARACTER(len = *),     INTENT(in) :: long_name
       CHARACTER(len = *),     INTENT(in) :: units
       REAL(wp), DIMENSION(2), INTENT(in) :: range
       REAL(wp), DIMENSION(:), INTENT(in) :: DATA
     END SUBROUTINE ropp_io_addvar_vlistD1d
     RECURSIVE SUBROUTINE ropp_io_addvar_vlistD2d(vlist, name, long_name, units, range, DATA)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       TYPE(VlisttypeD2d),       POINTER    :: vlist
       CHARACTER(len = *),       INTENT(in) :: name
       CHARACTER(len = *),       INTENT(in) :: long_name
       CHARACTER(len = *),       INTENT(in) :: units
       REAL(wp), DIMENSION(2),   INTENT(in) :: range
       REAL(wp), DIMENSION(:,:), INTENT(in) :: DATA
     END SUBROUTINE ropp_io_addvar_vlistD2d
  END INTERFACE

! 2.6 Occultation ID
! ------------------

  INTERFACE ropp_io_occid
     SUBROUTINE ropp_io_occid_1d(ROdata)
       USE ropp_io_types
       TYPE (ROprof), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_occid_1d
          SUBROUTINE ropp_io_occid_2d(ROdata)
       USE ropp_io_types
       TYPE (ROprof2d), INTENT(inout) :: ROdata
     END SUBROUTINE ropp_io_occid_2d
  END INTERFACE

! 2.7 Writing data
! ----------------

  INTERFACE ropp_io_write
     SUBROUTINE ropp_io_write_rodatass(ROdata, file, path, output_precision, append, rec, ranchk, ierr)
       USE ropp_io_types
       TYPE(ROprof),                     INTENT(inout) :: ROdata
       CHARACTER(len = *),               OPTIONAL      :: file
       CHARACTER(len = *),               OPTIONAL      :: path
       CHARACTER(len = *),               OPTIONAL      :: output_precision
       LOGICAL,                          OPTIONAL      :: append
       LOGICAL,                          OPTIONAL      :: ranchk
       INTEGER,                          OPTIONAL      :: rec
       INTEGER,                          OPTIONAL      :: ierr
     END SUBROUTINE ropp_io_write_rodatass
     SUBROUTINE ropp_io_write_rodatamm(ROdata, files, path, output_precision, multi, append, rec, ranchk, ierr)
       USE ropp_io_types
       TYPE(ROprof), DIMENSION(:),       INTENT(inout) :: ROdata
       CHARACTER(len = *), DIMENSION(:), OPTIONAL      :: files
       CHARACTER(len = *),               OPTIONAL      :: path
       CHARACTER(len = *),               OPTIONAL      :: output_precision
       LOGICAL,                          OPTIONAL      :: multi
       LOGICAL,                          OPTIONAL      :: append
       LOGICAL,                          OPTIONAL      :: ranchk
       INTEGER,                          OPTIONAL      :: rec
       INTEGER,                          OPTIONAL      :: ierr
     END SUBROUTINE ropp_io_write_rodatamm
     SUBROUTINE ropp_io_write_rodata2d(ROdata, file, path, output_precision, append, rec, ranchk, ierr)
       USE ropp_io_types
       TYPE(ROprof2d),                   INTENT(inout) :: ROdata
       CHARACTER(len = *),               OPTIONAL      :: file
       CHARACTER(len = *),               OPTIONAL      :: path
       CHARACTER(len = *),               OPTIONAL      :: output_precision
       LOGICAL,                          OPTIONAL      :: append
       LOGICAL,                          OPTIONAL      :: ranchk
       INTEGER,                          OPTIONAL      :: rec
       INTEGER,                          OPTIONAL      :: ierr
     END SUBROUTINE ropp_io_write_rodata2d
  END INTERFACE

  INTERFACE ropp_io_write_ncdf_def
     SUBROUTINE ropp_io_write_def_rodata(DATA, output_precision)
       USE ropp_io_types
       TYPE(ROprof)                 :: DATA
       CHARACTER(len = *), OPTIONAL :: output_precision
     END SUBROUTINE ropp_io_write_def_rodata
     SUBROUTINE ropp_io_write_def_rodata_2d(DATA, output_precision)
       USE ropp_io_types
       TYPE(ROprof2d)               :: DATA
       CHARACTER(len = *), OPTIONAL :: output_precision
     END SUBROUTINE ropp_io_write_def_rodata_2d
     RECURSIVE SUBROUTINE ropp_io_write_def_vlistD0d(vlist, output_precision)
       USE ropp_io_types
       TYPE(VlisttypeD0d)           :: vlist
       CHARACTER(len = *), OPTIONAL :: output_precision
     END SUBROUTINE ropp_io_write_def_vlistD0d
     RECURSIVE SUBROUTINE ropp_io_write_def_vlistD1d(vlist, output_precision)
       USE ropp_io_types
       TYPE(VlisttypeD1d)           :: vlist
       CHARACTER(len = *), OPTIONAL :: output_precision
     END SUBROUTINE ropp_io_write_def_vlistD1d
     RECURSIVE SUBROUTINE ropp_io_write_def_vlistD2d(vlist, output_precision)
       USE ropp_io_types
       TYPE(VlisttypeD2d)           :: vlist
       CHARACTER(len = *), OPTIONAL :: output_precision
     END SUBROUTINE ropp_io_write_def_vlistD2d
  END INTERFACE

  INTERFACE ropp_io_write_ncdf_put
     SUBROUTINE ropp_io_write_put_rodata(DATA, rec)
       USE ropp_io_types
       TYPE(ROprof), INTENT(in) :: DATA
       INTEGER,      OPTIONAL   :: rec
     END SUBROUTINE ropp_io_write_put_rodata
     SUBROUTINE ropp_io_write_put_rodata_2d(DATA, rec)
       USE ropp_io_types
       TYPE(ROprof2d), INTENT(in) :: DATA
       INTEGER,        OPTIONAL   :: rec
     END SUBROUTINE ropp_io_write_put_rodata_2d
     RECURSIVE SUBROUTINE ropp_io_write_put_vlistD0d(vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD0d), INTENT(in) :: vlist
       INTEGER,            OPTIONAL   :: rec
     END SUBROUTINE ropp_io_write_put_vlistD0d
     RECURSIVE SUBROUTINE ropp_io_write_put_vlistD1d(vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD1d), INTENT(in) :: vlist
       INTEGER,            OPTIONAL   :: rec
     END SUBROUTINE ropp_io_write_put_vlistD1d
     RECURSIVE SUBROUTINE ropp_io_write_put_vlistD2d(vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD2d), INTENT(in) :: vlist
       INTEGER,            OPTIONAL   :: rec
     END SUBROUTINE ropp_io_write_put_vlistD2d
  END INTERFACE

! 2.8 Reading data
! ----------------

  INTERFACE ropp_io_read
     SUBROUTINE ropp_io_read_rodatas(ROdata, file, rec, centre, ierr, ranchk, resolution, getlevel1a, getbufr)
       USE ropp_io_types
       TYPE(ROprof),         INTENT(inout) :: ROdata
       CHARACTER(len = *),   INTENT(in)    :: file
       INTEGER,              OPTIONAL      :: rec
       CHARACTER(len=20),    OPTIONAL      :: centre
       INTEGER,              OPTIONAL      :: ierr
       LOGICAL,              OPTIONAL      :: ranchk
       CHARACTER(len = 20),  OPTIONAL      :: resolution
       LOGICAL,              OPTIONAL      :: getlevel1a
       LOGICAL,              OPTIONAL      :: getbufr
     END SUBROUTINE ropp_io_read_rodatas
     SUBROUTINE ropp_io_read_rodatam(ROdata, file, ierr, ranchk)
       USE ropp_io_types
       TYPE(ROprof), DIMENSION(:),   POINTER     :: ROdata
       CHARACTER(len = *),           INTENT(in)  :: file
       INTEGER,            OPTIONAL, INTENT(out) :: ierr
       LOGICAL,            OPTIONAL, INTENT(in)  :: ranchk
     END SUBROUTINE ropp_io_read_rodatam
     SUBROUTINE ropp_io_read_rodata_2d(ROdata, file, rec, ierr, ranchk)
       USE ropp_io_types
       TYPE(ROprof2d),     INTENT(inout) :: ROdata
       CHARACTER(len = *), INTENT(in)    :: file
       INTEGER,            OPTIONAL      :: rec
       INTEGER,            OPTIONAL      :: ierr
       LOGICAL,            OPTIONAL      :: ranchk
     END SUBROUTINE ropp_io_read_rodata_2d
     SUBROUTINE ropp_io_read_rocorcov(ROdata, file, ierr)
       USE ropp_io_types
       TYPE(ROcorcov), DIMENSION(:), POINTER     :: ROdata
       CHARACTER(len = *),           INTENT(in)  :: file
       INTEGER,            OPTIONAL, INTENT(out) :: ierr
     END SUBROUTINE ropp_io_read_rocorcov
  END INTERFACE

  INTERFACE ropp_io_read_ncdf_get
     SUBROUTINE ropp_io_read_ncdf_get_rodata(DATA, rec)
       USE ropp_io_types
       TYPE(ROprof),       INTENT(inout) :: DATA
       INTEGER,            OPTIONAL      :: rec
     END SUBROUTINE ropp_io_read_ncdf_get_rodata
     SUBROUTINE ropp_io_read_ncdf_get_rodata_2d(DATA, rec)
       USE ropp_io_types
       TYPE(ROprof2d),     INTENT(inout) :: DATA
       INTEGER,            OPTIONAL      :: rec
     END SUBROUTINE ropp_io_read_ncdf_get_rodata_2d
     SUBROUTINE ropp_io_read_ncdf_get_rocorcov(DATA)
       USE ropp_io_types
       TYPE(ROcorcov), DIMENSION(:), POINTER :: DATA
     END SUBROUTINE ropp_io_read_ncdf_get_rocorcov
     SUBROUTINE ropp_io_read_ncdf_get_ucardata(DATA, file)
       USE ropp_io_types
       TYPE(ROprof),       INTENT(inout) :: DATA
       CHARACTER(len=*),   INTENT(in)    :: file
     END SUBROUTINE ropp_io_read_ncdf_get_ucardata
     SUBROUTINE ropp_io_read_ncdf_get_otherdata(DATA, file, centre, rec, resolution, getlevel1a, getbufr)
       USE ropp_io_types
       TYPE(ROprof),       INTENT(inout) :: DATA
       CHARACTER(len=*),   INTENT(in)    :: file
       CHARACTER(len=20),  INTENT(in)    :: centre
       LOGICAL,            OPTIONAL      :: getlevel1a
       LOGICAL,            OPTIONAL      :: getbufr
       CHARACTER(len=20),  OPTIONAL      :: resolution
       INTEGER,            OPTIONAL      :: rec
    END SUBROUTINE ropp_io_read_ncdf_get_otherdata
     SUBROUTINE ropp_io_read_ncdf_get_eumdata(DATA, file, resolution, getlevel1a, getbufr, ldummy)
       USE ropp_io_types
       TYPE(ROprof),       INTENT(inout) :: DATA
       CHARACTER(len=*),   INTENT(in)    :: file
       CHARACTER(len=20),  INTENT(in)    :: resolution
       LOGICAL,            INTENT(IN)    :: getlevel1a
       LOGICAL,            INTENT(IN)    :: getbufr
       LOGICAL,            INTENT(IN)    :: ldummy
     END SUBROUTINE ropp_io_read_ncdf_get_eumdata
  END INTERFACE

  INTERFACE ropp_io_vlist_read
     RECURSIVE SUBROUTINE ropp_io_read_ncdf_get_vlistD0d(varid, vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD0d), POINTER    :: vlist
       INTEGER,            INTENT(in) :: varid
       INTEGER,            INTENT(in) :: rec
     END SUBROUTINE ropp_io_read_ncdf_get_vlistD0d
     RECURSIVE SUBROUTINE ropp_io_read_ncdf_get_vlistD1d(varid, vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD1d), POINTER    :: vlist
       INTEGER,            INTENT(in) :: varid
       INTEGER,            INTENT(in) :: rec
     END SUBROUTINE ropp_io_read_ncdf_get_vlistD1d
     RECURSIVE SUBROUTINE ropp_io_read_ncdf_get_vlistD2d(varid, vlist, rec)
       USE ropp_io_types
       TYPE(VlisttypeD2d), POINTER    :: vlist
       INTEGER,            INTENT(in) :: varid
       INTEGER,            INTENT(in) :: rec
     END SUBROUTINE ropp_io_read_ncdf_get_vlistD2d
  END INTERFACE

! 2.9 Vlists
! ----------

  INTERFACE size
     FUNCTION size_vlist(vlist) RESULT(n)
       USE ropp_io_types
       TYPE(Vlisttype) :: vlist
       INTEGER         :: n
     END FUNCTION size_vlist
     RECURSIVE FUNCTION size_vlistD0d(vlist) RESULT(n)
       USE ropp_io_types
       TYPE(VlisttypeD0d), POINTER :: vlist
       INTEGER                     :: n
     END FUNCTION size_vlistD0d
     RECURSIVE FUNCTION size_vlistD1d(vlist) RESULT(n)
       USE ropp_io_types
       TYPE(VlisttypeD1d), POINTER :: vlist
       INTEGER                     :: n
     END FUNCTION size_vlistD1d
     RECURSIVE FUNCTION size_vlistD2d(vlist) RESULT(n)
       USE ropp_io_types
       TYPE(VlisttypeD2d), POINTER :: vlist
       INTEGER                     :: n
     END FUNCTION size_vlistD2d
  END INTERFACE

  INTERFACE PRINT
     SUBROUTINE print_vlist(vlist)
       USE ropp_io_types
       TYPE(Vlisttype)       :: vlist
     END SUBROUTINE print_vlist
     RECURSIVE SUBROUTINE print_vlistD0d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD0d), POINTER :: vlist
     END SUBROUTINE print_vlistD0d
     RECURSIVE SUBROUTINE print_vlistD1d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD1d), POINTER :: vlist
     END SUBROUTINE print_vlistD1d
     RECURSIVE SUBROUTINE print_vlistD2d(vlist)
       USE ropp_io_types
       TYPE(VlisttypeD2d), POINTER :: vlist
     END SUBROUTINE print_vlistD2d
  END INTERFACE

! 2.10 Thinning routines
! ---------------------

  INTERFACE
    SUBROUTINE ropp_io_thin ( rodata,   &
                              thinfile, &
                              impactalt, &
                              ranchk )
      USE ropp_io_types, ONLY: roprof
      TYPE(roprof)                            :: rodata
      CHARACTER (len=*),           INTENT(in) :: thinfile
      LOGICAL,           OPTIONAL, INTENT(in) :: impactalt
      LOGICAL,           OPTIONAL, INTENT(in) :: ranchk
    END SUBROUTINE ropp_io_thin
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_io_thin_select ( nlev,     lev,     val,     &
                                     nthinlev, thinlev, thinval, &
                                     method,   nsamp,   sigma )
      USE typesizes, ONLY: wp => eightbytereal
      CHARACTER (len=*), INTENT(in)    :: method
      INTEGER,           INTENT(in)    :: nlev
      INTEGER,           INTENT(in)    :: nthinlev
      REAL(wp),          INTENT(inout) :: lev(1:nlev)
      REAL(wp),          INTENT(inout) :: val(1:nlev)
      REAL(wp),          INTENT(inout) :: thinlev(1:nthinlev)
      REAL(wp),          INTENT(out)   :: thinval(1:nthinlev)
      INTEGER,           INTENT(out)   :: nsamp
      LOGICAL, OPTIONAL, INTENT(in)    :: sigma
    END SUBROUTINE ropp_io_thin_select
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_io_thin_fixed ( nlev,     lev,     val,     &
                                    nthinlev, thinlev, thinval, &
                                    method,   sigma)
      USE typesizes, ONLY: wp => eightbytereal
      CHARACTER (len=*), INTENT(in)  :: method
      INTEGER,           INTENT(in)  :: nlev
      INTEGER,           INTENT(in)  :: nthinlev
      REAL(wp),          INTENT(in)  :: lev(1:nlev)
      REAL(wp),          INTENT(in)  :: val(1:nlev)
      REAL(wp),          INTENT(in)  :: thinlev(1:nthinlev)
      REAL(wp),          INTENT(out) :: thinval(1:nthinlev)
      LOGICAL, OPTIONAL, INTENT(in)  :: sigma
    END SUBROUTINE ropp_io_thin_fixed
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_io_thin_skip ( maxpts, npoints,        &
                                   skip1,  skip2,   nsamp )
      INTEGER, INTENT(in)           :: maxpts
      INTEGER, INTENT(in)           :: npoints
      INTEGER, INTENT(out)          :: skip1, skip2
      INTEGER, INTENT(out)          :: nsamp
    END SUBROUTINE ropp_io_thin_skip
  END INTERFACE

  INTERFACE ropp_io_thin_sg
    SUBROUTINE ropp_io_thin_sg (nlev, val, npoints, order)
      USE typeSizes, ONLY: wp => EightByteReal
      INTEGER,           INTENT(in)    :: nlev
      REAL(wp),          INTENT(inout) :: val(1:nlev)
      INTEGER, OPTIONAL, INTENT(in)    :: npoints
      INTEGER, OPTIONAL, INTENT(in)    :: order
    END SUBROUTINE ropp_io_thin_sg
  END INTERFACE


! 2.11 Utility routines
! --------------------

  INTERFACE
     FUNCTION ropp_io_nrec(file) RESULT(n)
       CHARACTER(len = *)   :: file
       INTEGER              :: n
     END FUNCTION ropp_io_nrec
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_io_ascend ( ROdata )
      USE ropp_io_types, ONLY: roprof
      TYPE(roprof)        :: ROdata
    END SUBROUTINE ropp_io_ascend
  END INTERFACE

  INTERFACE isroppinrange
     FUNCTION ropp_io_isinrangeDT7(dt) RESULT(inrange)
       USE ropp_io_types
       TYPE(DT7type), INTENT(in) :: dt
       LOGICAL                   :: inrange
     END FUNCTION ropp_io_isinrangeDT7
  END INTERFACE

  INTERFACE ropp_io_version
    FUNCTION ropp_io_version() RESULT (version)
      CHARACTER (LEN=40) :: version
    END FUNCTION ropp_io_version
  END INTERFACE ropp_io_version

END MODULE ropp_io
