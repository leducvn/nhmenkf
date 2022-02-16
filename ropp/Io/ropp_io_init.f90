! $Id: ropp_io_init.f90 4060 2014-02-03 11:13:00Z idculv $

!****s* Initialisation/ropp_io_init *
!
! NAME
!    ropp_io_init - Initialise an RO derived type with dummy 'missing' data
!
! SYNOPSIS
!    use ropp_io
!    type(ROprof) :: ro_data
!      ...
!    call ropp_io_init(ro_data, NLev1a, NLev1b, NLev2a, NLev2b, Nlev2c, NLev2d)
!    call ropp_io_init(ro_data%Lev1a, NLev1a)
!    call ropp_io_init(ro_data%Lev1b, NLev1b)
!    call ropp_io_init(ro_data%Lev2a, NLev2a)
!    call ropp_io_init(ro_data%Lev2b, NLev2b)
!    call ropp_io_init(ro_data%Lev2c, NLev2c)
!    call ropp_io_init(ro_data%Lev2d, NLev2d)
!
! DESCRIPTION
!   This subroutine initialises a RO data structure or parts thereof with dummy
!   or 'missing' data values for all elements, appropriate to the parameter and
!   suitable for the ROPP data format.
!
! OUTPUT
!    ro_data   dtyp  RO data (derived type)
!
! SEE ALSO
!    ropp_io_types
!
! NOTES
!    The subroutine call ropp_io_init() for level2c data (i.e., surface parameters)
!    does not do much - it simply sets the number of elements to 1, regardless of
!    the argument. The subroutine is available for the sake of user interface
!    consistency only.
!
! REFERENCES
!    Format Definition for Radio Occultation Files -
!    CLIMAP Format Version 2.2a, Issue 1.6, 8 January 2004
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
! 1. L1atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l1atype(var, n)

! 1.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l1atype
  USE ropp_io_types, ONLY: L1atype

  IMPLICIT NONE

  TYPE(L1atype), INTENT(inout) :: var
  INTEGER                      :: n

! 1.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n > 0) THEN
    CALL callocate(var%dtime,      n,   ropp_MDFV)
    CALL callocate(var%snr_L1ca,   n,   ropp_MDFV)
    CALL callocate(var%snr_L1p,    n,   ropp_MDFV)
    CALL callocate(var%snr_L2p,    n,   ropp_MDFV)
    CALL callocate(var%phase_L1,   n,   ropp_MDFV)
    CALL callocate(var%phase_L2,   n,   ropp_MDFV)
    CALL callocate(var%phase_qual, n,   ropp_MDFV)
    CALL callocate(var%r_gns, (/n, 3/), ropp_MDFV)
    CALL callocate(var%v_gns, (/n, 3/), ropp_MDFV)
    CALL callocate(var%r_leo, (/n, 3/), ropp_MDFV)
    CALL callocate(var%v_leo, (/n, 3/), ropp_MDFV)
    var%Npoints = n
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l1atype


!-------------------------------------------------------------------------------
! 2. L1btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l1btype(var, n)

! 2.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l1btype
  USE ropp_io_types, ONLY: L1btype

  IMPLICIT NONE

  TYPE(L1btype), INTENT(inout) :: var
  INTEGER                      :: n

! 2.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n > 0) THEN
    CALL callocate(var%lat_tp,           n, ropp_MDFV)
    CALL callocate(var%lon_tp,           n, ropp_MDFV)
    CALL callocate(var%azimuth_tp,       n, ropp_MDFV)
    CALL callocate(var%impact_L1,        n, ropp_MDFV)
    CALL callocate(var%impact_L2,        n, ropp_MDFV)
    CALL callocate(var%impact,           n, ropp_MDFV)
    CALL callocate(var%impact_Opt,       n, ropp_MDFV)
    CALL callocate(var%bangle_L1,        n, ropp_MDFV)
    CALL callocate(var%bangle_L2,        n, ropp_MDFV)
    CALL callocate(var%bangle,           n, ropp_MDFV)
    CALL callocate(var%bangle_Opt,       n, ropp_MDFV)
    CALL callocate(var%bangle_L1_sigma,  n, ropp_MDFV)
    CALL callocate(var%bangle_L2_sigma,  n, ropp_MDFV)
    CALL callocate(var%bangle_sigma,     n, ropp_MDFV)
    CALL callocate(var%bangle_Opt_sigma, n, ropp_MDFV)
    CALL callocate(var%bangle_L1_qual,   n, ropp_MDFV)
    CALL callocate(var%bangle_L2_qual,   n, ropp_MDFV)
    CALL callocate(var%bangle_qual,      n, ropp_MDFV)
    CALL callocate(var%bangle_Opt_qual,  n, ropp_MDFV)
    var%Npoints = n
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l1btype


!-------------------------------------------------------------------------------
! 3. L2atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2atype(var, n)

! 3.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l2atype
  USE ropp_io_types, ONLY: L2atype

  IMPLICIT NONE

  TYPE(L2atype), INTENT(inout) :: var
  INTEGER                      :: n

! 3.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n > 0) THEN
    CALL callocate(var%alt_refrac,     n, ropp_MDFV)
    CALL callocate(var%geop_refrac,    n, ropp_MDFV)
    CALL callocate(var%refrac,         n, ropp_MDFV)
    CALL callocate(var%refrac_sigma,   n, ropp_MDFV)
    CALL callocate(var%refrac_qual,    n, ropp_MDFV)
    CALL callocate(var%dry_temp,       n, ropp_MDFV)
    CALL callocate(var%dry_temp_sigma, n, ropp_MDFV)
    CALL callocate(var%dry_temp_qual,  n, ropp_MDFV)
    var%Npoints = n
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l2atype


!-------------------------------------------------------------------------------
! 4a. L2btype (1d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2btype(var, n)

! 4.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l2btype
  USE ropp_io_types, ONLY: L2btype

  IMPLICIT NONE

  TYPE(L2btype), INTENT(inout) :: var
  INTEGER                      :: n

! 4.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n > 0) THEN
    CALL callocate(var%geop,        n, ropp_MDFV)
    CALL callocate(var%geop_sigma,  n, ropp_MDFV)
    CALL callocate(var%press,       n, ropp_MDFV)
    CALL callocate(var%press_sigma, n, ropp_MDFV)
    CALL callocate(var%temp,        n, ropp_MDFV)
    CALL callocate(var%temp_sigma,  n, ropp_MDFV)
    CALL callocate(var%shum,        n, ropp_MDFV)
    CALL callocate(var%shum_sigma,  n, ropp_MDFV)
    CALL callocate(var%meteo_qual,  n, ropp_MDFV)
    var%Npoints = n
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l2btype

!-------------------------------------------------------------------------------
! 4b. L2btype (2d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2btype_2d(var, n)

! 4b.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l2btype_2d
  USE ropp_io_types, ONLY: ropp_MDFV,  &
                           L2btype_2d

  IMPLICIT NONE

  TYPE(L2btype_2d), INTENT(inout) :: var
  INTEGER, DIMENSION(2)           :: n

! 4b.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n(1) > 0) THEN
    CALL callocate(var%geop,        n, ropp_MDFV)
    CALL callocate(var%geop_sigma,  n, ropp_MDFV)
    CALL callocate(var%press,       n, ropp_MDFV)
    CALL callocate(var%press_sigma, n, ropp_MDFV)
    CALL callocate(var%temp,        n, ropp_MDFV)
    CALL callocate(var%temp_sigma,  n, ropp_MDFV)
    CALL callocate(var%shum,        n, ropp_MDFV)
    CALL callocate(var%shum_sigma,  n, ropp_MDFV)
    CALL callocate(var%meteo_qual,  n, ropp_MDFV)
    var%Npoints = n(1)
    var%Nhoriz  = n(2)
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l2btype_2d

!-------------------------------------------------------------------------------
! 5a. L2ctype (1d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2ctype(var, n)

! 5.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io, not_this => ropp_io_init_l2ctype
  USE ropp_io_types, ONLY: L2ctype

  IMPLICIT NONE

  TYPE(L2ctype), INTENT(inout) :: var
  INTEGER                      :: n

! 5.2 Allocate memory for all structure elements
! ----------------------------------------------

  CALL ropp_io_free(var)

  IF (n > 0) THEN
    var%Npoints = 1
  ENDIF

END SUBROUTINE ropp_io_init_l2ctype

!-------------------------------------------------------------------------------
! 5b. L2ctype (2d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2ctype_2d(var, n)

! 5b.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io, not_this => ropp_io_init_l2ctype_2d
  USE ropp_io_types, ONLY: L2ctype_2d

  IMPLICIT NONE

  TYPE(L2ctype_2d), INTENT(inout) :: var
  INTEGER, DIMENSION(2)           :: n

! 5b.2 Allocate memory for all structure elements
! ----------------------------------------------

 IF (n(1) > 0) THEN
    CALL callocate(var%lat_2d,          n(2), ropp_MDFV)
    CALL callocate(var%lon_2d,          n(2), ropp_MDFV)
    CALL callocate(var%geop_sfc,        n(2), ropp_MDFV)
    CALL callocate(var%press_sfc,       n(2), ropp_MDFV)
    CALL callocate(var%press_sfc_sigma, n(2), ropp_MDFV)
    CALL callocate(var%press_sfc_qual,  n(2), ropp_MDFV)
    var%Npoints = 1
    var%Nhoriz  = n(2)
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l2ctype_2d

!-------------------------------------------------------------------------------
! 6. L2dtype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_l2dtype(var, n)

! 6.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_l2dtype
  USE ropp_io_types, ONLY: L2dtype

  IMPLICIT NONE

  TYPE(L2dtype), INTENT(inout) :: var
  INTEGER                      :: n

! 6.2 Allocate memory for all structure elements
! ----------------------------------------------

  IF (n > 0) THEN
    var%level_type = "UNKNOWN"
    CALL callocate(var%level_coeff_a, n, ropp_MDFV)
    CALL callocate(var%level_coeff_b, n, ropp_MDFV)
    var%Npoints = n
  ELSE
    CALL ropp_io_free(var)
  ENDIF

END SUBROUTINE ropp_io_init_l2dtype


!-------------------------------------------------------------------------------
! 7. Vlist
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_vlist(var)

! 6.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_init_vlist
  USE ropp_io_types, ONLY: Vlisttype

  IMPLICIT NONE

  TYPE(Vlisttype), INTENT(inout) :: var

! 6.2 Clear all previously existing additional variables
! -----------------------------------------------------

  CALL ropp_io_free(var)

END SUBROUTINE ropp_io_init_vlist


!-------------------------------------------------------------------------------
! 7. Joint RO data type
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_rotype(ROdata,         &
                               NLev1a, NLev1b, &
                               NLev2a, NLev2b, Nlev2c, NLev2d)

! 7.1 Declarations
! ----------------

  USE DateTimeProgs,  ONLY: Date_and_Time_UTC
  USE DateTimeTypes
! USE ropp_io,        not_this => ropp_io_init_rotype
  USE ropp_io_types,  ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof), INTENT(inout) :: ROdata
  INTEGER,      INTENT(in)    :: NLev1a
  INTEGER,      INTENT(in)    :: NLev1b
  INTEGER,      INTENT(in)    :: NLev2a
  INTEGER,      INTENT(in)    :: NLev2b
  INTEGER,      INTENT(in)    :: NLev2c
  INTEGER,      INTENT(in)    :: NLev2d

  INTEGER, DIMENSION(8)       :: DTnow

! 7.2 Time of processing is current date/time (UTC)
! -------------------------------------------------

  CALL Date_and_Time_UTC ( Values=DTnow )
  ROdata%DTpro%Year   = DTnow(IdxYear)
  ROdata%DTpro%Month  = DTnow(IdxMonth)
  ROdata%DTpro%Day    = DTnow(IdxDay)
  ROdata%DTpro%Hour   = DTnow(IdxHour)
  ROdata%DTpro%Minute = DTnow(IdxMinute)
  ROdata%DTpro%Second = DTnow(IdxSecond)
  ROdata%DTpro%Msec   = DTnow(IdxMSec)

! 7.3 Set PCD word to all 1's, including 'missing' flag bit
! ---------------------------------------------------------

  ROdata%PCD = 65535

! 7.4 Level 1a profile
! --------------------

  CALL ropp_io_init(ROdata%Lev1a, NLev1a)

! 7.5 Level 1b profile
! --------------------

  CALL ropp_io_init(ROdata%Lev1b, NLev1b)

! 7.6 Level 2a profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2a, NLev2a)

! 7.7 Level 2b profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2b, NLev2b)

! 7.8 Level 2c profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2c, NLev2c)

! 7.9 Level 2d profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2d, NLev2d)

! 7.10 Additional variables
! ------------------------

  CALL ropp_io_init(ROdata%vlist)

END SUBROUTINE ropp_io_init_rotype

!-------------------------------------------------------------------------------
! 8. Joint RO data type (two-dimensional meteorological data)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_init_rotype_2d(ROdata,         &
                                  NLev1a, NLev1b, &
                                  NLev2a, NLev2b, Nlev2c, NLev2d, NHoriz)

! 8.1 Declarations
! ----------------

  USE DateTimeProgs,  ONLY: Date_and_Time_UTC
  USE DateTimeTypes
! USE ropp_io,        not_this => ropp_io_init_rotype_2d
  USE ropp_io_types,  ONLY: ROprof2d

  IMPLICIT NONE

  TYPE(ROprof2d), INTENT(inout) :: ROdata
  INTEGER,        INTENT(in)    :: NLev1a
  INTEGER,        INTENT(in)    :: NLev1b
  INTEGER,        INTENT(in)    :: NLev2a
  INTEGER,        INTENT(in)    :: NLev2b
  INTEGER,        INTENT(in)    :: NLev2c
  INTEGER,        INTENT(in)    :: NLev2d
  INTEGER,        INTENT(in)    :: NHoriz

  INTEGER, DIMENSION(2)         :: ndim
  INTEGER, DIMENSION(8)         :: DTnow

! 8.2 Time of processing is current date/time (UTC)
! -------------------------------------------------

  CALL Date_and_Time_UTC ( Values=DTnow )
  ROdata%DTpro%Year   = DTnow(IdxYear)
  ROdata%DTpro%Month  = DTnow(IdxMonth)
  ROdata%DTpro%Day    = DTnow(IdxDay)
  ROdata%DTpro%Hour   = DTnow(IdxHour)
  ROdata%DTpro%Minute = DTnow(IdxMinute)
  ROdata%DTpro%Second = DTnow(IdxSecond)
  ROdata%DTpro%Msec   = DTnow(IdxMSec)

! 8.3 Set PCD word to all 1's, including 'missing' flag bit
! ---------------------------------------------------------

  ROdata%PCD = 65535

! 8.4 Level 1a profile
! --------------------

  CALL ropp_io_init(ROdata%Lev1a, NLev1a)

! 8.5 Level 1b profile
! --------------------

  CALL ropp_io_init(ROdata%Lev1b, NLev1b)

! 8.6 Level 2a profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2a, NLev2a)

! 8.7 Level 2b profile
! --------------------

  ndim(1) = NLev2b
  ndim(2) = NHoriz
  CALL ropp_io_init(ROdata%Lev2b, ndim)

! 8.8 Level 2c profile
! --------------------

  ndim(1) = NLev2c
  ndim(2) = NHoriz
  CALL ropp_io_init(ROdata%Lev2c, ndim)

! 8.9 Level 2d profile
! --------------------

  CALL ropp_io_init(ROdata%Lev2d, NLev2d)

! 8.10 Additional variables
! ------------------------

  CALL ropp_io_init(ROdata%vlist)

END SUBROUTINE ropp_io_init_rotype_2d

