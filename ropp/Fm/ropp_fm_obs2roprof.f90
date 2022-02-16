! $Id: ropp_fm_obs2roprof.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Copying2/ropp_fm_obs2roprof *
!
! NAME
!    ropp_fm_obs2roprof - Copy elements of an observation vector to an ROprof
!                            structure
!
! SYNOPSIS
!    type(<some obs vector type>) :: obs
!    type(ROprof)                 :: ro_data
!       ...
!    call ropp_fm_obs2roprof(obs, ro_data)
! 
! DESCRIPTION
!    This subroutine copies Level 1b (bending angle) or level 2a (refractivity)
!    data as contained in an observation vector to a radio occultation profile
!    data structure.
!
! INPUTS
!   obs      An observation state vector (Obs1dBangle, Obs1dRefrac)
!
! OUTPUT
!   ro_data  Radio occultation profile data.
!
! NOTES
!   Data is copied into the ROprof data structure without unit conversion;
!   thus, the units in the ro_data data structure for Level 1b and 2a must
!   be set to the units used by the observation vector variables. For the
!   units internally used within ropp_fm, this can be accomplished with
!   ropp_fm_set_units().
!
!   Currently, only 1-dimensional observation vectors for bending angle and
!   refractivity data (i.e., those of type(Obs1dbangle) or type(Obs1dRefrac))
!   are supported. In this case, the longitude and latitude coordinates of
!   the tangential points are taken from the georeferencing coordinates
!   contained in the header of the ROprof structure, assuming that they
!   reflect the profile's location properly.
!
!   The 1d bending angle observation vector contains a single vertical 
!   profile of bending angles only; it is assumed that these represent
!   neutral atmospheric bending. Thus, the data is copied into the generic
!   impact/bangle components of the ro_data structure. The earth and 
!   curvature radius information contained in the bending angle observation
!   structure is not further exploited; it is assumed that ro_data contains
!   the proper geolocation data in its header.
!
!   For refractivity, the calcultion of geopotential to altitude above the
!   reference ellipsoid is based on the latitude coordinate contained in the
!   header.
!
! SEE ALSO
!   Obs1dBangle
!   Obs1dRefrac
!   ropp_fm_roprof2obs
!   ropp_fm_set_units
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
! 1. 1D Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_obs1dbangle2roprof(y, ro_data)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dBangle

  IMPLICIT NONE

  TYPE(Obs1dBangle),    INTENT(in)    :: y          ! Bending angle structure
  TYPE(ROprof),         INTENT(inout) :: ro_data    ! RO data structure

  INTEGER                             :: i, n
  CHARACTER(len = 256)                :: routine

! 1.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_obs2roprof (1D bending angles)')

! 1.3 Copy data
! -------------

  CALL ropp_io_init(ro_data%Lev1b, int(SIZE(y%bangle)))

  ro_data%GEOref%lon   = y%lon
  ro_data%GEOref%lat   = y%lat
  ro_data%GEOref%azimuth      = y%azimuth
  ro_data%GEOref%undulation   = y%undulation
  ro_data%GEOref%roc          = y%r_curve

  ro_data%Lev1b%lon_tp = ro_data%GEOref%lon
  ro_data%Lev1b%lat_tp = ro_data%GEOref%lat

  ro_data%Lev1b%impact = y%impact
  ro_data%Lev1b%bangle = y%bangle

  IF(y%obs_ok)THEN
     ro_data%Lev1b%bangle_qual(:) = 100.0_wp
  ELSE
     ro_data%Lev1b%bangle_qual(:) = 0.0_wp
  ENDIF
  
! 1.4 Copy sigmas
! ---------------

  n = SIZE(y%bangle)

  IF (ASSOCIATED(y%cov%d) .AND. SIZE(y%cov%d) == n*(n+1)/2) THEN
     DO i = 1, SIZE(y%bangle)
       ! matrix_pp type, uplo = 'U'
        ro_data%Lev1b%bangle_sigma(i) = SQRT(y%cov%d(i + i*(i-1)/2)) 
     ENDDO
  ENDIF

! 1.5 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_obs1dbangle2roprof


!-------------------------------------------------------------------------------
! 2. 1D Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_obs1drefrac2roprof(y, ro_data)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dRefrac

  IMPLICIT NONE

  TYPE(Obs1dRefrac),    INTENT(in)    :: y        ! Refractivity structure
  TYPE(ROprof),         INTENT(inout) :: ro_data  ! RO data structure

  INTEGER                             :: i, n
  CHARACTER(len = 256)                :: routine

! 2.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_obs2roprof (1D refractivity)')

! 2.3 Copy data
! -------------

  CALL ropp_io_init(ro_data%Lev2a, int(SIZE(y%refrac)))

  ro_data%GEOref%lon   = y%lon
  ro_data%GEOref%lat   = y%lat

  ro_data%Lev2a%alt_refrac  = geopotential2geometric(ro_data%GEOref%lat, &
                                                     y%geop)
  ro_data%Lev2a%geop_refrac = y%geop
  ro_data%Lev2a%refrac      = y%refrac

  IF(y%obs_ok)THEN
     ro_data%Lev2a%refrac_qual(:) = 100.0_wp
  ELSE
     ro_data%Lev2a%refrac_qual(:) = 0.0_wp
  ENDIF

! 2.4 Copy sigmas
! ---------------

  n = SIZE(y%refrac)

  IF (ASSOCIATED(y%cov%d) .AND. SIZE(y%cov%d) == n*(n+1)/2) THEN
     DO i = 1, n
       ! matrix_pp type, uplo = 'U'
        ro_data%Lev2a%refrac_sigma(i) = SQRT(y%cov%d(i + i*(i-1)/2)) 
     ENDDO
  ENDIF

! 2.5 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_obs1drefrac2roprof

!-------------------------------------------------------------------------------
! 3. 1D Bending angles -> 2D roprof
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_obs2dbangle2roprof(y, ro_data)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof2d
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dBangle

  IMPLICIT NONE

  TYPE(Obs1dBangle),    INTENT(in)    :: y          ! Bending angle structure
  TYPE(ROprof2d),       INTENT(inout) :: ro_data    ! RO data structure

  INTEGER                             :: i, n
  CHARACTER(len = 256)                :: routine

! 1.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_obs2roprof (1D bending angles)')

! 1.3 Copy data
! -------------

  CALL ropp_io_init(ro_data%Lev1b, int(SIZE(y%bangle)))

  ro_data%GEOref%lon   = y%lon
  ro_data%GEOref%lat   = y%lat
  ro_data%GEOref%azimuth      = y%azimuth
  ro_data%GEOref%undulation   = y%undulation
  ro_data%GEOref%roc          = y%r_curve

  ro_data%Lev1b%lon_tp = ro_data%GEOref%lon
  ro_data%Lev1b%lat_tp = ro_data%GEOref%lat

  ro_data%Lev1b%impact = y%impact
  ro_data%Lev1b%bangle = y%bangle

  IF(y%obs_ok)THEN
     ro_data%Lev1b%bangle_qual(:) = 100.0_wp
  ELSE
     ro_data%Lev1b%bangle_qual(:) = 0.0_wp
  ENDIF


! 1.4 Copy sigmas
! ---------------

  n = SIZE(y%bangle)

  IF (ASSOCIATED(y%cov%d) .AND. SIZE(y%cov%d) == n*(n+1)/2) THEN
     DO i = 1, SIZE(y%bangle)
       ! matrix_pp type, uplo = 'U'
        ro_data%Lev1b%bangle_sigma(i) = SQRT(y%cov%d(i + i*(i-1)/2)) 
     ENDDO
  ENDIF

! 1.5 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_obs2dbangle2roprof
