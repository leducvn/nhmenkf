! $Id: ropp_fm_roprof2obs.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Copying2/ropp_fm_roprof2obs *
!
! NAME
!    ropp_fm_roprof2obs - Copy elements of an ROprof structure to an
!                            observation vector.
!
! SYNOPSIS
!    type(ROprof)                 :: ro_data
!    type(<some obs vector type>) :: obs
!       ...
!    call ropp_fm_roprof2obs(ro_data, obs)
! 
! DESCRIPTION
!    This subroutine copies Level 1b (bending angle) or level 2a 
!    (refractivity) data as contained in a radio occultation profile data
!    structure to an observation vector.
!
! INPUTS
!   ro_data  Radio occultation profile data.
!
! OUTPUT
!   obs      An observation state vector (Obs1dBangle, Obs1dRefrac).
!
! NOTES
!   Data is copied from the ROprof data structure without unit conversion;
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
!   neutral atmospheric bending. Thus, the data is copied from the generic
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

SUBROUTINE ropp_fm_roprof2obs1dbangle(ro_data, y)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dBangle
! USE ropp_fm_copy,  not_this => ropp_fm_roprof2obs1dbangle
  USE geodesy

  IMPLICIT NONE

  TYPE(ROprof),         INTENT(in)    :: ro_data     ! RO data structure
  TYPE(Obs1dBangle),    INTENT(inout) :: y           ! Bending angle structure

  INTEGER                             :: i, n, ipn
  INTEGER, DIMENSION(8)               :: DT8

  CHARACTER(len = 256)                :: routine
  CHARACTER(len = 10)                 :: err_val

! 1.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2obs (1D bending angles)')
  
  y%obs_ok = .TRUE.

! 1.3 Check and copy geolocation and time
! ---------------------------------------

  IF (isinrange(ro_data%georef%lon, ro_data%georef%range%lon)) THEN
    y%lon = ro_data%georef%lon
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lon
    CALL message(msg_warn, &
         "Longitude data for observations out of range or missing. " // &
         "(longitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

  IF (isinrange(ro_data%georef%lat, ro_data%georef%range%lat)) THEN
    y%lat = ro_data%georef%lat
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lat 
    CALL message(msg_warn, &
         "Latitude data for observations out of range or missing. " // &
         "(latitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

  IF ( isinrange(ro_data%DTocc%year,   ro_data%DTocc%range%year)   .AND. &
       isinrange(ro_data%DTocc%month,  ro_data%DTocc%range%month)  .AND. &
       isinrange(ro_data%DTocc%day,    ro_data%DTocc%range%day)    .AND. &
       isinrange(ro_data%DTocc%hour,   ro_data%DTocc%range%hour)   .AND. &
       isinrange(ro_data%DTocc%minute, ro_data%DTocc%range%minute) .AND. &
       isinrange(ro_data%DTocc%second, ro_data%DTocc%range%second) ) THEN
    DT8 = (/ro_data%DTocc%year,   ro_data%DTocc%month,  &
            ro_data%DTocc%day,    0,                    &
            ro_data%DTocc%hour,   ro_data%DTocc%minute, &
            ro_data%DTocc%second, ro_data%DTocc%msec/)
    CALL TimeSince ( DT8, y%time, 1, Base="JS2000" )
  ELSE
    CALL message(msg_warn, &
         "Time data for observations out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

! 1.4 Check that profiles are increasing in height - 1st element towards surface
! ------------------------------------------------

  CALL ropp_io_ascend(ro_data)

! 1.5 Check and copy observation data
! -----------------------------------

  IF (ro_data%Lev1b%Npoints == 0) THEN
    CALL message(msg_warn, &
         "RO data has no level 1b (bending angle) parameters.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
    RETURN
  ENDIF

  n = ro_data%Lev1b%Npoints

  IF (.NOT. ro_data%Lev2c%direct_ion) THEN

    ALLOCATE(y%impact(n), y%bangle(n), y%weights(n))
    y%impact           = ro_data%Lev1b%impact(1:n)
    y%bangle           = ro_data%Lev1b%bangle(1:n)
    y%weights          = 1.0_wp
    y%n_L1             = n

  ELSE

    ALLOCATE(y%impact(2*n), y%bangle(2*n), y%weights(2*n))

    IF (ALL(ro_data%Lev1b%impact_L1 < ropp_MDTV)) THEN ! Use neutral IP if no L1 IPs
      y%impact(1:n)      = ro_data%Lev1b%impact(1:n)
    ELSE
      y%impact(1:n)      = ro_data%Lev1b%impact_L1(1:n)
    END IF

    IF (ALL(ro_data%Lev1b%impact_L2 < ropp_MDTV)) THEN ! Use neutral IP if no L2 IPs
      y%impact(n+1:2*n)  = ro_data%Lev1b%impact(1:n)
    ELSE
      y%impact(n+1:2*n)  = ro_data%Lev1b%impact_L2(1:n)
    END IF

    y%bangle(1:n)      = ro_data%Lev1b%bangle_L1(1:n)
    y%bangle(n+1:2*n)  = ro_data%Lev1b%bangle_L2(1:n)
    y%weights          = 1.0_wp
    y%n_L1             = n

  END IF

  y%g_sfc              = gravity(ro_data%georef%lat, 0.0_wp)

  y%r_earth            = R_eff(ro_data%georef%lat)

  IF (isinrange(ro_data%georef%roc, ro_data%georef%range%roc)) THEN
    y%r_curve      = ro_data%georef%roc
  ELSE
    CALL message(msg_warn, "Radius of curvature out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

  IF (isinrange(ro_data%georef%undulation, ro_data%georef%range%undulation)) THEN
    y%undulation   = ro_data%georef%undulation 
  ELSE
    CALL message(msg_warn, "Undulation out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

! 1.6 Check and copy sigmas to diagonal error covariance
! ------------------------------------------------------

  y%cov_ok = .TRUE.

  IF (ASSOCIATED(y%cov%d)) DEALLOCATE(y%cov%d)

  IF (.NOT. ro_data%Lev2c%direct_ion) THEN

    CALL callocate(y%cov%d, n*(n+1)/2)

    DO i = 1, n
      IF (ro_data%Lev1b%bangle(i) > 0.0_wp) THEN
        IF (ro_data%Lev1b%bangle_sigma(i) > 0.0_wp) THEN
         ! matrix_pp type, uplo = 'U'
          y%cov%d(i + i*(i-1)/2) = ro_data%Lev1b%bangle_sigma(i)**2
        ELSE
          y%cov_ok = .FALSE.
        END IF
      ELSE
        y%cov%d(i + i*(i-1)/2) = 0.0003_wp  ! = 1 deg**2
      ENDIF
    END DO

  ELSE

    CALL callocate(y%cov%d, n*(2*n+1))

    DO i = 1, n

      IF (ro_data%Lev1b%bangle_L1(i) > 0.0_wp) THEN
        IF (ro_data%Lev1b%bangle_L1_sigma(i) > 0.0_wp) THEN
         ! matrix_pp type, uplo = 'U'
          y%cov%d(i + i*(i-1)/2) = ro_data%Lev1b%bangle_L1_sigma(i)**2
        ELSE
          y%cov_ok = .FALSE.
        END IF
      ELSE
        y%cov%d(i + i*(i-1)/2) = 0.0003_wp  ! = 1 deg**2
      END IF

      ipn = i + n
      IF (ro_data%Lev1b%bangle_L2(i) > 0.0_wp) THEN
        IF (ro_data%Lev1b%bangle_L2_sigma(i) > 0.0_wp) THEN
         ! matrix_pp type, uplo = 'U'
          y%cov%d(ipn + ipn*(ipn-1)/2) = ro_data%Lev1b%bangle_L2_sigma(i)**2
        ELSE
          y%cov_ok = .FALSE.
        END IF
      ELSE
        y%cov%d(ipn + ipn*(ipn-1)/2) = 0.0003_wp  ! = 1 deg**2
      END IF

    END DO

  END IF

  IF (ASSOCIATED(y%cov%e)) DEALLOCATE(y%cov%e)
  IF (ASSOCIATED(y%cov%f)) DEALLOCATE(y%cov%f)
  IF (ASSOCIATED(y%cov%s)) DEALLOCATE(y%cov%s)

  y%cov%fact_chol = .FALSE.
  y%cov%equi_chol = 'N'

! 1.7 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2obs1dbangle


!-------------------------------------------------------------------------------
! 2. 1D Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_roprof2obs1drefrac(ro_data, y)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dRefrac
! USE ropp_fm_copy,  not_this => ropp_fm_roprof2obs1drefrac

  IMPLICIT NONE

  TYPE(ROprof),          INTENT(in)    :: ro_data    ! RO data structure
  TYPE(Obs1dRefrac),     INTENT(inout) :: y          ! Refractivity structure

  INTEGER                              :: i, n
  INTEGER, DIMENSION(8)                :: DT8
  
  CHARACTER(len = 256)                 :: routine
  CHARACTER(len = 10)                  :: err_val

! 2.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2obs (1D refractivity)')

  y%obs_ok = .TRUE.

! 2.3 Check and copy geolocation and time
! ---------------------------------------

  IF (isinrange(ro_data%georef%lon, ro_data%georef%range%lon)) THEN
    y%lon = ro_data%georef%lon
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lon
    CALL message(msg_warn, &
         "Longitude data for observations out of range or missing. " // &
         "(longitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
  ENDIF

  IF (isinrange(ro_data%georef%lat, ro_data%georef%range%lat)) THEN
    y%lat = ro_data%georef%lat
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lat
    CALL message(msg_warn, &
         "Latitude data for observations out of range or missing. " // &
         "(latitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
  ENDIF

  IF ( isinrange(ro_data%DTocc%year,   ro_data%DTocc%range%year)   .AND. &
       isinrange(ro_data%DTocc%month,  ro_data%DTocc%range%month)  .AND. &
       isinrange(ro_data%DTocc%day,    ro_data%DTocc%range%day)    .AND. &
       isinrange(ro_data%DTocc%hour,   ro_data%DTocc%range%hour)   .AND. &
       isinrange(ro_data%DTocc%minute, ro_data%DTocc%range%minute) .AND. &
       isinrange(ro_data%DTocc%second, ro_data%DTocc%range%second) ) THEN
    DT8 = (/ro_data%DTocc%year,   ro_data%DTocc%month,  &
            ro_data%DTocc%day,    0,                    &
            ro_data%DTocc%hour,   ro_data%DTocc%minute, &
            ro_data%DTocc%second, ro_data%DTocc%msec/)
    CALL TimeSince ( DT8, y%time, 1, Base="JS2000" )
  ELSE
    CALL message(msg_warn, &
         "Time data for observations out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag state%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF
  
! 2.4 Check that profiles are increasing in height - 1st element towards surface
! ------------------------------------------------
  
  CALL ropp_io_ascend(ro_data)

! 2.5 Check and copy observation data
! -----------------------------------

  IF (ro_data%Lev2a%Npoints == 0) THEN
    CALL message(msg_warn, &
         "RO data has no Level 2a (refractivity) data.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag state%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF

  n = ro_data%Lev2a%Npoints
  
  ALLOCATE(y%geop(n), y%refrac(n), y%weights(n))
  y%geop    = ro_data%Lev2a%geop_refrac(1:n)
  y%refrac  = ro_data%Lev2a%refrac(1:n)
  y%weights = 1.0_wp

! 2.6 Check and copy sigmas to diagonal error covariance
! ------------------------------------------------------

  y%cov_ok = .TRUE.

  IF (ASSOCIATED(y%cov%d)) DEALLOCATE(y%cov%d)
  CALL callocate(y%cov%d, (n*(n+1)/2))

  DO i = 1, n
    IF (ro_data%Lev2a%refrac(i) > 0.0_wp) THEN
      IF (ro_data%Lev2a%refrac_sigma(i) > 0.0_wp) THEN
       ! matrix_pp type, uplo = 'U'
        y%cov%d(i + i*(i-1)/2) = ro_data%Lev2a%refrac_sigma(i)**2
      ELSE
        y%cov_ok = .FALSE.
      END IF
    ELSE 
      y%cov%d(i + i*(i-1)/2) = 0.0003_wp
    ENDIF
  END DO

  IF (ASSOCIATED(y%cov%e)) DEALLOCATE(y%cov%e)
  IF (ASSOCIATED(y%cov%f)) DEALLOCATE(y%cov%f)
  IF (ASSOCIATED(y%cov%s)) DEALLOCATE(y%cov%s)

  y%cov%fact_chol = .FALSE.
  y%cov%equi_chol = 'N'

! 2.7 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2obs1drefrac


!-------------------------------------------------------------------------------
! 3. 2D Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_roprof2obs2dbangle(ro_data, y)

! 3.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof2d
  USE ropp_fm
  USE ropp_fm_types, ONLY: Obs1dBangle
! USE ropp_fm_copy,  not_this => ropp_fm_roprof2obs2dbangle
  USE geodesy

  IMPLICIT NONE

  TYPE(ROprof2d),         INTENT(in)  :: ro_data
  TYPE(Obs1dBangle),    INTENT(inout) :: y

  INTEGER                             :: n
  INTEGER, DIMENSION(8)               :: DT8

  CHARACTER(len = 256)                :: routine
  CHARACTER(len = 10)                 :: err_val

! 3.2 Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2obs (2D bending angles)')
  
  y%obs_ok = .TRUE.

! 3.3 Check and copy geolocation and time
! ---------------------------------------

  IF (isinrange(ro_data%georef%lon, ro_data%georef%range%lon)) THEN
    y%lon = ro_data%georef%lon
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lon
    CALL message(msg_warn, &
         "Longitude data for observations out of range or missing. " // &
         "(longitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    y%obs_ok = .FALSE.
  ENDIF

  IF (isinrange(ro_data%georef%lat, ro_data%georef%range%lat)) THEN
    y%lat = ro_data%georef%lat
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lat
    CALL message(msg_warn, &
         "Latitude data for observations out of range or missing. " // &
         "(latitude value = " // TRIM(err_val) // ")")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    y%obs_ok = .FALSE.
  ENDIF

  IF ( isinrange(ro_data%DTocc%year,   ro_data%DTocc%range%year)   .AND. &
       isinrange(ro_data%DTocc%month,  ro_data%DTocc%range%month)  .AND. &
       isinrange(ro_data%DTocc%day,    ro_data%DTocc%range%day)    .AND. &
       isinrange(ro_data%DTocc%hour,   ro_data%DTocc%range%hour)   .AND. &
       isinrange(ro_data%DTocc%minute, ro_data%DTocc%range%minute) .AND. &
       isinrange(ro_data%DTocc%second, ro_data%DTocc%range%second)) THEN
    DT8 = (/ro_data%DTocc%year,   ro_data%DTocc%month,  &
            ro_data%DTocc%day,    0,                    &
            ro_data%DTocc%hour,   ro_data%DTocc%minute, &
            ro_data%DTocc%second, ro_data%DTocc%msec/)
    CALL TimeSince ( DT8, y%time, 1, Base="JS2000" )
  ELSE
    CALL message(msg_warn, &
         "Time data for observations is out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag state%ok to FALSE")
    y%obs_ok = .false.
  ENDIF

! 3.4 Check that profiles are increasing in height - 1st element towards surface
! ------------------------------------------------

!  call ropp_io_ascend(ro_data)

! 3.5 Check and copy observation data
! -----------------------------------

  IF (ro_data%Lev1b%Npoints == 0) THEN
    CALL message(msg_warn, &
         "RO data has no level 1b (bending angle) parameters.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag state%ok to FALSE")
    y%obs_ok = .FALSE.
    RETURN
  ENDIF
  
  n = ro_data%Lev1b%Npoints
  
! store number of bending angles
  
  y%nobs = n  
  
! allocate and set impact parameter and bending angle
  
  ALLOCATE(y%impact(n))
  y%impact  = ro_data%Lev1b%impact(1:n)
  ALLOCATE(y%bangle(n))
  y%bangle  = ro_data%Lev1b%bangle(1:n)
  ALLOCATE(y%weights(n))
  y%weights = 1.0_wp
  
  ALLOCATE(y%a_path(n,2))
  ALLOCATE(y%rtan(n))

! set the radius of curv, undulation and azimthuthal angle

  IF (isinrange(ro_data%georef%roc, ro_data%georef%range%roc)) THEN
    y%r_curve      = ro_data%georef%roc 
  ELSE
    CALL message(msg_warn, "Radius of curvature out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF
  IF (isinrange(ro_data%georef%undulation, ro_data%georef%range%undulation)) THEN
    y%undulation   = ro_data%georef%undulation 
  ELSE
    CALL message(msg_warn, "Undulation out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF
  IF (isinrange(ro_data%georef%azimuth, ro_data%georef%range%azimuth)) THEN
    y%azimuth = ro_data%georef%azimuth 
  ELSE
    CALL message(msg_warn, "Azimuth out of range or missing.")
    CALL message(msg_warn, "Check input data and valid_range attributes")
    CALL message(msg_warn, "Set status flag obs%ok to FALSE")
    y%obs_ok = .FALSE.
  ENDIF


! 3.7 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2obs2dbangle
