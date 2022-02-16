! $Id: ropp_1dvar_covar_bangle.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Error covariances/ropp_1dvar_covar_bangle *
!
! NAME
!    ropp_1dvar_covar_bangle - Set up covariance matrices for a bending angle
!                                observation profile.
!
! SYNOPSIS
!    call ropp_1dvar_covar(obs, obs_covar_config)
! 
! DESCRIPTION
!    This subroutine sets up an error covariance matrix for a vector of bending
!    angle observations. This particular routine provides the implementation for
!    an overloaded interface to various background and observation data types.
!
!
! INPUTS
!    type(Obs1dBangle)   :: obs  Bending angle observation vector
!
! OUTPUT
!    The observation vector's error covariance elements are
!    filled or updated according to the settings in the configuration
!    structures. 
!
! NOTES
!    Bending angle error covariances can be constructed using the following 
!    methods:
!
!       FSFC   Fixed Sigmas, Fixed correlations: Both error correlations and
!                and error standard deviations are read from an observation
!                error correlation file. The error correlation file must 
!                contain both the error correlation matrix as well as the
!                standard deviations (errors) for all observation vector
!                elements.
!
!       VSDC   Variable Sigmas, Diagonal Correlations: A diagonal error 
!                correlation structure (i.e., no error correlations) is assumed,
!                while per profile error estimates as contained in the 
!                observation data file are used. In this case, no error 
!                correlation / covariance data file is required.
!
!       VSFC   Variable Sigmas, Fixed Correlations: Error correlations are read
!                from an error correlation file, while per profile error 
!                estimates as contained in the observation data file are used.
!                In this case, the error correlation / covariance data files
!                only require to contain the error correlations.
!
!    Note that error correlation files may contain latitudinally binned error
!    correlations and standard deviations, allowing to have latitudinally 
!    varying error correlation structures and standard deviations in the FSFC  
!    and VSFC scenarios.
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_1dvar_covar_bg
!    ropp_1dvar_covar_refrac
!    ropp_1dvar_covar_bangle
!
! REFERENCES
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
! 1. Bending angle error covariances
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_covar_bangle(obs, config)

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROcorcov
  USE ropp_fm, Pi_loc => Pi
! USE ropp_1dvar, not_this => ropp_1dvar_covar_bangle
  USE ropp_1dvar

  IMPLICIT NONE

  TYPE(Obs1DBangle), INTENT(inout)      :: obs
  TYPE(VarConfig),   INTENT(in)         :: config

  TYPE(ROcorcov), DIMENSION(:),   POINTER :: corcov => null()
  REAL(wp),       DIMENSION(:,:), POINTER :: sigma => null()
  REAL(wp),       DIMENSION(:,:), POINTER :: covar => null()
  INTEGER                                 :: i, idx_lat, n
  LOGICAL                                 :: file_exists
  CHARACTER(len = 256)                    :: routine
  REAL(wp)                                :: season
  REAL(wp)                                :: season_scale_factor

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_covar_bangle')

!-------------------------------------------------------------------------------
! 2. Error correlations from file (if required)
!-------------------------------------------------------------------------------

  IF (config % obs_covar_method == "FSFC" .OR. &
      config % obs_covar_method == "VSFC") THEN

!    2.1 Read the error correlation / covariance data
!    ------------------------------------------------

     INQUIRE(file = config%obs_corr_file, exist = file_exists)
     IF (file_exists) THEN
        CALL ropp_io_read(corcov, config % obs_corr_file)
     ELSE
        CALL message(msg_error, &
          "Bending angle error correlation / covariance file does not exist:\n")
        CALL message(msg_cont,  &
             "   " // TRIM(config % obs_corr_file))
        obs%cov_ok = .false.
        RETURN
     ENDIF

!    2.2 Obtain latitude index
!    -------------------------

     IF (SIZE(corcov) == 1) THEN
        idx_lat = 1
     ELSE
        idx_lat = -1
        DO i = 1, SIZE(corcov)
           IF (obs%lat >= corcov(i)%lat_min .AND.    &
               obs%lat <= corcov(i)%lat_max) THEN
              idx_lat = i
           END IF
        END DO
        IF (idx_lat < 0) THEN
           CALL message(msg_error, &
                "No error correlation structure found for profile latitude.")
           obs%cov_ok = .false.
           RETURN
        END IF
     END IF

  END IF

!-------------------------------------------------------------------------------
! 3. Get standard error deviations and correlations based on specified method
!-------------------------------------------------------------------------------

  SELECT CASE (config % obs_covar_method)

! 3.1 Fixed sigmas, fixed correlations
! ------------------------------------

  CASE('FSFC')   ! Fixed sigmas and error correlations from correlation file

     IF (.NOT. ASSOCIATED(corcov(idx_lat)%sigma)) THEN
        CALL message(msg_error, &
          "Bending angle standard error deviations requested from \n" // &
          "error correlation file, but sigmas are not contained in data file.")
        obs%cov_ok = .false.
        RETURN
     END IF

     n = SIZE(corcov(idx_lat) % sigma)

     IF (SIZE(obs % bangle) > n) THEN
        CALL message(msg_error, &
             "Not enough elements in the bending angle observation " //      &
             "error correlation file.\n")
        obs%cov_ok = .FALSE.
        RETURN
     ENDIF

     IF (SIZE(obs % bangle) < n) THEN
        CALL message(msg_warn, &
             "Too many elements in the bending angle observation " //        &
             "error correlation file.\n" //                                 &
             "This will work, but is probably not what you wanted.\n")
     ENDIF

     CALL matrix_pp2full_alloc(corcov(idx_lat) % corr, covar)
     CALL matrix_bm2full_alloc(RESHAPE(corcov(idx_lat)%sigma,(/1,n/)),n,0,sigma)
     obs%cov_ok = .true.

! 3.2 Variable sigmas, no (diagonal) correlations
! -----------------------------------------------

  CASE('VSDC')   ! Variable sigmas and diagonal error correlations

     IF (.NOT. obs % cov_ok) THEN
        CALL message(msg_error, &
         "Bending angle standard error deviations requested on a \n" // &
         "per profile basis, but observation data file did not contain sigmas.")
        obs%cov_ok = .false.
        RETURN
     END IF

     n = SIZE(obs % bangle)

     ! Unit matrix
     CALL matrix_bm2full_alloc(RESHAPE((/ (1.0_wp, i=1,n)/),(/1,n/)),n,0,covar)
     CALL matrix_pp2full_alloc(obs % cov, sigma)
     sigma = SQRT(sigma)
        
! 3.3 Variable sigmas, fixed correlations
! ---------------------------------------

  CASE('VSFC')   ! Variable sigmas from profile, fixed corrns from error correlation file

     IF (.NOT. obs % cov_ok) THEN
        CALL message(msg_error, &
         "Bending angle standard error deviations requested on a \n" // &
         "per profile basis, but observation data file did not contain sigmas.")
        obs%cov_ok = .false.
        RETURN
     END IF

     n = NINT((SQRT(8*SIZE(corcov(idx_lat) % corr) + 1.0_wp)-1)/2.0_wp)

     IF (SIZE(obs % bangle) > n) THEN
        CALL message(msg_error, &
             "Not enough elements in the bending angle observation " //      &
             "error correlation file.\n")
        obs%cov_ok = .FALSE.
        RETURN
     ENDIF

     IF (SIZE(obs % bangle) < n) THEN
        CALL message(msg_warn, &
             "Too many elements in the bending angle observation " //        &
             "error correlation file.\n" //                                 &
             "This will work, but is probably not what you wanted.\n")

     ENDIF

     ! Variable sigmas
     CALL matrix_pp2full_alloc(obs % cov, sigma)
     sigma = SQRT(sigma)

     ! Select covariance matrix with same dimensions as variable sigma
     ALLOCATE(covar(SIZE(sigma,1),SIZE(sigma,1)))
     CALL matrix_pp2full_subset(corcov(idx_lat) % corr, covar)

! 3.4 Unknown error covariance method for observations
! ----------------------------------------------------

  CASE default

     CALL message(msg_error, &
          "Bending angle error covariance method " // &
          TRIM(config % obs_covar_method) // " unknown.")
     obs%cov_ok = .false.
     RETURN

  END SELECT

!-------------------------------------------------------------------------------
! 4. Perform seasonal scaling of observation errors (sigmas)
!-------------------------------------------------------------------------------

  ! Find fraction of way through year (obs%time is the number of seconds since
  ! 1/1/2000, 00Z)
  season = MOD(obs % time / (86400.0_wp), 365.25_wp) / 365.25_wp

  ! scale observation sigma values according to specified parameters
  season_scale_factor = 1.0_wp + config % season_amp * &
              COS(2.0_wp * Pi_loc * (season + config % season_phase)) + &
              config % season_offset

  IF (MINVAL(season_scale_factor * sigma) .GE. 0.0_wp) THEN
    sigma = season_scale_factor * sigma
  ELSE
    CALL message(msg_warn,     &
             'Observation sigma vector contains negative values,\n'// &
             'probably from seasonal scaling. Reverting to unscaled values. \n')
  END IF

!-------------------------------------------------------------------------------
! 5. Construct error covariance matrix
!-------------------------------------------------------------------------------

  CALL matrix_toast(covar, sigma)

  CALL matrix_full2pp(covar, obs % cov)

!-------------------------------------------------------------------------------
! 6. Check for invertibility
!-------------------------------------------------------------------------------

  covar = matrix_invert(obs % cov)

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(covar)
  DEALLOCATE(sigma)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_covar_bangle
