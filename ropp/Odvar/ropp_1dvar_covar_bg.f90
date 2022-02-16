! $Id: ropp_1dvar_covar_bg.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Error covariances/ropp_1dvar_covar_bg *
!
! NAME
!    ropp_1dvar_covar_bg - Set up covariance matrices for background data.
!
! SYNOPSIS
!    call ropp_1dvar_covar(bg,  bg_covar_config)
!    call ropp_1dvar_covar(obs, obs_covar_config)
! 
! DESCRIPTION
!    This subroutine sets up an error covariance matrix for a background
!    or an observation vector. It is an overloaded interface to various 
!    specialised error covariance routines.
!
! INPUTS
!    type(State1dFM) :: bg   Background state vector
!
! OUTPUT
!    The state vector's error covariance elements are
!    filled or updated according to the settings in the configuration
!    structures. 
!
! NOTES
!    Background covariances can be constructed using the following methods:
!
!       FSFC   Fixed Sigmas, Fixed correlations: Both error correlations and
!                and error standard deviations are read from a background
!                error correlation file. The error correlation file must 
!                contain both the error correlation matrix as well as the
!                standard deviations (errors) for all background state vector
!                elements.
!
!       VSFC   Variable Sigmas, Fixed Correlations: Error correlations are read
!                from an error correlation file, while per profile error 
!                estimates as contained in the background data file are used.
!                In this case, the error correlation / covariance data files
!                only require to contain the error correlations.
!
!    Note that error correlation files may contain latitudinally binned error
!    correlations and standard deviations, allowing to have latitudinally 
!    varying error correlation structures and standard deviations even in the 
!    FCFS scenario.
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
! 1. Background error covariances
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_covar_bg(bg, config)

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROcorcov
  USE ropp_fm
! USE ropp_1dvar, not_this => ropp_1dvar_covar_bg
  USE ropp_1dvar

  IMPLICIT NONE

  TYPE(State1dFM), INTENT(inout)        :: bg
  TYPE(VarConfig), INTENT(in)           :: config

  TYPE(ROcorcov), DIMENSION(:),   POINTER :: corcov => null()
  REAL(wp),       DIMENSION(:,:), POINTER :: sigma => null()
  REAL(wp),       DIMENSION(:,:), POINTER :: covar => null()
  INTEGER                                 :: i, idx_lat, n
  LOGICAL                                 :: file_exists
  CHARACTER(len = 256)                    :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_covar_bg')

!-------------------------------------------------------------------------------
! 2. Read error correlation and covariance data
!-------------------------------------------------------------------------------

  INQUIRE(file = config%bg_corr_file, exist = file_exists)
  IF (file_exists) THEN
     CALL ropp_io_read(corcov, config % bg_corr_file)
  ELSE
     CALL message(msg_error, &
          "Background error correlation / covariance file does not exist:\n")
     CALL message(msg_cont, &
        "   " // TRIM(config % bg_corr_file))
     bg%cov_ok = .FALSE.
     RETURN
  ENDIF

!-------------------------------------------------------------------------------
! 3. Get index of matching latitude bin
!-------------------------------------------------------------------------------

  IF (SIZE(corcov) == 1) THEN
     idx_lat = 1
  ELSE
     idx_lat = -1
     DO i = 1, SIZE(corcov)
        IF (bg%lat >= corcov(i)%lat_min .AND. bg%lat <= corcov(i)%lat_max) THEN
           idx_lat = i
        END IF
     END DO
     IF (idx_lat < 0) THEN
        CALL message(msg_error, &
             "No error correlation structure found for profile latitude.")
        bg%cov_ok = .FALSE.
        RETURN
     END IF
  END IF

!-------------------------------------------------------------------------------
! 4. Get standard error deviations and correlations based on specified method
!-------------------------------------------------------------------------------
  
  SELECT CASE (config % bg_covar_method)

! 4.1 Fixed sigmas, fixed correlations
! ------------------------------------

  CASE('FSFC')   ! Fixed correlations and standard deviations from error
                 !    correlation file

     IF (.NOT. ASSOCIATED(corcov(idx_lat)%sigma)) THEN
        CALL message(msg_fatal, &
             "Background standard error deviations requested from " //      &
             "error correlation file,\n" //                                 &
             "   but sigmas are not contained in data file.")
     END IF

     IF ( config%use_logp .OR. config%use_logq) THEN
       CALL message(msg_fatal, &
          "Support for FSFC bg error description with logq or logp " //     &
          "state vector elements not currently provided - use VSFC " //     &
          "method and define sigmas profile-by-profile")
     END IF

     n = SIZE(corcov(idx_lat) % sigma)

     IF (SIZE(bg % state) /= n) THEN
        CALL message(msg_error, &
             "Differing numbers of elements in background field " //        &
             "and background error correlation file.\n")
        bg%cov_ok = .FALSE.
        RETURN
     ENDIF

     CALL matrix_bm2full_alloc(RESHAPE(corcov(idx_lat) % sigma, (/1,n/)),   &
                               n, 0, sigma)

     ! Select covariance matrix with same dimensions as variable sigma
     CALL matrix_pp2full_alloc(corcov(idx_lat) % corr, covar)

     bg%cov_ok = .TRUE.
     
! 4.2 Variable sigmas, fixed correlations
! ---------------------------------------

  CASE('VSFC')   ! Fixed correlations from error correlation file, variable
                 !    standard deviations from background data

    IF (.NOT. bg % cov_ok) THEN
      CALL message(msg_error, &
                   "Background standard error deviations requested on a " // &
                   "per profile basis, \n" //  &
                   "   but background data file did not contain sigmas.")
      bg%cov_ok = .FALSE.
      RETURN
    END IF

    n = NINT((SQRT(8*SIZE(corcov(idx_lat) % corr) + 1.0_wp)-1)/2.0_wp)

    IF (SIZE(bg % state) /= n) THEN
      CALL message(msg_error, &
                   "Differing numbers of elements in background field " // &
                   "and background error correlation file.\n")
      bg%cov_ok = .FALSE.
      RETURN
    END IF

    CALL matrix_pp2full_alloc(bg % cov, sigma)
    sigma = SQRT(sigma)

    ! Select covariance matrix with same dimensions as variable sigma
    CALL matrix_pp2full_alloc(corcov(idx_lat) % corr, covar)

! 4.3 Variable sigmas, no (diagonal) correlations
! -----------------------------------------------

  CASE('VSDC')   ! Variable sigmas and diagonal error correlations

    IF (.NOT. bg % cov_ok) THEN
      CALL message(msg_error, & 
                   "Background standard error deviations requested on a " // &
                   "per profile basis,\n" // &
                   "   but background data file did not contain sigmas.")
      bg%cov_ok = .FALSE.
      RETURN
    END IF

    n = SIZE(bg % state)

    CALL matrix_bm2full_alloc(RESHAPE((/ (1.0_wp, i = 1,n)/), (/1,n/)), &
                              n, 0, covar) ! Unit matrix
    CALL matrix_pp2full_alloc(bg % cov, sigma)  
    sigma = SQRT(sigma)

! 4.4 Unknown error covariance method for background
! --------------------------------------------------

  CASE default

    CALL message(msg_error, &
                 "Background error covariance method " // &
                 TRIM(config % bg_covar_method) // " unknown.")
    bg%cov_ok = .FALSE.
    RETURN

  END SELECT

!-------------------------------------------------------------------------------
! 5. Construct error covariance matrix
!-------------------------------------------------------------------------------

  CALL matrix_toast(covar, sigma)

  CALL matrix_full2pp(covar, bg % cov)

!-------------------------------------------------------------------------------
! 6. Check for invertibility
!-------------------------------------------------------------------------------

  covar = matrix_invert(bg % cov)

!-------------------------------------------------------------------------------
! 7. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(covar)
  DEALLOCATE(sigma)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_covar_bg
