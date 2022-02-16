! $Id: ropp_1dvar_iono_repack_bangle.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_1dvar_iono_repack_bangle *
!
! NAME
!    ropp_1dvar_iono_repack_bangle - Define obs%cov for 2m levels 
!                                    if modelling L1 and L2 directly.
!
! SYNOPSIS
!    CALL ropp_1dvar_iono_repack_bangle(obs_data, obs, config)
!
! DESCRIPTION
!    This subroutine is invoked if the L1 and L2 bending angles are being 
!    modelled directly, via a model Chapman layer ionosphere 
!    (ie if -direct_ion is in force). It regenerates obs%cov%d and obs%cov%f
!    for 2m levels, as appropriate to the concatentated obs vector = (bangle_L1, bangle_L2).
!
! INPUTS
!    obs_data   ROprof structure
!    obs        Obs1dBangle structure containing neutral covariances
!    config     Configuration structure
!
! OUTPUT
!    obs        Modified Obs1dBangle structure containing L1 and L2 covariances
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_1dvar_iono_unpack_bangle
!    ropp_1dvar_iono_repack_bg
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

SUBROUTINE ropp_1dvar_iono_repack_bangle(obs_data, obs, config)

! 1.1 Declarations
! -----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_io_types
  USE ropp_fm
  USE ropp_fm_copy
  USE ropp_1dvar_types

  TYPE(ROprof),      INTENT(inout)      :: obs_data
  TYPE(Obs1dBangle), INTENT(inout)      :: obs
  TYPE(VarConfig),   INTENT(in)         :: config

  TYPE(Obs1dBangle)                     :: obs_neut
!  REAL(wp), PARAMETER                   :: sigma_L2_to_sigma_L1=3.0_wp
  REAL(wp), PARAMETER                   :: sigma_L1_min=10.0e-6_wp ! rad
  REAL(wp), PARAMETER                   :: sigma_L2_min=30.0e-6_wp ! rad

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: covar
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: sigma_n, sigma_L1, sigma_L2
  INTEGER                               :: i, j, jpm, jmin, jmax, k, m, m_tri
  CHARACTER(len=256)                    :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_iono_repack_bangle')

! 1.3 Save the neutral covs
! -------------------------

  obs_neut = obs

! 1.4 Extend obs%cov to size 2m
! -----------------------------

  obs_data%lev2c%direct_ion = .TRUE.

  CALL ropp_fm_roprof2obs(obs_data, obs) ! obs now has 2m elements

  obs%cov_ok = obs_neut%cov_ok

  m = obs_neut%n_L1  ;  m_tri = (m*(m+1)) / 2

! 1.5 Extract and store the sigmas
! --------------------------------

  ALLOCATE(sigma_n(m), sigma_L1(m), sigma_L2(m))

  DO j=1,m
    sigma_n(j)  = SQRT(obs_neut%cov%d((j*(j+1))/2))
    sigma_L1(j) = SQRT(obs%cov%d((j*(j+1))/2))
    jpm = j + m
    sigma_L2(j) = SQRT(obs%cov%d((jpm*(jpm+1))/2))
  END DO

! 1.6 Repack the neutral covs
! ---------------------------

  jmin =        1  ;  jmax = jmin + m_tri - 1
  obs%cov%d(jmin:jmax) = obs_neut%cov%d(jmin:jmax)

  jmin = jmax + 1
  DO j=1,m  ! Not exactly vectorisable, but it's only done once
    jmax = jmin + m - 1
    obs%cov%d(jmin:jmax) = 0.0_wp
    jmin = jmax + 1
    jmax = jmin + j - 1
    obs%cov%d(jmin:jmax) = obs_neut%cov%d(((j*(j-1))/2)+1:(j*(j+1))/2)
    jmin = jmax + 1
  END DO

! 1.7 Rescale these repacked neutral bangle covs to get L1 and L2 covs
! --------------------------------------------------------------------

  IF (config%obs_covar_method == 'FSFC') THEN ! Derive sigma_L1/2 from sigma_n with larger minmim

    sigma_L1 = sigma_n
    WHERE (sigma_L1 < sigma_L1_min)
      sigma_L1 = sigma_L1_min
    END WHERE

    sigma_L2 = sigma_n
    WHERE (sigma_L2 < sigma_L2_min)
      sigma_L2 = sigma_L2_min
    END WHERE

  END IF

  DO j=1,m

    jmin = (j*(j-1))/2
    DO i=1,j
      k = jmin + i
      obs%cov%d(k) = obs%cov%d(k) * &
                     (sigma_L1(i) * sigma_L1(j) / sigma_n(i) / sigma_n(j))
    END DO

    jmin = ((j+m)*(j+m-1))/2
    DO i=1,j
      k = jmin + m + i
      obs%cov%d(k) = obs%cov%d(k) * &
                     (sigma_L2(i) * sigma_L2(j) / sigma_n(i) / sigma_n(j))
    END DO

  END DO

! 1.8 Recalculate obs%cov%f, as ropp_1dvar_covar would do
! -------------------------------------------------------

  ALLOCATE(covar(2*m, 2*m))

  covar = matrix_invert(obs%cov)  ! This sets obs%cov%f behind the scenes

  DEALLOCATE(covar)

! 1.9 Clear up
! -------------

  CALL ropp_fm_free(obs_neut)

  DEALLOCATE(sigma_n, sigma_L1, sigma_L2)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_iono_repack_bangle

