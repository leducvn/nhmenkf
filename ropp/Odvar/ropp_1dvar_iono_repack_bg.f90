! $Id: ropp_1dvar_iono_repack_bg.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_1dvar_iono_repack_bg *
!
! NAME
!    ropp_1dvar_iono_repack_bg - Extend and repack state vector and bg%cov 
!                                with {Ne_max, H_peak, H_width} if 
!                                modelling L1 and L2 directly.
!
! SYNOPSIS
!    CALL ropp_1dvar_iono_repack_bg(bg_data, bg, config)
!
! DESCRIPTION
!    This subroutine is invoked if the L1 and L2 bending angles are being 
!    modelled directly, via a model Chapman layer ionosphere 
!    (ie if -direct_ion is in force). It extends the state vector by
!    appending {Ne_max, H_peak, H_width} to it, and regenerates the suitably 
!    modified bg%cov%d and bg%cov%f structures.
!
! INPUTS
!    bg_data    ROprof structure
!    bg         State1dFM structure containing neutral state
!    config     Configuration structure
!
! OUTPUT
!    bg         Modified State1dFM structure containing augmented state and
!               covariances.
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_1dvar_iono_repack_bangle
!    ropp_1dvar_iono_unpack_bangle
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

SUBROUTINE ropp_1dvar_iono_repack_bg(bg_data, bg, config)

! 1.1 Declarations
! -----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_io_types
  USE ropp_fm
  USE ropp_fm_copy
  USE ropp_fm_iono
  USE ropp_1dvar_types

  TYPE(ROprof),      INTENT(inout)      :: bg_data
  TYPE(State1dFM),   INTENT(inout)      :: bg
  TYPE(VarConfig),   INTENT(in)         :: config

  TYPE(State1dFM)                       :: bg_neut

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: covar
  INTEGER                               :: imin, imax, n, n_tri
  CHARACTER(len=256)                    :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_iono_repack_bg')

! 1.3 Save the neutral covs
! -------------------------

  bg_neut = bg

! 1.4 Default the iono params if necessary
! ----------------------------------------

  CALL ropp_fm_iono_set_default(bg_data)

! 1.5 Extend the bg state
! -----------------------

  bg_data%lev2c%direct_ion = .TRUE.

  CALL ropp_fm_roprof2state(bg_data, bg) ! bg now has 2*n_lev + 4 elements

! 1.6 Augment the neutral covs with the 3 extra diagonal elements
! ---------------------------------------------------------------

  bg%cov_ok = bg_neut%cov_ok

  bg%cov%fact_chol = bg_neut%cov%fact_chol

  n = SIZE(bg_neut%state)  ;  n_tri = (n*(n+1)) / 2

  imin =        1  ;  imax = imin + n_tri - 1
  bg%cov%d(imin:imax) = bg_neut%cov%d(imin:imax)

  imin = imax + 1  ;  imax = imin + n - 1
  bg%cov%d(imin:imax) = 0.0_wp
  imin = imax + 1  ;  imax = imin
  bg%cov%d(imin:imax) = bg_data%lev2c%ne_max_sigma**2 ! Should already be set in call to ropp_fm_roprof2state

  imin = imax + 1  ;  imax = imin + n + 0
  bg%cov%d(imin:imax) = 0.0_wp
  imin = imax + 1  ;  imax = imin
  bg%cov%d(imin:imax) = bg_data%lev2c%h_peak_sigma**2 ! Should already be set in call to ropp_fm_roprof2state

  imin = imax + 1  ;  imax = imin + n + 1
  bg%cov%d(imin:imax) = 0.0_wp
  imin = imax + 1  ;  imax = imin
  bg%cov%d(imin:imax) = bg_data%lev2c%h_width_sigma**2 ! Should already be set in call to ropp_fm_roprof2state

! 1.7 Recalculate bg%cov%f, as ropp_1dvar_covar would do
! ------------------------------------------------------

  ALLOCATE(covar(n+3, n+3))

  IF (ASSOCIATED(bg%cov%f)) DEALLOCATE(bg%cov%f)

  covar = matrix_invert(bg%cov)  ! This sets bg%cov%f behind the scenes

  DEALLOCATE(covar)

! 1.8 Clear up
! ------------

  CALL ropp_fm_free(bg_neut)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_1dvar_iono_repack_bg

