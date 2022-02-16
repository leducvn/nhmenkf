! $Id: ropp_fm_iono_set_default.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_fm_iono_set_default *
!
! NAME
!    ropp_fm_iono_set_default - Set default values of
!                               {Ne_max, H_peak, H_width} if
!                               modelling L1 and L2 directly.
!
! SYNOPSIS
!    CALL  ropp_fm_iono_set_default(ro_data)
!
! DESCRIPTION
!    This subroutine is invoked if the L1 and L2 bending angles are being
!    modelled directly, via a model Chapman layer ionosphere 
!    (ie if -direct_ion is in force). It defines {Ne_max, H_peak, H_width}
!    if they are missing.
!
! INPUTS
!    ro_data    ROprof structure
!
! OUTPUT
!    ro_data    Modified ROprof structure containing non-missing
!               ionospheric parameters.
!
! NOTES
!   Negative values of ne_max can be less than ropp_MDTV, so we need
!   another test for missing data in this case.
!
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_fm_bg2ro_1d
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

SUBROUTINE ropp_fm_iono_set_default(ro_data)

! 1.1 Declarations
! -----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io_types
  USE ropp_fm
  USE ropp_fm_copy

  TYPE(ROprof),      INTENT(inout)      :: ro_data

  REAL(wp), PARAMETER                   :: ne_max_default        = 300.0e9_wp ! m-3
  REAL(wp), PARAMETER                   :: ne_max_sigma_default  = 200.0e9_wp ! m-3
  REAL(wp), PARAMETER                   :: h_peak_default        = 300.0e3_wp ! m
  REAL(wp), PARAMETER                   :: h_peak_sigma_default  = 150.0e3_wp ! m
  REAL(wp), PARAMETER                   :: h_width_default       =  75.0e3_wp ! m
  REAL(wp), PARAMETER                   :: h_width_sigma_default =  25.0e3_wp ! m

  CHARACTER(len=5)                      :: str_const
  CHARACTER(len=256)                    :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_iono_set_default')

! 1.3 Default the iono params
! ---------------------------

! (Note that negative values of ne_max can be less than ropp_MDTV, so we need
!  another test for missing data in this case.)
  IF (ABS(ro_data%lev2c%ne_max - ropp_MDFV) < 1.0_wp) THEN
    WRITE (str_const, '(F5.1)') ne_max_default*1.0e-9_wp
    CALL message(msg_info, &
                 "-direct_ion requested but Ne_max not set ... " // &
                 "defaulting to " // str_const // "e9 m-3.")
    ro_data%lev2c%ne_max = ne_max_default
  END IF

  IF (ro_data%lev2c%ne_max_sigma < ropp_MDTV) THEN
    WRITE (str_const, '(F5.1)') ne_max_sigma_default*1.0e-9_wp
    CALL message(msg_info, &
                 "-direct_ion requested but sigma(Ne_max) not set ... " // &
                 "defaulting to " // str_const // "e9 m-3.")
    ro_data%lev2c%ne_max_sigma = ne_max_sigma_default
  END IF

  IF (ro_data%lev2c%h_peak < ropp_MDTV) THEN
    WRITE (str_const, '(F5.1)') h_peak_default*1.0e-3_wp
    CALL message(msg_info, &
                 "-direct_ion requested but H_peak not set ... " // &
                 "defaulting to " // str_const // " km.")
    ro_data%lev2c%h_peak = h_peak_default
  END IF

  IF (ro_data%lev2c%h_peak_sigma < ropp_MDTV) THEN
    WRITE (str_const, '(F5.1)') h_peak_sigma_default*1.0e-3_wp
    CALL message(msg_info, &
                 "-direct_ion requested but sigma(H_peak) not set ... " // &
                 "defaulting to " // str_const // " km.")
    ro_data%lev2c%h_peak_sigma = h_peak_sigma_default
  END IF

  IF (ro_data%lev2c%h_width < ropp_MDTV) THEN
    WRITE (str_const, '(F5.1)') h_width_default*1.0e-3_wp
    CALL message(msg_info, &
                 "-direct_ion requested but H_width not set ... " // &
                 "defaulting to " // str_const // " km.")
    ro_data%lev2c%h_width = h_width_default
  END IF

  IF (ro_data%lev2c%h_width_sigma < ropp_MDTV) THEN
    WRITE (str_const, '(F5.1)') h_width_sigma_default*1.0e-3_wp
    CALL message(msg_info, &
                 "-direct_ion requested but sigma(H_width) not set ... " // &
                 "defaulting to " // str_const // " km.")
    ro_data%lev2c%h_width_sigma = h_width_sigma_default
  END IF

! 1.4 Clear up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_iono_set_default

