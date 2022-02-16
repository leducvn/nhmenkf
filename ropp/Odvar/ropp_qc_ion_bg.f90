! $Id: ropp_qc_ion.f90 2280 2009-10-21 16:32:29Z frhl $

!****s* QC/ropp_qc_ion_bg *
!
! NAME
!    ropp_qc_ion_bg - Check whether the background ionospheric parameters 
!                     are appropriate for the L1 and L2 bending angle profiles,
!                     if we're modelling them directly using -direct_ion.
!                     If not, reset them.
!
! SYNOPSIS
!    CALL ropp_qc_ion_bg(bg, obs, config, diag)
!
! DESCRIPTION
!    This subroutine resets ionospheric parameters that are inconsistent 
!    with the bending angles.
!
! INPUTS
!    bg          State1dFM structure, comprising usual state augmented by
!                {Ne_max, H_peak, H_width}
!    obs         Obs1dBangle structure, containing {L1 and L2}
!    config      Configuration options (not currently used)
!
! OUTPUT
!    bg     Updated State1dFM structure
!    diag        Diagnostics structure (not currently used)
!
! NOTES
!    At present checking is only applied to Ne_max.
!    Config and diag are included for possible future use.
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_qc_bgqc
!    ropp_qc_ion_bangle
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
! 1. Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_ion_bg(bg, obs, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_1dvar_types
! USE ropp_qc, not_this => ropp_qc_ion_bg

  IMPLICIT NONE

  TYPE(State1dFM), INTENT(inout)      :: bg
  TYPE(Obs1dBangle), INTENT(in)       :: obs
  TYPE(VarConfig), INTENT(in)         :: config
  TYPE(VarDiag), INTENT(inout)        :: diag

  REAL(wp), PARAMETER                 :: bangle_diff_to_Ne_max=2.4e16_wp
  REAL(wp), PARAMETER                 :: min_fit_height=30.0e3_wp
  REAL(wp)                            :: avg_l2_minus_l1
  INTEGER                             :: j1, j2, m_both

  CHARACTER(len = 256)                :: routine
  CHARACTER(len = 7)                  :: str_const

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_ion_bg')

! 1.3 Fit background Ne_max to <L2-L1>
! ------------------------------------

  avg_l2_minus_l1 = 0.0_wp

  m_both = 0

  DO j1=1,obs%n_L1

    j2 = j1 + obs%n_L1

    IF ((obs%bangle(j1)  > ropp_MDTV) .AND. &
        (obs%bangle(j2)  > ropp_MDTV) .AND. &
        (obs%weights(j1) > ropp_ZERO) .AND. &
        (obs%weights(j2) > ropp_ZERO) .AND. &
        ((obs%impact(j1)-obs%r_curve) > min_fit_height)) THEN

      avg_l2_minus_l1 = avg_l2_minus_l1 + (obs%bangle(j2) - obs%bangle(j1))

      m_both = m_both + 1

    END IF

  END DO

  IF (m_both > 0) THEN

    avg_l2_minus_l1 = avg_l2_minus_l1 / REAL(m_both, wp)

    bg%ne_max = avg_l2_minus_l1 * bangle_diff_to_Ne_max

    WRITE (str_const, '(F7.2)') bg%ne_max*1.0e-9_wp

    CALL message(msg_info, "L2-L1 data cause background Ne_max to be " // &
                 "set to " // TRIM(ADJUSTL(str_const)) // "e9 m-3.\n")

  ELSE

    WRITE (str_const, '(F7.2)') bg%ne_max*1.0e-9_wp

    CALL message(msg_info, "Not enough L2-L1 data to intitialise Ne_max ... " // &
                 "leaving as " // TRIM(ADJUSTL(str_const)) // "e9 m-3.\n")

  END IF

! 1.4 Copy to state part of bg structure
! --------------------------------------

  bg%state(SIZE(bg%state)-2) = bg%ne_max

! 1.5 Recalculate OmB with the revised bg%ne_max
! ----------------------------------------------

  CALL ropp_qc_OmB(obs, bg, config, diag)

! 1.6 Return to sender
! --------------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_ion_bg
