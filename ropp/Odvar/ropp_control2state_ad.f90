! $Id: ropp_control2state_ad.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_control2state_ad(precon, control_ad, state_ad)

!****f* Preconditioning/ropp_control2state_ad *
!
! NAME
!    ropp_control2state_ad - Adjoint of ropp_control2state.
!
! SYNOPSIS
!    call ropp_control2state_ad(precon, control_ad, state_ad)
! 
! DESCRIPTION
!    This routine performs the variable transform from a control variable to
!    the original variable, given the preconditioner P.
!
! INPUTS
!    precon        preconditioner
!    state_ad      adjoint state variable
!
! OUTPUT
!    control_ad    adjoint control variable
!
! NOTES
!    The preconditioner must be set up previously by a call to preconditioner().
!
! EXAMPLE
!    use matrix
!    use ropp_1dvar
!      ...
!    type(matrix_pp) :: B
!    type(matrix_sq) :: precon
!      ...
!    call matrix_sqrt(B, precon)
!      ...
!    call ropp_state2control(precon, state, control)
!      ... [do something in control space] ...
!    call ropp_control2state(precon, control, state)
!
! SEE ALSO
!    matrix_sqrt
!    ropp_state2control
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
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix_types

  IMPLICIT NONE

  TYPE(matrix_sq),        INTENT(in)    :: precon
  REAL(wp), DIMENSION(:), INTENT(inout) :: control_ad
  REAL(wp), DIMENSION(:), INTENT(inout) :: state_ad

!-------------------------------------------------------------------------------
! 2. Do the variable transform
!-------------------------------------------------------------------------------

  control_ad = control_ad + MATMUL(TRANSPOSE(precon % L), state_ad)
  state_ad   = 0.0_wp

END SUBROUTINE ropp_control2state_ad
