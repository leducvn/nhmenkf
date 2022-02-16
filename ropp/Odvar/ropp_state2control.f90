! $Id: ropp_state2control.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_state2control(precon, state, control)

!****f* Matrix/ropp_state2control *
!
! NAME
!    ropp_state2control - Convert a state to a control variable.
!
! SYNOPSIS
!    call ropp_state2control(precon, state, control)
! 
! DESCRIPTION
!    This routine performs the variable transform from a control variable to
!    the original variable, given the preconditioner P.
!
! INPUTS
!    precon     preconditioner
!    state      variable in state space
!
! OUTPUT
!    control    variable in control space
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
!    ropp_controle2state
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

  TYPE(matrix_sq),        INTENT(in)  :: precon
  REAL(wp), DIMENSION(:), INTENT(in)  :: state
  REAL(wp), DIMENSION(:), INTENT(out) :: control

!-------------------------------------------------------------------------------
! 2. Do the variable transform
!-------------------------------------------------------------------------------

  control = MATMUL(precon % L_inv, state)

END SUBROUTINE ropp_state2control
