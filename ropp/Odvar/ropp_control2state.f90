! $Id: ropp_control2state.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_control2state(precon, control, state)

!****f* Preconditioning/ropp_control2state *
!
! NAME
!    control2state - Convert the control variable back to a state.
!
! SYNOPSIS
!    call ropp_control2state(precon, control, state)
! 
! DESCRIPTION
!    This routine performs the variable transform from a control variable to
!    the original variable, given the preconditioner P.
!
! INPUTS
!    precon     preconditioner
!    control    variable in control space
!
! OUTPUT
!    state      variable in state space
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

  TYPE(matrix_sq),        INTENT(in)  :: precon
  REAL(wp), DIMENSION(:), INTENT(in)  :: control
  REAL(wp), DIMENSION(:), INTENT(out) :: state

!-------------------------------------------------------------------------------
! 2. Do the variable transform
!-------------------------------------------------------------------------------

  state = MATMUL(precon % L, control)

END SUBROUTINE ropp_control2state
