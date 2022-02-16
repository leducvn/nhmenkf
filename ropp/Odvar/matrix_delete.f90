! $Id: matrix_delete.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrix/delete *
!
! NAME
!    delete - Delete the data contained in the matrix type.
!
! SYNOPSIS
!    call delete(matrix)
! 
! DESCRIPTION
!    This subroutine deletes all data / memory hold in the given matrix
!    type variable
!
! INPUTS
!    matrix
!
! OUTPUT
!
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!
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
! 1. General matrices
!-------------------------------------------------------------------------------

subroutine matrix_delete_ge(matrix)

! 1.1 Declarations
! ----------------

  use matrix_types

  implicit none

  type(matrix_ge), intent(inout) :: matrix

! 1.2 Deallocate memory
! ---------------------

  if (associated(matrix%d)) deallocate(matrix%d)
  if (associated(matrix%e)) deallocate(matrix%e)
  if (associated(matrix%f)) deallocate(matrix%f)
  if (associated(matrix%s)) deallocate(matrix%s)

  if (associated(matrix%i)) deallocate(matrix%i)
  if (associated(matrix%h)) deallocate(matrix%h)
  if (associated(matrix%g)) deallocate(matrix%g)
  if (associated(matrix%r)) deallocate(matrix%r)
  if (associated(matrix%c)) deallocate(matrix%c)

! 1.3 Set default values
! ----------------------

  matrix%fact_chol  = .false.
  matrix%equi_chol = 'N'

  matrix%fact_lu  = .false.
  matrix%equi_lu = 'N'

end subroutine matrix_delete_ge


!-------------------------------------------------------------------------------
! 2. Positive definite matrices
!-------------------------------------------------------------------------------

subroutine matrix_delete_pp(matrix)

! 2.1 Declarations
! ----------------

  use matrix_types

  implicit none

  type(matrix_pp), intent(inout) :: matrix

! 2.2 Deallocate memory
! ---------------------

  if (associated(matrix%d)) deallocate(matrix%d)
  if (associated(matrix%e)) deallocate(matrix%e)
  if (associated(matrix%f)) deallocate(matrix%f)
  if (associated(matrix%s)) deallocate(matrix%s)

! 2.3 Set default values
! ----------------------

  matrix%fact_chol = .false.
  matrix%equi_chol = 'N'

end subroutine matrix_delete_pp
