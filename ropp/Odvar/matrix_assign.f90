! $Id: matrix_assign.f90 3551 2013-02-25 09:51:28Z idculv $

!****m* Matrix/matrix_assign *
!
! NAME
!    matrix_assign - Implement matrix assignment (i.e., the operator =).
!
! SYNOPSIS
!
! 
! DESCRIPTION
!    This file contains implementations for the assignment operator (=) for
!    the various types in the matrix class.
!
! INPUTS
!
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
! 1. Full to general matrix type
!-------------------------------------------------------------------------------

subroutine matrix_full2get(ge, data)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_full2get

  implicit none

  type(matrix_ge),          intent(inout) :: ge
  real(wp), dimension(:,:), intent(in)    :: data

  integer                                 :: m, n

  call delete(ge)

  m = size(data, 1)
  n = size(data, 2)

  allocate(ge%d(m,n))

  ge%d(:,:) = data(:,:)

end subroutine matrix_full2get


!-------------------------------------------------------------------------------
! 2. Full to pp matrix type
!-------------------------------------------------------------------------------

subroutine matrix_full2ppt(pp, data)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_full2ppt

  implicit none

  type(matrix_pp),          intent(inout) :: pp
  real(wp), dimension(:,:), intent(in)    :: data

  integer                                 :: n

  call delete(pp)

  n = size(data, 1)
  if (n /= size(data, 2)) then
     print *, 'Only quadratic matrices are allowed for conversion into matrix_pp.'
     stop
  endif

  call matrix_full2pp_alloc(data, pp % d)

end subroutine matrix_full2ppt


!-------------------------------------------------------------------------------
! 3. pp to general matrix type
!-------------------------------------------------------------------------------

subroutine matrix_pp2get(ge, data)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_pp2get

  implicit none

  type(matrix_ge),        intent(inout) :: ge
  real(wp), dimension(:), intent(in)    :: data

  call delete(ge)

  call matrix_pp2full_alloc(data, ge % d)

end subroutine matrix_pp2get


!-------------------------------------------------------------------------------
! 4. pp to pp matrix type
!-------------------------------------------------------------------------------

subroutine matrix_pp2ppt(pp, data)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_pp2ppt

  implicit none

  type(matrix_pp),        intent(inout) :: pp
  real(wp), dimension(:), intent(in)    :: data

  call delete(pp)

  allocate(pp%d(size(data)))

  pp % d = data

end subroutine matrix_pp2ppt


!-------------------------------------------------------------------------------
! 5. General type to general type
!-------------------------------------------------------------------------------

subroutine matrix_get2get(oge, ige)

  use matrix_types
! use matrix, not_this => matrix_get2get

  implicit none

  type(matrix_ge),        intent(inout) :: oge
  type(matrix_ge),        intent(in)    :: ige

  call delete(oge)

  allocate(oge%d(size(ige%d,1), size(ige%d,2)))

  oge % d = ige % d

end subroutine matrix_get2get


!-------------------------------------------------------------------------------
! 6. General type to pp type
!-------------------------------------------------------------------------------

subroutine matrix_get2ppt(pp, ge)

  use matrix_types
! use matrix, not_this => matrix_get2ppt

  implicit none

  type(matrix_pp),        intent(inout) :: pp
  type(matrix_ge),        intent(in)    :: ge

  call delete(pp)

  call matrix_full2pp_alloc(ge % d, pp % d)

end subroutine matrix_get2ppt


!-------------------------------------------------------------------------------
! 7. pp type to general type
!-------------------------------------------------------------------------------

subroutine matrix_ppt2get(ge, pp)

  use matrix_types
! use matrix, not_this =>  matrix_ppt2get

  implicit none

  type(matrix_ge),        intent(inout) :: ge
  type(matrix_pp),        intent(in)    :: pp

  call delete(ge)

  call matrix_pp2full_alloc(pp % d, ge % d)

end subroutine matrix_ppt2get


!-------------------------------------------------------------------------------
! 8. pp type to pp type
!-------------------------------------------------------------------------------

subroutine matrix_ppt2ppt(opp, ipp)

  use matrix_types
  use arrays
! use matrix, not_this => matrix_ppt2ppt

  implicit none

  type(matrix_pp),        intent(inout) :: opp
  type(matrix_pp),        intent(in)    :: ipp

  call delete(opp)

  call copy_alloc(ipp % d, opp % d)

  if (size(ipp % e) == 0) then
     if (associated(opp % e)) deallocate(opp % e)
  else
     call copy_alloc(ipp % e, opp % e)
  endif

  if (size(ipp % f) == 0) then
     if (associated(opp % f)) deallocate(opp % f)
  else
     call copy_alloc(ipp % f, opp % f)
  endif

  if (size(ipp % s) == 0) then
     if (associated(opp % s)) deallocate(opp % s)
  else
     call copy_alloc(ipp % s, opp % s)
  endif

  opp % fact_chol = ipp % fact_chol
  opp % equi_chol = ipp % equi_chol

end subroutine matrix_ppt2ppt


!-------------------------------------------------------------------------------
! 9. General type to full matrix
!-------------------------------------------------------------------------------

subroutine matrix_get2full(data, ge)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_get2full
  use messages

  implicit none

  real(wp), dimension(:,:), intent(inout) :: data
  type(matrix_ge),          intent(in)    :: ge

  if ((size(data,1) == size(ge%d,1)) .and. &
       (size(data,1) == size(ge%d,1))) then
     data(:,:) = ge%d(:,:)
  else
     call message(msg_fatal, &
          "Matrix dimensions don't match.")
  endif

end subroutine matrix_get2full


!-------------------------------------------------------------------------------
! 10. pp type to full matrix
!-------------------------------------------------------------------------------

subroutine matrix_ppt2full(data, pp)

  use typesizes,   only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_ppt2full

  implicit none

  real(wp), dimension(:,:), intent(inout) :: data
  type(matrix_pp),          intent(in)    :: pp

  call matrix_pp2full(pp%d, data)

end subroutine matrix_ppt2full
