! $Id: matrix_invert.f90 3551 2013-02-25 09:51:28Z idculv $

!****f* Matrix/matrix_invert *
!
! NAME
!    matrix_invert - Invert a square matrix 
!
! SYNOPSIS
!    use matrix
!      ...
!    x = matrix_invert(A)
! 
! DESCRIPTION
!    This function calls matrix_solve to solve a linear equation of the form
!
!       A x = b
!
!    for matrix A. It is assumed that the matrix A is positive 
!    definite, and the solution is obtained using a Cholesky decomposition. 
!    The matrix b is defined zero everywhere except for one element equal to 1.
!
! INPUTS
!    A     Matrix to be inverted.
!
! OUTPUT
!    x     Inverted matrix.
!
! NOTES
!    The matrix A may be in either full or packed form. If the matrix is not 
!    positive definite, the attempt to solve will fail with an error message.
!    These routines are currently available in double precision only.
!
! EXAMPLE
!    To invert a positive definite matrix A,
!
!       A_inverted = matrix_invert(A)
!
! SEE ALSO
!    matrix_types
!    matrix_solve
!
! REFERENCES
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
! 1. Full matrices, double precision
!-------------------------------------------------------------------------------

function matrix_invert_gen(A) result(x)

! 1.1 Declarations

  use typesizes, only: wp => EightByteReal
  use messages
! use matrix, not_this => matrix_invert_gen
  use matrix

  implicit none

  real(wp), dimension(:,:),    intent(inout)   :: A
  real(wp), dimension(size(A,1), size(A,2)) :: x
  
  real(wp), dimension(size(A,1), size(A,2)) :: b
  integer                                   :: i
  character(len = 256)                      :: routine

! 1.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_invert')

! 1.3 Solve linear equation to find elements A^-1
  
  b(:,:) = 0.0_wp
  do i=1,size(b,1)
     b(i,i) = 1.0_wp
  enddo
  
  x =  matrix_solve(A, b)
  
! 1.4 Clean up

  call message_set_routine(routine)

end function matrix_invert_gen


!-------------------------------------------------------------------------------
! 2. Positive definite packed matrix, double precision
!-------------------------------------------------------------------------------

function matrix_invert_packed(A) result(x)
  
! 2.1 Declarations

  use typesizes, only: wp => EightByteReal
  use messages
! use matrix, not_this => matrix_invert_packed
  use matrix
  
  implicit none

  type(matrix_pp),    intent(inout)                          :: A
  real(wp), dimension((int(sqrt(8.*size(A%d))+1)-1)/2, &
                         (int(sqrt(8.*size(A%d))+1)-1)/2) :: x

  real(wp), dimension((int(sqrt(8.*size(A%d))+1)-1)/2, &
                         (int(sqrt(8.*size(A%d))+1)-1)/2) :: b
  integer                                                 :: i
  character(len = 256)                                    :: routine

! 2.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_invert')

! 2.3 Solve linear equation to find elements A^-1

  b(:,:) = 0.0_wp
  do i=1,size(b,1)
     b(i,i) = 1.0_wp
  enddo
  
  x =  matrix_solve(A, b)
  
! 2.4 Clean up

  call message_set_routine(routine)

end function matrix_invert_packed
  
  
