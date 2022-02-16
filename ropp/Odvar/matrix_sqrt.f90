! $Id: matrix_sqrt.f90 4010 2014-01-10 11:07:40Z idculv $

!****f* Matrix/matrix_sqrt *
!
! NAME
!    matrix_sqrt - Square root and it's inverse of a symmetric matrix
!
! SYNOPSIS
!    call matrix_sqrt(A, R)
! 
! DESCRIPTION
!    This subroutine calculates the (symmetric) square root of symmetric
!    matrix from it's SVD along with the inverse square root.
!
! INPUTS
!    A   Matrix to be used for calculating the preconditioner.
!
! OUTPUT
!    R   Symmetric square root and its inverse.
!
! NOTES
!   The routine works for regular (square) matrices as well as for positive
!   definite matrices in packed form. Matrix types matrix_ge and matrix_pp
!   are also supported.
!
!   If square roots have been calculated before, the matrices will silently
!   be deallocated, and all data will be lost.
!
! EXAMPLE
!   use matrix
!     ...
!   type(matrix_sq) :: R
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
! 1. General matrix
!-------------------------------------------------------------------------------

subroutine matrix_sqrt_gen(A, R)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
! use matrix,    not_this => matrix_sqrt_gen
  use matrix

  implicit none

  real(wp), dimension(:,:),   intent(in)    :: A
  type(matrix_sq),            intent(inout) :: R

  real(wp), dimension(size(A,1))            :: svalues
  real(wp), dimension(size(A,1), size(A,2)) :: U, diag, V

  integer                                   :: i, n
  character(len = 256)                      :: routine

! 1.2 Error handling
! ------------------

  call message_get_routine(routine)
  call message_set_routine('matrix_sqrt (general matrix)')

! 1.3 Prepare arrays
! ------------------

  n = size(A, 1)

  if (size(A, 2) /= n) then
     call message(msg_fatal, 'Matrix used for preconditionig must be quadratic.\n')
  endif

  if (associated(R % L))     deallocate(R % L)
  if (associated(R % L_inv)) deallocate(R % L_inv)

! 1.4 Set up the square root matrices
! -----------------------------------

  diag(:,:) = A  
  allocate(R % L(n,n), R % L_inv(n,n))

! 1.5 Compute singular value decomposition
! ----------------------------------------
  call matrix_svd(diag, svalues, U, V)

  diag(:,:) = 0.0_wp

  do i = 1, n  
     diag(i,i) = sqrt(svalues(i))  
  enddo
  R % L(:,:) = matmul(V, matmul(diag, transpose(V)))

  do i = 1, n  
     diag(i,i) = 1.0_wp / sqrt(svalues(i))  
  enddo
  R % L_inv(:,:) = matmul(V, matmul(diag, transpose(V)))

! 1.6 Clean up
! ------------

  call message_set_routine(routine)

end subroutine matrix_sqrt_gen


!-------------------------------------------------------------------------------
! 2. General matrix type
!-------------------------------------------------------------------------------

subroutine matrix_sqrt_get(A, R)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
! use matrix,    not_this => matrix_sqrt_get
  use matrix

  implicit none

  type(matrix_ge),                intent(in)    :: A
  type(matrix_sq),                intent(inout) :: R

  real(wp), dimension(size(A%d,1))              :: svalues
  real(wp), dimension(size(A%d,1), size(A%d,2)) :: U, diag, V

  integer                                       :: i, n
  character(len = 256)                          :: routine

! 2.2 Error handling
! ------------------

  call message_get_routine(routine)
  call message_set_routine('matrix_sqrt (general matrix)')

! 2.3 Prepare arrays
! ------------------

  n = size(A % d, 1)

  if (size(A % d, 2) /= n) then
     call message(msg_fatal, 'Matrix used for preconditionig must be quadratic.\n')
  endif

  if (associated(R % L))     deallocate(R % L)
  if (associated(R % L_inv)) deallocate(R % L_inv)

! 2.4 Set up the square root matrices
! -----------------------------------

  diag(:,:) = A  
  allocate(R % L(n,n), R % L_inv(n,n))

! 2.5 Compute singular value decomposition
! ----------------------------------------
  
  call matrix_svd(diag, svalues, U, V)

  diag(:,:) = 0.0_wp

  do i = 1, n  
     diag(i,i) = sqrt(svalues(i))  
  enddo
  R % L(:,:) = matmul(V, matmul(diag, transpose(V)))

  do i = 1, n  
     diag(i,i) = 1.0_wp / sqrt(svalues(i))  
  enddo
  R % L_inv(:,:) = matmul(V, matmul(diag, transpose(V)))

! 2.5 Clean up
! ------------

  call message_set_routine(routine)

end subroutine matrix_sqrt_get


!-------------------------------------------------------------------------------
! 3. Packed positive definite matrix
!-------------------------------------------------------------------------------

subroutine matrix_sqrt_pp(A, R)

! 3.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
! use matrix,    not_this => matrix_sqrt_pp
  use matrix

  implicit none

  real(wp), dimension(:),                intent(in)    :: A
  type(matrix_sq),                       intent(inout) :: R

  real(wp), dimension((int(sqrt(8.*size(A))+1.)-1)/2) :: svalues

  real(wp), dimension((int(sqrt(8.*size(A))+1)-1)/2, &
                      (int(sqrt(8.*size(A))+1)-1)/2) :: U, diag, V

  integer                                              :: i, n
  character(len = 256)                                 :: routine

! 3.2 Error handling
! ------------------

  call message_get_routine(routine)
  call message_set_routine('matrix_sqrt (packed matrix)')

! 3.3 Prepare arrays
! ------------------

  n = NINT(((sqrt(8.*size(A))+1.)-1.)/2.)

  if (associated(R % L))     deallocate(R % L)
  if (associated(R % L_inv)) deallocate(R % L_inv)

! 3.4 Set up the square root matrices
! -----------------------------------

  call matrix_pp2full(A, diag)  
  allocate(R % L(n,n), R % L_inv(n,n))

! 3.5 Compute singular value decomposition
! ----------------------------------------
  
  call matrix_svd(diag, svalues, U, V)

  diag(:,:) = 0.0_wp

  do i = 1, n  
     diag(i,i) = sqrt(svalues(i))  
  enddo
  R % L(:,:) = matmul(V, matmul(diag, transpose(V)))

  do i = 1, n  
     diag(i,i) = 1.0_wp / sqrt(svalues(i))  
  enddo
  R % L_inv(:,:) = matmul(V, matmul(diag, transpose(V)))

end subroutine matrix_sqrt_pp


!-------------------------------------------------------------------------------
! 4. Packed positive definite matrix
!-------------------------------------------------------------------------------

subroutine matrix_sqrt_ppt(A, R)

! 4.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
! use matrix,    not_this => matrix_sqrt_ppt
  use matrix

  implicit none

  type(matrix_pp),                         intent(in)    :: A
  type(matrix_sq),                         intent(inout) :: R

  real(wp), dimension((int(sqrt(8.*size(A%d))+1)-1)/2) :: svalues

  real(wp), dimension((int(sqrt(8.*size(A%d))+1)-1)/2, &
                      (int(sqrt(8.*size(A%d))+1)-1)/2) :: U, diag, V

  integer                                                :: i, n

! 4.2 Prepare arrays
! ------------------

  n = NINT((int(sqrt(8.*size(A%d))+1)-1)/2.)
  if (associated(R % L))     deallocate(R % L)
  if (associated(R % L_inv)) deallocate(R % L_inv)

! 4.3 Set up the square root matrices
! -----------------------------------

  call matrix_pp2full(A, diag)  
  allocate(R % L(n,n), R % L_inv(n,n))


! 4.4 Compute singular value decomposition
! ----------------------------------------
  
  call matrix_svd(diag, svalues, U, V)

  diag(:,:) = 0.0_wp

  do i = 1, n  
     diag(i,i) = sqrt(svalues(i))  
  enddo
  R % L(:,:) = matmul(V, matmul(diag, transpose(V)))

  do i = 1, n  
     diag(i,i) = 1.0_wp / sqrt(svalues(i))  
  enddo
  R % L_inv(:,:) = matmul(V, matmul(diag, transpose(V)))

end subroutine matrix_sqrt_ppt
