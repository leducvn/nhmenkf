! $Id: matrix_toast.f90 3551 2013-02-25 09:51:28Z idculv $

!****f* Matrix/matrix_toast *
!
! NAME
!    matrix_toast - Toast a matrix with a second one.
!
! SYNOPSIS
!    call matrix_toast(A, B [, BABt])
! 
! DESCRIPTION
!    This routine calculates the 'toast product' of two matrices. It is
!    possible to toast with a non-square matrix B if the optional result
!    matrix is provided in the call to toast_transpose and dimensioned
!    appropriately.
!
! INPUTS
!    A     Matrix to be toasted.
!    B     Toaster matrix.
!
! OUTPUT
!    BABt  Toasted matrix.
!
! NOTES
!    If no optional argument is given, calculations are done in place. The
!    result is then returned in A. This is only possible if the toaster is
!    a square matrix.
!
! EXAMPLE
!
!
! SEE ALSO
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
! 1. Single precision
!-------------------------------------------------------------------------------

subroutine matrix_toast_float(A, B, BABt)

  use typesizes, only: wp => FourByteReal

  implicit none

  real(wp), dimension(:,:), intent(inout) :: A
  real(wp), dimension(:,:), intent(in)    :: B
  real(wp), dimension(:,:), optional      :: BABt

  if (present(BABt)) then
     BABt = matmul(B, matmul(A, transpose(B)))
  else
     A = matmul(B, matmul(A, transpose(B)))
  endif

end subroutine matrix_toast_float


!-------------------------------------------------------------------------------
! 2. Double precision
!-------------------------------------------------------------------------------

subroutine matrix_toast_double(A, B, BABt)

  use typesizes, only: wp => EightByteReal

  implicit none

  real(wp), dimension(:,:), intent(inout) :: A
  real(wp), dimension(:,:), intent(in)    :: B
  real(wp), dimension(:,:), optional      :: BABt

  if (present(BABt)) then
     BABt = matmul(B, matmul(A, transpose(B)))
  else
     A = matmul(B, matmul(A, transpose(B)))
  endif

end subroutine matrix_toast_double
