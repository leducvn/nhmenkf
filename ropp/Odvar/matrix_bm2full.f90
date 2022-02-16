! $Id: matrix_bm2full.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_bm2full *
!
! NAME
!    matrix_bm2full - Convert a symmetric or hermitian matrix from
!                     Lapack's banded into full (general) form.
!
! SYNOPSIS
!    call matrix_bm2full(banded, ku, matrix)
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix in Lapack's packed form
!    into full (or general) form.
!
! INPUTS
!    banded    Matrix in banded form.
!    ku        Number of upper side diagonals.
!
! OUTPUT
!    matrix    Matrix in full (general) form.
!
! NOTES
!    The number of lower side diagonals will be computed from the shape
!    of the full matrix. If the dimensions of the full matrix do not equal
!    the shape of the original matrix, the conversion will give a wrong
!    result.
!
!    See the Lapack User Guide for details on matrix storage systems 
!    supported by Lapack. 
!
! EXAMPLE
!    To convert a general symmetric matrix with ku side diagonals into 
!    banded form, use
!
!       call matrix_full2bm(gen_matrix, ku, banded_matrix)
!
!    A full matrix can be obtained from a banded symmetric matrix via
!
!       call matrix_bm2full(banded_matrix, ku, gen_matrix)
!
! SEE ALSO
!    matrix_pp2full
!    matrix_pp2full_alloc
!
!    matrix_full2pp
!    matrix_full2pp_alloc
!
!    matrix_bm2full
!    matrix_bm2full_alloc
!
!    matrix_full2bm
!    matrix_full2bm_alloc
!
! REFERENCES
!    E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra,
!       J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, D. Sorensen,
!       LAPACK Users' Guide, 3rd Ed., SIAM, 1999.
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
! 1. Float
!-------------------------------------------------------------------------------

subroutine matrix_bm2full_float(banded, ku, matrix)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => FourByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: banded
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), intent(out) :: matrix

  integer                               :: i, j, k, m, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 1.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (size(matrix, 2) /= n) then
     call message(msg_fatal, &
          'Banded and full matrix must have the same number of columns.\n')
  endif

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Requested 1st dimension of full matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' m + ?.\n')
  endif

! 1.3 Copy matrix elements
! ------------------------

  matrix(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        matrix(i,j) = banded(k + i,j)
     enddo
  enddo

end subroutine matrix_bm2full_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_bm2full_double(banded, ku, matrix)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: banded
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), intent(out) :: matrix

  integer                               :: i, j, k, m, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 2.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (size(matrix, 2) /= n) then
     call message(msg_fatal, &
          'Banded and full matrix must have the same number of columns.\n')
  endif

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Requested 1st dimension of full matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' m + ?.\n')
  endif

! 2.3 Copy matrix elements
! ------------------------

  matrix(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        matrix(i,j) = banded(k + i,j)
     enddo
  enddo

end subroutine matrix_bm2full_double
