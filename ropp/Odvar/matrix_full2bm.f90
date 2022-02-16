! $Id: matrix_full2bm.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_full2bm *
!
! NAME
!    matrix_full2bm - Convert a symmetric or hermitian matrix in full
!                     (or general) form into Lapack's banded form.
!
! SYNOPSIS
!    call matrix_full2bm(matrix, ku, banded)
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix into Lapack's banded form.
!
! INPUTS
!    matrix    Matrix in full (general) form.
!    ku        Number of upper side diagonals.
!
! OUTPUT
!    banded    Matrix in banded form.
!
! NOTES
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

subroutine matrix_full2bm_float(matrix, ku, banded)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => FourByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: matrix
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), intent(out) :: banded

  integer                               :: i, j, k, m, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 1.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (size(banded, 2) /= n) then
     call message(msg_fatal, &
          'Full and banded matrix must have the same number of columns.\n')
  endif

  if (kl < 0) then
     write(size_str, '(i10)') mm     ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') ku + 1  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of banded array is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' + kl.\n')
  endif

! 1.3 Copy matrix elements
! ------------------------

  banded(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        banded(k + i,j) = matrix(i,j)
     enddo
  enddo

end subroutine matrix_full2bm_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_full2bm_double(matrix, ku, banded)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: matrix
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), intent(out) :: banded

  integer                               :: i, j, k, m, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 2.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (size(banded, 2) /= n) then
     call message(msg_fatal, &
          'Full and banded matrix must have the same number of columns.\n')
  endif

  if (kl < 0) then
     write(size_str, '(i10)') mm     ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') ku + 1  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of banded array is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' + kl.\n')
  endif

! 2.3 Copy matrix elements
! ------------------------

  banded(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        banded(k + i,j) = matrix(i,j)
     enddo
  enddo

end subroutine matrix_full2bm_double
