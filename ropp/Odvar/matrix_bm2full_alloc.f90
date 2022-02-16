! $Id: matrix_bm2full_alloc.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_bm2full_alloc *
!
! NAME
!    matrix_bm2full_alloc - Convert a symmetric or hermitian matrix from 
!                           Lapack's packed into full (general) form.
!
! SYNOPSIS
!    call matrix_bm2full_alloc(banded, m, ku, matrix)
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix in Lapack's banded form 
!    into full (or general) form. The general matrix to be filled is allocated.
!
! INPUTS
!    banded    Matrix in banded form.
!    m         Number rows of the full matrix.
!    ku        Number of upper side diagonals.
!
! OUTPUT
!    matrix    Pointer to a matrix in full (general) form.
!
! NOTES
!    The number of lower side diagonals will be computed from the number m
!    of rows of the full matrix. If the number of rows does not equal the
!    shape of the original matrix, the conversion will give a wrong result.
!
!    The matrix argument is a 2d Fortran pointer which is allocated within
!    this routine. If the pointer has been used before, all data will be lost.
!
!    See the Lapack User Guide for details on matrix storage systems supported 
!    by Lapack. 
!
! EXAMPLE
!    To convert a general symmetric matrix into packed form, use
!
!       call matrix_full2pp(gen_matrix, packed_matrix)
!
!    A full matrix can be obtained from a packed positive definitive matrix via
!
!       call matrix_pp2full(packed_matrix, gen_matrix)
!
! SEE ALSO
!    matrix_pp2full
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

subroutine matrix_bm2full_alloc_float(banded, m, ku, matrix)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => FourByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: banded
  integer,                  intent(in)  :: m
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), pointer     :: matrix

  integer                               :: i, j, k, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 1.2 Dimensions
! --------------

  n = size(banded, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (kl < 0) then
     write(size_str, '(i10)') mm     ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') ku + 1  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of banded matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' + kl.\n')
  endif

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Requested 1st dimension of full matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' m + ?.\n')
  endif

  if (associated(matrix)) then
     call message(msg_warn, &
          'Pointer to full matrix already associated - deallocating.\n')
     deallocate(matrix)
  endif

  allocate(matrix(m, n))

! 1.3 Copy matrix elements
! ------------------------

  matrix(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        matrix(i,j) = banded(k + i,j)
     enddo
  enddo

end subroutine matrix_bm2full_alloc_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_bm2full_alloc_double(banded, m, ku, matrix)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: banded
  integer,                  intent(in)  :: m
  integer,                  intent(in)  :: ku
  real(wp), dimension(:,:), pointer     :: matrix

  integer                               :: i, j, k, mm, kl, n
  character(len = 64)                   :: size_str, req_str

! 2.2 Dimensions
! --------------

  n = size(banded, 2)

  mm = size(banded, 1)
  kl = mm - ku - 1

  if (kl < 0) then
     write(size_str, '(i10)') mm    ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') ku + 1 ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of banded matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' + kl.\n')
  endif

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Requested 1st dimension of full matrix is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' m + ?.\n')
  endif

  if (associated(matrix)) then
     call message(msg_warn, &
          'Pointer to full matrix already associated - deallocating.\n')
     deallocate(matrix)
  endif

  allocate(matrix(m, n))

! 2.3 Copy matrix elements
! ------------------------

  matrix(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        matrix(i,j) = banded(k + i,j)
     enddo
  enddo

end subroutine matrix_bm2full_alloc_double
