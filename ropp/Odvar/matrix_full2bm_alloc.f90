! $Id: matrix_full2bm_alloc.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_full2bm_alloc *
!
! NAME
!    matrix_full2bm_alloc - Convert a symmetric or hermitian matrix in full
!                           (or general) form into Lapack's banded form.
!
! SYNOPSIS
!    call matrix_full2bm_alloc(matrix, ku, kl, banded)
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix into Lapack's banded form. 
!    The banded matrix to be filled is allocated.
!
! INPUTS
!    matrix    Matrix in full (general) form.
!    ku        Number of upper side diagonals.
!    kl        Number of lower side diagonals.
!
! OUTPUT
!    banded    Matrix in banded form.
!
! NOTES
!    The banded argument is a 2d Fortran pointer which is allocated within
!    this routine. If the pointer has been used before, all data will be lost.
!
!    See the Lapack User Guide for details on matrix storage systems 
!    supported by Lapack. 
!
! EXAMPLE
!    To convert a general symmetric matrix with ku side diagonals into 
!    banded form, use
!
!       call matrix_full2bm_alloc(gen_matrix, ku, kl, banded_matrix)
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

subroutine matrix_full2bm_alloc_float(matrix, ku, kl, banded)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => FourByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: matrix
  integer,                  intent(in)  :: ku
  integer,                  intent(in)  :: kl
  real(wp), dimension(:,:), pointer     :: banded

  integer                               :: i, j, k, m, mm, n
  character(len = 64)                   :: size_str, req_str

! 1.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = kl + ku + 1

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of full matrix is too small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' (ku + kl + 1).\n')
  endif

  if (associated(banded)) then
     call message(msg_warn, &
          'Pointer to banded matrix already associated - deallocating.\n')
     deallocate(banded)
  endif

  allocate(banded(mm, n))

! 1.3 Copy matrix elements
! ------------------------

  banded(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        banded(k + i,j) = matrix(i,j)
     enddo
  enddo

end subroutine matrix_full2bm_alloc_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_full2bm_alloc_double(matrix, ku, kl, banded)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages

  implicit none

  real(wp), dimension(:,:), intent(in)  :: matrix
  integer,                  intent(in)  :: ku
  integer,                  intent(in)  :: kl
  real(wp), dimension(:,:), pointer     :: banded

  integer                               :: i, j, k, m, mm, n
  character(len = 64)                   :: size_str, req_str

! 1.2 Dimensions
! --------------

  m = size(matrix, 1)
  n = size(matrix, 2)

  mm = kl + ku + 1

  if (m < mm) then
     write(size_str, '(i10)') m  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') mm  ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          '1st dimension of full matrix is too small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // ' (ku + kl + 1).\n')
  endif

  if (associated(banded)) then
     call message(msg_warn, &
          'Pointer to banded matrix already associated - deallocating.\n')
     deallocate(banded)
  endif

  allocate(banded(mm, n))

! 1.3 Copy matrix elements
! ------------------------

  banded(:,:) = 0.0_wp
  do j = 1, n
     k = ku + 1 - j
     do i = max(1,j - ku), min(m, j + kl)
        banded(k + i,j) = matrix(i,j)
     enddo
  enddo

end subroutine matrix_full2bm_alloc_double
