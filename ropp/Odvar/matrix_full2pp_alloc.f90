! $Id: matrix_full2pp_alloc.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_full2pp_alloc *
!
! NAME
!    matrix_full2pp_alloc - Convert a symmetric or hermitian matrix into full
!                     (or general) form into Lapack's packed form.
!
! SYNOPSIS
!    call matrix_full2pp_alloc(matrix, packed [, uplo])
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix into Lapack's packed form.
!    The packed matrix to be filled is allocated.
!
! INPUTS
!    matrix    Matrix in full (general) form.
!    uplo      'UPLO' parameter; determines if the packed matrix has been packed
!                 from an upper ('U') or lower ('L') full matrix. This parameter
!                 is optional; it defaults to 'U'.
!
! OUTPUT
!    packed    Matrix in packed form.
!
! NOTES
!    See the Lapack User Guide for details on matrix storage systems
!    supported by Lapack. 
!
! EXAMPLE
!    To convert a general symmetric matrix with ku side diagonals into 
!    banded form, use
!
!       call matrix_full2pp(gen_matrix, packed_matrix)
!
!    A full matrix can be obtained from a packed symmetric matrix via
!
!       call matrix_pp2full(packed_matrix, gen_matrix)
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

subroutine matrix_full2pp_alloc_float(array, packed, uplo)

! 1.1 Declarations
! ----------------

  use typesizes, only: sp => FourByteReal
  use messages
  use ropp_utils

  implicit none

  real(sp), dimension(:,:), intent(in)  :: array
  real(sp), dimension(:),   pointer     :: packed
  character(len = *),       optional    :: uplo

  integer                               :: i, j, n
  character                             :: ul

! 1.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 1.3 Allocate packed matrix
! --------------------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (associated(packed)) then
     call message(msg_warn, &
          'Pointer to packed matrix already associated - deallocating.\n')
     deallocate(packed)
  endif

  allocate(packed(n*(n+1)/2))

! 1.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           packed(i + j*(j-1)/2) = array(i, j)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           packed(i + (2*n-j)*(j-1)/2) = array(i, j)
        enddo
     enddo
  else
     call message(msg_fatal, "UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_full2pp_alloc_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_full2pp_alloc_double(array, packed, uplo)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use ropp_utils

  implicit none

  real(wp), dimension(:,:), intent(in)  :: array
  real(wp), dimension(:),   pointer     :: packed
  character(len = *),       optional    :: uplo

  integer                               :: i, j, n
  character                             :: ul

! 2.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 1.3 Allocate packed matrix
! --------------------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (associated(packed)) then
     call message(msg_warn, &
          'Pointer to packed matrix already associated - deallocating.\n')
     deallocate(packed)
  endif

  allocate(packed(n*(n+1)/2))

! 2.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           packed(i + j*(j-1)/2) = array(i, j)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           packed(i + (2*n-j)*(j-1)/2) = array(i, j)
        enddo
     enddo
  else
     call message(msg_fatal, "UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_full2pp_alloc_double


!-------------------------------------------------------------------------------
! 2. Double packed matrix type
!-------------------------------------------------------------------------------

subroutine matrix_full2mpp_alloc_double(array, packed, uplo)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use matrix_types
! use matrix, not_this => matrix_full2mpp_alloc_double
  use ropp_utils

  implicit none

  real(wp), dimension(:,:), intent(in)    :: array
  type(matrix_pp),          intent(inout) :: packed
  character(len = *),       optional      :: uplo

  integer                                 :: i, j, n
  character                               :: ul

! 2.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 1.3 Allocate packed matrix
! --------------------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (associated(packed%d)) then
     call message(msg_warn, &
          'Pointer to packed matrix already associated - deallocating.\n')
     call delete(packed)
  endif

! Deallocate the remaining elements silently

  allocate(packed%d(n*(n+1)/2))

! 2.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           packed%d(i + j*(j-1)/2) = array(i, j)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           packed%d(i + (2*n-j)*(j-1)/2) = array(i, j)
        enddo
     enddo
  else
     call message(msg_fatal, "UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_full2mpp_alloc_double
