! $Id: matrix_full2pp.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_full2pp *
!
! NAME
!    matrix_full2pp - Convert a symmetric or hermitian matrix into full
!                     (or general) form into Lapack's packed form.
!
! SYNOPSIS
!    call matrix_full2pp(matrix, packed [, uplo])
! 
! DESCRIPTION
!    This subroutine copies a symmetric matrix into Lapack's packed form.
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
!    To convert a general symmetric matrix into a packed form, use
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

subroutine matrix_full2pp_float(array, packed, uplo)

! 1.1 Declarations
! ----------------

  use typesizes, only: sp => FourByteReal
  use messages
  use ropp_utils

  implicit none

  real(sp), dimension(:,:), intent(in)  :: array
  real(sp), dimension(:),   intent(out) :: packed
  character(len = *),       optional    :: uplo

  integer                               :: i, j, n
  character                             :: ul
  character(len = 64)                   :: size_str, req_str

! 1.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 1.3 Dimensions
! --------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (size(packed) < n*(n+1)/2) then
     write(size_str, '(i10)') size(packed)  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') n*(n+1)/2      ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Size of packed array is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // '.\n')
  endif

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

end subroutine matrix_full2pp_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_full2pp_double(array, packed, uplo)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use ropp_utils

  implicit none

  real(wp), dimension(:,:), intent(in)  :: array
  real(wp), dimension(:),   intent(out) :: packed
  character(len = *),       optional    :: uplo

  integer                               :: i, j, n
  character                             :: ul
  character(len = 64)                   :: size_str, req_str

! 2.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 2.3 Dimensions
! --------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (size(packed) < n*(n+1)/2) then
     write(size_str, '(i10)') size(packed)  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') n*(n+1)/2      ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Size of packed array is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // '.\n')
  endif

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

end subroutine matrix_full2pp_double


!-------------------------------------------------------------------------------
! 3. Double packed matrix type
!-------------------------------------------------------------------------------

subroutine matrix_full2mpp_double(array, packed, uplo)

! 3.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use ropp_utils
  use matrix_types

  implicit none

  real(wp), dimension(:,:), intent(in)    :: array
  type(matrix_pp),          intent(inout) :: packed
  character(len = *),       optional      :: uplo

  integer                                 :: i, j, n
  character                               :: ul
  character(len = 64)                     :: size_str, req_str

! 3.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 3.3 Dimensions
! --------------

  n = size(array, 1)

  if (size(array, 2) /= n) then
     call message(msg_fatal, 'Input matrix must be quadratic.\n')
  endif

  if (size(packed%d) < n*(n+1)/2) then
     write(size_str, '(i10)') size(packed%d)  ;  size_str = adjustl(size_str)
     write(req_str, '(i10)') n*(n+1)/2        ;  req_str  = adjustl(req_str)
     call message(msg_fatal, &
          'Size of packed array is to small: ' // trim(size_str) // &
          ' instead of ' // trim(req_str) // '.\n')
  endif

! 3.4 Copy matrix elements
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

end subroutine matrix_full2mpp_double
