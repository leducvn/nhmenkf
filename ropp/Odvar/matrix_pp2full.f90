! $Id: matrix_pp2full.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_pp2full *
!
! NAME
!    matrix_pp2full - Convert a symmetric or hermitian matrix from Lapack's
!                        packed into full (general) form.
!
! SYNOPSIS
!    call matrix_pp2full(packed, matrix [, uplo])
! 
! DESCRIPTION
!    This subroutine copies a matrix in Lapack's packed form into full (or 
!    general) form.
!
! INPUTS
!    packed    Matrix in packed form.
!    uplo      'UPLO' parameter; determines if the packed matrix has been packed
!                 from an upper ('U') or lower ('L') full matrix. This parameter
!                 is optional; it defaults to 'U'.
!
! OUTPUT
!    matrix    Matrix in full (general) form.
!
! NOTES
!    Packed matrix storage in Lapack is reserved for symmetric or Hermitian
!    postive definite matrices. See the Lapack User Guide for details on this
!    and other matrix storage systems supported by Lapack. 
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

subroutine matrix_pp2full_float(packed, array, uplo)

! 1.1 Declarations
! ----------------

  use typesizes, only: sp => FourByteReal
  use messages
  use ropp_utils

  implicit none

  real(sp), dimension(:),   intent(in)  :: packed
  real(sp), dimension(:,:), intent(out) :: array
  character(len = *),       optional    :: uplo

  integer                               :: i, j
  integer                               :: n

  character                             :: ul
  character(len = 64)                   :: n_str

! 1.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 1.3 Size of full matrix
! -----------------------

  n = NINT((int(sqrt(8.0_sp * size(packed)) + 1.0_sp) - 1.0_sp) / 2.0_sp)

  if (size(array, 1) < n .or. size(array, 2) < n) then
     write(n_str, '(i10)') n  ;  n_str = adjustl(n_str)
     call message(msg_fatal, &
          'Full matrix must a least have dimension ' // trim(n_str) // &
          ' x ' // trim(n_str) // '.\n')
  endif

! 1.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           array(i, j) = packed(i + j*(j-1)/2)
        enddo
     enddo
     do i = 2, n
        do j = 1, i
           array(i, j) = array(j, i)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           array(i, j) = packed(i + (2*n-j)*(j-1)/2)
        enddo
     enddo
     do j = 2, n
        do i = 1, j
           array(i, j) = array(j, i)
        enddo
     enddo
  else
     call message(msg_fatal, "UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_pp2full_float


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_pp2full_double(packed, array, uplo)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use ropp_utils

  implicit none

  real(wp), dimension(:),   intent(in)  :: packed
  real(wp), dimension(:,:), intent(out) :: array
  character(len = *),       optional    :: uplo

  integer                               :: i, j
  integer                               :: n

  character                             :: ul
  character(len = 64)                   :: n_str

! 2.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 2.3 Size of full matrix
! -----------------------

  n = NINT((int(sqrt(8.0_wp * size(packed)) + 1.0_wp) - 1.0_wp) / 2.0_wp)

  if (size(array, 1) < n .or. size(array, 2) < n) then
     write(n_str, '(i10)') n  ;  n_str = adjustl(n_str)
     call message(msg_fatal, &
          'Full matrix must a least have dimension ' // trim(n_str) // &
          ' x ' // trim(n_str) // '.\n')
  endif

! 2.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           array(i, j) = packed(i + j*(j-1)/2)
        enddo
     enddo
     do i = 2, n
        do j = 1, i
           array(i, j) = array(j, i)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           array(i, j) = packed(i + (2*n-j)*(j-1)/2)
        enddo
     enddo
     do j = 2, n
        do i = 1, j
           array(i, j) = array(j, i)
        enddo
     enddo
  else
     call message(msg_fatal, "UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_pp2full_double


!-------------------------------------------------------------------------------
! 3. Double packed matrix type
!-------------------------------------------------------------------------------

subroutine matrix_mpp2full_double(packed, array, uplo)

! 3.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use messages
  use matrix_types
  use ropp_utils

  implicit none

  type(matrix_pp),          intent(in)  :: packed
  real(wp), dimension(:,:), intent(out) :: array
  character(len = *),       optional    :: uplo

  integer                               :: i, j
  integer                               :: n

  character                             :: ul
  character(len = 64)                   :: n_str

! 3.2 Default parameters
! ----------------------

  if (present(uplo)) then
     ul = adjustl(uplo)
     call To_Lower(ul)
  else
     ul = 'u'
  endif

! 3.3 Size of full matrix
! -----------------------

  n = NINT((int(sqrt(8.0_wp * size(packed%d)) + 1.0_wp) - 1.0_wp) / 2.0_wp)

  if (size(array, 1) < n .or. size(array, 2) < n) then
     write(n_str, '(i10)') n  ;  n_str = adjustl(n_str)
     call message(msg_fatal, &
          'Full matrix must a least have dimension ' // trim(n_str) // &
          ' x ' // trim(n_str) // '.\n')
  endif

! 3.4 Copy matrix elements
! ------------------------

  if (ul == 'u') then
     do j = 1, n
        do i = 1, j
           array(i, j) = packed%d(i + j*(j-1)/2)
        enddo
     enddo
     do i = 2, n
        do j = 1, i
           array(i, j) = array(j, i)
        enddo
     enddo
  else if (ul == 'l') then
     do i = 1, n
        do j = 1, i
           array(i, j) = packed%d(i + (2*n-j)*(j-1)/2)
        enddo
     enddo
     do j = 2, n
        do i = 1, j
           array(i, j) = array(j, i)
        enddo
     enddo
  else
     call message(msg_fatal, "If given, UPLO must be either 'U' or 'L'.\n")
  endif

end subroutine matrix_mpp2full_double
