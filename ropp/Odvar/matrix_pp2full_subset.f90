! $Id: matrix_pp2full_subset.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrices/matrix_pp2full_subset *
!
! NAME
!    matrix_pp2full_subset - Convert a subset of a symmetric or hermitian 
!                            matrix from Lapack's packed into full (general) 
!                            form.
!
! SYNOPSIS
!    call matrix_pp2full_subset(packed, matrix [, uplo])
! 
! DESCRIPTION
!    This subroutine copies part of a matrix in Lapack's packed form into full 
!    (or general) form.
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

subroutine matrix_pp2full_float_subset(packed, array, uplo)

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

  n = size(array, 1)

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

end subroutine matrix_pp2full_float_subset


!-------------------------------------------------------------------------------
! 2. Double
!-------------------------------------------------------------------------------

subroutine matrix_pp2full_double_subset(packed, array, uplo)

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

  n = size(array,1)

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

end subroutine matrix_pp2full_double_subset


!-------------------------------------------------------------------------------
! 3. Double packed matrix type
!-------------------------------------------------------------------------------

subroutine matrix_mpp2full_double_subset(packed, array, uplo)

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

  n = size(array, 1)

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

end subroutine matrix_mpp2full_double_subset
