! $Id: matrix.f90 3551 2013-02-25 09:51:28Z idculv $

module matrix

!****m* Modules/matrix *
!
! NAME
!    matrix - Matrix routines and functions.
!
! SYNOPSIS
!    use matrix
!
! DESCRIPTION
!    This Fortran module provides interfaces and data types required for
!    matrix operations.
!
! NOTES
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
! 1. Include derived data types
!-------------------------------------------------------------------------------

  use matrix_types

!-------------------------------------------------------------------------------
! 2. Public interfaces
!-------------------------------------------------------------------------------

! 2.1 Assignment
! --------------

  interface assignment(=)
     subroutine matrix_full2get(ge, data)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       type(matrix_ge),          intent(inout) :: ge
       real(wp), dimension(:,:), intent(in)    :: data
     end subroutine matrix_full2get
     subroutine matrix_full2ppt(pp, data)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),          intent(inout) :: pp
       real(wp), dimension(:,:), intent(in)    :: data
     end subroutine matrix_full2ppt
     subroutine matrix_pp2get(ge, data)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       type(matrix_ge),        intent(inout) :: ge
       real(wp), dimension(:), intent(in)    :: data
     end subroutine matrix_pp2get
     subroutine matrix_pp2ppt(pp, data)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),        intent(inout) :: pp
       real(wp), dimension(:), intent(in)    :: data
     end subroutine matrix_pp2ppt
     subroutine matrix_get2get(oge, ige)
       use matrix_types
       type(matrix_ge),        intent(inout) :: oge
       type(matrix_ge),        intent(in)    :: ige
     end subroutine matrix_get2get
     subroutine matrix_get2ppt(pp, ge)
       use matrix_types
       type(matrix_pp),        intent(inout) :: pp
       type(matrix_ge),        intent(in)    :: ge
     end subroutine matrix_get2ppt
     subroutine matrix_ppt2get(ge, pp)
       use matrix_types
       type(matrix_ge),        intent(inout) :: ge
       type(matrix_pp),        intent(in)    :: pp
     end subroutine matrix_ppt2get
     subroutine matrix_ppt2ppt(opp, ipp)
       use matrix_types
       type(matrix_pp),        intent(inout) :: opp
       type(matrix_pp),        intent(in)    :: ipp
     end subroutine matrix_ppt2ppt
     subroutine matrix_get2full(data, ge)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:,:), intent(inout) :: data
       type(matrix_ge),          intent(in)    :: ge
     end subroutine matrix_get2full
     subroutine matrix_ppt2full(data, pp)
       use typesizes,   only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:,:), intent(inout) :: data
       type(matrix_pp),          intent(in)    :: pp
     end subroutine matrix_ppt2full
  end interface

! 2.2 Addition
! ------------

  interface operator(+)
     function matrix_plus_get (A, B) result (C)
       use matrix_types
       implicit none
       type (matrix_ge), intent(in) :: A
       type (matrix_ge), intent(in) :: B
       type (matrix_ge) :: C
     end function matrix_plus_get
     function matrix_plus_ppt (A, B) result (C)
       use matrix_types
       implicit none
       type (matrix_pp), intent(in) :: A
       type (matrix_pp), intent(in) :: B
       type (matrix_pp) :: C
     end function matrix_plus_ppt
  end interface

! 2.3 Subtraction
! ---------------

  interface operator(-)
     function matrix_minus_get (A, B) result (C)
       use matrix_types
       implicit none
       type (matrix_ge), intent(in) :: A
       type (matrix_ge), intent(in) :: B
       type (matrix_ge) :: C
     end function matrix_minus_get
     function matrix_minus_ppt (A, B) result (C)
       use matrix_types
       implicit none
       type (matrix_pp), intent(in) :: A
       type (matrix_pp), intent(in) :: B
       type (matrix_pp) :: C
     end function matrix_minus_ppt
  end interface

! 2.4 Scalar multiplication
! -------------------------

  interface operator(*)
     function matrix_smul_get_float (s, A) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       real (wp), intent(in) :: s
       type (matrix_ge), intent(in) :: A
       type (matrix_ge) :: C
     end function matrix_smul_get_float
     function matrix_smul_get_dble (s, A) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       real (wp), intent(in) :: s
       type (matrix_ge), intent(in) :: A
       type (matrix_ge) :: C
     end function matrix_smul_get_dble
     function matrix_muls_get_float (A, s) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       type (matrix_ge), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_ge) :: C
     end function matrix_muls_get_float
     function matrix_muls_get_dble (A, s) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type (matrix_ge), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_ge) :: C
     end function matrix_muls_get_dble
     function matrix_smul_ppt_float (s, A) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_pp), intent(in) :: A
       type (matrix_pp) :: C
     end function matrix_smul_ppt_float
     function matrix_smul_ppt_dble (s, A) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_pp), intent(in) :: A
       type (matrix_pp) :: C
     end function matrix_smul_ppt_dble
     function matrix_muls_ppt_float (A, s) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       type (matrix_pp), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_pp) :: C
     end function matrix_muls_ppt_float
     function matrix_muls_ppt_dble (A, s) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type (matrix_pp), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_pp) :: C
     end function matrix_muls_ppt_dble
  end interface

! 2.5 Scalar division
! -------------------

  interface operator(/)
     function matrix_sdiv_get_float (s, A) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_ge), intent(in) :: A
       type (matrix_ge) :: C
     end function matrix_sdiv_get_float
     function matrix_sdiv_get_dble (s, A) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_ge), intent(in) :: A
       type (matrix_ge) :: C
     end function matrix_sdiv_get_dble
     function matrix_divs_get_float (A, s) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       implicit none
       type (matrix_ge), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_ge) :: C
     end function matrix_divs_get_float
     function matrix_divs_get_dble (A, s) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       implicit none
       type (matrix_ge), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_ge) :: C
     end function matrix_divs_get_dble
     function matrix_sdiv_ppt_float (s, A) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_pp), intent(in) :: A
       type (matrix_pp) :: C
     end function matrix_sdiv_ppt_float
     function matrix_sdiv_ppt_dble (s, A) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       implicit none
       real (wp), intent(in) :: s
       type (matrix_pp), intent(in) :: A
       type (matrix_pp) :: C
     end function matrix_sdiv_ppt_dble
     function matrix_divs_ppt_float (A, s) result (C)
       use typesizes, only: wp => FourByteReal
       use matrix_types
       implicit none
       type (matrix_pp), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_pp) :: C
     end function matrix_divs_ppt_float
     function matrix_divs_ppt_dble (A, s) result (C)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       implicit none
       type (matrix_pp), intent(in) :: A
       real (wp), intent(in) :: s
       type (matrix_pp) :: C
     end function matrix_divs_ppt_dble
  end interface

! 2.6 Destructor
! --------------

 interface delete
    subroutine matrix_delete_ge(matrix)
      use matrix_types
      type(matrix_ge), intent(inout) :: matrix
    end subroutine matrix_delete_ge
    subroutine matrix_delete_pp(matrix)
      use matrix_types
      type(matrix_pp), intent(inout) :: matrix
    end subroutine matrix_delete_pp
 end interface

! 2.7 Conversions
! ---------------

  interface matrix_full2bm
     subroutine matrix_full2bm_float (matrix, ku, banded)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent (in) :: matrix
       integer, intent (in) :: ku
       real (wp), dimension (:, :), intent (out) :: banded
     end subroutine matrix_full2bm_float
     subroutine matrix_full2bm_double (matrix, ku, banded)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent (in) :: matrix
       integer, intent (in) :: ku
       real (wp), dimension (:, :), intent (out) :: banded
     end subroutine matrix_full2bm_double
  end interface

  interface matrix_full2bm_alloc
     subroutine matrix_full2bm_alloc_float (matrix, ku, kl, banded)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent (in) :: matrix
       integer, intent (in) :: ku
       integer, intent (in) :: kl
       real (wp), dimension (:, :), pointer :: banded
     end subroutine matrix_full2bm_alloc_float
     subroutine matrix_full2bm_alloc_double (matrix, ku, kl, banded)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent (in) :: matrix
       integer, intent (in) :: ku
       integer, intent (in) :: kl
       real (wp), dimension (:, :), pointer :: banded
     end subroutine matrix_full2bm_alloc_double
  end interface

  interface matrix_bm2full
     subroutine matrix_bm2full_float (banded, ku, matrix)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent (in) :: banded
       integer, intent (in) :: ku
       real (wp), dimension (:, :), intent (out) :: matrix
     end subroutine matrix_bm2full_float
     subroutine matrix_bm2full_double (banded, ku, matrix)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent (in) :: banded
       integer, intent (in) :: ku
       real (wp), dimension (:, :), intent (out) :: matrix
     end subroutine matrix_bm2full_double
  end interface

  interface matrix_bm2full_alloc
      subroutine matrix_bm2full_alloc_float (banded, m, ku, matrix)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent (in) :: banded
       integer, intent (in) :: m
       integer, intent (in) :: ku
       real (wp), dimension (:, :), pointer :: matrix
     end subroutine matrix_bm2full_alloc_float
     subroutine matrix_bm2full_alloc_double (banded, m, ku, matrix)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent (in) :: banded
       integer, intent (in) :: m
       integer, intent (in) :: ku
       real (wp), dimension (:, :), pointer :: matrix
     end subroutine matrix_bm2full_alloc_double
  end interface

  interface matrix_full2pp
     subroutine matrix_full2pp_float (array, packed, uplo)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent(in)  :: array
       real (wp), dimension (:),    intent(out) :: packed
       character(len=*),            optional    :: uplo
     end subroutine matrix_full2pp_float
     subroutine matrix_full2pp_double (array, packed, uplo)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent(in)  :: array
       real (wp), dimension (:),    intent(out) :: packed
       character(len=*),            optional    :: uplo
     end subroutine matrix_full2pp_double
     subroutine matrix_full2mpp_double(array, packed, uplo)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:,:), intent(in)    :: array
       type(matrix_pp),          intent(inout) :: packed
       character(len = *),       optional      :: uplo
     end subroutine matrix_full2mpp_double
  end interface

  interface matrix_full2pp_alloc
     subroutine matrix_full2pp_alloc_float (array, packed, uplo)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:, :), intent(in) :: array
       real (wp), dimension (:),    pointer    :: packed
       character(len=*),            optional   :: uplo
     end subroutine matrix_full2pp_alloc_float
     subroutine matrix_full2pp_alloc_double (array, packed, uplo)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:, :), intent(in) :: array
       real (wp), dimension (:),    pointer    :: packed
       character(len=*),            optional   :: uplo
     end subroutine matrix_full2pp_alloc_double
     subroutine matrix_full2mpp_alloc_double(array, packed, uplo)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:,:), intent(in)    :: array
       type(matrix_pp),          intent(inout) :: packed
       character(len = *),       optional      :: uplo
     end subroutine matrix_full2mpp_alloc_double
  end interface

  interface matrix_pp2full
     subroutine matrix_pp2full_float (packed, array, uplo)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:),    intent(in)  :: packed
       real (wp), dimension (:, :), intent(out) :: array
       character(len=*),            optional    :: uplo
     end subroutine matrix_pp2full_float
     subroutine matrix_pp2full_double (packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:),    intent(in)  :: packed
       real (wp), dimension (:, :), intent(out) :: array
       character(len=*),            optional    :: uplo
     end subroutine matrix_pp2full_double
     subroutine matrix_mpp2full_double(packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),          intent(in)  :: packed
       real(wp), dimension(:,:), intent(out) :: array
       character(len = *),       optional    :: uplo
     end subroutine matrix_mpp2full_double
  end interface

  interface matrix_pp2full_alloc
     subroutine matrix_pp2full_alloc_float (packed, array, uplo)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:),    intent(in) :: packed
       real (wp), dimension (:, :), pointer    :: array
       character(len=*),            optional   :: uplo
     end subroutine matrix_pp2full_alloc_float
     subroutine matrix_pp2full_alloc_double (packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:),    intent(in) :: packed
       real (wp), dimension (:, :), pointer    :: array
       character(len=*),            optional   :: uplo
     end subroutine matrix_pp2full_alloc_double
     subroutine matrix_mpp2full_alloc_double(packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),          intent(in)  :: packed
       real(wp), dimension(:,:), pointer     :: array
       character(len = *),       optional    :: uplo
     end subroutine matrix_mpp2full_alloc_double
  end interface

  interface matrix_pp2full_subset
     subroutine matrix_pp2full_float_subset (packed, array, uplo)
       use typesizes, only: wp => FourByteReal
       real (wp), dimension (:),    intent(in)  :: packed
       real (wp), dimension (:, :), intent(out) :: array
       character(len=*),            optional    :: uplo
     end subroutine matrix_pp2full_float_subset
     subroutine matrix_pp2full_double_subset (packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       real (wp), dimension (:),    intent(in)  :: packed
       real (wp), dimension (:, :), intent(out) :: array
       character(len=*),            optional    :: uplo
     end subroutine matrix_pp2full_double_subset
     subroutine matrix_mpp2full_double_subset(packed, array, uplo)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),          intent(in)  :: packed
       real(wp), dimension(:,:), intent(out) :: array
       character(len = *),       optional    :: uplo
     end subroutine matrix_mpp2full_double_subset
  end interface

! 2.8 General matrix solution
! ---------------------------

  interface matrix_solve
     function matrix_solve_gen_1d (A, b) result(x)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:,:), intent(inout) :: A
       real(wp), dimension(:),   intent(in)    :: b
       real(wp), dimension(size(b))            :: x
     end function matrix_solve_gen_1d
     function matrix_solve_gen_2d (A, b) result(x)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:,:), intent(inout)  :: A
       real(wp), dimension(:,:),   intent(in)   :: b
       real(wp), dimension(size(b,1),size(b,2)) :: x
     end function matrix_solve_gen_2d
     function matrix_solve_packed_1d (A, b) result(x)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),        intent(inout) :: A
       real(wp), dimension(:), intent(in)    :: b
       real(wp), dimension(size(b))          :: x
     end function matrix_solve_packed_1d
     function matrix_solve_packed_2d (A, b) result(x)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),        intent(inout)    :: A
       real(wp), dimension(:,:), intent(in)     :: b
       real(wp), dimension(size(b,1),size(b,2)) :: x
     end function matrix_solve_packed_2d

  end interface

! 2.8 Positive definite matrix inverse
! ------------------------------------

  interface matrix_invert
     function matrix_invert_gen (A) result(x)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:,:),                 intent(inout) :: A
       real(wp), dimension(size(A,1),size(A,2))                :: x
     end function matrix_invert_gen
     function matrix_invert_packed (A) result(x)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       type(matrix_pp),                        intent(inout)  :: A
       real(wp), dimension((int(sqrt(8.*size(A%d))+1)-1)/2, &
                           (int(sqrt(8.*size(A%d))+1)-1)/2)   :: x
     end function matrix_invert_packed
  end interface


! 2.10 Symmetric square root and its inverse
! -----------------------------------------

  interface matrix_sqrt
     subroutine matrix_sqrt_gen(A, R)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:,:), intent(in)    :: A
       type(matrix_sq),          intent(inout) :: R
     end subroutine matrix_sqrt_gen
     subroutine matrix_sqrt_get(A, R)
       use matrix_types
       type(matrix_ge),          intent(in)    :: A
       type(matrix_sq),          intent(inout) :: R
     end subroutine matrix_sqrt_get
     subroutine matrix_sqrt_pp(A, R)
       use typesizes, only: wp => EightByteReal
       use matrix_types
       real(wp), dimension(:),   intent(in)    :: A
       type(matrix_sq),          intent(inout) :: R
     end subroutine matrix_sqrt_pp
     subroutine matrix_sqrt_ppt(A, R)
       use matrix_types
       type(matrix_pp),          intent(in)    :: A
       type(matrix_sq),          intent(inout) :: R
     end subroutine matrix_sqrt_ppt
  end interface

! 2.11 Toast product
! ------------------

  interface matrix_toast
     subroutine matrix_toast_float(A, B, BABt)
       use typesizes, only: wp => FourByteReal
       real(wp), dimension(:,:), intent(inout) :: A
       real(wp), dimension(:,:), intent(in)    :: B
       real(wp), dimension(:,:), optional      :: BABt
     end subroutine matrix_toast_float
     subroutine matrix_toast_double(A, B, BABt)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:,:), intent(inout) :: A
       real(wp), dimension(:,:), intent(in)    :: B
       real(wp), dimension(:,:), optional      :: BABt
     end subroutine matrix_toast_double
  end interface

! 2.12 Matrix SVD
! ---------------

  interface matrix_svd
     subroutine matrix_svd(A, svalues, U, V)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:,:)  ,               intent(in)  ::  A
       real(wp), dimension(size(A,1)),           intent(out) ::  svalues
       real(wp), dimension(size(A,1),size(A,2)), intent(out) ::  U
       real(wp), dimension(size(A,2),size(A,2)), intent(out) ::  V
     end subroutine matrix_svd
  end interface


end module matrix
