! $Id: matrix_operators.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Matrix/matrix_operators *
!
! NAME
!    matrix_operators - Implement matrix operators (i.e., the operators +, -, *, /).
!
! SYNOPSIS
!
! 
! DESCRIPTION
!    This file implements various operators for the types of the matrix clas.
!    In particular, these operations are supported:
!
!      - Addition
!      - Subtraction
!      - Scalar multiplication
!      - Scalar division
!
! INPUTS
!
!
! OUTPUT
!
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!
!
! REFERENCES
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
! 1. Addition
!-------------------------------------------------------------------------------

! 1.1 General matrix type
! -----------------------

function matrix_plus_get(A, B) result (C)

! 1.1.1 Declarations
  
  use matrix_types
  use messages

  implicit none

  type(matrix_ge), intent(in) :: A
  type(matrix_ge), intent(in) :: B
  type(matrix_ge)             :: C

  integer                     :: m, n

! 1.1.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  if ((size(B%d,1) == m) .and. (size(B%d,2) == n)) then
     allocate(C%d(m, n))
     C%d = A%d + B%d
  else
     call message(msg_fatal, "Matrix sizes do not conform.")
  endif

end function matrix_plus_get


! 1.2 General matrix type
! -----------------------

function matrix_plus_ppt(A, B) result (C)

! 1.2.1 Declarations
  
  use matrix_types
  use messages

  implicit none

  type(matrix_pp), intent(in) :: A
  type(matrix_pp), intent(in) :: B
  type(matrix_pp)             :: C

  integer                     :: n

! 1.2.2 Calculate

  n = size(A%d)

  if (size(B%d) == n) then
     allocate(C%d(n))
     C%d = A%d + B%d
  else
     call message(msg_fatal, "Matrix sizes do not conform.")
  endif

end function matrix_plus_ppt


!-------------------------------------------------------------------------------
! 2. Subtraction
!-------------------------------------------------------------------------------

! 2.1 General matrix type
! -----------------------

function matrix_minus_get(A, B) result (C)

! 2.1.1 Declarations
  
  use matrix_types
  use messages

  implicit none

  type(matrix_ge), intent(in) :: A
  type(matrix_ge), intent(in) :: B
  type(matrix_ge)             :: C

  integer                     :: m, n

! 2.1.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  if ((size(B%d,1) == m) .and. (size(B%d,2) == n)) then
     allocate(C%d(m, n))
     C%d = A%d - B%d
  else
     call message(msg_fatal, "Matrix sizes do not conform.")
  endif

end function matrix_minus_get


! 2.2 General matrix type
! -----------------------

function matrix_minus_ppt(A, B) result (C)

! 2.2.1 Declarations
  
  use matrix_types
  use messages

  implicit none

  type(matrix_pp), intent(in) :: A
  type(matrix_pp), intent(in) :: B
  type(matrix_pp)             :: C

  integer                     :: n

! 2.2.2 Calculate

  n = size(A%d)

  if (size(B%d) == n) then
     allocate(C%d(n))
     C%d = A%d - B%d
  else
     call message(msg_fatal, "Matrix sizes do not conform.")
  endif

end function matrix_minus_ppt


!-------------------------------------------------------------------------------
! 3. Scalar multiplication
!-------------------------------------------------------------------------------

! 3.1 General matrix type (float)
! -------------------------------

function matrix_smul_get_float(s, A) result (C)

! 3.1.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.1.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = s * A%d

end function matrix_smul_get_float

! 3.2 General matrix type (double)
! --------------------------------

function matrix_smul_get_dble(s, A) result (C)

! 3.2.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.2.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = s * A%d

end function matrix_smul_get_dble

! 3.3 General matrix type (float)
! -------------------------------

function matrix_muls_get_float(A, s) result (C)

! 3.3.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.3.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = s * A%d

end function matrix_muls_get_float

! 3.4 General matrix type (double)
! --------------------------------

function matrix_muls_get_dble(A, s) result (C)

! 3.4.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.4.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = s * A%d

end function matrix_muls_get_dble

! 3.5 General matrix type (float)
! -------------------------------

function matrix_smul_ppt_float(s, A) result (C)

! 3.5.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.5.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = s * A%d

end function matrix_smul_ppt_float

! 3.6 General matrix type (double)
! --------------------------------

function matrix_smul_ppt_dble(s, A) result (C)

! 3.6.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.6.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = s * A%d

end function matrix_smul_ppt_dble

! 3.7 General matrix type (float)
! -------------------------------

function matrix_muls_ppt_float(A, s) result (C)

! 3.7.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.7.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = s * A%d

end function matrix_muls_ppt_float

! 3.8 General matrix type (double)
! --------------------------------

function matrix_muls_ppt_dble(A, s) result (C)

! 3.8.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.8.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = s * A%d

end function matrix_muls_ppt_dble


!-------------------------------------------------------------------------------
! 3. Scalar division
!-------------------------------------------------------------------------------

! 3.1 General matrix type (float)
! -------------------------------

function matrix_sdiv_get_float(s, A) result (C)

! 3.1.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types
! use matrix, not_this => matrix_sdiv_get_float
  use matrix, only: matrix_invert

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C
  type(matrix_ge)             :: A_copy

  integer                     :: m, n

! 3.1.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)
  A_copy = A

  allocate(C%d(m, n))
  C%d = s * matrix_invert(A_copy%d)

end function matrix_sdiv_get_float

! 3.2 General matrix type (double)
! --------------------------------

function matrix_sdiv_get_dble(s, A) result (C)

! 3.2.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_sdiv_get_dble
  use matrix, only: matrix_invert

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C
  type(matrix_ge)             :: A_copy

  integer                     :: m, n

! 3.2.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)
  A_copy = A

  allocate(C%d(m, n))
  C%d = s * matrix_invert(A_copy%d)

end function matrix_sdiv_get_dble

! 3.3 General matrix type (float)
! ------------------------------

function matrix_divs_get_float(A, s) result (C)

! 3.3.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.3.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = A%d / s

end function matrix_divs_get_float

! 3.4 General matrix type (double)
! --------------------------------

function matrix_divs_get_dble(A, s) result (C)

! 3.4.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_ge), intent(in) :: A
  type(matrix_ge)             :: C

  integer                     :: m, n

! 3.4.2 Calculate

  m = size(A%d, 1)
  n = size(A%d, 2)

  allocate(C%d(m, n))
  C%d = A%d / s

end function matrix_divs_get_dble

! 3.5 General matrix type (float)
! -------------------------------

function matrix_sdiv_ppt_float(s, A) result (C)

! 3.5.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types
! use matrix, not_this => matrix_sdiv_ppt_float
  use matrix, only: operator(/), matrix_invert

  implicit none

  real(wp),        intent(in) :: s
! type(matrix_pp), intent(in) :: A      ! Changed at 21 July, 2016
! type(matrix_pp)             :: C      ! Changed at 21 July, 2016
! type(matrix_pp)             :: A_copy ! Changed at 21 July, 2016
  type(matrix_pb), intent(in) :: A
  type(matrix_pb)             :: C
  type(matrix_pb)             :: A_copy

  integer                     :: n

! 3.5.2 Calculate

! n = size(A%d) ! Changed at 21 July, 2016
  n = size(A%d,1)
  A_copy = A

! allocate(C%d(n)) ! Changed at 21 July, 2016
  allocate(C%d(n,n))
! C = s * matrix_invert(A_copy) ! Changed at 21 July, 2016
  C%d = s * matrix_invert(A_copy%d)

end function matrix_sdiv_ppt_float

! 3.6 General matrix type (double)
! --------------------------------

function matrix_sdiv_ppt_dble(s, A) result (C)

! 3.6.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types
! use matrix, not_this => matrix_sdiv_ppt_dble
  use matrix, ONLY: operator(/), matrix_invert

  implicit none

  real(wp),        intent(in) :: s
! type(matrix_pp), intent(in) :: A      ! Changed at 21 July, 2016
! type(matrix_pp)             :: C      ! Changed at 21 July, 2016
! type(matrix_pp)             :: A_copy ! Changed at 21 July, 2016
  type(matrix_pb), intent(in) :: A
  type(matrix_pb)             :: C
  type(matrix_pb)             :: A_copy

  integer                     :: n

! 3.6.2 Calculate

! n = size(A%d) ! Changed at 21 July, 2016
  n = size(A%d,1)
  A_copy = A

! allocate(C%d(n)) ! Changed at 21 July, 2016
  allocate(C%d(n,n))
! C = s * matrix_invert(A_copy) ! Changed at 21 July, 2016
  C%d = s * matrix_invert(A_copy%d)

end function matrix_sdiv_ppt_dble

! 3.7 General matrix type (float)
! -------------------------------

function matrix_divs_ppt_float(A, s) result (C)

! 3.7.1 Declarations
  
  use typesizes, only: wp => FourByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.7.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = A%d / s

end function matrix_divs_ppt_float

! 3.8 General matrix type (double)
! --------------------------------

function matrix_divs_ppt_dble(A, s) result (C)

! 3.8.1 Declarations
  
  use typesizes, only: wp => EightByteReal
  use matrix_types

  implicit none

  real(wp),        intent(in) :: s
  type(matrix_pp), intent(in) :: A
  type(matrix_pp)             :: C

  integer                     :: n

! 3.8.2 Calculate

  n = size(A%d)

  allocate(C%d(n))
  C%d = A%d / s

end function matrix_divs_ppt_dble
