! $Id: ropp_pp_msis.f90 3551 2013-02-25 09:51:28Z idculv $

!****m* Modules/ropp_pp_MSIS *
!
! NAME
!    ropp_pp_MSIS - Interface module for routines to obtain MSIS bending angle
!                   profiles in the ROPP pre-processor
!
! SYNOPSIS
!    use ropp_pp_MSIS
!
! DESCRIPTION
!    This module provides interfaces for all MSIS bending angle routines in the
!    ROPP Preprocessor library.
!
! SEE ALSO
!    ropp_pp
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

MODULE ropp_pp_MSIS_types

  USE typesizes, ONLY: wp => EightByteReal

!-------------------------------------------------------------------------------
! 1. MSIS spherical harmonic coefficients
!-------------------------------------------------------------------------------

  TYPE MSIScoeff_ref
     REAL(wp), DIMENSION(:,:), POINTER :: Ac0 => null()
     REAL(wp), DIMENSION(:,:), POINTER :: Ac1 => null()
     REAL(wp), DIMENSION(:,:), POINTER :: Bc1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Aa0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Aa1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ba1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ab0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ab1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Bb1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ad0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ad1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Bd1 => null()
  END TYPE MSIScoeff_ref

  TYPE MSIScoeff_ba
     REAL(wp), DIMENSION(:,:), POINTER :: Ac0 => null()
     REAL(wp), DIMENSION(:,:), POINTER :: Ac1 => null()
     REAL(wp), DIMENSION(:,:), POINTER :: Bc1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Aa0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Aa1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ba1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ab0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ab1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Bb1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ad0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ad1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Bd1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ae0 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Ae1 => null()
     REAL(wp), DIMENSION(:),   POINTER :: Be1 => null()
  END TYPE MSIScoeff_ba

  LOGICAL, SAVE, PUBLIC :: MSIS_read = .false.

END MODULE ropp_pp_MSIS_types

MODULE ropp_pp_MSIS

  USE ropp_pp_MSIS_types

!-------------------------------------------------------------------------------
! 2. Read MSIS file
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_read_MSIS
     SUBROUTINE ropp_pp_read_MSIS_ref(file, month, coeff)
       USE ropp_pp_MSIS_types
       CHARACTER(len=*),     INTENT(in)  :: file
       INTEGER,              INTENT(in)  :: month
       TYPE(MSIScoeff_ref),  INTENT(out) :: coeff
     END SUBROUTINE ropp_pp_read_MSIS_ref
          SUBROUTINE ropp_pp_read_MSIS_ba(file, month, coeff)
       USE ropp_pp_MSIS_types
       CHARACTER(len=*),     INTENT(in)  :: file
       INTEGER,              INTENT(in)  :: month
       TYPE(MSIScoeff_ba),   INTENT(out) :: coeff
     END SUBROUTINE ropp_pp_read_MSIS_ba
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_refrac_BG(bfile, month, lat, lon, alt, refrac)
       USE typesizes, ONLY: wp => EightByteReal
       CHARACTER(len=*),       INTENT(in)  :: bfile
       INTEGER,                INTENT(in)  :: month
       REAL(wp),               INTENT(in)  :: lat
       REAL(wp),               INTENT(in)  :: lon
       REAL(wp), DIMENSION(:), INTENT(in)  :: alt
       REAL(wp), DIMENSION(:), INTENT(out) :: refrac
     END SUBROUTINE ropp_pp_refrac_BG
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_refrac_MSIS(mfile, month, lat, lon, alt, refrac, grad)
       USE typesizes, ONLY: wp => EightByteReal
       CHARACTER(len=*),    INTENT(inout)  :: mfile
       INTEGER,                INTENT(in)  :: month
       REAL(wp),               INTENT(in)  :: lat
       REAL(wp),               INTENT(in)  :: lon
       REAL(wp), DIMENSION(:), INTENT(in)  :: alt
       REAL(wp), DIMENSION(:), INTENT(out) :: refrac
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: grad
     END SUBROUTINE ropp_pp_refrac_MSIS
  END INTERFACE

    INTERFACE
     SUBROUTINE ropp_pp_bangle_MSIS(mfile, month, lat, lon, alt, bangle)
       USE typesizes, ONLY: wp => EightByteReal
       CHARACTER(len=*),    INTENT(inout)  :: mfile
       INTEGER,                INTENT(in)  :: month
       REAL(wp),               INTENT(in)  :: lat
       REAL(wp),               INTENT(in)  :: lon
       REAL(wp), DIMENSION(:), INTENT(in)  :: alt
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle
     END SUBROUTINE ropp_pp_bangle_MSIS
  END INTERFACE

END MODULE ropp_pp_MSIS
