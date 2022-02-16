! $Id: ropp_fm_iono_unpack_bangle.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_fm_iono_unpack_bangle *
!
! NAME
!    ropp_fm_iono_unpack_bangle - Unpack (bangle_L1, bangle_L2) etc into 
!                                 {bangle_n, bangle_L1, bangle_L2} etc
!
! SYNOPSIS
!    CALL ropp_fm_iono_unpack_bangle(res_data)
!
! DESCRIPTION
!    This subroutine is invoked if the L1 and L2 bending angles are being 
!    modelled directly, via a model Chapman layer ionosphere 
!    (ie if -direct_ion is in force). It repacks the concatenated 
!    bangle vector (bangle_L1, bangle_L2) into bangle_n, bangle_L1 
!    and bangle_L2.  Similarly for impact parameters, sigmas and quals.
!
! INPUTS
!    res_data    ROprof structure
!
! OUTPUT
!    res_data    Modified ROprof structure
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!    ropp_fm_bg2ro_1d
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

SUBROUTINE ropp_fm_iono_unpack_bangle(res_data)

! 1.1 Declarations
! -----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types
  USE ropp_fm_constants, ONLY: f_L1,f_L2

  TYPE(ROprof), INTENT(inout)         :: res_data

  REAL(wp)                            :: f2_L1, f2_L2
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_n,impact_n,sigma_n,qual_n
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_1,impact_1,sigma_1,qual_1
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bangle_2,impact_2,sigma_2,qual_2
  INTEGER                             :: m, m2
  CHARACTER(len=256)                  :: routine

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_iono_unpack_bangle')

! 1.3 Unpack bangle fields
! -------------------------

  f2_L1 = f_L1**2  ;  f2_L2 = f_L2**2

  m2 = res_data%lev1b%Npoints  ;  m = m2 / 2

  ALLOCATE(bangle_1(m), impact_1(m), sigma_1(m), qual_1(m))
  ALLOCATE(bangle_2(m), impact_2(m), sigma_2(m), qual_2(m))
  ALLOCATE(bangle_n(m), impact_n(m), sigma_n(m), qual_n(m))

  bangle_1 = res_data%Lev1b%bangle(1:m)
  impact_1 = res_data%Lev1b%impact(1:m)
  sigma_1  = res_data%Lev1b%bangle_sigma(1:m)  ! Close enough
  qual_1   = res_data%Lev1b%bangle_qual(1:m)

  bangle_2 = res_data%Lev1b%bangle(m+1:m2)
  impact_2 = res_data%Lev1b%impact(m+1:m2)
  sigma_2  = res_data%Lev1b%bangle_sigma(m+1:m2)  ! Close enough
  qual_2   = res_data%Lev1b%bangle_qual(m+1:m2)

  bangle_n = (bangle_1*f2_L1 - bangle_2*f2_L2) / (f2_L1 - f2_L2) ! By construction
  impact_n = impact_1
  sigma_n  = (((sigma_1*f2_L1)**2)/((f2_L1 - f2_L2)**2)) + &
             (((sigma_2*f2_L2)**2)/((f2_L1 - f2_L2)**2))  ! Close enough
  sigma_n  = SQRT(sigma_n)
  qual_n   = qual_1  ! They're _all_ either 0 or 100 anyway, currently.

  CALL ropp_io_shrink(res_data%Lev1b, 1, m, 1)  ! Reduce res_data%Lev1b to size m

  res_data%Lev1b%bangle_L1       = bangle_1
  res_data%Lev1b%impact_L1       = impact_1
  res_data%Lev1b%bangle_L1_sigma = sigma_1
  res_data%Lev1b%bangle_L1_qual  = qual_1

  res_data%Lev1b%bangle_L2       = bangle_2
  res_data%Lev1b%impact_L2       = impact_2
  res_data%Lev1b%bangle_L2_sigma = sigma_2
  res_data%Lev1b%bangle_L2_qual  = qual_2

  res_data%Lev1b%bangle          = bangle_n
  res_data%Lev1b%impact          = impact_n
  res_data%Lev1b%bangle_sigma    = sigma_n
  res_data%Lev1b%bangle_qual     = qual_n

  DEALLOCATE(bangle_1, impact_1, sigma_1, qual_1)
  DEALLOCATE(bangle_2, impact_2, sigma_2, qual_2)
  DEALLOCATE(bangle_n, impact_n, sigma_n, qual_n)

! 1.4 Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_iono_unpack_bangle

