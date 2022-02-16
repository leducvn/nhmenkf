! $Id: ropp_fm_iono_bangle.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_fm_iono_bangle *
!
! NAME
!   ropp_fm_iono_bangle - Forward model bending caused by ionosphere.
!
! SYNOPSIS
!   CALL ropp_fm_iono_bangle(Ne_max, R_peak, H_width, n_L1, impact, bangle)
!
! DESCRIPTION
!   Forward model the ionospheric bending for the L1 and L2 bending angles
!   assuming Chapman layer profiles. Should only be called if
!   x%direct_ion=.TRUE.
!
! INPUTS
!   REAL(wp)                :: Ne_max     ! Peak electron density (m-3)
!   REAL(wp)                :: R_peak     ! Radius of peak value (m)
!   REAL(wp)                :: H_width    ! Width of Chapman Layer (m)
!   INTEGER                 :: n_L1       ! No. of L1 bending angles in array
!   REAL(wp), DIMENSION(:)  :: impact     ! Impact parameters
!
! INOUT
!   REAL(wp), DIMENSION(:)  :: bangle     ! Updated bending angles inc. ionospheric bending
!
! NOTES
!
! SEE ALSO
!
! AUTHOR
!   Met Office, Exeter, UK and ECMWF, Reading UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_fm_iono_bangle(Ne_max, R_peak, H_width, n_L1, impact, bangle)

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm,   not_this => ropp_fm_iono_bangle
! USE ropp_fm
  USE ropp_utils, ONLY: ropp_MDTV
! USE ropp_fm_constants, ONLY: k4, f_L1, f_L2

  IMPLICIT NONE

  REAL(wp), INTENT(in)                  :: Ne_max  ! Peak electron density (m-3)
  REAL(wp), INTENT(in)                  :: R_peak  ! Radius of peak value (m)
  REAL(wp), INTENT(in)                  :: H_width ! Width of Chapman Layer (m)
  INTEGER                               :: n_L1    ! Number of L1 bending angles
  REAL(wp), DIMENSION(:), INTENT(in)    :: impact  ! Impact parameters
  REAL(wp), DIMENSION(:), INTENT(inout) :: bangle  ! Updated bending angles

  INTEGER                               :: n_chap, n_impact, i
  REAL(wp)                              :: const_L1,const_L2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: lg, z, q, bangle_ion

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_chap = 1  ! For the moment
  n_impact = SIZE(impact)

  ALLOCATE(lg(n_impact), z(n_impact), q(n_impact), bangle_ion(n_impact))

  bangle_ion = 0.0_wp

!-------------------------------------------------------------------------------
! 3. Calculate the zorro function
!-------------------------------------------------------------------------------

  lg = (R_peak - impact) / H_width

  z = ropp_fm_zorro(lg)

!--------------------------------------------------------------------------------
! 4. Ionospheric bending BUT NOT YET including scaling by frequency
!--------------------------------------------------------------------------------  

  q = 2.0_wp * EXP(0.5_wp) * impact * R_peak / &
      SQRT(H_width*(R_peak + impact)**3)

  bangle_ion = bangle_ion + Ne_max * q * z

!--------------------------------------------------------------------------------
! 5. Now calculate total bending = (neutral + ionospheric)
!--------------------------------------------------------------------------------  

  const_L1 = k4 / f_L1**2
  const_L2 = k4 / f_L2**2

  n_L1 = n_impact / 2

  DO i = 1, n_impact

    IF (bangle(i) < ropp_MDTV) CYCLE ! neutral is missing

    IF (i <= n_L1) THEN ! the L1 bending angles

      bangle(i) = bangle(i) + const_L1 * bangle_ion(i)

    ELSE ! the L2 bending angles

      bangle(i) = bangle(i) + const_L2 * bangle_ion(i)

    END IF

  END DO

!--------------------------------------------------------------------------------
! 6. Deallocate/tidy up
!---------------------------------------------------------------------------------

  DEALLOCATE(lg, z, q, bangle_ion)

  RETURN

END SUBROUTINE ropp_fm_iono_bangle
