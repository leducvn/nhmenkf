! $Id: ropp_fm_iono_bangle_tl.f90 4010 2014-01-10 11:07:40Z idculv $

!****s* Ionosphere/ropp_fm_iono_bangle_tl *
!
! NAME
!   ropp_fm_iono_bangle_tl - Tangent linear of forward model bending caused by ionosphere.
!
! SYNOPSIS
!   CALL ropp_fm_iono_bangle_tl( &
!                               Ne_max,    R_peak,    H_width, &
!                               Ne_max_tl, R_peak_tl, H_width_tl, &
!                               n_L1, impact, bangle, bangle_tl)
!
! DESCRIPTION
!   Forward model the ionospheric bending for the L1 and L2 bending angles
!   assuming Chapman layer profiles. Should only be called if
!   x%direct_ion=.TRUE.
!
! INPUTS
!   REAL(wp), DIMENSION(:)  :: Ne_max     ! Peak electron density (m-3)
!   REAL(wp), DIMENSION(:)  :: R_peak     ! Radius of peak value (m)
!   REAL(wp), DIMENSION(:)  :: H_width    ! Width of Chapman Layer (m)
!   INTEGER,                :: n_L1       ! No. of L1 bending angles in array
!   REAL(wp), DIMENSION(:)  :: impact     ! Impact parameters
!   REAL(wp), DIMENSION(:)  :: bangle     ! Neutral bending angle
!
! INOUT
!   REAL(wp), DIMENSION(:) :: bangle_tl   ! Perturbed _total_ bending angle
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

SUBROUTINE ropp_fm_iono_bangle_tl( &
                                  Ne_max,    R_peak,    H_width, &
                                  Ne_max_tl, R_peak_tl, H_width_tl, &
                                  n_L1, impact, bangle, bangle_tl)

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm,   not_this => ropp_fm_iono_bangle_tl
! USE ropp_fm
  USE ropp_utils, ONLY: ropp_MDTV
! USE ropp_fm_constants, ONLY: k4, f_L1,f_L2

  IMPLICIT NONE

  REAL(wp), INTENT(in)                  :: Ne_max      ! Peak electron density (m-3)
  REAL(wp), INTENT(in)                  :: R_peak      ! Radius of peak value (m)
  REAL(wp), INTENT(in)                  :: H_width     ! Width of Chapman Layer (m)
  REAL(wp), INTENT(in)                  :: Ne_max_tl   ! Delta peak electron density (m-3)
  REAL(wp), INTENT(in)                  :: R_peak_tl   ! Delta Radius of peak value (m)
  REAL(wp), INTENT(in)                  :: H_width_tl  ! Delta width of Chapman Layer (m)
  INTEGER                               :: n_L1        ! Number of L1 bending angles
  REAL(wp), DIMENSION(:), INTENT(in)    :: impact      ! Impact parameters
  REAL(wp), DIMENSION(:), INTENT(in)    :: bangle      ! Neutral bending angle
  REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_tl   ! Delta _total_ bending angle

  INTEGER                               :: n_chap, n_impact, i
  REAL(wp)                              :: const_L1, const_L2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: lg,    z,    q, &
                                           lg_tl, z_tl, q_tl
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: bangle_ion_tl

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_chap = 1 ! For now
  n_impact = SIZE(impact)

  ALLOCATE(lg(n_impact),    z(n_impact),    q(n_impact))
  ALLOCATE(lg_tl(n_impact), z_tl(n_impact), q_tl(n_impact))
  ALLOCATE(bangle_ion_tl(n_impact))

  bangle_ion_tl = 0.0_wp

!-------------------------------------------------------------------------------
! 3. Calculate the Zorro function
!-------------------------------------------------------------------------------

  lg = (R_peak - impact) / H_width

  lg_tl = (R_peak_tl/H_width) - (lg/H_width)*H_width_tl

  z = ropp_fm_zorro(lg)

! We use the gradient of Zorro w.r.t lg here and in adjoint

  z_tl = ropp_fm_dzorro_dlg(lg)*lg_tl

!--------------------------------------------------------------------------------
! 4. Ionospheric bending BUT NOT YET including scaling by frequency
!--------------------------------------------------------------------------------  

  q = 2.0_wp * EXP(0.5_wp) * impact * R_peak / &
      SQRT(H_width*(R_peak + impact)**3)

  q_tl = -0.5_wp * (q/H_width) * H_width_tl + &
         q * R_peak_tl * &
         ((impact - 0.5_wp*R_peak) / R_peak / (R_peak + impact))

  bangle_ion_tl = bangle_ion_tl           +  &
                  Ne_max_tl * q    * z    +  &
                  Ne_max    * q_tl * z    +  &
                  Ne_max    * q    * z_tl

!--------------------------------------------------------------------------------
! 5. Now calculate total bending = (neutral + ionospheric)
!--------------------------------------------------------------------------------

  const_L1 = k4 / f_L1**2
  const_L2 = k4 / f_L2**2

  DO i = 1, n_impact

    IF (bangle(i) < ropp_MDTV) CYCLE  ! neutral is missing

    IF (i <= n_L1) THEN ! the L1 bending angles

      bangle_tl(i) = bangle_tl(i) + const_L1 * bangle_ion_tl(i)

    ELSE                ! the L2 bending angles

      bangle_tl(i) = bangle_tl(i) + const_L2 * bangle_ion_tl(i)

    END IF

  END DO

!--------------------------------------------------------------------------------
! 6. Deallocate/tidy up
!--------------------------------------------------------------------------------

  DEALLOCATE(lg,    z,    q)
  DEALLOCATE(lg_tl, z_tl, q_tl)
  DEALLOCATE(bangle_ion_tl)

  RETURN

END SUBROUTINE ropp_fm_iono_bangle_tl

