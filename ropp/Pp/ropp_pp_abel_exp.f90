! $Id: ropp_pp_abel_exp.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_abel_EXP(nr, refrac, impact, bangle)

!****s* BendingAngle/ropp_pp_abel_EXP *
!
! NAME
!    ropp_pp_abel_EXP - Calculate a one dimensional bending angle
!                       profile from refractivity / impact parameter profile at
!                       required output levels using a Fast Abel Transform
!                       Assume exponential variation of ln(n) with height 
!                       between successive impact parameter levels.
!
! SYNOPSIS
!    call ropp_pp_abel_EXP(nr, refrac, impact, bangle)
! 
! DESCRIPTION
!    This routine calculates bending angles at a given set of impact parameters
!    from a vertical profile of refractivity given at the observed set of 
!    x = nr levels.
!
! INPUTS
!    real(wp), dimension(:) :: nr          ! x=nr product
!    real(wp), dimension(:) :: refrac      ! Refractivity values
!    real(wp), dimension(:) :: impact      ! Impact parameters
!
! OUTPUT
!    real(wp), dimension(:) :: bangle      ! Calculated bending angles
!
! NOTES
!    The interpolation of bending angles calculated at the input data
!    geopotential levels to the output impact parameter levels is
!    carried out assuming that bending angle varies exponentially with
!    impact parameter.
!
! SEE ALSO
!     ropp_pp_invert_EXP
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
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils, ONLY: ropp_MDTV, ropp_MDFV, ropp_ZERO
  USE ropp_pp_constants, ONLY: pi

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)  :: nr            ! x=nr product
  REAL(wp), DIMENSION(:), INTENT(in)  :: refrac        ! refractivity
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact        ! impact parameter
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle        ! bending angle

  REAL(wp), DIMENSION(:), ALLOCATABLE :: kval          ! exponential decay rate
  REAL(wp)                            :: t_upper       ! upper bound on integral
  REAL(wp)                            :: t_lower       ! lower bound on integral
  REAL(wp)                            :: refrac_low    ! refrac at lower level
  REAL(wp)                            :: nr_low        ! x=nr at lower level

  REAL(wp)                            :: integral_erf  ! integral approximation
  REAL(wp)                            :: erf_up, erf_low, diff_ip
  REAL(wp)                            :: zt
  REAL(wp), PARAMETER                 :: a = 0.3480242_wp
  REAL(wp), PARAMETER                 :: b = 0.0958798_wp
  REAL(wp), PARAMETER                 :: c = 0.7478556_wp

  INTEGER                             :: n_lev, n_lower, n_impact
  INTEGER                             :: i, i_bot, i_top, l

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_lev = SIZE(nr)
  n_impact = SIZE(impact)

  ALLOCATE(kval(n_lev-1))

!-------------------------------------------------------------------------------
! 3. Calculate lowest usable level (because of superrefraction)
!-------------------------------------------------------------------------------

  n_lower = 1
  DO i = n_lev, 2, -1
     IF (nr(i) <= nr(i-1)) THEN
        n_lower = i
        EXIT
     ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! 4. Calculate exponential decay rate between levels
!-------------------------------------------------------------------------------

  DO i = 1, n_lev - 1

     IF (refrac(i) > ropp_MDTV .AND. refrac(i+1) > ropp_MDTV) THEN

       kval(i) = LOG(refrac(i)/refrac(i+1)) / MAX(1.0_wp, (nr(i+1)-nr(i)))

       kval(i) = MAX(1.0e-6_wp, kval(i))
       
        ! avoid problems with noisy data towards top of profile
        IF(kval(i) == 1.0e-6_wp)THEN
           kval(i:n_lev-1) = 1.0e-6_wp
           EXIT
        ENDIF
        
     ELSE
        
        kval(i) = 1.0e-6_wp

     ENDIF
       
  ENDDO

!-------------------------------------------------------------------------------
! 5. Calculate bending angles for observational heights
!-------------------------------------------------------------------------------

  bangle(:) = ropp_MDFV

  DO l = 1, n_impact
     
     IF (impact(l) < nr(n_lower) .OR. impact(l) >= nr(n_lev)) THEN
       CYCLE
     ENDIF
     
!    5.1 Find bottom state vector level
!    ----------------------------------

     i_bot = n_lower
     DO WHILE (impact(l) >= nr(i_bot + 1) .AND. i_bot < n_lev)
        i_bot = i_bot + 1
     ENDDO

     i_top = n_lev
     DO WHILE(refrac(i_top) <= ropp_MDTV)
        i_top = i_top - 1
     ENDDO

!    5.2 Loop over all levels above
!    ------------------------------

     bangle(l) = ropp_ZERO

     DO i = i_bot, i_top - 1

        IF(refrac(i) <= ropp_MDTV) CYCLE

!       5.2.1 Values of refractivity and impact parameter at lower level
!       ----------------------------------------------------------------

        IF (i == i_bot) THEN
           refrac_low = refrac(i_bot) * EXP(-kval(i_bot)*(impact(l)-nr(i_bot)))
           nr_low = impact(l)
        ELSE
           refrac_low = refrac(i) 
           nr_low = nr(i) 
        ENDIF

!       5.2.2 Upper (100 km above top end) and lower bounds of the integral
!       -------------------------------------------------------------------

        IF (i == n_lev - 1) THEN
           diff_ip = nr(i+1) + 1.0d5 - impact(l)
           IF (diff_ip >= 0.0_wp) THEN
              t_upper = SQRT(kval(i) * diff_ip)
           ELSE
              t_upper = ropp_MDFV
           ENDIF
        ELSE
           diff_ip = nr(i+1) - impact(l)
           IF (diff_ip >= 0.0_wp) THEN
              t_upper = SQRT(kval(i) * diff_ip)
           ELSE
              t_upper = ropp_MDFV
           ENDIF
        ENDIF

        IF (i == i_bot) THEN
           t_lower = 0.0_wp
        ELSE
           diff_ip = nr(i) - impact(l)
           IF (diff_ip >= 0.0_wp) THEN
              t_lower = SQRT(kval(i) * diff_ip)
           ELSE
              t_lower = ropp_MDFV
           ENDIF
        ENDIF

!       5.2.3 Integral
!       --------------

        ! Approximate error function with polynomial
        IF (t_lower > ropp_MDTV) THEN
           zt = 1.0_wp / (1.0_wp + 0.47047_wp * t_lower)
           erf_low = 1.0_wp - (a-(b-c*zt)*zt) * zt * EXP(-(t_lower*t_lower))
        ELSE
           erf_low = ropp_MDFV
        ENDIF

        IF (t_upper > ropp_MDTV) THEN
           zt = 1.0_wp / (1.0_wp + 0.47047_wp*t_upper)
           erf_up = 1.0_wp - (a-(b-c*zt)*zt) * zt * EXP(-(t_upper*t_upper))
        ELSE
           erf_up = ropp_MDFV
        ENDIF

        IF ((erf_low > ropp_MDTV) .AND. (erf_up > ropp_MDTV)) THEN
           integral_erf = erf_up - erf_low
        ELSE
           integral_erf = ropp_MDFV
        ENDIF

!       5.2.4 Bending angle value
!       -------------------------

       IF (integral_erf > ropp_MDTV) THEN
          bangle(l) = bangle(l) +   &
                     1.0e-6_wp * SQRT(2.0_wp*pi*impact(l)*kval(i))       &
                     * refrac_low * EXP(kval(i) * (nr_low - impact(l)))  &
                     * integral_erf
       ENDIF

     ENDDO
     
     IF(bangle(l) == ropp_ZERO) bangle(l) = ropp_MDFV

  ENDDO

  DEALLOCATE(kval)

END SUBROUTINE ropp_pp_abel_EXP
