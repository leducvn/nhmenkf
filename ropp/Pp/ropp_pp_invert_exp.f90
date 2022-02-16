! $Id: ropp_pp_invert_exp.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_invert_EXP(impact, bangle, nr, refrac)

!****s* BendingAngle/ropp_pp_invert_EXP *
!
! NAME
!    ropp_pp_invert_EXP - Calculate a one dimensional refractivity profile from
!                         bending angle / impact parameter profile using a 
!                         Fast Abel Transform. Assume exponential variation of 
!                         bending angle with height between successive impact 
!                         parameter levels.
!
! SYNOPSIS
!    call ropp_pp_invert_EXP(impact, bangle, nr, refrac)
! 
! DESCRIPTION
!    This routine calculates refractivity at a given set of x=nr levels
!    from a vertical profile of bending angles given at a set 
!    of impact parameters.
!
! INPUTS
!    real(wp), dimension(:) :: impact      ! Input impact parameters 
!    real(wp), dimension(:) :: bangle      ! Bending angles
!    real(wp), dimension(:) :: nr          ! x = nr product
!
! OUTPUT
!    real(wp), dimension(:) :: refrac      ! Refractivity values
!
! NOTES
!    The interpolation of refractivity calculated at the input data
!    impact height levels to the output geopotential levels is
!    carried out assuming that dln(n)/dx varies exponentially with x.
!
! SEE ALSO
!    ropp_pp_invert_LIN
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

  REAL(wp), DIMENSION(:), INTENT(in)  :: impact       ! impact parameter
  REAL(wp), DIMENSION(:), INTENT(in)  :: bangle       ! bending angle
  REAL(wp), DIMENSION(:), INTENT(in)  :: nr           ! x=nr product for output
  REAL(wp), DIMENSION(:), INTENT(out) :: refrac       ! refractivity

  REAL(wp), DIMENSION(:), ALLOCATABLE :: kval         ! exponential decay rate
  REAL(wp)                            :: t_upper      ! upper bound on integral
  REAL(wp)                            :: t_lower      ! lower bound on integral
  REAL(wp)                            :: bangle_low   ! bangle at lower level
  REAL(wp)                            :: impact_low   ! impact at lower level

  REAL(wp)                            :: integral_erf ! integral approximation
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
  
  ALLOCATE(kval(n_impact-1))

!-------------------------------------------------------------------------------
! 3. Calculate lowest usable level (because of superrefraction)
!-------------------------------------------------------------------------------

  n_lower = 1
  DO i = n_impact, 2, -1
     IF (impact(i) <= impact(i-1)) THEN
        n_lower = i
        EXIT
     ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! 4. Calculate exponential decay rate between levels
!-------------------------------------------------------------------------------

  DO i = 1, n_impact - 1

     IF (bangle(i) > ropp_MDTV .AND. bangle(i+1) > ropp_MDTV) THEN
     
        kval(i) = LOG(bangle(i)/bangle(i+1))/MAX(1.0_wp,(impact(i+1)-impact(i)))
        kval(i) = MAX(1.0e-6_wp, kval(i))
     
        ! avoid problems with noisy data towards top of profile
        IF(kval(i) == 1.0e-6_wp)THEN
           kval(i:n_impact-1) = 1.0e-6_wp
           EXIT
        ENDIF
        
     ELSE
        
        kval(i) = 1.0e-6_wp

     ENDIF
     
  ENDDO

!-------------------------------------------------------------------------------
! 5. Calculate bending angles for observational heights
!-------------------------------------------------------------------------------

  refrac(:) = ropp_MDFV

  DO l = 1, n_lev

     IF (bangle(l) <= ropp_MDTV) THEN
        CYCLE
     ENDIF

     IF (nr(l) < impact(n_lower) .OR. nr(l) >= impact(n_impact)) THEN
        CYCLE
     ENDIF

!    5.1 Find bottom state vector level
!    ----------------------------------

     i_bot = n_lower
     DO WHILE (nr(l) >= impact(i_bot + 1) .AND. i_bot < n_impact)
        i_bot = i_bot + 1
     ENDDO

     i_top = n_impact
     DO WHILE(bangle(i_top) <= ropp_MDTV)
        i_top = i_top - 1
     ENDDO

!    5.2 Loop over all levels above
!    ------------------------------

     refrac(l) = ropp_ZERO
     
     DO i = i_bot, i_top - 1
        
!       5.2.1 Values of refractivity and impact parameter at lower level
!       ----------------------------------------------------------------

        IF (i == i_bot) THEN
           bangle_low = bangle(i_bot)*EXP(-kval(i_bot)*(nr(l)-impact(i_bot)))
           impact_low = nr(l)
        ELSE
           bangle_low = bangle(i)
           impact_low = impact(i) 
        ENDIF

!       5.2.2 Upper (100 km above top end) and lower bounds of the integral
!       -------------------------------------------------------------------

        IF (i == i_top - 1) THEN
           diff_ip = impact(i+1) + 1.0d5 - nr(l)
           IF (diff_ip >= 0.0_wp) THEN
              t_upper = SQRT(kval(i) * diff_ip)
           ELSE
              t_upper = ropp_MDFV
           ENDIF
        ELSE
           diff_ip = impact(i+1) - nr(l)
           IF (diff_ip >= 0.0_wp) THEN
              t_upper = SQRT(kval(i) * diff_ip)
           ELSE
              t_upper = ropp_MDFV
           ENDIF
        ENDIF

        IF (i == i_bot) THEN
           t_lower = 0.0_wp
        ELSE
           diff_ip = impact(i) - nr(l)
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

!       5.2.4 Refractivity value
!       -------------------------

       IF (integral_erf > ropp_MDTV) THEN
          refrac(l) = refrac(l) + &
                      (bangle_low * EXP( kval(i) * (impact_low - nr(l))) *  &
                      integral_erf ) / SQRT(2.0_wp * pi * nr(l) * kval(i))
       ENDIF

     ENDDO

     IF(refrac(l) == ropp_ZERO) refrac(l) = ropp_MDFV
     IF(refrac(l) >  ropp_ZERO) refrac(l) = (1E6_wp)*(EXP(refrac(l))-1.0_wp)

  ENDDO

  DEALLOCATE(kval)

END SUBROUTINE ropp_pp_invert_EXP
