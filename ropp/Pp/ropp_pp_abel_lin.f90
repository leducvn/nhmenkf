! $Id: ropp_pp_abel_lin.f90 3551 2013-02-25 09:51:28Z idculv $

SUBROUTINE ropp_pp_abel_LIN(nr, refrac, impact, bangle, dln, scale)

!****s* BendingAngle/ropp_pp_abel_LIN *
!
! NAME
!    ropp_pp_abel_LIN - Calculate a one dimensional bending angle
!                       profile from refractivity / impact parameter profile at
!                       required output levels using a Fast Abel Transform
!                       Assume linear variation of ln(n) with impact height.  
!
! SYNOPSIS
!    call ropp_pp_abel(nr, refrac, impact, bangle, dlndx, scale)
! 
! DESCRIPTION
!    This routine calculates bending angles at a given set of impact parameters
!    from a vertical profile of refractivity given at the observation set of
!    impact parameters.
!
! INPUTS
!    real(wp), dimension(:) :: nr          ! x=nr product
!    real(wp), dimension(:) :: refrac      ! Refractivity values
!    real(wp), dimension(:) :: impact      ! Impact parameters
!    real(wp), dimension(:) :: dlndx       ! Refractive index gradient
!    real(wp),              :: scale       ! Scale height
!
! OUTPUT
!    real(wp), dimension(:) :: bangle       ! Calculated bending angles
!
! NOTES
!    The interpolation of bending angles calculated at the input data
!    geopotential levels to the output impact parameter levels is
!    carried out assuming that bending angle varies exponentially with
!    impact parameter.
!
! SEE ALSO
!    ropp_pp_abel_EXP
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
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
  REAL(wp), OPTIONAL, DIMENSION(:), INTENT(in)  :: dln ! d(ln[n])/dx gradient
  REAL(wp), OPTIONAL,     INTENT(in)  :: scale         ! scale height
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle        ! bending angle

  REAL(wp), DIMENSION(:), ALLOCATABLE :: dlndx         ! d(ln[n])/dx gradient
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ref_n         ! refractive index
  REAL(wp), DIMENSION(:), ALLOCATABLE :: kval          ! linear decay rate
  REAL(wp)                            :: t_upper       ! upper bound on integral
  REAL(wp)                            :: t_lower       ! lower bound on integral
  REAL(wp)                            :: refrac_low    ! refrac at lower level
  REAL(wp)                            :: nr_low        ! x=nr at lower level
  REAL(wp)                            :: delta_refrac
  REAL(wp)                            :: delta_nr
  REAL(wp)                            :: erf_up
  REAL(wp)                            :: zt

  REAL(wp), PARAMETER                 :: a = 0.3480242_wp
  REAL(wp), PARAMETER                 :: b = 0.0958798_wp
  REAL(wp), PARAMETER                 :: c = 0.7478556_wp

  INTEGER                             :: n_lev, n_lower, n_impact
  INTEGER                             :: i, i_bot, i_top, l, k, im1, ip1

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  ! 2.1 Array sizes
  
  n_lev = SIZE(nr)
  n_impact = SIZE(impact)

  ALLOCATE(dlndx(n_lev))
  ALLOCATE(ref_n(n_lev))
  ALLOCATE(kval(n_lev-1))
  
  ! 2.2 Compute gradient of refractive index (if required)
  
  IF(PRESENT(dln))THEN
     
     dlndx = dln
     
  ELSE

     WHERE(refrac > ropp_MDTV)
       ref_n = LOG(1.0_wp + 1.e-6_wp*refrac)
     ELSEWHERE
       ref_n = ropp_MDFV
     ENDWHERE
     
     DO i=1,n_lev
     
       ip1 = i+1
       IF(i == n_lev) ip1 = i
       
       im1 = i-1
       IF(i == 1) im1 = i

       IF (ref_n(ip1) > ropp_MDTV .AND. ref_n(im1) > ropp_MDTV) THEN
         dlndx(i) = (ref_n(ip1) - ref_n(im1)) / (nr(ip1) - nr(im1))
       ELSE
         dlndx(i) = ropp_MDFV
       ENDIF
       
     ENDDO
       
  ENDIF

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
! 4. Calculate linear decay rate between levels
!-------------------------------------------------------------------------------

  DO i = 1, n_lev - 1
     
     IF (dlndx(i) > ropp_ZERO .OR. dlndx(i+1) > ropp_ZERO &
         .OR. dlndx(i) < ropp_MDTV .OR. dlndx(i+1) < ropp_MDTV) THEN
        kval(i) = 1.0e-24_wp
     ELSE
! Should probably use linear decay rate in this routine
!        kval(i) = LOG(dlndx(i)/dlndx(i+1)) / MAX(1.0_wp, (nr(i+1)-nr(i)))
        kval(i) = -(dlndx(i)  - dlndx(i+1)) / MAX(1.0_wp, (nr(i+1)-nr(i))) !NB: dlndx < 0
        kval(i) = MAX(1.0e-24_wp, kval(i))
     ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! 5. Calculate bending angles for observational heights
!-------------------------------------------------------------------------------

  bangle(:) = ropp_MDFV

  DO l = 1, n_impact

     IF (impact(l) < nr(n_lower) .OR. impact(l) > nr(n_lev)) THEN
        CYCLE
     ENDIF

!    5.1 Find bottom state vector level
!    ----------------------------------

     i_bot = n_lower

     DO k=i_bot,n_lev-1  
       IF (impact(l) < nr(k + 1)) EXIT    
     ENDDO
     
     i_bot = k

     i_top = n_lev
     DO WHILE(refrac(i_top) <= ropp_MDTV .OR. dlndx(i_top) <= ropp_MDTV)
        i_top = i_top - 1
     ENDDO

!    5.2 Loop over all levels above
!    ------------------------------

     bangle(l) = ropp_ZERO
     
     DO i = i_bot, i_top - 1

        IF(refrac(i) <= ropp_MDTV .OR. refrac(i+1) <= ropp_MDTV   &
             .OR. dlndx(i) <= ropp_MDTV .OR. dlndx(i+1) <= ropp_MDTV ) &
             CYCLE

!       5.2.1 Values of refractivity and impact parameter at lower level
!       ----------------------------------------------------------------

        IF (i == i_bot) THEN
! Should probably use linear decay rate in this routine
!           refrac_low = dlndx(i_bot) * EXP(-kval(i_bot) * (impact(l)-nr(i_bot)))
           refrac_low = dlndx(i_bot) + kval(i_bot) * (impact(l)-nr(i_bot))
           nr_low = impact(l)
        ELSE
           refrac_low = dlndx(i) 
           nr_low = nr(i) 
        ENDIF

!       5.2.2 Analytical Abel integral solution
!       ---------------------------------------

        delta_refrac = dlndx(i+1) - refrac_low
        delta_nr     = nr(i+1) - nr_low
 
        t_upper = SQRT((nr(i+1)-impact(l))*(nr(i+1)+impact(l)))
        t_lower = SQRT((nr_low-impact(l))*(nr_low+impact(l)))

        bangle(l) = bangle(l) +                                      &
                     ((refrac_low*nr(i+1)-dlndx(i+1)*nr_low)    &
                     * LOG((t_upper+nr(i+1))/(t_lower+nr_low)) +     & 
                     delta_refrac*(t_upper-t_lower))/delta_nr

     ENDDO

!    5.3 Asymptotic correction
!    -------------------------

     IF(PRESENT(scale))THEN

        t_upper = SQRT((nr(i_top)-impact(l))/scale)
        zt = 1.0_wp / (1.0_wp + 0.47047_wp * t_upper)
        erf_up = 1.0_wp - (a-(b-c*zt)*zt) * zt * EXP(-(t_upper*t_upper))

        bangle(l) = bangle(l) + &
                     dlndx(i_top)*EXP((nr(i_top)-impact(l))/scale) *    &
                     SQRT(pi)*(1.0_wp-erf_up)/                         &
                     SQRT((nr(i_top)+impact(l))/scale)


      ENDIF
     
     bangle(l)  = -2.0_wp * impact(l) * bangle(l)  

     IF(bangle(l) == ropp_ZERO) bangle(l) = ropp_MDFV
     
  ENDDO

  DEALLOCATE(dlndx)
  DEALLOCATE(ref_n)
  DEALLOCATE(kval)
  
END SUBROUTINE ropp_pp_abel_LIN
