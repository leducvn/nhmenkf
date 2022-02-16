! $Id: ropp_fm_state2state_meto.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_meto(x)

!****s* Model_meto/ropp_fm_state2state_meto *
!
! NAME
!    ropp_fm_state2state_meto - Calculate pressure and temperature on B-levels 
!                               for MetOffice bg data and update state vector
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_meto(x)
! 
! DESCRIPTION
!    This subroutine calculates pressure and temperature on B-levels for a 
!    geopotential height vertical coordinate (e.g. Met Office UM). 
!
! INPUTS
!    type(State1DFM)  :: x       State vector structure
!
! OUTPUT
!    type(State1DFM)  :: x       Updated state vector structure
!                             
! NOTES
!
! SEE ALSO
!    ropp_fm_state2state_meto_tl
!    ropp_fm_state2state_meto_ad
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
! ---------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_constants
  USE ropp_fm_types

  IMPLICIT NONE

  TYPE(State1dFM), INTENT(inout)      :: x           ! State vector

  REAL(wp), DIMENSION(:), ALLOCATABLE :: pressA      ! Pressure on A-lvls
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geopA       ! Geop ht on A-levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerA      ! Exner pressure on A-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerB      ! Exner pressure on B-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv       ! Virtual temp on B-lvl
  REAL(wp)                            :: alpha       ! Interpolation factor
  INTEGER                             :: i           ! Counter
  REAL(wp), PARAMETER :: Pref=100000.0_wp            ! Reference pressure (Pa)

!-------------------------------------------------------------------------------
! 2. Memory allocation
!-------------------------------------------------------------------------------

  ALLOCATE(pressA(SIZE(x%geop)+1))
  ALLOCATE(geopA(SIZE(x%geop)+1))
  ALLOCATE(ExnerA(SIZE(x%geop)+1))
  ALLOCATE(ExnerB(SIZE(x%geop))) 
  ALLOCATE(Tvflv(SIZE(x%geop)))
  
!-------------------------------------------------------------------------------
! 3. Update state variables
!-------------------------------------------------------------------------------

  IF (x%use_logp) THEN
     pressA = EXP(x%state(1:x%n_lev+1)) * 100.0_wp
  ELSE
     pressA = x%state(1:x%n_lev+1)
  ENDIF

  IF (x%use_logq) THEN
    x%shum = EXP(x%state(x%n_lev + 2: 2*x%n_lev+1)) / 1000.0_wp
    IF (x%check_qsat) CALL check_qsat(x%shum, x%temp, x%pres)
  ELSE
    x%shum = x%state(x%n_lev + 2: 2*x%n_lev+1)
    IF (x%check_qsat) CALL check_qsat(x%shum, x%temp, x%pres)
    WHERE (x%shum <= 0.0_wp)
      x%shum = 1.0e-9_wp
    ENDWHERE
  ENDIF
  
  IF (x%use_logq) THEN
    x%state(x%n_lev + 2: 2*x%n_lev+1) = LOG(x%shum * 1000.0_wp)
  ELSE
    x%state(x%n_lev + 2: 2*x%n_lev+1) = x%shum
  ENDIF
  
  IF (x%direct_ion) THEN
    x%ne_max  = x%state(2*x%n_lev+2)
    x%h_peak  = x%state(2*x%n_lev+3)
    x%h_width = x%state(2*x%n_lev+4)
  ENDIF

!-------------------------------------------------------------------------------
! 4. Calculate geopotential height on A (density levels)
!-------------------------------------------------------------------------------

  geopA(1) = x%geop_sfc
  DO i=2,x%n_lev
     geopA(i) = x%geop(i-1) + ((x%geop(i)-x%geop(i-1))/2.0_wp)
  ENDDO
  geopA(x%n_lev+1)=x%geop(x%n_lev)+((x%geop(x%n_lev)-x%geop(x%n_lev-1))/2.0_wp)

!-------------------------------------------------------------------------------
! 5. Calculate Exner on A (density) levels
!-------------------------------------------------------------------------------

  ExnerA = (pressA/Pref)**(R_dry / C_p)

!-------------------------------------------------------------------------------
! 6. Calculate Exner on B (theta) levels (assume array index 1 towards sfc)
!-------------------------------------------------------------------------------

  DO i=1,x%n_lev    
     alpha = (geopA(i+1)-x%geop(i))/(geopA(i+1)-geopA(i)) 
     
     ExnerB(i) = alpha*ExnerA(i) + ((1.0_wp - alpha)*ExnerA(i+1)) 
  ENDDO

!-------------------------------------------------------------------------------
! 7. Calculate pressure on B (theta) levels
!-------------------------------------------------------------------------------

  x%pres = (ExnerB**(1.0_wp/(R_dry / C_p))) * Pref

!-------------------------------------------------------------------------------
! 8. Calculate mean layer virtual temperature
!-------------------------------------------------------------------------------

  DO i=1,x%n_lev
     Tvflv(i) = (g_wmo * (geopA(i+1)-geopA(i)) * ExnerB(i)) /   &
          (C_p * (ExnerA(i)-ExnerA(i+1)))
  ENDDO

!-------------------------------------------------------------------------------
! 9. Calculate temperature on B (theta) levels
!-------------------------------------------------------------------------------

  x%temp = Tvflv / ( 1.0_wp + (((1.0_wp/epsilon_water)-1.0_wp)*x%shum))
  
!-------------------------------------------------------------------------------
! 10. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(pressA)
  DEALLOCATE(geopA)
  DEALLOCATE(ExnerA)
  DEALLOCATE(ExnerB)
  DEALLOCATE(Tvflv)

CONTAINS
  
!-------------------------------------------------------------------------------
! 11. Compute Saturation Vapour Pressure - limit log(humidity) within saturation 
!-------------------------------------------------------------------------------

  SUBROUTINE check_qsat(shum, temp, press)

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(inout) :: shum
    REAL(wp), DIMENSION(:), INTENT(in)    :: temp
    REAL(wp), DIMENSION(:), INTENT(in)    :: press
    
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: Z
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: X
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: tfrac
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: es
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: qs

    ALLOCATE(es(SIZE(shum)))
    ALLOCATE(qs(SIZE(shum)))
    ALLOCATE(z(SIZE(shum)))
    ALLOCATE(x(SIZE(shum)))
    ALLOCATE(tfrac(SIZE(shum)))

! 11.1 Goff-Gratch equation (saturation vapour pressure over water)
    
    tfrac(:) = 373.16_wp / temp(:)
    Z = (-7.90298_wp * (tfrac - 1.0_wp)) +  &
       (5.02808_wp * LOG10(tfrac)) - (1.3816e-7_wp * (10**(11.344*(1.0_wp - (1.0_wp/tfrac))) - 1.0_wp)) + &
       (8.1328e-3_wp * (10**(-3.49149*(tfrac-1.0_wp)) - 1.0_wp))
    
! 11.2 Goff-Gratch equation (saturation vapour pressure over ice)

    tfrac(:) = 273.16_wp / temp(:)
    X = (-9.09718_wp * (tfrac - 1.0_wp)) + (-3.56654_wp*LOG10(tfrac)) + (0.876793_wp*(1.0_wp - (1.0_wp/tfrac)))

    WHERE(temp >= 273.16)
      es = 1013246.0_wp * (10**Z)
    ELSEWHERE
      es = 610.71_wp * (10**X)
    ENDWHERE
    
! 11.3 Compute saturation humidity

    qs = es * epsilon_water / (MAX(press,es) - ((1.0_wp - epsilon_water)*es))

    WHERE(qs == 1.0_wp)
      qs = 2.0_wp*shum(SIZE(shum))
    ENDWHERE

! 11.4 Limit humidity data below saturation

    WHERE(shum > qs)
      shum = qs
    ENDWHERE

    DEALLOCATE(es)
    DEALLOCATE(qs)
    DEALLOCATE(z)
    DEALLOCATE(x)
    DEALLOCATE(tfrac)

  END SUBROUTINE check_qsat

END SUBROUTINE ropp_fm_state2state_meto
