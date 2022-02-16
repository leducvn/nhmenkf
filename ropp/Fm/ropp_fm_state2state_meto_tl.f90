! $Id: ropp_fm_state2state_meto_tl.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_meto_tl(x, x_tl)

!****s* Model_meto/ropp_fm_state2state_meto_tl *
!
! NAME
!    ropp_fm_state2state_meto_tl - Tangent linear of ropp_fm_state2state_meto
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_meto_tl(x, x_tl)
! 
! DESCRIPTION
!    This routine is the tangent linear of ropp_fm_state2state_meto
!
! INPUTS
!    type(State1dFM)  :: x        State vector structure
!    type(State1dFM)  :: x_tl     Perturbation vector structure
!
! OUTPUT
!    type(State1dFM)  :: x_tl     Perturbation vector structure
!
! NOTES
!
! SEE ALSO
!    ropp_fm_state2state_meto
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

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_constants
  USE ropp_fm_types

  IMPLICIT NONE

  TYPE(State1dFM),        INTENT(in)  :: x           ! State vector
  TYPE(State1dFM),     INTENT(inout)  :: x_tl        ! Perturbation

  REAL(wp), DIMENSION(:), ALLOCATABLE :: pressA      ! Pressure on A-lvls
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pressA_tl   ! Pressure perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geopA       ! Geop ht on A-levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerA      ! Exner pressure on A-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerA_tl   ! ExnerA perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerB      ! Exner pressure on B-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerB_tl   ! ExnerB pertubation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv       ! Virtual temp on B-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv_tl    ! Virtual temp pertubation
  REAL(wp)                            :: alpha       ! Interpolation factor
  INTEGER                             :: i           ! Counter
  REAL(wp), PARAMETER :: Pref=100000.0_wp            ! Reference pressure (Pa)
 
!-------------------------------------------------------------------------------
! 2. Memory allocation
!-------------------------------------------------------------------------------

  ALLOCATE(pressA(SIZE(x%geop)+1))
  ALLOCATE(pressA_tl(SIZE(x%geop)+1))
  ALLOCATE(geopA(SIZE(x%geop)+1))
  ALLOCATE(ExnerA(SIZE(x%geop)+1))
  ALLOCATE(ExnerA_tl(SIZE(x%geop)+1))
  ALLOCATE(ExnerB(SIZE(x%geop))) 
  ALLOCATE(ExnerB_tl(SIZE(x%geop))) 
  ALLOCATE(Tvflv(SIZE(x%geop)))
  ALLOCATE(Tvflv_tl(SIZE(x%geop)))

!-------------------------------------------------------------------------------
! 3. Update state variables
!-------------------------------------------------------------------------------

  IF (x%use_logp) THEN
     pressA = EXP(x%state(1:x%n_lev+1)) * 100.0_wp
     pressA_tl = x_tl%state(1:x%n_lev+1) * EXP(x%state(1:x%n_lev+1)) * 100.0_wp
  ELSE
     pressA = x%state(1:x%n_lev+1)
     pressA_tl = x_tl%state(1:x%n_lev+1)
  ENDIF
  
  IF (x%use_logq) THEN
     x_tl%shum = x_tl%state(x%n_lev+2:2*x%n_lev+1) *    &
                    EXP(x%state(x%n_lev+2:2*x%n_lev+1)) / 1000.0_wp
  ELSE
     x_tl%shum = x_tl%state(x%n_lev+2:2*x%n_lev+1)
     WHERE(x%shum <= 0.0_wp)
        x_tl%shum = 0.0_wp
     endwhere
  ENDIF

  IF (x%direct_ion) THEN
    x_tl%ne_max  = x_tl%state(2*x_tl%n_lev+2)
    x_tl%h_peak  = x_tl%state(2*x_tl%n_lev+3)
    x_tl%h_width = x_tl%state(2*x_tl%n_lev+4)
  ENDIF

!-------------------------------------------------------------------------------
! 4. Calculate geopotential height on A (density levels)
!-------------------------------------------------------------------------------

  geopA(1) = x%geop_sfc
  DO i=2,x%n_lev
     geopA(i) = x%geop(i-1)+((x%geop(i)-x%geop(i-1))/2.0_wp)
  ENDDO
  geopA(x%n_lev+1)=x%geop(x%n_lev)+((x%geop(x%n_lev)-x%geop(x%n_lev-1))/2.0_wp)

!-------------------------------------------------------------------------------
! 5. Calculate Exner on A (density) levels
!-------------------------------------------------------------------------------

  ExnerA = (pressA/Pref)**(R_dry / C_p)
  ExnerA_tl = (R_dry / C_p) * ExnerA / pressA * pressA_tl

!-------------------------------------------------------------------------------
! 6. Calculate Exner on B (theta) levels (assume array index 1 towards sfc)
!-------------------------------------------------------------------------------

  DO i=1,x%n_lev    
     alpha = (geopA(i+1)-x%geop(i))/(geopA(i+1)-geopA(i)) 

     ExnerB(i) = alpha*ExnerA(i) + ((1.0_wp - alpha)*ExnerA(i+1)) 
     ExnerB_tl(i) = alpha*ExnerA_tl(i) + ((1.0_wp - alpha)*ExnerA_tl(i+1))
  ENDDO

!-------------------------------------------------------------------------------
! 7. Calculate pressure on B (theta) levels
!-------------------------------------------------------------------------------

!  x%pres = (ExnerB**(1.0_wp/(R_dry / C_p))) * Pref
  x_tl%pres = (C_p/R_dry) * Pref * (ExnerB**((C_p/R_dry)-1.0_wp)) * ExnerB_tl

!-------------------------------------------------------------------------------
! 8. Calculate mean layer virtual temperature
!-------------------------------------------------------------------------------

  DO i=1,x%n_lev  
     Tvflv(i) = (g_wmo * (geopA(i+1)-geopA(i)) * ExnerB(i)) /   &
                     (C_p * (ExnerA(i)-ExnerA(i+1)))
     Tvflv_tl(i) = Tvflv(i) / ExnerB(i) * ExnerB_tl(i) &
                     - Tvflv(i) / (ExnerA(i)-ExnerA(i+1)) * ExnerA_tl(i) &
                     +  Tvflv(i) / (ExnerA(i)-ExnerA(i+1)) * ExnerA_tl(i+1)
  ENDDO

!-------------------------------------------------------------------------------
! 9. Calculate temperature on B (theta) levels
!-------------------------------------------------------------------------------

!  x%temp = Tvflv / ( 1.0_wp + (((1.0_wp/epsilon_water)-1.0_wp)*x%shum))
  x_tl%temp = x%temp / Tvflv * Tvflv_tl -                                  &
              ((1.0_wp/epsilon_water)-1.0_wp) * Tvflv /                    &
              (1.0_wp+(((1.0_wp/epsilon_water)-1.0_wp)*x%shum))*x_tl%shum

!-------------------------------------------------------------------------------
! 10. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(pressA)
  DEALLOCATE(pressA_tl)
  DEALLOCATE(geopA)
  DEALLOCATE(ExnerA)
  DEALLOCATE(ExnerA_tl)
  DEALLOCATE(ExnerB)
  DEALLOCATE(ExnerB_tl)
  DEALLOCATE(Tvflv)
  DEALLOCATE(Tvflv_tl)

END SUBROUTINE ropp_fm_state2state_meto_tl
