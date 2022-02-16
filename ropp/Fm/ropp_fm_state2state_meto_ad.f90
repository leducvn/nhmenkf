! $Id: ropp_fm_state2state_meto_ad.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_meto_ad(x, x_ad)

!****s* Model_meto/ropp_fm_state2state_meto_ad *
!
! NAME
!    ropp_fm_state2state_meto_ad - Adjoint of ropp_fm_state2state_meto
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_meto_ad(x, x_ad)
! 
! DESCRIPTION
!    This routine is the adjoint of ropp_fm_state2state_meto
!
! INPUTS
!    type(State1DFM)  :: x       State vector structure
!    type(State1DFM)  :: x_ad    State vector adjoint structure
!
! OUTPUT
!    type(State1DFM)  :: x_ad    State vector adjoint structure
!                             
! NOTES
!
! SEE ALSO
!    ropp_fm_state2state_meto
!    ropp_fm_state2state_meto_tl
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

  TYPE(State1dFM), INTENT(in)         :: x           ! State vector
  TYPE(State1dFM), INTENT(inout)      :: x_ad        ! State adjoint

  REAL(wp), DIMENSION(:), ALLOCATABLE :: pressA      ! Pressure on A-lvls
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pressA_ad   ! Pressure perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geopA       ! Geop ht on A-levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerA      ! Exner pressure on A-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerA_ad   ! ExnerA perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerB      ! Exner pressure on B-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ExnerB_ad   ! ExnerB pertubation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv       ! Virtual temp on B-lvl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv_ad    ! Virtual temp pertubation
  REAL(wp)                            :: alpha       ! Interpolation factor
  INTEGER                             :: i           ! Counter
  REAL(wp), PARAMETER  :: Pref=100000.0_wp           ! Reference pressure (Pa)

!-------------------------------------------------------------------------------
! 2. Memory allocation
!-------------------------------------------------------------------------------

  ALLOCATE(pressA(SIZE(x%geop)+1))
  ALLOCATE(pressA_ad(SIZE(x%geop)+1))
  ALLOCATE(geopA(SIZE(x%geop)+1))
  ALLOCATE(ExnerA(SIZE(x%geop)+1))
  ALLOCATE(ExnerA_ad(SIZE(x%geop)+1))
  ALLOCATE(ExnerB(SIZE(x%geop))) 
  ALLOCATE(ExnerB_ad(SIZE(x%geop))) 
  ALLOCATE(Tvflv(SIZE(x%geop)))
  ALLOCATE(Tvflv_ad(SIZE(x%geop)))

!-------------------------------------------------------------------------------
! 3. Reset local adjoint variables
!-------------------------------------------------------------------------------

  pressA_ad(:) = 0.0_wp
  ExnerA_ad(:) = 0.0_wp
  ExnerB_ad(:) = 0.0_wp
  Tvflv_ad(:)  = 0.0_wp

!-------------------------------------------------------------------------------
! 4. Update state variables
!-------------------------------------------------------------------------------

  IF (x%use_logp) THEN
     pressA = EXP(x%state(1:x%n_lev+1)) * 100.0_wp
  ELSE
     pressA = x%state(1:x%n_lev+1)
  ENDIF
  
!-------------------------------------------------------------------------------
! 5. Recompute intermediate variables
!-------------------------------------------------------------------------------

! 5.1 Calculate geopotential height on A (density levels)

  geopA(1) = x%geop_sfc
  DO i=2,x%n_lev
     geopA(i) = x%geop(i-1)+((x%geop(i)-x%geop(i-1))/2.0_wp)
  ENDDO
  geopA(x%n_lev+1)=x%geop(x%n_lev)+((x%geop(x%n_lev)-x%geop(x%n_lev-1))/2.0_wp)

! 5.2 Calculate Exner on A (density) levels

  ExnerA = (pressA/Pref)**(R_dry / C_p)

! 5.3 Calculate Exner on B (theta) levels (assume array index 1 towards sfc)

  DO i=1,x%n_lev 
     alpha = (geopA(i+1)-x%geop(i))/(geopA(i+1)-geopA(i)) 

     ExnerB(i) = alpha*ExnerA(i) + ((1.0_wp - alpha)*ExnerA(i+1)) 
  ENDDO

! 5.4 Calculate mean layer virtual temperature

  DO i=1,x%n_lev    
     Tvflv(i) = (g_wmo * (geopA(i+1)-geopA(i)) * ExnerB(i)) /   &
                     (C_p * (ExnerA(i)-ExnerA(i+1)))
  ENDDO

!-------------------------------------------------------------------------------
! 6. Adjoint of computing temperature and pressure on full levels
!-------------------------------------------------------------------------------

! 6.1 Adjoint of computing temperature on B (theta) levels

  x_ad%shum = x_ad%shum - ((1.0_wp/epsilon_water)-1.0_wp) &
            * Tvflv / (1.0_wp+(((1.0_wp/epsilon_water)-1.0_wp)*x%shum)) &
            * x_ad%temp
  Tvflv_ad = Tvflv_ad + x%temp / Tvflv * x_ad%temp
  x_ad%temp = 0.0_wp

! 6.2 Adjoint of computing mean layer virtual temperature

  DO i=1,x%n_lev
     ExnerA_ad(i+1)=ExnerA_ad(i+1)+Tvflv(i)/(ExnerA(i)-ExnerA(i+1))*Tvflv_ad(i)
     ExnerA_ad(i) = ExnerA_ad(i) - Tvflv(i)/(ExnerA(i)-ExnerA(i+1))*Tvflv_ad(i)
     ExnerB_ad(i) = ExnerB_ad(i) + Tvflv(i) / ExnerB(i) * Tvflv_ad(i)
     Tvflv_ad(i) = 0.0_wp
  ENDDO

! 6.3 Adjoint of computing pressure on B (theta) levels
  
  ExnerB_ad = ExnerB_ad +      &
               (C_p/R_dry) * Pref * (ExnerB**((C_p/R_dry) - 1.0_wp)) * x_ad%pres
  x_ad%pres   = 0.0_wp

! 6.4 Adjoint of computing Exner on B (theta) levels 

  DO i=1,x%n_lev   
     alpha = (geopA(i+1)-x%geop(i))/(geopA(i+1)-geopA(i)) 
     ExnerA_ad(i+1) = ExnerA_ad(i+1) + (1.0_wp - alpha) * ExnerB_ad(i)
     ExnerA_ad(i)   = ExnerA_ad(i) + alpha * ExnerB_ad(i)
     ExnerB_ad(i)   = 0.0_wp
  ENDDO

! 6.5 Adjoint of computing Exner on A (density) levels

  pressA_ad = pressA_ad + (R_dry / C_p) * ExnerA / pressA * ExnerA_ad
  ExnerA_ad = 0.0_wp

!-------------------------------------------------------------------------------
! 7. Adjoint of state vector copying
!-------------------------------------------------------------------------------


  IF (x%direct_ion) THEN
    x_ad%state(2*x%n_lev+2) = x_ad%state(2*x%n_lev+2) + x_ad%ne_max
    x_ad%ne_max = 0.0_wp
    x_ad%state(2*x%n_lev+3) = x_ad%state(2*x%n_lev+3) + x_ad%h_peak
    x_ad%h_peak = 0.0_wp
    x_ad%state(2*x%n_lev+4) = x_ad%state(2*x%n_lev+4) + x_ad%h_width
    x_ad%h_width = 0.0_wp
  ENDIF

  IF (x%use_logq) THEN
     x_ad%state(x%n_lev+2:2*x%n_lev+1) = x_ad%state(x%n_lev+2:2*x%n_lev+1) +  &
                               x_ad%shum * EXP(x%state(x%n_lev+2:2*x%n_lev+1)) / 1000.0_wp 
                               
  ELSE
     WHERE (x%shum <= 0.0_wp)
        x_ad%shum = 0.0_wp
     endwhere
     x_ad%state(x%n_lev+2:2*x%n_lev+1) = x_ad%state(x%n_lev+2:2*x%n_lev+1) +  &
                                         x_ad%shum
  ENDIF
  x_ad%shum = 0.0_wp
  
  IF (x%use_logp) THEN
     x_ad%state(1:x%n_lev+1) = x_ad%state(1:x%n_lev+1) +           &
                               pressA_ad * EXP(x%state(1:x%n_lev+1)) * 100.0_wp
  ELSE
     x_ad%state(1:x%n_lev+1) = x_ad%state(1:x%n_lev+1) + pressA_ad
  ENDIF
  pressA_ad = 0.0_wp

!-------------------------------------------------------------------------------
! 8. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(pressA)
  DEALLOCATE(pressA_ad)
  DEALLOCATE(geopA)
  DEALLOCATE(ExnerA)
  DEALLOCATE(ExnerA_ad)
  DEALLOCATE(ExnerB)
  DEALLOCATE(ExnerB_ad)
  DEALLOCATE(Tvflv)
  DEALLOCATE(Tvflv_ad)

END SUBROUTINE ropp_fm_state2state_meto_ad
