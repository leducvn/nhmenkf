! $Id: ropp_fm_state2state_ecmwf_ad.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_ecmwf_ad(x, x_ad)

!****s* Model_ecmwf/ropp_fm_state2state_ecmwf_ad *
!
! NAME
!    ropp_fm_state2state_ecmwf_ad - Adjoint of ropp_fm_state2state_ecmwf
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_ecmwf_ad(x, x_ad)
! 
! DESCRIPTION
!    This routine is the adjoint of ropp_fm_state2state_ecmwf
!
! INPUTS
!    type(State1dFM)  :: x         State vector
!    type(State1dFM)  :: x_ad      State adjoint vector
!
! OUTPUT
!    type(State1dFM)  :: x_ad      State adjoint vector
!
! NOTES
!
! SEE ALSO
!   ropp_fm_state2state_ecmwf
!   ropp_fm_state2state_ecmwf_tl
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

  TYPE(State1dFM),    INTENT(in)      :: x        ! State vector
  TYPE(State1dFM),    INTENT(inout)   :: x_ad     ! State adjoint 

  REAL(wp)                            :: psfc     ! Sfc pressure
  REAL(wp)                            :: psfc_ad  ! Sfc pressure perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: p_hlv    ! Pressure on half levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: p_hlv_ad ! Pressure perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv    ! Virtual temperauture (K)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv_ad ! Temperauture perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: akk      ! a coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bkk      ! b coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_hlv ! Geopot height half level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_hlv_ad ! Geopot ht pertubation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ln_prflv ! log ratio of pressures
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ln_prflv_ad ! log ratio perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_p    ! Pressure difference (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_p_ad ! Pressure diff perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_geop ! Thickness (m)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_geop_ad ! Thickness perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha    ! Interpolation coefficient
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha_ad ! Coefficient perturbation
  INTEGER                             :: lvl
  INTEGER                             :: n_hlv
  INTEGER                             :: n_flv

!-------------------------------------------------------------------------------
! 2. Reset local adjoint variables
!-------------------------------------------------------------------------------

  n_hlv = SIZE(x%ak)
  n_flv = SIZE(x%pres)

  ALLOCATE(p_hlv(n_hlv))
  ALLOCATE(p_hlv_ad(n_hlv))
  ALLOCATE(Tvflv(SIZE(x%temp)))
  ALLOCATE(Tvflv_ad(SIZE(x%temp)))
  ALLOCATE(akk(n_flv))
  ALLOCATE(bkk(n_flv))
  ALLOCATE(geop_hlv(n_hlv))
  ALLOCATE(geop_hlv_ad(n_hlv))
  ALLOCATE(ln_prflv(n_flv))
  ALLOCATE(ln_prflv_ad(n_flv))
  ALLOCATE(del_p(n_flv))
  ALLOCATE(del_p_ad(n_flv))
  ALLOCATE(del_geop(n_flv))
  ALLOCATE(del_geop_ad(n_flv))
  ALLOCATE(alpha(n_flv))
  ALLOCATE(alpha_ad(n_flv))

  psfc_ad = 0.0_wp

  alpha_ad(:)    = 0.0_wp
  del_p_ad(:)    = 0.0_wp
  del_geop_ad(:)  = 0.0_wp
  ln_prflv_ad(:) = 0.0_wp
  geop_hlv_ad(:)  = 0.0_wp
  p_hlv_ad(:) = 0.0_wp
  Tvflv_ad(:) = 0.0_wp

!-------------------------------------------------------------------------------
! 3. Useful constants
!-------------------------------------------------------------------------------

  n_hlv = SIZE(p_hlv)
  
!-------------------------------------------------------------------------------
! 4. Update state variables
!-------------------------------------------------------------------------------

  IF (x%use_logp) THEN
     psfc = EXP(x%state(2*x%n_lev + 1)) * 100.0_wp
  ELSE
     psfc = x%state(2*x%n_lev + 1)
  ENDIF

!-------------------------------------------------------------------------------
! 5. Recompute intermediate variables
!-------------------------------------------------------------------------------

  bkk = 0.5_wp * ( x%bk(2:n_hlv) + x%bk(1:n_hlv-1) )
  p_hlv = x%ak + x%bk * psfc
  Tvflv = (1.0_wp + ((R_vap/R_dry)-1.d0) * x%shum) * x%temp
  
  del_p    = p_hlv(1:n_hlv-1) - p_hlv(2:n_hlv)
  ln_prflv = LOG(p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv))
  alpha    = 1.0_wp - p_hlv(2:n_hlv)/del_p*ln_prflv
  alpha(n_hlv-1) = LOG(2.0_wp)

!-------------------------------------------------------------------------------
! 6. Adjoint of computing geopotential height on full levels
!-------------------------------------------------------------------------------

! 6.1 Adjoint of interpolation onto full levels

  Tvflv_ad = Tvflv_ad + x_ad%geop * alpha * R_dry / g_wmo
  alpha_ad = alpha_ad + x_ad%geop * R_dry * Tvflv / g_wmo
  geop_hlv_ad(1:n_hlv-1) = geop_hlv_ad(1:n_hlv-1) + x_ad%geop
  x_ad%geop  = 0.0_wp

! 6.2 Adjoint of geopotential height integral

  DO lvl = 1, n_hlv
     del_geop_ad(1:lvl-1) = del_geop_ad(1:lvl-1) + geop_hlv_ad(lvl)
     geop_hlv_ad(lvl) = 0.0_wp
  END DO

! 6.3 Adjoint of function to integrated

  Tvflv_ad   = Tvflv_ad   + del_geop_ad * R_dry * ln_prflv / g_wmo
  ln_prflv_ad = ln_prflv_ad + del_geop_ad * R_dry * Tvflv / g_wmo
  del_geop_ad  = 0.0_wp

! 6.4 Adjoint of interpolation coefficients

  alpha_ad(n_hlv-1)   = 0.0_wp
  del_p_ad            = del_p_ad + alpha_ad * p_hlv(2:n_hlv)/(del_p*del_p) & 
                          * ln_prflv
  ln_prflv_ad         = ln_prflv_ad - alpha_ad * (p_hlv(2:n_hlv)/del_p)
  p_hlv_ad(1:n_hlv-1) = p_hlv_ad(1:n_hlv-1) - alpha_ad/del_p*ln_prflv
  alpha_ad            = 0.0_wp

! 6.5 Adjoint of log of pressure ratio

  p_hlv_ad(1:n_hlv-1) = p_hlv_ad(1:n_hlv-1) &
                        + ln_prflv_ad &
                        * (1./(p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv))/p_hlv(2:n_hlv))
  p_hlv_ad(2:n_hlv) = p_hlv_ad(2:n_hlv) &
                       - ln_prflv_ad &
                       * 1. / (p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv)) &
                       * (p_hlv(1:n_hlv-1) / (p_hlv(2:n_hlv) * p_hlv(2:n_hlv)))
  ln_prflv_ad         = 0.0_wp

! 6.6 Adjoint of pressure differences

  p_hlv_ad(1:n_hlv-1) = p_hlv_ad(1:n_hlv-1) + del_p_ad
  p_hlv_ad(2:n_hlv) = p_hlv_ad(2:n_hlv) - del_p_ad
  del_p_ad = 0.0_wp

! 6.7 Adjoint of virtual temperature calculation

  x_ad%shum  = x_ad%shum + 0.61_wp * Tvflv_ad * x%temp
  x_ad%temp  = x_ad%temp + Tvflv_ad * (1.0_wp + ((R_vap/R_dry)-1.d0) * x%shum)
  Tvflv_ad = 0.0_wp

! 6.8 Adjoint of calculating pressure levels

  psfc_ad = psfc_ad + SUM(x_ad%pres*bkk)
  x_ad%pres = 0.0_wp

! 6.9 Adjoint of calculating half pressure levels

  psfc_ad = psfc_ad + SUM(p_hlv_ad*x%bk)
  p_hlv_ad = 0.0_wp

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

  IF (x%use_logp) THEN
     x_ad%state(2*x%n_lev+1) = x_ad%state(2*x%n_lev+1) +          &
                                psfc_ad*EXP(x%state(2*x%n_lev+1)) * 100.0_wp
  ELSE
     x_ad%state(2*x%n_lev+1) = x_ad%state(2*x%n_lev+1) + psfc_ad
  ENDIF
  psfc_ad = 0.0_wp

  IF (x%use_logq) THEN
     x_ad%state(x%n_lev+1:2*x%n_lev) = x_ad%state(x%n_lev+1:2*x%n_lev) +     &
                                         x_ad%shum *                         &
                                         EXP(x%state(x%n_lev+1:2*x%n_lev)) / &
                                         1000.0_wp
   ELSE
     WHERE (x%shum <= 0.0_wp)
        x_ad%shum = 0.0_wp
     endwhere
     x_ad%state(x%n_lev+1:2*x%n_lev)=x_ad%state(x%n_lev+1:2*x%n_lev)+x_ad%shum
  ENDIF
  x_ad%shum = 0.0_wp

  x_ad%state(1:x%n_lev) = x_ad%state(1:x%n_lev) + x_ad%temp
  x_ad%temp = 0.0_wp
  
!-------------------------------------------------------------------------------
! 8. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(p_hlv)
  DEALLOCATE(p_hlv_ad)
  DEALLOCATE(Tvflv)
  DEALLOCATE(Tvflv_ad)
  DEALLOCATE(akk)
  DEALLOCATE(bkk)
  DEALLOCATE(geop_hlv)
  DEALLOCATE(geop_hlv_ad)
  DEALLOCATE(ln_prflv)
  DEALLOCATE(ln_prflv_ad)
  DEALLOCATE(del_p)
  DEALLOCATE(del_p_ad)
  DEALLOCATE(del_geop)
  DEALLOCATE(del_geop_ad)
  DEALLOCATE(alpha)
  DEALLOCATE(alpha_ad)


END SUBROUTINE ropp_fm_state2state_ecmwf_ad
