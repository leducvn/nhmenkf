! $Id: ropp_fm_state2state_ecmwf_tl.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_ecmwf_tl(x, x_tl)

!****s* Model_ecmwf/ropp_fm_state2state_ecmwf_tl *
!
! NAME
!    ropp_fm_state2state_ecmwf_tl - Tangent linear of ropp_fm_state2state_ecmwf
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_ecmwf_tl(x,x_tl)
! 
! DESCRIPTION
!    This routine is the tangent linear of ropp_fm_state2state_ecmwf
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
!   ropp_fm_state2state_ecmwf
!   ropp_fm_state2state_ecmwf_ad
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
  USE ropp_fm_types
  USE ropp_fm_constants

  IMPLICIT NONE

  TYPE(State1dFM),    INTENT(in)      :: x        ! State  
  TYPE(State1dFM),    INTENT(inout)   :: x_tl     ! Perturbation

  REAL(wp)                            :: psfc     ! Surface pressure
  REAL(wp)                            :: psfc_tl  ! psfc perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: p_hlv    ! Pressure on half levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: p_hlv_tl ! Pressure perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv    ! Virtual temperauture (K)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv_tl ! Temperauture perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: akk      ! a coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bkk      ! b coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_hlv ! Geopot height half level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_hlv_tl ! Geopot ht pertubation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ln_prflv ! log ratio of pressures
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ln_prflv_tl ! log ratio perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_p    ! Pressure difference (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_p_tl ! Pressure diff perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_geop ! Thickness (m)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_geop_tl ! Thickness perturbation
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha    ! Interpolation coefficient
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha_tl ! Coefficient perturbation
  INTEGER                             :: lvl      ! Index
  INTEGER                             :: n_hlv    ! Number of half levels
  INTEGER                             :: n_flv    ! Number of full levels

!-------------------------------------------------------------------------------
! 2. Useful constants
!-------------------------------------------------------------------------------

  n_hlv = SIZE(x%ak)
  n_flv = SIZE(x%pres)

  ALLOCATE(p_hlv(n_hlv))
  ALLOCATE(p_hlv_tl(n_hlv))
  ALLOCATE(Tvflv(SIZE(x%temp)))
  ALLOCATE(Tvflv_tl(SIZE(x%temp)))
  ALLOCATE(akk(n_flv))
  ALLOCATE(bkk(n_flv))
  ALLOCATE(geop_hlv(n_hlv))
  ALLOCATE(geop_hlv_tl(n_hlv))
  ALLOCATE(ln_prflv(n_flv))
  ALLOCATE(ln_prflv_tl(n_flv))
  ALLOCATE(del_p(n_flv))
  ALLOCATE(del_p_tl(n_flv))
  ALLOCATE(del_geop(n_flv))
  ALLOCATE(del_geop_tl(n_flv))
  ALLOCATE(alpha(n_flv))
  ALLOCATE(alpha_tl(n_flv))

!-------------------------------------------------------------------------------
! 3. Update state variables
!-------------------------------------------------------------------------------

  x_tl%temp = x_tl%state(1:x%n_lev)
  
  IF (x%use_logp) THEN
     psfc = EXP(x%state(2*x%n_lev+1)) * 100.0_wp
     psfc_tl = x_tl%state(2*x%n_lev+1) * EXP(x%state(2*x%n_lev+1)) * 100.0_wp
  ELSE
     psfc = x%state(2*x%n_lev+1)
     psfc_tl = x_tl%state(2*x%n_lev+1)
  ENDIF
  
  IF (x%use_logq) THEN
     x_tl%shum = x_tl%state(x%n_lev+1:2*x%n_lev) *     &
                   EXP(x%state(x%n_lev+1:2*x%n_lev)) / 1000.0_wp
  ELSE
     x_tl%shum = x_tl%state(x%n_lev+1:2*x%n_lev)
     WHERE(x%shum <= 0.0_wp)
        x_tl%shum = 0.0_wp
     endwhere
  ENDIF

  IF (x%direct_ion) THEN
    x_tl%ne_max  = x_tl%state(2*x%n_lev+2)
    x_tl%h_peak  = x_tl%state(2*x%n_lev+3)
    x_tl%h_width = x_tl%state(2*x%n_lev+4)
  ENDIF

!-------------------------------------------------------------------------------
! 4. Calculate half pressure levels
!-------------------------------------------------------------------------------

  p_hlv    = x%ak + x%bk * psfc
  p_hlv_tl = x%bk * psfc_tl 

!-------------------------------------------------------------------------------
! 5. Calculate pressure on full levels
!-------------------------------------------------------------------------------

  akk = 0.5_wp * ( x%ak(2:n_hlv) + x%ak(1:n_hlv-1) )
  bkk = 0.5_wp * ( x%bk(2:n_hlv) + x%bk(1:n_hlv-1) )
!  x%pres = akk + bkk * psfc
  x_tl%pres = bkk * psfc_tl

!-------------------------------------------------------------------------------
! 6. Calculate virtual temperature on full levels
!-------------------------------------------------------------------------------

  Tvflv    = (1.0_wp + ((R_vap/R_dry)-1.d0) * x%shum) * x%temp
  Tvflv_tl = Tvflv * x_tl%temp / x%temp + ((R_vap/R_dry)-1.d0)*x%temp*x_tl%shum

!-------------------------------------------------------------------------------
! 7. Calculate geopotential height on full levels (assume index 1 towards sfc)
!-------------------------------------------------------------------------------

! 7.1 Pressure differences

  del_p    = p_hlv(1:n_hlv-1) - p_hlv(2:n_hlv)
  del_p_tl = p_hlv_tl(1:n_hlv-1) - p_hlv_tl(2:n_hlv)

! 7.2 Log of pressure ratio

  ln_prflv    = LOG( p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv) )
  ln_prflv_tl = p_hlv_tl(1:n_hlv-1)/p_hlv(1:n_hlv-1) &
                  - p_hlv_tl(2:n_hlv)/p_hlv(2:n_hlv)

! 7.3 Interpolation coefficients

  alpha    = 1.0_wp - p_hlv(2:n_hlv)/del_p*ln_prflv
  alpha_tl = del_p_tl*p_hlv(2:n_hlv)/(del_p*del_p)*ln_prflv  &
              - ln_prflv_tl*(p_hlv(2:n_hlv)/del_p) &
                  - p_hlv_tl(2:n_hlv)/del_p*ln_prflv
  alpha(n_hlv-1)    = LOG(2.0_wp)
  alpha_tl(n_hlv-1) = 0.0_wp

! 7.4 Function to be integrated

  del_geop_tl = Tvflv_tl*R_dry*ln_prflv/g_wmo + ln_prflv_tl*R_dry*Tvflv/g_wmo

! 7.5 Calculate geopotential height integral 

  DO lvl = 1, n_hlv
     geop_hlv_tl(lvl) = SUM(del_geop_tl(1:lvl-1))
  END DO

! 7.6 Interpolate onto full levels 

  x_tl%geop = geop_hlv_tl(1:n_hlv-1) + Tvflv_tl*alpha*R_dry/g_wmo &
               + alpha_tl*R_dry*Tvflv/g_wmo 

!-------------------------------------------------------------------------------
! 8. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(p_hlv)
  DEALLOCATE(p_hlv_tl)
  DEALLOCATE(Tvflv)
  DEALLOCATE(Tvflv_tl)
  DEALLOCATE(akk)
  DEALLOCATE(bkk)
  DEALLOCATE(geop_hlv)
  DEALLOCATE(geop_hlv_tl)
  DEALLOCATE(ln_prflv)
  DEALLOCATE(ln_prflv_tl)
  DEALLOCATE(del_p)
  DEALLOCATE(del_p_tl)
  DEALLOCATE(del_geop)
  DEALLOCATE(del_geop_tl)
  DEALLOCATE(alpha)
  DEALLOCATE(alpha_tl)
        
END SUBROUTINE ropp_fm_state2state_ecmwf_tl
