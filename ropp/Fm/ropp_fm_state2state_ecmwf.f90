! $Id: ropp_fm_state2state_ecmwf.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_state2state_ecmwf_1d(x)

!****s* Model_ecmwf/ropp_fm_state2state_ecmwf *
!
! NAME
!    ropp_fm_state2state_ecmwf - Calculate pressure and geopotential height at 
!                                full levels for ECMWF hybrid-sigma vertical 
!                                level model type and update state vector
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_ecmwf(x)
! 
! DESCRIPTION
!    This subroutine calculates geopotential height and pressure at full levels
!    for a hybrid vertical coordinate.
!
! INPUTS
!    type(State1dFM)         :: x   State vector data structure
!
! OUTPUT
!    type(State1dFM)         :: x   State vector data structure
!
! NOTES
!    This function requires and returns SI units (i.e., Pa) for level 
!    coefficients, surface and level pressure.
!    Calculation of geopotential height assumes data arrays increasing with 
!    height (index 1 closest to surface).
!
! SEE ALSO
!   ropp_fm_state2state_ecmwf_ad
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

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_types, ONLY: State1dFM
  USE ropp_fm_constants

  IMPLICIT NONE

  TYPE(State1dFM),   INTENT(inout)    :: x        ! State vector structure

  REAL(wp)                            :: psfc     ! Surface pressure (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: p_hlv    ! Pressure on half levels
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tvflv    ! Virtual temperauture (K)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: akk      ! a coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bkk      ! b coefficient on full level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geop_hlv ! Geopot height half level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ln_prflv ! log ratio of pressures
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_p    ! Pressure difference (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: del_geop ! Thickness (m)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha    ! Interpolation coefficient
  INTEGER                             :: lvl      ! Index
  INTEGER                             :: n_hlv    ! Number of half levels
  INTEGER                             :: n_flv    ! Number of full levels

!-------------------------------------------------------------------------------
! 2. Useful constants
!-------------------------------------------------------------------------------

  n_hlv = SIZE(x%ak)
  n_flv = SIZE(x%pres)

  ALLOCATE(p_hlv(n_hlv))
  ALLOCATE(Tvflv(n_flv))
  ALLOCATE(akk(n_flv))
  ALLOCATE(bkk(n_flv))
  ALLOCATE(geop_hlv(n_hlv))
  ALLOCATE(ln_prflv(n_flv))
  ALLOCATE(del_p(n_flv))
  ALLOCATE(del_geop(n_flv))
  ALLOCATE(alpha(n_flv))

!-------------------------------------------------------------------------------
! 3. Update state variables
!-------------------------------------------------------------------------------

  x%temp = x%state(1:x%n_lev)

  IF (x%use_logq) THEN
    x%shum = EXP(x%state(x%n_lev + 1: 2*x%n_lev)) / 1000.0_wp
    IF (x%check_qsat) CALL check_qsat(x%shum, x%temp, x%pres)
  ELSE
    x%shum = x%state(x%n_lev + 1: 2*x%n_lev)
    IF (x%check_qsat) CALL check_qsat(x%shum, x%temp, x%pres)
    WHERE (x%shum <= 0.0_wp)
      x%shum = 1.0e-9_wp
    ENDWHERE
  ENDIF
  
  IF (x%use_logp) THEN
     psfc = EXP(x%state(2*x%n_lev + 1)) * 100.0_wp
  ELSE
     psfc = x%state(2*x%n_lev + 1)
  ENDIF
  x%ak(n_hlv) = 1.0e-32_wp

  IF (x%direct_ion) THEN
    x%ne_max  = x%state(2*x%n_lev+2)
    x%h_peak  = x%state(2*x%n_lev+3)
    x%h_width = x%state(2*x%n_lev+4)
  ENDIF

!-------------------------------------------------------------------------------
! 4. Calculate half pressure levels
!-------------------------------------------------------------------------------

  p_hlv = x%ak + x%bk * psfc

!-------------------------------------------------------------------------------
! 5. Calculate pressure on full levels
!-------------------------------------------------------------------------------

  akk = 0.5_wp * ( x%ak(1:n_hlv-1) + x%ak(2:n_hlv))
  bkk = 0.5_wp * ( x%bk(1:n_hlv-1) + x%bk(2:n_hlv))
  x%pres = akk + bkk*psfc

!-------------------------------------------------------------------------------
! 6. Calculate virtual temperature on full levels
!-------------------------------------------------------------------------------

  Tvflv  = (1.0_wp + 0.61_wp * x%shum) * x%temp

!-------------------------------------------------------------------------------
! 7. Calculate geopotential height on full levels (assume index 1 towards sfc)
!-------------------------------------------------------------------------------

! 7.1 Pressure differences 

  del_p = p_hlv(1:n_hlv-1) - p_hlv(2:n_hlv)

! 7.2 Log of pressure ratio

  ln_prflv = LOG(p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv))

! 7.3 Interpolation coefficients 

  alpha    = 1.0_wp - p_hlv(2:n_hlv)/del_p * ln_prflv
  alpha(n_hlv-1) = LOG(2.0_wp)

! 7.4 Function to be integrated 

  del_geop =   R_dry * Tvflv * ln_prflv / g_wmo

! 7.5 Calculate geopotential height integral

  DO lvl = 1, n_hlv

     geop_hlv(lvl) = x%geop_sfc + SUM(del_geop(1:lvl-1)) !! N.B:sum(X(1:0))=0

  ENDDO

! 7.6 Interpolate onto full levels 

  x%geop = geop_hlv(1:n_hlv-1) + alpha * R_dry * Tvflv / g_wmo
  
!-------------------------------------------------------------------------------
! 8. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(p_hlv)
  DEALLOCATE(Tvflv)
  DEALLOCATE(akk)
  DEALLOCATE(bkk)
  DEALLOCATE(geop_hlv)
  DEALLOCATE(ln_prflv)
  DEALLOCATE(del_p)
  DEALLOCATE(del_geop)
  DEALLOCATE(alpha)

CONTAINS
  
!-------------------------------------------------------------------------------
! 9. Compute Saturation Vapour Pressure - limit log(humidity) within saturation 
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

! 9.1 Goff-Gratch equation (saturation vapour pressure over water)
    
    tfrac(:) = 373.16_wp / temp(:)
    Z = (-7.90298_wp * (tfrac - 1.0_wp)) +  &
       (5.02808_wp * LOG10(tfrac)) - (1.3816e-7_wp * (10**(11.344*(1.0_wp - (1.0_wp/tfrac))) - 1.0_wp)) + &
       (8.1328e-3_wp * (10**(-3.49149*(tfrac-1.0_wp)) - 1.0_wp))
    
! 9.2 Goff-Gratch equation (saturation vapour pressure over ice)

    tfrac(:) = 273.16_wp / temp(:)
    X = (-9.09718_wp * (tfrac - 1.0_wp)) + (-3.56654_wp*LOG10(tfrac)) + (0.876793_wp*(1.0_wp - (1.0_wp/tfrac)))

    WHERE(temp >= 273.16)
      es = 1013246.0_wp * (10**Z)
    ELSEWHERE
      es = 610.71_wp * (10**X)
    ENDWHERE
    
! 9.3 Compute saturation humidity

    qs = es * epsilon_water / (MAX(press,es) - ((1.0_wp - epsilon_water)*es))

    WHERE(qs == 1.0_wp)
      qs = 2.0_wp*shum(SIZE(shum))
    ENDWHERE

! 9.4 Limit humidity data below saturation

    WHERE(shum > qs)
      shum = qs
    ENDWHERE

    DEALLOCATE(es)
    DEALLOCATE(qs)
    DEALLOCATE(z)
    DEALLOCATE(x)
    DEALLOCATE(tfrac)

  END SUBROUTINE check_qsat

END SUBROUTINE ropp_fm_state2state_ecmwf_1d


!****s* Model_ecmwf/ropp_fm_state2state_ecmwf_2d *
!
! NAME
!    ropp_fm_state2state_ecmwf_2d - Calculate pressure and geopotential height
!                                full levels for ECMWF hybrid-sigma vertical 
!                                level model type and update state vector
!                                TWO-DIMENSIONAL BACKGROUND VERSION
!
! SYNOPSIS
!    use ropp_fm
!      ...
!    call ropp_fm_state2state_ecmwf_2d(x)
! 
! DESCRIPTION
!    This subroutine calculates geopotential height and pressure at full levels
!    for a hybrid vertical coordinate.
!
! INPUTS
!    type(State2dFM)         :: x   State vector data structure
!
! OUTPUT
!    type(State2dFM)         :: x   State vector data structure
!
! NOTES
!    This function requires and returns SI units (i.e., Pa) for level 
!    coefficients, surface and level pressure.
!    Calculation of geopotential height assumes data arrays increasing with 
!    height (index 1 closest to surface).
!
! SEE ALSO
!   ropp_fm_state2state_ecmwf_ad
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


SUBROUTINE ropp_fm_state2state_ecmwf_2d(x)

  ! 1. Declarations
  ! ---------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_types, ONLY: State2dFM
  USE ropp_fm_constants

  IMPLICIT NONE

  TYPE(State2dFM),   INTENT(inout)  :: x

  REAL(wp), DIMENSION(x%n_horiz)              :: psfc     ! surface pressure
  REAL(wp), DIMENSION(SIZE(x%ak),x%n_horiz)   :: p_hlv    ! pressure on half levels
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz) :: Tvflv    ! virtual temperauture
  REAL(wp), DIMENSION(x%n_lev) :: akk      ! a coefficient on full levels
  REAL(wp), DIMENSION(x%n_lev) :: bkk      ! b coefficient on full levels
  INTEGER                           :: lvl  
  REAL(wp), DIMENSION(SIZE(x%ak),x%n_horiz)   :: geop_hlv ! geopotential height half lvl
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz) :: ln_prflv ! log ratio of pressures
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz) :: del_p    ! pressure difference
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz) :: del_geop ! thickness
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz) :: alpha
  INTEGER                           :: n_hlv
  INTEGER                           :: i

  ! 1. Useful constants
  ! -------------------

  n_hlv = SIZE(x%ak)

  ! 2. Update state variables
  ! -------------------------

  x%ak(n_hlv) = 1.0e-32_wp
  psfc(:) = x%pres_sfc(:)

  ! 3. Level coefficinets on full levels
  ! ---------------------------------

  akk = 0.5_wp * ( x%ak(1:n_hlv-1) + x%ak(2:n_hlv))
  bkk = 0.5_wp * ( x%bk(1:n_hlv-1) + x%bk(2:n_hlv))

  ! 4. Surface pressure and pressure on full levels
  !----------------------------------

  DO i = 1, x%n_horiz

   p_hlv(:,i) = x%ak(:) + x%bk(:) * psfc(i)
   x%pres(:,i) = akk(:) + bkk(:)*psfc(i)

  ENDDO

  ! 5. Calculate virtual temperature on full levels
  ! -----------------------------------------------

  Tvflv(:,:)  = (1.0_wp + 0.61_wp * x%shum(:,:)) * x%temp(:,:)

  ! 7. Calculate geopotential height on full levels (assume index 1 towards sfc)
  ! -----------------------------------------------

  ! 7.1 Pressure differences 

  del_p(:,:) = p_hlv(1:n_hlv-1,:) - p_hlv(2:n_hlv,:)

  ! 7.2 Log of pressure ratio

  ln_prflv(:,:) = LOG(p_hlv(1:n_hlv-1,:)/p_hlv(2:n_hlv,:))

  ! 7.3 Interpolation coefficients 

  alpha(:,:)   = 1.0_wp - p_hlv(2:n_hlv,:)/del_p(:,:) * ln_prflv(:,:)
  alpha(n_hlv-1,:) = LOG(2.0_wp)

  ! 7.4 Function to be integrated 

  del_geop(:,:) =   R_dry * Tvflv(:,:) * ln_prflv(:,:) / g_wmo

  ! 7.5 Calculate geopotential height integral 

  ! lowest half level is surface elevation

  geop_hlv(1,:) = x%geop_sfc(:)

  ! height of other half levels

  DO lvl = 2, n_hlv

   geop_hlv(lvl,:) = geop_hlv(lvl-1,:) + del_geop(lvl-1,:) ! 1D formulation doesn't work 

  ENDDO

  ! 7.6 Interpolate onto full levels 

  x%geop(:,:) = geop_hlv(1:n_hlv-1,:) + alpha(:,:) * R_dry * Tvflv(:,:) / g_wmo


END SUBROUTINE ropp_fm_state2state_ecmwf_2d
