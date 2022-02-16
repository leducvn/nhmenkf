! $Id: ropp_fm_bangle_1d_tl.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_bangle_1d_tl(x, x_tl, y, y_tl)

!****s* BendingAngle/ropp_fm_bangle_1d_tl *
!
! NAME
!    ropp_fm_bangle_1d_tl - Tangent linear of ropp_fm_bangle_1d().
!
! SYNOPSIS
!    call ropp_fm_bangle_1d_tl(x, x_tl, y, y_tl)
! 
! DESCRIPTION
!    This routine is the tangent linear of ropp_fm_bangle_1d.
!
! INPUTS
!    type(State1dFM)   :: x      ! State vector structure
!    type(State1dFM)   :: x_tl   ! Perturbation vector structure
!    type(Obs1dBangle) :: y      ! Bending angle observation vector
!
! OUTPUT
!    real(wp), dimension(:) :: y_tl  ! Bending angle perturbation
!
! NOTES
!    The obs vector is required only for the observation's impact parameter
!    levels; no forward simulated bending angle profile is returned.
!
!    The lengths of the arrays x_tl%state and y_tl must agree with the 
!    lengths of the x%state and y%bangle arrays, respectively.
!
! SEE ALSO
!    ropp_fm_types
!    ropp_fm_bangle_1d
!    ropp_fm_bangle_1d_ad
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
  USE ropp_utils
  USE ropp_fm,   not_this => ropp_fm_bangle_1d_tl
! USE ropp_fm
  
! USE ropp_fm_types
! USE ropp_fm_constants

  IMPLICIT NONE

  TYPE(State1dFM),   INTENT(in)                    :: x                ! State vector
  TYPE(State1dFM),   INTENT(in)                    :: x_tl             ! State perturbation
  TYPE(Obs1dBangle), INTENT(in)                    :: y                ! Obs vector
  REAL(wp), DIMENSION(SIZE(y%bangle)), INTENT(out) :: y_tl             ! Obs perturbation

  REAL(wp), DIMENSION(x%n_lev)                     :: pwvp             ! Partial water vapour pressure
  REAL(wp), DIMENSION(x%n_lev)                     :: pwvp_tl          ! Pwvp perturbation
  REAL(wp), DIMENSION(x%n_lev)                     :: pdry             ! Dry pressure
  REAL(wp), DIMENSION(x%n_lev)                     :: pdry_tl          ! Pdry perturbation
  REAL(wp), DIMENSION(x%n_lev)                     :: refrac           ! Refractivity on bg model levels
  REAL(wp), DIMENSION(x%n_lev)                     :: refrac_tl        ! Refractivity perturbation

  REAL(wp), DIMENSION(x%n_lev)                     :: z_geop           ! Geopotential height of model levels
  REAL(wp), DIMENSION(x%n_lev)                     :: z_geop_tl        ! GPH perturbation
  REAL(wp), DIMENSION(x%n_lev)                     :: zcomp_dry_inv    ! Dry compressibility
  REAL(wp), DIMENSION(x%n_lev)                     :: zcomp_dry_inv_tl ! Dry compressibility perturbation
  REAL(wp), DIMENSION(x%n_lev)                     :: zcomp_wet_inv    ! Wet compressibility
  REAL(wp), DIMENSION(x%n_lev)                     :: zcomp_wet_inv_tl ! Wet compressibility perturbation

  REAL(wp), DIMENSION(x%n_lev)                     :: h                ! Geometric height
  REAL(wp), DIMENSION(x%n_lev)                     :: h_tl             ! Geometric height perturbation
  REAL(wp), DIMENSION(x%n_lev)                     :: impact           ! Impact parameter
  REAL(wp), DIMENSION(x%n_lev)                     :: impact_tl        ! Impact perturbation
  REAL(wp), DIMENSION(SIZE(y%bangle))              :: bangle           ! Bending angle

  REAL(wp)                                         :: kap1,kap2,kap3   ! Refractivity coefficients used in routine
  REAL(wp)                                         :: R_peak           ! Radius of Chapman layer peak
  REAL(wp)                                         :: R_peak_tl        ! Radius of Chapman layer peak perturbation

!-------------------------------------------------------------------------------
! 2. Non ideal gas options 
!-------------------------------------------------------------------------------

! set inverse of compressibilities

  zcomp_dry_inv(:) = 1.0_wp
  zcomp_wet_inv(:) = 1.0_wp
 
  zcomp_dry_inv_tl(:) = 0.0_wp
  zcomp_wet_inv_tl(:) = 0.0_wp
    
! initialise geopotential heights
   
  z_geop(:) = x%geop(:)
  z_geop_tl(:) = x_tl%geop(:)

  IF (x%non_ideal) THEN

! if non ideal gas calculation, use adjusted coefficients

     kap1 = kappa1_comp
     kap2 = kappa2_comp
     kap3 = kappa3_comp   

!    calculate compressibilty and adjust geopotential heights in z_geop

     CALL ropp_fm_compress_tl &
     &(x,x_tl,z_geop,z_geop_tl,zcomp_dry_inv,zcomp_dry_inv_tl,&
     &zcomp_wet_inv,zcomp_wet_inv_tl)  

  ELSE

     kap1 = kappa1
     kap2 = kappa2
     kap3 = kappa3
     
  ENDIF  


!-------------------------------------------------------------------------------
! 2. Calculate water vapor  and dry pressure pressure
!-------------------------------------------------------------------------------
  
  pwvp = x%pres * x%shum / (epsilon_water + (1.0_wp - epsilon_water)*x%shum)
  
  pwvp_tl = pwvp * ( x_tl%pres/x%pres + x_tl%shum/x%shum &
                - x_tl%shum*pwvp*(1.0_wp - epsilon_water) / (x%pres*x%shum))

! dry pressure

  pdry = x%pres - pwvp
  pdry_tl =  x_tl%pres - pwvp_tl

!-------------------------------------------------------------------------------
! 3. Calculate refractivity
!-------------------------------------------------------------------------------

  refrac = kap1 * pdry * zcomp_dry_inv / x%temp    + &
           kap2 * pwvp * zcomp_wet_inv / x%temp**2 + &
           kap3 * pwvp * zcomp_wet_inv / x%temp

  refrac_tl = kap1 * pdry_tl * zcomp_dry_inv/ x%temp +  &
           kap2 * pwvp_tl * zcomp_wet_inv/ x%temp**2 + &
           kap3 * pwvp_tl * zcomp_wet_inv/ x%temp + &
           kap1 * pdry * zcomp_dry_inv_tl/ x%temp + &
           kap2 * pwvp * zcomp_wet_inv_tl/ x%temp**2 + &
           kap3 * pwvp * zcomp_wet_inv_tl/ x%temp - &
           (kap1 * pdry * zcomp_dry_inv/ x%temp**2 + &
           2.0_wp *kap2 * pwvp * zcomp_wet_inv/ x%temp**3 + &
           kap3 * pwvp * zcomp_wet_inv/ x%temp**2)* x_tl%temp

!-------------------------------------------------------------------------------
! 4. Calculate geometric height
!-------------------------------------------------------------------------------

  h    = y%r_earth * z_geop / (y%g_sfc / g_wmo * y%r_earth - z_geop)
  
  h_tl = z_geop_tl * (y%r_earth/(((y%g_sfc/g_wmo)*y%r_earth)-z_geop) &
               + y%r_earth*z_geop/((((y%g_sfc/g_wmo)*y%r_earth)-z_geop)**2))

!-------------------------------------------------------------------------------
! 5. Calculate impact parameter
!-------------------------------------------------------------------------------

  IF (y%undulation > ropp_MDTV) THEN
    impact    = (1.0_wp + 1.0e-6_wp * refrac) * (h + y%r_curve + y%undulation)
    impact_tl = 1.0e-6_wp * refrac_tl * (h + y%r_curve + y%undulation) &
                  + h_tl * (1.0_wp + 1.0e-6_wp*refrac)
  ELSE
    CALL message(msg_warn, "Undulation missing. " // &
                 "Will assume to be zero when calculating full and " // &
                 "perturbed impact parameters.")
    impact    = (1.0_wp + 1.0e-6_wp * refrac) * (h + y%r_curve)
    impact_tl = 1.0e-6_wp * refrac_tl * (h + y%r_curve) &
                  + h_tl * (1.0_wp + 1.0e-6_wp*refrac)
  END IF
  
!-------------------------------------------------------------------------------
! 6. Calculate neutral bending angle
!-------------------------------------------------------------------------------

  CALL ropp_fm_abel_tl( &
     impact, refrac, x%temp, x_tl%temp, y%r_curve, x%new_bangle_op, &
     y%impact, impact_tl, refrac_tl, y_tl)

!-------------------------------------------------------------------------------
! 7. Calculate the ionospheric bending if L1 and L2 are used in retrieval
!-------------------------------------------------------------------------------

  IF (x%direct_ion) THEN

! Need neutral bending angle for a test in ropp_fm_iono_bangle_tl

    CALL ropp_fm_abel(impact, refrac, x%temp, y%r_curve, &
                      x%new_bangle_op, y%impact, bangle)

    R_peak = x%H_peak + y%r_curve

    R_peak_tl = x_tl%H_peak

    CALL ropp_fm_iono_bangle_tl(x%Ne_max,    R_peak,    x%H_width, &
                                x_tl%Ne_max, R_peak_tl, x_tl%H_width, &
                                y%n_L1, y%impact, bangle, y_tl)

  END IF


END SUBROUTINE ropp_fm_bangle_1d_tl
