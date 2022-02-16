! $Id: ropp_fm_bangle_1d_ad.f90 4010 2014-01-10 11:07:40Z idculv $

SUBROUTINE ropp_fm_bangle_1d_ad(x, x_ad, y, y_ad)

!****s* BendingAngle/ropp_fm_bangle_1d_ad *
!
! NAME
!    ropp_fm_bangle_1d_ad - Adjoint of ropp_fm_bangle_1d().
!
! SYNOPSIS
!    call ropp_fm_bangle_1d_ad(x, x_ad, y, y_ad)
! 
! DESCRIPTION
!    This routine is the adjoint of ropp_fm_bangle_1d.
!
! INPUTS
!    type(State1dFM)        :: x        ! State vector
!    type(Obs1dBangle)      :: y        ! Observation vector
!    real(wp), dimension(:) :: y_ad     ! Adjoint forcing
!
! OUTPUT
!    type(State1dFM)        :: x_ad     ! State vector adjoint
!
! NOTES
!    The obs vector is required only for the observation's geopotential levels;
!    no forward simulated refractivity profile is returned.
!
!    The lengths of the arrays state_ad%state and obs_ad must agree with the 
!    lengths of the state%state and obs%obs arrays, respectively.
!
! SEE ALSO
!    ropp_fm_types
!    ropp_fm_bangle_1d
!    ropp_fm_bangle_1d_tl
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
  USE ropp_fm,   not_this => ropp_fm_bangle_1d_ad
! USE ropp_fm
! USE ropp_fm_types
! USE ropp_fm_iono
! USE ropp_fm_constants
  USE geodesy

  IMPLICIT NONE

  TYPE(State1dFM),        INTENT(in)    :: x                ! State vector
  TYPE(State1dFM),        INTENT(inout) :: x_ad             ! State vector adjoint
  TYPE(Obs1dBangle),      INTENT(in)    :: y                ! Observation vector
  REAL(wp), DIMENSION(:), INTENT(inout) :: y_ad             ! Observation adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: pwvp             ! Partial water vapour pressure
  REAL(wp), DIMENSION(x%n_lev)          :: pwvp_ad          ! Pwvp adjoint
 
  REAL(wp), DIMENSION(x%n_lev)          :: pdry             ! Dry pressure
  REAL(wp), DIMENSION(x%n_lev)          :: pdry_ad          ! Pdry adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: refrac           ! Refractivity
  REAL(wp), DIMENSION(x%n_lev)          :: refrac_ad        ! Refractivity adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: z_geop           ! |Geopotential height of model levels
  REAL(wp), DIMENSION(x%n_lev)          :: z_geop_ad        ! GPH adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: zcomp_dry_inv    ! Dry compressibility
  REAL(wp), DIMENSION(x%n_lev)          :: zcomp_dry_inv_ad ! Dry compressibility adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: zcomp_wet_inv    ! Wet compressibility
  REAL(wp), DIMENSION(x%n_lev)          :: zcomp_wet_inv_ad ! Wet compressibility adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: h                ! Geometric height
  REAL(wp), DIMENSION(x%n_lev)          :: h_ad             ! Geometric height adjoint

  REAL(wp), DIMENSION(x%n_lev)          :: impact           ! Impact parameter
  REAL(wp), DIMENSION(x%n_lev)          :: impact_ad        ! Impact parameter adjoint

  REAL(wp), DIMENSION(SIZE(y%bangle))   :: bangle

  REAL(wp)                              :: kap1,kap2,kap3   ! Refractivity coefficients used in routine

  REAL(wp)                              :: R_peak           ! Used in iono routines
  REAL(wp)                              :: R_peak_ad        ! Used in iono routines adjoint

!-------------------------------------------------------------------------------
! 1. Reset local adjoint variables
!-------------------------------------------------------------------------------
  
  pwvp_ad   = 0.0_wp
  refrac_ad = 0.0_wp
  h_ad      = 0.0_wp
  impact_ad = 0.0_wp
  
! for compressibility

  pdry_ad   = 0.0_wp
  z_geop_ad = 0.0_wp
  zcomp_dry_inv_ad = 0.0_wp
  zcomp_wet_inv_ad = 0.0_wp

! for iono

  R_peak_ad = 0.0_wp

!-------------------------------------------------------------------------------
! 2. Non ideal gas options 
!-------------------------------------------------------------------------------

! set inverse of compressibilities

  zcomp_dry_inv(:) = 1.0_wp
  zcomp_wet_inv(:) = 1.0_wp
   
! initialise geopotential heights
   
  z_geop(:) = x%geop(:)

  IF (x%non_ideal) THEN

! if non ideal gas calculation, use adjusted coefficients

     kap1 = kappa1_comp
     kap2 = kappa2_comp
     kap3 = kappa3_comp

!    calculate compressibilty and adjust geopotential heights in z_geop

     CALL ropp_fm_compress(x,z_geop,zcomp_dry_inv,zcomp_wet_inv)

  ELSE

     kap1 = kappa1
     kap2 = kappa2
     kap3 = kappa3

  ENDIF


!-------------------------------------------------------------------------------
! 3. Recompute forward model variables
!-------------------------------------------------------------------------------

! 3.1 Calculate water vapor pressure

  pwvp   = x%pres * x%shum / (epsilon_water + (1.0_wp-epsilon_water)*x%shum)

  pdry = x%pres - pwvp

! 3.2 Calculate refractivity

  refrac = kap1 * pdry * zcomp_dry_inv / x%temp    + &
           kap2 * pwvp * zcomp_wet_inv / x%temp**2 + &
           kap3 * pwvp * zcomp_wet_inv / x%temp

!-------------------------------------------------------------------------------
! 5. Calculate impact parameter
!-------------------------------------------------------------------------------

  h = geopotential2geometric(x%lat, z_geop)

  IF (y%undulation > ropp_MDTV) THEN
    impact = (1.0_wp + refrac*1.0e-6_wp) * (h + y%r_curve + y%undulation)
  ELSE
    CALL message(msg_warn, "Undulation missing. " // &
                 "Will assume to be zero when calculating full and " // &
                 "adjoint impact parameters.")
    impact = (1.0_wp + refrac*1.0e-6_wp) * (h + y%r_curve)
  END IF

!-------------------------------------------------------------------------------
! 6. Calculate bending angles
!-------------------------------------------------------------------------------

!  call ropp_fm_abel(impact, refrac, y%impact, bangle)  ! don't need this for neutral BAs
  
!-------------------------------------------------------------------------------
! 7. Adjoint of ionospheric bending angle computation **IF** we're using L1,L2
!-------------------------------------------------------------------------------

  IF (x%direct_ion) THEN

! need the neutral bending for a test in ropp_fm_iono

    CALL ropp_fm_abel &
      (impact, refrac, x%temp, y%r_curve, x%new_bangle_op, y%impact, bangle)

! now call adjoint

    R_peak = x%H_peak + y%r_curve

    CALL ropp_fm_iono_bangle_ad(x%Ne_max,    R_peak,    x%H_width, &
                                x_ad%Ne_max, R_peak_ad, x_ad%H_width, &
                                y%n_L1, y%impact, bangle, y_ad)

    x_ad%H_peak = x_ad%H_peak + R_peak_ad

    R_peak_ad = 0.0_wp

  END IF

!-------------------------------------------------------------------------------
! 8. Adjoint of bending angle computation
!-------------------------------------------------------------------------------

  CALL ropp_fm_abel_ad&
  (impact, refrac, x%temp, x_ad%temp, y%r_curve, x%new_bangle_op, &
   y%impact, impact_ad, refrac_ad, y_ad)

!-------------------------------------------------------------------------------
! 9. Adjoint of impact parameter calculation
!-------------------------------------------------------------------------------

  IF (y%undulation > ropp_MDTV) THEN
    refrac_ad = refrac_ad + 1.0e-6_wp * impact_ad * (h + y%r_curve + y%undulation)
  ELSE
    refrac_ad = refrac_ad + 1.0e-6_wp * impact_ad * (h + y%r_curve)
  END IF

  h_ad      = h_ad + impact_ad * (1.0_wp + 1.0e-6_wp*refrac)

  impact_ad = 0.0_wp

  z_geop_ad = z_geop_ad + h_ad*(y%r_earth/(((y%g_sfc/g_wmo)*y%r_earth)-z_geop) &
               + y%r_earth*z_geop/((((y%g_sfc/g_wmo)*y%r_earth)-z_geop)**2))

  h_ad      = 0.0_wp

! 10. Adjoint of refractivity calculation
!----------------------------------------------

  pdry_ad = pdry_ad +  refrac_ad * kap1 * zcomp_dry_inv/ x%temp
  
  pwvp_ad = pwvp_ad +refrac_ad * (kap2 * &
            zcomp_wet_inv/ x%temp**2 + kap3 * zcomp_wet_inv/ x%temp)
  
  zcomp_dry_inv_ad = zcomp_dry_inv_ad + refrac_ad * kap1 * pdry / x%temp
  
  zcomp_wet_inv_ad = zcomp_wet_inv_ad + refrac_ad * ( & 
                   kap2 * pwvp / x%temp**2 + kap3 * pwvp / x%temp) 
  
  x_ad%temp   = x_ad%temp - refrac_ad * (kap1*pdry*zcomp_dry_inv/(x%temp**2)     &
                                + 2.0_wp*kap2*pwvp*zcomp_wet_inv/(x%temp**3) &
                                  + kap3*pwvp*zcomp_wet_inv/x%temp**2)

  refrac_ad = 0.0_wp

! 11 Adjoint of water vapor and dry air pressure calculation
!-----------------------------------------------------
  x_ad%pres = x_ad%pres + pdry_ad 
  pwvp_ad = pwvp_ad - pdry_ad
  pdry_ad = 0.0_wp 
  
  x_ad%pres = x_ad%pres + pwvp_ad   &
              * (x%shum/(epsilon_water+(1.0_wp-epsilon_water)*x%shum))
  x_ad%shum = x_ad%shum + pwvp_ad   &
              * (x%pres/(epsilon_water+(1.0_wp-epsilon_water)*x%shum)         &
                 - x%pres * x%shum * (1.0_wp - epsilon_water) &
                   / ((epsilon_water + (1.0_wp-epsilon_water)*x%shum)         &
                    * (epsilon_water + (1.0_wp-epsilon_water)*x%shum)))
  pwvp_ad = 0.0_wp

! 11 Adjoint of the compressibilty calculation
!-----------------------------------------------------------
  
  IF (x%non_ideal) THEN

! call the compressibility routines

      CALL ropp_fm_compress_ad &
     &(x,x_ad,z_geop_ad,zcomp_dry_inv_ad,zcomp_wet_inv_ad)  

  ENDIF

  x_ad%geop = x_ad%geop + z_geop_ad
  z_geop_ad = 0.0_wp

  zcomp_dry_inv_ad(:) = 0.0_wp
  zcomp_wet_inv_ad(:) = 0.0_wp


END SUBROUTINE ropp_fm_bangle_1d_ad
