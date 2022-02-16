! $Id: ropp_pp_constants.f90 3551 2013-02-25 09:51:28Z idculv $

MODULE ropp_pp_constants

!****m* Modules/ropp_pp_constants *
!
! NAME
!    ropp_pp_constants - Module providing meteorological and physical constants.
!
! SYNOPSIS
!    use ropp_pp_constants
! 
! DESCRIPTION
!    This module provides numerical values for meteorological and other
!    physical constants which are used throughout the ROPP package.
!
! SEE ALSO
!
!    Thermodynamical constants:  R_dry, R_vap, C_p,
!                                mw_dry_air, mw_water, epsilon_water
!    GPSRO carrier frequencies:  f_L1, f_L2
!    Refractivity constants:     kappa1, kappa2
!    Light speed in vacuum:      c_light
!    Gravitational acceleration: g_wmo 
!    Mathematical constants:     pi
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

  USE typesizes, ONLY: wp => EightByteReal

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! 1. GPSRO carrier L1, L2 frequencies
!-------------------------------------------------------------------------------

!****t* GPSRO/L1L2
!
! DESCRIPTION
!    Carrier frequencies L1 and L2 as used in ROPP pre-processor library.
!
!****

!****p* GPSRO/L1 *
!
! NAME
!    f_L1 - Carrier frequency for L1 signal 
!
! SOURCE
!
  REAL(wp), PARAMETER :: f_L1 = 1.57542e9_wp         !   Hz
!
!****

!****p* GPSRO/L2 *
!
! NAME
!    f_L2 - Carrier frequency for L2 signal 
!
! SOURCE
!
  REAL(wp), PARAMETER :: f_L2 = 1.22760e9_wp         !   Hz
!
!****

!-------------------------------------------------------------------------------
! 2. Other physical constants
!-------------------------------------------------------------------------------

!****t* Constants/Earth
!
! DESCRIPTION
!    Physical constants related to the Earth as used in the ROPP Forward Model
!    library. 
!
! SEE ALSO
!    g_wmo
!    c_light
!    R_dry
!    R_vap
!    C_p
!    mw_dry_air
!    mw_water
!    epsilon_water
!
!****

!****p* Earth/g_wmo *
!
! NAME
!    g_wmo - Gravitational acceleration (WMO standard value)
!
! SOURCE
!
  REAL(wp), PARAMETER :: g_wmo = 9.80665_wp         !   m s^{-2}
!
!****

!****p* Physical/c_light *
!
! NAME
!    c_light - Light speed in vaccuum (m/s)
!
! SOURCE
!
  REAL(wp), PARAMETER :: c_light = 299792458.0_wp         !   m s^{-1}
!
!****

!****p* Thermodynamics/R_dry *
!
! NAME
!    R_dry - Gas constant of dry air
!
! SOURCE
!
!!  REAL(wp), PARAMETER :: R_dry = 287.0597_wp        !   J K^{-1} kg^{-1}
  REAL(wp), PARAMETER :: R_dry = 287.05_wp        !   J K^{-1} kg^{-1}
!
!****

!****p* Thermodynamics/R_vap *
!
! NAME
!    R_vap - Gas constant of water vapor
!
! SOURCE
!
!!  REAL(wp), PARAMETER :: R_vap = 461.5250_wp        !   J K^{-1} kg^{-1}
  REAL(wp), PARAMETER :: R_vap = 461.51_wp        !   J K^{-1} kg^{-1}
!
!****

!****p* Thermodynamics/C_p *
!
! NAME
!    C_p - Specific heat capacity of dry air at constant pressure
!
! SOURCE
!
!!  REAL(wp), PARAMETER :: C_p = 1005.7_wp       !   J K^{-1} kg^{-1}
  REAL(wp), PARAMETER :: C_p = 1004.6_wp       !   J K^{-1} kg^{-1}
!
!****

!****p* Thermodynamics/mw_dry_air *
!
! NAME
!    mw_dry_air - Molecular weight of dry air
!
! SOURCE
!
  REAL(wp), PARAMETER :: mw_dry_air = 28.9648e-3_wp !   kg mol^{-1}
!
!****

!****p* Thermodynamics/mw_water *
!
! NAME
!    mw_water - Molecular weight of water
!
! SOURCE
!
  REAL(wp), PARAMETER :: mw_water = 18.01528e-3_wp  !   kg mol^{-1}
!
!****

!****p* Thermodynamics/epsilon_water *
!
! NAME
!    epsilon_water - Ratio of molecular weight of water to that of dry air.
!
! SOURCE
!
  REAL(wp), PARAMETER :: epsilon_water = 287.05_wp / 461.51_wp  !mw_water / mw_dry_air
!
!****

!-------------------------------------------------------------------------------
! 3. RO specific constants
!-------------------------------------------------------------------------------

!****t* Constants/RadioOccultation
!
! DESCRIPTION
!    Physical constants specific for radio occultation measurement as used in
!    the ROPP Forward Model library.
!
! SEE ALSO
!    kappa1
!    kappa2
!
!****

!****p* RadioOccultation/kappa1 *
!
! NAME
!    kappa1 - Coefficient to calculate dry refractivity contribution.
!
! SOURCE
!
  REAL(wp), PARAMETER :: kappa1 = 77.60e-2_wp      !   K Pa^{-1}
!
!****

!****p* RadioOccultation/kappa2 *
!
! NAME
!    kappa2 - Coefficient to calculate moist refractivity contribution.
!
! SOURCE
!
  REAL(wp), PARAMETER :: kappa2 = 3.73e3_wp        !   K^2 Pa^{-1}
!
!****

!-------------------------------------------------------------------------------
! 4. Mathematical constants
!------------------------------------------------------------------------------

!****t* Constants/Maths
!
! DESCRIPTION
!    Mathematical constants as used in the ROPP Forward Model library. 
!
! SEE ALSO
!    pi
!
!****

!****p* Maths/pi *
!
! NAME
!    pi  
!
! SOURCE
!
  REAL(wp), PARAMETER :: pi = 3.141592653589793238_wp   
!
!****

END MODULE ropp_pp_constants
