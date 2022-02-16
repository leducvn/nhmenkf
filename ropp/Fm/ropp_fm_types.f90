! $Id: ropp_fm_types.f90 4010 2014-01-10 11:07:40Z idculv $

MODULE ropp_fm_types

!****t* ForwardModels/Datatypes
!
! DESCRIPTION
!    Data types / structures used by the forward models of ROPP.
!
! SEE ALSO
!    ropp_fm_types
!
!    ropp_fm_state_1d_type
!    ropp_fm_obs_type
!
!****

!****t* ForwardModels/Parameters
!
! DESCRIPTION
!    Parameters used by the forward models of ROPP.
!
! SEE ALSO
!    wp
!
!****

!****m* Modules/ropp_fm_types *
!
! NAME
!    ropp_fm_types - Type declarations and structure definitions for the ROPP
!                    Forward Model library.
!
! SYNOPSIS
!    use ropp_fm
!
! DESCRIPTION
!    This module provides dervived type / structure definitions for the
!    forward model routines in ROPP.
!
! NOTES
!    The ropp_fm_types module is loaded automatically with the ropp_fm module.
!
! SEE ALSO
!    ropp_fm_state_1d_type
!    ropp_fm_obs_type
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
! 1. Double precision variables
!-------------------------------------------------------------------------------

!****p* Parameters/wp *
!
! NAME
!    wp - Kind parameter for double precision floating point numbers.
!
! SOURCE
!
  USE typesizes, ONLY: wp => EightByteReal
!
!****

!-------------------------------------------------------------------------------
! 2. Matrix type
!-------------------------------------------------------------------------------
!****d* Datatypes/matrix_pp *
!
! NAME
!    matrix_pp - Type declaration for positive definite (packed) matrix.
!
! SOURCE
!
  USE matrix_types
!
!****

!-------------------------------------------------------------------------------
! 3. Observation vector for 1d refractivity
!-------------------------------------------------------------------------------

!****d* Datatypes/Obs1dRefrac *
!
! NAME
!    Obs1dRefrac - 1d Refractivity "observations" vector data type.
!
! SOURCE
!
  TYPE Obs1dRefrac
     REAL(wp), DIMENSION(:), POINTER :: refrac  => null() ! Refractivity values
     REAL(wp), DIMENSION(:), POINTER :: geop    => null() ! Geopotential height
     REAL(wp), DIMENSION(:), POINTER :: weights => null() ! Observation weights
     REAL(wp)                        :: lon               ! Longitude
     REAL(wp)                        :: lat               ! Latitude
     REAL(wp)                        :: time              ! JulianTime of obs
     TYPE(matrix_pp)                 :: cov               ! Error covariance
     LOGICAL                         :: obs_ok            ! Observations ok?
     LOGICAL                         :: cov_ok            ! Covariance ok?
  END TYPE Obs1dRefrac
!
!****

!-------------------------------------------------------------------------------
! 4. Observation vector for 1d bending angles
!-------------------------------------------------------------------------------

!****d* Datatypes/Obs1dBangle *
!
! NAME
!    Obs1dBangle - 1d bending angle "observations" vector data type.
!
! SOURCE
!
  TYPE Obs1dBangle
     REAL(wp), DIMENSION(:), POINTER :: bangle  => null() ! Bending angle values
     REAL(wp), DIMENSION(:), POINTER :: impact  => null() ! Impact parameter
     REAL(wp), DIMENSION(:), POINTER :: weights => null() ! Observation weights

     REAL(wp), DIMENSION(:), POINTER :: rtan    => null() ! Radius of tangent point (2D operator only)
     REAL(wp),DIMENSION(:,:),POINTER :: a_path  => null() ! nr sin(phi) at ray endpoints (2D operator only)
     REAL(wp)                        :: azimuth           ! Azimuthal angle (2D operator only)

     INTEGER                         :: nobs              ! Number of obs
     INTEGER                         :: n_L1              ! No. of L1 bangles when simulating L1,L2
     REAL(wp)                        :: lon               ! Longitude
     REAL(wp)                        :: lat               ! Latitude
     REAL(wp)                        :: time              ! JulianTime of obs
     REAL(wp)                        :: g_sfc             ! Surface gravity
     REAL(wp)                        :: r_earth           ! Effective radius
     REAL(wp)                        :: r_curve           ! Radius of curvature
     REAL(wp)                        :: undulation        ! Undulation
     TYPE(matrix_pp)                 :: cov               ! Error covariance
     LOGICAL                         :: obs_ok            ! Observations ok?
     LOGICAL                         :: cov_ok            ! Covariance ok?
  END TYPE Obs1dBangle
!
!****

!-------------------------------------------------------------------------------
! 5. Generic one dimensional state vector data type
!
!-------------------------------------------------------------------------------
!****d* Datatypes/State1dFM *
!
! NAME
!     - 1d state vector data type.
!
! SOURCE
!
  TYPE State1dFM

     REAL(wp), DIMENSION(:), POINTER :: state  => null()        ! State vector
     REAL(wp), DIMENSION(:), POINTER :: temp   => null()        ! Temperature (K)
     REAL(wp), DIMENSION(:), POINTER :: shum   => null()        ! Spec humidity (kg/kg)
     REAL(wp), DIMENSION(:), POINTER :: pres   => null()        ! Pressure (Pa)
     REAL(wp), DIMENSION(:), POINTER :: geop   => null()        ! Geopotential height (gpm)
     REAL(wp), DIMENSION(:), POINTER :: ak     => null()        ! Model level coefficients
     REAL(wp), DIMENSION(:), POINTER :: bk     => null()        ! Model level coefficients
     REAL(wp)                        :: Ne_max                  ! Peak elec.density
     REAL(wp)                        :: H_peak                  ! Height of peak elec. dens.
     REAL(wp)                        :: H_width                 ! Width of Chapman layer
     REAL(wp)                        :: geop_sfc                ! Surface geop height
     INTEGER                         :: n_lev                   ! Number of levels
     INTEGER                         :: n_chap = 0              ! No. of Chapman layers in iono calc
     REAL(wp)                        :: lon                     ! Longitude
     REAL(wp)                        :: lat                     ! Latitude
     REAL(wp)                        :: time                    ! JulianTime of bg
     TYPE(matrix_pp)                 :: cov                     ! Error covariance
     LOGICAL                         :: state_ok                ! State vector ok?
     LOGICAL                         :: cov_ok                  ! Covariance ok?
     LOGICAL                         :: use_logp      = .FALSE. ! Using log(pres)?
     LOGICAL                         :: use_logq      = .FALSE. ! Using log(shum)?
     LOGICAL                         :: non_ideal     = .FALSE. ! Non-ideal gas flag
     LOGICAL                         :: check_qsat    = .FALSE. ! Check against saturation
     LOGICAL                         :: new_ref_op    = .FALSE. ! New REF interpolation
     LOGICAL                         :: new_bangle_op = .FALSE. ! New BA computation
     LOGICAL                         :: direct_ion    = .FALSE. ! Calculate L1,L2 bangles

  END TYPE State1dFM
!
!****

!-------------------------------------------------------------------------------
! 6. Generic two dimensional state vector data type
!
!-------------------------------------------------------------------------------
!****d* Datatypes/State2dFM *
!
! NAME
!     - 2d state vector data type.
!
! SOURCE
!
  TYPE State2dFM

     REAL(wp), DIMENSION(:,:), POINTER :: temp   => null()     ! Temperature (K)
     REAL(wp), DIMENSION(:,:), POINTER :: shum   => null()     ! Specific humidity (kg/kg)
     REAL(wp), DIMENSION(:,:), POINTER :: pres   => null()     ! Pressure (Pa)
     REAL(wp), DIMENSION(:,:), POINTER :: geop   => null()     ! Geopotential height (gpm)
     REAL(wp),   DIMENSION(:), POINTER :: ak     => null()     ! Model level coefficients
     REAL(wp),   DIMENSION(:), POINTER :: bk     => null()     ! Model level coefficients
     REAL(wp),   DIMENSION(:), POINTER :: geop_sfc             ! Surface geop height
     REAL(wp),   DIMENSION(:), POINTER :: pres_sfc             ! Surface pressure
     REAL(wp), DIMENSION(:,:), POINTER :: refrac => null()     ! Refractivity on 2D grid
     REAL(wp), DIMENSION(:,:), POINTER :: nr     => null()     ! nr product on 2D grid
     INTEGER                           :: n_lev                ! Number of levels
     INTEGER                           :: n_horiz              ! Number of horizontal points
     REAL(wp)                          :: dtheta               ! Angle between profiles in plane
     REAL(wp), DIMENSION(:), POINTER   :: lon                  ! Longitude
     REAL(wp), DIMENSION(:), POINTER   :: lat                  ! Latitude
     REAL(wp)                          :: time                 ! JulianTime of bg
     LOGICAL                           :: state_ok             ! State vector ok?
     LOGICAL                           :: non_ideal  = .FALSE. ! Non-ideal gas flag
     LOGICAL                           :: check_qsat = .FALSE. ! Check against saturation

  END TYPE State2dFM
!
!****

END MODULE ropp_fm_types
