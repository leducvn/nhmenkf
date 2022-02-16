! $Id: ropp_pp_types.f90 4452 2015-01-29 14:42:02Z idculv $

MODULE ropp_pp_types

!****m* Modules/ropp_pp_types *
!
! NAME
!    ropp_pp_types - Type declarations for ROPP PP library.
!
! SYNOPSIS
!    use ropp_pp_types
!
! DESCRIPTION
!    This module provides dervived type / structure definitions for the
!    pre-processor routines in ROPP.
!
! NOTES
!    The ropp_pp_types module is loaded automatically with the ropp_pp module.

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
! 2. Pre-processor configuration data types
!-------------------------------------------------------------------------------

!****d* Datatypes/PPConfig
!
! NAME
!    PPConfig - Configuration options for excess phase to bending angle and 
!               ionospheric correction processing
!
! SYNOPSIS
!    use ropp_pp
!      ...
!    type(PPConfig) :: config
!
! SOURCE
!
  TYPE ppConfig

     LOGICAL  :: obs_ok   = .TRUE.      ! Observations QC flag
     
     LOGICAL  :: output_tdry = .TRUE.   ! Flag to output 'dry' parameters

     LOGICAL  :: output_diag = .FALSE.  ! Flag to output additional diagnostics

     !! Occultation (excess phase to bending angle) processing

     ! Excess phase to bangle method (GO = geometric optics, WO = wave optics)
     CHARACTER(len=10) :: occ_method = "WO"
     ! Filter type for phase differentiation (slpoly = sliding polynomial, optest = optimal estimation filter)
     CHARACTER(len=10) :: filter_method = "slpoly"
     ! Filter width for smoothed GO bending angles (m)
     REAL(wp) :: fw_go_smooth =  3000.0_wp
     ! Filter width for computing full resolution GO bending angles (m)
     REAL(wp) :: fw_go_full   =  3000.0_wp
     ! Filter width for wave optics bending angle above 7 km (m)
     REAL(wp) :: fw_wo        =  2000.0_wp
     ! Filter width for wave optics bending angle below 7 km (m)
     REAL(wp) :: fw_low       =  -1000.0_wp
     ! Maximum height for wave optics processing (m)
     REAL(wp) :: hmax_wo      = 25000.0_wp
     ! Fractional cut-off limit for amplitude
     REAL(wp) :: Acut     = 0.0_wp
     ! Cut-off limit for impact height (m)
     REAL(wp) :: Pcut     = -2000.0_wp
     ! Cut-off limit for bending angle (rad)
     REAL(wp) :: Bcut     = 0.1_wp
     ! Cut-off limit for straight-line tangent altitude (m)
     REAL(wp) :: Hcut     = -250000.0_wp
     ! Complex filter flag
     INTEGER  :: CFF = 3
     ! Shadow border width (m)
     REAL(wp) :: dsh = 200.0_wp
     ! Degraded L2 data flag
     LOGICAL  :: opt_DL2  = .TRUE.
     ! Calculate and output spectra flag
     LOGICAL  :: opt_spectra = .FALSE.
     ! Path to EGM96 model coefficients file
     CHARACTER(len=80) :: egm96 = "egm96.dat"
     ! Path to EGM96 model corrections file
     CHARACTER(len=80) :: corr_egm96 = "corrcoef.dat"
     ! Path to external navigation bit file
     CHARACTER(len=80) :: navbit_file = " "

     !! Inversion (bending angle to refractivity) processing

     ! Ionospheric correction method (NONE = linear combination, MSIS = full,
     !                                GMSIS = full with global MSIS search,
     !                                BG = full with background profile)
     CHARACTER(len=10) :: method = "GMSIS"
     ! Statistical optimization method (so = statistical optimisation (default)
     !                                  lcso = linear combination+stat opt)
     CHARACTER(len=10) :: so_method = "so"
     ! Abel integral algorithm (LIN = linear, EXP = exponential)
     CHARACTER(len=10) :: abel = "LIN"
     ! Model coefficients file path
     CHARACTER(len=80) :: mfile = "MSIS_coeff.nc"
     ! Background atmospheric profile file path
     CHARACTER(len=80) :: bfile = "n/a"
     ! Local radius of curvature (m)
     REAL(wp) :: r_curve
     ! Number of input data points
     INTEGER  :: npoints
     ! Step of standard impact parameter grid (m)
     REAL(wp) :: dpi      = 100.0_wp

     ! Minimum impact parameter (m)
     REAL(wp) :: pmin
     ! Maximum impact parameter (m)
     REAL(wp) :: pmax
     ! Polynomial degree for smoothing regression
     INTEGER  :: np_smooth = 3
     ! Filter width for smoothing profile
     REAL(wp) :: fw_smooth = 1000.0_wp

     ! Number of parameters used for model fit regression
     REAL(wp) :: nparm_fit = 2
     ! Lower limit for model fit regression (m)
     REAL(wp) :: hmin_fit = 20000.0_wp
     ! Upper limit for model fit regression (m)
     REAL(wp) :: hmax_fit = 70000.0_wp
     ! a priori standard dev of regression factor
     REAL(wp) :: omega_fit = 0.3_wp

     ! Ionospheric correction filter width (m)
     REAL(wp) :: f_width  = 2000.0_wp
     ! Step of homogeneous impact p. grid (m)
     REAL(wp) :: delta_p  =  20.0_wp
     ! External ionospheric smoothing scale (m)
     REAL(wp) :: s_smooth = 2000.0_wp
     ! Lower limit of ionospheric signal (m)
     REAL(wp) :: z_ion    = 50000.0_wp
     ! Upper limit of stratospheric signal (m)
     REAL(wp) :: z_str    = 35000.0_wp
     ! Upper limit of tropospheric signal (m)
     REAL(wp) :: z_ltr    = 12000.0_wp
     ! Number of points for smoothing (odd)
     INTEGER  :: n_smooth = 11
     ! A priori model error std.dev. (Negative ==> use dynamical estimate.)
     REAL(wp) :: model_err = -0.5_wp
     ! Height of atmosphere for inversion (m)
     REAL(wp) :: ztop_invert = 150000.0_wp
     ! Step of highest part inversion grid (m)
     REAL(wp) :: dzh_invert = 50.0_wp
     ! Interval for regression in inversion (m)
     REAL(wp) :: dzr_invert = 20000.0_wp

  END TYPE ppConfig
!
!****

!-------------------------------------------------------------------------------
! 3. Configuration read structure
!-------------------------------------------------------------------------------

!****d* Datatypes/KeyConfig
!
! NAME
!    KeyConfig - Key-value pairs defined in PP config files
!
! SYNOPSIS
!    use ropp_pp_types
!      ...
!    type(KeyConfig) :: key_value
!
! NOTES
!
! SEE ALSO
!
! SOURCE
!
  TYPE KeyConfig

     CHARACTER(len = 1024), DIMENSION(:), POINTER :: keys   => null()
     CHARACTER(len = 4096), DIMENSION(:), POINTER :: values => null()

  END TYPE KeyConfig
!
!****
!
!-------------------------------------------------------------------------------
! 4. Pre-processor output diagnostics data types
!-------------------------------------------------------------------------------

!****d* Datatypes/PPDiag
!
! NAME
!    PPDiag -
!
! SYNOPSIS
!    use ropp_pp
!      ...
!    type(PPDiag) :: diag
!
! SOURCE
!
  TYPE ppDiag
     REAL(wp) :: L2_badness    ! L2 phase correction badness score
     REAL(wp) :: sq            ! Statistical optimisation badness score
     ! CT processing impact parameter grid (m)
     REAL(wp), DIMENSION(:), POINTER :: CTimpact => null()
     ! CT processing amplitude 
     REAL(wp), DIMENSION(:), POINTER :: CTamplitude => null()
     ! CT processing smoothed amplitude 
     REAL(wp), DIMENSION(:), POINTER :: CTamplitude_smt => null()
     ! CT processing impact parameter grid (m)
     REAL(wp), DIMENSION(:), POINTER :: CTimpactL2 => null()
     ! CT processing amplitude 
     REAL(wp), DIMENSION(:), POINTER :: CTamplitudeL2 => null()
     ! CT processing smoothed amplitude 
     REAL(wp), DIMENSION(:), POINTER :: CTamplitudeL2_smt => null()
     ! Ionospheric bending angle in L1
     REAL(wp), DIMENSION(:), POINTER :: ba_ion => null()
     ! Error covariance of ionospheric bending angle (rad**2)
     REAL(wp), DIMENSION(:), POINTER :: err_ion => null()
     ! Error covariance of neutral bending angle (rad**2)
     REAL(wp), DIMENSION(:), POINTER :: err_neut => null()
     ! Weight of data (data:data+clim) in profile
     REAL(wp), DIMENSION(:), POINTER :: wt_data => null()

  END TYPE ppDiag
!
!****
!
!-------------------------------------------------------------------------------
! 5. Tropopause height diagnostics QC flag definitions
!-------------------------------------------------------------------------------

!****p* Parameters/TPH_QC_flag
!
! NAME
!    TPH QC flag parameters - Parameters that define the components 
!                             of the overall TPH QC flag.
!
! SOURCE
!
  INTEGER, PARAMETER :: TPH_QC_data_invalid     =  0  ! Input data validity check
  INTEGER, PARAMETER :: TPH_QC_prof_depth       =  1  ! Input data depth check
  INTEGER, PARAMETER :: TPH_QC_prof_height      =  2  ! Input data height check
  INTEGER, PARAMETER :: TPH_QC_CT_smooth_above  =  3  ! Cov trans sharpness above TPH check
  INTEGER, PARAMETER :: TPH_QC_CT_smooth_below  =  4  ! Cov trans sharpness below TPH check
  INTEGER, PARAMETER :: TPH_QC_double_trop      =  5  ! Double tropopause detected
  INTEGER, PARAMETER :: TPH_QC_too_low          =  6  ! TPH minimum height check
  INTEGER, PARAMETER :: TPH_QC_too_high         =  7  ! TPH maximum height check
  INTEGER, PARAMETER :: TPH_QC_reserved8        =  8  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved9        =  9  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved10       = 10  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved11       = 11  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved12       = 12  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved13       = 13  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved14       = 14  ! Reserved
  INTEGER, PARAMETER :: TPH_QC_reserved15       = 15  ! Reserved
!
!****

END MODULE ropp_pp_types
