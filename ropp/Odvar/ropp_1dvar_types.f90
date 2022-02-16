! $Id: ropp_1dvar_types.f90 4452 2015-01-29 14:42:02Z idculv $

!****m* Modules/ropp_1dvar_types *
!
! NAME
!    ropp_1dvar_types - Type declarations for the ROPP 1DVAR libary.
!
! SYNOPSIS
!    use ropp_1dvar
! 
! DESCRIPTION
!    This Fortran module supports the ROPP 1DVar library and provides
!    derived data types used by the ropp_1dvar library. It must be loaded
!    for all applications using the ropp_io library / package. Note that
!    loading the ropp_1dvar module includes the ropp_1dvar_types module.
!
! NOTES
!
! SEE ALSO
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

MODULE ropp_1dvar_types

!-------------------------------------------------------------------------------
! 1. Other modules
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix_types
  USE ropp_fm_types

  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax_fname = 2048

!-------------------------------------------------------------------------------
! 2. Configuration data types
!-------------------------------------------------------------------------------

!****d* Datatypes/minROPPConfig
!
! NAME
!    minROPPConfig - 
!
! SYNOPSIS
!    use ropp_1dvar_types
!      ...
!    type(VarConfig) :: config
!
! NOTES
!
! SEE ALSO
!    VarConfig
!
! SOURCE
!
  TYPE minROPPConfig
     CHARACTER(len = nmax_fname) :: method       = "MINROPP"
     CHARACTER(len = nmax_fname) :: log_file     = "screen"
     INTEGER                     :: impres       =    0
     INTEGER                     :: n_iter       = 1500
     INTEGER                     :: mode         =    0
     INTEGER                     :: n_updates    =   50
     REAL(wp)                    :: eps_grad     = 1.0e-8_wp
     REAL(wp)                    :: dx_min       = 1.0e-16_wp
  END TYPE minROPPConfig
!
!****


!****d* Datatypes/VarConfig
!
! NAME
!    VarConfig - 
!
! SYNOPSIS
!    use ropp_1dvar_types
!      ...
!    type(VarConfig) :: config
!
! NOTES
!    The VarConfigType structure is composed out of several other structures; 
!    see the user guide for a breakdown of the actual element names.
!
! SEE ALSO
!    minROPPConfig
!
! SOURCE
!
  TYPE VarConfig

!    Input and output files

     ! Input observation file
     CHARACTER(len = nmax_fname) :: obs_file         = "ropp_obs.nc"      
     ! Observation error covariance method
     CHARACTER(len = nmax_fname) :: obs_covar_method = "VSDC"             
     ! Observation error correlation file
     CHARACTER(len = nmax_fname) :: obs_corr_file    = "ropp_obs_corr.nc" 
     ! Input background file
     CHARACTER(len = nmax_fname) :: bg_file          = "ropp_bg.nc"       
     ! Output file
     CHARACTER(len = nmax_fname) :: out_file         = "ropp_out.nc"      
     ! Background error covariance method
     CHARACTER(len = nmax_fname) :: bg_covar_method  = "VSFC"             
     ! Background error correlation file
     CHARACTER(len = nmax_fname) :: bg_corr_file     = "ropp_bg_corr.nc"  

!    Observation cutoff range

     ! Minimum obs height to include in 1dVar analysis (km)
     REAL(wp)                    :: min_1dvar_height         = -10.0_wp 
     ! Maximum obs height to include in 1dVar analysis (km)
     REAL(wp)                    :: max_1dvar_height         = 60.0_wp   

!    Generic quality control

     ! Apply obs vs. bg colocation check?
     LOGICAL                     :: genqc_colocation_apply     = .TRUE.      
     ! Maximum obs vs. bg great circle distance (km)
     REAL(wp)                    :: genqc_max_distance         = 300.0_wp    
     ! Maximum obs vs. bg temporal seperation (sec)
     REAL(wp)                    :: genqc_max_time_sep         = 10800.0_wp    
     ! Minimium temperature value (K)
     REAL(wp)                    :: genqc_min_temperature      = 150.0_wp
     ! Maximium temperature value (K)
     REAL(wp)                    :: genqc_max_temperature      = 350.0_wp
     ! Minimium specific humidity value (g/kg)
     REAL(wp)                    :: genqc_min_spec_humidity    =  0.0_wp
     ! Maximium specific humidity value (g/kg)
     REAL(wp)                    :: genqc_max_spec_humidity    = 50.0_wp
     ! Minimium impact parameter value (m)
     REAL(wp)                    :: genqc_min_impact           = 6.2e6_wp   
     ! Maximium impact parameter value (m)
     REAL(wp)                    :: genqc_max_impact           = 6.6e6_wp  
     ! Minimium bending angle value (rad)
     REAL(wp)                    :: genqc_min_bangle           = -1.0e-4_wp
     ! Maximium bending angle value (rad)
     REAL(wp)                    :: genqc_max_bangle           =  0.1_wp 
     ! Minimium geopotential height for refractivity value (m)
     REAL(wp)                    :: genqc_min_geop_refrac      = -1.0e3_wp
     ! Maximum geopotential height for refractivity value (m)
     REAL(wp)                    :: genqc_max_geop_refrac      = 1.e5_wp  
     ! Minimium refractivity value (N-units)
     REAL(wp)                    :: genqc_min_refractivity     = 0.0_wp 
     ! Maximium refractivity value (N-units)
     REAL(wp)                    :: genqc_max_refractivity     = 500.0_wp  
     ! Minimum threshold observation height (m)
     REAL(wp)                    :: genqc_min_obheight         = 20000.0_wp

     !    Background quality control

     ! Apply background quality control?
     LOGICAL                     :: bgqc_apply                 = .TRUE.    
     ! Data rejected if O-B > factor * sigma
     REAL(wp)                    :: bgqc_reject_factor         = 10.0_wp    
     ! Maximum percentage of data rejected
     REAL(wp)                    :: bgqc_reject_max_percent    = 50.0_wp     

!    Probability of Gross Error

     ! Apply PGE for quality control?
     LOGICAL                     :: pge_apply                  = .FALSE.    
     ! First guess PGE
     REAL(wp)                    :: pge_fg                     = 0.001_wp  
     ! Width of gross error plateau
     REAL(wp)                    :: pge_d                      = 10.0_wp   

!    Minimiser

     TYPE(minROPPConfig)         :: minROPP
     LOGICAL                     :: use_precond                = .TRUE.

!    Additional convergence checks

     ! Apply additional convergence checks?
     LOGICAL                     :: conv_check_apply           = .TRUE.    
     ! Minimum number of iterations required
     INTEGER                     :: conv_check_n_previous      = 2         
     ! State vector must change less than this * bg error
     REAL(wp)                    :: conv_check_max_delta_state = 0.1_wp    
     ! Cost function must change less than this
     REAL(wp)                    :: conv_check_max_delta_J     = 0.1_wp     

!    Additional output

     LOGICAL                     :: extended_1dvar_diag     = .FALSE.

!    Log(pressure) and Log(humidity) calculations

     LOGICAL                     :: use_logp                = .FALSE.
     LOGICAL                     :: use_logq                = .FALSE.

!    compressibility 
     
     LOGICAL                     :: compress                = .FALSE. 

!    Parameters for simple seasonal scaling of observation errors (BA and REF)
     REAL                        :: season_amp              = 0.0_wp
     REAL                        :: season_offset           = 0.0_wp
     REAL                        :: season_phase            = 0.0_wp

  END TYPE VarConfig
!
!****
!
!****d* Datatypes/KeyConfig
!
! NAME
!    KeyConfig - Key-value pairs defined in 1dVar config files
!
! SYNOPSIS
!    use ropp_1dvar_types
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
! 3. Diagnostic data types
!-------------------------------------------------------------------------------

!****d* Datatypes/VarDiag
!
! NAME
!    VarDiag - Diagnostics of a 1DVar retrieval.
!
! SYNOPSIS
!    use ropp_1dvar_types
!      ...
!    type(VarDiag) :: diag
!
! NOTES
!    The VarDiag structure is composed out of several other structures; 
!    see the user guide for a breakdown of the actual element names.
!
! SEE ALSO
!
! SOURCE
!
  TYPE VarDiag
     INTEGER                         :: n_data
     INTEGER                         :: n_bgqc_reject
     INTEGER                         :: n_pge_reject
     TYPE(Obs1DBangle)               :: bg_bangle
     TYPE(Obs1DRefrac)               :: bg_refrac
     REAL(wp), DIMENSION(:), POINTER :: OmB         => null() ! O - B
     REAL(wp), DIMENSION(:), POINTER :: OmB_sigma   => null() ! O - B std dev
     REAL(wp)                        :: pge_gamma             ! PGE gamma value
     REAL(wp), DIMENSION(:), POINTER :: pge         => null() ! PGE of obs
     REAL(wp), DIMENSION(:), POINTER :: pge_weights => null() ! PGE weights
     LOGICAL                         :: ok             ! Overall quality flag
     REAL(wp)                        :: J              ! Value of cost function
     REAL(wp)                        :: J_scaled ! Value of cost (scaled, 2J/m)
     REAL(wp)                        :: J_init ! Initial cost function
     REAL(wp), DIMENSION(:), POINTER :: J_bgr => null() ! BG cost profile
     REAL(wp), DIMENSION(:), POINTER :: J_obs => null() ! OB cost profile
     REAL(wp), DIMENSION(:), POINTER :: B_sigma => null() ! Bg obs std dev
     INTEGER                         :: n_iter         ! Number of iterations
     INTEGER                         :: n_simul        ! Number of simulations
     INTEGER                         :: min_mode       ! Minimiser exit mode
     TYPE(Obs1DBangle)               :: res_bangle     ! Analysis bending angle
     TYPE(Obs1DRefrac)               :: res_refrac     ! Analysis refractivity
     REAL(wp), DIMENSION(:), POINTER :: OmA         => null() ! O - A
     REAL(wp), DIMENSION(:), POINTER :: OmA_sigma   => null() ! O - A std dev
  END TYPE VarDiag
!
!****

END MODULE ropp_1dvar_types
