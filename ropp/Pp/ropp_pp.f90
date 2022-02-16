! $Id: ropp_pp.f90 4060 2014-02-03 11:13:00Z idculv $

!****t* Interface/Modules
!
! SYNOPSIS
!    use ropp_pp
!    use ropp_pp_types
!    use ropp_pp_constants
!    use ropp_pp_utils
!    use ropp_pp_spline
!    use ropp_pp_MSIS_types
!
!****

!****m* Modules/ropp_pp *
!
! NAME
!    ropp_pp - Interface module for the ROPP pre-processor
!
! SYNOPSIS
!    use ropp_pp
!
! DESCRIPTION
!    This module provides interfaces for all pre-processor routines in the
!    ROPP Preprocessor library.
!
! NOTES
!
! SEE ALSO
!    ropp_pp_constants
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

MODULE ropp_pp

!-------------------------------------------------------------------------------
! 1. Other modules
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_pp_types
  USE ropp_pp_constants
  USE ropp_pp_utils
  USE ropp_pp_spline
  USE ropp_pp_msis_types

!-------------------------------------------------------------------------------
! 2. Compute bending angles - Geometric Optics processing
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_bending_angle_go(time, r_leo, r_gns, r_coc,       &
                                         phase_L1, phase_L2, w_smooth,    &
                                         filter, impact_L1, bangle_L1,    &
                                         impact_L2, bangle_L2)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Relative time
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! centre curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L1  ! L1 excess phase (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L2  ! L2 excess phase (m)
       INTEGER,                  INTENT(in)  :: w_smooth  ! smoothing window
       CHARACTER(len=*),         INTENT(in)  :: filter    ! filter type
       REAL(wp), DIMENSION(:),   INTENT(out) :: impact_L1 ! L1 impact
       REAL(wp), DIMENSION(:),   INTENT(out) :: bangle_L1 ! L1 bending angles
       REAL(wp), DIMENSION(:),   INTENT(out) :: impact_L2 ! L2 impact
       REAL(wp), DIMENSION(:),   INTENT(out) :: bangle_L2 ! L2 bending angles
     END SUBROUTINE ropp_pp_bending_angle_go
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_geometric_optics(r_leo, v_leo, r_gns, v_gns, doppler,  &
                                         impact, bangle)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo     ! LEO velocity
       REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns     ! GPS velocity
       REAL(wp),               INTENT(in)  :: doppler   ! Doppler shift
       REAL(wp),               INTENT(out) :: impact    ! Impact parameter
       REAL(wp),               INTENT(out) :: bangle    ! Bending angle
     END SUBROUTINE ropp_pp_geometric_optics
  END INTERFACE

    INTERFACE
     SUBROUTINE ropp_pp_geometric_optics_adj(r_leo, v_leo, r_gns, v_gns,     &
                                             doppler, impact, bangle,        &
                                             impact_dd, impact_dr,           &
                                             bangle_dd, bangle_dr)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo     ! LEO velocity
       REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns     ! GPS velocity
       REAL(wp),               INTENT(in)  :: doppler   ! Doppler shift
       REAL(wp),               INTENT(out) :: impact    ! Impact parameter
       REAL(wp),               INTENT(out) :: bangle    ! Bending angle
       REAL(wp),               INTENT(out) :: impact_dd ! d(IP)/d(d)
       REAL(wp), DIMENSION(6), INTENT(out) :: impact_dr ! d(IP)/d(r_leo,r_gns)
       REAL(wp),               INTENT(out) :: bangle_dd ! d(BA)/d(d)
       REAL(wp), DIMENSION(6), INTENT(out) :: bangle_dr ! d(BA)/d(r_leo,r_gns)
     END SUBROUTINE ropp_pp_geometric_optics_adj
  END INTERFACE

!-------------------------------------------------------------------------------
! 3. Compute bending angles - Wave Optics processing
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_bending_angle_wo(time, r_leo, r_gns, r_coc, roc,      &
                      phase_L1, phase_L2, snr_L1, snr_L2, w_ls, w_smooth,     &
                      w_low, hmax, filter, opt_DL2, cff, dsh,                 &
                      impact_L1, bangle_L1, ba_sigma_L1,                      &
                      impact_L2, bangle_L2, ba_sigma_L2, diag)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types, ONLY: PPdiag
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! Relative time
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! centre curvature
       REAL(wp),                 INTENT(in)  :: roc       ! radius curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L1  ! L1 excess phase (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_L2  ! L2 excess phase (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: snr_L1    ! L1 amplitude
       REAL(wp), DIMENSION(:),   INTENT(in)  :: snr_L2    ! L2 amplitude
       INTEGER,                  INTENT(in)  :: w_ls      ! smoothing window
       INTEGER,                  INTENT(in)  :: w_smooth  ! smoothing window
       INTEGER,                  INTENT(in)  :: w_low     ! smoothing window
       REAL(wp),                 INTENT(in)  :: hmax      ! max height for WO
       CHARACTER(len=*),         INTENT(in)  :: filter    ! Filter method
       LOGICAL,                  INTENT(in)  :: opt_DL2   ! Degraded L2 flag
       INTEGER,                  INTENT(in)    :: cff     ! Cmplx filter flag
       REAL(wp),                 INTENT(in)    :: dsh     ! Border width (m)
       REAL(wp), DIMENSION(:), INTENT(inout) :: impact_L1 ! L1 impact
       REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_L1 ! L1 bending angles
       REAL(wp), DIMENSION(:), INTENT(out)   :: ba_sigma_L1 ! L1 ba std dev
       REAL(wp), DIMENSION(:), INTENT(inout) :: impact_L2 ! L2 impact
       REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_L2 ! L2 bending angles
       REAL(wp), DIMENSION(:), INTENT(out)   :: ba_sigma_L2 ! L2 ba std dev
       TYPE(ppDiag), OPTIONAL, INTENT(inout) :: diag      ! Additional diags
     END SUBROUTINE ropp_pp_bending_angle_wo
  END INTERFACE

    INTERFACE
     SUBROUTINE ropp_pp_DCT(time,A,DS,r_leo,r_gns,r_coc,roc,w_ls,w_smooth,   &
                            w_low,hmax,filter,opt_DL2,cff,dsh, P, E, EC, diag)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types, ONLY: PPdiag
       REAL(wp), DIMENSION(:),   INTENT(in)    :: time     ! Relative time
       REAL(wp), DIMENSION(:,:), INTENT(in)    :: A        ! Amplitudes
       REAL(wp), DIMENSION(:,:), INTENT(in)    :: DS       ! Excess phase
       REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_leo    ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)    :: r_gns    ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)    :: r_coc    ! Centre curvature
       REAL(wp),                 INTENT(in)    :: roc      ! Radius curvature
       INTEGER,                  INTENT(in)    :: w_ls     ! Large smoothing
       INTEGER,                  INTENT(in)    :: w_smooth ! Smoothing > 7km
       INTEGER,                  INTENT(in)    :: w_low    ! Smoothing < 7km
       REAL(wp),                 INTENT(in)    :: hmax     ! max height for WO
       CHARACTER(len=*),         INTENT(in)    :: filter   ! Filter method
       LOGICAL,                  INTENT(in)    :: opt_DL2  ! Degraded L2 flag
       INTEGER,                  INTENT(in)    :: cff      ! Cmplx filter flag
       REAL(wp),                 INTENT(in)    :: dsh      ! Border width (m)
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: P        ! Impact parameters
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: E        ! Bending angles
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: EC       ! BA covariance
       TYPE(ppDiag), OPTIONAL, INTENT(inout) :: diag       ! Additional diags
     END SUBROUTINE ropp_pp_DCT
  END INTERFACE

!-------------------------------------------------------------------------------
! 4. Preprocessing - Model excess phase calculation and phase correction
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_modelphase(month, lat, lon, time, r_leo, r_gns, r_coc, &
                                   roc, phase_LM, impact_LM, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       INTEGER,                  INTENT(in)  :: month     ! month of year
       REAL(wp),                 INTENT(in)  :: lat       ! latitude
       REAL(wp),                 INTENT(in)  :: lon       ! longitude
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! time of samples
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! centre curvature
       REAL(wp),                 INTENT(in)  :: roc       ! radius curvature
       REAL(wp), DIMENSION(:),   INTENT(out) :: phase_LM  ! model excess phase
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: impact_LM ! model IP
       TYPE(ppConfig),           INTENT(inout) :: config  ! Configuration
     END SUBROUTINE ropp_pp_modelphase
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_bangle2phase(time, r_leo, r_gns, r_coc, impact,        &
                                     bangle, phase, dphi, impact_LM, IL, IU)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time   ! time of samples (s)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo  ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns  ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc  ! centre curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: impact ! impact parameter (m)
       REAL(wp), DIMENSION(:),   INTENT(inout) :: bangle ! bending angle (m)
       REAL(wp), DIMENSION(:),   INTENT(out) :: phase  ! excess phase (m)
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: dphi ! deriv phase(m/s)
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: impact_LM  ! model IP
       INTEGER,                OPTIONAL, INTENT(out) :: IL  ! model lower limit
       INTEGER,                OPTIONAL, INTENT(out) :: IU  ! model upper limit
     END SUBROUTINE ropp_pp_bangle2phase
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_impact2doppler(xleo, vleo, xgns, vgns, impact,         &
                                       doppler, bangle)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(3), INTENT(in)  :: xleo    ! Cartesian LEO coords
       REAL(wp), DIMENSION(3), INTENT(in)  :: vleo    ! Cartesian LEO velocity
       REAL(wp), DIMENSION(3), INTENT(in)  :: xgns    ! Cartesian GPS coords
       REAL(wp), DIMENSION(3), INTENT(in)  :: vgns    ! Cartesian GPS velocity
       REAL(wp),               INTENT(in)  :: impact  ! Impact parameter
       REAL(wp),               INTENT(out) :: doppler ! Doppler frequency shift
       REAL(wp),               INTENT(out) :: bangle  ! Bending angle
     END SUBROUTINE ropp_pp_impact2doppler
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_amplitude_go(time, r_leo, r_gns, r_coc, roc, impact,   &
                                     snr, w_smooth, snr_R)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! time of samples (s)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! centre curvature
       REAL(wp),                 INTENT(in)  :: roc       ! radius curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: impact    ! Impact parameters
       REAL(wp), DIMENSION(:),   INTENT(in)  :: snr       ! Observed amplitude
       INTEGER,                  INTENT(in)  :: w_smooth  ! smoothing window
       REAL(wp), DIMENSION(:),   INTENT(out) :: snr_R     ! Refract amplitude
     END SUBROUTINE ropp_pp_amplitude_go
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_correct_L2(time, r_leo, r_gns, r_coc, roc,             &
                                   impact_LM, phase_LM, phase_L1, phase_L2,    &
                                   snr_L1, snr_L2, LCF, hmid, L2Q)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time      ! time of samples (s)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! Centre curvature
       REAL(wp),                 INTENT(in)  :: roc       ! Radius curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: impact_LM ! Model impact (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM  ! Model excess phase
       REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L1  ! L1 excess phase (m)
       REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L2  ! L2 excess phase (m)
       REAL(wp), DIMENSION(:), INTENT(inout) :: snr_L1    ! L1 amplitude
       REAL(wp), DIMENSION(:), INTENT(inout) :: snr_L2    ! L2 amplitude
       INTEGER,  DIMENSION(:), INTENT(inout) :: lcf       ! Lost carrier flag
       REAL(wp),                 INTENT(in)  :: hmid      ! Border height (m)
       REAL(wp),                 INTENT(out) :: L2Q       ! L2 badness score
     END SUBROUTINE ropp_pp_correct_L2
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_spectra(time, phase_L1, phase_L2, phase_LM,   &
                                impact_LM, config, OutRO, filnam)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       REAL(wp), DIMENSION(:), INTENT(in)     :: time      ! time of samples (s)
       REAL(wp), DIMENSION(:), INTENT(in)     :: phase_L1  ! excess phase L1 (m)
       REAL(wp), DIMENSION(:), INTENT(in)     :: phase_L2  ! excess phase L2 (m)
       REAL(wp), DIMENSION(:), INTENT(in)     :: phase_LM  ! model excess phase (m)
       REAL(wp), DIMENSION(:), INTENT(in)     :: impact_LM ! model impact param (m)
       TYPE(PPConfig),         INTENT(in)     :: config    ! Configuration options
       LOGICAL,          OPTIONAL, INTENT(in) :: OutRO     ! Flag to output RO spectra
       CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: filnam    ! Output file name root
     END SUBROUTINE ropp_pp_spectra
  END INTERFACE

!-------------------------------------------------------------------------------
! 5. Preprocessing - Open loop preprocessing
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_openloop(time, phase_L1, phase_L2, phase_LM,           &
                                 r_leo, r_gns, r_coc, LCF)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: time      ! time of samples (s)
       REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L1  ! excess phase L1 (m)
       REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L2  ! excess phase L2 (m)
       REAL(wp), DIMENSION(:), INTENT(in)    :: phase_LM  ! model ex phase (m)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coords (ECF)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coords (ECF)
       REAL(wp), DIMENSION(:), INTENT(in)    :: r_coc     ! centre curv (ECF)
       INTEGER,  DIMENSION(:), INTENT(inout) :: LCF       ! lost carrier flag
     END SUBROUTINE ropp_pp_openloop
  END INTERFACE

!-------------------------------------------------------------------------------
! 6. Radioholographic filtering and analysis
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_radioholographic_filter(time, r_leo, r_gns, r_coc,     &
                                                roc, phase_LM, phase, snr)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time     ! time of samples (s)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo    ! LEO coordinates (m)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns    ! GPS coordinates (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc    ! centre curvature (m)
       REAL(wp),                 INTENT(in)  :: roc      ! radius of curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM ! model phase (m)
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: phase  ! L1,L2 excess phase
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: snr    ! L1,L2 amplitude
     END SUBROUTINE ropp_pp_radioholographic_filter
  END INTERFACE

    INTERFACE
     SUBROUTINE ropp_pp_radiooptic_analysis(time, r_leo, r_gns, r_coc, roc,    &
                                          phase_LM, phase, snr, PA, PD, OutRO, filnam)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: time     ! time of samples (s)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo    ! LEO coordinates (m)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns    ! GPS coordinates (m)
       REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc    ! centre curvature (m)
       REAL(wp),                 INTENT(in)  :: roc      ! radius of curvature
       REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM ! model phase (m)
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: phase    ! L1,L2 excess phase
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: snr      ! L1,L2 amplitude
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: PA  ! Average impact
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: PD  ! RMS impact
       LOGICAL,  OPTIONAL,       INTENT(in)  :: OutRO    ! Flag to output spectra
       CHARACTER(LEN=*), OPTIONAL,INTENT(in) :: filnam   ! Output file name root
     END SUBROUTINE ropp_pp_radiooptic_analysis
  END INTERFACE

!-------------------------------------------------------------------------------
! 7. Satellite location processing
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,   &
                                             xgns, vgns, abl, abg)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: time    ! time of samples
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo   ! LEO position
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns   ! GPS position
       REAL(wp), DIMENSION(:,:), INTENT(out) :: xleo    ! LEO position
       REAL(wp), DIMENSION(:,:), INTENT(out) :: vleo    ! LEO velocity
       REAL(wp), DIMENSION(:,:), INTENT(out) :: xgns    ! GPS position
       REAL(wp), DIMENSION(:,:), INTENT(out) :: vgns    ! GPS velocity
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abl ! LEO regrsn coeff
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abg ! GPS regrsn coeff
     END SUBROUTINE ropp_pp_satellite_velocities
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_interpolate_trajectory(time,cleo,cgns,r_coc,t_init,    &
                                               xleo, vleo, xgns, vgns, theta)
        USE typesizes, ONLY: wp => EightByteReal
        REAL(wp), DIMENSION(:),   INTENT(in)  :: time   ! time of samples (s)
        REAL(wp), DIMENSION(:,:), INTENT(in)  :: cleo   ! LEO regression coeff
        REAL(wp), DIMENSION(:,:), INTENT(in)  :: cgns   ! GPS regression coeff
        REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc  ! centre of curvature
        REAL(wp),                 INTENT(in)  :: t_init ! interpolation time
        REAL(wp), DIMENSION(:),   INTENT(out) :: xleo   ! LEO coordinates
        REAL(wp), DIMENSION(:),   INTENT(out) :: vleo   ! LEO velocity
        REAL(wp), DIMENSION(:),   INTENT(out) :: xgns   ! GPS coordinates
        REAL(wp), DIMENSION(:),   INTENT(out) :: vgns   ! GPS velocity
        REAL(wp), OPTIONAL,       INTENT(out) :: theta  ! GPS -> LEO angle
      END SUBROUTINE ropp_pp_interpolate_trajectory
END INTERFACE

!-------------------------------------------------------------------------------
! 8. Inverse Abel transform
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_invert_exp(impact, bangle, nr, refrac)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact  ! impact parameter
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle  ! bending angle obs
       REAL(wp), DIMENSION(:), INTENT(in)  :: nr      ! x=nr product
       REAL(wp), DIMENSION(:), INTENT(out) :: refrac  ! refractivity
     END SUBROUTINE ropp_pp_invert_exp
  END INTERFACE

  INTERFACE ropp_pp_invert
     SUBROUTINE ropp_pp_invert_lin(impact, bangle, nr, refrac, scale)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact  ! impact parameter
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle  ! bending angle obs
       REAL(wp), DIMENSION(:), INTENT(in)  :: nr      ! x=nr product
       REAL(wp), OPTIONAL,     INTENT(in)  :: scale   ! vertical scale
       REAL(wp), DIMENSION(:), INTENT(out) :: refrac  ! refractivity
     END SUBROUTINE ropp_pp_invert_lin
  END INTERFACE

!-------------------------------------------------------------------------------
! 9. Forward Abel transform
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_abel_exp(nr, refrac, impact, bangle)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: nr         ! x=nr product
       REAL(wp), DIMENSION(:), INTENT(in)  :: refrac     ! refractivity
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact     ! impact parameter
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle     ! bending angle
     END SUBROUTINE ropp_pp_abel_exp
  END INTERFACE

  INTERFACE ropp_pp_abel
     SUBROUTINE ropp_pp_abel_lin(nr, refrac, impact, bangle, dln, scale)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: nr         ! x=nr product
       REAL(wp), DIMENSION(:), INTENT(in)  :: refrac     ! refractivity
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact     ! impact parameter
       REAL(wp), OPTIONAL, DIMENSION(:), INTENT(in) :: dln ! ref index gradient
       REAL(wp), OPTIONAL,     INTENT(in)  :: scale      ! scale factor
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle     ! bending angle
     END SUBROUTINE ropp_pp_abel_lin
  END INTERFACE

!-------------------------------------------------------------------------------
! 10. Ionospheric correction
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_ionospheric_correction(impact_L1, bangle_L1,      &
          impact_L2, bangle_L2, impact_LM, bangle_LM, config,             &
          impact_LC, bangle_LC, diag)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L1  ! L1 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L1  ! L1 bending angles
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L2  ! L2 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L2  ! L2 bending angles
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_LM  ! Model bending angles
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_LM  ! Model bending angles
       TYPE(ppConfig),       INTENT(inout) :: config     ! Configuration
       REAL(wp), DIMENSION(:), INTENT(out) :: impact_LC  ! LC impact parameters
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_LC  ! Corrected bangle
       TYPE(ppDiag),         INTENT(inout) :: diag       ! Output diagnostics
     END SUBROUTINE ropp_pp_ionospheric_correction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_linear_combination(impact_L1, bangle_L1, impact_L2,    &
                                           bangle_L2, impact_LC, bangle_LC)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L1  ! L1 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L1  ! L1 bending angles
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L2  ! L2 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L2  ! L2 bending angles
       REAL(wp), DIMENSION(:), INTENT(out) :: impact_LC  ! LC impact parameters
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_LC  ! Corrected bangle
     END SUBROUTINE ropp_pp_linear_combination
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_merge_profile(impact_L1,bangle_L1,impact_L2,bangle_L2, &
                                      impact_I1,bangle_I1,impact_I2,bangle_I2, &
                                      Pmin, Pmax)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L1  ! L1 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L1  ! L1 bending angles
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact_L2  ! L2 impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle_L2  ! L2 bending angles
       REAL(wp), DIMENSION(:), INTENT(out) :: impact_I1  ! Interpolated impact
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_I1  ! Interpolated L1 ba
       REAL(wp), DIMENSION(:), INTENT(out) :: impact_I2  ! Interpolated impact
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_I2  ! Interpolated L2 ba
       REAL(wp), OPTIONAL, INTENT(in)      :: Pmin       ! Minimum impact
       REAL(wp), OPTIONAL, INTENT(in)      :: Pmax       ! Maximum impact
     END SUBROUTINE ropp_pp_merge_profile
  END INTERFACE

!-------------------------------------------------------------------------------
! 11. Fitting to model bending angle profile
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_invert_refraction(mfile, month, lat, lon, impact,      &
                                          bangle, geop, refrac, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       CHARACTER(len=*),    INTENT(inout)  :: mfile     ! Coefficients file
       INTEGER,                INTENT(in)  :: month     ! Month of year
       REAL(wp),               INTENT(in)  :: lat       ! Latitude
       REAL(wp),               INTENT(in)  :: lon       ! Longitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact    ! Impact parameter
       REAL(wp), DIMENSION(:), INTENT(in)  :: bangle    ! Bending angle
       REAL(wp), DIMENSION(:), INTENT(out) :: geop      ! Geopotential height
       REAL(wp), DIMENSION(:), INTENT(out) :: refrac    ! Refractivity
       TYPE(PPConfig),         INTENT(in)  :: config    ! Configuration
     END SUBROUTINE ropp_pp_invert_refraction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_model_refraction(mfile, month, lat, lon, impact,     &
                                         bangle_MSIS, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       CHARACTER(len=*),    INTENT(inout)  :: mfile        ! Coefficients file
       INTEGER,                INTENT(in)  :: month        ! Month of year
       REAL(wp),               INTENT(in)  :: lat          ! Latitude
       REAL(wp),               INTENT(in)  :: lon          ! Longitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact       ! Impact parameter
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_MSIS  ! Model bending angle
       TYPE(PPConfig),         INTENT(in)  :: config       ! Configuration
     END SUBROUTINE ropp_pp_model_refraction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_bg_refraction(bfile, month, lat, lon, impact,     &
                                         bangle_BG, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       CHARACTER(len=*),       INTENT(in)  :: bfile        ! Background file
       INTEGER,                INTENT(in)  :: month        ! Month of year
       REAL(wp),               INTENT(in)  :: lat          ! Latitude
       REAL(wp),               INTENT(in)  :: lon          ! Longitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact       ! Impact parameter
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_BG    ! BG bending angle
       TYPE(PPConfig),         INTENT(in)  :: config       ! Configuration
     END SUBROUTINE ropp_pp_bg_refraction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_search_model_refraction(mfile,in_month,in_lat,in_lon, &
                                                in_impact, in_bangle, &
                                                impact, bangle_MSIS, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       CHARACTER(len=*),    INTENT(inout)  :: mfile        ! Coefficients file
       INTEGER,                INTENT(in)  :: in_month     ! Month of year
       REAL(wp),               INTENT(in)  :: in_lat       ! Latitude
       REAL(wp),               INTENT(in)  :: in_lon       ! Longitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: in_impact    ! Input impact
       REAL(wp), DIMENSION(:), INTENT(in)  :: in_bangle    ! Input bangle
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact       ! Impact parameter
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle_MSIS  ! Model bangle
       TYPE(PPConfig),         INTENT(in)  :: config       ! Configuration
     END SUBROUTINE ropp_pp_search_model_refraction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_fit_model_refraction(impact_LC, bangle_LC,             &
                                             impact_model, bangle_model, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       REAL(wp), DIMENSION(:), INTENT(in)    :: impact_LC    ! Impact parameters
       REAL(wp), DIMENSION(:), INTENT(in)    :: bangle_LC    ! Bending angle
       REAL(wp), DIMENSION(:), INTENT(in)    :: impact_model ! Model impact
       REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_model ! Model bangle
       TYPE(PPConfig),         INTENT(in)    :: config       ! Configuration
     END SUBROUTINE ropp_pp_fit_model_refraction
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_smooth_profile(impact, bangle, smooth, config)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_pp_types
       REAL(wp), DIMENSION(:), INTENT(in)    :: impact    ! Impact parameter
       REAL(wp), DIMENSION(:), INTENT(in)    :: bangle    ! Bending angle
       REAL(wp), DIMENSION(:), INTENT(out)   :: smooth    ! Smoothed bangle
       TYPE(PPConfig),         INTENT(inout) :: config    ! Configuration
     END SUBROUTINE ropp_pp_smooth_profile
  END INTERFACE

!-------------------------------------------------------------------------------
! 12. Vertical interpolation
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_interpol
     SUBROUTINE ropp_pp_interpol_scl(x, newx, array, interp, Cext)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp),               INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp),               INTENT(out) :: interp
       LOGICAL,  OPTIONAL,     INTENT(in)  :: Cext
     END SUBROUTINE ropp_pp_interpol_scl
     SUBROUTINE ropp_pp_interpol_arr(x, newx, array, interp, Cext)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(out) :: interp
       LOGICAL,  OPTIONAL,     INTENT(in)  :: Cext
     END SUBROUTINE ropp_pp_interpol_arr
          SUBROUTINE ropp_pp_interpol_int(x, newx, array, interp)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp),               INTENT(in)  :: newx
       INTEGER,  DIMENSION(:), INTENT(in)  :: array
       INTEGER,                INTENT(out) :: interp
     END SUBROUTINE ropp_pp_interpol_int
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_interpol_log(x, newx, array, interp)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(out) :: interp
     END SUBROUTINE ropp_pp_interpol_log
  END INTERFACE

!-------------------------------------------------------------------------------
! 13. Monotonization
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_uni_monotonous(x, d)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(inout) :: x
       INTEGER, OPTIONAL                     :: d
     END SUBROUTINE ropp_pp_uni_monotonous
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_monotonous(x, d)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(inout) :: x
       INTEGER, OPTIONAL                     :: d
     END SUBROUTINE ropp_pp_monotonous
  END INTERFACE

!-------------------------------------------------------------------------------
! 14. Read Config file
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_read_config(file, config)
       USE ropp_pp_types
       CHARACTER(len=*), INTENT(in)  :: file
       TYPE(PPConfig), INTENT(inout) :: config
     END SUBROUTINE ropp_pp_read_config
  END INTERFACE

!-------------------------------------------------------------------------------
! 15. Dry temperature and pressure calculation
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_tdry(lat, alt, refrac, shum, t_dry, p_dry, Zmax)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp),               INTENT(in)  :: lat     ! Latitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: alt     ! Altitude
       REAL(wp), DIMENSION(:), INTENT(in)  :: refrac  ! Refractivity
       REAL(wp), DIMENSION(:), INTENT(in)  :: shum    ! Specific humidity
       REAL(wp), DIMENSION(:), INTENT(out) :: t_dry   ! Dry temperature
       REAL(wp), DIMENSION(:), INTENT(out) :: p_dry   ! Dry pressure
       REAL(wp), OPTIONAL,     INTENT(in)  :: Zmax    ! Maximum altitude
     END SUBROUTINE ropp_pp_tdry
  END INTERFACE

!-------------------------------------------------------------------------------
! 16. FFT
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_FFT
     SUBROUTINE ropp_pp_FFT_real(data, isign)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(INOUT) :: data       ! Complex signal
       INTEGER,                   INTENT(IN)    :: isign      ! FFT direction
     END SUBROUTINE ropp_pp_FFT_real
     SUBROUTINE ropp_pp_FFT_complex(data, isign)
       USE typesizes, ONLY: wp => EightByteReal
       COMPLEX(wp), DIMENSION(:), INTENT(INOUT) :: data       ! Complex signal
       INTEGER,                   INTENT(IN)    :: isign      ! FFT direction
     END SUBROUTINE ropp_pp_FFT_complex
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_fourier_filter(data, window)
       USE typesizes, ONLY: wp => EightByteReal
       COMPLEX(wp), DIMENSION(:), INTENT(INOUT) :: data       ! Complex signal
       REAL(wp),                  INTENT(IN)    :: window     ! Window width
     END SUBROUTINE ropp_pp_fourier_filter
  END INTERFACE

!-------------------------------------------------------------------------------
! 17. Signal Filtering (optimal estimation)
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_filter
     SUBROUTINE ropp_pp_filter_1d(DT, S, W, ND, FS, DS)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp),                 INTENT(in)  :: DT ! Time step
       REAL(wp), DIMENSION(:),   INTENT(in)  :: S  ! Signal samples [time]
       INTEGER,                  INTENT(in)  :: W  ! Window width [points]
       INTEGER,                  INTENT(in)  :: ND ! Differentiation points
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: FS ! Filtered signal
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal derivative
     END SUBROUTINE ropp_pp_filter_1d
     SUBROUTINE ropp_pp_filter_2d(DT, S, W, ND, FS, DS)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp),                 INTENT(in)  :: DT ! Time step
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: S  ! Signal samples [ch,time]
       INTEGER,                  INTENT(in)  :: W  ! Window width [points]
       INTEGER,                  INTENT(in)  :: ND ! Differentiation points
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: FS ! Filtered signal
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: DS ! Signal derivative
      END SUBROUTINE ropp_pp_filter_2d
  END INTERFACE

!-------------------------------------------------------------------------------
! 18. Sliding Polynomial Filtering
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_sliding_polynomial
     SUBROUTINE ropp_pp_sliding_poly_1d(t, s, w, np, fs, ds)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in) :: T  ! Time
       REAL(wp), DIMENSION(:),   INTENT(in) :: S  ! Signal samples (time)
       INTEGER,                  INTENT(in)  :: W  ! Window width [samples]
       INTEGER,                  INTENT(in)  :: NP ! Polynomial degree
       REAL(wp), DIMENSION(:),   INTENT(out) :: FS ! Filtered signal
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal deriv
     END SUBROUTINE ropp_pp_sliding_poly_1d
     SUBROUTINE ropp_pp_sliding_poly_2d(t, s, w, np, fs, ds)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:),   INTENT(in)  :: T  ! Time
       REAL(wp), DIMENSION(:,:), INTENT(in)  :: S  ! Signal [ch,t]
       INTEGER,                  INTENT(in)  :: W  ! Window width [samples]
       INTEGER,                  INTENT(in)  :: NP ! Polynomial degree
       REAL(wp), DIMENSION(:,:), INTENT(out) :: FS ! Filtered signal
       REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: DS ! Signal deriv

     END SUBROUTINE ropp_pp_sliding_poly_2d
     SUBROUTINE ropp_pp_sliding_poly_vec1d(t, s, w, np, fs, ds)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in) :: T  ! Time
       REAL(wp), DIMENSION(:), INTENT(in) :: S  ! Signal samples (time)
       INTEGER,  DIMENSION(:), INTENT(in) :: W  ! Window width (time)
       INTEGER,                INTENT(in) :: NP ! Polynomial degree
       REAL(wp), DIMENSION(:),           INTENT(out) :: FS ! Filtered signal
       REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS ! Signal derivative
     END SUBROUTINE ropp_pp_sliding_poly_vec1d
  END INTERFACE

!-------------------------------------------------------------------------------
! 19. Tropopause height
!-------------------------------------------------------------------------------

  INTERFACE
     SUBROUTINE ropp_pp_tph_bangle(ro_data, diag)
       USE ropp_io_types, ONLY: ROprof
       TYPE(ROprof), INTENT(INOUT)               :: ro_data
       LOGICAL, OPTIONAL, INTENT(IN)             :: diag
     END SUBROUTINE ropp_pp_tph_bangle
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_tph_refrac(ro_data, diag)
       USE ropp_io_types, ONLY: ROprof
       TYPE(ROprof), INTENT(INOUT)               :: ro_data
       LOGICAL, OPTIONAL, INTENT(IN)             :: diag
     END SUBROUTINE ropp_pp_tph_refrac
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_tph_tdry(ro_data, diag)
       USE ropp_io_types, ONLY: ROprof
       TYPE(ROprof), INTENT(INOUT)               :: ro_data
       LOGICAL, OPTIONAL, INTENT(IN)             :: diag
     END SUBROUTINE ropp_pp_tph_tdry
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_tph_temp(ro_data, diag)
       USE ropp_io_types, ONLY: ROprof
       TYPE(ROprof), INTENT(INOUT)               :: ro_data
       LOGICAL, OPTIONAL, INTENT(IN)             :: diag
     END SUBROUTINE ropp_pp_tph_temp
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_cov_transform(y, x, a, cov)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(IN)        :: y, x
       REAL(wp),               INTENT(IN)        :: a
       REAL(wp), DIMENSION(:), INTENT(INOUT)     :: cov
     END SUBROUTINE ropp_pp_cov_transform
  END INTERFACE

!-------------------------------------------------------------------------------
! 20. Common utilities
!-------------------------------------------------------------------------------

  INTERFACE ropp_pp_version
    FUNCTION ropp_pp_version() RESULT (version)
      CHARACTER (LEN=40) :: version
    END FUNCTION ropp_pp_version
  END INTERFACE ropp_pp_version

END MODULE ropp_pp
