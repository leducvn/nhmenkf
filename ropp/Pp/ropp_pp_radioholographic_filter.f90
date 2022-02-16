! $Id: ropp_pp_radioholographic_filter.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_radioholographic_filter(time, r_leo, r_gns, r_coc, roc,    &
                                           phase_LM, phase, snr)

!****s* Preprocessing/ropp_pp_radioholographic_filter *
!
! NAME
!    ropp_pp_radioholographic_filter - Radioholographic filtering of radio
!                                      occultation data
!                   
! SYNOPSIS
!    call ropp_pp_radioholographic_filter(time,r_leo,r_gns,r_coc,roc,phase_LM,
!                                         phase, snr)
! 
! DESCRIPTION
!     1) Multiplication of the radio occulutation data with a reference signal 
!     2) Fourier filter combined signal
!     3) Removal of reference signal
!
! INPUTS
!    real(wp), dimension(:)   :: time     ! relative time of samples (s)
!    real(wp), dimension(:,:) :: r_leo    ! cartesian LEO coordinates (m)
!    real(wp), dimension(:,:) :: r_gns    ! cartesian GPS coordinates (m)
!    real(wp), dimension(:)   :: r_coc    ! cartesian centre curvature (m)
!    real(wp)                 :: roc      ! radius of curvature (m)
!    real(wp), dimension(:)   :: phase_LM ! model excess phase (m)
!    real(wp), dimension(:,:) :: phase    ! L1 and L2 excess phase (m)
!    real(wp), dimension(:,:) :: snr      ! L1 and L2 amplitude
!
! OUTPUT
!    real(wp), dimension(:,:) :: phase    ! L1 and L2 excess phase (m)
!    real(wp), dimension(:,:) :: snr      ! L1 and L2 amplitude
!
! NOTES
!
! REFERENCES
!   Gorbunov M.E., Lauritsen K.B., Rhodin A., Tomassini M. and Kornblueh L. 
!   2006
!   Radio holographic filtering, error estimation, and quality control of 
!   radio occultation data
!   Journal of Geophysical Research (111) D10105
!
!   Gorbunov M.E., Lauritsen K.B., Rodin A., Tomassini M., Kornblueh L. 2005 
!   Analysis of the CHAMP experimental data on radio-occultation sounding of
!   the Earth's atmosphere.
!   Izvestiya Atmospheric and Oceanic Physics (41) 726-740.

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
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
! USE ropp_pp, not_this => ropp_pp_radioholographic_filter
  USE ropp_pp
  USE ropp_utils, ONLY: vector_angle
  
  IMPLICIT NONE

  REAL(wp), DIMENSION(:),   INTENT(in)  :: time     ! Time of samples (s)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo    ! Cartesian LEO coords (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns    ! Cartesian GPS coords (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc    ! Centre curvature (m)
  REAL(wp),                 INTENT(in)  :: roc      ! Radius of curvature (m)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: phase_LM ! Model excess phase (m)
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: phase  ! L1,L2 excess phase (m)
  REAL(wp), DIMENSION(:,:), INTENT(inout) :: snr    ! L1 and L2 amplitude
  
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xleo     ! LEO position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: xgns     ! GPS position by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vleo     ! LEO velocity by regression
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: vgns     ! GPS velocity by regression

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: k        ! Wave numbers
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: PA       ! Average impact parameter
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DA       ! Smoothed doppler
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: EA       ! Smoothed bending angle
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: SA       ! Smoothed excess phase
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: TH       ! Grid for complex field
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: PH       ! Phase of complex field
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: AH       ! Amplitude of complex field
  COMPLEX(wp), DIMENSION(:), ALLOCATABLE :: UH      ! Complex field 
  REAL(wp),    DIMENSION(:), ALLOCATABLE :: DS      ! phase - SA

  INTEGER             :: npoints      ! Number of data points
  INTEGER             :: nc           ! Number of channels
  INTEGER             :: i            ! Array index
  INTEGER             :: ic           ! Channel index
  INTEGER             :: nh           ! Number of UH samples
  REAL(wp)            :: Tmin         ! Minimum time
  REAL(wp)            :: Tmax         ! Maximum time
  REAL(wp)            :: ThetaMin     ! Minimum satellite-to-satellite angle
  REAL(wp)            :: ThetaMax     ! Maximum satellite-to-satellite angle
  REAL(wp)            :: ThetaDot     ! Angular velocity
  REAL(wp)            :: DTW          ! Time filter width
  REAL(wp)            :: W            ! Windows width

  REAL(wp), PARAMETER :: DES = 0.0025_wp ! Spectral width of bending angle (rad)
  REAL(wp), PARAMETER :: DPW = 500.0_wp  ! Filter width of impact parameter (m)
  COMPLEX(wp), PARAMETER :: Ci = (0.0_wp, 1.0_wp)  ! I = Sqrt(-1) 

!-------------------------------------------------------------------------------
! 2. Computation of occultation geometry
!-------------------------------------------------------------------------------

  ! 2.1. Determination of number of data
  
  npoints = SIZE(time)
  nc = SIZE(snr,1)
  
  ! 2.2 Calculation of GPS and LEO velocities 
  
  ALLOCATE(xgns(npoints,3))
  ALLOCATE(xleo(npoints,3)) 
  ALLOCATE(vgns(npoints,3))
  ALLOCATE(vleo(npoints,3))

  CALL ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo, xgns, vgns)

!-------------------------------------------------------------------------------
! 3. Computation of reference signal
!-------------------------------------------------------------------------------

  ! 3.1 Calculation of smoothed impact parameter

  ALLOCATE(PA(nc,npoints))
  
  CALL ropp_pp_radiooptic_analysis(time, r_leo, r_gns, r_coc, roc, phase_LM,  &
                                   phase, snr, PA)

  ! 3.2 Calculation of smoothed Doppler and bending angle

  ALLOCATE(DA(nc,npoints))
  ALLOCATE(EA(nc,npoints))
  ALLOCATE(SA(nc,npoints))

  DO ic=1,nc
     DO i=1,npoints
        CALL ropp_pp_impact2doppler(xleo(i,:)-r_coc, vleo(i,:),     &
                                    xgns(i,:)-r_coc, vgns(i,:),     &
                                    PA(ic,i), DA(ic,i), EA(ic,i))
     ENDDO
  ENDDO

  ! 3.3 Calculation of reference phase

  DO ic=1,nc
     
     SA(ic,1) = SQRT( SUM((xgns(1,:)-xleo(1,:))**2 ))
     
     DO i=2,npoints
        SA(ic,i) = SA(ic,i-1) -     &
                     C_Light*(DA(ic,i-1)+DA(ic,i))*(time(i)-time(i-1))/2.0_wp
     ENDDO
     
     DO i=1,npoints
        SA(ic,i) = SA(ic,i) - SQRT( SUM((xgns(i,:)-xleo(i,:))**2) )
     ENDDO
     
  ENDDO

  DEALLOCATE(PA) 
  DEALLOCATE(DA)
  DEALLOCATE(EA)

!-------------------------------------------------------------------------------
! 4. Filtering
!-------------------------------------------------------------------------------

  ALLOCATE(DS(npoints))
  ALLOCATE(k(nc))

  DO ic=1,nc

     IF (nc == 1) THEN
        k(ic) = 2.0_wp*Pi*f_L2/C_light
     ELSE
        IF (ic == 1) k = 2.0_wp*Pi*f_L1/C_Light
        IF (ic == 2) k = 2.0_wp*Pi*f_L2/C_Light
     ENDIF
     
     DS(:) = phase(ic,:) - SA(ic,:)
     
     ! 4.1 Determination of resolution
     
     nh = MAX(npoints, CEILING(npoints*k(ic)*    &
             ABS(MAXVAL(DS(1:npoints-1) - DS(2:npoints)))/(Pi/2)))
     nh = ropp_pp_nearest_power2(nh)

     ! 4.2. Sampling complex field
     
     ALLOCATE(TH(nh))
     ALLOCATE(AH(nh))
     ALLOCATE(PH(nh))
     ALLOCATE(UH(nh))

     Tmax = MAXVAL(time)
     Tmin = MINVAL(time)
     
     DO i=1,nh
        TH(i) = (Tmin*(nh-i) + Tmax*(i-1))/(nh-1)
     ENDDO
     
     CALL ropp_pp_fast_interpol(time, TH, snr(ic,:), AH)

     CALL ropp_pp_fast_interpol(time, TH, DS, PH)

     PH(:) = k(ic)*PH(:)

     UH(:) = AH(:)*EXP(Ci*(MODULO(PH(:), 2.0_wp*Pi))) 

     ! 4.3 Determination of filter width

     ThetaMin = Vector_Angle(xleo(1,:) - r_coc, xgns(1,:) - r_coc)
     ThetaMax = Vector_Angle(xleo(npoints,:) - r_coc, xgns(npoints,:) - r_coc)
     ThetaDot = ABS(ThetaMax - ThetaMin)/(Tmax - Tmin)
     
     DTW      = Pi/(k(ic)*ThetaDot*DPW)
     W        = DTW*(nh-1)/(Tmax - Tmin)

     ! 4.4 Filtering

     CALL ropp_pp_fourier_filter(UH, W)

     AH(:) = ABS(UH(:))

     WHERE (UH(:) /= 0.0_wp)
        PH(:) = ATAN2(AIMAG(UH(:)),REAL(UH(:)))
     ELSEWHERE
        PH(:) = 0.0_wp
     END WHERE
     
     CALL ropp_pp_accumulate_phase(ph(:))

     CALL ropp_pp_fast_interpol(TH, time, AH, snr(ic,:))
     
     CALL ropp_pp_fast_interpol(TH, time, PH, DS)
     DS(:) = DS(:)/k(ic)

     phase(ic,:) = SA(IC,:) + DS(:)

     ! 4.5 Memory deallocation
     
     DEALLOCATE(TH)
     DEALLOCATE(AH)
     DEALLOCATE(PH)
     DEALLOCATE(UH)
    
  ENDDO

!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(xgns)
  DEALLOCATE(xleo)
  DEALLOCATE(vgns)
  DEALLOCATE(vleo)
  
  DEALLOCATE(SA)
  DEALLOCATE(DS)
  DEALLOCATE(k)
  
CONTAINS

!-------------------------------------------------------------------------------
! 6. Transformation of phase to accumulated phase
!-------------------------------------------------------------------------------

  SUBROUTINE ropp_pp_accumulate_phase(ph)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp, ONLY: Pi
    IMPLICIT NONE
    
    ! 6.1 Declarations
    
    REAL(wp), DIMENSION(:), INTENT(inout)  :: ph   ! phase -> accumulated phase
    INTEGER                                :: i
    
    ! 6.2 Add +-2Pi where phase jumps from +-Pi to -+Pi

    DO i=2,SIZE(ph)
       ph(i) = ph(i-1) + MODULO(ph(i)-ph(i-1)+Pi, 2.0_wp*Pi) - Pi
    ENDDO

  END SUBROUTINE ropp_pp_accumulate_phase
  
!-------------------------------------------------------------------------------
! 7. Fast linear interpolation
!-------------------------------------------------------------------------------

  SUBROUTINE ropp_pp_fast_interpol(y1, y2, f1, f2)

    USE typesizes, ONLY: wp => EightByteReal
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(in)  :: y1
    REAL(wp), DIMENSION(:), INTENT(in)  :: y2
    REAL(wp), DIMENSION(:), INTENT(in)  :: f1
    REAL(wp), DIMENSION(:), INTENT(out) :: f2

    INTEGER :: n1, n2, i1, i2
    REAL(wp) :: dy1

    n1 = size(y1)
    n2 = size(y2)
    dy1 = y1(2) - y1(1)

    DO i2 = 1, n2
      i1 = 1 + FLOOR((y2(i2) - y1(1))/DY1)
      i1 = MAX(1, Min(N1-1, i1))
      f2(i2) = (f1(i1)*(y1(i1+1)-y2(i2)) + f1(i1+1)*(y2(i2)-y1(i1)))/dy1
    ENDDO

  END SUBROUTINE ropp_pp_fast_interpol

END SUBROUTINE ropp_pp_radioholographic_filter

