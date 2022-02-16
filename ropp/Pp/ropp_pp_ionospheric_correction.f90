! $Id: ropp_pp_ionospheric_correction.f90 4452 2015-01-29 14:42:02Z idculv $

SUBROUTINE ropp_pp_ionospheric_correction(impact_L1, ba_L1, impact_L2, ba_L2, &
                          impact_LM, ba_LM, config, impact_LC, ba_LC, diag_out)

!****s* IonosphericCorrection/ropp_pp_ionospheric_correction *
!
! NAME
!    ropp_pp_ionospheric_correction - 
!                   Calculate neutral and ionospheric bending angle 
!                   profile on L1 impact heights from L1/L2 bending angles 
!                   
! SYNOPSIS
!    call ropp_pp_ionospheric_correction(impact_L1, ba_L1, 
!                                        impact_L2, ba_L2,
!                                        impact_LM, ba_LM,
!                                        config, impact_LC, ba_LC, diag_out)
! 
! DESCRIPTION
!    This routine calculates bending angles at a given set of impact parameters
!    from vertical profiles of bending angles at the two measurement 
!    frequencies (channels) L1 and L2.
!
! INPUTS
!    real(wp), dimension(:) :: impact_L1   ! Impact parameters of channel L1 (m)
!    real(wp), dimension(:) :: ba_L1       ! Bending angles for channel L1 (rad)
!    real(wp), dimension(:) :: impact_L2   ! Impact parameters of channel L2 (m)
!    real(wp), dimension(:) :: ba_L2       ! Bending angles for channel L2 (rad)
!    real(wp), dimension(:) :: impact_LM   ! Model impact parameters (m)
!    real(wp), dimension(:) :: ba_LM       ! Model bending angles (rad)
!    type(ppConfig)         :: config      ! Configuration parameters
!
! OUTPUT
!    real(wp), dimension(:) :: impact_LC   ! Impact parameters of channel L1
!    real(wp), dimension(:) :: ba_LC       ! Corrected bending angles 
!    type(ppDiag)           :: diag_out    ! Output diagnostics
! 
! NOTES
!    Method:
!        1. Calculation of strongly smoothed ionospheric signal
!           (using the external scale). Further deviations from
!           this signal are calculated.
!        2. Estimation of ionospheric signal and noise covariances
!           using the highest part (> 50 km) of the occultation.
!        3. Calculation of relative mean deviation of neutral refraction
!           from the model refraction using signal at heights 12-35 km.
!        4. Optimal linear combination for the same impact parameters using the
!           covariances.
!
! REFERENCES
!     M.E. Gorbunov
!     Ionospheric correction and statistical optimization of radio 
!     occultation data
!     Radio Science, 37(5), 1084
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

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
! USE ropp_pp, not_this => ropp_pp_ionospheric_correction
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPConfig, PPDiag

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)    :: impact_L1 ! L1 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)    :: ba_L1     ! L1 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(in)    :: impact_L2 ! L2 impact parameters (m)
  REAL(wp), DIMENSION(:), INTENT(in)    :: ba_L2     ! L2 bending angles (rad)
  REAL(wp), DIMENSION(:), INTENT(in)    :: impact_LM ! Model impact parameters 
  REAL(wp), DIMENSION(:), INTENT(in)    :: ba_LM     ! Model bending angles 
  TYPE(ppConfig),         INTENT(inout) :: config    ! Configuration parameters
  REAL(wp), DIMENSION(:), INTENT(out)   :: impact_LC ! Output impact parameters
  REAL(wp), DIMENSION(:), INTENT(out)   :: ba_LC     ! Corrected bending angles
  TYPE(ppDiag),           INTENT(inout) :: diag_out  ! Output diagnostics

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: impact_LH ! Hi-res impact grid
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_L1H    ! L1 bangle on impact_LH
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_L2H    ! L2 bangle on impact_LH
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_LMH    ! LM bangle on impact_LH
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ba_diff   ! bangle differences
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ba_low    ! low-freq ba_ob-ba_model
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_is     ! smoothed ionospheric ba
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_ion    ! ionospheric bangle
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_neut   ! neutral bangle
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_nfilt  ! filtered neutral bangle
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_ifilt  ! filtered ionospheric ba
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: c_nfilt   ! error covariance ba_nfilt
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: c_ifilt   ! error covariance ba_ifilt
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: wlch      ! weight of filtered data in LC
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: ba_LI     ! ionospheric bangle ba_LI
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: err_LI    ! error covariance ba_LI
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: err_LC    ! error covariance ba_LC
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: WLC       ! weight of data in LC
  LOGICAL,  DIMENSION(:),   ALLOCATABLE :: m_diff    ! array mask
    
  REAL(wp), DIMENSION(2)   :: cin    ! covariance of ionospheric noise
  REAL(wp), DIMENSION(2)   :: cis    ! covariance of ionospheric signal
  REAL(wp)                 :: sih    ! covariance of L1 ionospheric signal
  REAL(wp)                 :: sem    ! mean square of relative model error
  REAL(wp), DIMENSION(2,2) :: K      ! system matrix
  REAL(wp), DIMENSION(2,2) :: cs     ! inverse signal covariance
  REAL(wp), DIMENSION(2,2) :: cn     ! inverse noise covariance
  REAL(wp), DIMENSION(2,2) :: KC     ! K^T*CN
  REAL(wp), DIMENSION(2,2) :: KCK    ! K^T*CN*K
  REAL(wp), DIMENSION(2,2) :: KK     ! K^T*CN*K + CS
  REAL(wp), DIMENSION(2,2) :: KI     ! (K^T*CN*K + CS)^(-1)
  REAL(wp), DIMENSION(2,2) :: KQI    ! (K^T*CN*K + CS)^(-1)*K^T*CN
  REAL(wp), DIMENSION(2)   :: ba_12  ! L1 and L2 bending angle
  REAL(wp), DIMENSION(2)   :: ba_ni  ! neut and ion bending angle

  REAL(wp) :: ba_diff0, ba_diffS
  REAL(wp) :: impact_lt
  REAL(wp) :: pmin, pmax
  INTEGER  :: iil, iiu, i_str, i_ltr
  INTEGER  :: i, n_obs, nh
  INTEGER  :: whi, wei
  INTEGER  :: w_smooth

!-------------------------------------------------------------------------------
! 2. Useful variables
!-------------------------------------------------------------------------------

  n_obs = SIZE(impact_L1)
  pmax = MAXVAL(impact_L1)
  pmin = MINVAL(impact_L1)
  
  ALLOCATE(ba_LI(n_obs))
  ALLOCATE(err_LI(n_obs))
  ALLOCATE(err_LC(n_obs))
  ALLOCATE(WLC(n_obs))

! IMPLEMENTATION FROM 'INVERT'
  IF (config%dpi < 50.) THEN
    w_smooth = Ceiling(config%f_width*(config%npoints-1) / &
               Abs(config%Pmax-config%Pmin))
  ELSE
! IMPLEMENTATION FROM 'OCC'
    w_smooth = CEILING(config%f_width/config%dpi)
  ENDIF

!-------------------------------------------------------------------------------
! 3. Array allocation
!-------------------------------------------------------------------------------

  nh = CEILING( ( pmax - pmin ) / config%delta_p ) + 1

  ALLOCATE(impact_LH(nh))
  ALLOCATE(ba_LMH(nh))  
  ALLOCATE(ba_is(nh))
  ALLOCATE(ba_ion(nh))
  ALLOCATE(ba_neut(nh))
  ALLOCATE(ba_ifilt(nh))
  ALLOCATE(ba_nfilt(nh))
  ALLOCATE(c_ifilt(nh))
  ALLOCATE(c_nfilt(nh))
  ALLOCATE(wlch(nh))
  ALLOCATE(ba_diff(2,nh))
  ALLOCATE(ba_low(2,nh))

!-------------------------------------------------------------------------------
! 4. Interpolation to homogeneous grid
!-------------------------------------------------------------------------------

  DO i=1,nh
     impact_LH(i) = pmin + (i-1.0_wp)*( pmax - pmin )/(nh-1.0_wp)
  ENDDO

  ALLOCATE(ba_L1H(nh))
  ALLOCATE(ba_L2H(nh)) 

  CALL ropp_pp_interpol(impact_L1, impact_LH, ba_L1, ba_L1H)
  CALL ropp_pp_interpol(impact_L2, impact_LH, ba_L2, ba_L2H)
  CALL ropp_pp_interpol(impact_LM, impact_LH, ba_LM, ba_LMH)

  ba_diff(1,:) = ba_L1H(:) - ba_LMH(:)
  ba_diff(2,:) = ba_L2H(:) - ba_LMH(:)
     
  DEALLOCATE(ba_L1H)
  DEALLOCATE(ba_L2H) 

!-------------------------------------------------------------------------------
! 5. Calculate smoothing windows
!-------------------------------------------------------------------------------
  
  whi = CEILING(w_smooth * (nh-1.0_wp)/( n_obs - 1.0_wp))
  wei = CEILING(config%s_smooth * (nh-1.0_wp)/( pmax - pmin))

!-------------------------------------------------------------------------------
! 6. Region separation
!-------------------------------------------------------------------------------

  ! 6.1 Stratospheric upper limit

  i_str = SUM(MINLOC(impact_LH(:), impact_LH(:)-config%r_curve > config%z_str))
  
  ! 6.2 Initial estimation of upper limit of lower-tropospheric L2 perturbations

  i_ltr = SUM(MINLOC(impact_LH(:), impact_LH(:)-config%r_curve > config%z_ltr))

  ! 6.3 Lower limit of ionospheric signal/noise

  iil = SUM(MINLOC(impact_LH(:), impact_LH(:)-config%r_curve > config%z_ion))
  IF(iil < 1 .OR. iil > nh) iil = 0

  ! 6.4 Dynamic upper limit for ionospheric noise estimate
  
  IF(iil == 0) THEN
    iiu = nh
  ELSE
    ba_diff0 = SUM(ba_diff(2,i_ltr:iil)-ba_diff(1,i_ltr:iil))/(iil-i_ltr+1.0_wp)
    ba_diffS = SQRT(SUM((ba_diff(2,i_ltr:iil) - ba_diff(1,i_ltr:iil) -    &
                    ba_diff0)**2) / (iil-i_ltr+1.0_wp))
    
    ALLOCATE(m_diff(nh))
    m_diff = ABS( ba_diff(2,:) - ba_diff(1,:) - ba_diff0) > 6.0_wp*ba_diffS
    IF (ANY(m_diff(iil:nh) .AND.   &
           impact_LH(iil:nh)-config%r_curve > 75000.0_wp)) THEN
      iiu = SUM(MaxLoc(impact_LH(:), impact_LH(:)-config%r_curve < 70000.0_wp))
    ELSE IF (ANY(m_diff(iil:nh) .AND. &
           impact_LH(iil:nh)-config%r_curve > 70000.0_wp)) THEN
      iiu = SUM(MaxLoc(impact_LH(:), impact_LH(:)-config%r_curve < 65000.0_wp))
    ELSE IF (ANY(m_diff(iil:nh) .AND. &
           impact_LH(iil:nh)-config%r_curve > 65000.0_wp)) THEN
      iiu = SUM(MaxLoc(impact_LH(:), impact_LH(:)-config%r_curve < 60000.0_wp))  
    ELSE
      iiu = SUM(MaxLoc(impact_LH(:), impact_LH(:)-config%r_curve < 80000.0_wp)) 
    ENDIF
    DEALLOCATE(m_diff)
  ENDIF
  
!-------------------------------------------------------------------------------
! 7. Covariance estimation
!-------------------------------------------------------------------------------

  ! 7.1 Smoothing with external scale 

  SELECT CASE (config%filter_method)
  CASE('optest')
    CALL ropp_pp_filter(impact_LH(2) - impact_LH(1), ba_diff(:, i_ltr:nh), wei, & 
                        config%n_smooth, ba_low(:, i_ltr:nh))
  CASE('slpoly')
    CALL ropp_pp_sliding_polynomial(impact_LH(:), ba_diff(:, i_ltr:nh), wei, & 
                                    config%np_smooth, ba_low(:, i_ltr:nh))
  END SELECT

  ! Calculate low-frequency component of ionospheric refraction for L1

  ba_is(i_ltr:nh) = (ba_low(1,i_ltr:nh) - ba_low(2,i_ltr:nh)) *  &
                          f_L2**2/(f_L2**2 - f_L1**2)
  
  ! 7.2 Smoothing with ionospheric scale

  SELECT CASE (config%filter_method)
  CASE('optest')
    CALL ropp_pp_filter(impact_LH(2) - impact_LH(1), ba_diff(:, i_ltr:nh), whi, & 
                        config%n_smooth, ba_low(:, i_ltr:nh))
  CASE('slpoly')
    CALL ropp_pp_sliding_polynomial(impact_LH(:), ba_diff(:, i_ltr:nh), whi, & 
                                    config%np_smooth, ba_low(:, i_ltr:nh))
  END SELECT

  ! 7.3 Estimation of ionospheric noise covariance

  ba_neut(i_ltr:nh) = ba_low(1,i_ltr:nh) -   &
                              ba_low(2,i_ltr:nh) * (f_L2/f_L1)**2       

  IF (iil > 0) THEN
     cin(1) = (SUM(ba_neut(iil:iiu)**2)/(iiu-iil+1))/2.0_wp
     cin(2) = ((f_L1/f_L2)**4)*(SUM(ba_neut(iil:iiu)**2)/(iiu-iil+1))/2.0_wp
  ELSE
     cin(:) = 0.0_wp
  ENDIF

  ! 7.4 Estimation of ionospheric signal covariance
 
  ba_ion(i_ltr:nh) = (ba_low(1,i_ltr:nh) +    &
                            ba_low(2,i_ltr:nh) * (f_L2/f_L1)**2)/2.0_wp -     &
                             ba_is(i_ltr:nh)
  
  IF (iil > 0) THEN
     cis(1) = SUM(ba_ion(iil:iiu)**2)/(iiu-iil+1)
     cis(2) = ((f_L1/f_L2)**4)*(SUM(ba_ion(iil:iiu)**2)/(iiu-iil+1))
     cis(:) = MAX(0.0_wp, cis(:) - cin(:)/2)
  ELSE
     cis(:) = 0.0_wp
  ENDIF

  ! 7.5 Estimation of neutral signal covariance

  IF (config%model_err < 0.0_wp) THEN
    sem = SUM((ba_diff(1, i_ltr:i_str)/ba_LMH(i_ltr:i_str))**2)/(i_str-i_ltr+1)  
  ELSE
    sem = config%model_err**2
  ENDIF

  ! 7.6 Final estimation of lower tropospheric height

  ba_ion(:) = (ba_diff(1, :) - ba_diff(2, :)) *    &
                    f_L2**2/(f_L2**2 - f_L1**2)

  sih = SUM((ba_ion(i_ltr:i_str)-ba_is(i_ltr:i_str))**2)/(i_str-i_ltr+1)

  i = SUM(MAXLOC(impact_LH(1:i_ltr), &
              (ba_ion(1:i_ltr) - ba_is(i_ltr+wei))**2 > 100*sih))
  IF ( i < 1 .OR. i > i_ltr ) THEN
     i = 0
  ENDIF

  i = MAX(1, i)
  impact_lt = impact_LH(i) + 1000.0_wp
  i_ltr = SUM(MINLOC(impact_LH(:), impact_LH(:) >= impact_lt))

!-------------------------------------------------------------------------------
! 8. Statistical optimization
!-------------------------------------------------------------------------------

  ! 8.1 Smoothing with external scale

  SELECT CASE(config%filter_method)
  CASE('optest')
    CALL ropp_pp_filter(impact_LH(2) - impact_LH(1), ba_diff(:, i_ltr:nh), wei, & 
                        config%n_smooth, ba_low(:, i_ltr:nh))
  CASE('slpoly')
    CALL ropp_pp_sliding_polynomial(impact_LH(:), ba_diff(:, i_ltr:nh), wei, & 
                                    config%np_smooth, ba_low(:, i_ltr:nh))
  END SELECT

  ba_is(i_ltr:nh) = (ba_low(1,i_ltr:nh) - ba_low(2,i_ltr:nh)) *  &
                          f_L2**2/(f_L2**2 - f_L1**2)

  ! 8.2 Smoothing with ionospheric scale

  SELECT CASE(config%filter_method)
  CASE('optest')
    CALL ropp_pp_filter(impact_LH(2) - impact_LH(1), ba_diff(:, i_ltr:nh), whi, & 
                        config%n_smooth, ba_low(:, i_ltr:nh))
  CASE('slpoly')
    CALL ropp_pp_sliding_polynomial(impact_LH(:), ba_diff(:, i_ltr:nh), whi, & 
                                    config%np_smooth, ba_low(:, i_ltr:nh))
  END SELECT

  ! 8.3 Border elimination

  i_ltr = i_ltr + wei
  impact_lt = impact_LH(i_ltr)

  ! 8.4 Matrix calculation
  
  K(:,:) = 1.0_wp
  K(2,2) = (f_L1/f_L2)**2

  IF (iil > 0) THEN
     cn(:,:) = Diag(1.0_wp/cin(:))
  ELSE
     cn(:,:) = Diag((/1.0_wp, 1.0_wp/))
  ENDIF
  
  KC(:,:)  = MATMUL(TRANSPOSE(K), CN)
  KCK(:,:) = MATMUL(KC, K)

  ! 8.5 Optimal linear combination of L1 and L2 bending angles

  DO i=i_ltr, nh
     
     IF (iil > 0) THEN
        cs(:,:) = Diag((/1.0_wp/(sem*ba_LMH(i)**2), 1.0_wp/cis(1) /))
     ELSE
        cs(:,:) = 0.0_wp
     ENDIF

     KK(:,:) = KCK(:,:) + CS(:,:)

     CALL ropp_pp_invert_matrix(KK, KI)

     ba_12 = (/ ba_diff(1,i) - ba_is(i), &
                ba_diff(2,i) - ba_is(i)*(f_L1/f_L2)**2 /)

     ba_ni = MATMUL(KI, MATMUL(KC, ba_12))

     KQI = MATMUL(KI, KC)
     wlch(i) = KQI(1,1) + KQI(1,2)

     ba_nfilt(i) = ba_ni(1) + ba_LMH(i)
     ba_ifilt(i) = ba_ni(2) + ba_is(i)

     c_nfilt(i) = KI(1,1)
     c_ifilt(i) = KI(2,2)

  ENDDO

  ! 8.6 Interpolation
    
  DO i=1,n_obs
 
     IF (impact_L1(i) > impact_lt) THEN
        
       CALL ropp_pp_interpol(impact_LH(i_ltr:nh),impact_L1(i),  &
                               ba_ifilt(i_ltr:nh),ba_LI(i))
       CALL ropp_pp_interpol(impact_LH(i_ltr:nh),impact_L1(i),  &
                               ba_nfilt(i_ltr:nh),ba_LC(i))
       CALL ropp_pp_interpol(impact_LH(i_ltr:nh),impact_L1(i),  &
                               c_ifilt(i_ltr:nh),err_LI(i))
       CALL ropp_pp_interpol(impact_LH(i_ltr:nh),impact_L1(i),  &
                               c_nfilt(i_ltr:nh),err_LC(i))
       CALL ropp_pp_interpol(impact_LH(i_ltr:nh),impact_L1(i),  &
                               wlch(i_ltr:nh),WLC(i))

     ELSE

       ba_LI(i) = ba_is(i_ltr)
       ba_LC(i) = ba_L1(i) - ba_LI(i)
       
       err_LI(i) = c_ifilt(i_ltr)
       err_LC(i) = c_nfilt(i_ltr)

       WLC(i) = 1.0_wp

     ENDIF

  ENDDO

  impact_LC(:) = impact_L1(:)

!-------------------------------------------------------------------------------
! 9. Linear combination plus statistical optimization (if required)
!-------------------------------------------------------------------------------

  SELECT CASE(config%so_method)
  CASE ( 'lcso' )  
    cin(1) = cin(1)*2.0_wp*f_L1**4 / (f_L1**2 - f_L2**2)**2
    ba_LC(:) = ba_L1(:) - ba_LI(:)
    ba_LC(:) = ba_LM(:) + (ba_LC(:)-ba_LM(:))*sem*ba_LM(:)**2 /   &
                  (sem*ba_LM(:)**2 + cin(1))
    
  END SELECT
  
!-------------------------------------------------------------------------------
! 10. Output diagnostics
!-------------------------------------------------------------------------------

  diag_out%sq = 100.0_wp * MAXVAL(SQRT(err_LC(:))/ba_LC(:))  ! 'Badness score'
  ALLOCATE(diag_out%ba_ion(n_obs))
  diag_out%ba_ion(:) =  ba_LI(:)         ! Ionospheric bending angle
  ALLOCATE(diag_out%err_ion(n_obs))
  diag_out%err_ion(:) = err_LI(:)        ! Error covariance ionospheric bending
  ALLOCATE(diag_out%err_neut(n_obs))
  diag_out%err_neut(:) = err_LC(:)       ! Error covariance neutral bending
  ALLOCATE(diag_out%wt_data(n_obs))
  diag_out%wt_data(:)  = WLC(:)          ! Ratio data:(data+clim) in profile

!-------------------------------------------------------------------------------
! 11. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(impact_LH)
  DEALLOCATE(ba_LMH)
  DEALLOCATE(ba_is)
  DEALLOCATE(ba_ion)
  DEALLOCATE(ba_neut)
  DEALLOCATE(ba_ifilt)
  DEALLOCATE(ba_nfilt)
  DEALLOCATE(c_ifilt)
  DEALLOCATE(c_nfilt)
  DEALLOCATE(ba_diff)
  DEALLOCATE(ba_low)
  DEALLOCATE(wlch)

  DEALLOCATE(ba_LI)
  DEALLOCATE(err_LI)
  DEALLOCATE(err_LC)
  DEALLOCATE(WLC)

CONTAINS

!-------------------------------------------------------------------------------
! 11. Generation of diagonal matrix
!-------------------------------------------------------------------------------

  FUNCTION diag(A)       
        
    IMPLICIT NONE
    
    REAL(wp), DIMENSION(:),  INTENT(in)  :: A     ! Array of diagonal elements
    REAL(wp), DIMENSION(SIZE(A),SIZE(A)) :: diag  ! Diagonal matrix
    
    INTEGER :: i   ! Diagonal index

    diag(:,:) = 0.0_wp
    
    DO i=1,SIZE(A)
       diag(i,i) = A(i)
    ENDDO
    
  END FUNCTION diag
  
END SUBROUTINE ropp_pp_ionospheric_correction
