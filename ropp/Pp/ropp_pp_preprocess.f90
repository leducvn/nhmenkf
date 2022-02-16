! $Id: ropp_pp_preprocess.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_preprocess(ro_data, config, diag)

!****s* Preprocessing/ropp_pp_preprocess *
!
! NAME
!    ropp_pp_preprocess - Level1a data preprocessing
!
! SYNOPSIS
!    call ropp_pp_preprocess(ro_data, config, diag)
!
! DESCRIPTION
!
! INPUTS
!    type(ROprof)   :: ro_data      ! Radio occultation data strucuture
!    type(PPConfig) :: config       ! Configuration options
!
! OUTPUT
!    type(ROprof)   :: ro_data      ! Corrected radio occultation data
!    type(PPConfig) :: config       ! Configuration options
!    type(PPDiag)   :: diag         ! Diagnostic output
!
! NOTES
!   Requires ROprof data structure type, defined in ropp_io module. This
!   routine therefore requires that the ropp_io module is pre-installed before
!   compilation.
!
! REFERENCES
!   Gorbunov M.E., Lauritsen K.B., Rhodin A., Tomassini M. and Kornblueh L.
!   2006
!   Radio holographic filtering, error estimation, and quality control of
!   radio occultation data
!   Journal of Geophysical Research (111) D10105
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
  USE ropp_io_types, ONLY: ROprof, PCD_rising, PCD_open_loop
! USE ropp_pp_preproc, not_this => ropp_pp_preprocess
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPConfig, PPDiag
  USE ropp_pp_constants
  USE ropp_utils

  IMPLICIT NONE

  TYPE(ROprof),   INTENT(inout) :: ro_data    ! Radio occultation data strucutre
  TYPE(PPconfig), INTENT(inout) :: config     ! Configuration options
  TYPE(PPdiag),   INTENT(inout) :: diag       ! Diagnostic output

  REAL(wp), DIMENSION(:), POINTER :: phase_LM => null() ! Model excess phase (m)
  REAL(wp), DIMENSION(:), POINTER :: impact_LM => null() ! Model impact (m)
  INTEGER,  DIMENSION(:), POINTER :: LCF => null()       ! Lost carrier flag
  INTEGER                         :: w_smooth   ! Smoothing window (points)
  REAL(wp)                        :: Pmax, Pmin ! Max/min impact parameter
  REAL(wp)                        :: ps1, psN   ! Start/end impact parameter
  REAL(wp)                        :: secs_past_hour ! No. of seconds past hour
  INTEGER                         :: n          ! Number of data points
  INTEGER                         :: ocd        ! Occ direction
  CHARACTER(len =   10)           :: nstr
  CHARACTER(len = 256)            :: routine
  REAL(wp), PARAMETER             :: DHS = 5000.0_wp  ! Safety border

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_preprocess')

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  config%obs_ok = .TRUE.

  n = ro_data%Lev1a%Npoints
  ALLOCATE(LCF(n))
  LCF(:) = 0

!-------------------------------------------------------------------------------
! 3. Data cut-off from amplitude (snr_L1p) and missing data flag (LCF)
!-------------------------------------------------------------------------------

  CALL ropp_pp_cutoff_amplitude(ro_data, LCF, config)

!-------------------------------------------------------------------------------
! 4. Mission-specific pre-processing
!-------------------------------------------------------------------------------

  IF (ANY(ro_data%Lev1a%r_leo == ropp_MDFV) .OR.  &
     ANY(ro_data%Lev1a%r_gns == ropp_MDFV)) THEN
    CALL message(msg_warn, 'Invalid coordinate values. Exiting processing.')
    config%obs_ok = .false.
    RETURN
  ENDIF

  SELECT CASE (ro_data%leo_id(1:2))

  ! 4.1 COSMIC

  CASE ('C0','CO')
     CALL message(msg_info, 'COSMIC data preprocessing')
     CALL ropp_pp_preprocess_COSMIC(ro_data, config, LCF)

  ! 4.2 CHAMP and GRACE

  CASE ('CH','GR')
     CALL message(msg_info, 'CHAMP data preprocessing')

     ! 3.2.2 Correct L2 amplitude if not defined

     IF ( ANY(ro_data%Lev1a%snr_L2p == 0.0_wp) ) THEN
        ro_data%Lev1a%snr_L2p(:) = ro_data%Lev1a%snr_L1ca(:)
     ENDIF

  CASE ('ME','MT')
    CALL message(msg_info, 'GRAS data preprocessing')
  
    IF (ASSOCIATED(ro_data%vlist%VlistD1d)) THEN
      IF (ro_data%vlist%VlistD1d%name == 'open_loop_lcf') THEN
        CALL ropp_pp_preprocess_GRASRS(ro_data, config, LCF)
        CALL ropp_pp_cutoff_amplitude(ro_data, LCF, config)
      ENDIF
    ENDIF

    WRITE(nstr, '(i6)') ro_data%lev1a%npoints
    CALL message(msg_diag, 'Merged RS+CL data size: ' // nstr)

  CASE default

    CALL message(msg_warn, 'Occultation LEO id '// TRIM(ro_data%leo_id) // &
        ' not recognised \n. ' // &
        ' No mission-specific data pre-processing conducted.')

  END SELECT

!-------------------------------------------------------------------------------
! 5. Compute model excess phase
!-------------------------------------------------------------------------------

  IF (ro_data%Lev1a%Npoints .lt. 100) THEN
    CALL message(msg_warn, "Error: Too few data")
    config%obs_ok = .FALSE.
    RETURN
  ENDIF

  n = ro_data%Lev1a%Npoints
  ALLOCATE(phase_LM(n))
  ALLOCATE(impact_LM(n))

  CALL ropp_pp_modelphase(ro_data%dtocc%month, ro_data%georef%lat,      &
                          ro_data%georef%lon,  ro_data%Lev1a%dtime,     &
                          ro_data%lev1a%r_leo, ro_data%lev1a%r_gns,     &
                          ro_data%georef%r_coc, ro_data%georef%roc,     &
                          phase_LM, impact_LM, config)

!-------------------------------------------------------------------------------
! 6. Open loop processing (GRAS and COSMIC missions only)
!-------------------------------------------------------------------------------

  SELECT CASE (ro_data%leo_id(1:2))

  CASE ('ME','MT')

    CALL message(msg_info,'GRAS data: openloop preprocessing')

    WHERE(MOD(LCF(:),2) /= 0)
      ro_data%Lev1a%phase_L2(:) = phase_LM(:)
    ENDWHERE

    WHERE(ro_data%Lev1a%phase_L2(:) < ropp_MDTV)
      ro_data%Lev1a%phase_L2(:) = -1.0_wp
      ro_data%Lev1a%phase_L2(:) = phase_LM(:)
    ENDWHERE

    secs_past_hour = 0.0_wp
    IF(ro_data%DTocc%Minute .LE.  59) &
      secs_past_hour = secs_past_hour + ro_data%DTocc%Minute*60.0_wp
    IF(ro_data%DTocc%Second .LE.  59) &
      secs_past_hour = secs_past_hour + ro_data%DTocc%Second*1.0_wp
    IF(ro_data%DTocc%Msec   .LT. 999) &  ! Sometimes msec not set, so defaults to 999
      secs_past_hour = secs_past_hour + ro_data%DTocc%Msec*1.0e-3_wp

    CALL ropp_pp_openloop(ro_data%Lev1a%dtime+secs_past_hour, &
                          ro_data%Lev1a%phase_L1, ro_data%Lev1a%phase_L2,     &
                          phase_LM, ro_data%lev1a%r_leo, ro_data%lev1a%r_gns, &
                          ro_data%georef%r_coc,   LCF)

    ! Set PCD flag
    ro_data%PCD = IBSET(ro_data%PCD, PCD_open_loop)

  CASE ('C0','CO')
    CALL message(msg_info,'COSMIC data: openloop preprocessing')

    CALL ropp_pp_openloop(ro_data%Lev1a%dtime,    ro_data%Lev1a%phase_L1,  &
                          ro_data%Lev1a%phase_L2, phase_LM,                &
                          ro_data%lev1a%r_leo,    ro_data%lev1a%r_gns,     &
                          ro_data%georef%r_coc,   LCF)

     ! Set PCD flag
     ro_data%PCD = IBSET(ro_data%PCD, PCD_open_loop)
  END SELECT

!-------------------------------------------------------------------------------
! 7. Data cutoff from bending angle and downsampling
!-------------------------------------------------------------------------------

    CALL ropp_pp_cutoff(ro_data, config, phase_LM, impact_LM, LCF)

!-------------------------------------------------------------------------------
! 8. Calculate spectra
!-------------------------------------------------------------------------------

  IF (config%opt_spectra) THEN
   CALL ropp_pp_spectra(ro_data%Lev1a%dtime, ro_data%Lev1a%phase_L1,     &
                        ro_data%Lev1a%phase_L2, phase_LM, impact_LM, config, &
                        OutRO=.true.)
  ENDIF

!-------------------------------------------------------------------------------
! 9. Degraded L2 amplitude and excess phase data correction
!-------------------------------------------------------------------------------

  Pmax = MAXVAL(impact_LM)
  Pmin = MAX(MINVAL(impact_LM), ro_data%georef%roc+2000.0_wp)
  w_smooth = CEILING(config%fw_go_smooth * (ro_data%Lev1a%Npoints - 1) /     &
                ABS(Pmax - Pmin))

  IF (config%opt_DL2) THEN

    CALL ropp_pp_amplitude_go(ro_data%Lev1a%dtime, ro_data%lev1a%r_leo,      &
                              ro_data%lev1a%r_gns, ro_data%georef%r_coc,     &
                              ro_data%georef%roc,  impact_LM,                &
                              ro_data%Lev1a%snr_L1ca, w_smooth,              &
                              ro_data%Lev1a%snr_L2p)

    IF (ALL(ro_data%Lev1a%snr_L2p == ropp_MDFV)) THEN
      CALL message(msg_warn, "Error: data unusable [L2 amplitude]")
      config%obs_ok = .FALSE.
      RETURN
    ENDIF

    CALL ropp_pp_correct_L2(ro_data%Lev1a%dtime,    ro_data%lev1a%r_leo,     &
                            ro_data%lev1a%r_gns,    ro_data%georef%r_coc,    &
                            ro_data%georef%roc,                              &
                            impact_LM,              phase_LM,                &
                            ro_data%Lev1a%phase_L1, ro_data%Lev1a%phase_L2,  &
                            ro_data%Lev1a%snr_L1ca, ro_data%Lev1a%snr_L2p,   &
                            LCF, config%hmax_wo-DHS, diag%L2_badness)

    IF (diag%L2_badness > Huge(diag%L2_badness)) THEN
      CALL message(msg_warn, "Error: data unusable [L2 badness]")
      config%obs_ok = .FALSE.
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
! 10. Set rising/setting flag
!-------------------------------------------------------------------------------

  n = ro_data%lev1a%npoints
  ps1 = impact_parameter(ro_data%lev1a%r_leo(1,:) - ro_data%georef%r_coc(:), &
                         ro_data%lev1a%r_gns(1,:) - ro_data%georef%r_coc)
  psN = impact_parameter(ro_data%lev1a%r_leo(n,:) - ro_data%georef%r_coc(:), &
                         ro_data%lev1a%r_gns(n,:) - ro_data%georef%r_coc)

  ocd = NINT(SIGN(1.0_wp, psN - ps1))
  IF (ocd == 1) ro_data%PCD = IBSET(ro_data%PCD, PCD_rising)

!-------------------------------------------------------------------------------
! 11. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(LCF)
  DEALLOCATE(impact_LM)
  DEALLOCATE(phase_LM)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_pp_preprocess
