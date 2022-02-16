! $Id: ropp_qc_genqc.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* QC/ropp_qc_genqc *
!
! NAME
!    ropp_qc_genqc - Generic quality control checks.
!
! SYNOPSIS
!    call ropp_qc_genqc(obs, bg, config, diag)
! 
! DESCRIPTION
!    This subroutine performs generic quality control checks, e.g. for the 
!    proper co-location of observation and background data in both space and 
!    time, for fullfilling monotonicity constraints, and also basic unit checks.
!
! INPUTS
!    obs         Observation vector
!    bg          Background vector
!    config      Configuration options
! 
! OUTPUT
!    diag        Diagnostics structure with updated status
!
! NOTES
!
!
! EXAMPLE
!
!
! SEE ALSO
!
!
! REFERENCES
!
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
! 1. Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_genqc_bangle(obs, bg, config, diag)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_genqc_bangle

  IMPLICIT NONE

  TYPE(Obs1DBangle)                     :: obs
  TYPE(State1DFM)                       :: bg
  TYPE(VarConfig)                       :: config
  TYPE(VarDiag)                         :: diag

  INTEGER                               :: cnt = 0
  INTEGER, DIMENSION(:), POINTER        :: idx => null()

  REAL(wp)                              :: z_sfc_guess
  REAL(wp)                              :: p_sfc_guess
  REAL(wp)                              :: distance
  REAL(wp)                              :: dtime

  CHARACTER(len =  256)                 :: routine
  CHARACTER(len =   15)                 :: cvalue1, cvalue2

! 1.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_genqc')

! 1.3 Background data checks
! --------------------------

! 1.3.1 Temperature

  IF (ANY(bg%temp < config%genqc_min_temperature) .OR. &
      ANY(bg%temp > config%genqc_max_temperature)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_temperature
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_temperature
     CALL message(msg_error, &
          "Background temperature values outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] K.")
     diag%ok = .FALSE.
  ENDIF

! 1.3.2 Humidity

  IF (ANY(bg%shum < config%genqc_min_spec_humidity/1000._wp) .OR. &  ! kg/kg
      ANY(bg%shum > config%genqc_max_spec_humidity/1000._wp)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_spec_humidity/1000._wp
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_spec_humidity/1000._wp
     CALL message(msg_error, &
          "Background humidity values outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] kg/kg.")
     diag%ok = .FALSE.
  ENDIF

! 1.3.3 Surface pressure    
  z_sfc_guess = ropp_MDFV
  p_sfc_guess = ropp_MDFV

  IF (bg%geop_sfc > ropp_MDTV) THEN
    z_sfc_guess = bg%geop_sfc        ! surface geopotential height
    p_sfc_guess = 1013.0_wp * EXP(- z_sfc_guess / 7000.0_wp) * 100.0_wp ! Pa
  ENDIF
  IF (ABS(bg%pres(1) - p_sfc_guess) > 100.e2_wp) THEN 
     CALL message(msg_error, &
          "Surface pressure and elevation / geopotential are inconsistent.")
     diag%ok = .FALSE.
  ENDIF

! 1.4 Observation data checks
! ---------------------------

! 1.4.1 Valid bending angle values

  idx => WHERE(obs%bangle > ropp_MDTV .AND. obs%weights > 0.0_wp, cnt)

! 1.4.1 Impact parameter

  IF (ANY(obs%impact(idx) < config%genqc_min_impact) .OR. &
      ANY(obs%impact(idx) > config%genqc_max_impact)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_impact
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_impact
     CALL message(msg_error, &
          "Impact parameter values outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] m.")
     diag%ok = .FALSE.
  ENDIF

! 1.4.2 Bending angles

  IF (ANY(obs%bangle(idx) < config%genqc_min_bangle) .OR. &
      ANY(obs%bangle(idx) > config%genqc_max_bangle)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_bangle
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_bangle
     CALL message(msg_error, &
          "Bending angle values outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] rad.")
     diag%ok = .FALSE.
  ENDIF

! 1.4.3 Radius of curvature

  IF (obs%r_curve < config%genqc_min_impact .OR. &
      obs%r_curve > config%genqc_max_impact) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_impact
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_impact
     CALL message(msg_error, &
          "Radius of curvature value outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] m.")
     diag%ok = .FALSE.
  ENDIF

! 1.4.4 Minimum required observation height

  IF (MINVAL(obs%impact(idx)-obs%r_curve) > config%genqc_min_obheight) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_obheight
     CALL message(msg_error, &
       "Observation profile does not reach below minimum required height: " &
       //TRIM(ADJUSTL(cvalue1))//" m.")
     diag%ok = .FALSE.
  ENDIF
  
! 1.5 Co-location
! ---------------

  IF (config%genqc_colocation_apply) THEN

     ! 1.5.1 Co-location in space (km)

     distance = great_circle_distance(obs%lon,obs%lat,bg%lon,bg%lat)/1000.0_wp
     
     IF (distance > config%genqc_max_distance) THEN
        WRITE(cvalue2, '(1pe15.5)') config%genqc_max_distance
        CALL message(msg_error, &
             "Great circle distance between observations and " // &
             "background profile exceeds upper limit: " &
             //TRIM(ADJUSTL(cvalue2))//" km.")
        diag%ok = .FALSE.
     ENDIF

     ! 1.5.2 Co-location in time

     dtime = ABS(bg%time - obs%time)   ! s
     IF(dtime > config%genqc_max_time_sep)THEN
        WRITE(cvalue2, '(1pe15.5)') config%genqc_max_time_sep
        CALL message(msg_error, &
             "Temporal separation between observations and " // &
             "background profile exceeds upper limit: " &
             //TRIM(ADJUSTL(cvalue2))//" s.")
        diag%ok = .FALSE.
     ENDIF

  ENDIF

! 1.6 Clean up
! ------------

  DEALLOCATE(idx)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_genqc_bangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_qc_genqc_refrac(obs, bg, config, diag)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE matrix
  USE ropp_utils
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_qc, not_this => ropp_qc_genqc_refrac

  IMPLICIT NONE

  TYPE(Obs1DRefrac)                     :: obs
  TYPE(State1DFM)                       :: bg
  TYPE(VarConfig)                       :: config
  TYPE(VarDiag)                         :: diag

  INTEGER                               :: cnt = 0
  INTEGER, DIMENSION(:), POINTER        :: idx => null()

  REAL(wp)                              :: z_sfc_guess
  REAL(wp)                              :: p_sfc_guess
  REAL(wp)                              :: distance
  REAL(wp)                              :: dtime

  CHARACTER(len =  256)                 :: routine
  CHARACTER(len =   15)                 :: cvalue1, cvalue2

! 2.2 Message handling
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_qc_genqc')

! 2.3 Background data checks
! --------------------------

! 2.3.1 Temperature
  IF (ANY(bg%temp < config%genqc_min_temperature) .OR. &
      ANY(bg%temp > config%genqc_max_temperature)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_temperature
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_temperature
     CALL message(msg_error, &
          "Background temperature values outside limits: [" &
          //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] K.")
     diag%ok = .FALSE.
  ENDIF

! 2.3.2 Humidity
  IF (ANY(bg%shum < config%genqc_min_spec_humidity/1000._wp) .OR. &  ! kg/kg
      ANY(bg%shum > config%genqc_max_spec_humidity/1000._wp)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_spec_humidity/1000._wp
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_spec_humidity/1000._wp
     CALL message(msg_error, &
           "Background humidity values outside limits: [" &
           //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] kg/kg.")
     diag%ok = .FALSE.
  ENDIF

  z_sfc_guess = ropp_MDFV
  p_sfc_guess = ropp_MDFV

! 2.3.3 Surface pressure
  IF (bg%geop_sfc > ropp_MDTV) THEN
    z_sfc_guess = bg%geop_sfc        ! surface geopotential height
    p_sfc_guess = 1013.0_wp * EXP(- z_sfc_guess / 7000.0_wp) * 100.0_wp ! Pa
  ENDIF
  IF (ABS(bg%pres(1) - p_sfc_guess) > 100.e2_wp) THEN
     CALL message(msg_error, &
          "Surface pressure and elevation / geopotential are inconsistent.")
     diag%ok = .FALSE.
  ENDIF

! 2.4 Observation data checks
! ---------------------------

! 2.4.1 Valid refractivity values

  idx => WHERE(obs%refrac > ropp_MDTV .AND. obs%weights > 0.0_wp, cnt)

! 2.4.2 Geopotential heights

  IF (ANY(obs%geop(idx) < config%genqc_min_geop_refrac) .OR. &
      ANY(obs%geop(idx) > config%genqc_max_geop_refrac)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_geop_refrac
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_geop_refrac
     CALL message(msg_error, &
          "Geopotential height values for refractivity outside limits: [" &
           //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] m.")
     diag%ok = .FALSE.
  ENDIF

  IF (ANY((obs%geop(idx(2:cnt)) - obs%geop(idx(1:cnt-1))) < 0.0_wp)) THEN
     CALL message(msg_error, &
      "Geopotential height values for refractivity must be strictly monotonic.")
     diag%ok = .FALSE.
  ENDIF

! 2.4.3 Refractivity

  IF (ANY(obs%refrac(idx) < config%genqc_min_refractivity) .OR. &
      ANY(obs%refrac(idx) > config%genqc_max_refractivity)) THEN
     WRITE(cvalue1, '(1pe15.5)') config%genqc_min_refractivity
     WRITE(cvalue2, '(1pe15.5)') config%genqc_max_refractivity
     CALL message(msg_error, &
          "Refractivity values outside limits: [" &
           //TRIM(ADJUSTL(cvalue1))//", "//TRIM(ADJUSTL(cvalue2))//"] N-units.")
     diag%ok = .FALSE.
  ENDIF

! 2.4.4 Minimum required observation height

  IF (MINVAL(obs%geop(idx)) > config%genqc_min_obheight) THEN
    WRITE(cvalue1, '(1pe15.5)') config%genqc_min_obheight
    CALL message(msg_error, &
       "Observation profile does not reach below minimum required height: " &
       //TRIM(ADJUSTL(cvalue1))//" m.")
    diag%ok = .FALSE.
  ENDIF

! 2.5 Co-location
! ---------------

  IF (config%genqc_colocation_apply) THEN

     ! 2.5.1 Co-location in space (km)

     distance = great_circle_distance(obs%lon,obs%lat,bg%lon,bg%lat)/1000.0_wp

     IF (distance > config%genqc_max_distance) THEN
        WRITE(cvalue2, '(1pe15.5)') config%genqc_max_distance
        CALL message(msg_error, &
             "Great circle distance between observations and " // &
             "background profile exceeds upper limit: " &
             //TRIM(ADJUSTL(cvalue2))//" km.")
        diag%ok = .FALSE.
     ENDIF

     ! 2.5.2 Co-location in time
     
     dtime = ABS(bg%time - obs%time)   ! s
     IF(dtime > config%genqc_max_time_sep)THEN
        WRITE(cvalue2, '(1pe15.5)') config%genqc_max_time_sep
        CALL message(msg_error, &
             "Temporal separation between observations and " // &
             "background profile exceeds upper limit: " &
             //TRIM(ADJUSTL(cvalue2))//" s.")
        diag%ok = .FALSE.
     ENDIF
     
  ENDIF

! 2.6 Clean up
! ------------

  DEALLOCATE(idx)

  CALL message_set_routine(routine)

END SUBROUTINE ropp_qc_genqc_refrac
