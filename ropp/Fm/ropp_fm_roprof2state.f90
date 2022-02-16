! $Id: ropp_fm_roprof2state.f90 4452 2015-01-29 14:42:02Z idculv $

SUBROUTINE ropp_fm_roprof2state1d(ro_data, x)

!****s* Copying2/ropp_fm_roprof2state *
!
! NAME
!    ropp_fm_roprof2state - Copy elements of an ROprof structure to a 
!                           state vector.
!
! SYNOPSIS
!    type(ROprof)                   :: ro_data
!    type(<some state vector type>) :: x
!       ...
!    call ropp_fm_roprof2state(ro_data, x)
! 
! DESCRIPTION
!    This subroutine copies Level 2b, c and d (if applicable) data from a
!    radio occultation profile data structure into a state vector. Data
!    is also checked for consistency.
!
! INPUTS
!   ro_data  Radio occultation profile data.
!
! OUTPUT
!   x        State vector structure.
!
! NOTES
!   Data is copied into the state vector structure without unit conversion;
!   thus, the units in the ROprof data structure for Level 2b,c, and (if
!   applicable) d must be set to the units used by the state vector 
!   variables. This can be accomplished with the ropp_fm_set_units()
!   subroutine.
!   The definition of the state vector and associated covariance matrix is 
!   model-dependent, set by level_type specified in input data.
!
! SEE ALSO
!   State1dFM
!   ropp_fm_state2roprof
!   ropp_fm_set_units
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
  USE ropp_io 
  USE ropp_io_types, ONLY: ROprof
  USE ropp_fm
  USE ropp_fm_types, ONLY: State1dFM

  IMPLICIT NONE

  TYPE(ROprof),     INTENT(in)       :: ro_data   ! Input ROprof data
  TYPE(State1dFM),  INTENT(inout)    :: x         ! Output state vector

  INTEGER                            :: i, j, n
  INTEGER, DIMENSION(8)              :: DT8
  CHARACTER(len = 256)               :: level_type
  CHARACTER(len = 256)               :: routine
  CHARACTER(len = 10)                :: err_val

!-------------------------------------------------------------------------------
! 2. Error handling
!-------------------------------------------------------------------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2state')

  x%state_ok = .TRUE.

!-------------------------------------------------------------------------------
! 3. Check and copy geolocation and time
!-------------------------------------------------------------------------------

  IF (isinrange(ro_data%georef%lon, ro_data%georef%range%lon)) THEN
    x%lon = ro_data%georef%lon
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lon
    CALL message(msg_warn, &
          "Longitude data for background data out of range or missing. " // &
          "(longitdude value = " // TRIM(err_val) // ")")
  ENDIF
  
  IF (isinrange(ro_data%georef%lat, ro_data%georef%range%lat)) THEN
    x%lat = ro_data%georef%lat
  ELSE
    WRITE(err_val, '(e8.1)') ro_data%georef%lat
    CALL message(msg_warn, &
          "Latitude data for background data out of range or missing. " // &
          "(latitude value = " // TRIM(err_val) // ")")
  ENDIF

  IF ( isinrange(ro_data%bg%year,   ro_data%bg%range%year)  .AND. &
       isinrange(ro_data%bg%month,  ro_data%bg%range%month) .AND. &
       isinrange(ro_data%bg%day,    ro_data%bg%range%day)   .AND. &
       isinrange(ro_data%bg%hour,   ro_data%bg%range%hour)  .AND. &
       isinrange(ro_data%bg%minute, ro_data%bg%range%minute) ) THEN
    DT8 = (/ro_data%bg%year, ro_data%bg%month,  ro_data%bg%day, 0, &
            ro_data%bg%hour, ro_data%bg%minute, 0, 0/)
    CALL TimeSince ( DT8, x%time, 1, Base="JS2000" )
! Fix to sort out DMI variable labelling issue - VT time vs analysis time (cf. 
! ROPP ticket #253).  Also increased max. temporal separation to 3 hours.
! This should no longer (ROPP8.0, 2014) be needed, as DMI now set bg%hour etc to the validity time.
!    IF ( isinrange(ro_data%bg%fcperiod,   ro_data%bg%range%fcperiod)  .AND. &
!         ro_data%processing_centre == 'DMI Copenhagen' .AND. &
!         ro_data%bg%source == 'ECMWF' ) THEN
!      x%time = x%time + 3600.0 * ro_data%bg%fcperiod
!    ENDIF

  ELSE
    CALL message(msg_warn, &
         "Time data for background is out of range or missing.")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Set status flag state%ok to FALSE")
    CALL message(msg_noin, '')
    x%state_ok = .FALSE.
  ENDIF
  
!-------------------------------------------------------------------------------
! 4. Check that profiles are increasing in height - 1st element towards surface
!-------------------------------------------------------------------------------

  CALL ropp_io_ascend(ro_data)

!-------------------------------------------------------------------------------
! 5. Check and copy meteorological parameters
!-------------------------------------------------------------------------------

! 5.1 Check and set level numbers

  IF (ro_data%Lev2b%Npoints == 0) THEN
    CALL message(msg_warn, &
         "RO data has no level 2b (free atmospheric) parameters.")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Set status flag state%ok to FALSE")
    CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
  ENDIF

  x%n_lev  = ro_data%Lev2b%Npoints

! 5.2 Check positive humidity and pressure data
  
  IF (x%use_logq) THEN
    IF (ANY(ro_data%Lev2b%shum <= 0.0_wp)) THEN
      WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%shum)
       CALL message(msg_warn,"One or more humidity values are negative. ")
       CALL message(msg_info,"(Minimum humidity value = "//TRIM(err_val)// ")")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Set status flag state%ok to FALSE")
       CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF
  ENDIF
 
  IF (ANY(ro_data%Lev2b%press <= 0.0_wp)) THEN
    WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%press)
    CALL message(msg_warn,"One or more pressure values are negative. ")
    CALL message(msg_info,"(Minimum pressure value = " // TRIM(err_val) // ")" )
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Set status flag state%ok to FALSE")
    CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
  ENDIF

! 5.3 Temperature data

  ALLOCATE(x%temp(x%n_lev))
  x%temp = ro_data%Lev2b%temp

! 5.4 Humidity data

  ALLOCATE(x%shum(x%n_lev))
  x%shum = ro_data%Lev2b%shum
  WHERE (x%shum <= 0.0_wp)
    x%shum = 1.0e-9_wp
  END WHERE

! 5.5 Pressure data

  ALLOCATE(x%pres(x%n_lev))
  x%pres = ro_data%lev2b%press

! 5.6 Geopotential height data

  ALLOCATE(x%geop(x%n_lev))
  x%geop = ro_data%lev2b%geop
  x%geop_sfc = ro_data%Lev2c%geop_sfc
 
!-------------------------------------------------------------------------------
! 6. Define state1dFM structure: ECMWF-type hybrid levels
!-------------------------------------------------------------------------------
 
  level_type = ro_data%Lev2d%level_type
  IF ( INDEX(level_type,'UNKNOWN') > 0) level_type = ro_data%bg%source
  CALL To_Upper(level_type)

  IF (INDEX(level_type,'HYBRID') > 0 .OR. &
      INDEX(level_type,'ECMWF') > 0) THEN

! 6.1 Hybrid level pressure coefficients

    ALLOCATE(x%ak(ro_data%Lev2d%Npoints))
    ALLOCATE(x%bk(ro_data%Lev2d%Npoints))
    x%ak = ro_data%Lev2d%level_coeff_a
    x%bk = ro_data%Lev2d%level_coeff_b
    IF (x%ak(ro_data%Lev2d%Npoints) == 0.0_wp)    &
         x%ak(ro_data%Lev2d%Npoints) = 1.0e-32_wp

    IF (ANY(x%ak < 0.0_wp) .OR. ALL(x%ak == 0.0_wp)) THEN
       WRITE(err_val, '(e10.2)') MINVAL(x%ak)
       CALL message(msg_warn, "Problem with level coefficients (ak) " //     &
                           "- either all zero, or one or more negative. ")
       CALL message(msg_info,"(Minimum ak value = " // TRIM(err_val) // ")")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Set status flag state%ok to FALSE")
       CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF

    IF (ANY(x%bk < 0.0_wp) .OR. ALL(x%bk == 0.0_wp)) THEN
      WRITE(err_val, '(e10.2)') MINVAL(x%bk)
       CALL message(msg_warn, "Problem with level coefficients (bk) " //     &
                           "- either all zero, or one or more negative. ")
       CALL message(msg_info, "(Minimum bk value = " // TRIM(err_val) // ")")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Set status flag state%ok to FALSE")
       CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF

    IF (ro_data%Lev2c%press_sfc <= 0.0_wp) THEN
      WRITE(err_val, '(e8.1)') ro_data%Lev2c%press_sfc
      CALL message(msg_warn,"Surface pressure value is negative. ")
      CALL message(msg_info,"(Surface pressure value = " // TRIM(err_val) // ")" )
      CALL message(msg_info, "Check input data and valid_range attributes")
      CALL message(msg_info, "Set status flag state%ok to FALSE")
      CALL message(msg_noin, '')
      x%state_ok = .FALSE.
      RETURN
    ENDIF

! 6.2 Define state vector - ECMWF

    IF (ro_data%Lev2c%direct_ion) x%direct_ion = .TRUE.

    n = 2*x%n_lev+1 ! Number of elements in the state vector

    IF (x%direct_ion) n = n + 3                            ! Append {Ne_max, H_peak, H_width} to state vector

    ALLOCATE(x%state(n))

    x%state(1:x%n_lev) = x%temp                            ! temperature on full levels

    IF(x%use_logq)THEN                                     ! log(shum [g/kg]) on full levels
       x%state(x%n_lev+1:2*x%n_lev) = LOG(x%shum*1000.0_wp)
    ELSE
       x%state(x%n_lev+1:2*x%n_lev) = x%shum               ! shum on full levels
    ENDIF

    IF (x%use_logp) THEN                                   ! log(sfc pressure [hPa])
       x%state(2*x%n_lev+1) = LOG(ro_data%Lev2c%press_sfc/100.0_wp)
    ELSE
       x%state(2*x%n_lev+1) = ro_data%Lev2c%press_sfc      ! surface pressure
    ENDIF

! 6.3 Compute pressure and geopotential height on full (humidity) levels
    
    CALL ropp_fm_state2state_ecmwf(x)

! 6.4 Include model ionospheric parameters

    x%Ne_max  = ro_data%Lev2c%Ne_max
    x%H_peak  = ro_data%Lev2c%H_peak
    x%H_width = ro_data%Lev2c%H_width

    IF (x%direct_ion) THEN

      x%state(2*x%n_lev+2) = x%Ne_max
      x%state(2*x%n_lev+3) = x%H_peak
      x%state(2*x%n_lev+4) = x%H_width

      x%n_chap = 1 ! for the moment.

    ENDIF

!-------------------------------------------------------------------------------
! 7. Define error covariances: ECMWF-type hybrid levels
!-------------------------------------------------------------------------------

    x%cov_ok = .TRUE.

! 7.1 Allocate memory

    IF (ASSOCIATED(x%cov%d)) DEALLOCATE(x%cov%d)
    CALL callocate(x%cov%d, n*(n+1)/2)

! 7.2 Temperature sigmas

    IF (ALL(ro_data%Lev2b%temp_sigma > 0.0_wp)) THEN

      DO i = 1, x%n_lev
        j = i
        x%cov%d((j*(j+1))/2) = ro_data%Lev2b%temp_sigma(i)**2
      END DO

    ELSE

      x%cov_ok = .FALSE.

    END IF

! 7.3 Humidity sigmas

    IF (ALL(ro_data%Lev2b%shum_sigma > 0.0_wp)) THEN

      IF (x%use_logq) THEN
        DO i = 1, x%n_lev
          j = x%n_lev + i
          x%cov%d((j*(j+1))/2) = ( ro_data%Lev2b%shum_sigma(i) /   &
                                   ro_data%Lev2b%shum(i) )**2
        END DO
      ELSE
        DO i = 1, x%n_lev
          j = x%n_lev + i
          x%cov%d((j*(j+1))/2) = ro_data%Lev2b%shum_sigma(i)**2
        END DO
      END IF

    ELSE

      x%cov_ok = .FALSE.

    END IF

! 7.4 Surface pressure sigmas

    IF (ro_data%Lev2c%press_sfc_sigma > 0.0_wp) THEN

      j = 2*x%n_lev + 1
      IF (x%use_logp) THEN
        x%cov%d((j*(j+1))/2) = ( ro_data%Lev2c%press_sfc_sigma /    &
                                 ro_data%Lev2c%press_sfc )**2
      ELSE
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%press_sfc_sigma**2
      END IF

    ELSE

      x%cov_ok = .FALSE.

    END IF

! 7.5 Ionospheric sigmas

    IF (x%direct_ion) THEN

      j = 2*x%n_lev + 2
      IF (ro_data%Lev2c%Ne_max_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%Ne_max_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

      j = 2*x%n_lev + 3
      IF (ro_data%Lev2c%H_peak_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%H_peak_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

      j = 2*x%n_lev + 4
      IF (ro_data%Lev2c%H_width_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%H_width_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

    END IF

  END IF

!-------------------------------------------------------------------------------
! 8. Define state1dFM structure: MetOffice-type geopotential height levels 
!-------------------------------------------------------------------------------

  IF ( INDEX(level_type,'METOFFICE') > 0) THEN
     
!     To simulate processing for MetOffice type model data it is assumed that 
!     the input ROprof data contain:
!          ro_data%geop_sfc - Z_A(1): geopotential height on surface A-level
!          ro_data%psfc - p_A(1) : pressure on surface A-level
!          ro_data%geop - Z_B: geopotential height on B-levels (humidity level)
!          ro_data%press - p_A: pressure on A-levels above surface
!          ro_data%shum - q_B: specific humidity on B-levels

! 8.1 Define state vector - METO 
    
    IF (ro_data%Lev2c%press_sfc <= 0.0_wp) THEN
      WRITE(err_val, '(e8.1)') ro_data%Lev2c%press_sfc
      CALL message(msg_warn,"Surface pressure value is negative. ")
      CALL message(msg_info,"(Surface pressure value = " // TRIM(err_val) // ")" )
      CALL message(msg_info, "Check input data and valid_range attributes")
      CALL message(msg_info, "Set status flag state%ok to FALSE")
      CALL message(msg_noin, '')
      x%state_ok = .FALSE.
      RETURN
    ENDIF

    IF (ro_data%Lev2c%direct_ion) x%direct_ion = .TRUE.
    
    n = 2*x%n_lev + 1                               ! Number of elements in the state vector

    IF (x%direct_ion) n = n + 3                     ! Append {Ne_max, H_peak, H_width} to state vector

    ALLOCATE(x%state(n))

    IF (x%use_logp) THEN                            ! log(pres [hPa]) on 'A'-levels
       x%state(1) = LOG(ro_data%Lev2c%press_sfc/100.0_wp)   
       x%state(2:x%n_lev+1) = LOG(ro_data%Lev2b%press/100.0_wp) 
    ELSE
       x%state(1) = ro_data%Lev2c%press_sfc         ! pressure on 'A'-level(1)
       x%state(2:x%n_lev+1) = ro_data%Lev2b%press   ! pressure on 'A'-levels
    ENDIF

    IF (x%use_logq) THEN                            ! log(shum [g/kg]) on 'B'-levels
       x%state(x%n_lev+2:2*x%n_lev+1) = LOG(x%shum*1000.0_wp) 
    ELSE
       x%state(x%n_lev+2:2*x%n_lev+1) = x%shum      ! shum on 'B'-levels
    ENDIF

    
! 8.2 Compute pressure and temperature on full (humidity) levels
    
    CALL ropp_fm_state2state_meto(x)

! 8.3 Include model ionospheric parameters

    x%Ne_max  = ro_data%Lev2c%Ne_max
    x%H_peak  = ro_data%Lev2c%H_peak
    x%H_width = ro_data%Lev2c%H_width

    IF (x%direct_ion) THEN

      x%state(2*x%n_lev+2) = x%Ne_max
      x%state(2*x%n_lev+3) = x%H_peak
      x%state(2*x%n_lev+4) = x%H_width

      x%n_chap = 1 ! for the moment.

    ENDIF

!-------------------------------------------------------------------------------
! 9. Define error covariances: MetOffice-type geopotential height levels 
!-------------------------------------------------------------------------------

    x%cov_ok = .TRUE.

! 9.1 Allocate memory

    IF (ASSOCIATED(x%cov%d)) DEALLOCATE(x%cov%d)
    CALL callocate(x%cov%d, n*(n+1)/2)

! 9.2 Pressure surface sigma

    IF (ro_data%Lev2c%press_sfc_sigma > 0.0_wp) THEN

      j = 1

      IF (x%use_logp) THEN
        x%cov%d((j*(j+1))/2) = (ro_data%Lev2c%press_sfc_sigma/ro_data%Lev2c%press_sfc)**2
      ELSE
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%press_sfc_sigma**2
      END IF

    ELSE

      x%cov_ok = .FALSE.

    ENDIF

! 9.3 Pressure sigmas

    IF (ALL(ro_data%Lev2b%press_sigma > 0.0_wp)) THEN

      IF (x%use_logp) THEN
        DO i = 1, x%n_lev
          j = 1 + i
          x%cov%d((j*(j+1))/2) = ( ro_data%Lev2b%press_sigma(i) /    &
                                   ro_data%Lev2b%press(i) )**2
        END DO
      ELSE
        DO i = 1, x%n_lev
          j = i + 1
          x%cov%d((j*(j+1))/2) = ro_data%Lev2b%press_sigma(i)**2
        END DO
      ENDIF

    ELSE

      x%cov_ok = .FALSE.

    ENDIF
  
! 9.4 Humidity sigmas

    IF (ALL(ro_data%Lev2b%shum_sigma > 0.0_wp)) THEN

      IF (x%use_logp) THEN
        DO i = 1, x%n_lev
          j = x%n_lev + 1 + i
          x%cov%d((j*(j+1))/2) = ( ro_data%Lev2b%shum_sigma(i) /    &
                                   ro_data%Lev2b%shum(i) )**2
        END DO
      ELSE
        DO i = 1, x%n_lev
          j = x%n_lev + 1 + i
          x%cov%d((j*(j+1))/2) = ro_data%Lev2b%shum_sigma(i)**2
        END DO
      ENDIF

    ELSE

      x%cov_ok = .FALSE.

    ENDIF

! 9.5 Ionospheric sigmas

    IF (x%direct_ion) THEN

      j = 2*x%n_lev + 2
      IF (ro_data%Lev2c%Ne_max_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%Ne_max_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

      j = 2*x%n_lev + 3
      IF (ro_data%Lev2c%H_peak_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%H_peak_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

      j = 2*x%n_lev + 4
      IF (ro_data%Lev2c%H_width_sigma > 0.0_wp) THEN
        x%cov%d((j*(j+1))/2) = ro_data%Lev2c%H_width_sigma**2
      ELSE
        x%cov_ok = .FALSE.
      END IF

    END IF

  ENDIF

!-------------------------------------------------------------------------------
! 10. Rest of the covariance marix
!-------------------------------------------------------------------------------

  IF (ASSOCIATED(x%cov%e)) DEALLOCATE(x%cov%e)
  IF (ASSOCIATED(x%cov%f)) DEALLOCATE(x%cov%f)
  IF (ASSOCIATED(x%cov%s)) DEALLOCATE(x%cov%s)

  x%cov%fact_chol = .FALSE.
  x%cov%equi_chol = 'N'

!-------------------------------------------------------------------------------
! 11. Clean up
!-------------------------------------------------------------------------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2state1d


SUBROUTINE ropp_fm_roprof2state2d(ro_data, x)

!****s* Copying2/ropp_fm_roprof2state_2d *
!
! NAME
!    ropp_fm_roprof2state_2d - Copy elements of an ROprof structure to a state
!                              vector.  TWO-DIMENSIONAL BACKGROUND DATA
!
! SYNOPSIS
!    type(ROprof2d)                 :: ro_data
!    type(<some state vector type>) :: x
!       ...
!    call ropp_fm_roprof2state(ro_data, x)
! 
! DESCRIPTION
!    This subroutine copies Level 2b, c and d (if applicable) data from a
!    radio occultation profile data structure into a state vector. Data is
!    is also checked for consistency. TWO-DIMENSIONAL BACKGROUND DATA.
!
! INPUTS
!   ro_data  Radio occultation profile data.
!
! OUTPUT
!   x        State vector structure.
!
! NOTES
!   Data is copied into the state vector structure without unit conversion;
!   thus, the units in the ROprof data structure for Level 2b,c, and (if
!   applicable) d must be set to the units used by the state vector 
!   variables. This can be accomplished with the ropp_fm_set_units()
!   subroutine.
!   The definition of the state vector and associated covariance matrix is 
!   model-dependent, set by level_type specified in input data.
!
! SEE ALSO
!   State2dFM
!   ropp_fm_state2roprof
!   ropp_fm_set_units
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

! 1. Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io 
  USE ropp_io_types, ONLY: ROprof2d
  USE ropp_fm
  USE ropp_fm_types, ONLY: State2dFM

  IMPLICIT NONE

  TYPE(ROprof2d),   INTENT(in)       :: ro_data   ! Input ROprof data
  TYPE(State2dFM),  INTENT(inout)    :: x         ! Output state vector

  INTEGER, DIMENSION(8)              :: DT8
  CHARACTER(len = 256)               :: level_type
  CHARACTER(len = 256)               :: routine
  CHARACTER(len = 10)                :: err_val

! 2. Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2state')

  x%state_ok = .TRUE.
  
! sizes of arrays  
  
  x%n_lev  = ro_data%Lev2b%Npoints
  x%n_horiz = ro_data%Lev2b%NHoriz

! angular spacing

  x%dtheta = ro_data%Lev2c%dtheta

! allocate the arrays containing the lats and lons
 
  ALLOCATE(x%lat(x%n_horiz))
  ALLOCATE(x%lon(x%n_horiz))

! 3. Copy geolocation and time
! --------------------------------------

  x%lon(:) = ro_data%Lev2c%lon_2d(:)
  x%lat(:) = ro_data%Lev2c%lat_2d(:)
  DT8 = (/ro_data%bg%year, ro_data%bg%month,  ro_data%bg%day, 0, &
          ro_data%bg%hour, ro_data%bg%minute, 0, 0/)
  CALL TimeSince ( DT8, x%time, 1, Base="JS2000" )
  
! 4. Check that profiles are increasing in height - 1st element towards surface
! -----------------------------------------------

!!!  call ropp_io_ascend(ro_data) 

! 5. Check and copy meteorological parameters
! --------------------------------------------

! 5.1 Check and set level numbers

  IF (ro_data%Lev2b%Npoints == 0 .OR. ro_data%Lev2b%NHoriz == 0) THEN
    CALL message(msg_warn, &
          "RO data has no level 2b (free atmospheric) parameters.")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Set status flag state%ok to FALSE")
    CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
  ENDIF
   

! 5.2 Allocate state vector variables
 
 ALLOCATE(x%geop_sfc(x%n_horiz)) 
 ALLOCATE(x%pres_sfc(x%n_horiz)) 
 ALLOCATE(x%pres(x%n_lev,x%n_horiz))
 ALLOCATE(x%temp(x%n_lev,x%n_horiz))
 ALLOCATE(x%shum(x%n_lev,x%n_horiz))
 ALLOCATE(x%geop(x%n_lev,x%n_horiz))
 ALLOCATE(x%refrac(x%n_lev,x%n_horiz))
 ALLOCATE(x%nr(x%n_lev,x%n_horiz))

! 5.3 Temperature data 
 x%temp(:,:) = ro_data%Lev2b%temp(:,:) 

! 5.4 Humidity data 
 IF (ANY(ro_data%Lev2b%shum < 0.0_wp)) THEN
   WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%shum)
   CALL message(msg_warn,"One or more humidity values are negative. ")
   CALL message(msg_info,"(Minimum humidity value = "//TRIM(err_val)// ")")
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Set status flag state%ok to FALSE")
   CALL message(msg_noin, '')
   x%state_ok = .FALSE.
   RETURN
 ENDIF
 x%shum(:,:) = ro_data%Lev2b%shum(:,:)

! 5.5 Pressure data
 IF (ANY(ro_data%Lev2b%press< 0.0_wp)) THEN
   WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%press)
   CALL message(msg_warn,"One or more pressure values are negative. ")
   CALL message(msg_info,"(Minimum pressure value = " // TRIM(err_val) // ")" )
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Set status flag state%ok to FALSE")
   CALL message(msg_noin, '')
   x%state_ok = .FALSE.
   RETURN
 ENDIF
 x%pres(:,:) = ro_data%lev2b%press(:,:)
 x%pres_sfc(:) = ro_data%Lev2c%press_sfc(:) ! surface pressure

! 5.6 Geopotential height data
 x%geop(:,:) = ro_data%lev2b%geop(:,:)
 
 x%geop_sfc(:) = ro_data%Lev2c%geop_sfc(:)
 
! 6. Define state2dFM structure: ECMWF-type hybrid levels
! -------------------------------------------------------
 
 level_type = ro_data%Lev2d%level_type
  CALL To_Upper(level_type)

 IF (INDEX(level_type,'HYBRID') > 0 .OR. &
      INDEX(level_type,'ECMWF') > 0) THEN

! 6.1 Hybrid level pressure coefficients

    ALLOCATE(x%ak(ro_data%Lev2d%Npoints)) 
    ALLOCATE(x%bk(ro_data%Lev2d%Npoints)) 
    x%ak = ro_data%Lev2d%level_coeff_a
    x%bk = ro_data%Lev2d%level_coeff_b
    IF (x%ak(ro_data%Lev2d%Npoints) == 0.0_wp)    &
         x%ak(ro_data%Lev2d%Npoints) = 1.0e-32_wp

    IF (ANY(x%ak < 0.0_wp) .OR. ALL(x%ak == 0.0_wp)) THEN
      WRITE(err_val, '(e10.2)') MINVAL(x%ak)
      CALL message(msg_warn, "Problem with level coefficients (ak) " //     &
                             "- either all zero, or one or more negative. ")
      CALL message(msg_info,"(Minimum ak value = " // TRIM(err_val) // ")")
      CALL message(msg_info, "Check input data and valid_range attributes")
      CALL message(msg_info, "Set status flag state%ok to FALSE")
      CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF
    
    IF (ANY(x%bk < 0.0_wp) .OR. ALL(x%bk == 0.0_wp)) THEN
      WRITE(err_val, '(e10.2)') MINVAL(x%bk)
       CALL message(msg_warn, "Problem with level coefficients (bk) " //     &
                           "- either all zero, or one or more negative. ")
       CALL message(msg_info, "(Minimum bk value = " // TRIM(err_val) // ")")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Set status flag state%ok to FALSE")
       CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF
        

! 6.3 Compute pressure and geopotential height on full (humidity) levels
    
    
    call ropp_fm_state2state_ecmwf(x)  

  ENDIF

! 7. Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2state2d


SUBROUTINE ropp_fm_roprof2d2state1d(ro_data, x)

!****s* Copying2/ropp_fm_roprof2d2state1d *
!
! NAME
!    ropp_fm_roprof2d2state1d - Copy elements of an 2D ROprof structure to a 
!                               1d state vector.  
!                               TWO-DIMENSIONAL BACKGROUND DATA
!
! SYNOPSIS
!    type(ROprof2d)                 :: ro_data
!    type(State1dFM)                :: x
!       ...
!    call ropp_fm_roprof2state(ro_data, x)
! 
! DESCRIPTION
!    This subroutine copies 2D Level 2b, c and d (if applicable) data from a
!    radio occultation profile data structure into a 1D state vector. Data is
!    is also checked for consistency. TWO-DIMENSIONAL BACKGROUND DATA.
!
! INPUTS
!   ro_data  Radio occultation profile data (2d).
!
! OUTPUT
!   x        State vector structure (1d).
!
! NOTES
!   Data is copied into the state vector structure without unit conversion;
!   thus, the units in the ROprof data structure for Level 2b,c, and (if
!   applicable) d must be set to the units used by the state vector 
!   variables. This can be accomplished with the ropp_fm_set_units()
!   subroutine.
!   The definition of the state vector and associated covariance matrix is 
!   model-dependent, set by level_type specified in input data.
!
! SEE ALSO
!   State1dFM
!   ropp_fm_state2roprof
!   ropp_fm_set_units
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

! 1. Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io 
  USE ropp_io_types, ONLY: ROprof2d
  USE ropp_fm
  USE ropp_fm_types, ONLY: State1dFM

  IMPLICIT NONE

  TYPE(ROprof2d),   INTENT(in)       :: ro_data   ! Input ROprof data
  TYPE(State1dFM),  INTENT(inout)    :: x         ! Output state vector

  INTEGER                            :: n_horiz
  INTEGER, DIMENSION(8)              :: DT8
  CHARACTER(len = 256)               :: level_type
  CHARACTER(len = 256)               :: routine
  CHARACTER(len = 10)                :: err_val

! 2. Error handling
! ------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_fm_roprof2state')

  x%state_ok = .TRUE.
  
! sizes of arrays  
  
  x%n_lev  = ro_data%Lev2b%Npoints
  n_horiz = ro_data%Lev2b%NHoriz

! 3. Copy geolocation and time
! --------------------------------------

  x%lon = ro_data%Lev2c%lon_2d(n_horiz/2)
  x%lat = ro_data%Lev2c%lat_2d(n_horiz/2)
  DT8 = (/ro_data%bg%year, ro_data%bg%month,  ro_data%bg%day, 0, &
          ro_data%bg%hour, ro_data%bg%minute, 0, 0/)
  CALL TimeSince ( DT8, x%time, 1, Base="JS2000" )
  
! 4. Check that profiles are increasing in height - 1st element towards surface
! -----------------------------------------------

!!!  call ropp_io_ascend(ro_data) 

! 5. Check and copy meteorological parameters
! --------------------------------------------

! 5.1 Check and set level numbers

  IF (ro_data%Lev2b%Npoints == 0 .OR. ro_data%Lev2b%NHoriz == 0) THEN
    CALL message(msg_warn, &
          "RO data has no level 2b (free atmospheric) parameters.")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Check input data and valid_range attributes")
    CALL message(msg_info, "Set status flag state%ok to FALSE")
    CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
   ENDIF

! 5.2 Allocate state vector variables

 ALLOCATE(x%pres(x%n_lev))
 ALLOCATE(x%temp(x%n_lev))
 ALLOCATE(x%shum(x%n_lev))
 ALLOCATE(x%geop(x%n_lev))

! 5.3 Temperature data 
 x%temp(:) = ro_data%Lev2b%temp(:,n_horiz/2) 

! 5.4 Humidity data 
 IF (ANY(ro_data%Lev2b%shum < 0.0_wp)) THEN
   WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%shum)
   CALL message(msg_warn,"One or more humidity values are negative. ")
   CALL message(msg_info,"(Minimum humidity value = "//TRIM(err_val)// ")")
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Set status flag state%ok to FALSE")
   CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
 ENDIF
 x%shum(:) = ro_data%Lev2b%shum(:,n_horiz/2)

! 5.5 Pressure data
 IF (ANY(ro_data%Lev2b%press< 0.0_wp)) THEN
    CALL message(msg_warn,"One or more pressure values are negative.")
    x%state_ok = .FALSE.
    RETURN
 ENDIF
 x%pres(:) = ro_data%lev2b%press(:,n_horiz/2)

! 5.6 Geopotential height data
 x%geop(:) = ro_data%lev2b%geop(:,n_horiz/2) 
 x%geop_sfc = ro_data%Lev2c%geop_sfc(n_horiz/2)
 
! 6. Define state2dFM structure: ECMWF-type hybrid levels
! -------------------------------------------------------
 
 level_type = ro_data%Lev2d%level_type
  CALL To_Upper(level_type)

 IF (INDEX(level_type,'HYBRID') > 0 .OR. &
      INDEX(level_type,'ECMWF') > 0) THEN

! 6.1 Hybrid level pressure coefficients

    ALLOCATE(x%ak(ro_data%Lev2d%Npoints)) 
    ALLOCATE(x%bk(ro_data%Lev2d%Npoints)) 
    x%ak = ro_data%Lev2d%level_coeff_a
    x%bk = ro_data%Lev2d%level_coeff_b
    IF (x%ak(ro_data%Lev2d%Npoints) == 0.0_wp)    &
         x%ak(ro_data%Lev2d%Npoints) = 1.0e-32_wp

    IF (ANY(x%ak < 0.0_wp) .OR. ALL(x%ak == 0.0_wp)) THEN
      WRITE(err_val, '(e10.2)') MINVAL(x%ak)
      CALL message(msg_warn, "Problem with level coefficients (ak) " //     &
                             "- either all zero, or one or more negative. ")
      CALL message(msg_info,"(Minimum ak value = " // TRIM(err_val) // ")")
      CALL message(msg_info, "Check input data and valid_range attributes")
      CALL message(msg_info, "Set status flag state%ok to FALSE")
      CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF
    
    IF (ANY(x%bk < 0.0_wp) .OR. ALL(x%bk == 0.0_wp)) THEN
      WRITE(err_val, '(e10.2)') MINVAL(x%bk)
       CALL message(msg_warn, "Problem with level coefficients (bk) " //     &
                           "- either all zero, or one or more negative. ")
       CALL message(msg_info, "(Minimum bk value = " // TRIM(err_val) // ")")
       CALL message(msg_info, "Check input data and valid_range attributes")
       CALL message(msg_info, "Set status flag state%ok to FALSE")
       CALL message(msg_noin, '')
       x%state_ok = .FALSE.
       RETURN
    ENDIF
        

! 6.3 Compute pressure and geopotential height on full (humidity) levels
    
    call ropp_fm_state2state_ecmwf(x)  

 IF (ANY(ro_data%Lev2b%press< 0.0_wp)) THEN
   WRITE(err_val, '(e8.1)') MINVAL(ro_data%Lev2b%press)
   CALL message(msg_warn,"One or more pressure values are negative. ")
   CALL message(msg_info,"(Minimum pressure value = " // TRIM(err_val) // ")" )
   CALL message(msg_info, "Check input data and valid_range attributes")
   CALL message(msg_info, "Set status flag state%ok to FALSE")
   CALL message(msg_noin, '')
    x%state_ok = .FALSE.
    RETURN
 ENDIF

  ENDIF

! 7. Clean up
! ------------

  CALL message_set_routine(routine)

END SUBROUTINE ropp_fm_roprof2d2state1d
