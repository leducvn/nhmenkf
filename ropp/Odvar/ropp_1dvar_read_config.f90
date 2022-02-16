! $Id: ropp_1dvar_read_config.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* 1DVar/ropp_1dvar_read_config *
!
! NAME
!    ropp_1dvar_read_config - Read a configuration file for a 1DVar retrieval.
!
! SYNOPSIS
!    call ropp_1dvar_read_config(file, config)
! 
! DESCRIPTION
!
! INPUTS
!      file        - configuration file name
!
! OUTPUT
!     config       - configuration options
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

SUBROUTINE ropp_1dvar_read_config(file, config)

!----------------------------------------------------------------------------
! 1. Declarations
!----------------------------------------------------------------------------

  USE ropp_1dvar
  USE messages

  IMPLICIT NONE

  CHARACTER(len = *), INTENT(in)    :: file
  TYPE(VarConfig),    INTENT(inout) :: config

  TYPE(KeyConfig)                   :: kv
  INTEGER                           :: i
  INTEGER                           :: season_range_error = 0

  CHARACTER(len = 256)              :: routine

!----------------------------------------------------------------------------
! 2. Error handling
!----------------------------------------------------------------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_1dvar_read_config')

!----------------------------------------------------------------------------
! 3. Read configuration file
!----------------------------------------------------------------------------

  CALL config_parsefile(file, kv)
 
!----------------------------------------------------------------------------
! 4. Loop over all configuration elements
!----------------------------------------------------------------------------

  DO i = 1, SIZE(kv % keys)

     SELECT CASE (kv % keys(i))

!----------------------------------------------------------------------------
! 5. Input and output data files
!----------------------------------------------------------------------------

        CASE ('bg_input_file')
           config % bg_file = kv % values(i)

        CASE ('bg_corr_file')
           config % bg_corr_file = kv % values(i)

        CASE ('obs_input_file')
           config % obs_file = kv % values(i)

        CASE ('obs_corr_file')
           config % obs_file = kv % values(i)

        CASE ('output_file')
           config % out_file = kv % values(i)

!----------------------------------------------------------------------------
! 6. Error covariance models / methods
!----------------------------------------------------------------------------

        CASE ('bg_covar_method')
           config % bg_covar_method  = kv % values(i)

        CASE ('obs_covar_method')
           config % obs_covar_method = kv % values(i)

!----------------------------------------------------------------------------
! 7. Observation cutoff range
!----------------------------------------------------------------------------
        
         CASE ('min_1dvar_height')
           READ(kv % values(i), *) config % min_1dvar_height

        CASE ('max_1dvar_height')
           READ(kv % values(i), *) config % max_1dvar_height

!----------------------------------------------------------------------------
! 8. Generic quality control
!----------------------------------------------------------------------------

        CASE ('genqc_colocation_apply')
           READ(kv % values(i), *) config % genqc_colocation_apply

        CASE ('genqc_max_distance')
           READ(kv % values(i), *) config % genqc_max_distance

        CASE ('genqc_max_time_sep')
           READ(kv % values(i), *) config % genqc_max_time_sep
         
        CASE ('genqc_min_temperature')
           READ(kv % values(i), *) config % genqc_min_temperature

        CASE ('genqc_max_temperature')
           READ(kv % values(i), *) config % genqc_max_temperature

        CASE ('genqc_min_spec_humidity')
           READ(kv % values(i), *) config % genqc_min_spec_humidity

        CASE ('genqc_max_spec_humidity')
           READ(kv % values(i), *) config % genqc_max_spec_humidity

        CASE ('genqc_min_impact')
           READ(kv % values(i), *) config % genqc_min_impact

        CASE ('genqc_max_impact')
           READ(kv % values(i), *) config % genqc_max_impact

        CASE ('genqc_min_bangle')
           READ(kv % values(i), *) config % genqc_min_bangle

        CASE ('genqc_max_bangle')
           READ(kv % values(i), *) config % genqc_max_bangle

        CASE ('genqc_min_geop_refrac')
           READ(kv % values(i), *) config % genqc_min_geop_refrac

        CASE ('genqc_max_geop_refrac')
           READ(kv % values(i), *) config % genqc_max_geop_refrac

        CASE ('genqc_min_refractivity')
           READ(kv % values(i), *) config % genqc_min_refractivity

        CASE ('genqc_max_refractivity')
           READ(kv % values(i), *) config % genqc_max_refractivity
           
         CASE ('genqc_min_obheight')
           READ(kv % values(i), *) config % genqc_min_obheight

!----------------------------------------------------------------------------
! 9. Background quality control
!----------------------------------------------------------------------------

        CASE ('bgqc_apply')
           READ(kv % values(i), *) config % bgqc_apply

        CASE ('bgqc_reject_factor')
           READ(kv % values(i), *) config % bgqc_reject_factor

        CASE ('bgqc_reject_max_percent')
           READ(kv % values(i), *) config % bgqc_reject_max_percent

!----------------------------------------------------------------------------
! 10. Probability of Gross Error quality control
!----------------------------------------------------------------------------

        CASE ('pge_apply')
           READ(kv % values(i), *) config % pge_apply

        CASE ('pge_fg')
           READ(kv % values(i), *) config % pge_fg

        CASE ('pge_d')
           READ(kv % values(i), *) config % pge_d

!----------------------------------------------------------------------------
! 11. Minimizer configuration 
!----------------------------------------------------------------------------

        CASE ('minropp_method')
           config % minropp % method = kv % values(i)

        CASE ('minropp_log_file')
           config % minropp % log_file = kv % values(i)

        CASE ('minropp_impres')
           READ(kv % values(i), *) config % minropp % impres

        CASE ('minropp_n_iter')
           READ(kv % values(i), *) config % minropp % n_iter

        CASE ('minropp_mode')
           READ(kv % values(i), *) config % minropp % mode

        CASE ('minropp_n_updates')
           READ(kv % values(i), *) config % minropp % n_updates

        CASE ('minropp_eps_grad')
           READ(kv % values(i), *) config % minropp % eps_grad

        CASE ('minropp_dx_min')
           READ(kv % values(i), *) config % minropp % dx_min

!----------------------------------------------------------------------------
! 12. Addtional convergence checks
!----------------------------------------------------------------------------

        CASE ('conv_check_apply')
           READ(kv % values(i), *) config % conv_check_apply

        CASE ('conv_check_n_previous')
           READ(kv % values(i), *) config % conv_check_n_previous

        CASE ('conv_check_max_delta_state')
           READ(kv % values(i), *) config % conv_check_max_delta_state

        CASE ('conv_check_max_delta_J')
           READ(kv % values(i), *) config % conv_check_max_delta_J

!----------------------------------------------------------------------------
! 13. Preconditioning 
!----------------------------------------------------------------------------

        CASE ('use_precond')
           READ(kv % values(i), *) config % use_precond

!----------------------------------------------------------------------------
! 14. Additional output
!----------------------------------------------------------------------------

        CASE ('extended_1dvar_diag')
           READ(kv % values(i), *) config % extended_1dvar_diag

!----------------------------------------------------------------------------
! 15. Log(pressure) and Log(humidity) options
!----------------------------------------------------------------------------

        CASE ('use_logq')
           READ(kv % values(i), *) config % use_logq

        CASE ('use_logp')
           READ(kv % values(i), *) config % use_logp

!----------------------------------------------------------------------------
! 16. Seasonal observation error scaling
!----------------------------------------------------------------------------

        CASE ('season_amp')
           READ(kv % values(i), *) config % season_amp
           IF (.NOT. (config % season_amp .GE.   0.0_wp .AND. &
                      config % season_amp .LE.  50.0_wp)) THEN
             season_range_error = season_range_error + 1
           END IF

        CASE ('season_offset')
           READ(kv % values(i), *) config % season_offset
           IF (.NOT. (config % season_offset .GE.   0.0_wp .AND. &
                      config % season_offset .LE. 500.0_wp)) THEN
             season_range_error = season_range_error + 1
           END IF

        CASE ('season_phase')
           READ(kv % values(i), *) config % season_phase
           IF (.NOT. (config % season_phase .GE. -1.0_wp .AND. &
                      config % season_phase .LE.  1.0_wp)) THEN
             season_range_error = season_range_error + 1
           END IF

!----------------------------------------------------------------------------
! 17. Unknown configuration item
!----------------------------------------------------------------------------

        CASE default
           CALL message(msg_warn,                  &
                "Unknown configuration item: " //  &
                TRIM(kv % keys(i)) // " - skipping.")

     END SELECT

  ENDDO

  ! If any specified seasonal parameters are out of range, use defaults.
  IF (season_range_error .GT. 0) THEN
    config % season_amp    = 0.0_wp
    config % season_offset = 0.0_wp
    config % season_phase  = 0.0_wp
    CALL message(msg_warn,     &
                 'Invalid seasonal observation error scaling parameters'// &
                 ' - switching off seasonal dependence. \n')
  END IF

!----------------------------------------------------------------------------
! 18. Clean up
!----------------------------------------------------------------------------

  CALL message_set_routine(routine)


CONTAINS


!----------------------------------------------------------------------------
! 19. Read and interpret 1dVar config text file
!----------------------------------------------------------------------------

  SUBROUTINE config_parsefile(file, keyvals)

    USE ropp_1dvar
    
    IMPLICIT NONE
    
    CHARACTER(len = *), INTENT(in)                :: file
    TYPE(KeyConfig),    INTENT(out)               :: keyvals
    
    CHARACTER(len = 1024), DIMENSION(:), POINTER  :: lines => null()
    CHARACTER(len = 1024), DIMENSION(2)           :: string_val  
    INTEGER                                       :: nlines, i
  
! 19.1 Read 1dVar config file, ignoring comment lines
    CALL config_readfile(file, lines, comments = '#!;*')
  
! 19.2 Split data lines into keys and values
    nlines = SIZE(lines)
    ALLOCATE(keyvals%keys(nlines), keyvals%values(nlines))
    
    DO i = 1, nlines
       CALL get_strings(lines(i), string_val)
       
       keyvals%keys(i)   = TRIM(ADJUSTL(string_val(1)))
       keyvals%values(i) = TRIM(ADJUSTL(string_val(2)))
       
    ENDDO
  
    DEALLOCATE(lines)
     
  END SUBROUTINE config_parsefile

!----------------------------------------------------------------------------
! 20. Read 1dVar config text file, excluding empty lines and comments
!----------------------------------------------------------------------------

  SUBROUTINE config_readfile(file, lines, comments)

    USE ropp_utils
    USE arrays 
    USE messages
    
    IMPLICIT NONE
    
    CHARACTER(len = *),               INTENT(in) :: file
    CHARACTER(len = *), DIMENSION(:), POINTER    :: lines 
    CHARACTER(len = *),               OPTIONAL   :: comments

    CHARACTER(len = 4096), DIMENSION(:), POINTER :: line_arr   => null()
    
    CHARACTER(len = 4096)                        :: line

    INTEGER                                      :: i, j, k, iostat
    INTEGER                                      :: ix, iln
    INTEGER, PARAMETER                           :: chunksize = 4096
    LOGICAL                                      :: exist
    
    INTEGER                                      :: unit, nlines
    INTEGER, DIMENSION(:), ALLOCATABLE           :: lengths
    INTEGER, DIMENSION(:), POINTER               :: idx
    
! 20.1 Open the data file
    
    INQUIRE(file = file, exist = exist)
    
    IF (exist) THEN
       unit = get_io_unit()
       OPEN(unit, file = file)
    ELSE
       CALL message(msg_fatal, 'File does not exist:' // TRIM(file) // '.\n')
    ENDIF
    
! 20.2 Loop over all lines

    i = 1 
    j = 1
    ALLOCATE(line_arr(chunksize))

    DO
       READ(unit, '(a)', iostat = iostat) line
       IF (iostat > 0) THEN
          CALL message(msg_warn,     &
                       'Error during config read - aborting the loop. \n')
          EXIT
       END IF
       IF (iostat < 0) THEN
          EXIT
     END IF

     IF (i <= chunksize) THEN
        line_arr(j) = line
        i = i + 1 ; j = j + 1
     ELSE
        CALL reallocate(line_arr, int(SIZE(line_arr)) + chunksize)
        line_arr(j) = line
        i = 2 ; j = j + 1
     END IF
     
  END DO
  
  nlines = j - 1
  
! The following line segfaults with the Intel v8.1 and v9.0 compilers:
!!$  call reallocate(line_arr, n)
  
! 20.3 Close data file
  
  CLOSE(unit)
  
! 20.4 Check for comments
  
  ALLOCATE(lengths(nlines))
  
  IF (PRESENT(comments)) THEN
     DO i = 1, nlines
        DO j = 1, LEN_TRIM(comments)
           ix = INDEX(line_arr(i), comments(j:j))
           IF (ix > 0) THEN
              iln = LEN(line_arr(i))
              DO k = ix, iln
                 line_arr(i)(k:k) = ' '
              ENDDO
           ENDIF
        ENDDO
        lengths(i) = LEN_TRIM(line_arr(i))
     ENDDO
  ELSE
     DO i = 1, nlines
        lengths(i) = LEN_TRIM(line_arr(i))
     ENDDO
     lengths = 1
  ENDIF

! 20.5 Copy found lines to output array

  idx => WHERE(lengths > 0, nlines)
  ALLOCATE(lines(nlines))
  DO i = 1,nlines
     lines(i) = line_arr(idx(i))
  ENDDO
  
! 20.6 Clean up

  DEALLOCATE(lengths, idx, line_arr)
  
END SUBROUTINE config_readfile


!----------------------------------------------------------------------------
! 21. Split strings into two sub-strings, separated by '='
!----------------------------------------------------------------------------

SUBROUTINE get_strings(str, output) 
  
  IMPLICIT NONE
  
  CHARACTER(len = *),               INTENT(in)  :: str
  CHARACTER(len = *), DIMENSION(:), INTENT(out) :: output
  INTEGER                                       :: idx

  idx = INDEX(str,'=')
  output(1) = str(1:idx-1)
  output(2) = TRIM(str(idx+1:LEN_TRIM(str)))

END SUBROUTINE get_strings



END SUBROUTINE ropp_1dvar_read_config
