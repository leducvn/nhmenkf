! $Id: ropp_pp_read_config.f90 4452 2015-01-29 14:42:02Z idculv $

SUBROUTINE ropp_pp_read_config(file, config)

!****s* Setup/ropp_pp_read_config *
!
! NAME
!    ropp_pp_read_config - Read a configuration file for ROPP_PP
!
! SYNOPSIS
!    call ropp_pp_read_config(file, config)
! 
! DESCRIPTION
!
! INPUTS
!      file
!
! OUTPUT
!     config
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

! USE ropp_pp, not_this => ropp_pp_read_config
  USE ropp_pp_types, ONLY: PPConfig, KeyConfig
  USE messages

  IMPLICIT NONE

  CHARACTER(len = *), INTENT(in)    :: file
  TYPE(PPConfig),    INTENT(inout) :: config

  TYPE(KeyConfig)                   :: kv
  INTEGER                           :: i

  CHARACTER(len = 256)              :: routine

!-------------------------------------------------------------------------------
! 2. Error handling
!-------------------------------------------------------------------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_read_config')

!-------------------------------------------------------------------------------
! 3. Read configuration file
!-------------------------------------------------------------------------------

  CALL config_parsefile(file, kv)

!-------------------------------------------------------------------------------
! 4. Loop over all configuration elements
!-------------------------------------------------------------------------------

  DO i = 1, SIZE(kv % keys)

     SELECT CASE (kv % keys(i))

        CASE ('output_tdry')
          READ(kv % values(i), *) config % output_tdry

        CASE ('output_diag')
          READ(kv % values(i), *) config % output_diag

        CASE ('occ_method')
           config % occ_method = kv % values(i)

        CASE ('filter_method')
           config % filter_method = kv % values(i)

        CASE ('fw_go_smooth')
           READ(kv % values(i), *) config % fw_go_smooth

        CASE ('fw_go_full')
           READ(kv % values(i), *) config % fw_go_full

        CASE ('fw_wo')
           READ(kv % values(i), *) config % fw_wo

        CASE ('fw_low')
           READ(kv % values(i), *) config % fw_low

        CASE ('hmax_wo')
           READ(kv % values(i), *) config % hmax_wo
           
        CASE ('Acut')
           READ(kv % values(i), *) config % Acut

        CASE ('Pcut')
           READ(kv % values(i), *) config % Pcut

        CASE ('Bcut')
           READ(kv % values(i), *) config % Bcut

        CASE ('Hcut')
           READ(kv % values(i), *) config % Hcut

        CASE ('CFF')
           READ(kv % values(i), *) config % cff

        CASE ('dsh')
           READ(kv % values(i), *) config % dsh

        CASE ('opt_DL2')
           READ(kv % values(i), *) config % opt_DL2

        CASE ('opt_spectra')
           READ(kv % values(i), *) config % opt_spectra

        CASE ('egm96')
           config % egm96 = kv % values(i)

        CASE ('corr_egm96')
           config % corr_egm96 = kv % values(i)

        CASE ('navbit_file')
           config % navbit_file = kv % values(i)

        CASE ('method')
           config % method = kv % values(i)

         CASE ('so_method')
           config % so_method = kv % values(i)

        CASE ('abel')
           config % abel = kv % values(i)

        CASE ('mfile')
           config % mfile = kv % values(i)

        CASE ('bfile')
           config % bfile = kv % values(i)

        CASE ('dpi')
           READ(kv % values(i), *) config % dpi

        CASE ('np_smooth')
           READ(kv % values(i), *) config % np_smooth

        CASE ('fw_smooth')
           READ(kv % values(i), *) config % fw_smooth

        CASE ('nparm_fit')
           READ(kv % values(i), *) config % nparm_fit

        CASE ('hmin_fit')
           READ(kv % values(i), *) config % hmin_fit

        CASE ('hmax_fit')
           READ(kv % values(i), *) config % hmax_fit

        CASE ('omega_fit')
           READ(kv % values(i), *) config % omega_fit

        CASE ('f_width')
           READ(kv % values(i), *) config % f_width

        CASE ('delta_p')
           READ(kv % values(i), *) config % delta_p

        CASE ('s_smooth')
           READ(kv % values(i), *) config % s_smooth

         CASE ('z_ion')
            READ(kv % values(i), *) config % z_ion

         CASE ('z_str')
            READ(kv % values(i), *) config % z_str

         CASE ('z_ltr')
            READ(kv % values(i), *) config % z_ltr

         CASE ('n_smooth')
           READ(kv % values(i), *) config % n_smooth

         CASE ('model_err')
           READ(kv % values(i), *) config % model_err

         CASE ('ztop_invert')
           READ(kv % values(i), *) config % ztop_invert

         CASE ('dzh_invert')
           READ(kv % values(i), *) config % dzh_invert

         CASE ('dzr_invert')
           READ(kv % values(i), *) config % dzr_invert

!-------------------------------------------------------------------------------
! 5. Unknown configuration item
!-------------------------------------------------------------------------------

        CASE default
           CALL message(msg_warn,                  &
                "Unknown configuration item: " //  &
                TRIM(kv % keys(i)) // " - skipping.")

     END SELECT

  ENDDO

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(kv%keys, kv%values)
  CALL message_set_routine(routine)


CONTAINS


!-------------------------------------------------------------------------------
! 7. Read and interpret PP config text file
!-------------------------------------------------------------------------------

  SUBROUTINE config_parsefile(file, keyvals)

    USE ropp_pp
    
    IMPLICIT NONE
    
    CHARACTER(len = *), INTENT(in)                :: file
    TYPE(KeyConfig),    INTENT(out)               :: keyvals
    
    CHARACTER(len = 1024), DIMENSION(:), POINTER  :: lines => null()
    CHARACTER(len = 1024), DIMENSION(2)           :: string_val  
    INTEGER                                       :: nlines, i
  
! 7.1 Read PP config file, ignoring comment lines
    CALL config_readfile(file, lines, comments = '#!;*')
  
! 7.2 Split data lines into keys and values
    nlines = SIZE(lines)
    ALLOCATE(keyvals%keys(nlines), keyvals%values(nlines))
    
    DO i = 1, nlines
       CALL get_strings(lines(i), string_val)
       
       keyvals%keys(i)   = TRIM(ADJUSTL(string_val(1)))
       keyvals%values(i) = TRIM(ADJUSTL(string_val(2)))
       
    ENDDO
  
    DEALLOCATE(lines)
     
  END SUBROUTINE config_parsefile

!-------------------------------------------------------------------------------
! 8. Read PP config text file, excluding empty lines and comments
!-------------------------------------------------------------------------------

  SUBROUTINE config_readfile(file, lines, comments)

    USE ropp_utils
    USE arrays, ONLY: WHERE, &
                    & reallocate
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
    
! 8.1 Open the data file
    
    INQUIRE(file = file, exist = exist)
    
    IF (exist) THEN
       unit = get_io_unit()
       OPEN(unit, file = file)
    ELSE
       CALL message(msg_fatal, 'File does not exist:' // TRIM(file) // '.\n')
    ENDIF
    
! 8.2 Loop over all lines

    i = 1 
    j = 1
    ALLOCATE(line_arr(chunksize))

    DO
       READ(unit, '(a)', iostat = iostat) line
       IF (iostat > 0) THEN
          CALL message(msg_warn,   &
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
  
! 8.3 Close data file
  
  CLOSE(unit)
  
! 8.4 Check for comments
  
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

! 8.5 Copy found lines to output array

  idx => WHERE(lengths > 0, nlines)
  ALLOCATE(lines(nlines))
  DO i = 1,nlines
     lines(i) = line_arr(idx(i))
  ENDDO
  
! 8.6 Clean up

  DEALLOCATE(lengths, idx, line_arr)
  
END SUBROUTINE config_readfile


!-------------------------------------------------------------------------------
! 9. Split strings into two sub-strings, separated by '='
!-------------------------------------------------------------------------------

SUBROUTINE get_strings(str, output) 
  
  IMPLICIT NONE
  
  CHARACTER(len = *),               INTENT(in)  :: str
  CHARACTER(len = *), DIMENSION(:), INTENT(out) :: output
  INTEGER                                       :: idx

  idx = INDEX(str,'=')
  output(1) = str(1:idx-1)
  output(2) = TRIM(str(idx+1:LEN_TRIM(str)))

END SUBROUTINE get_strings

END SUBROUTINE ropp_pp_read_config
