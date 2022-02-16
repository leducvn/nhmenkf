! $Id: ropp_io_success.f90 2197 2009-06-23 09:11:17Z idculv $

!****si* Test/ropp_io_success *
!
! NAME
!    ropp_io_success - Outputs test PASS/FAIL summary.
!
! SYNOPSIS
!    CALL ropp_io_success(pass, test_name, comment)
!
! DESCRIPTION
!    Writes PASS/FAIL test summary to stdout and a standard text file.
!    Called by various 'make test' routines, including ropp_<module>_compare.
!
! INPUTS
!    pass/fail, test_name, comment
!
! OUTPUT
!    Writes (to stdout and compare.txt) an overall PASS/FAIL message.
!
! AUTHOR
!    Met Office, Exeter, UK.
!    Any comments on this software should be given via the ROM SAF
!    Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_io_success(pass, test_name, comment)

  USE ropp_utils, ONLY: Get_IO_Unit

  LOGICAL,            INTENT(IN)  :: pass
  CHARACTER(LEN=*),   INTENT(IN)  :: test_name
  CHARACTER(LEN=*),   INTENT(IN)  :: comment

  CHARACTER(LEN=256), PARAMETER   :: sformat=&
                                     "('| ',a30,' | ',a30,' | ',a7,' | ',a6,' |')"

  CHARACTER(LEN=256), PARAMETER   :: compare_file='compare.txt'
  INTEGER                         :: compare_lun, iostatus

  LOGICAL                         :: exists
  CHARACTER(LEN=6)                :: spass

! 1. Write to stdout
! ------------------

  spass = '*FAIL*'

  IF ( pass ) spass = ' PASS '

  WRITE (*, '(a)') '****************************'
  WRITE (*, '(a)') '********** ' // spass // ' **********'
  WRITE (*, '(a)') '****************************'

! 2. Write same info to existing or new summary file
! --------------------------------------------------

  compare_lun = get_io_unit()  ! Find a free lun

  INQUIRE ( FILE=compare_file, EXIST=exists )

  IF ( .NOT. exists ) THEN

    OPEN ( compare_lun, FILE=compare_file, STATUS='NEW', &
           ACTION='WRITE', IOSTAT=iostatus )
    IF ( iostatus > 0 ) &
      CALL message ( msg_fatal, 'I/O error when opening summary file ' // &
                                 compare_file)

    WRITE ( compare_lun, '(a)' ) '************************* COMPARISON OF MODULE TEST RESULTS **************************'
    WRITE ( compare_lun, '(a)' ) '--------------------------------------------------------------------------------------'
    WRITE ( compare_lun, sformat ) 'Test name   ', 'Description      ', 'Run?', 'PASS?'
    WRITE ( compare_lun, '(a)' ) '--------------------------------------------------------------------------------------'

  ELSE

    OPEN ( compare_lun, FILE=compare_file, STATUS='OLD', &
           ACTION='WRITE', POSITION='APPEND', IOSTAT=iostatus )
    IF ( iostatus > 0 ) &
      CALL message ( msg_fatal, 'I/O error when opening summary file ' // &
                                 compare_file)

  ENDIF

  WRITE ( compare_lun, sformat ) TRIM(ADJUSTL(test_name)), &
                                 TRIM(ADJUSTL(comment)), 'Run', spass

  CLOSE ( compare_lun )


END SUBROUTINE ropp_io_success
