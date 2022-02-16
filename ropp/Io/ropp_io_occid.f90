! $Id: ropp_io_occid.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Initialisation/ropp_io_occid *
!
! NAME
!    ropp_io_occid - Generate an Occultation ID suitable for a file name.
!
! SYNOPSIS
!    use ropp_io
!    type(ROprof) :: ROdata
!     ...
!    call ropp_io_occid(ROdata)
!
! DESCRIPTION
!    This subroutine generates a unique occultation identifier of the form:
!
!       tt_yyyymmddhhmmss_llll_gggg_cccc
!
!    where:
!
!       tt             is 'OC' for measured or retrieved occultation
!                      or 'BG' for background data (e.g. from an NWP model)
!       yyyymmddhhmmss is the date & time of the occultation
!       llll           is the LEO  (Rx) 4-chr identifier
!       gggg           is the GNSS (Tx) 4-chr identifier
!       cccc           is the processing centre identifier (1st 4 chrs,
!                      padded with underscores if neccessary)
!
!    All these elements are as extracted from the RO derived type, and
!    the occultion ID is returned in that structure. Basic Q/C is applied
!    to the elements used for the occultaton ID only to ensure correct
!    formatting of the ID string; the elements themselves are not changed.
!
! INPUT
!   ROdata  dtyp  RO data structure (derived type)
!
! OUTPUT
!   ROdata  dtyp  RO data structure (derived type)
!
! USES
!  ropp_io
!
! CALLS
!   To_Upper
!
! NOTES
!   The subroutines requires that the PCD_occultation flag (specifying whether
!   the data originates from an occultation measurement and / or retrieval or
!   background data is set properly.
!
! REFERENCES
!    ROPP User Guide. Part I. Ref: SAF/ROM/METO/UG/ROPP/002
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
! 1. 1d version
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_occid_1d(ROdata)

! Modules

  USE ropp_utils, ONLY: To_Upper
! USE ropp_io, not_this => ropp_io_occid_1d
  USE ropp_io_types, ONLY: ROprof,          &
                           PCD_Occultation, &
                           PCD_Missing   

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN =  *), PARAMETER :: validchars = &
                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_" ! Valid ID chrs

! Argument list parameters

  TYPE(ROprof), INTENT(INOUT) :: ROdata


! Local variables

  CHARACTER (LEN=2)  :: Pty                  ! Profile type (OC or BG)
  CHARACTER (LEN=4)  :: Lid, Gid, Pid        ! Local LEO, GNSS & Proc.Centre IDs
  CHARACTER (LEN=14) :: CDT                  ! Occultation date time chr string
  INTEGER            :: Year, Month,  Day    ! Local occ date
  INTEGER            :: Hour, Minute, Second ! Local occ time
  INTEGER            :: i                    ! loop counter
  INTEGER            :: num                  ! Temporary number
  INTEGER            :: ierr                 ! I/O error

!-------------------------------------------------------------------------------
! 1.2 Extract PCD flag for profile type
!-------------------------------------------------------------------------------

   IF ( .NOT. BTEST(ROdata%PCD, PCD_missing) .AND. &
              BTEST(ROdata%PCD, PCD_occultation) ) THEN
      Pty = "BG"
   ELSE
      Pty = "OC"
   END IF

!-------------------------------------------------------------------------------
! 1.3 Extract occultation date time to yyyymmddhhmmss format with Q/C checks
!-------------------------------------------------------------------------------

  Year   = ROdata%DTocc%Year
  Month  = ROdata%DTocc%Month
  Day    = ROdata%DTocc%Day
  Hour   = ROdata%DTocc%Hour
  Minute = ROdata%DTocc%Minute
  Second = ROdata%DTocc%Second

  IF ( Year   >    0 .AND. Year   <  100 ) Year   = Year + 2000
  IF ( Year   < 1000 .OR.  Year   > 9999 ) Year   = 9999
  IF ( Month  <    1 .OR.  Month  >   12 ) Month  = 99
  IF ( Day    <    1 .OR.  Day    >   31 ) Day    = 99
  IF ( Hour   <    0 .OR.  Hour   >   23 ) Hour   = 99
  IF ( Minute <    0 .OR.  Minute >   59 ) Minute = 99
  IF ( Second <    0 .OR.  Second >   59 ) Second = 99

  WRITE ( CDT, fmt = "(I4.4,5I2.2)" ) Year, Month,  Day, &
                                      Hour, Minute, Second

!-------------------------------------------------------------------------------
! 1.4 Check LEO ID; make uppercase & replace any invalid chrs (not A-L,0-9)
!     with underscores
!-------------------------------------------------------------------------------

  Lid = ROdata%LEO_ID
  CALL To_Upper ( Lid )
  DO i = 1, 4
    IF ( INDEX(validchars,Lid(i:i)) == 0 ) Lid(i:i) = "_"
  END DO
  IF ( Lid == "____" ) Lid = "UNKN"

!-------------------------------------------------------------------------------
! 1.5 Check GNSS ID; check validity of
!     - GNSS Tx class is alphabetic  [A-Z] (subs, 'U' if not)
!     - SV no. is wholly numeric [001-999] (subs. '999' if not)
!-------------------------------------------------------------------------------

  Gid = ROdata%GNS_ID
  CALL To_Upper ( Gid )
  IF ( INDEX(validchars(1:26),Gid(1:1)) == 0 ) Gid(1:1) = "U"

  READ ( Gid(2:4), FMT="(I3)", IOSTAT=ierr ) num
  IF ( ierr /= 0 .OR. &
       num  < 001 .OR. num > 999 ) num = 999
  WRITE ( Gid(2:4), FMT="(I3.3)" ) num

!-------------------------------------------------------------------------------
! 1.6 Check Originator ID (first 4 chrs); make uppercase & replace
!     any invalid chrs (not A-L,0-9) with underscores
!-------------------------------------------------------------------------------

  Pid = ROdata%Processing_centre(1:4)
  CALL To_Upper ( Pid )
  DO i = 1, 4
    IF ( INDEX(validchars,Pid(i:i)) == 0 ) Pid(i:i) = "_"
  END DO
  IF ( Pid == "____" ) Pid = "UNKN"

!-------------------------------------------------------------------------------
! 1.7 Generate occultation ID string
!-------------------------------------------------------------------------------

  ROdata%occ_id = Pty // "_" // CDT // "_" // &
                  Lid // "_" // Gid // "_" // Pid

END SUBROUTINE ropp_io_occid_1d

!-------------------------------------------------------------------------------
! 2. 2d version
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_occid_2d(ROdata)

! Modules

  USE ropp_utils, ONLY: To_Upper
! USE ropp_io, not_this => ropp_io_occid_2d
  USE ropp_io_types, ONLY: ROprof2d,          &
                           PCD_Occultation, &
                           PCD_Missing   

  IMPLICIT NONE

! Fixed parameters

  CHARACTER (LEN =  *), PARAMETER :: validchars = &
                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_" ! Valid ID chrs

! Argument list parameters

  TYPE(ROprof2d), INTENT(INOUT) :: ROdata


! Local variables

  CHARACTER (LEN=2)  :: Pty                  ! Profile type (OC or BG)
  CHARACTER (LEN=4)  :: Lid, Gid, Pid        ! Local LEO, GNSS & Proc.Centre IDs
  CHARACTER (LEN=14) :: CDT                  ! Occultation date time chr string
  INTEGER            :: Year, Month,  Day    ! Local occ date
  INTEGER            :: Hour, Minute, Second ! Local occ time
  INTEGER            :: i                    ! loop counter
  INTEGER            :: num                  ! Temporary number
  INTEGER            :: ierr                 ! I/O error

!-------------------------------------------------------------------------------
! 2.1 Extract PCD flag for profile type
!-------------------------------------------------------------------------------

   IF ( .NOT. BTEST(ROdata%PCD, PCD_missing) .AND. &
              BTEST(ROdata%PCD, PCD_occultation) ) THEN
      Pty = "BG"
   ELSE
      Pty = "OC"
   END IF

!-------------------------------------------------------------------------------
! 2.2 Extract occultation date time to yyyymmddhhmmss format with Q/C checks
!-------------------------------------------------------------------------------

  Year   = ROdata%DTocc%Year
  Month  = ROdata%DTocc%Month
  Day    = ROdata%DTocc%Day
  Hour   = ROdata%DTocc%Hour
  Minute = ROdata%DTocc%Minute
  Second = ROdata%DTocc%Second

  IF ( Year   >    0 .AND. Year   <  100 ) Year   = Year + 2000
  IF ( Year   < 1000 .OR.  Year   > 9999 ) Year   = 9999
  IF ( Month  <    1 .OR.  Month  >   12 ) Month  = 99
  IF ( Day    <    1 .OR.  Day    >   31 ) Day    = 99
  IF ( Hour   <    0 .OR.  Hour   >   23 ) Hour   = 99
  IF ( Minute <    0 .OR.  Minute >   59 ) Minute = 99
  IF ( Second <    0 .OR.  Second >   59 ) Second = 99

  WRITE ( CDT, fmt = "(I4.4,5I2.2)" ) Year, Month,  Day, &
                                      Hour, Minute, Second

!-------------------------------------------------------------------------------
! 2.3 Check LEO ID; make uppercase & replace any invalid chrs (not A-L,0-9)
!     with underscores
!-------------------------------------------------------------------------------

  Lid = ROdata%LEO_ID
  CALL To_Upper ( Lid )
  DO i = 1, 4
    IF ( INDEX(validchars,Lid(i:i)) == 0 ) Lid(i:i) = "_"
  END DO
  IF ( Lid == "____" ) Lid = "UNKN"

!-------------------------------------------------------------------------------
! 2.4 Check GNSS ID; check validity of
!     - GNSS Tx class is alphabetic  [A-Z] (subs, 'U' if not)
!     - SV no. is wholly numeric [001-999] (subs. '999' if not)
!-------------------------------------------------------------------------------

  Gid = ROdata%GNS_ID
  CALL To_Upper ( Gid )
  IF ( INDEX(validchars(1:26),Gid(1:1)) == 0 ) Gid(1:1) = "U"

  READ ( Gid(2:4), FMT="(I3)", IOSTAT=ierr ) num
  IF ( ierr /= 0 .OR. &
       num  < 001 .OR. num > 999 ) num = 999
  WRITE ( Gid(2:4), FMT="(I3.3)" ) num

!-------------------------------------------------------------------------------
! 2.5 Check Originator ID (first 4 chrs); make uppercase & replace
!     any invalid chrs (not A-L,0-9) with underscores
!-------------------------------------------------------------------------------

  Pid = ROdata%Processing_centre(1:4)
  CALL To_Upper ( Pid )
  DO i = 1, 4
    IF ( INDEX(validchars,Pid(i:i)) == 0 ) Pid(i:i) = "_"
  END DO
  IF ( Pid == "____" ) Pid = "UNKN"

!-------------------------------------------------------------------------------
! 2.6 Generate occultation ID string
!-------------------------------------------------------------------------------

  ROdata%occ_id = Pty // "_" // CDT // "_" // &
                  Lid // "_" // Gid // "_" // Pid

END SUBROUTINE ropp_io_occid_2d
