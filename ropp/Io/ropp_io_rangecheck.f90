! $Id: ropp_io_rangecheck.f90 4452 2015-01-29 14:42:02Z idculv $
!
!****s* Initialisation/ropp_io_rangecheck *
!
! NAME
!   ropp_io_rangecheck - Check parameter ranges within an RO derived type
!
! SYNOPSIS
!   TYPE(ROprof) :: ROdata
!   CALL ropp_io_rangecheck ( ROdata )
!
!       - or -
!
!   CALL ropp_io_rangecheck ( ROdata%Lev1a)
!   CALL ropp_io_rangecheck ( ROdata%Lev1b)
!   CALL ropp_io_rangecheck ( ROdata%Lev2a)
!   CALL ropp_io_rangecheck ( ROdata%Lev2b)
!   CALL ropp_io_rangecheck ( ROdata%Lev2c)
!   CALL ropp_io_rangecheck ( ROdata%Lev2d)
!
! INPUT
!   ROdata     struc  RO profile structure
!
! OUTPUT
!   ROdata     struc  RO profile structure with invalid values
!                     set to appropriate 'missing data' values
!
! CALLS
!   ropp_io_occid
!   where
!
! USES
!   typesizes
!   ropp_io_types
!   ropp_io
!   ropp_utils
!
! DESCRIPTION
!   This subroutine check the values of all core numeric ROPP
!   parameters against their range attributes and if out-of-range,
!   substitutes an appropriate missing data flag value.
!   Longitudes in the range 180 to 360 deg. are treated as valid
!   but are converted to the range -180 to 0 deg.
!   Character parameters are made upper case, with some restricted
!   to alpha-numeric and underscore characters.
!   Profile levels having missing coordinates (time offset for
!   level 1a, Impact Parameter for L1b, geometric height for L2a or
!   geopotential height for L2c & L2d) are filtered out.
!   Profiles having at least one valid coordinate but no valid
!   observational data (e.g. SNR/phase, BA, REF or T/q/P) are flagged
!   as 'missing'.
!
! TODO
!   Set profile quality value missing (or zero) when critical profile
!   parameter(s) missing at any level
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

!----------------------------------------------------------------
! 1. 1d ROprof structure
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_1d ( ROdata ) ! (inout)

! Modules

  USE typesizes,     ONLY: wp => EightByteReal
  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_1d
  USE ropp_utils,    ONLY : To_Upper, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

  CHARACTER (LEN =  *), PARAMETER :: validchars = &
                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_- " ! Valid chr set

! Argument list parameters

  TYPE(ROprof), INTENT(INOUT) :: ROdata    ! RO data structure

! Local variables

  INTEGER  :: i    ! Loop indices
! INTEGER  :: ikl  ! Keyword length    ! Commented at 20 July, 2016
! INTEGER  :: ict  ! Invalid chr count ! Commented at 20 July, 2016
  INTEGER  :: num  ! Temporary number
  INTEGER  :: ierr ! I/O error
  REAL     :: tmp  ! Temporary number

! 1.1 Character variables; ensure left-justified and uppercase letters;
!     replace 'illegal' characters with underscores, or check numeric
!     representations are valid. If the whole field was invalid, set
!     as unknown

! The following strings are fixed 4-character keywords with
! restricted character set [A-Z],[0-9],'_','-' only.

  ROdata%LEO_ID = ADJUSTL(ROdata%LEO_ID)
  CALL To_Upper ( ROdata%LEO_ID )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%LEO_ID(i:i)) == 0 ) &
      ROdata%LEO_ID(i:i) = "_"
  END DO
  IF ( INDEX(ROdata%LEO_ID(1:3), "_") > 0 )  ROdata%LEO_ID = "UNKN"

  ROdata%GNS_ID = ADJUSTL(ROdata%GNS_ID)
  CALL To_Upper ( ROdata%GNS_ID(1:1) )
  IF ( INDEX(validchars(1:26),ROdata%GNS_ID(1:1)) == 0 ) &
    ROdata%GNS_ID(1:1) = "U"
  READ  ( ROdata%GNS_ID(2:4), FMT="(I3)", IOSTAT=ierr ) num
  IF ( ierr /= 0 .OR. &
       num  < 001 .OR. num > 999 ) num = 999
  WRITE ( ROdata%GNS_ID(2:4), FMT="(I3.3)" ) num

  ROdata%STN_ID = ADJUSTL(ROdata%STN_ID)
  CALL To_Upper ( ROdata%STN_ID )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%STN_ID(i:i)) == 0 ) &
      ROdata%STN_ID(i:i) = "_"
  END DO
  IF ( INDEX(ROdata%STN_ID(1:3), "_") > 0 )  ROdata%STN_ID = "UNKN"

  CALL keyword_check(ROdata%Processing_centre)

  CALL keyword_check(ROdata%Processing_software)

  CALL keyword_check(ROdata%POD_Method)

  CALL keyword_check(ROdata%Phase_Method)

  CALL keyword_check(ROdata%Bangle_Method)

  CALL keyword_check(ROdata%Refrac_Method)

  CALL keyword_check(ROdata%Meteo_Method)

  CALL keyword_check(ROdata%Thin_Method)

! Ensure software version string conforms to "Vnn.nnn" format

  ROdata%Software_Version = ADJUSTL(ROdata%Software_Version)
  CALL To_Upper ( ROdata%Software_Version )
  IF ( ROdata%Software_Version(2:2) == " " ) &
       ROdata%Software_Version(2:2) =  "0"
  ROdata%Software_Version(8:) = " "
  READ ( ROdata%Software_Version(2:7), FMT=*, IOSTAT=ierr ) tmp
  IF ( ROdata%Software_Version(1:1) /= "V" .OR. &
       ierr                         /=  0  .OR. &
       tmp < 0.0 .OR. tmp > 99.999 ) ROdata%Software_Version = "UNKNOWN"

! 1.2 Occultation date/time

  IF ( ROdata%DTocc%Year   < ROdata%DTocc%Range%Year(mn) .OR. &
       ROdata%DTocc%Year   > ROdata%DTocc%Range%Year(mx) )    &
       ROdata%DTocc%Year   = 9999

  IF ( ROdata%DTocc%Month  < ROdata%DTocc%Range%Month(mn) .OR. &
       ROdata%DTocc%Month  > ROdata%DTocc%Range%Month(mx) )    &
       ROdata%DTocc%Month  = 99

  IF ( ROdata%DTocc%Day    < ROdata%DTocc%Range%Day(mn) .OR. &
       ROdata%DTocc%Day    > ROdata%DTocc%Range%Day(mx) )    &
       ROdata%DTocc%Day    = 99

  IF ( ROdata%DTocc%Hour   < ROdata%DTocc%Range%Hour(mn) .OR. &
       ROdata%DTocc%Hour   > ROdata%DTocc%Range%Hour(mx) )    &
       ROdata%DTocc%Hour   = 99

  IF ( ROdata%DTocc%Minute < ROdata%DTocc%Range%Minute(mn) .OR. &
       ROdata%DTocc%Minute > ROdata%DTocc%Range%Minute(mx) )    &
       ROdata%DTocc%Minute = 99

  IF ( ROdata%DTocc%Second < ROdata%DTocc%Range%Second(mn) .OR. &
       ROdata%DTocc%Second > ROdata%DTocc%Range%Second(mx) )    &
       ROdata%DTocc%Second = 99

  IF ( ROdata%DTocc%Msec   < ROdata%DTocc%Range%Msec(mn) .OR. &
       ROdata%DTocc%Msec   > ROdata%DTocc%Range%Msec(mx) )    &
       ROdata%DTocc%Msec   = 9999

! 1.3 Processing date/time

  IF ( ROdata%DTpro%Year   < ROdata%DTpro%Range%Year(mn) .OR. &
       ROdata%DTpro%Year   > ROdata%DTpro%Range%Year(mx) )    &
       ROdata%DTpro%Year   = 9999

  IF ( ROdata%DTpro%Month  < ROdata%DTpro%Range%Month(mn) .OR. &
       ROdata%DTpro%Month  > ROdata%DTpro%Range%Month(mx) )    &
       ROdata%DTpro%Month  = 99

  IF ( ROdata%DTpro%Day    < ROdata%DTpro%Range%Day(mn) .OR. &
       ROdata%DTpro%Day    > ROdata%DTpro%Range%Day(mx) )    &
       ROdata%DTpro%Day    = 99

  IF ( ROdata%DTpro%Hour   < ROdata%DTpro%Range%Hour(mn) .OR. &
       ROdata%DTpro%Hour   > ROdata%DTpro%Range%Hour(mx) )    &
       ROdata%DTpro%Hour   = 99

  IF ( ROdata%DTpro%Minute < ROdata%DTpro%Range%Minute(mn) .OR. &
       ROdata%DTpro%Minute > ROdata%DTpro%Range%Minute(mx) )    &
       ROdata%DTpro%Minute = 99

  IF ( ROdata%DTpro%Second < ROdata%DTpro%Range%Second(mn) .OR. &
       ROdata%DTpro%Second > ROdata%DTpro%Range%Second(mx) )    &
       ROdata%DTpro%Second = 99

  IF ( ROdata%DTpro%Msec   < ROdata%DTpro%Range%Msec(mn) .OR. &
       ROdata%DTpro%Msec   > ROdata%DTpro%Range%Msec(mx) )    &
       ROdata%DTpro%Msec   = 999

! 1.4 Quality parameters

  IF ( ROdata%PCD < ROdata%Range%PCD(mn) .OR. &
       ROdata%PCD > ROdata%Range%PCD(mx) .OR. &
       BTEST ( ROdata%PCD, PCD_missing ) ) THEN
       ROdata%PCD = 65535
  ELSE IF ( BTEST ( ROdata%PCD, PCD_phase    ) .OR. &
            BTEST ( ROdata%PCD, PCD_bangle   ) .OR. &
            BTEST ( ROdata%PCD, PCD_refrac   ) .OR. &
            BTEST ( ROdata%PCD, PCD_met      ) .OR. &
            BTEST ( ROdata%PCD, PCD_bg       ) .OR. &
            BTEST ( ROdata%PCD, PCD_missing  ) ) THEN
       ROdata%PCD = IBSET ( ROdata%PCD, PCD_summary )
  END IF

  IF ( ROdata%Overall_Qual < ROdata%Range%Overall_Qual(mn) .OR. &
       ROdata%Overall_Qual > ROdata%Range%Overall_Qual(mx) )    &
       ROdata%Overall_Qual = ropp_MDFV

! 1.5 Re-generate the Occultation ID string

  CALL ropp_io_occid ( ROdata )

! 1.6 Geo-referencing parameters

  IF ( ROdata%GeoRef%Time_Offset < ROdata%GeoRef%Range%Time_Offset(mn) .OR. &
       ROdata%GeoRef%Time_Offset > ROdata%GeoRef%Range%Time_Offset(mx) )    &
       ROdata%GeoRef%Time_Offset = ropp_MDFV

  IF ( ROdata%GeoRef%Lat         < ROdata%GeoRef%Range%Lat(mn) .OR. &
       ROdata%GeoRef%Lat         > ROdata%GeoRef%Range%Lat(mx) )    &
       ROdata%GeoRef%Lat         = ropp_MDFV

  IF ( ROdata%GeoRef%Lon         > 180.0_wp .AND. &  ! throw 180-360 --> -180-0
       ROdata%GeoRef%Lon        <= 360.0_wp )     &
       ROdata%GeoRef%Lon         = ROdata%GeoRef%Lon - 360.0_wp

  IF ( ROdata%GeoRef%Lon         < ROdata%GeoRef%Range%Lon(mn) .OR. &
       ROdata%GeoRef%Lon         > ROdata%GeoRef%Range%Lon(mx) )    &
       ROdata%GeoRef%Lon         = ropp_MDFV

  IF ( ROdata%GeoRef%RoC         < ROdata%GeoRef%Range%RoC(mn) .OR. &
       ROdata%GeoRef%RoC         > ROdata%GeoRef%Range%RoC(mx) )    &
       ROdata%GeoRef%RoC         = ropp_MDFV

  WHERE ( ROdata%GeoRef%R_CoC    < ROdata%GeoRef%Range%R_CoC(mn) .OR. &
          ROdata%GeoRef%R_CoC    > ROdata%GeoRef%Range%R_CoC(mx) )    &
          ROdata%GeoRef%R_CoC    = ropp_MDFV

  IF ( ABS(ROdata%GeoRef%R_CoC(1)-ropp_MDTV) < 1.0 .AND. &   ! all are -9999
       ABS(ROdata%GeoRef%R_CoC(2)-ropp_MDTV) < 1.0 .AND. &
       ABS(ROdata%GeoRef%R_CoC(3)-ropp_MDTV) < 1.0 )     &
           ROdata%GeoRef%R_CoC   = ropp_MDFV

  IF ( ABS(ROdata%GeoRef%R_CoC(1)-ropp_MDFV) < 1.0 .AND. &   ! all are -9999999
       ABS(ROdata%GeoRef%R_CoC(2)-ropp_MDFV) < 1.0 .AND. &
       ABS(ROdata%GeoRef%R_CoC(3)-ropp_MDFV) < 1.0 )     &
           ROdata%GeoRef%R_CoC   = ropp_MDFV

  IF ( ROdata%GeoRef%Azimuth     < ROdata%GeoRef%Range%Azimuth(mn) .OR. &
       ROdata%GeoRef%Azimuth     > ROdata%GeoRef%Range%Azimuth(mx) )    &
       ROdata%GeoRef%Azimuth     = ropp_MDFV

  IF ( ROdata%GeoRef%Undulation  < ROdata%GeoRef%Range%Undulation(mn) .OR. &
       ROdata%GeoRef%Undulation  > ROdata%GeoRef%Range%Undulation(mx) )    &
       ROdata%GeoRef%Undulation  = ropp_MDFV

! 1.7 Meteo background parameters

  CALL To_Upper ( ROdata%BG%Source )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%BG%Source(i:i)) == 0 ) &
      ROdata%BG%Source(i:i) = "_"
  END DO
  IF ( ROdata%BG%Source(1:4) == "____" ) ROdata%BG%Source = "UNKNOWN"

  IF ( ROdata%STN_ID == "____" ) ROdata%STN_ID = "UNKN"

  IF ( ROdata%BG%Year     < ROdata%BG%Range%Year(mn) .OR. &
       ROdata%BG%Year     > ROdata%BG%Range%Year(mx) )    &
       ROdata%BG%Year     = 9999

  IF ( ROdata%BG%Month    < ROdata%BG%Range%Month(mn) .OR. &
       ROdata%BG%Month    > ROdata%BG%Range%Month(mx) )    &
       ROdata%BG%Month    = 99

  IF ( ROdata%BG%Day      < ROdata%BG%Range%Day(mn) .OR. &
       ROdata%BG%Day      > ROdata%BG%Range%Day(mx) )    &
       ROdata%BG%Day      = 99

  IF ( ROdata%BG%Hour     < ROdata%BG%Range%Hour(mn) .OR. &
       ROdata%BG%Hour     > ROdata%BG%Range%Hour(mx) )    &
       ROdata%BG%Hour     = 99

  IF ( ROdata%BG%Minute   < ROdata%BG%Range%Minute(mn) .OR. &
       ROdata%BG%Minute   > ROdata%BG%Range%Minute(mx) )    &
       ROdata%BG%Minute   = 99

  IF ( ROdata%BG%FCperiod < ROdata%BG%Range%FCperiod(mn) .OR. &
       ROdata%BG%FCperiod > ROdata%BG%Range%FCperiod(mx) )    &
       ROdata%BG%FCperiod = 999.9

! 1.8 Level 1a profile

  CALL ropp_io_rangecheck (ROdata%Lev1a)

! 1.9 Level 1b profile

  CALL ropp_io_rangecheck (ROdata%Lev1b)

! 1.10 Level 2a profile

  CALL ropp_io_rangecheck (ROdata%Lev2a, ROdata%GeoRef%Lat)

! 1.11 Level 2b profile

  CALL ropp_io_rangecheck (ROdata%Lev2b)

! 1.12 Level 2c profile

  CALL ropp_io_rangecheck (ROdata%Lev2c)

! 1.13 Level 2d profile

  CALL ropp_io_rangecheck (ROdata%Lev2d)


END SUBROUTINE ropp_io_rangecheck_1d

!----------------------------------------------------------------
! 2. 2d ROprof structure
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_2d ( ROdata ) ! (inout)

! Modules

  USE typesizes,     ONLY: wp => EightByteReal
  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_2d
  USE ropp_utils,    ONLY : To_Upper, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

  CHARACTER (LEN =  *), PARAMETER :: validchars = &
                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_- " ! Valid chr set

! Argument list parameters

  TYPE(ROprof2d), INTENT(INOUT) :: ROdata    ! RO data structure

! Local variables

  INTEGER  :: i    ! Loop indices
! INTEGER  :: ikl  ! Keyword length    ! Commented at 20 July, 2016
! INTEGER  :: ict  ! Invalid chr count ! Commented at 20 July, 2016
  INTEGER  :: num  ! Temporary number
  INTEGER  :: ierr ! I/O error
  REAL     :: tmp  ! Temporary number

! 2.1 Character variables; ensure left-justified and uppercase letters;
!     replace 'illegal' characters with underscores, or check numeric
!     representations are valid. If the whole field was invalid, set
!     as unknown

! The following strings are fixed 4-character keywords with
! restricted character set [A-Z],[0-9],'_','-' only.

  ROdata%LEO_ID = ADJUSTL(ROdata%LEO_ID)
  CALL To_Upper ( ROdata%LEO_ID )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%LEO_ID(i:i)) == 0 ) &
      ROdata%LEO_ID(i:i) = "_"
  END DO
  IF ( INDEX(ROdata%LEO_ID(1:3), "_") > 0 )  ROdata%LEO_ID = "UNKN"

  ROdata%GNS_ID = ADJUSTL(ROdata%GNS_ID)
  CALL To_Upper ( ROdata%GNS_ID(1:1) )
  IF ( INDEX(validchars(1:26),ROdata%GNS_ID(1:1)) == 0 ) &
    ROdata%GNS_ID(1:1) = "U"
  READ  ( ROdata%GNS_ID(2:4), FMT="(I3)", IOSTAT=ierr ) num
  IF ( ierr /= 0 .OR. &
       num  < 001 .OR. num > 999 ) num = 999
  WRITE ( ROdata%GNS_ID(2:4), FMT="(I3.3)" ) num

  ROdata%STN_ID = ADJUSTL(ROdata%STN_ID)
  CALL To_Upper ( ROdata%STN_ID )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%STN_ID(i:i)) == 0 ) &
      ROdata%STN_ID(i:i) = "_"
  END DO
  IF ( INDEX(ROdata%STN_ID(1:3), "_") > 0 )  ROdata%STN_ID = "UNKN"

  CALL keyword_check(ROdata%Processing_centre)

  CALL keyword_check(ROdata%Processing_software)

  CALL keyword_check(ROdata%POD_Method)

  CALL keyword_check(ROdata%Phase_Method)

  CALL keyword_check(ROdata%Bangle_Method)

  CALL keyword_check(ROdata%Refrac_Method)

  CALL keyword_check(ROdata%Meteo_Method)

  CALL keyword_check(ROdata%Thin_Method)

! Ensure software version string conforms to "Vnn.nnn" format

  ROdata%Software_Version = ADJUSTL(ROdata%Software_Version)
  CALL To_Upper ( ROdata%Software_Version )
  IF ( ROdata%Software_Version(2:2) == " " ) &
       ROdata%Software_Version(2:2) =  "0"
  ROdata%Software_Version(8:) = " "
  READ ( ROdata%Software_Version(2:7), FMT=*, IOSTAT=ierr ) tmp
  IF ( ROdata%Software_Version(1:1) /= "V" .OR. &
       ierr                         /=  0  .OR. &
       tmp < 0.0 .OR. tmp > 99.999 ) ROdata%Software_Version = "UNKNOWN"

! 2.2 Occultation date/time

  IF ( ROdata%DTocc%Year   < ROdata%DTocc%Range%Year(mn) .OR. &
       ROdata%DTocc%Year   > ROdata%DTocc%Range%Year(mx) )    &
       ROdata%DTocc%Year   = 9999

  IF ( ROdata%DTocc%Month  < ROdata%DTocc%Range%Month(mn) .OR. &
       ROdata%DTocc%Month  > ROdata%DTocc%Range%Month(mx) )    &
       ROdata%DTocc%Month  = 99

  IF ( ROdata%DTocc%Day    < ROdata%DTocc%Range%Day(mn) .OR. &
       ROdata%DTocc%Day    > ROdata%DTocc%Range%Day(mx) )    &
       ROdata%DTocc%Day    = 99

  IF ( ROdata%DTocc%Hour   < ROdata%DTocc%Range%Hour(mn) .OR. &
       ROdata%DTocc%Hour   > ROdata%DTocc%Range%Hour(mx) )    &
       ROdata%DTocc%Hour   = 99

  IF ( ROdata%DTocc%Minute < ROdata%DTocc%Range%Minute(mn) .OR. &
       ROdata%DTocc%Minute > ROdata%DTocc%Range%Minute(mx) )    &
       ROdata%DTocc%Minute = 99

  IF ( ROdata%DTocc%Second < ROdata%DTocc%Range%Second(mn) .OR. &
       ROdata%DTocc%Second > ROdata%DTocc%Range%Second(mx) )    &
       ROdata%DTocc%Second = 99

  IF ( ROdata%DTocc%Msec   < ROdata%DTocc%Range%Msec(mn) .OR. &
       ROdata%DTocc%Msec   > ROdata%DTocc%Range%Msec(mx) )    &
       ROdata%DTocc%Msec   = 999

! 2.3 Processing date/time

  IF ( ROdata%DTpro%Year   < ROdata%DTpro%Range%Year(mn) .OR. &
       ROdata%DTpro%Year   > ROdata%DTpro%Range%Year(mx) )    &
       ROdata%DTpro%Year   = 9999

  IF ( ROdata%DTpro%Month  < ROdata%DTpro%Range%Month(mn) .OR. &
       ROdata%DTpro%Month  > ROdata%DTpro%Range%Month(mx) )    &
       ROdata%DTpro%Month  = 99

  IF ( ROdata%DTpro%Day    < ROdata%DTpro%Range%Day(mn) .OR. &
       ROdata%DTpro%Day    > ROdata%DTpro%Range%Day(mx) )    &
       ROdata%DTpro%Day    = 99

  IF ( ROdata%DTpro%Hour   < ROdata%DTpro%Range%Hour(mn) .OR. &
       ROdata%DTpro%Hour   > ROdata%DTpro%Range%Hour(mx) )    &
       ROdata%DTpro%Hour   = 99

  IF ( ROdata%DTpro%Minute < ROdata%DTpro%Range%Minute(mn) .OR. &
       ROdata%DTpro%Minute > ROdata%DTpro%Range%Minute(mx) )    &
       ROdata%DTpro%Minute = 99

  IF ( ROdata%DTpro%Second < ROdata%DTpro%Range%Second(mn) .OR. &
       ROdata%DTpro%Second > ROdata%DTpro%Range%Second(mx) )    &
       ROdata%DTpro%Second = 99

  IF ( ROdata%DTpro%Msec   < ROdata%DTpro%Range%Msec(mn) .OR. &
       ROdata%DTpro%Msec   > ROdata%DTpro%Range%Msec(mx) )    &
       ROdata%DTpro%Msec   = 999

! 2.4 Quality parameters

  IF ( ROdata%PCD < ROdata%Range%PCD(mn) .OR. &
       ROdata%PCD > ROdata%Range%PCD(mx) .OR. &
       BTEST ( ROdata%PCD, PCD_missing ) ) THEN
       ROdata%PCD = 65535
  ELSE IF ( BTEST ( ROdata%PCD, PCD_phase    ) .OR. &
            BTEST ( ROdata%PCD, PCD_bangle   ) .OR. &
            BTEST ( ROdata%PCD, PCD_refrac   ) .OR. &
            BTEST ( ROdata%PCD, PCD_met      ) .OR. &
            BTEST ( ROdata%PCD, PCD_bg       ) .OR. &
            BTEST ( ROdata%PCD, PCD_missing  ) ) THEN
       ROdata%PCD = IBSET ( ROdata%PCD, PCD_summary )
  END IF

  IF ( ROdata%Overall_Qual < ROdata%Range%Overall_Qual(mn) .OR. &
       ROdata%Overall_Qual > ROdata%Range%Overall_Qual(mx) )    &
       ROdata%Overall_Qual = ropp_MDFV

! 2.5 Re-generate the Occultation ID string

  CALL ropp_io_occid ( ROdata )

! 2.6 Geo-referencing parameters

  IF ( ROdata%GeoRef%Time_Offset < ROdata%GeoRef%Range%Time_Offset(mn) .OR. &
       ROdata%GeoRef%Time_Offset > ROdata%GeoRef%Range%Time_Offset(mx) )    &
       ROdata%GeoRef%Time_Offset = ropp_MDFV

  IF ( ROdata%GeoRef%Lat         < ROdata%GeoRef%Range%Lat(mn) .OR. &
       ROdata%GeoRef%Lat         > ROdata%GeoRef%Range%Lat(mx) )    &
       ROdata%GeoRef%Lat         = ropp_MDFV

  IF ( ROdata%GeoRef%Lon         > 180.0_wp .AND. &  ! throw 180-360 --> -180-0
       ROdata%GeoRef%Lon        <= 360.0_wp )     &
       ROdata%GeoRef%Lon         = ROdata%GeoRef%Lon - 360.0_wp

  IF ( ROdata%GeoRef%Lon         < ROdata%GeoRef%Range%Lon(mn) .OR. &
       ROdata%GeoRef%Lon         > ROdata%GeoRef%Range%Lon(mx) )    &
       ROdata%GeoRef%Lon         = ropp_MDFV

  IF ( ROdata%GeoRef%RoC         < ROdata%GeoRef%Range%RoC(mn) .OR. &
       ROdata%GeoRef%RoC         > ROdata%GeoRef%Range%RoC(mx) )    &
       ROdata%GeoRef%RoC         = ropp_MDFV

  WHERE ( ROdata%GeoRef%R_CoC    < ROdata%GeoRef%Range%R_CoC(mn) .OR. &
          ROdata%GeoRef%R_CoC    > ROdata%GeoRef%Range%R_CoC(mx) )    &
          ROdata%GeoRef%R_CoC    = ropp_MDFV

  IF ( ABS(ROdata%GeoRef%R_CoC(1)-ropp_MDTV) < 1.0 .AND. &   ! all are -9999
       ABS(ROdata%GeoRef%R_CoC(2)-ropp_MDTV) < 1.0 .AND. &
       ABS(ROdata%GeoRef%R_CoC(3)-ropp_MDTV) < 1.0 )     &
           ROdata%GeoRef%R_CoC   = ropp_MDFV

  IF ( ABS(ROdata%GeoRef%R_CoC(1)-ropp_MDFV) < 1.0 .AND. &   ! all are -9999999
       ABS(ROdata%GeoRef%R_CoC(2)-ropp_MDFV) < 1.0 .AND. &
       ABS(ROdata%GeoRef%R_CoC(3)-ropp_MDFV) < 1.0 )     &
           ROdata%GeoRef%R_CoC   = ropp_MDFV

  IF ( ROdata%GeoRef%Azimuth     < ROdata%GeoRef%Range%Azimuth(mn) .OR. &
       ROdata%GeoRef%Azimuth     > ROdata%GeoRef%Range%Azimuth(mx) )    &
       ROdata%GeoRef%Azimuth     = ropp_MDFV

  IF ( ROdata%GeoRef%Undulation  < ROdata%GeoRef%Range%Undulation(mn) .OR. &
       ROdata%GeoRef%Undulation  > ROdata%GeoRef%Range%Undulation(mx) )    &
       ROdata%GeoRef%Undulation  = ropp_MDFV

! 2.7 Meteo background parameters

  CALL To_Upper ( ROdata%BG%Source )
  DO i = 1, 4
    IF ( INDEX(validchars,ROdata%BG%Source(i:i)) == 0 ) &
      ROdata%BG%Source(i:i) = "_"
  END DO
  IF ( ROdata%BG%Source(1:4) == "____" ) ROdata%BG%Source = "UNKNOWN"

  IF ( ROdata%STN_ID == "____" ) ROdata%STN_ID = "UNKN"

  IF ( ROdata%BG%Year     < ROdata%BG%Range%Year(mn) .OR. &
       ROdata%BG%Year     > ROdata%BG%Range%Year(mx) )    &
       ROdata%BG%Year     = 9999

  IF ( ROdata%BG%Month    < ROdata%BG%Range%Month(mn) .OR. &
       ROdata%BG%Month    > ROdata%BG%Range%Month(mx) )    &
       ROdata%BG%Month    = 99

  IF ( ROdata%BG%Day      < ROdata%BG%Range%Day(mn) .OR. &
       ROdata%BG%Day      > ROdata%BG%Range%Day(mx) )    &
       ROdata%BG%Day      = 99

  IF ( ROdata%BG%Hour     < ROdata%BG%Range%Hour(mn) .OR. &
       ROdata%BG%Hour     > ROdata%BG%Range%Hour(mx) )    &
       ROdata%BG%Hour     = 99

  IF ( ROdata%BG%Minute   < ROdata%BG%Range%Minute(mn) .OR. &
       ROdata%BG%Minute   > ROdata%BG%Range%Minute(mx) )    &
       ROdata%BG%Minute   = 99

  IF ( ROdata%BG%FCperiod < ROdata%BG%Range%FCperiod(mn) .OR. &
       ROdata%BG%FCperiod > ROdata%BG%Range%FCperiod(mx) )    &
       ROdata%BG%FCperiod = 999.9

! 2.8 Level 1a profile

  CALL ropp_io_rangecheck (ROdata%Lev1a)

! 2.9 Level 1b profile

  CALL ropp_io_rangecheck (ROdata%Lev1b)

! 2.10 Level 2a profile

  CALL ropp_io_rangecheck (ROdata%Lev2a, ROdata%GeoRef%Lat)

! 2.11 Level 2b profile

  CALL ropp_io_rangecheck (ROdata%Lev2b)

! 2.12 Level 2c profile

  CALL ropp_io_rangecheck (ROdata%Lev2c)

! 2.13 Level 2d profile

  CALL ropp_io_rangecheck (ROdata%Lev2d)


END SUBROUTINE ropp_io_rangecheck_2d


!----------------------------------------------------------------
! 3. Level 1a profile parameters
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l1atype ( Lev1a ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l1atype
  USE ropp_utils,    ONLY : WHERE, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L1atype), INTENT(INOUT) :: Lev1a    ! RO data structure

! Local variables

  INTEGER  :: i, j ! Loop indices

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0
  INTEGER                        :: midx = 0

  IF ( Lev1a%Npoints <= 0 )THEN
       Lev1a%Npoints =  0

  ELSE

    WHERE ( Lev1a%Dtime      < Lev1a%Range%Dtime(mn) .OR.      &
            Lev1a%Dtime      > Lev1a%Range%Dtime(mx) )         &
            Lev1a%Dtime      = ropp_MDFV

    WHERE ( Lev1a%SNR_L1ca   < Lev1a%Range%SNR(mn) .OR.        &
            Lev1a%SNR_L1ca   > Lev1a%Range%SNR(mx) )           &
            Lev1a%SNR_L1ca   = ropp_MDFV

    WHERE ( Lev1a%SNR_L1p    < Lev1a%Range%SNR(mn) .OR.        &
            Lev1a%SNR_L1p    > Lev1a%Range%SNR(mx) )           &
            Lev1a%SNR_L1p    = ropp_MDFV

    WHERE ( Lev1a%SNR_L2p    < Lev1a%Range%SNR(mn) .OR.        &
            Lev1a%SNR_L2p    > Lev1a%Range%SNR(mx) )           &
            Lev1a%SNR_L2p    = ropp_MDFV

    WHERE ( Lev1a%Phase_L1   < Lev1a%Range%Phase(mn) .OR.      &
            Lev1a%Phase_L1   > Lev1a%Range%Phase(mx) )         &
            Lev1a%Phase_L1   = ropp_MDFV

    WHERE ( Lev1a%Phase_L2   < Lev1a%Range%Phase(mn) .OR.      &
            Lev1a%Phase_L2   > Lev1a%Range%Phase(mx) )         &
            Lev1a%Phase_L2   = ropp_MDFV

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%r_GNS(i,:)**2) > Lev1a%Range%r_GNS(mx)**2 ) &
           Lev1a%r_GNS(i,:) = ropp_MDFV
    END DO

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%r_GNS(i,:)**2) < 6.2e6_wp**2 ) & !6.2e6 = GeoRange%roc(mn)
           Lev1a%r_GNS(i,:) = ropp_MDFV
    END DO

    DO j = 1, 3
      DO i = 1, Lev1a%Npoints
        IF ( Lev1a%r_GNS(i,j) < Lev1a%Range%r_GNS(mn) .OR.    &
             Lev1a%r_GNS(i,j) > Lev1a%Range%r_GNS(mx) )       &
             Lev1a%r_GNS(i,j) = ropp_MDFV
      END DO
    END DO

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%v_GNS(i,:)**2) > Lev1a%Range%v_GNS(mx)**2 ) &
           Lev1a%v_GNS(i,:) = ropp_MDFV
    END DO

    DO j = 1, 3
      DO i = 1, Lev1a%Npoints
        IF ( Lev1a%v_GNS(i,j) < Lev1a%Range%v_GNS(mn) .OR.    &
             Lev1a%v_GNS(i,j) > Lev1a%Range%v_GNS(mx) )       &
             Lev1a%v_GNS(i,j) = ropp_MDFV
      END DO
    END DO

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%r_LEO(i,:)**2) > Lev1a%Range%r_LEO(mx)**2 ) &
           Lev1a%r_LEO(i,:) = ropp_MDFV
    END DO

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%r_LEO(i,:)**2) < 6.2e6_wp**2 ) & !6.2e6 = GeoRange%roc(mn)
           Lev1a%r_LEO(i,:) = ropp_MDFV
    END DO

    DO j = 1, 3
      DO i = 1, Lev1a%Npoints
        IF ( Lev1a%r_LEO(i,j) < Lev1a%Range%r_LEO(mn) .OR.    &
             Lev1a%r_LEO(i,j) > Lev1a%Range%r_LEO(mx) )       &
             Lev1a%r_LEO(i,j) = ropp_MDFV
      END DO
    END DO

    DO i = 1, Lev1a%Npoints
      IF ( SUM(Lev1a%v_LEO(i,:)**2) > Lev1a%Range%v_LEO(mx)**2 ) &
           Lev1a%v_LEO(i,:) = ropp_MDFV
    END DO

    DO j = 1, 3
      DO i = 1, Lev1a%Npoints
        IF ( Lev1a%v_LEO(i,j) < Lev1a%Range%v_LEO(mn) .OR.    &
             Lev1a%v_LEO(i,j) > Lev1a%Range%v_LEO(mx) )       &
             Lev1a%v_LEO(i,j) = ropp_MDFV
      END DO
    END DO

    WHERE ( Lev1a%Phase_Qual < Lev1a%Range%Phase_Qual(mn) .OR. &
            Lev1a%Phase_Qual > Lev1a%Range%Phase_Qual(mx) )    &
            Lev1a%Phase_Qual = ropp_MDFV

! All samples must have a valid time offset value.
! If all are missing, then set whole L1a profile to zero length.
! If only some times are missing, filter them out

    idx => WHERE ( Lev1a%DTime > ropp_MDTV, nidx )

    IF ( nidx == 0 ) THEN
      Lev1a%Npoints = 0
    ELSE IF ( nidx < Lev1a%Npoints ) THEN
      Lev1a%Npoints              = nidx
      Lev1a%DTime(1:nidx)        = Lev1a%DTime(idx)
      Lev1a%SNR_L1ca(1:nidx)     = Lev1a%SNR_L1ca(idx)
      Lev1a%SNR_L1p(1:nidx)      = Lev1a%SNR_L1p(idx)
      Lev1a%SNR_L2p(1:nidx)      = Lev1a%SNR_L2p(idx)
      Lev1a%Phase_L1(1:nidx)     = Lev1a%Phase_L1(idx)
      Lev1a%Phase_L2(1:nidx)     = Lev1a%Phase_L2(idx)
      DO j = 1, 3
         Lev1a%r_GNS(1:nidx,j)   = Lev1a%r_GNS(idx,j)
         Lev1a%v_GNS(1:nidx,j)   = Lev1a%v_GNS(idx,j)
         Lev1a%r_LEO(1:nidx,j)   = Lev1a%r_LEO(idx,j)
         Lev1a%v_LEO(1:nidx,j)   = Lev1a%v_LEO(idx,j)
      END DO
      Lev1a%Phase_Qual(1:nidx)   = Lev1a%Phase_Qual(idx)
    END IF

! If there are no valid SNR or Phase or POD values of any type,
! set profile missing flag to indicate no valid Level 1a data

    IF ( Lev1a%Npoints > 0 ) THEN
      idx => WHERE ( Lev1a%SNR_L1ca   > ropp_MDTV .OR. &
                     Lev1a%SNR_L1p    > ropp_MDTV .OR. &
                     Lev1a%SNR_L2p    > ropp_MDTV .OR. &
                     Lev1a%Phase_L1   > ropp_MDTV .OR. &
                     Lev1a%Phase_L2   > ropp_MDTV, nidx )
      idx => WHERE ( Lev1a%r_LEO(:,1) > ropp_MDFV .OR. &
                     Lev1a%v_LEO(:,1) > ropp_MDFV .OR. &
                     Lev1a%r_GNS(:,1) > ropp_MDFV .OR. &
                     Lev1a%v_GNS(:,1) > ropp_MDFV, midx )
      IF ( nidx+midx == 0 ) Lev1a%Missing = .TRUE.
    END IF

  END IF

  IF ( Lev1a%Npoints == 0 ) Lev1a%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l1atype

!----------------------------------------------------------------
! 4. Level 1b profile parameters
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l1btype ( Lev1b ) ! (inout)

  USE typesizes,     ONLY: wp => EightByteReal
  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l1btype
  USE ropp_utils,    ONLY : WHERE, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L1btype), INTENT(INOUT) :: Lev1b    ! RO data structure

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0

  IF ( Lev1b%Npoints <= 0 ) THEN
       Lev1b%Npoints =  0

  ELSE

    WHERE ( Lev1b%Lat_tp           < Lev1b%Range%Lat_tp(mn) .OR. &
            Lev1b%Lat_tp           > Lev1b%Range%Lat_tp(mx) )    &
            Lev1b%Lat_tp           = ropp_MDFV

    WHERE ( Lev1b%Lon_tp           > 180.0_wp .AND. &  ! throw 180-360 --> -180-0
            Lev1b%Lon_tp          <= 360.0_wp )     &
            Lev1b%Lon_tp           = Lev1b%Lon_tp - 360.0_wp

    WHERE ( Lev1b%Lon_tp           < Lev1b%Range%Lon_tp(mn) .OR. &
            Lev1b%Lon_tp           > Lev1b%Range%Lon_tp(mx) )    &
            Lev1b%Lon_tp           = ropp_MDFV

    WHERE ( Lev1b%Azimuth_tp       < Lev1b%Range%Azimuth_tp(mn) .OR. &
            Lev1b%Azimuth_tp       > Lev1b%Range%Azimuth_tp(mx) )    &
            Lev1b%Azimuth_tp       = ropp_MDFV

    WHERE ( Lev1b%Impact_L1        < Lev1b%Range%Impact(mn) .OR. &
            Lev1b%Impact_L1        > Lev1b%Range%Impact(mx) )    &
            Lev1b%Impact_L1        = ropp_MDFV

    WHERE ( Lev1b%Impact_L2        < Lev1b%Range%Impact(mn) .OR. &
            Lev1b%Impact_L2        > Lev1b%Range%Impact(mx) )    &
            Lev1b%Impact_L2        = ropp_MDFV

    WHERE ( Lev1b%Impact           < Lev1b%Range%Impact(mn) .OR. &
            Lev1b%Impact           > Lev1b%Range%Impact(mx) )    &
            Lev1b%Impact           = ropp_MDFV

    WHERE ( Lev1b%Impact_Opt       < Lev1b%Range%Impact(mn) .OR. &
            Lev1b%Impact_Opt       > Lev1b%Range%Impact(mx) )    &
            Lev1b%Impact_Opt       = ropp_MDFV

    WHERE ( Lev1b%Bangle_L1        < Lev1b%Range%Bangle(mn) .OR. &
            Lev1b%Bangle_L1        > Lev1b%Range%Bangle(mx) )    &
            Lev1b%Bangle_L1        = ropp_MDFV

    WHERE ( Lev1b%Bangle_L2        < Lev1b%Range%Bangle(mn) .OR. &
            Lev1b%Bangle_L2        > Lev1b%Range%Bangle(mx) )    &
            Lev1b%Bangle_L2        = ropp_MDFV

    WHERE ( Lev1b%Bangle           < Lev1b%Range%Bangle(mn) .OR. &
            Lev1b%Bangle           > Lev1b%Range%Bangle(mx) )    &
           Lev1b%Bangle           = ropp_MDFV

    WHERE ( Lev1b%Bangle_Opt       < Lev1b%Range%Bangle(mn) .OR. &
            Lev1b%Bangle_Opt       > Lev1b%Range%Bangle(mx) )    &
            Lev1b%Bangle_Opt       = ropp_MDFV

    WHERE ( Lev1b%Bangle_L1_Sigma  < Lev1b%Range%Bangle_Sigma(mn) .OR. &
            Lev1b%Bangle_L1_Sigma  > Lev1b%Range%Bangle_Sigma(mx) )    &
            Lev1b%Bangle_L1_Sigma  = ropp_MDFV

    WHERE ( Lev1b%Bangle_L2_Sigma  < Lev1b%Range%Bangle_Sigma(mn) .OR. &
            Lev1b%Bangle_L2_Sigma  > Lev1b%Range%Bangle_Sigma(mx) )    &
            Lev1b%Bangle_L2_Sigma  = ropp_MDFV

    WHERE ( Lev1b%Bangle_Sigma     < Lev1b%Range%Bangle_Sigma(mn) .OR. &
            Lev1b%Bangle_Sigma     > Lev1b%Range%Bangle_Sigma(mx) )    &
            Lev1b%Bangle_Sigma     = ropp_MDFV

    WHERE ( Lev1b%Bangle_Opt_Sigma < Lev1b%Range%Bangle_Sigma(mn) .OR. &
            Lev1b%Bangle_Opt_Sigma > Lev1b%Range%Bangle_Sigma(mx) )    &
            Lev1b%Bangle_Opt_Sigma = ropp_MDFV

    WHERE ( Lev1b%Bangle_L1_Qual   < Lev1b%Range%Bangle_Qual(mn) .OR. &
            Lev1b%Bangle_L1_Qual   > Lev1b%Range%Bangle_Qual(mx) )    &
            Lev1b%Bangle_L1_Qual   = ropp_MDFV

    WHERE ( Lev1b%Bangle_L2_Qual   < Lev1b%Range%Bangle_Qual(mn) .OR. &
            Lev1b%Bangle_L2_Qual   > Lev1b%Range%Bangle_Qual(mx) )    &
            Lev1b%Bangle_L2_Qual   = ropp_MDFV

    WHERE ( Lev1b%Bangle_Qual      < Lev1b%Range%Bangle_Qual(mn) .OR. &
            Lev1b%Bangle_Qual      > Lev1b%Range%Bangle_Qual(mx) )    &
            Lev1b%Bangle_Qual      = ropp_MDFV

    WHERE ( Lev1b%Bangle_Opt_Qual  < Lev1b%Range%Bangle_Qual(mn) .OR. &
            Lev1b%Bangle_Opt_Qual  > Lev1b%Range%Bangle_Qual(mx) )    &
            Lev1b%Bangle_Opt_Qual  = ropp_MDFV

! Every level must have at least one valid impact parameter.
! If all four are missing on all levels, then set whole L1a profile to
! zero length. If only some impact height levels are missing, filter them out

    idx => WHERE ( Lev1b%Impact_L1  > ropp_MDTV .OR. &
                   Lev1b%Impact_L2  > ropp_MDTV .OR. &
                   Lev1b%Impact     > ropp_MDTV .OR. &
                   Lev1b%Impact_Opt > ropp_MDTV, nidx )

    IF ( nidx == 0 ) THEN
      Lev1b%Npoints = 0
    ELSE IF ( nidx < Lev1b%Npoints ) THEN
      Lev1b%Npoints                   = nidx
      Lev1b%Lat_tp(1:nidx)            = Lev1b%Lat_tp(idx)
      Lev1b%Lon_tp(1:nidx)            = Lev1b%Lon_tp(idx)
      Lev1b%Azimuth_tp(1:nidx)        = Lev1b%Azimuth_tp(idx)
      Lev1b%Impact_L1(1:nidx)         = Lev1b%Impact_L1(idx)
      Lev1b%Impact_L2(1:nidx)         = Lev1b%Impact_L2(idx)
      Lev1b%Impact(1:nidx)            = Lev1b%Impact(idx)
      Lev1b%Impact_Opt(1:nidx)        = Lev1b%Impact_Opt(idx)
      Lev1b%Bangle_L1(1:nidx)         = Lev1b%Bangle_L1(idx)
      Lev1b%Bangle_L2(1:nidx)         = Lev1b%Bangle_L2(idx)
      Lev1b%Bangle(1:nidx)            = Lev1b%Bangle(idx)
      Lev1b%Bangle_Opt(1:nidx)        = Lev1b%Bangle_Opt(idx)
      Lev1b%Bangle_L1_Sigma(1:nidx)   = Lev1b%Bangle_L1_Sigma(idx)
      Lev1b%Bangle_L2_Sigma(1:nidx)   = Lev1b%Bangle_L2_Sigma(idx)
      Lev1b%Bangle_Sigma(1:nidx)      = Lev1b%Bangle_Sigma(idx)
      Lev1b%Bangle_Opt_Sigma(1:nidx)  = Lev1b%Bangle_Opt_Sigma(idx)
      Lev1b%Bangle_L1_Qual(1:nidx)    = Lev1b%Bangle_L1_Qual(idx)
      Lev1b%Bangle_L2_Qual(1:nidx)    = Lev1b%Bangle_L2_Qual(idx)
      Lev1b%Bangle_Qual(1:nidx)       = Lev1b%Bangle_Qual(idx)
      Lev1b%Bangle_Opt_Qual(1:nidx)   = Lev1b%Bangle_Opt_Qual(idx)
    END IF

! If there are no valid bending angles of any type, set profile
! missing flag to indicate no valid Level 1b data

    IF ( Lev1b%Npoints > 0 ) THEN
      idx => where ( Lev1b%Bangle_L1  > ropp_MDTV .OR. &
                     Lev1b%Bangle_L2  > ropp_MDTV .OR. &
                     Lev1b%Bangle     > ropp_MDTV .OR. &
                     Lev1b%Bangle_Opt > ropp_MDTV, nidx )
      IF ( nidx == 0 ) Lev1b%Missing = .TRUE.
    END IF

  END IF

  IF ( Lev1b%Npoints == 0 ) Lev1b%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l1btype

!----------------------------------------------------------------
! 5. Level 2a profile parameters
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2atype ( Lev2a, Lat ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2atype
  USE ropp_utils,    ONLY : WHERE, ropp_MDFV, ropp_MDTV
  USE geodesy,       ONLY: geometric2geopotential

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L2atype), INTENT(INOUT) :: Lev2a    ! RO data structure
  REAL(wp),      INTENT(IN)    :: Lat      ! Occultation latitude

  INTEGER  :: n    ! No. of points in profile

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0
  REAL(wp)                       :: occlat

  IF ( Lev2a%Npoints <= 0 ) THEN
       Lev2a%Npoints =  0

  ELSE

    WHERE ( Lev2a%Alt_Refrac   < Lev2a%Range%Alt_Refrac(mn) .OR. &
            Lev2a%Alt_Refrac   > Lev2a%Range%Alt_Refrac(mx) )    &
            Lev2a%Alt_Refrac   = ropp_MDFV

    WHERE ( Lev2a%Geop_Refrac  < Lev2a%Range%Geop_Refrac(mn) .OR. &
            Lev2a%Geop_Refrac  > Lev2a%Range%Geop_Refrac(mx) )    &
            Lev2a%Geop_Refrac  = ropp_MDFV

    WHERE ( Lev2a%Refrac       < Lev2a%Range%Refrac(mn) .OR. &
            Lev2a%Refrac       > Lev2a%Range%Refrac(mx) )    &
            Lev2a%Refrac       = ropp_MDFV

    WHERE ( Lev2a%Refrac_Sigma < Lev2a%Range%Refrac_Sigma(mn) .OR. &
            Lev2a%Refrac_Sigma > Lev2a%Range%Refrac_Sigma(mx) )    &
            Lev2a%Refrac_Sigma = ropp_MDFV

    WHERE ( Lev2a%Refrac_Qual  < Lev2a%Range%Refrac_Qual(mn) .OR. &
            Lev2a%Refrac_Qual  > Lev2a%Range%Refrac_Qual(mx) )    &
            Lev2a%Refrac_Qual  = ropp_MDFV

    WHERE ( Lev2a%Dry_Temp       < Lev2a%Range%Dry_Temp(mn) .OR. &
            Lev2a%Dry_Temp       > Lev2a%Range%Dry_Temp(mx) )    &
            Lev2a%Dry_Temp       = ropp_MDFV

    WHERE ( Lev2a%Dry_Temp_Sigma < Lev2a%Range%Dry_Temp_Sigma(mn) .OR. &
            Lev2a%Dry_Temp_Sigma > Lev2a%Range%Dry_Temp_Sigma(mx) )    &
            Lev2a%Dry_Temp_Sigma = ropp_MDFV

    WHERE ( Lev2a%Dry_Temp_Qual  < Lev2a%Range%Dry_Temp_Qual(mn) .OR. &
            Lev2a%Dry_Temp_Qual  > Lev2a%Range%Dry_Temp_Qual(mx) )    &
            Lev2a%Dry_Temp_Qual  = ropp_MDFV

! Every level must have a valid geometric altitude value;
! If all are missing, then set whole L1a profile to zero length.
! If only some altitudes are missing, filter them out

    idx => WHERE ( Lev2a%Alt_Refrac > ropp_MDTV, nidx )
    IF ( nidx == 0 ) THEN
      Lev2a%Npoints = 0
    ELSE IF ( nidx < Lev2a%Npoints ) THEN
      Lev2a%Npoints                = nidx
      Lev2a%Alt_Refrac(1:nidx)     = Lev2a%Alt_Refrac(idx)
      Lev2a%Geop_Refrac(1:nidx)    = Lev2a%Geop_Refrac(idx)
      Lev2a%Refrac(1:nidx)         = Lev2a%Refrac(idx)
      Lev2a%Refrac_Sigma(1:nidx)   = Lev2a%Refrac_Sigma(idx)
      Lev2a%Refrac_Qual(1:nidx)    = Lev2a%Refrac_Qual(idx)
      Lev2a%Dry_Temp(1:nidx)       = Lev2a%Dry_Temp(idx)
      Lev2a%Dry_Temp_Sigma(1:nidx) = Lev2a%Dry_Temp_Sigma(idx)
      Lev2a%Dry_Temp_Qual(1:nidx)  = Lev2a%Dry_Temp_Qual(idx)
    END IF

! If any missing geopotential heights, generate whole profile
! from geometric heights for consistency (if occ. Latitude is
! invalid, assume 45deg)

    IF ( Lev2a%Npoints > 0 ) THEN
      n   = Lev2a%Npoints
      occlat = Lat
      IF ( occlat < ropp_MDTV ) occlat = 45.0_wp
      idx => WHERE ( Lev2a%Geop_Refrac(1:n) < ropp_MDTV, nidx )
      IF ( nidx > 0 ) THEN
      Lev2a%Geop_Refrac(idx) = &
             geometric2geopotential( occlat, Lev2a%Alt_Refrac(idx) )
      END IF
    END IF

! If there are no valid refractivities, set profile
! missing flag to indicate no valid Level 2a data

    IF ( Lev2a%Npoints > 0 ) THEN
      idx => where ( Lev2a%Refrac > ropp_MDTV, nidx )
      IF ( nidx == 0 ) Lev2a%Missing = .TRUE.
    END IF

  END IF

  IF ( Lev2a%Npoints == 0 ) Lev2a%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2atype

!----------------------------------------------------------------
! 6. Level 2b profile parameters (1d version)
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2btype ( Lev2b ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2btype
  USE ropp_utils,    ONLY : WHERE, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L2btype), INTENT(INOUT) :: Lev2b    ! RO data structure

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0

  IF ( Lev2b%Npoints <= 0 ) THEN
       Lev2b%Npoints =  0

  ELSE

    WHERE ( Lev2b%Geop        < Lev2b%Range%Geop(mn) .OR. &
            Lev2b%Geop        > Lev2b%Range%Geop(mx) )    &
            Lev2b%Geop        = ropp_MDFV

    WHERE ( Lev2b%Geop_Sigma  < Lev2b%Range%Geop_Sigma(mn) .OR. &
            Lev2b%Geop_Sigma  > Lev2b%Range%Geop_Sigma(mx) )    &
            Lev2b%Geop_Sigma  = ropp_MDFV

    WHERE ( Lev2b%Press       < Lev2b%Range%Press(mn) .OR. &
            Lev2b%Press       > Lev2b%Range%Press(mx) )    &
            Lev2b%Press       = ropp_MDFV

    WHERE ( Lev2b%Press_Sigma < Lev2b%Range%Press_Sigma(mn) .OR. &
            Lev2b%Press_Sigma > Lev2b%Range%Press_Sigma(mx) )    &
            Lev2b%Press_Sigma = ropp_MDFV

    WHERE ( Lev2b%Temp        < Lev2b%Range%Temp(mn) .OR. &
            Lev2b%Temp        > Lev2b%Range%Temp(mx) )    &
            Lev2b%Temp        = ropp_MDFV

    WHERE ( Lev2b%Temp_Sigma  < Lev2b%Range%Temp_Sigma(mn) .OR. &
            Lev2b%Temp_Sigma  > Lev2b%Range%Temp_Sigma(mx) )    &
            Lev2b%Temp_Sigma  = ropp_MDFV

    WHERE ( Lev2b%SHum        < Lev2b%Range%SHum(mn) .OR. &
            Lev2b%SHum        > Lev2b%Range%SHum(mx) )    &
            Lev2b%SHum        = ropp_MDFV

    WHERE ( Lev2b%SHum_Sigma  < Lev2b%Range%SHum_Sigma(mn) .OR. &
            Lev2b%SHum_Sigma  > Lev2b%Range%SHum_Sigma(mx) )    &
            Lev2b%SHum_Sigma  = ropp_MDFV

    WHERE ( Lev2b%Meteo_Qual  < Lev2b%Range%Meteo_Qual(mn) .OR. &
            Lev2b%Meteo_Qual  > Lev2b%Range%Meteo_Qual(mx) )    &
            Lev2b%Meteo_Qual  = ropp_MDFV

! Every level must have a valid geopotential height value;
! if all are missing, then set whole L2b profile to zero length.
! If only some heights are missing, filter them out

    idx => WHERE ( Lev2b%Geop > ropp_MDTV, nidx )
    IF ( nidx == 0 ) THEN
      Lev2b%Npoints = 0
    ELSE IF ( nidx < Lev2b%Npoints ) THEN
      Lev2b%Npoints             = nidx
      Lev2b%Geop(1:nidx)        = Lev2b%Geop(idx)
      Lev2b%Geop_Sigma(1:nidx)  = Lev2b%Geop_Sigma(idx)
      Lev2b%Press(1:nidx)       = Lev2b%Press(idx)
      Lev2b%Press_Sigma(1:nidx) = Lev2b%Press_Sigma(idx)
      Lev2b%Temp(1:nidx)        = Lev2b%Temp(idx)
      Lev2b%Temp_Sigma(1:nidx)  = Lev2b%Temp_Sigma(idx)
      Lev2b%SHum(1:nidx)        = Lev2b%SHum(idx)
      Lev2b%SHum_Sigma(1:nidx)  = Lev2b%SHum_Sigma(idx)
      Lev2b%Meteo_Qual(1:nidx)  = Lev2b%Meteo_Qual(idx)
    END IF

! If there are no valid T,q,P values, set profile missing flag
! to indicate no valid Level 2b data

    IF ( Lev2b%Npoints > 0 ) THEN
      idx => where ( Lev2b%Press > ropp_MDTV .OR. &
                     Lev2b%Temp  > ropp_MDTV .OR. &
                     Lev2b%SHum  > ropp_MDTV, nidx )
      IF ( nidx == 0 ) Lev2b%Missing = .TRUE.
    END IF

  END IF

  IF ( Lev2b%Npoints == 0 ) Lev2b%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2btype

!----------------------------------------------------------------
! 7. Level 2b profile parameters (2d version)
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2btype_2d ( Lev2b ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2btype_2d
  USE ropp_utils,    ONLY : WHERE, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L2btype_2d), INTENT(INOUT) :: Lev2b    ! RO data structure

! Local variables

  INTEGER  :: i ! Loop indices

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0

  IF ( Lev2b%Nhoriz <= 0 ) THEN
    Lev2b%Npoints =  0
    Lev2b%Nhoriz  =  0
  ENDIF

  IF ( Lev2b%Npoints <= 0 ) THEN
       Lev2b%Npoints =  0
       Lev2b%Nhoriz  =  0

  ELSE

    WHERE ( Lev2b%Geop        < Lev2b%Range%Geop(mn) .OR. &
            Lev2b%Geop        > Lev2b%Range%Geop(mx) )    &
            Lev2b%Geop        = ropp_MDFV

    WHERE ( Lev2b%Geop_Sigma  < Lev2b%Range%Geop_Sigma(mn) .OR. &
            Lev2b%Geop_Sigma  > Lev2b%Range%Geop_Sigma(mx) )    &
            Lev2b%Geop_Sigma  = ropp_MDFV

    WHERE ( Lev2b%Press       < Lev2b%Range%Press(mn) .OR. &
            Lev2b%Press       > Lev2b%Range%Press(mx) )    &
            Lev2b%Press       = ropp_MDFV

    WHERE ( Lev2b%Press_Sigma < Lev2b%Range%Press_Sigma(mn) .OR. &
            Lev2b%Press_Sigma > Lev2b%Range%Press_Sigma(mx) )    &
            Lev2b%Press_Sigma = ropp_MDFV

    WHERE ( Lev2b%Temp        < Lev2b%Range%Temp(mn) .OR. &
            Lev2b%Temp        > Lev2b%Range%Temp(mx) )    &
            Lev2b%Temp        = ropp_MDFV

    WHERE ( Lev2b%Temp_Sigma  < Lev2b%Range%Temp_Sigma(mn) .OR. &
            Lev2b%Temp_Sigma  > Lev2b%Range%Temp_Sigma(mx) )    &
            Lev2b%Temp_Sigma  = ropp_MDFV

    WHERE ( Lev2b%SHum        < Lev2b%Range%SHum(mn) .OR. &
            Lev2b%SHum        > Lev2b%Range%SHum(mx) )    &
            Lev2b%SHum        = ropp_MDFV

    WHERE ( Lev2b%SHum_Sigma  < Lev2b%Range%SHum_Sigma(mn) .OR. &
            Lev2b%SHum_Sigma  > Lev2b%Range%SHum_Sigma(mx) )    &
            Lev2b%SHum_Sigma  = ropp_MDFV

    WHERE ( Lev2b%Meteo_Qual  < Lev2b%Range%Meteo_Qual(mn) .OR. &
            Lev2b%Meteo_Qual  > Lev2b%Range%Meteo_Qual(mx) )    &
            Lev2b%Meteo_Qual  = ropp_MDFV

! Every level must have a valid geopotential height value;
! if all are missing, then set whole L2b profile to zero length.
! If only some heights are missing, filter them out

    DO i = 1, Lev2b%Nhoriz
      idx => WHERE ( Lev2b%Geop(:,i) > ropp_MDTV, nidx )
      IF ( nidx < Lev2b%Npoints ) THEN
        Lev2b%Npoints             = nidx
        Lev2b%Geop(1:nidx,i)        = Lev2b%Geop(idx,i)
        Lev2b%Geop_Sigma(1:nidx,i)  = Lev2b%Geop_Sigma(idx,i)
        Lev2b%Press(1:nidx,i)       = Lev2b%Press(idx,i)
        Lev2b%Press_Sigma(1:nidx,i) = Lev2b%Press_Sigma(idx,i)
        Lev2b%Temp(1:nidx,i)        = Lev2b%Temp(idx,i)
        Lev2b%Temp_Sigma(1:nidx,i)  = Lev2b%Temp_Sigma(idx,i)
        Lev2b%SHum(1:nidx,i)        = Lev2b%SHum(idx,i)
        Lev2b%SHum_Sigma(1:nidx,i)  = Lev2b%SHum_Sigma(idx,i)
        Lev2b%Meteo_Qual(1:nidx,i)  = Lev2b%Meteo_Qual(idx,i)
      END IF
    ENDDO

  END IF

  IF ( Lev2b%Npoints == 0 ) Lev2b%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2btype_2d

!----------------------------------------------------------------
! 8. Level 2c profile parameters (1d version)
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2ctype ( Lev2c ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2ctype
  USE ropp_utils,    ONLY : ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L2ctype), INTENT(INOUT) :: Lev2c    ! RO data structure

  Lev2c%Npoints = MIN ( MAX ( Lev2c%Npoints, 0 ), 1 )

  IF ( Lev2c%Npoints > 0 ) THEN

    IF ( Lev2c%Geop_Sfc          < Lev2c%Range%Geop_Sfc(mn) .OR. &
         Lev2c%Geop_Sfc          > Lev2c%Range%Geop_Sfc(mx) )    &
         Lev2c%Geop_Sfc          = ropp_MDFV

    IF ( Lev2c%Press_Sfc         < Lev2c%Range%Press_Sfc(mn) .OR. &
         Lev2c%Press_Sfc         > Lev2c%Range%Press_Sfc(mx) )    &
         Lev2c%Press_Sfc         = ropp_MDFV

    IF ( Lev2c%Press_Sfc_Sigma   < Lev2c%Range%Press_Sfc_Sigma(mn) .OR. &
         Lev2c%Press_Sfc_Sigma   > Lev2c%Range%Press_Sfc_Sigma(mx) )    &
         Lev2c%Press_Sfc_Sigma   = ropp_MDFV

    IF ( Lev2c%Press_Sfc_Qual    < Lev2c%Range%Press_Sfc_Qual(mn) .OR. &
         Lev2c%Press_Sfc_Qual    > Lev2c%Range%Press_Sfc_Qual(mx) )    &
         Lev2c%Press_Sfc_Qual    = ropp_MDFV

    IF ( Lev2c%Ne_max            < Lev2c%Range%Ne_max(mn) .OR. &
         Lev2c%Ne_max            > Lev2c%Range%Ne_max(mx) )    &
         Lev2c%Ne_max            = ropp_MDFV

    IF ( Lev2c%Ne_max_sigma      < Lev2c%Range%Ne_max_sigma(mn) .OR. &
         Lev2c%Ne_max_sigma      > Lev2c%Range%Ne_max_sigma(mx) )    &
         Lev2c%Ne_max_sigma      = ropp_MDFV

    IF ( Lev2c%H_peak            < Lev2c%Range%H_peak(mn) .OR. &
         Lev2c%H_peak            > Lev2c%Range%H_peak(mx) )    &
         Lev2c%H_peak            = ropp_MDFV

    IF ( Lev2c%H_peak_sigma      < Lev2c%Range%H_peak_sigma(mn) .OR. &
         Lev2c%H_peak_sigma      > Lev2c%Range%H_peak_sigma(mx) )    &
         Lev2c%H_peak_sigma      = ropp_MDFV

    IF ( Lev2c%H_width           < Lev2c%Range%H_width(mn) .OR. &
         Lev2c%H_width           > Lev2c%Range%H_width(mx) )    &
         Lev2c%H_width           = ropp_MDFV

    IF ( Lev2c%H_width_sigma     < Lev2c%Range%H_width_sigma(mn) .OR. &
         Lev2c%H_width_sigma     > Lev2c%Range%H_width_sigma(mx) )    &
         Lev2c%H_width_sigma     = ropp_MDFV

    IF ( Lev2c%TPH_Bangle        < Lev2c%Range%TPH_Bangle(mn) .OR. &
         Lev2c%TPH_Bangle        > Lev2c%Range%TPH_Bangle(mx) )    &
         Lev2c%TPH_Bangle        = ropp_MDFV

    IF ( Lev2c%TPA_Bangle        < Lev2c%Range%TPA_Bangle(mn) .OR. &
         Lev2c%TPA_Bangle        > Lev2c%Range%TPA_Bangle(mx) )    &
         Lev2c%TPA_Bangle        = ropp_MDFV

    IF ( Lev2c%TPH_Bangle_Flag   < Lev2c%Range%TPH_Bangle_Flag(mn) .OR. &
         Lev2c%TPH_Bangle_Flag   > Lev2c%Range%TPH_Bangle_Flag(mx) )    &
         Lev2c%TPH_Bangle_Flag   = ropp_MIFV

    IF ( Lev2c%TPH_Refrac        < Lev2c%Range%TPH_Refrac(mn) .OR. &
         Lev2c%TPH_Refrac        > Lev2c%Range%TPH_Refrac(mx) )    &
         Lev2c%TPH_Refrac        = ropp_MDFV

    IF ( Lev2c%TPN_Refrac        < Lev2c%Range%TPN_Refrac(mn) .OR. &
         Lev2c%TPN_Refrac        > Lev2c%Range%TPN_Refrac(mx) )    &
         Lev2c%TPN_Refrac        = ropp_MDFV

    IF ( Lev2c%TPH_Refrac_Flag   < Lev2c%Range%TPH_Refrac_Flag(mn) .OR. &
         Lev2c%TPH_Refrac_Flag   > Lev2c%Range%TPH_Refrac_Flag(mx) )    &
         Lev2c%TPH_Refrac_Flag   = ropp_MIFV

    IF ( Lev2c%TPH_Tdry_LRT      < Lev2c%Range%TPH_Tdry_LRT(mn) .OR. &
         Lev2c%TPH_Tdry_LRT      > Lev2c%Range%TPH_Tdry_LRT(mx) )    &
         Lev2c%TPH_Tdry_LRT      = ropp_MDFV

    IF ( Lev2c%TPT_Tdry_LRT      < Lev2c%Range%TPT_Tdry_LRT(mn) .OR. &
         Lev2c%TPT_Tdry_LRT      > Lev2c%Range%TPT_Tdry_LRT(mx) )    &
         Lev2c%TPT_Tdry_LRT      = ropp_MDFV

    IF ( Lev2c%TPH_Tdry_LRT_Flag < Lev2c%Range%TPH_Tdry_LRT_Flag(mn) .OR. &
         Lev2c%TPH_Tdry_LRT_Flag > Lev2c%Range%TPH_Tdry_LRT_Flag(mx) )    &
         Lev2c%TPH_Tdry_LRT_Flag = ropp_MIFV

    IF ( Lev2c%TPH_Tdry_CPT      < Lev2c%Range%TPH_Tdry_CPT(mn) .OR. &
         Lev2c%TPH_Tdry_CPT      > Lev2c%Range%TPH_Tdry_CPT(mx) )    &
         Lev2c%TPH_Tdry_CPT      = ropp_MDFV

    IF ( Lev2c%TPT_Tdry_CPT      < Lev2c%Range%TPT_Tdry_CPT(mn) .OR. &
         Lev2c%TPT_Tdry_CPT      > Lev2c%Range%TPT_Tdry_CPT(mx) )    &
         Lev2c%TPT_Tdry_CPT      = ropp_MDFV

    IF ( Lev2c%TPH_Tdry_CPT_Flag < Lev2c%Range%TPH_Tdry_CPT_Flag(mn) .OR. &
         Lev2c%TPH_Tdry_CPT_Flag > Lev2c%Range%TPH_Tdry_CPT_Flag(mx) )    &
         Lev2c%TPH_Tdry_CPT_Flag = ropp_MIFV

    IF ( Lev2c%PRH_Tdry_CPT      < Lev2c%Range%PRH_Tdry_CPT(mn) .OR. &
         Lev2c%PRH_Tdry_CPT      > Lev2c%Range%PRH_Tdry_CPT(mx) )    &
         Lev2c%PRH_Tdry_CPT      = ropp_MDFV

    IF ( Lev2c%PRT_Tdry_CPT      < Lev2c%Range%PRT_Tdry_CPT(mn) .OR. &
         Lev2c%PRT_Tdry_CPT      > Lev2c%Range%PRT_Tdry_CPT(mx) )    &
         Lev2c%PRT_Tdry_CPT      = ropp_MDFV

    IF ( Lev2c%PRH_Tdry_CPT_Flag < Lev2c%Range%PRH_Tdry_CPT_Flag(mn) .OR. &
         Lev2c%PRH_Tdry_CPT_Flag > Lev2c%Range%PRH_Tdry_CPT_Flag(mx) )    &
         Lev2c%PRH_Tdry_CPT_Flag = ropp_MIFV

    IF ( Lev2c%TPH_Temp_LRT      < Lev2c%Range%TPH_Temp_LRT(mn) .OR. &
         Lev2c%TPH_Temp_LRT      > Lev2c%Range%TPH_Temp_LRT(mx) )    &
         Lev2c%TPH_Temp_LRT      = ropp_MDFV

    IF ( Lev2c%TPT_Temp_LRT      < Lev2c%Range%TPT_Temp_LRT(mn) .OR. &
         Lev2c%TPT_Temp_LRT      > Lev2c%Range%TPT_Temp_LRT(mx) )    &
         Lev2c%TPT_Temp_LRT      = ropp_MDFV

    IF ( Lev2c%TPH_Temp_LRT_Flag < Lev2c%Range%TPH_Temp_LRT_Flag(mn) .OR. &
         Lev2c%TPH_Temp_LRT_Flag > Lev2c%Range%TPH_Temp_LRT_Flag(mx) )    &
         Lev2c%TPH_Temp_LRT_Flag = ropp_MIFV

    IF ( Lev2c%TPH_Temp_CPT      < Lev2c%Range%TPH_Temp_CPT(mn) .OR. &
         Lev2c%TPH_Temp_CPT      > Lev2c%Range%TPH_Temp_CPT(mx) )    &
         Lev2c%TPH_Temp_CPT      = ropp_MDFV

    IF ( Lev2c%TPT_Temp_CPT      < Lev2c%Range%TPT_Temp_CPT(mn) .OR. &
         Lev2c%TPT_Temp_CPT      > Lev2c%Range%TPT_Temp_CPT(mx) )    &
         Lev2c%TPT_Temp_CPT      = ropp_MDFV

    IF ( Lev2c%TPH_Temp_CPT_Flag < Lev2c%Range%TPH_Temp_CPT_Flag(mn) .OR. &
         Lev2c%TPH_Temp_CPT_Flag > Lev2c%Range%TPH_Temp_CPT_Flag(mx) )    &
         Lev2c%TPH_Temp_CPT_Flag = ropp_MIFV

    IF ( Lev2c%PRH_Temp_CPT      < Lev2c%Range%PRH_Temp_CPT(mn) .OR. &
         Lev2c%PRH_Temp_CPT      > Lev2c%Range%PRH_Temp_CPT(mx) )    &
         Lev2c%PRH_Temp_CPT      = ropp_MDFV

    IF ( Lev2c%PRT_Temp_CPT      < Lev2c%Range%PRT_Temp_CPT(mn) .OR. &
         Lev2c%PRT_Temp_CPT      > Lev2c%Range%PRT_Temp_CPT(mx) )    &
         Lev2c%PRT_Temp_CPT      = ropp_MDFV

    IF ( Lev2c%PRH_Temp_CPT_Flag < Lev2c%Range%PRH_Temp_CPT_Flag(mn) .OR. &
         Lev2c%PRH_Temp_CPT_Flag > Lev2c%Range%PRH_Temp_CPT_Flag(mx) )    &
         Lev2c%PRH_Temp_CPT_Flag = ropp_MIFV

! If geopotential is missing set 'profile' to zero points

    IF ( Lev2c%Geop_Sfc < ropp_MDTV ) Lev2c%Npoints = 0

! If there are no valid surface P values, set surface missing flag
! to indicate no valid Level 2c data

    IF ( Lev2c%Npoints > 0                 .AND. &
         Lev2c%Press_Sfc       < ropp_MDTV .AND. &
         Lev2c%Press_Sfc_Sigma < ropp_MDTV .AND. &
         Lev2c%Press_Sfc_Qual  < ropp_MDTV )     &
         Lev2c%Missing = .TRUE.

  END IF

  IF ( Lev2c%Npoints == 0 ) Lev2c%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2ctype

!----------------------------------------------------------------
! 9. Level 2c profile parameters (2d version)
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2ctype_2d ( Lev2c ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2ctype_2d
  USE ropp_utils,    ONLY : ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

! Argument list parameters

  TYPE(L2ctype_2d), INTENT(INOUT) :: Lev2c    ! RO data structure

! Local variables

  INTEGER  :: i    ! Loop indices

  Lev2c%Npoints = MIN ( MAX ( Lev2c%Npoints, 0 ), 1 )

  IF ( Lev2c%Nhoriz <= 0 ) THEN
    Lev2c%Npoints = 0
  ENDIF

  IF ( Lev2c%Npoints > 0 ) THEN

    DO i = 1, Lev2c%Nhoriz

      IF ( Lev2c%Geop_Sfc(i)        < Lev2c%Range%Geop_Sfc(mn) .OR. &
           Lev2c%Geop_Sfc(i)        > Lev2c%Range%Geop_Sfc(mx) )    &
           Lev2c%Geop_Sfc(i)        = ropp_MDFV

      IF ( Lev2c%Press_Sfc(i)       < Lev2c%Range%Press_Sfc(mn) .OR. &
           Lev2c%Press_Sfc(i)       > Lev2c%Range%Press_Sfc(mx) )    &
           Lev2c%Press_Sfc(i)       = ropp_MDFV

      IF ( Lev2c%Press_Sfc_Sigma(i) < Lev2c%Range%Press_Sfc_Sigma(mn) .OR. &
           Lev2c%Press_Sfc_Sigma(i) > Lev2c%Range%Press_Sfc_Sigma(mx) )    &
           Lev2c%Press_Sfc_Sigma(i) = ropp_MDFV

      IF ( Lev2c%Press_Sfc_Qual(i)  < Lev2c%Range%Press_Sfc_Qual(mn) .OR. &
           Lev2c%Press_Sfc_Qual(i)  > Lev2c%Range%Press_Sfc_Qual(mx) )    &
           Lev2c%Press_Sfc_Qual(i)  = ropp_MDFV

! If geopotential is missing set 'profile' to zero points
      IF ( Lev2c%Geop_Sfc(i) < ropp_MDTV ) Lev2c%Npoints = 0

! If there are no valid surface P values, set surface missing flag
! to indicate no valid Level 2c data

      IF ( Lev2c%Geop_Sfc(i)        < ropp_MDTV .AND. &
           Lev2c%Press_Sfc(i)       < ropp_MDTV .AND. &
           Lev2c%Press_Sfc_Sigma(i) < ropp_MDTV .AND. &
           Lev2c%Press_Sfc_Qual(i)  < ropp_MDTV )     &
           Lev2c%Missing = .TRUE.
    ENDDO
  END IF

  IF ( Lev2c%Npoints == 0 ) Lev2c%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2ctype_2d

!----------------------------------------------------------------
! 10. Level 2d profile parameters
!----------------------------------------------------------------

SUBROUTINE ropp_io_rangecheck_l2dtype ( Lev2d ) ! (inout)

  USE ropp_io_types
! USE ropp_io,       not_this => ropp_io_rangecheck_l2dtype
  USE ropp_utils,    ONLY : WHERE, To_Upper, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed parameters

  INTEGER, PARAMETER :: mn   = 1    ! Range attributes: min value in (mn)
  INTEGER, PARAMETER :: mx   = 2    !               and max value in (mx)

  CHARACTER (LEN =  *), PARAMETER :: validchars = &
                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_- " ! Valid chr set

! Argument list parameters

  TYPE(L2dtype), INTENT(INOUT) :: Lev2d    ! RO data structure

! Local variables

  INTEGER  :: i    ! Loop indices

! Holds WHERE output (NB this is NOT the F90 intrinsic)

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx = 0

  CALL To_Upper ( Lev2d%Level_Type )
  DO i = 1, 4
    IF ( INDEX(validchars,Lev2d%Level_Type(i:i)) == 0 ) &
      Lev2d%Level_Type(i:i) = "_"
  END DO
  IF ( Lev2d%Level_Type(1:4) == "____" ) Lev2d%Level_Type = "UNKNOWN"

  IF ( Lev2d%Npoints <= 0 ) THEN
       Lev2d%Npoints =  0

  ELSE

    WHERE ( Lev2d%Level_Coeff_A  < Lev2d%Range%Level_Coeff_A(mn) .OR. &
            Lev2d%Level_Coeff_A  > Lev2d%Range%Level_Coeff_A(mx) )    &
            Lev2d%Level_Coeff_A  = ropp_MDFV

    WHERE ( Lev2d%Level_Coeff_B  < Lev2d%Range%Level_Coeff_B(mn) .OR. &
            Lev2d%Level_Coeff_B  > Lev2d%Range%Level_Coeff_B(mx) )    &
            Lev2d%Level_Coeff_B  = ropp_MDFV

! If there are no valid coefficent values, set coeffs missing flag
! to indicate no valid Level 2d data

    idx => WHERE ( Lev2d%Level_Coeff_A > ropp_MDTV .OR. &
                   Lev2d%Level_Coeff_B > ropp_MDTV, nidx )
    IF ( nidx == 0 ) Lev2d%Missing = .TRUE.

  END IF

  IF ( Lev2d%Npoints == 0 ) Lev2d%Missing = .TRUE.

END SUBROUTINE ropp_io_rangecheck_l2dtype


!----------------------------------------------------------------
! 11. String checking routine
!----------------------------------------------------------------

SUBROUTINE keyword_check ( str ) ! (inout)

! The input string consists of a variable-length keyword
! followed by a space then free-form text. We just check the keyword;
! if all characters are invalid, set whole string as unknown.

  IMPLICIT NONE

! In/out
  CHARACTER (LEN=*), INTENT(inout)    :: str

! Fixed parameters

  INTEGER                             :: i, ikl, ict
  CHARACTER (LEN =*), PARAMETER       :: validchars = &
   "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_- " ! Valid chr set

  str = ADJUSTL(str)
  ikl = INDEX ( str, " " ) - 1
  IF ( ikl > 0 ) THEN
!    CALL To_Upper ( str(1:ikl) )
    ict = 0
    DO i = 1, ikl
      IF ( INDEX(validchars, str(i:i)) == 0 ) THEN
        str(i:i) = "_"
        ict = ict + 1
      END IF
    END DO
    IF ( ict == ikl ) str = "UNKNOWN"
  END IF

END SUBROUTINE keyword_check
