! $Id: datetimeprogs.f90 3696 2013-06-17 08:48:37Z idculv $

!****m* DateTime/DateTimeProgs *
!
! NAME
!    DateTimeProgs    (datetimeprogs.f90)
!
! SYNOPSIS
!   Interface module for date & time manipulation routines
!
!    USE DateTimeProgs
!
! DESCRIPTION
!    This module provides interfaces to some date and time manipulation 
!    routines
!
! SEE ALSO
!   DateTimeTypes
!   CalToJul
!   ConvertDT
!   Date_and_Time_UTC
!   DateTimeOffset
!   MonthOfYear
!   TimeSince
!
! AUTHOR
!    D. Offiler, Met Office
!
! COPYRIGHT
!
!    Crown copyright 2013, Met Office. All rights reserved.
!
!    Use, duplication or disclosure of this code is subject to the restrictions
!    as set forth in the contract. If no contract has been raised with this
!    copy of the code, the use, duplication or disclosure of it is strictly
!    prohibited. Permission to do so must first be obtained in writing from
!    the Head of Satellite Applications at the following address:
!       Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB  United Kingdom
!
!****

MODULE DateTimeProgs

  use DateTimeTypes

! Interfaces

  interface
    subroutine CalToJul (CDT,JDF,ind)
      use datetimetypes
      implicit none
      integer,  intent(inout) :: CDT(:)
      real(dp), intent(inout) :: JDF
      integer,  intent(in)    :: ind
    end subroutine CalToJul
  end interface

  interface
    subroutine ConvertDT (DTlong,DTshort,ind,LongFmt,MSec)
      implicit none
      character (len=*), intent(inout) :: DTlong
      character (len=*), intent(inout) :: DTshort
      integer,           intent(in)    :: ind
      integer, optional, intent(in)    :: LongFmt
      integer, optional, intent(inout) :: Msec
    end subroutine ConvertDT
  end interface 

  interface
    subroutine Date_and_Time_UTC (Date,Time,Zone,Values)
      implicit none
      character (len=*), optional, intent(out)  :: Date
      character (len=*), optional, intent(out)  :: Time
      character (len=*), optional, intent(out)  :: Zone
      integer,           optional, intent(out)  :: Values(:)
    end subroutine Date_and_Time_UTC
  end interface

  interface
    subroutine DateTimeOffset (CDT1,oper,CDT2)
      implicit none
      character,         intent(in)    :: oper
      integer,           intent(inout) :: CDT1(:)
      integer,           intent(inout) :: CDT2(:)
    end subroutine DateTimeOffset
  end interface

  interface
    subroutine MonthOfYear (MonthNumber,MonthString,inv)
      implicit none
      integer,           intent(inout) :: MonthNumber
      character (len=*), intent(inout) :: MonthString
      integer,           intent(in)    :: inv
    end subroutine MonthOfYear
  end interface

  interface
    subroutine TimeSince (CDT, Tsince, inv, Base)
      use datetimetypes
      implicit none
      integer,                     intent(inout) :: CDT(:)
      real(dp),                    intent(inout) :: Tsince
      integer,                     intent(in)    :: inv
      character (len=*), optional, intent(in)    :: Base
    end subroutine TimeSince
  end interface
  
END MODULE DateTimeProgs
