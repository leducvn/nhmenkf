module rttov_coef_io_mod
! Description:
!   Useful functions for I/O
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
use parkind1, Only : jpim, jprb, jplm

#include "throw.h"

Implicit None

#include "rttov_errorreport.h"
#include "rttov_cmpuc.h"


contains


!
! open a coefficient file; if lun was passed as input, use it, otherwise open the file
! whose name was passed in f, otherwise create the coefficient file name from the instrument
! definition
!
subroutine getlun( err, lun1, form1, loutput, t, lun, f, form, instrument )

  USE rttov_const, ONLY :  &
       & rttov_magic_string

integer(kind=jpim), intent(out) :: err
integer(Kind=jpim), intent(out) :: lun1    ! logical unit the file was opened with
character(len=*),   intent(out) :: form1   ! form (= unformatted, formatted) the file was opened with
logical(Kind=jplm), intent(in)  :: loutput 
character(len=*),   intent(in)  :: t       ! type = "rtcoef", "pccoef", ...
integer(Kind=jpim), intent(in), optional :: lun
character(len=*),   intent(in), optional :: f
character(len=*),   intent(in), optional :: form
integer(kind=jpim), intent(in), optional :: instrument(3)
!

#include "rttov_coeffname.h"

Logical(Kind=jplm) :: lopen, lexist
character(256) :: f1

CHARACTER(LEN = 16) :: bin_check_string

TRY

form1 = ""
if(present(form))then
  if(form.ne."") then
    form1 = form
  endif
endif

!
! if lun was passed as argument, then take it and return
!
if(present(lun)) then
  lun1 = lun
!
! if form was passed as argument, then take it,
! otherwise make a guess
!
  if(form.ne."")then
    form1 = form
  else
    inquire(lun1, form = form1)
  endif
  return
endif

!
! Try to guess the filename
!
!
! If f was passed as argument, then take it
!
if(present(f)) then
  f1 = f
  if (loutput.and.(form1.eq."")) then
    err = errorstatus_fatal
    THROWM(err.ne.0,"Format argument missing for output file")
  endif

!
! Guess the filename from the instrument triplet
!
else if(present(instrument)) then

  if (loutput) then

!
! Formatted output is the default
!
    if ((.not.rttov_cmpuc ("unformatted", form1)).and.(.not.rttov_cmpuc ("formatted", form1))) then
      form1 = "formatted"
    endif
    call rttov_coeffname(err, instrument, f1, t, rttov_cmpuc ("unformatted", form1))
    THROW(err.ne.0)

  else if (form1 .ne. "") then

    call rttov_coeffname(err, instrument, f1, t, rttov_cmpuc ("unformatted", form1))
    THROW(err.ne.0)

  else

!
! Guess the input format (check for file existence)
!

   form1 = "unformatted"
   call rttov_coeffname(err, instrument, f1, t, .true._jplm)
   THROW(err.ne.0)
   Inquire (file = f1, exist = lexist)

   If(.not.lexist) then
     form1 = "formatted"
     call rttov_coeffname(err, instrument, f1, t, .false._jplm)
     THROW(err.ne.0)
   endif

  endif

else
  err = errorstatus_fatal
  THROWM(err.ne.0,"Opening the coefficient file requires the filename or the instrument id")
endif

!
! Look for a free logical unit (not thread-safe)
!
do lun1 = 9, 99
  Inquire( lun1, opened = lopen )
  If(.Not. lopen) Exit
End do
If(lopen) then
  err = errorstatus_fatal
  THROWM(err.ne.0,"Cannot find a free lun")
endif


if (loutput) then
  open (lun1, file = f1, form = form1, status = 'replace', action = 'write', iostat = err)
  THROWM(err.ne.0,"Cannot open "//Trim(f1))
else

  if (form1.eq."") then

!
! Guess the input file format with the embedded magic string
!

    open (lun1, file = f1, form = 'unformatted', status = 'old', action = 'read', iostat = err)
    THROWM(err.ne.0,"Cannot open "//Trim(f1))
    READ (lun1, iostat=err) bin_check_string
    close (lun1, iostat = err)
    THROW(err.ne.0)

    if (bin_check_string /= rttov_magic_string) then
      form1 = "formatted"
    else
      form1 = "unformatted"
    endif

  endif

  open (lun1, file = f1, form = form1, status = 'old', action = 'read', iostat = err)
  THROWM(err.ne.0,"Cannot open "//Trim(f1))
endif

CATCH

end subroutine

subroutine closelun(err, lun1, lun)
integer(kind=jpim), intent(out) :: err
integer(kind=jpim), intent(in) :: lun1
integer(kind=jpim), intent(in), optional :: lun


TRY

if(.not.present(lun)) then
  close (lun1, iostat = err)
  THROW(err.ne.0)
endif


CATCH

end subroutine


end module
