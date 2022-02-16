subroutine rttov_write_coefs( err, coefs, opts, &
                             form_coef,     &
                             form_scaer,    &
                             form_sccld,    &
                             form_pccoef,   &
                             file_coef,     &
                             file_scaer,    &
                             file_sccld,    &
                             file_pccoef,   &
                             file_id_coef,  &
                             file_id_scaer, &
                             file_id_sccld, &
                             file_id_pccoef,&
                             instrument)
! Description:
!   Write the coefficients structure
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
#include "throw.h"

Use rttov_types, Only : rttov_coefs, rttov_options
Use parkind1, Only : jpim
!INTF_OFF
Use parkind1, Only : jprb, jplm
Use rttov_coef_io_mod, only : getlun, closelun
Use yomhook, only : LHOOK, DR_HOOK
!INTF_ON

Implicit None

Integer(Kind=jpim), Intent(out) :: err
Type(rttov_coefs),  Intent(in)  :: coefs
Type(rttov_options),Intent(in)  :: opts
Character(len=*),   Intent(in), Optional :: form_coef
Character(len=*),   Intent(in), Optional :: form_scaer
Character(len=*),   Intent(in), Optional :: form_sccld
Character(len=*),   Intent(in), Optional :: form_pccoef
Character(len=*),   Intent(in), Optional :: file_coef
Character(len=*),   Intent(in), Optional :: file_scaer
Character(len=*),   Intent(in), Optional :: file_sccld
Character(len=*),   Intent(in), Optional :: file_pccoef
Integer(Kind=jpim), Intent(in), Optional :: file_id_coef
Integer(Kind=jpim), Intent(in), Optional :: file_id_scaer
Integer(Kind=jpim), Intent(in), Optional :: file_id_sccld
Integer(Kind=jpim), Intent(in), Optional :: file_id_pccoef
Integer(Kind=jpim), Intent(in), Optional :: instrument(3)

!INTF_END

#include "rttov_errorreport.h"
#include "rttov_cmpuc.h"


#include "rttov_write_ascii_coef.h"
#include "rttov_write_ascii_pccoef.h"
#include "rttov_write_ascii_scaercoef.h"
#include "rttov_write_ascii_sccldcoef.h"
#include "rttov_write_binary_coef.h"
#include "rttov_write_binary_pccoef.h"
#include "rttov_write_binary_scaercoef.h"
#include "rttov_write_binary_sccldcoef.h"


Integer(Kind=jpim) :: file_id_coef1
Integer(Kind=jpim) :: file_id_scaer1
Integer(Kind=jpim) :: file_id_sccld1
Integer(Kind=jpim) :: file_id_pccoef1

character(len=32) :: form_coef1
character(len=32) :: form_scaer1
character(len=32) :: form_sccld1
character(len=32) :: form_pccoef1

Real(Kind=jprb)   :: ZHOOK_HANDLE

TRY

!
IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 0_jpim, ZHOOK_HANDLE)

call getlun(err, file_id_coef1, form_coef1, .true._jplm, "rtcoef", file_id_coef, file_coef, form_coef, instrument )
THROW(err.ne.0)

if (rttov_cmpuc(form_coef1,"unformatted")) then
  call rttov_write_binary_coef (err, coefs%coef, file_id_coef1)
  THROW(err.ne.0)
else if (rttov_cmpuc(form_coef1,"formatted")) then
  call rttov_write_ascii_coef (err, coefs%coef, file_id_coef1)
  THROW(err.ne.0)
else
  err = errorstatus_fatal
  THROWM(err.ne.0,"Unknown format "//Trim(form_coef1))
endif

call closelun(err, file_id_coef1, file_id_coef)
THROW(err.ne.0)

!


if(opts%addaerosl) then
  call getlun(err, file_id_scaer1, form_scaer1, .true._jplm, "scaer", file_id_scaer, &
              file_scaer, form_scaer, instrument )
  THROW(err.ne.0)
  if (rttov_cmpuc(form_scaer1,"unformatted")) then
    call rttov_write_binary_scaercoef (err, coefs%coef, coefs%coef_scatt_ir, &
                                      coefs%optp, file_id_scaer1)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_scaer1,"formatted")) then
    call rttov_write_ascii_scaercoef (err, coefs%coef, coefs%coef_scatt_ir, &
                                     coefs%optp, file_id_scaer1)
    THROW(err.ne.0)
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_scaer1))
  endif

  call closelun(err, file_id_scaer1, file_id_scaer)
  THROW(err.ne.0)

endif

!

if(opts%addclouds) then
  call getlun(err, file_id_sccld1, form_sccld1, .true._jplm, "sccld", file_id_sccld, file_sccld, form_sccld, instrument)
  THROW(err.ne.0)
  if (rttov_cmpuc(form_sccld1,"unformatted")) then
    call rttov_write_binary_sccldcoef (err, coefs%coef, coefs%coef_scatt_ir, &                            
                                      coefs%optp, file_id_sccld1)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_sccld1,"formatted")) then
    call rttov_write_ascii_sccldcoef (err, coefs%coef, coefs%coef_scatt_ir, &                             
                                     coefs%optp, file_id_sccld1)
    THROW(err.ne.0)                                                                                      
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_sccld1))
  endif

  call closelun(err, file_id_sccld1, file_id_sccld)
  THROW(err.ne.0)

endif

!

if(opts%addpc) then
  call getlun(err, file_id_pccoef1, form_pccoef1, .true._jplm, "pccoef", file_id_pccoef, file_pccoef, form_pccoef, &
    instrument)
  THROW(err.ne.0)
  if (rttov_cmpuc(form_pccoef1,"unformatted")) then
    call rttov_write_binary_pccoef (err, coefs%coef_pccomp, file_id_pccoef1)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_pccoef1,"formatted")) then
    call rttov_write_ascii_pccoef (err,  coefs%coef_pccomp, file_id_pccoef1)
    THROW(err.ne.0)    
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_pccoef1))
  endif

  call closelun(err, file_id_pccoef1, file_id_pccoef)
  THROW(err.ne.0)

endif

IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)
CATCH
IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)

end subroutine
