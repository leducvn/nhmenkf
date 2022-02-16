Interface
Subroutine rttov_readscattcoeffs  (&
      & errorstatus,   &! out
      & coef_rttov,    &! in
      & coef_scatt,    &! inout
      & kmyproc,       &! my proc
      & knproc,        &! number of processors
      & kioproc,       &! proc for io
      & file_id       ) ! in  
  Use rttov_types, Only : &
       & rttov_coef, &
       & rttov_scatt_coef 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  Use rttov_const, Only :   &
       & inst_name           ,&
       & platform_name       ,&
       & errorstatus_success ,&
       & errorstatus_fatal   ,&
       & lensection  
  Use parkind1, Only : jpim, jprb, jplm
  Implicit None
  Integer(Kind=jpim), Optional, Intent(in) :: kmyproc  ! logical processor id
  Integer(Kind=jpim), Optional, Intent(in) :: knproc   ! Number of processors
  Integer(Kind=jpim), Optional, Intent(in) :: kioproc  ! processor  dedicated for io
  Integer(Kind=jpim), Optional, Intent(in) :: file_id  ! file logical unit number
  Integer(Kind=jpim), Intent (out) :: errorstatus      ! return code
  Type( rttov_coef ), Intent (in) :: coef_rttov        ! clear-sky coefficients
  Type( rttov_scatt_coef ), Intent (inout) :: coef_scatt ! coefficients
End Subroutine
End Interface
