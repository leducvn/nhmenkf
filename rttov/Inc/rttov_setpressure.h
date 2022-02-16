Interface
  subroutine rttov_setpressure (p_sfc, p, ph)
  Use parkind1, Only : jpim     ,jprb
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  implicit none
  Integer (Kind=jpim), parameter :: nlev = 60
  Real (Kind=jprb) :: p_sfc
  Real (Kind=jprb) :: p   (nlev)  , ph  (nlev+1)
End Subroutine
End Interface
