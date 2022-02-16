Interface
SUBROUTINE rttov_dealloc_coefs(  ERR, coefs )
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coefs
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_coefs         ), INTENT(INOUT)            :: coefs        ! coefficients
End Subroutine
End Interface
