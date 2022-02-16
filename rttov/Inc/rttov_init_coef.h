Interface
SUBROUTINE rttov_init_coef(  ERR, coef )
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef         ! coefficients
End Subroutine
End Interface
