Interface
SUBROUTINE rttov_init_coef_pccomp( &
            & ERR,           &
            & coef,          &
            & coef_pccomp)
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef         ! coefficients
  TYPE(rttov_coef_pccomp  ), INTENT(INOUT)             :: coef_pccomp
End Subroutine
End Interface
