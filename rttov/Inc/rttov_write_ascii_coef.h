Interface
SUBROUTINE rttov_write_ascii_coef(err, coef, file_id)
#include "throw.h"
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(IN)  :: file_id! file logical unit number
  TYPE(rttov_coef)  , INTENT(IN)  :: coef   ! coefficients
  INTEGER(KIND=jpim), INTENT(OUT) :: err    ! return code
End Subroutine
End Interface
