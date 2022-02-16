Interface
SUBROUTINE rttov_write_ascii_pccoef( &
            & err,         &
            & coef_pccomp, &
            & file_lu)
#include "throw.h"
  USE rttov_types, ONLY : rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)     , INTENT(IN)  :: file_lu    ! file logical unit number
  TYPE(rttov_coef_pccomp), INTENT(IN)  :: coef_pccomp
  INTEGER(KIND=jpim)     , INTENT(OUT) :: err        ! return code
End Subroutine
End Interface
