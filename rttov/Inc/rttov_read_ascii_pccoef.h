Interface
SUBROUTINE rttov_read_ascii_pccoef( &
            & err,           &
            & coef,          &
            & coef_pccomp,   &
            & file_lu,       &
            & channels,      &
            & channels_rec)! in Optional
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu        ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels    (:)! list of channels to extract
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels_rec(:)
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  TYPE(rttov_coef_pccomp  ), INTENT(INOUT)             :: coef_pccomp
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: err           ! return code
End Subroutine
End Interface
