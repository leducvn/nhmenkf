Interface
SUBROUTINE rttov_read_binary_pccoef( &
            & ERR,           &
            & coef,          &
            & coef_pccomp,   &
            & file_lu,       &
            & channels,      &
            & channels_rec)
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_optpar_ir,     &
       & rttov_coef_scatt_ir, &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu        ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels    (:)! list of channels to extract
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels_rec(:)
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  TYPE(rttov_coef_pccomp  ), INTENT(INOUT)             :: coef_pccomp
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
End Subroutine
End Interface
