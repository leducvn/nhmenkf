Interface
SUBROUTINE rttov_reconstruct_k( &
            & chanprof_in, &
            & chanprof_pc, &
            & pccomp,      &
            & pccomp_k,    &
            & pcscores_k,  &
            & coef_pccomp)
  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp
  USE parkind1, ONLY : jprb
  USE PARKIND1, ONLY : JPRB
  IMPLICIT NONE
  TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof_in(:)      ! Channel indices
  TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof_pc(:)
  REAL(KIND=jprb)        , INTENT(INOUT) :: pcscores_k (:, :, :)
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_k
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
End Subroutine
End Interface
