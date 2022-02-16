Interface
SUBROUTINE rttov_pcscores_rec_k( &
            & opts,        &
            & chanprof,    &
            & chanprof_pc, &
            & pccomp,      &
            & pcscores_k,  &
            & coef_pccomp, &
            & clear_k_pc)
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_pccomp,      &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(rttov_options )   , INTENT(IN)    :: opts
  TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof   (:)      ! Channel indices
  TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof_pc(:)      ! Channel indices
  REAL(KIND=jprb)        , INTENT(INOUT) :: pcscores_k (:, :, :)
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  REAL(KIND=jprb)        , INTENT(INOUT) :: clear_k_pc(:, :, :)
End Subroutine
End Interface
