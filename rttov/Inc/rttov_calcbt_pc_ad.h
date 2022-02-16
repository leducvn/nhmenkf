Interface
SUBROUTINE rttov_calcbt_pc_ad( &
            & chanprof_in, &
            & coef_pccomp, &
            & pccomp,      &
            & pccomp_ad)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef_pccomp, rttov_pccomp
  IMPLICIT NONE
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_pccomp     ), INTENT(IN)    :: pccomp        ! rad%total rad%clear
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_ad
End Subroutine
End Interface
