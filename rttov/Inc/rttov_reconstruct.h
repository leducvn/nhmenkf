Interface
SUBROUTINE rttov_reconstruct( &
            & chanprof_in, &
            & chanprof_pc, &
            & pccomp,      &
            & coef_pccomp)
  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp
  IMPLICIT NONE
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
End Subroutine
End Interface
