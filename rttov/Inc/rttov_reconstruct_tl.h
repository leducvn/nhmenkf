Interface
SUBROUTINE rttov_reconstruct_tl( &
            & chanprof_in, &
            & chanprof_pc, &
            & pccomp,      &
            & pccomp_tl,   &
            & coef_pccomp)
  USE rttov_types, ONLY : rttov_chanprof, rttov_pccomp, rttov_coef_pccomp
  IMPLICIT NONE
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_in(:)! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_tl
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
End Subroutine
End Interface
