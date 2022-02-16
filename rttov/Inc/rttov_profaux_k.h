Interface
SUBROUTINE rttov_profaux_k( &
            & opts,       &
            & chanprof,   &
            & profiles,   &
            & profiles_k, &
            & coef,       &
            & aux_prof,   &
            & aux_prof_k)
  USE rttov_types, ONLY :  &
       & rttov_options,  &
       & rttov_chanprof, &
       & rttov_coef,     &
       & profile_Type,   &
       & profile_aux
  IMPLICIT NONE
  TYPE(rttov_options ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(profile_aux   ), INTENT(IN)    :: aux_prof
  TYPE(profile_aux   ), INTENT(INOUT) :: aux_prof_k
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
End Subroutine
End Interface
