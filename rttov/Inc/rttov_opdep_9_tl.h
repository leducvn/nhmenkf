Interface
SUBROUTINE rttov_opdep_9_tl( &
            & nlayers,       &
            & chanprof,      &
            & predictors,    &
            & predictors_tl, &
            & aux,           &
            & aux_tl,        &
            & coef,          &
            & opdp_path,     &
            & opdp_path_tl,  &
            & opdp_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  TYPE(rttov_chanprof) , INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(IN)    :: predictors_tl
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(IN)    :: aux_tl
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path                        ! optical depths
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_tl
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, size(chanprof))
End Subroutine
End Interface
