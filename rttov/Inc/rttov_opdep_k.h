Interface
SUBROUTINE rttov_opdep_k( &
            & nlayers,      &
            & chanprof,     &
            & predictors,   &
            & predictors_k, &
            & aux,          &
            & aux_k,        &
            & coef,         &
            & opdp_path,    &
            & opdp_path_k,  &
            & opdp_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(INOUT) :: aux_k
  TYPE(opdp_path_Type ), INTENT(IN)    :: opdp_path                        ! optical depths
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_k
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, size(chanprof))
End Subroutine
End Interface
