Interface
SUBROUTINE rttov_profaux_tl( &
            & opts,    &
            & prof,    &
            & prof_tl, &
            & coef,    &
            & aux,     &
            & aux_tl)
  USE rttov_types, ONLY :  &
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(IN)    :: prof_tl(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_tl
End Subroutine
End Interface
