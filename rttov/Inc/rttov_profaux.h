Interface
SUBROUTINE rttov_profaux( &
            & opts, &
            & prof, &
            & coef, &
            & aux)
  USE rttov_types, ONLY :  &
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof(:)! profile
  TYPE(rttov_coef   ), INTENT(IN)    :: coef   ! coefficients
  TYPE(profile_aux  ), INTENT(INOUT) :: aux    ! auxilary profile info
End Subroutine
End Interface
