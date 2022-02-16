Interface
SUBROUTINE rttov_fresnel( &
            & chanprof,  &
            & profiles,  &
            & coef,      &
            & sunglint,  &
            & fresnrefl)
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_Type,   &
       & rttov_coef,     &
       & sunglint_type
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type  ), INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
  TYPE(sunglint_type ), INTENT(INOUT) :: sunglint
  REAL(KIND=jprb)     , INTENT(OUT)   :: fresnrefl(size(chanprof))
End Subroutine
End Interface
