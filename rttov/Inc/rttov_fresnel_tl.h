Interface
SUBROUTINE rttov_fresnel_tl( &
            & chanprof,     &
            & profiles,     &
            & coef,         &
            & sunglint,     &
            & sunglint_tl,  &
            & fresnrefl,    &
            & fresnrefl_tl)
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_Type,   &
       & rttov_coef,     &
       & sunglint_type
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof(:)
  TYPE(profile_type  ), INTENT(IN)  :: profiles(:)
  TYPE(rttov_coef    ), INTENT(IN)  :: coef
  TYPE(sunglint_type ), INTENT(IN)  :: sunglint
  TYPE(sunglint_type ), INTENT(IN)  :: sunglint_tl
  REAL(KIND=jprb)     , INTENT(IN)  :: fresnrefl   (size(chanprof))
  REAL(KIND=jprb)     , INTENT(OUT) :: fresnrefl_tl(size(chanprof))
End Subroutine
End Interface
