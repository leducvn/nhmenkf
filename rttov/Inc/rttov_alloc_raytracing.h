Interface
SUBROUTINE rttov_alloc_raytracing( &
            & ERR,          &
            & nraytracings, &
            & raytracings,  &
            & nlevels,      &
            & asw,          &
            & init)
  USE rttov_types, ONLY : raytracing_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(OUT)             :: ERR         ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nraytracings
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nlevels
  TYPE(raytracing_Type), INTENT(INOUT)           :: raytracings
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
