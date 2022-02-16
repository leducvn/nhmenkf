Interface
SUBROUTINE rttov_alloc_rad( &
            & ERR,       &
            & nchannels, &
            & radiance,  &
            & nlayers,   &
            & asw,       &
            & init)
  USE rttov_types, ONLY : radiance_type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim) , INTENT(OUT)          :: ERR      ! return code
  INTEGER(KIND=jpim) , INTENT(IN)           :: nchannels! number of channels
  INTEGER(KIND=jpim) , INTENT(IN)           :: nlayers  ! number of levels
  TYPE(radiance_Type), INTENT(INOUT)        :: radiance ! radiances
  INTEGER(KIND=jpim) , INTENT(IN)           :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm) , OPTIONAL, INTENT(IN) :: init
End Subroutine
End Interface
