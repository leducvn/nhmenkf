Interface
SUBROUTINE rttov_alloc_scatt_prof ( nprof, cld_profiles, nlev, use_totalice, asw, init)
USE YOMHOOK,  ONLY: LHOOK , DR_HOOK
use parkind1, only: jpim, jprb, jplm
use rttov_types, only : profile_cloud_type 
IMPLICIT NONE
integer(kind=jpim), intent(in) :: nlev         ! number of levels
integer(kind=jpim), intent(in) :: nprof        ! number of profiles
integer(kind=jpim), intent(in) :: asw          ! 1=allocate,      0=deallocate
logical(kind=jplm), optional,  intent(in) :: init         ! true=zero contents, false=don't bother
logical(kind=jplm), intent(in) :: use_totalice ! Choose separate ciw and snow, or totalice
type(profile_cloud_type), intent (inout) :: cld_profiles (nprof)
End Subroutine
End Interface
