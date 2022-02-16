!
SUBROUTINE rttov_alloc_scatt_prof ( nprof, cld_profiles, nlev, use_totalice, asw, init)
! Description:
! allocation/deallocation of a RTTOV_SCATT cld_profile structure
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!  1.0      28/04/2010  New F90 code (Alan Geer)
!

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
!INTF_END
  
integer(kind=jpim) :: iprof
logical(kind=jplm) :: init1
real(kind=jprb)    :: ZHOOK_HANDLE

if (lhook) call dr_hook('RTTOV_ALLOC_SCATT_PROF',0_jpim,zhook_handle)

init1 = .false.
if(present(init)) init1 = init

if( asw .eq. 1) then 
  do iprof = 1, nprof

    cld_profiles (iprof) % nlevels      = nlev 
    cld_profiles (iprof) % use_totalice = use_totalice
    cld_profiles (iprof) % cfrac        = 0.0_JPRB    
   
    allocate( cld_profiles (iprof) % ph (nlev+1) )
    allocate( cld_profiles (iprof) % cc (nlev) )
    allocate( cld_profiles (iprof) % clw (nlev) ) 
    allocate( cld_profiles (iprof) % rain (nlev) )
    if (use_totalice) then 
      allocate( cld_profiles (iprof) % totalice (nlev) )
    else
      allocate( cld_profiles (iprof) % sp (nlev) )
      allocate( cld_profiles (iprof) % ciw (nlev) )
    endif

    if( init1) then

      cld_profiles (iprof) % ph   (:) = 0.0_JPRB  
      cld_profiles (iprof) % cc   (:) = 0.0_JPRB  
      cld_profiles (iprof) % clw  (:) = 0.0_JPRB  
      cld_profiles (iprof) % rain (:) = 0.0_JPRB  
      if (use_totalice) then 
        cld_profiles (iprof) % totalice (:) = 0.0_JPRB  
      else
        cld_profiles (iprof) % sp       (:) = 0.0_JPRB  
        cld_profiles (iprof) % ciw      (:) = 0.0_JPRB  
      endif

    endif
  enddo
else
  do iprof = 1, nprof

    deallocate( cld_profiles (iprof) % ph )
    deallocate( cld_profiles (iprof) % cc )
    deallocate( cld_profiles (iprof) % clw )
    deallocate( cld_profiles (iprof) % rain )
    if (cld_profiles (iprof) % use_totalice) then 
      deallocate( cld_profiles (iprof) % totalice )
    else
      deallocate( cld_profiles (iprof) % ciw )
      deallocate( cld_profiles (iprof) % sp )
    endif
    
  enddo  
endif

if (lhook) call dr_hook('RTTOV_ALLOC_SCATT_PROF',1_jpim,zhook_handle)

end subroutine

!- End of header --------------------------------------------------------
