SUBROUTINE rttov_init_prof(profiles, p)
! Description:
!   Initialise profile structure
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
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : profile_Type
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(profile_Type), INTENT(INOUT)           :: profiles(:)
  REAL(KIND=jprb)   , INTENT(IN)   , OPTIONAL :: p(:)
!INTF_END
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: nprofiles! Number of profiles
  nprofiles = size(profiles)
  DO j = 1, nprofiles
    IF (Present(p)) THEN
      profiles(j)%p = p
    ELSE
      profiles(j)%p = 0._jprb
    ENDIF
    profiles(j)%t = 0._jprb
    profiles(j)%q = 0._jprb
    IF (Associated(profiles(j)%co2)     ) profiles(j)%co2      = 0._jprb
    IF (Associated(profiles(j)%o3)      ) profiles(j)%o3       = 0._jprb
    IF (Associated(profiles(j)%n2o)     ) profiles(j)%n2o      = 0._jprb
    IF (Associated(profiles(j)%co)      ) profiles(j)%co       = 0._jprb
    IF (Associated(profiles(j)%ch4)     ) profiles(j)%ch4      = 0._jprb
    IF (Associated(profiles(j)%clw)     ) profiles(j)%clw      = 0._jprb
    IF (associated(profiles(j)%aerosols)) profiles(j)%aerosols = 0._jprb
    IF (associated(profiles(j)%cloud)   ) profiles(j)%cloud    = 0._jprb
    IF (associated(profiles(j)%cfrac)   ) profiles(j)%cfrac    = 0._jprb
    IF (associated(profiles(j)%icede)   ) profiles(j)%icede    = 0._jprb
    profiles(j)%id             = " "
    profiles(j)%date           = (/ 1950_jpim, 1_jpim, 1_jpim /)
    profiles(j)%time           = (/ 0_jpim, 0_jpim, 0_jpim /)
    profiles(j)%zenangle       = 0._jprb
    profiles(j)%azangle        = 0._jprb
    profiles(j)%sunzenangle    = 0._jprb
    profiles(j)%sunazangle     = 0._jprb
    profiles(j)%elevation      = 0._jprb
    profiles(j)%latitude       = 0._jprb
    profiles(j)%longitude      = 0._jprb
    profiles(j)%snow_frac      = 0._jprb
    profiles(j)%soil_moisture  = 0._jprb
    profiles(j)%ctp            = 0._jprb
    profiles(j)%cfraction      = 0._jprb
    profiles(j)%s2m%t          = 0._jprb
    profiles(j)%s2m%q          = 0._jprb
    profiles(j)%s2m%o          = 0._jprb
    profiles(j)%s2m%p          = 0._jprb
    profiles(j)%s2m%u          = 0._jprb
    profiles(j)%s2m%v          = 0._jprb
    profiles(j)%s2m%wfetc      = 0._jprb
    profiles(j)%skin%t         = 0._jprb
    profiles(j)%skin%fastem    = 0._jprb
    profiles(j)%skin%surftype  = 0._jprb
    profiles(j)%skin%watertype = 0._jprb 
    profiles(j)%skin%salinity  = 0._jprb
    profiles(j)%idg            =  - 1_jpim
    profiles(j)%ish            =  - 1_jpim
    profiles(j)%Be             = 0._jprb
    profiles(j)%cosbk          = 0._jprb
  ENDDO
END SUBROUTINE
