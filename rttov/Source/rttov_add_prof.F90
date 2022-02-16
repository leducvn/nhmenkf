SUBROUTINE rttov_add_prof( &
            & profiles,  &
            & profiles1, &
            & profiles2, &
            & lair,      &
            & lground)
! Description:
!   Adds two profile structures
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
  USE rttov_types, ONLY : profile_type
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  TYPE(profile_type), INTENT(INOUT)           :: profiles (:)
  TYPE(profile_type), INTENT(IN)              :: profiles1(size(profiles))
  TYPE(profile_type), INTENT(IN)              :: profiles2(size(profiles))
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: lair
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: lground
!INTF_END
  LOGICAL(KIND=jplm) :: lair1
  LOGICAL(KIND=jplm) :: lground1
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  lground1 = .TRUE.
  IF (Present(lground)) lground1 = lground
  lair1 = .TRUE.
  IF (Present(lair)) lair1 = lair
  nprofiles = size(profiles1)
! do not add 
! profiles%id %date %time
!
  DO iprof = 1, nprofiles
    IF (lground1) THEN
      profiles(iprof)%zenangle       = profiles1(iprof)%zenangle + profiles2(iprof)%zenangle
      profiles(iprof)%azangle        = profiles1(iprof)%azangle + profiles2(iprof)%azangle
      profiles(iprof)%sunzenangle    = profiles1(iprof)%sunzenangle + profiles2(iprof)%sunzenangle
      profiles(iprof)%sunazangle     = profiles1(iprof)%sunazangle + profiles2(iprof)%sunazangle
      profiles(iprof)%latitude       = profiles1(iprof)%latitude + profiles2(iprof)%latitude
      profiles(iprof)%longitude      = profiles1(iprof)%longitude + profiles2(iprof)%longitude
      profiles(iprof)%snow_frac      = profiles1(iprof)%snow_frac + profiles2(iprof)%snow_frac
      profiles(iprof)%soil_moisture  = profiles1(iprof)%soil_moisture + profiles2(iprof)%soil_moisture
      profiles(iprof)%ctp            = profiles1(iprof)%ctp + profiles2(iprof)%ctp
      profiles(iprof)%cfraction      = profiles1(iprof)%cfraction + profiles2(iprof)%cfraction
      profiles(iprof)%elevation      = profiles1(iprof)%elevation + profiles2(iprof)%elevation
      profiles(iprof)%S2m%t          = profiles1(iprof)%S2m%t + profiles2(iprof)%S2m%t
      profiles(iprof)%S2m%q          = profiles1(iprof)%S2m%q + profiles2(iprof)%S2m%q
      profiles(iprof)%S2m%o          = profiles1(iprof)%S2m%o + profiles2(iprof)%S2m%o
      profiles(iprof)%S2m%p          = profiles1(iprof)%S2m%p + profiles2(iprof)%S2m%p
      profiles(iprof)%S2m%u          = profiles1(iprof)%S2m%u + profiles2(iprof)%S2m%u
      profiles(iprof)%S2m%v          = profiles1(iprof)%S2m%v + profiles2(iprof)%S2m%v
      profiles(iprof)%S2m%wfetc      = profiles1(iprof)%S2m%wfetc + profiles2(iprof)%S2m%wfetc
      profiles(iprof)%skin%t         = profiles1(iprof)%skin%t + profiles2(iprof)%skin%t
      profiles(iprof)%skin%fastem(:) = profiles1(iprof)%skin%fastem(:) + profiles2(iprof)%skin%fastem(:)
      profiles(iprof)%skin%salinity  = profiles1(iprof)%skin%salinity + profiles2(iprof)%skin%salinity
      profiles(iprof)%be             = profiles1(iprof)%be + profiles2(iprof)%be
      profiles(iprof)%cosbk          = profiles1(iprof)%cosbk + profiles2(iprof)%cosbk
    ENDIF
    IF (lair1) THEN
      profiles(iprof)%p = profiles1(iprof)%p + profiles2(iprof)%p
      profiles(iprof)%t = profiles1(iprof)%t + profiles2(iprof)%t
      profiles(iprof)%q = profiles1(iprof)%q + profiles2(iprof)%q
      IF (Associated(profiles(iprof)%o3) .AND. Associated(profiles1(iprof)%o3) .AND. Associated(profiles2(iprof)%o3)     &
        &                   ) profiles(iprof)%o3       = profiles1(iprof)%o3 + profiles2(iprof)%o3
      IF (Associated(profiles(iprof)%co2) .AND. Associated(profiles1(iprof)%co2) .AND. Associated(profiles2(iprof)%co2)     &
        &                ) profiles(iprof)%co2      = profiles1(iprof)%co2 + profiles2(iprof)%co2
      IF (Associated(profiles(iprof)%n2o) .AND. Associated(profiles1(iprof)%n2o) .AND. Associated(profiles2(iprof)%n2o)     &
        &                ) profiles(iprof)%n2o      = profiles1(iprof)%n2o + profiles2(iprof)%n2o
      IF (Associated(profiles(iprof)%co) .AND. Associated(profiles1(iprof)%co) .AND. Associated(profiles2(iprof)%co)     &
        &                   ) profiles(iprof)%co       = profiles1(iprof)%co + profiles2(iprof)%co
      IF (Associated(profiles(iprof)%ch4) .AND. Associated(profiles1(iprof)%ch4) .AND. Associated(profiles2(iprof)%ch4)     &
        &                ) profiles(iprof)%ch4      = profiles1(iprof)%ch4 + profiles2(iprof)%ch4
      IF (Associated(profiles(iprof)%clw) .AND. Associated(profiles1(iprof)%clw) .AND. Associated(profiles2(iprof)%clw)     &
        &                ) profiles(iprof)%clw      = profiles1(iprof)%clw + profiles2(iprof)%clw
      IF (Associated(profiles(iprof)%aerosols) .AND. Associated(profiles1(iprof)%aerosols) .AND.      &
        & Associated(profiles2(iprof)%aerosols))                                                      &
        &  profiles(iprof)%aerosols = profiles1(iprof)%aerosols + profiles2(iprof)%aerosols
      IF (Associated(profiles(iprof)%cloud) .AND. Associated(profiles1(iprof)%cloud) .AND.      &
        & Associated(profiles2(iprof)%cloud)         )                                          &
        &  profiles(iprof)%cloud    = profiles1(iprof)%cloud + profiles2(iprof)%cloud
      IF (Associated(profiles(iprof)%cfrac) .AND. Associated(profiles1(iprof)%cfrac) .AND.      &
        & Associated(profiles2(iprof)%cfrac)         )                                          &
        &  profiles(iprof)%cfrac    = profiles1(iprof)%cfrac + profiles2(iprof)%cfrac
      IF (Associated(profiles(iprof)%icede) .AND. Associated(profiles1(iprof)%icede) .AND.      &
        & Associated(profiles2(iprof)%icede)         )                                          &
        &  profiles(iprof)%icede    = profiles1(iprof)%icede + profiles2(iprof)%icede
    ENDIF
  ENDDO
END SUBROUTINE 
