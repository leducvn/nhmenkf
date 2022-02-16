SUBROUTINE rttov_mult_profiles_k(profiles_k_rec, profiles_k, clear_k_pc)
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
  USE rttov_types, ONLY : profile_type
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(profile_type), INTENT(INOUT) :: profiles_k_rec(:)
  TYPE(profile_type), INTENT(IN)    :: profiles_k   (:)
  REAL(KIND=jprb)   , INTENT(IN)    :: clear_k_pc   (:, :, :)! dBT(PC)/dRadiance(RTTOV)
!INTF_END
  INTEGER(KIND=jpim) :: nchannels_rec_p
  INTEGER(KIND=jpim) :: nchannels_p
  INTEGER(KIND=jpim) :: ichan_in1
  INTEGER(KIND=jpim) :: ichan1
  INTEGER(KIND=jpim) :: ichan_in
  INTEGER(KIND=jpim) :: ichan
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: nprofiles
  nchannels_rec_p = size(clear_k_pc, 1)
  nchannels_p     = size(clear_k_pc, 2)
  nprofiles       = size(clear_k_pc, 3)

  DO iprof = 1, nprofiles
    DO ichan_in = 1, nchannels_rec_p
      ichan_in1 = ichan_in + (iprof - 1) * nchannels_rec_p
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        profiles_k_rec(ichan_in1)%p   =      &
          & profiles_k_rec(ichan_in1)%p + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%p
        profiles_k_rec(ichan_in1)%t   =      &
          & profiles_k_rec(ichan_in1)%t + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%t
        profiles_k_rec(ichan_in1)%q   =      &
          & profiles_k_rec(ichan_in1)%q + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%q
        IF (associated(profiles_k_rec(ichan_in1)%co2) ) profiles_k_rec(ichan_in1)%co2 =     &
   & profiles_k_rec(ichan_in1)%co2 + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%co2
        IF (associated(profiles_k_rec(ichan_in1)%o3) ) profiles_k_rec(ichan_in1)%o3  =      &
   & profiles_k_rec(ichan_in1)%o3 + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%o3
        IF (associated(profiles_k_rec(ichan_in1)%n2o)) profiles_k_rec(ichan_in1)%n2o =      &
   & profiles_k_rec(ichan_in1)%n2o + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%n2o
        IF (associated(profiles_k_rec(ichan_in1)%co) ) profiles_k_rec(ichan_in1)%co  =      &
   & profiles_k_rec(ichan_in1)%co + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%co
        IF (associated(profiles_k_rec(ichan_in1)%ch4)) profiles_k_rec(ichan_in1)%ch4 =      &
   & profiles_k_rec(ichan_in1)%ch4 + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%ch4
        IF (associated(profiles_k_rec(ichan_in1)%clw)) profiles_k_rec(ichan_in1)%clw =      &
   & profiles_k_rec(ichan_in1)%clw + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%clw
        profiles_k_rec(ichan_in1)%ctp       =      &
          & profiles_k_rec(ichan_in1)%ctp + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%ctp
        profiles_k_rec(ichan_in1)%cfraction =      &
          & profiles_k_rec(ichan_in1)%cfraction + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%cfraction
        profiles_k_rec(ichan_in1)%s2m%t     =      &
          & profiles_k_rec(ichan_in1)%s2m%t + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%t
        profiles_k_rec(ichan_in1)%s2m%q     =      &
          & profiles_k_rec(ichan_in1)%s2m%q + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%q
        profiles_k_rec(ichan_in1)%s2m%o     =      &
          & profiles_k_rec(ichan_in1)%s2m%o + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%o
        profiles_k_rec(ichan_in1)%s2m%p     =      &
          & profiles_k_rec(ichan_in1)%s2m%p + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%p
        profiles_k_rec(ichan_in1)%s2m%u     =      &
          & profiles_k_rec(ichan_in1)%s2m%u + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%u
        profiles_k_rec(ichan_in1)%s2m%v     =      &
          & profiles_k_rec(ichan_in1)%s2m%v + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%v
        profiles_k_rec(ichan_in1)%s2m%wfetc =      &
          & profiles_k_rec(ichan_in1)%s2m%wfetc + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%s2m%wfetc
        profiles_k_rec(ichan_in1)%skin%t    =      &
          & profiles_k_rec(ichan_in1)%skin%t + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%skin%t
        IF (associated(profiles_k_rec(ichan_in1)%aerosols)) THEN
          profiles_k_rec(ichan_in1)%aerosols =      &
            & profiles_k_rec(ichan_in1)%aerosols + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%aerosols
        ENDIF
        IF (associated(profiles_k_rec(ichan_in1)%cloud)) THEN
          profiles_k_rec(ichan_in1)%cloud =      &
            & profiles_k_rec(ichan_in1)%cloud + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%cloud
        ENDIF
        IF (associated(profiles_k_rec(ichan_in1)%cfrac)) THEN
          profiles_k_rec(ichan_in1)%cfrac =      &
            & profiles_k_rec(ichan_in1)%cfrac + clear_k_pc(ichan_in, ichan, iprof) * profiles_k(ichan1)%cfrac
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE 
