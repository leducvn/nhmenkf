!
SUBROUTINE rttov_calcemis_ir_ad( &
            & profiles,      &
            & profiles_ad,   &
            & coef,          &
            & addpc,         &
            & coef_pccomp,   &
            & chanprof,      &
            & calcemis,      &
            & emissivity_ad)
! Description:
! To compute IR surface emissivities for all channels and all
! profiles if desired
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method:
! RTTOV-6 IR surface emissivity report, V. Sherlock at:
! http://www.metoffice.com/research/interproj/nwpsaf/rtm/papers/isem6.pdf
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       02/01/2003  Comments added (R Saunders)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       02/12/2009  Introduced RTIASI emissivity model to be used with
!                        principal components. Marco Matricardi. ECMWF
!  1.4       18/10/2010  Always compute emissivity for PC-RTTOV regardless
!                        of calcemis (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_coef_pccomp, &
       & profile_Type
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : surftype_sea
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: addpc
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:)
  TYPE(profile_Type     ), INTENT(INOUT) :: profiles_ad(:)
  TYPE(profile_Type     ), INTENT(IN)    :: profiles   (:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis     (size(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity_ad(size(chanprof))
!INTF_END
!local variables:
! Type(profile_Type),  Pointer   :: prof
  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL   (KIND=jprb) :: windsp
  REAL   (KIND=jprb) :: windsp_ad
  REAL   (KIND=jprb) :: aems(size(chanprof))
  REAL   (KIND=jprb) :: bems(size(chanprof))
  REAL   (KIND=jprb) :: cems(size(chanprof))
  REAL   (KIND=jprb) :: expf
  REAL   (KIND=jprb) :: fac
  REAL   (KIND=jprb) :: aems_ad(size(chanprof))
  REAL   (KIND=jprb) :: bems_ad(size(chanprof))
  REAL   (KIND=jprb) :: cems_ad(size(chanprof))
  INTEGER(KIND=jpim) :: nchannels              ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! Loop on all channels
  nchannels = size(chanprof)
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_AD', 0_jpim, ZHOOK_HANDLE)
  aems_ad(:) = 0._jprb
  bems_ad(:) = 0._jprb
  cems_ad(:) = 0._jprb
  windsp_ad  = 0._jprb
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    iprof = chanprof(j)%prof

    IF (addpc) THEN
    
      IF (profiles(iprof)%skin%surftype == surftype_sea) THEN
        !------------------------------------------------------------------
        !1. Over sea, emissivity is computed using RTIASI model
        !------------------------------------------------------------------
        windsp     = sqrt(profiles(iprof)%s2m%u ** 2 + profiles(iprof)%s2m%v ** 2)
        aems(j)    =      &
          & coef_pccomp%emiss_c1(chan) + coef_pccomp%emiss_c2(chan) * windsp + coef_pccomp%emiss_c3(chan) * windsp ** 2
        bems(j)    =      &
          & coef_pccomp%emiss_c4(chan) + coef_pccomp%emiss_c5(chan) * windsp + coef_pccomp%emiss_c6(chan) * windsp ** 2
        cems(j)    = coef_pccomp%emiss_c7(chan) + coef_pccomp%emiss_c8(chan) * windsp
        expf       = exp(                                                                                                 &
          & ((coef_pccomp%emiss_c9(chan) - 60.) ** 2_jpim - (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) /  &
          & cems(j))
        fac        =  -                                                                                                   &
          & ((coef_pccomp%emiss_c9(chan) - 60.) ** 2_jpim - (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) /  &
          & cems(j) ** 2
        aems_ad(j) = aems_ad(j) + emissivity_ad(J)
        bems_ad(j) = bems_ad(j) + emissivity_ad(J) * expf
        cems_ad(j) = cems_ad(j) + bems(j) * expf * emissivity_ad(J) * fac
        aems_ad(j) = aems_ad(j) - emissivity_ad(J) * expf
        cems_ad(j) = cems_ad(j) - aems(j) * expf * emissivity_ad(J) * fac
        windsp_ad  = windsp_ad + coef_pccomp%emiss_c8(chan) * Cems_ad(j)
        windsp_ad  = windsp_ad + (coef_pccomp%emiss_c5(chan) + 2 * coef_pccomp%emiss_c6(chan) * windsp) * Bems_ad(j)
        windsp_ad  = windsp_ad + (coef_pccomp%emiss_c2(chan) + 2 * coef_pccomp%emiss_c3(chan) * windsp) * aems_ad(j)
        IF (profiles(iprof)%s2m%u /= 0. .AND. profiles(iprof)%s2m%v /= 0.) THEN
          profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + profiles(iprof)%s2m%u * windsp_AD / windsp
          profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + profiles(iprof)%s2m%v * windsp_AD / windsp
          windsp_ad                = 0.
        ELSE IF (profiles(iprof)%s2m%u == 0. .AND. profiles(iprof)%s2m%v == 0.) THEN
          profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + windsp_ad / sqrt(2.)
          profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + windsp_AD / sqrt(2.)
          windsp_ad                = 0.
        ELSE IF (profiles(iprof)%s2m%u == 0. .AND. profiles(iprof)%s2m%v /= 0.) THEN
          profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + WINDSP_ad
          windsp_ad                = 0
        ELSE IF (profiles(iprof)%s2m%u /= 0. .AND. profiles(iprof)%s2m%v == 0.) THEN
          profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + windsp_ad
          windsp_ad                = 0
        ENDIF
      ENDIF
    
    ELSE ! Not PC-RTTOV
    
      IF (.NOT. calcemis(j)) CYCLE
    
    END IF
    
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir_ad
