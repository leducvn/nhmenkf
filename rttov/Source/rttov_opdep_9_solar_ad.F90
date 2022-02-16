SUBROUTINE rttov_opdep_9_solar_ad( &
            & nlayers,       &
            & chanprof,      &
            & profiles,      &
            & sun,           &
            & predictors,    &
            & predictors_ad, &
            & coef,          &
            & opdp_path,     &
            & opdp_path_ad,  &
            & opdpsun_ref)
!
! Description:
! Adjoint of rttov_transmit_ad
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPAD and RTTAUAD from previous versions of RTTOV
! Only one profile per call
!
! Adjoint variables
! input transmission_ad % tau_surf and transmission_ad % tau_level
! set inside integrate_ad
!
! output predictors_ad initialised inside rttov_ad (need input
!    intent for memory allocation in calling routine)
!
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_ad.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    06/02/2007  Removed polarisation index R Saunders
!  1.2    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.3    30/06/2009  Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
!  1.4    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.5    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.6    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  1.7    14/12/2010  Use traj0_sta%sun array to flag channels for which solar calculations
!                     should be performed (J Hocking)
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
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov9_wv2000_00, &
       & rttov9_wv2295_25, &
       & rttov9_wv2360_00, &
       & max_sol_zen
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type   ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: sun(size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_ad
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  REAL(KIND=jprb)      , INTENT(IN)    :: opdpsun_ref(nlayers, size(chanprof))
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: opticaldepth_ad(nlayers, size(chanprof))
  REAL   (KIND=jprb) :: chanx
  INTEGER(KIND=jpim) :: lev         , lay, chan, j, nlevels
  INTEGER(KIND=jpim) :: prof
! cloud liquid water local variables
  INTEGER(KIND=jpim) :: II
  INTEGER(KIND=jpim) :: nchannels                               ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels            = size(chanprof)
  nlevels              = nlayers + 1
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  opticaldepth_ad(:,:) = 0._JPRB
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = nlevels, 2,  - 1
        lay = lev - 1
        opdp_path_ad%sun_level(lev - 1, j) = opdp_path_ad%sun_level(lev - 1, j) + opdp_path_ad%sun_level(lev, j)
        opticaldepth_ad(lay, j)            = opticaldepth_ad(lay, j) + opdp_path_ad%sun_level(lev, j)
        opdp_path_ad%sun_level(lev, j)     = 0.0_JPRB
      ENDDO
      opdp_path_ad%sun_level(1, j) = 0.0_jprb
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lay = 1, nlayers
        IF (opdpsun_ref(lay, j) > 0.0_JPRB) THEN
          opticaldepth_ad(lay, j) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDIF
  ENDDO
!-----------
!1.8 add CH4
!-----------
  IF (coef%nch4 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          DO ii = 1, coef%nch4
            predictors_ad%ch4_sun(ii, lay, prof) =      &
              & predictors_ad%ch4_sun(ii, lay, prof) + coef%ch4(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.7 add CO
!-----------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_ad%co_sun(1:12, lay, prof) =      &
            & predictors_ad%co_sun(1:12, lay, prof) + coef%co(lay, chan, 1:12) * opticaldepth_ad(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.6 add N2O
!-----------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o_sun(12, lay, prof) =      &
              & predictors_ad%n2o_sun(12, lay, prof) + coef%n2o(lay, chan, 14) * opticaldepth_ad(lay, j)
            predictors_ad%n2o_sun(11, lay, prof) =      &
              & predictors_ad%n2o_sun(11, lay, prof) + coef%n2o(lay, chan, 13) * opticaldepth_ad(lay, j)
            IF (coef%nco > 0) THEN
              predictors_ad%co_sun(8, lay, prof) = predictors_ad%co_sun(8, lay, prof) -                       &
                & coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * predictors%co_sun(1, lay, prof) ** 3 /  &
                & predictors%co_sun(8, lay, prof) ** 2
              predictors_ad%co_sun(1, lay, prof) = predictors_ad%co_sun(1, lay, prof) +                           &
                & coef%n2o(lay, chan, 12) * opticaldepth_ad(lay, j) * 3 * predictors%co_sun(1, lay, prof) ** 2 /  &
                & predictors%co_sun(8, lay, prof)
              predictors_ad%co_sun(1, lay, prof) =      &
                & predictors_ad%co_sun(1, lay, prof) + opticaldepth_ad(lay, j) * coef%n2o(lay, chan, 11)
            ENDIF
            predictors_ad%n2o_sun(1:10, lay, prof) =      &
              & predictors_ad%n2o_sun(1:10, lay, prof) + coef%n2o(lay, chan, 1:10) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_ad%n2o_sun(1:12, lay, prof) =      &
              & predictors_ad%n2o_sun(1:12, lay, prof) + coef%n2o(lay, chan, 1:12) * opticaldepth_ad(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      chanx = coef%ff_cwn(chan)
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2_sun(1:14, lay, prof) =      &
              & predictors_ad%co2_sun(1:14, lay, prof) + coef%co2(lay, chan, 1:14) * opticaldepth_ad(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_ad%co2_sun(14, lay, prof)  =      &
              & predictors_ad%co2_sun(14, lay, prof) + coef%co2(lay, chan, 13) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(13, lay, prof)  =      &
              & predictors_ad%co2_sun(13, lay, prof) + coef%co2(lay, chan, 12) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(12, lay, prof)  =      &
              & predictors_ad%co2_sun(12, lay, prof) + coef%co2(lay, chan, 11) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(11, lay, prof)  =      &
              & predictors_ad%co2_sun(11, lay, prof) + coef%co2(lay, chan, 10) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(10, lay, prof)  =      &
              & predictors_ad%co2_sun(10, lay, prof) + coef%co2(lay, chan, 9) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(9, lay, prof)   =      &
              & predictors_ad%co2_sun(9, lay, prof) + coef%co2(lay, chan, 8) * opticaldepth_ad(lay, j)
            predictors_ad%co2_sun(1:7, lay, prof) =      &
              & predictors_ad%co2_sun(1:7, lay, prof) + coef%co2(lay, chan, 1:7) * opticaldepth_ad(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_ad%wvcont_sun(:, lay, prof) =      &
            & predictors_ad%wvcont_sun(:, lay, prof) + coef%wvcont(lay, chan, :) * opticaldepth_ad(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          DO ii = 1, 13
            predictors_ad%ozone_sun(ii, lay, prof) =      &
              & predictors_ad%ozone_sun(ii, lay, prof) + coef%ozone(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    chanx = coef%ff_cwn(chan)
    prof  = chanprof(j)%prof
    IF (sun(j)) THEN
      IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_ad%co_sun(1, lay, prof) =      &
              & predictors_ad%co_sun(1, lay, prof) + coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_ad%co2_sun(1, lay, prof) =      &
              & predictors_ad%co2_sun(1, lay, prof) + coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j)
          ENDIF
          DO ii = 1, 13
            predictors_ad%watervapour_sun(ii, lay, prof) =      &
              & predictors_ad%watervapour_sun(ii, lay, prof) + coef%watervapour(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2360_00) THEN
        DO lay = nlayers, 1,  - 1
          predictors_ad%watervapour_sun(13, lay, prof) =      &
            & predictors_ad%watervapour_sun(13, lay, prof) + coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(12, lay, prof) =      &
            & predictors_ad%watervapour_sun(12, lay, prof) + coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(11, lay, prof) =      &
            & predictors_ad%watervapour_sun(11, lay, prof) + coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(10, lay, prof) =      &
            & predictors_ad%watervapour_sun(10, lay, prof) + coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          DO ii = 1, 9
            predictors_ad%watervapour_sun(ii, lay, prof) =      &
              & predictors_ad%watervapour_sun(ii, lay, prof) + coef%watervapour(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ELSE
        DO lay = nlayers, 1,  - 1
          IF (coef%nch4 > 0) THEN
            predictors_ad%ch4_sun(3, lay, prof) = predictors_ad%ch4_sun(3, lay, prof) +      &
              & coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j) * predictors%ch4_sun(1, lay, prof) ** 0.25
            predictors_ad%ch4_sun(1, lay, prof) = predictors_ad%ch4_sun(1, lay, prof) +      &
              & coef%watervapour(lay, chan, 15) * opticaldepth_ad(lay, j) * 0.25 *           &
              & predictors%ch4_sun(1, lay, prof) ** ( - 0.75) * predictors%ch4_sun(3, lay, prof)
            predictors_ad%ch4_sun(1, lay, prof) = predictors_ad%ch4_sun(1, lay, prof) +      &
              & coef%watervapour(lay, chan, 14) * opticaldepth_ad(lay, j) * 1.25 *           &
              & predictors%ch4_sun(1, lay, prof) ** 0.25
          ENDIF
          predictors_ad%watervapour_sun(15, lay, prof) =      &
            & predictors_ad%watervapour_sun(15, lay, prof) + coef%watervapour(lay, chan, 13) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(12, lay, prof) =      &
            & predictors_ad%watervapour_sun(12, lay, prof) + coef%watervapour(lay, chan, 12) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(11, lay, prof) =      &
            & predictors_ad%watervapour_sun(11, lay, prof) + coef%watervapour(lay, chan, 11) * opticaldepth_ad(lay, j)
          predictors_ad%watervapour_sun(14, lay, prof) =      &
            & predictors_ad%watervapour_sun(14, lay, prof) + coef%watervapour(lay, chan, 10) * opticaldepth_ad(lay, j)
          DO ii = 1, 9
            predictors_ad%watervapour_sun(ii, lay, prof) =      &
              & predictors_ad%watervapour_sun(ii, lay, prof) + coef%watervapour(lay, chan, ii) * opticaldepth_ad(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    prof  = chanprof(j)%prof
    IF (sun(j)) THEN
      DO lay = nlayers, 1,  - 1
        predictors_ad%mixedgas_sun(10, lay, prof) =      &
          & predictors_ad%mixedgas_sun(10, lay, prof) + coef%mixedgas(lay, chan, 10) * opticaldepth_ad(lay, j)
        predictors_ad%mixedgas_sun(9, lay, prof)  =      &
          & predictors_ad%mixedgas_sun(9, lay, prof) + coef%mixedgas(lay, chan, 9) * opticaldepth_ad(lay, j)
        DO ii = 2, 8
          predictors_ad%mixedgas_sun(ii, lay, prof) =      &
            & predictors_ad%mixedgas_sun(ii, lay, prof) + coef%mixedgas(lay, chan, ii) * opticaldepth_ad(lay, j)
        ENDDO
        predictors_ad%mixedgas_sun(1, lay, prof) =      &
          & predictors_ad%mixedgas_sun(1, lay, prof) + coef%mixedgas(lay, chan, 1) * opticaldepth_ad(lay, j)
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_solar_ad
