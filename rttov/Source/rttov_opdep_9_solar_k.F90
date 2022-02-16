SUBROUTINE rttov_opdep_9_solar_k( &
            & nlayers,      &
            & chanprof,     &
            & profiles,     &
            & sun,          &
            & predictors,   &
            & predictors_k, &
            & coef,         &
            & opdp_path,    &
            & opdp_path_k,  &
            & opdpsun_ref)
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
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_k.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  2.0    01/07/2006  Marco Matricardi (ECMWF):
!            --       The contribution of aerosols and clouds
!                     has been added to the total transmission.
!  3.0   06/02/2007   removed polarisation R Saunders
!  3.1   30/06/2009   Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
!  4.0   15/09/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!  5.0   03/11/2009   Transmittances / optical depths on levels (A Geer)
!  5.1   05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  5.2   04/08/2010   Move addsolar check to calling routine (J Hocking)
!  5.3   14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                     should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Adjoint variables
! input transmission_k% tau_surf and transmission_k% tau_level set inside integrate_k
!
! output predictors_k initialised inside rttov_k (need input
!    intent for memory allocation in calling routine)
!
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
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  REAL(KIND=jprb)      , INTENT(IN)    :: opdpsun_ref(nlayers, size(chanprof))
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: opticaldepth_k(nlayers, size(chanprof))
  REAL   (KIND=jprb) :: chanx
  INTEGER(KIND=jpim) :: lev         , lay, chan, i, j, ii, nlevels
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: nchannels                                 ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels           = size(chanprof)
  nlevels             = nlayers + 1
!----------------------------------------
!2.Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
  opticaldepth_k(:,:) = 0._JPRB
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = nlevels, 2,  - 1
        lay = lev - 1
        opdp_path_k%sun_level(lev - 1, j) = opdp_path_k%sun_level(lev - 1, j) + opdp_path_k%sun_level(lev, j)
        opticaldepth_k(lay, j)            = opticaldepth_k(lay, j) + opdp_path_k%sun_level(lev, j)
        opdp_path_k%sun_level(lev, j)     = 0.0_JPRB
      ENDDO
      opdp_path_k%sun_level(1, j) = 0.0_jprb
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO i = 1, nlayers
        IF (opdpsun_ref(i, j) > 0.0_JPRB) THEN
          opticaldepth_k(i, j) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!1.8 add CH4
!---------------------------------------------------------------------------------------
  IF (coef%nch4 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_k%ch4_sun(:, lay, j) =      &
            & predictors_k%ch4_sun(:, lay, j) + coef%ch4(lay, chan, :) * opticaldepth_k(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.7 add CO
!---------------------------------------------------------------------------------------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_k%co_sun(1:12, lay, j) =      &
            & predictors_k%co_sun(1:12, lay, j) + coef%co(lay, chan, 1:12) * opticaldepth_k(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.6 add N2O
!---------------------------------------------------------------------------------------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o_sun(12, lay, j) =      &
              & predictors_k%n2o_sun(12, lay, j) + coef%n2o(lay, chan, 14) * opticaldepth_k(lay, j)
            predictors_k%n2o_sun(11, lay, j) =      &
              & predictors_k%n2o_sun(11, lay, j) + coef%n2o(lay, chan, 13) * opticaldepth_k(lay, j)
            IF (coef%nco > 0) THEN
              predictors_k%co_sun(8, lay, j) = predictors_k%co_sun(8, lay, j) -                              &
                & coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * predictors%co_sun(1, lay, prof) ** 3 /  &
                & predictors%co_sun(8, lay, prof) ** 2
              predictors_k%co_sun(1, lay, j) = predictors_k%co_sun(1, lay, j) +                                  &
                & coef%n2o(lay, chan, 12) * opticaldepth_k(lay, j) * 3 * predictors%co_sun(1, lay, prof) ** 2 /  &
                & predictors%co_sun(8, lay, prof)
              predictors_k%co_sun(1, lay, j) =      &
                & predictors_k%co_sun(1, lay, j) + opticaldepth_k(lay, j) * coef%n2o(lay, chan, 11)
            ENDIF
            predictors_k%n2o_sun(1:10, lay, j) =      &
              & predictors_k%n2o_sun(1:10, lay, j) + coef%n2o(lay, chan, 1:10) * opticaldepth_k(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_k%n2o_sun(1:12, lay, j) =      &
              & predictors_k%n2o_sun(1:12, lay, j) + coef%n2o(lay, chan, 1:12) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.5 add CO2
!---------------------------------------------------------------------------------------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = nlayers, 1,  - 1
            predictors_k%co2_sun(1:14, lay, j) =      &
              & predictors_k%co2_sun(1:14, lay, j) + coef%co2(lay, chan, 1:14) * opticaldepth_k(lay, j)
          ENDDO
        ELSE
          DO lay = nlayers, 1,  - 1
            predictors_k%co2_sun(14, lay, j)  =      &
              & predictors_k%co2_sun(14, lay, j) + coef%co2(lay, chan, 13) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(13, lay, j)  =      &
              & predictors_k%co2_sun(13, lay, j) + coef%co2(lay, chan, 12) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(12, lay, j)  =      &
              & predictors_k%co2_sun(12, lay, j) + coef%co2(lay, chan, 11) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(11, lay, j)  =      &
              & predictors_k%co2_sun(11, lay, j) + coef%co2(lay, chan, 10) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(10, lay, j)  =      &
              & predictors_k%co2_sun(10, lay, j) + coef%co2(lay, chan, 9) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(9, lay, j)   =      &
              & predictors_k%co2_sun(9, lay, j) + coef%co2(lay, chan, 8) * opticaldepth_k(lay, j)
            predictors_k%co2_sun(1:7, lay, j) =      &
              & predictors_k%co2_sun(1:7, lay, j) + coef%co2(lay, chan, 1:7) * opticaldepth_k(lay, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.4 add Water Vapour Continuum
!---------------------------------------------------------------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          predictors_k%wvcont_sun(:, lay, j) =      &
            & predictors_k%wvcont_sun(:, lay, j) + coef%wvcont(lay, chan, :) * opticaldepth_k(lay, j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.3 add ozone
!---------------------------------------------------------------------------------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      IF (sun(j)) THEN
        DO lay = nlayers, 1,  - 1
          DO ii = 1, 13
            predictors_k%ozone_sun(ii, lay, j) =      &
              & predictors_k%ozone_sun(ii, lay, j) + coef%ozone(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!1.2 add water vapour
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    chanx = coef%ff_cwn(chan)
    prof  = chanprof(j)%prof
    IF (sun(j)) THEN
      IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = nlayers, 1,  - 1
          IF (coef%nco > 0) THEN
            predictors_k%co_sun(1, lay, j) =      &
              & predictors_k%co_sun(1, lay, j) + coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j)
          ENDIF
          IF (coef%nco2 > 0) THEN
            predictors_k%co2_sun(1, lay, j) =      &
              & predictors_k%co2_sun(1, lay, j) + coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j)
          ENDIF
          DO ii = 1, 13
            predictors_k%watervapour_sun(ii, lay, j) =      &
              & predictors_k%watervapour_sun(ii, lay, j) + coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2360_00) THEN
        DO lay = nlayers, 1,  - 1
          predictors_k%watervapour_sun(13, lay, j) =      &
            & predictors_k%watervapour_sun(13, lay, j) + coef%watervapour(lay, chan, 13) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(12, lay, j) =      &
            & predictors_k%watervapour_sun(12, lay, j) + coef%watervapour(lay, chan, 12) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(11, lay, j) =      &
            & predictors_k%watervapour_sun(11, lay, j) + coef%watervapour(lay, chan, 11) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(10, lay, j) =      &
            & predictors_k%watervapour_sun(10, lay, j) + coef%watervapour(lay, chan, 10) * opticaldepth_k(lay, j)
          DO ii = 1, 9
            predictors_k%watervapour_sun(ii, lay, j) =      &
              & predictors_k%watervapour_sun(ii, lay, j) + coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ELSE
        DO lay = nlayers, 1,  - 1
          predictors_k%ch4_sun(3, lay, j)          = predictors_k%ch4_sun(3, lay, j) +      &
            & coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j) * predictors%ch4_sun(1, lay, prof) ** 0.25
          predictors_k%ch4_sun(1, lay, j)          = predictors_k%ch4_sun(1, lay, j) +      &
            & coef%watervapour(lay, chan, 15) * opticaldepth_k(lay, j) * 0.25 *             &
            & predictors%ch4_sun(1, lay, prof) ** ( - 0.75) * predictors%ch4_sun(3, lay, prof)
          predictors_k%ch4_sun(1, lay, j)          = predictors_k%ch4_sun(1, lay, j) +      &
            & coef%watervapour(lay, chan, 14) * opticaldepth_k(lay, j) * 1.25 * predictors%ch4_sun(1, lay, prof) ** 0.25
          predictors_k%watervapour_sun(15, lay, j) =      &
            & predictors_k%watervapour_sun(15, lay, j) + coef%watervapour(lay, chan, 13) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(12, lay, j) =      &
            & predictors_k%watervapour_sun(12, lay, j) + coef%watervapour(lay, chan, 12) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(11, lay, j) =      &
            & predictors_k%watervapour_sun(11, lay, j) + coef%watervapour(lay, chan, 11) * opticaldepth_k(lay, j)
          predictors_k%watervapour_sun(14, lay, j) =      &
            & predictors_k%watervapour_sun(14, lay, j) + coef%watervapour(lay, chan, 10) * opticaldepth_k(lay, j)
          DO ii = 1, 9
            predictors_k%watervapour_sun(ii, lay, j) =      &
              & predictors_k%watervapour_sun(ii, lay, j) + coef%watervapour(lay, chan, ii) * opticaldepth_k(lay, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!1.1 start with mixed gases
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    IF (sun(j)) THEN
      DO lay = nlayers, 1,  - 1
        predictors_k%mixedgas_sun(10, lay, j) =      &
          & predictors_k%mixedgas_sun(10, lay, j) + coef%mixedgas(lay, chan, 10) * opticaldepth_k(lay, j)
        predictors_k%mixedgas_sun(9, lay, j)  =      &
          & predictors_k%mixedgas_sun(9, lay, j) + coef%mixedgas(lay, chan, 9) * opticaldepth_k(lay, j)
        DO ii = 2, 8
          predictors_k%mixedgas_sun(ii, lay, j) =      &
            & predictors_k%mixedgas_sun(ii, lay, j) + coef%mixedgas(lay, chan, ii) * opticaldepth_k(lay, j)
        ENDDO
        predictors_k%mixedgas_sun(1, lay, j) =      &
          & predictors_k%mixedgas_sun(1, lay, j) + coef%mixedgas(lay, chan, 1) * opticaldepth_k(lay, j)
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_solar_k
