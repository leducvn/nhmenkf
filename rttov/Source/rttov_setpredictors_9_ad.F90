!
SUBROUTINE rttov_setpredictors_9_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef_pccomp,   &
            & coef,          &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
! Description
! RTTOV-8 Model
! AD of rttov_setpredictors_8
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
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
! see RTTOV7 science and validation report pages 18/19
! variable names are close to the documentation
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0   01/06/2005  Marco Matricardi (ECMWF):
!           --       New routine based on rttov_setpredictors_ad.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.2   02/12/2009  Introduced principal component capability. Pathsat, Pathsun and
!                    related quantities are now layer arrays (Marco Matricardi).
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
       & rttov_coef,        &
       & rttov_options,     &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type,     &
       & predictors_Type,   &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE parkind1, ONLY : JPRB
  USE parkind1, ONLY : JPRB
  USE parkind1, ONLY : JPRB
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(profile_Type     ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type     ), INTENT(INOUT) :: prof_ad(size(prof))
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type  ), INTENT(IN)    :: predictors
  TYPE(predictors_Type  ), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type  ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type  ), INTENT(INOUT) :: raytracing_ad
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer, iprof
! user profile
  REAL   (KIND=Jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2o  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4  (prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tro  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2r (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2or (prof(1)%nlayers)
  REAL   (KIND=jprb) :: cor  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr  (prof(1)%nlayers)
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuw  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuwr (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2w (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2ow (prof(1)%nlayers)
  REAL   (KIND=jprb) :: cow  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4w (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2owr(prof(1)%nlayers)
  REAL   (KIND=jprb) :: cowr (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4wr(prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: deltac     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ww    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_wwr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ow    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_twr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_tuw   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_co2w  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ch4w  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_ch4wr (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_n2ow  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_n2owr (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_cow   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sum2_cowr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sq      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_4       (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wrwr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr     (prof(1)%nlayers)
! AD variables
  REAL   (KIND=Jprb) :: t_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o_ad       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2o_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tro_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2or_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cor_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4r_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuw_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuwr_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dt_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tw_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow_ad      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2ow_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cow_ad     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4w_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2owr_ad   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cowr_ad    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4wr_ad   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_or_ad  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr_ad  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr_ad(prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles, nlayers          ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! Recompute Direct variables
!-------------------------------------------------------------------------------
!1) Profile layer quantities
!-------------------------------------------------------------------------------
!-Temperature--------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nlayers = prof(1)%nlayers

    IF (opts%addpc) THEN
!-CO2--------------------------------------------------------------------------

      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
          co2(layer) = (coef_pccomp%co2_pc_ref(level - 1) + coef_pccomp%co2_pc_ref(level)) * 0.5_JPRB
        ENDIF
      ENDDO

!-N2O--------------------------------------------------------------------------
      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
          n2o(layer) = (coef_pccomp%n2o_pc_ref(level - 1) + coef_pccomp%n2o_pc_ref(level)) * 0.5_JPRB
        ENDIF
      ENDDO

!-CO---------------------------------------------------------------------------
      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        IF (opts%co_Data .AND. coef%nco > 0) THEN
          co(layer) = (coef_pccomp%co_pc_ref(level - 1) + coef_pccomp%co_pc_ref(level)) * 0.5_JPRB
        ENDIF
      ENDDO

!-CH4--------------------------------------------------------------------------
      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
          ch4(layer) = (coef_pccomp%ch4_pc_ref(level - 1) + coef_pccomp%ch4_pc_ref(level)) * 0.5_JPRB
        ENDIF
      ENDDO

    ENDIF
! addpc


  DO iprof = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
!-Temperature  ----------------------------------------------------------------
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) * 0.5_JPRB
!-H2O--------------------------------------------------------------------------
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) * 0.5_JPRB
!-O3---------------------------------------------------------------------------
      IF (opts%ozone_Data .AND. coef%nozone > 0)     &
        &  o(layer) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) * 0.5_JPRB
    ENDDO

    IF (.not. opts%addpc) THEN

      DO layer = 1, prof(1)%nlayers
        level    = layer + 1

!-CO2--------------------------------------------------------------------------

        IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
          co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) * 0.5_JPRB
        ENDIF

!-N2O--------------------------------------------------------------------------

        IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
          n2o(layer) = (prof(iprof)%n2o(level - 1) + prof(iprof)%n2o(level)) * 0.5_JPRB
        ENDIF

!-CO---------------------------------------------------------------------------

        IF (opts%co_Data .AND. coef%nco > 0) THEN
          co(layer) = (prof(iprof)%co(level - 1) + prof(iprof)%co(level)) * 0.5_JPRB
        ENDIF

!-CH4--------------------------------------------------------------------------

        IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
          ch4(layer) = (prof(iprof)%ch4(level - 1) + prof(iprof)%ch4(level)) * 0.5_JPRB
        ENDIF

      ENDDO
! layers

    ENDIF
! not addpc


!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:)      = t(:) - coef%tstar(:)
    dtabs(:)   = Abs(dt(:))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:)      = t(:) / coef%tstar(:)
    tr_sq(:)   = tr(:) * tr(:)
    tr_4(:)    = tr_sq(:) * tr_sq(:)
    wr(:)      = w(:) / coef%wstar(:)

    IF (coef%nozone > 0) THEN
      dto(:) = t(:) - coef%to3star(:)
      tro(:) = t(:) / coef%to3star(:)
      IF (opts%ozone_Data) THEN
        or(:)  = o(:) / coef%ostar(:)
      ELSE
        or(:)  = 1._JPRB
      ENDIF
    ENDIF

!-CO2----------------------------------------------------------------------------

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:) = co2(:) / coef%co2star(:)
    ELSE
      co2r(:) = 1._JPRB
    ENDIF

!-N2O----------------------------------------------------------------------------

    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:) = n2o(:) / coef%n2ostar(:)
    ELSE
      n2or(:) = 1._JPRB
    ENDIF

!-CO----------------------------------------------------------------------------

    IF (opts%co_Data .AND. coef%nco > 0) THEN
      cor(:) = co(:) / coef%costar(:)
    ELSE
      cor(:) = 1._JPRB
    ENDIF

!-CH4---------------------------------------------------------------------------

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:) = ch4(:) / coef%ch4star(:)
    ELSE
      ch4r(:) = 1._JPRB
    ENDIF

!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr
!--------------------------------------------------------------------
    tw(1) = 0._jprb

    DO layer = 2, prof(1)%nlayers
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1            = sum1 + t(layer)
      sum2            = sum2 + coef%tstar(layer)
      sum2_tuw(layer) = sum2
      tuw(layer)      = sum1 / sum2
    ENDDO

    tuwr(1)                 = coef%dpp(0) * t(1) / (coef%dpp(0) * coef%tstar(1))
    tuwr(2:prof(1)%nlayers) = tuw(2:prof(1)%nlayers)

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1            = sum1 + coef%dpp(layer - 1) * t(layer)
        sum2            = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        sum2_twr(layer) = sum2
        twr(layer)      = sum1 / sum2
      ENDDO
    ENDIF

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1            = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2            = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      sum2_wwr(layer) = sum2
      wwr(layer)      = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1           = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2           = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(layer) = sum2
      ww(layer)      = sum1 / sum2
    ENDDO


    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1           = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2           = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        sum2_ow(layer) = sum2
        ow(layer)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_ow(:) = 0._JPRB
      ow(:)      = 1._JPRB
    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1             = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2             = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        sum2_co2w(layer) = sum2
        co2w(layer)      = sum1 / sum2
      ENDDO

    ELSE
      sum2_co2w(:) = 0._JPRB
      co2w(:)      = 1._JPRB
    ENDIF

!-N2O---------------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * n2o(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          sum2_n2ow(layer) = sum2
          n2ow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * n2o(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer) = sum2
          n2owr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_n2ow(:) = 0._JPRB
        n2ow(:)  = 1._JPRB

        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer) = sum2
          n2owr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF

!-CO---------------------------------------------------------------------------

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1       = sum1 + coef%dpp(layer - 1) * co(layer)
          sum2       = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          sum2_cow(layer) = sum2
          cow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * co(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer) = sum2
          cowr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_cow(:) = 0._JPRB
        cow(:)  = 1._JPRB
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer) = sum2
          cowr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF

!-CH4---------------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * ch4(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          sum2_ch4w(layer) = sum2
          ch4w(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * ch4(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer) = sum2
          ch4wr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        sum2_ch4w(:) = 0._JPRB
        ch4w(:)  = 1._JPRB

        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer) = sum2
          ch4wr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF

!-------------------------------------------------------------------------
! Adjoint code
!-------------------------------------------------------------------------
! DAR - 1:nlayers is seemingly necessary to stop ifort from optimising these lines out at -O3

    w_ad(1:nlayers)        = 0._JPRB
    wr_ad(1:nlayers)       = 0._JPRB
    ww_ad(1:nlayers)       = 0._JPRB
    wwr_ad(1:nlayers)      = 0._JPRB
    sec_wr_ad(1:nlayers)   = 0._JPRB
    sec_wrwr_ad(1:nlayers) = 0._JPRB
    dt_ad(1:nlayers)       = 0._JPRB
    dto_ad(1:nlayers)      = 0._JPRB
    t_ad(1:nlayers)        = 0._JPRB
    tr_ad(1:nlayers)       = 0._JPRB
    tro_ad(1:nlayers)      = 0._JPRB
    tw_ad(1:nlayers)       = 0._JPRB
    twr_ad(1:nlayers)      = 0._JPRB
    tuw_ad(1:nlayers)      = 0._JPRB
    tuwr_ad(1:nlayers)     = 0._JPRB
    ch4_ad(1:nlayers)      = 0._JPRB
    ch4r_ad(1:nlayers)     = 0._JPRB
    ch4w_ad(1:nlayers)     = 0._JPRB
    ch4wr_ad(1:nlayers)    = 0._JPRB
    co_ad(1:nlayers)       = 0._JPRB
    cor_ad(1:nlayers)      = 0._JPRB
    cow_ad(1:nlayers)      = 0._JPRB
    cowr_ad(1:nlayers)     = 0._JPRB
    n2o_ad(1:nlayers)      = 0._JPRB
    n2or_ad(1:nlayers)     = 0._JPRB
    n2ow_ad(1:nlayers)     = 0._JPRB
    n2owr_ad(1:nlayers)    = 0._JPRB
    co2_ad(1:nlayers)      = 0._JPRB
    co2r_ad(1:nlayers)     = 0._JPRB
    co2w_ad(1:nlayers)     = 0._JPRB
!5.9 ch4            transmittance based on RTIASI
!-------------------------------------------------
!

    DO layer = 1, prof(1)%nlayers
      level = layer + 1

      IF (coef%nch4 > 0) THEN
        ch4r_ad(layer)                      =      &
          & ch4r_ad(layer) + predictors_ad%ch4(11, layer, iprof) * 1.5 * predictors%ch4(2, layer, iprof) / ch4w(layer)
        ch4w_ad(layer)                      = ch4w_ad(layer) -                            &
          & predictors_ad%ch4(11, layer, iprof) * predictors%ch4(2, layer, iprof) ** 3 /  &
          & (raytracing%pathsat(layer, iprof) * ch4w(layer) ** 2)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ch4(11, layer, iprof) * 0.5 * ch4r(layer) ** 1.5 /             &
          & (raytracing%pathsat(layer, iprof) ** 0.5 * ch4w(layer))
        ch4w_ad(layer)                      =      &
          & ch4w_ad(layer) + raytracing%pathsat(layer, iprof) * predictors_ad%ch4(10, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ch4(10, layer, iprof) * ch4w(layer)
        ch4w_ad(layer)                      = ch4w_ad(layer) +      &
          & predictors_ad%ch4(9, layer, iprof) * predictors%ch4(10, layer, iprof) * 2 * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ch4(9, layer, iprof) * 2 * predictors%ch4(10, layer, iprof) * ch4w(layer)
        ch4wr_ad(layer)                     = ch4wr_ad(layer) + predictors_ad%ch4(8, layer, iprof)
        ch4wr_ad(layer)                     =      &
          & ch4wr_ad(layer) + predictors_ad%ch4(7, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ch4(7, layer, iprof) * ch4wr(layer)
        ch4r_ad(layer)                      = ch4r_ad(layer) +                              &
          & predictors_ad%ch4(6, layer, iprof) * raytracing%pathsat(layer, iprof) * 0.25 /  &
          & predictors%ch4(6, layer, iprof) ** 3
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ch4(6, layer, iprof) * ch4r(layer) * 0.25 / predictors%ch4(6, layer, iprof) ** 3
        dt_ad(layer)                        = dt_ad(layer) + predictors_ad%ch4(5, layer, iprof) * ch4r(layer)
        ch4r_ad(layer)                      = ch4r_ad(layer) + predictors_ad%ch4(5, layer, iprof) * dt(layer)
        ch4r_ad(layer)                      = ch4r_ad(layer) +      &
          & predictors_ad%ch4(4, layer, iprof) * predictors%ch4(1, layer, iprof) * raytracing%pathsat(layer, iprof) * 2
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ch4(4, layer, iprof) * 2 * predictors%ch4(1, layer, iprof) * ch4r(layer)
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%ch4(3, layer, iprof) * predictors%ch4(1, layer, iprof)
        ch4r_ad(layer)                      =      &
          & ch4r_ad(layer) + predictors_ad%ch4(3, layer, iprof) * raytracing%pathsat(layer, iprof) * dt(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ch4(3, layer, iprof) * ch4r(layer) * dt(layer)
        ch4r_ad(layer)                      = ch4r_ad(layer) +      &
          & predictors_ad%ch4(2, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) / predictors%ch4(2, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ch4(2, layer, iprof) * 0.5 * ch4r(layer) / predictors%ch4(2, layer, iprof)
        ch4r_ad(layer)                      =      &
          & ch4r_ad(layer) + predictors_ad%ch4(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ch4(1, layer, iprof) * ch4r(layer)
      ENDIF

!5.8 co             transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nco > 0) THEN
        cowr_ad(layer)                      = cowr_ad(layer) +                              &
          & predictors_ad%co(13, layer, iprof) * 0.25 * raytracing%pathsat(layer, iprof) *  &
          & (raytracing%pathsat(layer, iprof) * cowr(layer)) ** ( - 0.75)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(13, layer, iprof) * 0.25 * cowr(layer) *                    &
          & (raytracing%pathsat(layer, iprof) * cowr(layer)) ** ( - 0.75)
        cowr_ad(layer)                      = cowr_ad(layer) +                             &
          & predictors_ad%co(12, layer, iprof) * raytracing%pathsat(layer, iprof) * 0.4 *  &
          & (raytracing%pathsat(layer, iprof) * cowr(layer)) ** ( - 0.6)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(12, layer, iprof) * cowr(layer) * 0.4 *                     &
          & (raytracing%pathsat(layer, iprof) * cowr(layer)) ** ( - 0.6)
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(11, layer, iprof) * 2 * predictors%co(1, layer, iprof) / cow(layer) ** 0.25
        cow_ad(layer)                       =      &
          & cow_ad(layer) - predictors_ad%co(11, layer, iprof) * 0.25 * predictors%co(11, layer, iprof) / cow(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(11, layer, iprof) * cor(layer) ** 2 / cow(layer) ** 0.25
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(10, layer, iprof) * 2 * predictors%co(1, layer, iprof) / cow(layer) ** 0.5
        cow_ad(layer)                       =      &
          & cow_ad(layer) - predictors_ad%co(10, layer, iprof) * 0.5 * predictors%co(10, layer, iprof) / cow(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(10, layer, iprof) * cor(layer) ** 2 / SQRT(cow(layer))
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(9, layer, iprof) * 1.5 * predictors%co(2, layer, iprof) / cow(layer)
        cow_ad(layer)                       = cow_ad(layer) -                          &
          & predictors_ad%co(9, layer, iprof) * predictors%co(2, layer, iprof) ** 3 /  &
          & (raytracing%pathsat(layer, iprof) * cow(layer) ** 2)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(9, layer, iprof) * 0.5 * cor(layer) ** 1.5 /                &
          & (raytracing%pathsat(layer, iprof) ** 0.5 * cow(layer))
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(8, layer, iprof) * 2 * predictors%co(1, layer, iprof) / cow(layer)
        cow_ad(layer)                       =      &
          & cow_ad(layer) - predictors_ad%co(8, layer, iprof) * predictors%co(8, layer, iprof) / cow(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(8, layer, iprof) * cor(layer) ** 2 / cow(layer)
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(7, layer, iprof) * dtabs(layer) * raytracing%pathsat(layer, iprof) * dt(layer)
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%co(7, layer, iprof) * dtabs(layer) * predictors%co(1, layer, iprof) * 2
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(7, layer, iprof) * dtabs(layer) * cor(layer) * dt(layer)
        cor_ad(layer)                       = cor_ad(layer) +                              &
          & predictors_ad%co(6, layer, iprof) * raytracing%pathsat(layer, iprof) * 0.25 /  &
          & predictors%co(6, layer, iprof) ** 3
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(6, layer, iprof) * cor(layer) * 0.25 / predictors%co(6, layer, iprof) ** 3
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%co(5, layer, iprof) * predictors%co(2, layer, iprof)
        cor_ad(layer)                       = cor_ad(layer) +                                         &
          & predictors_ad%co(5, layer, iprof) * 0.5 * dt(layer) * raytracing%pathsat(layer, iprof) /  &
          & predictors%co(2, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(5, layer, iprof) * 0.5 * dt(layer) * cor(layer) / predictors%co(2, layer, iprof)
        cor_ad(layer)                       = cor_ad(layer) +      &
          & predictors_ad%co(4, layer, iprof) * 2 * predictors%co(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(4, layer, iprof) * 2 * predictors%co(1, layer, iprof) * cor(layer)
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%co(3, layer, iprof) * predictors%co(1, layer, iprof)
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(3, layer, iprof) * dt(layer) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(3, layer, iprof) * cor(layer) * dt(layer)
        cor_ad(layer)                       = cor_ad(layer) +      &
          & predictors_ad%co(2, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) / predictors%co(2, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co(2, layer, iprof) * 0.5 * cor(layer) / predictors%co(2, layer, iprof)
        cor_ad(layer)                       =      &
          & cor_ad(layer) + predictors_ad%co(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co(1, layer, iprof) * cor(layer)
      ENDIF

!5.7 n2o            transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nn2o > 0) THEN
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%n2o(13, layer, iprof) * raytracing%pathsat(layer, iprof) ** 2 * n2owr(layer)
        n2owr_ad(layer)                     =      &
          & n2owr_ad(layer) + predictors_ad%n2o(13, layer, iprof) * raytracing%pathsat(layer, iprof) ** 2 * dt(layer)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(13, layer, iprof) * 2 * raytracing%pathsat(layer, iprof) * n2owr(layer) * dt(layer)
        n2owr_ad(layer)                     = n2owr_ad(layer) +                               &
          & predictors_ad%n2o(12, layer, iprof) * 3 * predictors%n2o(8, layer, iprof) ** 2 *  &
          & raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(12, layer, iprof) * 3 * predictors%n2o(8, layer, iprof) ** 2 * n2owr(layer)
        n2owr_ad(layer)                     = n2owr_ad(layer) +      &
          & predictors_ad%n2o(11, layer, iprof) * 2 * predictors%n2o(8, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(11, layer, iprof) * predictors%n2o(8, layer, iprof) * 2 * n2owr(layer)
        n2or_ad(layer)                      =      &
          & n2or_ad(layer) + predictors_ad%n2o(10, layer, iprof) * 1.5 * predictors%n2o(2, layer, iprof) / n2ow(layer)
        n2ow_ad(layer)                      = n2ow_ad(layer) -                            &
          & predictors_ad%n2o(10, layer, iprof) * predictors%n2o(2, layer, iprof) ** 3 /  &
          & (raytracing%pathsat(layer, iprof) * n2ow(layer) ** 2)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(10, layer, iprof) * 0.5 * n2or(layer) ** 1.5 /             &
          & (raytracing%pathsat(layer, iprof) ** 0.5 * n2ow(layer))
        n2owr_ad(layer)                     = n2owr_ad(layer) + predictors_ad%n2o(9, layer, iprof)
        n2owr_ad(layer)                     =      &
          & n2owr_ad(layer) + predictors_ad%n2o(8, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%n2o(8, layer, iprof) * n2owr(layer)
        n2ow_ad(layer)                      =      &
          & n2ow_ad(layer) + predictors_ad%n2o(7, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%n2o(7, layer, iprof) * n2ow(layer)
        n2or_ad(layer)                      = n2or_ad(layer) +                              &
          & predictors_ad%n2o(6, layer, iprof) * raytracing%pathsat(layer, iprof) * 0.25 /  &
          & predictors%n2o(6, layer, iprof) ** 3
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(6, layer, iprof) * n2or(layer) * 0.25 / predictors%n2o(6, layer, iprof) ** 3
        dt_ad(layer)                        = dt_ad(layer) + predictors_ad%n2o(5, layer, iprof) * n2or(layer)
        n2or_ad(layer)                      = n2or_ad(layer) + predictors_ad%n2o(5, layer, iprof) * dt(layer)
        n2or_ad(layer)                      = n2or_ad(layer) +      &
          & predictors_ad%n2o(4, layer, iprof) * 2 * predictors%n2o(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(4, layer, iprof) * 2 * predictors%n2o(1, layer, iprof) * n2or(layer)
        dt_ad(layer)                        =      &
          & dt_ad(layer) + predictors_ad%n2o(3, layer, iprof) * predictors%n2o(1, layer, iprof)
        n2or_ad(layer)                      =      &
          & n2or_ad(layer) + predictors_ad%n2o(3, layer, iprof) * raytracing%pathsat(layer, iprof) * dt(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%n2o(3, layer, iprof) * n2or(layer) * dt(layer)
        n2or_ad(layer)                      = n2or_ad(layer) +      &
          & predictors_ad%n2o(2, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) / predictors%n2o(2, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%n2o(2, layer, iprof) * 0.5 * n2or(layer) / predictors%n2o(2, layer, iprof)
        n2or_ad(layer)                      =      &
          & n2or_ad(layer) + predictors_ad%n2o(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%n2o(1, layer, iprof) * n2or(layer)
      ENDIF

!5.6 CO2
!-------

      IF (coef%nco2 > 0) THEN
!    co2r_ad(:)    = 0._JPRB
!    co2w_ad(:)    = 0._JPRB
!    twr_ad(:)     = 0._JPRB
!    co2_ad(:)     = 0._JPRB
        tr_ad(layer)                        =      &
          & tr_ad(layer) + predictors_ad%co2(15, layer, iprof) * 2 * SQRT(predictors%co2(15, layer, iprof)) * twr(layer)
        twr_ad(layer)                       =      &
          & twr_ad(layer) + predictors_ad%co2(15, layer, iprof) * 2 * SQRT(predictors%co2(15, layer, iprof)) * tr(layer)
        twr_ad(layer)                       = twr_ad(layer) +                             &
          & predictors_ad%co2(14, layer, iprof) * 3 * twr(layer) ** 2 * tr(layer) ** 2 *  &
          & raytracing%pathsat(layer, iprof) ** 0.5
        tr_ad(layer)                        = tr_ad(layer) +      &
          & predictors_ad%co2(14, layer, iprof) * 2 * tr(layer) * twr(layer) ** 3 * SQRT(raytracing%pathsat(layer, iprof))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +         &
          & predictors_ad%co2(14, layer, iprof) * tr(layer) ** 2 * twr(layer) ** 3 * 0.5 /  &
          & SQRT(raytracing%pathsat(layer, iprof))
        tr_ad(layer)                        =      &
          & tr_ad(layer) + predictors_ad%co2(13, layer, iprof) * 3 * tr(layer) ** 2 * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(13, layer, iprof) * tr(layer) ** 3
        tr_ad(layer)                        = tr_ad(layer) + predictors_ad%co2(12, layer, iprof) * 3 * tr(layer) ** 2
        co2r_ad(layer)                      = co2r_ad(layer) +                              &
          & predictors_ad%co2(11, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) /  &
          & SQRT(predictors%co2(1, layer, iprof))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co2(11, layer, iprof) * 0.5 * co2r(layer) / SQRT(predictors%co2(1, layer, iprof))
        twr_ad(layer)                       = twr_ad(layer) +      &
          & predictors_ad%co2(10, layer, iprof) * raytracing%pathsat(layer, iprof) * SQRT(predictors%co2(5, layer, iprof))
        tr_ad(layer)                        = tr_ad(layer) +                               &
          & predictors_ad%co2(10, layer, iprof) * 0.5 * predictors%co2(7, layer, iprof) /  &
          & SQRT(predictors%co2(5, layer, iprof))
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(10, layer, iprof) * twr(layer) * tr(layer) ** 0.5
        twr_ad(layer)                       = twr_ad(layer) + predictors_ad%co2(9, layer, iprof) * 3 * twr(layer) ** 2
        co2w_ad(layer)                      = co2w_ad(layer) +                           &
          & predictors_ad%co2(8, layer, iprof) * 2 * raytracing%pathsat(layer, iprof) *  &
          & SQRT(predictors%co2(8, layer, iprof))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%co2(8, layer, iprof) * 2 * co2w(layer) * SQRT(predictors%co2(8, layer, iprof))
        twr_ad(layer)                       =      &
          & twr_ad(layer) + predictors_ad%co2(7, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(7, layer, iprof) * twr(layer)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(6, layer, iprof)
        tr_ad(layer)                        = tr_ad(layer) + predictors_ad%co2(5, layer, iprof)
        tr_ad(layer)                        =      &
          & tr_ad(layer) + predictors_ad%co2(4, layer, iprof) * 2 * predictors%co2(3, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(4, layer, iprof) * tr(layer) ** 2
        tr_ad(layer)                        =      &
          & tr_ad(layer) + predictors_ad%co2(3, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(3, layer, iprof) * tr(layer)
        tr_ad(layer)                        =      &
          & tr_ad(layer) + predictors_ad%co2(2, layer, iprof) * 2 * predictors%co2(5, layer, iprof)
        co2r_ad(layer)                      =      &
          & co2r_ad(layer) + predictors_ad%co2(1, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%co2(1, layer, iprof) * co2r(layer)
      ENDIF

    ENDDO

!5.5 cloud
!---------

    IF (opts%clw_Data .AND. coef%id_sensor == sensor_id_mw) THEN

      DO layer = 1, prof(1)%nlayers
        deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
      ENDDO


      DO layer = 2, prof(1)%nlayers
        level = layer + 1
        prof_ad(iprof)%clw(level - 1)   =      &
          & prof_ad(iprof)%clw(level - 1) + 0.5_JPRB * predictors_ad%clw(layer, iprof) * deltac(layer) * geom(iprof)%seczen
        predictors_ad%clw(layer, iprof) = 0.5_JPRB * predictors_ad%clw(layer, iprof)
      ENDDO


      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        prof_ad(iprof)%clw(level) =      &
          & prof_ad(iprof)%clw(level) + predictors_ad%clw(layer, iprof) * deltac(layer) * geom(iprof)%seczen
      ENDDO

    ENDIF

!5.4 ozone
!---------

    IF (coef%nozone > 0) THEN
      o_ad(:)      = 0._JPRB
      or_ad(:)     = 0._JPRB
      wr_ad(:)     = 0._JPRB
      ow_ad(:)     = 0._JPRB
      dto_ad(:)    = 0._JPRB
      sec_or_ad(:) = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
! One can pack all ow_ad lines in one longer statement
! same for sec_or_ad and dto_ad
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(15, layer, iprof) * tro(layer) ** 3
        tro_ad(layer)                       =      &
          & tro_ad(layer) + predictors_ad%ozone(15, layer, iprof) * 3 * tro(layer) ** 2 * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +                                        &
          & predictors_ad%ozone(14, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) ** ( - 0.5) * ow(layer) ** 2 *  &
          & dto(layer)
        ow_ad(layer)                        = ow_ad(layer) +      &
          & predictors_ad%ozone(14, layer, iprof) * 2 * ow(layer) * raytracing%pathsat(layer, iprof) ** 0.5 * dto(layer)
        dto_ad(layer)                       =      &
          & dto_ad(layer) + predictors_ad%ozone(14, layer, iprof) * raytracing%pathsat(layer, iprof) ** 0.5 * ow(layer) ** 2
        ow_ad(layer)                        = ow_ad(layer) +                                   &
          & predictors_ad%ozone(13, layer, iprof) * 1.75 * raytracing%pathsat(layer, iprof) *  &
          & (raytracing%pathsat(layer, iprof) * ow(layer)) ** 0.75
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(13, layer, iprof) * 1.75 * ow(layer) *                   &
          & (raytracing%pathsat(layer, iprof) * ow(layer)) ** 0.75
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(12, layer, iprof) * or(layer) / ow(layer)
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(12, layer, iprof) * raytracing%pathsat(layer, iprof) / ow(layer)
        ow_ad(layer)                        = ow_ad(layer) -      &
          & predictors_ad%ozone(12, layer, iprof) * or(layer) * raytracing%pathsat(layer, iprof) / ow(layer) ** 2
        ow_ad(layer)                        = ow_ad(layer) +                                &
          & predictors_ad%ozone(11, layer, iprof) * 2 * raytracing%pathsat(layer, iprof) *  &
          & predictors%ozone(10, layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(11, layer, iprof) * 2 * ow(layer) * predictors%ozone(10, layer, iprof)
        ow_ad(layer)                        =      &
          & ow_ad(layer) + predictors_ad%ozone(10, layer, iprof) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(10, layer, iprof) * ow(layer)
        or_ad(layer)                        = or_ad(layer) +                                             &
          & predictors_ad%ozone(9, layer, iprof) * SQRT(raytracing%pathsat(layer, iprof) * ow(layer)) *  &
          & raytracing%pathsat(layer, iprof)
        ow_ad(layer)                        = ow_ad(layer) +                                  &
          & predictors_ad%ozone(9, layer, iprof) * predictors%ozone(1, layer, iprof) * 0.5 *  &
          & raytracing%pathsat(layer, iprof) / SQRT(raytracing%pathsat(layer, iprof) * ow(layer))
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(9, layer, iprof) * 1.5 * or(layer) * SQRT(predictors%ozone(10, layer, iprof))
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(8, layer, iprof) * predictors%ozone(10, layer, iprof)
        ow_ad(layer)                        =      &
          & ow_ad(layer) + predictors_ad%ozone(8, layer, iprof) * raytracing%pathsat(layer, iprof) * or(layer)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(8, layer, iprof) * or(layer) * ow(layer)
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(7, layer, iprof) * 1.5 * predictors%ozone(2, layer, iprof) / ow(layer)
        ow_ad(layer)                        =      &
          & ow_ad(layer) - predictors_ad%ozone(7, layer, iprof) * predictors%ozone(7, layer, iprof) / ow(layer)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
          & predictors_ad%ozone(7, layer, iprof) * 0.5 * or(layer) ** 1.5 /              &
          & (raytracing%pathsat(layer, iprof) ** 0.5 * ow(layer))
        or_ad(layer)                        =      &
          & or_ad(layer) + predictors_ad%ozone(6, layer, iprof) * 2 * predictors%ozone(1, layer, iprof) * ow(layer)
        ow_ad(layer)                        = ow_ad(layer) +      &
          & predictors_ad%ozone(6, layer, iprof) * predictors%ozone(4, layer, iprof) / raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) =      &
          & raytracing_ad%pathsat(layer, iprof) + predictors_ad%ozone(6, layer, iprof) * or(layer) ** 2 * ow(layer)
        sec_or_ad(layer)                    = sec_or_ad(layer) +                                   &
          & predictors_ad%ozone(5, layer, iprof) * 0.5_JPRB * predictors%ozone(3, layer, iprof) /  &
          & (predictors%ozone(1, layer, iprof) * predictors%ozone(2, layer, iprof))
        dto_ad(layer)                       =      &
          & dto_ad(layer) + predictors_ad%ozone(5, layer, iprof) * predictors%ozone(2, layer, iprof)
        sec_or_ad(layer)                    =      &
          & sec_or_ad(layer) + predictors_ad%ozone(4, layer, iprof) * 2 * predictors%ozone(1, layer, iprof)
        sec_or_ad(layer)                    = sec_or_ad(layer) +      &
          & predictors_ad%ozone(3, layer, iprof) * predictors%ozone(3, layer, iprof) / predictors%ozone(1, layer, iprof)
        dto_ad(layer)                       =      &
          & dto_ad(layer) + predictors_ad%ozone(3, layer, iprof) * predictors%ozone(1, layer, iprof)
        sec_or_ad(layer)                    =      &
          & sec_or_ad(layer) + predictors_ad%ozone(2, layer, iprof) * 0.5_JPRB / predictors%ozone(2, layer, iprof)
        sec_or_ad(layer)                    = sec_or_ad(layer) + predictors_ad%ozone(1, layer, iprof)
        or_ad(layer)                        = or_ad(layer) + sec_or_ad(layer) * raytracing%pathsat(layer, iprof)
        raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + sec_or_ad(layer) * or(layer)
      ENDDO

    ENDIF

!5.3 Water Vapour Continuum based on RTIASI
!------------------------------------------

    IF (coef%nwvcont > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level              = layer + 1
        sec_wr_ad(layer)   = sec_wr_ad(layer) + predictors_ad%wvcont(4, layer, iprof) / tr_sq(layer)
        tr_ad(layer)       = tr_ad(layer) -      &
          & 2 * predictors_ad%wvcont(4, layer, iprof) * predictors%watervapour(7, layer, iprof) / (tr_sq(layer) * tr(layer))
        sec_wr_ad(layer)   = sec_wr_ad(layer) + predictors_ad%wvcont(3, layer, iprof) / tr(layer)
        tr_ad(layer)       =      &
          & tr_ad(layer) - predictors_ad%wvcont(3, layer, iprof) * predictors%watervapour(7, layer, iprof) / tr_sq(layer)
        sec_wrwr_ad(layer) = sec_wrwr_ad(layer) + predictors_ad%wvcont(2, layer, iprof) / tr_4(layer)
        tr_ad(layer)       =      &
          & tr_ad(layer) - 4 * predictors_ad%wvcont(2, layer, iprof) * predictors%wvcont(1, layer, iprof) / tr_4(layer)
        sec_wrwr_ad(layer) = sec_wrwr_ad(layer) + predictors_ad%wvcont(1, layer, iprof) / tr(layer)
        tr_ad(layer)       =      &
          & tr_ad(layer) - predictors_ad%wvcont(1, layer, iprof) * predictors%wvcont(1, layer, iprof) / tr(layer)
      ENDDO

    ENDIF

!
!5.2 water vapour based on RTIASI
!--------------------------------

    DO layer = 1, prof(1)%nlayers
      level = layer + 1
      sec_wr = raytracing%pathsat(layer, iprof) * wr(layer)
      sec_wrwr = sec_wr(layer) * wr(layer)
      wr_ad(layer)                        = wr_ad(layer) +      &
        & predictors_ad%watervapour(19, layer, iprof) * 2 * predictors%watervapour(7, layer, iprof) / ww(layer)
      ww_ad(layer)                        = ww_ad(layer) -                                         &
        & predictors_ad%watervapour(19, layer, iprof) * predictors%watervapour(1, layer, iprof) /  &
        & (raytracing%pathsat(layer, iprof) * ww(layer) ** 2)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%watervapour(19, layer, iprof) * wr(layer) ** 2 / ww(layer)
      dt_ad(layer)                        =      &
        & dt_ad(layer) + predictors_ad%watervapour(18, layer, iprof) * predictors%watervapour(7, layer, iprof) ** 1.5
      wr_ad(layer)                        = wr_ad(layer) +                                               &
        & predictors_ad%watervapour(18, layer, iprof) * 1.5 * predictors%watervapour(5, layer, iprof) *  &
        & raytracing%pathsat(layer, iprof) * dt(layer)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +                                    &
        & predictors_ad%watervapour(18, layer, iprof) * 1.5 * predictors%watervapour(5, layer, iprof) * wr(layer) *  &
        & dt(layer)
      wr_ad(layer)                        =      &
        & wr_ad(layer) + predictors_ad%watervapour(17, layer, iprof) * 1.5 * predictors%watervapour(5, layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(17, layer, iprof) * 0.5 * wr(layer) ** 1.5 / raytracing%pathsat(layer, iprof) ** 0.5
      ww_ad(layer)                        = ww_ad(layer) +                                         &
        & predictors_ad%watervapour(16, layer, iprof) * 1.25 * raytracing%pathsat(layer, iprof) *  &
        & (predictors%watervapour(2, layer, iprof)) ** 0.25
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(16, layer, iprof) * 1.25 * ww(layer) * (predictors%watervapour(2, layer, iprof)) ** 0.25
      wr_ad(layer)                        = wr_ad(layer) +                                        &
        & predictors_ad%watervapour(15, layer, iprof) * 1.5 * raytracing%pathsat(layer, iprof) *  &
        & SQRT(predictors%watervapour(7, layer, iprof))
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(15, layer, iprof) * 1.5 * wr(layer) * SQRT(predictors%watervapour(7, layer, iprof))
      ww_ad(layer)                        = ww_ad(layer) +                                        &
        & predictors_ad%watervapour(14, layer, iprof) * 1.5 * raytracing%pathsat(layer, iprof) *  &
        & SQRT(predictors%watervapour(2, layer, iprof))
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(14, layer, iprof) * 1.5 * ww(layer) * SQRT(predictors%watervapour(2, layer, iprof))
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +                       &
        & predictors_ad%watervapour(13, layer, iprof) * 0.5 / Sqrt(raytracing%pathsat(layer, iprof)) *  &
        & (wr(layer) ** 1.5_JPRB / wwr(layer))
      wr_ad(layer)                        = wr_ad(layer) +                                                                 &
        & predictors_ad%watervapour(13, layer, iprof) * sqrt(raytracing%pathsat(layer, iprof)) * 1.5 * wr(layer) ** 0.5 /  &
        & wwr(layer)
      wwr_ad(layer)                       = wwr_ad(layer) -                                                               &
        & predictors_ad%watervapour(13, layer, iprof) * sqrt(raytracing%pathsat(layer, iprof)) * wr(layer) ** 1.5_JPRB /  &
        & wwr(layer) ** 2
      sec_wrwr_ad(layer)                  =      &
        & sec_wrwr_ad(layer) + predictors_ad%watervapour(12, layer, iprof) / wwr(layer)
      wwr_ad(layer)                       =      &
        & wwr_ad(layer) - predictors_ad%watervapour(12, layer, iprof) * sec_wrwr(layer) / wwr(layer) ** 2
      dt_ad(layer)                        =      &
        & dt_ad(layer) + predictors_ad%watervapour(11, layer, iprof) * predictors%watervapour(5, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) +      &
        & 0.5_JPRB * predictors_ad%watervapour(11, layer, iprof) * dt(layer) / predictors%watervapour(5, layer, iprof)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(10, layer, iprof) * dtabs(layer) * dt(layer)
      dt_ad(layer)                        = dt_ad(layer) +      &
        & 2 * predictors_ad%watervapour(10, layer, iprof) * predictors%watervapour(7, layer, iprof) * dtabs(layer)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(9, layer, iprof) * 4 * predictors%watervapour(8, layer, iprof)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 3 * predictors_ad%watervapour(8, layer, iprof) * predictors%watervapour(1, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + predictors_ad%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) +      &
        & 0.25_JPRB * predictors_ad%watervapour(6, layer, iprof) / predictors%watervapour(6, layer, iprof) ** 3
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 0.5_JPRB * predictors_ad%watervapour(5, layer, iprof) / predictors%watervapour(5, layer, iprof)
      dt_ad(layer)                        =      &
        & dt_ad(layer) + predictors_ad%watervapour(4, layer, iprof) * predictors%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + predictors_ad%watervapour(4, layer, iprof) * dt(layer)
      ww_ad(layer)                        = ww_ad(layer) +                                            &
        & predictors_ad%watervapour(3, layer, iprof) * 2 * predictors%watervapour(2, layer, iprof) *  &
        & raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%watervapour(3, layer, iprof) * 2 * predictors%watervapour(2, layer, iprof) * ww(layer)
      ww_ad(layer)                        =      &
        & ww_ad(layer) + predictors_ad%watervapour(2, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%watervapour(2, layer, iprof) * ww(layer)
      sec_wr_ad(layer)                    =      &
        & sec_wr_ad(layer) + 2 * predictors_ad%watervapour(1, layer, iprof) * predictors%watervapour(7, layer, iprof)
      sec_wr_ad(layer)                    = sec_wr_ad(layer) + sec_wrwr_ad(layer) * wr(layer)
      wr_ad(layer)                        = wr_ad(layer) + sec_wrwr_ad(layer) * predictors%watervapour(7, layer, iprof)
      wr_ad(layer)                        = wr_ad(layer) + sec_wr_ad(layer) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) + sec_wr_ad(layer) * wr(layer)
!5.1 mixed gases
!---------------
! X10
      tr_ad(layer)                        = tr_ad(layer) +      &
        & predictors_ad%mixedgas(10, layer, iprof) * 0.5 * raytracing%pathsat(layer, iprof) ** 1.5 / SQRT(tr(layer))
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(10, layer, iprof) * 1.5 * SQRT(predictors%mixedgas(3, layer, iprof))
! X9
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(9, layer, iprof) * 3 * raytracing%pathsat(layer, iprof) * tr(layer) ** 2
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(9, layer, iprof) * tr(layer) ** 3
! X8
      tuwr_ad(layer)                      =      &
        & tuwr_ad(layer) + predictors_ad%mixedgas(8, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(8, layer, iprof) * tuwr(layer)
! X7
      tuw_ad(layer)                       =      &
        & tuw_ad(layer) + predictors_ad%mixedgas(7, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(7, layer, iprof) * tuw(layer)
! X6
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(6, layer, iprof) * 2 * predictors%mixedgas(5, layer, iprof)
! X5
      tr_ad(layer)                        = tr_ad(layer) + predictors_ad%mixedgas(5, layer, iprof)
! X4
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(4, layer, iprof) * 2 * predictors%mixedgas(3, layer, iprof)
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(4, layer, iprof) * predictors%mixedgas(5, layer, iprof) ** 2
! X3
      tr_ad(layer)                        =      &
        & tr_ad(layer) + predictors_ad%mixedgas(3, layer, iprof) * raytracing%pathsat(layer, iprof)
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(3, layer, iprof) * tr(layer)
! X2
      raytracing_ad%pathsat(layer, iprof) = raytracing_ad%pathsat(layer, iprof) +      &
        & predictors_ad%mixedgas(2, layer, iprof) * 2 * raytracing%pathsat(layer, iprof)
! X1
      raytracing_ad%pathsat(layer, iprof) =      &
        & raytracing_ad%pathsat(layer, iprof) + predictors_ad%mixedgas(1, layer, iprof)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile/reference sums
!-------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + ch4w_ad(layer) / sum2_ch4w(layer)
          ch4_ad(layer) = ch4_ad(layer) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + ch4wr_ad(layer) / sum2_ch4wr(layer)
          ch4_ad(layer) = ch4_ad(layer) + sum1 * coef%dpp(layer - 1) * t(layer)
          t_ad(layer)   = t_ad(layer) + sum1 * coef%dpp(layer - 1) * ch4(layer)
        ENDDO
  
      ELSE
        ch4_ad(:) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + ch4wr_ad(layer) / sum2_ch4wr(layer)
          t_ad(layer)   = t_ad(layer) + sum1 * coef%dpp(layer - 1) * coef%ch4star(layer)
        ENDDO
        
      ENDIF
    ENDIF

    
    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1         = sum1 + cow_ad(layer) / sum2_cow(layer)
          co_ad(layer) = co_ad(layer) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1         = sum1 + cowr_ad(layer) / sum2_cowr(layer)
          co_ad(layer) = co_ad(layer) + sum1 * coef%dpp(layer - 1) * t(layer)
          t_ad(layer)  = t_ad(layer) + sum1 * coef%dpp(layer - 1) * co(layer)
        ENDDO
  
      ELSE
        co_ad(:) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1         = sum1 + cowr_ad(layer) / sum2_cowr(layer)
          t_ad(layer)  = t_ad(layer) + sum1 * coef%dpp(layer - 1) * coef%costar(layer)
        ENDDO

      ENDIF
    ENDIF
    

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + n2ow_ad(layer) / sum2_n2ow(layer)
          n2o_ad(layer) = n2o_ad(layer) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + n2owr_ad(layer) / sum2_n2owr(layer)
          n2o_ad(layer) = n2o_ad(layer) + sum1 * coef%dpp(layer - 1) * t(layer)
          t_ad(layer)   = t_ad(layer) + sum1 * coef%dpp(layer - 1) * n2o(layer)
        ENDDO
  
      ELSE
        n2o_ad(:) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof(1)%nlayers, 1,  - 1
          sum1          = sum1 + n2owr_ad(layer) / sum2_n2owr(layer)
          t_ad(layer)   = t_ad(layer) + sum1 * coef%dpp(layer - 1) * coef%n2ostar(layer)
        ENDDO

      ENDIF
    ENDIF

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof(1)%nlayers, 1,  - 1
        sum1          = sum1 + co2w_ad(layer) / sum2_co2w(layer)
        co2_ad(layer) = co2_ad(layer) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      co2_ad(:) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = prof(1)%nlayers, 1,  - 1
        sum1        = sum1 + ow_ad(layer) / sum2_ow(layer)
        o_ad(layer) = o_ad(layer) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      o_ad(:) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    DO layer = prof(1)%nlayers, 1,  - 1
      sum1        = sum1 + wwr_ad(layer) / sum2_wwr(layer)
      w_ad(layer) = w_ad(layer) + sum1 * coef%dpp(layer - 1) * t(layer)
      t_ad(layer) = t_ad(layer) + sum1 * coef%dpp(layer - 1) * w(layer)
    ENDDO

!
    sum1 = 0._JPRB

    DO layer = prof(1)%nlayers, 1,  - 1
      sum1        = sum1 + ww_ad(layer) / sum2_ww(layer)
      w_ad(layer) = w_ad(layer) + sum1 * coef%dpp(layer - 1)
    ENDDO

    sum1 = 0._JPRB

    IF (coef%nco2 > 0) THEN
      DO layer = prof(1)%nlayers, 1,  - 1
        sum1        = sum1 + twr_ad(layer) / sum2_twr(layer)
        t_ad(layer) = t_ad(layer) + sum1 * coef%dpp(layer - 1)
      ENDDO
    ENDIF

    tuw_ad(2:prof(1)%nlayers) = tuw_ad(2:prof(1)%nlayers) + tuwr_ad(2:prof(1)%nlayers)
    t_ad(1)                   = t_ad(1) + tuwr_ad(1) * coef%dpp(0) / (coef%dpp(0) * coef%tstar(1))
    sum1 = 0._JPRB

    DO layer = prof(1)%nlayers, 1,  - 1
      sum1        = sum1 + tuw_ad(layer) / sum2_tuw(layer)
      t_ad(layer) = t_ad(layer) + sum1
    ENDDO


    DO layer = prof(1)%nlayers, 2,  - 1
      tw_ad(layer - 1) = tw_ad(layer - 1) + tw_ad(layer)
      tr_ad(layer - 1) = tr_ad(layer - 1) + tw_ad(layer) * coef%dpp(layer - 1)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile deviations
!-------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers

    DO layer = 1, prof(1)%nlayers

      IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
        ch4_ad(layer) = ch4_ad(layer) + ch4r_ad(layer) / coef%ch4star(layer)
      ELSE
        ch4_ad(layer) = 0
      ENDIF


      IF (opts%co_Data .AND. coef%nco > 0) THEN
        co_ad(layer) = co_ad(layer) + cor_ad(layer) / coef%costar(layer)
      ELSE
        co_ad(layer) = 0
      ENDIF


      IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
        n2o_ad(layer) = n2o_ad(layer) + n2or_ad(layer) / coef%n2ostar(layer)
      ELSE
        n2o_ad(layer) = 0
      ENDIF


      IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
        co2_ad(layer) = co2_ad(layer) + co2r_ad(layer) / coef%co2star(layer)
      ELSE
        co2_ad(layer) = 0
      ENDIF

      IF (coef%nozone > 0) THEN
        t_ad(layer) = t_ad(layer) + tro_ad(layer) / coef%to3star(layer)
        IF (opts%ozone_Data) THEN
          o_ad(layer) = o_ad(layer) + or_ad(layer) / coef%ostar(layer)
        ELSE
          o_ad(layer) = 0
        ENDIF
      ENDIF

      w_ad(layer) = w_ad(layer) + wr_ad(layer) / coef%wstar(layer)
      t_ad(layer) = t_ad(layer) + tr_ad(layer) / coef%tstar(layer)

      IF (coef%nozone > 0) t_ad(layer) = t_ad(layer) + dto_ad(layer)

      t_ad(layer) = t_ad(layer) + dt_ad(layer)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile layer means
!-------------------------------------------------------------------

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN

      DO level = 2, prof(1)%nlevels

        IF (opts%addpc) THEN
          layer         = level - 1
          ch4_ad(layer) = 0._jprb
        ELSE
          layer = level - 1
          prof_ad(iprof)%ch4(level - 1) = prof_ad(iprof)%ch4(level - 1) + 0.5_JPRB * ch4_ad(layer)
          prof_ad(iprof)%ch4(level)     = prof_ad(iprof)%ch4(level) + 0.5_JPRB * ch4_ad(layer)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN

      DO level = 2, prof(1)%nlevels

        IF (opts%addpc) THEN
          layer         = level - 1
          n2o_ad(layer) = 0._jprb
        ELSE
          layer = level - 1
          prof_ad(iprof)%n2o(level - 1) = prof_ad(iprof)%n2o(level - 1) + 0.5_JPRB * n2o_ad(layer)
          prof_ad(iprof)%n2o(level)     = prof_ad(iprof)%n2o(level) + 0.5_JPRB * n2o_ad(layer)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%co_Data .AND. coef%nco > 0) THEN

      DO level = 2, prof(1)%nlevels

        IF (opts%addpc) THEN
          layer        = level - 1
          co_ad(layer) = 0._jprb
        ELSE
          layer = level - 1
          prof_ad(iprof)%co(level - 1) = prof_ad(iprof)%co(level - 1) + 0.5_JPRB * co_ad(layer)
          prof_ad(iprof)%co(level)     = prof_ad(iprof)%co(level) + 0.5_JPRB * co_ad(layer)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN

      DO level = 2, prof(1)%nlevels

        IF (opts%addpc) THEN
          layer         = level - 1
          co2_ad(layer) = 0._jprb
        ELSE
          layer = level - 1
          prof_ad(iprof)%co2(level - 1) = prof_ad(iprof)%co2(level - 1) + 0.5_JPRB * co2_ad(layer)
          prof_ad(iprof)%co2(level)     = prof_ad(iprof)%co2(level) + 0.5_JPRB * co2_ad(layer)
        ENDIF

      ENDDO

    ENDIF


    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO level = 2, prof(1)%nlevels
        layer = level - 1
        prof_ad(iprof)%o3(level - 1) = prof_ad(iprof)%o3(level - 1) + 0.5_JPRB * o_ad(layer)
        prof_ad(iprof)%o3(level)     = prof_ad(iprof)%o3(level) + 0.5_JPRB * o_ad(layer)
      ENDDO

    ENDIF


    DO level = 2, prof(1)%nlevels
      layer = level - 1
      prof_ad(iprof)%q(level - 1) = prof_ad(iprof)%q(level - 1) + 0.5_JPRB * w_ad(layer)
      prof_ad(iprof)%q(level)     = prof_ad(iprof)%q(level) + 0.5_JPRB * w_ad(layer)
      prof_ad(iprof)%t(level - 1) = prof_ad(iprof)%t(level - 1) + 0.5_JPRB * t_ad(layer)
      prof_ad(iprof)%t(level)     = prof_ad(iprof)%t(level) + 0.5_JPRB * t_ad(layer)
    ENDDO

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9_ad
