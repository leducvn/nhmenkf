!
SUBROUTINE rttov_setpredictors_9_solar( &
            & opts,       &
            & prof,       &
            & coef,       &
            & predictors, &
            & raytracing)
! Description
! RTTOV-8 Model
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
!           --       New routine based on rttov_setpredictors.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   27/02/2009  Profile levels to include ToA. Distinguish arrays
!                    in raytracing (on levels) from all others (on
!                    layers). Predictors prepared to maintain agreement
!                    with the RTTOV-9 prediction scheme (P. Rayer)
!  1.2   02/12/2009  Pathsat, Pathsun and related quantities are now
!                    layer arrays (Marco Matricardi).
!  1.3   05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.4   04/08/2010  Move addsolar check to calling routine (J Hocking)
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
       & rttov_options,   &
       & rttov_coef,      &
       & profile_Type,    &
       & predictors_Type, &
       & raytracing_type
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : max_sol_zen
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(predictors_Type), INTENT(INOUT) :: predictors      ! predictors
  TYPE(raytracing_Type), INTENT(INOUT) :: raytracing      !
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer, iprof
! user profile
  REAL   (KIND=jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2o  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4  (prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tro  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cor  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2or (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuw  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuwr (prof(1)%nlayers)
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cow  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2ow (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4w (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2owr(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cowr (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4wr(prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: sec_or  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sq   (prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles                ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
! 1 profile layer quantities
!------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO iprof = 1, nprofiles
    IF (prof(iprof)%sunzenangle < 0.0 .OR. &
        prof(iprof)%sunzenangle >= max_sol_zen) CYCLE

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
!-Temperature-----------------------------------------------------------------
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) / 2._JPRB
!-H2O--------------------------------------------------------------------------
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) / 2._JPRB
!-O3---------------------------------------------------------------------------
      IF (opts%ozone_Data .AND. coef%nozone > 0)     &
        &  o(layer)   = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) / 2._JPRB
!-CO2--------------------------------------------------------------------------
      IF (opts%co2_Data .AND. coef%nco2 > 0    )     &
        &  co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) / 2._JPRB
!-N2O--------------------------------------------------------------------------
      IF (opts%n2o_Data .AND. coef%nn2o > 0    )     &
        &  n2o(layer) = (prof(iprof)%n2o(level - 1) + prof(iprof)%n2o(level)) / 2._JPRB
!-CO---------------------------------------------------------------------------
      IF (opts%co_Data .AND. coef%nco > 0      )     &
        &  co(layer)  = (prof(iprof)%co(level - 1) + prof(iprof)%co(level)) / 2._JPRB
!-CH4--------------------------------------------------------------------------
      IF (opts%ch4_Data .AND. coef%nch4 > 0    )     &
        &  ch4(layer) = (prof(iprof)%ch4(level - 1) + prof(iprof)%ch4(level)) / 2._JPRB
    ENDDO
! layers
!------------------------------------------------------------------------------
! 2 calculate deviations from reference profile (layers)
! if no input O3 profile we still use the input temperature profile for dto
!-----------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:)    = t(:) - coef%tstar(:)
    dtabs(:) = Abs(dt(:))
    IF (coef%nozone > 0) dto(:) = t(:) - coef%to3star(:)
    
!------------------------------------------------------------------------------
! 3 calculate (profile / reference profile) ratios; tr,wr,or,co2r etc.
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
!-Temperature------------------------------------------------------------------------
    tr(:) = t(:) / coef%tstar(:)
!-H2O-----------------------------------------------------------------------------
    wr(:) = w(:) / coef%wstar(:)
!-Ozone------------------------------------------------------------------------
! if no input O3 profile, set to reference value for or, but still use the input
! temperature profile for tro.

    IF (coef%nozone > 0) THEN
      tro(:) = t(:) / coef%to3star(:)
      IF (opts%ozone_Data) THEN
        or(:)  = o(:) / coef%ostar(:)
      ELSE
        or(:)  = 1._JPRB
      ENDIF
    ENDIF

!-CO2--------------------------------------------------------------------------
! if no input CO2 profile, set to reference value (co2r=1)

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:) = co2(:) / coef%co2star(:)
    ELSE
      co2r(:) = 1._JPRB
    ENDIF

!-N2O-------------------------------------------------------------------------
! if no input N2O profile, set to reference value (n2or=1)

    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:) = n2o(:) / coef%n2ostar(:)
    ELSE
      n2or(:) = 1._JPRB
    ENDIF

!-CO--------------------------------------------------------------------------
! if no input CO profile, set to reference value (cor=1)

    IF (opts%co_Data .AND. coef%nco > 0) THEN
      cor(:) = co(:) / coef%costar(:)
    ELSE
      cor(:) = 1._JPRB
    ENDIF

!-CH4-------------------------------------------------------------------------
! if no input CH4 profile, set to reference value (ch4r=1)

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:) = ch4(:) / coef%ch4star(:)
    ELSE
      ch4r(:) = 1._JPRB
    ENDIF

!------------------------------------------------------------------------------
! 4 calculate profile / reference profile sums: tw,ww,ow,co2w,twr etc.
!------------------------------------------------------------------------------
!-Temperature-----------------------------------------------------------------
    tw(1) = 0._JPRB

    DO layer = 2, prof(1)%nlayers
! cumulate overlying layers (tr relates to same layer as dpp)
! do not need dpp(0) to start
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
! cumulating
      sum1       = sum1 + t(layer)
      sum2       = sum2 + coef%tstar(layer)
      tuw(layer) = sum1 / sum2
    ENDDO

! special case
    tuwr(1)                 = coef%dpp(0) * t(1) / (coef%dpp(0) * coef%tstar(1))
! otherwise
    tuwr(2:prof(1)%nlayers) = tuw(2:prof(1)%nlayers)

    ! twr used only in CO2 predictors
    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (t,tstar relate to layer below dpp)
! need dpp(0) to start
        sum1       = sum1 + coef%dpp(layer - 1) * t(layer)
        sum2       = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        twr(layer) = sum1 / sum2
      ENDDO
    ENDIF

!-H2O--------------------------------------------------------------------------
    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
! cumulate overlying layers: (w,wstar relate to layer below dpp)
! need dpp(0) to start
      sum1      = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2      = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(layer) = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (w,wstar relate to layer below dpp)
! need dpp(0) to start
      sum1       = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2       = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
    ENDDO

!-O3---------------------------------------------------------------------------
! if no input O3 profile, set to reference value (ow =1)

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (o,ostar relate to layer below dpp)
! need dpp(0) to start
        sum1      = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2      = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow(:) = 1._JPRB
    ENDIF

!-CO2-------------------------------------------------------------------------
! if no input co2 profile, set to reference value (co2w=1)

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co2,co2star relate to layer below dpp)
! need dpp(0) to start
        sum1        = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2        = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w(:) = 1._JPRB
    ENDIF

!-N2O-------------------------------------------------------------------------
! if no input n2o profile, set to reference value (n2ow=1)

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (n2o,n2ostar relate to layer below dpp)
! need dpp(0) to start
          sum1        = sum1 + coef%dpp(layer - 1) * n2o(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          n2ow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (n2o,n2ostar relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * n2o(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        n2ow(:)  = 1._JPRB

! no input n2o profile so use the n2o reference profile, 
! but still use the input temperature profile to calculate n2owr
        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (n2o,n2ostar relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF
    
!-CO--------------------------------------------------------------------------
! if no input co profile, set to reference value (cow=1 )

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co,costar relate to layer below dpp)
! need dpp(0) to start
          sum1       = sum1 + coef%dpp(layer - 1) * co(layer)
          sum2       = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          cow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co,costar relate to layer below dpp)
! need dpp(0) to start
          sum1        = sum1 + coef%dpp(layer - 1) * co(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        cow(:)  = 1._JPRB
  
! no input co profile so use the co reference profile, 
! but still use the input temperature profile to calculate cowr
        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co,costar relate to layer below dpp)
! need dpp(0) to start
          sum1        = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF
    
!-CH4-------------------------------------------------------------------------
! if no input ch4 profile, set to reference value (ch4w=1 )

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (ch4,ch4star relate to layer below dpp)
! need dpp(0) to start
          sum1        = sum1 + coef%dpp(layer - 1) * ch4(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          ch4w(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (ch4,ch4star relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * ch4(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        ch4w(:)  = 1._JPRB
  
! no input ch4 profile so use the ch4 reference profile, 
! but still use the input temperature profile to calculate ch4wr
        sum1 = 0._JPRB
        sum2 = 0._JPRB
    
        DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (ch4,ch4star relate to layer below dpp)
! need dpp(0) to start
          sum1         = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr(layer) = sum1 / sum2
        ENDDO
      ENDIF
    ENDIF
  
!5) set predictors for RTTOV-9 options (same for RTTOV-10)
!--

    DO layer = 1, prof(1)%nlayers
      level = layer + 1! raytracing % pathsat(level,i) refers to angle at lower boundary of layer
      raytracing%patheff(layer, iprof)             = raytracing%pathsat(layer, iprof) + raytracing%pathsun(layer, iprof)
!5.1 mixed gases
!---
      tr_sq(layer)                                 = tr(layer) * tr(layer)
      predictors%mixedgas_sun(1, layer, iprof)     = raytracing%patheff(layer, iprof)
      predictors%mixedgas_sun(2, layer, iprof)     = raytracing%patheff(layer, iprof) * raytracing%patheff(layer, iprof)
      predictors%mixedgas_sun(3, layer, iprof)     = raytracing%patheff(layer, iprof) * tr(layer)
      predictors%mixedgas_sun(4, layer, iprof)     = raytracing%patheff(layer, iprof) * tr_sq(layer)
      predictors%mixedgas_sun(5, layer, iprof)     = tr(layer)
      predictors%mixedgas_sun(6, layer, iprof)     = tr_sq(layer)
      predictors%mixedgas_sun(7, layer, iprof)     = raytracing%patheff(layer, iprof) * tuw(layer)
      predictors%mixedgas_sun(8, layer, iprof)     = raytracing%patheff(layer, iprof) * tuwr(layer)
      predictors%mixedgas_sun(9, layer, iprof)     = raytracing%patheff(layer, iprof) * tr(layer) ** 3
      predictors%mixedgas_sun(10, layer, iprof)    =      &
        & raytracing%patheff(layer, iprof) * (raytracing%patheff(layer, iprof) * tr(layer)) ** 0.5
!5.2 water vapour line transmittance based on RTIASI but with pred 9 removed
!----------------
      sec_wr(layer)                                = raytracing%patheff(layer, iprof) * wr(layer)
      sec_wrwr(layer)                              = sec_wr(layer) * wr(layer)
!predictors % watervapour_sun(:,:,iprof) = 0._JPRB
      predictors%watervapour_sun(1, layer, iprof)  = sec_wr(layer) * sec_wr(layer)
      predictors%watervapour_sun(2, layer, iprof)  = raytracing%patheff(layer, iprof) * ww(layer)
      predictors%watervapour_sun(3, layer, iprof)  = (raytracing%patheff(layer, iprof) * ww(layer)) ** 2
      predictors%watervapour_sun(4, layer, iprof)  = sec_wr(layer) * dt(layer)
      predictors%watervapour_sun(5, layer, iprof)  = Sqrt(sec_wr(layer))
      predictors%watervapour_sun(6, layer, iprof)  = sec_wr(layer) ** 0.25_JPRB
      predictors%watervapour_sun(7, layer, iprof)  = sec_wr(layer)
      predictors%watervapour_sun(8, layer, iprof)  = (raytracing%patheff(layer, iprof) * ww(layer)) ** 1.5
      predictors%watervapour_sun(9, layer, iprof)  = sec_wr(layer) ** 1.5
      predictors%watervapour_sun(10, layer, iprof) = sec_wr(layer) * dt(layer) * dtabs(layer)
      predictors%watervapour_sun(11, layer, iprof) = Sqrt(sec_wr(layer)) * dt(layer)
      predictors%watervapour_sun(12, layer, iprof) = (raytracing%patheff(layer, iprof) * ww(layer)) ** 1.25
      predictors%watervapour_sun(13, layer, iprof) = sec_wr(layer) ** 0.5 * wr(layer)
      predictors%watervapour_sun(14, layer, iprof) = sec_wr(layer) ** 1.5 * dt(layer)
      predictors%watervapour_sun(15, layer, iprof) = raytracing%patheff(layer, iprof) * wr(layer) ** 2 / ww(layer)
!5.3 water vapour continuum transmittance based on RTIASI
!----------------
!

      IF (coef%nwvcont > 0) THEN
!predictors % wvcont_sun(:,:,iprof)  = 0._JPRB
        predictors%wvcont_sun(1, layer, iprof) = sec_wrwr(layer) / tr(layer)
        predictors%wvcont_sun(2, layer, iprof) = sec_wrwr(layer) / (tr_sq(layer) * tr_sq(layer))
        predictors%wvcont_sun(3, layer, iprof) = sec_wr(layer) / tr(layer)
        predictors%wvcont_sun(4, layer, iprof) = sec_wr(layer) / tr_sq(layer)
      ENDIF

!5.4 ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

      IF (coef%nozone > 0) THEN
        sec_or(layer)                          = or(layer) * raytracing%patheff(layer, iprof)
        predictors%ozone_sun(1, layer, iprof)  = sec_or(layer)
        predictors%ozone_sun(2, layer, iprof)  = Sqrt(sec_or(layer))
        predictors%ozone_sun(3, layer, iprof)  = sec_or(layer) / ow(layer)
        predictors%ozone_sun(4, layer, iprof)  = sec_or(layer) * sec_or(layer)
        predictors%ozone_sun(5, layer, iprof)  = Sqrt(sec_or(layer)) * dto(layer)
        predictors%ozone_sun(6, layer, iprof)  = sec_or(layer) * or(layer) * ow(layer)
        predictors%ozone_sun(7, layer, iprof)  = Sqrt(sec_or(layer)) * or(layer) / ow(layer)
        predictors%ozone_sun(8, layer, iprof)  = (raytracing%patheff(layer, iprof) * ow(layer)) ** 1.75
        predictors%ozone_sun(9, layer, iprof)  = sec_or(layer) * Sqrt(raytracing%patheff(layer, iprof) * ow(layer))
        predictors%ozone_sun(10, layer, iprof) = raytracing%patheff(layer, iprof) * ow(layer)
        predictors%ozone_sun(11, layer, iprof) =      &
          & raytracing%patheff(layer, iprof) * ow(layer) * raytracing%patheff(layer, iprof) * ow(layer)
        predictors%ozone_sun(12, layer, iprof) = raytracing%patheff(layer, iprof) ** 0.5 * ow(layer) ** 2 * dto(layer)
        predictors%ozone_sun(13, layer, iprof) = raytracing%patheff(layer, iprof) * tro(layer) ** 3
      ENDIF

!5.5 carbon dioxide transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nco2 > 0) THEN
        predictors%co2_sun(1, layer, iprof)  = raytracing%patheff(layer, iprof) * co2r(layer)
        predictors%co2_sun(2, layer, iprof)  = tr_sq(layer)
        predictors%co2_sun(3, layer, iprof)  = raytracing%patheff(layer, iprof) * tr(layer)
        predictors%co2_sun(4, layer, iprof)  = raytracing%patheff(layer, iprof) * tr_sq(layer)
        predictors%co2_sun(5, layer, iprof)  = tr(layer)
        predictors%co2_sun(6, layer, iprof)  = raytracing%patheff(layer, iprof) * twr(layer)
        predictors%co2_sun(7, layer, iprof)  = (raytracing%patheff(layer, iprof) * co2w(layer)) ** 2
        predictors%co2_sun(8, layer, iprof)  = raytracing%patheff(layer, iprof) * cor(layer)
        predictors%co2_sun(9, layer, iprof)  = raytracing%patheff(layer, iprof) * twr(layer) * Sqrt(tr(layer))
        predictors%co2_sun(10, layer, iprof) = (raytracing%patheff(layer, iprof) * co2r(layer)) ** 0.5
        predictors%co2_sun(11, layer, iprof) = tr(layer) ** 3
        predictors%co2_sun(12, layer, iprof) = raytracing%patheff(layer, iprof) * tr(layer) ** 3
        predictors%co2_sun(13, layer, iprof) =      &
          & raytracing%patheff(layer, iprof) ** 0.5 * tr(layer) ** 2 * twr(layer) ** 3
        predictors%co2_sun(14, layer, iprof) = tr(layer) ** 2 * twr(layer) ** 2
      ENDIF

!5.6 n2o transmittance based on RTIASI
!-------------------------------------
!

      IF (coef%nn2o > 0) THEN
        predictors%n2o_sun(1, layer, iprof)  = raytracing%patheff(layer, iprof) * n2or(layer)
        predictors%n2o_sun(2, layer, iprof)  = (raytracing%patheff(layer, iprof) * n2or(layer)) ** 0.5_JPRB
        predictors%n2o_sun(3, layer, iprof)  = raytracing%patheff(layer, iprof) * n2or(layer) * dt(layer)
        predictors%n2o_sun(4, layer, iprof)  = (raytracing%patheff(layer, iprof) * n2or(layer)) ** 2
        predictors%n2o_sun(5, layer, iprof)  = n2or(layer) * dt(layer)
        predictors%n2o_sun(6, layer, iprof)  = (raytracing%patheff(layer, iprof) * n2or(layer)) ** 0.25_JPRB
        predictors%n2o_sun(7, layer, iprof)  = raytracing%patheff(layer, iprof) * n2ow(layer)
        predictors%n2o_sun(8, layer, iprof)  = raytracing%patheff(layer, iprof) * n2owr(layer)
        predictors%n2o_sun(9, layer, iprof)  = (raytracing%patheff(layer, iprof) * n2owr(layer)) ** 2
        predictors%n2o_sun(10, layer, iprof) =      &
          & (raytracing%patheff(layer, iprof) * n2or(layer)) ** 0.5_JPRB * n2or(layer) / n2ow(layer)
        predictors%n2o_sun(11, layer, iprof) = (raytracing%patheff(layer, iprof) * n2owr(layer)) ** 3
        predictors%n2o_sun(12, layer, iprof) = raytracing%patheff(layer, iprof) ** 2 * n2owr(layer) * dt(layer)
      ENDIF

!5.7 co transmittance based on RTIASI
!------------------------------------

      IF (coef%nco > 0) THEN
        predictors%co_sun(1, layer, iprof)  = raytracing%patheff(layer, iprof) * cor(layer)
        predictors%co_sun(2, layer, iprof)  = (raytracing%patheff(layer, iprof) * cor(layer)) ** 0.5
        predictors%co_sun(3, layer, iprof)  = raytracing%patheff(layer, iprof) * cor(layer) * dt(layer)
        predictors%co_sun(4, layer, iprof)  = (raytracing%patheff(layer, iprof) * cor(layer)) ** 2
        predictors%co_sun(5, layer, iprof)  = (raytracing%patheff(layer, iprof) * cor(layer)) ** 0.5 * dt(layer)
        predictors%co_sun(6, layer, iprof)  = (raytracing%patheff(layer, iprof) * cor(layer)) ** 0.25
        predictors%co_sun(7, layer, iprof)  = raytracing%patheff(layer, iprof) * cor(layer) * dt(layer) * ABS(dt(layer))
        predictors%co_sun(8, layer, iprof)  = raytracing%patheff(layer, iprof) * cor(layer) ** 2 / cow(layer)
        predictors%co_sun(9, layer, iprof)  =      &
          & (raytracing%patheff(layer, iprof) * cor(layer)) ** 0.5 * cor(layer) / cow(layer)
        predictors%co_sun(10, layer, iprof) = raytracing%patheff(layer, iprof) * cor(layer) ** 2 / SQRT(cow(layer))
        predictors%co_sun(11, layer, iprof) = (raytracing%patheff(layer, iprof) * cowr(layer)) ** 0.4
        predictors%co_sun(12, layer, iprof) = (raytracing%patheff(layer, iprof) * cowr(layer)) ** 0.25
      ENDIF

!5.8 ch4 transmittance based on RTIASI
!-------------------------------------
!

      IF (coef%nch4 > 0) THEN
        predictors%ch4_sun(1, layer, iprof)  = raytracing%patheff(layer, iprof) * ch4r(layer)
        predictors%ch4_sun(2, layer, iprof)  = (raytracing%patheff(layer, iprof) * ch4r(layer)) ** 0.5
        predictors%ch4_sun(3, layer, iprof)  = raytracing%patheff(layer, iprof) * ch4r(layer) * dt(layer)
        predictors%ch4_sun(4, layer, iprof)  = (raytracing%patheff(layer, iprof) * ch4r(layer)) ** 2
        predictors%ch4_sun(5, layer, iprof)  = ch4r(layer) * dt(layer)
        predictors%ch4_sun(6, layer, iprof)  = (raytracing%patheff(layer, iprof) * ch4r(layer)) ** 0.25
        predictors%ch4_sun(7, layer, iprof)  = raytracing%patheff(layer, iprof) * ch4wr(layer)
        predictors%ch4_sun(8, layer, iprof)  = ch4wr(layer)
        predictors%ch4_sun(9, layer, iprof)  = (raytracing%patheff(layer, iprof) * ch4w(layer)) ** 2
        predictors%ch4_sun(10, layer, iprof) = raytracing%patheff(layer, iprof) * ch4w(layer)
        predictors%ch4_sun(11, layer, iprof) =      &
          & (raytracing%patheff(layer, iprof) * ch4r(layer)) ** 0.5 * ch4r(layer) / ch4w(layer)
      ENDIF

    ENDDO
! layers
  ENDDO
! profiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9_solar
