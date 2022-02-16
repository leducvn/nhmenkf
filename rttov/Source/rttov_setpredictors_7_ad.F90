!
SUBROUTINE rttov_setpredictors_7_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
!
! Description
! RTTOV-7 Model
! AD of rttov_setpredictors
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
! Code based on PRFTAU from previous versions of RTTOV
! Only one profile per call
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       07/12/2005  Add surface humidity (R. Saunders)
!  1.4       03/03/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.6       15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.7       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                        (Marco Matricardi).
!  1.8       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
!                        use with single coefficient file  (P Rayer)                       (Marco Matricardi).
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
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : gravity
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE parkind1, ONLY : JPRB
  USE parkind1, ONLY : JPRB
  USE parkind1, ONLY : JPRB
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_ad(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(profile_aux    ), INTENT(IN)    :: aux                ! auxillary profiles info.
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer
  INTEGER(KIND=jpim) :: iv2lev, iv3lev(size(prof))
  INTEGER(KIND=jpim) :: iv2lay, i
  REAL   (KIND=jprb) :: ztemp
! user profile
  REAL   (KIND=jprb) :: t(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: w(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(size(prof)     , prof(1)%nlayers)
! reference profile
  REAL   (KIND=jprb) :: tr(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: wr(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or(size(prof)     , prof(1)%nlayers)
! user - reference
  REAL   (KIND=jprb) :: dt(size(prof)     , prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=jprb) :: tw(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: sum1     (size(prof)                      ), sum2  (size(prof))
  REAL   (KIND=jprb) :: deltac   (prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_wr   (prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sum2_ww  (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: sum2_ow  (size(prof)     , prof(1)%nlayers)
! TL variables
  REAL   (KIND=jprb) :: t_ad     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: w_ad     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: o_ad     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: tr_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: wr_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: or_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: dt_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: dto_ad   (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: tw_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow_ad    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_or_ad(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_wr_ad(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: zsqrt       , zrecip
  INTEGER(KIND=jpim) :: nprofiles                                                      ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
 
!- End of header --------------------------------------------------------
! profile layer quantities
! Direct variables
!CDIR NOLOOPCHG
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      t(i, layer) = (prof(i)%t(level - 1) + prof(i)%t(level)) * 0.5_JPRB
      w(i, layer) = (prof(i)%q(level - 1) + prof(i)%q(level)) * 0.5_JPRB
    ENDDO

  ENDDO


  IF (opts%use_q2m) THEN

    DO i = 1, nprofiles
! include surface humidity
      iv3lev(i) = aux%s(i)%nearestlev_surf - 1
      iv2lev    = aux%s(i)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay    = iv2lev - 1
        w(i, iv2lay) = (prof(i)%s2m%q + prof(i)%q(iv3lev(i))) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDIF

!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        o(i, layer) = (prof(i)%o3(level - 1) + prof(i)%o3(level)) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDDO


! 2 calculate deviations from reference profile (layers)
! if no input O3 profile, set to reference value (dto =0)
!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers

    DO i = 1, nprofiles
      dt(i, layer) = t(i, layer) - coef%tstar(layer)
! 3 calculate (profile / reference profile) ratios; tr wr or
      tr(i, layer) = t(i, layer) / coef%tstar(layer)
      wr(i, layer) = w(i, layer) / coef%wstar(layer)
! if no input O3 profile, set to reference value (or =1)

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        or(i, layer) = o(i, layer) / coef%ostar(layer)
      ELSE
        or(i, layer) = 1._JPRB
      ENDIF

    ENDDO

  ENDDO

! 4 calculate profile / reference profile sums: tw ww ow

  DO i = 1, nprofiles
    tw(i, 1) = 0._JPRB
  ENDDO


  DO layer = 2, prof(1)%nlayers

    DO i = 1, nprofiles
      tw(i, layer) = tw(i, layer - 1) + coef%dpp(layer - 1) * tr(i, layer - 1)
    ENDDO

  ENDDO

  sum1 = 0._JPRB
  sum2 = 0._JPRB

  DO layer = 1, prof(1)%nlayers

    DO i = 1, nprofiles
      sum1(i)           = sum1(i) + coef%dpp(layer - 1) * w(i, layer)
      sum2(i)           = sum2(i) + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(i, layer) = sum2(i)
      ww(i, layer)      = sum1(i) / sum2(i)
    ENDDO

  ENDDO

! if no input O3 profile, set to reference value (ow =1)
  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
    sum2 = 0._JPRB
  
    DO layer = 1, prof(1)%nlayers
  
      DO i = 1, nprofiles
  
        IF (opts%ozone_Data) THEN
          sum1(i)           = sum1(i) + coef%dpp(layer - 1) * o(i, layer)
          sum2(i)           = sum2(i) + coef%dpp(layer - 1) * coef%ostar(layer)
          sum2_ow(i, layer) = sum2(i)
          ow(i, layer)      = sum1(i) / sum2(i)
        ELSE
          sum2_ow(i, layer) = 0._JPRB
          ow(i, layer)      = 1._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF
  
! Ajoint code
!-------------
  w_ad(:,:)  = 0._JPRB
  wr_ad(:,:) = 0._JPRB
  ww_ad(:,:) = 0._JPRB
  dt_ad(:,:) = 0._JPRB
  t_ad(:,:)  = 0._JPRB
  tr_ad(:,:) = 0._JPRB
  tw_ad(:,:) = 0._JPRB

  IF (coef%nozone > 0) THEN
    o_ad(:,:)   = 0._JPRB
    or_ad(:,:)  = 0._JPRB
    ow_ad(:,:)  = 0._JPRB
    dto_ad(:,:) = 0._JPRB
  ENDIF

!5.4 cloud
!---------
  ztemp = 1.0_JPRB / (4.3429_JPRB * gravity)

  DO layer = 1, prof(1)%nlayers
    deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) * ztemp
  ENDDO

!CDIR NOLOOPCHG

  DO layer = 2, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles

      IF (opts%clw_Data) THEN
        prof_ad(i)%clw(level - 1)   =      &
          & prof_ad(i)%clw(level - 1) + 0.5_JPRB * predictors_ad%clw(layer, i) * deltac(layer) * geom(i)%seczen
        predictors_ad%clw(layer, i) = 0.5_JPRB * predictors_ad%clw(layer, i)
      ENDIF

    ENDDO

  ENDDO

!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles

      IF (opts%clw_Data) THEN
        prof_ad(i)%clw(level) = prof_ad(i)%clw(level) + predictors_ad%clw(layer, i) * deltac(layer) * geom(i)%seczen
      ENDIF

    ENDDO

  ENDDO

!5.3 ozone
!---------

  IF (coef%nozone > 0) THEN

    DO layer = 1, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles
        sec_or_ad(layer)                = 0._JPRB
! One can pack all ow_ad lines in one longer statement
! same for sec_or_ad and dto_ad
        ow_ad(i, layer)                 = ow_ad(i, layer) +      &
          & predictors_ad%ozone(11, layer, i) * 2 * raytracing%pathsat(layer, i) * predictors%ozone(10, layer, i)
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%ozone(11, layer, i) * 2 * ow(i, layer) * predictors%ozone(10, layer, i)
        ow_ad(i, layer)                 =      &
          & ow_ad(i, layer) + predictors_ad%ozone(10, layer, i) * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%ozone(10, layer, i) * ow(i, layer)
        or_ad(i, layer)                 = or_ad(i, layer) +                                         &
          & predictors_ad%ozone(9, layer, i) * SQRT(raytracing%pathsat(layer, i) * ow(i, layer)) *  &
          & raytracing%pathsat(layer, i)
        ow_ad(i, layer)                 = ow_ad(i, layer) +                                                          &
          & predictors_ad%ozone(9, layer, i) * predictors%ozone(1, layer, i) * 0.5 * raytracing%pathsat(layer, i) /  &
          & SQRT(raytracing%pathsat(layer, i) * ow(i, layer))
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%ozone(9, layer, i) * 1.5 * or(i, layer) * SQRT(predictors%ozone(10, layer, i))
        or_ad(i, layer)                 =      &
          & or_ad(i, layer) + predictors_ad%ozone(8, layer, i) * predictors%ozone(10, layer, i)
        ow_ad(i, layer)                 =      &
          & ow_ad(i, layer) + predictors_ad%ozone(8, layer, i) * raytracing%pathsat(layer, i) * or(i, layer)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%ozone(8, layer, i) * or(i, layer) * ow(i, layer)
        or_ad(i, layer)                 =      &
          & or_ad(i, layer) + predictors_ad%ozone(7, layer, i) * 1.5 * predictors%ozone(2, layer, i) / ow(i, layer)
        ow_ad(i, layer)                 =      &
          & ow_ad(i, layer) - predictors_ad%ozone(7, layer, i) * predictors%ozone(7, layer, i) / ow(i, layer)
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%ozone(7, layer, i) * 0.5 * or(i, layer) * sqrt(or(i, layer)) /       &
          & (sqrt(raytracing%pathsat(layer, i)) * ow(i, layer))
        or_ad(i, layer)                 =      &
          & or_ad(i, layer) + predictors_ad%ozone(6, layer, i) * 2 * predictors%ozone(1, layer, i) * ow(i, layer)
        ow_ad(i, layer)                 = ow_ad(i, layer) +      &
          & predictors_ad%ozone(6, layer, i) * predictors%ozone(4, layer, i) / raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%ozone(6, layer, i) * or(i, layer) ** 2 * ow(i, layer)
        sec_or_ad(layer)                = sec_or_ad(layer) +                               &
          & predictors_ad%ozone(5, layer, i) * 0.5_JPRB * predictors%ozone(3, layer, i) /  &
          & (predictors%ozone(1, layer, i) * predictors%ozone(2, layer, i))
        dto_ad(i, layer)                =      &
          & dto_ad(i, layer) + predictors_ad%ozone(5, layer, i) * predictors%ozone(2, layer, i)
        sec_or_ad(layer)                =      &
          & sec_or_ad(layer) + predictors_ad%ozone(4, layer, i) * 2 * predictors%ozone(1, layer, i)
        sec_or_ad(layer)                = sec_or_ad(layer) +      &
          & predictors_ad%ozone(3, layer, i) * predictors%ozone(3, layer, i) / predictors%ozone(1, layer, i)
        dto_ad(i, layer)                =      &
          & dto_ad(i, layer) + predictors_ad%ozone(3, layer, i) * predictors%ozone(1, layer, i)
        sec_or_ad(layer)                =      &
          & sec_or_ad(layer) + predictors_ad%ozone(2, layer, i) * 0.5_JPRB / predictors%ozone(2, layer, i)
        sec_or_ad(layer)                = sec_or_ad(layer) + predictors_ad%ozone(1, layer, i)
        or_ad(i, layer)                 = or_ad(i, layer) + sec_or_ad(layer) * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) + sec_or_ad(layer) * or(i, layer)
      ENDDO

    ENDDO

  ENDIF

!5.2 water vapour  ( numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO layer = 1, prof_ad(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      sec_wr_ad(layer)                       = 0._JPRB
      sec_wr(layer)                          = raytracing%pathsat(layer, i) * wr(i, layer)
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(15, layer, i) * wr(i, layer) / tr(i, layer) ** 4
      wr_ad(i, layer)                        =      &
        & wr_ad(i, layer) + predictors_ad%watervapour(15, layer, i) * sec_wr(layer) / tr(i, layer) ** 4
      tr_ad(i, layer)                        =      &
        & tr_ad(i, layer) - predictors_ad%watervapour(15, layer, i) * 4 * sec_wr(layer) * wr(i, layer) / tr(i, layer) ** 5
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(14, layer, i) * wr(i, layer) / tr(i, layer)
      wr_ad(i, layer)                        =      &
        & wr_ad(i, layer) + predictors_ad%watervapour(14, layer, i) * sec_wr(layer) / tr(i, layer)
      tr_ad(i, layer)                        =      &
        & tr_ad(i, layer) - predictors_ad%watervapour(14, layer, i) * sec_wr(layer) * wr(i, layer) / tr(i, layer) ** 2
      zsqrt = Sqrt(predictors%watervapour(13, layer, i))
      ww_ad(i, layer)                        =      &
        & ww_ad(i, layer) + predictors_ad%watervapour(13, layer, i) * 2 * raytracing%pathsat(layer, i) * zsqrt
      raytracing_ad%pathsat(layer, i)        =      &
        & raytracing_ad%pathsat(layer, i) + predictors_ad%watervapour(13, layer, i) * 2 * ww(i, layer) * zsqrt
      ww_ad(i, layer)                        = ww_ad(i, layer) +                                                            &
        & predictors_ad%watervapour(12, layer, i) * 4 * raytracing%pathsat(layer, i) * predictors%watervapour(12, layer, i) &
        &  / zsqrt
      raytracing_ad%pathsat(layer, i)        = raytracing_ad%pathsat(layer, i) +      &
        & predictors_ad%watervapour(12, layer, i) * 4 * ww(i, layer) * predictors%watervapour(12, layer, i) / zsqrt
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(11, layer, i) * Abs(dt(i, layer)) * dt(i, layer)
      dt_ad(i, layer)                        =      &
        & dt_ad(i, layer) + predictors_ad%watervapour(11, layer, i) * Abs(dt(i, layer)) * 2 * sec_wr(layer)
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(10, layer, i) * 4 * predictors%watervapour(9, layer, i)
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(9, layer, i) * 3 * predictors%watervapour(5, layer, i)
      predictors_ad%watervapour(2, layer, i) =      &
        & predictors_ad%watervapour(2, layer, i) + predictors_ad%watervapour(8, layer, i) * wr(i, layer) / ww(i, layer)
      wr_ad(i, layer)                        =      &
        & wr_ad(i, layer) + predictors_ad%watervapour(8, layer, i) * predictors%watervapour(2, layer, i) / ww(i, layer)
      ww_ad(i, layer)                        = ww_ad(i, layer) -      &
        & predictors_ad%watervapour(8, layer, i) * predictors%watervapour(2, layer, i) * wr(i, layer) / ww(i, layer) ** 2
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(7, layer, i) * 0.25_JPRB / predictors%watervapour(7, layer, i) ** 3
      zrecip = 1.0_JPRB / predictors%watervapour(1, layer, i)
      dt_ad(i, layer)                        =      &
        & dt_ad(i, layer) + predictors_ad%watervapour(6, layer, i) * predictors%watervapour(2, layer, i)
      sec_wr_ad(layer)                       = sec_wr_ad(layer) +      &
        & predictors_ad%watervapour(6, layer, i) * 0.5_JPRB * predictors%watervapour(6, layer, i) * zrecip
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(5, layer, i) * 2 * predictors%watervapour(1, layer, i)
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(4, layer, i) * predictors%watervapour(4, layer, i) * zrecip
      dt_ad(i, layer)                        =      &
        & dt_ad(i, layer) + predictors_ad%watervapour(4, layer, i) * predictors%watervapour(1, layer, i)
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(3, layer, i) * wr(i, layer) / ww(i, layer)
      wr_ad(i, layer)                        =      &
        & wr_ad(i, layer) + predictors_ad%watervapour(3, layer, i) * sec_wr(layer) / ww(i, layer)
      ww_ad(i, layer)                        =      &
        & ww_ad(i, layer) - predictors_ad%watervapour(3, layer, i) * wr(i, layer) * sec_wr(layer) / ww(i, layer) ** 2
      sec_wr_ad(layer)                       =      &
        & sec_wr_ad(layer) + predictors_ad%watervapour(2, layer, i) * 0.5_JPRB / predictors%watervapour(2, layer, i)
      sec_wr_ad(layer)                       = sec_wr_ad(layer) + predictors_ad%watervapour(1, layer, i)
      raytracing_ad%pathsat(layer, i)        = raytracing_ad%pathsat(layer, i) + sec_wr_ad(layer) * wr(i, layer)
      wr_ad(i, layer)                        = wr_ad(i, layer) + sec_wr_ad(layer) * raytracing%pathsat(layer, i)
    ENDDO

  ENDDO

!5.1 mixed gases
!---------------

 IF (coef%id_inst == 10 .and. coef%IncZeeman) THEN
 ! SSMIS with Zeeman coefficient file
 ! geomagnetic field variables (Be, cosbk) are part of the user input

 ! X21 -> X11
    DO layer = 1, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles

 ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        Predictors_ad%mixedgas(20, layer, i)   =      &
          & Predictors_ad%mixedgas(20, layer, i) + Predictors_ad%mixedgas(21, layer, i) * prof(i)%Be
        Predictors_ad%mixedgas(13, layer, i)    =      &
          & Predictors_ad%mixedgas(13, layer, i) + Predictors_ad%mixedgas(20, layer, i) * prof(i)%Be
        raytracing_ad%pathsat(layer, i)        =                                                       &
          & raytracing_ad%pathsat(layer, i) + Predictors_ad%mixedgas(19, layer, i) * prof(i)%Be ** 3 +  &
          & Predictors_ad%mixedgas(18, layer, i) * prof(i)%Be
        Predictors_ad%mixedgas(16, layer, i)    =      &
          & Predictors_ad%mixedgas(16, layer, i) + Predictors_ad%mixedgas(17, layer, i) / prof(i)%Be
        raytracing_ad%pathsat(layer, i)        =      &
          & raytracing_ad%pathsat(layer, i) + Predictors_ad%mixedgas(16, layer, i) / prof(i)%Be
        Predictors_ad%mixedgas(12, layer, i)    =      &
          & Predictors_ad%mixedgas(12, layer, i) + Predictors_ad%mixedgas(15, layer, i) * prof(i)%cosbk ** 2
        Predictors_ad%mixedgas(12, layer, i)    =      &
          & Predictors_ad%mixedgas(12, layer, i) + Predictors_ad%mixedgas(14, layer, i) / prof(i)%Be
        raytracing_ad%pathsat(layer, i)        =                                                          &
          & raytracing_ad%pathsat(layer, i) + Predictors_ad%mixedgas(13, layer, i) * prof(i)%cosbk ** 2 +  &
          & Predictors_ad%mixedgas(12, layer, i) * (300.0_JPRB / t(i, layer))
        t_ad(i, layer)                         =      &
          & t_ad(i, layer) - (Predictors%mixedgas(12, layer, i) / t(i, layer)) * Predictors_ad%mixedgas(12, layer, i)
        raytracing_ad%pathsat(layer, i)        = raytracing_ad%pathsat(layer, i) + Predictors_ad%mixedgas(11, layer, i)

      ENDDO

    ENDDO

 ! X10 -  NB tw(i,1) set to zero in direct

    DO layer = 2, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles

 ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%mixedgas(10, layer, i) * 0.5_JPRB * sqrt(sqrt(tw(i, layer))) / predictors%mixedgas(9, layer, i)
        tw_ad(i, layer)                 = tw_ad(i, layer) +      &
          & predictors_ad%mixedgas(10, layer, i) * 0.25_JPRB * predictors%mixedgas(9, layer, i) * &
          sqrt(sqrt(tw(i, layer)))**(-3_jpim)
      ENDDO
    ENDDO


 ! X9 -> X1
    DO layer = 1, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles
 
 ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%mixedgas(9, layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(8, layer, i) * tw(i, layer) / tr(i, layer)
        tw_ad(i, layer)                 =      &
          & tw_ad(i, layer) + predictors_ad%mixedgas(8, layer, i) * raytracing%pathsat(layer, i) / tr(i, layer)
        tr_ad(i, layer)                 = tr_ad(i, layer) -      &
          & predictors_ad%mixedgas(8, layer, i) * tw(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) ** 2._JPIM
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(7, layer, i) * tw(i, layer)
        tw_ad(i, layer)                 =      &
          & tw_ad(i, layer) + predictors_ad%mixedgas(7, layer, i) * raytracing%pathsat(layer, i)
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(6, layer, i) * 2 * predictors%mixedgas(5, layer, i)
        tr_ad(i, layer)                 = tr_ad(i, layer) + predictors_ad%mixedgas(5, layer, i)
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(4, layer, i) * 2 * predictors%mixedgas(3, layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, i) ** 2
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(3, layer, i) * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(3, layer, i) * tr(i, layer)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(2, layer, i) * 2 * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(1, layer, i)
        Predictors_ad%mixedgas(1:21, layer, i) = 0.0_JPRB
      ENDDO

    ENDDO


  ELSE ! all sensors, except when running SSMIS with Zeeman coefficient file


! X10 -  NB tw(i,1) set to zero in direct

    DO layer = 2, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles

 ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%mixedgas(10, layer, i) * 0.5_JPRB * sqrt(sqrt(tw(i, layer))) / predictors%mixedgas(9, layer, i)
        tw_ad(i, layer)                 = tw_ad(i, layer) +      &
          predictors_ad%mixedgas(10, layer, i) * 0.25_JPRB * predictors%mixedgas(9, layer, i) * &
          sqrt(sqrt(tw(i, layer)))**(-3_jpim)
      ENDDO

    ENDDO


! X14 -> X11
    DO layer = 1, prof_ad(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles

        IF (coef%id_inst == 3 .AND. coef%IncZeeman) THEN
	! AMSU-A with Zeeman coefficient file
 ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
          raytracing_ad%pathsat(layer, i)         = raytracing_ad%pathsat(layer, i) +                          &
            & 2.0_JPRB * (prof(i)%cosbk * prof(i)%Be) ** 2 * raytracing%pathsat(layer, i) *                    &
            & predictors_ad%mixedgas(14, layer, i) + prof(i)%Be ** 3 * predictors_ad%mixedgas(13, layer, i) +  &
            & 2.0_JPRB * prof(i)%Be * raytracing%pathsat(layer, i) * predictors_ad%mixedgas(12, layer, i) +    &
            & prof(i)%cosbk ** 2 * predictors_ad%mixedgas(11, layer, i)
          predictors_ad%mixedgas(11:14, layer, i) = 0.0_JPRB
 ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
        ENDIF

 ! X9 -> X1...
 ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) +      &
          & predictors_ad%mixedgas(9, layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(8, layer, i) * tw(i, layer) / tr(i, layer)
        tw_ad(i, layer)                 =      &
          & tw_ad(i, layer) + predictors_ad%mixedgas(8, layer, i) * raytracing%pathsat(layer, i) / tr(i, layer)
        tr_ad(i, layer)                 = tr_ad(i, layer) -      &
          & predictors_ad%mixedgas(8, layer, i) * tw(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) ** 2._JPIM
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(7, layer, i) * tw(i, layer)
        tw_ad(i, layer)                 =      &
          & tw_ad(i, layer) + predictors_ad%mixedgas(7, layer, i) * raytracing%pathsat(layer, i)
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(6, layer, i) * 2 * predictors%mixedgas(5, layer, i)
        tr_ad(i, layer)                 = tr_ad(i, layer) + predictors_ad%mixedgas(5, layer, i)
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(4, layer, i) * 2 * predictors%mixedgas(3, layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, i) ** 2
        tr_ad(i, layer)                 =      &
          & tr_ad(i, layer) + predictors_ad%mixedgas(3, layer, i) * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(3, layer, i) * tr(i, layer)
        raytracing_ad%pathsat(layer, i) =      &
          & raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(2, layer, i) * 2 * raytracing%pathsat(layer, i)
        raytracing_ad%pathsat(layer, i) = raytracing_ad%pathsat(layer, i) + predictors_ad%mixedgas(1, layer, i)
      ENDDO

    ENDDO


  ENDIF

  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
  
    DO layer = prof_ad(1)%nlayers, 1,  - 1
  
      DO i = 1, nprofiles
  
        IF (opts%ozone_Data) THEN
          sum1(i)        = sum1(i) + ow_ad(i, layer) / sum2_ow(i, layer)
          o_ad(i, layer) = o_ad(i, layer) + sum1(i) * coef%dpp(layer - 1)
        ELSE
          o_ad(i, layer) = 0._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF
  
  sum1 = 0._JPRB

  DO layer = prof_ad(1)%nlayers, 1,  - 1

    DO i = 1, nprofiles
      sum1(i)        = sum1(i) + ww_ad(i, layer) / sum2_ww(i, layer)
      w_ad(i, layer) = w_ad(i, layer) + sum1(i) * coef%dpp(layer - 1)
    ENDDO

  ENDDO


  DO layer = prof_ad(1)%nlayers, 2,  - 1

    DO i = 1, nprofiles
      tw_ad(i, layer - 1) = tw_ad(i, layer - 1) + tw_ad(i, layer)
      tr_ad(i, layer - 1) = tr_ad(i, layer - 1) + tw_ad(i, layer) * coef%dpp(layer - 1)
    ENDDO

  ENDDO

!CDIR NOLOOPCHG

  IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
    DO layer = 1, prof(1)%nlayers
  
      DO i = 1, nprofiles

        o_ad(i, layer) = o_ad(i, layer) + or_ad(i, layer) / coef%ostar(layer)

      ENDDO
  
    ENDDO
  ENDIF

!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers

    DO i = 1, nprofiles
      w_ad(i, layer) = w_ad(i, layer) + wr_ad(i, layer) / coef%wstar(layer)
      t_ad(i, layer) = t_ad(i, layer) + tr_ad(i, layer) / coef%tstar(layer)
    ENDDO

  ENDDO


  DO layer = 1, prof(1)%nlayers

    DO i = 1, nprofiles

      IF (coef%nozone > 0) t_ad(i, layer) = t_ad(i, layer) + dto_ad(i, layer)

      t_ad(i, layer) = t_ad(i, layer) + dt_ad(i, layer)
    ENDDO

  ENDDO


  IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
    DO level = 2, prof(1)%nlevels
      layer = level - 1
  
      DO i = 1, nprofiles
  
        
        prof_ad(i)%o3(level - 1) = prof_ad(i)%o3(level - 1) + 0.5_JPRB * o_ad(i, layer)
        prof_ad(i)%o3(level)     = prof_ad(i)%o3(level) + 0.5_JPRB * o_ad(i, layer)
  
      ENDDO
  
    ENDDO
  ENDIF

!
! include (or not) adjoint surface humidity
!

  IF (opts%use_q2m) THEN

    DO i = 1, nprofiles
! include surface humidity
      iv2lev                  = aux%s(i)%nearestlev_surf
      iv2lay                  = iv2lev - 1
      prof_ad(i)%s2m%q        = prof_ad(i)%s2m%q + 0.5_JPRB * w_ad(i, iv2lay)
      prof_ad(i)%q(iv3lev(i)) = prof_ad(i)%q(iv3lev(i)) + 0.5_JPRB * w_ad(i, iv2lay)
    ENDDO


    DO i = 1, nprofiles

      DO level = 2, iv3lev(i)
        layer = level - 1
        prof_ad(i)%q(level - 1) = prof_ad(i)%q(level - 1) + 0.5_JPRB * w_ad(i, layer)
        prof_ad(i)%q(level)     = prof_ad(i)%q(level) + 0.5_JPRB * w_ad(i, layer)
      ENDDO

    ENDDO

  ELSE

    DO level = 2, prof_ad(1)%nlevels
      layer = level - 1

      DO i = 1, nprofiles
        prof_ad(i)%q(level - 1) = prof_ad(i)%q(level - 1) + 0.5_JPRB * w_ad(i, layer)
        prof_ad(i)%q(level)     = prof_ad(i)%q(level) + 0.5_JPRB * w_ad(i, layer)
      ENDDO

    ENDDO

  ENDIF


  DO level = 2, prof_ad(1)%nlevels
    layer = level - 1

    DO i = 1, nprofiles
      prof_ad(i)%t(level - 1) = prof_ad(i)%t(level - 1) + 0.5_JPRB * t_ad(i, layer)
      prof_ad(i)%t(level)     = prof_ad(i)%t(level) + 0.5_JPRB * t_ad(i, layer)
    ENDDO

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_ad
