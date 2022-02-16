!
SUBROUTINE rttov_setpredictors_7_k( &
            & opts,         &
            & nlayers,      &
            & angles,       &
            & chanprof,     &
            & profiles,     &
            & profiles_k,   &
            & coef,         &
            & aux_prof,     &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
! Description
! RTTOV-9 Model
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
!    Copyright 2005, EUMETSAT, All Rights Reserved.
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
!  1.0       22/06/2005  initial (P Brunel)
!                        based on version 1.2 (29/03/05) of AD code
!  1.1       07/12/2005  Add surface humidity (R. Saunders)
!  1.2       06/03/2006  Add RTIASI features (RTTOV-8M) (R. Saunders)
!  1.3       12/02/2007  Removed polarisation indexing ( R Saunders)
!  1.4       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.5       15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.6       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                       (Marco Matricardi).
!  1.7       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
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
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options)  , INTENT(IN)    :: opts
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                   ! Number of layers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type   ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(geometry_Type  ), INTENT(IN)    :: angles    (size(profiles))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(profile_aux    ), INTENT(IN)    :: aux_prof                  ! auxillary profiles info.
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer
  INTEGER(KIND=jpim) :: i                                                        ! channel indice
  INTEGER(KIND=jpim) :: j                                                        ! profile indice
  INTEGER(KIND=jpim) :: iv2lev , iv3lev(size(profiles))
  INTEGER(KIND=jpim) :: iv2lay
! user profile
  REAL   (KIND=jprb) :: t(size(profiles), nlayers)
  REAL   (KIND=jprb) :: w(size(profiles), nlayers)
  REAL   (KIND=jprb) :: o(size(profiles), nlayers)
! reference profile
  REAL   (KIND=jprb) :: tr      (size(profiles), nlayers)
  REAL   (KIND=jprb) :: wr      (size(profiles), nlayers)
  REAL   (KIND=jprb) :: or      (size(profiles), nlayers)
! user - reference
  REAL   (KIND=jprb) :: dt      (size(profiles), nlayers)
! pressure weighted
  REAL   (KIND=jprb) :: tw      (size(profiles), nlayers)
  REAL   (KIND=jprb) :: ww      (size(profiles), nlayers)
  REAL   (KIND=jprb) :: ow      (size(profiles), nlayers)
  REAL   (KIND=jprb) :: sum1    (size(chanprof)         ), sum2  (size(chanprof))
  REAL   (KIND=jprb) :: deltac  (nlayers                )
  REAL   (KIND=jprb) :: sec_wr  (nlayers                )
  REAL   (KIND=jprb) :: sum2_ww (size(profiles), nlayers)
  REAL   (KIND=jprb) :: sum2_ow (size(profiles), nlayers)
! K variables
  REAL   (KIND=jprb) :: t_k     (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: w_k     (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: o_k     (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: tr_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: wr_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: or_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: dt_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: dto_k   (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: tw_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: ww_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: ow_k    (size(chanprof), nlayers)
  REAL   (KIND=jprb) :: sec_or_k(nlayers                )
  REAL   (KIND=jprb) :: sec_wr_k(nlayers                )
  REAL   (KIND=jprb) :: zsqrt       , zrecip
  INTEGER(KIND=jpim) :: nprofiles                                                ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!
! Keep use of prof%nlevels in the code instead of input argument nlevels
! This is to allow profiles on variable levels, for future version.
!!!profiles_k(:) % nlevels =nlevels
!
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_K', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(profiles)
  nchannels = size(chanprof)
! profile layer quantities
! Direct variables
!CDIR NOLOOPCHG

  DO layer = 1, profiles(1)%nlayers
    level = layer + 1

    DO j = 1, nprofiles
      t(j, layer) = (profiles(j)%t(level - 1) + profiles(j)%t(level)) * 0.5_JPRB
      w(j, layer) = (profiles(j)%q(level - 1) + profiles(j)%q(level)) * 0.5_JPRB
    ENDDO

  ENDDO

!

  IF (opts%use_q2m) THEN

    DO j = 1, nprofiles
! include surface humidity
      iv3lev(j) = aux_prof%s(j)%nearestlev_surf - 1
      iv2lev    = aux_prof%s(j)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay    = iv2lev - 1
        w(j, iv2lay) = (profiles(j)%s2m%q + profiles(j)%q(iv3lev(j))) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDIF

!
!CDIR NOLOOPCHG

  DO layer = 1, profiles(1)%nlayers
    level = layer + 1

    DO j = 1, nprofiles

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        o(j, layer) = (profiles(j)%o3(level - 1) + profiles(j)%o3(level)) * 0.5_JPRB
      ELSE
        o(j, layer) = 0.0_JPRB
      ENDIF

    ENDDO

  ENDDO


!3) calculate deviations from reference profile (layers)
! if no input O3 profile, set to reference value (dto =0)
! Direct variables
!CDIR NOLOOPCHG

  DO layer = 1, profiles(1)%nlayers

    DO j = 1, nprofiles
      dt(j, layer) = t(j, layer) - coef%tstar(layer)
!2) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
! Direct variables
      tr(j, layer) = t(j, layer) / coef%tstar(layer)
      wr(j, layer) = w(j, layer) / coef%wstar(layer)
! if no input O3 profile, set to reference value (or =1)

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        or(j, layer) = o(j, layer) / coef%ostar(layer)
      ELSE
        or(j, layer) = 1._JPRB
      ENDIF

    ENDDO

  ENDDO

! calculate profile / reference profile sums: tw ww ow
! if no input O3 profile, set to reference value (ow =1)
! Direct variables

  DO j = 1, nprofiles
    tw(j, 1) = 0._JPRB
  ENDDO


  DO layer = 2, profiles(1)%nlayers

    DO j = 1, nprofiles
      tw(j, layer) = tw(j, layer - 1) + coef%dpp(layer - 1) * tr(j, layer - 1)
    ENDDO

  ENDDO

  sum1 = 0._JPRB
  sum2 = 0._JPRB

  DO layer = 1, profiles(1)%nlayers

    DO j = 1, nprofiles
      sum1(j)           = sum1(j) + coef%dpp(layer - 1) * w(j, layer)
      sum2(j)           = sum2(j) + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(j, layer) = sum2(j)
      ww(j, layer)      = sum1(j) / sum2(j)
    ENDDO

  ENDDO

! if no input O3 profile, set to reference value (ow =1)
  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
    sum2 = 0._JPRB
  
    DO layer = 1, profiles(1)%nlayers
  
      DO j = 1, nprofiles
  
        IF (opts%ozone_Data) THEN
          sum1(j)           = sum1(j) + coef%dpp(layer - 1) * o(j, layer)
          sum2(j)           = sum2(j) + coef%dpp(layer - 1) * coef%ostar(layer)
          sum2_ow(j, layer) = sum2(j)
          ow(j, layer)      = sum1(j) / sum2(j)
        ELSE
          sum2_ow(j, layer) = 0.0_JPRB
          ow(j, layer)      = 1.0_JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF
  
! Ajoint code
!-------------
  w_k(:,:)  = 0._JPRB
  wr_k(:,:) = 0._JPRB
  ww_k(:,:) = 0._JPRB
  dt_k(:,:) = 0._JPRB
  t_k(:,:)  = 0._JPRB
  tr_k(:,:) = 0._JPRB
  tw_k(:,:) = 0._JPRB

  IF (coef%nozone > 0) THEN
    o_k(:,:)   = 0._JPRB
    or_k(:,:)  = 0._JPRB
    ow_k(:,:)  = 0._JPRB
    dto_k(:,:) = 0._JPRB
  ENDIF

!5.4 cloud
!---------

  IF (coef%id_sensor == sensor_id_mw) THEN

    DO layer = 1, profiles(1)%nlayers
      deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
    ENDDO

!CDIR NOLOOPCHG

    DO layer = 2, profiles(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

        IF (opts%clw_Data) THEN
          profiles_k(i)%clw(level - 1) =      &
            & profiles_k(i)%clw(level - 1) + 0.5_JPRB * predictors_k%clw(layer, i) * deltac(layer) * angles(j)%seczen
          predictors_k%clw(layer, i)   = 0.5_JPRB * predictors_k%clw(layer, i)
        ENDIF

      ENDDO

    ENDDO

!CDIR NOLOOPCHG

    DO layer = 1, profiles(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

        IF (opts%clw_Data) THEN
          profiles_k(i)%clw(level) =      &
            & profiles_k(i)%clw(level) + predictors_k%clw(layer, i) * deltac(layer) * angles(j)%seczen
        ENDIF

      ENDDO

    ENDDO

  ENDIF

!5.3 ozone
!---------

  IF (coef%nozone > 0) THEN

    DO layer = 1, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof
        sec_or_k(layer)                = 0._JPRB
! One can pack all ow_k lines in one longer statement
! same for sec_or_k and dto_k
        ow_k(i, layer)                 = ow_k(i, layer) +      &
          & predictors_k%ozone(11, layer, i) * 2 * raytracing%pathsat(layer, j) * predictors%ozone(10, layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(11, layer, i) * 2 * ow(j, layer) * predictors%ozone(10, layer, j)
        ow_k(i, layer)                 =      &
          & ow_k(i, layer) + predictors_k%ozone(10, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(10, layer, i) * ow(j, layer)
        or_k(i, layer)                 = or_k(i, layer) +                                          &
          & predictors_k%ozone(9, layer, i) * SQRT(raytracing%pathsat(layer, j) * ow(j, layer)) *  &
          & raytracing%pathsat(layer, j)
        ow_k(i, layer)                 = ow_k(i, layer) +                                                           &
          & predictors_k%ozone(9, layer, i) * predictors%ozone(1, layer, j) * 0.5 * raytracing%pathsat(layer, j) /  &
          & SQRT(raytracing%pathsat(layer, j) * ow(j, layer))
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(9, layer, i) * 1.5 * or(j, layer) * SQRT(predictors%ozone(10, layer, j))
        or_k(i, layer)                 =      &
          & or_k(i, layer) + predictors_k%ozone(8, layer, i) * predictors%ozone(10, layer, j)
        ow_k(i, layer)                 =      &
          & ow_k(i, layer) + predictors_k%ozone(8, layer, i) * raytracing%pathsat(layer, j) * or(j, layer)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(8, layer, i) * or(j, layer) * ow(j, layer)
        or_k(i, layer)                 =      &
          & or_k(i, layer) + predictors_k%ozone(7, layer, i) * 1.5 * predictors%ozone(2, layer, j) / ow(j, layer)
        ow_k(i, layer)                 =      &
          & ow_k(i, layer) - predictors_k%ozone(7, layer, i) * predictors%ozone(7, layer, j) / ow(j, layer)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%ozone(7, layer, i) * 0.5 * or(j, layer) ** 1.5 /      &
          & (raytracing%pathsat(layer, j) ** 0.5 * ow(j, layer))
        or_k(i, layer)                 =      &
          & or_k(i, layer) + predictors_k%ozone(6, layer, i) * 2 * predictors%ozone(1, layer, j) * ow(j, layer)
        ow_k(i, layer)                 =      &
          & ow_k(i, layer) + predictors_k%ozone(6, layer, i) * predictors%ozone(4, layer, j) / raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%ozone(6, layer, i) * or(j, layer) ** 2 * ow(j, layer)
        sec_or_k(layer)                = sec_or_k(layer) +                                &
          & predictors_k%ozone(5, layer, i) * 0.5_JPRB * predictors%ozone(3, layer, j) /  &
          & (predictors%ozone(1, layer, j) * predictors%ozone(2, layer, j))
        dto_k(i, layer)                =      &
          & dto_k(i, layer) + predictors_k%ozone(5, layer, i) * predictors%ozone(2, layer, j)
        sec_or_k(layer)                =      &
          & sec_or_k(layer) + predictors_k%ozone(4, layer, i) * 2 * predictors%ozone(1, layer, j)
        sec_or_k(layer)                = sec_or_k(layer) +      &
          & predictors_k%ozone(3, layer, i) * predictors%ozone(3, layer, j) / predictors%ozone(1, layer, j)
        dto_k(i, layer)                =      &
          & dto_k(i, layer) + predictors_k%ozone(3, layer, i) * predictors%ozone(1, layer, j)
        sec_or_k(layer)                =      &
          & sec_or_k(layer) + predictors_k%ozone(2, layer, i) * 0.5_JPRB / predictors%ozone(2, layer, j)
        sec_or_k(layer)                = sec_or_k(layer) + predictors_k%ozone(1, layer, i)
        or_k(i, layer)                 = or_k(i, layer) + sec_or_k(layer) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + sec_or_k(layer) * or(j, layer)
      ENDDO

    ENDDO

 ENDIF

!5.2 water vapour  ( numbers in right hand are predictorsictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO layer = 1, profiles_k(1)%nlayers
    level = layer + 1

    DO i = 1, nchannels
      j = chanprof(i)%prof
      sec_wr_k(layer)                       = 0._JPRB
      sec_wr(layer)                         = raytracing%pathsat(layer, j) * wr(j, layer)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(15, layer, i) * wr(j, layer) / tr(j, layer) ** 4
      wr_k(i, layer)                        =      &
        & wr_k(i, layer) + predictors_k%watervapour(15, layer, i) * sec_wr(layer) / tr(j, layer) ** 4
      tr_k(i, layer)                        =      &
        & tr_k(i, layer) - predictors_k%watervapour(15, layer, i) * 4 * sec_wr(layer) * wr(j, layer) / tr(j, layer) ** 5
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(14, layer, i) * wr(j, layer) / tr(j, layer)
      wr_k(i, layer)                        =      &
        & wr_k(i, layer) + predictors_k%watervapour(14, layer, i) * sec_wr(layer) / tr(j, layer)
      tr_k(i, layer)                        =      &
        & tr_k(i, layer) - predictors_k%watervapour(14, layer, i) * sec_wr(layer) * wr(j, layer) / tr(j, layer) ** 2
      zsqrt = Sqrt(predictors%watervapour(13, layer, j))
      ww_k(i, layer)                        =      &
        & ww_k(i, layer) + predictors_k%watervapour(13, layer, i) * 2 * raytracing%pathsat(layer, j) * zsqrt
      raytracing_k%pathsat(layer, i)        =      &
        & raytracing_k%pathsat(layer, i) + predictors_k%watervapour(13, layer, i) * 2 * ww(j, layer) * zsqrt
      ww_k(i, layer)                        = ww_k(i, layer) +                                                             &
        & predictors_k%watervapour(12, layer, i) * 4 * raytracing%pathsat(layer, j) * predictors%watervapour(12, layer, j) &
        &  / zsqrt
      raytracing_k%pathsat(layer, i)        = raytracing_k%pathsat(layer, i) +      &
        & predictors_k%watervapour(12, layer, i) * 4 * ww(j, layer) * predictors%watervapour(12, layer, j) / zsqrt
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(11, layer, i) * Abs(dt(j, layer)) * dt(j, layer)
      dt_k(i, layer)                        =      &
        & dt_k(i, layer) + predictors_k%watervapour(11, layer, i) * Abs(dt(j, layer)) * 2 * sec_wr(layer)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(10, layer, i) * 4 * predictors%watervapour(9, layer, j)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(9, layer, i) * 3 * predictors%watervapour(5, layer, j)
      predictors_k%watervapour(2, layer, i) =      &
        & predictors_k%watervapour(2, layer, i) + predictors_k%watervapour(8, layer, i) * wr(j, layer) / ww(j, layer)
      wr_k(i, layer)                        =      &
        & wr_k(i, layer) + predictors_k%watervapour(8, layer, i) * predictors%watervapour(2, layer, j) / ww(j, layer)
      ww_k(i, layer)                        = ww_k(i, layer) -      &
        & predictors_k%watervapour(8, layer, i) * predictors%watervapour(2, layer, j) * wr(j, layer) / ww(j, layer) ** 2
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(7, layer, i) * 0.25_JPRB / predictors%watervapour(7, layer, j) ** 3
      zrecip = 1.0_JPRB / predictors%watervapour(1, layer, j)
      dt_k(i, layer)                        =      &
        & dt_k(i, layer) + predictors_k%watervapour(6, layer, i) * predictors%watervapour(2, layer, j)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(6, layer, i) * 0.5_JPRB * predictors%watervapour(6, layer, j) * zrecip
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(5, layer, i) * 2 * predictors%watervapour(1, layer, j)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(4, layer, i) * predictors%watervapour(4, layer, j) * zrecip
      dt_k(i, layer)                        =      &
        & dt_k(i, layer) + predictors_k%watervapour(4, layer, i) * predictors%watervapour(1, layer, j)
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(3, layer, i) * wr(j, layer) / ww(j, layer)
      wr_k(i, layer)                        =      &
        & wr_k(i, layer) + predictors_k%watervapour(3, layer, i) * sec_wr(layer) / ww(j, layer)
      ww_k(i, layer)                        =      &
        & ww_k(i, layer) - predictors_k%watervapour(3, layer, i) * wr(j, layer) * sec_wr(layer) / ww(j, layer) ** 2
      sec_wr_k(layer)                       =      &
        & sec_wr_k(layer) + predictors_k%watervapour(2, layer, i) * 0.5_JPRB / predictors%watervapour(2, layer, j)
      sec_wr_k(layer)                       = sec_wr_k(layer) + predictors_k%watervapour(1, layer, i)
      raytracing_k%pathsat(layer, i)        = raytracing_k%pathsat(layer, i) + sec_wr_k(layer) * wr(j, layer)
      wr_k(i, layer)                        = wr_k(i, layer) + sec_wr_k(layer) * raytracing%pathsat(layer, j)
    ENDDO

  ENDDO


!5.1 mixed gases
!---------------

  IF (coef%id_inst == 10 .and. coef%IncZeeman) THEN
  ! SSMIS with Zeeman coefficient file
  ! geomagnetic field variables (Be, cosbk) are part of the user input

  ! X21 -> X11
    DO layer = 1, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

  ! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        Predictors_k%mixedgas(20, layer, i)   =      &
          & Predictors_k%mixedgas(20, layer, i) + Predictors_k%mixedgas(21, layer, i) * profiles(j)%Be
        Predictors_k%mixedgas(13, layer, i)    =      &
          & Predictors_k%mixedgas(13, layer, i) + Predictors_k%mixedgas(20, layer, i) * profiles(j)%Be
        raytracing_k%pathsat(layer, i)        =                                                          &
          & raytracing_k%pathsat(layer, i) + Predictors_k%mixedgas(19, layer, i) * profiles(j)%Be ** 3 +  &
          & Predictors_k%mixedgas(18, layer, i) * profiles(j)%Be
        Predictors_k%mixedgas(16, layer, i)    =      &
          & Predictors_k%mixedgas(16, layer, i) + Predictors_k%mixedgas(17, layer, i) / profiles(j)%Be
        raytracing_k%pathsat(layer, i)        =      &
          & raytracing_k%pathsat(layer, i) + Predictors_k%mixedgas(16, layer, i) / profiles(j)%Be
        Predictors_k%mixedgas(12, layer, i)    =      &
          & Predictors_k%mixedgas(12, layer, i) + Predictors_k%mixedgas(15, layer, i) * profiles(j)%cosbk ** 2
        Predictors_k%mixedgas(12, layer, i)    =      &
          & Predictors_k%mixedgas(12, layer, i) + Predictors_k%mixedgas(14, layer, i) / profiles(j)%Be
        raytracing_k%pathsat(layer, i)        =                                                             &
          & raytracing_k%pathsat(layer, i) + Predictors_k%mixedgas(13, layer, i) * profiles(j)%cosbk ** 2 +  &
          & Predictors_k%mixedgas(12, layer, i) * (300.0_JPRB / t(j, layer))
        t_k(i, layer)                         =      &
          & t_k(i, layer) - (Predictors%mixedgas(12, layer, j) / t(j, layer)) * Predictors_k%mixedgas(12, layer, i)
        raytracing_k%pathsat(layer, i)        = raytracing_k%pathsat(layer, i) + Predictors_k%mixedgas(11, layer, i)

      ENDDO

    ENDDO

  ! X10 -  NB tw(i,1) set to zero in direct

    DO layer = 2, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof
 
  ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%mixedgas(10, layer, i) * 0.5_JPRB * tw(j, layer) ** 0.25_JPRB / predictors%mixedgas(9, layer, j)
        tw_k(i, layer)                 = tw_k(i, layer) +      &
          & predictors_k%mixedgas(10, layer, i) * 0.25_JPRB * predictors%mixedgas(9, layer, j) / tw(j, layer) ** 0.75
      ENDDO

    ENDDO


  ! X9 -> X1
    DO layer = 1, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

  ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%mixedgas(9, layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(8, layer, i) * tw(j, layer) / tr(j, layer)
        tw_k(i, layer)                 =      &
          & tw_k(i, layer) + predictors_k%mixedgas(8, layer, i) * raytracing%pathsat(layer, j) / tr(j, layer)
        tr_k(i, layer)                 = tr_k(i, layer) -      &
          & predictors_k%mixedgas(8, layer, i) * tw(j, layer) * raytracing%pathsat(layer, j) / tr(j, layer) ** 2._JPRB
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(7, layer, i) * tw(j, layer)
        tw_k(i, layer)                 =      &
          & tw_k(i, layer) + predictors_k%mixedgas(7, layer, i) * raytracing%pathsat(layer, j)
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(6, layer, i) * 2 * predictors%mixedgas(5, layer, j)
        tr_k(i, layer)                 = tr_k(i, layer) + predictors_k%mixedgas(5, layer, i)
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(4, layer, i) * 2 * predictors%mixedgas(3, layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, j) ** 2
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(3, layer, i) * raytracing%pathsat(layer, j)
      
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(3, layer, i) * tr(j, layer)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(2, layer, i) * 2 * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(1, layer, i)

        Predictors_k%mixedgas(1:21, layer, i) = 0.0_JPRB
      ENDDO

    ENDDO

  ELSE ! no Zeeman effect
  ! all sensors, except when running SSMIS with Zeeman coefficient file

  ! X14 -> X11
    IF (coef%id_inst == 3 .AND. coef%IncZeeman) THEN
    ! AMSU-A with Zeeman coefficient file

      DO layer = 1, profiles_k(1)%nlayers
        level = layer + 1

        DO i = 1, nchannels
          j = chanprof(i)%prof

  ! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
          raytracing_k%pathsat(layer, i)         = raytracing_k%pathsat(layer, i) +                              &
            & 2.0_JPRB * (profiles(j)%cosbk * profiles(j)%Be) ** 2 * raytracing%pathsat(layer, j) *              &
            & predictors_k%mixedgas(14, layer, i) + profiles(j)%Be ** 3 * predictors_k%mixedgas(13, layer, i) +  &
            & 2.0_JPRB * profiles(j)%Be * raytracing%pathsat(layer, j) * predictors_k%mixedgas(12, layer, i) +   &
            & profiles(j)%cosbk ** 2 * predictors_k%mixedgas(11, layer, i)
          predictors_k%mixedgas(11:14, layer, i) = 0.0_JPRB
  ! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
        ENDDO

      ENDDO

    ENDIF

  ! X10 -  NB tw(i,1) set to zero in direct
    DO layer = 2, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

  ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%mixedgas(10, layer, i) * 0.5_JPRB * tw(j, layer) ** 0.25_JPRB / predictors%mixedgas(9, layer, j)
        tw_k(i, layer)                 = tw_k(i, layer) +      &
          & predictors_k%mixedgas(10, layer, i) * 0.25_JPRB * predictors%mixedgas(9, layer, j) / tw(j, layer) ** 0.75
      ENDDO

    ENDDO

  ! X9 -> X1...

    DO layer = 1, profiles_k(1)%nlayers
      level = layer + 1

      DO i = 1, nchannels
        j = chanprof(i)%prof

  ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) +      &
          & predictors_k%mixedgas(9, layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(8, layer, i) * tw(j, layer) / tr(j, layer)
        tw_k(i, layer)                 =      &
          & tw_k(i, layer) + predictors_k%mixedgas(8, layer, i) * raytracing%pathsat(layer, j) / tr(j, layer)
        tr_k(i, layer)                 = tr_k(i, layer) -      &
          & predictors_k%mixedgas(8, layer, i) * tw(j, layer) * raytracing%pathsat(layer, j) / tr(j, layer) ** 2._JPRB
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(7, layer, i) * tw(j, layer)
        tw_k(i, layer)                 =      &
          & tw_k(i, layer) + predictors_k%mixedgas(7, layer, i) * raytracing%pathsat(layer, j)
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(6, layer, i) * 2 * predictors%mixedgas(5, layer, j)
        tr_k(i, layer)                 = tr_k(i, layer) + predictors_k%mixedgas(5, layer, i)
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(4, layer, i) * 2 * predictors%mixedgas(3, layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(4, layer, i) * predictors%mixedgas(5, layer, j) ** 2
        tr_k(i, layer)                 =      &
          & tr_k(i, layer) + predictors_k%mixedgas(3, layer, i) * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(3, layer, i) * tr(j, layer)
        raytracing_k%pathsat(layer, i) =      &
          & raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(2, layer, i) * 2 * raytracing%pathsat(layer, j)
        raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + predictors_k%mixedgas(1, layer, i)
      ENDDO

    ENDDO

  ENDIF

  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
  
    DO layer = profiles_k(1)%nlayers, 1,  - 1
  
      DO i = 1, nchannels
        j = chanprof(i)%prof
  
        IF (opts%ozone_Data) THEN
          sum1(i)       = sum1(i) + ow_k(i, layer) / sum2_ow(j, layer)
          o_k(i, layer) = o_k(i, layer) + sum1(i) * coef%dpp(layer - 1)
        ELSE
          o_k(i, layer) = 0._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF
  
  sum1 = 0._JPRB

  DO layer = profiles_k(1)%nlayers, 1,  - 1

    DO i = 1, nchannels
      j = chanprof(i)%prof
      sum1(i)       = sum1(i) + ww_k(i, layer) / sum2_ww(j, layer)
      w_k(i, layer) = w_k(i, layer) + sum1(i) * coef%dpp(layer - 1)
    ENDDO

  ENDDO


  DO layer = profiles_k(1)%nlayers, 2,  - 1

    DO i = 1, nchannels
      tw_k(i, layer - 1) = tw_k(i, layer - 1) + tw_k(i, layer)
      tr_k(i, layer - 1) = tr_k(i, layer - 1) + tw_k(i, layer) * coef%dpp(layer - 1)
    ENDDO

  ENDDO
!CDIR NOLOOPCHG

  IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
    DO layer = 1, profiles(1)%nlayers
  
      DO i = 1, nchannels
        j = chanprof(i)%prof

        o_k(i, layer) = o_k(i, layer) + or_k(i, layer) / coef%ostar(layer)
  
      ENDDO
  
    ENDDO
  ENDIF
  
!CDIR NOLOOPCHG

  DO layer = 1, profiles(1)%nlayers

    DO i = 1, nchannels
      w_k(i, layer) = w_k(i, layer) + wr_k(i, layer) / coef%wstar(layer)
      t_k(i, layer) = t_k(i, layer) + tr_k(i, layer) / coef%tstar(layer)
    ENDDO

  ENDDO

  DO layer = 1, profiles(1)%nlayers

    DO i = 1, nchannels
      j = chanprof(i)%prof

      IF (coef%nozone > 0) t_k(i, layer) = t_k(i, layer) + dto_k(i, layer)

      t_k(i, layer) = t_k(i, layer) + dt_k(i, layer)
    ENDDO

  ENDDO

  IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
    DO level = 2, profiles(1)%nlevels
      layer = level - 1
  
      DO i = 1, nchannels
        j = chanprof(i)%prof
        
        profiles_k(i)%o3(level - 1) = profiles_k(i)%o3(level - 1) + 0.5_JPRB * o_k(i, layer)
        profiles_k(i)%o3(level)     = profiles_k(i)%o3(level) + 0.5_JPRB * o_k(i, layer)
  
      ENDDO
  
    ENDDO
  ENDIF

!
! include (or not) adjoint surface humidity
!

  IF (opts%use_q2m) THEN

    DO i = 1, nchannels
      j = chanprof(i)%prof
! include surface humidity
      iv2lev = aux_prof%s(j)%nearestlev_surf
      iv2lay = iv2lev - 1
      profiles_k(i)%s2m%q        = profiles_k(i)%s2m%q + 0.5_JPRB * w_k(i, iv2lay)
      profiles_k(i)%q(iv3lev(j)) = profiles_k(i)%q(iv3lev(j)) + 0.5_JPRB * w_k(i, iv2lay)
    ENDDO


    DO i = 1, nchannels
      j = chanprof(i)%prof

      DO level = 2, iv3lev(j)
        layer = level - 1
        profiles_k(i)%q(level - 1) = profiles_k(i)%q(level - 1) + 0.5_JPRB * w_k(i, layer)
        profiles_k(i)%q(level)     = profiles_k(i)%q(level) + 0.5_JPRB * w_k(i, layer)
      ENDDO

    ENDDO

  ELSE

    DO level = 2, profiles_k(1)%nlevels
      layer = level - 1

      DO i = 1, nchannels
        profiles_k(i)%q(level - 1) = profiles_k(i)%q(level - 1) + 0.5_JPRB * w_k(i, layer)
        profiles_k(i)%q(level)     = profiles_k(i)%q(level) + 0.5_JPRB * w_k(i, layer)
      ENDDO

    ENDDO

  ENDIF


  DO level = 2, profiles_k(1)%nlevels
    layer = level - 1

    DO i = 1, nchannels
      profiles_k(i)%t(level - 1) = profiles_k(i)%t(level - 1) + 0.5_JPRB * t_k(i, layer)
      profiles_k(i)%t(level)     = profiles_k(i)%t(level) + 0.5_JPRB * t_k(i, layer)
    ENDDO

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_k
