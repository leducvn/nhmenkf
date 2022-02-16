!
SUBROUTINE rttov_setpredictors_8( &
            & opts,       &
            & prof,       &
            & geom,       &
            & coef,       &
            & aux,        &
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
!  1.0   29/01/2003  Original - copy of RTTOV7 model (P Brunel)
!  1.1   11/09/2003  Added predictors for wv line and continuum and CO2 (R Saunders)
!  1.2   03/06/2004  Parkind parametrisation (P. Brunel)
!  1.3   23/02/2005  Correction of Twr definition (P. Brunel)
!  1.4   29/03/2005  Add end of header comment (J. Cameron)
!  1.5   07/12/2005  Add surface humidity (R. Saunders)
!  1.6   17/01/2006  Marco Matricardi (ECMWF):
!           --       Introduced altitude dependent local zenith angle
!  1.7   14/03/2007  Corrected CO2 profile logic+added raytracing  (R Saunders)
!  1.8   27/02/2009  Profile levels to include ToA. Distinguish arrays
!                    in raytracing (on levels) from all others (on
!                    layers). Predictors prepared to maintain agreement
!                    with the RTTOV-8 scheme in RTTOV-9 (P. Rayer)
!  1.9   02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                    (Marco Matricardi).
!  1.10  10/11/2010  Remove rttov9_compat flag from code (J Hocking)
!  1.11  25/10/2011  Fix bug so that omitted trace gas profiles are treated 
!                    correctly (J Hocking)
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
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE parkind1, ONLY : JPRB
  USE parkind1, ONLY : JPRB
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))! geometry
  TYPE(predictors_Type), INTENT(INOUT) :: predictors      ! predictors
  TYPE(profile_aux    ), INTENT(IN)    :: aux             ! auxillary profiles info.
  TYPE(raytracing_Type), INTENT(INOUT) :: raytracing      !
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer , iprof
  INTEGER(KIND=jpim) :: iv2lev, iv3lev, iv2lay
! user profile
  REAL   (KIND=jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2  (prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers)
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w (prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: deltac  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_or  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_sq   (prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles                     ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!-------------------------------------------------------------------------------
! 1 profile layer quantities
!-------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO iprof = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) / 2._JPRB
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) / 2._JPRB
    ENDDO


    IF (opts%use_q2m) THEN
! include surface humidity
      iv3lev = aux%s(iprof)%nearestlev_surf - 1
      iv2lev = aux%s(iprof)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay    = iv2lev - 1
        w(iv2lay) = (prof(iprof)%s2m%q + prof(iprof)%q(iv3lev)) / 2._JPRB
      ENDIF

    ENDIF

!

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level    = layer + 1
        o(layer) = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) / 2._JPRB
      ENDDO

    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level      = layer + 1
        co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) / 2._JPRB
      ENDDO

    ENDIF


!------------------------------------------------------------------------------
! 2 calculate deviations from reference profile (layers)
! if no input O3 profile we still use the input temperature profile for dto
!-----------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:)    = t(:) - coef%tstar(:)
    dtabs(:) = Abs(dt(:))
    IF (coef%nozone > 0) dto(:) = t(:) - coef%to3star(:)

!------------------------------------------------------------------------------
! 3 calculate (profile / reference profile) ratios; tr wr or co2r
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:) = t(:) / coef%tstar(:)
    wr(:) = w(:) / coef%wstar(:)
! if no input O3 profile, set to reference value (or =1)

    IF (coef%nozone > 0) THEN
      IF (opts%ozone_Data) THEN
        or(:) = o(:) / coef%ostar(:)
      ELSE
        or(:) = 1._JPRB
      ENDIF
    ENDIF
    
! if no input CO2 profile, set to reference value (co2r=1)

    IF (coef%nco2 > 0) THEN
      IF (opts%co2_Data) THEN
        co2r(:) = co2(:) / coef%co2star(:)
      ELSE
        co2r(:) = 1._JPRB
      ENDIF
    ENDIF
    
!------------------------------------------------------------------------------
! 4 calculate profile / reference profile sums: tw ww ow co2w twr
!------------------------------------------------------------------------------
    tw(1) = 0._JPRB

    DO layer = 2, prof(1)%nlayers
! cumulate overlying layers (tr relates to same layer as dpp)
! do not need dpp(0) to start
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB
! cumulating overlying layer and the layer itself

    DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (w, wstar relate to layer below dpp)
! need dpp(0) to start
      sum1      = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2      = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(layer) = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1       = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2       = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      wwr(layer) = sum1 / sum2
    ENDDO

! if no input O3 profile, set to reference value (ow =1)

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (o, ostar relate to layer below dpp)
! need dpp(0) to start
        sum1      = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2      = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow(:) = 1._JPRB
    ENDIF

! if no input co2 profile, set to reference value (co2w=1) 
! but twr calculation still uses the input temperature profile

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
! cumulate overlying layers (co2, co2star relate to layer below dpp)
! need dpp(0) to start
        sum1        = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2        = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w(:) = 1._JPRB
    ENDIF
    
    ! twr only used in CO2 predictors
    IF (coef%nco2 > 0) THEN
      sum1   = 0._JPRB
      sum2   = 0._JPRB
      twr(1) = 1._JPRB
  
      DO layer = 2, prof(1)%nlayers
! cumulate overlying layers (t, tstar relate to same layer as dpp)
! do not need dpp(0) to start
        sum1       = sum1 + coef%dpp(layer - 1) * t(layer - 1)
        sum2       = sum2 + coef%dpp(layer - 1) * coef%tstar(layer - 1)
        twr(layer) = sum1 / sum2
      ENDDO
    ENDIF

!5) set predictors for RTTOV-8 options
!--
!5.1 mixed gases
!---

    DO layer = 1, prof(1)%nlayers
      level = layer + 1! raytracing % pathsat(layer,i) is angle at lower boundary of layer
      tr_sq(layer)                             = tr(layer) * tr(layer)
      predictors%mixedgas(1, layer, iprof)     = raytracing%pathsat(layer, iprof)
      predictors%mixedgas(2, layer, iprof)     = raytracing%pathsat(layer, iprof) * raytracing%pathsat(layer, iprof)
      predictors%mixedgas(3, layer, iprof)     = raytracing%pathsat(layer, iprof) * tr(layer)
      predictors%mixedgas(4, layer, iprof)     = raytracing%pathsat(layer, iprof) * tr_sq(layer)
      predictors%mixedgas(5, layer, iprof)     = tr(layer)
      predictors%mixedgas(6, layer, iprof)     = tr_sq(layer)
      predictors%mixedgas(7, layer, iprof)     = raytracing%pathsat(layer, iprof) * tw(layer)
      predictors%mixedgas(8, layer, iprof)     = raytracing%pathsat(layer, iprof) * tw(layer) / tr(layer)
      predictors%mixedgas(9, layer, iprof)     = raytracing%pathsat(layer, iprof) ** 0.5_JPRB
      predictors%mixedgas(10, layer, iprof)    = raytracing%pathsat(layer, iprof) ** 0.5_JPRB * tw(layer) ** 0.25_JPRB
!5.2 water vapour line transmittance based on RTIASI but with pred 9 removed
!----------------
      sec_wr(layer)                            = raytracing%pathsat(layer, iprof) * wr(layer)
      sec_wrwr(layer)                          = sec_wr(layer) * wr(layer)
!predictors % watervapour(:,:,iprof) = 0._JPRB
      predictors%watervapour(1, layer, iprof)  = sec_wr(layer) * sec_wr(layer)
      predictors%watervapour(2, layer, iprof)  = raytracing%pathsat(layer, iprof) * ww(layer)
      predictors%watervapour(3, layer, iprof)  = (raytracing%pathsat(layer, iprof) * ww(layer)) ** 2
      predictors%watervapour(4, layer, iprof)  = sec_wr(layer) * dt(layer)
      predictors%watervapour(5, layer, iprof)  = Sqrt(sec_wr(layer))
      predictors%watervapour(6, layer, iprof)  = sec_wr(layer) ** 0.25_JPRB
      predictors%watervapour(7, layer, iprof)  = sec_wr(layer)
      predictors%watervapour(8, layer, iprof)  = sec_wr(layer) ** 3
      predictors%watervapour(9, layer, iprof)  = sec_wr(layer) * dt(layer) * dtabs(layer)
      predictors%watervapour(10, layer, iprof) = Sqrt(sec_wr(layer)) * dt(layer)
      predictors%watervapour(11, layer, iprof) = sec_wrwr(layer) / wwr(layer)
      predictors%watervapour(12, layer, iprof) =      &
        & Sqrt(raytracing%pathsat(layer, iprof)) * wr(layer) ** 1.5_JPRB / wwr(layer)
!5.3 water vapour continuum transmittance based on RTIASI
!----------------
!

      IF (coef%nwvcont > 0) THEN
!predictors % wvcont(:,:,iprof)  = 0._JPRB
        predictors%wvcont(1, layer, iprof) = sec_wrwr(layer) / tr(layer)
        predictors%wvcont(2, layer, iprof) = sec_wrwr(layer) / (tr_sq(layer) * tr_sq(layer))
        predictors%wvcont(3, layer, iprof) = sec_wr(layer) / tr(layer)
        predictors%wvcont(4, layer, iprof) = sec_wr(layer) / tr_sq(layer)
      ENDIF

!5.4 ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

      IF (coef%nozone > 0) THEN
        sec_or(layer)                      = or(layer) * raytracing%pathsat(layer, iprof)
        predictors%ozone(1, layer, iprof)  = sec_or(layer)
        predictors%ozone(2, layer, iprof)  = Sqrt(sec_or(layer))
        predictors%ozone(3, layer, iprof)  = sec_or(layer) * dto(layer)
        predictors%ozone(4, layer, iprof)  = sec_or(layer) * sec_or(layer)
        predictors%ozone(5, layer, iprof)  = Sqrt(sec_or(layer)) * dto(layer)
        predictors%ozone(6, layer, iprof)  = sec_or(layer) * or(layer) * ow(layer)
        predictors%ozone(7, layer, iprof)  = Sqrt(sec_or(layer)) * or(layer) / ow(layer)
        predictors%ozone(8, layer, iprof)  = sec_or(layer) * ow(layer)
        predictors%ozone(9, layer, iprof)  = sec_or(layer) * Sqrt(raytracing%pathsat(layer, iprof) * ow(layer))
        predictors%ozone(10, layer, iprof) = raytracing%pathsat(layer, iprof) * ow(layer)
        predictors%ozone(11, layer, iprof) =      &
          & raytracing%pathsat(layer, iprof) * ow(layer) * raytracing%pathsat(layer, iprof) * ow(layer)
      ENDIF

    ENDDO

!5.5 cloud
!---------

    IF (opts%clw_Data .AND. coef%id_sensor == sensor_id_mw) THEN

      DO layer = 1, prof(1)%nlayers
! NB in RTTOV-10, layer 1 is the layer above level 2
        level = layer + 1
        deltac(layer)                = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)
        predictors%clw(layer, iprof) = deltac(layer) * prof(iprof)%clw(level) * geom(iprof)%seczen
      ENDDO


      DO layer = 2, prof(1)%nlayers
        level = layer + 1
        predictors%clw(layer, iprof) =      &
          & 0.5_JPRB * (predictors%clw(layer, iprof) + deltac(layer) * prof(iprof)%clw(level - 1) * geom(iprof)%seczen)
      ENDDO

      predictors%ncloud = 1
    ELSE
      predictors%clw = 0._jprb
    ENDIF

!5.6 carbon dioxide transmittance based on RTIASI
!-------------------------------------------------
!

    IF (coef%nco2 > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level = layer + 1
        predictors%co2(1, layer, iprof)  = raytracing%pathsat(layer, iprof) * co2r(layer)
        predictors%co2(2, layer, iprof)  = tr_sq(layer)
        predictors%co2(3, layer, iprof)  = raytracing%pathsat(layer, iprof) * tr(layer)
        predictors%co2(4, layer, iprof)  = raytracing%pathsat(layer, iprof) * tr_sq(layer)
        predictors%co2(5, layer, iprof)  = tr(layer)
        predictors%co2(6, layer, iprof)  = raytracing%pathsat(layer, iprof)
        predictors%co2(7, layer, iprof)  = raytracing%pathsat(layer, iprof) * twr(layer)
        predictors%co2(8, layer, iprof)  = (raytracing%pathsat(layer, iprof) * co2w(layer)) ** 2
        predictors%co2(9, layer, iprof)  = twr(layer) * twr(layer) * twr(layer)
        predictors%co2(10, layer, iprof) = raytracing%pathsat(layer, iprof) * twr(layer) * Sqrt(tr(layer))
      ENDDO

    ENDIF

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_8', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_8
