SUBROUTINE rttov_transmit_9_solar_ad( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & chanprof,                        &
            & profiles,                        &
            & sun,                             &
            & aux,                             &
            & aux_ad,                          &
            & coef,                            &
            & raytracing,                      &
            & raytracing_ad,                   &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_ad,                    &
            & odsun_level,                     &
            & odsun_singlelayer,               &
            & od_frac,                         &
            & transmission_aux,                &
            & transmission_aux_ad,             &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_ad, &
            & tausun_ref,                      &
            & tausun_ref_surf,                 &
            & tausun_level,                    &
            & tausun_surf)
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
! input transmission_aux_ad % tau_surf and transmission_aux_ad % tau_level
! set inside integrate_ad
!
! input/output aux_ad
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
!  1.3    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.4    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
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
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & profile_type,               &
       & ircld_type,                 &
       & raytracing_type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_const, ONLY : max_sol_zen
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type              ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_ad
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_ad
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_ad
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing_ad
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_ad
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_level      (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_singlelayer(nlayers       , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref       (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref_surf  (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_level     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_surf      (size(chanprof)                )
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_ad(nlayers       , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_ad    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_level_ad      (nlayers + 1   , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_ad       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_frac_ad       (size(chanprof)                )
  REAL   (KIND=jprb) :: tausun_level_ad  (nlayers + 1   , size(chanprof))            ! sat to level transmission_aux at each frequency
  REAL   (KIND=jprb) :: tausun_surf_ad   (size(chanprof)                )
  INTEGER(KIND=jpim) :: lev         , lay, chan, j, ist, levsurf, nlevels
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: nchannels                                                    ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels              = size(chanprof)
  nlevels                = nlayers + 1
!---------------------------------------------------------------------------------------
!AD of store transmittances for other polarisations
!---------------------------------------------------------------------------------------
  od_level_ad(:,:)       = 0._JPRB
  tausun_level_ad(:,:)   = 0.0_JPRB
  od_singlelayer_ad(:,:) = 0.0_JPRB
  tausun_surf_ad(:)      = 0.0_JPRB
  od_frac_ad(:)          = 0.0_JPRB
  od_surf_ad(:)          = 0.0_JPRB
  od_surf_ac_ad(:)       = 0.0_JPRB
  DO j = nchannels, 1,  - 1
    prof = chanprof(j)%prof! Profile index
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      DO ist = ircld%nstream(prof), 0,  - 1
        IF (addaerosl .OR. addclouds) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
            IF (transmission_scatt_ir_stream%opdpext(ist, j, lay) /= 0._jprb) THEN
              transmission_scatt_ir_stream_ad%opdpsca(ist, j, lay) =      &
                & transmission_scatt_ir_stream_ad%opdpsca(ist, j, lay) +  &
                & transmission_scatt_ir_stream_ad%ssa(ist, j, lay) / transmission_scatt_ir_stream%opdpext(ist, j, lay)
              transmission_scatt_ir_stream_ad%opdpext(ist, j, lay) =                                                      &
                & transmission_scatt_ir_stream_ad%opdpext(ist, j, lay) -                                                  &
                & transmission_scatt_ir_stream_ad%ssa(ist, j, lay) * transmission_scatt_ir_stream%opdpsca(ist, j, lay) /  &
                & transmission_scatt_ir_stream%opdpext(ist, j, lay) ** 2
            ENDIF
            od_singlelayer_ad(lay, j)                               =      &
              & od_singlelayer_ad(lay, j) + transmission_scatt_ir_stream_ad%opdpext(ist, j, lay)
            transmission_scatt_ir_stream_ad%opdpabs(ist, j, lay)    =      &
              & transmission_scatt_ir_stream_ad%opdpabs(ist, j, lay) +     &
              & transmission_scatt_ir_stream_ad%opdpext(ist, j, lay)
            transmission_scatt_ir_stream_ad%opdpsca(ist, j, lay)    =      &
              & transmission_scatt_ir_stream_ad%opdpsca(ist, j, lay) +     &
              & transmission_scatt_ir_stream_ad%opdpext(ist, j, lay)
            od_singlelayer_ad(lay, j)                               = od_singlelayer_ad(lay, j) +      &
              & transmission_aux_ad%odsun_singlelayer(lay, ist, j) * raytracing%pathsun(lay, prof) /   &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            transmission_scatt_ir_stream_ad%opdpaclsun(ist, j, lay) =                                 &
              & transmission_scatt_ir_stream_ad%opdpaclsun(ist, j, lay) +                             &
              & transmission_aux_ad%odsun_singlelayer(lay, ist, j) * raytracing%pathsun(lay, prof) /  &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            raytracing_ad%pathsat(lay, prof)                        = raytracing_ad%pathsat(lay, prof) -      &
              & transmission_aux_ad%odsun_singlelayer(lay, ist, j) *                                          &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *          &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2
            raytracing_ad%pathsun(lay, prof)                        = raytracing_ad%pathsun(lay, prof) +      &
              & transmission_aux_ad%odsun_singlelayer(lay, ist, j) *                                          &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) /          &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            raytracing_ad%pathsun(lay, prof)                        = raytracing_ad%pathsun(lay, prof) -      &
              & transmission_aux_ad%odsun_singlelayer(lay, ist, j) *                                          &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *          &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2
          ENDDO
          transmission_aux_ad%odsun_frac_t(ist, J)  =      &
            & transmission_aux_ad%odsun_frac_t(ist, J) + transmission_aux_ad%odsun_sfrac(ist, j)
          transmission_aux_ad%tausun_surf_t(ist, j) =      &
            & transmission_aux_ad%tausun_surf_t(ist, j) + transmission_aux_ad%tausun_surf(ist, j)
          DO lev = nlevels, 1,  - 1
            IF (tausun_level(lev, j) >= 0) THEN
              tausun_level_ad(lev, j)                                = tausun_level_ad(lev, j) +      &
                & transmission_aux_ad%tausun_level(lev, ist, j) *                                     &
                & exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
              transmission_scatt_ir_stream_ad%opdpacsun(ist, j, lev) =                    &
                & transmission_scatt_ir_stream_ad%opdpacsun(ist, j, lev) -                &
                & transmission_aux_ad%tausun_level(lev, ist, j) * tausun_level(lev, j) *  &
                & exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
            ELSE
              tausun_level_ad(lev, j) = tausun_level_ad(lev, j) + transmission_aux_ad%tausun_level(lev, ist, j)
            ENDIF
          ENDDO
        ELSE
          tausun_level_ad(:, j) = tausun_level_ad(:, j) + transmission_aux_ad%tausun_level(:, ist, j)
          tausun_surf_ad(j)     = tausun_surf_ad(j) + transmission_aux_ad%tausun_surf(ist, j)
        ENDIF
      ENDDO
!     ENDIF
!     IF (profiles(prof)%sunzenangle >= 0.0 .AND. &
!         profiles(prof)%sunzenangle < max_sol_zen) THEN
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_ad%tausun_level(:, ist, j)      = 0.0_JPRB
        transmission_aux_ad%tausun_surf(ist, j)          = 0.0_JPRB
        transmission_aux_ad%odsun_singlelayer(:, ist, j) = 0.0_JPRB
        transmission_aux_ad%odsun_sfrac(ist, j)          = 0.0_JPRB
      ENDDO
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!AD of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (sun(j)) THEN
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = ircld%nstream(prof), 0,  - 1
          IF (tausun_surf(j) >= 0) THEN
            tausun_surf_ad(j)                          =      &
              & tausun_surf_ad(j) + transmission_aux_ad%tausun_surf_t(ist, j) * transmission_aux%tau_surf_acsun(ist, J)
            transmission_aux_ad%tau_surf_acsun(ist, J) =      &
              & transmission_aux_ad%tau_surf_acsun(ist, J) + transmission_aux_ad%tausun_surf_t(ist, j) * tausun_surf(j)
          ELSE
            tausun_surf_ad(j) = tausun_surf_ad(j) + transmission_aux_ad%tausun_surf_t(ist, j)
          ENDIF
          od_frac_ad(j)                             = od_frac_ad(j) -                             &
            & transmission_aux_ad%odsun_frac_t(ist, J) * raytracing%pathsun(levsurf - 1, prof) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          transmission_aux_ad%odsun_frac_ac(ist, J) = transmission_aux_ad%odsun_frac_ac(ist, J) +      &
            & transmission_aux_ad%odsun_frac_t(ist, J) * raytracing%pathsun(levsurf - 1, prof) /       &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          raytracing_ad%pathsat(levsurf - 1, prof)  = raytracing_ad%pathsat(levsurf - 1, prof) -                     &
            & transmission_aux_ad%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *  &
            & raytracing%pathsun(levsurf - 1, prof) /                                                                &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2
          raytracing_ad%pathsun(levsurf - 1, prof)  = raytracing_ad%pathsun(levsurf - 1, prof) +                     &
            & transmission_aux_ad%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          raytracing_ad%pathsun(levsurf - 1, prof)  = raytracing_ad%pathsun(levsurf - 1, prof) -                     &
            & transmission_aux_ad%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *  &
            & raytracing%pathsun(levsurf - 1, prof) /                                                                &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2
          od_surf_ac_ad(j)                          =      &
            & od_surf_ac_ad(j) - transmission_aux_ad%tau_surf_acsun(ist, j) * transmission_aux%tau_ref_surf_acsun(ist, j)
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_ad(j)                                               =      &
              & od_surf_ac_ad(j) + transmission_aux_ad%odsun_frac_ac(ist, J)
            transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf - 1) =      &
              & transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf - 1) - transmission_aux_ad%odsun_frac_ac(ist, J)
          ELSE
            od_surf_ac_ad(j)                                           =      &
              & od_surf_ac_ad(j) + transmission_aux_ad%odsun_frac_ac(ist, J)
            transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf) =      &
              & transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf) - transmission_aux_ad%odsun_frac_ac(ist, J)
          ENDIF
          transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf)     =      &
            & transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf) +      &
            & od_surf_ac_ad(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf - 1) =      &
            & transmission_scatt_ir_stream_ad%opdpacsun(ist, j, levsurf - 1) +  &
            & od_surf_ac_ad(j) * aux%s(prof)%pfraction_surf
          aux_ad%s(prof)%pfraction_surf                                  = aux_ad%s(prof)%pfraction_surf +      &
            & od_surf_ac_ad(j) * (transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf - 1) -                 &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf))
          od_surf_ac_ad(j)                                               = 0._JPRB
        ENDDO
      ENDIF
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tausun_surf(j) < coef%tt_a0(chan)) THEN
          tausun_surf_ad(j) = 0._jprb
        ENDIF
      ENDIF
      od_surf_ad(j) = od_surf_ad(j) + tausun_surf_ad(j) * tausun_ref_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_ad(j)               = od_surf_ad(j) + od_frac_ad(J)
        od_level_ad(levsurf - 1, j) = od_level_ad(levsurf - 1, j) - od_frac_ad(J)
      ELSE
        od_surf_ad(j)           = od_surf_ad(j) + od_frac_ad(J)
        od_level_ad(levsurf, j) = od_level_ad(levsurf, j) - od_frac_ad(J)
      ENDIF
      od_level_ad(levsurf, j)       = od_level_ad(levsurf, j) + od_surf_ad(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
      od_level_ad(levsurf - 1, j)   = od_level_ad(levsurf - 1, j) + od_surf_ad(j) * aux%s(prof)%pfraction_surf
      aux_ad%s(prof)%pfraction_surf =      &
        & aux_ad%s(prof)%pfraction_surf + od_surf_ad(j) * (odsun_level(levsurf - 1, j) - odsun_level(levsurf, j))
      od_surf_ad(j)                 = 0._JPRB
    ENDIF
  ENDDO
!-------------------------------------------
!AD of assemble layer optical depths
!-------------------------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tausun_ref(lev, j) < coef%tt_a0(chan)) THEN
            tausun_level_ad(lev, j) = 0._jprb
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level_ad(lev, j) = od_level_ad(lev, j) + tausun_level_ad(lev, j) * tausun_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      od_singlelayer_ad(:, j) = coef%ff_gam(chan) * od_singlelayer_ad(:, j)
      od_level_ad(:, j)       = coef%ff_gam(chan) * od_level_ad(:, j)
    ENDIF
  ENDDO
! ad of level to space optical depths
  opdp_path_ad%sun_level(:,:) = opdp_path_ad%sun_level(:,:) + od_level_ad(:,:)
! ad of single layer optical depths
  DO lay = nlayers, 1,  - 1
    opdp_path_ad%sun_level(lay, :)     = opdp_path_ad%sun_level(lay, :) + od_singlelayer_ad(lay, :)
    opdp_path_ad%sun_level(lay + 1, :) = opdp_path_ad%sun_level(lay + 1, :) - od_singlelayer_ad(lay, :)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_ad
