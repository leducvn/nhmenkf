SUBROUTINE rttov_transmit_9_solar_k( &
            & addaerosl,                      &
            & addclouds,                      &
            & nlayers,                        &
            & chanprof,                       &
            & profiles,                       &
            & sun,                            &
            & aux,                            &
            & aux_k,                          &
            & coef,                           &
            & raytracing,                     &
            & raytracing_k,                   &
            & ircld,                          &
            & opdp_path,                      &
            & opdp_path_k,                    &
            & odsun_level,                    &
            & odsun_singlelayer,              &
            & od_frac,                        &
            & transmission_aux,               &
            & transmission_aux_k,             &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & tausun_ref,                     &
            & tausun_ref_surf,                &
            & tausun_level,                   &
            & tausun_surf)
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
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
!                     has been added to the total transmission_aux.
!  3.0    06/02/2007  removed polarisation R Saunders
!  4.0    15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  5.0    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  6.0    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
!  6.1    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  6.2    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  6.3    14/12/2010  Use traj0_sta%sun array to flag channels for which solar calculations
!                     should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Adjoint variables
! input transmission_aux_k% tau_surf and transmission_aux_k% tau_level set inside integrate_k
!
! input/output aux_k
!
! output predictors_k initialised inside rttov_k (need input
!    intent for memory allocation in calling routine)
!
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
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_k
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_k
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing_k
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_level      (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_singlelayer(nlayers       , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref       (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref_surf  (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_level     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_surf      (size(chanprof)                )
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_k(nlayers       , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_k    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_level_k      (nlayers + 1   , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_k       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_frac_k       (size(chanprof)                )
  REAL   (KIND=jprb) :: tausun_level_k  (nlayers + 1   , size(chanprof))  ! sat to level transmission_aux at each frequency
  REAL   (KIND=jprb) :: tausun_surf_k   (size(chanprof)                )
  INTEGER(KIND=jpim) :: lev         , lay    , chan   , j, ist, nlevels
  INTEGER(KIND=jpim) :: prof        , levsurf
  INTEGER(KIND=jpim) :: nchannels                                         ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels             = size(chanprof)
  nlevels               = nlayers + 1
!---------------------------------------------------------------------------------------
!K of store transmittances for other polarisations
!---------------------------------------------------------------------------------------
  od_level_k(:,:)       = 0._JPRB
  tausun_level_k(:,:)   = 0.0_JPRB
  od_singlelayer_k(:,:) = 0.0_JPRB
  tausun_surf_k(:)      = 0.0_JPRB
  od_frac_k(:)          = 0.0_JPRB
  od_surf_k(:)          = 0.0_JPRB
  od_surf_ac_k(:)       = 0.0_JPRB
  DO j = nchannels, 1,  - 1
    prof = chanprof(j)%prof! Profile index
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      DO ist = ircld%nstream(prof), 0,  - 1
        IF (addaerosl .OR. addclouds) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
            IF (transmission_scatt_ir_stream%opdpext(ist, j, lay) /= 0._jprb) THEN
              transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) =      &
                & transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) +  &
                & transmission_scatt_ir_stream_k%ssa(ist, j, lay) / transmission_scatt_ir_stream%opdpext(ist, j, lay)
              transmission_scatt_ir_stream_k%opdpext(ist, j, lay) =                                                      &
                & transmission_scatt_ir_stream_k%opdpext(ist, j, lay) -                                                  &
                & transmission_scatt_ir_stream_k%ssa(ist, j, lay) * transmission_scatt_ir_stream%opdpsca(ist, j, lay) /  &
                & transmission_scatt_ir_stream%opdpext(ist, j, lay) ** 2
            ENDIF
            od_singlelayer_k(lay, j)                               =      &
              & od_singlelayer_k(lay, j) + transmission_scatt_ir_stream_k%opdpext(ist, j, lay)
            transmission_scatt_ir_stream_k%opdpabs(ist, j, lay)    =      &
              & transmission_scatt_ir_stream_k%opdpabs(ist, j, lay) + transmission_scatt_ir_stream_k%opdpext(ist, j, lay)
            transmission_scatt_ir_stream_k%opdpsca(ist, j, lay)    =      &
              & transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) + transmission_scatt_ir_stream_k%opdpext(ist, j, lay)
            od_singlelayer_k(lay, j)                               = od_singlelayer_k(lay, j) +      &
              & transmission_aux_k%odsun_singlelayer(lay, ist, j) * raytracing%pathsun(lay, prof) /  &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) =                                 &
              & transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) +                             &
              & transmission_aux_k%odsun_singlelayer(lay, ist, j) * raytracing%pathsun(lay, prof) /  &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            raytracing_k%pathsat(lay, j)                           = raytracing_k%pathsat(lay, j) -      &
              & transmission_aux_k%odsun_singlelayer(lay, ist, j) *                                      &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *     &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2
            raytracing_k%pathsun(lay, j)                           = raytracing_k%pathsun(lay, j) +      &
              & transmission_aux_k%odsun_singlelayer(lay, ist, j) *                                      &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) /     &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            raytracing_k%pathsun(lay, j)                           = raytracing_k%pathsun(lay, j) -      &
              & transmission_aux_k%odsun_singlelayer(lay, ist, j) *                                      &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *     &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2
          ENDDO
          transmission_aux_k%odsun_frac_t(ist, J)  =      &
            & transmission_aux_k%odsun_frac_t(ist, J) + transmission_aux_k%odsun_sfrac(ist, j)
          transmission_aux_k%tausun_surf_t(ist, j) =      &
            & transmission_aux_k%tausun_surf_t(ist, j) + transmission_aux_k%tausun_surf(ist, j)
          DO lev = nlevels, 1,  - 1
            IF (tausun_level(lev, j) >= 0) THEN
              tausun_level_k(lev, j)                                = tausun_level_k(lev, j) +      &
                & transmission_aux_k%tausun_level(lev, ist, j) *                                    &
                & exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
              transmission_scatt_ir_stream_k%opdpacsun(ist, j, lev) =                    &
                & transmission_scatt_ir_stream_k%opdpacsun(ist, j, lev) -                &
                & transmission_aux_k%tausun_level(lev, ist, j) * tausun_level(lev, j) *  &
                & exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
            ELSE
              tausun_level_k(lev, j) = tausun_level_k(lev, j) + transmission_aux_k%tausun_level(lev, ist, j)
            ENDIF
          ENDDO
        ELSE
          tausun_level_k(:, j) = tausun_level_k(:, j) + transmission_aux_k%tausun_level(:, ist, j)
          tausun_surf_k(j)     = tausun_surf_k(j) + transmission_aux_k%tausun_surf(ist, j)
        ENDIF
      ENDDO
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_k%tausun_level(:, ist, j)      = 0.0_JPRB
        transmission_aux_k%tausun_surf(ist, j)          = 0.0_JPRB
        transmission_aux_k%odsun_singlelayer(:, ist, j) = 0.0_JPRB
        transmission_aux_k%odsun_sfrac(ist, j)          = 0.0_JPRB
      ENDDO
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!K of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof           ! Profile index
    chan    = chanprof(j)%chan
    levsurf = aux%s(prof)%nearestlev_surf
    IF (sun(j)) THEN
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = ircld%nstream(prof), 0,  - 1
          IF (tausun_surf(j) >= 0) THEN
            tausun_surf_k(j)                          =      &
              & tausun_surf_k(j) + transmission_aux_k%tausun_surf_t(ist, j) * transmission_aux%tau_surf_acsun(ist, J)
            transmission_aux_k%tau_surf_acsun(ist, J) =      &
              & transmission_aux_k%tau_surf_acsun(ist, J) + transmission_aux_k%tausun_surf_t(ist, j) * tausun_surf(j)
          ELSE
            tausun_surf_k(j) = tausun_surf_k(j) + transmission_aux_k%tausun_surf_t(ist, j)
          ENDIF
          od_frac_k(j)                             = od_frac_k(j) -                              &
            & transmission_aux_k%odsun_frac_t(ist, J) * raytracing%pathsun(levsurf - 1, prof) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          transmission_aux_k%odsun_frac_ac(ist, J) = transmission_aux_k%odsun_frac_ac(ist, J) +      &
            & transmission_aux_k%odsun_frac_t(ist, J) * raytracing%pathsun(levsurf - 1, prof) /      &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          raytracing_k%pathsat(levsurf - 1, j)     = raytracing_k%pathsat(levsurf - 1, j) -                         &
            & transmission_aux_k%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *  &
            & raytracing%pathsun(levsurf - 1, prof) /                                                               &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2
          raytracing_k%pathsun(levsurf - 1, j)     = raytracing_k%pathsun(levsurf - 1, j) +                         &
            & transmission_aux_k%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          raytracing_k%pathsun(levsurf - 1, j)     = raytracing_k%pathsun(levsurf - 1, j) -                         &
            & transmission_aux_k%odsun_frac_t(ist, J) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *  &
            & raytracing%pathsun(levsurf - 1, prof) /                                                               &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2
          od_surf_ac_k(j)                          =      &
            & od_surf_ac_k(j) - transmission_aux_k%tau_surf_acsun(ist, j) * transmission_aux%tau_ref_surf_acsun(ist, j)
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_k(j) = od_surf_ac_k(j) + transmission_aux_k%odsun_frac_ac(ist, J)
            transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf - 1) =      &
              & transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf - 1) - transmission_aux_k%odsun_frac_ac(ist, J)
          ELSE
            od_surf_ac_k(j)                                           =      &
              & od_surf_ac_k(j) + transmission_aux_k%odsun_frac_ac(ist, J)
            transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf) =      &
              & transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf) - transmission_aux_k%odsun_frac_ac(ist, J)
          ENDIF
          transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf)     =      &
            & transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf) +      &
            & od_surf_ac_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf - 1) =      &
            & transmission_scatt_ir_stream_k%opdpacsun(ist, j, levsurf - 1) + od_surf_ac_k(j) * aux%s(prof)%pfraction_surf
          aux_k%s(j)%pfraction_surf                                     = aux_k%s(j)%pfraction_surf + od_surf_ac_k(j)     &
            &  * (transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf - 1) -                                           &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf))
          od_surf_ac_k(j) = 0._JPRB
        ENDDO
      ENDIF
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tausun_surf(j) < coef%tt_a0(chan)) THEN
          tausun_surf_k(j) = 0._jprb
        ENDIF
      ENDIF
      od_surf_k(j) = od_surf_k(j) + tausun_surf_k(j) * tausun_ref_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_k(j)               = od_surf_k(j) + od_frac_k(J)
        od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) - od_frac_k(J)
      ELSE
        od_surf_k(j)           = od_surf_k(j) + od_frac_k(J)
        od_level_k(levsurf, j) = od_level_k(levsurf, j) - od_frac_k(J)
      ENDIF
      od_level_k(levsurf, j)     = od_level_k(levsurf, j) + od_surf_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
      od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) + od_surf_k(j) * aux%s(prof)%pfraction_surf
      aux_k%s(j)%pfraction_surf  =      &
        & aux_k%s(j)%pfraction_surf + od_surf_k(j) * (odsun_level(levsurf - 1, j) - odsun_level(levsurf, j))
      od_surf_k(j)               = 0._JPRB
    ENDIF
  ENDDO
!-------------------------------------------
!K of assemble layer optical depths
!-------------------------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tausun_ref(lev, j) < coef%tt_a0(chan)) THEN
            tausun_level_k(lev, j) = 0._jprb
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level_k(lev, j) = od_level_k(lev, j) + tausun_level_k(lev, j) * tausun_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      od_singlelayer_k(:, j) = coef%ff_gam(chan) * od_singlelayer_k(:, j)
      od_level_k(:, j)       = coef%ff_gam(chan) * od_level_k(:, j)
    ENDIF
  ENDDO
! k of level to space optical depths
  opdp_path_k%sun_level(:,:) = opdp_path_k%sun_level(:,:) + od_level_k(:,:)
! k of single layer optical depths
  DO lay = nlayers, 1,  - 1
    opdp_path_k%sun_level(lay, :)     = opdp_path_k%sun_level(lay, :) + od_singlelayer_k(lay, :)
    opdp_path_k%sun_level(lay + 1, :) = opdp_path_k%sun_level(lay + 1, :) - od_singlelayer_k(lay, :)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_k
