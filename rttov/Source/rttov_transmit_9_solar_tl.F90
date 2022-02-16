SUBROUTINE rttov_transmit_9_solar_tl( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & chanprof,                        &
            & profiles,                        &
            & sun,                             &
            & aux,                             &
            & aux_tl,                          &
            & coef,                            &
            & raytracing,                      &
            & raytracing_tl,                   &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_tl,                    &
            & odsun_level,                     &
            & odsun_singlelayer,               &
            & od_frac,                         &
            & transmission_aux,                &
            & transmission_aux_tl,             &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & tausun_ref,                      &
            & tausun_ref_surf,                 &
            & tausun_level,                    &
            & tausun_surf)
! Description:
! Tangent linear of rttov_transmit
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPTL and RTTAUTL from previous versions of RTTOV
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_tl.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    06/02/2007  Modified to remove polarisation (R Saunders)
!  1.2    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
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
       & profile_Type,               &
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
  TYPE(rttov_chanprof)            , INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type  )            , INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(IN)    :: aux_tl
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_tl
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_tl
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing_tl
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_tl
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_level      (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_singlelayer(nlayers       , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref_surf  (size(chanprof)                )! sat to surface transmittance
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref       (nlayers + 1   , size(chanprof))! sat to layer transmittance
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_level     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_surf      (size(chanprof)                )
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_singlelayer_tl(nlayers       , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_tl    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_level_tl      (nlayers + 1   , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_frac_tl       (size(chanprof)                )
  REAL   (KIND=jprb) :: tausun_level_tl  (nlayers + 1   , size(chanprof))! sat to layer transmission_aux at each frequency
  REAL   (KIND=jprb) :: tausun_surf_tl   (size(chanprof)                )
  INTEGER(KIND=jpim) :: lev         , lay    , chan   , j, nlevels
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: ist         , levsurf
  INTEGER(KIND=jpim) :: nchannels                                        ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!----------------------------------------
!2. Compute layer to space optical depths
!----------------------------------------
! notes: apply gamma correction; check value is sensible and constrain
! if necessary.
! note that optical depth in the calculations is negative
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lay = 1, nlayers
        od_singlelayer_tl(lay, j) =  - (opdp_path_tl%sun_level(lay + 1, j) - opdp_path_tl%sun_level(lay, j))
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level_tl(lev, j) = opdp_path_tl%sun_level(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level_tl(lev, j) = coef%ff_gam(chan) * od_level_tl(lev, j)
      ENDDO
      DO lay = 1, nlayers
        od_singlelayer_tl(lay, j) = coef%ff_gam(chan) * od_singlelayer_tl(lay, j)
      ENDDO
    ENDIF
  ENDDO
! associated transmittances
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        tausun_level_tl(lev, j) = od_level_tl(lev, j) * tausun_ref(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tausun_ref(lev, j) < coef%tt_a0(chan)) THEN
            tausun_level_tl(lev, j) = 0._jprb
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (sun(j)) THEN
      od_surf_tl(j) = od_level_tl(levsurf, j) +                                                      &
        & aux_tl%s(prof)%pfraction_surf * (odsun_level(levsurf - 1, j) - odsun_level(levsurf, j)) +  &
        & aux%s(prof)%pfraction_surf * (od_level_tl(levsurf - 1, j) - od_level_tl(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac_tl(J) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
      ELSE
        od_frac_tl(J) = od_surf_tl(j) - od_level_tl(levsurf, j)
      ENDIF
      tausun_surf_tl(j) = od_surf_tl(j) * tausun_ref_surf(j)
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tausun_ref_surf(j) < coef%tt_a0(chan)) THEN
          tausun_surf_tl(j) = 0._jprb
        ENDIF
      ENDIF
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac_tl(j) = transmission_scatt_ir_stream_tl%opdpacsun(ist, j, levsurf) + aux%s(prof)%pfraction_surf     &
            &  * (transmission_scatt_ir_stream_tl%opdpacsun(ist, j, levsurf - 1) -                                       &
            & transmission_scatt_ir_stream_tl%opdpacsun(ist, j, levsurf)) + aux_tl%s(prof)%pfraction_surf * (            &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf - 1) -                                              &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            transmission_aux_tl%odsun_frac_ac(ist, J) =      &
              & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpacsun(ist, j, levsurf - 1)
          ELSE
            transmission_aux_tl%odsun_frac_ac(ist, J) =      &
              & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpacsun(ist, j, levsurf)
          ENDIF
          transmission_aux_tl%tau_surf_acsun(ist, j) =      &
            &  - od_surf_ac_tl(j) * transmission_aux%tau_ref_surf_acsun(ist, j)
          transmission_aux_tl%odsun_frac_t(ist, J)   =                                                                  &
            & ( - od_frac_tl(j) + transmission_aux_tl%odsun_frac_ac(ist, J)) * raytracing%pathsun(levsurf - 1, prof) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) -                         &
            & raytracing_tl%pathsat(levsurf - 1, prof) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *     &
            & raytracing%pathsun(levsurf - 1, prof) /                                                                   &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2 +                    &
            & raytracing_tl%pathsun(levsurf - 1, prof) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) /     &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) -                         &
            & raytracing_tl%pathsun(levsurf - 1, prof) * ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) *     &
            & raytracing%pathsun(levsurf - 1, prof) /                                                                   &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof)) ** 2
          IF (tausun_surf(j) >= 0) THEN
            transmission_aux_tl%tausun_surf_t(ist, j) = tausun_surf_tl(j) * transmission_aux%tau_surf_acsun(ist, J) +      &
              & transmission_aux_tl%tau_surf_acsun(ist, J) * tausun_surf(j)
          ELSE
            transmission_aux_tl%tausun_surf_t(ist, j) = tausun_surf_tl(j)
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------------------------
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_tl%tausun_level(:, ist, j)      = 0.0_JPRB
        transmission_aux_tl%tausun_surf(ist, j)          = 0.0_JPRB
        transmission_aux_tl%odsun_singlelayer(:, ist, j) = 0.0_JPRB
        transmission_aux_tl%odsun_sfrac(ist, j)          = 0.0_JPRB
      ENDDO
      DO ist = 0, ircld%nstream(prof)
        IF (addaerosl .OR. addclouds) THEN
          DO lev = 1, nlevels
            IF (tausun_level(lev, j) >= 0) THEN
              transmission_aux_tl%tausun_level(lev, ist, j) =                                              &
                & tausun_level_tl(lev, j) * exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev)) -  &
                & transmission_scatt_ir_stream_tl%opdpacsun(ist, j, lev) * tausun_level(lev, j) *          &
                & exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
            ELSE
              transmission_aux_tl%tausun_level(lev, ist, j) = tausun_level_tl(lev, j)
            ENDIF
          ENDDO
          transmission_aux_tl%tausun_surf(ist, j) = transmission_aux_tl%tausun_surf_t(ist, j)
          transmission_aux_tl%odsun_sfrac(ist, j) = transmission_aux_tl%odsun_frac_t(ist, J)
          DO lay = 1, nlayers
            lev = lay + 1
            transmission_aux_tl%odsun_singlelayer(lay, ist, j)   =                                                      &
              & (od_singlelayer_tl(lay, j) + transmission_scatt_ir_stream_tl%opdpaclsun(ist, j, lay)) *                 &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) -       &
              & raytracing_tl%pathsat(lay, prof) *                                                                      &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *                    &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2 +  &
              & raytracing_tl%pathsun(lay, prof) *                                                                      &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) /                    &
              & (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) - raytracing_tl%pathsun(lay, prof) *    &
              & (odsun_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *                    &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof)) ** 2
            transmission_scatt_ir_stream_tl%opdpext(ist, j, lay) =                                  &
              & od_singlelayer_tl(lay, j) + transmission_scatt_ir_stream_tl%opdpabs(ist, j, lay) +  &
              & transmission_scatt_ir_stream_tl%opdpsca(ist, j, lay)
            IF (transmission_scatt_ir_stream%opdpext(ist, j, lay) /= 0._jprb) THEN
              transmission_scatt_ir_stream_tl%ssa(ist, j, lay) =                                                           &
                & transmission_scatt_ir_stream_tl%opdpsca(ist, j, lay) / transmission_scatt_ir_stream%opdpext(ist, j, lay) &
                &  -                                                                                                       &
                & transmission_scatt_ir_stream_tl%opdpext(ist, j, lay) * transmission_scatt_ir_stream%opdpsca(ist, j, lay) &
                &  / transmission_scatt_ir_stream%opdpext(ist, j, lay) ** 2
            ENDIF
          ENDDO
        ELSE
          transmission_aux_tl%tausun_level(:, ist, j) = tausun_level_tl(:, j)
          transmission_aux_tl%tausun_surf(ist, j)     = tausun_surf_tl(j)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar_tl
