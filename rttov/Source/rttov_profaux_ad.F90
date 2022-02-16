SUBROUTINE rttov_profaux_ad( &
            & opts,    &
            & prof,    &
            & prof_ad, &
            & coef,    &
            & aux,     &
            & aux_ad)
!
! Description:
! AD of rttov_profaux
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
!  1.0    07/10/2004 Added history
!  1.1    29/03/2005 Add end of header comment (J. Cameron)
!  1.2    04/02/2008 opts%lgradp option for TL/AD of pressure levels (N. Bormann)
!  1.3    15/08/2009 User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4    02/12/2009 Fixed a number of bugs due to the wrong assumption that cloud
!                    related quantities are defined on levels (thay are layer
!                    average quantities). Marco Matricardi
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
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, dcoeff
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(INOUT) :: prof_ad(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_ad
!INTF_END
! local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: dp_ad(prof(1)%nlayers)
  REAL   (KIND=jprb) :: v                                                       ! temperature ratio
  REAL   (KIND=jprb) :: v_ad                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_ad
  REAL   (KIND=jprb) :: ztempc      , ztempc_ad
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_ad   , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  DO iprof = 1, nprofiles
    zmcfarq_ad = 0._jprb
    ztempc_ad  = 0._jprb
    bwyser_ad  = 0._jprb
    dp_ad      = 0._jprb
!-----------------------------------------------------------------------------------------
! AD for ice water content into effective generalized diameter
!-----------------------------------------------------------------------------------------
    IF (opts%addclouds) THEN
! Calculate upper and lower limits for Ou-Liou effective size
!
      zradipou_upp = 326.3_JPRB + rtou_upp * (12.42_JPRB + rtou_upp * (0.197_JPRB + rtou_upp * 0.0012_JPRB))
      zradipou_low = 326.3_JPRB + rtou_low * (12.42_JPRB + rtou_low * (0.197_JPRB + rtou_low * 0.0012_JPRB))
!
! and convert these to the "generalized" effective size used here (using McFarquhar et al 2003 equation),
! not forgetting the factor of 2 to convert from McFarquhar's radius to a diameter
!
      zradipou_upp =  - 1.56_JPRB + zradipou_upp * (0.388_JPRB + zradipou_upp * 0.00051_JPRB)
      zradipou_upp = 2.0_JPRB * zradipou_upp
      zradipou_low =  - 1.56_JPRB + zradipou_low * (0.388_JPRB + zradipou_low * 0.00051_JPRB)
      zradipou_low = 2.0_JPRB * zradipou_low
      DO lay = prof(iprof)%nlayers, 1,  - 1
        lev = lay + 1
        IF (prof(iprof)%cloud(6, lay) /= 0._jprb) THEN
          ! Use effective diameter from input profile if specified
          IF (prof(iprof)%icede(lay) /= 0._jprb) THEN
            prof_ad(iprof)%icede(lay) = prof_ad(iprof)%icede(lay) + aux_ad%dg(lay, iprof)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = prof(iprof)%t(lev) - rtt
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux%fac3_dg(lay, iprof) < zradipou_low .OR. aux%fac3_dg(lay, iprof) > zradipou_upp) THEN
                aux_ad%dg(lay, iprof) = 0._JPRB
              ELSE
                aux_ad%fac3_dg(lay, iprof) = aux_ad%fac3_dg(lay, iprof) + aux_ad%dg(lay, iprof)
              ENDIF
              aux_ad%fac2_dg(lay, iprof) = aux_ad%fac2_dg(lay, iprof) + 2.0_JPRB * aux_ad%fac3_dg(lay, iprof)
              aux_ad%fac1_dg(lay, iprof) = aux_ad%fac1_dg(lay, iprof) + 0.388_JPRB * aux_ad%fac2_dg(lay, iprof)
              aux_ad%fac1_dg(lay, iprof) =      &
                & aux_ad%fac1_dg(lay, iprof) + aux_ad%fac2_dg(lay, iprof) * 0.00051_JPRB * 2 * aux%fac1_dg(lay, iprof)
              ztempc_ad                  = ztempc_ad + aux_ad%fac1_dg(lay, iprof) *      &
                & (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              prof_ad(iprof)%t(lev)      = prof_ad(iprof)%t(lev) + ztempc_ad
              ztempc_ad                  = 0._jprb
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * &
                  ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
              ENDIF
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_ad%fac2_dg(lay, iprof) =      &
                & aux_ad%fac2_dg(lay, iprof) + aux_ad%dg(lay, iprof) * 2._JPRB * 4._JPRB * sqrt(3._JPRB) / 9._JPRB
              aux_ad%fac1_dg(lay, iprof) = aux_ad%fac1_dg(lay, iprof) + aux_ad%fac2_dg(lay, iprof) / nft
              bwyser_ad                  = bwyser_ad +      &
                & aux_ad%fac1_dg(lay, iprof) * (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB)
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                prof_ad(iprof)%t(lev)        = prof_ad(iprof)%t(lev) -                                    &
                  & bwyser_ad * 0.001_JPRB * 1.5_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 0.5_JPRB) *  &
                  & Log10(prof(iprof)%cloud(6, lay) / 50._JPRB)
                prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) + &
                  bwyser_ad * 0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) / &
                  prof(iprof)%cloud(6, lay) / log(10._JPRB)
              ENDIF
              bwyser_ad = 0._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc = prof(iprof)%t(lev) - rtt
              prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) +                               &
                & aux_ad%dg(lay, iprof) * 53.005_JPRB * (prof(iprof)%cloud(6, lay) ** (0.06_JPRB - 1)) *  &
                & exp(0.013_JPRB * ztempc) * 0.06_JPRB
              ztempc_ad                    = ztempc_ad +                                                           &
                & aux_ad%dg(lay, iprof) * 53.005_JPRB * ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB *  &
                & exp(0.013_JPRB * ztempc)
              prof_ad(iprof)%t(lev)        = prof_ad(iprof)%t(lev) + ztempc_ad
              ztempc_ad                    = 0._JPRB
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                      = 1.78449_JPRB
              bmcfarq                      = 0.281301_JPRB
              cmcfarq                      = 0.0177166_JPRB
              zmcfarq                      = prof(iprof)%cloud(6, lay)
              aux_ad%fac1_dg(lay, iprof)   = aux_ad%fac1_dg(lay, iprof) + aux_ad%dg(lay, iprof) * 2.0_JPRB
              zmcfarq_ad                   = zmcfarq_ad + aux_ad%fac1_dg(lay, iprof) *                             &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq
              prof_ad(iprof)%cloud(6, lay) = prof_ad(iprof)%cloud(6, lay) + zmcfarq_ad
              zmcfarq_ad                   = 0._jprb
            ELSE
              aux_ad%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
!--------------------------------------------------
! AD for debye terms for pewrmittivity calculations
!--------------------------------------------------
    IF ((coef % id_sensor == sensor_id_mw .or. coef % id_sensor == sensor_id_po) .and. &
         opts%clw_data) THEN
      DO lev = prof(iprof)%nlevels, coef%mwcldtop,  - 1
        v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
        aux_ad%debye_prof(3, lev, iprof) =      &
          & aux_ad%debye_prof(3, lev, iprof) + aux_ad%debye_prof(4, lev, iprof) * dcoeff(7)
        v_ad = aux_ad%debye_prof(3, lev, iprof) * dcoeff(5)
        aux_ad%debye_prof(1, lev, iprof) =      &
          & aux_ad%debye_prof(1, lev, iprof) + aux_ad%debye_prof(2, lev, iprof) * dcoeff(4)
        v_ad = v_ad + aux_ad%debye_prof(1, lev, iprof) * ( - dcoeff(2) + 2 * dcoeff(3) * v)
        prof_ad(iprof)%t(lev)            = prof_ad(iprof)%t(lev) + v_ad * ( - 300.0_JPRB / prof(iprof)%t(lev) ** 2)
      ENDDO
!aux_ad % debye_prof(:,:) = 0.
    ENDIF
!-----------------------------------------
! AD for cloud top and surface levels
!-----------------------------------------
    dp(1:prof(iprof)%nlayers) = prof(iprof)%p(2:prof(iprof)%nlevels) - prof(iprof)%p(1:prof(iprof)%nlevels - 1)

    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      prof_ad(iprof)%cfraction = prof_ad(iprof)%cfraction + aux_ad%s(iprof)%cfraction
      prof_ad(iprof)%ctp       =      &
        & prof_ad(iprof)%ctp - aux_ad%s(iprof)%pfraction_ctp / dp(aux%s(iprof)%nearestlev_ctp - 1)
      IF (opts%lgradp) THEN
        prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_ctp) +      &
          & aux_ad%s(iprof)%pfraction_ctp / dp(aux%s(iprof)%nearestlev_ctp - 1)
        dp_ad(aux%s(iprof)%nearestlev_ctp - 1)        = dp_ad(aux%s(iprof)%nearestlev_ctp - 1) -                         &
          & (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp(aux%s(iprof)%nearestlev_ctp - 1) ** 2 *  &
          & aux_ad%s(iprof)%pfraction_ctp
      ENDIF
!Else
! for micro waves do not consider clouds in the RTTOV basis routines
    ENDIF
!nearest level above surface
    prof_ad(iprof)%s2m%p = prof_ad(iprof)%s2m%p - aux_ad%s(iprof)%pfraction_surf / dp(aux%s(iprof)%nearestlev_surf - 1)
    IF (opts%lgradp) THEN
      prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) = prof_ad(iprof)%p(aux%s(iprof)%nearestlev_surf) +      &
        & aux_ad%s(iprof)%pfraction_surf / dp(aux%s(iprof)%nearestlev_surf - 1)
      dp_ad(aux%s(iprof)%nearestlev_surf - 1)        = dp_ad(aux%s(iprof)%nearestlev_surf - 1) -                           &
        & (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp(aux%s(iprof)%nearestlev_surf - 1) ** 2 *  &
        & aux_ad%s(iprof)%pfraction_surf

!      prof_ad(iprof)%p(3:prof(iprof)%nlevels)        =      &
!        & prof_ad(iprof)%p(3:prof(iprof)%nlevels) + dp_ad(2:prof(iprof)%nlayers)
!      prof_ad(iprof)%p(2:prof(iprof)%nlevels - 1)    =      &
!        & prof_ad(iprof)%p(2:prof(iprof)%nlevels - 1) - dp_ad(2:prof(iprof)%nlayers)
!      prof_ad(iprof)%p(2)                            = prof_ad(iprof)%p(2) + dp_ad(1)

      prof_ad(iprof)%p(2:prof(iprof)%nlevels)        =      &
        & prof_ad(iprof)%p(2:prof(iprof)%nlevels) + dp_ad(1:prof(iprof)%nlayers)
      prof_ad(iprof)%p(1:prof(iprof)%nlevels - 1)        =      &
        & prof_ad(iprof)%p(1:prof(iprof)%nlevels - 1) - dp_ad(1:prof(iprof)%nlayers)

    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_ad
