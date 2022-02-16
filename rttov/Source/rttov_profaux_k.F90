SUBROUTINE rttov_profaux_k( &
            & opts,       &
            & chanprof,   &
            & profiles,   &
            & profiles_k, &
            & coef,       &
            & aux_prof,   &
            & aux_prof_k)
!
! Description:
! K of rttov_profaux
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       22/06/2005  initial (P Brunel)
!                        based on version 1.0 (07/10/04) of AD code
!  1.1       12/02/2007  Removed polarisation index (R Saunders)
!  1.2       04/02/2008  opts%lgradp option for TL/AD of pressure levels (N Bormann)
!  1.3       15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4       02/12/2009  Fixed a number of bugs due to the wrong assumption that cloud
!                        related quantities are defined on levels (thay are layer
!                        average quantities). Marco Matricardi
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
       & rttov_options,  &
       & rttov_chanprof, &
       & rttov_coef,     &
       & profile_Type,   &
       & profile_aux
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, dcoeff
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(profile_aux   ), INTENT(IN)    :: aux_prof
  TYPE(profile_aux   ), INTENT(INOUT) :: aux_prof_k
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: dp  (profiles(1)%nlayers)
  REAL   (KIND=jprb) :: dp_k(profiles(1)%nlayers)
  INTEGER(KIND=jpim) :: lev         , lay
  INTEGER(KIND=jpim) :: i
  INTEGER(KIND=jpim) :: j
  REAL   (KIND=jprb) :: v                                                      ! temperature ratio
  REAL   (KIND=jprb) :: v_k                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_k
  REAL   (KIND=jprb) :: ztempc      , ztempc_k
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_k    , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  INTEGER(KIND=jpim) :: nchannels   ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  zmcfarq_k = 0._jprb
  ztempc_k  = 0._jprb
  bwyser_k  = 0._jprb
!------------------------------------------------------------
! K for ice water content into effective generalized diameter
!------------------------------------------------------------
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
    DO i = 1, nchannels
      j = chanprof(i)%prof
      DO lay = profiles(j)%nlayers, 1,  - 1
        lev = lay + 1
        IF (profiles(j)%cloud(6, lay) /= 0._jprb) THEN
          ! Use effective diameter from input profile if specified
          IF (profiles(j)%icede(lay) /= 0._jprb) THEN
            profiles_k(i)%icede(lay) = profiles_k(i)%icede(lay) + aux_prof_k%dg(lay, i)
          ELSE
            IF (profiles(j)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = profiles(j)%t(lev) - rtt
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux_prof%fac3_dg(lay, j) < zradipou_low .OR. aux_prof%fac3_dg(lay, j) > zradipou_upp) THEN
                aux_prof_k%dg(lay, i) = 0._JPRB
              ELSE
                aux_prof_k%fac3_dg(lay, i) = aux_prof_k%fac3_dg(lay, i) + aux_prof_k%dg(lay, i)
              ENDIF
              aux_prof_k%fac2_dg(lay, i) = aux_prof_k%fac2_dg(lay, i) + 2.0_JPRB * aux_prof_k%fac3_dg(lay, i)
              aux_prof_k%fac1_dg(lay, i) = aux_prof_k%fac1_dg(lay, i) + 0.388_JPRB * aux_prof_k%fac2_dg(lay, i)
              aux_prof_k%fac1_dg(lay, i) =      &
                & aux_prof_k%fac1_dg(lay, i) + aux_prof_k%fac2_dg(lay, i) * 0.00051_JPRB * 2 * aux_prof%fac1_dg(lay, j)
              ztempc_k                   = ztempc_k + aux_prof_k%fac1_dg(lay, i) *      &
                & (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              profiles_k(i)%t(lev)       = profiles_k(i)%t(lev) + ztempc_k
              ztempc_k                   = 0._jprb
            ELSE IF (profiles(j)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (profiles(j)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * &
                  ((273._JPRB - profiles(j)%t(lev)) ** 1.5_JPRB) * Log10(profiles(j)%cloud(6, lay) / 50._JPRB))
              ENDIF
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_prof_k%fac2_dg(lay, i) =      &
                & aux_prof_k%fac2_dg(lay, i) + aux_prof_k%dg(lay, i) * 2._JPRB * 4._JPRB * sqrt(3._JPRB) / 9._JPRB
              aux_prof_k%fac1_dg(lay, i) = aux_prof_k%fac1_dg(lay, i) + aux_prof_k%fac2_dg(lay, i) / nft
              bwyser_k                   = bwyser_k +      &
                & aux_prof_k%fac1_dg(lay, i) * (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB)
              IF (profiles(j)%t(lev) < 273._JPRB) THEN
                profiles_k(i)%t(lev)        = profiles_k(i)%t(lev) -                                     &
                  & bwyser_k * 0.001_JPRB * 1.5_JPRB * ((273._JPRB - profiles(j)%t(lev)) ** 0.5_JPRB) *  &
                  & Log10(profiles(j)%cloud(6, lay) / 50._JPRB)
                profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) +                                               &
                  & bwyser_k * 0.001_JPRB * ((273._JPRB - profiles(j)%t(lev)) ** 1.5_JPRB) / profiles(j)%cloud(6, lay) /  &
                  & log(10._JPRB)
              ENDIF
              bwyser_k = 0._JPRB
            ELSE IF (profiles(j)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc = profiles(j)%t(lev) - rtt
              profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) +                                 &
                & aux_prof_k%dg(lay, i) * 53.005_JPRB * (profiles(j)%cloud(6, lay) ** (0.06_JPRB - 1)) *  &
                & exp(0.013_JPRB * ztempc) * 0.06_JPRB
              ztempc_k                    = ztempc_k +                                                             &
                & aux_prof_k%dg(lay, i) * 53.005_JPRB * ((profiles(j)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB *  &
                & exp(0.013_JPRB * ztempc)
              profiles_k(i)%t(lev)        = profiles_k(i)%t(lev) + ztempc_k
              ztempc_k                    = 0._JPRB
            ELSE IF (profiles(j)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                     = 1.78449_JPRB
              bmcfarq                     = 0.281301_JPRB
              cmcfarq                     = 0.0177166_JPRB
              zmcfarq                     = profiles(j)%cloud(6, lay)
              aux_prof_k%fac1_dg(lay, i)  = aux_prof_k%fac1_dg(lay, i) + aux_prof_k%dg(lay, i) * 2.0_JPRB
              zmcfarq_k                   = zmcfarq_k + aux_prof_k%fac1_dg(lay, i) *                               &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq
              profiles_k(i)%cloud(6, lay) = profiles_k(i)%cloud(6, lay) + zmcfarq_k
              zmcfarq_k                   = 0._jprb
            ELSE
              aux_prof_k%dg(lay, i) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-------------------------------------------------
! K for debye terms for pewrmittivity calculations
!-------------------------------------------------
  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (opts%lgradp) dp_k = 0._jprb
    IF ((coef % id_sensor == sensor_id_mw .or. coef % id_sensor == sensor_id_po) .and. &
         opts%clw_data) THEN
      DO lev = profiles(j)%nlevels, coef%mwcldtop,  - 1
        v = 300.0_JPRB / profiles(j)%t(lev) - 1.0_JPRB
        aux_prof_k%debye_prof(3, lev, i) =      &
          & aux_prof_k%debye_prof(3, lev, i) + aux_prof_k%debye_prof(4, lev, i) * dcoeff(7)
        v_k = aux_prof_k%debye_prof(3, lev, i) * dcoeff(5)
        aux_prof_k%debye_prof(1, lev, i) =      &
          & aux_prof_k%debye_prof(1, lev, i) + aux_prof_k%debye_prof(2, lev, i) * dcoeff(4)
        v_k = v_k + aux_prof_k%debye_prof(1, lev, i) * ( - dcoeff(2) + 2 * dcoeff(3) * v)
        profiles_k(i)%t(lev)             = profiles_k(i)%t(lev) + v_k * ( - 300.0_JPRB / profiles(j)%t(lev) ** 2)
      ENDDO
!aux_prof_k(i) % debye_prof(:,:) = 0.
    ENDIF
!-----------------------------------
! K for cloud top and surface levels
!-----------------------------------
    !dp(1) = profiles(j)%p(2)
    !dp(2:profiles(j)%nlayers) = profiles(j)%p(3:profiles(j)%nlevels) - profiles(j)%p(2:profiles(j)%nlevels - 1)
    dp(1:profiles(j)%nlayers) = profiles(j)%p(2:profiles(j)%nlevels) - profiles(j)%p(1:profiles(j)%nlevels - 1)
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      profiles_k(i)%cfraction = profiles_k(i)%cfraction + aux_prof_k%s(i)%cfraction
      profiles_k(i)%ctp       = profiles_k(i)%ctp - aux_prof_k%s(i)%pfraction_ctp / dp(aux_prof%s(j)%nearestlev_ctp - 1)
      IF (opts%lgradp) THEN
        profiles_k(i)%p(aux_prof%s(j)%nearestlev_ctp) = profiles_k(i)%p(aux_prof%s(j)%nearestlev_ctp) +      &
          & aux_prof_k%s(i)%pfraction_ctp / dp(aux_prof%s(j)%nearestlev_ctp - 1)
        dp_k(aux_prof%s(j)%nearestlev_ctp - 1)        = dp_k(aux_prof%s(j)%nearestlev_ctp - 1) -                           &
          & (profiles(j)%p(aux_prof%s(j)%nearestlev_ctp) - profiles(j)%ctp) / dp(aux_prof%s(j)%nearestlev_ctp - 1) ** 2 *  &
          & aux_prof_k%s(i)%pfraction_ctp
      ENDIF
!Else
! for micro waves do not consider clouds in the RTTOV basis routines
    ENDIF
!nearest level above surface
    profiles_k(i)%s2m%p = profiles_k(i)%s2m%p - aux_prof_k%s(i)%pfraction_surf / dp(aux_prof%s(j)%nearestlev_surf - 1)
    IF (opts%lgradp) THEN
      profiles_k(i)%p(aux_prof%s(j)%nearestlev_surf) = profiles_k(i)%p(aux_prof%s(j)%nearestlev_surf) +      &
        & aux_prof_k%s(i)%pfraction_surf / dp(aux_prof%s(j)%nearestlev_surf - 1)
      dp_k(aux_prof%s(j)%nearestlev_surf - 1)        = dp_k(aux_prof%s(j)%nearestlev_surf - 1) -                             &
        & (profiles(j)%p(aux_prof%s(j)%nearestlev_surf) - profiles(j)%s2m%p) / dp(aux_prof%s(j)%nearestlev_surf - 1) ** 2 *  &
        & aux_prof_k%s(i)%pfraction_surf
      !profiles_k(i)%p(3:profiles(j)%nlevels)         =      &
      !  & profiles_k(i)%p(3:profiles(j)%nlevels) + dp_k(2:profiles(j)%nlayers)
      !profiles_k(i)%p(2:profiles(j)%nlevels - 1)     =      &
      !  & profiles_k(i)%p(2:profiles(j)%nlevels - 1) - dp_k(2:profiles(j)%nlayers)
      !profiles_k(i)%p(2)                             = profiles_k(i)%p(2) + dp_k(1)

      profiles_k(i)%p(2:profiles(j)%nlevels)         =      &
        & profiles_k(i)%p(2:profiles(j)%nlevels) + dp_k(1:profiles(j)%nlayers)
      profiles_k(i)%p(1:profiles(j)%nlevels - 1)     =      &
        & profiles_k(i)%p(1:profiles(j)%nlevels - 1) - dp_k(1:profiles(j)%nlayers)

    ENDIF
  ENDDO
! Channels loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_k
