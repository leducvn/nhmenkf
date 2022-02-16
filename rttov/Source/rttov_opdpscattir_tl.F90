!     Compute Optical parameters for aerosols and clouds
SUBROUTINE rttov_opdpscattir_tl( &
            & nlayers,                         &
            & chanprof,                        &
            & opts,                            &
            & aux,                             &
            & aux_tl,                          &
            & profiles,                        &
            & profiles_tl,                     &
            & sun,                             &
            & coef,                            &
            & coef_scatt_ir,                   &
            & raytracing,                      &
            & raytracing_tl,                   &
            & transmission_scatt_ir,           &
            & transmission_scatt_ir_tl,        &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & optp,                            &
            & ircld,                           &
            & ircld_tl)
!     Description:
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
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
! Current Code Owner: SAF NWP
!
!     Method:
!     See: User's manual and scientific report for RTIASI-5
!          (Available from EUMETSAT)
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           30/07/2004   Marco Matricardi. ECMWF.
!     1.1         05/02/2007   Removed polarisation R Saunders
!     1.2         15/07/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     1.3         03/11/2009   Transmittances / optical depths on levels (A Geer)
!     1.4         02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     1.5         02/12/2009   Introduced multiple cloud types in a single layer. Pathsun, Pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     1.6         05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     1.7         02/08/2010   Allow variable pressure levels (lgradp) (N Bormann)
!     1.8         14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                              should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_options,              &
       & rttov_coef,                 &
       & profile_Type,               &
       & raytracing_type,            &
       & transmission_scatt_ir_type, &
       & rttov_coef_scatt_ir,        &
       & rttov_optpar_ir,            &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       & pi,      &
       & deg2rad, &
       & jpazn,   &
       & t00,     &
       & ti,      &
       & ncldtyp, &
       & max_sol_zen
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof   (:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type              ), INTENT(IN)    :: profiles_tl(size(profiles))
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_tl
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_tl
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_tl
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type),                 INTENT(IN)    :: ircld
  TYPE(ircld_type),                 INTENT(INOUT) :: ircld_tl
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_tl
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: j, ich   , ich1   , i, ipf, ish, k, kk, ityp, iconf, iae, ikk
  INTEGER(KIND=jpim) :: lev         , lay   , lctyp
  INTEGER(KIND=jpim) :: iang        , iend  , ist
  REAL   (KIND=jprb) :: opd         , opdsun
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: scattangle
  REAL   (KIND=jprb) :: deltapice
  REAL   (KIND=jprb) :: deltapice_d
  REAL   (KIND=jprb) :: deltapaer
  REAL   (KIND=jprb) :: deltapaer_d
  REAL   (KIND=jprb) :: deltapwcl_d
  REAL   (KIND=jprb) :: deltaaer
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: deltawcl
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: phasint
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
  REAL   (KIND=jprb) :: phup
  REAL   (KIND=jprb) :: phdo
  REAL   (KIND=jprb) :: scattangle_d
  REAL   (KIND=jprb) :: musat_d
  REAL   (KIND=jprb) :: musun_d
  REAL   (KIND=jprb) :: frach_d
  REAL   (KIND=jprb) :: absch_d
  REAL   (KIND=jprb) :: scach_d
  REAL   (KIND=jprb) :: bparh_d
!       Local arrays:
!      REAL OPDPE      (KNCHPF,JPLAY)
!      REAL OPDPEAER   (KNCHPF,JPLAY)
  REAL   (KIND=jprb) :: opdpaerl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpaerlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: phasice_d  (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: phash_d    (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=JPRB) :: ZPI
  INTEGER(KIND=jpim) :: nprofiles                                                    ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                    ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles   = size(profiles)
  nchannels   = size(chanprof)
  OPDPCLDL    = 0._jprb
  OPDPAERL    = 0._jprb
  OPDPCLDLSUN = 0._jprb
  OPDPAERLSUN = 0._jprb
  bparh       = 0._jprb
  ZPI         = deg2rad / (2._jprb * PI)
  IANG        = 360_jpim / (JPAZN + 1)
  IEND        = IANG * JPAZN
!-----------------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-----------------------------------------------------------------------------------------
!-----Compute relative humidity-----------------------------------------------------------
  IF (opts%addaerosl) THEN
    DO j = 1, nprofiles
      DO lay = 1, nlayers
        lev = lay + 1
        ircld_tl%tave(lay, j)     = (profiles_tl(j)%t(lev - 1) + profiles_tl(j)%t(lev)) / 2._jprb
        ircld_tl%wmixave(lay, j)  = (profiles_tl(j)%q(lev - 1) + profiles_tl(j)%q(lev)) / 2._jprb
        If(opts%lgradp) ircld_tl%xpresave(lay, j) = (profiles_tl(j)%p(lev - 1) + profiles_tl(j)%p(lev)) / 2._jprb
!-----------Compute vater vapour partial pressure-----------------------------------------
        ircld_tl%ESW(lay, j)     = ircld_tl%TAVE(lay, j) * ircld%ESW(lay, j) * 17.502_jprb * (T00 - 32.19_jprb) /      &
          & (ircld%TAVE(lay, j) - 32.19_jprb) ** 2
        ircld_tl%ESI(lay, j)     = ircld_tl%TAVE(lay, j) * ircld%ESI(lay, j) * 22.587_jprb * (T00 + 0.7_jprb) /      &
          & (ircld%TAVE(lay, j) + 0.7_jprb) ** 2
        IF (ircld%TAVE(lay, j) > T00) THEN
          ircld_tl%PPV(lay, j) = ircld_tl%ESW(lay, j)! Water phase
        ELSE IF (ircld%TAVE(lay, j) > TI .AND. ircld%TAVE(lay, j) <= T00) THEN
          ircld_tl%PPV(lay, j) =                                                                             &
            & ircld_tl%ESI(lay, j) + ircld_tl%ESW(lay, j) * ((ircld%TAVE(lay, j) - TI) / (T00 - TI)) ** 2 -  &
            & ircld_tl%ESI(lay, j) * ((ircld%TAVE(lay, j) - TI) / (T00 - TI)) ** 2 +                         &
            & ircld_tl%TAVE(lay, j) * (ircld%ESW(lay, j) - ircld%ESI(lay, j)) * 2 *                          &
            & ((ircld%TAVE(lay, j) - TI) / (T00 - TI) ** 2)! Mixed phase
        ELSE IF (ircld%TAVE(lay, j) <= TI) THEN
          ircld_tl%PPV(lay, j) = ircld_tl%ESI(lay, j)! Ice phase
        ENDIF
        ircld_tl%PPV(lay, j)  = ircld_tl%PPV(lay, j) / 100._jprb
!-----------------------------------------------------------------------------------------
! layer average relative humidity
        aux_tl%RELHUM(lay, j) = ircld_tl%wmixave(lay, j) * 100._jprb * 1e-6_jprb * 0.622_jprb * ircld%xpresave(lay, j)     &
          &  / (ircld%PPV(lay, j) * (0.622_jprb + ircld%wmixave(lay, j) * 1e-6_jprb * 0.622_jprb)) -                       &
          & ircld_tl%wmixave(lay, j) * 100._jprb * (1e-6_jprb * 0.622_jprb) ** 2 * ircld%xpresave(lay, j) *                &
          & ircld%wmixave(lay, j) * ircld%PPV(lay, j) /                                                                    &
          & (ircld%PPV(lay, j) * (0.622_jprb + ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb)) ** 2 -                    &
          & ircld_tl%PPV(lay, j) * 100._jprb * ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb * ircld%xpresave(lay, j) *  &
          & (0.622_jprb + ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb) /                                               &
          & (ircld%PPV(lay, j) * (0.622_jprb + ircld%wmixave(lay, j) * 1e-6_jprb * 0.622_jprb)) ** 2
        IF(opts%lgradp) aux_tl%RELHUM(lay, j) = aux_tl%RELHUM(lay, j) + &
          &   100._jprb * ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb * ircld_tl%xpresave(lay, j) /      &
          &   (ircld%ppv(lay, j) * (0.622_jprb + ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb))
        IF (aux%RELHUMREF(lay, j) > 99._jprb) THEN
          aux_tl%RELHUM(lay, j) = 0._jprb
        ENDIF
      ENDDO
! layers
    ENDDO
! profiles
  ENDIF
! opts%addaerosl
!-----------------------------------------------------------------------------------------
  IF (opts%addaerosl .OR. opts%addclouds) THEN
    transmission_scatt_ir_tl%OPDPAAER(:,:)    = 0._jprb
    transmission_scatt_ir_tl%OPDPSAER(:,:)    = 0._jprb
    transmission_scatt_ir_tl%OPDPA(:,:)       = 0._jprb
!     transmission_scatt_ir_tl%OPDPACLS(:,:,:)  = 0._jprb
    transmission_scatt_ir_tl%OPDPS(:,:)       = 0._jprb
!     transmission_scatt_ir_tl%OPDPSCLS(:,:,:)  = 0._jprb
    transmission_scatt_ir_tl%GPARAER(:,:)     = 0._jprb
    transmission_scatt_ir_tl%GPAR(:,:)        = 0._jprb
    transmission_scatt_ir_tl%GPARTOT(:,:)     = 0._jprb
!     transmission_scatt_ir_tl%GPARCLS(:,:,:)   = 0._jprb
    transmission_scatt_ir_tl%OPDPAERLA(:,:)   = 0._jprb
    transmission_scatt_ir_tl%GPARAERA(:,:)    = 0._jprb
    transmission_scatt_ir_tl%azphaerup(:,:)   = 0._jprb
    transmission_scatt_ir_tl%azphaerupa(:,:)  = 0._jprb
    transmission_scatt_ir_tl%azphaerdo(:,:)   = 0._jprb
    transmission_scatt_ir_tl%azphaerdoa(:,:)  = 0._jprb
    transmission_scatt_ir_tl%azphupcls(:,:,:) = 0._jprb
    transmission_scatt_ir_tl%azphdocls(:,:,:) = 0._jprb
    transmission_scatt_ir_tl%azphuptot(:,:)   = 0._jprb
    transmission_scatt_ir_tl%azphdotot(:,:)   = 0._jprb
  ENDIF
  IF (opts%addaerosl) THEN
    DO j = 1, nchannels
      ich = chanprof(J)%chan
      ipf = chanprof(J)%prof
      DO lay = 1, nlayers
        lev = lay + 1
        DO i = 1, aux%iaernum(lay, ipf)
          iae = aux%iaertyp(I, lay, ipf)
          IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
!---------------Interpolate scattering parameters to actual value of relative-------------
!               humidity.
            DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
              IF (aux%relhum(lay, ipf) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                & aux%relhum(lay, ipf) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                delth   = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                frach   = aux_tl%relhum(lay, ipf)
                afac    = (optp%optpaer(iae)%abs(ich, k + 1) - optp%optpaer(iae)%abs(ich, k)) / delth
                sfac    = (optp%optpaer(iae)%sca(ich, k + 1) - optp%optpaer(iae)%sca(ich, k)) / delth
!!                    efac=(EXTC(ICH,IAE,K+1)-EXTC(ICH,IAE,K))/delth
                gfac    = (optp%optpaer(iae)%bpr(ich, k + 1) - optp%optpaer(iae)%bpr(ich, k)) / delth
                frach_d = (aux%relhum(lay, ipf) - optp%optpaer(iae)%fmv_aer_rh_val(k))
                absch_d = optp%optpaer(iae)%abs(ich, k) + afac * frach_d
                scach_d = optp%optpaer(iae)%sca(ich, k) + sfac * frach_d
!!                    EXTCH_D=EXTC(ICH,IAE,K)+EFAC*FRACH_D
                bparh_d = optp%optpaer(iae)%bpr(ich, k) + gfac * frach_d
                absch   = afac * frach
                scach   = sfac * frach
                bparh   = gfac * frach
                IF (sun(j)) THEN
                  ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
                  pfac(1:coef_scatt_ir%fmv_aer_ph)    = (                               &
                    & optp%optpaer(iae)%pha(ich1, K + 1, 1:coef_scatt_ir%fmv_aer_ph) -  &
                    & optp%optpaer(iae)%pha(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph)) / delth
                  phash_d(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph) +      &
                    & pfac(1:coef_scatt_ir%fmv_aer_ph) * frach_d
                  phash(1:coef_scatt_ir%fmv_aer_ph)   = pfac(1:coef_scatt_ir%fmv_aer_ph) * frach
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ELSE
            absch_d = optp%optpaer(iae)%abs(ich, 1)
            scach_d = optp%optpaer(iae)%sca(ich, 1)
            bparh_d = optp%optpaer(iae)%bpr(ich, 1)
            absch   = 0._jprb
            scach   = 0._jprb
            bparh   = 0._jprb
            IF (sun(j)) THEN
              ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
              phash_d(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(ich1, 1, 1:coef_scatt_ir%fmv_aer_ph)
              phash(1:coef_scatt_ir%fmv_aer_ph)   = 0._jprb
            ENDIF
          ENDIF
!-------------Compute optical parameters considering the contribution of------------------
!             all the aerosol components present in the layer
          IF (i == 1) THEN
! NB lev=lay+1, so lev labels lower level of the layer
            transmission_scatt_ir_tl%opdpaaer(j, lay) =                                                     &
              & profiles_tl(ipf)%aerosols(iae, lay) * absch_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & absch * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * absch_d * profiles(ipf)%aerosols(iae, lay)
            transmission_scatt_ir_tl%opdpsaer(j, lay) =                                                     &
              & profiles_tl(ipf)%aerosols(iae, lay) * scach_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(iae, lay)
            transmission_scatt_ir_tl%gparaera(j, lay) =                                                               &
              & profiles_tl(ipf)%aerosols(iae, lay) * scach_d * bparh_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) * bparh_d +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(iae, lay) * bparh_d +  &
              & bparh * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) * scach_d
          ELSE
            transmission_scatt_ir_tl%opdpaaer(j, lay) = transmission_scatt_ir_tl%opdpaaer(j, lay) +         &
              & profiles_tl(ipf)%aerosols(iae, lay) * absch_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & absch * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * absch_d * profiles(ipf)%aerosols(iae, lay)
            transmission_scatt_ir_tl%opdpsaer(j, lay) = transmission_scatt_ir_tl%opdpsaer(j, lay) +         &
              & profiles_tl(ipf)%aerosols(iae, lay) * scach_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(iae, lay)
            transmission_scatt_ir_tl%gparaera(j, lay) = transmission_scatt_ir_tl%gparaera(j, lay) +                   &
              & profiles_tl(ipf)%aerosols(iae, lay) * scach_d * bparh_d * raytracing%ltick(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) * bparh_d +       &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(iae, lay) * bparh_d +  &
              & bparh * profiles(ipf)%aerosols(iae, lay) * raytracing%ltick(nlayers + 1 - lay, ipf) * scach_d
          ENDIF
!-------------If solar radiation is present,compute the azimuthally averaged--------------
!             value of the phase function for the given value of the viewing--------------
!             angle and solar zenith angle.
          IF (sun(j)) THEN
            ich1    = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
!-------------Average phase function for the upward scattered solar beam------------------
            phasint = 0._jprb
!                 cosan_aer(1:coef_scatt_ir%fmv_aer_ph)  =                               &
!                   & cos(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph)*deg2rad)
            musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
            musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
!                  musat           = -raytracing_tl%pathsat(lay,ipf)                    / &
!                                   raytracing%pathsat(lay,ipf)*raytracing%pathsat(lay,ipf)
!                  musun           =  raytracing_tl%pathsun(lay,ipf)                    / &
!                                   raytracing%pathsun(lay,ipf)*raytracing%pathsun(lay,ipf)
            musat   =  - raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2
            musun   = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2
            loop1 : DO K = 0, iend, iang
              scattangle_d = musat_d * musun_d + sqrt((1 - musat_d ** 2) * (1 - musun_d ** 2)) * cos(k * deg2rad)
              IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                scattangle = musat * musun_d + musun * musat_d
              ELSE
                scattangle = musat * musun_d + musun * musat_d -                                                        &
                  & musat * musat_d * sqrt(1._jprb - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) -  &
                  & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
              ENDIF
!                     do kk=1,coef_scatt_ir%fmv_aer_ph-1
!                       if(scattangle_d>=coef_scatt_ir%fmv_aer_ph_val_cos(kk+1).and.   &
!                                                  & scattangle_d<=coef_scatt_ir%fmv_aer_ph_val_cos(kk))then
              ikk         = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad), jpim))
              KK = INT (coef_scatt_ir%ifmv_aer_ph_val(ikk) - 1._jprb, jpim)
              deltapaer   = (phash(kk + 1) - phash(kk))
              deltapaer_d = (phash_d(kk + 1) - phash_d(kk))
              deltaaer    = (coef_scatt_ir%fmv_aer_ph_val_cos(kk) - coef_scatt_ir%fmv_aer_ph_val_cos(kk + 1))
              phup        = phash(kk) + deltapaer * (coef_scatt_ir%fmv_aer_ph_val_cos(kk) - scattangle_d) / deltaaer     &
                &  - scattangle * deltapaer_d / deltaaer
!                         exit
!                       endif
!                     enddo
              phasint     = phasint + phup * iang * ZPI
            ENDDO loop1
            transmission_scatt_ir_tl%azphaerup(j, lay) = transmission_scatt_ir_tl%azphaerup(j, lay) +              &
              & profiles_tl(ipf)%aerosols(IAE, lay) * transmission_scatt_ir%phasintupref(j, lay, i) * scach_d *    &
              & raytracing%LTICK(nlayers + 1 - lay, ipf) +                                                         &
              & phasint * profiles(ipf)%aerosols(iae, lay) * scach_d * raytracing%LTICK(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(iae, lay) * transmission_scatt_ir%phasintupref(j, lay, i) *         &
              & raytracing%LTICK(nlayers + 1 - lay, ipf) +                                                         &
              & raytracing_tl%LTICK(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(iae, lay) *         &
              & transmission_scatt_ir%phasintupref(j, lay, i)
!-------------Average phase function for the downward scattered solar beam------
            phasint = 0._jprb
!                 cosan_aer(1:coef_scatt_ir%fmv_aer_ph)  =                               &
!                    & cos(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph)*deg2rad)
            musat_d =  - 1._jprb / raytracing%pathsat(lay, ipf)
            musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
            musat = raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2!*
!                                    raytracing%pathsat(lay,ipf)
            musun = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2!*
!                                    raytracing%pathsun(lay,ipf)
            loop2 : DO k = 0, iend, iang
              scattangle_d = musat_d * musun_d + sqrt((1 - musat_d ** 2) * (1 - musun_d ** 2)) * cos(k * deg2rad)
              IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                scattangle = musat * musun_d + musun * musat_d
              ELSE
                scattangle = musat * musun_d + musun * musat_d -                                                  &
                  & musat * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) -  &
                  & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
              ENDIF
!                     do kk=1,coef_scatt_ir%fmv_aer_ph-1
!                       if(scattangle_d>=coef_scatt_ir%fmv_aer_ph_val_cos(kk+1).and.                            &
!                                                  & scattangle_d<=coef_scatt_ir%fmv_aer_ph_val_cos(kk))then
              ikk         = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad), jpim))
              KK = INT (coef_scatt_ir%ifmv_aer_ph_val(ikk) - 1._jprb, jpim)
              deltapaer_d = (phash_d(kk + 1) - phash_d(kk))
              deltapaer   = (phash(kk + 1) - phash(kk))
              deltaaer    = (coef_scatt_ir%fmv_aer_ph_val_cos(kk) - coef_scatt_ir%fmv_aer_ph_val_cos(kk + 1))
              phdo        = phash(kk) + deltapaer * (coef_scatt_ir%fmv_aer_ph_val_cos(kk) - scattangle_d) / deltaaer     &
                &  - scattangle * deltapaer_d / deltaaer
!                         exit
!                       endif
!                     enddo
              phasint     = phasint + phdo * iang * ZPI
            ENDDO loop2
            transmission_scatt_ir_tl%azphaerdo(j, lay) = transmission_scatt_ir_tl%azphaerdo(j, lay) +              &
              & profiles_tl(ipf)%aerosols(IAE, lay) * transmission_scatt_ir%phasintdoref(j, lay, i) * scach_d *    &
              & raytracing%LTICK(nlayers + 1 - lay, ipf) +                                                         &
              & phasint * profiles(ipf)%aerosols(IAE, lay) * scach_d * raytracing%LTICK(nlayers + 1 - lay, ipf) +  &
              & scach * profiles(ipf)%aerosols(IAE, lay) * transmission_scatt_ir%phasintdoref(j, lay, i) *         &
              & raytracing%LTICK(nlayers + 1 - lay, ipf) +                                                         &
              & raytracing_tl%ltick(nlayers + 1 - lay, ipf) * scach_d * profiles(ipf)%aerosols(IAE, lay) *         &
              & transmission_scatt_ir%phasintdoref(j, lay, i)
          ENDIF
        ENDDO
      ENDDO
!---------Compute final values for optical parameters-------------------------------------
      DO lay = 1, nlayers
        lev = lay + 1
        IF (transmission_scatt_ir%opdpsaer(j, lay) /= 0._jprb) THEN
          transmission_scatt_ir_tl%gparaer(j, lay) =                                                &
            & transmission_scatt_ir_tl%gparaera(j, lay) / transmission_scatt_ir%opdpsaer(J, lay) -  &
            & transmission_scatt_ir_tl%opdpsaer(j, lay) * transmission_scatt_ir%gparaera(J, lay) /  &
            & transmission_scatt_ir%opdpsaer(j, lay) ** 2
          IF (sun(j)) THEN
            transmission_scatt_ir_tl%azphaerupa(j, lay) =                                              &
              & transmission_scatt_ir_tl%azphaerup(j, lay) / transmission_scatt_ir%opdpsaer(j, lay) -  &
              & transmission_scatt_ir_tl%opdpsaer(j, lay) * transmission_scatt_ir%azphaerup(j, lay) /  &
              & transmission_scatt_ir%opdpsaer(j, lay) ** 2
            transmission_scatt_ir_tl%azphaerdoa(j, lay) =                                              &
              & transmission_scatt_ir_tl%azphaerdo(j, lay) / transmission_scatt_ir%opdpsaer(j, lay) -  &
              & transmission_scatt_ir_tl%opdpsaer(j, lay) * transmission_scatt_ir%azphaerdo(j, lay) /  &
              & transmission_scatt_ir%opdpsaer(J, lay) ** 2
          ENDIF
        ENDIF
        transmission_scatt_ir_tl%opdpaerla(j, lay) = transmission_scatt_ir_tl%opdpaaer(j, lay) +      &
          & transmission_scatt_ir_tl%opdpsaer(j, lay) * transmission_scatt_ir%gparaer(j, lay) +       &
          & transmission_scatt_ir_tl%gparaer(j, lay) * transmission_scatt_ir%opdpsaer(j, lay)
        opdpaerl(j, lay)                           =                                                        &
          & transmission_scatt_ir_tl%opdpaerla(j, lay) * raytracing%pathsat(lay, ipf) * coef%ff_gam(ich) +  &
          & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%opdpaerla(j, lay) * coef%ff_gam(ich)
        IF (sun(j)) THEN
          opdpaerlsun(j, lay) =                                                                                          &
            & transmission_scatt_ir_tl%opdpaerla(j, lay) * (raytracing%pathsun(lay, ipf) + raytracing%pathsat(lay, ipf)) &
            &  * coef%ff_gam(ich) +                                                                                      &
            & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%opdpaerla(j, lay) * coef%ff_gam(ich) +             &
            & raytracing_TL%pathsun(lay, ipf) * transmission_scatt_ir%opdpaerla(j, lay) * coef%ff_gam(ich)
        ENDIF
      ENDDO
! layers
    ENDDO
! channels
  ENDIF
! opts%addaerosl
!-------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-------------------------------------------------------------------------------
  IF (opts%addclouds) THEN
    DO j = 1, nchannels
      ich = chanprof(j)%chan
      ipf = chanprof(j)%prof
      ish = profiles(ipf)%ish
!           IOFF=NUMAE
      DO lay = 1, nlayers
        lev = lay + 1
        transmission_scatt_ir_tl%azphup(j, lay) = 0._jprb
        transmission_scatt_ir_tl%azphdo(j, lay) = 0._jprb
        transmission_scatt_ir_tl%OPDPA(j, lay)  = 0._jprb
        transmission_scatt_ir_tl%OPDPS(j, lay)  = 0._jprb
        transmission_scatt_ir_tl%GPAR(j, lay)   = 0._jprb
!               transmission_scatt_ir_tl%OPDPAAER(j,lay)=0._jprb
!               transmission_scatt_ir_tl%OPDPSAER(j,lay)=0._jprb
!               transmission_scatt_ir_tl%GPARAER(j,lay) =0._jprb
        DO lctyp = 1, ncldtyp
          IF (ircld%cldtyp(lctyp, lay, ipf) /= 0) THEN
            ityp  = ircld%cldtyp(lctyp, lay, ipf)
            iconf = ityp
!                   ITMP = ITEMP  (lay,IPF)
!                   ICONF= CLDTYP (lay,IPF)+ITEMP(lay,IPF)-1
!-----------------Compute cloud  optical parameters ----------------------------------------
            IF (ityp <= 5_jpim) THEN
!-----------------------------------------------------------------------------------------
!                   For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
              transmission_scatt_ir_tl%opdpacls(J, lay, ityp) =                                                              &
                & profiles_tl(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%abs(ich, 1) *         &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) +                                                                 &
                & raytracing_tl%LTICK(nlayers + 1 - lay, IPF) * profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) &
                &  * optp%optpwcl(ityp)%abs(ich, 1)
              transmission_scatt_ir_tl%OPDPA(J, lay)          =      &
                & transmission_scatt_ir_tl%OPDPA(J, lay) + transmission_scatt_ir_tl%OPDPACLS(J, lay, ityp)
              transmission_scatt_ir_tl%opdpscls(j, lay, ityp) =                                                              &
                & profiles_tl(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%sca(ich, 1) *         &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) +                                                                 &
                & raytracing_tl%LTICK(nlayers + 1 - lay, IPF) * profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) &
                &  * optp%optpwcl(ityp)%sca(ich, 1)
              transmission_scatt_ir_tl%OPDPS(J, lay)          =      &
                & transmission_scatt_ir_tl%OPDPS(J, lay) + transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
!                OPDPE(J,jl)=PROFCLD(ITYP,jl,IPF)*CONFAC(ICONF)                            * &
!                              EXTC(ICH,IOFF+ITYP,ITMP)*LTICK(90-jl,IPF)                   + &
!                            LTICK(90-lay,IPF)*PROFCLD_D(ITYP,lay,IPF)                     * &
!                              CONFAC(ICONF)*EXTC(ICH,IOFF+ITYP,ITMP)
              transmission_scatt_ir_tl%GPARCLS(J, lay, ityp)  = 0._jprb
              transmission_scatt_ir_tl%GPARTOT(J, lay)        = transmission_scatt_ir_tl%GPARTOT(J, lay) +         &
                & transmission_scatt_ir_tl%GPARCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                & transmission_scatt_ir%GPARCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
            ELSE
!-----------------------------------------------------------------------------------------
!                   For ice clouds optical parameters are computed using regression
!                   coefficients
!-----------------------------------------------------------------------------------------
              transmission_scatt_ir_tl%OPDPACLS(J, lay, ityp) = profiles_tl(IPF)%cloud(Ityp, lay) * (                       &
                & optp%optpicl(ish)%abs(ich, 1) + optp%optpicl(ish)%abs(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%abs(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%abs(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF) +                                                             &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%abs(ich, 2) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) -                                                                &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%abs(ich, 3) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) / aux%dg(lay, ipf) ** 2 -                                        &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%abs(ich, 4) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) * 2._jprb / aux%dg(lay, ipf) ** 3 +                              &
                & raytracing_tl%LTICK(nlayers + 1 - lay, IPF) * profiles(IPF)%cloud(Ityp, lay) * (                          &
                & optp%optpicl(ish)%abs(ich, 1) + optp%optpicl(ish)%abs(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%abs(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%abs(ich, 4) / aux%dg(lay, ipf) ** 2)
              transmission_scatt_ir_tl%OPDPA(J, lay)          =      &
                & transmission_scatt_ir_tl%OPDPA(J, lay) + transmission_scatt_ir_tl%OPDPACLS(J, lay, ityp)
              transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp) = profiles_tl(IPF)%cloud(Ityp, lay) * (                       &
                & optp%optpicl(ish)%sca(ich, 1) + optp%optpicl(ish)%sca(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%sca(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%sca(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF) +                                                             &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%sca(ich, 2) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) -                                                                &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%sca(ich, 3) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) / aux%dg(lay, ipf) ** 2 -                                        &
                & aux_tl%dg(lay, ipf) * optp%optpicl(ish)%sca(ich, 4) * profiles(IPF)%cloud(Ityp, lay) *                    &
                & raytracing%LTICK(nlayers + 1 - lay, IPF) * 2._jprb / aux%dg(lay, ipf) ** 3 +                              &
                & raytracing_tl%LTICK(nlayers + 1 - lay, IPF) * profiles(IPF)%cloud(Ityp, lay) * (                          &
                & optp%optpicl(ish)%sca(ich, 1) + optp%optpicl(ish)%sca(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%sca(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%sca(ich, 4) / aux%dg(lay, ipf) ** 2)
              transmission_scatt_ir_tl%OPDPS(J, lay)          =      &
                & transmission_scatt_ir_tl%OPDPS(J, lay) + transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
              transmission_scatt_ir_tl%GPARCLS(J, lay, ityp)  = aux_tl%dg(lay, ipf) * optp%optpicl(ish)%bpr(ich, 2) +      &
                & aux_tl%dg(lay, ipf) * 2._jprb * aux%dg(lay, ipf) * optp%optpicl(ish)%bpr(ich, 3) +                       &
                & aux_tl%dg(lay, ipf) * 3._jprb * aux%dg(lay, ipf) ** 2 * optp%optpicl(ish)%bpr(ich, 4)
              transmission_scatt_ir_tl%GPARTOT(J, lay)        = transmission_scatt_ir_tl%GPARTOT(J, lay) +         &
                & transmission_scatt_ir_tl%GPARCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                & transmission_scatt_ir%GPARCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
            ENDIF
!ice/water
!-------------------If solar radiation is present,compute the azimuthally averaged----
!                   value of the phase function for the given value of the viewing----
!                   angle and solar zenith angle.
            IF (sun(j)) THEN
!-------------------If ice clouds are present the phase function for for the current value
!                   of the effective generalized diameter is obtained by linear
!                   interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
              IF (ITYP == 6_jpim) THEN
                ich1 = ich - coef_scatt_ir%fmv_icl_pha_ioff + 1
                DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                  IF (aux%dg(lay, ipf) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                    & aux%dg(lay, ipf) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                    deltap       = (optp%optpicl(ish)%PHA(ICH1, k + 1, :) - optp%optpicl(ish)%PHA(ICH1, k, :))
                    DELTAICE     = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                    phasice_d(:) = optp%optpicl(ish)%PHA(ICH1, k, :) +      &
                      & deltap(:) * (aux%dg(lay, ipf) - coef_scatt_ir%fmv_icl_dg(k, ish)) / DELTAICE
                    DO kk = 1, coef_scatt_ir%fmv_icl_ph
                      phasice(kk) = deltap(kk) * aux_tl%dg(lay, ipf) / DELTAICE
                    ENDDO
                    EXIT
                  ENDIF
                ENDDO
!-----------------------Average phase function for the upward scattered solar beam----------
                phasint = 0._jprb
!                       cosan_icl(1:coef_scatt_ir%fmv_icl_ph)  =                               &
!                       & cos(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph)*deg2rad)
                musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                musat   =  - raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2
                musun   = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2
                loop3 : DO k = 0, iend, iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt((1._jprb - musat_d ** 2) * (1._jprb - musun_d ** 2)) * cos(k * deg2rad)
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    scattangle = musat * musun_d + musun * musat_d
                  ELSE
                    scattangle = musat * musun_d + musun * musat_d -                                                     &
                      & musat * musat_d * sqrt(1._jprb - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) &
                      &  -                                                                                               &
                      & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
                  ENDIF
!                   do kk=1,coef_scatt_ir%fmv_icl_ph-1
!                     if(scattangle_d>=coef_scatt_ir%fmv_icl_ph_val_cos(kk+1) &
!                       &.and.scattangle_d<=coef_scatt_ir%fmv_icl_ph_val_cos(kk))then
                  ikk = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_icl_ph_val(ikk) - 1._jprb, jpim)
                  deltapice = (PHAsice(KK + 1) - PHAsice(KK))
                  deltapice_d = (PHAsice_d(KK + 1) - PHAsice_d(KK))
                  deltaice = (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - coef_scatt_ir%fmv_icl_ph_val_cos(KK + 1))
                  phasint = PHAsice(KK) + DELTAPice * (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - SCATTANGLE_d) / DELTAice     &
                    &  - SCATTANGLE * DELTAPice_d / DELTAice
!                         exit
!                     endif
!                   enddo
                  transmission_scatt_ir_tl%azphupcls(j, lay, ityp) =      &
                    & transmission_scatt_ir_tl%AZPHUPCLS(J, lay, ityp) + phasint * iang * ZPI
                ENDDO loop3
                transmission_scatt_ir_tl%AZPHUPTOT(J, lay) = transmission_scatt_ir_tl%AZPHUPTOT(J, lay) +              &
                  & transmission_scatt_ir_tl%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                  & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
!---------------------Average phase function for the downward scattered solar beam--------
                phasint = 0._jprb
!                      cosan_icl(1:coef_scatt_ir%fmv_icl_ph)  =                               &
!                      & cos(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph)*deg2rad)
                musat_d =  - 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                musat = raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2
                musun = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2
                loop4 : DO k = 0, iend, iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt((1._jprb - musat_d ** 2) * (1._jprb - musun_d ** 2)) * cos(k * deg2rad)
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    scattangle = musat * musun_d + musun * musat_d
                  ELSE
                    scattangle = musat * musun_d + musun * musat_d -                                                     &
                      & musat * musat_d * sqrt(1._jprb - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) &
                      &  -                                                                                               &
                      & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
                  ENDIF
!                   do kk=1,coef_scatt_ir%fmv_icl_ph-1
!                     if(scattangle_d>=coef_scatt_ir%fmv_icl_ph_val_cos(kk+1) &
!                       &.and.scattangle_d<=coef_scatt_ir%fmv_icl_ph_val_cos(kk))then
                  ikk = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_icl_ph_val(ikk) - 1._jprb, jpim)
                  deltapice_d = (PHAsice_d(KK + 1) - PHAsice_d(KK))
                  deltapice = (PHAsice(KK + 1) - PHAsice(KK))
                  deltaice = (coef_scatt_ir%fmv_icl_ph_val_cos(kk) - coef_scatt_ir%fmv_icl_ph_val_cos(kk + 1))
                  phasint = PHAsice(KK) + DELTAPice * (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - SCATTANGLE_d) / DELTAice     &
                    &  - SCATTANGLE * DELTAPice_d / DELTAice
!                        exit
!                     endif
!                   enddo
                  transmission_scatt_ir_tl%azphdocls(j, lay, ityp) =      &
                    & transmission_scatt_ir_tl%azphdocls(J, lay, ityp) + phasint * iang * ZPI
                ENDDO loop4
                transmission_scatt_ir_tl%AZPHDOTOT(J, lay) = transmission_scatt_ir_tl%AZPHDOTOT(J, lay) +              &
                  & transmission_scatt_ir_tl%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                  & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
              ELSE
!-----------------------------------------------------------------------------------------
!                       Water clouds
!-----------------------------------------------------------------------------------------
                ich1    = ich - coef_scatt_ir%fmv_wcl_pha_ioff + 1
!-----------------Average phase function for the upward scattered solar beam--------------
                phasint = 0._jprb
!                       cosan_wcl(1:coef_scatt_ir%fmv_wcl_ph)  =                               &
!                        & cos(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph)*deg2rad)
                musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                musat   =  - raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2
                musun   = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2
                loop5 : DO k = 0, iend, iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt((1._jprb - musat_d ** 2) * (1._jprb - musun_d ** 2)) * cos(k * deg2rad)
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    scattangle = musat * musun_d + musun * musat_d
                  ELSE
                    scattangle = musat * musun_d + musun * musat_d -                                                     &
                      & musat * musat_d * sqrt(1._jprb - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) &
                      &  -                                                                                               &
                      & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
                  ENDIF
!                   do kk=1,coef_scatt_ir%fmv_wcl_ph-1
!                     if(scattangle_d>=coef_scatt_ir%fmv_wcl_ph_val_cos(kk+1) &
!                       &.and.scattangle_d<=coef_scatt_ir%fmv_wcl_ph_val_cos(kk))then
                  ikk = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_wcl_ph_val(ikk) - 1._jprb, jpim)
                  deltapwcl_d = (optp%optpwcl(ityp)%PHA(ICH1, 1, KK + 1) - optp%optpwcl(ityp)%PHA(ICH1, 1, KK))
                  deltawcl = (coef_scatt_ir%fmv_wcl_ph_val_cos(KK) - coef_scatt_ir%fmv_wcl_ph_val_cos(KK + 1))
                  phasint =  - scattangle * deltapwcl_d / deltawcl
!                       exit
!                     endif
!                   enddo
                  transmission_scatt_ir_tl%azphupcls(j, lay, ityp) =      &
                    & transmission_scatt_ir_tl%AZPHUPCLS(J, lay, ityp) + phasint * iang * ZPI
                ENDDO loop5
                transmission_scatt_ir_tl%AZPHUPTOT(J, lay) = transmission_scatt_ir_tl%AZPHUPTOT(J, lay) +              &
                  & transmission_scatt_ir_tl%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                  & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
!-----------------------Average phase function for the downward scattered solar beam------------
                phasint = 0._jprb
!                       cosan_wcl(1:coef_scatt_ir%fmv_wcl_ph)  =                               &
!                       & cos(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph)*deg2rad)
                musat_d =  - 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                musat = raytracing_tl%pathsat(lay, ipf) / raytracing%pathsat(lay, ipf) ** 2
                musun = raytracing_tl%pathsun(lay, ipf) / raytracing%pathsun(lay, ipf) ** 2
                loop6 : DO k = 0, iend, iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt((1._jprb - musat_d ** 2) * (1._jprb - musun_d ** 2)) * cos(k * deg2rad)
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    scattangle = musat * musun_d + musun * musat_d
                  ELSE
                    scattangle = musat * musun_d + musun * musat_d -                                                     &
                      & musat * musat_d * sqrt(1._jprb - musun_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musat_d ** 2) &
                      &  -                                                                                               &
                      & musun * musun_d * sqrt(1._jprb - musat_d ** 2) * cos(k * deg2rad) / sqrt(1._jprb - musun_d ** 2)
                  ENDIF
!                   do kk=1,coef_scatt_ir%fmv_wcl_ph-1
!                     if(scattangle_d>=coef_scatt_ir%fmv_wcl_ph_val_cos(kk+1) &
!                       &.and.scattangle_d<=coef_scatt_ir%fmv_wcl_ph_val_cos(kk))then
                  ikk = max(1_jpim, int(acos(scattangle_d) / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_wcl_ph_val(ikk) - 1._jprb, jpim)
                  deltapwcl_d = (optp%optpwcl(ityp)%pha(ich1, 1, kk + 1) - optp%optpwcl(ityp)%pha(ich1, 1, kk))
                  deltawcl = (coef_scatt_ir%fmv_wcl_ph_val_cos(kk) - coef_scatt_ir%fmv_wcl_ph_val_cos(kk + 1))
                  phasint =  - scattangle * deltapwcl_d / deltawcl
!                       exit
!                     endif
!                   enddo
                  transmission_scatt_ir_tl%azphdocls(j, lay, ityp) =      &
                    & transmission_scatt_ir_tl%azphdocls(J, lay, ityp) + phasint * iang * ZPI
                ENDDO loop6
                transmission_scatt_ir_tl%AZPHDOTOT(J, lay) = transmission_scatt_ir_tl%AZPHDOTOT(J, lay) +              &
                  & transmission_scatt_ir_tl%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp) +  &
                  & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir_tl%OPDPSCLS(J, lay, ityp)
              ENDIF
!ityp
            ENDIF
!addsun
          ENDIF
!cldtyp
        ENDDO
!ncldtyp
        IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
          transmission_scatt_ir_tl%GPAR(J, lay) =                                               &
            & transmission_scatt_ir_tl%GPARTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay) -  &
            & transmission_scatt_ir_tl%OPDPS(J, lay) * transmission_scatt_ir%GPARTOT(J, lay) /  &
            & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
        ENDIF
        transmission_scatt_ir_tl%opdpcldla(j, lay) = transmission_scatt_ir_tl%opdpa(j, lay) +      &
          & transmission_scatt_ir_tl%opdps(j, lay) * transmission_scatt_ir%gpar(J, lay) +          &
          & transmission_scatt_ir_tl%gpar(j, lay) * transmission_scatt_ir%opdps(J, lay)
        OPDPCLDL(J, lay)                           =                                                        &
          & transmission_scatt_ir_tl%opdpcldla(j, lay) * raytracing%pathsat(lay, ipf) * coef%ff_gam(ICH) +  &
          & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%opdpcldla(j, lay) * coef%ff_gam(ICH)
        IF (sun(j)) THEN
          IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
            transmission_scatt_ir_tl%AZPHUP(J, lay) =                                               &
              & transmission_scatt_ir_tl%AZPHUPTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay) -  &
              & transmission_scatt_ir%AZPHUPTOT(J, lay) * transmission_scatt_ir_tl%OPDPS(J, lay) /  &
              & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
            transmission_scatt_ir_tl%AZPHDO(J, lay) =                                               &
              & transmission_scatt_ir_tl%AZPHDOTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay) -  &
              & transmission_scatt_ir%AZPHDOTOT(J, lay) * transmission_scatt_ir_tl%OPDPS(J, lay) /  &
              & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
          ENDIF
          OPDPCLDLSUN(J, lay) =                                                                                          &
            & transmission_scatt_ir_tl%OPDPCLDLA(J, lay) * (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) &
            &  * coef%ff_gam(ICH) +                                                                                      &
            & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%OPDPCLDLA(J, lay) * coef%ff_gam(ICH) +             &
            & raytracing_tl%pathsun(lay, ipf) * transmission_scatt_ir%OPDPCLDLA(J, lay) * coef%ff_gam(ICH)
        ENDIF
      ENDDO
!nlayers
    ENDDO
!nchannels
  ENDIF
!opts%addclouds
!-----Compute optical parameters for each stream--------------------------------
  transmission_scatt_ir_stream_tl%SSA = 0._jprb
  DO j = 1, nchannels
    ich = chanprof(J)%chan
    ipf = chanprof(J)%prof
    DO ist = 0, ircld%nstream(ipf)
      IF (ist == 0) THEN
        IF (opts%addaerosl) THEN
          OPD = 0._jprb
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPAC(IST, J, 1)    = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPACSUN(IST, J, 1) = 0._jprb
          DO lay = 1, nlayers
            lev = lay + 1
            IF (sun(j)) THEN
              transmission_scatt_ir_stream_tl%azphacup(ist, j, lay)   = transmission_scatt_ir_tl%azphaerupa(j, lay)
              transmission_scatt_ir_stream_tl%azphacdo(ist, j, lay)   = transmission_scatt_ir_tl%azphaerdoa(j, lay)
              transmission_scatt_ir_stream_tl%bcksp(ist, j, lay)      = transmission_scatt_ir_tl%gparaer(j, lay)
              transmission_scatt_ir_stream_tl%opdpabs(ist, j, lay)    = transmission_scatt_ir_tl%opdpaaer(j, lay) *      &
                & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +                     &
                & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%OPDPAAER(J, lay) * coef%ff_gam(ICH) +          &
                & raytracing_tl%pathsun(lay, ipf) * transmission_scatt_ir%OPDPAAER(J, lay) * coef%ff_gam(ICH)
              transmission_scatt_ir_stream_tl%opdpsca(ist, j, lay)    = transmission_scatt_ir_tl%opdpsaer(j, lay) *      &
                & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +                     &
                & raytracing_tl%pathsat(lay, ipf) * transmission_scatt_ir%opdpsaer(j, lay) * coef%ff_gam(ICH) +          &
                & raytracing_tl%pathsun(lay, ipf) * transmission_scatt_ir%opdpsaer(j, lay) * coef%ff_gam(ICH)
              transmission_scatt_ir_stream_tl%opdpaclsun(ist, j, lay) = opdpaerlsun(j, lay)
              opdsun = opdsun + transmission_scatt_ir_stream_tl%OPDPACLSUN(IST, J, lay)
              transmission_scatt_ir_stream_tl%OPDPACSUN(IST, J, lev)  = OPDSUN
            ENDIF
            transmission_scatt_ir_stream_tl%OPDPACL(IST, J, lay) = OPDPAERL(J, lay)
            OPD = OPD + transmission_scatt_ir_stream_tl%OPDPACL(IST, J, lay)
            transmission_scatt_ir_stream_tl%opdpac(ist, j, lev)  = OPD
          ENDDO
! layers
        ELSE IF ((.NOT. opts%addaerosl) .AND. (opts%addclouds)) THEN
          transmission_scatt_ir_stream_tl%opdpaclsun(ist, j, :) = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPACSUN(IST, J, :)  = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPACL(IST, J, :)    = 0._jprb
          transmission_scatt_ir_stream_tl%opdpac(ist, j, :)     = 0._jprb
          transmission_scatt_ir_stream_tl%azphacup(ist, j, :)   = 0._jprb
          transmission_scatt_ir_stream_tl%azphacdo(ist, j, :)   = 0._jprb
          IF (sun(j)) THEN
            transmission_scatt_ir_stream_tl%opdpabs(ist, j, :)    = 0._jprb
            transmission_scatt_ir_stream_tl%opdpsca(ist, j, :)    = 0._jprb
          ENDIF
        ENDIF
! opts%addaerosl
      ELSE
        IF (opts%addclouds) THEN
          opd = 0._jprb
          opdsun = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPAC(IST, J, 1)    = 0._jprb
          transmission_scatt_ir_stream_tl%OPDPACSUN(IST, J, 1) = 0._jprb
          DO lay = 1, nlayers
            lev = lay + 1
            IF (sun(j)) THEN
              IF ((ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +      &
                & transmission_scatt_ir%OPDPSAER(J, lay)) /= 0._jprb) THEN
                transmission_scatt_ir_stream_tl%AZPHACUP(IST, J, lay) =                                                   &
                  & transmission_scatt_ir_tl%AZPHAERUPA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHAERUPA(J, lay) / (              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHAERUPA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (                 &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHAERUPA(J, lay) *                &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                            &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 +                                                        &
                  & transmission_scatt_ir_tl%AZPHUP(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *                       &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHUP(J, lay) / (                                                              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHUP(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (                        &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHUP(J, lay) *                    &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                                &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_stream_tl%AZPHACDO(IST, J, lay) =                                                   &
                  & transmission_scatt_ir_tl%AZPHAERDOA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHAERDOA(J, lay) / (              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHAERDOA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (                 &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHAERDOA(J, lay) *                &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                            &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 +                                                        &
                  & transmission_scatt_ir_tl%AZPHDO(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *                       &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHDO(J, lay) / (                                                              &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%AZPHDO(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (                        &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%AZPHDO(J, lay) *                    &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                                &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_stream_tl%BCKSP(IST, J, lay)    =                                                   &
                  & transmission_scatt_ir_tl%GPARAER(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (                 &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%GPARAER(J, lay) / (                 &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%GPARAER(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (                    &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%GPARAER(J, lay) *                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                            &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 +                                                        &
                  & transmission_scatt_ir_tl%GPAR(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *                         &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%GPAR(J, lay) / (                                                                &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) -                                                             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                               &
                  & transmission_scatt_ir%GPAR(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (                          &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2 -                                                        &
                  & transmission_scatt_ir_tl%OPDPSAER(J, lay) * transmission_scatt_ir%GPAR(J, lay) *                      &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                                &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_stream_tl%OPDPABS(IST, J, lay)  =                                            &
                  & raytracing_tl%pathsat(lay, ipf) * coef%ff_gam(ICH) * transmission_scatt_ir%OPDPAAER(J, lay) +  &
                  & raytracing_tl%pathsat(lay, ipf) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) *            &
                  & transmission_scatt_ir%OPDPA(J, lay) + transmission_scatt_ir_tl%OPDPAAER(J, lay) *              &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +             &
                  & transmission_scatt_ir_tl%OPDPA(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                        &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +             &
                  & raytracing_tl%pathsun(lay, ipf) * coef%ff_gam(ICH) * transmission_scatt_ir%OPDPAAER(J, lay) +  &
                  & raytracing_tl%pathsun(lay, ipf) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) *            &
                  & transmission_scatt_ir%OPDPA(J, lay)
                transmission_scatt_ir_stream_tl%OPDPSCA(IST, J, lay)  =                                            &
                  & raytracing_tl%pathsat(lay, ipf) * coef%ff_gam(ICH) * transmission_scatt_ir%OPDPSAER(J, lay) +  &
                  & raytracing_tl%pathsat(lay, ipf) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) *            &
                  & transmission_scatt_ir%OPDPS(J, lay) + transmission_scatt_ir_tl%OPDPSAER(J, lay) *              &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +             &
                  & transmission_scatt_ir_tl%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) *                        &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH) +             &
                  & raytracing_tl%pathsun(lay, ipf) * coef%ff_gam(ICH) * transmission_scatt_ir%OPDPSAER(J, lay) +  &
                  & raytracing_tl%pathsun(lay, ipf) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) *            &
                  & transmission_scatt_ir%OPDPS(J, lay)
              ELSE
                transmission_scatt_ir_stream_tl%AZPHACUP(IST, J, lay) = 0._jprb
                transmission_scatt_ir_stream_tl%AZPHACDO(IST, J, lay) = 0._jprb
                transmission_scatt_ir_stream_tl%BCKSP(IST, J, lay)    = 0._jprb
                transmission_scatt_ir_stream_tl%OPDPSCA(IST, J, lay)  = 0._jprb
                transmission_scatt_ir_stream_tl%OPDPABS(IST, J, lay)  = 0._jprb
              ENDIF
              transmission_scatt_ir_stream_tl%OPDPACLSUN(IST, J, lay) =      &
                & OPDPCLDLSUN(J, lay) * ircld%ICLDARR(IST, lay, ipf) + OPDPAERLSUN(J, lay)
              OPDSUN = OPDSUN + transmission_scatt_ir_stream_tl%OPDPACLSUN(IST, J, lay)
              transmission_scatt_ir_stream_tl%OPDPACSUN(IST, J, lev)  = OPDSUN
            ENDIF
            transmission_scatt_ir_stream_tl%OPDPACL(IST, J, lay) =      &
              & OPDPCLDL(J, lay) * ircld%ICLDARR(IST, lay, ipf) + OPDPAERL(J, lay)
            OPD = OPD + transmission_scatt_ir_stream_tl%OPDPACL(IST, J, lay)
!                OPDSUN                   = OPDSUN                                      + &
!                          transmission_scatt_ir_stream_tl%OPDPACLSUN(IST,J,lay)
            transmission_scatt_ir_stream_tl%OPDPAC(IST, J, lev)  = OPD
!                transmission_scatt_ir_stream_tl%OPDPACSUN(IST,J,lev)  = OPDSUN
          ENDDO
! layers
        ENDIF
! opts%addclouds
      ENDIF
! ist=0
    ENDDO
! istream
  ENDDO
! channels/j
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdpscattir_tl
