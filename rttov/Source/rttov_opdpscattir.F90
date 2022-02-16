!!    Compute Optical parameters for aerosols and clouds
SUBROUTINE rttov_opdpscattir( &
            & nlayers,                      &
            & chanprof,                     &
            & opts,                         &
            & aux,                          &
            & profiles,                     &
            & sun,                          &
            & coef,                         &
            & coef_scatt_ir,                &
            & raytracing,                   &
            & transmission_scatt_ir,        &
            & transmission_scatt_ir_stream, &
            & optp,                         &
            & ircld)
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
!     1           30/7/2004    RTIASI-5. Marco Matricardi. ECMWF.
!     2           27/02/2009   Profile levels to include ToA. Distinguish
!                              arrays in raytracing and profiles (on levels)
!                              from all others (on layers) - size, index
!                              labels, looping (P. Rayer)
!     3           09/06/2009   Corrected bug - to now assign ircld%xpresave from
!                              user pressure array, not refprf (R.Saunders)
!     4           03/11/2009   Transmittances / optical depths on levels (A Geer)
!     5           02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     6           02/12/2009   Introduced multiple cloud types in a single layer. Pathsun, Pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     7           05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     8           14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                              should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
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
       & e00,     &
       & t00,     &
       & ti,      &
       & ncldtyp, &
       & max_sol_zen
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(profile_aux               ), INTENT(INOUT) :: aux
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type                ), INTENT(INOUT) :: ircld
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: j, ich   , ich1   , i, ipf, ish, istream, k, kk, ityp, iconf, iae, ikk
  INTEGER(KIND=jpim) :: lev       , lay   , lctyp
  INTEGER(KIND=jpim) :: iang      , iend  , icount , icl
  REAL   (KIND=jprb) :: opd       , opdsun
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: sum
  REAL   (KIND=jprb) :: sum1
  REAL   (KIND=jprb) :: sum2
  REAL   (KIND=jprb) :: sum3
  REAL   (KIND=jprb) :: scattangle
  REAL   (KIND=jprb) :: deltapaer
  REAL   (KIND=jprb) :: deltapice
  REAL   (KIND=jprb) :: deltapwcl
  REAL   (KIND=jprb) :: deltaaer
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: deltawcl
  REAL   (KIND=jprb) :: deltadg
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: zdelth
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: phasint
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
  REAL   (KIND=jprb) :: phup
  REAL   (KIND=jprb) :: phdo
!       Local arrays:
  REAL   (KIND=jprb) :: opdpaerl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpaerlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=JPRB) :: ZPI         , ztmpx , ztmpy
  INTEGER(KIND=jpim) :: nprofiles                                                             ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                                                             ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 0_jpim, ZHOOK_HANDLE)
  nprofiles   = size(profiles)
  nchannels   = size(chanprof)
  bparh       = 0._jprb
  opdpaerl    = 0._jprb
  opdpaerlsun = 0._jprb
  opdpcldl    = 0._jprb
  opdpcldlsun = 0._jprb
  ZPI         = deg2rad / (2._jprb * PI)
  iang        = 360_jpim / (jpazn + 1)
  iend        = iang * jpazn
!-----------------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-----------------------------------------------------------------------------------------
  IF (opts%addaerosl) THEN
!       opdpaerl(:,:)     =0._jprb
!       opdpaerlsun(:,:)  =0._jprb
    DO j = 1, nprofiles
      DO lay = 1, nlayers
        lev    = lay + 1
        icount = 0
        DO icl = 1, coef_scatt_ir%fmv_aer_comp
          IF (profiles(j)%aerosols(icl, lay) /= 0._jprb) THEN
            icount = icount + 1
            aux%iaertyp(icount, lay, j) = icl
          ENDIF
        ENDDO
        aux%iaernum(lay, j)    = icount
!-----Compute relative humidity-----------------------------------------------------------
        ircld%tave(lay, j)     = (profiles(j)%t(lev - 1) + profiles(j)%t(lev)) / 2._jprb
        ircld%wmixave(lay, j)  = (profiles(j)%q(lev - 1) + profiles(j)%q(lev)) / 2._jprb
        ircld%xpresave(lay, j) = (profiles(j)%p(lev - 1) + profiles(j)%p(lev)) / 2._jprb
! saturated vapour pressure
        ircld%esw(lay, j) = e00 * exp(17.502_jprb * (ircld%tave(lay, j) - T00) / (ircld%tave(lay, j) - 32.19_jprb))
        ircld%esi(lay, j) = E00 * exp(22.587_jprb * (ircld%tave(lay, j) - T00) / (ircld%tave(lay, j) + 0.7_jprb))
        IF (ircld%tave(lay, j) > t00) THEN
          ircld%ppv(lay, j) = ircld%esw(lay, j)! Water phase
        ELSE IF (ircld%tave(lay, j) > ti .AND. ircld%tave(lay, j) <= t00) THEN
          ircld%ppv(lay, j) =      &
            & ircld%esi(lay, j) + (ircld%esw(lay, j) - ircld%esi(lay, j)) * ((ircld%tave(lay, j) - ti) / (t00 - ti)) ** 2! Mixed phase
        ELSE IF (ircld%tave(lay, j) <= ti) THEN
          ircld%ppv(lay, j) = ircld%esi(lay, j)! Ice phase
        ENDIF
        ircld%ppv(lay, j)     = ircld%ppv(lay, j) / 100._jprb
! layer average relative humidity
        aux%relhum(lay, j)    = 100._jprb * ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb * ircld%xpresave(lay, j) /      &
          & (ircld%ppv(lay, j) * (0.622_jprb + ircld%wmixave(lay, j) * 1.e-6_jprb * 0.622_jprb))
        aux%relhumref(lay, j) = aux%relhum(lay, j)
        IF (aux%relhum(lay, j) > 99._jprb) THEN
          aux%relhum(lay, j) = 99._jprb
        ENDIF
!----------------------------------------------------------------------------------------
      ENDDO
! layers
    ENDDO
! profiles
  ENDIF
! opts%addaerosl
  IF (opts%addaerosl .OR. opts%addclouds) THEN
    transmission_scatt_ir%OPDPAAER(:,:)    = 0._jprb
    transmission_scatt_ir%OPDPSAER(:,:)    = 0._jprb
    transmission_scatt_ir%OPDPA(:,:)       = 0._jprb
!     transmission_scatt_ir%OPDPACLS(:,:,:)  = 0._jprb
    transmission_scatt_ir%OPDPS(:,:)       = 0._jprb
!     transmission_scatt_ir%OPDPSCLS(:,:,:)  = 0._jprb
    transmission_scatt_ir%GPARAER(:,:)     = 0._jprb
    transmission_scatt_ir%GPAR(:,:)        = 0._jprb
    transmission_scatt_ir%GPARTOT(:,:)     = 0._jprb
!     transmission_scatt_ir%GPARCLS(:,:,:)   = 0._jprb
    transmission_scatt_ir%OPDPAERLA(:,:)   = 0._jprb
    transmission_scatt_ir%GPARAERA(:,:)    = 0._jprb
    transmission_scatt_ir%azphaerup(:,:)   = 0._jprb
    transmission_scatt_ir%azphaerupa(:,:)  = 0._jprb
    transmission_scatt_ir%azphaerdo(:,:)   = 0._jprb
    transmission_scatt_ir%azphaerdoa(:,:)  = 0._jprb
    transmission_scatt_ir%azphupcls(:,:,:) = 0._jprb
    transmission_scatt_ir%azphdocls(:,:,:) = 0._jprb
    transmission_scatt_ir%azphuptot(:,:)   = 0._jprb
    transmission_scatt_ir%azphdotot(:,:)   = 0._jprb
  ENDIF
  IF (opts%addaerosl) THEN
    frequency : DO j = 1, nchannels
      ich  = chanprof(J)%chan
      ipf  = chanprof(J)%prof
      sum  = 0._jprb
      sum1 = 0._jprb
      sum2 = 0._jprb
      sum3 = 0._jprb
      layers : DO lay = 1, nlayers
        lev = lay + 1
        numaer : DO I = 1, aux%iaernum(lay, ipf)
          iae = aux%iaertyp(I, lay, ipf)
          IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
!
!---------------Interpolate scattering parameters to actual value of relative-------------
!               humidity.
            DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
              IF (aux%relhum(lay, ipf) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                & aux%relhum(lay, ipf) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                delth  = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                zdelth = 1.0_JPRB / delth
                frach  = (aux%relhum(lay, ipf) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                afac   = (optp%optpaer(iae)%abs(ich, k + 1) - optp%optpaer(iae)%abs(ich, k)) * zdelth
                sfac   = (optp%optpaer(iae)%sca(ich, k + 1) - optp%optpaer(iae)%sca(ich, k)) * zdelth
                gfac   = (optp%optpaer(iae)%bpr(ich, k + 1) - optp%optpaer(iae)%bpr(ich, k)) * zdelth
                absch  = optp%optpaer(iae)%abs(ich, k) + afac * frach
                scach  = optp%optpaer(iae)%sca(ich, k) + sfac * frach
                bparh  = optp%optpaer(iae)%bpr(ich, k) + gfac * frach
                IF (sun(j)) THEN
                  ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
                  PFAC(1:coef_scatt_ir%fmv_aer_ph)  = (optp%optpaer(iae)%PHA(ICH1, K + 1, 1:coef_scatt_ir%fmv_aer_ph)     &
                    &  - optp%optpaer(iae)%PHA(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph)) * zdelth
                  PHASH(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%PHA(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph) +      &
                    & PFAC(1:coef_scatt_ir%fmv_aer_ph) * FRACH
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ELSE
            absch = optp%optpaer(iae)%abs(ich, 1)
            scach = optp%optpaer(iae)%sca(ich, 1)
!!                EXTCH=EXTC(ICH,IAE,1)
            bparh = optp%optpaer(iae)%bpr(ich, 1)
            IF (sun(j)) THEN
              ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
              PHASH(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(ich1, 1, 1:coef_scatt_ir%fmv_aer_ph)
            ENDIF
          ENDIF
!-------------Compute optical parameters considering the contribution of------------------
!             all the aerosol components present in the layer
          IF (i == 1) THEN
! NB lev=lay+1, so lev labels lower level of the layer
            transmission_scatt_ir%OPDPAAER(J, lay) =      &
              & profiles(IPF)%aerosols(IAE, lay) * ABSCH * raytracing%LTICK(nlayers + 1 - lay, IPF)
! = (aerosol amount at bottom of layer) * (aerosol abs coeff for given layer RH) * (layer depth as height diff of two upwardly successive levels):
            transmission_scatt_ir%OPDPSAER(J, lay) =      &
              & profiles(IPF)%aerosols(IAE, lay) * SCACH * raytracing%LTICK(nlayers + 1 - lay, IPF)
!!                OPDPEAER(J,lay)=OPDPEAER(J,lay)+PROFAER(IAE,lay,IPF)                    * &
!!                                 EXTCH*LTICK(90-lay,IPF)
            transmission_scatt_ir%GPARAERA(J, lay) =      &
              & profiles(IPF)%aerosols(IAE, lay) * SCACH * BPARH * raytracing%LTICK(nlayers + 1 - lay, IPF)
          ELSE
            transmission_scatt_ir%OPDPAAER(J, lay) = transmission_scatt_ir%OPDPAAER(J, lay) +      &
              & profiles(IPF)%aerosols(IAE, lay) * ABSCH * raytracing%LTICK(nlayers + 1 - lay, IPF)
            transmission_scatt_ir%OPDPSAER(J, lay) = transmission_scatt_ir%OPDPSAER(J, lay) +      &
              & profiles(IPF)%aerosols(IAE, lay) * SCACH * raytracing%LTICK(nlayers + 1 - lay, IPF)
            transmission_scatt_ir%GPARAERA(J, lay) = transmission_scatt_ir%GPARAERA(J, lay) +      &
              & profiles(IPF)%aerosols(IAE, lay) * SCACH * BPARH * raytracing%LTICK(nlayers + 1 - lay, IPF)
          ENDIF
          IF (j == 1421) THEN
            sum  = sum + profiles(IPF)%aerosols(IAE, lay) * ABSCH * raytracing%LTICK(nlayers + 1 - lay, IPF)
            sum1 = sum1 + profiles(IPF)%aerosols(IAE, lay) * SCACH * raytracing%LTICK(nlayers + 1 - lay, IPF)
          ENDIF
!-------------If solar radiation is present,compute the azimuthally averaged--------------
!             value of the phase function for the given value of the viewing--------------
!             angle and solar zenith angle.
          IF (sun(j)) THEN
            ICH1    = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
!-------------Average phase function for the upward scattered solar beam------------------
            PHASINT = 0._jprb
!                 COSAN_aer(1:coef_scatt_ir%fmv_aer_ph)  =                               &
!                    & COS(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph)*deg2rad)
            MUSAT   = 1._jprb / raytracing%PATHSAT(lay, ipf)
            MUSUN   =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
!if(ich==8461)print*,raytracing%PATHSAT(lay,ipf),raytracing%PATHSUN(lay,ipf),lay
            ztmpx   = SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2))
            ztmpy   = 1.E0_jprb / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad)
            LOOP1 : DO k = 0, IEND, IANG
              SCATTANGLE = MUSAT * MUSUN + ztmpx * COS(K * deg2rad)
!                     DO kk=1,coef_scatt_ir%fmv_aer_ph-1
!                       IF(SCATTANGLE>=coef_scatt_ir%fmv_aer_ph_val_cos(KK+1) &
!                         &  .AND.SCATTANGLE<=coef_scatt_ir%fmv_aer_ph_val_cos(KK))THEN
              ikk        = max(1_jpim, int(acos(SCATTANGLE) * ztmpy, jpim))
              KK = INT (coef_scatt_ir%ifmv_aer_ph_val(ikk) - 1._jprb, jpim)
              DELTAPAER  = (PHASH(KK + 1) - PHASH(KK))
              DELTAAER   = (coef_scatt_ir%fmv_aer_ph_val_cos(KK) - coef_scatt_ir%fmv_aer_ph_val_cos(KK + 1))
              PHUP       = PHASH(KK) + DELTAPAER * (coef_scatt_ir%fmv_aer_ph_val_cos(KK) - SCATTANGLE) / DELTAAER
!                         EXIT
!                       ENDIF
!                     ENDDO
              PHASINT    = PHASINT + PHUP * IANG * ZPI
            ENDDO LOOP1
            transmission_scatt_ir%PHASINTUPREF(J, lay, I) = PHASINT
            transmission_scatt_ir%AZPHAERUP(J, lay)       = transmission_scatt_ir%AZPHAERUP(J, lay) +      &
              & profiles(ipf)%aerosols(IAE, lay) * PHASINT * SCACH * raytracing%LTICK(nlayers + 1 - lay, ipf)
!-------------Average phase function for the downward scattered solar beam----------------
            PHASINT = 0._jprb
!                 COSAN_aer(1:coef_scatt_ir%fmv_aer_ph)  =                               &
!                    & COS(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph)*deg2rad)
            MUSAT =  - 1._jprb / raytracing%PATHSAT(lay, ipf)
            MUSUN =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
            ztmpx = SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2))
            ztmpy = 1.E0_jprb / (coef_scatt_ir%fmv_aer_ph_val_min * deg2rad)
            LOOP2 : DO k = 0, IEND, IANG
              SCATTANGLE = MUSAT * MUSUN + ztmpx * COS(K * deg2rad)
!                     DO kk=1,coef_scatt_ir%fmv_aer_ph-1
!                       IF(SCATTANGLE>=coef_scatt_ir%fmv_aer_ph_val_cos(KK+1) &
!                         & .AND.SCATTANGLE<=coef_scatt_ir%fmv_aer_ph_val_cos(KK))THEN
              ikk        = max(1_jpim, int(acos(SCATTANGLE) * ztmpy, jpim))
              KK = INT (coef_scatt_ir%ifmv_aer_ph_val(ikk) - 1._jprb, jpim)
              DELTAPAER  = (PHASH(KK + 1) - PHASH(KK))
              DELTAAER   = (coef_scatt_ir%fmv_aer_ph_val_cos(KK) - coef_scatt_ir%fmv_aer_ph_val_cos(KK + 1))
              PHDO       = PHASH(KK) + DELTAPAER * (coef_scatt_ir%fmv_aer_ph_val_cos(KK) - SCATTANGLE) / DELTAAER
!                         EXIT
!                       ENDIF
!                     ENDDO
              PHASINT    = PHASINT + PHDO * IANG * ZPI
            ENDDO LOOP2
            transmission_scatt_ir%PHASINTDOREF(J, lay, I) = PHASINT
            transmission_scatt_ir%AZPHAERDO(J, lay)       = transmission_scatt_ir%AZPHAERDO(J, lay) +      &
              & profiles(ipf)%aerosols(IAE, lay) * PHASINT * SCACH * raytracing%LTICK(nlayers + 1 - lay, ipf)
          ENDIF
        ENDDO numaer
      ENDDO layers
!---------Compute final values for optical parameters-------------------------------------
      DO lay = 1, nlayers
        lev = lay + 1
        IF (transmission_scatt_ir%OPDPSAER(J, lay) /= 0._JPRB) THEN
          transmission_scatt_ir%GPARAER(J, lay) =      &
            & transmission_scatt_ir%GPARAERA(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
          IF (sun(j)) THEN
            transmission_scatt_ir%AZPHAERUPA(J, lay) =      &
              & transmission_scatt_ir%AZPHAERUP(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
            transmission_scatt_ir%AZPHAERDOA(J, lay) =      &
              & transmission_scatt_ir%AZPHAERDO(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
          ENDIF
        ENDIF
        transmission_scatt_ir%OPDPAERLA(J, lay) = transmission_scatt_ir%OPDPAAER(J, lay) +      &
          & transmission_scatt_ir%OPDPSAER(J, lay) * transmission_scatt_ir%GPARAER(J, lay)
        OPDPAERL(J, lay)                        =      &
          & transmission_scatt_ir%OPDPAERLA(J, lay) * RAYTRACING%PATHSAT(lay, IPF) * coef%ff_gam(ICH)
        IF (sun(j)) THEN
          OPDPAERLSUN(J, lay) =                                                                                          &
            & transmission_scatt_ir%OPDPAERLA(J, lay) * (RAYTRACING%PATHSUN(lay, IPF) + RAYTRACING%PATHSAT(lay, IPF)) *  &
            & coef%ff_gam(ICH)
        ENDIF
      ENDDO
! layers
!        if(j==1421)print*,sum+sum1,'EXTOPD_AER'
    ENDDO frequency
  ENDIF
! opts%addaerosl
!-----------------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-----------------------------------------------------------------------------------------
  IF (opts%addclouds) THEN
!          opdpcldl(:,:)      =0._jprb
!          opdpcldlsun(:,:)   =0._jprb
    DO j = 1, nchannels
      ich  = chanprof(j)%chan
      ipf  = chanprof(j)%prof
      ish  = profiles(ipf)%ish
      sum  = 0._jprb
      sum1 = 0._jprb
      sum2 = 0._jprb
      sum3 = 0._jprb
!            IOFF=NUMAE
      DO lay = 1, nlayers
        lev = lay + 1
        transmission_scatt_ir%azphup(j, lay) = 0._jprb
        transmission_scatt_ir%azphdo(j, lay) = 0._jprb
        transmission_scatt_ir%OPDPA(j, lay)  = 0._jprb
        transmission_scatt_ir%OPDPS(j, lay)  = 0._jprb
        transmission_scatt_ir%GPAR(j, lay)   = 0._jprb
!               transmission_scatt_ir%OPDPAAER(j,lay)=0._jprb
!               transmission_scatt_ir%OPDPSAER(j,lay)=0._jprb
!               transmission_scatt_ir%GPARAER(j,lay) =0._jprb
        DO lctyp = 1, ncldtyp
          IF (ircld%CLDTYP(lctyp, lay, ipf) /= 0) THEN
            ITYP  = ircld%CLDTYP(lctyp, lay, ipf)
            ICONF = ITYP
!                  ITMP = ITEMP  (lay,IPF)
!                  ICONF= CLDTYP (lay,IPF)+ITEMP(lay,IPF)-1
!---------------Compute cloud  optical parameters ----------------------------------------
            IF (ityp <= 5_jpim) THEN
!-----------------------------------------------------------------------------------------
!                 For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
              transmission_scatt_ir%OPDPACLS(J, lay, ityp) =                                                       &
                & profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%abs(ich, 1) *  &
                & raytracing%LTICK(nlayers + 1 - lay, IPF)
              transmission_scatt_ir%OPDPA(J, lay)          =      &
                & transmission_scatt_ir%OPDPA(J, lay) + transmission_scatt_ir%OPDPACLS(J, lay, ityp)
              transmission_scatt_ir%OPDPSCLS(J, lay, ityp) =  +                                                    &
                & profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%sca(ich, 1) *  &
                & raytracing%LTICK(nlayers + 1 - lay, IPF)
              transmission_scatt_ir%OPDPS(J, lay)          =      &
                & transmission_scatt_ir%OPDPS(J, lay) + transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              sum2 = sum2 +                                                                                        &
                & profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%abs(ich, 1) *  &
                & raytracing%LTICK(nlayers + 1 - lay, IPF)
              sum3 = sum3 +                                                                                        &
                & profiles(IPF)%cloud(Ityp, lay) * coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%sca(ich, 1) *  &
                & raytracing%LTICK(nlayers + 1 - lay, IPF)
!                  OPDPE(J,jl)=PROFCLD(ITYP,jl,IPF)*coef_scatt_ir%CONFAC(ICONF)          * &
!                                EXTC(ICH,IOFF+ITYP,ITMP)*LTICK(90-jl,IPF)
              transmission_scatt_ir%GPARCLS(J, lay, ityp)  = optp%optpwcl(ityp)%bpr(ich, 1)
              transmission_scatt_ir%GPARTOT(J, lay)        = transmission_scatt_ir%GPARTOT(J, lay) +      &
                & transmission_scatt_ir%GPARCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
            ELSE
!-----------------------------------------------------------------------------------------
!                 For ice clouds optical parameters are computed using regression
!                 coefficients
!-----------------------------------------------------------------------------------------
              transmission_scatt_ir%OPDPACLS(J, lay, ityp) =  + profiles(IPF)%cloud(Ityp, lay) * (                          &
                & optp%optpicl(ish)%abs(ich, 1) + optp%optpicl(ish)%abs(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%abs(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%abs(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF)
              transmission_scatt_ir%OPDPA(J, lay)          =      &
                & transmission_scatt_ir%OPDPA(J, lay) + transmission_scatt_ir%OPDPACLS(J, lay, ityp)
              transmission_scatt_ir%OPDPSCLS(J, lay, ityp) =  + profiles(IPF)%cloud(Ityp, lay) * (                          &
                & optp%optpicl(ish)%sca(ich, 1) + optp%optpicl(ish)%sca(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%sca(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%sca(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF)
              transmission_scatt_ir%OPDPS(J, lay)          =      &
                & transmission_scatt_ir%OPDPS(J, lay) + transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              transmission_scatt_ir%GPARCLS(J, lay, ityp)  =  + (                                     &
                & optp%optpicl(ish)%bpr(ich, 1) + optp%optpicl(ish)%bpr(ich, 2) * aux%dg(lay, ipf) +  &
                & optp%optpicl(ish)%bpr(ich, 3) * aux%dg(lay, ipf) ** 2 +                             &
                & optp%optpicl(ish)%bpr(ich, 4) * aux%dg(lay, ipf) ** 3)
              transmission_scatt_ir%GPARTOT(J, lay)        = transmission_scatt_ir%GPARTOT(J, lay) +      &
                & transmission_scatt_ir%GPARCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
            ENDIF
!---------------If solar radiation is present,compute the azimuthally averaged------------
!               value of the phase function for the given value of the viewing------------
!               angle and solar zenith angle.
!-----------------------------------------------------------------------------------------
            IF (sun(j)) THEN
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
              IF (ITYP == 6_jpim) THEN
                ich1 = ich - coef_scatt_ir%fmv_icl_pha_ioff + 1
                DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                  IF (aux%dg(lay, ipf) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                    & aux%dg(lay, ipf) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                    deltap     = (optp%optpicl(ish)%PHA(ICH1, k + 1, :) - optp%optpicl(ish)%PHA(ICH1, k, :))
                    DELTADG    = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                    phasice(:) = optp%optpicl(ish)%PHA(ICH1, k, :) +      &
                      & deltap(:) * (aux%dg(lay, ipf) - coef_scatt_ir%fmv_icl_dg(k, ish)) / DELTADG
                    EXIT
                  ENDIF
                ENDDO
!-----------------------Average phase function for the upward scattered solar beam----------
                PHASINT = 0._jprb
!                       COSAN_icl(1:coef_scatt_ir%fmv_icl_ph)  =                           &
!                       & COS(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph)*deg2rad)
                MUSAT   = 1._jprb / raytracing%PATHSAT(lay, ipf)
                MUSUN   =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
                LOOP3 : DO K = 0, IEND, IANG
                  SCATTANGLE =      &
                    & MUSAT * MUSUN + SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2)) * COS(K * deg2rad)
!                         DO KK=1,coef_scatt_ir%fmv_icl_ph-1
!                           IF(SCATTANGLE>=coef_scatt_ir%fmv_icl_ph_val_cos(KK+1) &
!                             &  .AND.SCATTANGLE<=coef_scatt_ir%fmv_icl_ph_val_cos(KK))THEN
                  ikk = max(1_jpim, int(acos(SCATTANGLE) / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_icl_ph_val(ikk) - 1._jprb, jpim)
                  DELTAPICE = (phasice(KK + 1) - PHAsice(KK))
                  DELTAICE = (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - coef_scatt_ir%fmv_icl_ph_val_cos(KK + 1))
                  PHASINT = PHAsice(KK) + DELTAPICE * (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - SCATTANGLE) / DELTAICE
!                             EXIT
!                           ENDIF
!                         ENDDO
                  transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) =      &
                    & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) + PHASINT * IANG * ZPI
                ENDDO LOOP3
                transmission_scatt_ir%AZPHUPTOT(J, lay) = transmission_scatt_ir%AZPHUPTOT(J, lay) +      &
                  & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
!---------------------Average phase function for the downward scattered solar beam--------
                PHASINT = 0._jprb
!                       COSAN_icl(1:coef_scatt_ir%fmv_icl_ph)  =                           &
!                       & COS(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph)*deg2rad)
                MUSAT =  - 1._jprb / raytracing%PATHSAT(lay, ipf)
                MUSUN =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
                LOOP4 : DO K = 0, IEND, IANG
                  SCATTANGLE =      &
                    & MUSAT * MUSUN + SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2)) * COS(K * deg2rad)
!                       DO KK=1,coef_scatt_ir%fmv_icl_ph-1
!                         IF(SCATTANGLE>=coef_scatt_ir%fmv_icl_ph_val_cos(KK+1) &
!                           &  .AND.SCATTANGLE<=coef_scatt_ir%fmv_icl_ph_val_cos(KK))THEN
                  ikk = max(1_jpim, int(acos(SCATTANGLE) / (coef_scatt_ir%fmv_icl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_icl_ph_val(ikk) - 1._jprb, jpim)
                  DELTAPICE = (PHAsice(KK + 1) - PHAsice(KK))
                  DELTAICE = (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - coef_scatt_ir%fmv_icl_ph_val_cos(KK + 1))
                  PHASINT = PHAsice(KK) + DELTAPICE * (coef_scatt_ir%fmv_icl_ph_val_cos(KK) - SCATTANGLE) / DELTAICE
!                           EXIT
!                         ENDIF
!                       ENDDO
                  transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) =      &
                    & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) + PHASINT * IANG * ZPI
                ENDDO LOOP4
                transmission_scatt_ir%AZPHDOTOT(J, lay) = transmission_scatt_ir%AZPHDOTOT(J, lay) +      &
                  & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              ELSE
!-----------------------------------------------------------------------------------------
!                     Water clouds
!-----------------------------------------------------------------------------------------
!---------------------Average phase function for the upward scattered solar beam----------
                ich1    = ich - coef_scatt_ir%fmv_wcl_pha_ioff + 1
                PHASINT = 0._jprb
!                       COSAN_wcl(1:coef_scatt_ir%fmv_wcl_ph)  =                           &
!                       & COS(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph)*deg2rad)
                MUSAT   = 1._jprb / raytracing%PATHSAT(lay, ipf)
                MUSUN   =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
                LOOP5 : DO K = 0, IEND, IANG
                  SCATTANGLE =      &
                    & MUSAT * MUSUN + SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2)) * COS(K * deg2rad)
!                       DO KK=1,coef_scatt_ir%fmv_wcl_ph-1
!                         IF(SCATTANGLE>=coef_scatt_ir%fmv_wcl_ph_val_cos(KK+1) &
!                           & .AND.SCATTANGLE<=coef_scatt_ir%fmv_wcl_ph_val_cos(KK))THEN
                  ikk = max(1_jpim, int(acos(SCATTANGLE) / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_wcl_ph_val(ikk) - 1._jprb, jpim)
                  DELTAPWCL = (optp%optpwcl(ityp)%PHA(ICH1, 1, KK + 1) - optp%optpwcl(ityp)%PHA(ICH1, 1, KK))
                  DELTAWCL = (coef_scatt_ir%fmv_wcl_ph_val_cos(KK) - coef_scatt_ir%fmv_wcl_ph_val_cos(KK + 1))
                  PHASINT = optp%optpwcl(ityp)%PHA(ICH1, 1, KK) +      &
                    & DELTAPWCL * (coef_scatt_ir%fmv_wcl_ph_val_cos(KK) - SCATTANGLE) / DELTAWCL
!                         EXIT
!                         ENDIF
!                       ENDDO
                  transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) =      &
                    & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) + PHASINT * IANG * ZPI
                ENDDO LOOP5
                transmission_scatt_ir%AZPHUPTOT(J, lay) = transmission_scatt_ir%AZPHUPTOT(J, lay) +      &
                  & transmission_scatt_ir%AZPHUPCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
!---------------------Average phase function for the downward scattered solar beam--------
                PHASINT = 0._jprb
!                       COSAN_wcl(1:coef_scatt_ir%fmv_wcl_ph)  =                           &
!                       & COS(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph)*deg2rad)
                MUSAT =  - 1._jprb / raytracing%PATHSAT(lay, ipf)
                MUSUN =  - 1._jprb / raytracing%PATHSUN(lay, ipf)
                LOOP6 : DO K = 0, IEND, IANG
                  SCATTANGLE =      &
                    & MUSAT * MUSUN + SQRT((1._jprb - MUSAT ** 2) * (1._jprb - MUSUN ** 2)) * COS(K * deg2rad)
!                       DO KK=1,coef_scatt_ir%fmv_wcl_ph-1
!                         IF(SCATTANGLE>=coef_scatt_ir%fmv_wcl_ph_val_cos(KK+1) &
!                           &   .AND.SCATTANGLE<=coef_scatt_ir%fmv_wcl_ph_val_cos(KK))THEN
                  ikk = max(1_jpim, int(acos(SCATTANGLE) / (coef_scatt_ir%fmv_wcl_ph_val_min * deg2rad), jpim))
                  KK = INT (coef_scatt_ir%ifmv_wcl_ph_val(ikk) - 1._jprb, jpim)
                  DELTAPWCL = (optp%optpwcl(ityp)%PHA(ICH1, 1, KK + 1) - optp%optpwcl(ityp)%PHA(ICH1, 1, KK))
                  DELTAWCL = (coef_scatt_ir%fmv_wcl_ph_val_cos(KK) - coef_scatt_ir%fmv_wcl_ph_val_cos(KK + 1))
                  PHASINT = optp%optpwcl(ityp)%PHA(ICH1, 1, KK) +      &
                    & DELTAPWCL * (coef_scatt_ir%fmv_wcl_ph_val_cos(KK) - SCATTANGLE) / DELTAWCL
!                           EXIT
!                         ENDIF
!                       ENDDO
                  transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) =      &
                    & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) + PHASINT * IANG * ZPI
                ENDDO LOOP6
                transmission_scatt_ir%AZPHDOTOT(J, lay) = transmission_scatt_ir%AZPHDOTOT(J, lay) +      &
                  & transmission_scatt_ir%AZPHDOCLS(J, lay, ityp) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
          transmission_scatt_ir%GPAR(J, lay) =      &
            & transmission_scatt_ir%GPARTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
        ENDIF
        transmission_scatt_ir%OPDPCLDLA(J, lay) =      &
          & transmission_scatt_ir%OPDPA(J, lay) + transmission_scatt_ir%GPAR(J, lay) * transmission_scatt_ir%OPDPS(J, lay)
        OPDPCLDL(J, lay)                        =      &
          & transmission_scatt_ir%OPDPCLDLA(J, lay) * raytracing%PATHSAT(lay, ipf) * coef%ff_gam(ICH)
        IF (sun(j)) THEN
          IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
            transmission_scatt_ir%AZPHUP(J, lay) =      &
              & transmission_scatt_ir%AZPHUPTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
            transmission_scatt_ir%AZPHDO(J, lay) =      &
              & transmission_scatt_ir%AZPHDOTOT(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
          ENDIF
          OPDPCLDLSUN(J, lay) =                                                                                          &
            & transmission_scatt_ir%OPDPCLDLA(J, lay) * (raytracing%PATHSAT(lay, ipf) + raytracing%PATHSUN(lay, ipf)) *  &
            & coef%ff_gam(ICH)
        ENDIF
      ENDDO
!      if(j==1421)print*,sum2+sum3,'EXTOPD_CLD'
    ENDDO
  ENDIF
!-----Compute optical parameters for each stream------------------------------------------
!      if(opts%addaerosl.or.opts%addclouds )then
  transmission_scatt_ir_stream%SSA = 0._jprb
  DO J = 1, nchannels
    ICH = chanprof(J)%chan
!opd=0
    IPF = chanprof(J)%prof
    DO ISTREAM = 0, ircld%NSTREAM(IPF)
      IF (ISTREAM == 0) THEN
!              OPD=0._jprb
!              OPDSUN=0._jprb
        IF (opts%addaerosl) THEN
          OPD = 0._jprb
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream%OPDPAC(ISTREAM, J, 1)    = 0._jprb
          transmission_scatt_ir_stream%OPDPACSUN(ISTREAM, J, 1) = 0._jprb
          DO lay = 1, nlayers
            lev = lay + 1
            IF (sun(j)) THEN
              transmission_scatt_ir_stream%AZPHACUP(ISTREAM, J, lay)   = transmission_scatt_ir%AZPHAERUPA(J, lay)
              transmission_scatt_ir_stream%AZPHACDO(ISTREAM, J, lay)   = transmission_scatt_ir%AZPHAERDOA(J, lay)
              transmission_scatt_ir_stream%BCKSP(ISTREAM, J, lay)      = transmission_scatt_ir%GPARAER(J, lay)
              transmission_scatt_ir_stream%OPDPABS(ISTREAM, J, lay)    =                                                 &
                & transmission_scatt_ir%OPDPAAER(J, lay) * (raytracing%PATHSAT(lay, ipf) + raytracing%PATHSUN(lay, ipf)) &
                &  * coef%ff_gam(ICH)
              transmission_scatt_ir_stream%OPDPSCA(ISTREAM, J, lay)    =                                                 &
                & transmission_scatt_ir%OPDPSAER(J, lay) * (raytracing%PATHSAT(lay, ipf) + raytracing%PATHSUN(lay, ipf)) &
                &  * coef%ff_gam(ICH)
              transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM, J, lay) = OPDPAERLSUN(J, lay)
              OPDSUN = OPDSUN + transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM, J, lay)
              transmission_scatt_ir_stream%OPDPACSUN(ISTREAM, J, lev)  = OPDSUN
            ENDIF
            transmission_scatt_ir_stream%OPDPACL(ISTREAM, J, lay) = OPDPAERL(J, lay)
            OPD = OPD + transmission_scatt_ir_stream%OPDPACL(ISTREAM, J, lay)
            transmission_scatt_ir_stream%OPDPAC(ISTREAM, J, lev)  = OPD
          ENDDO
! layers
        ELSE IF ((.NOT. opts%addaerosl) .AND. (opts%addclouds)) THEN
          transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM, J, :) = 0._jprb
          transmission_scatt_ir_stream%OPDPACSUN(ISTREAM, J, :)  = 0._jprb
          transmission_scatt_ir_stream%OPDPACL(ISTREAM, J, :)    = 0._jprb
          transmission_scatt_ir_stream%OPDPAC(ISTREAM, J, :)     = 0._jprb
          transmission_scatt_ir_stream%AZPHACUP(ISTREAM, J, :)   = 0._jprb
          transmission_scatt_ir_stream%AZPHACDO(ISTREAM, J, :)   = 0._jprb
          IF (sun(j)) THEN
            transmission_scatt_ir_stream%OPDPABS(ISTREAM, J, :)    = 0._jprb
            transmission_scatt_ir_stream%OPDPSCA(ISTREAM, J, :)    = 0._jprb
          ENDIF
        ENDIF
! opts%addaerosl
      ELSE
        IF (opts%addclouds) THEN
          OPD = 0._jprb
          OPDSUN = 0._jprb
          transmission_scatt_ir_stream%OPDPAC(ISTREAM, J, 1)    = 0._jprb
          transmission_scatt_ir_stream%OPDPACSUN(ISTREAM, J, 1) = 0._jprb
          DO lay = 1, nlayers
            lev = lay + 1
            IF (sun(j)) THEN
              IF ((ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +      &
                & transmission_scatt_ir%OPDPSAER(J, lay)) /= 0) THEN
                transmission_scatt_ir_stream%AZPHACUP(ISTREAM, J, lay) =                                  &
                  & transmission_scatt_ir%AZPHAERUPA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +              &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                             &
                  & transmission_scatt_ir%AZPHUP(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *          &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) / (                                                  &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +              &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_stream%AZPHACDO(ISTREAM, J, lay) =                                  &
                  & transmission_scatt_ir%AZPHAERDOA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +              &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                             &
                  & transmission_scatt_ir%AZPHDO(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *          &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) / (                                                  &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +              &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_stream%BCKSP(ISTREAM, J, lay)    =                               &
                  & transmission_scatt_ir%GPARAER(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +           &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) +                                          &
                  & transmission_scatt_ir%GPAR(J, lay) * transmission_scatt_ir%OPDPS(J, lay) *         &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) / (                                               &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +           &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_stream%OPDPABS(ISTREAM, J, lay)  = (transmission_scatt_ir%OPDPAAER(J, lay) +      &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPA(J, lay)) *                           &
                  & (raytracing%PATHSAT(lay, ipf) + raytracing%PATHSUN(lay, ipf)) * coef%ff_gam(ICH)
                transmission_scatt_ir_stream%OPDPSCA(ISTREAM, J, lay)  = (transmission_scatt_ir%OPDPSAER(J, lay) +      &
                  & ircld%ICLDARR(ISTREAM, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay)) *                           &
                  & (raytracing%PATHSAT(lay, ipf) + raytracing%PATHSUN(lay, ipf)) * coef%ff_gam(ICH)
              ELSE
                transmission_scatt_ir_stream%AZPHACUP(ISTREAM, J, lay) = 0._jprb
                transmission_scatt_ir_stream%AZPHACDO(ISTREAM, J, lay) = 0._jprb
                transmission_scatt_ir_stream%BCKSP(ISTREAM, J, lay)    = 0._jprb
                transmission_scatt_ir_stream%OPDPABS(ISTREAM, J, lay)  = 0._jprb
                transmission_scatt_ir_stream%OPDPSCA(ISTREAM, J, lay)  = 0._jprb
              ENDIF
              transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM, J, lay) =      &
                & OPDPCLDLSUN(J, lay) * ircld%ICLDARR(ISTREAM, lay, ipf) + OPDPAERLSUN(J, lay)
              OPDSUN = OPDSUN + transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM, J, lay)
              transmission_scatt_ir_stream%OPDPACSUN(ISTREAM, J, lev)  = OPDSUN
            ENDIF
            transmission_scatt_ir_stream%OPDPACL(ISTREAM, J, lay) =      &
              & OPDPCLDL(J, lay) * ircld%ICLDARR(ISTREAM, lay, ipf) + OPDPAERL(J, lay)
!                transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM,J,lay) =             &
!                   OPDPCLDLSUN(J,lay)*ircld%ICLDARR(ISTREAM,lay,ipf)+OPDPAERLSUN(J,lay)
            OPD = OPD + transmission_scatt_ir_stream%OPDPACL(ISTREAM, J, lay)
!                OPDSUN= OPDSUN+transmission_scatt_ir_stream%OPDPACLSUN(ISTREAM,J,lay)
            transmission_scatt_ir_stream%OPDPAC(ISTREAM, J, lev)  = OPD
!                transmission_scatt_ir_stream%OPDPACSUN(ISTREAM,J,lev)  = OPDSUN
          ENDDO
! layers
        ENDIF
! opts%addclouds
      ENDIF
! istream
!if(j==1421)print*,opd,istream
    ENDDO
! istream
  ENDDO
! channels/j
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR', 1_jpim, ZHOOK_HANDLE)
!      endif
END SUBROUTINE rttov_opdpscattir
