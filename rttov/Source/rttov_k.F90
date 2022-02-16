!
SUBROUTINE rttov_k( &
            & errorstatus,      &
            & chanprof,         &
            & opts,             &
            & profiles,         &
            & profiles_k,       &
            & coefs,            &
            & calcemis,         &
            & emissivity,       &
            & emissivity_k,     &
            & emissivity_out,   &
            & emissivity_out_k, &
            & transmission,     &
            & transmission_k,   &
            & radiancedata,     &
            & radiancedata_k,   &
            & traj,             &
            & traj_k,           &
            & pccomp,           &
            & pccomp_k,         &
            & profiles_k_pc,    &
            & profiles_k_rec,   &
            & channels_rec)
!
! Description:
! K matrix of rttov_direct
! to compute multi-channel level to space transmittances,
! top of atmosphere and level to space radiances and brightness
! temperatures and optionally surface emissivities, for many
! profiles in a single call, for satellite
! infrared or microwave sensors. The code requires a coefficient file
! for each sensor for which simulated radiances are requested.
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
! Method: The methodology is described in the following:
!
! Eyre J.R. and H.M. Woolf  1988 Transmittance of atmospheric gases
! in the microwave region: a fast model. Applied Optics 27  3244-3249
!
! Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
! systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
! librarian at ECMWF).
!
! Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
! Transfer Model for Assimilation of Satellite Radiance Observations.
! QJRMS, 125, 1407-1425.
!
! Matricardi, M., F. Chevallier and S. Tjemkes 2001 An improved general
! fast radiative transfer model for the assimilation of radiance
! observations. ECMWF Research Dept. Tech. Memo. 345
! (available from the librarian at ECMWF).
!
! Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
! transfer model for the infrared atmospheric sounding interferometer.
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
!
! Matricardi, M. 2009: An Observation operator for the assimilation of
! principal component scores into a NWP system. Available from EUMETSAT
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.1   01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.2   02/01/2003  More comments added (R Saunders)
!  1.3   24/01/2003  Error return code by input profile (P Brunel)
!                    Add WV Continuum and CO2 capability
!  1.4   02/06/2004  Change tests on id_comp_lvl == 7 by tests on
!                    fmv_model_ver (P. Brunel)
!  1.5   17/02/2005  Changed to allow calls from RTTOV_SCATT_K and
!                    RTTOV_CLD_K.   (A. Collard)
!  1.6   24/08/2005  More changes so routine still works for clear
!                    sky RTTOV (R.Saunders)
!  1.7   08/12/2005  Added surface humidity to lowest mean layer q (R Saunders)
!  1.8   13/03/2006  Marco Matricardi (ECMWF):
!              --    IASI capability added.
!              --    Variable trace gases added for IASI and AIRS.
!              --    Solar radiation added for IASI and AIRS
!              --    Altitude dependent local zenith angle added.
!              --    Linear in tau approximation for RT equation added.
!  1.9   08/02/2007  Removed polarisation index (Roger Saunders)
!  1.10  11/10/2007  Profile allocation though rttov_alloc_prof
!                    Move iaernum & iaertyp members to profile_aux P.Marguinaud
!  1.11   16/01/2008 Facility to apply regression limits (N. Bormann)
!  1.12   26/09/2008 Nullify temp pointers; this is for the NEC
!                    under -Chopt. P. Marguinaud
!
!  1.12   03/10/2008 Fix bug in checking addsolar (A. Geer)
!  1.13   08/12/2008 Fixed bug in addcosmic for polarim (R Saunders)
!  1.14   15/09/2009 User defined ToA. Layers distinct from levels (P.Rayer)
!  1.15   03/11/2009 Transmittances on levels (A Geer)
!  1.16   02/12/2009 Introduced principal component capability. Marco Matricardi. ECMWF
!  1.17   17/06/2010 Introduced optional spacetop flag to zero opdeps at user's
!                    model-top in Zeeman channels (P Rayer)
!  1.18   05/07/2010 Remove addsolar flag from profiles structure (J Hocking)
!  1.19   14/10/2010 Remove rt8_mode (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
! A user guide and technical documentation is available at
! http://www.metoffice.com/research/interproj/nwpsaf/rtm/index.html
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & transmission_Type, &
       & profile_Type,      &
       & radiance_Type,     &
       & rttov_chanprof,    &
       & rttov_traj
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw, &
       & sensor_id_ir, &
       & sensor_id_po, &
       & sensor_id_hi, &
       & ncldtyp
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_types, ONLY : &
       & rttov_traj_sta,    &
       & rttov_traj_dyn
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_k(size(chanprof))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis        (size(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  REAL   (KIND=jprb)     , INTENT(IN)                           :: emissivity      (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_k    (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(OUT)                          :: emissivity_out  (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_out_k(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                    ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_k                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                    ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_k
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj          , traj_k          ! Target is needed here (see rttov_check_temp)
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_k
  TYPE(profile_Type     ), OPTIONAL     , INTENT(INOUT)         :: profiles_k_pc(:)
  TYPE(profile_Type     ), OPTIONAL     , INTENT(INOUT)         :: profiles_k_rec(:)
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_checkinput_k.h"
#include "rttov_setgeometry_k.h"
#include "rttov_calcemis_ir_k.h"
#include "rttov_profaux_k.h"
#include "rttov_setpredictors_7_k.h"
#include "rttov_setpredictors_8_k.h"
#include "rttov_setpredictors_9_k.h"
#include "rttov_opdep_k.h"
#include "rttov_opdep_9_k.h"
#include "rttov_opdep_9_solar_k.h"
#include "rttov_transmit_k.h"
#include "rttov_transmit_9_solar_k.h"
#include "rttov_calcemis_mw_k.h"
#include "rttov_integrate_k.h"
#include "rttov_intavg_prof_k.h"
#include "rttov_intavg_chan_k.h"
#include "rttov_refsun_k.h"
#include "rttov_init_prof.h"
#include "rttov_init_predictor.h"
#include "rttov_init_raytracing.h"
#include "rttov_check_traj.h"
#include "rttov_setpredictors_9_solar_k.h"
#include "rttov_cldstr_k.h"
#include "rttov_opdpscattir_k.h"
#include "rttov_fresnel_k.h"
#include "rttov_mult_profiles_k.h"
#include "rttov_pcscores.h"
#include "rttov_pcscores_k.h"
#include "rttov_pcscores_rec_k.h"
#include "rttov_reconstruct.h"
#include "rttov_reconstruct_k.h"
#include "rttov_calcbt_pc.h"
#include "rttov_calcbt_pc_ad.h"
#include "rttov_add_raytracing.h"
#include "rttov_add_aux_prof.h"
#include "rttov_add_prof.h"
#include "rttov_init_rad.h"
#include "rttov_init_sunglint.h"
#include "rttov_init_ircld.h"
#include "rttov_init_opdp_path.h"
#include "rttov_add_opdp_path.h"
#include "rttov_init_aux_prof.h"
#include "rttov_init_auxrad_stream.h"
#include "rttov_init_trans_scatt_ir.h"
#include "rttov_init_transmission_aux.h"
#include "rttov_alloc_traj_dyn.h"
#include "rttov_alloc_traj_sta.h"
#include "rttov_direct.h"
!local variables:
  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_k_dyn

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_k
  TYPE(rttov_traj), POINTER :: traj0_k

  TYPE(rttov_options) :: opts_COEF
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nprofiles                     ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                     ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim) :: nchannels_rec
  REAL(KIND=jprb)     , ALLOCATABLE :: clear_k_pc (:, :, :), pcscores_k(:, :, :)
  INTEGER(KIND=jpim)  :: ERR
  REAL   (KIND=JPRB)  :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 0_jpim, ZHOOK_HANDLE)
!-------------
!0. initialize
!-------------
  nprofiles           = size(profiles)
  nchannels           = size(chanprof)
  nlevels             = profiles(1)%nlevels
  opts_COEF           = opts
  opts_COEF%addclouds = .FALSE.
  opts_COEF%addaerosl = .FALSE.
  errorstatus         = errorstatus_success

  IF (opts%addpc) THEN
    npcscores = size(pccomp%pcscores)
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_k)
!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------
! do immediately - for quick escape if target values are out of bounds
  CALL rttov_check_traj( &
        & ERR,               &
        & nprofiles,         &
        & nchannels,         &
        & opts,              &
        & nlevels,           &
        & coefs,             &
        & 1_jpim,            &
        & traj0 = traj0,     &
        & traj0_k = traj0_k, &
        & traj1 = traj1,     &
        & traj1_k = traj1_k, &
        & traj2 = traj,      &
        & traj2_k = traj_k)

  THROWM( ERR .NE. 0 , "rttov_check_traj fatal error")

  CALL rttov_direct( &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & coefs,                     &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_out,            &
            & transmission,              &
            & radiancedata,              &
            & traj         = traj0,      &
            & traj_dyn     = traj0_dyn,  &
            & traj_sta     = traj0_sta,  &
            & pccomp       = pccomp,     &
            & channels_rec = channels_rec)
  
  IF (errorstatus == errorstatus_fatal) THEN
    err = errorstatus_fatal
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_traj_dyn (err, traj0_k_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
  
! K matrix
!----------------
!---------------------------------------------
!0. allocate and initialize local AD variables
!---------------------------------------------
!0.1 auxillary profile K-variables
!----------------------------------
!0.6 path opdep K-variables on COEF levels
!-----------------------------------------
  CALL rttov_init_opdp_path(traj0_k%opdp_path)
  CALL rttov_init_opdp_path(traj0_k%opdp_path_COEF)

  CALL rttov_init_ircld(traj0_k%ircld)
  IF (opts%addaerosl .OR. opts%addclouds) THEN
    CALL rttov_init_trans_scatt_ir(traj0_k%transmission_scatt_ir)
    CALL rttov_init_trans_scatt_ir(traj0_k_dyn%transmission_scatt_ir_stream)
  ENDIF

  CALL rttov_init_auxrad_stream (traj0_k_dyn%auxrad_stream)

  CALL rttov_init_transmission_aux (traj0_k_dyn%transmission_aux)

  CALL rttov_init_prof(traj0_k%profiles_COEF)

! on USER levels
  CALL rttov_init_aux_prof(traj0_k%aux_prof)
! on COEF levels
  CALL rttov_init_aux_prof(traj0_k%aux_prof_COEF)
!0.2 profile_all K-variables on USER levels
!------------------------------------------
!0.2.1 USER LEVELS (all channels)
!--------------------------------

  CALL rttov_init_predictor(traj0_k%predictors)

!0.4 raytracing K-variables
!--------------------------
  CALL rttov_init_raytracing(traj0_k%raytracing)
  CALL rttov_init_raytracing(traj0_k%raytracing_COEF)

  traj0_k%reflectivity(:) = 0._JPRB
  traj0_k%fresnrefl(:)    = 0._JPRB


!0.8 auxrad_stream K-variables on USER levels
!--------------------------------------------

  IF (opts%addpc) THEN
    CALL rttov_init_rad(radiancedata_k)
    radiancedata_k%clear = 1._jprb
  ENDIF


  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_init_sunglint(traj0_k%sunglint)
  ENDIF

!--------------------------------------------------
!1. K of radiative transfer integration - USER levs
!--------------------------------------------------
! use the subroutine argument for K radiances
  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)
  CALL rttov_integrate_k( &
        & addcosmic,                                &
        & opts,                                     &
        & traj0_dyn%nstreams,                       &
        & chanprof,                                 &
        & emissivity_out,                           &
        & emissivity_out_k,                         &
        & traj0%reflectivity,                       &
        & traj0_k%reflectivity,                     &
        & traj0%fresnrefl,                          &
        & traj0_k%fresnrefl,                        &
        & traj0%sunglint,                           &
        & traj0_k%sunglint,                         &
        & traj0_sta%sun,                            &
        & traj0_dyn%transmission_aux,               &
        & traj0_k_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream,   &
        & traj0_k_dyn%transmission_scatt_ir_stream, &
        & profiles,                                 &
        & profiles_k,                               &
        & traj0%aux_prof,                           &
        & traj0_k%aux_prof,                         &
        & coefs%coef,                               &
        & traj0%raytracing,                         &
        & traj0_k%raytracing,                       &
        & traj0%ircld,                              &
        & traj0_k%ircld,                            &
        & radiancedata,                             &
        & traj0_sta%auxrad,                         &
        & traj0_dyn%auxrad_stream,                  &
        & traj0_k_dyn%auxrad_stream,                &
        & radiancedata_k)
!--------------------------------------
!2. K of channel emissivities - SURFACE
!--------------------------------------

  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_fresnel_k( &
          & chanprof,          &
          & profiles,          &
          & coefs%coef,        &
          & traj0%sunglint,    &
          & traj0_k%sunglint,  &
          & traj0%fresnrefl,   &
          & traj0_k%fresnrefl)
    CALL rttov_refsun_k( &
          & chanprof,           &
          & profiles,           &
          & profiles_k,         &
          & coefs%coef,         &
          & traj0%aux_prof,     &
          & traj0%sunglint,     &
          & traj0_k%sunglint,   &
          & traj0%raytracing,   &
          & traj0_k%raytracing)
  ENDIF


  IF (Any(calcemis)) THEN
! calculate surface emissivity for selected channels
! and traj0%reflectivity

    IF (coefs%coef%id_sensor == sensor_id_ir) THEN
!Infrared
      emissivity_out_k(:) =  - traj0_k%reflectivity(:) + emissivity_out_k(:)
    ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
!Microwave
      CALL rttov_calcemis_mw_k( &
            & opts,                           &
            & profiles,                       &
            & profiles_k,                     &
            & traj0_sta%angles,               &
            & coefs%coef,                     &
            & chanprof,                       &
            & traj0_dyn%transmission_aux,     &
            & traj0_k_dyn%transmission_aux,   &
            & calcemis,                       &
            & emissivity_out_k,               &
            & traj0_k%reflectivity)
    ELSE
! Hires
      emissivity_out_k(:) =  - traj0_k%reflectivity(:) + emissivity_out_k(:)
      CALL rttov_calcemis_ir_k( &
            & profiles,          &
            & profiles_k,        &
            & coefs%coef,        &
            & opts%addpc,        &
            & coefs%coef_pccomp, &
            & chanprof,          &
            & calcemis,          &
            & emissivity_out_k)
    ENDIF

  ENDIF


  WHERE (.NOT. calcemis)
    emissivity_out_k = emissivity_out_k - traj0_k%reflectivity
  ENDWHERE

! Where( .not. calcemis )
!   emissivity_k = emissivity_k + emissivity_out_k
! End Where
! Where( .not. calcemis )
!   emissivity_out_k = 0._jprb
! End Where
  emissivity_k = emissivity_k + emissivity_out_k
!----------------------------------
!3. K of transmittances - USER levs
!----------------------------------

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_transmit_9_solar_k( &
            & opts%addaerosl,                           &
            & opts%addclouds,                           &
            & profiles(1)%nlayers,                      &
            & chanprof,                                 &
            & profiles,                                 &
            & traj0_sta%sun,                            &
            & traj0%aux_prof,                           &
            & traj0_k%aux_prof,                         &
            & coefs%coef,                               &
            & traj0%raytracing,                         &
            & traj0_k%raytracing,                       &
            & traj0%ircld,                              &
            & traj0%opdp_path,                          &
            & traj0_k%opdp_path,                        &
            & traj0_sta%odsun_level,                    &
            & traj0_sta%odsun_singlelayer,              &
            & traj0_sta%od_frac,                        &
            & traj0_dyn%transmission_aux,               &
            & traj0_k_dyn%transmission_aux,             &
            & traj0_dyn%transmission_scatt_ir_stream,   &
            & traj0_k_dyn%transmission_scatt_ir_stream, &
            & traj0_sta%tausun_ref,                     &
            & traj0_sta%tausun_ref_surf,                &
            & traj0_sta%tausun_level,                   &
            & traj0_sta%tausun_surf)
    ENDIF
  ENDIF

  CALL rttov_transmit_k( &
        & opts%addaerosl,                           &
        & opts%addclouds,                           &
        & profiles(1)%nlayers,                      &
        & chanprof,                                 &
        & traj0%aux_prof,                           &
        & traj0_k%aux_prof,                         &
        & coefs%coef,                               &
        & traj0%ircld,                              &
        & traj0%opdp_path,                          &
        & traj0_k%opdp_path,                        &
        & traj0_sta%od_level,                       &
        & transmission,                             &
        & transmission_k,                           &
        & traj0_dyn%transmission_aux,               &
        & traj0_k_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream,   &
        & traj0_k_dyn%transmission_scatt_ir_stream, &
        & traj0_sta%tau_ref,                        &
        & traj0_sta%tau_ref_surf,                   &
        & traj0_sta%tau_surf,                       &
        & traj0_sta%tau_level)
!------------------------------------------------------------
!4. K of optical depths of aerosols and/or clouds - USER levs
!------------------------------------------------------------

  IF (opts%addaerosl .OR. opts%addclouds) THEN
    CALL rttov_opdpscattir_k( &
          & profiles(1)%nlayers,                       &
          & chanprof,                                  &
          & opts,                                      &
          & traj0%aux_prof,                            &
          & traj0_k%aux_prof,                          &
          & profiles,                                  &
          & profiles_k,                                &
          & traj0_sta%sun,                             &
          & coefs%coef,                                &
          & coefs%coef_scatt_ir,                       &
          & traj0%raytracing,                          &
          & traj0_k%raytracing,                        &
          & traj0%transmission_scatt_ir,               &
          & traj0_k%transmission_scatt_ir,             &
          & traj0_dyn%transmission_scatt_ir_stream,    &
          & traj0_k_dyn%transmission_scatt_ir_stream,  &
          & coefs%optp,                                &
          & traj0%ircld,                               &
          & traj0_k%ircld)
  ENDIF

!---------------------------------------------------
!5. K of cloud streams and distributions - USER levs
!---------------------------------------------------

  IF (opts%addclouds) THEN
    CALL rttov_cldstr_k( &
          & chanprof,       &
          & profiles,       &
          & profiles_k,     &
          & traj0%ircld,    &
          & traj0_k%ircld)
  ENDIF

!-----------------------------------
!6. K of optical depth interpolation
!-----------------------------------

  if(opts%SpaceTop) traj0_k%opdp_path%atm_level(1,:) = 0._jprb
  if(opts%SpaceTop) traj0_k%opdp_path%sun_level(1,:) = 0._jprb

  IF (opts%addinterp) THEN
    CALL rttov_intavg_chan_k( &
          & opts%lgradp,                    &
          & traj0_sta%sun,                  &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_COEF,            &
          & profiles,                       &
          & profiles_k,                     &
          & traj0%opdp_path_COEF,           &
          & traj0_k%opdp_path_COEF,         &
          & traj0%opdp_path,                &
          & traj0_k%opdp_path)
  ELSE
    CALL rttov_add_opdp_path(traj0_k%opdp_path_COEF, traj0_k%opdp_path_COEF, traj0_k%opdp_path)
  ENDIF

!--------------------------------------------------------
!7. K of atmospheric and solar optical depths - COEF levs
!--------------------------------------------------------
! optical depth arrays allocated at start of K-code

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_opdep_9_solar_k( &
            & coefs%coef%nlayers,     &
            & chanprof,               &
            & profiles,               &
            & traj0_sta%sun,          &
            & traj0%predictors,       &
            & traj0_k%predictors,     &
            & coefs%coef,             &
            & traj0%opdp_path_COEF,   &
            & traj0_k%opdp_path_COEF, &
            & traj0_sta%opdpsun_ref_COEF)
    ENDIF
    CALL rttov_opdep_9_k( &
          & coefs%coef%nlayers,     &
          & chanprof,               &
          & traj0%predictors,       &
          & traj0_k%predictors,     &
          & traj0%aux_prof_COEF,    &
          & traj0_k%aux_prof_COEF,  &
          & coefs%coef,             &
          & traj0%opdp_path_COEF,   &
          & traj0_k%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ELSE
    CALL rttov_opdep_k( &
          & coefs%coef%nlayers,     &
          & chanprof,               &
          & traj0%predictors,       &
          & traj0_k%predictors,     &
          & traj0%aux_prof_COEF,    &
          & traj0_k%aux_prof_COEF,  &
          & coefs%coef,             &
          & traj0%opdp_path_COEF,   &
          & traj0_k%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ENDIF

!-------------------------------------------------------
!8. K of RTTOV-7 RTTOV-8 RTTOV-9 predictors - COEF levs
!-------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_k( &
          & opts,                    &
          & coefs%coef%nlayers,      &
          & traj0_sta%angles_COEF,   &
          & chanprof,                &
          & traj0%profiles_COEF,     &
          & traj0_k%profiles_COEF,   &
          & coefs%coef,              &
          & traj0%aux_prof_COEF,     &
          & traj0%predictors,        &
          & traj0_k%predictors,      &
          & traj0%raytracing_COEF,   &
          & traj0_k%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_k( &
          & opts,                    &
          & coefs%coef%nlayers,      &
          & traj0_sta%angles_COEF,   &
          & chanprof,                &
          & traj0%profiles_COEF,     &
          & traj0_k%profiles_COEF,   &
          & coefs%coef,              &
          & traj0%aux_prof_COEF,     &
          & traj0%predictors,        &
          & traj0_k%predictors,      &
          & traj0%raytracing_COEF,   &
          & traj0_k%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_setpredictors_9_solar_k( &
            & opts,                    &
            & chanprof,                &
            & traj0%profiles_COEF,     &
            & traj0_k%profiles_COEF,   &
            & coefs%coef,              &
            & traj0%predictors,        &
            & traj0_k%predictors,      &
            & traj0%raytracing_COEF,   &
            & traj0_k%raytracing_COEF)
    ENDIF
    CALL rttov_setpredictors_9_k( &
          & opts,                    &
          & coefs%coef%nlayers,      &
          & traj0_sta%angles_COEF,   &
          & coefs%coef_pccomp,       &
          & chanprof,                &
          & traj0%profiles_COEF,     &
          & traj0_k%profiles_COEF,   &
          & coefs%coef,              &
          & traj0%predictors,        &
          & traj0_k%predictors,      &
          & traj0%raytracing_COEF,   &
          & traj0_k%raytracing_COEF)
  ENDIF

!------------------------------------------------------
!9. K of common geometric set-up for the RT integration
!------------------------------------------------------
!9.1 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_setgeometry_k( &
          & opts,                    &
          & chanprof,                &
          & traj0%profiles_COEF,     &
          & traj0_k%profiles_COEF,   &
          & traj0%aux_prof_COEF,     &
          & coefs%coef,              &
          & traj0_sta%angles_COEF,   &
          & traj0%raytracing_COEF,   &
          & traj0_k%raytracing_COEF)
  ELSE
    CALL rttov_add_raytracing(traj0_k%raytracing, traj0_k%raytracing, traj0_k%raytracing_COEF)
  ENDIF

!9.1 USER levs
!-------------
  CALL rttov_setgeometry_k( &
        & opts,               &
        & chanprof,           &
        & profiles,           &
        & profiles_k,         &
        & traj0%aux_prof,     &
        & coefs%coef,         &
        & traj0_sta%angles,   &
        & traj0%raytracing,   &
        & traj0_k%raytracing)
!---------------------------------------------------------
!10. AD of cloud top, surface levels, ice cloud parameters
!---------------------------------------------------------
!10.1 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_profaux_k( &
          & opts_COEF,             &
          & chanprof,              &
          & traj0%profiles_COEF,   &
          & traj0_k%profiles_COEF, &
          & coefs%coef,            &
          & traj0%aux_prof_COEF,   &
          & traj0_k%aux_prof_COEF)
  ELSE
    CALL rttov_add_aux_prof(traj0_k%aux_prof, traj0_k%aux_prof, traj0_k%aux_prof_COEF)
  ENDIF

!10.2 USER levs
!-------------
  CALL rttov_profaux_k( &
        & opts,             &
        & chanprof,         &
        & profiles,         &
        & profiles_k,       &
        & coefs%coef,       &
        & traj0%aux_prof,   &
        & traj0_k%aux_prof)
!------------------------------------------------------------------
!11 K of check input data is within suitable physical limits - COEF levs
!------------------------------------------------------------------

  IF (opts%apply_reg_limits) THEN

    CALL rttov_checkinput_k( &
          & opts,                         &
          & chanprof,                     &
          & traj0_sta%profiles_COEF_ref,  &
          & traj0_k%profiles_COEF,        &
          & coefs%coef)

  ENDIF

!----------------------------------------------------------------
!12. K of profile variable interpolation - COEF levs to USER levs
!----------------------------------------------------------------
! intvar logical array set by direct code
! loop over chanels inside subroutine

  IF (opts%addinterp) THEN
    CALL rttov_intavg_prof_k( &
          & opts,                           &
          & chanprof,                       &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_k,                     &
          & traj0%profiles_COEF,            &
          & traj0_k%profiles_COEF)!inout  target variables
  ELSE
    CALL rttov_add_prof( &
          & profiles_k,            &
          & profiles_k,            &
          & traj0_k%profiles_COEF, &
          & lair = .TRUE._jplm,    &
          & lground = .FALSE._jplm)
  ENDIF

  CALL rttov_add_prof( &
        & profiles_k,            &
        & profiles_k,            &
        & traj0_k%profiles_COEF, &
        & lair = .FALSE._jplm,   &
        & lground = .TRUE._jplm)

  IF (opts%addpc) THEN

    IF (opts%addradrec) THEN
      nchannels_rec = size(traj0_sta%chanprof_in)
      npcscores     = size(traj0_sta%chanprof_pc)
      ALLOCATE (pcscores_k(nchannels_rec / nprofiles, npcscores / nprofiles, nprofiles),      &
        & clear_k_pc(nchannels_rec / nprofiles, nchannels / nprofiles, nprofiles), STAT = ERR)

      THROWM( ERR .NE. 0 , "allocation pc k arrays")


      IF (opts%switchrad) THEN
        pccomp_k%bt_pccomp    = 1._jprb
        clear_k_pc            = 0._jprb
        pcscores_k            = 0._jprb
        pccomp_k%clear_pccomp = 0._jprb
        CALL rttov_calcbt_pc_ad( &
              & traj0_sta%chanprof_in,       &
              & coefs%coef_pccomp,           &
              & pccomp,                      &
              & pccomp_k)
      ELSE
        pccomp_k%bt_pccomp    = 0._jprb
        clear_k_pc            = 0._jprb
        pccomp_k%clear_pccomp = 1._jprb
        pcscores_k            = 0._jprb
      ENDIF

      CALL rttov_reconstruct_k( &
            & traj0_sta%chanprof_in,       &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & pccomp_k,                    &
            & pcscores_k,                  &
            & coefs%coef_pccomp)
      CALL rttov_pcscores_rec_k(           &
            & opts,                        &
            & chanprof,                    &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & pcscores_k,                  &
            & coefs%coef_pccomp,           &
            & clear_k_pc)
      CALL rttov_init_prof(profiles_k_rec)
! Pas de gradient d'emissivite pour les PC (calcemis==false)
      CALL rttov_mult_profiles_k(profiles_k_rec, profiles_k, clear_k_pc)
      DEALLOCATE (pcscores_k, clear_k_pc, STAT = ERR)

      THROWM( ERR .NE. 0 , "deallocation pc k arrays")

    ELSE
      npcscores = size(traj0_sta%chanprof_pc)
      ALLOCATE (clear_k_pc(npcscores / nprofiles, nchannels / nprofiles, nprofiles), STAT = ERR)

      THROWM( ERR .NE. 0 , "allocation pc k arrays")

      pccomp_k%pcscores = 1._jprb
      clear_k_pc        = 0._jprb
      CALL rttov_pcscores_k( &
            & opts,                        &
            & chanprof,                    &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & pccomp_k,                    &
            & coefs%coef_pccomp,           &
            & clear_k_pc)
      CALL rttov_init_prof(profiles_k_pc)
      CALL rttov_mult_profiles_k(profiles_k_pc, profiles_k, clear_k_pc)
      DEALLOCATE (clear_k_pc, STAT = ERR)

      THROWM( ERR .NE. 0 , "deallocation pc k arrays")

    ENDIF

  ENDIF

!
!-----------------------------------------
!13. deallocate memory for local variables
!-----------------------------------------

  CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
  THROW(err.ne.0)
!
  CALL rttov_alloc_traj_dyn (err, traj0_k_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
  THROW(err.ne.0)
!
  CALL rttov_alloc_traj_sta (err, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, 0_jpim, npcscores, channels_rec)
  THROW(err.ne.0)
!


  CALL rttov_check_traj( &
        & ERR,               &
        & nprofiles,         &
        & nchannels,         &
        & opts,              &
        & nlevels,           &
        & coefs,             &
        & 0_jpim,            &
        & traj0 = traj0,     &
        & traj0_k = traj0_k, &
        & traj1 = traj1,     &
        & traj1_k = traj1_k, &
        & traj2 = traj,      &
        & traj2_k = traj_k)

  THROWM( ERR .NE. 0 , "deallocation check_traj")

  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)
  CATCH
  errorstatus = ERR
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_k
