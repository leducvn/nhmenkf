!
SUBROUTINE rttov_direct( &
            & errorstatus,    &
            & chanprof,       &
            & opts,           &
            & profiles,       &
            & coefs,          &
            & calcemis,       &
            & emissivity,     &
            & emissivity_out, &
            & transmission,   &
            & radiancedata,   &
            & traj,           &
            & traj_dyn,       &
            & traj_sta,       &
            & pccomp,         &
            & channels_rec,   &
            & lbl_check)
!
! Description:
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
! Rochon Y.J., L. Garand, D.S. Turner and S. Polavarapu Jacobian mapping
! between  vertical coordinate systems in data assimilation (submitted QJRMS
! June 2006)
!
! Matricardi, M. 2009: An Observation operator for the assimilation of
! principal component scores into a NWP system. Available from EUMETSAT
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!          13/8/92.  For version 2.
!                    ksat added to argument list; ssu included;
!                    internal changes to move big arrays from commons to
!                    arguments and to introduce taskcommons
!           8/7/97   added ozone and extended water vapour in control vector
!        01/05/2000  F90 code
!        21/08/2000  Interface to rtint changed to include pref (surface traj0%reflectivity).
!                    (Stephen English)
!        31/01/2001  More cloud computations. stored in radov (F. Chevallier)
!        6/2/2001    pgrody and knav etc arrays removed from call (R Saunders)
!        18/01/2002  Thread safe (D.Salmond)
!        01/12/2002  New F90 code with structures (P Brunel A Smith)
!        02/01/2003  More comments added (R Saunders)
!        24/01/2003  Error return code by input profile (P Brunel)
!                    Add WV Continuum and CO2 capability
!        02/06/2004  Change tests on id_comp_lvl == 7 by tests on fmv_model_ver (P. Brunel)
!        08/12/2005  Added surface humidity to lowest mean layer q (R Saunders)
!        16/01/2006  Marco Matricardi (ECMWF):
!              --    IASI capability added.
!              --    Variable trace gases added for IASI and AIRS.
!              --    Solar radiation added for IASI and AIRS
!              --    Altitude dependent local zenith angle added.
!              --    Linear in tau approximation for RT equation added.
!         1/07/2006  Marco Matricardi (ECMWF):
!              --    Parameterization of multiple scattering for aerosols
!                    and clouds added
!        22/01/2007  Removed polarisation indexing R. Saunders.
!        29/06/2007  Introduce interpolator (P.J.Rayer)
!              --    Profiles on USER levs -> COEF levs
!              --    Duplicate profaux and setgeometry (USER and COEF levs)
!              --    Predict atm/sol optical depths on COEF levs -> USER levs
!              --    Calculate transmittance and integrate RT on USER levs
!        22/08/2007  Remove Allocates (D Salmond)
!        11/10/2007  Move iaernum and iaertyp to profile_aux P.Marguinaud
!        26/09/2008  Nullify temp pointers; this is for the NEC
!                    under -Chopt. P. Marguinaud
!        03/10/2008  Fix bug in checking addsolar
!        20/01/2009  Deallocate profile memory if profile outside limits R Saunders
!        27/02/2009  Profile levels to include ToA. Distinguish between layer
!                    arrays  and level arrays - size, index labels, looping
!                    (P. Rayer)
!        03/11/2009  Transmittances on levels (A Geer)
!        02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
!        17/06/2010  Introduced optional spacetop flag to zero opdeps at user's
!                    model-top in Zeeman channels (P Rayer)
!        05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!        14/10/2010  Remove rt8_mode (J Hocking)
!        18/10/2010  Add control on regression channels for PC calcs (P Brunel)
!        14/12/2010  Use traj0_sta%sun array to flag channels for which solar calculations
!                    should be performed (J Hocking)
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
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_Type,      &
       & transmission_type, &
       & radiance_Type,     &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_traj,        &
       & rttov_traj_dyn,    &
       & rttov_traj_sta,    &
       & rttov_lbl_check
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw, &
       & sensor_id_ir, &
       & sensor_id_hi, &
       & sensor_id_po, &
       & ncldtyp,      &
       & max_sol_zen
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)                   ! Atmospheric profiles
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)                   ! Channel indices (nchannels)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
! (supplied on User levels) (nprofiles)
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs                         ! It is necessary to have "Target" attribute here
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis      (size(chanprof))! switch for emmissivity calc.
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity    (size(chanprof))! surface emmissivity
  REAL(KIND=jprb)        , INTENT(INOUT)                        :: emissivity_out(size(chanprof))! surface emmissivity
  TYPE(transmission_type), INTENT(INOUT)                        :: transmission                  ! transmittances and singlelayer optical depths (on User levels)
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                  ! radiances (mw/cm-1/ster/sq.m) and degK
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus                   ! return flag
  TYPE(rttov_traj  )     , INTENT(INOUT), OPTIONAL     , TARGET :: traj                          ! Target is *NEEDED* here (see rttov_check_temp)
  TYPE(rttov_traj_dyn)   , INTENT(INOUT), OPTIONAL     , TARGET :: traj_dyn                      
  TYPE(rttov_traj_sta)   , INTENT(INOUT), OPTIONAL     , TARGET :: traj_sta
  TYPE(rttov_pccomp)     , OPTIONAL     , INTENT(INOUT)         :: pccomp
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
  TYPE(rttov_lbl_check)  , OPTIONAL     , INTENT(IN)            :: lbl_check
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_cldstr.h"
#include "rttov_checkinput.h"
#include "rttov_profaux.h"
#include "rttov_setgeometry.h"
#include "rttov_setpredictors_7.h"
#include "rttov_setpredictors_8.h"
#include "rttov_setpredictors_9.h"
#include "rttov_setpredictors_9_solar.h"
#include "rttov_opdpscattir.h"
#include "rttov_fresnel.h"
#include "rttov_opdep.h"
#include "rttov_opdep_9.h"
#include "rttov_opdep_9_solar.h"
#include "rttov_transmit.h"
#include "rttov_transmit_9_solar.h"
#include "rttov_calcemis_ir.h"
#include "rttov_calcemis_mw.h"
#include "rttov_integrate.h"
#include "rttov_intavg_chan.h"
#include "rttov_intavg_prof.h"
#include "rttov_refsun.h"
#include "rttov_check_traj.h"
#include "rttov_init_prof.h"
#include "rttov_init_raytracing.h"
#include "rttov_reconstruct.h"
#include "rttov_pcscores.h"
#include "rttov_calcbt_pc.h"
#include "rttov_copy_raytracing.h"
#include "rttov_copy_prof.h"
#include "rttov_copy_aux_prof.h"
#include "rttov_copy_opdp_path.h"
#include "rttov_init_ircld.h"
#include "rttov_init_aux_prof.h"
#include "rttov_alloc_traj_dyn.h"
#include "rttov_alloc_traj_sta.h"
#include "rttov_checkpcchan.h"

!local variables:
  INTEGER(KIND=jpim)               :: i, j, prof, chan                                      ! loop index
  INTEGER(KIND=jpim)               :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_options)              :: opts_COEF
!
  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj_dyn), TARGET  :: traj1_dyn
  TYPE(rttov_traj_dyn), POINTER :: traj0_dyn
  TYPE(rttov_traj_sta), TARGET  :: traj1_sta
  TYPE(rttov_traj_sta), POINTER :: traj0_sta
  LOGICAL(KIND=jplm) :: ltraj_dyn_dealloc
  LOGICAL(KIND=jplm) :: ltraj_sta_dealloc
  INTEGER(KIND=jpim) :: nprofiles                    ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                    ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim)  :: npcscores
  INTEGER(KIND=jpim)  :: ERR
  INTEGER(KIND=JPIM)  :: dc(size(chanprof)/size(profiles))

  REAL(KIND=JPRB)     :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

  TRY
!-------------
!0. initialize
!-------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 0_jpim, ZHOOK_HANDLE)
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

  IF(opts%addpc) THEN
    Call rttov_checkpcchan( & 
       & nprofiles,         &
       & nchannels,         &
       & opts,              &
       & chanprof,          &
       & coefs,             &
       & ERR                )
    THROWM( ERR .NE. 0 , "rttov_checkpcchan fatal error")
  ENDIF

  NULLIFY (traj0, traj0_dyn, traj0_sta)
!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------
! do immediately - for quick escape if target values are out of bounds
  CALL rttov_check_traj( &
        & ERR,           &
        & nprofiles,     &
        & nchannels,     &
        & opts,          &
        & nlevels,       &
        & coefs,         &
        & 1_jpim,        &
        & traj0 = traj0, &
        & traj1 = traj1, &
        & traj2 = traj)
  THROWM( ERR .NE. 0 , "rttov_check_traj fatal error")

!
!
  IF (Present (traj_dyn)) THEN
    traj0_dyn => traj_dyn
  ELSE
    traj0_dyn => traj1_dyn
  ENDIF
  ltraj_dyn_dealloc = .not. Present(traj_dyn)

!
!
  IF (Present (traj_sta)) THEN
    traj0_sta => traj_sta
  ELSE
    traj0_sta => traj1_sta
  ENDIF
  ltraj_sta_dealloc = .not. Present(traj_sta)


!

  CALL rttov_alloc_traj_sta (err, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                            & 1_jpim, npcscores, channels_rec)
  THROW(err.ne.0)


  CALL rttov_init_prof(traj0%profiles_COEF, p = coefs%coef%ref_prfl_p)

  IF (opts%addinterp) THEN
! set up profile variables for interpolation
! easy extension to more trace following gas id codes in rttov_const
!    nvar=ngases_max
! add more trace gases here if appropriate
    CALL rttov_intavg_prof( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & traj0%profiles_COEF)!inout  target variables
!
  ELSE
    CALL rttov_copy_prof( &
          & traj0%profiles_COEF,  &
          & profiles,             &
          & larray = .TRUE._jplm, &
          & lscalar = .FALSE._jplm)
  ENDIF

  CALL rttov_copy_prof( &
        & traj0%profiles_COEF,   &
        & profiles,              &
        & larray = .FALSE._jplm, &
        & lscalar = .TRUE._jplm)
  CALL rttov_init_raytracing(traj0%raytracing)
!------------------------------------------------------------------
!2. check input data is within suitable physical limits - COEF levs
!------------------------------------------------------------------


  IF (opts%do_checkinput) THEN

    IF (opts%apply_reg_limits) THEN
      CALL rttov_copy_prof (traj0_sta%profiles_coef_ref, traj0%profiles_coef)
    ENDIF

     CALL rttov_checkinput( &
           & opts,                   &
           & traj0%profiles_COEF,    &
           & coefs%coef,             &
           & coefs%coef_pccomp,      &
           & err)
    IF (err.eq.errorstatus_fatal) THEN
!  nothing processed so all profiles get the fatal error code
!  user will know which profile
!  Deallocate profile memory
      CALL rttov_check_traj( &
            & ERR,           &
            & nprofiles,     &
            & nchannels,     &
            & opts,          &
            & nlevels,       &
            & coefs,         &
            & 0_jpim,        &
            & traj0 = traj0, &
            & traj1 = traj1, &
            & traj2 = traj)
  
      THROWM( ERR .NE. 0 , "deallocation check_traj")
  
      err = errorstatus_fatal
      THROW(err.ne.0)
  
  
    ENDIF
  ENDIF

! Set flags to indicate if solar calculations should be performed for each channel.
  IF (ASSOCIATED(coefs%coef%ss_val_chn)) THEN
    DO j = 1, nchannels
      prof = chanprof(j)%prof
      chan = chanprof(j)%chan
    
      IF (opts%addsolar .AND. profiles(prof)%sunzenangle >= 0.0 .AND. &
                              profiles(prof)%sunzenangle < max_sol_zen) THEN  
        IF (coefs%coef%ss_val_chn(chan) .NE. 0_jpim) THEN
          traj0_sta%sun(j) = .TRUE.
        ELSE
          traj0_sta%sun(j) = .FALSE.
        ENDIF
      ELSE
        traj0_sta%sun(j) = .FALSE.
      END IF
    ENDDO
  ELSE
    traj0_sta%sun(:) = .FALSE.
  ENDIF



!------------------------------------------------------------------------
!3. determine cloud top, surface levels, ice cloud parameters
!------------------------------------------------------------------------
!3.1 USER levs
!-------------
  CALL rttov_profaux( &
        & opts,           &
        & profiles,       &
        & coefs%coef,     &
        & traj0%aux_prof)
!3.2 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_profaux( &
          & opts_COEF,           &
          & traj0%profiles_COEF, &
          & coefs%coef,          &
          & traj0%aux_prof_COEF)
  ELSE
    CALL rttov_copy_aux_prof(traj0%aux_prof_COEF, traj0%aux_prof)
  ENDIF

!------------------------------------------------------------------------
!4. set up common geometric variables for the RT integration
!------------------------------------------------------------------------
!4.1 USER levs
!-------------
  CALL rttov_setgeometry( &
        & opts,             &
        & profiles,         &
        & traj0%aux_prof,   &
        & coefs%coef,       &
        & traj0_sta%angles, &
        & traj0%raytracing)
!4.2 COEF levs
!-------------


  IF (Present (lbl_check)) THEN
!
! For lbl checks, we may want to test in a plane geometry
!
    IF( lbl_check%plane_geometry ) THEN
      DO i = 1, nprofiles
        traj0%raytracing%PATHSAT(:,i) = traj0_sta%angles(i)%seczen
      ENDDO
    ENDIF
  ENDIF


  IF (opts%addinterp) THEN
    CALL rttov_setgeometry( &
          & opts,                  &
          & traj0%profiles_COEF,   &
          & traj0%aux_prof_COEF,   &
          & coefs%coef,            &
          & traj0_sta%angles_COEF, &
          & traj0%raytracing_COEF)
  ELSE
    traj0_sta%angles_COEF = traj0_sta%angles
    CALL rttov_copy_raytracing(traj0%raytracing_COEF, traj0%raytracing)
  ENDIF

!---------------------------------------------------
!5. calculate atm/sol predictors - atm/sol COEF levs
!---------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7( &
          & opts,                  &
          & traj0%profiles_COEF,   &
          & traj0_sta%angles_COEF, &
          & coefs%coef,            &
          & traj0%aux_prof_COEF,   &
          & traj0%predictors,      &
          & traj0%raytracing_COEF)! inout  (in because of mem allocation)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8( &
          & opts,                  &
          & traj0%profiles_COEF,   &
          & traj0_sta%angles_COEF, &
          & coefs%coef,            &
          & traj0%aux_prof_COEF,   &
          & traj0%predictors,      &
          & traj0%raytracing_COEF)! inout  (in because of mem allocation)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    CALL rttov_setpredictors_9( &
          & opts,                  &
          & traj0%profiles_COEF,   &
          & traj0_sta%angles_COEF, &
          & coefs%coef_pccomp,     &
          & coefs%coef,            &
          & traj0%predictors,      &
          & traj0%raytracing_COEF)! inout  (in because of mem allocation)
    IF (opts%addsolar) THEN
      CALL rttov_setpredictors_9_solar( &
            & opts,                  &
            & traj0%profiles_COEF,   &
            & coefs%coef,            &
            & traj0%predictors,      &
            & traj0%raytracing_COEF)
    ENDIF
  ELSE
    err = errorstatus_fatal

    THROWM( ERR .NE. 0 , "Unexpected RTTOV compatibility version number")

  ENDIF

!  End Do ! Profile loop
!--------------------------------------------------------------------
!6. predict atmospheric and solar optical depths - atm/sol COEF levs
!--------------------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 9) THEN
    CALL rttov_opdep_9( &
          & coefs%coef%nlayers,   &
          & chanprof,             &
          & traj0%predictors,     &
          & traj0%aux_prof_COEF,  &
          & coefs%coef,           &
          & traj0%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
    IF (opts%addsolar) THEN
      CALL rttov_opdep_9_solar( &
            & coefs%coef%nlayers,   &
            & chanprof,             &
            & profiles,             &
            & traj0_sta%sun,        &
            & traj0%predictors,     &
            & coefs%coef,           &
            & traj0%opdp_path_COEF, &
            & traj0_sta%opdpsun_ref_COEF)
    ENDIF
  ELSE
    CALL rttov_opdep( &
          & coefs%coef%nlayers,   &
          & chanprof,             &
          & traj0%predictors,     &
          & traj0%aux_prof_COEF,  &
          & coefs%coef,           &
          & traj0%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ENDIF

!--------------------------------------------------------------------------
!7. interpolator second  call - optical depths from COEF levs to USER levs
!--------------------------------------------------------------------------

  IF (opts%addinterp) THEN
    CALL rttov_intavg_chan( &
          & opts%lgradp,                    &
          & traj0_sta%sun,                  &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_COEF,            &
          & profiles,                       &
          & traj0%opdp_path_COEF,           &
          & traj0%opdp_path)
  ELSE
    CALL rttov_copy_opdp_path(traj0%opdp_path, traj0%opdp_path_COEF)
  ENDIF
  if(opts%SpaceTop) traj0%opdp_path%atm_level(1,:) = 0._jprb
  if(opts%SpaceTop) traj0%opdp_path%sun_level(1,:) = 0._jprb
 
   !
   ! This is for coefficient testing; we replace the optical depths
   ! with those from the line-by-line convoluated with ISRF
   !

  IF (Present (lbl_check)) THEN
    IF (Associated (lbl_check%atm_layer)) THEN
      traj0%opdp_path%atm_level(1,:) = 0._jprb
      traj0%opdp_path%atm_level(2:,:) = lbl_check%atm_layer(:,:)
    ENDIF
  ENDIF


!--------------------------------------------------------------------
!8. If clouds are present,calculate the number of streams and
!   the cloud distribution  in each stream - remains on USER levs
!--------------------------------------------------------------------
  traj0_dyn%nstreams = 0_jpim

  IF (opts%addclouds) THEN
    CALL rttov_cldstr( &
          & profiles,              &
          & opts%cldstr_threshold, &
          & traj0%ircld,           &
          & traj0_dyn%nstreams)
  ELSE
    traj0%ircld%XSTRCLR = 1._jprb
    traj0%ircld%NSTREAM = 0_jpim
  ENDIF

!----------------------------------------------------------------------------
!9. Calculate optical depths of aerosols and/or clouds - USER levs
!----------------------------------------------------------------------------

  CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(ERR .NE. 0)

  IF (opts%addaerosl .OR. opts%addclouds) THEN
    CALL rttov_opdpscattir( &
          & profiles(1)%nlayers,                    &
          & chanprof,                               &
          & opts,                                   &
          & traj0%aux_prof,                         &
          & profiles,                               &
          & traj0_sta%sun,                          &
          & coefs%coef,                             &
          & coefs%coef_scatt_ir,                    &
          & traj0%raytracing,                       &
          & traj0%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream, &
          & coefs%optp,                             &
          & traj0%ircld)
  ENDIF

!----------------------------------------
!10. calculate transmittances - USER levs
!----------------------------------------
  CALL rttov_transmit( &
        & opts%addaerosl,                         &
        & opts%addclouds,                         &
        & profiles(1)%nlayers,                    &
        & chanprof,                               &
        & traj0%aux_prof,                         &
        & coefs%coef,                             &
        & traj0%ircld,                            &
        & traj0%opdp_path,                        &
        & traj0_sta%od_level,                     &
        & transmission,                           &
        & traj0_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream, &
        & traj0_sta%tau_ref,                      &
        & traj0_sta%tau_ref_surf,                 &
        & traj0_sta%tau_surf,                     &
        & traj0_sta%tau_level)

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_transmit_9_solar( &
            & opts%addaerosl,                         &
            & opts%addclouds,                         &
            & profiles(1)%nlayers,                    &
            & chanprof,                               &
            & profiles,                               &
            & traj0_sta%sun,                          &
            & traj0%aux_prof,                         &
            & coefs%coef,                             &
            & traj0%raytracing,                       &
            & traj0%ircld,                            &
            & traj0%opdp_path,                        &
            & traj0_sta%odsun_level,                  &
            & traj0_sta%odsun_singlelayer,            &
            & traj0_sta%od_frac,                      &
            & traj0_dyn%transmission_aux,             &
            & traj0_dyn%transmission_scatt_ir_stream, &
            & traj0_sta%tausun_ref,                   &
            & traj0_sta%tausun_ref_surf,              &
            & traj0_sta%tausun_surf,                  &
            & traj0_sta%tausun_level)
    ENDIF
  ENDIF

!-------------------------------------------------
!11. calculate channel emissivity values - SURFACE
!-------------------------------------------------

  WHERE (.NOT. calcemis)
    emissivity_out = emissivity
  ELSEWHERE
    emissivity_out = 0._jprb
  ENDWHERE


  WHERE (.NOT. calcemis)
    traj0%reflectivity = 1 - emissivity_out
  ENDWHERE


  IF (Any(calcemis)) THEN
! calculate surface emissivity for selected channels
! and traj0%reflectivity

    IF (coefs%coef%id_sensor == sensor_id_ir) THEN
!Infrared
      CALL rttov_calcemis_ir( &
            & profiles,          &
            & traj0_sta%angles,  &
            & coefs%coef,        &
            & opts%addpc,        &
            & coefs%coef_pccomp, &
            & chanprof,          &
            & calcemis,          &
            & emissivity_out,    &
            & errorstatus)
      IF (errorstatus == errorstatus_fatal) err = errorstatus_fatal
      THROWM( ERR .NE. 0 , "calcemis_ir")      
      traj0%reflectivity(:) = 1 - emissivity_out(:)
    ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
!Microwave
      CALL rttov_calcemis_mw( &
            & opts,                          &
            & profiles,                      &
            & traj0_sta%angles,              &
            & coefs%coef,                    &
            & chanprof,                      &
            & traj0_dyn%transmission_aux,    &
            & calcemis,                      &
            & emissivity_out,                &
            & traj0%reflectivity,            &
            & errorstatus)
      IF (errorstatus == errorstatus_fatal) err = errorstatus_fatal

      THROWM( ERR .NE. 0 , "calcemis_mw")

    ELSE
! Hires
      CALL rttov_calcemis_ir( &
            & profiles,          &
            & traj0_sta%angles,  &
            & coefs%coef,        &
            & opts%addpc,        &
            & coefs%coef_pccomp, &
            & chanprof,          &
            & calcemis,          &
            & emissivity_out,    &
            & errorstatus)
      IF (errorstatus == errorstatus_fatal) err = errorstatus_fatal
      THROWM( ERR .NE. 0 , "calcemis_ir")
      traj0%reflectivity(:) = 1 - emissivity_out(:)
    ENDIF

  ENDIF

!-------------------------------------------------------
!12. Compute Fresnel reflectance and sun glint - SURFACE
!-------------------------------------------------------

  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_refsun( &
          & profiles,         &
          & coefs%coef,       &
          & traj0%aux_prof,   &
          & traj0%sunglint,   &
          & traj0%raytracing)
    CALL rttov_fresnel( &
          & chanprof,        &
          & profiles,        &
          & coefs%coef,      &
          & traj0%sunglint,  &
          & traj0%fresnrefl)
  ENDIF

!---------------------------------------------------------
!13. integrate the radiative transfer equation - USER levs
!---------------------------------------------------------
! for new modeltop - '% layer' renamed '% air' (use of 'layer' now misleading)

  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)
  CALL rttov_integrate( &
        & addcosmic,                              &
        & opts,                                   &
        & traj0_dyn%nstreams,                     &
        & chanprof,                               &
        & emissivity_out,                         &
        & traj0%reflectivity,                     &
        & traj0%fresnrefl,                        &
        & traj0%sunglint,                         &
        & traj0_sta%sun,                          &
        & traj0_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream, &
        & profiles,                               &
        & traj0%aux_prof,                         &
        & coefs%coef,                             &
        & traj0%raytracing,                       &
        & traj0%ircld,                            &
        & radiancedata,                           &
        & traj0_sta%auxrad,                       &
        & traj0_dyn%auxrad_stream)

  IF (opts%addpc) THEN
    CALL rttov_pcscores( &
          & opts,                        &
          & chanprof,                    &
          & traj0_sta%chanprof_pc,       &
          & pccomp,                      &
          & coefs%coef_pccomp,           &
          & radiancedata)

    IF (opts%addradrec) THEN
      CALL rttov_reconstruct( &
            & traj0_sta%chanprof_in,       &
            & traj0_sta%chanprof_pc,       &
            & pccomp,                      &
            & coefs%coef_pccomp)
      CALL rttov_calcbt_pc(traj0_sta%chanprof_in, coefs%coef_pccomp, pccomp)
    ENDIF

  ENDIF

!---------------------
!14. deallocate memory
!---------------------


  IF (ltraj_dyn_dealloc) THEN

    CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
    THROW(ERR .NE. 0)

  ENDIF

  IF (ltraj_sta_dealloc) THEN

    CALL rttov_alloc_traj_sta (err, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, &
                              & 0_jpim, npcscores, channels_rec)
    THROW(err.ne.0)

  ENDIF

  CALL rttov_check_traj( &
        & ERR,           &
        & nprofiles,     &
        & nchannels,     &
        & opts,          &
        & nlevels,       &
        & coefs,         &
        & 0_jpim,        &
        & traj0 = traj0, &
        & traj1 = traj1, &
        & traj2 = traj)

  THROWM( ERR .NE. 0 , "deallocation check_traj")

  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)
  CATCH
  errorstatus = ERR
  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_direct
