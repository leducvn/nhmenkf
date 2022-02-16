SUBROUTINE rttov_alloc_traj( &
            & ERR,       &
            & nprofiles, &
            & nchannels, &
            & opts,      &
            & nlevels,   &
            & coefs,     &
            & asw,       &
            & traj,      &
            & traj_tl,   &
            & traj_ad,   &
            & traj_k)
! Description:
!   Allocates/deallocates trajectory data for rttov_direct
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY :  &
       & rttov_options,  &
       & rttov_coefs,    &
       & rttov_traj
!INTF_OFF
  USE rttov_const, ONLY : ncldtyp, sensor_id_hi
  USE parkind1, ONLY : jplm
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nprofiles
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nchannels
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nlevels
  TYPE(rttov_options) , INTENT(IN)                 :: opts
  TYPE(rttov_coefs  ) , INTENT(IN) , TARGET        :: coefs                             ! Target attribute needed
  INTEGER(KIND=jpim)  , INTENT(OUT)                :: ERR
  TYPE(rttov_traj)    , OPTIONAL   , INTENT(INOUT) :: traj    , traj_tl, traj_ad, traj_k
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: asw
!INTF_END

#include "rttov_errorreport.h"


  TRY

  IF (Present(traj)) CALL alloc_traj( &
       & ERR,       &
       & nprofiles, &
       & nchannels, &
       & opts,      &
       & nlevels,   &
       & coefs,     &
       & asw,       &
       & traj)
  THROW(ERR.NE.0)
  IF (Present(traj_tl))     &
    &  CALL alloc_traj( &
       & ERR,       &
       & nprofiles, &
       & nchannels, &
       & opts,      &
       & nlevels,   &
       & coefs,     &
       & asw,       &
       & traj_tl)
  THROW(ERR.NE.0)
  IF (Present(traj_ad))     &
    &  CALL alloc_traj( &
       & ERR,       &
       & nprofiles, &
       & nchannels, &
       & opts,      &
       & nlevels,   &
       & coefs,     &
       & asw,       &
       & traj_ad)
  THROW(ERR.NE.0)
  IF (Present(traj_k)) THEN
    CALL alloc_traj( &
          & ERR,       &
          & nchannels, &
          & nchannels, &
          & opts,      &
          & nlevels,   &
          & coefs,     &
          & asw,       &
          & traj_k)
    THROW(ERR.NE.0)
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE alloc_traj( &
              & ERR,       &
              & nprofiles, &
              & nchannels, &
              & opts,      &
              & nlevels,   &
              & coefs,     &
              & asw,       &
              & traj)
    INTEGER(KIND=jpim)  , INTENT(IN)            :: nprofiles
    INTEGER(KIND=jpim)  , INTENT(IN)            :: nchannels
    TYPE(rttov_options) , INTENT(IN)            :: opts
    TYPE(rttov_coefs  ) , INTENT(IN)   , TARGET :: coefs
    INTEGER(KIND=jpim)  , INTENT(IN)            :: nlevels
    INTEGER(KIND=jpim)  , INTENT(OUT)           :: ERR
    TYPE(rttov_traj)    , INTENT(INOUT)         :: traj
    INTEGER(KIND=jpim)  , INTENT(IN)            :: asw
#include "rttov_alloc_prof.h"
#include "rttov_alloc_raytracing.h"
#include "rttov_alloc_predictor.h"
#include "rttov_alloc_opdp_path.h"
#include "rttov_alloc_ircld.h"
#include "rttov_alloc_aux_prof.h"
#include "rttov_alloc_sunglint.h"
#include "rttov_alloc_trans_scatt_ir.h"
    INTEGER(KIND=jpim)  :: nlayers
    TYPE(rttov_options) :: opts_COEF
    TRY
    nlayers   = nlevels - 1
    IF (asw .EQ. 1_jpim) THEN
      opts_COEF           = opts
      opts_COEF%addclouds = .FALSE.
      opts_COEF%addaerosl = .FALSE.
      ALLOCATE (traj%profiles_COEF(nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of traj%profiles_COEF")
! ON LEVELS
! ---------
      CALL rttov_alloc_prof( &
            & ERR,                            &
            & nprofiles,                      &
            & traj%profiles_COEF,             &
            & coefs%coef%nlevels,             &
            & opts_COEF,                      &
            & 1_jpim,                         &
            & coefs = coefs,                  &
            & blob = traj%profiles_COEF_blob)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_raytracing( &
            & ERR,                  &
            & nprofiles,            &
            & traj%raytracing_COEF, &
            & coefs%coef%nlevels,   &
            & 1_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_raytracing( &
            & ERR,             &
            & nprofiles,       &
            & traj%raytracing, &
            & nlevels,         &
            & 1_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_opdp_path( &
            & ERR,            &
            & traj%opdp_path, &
            & nlevels,        &
            & nchannels,      &
            & 1_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_opdp_path( &
            & ERR,                 &
            & traj%opdp_path_COEF, &
            & coefs%coef%nlevels,  &
            & nchannels,           &
            & 1_jpim,              &
            & init=.TRUE._jplm)
      THROW( ERR .NE. 0 )
! ON LAYERS
! ---------
!
      CALL rttov_alloc_predictor( &
            & ERR,             &
            & nprofiles,       &
            & traj%predictors, &
            & coefs%coef,      &
            & 1_jpim,          &
            & opts%addsolar)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_ircld( &
            & ERR,            &
            & nprofiles,      &
            & traj%ircld,     &
            & nlayers,        &
            & 1_jpim,         &
            & opts%addaerosl, &
            & opts%addclouds)
      THROW( ERR .NE. 0 )
      ALLOCATE (traj%reflectivity(nchannels), STAT = err)
      THROWM(err.ne.0,"Allocation of reflectivity failed")
      ALLOCATE (traj%fresnrefl(nchannels), STAT = err)
      THROWM(err.ne.0,"Allocation of fresnrefl failed")
      CALL rttov_alloc_aux_prof( &
            & err,                  &
            & nprofiles,            &
            & nlevels,              &
            & coefs%coef%id_sensor, &
            & traj%aux_prof,        &
            & opts,                 &
            & 1_jpim)
      THROW(err.ne.0)
      CALL rttov_alloc_aux_prof( &
            & err,                  &
            & nprofiles,            &
            & coefs%coef%nlevels,   &
            & coefs%coef%id_sensor, &
            & traj%aux_prof_COEF,   &
            & opts_COEF,            &
            & 1_jpim)
      THROW(err.ne.0)
      IF (opts%addaerosl .OR. opts%addclouds) THEN
        CALL rttov_alloc_trans_scatt_ir( &
              & err,                        &
              & traj%transmission_scatt_ir, &
              & nchannels,                  &
              & ncldtyp,                    &
              & nlayers,                    &
              & 1_jpim)
        THROWM( ERR .NE. 0 , "allocation of trans_scatt_ir")
      ENDIF
      IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
        CALL rttov_alloc_sunglint( &
              & err,                  &
              & traj%sunglint,        &
              & nprofiles,            &
              & coefs%coef%ws_nomega, &
              & 1_jpim)
        THROW(err.ne.0)
      ENDIF
      traj%coefs => coefs
      traj%nchannels = nchannels
      traj%opts      = opts
      traj%nlevels   = nlevels
      traj%nlayers   = traj%nlevels - 1
    ELSE
      IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
        CALL rttov_alloc_sunglint( &
              & err,                  &
              & traj%sunglint,        &
              & nprofiles,            &
              & coefs%coef%ws_nomega, &
              & 0_jpim)
        THROW(err.ne.0)
      ENDIF
      IF (opts%addaerosl .OR. opts%addclouds) THEN
        CALL rttov_alloc_trans_scatt_ir( &
              & err,                        &
              & traj%transmission_scatt_ir, &
              & nchannels,                  &
              & ncldtyp,                    &
              & nlayers,                    &
              & 0_jpim)
        THROWM( ERR .NE. 0 , "allocation of trans_scatt_ir")
      ENDIF
      CALL rttov_alloc_aux_prof( &
            & err,                  &
            & nprofiles,            &
            & coefs%coef%nlevels,   &
            & coefs%coef%id_sensor, &
            & traj%aux_prof_COEF,   &
            & opts_COEF,            &
            & 0_jpim)
      THROW(err.ne.0)
      CALL rttov_alloc_aux_prof( &
            & err,                  &
            & nprofiles,            &
            & nlevels,              &
            & coefs%coef%id_sensor, &
            & traj%aux_prof,        &
            & opts,                 &
            & 0_jpim)
      THROW(err.ne.0)
      DEALLOCATE (traj%fresnrefl, STAT = err)
      THROWM(err.ne.0,"DeAllocation of fresnrefl failed")
      DEALLOCATE (traj%reflectivity, STAT = err)
      THROWM(err.ne.0,"DeAllocation of reflectivity failed")
! ON LAYERS
! ---------
!
      CALL rttov_alloc_ircld( &
            & ERR,            &
            & nprofiles,      &
            & traj%ircld,     &
            & nlayers,        &
            & 0_jpim,         &
            & opts%addaerosl, &
            & opts%addclouds)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_opdp_path( &
            & ERR,                 &
            & traj%opdp_path_COEF, &
            & coefs%coef%nlayers,  &
            & nchannels,           &
            & 0_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_opdp_path( &
            & ERR,            &
            & traj%opdp_path, &
            & nlayers,        &
            & nchannels,      &
            & 0_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_predictor( &
            & ERR,             &
            & nprofiles,       &
            & traj%predictors, &
            & coefs%coef,      &
            & 0_jpim,          &
            & opts%addsolar)
      THROW( ERR .NE. 0 )
!
! ON LEVELS
! ---------
!
!
      CALL rttov_alloc_raytracing( &
            & ERR,             &
            & nprofiles,       &
            & traj%raytracing, &
            & nlevels,         &
            & 0_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_raytracing( &
            & ERR,                  &
            & nprofiles,            &
            & traj%raytracing_COEF, &
            & coefs%coef%nlevels,   &
            & 0_jpim)
      THROW( ERR .NE. 0 )
      CALL rttov_alloc_prof( &
            & ERR,                            &
            & nprofiles,                      &
            & traj%profiles_COEF,             &
            & coefs%coef%nlevels,             &
            & opts_COEF,                      &
            & 0_jpim,                         &
            & blob = traj%profiles_COEF_blob)
      THROW( ERR .NE. 0 )
      DEALLOCATE (traj%profiles_COEF, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of traj%profiles_COEF")
    ENDIF
    CATCH
  END SUBROUTINE
END SUBROUTINE 
