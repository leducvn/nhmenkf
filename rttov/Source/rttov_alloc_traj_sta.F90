SUBROUTINE rttov_alloc_traj_sta (err, traj_sta, opts, coef, nlevels, nchannels, nprofiles, asw, npcscores, channels_rec)
! Description:
!   Allocates/deallocates static trajectory data
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
#include "throw.h"

USE rttov_types, ONLY : &
 & rttov_traj_sta, &
 & rttov_coef,     &
 & rttov_options

USE parkind1, ONLY : &
  & jpim, jprb

IMPLICIT NONE


INTEGER(KIND=jpim),   INTENT(OUT)   :: err
TYPE(rttov_traj_sta), INTENT(INOUT) :: traj_sta
TYPE(rttov_options),  INTENT(IN)    :: opts
TYPE(rttov_coef),     INTENT(IN)    :: coef
INTEGER(KIND=jpim),   INTENT(IN)    :: nlevels
INTEGER(KIND=jpim),   INTENT(IN)    :: nchannels
INTEGER(KIND=jpim),   INTENT(IN)    :: nprofiles
INTEGER(KIND=jpim),   INTENT(IN)    :: asw
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: npcscores
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: channels_rec(:)

!INTF_END

#include "rttov_alloc_prof.h"
#include "rttov_alloc_auxrad.h"
#include "rttov_alloc_pc_dimensions.h"

INTEGER(KIND=jpim)  :: nlayers
TYPE(rttov_options) :: opts_COEF


TRY

  nlayers = nlevels - 1

  opts_COEF = opts
  opts_COEF%addclouds = .false.
  opts_COEF%addaerosl = .false.

  IF (asw .eq. 1) THEN

    IF (opts%addpc) THEN
   
      CALL rttov_alloc_pc_dimensions( &
            & err,                   &
            & opts,                  &
            & npcscores,             &
            & nprofiles,             &
            & traj_sta%chanprof_in,  &
            & traj_sta%chanprof_pc,  &
            & 1_jpim,                &
            & channels_rec = channels_rec)

     THROW(err.ne.0)

    ENDIF

    ALLOCATE (    &
            & traj_sta%sun(nchannels),                                &
            & traj_sta%tau_ref(nlevels, nchannels),                   &
            & traj_sta%tau_ref_surf(nchannels),                       &
            & traj_sta%tau_surf(nchannels),                           &
            & traj_sta%tausun_ref(nlevels, nchannels),                &
            & traj_sta%tausun_ref_surf(nchannels),                    &
            & traj_sta%tausun_level(nlevels, nchannels),              &
            & traj_sta%tausun_surf(nchannels),                        &
            & traj_sta%opdp_ref_COEF(coef%nlayers, nchannels),        &
            & traj_sta%od_level(nlevels, nchannels),                  &
            & traj_sta%tau_level(nlevels, nchannels),                 &
            & traj_sta%opdpsun_ref_COEF(coef%nlayers , nchannels),    &
            & traj_sta%odsun_level(nlevels, nchannels),               &
            & traj_sta%odsun_singlelayer(nlayers, nchannels),         &
            & traj_sta%od_frac(nchannels),                            &
            & traj_sta%angles(nprofiles),                             & 
            & traj_sta%angles_COEF(nprofiles),                        &
            & traj_sta%profiles_COEF_ref(nprofiles),                  &
            & STAT = err)
    THROWM(err.ne.0,"Allocation of traj_sta failed")

    IF (opts%apply_reg_limits) THEN

      CALL rttov_alloc_prof (err, nprofiles, traj_sta%profiles_COEF_ref, coef % nlevels, opts_COEF, &
                            & 1_jpim, blob = traj_sta%profiles_COEF_blob_ref)
      THROW(err.ne.0)

    ENDIF

    CALL rttov_alloc_auxrad( &
          & err,                 &
          & traj_sta%auxrad,     &
          & nlevels,             &
          & nchannels,           &
          & 1_jpim)
    THROW(err.ne.0)
  
  ENDIF

  IF (asw .eq. 0) THEN

    CALL rttov_alloc_auxrad( &
          & err,                 &
          & traj_sta%auxrad,     &
          & nlevels,             &
          & nchannels,           &
          & 0_jpim)
    THROW(err.ne.0)
  
    IF (opts%apply_reg_limits) THEN

      CALL rttov_alloc_prof (err, nprofiles, traj_sta%profiles_COEF_ref, coef % nlevels, opts_COEF, &
                            & 0_jpim, blob = traj_sta%profiles_COEF_blob_ref)
      THROW(err.ne.0)

    ENDIF
  
    DEALLOCATE (    &
            & traj_sta%sun,                  &
            & traj_sta%tau_ref,              &
            & traj_sta%tau_ref_surf,         &
            & traj_sta%tau_surf,             &
            & traj_sta%tausun_ref,           &
            & traj_sta%tausun_ref_surf,      &
            & traj_sta%tausun_level,         &
            & traj_sta%tausun_surf,          &
            & traj_sta%opdp_ref_COEF,        &
            & traj_sta%od_level,             &
            & traj_sta%tau_level,            &
            & traj_sta%opdpsun_ref_COEF,     &
            & traj_sta%odsun_level,          &
            & traj_sta%odsun_singlelayer,    &
            & traj_sta%od_frac,              &
            & traj_sta%angles,               &
            & traj_sta%angles_COEF,          &
            & traj_sta%profiles_COEF_ref,    &
            & STAT = err)
    THROWM(err.ne.0,"DeAllocation of traj_sta failed")

    IF (opts%addpc) THEN
   
      CALL rttov_alloc_pc_dimensions( &
            & err,                   &
            & opts,                  &
            & npcscores,             &
            & nprofiles,             &
            & traj_sta%chanprof_in,  &
            & traj_sta%chanprof_pc,  &
            & 0_jpim,                &
            & channels_rec = channels_rec)

     THROW(err.ne.0)

    ENDIF

  ENDIF

CATCH

END SUBROUTINE
