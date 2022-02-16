!
SUBROUTINE rttov_setgeometry_k( &
            & opts,         &
            & chanprof,     &
            & profiles,     &
            & profiles_k,   &
            & aux,          &
            & coef,         &
            & angles,       &
            & raytracing,   &
            & raytracing_k)
! Description:
! compute all profile related viewing geometry
! The only profile input value is profile%zenangle (zenith angle)
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
!  1.0       01/06/2005  New code. Marco Matricardi (ECMWF).
!  1.1       14/10/2010  Remove rt8_mode (J Hocking)
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
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & profile_Type,    &
       & profile_aux,     &
       & geometry_Type,   &
       & raytracing_type, &
       & rttov_options
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE parkind1, ONLY : JPRB
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type   ), INTENT(IN)    :: profiles  (:)             ! profile
  TYPE(profile_Type   ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                      ! coefficient
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_Type  ), INTENT(IN)    :: angles(size(profiles))    ! angles
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
!INTF_END
#include "rttov_locpat_k.h"
!local arguments
  INTEGER(KIND=jpim) :: i           ! loop index
  INTEGER(KIND=jpim) :: nchannels   ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  CALL RTTOV_LOCPAT_K( &
        & OPTS,         &
        & CHANPROF,     &
        & PROFILES,     &
        & PROFILES_K,   &
        & AUX,          &
        & COEF,         &
        & ANGLES,       &
        & RAYTRACING,   &
        & RAYTRACING_K)
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setgeometry_k
