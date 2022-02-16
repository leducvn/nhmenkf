!
SUBROUTINE rttov_setgeometry( &
            & opts,       &
            & profiles,   &
            & aux,        &
            & coef,       &
            & angles,     &
            & raytracing)! inout, optional
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       02/01/2003  Added more comments (R Saunders)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       01/06/2005  Marco Matricardi (ECMWF):
!               --       Added raytracing_type to imported type definitions
!               --       Added subroutine LOCPAT for the computation of the
!               --       altitude dependent local zenith angle
!  1.4       16/02/2006  Modified to make backward compatible with RTTOV-8
!                        (R Saunders)
!  1.5       11/11/2007  Made 'raytracing' and 'aux' optional to allow (if
!                        needed, as by RTTOV-SCATT) only 'angles' output
!  1.6       14/10/2010  Remove rt8_mode (J Hocking)
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
       & rttov_coef,      &
       & profile_Type,    &
       & profile_aux,     &
       & geometry_Type,   &
       & raytracing_type, &
       & rttov_options
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE parkind1, ONLY : JPRB
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)              :: opts
  TYPE(profile_Type   ), INTENT(IN)              :: profiles(:)           ! profile
  TYPE(rttov_coef     ), INTENT(IN)              :: coef                  ! coefficient
  TYPE(profile_aux    ), INTENT(IN)   , OPTIONAL :: aux
  TYPE(geometry_Type  ), INTENT(OUT)             :: angles(size(profiles))! angles
  TYPE(raytracing_type), INTENT(INOUT), OPTIONAL :: raytracing
!INTF_END
#include "rttov_locpat.h"
!local arguments
  INTEGER(KIND=jpim) :: i           ! loop index
  INTEGER(KIND=jpim) :: nprofiles   ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
! local
!- End of header --------------------------------------------------------
!Notes on notation:
! zen  => zenith angle
!   (definition: angle at surface between view path to satellite and zenith)
! view => view angle
!   (definition: angle at the satellite between view path and nadir)
! _sq = square of given value
! _sqrt = square root of given value
! _minus1 = given value - 1
! trigonometric function abbreviations have their usual meanings
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(profiles)
  DO i = 1, nprofiles
    angles(i)%sinzen           = Sin(profiles(i)%zenangle * deg2rad)
    angles(i)%sinzen_sq        = angles(i)%sinzen * angles(i)%sinzen
    angles(i)%coszen           = Cos(profiles(i)%zenangle * deg2rad)
    angles(i)%coszen_sq        = angles(i)%coszen * angles(i)%coszen
    angles(i)%seczen           = 1.0_JPRB / Abs(angles(i)%coszen)
    angles(i)%seczen_sq        = angles(i)%seczen * angles(i)%seczen
    angles(i)%seczen_sqrt      = Sqrt(angles(i)%seczen)
    angles(i)%seczen_minus1    = angles(i)%seczen - 1.0_JPRB
    angles(i)%seczen_minus1_sq = angles(i)%seczen_minus1 * angles(i)%seczen_minus1
    angles(i)%sinview          = angles(i)%sinzen / coef%ratoe
    angles(i)%sinview_sq       = angles(i)%sinview * angles(i)%sinview
    angles(i)%cosview_sq       = 1.0_JPRB - angles(i)%sinview_sq
    angles(i)%normzen          = profiles(i)%zenangle / 60.0_JPRB                 !normalized zenith angle
    angles(i)%viewang          = asin(angles(i)%sinview) / deg2rad
    IF (Present(raytracing)) raytracing%pathsat(:, i) = angles(i)%seczen
  ENDDO
  IF (Present(raytracing) .AND. Present(aux)) THEN
    CALL rttov_locpat( &
          & opts,       &
          & profiles,   &
          & aux,        &
          & coef,       &
          & angles,     &
          & raytracing)
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setgeometry
