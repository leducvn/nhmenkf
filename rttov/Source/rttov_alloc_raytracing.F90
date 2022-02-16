!
SUBROUTINE rttov_alloc_raytracing( &
            & ERR,          &
            & nraytracings, &
            & raytracings,  &
            & nlevels,      &
            & asw,          &
            & init)
! Description:
! allocation/deallocation of a raytracing structure
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
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code
!  1.1       11/10/2007  Add addclouds, addaerosl, init logicals
!                        nullify unused pointers P.Marguinaud
!  1.2       02/12/2009  Pathsat, Patsun, Patheff and Ltick are allocated
!                        on number of layers (Marco Matricardi).
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : raytracing_Type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)   , INTENT(OUT)             :: ERR         ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nraytracings
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nlevels
  TYPE(raytracing_Type), INTENT(INOUT)           :: raytracings
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_raytracing.h"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (Present(init)) init1 = init
!  Allocate section
  IF (asw .EQ. 1) THEN
    ALLOCATE (raytracings%HGPL(nlevels + 1, nraytracings), raytracings%DAIR(nlevels, nraytracings),      &
      & raytracings%DMAIR(nlevels, nraytracings), raytracings%R(nlevels, nraytracings),                  &
      & raytracings%RATOESUN(nlevels, nraytracings), raytracings%RATOESAT(nlevels, nraytracings),        &
      & raytracings%ZASUN(nlevels, nraytracings), raytracings%ZASAT(nlevels, nraytracings),              &
      & raytracings%INT(nlevels, nraytracings), raytracings%HL(nlevels, nraytracings),                   &
      & raytracings%PPW(nlevels, nraytracings), raytracings%DISPCO2(nlevels, nraytracings),              &
      & raytracings%PATHSAT(nlayers, nraytracings), raytracings%PATHSUN(nlayers, nraytracings),          &
      & raytracings%PATHEFF(nlayers, nraytracings), raytracings%LTICK(nlayers, nraytracings), STAT = err)
    THROWM(err.ne.0,"Allocation of raytracings failed")
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (raytracings%HGPL, raytracings%DAIR, raytracings%DMAIR, raytracings%R, raytracings%RATOESUN,            &
      & raytracings%RATOESAT, raytracings%ZASUN, raytracings%ZASAT, raytracings%INT, raytracings%HL, raytracings%PPW,  &
      & raytracings%DISPCO2, raytracings%PATHSAT, raytracings%PATHSUN, raytracings%PATHEFF, raytracings%LTICK, STAT = err)
    THROWM(err.ne.0,"DeAllocation of raytracings failed")
  ENDIF
  IF (asw .EQ. 0) THEN
    NULLIFY (raytracings%LTICK)
    NULLIFY (raytracings%HGPL)
    NULLIFY (raytracings%DAIR)
    NULLIFY (raytracings%DMAIR)
    NULLIFY (raytracings%R)
    NULLIFY (raytracings%RATOESUN)
    NULLIFY (raytracings%RATOESAT)
    NULLIFY (raytracings%ZASUN)
    NULLIFY (raytracings%ZASAT)
    NULLIFY (raytracings%INT)
    NULLIFY (raytracings%HL)
    NULLIFY (raytracings%PPW)
    NULLIFY (raytracings%DISPCO2)
    NULLIFY (raytracings%PATHSAT)
    NULLIFY (raytracings%PATHSUN)
    NULLIFY (raytracings%PATHEFF)
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    CALL rttov_init_raytracing(raytracings)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_raytracing
