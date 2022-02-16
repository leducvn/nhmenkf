SUBROUTINE rttov_add_raytracing(raytracing, raytracing1, raytracing2)
! Description:
!   Adds two raytracing structures
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
  USE rttov_types, ONLY : raytracing_type
  IMPLICIT NONE
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing1
  TYPE(raytracing_type), INTENT(IN)    :: raytracing2
!INTF_END
  raytracing%HGPL     = raytracing1%HGPL + raytracing2%HGPL
  raytracing%DAIR     = raytracing1%DAIR + raytracing2%DAIR
  raytracing%DMAIR    = raytracing1%DMAIR + raytracing2%DMAIR
  raytracing%R        = raytracing1%R + raytracing2%R
  raytracing%RATOESUN = raytracing1%RATOESUN + raytracing2%RATOESUN
  raytracing%RATOESAT = raytracing1%RATOESAT + raytracing2%RATOESAT
  raytracing%ZASUN    = raytracing1%ZASUN + raytracing2%ZASUN
  raytracing%ZASAT    = raytracing1%ZASAT + raytracing2%ZASAT
  raytracing%INT      = raytracing1%INT + raytracing2%INT
  raytracing%HL       = raytracing1%HL + raytracing2%HL
  raytracing%PPW      = raytracing1%PPW + raytracing2%PPW
  raytracing%DISPCO2  = raytracing1%DISPCO2 + raytracing2%DISPCO2
  raytracing%PATHSAT  = raytracing1%PATHSAT + raytracing2%PATHSAT
  raytracing%PATHSUN  = raytracing1%PATHSUN + raytracing2%PATHSUN
  raytracing%PATHEFF  = raytracing1%PATHEFF + raytracing2%PATHEFF
  raytracing%LTICK    = raytracing1%LTICK + raytracing2%LTICK
END SUBROUTINE 
