SUBROUTINE rttov_copy_raytracing(raytracing1, raytracing2)
! Description:
!   Copy raytracing structure
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
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing1
  TYPE(raytracing_type), INTENT(IN)    :: raytracing2
!INTF_END
  raytracing1%HGPL     = raytracing2%HGPL
  raytracing1%DAIR     = raytracing2%DAIR
  raytracing1%DMAIR    = raytracing2%DMAIR
  raytracing1%R        = raytracing2%R
  raytracing1%RATOESUN = raytracing2%RATOESUN
  raytracing1%RATOESAT = raytracing2%RATOESAT
  raytracing1%ZASUN    = raytracing2%ZASUN
  raytracing1%ZASAT    = raytracing2%ZASAT
  raytracing1%INT      = raytracing2%INT
  raytracing1%HL       = raytracing2%HL
  raytracing1%PPW      = raytracing2%PPW
  raytracing1%DISPCO2  = raytracing2%DISPCO2
  raytracing1%PATHSAT  = raytracing2%PATHSAT
  raytracing1%PATHSUN  = raytracing2%PATHSUN
  raytracing1%PATHEFF  = raytracing2%PATHEFF
  raytracing1%LTICK    = raytracing2%LTICK
END SUBROUTINE 
