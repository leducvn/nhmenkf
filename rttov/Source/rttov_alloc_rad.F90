!
SUBROUTINE rttov_alloc_rad( &
            & ERR,       &
            & nchannels, &
            & radiance,  &
            & nlayers,   &
            & asw,       &
            & init)
! Description:
! allocation/deallocation of a radiance structure
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
!  2.0       02/12/2009  Marco matricardi:Added Principal component capability
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
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : radiance_type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim) , INTENT(OUT)          :: ERR      ! return code
  INTEGER(KIND=jpim) , INTENT(IN)           :: nchannels! number of channels
  INTEGER(KIND=jpim) , INTENT(IN)           :: nlayers  ! number of levels
  TYPE(radiance_Type), INTENT(INOUT)        :: radiance ! radiances
  INTEGER(KIND=jpim) , INTENT(IN)           :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm) , OPTIONAL, INTENT(IN) :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_rad.h"
! Local arrays and scalars:
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1   = .FALSE.
  IF (Present(init)) init1 = init
!  Allocate section
  IF (asw .EQ. 1) THEN
!
    ALLOCATE (radiance%clear(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%clear")
    ALLOCATE (radiance%cloudy(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%cloudy")
    ALLOCATE (radiance%total(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%total")
    ALLOCATE (radiance%bt(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%bt")
    ALLOCATE (radiance%bt_clear(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%bt_clear")
    ALLOCATE (radiance%upclear(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%upclear")
    ALLOCATE (radiance%dnclear(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%dnclear")
    ALLOCATE (radiance%reflclear(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%reflclear")
    ALLOCATE (radiance%overcast(nlayers, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%overcast")
    ALLOCATE (radiance%up(nlayers, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%up")
    ALLOCATE (radiance%down(nlayers, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%down")
    ALLOCATE (radiance%surf(nlayers, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of radiance%surf")
    
    IF (init1) THEN
      CALL rttov_init_rad(radiance)
    ENDIF
  ELSE
! deallocate radiance results arrays with number of channels
    DEALLOCATE (radiance%clear, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%clear")
    DEALLOCATE (radiance%cloudy, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%cloudy")
    DEALLOCATE (radiance%total, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%total")
    DEALLOCATE (radiance%bt, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%bt")
    DEALLOCATE (radiance%bt_clear, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%bt_clear")
    DEALLOCATE (radiance%upclear, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%upclear")
    DEALLOCATE (radiance%dnclear, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%dnclear")
    DEALLOCATE (radiance%reflclear, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%reflclear")
    DEALLOCATE (radiance%overcast, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%overcast")
    DEALLOCATE (radiance%up, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%up")
    DEALLOCATE (radiance%down, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%down")
    DEALLOCATE (radiance%surf, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of radiance%surf")
    NULLIFY (radiance%clear)
    NULLIFY (radiance%cloudy)
    NULLIFY (radiance%total)
    NULLIFY (radiance%bt)
    NULLIFY (radiance%bt_clear)
    NULLIFY (radiance%upclear)
    NULLIFY (radiance%dnclear)
    NULLIFY (radiance%reflclear)
    NULLIFY (radiance%overcast)
    NULLIFY (radiance%up)
    NULLIFY (radiance%down)
    NULLIFY (radiance%surf)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_rad
