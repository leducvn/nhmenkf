!
SUBROUTINE rttov_alloc_ircld( &
            & ERR,       &
            & nircld,    &
            & irclds,    &
            & nlayers,   &
            & asw,       &
            & addaerosl, &
            & addclouds, &
            & init)
! Description:
! allocation/deallocation of a ircld structure
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
!  1.0       22/20/2002  Creation
!  1.1       02/12/2009  Increased the dimensions of cldtyp (Marco Matricardi)
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
  USE rttov_types, ONLY : ircld_Type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : ncldtyp
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim), INTENT(OUT)             :: ERR      ! return code
  INTEGER(KIND=jpim), INTENT(IN)              :: nircld   !
  INTEGER(KIND=jpim), INTENT(IN)              :: nlayers  ! number of layers
  TYPE(ircld_Type)  , INTENT(INOUT)           :: irclds
  INTEGER(KIND=jpim), INTENT(IN)              :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: addclouds
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_ircld.h"
! Local Arrays and Scalars:
  LOGICAL(KIND=jplm) :: addclouds1
  LOGICAL(KIND=jplm) :: addaerosl1
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  addclouds1 = .FALSE.
  addaerosl1 = .FALSE.
  init1      = .FALSE.
  IF (Present(addclouds)) addclouds1 = addclouds
  IF (Present(addaerosl)) addaerosl1 = addaerosl
  IF (Present(init)     ) init1      = init
  IF (asw .EQ. 1) THEN
    NULLIFY (irclds%tave)
    NULLIFY (irclds%wmixave)
    NULLIFY (irclds%xpresave)
    NULLIFY (irclds%esw)
    NULLIFY (irclds%esi)
    NULLIFY (irclds%ppv)
    NULLIFY (irclds%icldarr)
    NULLIFY (irclds%xstrref1)
    NULLIFY (irclds%xstrref2)
    NULLIFY (irclds%cldtyp)
    NULLIFY (irclds%xstr)
    NULLIFY (irclds%xstrminref)
    NULLIFY (irclds%xstrref)
    NULLIFY (irclds%cldcfr)
    NULLIFY (irclds%maxcov)
    NULLIFY (irclds%xstrmax)
    NULLIFY (irclds%xstrmin)
    NULLIFY (irclds%a)
    NULLIFY (irclds%ntotref)
    NULLIFY (irclds%indexstr)
    NULLIFY (irclds%icount1ref)
    NULLIFY (irclds%iloopin)
    NULLIFY (irclds%flag)
    NULLIFY (irclds%iflag)
    NULLIFY (irclds%nstream)
    NULLIFY (irclds%nstreamref)
    NULLIFY (irclds%iloop)
    NULLIFY (irclds%icount)
    NULLIFY (irclds%icounstr)
    NULLIFY (irclds%icount1)
    NULLIFY (irclds%xstrclr)
  ENDIF
  IF (asw .EQ. 1) THEN
!
    IF (addaerosl1) THEN
      ALLOCATE (irclds%tave(nlayers, nircld), irclds%wmixave(nlayers, nircld), irclds%xpresave(nlayers, nircld),      &
        & irclds%esw(nlayers, nircld), irclds%esi(nlayers, nircld), irclds%ppv(nlayers, nircld), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of irclds")
    ENDIF
    IF (addclouds1) THEN
      ALLOCATE (irclds%icldarr(2*nlayers, nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % icldarr ")
      ALLOCATE (irclds%xstrref1(2*nlayers, 2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrref1 ")
      ALLOCATE (irclds%xstrref2(2*nlayers, 2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrref2 ")
      ALLOCATE (irclds%cldtyp(ncldtyp, nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % cldtyp ")
      ALLOCATE (irclds%xstr(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstr ")
      ALLOCATE (irclds%xstrminref(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrminref ")
      ALLOCATE (irclds%xstrref(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrref ")
      ALLOCATE (irclds%cldcfr(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % cldcfr ")
      ALLOCATE (irclds%maxcov(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % maxcov ")
      ALLOCATE (irclds%xstrmax(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrmax ")
      ALLOCATE (irclds%xstrmin(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % xstrmin ")
      ALLOCATE (irclds%a(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % a ")
      ALLOCATE (irclds%ntotref(nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % ntotref ")
      ALLOCATE (irclds%indexstr(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % indexstr ")
      ALLOCATE (irclds%icount1ref(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % icount1ref ")
      ALLOCATE (irclds%iloopin(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % iloopin ")
      ALLOCATE (irclds%flag(2*nlayers, nircld), STAT = err)
      THROWM( ERR .NE. 0 , "allocation of irclds % flag ")
      ALLOCATE (irclds%iflag(2*nlayers, nircld), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of irclds % iflag")
    ENDIF
    ALLOCATE (irclds%nstream(nircld), irclds%nstreamref(nircld), irclds%iloop(nircld), irclds%icount(nircld),      &
      & irclds%icounstr(nircld), irclds%icount1(nircld), irclds%xstrclr(nircld), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of irclds % ")
    irclds%NSTREAM = 0_jpim
    irclds%XSTRCLR = 1._jprb
  ENDIF
  IF (asw .EQ. 0) THEN
    IF (Associated(irclds%tave)) THEN
      DEALLOCATE (irclds%tave, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % tave")
    ENDIF
    IF (Associated(irclds%wmixave)) THEN
      DEALLOCATE (irclds%wmixave, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % wmixave")
    ENDIF
    IF (Associated(irclds%xpresave)) THEN
      DEALLOCATE (irclds%xpresave, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xpresave")
    ENDIF
    IF (Associated(irclds%esw)) THEN
      DEALLOCATE (irclds%esw, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % esw")
    ENDIF
    IF (Associated(irclds%esi)) THEN
      DEALLOCATE (irclds%esi, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % esi")
    ENDIF
    IF (Associated(irclds%ppv)) THEN
      DEALLOCATE (irclds%ppv, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % ppv")
    ENDIF
    IF (Associated(irclds%icldarr)) THEN
      DEALLOCATE (irclds%icldarr, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % icldarr")
    ENDIF
    IF (Associated(irclds%xstrref1)) THEN
      DEALLOCATE (irclds%xstrref1, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrref1")
    ENDIF
    IF (Associated(irclds%xstrref2)) THEN
      DEALLOCATE (irclds%xstrref2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrref2")
    ENDIF
    IF (Associated(irclds%cldtyp)) THEN
      DEALLOCATE (irclds%cldtyp, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % cldtyp")
    ENDIF
    IF (Associated(irclds%xstr)) THEN
      DEALLOCATE (irclds%xstr, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstr")
    ENDIF
    IF (Associated(irclds%xstrminref)) THEN
      DEALLOCATE (irclds%xstrminref, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrminref")
    ENDIF
    IF (Associated(irclds%xstrref)) THEN
      DEALLOCATE (irclds%xstrref, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrref")
    ENDIF
    IF (Associated(irclds%cldcfr)) THEN
      DEALLOCATE (irclds%cldcfr, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % cldcfr")
    ENDIF
    IF (Associated(irclds%maxcov)) THEN
      DEALLOCATE (irclds%maxcov, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % maxcov")
    ENDIF
    IF (Associated(irclds%xstrmax)) THEN
      DEALLOCATE (irclds%xstrmax, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrmax")
    ENDIF
    IF (Associated(irclds%xstrmin)) THEN
      DEALLOCATE (irclds%xstrmin, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % xstrmin")
    ENDIF
    IF (Associated(irclds%a)) THEN
      DEALLOCATE (irclds%a, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % a")
    ENDIF
    IF (Associated(irclds%ntotref)) THEN
      DEALLOCATE (irclds%ntotref, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % ntotref")
    ENDIF
    IF (Associated(irclds%indexstr)) THEN
      DEALLOCATE (irclds%indexstr, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % indexstr")
    ENDIF
    IF (Associated(irclds%icount1ref)) THEN
      DEALLOCATE (irclds%icount1ref, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % icount1ref")
    ENDIF
    IF (Associated(irclds%iloopin)) THEN
      DEALLOCATE (irclds%iloopin, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % iloopin")
    ENDIF
    IF (Associated(irclds%flag)) THEN
      DEALLOCATE (irclds%flag, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % flag")
    ENDIF
    IF (Associated(irclds%iflag)) THEN
      DEALLOCATE (irclds%iflag, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of irclds % iflag")
    ENDIF
    DEALLOCATE (                                                                                                        &
      & irclds%nstream, irclds%nstreamref, irclds%iloop, irclds%icount, irclds%icounstr, irclds%icount1, irclds%xstrclr &
      & , STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of irclds")
!
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    irclds%nstream    = 0_jpim
    irclds%nstreamref = 0_jpim
    irclds%iloop      = 0_jpim
    irclds%icount     = 0_jpim
    irclds%icounstr   = 0_jpim
    irclds%icount1    = 0_jpim
    CALL rttov_init_ircld(irclds)
  ENDIF
  IF (asw .EQ. 0) THEN
    NULLIFY (irclds%tave)
    NULLIFY (irclds%wmixave)
    NULLIFY (irclds%xpresave)
    NULLIFY (irclds%esw)
    NULLIFY (irclds%esi)
    NULLIFY (irclds%ppv)
    NULLIFY (irclds%icldarr)
    NULLIFY (irclds%xstrref2)
    NULLIFY (irclds%cldtyp)
    NULLIFY (irclds%xstr)
    NULLIFY (irclds%xstrminref)
    NULLIFY (irclds%xstrref)
    NULLIFY (irclds%cldcfr)
    NULLIFY (irclds%maxcov)
    NULLIFY (irclds%xstrmax)
    NULLIFY (irclds%xstrmin)
    NULLIFY (irclds%a)
    NULLIFY (irclds%ntotref)
    NULLIFY (irclds%indexstr)
    NULLIFY (irclds%icount1ref)
    NULLIFY (irclds%iloopin)
    NULLIFY (irclds%flag)
    NULLIFY (irclds%iflag)
    NULLIFY (irclds%nstream)
    NULLIFY (irclds%nstreamref)
    NULLIFY (irclds%iloop)
    NULLIFY (irclds%icount)
    NULLIFY (irclds%icounstr)
    NULLIFY (irclds%icount1)
    NULLIFY (irclds%xstrclr)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_ircld
