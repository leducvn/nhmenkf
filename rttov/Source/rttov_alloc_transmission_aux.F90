!
SUBROUTINE rttov_alloc_transmission_aux( &
            & ERR,              &
            & transmission_aux, &
            & nlayers,          &
            & nchannels,        &
            & asw,              &
            & nstreams,         &
            & init)
! Description:
! allocation/deallocation of a profile structure
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
!  1.2       03/11/2009  Transmittances / optical depths on levels (A Geer)
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
  USE rttov_types, ONLY : transmission_Type_aux
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)         , INTENT(OUT)             :: ERR             ! return code
  TYPE(transmission_type_aux), INTENT(INOUT)           :: transmission_aux
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nlayers         ! number of levels
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)         , INTENT(IN)              :: asw             ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nstreams
  LOGICAL(KIND=jplm)         , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_transmission_aux.h"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    NULLIFY (transmission_aux%fac)
    NULLIFY (transmission_aux%surf_fac)
    NULLIFY (transmission_aux%tau_surf)
    NULLIFY (transmission_aux%tau_surf_r)
    NULLIFY (transmission_aux%tau_level_r)
    NULLIFY (transmission_aux%tau_level)
    NULLIFY (transmission_aux%od_singlelayer)
    NULLIFY (transmission_aux%od_singlelayer_r)
    NULLIFY (transmission_aux%od_level)
    NULLIFY (transmission_aux%odsun_singlelayer)
    NULLIFY (transmission_aux%tausun_surf)
    NULLIFY (transmission_aux%tausun_level)
    NULLIFY (transmission_aux%od_sfrac)
    NULLIFY (transmission_aux%od_sfrac_r)
    NULLIFY (transmission_aux%odsun_sfrac)
    NULLIFY (transmission_aux%od_frac_ac)
    NULLIFY (transmission_aux%odsun_frac_ac)
    NULLIFY (transmission_aux%tau_surf_ac)
    NULLIFY (transmission_aux%tau_surf_acsun)
    NULLIFY (transmission_aux%tau_ref_surf_ac)
    NULLIFY (transmission_aux%tau_ref_surf_acsun)
    NULLIFY (transmission_aux%od_frac_t)
    NULLIFY (transmission_aux%odsun_frac_t)
    NULLIFY (transmission_aux%tau_surf_t)
    NULLIFY (transmission_aux%tausun_surf_t)
    NULLIFY (transmission_aux%tau_ref_surf_t)
    NULLIFY (transmission_aux%tausun_ref_surf_t)
    NULLIFY (transmission_aux%refl_norm)
  ENDIF
  IF (asw .EQ. 1) THEN
    ALLOCATE (transmission_aux%fac(2, nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % fac")
    ALLOCATE (transmission_aux%surf_fac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % surf_fac")
    ALLOCATE (transmission_aux%tau_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_level")
    ALLOCATE (transmission_aux%tau_level_r(nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_level_r")
    ALLOCATE (transmission_aux%od_singlelayer(nlayers, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_singlelayer")
    ALLOCATE (transmission_aux%od_singlelayer_r(nlayers, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_singlelayer_r")
    ALLOCATE (transmission_aux%od_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_level")
    ALLOCATE (transmission_aux%odsun_singlelayer(nlayers, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % odsun_singlelayer")
    ALLOCATE (transmission_aux%tausun_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tausun_level")
    ALLOCATE (transmission_aux%od_sfrac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_sfrac")
    ALLOCATE (transmission_aux%od_sfrac_r(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_sfrac_r")
    ALLOCATE (transmission_aux%odsun_sfrac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % odsun_sfrac")
    ALLOCATE (transmission_aux%tau_surf(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_surf")
    ALLOCATE (transmission_aux%tau_surf_r(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_surf_r")
    ALLOCATE (transmission_aux%tausun_surf(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tausun_surf")
    ALLOCATE (transmission_aux%od_frac_ac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_frac_ac")
    ALLOCATE (transmission_aux%odsun_frac_ac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % odsun_frac_ac")
    ALLOCATE (transmission_aux%tau_surf_ac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_surf_ac")
    ALLOCATE (transmission_aux%tau_surf_acsun(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_surf_acsun")
    ALLOCATE (transmission_aux%tau_ref_surf_ac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_ref_surf_ac")
    ALLOCATE (transmission_aux%tau_ref_surf_acsun(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_ref_surf_acsun")
    ALLOCATE (transmission_aux%od_frac_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % od_frac_t")
    ALLOCATE (transmission_aux%odsun_frac_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % odsun_frac_t")
    ALLOCATE (transmission_aux%tau_surf_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_surf_t")
    ALLOCATE (transmission_aux%tausun_surf_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tausun_surf_t")
    ALLOCATE (transmission_aux%tau_ref_surf_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tau_ref_surf_t")
    ALLOCATE (transmission_aux%tausun_ref_surf_t(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % tausun_ref_surf_t")
    ALLOCATE (transmission_aux%refl_norm(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % refl_norm")
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (transmission_aux%fac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%fac")
    DEALLOCATE (transmission_aux%surf_fac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%fac")
    DEALLOCATE (transmission_aux%tau_level_r, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_level")
    DEALLOCATE (transmission_aux%tau_level, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_level")
    DEALLOCATE (transmission_aux%od_singlelayer, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_singlelayer")
    DEALLOCATE (transmission_aux%od_singlelayer_r, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_singlelayer")
    DEALLOCATE (transmission_aux%od_level, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_level")
    DEALLOCATE (transmission_aux%odsun_singlelayer, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%odsun_singlelayer")
    DEALLOCATE (transmission_aux%tausun_level, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tausun_level")
    DEALLOCATE (transmission_aux%od_sfrac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_sfrac")
    DEALLOCATE (transmission_aux%od_sfrac_r, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_sfrac_r")
    DEALLOCATE (transmission_aux%odsun_sfrac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%odsun_sfrac")
    DEALLOCATE (transmission_aux%tau_surf, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_surf")
    DEALLOCATE (transmission_aux%tau_surf_r, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_surf_r")
    DEALLOCATE (transmission_aux%tausun_surf, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tausun_surf")
    DEALLOCATE (transmission_aux%od_frac_ac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_frac_ac ")
    DEALLOCATE (transmission_aux%odsun_frac_ac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%odsun_frac_ac")
    DEALLOCATE (transmission_aux%tau_surf_ac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_surf_ac")
    DEALLOCATE (transmission_aux%tau_surf_acsun, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_surf_acsun")
    DEALLOCATE (transmission_aux%tau_ref_surf_ac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_ref_surf_ac")
    DEALLOCATE (transmission_aux%tau_ref_surf_acsun, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_ref_surf_acsun")
    DEALLOCATE (transmission_aux%od_frac_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%od_frac_t")
    DEALLOCATE (transmission_aux%odsun_frac_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%odsun_frac_t")
    DEALLOCATE (transmission_aux%tau_surf_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_surf_t")
    DEALLOCATE (transmission_aux%tausun_surf_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tausun_surf_t")
    DEALLOCATE (transmission_aux%tau_ref_surf_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tau_ref_surf_t")
    DEALLOCATE (transmission_aux%tausun_ref_surf_t, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%tausun_ref_surf_t")
    DEALLOCATE (transmission_aux%refl_norm, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux% refl_norm")
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    CALL rttov_init_transmission_aux(transmission_aux)
  ENDIF
  IF (asw .EQ. 0) THEN
    NULLIFY (transmission_aux%fac)
    NULLIFY (transmission_aux%surf_fac)
    NULLIFY (transmission_aux%tau_surf)
    NULLIFY (transmission_aux%tau_surf_r)
    NULLIFY (transmission_aux%tau_level_r)
    NULLIFY (transmission_aux%tau_level)
    NULLIFY (transmission_aux%od_level)
    NULLIFY (transmission_aux%od_singlelayer)
    NULLIFY (transmission_aux%od_singlelayer_r)
    NULLIFY (transmission_aux%odsun_singlelayer)
    NULLIFY (transmission_aux%tausun_surf)
    NULLIFY (transmission_aux%tausun_level)
    NULLIFY (transmission_aux%od_sfrac)
    NULLIFY (transmission_aux%od_sfrac_r)
    NULLIFY (transmission_aux%odsun_sfrac)
    NULLIFY (transmission_aux%od_frac_ac)
    NULLIFY (transmission_aux%odsun_frac_ac)
    NULLIFY (transmission_aux%tau_surf_ac)
    NULLIFY (transmission_aux%tau_surf_acsun)
    NULLIFY (transmission_aux%tau_ref_surf_ac)
    NULLIFY (transmission_aux%tau_ref_surf_acsun)
    NULLIFY (transmission_aux%od_frac_t)
    NULLIFY (transmission_aux%odsun_frac_t)
    NULLIFY (transmission_aux%tau_surf_t)
    NULLIFY (transmission_aux%tausun_surf_t)
    NULLIFY (transmission_aux%tau_ref_surf_t)
    NULLIFY (transmission_aux%tausun_ref_surf_t)
    NULLIFY (transmission_aux%refl_norm)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_transmission_aux
