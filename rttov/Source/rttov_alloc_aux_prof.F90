SUBROUTINE rttov_alloc_aux_prof( &
            & err,       &
            & nprofiles, &
            & nlevels,   &
            & id_sensor, &
            & aux_prof,  &
            & opts,      &
            & asw,       &
            & init)
! Description:
!   Allocates/deallocates the auxiliary profile structure
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
  USE rttov_types, ONLY : rttov_options, profile_aux
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po
  USE parkind1, ONLY : jprb
!INTF_ON
  INTEGER(KIND=jpim) , INTENT(OUT)               :: err
  INTEGER(KIND=jpim) , INTENT(IN)                :: nprofiles
  INTEGER(KIND=jpim) , INTENT(IN)                :: nlevels
  INTEGER(KIND=jpim) , INTENT(IN)                :: id_sensor
  TYPE(profile_aux  ), INTENT(INOUT)             :: aux_prof
  TYPE(rttov_options), INTENT(IN)                :: opts
  INTEGER(KIND=jpim) , INTENT(IN)                :: asw
  LOGICAL(KIND=jplm) , OPTIONAL     , INTENT(IN) :: init
!INTF_END

#include "rttov_errorreport.h"
#include "rttov_init_aux_prof.h"

  LOGICAL(KIND=jplm) :: init1
  TRY
  init1 = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    NULLIFY (aux_prof%s)
    NULLIFY (aux_prof%debye_prof)
    NULLIFY (aux_prof%dg)
    NULLIFY (aux_prof%fac1_dg)
    NULLIFY (aux_prof%fac2_dg)
    NULLIFY (aux_prof%fac3_dg)
    NULLIFY (aux_prof%iaernum)
    NULLIFY (aux_prof%iaertyp)
    NULLIFY (aux_prof%relhum)
    NULLIFY (aux_prof%relhumref)
    ALLOCATE (aux_prof%s(nprofiles), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of aux_prof%s")
    IF (id_sensor == sensor_id_mw .OR. id_sensor == sensor_id_po) THEN
      ALLOCATE (aux_prof%debye_prof(5, nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%debye_prof")
    ENDIF
    IF (opts%addclouds) THEN
      ALLOCATE (aux_prof%dg(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%dg")
      ALLOCATE (aux_prof%fac1_dg(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac1_dg")
      ALLOCATE (aux_prof%fac2_dg(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac2_dg")
      ALLOCATE (aux_prof%fac3_dg(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac3_dg")
    ENDIF
    IF (opts%addaerosl) THEN
      ALLOCATE (aux_prof%iaernum(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%iaernum")
      ALLOCATE (aux_prof%iaertyp(11, nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%iaertyp")
      ALLOCATE (aux_prof%relhum(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%relhum")
      ALLOCATE (aux_prof%relhumref(nlevels, nprofiles), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%relhumref")
    ENDIF
    IF (init1) THEN
      CALL rttov_init_aux_prof(aux_prof)
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (aux_prof%s, STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of aux_prof%s")
    IF (id_sensor == sensor_id_mw .OR. id_sensor == sensor_id_po) THEN
      DEALLOCATE (aux_prof%debye_prof, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%debye_prof")
    ENDIF
    IF (opts%addclouds) THEN
      DEALLOCATE (aux_prof%dg, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%dg")
      DEALLOCATE (aux_prof%fac1_dg, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac1_dg")
      DEALLOCATE (aux_prof%fac2_dg, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac2_dg")
      DEALLOCATE (aux_prof%fac3_dg, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%fac3_dg")
    ENDIF
    IF (opts%addaerosl) THEN
      DEALLOCATE (aux_prof%iaernum, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%iaernum")
      DEALLOCATE (aux_prof%iaertyp, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%iaertyp")
      DEALLOCATE (aux_prof%relhum, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%relhum")
      DEALLOCATE (aux_prof%relhumref, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of aux_prof%relhumref")
    ENDIF
    NULLIFY (aux_prof%s)
    NULLIFY (aux_prof%debye_prof)
    NULLIFY (aux_prof%dg)
    NULLIFY (aux_prof%fac1_dg)
    NULLIFY (aux_prof%fac2_dg)
    NULLIFY (aux_prof%fac3_dg)
    NULLIFY (aux_prof%iaernum)
    NULLIFY (aux_prof%iaertyp)
    NULLIFY (aux_prof%relhum)
    NULLIFY (aux_prof%relhumref)
  ENDIF
  CATCH
END SUBROUTINE 
