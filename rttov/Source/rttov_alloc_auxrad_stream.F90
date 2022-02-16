SUBROUTINE rttov_alloc_auxrad_stream( &
            & ERR,           &
            & auxrad_stream, &
            & opts,          &
            & nstreams,      &
            & nlayers,       &
            & nchannels,     &
            & asw,           &
            & init)
! Description:
!   Allocates/deallocates the auxiliary radiance structure 
!   for streams
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

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_options, radiance_aux
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
!INTF_OFF
#include "throw.h"
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim) , INTENT(OUT)             :: ERR
  TYPE(radiance_aux ), INTENT(INOUT)           :: auxrad_stream
  TYPE(rttov_options), INTENT(IN)              :: opts
  INTEGER(KIND=jpim) , INTENT(IN)              :: nstreams
  INTEGER(KIND=jpim) , INTENT(IN)              :: nlayers
  INTEGER(KIND=jpim) , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim) , INTENT(IN)              :: asw
  LOGICAL(KIND=jplm) , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_auxrad_stream.h"
  LOGICAL(KIND=jplm) :: init1
  TRY
  init1 = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    ALLOCATE (auxrad_stream%up(nlayers, 0:nstreams, nchannels), auxrad_stream%down(nlayers, 0:nstreams, nchannels),      &
      & auxrad_stream%down_ref(nlayers, 0:nstreams, nchannels), auxrad_stream%cloudy(0:nstreams, nchannels), STAT = ERR)
    ALLOCATE (auxrad_stream%meanrad_up(0:nstreams, nchannels), auxrad_stream%meanrad_down(0:nstreams, nchannels), STAT = ERR)

    THROWM( ERR .NE. 0 , "allocation of auxrad_stream")
    IF (opts%addsolar .AND. (opts%addaerosl .OR. opts%addclouds)) THEN
      ALLOCATE (auxrad_stream%FAC1_1(0:nstreams, nchannels), auxrad_stream%FAC2_1(0:nstreams, nchannels),              &
                auxrad_stream%FAC3_1(0:nstreams, nchannels), auxrad_stream%FAC4_1(0:nstreams, nchannels),              &
                auxrad_stream%FAC5_1(0:nstreams, nchannels), auxrad_stream%FAC6_1(0:nstreams, nchannels),              &
                auxrad_stream%FAC7_1(0:nstreams, nchannels), &
                auxrad_stream%FAC1_2(nlayers, 0:nstreams, nchannels), auxrad_stream%FAC2_2(nlayers, 0:nstreams, nchannels),  &
                auxrad_stream%FAC3_2(nlayers, 0:nstreams, nchannels), auxrad_stream%FAC4_2(nlayers, 0:nstreams, nchannels),  &
                auxrad_stream%FAC5_2(nlayers, 0:nstreams, nchannels), auxrad_stream%FAC6_2(nlayers, 0:nstreams, nchannels),  &
                auxrad_stream%FAC7_2(nlayers, 0:nstreams, nchannels), &
                auxrad_stream%FAC1_3(0:nstreams, nchannels), auxrad_stream%FAC2_3(0:nstreams, nchannels),                    &
                auxrad_stream%FAC3_3(0:nstreams, nchannels), auxrad_stream%FAC4_3(0:nstreams, nchannels),                    &
                auxrad_stream%FAC5_3(0:nstreams, nchannels), auxrad_stream%FAC6_3(0:nstreams, nchannels),                    &
                auxrad_stream%FAC7_3(0:nstreams, nchannels), STAT = ERR)
    ELSE
      NULLIFY (auxrad_stream%FAC1_1, auxrad_stream%FAC2_1, auxrad_stream%FAC3_1, auxrad_stream%FAC4_1,                   &
        & auxrad_stream%FAC5_1, auxrad_stream%FAC6_1, auxrad_stream%FAC7_1, auxrad_stream%FAC1_2, auxrad_stream%FAC2_2,  &
        & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, auxrad_stream%FAC6_2, auxrad_stream%FAC7_2,  &
        & auxrad_stream%FAC1_3, auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, auxrad_stream%FAC5_3,  &
        & auxrad_stream%FAC6_3, auxrad_stream%FAC7_3)
    ENDIF
    IF (init1) THEN
      CALL rttov_init_auxrad_stream(auxrad_stream)
    ENDIF
    THROWM( ERR .NE. 0 , "allocation of auxrad_stream%FAC")
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (auxrad_stream%up, auxrad_stream%down, auxrad_stream%down_ref, auxrad_stream%cloudy, STAT = ERR)
    DEALLOCATE (auxrad_stream%meanrad_up, auxrad_stream%meanrad_down, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of auxrad_stream")
    IF (opts%addsolar) THEN
      IF (opts%addaerosl .OR. opts%addclouds) THEN
        DEALLOCATE (auxrad_stream%FAC1_1, auxrad_stream%FAC2_1, auxrad_stream%FAC3_1, auxrad_stream%FAC4_1,                &
          & auxrad_stream%FAC5_1, auxrad_stream%FAC6_1, auxrad_stream%FAC7_1, auxrad_stream%FAC1_2, auxrad_stream%FAC2_2,  &
          & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, auxrad_stream%FAC6_2, auxrad_stream%FAC7_2,  &
          & auxrad_stream%FAC1_3, auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, auxrad_stream%FAC5_3,  &
          & auxrad_stream%FAC6_3, auxrad_stream%FAC7_3, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of auxrad_stream%FAC")
      ENDIF
    ENDIF
    NULLIFY (auxrad_stream%up, auxrad_stream%down, auxrad_stream%down_ref, auxrad_stream%cloudy)
    NULLIFY (auxrad_stream%meanrad_up, auxrad_stream%meanrad_down)
    IF (opts%addsolar) THEN
      IF (opts%addaerosl .OR. opts%addclouds) THEN
        NULLIFY (auxrad_stream%FAC1_1, auxrad_stream%FAC2_1, auxrad_stream%FAC3_1, auxrad_stream%FAC4_1,                   &
          & auxrad_stream%FAC5_1, auxrad_stream%FAC6_1, auxrad_stream%FAC7_1, auxrad_stream%FAC1_2, auxrad_stream%FAC2_2,  &
          & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, auxrad_stream%FAC6_2, auxrad_stream%FAC7_2,  &
          & auxrad_stream%FAC1_3, auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, auxrad_stream%FAC5_3,  &
          & auxrad_stream%FAC6_3, auxrad_stream%FAC7_3)
      ENDIF
    ENDIF
  ENDIF
  CATCH
END SUBROUTINE 
