!
SUBROUTINE rttov_alloc_prof( &
            & ERR,      &
            & nprof,    &
            & profiles, &
            & nlevels,  &
            & opts,     &
            & asw,      &
            & coefs,    &
            & init,     &
            & blob)
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
!  1.2       26/09/2008  Disable vectorization of the pointer association
!                        loop; this is for the NEC under -Chopt
!                        P. Marguinaud
!  1.3       02/12/2009  Allocate aerosols,cloud and cfrac variables on number of layers.
!                        Marco Matricardi
!  1.4       05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.5       14/07/2010  Replace iaer and icld arguments with optional coefs argument to
!                        obtain size of clouds and aerosol arrays from coef file (J Hocking)
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
  USE rttov_types, ONLY : rttov_options, profile_Type, rttov_coefs, blob_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim) , INTENT(OUT)             :: ERR            ! return code
  INTEGER(KIND=jpim) , INTENT(IN)              :: nprof          ! number of profiles
  INTEGER(KIND=jpim) , INTENT(IN)              :: nlevels        ! number of levels
  TYPE(profile_Type ), INTENT(INOUT)           :: profiles(nprof)! profiles
  TYPE(rttov_options), INTENT(IN)              :: opts
  INTEGER(KIND=jpim) , INTENT(IN)              :: asw            ! 1=allocate, 0=deallocate
  TYPE(rttov_coefs)  , INTENT(IN)   , OPTIONAL :: coefs
  LOGICAL(KIND=jplm) , INTENT(IN)   , OPTIONAL :: init
  TYPE(blob_type)    , INTENT(INOUT), OPTIONAL :: blob
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_prof.h"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: nlayers
  INTEGER(KIND=jpim) :: natm
  INTEGER(KIND=jpim) :: iatm
  INTEGER(KIND=jpim) :: ncldtyp
  INTEGER(KIND=jpim) :: naertyp
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (Present(init)) init1 = init
  IF (opts%addclouds .AND. .NOT. Present(coefs) .AND. (asw .EQ. 1_jpim)) THEN
    err = errorstatus_fatal
    THROWM( ERR .NE. 0 , "Dummy argument coefs is required")
  ENDIF
  IF (opts%addaerosl .AND. .NOT. Present(coefs) .AND. (asw .EQ. 1_jpim)) THEN
    err = errorstatus_fatal
    THROWM( ERR .NE. 0 , "Dummy argument coefs is required")
  ENDIF
!  Allocate section
  IF (opts%addclouds .AND. (asw .EQ. 1_jpim)) &
      & ncldtyp = max(1_jpim, coefs % coef_scatt_ir % fmv_wcl_comp+1_jpim)
  IF (opts%addaerosl .AND. (asw .EQ. 1_jpim)) &
      & naertyp = max(1_jpim, coefs % coef_scatt_ir % fmv_aer_comp)
  IF (asw .EQ. 1) THEN
    DO j = 1, nprof
      NULLIFY (profiles(j)%p)
      NULLIFY (profiles(j)%t)
      NULLIFY (profiles(j)%q)
      NULLIFY (profiles(j)%o3)
      NULLIFY (profiles(j)%co2)
      NULLIFY (profiles(j)%n2o)
      NULLIFY (profiles(j)%co)
      NULLIFY (profiles(j)%ch4)
      NULLIFY (profiles(j)%clw)
      NULLIFY (profiles(j)%aerosols)
      NULLIFY (profiles(j)%cloud)
      NULLIFY (profiles(j)%cfrac)
      NULLIFY (profiles(j)%icede)
    ENDDO
  ENDIF
  IF (Present(blob) .AND. (asw .EQ. 1)) THEN
    NULLIFY (blob%rlvp)
    NULLIFY (blob%ralp)
    NULLIFY (blob%rclp)
    NULLIFY (blob%rtlp)
    natm = 3_jpim! p, t, q
    IF (opts%ozone_data) natm = natm + 1
    IF (opts%co2_data  ) natm = natm + 1
    IF (opts%n2o_data  ) natm = natm + 1
    IF (opts%co_data   ) natm = natm + 1
    IF (opts%ch4_data  ) natm = natm + 1
    IF (opts%clw_data  ) natm = natm + 1
    IF (opts%addclouds ) natm = natm + 1  ! icede
    ALLOCATE (blob%rlvp(nlevels, natm, nprof), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of blob%rlvp")
    IF (opts%addaerosl) ALLOCATE (blob%ralp(naertyp, nlayers, nprof), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of blob%ralp")
    IF (opts%addclouds) THEN
      ALLOCATE (blob%rclp(ncldtyp, nlayers, nprof), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of blob%rclp")
      ALLOCATE (blob%rtlp(ncldtyp, nlayers, nprof), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of blob%rtlp")
    ENDIF
!CDIR NOVECTOR
    DO j = 1, nprof
      profiles(j)%nlevels = nlevels
      profiles(j)%nlayers = nlevels - 1
      profiles(j)%p => blob%rlvp(:, 1, j)
      profiles(j)%t => blob%rlvp(:, 2, j)
      profiles(j)%q => blob%rlvp(:, 3, j)
      iatm = 3_jpim
      IF (opts%co2_data) THEN
        iatm = iatm + 1
        profiles(j)%co2 => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%ozone_data) THEN
        iatm = iatm + 1
        profiles(j)%o3 => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%n2o_data) THEN
        iatm = iatm + 1
        profiles(j)%n2o => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%co_data) THEN
        iatm = iatm + 1
        profiles(j)%co => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%ch4_data) THEN
        iatm = iatm + 1
        profiles(j)%ch4 => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%clw_data) THEN
        iatm = iatm + 1
        profiles(j)%clw => blob%rlvp(:, iatm, j)
      ENDIF
      IF (opts%addaerosl) profiles(j)%aerosols => blob%ralp(:, :, j)
      IF (opts%addclouds) THEN
        profiles(j)%cloud => blob%rtlp(:, :, j)
        profiles(j)%cfrac => blob%rclp(:, :, j)
        iatm = iatm + 1
        profiles(j)%icede => blob%rlvp(:, iatm, j)
      ENDIF
    ENDDO
  ENDIF
  IF (Present(blob) .AND. (asw .EQ. 0)) THEN
    IF (Associated(blob%rlvp)) DEALLOCATE (blob%rlvp, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (Associated(blob%ralp)) DEALLOCATE (blob%ralp, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (Associated(blob%rclp)) DEALLOCATE (blob%rclp, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (Associated(blob%rtlp)) DEALLOCATE (blob%rtlp, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    NULLIFY (blob%rlvp)
    NULLIFY (blob%ralp)
    NULLIFY (blob%rclp)
    NULLIFY (blob%rtlp)
  ENDIF
  IF ((.NOT. Present(blob)) .AND. (asw .EQ. 1)) THEN
    DO j = 1, nprof
      profiles(j)%nlevels = nlevels
      profiles(j)%nlayers = nlevels - 1
      ALLOCATE (profiles(j)%p(nlevels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of profiles%p")
      ALLOCATE (profiles(j)%t(nlevels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of profiles%t")
      ALLOCATE (profiles(j)%q(nlevels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of profiles%q")
      IF (opts%co2_data) THEN
        ALLOCATE (profiles(j)%co2(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%co2")
      ENDIF
      IF (opts%ozone_data) THEN
        ALLOCATE (profiles(j)%o3(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%o3 ")
      ENDIF
      IF (opts%n2o_data) THEN
        ALLOCATE (profiles(j)%n2o(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%n2o")
      ENDIF
      IF (opts%co_data) THEN
        ALLOCATE (profiles(j)%co(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%co")
      ENDIF
      IF (opts%ch4_data) THEN
        ALLOCATE (profiles(j)%ch4(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%ch4")
      ENDIF
      IF (opts%clw_data) THEN
        ALLOCATE (profiles(j)%clw(nlevels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%clw")
      ENDIF
      IF (opts%addaerosl) ALLOCATE (profiles(j)%aerosols(naertyp, nlayers), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of profiles%aerosols")
      IF (opts%addclouds) THEN
        ALLOCATE (profiles(j)%cloud(ncldtyp, nlayers), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%cloud")
        ALLOCATE (profiles(j)%cfrac(ncldtyp, nlayers), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%cfrac")
        ALLOCATE (profiles(j)%icede(nlayers), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of profiles%icede")
      ENDIF
    ENDDO
  ENDIF
  IF ((.NOT. Present(blob)) .AND. (asw .EQ. 0)) THEN
    DO j = 1, nprof
! deallocate model profiles atmospheric arrays
      DEALLOCATE (profiles(j)%p, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of profiles%p")
      DEALLOCATE (profiles(j)%t, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of profiles%t")
      DEALLOCATE (profiles(j)%q, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of profiles%q")
      IF (associated(profiles(j)%o3)) THEN
        DEALLOCATE (profiles(j)%o3, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%o3")
      ENDIF
      IF (associated(profiles(j)%co2)) THEN
        DEALLOCATE (profiles(j)%co2, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%co2")
      ENDIF
      IF (associated(profiles(j)%n2o)) THEN
        DEALLOCATE (profiles(j)%n2o, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%n2o")
      ENDIF
      IF (associated(profiles(j)%co)) THEN
        DEALLOCATE (profiles(j)%co, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%co")
      ENDIF
      IF (associated(profiles(j)%ch4)) THEN
        DEALLOCATE (profiles(j)%ch4, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%ch4")
      ENDIF
      IF (associated(profiles(j)%clw)) THEN
        DEALLOCATE (profiles(j)%clw, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%clw")
      ENDIF
      IF (associated(profiles(j)%aerosols)) THEN
        DEALLOCATE (profiles(j)%aerosols, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%aerosols")
      ENDIF
      IF (associated(profiles(j)%cloud)) THEN
        DEALLOCATE (profiles(j)%cloud, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%cloud")
      ENDIF
      IF (associated(profiles(j)%cfrac)) THEN
        DEALLOCATE (profiles(j)%cfrac, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%cfrac")
      ENDIF
      IF (associated(profiles(j)%icede)) THEN
        DEALLOCATE (profiles(j)%icede, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of profiles%icede")
      ENDIF
      NULLIFY (profiles(j)%p)
      NULLIFY (profiles(j)%t)
      NULLIFY (profiles(j)%q)
      NULLIFY (profiles(j)%o3)
      NULLIFY (profiles(j)%co2)
      NULLIFY (profiles(j)%n2o)
      NULLIFY (profiles(j)%co)
      NULLIFY (profiles(j)%ch4)
      NULLIFY (profiles(j)%clw)
      NULLIFY (profiles(j)%aerosols)
      NULLIFY (profiles(j)%cloud)
      NULLIFY (profiles(j)%cfrac)
      NULLIFY (profiles(j)%icede)
    ENDDO
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    CALL rttov_init_prof(profiles)
  ENDIF
  IF (asw .EQ. 0) THEN
    DO j = 1, nprof
      NULLIFY (profiles(j)%p)
      NULLIFY (profiles(j)%t)
      NULLIFY (profiles(j)%q)
      NULLIFY (profiles(j)%o3)
      NULLIFY (profiles(j)%co2)
      NULLIFY (profiles(j)%n2o)
      NULLIFY (profiles(j)%co)
      NULLIFY (profiles(j)%ch4)
      NULLIFY (profiles(j)%clw)
      NULLIFY (profiles(j)%aerosols)
      NULLIFY (profiles(j)%cloud)
      NULLIFY (profiles(j)%cfrac)
      NULLIFY (profiles(j)%icede)
    ENDDO
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_prof
