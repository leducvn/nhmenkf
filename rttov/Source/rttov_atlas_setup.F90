!
SUBROUTINE rttov_atlas_setup( &
             &  err,          &! out
             &  imonth,       &! in
             &  coef,         &! in
             &  path,         &! in, optional
             &  ir_atlas_ver, &! in, optional
             &  mw_atlas_ver)  ! in, optional
             
  ! Description:
  !   Sets up the emissivity atlas appropriate to the sensor
  !   type. Data is loaded for the specified month.
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
  !  1.0      02/06/2010  Created (J. Hocking)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
#include "throw.h"

  USE parkind1, ONLY : jpim
  
  USE rttov_types, ONLY : &
        & rttov_coef

  USE rttov_const, ONLY :      &
        & sensor_id_mw,        &
        & sensor_id_po,        &
        & errorstatus_success, &
        & errorstatus_fatal
        
  USE mod_iratlas, ONLY :      &
        & rttov_uwiremis_init, &
        & uw_ir_atlas_version => ir_atlas_version
        
  USE mod_mwatlas, ONLY :     &
        & rttov_readmw_atlas, &
        & telsem_mw_atlas_version => mw_atlas_version 

  USE mod_cnrm_mw_atlas, ONLY :     &
        & rttov_cnrmmwemis_init,    &
        & cnrm_mw_atlas_version => mw_atlas_version 

  USE mod_rttov_emis_atlas, ONLY : &
        & ir_atlas_version,        &
        & mw_atlas_version,        &
        & ir_atlas_init,           &
        & mw_atlas_init
        
  IMPLICIT NONE
  
  INTEGER(KIND=jpim), INTENT(IN)           :: imonth       ! month for which atlas data required
  TYPE(rttov_coef),   INTENT(IN)           :: coef         ! RTTOV instrument coefficients
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path         ! path to atlas data (if not the default)
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: ir_atlas_ver ! version of IR atlas to use
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: mw_atlas_ver ! version of MW atlas to use
  
  INTEGER(KIND=jpim), INTENT(OUT)          :: err          ! output error status
  
!INTF_END

  CHARACTER(LEN=300) :: fpath
  CHARACTER(LEN=300) :: msg
  
!----------------------------------------------------------------------------
  TRY
  
  IF (PRESENT(path)) THEN
    fpath = TRIM(path)//'/'
  ELSE
    fpath = './'
  END IF
  
  ! default values for Atlas versions
  mw_atlas_version = telsem_mw_atlas_version
  ir_atlas_version = uw_ir_atlas_version

  IF (PRESENT(ir_atlas_ver)) ir_atlas_version = ir_atlas_ver
  IF (PRESENT(mw_atlas_ver)) mw_atlas_version = mw_atlas_ver
  
  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
  
    ! MW atlas
    
    IF (mw_atlas_init) THEN
      WARN('MW emissivity atlas already initialised.')
    ELSE
      IF (mw_atlas_version == telsem_mw_atlas_version) THEN
        CALL rttov_readmw_atlas(TRIM(fpath),imonth,err)
        THROWM(err .NE. 0, 'Error initialising TELSEM MW emissivity atlas.')
      ELSEIF (mw_atlas_version == cnrm_mw_atlas_version) THEN
        CALL rttov_cnrmmwemis_init(TRIM(fpath),imonth,coef,err)
        THROWM(err .NE. 0, 'Error initialising CNRM MW emissivity atlas.')
      ELSE
        WRITE(msg,'(a,i5)') 'Unknown MW atlas version: ', mw_atlas_version
        err = errorstatus_fatal
        THROWM(err .NE. 0, msg)
      END IF
      IF (err == 0) mw_atlas_init = .TRUE.
    END IF
  ELSE
    
    ! IR atlas
    
    IF (ir_atlas_init) THEN
      WARN('IR emissivity atlas already initialised.')
    ELSE
      IF (ir_atlas_version == uw_ir_atlas_version) THEN
        CALL rttov_uwiremis_init(TRIM(fpath),imonth,err)
        THROWM(err .NE. 0, 'Error initialising IR emissivity atlas.')
      ELSE
        WRITE(msg,'(a,i5)') 'Unknown IR atlas version: ', ir_atlas_version
        err = errorstatus_fatal
        THROWM(err .NE. 0, msg)
      END IF
      IF (err == 0) ir_atlas_init = .TRUE.
    END IF
  END IF
  
  CATCH
  
END SUBROUTINE rttov_atlas_setup
