!
SUBROUTINE rttov_deallocate_atlas( &
                              coef)
  ! Description:
  !   Deallocate any arrays used by emissivity atlases.
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
  USE parkind1, ONLY : jpim
  
  USE rttov_types, ONLY : &
        & rttov_coef

  USE rttov_const, ONLY :      &
        & sensor_id_mw,        &
        & sensor_id_po
        
  USE mod_iratlas, ONLY :             &
        & rttov_uwiremis_close_atlas, &
        & uw_ir_atlas_version => ir_atlas_version
  
  USE mod_mwatlas, ONLY :      &
        & rttov_closemw_atlas, &
        & telsem_mw_atlas_version => mw_atlas_version 
        
  USE mod_cnrm_mw_atlas, ONLY :     &
        & cnrm_mw_atlas_version => mw_atlas_version 
        
  USE mod_rttov_emis_atlas, ONLY : &
        & ir_atlas_version,       &
        & mw_atlas_version,       & 
        & ir_atlas_init,          &
        & mw_atlas_init
        
  IMPLICIT NONE

  TYPE(rttov_coef),  INTENT(IN) :: coef ! RTTOV instrument coefficients
  
!INTF_END

  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
  
    ! MW atlas
    
    IF (mw_atlas_init) THEN
      IF (mw_atlas_version == telsem_mw_atlas_version) THEN
        CALL rttov_closemw_atlas
        mw_atlas_init = .FALSE.
      END IF
      
      IF (mw_atlas_version == cnrm_mw_atlas_version) THEN
        ! nothing to do for CNRM MW atlas
        mw_atlas_init = .FALSE.
      END IF
    END IF
        
  ELSE    
  
    ! IR atlas
    
    IF (ir_atlas_init) THEN
      CALL rttov_uwiremis_close_atlas
      ir_atlas_init = .FALSE.
    END IF
    
  END IF
  
END SUBROUTINE rttov_deallocate_atlas
