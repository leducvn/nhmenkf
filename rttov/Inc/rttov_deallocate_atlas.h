Interface
SUBROUTINE rttov_deallocate_atlas( &
                              coef)
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
End Subroutine
End Interface
