Interface
SUBROUTINE rttov_atlas_setup( &
             &  err,          &! out
             &  imonth,       &! in
             &  coef,         &! in
             &  path,         &! in, optional
             &  ir_atlas_ver, &! in, optional
             &  mw_atlas_ver)  ! in, optional
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
End Subroutine
End Interface
