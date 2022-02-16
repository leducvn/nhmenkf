Interface
SUBROUTINE rttov_get_emis(  &
              & err,        & ! out
              & chanprof,   & ! in
              & profiles,   & ! in
              & coef,       & ! in
              & resolution, & ! in, optional (MW atlas only)
              & emissivity, & ! out
              & emis_std,   & ! out, optional (not CNRM atlas)
              & emis_cov,   & ! out, optional (MW atlas only)
              & emis_flag,  & ! out, optional (IR atlas only)
              & pbats_veg)    ! out, optional (MW CNRM atlas only)
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : &
          rttov_chanprof, &
          profile_type,   &
          rttov_coef
  USE rttov_const, ONLY :      &
          deg2rad,             &
          sensor_id_mw,        &
          sensor_id_po,        &
          surftype_land,       &
          surftype_seaice,     &
          pol_v,               &
          pol_h,               &
          default_err_unit,    &
          errorstatus_success, &
          errorstatus_fatal
  USE mod_iratlas, ONLY :  &
          rttov_uwiremis,  &
          uw_ir_atlas_version => ir_atlas_version
  USE mod_mwatlas, ONLY :       &
          atlas,                &
          emis_interp_ind_mult, &
          emis_interp_int_mult, &
          telsem_mw_atlas_version =>  mw_atlas_version
  USE mod_cnrm_mw_atlas, ONLY :     &
          rttov_cnrmmwemis,         &
          cnrm_mw_atlas_version => mw_atlas_version
  USE mod_rttov_emis_atlas, ONLY : &
          ir_atlas_version,        &
          mw_atlas_version,        & 
          ir_atlas_init,           &
          mw_atlas_init
  USE yomhook, ONLY : &
          LHOOK,      &
          DR_HOOK
  IMPLICIT NONE
  INTEGER(KIND=jpim),   INTENT(OUT)           :: err
  TYPE(rttov_chanprof), INTENT(IN)            :: chanprof(:)
  TYPE(profile_type),   INTENT(IN)            :: profiles(:)
  TYPE(rttov_coef),     INTENT(IN)            :: coef
  REAL(KIND=jprb),      INTENT(IN),  OPTIONAL :: resolution
  REAL(KIND=jprb),      INTENT(OUT)           :: emissivity(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: emis_std(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: emis_cov(:,:,:)
  INTEGER(KIND=jpim),   INTENT(OUT), OPTIONAL :: emis_flag(size(chanprof))
  REAL(KIND=jprb),      INTENT(OUT), OPTIONAL :: pbats_veg(size(chanprof))
End Subroutine
End Interface
