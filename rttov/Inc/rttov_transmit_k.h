Interface
SUBROUTINE rttov_transmit_k( &
            & addaerosl,                      &
            & addclouds,                      &
            & nlayers,                        &
            & chanprof,                       &
            & aux,                            &
            & aux_k,                          &
            & coef,                           &
            & ircld,                          &
            & opdp_path,                      &
            & opdp_path_k,                    &
            & od_level,                       &
            & transmission,                   &
            & transmission_k,                 &
            & transmission_aux,               &
            & transmission_aux_k,             &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & tau_ref,                        &
            & tau_ref_surf,                   &
            & tau_surf,                       &
            & tau_level)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type,          &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path                                   ! path  optical depths
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_k
  TYPE(transmission_Type         ), INTENT(IN)    :: transmission
  TYPE(transmission_Type         ), INTENT(INOUT) :: transmission_k
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_k
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref_surf(size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_level    (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_surf    (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_level   (nlayers + 1   , size(chanprof))
End Subroutine
End Interface
