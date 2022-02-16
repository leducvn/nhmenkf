Interface
SUBROUTINE rttov_opdpscattir_k( &
            & nlayers,                        &
            & chanprof,                       &
            & opts,                           &
            & aux,                            &
            & aux_k,                          &
            & profiles,                       &
            & profiles_k,                     &
            & sun,                            &
            & coef,                           &
            & coef_scatt_ir,                  &
            & raytracing,                     &
            & raytracing_k,                   &
            & transmission_scatt_ir,          &
            & transmission_scatt_ir_k,        &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & optp,                           &
            & ircld,                          &
            & ircld_k)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_options,              &
       & rttov_coef,                 &
       & profile_Type,               &
       & raytracing_type,            &
       & transmission_scatt_ir_type, &
       & rttov_coef_scatt_ir,        &
       & rttov_optpar_ir,            &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof  (:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_k(size(chanprof))
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type),                 INTENT(IN)    :: ircld
  TYPE(ircld_type),                 INTENT(INOUT) :: ircld_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
End Subroutine
End Interface
