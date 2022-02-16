Interface
SUBROUTINE rttov_opdpscattir( &
            & nlayers,                      &
            & chanprof,                     &
            & opts,                         &
            & aux,                          &
            & profiles,                     &
            & sun,                          &
            & coef,                         &
            & coef_scatt_ir,                &
            & raytracing,                   &
            & transmission_scatt_ir,        &
            & transmission_scatt_ir_stream, &
            & optp,                         &
            & ircld)
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
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(profile_aux               ), INTENT(INOUT) :: aux
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type                ), INTENT(INOUT) :: ircld
  TYPE(raytracing_type           ), INTENT(INOUT) :: raytracing
End Subroutine
End Interface
