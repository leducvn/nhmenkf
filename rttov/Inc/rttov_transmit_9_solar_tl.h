Interface
SUBROUTINE rttov_transmit_9_solar_tl( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & chanprof,                        &
            & profiles,                        &
            & sun,                             &
            & aux,                             &
            & aux_tl,                          &
            & coef,                            &
            & raytracing,                      &
            & raytracing_tl,                   &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_tl,                    &
            & odsun_level,                     &
            & odsun_singlelayer,               &
            & od_frac,                         &
            & transmission_aux,                &
            & transmission_aux_tl,             &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & tausun_ref,                      &
            & tausun_ref_surf,                 &
            & tausun_level,                    &
            & tausun_surf)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & profile_Type,               &
       & ircld_type,                 &
       & raytracing_type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  TYPE(rttov_chanprof)            , INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type  )            , INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(IN)    :: aux_tl
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_tl
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_tl
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing_tl
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_tl
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_level      (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: odsun_singlelayer(nlayers       , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref_surf  (size(chanprof)                )! sat to surface transmittance
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_ref       (nlayers + 1   , size(chanprof))! sat to layer transmittance
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_level     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tausun_surf      (size(chanprof)                )
End Subroutine
End Interface
