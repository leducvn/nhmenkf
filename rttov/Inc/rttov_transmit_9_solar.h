Interface
SUBROUTINE rttov_transmit_9_solar( &
            & addaerosl,                    &
            & addclouds,                    &
            & nlayers,                      &
            & chanprof,                     &
            & profiles,                     &
            & sun,                          &
            & aux,                          &
            & coef,                         &
            & raytracing,                   &
            & ircld,                        &
            & opdp_path,                    &
            & odsun_level,                  &
            & odsun_singlelayer,            &
            & od_frac,                      &
            & transmission_aux,             &
            & transmission_scatt_ir_stream, &
            & tausun_ref,                   &
            & tausun_ref_surf,              &
            & tausun_surf,                  &
            & tausun_level)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type,                 &
       & raytracing_type,            &
       & profile_Type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers                                          ! Number of pressure levels
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)                                      ! Channel indices
  TYPE(profile_Type              ), INTENT(IN)    :: profiles(:)                                      ! Atmospheric profiles
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef                                             ! Coefficients
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux                                 ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux                                              ! auxillary profiles informations
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  REAL(KIND=jprb)                 , INTENT(OUT)   :: odsun_level      (nlayers + 1   , size(chanprof))! sat to level optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: odsun_singlelayer(nlayers       , size(chanprof))! single layer optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_ref_surf  (size(chanprof)                )! sat to surface transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_ref       (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_level     (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_surf      (size(chanprof)                )! sat to surface transmittance
End Subroutine
End Interface
