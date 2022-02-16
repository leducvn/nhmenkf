Interface
SUBROUTINE rttov_transmit( &
            & addaerosl,                    &
            & addclouds,                    &
            & nlayers,                      &
            & chanprof,                     &
            & aux,                          &
            & coef,                         &
            & ircld,                        &
            & opdp_path,                    &
            & od_level,                     &
            & transmission,                 &
            & transmission_aux,             &
            & transmission_scatt_ir_stream, &
            & tau_ref,                      &
            & tau_ref_surf,                 &
            & tau_surf,                     &
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
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers                                     ! Number of pressure levels
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)                                 ! Channel indices
  TYPE(rttov_coef                ), INTENT(IN)    :: coef                                        ! Coefficients
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(transmission_Type         ), INTENT(INOUT) :: transmission                                ! Transmittances and single-layer od
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux                            ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux                                         ! auxillary profiles informations
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_ref     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_ref_surf(size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(OUT)   :: od_level    (nlayers + 1   , size(chanprof))! sat to level optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_level   (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_surf    (size(chanprof)                )! sat to surfacetransmittance
End Subroutine
End Interface
