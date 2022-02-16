Interface
SUBROUTINE rttov_alloc_trans_scatt_ir( &
            & ERR,                   &
            & transmission_scatt_ir, &
            & nchannels,             &
            & ncldtyp,               &
            & nlayers,               &
            & asw,                   &
            & nstreams,              &
            & stream,                &
            & init)
  USE rttov_types, ONLY : transmission_scatt_ir_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(OUT)             :: ERR                  ! return code
  INTEGER(KIND=jpim)              , INTENT(IN)              :: nlayers              ! number of levels
  INTEGER(KIND=jpim)              , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)              , INTENT(IN)              :: ncldtyp
  TYPE(transmission_scatt_ir_Type), INTENT(INOUT)           :: transmission_scatt_ir
  INTEGER(KIND=jpim)              , INTENT(IN)              :: asw                  ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim)              , INTENT(IN)   , OPTIONAL :: nstreams
  LOGICAL(KIND=jplm)              , INTENT(IN)   , OPTIONAL :: stream
  LOGICAL(KIND=jplm)              , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
