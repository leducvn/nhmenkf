Interface
Subroutine rttov_q2v (&
       & h2o_unit,  &! in
       & h2o,       &! in
       & gaz_id,    &! in
       & q_gaz,     &! in
       & v_gaz     ) ! inout
  Use parkind1, Only : jpim     ,jprb
USE PARKIND1  ,ONLY : JPRB
  Implicit None
  Integer(Kind=jpim) , Intent (in) :: h2o_unit ! Water vapour input unit
  Real(Kind=jprb)    , Intent (in) :: h2o      ! Water Vapour content in unit h2o_unit
  Integer(Kind=jpim) , Intent (in) :: gaz_id   ! Gaz identification number
  Real(Kind=jprb)    , Intent (in)   :: q_gaz  ! specific concentration for gaz (kg/kg)
  Real(Kind=jprb)    , Intent (inout):: v_gaz  ! volume mixing ratio for gaz (ppmv)
End Subroutine
End Interface
