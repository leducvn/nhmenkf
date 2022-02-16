Interface
SUBROUTINE rttov_copy_opdp_path(opdp_path1, opdp_path2)
  USE rttov_types, ONLY : opdp_path_Type
  IMPLICIT NONE
  TYPE(opdp_path_Type), INTENT(INOUT) :: opdp_path1
  TYPE(opdp_path_Type), INTENT(IN)    :: opdp_path2
End Subroutine
End Interface
