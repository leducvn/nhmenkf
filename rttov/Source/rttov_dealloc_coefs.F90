!
SUBROUTINE rttov_dealloc_coefs(  ERR, coefs )
! Description:
!   Deallocate coefficients structure
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
! Imported Parameters:
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coefs
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_coefs         ), INTENT(INOUT)            :: coefs        ! coefficients
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_dealloc_optpar_ir.h"
#include "rttov_dealloc_coef_scatt_ir.h"
#include "rttov_dealloc_coef_pccomp.h"
#include "rttov_dealloc_coef.h"
! Local Scalars:
!- End of header --------------------------------------------------------
  TRY

  IF (Associated (coefs%optp%optpaer)) THEN
    CALL rttov_dealloc_optpar_ir(ERR, coefs%optp)
    THROW(err.ne.0)
  End If

  IF( Associated(coefs%coef_scatt_ir%fmv_aer_rh)) Then
    CALL rttov_dealloc_coef_scatt_ir (ERR, coefs%coef_scatt_ir)
    THROW(err.ne.0)
  End If

  IF( ASSOCIATED(coefs%coef_pccomp%pcreg) ) Then
    CALL rttov_dealloc_coef_pccomp(ERR, coefs%coef_pccomp)
    THROW(err.ne.0)
  End If

  CALL rttov_dealloc_coef(  ERR, coefs%coef )
  THROW(err.ne.0)

  CATCH
END SUBROUTINE 
