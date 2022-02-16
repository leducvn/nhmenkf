!
SUBROUTINE rttov_dealloc_coef (ERR, coef)
! Description:
! de-allocation of a coefficient structure
! The allocation is done by the readcoef subroutine called by the user
! this subroutine should be called once per coef structure when
! all rttov calls are completed.
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       03/05/2004  Add specific RTTOV8 CO2 variable (P. Brunel)
!  1.2       02/06/2004  Update for RTTOV8 coefS (P. Brunel)
!  1.3       01/06/2005  Marco Matricardi (ECMWF):
!               --       N2O,CO and CH4 variables added
!  1.4       26/04/2007  Cloud/aerosol variables added (R Saunders)
!  1.5       11/10/2007  Nullify unused coef pointers P.Marguinaud
!  1.6       23/94/2008  Add some initialisation of scalars (P. Brunel)
!  1.7       06/03/2009  Deafault now coef % IncTop = .true. (P.Rayer)
!  1.8       02/12/2009  Add principal components (M. Matricardi, ECMWF)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
#include "throw.h"
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)   :: ERR          ! return code
  TYPE(rttov_coef         ), INTENT(INOUT) :: coef         ! coefficients
!INTF_END
#include "rttov_nullify_coef.h"
#include "rttov_errorreport.h"
! Local Arrays and Scalars:
!- End of header --------------------------------------------------------

  TRY

  IF (Associated(coef%fmv_gas_id)) DEALLOCATE (coef%fmv_gas_id, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%fmv_gas_pos)) DEALLOCATE (coef%fmv_gas_pos, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%fmv_var)) DEALLOCATE (coef%fmv_var, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%fmv_coe)) DEALLOCATE (coef%fmv_coe, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%fmv_lvl)) DEALLOCATE (coef%fmv_lvl, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_ori_chn)) DEALLOCATE (coef%ff_ori_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_val_chn)) DEALLOCATE (coef%ff_val_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_cwn)) DEALLOCATE (coef%ff_cwn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_bco)) DEALLOCATE (coef%ff_bco, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_bcs)) DEALLOCATE (coef%ff_bcs, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ff_gam)) DEALLOCATE (coef%ff_gam, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%gaz_units)) DEALLOCATE (coef%gaz_units, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%fastem_polar)) DEALLOCATE (coef%fastem_polar, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%ssirem_chn)) DEALLOCATE (coef%ssirem_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ssirem_a0)) DEALLOCATE (coef%ssirem_a0, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ssirem_a1)) DEALLOCATE (coef%ssirem_a1, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ssirem_a2)) DEALLOCATE (coef%ssirem_a2, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ssirem_xzn1)) DEALLOCATE (coef%ssirem_xzn1, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ssirem_xzn2)) DEALLOCATE (coef%ssirem_xzn2, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%ref_prfl_p)) DEALLOCATE (coef%ref_prfl_p, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ref_prfl_t)) DEALLOCATE (coef%ref_prfl_t, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ref_prfl_mr)) DEALLOCATE (coef%ref_prfl_mr, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%lim_prfl_p)) DEALLOCATE (coef%lim_prfl_p, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%lim_prfl_tmax)) DEALLOCATE (coef%lim_prfl_tmax, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%lim_prfl_tmin)) DEALLOCATE (coef%lim_prfl_tmin, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%lim_prfl_gmin)) DEALLOCATE (coef%lim_prfl_gmin, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%lim_prfl_gmax)) DEALLOCATE (coef%lim_prfl_gmax, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%mixedgas)) DEALLOCATE (coef%mixedgas, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%mixedgasint)) DEALLOCATE (coef%mixedgasint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%watervapour)) DEALLOCATE (coef%watervapour, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%watervapourint)) DEALLOCATE (coef%watervapourint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%ozone)) DEALLOCATE (coef%ozone, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ozoneint)) DEALLOCATE (coef%ozoneint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%wvcont)) DEALLOCATE (coef%wvcont, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%wvcontint)) DEALLOCATE (coef%wvcontint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%co2)) DEALLOCATE (coef%co2, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%co2int)) DEALLOCATE (coef%co2int, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%n2o)) DEALLOCATE (coef%n2o, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%n2oint)) DEALLOCATE (coef%n2oint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%co)) DEALLOCATE (coef%co, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%coint)) DEALLOCATE (coef%coint, STAT = ERR)

  THROW( ERR .NE. 0 )



  IF (Associated(coef%ch4)) DEALLOCATE (coef%ch4, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ch4int)) DEALLOCATE (coef%ch4int, STAT = ERR)

  THROW( ERR .NE. 0 )

! planck variables
  IF (Associated(coef%planck1)) DEALLOCATE (coef%planck1, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%planck2)) DEALLOCATE (coef%planck2, STAT = ERR)

  THROW( ERR .NE. 0 )

! frequency in GHz for MicroWaves

  IF (Associated(coef%frequency_ghz)) DEALLOCATE (coef%frequency_ghz, STAT = ERR)

  THROW( ERR .NE. 0 )


  IF (Associated(coef%dp)) DEALLOCATE (coef%dp, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%dpp)) DEALLOCATE (coef%dpp, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tstar)) DEALLOCATE (coef%tstar, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%to3star)) DEALLOCATE (coef%to3star, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%wstar)) DEALLOCATE (coef%wstar, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ostar)) DEALLOCATE (coef%ostar, STAT = ERR)

  THROW( ERR .NE. 0 )

! Specific variables for RTTOV8

  IF (Associated(coef%co2star)) DEALLOCATE (coef%co2star, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%co2star)) DEALLOCATE (coef%co2star, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%n2ostar)) DEALLOCATE (coef%n2ostar, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%costar)) DEALLOCATE (coef%costar, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ch4star)) DEALLOCATE (coef%ch4star, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tt_chn)) DEALLOCATE (coef%tt_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tt_val_chn)) DEALLOCATE (coef%tt_val_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tt_cwn)) DEALLOCATE (coef%tt_cwn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tt_a0)) DEALLOCATE (coef%tt_a0, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%tt_a1)) DEALLOCATE (coef%tt_a1, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ss_chn)) DEALLOCATE (coef%ss_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ss_val_chn)) DEALLOCATE (coef%ss_val_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ss_cwn)) DEALLOCATE (coef%ss_cwn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ss_solar_spectrum)) DEALLOCATE (coef%ss_solar_spectrum, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%woc_chn)) DEALLOCATE (coef%woc_chn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%woc_cwn)) DEALLOCATE (coef%woc_cwn, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%woc_waopc_ow)) DEALLOCATE (coef%woc_waopc_ow, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%woc_waopc_fw)) DEALLOCATE (coef%woc_waopc_fw, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ws_k_omega)) DEALLOCATE (coef%ws_k_omega, STAT = ERR)

  THROW( ERR .NE. 0 )

  IF (Associated(coef%ws_npoint)) DEALLOCATE (coef%ws_npoint, STAT = ERR)

  THROW( ERR .NE. 0 )

  Call rttov_nullify_coef (coef)

  CATCH
END SUBROUTINE rttov_dealloc_coef
