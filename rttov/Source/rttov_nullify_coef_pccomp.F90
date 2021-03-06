SUBROUTINE rttov_nullify_coef_pccomp(coef_pccomp)
! Description:
!   Nullify the PC coefficient structure
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
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_coef_pccomp
  IMPLICIT NONE
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
!INTF_END

  coef_pccomp%fmv_pc_comp_pc= 0_JPIM
  coef_pccomp%fmv_pc_sets   = 0_JPIM
  coef_pccomp%fmv_pc_mnum   = 0_JPIM
  coef_pccomp%fmv_pc_nchn   = 0_JPIM
  coef_pccomp%fmv_pc_nchn_noise = 0_JPIM
  coef_pccomp%fmv_pc_nche   = 0_JPIM
  coef_pccomp%fmv_pc_gas    = 0_JPIM
  coef_pccomp%fmv_pc_nlev   = 0_JPIM

  coef_pccomp%lim_pc_prfl_pmin     = 0_JPRB
  coef_pccomp%lim_pc_prfl_pmax     = 0_JPRB
  coef_pccomp%lim_pc_prfl_tsmin    = 0_JPRB
  coef_pccomp%lim_pc_prfl_tsmax    = 0_JPRB
  coef_pccomp%lim_pc_prfl_skmin    = 0_JPRB
  coef_pccomp%lim_pc_prfl_skmax    = 0_JPRB
  coef_pccomp%lim_pc_prfl_wsmin    = 0_JPRB
  coef_pccomp%lim_pc_prfl_wsmax    = 0_JPRB

  NULLIFY (coef_pccomp%eigenvectors)
  NULLIFY (coef_pccomp%emiss_chn)
  NULLIFY (coef_pccomp%emiss_c1)
  NULLIFY (coef_pccomp%emiss_c2)
  NULLIFY (coef_pccomp%emiss_c3)
  NULLIFY (coef_pccomp%emiss_c4)
  NULLIFY (coef_pccomp%emiss_c5)
  NULLIFY (coef_pccomp%emiss_c6)
  NULLIFY (coef_pccomp%emiss_c7)
  NULLIFY (coef_pccomp%emiss_c8)
  NULLIFY (coef_pccomp%emiss_c9)
  NULLIFY (coef_pccomp%ref_pc_prfl_p)
  NULLIFY (coef_pccomp%ref_pc_prfl_mr)
  NULLIFY (coef_pccomp%lim_pc_prfl_tmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_tmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_qmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_qmax)
  NULLIFY (coef_pccomp%lim_pc_prfl_ozmin)
  NULLIFY (coef_pccomp%lim_pc_prfl_ozmax)
  NULLIFY (coef_pccomp%co2_pc_ref)
  NULLIFY (coef_pccomp%n2o_pc_ref)
  NULLIFY (coef_pccomp%co_pc_ref)
  NULLIFY (coef_pccomp%ch4_pc_ref)
  NULLIFY (coef_pccomp%noise_in)
  NULLIFY (coef_pccomp%noise)
  NULLIFY (coef_pccomp%ff_ori_chn_in)
  NULLIFY (coef_pccomp%ff_cwn_in)
  NULLIFY (coef_pccomp%ff_bco_in)
  NULLIFY (coef_pccomp%ff_bcs_in)
  NULLIFY (coef_pccomp%planck1_in)
  NULLIFY (coef_pccomp%planck2_in)
  NULLIFY (coef_pccomp%pcreg)
END SUBROUTINE 
