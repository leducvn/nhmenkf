!
SUBROUTINE rttov_init_coef(  ERR, coef )
! Description:
!
!   coef arrays initialisation for all pes
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
!  1.0       17/05/2004  Original (based on rttov_readcoeffs)
!  1.1       01/06/2005  Marco Matricardi (ECMWF):
!               --       Added variables for CO2,N2O, CO and CH4.
!  1.2       2/11/2007   Made mwcldtop variable (R Saunders)
!  1.3       2/11/2007   Included polarimetric type for GHz calc (RS)
!  1.4       27/02/2009  Profile levels to include ToA. Distinguish between
!                        layer arrays and level arrays - size, index labels,
!                        looping. Prepares predictor variables to agree
!                        with RTTOV-9 (P. Rayer)
!  1.5       16/7/2009   Remove trace gases for RTTOV-8 and 9 which are not supported (RS)
!  1.6       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw,       &
       & sensor_id_po,       &
       & gas_id_mixed,       &
       & gas_id_watervapour, &
       & gas_id_ozone,       &
       & gas_id_co2,         &
       & gas_id_n2o,         &
       & gas_id_co,          &
       & gas_id_ch4,         &
       & gas_unit_specconc,  &
       & gas_unit_ppmv,      &
       & earthradius,        &
       & mwcldtp
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR          ! return code
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef         ! coefficients
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_q2v.h"
! Local Scalars:
  INTEGER(KIND=jpim) :: i, n, h
!- End of header --------------------------------------------------------
  TRY
! 0 Initialise variables
!---------------------------------------------

! 5 Now compute auxillary variables and change unit for mixing ratios
!--------------------------------------------------------------------
! ratio satellite altitude to earth radius
  coef%ratoe = (earthradius + coef%fc_sat_height) / earthradius
! planck variables
  ALLOCATE (coef%planck1(coef%fmv_chn), STAT = ERR)

  THROWM( ERR .NE. 0, "allocation of Planck1" )

  ALLOCATE (coef%planck2(coef%fmv_chn), STAT = ERR)

  THROWM( ERR .NE. 0, "allocation of Planck2" )


  coef%planck1(:) = coef%fc_planck_c1 * coef%ff_cwn(:) ** 3
  coef%planck2(:) = coef%fc_planck_c2 * coef%ff_cwn(:)

! frequency in GHz for MicroWaves

  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
    ALLOCATE (coef%frequency_ghz(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of frequency_ghz" )

    coef%frequency_ghz(:) = coef%fc_speedl * 1.0e-09_JPRB * coef%ff_cwn(:)
  ENDIF

! Conversion of gas mixing ratio units, if needed
! The correct unit for RTTOV calculations is ppmv
! Take care this conversion is only valid
! if all gases have the same unit
  h = coef%fmv_gas_pos(gas_id_watervapour)

  DO n = 1, coef%fmv_gas

    IF (coef%gaz_units(n) == gas_unit_specconc) THEN
! Unit of gaz n is specific concentration

      DO i = 1, coef%nlevels
! Convert reference profile mixing ratio
        CALL rttov_q2v( &
              & coef%gaz_units(h),      &
              & coef%ref_prfl_mr(i, h), &
              & coef%fmv_gas_id(n),     &
              & coef%ref_prfl_mr(i, n), &
              & coef%ref_prfl_mr(i, n))
! Now unit of gaz n is ppmv for the reference
! in particular for H2O
        coef%gaz_units(n) = gas_unit_ppmv
! Convert profile minimum limit mixing ratio
        CALL rttov_q2v( &
              & coef%gaz_units(h),        &
              & coef%ref_prfl_mr(i, h),   &
              & coef%fmv_gas_id(n),       &
              & coef%lim_prfl_gmin(i, n), &
              & coef%lim_prfl_gmin(i, n))
! Convert profile maximum limit mixing ratio
        CALL rttov_q2v( &
              & coef%gaz_units(h),        &
              & coef%ref_prfl_mr(i, h),   &
              & coef%fmv_gas_id(n),       &
              & coef%lim_prfl_gmax(i, n), &
              & coef%lim_prfl_gmax(i, n))
      ENDDO

    ENDIF

  ENDDO

! 6 Construct variables from reference profile
! --------------------------------------------
! 6.1   Specific variables for RTTOV7
! -----------------------------------
! Test on model and coeff versions

  IF (coef%fmv_model_ver == 7) THEN
    ALLOCATE (coef%dp(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dp" )

    ALLOCATE (coef%dpp(0:coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dpp" )

    ALLOCATE (coef%tstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of tstar" )

    ALLOCATE (coef%to3star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of to3star" )

    ALLOCATE (coef%wstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of wstar" )

    ALLOCATE (coef%ostar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ostar" )

! pressure intervals dp
! ------------------
! dp(layer=N), where layer N lies between levels N and N+1
    coef%dp(1:coef%nlayers)      = coef%ref_prfl_p(2:coef%nlevels) - coef%ref_prfl_p(1:coef%nlevels - 1)
! pressure quantity dpp
! -----------------
! special coef % dpp(0) replaces RTTOV-9 coef % dpp(1)
! needed for predictor variables ww and ow (but not tw)
    coef%dpp(0)                  = coef%dp(1) * coef%ref_prfl_p(1)
! NB coef % ref_prfl_p(1) replaces RTTOV-9 parameter pressure_top
! coef % dpp(layer=N), where layer N bounded by levels N and N+1
! coef % dpp(N) replaces RTTOV-9 coef % dpp(N+1) -  NB coef % dpp(43) not used
    coef%dpp(1:coef%nlayers - 1) = coef%dp(1:coef%nlayers - 1) * coef%ref_prfl_p(2:coef%nlevels - 1)

    DO i = 1, coef%nlevels
      IF (coef%ref_prfl_p(i) < mwcldtp) coef%mwcldtop = i
    ENDDO

! reference layer quantities
! --------------------------
! temperature
    n = coef%fmv_gas_pos(gas_id_mixed)
    coef%tstar(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
! temperature for O3 profiles

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%to3star(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
    ENDIF

! water vapour
    n = coef%fmv_gas_pos(gas_id_watervapour)
    coef%wstar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
! ozone

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%ostar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
    ENDIF

  ENDIF

! 6.2   Specific variables for RTTOV8
! AT PRESENT SAME VARIABLES AS RTTOV7
! -----------------------------------
! Test on model and coeff versions

  IF (coef%fmv_model_ver == 8) THEN
    ALLOCATE (coef%dp(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dp" )

    ALLOCATE (coef%dpp(0:coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dpp" )

    ALLOCATE (coef%tstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of tstar" )

    ALLOCATE (coef%to3star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of to3star" )

    ALLOCATE (coef%wstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of wstar" )

    ALLOCATE (coef%ostar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ostar" )

    ALLOCATE (coef%co2star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of co2star" )

! pressure intervals dp
! ------------------
! dp(layer=N), where layer N lies between levels N and N+1
    coef%dp(1:coef%nlayers)      = coef%ref_prfl_p(2:coef%nlevels) - coef%ref_prfl_p(1:coef%nlevels - 1)
! pressure quantity dpp
! -----------------
! special coef % dpp(0) replaces RTTOV-9 coef % dpp(1)
! needed for predictor variables ww and ow (but not tw)
    coef%dpp(0)                  = coef%dp(1) * coef%ref_prfl_p(1)
! NB coef % ref_prfl_p(1) replaces RTTOV-9 parameter pressure_top
! coef % dpp(layer=N), where layer N bounded by levels N and N+1
! coef % dpp(N) replaces RTTOV-9 coef % dpp(N+1) -  NB coef % dpp(43) not used
    coef%dpp(1:coef%nlayers - 1) = coef%dp(1:coef%nlayers - 1) * coef%ref_prfl_p(2:coef%nlevels - 1)

    DO i = 1, coef%nlevels
      IF (coef%ref_prfl_p(i) < mwcldtp) coef%mwcldtop = i
    ENDDO

! reference layer quantities
! --------------------------
! temperature
    n = coef%fmv_gas_pos(gas_id_mixed)
    coef%tstar(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
! temperature for O3 profiles

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%to3star(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
    ENDIF

! water vapour
    n = coef%fmv_gas_pos(gas_id_watervapour)
    coef%wstar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
! ozone

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%ostar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
    ENDIF

! CO2

    IF (coef%nco2 > 0) THEN
      n = coef%fmv_gas_pos(gas_id_co2)
      coef%co2star(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
    ENDIF

  ENDIF

! 6.3   Specific variables for RTTOV9
! AT PRESENT SAME VARIABLES AS RTTOV8
! -----------------------------------
! Test on model and coeff versions

  IF (coef%fmv_model_ver == 9) THEN
    ALLOCATE (coef%dp(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dp" )

    ALLOCATE (coef%dpp(0:coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of dpp" )

    ALLOCATE (coef%tstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of tstar" )

    ALLOCATE (coef%to3star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of to3star" )

    ALLOCATE (coef%wstar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of wstar" )

    ALLOCATE (coef%ostar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ostar" )

    ALLOCATE (coef%co2star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of co2star" )

    ALLOCATE (coef%n2ostar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of n2ostar" )

    ALLOCATE (coef%costar(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of costar" )

    ALLOCATE (coef%ch4star(coef%nlayers), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ch4star" )


! pressure intervals dp
! ------------------
! dp(layer=N), where layer N lies between levels N and N+1
    coef%dp(1:coef%nlayers)      = coef%ref_prfl_p(2:coef%nlevels) - coef%ref_prfl_p(1:coef%nlevels - 1)
! pressure quantity dpp
! -----------------
! special coef % dpp(0) replaces RTTOV-9 coef % dpp(1)
! needed for predictor variables ww and ow (but not tw)
    coef%dpp(0)                  = coef%dp(1) * coef%ref_prfl_p(1)
! NB coef % ref_prfl_p(1) replaces RTTOV-9 parameter pressure_top
! coef % dpp(layer=N), where layer N bounded by levels N and N+1
! coef % dpp(N) replaces RTTOV-9 coef % dpp(N+1) -  NB coef % dpp(43) not used
    coef%dpp(1:coef%nlayers - 1) = coef%dp(1:coef%nlayers - 1) * coef%ref_prfl_p(2:coef%nlevels - 1)

    DO i = 1, coef%nlevels
      IF (coef%ref_prfl_p(i) < mwcldtp) coef%mwcldtop = i
    ENDDO

! reference layer quantities
! --------------------------
! temperature
    n = coef%fmv_gas_pos(gas_id_mixed)
    coef%tstar(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
! temperature for O3 profiles

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%to3star(1:coef%nlayers) = (coef%ref_prfl_t(1:coef%nlevels - 1, n) + coef%ref_prfl_t(2:coef%nlevels, n)) / 2
    ENDIF

! water vapour
    n = coef%fmv_gas_pos(gas_id_watervapour)
    coef%wstar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
! ozone

    IF (coef%nozone > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ozone)
      coef%ostar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
    ENDIF

! CO2

    IF (coef%nco2 > 0) THEN
      n = coef%fmv_gas_pos(gas_id_co2)
      coef%co2star(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2

    ENDIF

! N2O

    IF (coef%nn2o > 0) THEN
      n = coef%fmv_gas_pos(gas_id_n2o)
      coef%n2ostar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2

    ENDIF

! CO

    IF (coef%nco > 0) THEN
      n = coef%fmv_gas_pos(gas_id_co)
      coef%costar(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2

    ENDIF

! CH4

    IF (coef%nch4 > 0) THEN
      n = coef%fmv_gas_pos(gas_id_ch4)
      coef%ch4star(1:coef%nlayers) = (coef%ref_prfl_mr(1:coef%nlevels - 1, n) + coef%ref_prfl_mr(2:coef%nlevels, n)) / 2
    ENDIF

  ENDIF

  CATCH
END SUBROUTINE 
