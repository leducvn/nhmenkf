SUBROUTINE rttov_write_ascii_coef(err, coef, file_id)
! Description:
! write on unit file_id the coef structure.
! If lbinary is false or not present the file is assumed as
! an ASCII sequential formatted, in other case it is sequential unformatted.
! I/O write status are only tested at the end of the code
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
!  1.1       24/01/2003  insert I/O status (P Brunel)
!                        one record per channel for coefficients in binary format
!                        New header to allow checking R4<->R8
!  1.2       02/06/2004  Update for RTTOV8 coefs (P. Brunel)
!  1.3       02/08/2006  Change format for pressure levels f8.3 -> f9.4 (P. Brunel)
!  1.4       14/05/2007  Updated for RTTOV-9 (P Brunel)
!  1.5       19/02/2008  Another update for RTTOV-9 for inc_top (R Saunders)
!  1.6       27/06/2008  Introduced the case where no channels are available for
!                        the phase function in the solar range (M. Matricardi)
!  1.7       06/03/2009  Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.8       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :    &
       & version,            &
       & release,            &
       & minor_version,      &
       & rttov_magic_string, &
       & rttov_magic_number, &
       & gas_id_mixed,       &
       & gas_id_watervapour, &
       & gas_id_ozone,       &
       & gas_id_wvcont,      &
       & gas_id_co2,         &
       & gas_id_n2o,         &
       & gas_id_co,          &
       & gas_id_ch4,         &
       & gas_name,           &
       & gas_unit_name
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim), INTENT(IN)  :: file_id! file logical unit number
  TYPE(rttov_coef)  , INTENT(IN)  :: coef   ! coefficients
! scalar arguments with intent(in):
  INTEGER(KIND=jpim), INTENT(OUT) :: err    ! return code
!INTF_END
#include "rttov_errorreport.h"
! local scalars
  INTEGER(KIND=jpim)  :: i, j, l, k
  CHARACTER(LEN = 2 ) :: sensor
  CHARACTER(LEN = 32) :: section
  CHARACTER(LEN = 80) :: errMessage
  CHARACTER(LEN = 80) :: version_name
  ! Integer(Kind=jpim) :: Err_Unit         Logical error unit 
  ! Integer(Kind=jpim) :: verbosity_level
  INTEGER(KIND=jpim)  :: IncZeeman
!- End of header --------------------------------------------------------
  TRY
  !  Err_unit=-1_jpim 
  !  verbosity_level=3_jpim 
  !  Call rttov_errorhandling (Err_unit, verbosity_level) 
  
! Consider lbinary option to create the option
!ASCII file
  WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
  INFO(errMessage)

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV coefficient file '//Trim(coef%id_Common_name)

  THROW(err.ne.0)
 
  WRITE (file_id, '(a)', iostat=err)' ! automatic creation by subroutine rttov_writecoef'
  
  THROW(err.ne.0)
  
  IF (release < 10 .and. minor_version < 10 ) THEN
    WRITE (version_name, '(I2.2,".",i1,".",i1)', iostat=err) version, release, minor_version
  ELSE
    WRITE (version_name, '(I2.2,".",i2.2,".",i2.2)', iostat=err) version, release, minor_version
  ENDIF

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV library version '//TRIM(version_name)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

! COEF structure (V10)
! IDENTIFICATION
  section = 'IDENTIFICATION'

  SELECT CASE (coef%id_sensor)
  CASE (1_jpim)
    sensor = 'ir'
  CASE (2_jpim)
    sensor = 'mw'
  CASE (3_jpim)
    sensor = 'hi'
  CASE (4_jpim)
    sensor = 'po'
  END SELECT


  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(3i3,T20,a)', iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, '! platform sat_id instrument'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_Common_name)

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a,T20,a)', iostat=err)sensor, '! sensor type [ir,mw,hi]'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i2,T20,a)', iostat=err)coef%id_comp_lvl, '! RTTOV coefficient file version number'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_creation)

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,1x,i2.2,1x,i2.2,t20,a)', iostat=err)coef%id_creation_date, '! creation date'

  THROW(err.ne.0)

! No LINE-BY-LINE section

  IF (coef%line_by_line(1) .NE. 'xxxx') THEN
    section = 'LINE-BY-LINE'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)


    DO i = 1, Size(coef%line_by_line)
      IF (coef%line_by_line(i) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)Trim(coef%line_by_line(i))

      THROW(err.ne.0)

    ENDDO

  ENDIF

! FAST_MODEL_VARIABLES
  section = 'FAST_MODEL_VARIABLES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a,t20,a)', iostat=err)coef%fmv_model_def, '! fast model name'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_model_ver, '! fast model version compatibility level'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_chn, '! Number of channels described in the coef file'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_gas, '! Number of gases described in the coef file'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)' , iostat=ERR)coef%id_comp_pc, '! PC compatibility level'

  THROW(err.ne.0)

  IncZeeman = 0_jpim
  IF (coef%IncZeeman) IncZeeman = 1_jpim
  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)IncZeeman, '! Zeeman flag'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(coef%fmv_gas_id(i))), '! gas identification'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,3i4,t20,a)', iostat=err)coef%fmv_var(i), coef%fmv_coe(i), coef%fmv_lvl(i),&
            & '! variables/predictors  levels (pressure/absorber)'

    THROW(err.ne.0)

  ENDDO


! GAZ_UNITS
  section = 'GAZ_UNITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gaz concentrations can be expressed in '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! volume mixing ratio (ppmv)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! specific concentration (kg/kg)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i4,t20,"! ",a)', iostat=err)coef%gaz_units(i), gas_unit_name(coef%gaz_units(i))

    THROW(err.ne.0)

  ENDDO


  IF ((coef%nintmixed > 0) .AND. (coef%nintwater > 0) .AND. (coef%nintozone > 0) .AND. (coef%nintwvcont > 0) .AND. &
      (coef%nintco2 > 0) .AND. (coef%nintn2o > 0) .AND. (coef%nintco > 0) .AND. (coef%nintch4 > 0)) THEN!pb
    section = 'GAS_SPECTRAL_INTERVAL'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)


    IF (coef%nintmixed > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_mixed)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintmixed

      THROW(err.ne.0)


      DO i = 1, coef%nintmixed
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%mixedgasint(1, i), coef%mixedgasint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintwater > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_watervapour)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintwater

      THROW(err.ne.0)


      DO i = 1, coef%nintwater
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%watervapourint(1, i), coef%watervapourint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintozone > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_ozone)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintozone

      THROW(err.ne.0)


      DO i = 1, coef%nintozone
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%ozoneint(1, i), coef%ozoneint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintwvcont > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_wvcont)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintwvcont

      THROW(err.ne.0)


      DO i = 1, coef%nintwvcont
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%wvcontint(1, i), coef%wvcontint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintco2 > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_co2)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintco2

      THROW(err.ne.0)


      DO i = 1, coef%nintco2
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%co2int(1, i), coef%co2int(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintn2o > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_n2o)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintn2o

      THROW(err.ne.0)


      DO i = 1, coef%nintn2o
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%n2oint(1, i), coef%n2oint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintco > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_co)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintco

      THROW(err.ne.0)


      DO i = 1, coef%nintco
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%coint(1, i), coef%coint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintch4 > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_ch4)), '! gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintch4

      THROW(err.ne.0)


      DO i = 1, coef%nintch4
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%ch4int(1, i), coef%ch4int(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDIF


  IF (Associated(coef%tt_chn)) THEN!pb
    section = 'TRANSMITTANCE_TRESHOLD'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! chan number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! validity of channel '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! central wave number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! transmittance treshold'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! transmittance value   '

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5,3e19.10)', iostat=err)     &
        & coef%tt_chn(i), coef%tt_val_chn(i), coef%tt_cwn(i), coef%tt_a0(i), coef%tt_a1(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb

  IF (Associated(coef%ss_chn)) THEN!pb
    section = 'SOLAR_SPECTRUM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! chan number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! validity of channel'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! central wave number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! solar spectrum'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5,2e19.10)', iostat=err)     &
        & coef%ss_chn(i), coef%ss_val_chn(i), coef%ss_cwn(i), coef%ss_solar_spectrum(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb

  IF (associated(coef%woc_chn)) THEN!pb
    section = 'WATER_OPTICAL_CONSTANT'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! chan number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! central wave number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! ocean water optical constants(real and imaginary part)'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! fresh water optical constants(real and imaginary part)'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(i5,e19.10,2(" (",e17.10,",",e17.10,")"))', iostat=err)     &
        & coef%woc_chn(i), coef%woc_cwn(i), coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb

  IF (Associated(coef%ws_npoint)) THEN!pb
    section = 'WAVE_SPECTRUM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Number of points'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Point number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! WAve spectrum'

    THROW(err.ne.0)

    WRITE (file_id,  * , iostat=err)coef%ws_nomega

    THROW(err.ne.0)


    DO i = 1, coef%ws_nomega
      WRITE (file_id, '(f10.3,f12.5)', iostat=err)coef%ws_npoint(i), coef%ws_k_omega(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb
  section = 'FILTER_FUNCTIONS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel Number (from instrument original description)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel status '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Central Wavenumber'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Band Correction coefficients(Offset,Slope)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gamma correction factor'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_chn
    WRITE (file_id, '(1x,i4,1x,i4,4(1x,e18.10))', iostat=err)     &
      & coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)

    THROW(err.ne.0)

  ENDDO

  section = 'FUNDAMENTAL_CONSTANTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! units of constants for spectral radiance'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! first radiation constant(mW/(m2.sr.cm-4))'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! second radiation constant (cm.K)'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,f14.1,t30,a)', iostat=err)coef%fc_speedl, '! speed of light (cm/s)'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,1p,e15.8,0p,f10.6,t30,a)', iostat=err)coef%fc_planck_c1, coef%fc_planck_c2, '! Planck constants'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,f8.1,t30,a)', iostat=err)coef%fc_sat_height, '! nominal satellite height (km)'

  THROW(err.ne.0)


  IF (coef%fastem_ver >= 1) THEN
    section = 'FASTEM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! S. English fast generic millimetre wave ocean emissivity model'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Polarisation of each channel', &             
              & ' !       MPOL=0: Average of vertical and horizontal polarisation ie 0.5(H+V)', &
              & ' !       MPOL=1: Nominal vertical at nadir rotating with view angle QV', &
              & ' !       MPOL=2: Nominal horizontal at nadir rotating with view angle QH', &
              & ' !       MPOL=3: Vertical V', &
              & ' !       MPOL=4: Horizontal H', &
              & ' !       MPOL=5: +45 minus -45 (3rd stokes vector) S3', &
              & ' !       MPOL=6: Left circular - right circular (4th stokes vector) S4'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)coef%fastem_ver, '   ! version number'

    THROW(err.ne.0)

    WRITE (file_id, '(20i3)', iostat=err)(coef%fastem_polar(i), i = 1, coef%fmv_chn)

    THROW(err.ne.0)

  ENDIF


  IF (coef%ssirem_ver >= 1) THEN
    section = 'SSIREM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel ID'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! 5 coefficients for emissivity model ssirem'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)coef%ssirem_ver, '   ! version number'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(1x,i4,3f12.7,2f4.1)', iostat=err)coef%ssirem_chn(i), coef%ssirem_a0(i), coef%ssirem_a1(i),      &
        & coef%ssirem_a2(i), coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF

  section = 'REFERENCE_PROFILE'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Ref.pressure (hPa)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Ref.Temp (K) Ref.Volume Mixing Ratio [ppmv] for each gas'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Note for MxG that mixing ratio is "missing"'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)


    DO l = 1, coef%fmv_lvl(i)
      WRITE (file_id, '(1x,f9.4,2x,f7.3,1x,e13.6)')coef%ref_prfl_p(l), coef%ref_prfl_t(l, i), coef%ref_prfl_mr(l, i)
!!$                & coef % ref_prfl_p(l), coef % ref_prfl_t(l,i), ref_mr(l,i)

      THROW(err.ne.0)

    ENDDO

  ENDDO

  section = 'PROFILE_LIMITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Ref.pressure (hPa)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Temp Max (K) Temp Min (K)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Volume Mixing Ratio for  Max and Min [ppmv] for each gas'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !      Temperature'

  THROW(err.ne.0)


  DO l = 1, coef%fmv_lvl(1)
    WRITE (file_id, '(1x,f9.4,2(1x,f7.2))', iostat=err)coef%lim_prfl_p(l), coef%lim_prfl_tmax(l), coef%lim_prfl_tmin(l)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)


    DO l = 1, coef%fmv_lvl(i)
      WRITE (file_id, '(1x,f9.4,2x,e12.4,e12.4)', iostat=err)     &
        & coef%lim_prfl_p(l), coef%lim_prfl_gmax(l, i), coef%lim_prfl_gmin(l, i)

      THROW(err.ne.0)

    ENDDO

  ENDDO

  section = 'FAST_COEFFICIENTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! transmission coefficients'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Order of the gases:'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)

  ENDDO


  DO l = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)gas_name(coef%fmv_gas_id(l))

    THROW(err.ne.0)


    SELECT CASE (coef%fmv_gas_id(l))
    CASE (gas_id_mixed)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%mixedgas(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncmixed)

      THROW(err.ne.0)

    CASE (gas_id_watervapour)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%watervapour(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwater)

      THROW(err.ne.0)

    CASE (gas_id_ozone)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%ozone(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncozone)

      THROW(err.ne.0)

    CASE (gas_id_wvcont)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%wvcont(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwvcont)

      THROW(err.ne.0)

    CASE (gas_id_co2)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%co2(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco2)

      THROW(err.ne.0)

    CASE (gas_id_n2o)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%n2o(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncn2o)

      THROW(err.ne.0)

    CASE (gas_id_co)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%co(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco)

      THROW(err.ne.0)

    CASE (gas_id_ch4)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%ch4(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncch4)

      THROW(err.ne.0)

    END SELECT

  ENDDO

  section = 'END'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
