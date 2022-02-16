!
SUBROUTINE rttov_read_binary_coef( &
            & ERR,           &
            & coef,          &
            & file_lu,       &
            & channels)
! Description:
!
! Read an binary coefficient file and fills coeff structure
!   arrays according to the optional list of channels.
!
! The user can provide an optional list of channels in "channels" argument
!  array to reduce the output coefficient structure to this list. This
! can be important for reducing the memory allocation required when running
! with advanced IR sounders (e.g. AIRS or IASI). If the user
!  wants all channels the "channels" argument shall not be present.
!
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
!  1.1       02/01/2003  A few comments added (R Saunders)
!  1.2       24/01/2003  add tests on all read statements (P Brunel)
!                        one record per channel for coefficients in binary format
!                        New header to allow checking R4<->R8
!  1.3       06/05/2003  Remove "optional" attribute of argument file_lu (P Brunel)
!  1.4       02/06/2004  New format for FMV section with RTTOV8  (P. Brunel)
!  1.5       15/06/2004  Corrected array dimension for coef % fmv_gas_pos (R Saunders)
!  1.6       14/05/2007  Updated for RTTOV-9 ( P Brunel)
!  1.7       12/12/2007  Option of additional layer at top flag read (R Saunders)
!  1.8       27/06/2008  Introduced the case when no channels are present for the
!                        phase function in the solar range (M. Matricardi)
!  1.9       27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.10      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.11      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.12      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
       & rttov_magic_string,     &
       & rttov_magic_number,     &
       & version_compatible_min, &
       & version_compatible_max, &
       & sensor_id_hi,           &
       & gas_id_mixed,           &
       & gas_id_watervapour,     &
       & gas_id_ozone,           &
       & gas_id_wvcont,          &
       & gas_id_co2,             &
       & gas_id_n2o,             &
       & gas_id_co,              &
       & gas_id_ch4,             &
       & ngases_max,             &
       & deg2rad
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu       ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels   (:)! list of channels to extract
! scalar arguments with intent(inout):
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
!INTF_END
#include "rttov_errorreport.h"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: n
  INTEGER(KIND=jpim) :: chn
  INTEGER(KIND=jpim) :: i, IncZeeman
! pointers for generic inputs
  COMPLEX(KIND=jprb), POINTER :: values_c_0   (:)
  COMPLEX(KIND=jprb), POINTER :: values_c_1   (:)
  REAL   (KIND=jprb), POINTER     :: values0         (:)
  REAL   (KIND=jprb), POINTER     :: values1         (:)
  REAL   (KIND=jprb), POINTER     :: values2         (:)
  REAL   (KIND=jprb), POINTER     :: values3         (:)
  REAL   (KIND=jprb), POINTER     :: values4         (:)
  INTEGER(KIND=jpim), POINTER     :: ivalues0        (:)
  INTEGER(KIND=jpim), POINTER     :: ivalues1        (:)
  CHARACTER(LEN = 16) :: bin_check_string
  REAL(KIND=jprb)     :: bin_check_number
  REAL(KIND=jprb)     :: bin_check_value
  CHARACTER(LEN = 80) :: errMessage
  LOGICAL(KIND=jplm)  :: section_present
!- End of header --------------------------------------------------------
  TRY
! 0 Initialise variables
!---------------------------------------------
! test presence of channels argument

  IF (Present(channels)) THEN
    all_channels = .FALSE.
  ELSE
    all_channels = .TRUE.
  ENDIF

! 3 Read binary file
!-------------------
! Binary file
  READ (file_lu, iostat=ERR)bin_check_string, bin_check_number

  THROWM( ERR .NE. 0, 'io status while reading header')

! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'Wrong header string in file')

! Verification of single/double precision using a 5 digit number
! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._JPRB - abs(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_JPRB .OR. bin_check_value < 0.99_JPRB) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'File created with a different real precision (R4<->R8)')

! COEF structure (V10)

    errMessage = 'io status while reading IDENTIFICATION'
    READ (file_lu, iostat=ERR)coef%id_platform, coef%id_sat, coef%id_inst, coef%id_sensor

    THROWM(err.ne.0,errMessage)

    READ (file_lu, iostat=ERR)coef%id_comp_lvl, coef%id_creation_date, coef%id_creation, coef%id_Common_name

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading FAST_MODEL_VARIABLES'

    READ (file_lu, iostat=ERR)coef%fmv_model_def, coef%fmv_model_ver, coef%fmv_chn, coef%fmv_gas


    THROWM(err.ne.0,errMessage)

    coef%id_comp_pc = 0_jpim

    READ (file_lu, iostat=ERR)coef%id_comp_pc

    THROWM(err.ne.0,errMessage)

    
    coef%IncZeeman = .FALSE.

    READ (file_lu, iostat=ERR)IncZeeman

    THROWM(err.ne.0,errMessage)

    IF (IncZeeman == 0) coef%IncZeeman = .FALSE.
    IF (IncZeeman == 1) coef%IncZeeman = .TRUE.

! Error if the compatibility version of the coefficient file
! is not in the range defined by the constant module

    IF (coef%id_comp_lvl < version_compatible_min .OR. coef%id_comp_lvl > version_compatible_max) THEN
      err = errorstatus_fatal

      THROWM(err.ne.0,"Version of coefficient file is incompatible with RTTOV library")

    ENDIF

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests
    file_channels = coef%fmv_chn

    IF (.NOT. all_channels) THEN
      coef%fmv_chn = Size(channels)
    ENDIF

    ALLOCATE (coef%fmv_gas_id(coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_gas_id")

    ALLOCATE (coef%fmv_gas_pos(ngases_max), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_gas_pos")

    ALLOCATE (coef%fmv_var(coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_var")

    ALLOCATE (coef%fmv_coe(coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_coe")
!pb
    ALLOCATE (coef%fmv_lvl(coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_lvl")

    ALLOCATE (coef%ff_ori_chn(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_ori_chn")

    ALLOCATE (coef%ff_val_chn(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_val_chn")

    ALLOCATE (coef%ff_cwn(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_cwn")

    ALLOCATE (coef%ff_bco(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_bco")

    ALLOCATE (coef%ff_bcs(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_bcs")

    ALLOCATE (coef%ff_gam(coef%fmv_chn), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ff_gam")

    ALLOCATE (coef%gaz_units(coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of gaz_units")

    coef%fmv_gas_id(:)  = 0_jpim
    coef%fmv_gas_pos(:) = 0_jpim
    coef%fmv_var(:)     = 0_jpim
    coef%fmv_lvl(:)     = 0_jpim
    coef%fmv_coe(:)     = 0_jpim
    READ (file_lu, iostat=ERR)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var, coef%fmv_lvl

    THROWM(err.ne.0,errMessage)


    READ (file_lu, iostat=ERR)coef%fmv_coe!pb

    THROWM(err.ne.0,errMessage)

    DO n = 1, coef%fmv_gas

      SELECT CASE (coef%fmv_gas_id(n))
      CASE (gas_id_mixed)
        coef%nmixed  = coef%fmv_var(n)
        coef%nlevels = coef%fmv_lvl(n)
        coef%ncmixed = coef%fmv_coe(n)!pb
      CASE (gas_id_watervapour)
        coef%nwater  = coef%fmv_var(n)
        coef%ncwater = coef%fmv_coe(n)!pb
      CASE (gas_id_ozone)
        coef%nozone  = coef%fmv_var(n)
        coef%ncozone = coef%fmv_coe(n)!pb
      CASE (gas_id_wvcont)
        coef%nwvcont  = coef%fmv_var(n)
        coef%ncwvcont = coef%fmv_coe(n)!pb
      CASE (gas_id_co2)
        coef%nco2  = coef%fmv_var(n)
        coef%ncco2 = coef%fmv_coe(n)!pb
      CASE (gas_id_n2o)
        coef%nn2o  = coef%fmv_var(n)
        coef%ncn2o = coef%fmv_coe(n)!pb
      CASE (gas_id_co)
        coef%nco  = coef%fmv_var(n)
        coef%ncco = coef%fmv_coe(n)!pb
      CASE (gas_id_ch4)
        coef%nch4  = coef%fmv_var(n)
        coef%ncch4 = coef%fmv_coe(n)!pb
      END SELECT

    ENDDO

    coef%nlayers = coef%nlevels - 1
    errMessage   = 'io status while reading GAZ_UNITS'
    READ (file_lu, iostat=ERR)coef%gaz_units

    THROWM(err.ne.0,errMessage)


! GAS_SPECTRAL_INTERVAL                  !pb
    READ (file_lu, iostat=ERR)coef%nintmixed, coef%nintwater, coef%nintozone, coef%nintwvcont, coef%nintco2,      &
      & coef%nintn2o, coef%nintco, coef%nintch4!pb

    THROWM(err.ne.0,errMessage)


    IF (coef%nintmixed > 0) THEN!pb
      ALLOCATE (coef%mixedgasint(2, coef%nintmixed), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of mixedgasint arra")

      READ (file_lu, iostat=ERR)coef%mixedgasint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintwater > 0) THEN!pb
      ALLOCATE (coef%watervapourint(2, coef%nintwater), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of watervapourint array")

      READ (file_lu, iostat=ERR)coef%watervapourint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintozone > 0) THEN!pb
      ALLOCATE (coef%ozoneint(2, coef%nintozone), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ozoneint array")

      READ (file_lu, iostat=ERR)coef%ozoneint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintwvcont > 0) THEN!pb
      ALLOCATE (coef%wvcontint(2, coef%nintwvcont), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of wvcontint array")

      READ (file_lu, iostat=ERR)coef%wvcontint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintco2 > 0) THEN!pb
      ALLOCATE (coef%co2int(2, coef%nintco2), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co2int array")

      READ (file_lu, iostat=ERR)coef%co2int!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintn2o > 0) THEN!pb
      ALLOCATE (coef%n2oint(2, coef%nintn2o), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of n2oint array")

      READ (file_lu, iostat=ERR)coef%n2oint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintco > 0) THEN!pb
      ALLOCATE (coef%coint(2, coef%nintco), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of coint array")

      READ (file_lu, iostat=ERR)coef%coint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintch4 > 0) THEN!pb
      ALLOCATE (coef%ch4int(2, coef%nintch4), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ch4int array")

      READ (file_lu, iostat=ERR)coef%ch4int!pb

      THROWM(err.ne.0,errMessage)

    ENDIF

      
    READ (file_lu)section_present
    errMessage = 'io status while reading TRANSMITTANCE_TRESHOLD'

    IF (section_present) THEN
!TRANSMITTANCE_TRESHOLD                  !pb
      ALLOCATE (coef%tt_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_chn")

      ALLOCATE (coef%tt_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_val_chn")

      ALLOCATE (coef%tt_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_cwn")

      ALLOCATE (coef%tt_a0(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_a0")

      ALLOCATE (coef%tt_a1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_a1")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%tt_chn, coef%tt_val_chn, coef%tt_cwn, coef%tt_a0, coef%tt_a1!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (ivalues1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues1")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values2(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values2 ")

        READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1, values2!pb

        THROWM(err.ne.0,errMessage)

        coef%tt_chn(:)     = ivalues0(channels(:))
        coef%tt_val_chn(:) = ivalues1(channels(:))
        coef%tt_cwn(:)     = values0(channels(:))
        coef%tt_a0(:)      = values1(channels(:))
        coef%tt_a1(:)      = values2(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values2")

      ENDIF

    ENDIF

    
    READ (file_lu)section_present
    errMessage = 'io status while reading SOLAR_SPECTRUM'

    IF (section_present) THEN
!SOLAR_SPECTRUM                          !pb
      ALLOCATE (coef%ss_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_chn")

      ALLOCATE (coef%ss_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_val_chn")

      ALLOCATE (coef%ss_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_cwn")

      ALLOCATE (coef%ss_solar_spectrum(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_solar_spectrum")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%ss_chn, coef%ss_val_chn, coef%ss_cwn, coef%ss_solar_spectrum!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (ivalues1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues1")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1

        THROWM(err.ne.0,errMessage)

        coef%ss_chn(:)            = ivalues0(channels(:))
        coef%ss_val_chn(:)        = ivalues1(channels(:))
        coef%ss_cwn(:)            = values0(channels(:))
        coef%ss_solar_spectrum(:) = values1(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

      ENDIF

    ENDIF

    
    READ (file_lu)section_present
    errMessage = 'io status while reading WATER_OPTICAL_CONSTANT'

    IF (section_present) THEN
!WATER_OPTICAL_CONSTANT                 !pb
      ALLOCATE (coef%woc_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_chn")

      ALLOCATE (coef%woc_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_cwn")

      ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_ow")

      ALLOCATE (coef%woc_waopc_fw(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_fw")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%woc_chn, coef%woc_cwn, coef%woc_waopc_ow, coef%woc_waopc_fw!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values_c_0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_0")

        ALLOCATE (values_c_1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_1")

        READ (file_lu, iostat=ERR)ivalues0, values1, values_c_0, values_c_1

        THROWM(err.ne.0,errMessage)

        coef%woc_chn(:)      = ivalues0(channels(:))
        coef%woc_cwn(:)      = values1(channels(:))
        coef%woc_waopc_ow(:) = values_c_0(channels(:))
        coef%woc_waopc_fw(:) = values_c_1(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values_c_0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values_c_0")

        DEALLOCATE (values_c_1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values_c_1")

      ENDIF

    ENDIF


    READ (file_lu)section_present
    errMessage = 'io status while reading WAVE_SPECTRUM'

    IF (section_present) THEN
!WAVE_SPECTRUM                          !pb
      READ (file_lu, iostat=ERR)coef%ws_nomega!pb

      THROWM(err.ne.0,errMessage)

      ALLOCATE (coef%ws_k_omega(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_k_omega")

      ALLOCATE (coef%ws_npoint(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_npoint")

      READ (file_lu, iostat=ERR)coef%ws_k_omega, coef%ws_npoint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF
    
    errMessage = 'io status while reading FILTER_FUNCTIONS'

    IF (all_channels) THEN
      READ (file_lu, iostat=ERR)coef%ff_ori_chn, coef%ff_val_chn, coef%ff_cwn, coef%ff_bco, coef%ff_bcs, coef%ff_gam

      THROWM(err.ne.0,errMessage)

    ELSE
      ALLOCATE (ivalues0(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ivalues0")

      ALLOCATE (ivalues1(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ivalues1")

      ALLOCATE (values0(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values0")

      ALLOCATE (values1(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values1")

      ALLOCATE (values2(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values2")

      ALLOCATE (values3(file_channels), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values3")

      READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1, values2, values3

      THROWM(err.ne.0,errMessage)

      coef%ff_ori_chn(:) = ivalues0(channels(:))
      coef%ff_val_chn(:) = ivalues1(channels(:))
      coef%ff_cwn(:)     = values0(channels(:))
      coef%ff_bco(:)     = values1(channels(:))
      coef%ff_bcs(:)     = values2(channels(:))
      coef%ff_gam(:)     = values3(channels(:))
      DEALLOCATE (ivalues0, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of ivalues0")

      DEALLOCATE (ivalues1, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of ivalues1")

      DEALLOCATE (values0, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values0")

      DEALLOCATE (values1, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values1")

      DEALLOCATE (values2, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values2")

      DEALLOCATE (values3, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values3")

    ENDIF

    errMessage = 'io status while reading FUNDAMENTAL_CONSTANTS'
    READ (file_lu, iostat=ERR)coef%fc_speedl, coef%fc_planck_c1, coef%fc_planck_c2, coef%fc_sat_height

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading EMISSIVITY model versions'
    READ (file_lu, iostat=ERR)coef%fastem_ver, coef%ssirem_ver

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading FASTEM'

    IF (coef%fastem_ver >= 1) THEN
      ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fastem_polar")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%fastem_polar

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        READ (file_lu, iostat=ERR)ivalues0

        THROWM(err.ne.0,errMessage)

        coef%fastem_polar(:) = ivalues0(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

      ENDIF

    ENDIF

    errMessage = 'io status while reading SSIREM'

    IF (coef%ssirem_ver >= 1) THEN
      ALLOCATE (coef%ssirem_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_chn")

      ALLOCATE (coef%ssirem_a0(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a0")

      ALLOCATE (coef%ssirem_a1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a1")

      ALLOCATE (coef%ssirem_a2(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a2")

      ALLOCATE (coef%ssirem_xzn1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_xzn1")

      ALLOCATE (coef%ssirem_xzn2(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_xzn2")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)     &
          & coef%ssirem_chn, coef%ssirem_a0, coef%ssirem_a1, coef%ssirem_a2, coef%ssirem_xzn1, coef%ssirem_xzn2

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values2(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values2")

        ALLOCATE (values3(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values3")

        ALLOCATE (values4(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values4")

        READ (file_lu, iostat=ERR)ivalues0, values0, values1, values2, values3, values4

        THROWM(err.ne.0,errMessage)

        coef%ssirem_chn(:)  = ivalues0(channels(:))
        coef%ssirem_a0(:)   = values0(channels(:))
        coef%ssirem_a1(:)   = values1(channels(:))
        coef%ssirem_a2(:)   = values2(channels(:))
        coef%ssirem_xzn1(:) = values3(channels(:))
        coef%ssirem_xzn2(:) = values4(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values2")

        DEALLOCATE (values3, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values3")

        DEALLOCATE (values4, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values4")

      ENDIF

    ENDIF

    ALLOCATE (coef%ref_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_p")

    ALLOCATE (coef%ref_prfl_t(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_t")

    ALLOCATE (coef%ref_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_mr")

    ALLOCATE (coef%lim_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_p")

    ALLOCATE (coef%lim_prfl_tmax(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_tmax")

    ALLOCATE (coef%lim_prfl_tmin(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_tmin")

    ALLOCATE (coef%lim_prfl_gmin(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_gmin")

    ALLOCATE (coef%lim_prfl_gmax(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_gmax")

    errMessage = 'io status while reading REFERENCE PROFILE'
    READ (file_lu, iostat=ERR)coef%ref_prfl_p, coef%ref_prfl_t, coef%ref_prfl_mr

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading PROFILE LIMITS'
    READ (file_lu, iostat=ERR)     &
      & coef%lim_prfl_p, coef%lim_prfl_tmax, coef%lim_prfl_tmin, coef%lim_prfl_gmax, coef%lim_prfl_gmin

    THROWM(err.ne.0,errMessage)

! FAST COEFFICIENT section
    errMessage = 'io status while reading Mixed gases coefs'

    IF (coef%ncmixed > 0) THEN
      ALLOCATE (coef%mixedgas(coef%nlayers, coef%fmv_chn, coef%ncmixed), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of mixedgas")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%mixedgas(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%mixedgas(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading Water vapour coefs'

    IF (coef%ncwater > 0) THEN
      ALLOCATE (coef%watervapour(coef%nlayers, coef%fmv_chn, coef%ncwater), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of watervapour")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%watervapour(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%watervapour(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading Ozone coefs'

    IF (coef%ncozone > 0) THEN
      ALLOCATE (coef%ozone(coef%nlayers, coef%fmv_chn, coef%ncozone), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ozone")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%ozone(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%ozone(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading WV continuum coefs'

    IF (coef%ncwvcont > 0) THEN
      ALLOCATE (coef%wvcont(coef%nlayers, coef%fmv_chn, coef%ncwvcont), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of wvcont")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%wvcont(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%wvcont(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading CO2 coefs'

    IF (coef%ncco2 > 0) THEN
      ALLOCATE (coef%co2(coef%nlayers, coef%fmv_chn, coef%ncco2), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co2 ")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%co2(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%co2(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading N2O coefs'

    IF (coef%ncn2o > 0) THEN
      ALLOCATE (coef%n2o(coef%nlayers, coef%fmv_chn, coef%ncn2o), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of n2o")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%n2o(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%n2o(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading CO coefs'

    IF (coef%ncco > 0) THEN
      ALLOCATE (coef%co(coef%nlayers, coef%fmv_chn, coef%ncco), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%co(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%co(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

    ENDIF

    errMessage = 'io status while reading CH4 coefs'

    IF (coef%ncch4 > 0) THEN
      ALLOCATE (coef%ch4(coef%nlayers, coef%fmv_chn, coef%ncch4), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ch4")

      i = 1

      DO chn = 1, file_channels

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%ch4(:, chn, :)
        ELSE IF (chn == channels(i)) THEN
          READ (file_lu, iostat=ERR)coef%ch4(:, i, :)

          IF (i < coef%fmv_chn) THEN
            i = i + 1
          ENDIF

        ELSE
          READ (file_lu, iostat=ERR)
        ENDIF


        THROWM(err.ne.0,errMessage)

      ENDDO

  ENDIF

!
! Here add reading of new sections for binary format in order to keep compatibility with
! previous versions
!
  CATCH
END SUBROUTINE rttov_read_binary_coef
