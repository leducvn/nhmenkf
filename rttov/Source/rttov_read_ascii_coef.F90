!
SUBROUTINE rttov_read_ascii_coef( &
            & err,           &
            & coef,          &
            & file_lu,       &
            & channels)
! Description:
!
! Read an ASCII coefficient file and fills coeff structure
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
!  1.2       24/01/2003  Add return when section END encountered (P Brunel)
!                        any I/O error is coded as fatal
!                        Add GAZ_UNITS section
!  1.3       02/06/2004  New format for FMV section with RTTOV8 (P. Brunel)
!  1.4       15/06/2004  Corrected array dimension for coef % fmv_gas_pos (R Saunders)
!  1.5        June 2005  Added Additional arrays for RTTOV8M Marco Matricardi
!  1.6       05/12/2005  Corrected some array types for optical  const (R Saunders)
!  1.7       19/03/2007  Reduced no of continuation lines to below 39 (R Saunders)
!  1.8       12/12/2007  Add option of using top level (R Saunders)
!  1.9       01/11/2007  Removed hardcoded section length (A Geer)
!  1.10      26/06/2008  Introduced the case where no channels are available for the
!                        phase function in the solar range (M. Matricardi)
!  1.11      27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.12      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.13      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.14      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
       & version_compatible_min, &
       & version_compatible_max, &
       & sensor_id_hi,           &
       & sensor_id_mw,           &
       & sensor_id_ir,           &
       & sensor_id_po,           &
       & gas_id_mixed,           &
       & gas_id_watervapour,     &
       & gas_id_ozone,           &
       & gas_id_wvcont,          &
       & gas_id_co2,             &
       & gas_id_n2o,             &
       & gas_id_co,              &
       & gas_id_ch4,             &
       & sensor_name,            &
       & ngases_max,             &
       & gas_name,               &
       & gas_unit_specconc,      &
       & deg2rad,                &
       & lensection
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
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: err           ! return code
!INTF_END
#include "rttov_opencoeff.h"
#include "rttov_errorreport.h"
#include "rttov_skipcommentline.h"
#include "rttov_deletecomment.h"
#include "rttov_cmpuc.h"
#include "rttov_findnextsection.h"


#include "rttov_nullify_coef.h"
#include "rttov_nullify_coef_pccomp.h"
#include "rttov_nullify_coef_scatt_ir.h"
#include "rttov_nullify_optpar_ir.h"

! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  LOGICAL(KIND=jplm) :: for_output
  INTEGER(KIND=jpim) :: file_lu_coef
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: io_status
  REAL   (KIND=jprb) :: pres
  INTEGER(KIND=jpim) :: i, j, k, l, n, IncZeeman
  INTEGER(KIND=jpim) :: gas_id
  INTEGER(KIND=jpim) :: nbint
! pointers for generic inputs
  INTEGER(KIND=jpim) :: nvalues
  COMPLEX(KIND=jprb), POINTER :: values_c_0   (:)
  COMPLEX(KIND=jprb), POINTER :: values_c_1   (:)
  REAL   (KIND=jprb), POINTER :: values0      (:)
  REAL   (KIND=jprb), POINTER :: values1      (:)
  REAL   (KIND=jprb), POINTER :: values2      (:)
  REAL   (KIND=jprb), POINTER :: values3      (:)
  REAL   (KIND=jprb), POINTER :: values4      (:)
  INTEGER(KIND=jpim), POINTER :: ivalues0     (:)
  INTEGER(KIND=jpim), POINTER :: ivalues1     (:)
  REAL   (KIND=jprb), POINTER :: coeffsarray  (:, :, :)

  CHARACTER(LEN = 36)         :: input_string
  CHARACTER(LEN = 32)         :: gas_Type
  CHARACTER(LEN = lensection) :: section
  CHARACTER(LEN = 100       ) :: errMessage
  LOGICAL(KIND=jplm)          :: found
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

  CALL rttov_nullify_coef(coef)

!read the file

  readfile : DO
    CALL rttov_findnextsection(file_lu, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (Trim(section))

    CASE ('IDENTIFICATION')
! Identification section
! 6 lines
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef%id_platform, coef%id_sat, coef%id_inst

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu, '(a)', iostat=ERR)coef%id_Common_name

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)input_string

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


      SELECT CASE (input_string)
      CASE (sensor_name(sensor_id_ir))
        coef%id_sensor = sensor_id_ir! Infrared
      CASE (sensor_name(sensor_id_mw))
        coef%id_sensor = sensor_id_mw! Micro Wave
      CASE (sensor_name(sensor_id_hi))
        coef%id_sensor = sensor_id_hi! High resolution
      CASE (sensor_name(sensor_id_po))
        coef%id_sensor = sensor_id_po! Polarimetric
      CASE DEFAULT
        coef%id_sensor = sensor_id_ir
      END SELECT

      READ (file_lu,  * , iostat=ERR)coef%id_comp_lvl

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

! Error if the compatibility version of the coefficient file
! is not in the range defined by the constant module

      IF (coef%id_comp_lvl < version_compatible_min .OR. coef%id_comp_lvl > version_compatible_max) THEN
        err = errorstatus_fatal

        THROWM( ERR .NE. 0, "Version of coefficient file is incompatible with RTTOV library")

      ENDIF

      READ (file_lu, '(a)', iostat=ERR)coef%id_creation

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_lu,  * , iostat=ERR)coef%id_creation_date

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

    CASE ('LINE-BY-LINE')
      CALL rttov_skipcommentline(file_lu, ERR)

      LoopLBL : DO i = 1, Size(coef%line_by_line)
        READ (file_lu, '(a)', iostat=ERR)coef%line_by_line(i)

        DO j = 1, Len(coef%line_by_line(i))

          SELECT CASE (coef%line_by_line(i)(j:j))
          CASE ('!')
            EXIT LoopLBL
          CASE (' ')
            CYCLE
          CASE DEFAULT
            EXIT
          END SELECT

        ENDDO

      ENDDO LoopLBL

      IF (i .LE. Size(coef%line_by_line)) coef%line_by_line(i) = 'xxxx'
    CASE ('FAST_MODEL_VARIABLES')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

! fast model variables definition
      READ (file_lu,  * , iostat=ERR)coef%fmv_model_def

      THROWM( ERR .NE. 0, 'io status while reading section '//section)


! fast model variables version
      READ (file_lu,  * , iostat=ERR)coef%fmv_model_ver

      THROWM(err.ne.0,"io status while reading section "//section)

! number of channels stored
      READ (file_lu,  * , iostat=ERR)coef%fmv_chn

      THROWM(err.ne.0,"io status while reading section "//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests
      file_channels = coef%fmv_chn

      IF (.NOT. all_channels) THEN
        coef%fmv_chn = Size(channels)
      ENDIF

! number of gases in file
      READ (file_lu,  * , iostat=ERR)coef%fmv_gas

      THROWM(err.ne.0,"io status while reading section "//section)

! Flag to check compatibility with principal component regression file
      coef%id_comp_pc = 0_jpim

      READ (file_lu,  * , iostat=ERR)coef%id_comp_pc

      THROWM(err.ne.0,"io status while reading section "//section)

! Flag requesting full Zeeman code where appropriate
      coef%IncZeeman = .FALSE.

! Include Zeeman flag
      READ (file_lu,  * , iostat=ERR)IncZeeman

      THROWM(err.ne.0,"io status while reading section "//section)

      IF (IncZeeman == 0) coef%IncZeeman = .FALSE.
      IF (incZeeman == 1) coef%IncZeeman = .TRUE.

! Flag indicating that profile is to include ToA (new for RTTOV-9 and later files)
! now separate from IncZeeman. New default is .True.
! allocate arrays of FAST_MODEL_VARIABLES section
      ALLOCATE (coef%fmv_gas_id(coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_gas_id")

      ALLOCATE (coef%fmv_gas_pos(ngases_max), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_gas_pos")

      ALLOCATE (coef%fmv_var(coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_var")

      ALLOCATE (coef%fmv_lvl(coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_lvl")

      ALLOCATE (coef%fmv_coe(coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_coe")

      coef%fmv_gas_id(:)  = 0_jpim
      coef%fmv_gas_pos(:) = 0_jpim
      coef%fmv_var(:)     = 0_jpim
      coef%fmv_lvl(:)     = 0_jpim
      coef%fmv_coe(:)     = 0_jpim

      DO n = 1, coef%fmv_gas
! gas id. number i gas_id list (fmv_gas)
        READ (file_lu, '(a)', iostat=ERR)gas_Type

        THROWM(err.ne.0,"io status while reading section "//section)

        CALL Rttov_deletecomment(gas_Type)
        found = .FALSE.

        DO i = 1, ngases_max

          IF (rttov_cmpuc(gas_Type, gas_name(i))) THEN
            coef%fmv_gas_id(n) = i
            found              = .TRUE.
            EXIT
          ENDIF

        ENDDO


        IF (.NOT. found) THEN
          err = errorstatus_fatal

          THROWM(err.ne.0,Trim(gas_Type))

          RETURN
        ENDIF

! store also the indice of this gas in the
! identification list
! so fmv_gas_pos(1) will give position of MxG in the file
        coef%fmv_gas_pos(coef%fmv_gas_id(n)) = n
! number of variables/predictors by gaz
! number of levels(pres/absorber) by gaz (fmv_gas

        READ (file_lu,  * , iostat=ERR)coef%fmv_var(n), coef%fmv_coe(n), coef%fmv_lvl(n)

        THROWM(err.ne.0,"io status while reading section "//section)

! Transfer information to some "classical" variables
! with more common names
! Note that the number of levels is taken from the Mixed Gases line

        SELECT CASE (coef%fmv_gas_id(n))
        CASE (gas_id_mixed)
          coef%nmixed  = coef%fmv_var(n)
          coef%ncmixed = coef%fmv_coe(n)
          coef%nlevels = coef%fmv_lvl(n)
        CASE (gas_id_watervapour)
          coef%nwater  = coef%fmv_var(n)
          coef%ncwater = coef%fmv_coe(n)
        CASE (gas_id_ozone)
          coef%nozone  = coef%fmv_var(n)
          coef%ncozone = coef%fmv_coe(n)
        CASE (gas_id_wvcont)
          coef%nwvcont  = coef%fmv_var(n)
          coef%ncwvcont = coef%fmv_coe(n)
        CASE (gas_id_co2)
          coef%nco2  = coef%fmv_var(n)
          coef%ncco2 = coef%fmv_coe(n)
        CASE (gas_id_n2o)
          coef%nn2o  = coef%fmv_var(n)
          coef%ncn2o = coef%fmv_coe(n)
        CASE (gas_id_co)
          coef%nco  = coef%fmv_var(n)
          coef%ncco = coef%fmv_coe(n)
        CASE (gas_id_ch4)
          coef%nch4  = coef%fmv_var(n)
          coef%ncch4 = coef%fmv_coe(n)
        END SELECT

      ENDDO

      coef%nlayers = coef%nlevels - 1
! Initialise the gaz units array with defaults values
! (specific concetration kg/kg)
      ALLOCATE (coef%gaz_units(coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of gaz_units")

      coef%gaz_units(:) = gas_unit_specconc
    CASE ('GAS_SPECTRAL_INTERVAL')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO n = 1, coef%fmv_gas
! gas id. number i gas_id list (fmv_gas)
        READ (file_lu, '(a)', iostat=ERR)gas_Type

        THROWM(err.ne.0,"io status while reading section "//section)

        CALL Rttov_deletecomment(gas_Type)
        gas_id = 0

        DO i = 1, ngases_max

          IF (rttov_cmpuc(gas_Type, gas_name(i))) THEN
            gas_id = i
            EXIT
          ENDIF

        ENDDO

        IF (gas_id == 0) WRITE ( * , '(a)')'Error: gas type '//Trim(gas_Type)//' not recognised'
        READ (file_lu,  * , iostat=ERR)nbint

        THROWM(err.ne.0,"io status while reading section "//section)

! Transfer information to some "classical" variables
! with more common names

        SELECT CASE (gas_id)
        CASE (gas_id_mixed)
          coef%nintmixed = nbint
          ALLOCATE (coef%mixedgasint(2, coef%nintmixed), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of mixedgasint")


          DO i = 1, coef%nintmixed
            READ (file_lu,  * , iostat=ERR)coef%mixedgasint(1, i), coef%mixedgasint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_watervapour)
          coef%nintwater = nbint
          ALLOCATE (coef%watervapourint(2, coef%nintwater), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of watervapourint")


          DO i = 1, coef%nintwater
            READ (file_lu,  * , iostat=ERR)coef%watervapourint(1, i), coef%watervapourint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_ozone)
          coef%nintozone = nbint
          ALLOCATE (coef%ozoneint(2, coef%nintozone), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of ozoneint")


          DO i = 1, coef%nintozone
            READ (file_lu,  * , iostat=ERR)coef%ozoneint(1, i), coef%ozoneint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_wvcont)
          coef%nintwvcont = nbint
          ALLOCATE (coef%wvcontint(2, coef%nintwvcont), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of wvcontint")


          DO i = 1, coef%nintwvcont
            READ (file_lu,  * , iostat=ERR)coef%wvcontint(1, i), coef%wvcontint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_co2)
          coef%nintco2 = nbint
          ALLOCATE (coef%co2int(2, coef%nintco2), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of co2int")


          DO i = 1, coef%nintco2
            READ (file_lu,  * , iostat=ERR)coef%co2int(1, i), coef%co2int(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_n2o)
          coef%nintn2o = nbint
          ALLOCATE (coef%n2oint(2, coef%nintn2o), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of n2oint")


          DO i = 1, coef%nintn2o
            READ (file_lu,  * , iostat=ERR)coef%n2oint(1, i), coef%n2oint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_co)
          coef%nintco = nbint
          ALLOCATE (coef%coint(2, coef%nintco), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of coint")


          DO i = 1, coef%nintco
            READ (file_lu,  * , iostat=ERR)coef%coint(1, i), coef%coint(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        CASE (gas_id_ch4)
          coef%nintch4 = nbint
          ALLOCATE (coef%ch4int(2, coef%nintch4), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of ch4int")


          DO i = 1, coef%nintch4
            READ (file_lu,  * , iostat=ERR)coef%ch4int(1, i), coef%ch4int(2, i)

            THROWM(err.ne.0,"io status while reading section "//section)

          ENDDO

        END SELECT

      ENDDO

    CASE ('FILTER_FUNCTIONS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

!allocate FILTER_FUNCTIONS section  array size is fmv_chn
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


      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
          READ (file_lu,  * , iostat=ERR)     &
            & coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv0")

        ALLOCATE (ivalues1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv1")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v1")

        ALLOCATE (values2(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v2")

        ALLOCATE (values3(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v3")


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), ivalues1(i), values0(i), values1(i), values2(i), values3(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

        coef%ff_ori_chn(:) = ivalues0(channels(:))
        coef%ff_val_chn(:) = ivalues1(channels(:))
        coef%ff_cwn(:)     = values0(channels(:))
        coef%ff_bco(:)     = values1(channels(:))
        coef%ff_bcs(:)     = values2(channels(:))
        coef%ff_gam(:)     = values3(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v2")

        DEALLOCATE (values3, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v3")

      ENDIF

!       coef % ff_gam(:)=1.
    CASE ('TRANSMITTANCE_TRESHOLD')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

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

        DO i = 1, coef%fmv_chn
! chan number
! validity of channel
! central wave number
! transmittance treshold
! transmittance value
          READ (file_lu,  * , iostat=ERR)     &
            & coef%tt_chn(i), coef%tt_val_chn(i), coef%tt_cwn(i), coef%tt_a0(i), coef%tt_a1(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv0")

        ALLOCATE (ivalues1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv1")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v1")

        ALLOCATE (values2(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v2")


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), ivalues1(i), values0(i), values1(i), values2(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

        coef%tt_chn(:)     = ivalues0(channels(:))
        coef%tt_val_chn(:) = ivalues1(channels(:))
        coef%tt_cwn(:)     = values0(channels(:))
        coef%tt_a0(:)      = values1(channels(:))
        coef%tt_a1(:)      = values2(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v2")

      ENDIF

    CASE ('SOLAR_SPECTRUM')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef%ss_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_chn")

      ALLOCATE (coef%ss_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_val_chn")

      ALLOCATE (coef%ss_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_cwn")

      ALLOCATE (coef%ss_solar_spectrum(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_solar_spectrum")


      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
! chan number
! validity of channel
! central wave number
! solar spectrum
          READ (file_lu,  * , iostat=ERR)coef%ss_chn(i), coef%ss_val_chn(i), coef%ss_cwn(i), coef%ss_solar_spectrum(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv0")

        ALLOCATE (ivalues1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of iv1")

        ALLOCATE (values0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of v1")


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), ivalues1(i), values0(i), values1(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

        coef%ss_chn(:)            = ivalues0(channels(:))
        coef%ss_val_chn(:)        = ivalues1(channels(:))
        coef%ss_cwn(:)            = values0(channels(:))
        coef%ss_solar_spectrum(:) = values1(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of iv1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of v1")

      ENDIF

    CASE ('WATER_OPTICAL_CONSTANT')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef%woc_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_chn")

      ALLOCATE (coef%woc_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_cwn")

      ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_ow")

      ALLOCATE (coef%woc_waopc_fw(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_fw")


      IF (all_channels) THEN

        DO i = 1, coef%fmv_chn
! chan number
! central wave number
! ocean water optical constants(real and imaginary part)
! fresh water optical constants(real and imaginary part)
          READ (file_lu,  * , iostat=ERR)coef%woc_chn(i), coef%woc_cwn(i), coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (values1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values_c_0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_0")

        ALLOCATE (values_c_1(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_1")


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), values1(i), values_c_0(i), values_c_1(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

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

    CASE ('WAVE_SPECTRUM')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

      READ (file_lu,  * , iostat=ERR)coef%ws_nomega
      ALLOCATE (coef%ws_k_omega(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_k_omega")

      ALLOCATE (coef%ws_npoint(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_npoint")


      DO i = 1, coef%ws_nomega
        READ (file_lu,  * , iostat=ERR)coef%ws_npoint(i), coef%ws_k_omega(i)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO

    CASE ('FUNDAMENTAL_CONSTANTS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! speed of light (cm/s)
      READ (file_lu,  * , iostat=ERR)coef%fc_speedl

      THROWM(err.ne.0,"io status while reading section "//section)

! first radiation constant (mW/(m2*sr*cm-4))
! second radiation constant (cm*K)
      READ (file_lu,  * , iostat=ERR)coef%fc_planck_c1, coef%fc_planck_c2

      THROWM(err.ne.0,"io status while reading section "//section)

! satellite nominal altitude (km)
      READ (file_lu,  * , iostat=ERR)coef%fc_sat_height

      THROWM(err.ne.0,"io status while reading section "//section)

    CASE ('FASTEM')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! fastem version number
      READ (file_lu,  * , iostat=ERR)coef%fastem_ver

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fastem_polar")

! polarisation of each channel

      IF (all_channels) THEN
        READ (file_lu,  * , iostat=ERR)(coef%fastem_polar(i), i = 1, coef%fmv_chn)
      ELSE
        ALLOCATE (ivalues0(file_channels), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        READ (file_lu,  * , iostat=ERR)(ivalues0(i), i = 1, file_channels)
        coef%fastem_polar(:) = ivalues0(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

      ENDIF


      THROWM(err.ne.0,"io status while reading section "//section)

!-------------------------------------------------------
    CASE ('SSIREM')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! version number
      READ (file_lu,  * , iostat=ERR)coef%ssirem_ver

      THROWM(err.ne.0,"io status while reading section "//section)

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

        DO i = 1, coef%fmv_chn
! original chan number
! constant coef
! first order coef
! second order coef
! 1st exponent on zenith angle
! 2nd exponent on zenith angle
          READ (file_lu,  * , iostat=ERR)coef%ssirem_chn(i), coef%ssirem_a0(i), coef%ssirem_a1(i), coef%ssirem_a2(i),      &
            & coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

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


        DO i = 1, file_channels
          READ (file_lu,  * , iostat=ERR)ivalues0(i), values0(i), values1(i), values2(i), values3(i), values4(i)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

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

!-------------------------------------------------------
    CASE ('GAZ_UNITS')
! the array has already been allocated and initialised
! to specific concentration (kg/kg)
!
! This section needs one input line per gaz
! in the same order as the gaz list defined inside
!
! This is defining the units used for the sections
! REFERENCE_PROFILE and PROFILE_LIMITS
!
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  * , iostat=ERR)coef%gaz_units(n)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO

!-------------------------------------------------------
    CASE ('REFERENCE_PROFILE')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef%ref_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ref_prfl_p")

      ALLOCATE (coef%ref_prfl_t(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ref_prfl_t")

      ALLOCATE (coef%ref_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ref_prfl_mr")


      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)

! units for reference gaz concentration is
! specified in GAZ_UNITS section (default is specific concentration (kg/kg))
!

        DO i = 1, coef%nlevels
          READ (file_lu,  * , iostat=ERR)pres, coef%ref_prfl_t(i, n), coef%ref_prfl_mr(i, n)

          THROWM(err.ne.0,"io status while reading section "//section)

          IF (coef%fmv_gas_id(n) == gas_id_mixed) coef%ref_prfl_p(i) = pres
        ENDDO

      ENDDO

!-------------------------------------------------------
    CASE ('PROFILE_LIMITS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

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


      DO l = 1, coef%nlevels
! pressure  (hPa)       (levels)
! max temperature (K)   (levels)
! min temperature (K)   (levels)
        READ (file_lu,  * , iostat=ERR)coef%lim_prfl_p(l), coef%lim_prfl_tmax(l), coef%lim_prfl_tmin(l)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)


        DO i = 1, coef%nlevels
! max specific concentration (kg/kg) (levels, gases)
! min specific concentration (kg/kg) (levels, gases)
! or
! max volume mixing r (ppmv) (levels, gases)
! min volume mixing r (ppmv) (levels, gases)
! according to
! units specified in GAZ_UNITS section (default is specific concentration (kg/kg))
          READ (file_lu,  * , iostat=ERR)pres, coef%lim_prfl_gmax(i, n), coef%lim_prfl_gmin(i, n)

          THROWM(err.ne.0,"io status while reading section "//section)

        ENDDO

      ENDDO

!-------------------------------------------------------
    CASE ('FAST_COEFFICIENTS', 'COEF_SUB_FILES')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! If section is COEF_SUB_FILES then coefficients for each gaz is stored
! in separate files.
! This possibility could be used to store very large coefficient files
! (large number of channels or gases)
! Section contains 1 line per gaz in the same order as the
! FAST_MODEL_VARIABLES section
! line indicates the filename of the coefficient for that gas
!
! No verification is done on the file.
! header lines starting with "!" sign are skipped
! loop on gases

      DO n = 1, coef%fmv_gas
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)

! read dummy string of gas name or filename of the sub_coefficient file
        READ (file_lu,  * , iostat=ERR)input_string

        THROWM(err.ne.0,"io status while reading section "//section)

! Case of Sub coefficient files
! Open the file and skip the header

        IF (Trim(section) == 'COEF_SUB_FILES') THEN
          file_lu_coef = 0
          for_output   = .FALSE.
          CALL rttov_opencoeff( &
                & err,          &
                & input_string, &
                & file_lu_coef, &
                & for_output)

          IF (err /= 0) THEN
            err = errorstatus_fatal
            WRITE (errMessage, '( "Error opening sub_coef file" )')

            THROWM(err.ne.0,errMessage)

            RETURN
          ENDIF

          CALL rttov_skipcommentline(file_lu_coef, ERR)

          THROWM(err.ne.0,"io status while reading section "//section)

        ELSE
          file_lu_coef = file_lu
        ENDIF


        SELECT CASE (coef%fmv_gas_id(n))
        CASE (gas_id_mixed)
          nvalues = coef%ncmixed
          ALLOCATE (coef%mixedgas(coef%nlayers, coef%fmv_chn, coef%ncmixed), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of mixedgas")


          IF (all_channels) THEN
            coeffsarray => coef%mixedgas
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncmixed), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_watervapour)
          nvalues = coef%ncwater
          ALLOCATE (coef%watervapour(coef%nlayers, coef%fmv_chn, coef%ncwater), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of watervapour")


          IF (all_channels) THEN
            coeffsarray => coef%watervapour
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncwater), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_ozone)
          nvalues = coef%ncozone
          ALLOCATE (coef%ozone(coef%nlayers, coef%fmv_chn, coef%ncozone), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of ozone")


          IF (all_channels) THEN
            coeffsarray => coef%ozone
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncozone), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_wvcont)
          nvalues = coef%ncwvcont
          ALLOCATE (coef%wvcont(coef%nlayers, coef%fmv_chn, coef%ncwvcont), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of wvcont")


          IF (all_channels) THEN
            coeffsarray => coef%wvcont
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncwvcont), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_co2)
          nvalues = coef%ncco2
          ALLOCATE (coef%co2(coef%nlayers, coef%fmv_chn, coef%ncco2), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of co2")


          IF (all_channels) THEN
            coeffsarray => coef%co2
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncco2), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_n2o)
          nvalues = coef%ncn2o
          ALLOCATE (coef%n2o(coef%nlayers, coef%fmv_chn, coef%ncn2o), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of n2o")


          IF (all_channels) THEN
            coeffsarray => coef%n2o
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncn2o), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_co)
          nvalues = coef%ncco
          ALLOCATE (coef%co(coef%nlayers, coef%fmv_chn, coef%ncco), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of co")


          IF (all_channels) THEN
            coeffsarray => coef%co
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncco), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        CASE (gas_id_ch4)
          nvalues = coef%ncch4
          ALLOCATE (coef%ch4(coef%nlayers, coef%fmv_chn, coef%ncch4), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of ch4")


          IF (all_channels) THEN
            coeffsarray => coef%ch4
          ELSE
            ALLOCATE (coeffsarray(coef%nlayers, file_channels, coef%ncch4), STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of coeffsarray")

          ENDIF

        END SELECT

        READ (file_lu_coef,  * , iostat=ERR)     &
          & (((coeffsarray(i, j, k), i = 1, coef%nlayers), j = 1, file_channels), k = 1, nvalues)

        THROWM(err.ne.0,"io status while reading section "//section)


        IF (.NOT. all_channels) THEN

          SELECT CASE (coef%fmv_gas_id(n))
          CASE (gas_id_mixed)
            coef%mixedgas(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_watervapour)
            coef%watervapour(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_ozone)
            coef%ozone(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_wvcont)
            coef%wvcont(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_co2)
            coef%co2(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_n2o)
            coef%n2o(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_co)
            coef%co(:,:,:) = coeffsarray(:, channels(:), :)
          CASE (gas_id_ch4)
            coef%ch4(:,:,:) = coeffsarray(:, channels(:), :)
!
          END SELECT

          DEALLOCATE (coeffsarray, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of coeffsarray")

        ENDIF

! For COEF_SUB_FILES close the intermediate coef file

        IF (Trim(section) == 'COEF_SUB_FILES') THEN
          CLOSE (unit=file_lu_coef)
        ENDIF

      ENDDO

    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_coef
