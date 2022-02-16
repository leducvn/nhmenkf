!
SUBROUTINE rttov_read_ascii_scaercoef( &
            & err,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
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
       & rttov_coef,          &
       & rttov_optpar_ir,     &
       & rttov_coef_scatt_ir, &
       & rttov_coef_pccomp
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
       & version_compatible_min, &
       & version_compatible_max, &
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
  TYPE(rttov_optpar_ir    ), INTENT(INOUT)             :: optp
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT)             :: coef_scatt_ir
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: err           ! return code
!INTF_END
#include "rttov_opencoeff.h"
#include "rttov_errorreport.h"
#include "rttov_skipcommentline.h"
#include "rttov_deletecomment.h"
#include "rttov_cmpuc.h"
#include "rttov_findnextsection.h"
#include "rttov_nullify_coef_scatt_ir.h"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  INTEGER(KIND=jpim) :: pha_channels
  INTEGER(KIND=jpim) :: pha_1           , pha_2
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: nbchannels_solar
  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, j, k, n, nrh, is, icount
  INTEGER(KIND=jpim), ALLOCATABLE :: list_of_channels(:)
! pointers for generic inputs
  REAL   (KIND=jprb), POINTER :: abs_aer_array(:, :   )
  REAL   (KIND=jprb), POINTER :: sca_aer_array(:, :   )
  REAL   (KIND=jprb), POINTER :: bpr_aer_array(:, :   )
  REAL   (KIND=jprb), POINTER :: pha_aer_array(:, :, :)
  CHARACTER(LEN = 36)         :: input_string
  CHARACTER(LEN = 32)         :: aer_comp
  CHARACTER(LEN = lensection) :: section
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


!read the file

  readfile : DO
    CALL rttov_findnextsection(file_lu, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (Trim(section))

!-------------------------------------------------------------------------------
    CASE ('AEROSOLS_COMPONENTS')

      IF (.NOT. all_channels) THEN
        ALLOCATE (list_of_channels(size(channels)))
        list_of_channels = channels
      ELSE
        ALLOCATE (list_of_channels(coef%fmv_chn))
        list_of_channels = (/(i, i = 1, coef%fmv_chn)/)
      ENDIF

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! number of channels for which optical parameters are stored
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_chn

      THROWM(err.ne.0,"io status while reading section "//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_aer_chn is the number of channels that the user requests
      file_channels = coef_scatt_ir%fmv_aer_chn

      IF (.NOT. all_channels) THEN
        coef_scatt_ir%fmv_aer_chn = Size(channels)
      ENDIF

! Number of channels for which phase function values are stored
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_pha_chn

      THROWM(err.ne.0,"io status while reading section "//section)

! index of first channel for which phase function values are available
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_sun_chn

      THROWM(err.ne.0,"io status while reading section "//section)

! number of aerosols components
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_comp

      THROWM(err.ne.0,"io status while reading section "//section)

! number of angles for phase function
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_ph

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef_scatt_ir%fmv_aer_ph_val(coef_scatt_ir%fmv_aer_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val")

      ALLOCATE (coef_scatt_ir%fmv_aer_ph_val_cos(coef_scatt_ir%fmv_aer_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val_cos")

! angles at which values of the phase function are given
      READ (file_lu,  * , iostat=ERR)(coef_scatt_ir%fmv_aer_ph_val(i), i = 1, coef_scatt_ir%fmv_aer_ph)

      DO i = 1, coef_scatt_ir%fmv_aer_ph
        coef_scatt_ir%fmv_aer_ph_val_cos(i) = cos(coef_scatt_ir%fmv_aer_ph_val(i) * deg2rad)
      ENDDO

      coef_scatt_ir%fmv_aer_ph_val_min = 99999.0_JPRB

      DO i = 1, coef_scatt_ir%fmv_aer_ph - 1

        IF (coef_scatt_ir%fmv_aer_ph_val_min > coef_scatt_ir%fmv_aer_ph_val(i + 1) - coef_scatt_ir%fmv_aer_ph_val(i)     &
          & ) THEN
          coef_scatt_ir%fmv_aer_ph_val_min = coef_scatt_ir%fmv_aer_ph_val(i + 1) - coef_scatt_ir%fmv_aer_ph_val(i)
        ENDIF

      ENDDO

      coef_scatt_ir%fmv_aer_ph_val_min = coef_scatt_ir%fmv_aer_ph_val_min / 2.0_JPRB
      icount = coef_scatt_ir%fmv_aer_ph_val(coef_scatt_ir%fmv_aer_ph) / coef_scatt_ir%fmv_aer_ph_val_min
      ALLOCATE (coef_scatt_ir%ifmv_aer_ph_val(icount), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ifmv_aer_ph_val")

      k = 1

      DO i = 1, icount

        IF (coef_scatt_ir%fmv_aer_ph_val(k) >= (i + 1) * coef_scatt_ir%fmv_aer_ph_val_min) THEN
          coef_scatt_ir%ifmv_aer_ph_val(i) = k
        ELSE
          k = k + 1
          coef_scatt_ir%ifmv_aer_ph_val(i) = k
        ENDIF

      ENDDO


      THROWM(err.ne.0,"io status while reading section "//section)

      is = 0

      DO i = 1, size(list_of_channels)

        IF ((list_of_channels(i) >= coef_scatt_ir%fmv_aer_sun_chn) .AND.      &
          & (list_of_channels(i) <= coef_scatt_ir%fmv_aer_sun_chn + coef_scatt_ir%fmv_aer_pha_chn - 1)) THEN
          is = is + 1
        ENDIF

      ENDDO

!--------------a dirty fix for imagers with chans ordered in wavelength-----
      IF (coef_scatt_ir%fmv_aer_sun_chn == 1) is = 1
!---------------------------------------------------------------------------
      pha_channels = coef_scatt_ir%fmv_aer_pha_chn

      IF (is /= coef_scatt_ir%fmv_aer_pha_chn) THEN
        coef_scatt_ir%fmv_aer_pha_chn = is
      ENDIF

! allocate arrays of FAST_MODEL_VARIABLES section
      ALLOCATE (coef_scatt_ir%fmv_aer_rh(coef_scatt_ir%fmv_aer_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_aer_rh")

      IF (is .lt. 1) is = 1
      ALLOCATE (coef_scatt_ir%channels_solar(is), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of channels_solar")

      ALLOCATE (optp%optpaer(coef_scatt_ir%fmv_aer_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_aer_comp")

      is = 0
      coef_scatt_ir%fmv_aer_pha_ioff = 0

      DO i = 1, size(list_of_channels)

        IF (list_of_channels(i) >= coef_scatt_ir%fmv_aer_sun_chn) THEN
          coef_scatt_ir%fmv_aer_pha_ioff = i
          EXIT
        ENDIF

      ENDDO

      is = 0

      IF (coef_scatt_ir%fmv_aer_sun_chn > 1) THEN

        DO i = 1, size(list_of_channels)

          IF ((list_of_channels(i) >= coef_scatt_ir%fmv_aer_sun_chn) .AND.      &
            & (list_of_channels(i) <= coef_scatt_ir%fmv_aer_sun_chn + pha_channels - 1)) THEN
            is = is + 1
            coef_scatt_ir%channels_solar(is) = list_of_channels(i)
          ENDIF

        ENDDO
        IF (is < 1) coef_scatt_ir%channels_solar(1) = 1

      ELSE
        coef_scatt_ir%channels_solar(1) = 1
      ENDIF
      nbchannels_solar = is


      DO n = 1, coef_scatt_ir%fmv_aer_comp
        CALL rttov_nullify_coef_scatt_ir (optp%optpaer(n))
! aerosol id. number i aerosol_id list (fmv_gas)
        READ (file_lu, '(a)', iostat=ERR)aer_comp

        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_aer_rh(n)

        THROWM(err.ne.0,"io status while reading section "//section)

        ALLOCATE (optp%optpaer(n)%fmv_aer_rh_val(coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpaer(n)% fmv_aer_rh_val")

        READ (file_lu,  * , iostat=ERR)(optp%optpaer(n)%fmv_aer_rh_val(i), i = 1, coef_scatt_ir%fmv_aer_rh(n))

        THROWM(err.ne.0,"io status while reading section "//section)

! Transfer information to some "classical" variables
! with more common names
      ENDDO

      DEALLOCATE (list_of_channels)
    CASE ('AEROSOLS_PARAMETERS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

      pha_1 = coef_scatt_ir%fmv_aer_sun_chn
      pha_2 = coef_scatt_ir%fmv_aer_sun_chn + pha_channels - 1
!        print*,pha_1,pha_2
! loop on aerosol components

      DO n = 1, coef_scatt_ir%fmv_aer_comp
        ALLOCATE (optp%optpaer(n)%abs(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of fmv_aer_rh(n)")

        ALLOCATE (optp%optpaer(n)%sca(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of fmv_aer_rh(n)")

        ALLOCATE (optp%optpaer(n)%bpr(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of fmv_aer_rh(n)")

        ALLOCATE (                                                                                                    &
          & optp%optpaer(n)%pha(coef_scatt_ir%fmv_aer_pha_chn, coef_scatt_ir%fmv_aer_rh(n), coef_scatt_ir%fmv_aer_ph) &
          & , STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of fmv_aer_rh(n)")


        IF (all_channels) THEN
          abs_aer_array => optp%optpaer(n)%abs
          sca_aer_array => optp%optpaer(n)%sca
          bpr_aer_array => optp%optpaer(n)%bpr
          pha_aer_array => optp%optpaer(n)%pha
        ELSE
          ALLOCATE (abs_aer_array(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of abs_aer_array")

          ALLOCATE (sca_aer_array(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of sca_aer_array")

          ALLOCATE (bpr_aer_array(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of bpr_aer_array")

          ALLOCATE (pha_aer_array(pha_1:pha_2, coef_scatt_ir%fmv_aer_rh(n), coef_scatt_ir%fmv_aer_ph), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pha_aer_array")

        ENDIF


        DO nrh = 1, coef_scatt_ir%fmv_aer_rh(n)
          CALL rttov_skipcommentline(file_lu, ERR)

          THROWM(err.ne.0,"io status while reading section "//section)

! read dummy string of gas name or filename of the sub_coefficient file
          READ (file_lu,  * , iostat=ERR)input_string

          THROWM(err.ne.0,"io status while reading section "//section)

!write(*,*)'Optical parameters available for: ',input_string
          READ (file_lu,  * , iostat=ERR)(abs_aer_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)

          READ (file_lu,  * , iostat=ERR)(sca_aer_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)

          READ (file_lu,  * , iostat=ERR)(bpr_aer_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)


          IF (pha_channels > 0) THEN

            IF (all_channels) THEN
              READ (file_lu,  * , iostat=ERR)     &
                & ((pha_aer_array(i, nrh, j), j = 1, coef_scatt_ir%fmv_aer_ph), i = 1, coef_scatt_ir%fmv_aer_pha_chn)
            ELSE
              READ (file_lu,  * , iostat=ERR)     &
                & ((pha_aer_array(i, nrh, j), j = 1, coef_scatt_ir%fmv_aer_ph), i = pha_1, pha_2)
            ENDIF

          ENDIF


          THROWM(err.ne.0,"io status while reading pha_aer_array section "//section)

!abs_aer_array=0.
!sca_aer_array=0.
!bpr_aer_array=0.
        ENDDO


        IF (.NOT. all_channels) THEN
          optp%optpaer(n)%abs(:,:)   = abs_aer_array(channels(:), :)
          optp%optpaer(n)%sca(:,:)   = sca_aer_array(channels(:), :)
          optp%optpaer(n)%bpr(:,:)   = bpr_aer_array(channels(:), :)
          IF ( nbchannels_solar  > 0 ) THEN
            optp%optpaer(n)%pha(:,:,:) = pha_aer_array(coef_scatt_ir%channels_solar(:), :, :)
          ENDIF
          DEALLOCATE (abs_aer_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of abs_aer_array")

          DEALLOCATE (sca_aer_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of sca_aer_array")

          DEALLOCATE (bpr_aer_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of bpr_aer_array")

          DEALLOCATE (pha_aer_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of pha_aer_array")

        ENDIF

      ENDDO


    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_scaercoef
