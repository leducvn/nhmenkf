SUBROUTINE rttov_write_ascii_sccldcoef ( &
            & err,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
            & file_id)
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
  USE rttov_types, ONLY :  &
       & rttov_coef,          &
       & rttov_optpar_ir,     &
       & rttov_coef_scatt_ir
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
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
  INTEGER(KIND=jpim)       , INTENT(IN)              :: file_id      ! file logical unit number
  TYPE(rttov_coef         ), INTENT(IN)              :: coef         ! coefficients
  TYPE(rttov_optpar_ir    ), INTENT(IN)              :: optp
  TYPE(rttov_coef_scatt_ir), INTENT(IN)              :: coef_scatt_ir
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(OUT)             :: err          ! return code
!INTF_END
#include "rttov_errorreport.h"
! local scalars
  INTEGER(KIND=jpim)  :: i, n, nr, nrh
  CHARACTER(LEN = 32) :: section
  CHARACTER(LEN = 80) :: errMessage
!- End of header --------------------------------------------------------
  TRY

  
!ASCII file
  WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
  INFO(errMessage)

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV coefficient file '//Trim(coef%id_Common_name)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! automatic creation by subroutine Rttov_writecoef '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

! COEF structure (V10)

  section = 'WATERCLOUD_TYPES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_chn, &
      &'! number of channels for which optical parameters are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_pha_chn,&
      &'! Number of channels for which phase function values are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_pha_ioff,&
      &'! index of channel for solar term'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_comp, &
      &'! number of water cloud types'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_ph, &
      &'! number of angles for phase function for water cloud types'

  THROW(err.ne.0)

  WRITE (file_id, '(10f7.2)', iostat=err)coef_scatt_ir%fmv_wcl_ph_val

  THROW(err.ne.0)


  DO n = 1, coef_scatt_ir%fmv_wcl_comp
    WRITE (file_id, '(a)', iostat=err)' cloud   ! default name for rttov_writecoef'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_rh(n),&
        &'!RH values for which parameters are available'

    THROW(err.ne.0)

    WRITE (file_id, '(10f7.2)', iostat=err)optp%optpwcl(n)%fmv_wcl_rh_val

    THROW(err.ne.0)

    WRITE (file_id, '(1x,f12.6,t20,a)', iostat=err)coef_scatt_ir%confac(n), &
        & '!Conversion from LWC to particle density'

    THROW(err.ne.0)

  ENDDO

  section = 'WATERCLOUD_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO n = 1, coef_scatt_ir%fmv_wcl_comp

    DO nrh = 1, coef_scatt_ir%fmv_wcl_rh(n)
      WRITE (file_id, '(a)', iostat=err)' cloud   ! default name for rttov_writecoef'

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%abs(:, nrh)

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%sca(:, nrh)

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%bpr(:, nrh)

      THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!& ((optp%optpwcl(n) % pha(i,nrh,j),&
!& j = 1, coef_scatt_ir % fmv_wcl_ph),i=1,coef_scatt_ir % fmv_wcl_pha_chn)

      IF (coef_scatt_ir%fmv_wcl_pha_chn > 0) THEN

        DO i = 1, coef_scatt_ir%fmv_wcl_pha_chn
          WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%pha(i, nrh, :)

          THROW(err.ne.0)

        ENDDO

      ENDIF

    ENDDO
  
  ENDDO

  section = 'ICECLOUD_TYPES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_chn, &
      & '! number of channels for which regression coefficients are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_pha_chn,&
      & '! Number of channels for which phase function values are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_pha_ioff,&
      & '! index of channel for solar term'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nabs,&
      & '! number of coefficients used in the regression for absorption optical depth'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nsca, &
      & '! number of coefficients used in the regression for scattering optical depth'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nbpr,&
      & '! number of coefficients used in the regression for backscattering parameter' 

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_comp, &
      & '! number of size distributions used in the regression' 

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_ishp,&
      & '! number of ice crystal shapes for which parameters are available'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_ph,&
      & '! number of angles for phase function for ice clouds' 

  THROW(err.ne.0)

  WRITE (file_id, '(10f7.2)', iostat=err)coef_scatt_ir%fmv_icl_ph_val

  THROW(err.ne.0)

  section = 'HEXAGONAL_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Effective diameter for each size distribution'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(f10.4)', iostat=err)coef_scatt_ir%fmv_icl_dg(:, 1)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Regression coefficients for ice clouds'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%abs(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%sca(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%bpr(i, :)

    THROW(err.ne.0)

  ENDDO


  DO nr = 1, coef_scatt_ir%fmv_icl_comp
    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Phase function values'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' Hexagonal   ! default name for rttov_writecoef'

    THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!  & ((optp%optpicl(1) % pha(i,nr,j),&
!  & j = 1, coef_scatt_ir % fmv_icl_ph),i=1,coef_scatt_ir % fmv_icl_pha_chn)

    IF (coef_scatt_ir%fmv_icl_pha_chn > 0) THEN

      DO i = 1, coef_scatt_ir%fmv_icl_pha_chn
        WRITE (file_id, '(5e16.8)', iostat=err)optp%optpicl(1)%pha(i, nr, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDDO

  section = 'AGGREGATE_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Effective diameter for each size distribution'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(f10.4)', iostat=err)coef_scatt_ir%fmv_icl_dg(:, 2)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Regression coefficients for ice clouds'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%abs(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%sca(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%bpr(i, :)

    THROW(err.ne.0)

  ENDDO


  DO nr = 1, coef_scatt_ir%fmv_icl_comp
    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Phase function values'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' aggregate   ! default name for rttov_writecoef'

    THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!  & ((optp%optpicl(2) % pha(i,nr,j),&
!  & j = 1, coef_scatt_ir % fmv_icl_ph),i=1,coef_scatt_ir % fmv_icl_pha_chn)

    IF (coef_scatt_ir%fmv_icl_pha_chn > 0) THEN

      DO i = 1, coef_scatt_ir%fmv_icl_pha_chn
        WRITE (file_id, '(5e16.8)', iostat=err)optp%optpicl(2)%pha(i, nr, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDDO

  INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
