SUBROUTINE rttov_write_ascii_pccoef( &
            & err,         &
            & coef_pccomp, &
            & file_lu)
! Description:
!   Write the PC coefficient structure
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
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)     , INTENT(IN)  :: file_lu    ! file logical unit number
  TYPE(rttov_coef_pccomp), INTENT(IN)  :: coef_pccomp
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)     , INTENT(OUT) :: err        ! return code
!INTF_END
#include "rttov_errorreport.h"
  
  CHARACTER(LEN=32)   :: section
  INTEGER(KIND=jpim)  :: n, i, j
!
  TRY

     
    section = 'PRINCOMP_PREDICTORS'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_comp_pc, '!Principal components coefficient file version number'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_sets, '!Number of regression sets'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

! loop on predictor sets

    DO n = 1, coef_pccomp%fmv_pc_sets

      WRITE (file_lu,  '(i5,1x,a,i5)', iostat=ERR)coef_pccomp%pcreg(n)%fmv_pc_npred, '!Number of predictors in set', n
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

      WRITE (file_lu, '(a,i3)', iostat = err) '!Channel indices for set', n
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

      WRITE (file_lu,  '(10i6)', iostat=ERR)(coef_pccomp%pcreg(n)%predictindex(i), i = 1, coef_pccomp%pcreg(n)%fmv_pc_npred)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    section = 'PRINCOMP_EIGENVECTORS'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_mnum, '!Maximum number of eigenvectors'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_nchn, '!Number of channels for which eigenvectors are available'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests

    DO n = 1, coef_pccomp%fmv_pc_mnum

      WRITE (file_lu,  '(5e23.14)', iostat=ERR)(coef_pccomp%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    section = 'PRINCOMP_COEFFICIENTS'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)


    DO n = 1, coef_pccomp%fmv_pc_sets

      WRITE (file_lu, '(a,1x,i5)', iostat = err) '! Predictor set ', n
      THROW(err.ne.0)

      DO j = 1, coef_pccomp%fmv_pc_mnum

        WRITE (file_lu, '(5e23.14)', iostat=ERR)     &
          & (coef_pccomp%pcreg(n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(n)%fmv_pc_npred)
        THROWM( ERR .NE. 0, 'io status while writing section '//section)

      ENDDO

    ENDDO

    section = 'EMISSIVITY_COEFFICIENTS'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    WRITE (file_lu, "(a)", iostat = err) '! Emissivity coefficients for RTIASI routine'
    THROW(err.ne.0)


    WRITE (file_lu, '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_nche, '!Number of channels for which coefficients are available'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests


    DO i = 1, coef_pccomp%fmv_pc_nche

      WRITE (file_lu,  '(i5,5e17.9,/,5x,4e17.9)', iostat=ERR)&
        & coef_pccomp%emiss_chn(i), coef_pccomp%emiss_c1(i), coef_pccomp%emiss_c2(i),                                 &
        & coef_pccomp%emiss_c3(i), coef_pccomp%emiss_c4(i), coef_pccomp%emiss_c5(i), coef_pccomp%emiss_c6(i),         &
        & coef_pccomp%emiss_c7(i), coef_pccomp%emiss_c8(i), coef_pccomp%emiss_c9(i)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    section = 'PC_REFERENCE_PROFILE'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_gas, '! Number of gases for which reference profiles must be used'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(i5,1x,a)', iostat=ERR)coef_pccomp%fmv_pc_nlev, '! Number of levels'
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    DO n = 1, coef_pccomp%fmv_pc_gas

      WRITE (file_lu,  '("! ",i5)', iostat=ERR) n
      THROWM( ERR .NE. 0, 'io status while writing section '//section)
      
      DO i = 1, coef_pccomp%fmv_pc_nlev

        WRITE (file_lu,  '(2g26.16)', iostat=ERR)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i, n)
        THROWM( ERR .NE. 0, 'io status while writing section '//section)

      ENDDO

    ENDDO

    section = 'PC_PROFILE_LIMITS'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    WRITE (file_lu, "(a)", iostat = err) '! Ref.pressure (hPa)'
    THROW(err.ne.0)
    WRITE (file_lu, "(a)", iostat = err) '! Temperature Max and Min [K]'
    THROW(err.ne.0)

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu,  '(3g26.16)', iostat=ERR)     &
        & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_tmin(i), coef_pccomp%lim_pc_prfl_tmax(i)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    WRITE (file_lu, "(a)", iostat = err) '! Ref.pressure (hPa)'
    THROW(err.ne.0)
    WRITE (file_lu, "(a)", iostat = err) '! H2O Max and Min [ppmv]'
    THROW(err.ne.0)

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu, '(3g26.16)', iostat=ERR)     &
        & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_qmin(i), coef_pccomp%lim_pc_prfl_qmax(i)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    WRITE (file_lu, "(a)", iostat = err) '! Ref.pressure (hPa)'
    THROW(err.ne.0)
    WRITE (file_lu, "(a)", iostat = err) '! O3 Max and Min [ppmv]'
    THROW(err.ne.0)

    DO i = 1, coef_pccomp%fmv_pc_nlev

      WRITE (file_lu,  '(3g26.16)', iostat=ERR)     &
        & coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_ozmin(i), coef_pccomp%lim_pc_prfl_ozmax(i)
      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

   
    WRITE (file_lu,  '(a)', iostat = err) '! -------------------------------------------------------'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Surface pressure'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  Max and Min [hPa]'
    THROW(err.ne.0)

    WRITE (file_lu,  '(2g26.16)', iostat=ERR)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(a)', iostat = err) '! -------------------------------------------------------'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Surface temperature'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  Max and Min [K]'
    THROW(err.ne.0)

    WRITE (file_lu,  '(2g26.16)', iostat=ERR)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(a)', iostat = err) '! -------------------------------------------------------'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Skin temperature'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  Max and Min [K]'
    THROW(err.ne.0)

    WRITE (file_lu,  '(2g26.16)', iostat=ERR)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax
    THROWM( ERR .NE. 0, 'io status while writing section '//section)

    WRITE (file_lu,  '(a)', iostat = err) '! -------------------------------------------------------'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS 10m wind speed'
    THROW(err.ne.0)
    WRITE (file_lu,  '(a)', iostat = err) '!  Max and Min [m/s]'
    THROW(err.ne.0)

    WRITE (file_lu,  '(2g26.16)', iostat=ERR)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax
    THROWM( ERR .NE. 0, 'io status while writing section '//section)


    section = 'INSTRUMENT_NOISE'
    WRITE (file_lu, "(a)", iostat = err) Trim(section)
    THROW(err.ne.0)

    DO n = 1, coef_pccomp%fmv_pc_nchn

      WRITE (file_lu,  '(i8,4g26.16)', iostat=ERR)coef_pccomp%ff_ori_chn_in(n), coef_pccomp%ff_cwn_in(n),      &
        & coef_pccomp%ff_bco_in(n), coef_pccomp%ff_bcs_in(n), coef_pccomp%noise_in(n)

      THROWM( ERR .NE. 0, 'io status while writing section '//section)

    ENDDO

    WRITE(file_lu, '(a)', iostat = err) 'END'
    THROW(err.ne.0)

    INFO("end of write coefficient")
    
  CATCH
END SUBROUTINE 
