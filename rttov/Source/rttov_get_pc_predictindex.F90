!
SUBROUTINE rttov_get_pc_predictindex( &
            & err,           &
            & opts,          &
            & predictindex,  &
            & form_pccoef,   &
            & file_pccoef,   &
            & file_id_pccoef,& 
            & instrument)
! Description:
!   Read channel indices of PC predictor set
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
  USE rttov_types, ONLY : rttov_options

  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :    & 
       & rttov_magic_string, &
       & rttov_magic_number, &
       & lensection
  USE parkind1, ONLY : jprb, jplm
  use rttov_coef_io_mod, Only : getlun, closelun
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
Integer(Kind=jpim), Intent(out) :: err
Type(rttov_options),Intent(in)  :: opts
Integer (Kind=jpim), Pointer    ::predictindex(:)
Character(len=*),   Intent(in), Optional :: form_pccoef
Character(len=*),   Intent(in), Optional :: file_pccoef
Integer(Kind=jpim), Intent(in), Optional :: file_id_pccoef
Integer(Kind=jpim), Intent(in), Optional :: instrument(3)

!INTF_END
#include "rttov_errorreport.h"
#include "rttov_skipcommentline.h"
#include "rttov_deletecomment.h"
#include "rttov_cmpuc.h"
#include "rttov_findnextsection.h"

! Local Scalars:


  Integer(Kind=jpim) :: file_id_pccoef1
  character(len=32)  :: form_pccoef1

  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, n

  CHARACTER(LEN = 16) :: bin_check_string
  REAL(KIND=jprb)     :: bin_check_number
  REAL(KIND=jprb)     :: bin_check_value

  INTEGER(KIND=jpim) ::fmv_pc_comp_pc
  INTEGER(KIND=jpim) ::fmv_pc_sets
  INTEGER(KIND=jpim) ::fmv_pc_npred


  CHARACTER(LEN = lensection) :: section
!- End of header --------------------------------------------------------
  TRY

  call getlun(err, file_id_pccoef1, form_pccoef1, .false._jplm, "pccoef", file_id_pccoef, file_pccoef, form_pccoef, &
    instrument )
  THROW(err.ne.0)

!read the file
  if (rttov_cmpuc(form_pccoef1,"formatted")) then
  readascii : DO
    CALL rttov_findnextsection(file_id_pccoef1, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (Trim(section))
    CASE ('PRINCOMP_PREDICTORS')
      CALL rttov_skipcommentline(file_id_pccoef1, ERR)

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_id_pccoef1,  * , iostat=ERR)fmv_pc_comp_pc

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

      READ (file_id_pccoef1,  * , iostat=ERR)fmv_pc_sets

      THROWM( ERR .NE. 0, 'io status while reading section '//section)

! loop on predictor sets

      DO n = 1, fmv_pc_sets

        READ (file_id_pccoef1,  * , iostat=ERR) fmv_pc_npred

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

        CALL rttov_skipcommentline(file_id_pccoef1, ERR)

        THROWM( ERR .NE. 0, 'io status while reading section '//section)

        ALLOCATE (predictindex(fmv_pc_npred), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of predictindex arrays")

        READ (file_id_pccoef1,  * , iostat=ERR)(predictindex(i), i = 1, fmv_pc_npred)
        THROW(err.ne.0)

        if ( opts%ipcreg == n ) EXIT readascii

        DEALLOCATE (predictindex, STAT = ERR)
        THROW(err.ne.0)

      ENDDO
      EXIT readascii

    CASE ('END')
      EXIT readascii
    CASE DEFAULT
      CYCLE readascii
    END SELECT

  ENDDO readascii

  Else

    READ (file_id_pccoef1, iostat=ERR)bin_check_string, bin_check_number
    THROW( ERR .NE. 0)

    ! Verification of header string
    IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal
    THROWM( ERR .NE. 0,'Wrong header string in file')

    ! Verification of single/double precision using a 5 digit number
    ! with exponent 12, which is always Ok for single precision
    bin_check_value = 1._JPRB - abs(bin_check_number - rttov_magic_number)
    IF (bin_check_value > 1.01_JPRB .OR. bin_check_value < 0.99_JPRB) err = errorstatus_fatal

    THROWM( ERR .NE. 0,'File created with a different real precision (R4<->R8)')


    READ (file_id_pccoef1, iostat=ERR)fmv_pc_comp_pc

    THROWM( ERR .NE. 0, "reading coef_pccomp % fmv_pc_comp_pc")

    READ (file_id_pccoef1, iostat=ERR) fmv_pc_sets

    THROWM( ERR .NE. 0, "reading coef_pccomp % fmv_pc_sets")

    readbin: DO n = 1, fmv_pc_sets
      READ (file_id_pccoef1, iostat=ERR) fmv_pc_npred

      THROWM(err.ne.0,"reading fmv_pc_npred")

      ALLOCATE ( predictindex(fmv_pc_npred), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of predictindex")

      READ (file_id_pccoef1, iostat=ERR)(predictindex(i), i = 1, fmv_pc_npred)

      THROWM(err.ne.0,"reading predictindex")

      if ( opts%ipcreg == n ) Exit readbin

      DEALLOCATE (predictindex, STAT = ERR)
      THROW(err.ne.0)
    ENDDO readbin

  EndIf

  call closelun(err, file_id_pccoef1, file_id_pccoef)
  THROW(err.ne.0)

  If(.not. Associated( predictindex ) ) ERR = errorstatus_fatal
  THROWM(err.ne.0,"PC regression set not found")

  CATCH
END SUBROUTINE
