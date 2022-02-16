!+ routine to read Mie coeficient file
!
Subroutine rttov_readscattcoeffs  (&
      & errorstatus,   &! out
      & coef_rttov,    &! in
      & coef_scatt,    &! inout
      & kmyproc,       &! my proc
      & knproc,        &! number of processors
      & kioproc,       &! proc for io
      & file_id       ) ! in  

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
  ! Description:
  ! to initialise Mie look-up table
  !
  ! Method:
  !
  ! Current code owner: saf nwp
  !
  ! History:
  ! version   date        comment
  ! -------   ----        -------
  !   1.0    09/2002      RTTOV7 compatible  (ECMWF)
  !   1.1    05/2003      RTTOV7.3 compatible (ECMWF)
  !   1.2    10/2004      Change stop to return (J Cameron)
  !   1.3    10/2004      Make file_id optional in analogy with rttov_readcoeffs (J Cameron)
  !   1.4    11/2007      Merge with IFS version for RTTOV9 (A Geer)  
  !   1.5    03/2010      Stop adding 1 to nhydro (A Geer)

  ! Imported Type Definitions:
  Use rttov_types, Only : &
       & rttov_coef, &
       & rttov_scatt_coef 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

  Use rttov_const, Only :   &
       & inst_name           ,&
       & platform_name       ,&
       & errorstatus_success ,&
       & errorstatus_fatal   ,&
       & lensection  

  Use parkind1, Only : jpim, jprb, jplm
  Implicit None


  ! subroutine arguments
  ! scalar arguments with intent(in):
  Integer(Kind=jpim), Optional, Intent(in) :: kmyproc  ! logical processor id
  Integer(Kind=jpim), Optional, Intent(in) :: knproc   ! Number of processors
  Integer(Kind=jpim), Optional, Intent(in) :: kioproc  ! processor  dedicated for io
  Integer(Kind=jpim), Optional, Intent(in) :: file_id  ! file logical unit number

  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent (out) :: errorstatus      ! return code

  ! array arguments with intent(in):
  Type( rttov_coef ), Intent (in) :: coef_rttov        ! clear-sky coefficients

  ! array arguments with intent(inout):
  Type( rttov_scatt_coef ), Intent (inout) :: coef_scatt ! coefficients

!INTF_END

#include "rttov_errorreport.h"
#include "rttov_findnextsection.h"
#include "rttov_skipcommentline.h"
#include "rttov_opencoeff.h"

! local variables
  Integer(Kind=jpim)    :: file_lu, inst, platform, i, j, k
  Integer(Kind=jpim)    :: io_status
  Logical(Kind=jplm)    :: existence
  Logical(Kind=jplm)    :: file_toclose
  Logical(Kind=jplm)    :: file_open
  Character (len=32)    :: NameOfRoutine = 'rttov_readscattcoeffs' ! name for error message
  Character (len=132)   :: ErrMessage    ! error message string
  Character (len=256)   :: coeffname     ! file name for coefficient file
  Character (len=lensection) :: section
  Integer(Kind=jpim)    :: ITAG,iioproc,inproc,imyproc
REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',0_jpim,ZHOOK_HANDLE)

  If ( .Not. Present (kmyproc) ) Then
     imyproc = 1
  Else
     imyproc = kmyproc
  Endif

  If ( .Not. Present (kioproc) ) Then
     iioproc = 1
  Else
     iioproc = kioproc
  Endif

  If ( .Not. Present (knproc) ) Then
     inproc = 1
  Else
     inproc = knproc
  Endif

#ifndef _RTTOV_DO_DISTRIBCOEF
  If (inproc /= 1) Then
    errorstatus = errorstatus_fatal
    Write( errMessage, '( "Multiple processors need routine compiled with _RTTOV_DO_DISTRIBCOEF" )' ) 
    Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
    IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
    Return  
  Endif
#endif

  ITAG=1

  errorstatus     = errorstatus_success

  If ( Present (file_id) ) Then
     ! Scatt coefficient file has been opened externally
     file_lu = file_id
     file_toclose = .FALSE.

     Inquire( file_lu, OPENED = file_open )
     If ( .NOT. file_open ) Then
        errorstatus = errorstatus_fatal
        Write( errMessage, '( "File is not open on unit: ",i5 )' ) file_lu
        Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
        Return
     End If
  Else
     ! Open the scatt coefficients internally
     file_lu = 0
     file_toclose = .TRUE.

     platform = coef_rttov % id_platform
     inst     = coef_rttov % id_inst
     coeffname = 'mietable_'//Trim(platform_name(platform))//'_'//Trim(inst_name(inst))//'.dat'

     Inquire( FILE = coeffname, EXIST = existence )
     If ( .Not. existence ) Then
          errorstatus = errorstatus_fatal
          Write( errMessage, '( "Coefficient file, ", a, " not found." )' ) &
               & Trim( coeffname ) 
          Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
          IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
          Return
     End If

     if(imyproc == iioproc .or. inproc == 1)then
       Call rttov_opencoeff (errorstatus, coeffname, file_lu)
     endif
     if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
        Call broadcint(errorstatus, 1,iioproc,ITAG)
#endif        
     endif

     If (errorstatus /= errorstatus_success) Then
       ! rttov_opencoeff will have already reported an error
       errorstatus = errorstatus_fatal
       write(0,*) 'RTTOV_READSCATTCOEFFS: Error Opening File ' // coeffname
       IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
       Return
     Endif
  Endif

  readfile: Do
   if(imyproc == iioproc .or. inproc == 1)then
     Call rttov_findnextsection( file_lu, io_status, section )
   endif
   if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
      Call broadcint(io_status, 1,iioproc,ITAG)
      Call broadcchar(section, lensection, iioproc,ITAG)
#endif
   endif
     If ( io_status < 0 ) Exit !end-of-file

     ! error message if any problem when reading
     errMessage = 'io status while reading section '//section
   if(imyproc == iioproc .or. inproc == 1)then
     Call rttov_skipcommentline ( file_lu, io_status )
   endif
   if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
      Call broadcint(io_status, 1,iioproc,ITAG)
#endif
   endif
     If(io_status /= 0) Then
           Call rttov_errorreport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
           Return
     Endif

     Select Case( Trim(section) )


     Case( 'IDENTIFICATION' )
        if(imyproc == iioproc .or. inproc == 1)then
          Read(file_lu,*)  ! platform instrument in id
          Read(file_lu,*)  ! platform instrument in letters
          Read(file_lu,*)  ! sensor type [ir,mw,hi]
          Read(file_lu,*)  ! RTTOV compatibility version
          Read(file_lu,*)  ! version
          Read(file_lu,*)  ! creation date
        endif

     Case( 'DIMENSIONS')

        if(imyproc == iioproc .or. inproc == 1)then
          Read(file_lu,*)  coef_scatt%mfreqm,  coef_scatt%mtype,  coef_scatt%mtemp,  coef_scatt%mwc 
        endif
        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcint(coef_scatt % mfreqm, 1,iioproc,ITAG)
          Call broadcint(coef_scatt % mtype , 1,iioproc,ITAG)
          Call broadcint(coef_scatt % mtemp , 1,iioproc,ITAG)
          Call broadcint(coef_scatt % mwc   , 1,iioproc,ITAG)
#endif          
        endif

        If (coef_scatt%mtype /= 5) Then
          errorstatus = errorstatus_fatal
          errMessage = 'Wrong no of hydrometeors in parameter file (should be 5)'
              ! liquid prec., solid prec., ice water, liquid water, total ice
          Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
          IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
          Return
        Endif
        coef_scatt % nhydro = coef_scatt%mtype

     Case( 'FREQUENCIES')

        Allocate (coef_scatt % mie_freq(coef_scatt%mfreqm))
        if(imyproc == iioproc .or. inproc == 1)then
          Read(file_lu,*)  coef_scatt%mie_freq (:)
        endif
        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcreal(coef_scatt % mie_freq, &
           &coef_scatt%mfreqm,&
           &iioproc,ITAG)
#endif      
        endif

     Case( 'HYDROMETEOR')

        if(imyproc == iioproc .or. inproc == 1)then
          Read(file_lu,*)  
        endif

     Case( 'CONVERSIONS')

        if(imyproc == iioproc .or. inproc == 1)then
          Read(file_lu,*) coef_scatt%conv_rain(:)
          Read(file_lu,*) coef_scatt%conv_sp  (:)
          Read(file_lu,*) coef_scatt%conv_liq(:)
          Read(file_lu,*) coef_scatt%conv_ice(:)
          Read(file_lu,*) coef_scatt%conv_totalice(:)
          Read(file_lu,*) 
          Read(file_lu,*) coef_scatt%offset_temp_rain
          Read(file_lu,*) coef_scatt%offset_temp_sp
          Read(file_lu,*) coef_scatt%offset_temp_liq
          Read(file_lu,*) coef_scatt%offset_temp_ice
          Read(file_lu,*) coef_scatt%offset_temp_totalice
          Read(file_lu,*)
          Read(file_lu,*) coef_scatt%scale_water, coef_scatt%offset_water
        endif

        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcreal(coef_scatt % conv_rain, 2,iioproc,ITAG)
          Call broadcreal(coef_scatt % conv_sp, 2,iioproc,ITAG)
          Call broadcreal(coef_scatt % conv_liq, 2,iioproc,ITAG)
          Call broadcreal(coef_scatt % conv_ice, 2,iioproc,ITAG)
          Call broadcreal(coef_scatt % conv_totalice, 2,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_temp_rain, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_temp_sp, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_temp_liq, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_temp_ice, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_temp_totalice, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % scale_water, 1,iioproc,ITAG)
          Call broadcreal(coef_scatt % offset_water, 1,iioproc,ITAG)
#endif
        endif

        coef_scatt%conv_rain(:) = 1._JPRB/coef_scatt%conv_rain(:)
        coef_scatt%conv_sp  (:) = 1._JPRB/coef_scatt%conv_sp  (:)
        coef_scatt%scale_water = 1._JPRB/coef_scatt%scale_water
        coef_scatt%offset_water = - coef_scatt%offset_water
        coef_scatt%from_scale_water = 10**( 1._JPRB / coef_scatt%scale_water )
        
     Case( 'EXTINCTION')

        Allocate (coef_scatt % ext(coef_scatt%mfreqm, coef_scatt%mtype, coef_scatt%mtemp, coef_scatt%mwc))
        ! The loops should be inverted for better efficiency, but generation program currently not appropriate
        if(imyproc == iioproc .or. inproc == 1)then
          Do i = 1, coef_scatt%mfreqm
            Do j = 1, coef_scatt%mtype
              Do k = 1, coef_scatt%mtemp
                Read(file_lu,'(5(1x,e23.16))') coef_scatt % ext(i,j,k,:)
              Enddo
            Enddo
          Enddo
        endif
        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcreal(coef_scatt % ext, &
           &coef_scatt%mfreqm*coef_scatt%mtype*coef_scatt%mtemp*coef_scatt%mwc,&
           &iioproc,ITAG)
#endif
        endif

     Case( 'ALBEDO')

        Allocate (coef_scatt % ssa(coef_scatt%mfreqm, coef_scatt%mtype, coef_scatt%mtemp, coef_scatt%mwc))
        if(imyproc == iioproc .or. inproc == 1)then
          Do i = 1, coef_scatt%mfreqm
            Do j = 1, coef_scatt%mtype
              Do k = 1, coef_scatt%mtemp
                Read(file_lu,'(5(1x,e23.16))') coef_scatt % ssa(i,j,k,:)
              Enddo
            Enddo
          Enddo
        endif
        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcreal(coef_scatt % ssa, &
           &coef_scatt%mfreqm*coef_scatt%mtype*coef_scatt%mtemp*coef_scatt%mwc,&
           &iioproc,ITAG)
#endif
        endif

     Case( 'ASYMMETRY')

        Allocate (coef_scatt % asp(coef_scatt%mfreqm, coef_scatt%mtype, coef_scatt%mtemp, coef_scatt%mwc))
        if(imyproc == iioproc .or. inproc == 1)then
          Do i = 1, coef_scatt%mfreqm
            Do j = 1, coef_scatt%mtype
              Do k = 1, coef_scatt%mtemp
                Read(file_lu,'(5(1x,e23.16))') coef_scatt % asp(i,j,k,:)
              Enddo
            Enddo
          Enddo
        endif
        if(inproc /= 1)then
#ifdef _RTTOV_DO_DISTRIBCOEF
          Call broadcreal(coef_scatt % asp, &
           &coef_scatt%mfreqm*coef_scatt%mtype*coef_scatt%mtemp*coef_scatt%mwc,&
           &iioproc,ITAG)
#endif
        endif
  
     Case default
        Cycle readfile

     End Select

  End Do readfile

    
  if(imyproc == iioproc .or. inproc == 1)then
    If ( file_toclose ) Then
      Close ( unit = file_lu )
    Endif
  endif

  IF (LHOOK) CALL DR_HOOK('RTTOV_READSCATTCOEFFS',1_jpim,ZHOOK_HANDLE)
  End Subroutine rttov_readscattcoeffs
