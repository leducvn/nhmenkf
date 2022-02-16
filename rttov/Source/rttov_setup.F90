!
Subroutine rttov_setup (&
       & ERR,      &! out
       & Err_unit,         &! in
       & verbosity_level,  &! in
       & opts,             &
       & coefs,            &
       & instrument,       &! in
       & channels  ,       &! in optional
       & channels_rec)      ! in optional
  !
  ! Description:
  !
  ! Setup routine for RTTOV
  ! Handling of error messages. (rttov_errorhandling)
  ! Read coefficients (rttov_readcoeffs)
  !
  ! Error messages will be sent on the optional unit number errunit.
  !      Default is the value defined in the module for constants.
  ! 
  ! The levels of verbosity are 
  !  0 = no error messages output 
  !  1 = FATAL errors only printed. these are errors which 
  !      mean that profile should be aborted (e.g. unphysical 
  !      profile input) 
  !  2 = WARNING errors only printed. Errors which can allow 
  !      the computation to continue but the results may be 
  !      suspect (e.g. profile outside basis profile limits) 
  !  3 = INFORMATION messages which inform the user about 
  !      the computation 
  !
  ! For each instrument:
  ! Read an ASCII or binary coefficient file and allocate coeff structure
  !   arrays according to the optional list of channels.
  ! The user can provide an optional list of channels in "channels" argument
  !  array to reduce the output coefficient structure to this list. This
  ! can be important for reducing the memory allocation required when running
  ! with advanced IR sounders (e.g. AIRS or IASI). If the user
  !  wants all channels the "channels" argument shall not be present.
  !
  !
  ! Copyright:
  !
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
  !
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version    Date       Comment
  !  1.0    10/03/2003   Original code (P Brunel)
  !  1.1    22/03/2007   Added check on read status (R Saunders)
  !  1.2    02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !
  ! Code Description:
  !   FORTRAN 90, following AAPP standards
  !
  ! Declarations
  !
  ! Global variables:
  ! Modules used:
  !
#include "throw.h"

  ! Imported Type Definitions:
  Use rttov_types, Only : &
        & rttov_coefs,           &
        & rttov_options 

  Use parkind1, Only : jpim
!INTF_OFF
  Use parkind1, Only : jprb
  Use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON
  Implicit None

  !
  ! Subroutine arguments
  !   Scalar arguments with intent(in):
  Type(rttov_options), Intent(in) :: opts
  Integer(Kind=jpim), Intent (in) :: Err_Unit        ! Logical error unit (<0 for default) 
  Integer(Kind=jpim), Intent (in) :: verbosity_level ! (<0 for default)
  Integer(Kind=jpim), Intent (in) :: instrument(3) ! Instrument triplet
         ! first dimension  : (platform, satellite identification, instrument) number
         ! second dimension : nsat
  Integer(Kind=jpim), Optional, Intent (in) :: channels    (:)   ! list of channels to extract (channels,msat)
  Integer(Kind=jpim), Optional, Intent (in) :: channels_rec(:)

  ! scalar arguments with intent(out):
  Integer(Kind=jpim),  Intent (out) :: ERR
  Type( rttov_coefs ), Intent (out) :: coefs

!INTF_END

#include "rttov_errorhandling.h"
#include "rttov_errorreport.h"
#include "rttov_read_coefs.h"
#include "rttov_init_coefs.h"

  ! Local scalars/arrays
  Integer(Kind=jpim) :: dimchans ! size of array channels for channels dimension
  Integer(Kind=jpim) :: dimchans_in ! size of array channels for channels dimension
  Integer(Kind=jpim) :: nchans   ! number of requested channels per instrument (0 = all)
  Integer(Kind=jpim) :: nchans_in   ! number of requested channels per instrument (0 = all)
  Integer(Kind=jpim) :: i        ! loop index
  Integer(Kind=jpim), allocatable :: channels_list(:) ! list of requested channels
  Integer(Kind=jpim), allocatable :: channels_rec_list(:) ! list of requested channels
  Real(Kind=jprb)    :: ZHOOK_HANDLE
  !- End of header --------------------------------------------------------
TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETUP', 0_jpim, ZHOOK_HANDLE)
  
  ! Error Handling setup routine
  call rttov_errorhandling(Err_Unit, verbosity_level)

  ! Check optional argument channels
  If( Present ( channels ) ) Then
     dimchans = Size( channels, dim=1 )
  Else
     dimchans = 0
  End If

  If( Present ( channels_rec ) ) Then
     dimchans_in = Size( channels_rec, dim=1 )
  Else
     dimchans_in = 0
  End If


     ! Finds the last non null channel for ninst
     nchans = 0
     if( Present ( channels ) ) Then
        do i = 1, dimchans
           if( channels(i) > 0 ) then
              nchans  = nchans + 1
           endif
        End Do
     Endif

     nchans_in = 0
     if( Present ( channels_rec ) ) Then
        do i = 1, dimchans_in
           if( channels_rec(i) > 0 ) then
              nchans_in  = nchans_in + 1
           endif
        End Do
     Endif


     If( nchans > 0 ) Then
        ! Some channels wanted, create a list of the
        ! selected channels without O values

        ! Allocate intermediate channels list
        Allocate ( channels_list ( nchans ), stat= ERR)
        THROWM( ERR .NE. 0, "allocation of channels_list")
        Allocate ( channels_rec_list ( nchans_in ), stat= ERR)
        THROWM( ERR .NE. 0, "allocation of channels_rec_list")

        ! Allocate intermediate channels list
        ! Create intermediate channels list (use nchans var. again)
        nchans = 0
        do i = 1, dimchans
           if( channels(i) > 0 ) then
              nchans  = nchans + 1
              channels_list(nchans) = channels(i)
           endif
        End Do


        nchans_in=0
        do i = 1, dimchans_in
           if( channels_rec(i) > 0 ) then
              nchans_in  = nchans_in + 1
              channels_rec_list(nchans_in) = channels_rec(i)
           endif
        End Do

        ! default value 

        Call rttov_read_coefs( err, coefs, opts,                   &
                           channels     = channels_list,           &
                           channels_rec = channels_rec_list,       &
                           instrument   = instrument)
        THROWM( ERR .NE. 0, "readcoeffs list channels with addpc")

        Deallocate ( channels_list , stat=ERR)
        THROWM( ERR .NE. 0, "deallocation of channels_list")

        Deallocate ( channels_rec_list , stat=ERR)
        THROWM( ERR .NE. 0, "deallocation of channels_rec_list")


     Else 
        Call rttov_read_coefs( err, coefs, opts, instrument  = instrument)
        THROWM( ERR .NE. 0, "readcoeffs all channels")
     Endif

     Call rttov_init_coefs  ( err, opts, coefs )
     THROWM( ERR .NE. 0, "initcoeffs error")

IF (LHOOK) CALL DR_HOOK('RTTOV_SETUP', 1_jpim, ZHOOK_HANDLE)
CATCH
IF (LHOOK) CALL DR_HOOK('RTTOV_SETUP', 1_jpim, ZHOOK_HANDLE)
End Subroutine rttov_setup
