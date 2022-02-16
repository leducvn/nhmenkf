  !
  ! Description:
  !    Parallel ( openmp ) version of RTTOV_DIRECT, _TL, _AD and _K.   
  !    This file contains the definition of RTTOV_PARALLEL_DIRECT, _TL, _AD and _K. 
  !    RTTOV_PARALLEL_DIRECT is activated when _RTTOV_PARALLEL_DIRECT is defined, etc...
  !    The interfaces of these parallel routines are identical to the original ones,
  !    plus the optional parameter nthreads.
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
  !    Copyright 2007, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !          09/2007    Creation P. Marguinaud & P. Brunel
  !          10/2008    Fix for po sounders P. Marguinaud
  !          11/2009    Add strategy and debug parameters P. Marguinaud
  !


Subroutine &
rttov_parallel_ad      &
      & ( errorstatus      &
      & , chanprof         &
      & , opts             &
      & , profiles         &
      & , profiles_ad      &
      & , coefs            &
      & , calcemis         &
      & , emissivity       &
      & , emissivity_ad    &
      & , emissivity_out   &
      & , emissivity_out_ad&
      & , transmission     &
      & , transmission_ad  &
      & , radiancedata     &
      & , radiancedata_ad  &
      & , traj             &
      & , traj_ad          &
      &,  pccomp           &
      &,  pccomp_ad        &
      &,  channels_rec     &
      &,  nthreads         &    
      &,  strategy         &    
      &,  debug            &    
& )   


  Use rttov_const, Only :  &
       & errorstatus_success,&
       & errorstatus_warning,&
       & errorstatus_fatal,  & 
       & sensor_id_po
  
! Imported Type Definitions:
  Use rttov_types, Only :            &
        & rttov_options,             &
        & rttov_chanprof,            &
        & rttov_pccomp,              &
        & rttov_coefs,               &
        & profile_Type,              &
        & transmission_type,         &
        & radiance_Type,             &
        & rttov_traj
                      
  
  Use parkind1, Only : jpim,jprb,jplm
  Implicit None



!subroutine arguments:
  Integer(Kind=jpim),            Intent(out)   :: errorstatus ! return flag
  Type(rttov_chanprof),          Intent(in)    :: chanprof(:)
  Type(rttov_options),           Intent(in)    :: opts
  Type(profile_Type),            Intent(in)    :: profiles( : ) ! Atmospheric profiles
  Type(rttov_coefs),             Intent(in)    :: coefs
  Logical(Kind=jplm),            Intent(in)    :: calcemis( : )   ! switch for emmissivity calc.
  Real(Kind=jprb),               Intent(inout) :: emissivity( : ) ! surface emmissivity
  Real(Kind=jprb),               Intent(inout) :: emissivity_out( : ) ! surface emmissivity
  Type(transmission_type),       Intent(inout) :: transmission   ! transmittances and layer optical depths
  Type(radiance_Type),           Intent(inout) :: radiancedata   ! radiances (mw/cm-1/ster/sq.m) and degK


  Type(profile_Type),            Intent(inout) :: profiles_ad( : ) 
  Real(Kind=jprb),               Intent(inout) :: emissivity_ad( : ) 
  Real(Kind=jprb),               Intent(inout) :: emissivity_out_ad( : ) 
  Type(transmission_type),       Intent(inout) :: transmission_ad
  Type(radiance_Type),           Intent(inout) :: radiancedata_ad 

  
  Type(rttov_traj),    Optional, Intent(inout) :: traj
  Type(rttov_traj),    Optional, Intent(inout) :: traj_ad
  Type(rttov_pccomp),  Optional, Intent(inout) :: pccomp
  Type(rttov_pccomp),  Optional, Intent(inout) :: pccomp_ad
  Integer(Kind=jpim),  Optional, Intent(in)    :: channels_rec( : )
  

  Integer(Kind=jpim),  Optional, Intent(in)    :: nthreads
  Integer(Kind=jpim),  Optional, Intent(in)    :: strategy ! 0 = no strategy (RTTOV computes band limits), 1 = profile-wise, 2 = channel-wise
  Logical(Kind=jplm),  Optional, Intent(in)    :: debug

!INTF_END

#include "rttov_ad.h"

#include "rttov_alloc_prof.h"
#include "rttov_alloc_rad.h"
#include "rttov_errorreport.h"
#include "rttov_init_prof.h"
#include "rttov_add_prof.h"




  Character(len=*), Parameter :: NameOfRoutine = &
    "rttov_parallel_ad"



  Integer(Kind=jpim) :: nprofiles12max
  Integer(Kind=jpim) :: k
  Integer(Kind=jpim) :: iband, nbands, ibandk
  Integer(Kind=jpim) :: ithread, mthreads
  Integer(Kind=jpim) :: nsplits
  Integer(Kind=jpim) :: ichannelk, iprofilek
  Integer(Kind=jpim) :: nlevels
  Integer(Kind=jpim) :: nchannels_po
  Integer(Kind=jpim) :: nchannels
  Integer(Kind=jpim) :: nprofiles
  Integer(Kind=jpim) :: strategy1
  Logical(Kind=jplm) :: debug1
  Integer            :: mstat
  Integer            :: npcscores, nchan_rec
  Integer(Kind=jpim),      Allocatable :: errorstatus_x( : )               ! errorstatus/bands
  Integer(Kind=jpim),      Allocatable :: ichannel1( : ), ichannel2( : )   ! channel limits
  Integer(Kind=jpim),      Allocatable :: iprofile1( : ), iprofile2( : )   ! profile limits
  Integer(Kind=jpim),      Allocatable :: ipc1( : ), ipc2( : )             ! PCscores limits
  Integer(Kind=jpim),      Allocatable :: ichannelrec1( : ), ichannelrec2( : ) ! PC reconstructed channel limits
  Type(rttov_chanprof),    Allocatable :: chanprof_x( :, : )               ! chanprof/bands

  Type(transmission_type), Allocatable :: transmission_x( : )         ! transmission/bands ( pointer assoc )
  Type(radiance_Type),     Allocatable :: radiancedata_x( : )         ! radiancedata/bands ( pointer assoc )
  Type(rttov_pccomp),      Allocatable :: pccomp_x( : )               ! pccomp/bands ( pointer assoc )


  Integer(Kind=jpim)                   :: errorstatus_ad_x
  Type(transmission_type), Allocatable :: transmission_ad_x( : )  
  Type(radiance_Type),     Allocatable :: radiancedata_ad_x( : )
  Type(rttov_pccomp),      Allocatable :: pccomp_ad_x( : )
  Type(profile_Type),      Allocatable :: profiles_ad_x( :, : ) 
  


! Activated if openmp
!  
!$  Integer, External :: omp_get_thread_num
!$  Integer, External :: omp_get_max_threads
!

  nchannels = size(chanprof)
  nprofiles = size(profiles)

  debug1 = .FALSE.
  If( Present( debug ) ) debug1 = debug

  If( Present( traj ) ) Then
    Call rttov_errorreport( errorstatus_warning, "traj argument will not be used !!", NameOfRoutine )
  EndIf
  If( Present( traj_ad ) ) Then
    Call rttov_errorreport( errorstatus_warning, "traj_ad argument will not be used !!", NameOfRoutine )
  EndIf

  nlevels = profiles( 1 ) % nlevels

  strategy1 = 0
  If( Present( strategy ) ) Then
    strategy1 = strategy
  EndIf
  ! NB Must distribute work profile-wise among threads for PC-RTTOV
  If ( opts % addpc ) strategy1 = 1

  If( Present( nthreads ) ) Then
    mthreads = nthreads
  Else
! Default value for mthreads
    mthreads = 1_jpim
! Activated if openmp
!  
!$  mthreads = omp_get_max_threads()
!
  EndIf

! Check we do not have more bands than channels

  nsplits = 1_jpim
  nbands = mthreads * nsplits


  If( strategy1 .eq. 0 ) Then
    If( coefs%coef%id_sensor .eq. sensor_id_po ) Then
      If( nbands .gt. nprofiles ) nbands = nprofiles
    Else
      If( nbands .gt. nchannels ) nbands = nchannels
    EndIf
  Else If( strategy1 .eq. 1 ) Then
      nbands = nprofiles
  Else If( strategy1 .eq. 2 ) Then
      If( coefs%coef%id_sensor .eq. sensor_id_po ) Then
        errorstatus = errorstatus_fatal
        Call rttov_errorreport( errorstatus_fatal, "Polarized sounders cannot be run channel by channel", NameOfRoutine )
        Return
      EndIf
      nbands = nchannels
  Else
    errorstatus = errorstatus_fatal
    Call rttov_errorreport( errorstatus_fatal, "Unknown strategy", NameOfRoutine )
    Return
  Endif
  

  If( debug1 ) Then
  Write( *, '(A20," = ",I8)' ) "nprofiles", nprofiles
  Write( *, '(A20," = ",I8)' ) "nchannels", nchannels
  Write( *, '(A20," = ",I5)' ) "nbands", nbands
  Endif
  
  Allocate( ichannel1( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  Allocate( ichannel2( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  
  If( opts%addpc ) Then
    Allocate( ipc1( nbands ), stat = mstat )
    If( mstat .ne. 0 ) Goto 100
    Allocate( ipc2( nbands ), stat = mstat )
    If( mstat .ne. 0 ) Goto 100
    If( opts%addradrec ) Then
      Allocate( ichannelrec1( nbands ), stat = mstat )
      If( mstat .ne. 0 ) Goto 100
      Allocate( ichannelrec2( nbands ), stat = mstat )
      If( mstat .ne. 0 ) Goto 100
    EndIf
  EndIf
  
  Allocate( iprofile1( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  Allocate( iprofile2( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  
  Allocate( transmission_x( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  Allocate( radiancedata_x( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  If( opts%addpc ) Then
    Allocate( pccomp_x( nbands ), stat = mstat )
    If( mstat .ne. 0 ) Goto 100
  EndIf
  

  Allocate( transmission_ad_x( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  Allocate( radiancedata_ad_x( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  If( opts%addpc ) Then
    Allocate( pccomp_ad_x( nbands ), stat = mstat )
    If( mstat .ne. 0 ) Goto 100
  EndIf
  
  
  If ( coefs%coef%id_sensor .eq. sensor_id_po ) Then
    ! The parallel work for polarised sensors is divided by profile. In the extreme 
    ! case, one profile may be run for all channels, and another may be run for just 
    ! one channel. Therefore we ensure chanprof_x is large enough to handle every 
    ! possible contingency (no matter how unlikely)
    Allocate( chanprof_x( (nprofiles / nbands + 1_jpim) * coefs%coef%fmv_chn, nbands ), stat = mstat )
  Else
    Allocate( chanprof_x( nchannels / nbands + nbands, nbands ), stat = mstat )
  EndIf
  If( mstat .ne. 0 ) Goto 100
  Allocate( errorstatus_x( nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  
  
  ichannel1( : ) = 0_jpim
  ichannel2( : ) = 0_jpim
  
  If( opts%addpc ) Then
    ipc1( : ) = 0_jpim
    ipc2( : ) = 0_jpim
    If( opts%addradrec ) Then
      ichannelrec1( : ) = 0_jpim
      ichannelrec2( : ) = 0_jpim
    EndIf
  EndIf
  
  iprofile1( : ) = 0_jpim
  iprofile2( : ) = 0_jpim
  
  !
  ! Build the array of lower/upper bounds of bands
  !

  If( ( strategy1 .eq. 0 ) .and. ( coefs%coef%id_sensor .eq. sensor_id_po ) ) Then
    !
    ! The approach is different for PO sensors
    ! because some channels cannot be separated in the
    ! computations; hence we must split the workload
    ! profile-wise
    !
    ichannelk = 1_jpim
    iprofilek = 1_jpim
    Do iband = 1, nbands
      nchannels_po = 0_jpim
      ichannel1( iband ) = ichannelk
      Do While( ( iprofilek .le. nprofiles ) &
          .and. ( ichannelk .lt. nchannels ) )
        Do While( ( iprofilek .le. nprofiles ) &
            .and. ( ichannelk .lt. nchannels ) )
          ichannelk = ichannelk + 1_jpim
          nchannels_po = nchannels_po + 1_jpim
          If( chanprof( ichannelk )%prof .ne. iprofilek ) Then
            iprofilek = iprofilek + 1_jpim
            Exit
          End If
        End Do
        ! This logic splits the profiles more evenly among the bands
        If (iband .eq. 1) Then
          If( nchannels_po .ge. nchannels / nbands ) Exit
        Else
          If( nchannels_po .ge. (nchannels - ichannel2( iband-1 )) / (nbands - (iband-1))) Exit
        EndIf
      End Do
      ichannel2( iband ) = ichannelk - 1_jpim
    End Do
    ! What remains is given to the last band
    ichannel2( nbands ) = nchannels
  Else If( strategy1 .eq. 1 ) Then
    If( opts%addpc ) Then
      npcscores = size(pccomp%pcscores) / size(profiles)
      If( opts%addradrec ) nchan_rec = size(channels_rec)
    EndIf
    ichannelk = 1_jpim
    Do iband = 1, nbands
      ichannel1( iband ) = ichannelk
      iprofilek = chanprof(ichannelk)%prof
      Do
        If( ichannelk .gt. nchannels ) Exit
        If( iprofilek .gt. nprofiles ) Exit
        If( iprofilek .ne. chanprof(ichannelk)%prof ) Exit
        ichannelk = ichannelk + 1_jpim
      EndDo
      ichannel2( iband ) = ichannelk-1_jpim
      If( opts%addpc ) Then
        ! With one profile per band, splitting the pcscores and reconstructed channels
        ! among bands is quite simple
        ipc1( iband ) = (iband-1_jpim) * npcscores + 1_jpim
        ipc2( iband ) = ipc1( iband ) + npcscores - 1_jpim        
        If( opts%addradrec ) Then
          ichannelrec1( iband ) = (iband-1_jpim) * nchan_rec + 1_jpim
          ichannelrec2( iband ) = ichannelrec1( iband ) + nchan_rec - 1_jpim
        EndIf
      EndIf
    EndDo
  Else
    !
    ! For IR & MW sensors, we try to make bands
    ! with the same number of channels
    !
    ichannelk = 1_jpim
    Do iband = 1, nbands
      ichannel1( iband ) = ichannelk
      ichannelk = ichannelk + nchannels / nbands - 1_jpim
      If( iband .eq. nbands ) ichannelk = nchannels
      ichannel2( iband ) = ichannelk
      ichannelk = ichannelk + 1_jpim
    End Do
  End If


 
  chanprof_x( :, : )%chan = 0_jpim
  chanprof_x( :, : )%prof = 0_jpim
  
  Do iband = 1, nbands
    iprofile1( iband ) = chanprof( ichannel1( iband ) )%prof
    iprofile2( iband ) = chanprof( ichannel2( iband ) )%prof
    chanprof_x( 1 : ichannel2( iband ) - ichannel1( iband ) + 1, iband )%chan = &
          & chanprof( ichannel1( iband ) : ichannel2( iband ) )%chan
    Do ichannelk = ichannel1( iband ), ichannel2( iband )
      chanprof_x( ichannelk - ichannel1( iband ) + 1, iband )%prof &
        & = chanprof( ichannelk )%prof - iprofile1( iband ) + 1
    EndDo
  EndDo 
  
  !
  ! for rttov_ad, we have to allocate a temporary array holding
  ! results from every band, and add the contribution of every band
  ! after the parallel loop
  !
  nprofiles12max = maxval( iprofile2( 1 : nbands ) - iprofile1( 1 : nbands ) + 1 )
  Allocate( profiles_ad_x( nprofiles12max, nbands ), stat = mstat )
  If( mstat .ne. 0 ) Goto 100
  Do iband = 1, nbands
    call rttov_alloc_prof( errorstatus_ad_x, nprofiles12max, profiles_ad_x( :, iband ), nlevels, opts, 1_jpim, & 
      init = .true._jplm, coefs=coefs )
    If( errorstatus_ad_x .ne. 0 ) Goto 100
    Call rttov_init_prof(profiles_ad_x( :, iband ))
  EndDo
  
  
  If( debug1 ) Then

  Do iband = 1, nbands
    Write( *, * )
    Write( *, '(A20," = ",I5)' ) "band", iband
    Write( *, '(A20," = ",100(I5))' ) "channels", chanprof( ichannel1( iband ) : ichannel2( iband ) )%chan
    Write( *, '(A20," = ",100(I5))' ) "lprofiles", chanprof( ichannel1( iband ) : ichannel2( iband ) )%prof
    Write( *, '(A20," = ",100(I5))' ) "chanprof_x%chan", chanprof_x( :, iband )%chan
    Write( *, '(A20," = ",100(I5))' ) "chanprof_x%prof", chanprof_x( :, iband )%prof
    Write( *, * )
  EndDo

  Endif
  
  
  !
  ! we make the correct pointers associations between temporary transmissions
  ! and radiances and arrays passed as input
  !
  Do iband = 1, nbands
  
    Call AssociateTransmission( ichannel1( iband ), ichannel2( iband ), &
     & transmission_x( iband ), transmission )
    Call AssociateRadiance( ichannel1( iband ), ichannel2( iband ), &
     & radiancedata_x( iband ), radiancedata )

    If( opts%addpc ) Then
      If( opts%addradrec ) Then
        Call AssociatePCcomp( ipc1( iband ), ipc2( iband ), &
        & pccomp_x( iband ), pccomp, ichannelrec1( iband ), ichannelrec2( iband ) )
      Else
        Call AssociatePCcomp( ichannel1( iband ), ichannel2( iband ), &
        & pccomp_x( iband ), pccomp )
      EndIf
    EndIf 

    
    Call AssociateTransmission( ichannel1( iband ), ichannel2( iband ), &
     & transmission_ad_x( iband ), transmission_ad )
    Call AssociateRadiance( ichannel1( iband ), ichannel2( iband ), &
     & radiancedata_ad_x( iband ), radiancedata_ad )
     
    If( opts%addpc ) Then
      If( opts%addradrec ) Then
        Call AssociatePCcomp( ipc1( iband ), ipc2( iband ), &
        & pccomp_ad_x( iband ), pccomp_ad, ichannelrec1( iband ), ichannelrec2( iband ) )
      Else
        Call AssociatePCcomp( ichannel1( iband ), ichannel2( iband ), &
        & pccomp_ad_x( iband ), pccomp_ad )
      EndIf
    EndIf 
    
    
    
  End Do

  
  
  
!$OMP PARALLEL DO PRIVATE( iband ) NUM_THREADS( mthreads ) SCHEDULE( DYNAMIC )
  Do iband = 1, nbands

  If( debug1 ) Then
!    Print *, "================================================"
!    Print *, "THREAD = ", OMP_GET_THREAD_NUM(), " BAND = ", iband, &
!  & " ichannel1 = ", ichannel1( iband ), " ichannel2 = ", ichannel2( iband )

    Write( *, * )
    Write( *, '(A)' ) "=========================== rttov_ad ==========================="
    Write( *, '(A20)' ) "errorstatus"
    Write( *, '(A20," = ",I5)' ) "nprofiles", iprofile2( iband ) - iprofile1( iband ) + 1
    Write( *, '(A20," = ",I5)' ) "nchannels", ichannel2( iband ) - ichannel1( iband ) + 1
    Write( *, '(A20," = ",100(I5))' ) "channels", chanprof_x( :, iband )%chan
    Write( *, '(A20," = ",100(I5))' ) "lprofiles", chanprof_x( :, iband )%prof
    Write( *, '(A20," = ",L5)' ) "opts", opts
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "profiles",          &
    &"profiles", iprofile1( iband ), iprofile2( iband )
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "profiles_ad",       &
    &"profiles_ad", iprofile1( iband ), iprofile2( iband )
    Write( *, '(A20)' ) "coefs"
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "calcemis",          &
    &"calcemis", ichannel1( iband ), ichannel2( iband )
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "emissivity",        &
    &"emissivity", ichannel1( iband ), ichannel2( iband )
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "emissivity_ad",     &
    &"emissivity_ad", ichannel1( iband ), ichannel2( iband )
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "emissivity_out",    &
    &"emissivity_out", ichannel1( iband ), ichannel2( iband )
    Write( *, '(A20," = ",A20," ( ",I5," : ",I5," ) ")' ) "emissivity_out_ad", &
    &"emissivity_out_ad", ichannel1( iband ), ichannel2( iband )
    Write( *, '(A20)' ) "transmission"
    Write( *, '(A20)' ) "transmission_ad"
    Write( *, '(A20)' ) "radiancedata"
    Write( *, '(A20)' ) "radiancedata_ad"
    Write( *, * )
    
  Endif


  If( opts%addpc ) Then
      Call &
        & rttov_ad                                                               &
        & ( errorstatus_x( iband )                                               &
        & , chanprof_x( 1 : ichannel2( iband ) - ichannel1( iband ) + 1, iband ) &
        & , opts                                                                 &
        & , profiles( iprofile1( iband ) : iprofile2( iband ) )                  &
        & , profiles_ad_x( :, iband )                                            &
        & , coefs                                                                &
        & , calcemis( ichannel1( iband ) : ichannel2( iband ) )                  &
        & , emissivity( ichannel1( iband ) : ichannel2( iband ) )                &
        & , emissivity_ad( ichannel1( iband ) : ichannel2( iband ) )             &
        & , emissivity_out( ichannel1( iband ) : ichannel2( iband ) )            &
        & , emissivity_out_ad( ichannel1( iband ) : ichannel2( iband ) )         &
        & , transmission_x( iband )                                              &
        & , transmission_ad_x( iband )                                           &
        & , radiancedata_x( iband )                                              &
        & , radiancedata_ad_x( iband )                                           &
        &,  pccomp         = pccomp_x( iband )                                   &
        &,  pccomp_ad      = pccomp_ad_x( iband )                                &
        &,  channels_rec   = channels_rec                                        &
        & )
   Else
      Call &
        & rttov_ad                                                               &
        & ( errorstatus_x( iband )                                               &
        & , chanprof_x( 1 : ichannel2( iband ) - ichannel1( iband ) + 1, iband ) &
        & , opts                                                                 &
        & , profiles( iprofile1( iband ) : iprofile2( iband ) )                  &
        & , profiles_ad_x( :, iband )                                            &
        & , coefs                                                                &
        & , calcemis( ichannel1( iband ) : ichannel2( iband ) )                  &
        & , emissivity( ichannel1( iband ) : ichannel2( iband ) )                &
        & , emissivity_ad( ichannel1( iband ) : ichannel2( iband ) )             &
        & , emissivity_out( ichannel1( iband ) : ichannel2( iband ) )            &
        & , emissivity_out_ad( ichannel1( iband ) : ichannel2( iband ) )         &
        & , transmission_x( iband )                                              &
        & , transmission_ad_x( iband )                                           &
        & , radiancedata_x( iband )                                              &
        & , radiancedata_ad_x( iband )                                           &
        & )
   EndIf
!     Print *, iband, errorstatus_x(iband )
   
  End Do
!$OMP END PARALLEL DO


 !
 ! add every contribution in output profile_k array
 !
  Do iband = 1, nbands
    Do iprofilek = iprofile1( iband ), iprofile2( iband )
!      Print *, "profiles_ad( ", iprofilek, " ) = ", &
!      & "profiles_ad( ", iprofilek, " ) + profiles_ad_x( ", &
!      & iprofilek - iprofile1( iband ) + 1, ", ", iband, " )"
      k = iprofilek - iprofile1( iband ) + 1
      Call rttov_add_prof( profiles_ad( iprofilek:iprofilek ), &
                         & profiles_ad( iprofilek:iprofilek ), &
                         & profiles_ad_x( k:k, iband ) )
    EndDo
  EndDo




!  Do iband = 1, nbands
!    Call CopyTransmission( ichannel1( iband ), ichannel2( iband ), transmission_x( iband ), transmission )
!  EndDo

  !
  ! we have to construct the global errorstatus array 
  ! from what we got from every thread
  !
  errorstatus = errorstatus_success
  If (Any(errorstatus_x( : ) /= errorstatus_success)) errorstatus = errorstatus_fatal

!  Do iband = 1, nbands
!    Do iprofilek = iprofile1( iband ), iprofile2( iband )
!      Print *, iprofilek, " <= ", iband, iprofilek - iprofile1( iband ) + 1, &
!      & errorstatus_x( iprofilek - iprofile1( iband ) + 1, iband )
!    EndDo
!  EndDo

  
  DeAllocate( chanprof_x )
  DeAllocate( errorstatus_x )

  DeAllocate( transmission_x )
  DeAllocate( radiancedata_x )

  If( opts%addpc ) DeAllocate( pccomp_x )

  
  DeAllocate( transmission_ad_x )
  DeAllocate( radiancedata_ad_x )
  If( opts%addpc ) DeAllocate( pccomp_ad_x )



  Do iband = 1, nbands
    call rttov_alloc_prof( errorstatus_ad_x, nprofiles12max, profiles_ad_x( :, iband ), nlevels, opts, 0_jpim, coefs = coefs)
  EndDo
  DeAllocate( profiles_ad_x )
  
  DeAllocate( iprofile2 )
  DeAllocate( iprofile1 )

  If( opts%addpc ) Then
    DeAllocate( ipc2 )
    DeAllocate( ipc1 )
    If( opts%addradrec )Then
      DeAllocate( ichannelrec2 )
      DeAllocate( ichannelrec1 )
    EndIf
  EndIf

  DeAllocate( ichannel2 )
  DeAllocate( ichannel1 )

  
  Return
  
100 Continue
    
    Call rttov_errorreport( errorstatus_fatal, "Memory allocation failed", NameOfRoutine )

    Stop  
  
Contains

  ! ancillary routines

  Subroutine AssociateRadiance( ichannel1, ichannel2, radiancedata_x, radiancedata )
    Integer(Kind=jpim), Intent(in) :: ichannel1, ichannel2
    Type(radiance_Type), Intent(in) :: radiancedata
    Type(radiance_Type), Intent(out) :: radiancedata_x
    
    radiancedata_x % clear &
      => radiancedata % clear( ichannel1 : ichannel2 )

    radiancedata_x % cloudy &
      => radiancedata % cloudy( ichannel1 : ichannel2 )
    
    radiancedata_x % total &
      => radiancedata % total( ichannel1 : ichannel2 )
    
    radiancedata_x % bt &
      => radiancedata % bt( ichannel1 : ichannel2 )
    
    radiancedata_x % bt_clear &
      => radiancedata % bt_clear( ichannel1 : ichannel2 )
    
    radiancedata_x % upclear &
      => radiancedata % upclear( ichannel1 : ichannel2 )
    
    radiancedata_x % dnclear &
      => radiancedata % dnclear( ichannel1 : ichannel2 )
    
    radiancedata_x % reflclear &
      => radiancedata % reflclear( ichannel1 : ichannel2 )
    
    radiancedata_x % overcast &
      => radiancedata % overcast( :, ichannel1 : ichannel2 )
  
    radiancedata_x % up &
      => radiancedata % up( :, ichannel1 : ichannel2 )
  
    radiancedata_x % down &
      => radiancedata % down( :, ichannel1 : ichannel2 )
  
    radiancedata_x % surf &
      => radiancedata % surf( :, ichannel1 : ichannel2 )
  
  End Subroutine

  Subroutine AssociateTransmission( ichannel1, ichannel2, transmission_x, transmission )
    Integer(Kind=jpim), Intent(in) :: ichannel1, ichannel2
    Type(transmission_type), Intent(in) :: transmission
    Type(transmission_type), Intent(out) :: transmission_x
    
    transmission_x % tau_levels & 
      => transmission % tau_levels( :, ichannel1 : ichannel2 )

    transmission_x % tau_total &
      => transmission % tau_total( ichannel1 : ichannel2 )
      
  End Subroutine
  
  Subroutine AssociatePCcomp( ipc1, ipc2, pccomp_x, pccomp, ichannelrec1, ichannelrec2 )
    Integer(Kind=jpim), Intent(in) :: ipc1, ipc2
    Type(rttov_pccomp), Intent(in) :: pccomp
    Type(rttov_pccomp), Intent(out) :: pccomp_x
    Integer(Kind=jpim), Intent(in), Optional :: ichannelrec1, ichannelrec2
    
    pccomp_x % pcscores & 
      => pccomp % pcscores( ipc1 : ipc2 )

    If( Present(ichannelrec1) ) Then
      pccomp_x % bt_pccomp &
        => pccomp % bt_pccomp( ichannelrec1 : ichannelrec2 )
        
      pccomp_x % clear_pccomp &
        => pccomp % clear_pccomp( ichannelrec1 : ichannelrec2 )
    EndIf

  End Subroutine
  
End Subroutine






