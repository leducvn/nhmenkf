!
Subroutine rttov_user_profile_checkinput( &
       & opts,               &
       & prof,               &
       & coefs,              &
       & ERR                 )
  ! Description:
  ! Check input profile/angles
  ! (i)  Are physically realistic
  ! (ii) Profile values are within the basis set used to
  !      generate the coefficients
  ! Unphysical values return a fatal error status
  ! Profile values outside the basis set return a warning status
  !
  ! Valid for one profile on any pressure levels
  ! profile levels are tested against the nearest higher pressure coefficient level
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
  ! Method: Check input profiles with fixed limits specified
  !         in constants and coeff file.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       16/07/2010  P Brunel
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:

!INTF_OFF
#include "throw.h"
!INTF_ON

  ! Imported Type Definitions:
  Use rttov_types, Only :      &
         & rttov_coefs       , &
         & rttov_options     , &
         & profile_Type

  Use parkind1, Only : jpim

!INTF_OFF
Use rttov_const, Only : &
    nsurftype,&
    nwatertype,&
    gas_id_watervapour,&
    gas_id_ozone,&
    gas_id_co2,&
    gas_id_co,&
    gas_id_ch4,&
    gas_id_n2o,&
    tmax,&
    tmin,&
    qmax,&
    qmin,&
    o3max,&
    o3min,&
    co2max,&
    co2min,&
    comax,&
    comin,&
    n2omax,&
    n2omin,&
    ch4max,&
    ch4min,&
    clwmax,&
    clwmin,&
    pmax,&
    pmin,&
    wmax,&
    zenmax,&
    ctpmax,&
    ctpmin,&
    bemax,&
    bemin
Use yomhook, Only : &
    LHOOK,&
    DR_HOOK
Use parkind1, Only : &
    jprb, &
    jplm
!INTF_ON

  Implicit None

  ! subroutine arguments
  ! scalar arguments with intent(in):
  Type(rttov_options),Intent (in) :: opts    ! rttov options
  Type(profile_Type), Intent (in) :: prof    ! input profiles
  Type( rttov_coefs), Intent (in) :: coefs   ! coefficients

  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent (out) :: ERR    ! return code

!INTF_END

#include "rttov_errorreport.h"


  !local variables:


  Real(Kind=jprb)    :: wind
  Integer(Kind=jpim) :: firstlevel, ilev, jlev
  Integer(Kind=jpim) :: ilay
  Integer(Kind=jpim) :: ig
  Logical(Kind=jplm) :: cfracflag
  Logical(Kind=jplm) :: OK
  Character(len=256) :: msg
  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
TRY

  !-------------
  !0. Initialize
  !-------------

IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

  !determine first pressure level above the surface (note levels are top down)
  Do firstlevel = prof % nlevels, 1, -1
     If ( prof % p(firstlevel) <= prof % s2m % p ) Exit
  End Do
 

  !------------------------------
  !1. Check for unphysical values
  !------------------------------
  ! zenith angle
  If ( prof % zenangle > zenmax .Or. &
     & prof % zenangle < 0._JPRB          ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid zenith angle")
  End If
  
  ! Cloud Top Pressure
  If ( prof % ctp > ctpmax .Or. &
     & prof % ctp < ctpmin      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid cloud top pressure")
  End If

  ! Cloud Fraction
  If ( prof % cfraction > 1._JPRB .Or. &
     & prof % cfraction < 0._JPRB      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid cloud fraction")
  End If
  
  ! Zeeman variables
  If (coefs%coef % inczeeman) Then
    ! Magnetic field strength
    If ( prof % be > bemax .Or. &
       & prof % be < bemin      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid magnetic field strength")
    End If

    ! Cosine of angle between path and mag. field
    If ( prof % cosbk > 1._JPRB  .Or. &
       & prof % cosbk < -1._JPRB      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid cosbk")
    End If    
  End If

  !1.1 surface variables
  !---------------------

  ! Pressure
  If ( prof % s2m % p > pmax .Or. &
     & prof % s2m % p < pmin      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid surface pressure")
  End If

  ! 2m air temperature
  If ( prof % s2m % t    > tmax .Or. &
     & prof % s2m % t    < tmin ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid 2m air temperature")
  End If

  ! 2m water vapour - only used if opts % use_q2m is TRUE
  If ( opts % use_q2m ) Then
    If ( prof % s2m % q > qmax .Or. &
       & prof % s2m % q < qmin      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid 2m water vapour")
    End If
  End If

  !  surface wind speed
  wind = sqrt(&
          & prof % s2m % u * prof % s2m % u + &
          & prof % s2m % v * prof % s2m % v   )
  If ( wind > wmax .Or. &
     & wind < 0._JPRB      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid 10m wind speed")
  End If

  ! surface skin temperature
  If ( prof % skin % t > tmax .Or. &
     & prof % skin % t < tmin      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid skin surface temperature")
  End If

  ! surface type
  If ( prof % skin % surftype < 0 .Or. &
     & prof % skin % surftype > nsurftype ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid surface type")
  End If
  
  ! water type
  If ( prof % skin % watertype < 0 .Or. &
     & prof % skin % watertype > nwatertype ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid water type")
  End If

  
  !1.2 atmospheric variables
  !-------------------------

  ! Predictors are calculated on *all* levels so check
  ! the hard limits on every level to avoid errors.

  ! temperature
  If ( Any( prof % t(:) > tmax ) .Or. &
     & Any( prof % t(:) < tmin )      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid atmospheric temperature")
  End If

  ! water vapour
  If ( Any( prof % q(:) > qmax ) .Or. &
     & Any( prof % q(:) < qmin )      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid atmospheric water vapour")
  End If

  ! ozone
  If ( opts%ozone_Data .And. coefs%coef % nozone > 0) Then
     If ( Any( prof % o3(:) > o3max ) .Or. &
        & Any( prof % o3(:) < o3min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric ozone")
    End If
  End If

  ! CO2
  If ( opts%co2_Data .And. coefs%coef % nco2 > 0) Then
     If ( Any( prof % co2(:) > co2max ) .Or. &
        & Any( prof % co2(:) < co2min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CO2")
    End If
  End If

  ! CO
  If ( opts%co_Data .And. coefs%coef % nco > 0) Then
     If ( Any( prof % co(:) > comax ) .Or. &
        & Any( prof % co(:) < comin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CO")
    End If
  End If

  ! N2O
  If ( opts%n2o_Data .And. coefs%coef % nn2o > 0) Then
     If ( Any( prof % n2o(:) > n2omax ) .Or. &
        & Any( prof % n2o(:) < n2omin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric N2O")
    End If
  End If

  ! CH4
  If ( opts%ch4_Data .And. coefs%coef % nch4 > 0) Then
     If ( Any( prof % ch4(:) > ch4max ) .Or. &
        & Any( prof % ch4(:) < ch4min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CH4")
    End If
  End If

  ! cloud liquid water
  If ( opts%clw_Data ) Then
     If ( Any( prof % clw(:) > clwmax ) .Or. &
        & Any( prof % clw(:) < clwmin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid cloud liquid water")
    End If
  Endif

  ! Cloud input profile
  If ( opts%addclouds ) Then
     If ( Any( prof % cfrac(:,:) > 1._JPRB ) .Or. &
        & Any( prof % cfrac(:,:) < 0._JPRB )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid cloud profile fraction (cfrac)")
    End If
    
    If ( Any( prof % cloud(:,:) < 0._JPRB ) ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid cloud concentration")
    End If
    
    If ( Any( prof % icede(:) < 0._JPRB ) ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid ice effective diameter")
    End If    
  End If

  ! Check whether input cloud profile has more than one non-zero cfrac on each layer
  If ( opts%addclouds ) Then
    cfracflag = .FALSE.
    Do ilay = 1, prof % nlayers
      If ( SUM(prof % cfrac(:,ilay)) > MAXVAL(prof % cfrac(:,ilay)) ) cfracflag = .TRUE.
    End Do
    If ( cfracflag ) Then
      ERR = ERRORSTATUS_FATAL
      msg = "there must be at most one non-zero cfrac value specified per layer"
      THROWM( ERR .NE. 0 , msg)
    End If        
  End If

  ! Aerosol input profile
  If ( opts%addaerosl ) Then
    If ( Any( prof % aerosols(:,:) < 0._JPRB ) ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid aerosol concentration")
    End If
  End If


  !-----------------------------
  !2. Check against basis values
  !-----------------------------

  if(.not.opts%addpc)then

    jlev = 1
    Do ilev = 1, firstlevel
      ! recherche des niveaux RTTOV qui encadrent le niveau user, les pressions sont
      ! arondies 0.1Pa
      If (jlev < coefs%coef % nlevels) Then
        OK =  nint(prof % p(ilev)*1000) >= nint(coefs%coef % lim_prfl_p(jlev)*1000)   .and. &
           &  nint(prof % p(ilev)*1000) <  nint(coefs%coef % lim_prfl_p(jlev+1)*1000)
      End If
      
      Do while (  .not. OK .and. jlev < coefs%coef % nlevels)
        jlev = jlev+1
        If ( jlev < coefs%coef % nlevels ) Then
          OK =  nint(prof % p(ilev)*1000) >= nint(coefs%coef % lim_prfl_p(jlev)*1000)   .and. &
             &  nint(prof % p(ilev)*1000) <  nint(coefs%coef % lim_prfl_p(jlev+1)*1000)
        Else
          ! on est sur le dernier niveau des fichiers de coefs
          OK =  .true.
        End If
      End Do

      If  ( ( prof % t(ilev) > coefs%coef % lim_prfl_tmax(jlev) ) .or. &
         &  ( prof % t(ilev) < coefs%coef % lim_prfl_tmin(jlev) ) ) Then
        ERR = errorstatus_warning
        If ( opts%verbose_checkinput_warnings ) Then
          msg = "some atmospheric temperature outside coef. limits"
          WARN(msg)
        End If
      End If

      ig = coefs%coef % fmv_gas_pos( gas_id_watervapour )
      If  ( ( prof % q(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
         &  ( prof % q(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
        ERR = errorstatus_warning
        If ( opts%verbose_checkinput_warnings ) Then
          msg = "some atmospheric water vapour outside coef. limits"
          WARN(msg)
        End If
      End If

      If ( opts%ozone_Data .And. coefs%coef % nozone > 0) Then
        ig = coefs%coef % fmv_gas_pos( gas_id_ozone )
        If  ( ( prof % o3(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
           &  ( prof % o3(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric ozone outside coef. limits"
            WARN(msg)
          End If
        End If
      End If

      If ( opts%co2_Data .And. coefs%coef % nco2 > 0) Then
        ig = coefs%coef % fmv_gas_pos( gas_id_co2 )
        If  ( ( prof % co2(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
           &  ( prof % co2(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric co2 outside coef. limits"
            WARN(msg)
          End If
        End If
      End If

      If ( opts%co_Data .And. coefs%coef % nco > 0) Then
        ig = coefs%coef % fmv_gas_pos( gas_id_co )
        If  ( ( prof % co(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
           &  ( prof % co(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric co outside coef. limits"
            WARN(msg)
          End If
        End If
      End If

      If ( opts%n2o_Data .And. coefs%coef % nn2o > 0) Then
        ig = coefs%coef % fmv_gas_pos( gas_id_n2o )
        If  ( ( prof % n2o(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
           &  ( prof % n2o(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric n2o outside coef. limits"
            WARN(msg)
          End If
        End If
      End If

      If ( opts%ch4_Data .And. coefs%coef % nch4 > 0) Then
        ig = coefs%coef % fmv_gas_pos( gas_id_ch4 )
        If  ( ( prof % ch4(ilev) > coefs%coef % lim_prfl_gmax(jlev, ig) ) .or. &
           &  ( prof % ch4(ilev) < coefs%coef % lim_prfl_gmin(jlev, ig) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric ch4 outside coef. limits"
            WARN(msg)
          End If
        End If
      End If


    End Do
    
  ELSE IF(opts%addpc)then


    if( (prof % s2m % p < coefs%coef_pccomp % lim_pc_prfl_pmin).or. &
      & (prof % s2m % p > coefs%coef_pccomp % lim_pc_prfl_pmax)) Then
      ERR = errorstatus_warning
      If ( opts%verbose_checkinput_warnings ) Then
        msg = "PC-RTTOV: surface pressure outside limits"
        WARN(msg)
      End If
    endif

    if( (prof % s2m % t < coefs%coef_pccomp % lim_pc_prfl_tsmin).or. &
      & (prof % s2m % t > coefs%coef_pccomp % lim_pc_prfl_tsmax)) Then
      ERR = errorstatus_warning
      If ( opts%verbose_checkinput_warnings ) Then
        msg = "PC-RTTOV: surface temperature outside limits"
        WARN(msg)
      End If
    endif

    if( (prof % skin % t < coefs%coef_pccomp % lim_pc_prfl_skmin).or.   &
      & (prof % skin % t > coefs%coef_pccomp % lim_pc_prfl_skmax)) Then
      ERR = errorstatus_warning
      If ( opts%verbose_checkinput_warnings ) Then
        msg = "PC-RTTOV: skin temperature outside limits"
        WARN(msg)
      End If
    endif

    wind = sqrt(&
       & prof % s2m % u * prof % s2m % u + &
       & prof % s2m % v * prof % s2m % v   )

    if( (wind < coefs%coef_pccomp % lim_pc_prfl_wsmin).or.  &
      & (wind > coefs%coef_pccomp % lim_pc_prfl_wsmax)) Then
      ERR = errorstatus_warning
      If ( opts%verbose_checkinput_warnings ) Then
        msg = "PC-RTTOV: 10m wind speed outside limits"
        WARN(msg)
      End If
    endif


    jlev = 1
    Do ilev = 1, firstlevel
      ! recherche des niveaux RTTOV qui encadrent le niveau user, les pressions sont
      ! arondies 0.1Pa
      OK =  nint(prof % p(ilev)*1000) >= nint(coefs%coef % lim_prfl_p(jlev)*1000)   .and. &
         &  nint(prof % p(ilev)*1000) <  nint(coefs%coef % lim_prfl_p(jlev+1)*1000)

      Do while (  .not. OK .and. jlev < coefs%coef % nlevels)
        jlev = jlev+1
        If ( jlev < coefs%coef % nlevels ) Then
          OK =  nint(prof % p(ilev)*1000) >= nint(coefs%coef % lim_prfl_p(jlev)*1000)   .and. &
             &  nint(prof % p(ilev)*1000) <  nint(coefs%coef % lim_prfl_p(jlev+1)*1000)
        Else
          ! on est sur le dernier niveau des fichiers de coefs
          OK =  .true.
        End If
      End Do

      If  ( ( prof % t(ilev) > coefs%coef_pccomp % lim_pc_prfl_tmax(jlev) ) .or. &
         &  ( prof % t(ilev) < coefs%coef_pccomp % lim_pc_prfl_tmin(jlev) ) ) Then
        ERR = errorstatus_warning
        If ( opts%verbose_checkinput_warnings ) Then
          msg = "PC-RTTOV: some atmospheric temperature outside coef. limits"
          WARN(msg)
        End If
      End If

      If  ( ( prof % q(ilev) > coefs%coef_pccomp % lim_pc_prfl_qmax(jlev) ) .or. &
         &  ( prof % q(ilev) < coefs%coef_pccomp % lim_pc_prfl_qmin(jlev) ) ) Then
        ERR = errorstatus_warning
        If ( opts%verbose_checkinput_warnings ) Then
          msg = "PC-RTTOV: some atmospheric water vapour outside coef. limits"
          WARN(msg)
        End If
      End If

      If ( opts%ozone_Data .And. coefs%coef % nozone > 0) Then
        If  ( ( prof % o3(ilev) > coefs%coef_pccomp % lim_pc_prfl_ozmax(jlev) ) .or. &
           &  ( prof % o3(ilev) < coefs%coef_pccomp % lim_pc_prfl_ozmin(jlev) ) ) Then
          ERR = errorstatus_warning
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: some atmospheric ozone outside coef. limits"
            WARN(msg)
          End If
        End If
      End If

    End Do
 
  Endif


IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CATCH

IF (LHOOK) CALL DR_HOOK('RTTOV_USER_PROFILE_CHECKINPUT',1_jpim,ZHOOK_HANDLE)
End Subroutine rttov_user_profile_checkinput
