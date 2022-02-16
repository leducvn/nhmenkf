!
Subroutine rttov_checkinput( &
       & opts,        &
       & prof,        &
       & coef,        &! in
       & coef_pccomp, &
       & ERR          ) ! out
  ! Description:
  ! Check input profile/angles
  ! (i)  Are physically realistic
  ! (ii) Profile values are within the basis set used to
  !      generate the coefficients
  ! Unphysical values return a fatal error status
  ! Profile values outside the basis set return a warning status
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
  !  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1       02/01/2003  More comments added (R Saunders)
  !  1.2       29/01/2003  More tests and add CO2 (P Brunel)
  !  1.3       27/06/2005  Uncommented water vapor and ozone profile checks  (R Saunders)
  !  1.4       23/01/2006  Changes from Marco (R Saunders)
  !  1.5       25/09/2007  Introduce a maximum number of warning messages (P Brunel)
  !  1.6       07/12/2007  Remove above; replace with a logical global variable... P.B.
  !  1.7       12/12/2007  Added limit checking for CH4, N2O, CO (R Saunders)
  !  1.8       16/01/2008  Added facility to apply regression limits (N Bormann)
  !  1.9       02/12/2009  Marco Matricardi: Added principal components
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
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

  ! Imported Parameters:




  ! Imported Type Definitions:
  Use rttov_types, Only : &
         & rttov_coef     ,&
         & rttov_options,  &
         & rttov_coef_pccomp,   &
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
    jprb
!INTF_ON



  Implicit None


  ! subroutine arguments
  ! scalar arguments with intent(in):
  Type(rttov_options),Intent(in)  :: opts
  Type(profile_Type), Intent (inout) :: prof(:)    ! input profiles
  Type( rttov_coef ), Intent (in)    :: coef    ! coefficients
  Type( rttov_coef_pccomp ), Intent (in)  :: coef_pccomp

  ! scalar arguments with intent(out):
  Integer(Kind=jpim), Intent (out) :: ERR       ! return code

!INTF_END

#include "rttov_errorreport.h"


  !local variables:
  Real(Kind=jprb)    :: wind
  Real(Kind=jprb)    :: dp( coef % nlevels )
  Integer(Kind=jpim) :: firstlevel,ilev
  Integer(Kind=jpim) :: ig          ! gas number
  Integer(kind=jpim) :: nprofiles, iprof
  character(32)      :: sprof
  character(128)     :: msg
REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
TRY

  !-------------
  !0. Initialize
  !-------------

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

nprofiles = size(prof)

Do iprof = 1, nprofiles
! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
  write(sprof,'(" (profile number = ",I8,")")') iprof
!$OMP END CRITICAL

  !determine first pressure level above the surface (note levels are top down)
  Do firstlevel = coef % nlevels, 1, -1
     If ( coef % ref_prfl_p(firstlevel) <= prof(iprof) % s2m % p ) Exit
  End Do

  ! Compare Profile Levels and Model Levels
  If ( prof(iprof) % nlevels /= coef % nlevels ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid profile number of levels"//sprof)
  End If
  
  dp(:) = abs ( coef % ref_prfl_p(:) - prof(iprof) % p(:) ) / coef % ref_prfl_p(:)
  If ( Any( dp > 0.01_JPRB ) ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid profile pressure levels"//sprof)
  End If
  

  !------------------------------
  !1. Check for unphysical values
  !------------------------------
  ! zenith angle
  If ( prof(iprof) % zenangle > zenmax .Or. &
     & prof(iprof) % zenangle < 0._JPRB          ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid zenith angle"//sprof)
  End If
  
  ! Cloud Fraction
  If ( prof(iprof) % cfraction > 1._JPRB .Or. &
     & prof(iprof) % cfraction < 0._JPRB      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid cloud fraction"//sprof)
  End If

  ! Cloud Top Pressure
  If ( prof(iprof) % cfraction .ne. 0 ) Then
    If ( prof(iprof) % ctp > ctpmax .Or. &
       & prof(iprof) % ctp < ctpmin      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid cloud top pressure"//sprof)
    End If
  End If
  
  ! Zeeman variables
  If (coef % inczeeman) Then
    ! Magnetic field strength
    If ( prof(iprof) % be > bemax .Or. &
       & prof(iprof) % be < bemin      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid magnetic field strength"//sprof)
    End If

    ! Cosine of angle between path and mag. field
    If ( prof(iprof) % cosbk > 1._JPRB  .Or. &
       & prof(iprof) % cosbk < -1._JPRB      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid cosbk"//sprof)
    End If    
  End If

  !1.1 surface variables
  !---------------------

  ! Pressure
  If ( prof(iprof) % s2m % p > pmax .Or. &
     & prof(iprof) % s2m % p < pmin      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid surface pressure"//sprof)
  End If

  ! 2m air temperature
  If ( prof(iprof) % s2m % t    > tmax .Or. &
     & prof(iprof) % s2m % t    < tmin ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid 2m air temperature"//sprof)
  End If

  ! 2m water vapour - only used if opts % use_q2m is TRUE
  If ( opts % use_q2m ) Then
    If ( prof(iprof) % s2m % q > qmax .Or. &
       & prof(iprof) % s2m % q <= qmin      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "invalid 2m water vapour"//sprof)
    End If
  End If

  !  surface wind speed
  wind = sqrt(&
          & prof(iprof) % s2m % u * prof(iprof) % s2m % u + &
          & prof(iprof) % s2m % v * prof(iprof) % s2m % v   )
  If ( wind > wmax .Or. &
     & wind < 0._JPRB      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid 10m wind speed"//sprof)
  End If

  ! surface skin temperature
  If ( prof(iprof) % skin % t > tmax .Or. &
     & prof(iprof) % skin % t < tmin      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "invalid skin surface temperature"//sprof)
  End If

  ! surface type
  If ( prof(iprof) % skin % surftype < 0 .Or. &
     & prof(iprof) % skin % surftype > nsurftype ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid surface type"//sprof)
  End If
  
  ! water type
  If ( prof(iprof) % skin % watertype < 0 .Or. &
     & prof(iprof) % skin % watertype > nwatertype ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid water type"//sprof)
  End If

  
  !1.2 atmospheric variables
  !-------------------------

  ! Predictors are calculated on *all* levels so check
  ! the hard limits on every level to avoid errors.

  ! temperature
  If ( Any( prof(iprof) % t(:) > tmax ) .Or. &
     & Any( prof(iprof) % t(:) < tmin )      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid atmospheric temperature"//sprof)
  End If

  ! water vapour
  If ( Any( prof(iprof) % q(:) > qmax ) .Or. &
     & Any( prof(iprof) % q(:) <= qmin )      ) Then
    ERR = ERRORSTATUS_FATAL
    THROWM( ERR .NE. 0 , "some invalid atmospheric water vapour"//sprof)
  End If

  ! ozone
  If ( opts%ozone_Data .And. coef % nozone > 0) Then
     If ( Any( prof(iprof) % o3(:) > o3max ) .Or. &
        & Any( prof(iprof) % o3(:) < o3min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric ozone"//sprof)
    End If
  End If

  ! CO2
  If ( opts%co2_Data .And. coef % nco2 > 0) Then
     If ( Any( prof(iprof) % co2(:) > co2max ) .Or. &
        & Any( prof(iprof) % co2(:) < co2min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CO2"//sprof)
    End If
  End If

  ! CO
  If ( opts%co_Data .And. coef % nco > 0) Then
     If ( Any( prof(iprof) % co(:) > comax ) .Or. &
        & Any( prof(iprof) % co(:) < comin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CO"//sprof)
    End If
  End If

  ! N2O
  If ( opts%n2o_Data .And. coef % nn2o > 0) Then
     If ( Any( prof(iprof) % n2o(:) > n2omax ) .Or. &
        & Any( prof(iprof) % n2o(:) < n2omin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric N2O"//sprof)
    End If
  End If

  ! CH4
  If ( opts%ch4_Data .And. coef % nch4 > 0) Then
     If ( Any( prof(iprof) % ch4(:) > ch4max ) .Or. &
        & Any( prof(iprof) % ch4(:) < ch4min )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid atmospheric CH4"//sprof)
    End If
  End If

  ! cloud liquid water
  If ( opts%clw_Data ) Then
     If ( Any( prof(iprof) % clw(:) > clwmax ) .Or. &
        & Any( prof(iprof) % clw(:) < clwmin )      ) Then
      ERR = ERRORSTATUS_FATAL
      THROWM( ERR .NE. 0 , "some invalid cloud liquid water"//sprof)
    End If
  Endif

  !-----------------------------
  !2. Check against basis values
  !-----------------------------


    if(.not.opts%addpc)then

      If ( opts%apply_reg_limits ) Then
        Do ilev = 1, firstlevel
          If ( prof(iprof) % t(ilev) > coef % lim_prfl_tmax(ilev) ) Then
            prof(iprof) % t(ilev) = coef % lim_prfl_tmax(ilev)
          Else If ( prof(iprof) % t(ilev) < coef % lim_prfl_tmin(ilev) ) Then
            prof(iprof) % t(ilev) = coef % lim_prfl_tmin(ilev)
          End If

          ig = coef % fmv_gas_pos( gas_id_watervapour )
          If ( prof(iprof) % q(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
            prof(iprof) % q(ilev) = coef % lim_prfl_gmax(ilev, ig)
          Else If ( prof(iprof) % q(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
            prof(iprof) % q(ilev) = coef % lim_prfl_gmin(ilev, ig)
          End If
        End Do

        If ( opts%ozone_Data .And. coef % nozone > 0) Then
          ig = coef % fmv_gas_pos( gas_id_ozone )
          Do ilev = 1, firstlevel
            If ( prof(iprof) % o3(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
              prof(iprof) % o3(ilev) = coef % lim_prfl_gmax(ilev, ig)
            Else If ( prof(iprof) % o3(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
              prof(iprof) % o3(ilev) = coef % lim_prfl_gmin(ilev, ig)
            End If
          End Do
        End If

        If ( opts%co2_Data .And. coef % nco2 > 0) Then
          ig = coef % fmv_gas_pos( gas_id_co2 )
          Do ilev = 1, firstlevel
            If ( prof(iprof) % co2(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
              prof(iprof) % co2(ilev) = coef % lim_prfl_gmax(ilev, ig)
            Else If ( prof(iprof) % co2(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
              prof(iprof) % co2(ilev) = coef % lim_prfl_gmin(ilev, ig)
            End If
          End Do
        End If

        If ( opts%co_Data .And. coef % nco > 0) Then
          ig = coef % fmv_gas_pos( gas_id_co )
          Do ilev = 1, firstlevel
            If ( prof(iprof) % co(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
              prof(iprof) % co(ilev) = coef % lim_prfl_gmax(ilev, ig)
            Else If ( prof(iprof) % co(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
              prof(iprof) % co(ilev) = coef % lim_prfl_gmin(ilev, ig)
            End If
          End Do
        End If

        If ( opts%n2o_Data .And. coef % nn2o > 0) Then
          ig = coef % fmv_gas_pos( gas_id_n2o )
          Do ilev = 1, firstlevel
            If ( prof(iprof) % n2o(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
              prof(iprof) % n2o(ilev) = coef % lim_prfl_gmax(ilev, ig)
            Else If ( prof(iprof) % n2o(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
              prof(iprof) % n2o(ilev) = coef % lim_prfl_gmin(ilev, ig)
            End If
          End Do
        End If

        If ( opts%ch4_Data .And. coef % nch4 > 0) Then
          ig = coef % fmv_gas_pos( gas_id_ch4 )
          Do ilev = 1, firstlevel
            If ( prof(iprof) % ch4(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
              prof(iprof) % ch4(ilev) = coef % lim_prfl_gmax(ilev, ig)
            Else If ( prof(iprof) % ch4(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
              prof(iprof) % ch4(ilev) = coef % lim_prfl_gmin(ilev, ig)
            End If
          End Do
        End If

      Else

        If ( Any( prof(iprof) % t(1:firstlevel) > coef % lim_prfl_tmax(1:firstlevel) ) .Or. &
              & Any( prof(iprof) % t(1:firstlevel) < coef % lim_prfl_tmin(1:firstlevel) ) ) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric temperature outside coef. limits"//sprof 
            WARN(msg)
          End If
        End If

        ig = coef % fmv_gas_pos( gas_id_watervapour )
!        If ( Any( prof(iprof) % q(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
!             & Any( prof(iprof) % q(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
        If ( Any( prof(iprof) % q(1:coef % nlevels) > coef % lim_prfl_gmax(1:coef % nlevels, ig) ) .Or. &
             & Any( prof(iprof) % q(1:coef % nlevels) < coef % lim_prfl_gmin(1:coef % nlevels, ig) ) ) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "some atmospheric water vapour outside coef. limits"//sprof
            WARN(msg)
          End If
        End If

        If ( opts%ozone_Data .And. coef % nozone > 0) Then
          ig = coef % fmv_gas_pos( gas_id_ozone )
          If ( Any( prof(iprof) % o3(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
                & Any( prof(iprof) % o3(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
               msg = "some atmospheric ozone outside coef. limits"//sprof
               WARN(msg)
             End If
          End If
        End If

        If ( opts%co2_Data .And. coef % nco2 > 0) Then
          ig = coef % fmv_gas_pos( gas_id_co2 )
          If ( Any( prof(iprof) % co2(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
                 & Any( prof(iprof) % co2(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
              msg = "some atmospheric CO2 outside coef. limits"//sprof
              WARN(msg)
            End If
          End If
        End If

        If ( opts%co_Data .And. coef % nco > 0) Then
          ig = coef % fmv_gas_pos( gas_id_co )
          If ( Any( prof(iprof) % co(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
                 & Any( prof(iprof) % co(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
              msg = "some atmospheric CO outside coef. limits"//sprof
              WARN(msg)
            End If
          End If
        End If

        If ( opts%n2o_Data .And. coef % nn2o > 0) Then
          ig = coef % fmv_gas_pos( gas_id_n2o )
          If ( Any( prof(iprof) % n2o(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
                 & Any( prof(iprof) % n2o(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
              msg = "some atmospheric N2O outside coef. limits"//sprof
              WARN(msg)
            End If
          End If
        End If

        If ( opts%ch4_Data .And. coef % nch4 > 0) Then
          ig = coef % fmv_gas_pos( gas_id_ch4 )
          If ( Any( prof(iprof) % ch4(1:firstlevel) > coef % lim_prfl_gmax(1:firstlevel, ig) ) .Or. &
                 & Any( prof(iprof) % ch4(1:firstlevel) < coef % lim_prfl_gmin(1:firstlevel, ig) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
              msg = "some atmospheric CH4 outside coef. limits"//sprof
              WARN(msg)
            End If
          End If
        End If

      End If
    else if(opts%addpc)then
      If ( opts%apply_reg_limits ) Then
        Do ilev = 1, firstlevel

          If ( prof(iprof) % t(ilev) > coef_pccomp % lim_pc_prfl_tmax(ilev) ) Then
            prof(iprof) % t(ilev) = coef_pccomp % lim_pc_prfl_tmax(ilev)
          Else If ( prof(iprof) % t(ilev) < coef_pccomp % lim_pc_prfl_tmin(ilev) ) Then
            prof(iprof) % t(ilev) = coef_pccomp % lim_pc_prfl_tmin(ilev)
          End If

          If ( prof(iprof) % q(ilev) > coef_pccomp % lim_pc_prfl_qmax(ilev) ) Then
            prof(iprof) % q(ilev) = coef_pccomp % lim_pc_prfl_qmax(ilev)
          Else If ( prof(iprof) % q(ilev) < coef_pccomp % lim_pc_prfl_qmin(ilev) ) Then
            prof(iprof) % q(ilev) = coef_pccomp % lim_pc_prfl_qmin(ilev)
          End If
        End Do

        If ( opts%ozone_Data .And. coef % nozone > 0) Then
          Do ilev = 1, firstlevel
            If ( prof(iprof) % o3(ilev) > coef_pccomp % lim_pc_prfl_ozmax(ilev) ) Then
              prof(iprof) % o3(ilev) = coef_pccomp % lim_pc_prfl_ozmax(ilev)
            Else If ( prof(iprof) % o3(ilev) < coef_pccomp % lim_pc_prfl_ozmin(ilev) ) Then
              prof(iprof) % o3(ilev) = coef_pccomp % lim_pc_prfl_ozmin(ilev)
            End If
          End Do
        End If

        if( prof(iprof) % s2m % p < coef_pccomp % lim_pc_prfl_pmin) Then
          prof(iprof) % s2m % p=coef_pccomp % lim_pc_prfl_pmin
        else if(prof(iprof) % s2m % p > coef_pccomp % lim_pc_prfl_pmax) Then
          prof(iprof) % s2m % p=coef_pccomp % lim_pc_prfl_pmax
        endif

        if( prof(iprof) % s2m % t < coef_pccomp % lim_pc_prfl_tsmin) Then
          prof(iprof) % s2m % t=coef_pccomp % lim_pc_prfl_tsmin
        else if(prof(iprof) % s2m % t > coef_pccomp % lim_pc_prfl_tsmax) Then
          prof(iprof) % s2m % t=coef_pccomp % lim_pc_prfl_tsmax
        endif

        if( prof(iprof) % skin % t < coef_pccomp % lim_pc_prfl_skmin) Then
          prof(iprof) % skin % t=coef_pccomp % lim_pc_prfl_skmin
        else if(prof(iprof) % skin % t > coef_pccomp % lim_pc_prfl_skmax) Then
          prof(iprof) % skin % t=coef_pccomp % lim_pc_prfl_skmax
        endif

        wind = sqrt(&
          & prof(iprof) % s2m % u * prof(iprof) % s2m % u + &
          & prof(iprof) % s2m % v * prof(iprof) % s2m % v   )

        if( wind < coef_pccomp % lim_pc_prfl_wsmin) Then
          prof(iprof) % s2m % u=sqrt(coef_pccomp % lim_pc_prfl_wsmin**2/2._jprb)
          prof(iprof) % s2m % v=sqrt(coef_pccomp % lim_pc_prfl_wsmin**2/2._jprb)
        else if(wind > coef_pccomp % lim_pc_prfl_wsmax) Then
          prof(iprof) % s2m % u=sqrt(coef_pccomp % lim_pc_prfl_wsmax**2/2._jprb)
          prof(iprof) % s2m % v=sqrt(coef_pccomp % lim_pc_prfl_wsmax**2/2._jprb)
        endif

      else
        If ( Any( prof(iprof) % t(1:firstlevel) > coef_pccomp % lim_pc_prfl_tmax(1:firstlevel) ) .Or. &
              & Any( prof(iprof) % t(1:firstlevel) < coef_pccomp % lim_pc_prfl_tmin(1:firstlevel) ) ) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: some atmospheric temperature outside coef. limits"//sprof
            WARN(msg)
          End If
        End If

        If ( Any( prof(iprof) % q(1:firstlevel) > coef_pccomp % lim_pc_prfl_qmax(1:firstlevel) ) .Or. &
              & Any( prof(iprof) % q(1:firstlevel) < coef_pccomp % lim_pc_prfl_qmin(1:firstlevel) ) ) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
             msg = "PC-RTTOV: some atmospheric water vapour outside coef. limits"//sprof
             WARN(msg)
           End If
        End If

        If ( opts%ozone_Data .And. coef % nozone > 0) Then
           If ( Any( prof(iprof) % o3(1:firstlevel) > coef_pccomp % lim_pc_prfl_ozmax(1:firstlevel) ) .Or. &
                 & Any( prof(iprof) % o3(1:firstlevel) < coef_pccomp % lim_pc_prfl_ozmin(1:firstlevel) ) ) Then
            ERR = errorstatus_warning 
            If ( opts%verbose_checkinput_warnings ) Then
                msg = "PC-RTTOV: some atmospheric ozone outside coef. limits" //sprof
                WARN(msg)
              End If
           End If
        End If

        if( (prof(iprof) % s2m % p < coef_pccomp % lim_pc_prfl_pmin).or. &
        & (prof(iprof) % s2m % p > coef_pccomp % lim_pc_prfl_pmax)) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: surface pressure outside limits"//sprof
            WARN(msg)
          End If
        endif

        if( (prof(iprof) % s2m % t < coef_pccomp % lim_pc_prfl_tsmin).or. &
            & (prof(iprof) % s2m % t > coef_pccomp % lim_pc_prfl_tsmax)) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: surface temperature outside limits"//sprof
            WARN(msg)
          End If
        endif

        if( (prof(iprof) % skin % t < coef_pccomp % lim_pc_prfl_skmin).or.   &
        & (prof(iprof) % skin % t > coef_pccomp % lim_pc_prfl_skmax)) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: skin temperature outside limits"//sprof
            WARN(msg)
          End If
        endif

        wind = sqrt(&
          & prof(iprof) % s2m % u * prof(iprof) % s2m % u + &
          & prof(iprof) % s2m % v * prof(iprof) % s2m % v   )

        if( (wind < coef_pccomp % lim_pc_prfl_wsmin).or.  &
        & (wind > coef_pccomp % lim_pc_prfl_wsmax)) Then
          ERR = errorstatus_warning 
          If ( opts%verbose_checkinput_warnings ) Then
            msg = "PC-RTTOV: 10m wind speed outside limits"//sprof
            WARN(msg)
          End If
        endif

      endif

    Endif
enddo

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

CATCH

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT',1_jpim,ZHOOK_HANDLE)
End Subroutine rttov_checkinput
