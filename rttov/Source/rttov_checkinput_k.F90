!
Subroutine rttov_checkinput_k( &
       & opts,        &
       & chanprof,    &
       & prof,        &
       & prof_k,      &
       & coef) ! out
  ! Description:
  !    K of resetting profile variables to regression limits
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
  !  1.0       01/12/2002 New code (N. Bormann)
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:

  ! Imported Parameters:




  ! Imported Type Definitions:
  Use rttov_types, Only : &
         & rttov_chanprof, &
         & rttov_coef     ,&
         & rttov_options,  &
         & profile_Type

!INTF_OFF
Use rttov_const, Only : &
    gas_id_watervapour,&
    gas_id_ozone,&
    gas_id_co2,&
    gas_id_co,&
    gas_id_n2o,&
    gas_id_ch4
Use yomhook, Only : &
    LHOOK,&
    DR_HOOK
Use parkind1, Only : &
    jpim,&
    jprb
!INTF_ON





  Implicit None


  ! subroutine arguments
  ! scalar arguments with intent(in):
  Type(rttov_options),Intent(in)  :: opts
  Type(rttov_chanprof),Intent(in) :: chanprof(:)
  Type(profile_Type), Intent (in) :: prof(:)    ! input profiles
  Type(profile_Type), Intent (inout) :: prof_k(:)    ! input profiles
  Type( rttov_coef ), Intent (in) :: coef    ! coefficients

!INTF_END

  !local variables:
  Integer(Kind=jpim) :: firstlevel,ilev
  Integer(Kind=jpim) :: ig          ! gas number
  Integer(kind=jpim) :: nchannels, ichan, iprof
REAL(KIND=JPRB) :: ZHOOK_HANDLE


  !- End of header --------------------------------------------------------

  !-------------
  !0. Initialize
  !-------------


  !determine first pressure level above the surface (note levels are top down)
IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT_K',0_jpim,ZHOOK_HANDLE)
nchannels = size(chanprof)

Do ichan = 1, nchannels
  iprof = chanprof(ichan)%prof

  Do firstlevel = coef % nlevels, 1, -1
     If ( coef % ref_prfl_p(firstlevel) <= prof(iprof) % s2m % p ) Exit
  End Do

  !-----------------------------
  !2. Check against basis values
  !-----------------------------

  If ( opts%apply_reg_limits ) Then
    Do ilev = 1, firstlevel
      If ( prof(iprof) % t(ilev) > coef % lim_prfl_tmax(ilev) ) Then
        prof_k(ichan) % t(ilev) = 0.0_JPRB
      Else If ( prof(iprof) % t(ilev) < coef % lim_prfl_tmin(ilev) ) Then
        prof_k(ichan) % t(ilev) = 0.0_JPRB
      End If

      ig = coef % fmv_gas_pos( gas_id_watervapour )
      If ( prof(iprof) % q(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
        prof_k(ichan) % q(ilev) = 0.0_JPRB
      Else If ( prof(iprof) % q(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
        prof_k(ichan) % q(ilev) = 0.0_JPRB
      End If
    End Do

    If ( opts%ozone_Data .And. coef % nozone > 0) Then
      ig = coef % fmv_gas_pos( gas_id_ozone )
      Do ilev = 1, firstlevel
        If ( prof(iprof) % o3(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
          prof_k(ichan) % o3(ilev) = 0.0_JPRB
        Else If ( prof(iprof) % o3(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
          prof_k(ichan) % o3(ilev) =  0.0_JPRB
        End If
      End Do
    End If

    If ( opts%co2_Data .And. coef % nco2 > 0) Then
      ig = coef % fmv_gas_pos( gas_id_co2 )
      Do ilev = 1, firstlevel
        If ( prof(iprof) % co2(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
          prof_k(ichan) % co2(ilev) =  0.0_JPRB
        Else If ( prof(iprof) % co2(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
          prof_k(ichan) % co2(ilev) =  0.0_JPRB
        End If
      End Do
    End If

    If ( opts%co_Data .And. coef % nco > 0) Then
      ig = coef % fmv_gas_pos( gas_id_co )
      Do ilev = 1, firstlevel
        If ( prof(iprof) % co(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
          prof_k(ichan) % co(ilev) =  0.0_JPRB
        Else If ( prof(iprof) % co(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
          prof_k(ichan) % co(ilev) =  0.0_JPRB
        End If
      End Do
    End If

    If ( opts%n2o_Data .And. coef % nn2o > 0) Then
      ig = coef % fmv_gas_pos( gas_id_n2o )
      Do ilev = 1, firstlevel
        If ( prof(iprof) % n2o(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
          prof_k(ichan) % n2o(ilev) =  0.0_JPRB
        Else If ( prof(iprof) % n2o(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
          prof_k(ichan) % n2o(ilev) =  0.0_JPRB
        End If
      End Do
    End If

    If ( opts%ch4_Data .And. coef % nch4 > 0) Then
      ig = coef % fmv_gas_pos( gas_id_ch4 )
      Do ilev = 1, firstlevel
        If ( prof(iprof) % ch4(ilev) > coef % lim_prfl_gmax(ilev, ig) ) Then
          prof_k(ichan) % ch4(ilev) = 0.0_JPRB
        Else If ( prof(iprof) % ch4(ilev) < coef % lim_prfl_gmin(ilev, ig) ) Then
          prof_k(ichan) % ch4(ilev) = 0.0_JPRB
        End If
      End Do
    End If

  End If
enddo

IF (LHOOK) CALL DR_HOOK('RTTOV_CHECKINPUT_K',1_jpim,ZHOOK_HANDLE)



End Subroutine rttov_checkinput_k
