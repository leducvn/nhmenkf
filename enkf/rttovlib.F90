module rttovlib
! Author: Le Duc
! Created date: 10 Apr 2016
   use variable, only : r_size, r_dble
   use rttov_types
   use parkind1, only : jpim, jprb
   !
   integer(4),parameter :: jppf   = 1        ! Max no. profiles per RTTOV call
   integer(4),parameter :: jpnsat = 60       ! Max no. of used sensors
   integer(4),parameter :: jpchus =100       ! Max no. of used channels
   integer(4),parameter :: jplevo = 68
   integer(4),parameter :: jplev  = 40       ! No. of pressure levels        ! org
   integer(4),parameter :: jpnav  =  4       ! No. of profile variables       ! org
   integer(4),parameter :: jpnsav =  5       ! No. of surface air variables   ! org
   integer(4),parameter :: jpnssv =  6       ! No. of skin variables            ! modified
   integer(4),parameter :: jpncv  =  2       ! No. of cloud variables
   integer(4),parameter :: jplenpf=jplev*jpnav+jpnsav+jpnssv+jpncv+jpchus
   integer(4),parameter :: jpchpf = jppf*jpchus ! Max no. of profs * chans used
   integer(4),parameter :: rttv_chan_hyperch = 8461
   !
   integer(4),parameter :: &
     &    id1000_tovs =15,         &
     &    id1000_mwpw =17,         &
     &    id1000_gtbb =21,         &
     &    idtype_gpspw=14100,      &
     &    idtype_gleo =14200,      &
     &    idtype_tovs =15200,      &! TOVS
     &    idtype_rtovs=15210,      &! RTOVS
     &    idtype_atovs=15215,      &! ATOVS(L1D,HIRS+AMSUA)
     &    idtype_amsub=15300,      &! AMSUB(L1D/B)
     &    idtype_amsua=15310,      &! AMSUA(L1B)
     &    idtype_hirs =15320,      &! HIRS (L1B)
     &    idtype_airs =15330,      &! AIRS
     &    idtype_iasi =15410,      &! IASI
     &    idtype_cris =15510,      &! CrIS
     &    idtype_mhs  =15350,      &! MHS  (L1B)
     &    idtype_atms =15500,      &! ATMS
     &    idtype_mwhs =15600,      &! FY3 MWHS
     &    idtype_mwts =15610,      &! FY3 MWTS
     &    idtype_mwhs2=15620,      &! FY3 MWHS2
     &    idtype_mwts2=15630,      &! FY3 MWTS2
     &    idtype_saphir=15700,     &! SAPHIR
     &    idtype_ssmi =17100,      &! SSM/I
     &    idtype_tmi  =17200,      &! TMI
     &    idtype_gmi  =17250,      &! GMI
     &    idtype_amsr =17300,      &! AMSR
     &    idtype_amse =17350,      &! AMSR-E
     &    idtype_amsr2=17360,      &! AMSR2
     &    idtype_ssmis_env =17400, &! SSMIS-ENV
     &    idtype_ssmis_img =17410, &! SSMIS-IMG
     &    idtype_ssmis_uas =17420, &! SSMIS-UAS
     &    idtype_ssmis_las =17430, &! SSMIS-LAS
     &    idtype_ssmis_uk =17440,  &! SSMIS-UK_LAS
     &    idtype_ssmis_lasenv =17445,  &! SSMIS-LASENV
     &    idtype_ssmis_upp    =17451,  &! SSMIS-UPP(DCDH)
     &    idtype_windsat =17510,   &! WindSat
     &    idtype_mwri =17600        ! FY3 MWRI
   integer(4),parameter :: &
     &    idtype_goes    =21500,   &! GOES
     &    idtype_meteosat=21200,   &! METEOSAT
     &    idtype_mtsat   =21300,   &! MTSAT
     &    idtype_msg     =21400,   &! MSG
     &    idtype_himawari=21600,   &! Himawari
     &    idtype_mtg     =21700     ! MTG
   !
   integer(4), save :: rttv_chidx(jpchus+1,jpnsat)
   integer(4), save :: rttv_ch2378tidx(rttv_chan_hyperch*jppf,jpnsat)  ! convert org.ch index to Channel Index
   !
   type type_rt
      integer(kind=jpim) :: n_coef
      type(rttov_coefs), pointer :: coefs(:)
      type(rttov_scatt_coef), pointer :: coef_scatt(:)
      type(rttov_options) :: opts
   end type type_rt
   type(type_rt), save :: rt
   !
   type type_rt_io
      integer(kind=jpim), pointer :: errorstatus(:) ! nprofiles
      integer(kind=jpim) :: nprofiles
      integer(kind=jpim) :: nchannels
      type(profile_type), pointer :: profiles(:) ! nprofiles
      type(profile_type), pointer :: profiles_k(:)
      type(profile_type), pointer :: profiles_ad(:)
      !
      type(profile_cloud_type), pointer :: cld_profiles(:) ! nprofiles
      type(profile_cloud_type), pointer :: cld_profiles_k(:)
      type(profile_cloud_type), pointer :: cld_profiles_ad(:)
      !
      logical, pointer :: calcemiss(:) ! nchannels
      real(kind=jprb), pointer :: emissivity_in(:) ! nchannels
      real(kind=jprb), pointer :: emissivity_in_k(:) ! nchannels
      real(kind=jprb), pointer :: emissivity_out(:) ! nchannels
      real(kind=jprb), pointer :: emissivity_out_k(:) ! nchannels
      real(kind=jprb), pointer :: emissivity_ad(:) ! nchannels
      !
      type(transmission_type) :: transmission
      type(transmission_type) :: transmission_k
      type(transmission_type) :: transmission_ad
      type(radiance_type) :: radiance
      type(radiance_type) :: radiance_k
      type(radiance_type) :: radiance_ad
   end type type_rt_io
   !
   !
   !
contains
   !
   !
   !
   subroutine initialize_rttov(myid)
      use rttov_const, only :         &
     & errorstatus_success,           &
     & nplatforms,ninst,sensor_id_mw, &
     & inst_id_goesim,inst_id_gmsim,  &
     & platform_name,inst_name,       &
     & mh2o, mair
      implicit none
      !
      integer, intent(in) :: myid
      integer(kind=jpim), allocatable :: input_instrument(:) ! instrument id
      integer(kind=jpim), allocatable :: input_channels(:)
      integer(kind=jpim) :: input_errorstatus
      integer(Kind=jpim), parameter  :: nprofiles = 1                ! Number of profiles
      ! AD variables for rttov_k calls
      integer(kind=jpim), allocatable :: errorstatus(:)  ! rttov error return code
      integer(kind=jpim), allocatable :: setup_errorstatus(:) ! setup return code
      !
      integer, parameter :: ERR_UNIT=6
      integer, parameter :: VERBOSITY_LEVEL=1 ! 0: no error messages output
                                              ! 1: fatal errors only printed
                                              ! 2: warning erros only printed
                                              ! 3: information message
      !
      integer :: num_sensor      ! number of the sensors used
      integer :: r_satn          ! satn : satellite number
      integer :: r_nidx          ! nidx : number of channel index
      character(len=8) :: r_plat ! plat : platform
      character(len=8) :: r_inst ! inst : instrument
      integer :: tvsinst(3), r_cidx(10)
      integer, allocatable :: tvschan(:)
      integer(4) :: ich, chidx(jpchus+1)
      integer :: i_coef, c_plat, c_inst, c_chan,ctchan
      integer :: md
      !integer :: rttv_varbc, rttv_varqc, var_on
      !character(len=256) :: dmmy
      !
#include "rttov_errorhandling.h"
#include "rttov_setup.h"
      !
      rt%opts%addinterp = .true.
      rt%opts%addsolar  = .false.
      rt%opts%addaerosl = .false.
      rt%opts%addclouds = .false.
      rt%opts%switchrad =.true.
      rt%opts%spacetop  =.true.
      rt%opts%lgradp    =.false.
      rt%opts%use_q2m   = .false.
      rt%opts%apply_reg_limits = .true.
      rt%opts%verbose_checkinput_warnings = .true.
      rt%opts%ozone_data    = .false.
      rt%opts%co2_data      = .false.
      rt%opts%n2o_data      = .false.
      rt%opts%co_data       = .false.
      rt%opts%ch4_data      = .false.
      rt%opts%clw_data      = .true.
      rt%opts%addrefrac     = .false.
      rt%opts%do_checkinput = .true.
      !
      md = 99
      open(md, file="satellites.txt")
      read(md,'(i2)') num_sensor
      if (myid == 0) write(6,'(2x,a,i4)') 'num_sensor=',num_sensor
      allocate(rt%coefs(num_sensor))
      allocate(input_instrument(3))
      call rttov_errorhandling(err_unit, verbosity_level)
      input_errorstatus = 0
      rttv_chidx(:,:) = 0
      rttv_ch2378tidx(:,:) = 0
      !
      do i_coef = 1,num_sensor
         read(md,'(a8,i2,2x,a8,i3)') r_plat,r_satn,r_inst,r_nidx
         ! ch number to be processed: 0 is equivalent to jpchus, or all available ch are processed.
         ! The negative value allows smaller # of ch than those set in preprocessing (QCRTX) to be processed
	 chidx(jpchus+1) = r_nidx  
         if(r_nidx<0) r_nidx=-r_nidx
         ! platform id
         do c_plat = 1, nplatforms
            if (r_plat == platform_name(c_plat)) tvsinst(1) = c_plat
         end do
         ! satellite number id
         tvsinst(2)=r_satn
         ! instrument id
         if (r_inst == 'imager  ') then
            select case (r_plat)
            case('goes    ')
               tvsinst(3)=inst_id_goesim
            case('mtsat   ')
               tvsinst(3)=inst_id_gmsim
            end select
         else
            do c_inst = 0,ninst-1
               if(r_inst==inst_name(c_inst)) tvsinst(3)=c_inst
            end do
         end if
         !
         if (r_nidx == 0) then ! non.hyperCH
            do ich=1,jpchus
               chidx(ich) = ich
            end do
            rttv_chidx(1:(jpchus+1),i_coef) = chidx(1:(jpchus+1))
         else ! hyperCH
            allocate(tvschan(r_nidx))
            do ctchan = 0, (r_nidx-1)/10
               read(md,'(10i5)') r_cidx
               do c_chan = 1, min(10,r_nidx-ctchan*10)
                  tvschan(ctchan*10+c_chan) = r_cidx(c_chan)
               end do
            end do
            allocate(input_channels(r_nidx))
            input_channels(1:r_nidx) = tvschan(1:r_nidx)
            chidx(1:r_nidx) = tvschan(1:r_nidx)
            rttv_chidx(1:(jpchus+1),i_coef) = chidx(1:(jpchus+1))
            do ich = 1, r_nidx
               rttv_ch2378tidx(chidx(ich),i_coef) = ich
            end do
            deallocate(tvschan)
         end if ! hyperCH or not
         !
         input_instrument(1) = tvsinst(1)
         input_instrument(2) = tvsinst(2)
         input_instrument(3) = tvsinst(3)
         if (myid == 0) then
            if (r_nidx /= 0)then
               write(6,'(i4,2x,a8,i2,2x,a8,4i4)') i_coef, r_plat,r_satn,r_inst,input_instrument(1:3),r_nidx
               write(6,'(22i5)') input_channels(1:r_nidx)
            else
               write(6,'(i4,2x,a8,i2,2x,a8,4i4)') i_coef, r_plat,r_satn,r_inst,input_instrument(1:3),0
            endif
         end if
         !
         if (r_nidx /= 0)then
            call rttov_setup (&
          & input_errorstatus,        &! out
          & err_unit,                 &! in
          & verbosity_level,          &! in
          & rt%opts,                  &! in
          & rt%coefs(i_coef),         &! out
          & input_instrument,         &! in
          & channels = input_channels , & ! in Optional 
          & channels_rec = input_channels  ) ! in Optional 
            deallocate(input_channels)
         else
            call rttov_setup (&
          & input_errorstatus,        &! out
          & err_unit,                 &! in
          & verbosity_level,          &! in
          & rt%opts,                  &! in
          & rt%coefs(i_coef),         &! out
          & input_instrument )         ! in
         endif
         if (input_errorstatus /= errorstatus_success ) then 
            write ( 6, * ) 'rttov_setup fatal error'
            stop 
         endif
      end do  ! i_coef=1,num_sensor
      !
      if (myid == 0) then
         do i_coef = 1,num_sensor
            write(6,*) &
        &   i_coef, &
        &   rt%coefs(i_coef)%coef%id_platform, &
        &   rt%coefs(i_coef)%coef%id_sat, &
        &   rt%coefs(i_coef)%coef%id_inst, &
        &   input_instrument(1:3),rt%coefs(i_coef)%coef%fmv_chn
            write(6,'(22i5)') &
        &   rt%coefs(i_coef)%coef%ff_ori_chn(1:rt%coefs(i_coef)%coef%fmv_chn)
            write(6,'(22i5)') &
        &   rt%coefs(i_coef)%coef%ff_val_chn(1:rt%coefs(i_coef)%coef%fmv_chn)
            !if (rt%coefs(i_coef)%coef%id_sensor == sensor_id_mw) then !MW
               !write(6,'(22i5)') rt%coefs(i_coef)%coef%fastem_polar(:)
               !write(6,'(i5)') rt%coefs(i_coef)%coef%fastem_ver
            !end if
            write(6,*) "rt%coefs(i_coef)%coef%fmv_model_ver = ", i_coef,rt%coefs(i_coef)%coef%fmv_model_ver
         end do  ! i_coef=1,num_sensor
      end if
      close(md)
      !write(6,*) "PARM_RTTOV FILE READ END"
      rt%n_coef = num_sensor
      deallocate(input_instrument)
      !
      return
   end subroutine initialize_rttov
   !
   !
   !
   subroutine finalize_rttov()
      implicit none
      deallocate(rt%coefs)
      !deallocate(rt%coef_scatt)
      return
   end subroutine finalize_rttov
   !
   !
   !
   subroutine get_satid(idname, idtype, ksat)
      use rttov_const, only : platform_id_goes,&
     &  platform_id_meteosat,&
     &  platform_id_msg, platform_id_mtsat,platform_id_noaa,&
     &  platform_id_mtg, platform_id_himawari,&
     &  platform_id_dmsp,&
     &  platform_id_eos, platform_id_trmm,platform_id_metop,&
     &  platform_id_coriolis, platform_id_fy3, platform_id_jpss,&
     &  platform_id_gcomw, platform_id_gpm, &
     &  platform_id_meghatr, &
     &  inst_id_goesim,inst_id_mviri,inst_id_seviri,&
     &  inst_id_gmsim,&
     &  inst_id_irs, inst_id_ahi,&
     &  inst_id_hirs, inst_id_amsua, inst_id_amsub,&
     &  inst_id_mhs, inst_id_airs, inst_id_ssmi, inst_id_tmi,&
     &  inst_id_gmi, inst_id_saphir,&
     &  inst_id_amsr, inst_id_ssmis, inst_id_iasi, inst_id_cris,&
     &  inst_id_windsat, inst_id_mwts, inst_id_mwhs, inst_id_mwri, &
     &  inst_id_mwts2, inst_id_mwhs2, &
     &  inst_id_atms
      implicit none
      !
      character(len=8), intent(in) :: idname
      integer(4), intent(in) :: idtype
      integer(4), intent(out) :: ksat
      integer :: plat, satn, inst, isat
      !
      ksat= 0
      plat=0
      satn=0
      inst=0
      !
      ! Platform and satellite number
      if(     idname(1:4)=='0253')then
      plat=platform_id_goes
      satn=9
      else if(idname(1:4)=='0254')then
      plat=platform_id_goes
      satn=10
      else if(idname(1:4)=='0255')then
      plat=platform_id_goes
      satn=11
      else if(idname(1:4)=='0256')then
      plat=platform_id_goes
      satn=12
      else if(idname(1:4)=='0257')then
      plat=platform_id_goes
      satn=13
      else if(idname(1:4)=='0258')then
      plat=platform_id_goes
      satn=14
      else if(idname(1:4)=='0259')then
      plat=platform_id_goes
      satn=15
      !
      else if(idname(1:4)=='0052')then
      plat=platform_id_meteosat
      satn=5
      else if(idname(1:4)=='0054')then
      plat=platform_id_meteosat
      satn=7
      else if(idname(1:4)=='0055')then
      plat=platform_id_msg
      satn=1
      else if(idname(1:4)=='0056')then
      plat=platform_id_msg
      satn=2
      else if(idname(1:4)=='0057')then
      plat=platform_id_msg
      satn=3
      else if(idname(1:4)=='0070')then
      plat=platform_id_msg
      satn=4
      else if(idname(1:4)=='0071')then
      plat=platform_id_mtg
      satn=1
      !
      else if(idname(1:4)=='0171')then
      plat=platform_id_mtsat
      satn=1
      else if(idname(1:4)=='0172')then
      plat=platform_id_mtsat
      satn=2
      else if(idname(1:4)=='0173')then
      plat=platform_id_himawari
      satn=8
      else if(idname(1:4)=='0174')then
      plat=platform_id_himawari
      satn=9
      !
      else if(idname(1:4)=='0206')then
      plat=platform_id_noaa
      satn=15
      else if(idname(1:4)=='0207')then
      plat=platform_id_noaa
      satn=16
      else if(idname(1:4)=='0208')then
      plat=platform_id_noaa
      satn=17
      else if(idname(1:4)=='0209')then
      plat=platform_id_noaa
      satn=18
      else if(idname(1:4)=='0223')then
      plat=platform_id_noaa
      satn=19
      !
      else if(idname(1:4)=='0246')then
      plat=platform_id_dmsp
      satn=13
      else if(idname(1:4)=='0247')then
      plat=platform_id_dmsp
      satn=14
      else if(idname(1:4)=='0248')then
      plat=platform_id_dmsp
      satn=15
      else if(idname(1:4)=='0249')then
      plat=platform_id_dmsp
      satn=16
      else if(idname(1:4)=='0285')then
      plat=platform_id_dmsp
      satn=17
      else if(idname(1:4)=='0286')then
      plat=platform_id_dmsp
      satn=18
      else if(idname(1:4)=='0287')then
      plat=platform_id_dmsp
      satn=19
      !
      else if(idname(1:4)=='0783')then
      plat=platform_id_eos
      satn=1
      else if(idname(1:4)=='0784')then
      plat=platform_id_eos
      satn=2
      else if(idname(1:4)=='0122')then
      plat=platform_id_gcomw
      satn=1
      else if(idname(1:4)=='0288')then
      plat=platform_id_gpm
      satn=1
      else if(idname(1:4)=='0440')then
      plat=platform_id_meghatr
      satn=1
      !
      else if(idname(1:4)=='TRMM')then
      plat=platform_id_trmm
      satn=1
      else if(idname(1:4)=='0003')then
      plat=platform_id_metop
      satn=1
      else if(idname(1:4)=='0004')then
      plat=platform_id_metop
      satn=2
      else if(idname(1:4)=='0005')then
      plat=platform_id_metop
      satn=3
      !
      else if(idname(1:4)=='0520')then
      plat=platform_id_fy3
      satn=1
      else if(idname(1:4)=='0521')then
      plat=platform_id_fy3
      satn=2
      else if(idname(1:4)=='0522')then
      plat=platform_id_fy3
      satn=3
      else if(idname(1:4)=='0523')then
      plat=platform_id_fy3
      satn=4
      else if(idname(1:4)=='0283')then
      plat=platform_id_coriolis
      satn=1
      else if(idname(1:4)=='0224')then
      plat=platform_id_jpss
      satn=0
      end if
      !
      ! Instrument
      if(idtype==idtype_goes)          then
      inst=inst_id_goesim
      else if(idtype==idtype_meteosat) then
      inst=inst_id_mviri
      else if(idtype==idtype_msg)      then
      inst=inst_id_seviri
      else if(idtype==idtype_mtg)      then
      inst=inst_id_irs
      else if(idtype==idtype_mtsat)    then
      inst=inst_id_gmsim
      else if(idtype==idtype_himawari) then
      inst=inst_id_ahi
      !
      else if(idtype==idtype_hirs)  then
      inst=inst_id_hirs
      else if(idtype==idtype_amsua) then
      inst=inst_id_amsua
      else if(idtype==idtype_amsub) then
      inst=inst_id_amsub
      else if(idtype==idtype_mhs)   then
      inst=inst_id_mhs
      else if(idtype==idtype_airs)  then
      inst=inst_id_airs
      else if(idtype==idtype_iasi)  then
      inst=inst_id_iasi
      else if(idtype==idtype_cris)  then
      inst=inst_id_cris
      else if(idtype==idtype_ssmi)  then
      inst=inst_id_ssmi
      else if(idtype==idtype_tmi)   then
      inst=inst_id_tmi
      else if(idtype==idtype_gmi)   then
      inst=inst_id_gmi
      else if(idtype==idtype_saphir)   then
      inst=inst_id_saphir
      else if(idtype==idtype_amse)  then
      inst=inst_id_amsr
      else if(idtype==idtype_amsr2)  then
      inst=inst_id_amsr
      else if(idtype==idtype_windsat)  then
      inst=inst_id_windsat
      else if(idtype==idtype_mwts)  then
      inst=inst_id_mwts
      else if(idtype==idtype_mwhs)  then
      inst=inst_id_mwhs
      else if(idtype==idtype_mwts2) then
      inst=inst_id_mwts2
      else if(idtype==idtype_mwhs2) then
      inst=inst_id_mwhs2
      else if(idtype==idtype_mwri)  then
      inst=inst_id_mwri
      else if((idtype==idtype_ssmis_env) .or. &
     &       (idtype==idtype_ssmis_img) .or. &
     &       (idtype==idtype_ssmis_uas) .or. &
     &       (idtype==idtype_ssmis_las) .or. &
     &       (idtype==idtype_ssmis_uk)  .or. &
     &       (idtype==idtype_ssmis_upp) .or. &
     &       (idtype==idtype_ssmis_lasenv)) then
      inst=inst_id_ssmis
      else if(idtype==idtype_atms)  then
      inst=inst_id_atms
      end if
      !
      ! RTTOV INDEX CHECK
      do isat=1,rt%n_coef
         if (rt%coefs(isat)%coef%id_platform == plat .and. &
     &       rt%coefs(isat)%coef%id_sat == satn       .and. &
     &       rt%coefs(isat)%coef%id_inst == inst ) then
            ksat=isat
            return
         end if
      end do
      ksat=0
      !
      return
   end subroutine get_satid
   !
   !
   !
end module rttovlib
