!
module rttov_types
  ! Description:
  ! defines all derived types for RTTOV
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
  !  1.0   01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1   29/01/2003  Add CO2 variable gaz to profile sturcture (P Brunel)
  !                    Add rain and solid precip. to profile cloud structure
  !  1.2   13/05/2003  Add structure for transmissions and optical depths (F Chevallier)
  !  1.3      08/2003  Add scattering facility (F Chevallier)
  !  1.4   18/09/2003  Add kice and kradip to profile_cloud_type (P Francis)
  !  1.5   09/12/2003  Change type for mclayer to INTEGER (R Saunders)
  !  1.6   06/01/2004  Add CO2 to ref profile (R Saunders)
  !  1.7   02/06/2004  Add fast model version compatibility level in coef type (P. Brunel)
  !  1.8   17/05/2005  Add q to profile_cloud_type ( U O'Keeffe)
  !  1.9   12/01/2006  Marco Matricardi (ECMWF):
  !           --       Added raytracing_type structure  for computation of variable
  !           --       local zenith angle.
  !           --       Added N2O,CO and CH4 to profile_type structure.
  !           --       Modified type rttov_coef.
  ! 1.10   30/01/2007  Removed frequency tagged variables R Saunders
  ! 1.11   15/03/2007  Added single stream transmittances for RTTOV
  !                    o/p R Saunders
  ! 1.12   11/10/2007  Move iaernum & iaertyp profile members to profile_aux
  !                    P.Marguinaud
  ! 1.13   05/07/2007  Added extra radiance outputs from cloud (R Saunders)
  ! 1.14   22/11/2007  RTTOV-SCATT version 9 - remove cld_radiance_type; slim
  !                    down profile_cloud_type ( A Geer )
  ! 1.15   12/12/2007  add IncTop (R Saunders)
  ! 1.16   12/08/2008  add IncZeeman, Be, cosbk (P. Rayer)
  ! 1.16   03/10/2008  RTTOV-SCATT revised cloud partitioning (A Geer)
  ! 1.17   01/09/2009  Added salinity for FASTEM-4 (Mark Liu)
  ! 1.19   03/11/2009  Transmittances / optical depths on levels (A Geer)
  ! 1.20   02/12/2009  Added principal component capability (Marco Matricard)
  ! 1.21   02/12/2009  Introduced new variables for the mixed cloud scheme (Marco Matricardi)
  ! 1.22   17/06/2010  Introduced spacetop flag to zero opdeps at user's
  !                    model-top in Zeeman channels (P Rayer)
  ! 1.23   05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
  ! 1.24   14/10/2010  Remove rt8_mode (J Hocking)
  !
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
  use rttov_const, only: &
        & fastem_sp,      &
        & ncldtyp
  Use parkind1, Only : jpim, jprb, jplm
  Implicit None

  ! Surface skin
  Type sskin_type
  Integer(Kind=jpim) :: surftype        ! 0=land, 1=sea, 2=sea-ice
  Integer(Kind=jpim) :: watertype       ! 0=fresh water, 1=ocean water
  Real(Kind=jprb)    :: t               ! radiative skin temperature (K)
  Real(Kind=jprb)    :: salinity        ! practical salinity unit %o
  Real(Kind=jprb)    :: fastem(fastem_sp)  ! land/sea-ice surface parameters for fastem-2
  End Type sskin_type

  ! Surface 2m
  Type s2m_type
  Real(Kind=jprb) :: t                  ! temperature (K)
  Real(Kind=jprb) :: q                  ! water vapour (ppmv)
  Real(Kind=jprb) :: o                  ! ozone (ppmv)
  Real(Kind=jprb) :: p                  ! surface pressure (hPa)
  Real(Kind=jprb) :: u                  ! U 10m wind component (m/s)
  Real(Kind=jprb) :: v                  ! V 10m wind component (m/s)
  Real(Kind=jprb) :: wfetc              ! Wind fetch (metres)
  End Type s2m_type


  ! structure for atmospheric profiles on model pressure levels
  Type profile_type
  Character(Len=128) :: id
  Integer(Kind=jpim) :: date(3)  ! Year, Month, Day
  Integer(Kind=jpim) :: time(3)  ! Hour, Minute, Second
     ! number of atmospheric levels/layers
  Integer(Kind=jpim) :: nlevels
  Integer(Kind=jpim) :: nlayers
     ! ozone, CO2 and cloud liquid water profiles available
     ! atmosphere defined on nlevels
  Real(Kind=jprb), Pointer :: p(:)      ! pressure (hPa)
  Real(Kind=jprb), Pointer :: t(:)      ! temperature (K)
  Real(Kind=jprb), Pointer :: q(:)      ! water vapour (ppmv)
  Real(Kind=jprb), Pointer :: o3(:)     ! ozone (ppmv)
  Real(Kind=jprb), Pointer :: co2(:)    ! carbon dioxide (ppmv)
  Real(Kind=jprb), Pointer :: n2o(:)    ! n2o(ppmv)
  Real(Kind=jprb), Pointer :: co(:)     ! co(ppmv)
  Real(Kind=jprb), Pointer :: ch4(:)    ! ch4(ppmv)
  Real(Kind=jprb), Pointer :: clw(:)    ! cloud liquid water (kg/kg)
  Real(Kind=jprb), Pointer :: aerosols(:,:)
  Real(Kind=jprb), Pointer :: cloud(:,:)
  Real(Kind=jprb), Pointer :: cfrac(:,:)
  Real(Kind=jprb), Pointer :: icede(:)  ! ice particle effective diameter (microns)
  Integer(Kind=jpim)       :: idg
  Integer(Kind=jpim)       :: ish
     ! surface
     Type(sskin_type) :: skin
     Type(s2m_type)   :: s2m
     !angles
  Real(Kind=jprb) :: zenangle
  Real(Kind=jprb) :: azangle
  Real(Kind=jprb) :: sunzenangle
  Real(Kind=jprb) :: sunazangle
  Real(Kind=jprb) :: elevation
  Real(Kind=jprb) :: latitude
  Real(Kind=jprb) :: longitude
  Real(Kind=jprb) :: snow_frac     ! snow coverage fraction for IR emissivity atlas (0 - 1)
  Real(Kind=jprb) :: soil_moisture ! soil moisture (m^3/m^3)
  ! Be - Earth magnetic field strength (Gauss)
  ! cosbk - cosine of the angle between the Earth magnetic field and wave
  !         propagation direction
  Real(Kind=jprb) :: Be
  Real(Kind=jprb) :: cosbk
     ! Black body cloud
  Real(Kind=jprb) :: ctp               ! cloud top pressure  (hPa)
  Real(Kind=jprb) :: cfraction         ! cloud fraction (0 - 1) 1 for 100% cloud cover
  End Type profile_type

  ! structure for atmospheric profile additional input for RTTOV-SCATT,
  ! with information on clouds and rain for each level. Full level pressures,
  ! t, q etc. should be placed in the profile_type.

  Type profile_cloud_type

  Integer(Kind=jpim) :: nlevels ! number of atmospheric levels (nlevels+1 for ph)
  logical(Kind=jplm) :: use_totalice ! False => separate ice and snow  True => total ice
  Real(Kind=jprb) :: cfrac      ! Average cloud fraction (only used if lusercfrac = .true.)

  Real(Kind=jprb), Pointer :: ph(:)        ! nlevels+1 of half-level model pressures (hPa)
  Real(Kind=jprb), Pointer :: cc(:)        ! nlevels of cloud cover
  Real(Kind=jprb), Pointer :: clw(:)       ! nlevels of cloud liquid water (kg/kg)
  Real(Kind=jprb), Pointer :: ciw(:)       ! nlevels of cloud ice water (kg/kg)
  Real(Kind=jprb), Pointer :: totalice(:)  ! nlevels of total ice (kg/kg)
  Real(Kind=jprb), Pointer :: rain(:)      ! nlevels of rain (kg/m2/s)
  Real(Kind=jprb), Pointer :: sp(:)        ! nlevels of solid precipitation (kg/m2/s)

  End Type profile_cloud_type

  ! satellite geometry
  Type geometry_type
  Real(Kind=jprb) :: sinzen
  Real(Kind=jprb) :: sinzen_sq
  Real(Kind=jprb) :: coszen
  Real(Kind=jprb) :: coszen_sq
  Real(Kind=jprb) :: seczen
  Real(Kind=jprb) :: seczen_sq
  Real(Kind=jprb) :: seczen_sqrt
  Real(Kind=jprb) :: seczen_minus1
  Real(Kind=jprb) :: seczen_minus1_sq
  Real(Kind=jprb) :: sinview
  Real(Kind=jprb) :: sinview_sq
  Real(Kind=jprb) :: cosview_sq
  Real(Kind=jprb) :: normzen
  Real(Kind=jprb) :: viewang
  End Type geometry_type

  ! Predictors
  Type predictors_type
     ! the nxxxx could be set to 0 to indicate the abscence
     ! of the predictor, in that case there is no need to
     ! allocate the corresponding predictor
  Integer(Kind=jpim) :: nlevels   ! number of levels for predictors (all same)
  Integer(Kind=jpim) :: nmixed    ! number of variables for Mixed Gases
  Integer(Kind=jpim) :: nwater    ! number of variables for Water Vapour
  Integer(Kind=jpim) :: nozone    ! number of variables for Ozone
  Integer(Kind=jpim) :: nwvcont   ! number of variables for WV Continuum
  Integer(Kind=jpim) :: nco2      ! number of variables for CO2
  Integer(Kind=jpim) :: nn2o      ! number of variables for N2O
  Integer(Kind=jpim) :: nco       ! number of variables for CO
  Integer(Kind=jpim) :: nch4      ! number of variables for CH4
  Integer(Kind=jpim) :: ncloud    ! number of variables for MW Cloud
  Real(Kind=jprb), Pointer     :: mixedgas(:,:,:)          ! (nmixed,  nlevels, nprofiles )
  Real(Kind=jprb), Pointer     :: watervapour(:,:,:)       ! (nwater,  nlevels, nprofiles)
  Real(Kind=jprb), Pointer     :: ozone(:,:,:)             ! (nozone,  nlevels)
  Real(Kind=jprb), Pointer     :: wvcont(:,:,:)            ! (nwvcont, nlevels)
  Real(Kind=jprb), Pointer     :: co2(:,:,:)               ! (nco2,    nlevels)
  Real(Kind=jprb), Pointer     :: n2o(:,:,:)               ! (nn2o,    nlevels)
  Real(Kind=jprb), Pointer     :: co(:,:,:)                ! (nco,     nlevels)
  Real(Kind=jprb), Pointer     :: ch4(:,:,:)               ! (nch4,    nlevels)
  Real(Kind=jprb), Pointer     :: clw(:,:)                 ! (         nlevels)
  Real(Kind=jprb), Pointer     :: mixedgas_sun(:,:,:)      ! (nmixed,  nlevels)
  Real(Kind=jprb), Pointer     :: watervapour_sun(:,:,:)   ! (nwater,  nlevels)
  Real(Kind=jprb), Pointer     :: ozone_sun(:,:,:)         ! (nozone,  nlevels)
  Real(Kind=jprb), Pointer     :: wvcont_sun(:,:,:)        ! (nwvcont, nlevels)
  Real(Kind=jprb), Pointer     :: co2_sun(:,:,:)           ! (nco2,    nlevels)
  Real(Kind=jprb), Pointer     :: n2o_sun(:,:,:)           ! (nn2o,    nlevels)
  Real(Kind=jprb), Pointer     :: co_sun(:,:,:)            ! (nco,     nlevels)
  Real(Kind=jprb), Pointer     :: ch4_sun(:,:,:)           ! (nch4,    nlevels)
  End Type predictors_type


  Type rttov_coef
     ! Structure for the storage of RTTOV coefficients
     ! this may differ from what is stored in the coefficient files especially
     ! for the units (ie kg/kg to ppmv)
     ! Gases are separated in MxG WV O3
     ! Number of levels is the same for all gases (taken from MxG).
     !
  Integer(Kind=jpim) :: id_platform  ! platform   (see documentation or MOD_CPARAM)
  Integer(Kind=jpim) :: id_sat    ! satellite  (.....)
  Integer(Kind=jpim) :: id_inst    ! instrument (.....)
  Integer(Kind=jpim) :: id_sensor  ! sensor
     !  1 = Infrared
     !  2 = Micro Wave
     !  3 = High resolution
  Integer(Kind=jpim) :: id_comp_lvl  ! RTTOV coefficient file version number
  Integer(Kind=jpim) :: id_comp_pc   ! Principal component coefficient file version number
  Integer(Kind=jpim) ,Dimension(3) :: id_creation_date  ! YYYY MM DD
  Character (len=80)    :: id_creation    ! Creation comment
  Character (len=32)    :: id_Common_name  ! usual name of the satellite
  Character (len=80)    :: line_by_line(20)


     !FAST_MODEL_VARIABLES section
  Character (len=32)    :: fmv_model_def  ! FMV definition (RTTOV6 OPTRAN RTTOV7)
  Integer(Kind=jpim)               :: fmv_model_ver  ! fast model version compatibility level
  Integer(Kind=jpim)               :: fmv_chn        ! number of channels in file
  Integer(Kind=jpim)               :: fmv_gas        ! number of gases in file
  Integer(Kind=jpim), pointer      :: fmv_gas_id(:)    ! gas id. number i gas_id list (fmv_gas)
  Integer(Kind=jpim), Pointer      :: fmv_gas_pos(:)   ! respective position of each gas of gas_id list (ngases_max)
  Integer(Kind=jpim), Pointer      :: fmv_var(:)       ! number of variables/predictors by gaz (fmv_gas)
  Integer(Kind=jpim), Pointer      :: fmv_coe(:)       ! number of coefficients by gaz (fmv_gas)
  Integer(Kind=jpim), Pointer      :: fmv_int(:)       ! number of spectral intervals by gaz (fmv_gas)
  Integer(Kind=jpim), Pointer      :: fmv_lvl(:)       ! number of levels(pres/absorber) by gaz (fmv_gas)

  Integer(Kind=jpim)               :: nmixed         ! number of variables/predictors for Mixed Gases
  Integer(Kind=jpim)               :: nwater         ! number of variables/predictors for Water Vapour
  Integer(Kind=jpim)               :: nozone         ! number of variables/predictors for Ozone
  Integer(Kind=jpim)               :: nwvcont        ! number of variables/predictors for WV continuum
  Integer(Kind=jpim)               :: nco2           ! number of variables/predictors for CO2
  Integer(Kind=jpim)               :: nn2o           ! number of variables/predictors for N2O
  Integer(Kind=jpim)               :: nco            ! number of variables/predictors for CO
  Integer(Kind=jpim)               :: nch4           ! number of variables/predictors for CH4
  Integer(Kind=jpim)               :: nlevels        ! number of levels(pres/absorber) same for all gases
  Integer(Kind=jpim)               :: nlayers        ! number of layers(pres/absorber) nlevels-1
  Logical(Kind=jplm)               :: IncZeeman      ! Flag to include Zeeman effect for this sensor
  Integer(Kind=jpim)               :: ncmixed        ! number of coefficients for Mixed Gases
  Integer(Kind=jpim)               :: ncwater        ! number of coefficients for Water Vapour
  Integer(Kind=jpim)               :: ncozone        ! number of coefficients for Ozone
  Integer(Kind=jpim)               :: ncwvcont       ! number of coefficients for WV continuum
  Integer(Kind=jpim)               :: ncco2          ! number of coefficients for CO2
  Integer(Kind=jpim)               :: ncn2o          ! number of coefficients for N2O
  Integer(Kind=jpim)               :: ncco           ! number of coefficients for CO
  Integer(Kind=jpim)               :: ncch4          ! number of coefficients for CH4
  Integer(Kind=jpim)               :: nintmixed
  Integer(Kind=jpim)               :: nintwater
  Integer(Kind=jpim)               :: nintozone
  Integer(Kind=jpim)               :: nintwvcont
  Integer(Kind=jpim)               :: nintco2
  Integer(Kind=jpim)               :: nintn2o
  Integer(Kind=jpim)               :: nintco
  Integer(Kind=jpim)               :: nintch4


     !GAZ_UNITS section
     ! gases are in the order of gas id codes
  Integer(Kind=jpim), Pointer      :: gaz_units(:)  ! unit of gaz concentration for each gaz
                                             ! default value is specific conc. (kg/kg)
                                             ! value inside RTTOV calculations (ppmv)
     !FILTER_FUNCTIONS section  array size is fmv_chn
  Integer(Kind=jpim) ,Pointer :: ff_ori_chn(:)   ! original chan number
  Integer(Kind=jpim) ,Pointer :: ff_val_chn(:)   ! validity of the channel (1=OK)
  Real(Kind=jprb) ,Pointer :: ff_cwn (:)         ! cental wave number (cm-1)
  Real(Kind=jprb) ,Pointer :: ff_bco (:)         ! band correction offset (K)
  Real(Kind=jprb) ,Pointer :: ff_bcs (:)         ! band correction slope (K/K)
  Real(Kind=jprb) ,Pointer :: ff_gam (:)         ! gamma factor transm. correction

     !TRANSMITTANCE_TRESHOLD section  array size is fmv_chn
  Integer(Kind=jpim) ,Pointer :: tt_chn(:)
  Integer(Kind=jpim) ,Pointer :: tt_val_chn(:)
  Real(Kind=jprb)    ,Pointer :: tt_cwn (:)
  Real(Kind=jprb)    ,Pointer :: tt_a0(:)
  Real(Kind=jprb)    ,Pointer :: tt_a1(:)

     !SOLAR_SPECTRUM section array size is fmv_chn
  Integer(Kind=jpim) ,Pointer :: ss_chn(:)
  Integer(Kind=jpim) ,Pointer :: ss_val_chn(:)
  Real(Kind=jprb)    ,Pointer :: ss_cwn (:)
  Real(Kind=jprb)    ,Pointer :: ss_solar_spectrum(:)

     !WATER_OPTICAL_CONSTANT section array size is fmv_chn
  Integer(Kind=jpim) ,Pointer :: woc_chn(:)
  Real(Kind=jprb)    ,Pointer :: woc_cwn (:)
  Complex(Kind=jprb) ,Pointer :: woc_waopc_ow(:)
  Complex(Kind=jprb) ,Pointer :: woc_waopc_fw(:)

     !WAVE_SPECTRUM section array size is ws_nomega
     !Data used to compute the frequency spectrum of the JONSWAP
     !wave model surface wave.
  Integer(Kind=jpim)          :: ws_nomega
  Real(Kind=jprb)  ,Pointer   :: ws_npoint(:)
  Real(Kind=jprb)  ,Pointer   :: ws_k_omega(:)

     !FUNDAMENTAL_CONSTANTS section
  Real(Kind=jprb) :: fc_speedl         ! speed of light (cm/s)
  Real(Kind=jprb) :: fc_planck_c1      ! first radiation constant (mW/(m2*sr*cm-4))
  Real(Kind=jprb) :: fc_planck_c2      ! second radiation constant (cm*K)
  Real(Kind=jprb) :: fc_sat_height     ! satellite nominal altitude (km)

     !FASTEM section
  Integer(Kind=jpim) :: fastem_ver      ! fastem version number
  Integer(Kind=jpim), Pointer :: fastem_polar(:)  ! polarisation of each channel
     ! 0 = 0.5 V+H
     ! 1 = 90 - incident angle
     ! 2 = incident angle
     ! 3 = vertical
     ! 4 = horizontal
     ! 5 = V+H
     ! Full stokes vector

     !SSIREM section     array size is fmv_chn
     ! ems =   ssirem_a0
     !       - ssirem_a1*(zen**ssirem_xzn1)
     !       - ssirem_a2*(zen**ssirem_xzn2)
     ! where zen is satellite zenith angle in degrees, divided by 60.
  Integer(Kind=jpim) :: ssirem_ver                ! version number
  Integer(Kind=jpim),  Pointer  :: ssirem_chn(:)   ! original chan number
  Real(Kind=jprb),  Pointer     :: ssirem_a0(:)    ! constant coef
  Real(Kind=jprb),  Pointer     :: ssirem_a1(:)    ! first coef
  Real(Kind=jprb),  Pointer     :: ssirem_a2(:)    ! second coef
  Real(Kind=jprb),  Pointer     :: ssirem_xzn1(:)  ! 1st exponent on zenith angle
  Real(Kind=jprb),  Pointer     :: ssirem_xzn2(:)  ! 2nd exponent on zenith angle

     !REFERENCE_PROFILE section  defined on Mixed gases pressure levels
     ! Not working for OPTRAN gas absorber levels
     ! gases are in the order of gas id codes
     ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
     ! unit for mr for optical depth calculations is ppmv
  Real(Kind=jprb), Pointer      :: ref_prfl_p(:)     ! pressure  (hPa)       (levels)
  Real(Kind=jprb), Pointer      :: ref_prfl_t(:,:)   ! temperature (K)       (levels, gases)
  Real(Kind=jprb), Pointer      :: ref_prfl_mr(:,:)  ! mixing ratio (ppmv)   (levels, gases)
     !PROFILE_LIMITS section
     ! gases are in the order of gas id codes
     ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
     ! unit for mr for optical depth calculations is ppmv
  Real(Kind=jprb), Pointer      :: lim_prfl_p(:)       ! pressure  (hPa)       (levels)
  Real(Kind=jprb), Pointer      :: lim_prfl_tmax(:)    ! max temperature (K)   (levels)
  Real(Kind=jprb), Pointer      :: lim_prfl_tmin(:)    ! min temperature (K)   (levels)
  Real(Kind=jprb), Pointer      :: lim_prfl_gmax(:,:)  ! max mixing r (ppmv) (levels, gases)
  Real(Kind=jprb), Pointer      :: lim_prfl_gmin(:,:)  ! min mixing r (ppmv) (levels, gases)


     !FAST_COEFFICIENTS section
     ! separate arrays to allow different number of variables for each gaz
  Real(Kind=jprb), Pointer      :: mixedgas(:,:,:)     ! Mixed gases coefs  (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: watervapour(:,:,:)  ! Water vapour coefs (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: ozone(:,:,:)        ! Ozone coefs        (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: wvcont(:,:,:)       ! WV Cont coefs      (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: co2(:,:,:)          ! CO2 coefs          (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: n2o(:,:,:)          ! N2O coefs          (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: co(:,:,:)           ! CO coefs           (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: ch4(:,:,:)          ! CH4 coefs          (levels, channels, variables)
  Real(Kind=jprb), Pointer      :: mixedgasint(:,:)
  Real(Kind=jprb), Pointer      :: watervapourint(:,:)
  Real(Kind=jprb), Pointer      :: ozoneint(:,:)
  Real(Kind=jprb), Pointer      :: wvcontint(:,:)
  Real(Kind=jprb), Pointer      :: co2int(:,:)
  Real(Kind=jprb), Pointer      :: n2oint(:,:)
  Real(Kind=jprb), Pointer      :: coint(:,:)
  Real(Kind=jprb), Pointer      :: ch4int(:,:)

     ! Auxillary variables
  Real(Kind=jprb)               :: ratoe       ! ratio (H+R)/R  H=sat height, R=Earth radius
  Integer(Kind=jpim)            :: mwcldtop    ! Upper layer for MW LWP calcs
  Real(Kind=jprb), pointer      :: planck1(:)        ! C1 * Nu**3
  Real(Kind=jprb), pointer      :: planck2(:)        ! C2 * Nu
  Real(Kind=jprb), pointer      :: frequency_ghz(:)  ! frequency in GHz

     ! other predictor variables see Science and Validation report
  Real(Kind=jprb), pointer      :: dp(:)        ! interval between standard p levels (hPa)
  Real(Kind=jprb), pointer      :: dpp(:)       ! pressure based variable (hPa**2)
  Real(Kind=jprb), pointer      :: tstar(:)     ! layer temp (K)
  Real(Kind=jprb), pointer      :: to3star(:)   ! layer temp for O3 calculations (K)
  Real(Kind=jprb), pointer      :: wstar(:)     ! layer WV  (ppmv)
  Real(Kind=jprb), pointer      :: ostar(:)     ! layer O3  (ppmv)
  Real(Kind=jprb), pointer      :: co2star(:)   ! layer co2 (ppmv)
  Real(Kind=jprb), pointer      :: n2ostar(:)   ! layer n2o (ppmv)
  Real(Kind=jprb), pointer      :: costar(:)    ! layer co  (ppmv)
  Real(Kind=jprb), pointer      :: ch4star(:)   ! layer ch4 (ppmv)
  End Type rttov_coef

  Type rttov_scatt_coef
     ! Structure for the storage of RTTOV_SCATT coefficients
  Integer(Kind=jpim) :: nhydro ! Number of hydrometeors in computation
  Integer(Kind=jpim) :: mtype  ! Number of hydrometeors     in Mie tables
  Integer(Kind=jpim) :: mfreqm ! Number of frequencies      in Mie tables
  Integer(Kind=jpim) :: mtemp  ! Number of temperature bins in Mie tables
  Integer(Kind=jpim) :: mwc    ! Number of water bins       in Mie tables
  Real(Kind=jprb)    :: offset_temp_rain       ! temperature offset in table for rain type
  Real(Kind=jprb)    :: offset_temp_sp         ! temperature offset in table for solid prec. type
  Real(Kind=jprb)    :: offset_temp_liq        ! temperature offset in table for cloud water type
  Real(Kind=jprb)    :: offset_temp_ice        ! temperature offset in table for cloud ice type
  Real(Kind=jprb)    :: offset_temp_totalice   ! temperature offset in table for total ice type
  Real(Kind=jprb)    :: offset_water           ! liquid/ice water offset in table
  Real(Kind=jprb)    :: scale_water            ! log10(liquid/ice water) scaling factor in table
  Real(Kind=jprb)    :: from_scale_water       ! 10**(1._JPRB/scale_water)
  Real(Kind=jprb)    :: conv_rain(2)           ! coefficients for rain unit conversion (mm.h-1 to g.m-3)
  Real(Kind=jprb)    :: conv_sp  (2)           ! coefficients for solid prec. unit conversion (mm.h-1 to g.m-3)
  Real(Kind=jprb)    :: conv_liq (2)           ! coefficients for cloud water conversion (not used)
  Real(Kind=jprb)    :: conv_ice (2)           ! coefficients for cloud ice conversion   (not used)
  Real(Kind=jprb)    :: conv_totalice (2)      ! coefficients for total ice conversion   (not used)
  Real(Kind=jprb), pointer :: mie_freq(:)      ! list of frequencies in Mie table
  Real(Kind=jprb), pointer :: ext(:,:,:,:)     ! extinction coefficent table
  Real(Kind=jprb), pointer :: ssa(:,:,:,:)     ! single scattering albedo table
  Real(Kind=jprb), pointer :: asp(:,:,:,:)     ! assymetry parameter table

  End Type rttov_scatt_coef


  Type profile_aux_s
  Integer(Kind=jpim) :: nearestlev_surf ! nearest model level above surface
  Real(Kind=jprb)    :: pfraction_surf  ! pressure fraction of surface in model layer (hPa)
  Integer(Kind=jpim) :: nearestlev_ctp  ! nearest model level above cloud top
  Real(Kind=jprb) :: pfraction_ctp   ! pressure fraction of cloud top pressure in layer (hPa)
  Real(Kind=jprb) :: cfraction       ! cloud fraction (0 - 1) 1 for 100% cloud cover
  End Type

  ! Auxillary profile variables
  ! variables calculated by the model from profile
  type profile_aux
  type(profile_aux_s), pointer :: s(:)
  Real(Kind=jprb), pointer :: debye_prof(:,:,:)  ! Debye terms
  Real(Kind=jprb), pointer :: relhum(:,:)        !Relative humidity
  Real(Kind=jprb), pointer :: relhumref(:,:)
  Real(Kind=jprb), pointer :: dg(:,:)            !Generalized effective diameter
  Real(Kind=jprb), pointer :: fac1_dg(:,:)       !Intermediate variables used to compute the
  Real(Kind=jprb), pointer :: fac2_dg(:,:)       !generalized diameter
  Real(Kind=jprb), pointer :: fac3_dg(:,:)
  Integer(Kind=jpim), Pointer :: iaertyp(:,:,:)
  Integer(Kind=jpim), Pointer :: iaernum(:,:)
  end type profile_aux

  ! Auxillary profile variables for RTTOV_SCATT
  ! variables calculated by the model from profile
  Type profile_scatt_aux
  Real(Kind=jprb), pointer :: cfrac(:)        ! horizontal cloud fraction (one value used for all layers)
  Real(Kind=jprb), pointer :: ems_bnd(:)      ! surface emissivity for boundary conditions
  Real(Kind=jprb), pointer :: ref_bnd(:)      ! surface emissivity for boundary conditions
  Real(Kind=jprb), pointer :: ems_cld(:)      ! surface emissivity taking into account cloud/rain impact on od
  Real(Kind=jprb), pointer :: ref_cld(:)      ! surface reflectivity taking into account cloud/rain impact on od
  Real(Kind=jprb), pointer :: dz(:,:)         ! layer depth   [km]
  Real(Kind=jprb), pointer :: tbd(:,:)        ! temperature at layer boundary [K]
  Real(Kind=jprb), Pointer :: clw(:,:)        ! cloud liquid water (g/m3)
  Real(Kind=jprb), Pointer :: ciw(:,:)        ! cloud ice water (g/m3)
  Real(Kind=jprb), Pointer :: totalice(:,:)   ! total ice (g/m3)
  Real(Kind=jprb), Pointer :: rain(:,:)       ! rain (g/m3)
  Real(Kind=jprb), Pointer :: sp(:,:)         ! solid precipitation (g/m3)
!RWS  Real(Kind=jprb), pointer :: mclayer(:)  ! upper level cloud layer
  Integer(Kind=jpim), pointer :: mclayer(:)   ! upper level cloud layer
  Real(Kind=jprb), pointer :: delta(:,:)      ! (= ext*dz/coszen)
  Real(Kind=jprb), pointer :: tau(:,:)        ! optical depths (= exp(-delta))
  Real(Kind=jprb), pointer :: ext(:,:)        ! extinction coefficient integreated over hydrometeor types
  Real(Kind=jprb), pointer :: ssa(:,:)        ! single scattering albedo integreated over hydrometeor types
  Real(Kind=jprb), pointer :: asm(:,:)        ! asymetry parameter integreated over hydrometeor types [-1,1]
  Real(Kind=jprb), pointer :: lambda(:,:)     ! eddington approx. variable
                                  ! (= sqrt( 3*ext*ext*(1-ssa)*(1-ssa*asm) )
  Real(Kind=jprb), pointer :: h (:,:)         ! boundary condition variable (= 1.5_JPRB*ext(1-ssa*asm))
  Real(Kind=jprb), pointer :: b0(:,:)         ! lower level temperature
  Real(Kind=jprb), pointer :: b1(:,:)         ! temperature gradient
  Real(Kind=jprb), pointer :: bn(:,:)         ! upper level temperature
  end type profile_scatt_aux

  type opdp_path_type
     ! path optical depths as predicted or interpolated (unitless)
  Real(Kind=jprb), pointer :: atm_level(:,:)   ! neg optical depth for thermal radiation (levels to space), size (levels, channels)
  Real(Kind=jprb), pointer :: sun_level(:,:)   ! neg optical depth for solar radiation (levels to space), size (levels, channels)
  end type opdp_path_type

  type transmission_type
  Real(Kind=jprb), pointer  :: tau_total(:)       ! transmittance from surface (array size is of size nchannels)
  Real(Kind=jprb), pointer  :: tau_levels(:,:)    ! transmittance from each standard pressure level array (levels,channels)
  end type

  type transmission_type_aux
     ! Transmissions and optical depths (unitless)
  Real(Kind=jprb), pointer  :: fac(:,:,:,:)        ! Mask for integration calculation ! leading dimension size 2!
  Real(Kind=jprb), pointer  :: surf_fac(:,:)       ! Mask for integration calculation on surface
  Real(Kind=jprb), pointer  :: tau_surf(:,:)       ! transmittance from surface (array size is of size streams,nchannels)
  Real(Kind=jprb), pointer  :: tau_surf_r(:,:)     ! reciprocal transmittance from surface (array size is of size streams,nchannels)
  Real(Kind=jprb), pointer  :: tau_level(:,:,:)    ! transmittance from each standard pressure level array (levels,streams,channels)
  Real(Kind=jprb), pointer  :: tau_level_r(:,:,:)  ! reciprocal transmittance from each standard pressure level array
  Real(Kind=jprb), pointer  :: od_singlelayer(:,:,:)          ! single-layer optical depth
  Real(Kind=jprb), pointer  :: od_singlelayer_r(:,:,:)        ! reciprocal single-layer optical depth
  Real(Kind=jprb), pointer  :: od_level(:,:,:)                ! op dep from each standard pressure level array 
                                                              ! (levels,streams,channels)
  Real(Kind=jprb), pointer  :: odsun_singlelayer(:,:,:)
  Real(Kind=jprb), pointer  :: tausun_surf(:,:)               ! transmittance from surface (array size is of size nchannels)
  Real(Kind=jprb), pointer  :: tausun_level(:,:,:)            ! transmittance from each standard pressure level
                                                    !   (array size is of size (nlevels,nchannels))
  Real(Kind=jprb), pointer  :: od_sfrac(:,:)
  Real(Kind=jprb), pointer  :: od_sfrac_r(:,:)
  Real(Kind=jprb), pointer  :: odsun_sfrac(:,:)
  Real(Kind=jprb), pointer  :: od_frac_ac(:,:)
  Real(Kind=jprb), pointer  :: odsun_frac_ac(:,:)
  Real(Kind=jprb), pointer  :: tau_surf_ac(:,:)
  Real(Kind=jprb), pointer  :: tau_surf_acsun(:,:)
  Real(Kind=jprb), pointer  :: tau_ref_surf_ac(:,:)
  Real(Kind=jprb), pointer  :: tau_ref_surf_acsun(:,:)
  Real(Kind=jprb), pointer  :: od_frac_t(:,:)
  Real(Kind=jprb), pointer  :: odsun_frac_t(:,:)
  Real(Kind=jprb), pointer  :: tau_surf_t(:,:)
  Real(Kind=jprb), pointer  :: tausun_surf_t(:,:)
  Real(Kind=jprb), pointer  :: tau_ref_surf_t(:,:)
  Real(Kind=jprb), pointer  :: tausun_ref_surf_t(:,:)

  Real(kind=jprb), pointer  :: refl_norm(:)
  Real(Kind=jprb) :: anynegtau ! used to store information about presence of any negative transmittances
  end type transmission_type_aux

  type radiance_type
     ! Radiance and corresponding brightness temperature
     ! Array size is of size nchannels
     ! except for cloudy calculations (nlevels, nchannels)
     ! unit for radiance is mw/cm-1/ster/sq.m
     ! unit for temperature is Kelvin
     !
  Real(Kind=jprb), pointer  :: clear(:)       ! clear sky radiance
  Real(Kind=jprb), pointer  :: cloudy(:)      ! 100% cloudy radiance for given cloud
  Real(Kind=jprb), pointer  :: total(:)       ! cloudy radiance for given cloud
  Real(Kind=jprb), pointer  :: bt(:)          ! Brightness temp equivalent to total radiance
  Real(Kind=jprb), pointer  :: bt_clear(:)    ! Brightness temp equivalent to clear radiance
  Real(Kind=jprb), pointer  :: upclear(:)     ! clear sky radiance without reflection term
  Real(Kind=jprb), pointer  :: dnclear(:)     ! clear sky downwelling radiance
  Real(Kind=jprb), pointer  :: reflclear(:)   ! reflected clear sky downwelling radiance
  Real(Kind=jprb), pointer  :: overcast(:,:)  ! overcast radiance at given cloud
                                                       !   top  (levels,channels)
  Real(Kind=jprb), pointer  :: up(:,:)        ! sum( B * dT ) above cloud upwelling radiance for each pressure level
  Real(Kind=jprb), pointer  :: down(:,:)      ! sum ( B / T**2 dT ) above cloud downwelling radiance for each pressure level
  Real(Kind=jprb), pointer  :: surf(:,:)      !radiance at surface emitted from a black cloud (levels,channels).
  end type radiance_type

  type radiance_aux
     ! auxillary calculation arrays for RTE integration
     ! Direct model arrays need to be passed to TL AD and K codes
     ! array size is of (nchannels) or (nlevels, nchannels)
  Real(Kind=jprb), pointer :: air(:,:)
  Real(Kind=jprb), pointer :: surfair(:)
  Real(Kind=jprb), pointer :: skin(:)
  Real(Kind=jprb), pointer :: cosmic(:)
  Real(Kind=jprb), pointer :: up(:,:,:)               ! sum( B * dT )
  Real(Kind=jprb), pointer :: down(:,:,:)             ! sum ( B / T**2 dT )
  Real(Kind=jprb), pointer :: meanrad_up(:,:) 
  Real(Kind=jprb), pointer :: meanrad_down(:,:)
  Real(Kind=jprb), pointer :: down_ref(:,:,:)
  Real(Kind=jprb), pointer :: FAC1_1(:,:)
  Real(Kind=jprb), pointer :: FAC2_1(:,:)
  Real(Kind=jprb), pointer :: FAC3_1(:,:)
  Real(Kind=jprb), pointer :: FAC4_1(:,:)
  Real(Kind=jprb), pointer :: FAC5_1(:,:)
  Real(Kind=jprb), pointer :: FAC6_1(:,:)
  Real(Kind=jprb), pointer :: FAC7_1(:,:)
  Real(Kind=jprb), pointer :: FAC1_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC2_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC3_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC4_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC5_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC6_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC7_2(:,:,:)
  Real(Kind=jprb), pointer :: FAC1_3(:,:)
  Real(Kind=jprb), pointer :: FAC2_3(:,:)
  Real(Kind=jprb), pointer :: FAC3_3(:,:)
  Real(Kind=jprb), pointer :: FAC4_3(:,:)
  Real(Kind=jprb), pointer :: FAC5_3(:,:)
  Real(Kind=jprb), pointer :: FAC6_3(:,:)
  Real(Kind=jprb), pointer :: FAC7_3(:,:)
  Real(Kind=jprb), pointer :: cloudy(:,:)
  end type radiance_aux

  type raytracing_type
  Real(Kind=jprb), pointer :: LTICK   (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: HGPL    (:,:)   !(levels+1)
  Real(Kind=jprb), pointer :: DAIR    (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: DMAIR   (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: R       (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: RATOESUN(:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: RATOESAT(:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: ZASUN   (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: ZASAT   (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: INT     (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: HL      (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: PPW     (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: DISPCO2 (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: PATHSAT (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: PATHSUN (:,:)   !(levels,profiles)
  Real(Kind=jprb), pointer :: PATHEFF (:,:)
  end type raytracing_type

  type sunglint_type_s
  Real(Kind=jprb)          :: CSI
  Real(Kind=jprb)          :: ALFA
  Real(Kind=jprb)          :: C_SHAD
  Real(Kind=jprb)          :: P_PRIME
  Real(Kind=jprb)          :: PXY_GAMMAXY
  Real(Kind=jprb)          :: GAMMA_O
  Real(Kind=jprb)          :: GAMMA_P
  Real(Kind=jprb)          :: G_SHAD
  Real(Kind=jprb)          :: GAMMAX
  Real(Kind=jprb)          :: Q_SHAD
  Real(Kind=jprb)          :: ZENSAT
  Real(Kind=jprb)          :: ZENSUN
  Real(Kind=jprb)          :: DAZNG
  Real(Kind=jprb)          :: FAC1
  Real(Kind=jprb)          :: A_SHAD
  Real(Kind=jprb)          :: B_SHAD
  Real(Kind=jprb)          :: LAMBDA_A
  Real(Kind=jprb)          :: LAMBDA_B
  Real(Kind=jprb)          :: X_U
  Real(Kind=jprb)          :: ALFA1
  Real(Kind=jprb)          :: OMEGA_M
  Real(Kind=jprb)          :: WINDSP
  Real(Kind=jprb)          :: WANGL
  Real(Kind=jprb)          :: GAMMA_SQ
  Real(Kind=jprb)          :: GLINT
  Real(Kind=jprb)          :: OMEGA
  end type sunglint_type_s

  type sunglint_type
  type(sunglint_type_s), pointer :: s(:)
  Real(Kind=jprb),pointer  :: BETA (:,:)
  Real(Kind=jprb),pointer  :: PSI  (:,:)
  end type sunglint_type

  type transmission_scatt_ir_type
  Real(Kind=jprb)   , pointer  :: opdps         (:,:)
  Real(Kind=jprb)   , pointer  :: opdpa         (:,:)
  Real(Kind=jprb)   , pointer  :: gpar          (:,:)
  Real(Kind=jprb)   , pointer  :: gpartot       (:,:)
  Real(Kind=jprb)   , pointer  :: opdpscls      (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpacls      (:,:,:)
  Real(Kind=jprb)   , pointer  :: gparcls       (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpaerla     (:,:)
  Real(Kind=jprb)   , pointer  :: opdpcldla     (:,:)
  Real(Kind=jprb)   , pointer  :: opdpsaer      (:,:)
  Real(Kind=jprb)   , pointer  :: opdpaaer      (:,:)
  Real(Kind=jprb)   , pointer  :: gparaera      (:,:)
  Real(Kind=jprb)   , pointer  :: gparaer       (:,:)
  Real(Kind=jprb)   , pointer  :: azphup        (:,:)
  Real(Kind=jprb)   , pointer  :: azphdo        (:,:)
  Real(Kind=jprb)   , pointer  :: azphupcls     (:,:,:)
  Real(Kind=jprb)   , pointer  :: azphdocls     (:,:,:)
  Real(Kind=jprb)   , pointer  :: azphuptot     (:,:)
  Real(Kind=jprb)   , pointer  :: azphdotot     (:,:)
  Real(Kind=jprb)   , pointer  :: azphaerup     (:,:)
  Real(Kind=jprb)   , pointer  :: azphaerdo     (:,:)
  Real(Kind=jprb)   , pointer  :: azphaerupa    (:,:)
  Real(Kind=jprb)   , pointer  :: azphaerdoa    (:,:)
  Real(Kind=jprb)   , pointer  :: phasintupref  (:,:,:)
  Real(Kind=jprb)   , pointer  :: phasintdoref  (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpabs       (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpsca       (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpac        (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpacl       (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpacsun     (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpaclsun    (:,:,:)
  Real(Kind=jprb)   , pointer  :: azphacup      (:,:,:)
  Real(Kind=jprb)   , pointer  :: azphacdo      (:,:,:)
  Real(Kind=jprb)   , pointer  :: bcksp         (:,:,:)
  Real(Kind=jprb)   , pointer  :: opdpext       (:,:,:)
  Real(Kind=jprb)   , pointer  :: ssa           (:,:,:)
  end type transmission_scatt_ir_type


  type rttov_pccomp
  Real(Kind=jprb), pointer  :: pcscores(:)    ! Principal component scores
  Real(Kind=jprb), pointer  :: bt_pccomp(:)   ! Brightness temp equivalent to radiances
                                              ! reconstructed using principal components
  Real(Kind=jprb), pointer  :: clear_pccomp(:)! Clear radiances reconstructed using principal
                                              ! components
  end type rttov_pccomp

  type rttov_coef_pccomp1
  Integer(Kind=jpim)           :: fmv_pc_npred          ! Number of predictors in the regression set
  Integer(Kind=jpim), pointer  :: predictindex  (:)     ! Precitors channel indices
  Real(Kind=jprb)   , pointer  :: coefficients  (:,:)   ! Regression coefficients
  end type rttov_coef_pccomp1


  type rttov_coef_pccomp
    Integer(Kind=jpim)           :: fmv_pc_comp_pc
    Integer(Kind=jpim)           :: fmv_pc_sets           ! Number of regression sets
    Integer(Kind=jpim)           :: fmv_pc_mnum           ! Maximum number of eigenvectors
    Integer(Kind=jpim)           :: fmv_pc_nchn           ! Number of channels for which eigenvectors are available
    Integer(Kind=jpim)           :: fmv_pc_nchn_noise     ! Number of channels for which instrument noise is available
    Integer(Kind=jpim)           :: fmv_pc_nche           ! Number of channels for which emissisity coefs are available
    Integer(Kind=jpim)           :: fmv_pc_gas            ! Number of gases for which a reference profile is given
    Real(Kind=jprb)   , pointer  :: eigenvectors  (:,:)   ! Eigenvectors
    Integer(Kind=jpim), pointer  :: emiss_chn     (:)     ! Number of channels for which emissivity coefficients are
    Real   (Kind=jprb), pointer  :: emiss_c1      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c2      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c3      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c4      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c5      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c6      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c7      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c8      (:)     ! Emissivity coefficient
    Real   (Kind=jprb), pointer  :: emiss_c9      (:)     ! Emissivity coefficient
    Integer(Kind=jpim)           :: fmv_pc_nlev           ! Number of reference profile levels
    Real(Kind=jprb), Pointer     :: ref_pc_prfl_p (:)     ! pressure  (hPa)       (levels)
    Real(Kind=jprb), Pointer     :: ref_pc_prfl_mr(:,:)   ! mixing ratio (ppmv)   (levels)
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_tmin(:)   ! Profile limit :temperature
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_tmax(:)   ! Profile limit :temperature
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_qmin(:)   ! Profile limit :water vapour
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_qmax(:)   ! Profile limit :water vapour
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_ozmin(:)  ! Profile limit :ozone
    Real(Kind=jprb), Pointer     :: lim_pc_prfl_ozmax(:)  ! Profile limit :ozone
    Real(Kind=jprb)              :: lim_pc_prfl_pmin      ! Surface pressure
    Real(Kind=jprb)              :: lim_pc_prfl_pmax      ! Surface pressure
    Real(Kind=jprb)              :: lim_pc_prfl_tsmin     ! Surface temperature
    Real(Kind=jprb)              :: lim_pc_prfl_tsmax     ! Surface temperature
    Real(Kind=jprb)              :: lim_pc_prfl_skmin     ! Skin temperature
    Real(Kind=jprb)              :: lim_pc_prfl_skmax     ! Skin temperature
    Real(Kind=jprb)              :: lim_pc_prfl_wsmin     ! 10m wind speed
    Real(Kind=jprb)              :: lim_pc_prfl_wsmax     ! 10m wind speed
    Real(Kind=jprb), Pointer     :: co2_pc_ref    (:)     ! Fixed co2 profile to be used in the computation of PC's
    Real(Kind=jprb), Pointer     :: n2o_pc_ref    (:)     ! Fixed n2o profile to be used in the computation of PC's
    Real(Kind=jprb), Pointer     :: co_pc_ref     (:)     ! Fixed co  profile to be used in the computation of PC's
    Real(Kind=jprb), Pointer     :: ch4_pc_ref    (:)     ! Fixed ch4 profile to be used in the computation of PC's
    Real(Kind=jprb), Pointer     :: noise_in      (:)     ! Noise values for the channels whose radiances are
                                                          ! reconstrucetd using principal components
    Real(Kind=jprb), Pointer     :: noise         (:)     ! Noise values for the channels whose radiances are
                                                          ! used as predictors in the computation of principal components
    Integer(Kind=jpim), Pointer  :: ff_ori_chn_in(:)
    Real(Kind=jprb),    Pointer  :: ff_cwn_in(:)          ! central wave number of reconstructed radiances
    Real(Kind=jprb),    Pointer  :: ff_bco_in (:)         ! band correction offset (K)
    Real(Kind=jprb),    Pointer  :: ff_bcs_in (:)         ! band correction slope (K/K)
    Real(Kind=jprb),    pointer  :: planck1_in(:)         ! C1 * Nu**3
    Real(Kind=jprb),    pointer  :: planck2_in(:)         ! C2 * Nu
    type(rttov_coef_pccomp1), pointer:: pcreg         (:)
  end type rttov_coef_pccomp

  type rttov_coef_scatt_ir
  Integer(Kind=jpim)           :: dim
  Integer(Kind=jpim)           :: fmv_aer_chn    ! number of channels for which aerosol optical parameters are stored
  Integer(Kind=jpim)           :: fmv_wcl_chn
  Integer(Kind=jpim)           :: fmv_icl_chn
  Integer(Kind=jpim)           :: fmv_aer_pha_chn
  Integer(Kind=jpim)           :: fmv_wcl_pha_chn
  Integer(Kind=jpim)           :: fmv_icl_pha_chn
  Integer(Kind=jpim)           :: fmv_aer_sun_chn
  Integer(Kind=jpim)           :: fmv_wcl_sun_chn
  Integer(Kind=jpim)           :: fmv_icl_sun_chn
  Integer(Kind=jpim)           :: fmv_aer_comp
  Integer(Kind=jpim)           :: fmv_wcl_comp
  Integer(Kind=jpim)           :: fmv_icl_comp
  Integer(Kind=jpim)           :: fmv_icl_ishp
  Integer(Kind=jpim)           :: fmv_aer_pha_ioff
  Integer(Kind=jpim)           :: fmv_wcl_pha_ioff
  Integer(Kind=jpim)           :: fmv_icl_pha_ioff
  Integer(Kind=jpim)           :: fmv_aer_ph
  Integer(Kind=jpim)           :: fmv_wcl_ph
  Integer(Kind=jpim)           :: fmv_icl_ph
  Integer(Kind=jpim)           :: icl_nabs
  Integer(Kind=jpim)           :: icl_nsca
  Integer(Kind=jpim)           :: icl_nbpr
  Integer(Kind=jpim), Pointer  :: fmv_aer_rh    (:)
  Integer(Kind=jpim), Pointer  :: fmv_wcl_rh    (:)
  Real   (Kind=jprb), Pointer  :: fmv_aer_rh_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_wcl_rh_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_wcl_ph_val_cos(:)
  Real   (Kind=jprb)           :: fmv_wcl_ph_val_min
  Integer(Kind=jpim), Pointer  :: ifmv_wcl_ph_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_aer_ph_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_aer_ph_val_cos(:)
  Integer(Kind=jpim), Pointer  :: ifmv_aer_ph_val(:)
  Real   (Kind=jprb)           :: fmv_aer_ph_val_min
  Real   (Kind=jprb), Pointer  :: fmv_wcl_ph_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_icl_ph_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_icl_ph_val_cos(:)
  Real   (Kind=jprb)           :: fmv_icl_ph_val_min
  Integer(Kind=jpim), Pointer  :: ifmv_icl_ph_val(:)
  Real   (Kind=jprb), Pointer  :: fmv_icl_dg    (:,:)
  Integer(Kind=jpim), Pointer  :: channels_solar(:)
  Real(Kind=jprb)   , pointer  :: abs           (:,:)
  Real(Kind=jprb)   , pointer  :: sca           (:,:)
  Real(Kind=jprb)   , pointer  :: bpr           (:,:)
  Real(Kind=jprb)   , pointer  :: pha           (:,:,:)
  Real(Kind=jprb)   , pointer  :: confac        (:)
  end type rttov_coef_scatt_ir

  type rttov_optpar_ir
    type(rttov_coef_scatt_ir), pointer :: optpaer(:)
    type(rttov_coef_scatt_ir), pointer :: optpwcl(:)
    type(rttov_coef_scatt_ir), pointer :: optpicl(:)
  end type rttov_optpar_ir

  type ircld_type
  Integer(Kind=jpim), pointer  :: nstream(:)
  Integer(Kind=jpim), pointer  :: nstreamref(:)
  Integer(Kind=jpim), pointer  :: iloop(:)
  Integer(Kind=jpim), pointer  :: icount(:)
  Integer(Kind=jpim), pointer  :: icounstr(:)
  Integer(Kind=jpim), pointer  :: icount1(:)
  Real(Kind=jprb)   , pointer  :: xstrclr(:)
  Integer(Kind=jpim), pointer  :: icldarr   (:,:,:)
  Real(Kind=jprb)   , pointer  :: xstrref1  (:,:,:)
  Real(Kind=jprb)   , pointer  :: xstrref2  (:,:,:)
  Integer(Kind=jpim), pointer  :: cldtyp    (:,:,:)
  Integer(Kind=jpim), pointer  :: indexstr  (:,:)
  Integer(Kind=jpim), pointer  :: icount1ref(:,:)
  Integer(Kind=jpim), pointer  :: iloopin   (:,:)
  Integer(Kind=jpim), pointer  :: iflag     (:,:)
  Real(Kind=jprb)   , pointer  :: xstr      (:,:)
  Real(Kind=jprb)   , pointer  :: xstrminref(:,:)
  Real(Kind=jprb)   , pointer  :: xstrref   (:,:)
  Real(Kind=jprb)   , pointer  :: cldcfr    (:,:)
  Real(Kind=jprb)   , pointer  :: maxcov    (:,:)
  Real(Kind=jprb)   , pointer  :: xstrmax   (:,:)
  Real(Kind=jprb)   , pointer  :: xstrmin   (:,:)
  Real(Kind=jprb)   , pointer  :: a         (:,:)
  Real(Kind=jprb)   , pointer  :: ntotref   (:,:)

  Real(Kind=jprb)   , pointer  :: tave      (:,:)
  Real(Kind=jprb)   , pointer  :: wmixave   (:,:)
  Real(Kind=jprb)   , pointer  :: xpresave  (:,:)
  Real(Kind=jprb)   , pointer  :: ppv       (:,:)
  Real(Kind=jprb)   , pointer  :: esw       (:,:)
  Real(Kind=jprb)   , pointer  :: esi       (:,:)
  Logical(Kind=jplm), pointer  :: flag      (:,:)
  end type ircld_type


  Type blob_type

    ! l = level
    ! a = aerosol
    ! p = profile
    ! c = cloud
    ! t = cloudtype
    ! i = ircld
    ! v = var
    ! x = ??
    ! s = stream
    ! h = channel
    ! r = raytracing
    ! o = coef
    ! e = predictor
    ! lp1 = l + 1
    !

    Real(Kind=jprb), pointer :: rlvp(:,:,:)
    Real(Kind=jprb), Pointer :: ralp(:,:,:)
    Real(Kind=jprb), Pointer :: rclp(:,:,:)
    Real(Kind=jprb), Pointer :: rtlp(:,:,:)

  End Type blob_type


  Type rttov_coefs
    Type(rttov_coef)           :: coef
    Type(rttov_coef_scatt_ir)  :: coef_scatt_ir
    Type(rttov_optpar_ir)      :: optp
    type(rttov_coef_pccomp)    :: coef_pccomp
  End Type

  Type rttov_options
    Integer(Kind=jpim) :: ipcreg = -1_jpim
    Real(Kind=jprb)    :: cldstr_threshold=-1.0_jprb ! Recommended to set this negative
    Integer(Kind=jpim) :: fastem_version = 0_jpim   ! Valid range: 1-5. Otherwise version taken from coef file.
    Logical(Kind=jplm) :: addinterp  = .false.
    Logical(Kind=jplm) :: addpc      = .false.
    Logical(Kind=jplm) :: addradrec  = .false.
    Logical(Kind=jplm) :: addsolar   = .false.
    Logical(Kind=jplm) :: addaerosl  = .false.
    Logical(Kind=jplm) :: addclouds  = .false.
    Logical(Kind=jplm) :: switchrad  = .false.
    Logical(Kind=jplm) :: spacetop   = .true.  ! treat user's model-top as space boundary
    Logical(Kind=jplm) :: lgradp     = .false.
    Logical(Kind=jplm) :: use_q2m    = .false. ! set to true to activate use of surface humidity
    Logical(Kind=jplm) :: apply_reg_limits = .false. ! set to true makes rttov_checkinput reset the profiles
                                                     ! variables to the regression limits if outside
    Logical(Kind=jplm) :: verbose_checkinput_warnings    = .true.
    Logical(Kind=jplm) :: ozone_data = .false.
    Logical(Kind=jplm) :: co2_data   = .false.
    Logical(Kind=jplm) :: n2o_data   = .false.
    Logical(Kind=jplm) :: co_data    = .false.
    Logical(Kind=jplm) :: ch4_data   = .false.
    Logical(Kind=jplm) :: clw_data   = .false.
    logical(Kind=jplm) :: addrefrac  = .false.
    Logical(Kind=jplm) :: do_checkinput = .true.
  End Type


  Type rttov_traj
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! and their dimensions are known before running RTTOV (nlevels, nprofiles, nchannels)
! it is possible to allocate these variables from outside RTTOV
!
    Type(profile_Type), Pointer :: profiles_COEF(:)
    Type(blob_type)             :: profiles_COEF_blob
    Type(predictors_Type)       :: predictors
    Type(raytracing_type)       :: raytracing
    Type(raytracing_type)       :: raytracing_COEF
    Type(ircld_type)            :: ircld
    Type(opdp_path_type)        :: opdp_path
    Type(opdp_path_type)        :: opdp_path_COEF

    Real(Kind=jprb), Pointer :: reflectivity(:) ! channels
    Real(Kind=jprb), Pointer :: fresnrefl(:)    ! channels
    Type(profile_aux)  :: aux_prof
    Type(profile_aux)  :: aux_prof_COEF
    Type(transmission_scatt_ir_type)  :: transmission_scatt_ir
    Type(sunglint_type):: sunglint

    Type(rttov_coefs),         Pointer :: coefs
    Integer(Kind=jpim)   :: nchannels
    Integer(Kind=jpim)   :: nlevels
    Integer(Kind=jpim)   :: nlayers
    Type(rttov_options)  :: opts
  End Type

  Type rttov_traj_dyn
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! but their dimensions are known when RTTOV starts running (nstreams)
!
    INTEGER(KIND=jpim)               :: nstreams
    TYPE(radiance_aux              ) :: auxrad_stream
    TYPE(transmission_scatt_ir_type) :: transmission_scatt_ir_stream
    TYPE(transmission_type_aux     ) :: transmission_aux
  End Type

  Type rttov_traj_sta
!
! Hold RTTOV trajectory; these variables do not have counterparts in TL, AD, K
!
    LOGICAL(KIND=jplm),  POINTER :: sun(:)              ! solar switch located in channel sequence
    REAL(KIND=jprb),     POINTER :: tau_ref          (:,:)
    REAL(KIND=jprb),     POINTER :: tau_ref_surf     (:)
    REAL(KIND=jprb),     POINTER :: tau_surf         (:)
    REAL(KIND=jprb),     POINTER :: tausun_ref       (:,:)
    REAL(KIND=jprb),     POINTER :: tausun_ref_surf  (:)
    REAL(KIND=jprb),     POINTER :: tausun_level     (:,:)
    REAL(KIND=jprb),     POINTER :: tausun_surf      (:)
    REAL(KIND=jprb),     POINTER :: opdp_ref_COEF    (:,:)! layer optical depth before threshold
    REAL(KIND=jprb),     POINTER :: od_level         (:,:)! sat to level optical depth
    REAL(KIND=jprb),     POINTER :: tau_level        (:,:)! sat to level transmittance
    REAL(KIND=jprb),     POINTER :: opdpsun_ref_COEF (:,:)! layer optical depth
    REAL(KIND=jprb),     POINTER :: odsun_level      (:,:)! sat to level optical depth
    REAL(KIND=jprb),     POINTER :: odsun_singlelayer(:,:)! single layer optical depth
    REAL(KIND=jprb),     POINTER :: od_frac          (:)
    TYPE(geometry_Type), POINTER :: angles           (:)! geometry angles
    TYPE(geometry_Type), POINTER :: angles_COEF      (:)! geometry angles
    TYPE(profile_Type),  POINTER :: profiles_COEF_ref(:)
    TYPE(blob_type)              :: profiles_COEF_blob_ref
    TYPE(radiance_aux)           :: auxrad
    TYPE(rttov_chanprof), POINTER :: chanprof_in(:)
    TYPE(rttov_chanprof), POINTER :: chanprof_pc(:)
  End Type

!
  Type rttov_chanprof
    Integer(Kind=jpim) :: chan
    Integer(Kind=jpim) :: prof
  End Type

  Type rttov_lbl_check
    Real(Kind=jprb), Pointer :: atm_layer(:,:)
    Logical(Kind=jplm) :: plane_geometry
  End Type


End Module rttov_types
