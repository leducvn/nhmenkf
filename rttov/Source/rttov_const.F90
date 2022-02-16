!
Module rttov_const
  ! Description:
  ! Definition of all parameters (constants) for RTTOV
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
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0   01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1   29/01/2003  New platforms and instruments (P Brunel)
  !                    Hard limits for input profiles
  !  1.2   19/02/2003  Some changes to limits and comments (R Saunders)
  !  1.3   06/05/2003  Change version number to 7.3.1
  !                    and add references for physical constants (P Brunel)
  !  1.4      08/2003  Added variables for MW scattering (F Chevallier)
  !  1.5   18/09/2003  Added coefficients for cloud absorption properties (P Francis)
  !  1.6   15/10/2003  Added new sections in parameter files for scatt   (F Chevallier)
  !  1.7   23/11/2003  Added new definitions of polarisations 2.1 (S English)
  !  1.8   25/08/2005  Made inst_name a parameter (R Saunders)
  !  1.9   11/01/2006  Added logical flag for surface humidity use (R Saunders)
  !  1.10  12/01/2006  Marco Matricardi (ECMWF):
  !           --       Added variables for CO2,CO,N2O and CH4 molecules.
  !           --       Added parameters for the computation of the refractive index
  !           --       of air.
  !  1.11  06/02/2006  Added logical flag for linear in tau approx (R Saunders)
  !  1.12  06/04/2006  Added Meghatropiques (R. Saunders)
  !  1.13  14/03/2007  Added units conversion constants
  !  1.14  16/05/2007  Added polarimetric sensor type (R Saunders)
  !  1.15  25/09/2007  Added maximum number of warnings for checkinput (P Brunel)
  !  1.16  11/10/2007  Remove zhusta* and zice* constants ( P.Marguinaud )
  !  1.17  07/12/2007  Remove maximum number of warnings for checkinput (P Brunel)
  !  1.18  12/12/2007  Added hard limits for trace gases (R Saunders)
  !  1.19  13/12/2007  Renamed linear_tau (R Saunders)
  !  1.20  01/11/2007  Added parameters for section length and AD/K code (A. Geer)
  !  1.21  16/01/2008  Facility to apply regression limits  (N. Bormann)
  !  1.22  04/03/2008  Made min hard limit > zero (R Saunders)
  !  1.23  14/04/2008  Added SSM/T2 (R Saunders)
  !  1.24  02/06/2008  Changed mixing ratio for CO (R Saunders)
  !  1.25  12/08/2008  Added SSMISZ for SSMIS chan19-22 - Zeeman (P. Rayer)
  !  1.26  29/01/2009  Add Kalpana and FY-3 (R Saunders)
  !  1.27  26/05/2009  Add more platforms and sensors (R Saunders)
  !  1.28  02/12/2009  Add principal component capability (Marco matricardi)
  !  1.29  15/01/2010  Add rttov9 intervals constants (P Marguinaud)
  !  1.30  05/07/2010  Add maximum solar zenith angle constant (J Hocking)
  !  1.31  01/02/2011  Updates to platform and sensor lists (J Hocking)
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !

  Use parkind1, Only : jpim     ,jprb
  Implicit None

  !1.0 Precision and numerical constants
  ! Try to ensure this is large enough to avoid overflows in reciprocals but small enough to not affect the calculations.
  ! these parameters are defined at bottom of module (because they use 
  Real(jprb), parameter :: max_exp_exponent = 50_jprb ! approx 1e22 which should be sufficiently big for most purposes
  Real(jprb), parameter :: min_exponent     = 1e-16_jprb ! approx log_10(1+2^-52) - anything raised to this power or smaller 
                                                         ! should be approx equal to 1
! small_val is defined in rttov_transmit to avoid compiler incompatibility
!   ! small_val is used in rttov_transmit to ensure small values do not result in underflows. In subsequent calculations
!   ! these values are multiplied together hence the exponent of 1/3.
!   Real(jprb)            :: small_val = (tiny(min_exponent)) ** (0.333333_jprb) ! XLF doesn't like 1/3

  !1.1 general
  !-----------
  ! Version number of the current code

  Integer(Kind=jpim), Parameter :: version = 10
  Integer(Kind=jpim), Parameter :: release = 2
  Integer(Kind=jpim), Parameter :: minor_version = 0

  Integer(Kind=jpim), Parameter :: version_compatible_min = 10 ! minimum version number
  Integer(Kind=jpim), Parameter :: version_compatible_max = 10 ! maximum version number
          ! compatible for coefficients.
          ! coef files with "id_comp_lvl" outside range will be rejected

  Character (len=16), Parameter :: rttov_magic_string = '%RTTOV_COEFF    '
  Real(Kind=jprb),    Parameter :: rttov_magic_number = 1.2345E+12_JPRB

  Integer(Kind=jpim), Parameter :: default_err_unit = 0  ! standard error unit number
                              ! standard error unit number is 7 for HPUX

  !1.2 physical constants
  !----------------------
  ! Molecular weights  (g/mole) are calculated by adding NIST Standard Atomic Weights
  ! Molecular weight of dry air refers to US standard atmosphere 1976
  ! NIST  Standard Atomic Weight are:
  ! H    1.00794   (7)
  ! C   12.0107    (8)
  ! N   14.0067    (2)
  ! O   15.9994    (3)
  Real(Kind=jprb), Parameter :: mair = 28.9644_JPRB
  Real(Kind=jprb), Parameter :: mh2o = 18.01528_JPRB
  Real(Kind=jprb), Parameter :: mo3  = 47.9982_JPRB
  Real(Kind=jprb), Parameter :: mco2 = 44.0095_JPRB
  Real(Kind=jprb), Parameter :: mch4 = 16.04246_JPRB
  Real(Kind=jprb), Parameter :: mn2o = 44.0128_JPRB
  Real(Kind=jprb), Parameter :: mco  = 28.0101_JPRB

  ! Gravity from NIST 9.80665 ms-1 (exact)
  Real(Kind=jprb), Parameter :: gravity = 9.80665_JPRB

  !
  ! Kaye & Laby latest library edition is 16e 1995, and gives
  ! * standard value  g = 9.80665 ms-1 exactly (p.191)
  ! * earth mean radius r= 6371.00 km (p191)
  !    [defined as [(r_equator)^2 (r_pole)]^1/3]
  Real(Kind=jprb), Parameter :: pi      = 3.1415926535_JPRB
  Real(Kind=jprb), Parameter :: deg2rad = pi/180.0_JPRB
  Real(Kind=jprb), Parameter :: earthradius = 6371.00_JPRB
  Real(Kind=jprb), Parameter :: flatt       = 3.3528107E-3_JPRB
  Real(Kind=jprb), Parameter :: omega       = 7292115E-11_JPRB
  Real(Kind=jprb), Parameter :: eqrad       = 6378.137_JPRB
  Real(Kind=jprb), Parameter :: grave       = 9.7803267715_JPRB
  Real(Kind=jprb), Parameter :: z4pi_r      = 0.0795774715_JPRB
  Real(Kind=jprb), Parameter :: pi_r        = 0.3183098862_JPRB

  ! The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set
  ! Fixsen D.J. et all
  ! Astrophysical Journal v.473, p.576 December 1996
  ! CMBR = 2.728 +- 0.004K
  Real(Kind=jprb), Parameter :: tcosmic     = 2.728_JPRB
  !  Real(Kind=jprb), Parameter :: tcosmic     = 0.1_JPRB !used for ECMWF tests

  ! Universal gas constant R = 8.314510 J/mol/K
  Real(Kind=jprb), Parameter :: rgp = 8.314510_JPRB
  Real(Kind=jprb), Parameter :: rgc = 8.314472_JPRB

  ! mean molar mass of dry air rm = 0.0289644 kg.mol^-1
  Real(Kind=jprb), Parameter :: rm = 0.0289644_JPRB

  ! units conversion from  mixing ratio to ppmv
  Real(Kind=jprb), Parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB
  Real(Kind=jprb), Parameter :: o3_mixratio_to_ppmv = 6.03504e+5_JPRB
  Real(Kind=jprb), Parameter :: co2_mixratio_to_ppmv= 6.58114e+5_JPRB
  Real(Kind=jprb), Parameter :: co_mixratio_to_ppmv = 1.0340699e+6_JPRB
  Real(Kind=jprb), Parameter :: n2o_mixratio_to_ppmv= 6.58090e+5_JPRB
  Real(Kind=jprb), Parameter :: ch4_mixratio_to_ppmv= 1.80548e+6_JPRB

  ! zero temperature(K)
  Real(Kind=jprb), Parameter :: t0 =273.15

  !1.3 satellite and instrument information
  !----------------------------------------

  !platform id codes
  Integer(Kind=jpim), Parameter :: nplatforms = 37
  Integer(Kind=jpim), Parameter :: &
       & platform_id_noaa      = 1, &
       & platform_id_dmsp      = 2, &
       & platform_id_meteosat  = 3, &
       & platform_id_goes      = 4, &
       & platform_id_gms       = 5, &
       & platform_id_fy2       = 6, &
       & platform_id_trmm      = 7, &
       & platform_id_ers       = 8, &
       & platform_id_eos       = 9, &
       & platform_id_metop     = 10, &
       & platform_id_envisat   = 11, &
       & platform_id_msg       = 12, &
       & platform_id_fy1       = 13, &
       & platform_id_adeos     = 14, &
       & platform_id_mtsat     = 15, &
       & platform_id_coriolis  = 16, &
       & platform_id_jpss      = 17, &
       & platform_id_gifts     = 18, &
       & platform_id_sentinel3 = 19, &
       & platform_id_meghatr   = 20, &
       & platform_id_kalpana   = 21, &
       & platform_id_insat_3d  = 22, &
       & platform_id_fy3       = 23, &
       & platform_id_coms      = 24, &
       & platform_id_meteorm   = 25, &
       & platform_id_gosat     = 26, &
       & platform_id_calipso   = 27, &
       & platform_id_dummy     = 28, &
       & platform_id_gcomw     = 29, &
       & platform_id_nimbus    = 30, &
       & platform_id_himawari  = 31, &
       & platform_id_mtg       = 32, &
       & platform_id_saral     = 33, &
       & platform_id_metopsg   = 34, &
       & platform_id_landsat   = 35, &
       & platform_id_jason     = 36, &
       & platform_id_gpm       = 37

  !platform names
  Character (len=9), Parameter :: platform_name(nplatforms) = &
       & (/ 'noaa     ', 'dmsp     ', 'meteosat ', 'goes     ', 'gms      ', &
          & 'fy2      ', 'trmm     ', 'ers      ', 'eos      ', 'metop    ', &
          & 'envisat  ', 'msg      ', 'fy1      ', 'adeos    ', 'mtsat    ', &
          & 'coriolis ', 'jpss     ', 'gifts    ', 'sentinel3', 'meghatr  ', &
          & 'kalpana  ', 'insat_3d ', 'fy3      ', 'coms     ', 'meteor-m ', &
          & 'gosat    ', 'calipso  ', 'dummy    ', 'gcom-w   ', 'nimbus   ', &
          & 'himawari ', 'mtg      ', 'saral    ', 'metopsg  ', 'landsat  ', &
          & 'jason    ', 'gpm      '  /)

  !instrument id codes
  Integer(Kind=jpim), Parameter :: &
       & inst_id_hirs   =  0, inst_id_msu    =  1, inst_id_ssu    =  2, inst_id_amsua   =  3, &
       & inst_id_amsub  =  4, inst_id_avhrr  =  5, inst_id_ssmi   =  6, inst_id_vtpr1   =  7, &
       & inst_id_vtpr2  =  8, inst_id_tmi    =  9, inst_id_ssmis  = 10, inst_id_airs    = 11, &
       & inst_id_hsb    = 12, inst_id_modis  = 13, inst_id_atsr   = 14, inst_id_mhs     = 15, &
       & inst_id_iasi   = 16, inst_id_amsr   = 17, inst_id_mtsatim= 18, inst_id_atms    = 19, &
       & inst_id_mviri  = 20, inst_id_seviri = 21, inst_id_goesim = 22, inst_id_goessd  = 23, &
       & inst_id_gmsim  = 24, inst_id_vissr  = 25, inst_id_mvisr  = 26, inst_id_cris    = 27, &
       & inst_id_cmis   = 28, inst_id_viirs  = 29, inst_id_windsat= 30, inst_id_gifts   = 31, &
       & inst_id_ssmt1  = 32, inst_id_ssmt2  = 33, inst_id_saphir = 34, inst_id_madras  = 35, &
       & inst_id_ssmisz = 36, inst_id_kavhrr = 37, inst_id_iimager= 38, inst_id_isoundr = 39, &
       & inst_id_mwts   = 40, inst_id_mwhs   = 41, inst_id_iras   = 42, inst_id_mwri    = 43, &
       & inst_id_abi    = 44, inst_id_mi     = 45, inst_id_msumr  = 46, inst_id_tansofts= 47, &
       & inst_id_iir    = 48, inst_id_mwr    = 49, inst_id_dummyir= 50, inst_id_dummymw = 51, &
       & inst_id_dummyhi= 52, inst_id_dummypo= 53, inst_id_scams  = 54, inst_id_smmr    = 55, &
       & inst_id_ahi    = 56, inst_id_irs    = 57, inst_id_altika = 58, inst_id_iasing  = 59, &
       & inst_id_tm     = 60, inst_id_fci    = 61, inst_id_amsr1  = 62, inst_id_amsr2   = 63, &
       & inst_id_vissr2 = 64, inst_id_slstr  = 65, inst_id_tirs   = 66, inst_id_amr     = 67, &
       & inst_id_oli    = 68, inst_id_iris   = 69, inst_id_ici    = 70, inst_id_gmi     = 71, &
       & inst_id_mwts2  = 72, inst_id_mwhs2  = 73

  Integer(Kind=jpim), Parameter :: ninst = 74
  ! List of instruments  !!!! HIRS is number 0
  Character (len=8), Dimension(0:ninst-1),parameter :: inst_name =       &
        & (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ',  &
           & 'avhrr   ', 'ssmi    ', 'vtpr1   ', 'vtpr2   ', 'tmi     ',  &
           & 'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ',  &
           & 'mhs     ', 'iasi    ', 'amsr    ', 'imager  ', 'atms    ',  &
           & 'mviri   ', 'seviri  ', 'imager  ', 'sounder ', 'imager  ',  &
           & 'vissr   ', 'mvisr   ', 'cris    ', 'cmis    ', 'viirs   ',  &
           & 'windsat ', 'gifts   ', 'ssmt1   ', 'ssmt2   ', 'saphir  ',  &
           & 'madras  ', 'ssmisz  ', 'kavhrr  ', 'iimager ', 'isoundr ',  &
           & 'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ',  &
           & 'mi      ', 'msumr   ', 'tansofts', 'iir     ', 'mwr     ',  &
           & 'dummyir ', 'dummymw ', 'dummyhi ', 'dummypo ', 'scams   ',  &
           & 'smmr    ', 'ahi     ', 'irs     ', 'altika  ', 'iasing  ',  &
           & 'tm      ', 'fci     ', 'amsr    ', 'amsr2   ', 'vissr   ',  &
           & 'slstr   ', 'tirs    ', 'amr     ', 'oli     ', 'iris    ',  &
           & 'ici     ', 'gmi     ', 'mwts2   ', 'mwhs2   ' /)


  !1.4 Coefficient file Section names
  !----------------------------------
  Integer(Kind=jpim), Parameter :: nsections = 38
  Integer(Kind=jpim), Parameter :: lensection = 23
  Character(len=lensection), Parameter :: section_types(nsections) = &
    & (/ 'IDENTIFICATION         ', 'LINE-BY-LINE           ', &
       & 'FAST_MODEL_VARIABLES   ', 'FILTER_FUNCTIONS       ', &
       & 'FUNDAMENTAL_CONSTANTS  ', 'SSIREM                 ', &
       & 'FASTEM                 ', 'REFERENCE_PROFILE      ', &
       & 'PROFILE_LIMITS         ', 'FAST_COEFFICIENTS      ', &
       & 'COEF_SUB_FILES         ', 'GAZ_UNITS              ', &
       & 'DIMENSIONS             ', 'FREQUENCIES            ', &
       & 'HYDROMETEOR            ', 'CONVERSIONS            ', &
       & 'EXTINCTION             ', 'ALBEDO                 ', &
       & 'ASYMMETRY              ', 'GAS_SPECTRAL_INTERVAL  ', &
       & 'TRANSMITTANCE_TRESHOLD ', 'SOLAR_SPECTRUM         ', &
       & 'WATER_OPTICAL_CONSTANT ', 'WAVE_SPECTRUM          ', &
       & 'AEROSOLS_PARAMETERS    ', 'AEROSOLS_COMPONENTS    ', &
       & 'WATERCLOUD_TYPES       ', 'WATERCLOUD_PARAMETERS  ', &
       & 'ICECLOUD_TYPES         ', 'HEXAGONAL_PARAMETERS   ', &
       & 'AGGREGATE_PARAMETERS   ', 'PRINCOMP_PREDICTORS    ', &
       & 'PRINCOMP_EIGENVECTORS  ', 'PRINCOMP_COEFFICIENTS  ', &
       & 'EMISSIVITY_COEFFICIENTS', 'PC_REFERENCE_PROFILE   ', &
       & 'PC_PROFILE_LIMITS      ', 'INSTRUMENT_NOISE       '/)

  !sensors id codes
  Integer(Kind=jpim), Parameter :: nsensors = 4
  Integer(Kind=jpim), Parameter :: &
       & sensor_id_ir     = 1, &
       & sensor_id_mw     = 2, &
       & sensor_id_hi     = 3, &
       & sensor_id_po     = 4

  !sensors names
  Character (len=2), Parameter :: sensor_name(nsensors) = &
       & (/ 'ir', 'mw', 'hi', 'po' /)

  ! these codes are for the instrument from the inst_name array
  Integer(Kind=jpim), Parameter :: sensor_id(0:ninst-1) = (/ &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_mw,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_ir,  &
    sensor_id_po, sensor_id_hi, sensor_id_mw, sensor_id_mw, sensor_id_mw,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_hi, sensor_id_po, sensor_id_mw,  &
    sensor_id_mw, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_hi,  &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_hi,  &
    sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_mw /)

  !gas id codes
  Integer(Kind=jpim), Parameter :: ngases_max = 8
  Integer(Kind=jpim), Parameter :: &
        & gas_id_mixed       = 1, &
        & gas_id_watervapour = 2, &
        & gas_id_ozone       = 3, &
        & gas_id_wvcont      = 4, &
        & gas_id_co2         = 5, &
        & gas_id_n2o         = 6, &
        & gas_id_co          = 7, &
        & gas_id_ch4         = 8

  !gas names
  Character (len=12), Parameter :: gas_name(ngases_max) = &
        & (/ 'Mixed_gases ', &
           & 'Water_vapour', &
           & 'Ozone       ', &
           & 'WV_Continuum', &
           & 'CO2         ', &
           & 'N2O         ', &
           & 'CO          ', &
           & 'CH4         ' /)

  !gas units
  Integer(Kind=jpim), Parameter :: ngases_unit = 2
  Integer(Kind=jpim), Parameter :: &
        & gas_unit_specconc  = 1, &
        & gas_unit_ppmv      = 2
  Character (len=12), Parameter :: gas_unit_name(ngases_unit) = &
        & (/ 'spec. concen', &
           & 'ppmv        '  /)


  !1.5 error reporting
  !-------------------
  !error status values
  Integer(Kind=jpim), Parameter :: nerrorstatus = 3
  Integer(Kind=jpim), Parameter :: errorstatus_success = 0
  Integer(Kind=jpim), Parameter :: errorstatus_warning = 1
  Integer(Kind=jpim), Parameter :: errorstatus_fatal   = 2
  Integer(Kind=jpim), Parameter :: errorstatus_info    = 3 
  Character(len=*), Parameter :: errorstatus_text(0:nerrorstatus) = & 
       & (/ 'success', & 
       & 'warning', & 
       & 'fatal  ', & 
       & 'info   '  /)


  !1.6 surface types
  !-----------------
  Integer(Kind=jpim), Parameter :: nsurftype = 2
  Integer(Kind=jpim), Parameter :: surftype_land = 0
  Integer(Kind=jpim), Parameter :: surftype_sea = 1
  Integer(Kind=jpim), Parameter :: surftype_seaice = 2

  !1.7 water types
  !---------------
  Integer(Kind=jpim), Parameter :: nwatertype = 1
  Integer(Kind=jpim), Parameter :: watertype_fresh_water = 0
  Integer(Kind=jpim), Parameter :: watertype_ocean_water = 1

  !1.8 cloud emissivity
  !---------------------
  Integer(Kind=jpim), Parameter :: overlap_scheme = 2    ! overlap scheme
  ! 1 => Geleyn and Hollingsworth (1979)
  ! 2 => Raisanen (1998)

  !
  !1.9 Hard limits for control of input profile
  !--------------------------------------------
  ! Temperature
  Real(Kind=jprb), Parameter :: tmax   = 400.0_JPRB       ! degK
  Real(Kind=jprb), Parameter :: tmin   = 90.0_JPRB        ! degK
  ! Water Vapour
  Real(Kind=jprb), Parameter :: qmax   = 0.60E+06_JPRB    ! ppmv 0.373_JPRB kg/kg
  Real(Kind=jprb), Parameter :: qmin   = 0.1E-10_JPRB     ! ppmv
  ! Ozone
  Real(Kind=jprb), Parameter :: o3max  = 1000.0_JPRB      ! ppmv  1.657E-3_JPRB kg/kg
  Real(Kind=jprb), Parameter :: o3min  = 0.1E-10_JPRB     ! ppmv
  ! CO2
  Real(Kind=jprb), Parameter :: co2max = 1000.0_JPRB      ! ppmv
  Real(Kind=jprb), Parameter :: co2min = 0.1E-10_JPRB     ! ppmv
  ! CO
  Real(Kind=jprb), Parameter :: comax  = 10.0_JPRB        ! ppmv
  Real(Kind=jprb), Parameter :: comin  = 0.1E-10_JPRB     ! ppmv
  ! N2O
  Real(Kind=jprb), Parameter :: n2omax = 10.0_JPRB        ! ppmv
  Real(Kind=jprb), Parameter :: n2omin = 0.1E-10_JPRB     ! ppmv
  ! CH4
  Real(Kind=jprb), Parameter :: ch4max = 50.0_JPRB        ! ppmv
  Real(Kind=jprb), Parameter :: ch4min = 0.1E-10_JPRB     ! ppmv
  ! Cloud Liquid Water
  Real(Kind=jprb), Parameter :: clwmax = 1.0_JPRB         ! kg/kg
  Real(Kind=jprb), Parameter :: clwmin = 0.0_JPRB         ! kg/kg
  ! Surface Pressure
  Real(Kind=jprb), Parameter :: pmax   = 1100.0_JPRB      ! surface pressure hPa
  Real(Kind=jprb), Parameter :: pmin   = 400.0_JPRB       ! hPa
  ! Surface Wind
  Real(Kind=jprb), Parameter :: wmax   =  100.0_JPRB      ! surface wind speed (m/s)
  ! Zenith Angle
  Real(Kind=jprb), Parameter :: zenmax = 75.0_JPRB        ! zenith angle (Deg) = secant 3.86_JPRB
  ! Cloud Top Pressure
  Real(Kind=jprb), Parameter :: ctpmax = 1100.0_JPRB      ! (hPa)
  Real(Kind=jprb), Parameter :: ctpmin =   50.0_JPRB      ! (hPa)
  ! Magnetic field strength
  Real(Kind=jprb), Parameter :: bemax = 0.7_JPRB          ! (Gauss)
  Real(Kind=jprb), Parameter :: bemin = 0.2_JPRB          ! (Guass)


  !1.10  Maximum Optical Depth
  !--------------------------
  ! maximum value of optical depth for transmittance calculation
  ! e(-30) -> 10**-14
  ! e(-50) -> 10**-22
  Real(Kind=jprb), Parameter  :: max_optical_depth = 50._JPRB

  !1.11  Maximum solar zenith angle for which to apply solar calculation
  !---------------------------------------------------------------------
  Real(Kind=jprb), Parameter  :: max_sol_zen = 84._JPRB
  
  
  !2 RTTOV7 aux parameters
  !-------------------------
  Integer(Kind=jpim), Parameter :: fastem_sp = 5  ! max. number of fastem surface parameters
  Real(Kind=jprb), Parameter    :: mwcldtp = 322.0_JPRB  ! Upper pressure level (HPa) for lwp calcs
  Real(Kind=jprb), Parameter    :: pressure_top = 0.004985_JPRB ! Pressure of top level for
                                                ! Line/Line calculations (hPa)
  Real(Kind=jprb) , Dimension(8), Parameter :: dcoeff =        &! Debye coefs
        & (/ 17.1252_JPRB, 134.2450_JPRB, 310.2125_JPRB,  5.667_JPRB,   &
          & 188.7979_JPRB,  80.5419_JPRB,   0.1157_JPRB,  4.8417_JPRB/)

  !2.1 Polarisation definitions
  !----------------------------
  ! == pol_id +1
  !   1 average of vertical and horizontal
  !   2 nominal vertical at nadir, rotating
  !      with view angle
  !   3 nominal horizontal at nadir, rotating
  !      with view angle
  !   4 vertical
  !   5 horizontal
  !   6 + 45 minus -45 (3rd stokes vector)
  !   7 left circular - right circular (4th stokes vector)
  Integer(Kind=jpim), Dimension(7), Parameter :: npolar_compute = &
   & (/ 2, 2, 2, 1, 1, 2, 4/)
  Integer(Kind=jpim), Dimension(7), Parameter :: npolar_return = &
   & (/ 1, 1, 1, 1, 1, 2, 4/)

  ! pol_v and pol_h give proportion of v and h pol to use in emissivity calculation
  ! pol_s3 adds the 3rd/4th stokes vectors
  Real(Kind=jprb), Parameter :: pol_v(3,7) = Reshape( &
    & (/ 0.5_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 1.0_JPRB, &
       & 0.0_JPRB, 1.0_JPRB, 0.0_JPRB, &
       & 1.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB  /), (/3,7/) )
  Real(Kind=jprb), Parameter :: pol_h(3,7) = Reshape( &
    & (/ 0.5_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 1.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 1.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 1.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, 0.0_JPRB  /), (/3,7/) )
  Real(Kind=jprb), Parameter :: pol_s3(0:1,7) = Reshape( &
    & (/ 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 0.0_JPRB, &
       & 1.0_JPRB, 0.0_JPRB, &
       & 0.0_JPRB, 1.0_JPRB  /), (/2,7/) )

  !3 RTTOVSCATT aux parameters
  !---------------------------
  ! Minimum cloud cover processed by rttov_scatt
  Real(Kind=jprb), Parameter :: ccthres = 0.05_JPRB
  ! Minimum single scattering albedo processed by rttov_scatt
  Real(Kind=jprb), Parameter :: min_ssa = 1.0E-03_JPRB
  ! Rain density (g.cm-3)
  Real(Kind=jprb), Parameter :: rho_rain = 1.0_JPRB
  ! Snow density (g.cm-3)
  Real(Kind=jprb), Parameter :: rho_snow = 0.1_JPRB

  ! Flags to identify function in shared K/Adjoint routines
  Integer(Kind=jpim), Parameter :: adk_adjoint = 0
  Integer(Kind=jpim), Parameter :: adk_k       = 1

  !4 Parameters to compute refractive index of air
  !--------------------------------------------------------
  Real(Kind=jprb), Parameter :: D1   =8341.87_JPRB
  Real(Kind=jprb), Parameter :: D2   =2405955.0_JPRB
  Real(Kind=jprb), Parameter :: D3   =130.0_JPRB
  Real(Kind=jprb), Parameter :: D4   =15996.0_JPRB
  Real(Kind=jprb), Parameter :: D5   =38.9_JPRB
  Real(Kind=jprb), Parameter :: DCO2 =0.540_JPRB
  Real(Kind=jprb), Parameter :: ED1  =96095.43_JPRB
  Real(Kind=jprb), Parameter :: ED2  =0.601_JPRB
  Real(Kind=jprb), Parameter :: ED3  =0.00972_JPRB
  Real(Kind=jprb), Parameter :: ED4  =0.003661_JPRB
  Real(Kind=jprb), Parameter :: EW1  =3.7345_JPRB
  Real(Kind=jprb), Parameter :: EW2  =0.0401_JPRB
  Real(Kind=jprb), Parameter :: HTOP =100.0_JPRB
  Real(Kind=jprb), Parameter :: CTOM =1.0E-4_JPRB
  Real(Kind=jprb), Parameter :: WAVER=1700.0_JPRB

  !5 RTTOV8_M_SCATT
  !--------------------------------------------------------
  Integer(Kind=jpim), Parameter :: naer_max = 11
  Integer(Kind=jpim), Parameter :: naer_cl  = 10
  Integer(Kind=jpim), Parameter :: nhumaer(naer_max)=                     &
       & (/1,8,1,8,8,1,1,1,1,8,1/)

  Integer(Kind=jpim), Parameter :: &
        & aer_id_inso       = 1, &
        & aer_id_waso       = 2, &
        & aer_id_soot       = 3, &
        & aer_id_ssam       = 4, &
        & aer_id_sscm       = 5, &
        & aer_id_minm       = 6, &
        & aer_id_miam       = 7, &
        & aer_id_micm       = 8, &
        & aer_id_mitr       = 9, &
        & aer_id_suso       =10, &
        & aer_id_vola       =11

  Character (len=4), Parameter :: aer_name(naer_max) = &
        & (/ 'inso', &
           & 'waso', &
           & 'soot', &
           & 'ssam', &
           & 'sscm', &
           & 'minm', &
           & 'miam', &
           & 'micm', &
           & 'mitr', &
           & 'suso', &
           & 'vola' /)

  Integer(Kind=jpim), Parameter :: nwcl_max = 5
  Integer(Kind=jpim), Parameter :: nhumwcl(nwcl_max)=                     &
       & (/1,1,1,1,1/)

  Integer(Kind=jpim), Parameter :: &
        & wcl_id_stco       = 1, &
        & wcl_id_stma       = 2, &
        & wcl_id_cucc       = 3, &
        & wcl_id_cucp       = 4, &
        & wcl_id_cuma       = 5

  Character (len=4), Parameter :: wcl_name(nwcl_max) = &
        & (/ 'stco', &
           & 'stma', &
           & 'cucc', &
           & 'cucp', &
           & 'cuma' /)

  Integer(Kind=jpim), Parameter:: ncldtyp=6

  Integer(Kind=jpim), Parameter:: jpazn=11

  Real(Kind=jprb), Parameter :: E00       = 611.21_JPRB
  Real(Kind=jprb), Parameter :: T00       = 273.16_JPRB
  Real(Kind=jprb), Parameter :: TI        = T00 - 23.0_JPRB

  Real(Kind=jprb), Parameter :: min_tau = 1.0e-8_JPRB
  Real(Kind=jprb), Parameter :: min_od  = 1.0e-5_JPRB

!
! These are the RTTOV9 wavenumbers that make intervals
!
  Real(Kind=jprb), Parameter :: rttov9_wv0690_50 =  690.50_JPRB, &
                                rttov9_wv1050_00 = 1050.00_JPRB, &
                                rttov9_wv1095_25 = 1095.25_JPRB, &
                                rttov9_wv1100_25 = 1100.25_JPRB, &
                                rttov9_wv1350_25 = 1350.25_JPRB, &
                                rttov9_wv1750_25 = 1750.25_JPRB, &
                                rttov9_wv1900_25 = 1900.25_JPRB, &
                                rttov9_wv1995_00 = 1995.00_JPRB, &
                                rttov9_wv2000_00 = 2000.00_JPRB, &
                                rttov9_wv2250_00 = 2250.00_JPRB, &
                                rttov9_wv2295_25 = 2295.25_JPRB, &
                                rttov9_wv2360_00 = 2360.00_JPRB, &
                                rttov9_wv2380_25 = 2380.25_JPRB, &
                                rttov9_wv2660_25 = 2660.25_JPRB, &
                                rttov9_wv2760_25 = 2760.25_JPRB  
End Module rttov_const
