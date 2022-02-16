!
MODULE mod_iratlas
  ! Description:
  !   Data and routines for IR emissivity atlas.
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
  !  1.0      02/06/2010  Based on UW IR atlas code (E. Borbas, B. Ruston, J. Hocking)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
#include "throw.h"

  USE parkind1, Only : jpim, jprb, jplm
  
  ! Disable implicit typing
  IMPLICIT NONE

#include "rttov_errorreport.h"

  INCLUDE 'netcdf.inc'

  ! Specify kinds for NetCDF interface
  INTEGER, PARAMETER :: ncint32  = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: ncint16  = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: ncint8   = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: ncreal32 = SELECTED_REAL_KIND(6,37)
  
  ! User can specify atlas version at set-up. Version 100 is the default.
  INTEGER(KIND=jpim) :: ir_atlas_version=100   ! Version of atlas
  
  ! Atlas constants
    
  INTEGER(KIND=jpim), PARAMETER :: numpcs=6
  INTEGER(KIND=jpim), PARAMETER :: hngpnts=10
  INTEGER(KIND=jpim), PARAMETER :: numwave=416

  INTEGER(KIND=jpim), PARAMETER :: nb_lats=1800
  INTEGER(KIND=jpim), PARAMETER :: nb_lons=3600
  ! INTEGER(KIND=jpim), PARAMETER :: nb_pack=2298394  !ver2.1  2006
  INTEGER(KIND=jpim), PARAMETER :: nb_pack=2250931  !ver2.1  2007

  INTEGER(KIND=jpim), PARAMETER :: cv_lats=360
  INTEGER(KIND=jpim), PARAMETER :: cv_lons=720
  INTEGER(KIND=jpim), PARAMETER :: cv_pack=98008

  INTEGER(KIND=jpim), PARAMETER :: db_ver_year=2007
  
  INTEGER(KIND=jpim), PARAMETER :: seaice_flag=70          ! flag value returned for sea-ice
  REAL(KIND=jprb),    PARAMETER :: default_std = 0.05_JPRB ! default standard deviation

  REAL(KIND=jprb),    PARAMETER :: bfemis_gridres = 0.1_JPRB
  REAL(KIND=jprb),    PARAMETER :: bfemis_ygrid1 = 89.95_JPRB
  REAL(KIND=jprb),    PARAMETER :: bfemis_xgrid1 = -179.95_JPRB

  REAL(KIND=jprb),    PARAMETER :: cov_emis_gridres = 0.5_JPRB
  REAL(KIND=jprb),    PARAMETER :: cov_emis_ygrid1 = 89.75_JPRB
  REAL(KIND=jprb),    PARAMETER :: cov_emis_xgrid1 = -179.75_JPRB
  
  
  ! Atlas data loaded by initialisation routine
  
  INTEGER(KIND=ncint16), ALLOCATABLE :: bfemis_flag(:,:)  ! dims are (nb_lats,nb_lons)
  INTEGER(KIND=jpim),    ALLOCATABLE :: bfemis_lut(:,:)   ! dims are (nb_lats,nb_lons)
  REAL(KIND=jprb),       ALLOCATABLE :: bfemis(:,:)       ! dims are (nb_pack,hngpnts)

  INTEGER(KIND=jpim),    ALLOCATABLE :: cov_emis_lut(:,:) ! dims are (cv_lats,cv_lons)
  REAL(KIND=jprb),       ALLOCATABLE :: cov_emis(:,:)     ! dims are (cv_pack,numwave)

  REAL(KIND=jprb)     :: D(hngpnts,numpcs)
  REAL(KIND=ncreal32) :: pcu(numwave,numwave)              ! pcu need to be read in from flat file  
  
  
  ! Atlas data contained in this file
  
  REAL(KIND=jprb) :: sice_em(numwave), snow_em(numwave)
  REAL(KIND=jprb) :: sice_stdv, snow_stdv
  REAL(KIND=jprb) :: pcm(numwave)
  REAL(KIND=jprb) :: hsr_wavenum(numwave)
  REAL(KIND=jprb) :: pcu_modres(hngpnts,hngpnts)
  REAL(KIND=jprb) :: pcm_modres(hngpnts)


  DATA pcm_modres / &
        0.9782182_JPRB,   0.9759824_JPRB,   0.9681804_JPRB,   0.9187206_JPRB,   0.9265846_JPRB,  &
        0.9733831_JPRB,   0.9684200_JPRB,   0.9459322_JPRB,   0.8725729_JPRB,   0.8704424_JPRB/

!   number of PC is the first dimension and spectra is the second for pcu_modres (eigenvectors) 
  DATA pcu_modres / &
       -0.0036189_JPRB,  -0.0107326_JPRB,  -0.0357575_JPRB,  -0.0199025_JPRB,  -0.0268988_JPRB,  &
       -0.0458035_JPRB,  -0.1067330_JPRB,  -0.1013470_JPRB,  -0.0089873_JPRB,   0.0614024_JPRB,  &
       -0.0024519_JPRB,  -0.0111489_JPRB,  -0.0552112_JPRB,   0.0015246_JPRB,  -0.0248343_JPRB,  &
       -0.0227242_JPRB,  -0.0621471_JPRB,  -0.0977514_JPRB,  -0.0291845_JPRB,   0.0449263_JPRB,  &
        0.0030866_JPRB,   0.0003667_JPRB,  -0.0769359_JPRB,   0.0242423_JPRB,   0.0167816_JPRB,  &
       -0.0266578_JPRB,  -0.0145453_JPRB,  -0.0207119_JPRB,  -0.0309243_JPRB,  -0.0287100_JPRB,  &
        0.0362763_JPRB,   0.1269156_JPRB,  -0.0941633_JPRB,  -0.0524558_JPRB,   0.0621878_JPRB,  &
       -0.0331482_JPRB,   0.0164898_JPRB,  -0.0062478_JPRB,   0.0276787_JPRB,   0.0014637_JPRB,  &
        0.0393517_JPRB,   0.1177192_JPRB,  -0.0421300_JPRB,  -0.1202559_JPRB,  -0.0664961_JPRB,  &
        0.0741313_JPRB,   0.0287669_JPRB,   0.0501746_JPRB,   0.0506316_JPRB,   0.0914001_JPRB,  &
        0.0002204_JPRB,  -0.0087577_JPRB,  -0.0449563_JPRB,   0.0327884_JPRB,  -0.0167350_JPRB,  &
        0.0526362_JPRB,  -0.0313638_JPRB,   0.0499965_JPRB,   0.0695323_JPRB,   0.0616089_JPRB,  &
        0.0112971_JPRB,  -0.0362747_JPRB,  -0.0328459_JPRB,   0.0031088_JPRB,  -0.0465960_JPRB,  &
       -0.0110655_JPRB,  -0.0633873_JPRB,   0.0473527_JPRB,   0.0344864_JPRB,  -0.0452882_JPRB,  &
        0.0207165_JPRB,  -0.0652447_JPRB,  -0.0579818_JPRB,   0.0071345_JPRB,  -0.0130217_JPRB,  &
       -0.0274127_JPRB,   0.0848268_JPRB,   0.0065488_JPRB,   0.0052432_JPRB,  -0.0071355_JPRB,  &
        0.0836001_JPRB,  -0.0290997_JPRB,  -0.0019292_JPRB,  -0.0057489_JPRB,   0.0253799_JPRB,  &
        0.0964016_JPRB,   0.0197362_JPRB,   0.0286795_JPRB,  -0.1536087_JPRB,   0.0228268_JPRB,  &
        0.0887494_JPRB,   0.0377112_JPRB,   0.0174746_JPRB,   0.1030488_JPRB,  -0.0500023_JPRB,  &
        0.0222185_JPRB,   0.0446731_JPRB,  -0.0709027_JPRB,   0.0504857_JPRB,   0.0019073_JPRB/

  DATA hsr_wavenum / &
    699.3_JPRB,  704.3_JPRB,  709.3_JPRB,  714.3_JPRB,  719.3_JPRB,  724.3_JPRB,  729.3_JPRB,  734.3_JPRB,  &
    739.3_JPRB,  744.3_JPRB,  749.3_JPRB,  754.3_JPRB,  759.3_JPRB,  764.3_JPRB,  769.3_JPRB,  774.3_JPRB,  &
    779.3_JPRB,  784.3_JPRB,  789.3_JPRB,  794.3_JPRB,  799.3_JPRB,  804.3_JPRB,  809.3_JPRB,  814.3_JPRB,  &
    819.3_JPRB,  824.3_JPRB,  829.3_JPRB,  834.3_JPRB,  839.3_JPRB,  844.3_JPRB,  849.3_JPRB,  854.3_JPRB,  &
    859.3_JPRB,  864.3_JPRB,  869.3_JPRB,  874.3_JPRB,  879.3_JPRB,  884.3_JPRB,  889.3_JPRB,  894.3_JPRB,  &
    899.3_JPRB,  904.3_JPRB,  909.3_JPRB,  914.3_JPRB,  919.3_JPRB,  924.3_JPRB,  929.3_JPRB,  934.3_JPRB,  &
    939.3_JPRB,  944.3_JPRB,  949.3_JPRB,  954.3_JPRB,  959.3_JPRB,  964.3_JPRB,  969.3_JPRB,  974.3_JPRB,  &
    979.3_JPRB,  984.3_JPRB,  989.3_JPRB,  994.3_JPRB,  999.3_JPRB, 1004.3_JPRB, 1009.3_JPRB, 1014.3_JPRB,  &
    1019.3_JPRB, 1024.3_JPRB, 1029.3_JPRB, 1034.3_JPRB, 1039.3_JPRB, 1044.3_JPRB, 1049.3_JPRB, 1054.3_JPRB,  &
    1059.3_JPRB, 1064.3_JPRB, 1069.3_JPRB, 1074.3_JPRB, 1079.3_JPRB, 1084.3_JPRB, 1089.3_JPRB, 1094.3_JPRB,  &
    1099.3_JPRB, 1104.3_JPRB, 1109.3_JPRB, 1114.3_JPRB, 1119.3_JPRB, 1124.3_JPRB, 1129.3_JPRB, 1134.3_JPRB,  &
    1139.3_JPRB, 1144.3_JPRB, 1149.3_JPRB, 1154.3_JPRB, 1159.3_JPRB, 1164.3_JPRB, 1169.3_JPRB, 1174.3_JPRB,  &
    1179.3_JPRB, 1184.3_JPRB, 1189.3_JPRB, 1194.3_JPRB, 1199.3_JPRB, 1204.3_JPRB, 1209.3_JPRB, 1214.3_JPRB,  &
    1219.3_JPRB, 1224.3_JPRB, 1229.3_JPRB, 1234.3_JPRB, 1239.3_JPRB, 1244.3_JPRB, 1249.3_JPRB, 1254.3_JPRB,  &
    1259.3_JPRB, 1264.3_JPRB, 1269.3_JPRB, 1274.3_JPRB, 1279.3_JPRB, 1284.3_JPRB, 1289.3_JPRB, 1294.3_JPRB,  &
    1299.3_JPRB, 1304.3_JPRB, 1309.3_JPRB, 1314.3_JPRB, 1319.3_JPRB, 1324.3_JPRB, 1329.3_JPRB, 1334.3_JPRB,  &
    1339.3_JPRB, 1344.3_JPRB, 1349.3_JPRB, 1354.3_JPRB, 1359.3_JPRB, 1364.3_JPRB, 1369.3_JPRB, 1374.3_JPRB,  &
    1379.3_JPRB, 1384.3_JPRB, 1389.3_JPRB, 1394.3_JPRB, 1399.3_JPRB, 1404.3_JPRB, 1409.3_JPRB, 1414.3_JPRB,  &
    1419.3_JPRB, 1424.3_JPRB, 1429.3_JPRB, 1434.3_JPRB, 1439.3_JPRB, 1444.3_JPRB, 1449.3_JPRB, 1454.3_JPRB,  &
    1459.3_JPRB, 1464.3_JPRB, 1469.3_JPRB, 1474.3_JPRB, 1479.3_JPRB, 1484.3_JPRB, 1489.3_JPRB, 1494.3_JPRB,  &
    1499.3_JPRB, 1504.3_JPRB, 1509.3_JPRB, 1514.3_JPRB, 1519.3_JPRB, 1524.3_JPRB, 1529.3_JPRB, 1534.3_JPRB,  &
    1539.3_JPRB, 1544.3_JPRB, 1549.3_JPRB, 1554.3_JPRB, 1559.3_JPRB, 1564.3_JPRB, 1569.3_JPRB, 1574.3_JPRB,  &
    1579.3_JPRB, 1584.3_JPRB, 1589.3_JPRB, 1594.3_JPRB, 1599.3_JPRB, 1604.3_JPRB, 1609.3_JPRB, 1614.3_JPRB,  &
    1619.3_JPRB, 1624.3_JPRB, 1629.3_JPRB, 1634.3_JPRB, 1639.3_JPRB, 1644.3_JPRB, 1649.3_JPRB, 1654.3_JPRB,  &
    1659.3_JPRB, 1664.3_JPRB, 1669.3_JPRB, 1674.3_JPRB, 1679.3_JPRB, 1684.3_JPRB, 1689.3_JPRB, 1694.3_JPRB,  &
    1699.3_JPRB, 1704.3_JPRB, 1709.3_JPRB, 1714.3_JPRB, 1719.3_JPRB, 1724.3_JPRB, 1729.3_JPRB, 1734.3_JPRB,  &
    1739.3_JPRB, 1744.3_JPRB, 1749.3_JPRB, 1754.3_JPRB, 1759.3_JPRB, 1764.3_JPRB, 1769.3_JPRB, 1774.3_JPRB,  &
    1779.3_JPRB, 1784.3_JPRB, 1789.3_JPRB, 1794.3_JPRB, 1799.3_JPRB, 1804.3_JPRB, 1809.3_JPRB, 1814.3_JPRB,  &
    1819.3_JPRB, 1824.3_JPRB, 1829.3_JPRB, 1834.3_JPRB, 1839.3_JPRB, 1844.3_JPRB, 1849.3_JPRB, 1854.3_JPRB,  &
    1859.3_JPRB, 1864.3_JPRB, 1869.3_JPRB, 1874.3_JPRB, 1879.3_JPRB, 1884.3_JPRB, 1889.3_JPRB, 1894.3_JPRB,  &
    1899.3_JPRB, 1904.3_JPRB, 1909.3_JPRB, 1914.3_JPRB, 1919.3_JPRB, 1924.3_JPRB, 1929.3_JPRB, 1934.3_JPRB,  &
    1939.3_JPRB, 1944.3_JPRB, 1949.3_JPRB, 1954.3_JPRB, 1959.3_JPRB, 1964.3_JPRB, 1969.3_JPRB, 1974.3_JPRB,  &
    1979.3_JPRB, 1984.3_JPRB, 1989.3_JPRB, 1994.3_JPRB, 1999.3_JPRB, 2004.3_JPRB, 2009.3_JPRB, 2014.3_JPRB,  &
    2019.3_JPRB, 2024.3_JPRB, 2029.3_JPRB, 2034.3_JPRB, 2039.3_JPRB, 2044.3_JPRB, 2049.3_JPRB, 2054.3_JPRB,  &
    2059.3_JPRB, 2064.3_JPRB, 2069.3_JPRB, 2074.3_JPRB, 2079.3_JPRB, 2084.3_JPRB, 2089.3_JPRB, 2094.3_JPRB,  &
    2099.3_JPRB, 2104.3_JPRB, 2109.3_JPRB, 2114.3_JPRB, 2119.3_JPRB, 2124.3_JPRB, 2129.3_JPRB, 2134.3_JPRB,  &
    2139.3_JPRB, 2144.3_JPRB, 2149.3_JPRB, 2154.3_JPRB, 2159.3_JPRB, 2164.3_JPRB, 2169.3_JPRB, 2174.3_JPRB,  &
    2179.3_JPRB, 2184.3_JPRB, 2189.3_JPRB, 2194.3_JPRB, 2199.3_JPRB, 2204.3_JPRB, 2209.3_JPRB, 2214.3_JPRB,  &
    2219.3_JPRB, 2224.3_JPRB, 2229.3_JPRB, 2234.3_JPRB, 2239.3_JPRB, 2244.3_JPRB, 2249.3_JPRB, 2254.3_JPRB,  &
    2259.3_JPRB, 2264.3_JPRB, 2269.3_JPRB, 2274.3_JPRB, 2279.3_JPRB, 2284.3_JPRB, 2289.3_JPRB, 2294.3_JPRB,  &
    2299.3_JPRB, 2304.3_JPRB, 2309.3_JPRB, 2314.3_JPRB, 2319.3_JPRB, 2324.3_JPRB, 2329.3_JPRB, 2334.3_JPRB,  &
    2339.3_JPRB, 2344.3_JPRB, 2349.3_JPRB, 2354.3_JPRB, 2359.3_JPRB, 2364.3_JPRB, 2369.3_JPRB, 2374.3_JPRB,  &
    2379.3_JPRB, 2384.3_JPRB, 2389.3_JPRB, 2394.3_JPRB, 2399.3_JPRB, 2404.3_JPRB, 2409.3_JPRB, 2414.3_JPRB,  &
    2419.3_JPRB, 2424.3_JPRB, 2429.3_JPRB, 2434.3_JPRB, 2439.3_JPRB, 2444.3_JPRB, 2449.3_JPRB, 2454.3_JPRB,  &
    2459.3_JPRB, 2464.3_JPRB, 2469.3_JPRB, 2474.3_JPRB, 2479.3_JPRB, 2484.3_JPRB, 2489.3_JPRB, 2494.3_JPRB,  &
    2499.3_JPRB, 2504.3_JPRB, 2509.3_JPRB, 2514.3_JPRB, 2519.3_JPRB, 2524.3_JPRB, 2529.3_JPRB, 2534.3_JPRB,  &
    2539.3_JPRB, 2544.3_JPRB, 2549.3_JPRB, 2554.3_JPRB, 2559.3_JPRB, 2564.3_JPRB, 2569.3_JPRB, 2574.3_JPRB,  &
    2579.3_JPRB, 2584.3_JPRB, 2589.3_JPRB, 2594.3_JPRB, 2599.3_JPRB, 2604.3_JPRB, 2609.3_JPRB, 2614.3_JPRB,  &
    2619.3_JPRB, 2624.3_JPRB, 2629.3_JPRB, 2634.3_JPRB, 2639.3_JPRB, 2644.3_JPRB, 2649.3_JPRB, 2654.3_JPRB,  &
    2659.3_JPRB, 2664.3_JPRB, 2669.3_JPRB, 2674.3_JPRB, 2679.3_JPRB, 2684.3_JPRB, 2689.3_JPRB, 2694.3_JPRB,  &
    2699.3_JPRB, 2704.3_JPRB, 2709.3_JPRB, 2714.3_JPRB, 2719.3_JPRB, 2724.3_JPRB, 2729.3_JPRB, 2734.3_JPRB,  &
    2739.3_JPRB, 2744.3_JPRB, 2749.3_JPRB, 2754.3_JPRB, 2759.3_JPRB, 2764.3_JPRB, 2769.3_JPRB, 2774.3_JPRB /

  DATA pcm / &
        0.9782182_JPRB,   0.9770744_JPRB,   0.9763290_JPRB,   0.9763215_JPRB,   0.9760258_JPRB,  &
        0.9763704_JPRB,   0.9767076_JPRB,   0.9763077_JPRB,   0.9758835_JPRB,   0.9753462_JPRB,  &
        0.9748067_JPRB,   0.9734465_JPRB,   0.9721510_JPRB,   0.9717180_JPRB,   0.9714773_JPRB,  &
        0.9706340_JPRB,   0.9710826_JPRB,   0.9722888_JPRB,   0.9731166_JPRB,   0.9732918_JPRB,  &
        0.9736975_JPRB,   0.9751787_JPRB,   0.9770049_JPRB,   0.9773170_JPRB,   0.9765164_JPRB,  &
        0.9759824_JPRB,   0.9750199_JPRB,   0.9746831_JPRB,   0.9738413_JPRB,   0.9731615_JPRB,  &
        0.9720387_JPRB,   0.9716908_JPRB,   0.9708628_JPRB,   0.9705366_JPRB,   0.9697853_JPRB,  &
        0.9694459_JPRB,   0.9688896_JPRB,   0.9688236_JPRB,   0.9689180_JPRB,   0.9692774_JPRB,  &
        0.9693237_JPRB,   0.9692513_JPRB,   0.9689918_JPRB,   0.9686664_JPRB,   0.9684489_JPRB,  &
        0.9681804_JPRB,   0.9672847_JPRB,   0.9667084_JPRB,   0.9661347_JPRB,   0.9655386_JPRB,  &
        0.9650131_JPRB,   0.9641176_JPRB,   0.9628995_JPRB,   0.9620982_JPRB,   0.9605948_JPRB,  &
        0.9590283_JPRB,   0.9572537_JPRB,   0.9552648_JPRB,   0.9529146_JPRB,   0.9505763_JPRB,  &
        0.9486620_JPRB,   0.9468448_JPRB,   0.9446425_JPRB,   0.9428397_JPRB,   0.9415421_JPRB,  &
        0.9398234_JPRB,   0.9378662_JPRB,   0.9358756_JPRB,   0.9338515_JPRB,   0.9317511_JPRB,  &
        0.9296144_JPRB,   0.9274116_JPRB,   0.9248639_JPRB,   0.9219664_JPRB,   0.9197029_JPRB,  &
        0.9187206_JPRB,   0.9195539_JPRB,   0.9211251_JPRB,   0.9227578_JPRB,   0.9242273_JPRB,  &
        0.9256495_JPRB,   0.9265392_JPRB,   0.9276078_JPRB,   0.9279289_JPRB,   0.9282181_JPRB,  &
        0.9284544_JPRB,   0.9289097_JPRB,   0.9299400_JPRB,   0.9314128_JPRB,   0.9329405_JPRB,  &
        0.9349486_JPRB,   0.9377099_JPRB,   0.9380918_JPRB,   0.9354525_JPRB,   0.9330018_JPRB,  &
        0.9316696_JPRB,   0.9308965_JPRB,   0.9296793_JPRB,   0.9282659_JPRB,   0.9273711_JPRB,  &
        0.9268156_JPRB,   0.9265846_JPRB,   0.9264724_JPRB,   0.9278417_JPRB,   0.9298262_JPRB,  &
        0.9342009_JPRB,   0.9397170_JPRB,   0.9451398_JPRB,   0.9501663_JPRB,   0.9547508_JPRB,  &
        0.9586911_JPRB,   0.9618842_JPRB,   0.9649577_JPRB,   0.9675525_JPRB,   0.9696881_JPRB,  &
        0.9708689_JPRB,   0.9717879_JPRB,   0.9722518_JPRB,   0.9724457_JPRB,   0.9728941_JPRB,  &
        0.9731293_JPRB,   0.9731925_JPRB,   0.9730867_JPRB,   0.9733831_JPRB,   0.9735166_JPRB,  &
        0.9740434_JPRB,   0.9742066_JPRB,   0.9746855_JPRB,   0.9748268_JPRB,   0.9749292_JPRB,  &
        0.9751188_JPRB,   0.9752902_JPRB,   0.9751062_JPRB,   0.9751985_JPRB,   0.9752622_JPRB,  &
        0.9750626_JPRB,   0.9755121_JPRB,   0.9755228_JPRB,   0.9760818_JPRB,   0.9759580_JPRB,  &
        0.9758280_JPRB,   0.9755163_JPRB,   0.9754220_JPRB,   0.9750829_JPRB,   0.9743836_JPRB,  &
        0.9745844_JPRB,   0.9742978_JPRB,   0.9740397_JPRB,   0.9744191_JPRB,   0.9745796_JPRB,  &
        0.9749123_JPRB,   0.9750853_JPRB,   0.9746974_JPRB,   0.9747824_JPRB,   0.9746920_JPRB,  &
        0.9735873_JPRB,   0.9733123_JPRB,   0.9725510_JPRB,   0.9718717_JPRB,   0.9713586_JPRB,  &
        0.9706160_JPRB,   0.9701124_JPRB,   0.9698699_JPRB,   0.9698430_JPRB,   0.9694992_JPRB,  &
        0.9691019_JPRB,   0.9690002_JPRB,   0.9678345_JPRB,   0.9668854_JPRB,   0.9659764_JPRB,  &
        0.9666998_JPRB,   0.9669611_JPRB,   0.9665817_JPRB,   0.9679645_JPRB,   0.9695909_JPRB,  &
        0.9711555_JPRB,   0.9724632_JPRB,   0.9737635_JPRB,   0.9746142_JPRB,   0.9748497_JPRB,  &
        0.9752109_JPRB,   0.9752749_JPRB,   0.9754022_JPRB,   0.9753313_JPRB,   0.9746057_JPRB,  &
        0.9745884_JPRB,   0.9747860_JPRB,   0.9752877_JPRB,   0.9753085_JPRB,   0.9759305_JPRB,  &
        0.9752344_JPRB,   0.9748027_JPRB,   0.9757417_JPRB,   0.9751943_JPRB,   0.9748128_JPRB,  &
        0.9743713_JPRB,   0.9741939_JPRB,   0.9725359_JPRB,   0.9723988_JPRB,   0.9716700_JPRB,  &
        0.9708291_JPRB,   0.9705051_JPRB,   0.9699901_JPRB,   0.9689955_JPRB,   0.9683419_JPRB,  &
        0.9684200_JPRB,   0.9672046_JPRB,   0.9660766_JPRB,   0.9658424_JPRB,   0.9648336_JPRB,  &
        0.9640325_JPRB,   0.9642861_JPRB,   0.9636880_JPRB,   0.9638920_JPRB,   0.9638573_JPRB,  &
        0.9641714_JPRB,   0.9648057_JPRB,   0.9648220_JPRB,   0.9639065_JPRB,   0.9635883_JPRB,  &
        0.9626419_JPRB,   0.9616417_JPRB,   0.9600965_JPRB,   0.9587714_JPRB,   0.9576451_JPRB,  &
        0.9557189_JPRB,   0.9545730_JPRB,   0.9550443_JPRB,   0.9551759_JPRB,   0.9560625_JPRB,  &
        0.9576327_JPRB,   0.9587138_JPRB,   0.9594474_JPRB,   0.9598546_JPRB,   0.9601094_JPRB,  &
        0.9601356_JPRB,   0.9597549_JPRB,   0.9590299_JPRB,   0.9581512_JPRB,   0.9572046_JPRB,  &
        0.9557602_JPRB,   0.9538486_JPRB,   0.9521495_JPRB,   0.9503905_JPRB,   0.9491790_JPRB,  &
        0.9485527_JPRB,   0.9479896_JPRB,   0.9475234_JPRB,   0.9468080_JPRB,   0.9469628_JPRB,  &
        0.9469683_JPRB,   0.9465806_JPRB,   0.9468755_JPRB,   0.9466828_JPRB,   0.9471480_JPRB,  &
        0.9470276_JPRB,   0.9470209_JPRB,   0.9468378_JPRB,   0.9464890_JPRB,   0.9462101_JPRB,  &
        0.9459322_JPRB,   0.9449111_JPRB,   0.9435923_JPRB,   0.9416961_JPRB,   0.9401403_JPRB,  &
        0.9387150_JPRB,   0.9374595_JPRB,   0.9347988_JPRB,   0.9319339_JPRB,   0.9295776_JPRB,  &
        0.9268476_JPRB,   0.9243815_JPRB,   0.9224647_JPRB,   0.9208075_JPRB,   0.9195780_JPRB,  &
        0.9183103_JPRB,   0.9171674_JPRB,   0.9164810_JPRB,   0.9160877_JPRB,   0.9151877_JPRB,  &
        0.9148492_JPRB,   0.9142842_JPRB,   0.9142084_JPRB,   0.9138089_JPRB,   0.9137760_JPRB,  &
        0.9137531_JPRB,   0.9141592_JPRB,   0.9136598_JPRB,   0.9125727_JPRB,   0.9108481_JPRB,  &
        0.9093652_JPRB,   0.9080561_JPRB,   0.9062355_JPRB,   0.9046820_JPRB,   0.9028210_JPRB,  &
        0.9018152_JPRB,   0.9008504_JPRB,   0.9000632_JPRB,   0.8995758_JPRB,   0.8989593_JPRB,  &
        0.8987811_JPRB,   0.8992507_JPRB,   0.8999549_JPRB,   0.9013391_JPRB,   0.9020863_JPRB,  &
        0.9025120_JPRB,   0.9023982_JPRB,   0.9015658_JPRB,   0.9008633_JPRB,   0.8996401_JPRB,  &
        0.8981582_JPRB,   0.8969440_JPRB,   0.8946483_JPRB,   0.8925536_JPRB,   0.8906261_JPRB,  &
        0.8889833_JPRB,   0.8870751_JPRB,   0.8845615_JPRB,   0.8825631_JPRB,   0.8811586_JPRB,  &
        0.8796447_JPRB,   0.8779839_JPRB,   0.8765292_JPRB,   0.8754975_JPRB,   0.8739760_JPRB,  &
        0.8725729_JPRB,   0.8714029_JPRB,   0.8706908_JPRB,   0.8710466_JPRB,   0.8699325_JPRB,  &
        0.8697992_JPRB,   0.8718969_JPRB,   0.8713725_JPRB,   0.8701416_JPRB,   0.8695096_JPRB,  &
        0.8698574_JPRB,   0.8700698_JPRB,   0.8694080_JPRB,   0.8693934_JPRB,   0.8693246_JPRB,  &
        0.8698239_JPRB,   0.8696592_JPRB,   0.8681608_JPRB,   0.8656288_JPRB,   0.8654716_JPRB,  &
        0.8640761_JPRB,   0.8639477_JPRB,   0.8635154_JPRB,   0.8630069_JPRB,   0.8623275_JPRB,  &
        0.8623751_JPRB,   0.8627441_JPRB,   0.8630516_JPRB,   0.8638958_JPRB,   0.8644919_JPRB,  &
        0.8655882_JPRB,   0.8666160_JPRB,   0.8676174_JPRB,   0.8692035_JPRB,   0.8695340_JPRB,  &
        0.8703975_JPRB,   0.8714244_JPRB,   0.8715467_JPRB,   0.8713564_JPRB,   0.8712272_JPRB,  &
        0.8714187_JPRB,   0.8701625_JPRB,   0.8697796_JPRB,   0.8688766_JPRB,   0.8682391_JPRB,  &
        0.8680181_JPRB,   0.8676605_JPRB,   0.8672657_JPRB,   0.8679592_JPRB,   0.8675538_JPRB,  &
        0.8686572_JPRB,   0.8682060_JPRB,   0.8688578_JPRB,   0.8693632_JPRB,   0.8689557_JPRB,  &
        0.8681611_JPRB,   0.8684876_JPRB,   0.8680010_JPRB,   0.8675498_JPRB,   0.8675414_JPRB,  &
        0.8677824_JPRB,   0.8665875_JPRB,   0.8668503_JPRB,   0.8665696_JPRB,   0.8671130_JPRB,  &
        0.8669835_JPRB,   0.8671956_JPRB,   0.8683699_JPRB,   0.8685648_JPRB,   0.8682314_JPRB,  &
        0.8683055_JPRB,   0.8694246_JPRB,   0.8689486_JPRB,   0.8693868_JPRB,   0.8694460_JPRB,  &
        0.8701811_JPRB,   0.8704424_JPRB,   0.8709887_JPRB,   0.8712862_JPRB,   0.8721344_JPRB,  &
        0.8724745_JPRB,   0.8727338_JPRB,   0.8740577_JPRB,   0.8748575_JPRB,   0.8747587_JPRB,  &
        0.8762293_JPRB,   0.8772818_JPRB,   0.8779803_JPRB,   0.8791369_JPRB,   0.8807610_JPRB,  &
        0.8813813_JPRB/

  DATA sice_stdv / 0.015_JPRB /
  DATA sice_em / &
      0.9370_JPRB, 0.9370_JPRB, 0.9370_JPRB, 0.9370_JPRB, 0.9367_JPRB, 0.9367_JPRB, 0.9366_JPRB, 0.9365_JPRB, &
      0.9365_JPRB, 0.9365_JPRB, 0.9365_JPRB, 0.9366_JPRB, 0.9367_JPRB, 0.9370_JPRB, 0.9374_JPRB, 0.9381_JPRB, &
      0.9386_JPRB, 0.9393_JPRB, 0.9401_JPRB, 0.9408_JPRB, 0.9415_JPRB, 0.9427_JPRB, 0.9440_JPRB, 0.9452_JPRB, &
      0.9464_JPRB, 0.9481_JPRB, 0.9496_JPRB, 0.9511_JPRB, 0.9525_JPRB, 0.9544_JPRB, 0.9563_JPRB, 0.9582_JPRB, &
      0.9602_JPRB, 0.9620_JPRB, 0.9640_JPRB, 0.9658_JPRB, 0.9678_JPRB, 0.9702_JPRB, 0.9725_JPRB, 0.9748_JPRB, &
      0.9770_JPRB, 0.9792_JPRB, 0.9814_JPRB, 0.9836_JPRB, 0.9856_JPRB, 0.9872_JPRB, 0.9885_JPRB, 0.9897_JPRB, &
      0.9905_JPRB, 0.9911_JPRB, 0.9913_JPRB, 0.9913_JPRB, 0.9912_JPRB, 0.9910_JPRB, 0.9907_JPRB, 0.9904_JPRB, &
      0.9901_JPRB, 0.9897_JPRB, 0.9893_JPRB, 0.9889_JPRB, 0.9885_JPRB, 0.9880_JPRB, 0.9876_JPRB, 0.9871_JPRB, &
      0.9867_JPRB, 0.9864_JPRB, 0.9861_JPRB, 0.9858_JPRB, 0.9854_JPRB, 0.9852_JPRB, 0.9849_JPRB, 0.9846_JPRB, &
      0.9844_JPRB, 0.9842_JPRB, 0.9840_JPRB, 0.9838_JPRB, 0.9836_JPRB, 0.9834_JPRB, 0.9832_JPRB, 0.9831_JPRB, &
      0.9829_JPRB, 0.9828_JPRB, 0.9826_JPRB, 0.9824_JPRB, 0.9822_JPRB, 0.9821_JPRB, 0.9820_JPRB, 0.9819_JPRB, &
      0.9817_JPRB, 0.9816_JPRB, 0.9814_JPRB, 0.9813_JPRB, 0.9811_JPRB, 0.9810_JPRB, 0.9808_JPRB, 0.9807_JPRB, &
      0.9805_JPRB, 0.9804_JPRB, 0.9803_JPRB, 0.9801_JPRB, 0.9799_JPRB, 0.9797_JPRB, 0.9796_JPRB, 0.9794_JPRB, &
      0.9792_JPRB, 0.9791_JPRB, 0.9789_JPRB, 0.9787_JPRB, 0.9786_JPRB, 0.9785_JPRB, 0.9784_JPRB, 0.9784_JPRB, &
      0.9783_JPRB, 0.9782_JPRB, 0.9782_JPRB, 0.9782_JPRB, 0.9781_JPRB, 0.9781_JPRB, 0.9781_JPRB, 0.9781_JPRB, &
      0.9781_JPRB, 0.9781_JPRB, 0.9781_JPRB, 0.9780_JPRB, 0.9780_JPRB, 0.9780_JPRB, 0.9780_JPRB, 0.9780_JPRB, &
      0.9780_JPRB, 0.9780_JPRB, 0.9779_JPRB, 0.9779_JPRB, 0.9778_JPRB, 0.9778_JPRB, 0.9777_JPRB, 0.9777_JPRB, &
      0.9777_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, &
      0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9775_JPRB, 0.9775_JPRB, 0.9775_JPRB, &
      0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9777_JPRB, &
      0.9777_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9776_JPRB, 0.9776_JPRB, 0.9776_JPRB, &
      0.9775_JPRB, 0.9775_JPRB, 0.9774_JPRB, 0.9773_JPRB, 0.9773_JPRB, 0.9773_JPRB, 0.9773_JPRB, 0.9773_JPRB, &
      0.9774_JPRB, 0.9774_JPRB, 0.9775_JPRB, 0.9776_JPRB, 0.9777_JPRB, 0.9778_JPRB, 0.9779_JPRB, 0.9780_JPRB, &
      0.9781_JPRB, 0.9782_JPRB, 0.9783_JPRB, 0.9785_JPRB, 0.9786_JPRB, 0.9788_JPRB, 0.9790_JPRB, 0.9792_JPRB, &
      0.9793_JPRB, 0.9795_JPRB, 0.9797_JPRB, 0.9799_JPRB, 0.9801_JPRB, 0.9802_JPRB, 0.9803_JPRB, 0.9805_JPRB, &
      0.9806_JPRB, 0.9807_JPRB, 0.9808_JPRB, 0.9809_JPRB, 0.9810_JPRB, 0.9811_JPRB, 0.9811_JPRB, 0.9811_JPRB, &
      0.9810_JPRB, 0.9810_JPRB, 0.9810_JPRB, 0.9809_JPRB, 0.9808_JPRB, 0.9808_JPRB, 0.9807_JPRB, 0.9807_JPRB, &
      0.9806_JPRB, 0.9805_JPRB, 0.9805_JPRB, 0.9804_JPRB, 0.9803_JPRB, 0.9802_JPRB, 0.9802_JPRB, 0.9801_JPRB, &
      0.9800_JPRB, 0.9799_JPRB, 0.9798_JPRB, 0.9797_JPRB, 0.9797_JPRB, 0.9795_JPRB, 0.9795_JPRB, 0.9794_JPRB, &
      0.9793_JPRB, 0.9792_JPRB, 0.9791_JPRB, 0.9791_JPRB, 0.9789_JPRB, 0.9789_JPRB, 0.9788_JPRB, 0.9787_JPRB, &
      0.9786_JPRB, 0.9785_JPRB, 0.9785_JPRB, 0.9783_JPRB, 0.9783_JPRB, 0.9782_JPRB, 0.9781_JPRB, 0.9781_JPRB, &
      0.9780_JPRB, 0.9779_JPRB, 0.9779_JPRB, 0.9778_JPRB, 0.9777_JPRB, 0.9777_JPRB, 0.9776_JPRB, 0.9775_JPRB, &
      0.9774_JPRB, 0.9774_JPRB, 0.9773_JPRB, 0.9772_JPRB, 0.9771_JPRB, 0.9771_JPRB, 0.9770_JPRB, 0.9769_JPRB, &
      0.9769_JPRB, 0.9768_JPRB, 0.9767_JPRB, 0.9766_JPRB, 0.9765_JPRB, 0.9765_JPRB, 0.9764_JPRB, 0.9764_JPRB, &
      0.9763_JPRB, 0.9762_JPRB, 0.9762_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9760_JPRB, 0.9759_JPRB, 0.9758_JPRB, &
      0.9757_JPRB, 0.9757_JPRB, 0.9756_JPRB, 0.9756_JPRB, 0.9755_JPRB, 0.9755_JPRB, 0.9754_JPRB, 0.9754_JPRB, &
      0.9754_JPRB, 0.9754_JPRB, 0.9754_JPRB, 0.9753_JPRB, 0.9753_JPRB, 0.9753_JPRB, 0.9753_JPRB, 0.9752_JPRB, &
      0.9752_JPRB, 0.9752_JPRB, 0.9753_JPRB, 0.9754_JPRB, 0.9755_JPRB, 0.9756_JPRB, 0.9757_JPRB, 0.9757_JPRB, &
      0.9758_JPRB, 0.9758_JPRB, 0.9759_JPRB, 0.9759_JPRB, 0.9760_JPRB, 0.9760_JPRB, 0.9760_JPRB, 0.9760_JPRB, &
      0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, 0.9761_JPRB, &
      0.9761_JPRB, 0.9760_JPRB, 0.9760_JPRB, 0.9759_JPRB, 0.9759_JPRB, 0.9758_JPRB, 0.9758_JPRB, 0.9757_JPRB, &
      0.9757_JPRB, 0.9757_JPRB, 0.9756_JPRB, 0.9756_JPRB, 0.9755_JPRB, 0.9754_JPRB, 0.9753_JPRB, 0.9753_JPRB, &
      0.9752_JPRB, 0.9751_JPRB, 0.9751_JPRB, 0.9750_JPRB, 0.9750_JPRB, 0.9749_JPRB, 0.9749_JPRB, 0.9748_JPRB, &
      0.9747_JPRB, 0.9746_JPRB, 0.9746_JPRB, 0.9746_JPRB, 0.9745_JPRB, 0.9744_JPRB, 0.9743_JPRB, 0.9742_JPRB, &
      0.9742_JPRB, 0.9741_JPRB, 0.9740_JPRB, 0.9739_JPRB, 0.9739_JPRB, 0.9739_JPRB, 0.9738_JPRB, 0.9737_JPRB, &
      0.9736_JPRB, 0.9735_JPRB, 0.9735_JPRB, 0.9734_JPRB, 0.9733_JPRB, 0.9732_JPRB, 0.9731_JPRB, 0.9731_JPRB, &
      0.9730_JPRB, 0.9729_JPRB, 0.9728_JPRB, 0.9727_JPRB, 0.9726_JPRB, 0.9725_JPRB, 0.9724_JPRB, 0.9723_JPRB, &
      0.9723_JPRB, 0.9722_JPRB, 0.9721_JPRB, 0.9720_JPRB, 0.9719_JPRB, 0.9718_JPRB, 0.9717_JPRB, 0.9716_JPRB, &
      0.9715_JPRB, 0.9714_JPRB, 0.9713_JPRB, 0.9712_JPRB, 0.9711_JPRB, 0.9709_JPRB, 0.9708_JPRB, 0.9706_JPRB, &
      0.9705_JPRB, 0.9704_JPRB, 0.9703_JPRB, 0.9702_JPRB, 0.9700_JPRB, 0.9699_JPRB, 0.9698_JPRB, 0.9696_JPRB, &
      0.9695_JPRB, 0.9693_JPRB, 0.9691_JPRB, 0.9690_JPRB, 0.9688_JPRB, 0.9686_JPRB, 0.9685_JPRB, 0.9683_JPRB, &
      0.9682_JPRB, 0.9681_JPRB, 0.9679_JPRB, 0.9677_JPRB, 0.9676_JPRB, 0.9674_JPRB, 0.9671_JPRB, 0.9669_JPRB/

  DATA snow_stdv / 0.015_JPRB /
  DATA snow_em / &
      0.9716_JPRB, 0.9716_JPRB, 0.9716_JPRB, 0.9716_JPRB, 0.9713_JPRB, 0.9710_JPRB, 0.9708_JPRB, 0.9706_JPRB, &
      0.9705_JPRB, 0.9705_JPRB, 0.9705_JPRB, 0.9703_JPRB, 0.9701_JPRB, 0.9700_JPRB, 0.9699_JPRB, 0.9700_JPRB, &
      0.9702_JPRB, 0.9703_JPRB, 0.9705_JPRB, 0.9707_JPRB, 0.9710_JPRB, 0.9714_JPRB, 0.9717_JPRB, 0.9722_JPRB, &
      0.9728_JPRB, 0.9734_JPRB, 0.9740_JPRB, 0.9746_JPRB, 0.9753_JPRB, 0.9759_JPRB, 0.9765_JPRB, 0.9771_JPRB, &
      0.9778_JPRB, 0.9784_JPRB, 0.9792_JPRB, 0.9798_JPRB, 0.9806_JPRB, 0.9814_JPRB, 0.9824_JPRB, 0.9833_JPRB, &
      0.9842_JPRB, 0.9852_JPRB, 0.9863_JPRB, 0.9873_JPRB, 0.9882_JPRB, 0.9891_JPRB, 0.9901_JPRB, 0.9908_JPRB, &
      0.9914_JPRB, 0.9920_JPRB, 0.9925_JPRB, 0.9926_JPRB, 0.9928_JPRB, 0.9927_JPRB, 0.9926_JPRB, 0.9926_JPRB, &
      0.9923_JPRB, 0.9920_JPRB, 0.9918_JPRB, 0.9916_JPRB, 0.9915_JPRB, 0.9913_JPRB, 0.9911_JPRB, 0.9907_JPRB, &
      0.9905_JPRB, 0.9903_JPRB, 0.9902_JPRB, 0.9900_JPRB, 0.9897_JPRB, 0.9896_JPRB, 0.9894_JPRB, 0.9892_JPRB, &
      0.9890_JPRB, 0.9889_JPRB, 0.9886_JPRB, 0.9884_JPRB, 0.9883_JPRB, 0.9884_JPRB, 0.9885_JPRB, 0.9885_JPRB, &
      0.9884_JPRB, 0.9883_JPRB, 0.9881_JPRB, 0.9880_JPRB, 0.9880_JPRB, 0.9880_JPRB, 0.9880_JPRB, 0.9879_JPRB, &
      0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9878_JPRB, 0.9877_JPRB, &
      0.9876_JPRB, 0.9876_JPRB, 0.9877_JPRB, 0.9876_JPRB, 0.9875_JPRB, 0.9874_JPRB, 0.9873_JPRB, 0.9873_JPRB, &
      0.9873_JPRB, 0.9874_JPRB, 0.9875_JPRB, 0.9875_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, &
      0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9873_JPRB, &
      0.9873_JPRB, 0.9873_JPRB, 0.9873_JPRB, 0.9873_JPRB, 0.9873_JPRB, 0.9872_JPRB, 0.9871_JPRB, 0.9872_JPRB, &
      0.9871_JPRB, 0.9870_JPRB, 0.9870_JPRB, 0.9870_JPRB, 0.9870_JPRB, 0.9869_JPRB, 0.9868_JPRB, 0.9868_JPRB, &
      0.9867_JPRB, 0.9866_JPRB, 0.9866_JPRB, 0.9865_JPRB, 0.9865_JPRB, 0.9865_JPRB, 0.9866_JPRB, 0.9866_JPRB, &
      0.9865_JPRB, 0.9865_JPRB, 0.9865_JPRB, 0.9866_JPRB, 0.9866_JPRB, 0.9866_JPRB, 0.9866_JPRB, 0.9867_JPRB, &
      0.9868_JPRB, 0.9868_JPRB, 0.9868_JPRB, 0.9867_JPRB, 0.9867_JPRB, 0.9867_JPRB, 0.9866_JPRB, 0.9867_JPRB, &
      0.9867_JPRB, 0.9867_JPRB, 0.9867_JPRB, 0.9867_JPRB, 0.9867_JPRB, 0.9868_JPRB, 0.9868_JPRB, 0.9868_JPRB, &
      0.9869_JPRB, 0.9869_JPRB, 0.9870_JPRB, 0.9872_JPRB, 0.9873_JPRB, 0.9874_JPRB, 0.9874_JPRB, 0.9874_JPRB, &
      0.9875_JPRB, 0.9875_JPRB, 0.9875_JPRB, 0.9875_JPRB, 0.9876_JPRB, 0.9876_JPRB, 0.9876_JPRB, 0.9876_JPRB, &
      0.9877_JPRB, 0.9877_JPRB, 0.9877_JPRB, 0.9877_JPRB, 0.9878_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, &
      0.9878_JPRB, 0.9878_JPRB, 0.9878_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9879_JPRB, 0.9878_JPRB, &
      0.9877_JPRB, 0.9876_JPRB, 0.9876_JPRB, 0.9877_JPRB, 0.9877_JPRB, 0.9876_JPRB, 0.9876_JPRB, 0.9876_JPRB, &
      0.9876_JPRB, 0.9876_JPRB, 0.9876_JPRB, 0.9875_JPRB, 0.9874_JPRB, 0.9873_JPRB, 0.9873_JPRB, 0.9872_JPRB, &
      0.9870_JPRB, 0.9869_JPRB, 0.9869_JPRB, 0.9869_JPRB, 0.9869_JPRB, 0.9869_JPRB, 0.9868_JPRB, 0.9867_JPRB, &
      0.9867_JPRB, 0.9866_JPRB, 0.9866_JPRB, 0.9865_JPRB, 0.9865_JPRB, 0.9865_JPRB, 0.9864_JPRB, 0.9863_JPRB, &
      0.9862_JPRB, 0.9862_JPRB, 0.9862_JPRB, 0.9862_JPRB, 0.9862_JPRB, 0.9861_JPRB, 0.9860_JPRB, 0.9860_JPRB, &
      0.9860_JPRB, 0.9859_JPRB, 0.9859_JPRB, 0.9859_JPRB, 0.9858_JPRB, 0.9858_JPRB, 0.9857_JPRB, 0.9857_JPRB, &
      0.9857_JPRB, 0.9856_JPRB, 0.9855_JPRB, 0.9855_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9853_JPRB, 0.9853_JPRB, &
      0.9852_JPRB, 0.9852_JPRB, 0.9852_JPRB, 0.9852_JPRB, 0.9852_JPRB, 0.9852_JPRB, 0.9851_JPRB, 0.9850_JPRB, &
      0.9850_JPRB, 0.9850_JPRB, 0.9850_JPRB, 0.9851_JPRB, 0.9850_JPRB, 0.9850_JPRB, 0.9849_JPRB, 0.9849_JPRB, &
      0.9850_JPRB, 0.9851_JPRB, 0.9851_JPRB, 0.9850_JPRB, 0.9850_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9849_JPRB, &
      0.9849_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9849_JPRB, &
      0.9849_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9850_JPRB, &
      0.9850_JPRB, 0.9850_JPRB, 0.9851_JPRB, 0.9851_JPRB, 0.9852_JPRB, 0.9852_JPRB, 0.9853_JPRB, 0.9853_JPRB, &
      0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, 0.9854_JPRB, &
      0.9855_JPRB, 0.9856_JPRB, 0.9856_JPRB, 0.9856_JPRB, 0.9856_JPRB, 0.9856_JPRB, 0.9856_JPRB, 0.9856_JPRB, &
      0.9856_JPRB, 0.9855_JPRB, 0.9855_JPRB, 0.9855_JPRB, 0.9854_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, &
      0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, 0.9853_JPRB, &
      0.9852_JPRB, 0.9851_JPRB, 0.9851_JPRB, 0.9850_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9849_JPRB, 0.9848_JPRB, &
      0.9848_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9848_JPRB, 0.9847_JPRB, 0.9846_JPRB, &
      0.9846_JPRB, 0.9846_JPRB, 0.9847_JPRB, 0.9846_JPRB, 0.9845_JPRB, 0.9844_JPRB, 0.9844_JPRB, 0.9843_JPRB, &
      0.9842_JPRB, 0.9842_JPRB, 0.9842_JPRB, 0.9842_JPRB, 0.9841_JPRB, 0.9841_JPRB, 0.9840_JPRB, 0.9839_JPRB, &
      0.9838_JPRB, 0.9838_JPRB, 0.9837_JPRB, 0.9837_JPRB, 0.9837_JPRB, 0.9836_JPRB, 0.9836_JPRB, 0.9835_JPRB, &
      0.9835_JPRB, 0.9834_JPRB, 0.9833_JPRB, 0.9832_JPRB, 0.9832_JPRB, 0.9832_JPRB, 0.9831_JPRB, 0.9831_JPRB, &
      0.9830_JPRB, 0.9829_JPRB, 0.9828_JPRB, 0.9828_JPRB, 0.9827_JPRB, 0.9827_JPRB, 0.9826_JPRB, 0.9826_JPRB, &
      0.9825_JPRB, 0.9824_JPRB, 0.9824_JPRB, 0.9823_JPRB, 0.9821_JPRB, 0.9821_JPRB, 0.9820_JPRB, 0.9821_JPRB, &
      0.9820_JPRB, 0.9820_JPRB, 0.9818_JPRB, 0.9817_JPRB, 0.9817_JPRB, 0.9816_JPRB, 0.9815_JPRB, 0.9815_JPRB, &
      0.9815_JPRB, 0.9814_JPRB, 0.9813_JPRB, 0.9813_JPRB, 0.9812_JPRB, 0.9811_JPRB, 0.9811_JPRB, 0.9810_JPRB/  

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------
  
  SUBROUTINE rttov_uwiremis_init( &
        &             path,       &! in
        &             imonth,     &! in
        &             err         )! out
          
    ! Description:
    ! initialize the rttov_uwiremis algorithm by (1) reading in the UW BF IR Global 
    ! Emissivty data and (2) the eigenvectors of the laboratory spectra, and (2) make some
    ! precalculations for the PCA regression
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9    03/31/2009   origianl code E. Borbas UW-Madison/CIMSS 
    !  1.0    03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
    
    Use rttov_const, Only :      &
            errorstatus_success, &
            errorstatus_fatal
                      
    Implicit None

    Character (len=*),  Intent(in)  :: path
    Integer(Kind=jpim), Intent(in)  :: imonth
    Integer(Kind=jpim), Intent(out) :: err

    
    Integer(Kind=jpim) :: i,info
    Integer(Kind=jpim) :: ipiv(numpcs)
    Integer(Kind=jpim), Parameter :: nmonth=12_jpim

    Real(Kind=jprb)  :: A(numpcs,hngpnts),AT(hngpnts,numpcs)
    Double Precision :: B(numpcs,numpcs),C(numpcs,numpcs)
    
    Character (len=300) :: fn
    Character (len=4) :: cyear
    Character (len=2) :: cmonth

    Logical(Kind=jplm) :: file_exists
    
    TRY
    
    err = errorstatus_success
    Write(cyear,'(i4)') db_ver_year
    Write(cmonth,'(i2.2)') imonth

    Allocate(bfemis_flag  (nb_lats,nb_lons))
    Allocate(bfemis_lut   (nb_lats,nb_lons))
    Allocate(bfemis       (nb_pack,hngpnts))
    
    Allocate(cov_emis_lut (cv_lats,cv_lons))
    Allocate(cov_emis     (cv_pack,numwave))
      
    !----------------------------------------------------------------------------
    ! reading the 0.1 degree resolution UW BF IR Land Surface Global Emissivty 
    !----------------------------------------------------------------------------
    fn=TRIM(path)//'UWirbfemis_V2.1_0.1deg_'//cyear//cmonth//'_mask.nc'
    Inquire(FILE=fn, EXIST=file_exists)
    If (.not. file_exists) Then
      err = errorstatus_fatal
      THROWM(err .ne. errorstatus_success, 'UWiremis means file not found : '//Trim(fn))
    End If
    
    INFO('Using UWiremis means fn=   '//Trim(fn))
    Call rttov_uwiremis_read_db(Trim(fn))
    
    !----------------------------------------------------------------------------
    ! reading the 0.5 degree resolution UW IR Land Surface Global Emissivty STD DEV
    !----------------------------------------------------------------------------
    fn=TRIM(path)//'UWiremis_hsremis_covmat_V1.0_deg0.5_month'//cmonth//'_mask.nc'
    Inquire(FILE=fn, EXIST=file_exists)
    If (.not. file_exists) Then
      err = errorstatus_fatal
      THROWM(err .ne. errorstatus_success, 'UWiremis covariances file not found : '//Trim(fn))
    End If

    INFO('Using UWiremis covariances fn=   '//Trim(fn))
    Call rttov_uwiremis_read_cov(Trim(fn))
          
    !-------------------------------------------------------------------
    !  reading the eigienvectors of the 128 selected laboratory spectra 
    !-------------------------------------------------------------------
    fn=TRIM(path)//'UWiremis_labeigvects.nc'
    Inquire(FILE=fn, EXIST=file_exists)
    If (.not. file_exists) Then
      err = errorstatus_fatal
      THROWM(err .ne. errorstatus_success, 'UWiremis file not found : '//Trim(fn))
    End If

    Call rttov_uwiremis_read_labeigvects(Trim(fn))
          
    !---------------------------------------------------
    ! compute D matrix for HSR emissivity construction
    !---------------------------------------------------
  
    ! define matrix A
    
    A=pcu_modres(1:numpcs,1:hngpnts)
    
    ! calculate transpose matrix A : AT [6X4]         
    
    AT=transpose(A) 
  
    !  B=A*A' [4X6]X[6X4]   first: row second: column
  
    B(:,:)=0.0_jprb
    B=MATMUL(A,AT)
  
    ! computing inverse of B
  
    C(:,:)=0.0_jprb
    Do i=1,numpcs
      C(i,i)=1.0_jprb
    Enddo
    
    Call DGESV( numpcs, numpcs, B, numpcs, IPIV, C, numpcs, info )
                  
    ! compute D=A'*inv(A*A')=A'*inv(B)= AT * C [6X4] * [4X4]
  
    D=MATMUL(AT,C)
    
    CATCH
  END SUBROUTINE rttov_uwiremis_init
  
  
  SUBROUTINE  rttov_uwiremis_read_cov( &
        &                          fn) !in
  
    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data   
    ! from the netCDF file into memory. 
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston 
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
              
    Implicit None
  
    ! Define Variables.
  
    Character(len=*), Intent(in) :: fn    ! filename including full path?
  

    Integer(Kind=ncint32) :: nvars     ! number of variables
    Integer(Kind=ncint32) :: ndims     ! number of dimensions
    Integer(Kind=ncint32) :: errstat   ! error code
    Integer(Kind=ncint32) :: recdim    ! record dimension
    Integer(Kind=ncint32) :: nc_dim(4)
  
    Integer(Kind=jpim)    :: i, j
    Integer(Kind=ncint32) :: ncid, ngatts, nrecs, varid
  
    Real(Kind=ncreal32)   ::  sfac
    Integer(Kind=ncint16), Allocatable, Dimension(:,:) :: emis_cov
    Integer(Kind=ncint8),  Allocatable, Dimension(:,:) :: pack_cov
    
    Character (len=1024) :: strbuf ! string buffer for var
    Integer(Kind=jpim)   :: indexlut
    
  
    ! Open netCDF file.
    errstat = nf_open(Trim(fn),nf_nowrite,ncid)
  
    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid,ndims,nvars,ngatts,recdim)
  
    Do recdim=1,ndims
      errstat = nf_inq_dim(ncid,recdim,strbuf,nrecs)
      nc_dim(recdim) = nrecs
    End Do
        
    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    Allocate(pack_cov(nc_dim(2),nc_dim(3)))
    errstat = nf_inq_varid (ncid,'mask',varid)
    errstat = nf_get_var_int1(ncid,varid,pack_cov)
  
    ! Generate the look-up table into the covariance data
    cov_emis_lut(:,:) = -1
    indexlut = 1
    Do i=1,nc_dim(3)
      Do j=1,nc_dim(2)
        If (pack_cov(j,i) > 0) then
          cov_emis_lut(j,i) = indexlut
          indexlut = indexlut + 1
        End If
      End Do
    End Do
  
    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    Allocate(emis_cov(nc_dim(1),nc_dim(4)))
    errstat = nf_inq_varid (ncid,'emis_diagCov',varid)
    errstat = nf_get_var_int2(ncid,varid,emis_cov)
    errstat = nf_get_att_real(ncid,varid,'scale_factor',sfac)
  
    Do i=1,nc_dim(1)
      cov_emis(:,i) = Real(emis_cov(i,:),Kind=jprb)*sfac 
    End Do
  
    Deallocate(pack_cov,emis_cov)
          
    errstat = nf_close(ncid)

  END SUBROUTINE rttov_uwiremis_read_cov
  
  
  SUBROUTINE  rttov_uwiremis_read_db( &
        &                          fn) !in
  
    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data   
    ! from the netCDF file into memory. 
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston 
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
  
    Implicit None
  
    ! Define Variables.
  
    Character(len=*), Intent(in) :: fn    ! filename including full path?


    Integer(Kind=ncint32) :: nvars     ! number of variables
    Integer(Kind=ncint32) :: ndims     ! number of dimensions
    Integer(Kind=ncint32) :: errstat   ! error code
    Integer(Kind=ncint32) :: recdim    ! record dimension
    Integer(Kind=ncint32) :: nc_dim(4) ! hng_pnt, lats, lons, pack_len
  
    Integer(Kind=jpim)    :: i,j,k
    Integer(Kind=ncint32) :: ncid, ngatts, nrecs, varid
  
    Real(Kind=ncreal32)   :: sfac, offs
    Integer(Kind=ncint16), Allocatable, Dimension(:) :: emis_ch
  
    Character (len=1024) :: strbuf ! string buffer for var
    Character (len=6)    :: cfld 
    Integer(Kind=jpim)   :: indexlut

    ! Open netCDF file.
    errstat = nf_open(Trim(fn),nf_nowrite,ncid)
    
    ! Get info on the record dimension for this file. 
    errstat = nf_inq(ncid,ndims,nvars,ngatts,recdim)
    
    Do recdim=1,ndims
      errstat = nf_inq_dim(ncid,recdim,strbuf,nrecs)
      nc_dim(recdim) = nrecs
    End Do
  
    Allocate(emis_ch(nc_dim(4)))
  
    ! Retrieve emissivity database flag value
    errstat = nf_inq_varid (ncid, 'emis_flag', varid)
    errstat = nf_get_var_int2(ncid,varid,bfemis_flag)
  
    ! Generate the look-up table into the emissivity data
    bfemis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    Do i=1,nc_dim(3)
      Do j=1,nc_dim(2)
        If (bfemis_flag(j,i) > 0) then
          bfemis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        End If
      End Do
    End Do
  
    ! Retrieve database of 10 MODIS hinge-point emissivity values
    Do k=1,nc_dim(1)
      If (k < 10) Then
        write(cfld,'("emis",i1," ")') k
      Else
        write(cfld,'("emis",i2.2)') k
      End if
      errstat = nf_inq_varid (ncid,cfld,varid)
      errstat = nf_get_var_int2(ncid,varid,emis_ch)
      errstat = nf_get_att_real(ncid,varid,'scale_factor',sfac)
      errstat = nf_get_att_real(ncid,varid,'add_offset',offs)

      bfemis(:,11-k) = Real(emis_ch,Kind=jprb)*sfac + offs
    End Do
  
    Deallocate(emis_ch)

    errstat = nf_close(ncid)
  
  END SUBROUTINE rttov_uwiremis_read_db
  
  
  SUBROUTINE  rttov_uwiremis_read_labeigvects( &
        &                                   fn) !in
  
    ! Description:
    ! read the eigienvectors of the 128 selected laboratory spectra 
    !(created by E borbas) from the netCDF file into memory. 
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston 
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
  
    Implicit None
  
    ! Define Variables.
    Character(len=*), Intent(in) :: fn    ! filename including full path?
  

    Integer(Kind=ncint32) :: nvars     ! number of variables
    Integer(Kind=ncint32) :: ndims     ! parameter (ndims = 1) ! number of dimensions
    Integer(Kind=ncint32) :: errstat   ! error code
    Integer(Kind=ncint32) :: recdim    ! record dimension
    Integer(Kind=ncint32) :: xy_dim(2)
  
    Integer(Kind=ncint32) :: ncid, ngatts, nrecs, varid
  
    Character (len=1024)  :: strbuf  ! string buffer for var
  
  
    ! Open netCDF file.
    errstat = nf_open(Trim(fn),nf_nowrite,ncid)
  
    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid,ndims,nvars,ngatts,recdim)
    
    Do recdim=1,ndims
      errstat = nf_inq_dim(ncid,recdim,strbuf,nrecs)
      xy_dim(recdim) = nrecs
    End Do
  
    ! Retrieve the laboratory eigenvalues into array pcu 
    errstat = nf_inq_varid (ncid, 'PC_scores', varid)
    errstat = nf_get_var_real(ncid,varid,pcu)

    errstat = nf_close(ncid)
  
  END SUBROUTINE rttov_uwiremis_read_labeigvects
  
  
   
!------------------------------------------
! Routines for returning emissivity values
!------------------------------------------
  
  SUBROUTINE rttov_uwiremis( &
        & nchs,              &! in 
        & lat,               &! in 
        & lon,               &! in
        & surfacetype,       &! in
        & snowfrac,          &! in 
        & instr_wavenum,     &! in 
        & instr_emis,        &! out 
        & instr_emis_cov,    &! out 
        & instr_emis_flag)    ! out 
  
    ! Description:
    ! To compute IR emissivty for a given location and frequency 
    ! from the 0.1 degree resolution UW BF IR Global Emissivty data 
    ! (at 10 hinge points) (http://cimss.ssec.wisc.edu/iremis/)
    ! and labratory measurements using principal component analyses
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
  
    Use rttov_const, Only : surftype_land, surftype_seaice
  
    Implicit None
  
    Integer(Kind=jpim), intent(in)  :: nchs
    Integer(Kind=jpim), intent(in)  :: surfacetype
    Integer(Kind=jpim), Intent(out) :: instr_emis_flag
  
    Real(Kind=jprb),    Intent(in)  :: lat, lon
    Real(Kind=jprb),    Intent(in)  :: snowfrac
    Real(Kind=jprb),    Intent(in)  :: instr_wavenum(nchs)
    Real(Kind=jprb),    Intent(out) :: instr_emis(nchs)
    Real(Kind=jprb),    Intent(out) :: instr_emis_cov(nchs)
  

    Real(Kind=jprb) :: bfem(hngpnts)
    Real(Kind=jprb) :: hsremis(numwave)
    Real(Kind=jprb) :: emis_cov(numwave)
    
    Real(Kind=jprb) :: long
    Integer(Kind=jpim) :: gridy, gridx, rnd_x, rnd_y
          
    Real(Kind=jprb) :: hkod = -999._JPRB
    
    Real(Kind=jprb) :: hsr_ir_emis(nchs)
    Real(Kind=jprb) :: hsr_ir_emis_cov(nchs)
    Real(Kind=jprb) :: cov_buff(numwave)
    
    Character(len=128) :: msg
    
    instr_emis(:) = hkod  
        
    If ( surfacetype == surftype_land ) Then
        
    !------------------------------------------------------------
    ! find the closest grid point from the uwiremis database
    !------------------------------------------------------------
      
      long = modulo(lon,360.0_jprb)
      If (long > 180.0) Then
        long = long - 360.0_jprb
      End If
    
      gridy = NInt(Abs(bfemis_ygrid1-lat)/bfemis_gridres,Kind=jpim)+1_jpim
      gridx = NInt(Abs(bfemis_xgrid1-long)/bfemis_gridres,Kind=jpim)+1_jpim
                                      
      rnd_y = NInt(Abs(cov_emis_ygrid1-lat)/cov_emis_gridres,Kind=jpim)+1_jpim
      rnd_x = NInt(Abs(cov_emis_xgrid1-long)/cov_emis_gridres,Kind=jpim)+1_jpim
                                      
      instr_emis_flag = bfemis_flag(gridy,gridx)

    !------------------------------
    ! check if it is a land pixel
    !------------------------------ 
            
      If ( instr_emis_flag > 0 ) Then
        
        ! Find the emissivity and covariances
        
        If (bfemis_lut(gridy,gridx) > 0) Then
          bfem(:) = Real(bfemis(bfemis_lut(gridy,gridx),:),Kind=JPRB)
        Else
          bfem(:) = 0.0_jprb
        End If
        If (cov_emis_lut(rnd_y,rnd_x) > 0) Then
          emis_cov(:) = SQRT(Real(cov_emis(cov_emis_lut(rnd_y,rnd_x),:),Kind=JPRB))
        Else
          emis_cov(:) = default_std
        End If
          
    !--------------------------------------------------------------------------------------------------------
    ! compute the hsr emissivity spectra at 416 wavenumber points from the 10 BF emissivity hinge points
    !-------------------------------------------------------------------------------------------------------
    
        Call rttov_uwiremis_recon_hsremis( &
          & bfem,                          &! in
          & hsremis)                        ! out
      
    !--------------------------------------------------------------------------------
    ! create instrument specific emis/stdv by finidng the closest wavenumber value
    !--------------------------------------------------------------------------------
    
        Call rttov_uwiremis_select_wavenum( &
          & hsremis,                        &! in  
          & emis_cov,                       &! in  
          & nchs,                           &! in 
          & instr_wavenum(1:nchs),          &! in 
          & instr_emis,                     &! out 
          & instr_emis_cov       )           ! out 
            
    !---------------------------------------------------
    ! Linearly blend avg snow emis with snow cover frac
    !---------------------------------------------------
        If ( snowfrac > 0.0 ) Then

          ! snow_stdv is a fixed stdv
          cov_buff(:) = snow_stdv
          Call rttov_uwiremis_select_wavenum( &
                  & snow_em,                  &! in  
                  & cov_buff,                 &! in  
                  & nchs,                     &! in 
                  & instr_wavenum,            &! in 
                  & hsr_ir_emis,              &! out 
                  & hsr_ir_emis_cov       )    ! out 
      
          If ( snowfrac > 1.0 ) Then
            instr_emis(:) = hsr_ir_emis(:)
            ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
            instr_emis_cov(:) = hsr_ir_emis_cov(:)
          Else
            instr_emis(:) = snowfrac*hsr_ir_emis(:) + (1.0_jprb-snowfrac)*instr_emis(:)
            ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
            instr_emis_cov(:) = snowfrac*hsr_ir_emis_cov(:) + (1.0_jprb-snowfrac)*instr_emis_cov(:)
          End If
          
        End If  ! snow chk
        
      End If  ! emis flag > 0
    
    Else If ( surfacetype == surftype_seaice) Then
    
    !---------------------------------------
    ! Return emissivity and stdv for seaice
    !---------------------------------------
  
      ! sice_stdv is a fixed stdv
      cov_buff(:) = sice_stdv
      Call rttov_uwiremis_select_wavenum( &
                & sice_em,                &! in  
                & cov_buff,               &! in  
                & nchs,                   &! in 
                & instr_wavenum,          &! out 
                & instr_emis,             &! out 
                & instr_emis_cov        )  ! out 
        
      instr_emis_flag = seaice_flag
    Else
      Write (msg,'(a)') 'Warning: IR emissivity atlas should only be called for land and seaice surface types'
      Call rttov_errorreport(errorstatus_success, msg, '')
    End If
    
  END SUBROUTINE rttov_uwiremis
  
  
  SUBROUTINE rttov_uwiremis_recon_hsremis( &
        & bfem,                            &! in 
        & hsremis)                          ! out 
          
    ! Description:
    ! To creates high spectra resolution emissivties at 416 wavenumbers
    ! from the UW BF IR Global Emissivty data (at 10 hinge points) 
    ! and labratory measurements using principal component analyses
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
          
  
    Implicit None
  
    
    Real(Kind=jprb), Intent(in)  :: bfem(hngpnts)
    Real(Kind=jprb), Intent(out) :: hsremis(numwave)
  

    Integer(Kind=jpim) :: j,k
  
    Real(Kind=jprb) :: col(hngpnts)     
    Real(Kind=jprb) :: coef(numpcs) 
    Real(Kind=jprb) :: emis(numwave) 
  
    !-------------------------
    ! calculate the regcoef
    !-------------------------
  
    coef(:) = 0._JPRB
      
    If ( bfem(1) > 0._JPRB ) Then
      
      Do j=1,hngpnts
        col(j) = bfem(j) - pcm_modres(j)
      End Do
    
      Do k=1,numpcs
        Do j=1,hngpnts
          coef(k) = coef(k) + col(j) * D(j,k)
        End Do
      End Do

    !-----------------------------------
    ! apply regcoef to get the hsr dataset
    !-----------------------------------
  
      emis=0._JPRB
          
      If ( coef(1) .ne. -999._JPRB ) Then
  
        Do k=1,numwave
          Do j=1,numpcs
            emis(k) = emis(k) + coef(j) * Real(pcu(j,k),Kind=JPRB)
          End Do
  
          emis(k) = emis(k) + pcm(k)
        End Do
  
    !--------------------------------------------
    !  filter out emissivty values larger then 1.
    !--------------------------------------------
  
        Do j=1,numwave
          If ( emis(j) > 1._JPRB ) emis(j) = 1._JPRB                        
          hsremis(j) = emis(j)
        End Do
          
      Else
                  
        coef = -999._JPRB
        hsremis = -999._JPRB
              
      End If

    End If
                  
  END SUBROUTINE rttov_uwiremis_recon_hsremis
  
  
  SUBROUTINE rttov_uwiremis_select_wavenum ( &
        & hsremis,                           &! in
        & emis_cov,                          &! in
        & nchs,                              &! in
        & instr_wavenum,                     &! in
        & instr_emis,                        &! out
        & instr_emis_cov)                     ! out
  
    ! Description:
    ! subroutine to find the closest wavenumber from the UW HSR emissivity spectra 
    ! for the instrument frequency and assign the instrument emissivity by choosing the 
    ! closest spectral point value or bilinear interpolating  between the two 
    ! closest spectral point values
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
              
    Integer(Kind=jpim), Intent(in) :: nchs
  
    Real(Kind=jprb), Intent(in) :: hsremis(numwave)
    Real(Kind=jprb), Intent(in) :: emis_cov(numwave)
    Real(Kind=jprb), Intent(in) :: instr_wavenum(nchs)
  
    Real(Kind=jprb), Intent(out) :: instr_emis(nchs)
    Real(Kind=jprb), Intent(out) :: instr_emis_cov(nchs)
  

    Integer(Kind=jpim) :: j,k
  
    Real(Kind=jprb) :: dist(numwave)
    Real(Kind=jprb) :: mindist(nchs)
    Integer(Kind=jpim) :: ind_mindist(nchs)
          
    Real(Kind=jprb) :: dwvnum1,dwvnum2,dwvsum
    Real(Kind=jprb) :: hsremis1,hsremis2,emis_cov1,emis_cov2
    Logical(Kind=jplm) :: lcpu_emis, lcpu_cov
  
    Real(Kind=jprb) :: hkod = -999._JPRB
  
    ! initialize instrument emissivity
    
    !---------------------------------------------------------------
    ! finindg the closest frequency from the hsr emissivity spectr
    !--------------------------------------------------------------
    lcpu_emis = .True.
    lcpu_cov  = .True.
    If ( All(hsremis == hsremis(1)) ) lcpu_emis = .False.
    If ( All(emis_cov == emis_cov(1)) ) lcpu_cov  = .False.
              
    If (lcpu_emis .or. lcpu_cov) Then           
      instr_emis(:)       = hkod
      instr_emis_cov(:)   = hkod
      Do j=1,nchs
      
        If(instr_wavenum(j) <= hsr_wavenum(1)) Then

          instr_emis(j)=hsremis(1)
          instr_emis_cov(j)=emis_cov(1)

        Elseif(instr_wavenum(j) >= hsr_wavenum(numwave)) Then     

          instr_emis(j)=hsremis(numwave)
          instr_emis_cov(j)=emis_cov(numwave)

        Else ! within wavenumber compute range

          mindist(j) = 100._JPRB
          ind_mindist(j) = 100000_jpim
                                                  
          Do k=1,numwave
  
            ! calucalte distances between the instr freq end hsr emissivities                             
            dist(k) = abs( instr_wavenum(j) - hsr_wavenum(k) )
        
            ! finding the closest frequency from the hsr emissivity 
            If( dist(k) <=  mindist(j) ) Then
        
              mindist(j) = dist(k)
              ind_mindist(j) = k
        
            End If
          End do
  
    !--------------------------------
    ! assign instrument emissivity
    !-------------------------------- 
    !  closest spectral point
    !                       instr_emis(j)=hsremis(ind_mindist(j))
    !                       instr_emis_cov(j)=cov_emis(ind_mindist(j))
    
    ! or bilinear mean of the two closest spectral points
  
          k=1 
          dwvnum1=0._jprb
          dwvnum2=0._jprb
          hsremis1=0._jprb
          hsremis2=0._jprb
          emis_cov1=0._jprb
          emis_cov2=0._jprb

          If ( instr_wavenum(j) <= hsr_wavenum(ind_mindist(j)) ) k=-1

            dwvnum1 = dist( ind_mindist(j) )
            dwvnum2 = dist( ind_mindist(j) + k )
            dwvsum = dwvnum1 + dwvnum2

          If ( lcpu_emis ) Then
            hsremis1 = dwvnum1 * hsremis( ind_mindist(j) + k )
            hsremis2 = dwvnum2 * hsremis( ind_mindist(j) )
            instr_emis(j) = ( hsremis1 + hsremis2 ) / dwvsum
          Else
            instr_emis(j) = hsremis(1)
          End If

          If ( lcpu_cov ) Then
            emis_cov1 = dwvnum1 * emis_cov( ind_mindist(j) + k )
            emis_cov2 = dwvnum2 * emis_cov( ind_mindist(j) )
            instr_emis_cov(j) = ( emis_cov1 + emis_cov2 ) / dwvsum
          Else
            instr_emis_cov(j) = emis_cov(1)
          End If
  
        End If    !==  (instr_wavenum(j) <= hsr_wavenum(1)) 
  
      End Do
  
    Else  ! all logical computes (lcpu_emis, lcpu_cov) are false
      instr_emis(:)=hsremis(1)
      instr_emis_cov(:)=emis_cov(1)
    End If
  
  END SUBROUTINE rttov_uwiremis_select_wavenum
  
  
!------------------------------------
! Routine to deallocate atlas arrays
!------------------------------------
  SUBROUTINE rttov_uwiremis_close_atlas
    IF ( ALLOCATED(bfemis_flag)  ) DEALLOCATE(bfemis_flag)
    IF ( ALLOCATED(bfemis_lut)   ) DEALLOCATE(bfemis_lut)
    IF ( ALLOCATED(bfemis)       ) DEALLOCATE(bfemis)
    IF ( ALLOCATED(cov_emis_lut) ) DEALLOCATE(cov_emis_lut)
    IF ( ALLOCATED(cov_emis)     ) DEALLOCATE(cov_emis)
  END SUBROUTINE
  
END MODULE mod_iratlas
