! $Id: ropp_io_free.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Initialisation/ropp_io_free *
!
! NAME
!    ropp_io_free - Free arrays within an an RO derived type
!
! SYNOPSIS
!    use ropp_io
!    type(ROprof) :: ro_data
!      ...
!    call ropp_io_free(ro_data)
!
!      - or -
!
!    call ropp_io_free(ro_data%Lev1a)
!    call ropp_io_free(ro_data%Lev1b)
!    call ropp_io_free(ro_data%Lev2a)
!    call ropp_io_free(ro_data%Lev2b)
!    call ropp_io_free(ro_data%Lev2c)
!    call ropp_io_free(ro_data%Lev2d)
!
! DESCRIPTION
!   This subroutine frees memory from a previously initialised RO data structure
!   or parts thereof; scalar values are set to their 'missing' data status again.
!
! OUTPUT
!    ro_data   dtyp  RO data (derived type)
!
! SEE ALSO
!    ropp_io_types
!
! NOTES
!    Calling the deallocation / freeing routines for individual sublevels of
!    the input data structure is NOT equivalent to freeing the same data
!    structure in one go, as the header data is only reset to its default
!    (missing) values in the second case.
!
!    Meta data information like units and ranges are NOT changed back to their
!    default values, but left with the definition as set by the user.
!
!    The current interface only allows the deallocation of scalar structures.
!    To free an array of structures (e.g., as obtained from reading a multifile),
!    a loop over all elements of the array of structures is required:
!
!       do i = 1, size(ro_data)
!          call ropp_io_free(ro_data(i))
!       enddo
!
! REFERENCES
!    Format Definition for Radio Occultation Files -
!    CLIMAP Format Version 2.2a, Issue 1.6, 8 January 2004
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. L1atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l1atype(var)

! 1.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_free_l1atype
  USE ropp_io_types, ONLY: L1atype

  IMPLICIT NONE

  TYPE(L1atype), INTENT(inout) :: var

! 1.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%dtime))      DEALLOCATE(var%dtime)
  IF (ASSOCIATED(var%snr_L1ca))   DEALLOCATE(var%snr_L1ca)
  IF (ASSOCIATED(var%snr_L1p))    DEALLOCATE(var%snr_L1p)
  IF (ASSOCIATED(var%snr_L2p))    DEALLOCATE(var%snr_L2p)
  IF (ASSOCIATED(var%phase_L1))   DEALLOCATE(var%phase_L1)
  IF (ASSOCIATED(var%phase_L2))   DEALLOCATE(var%phase_L2)
  IF (ASSOCIATED(var%phase_qual)) DEALLOCATE(var%phase_qual)

  IF (ASSOCIATED(var%r_gns))      DEALLOCATE(var%r_gns)
  IF (ASSOCIATED(var%v_gns))      DEALLOCATE(var%v_gns)
  IF (ASSOCIATED(var%r_leo))      DEALLOCATE(var%r_leo)
  IF (ASSOCIATED(var%v_leo))      DEALLOCATE(var%v_leo)

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l1atype


!-------------------------------------------------------------------------------
! 2. L1btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l1btype(var)

! 2.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L1btype

  IMPLICIT NONE

  TYPE(L1btype), INTENT(inout) :: var

! 2.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%lat_tp))           DEALLOCATE(var%lat_tp)
  IF (ASSOCIATED(var%lon_tp))           DEALLOCATE(var%lon_tp)
  IF (ASSOCIATED(var%azimuth_tp))       DEALLOCATE(var%azimuth_tp)
  IF (ASSOCIATED(var%impact_L1))        DEALLOCATE(var%impact_L1)
  IF (ASSOCIATED(var%impact_L2))        DEALLOCATE(var%impact_L2)
  IF (ASSOCIATED(var%impact))           DEALLOCATE(var%impact)
  IF (ASSOCIATED(var%impact_opt))       DEALLOCATE(var%impact_opt)
  IF (ASSOCIATED(var%bangle_L1))        DEALLOCATE(var%bangle_L1)
  IF (ASSOCIATED(var%bangle_L2))        DEALLOCATE(var%bangle_L2)
  IF (ASSOCIATED(var%bangle))           DEALLOCATE(var%bangle)
  IF (ASSOCIATED(var%bangle_opt))       DEALLOCATE(var%bangle_opt)
  IF (ASSOCIATED(var%bangle_L1_sigma))  DEALLOCATE(var%bangle_L1_sigma)
  IF (ASSOCIATED(var%bangle_L2_sigma))  DEALLOCATE(var%bangle_L2_sigma)
  IF (ASSOCIATED(var%bangle_sigma))     DEALLOCATE(var%bangle_sigma)
  IF (ASSOCIATED(var%bangle_opt_sigma)) DEALLOCATE(var%bangle_opt_sigma)
  IF (ASSOCIATED(var%bangle_L1_qual))   DEALLOCATE(var%bangle_L1_qual)
  IF (ASSOCIATED(var%bangle_L2_qual))   DEALLOCATE(var%bangle_L2_qual)
  IF (ASSOCIATED(var%bangle_qual))      DEALLOCATE(var%bangle_qual)
  IF (ASSOCIATED(var%bangle_opt_qual))  DEALLOCATE(var%bangle_opt_qual)

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l1btype


!-------------------------------------------------------------------------------
! 3. L2atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2atype(var)

! 3.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2atype

  IMPLICIT NONE

  TYPE(L2atype), INTENT(inout) :: var

! 3.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%alt_refrac))     DEALLOCATE(var%alt_refrac)
  IF (ASSOCIATED(var%geop_refrac))    DEALLOCATE(var%geop_refrac)
  IF (ASSOCIATED(var%refrac))         DEALLOCATE(var%refrac)
  IF (ASSOCIATED(var%refrac_sigma))   DEALLOCATE(var%refrac_sigma)
  IF (ASSOCIATED(var%refrac_qual))    DEALLOCATE(var%refrac_qual)
  IF (ASSOCIATED(var%dry_temp))       DEALLOCATE(var%dry_temp)
  IF (ASSOCIATED(var%dry_temp_sigma)) DEALLOCATE(var%dry_temp_sigma)
  IF (ASSOCIATED(var%dry_temp_qual))  DEALLOCATE(var%dry_temp_qual)

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l2atype


!-------------------------------------------------------------------------------
! 4a. L2btype (1d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2btype(var)

! 4.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2btype

  IMPLICIT NONE

  TYPE(L2btype), INTENT(inout) :: var

! 4.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%geop))        DEALLOCATE(var%geop)
  IF (ASSOCIATED(var%geop_sigma))  DEALLOCATE(var%geop_sigma)
  IF (ASSOCIATED(var%press))       DEALLOCATE(var%press)
  IF (ASSOCIATED(var%press_sigma)) DEALLOCATE(var%press_sigma)
  IF (ASSOCIATED(var%temp))        DEALLOCATE(var%temp)
  IF (ASSOCIATED(var%temp_sigma))  DEALLOCATE(var%temp_sigma)
  IF (ASSOCIATED(var%shum))        DEALLOCATE(var%shum)
  IF (ASSOCIATED(var%shum_sigma))  DEALLOCATE(var%shum_sigma)
  IF (ASSOCIATED(var%meteo_qual))  DEALLOCATE(var%meteo_qual)

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l2btype

!-------------------------------------------------------------------------------
! 4b. L2btype (2d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2btype_2d(var)

! 4.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2btype_2d

  IMPLICIT NONE

  TYPE(L2btype_2d), INTENT(inout) :: var

! 4.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%geop))        DEALLOCATE(var%geop)
  IF (ASSOCIATED(var%geop_sigma))  DEALLOCATE(var%geop_sigma)
  IF (ASSOCIATED(var%press))       DEALLOCATE(var%press)
  IF (ASSOCIATED(var%press_sigma)) DEALLOCATE(var%press_sigma)
  IF (ASSOCIATED(var%temp))        DEALLOCATE(var%temp)
  IF (ASSOCIATED(var%temp_sigma))  DEALLOCATE(var%temp_sigma)
  IF (ASSOCIATED(var%shum))        DEALLOCATE(var%shum)
  IF (ASSOCIATED(var%shum_sigma))  DEALLOCATE(var%shum_sigma)
  IF (ASSOCIATED(var%meteo_qual))  DEALLOCATE(var%meteo_qual)

  var%Npoints = 0
  var%Nhoriz = 0

END SUBROUTINE ropp_io_free_l2btype_2d

!-------------------------------------------------------------------------------
! 5a. L2ctype (1d version)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2ctype(var)

! 5.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2ctype

  IMPLICIT NONE

  TYPE(L2ctype), INTENT(inout) :: var

! 5.2 Reset scalar values for all structure elements
! --------------------------------------------------

  var%geop_sfc        = ropp_MDFV
  var%press_sfc       = ropp_MDFV
  var%press_sfc_sigma = ropp_MDFV
  var%press_sfc_qual  = ropp_MDFV

  var%Ne_max          = ropp_MDFV
  var%Ne_max_sigma    = ropp_MDFV
  var%H_peak          = ropp_MDFV
  var%H_peak_sigma    = ropp_MDFV
  var%H_width         = ropp_MDFV
  var%H_width_sigma   = ropp_MDFV
  var%direct_ion      = .FALSE.

  var%tph_bangle      = ropp_MDFV
  var%tpa_bangle      = ropp_MDFV
  var%tph_bangle_flag = ropp_MIFV

  var%tph_refrac      = ropp_MDFV
  var%tpn_refrac      = ropp_MDFV
  var%tph_refrac_flag = ropp_MIFV

  var%tph_tdry_lrt        = ropp_MDFV
  var%tpt_tdry_lrt        = ropp_MDFV
  var%tph_tdry_lrt_flag   = ropp_MIFV

  var%tph_tdry_cpt        = ropp_MDFV
  var%tpt_tdry_cpt        = ropp_MDFV
  var%tph_tdry_cpt_flag   = ropp_MIFV

  var%prh_tdry_cpt        = ropp_MDFV
  var%prt_tdry_cpt        = ropp_MDFV
  var%prh_tdry_cpt_flag   = ropp_MIFV

  var%tph_temp_lrt        = ropp_MDFV
  var%tpt_temp_lrt        = ropp_MDFV
  var%tph_temp_lrt_flag   = ropp_MIFV

  var%tph_temp_cpt        = ropp_MDFV
  var%tpt_temp_cpt        = ropp_MDFV
  var%tph_temp_cpt_flag   = ropp_MIFV

  var%prh_temp_cpt        = ropp_MDFV
  var%prt_temp_cpt        = ropp_MDFV
  var%prh_temp_cpt_flag   = ropp_MIFV

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l2ctype

!------------------------------------------------------------------------------
! 5b. L2ctype (2d version)
!------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2ctype_2d(var)

! 5b.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2ctype_2d

  IMPLICIT NONE

  TYPE(L2ctype_2d), INTENT(inout) :: var

! 5b.2 Reset scalar values for all structure elements
! --------------------------------------------------

  IF (ASSOCIATED(var%lat_2d))          DEALLOCATE(var%lat_2d)
  IF (ASSOCIATED(var%lon_2d))          DEALLOCATE(var%lon_2d)
  IF (ASSOCIATED(var%geop_sfc))        DEALLOCATE(var%geop_sfc)
  IF (ASSOCIATED(var%press_sfc))       DEALLOCATE(var%press_sfc)
  IF (ASSOCIATED(var%press_sfc_sigma)) DEALLOCATE(var%press_sfc_sigma)
  IF (ASSOCIATED(var%press_sfc_qual))  DEALLOCATE(var%press_sfc_qual)

  var%Npoints = 0
  var%Nhoriz  = 0

END SUBROUTINE ropp_io_free_l2ctype_2d

!-------------------------------------------------------------------------------
! 6. L2dtype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_l2dtype(var)

! 6.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io_types, ONLY: L2dtype

  IMPLICIT NONE

  TYPE(L2dtype), INTENT(inout) :: var

! 6.2 Deallocate memory for all structure elements
! ------------------------------------------------

  IF (ASSOCIATED(var%level_coeff_a)) DEALLOCATE(var%level_coeff_a)
  IF (ASSOCIATED(var%level_coeff_b)) DEALLOCATE(var%level_coeff_b)

  var%Npoints = 0

END SUBROUTINE ropp_io_free_l2dtype


!-------------------------------------------------------------------------------
! 7. Joint RO data type
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_rotype(ROdata)

! 7.1 Declarations
! ----------------

  USE ropp_utils,     ONLY: ropp_MDFV
! USE ropp_io,        not_this => ropp_io_free_rotype
  USE ropp_io_types,  ONLY: ROprof,       &
                            PCD_missing

  IMPLICIT NONE

  TYPE(ROprof), INTENT(inout) :: ROdata

! 7.2 Reset scalar values
! -----------------------

  ROdata%FmtVersion        = "UNKNOWN" ! File format version ID
  ROdata%occ_id            = "UNKNOWN" ! Occultation ID
  ROdata%leo_id            = "UNKN"    ! LEO identifier
  ROdata%gns_id            = "UNKN"    ! GNSS identifier
  ROdata%stn_id            = "UNKN"    ! GSN station identifier
  ROdata%processing_centre = "UNKNOWN" ! Processing centre
  ROdata%processing_software = "UNKNOWN" ! Processing centre software
  ROdata%pod_method        = "UNKNOWN" ! POD processing method
  ROdata%bangle_method     = "UNKNOWN" ! Bending angle processing method
  ROdata%refrac_method     = "UNKNOWN" ! Refractivity processing method
  ROdata%meteo_method      = "UNKNOWN" ! Meteorological processing method
  ROdata%thin_method       = "UNKNOWN" ! Profile thinning method
  ROdata%software_version  = "UNKNOWN" ! Software version ID

  ROdata%DTpro%Year   = 9999
  ROdata%DTpro%Month  =   99
  ROdata%DTpro%Day    =   99
  ROdata%DTpro%Hour   =   99
  ROdata%DTpro%Minute =   99
  ROdata%DTpro%Second =   99
  ROdata%DTpro%Msec   =  999

  ROdata%DTocc%Year   = 9999
  ROdata%DTocc%Month  =   99
  ROdata%DTocc%Day    =   99
  ROdata%DTocc%Hour   =   99
  ROdata%DTocc%Minute =   99
  ROdata%DTocc%Second =   99
  ROdata%DTocc%Msec   =  999

  ROdata%PCD          = 0
  ROdata%PCD          = IBSET ( ROdata%PCD, PCD_missing )

  ROdata%overall_qual = ropp_MDFV   ! Overall quality value

  ROdata%georef%time_offset =  ropp_MDFV ! Time since start (s)
  ROdata%georef%lat         =  ropp_MDFV ! Latitude  (deg)
  ROdata%georef%lon         =  ropp_MDFV ! Longitude (deg)
  ROdata%georef%roc         =  ropp_MDFV ! RoC value (m)
  ROdata%georef%r_coc       =  ropp_MDFV ! RoC offset X,Y,Z vector (m)
  ROdata%georef%azimuth     =  ropp_MDFV ! GNSS->LEO line of sight angle (degT)
  ROdata%georef%undulation  =  ropp_MDFV ! Geoid undulation (EGM96-WGS84) (m)

  ROdata%bg%source   = 'NONE' ! Source of b/g profile
  ROdata%bg%year     = 9999   ! VT year   of b/g
  ROdata%bg%month    =   99   ! VT month  of b/g
  ROdata%bg%day      =   99   ! VT day    of b/g
  ROdata%bg%hour     =   99   ! VT hour   of b/g
  ROdata%bg%minute   =   99   ! VT minute of b/g
  ROdata%bg%fcperiod =  999   ! F/c period (hrs)

! 7.3 Level 1a profile
! --------------------

  CALL ropp_io_free(ROdata%Lev1a)

! 7.4 Level 1b profile
! --------------------

  CALL ropp_io_free(ROdata%Lev1b)

! 7.5 Level 2a profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2a)

! 7.6 Level 2b profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2b)

! 7.7 Level 2c profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2c)

! 7.8 Level 2d profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2d)

! 7.9 Additional variables
! ------------------------

  CALL ropp_io_free(ROdata%vlist)

END SUBROUTINE ropp_io_free_rotype

!-------------------------------------------------------------------------------
! 8. Joint RO data type (two-dimensional meteorological data)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_rotype_2d(ROdata)

! 8.1 Declarations
! ----------------

! USE ropp_io,        not_this => ropp_io_free_rotype_2d
  USE ropp_io_types,  ONLY: ROprof2d,     &
                            ropp_MDFV,    &
                            PCD_missing

  IMPLICIT NONE

  TYPE(ROprof2d), INTENT(inout) :: ROdata

! 8.2 Reset scalar values
! -----------------------

  ROdata%FmtVersion        = "UNKNOWN" ! File format version ID
  ROdata%occ_id            = "UNKNOWN" ! Occultation ID
  ROdata%leo_id            = "UNKN"    ! LEO identifier
  ROdata%gns_id            = "UNKN"    ! GNSS identifier
  ROdata%stn_id            = "UNKN"    ! GSN station identifier
  ROdata%processing_centre = "UNKNOWN" ! Processing centre
  ROdata%processing_software = "UNKNOWN" ! Processing centre software
  ROdata%pod_method        = "UNKNOWN" ! POD processing method
  ROdata%bangle_method     = "UNKNOWN" ! Bending angle processing method
  ROdata%refrac_method     = "UNKNOWN" ! Refractivity processing method
  ROdata%meteo_method      = "UNKNOWN" ! Meteorological processing method
  ROdata%thin_method       = "UNKNOWN" ! Profile thinning method
  ROdata%software_version  = "UNKNOWN" ! Software version ID

  ROdata%DTpro%Year   = 9999
  ROdata%DTpro%Month  =   99
  ROdata%DTpro%Day    =   99
  ROdata%DTpro%Hour   =   99
  ROdata%DTpro%Minute =   99
  ROdata%DTpro%Second =   99
  ROdata%DTpro%Msec   =  999

  ROdata%DTocc%Year   = 9999
  ROdata%DTocc%Month  =   99
  ROdata%DTocc%Day    =   99
  ROdata%DTocc%Hour   =   99
  ROdata%DTocc%Minute =   99
  ROdata%DTocc%Second =   99
  ROdata%DTocc%Msec   =  999

  ROdata%PCD          = 0
  ROdata%PCD          = IBSET ( ROdata%PCD, PCD_missing )

  ROdata%overall_qual = ropp_MDFV         ! Overall quality value

  ROdata%georef%time_offset =  ropp_MDFV  ! Time since start (s)
  ROdata%georef%lat         =  ropp_MDFV  ! Latitude  (deg)
  ROdata%georef%lon         =  ropp_MDFV  ! Longitude (deg)
  ROdata%georef%roc         =  ropp_MDFV  ! RoC value (m)
  ROdata%georef%r_coc       =  ropp_MDFV  ! RoC offset X,Y,Z vector (m)
  ROdata%georef%azimuth     =  ropp_MDFV  ! GNSS->LEO line of sight angle (degT)
  ROdata%georef%undulation  =  ropp_MDFV  ! Geoid undulation (EGM96-WGS84) (m)

  ROdata%bg%source   = 'NONE' ! Source of b/g profile
  ROdata%bg%year     = 9999   ! VT year   of b/g
  ROdata%bg%month    =   99   ! VT month  of b/g
  ROdata%bg%day      =   99   ! VT day    of b/g
  ROdata%bg%hour     =   99   ! VT hour   of b/g
  ROdata%bg%minute   =   99   ! VT minute of b/g
  ROdata%bg%fcperiod =  999   ! F/c period (hrs)

! 8.3 Level 1a profile
! --------------------

  CALL ropp_io_free(ROdata%Lev1a)

! 8.4 Level 1b profile
! --------------------

  CALL ropp_io_free(ROdata%Lev1b)

! 8.5 Level 2a profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2a)

! 8.6 Level 2b profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2b)

! 8.7 Level 2c profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2c)

! 8.8 Level 2d profile
! --------------------

  CALL ropp_io_free(ROdata%Lev2d)

! 8.9 Additional variables
! ------------------------

  CALL ropp_io_free(ROdata%vlist)

END SUBROUTINE ropp_io_free_rotype_2d

!-------------------------------------------------------------------------------
! 9. Vlist (D0d)
!-------------------------------------------------------------------------------

RECURSIVE SUBROUTINE ropp_io_free_vlistD0d(vlist)

! 9.1 Declarations
! ----------------

! USE ropp_io,       not_this => ropp_io_free_vlistD0d
  USE ropp_io_types, ONLY: VlisttypeD0d

  IMPLICIT NONE

  TYPE(VlisttypeD0d), POINTER :: vlist

! 9.2 Free next member of the list
! --------------------------------

  IF (ASSOCIATED(vlist)) THEN

     IF (ASSOCIATED(vlist%next)) THEN
        CALL ropp_io_free_vlistD0d(vlist%next)
     ENDIF

! 9.3 Free this member of the list
! --------------------------------

     DEALLOCATE(vlist)

  ENDIF

END SUBROUTINE ropp_io_free_vlistD0d


!-------------------------------------------------------------------------------
! 10. Vlist (D1d)
!-------------------------------------------------------------------------------

RECURSIVE SUBROUTINE ropp_io_free_vlistD1d(vlist)

! 10.1 Declarations
! -----------------

! USE ropp_io,       not_this => ropp_io_free_vlistD1d
  USE ropp_io_types, ONLY: VlisttypeD1d

  IMPLICIT NONE

  TYPE(VlisttypeD1d), POINTER :: vlist

! 10.2 Free next member of the list
! ---------------------------------

  IF (ASSOCIATED(vlist)) THEN

     IF (ASSOCIATED(vlist%next)) THEN
        CALL ropp_io_free_vlistD1d(vlist%next)
     ENDIF

! 10.3 Free this member of the list
! ---------------------------------

     IF (ASSOCIATED(vlist%data)) THEN
        DEALLOCATE(vlist%data)
     ENDIF

     DEALLOCATE(vlist)

  ENDIF

END SUBROUTINE ropp_io_free_vlistD1d


!-------------------------------------------------------------------------------
! 11. Vlist (D2d)
!-------------------------------------------------------------------------------

RECURSIVE SUBROUTINE ropp_io_free_vlistD2d(vlist)

! 11.1 Declarations
! -----------------

! USE ropp_io,       not_this => ropp_io_free_vlistD2d
  USE ropp_io_types, ONLY: VlisttypeD2d

  IMPLICIT NONE

  TYPE(VlisttypeD2d), POINTER :: vlist

! 11.2 Free next member of the list
! ---------------------------------

  IF (ASSOCIATED(vlist)) THEN

     IF (ASSOCIATED(vlist%next)) THEN
        CALL ropp_io_free_vlistD2d(vlist%next)
     ENDIF

! 11.3 Free this member of the list
! ---------------------------------

     IF (ASSOCIATED(vlist%data)) THEN
        DEALLOCATE(vlist%data)
     ENDIF

     DEALLOCATE(vlist)

  ENDIF

END SUBROUTINE ropp_io_free_vlistD2d


!-------------------------------------------------------------------------------
! 12. Vlist (generic)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_vlist(vlist)

! 12.1 Declarations
! -----------------

! USE ropp_io,       not_this => ropp_io_free_vlist
  USE ropp_io_types, ONLY: Vlisttype

  IMPLICIT NONE

  TYPE(Vlisttype), INTENT(inout) :: vlist

! 12.2 Free members of the list
! -----------------------------

  IF (ASSOCIATED(vlist%VlistD0d)) THEN
     CALL ropp_io_free_vlistD0d(vlist%VlistD0d)
  ENDIF
  IF (ASSOCIATED(vlist%VlistD1d)) THEN
     CALL ropp_io_free_vlistD1d(vlist%VlistD1d)
  ENDIF
  IF (ASSOCIATED(vlist%VlistD2d)) THEN
     CALL ropp_io_free_vlistD2d(vlist%VlistD2d)
  ENDIF

END SUBROUTINE ropp_io_free_vlist


!-------------------------------------------------------------------------------
! 13. Error correlation / covariance structure (scalar)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_corcov_sca(corcov)

! 13.1 Declarations
! -----------------

  USE typesizes,     ONLY: wp => EightByteReal
! USE ropp_io,       not_this => ropp_io_free_corcov_sca
  USE ropp_io_types, ONLY: ROcorcov

  IMPLICIT NONE

  TYPE(ROcorcov), INTENT(inout) :: corcov

! 13.2 Reinitialise latitude boundaries
! -------------------------------------

  corcov%lat_min = -90.0_wp
  corcov%lat_max =  90.0_wp

! 13.3 Free members of the structure
! ----------------------------------

  IF (ASSOCIATED(corcov%corr)) THEN
     DEALLOCATE(corcov%corr)
  ENDIF

  IF (ASSOCIATED(corcov%sigma)) THEN
     DEALLOCATE(corcov%sigma)
  ENDIF

END SUBROUTINE ropp_io_free_corcov_sca


!-------------------------------------------------------------------------------
! 14. Error correlation / covariance structure (array)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_free_corcov_arr(corcov)

! 14.1 Declarations
! -----------------

! USE ropp_io,       not_this => ropp_io_free_corcov_arr
  USE ropp_io_types, ONLY: ROcorcov

  IMPLICIT NONE

  TYPE(ROcorcov), DIMENSION(:), POINTER :: corcov

  INTEGER                               :: i

! 14.2 Free memeory in all elements
! ---------------------------------

  IF (ASSOCIATED(corcov)) THEN
     IF (SIZE(corcov) > 0) THEN
        DO i = 1, SIZE(corcov)
           CALL ropp_io_free(corcov(i))
        ENDDO
        DEALLOCATE(corcov)
     ENDIF
  ENDIF

END SUBROUTINE ropp_io_free_corcov_arr
