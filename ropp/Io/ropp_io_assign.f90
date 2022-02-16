! $Id: ropp_io_assign.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Initialisation/ropp_io_assign *
!
! NAME
!    ropp_io_assign - Initialise an RO derived type with another RO derived type
!
! SYNOPSIS
!    use ropp_io
!    type(ROprof) :: ro_data
!      ...
!    call ropp_io_assign(ro_data_old, ro_data_new) 
!    call ropp_io_assign(ro_data%Lev1a_old, ro_data%Lev1a_old)
!    call ropp_io_assign(ro_data%Lev1b_old, ro_data%Lev1b_old)
!    call ropp_io_assign(ro_data%Lev2a_old, ro_data%Lev2a_old)
!    call ropp_io_assign(ro_data%Lev2b_old, ro_data%Lev2b_old)
!    call ropp_io_assign(ro_data%Lev2c_old, ro_data%Lev2c_old)
!    call ropp_io_assign(ro_data%Lev2d_old, ro_data%Lev2d_old)
!
! DESCRIPTION
!   This subroutine initialises a RO data structure or parts thereof with the 
!   elements of an exisitng RO data structure.
!   These subroutines can be used to overloading the assign operator (=)
!   for RO data structures
!
! INPUT
!    ro_data_old   dtyp  RO data (derived type)
!
! OUTPUT
!    ro_data_new   dtyp  RO data (derived type)
!
! SEE ALSO
!    ropp_io_types
!    ropp_io_init
!
! NOTES
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

SUBROUTINE ropp_io_l1atype_l1atype(from_var, to_var)

! 1.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_l1atype_l1atype
  USE ropp_io_types, ONLY: L1atype

  IMPLICIT NONE

  TYPE(L1atype), INTENT(in)    :: from_var
  TYPE(L1atype), INTENT(inout) :: to_var

! 1.2 Copy contents
! -----------------  

  to_var % Npoints = from_var % Npoints

  to_var % Missing = from_var % Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var % dtime(1:from_var%Npoints), to_var % dtime )
  CALL copy_alloc( from_var % snr_L1ca(1:from_var % Npoints)   , to_var % snr_L1ca   )
  CALL copy_alloc( from_var % snr_L1p(1:from_var % Npoints)    , to_var % snr_L1p    )
  CALL copy_alloc( from_var % snr_L2p(1:from_var % Npoints)   , to_var % snr_L2p    )
  CALL copy_alloc( from_var % phase_L1(1:from_var % Npoints)   , to_var % phase_L1   )
  CALL copy_alloc( from_var % phase_L2(1:from_var % Npoints)   , to_var % phase_L2   )
  CALL copy_alloc( from_var % phase_qual(1:from_var % Npoints) , to_var % phase_qual )
  CALL copy_alloc( from_var % r_gns(1:from_var % Npoints, :)      , to_var % r_gns      )
  CALL copy_alloc( from_var % v_gns(1:from_var % Npoints, :)      , to_var % v_gns      )
  CALL copy_alloc( from_var % r_leo(1:from_var % Npoints, :)      , to_var % r_leo      )
  CALL copy_alloc( from_var % v_leo(1:from_var % Npoints, :)      , to_var % v_leo      )

  ENDIF

END SUBROUTINE ropp_io_l1atype_l1atype

!-------------------------------------------------------------------------------
! 2. L1atype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l1atype(to_var, from_var)

! 2.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_assign_l1atype
  USE ropp_io_types, ONLY: L1atype

  IMPLICIT NONE

  TYPE(L1atype), INTENT(in)    :: from_var
  TYPE(L1atype), INTENT(inout) :: to_var

! 2.2 Copy contents
! -----------------  

  to_var % Npoints = from_var % Npoints

  to_var % Missing = from_var % Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var % dtime(1:from_var % Npoints)      , to_var % dtime      )
  CALL copy_alloc( from_var % snr_L1ca(1:from_var % Npoints)   , to_var % snr_L1ca   )
  CALL copy_alloc( from_var % snr_L1p(1:from_var % Npoints)    , to_var % snr_L1p    )
  CALL copy_alloc( from_var % snr_L2p(1:from_var % Npoints)    , to_var % snr_L2p    )
  CALL copy_alloc( from_var % phase_L1(1:from_var % Npoints)   , to_var % phase_L1   )
  CALL copy_alloc( from_var % phase_L2(1:from_var % Npoints)   , to_var % phase_L2   )
  CALL copy_alloc( from_var % phase_qual(1:from_var % Npoints) , to_var % phase_qual )
  CALL copy_alloc( from_var % r_gns(1:from_var % Npoints, :)      , to_var % r_gns   )   
  CALL copy_alloc( from_var % v_gns(1:from_var % Npoints, :)      , to_var % v_gns      )
  CALL copy_alloc( from_var % r_leo(1:from_var % Npoints, :)      , to_var % r_leo      )
  CALL copy_alloc( from_var % v_leo(1:from_var % Npoints, :)      , to_var % v_leo      )

  ENDIF

END SUBROUTINE ropp_io_assign_l1atype


!-------------------------------------------------------------------------------
! 3. L1btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_l1btype_l1btype(from_var, to_var)

! 3.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_l1btype_l1btype
  USE ropp_io_types, ONLY: L1btype

  IMPLICIT NONE

  TYPE(L1btype), INTENT(in)    :: from_var
  TYPE(L1btype), INTENT(inout) :: to_var

! 3.2 Copy contents
! -----------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%lat_tp(1:from_var % Npoints)           , to_var%lat_tp           )   
  CALL copy_alloc( from_var%lon_tp(1:from_var % Npoints)           , to_var%lon_tp           ) 
  CALL copy_alloc( from_var%azimuth_tp(1:from_var % Npoints)       , to_var%azimuth_tp       )  
  CALL copy_alloc( from_var%impact_L1(1:from_var % Npoints)        , to_var%impact_L1        )     
  CALL copy_alloc( from_var%impact_L2 (1:from_var % Npoints)       , to_var%impact_L2        )     
  CALL copy_alloc( from_var%impact(1:from_var % Npoints)           , to_var%impact           )     
  CALL copy_alloc( from_var%impact_Opt(1:from_var % Npoints)       , to_var%impact_Opt       )     
  CALL copy_alloc( from_var%bangle_L1(1:from_var % Npoints)        , to_var%bangle_L1        )  
  CALL copy_alloc( from_var%bangle_L2(1:from_var % Npoints)        , to_var%bangle_L2        )    
  CALL copy_alloc( from_var%bangle(1:from_var % Npoints)           , to_var%bangle           )    
  CALL copy_alloc( from_var%bangle_Opt(1:from_var % Npoints)       , to_var%bangle_Opt       )
  CALL copy_alloc( from_var%bangle_L1_sigma(1:from_var % Npoints)  , to_var%bangle_L1_sigma  ) 
  CALL copy_alloc( from_var%bangle_L2_sigma(1:from_var % Npoints)  , to_var%bangle_L2_sigma  )  
  CALL copy_alloc( from_var%bangle_sigma(1:from_var % Npoints)     , to_var%bangle_sigma     )
  CALL copy_alloc( from_var%bangle_Opt_sigma(1:from_var % Npoints) , to_var%bangle_Opt_sigma )
  CALL copy_alloc( from_var%bangle_L1_qual(1:from_var % Npoints)   , to_var%bangle_L1_qual   )
  CALL copy_alloc( from_var%bangle_L2_qual(1:from_var % Npoints)   , to_var%bangle_L2_qual   ) 
  CALL copy_alloc( from_var%bangle_qual(1:from_var % Npoints)      , to_var%bangle_qual      )  
  CALL copy_alloc( from_var%bangle_Opt_qual(1:from_var % Npoints)  , to_var%bangle_Opt_qual  )  
  
  ENDIF

END SUBROUTINE ropp_io_l1btype_l1btype

!-------------------------------------------------------------------------------
! 4. L1btype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l1btype(to_var, from_var)

! 4.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_assign_l1btype
  USE ropp_io_types, ONLY: L1btype

  IMPLICIT NONE

  TYPE(L1btype), INTENT(in)    :: from_var
  TYPE(L1btype), INTENT(inout) :: to_var

! 4.2 Copy contents
! -----------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%lat_tp(1:from_var % Npoints)           , to_var%lat_tp           )   
  CALL copy_alloc( from_var%lon_tp(1:from_var % Npoints)           , to_var%lon_tp           ) 
  CALL copy_alloc( from_var%azimuth_tp(1:from_var % Npoints)       , to_var%azimuth_tp       )  
  CALL copy_alloc( from_var%impact_L1(1:from_var % Npoints)       , to_var%impact_L1        )     
  CALL copy_alloc( from_var%impact_L2(1:from_var % Npoints)        , to_var%impact_L2        )     
  CALL copy_alloc( from_var%impact(1:from_var % Npoints)           , to_var%impact           )     
  CALL copy_alloc( from_var%impact_Opt(1:from_var % Npoints)       , to_var%impact_Opt       )     
  CALL copy_alloc( from_var%bangle_L1(1:from_var % Npoints)        , to_var%bangle_L1        )  
  CALL copy_alloc( from_var%bangle_L2(1:from_var % Npoints)       , to_var%bangle_L2        )    
  CALL copy_alloc( from_var%bangle(1:from_var % Npoints)           , to_var%bangle           )    
  CALL copy_alloc( from_var%bangle_Opt(1:from_var % Npoints)       , to_var%bangle_Opt       )
  CALL copy_alloc( from_var%bangle_L1_sigma(1:from_var % Npoints)  , to_var%bangle_L1_sigma  ) 
  CALL copy_alloc( from_var%bangle_L2_sigma(1:from_var % Npoints)  , to_var%bangle_L2_sigma  )  
  CALL copy_alloc( from_var%bangle_sigma(1:from_var % Npoints)    , to_var%bangle_sigma     )
  CALL copy_alloc( from_var%bangle_Opt_sigma(1:from_var % Npoints) , to_var%bangle_Opt_sigma )
  CALL copy_alloc( from_var%bangle_L1_qual(1:from_var % Npoints)   , to_var%bangle_L1_qual   )
  CALL copy_alloc( from_var%bangle_L2_qual(1:from_var % Npoints)   , to_var%bangle_L2_qual   ) 
  CALL copy_alloc( from_var%bangle_qual(1:from_var % Npoints)      , to_var%bangle_qual      )  
  CALL copy_alloc( from_var%bangle_Opt_qual(1:from_var % Npoints)  , to_var%bangle_Opt_qual  )  
  
  ENDIF

END SUBROUTINE ropp_io_assign_l1btype

!-------------------------------------------------------------------------------
! 5. L2atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_l2atype_l2atype(from_var, to_var)

! 5.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_l2atype_l2atype
  USE ropp_io_types, ONLY: L2atype

  IMPLICIT NONE

  TYPE(L2atype), INTENT(in)    :: from_var
  TYPE(L2atype), INTENT(inout) :: to_var

! 5.2 Allocate memory for all structure elements
! ----------------------------------------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%alt_refrac(1:from_var % Npoints)     , to_var%alt_refrac     )   
  CALL copy_alloc( from_var%geop_refrac(1:from_var % Npoints)    , to_var%geop_refrac    ) 
  CALL copy_alloc( from_var%refrac(1:from_var % Npoints)         , to_var%refrac         ) 
  CALL copy_alloc( from_var%refrac_sigma(1:from_var % Npoints)   , to_var%refrac_sigma   ) 
  CALL copy_alloc( from_var%refrac_qual(1:from_var % Npoints)    , to_var%refrac_qual    ) 
  CALL copy_alloc( from_var%dry_temp(1:from_var % Npoints)       , to_var%dry_temp       ) 
  CALL copy_alloc( from_var%dry_temp_sigma(1:from_var % Npoints) , to_var%dry_temp_sigma ) 
  CALL copy_alloc( from_var%dry_temp_qual(1:from_var % Npoints)  , to_var%dry_temp_qual  ) 

  ENDIF
  
END SUBROUTINE ropp_io_l2atype_l2atype

!-------------------------------------------------------------------------------
! 6. L2atype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l2atype(to_var, from_var)

! 6.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_assign_l2atype
  USE ropp_io_types, ONLY: L2atype

  IMPLICIT NONE

  TYPE(L2atype), INTENT(in)    :: from_var
  TYPE(L2atype), INTENT(inout) :: to_var

! 6.2 Allocate memory for all structure elements
! ----------------------------------------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%alt_refrac(1:from_var % Npoints)     , to_var%alt_refrac     )   
  CALL copy_alloc( from_var%geop_refrac(1:from_var % Npoints)    , to_var%geop_refrac    ) 
  CALL copy_alloc( from_var%refrac(1:from_var % Npoints)         , to_var%refrac         )       
  CALL copy_alloc( from_var%refrac_sigma(1:from_var % Npoints)   , to_var%refrac_sigma   ) 
  CALL copy_alloc( from_var%refrac_qual(1:from_var % Npoints)    , to_var%refrac_qual    )  
  CALL copy_alloc( from_var%dry_temp(1:from_var % Npoints)       , to_var%dry_temp       )       
  CALL copy_alloc( from_var%dry_temp_sigma(1:from_var % Npoints) , to_var%dry_temp_sigma ) 
  CALL copy_alloc( from_var%dry_temp_qual(1:from_var % Npoints)  , to_var%dry_temp_qual  )  

  ENDIF

END SUBROUTINE ropp_io_assign_l2atype

!-------------------------------------------------------------------------------
! 7. L2btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_l2btype_l2btype(from_var, to_var)

! 7.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_l2btype_l2btype
  USE ropp_io_types, ONLY: L2btype

  IMPLICIT NONE

  TYPE(L2btype), INTENT(in)    :: from_var
  TYPE(L2btype), INTENT(inout) :: to_var

! 7.2 Copy contents
! -----------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%geop(1:from_var % Npoints)        , to_var%geop        )
  CALL copy_alloc( from_var%geop_sigma(1:from_var % Npoints)  , to_var%geop_sigma  )
  CALL copy_alloc( from_var%press(1:from_var % Npoints)       , to_var%press       )
  CALL copy_alloc( from_var%press_sigma(1:from_var % Npoints) , to_var%press_sigma )
  CALL copy_alloc( from_var%temp(1:from_var % Npoints)        , to_var%temp        )
  CALL copy_alloc( from_var%temp_sigma(1:from_var % Npoints)  , to_var%temp_sigma  )
  CALL copy_alloc( from_var%shum(1:from_var % Npoints)        , to_var%shum        )
  CALL copy_alloc( from_var%shum_sigma(1:from_var % Npoints)  , to_var%shum_sigma  ) 
  CALL copy_alloc( from_var%meteo_qual(1:from_var % Npoints)  , to_var%meteo_qual  )  

  ENDIF

END SUBROUTINE ropp_io_l2btype_l2btype

!-------------------------------------------------------------------------------
! 8. L2btype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l2btype(to_var, from_var)

! 8.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_assign_l2btype
  USE ropp_io_types, ONLY: L2btype

  IMPLICIT NONE

  TYPE(L2btype), INTENT(in)    :: from_var
  TYPE(L2btype), INTENT(inout) :: to_var

! 8.2 Copy contents
! -----------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%geop(1:from_var % Npoints)        , to_var%geop        )
  CALL copy_alloc( from_var%geop_sigma(1:from_var % Npoints)  , to_var%geop_sigma  )
  CALL copy_alloc( from_var%press(1:from_var % Npoints)       , to_var%press       )
  CALL copy_alloc( from_var%press_sigma(1:from_var % Npoints) , to_var%press_sigma )
  CALL copy_alloc( from_var%temp(1:from_var % Npoints)        , to_var%temp        )
  CALL copy_alloc( from_var%temp_sigma(1:from_var % Npoints)  , to_var%temp_sigma  )
  CALL copy_alloc( from_var%shum(1:from_var % Npoints)        , to_var%shum        )
  CALL copy_alloc( from_var%shum_sigma(1:from_var % Npoints)  , to_var%shum_sigma  ) 
  CALL copy_alloc( from_var%meteo_qual(1:from_var % Npoints)  , to_var%meteo_qual  )  

  ENDIF

END SUBROUTINE ropp_io_assign_l2btype


!-------------------------------------------------------------------------------
! 9. L2ctype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_l2ctype_l2ctype(from_var, to_var)

! 9.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io, not_this => ropp_io_l2ctype_l2ctype
  USE ropp_io_types, ONLY: L2ctype

  IMPLICIT NONE

  TYPE(L2ctype), INTENT(in)    :: from_var
  TYPE(L2ctype), INTENT(inout) :: to_var

! 9.2 Copy contents
! -----------------

  to_var%Npoints           = from_var%Npoints
  to_var%Missing           = from_var%Missing
  to_var%geop_sfc          = from_var%geop_sfc
  to_var%press_sfc         = from_var%press_sfc
  to_var%press_sfc_sigma   = from_var%press_sfc_sigma
  to_var%press_sfc_qual    = from_var%press_sfc_qual

  to_var%Ne_max            = from_var%Ne_max
  to_var%Ne_max_sigma      = from_var%Ne_max_sigma

  to_var%H_peak            = from_var%H_peak
  to_var%H_peak_sigma      = from_var%H_peak_sigma

  to_var%H_width           = from_var%H_width
  to_var%H_width_sigma     = from_var%H_width_sigma

  to_var%tph_bangle        = from_var%tph_bangle
  to_var%tpa_bangle        = from_var%tpa_bangle
  to_var%tph_bangle_flag   = from_var%tph_bangle_flag

  to_var%tph_refrac        = from_var%tph_refrac
  to_var%tpn_refrac        = from_var%tpn_refrac
  to_var%tph_refrac_flag   = from_var%tph_refrac_flag

  to_var%tph_tdry_lrt      = from_var%tph_tdry_lrt
  to_var%tpt_tdry_lrt      = from_var%tpt_tdry_lrt
  to_var%tph_tdry_lrt_flag = from_var%tph_tdry_lrt_flag

  to_var%tph_tdry_cpt      = from_var%tph_tdry_cpt
  to_var%tpt_tdry_cpt      = from_var%tpt_tdry_cpt
  to_var%tph_tdry_cpt_flag = from_var%tph_tdry_cpt_flag

  to_var%prh_tdry_cpt      = from_var%prh_tdry_cpt
  to_var%prt_tdry_cpt      = from_var%prt_tdry_cpt
  to_var%prh_tdry_cpt_flag = from_var%prh_tdry_cpt_flag

  to_var%tph_temp_lrt      = from_var%tph_temp_lrt
  to_var%tpt_temp_lrt      = from_var%tpt_temp_lrt
  to_var%tph_temp_lrt_flag = from_var%tph_temp_lrt_flag

  to_var%tph_temp_cpt      = from_var%tph_temp_cpt
  to_var%tpt_temp_cpt      = from_var%tpt_temp_cpt
  to_var%tph_temp_cpt_flag = from_var%tph_temp_cpt_flag

  to_var%prh_temp_cpt      = from_var%prh_temp_cpt
  to_var%prt_temp_cpt      = from_var%prt_temp_cpt
  to_var%prh_temp_cpt_flag = from_var%prh_temp_cpt_flag

END SUBROUTINE ropp_io_l2ctype_l2ctype

!-------------------------------------------------------------------------------
! 10. L2ctype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l2ctype(to_var, from_var)

! 10.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io, not_this => ropp_io_assign_l2ctype
  USE ropp_io_types, ONLY: L2ctype

  IMPLICIT NONE

  TYPE(L2ctype), INTENT(in)    :: from_var
  TYPE(L2ctype), INTENT(inout) :: to_var

! 10.2 Copy contents
! -----------------

  to_var%Npoints           = from_var%Npoints
  to_var%Missing           = from_var%Missing
  to_var%geop_sfc          = from_var%geop_sfc
  to_var%press_sfc         = from_var%press_sfc
  to_var%press_sfc_sigma   = from_var%press_sfc_sigma
  to_var%press_sfc_qual    = from_var%press_sfc_qual

  to_var%Ne_max            = from_var%Ne_max
  to_var%Ne_max_sigma      = from_var%Ne_max_sigma

  to_var%H_peak            = from_var%H_peak
  to_var%H_peak_sigma      = from_var%H_peak_sigma

  to_var%H_width           = from_var%H_width
  to_var%H_width_sigma     = from_var%H_width_sigma

  to_var%tph_bangle        = from_var%tph_bangle
  to_var%tpa_bangle        = from_var%tpa_bangle
  to_var%tph_bangle_flag   = from_var%tph_bangle_flag

  to_var%tph_refrac        = from_var%tph_refrac
  to_var%tpn_refrac        = from_var%tpn_refrac
  to_var%tph_refrac_flag   = from_var%tph_refrac_flag

  to_var%tph_tdry_lrt      = from_var%tph_tdry_lrt
  to_var%tpt_tdry_lrt      = from_var%tpt_tdry_lrt
  to_var%tph_tdry_lrt_flag = from_var%tph_tdry_lrt_flag

  to_var%tph_tdry_cpt      = from_var%tph_tdry_cpt
  to_var%tpt_tdry_cpt      = from_var%tpt_tdry_cpt
  to_var%tph_tdry_cpt_flag = from_var%tph_tdry_cpt_flag

  to_var%prh_tdry_cpt      = from_var%prh_tdry_cpt
  to_var%prt_tdry_cpt      = from_var%prt_tdry_cpt
  to_var%prh_tdry_cpt_flag = from_var%prh_tdry_cpt_flag

  to_var%tph_temp_lrt      = from_var%tph_temp_lrt
  to_var%tpt_temp_lrt      = from_var%tpt_temp_lrt
  to_var%tph_temp_lrt_flag = from_var%tph_temp_lrt_flag

  to_var%tph_temp_cpt      = from_var%tph_temp_cpt
  to_var%tpt_temp_cpt      = from_var%tpt_temp_cpt
  to_var%tph_temp_cpt_flag = from_var%tph_temp_cpt_flag

  to_var%prh_temp_cpt      = from_var%prh_temp_cpt
  to_var%prt_temp_cpt      = from_var%prt_temp_cpt
  to_var%prh_temp_cpt_flag = from_var%prh_temp_cpt_flag

END SUBROUTINE ropp_io_assign_l2ctype


!-------------------------------------------------------------------------------
! 11. L2dtype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_l2dtype_l2dtype(from_var, to_var)

! 11.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_l2dtype_l2dtype
  USE ropp_io_types, ONLY: L2dtype

  IMPLICIT NONE

  TYPE(L2dtype), INTENT(in)    :: from_var
  TYPE(L2dtype), INTENT(inout) :: to_var

! 11.2 Copy contents
! ------------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  to_var%level_type = from_var%level_type

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%level_coeff_a(1:from_var % Npoints) , to_var%level_coeff_a )
  CALL copy_alloc( from_var%level_coeff_b(1:from_var % Npoints) , to_var%level_coeff_b )

  ENDIF

END SUBROUTINE ropp_io_l2dtype_l2dtype

!-------------------------------------------------------------------------------
! 12. L2dtype (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_l2dtype(to_var, from_var)

! 12.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_assign_l2dtype
  USE ropp_io_types, ONLY: L2dtype

  IMPLICIT NONE

  TYPE(L2dtype), INTENT(in)    :: from_var
  TYPE(L2dtype), INTENT(inout) :: to_var

! 12.2 Copy contents
! ------------------

  to_var%Npoints = from_var%Npoints

  to_var%Missing = from_var%Missing

  to_var%level_type = from_var%level_type

  IF ( to_var%Npoints > 0 ) THEN

  CALL copy_alloc( from_var%level_coeff_a(1:from_var % Npoints) , to_var%level_coeff_a )
  CALL copy_alloc( from_var%level_coeff_b(1:from_var % Npoints) , to_var%level_coeff_b )

  ENDIF

END SUBROUTINE ropp_io_assign_l2dtype

!-------------------------------------------------------------------------------
! 13. Joint RO data type
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_rotype_rotype(from_ROdata, to_ROdata)

! 13.1 Declarations
! ----------------

! USE ropp_io,        not_this => ropp_io_rotype_rotype
  USE ropp_io_types,  ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof), INTENT(in)    :: from_ROdata
  TYPE(ROprof), INTENT(inout) :: to_ROdata

! 13.2 Copy attributes
! --------------------
  
  to_ROdata%FmtVersion = from_ROdata%FmtVersion
  to_ROdata%occ_id = from_ROdata%occ_id
  to_ROdata%leo_id = from_ROdata%leo_id
  to_ROdata%gns_id = from_ROdata%gns_id
  to_ROdata%stn_id = from_ROdata%stn_id
  to_ROdata%processing_centre = from_ROdata%processing_centre
  to_ROdata%processing_software = from_ROdata%processing_software
  to_ROdata%pod_method = from_ROdata%pod_method
  to_ROdata%phase_method = from_ROdata%phase_method
  to_ROdata%bangle_method = from_ROdata%bangle_method
  to_ROdata%refrac_method = from_ROdata%refrac_method
  to_ROdata%meteo_method = from_ROdata%meteo_method
  to_ROdata%thin_method = from_ROdata%thin_method
  to_ROdata%software_version = from_ROdata%software_version

! 13.3 Copy time contents
! -----------------------

  to_ROdata%DTocc%Year   = from_ROdata%DTocc%Year
  to_ROdata%DTocc%Month  = from_ROdata%DTocc%Month
  to_ROdata%DTocc%Day    = from_ROdata%DTocc%Day
  to_ROdata%DTocc%Hour   = from_ROdata%DTocc%Hour 
  to_ROdata%DTocc%Minute = from_ROdata%DTocc%Minute
  to_ROdata%DTocc%Second = from_ROdata%DTocc%Second
  to_ROdata%DTocc%Msec   = from_ROdata%DTocc%Msec

  to_ROdata%DTpro%Year   = from_ROdata%DTpro%Year
  to_ROdata%DTpro%Month  = from_ROdata%DTpro%Month
  to_ROdata%DTpro%Day    = from_ROdata%DTpro%Day
  to_ROdata%DTpro%Hour   = from_ROdata%DTpro%Hour 
  to_ROdata%DTpro%Minute = from_ROdata%DTpro%Minute
  to_ROdata%DTpro%Second = from_ROdata%DTpro%Second
  to_ROdata%DTpro%Msec   = from_ROdata%DTpro%Msec

! 13.4 Copy PCD and quality 
! -------------------------

  to_ROdata%PCD = from_ROdata%PCD 
  to_ROdata%overall_qual = from_ROdata%overall_qual 

! 13.5 Copy georeferencing
! ------------------------

  to_ROdata%GEOref%time_offset = from_ROdata%GEOref%time_offset 
  to_ROdata%GEOref%lat = from_ROdata%GEOref%lat
  to_ROdata%GEOref%lon = from_ROdata%GEOref%lon
  to_ROdata%GEOref%roc = from_ROdata%GEOref%roc 
  to_ROdata%GEOref%r_coc(1:3) = from_ROdata%GEOref%r_coc(1:3) 
  to_ROdata%GEOref%azimuth = from_ROdata%GEOref%azimuth
  to_ROdata%GEOref%undulation = from_ROdata%GEOref%undulation

! 13.6 Level 1a profile
! --------------------

  to_ROdata%Lev1a = from_ROdata%Lev1a

! 13.7 Level 1b profile
! --------------------

  to_ROdata%Lev1b = from_ROdata%Lev1b

! 13.8 Level 2a profile
! --------------------

  to_ROdata%Lev2a = from_ROdata%Lev2a

! 13.9 Level 2b profile
! --------------------

  to_ROdata%Lev2b = from_ROdata%Lev2b

! 13.10 Level 2c profile
! --------------------

  to_ROdata%Lev2c = from_ROdata%Lev2c

! 13.11 Level 2d profile
! --------------------

  to_ROdata%Lev2d = from_ROdata%Lev2d

END SUBROUTINE ropp_io_rotype_rotype


!-------------------------------------------------------------------------------
! 14. Joint RO data type (assignment)
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_assign_rotype(to_ROdata, from_ROdata)

! 14.1 Declarations
! ----------------

! USE ropp_io,        not_this => ropp_io_assign_rotype
  USE ropp_io_types,  ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof), INTENT(in)    :: from_ROdata
  TYPE(ROprof), INTENT(inout) :: to_ROdata

! 14.2 Copy attributes
! --------------------
  
  to_ROdata%FmtVersion = from_ROdata%FmtVersion
  to_ROdata%occ_id = from_ROdata%occ_id
  to_ROdata%leo_id = from_ROdata%leo_id
  to_ROdata%gns_id = from_ROdata%gns_id
  to_ROdata%stn_id = from_ROdata%stn_id
  to_ROdata%processing_centre = from_ROdata%processing_centre
  to_ROdata%processing_software = from_ROdata%processing_software
  to_ROdata%pod_method = from_ROdata%pod_method
  to_ROdata%phase_method = from_ROdata%phase_method
  to_ROdata%bangle_method = from_ROdata%bangle_method
  to_ROdata%refrac_method = from_ROdata%refrac_method
  to_ROdata%meteo_method = from_ROdata%meteo_method
  to_ROdata%thin_method = from_ROdata%thin_method
  to_ROdata%software_version = from_ROdata%software_version

! 14.3 Copy time contents
! -----------------------

  to_ROdata%DTocc%Year   = from_ROdata%DTocc%Year
  to_ROdata%DTocc%Month  = from_ROdata%DTocc%Month
  to_ROdata%DTocc%Day    = from_ROdata%DTocc%Day
  to_ROdata%DTocc%Hour   = from_ROdata%DTocc%Hour 
  to_ROdata%DTocc%Minute = from_ROdata%DTocc%Minute
  to_ROdata%DTocc%Second = from_ROdata%DTocc%Second
  to_ROdata%DTocc%Msec   = from_ROdata%DTocc%Msec

  to_ROdata%DTpro%Year   = from_ROdata%DTpro%Year
  to_ROdata%DTpro%Month  = from_ROdata%DTpro%Month
  to_ROdata%DTpro%Day    = from_ROdata%DTpro%Day
  to_ROdata%DTpro%Hour   = from_ROdata%DTpro%Hour 
  to_ROdata%DTpro%Minute = from_ROdata%DTpro%Minute
  to_ROdata%DTpro%Second = from_ROdata%DTpro%Second
  to_ROdata%DTpro%Msec   = from_ROdata%DTpro%Msec

! 14.4 Copy PCD and quality 
! -------------------------

  to_ROdata%PCD = from_ROdata%PCD 
  to_ROdata%overall_qual = from_ROdata%overall_qual 

! 14.5 Copy georeferencing
! ------------------------

  to_ROdata%GEOref%time_offset = from_ROdata%GEOref%time_offset 
  to_ROdata%GEOref%lat = from_ROdata%GEOref%lat
  to_ROdata%GEOref%lon = from_ROdata%GEOref%lon
  to_ROdata%GEOref%roc = from_ROdata%GEOref%roc 
  to_ROdata%GEOref%r_coc(1:3) = from_ROdata%GEOref%r_coc(1:3) 
  to_ROdata%GEOref%azimuth = from_ROdata%GEOref%azimuth
  to_ROdata%GEOref%undulation = from_ROdata%GEOref%undulation

! 14.6 Level 1a profile
! --------------------

  to_ROdata%Lev1a = from_ROdata%Lev1a

! 14.7 Level 1b profile
! --------------------

  to_ROdata%Lev1b = from_ROdata%Lev1b

! 14.8 Level 2a profile
! --------------------

  to_ROdata%Lev2a = from_ROdata%Lev2a

! 14.9 Level 2b profile
! --------------------

  to_ROdata%Lev2b = from_ROdata%Lev2b

! 14.10 Level 2c profile
! --------------------

  to_ROdata%Lev2c = from_ROdata%Lev2c

! 14.11 Level 2d profile
! --------------------

  to_ROdata%Lev2d = from_ROdata%Lev2d

END SUBROUTINE ropp_io_assign_rotype


