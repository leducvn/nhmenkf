! $Id: ropp_io_shrink.f90 1959 2008-11-13 12:15:18Z frhl $

!****s* Initialisation/ropp_io_shrink *
!
! NAME
!    ropp_io_shrink - Copy subset of ROprof structure to itself
!
! SYNOPSIS
!    use ropp_io
!    type(ROprof) :: ro_data
!      ...
!    call ropp_io_shrink(ro_data, imin, imax)
!
!      - or -
!
!    call ropp_io_shrink(ro_data%Lev1a, imin, imax, stride)
!    call ropp_io_shrink(ro_data%Lev1b, imin, imax, stride)
!    call ropp_io_shrink(ro_data%Lev2a, imin, imax, stride)
!    call ropp_io_shrink(ro_data%Lev2b, imin, imax, stride)
!    call ropp_io_shrink(ro_data%Lev2c, imin, imax, stride)
!    call ropp_io_shrink(ro_data%Lev2d, imin, imax, stride)
!
! DESCRIPTION
!    This subroutine copies data to be kept from the ROprof structure to a
!    local temporary storage, reallocates the array and copies data back.
!    Only single chunks of data are copied (i.e. all data between imin and
!    imax inclusively are kept). All elements within a data Level are copied
!    between the same limits.
!
! INPUTS
!    ro_data   dtyp  RO data (derived type)
!    imin            Initial index
!    imax            Final index
!    stride          Sampling interval
!
! OUTPUT
!    ro_data   dtyp  RO data (derived type)
!
! SEE ALSO
!    ropp_io_types
!
! NOTES
!    The current interface only allows the shrinking of scalar structures.
!    To shrink an array of structures (e.g., as obtained from reading a multifile),
!    a loop over all elements of the array of structures is required:
!
!       do i = 1, size(ro_data)
!          call ropp_io_shrink(ro_data(i), imin, imax, stride)
!       enddo
!
! REFERENCES
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

SUBROUTINE ropp_io_shrink_l1atype(var, imin, imax, stride)

! 1.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l1atype
  USE ropp_io_types, ONLY: L1atype

  IMPLICIT NONE

  TYPE(L1atype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  TYPE(L1atype)                :: work
  INTEGER                      :: n

  IF (var%Npoints < imax) RETURN
  n = size(var%dtime(imin:imax:stride))

  IF (n < var%Npoints) THEN

! 1.2 Allocate new structure
! --------------------------

  CALL ropp_io_init(work, n)

! 1.3 Copy old structure to new structure
! ---------------------------------------

  CALL copy_alloc( var%dtime(imin:imax:stride)      , work%dtime      )
  CALL copy_alloc( var%snr_L1ca(imin:imax:stride)   , work%snr_L1ca   )
  CALL copy_alloc( var%snr_L1p(imin:imax:stride)    , work%snr_L1p    )
  CALL copy_alloc( var%snr_L2p(imin:imax:stride)    , work%snr_L2p    )
  CALL copy_alloc( var%phase_L1(imin:imax:stride)   , work%phase_L1   )
  CALL copy_alloc( var%phase_L2(imin:imax:stride)   , work%phase_L2   )
  CALL copy_alloc( var%phase_qual(imin:imax:stride) , work%phase_qual )
  CALL copy_alloc( var%r_gns(imin:imax:stride, :)    , work%r_gns     )
  CALL copy_alloc( var%v_gns(imin:imax:stride, :)    , work%v_gns     )
  CALL copy_alloc( var%r_leo(imin:imax:stride, :)    , work%r_leo     )
  CALL copy_alloc( var%v_leo(imin:imax:stride, :)    , work%v_leo     )

! 1.4 Copy work structure to output structure
! -------------------------------------------

  CALL ropp_io_free( var )
  CALL ropp_io_init( var, n )
  var = work

  ENDIF

END SUBROUTINE ropp_io_shrink_l1atype


!-------------------------------------------------------------------------------
! 2. L1btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_l1btype(var, imin, imax, stride)

! 2.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l1btype
  USE ropp_io_types, ONLY: L1btype

  IMPLICIT NONE

  TYPE(L1btype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  TYPE(L1btype)                :: work
  INTEGER                      :: n

  IF (var%Npoints < imax) RETURN
  n = size(var%lat_tp(imin:imax:stride))

  IF (n < var%Npoints) THEN

! 2.2 Allocate new structure
! --------------------------

  CALL ropp_io_init(work, n)

! 2.3 Copy old structure to new structure
! ---------------------------------------

  CALL copy_alloc( var%lat_tp(imin:imax:stride)           , work%lat_tp           )
  CALL copy_alloc( var%lon_tp(imin:imax:stride)           , work%lon_tp           )
  CALL copy_alloc( var%azimuth_tp(imin:imax:stride)       , work%azimuth_tp       )
  CALL copy_alloc( var%impact_L1(imin:imax:stride)        , work%impact_L1        )
  CALL copy_alloc( var%impact_L2 (imin:imax:stride)       , work%impact_L2        )
  CALL copy_alloc( var%impact(imin:imax:stride)           , work%impact           )
  CALL copy_alloc( var%impact_Opt(imin:imax:stride)       , work%impact_Opt       )
  CALL copy_alloc( var%bangle_L1(imin:imax:stride)        , work%bangle_L1        )
  CALL copy_alloc( var%bangle_L2(imin:imax:stride)        , work%bangle_L2        )
  CALL copy_alloc( var%bangle(imin:imax:stride)           , work%bangle           )
  CALL copy_alloc( var%bangle_Opt(imin:imax:stride)       , work%bangle_Opt       )
  CALL copy_alloc( var%bangle_L1_sigma(imin:imax:stride)  , work%bangle_L1_sigma  )
  CALL copy_alloc( var%bangle_L2_sigma(imin:imax:stride)  , work%bangle_L2_sigma  )
  CALL copy_alloc( var%bangle_sigma(imin:imax:stride)     , work%bangle_sigma     )
  CALL copy_alloc( var%bangle_Opt_sigma(imin:imax:stride) , work%bangle_Opt_sigma )
  CALL copy_alloc( var%bangle_L1_qual(imin:imax:stride)   , work%bangle_L1_qual   )
  CALL copy_alloc( var%bangle_L2_qual(imin:imax:stride)   , work%bangle_L2_qual   )
  CALL copy_alloc( var%bangle_qual(imin:imax:stride)      , work%bangle_qual      )
  CALL copy_alloc( var%bangle_Opt_qual(imin:imax:stride)  , work%bangle_Opt_qual  )

! 2.4 Copy work structure to output structure
! -------------------------------------------

  CALL ropp_io_free( var )
  CALL ropp_io_init( var, n )
  var = work

  ENDIF

END SUBROUTINE ropp_io_shrink_l1btype


!-------------------------------------------------------------------------------
! 3. L2atype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_l2atype(var, imin, imax, stride)

! 3.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l2atype
  USE ropp_io_types, ONLY: L2atype

  IMPLICIT NONE

  TYPE(L2atype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  TYPE(L2atype)                :: work
  INTEGER                      :: n

  IF (var%Npoints < imax) RETURN
  n = size(var%alt_refrac(imin:imax:stride))

  IF (n < var%Npoints) THEN

! 3.2 Allocate new structure
! --------------------------

  CALL ropp_io_init(work, n)

! 3.3 Copy old structure to new structure
! ---------------------------------------

  CALL copy_alloc( var%alt_refrac(imin:imax:stride)     , work%alt_refrac     )
  CALL copy_alloc( var%geop_refrac(imin:imax:stride)    , work%geop_refrac    )
  CALL copy_alloc( var%refrac(imin:imax:stride)         , work%refrac         )
  CALL copy_alloc( var%refrac_sigma(imin:imax:stride)   , work%refrac_sigma   )
  CALL copy_alloc( var%refrac_qual(imin:imax:stride)    , work%refrac_qual    )
  CALL copy_alloc( var%dry_temp(imin:imax:stride)       , work%dry_temp       )
  CALL copy_alloc( var%dry_temp_sigma(imin:imax:stride) , work%dry_temp_sigma )
  CALL copy_alloc( var%dry_temp_qual(imin:imax:stride)  , work%dry_temp_qual  )

! 3.4 Copy work structure to output structure
! -------------------------------------------

  CALL ropp_io_free( var )
  CALL ropp_io_init( var, n )
  var = work

  ENDIF

END SUBROUTINE ropp_io_shrink_l2atype


!-------------------------------------------------------------------------------
! 4. L2btype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_l2btype(var, imin, imax, stride)

! 4.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l2btype
  USE ropp_io_types, ONLY: L2btype

  IMPLICIT NONE

  TYPE(L2btype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  TYPE(L2btype)                :: work
  INTEGER                      :: n

  IF (var%Npoints < imax) RETURN
  n = size(var%geop(imin:imax:stride))

  IF (n < var%Npoints) THEN

! 4.2 Allocate new structure
! --------------------------

  CALL ropp_io_init(work, n)

! 4.3 Copy old structure to new structure
! ---------------------------------------

  CALL copy_alloc( var%geop(imin:imax:stride)        , work%geop        )
  CALL copy_alloc( var%geop_sigma(imin:imax:stride)  , work%geop_sigma  )
  CALL copy_alloc( var%press(imin:imax:stride)       , work%press       )
  CALL copy_alloc( var%press_sigma(imin:imax:stride) , work%press_sigma )
  CALL copy_alloc( var%temp(imin:imax:stride)        , work%temp        )
  CALL copy_alloc( var%temp_sigma(imin:imax:stride)  , work%temp_sigma  )
  CALL copy_alloc( var%shum(imin:imax:stride)        , work%shum        )
  CALL copy_alloc( var%shum_sigma(imin:imax:stride)  , work%shum_sigma  )
  CALL copy_alloc( var%meteo_qual(imin:imax:stride)  , work%meteo_qual  )

! 4.4 Copy work structure to output structure
! -------------------------------------------

  CALL ropp_io_free( var )
  CALL ropp_io_init( var, n )
  var = work

  ENDIF

END SUBROUTINE ropp_io_shrink_l2btype


!-------------------------------------------------------------------------------
! 5. L2ctype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_l2ctype(var, imin, imax, stride)

! 5.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l2ctype
  USE ropp_io_types, ONLY: L2ctype

  IMPLICIT NONE

  TYPE(L2ctype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  INTEGER :: Dummy

! 5.2 Cannot shrink scalar variables
! --------------------------------------------------

  var = var
  IF (imax > 0 .AND. imin > 0) THEN
    CONTINUE
    Dummy = stride
  ENDIF

END SUBROUTINE ropp_io_shrink_l2ctype


!-------------------------------------------------------------------------------
! 6. L2dtype
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_l2dtype(var, imin, imax, stride)

! 6.1 Declarations
! ----------------

  USE ropp_utils
! USE ropp_io,       not_this => ropp_io_shrink_l2dtype
  USE ropp_io_types, ONLY: L2dtype

  IMPLICIT NONE

  TYPE(L2dtype), INTENT(inout) :: var
  INTEGER,       INTENT(in)    :: imin
  INTEGER,       INTENT(in)    :: imax
  INTEGER,       INTENT(in)    :: stride

  TYPE(L2dtype)                :: work
  INTEGER                      :: n

  IF (var%Npoints < imax) RETURN
  n = size(var%level_coeff_a(imin:imax:stride))

  IF (n < var%Npoints) THEN

! 6.2 Allocate new structure
! --------------------------

  CALL ropp_io_init(work, n)

! 6.3 Copy old structure to new structure
! ---------------------------------------

  CALL copy_alloc( var%level_coeff_a(imin:imax:stride) , work%level_coeff_a )
  CALL copy_alloc( var%level_coeff_b(imin:imax:stride) , work%level_coeff_b )

! 6.4 Copy work structure to output structure
! -------------------------------------------

  CALL ropp_io_free( var )
  CALL ropp_io_init( var, n )
  var = work

  ENDIF

END SUBROUTINE ropp_io_shrink_l2dtype


!-------------------------------------------------------------------------------
! 7. Joint RO data type
!-------------------------------------------------------------------------------

SUBROUTINE ropp_io_shrink_rotype(ROdata, imin, imax, stride)

! 7.1 Declarations
! ----------------

! USE ropp_io,        not_this => ropp_io_shrink_rotype
  USE ropp_io_types,  ONLY: ROprof

  IMPLICIT NONE

  TYPE(ROprof), INTENT(inout) :: ROdata
  INTEGER,       INTENT(in)   :: imin
  INTEGER,       INTENT(in)   :: imax
  INTEGER,       INTENT(in)   :: stride


! 7.2 Level 1a profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev1a, imin, imax, stride)

! 7.3 Level 1b profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev1b, imin, imax, stride)

! 7.4 Level 2a profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev2a, imin, imax, stride)

! 7.5 Level 2b profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev2b, imin, imax, stride)

! 7.6 Level 2c profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev2c, imin, imax, stride)

! 7.7 Level 2d profile
! --------------------

  CALL ropp_io_shrink(ROdata%Lev2d, imin, imax, stride)

END SUBROUTINE ropp_io_shrink_rotype











