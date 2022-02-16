! $Id: ropp_pp_preprocess_GRASRS.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_preprocess_GRASRS(ro_data, config, LCF)

!****s* Preprocessing/ropp_pp_preprocess_GRASRS *
!
! NAME
!    ropp_pp_preprocess_GRASRS - Mission-specific Level1a data preprocessing
!                                for GRAS Raw Sampling data
!
! SYNOPSIS
!    call ropp_pp_preprocess_GRASRS(ro_data, config, LCF)
!
! DESCRIPTION
!    Merge and upsample CL and RS data
!    1. Select CL and RS records by LCF flag
!    2. Generate merged time grid anchored at highest point of RS record
!    3. Interpolate CL and RS data to merged time grid
!    4. Merge data
!    5. Generate phase model and connecting phase
!    6. Interpolate data on merged time grid within small gaps
!    7. Restore phase variation
!
! INPUTS
!    type(ROprof)   :: ro_data      ! Radio occultation data strucuture
!    type(PPConfig) :: config       ! Configuration options
!    integer        :: LCF          ! Lost carrier flag
!
! OUTPUT
!    type(ROprof)   :: ro_data      ! Corrected radio occultation data
!    type(PPConfig) :: config       ! Configuration options
!    integer        :: LCF          ! Lost carrier flag
!
! NOTES
!   Requires ROprof data structure type, defined in ropp_io module. This
!   routine therefore requires that the ropp_io module is pre-installed before
!   compilation.
!
! REFERENCES
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils, ONLY: impact_parameter, ropp_MDFV, occ_point
  USE ropp_io, ONLY: ropp_io_init, ropp_io_free
  USE ropp_io_types, ONLY: ROprof, L1atype
! USE ropp_pp_preproc, not_this => ropp_pp_preprocess_GRASRS
  USE ropp_pp
  USE ropp_pp_spline
  USE ropp_pp_constants
  USE ropp_pp_types, ONLY: PPConfig

  IMPLICIT NONE

  TYPE(ROprof),          INTENT(inout) :: ro_data    ! Occultation data struct
  TYPE(PPconfig),        INTENT(inout) :: config     ! Configuration options
  INTEGER, DIMENSION(:), POINTER       :: LCF        ! Lost carrier flag

  INTEGER, PARAMETER                   :: np = 500   ! Dimension phase model
  INTEGER, PARAMETER                   :: nt = 300   ! RS grid interval
  INTEGER, PARAMETER                   :: nv = 5     ! Regression order
  REAL(wp), PARAMETER                  :: gap=0.04_wp ! Maximum data gap size

  TYPE(L1atype)                        :: mg_data    ! Merged data
  TYPE(L1atype)                        :: rs_data    ! Raw sampling data
  TYPE(L1atype)                        :: cl_data    ! Closed loop data
  INTEGER,  DIMENSION(:), ALLOCATABLE  :: mg_lcf     ! Merged LCF
  INTEGER,  DIMENSION(:), ALLOCATABLE  :: rs_lcf     ! Raw sampling LCF
  INTEGER,  DIMENSION(:), ALLOCATABLE  :: cl_lcf     ! Closed loop LCF

  REAL(wp), DIMENSION(:),  ALLOCATABLE :: mg_msis    ! MSIS phase model (merged)
  REAL(wp), DIMENSION(:),  ALLOCATABLE :: mg_impact  ! MSIS phase model (merged)
  REAL(wp), DIMENSION(:),  ALLOCATABLE :: mg_ds      ! phase deviation (merged)
  REAL(wp), DIMENSION(:),  ALLOCATABLE :: t_norm     ! Normalised time
  REAL(wp), DIMENSION(:,:),ALLOCATABLE :: KV         ! Regression matrix
  REAL(wp), DIMENSION(0:nv,3)          :: coeff_vleo ! Regression coeffs
  REAL(wp), DIMENSION(0:nv,3)          :: coeff_vgns ! Regression coeffs

  INTEGER                              :: ocd
  INTEGER                              :: i, j, n, n1, n_cl, n_rs
  INTEGER                              :: icl_min, icl_max
  INTEGER                              :: irs_min, irs_max
  INTEGER                              :: ig1, ig2, igl
  INTEGER                              :: i_int
  LOGICAL                              :: mrg_Mode  ! Merging mode

  REAL(wp)                             :: p1, pN, sb
  REAL(wp)                             :: ts, t1, tgl, tmin, tmax

!-------------------------------------------------------------------------------
! 2. Retrieve lost carrier flag information from input data file
!-------------------------------------------------------------------------------

  IF (ASSOCIATED(ro_data%vlist%VlistD1d)) THEN
    LCF(:) = NINT(ro_data%vlist%VlistD1d%data(1:size(LCF)))
    DEALLOCATE(ro_data%vlist%VlistD1d%data)
  ENDIF

!-------------------------------------------------------------------------------
! 3. Data quality checks
!-------------------------------------------------------------------------------

  WHERE(ro_data%Lev1a%phase_L1(:) == ropp_MDFV)
    ro_data%Lev1a%phase_L1(:) = -1.0_wp
  ENDWHERE
  WHERE(ro_data%Lev1a%phase_L2(:) == ropp_MDFV)
    ro_data%Lev1a%phase_L2(:) = -1.0_wp
  ENDWHERE

!-------------------------------------------------------------------------------
! 4. Initialisation
!-------------------------------------------------------------------------------

  mrg_Mode = .TRUE.

  ! 4.1 Select CL data

  icl_min = SUM(MINLOC(ro_data%lev1a%dtime(:), MASK = .NOT. BTEST(LCF(:),0)))
  icl_max = SUM(MAXLOC(ro_data%lev1a%dtime(:), MASK = .NOT. BTEST(LCF(:),0)))

  n_cl = icl_max-icl_min+1
  CALL ropp_io_init(cl_data, n_cl)
  ALLOCATE(cl_lcf(n_cl))

  cl_data%dtime = ro_data%lev1a%dtime(icl_min:icl_max)
  DO j=1,3
    cl_data%r_gns(:,j) = ro_data%lev1a%r_gns(icl_min:icl_max,j)
    cl_data%v_gns(:,j) = ro_data%lev1a%v_gns(icl_min:icl_max,j)
    cl_data%r_leo(:,j) = ro_data%lev1a%r_leo(icl_min:icl_max,j)
    cl_data%v_leo(:,j) = ro_data%lev1a%v_leo(icl_min:icl_max,j)
  ENDDO
  cl_data%snr_L1ca = ro_data%lev1a%snr_L1ca(icl_min:icl_max)
  cl_data%snr_L1p  = ro_data%lev1a%snr_L1p(icl_min:icl_max)
  cl_data%snr_L2p  = ro_data%lev1a%snr_L2p(icl_min:icl_max)
  cl_data%phase_L1 = ro_data%lev1a%phase_L1(icl_min:icl_max)
  cl_data%phase_L2 = ro_data%lev1a%phase_L2(icl_min:icl_max)
  cl_lcf           = lcf(icl_min:icl_max)

  ! 4.2 Select RS data

  irs_min = SUM(MINLOC(ro_data%lev1a%dtime(:), MASK = BTEST(LCF(:),0)))
  irs_max = SUM(MAXLOC(ro_data%lev1a%dtime(:), MASK = BTEST(LCF(:),0)))

  IF (irs_min >= 1 .AND. irs_min <= ro_data%lev1a%npoints) THEN

    n_rs = irs_max-irs_min+1
    CALL ropp_io_init(rs_data, n_rs)
    ALLOCATE(rs_lcf(n_rs))
    rs_data%dtime = ro_data%lev1a%dtime(irs_min:irs_max)
    DO j=1,3
      rs_data%r_gns(:,j) = ro_data%lev1a%r_gns(irs_min:irs_max,j)
      rs_data%v_gns(:,j) = ro_data%lev1a%v_gns(irs_min:irs_max,j)
      rs_data%r_leo(:,j) = ro_data%lev1a%r_leo(irs_min:irs_max,j)
      rs_data%v_leo(:,j) = ro_data%lev1a%v_leo(irs_min:irs_max,j)
    ENDDO
    rs_data%snr_L1ca = ro_data%lev1a%snr_L1ca(irs_min:irs_max)
    rs_data%snr_L1p  = ro_data%lev1a%snr_L1p(irs_min:irs_max)
    rs_data%snr_L2p  = ro_data%lev1a%snr_L2p(irs_min:irs_max)
    rs_data%phase_L1 = ro_data%lev1a%phase_L1(irs_min:irs_max)
    rs_data%phase_L2 = ro_data%lev1a%phase_L2(irs_min:irs_max)
    rs_lcf           = lcf(irs_min:irs_max)

  ELSE

    n_rs     = 0
    mrg_Mode = .FALSE.     ! Closed loop data only

  ENDIF

!-------------------------------------------------------------------------------
! 5. Generate merged time grid
!-------------------------------------------------------------------------------

  ! 5.1 Determine occultation direction

  p1 = impact_parameter(cl_data%r_leo(1,:), cl_data%r_gns(1,:))
  pN = impact_parameter(cl_data%r_leo(n_cl,:), cl_data%r_gns(n_cl,:))
  ocd = NINT(SIGN(1.0_wp, pN - p1))

  ! 5.2 Determine reference points and grid time step

  IF (mrg_Mode) THEN      ! Merge RS+CL

    ts = MINVAL(ABS(rs_data%dtime(1+nt:n_rs)-rs_data%dtime(1:n_rs-nt)))/nt
    IF (ocd < 0) THEN
      tmin = MINVAL(cl_data%dtime(:))
      tmax = MAXVAL(rs_data%dtime(:))
      t1 = MINVAL(rs_data%dtime(:))
    ELSE
      tmin = MINVAL(rs_data%dtime(:))
      tmax = MAXVAL(cl_data%dtime(:))
      t1 = MAXVAL(rs_data%dtime(:))
    ENDIF

  ELSE

    ts = MINVAL(ABS(cl_data%dtime(2:n_cl)-cl_data%dtime(1:n_cl-1)))
    tmin = MINVAL(cl_data%dtime(:))
    tmax = MAXVAL(cl_data%dtime(:))
    IF (ocd < 0) THEN
      t1 = MINVAL(cl_data%dtime(:))
    ELSE
      t1 = MAXVAL(cl_data%dtime(:))
    ENDIF

  ENDIF

  ! 5.3 Determine anchor point index and grid dimension

  n1 = FLOOR((t1 - tmin)/ts) + 1
  n  = FLOOR((tmax - t1)/ts) + n1

  ! 5.4 Generate merged time grid
  CALL ropp_io_init(mg_data, n)
  ALLOCATE(mg_LCF(n))
  mg_data%reference_frame%r_leo = ro_data%lev1a%reference_frame%r_leo
  mg_data%reference_frame%r_gns = ro_data%lev1a%reference_frame%r_gns

  IF (ocd < 0) THEN
    DO i=1,n
      mg_data%dtime(i) = ((n-i)*t1 + (i-n1)*tmax)/(n-n1)
    ENDDO
  ELSE
    DO i=1,n
      mg_data%dtime(i) = ((n1-i)*tmin + (i-1)*t1)/(n1-1)
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! 6. Interpolate CL and RS data to merged time grid
!-------------------------------------------------------------------------------

  ! 6.1 Interpolate trajectories

  ALLOCATE(t_norm(ro_data%Lev1a%npoints))
  ALLOCATE(KV(ro_data%Lev1a%npoints, 0:nv))

  t_norm(:) = (ro_data%Lev1a%dtime(:) - ro_data%Lev1a%dtime(1))/   &
          (ro_data%Lev1a%dtime(ro_data%Lev1a%npoints) - ro_data%Lev1a%dtime(1))

  CALL ropp_pp_init_polynomial(t_norm, KV)

  ! 6.2 Perform regression on positions

  DO j=1,3
    CALL ropp_pp_regression(KV,ro_data%Lev1a%r_leo(:,j),coeff_vleo(:,j))
    CALL ropp_pp_regression(KV,ro_data%Lev1a%r_gns(:,j),coeff_vgns(:,j))
  ENDDO

  ! 6.3 Perform regression on residual positions to gain higher accuracy

  DO j=1,3
    CALL ropp_pp_residual_regression(KV,t_norm,  &
                                     ro_data%Lev1a%r_leo(:,j),coeff_vleo(:,j))
    CALL ropp_pp_residual_regression(KV,t_norm,  &
                                     ro_data%Lev1a%r_gns(:,j),coeff_vgns(:,j))
  ENDDO

  DEALLOCATE(t_norm)
  DEALLOCATE(KV)

  DO i=1,n
    CALL ropp_pp_interpolate_trajectory(ro_data%Lev1a%dtime,         &
                                        coeff_vleo, coeff_vgns,      &
                                        ro_data%georef%r_coc,        &
                                        mg_data%dtime(i),            &
                                        mg_data%r_leo(i,:),          &
                                        mg_data%v_leo(i,:),          &
                                        mg_data%r_gns(i,:),          &
                                        mg_data%v_gns(i,:))
  ENDDO

  ! 6.4 Interpolate LCF, phase and amplitude

  DO i=1,n

    IF (mrg_Mode) THEN   ! Merge LC + RS data

      i_int = ropp_pp_seek_index(rs_data%dtime, mg_data%dtime(i))

      IF (i_int == 0 .OR. i_int == rs_data%npoints) THEN

        ! RS data range outside current time period
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_LCF, mg_LCF(i))
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_data%snr_L1ca, mg_data%snr_L1ca(i))
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_data%snr_L1p, mg_data%snr_L1p(i))
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_data%snr_L2p, mg_data%snr_L2p(i))
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_data%phase_L1, mg_data%phase_L1(i))
        CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                 cl_data%phase_L2, mg_data%phase_L2(i))

      ELSE

        ! RS data range within current time period
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_LCF, mg_LCF(i))
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_data%snr_L1ca, mg_data%snr_L1ca(i))
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_data%snr_L1p, mg_data%snr_L1p(i))
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_data%snr_L2p, mg_data%snr_L2p(i))
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_data%phase_L1, mg_data%phase_L1(i))
        CALL ropp_pp_interpol(rs_data%dtime, mg_data%dtime(i),   &
                                 rs_data%phase_L2, mg_data%phase_L2(i))

      ENDIF

    ELSE      ! Use only CL data

      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),         &
                               cl_LCF, mg_LCF(i))
      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                               cl_data%snr_L1ca, mg_data%snr_L1ca(i))
      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                               cl_data%snr_L1p, mg_data%snr_L1p(i))
      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                               cl_data%snr_L2p, mg_data%snr_L2p(i))
      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                                cl_data%phase_L1, mg_data%phase_L1(i))
      CALL ropp_pp_interpol(cl_data%dtime, mg_data%dtime(i),   &
                               cl_data%phase_L2, mg_data%phase_L2(i))

    ENDIF


  ENDDO

!-------------------------------------------------------------------------------
! 7. Connect phase
!-------------------------------------------------------------------------------

  ! 7.1 Calculate MSIS excess phase on merged time grid
  CALL occ_point( mg_data%r_leo,  mg_data%r_gns,       &
                  ro_data%georef%lat,   ro_data%georef%lon,        &
                  ro_data%georef%r_coc, ro_data%georef%roc,        &
                  ro_data%georef%azimuth,                          &
                  ro_data%georef%undulation,                       &
                  config%egm96, config%corr_egm96)

  ALLOCATE(mg_msis(n))
  ALLOCATE(mg_impact(n))
  CALL ropp_pp_modelphase(ro_data%dtocc%month, ro_data%georef%lat,  &
                          ro_data%georef%lon,  mg_data%dtime,       &
                          mg_data%r_leo,       mg_data%r_gns,       &
                          ro_data%georef%r_coc, ro_data%georef%roc, &
                          mg_msis, mg_impact, config)
  DEALLOCATE(mg_impact)

  ! 7.2 Phase re-accumulation

  ALLOCATE(mg_ds(n))
  WHERE (.NOT. BTEST(mg_LCF(:),3))
    mg_ds(:) = (2.0_wp*Pi*f_L1/c_light)*(mg_data%phase_L1(:)-mg_msis(:))
  ELSEWHERE
    mg_ds(:) = 0.0_wp
  ENDWHERE

  CALL Accumulate_Phase(mg_DS)

!-------------------------------------------------------------------------------
! 8. Fill in data gaps
!-------------------------------------------------------------------------------

  ts = MINVAL(ABS(cl_data%dtime(2:n_cl) - cl_data%dtime(1:n_cl-1)))

  ig1 = 1

  DO i=2,n

    IF (BTEST(mg_LCF(i),3) .AND. .NOT. BTEST(mg_LCF(i-1),3)) THEN
      ig1 = i-1
    ENDIF

    IF (.NOT. BTEST(mg_LCF(i),3) .AND. BTEST(mg_LCF(i-1),3)) THEN

      ig2 = i
      igl = ig2-ig1
      tgl = mg_data%dtime(ig2) - mg_data%dtime(ig1)

      IF (tgl < gap) THEN     ! Fill in data gap

        DO j=ig1+1,ig2-1
          CALL ropp_pp_interpol(mg_data%dtime(ig1:ig2:igl), mg_data%dtime(j), &
                                mg_ds(ig1:ig2:igl), mg_ds(j))
          mg_LCF(j) = IBCLR(mg_LCF(j),3)
        ENDDO

      ENDIF

    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
! 9. Restore and normalise phase
!-------------------------------------------------------------------------------

  ! 9.1 Restore phase variation

  mg_data%phase_L1(:) = mg_msis(:) + (c_light/(2.0_wp*pi*f_L1))*mg_ds(:)

  ! 9.2 Normalise to zero minimum value

  sb = MINVAL(mg_data%phase_L1(:))
  mg_data%phase_L1(:) = mg_data%phase_L1(:) - sb

!-------------------------------------------------------------------------------
! 7. Setup merged output data
!-------------------------------------------------------------------------------

  CALL ropp_io_free(ro_data%Lev1a)
  DEALLOCATE(LCF)

  CALL ropp_io_init(ro_data%Lev1a, n)
  ALLOCATE(LCF(n))

  ro_data%Lev1a = mg_data
  LCF(:) = mg_LCF(:)

!-------------------------------------------------------------------------------
! 8. Re-compute occultation point
!-------------------------------------------------------------------------------

    CALL occ_point( ro_data%Lev1a%r_leo,  ro_data%Lev1a%r_gns,       &
                  ro_data%georef%lat,   ro_data%georef%lon,        &
                  ro_data%georef%r_coc, ro_data%georef%roc,        &
                  ro_data%georef%azimuth,                          &
                  ro_data%georef%undulation,                       &
                  config%egm96, config%corr_egm96)

!-------------------------------------------------------------------------------
! 9. Quality control
!-------------------------------------------------------------------------------

  WHERE(ro_data%Lev1a%phase_L1(:) == ropp_MDFV)
    ro_data%Lev1a%phase_L1(:) = -1.0_wp
  ENDWHERE
  WHERE(ro_data%Lev1a%phase_L2(:) == ropp_MDFV)
    ro_data%Lev1a%phase_L2(:) = -1.0_wp
  ENDWHERE

!-------------------------------------------------------------------------------
! 10. Clean up
!-------------------------------------------------------------------------------

  DEALLOCATE(mg_msis)
  CALL ropp_io_free(cl_data)
  CALL ropp_io_free(rs_data)

CONTAINS

!-------------------------------------------------------------------------------
! 11. Transform phase to accumulated phase
!-------------------------------------------------------------------------------

  SUBROUTINE Accumulate_Phase(Ph, Sign)   ! (Array of (accumulated) phase, dir)

! Method:
!   Sign = 0 or no Sign:
!      Adding +-2*Pi where phase jumps from
!      +-Pi to -+Pi,
!   Sign > 0:
!      Adding +2*Pi where phase jumps from
!      - to +
!   Sign < 0
!      Adding -2*Pi where phase jumps from
!      + to -

    ! 11.1 Declarations

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_constants, ONLY: pi
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(inout) :: Ph   ! Phase --> accumulated phase
    INTEGER, OPTIONAL,      INTENT(in)    :: Sign ! Phase change sign

    INTEGER  :: i     ! Array index
    INTEGER  :: PSign ! Phase change sign

    ! 11.2 Determine phase change sign

    IF (.NOT. PRESENT(Sign)) THEN
      PSign = 0
    ELSE
      PSign = Sign
    ENDIF

    ! 11.3 Accumulate phase

    IF (PSign == 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1)+pi, 2*pi) - pi
      ENDDO
    ELSEIF (PSign > 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1), 2*pi)
      ENDDO
    ELSEIF (PSign < 0) THEN
      DO i=2,SIZE(Ph)
        Ph(i) = Ph(i-1) + MODULO(Ph(i)-Ph(i-1)+2*pi, 2*pi) - 2*pi
      ENDDO
    ENDIF

  END SUBROUTINE Accumulate_Phase


END SUBROUTINE ropp_pp_preprocess_GRASRS
