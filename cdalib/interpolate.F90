#ifndef R_KIND
#define R_KIND 4
#endif
#define DEBUG 0
#define TEST 0

!#-- 2005.01.17 NISHIJIMA start ---
!# define macoro "R_KIND2" and change declaration : real(4) -> real(R_KIND2)
!# to vary accuracy according to situation (routine or TL/AD code check).
#ifndef R_KIND2
! for routine
#  define R_KIND2 4
! for TL/AD code check
!#  define R_KIND2 8    
#endif
!#-- 2005.01.17 NISHIJIMA end ---

module interpolate
!  === Interpolation Package Module ===
  implicit none
  save
#if TEST == 1
  public
  private :: Miss
#else
  private
  public :: setGridInf_fromNuSDaS, setGridInf, getGridInf, &
          & PolarInterpolation3D, &
          & NPInterpolation3D, SPInterpolation3D, &
          & PolarInterpolation2D, &
          & PolarInterpolation3D_Wind, &
          & NPInterpolation3D_Wind, SPInterpolation3D_Wind, &
          & NPInterpolation2D_Wind, SPInterpolation2D_Wind, &
          & interpolation2D, interpolation3D, &
          & interpolation3D_T, interpolation3D_Z, &
          & convert_LatLon_to_Grid, convert_LatLon_to_GridPos, &
          & convert_GridPos_offset, convert_GridPos_4P, &
          & convert_Lev_to_GridPos, convert_Lev_to_GridPos_T, &
          & convert_Lev_to_GridPos_Z, &
          & makeGroundT, vertInterpolation, &
          & PosGrid, gsize, proj, distance, basepoint, stand_ll, &
          & glat, glon, &
          & FlagOverTop, FlagUnderGround, FlagRange, &
          & rad
#endif
!  --- Parameter ---
  real(R_KIND), parameter   :: Miss = -32768.d0
  real(R_KIND), parameter   :: pi = 3.141593d0, rad = pi / 180.d0
!  --- Flag of flag of under the ground or over top ---
  integer(4), parameter     :: FlagOverTop = 1, &
                             & FlagUnderGround = -1, &
                             & FlagRange = 0
!  --- Position in the grid ---
  type PosGrid
#if DEBUG == 2
    sequence
#endif
    real(R_KIND2)           :: px, py, rx, ry
    integer(4)              :: i1, i2, j1, j2         ! grid No.
    real(R_KIND2)           :: rz11, rz21, rz12, rz22
    integer(4)              :: k11, k21, k12, k22     ! Level No.
!     flag of under the ground or over top
    integer(4)              :: range11, range21, range12, range22
!     height from surfce that are used under the ground only
    real(R_KIND2)           :: dlev11, dlev21, dlev12, dlev22
  end type PosGrid
!  --- Grid Information ---
  integer(4)                :: gsize(3) ! GPV size imax, jmax, kmax
  character(len=4)          :: proj     ! map projection
  real(R_KIND2)             :: distance(2), basepoint(4), stand_ll(4)
  logical                   :: circulate
!  --- for GAUSSIAN Grid ---
  real(R_KIND2), pointer    :: glat(:,:), glon(:,:)    ! Lat/Lon at Grid Points
!------------------------------------------------------------------------
!  Interface
  interface interpolation2D
    module procedure interpolation2D_Gauss, interpolation2D_Standard
  end interface
  interface interpolation3D
    module procedure interpolation3D_Gauss, interpolation3D_Standard
  end interface
  interface interpolation3D_T
    module procedure interpolation3D_T_Gauss, interpolation3D_T_Standard
  end interface
  interface interpolation3D_Z
    module procedure interpolation3D_Z_Gauss, interpolation3D_Z_Standard
  end interface
!------------------------------------------------------------------------
!************************************************************************
contains
!************************************************************************
subroutine setGridInf_fromNuSDaS(type1, type2, type3)
!  === set gsize, proj, distance, basepoint, stand_ll ===
!  === set glat, glon, circulate for Global Model ===
  character(len=8), intent(in) :: type1
  character(len=4), intent(in) :: type2, type3
!
  include 'nusdas_fort.h'
  integer(4)       :: rt
!
  call nusdas_inq_def(type1, type2, type3, N_GRID_SIZE, gsize(1:2), 2, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_GRID_SIZE) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 1
  end if
  call nusdas_inq_def(type1, type2, type3, N_PROJECTION, proj, 1, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_PROJECTION) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 2
  end if
  call nusdas_inq_def(type1, type2, type3, N_GRID_DISTANCE, distance, 2, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_GRID_DISTANCE) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 3
  end if
  call nusdas_inq_def(type1, type2, type3, N_GRID_BASEPOINT, basepoint, 4, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_GRID_BASEPOINT) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 4
  end if
  call nusdas_inq_def(type1, type2, type3, N_STAND_LATLON, stand_ll, 4, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_STAND_LATLON) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 5
  end if
  call nusdas_inq_def(type1, type2, type3, N_PLANE_NUM, gsize(3), 1, rt)
  if (rt < 0) then
    print *, "nusdas_inq_def(N_PLANE_NUM) error!! rt = ", rt
    print *, "type1, type2, type3 = ", type1, type2, type3
    stop 6
  end if
  if (proj == "GS  ") then
    if (associated(glat)) then
      deallocate(glat) ; deallocate(glon)
    end if
    allocate(glat(gsize(1),gsize(2))) ; allocate(glon(gsize(1),gsize(2)))
    call gausltln(glat, glon, gsize(1), gsize(2))
    circulate = .true.
  else
    circulate = .false.
  end if
#if DEBUG == 1
  print *, "in setGridInf_fromNuSDaS:", type1, type2, type3
  print *, "gsize:", gsize
  print *, "proj:", proj
  print *, "distance:", distance
  print *, "basepoint:", basepoint
  print *, "stand_ll:", stand_ll
#endif
end subroutine setGridInf_fromNuSDaS
!************************************************************************
subroutine getGridInf(o_gsize, o_proj, o_distance, o_basepoint, o_stand_ll, &
                    & o_glat, o_glon)
!  === get gsize, proj, distance, basepoint, stand_ll ===
!  === get glat, glon for Global Model ===
!  === this subroutine must be called after setGridInf ===
  integer(4), intent(out), optional  :: o_gsize(3) ! GPV size imax, jmax, kmax
  character(len=4), intent(out), optional  :: o_proj     ! map projection
  real(R_KIND2), intent(out), optional     :: o_distance(2), &
                                            & o_basepoint(4), o_stand_ll(4)
  real(R_KIND2), pointer, optional         :: o_glat(:,:), o_glon(:,:)
!
  if (present(o_gsize)) o_gsize(1:3) = gsize(1:3)
  if (present(o_proj)) o_proj = proj
  if (present(o_distance)) o_distance(1:2) = distance(1:2)
  if (present(o_basepoint)) o_basepoint(1:4) = basepoint(1:4)
  if (present(o_stand_ll)) o_stand_ll(1:4) = stand_ll(1:4)
  if (present(o_glat)) o_glat => glat
  if (present(o_glon)) o_glon => glon
!
end subroutine getGridInf
!************************************************************************
subroutine setGridInf(im, jm, km, npro, dels, &
                    & slon, xi, xj, xlat, xlon, slata, slatb)
!  === set gsize, proj, distance, basepoint, stand_ll ===
!  === set glat, glon, circulate for Global Model ===
  integer(4), intent(in), optional       :: im, jm, km
  character(len=4), intent(in), optional :: npro
  real(R_KIND2), intent(in), optional    :: dels, slon, xi, xj, xlat, xlon, &
                                          & slata, slatb
!
  if (present(im)) gsize(1) = im
  if (present(jm)) gsize(2) = jm
  if (present(km)) gsize(3) = km
  if (present(npro)) proj = npro
  if (present(dels)) distance(1:2) = dels
  if (present(xi)) basepoint(1) = xi
  if (present(xj)) basepoint(2) = xj
  if (present(xlat)) basepoint(3) = xlat
  if (present(xlon)) basepoint(4) = xlon
  if (present(slon)) then; stand_ll(2) = slon; stand_ll(4) = slon; end if
  if (present(slata)) stand_ll(1) = slata
  if (present(slatb)) stand_ll(3) = slatb
  if (present(im) .or. present(jm)) then
    if (proj == "GS  ") then
      if (associated(glat)) then
        deallocate(glat) ; deallocate(glon)
      end if
      allocate(glat(gsize(1),gsize(2))) ; allocate(glon(gsize(1),gsize(2)))
      call gausltln(glat, glon, gsize(1), gsize(2))
    end if
  end if
  if (present(npro)) then
    if (proj == "GS  ") then
      circulate = .true.
    else
      circulate = .false.
    end if
  end if
end subroutine setGridInf
!************************************************************************
subroutine PolarInterpolation3D(NValue, SValue, gpv)
  real(R_KIND), intent(out) :: NValue(:,:), SValue(:,:)
  real(R_KIND), intent(in)  :: gpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call PolarInterpolation2D(NValue(1:gsize(1),k), &
                              & gpv(1:gsize(1),1,k))
      call PolarInterpolation2D(SValue(1:gsize(1),k), &
                              & gpv(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine PolarInterpolation3D
!************************************************************************
subroutine NPInterpolation3D(NValue, gpv)
  real(R_KIND), intent(out) :: NValue(:,:)
  real(R_KIND), intent(in)  :: gpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call PolarInterpolation2D(NValue(1:gsize(1),k), &
                              & gpv(1:gsize(1),1,k))
    end do
  end if
end subroutine NPInterpolation3D
!************************************************************************
subroutine SPInterpolation3D(SValue, gpv)
  real(R_KIND), intent(out) :: SValue(:,:)
  real(R_KIND), intent(in)  :: gpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call PolarInterpolation2D(SValue(1:gsize(1),k), &
                              & gpv(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine SPInterpolation3D
!************************************************************************
subroutine PolarInterpolation2D(PValue, gpv)
!  === get interpolated vaule at the North/South Pole ===
!      for GAUSSIAN GRID only
  real(R_KIND), intent(out) :: PValue(:)
  real(R_KIND), intent(in)  :: gpv(:)
!
  if (proj == "GS  ") then
!   avarage of gpv near the pole
    PValue(1:gsize(1)) = sum(gpv(1:gsize(1))) / gsize(1)
  end if
end subroutine PolarInterpolation2D
!************************************************************************
subroutine PolarInterpolation3D_Wind(NValueU, NValueV, SValueU, SValueV, &
                                  & gpvU, gpvV)
  real(R_KIND), intent(out) :: NValueU(:,:), NValueV(:,:), &
                             & SValueU(:,:), SValueV(:,:)
  real(R_KIND), intent(in)  :: gpvU(:,:,:), gpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call NPInterpolation2D_Wind(NValueU(1:gsize(1),k), &
                               & NValueV(1:gsize(1),k), &
                               & gpvU(1:gsize(1),1,k), &
                               & gpvV(1:gsize(1),1,k))
      call SPInterpolation2D_Wind(SValueU(1:gsize(1),k), &
                               & SValueV(1:gsize(1),k), &
                               & gpvU(1:gsize(1),gsize(2),k), &
                               & gpvV(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine PolarInterpolation3D_Wind
!************************************************************************
subroutine NPInterpolation3D_Wind(NValueU, NValueV, gpvU, gpvV)
  real(R_KIND), intent(out) :: NValueU(:,:), NValueV(:,:)
  real(R_KIND), intent(in)  :: gpvU(:,:,:), gpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call NPInterpolation2D_Wind(NValueU(1:gsize(1),k), &
                               & NValueV(1:gsize(1),k), &
                               & gpvU(1:gsize(1),1,k), &
                               & gpvV(1:gsize(1),1,k))
    end do
  end if
end subroutine NPInterpolation3D_Wind
!************************************************************************
subroutine SPInterpolation3D_Wind(SValueU, SValueV, gpvU, gpvV)
  real(R_KIND), intent(out) :: SValueU(:,:), SValueV(:,:)
  real(R_KIND), intent(in)  :: gpvU(:,:,:), gpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call SPInterpolation2D_Wind(SValueU(1:gsize(1),k), &
                               & SValueV(1:gsize(1),k), &
                               & gpvU(1:gsize(1),gsize(2),k), &
                               & gpvV(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine SPInterpolation3D_Wind
!************************************************************************
subroutine NPInterpolation2D_Wind(NValueU, NValueV, gpvU, gpvV)
!  === get interpolated vaule (Wind) at the North Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(out) :: NValueU(:), NValueV(:)
  real(R_KIND), intent(in)  :: gpvU(:), gpvV(:)
!
  real(R_KIND) :: NPU, NPV
!
  if (proj == "GS  ") then
!  --- Œo“x0“x‚Ì•ûŒü‚ðŠî€‚É‚µ‚Ä‹É“_‚ÌŽü‚è‚Ì•—(U,V)‚ð
!  --- Œo“x•ª‚¾‚¯‰ñ“](–k‹É),‹t‰ñ“](“ì‹É)‚³‚¹‚é
!        U = ucosA-vsinA, V =  usinA+vcosA    (‰ñ“])
!        U = ucosA+vsinA, V = -usinA+vcosA    (‹t‰ñ“])
!  --- ‚»‚ê‚ç‚Ì•½‹Ï‚ð‚Æ‚Á‚ÄŒo“x0“x‚Ì•ûŒü‚ðŠî€‚Æ‚·‚é‹É“_‚Ì•—‚Æ‚·‚é
    NPU = sum( gpvU(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) &
           & - gpvV(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
    NPV = sum( gpvU(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) &
           & + gpvV(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
!  --- ‹É“_‚Ì•—‚ðŠeŒo“x‚Ì•ûŒü‚É‡‚¤‚æ‚¤‚É‰ñ“]‚³‚¹‚é
!      ‹t‰ñ“](–k‹É),‰ñ“](“ì‹É)‚³‚¹‚é
    NValueU(1:gsize(1)) =  NPU * cos(glon(1:gsize(1),1) * rad) &
                       & + NPV * sin(glon(1:gsize(1),1) * rad)
    NValueV(1:gsize(1)) = -NPU * sin(glon(1:gsize(1),1) * rad) &
                       & + NPV * cos(glon(1:gsize(1),1) * rad)
  end if
end subroutine NPInterpolation2D_Wind
!************************************************************************
subroutine SPInterpolation2D_Wind(SValueU, SValueV, gpvU, gpvV)
!  === get interpolated vaule (Wind) at the South Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(out) :: SValueU(:), SValueV(:)
  real(R_KIND), intent(in)  :: gpvU(:), gpvV(:)
!
  real(R_KIND) :: SPU, SPV
!
  if (proj == "GS  ") then
!  --- Œo“x0“x‚Ì•ûŒü‚ðŠî€‚É‚µ‚Ä‹É“_‚ÌŽü‚è‚Ì•—(U,V)‚ð
!  --- Œo“x•ª‚¾‚¯‰ñ“](–k‹É),‹t‰ñ“](“ì‹É)‚³‚¹‚é
!        U = ucosA-vsinA, V =  usinA+vcosA    (‰ñ“])
!        U = ucosA+vsinA, V = -usinA+vcosA    (‹t‰ñ“])
!  --- ‚»‚ê‚ç‚Ì•½‹Ï‚ð‚Æ‚Á‚ÄŒo“x0“x‚Ì•ûŒü‚ðŠî€‚Æ‚·‚é‹É“_‚Ì•—‚Æ‚·‚é
    SPU = sum( gpvU(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) &
           & + gpvV(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
    SPV = sum(-gpvU(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) &
           & + gpvV(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
!  --- ‹É“_‚Ì•—‚ðŠeŒo“x‚Ì•ûŒü‚É‡‚¤‚æ‚¤‚É‰ñ“]‚³‚¹‚é
!      ‹t‰ñ“](–k‹É),‰ñ“](“ì‹É)‚³‚¹‚é
    SValueU(1:gsize(1)) =  SPU * cos(glon(1:gsize(1),gsize(2)) * rad) &
                       & - SPV * sin(glon(1:gsize(1),gsize(2)) * rad)
    SValueV(1:gsize(1)) =  SPU * sin(glon(1:gsize(1),gsize(2)) * rad) &
                       & + SPV * cos(glon(1:gsize(1),gsize(2)) * rad)
  end if
end subroutine SPInterpolation2D_Wind
!************************************************************************
subroutine interpolation2D_Gauss(value, pg, gpv, NPgpv, SPgpv)
  real(R_KIND), intent(out) :: value
  type(PosGrid), intent(in) :: pg
  real(R_KIND), intent(in)  :: gpv(:,:)
  real(R_KIND), intent(in)  :: NPgpv(:), SPgpv(:)
!  --- near the North Pole ---
  if (pg%j1 == 0 .and. proj == "GS  ") then
    call bilinearInterpolation(value, &
    & NPgpv(pg%i1), NPgpv(pg%i2), gpv(pg%i1,pg%j2), gpv(pg%i2,pg%j2), &
    & pg%rx, pg%ry)
!  --- near the South Pole ---
  else if (pg%j1 == gsize(2) .and. proj == "GS  ") then
    call bilinearInterpolation(value, &
    & gpv(pg%i1,pg%j1), gpv(pg%i2,pg%j1), SPgpv(pg%i1), SPgpv(pg%i2), &
    & pg%rx, pg%ry)
!  --- other area ---
  else if (pg%j1 > 0) then
    call bilinearInterpolation(value, &
    & gpv(pg%i1,pg%j1), gpv(pg%i2,pg%j1), gpv(pg%i1,pg%j2), gpv(pg%i2,pg%j2), &
    & pg%rx, pg%ry)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation2D_Gauss
!************************************************************************
subroutine interpolation2D_Standard(value, pg, gpv)
  real(R_KIND), intent(out) :: value
  type(PosGrid), intent(in) :: pg
  real(R_KIND), intent(in)  :: gpv(:,:)
!
  if (pg%j1 > 0) then
    call bilinearInterpolation(value, &
    & gpv(pg%i1,pg%j1), gpv(pg%i2,pg%j1), gpv(pg%i1,pg%j2), gpv(pg%i2,pg%j2), &
    & pg%rx, pg%ry)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation2D_Standard
!************************************************************************
subroutine interpolation3D_Gauss &
         & (value, pg, lev, gpv, lgpv, &
          & NPgpv, NPlgpv, SPgpv, SPlgpv)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
  real(R_KIND), target           :: NPgpv(:,:), NPlgpv(:,:)
  real(R_KIND), target           :: SPgpv(:,:), SPlgpv(:,:)
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  if (gsize(3) <= 1) then
    call interpolation2D_Gauss(value, pg, gpv(:,:,1), NPgpv(:,1), SPgpv(:,1))
    return
  end if
!  --- near the North Pole ---
  if (pg%j1 == 0 .and. proj == "GS  ") then
    p11 => NPlgpv(pg%i1,:)
    p21 => NPlgpv(pg%i2,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos(pg, p11, p21, p12, p22, lev)
    p11 => NPgpv(pg%i1,:)
    p21 => NPgpv(pg%i2,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- near the South Pole ---
  else if (pg%j1 == gsize(2) .and. proj == "GS  ") then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => SPlgpv(pg%i1,:)
    p22 => SPlgpv(pg%i2,:)
    call convert_Lev_to_GridPos(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => SPgpv(pg%i1,:)
    p22 => SPgpv(pg%i2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- other area ---
  else if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_Gauss
!************************************************************************
subroutine interpolation3D_Standard &
         & (value, pg, lev, gpv, lgpv)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  if (gsize(3) <= 1) then
    call interpolation2D_Standard(value, pg, gpv(:,:,1))
    return
  end if
  if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_Standard
!************************************************************************
subroutine interpolation3D_T_Gauss &
         & (value, pg, lev, gpv, lgpv, lElement, &
          & NPgpv, NPlgpv, SPgpv, SPlgpv)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
  character(len=6), intent(in)   :: lElement  ! element name of level
  real(R_KIND), target           :: NPgpv(:,:), NPlgpv(:,:)
  real(R_KIND), target           :: SPgpv(:,:), SPlgpv(:,:)
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- temperature under the ground ---
  real(R_KIND), target      :: t11(2), t21(2), t12(2), t22(2)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  if (gsize(3) <= 1) then
    call interpolation2D_Gauss(value, pg, gpv(:,:,1), NPgpv(:,1), SPgpv(:,1))
    return
  end if
!  --- near the North Pole ---
  if (pg%j1 == 0 .and. proj == "GS  ") then
    p11 => NPlgpv(pg%i1,:)
    p21 => NPlgpv(pg%i2,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_T(pg, p11, p21, p12, p22, lev)
    if (pg%range11 /= FlagUnderGround) then
      p11 => NPgpv(pg%i1,:)
    else
      call makeGroundT(t11(1), NPgpv(pg%i1,1), pg%dlev11, lElement)
      p11 => t11
    end if
    if (pg%range21 /= FlagUnderGround) then
      p21 => NPgpv(pg%i2,:)
    else
      call makeGroundT(t21(1), NPgpv(pg%i2,1), pg%dlev21, lElement)
      p21 => t21
    end if
    if (pg%range12 /= FlagUnderGround) then
      p12 => gpv(pg%i1,pg%j2,:)
    else
      call makeGroundT(t12(1), gpv(pg%i1,pg%j2,1), pg%dlev12, lElement)
      p12 => t12
    end if
    if (pg%range22 /= FlagUnderGround) then
      p22 => gpv(pg%i2,pg%j2,:)
    else
      call makeGroundT(t22(1), gpv(pg%i2,pg%j2,1), pg%dlev22, lElement)
      p12 => t12
    end if
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- near the South Pole ---
  else if (pg%j1 == gsize(2) .and. proj == "GS  ") then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => SPlgpv(pg%i1,:)
    p22 => SPlgpv(pg%i2,:)
    call convert_Lev_to_GridPos_T(pg, p11, p21, p12, p22, lev)
    if (pg%range11 /= FlagUnderGround) then
      p11 => gpv(pg%i1,pg%j1,:)
    else
      call makeGroundT(t11(1), gpv(pg%i1,pg%j1,1), pg%dlev11, lElement)
      p11 => t11
    end if
    if (pg%range21 /= FlagUnderGround) then
      p21 => gpv(pg%i2,pg%j1,:)
    else
      call makeGroundT(t21(1), gpv(pg%i2,pg%j1,1), pg%dlev21, lElement)
      p21 => t21
    end if
    if (pg%range12 /= FlagUnderGround) then
      p12 => SPgpv(pg%i1,:)
    else
      call makeGroundT(t12(1), SPgpv(pg%i1,1), pg%dlev12, lElement)
      p12 => t12
    end if
    if (pg%range22 /= FlagUnderGround) then
      p22 => SPgpv(pg%i2,:)
    else
      call makeGroundT(t22(1), SPgpv(pg%i2,1), pg%dlev22, lElement)
      p12 => t12
    end if
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- other area ---
  else if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_T(pg, p11, p21, p12, p22, lev)
    if (pg%range11 /= FlagUnderGround) then
      p11 => gpv(pg%i1,pg%j1,:)
    else
      call makeGroundT(t11(1), gpv(pg%i1,pg%j1,1), pg%dlev11, lElement)
      p11 => t11
    end if
    if (pg%range21 /= FlagUnderGround) then
      p21 => gpv(pg%i2,pg%j1,:)
    else
      call makeGroundT(t21(1), gpv(pg%i2,pg%j1,1), pg%dlev21, lElement)
      p21 => t21
    end if
    if (pg%range12 /= FlagUnderGround) then
      p12 => gpv(pg%i1,pg%j2,:)
    else
      call makeGroundT(t12(1), gpv(pg%i1,pg%j2,1), pg%dlev12, lElement)
      p12 => t12
    end if
    if (pg%range22 /= FlagUnderGround) then
      p22 => gpv(pg%i2,pg%j2,:)
    else
      call makeGroundT(t22(1), gpv(pg%i2,pg%j2,1), pg%dlev22, lElement)
      p22 => t22
    end if
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_T_Gauss
!************************************************************************
subroutine interpolation3D_T_Standard &
         & (value, pg, lev, gpv, lgpv, lElement)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
  character(len=6), intent(in)   :: lElement  ! element name of level
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- temperature under the ground ---
  real(R_KIND), target      :: t11(2), t21(2), t12(2), t22(2)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  ! initialize for MPI auto
  t11(:) = 0.0d0; t12(:) = 0.0d0; t21(:) = 0.0d0; t22(:) = 0.0d0

  if (gsize(3) <= 1) then
    call interpolation2D_Standard(value, pg, gpv(:,:,1))
    return
  end if
  if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_T(pg, p11, p21, p12, p22, lev)
    if (pg%range11 /= FlagUnderGround) then
      p11 => gpv(pg%i1,pg%j1,:)
    else
      call makeGroundT(t11(1), gpv(pg%i1,pg%j1,1), pg%dlev11, lElement)
      p11 => t11
    end if
    if (pg%range21 /= FlagUnderGround) then
      p21 => gpv(pg%i2,pg%j1,:)
    else
      call makeGroundT(t21(1), gpv(pg%i2,pg%j1,1), pg%dlev21, lElement)
      p21 => t21
    end if
    if (pg%range12 /= FlagUnderGround) then
      p12 => gpv(pg%i1,pg%j2,:)
    else
      call makeGroundT(t12(1), gpv(pg%i1,pg%j2,1), pg%dlev12, lElement)
      p12 => t12
    end if
    if (pg%range22 /= FlagUnderGround) then
      p22 => gpv(pg%i2,pg%j2,:)
    else
      call makeGroundT(t22(1), gpv(pg%i2,pg%j2,1), pg%dlev22, lElement)
      p22 => t22
    end if
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_T_Standard
!************************************************************************
subroutine interpolation3D_Z_Gauss &
         & (value, pg, lev, gpv, lgpv, &
          & NPgpv, NPlgpv, SPgpv, SPlgpv)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
  real(R_KIND), target           :: NPgpv(:,:), NPlgpv(:,:)
  real(R_KIND), target           :: SPgpv(:,:), SPlgpv(:,:)
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  if (gsize(3) <= 1) then
    call interpolation2D_Gauss(value, pg, gpv(:,:,1), NPgpv(:,1), SPgpv(:,1))
    return
  end if
!  --- near the North Pole ---
  if (pg%j1 == 0 .and. proj == "GS  ") then
    p11 => NPlgpv(pg%i1,:)
    p21 => NPlgpv(pg%i2,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_Z(pg, p11, p21, p12, p22, lev)
    p11 => NPgpv(pg%i1,:)
    p21 => NPgpv(pg%i2,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- near the South Pole ---
  else if (pg%j1 == gsize(2) .and. proj == "GS  ") then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => SPlgpv(pg%i1,:)
    p22 => SPlgpv(pg%i2,:)
    call convert_Lev_to_GridPos_Z(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => SPgpv(pg%i1,:)
    p22 => SPgpv(pg%i2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- other area ---
  else if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_Z(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_Z_Gauss
!************************************************************************
subroutine interpolation3D_Z_Standard &
         & (value, pg, lev, gpv, lgpv)
  real(R_KIND), intent(out)      :: value
  type(PosGrid), intent(inout)   :: pg
  real(R_KIND2), intent(in)      :: lev
  real(R_KIND), target           :: gpv(:,:,:), lgpv(:,:,:)
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND), pointer     :: p11(:), p21(:), p12(:), p22(:)
!  --- in the case num level == 1 then as same as interpolate 2D ---
  if (gsize(3) <= 1) then
    call interpolation2D_Standard(value, pg, gpv(:,:,1))
    return
  end if
  if (pg%j1 > 0) then
    p11 => lgpv(pg%i1,pg%j1,:)
    p21 => lgpv(pg%i2,pg%j1,:)
    p12 => lgpv(pg%i1,pg%j2,:)
    p22 => lgpv(pg%i2,pg%j2,:)
    call convert_Lev_to_GridPos_Z(pg, p11, p21, p12, p22, lev)
    p11 => gpv(pg%i1,pg%j1,:)
    p21 => gpv(pg%i2,pg%j1,:)
    p12 => gpv(pg%i1,pg%j2,:)
    p22 => gpv(pg%i2,pg%j2,:)
    call isobaricInterpolation(value, p11, p21, p12, p22, pg)
!  --- out of range ---
  else
    value = Miss
  end if
end subroutine interpolation3D_Z_Standard
!************************************************************************
subroutine convert_LatLon_to_Grid(fi, fj, lat, lon, rev_j)
  real(R_KIND2), intent(out)         :: fi, fj
  real(R_KIND2), intent(in)          :: lat, lon
  integer(4), intent(in), optional :: rev_j
!
  type(PosGrid)                :: pg
!
  call convert_LatLon_to_GridPos(pg, lat, lon, rev_j)
  fi = pg%px
  fj = pg%py
end subroutine convert_LatLon_to_Grid
!************************************************************************
subroutine convert_LatLon_to_GridPos(pg, lat, lon, rev_j)
  type(PosGrid), intent(inout) :: pg
  real(R_KIND2), intent(in)          :: lat, lon
  integer(4), intent(in), optional :: rev_j
!
  real(R_KIND2)                      :: wk_lat, wk_lon
#if DEBUG == 9 || DEBUG == 8
  integer(4)                   :: ii, jj
#endif
!
  wk_lat = lat
  if (lon >= 0) then
    wk_lon = lon
  else
    wk_lon = lon + 360.
  end if
  if (proj == "GS  ") then
    call lltoij_gauss(pg%px, pg%py, pg%rx, pg%ry, &
              & wk_lat, wk_lon, glat, gsize(1), gsize(2))
#if DEBUG == 9
    call posdet(wk_lat, wk_lon, glat, gsize(1), gsize(2), &
              & ii, jj, pg%rx, pg%ry)
    pg%px = ii + pg%rx
    pg%py = jj + pg%ry
#endif
#if DEBUG == 8
    pg%px = wk_lon * gsize(1) / 360. + 1.
    pg%rx = pg%px - int(pg%px)
    if (wk_lat <= glat(1,1) .and. wk_lat > glat(1,gsize(2))) then
      do j = 2, gsize(2)
        if (wk_lat > glat(1,j)) then
          pg%py = (j - 1) + (glat(1,j-1) - wk_lat) / (glat(1,j-1) - glat(1,j))
          exit
        end if
      end do
      pg%ry = pg%py - int(pg%py)
    else if (wk_lat > glat(1,1)) then  ! near the north pole
      pg%py = (90. - wk_lat) / (90. - glat(1,1))
      if (pg%py < 0.) then
        pg%py = 0.
      end if
      pg%ry = pg%py - int(pg%py)
    else if (wk_lat <= glat(1,gsize(2))) then ! near the south pole
      pg%py = gsize(2) &
          & + (glat(1,gsize(2)) - wk_lat) / (glat(1,gsize(2)) - (-90.))
      if (pg%py >= gsize(2) + 1) then ! at the south pole
        pg%py = gsize(2)
        pg%ry = 1.
      else
        pg%ry = pg%py - int(pg%py)
      end if
    end if
#endif
  else if (proj == "LL  ") then
    pg%px = basepoint(1) + (wk_lon - basepoint(4)) / distance(1)
    pg%py = basepoint(2) + (wk_lat - basepoint(3)) / distance(1)
    pg%rx = pg%px - int(pg%px)
    pg%ry = pg%py - int(pg%py)
  else
    call rltln(pg%px, pg%py, wk_lat, wk_lon, &
              & proj, distance(1), stand_ll(2), &
              & basepoint(1), basepoint(2), basepoint(3), basepoint(4))
    pg%rx = pg%px - int(pg%px)
    if (present(rev_j)) pg%py = gsize(2) - pg%py + 1.
    pg%ry = pg%py - int(pg%py)
  end if
  if (circulate .and. pg%px < 1) then
    pg%px = pg%px + gsize(1)
  end if
!
  call getGridNo(pg)
end subroutine convert_LatLon_to_GridPos
!************************************************************************
subroutine convert_GridPos_offset(pg, i_off, j_off)
  type(PosGrid), intent(inout)  :: pg
  integer, intent(in), optional :: i_off, j_off
!
  if (present(i_off)) then
    pg%i1 = pg%i1 - i_off
    pg%i2 = pg%i2 - i_off
  end if
  if (present(j_off)) then
    pg%j1 = pg%j1 - j_off
    pg%j2 = pg%j2 - j_off
  end if
end subroutine convert_GridPos_offset
!************************************************************************
subroutine convert_GridPos_4P(pg)
  type(PosGrid), intent(inout)  :: pg
!
  pg%i1 = 1
  pg%i2 = 2
  pg%j1 = 1
  pg%j2 = 2
end subroutine convert_GridPos_4P
!************************************************************************
subroutine convert_Lev_to_GridPos(pg, p11, p21, p12, p22, lev)
  type(PosGrid), intent(inout) :: pg
  real(R_KIND), intent(in)     :: p11(:), p21(:), p12(:), p22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call getLevNo(pg%k11, pg%range11, p11, lev)
  call getLevNo(pg%k21, pg%range21, p21, lev)
  call getLevNo(pg%k12, pg%range12, p12, lev)
  call getLevNo(pg%k22, pg%range22, p22, lev)
  call getZrate(pg%rz11, p11, lev, pg%k11, pg%range11)
  call getZrate(pg%rz21, p21, lev, pg%k21, pg%range21)
  call getZrate(pg%rz12, p12, lev, pg%k12, pg%range12)
  call getZrate(pg%rz22, p22, lev, pg%k22, pg%range22)
end subroutine convert_Lev_to_GridPos
!************************************************************************
subroutine convert_Lev_to_GridPos_T(pg, p11, p21, p12, p22, lev)
  type(PosGrid), intent(inout) :: pg
  real(R_KIND), intent(in)     :: p11(:), p21(:), p12(:), p22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call getLevNo(pg%k11, pg%range11, p11, lev)
  call getLevNo(pg%k21, pg%range21, p21, lev)
  call getLevNo(pg%k12, pg%range12, p12, lev)
  call getLevNo(pg%k22, pg%range22, p22, lev)
  call getZrate_T(pg%rz11, pg%dlev11, p11, lev, pg%k11, pg%range11)
  call getZrate_T(pg%rz21, pg%dlev21, p21, lev, pg%k21, pg%range21)
  call getZrate_T(pg%rz12, pg%dlev12, p12, lev, pg%k12, pg%range12)
  call getZrate_T(pg%rz22, pg%dlev22, p22, lev, pg%k22, pg%range22)
end subroutine convert_Lev_to_GridPos_T
!************************************************************************
subroutine convert_Lev_to_GridPos_Z(pg, p11, p21, p12, p22, lev)
  type(PosGrid), intent(inout) :: pg
  real(R_KIND), intent(in)     :: p11(:), p21(:), p12(:), p22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call getLevNo(pg%k11, pg%range11, p11, lev)
  call getLevNo(pg%k21, pg%range21, p21, lev)
  call getLevNo(pg%k12, pg%range12, p12, lev)
  call getLevNo(pg%k22, pg%range22, p22, lev)
  call getZrate_Z(pg%rz11, p11, lev, pg%k11, pg%range11)
  call getZrate_Z(pg%rz21, p21, lev, pg%k21, pg%range21)
  call getZrate_Z(pg%rz12, p12, lev, pg%k12, pg%range12)
  call getZrate_Z(pg%rz22, p22, lev, pg%k22, pg%range22)
end subroutine convert_Lev_to_GridPos_Z
!************************************************************************
subroutine getGridNo(pg)
  type(PosGrid), intent(inout) :: pg
!
  pg%i1 = int(pg%px) ; pg%j1 = int(pg%py)
  if (pg%px /= gsize(1)) then
    pg%i2 = pg%i1 + 1
  else
    pg%i2 = pg%i1
  end if
  if (pg%py /= gsize(2)) then
    pg%j2 = pg%j1 + 1
  else
    pg%j2 = pg%j1
  end if
  if (circulate .and. pg%i2 > gsize(1)) pg%i2 = 1
!  --- set j1 = -1 in the case out of range ---
  if (proj == "GS  ") then
    if (pg%j1 < 0 .or. pg%j2 > gsize(2)+1 .or. &
      & pg%i1 < 1 .or. pg%i2 > gsize(1)) pg%j1 = -1
  else
    if (pg%j1 < 1 .or. pg%j2 > gsize(2) .or. &
      & pg%i1 < 1 .or. pg%i2 > gsize(1)) pg%j1 = -1
  end if
end subroutine getGridNo
!************************************************************************
subroutine getLevNo(k, flag, levs, lev)
  integer(4), intent(out)   :: k
  integer(4), intent(out)   :: flag ! flag of under the ground or over top
  real(R_KIND), intent(in)  :: levs(:)
  real(R_KIND2), intent(in) :: lev
!
  integer(4) :: nk
!
  k = -1
  do nk = 1, gsize(3) - 1
    if ((lev - levs(nk)) * (lev - levs(nk+1)) <= 0) then
      k = nk
      flag = FlagRange
      exit
    end if
  end do
  if (k < 0) then
!  --- over the top of model level ---
    if ((lev - levs(1)) * (levs(2) - levs(1)) > 0) then
      k = gsize(3) - 1
      flag = FlagOverTop
    else
      k = 1
      flag = FlagUnderGround
    end if
  end if
end subroutine getLevNo
!************************************************************************
subroutine getZrate(z, levs, lev, k, flag)
  real(R_KIND2), intent(out):: z
  real(R_KIND), intent(in)  :: levs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
  if (flag == FlagOverTop) then
!  --- over the top of model level ---
    z = 1.            ! same value of the top of the model
  else if (flag == FlagUnderGround) then
!  --- under the ground ---
    z = 0.            ! same value of the ground
  else
    z = (lev - levs(k)) / (levs(k+1) - levs(k))
  end if
end subroutine getZrate
!************************************************************************
subroutine getZrate_T(z, dlev, levs, lev, k, flag)
  real(R_KIND2), intent(out):: z
  real(R_KIND2), intent(out):: dlev ! height from surfce
  real(R_KIND), intent(in)  :: levs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
  if (flag == FlagUnderGround) then
    z = 0.            ! exception, continue to makeGroundT
    dlev = lev - levs(1) ! difference of levels
  else ! including "if (flag == FlagOverTop)", exterpolation
    z = (lev - levs(k)) / (levs(k+1) - levs(k))
  end if
end subroutine getZrate_T
!************************************************************************
subroutine getZrate_Z(z, levs, lev, k, flag)
  real(R_KIND2), intent(out):: z
  real(R_KIND), intent(in)  :: levs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
! including "if (flag == FlagOverTop .and. FdlagUnderGround)", exterpolation
  z = (lev - levs(k)) / (levs(k+1) - levs(k))
end subroutine getZrate_Z
!************************************************************************
subroutine makeGroundT(t, t1, dlev, lElement)
  real(R_KIND), intent(out) :: t    ! temperature under the ground
  real(R_KIND), intent(in)  :: t1   ! temperature at the surface
  real(R_KIND2), intent(in) :: dlev ! difference of levels
  character(len=6), intent(in)   :: lElement  ! element name of level
!
  real(R_KIND), parameter :: &
  & gamma = 0.005d0, &  ! temperature lapse rate (K/-m)
  & rd = 287.05d0, &    ! dry gas constant
  & grv = 9.80665d0     ! gravitational acceleration
!
  if (lElement == "logP  ") then      ! dlev > 1
!   Nomary t1 must be virtual temperature, but supposed as q = 0
    t = exp(dlev*gamma*rd/grv) * t1
  else if (lElement == "Z     ") then ! dlev < 0
    t = t1 - gamma * dlev
  else
    t = t1
  end if
end subroutine makeGroundT
!************************************************************************
subroutine bilinearInterpolation(v, v11, v21, v12, v22, x, y)
  real(R_KIND), intent(out) :: v
  real(R_KIND), intent(in)  :: v11, v21, v12, v22
  real(R_KIND2), intent(in) :: x, y
!
  v = (x * v21 + (1. - x) * v11) * (1. - y) &
  & + (x * v22 + (1. - x) * v12) * y
end subroutine bilinearInterpolation
!************************************************************************
subroutine isobaricInterpolation(value, p11, p21, p12, p22, pg)
  real(R_KIND), intent(out) :: value
  real(R_KIND), intent(in)  :: p11(:), p21(:), p12(:), p22(:)
  type(PosGrid), intent(in) :: pg
!
  real(R_KIND)              :: v11, v21, v12, v22
!
  call vertInterpolation(v11, p11, pg%rz11, pg%k11)
  call vertInterpolation(v21, p21, pg%rz21, pg%k21)
  call vertInterpolation(v12, p12, pg%rz12, pg%k12)
  call vertInterpolation(v22, p22, pg%rz22, pg%k22)
#if DEBUG == 4
  print *, "in isobaricInterpolation: v11, v21, v12, v22:", v11, v21, v12, v22
#endif
  call bilinearInterpolation(value, v11, v21, v12, v22, &
                           & pg%rx, pg%ry)
!
end subroutine isobaricInterpolation
!************************************************************************
subroutine vertInterpolation(v, p, rz, k)
  real(R_KIND), intent(out)  :: v
  real(R_KIND), intent(in)   :: p(:)
  real(R_KIND2), intent(in)  :: rz
  integer, intent(in)        :: k
!
  v = p(k) + rz * (p(k+1) - p(k))
!
end subroutine vertInterpolation
!************************************************************************


end module interpolate
