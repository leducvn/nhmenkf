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

module interpolate_TLAD
!  === Interpolation Package Module (Tangent Linear & Adjoint) ===
  use interpolate, only : makeGroundT, PosGrid, gsize, proj, &
                        & glat, glon, &
                        & FlagOverTop, FlagUnderGround, FlagRange, &
                        & rad

  implicit none
  save
#if TEST == 1
  public
#else
  private
  public :: TL_interpolation2D, AD_interpolation2D, &
          & TL_interpolation3D, AD_interpolation3D, &
          & TL_interpolation3D_T, AD_interpolation3D_T, &
          & TL_interpolation3D_Z, AD_interpolation3D_Z, &
          & TL_PolarInterpolation2D, AD_PolarInterpolation2D, &
          & TL_NPInterpolation2D_Wind, AD_NPInterpolation2D_Wind, &
          & TL_SPInterpolation2D_Wind, AD_SPInterpolation2D_Wind, &
!          & TL_NPInterpolation2D, AD_NPInterpolation2D, &
!          & TL_SPInterpolation2D, AD_SPInterpolation2D, &
          & TL_NPInterpolation3D_Wind, AD_NPInterpolation3D_Wind, &
          & TL_SPInterpolation3D_Wind, AD_SPInterpolation3D_Wind
!#-- 2004.12.11 NISHIJIMA start ---
!#-- need for SR8000
  public :: PosGrid
!#-- 2004.12.11 NISHIJIMA end   ---
#endif
!------------------------------------------------------------------------
!  Interface
  interface TL_interpolation2D
    module procedure TL_interpolation2D_Gauss, TL_interpolation2D_Standard
  end interface
  interface AD_interpolation2D
    module procedure AD_interpolation2D_Gauss, AD_interpolation2D_Standard
  end interface
  interface TL_interpolation3D
    module procedure TL_interpolation3D_Gauss, TL_interpolation3D_Standard
  end interface
  interface AD_interpolation3D
    module procedure AD_interpolation3D_Gauss, AD_interpolation3D_Standard
  end interface
  interface TL_interpolation3D_T
    module procedure TL_interpolation3D_T_Gauss, TL_interpolation3D_T_Standard
  end interface
  interface AD_interpolation3D_T
    module procedure AD_interpolation3D_T_Gauss, AD_interpolation3D_T_Standard
  end interface
  interface TL_interpolation3D_Z
    module procedure TL_interpolation3D_Z_Gauss, TL_interpolation3D_Z_Standard
  end interface
  interface AD_interpolation3D_Z
    module procedure AD_interpolation3D_Z_Gauss, AD_interpolation3D_Z_Standard
  end interface
!------------------------------------------------------------------------
!************************************************************************
contains
!************************************************************************
subroutine TL_bilinearInterpolation(Dv, Dv11, Dv21, Dv12, Dv22, x, y)
  real(R_KIND), intent(out) :: Dv
  real(R_KIND), intent(in)  :: Dv11, Dv21, Dv12, Dv22
  real(R_KIND2), intent(in) :: x, y
!
!#-- 2005.01.17 NISHIJIMA start ---
!  Dv = (1. - x) * (1. - y) * Dv11 &
!   & +       x  * (1. - y) * Dv21 &
!   & + (1. - x) *       y  * Dv12 &
!   & +       x  *       y  * Dv22
  Dv = (1.0d0 - x) * (1.0d0 - y) * Dv11 &
   & +          x  * (1.0d0 - y) * Dv21 &
   & + (1.0d0 - x) *          y  * Dv12 &
   & +          x  *          y  * Dv22
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine TL_bilinearInterpolation
!************************************************************************
subroutine AD_bilinearInterpolation(Dv, Dv11, Dv21, Dv12, Dv22, x, y)
  real(R_KIND), intent(in)    :: Dv
  real(R_KIND), intent(inout) :: Dv11, Dv21, Dv12, Dv22
  real(R_KIND2), intent(in)   :: x, y
!
!#-- 2005.01.17 NISHIJIMA start ---
!  Dv11 = (1. - x) * (1. - y) * Dv
!  Dv21 =       x  * (1. - y) * Dv
!  Dv12 = (1. - x) *       y  * Dv
!  Dv22 =       x  *       y  * Dv
  Dv11 = (1.0d0 - x) * (1.0d0 - y) * Dv
  Dv21 =          x  * (1.0d0 - y) * Dv
  Dv12 = (1.0d0 - x) *          y  * Dv
  Dv22 =          x  *          y  * Dv
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine AD_bilinearInterpolation
!************************************************************************
subroutine TL_vertInterpolation(Dv, Dp, Drz, Bp, Brz, k)
  real(R_KIND), intent(out)  :: Dv
  real(R_KIND), intent(in)   :: Dp(:)
  real(R_KIND2), intent(in)        :: Drz
  real(R_KIND), intent(in)   :: Bp(:)
  real(R_KIND2), intent(in)  :: Brz
  integer, intent(in)        :: k
!
!#-- 2005.01.17 NISHIJIMA start ---
!  Dv =       Brz         * Dp(k+1) &
!   & + (1. - Brz)        * Dp(k) &
!   & + (Bp(k+1) - Bp(k)) * Drz
  Dv =          Brz         * Dp(k+1) &
   & + (1.0d0 - Brz)        * Dp(k) &
   & + (Bp(k+1) - Bp(k))    * Drz
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine TL_vertInterpolation
!************************************************************************
subroutine AD_vertInterpolation(Dv, Dp, Drz, Bp, Brz, k)
  real(R_KIND), intent(in)     :: Dv
  real(R_KIND), intent(inout)  :: Dp(:)
  real(R_KIND2), intent(out)   :: Drz
  real(R_KIND), intent(in)     :: Bp(:)
  real(R_KIND2), intent(in)    :: Brz
  integer, intent(in)          :: k
!
!#-- 2005.01.17 NISHIJIMA start ---
!  Dp(k+1) =       Brz         * Dv
!  Dp(k)   = (1. - Brz)        * Dv
!  Drz     = (Bp(k+1) - Bp(k)) * Dv
  Dp(k+1) =          Brz         * Dv
  Dp(k)   = (1.0d0 - Brz)        * Dv
  Drz     = (Bp(k+1) - Bp(k))    * Dv
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine AD_vertInterpolation
!************************************************************************
subroutine TL_isobaricInterpolation &
         & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  real(R_KIND), intent(out)  :: Dvalue
  real(R_KIND), intent(in)   :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)  :: Dpg
  real(R_KIND), intent(in)   :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  type(PosGrid), intent(in)  :: Bpg
!
  real(R_KIND)               :: Dv11, Dv21, Dv12, Dv22
!
  call TL_vertInterpolation(Dv11, Dp11, Dpg%rz11, Bp11, Bpg%rz11, 1)
  call TL_vertInterpolation(Dv21, Dp21, Dpg%rz21, Bp21, Bpg%rz21, 1)
  call TL_vertInterpolation(Dv12, Dp12, Dpg%rz12, Bp12, Bpg%rz12, 1)
  call TL_vertInterpolation(Dv22, Dp22, Dpg%rz22, Bp22, Bpg%rz22, 1)
  call TL_bilinearInterpolation(Dvalue, Dv11, Dv21, Dv12, Dv22, &
                              & Bpg%rx, Bpg%ry)
!
end subroutine TL_isobaricInterpolation
!************************************************************************
subroutine AD_isobaricInterpolation &
         & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(inout) :: Dpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  type(PosGrid), intent(in)    :: Bpg
!
  real(R_KIND)               :: Dv11, Dv21, Dv12, Dv22
!
  call AD_bilinearInterpolation(Dvalue, Dv11, Dv21, Dv12, Dv22, &
                              & Bpg%rx, Bpg%ry)
  call AD_vertInterpolation(Dv11, Dp11, Dpg%rz11, Bp11, Bpg%rz11, 1)
  call AD_vertInterpolation(Dv21, Dp21, Dpg%rz21, Bp21, Bpg%rz21, 1)
  call AD_vertInterpolation(Dv12, Dp12, Dpg%rz12, Bp12, Bpg%rz12, 1)
  call AD_vertInterpolation(Dv22, Dp22, Dpg%rz22, Bp22, Bpg%rz22, 1)
!
end subroutine AD_isobaricInterpolation
!************************************************************************
subroutine TL_makeGroundT(Dt, Dt1, Ddlev, Bt1, Bdlev, lElement)
  real(R_KIND), intent(out) :: Dt     ! temperature under the ground
  real(R_KIND), intent(in)  :: Dt1    ! temperature at the surface
  real(R_KIND2), intent(in) :: Ddlev  ! difference of levels
  real(R_KIND), intent(in)  :: Bt1    ! temperature at the surface
  real(R_KIND2), intent(in) :: Bdlev  ! difference of levels
  character(len=6), intent(in)   :: lElement  ! element name of level
!
  real(R_KIND), parameter :: &
  & gamma = 0.005d0, &  ! temperature lapse rate (K/-m)
  & rd = 287.05d0, &    ! dry gas constant
  & grv = 9.80665d0     ! gravitational acceleration
!
  if (lElement == "logP  ") then      ! dlev > 1
!   Nomary t1 must be virtual temperature, but supposed as q = 0
    Dt = exp(Bdlev*gamma*rd/grv) * Dt1 &
    &  + gamma*rd/grv * exp(Bdlev*gamma*rd/grv) * Bt1 * Ddlev
  else if (lElement == "Z     ") then ! dlev < 0
    Dt = Dt1 &
    &  - gamma * Ddlev
  else
    Dt = Dt1
  end if
end subroutine TL_makeGroundT
!************************************************************************
subroutine AD_makeGroundT(Dt, Dt1, Ddlev, Bt1, Bdlev, lElement)
  real(R_KIND), intent(in)    :: Dt     ! temperature under the ground
  real(R_KIND), intent(inout) :: Dt1    ! temperature at the surface
  real(R_KIND2), intent(inout):: Ddlev  ! difference of levels
  real(R_KIND), intent(in)    :: Bt1    ! temperature at the surface
  real(R_KIND2), intent(in)   :: Bdlev  ! difference of levels
  character(len=6), intent(in)   :: lElement  ! element name of level
!
  real(R_KIND), parameter :: &
  & gamma = 0.005d0, &  ! temperature lapse rate (K/-m)
  & rd = 287.05d0, &    ! dry gas constant
  & grv = 9.80665d0     ! gravitational acceleration
!
  if (lElement == "logP  ") then      ! dlev > 1
!   Nomary t1 must be virtual temperature, but supposed as q = 0
    Dt1   = exp(Bdlev*gamma*rd/grv) * Dt
    Ddlev = gamma*rd/grv * exp(Bdlev*gamma*rd/grv) * Bt1 * Dt
  else if (lElement == "Z     ") then ! dlev < 0
    Dt1   = Dt
    Ddlev = -gamma * Dt
  else
    Dt1 = Dt
  end if
end subroutine AD_makeGroundT
!************************************************************************
subroutine TL_getZrate(Dz, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(out):: Dz
  real(R_KIND), intent(in)  :: Dlevs(:)
  real(R_KIND), intent(in)  :: Blevs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
!
!#-- 2005.01.17 NISHIJIMA start ---
!#  if (flag == FlagOverTop) then
!#!  --- over the top of model level ---
!#    Dz = 0.
!#  else if (flag == FlagUnderGround) then
!#!  --- under the ground ---
!#    Dz = 0.
!#  else
!#    Bdlev2inv = 1 / ((Blevs(k+1) - Blevs(k)) ** 2)
!#  end if
  if (flag == FlagOverTop) then
!  --- over the top of model level ---
    Dz = 0.0d0
  else if (flag == FlagUnderGround) then
!  --- under the ground ---
    Dz = 0.0d0
  else
    Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
    Dz = (Blevs(k) - lev)   * Bdlev2inv * Dlevs(k+1) &
     & + (lev - Blevs(k+1)) * Bdlev2inv * Dlevs(k)
  end if
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine TL_getZrate
!************************************************************************
subroutine AD_getZrate(Dz, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(in)   :: Dz
  real(R_KIND), intent(inout) :: Dlevs(:)
  real(R_KIND), intent(in)    :: Blevs(:)
  real(R_KIND2), intent(in)   :: lev
  integer(4), intent(in)      :: k
  integer(4), intent(in)      :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
!
!#-- 2005.01.17 NISHIJIMA start ---
!#  if (flag == FlagOverTop) then
!#!  --- over the top of model level ---
!#    Dlevs(k:k+1) = 0.
!#  else if (flag == FlagUnderGround) then
!#!  --- under the ground ---
!#    Dlevs(k:k+1) = 0.
!#  else
!#    Bdlev2inv = 1 / ((Blevs(k+1) - Blevs(k)) ** 2)
!#  end if
!#
  if (flag == FlagOverTop) then
!  --- over the top of model level ---
    Dlevs(k:k+1) = 0.0d0
  else if (flag == FlagUnderGround) then
!  --- under the ground ---
    Dlevs(k:k+1) = 0.0d0
  else
    Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
    Dlevs(k+1) = (Blevs(k) - lev)   * Bdlev2inv * Dz
    Dlevs(k)   = (lev - Blevs(k+1)) * Bdlev2inv * Dz
  end if
!#-- 2005.01.17 NISHIJIMA end   ---
end subroutine AD_getZrate
!************************************************************************
subroutine TL_getZrate_T(Dz, Ddlev, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(out):: Dz
  real(R_KIND2), intent(out):: Ddlev ! height from surfce
  real(R_KIND), intent(in)  :: Dlevs(:)
  real(R_KIND), intent(in)  :: Blevs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
!
  if (flag == FlagUnderGround) then
!#-- 2005.01.17 NISHIJIMA start ---
!#    Dz = 0.
    Dz = 0.0d0
!#-- 2005.01.17 NISHIJIMA end   ---
    Ddlev = - Dlevs(1)
  else ! including "if (flag == FlagOverTop)", exterpolation
    Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
    Dz = (Blevs(k) - lev)   * Bdlev2inv * Dlevs(k+1) &
     & + (lev - Blevs(k+1)) * Bdlev2inv * Dlevs(k)
  end if
end subroutine TL_getZrate_T
!************************************************************************
subroutine AD_getZrate_T(Dz, Ddlev, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(in)   :: Dz
  real(R_KIND2), intent(in)   :: Ddlev ! height from surfce
  real(R_KIND), intent(inout) :: Dlevs(:)
  real(R_KIND), intent(in)    :: Blevs(:)
  real(R_KIND2), intent(in)   :: lev
  integer(4), intent(in)      :: k
  integer(4), intent(in)      :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
!
  if (flag == FlagUnderGround) then
    Dlevs(1:2) = 0.
    Dlevs(1) = - Ddlev
  else ! including "if (flag == FlagOverTop)", exterpolation
    Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
    Dlevs(k+1) = (Blevs(k) - lev)   * Bdlev2inv * Dz
    Dlevs(k)   = (lev - Blevs(k+1)) * Bdlev2inv * Dz
  end if
end subroutine AD_getZrate_T
!************************************************************************
subroutine TL_getZrate_Z(Dz, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(out):: Dz
  real(R_KIND), intent(in)  :: Dlevs(:)
  real(R_KIND), intent(in)  :: Blevs(:)
  real(R_KIND2), intent(in) :: lev
  integer(4), intent(in)    :: k
  integer(4), intent(in)    :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
! including "if (flag == FlagOverTop .and. FlagUnderGround)", exterpolation
  Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
  Dz = (Blevs(k) - lev)   * Bdlev2inv * Dlevs(k+1) &
   & + (lev - Blevs(k+1)) * Bdlev2inv * Dlevs(k)
end subroutine TL_getZrate_Z
!************************************************************************
subroutine AD_getZrate_Z(Dz, Dlevs, Blevs, lev, k, flag)
  real(R_KIND2), intent(in)   :: Dz
  real(R_KIND), intent(inout) :: Dlevs(:)
  real(R_KIND), intent(in)    :: Blevs(:)
  real(R_KIND2), intent(in)   :: lev
  integer(4), intent(in)      :: k
  integer(4), intent(in)      :: flag ! flag of under the ground or over top
!
  real(R_KIND)              :: Bdlev2inv
! including "if (flag == FlagOverTop .and. FlagUnderGround)", exterpolation
  Bdlev2inv = 1.0d0 / ((Blevs(k+1) - Blevs(k)) ** 2)
  Dlevs(k+1) = (Blevs(k) - lev)   * Bdlev2inv * Dz
  Dlevs(k)   = (lev - Blevs(k+1)) * Bdlev2inv * Dz
end subroutine AD_getZrate_Z
!************************************************************************
subroutine TL_convert_Lev_to_GridPos &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(out)   :: Dpg
  real(R_KIND), intent(in)     :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call TL_getZrate(Dpg%rz11, Dp11, Bp11, lev, 1, Bpg%range11)
  call TL_getZrate(Dpg%rz21, Dp21, Bp21, lev, 1, Bpg%range21)
  call TL_getZrate(Dpg%rz12, Dp12, Bp12, lev, 1, Bpg%range12)
  call TL_getZrate(Dpg%rz22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine TL_convert_Lev_to_GridPos
!************************************************************************
subroutine AD_convert_Lev_to_GridPos &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(in)    :: Dpg
  real(R_KIND), intent(inout)  :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call AD_getZrate(Dpg%rz11, Dp11, Bp11, lev, 1, Bpg%range11)
  call AD_getZrate(Dpg%rz21, Dp21, Bp21, lev, 1, Bpg%range21)
  call AD_getZrate(Dpg%rz12, Dp12, Bp12, lev, 1, Bpg%range12)
  call AD_getZrate(Dpg%rz22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine AD_convert_Lev_to_GridPos
!************************************************************************
subroutine TL_convert_Lev_to_GridPos_T &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(out)   :: Dpg
  real(R_KIND), intent(in)     :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call TL_getZrate_T(Dpg%rz11, Dpg%dlev11, Dp11, Bp11, lev, 1, Bpg%range11)
  call TL_getZrate_T(Dpg%rz21, Dpg%dlev21, Dp21, Bp21, lev, 1, Bpg%range21)
  call TL_getZrate_T(Dpg%rz12, Dpg%dlev12, Dp12, Bp12, lev, 1, Bpg%range12)
  call TL_getZrate_T(Dpg%rz22, Dpg%dlev22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine TL_convert_Lev_to_GridPos_T
!************************************************************************
subroutine AD_convert_Lev_to_GridPos_T &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(in)    :: Dpg
  real(R_KIND), intent(inout)  :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call AD_getZrate_T(Dpg%rz11, Dpg%dlev11, Dp11, Bp11, lev, 1, Bpg%range11)
  call AD_getZrate_T(Dpg%rz21, Dpg%dlev21, Dp21, Bp21, lev, 1, Bpg%range21)
  call AD_getZrate_T(Dpg%rz12, Dpg%dlev12, Dp12, Bp12, lev, 1, Bpg%range12)
  call AD_getZrate_T(Dpg%rz22, Dpg%dlev22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine AD_convert_Lev_to_GridPos_T
!************************************************************************
subroutine TL_convert_Lev_to_GridPos_Z &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(out)   :: Dpg
  real(R_KIND), intent(in)     :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call TL_getZrate_Z(Dpg%rz11, Dp11, Bp11, lev, 1, Bpg%range11)
  call TL_getZrate_Z(Dpg%rz21, Dp21, Bp21, lev, 1, Bpg%range21)
  call TL_getZrate_Z(Dpg%rz12, Dp12, Bp12, lev, 1, Bpg%range12)
  call TL_getZrate_Z(Dpg%rz22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine TL_convert_Lev_to_GridPos_Z
!************************************************************************
subroutine AD_convert_Lev_to_GridPos_Z &
         & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
  type(PosGrid), intent(in)    :: Dpg
  real(R_KIND), intent(inout)  :: Dp11(:), Dp21(:), Dp12(:), Dp22(:)
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bp11(:), Bp21(:), Bp12(:), Bp22(:)
  real(R_KIND2), intent(in)    :: lev
!
  call AD_getZrate_Z(Dpg%rz11, Dp11, Bp11, lev, 1, Bpg%range11)
  call AD_getZrate_Z(Dpg%rz21, Dp21, Bp21, lev, 1, Bpg%range21)
  call AD_getZrate_Z(Dpg%rz12, Dp12, Bp12, lev, 1, Bpg%range12)
  call AD_getZrate_Z(Dpg%rz22, Dp22, Bp22, lev, 1, Bpg%range22)
end subroutine AD_convert_Lev_to_GridPos_Z
!************************************************************************
subroutine TL_interpolation3D_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(in)     :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(in)     :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Dp11(1:2) = DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_Gauss
!************************************************************************
subroutine AD_interpolation3D_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(inout)  :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(inout)  :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)  
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_Gauss
!************************************************************************
subroutine TL_interpolation3D_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_Standard
!************************************************************************
subroutine AD_interpolation3D_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_Standard
!************************************************************************
subroutine TL_interpolation3D_T_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv, lElement)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(in)     :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(in)     :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
  character(len=6), intent(in) :: lElement  ! element name of level
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Dp11(1:2) = DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    if (Bpg%range11 /= FlagUnderGround) then
      Dp11(1:2) = DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
      Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    else
      call TL_makeGroundT(Dp11(1), DNPgpv(Bpg%i1,1), Dpg%dlev11, &
                        & BNPgpv(Bpg%i1,1), Bpg%dlev11, lElement)
      call makeGroundT(Bp11(1), BNPgpv(Bpg%i1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dp21(1:2) = DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
      Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    else
      call TL_makeGroundT(Dp21(1), DNPgpv(Bpg%i2,1), Dpg%dlev21, &
                        & BNPgpv(Bpg%i2,1), Bpg%dlev21, lElement)
      call makeGroundT(Bp21(1), BNPgpv(Bpg%i2,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else
      call TL_makeGroundT(Dp12(1), Dgpv(Bpg%i1,Bpg%j2,1), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call TL_makeGroundT(Dp22(1), Dgpv(Bpg%i2,Bpg%j2,1), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    if (Bpg%range11 /= FlagUnderGround) then
      Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call TL_makeGroundT(Dp11(1), Dgpv(Bpg%i1,Bpg%j1,1), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call TL_makeGroundT(Dp21(1), Dgpv(Bpg%i2,Bpg%j1,1), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dp12(1:2) = DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
      Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    else
      call TL_makeGroundT(Dp12(1), DSPgpv(Bpg%i1,1), Dpg%dlev12, &
                        & BSPgpv(Bpg%i1,1), Bpg%dlev12, lElement)
      call makeGroundT(Bp12(1), BSPgpv(Bpg%i1,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dp22(1:2) = DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
      Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    else
      call TL_makeGroundT(Dp22(1), DSPgpv(Bpg%i2,1), Dpg%dlev22, &
                        & BSPgpv(Bpg%i2,1), Bpg%dlev22, lElement)
      call makeGroundT(Bp22(1), BSPgpv(Bpg%i2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    if (Bpg%range11 /= FlagUnderGround) then
      Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call TL_makeGroundT(Dp11(1), Dgpv(Bpg%i1,Bpg%j1,1), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call TL_makeGroundT(Dp21(1), Dgpv(Bpg%i2,Bpg%j1,1), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else                             
      call TL_makeGroundT(Dp12(1), Dgpv(Bpg%i1,Bpg%j2,1), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call TL_makeGroundT(Dp22(1), Dgpv(Bpg%i2,Bpg%j2,1), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_T_Gauss
!************************************************************************
subroutine AD_interpolation3D_T_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv, lElement)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(inout)  :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(inout)  :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
  character(len=6), intent(in) :: lElement  ! element name of level
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    if (Bpg%range11 /= FlagUnderGround) then
      Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    else
      call makeGroundT(Bp11(1), BNPgpv(Bpg%i1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    else
      call makeGroundT(Bp21(1), BNPgpv(Bpg%i2,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    if (Bpg%range11 /= FlagUnderGround) then
      DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    else
      call AD_makeGroundT(Dp11(1), Dp11(2), Dpg%dlev11, &
                        & BNPgpv(Bpg%i1,1), Bpg%dlev11, lElement)
      DNPgpv(Bpg%i1,1) = Dp11(2) + DNPgpv(Bpg%i1,1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    else
      call AD_makeGroundT(Dp21(1), Dp21(2), Dpg%dlev21, &
                        & BNPgpv(Bpg%i2,1), Bpg%dlev21, lElement)
      DNPgpv(Bpg%i2,1) = Dp21(2) + DNPgpv(Bpg%i2,1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)  
    else
      call AD_makeGroundT(Dp12(1), Dp12(2), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Dgpv(Bpg%i1,Bpg%j2,1) = Dp12(2) + Dgpv(Bpg%i1,Bpg%j2,1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call AD_makeGroundT(Dp22(1), Dp22(2), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Dgpv(Bpg%i2,Bpg%j2,1) = Dp22(2) + Dgpv(Bpg%i2,Bpg%j2,1)
    end if
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    if (Bpg%range11 /= FlagUnderGround) then
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    else
      call makeGroundT(Bp12(1), BSPgpv(Bpg%i1,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    else
      call makeGroundT(Bp22(1), BSPgpv(Bpg%i2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    if (Bpg%range11 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call AD_makeGroundT(Dp11(1), Dp11(2), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Dgpv(Bpg%i1,Bpg%j1,1) = Dp11(2) + Dgpv(Bpg%i1,Bpg%j1,1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call AD_makeGroundT(Dp21(1), Dp21(2), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Dgpv(Bpg%i2,Bpg%j1,1) = Dp21(2) + Dgpv(Bpg%i2,Bpg%j1,1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
     DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    else
      call AD_makeGroundT(Dp12(1), Dp12(2), Dpg%dlev12, &
                        & BSPgpv(Bpg%i1,1), Bpg%dlev12, lElement)
      DSPgpv(Bpg%i1,1) = Dp12(2) + DSPgpv(Bpg%i1,1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    else
      call AD_makeGroundT(Dp22(1), Dp22(2), Dpg%dlev22, &
                        & BSPgpv(Bpg%i2,1), Bpg%dlev22, lElement)
      DSPgpv(Bpg%i2,1) = Dp22(2) + DSPgpv(Bpg%i2,1)
    end if
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    if (Bpg%range11 /= FlagUnderGround) then
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else                             
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    if (Bpg%range11 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call AD_makeGroundT(Dp11(1), Dp11(2), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Dgpv(Bpg%i1,Bpg%j1,1) = Dp11(2) + Dgpv(Bpg%i1,Bpg%j1,1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call AD_makeGroundT(Dp21(1), Dp21(2), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Dgpv(Bpg%i2,Bpg%j1,1) = Dp21(2) + Dgpv(Bpg%i2,Bpg%j1,1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)  
    else
      call AD_makeGroundT(Dp12(1), Dp12(2), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Dgpv(Bpg%i1,Bpg%j2,1) = Dp12(2) + Dgpv(Bpg%i1,Bpg%j2,1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call AD_makeGroundT(Dp22(1), Dp22(2), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Dgpv(Bpg%i2,Bpg%j2,1) = Dp22(2) + Dgpv(Bpg%i2,Bpg%j2,1)
    end if
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_T_Gauss
!************************************************************************
subroutine TL_interpolation3D_T_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv, lElement)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  character(len=6), intent(in) :: lElement  ! element name of level
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    if (Bpg%range11 /= FlagUnderGround) then
      Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call TL_makeGroundT(Dp11(1), Dgpv(Bpg%i1,Bpg%j1,1), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call TL_makeGroundT(Dp21(1), Dgpv(Bpg%i2,Bpg%j1,1), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else                             
      call TL_makeGroundT(Dp12(1), Dgpv(Bpg%i1,Bpg%j2,1), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call TL_makeGroundT(Dp22(1), Dgpv(Bpg%i2,Bpg%j2,1), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_T_Standard
!************************************************************************
subroutine AD_interpolation3D_T_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv, lElement)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  character(len=6), intent(in) :: lElement  ! element name of level
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    if (Bpg%range11 /= FlagUnderGround) then
      Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call makeGroundT(Bp11(1), Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Bp11(2) = Bp11(1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call makeGroundT(Bp21(1), Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Bp21(2) = Bp21(1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    else                             
      call makeGroundT(Bp12(1), Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Bp12(2) = Bp12(1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call makeGroundT(Bp22(1), Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Bp22(2) = Bp22(1)
    end if
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    if (Bpg%range11 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    else
      call AD_makeGroundT(Dp11(1), Dp11(2), Dpg%dlev11, &
                        & Bgpv(Bpg%i1,Bpg%j1,1), Bpg%dlev11, lElement)
      Dgpv(Bpg%i1,Bpg%j1,1) = Dp11(2) + Dgpv(Bpg%i1,Bpg%j1,1)
    end if
    if (Bpg%range21 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    else
      call AD_makeGroundT(Dp21(1), Dp21(2), Dpg%dlev21, &
                        & Bgpv(Bpg%i2,Bpg%j1,1), Bpg%dlev21, lElement)
      Dgpv(Bpg%i2,Bpg%j1,1) = Dp21(2) + Dgpv(Bpg%i2,Bpg%j1,1)
    end if
    if (Bpg%range12 /= FlagUnderGround) then
      Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)  
    else
      call AD_makeGroundT(Dp12(1), Dp12(2), Dpg%dlev12, &
                        & Bgpv(Bpg%i1,Bpg%j2,1), Bpg%dlev12, lElement)
      Dgpv(Bpg%i1,Bpg%j2,1) = Dp12(2) + Dgpv(Bpg%i1,Bpg%j2,1)
    end if
    if (Bpg%range22 /= FlagUnderGround) then
      Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    else
      call AD_makeGroundT(Dp22(1), Dp22(2), Dpg%dlev22, &
                        & Bgpv(Bpg%i2,Bpg%j2,1), Bpg%dlev22, lElement)
      Dgpv(Bpg%i2,Bpg%j2,1) = Dp22(2) + Dgpv(Bpg%i2,Bpg%j2,1)
    end if
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_T &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_T_Standard
!************************************************************************
subroutine TL_interpolation3D_Z_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(in)     :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(in)     :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Dp11(1:2) = DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_Z_Gauss
!************************************************************************
subroutine AD_interpolation3D_Z_Gauss &
         & (Dvalue, Dgpv, Dlgpv, &
          & DNPgpv, DNPlgpv, DSPgpv, DSPlgpv, &
          & Bpg, lev, Bgpv, Blgpv, &
          & BNPgpv, BNPlgpv, BSPgpv, BSPlgpv)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND), intent(inout)  :: DNPgpv(:,:), DNPlgpv(:,:)
  real(R_KIND), intent(inout)  :: DSPgpv(:,:), DSPlgpv(:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
  real(R_KIND), intent(in)     :: BNPgpv(:,:), BNPlgpv(:,:)
  real(R_KIND), intent(in)     :: BSPgpv(:,:), BSPlgpv(:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Bp11(1:2) = BNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)  
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = BNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = BNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + DNPlgpv(Bpg%i1,Bpg%k11:Bpg%k11+1)
    DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + DNPlgpv(Bpg%i2,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = BSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = BSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + DSPlgpv(Bpg%i1,Bpg%k12:Bpg%k12+1)
    DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + DSPlgpv(Bpg%i2,Bpg%k22:Bpg%k22+1)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_Z_Gauss
!************************************************************************
subroutine TL_interpolation3D_Z_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv)
  real(R_KIND), intent(out)    :: Dvalue
  real(R_KIND), intent(in)     :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    Dp11(1:2) = Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dp11(1:2) = Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dp21(1:2) = Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dp12(1:2) = Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dp22(1:2) = Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call TL_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
  end if
end subroutine TL_interpolation3D_Z_Standard
!************************************************************************
subroutine AD_interpolation3D_Z_Standard &
         & (Dvalue, Dgpv, Dlgpv, &
          & Bpg, lev, Bgpv, Blgpv)
  real(R_KIND), intent(in)     :: Dvalue
  real(R_KIND), intent(inout)  :: Dgpv(:,:,:), Dlgpv(:,:,:)
  real(R_KIND2), intent(in)    :: lev
  type(PosGrid), intent(in)    :: Bpg
  real(R_KIND), intent(in)     :: Bgpv(:,:,:), Blgpv(:,:,:)
!
  type(PosGrid)             :: Dpg   ! Dpg%rz?? are only used
!  --- the pillar of the grid in the surroundings ---
  real(R_KIND)              :: Dp11(2), Dp21(2), Dp12(2), Dp22(2)
  real(R_KIND)              :: Bp11(2), Bp21(2), Bp12(2), Bp22(2)
!
  if (Bpg%j1 > 0) then
    Bp11(1:2) = Bgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Bgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Bgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Bgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_isobaricInterpolation &
    & (Dvalue, Dp11, Dp21, Dp12, Dp22, Dpg, Bp11, Bp21, Bp12, Bp22, Bpg)
    Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    Bp11(1:2) = Blgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Bp21(1:2) = Blgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Bp12(1:2) = Blgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Bp22(1:2) = Blgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
    call AD_convert_Lev_to_GridPos_Z &
    & (Dpg, Dp11, Dp21, Dp12, Dp22, Bpg, Bp11, Bp21, Bp12, Bp22, lev)
    Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1) = Dp11(1:2) + Dlgpv(Bpg%i1,Bpg%j1,Bpg%k11:Bpg%k11+1)
    Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1) = Dp21(1:2) + Dlgpv(Bpg%i2,Bpg%j1,Bpg%k21:Bpg%k21+1)
    Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1) = Dp12(1:2) + Dlgpv(Bpg%i1,Bpg%j2,Bpg%k12:Bpg%k12+1)
    Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1) = Dp22(1:2) + Dlgpv(Bpg%i2,Bpg%j2,Bpg%k22:Bpg%k22+1)
  end if
end subroutine AD_interpolation3D_Z_Standard
!************************************************************************
subroutine TL_interpolation2D_Gauss(Dvalue, Dgpv, DNPgpv, DSPgpv, Bpg)
  real(R_KIND), intent(out) :: Dvalue
  real(R_KIND), intent(in)  :: Dgpv(:,:)
  real(R_KIND), intent(in)  :: DNPgpv(:), DSPgpv(:)
  type(PosGrid), intent(in) :: Bpg
!
  real(R_KIND)              :: Dv(4)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    Dv(1) = DNPgpv(Bpg%i1)
    Dv(2) = DNPgpv(Bpg%i2)
    Dv(3) = Dgpv(Bpg%i1,Bpg%j2)
    Dv(4) = Dgpv(Bpg%i2,Bpg%j2)
    call TL_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    Dv(1) = Dgpv(Bpg%i1,Bpg%j1)
    Dv(2) = Dgpv(Bpg%i2,Bpg%j1)
    Dv(3) = DSPgpv(Bpg%i1)
    Dv(4) = DSPgpv(Bpg%i2)
    call TL_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    Dv(1) = Dgpv(Bpg%i1,Bpg%j1)
    Dv(2) = Dgpv(Bpg%i2,Bpg%j1)
    Dv(3) = Dgpv(Bpg%i1,Bpg%j2)
    Dv(4) = Dgpv(Bpg%i2,Bpg%j2)
    call TL_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
  end if
end subroutine TL_interpolation2D_Gauss
!************************************************************************
subroutine AD_interpolation2D_Gauss(Dvalue, Dgpv, DNPgpv, DSPgpv, Bpg)
  real(R_KIND), intent(in)    :: Dvalue
  real(R_KIND), intent(inout) :: Dgpv(:,:)
  real(R_KIND), intent(inout) :: DNPgpv(:), DSPgpv(:)
  type(PosGrid), intent(in)   :: Bpg
!
  real(R_KIND)              :: Dv(4)
!  --- near the North Pole ---
  if (Bpg%j1 == 0 .and. proj == "GS  ") then
    call AD_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
    DNPgpv(Bpg%i1) = Dv(1) + DNPgpv(Bpg%i1)
    DNPgpv(Bpg%i2) = Dv(2) + DNPgpv(Bpg%i2)
    Dgpv(Bpg%i1,Bpg%j2) = Dv(3) + Dgpv(Bpg%i1,Bpg%j2)
    Dgpv(Bpg%i2,Bpg%j2) = Dv(4) + Dgpv(Bpg%i2,Bpg%j2)
!  --- near the South Pole ---
  else if (Bpg%j1 == gsize(2) .and. proj == "GS  ") then
    call AD_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
    Dgpv(Bpg%i1,Bpg%j1) = Dv(1) + Dgpv(Bpg%i1,Bpg%j1)
    Dgpv(Bpg%i2,Bpg%j1) = Dv(2) + Dgpv(Bpg%i2,Bpg%j1)
    DSPgpv(Bpg%i1) = Dv(3) + DSPgpv(Bpg%i1)
    DSPgpv(Bpg%i2) = Dv(4) + DSPgpv(Bpg%i2)
!  --- other area ---
  else if (Bpg%j1 > 0) then
    call AD_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
    Dgpv(Bpg%i1,Bpg%j1) = Dv(1) + Dgpv(Bpg%i1,Bpg%j1)
    Dgpv(Bpg%i2,Bpg%j1) = Dv(2) + Dgpv(Bpg%i2,Bpg%j1)
    Dgpv(Bpg%i1,Bpg%j2) = Dv(3) + Dgpv(Bpg%i1,Bpg%j2)
    Dgpv(Bpg%i2,Bpg%j2) = Dv(4) + Dgpv(Bpg%i2,Bpg%j2)
  end if
end subroutine AD_interpolation2D_Gauss
!************************************************************************
subroutine TL_interpolation2D_Standard(Dvalue, Dgpv, Bpg)
  real(R_KIND), intent(out) :: Dvalue
  real(R_KIND), intent(in)  :: Dgpv(:,:)
  type(PosGrid), intent(in) :: Bpg
!
  real(R_KIND)              :: Dv(4)
!
  if (Bpg%j1 > 0) then
    Dv(1) = Dgpv(Bpg%i1,Bpg%j1)
    Dv(2) = Dgpv(Bpg%i2,Bpg%j1)
    Dv(3) = Dgpv(Bpg%i1,Bpg%j2)
    Dv(4) = Dgpv(Bpg%i2,Bpg%j2)
    call TL_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
  end if
end subroutine TL_interpolation2D_Standard
!************************************************************************
subroutine AD_interpolation2D_Standard(Dvalue, Dgpv, Bpg)
  real(R_KIND), intent(in)    :: Dvalue
  real(R_KIND), intent(inout) :: Dgpv(:,:)
  type(PosGrid), intent(in)   :: Bpg
!
  real(R_KIND)              :: Dv(4)
!
  if (Bpg%j1 > 0) then
    call AD_bilinearInterpolation(Dvalue, &
    & Dv(1), Dv(2), Dv(3), Dv(4), &
    & Bpg%rx, Bpg%ry)
    Dgpv(Bpg%i1,Bpg%j1) = Dv(1) + Dgpv(Bpg%i1,Bpg%j1)
    Dgpv(Bpg%i2,Bpg%j1) = Dv(2) + Dgpv(Bpg%i2,Bpg%j1)
    Dgpv(Bpg%i1,Bpg%j2) = Dv(3) + Dgpv(Bpg%i1,Bpg%j2)
    Dgpv(Bpg%i2,Bpg%j2) = Dv(4) + Dgpv(Bpg%i2,Bpg%j2)
  end if
end subroutine AD_interpolation2D_Standard
!************************************************************************
subroutine TL_PolarInterpolation2D(DPValue, Dgpv)
!  === get interpolated vaule at the North/South Pole ===
!      for GAUSSIAN GRID only
  real(R_KIND), intent(out)    :: DPValue(:)
  real(R_KIND), intent(in)     :: Dgpv(:)
!
  real(R_KIND)                 :: DPV
!
  if (proj == "GS  ") then
!   avarage of gpv near the pole
    DPV = sum(Dgpv(1:gsize(1))) / gsize(1)
    DPValue(1:gsize(1)) = DPV
  end if
end subroutine TL_PolarInterpolation2D
!************************************************************************
subroutine AD_PolarInterpolation2D(DPValue, Dgpv)
!  === get interpolated vaule at the North/South Pole ===
!      for GAUSSIAN GRID only
  real(R_KIND), intent(in)     :: DPValue(:)
  real(R_KIND), intent(inout)  :: Dgpv(:)
!
  real(R_KIND)                 :: DPV
!
  if (proj == "GS  ") then
!   avarage of gpv near the pole
    DPV = sum(DPValue(1:gsize(1))) / gsize(1)
    Dgpv(1:gsize(1)) = DPV + Dgpv(1:gsize(1))
  end if
end subroutine AD_PolarInterpolation2D
!************************************************************************
subroutine TL_NPInterpolation2D_Wind(DNValueU, DNValueV, DgpvU, DgpvV)
!  === get interpolated vaule (Wind) at the North Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(out)   :: DNValueU(:), DNValueV(:)
  real(R_KIND), intent(in)    :: DgpvU(:), DgpvV(:)
!
  real(R_KIND) :: NPU, NPV
!
  if (proj == "GS  ") then
    NPU = sum( DgpvU(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) &
           & - DgpvV(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
    NPV = sum( DgpvU(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) &
           & + DgpvV(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
    DNValueU(1:gsize(1)) =  NPU * cos(glon(1:gsize(1),1) * rad) &
                       & + NPV * sin(glon(1:gsize(1),1) * rad)
    DNValueV(1:gsize(1)) = -NPU * sin(glon(1:gsize(1),1) * rad) &
                       & + NPV * cos(glon(1:gsize(1),1) * rad)
  end if
end subroutine TL_NPInterpolation2D_Wind
!************************************************************************
subroutine AD_NPInterpolation2D_Wind(DNValueU, DNValueV, DgpvU, DgpvV)
!  === get interpolated vaule (Wind) at the North Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(in)    :: DNValueU(:), DNValueV(:)
  real(R_KIND), intent(inout) :: DgpvU(:), DgpvV(:)
!
  real(R_KIND) :: NPU, NPV
!
  if (proj == "GS  ") then
    NPU = sum( DNValueU(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) &
           & - DNValueV(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
    NPV = sum( DNValueU(1:gsize(1)) * sin(glon(1:gsize(1),1) * rad) &
           & + DNValueV(1:gsize(1)) * cos(glon(1:gsize(1),1) * rad) ) &
        & / gsize(1)
    DgpvU(1:gsize(1)) =  NPU * cos(glon(1:gsize(1),1) * rad) &
                    & + NPV * sin(glon(1:gsize(1),1) * rad) &
                    & + DgpvU(1:gsize(1))
    DgpvV(1:gsize(1)) = -NPU * sin(glon(1:gsize(1),1) * rad) &
                    & + NPV * cos(glon(1:gsize(1),1) * rad) &
                    & + DgpvV(1:gsize(1))
  end if
end subroutine AD_NPInterpolation2D_Wind
!************************************************************************
subroutine TL_SPInterpolation2D_Wind(DSValueU, DSValueV, DgpvU, DgpvV)
!  === get interpolated vaule (Wind) at the South Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(out)   :: DSValueU(:), DSValueV(:)
  real(R_KIND), intent(in)    :: DgpvU(:), DgpvV(:)
!
  real(R_KIND) :: SPU, SPV
!
  if (proj == "GS  ") then
    SPU = sum( DgpvU(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) &
           & + DgpvV(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
    SPV = sum(-DgpvU(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) &
           & + DgpvV(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
    DSValueU(1:gsize(1)) =  SPU * cos(glon(1:gsize(1),gsize(2)) * rad) &
                       & - SPV * sin(glon(1:gsize(1),gsize(2)) * rad)
    DSValueV(1:gsize(1)) =  SPU * sin(glon(1:gsize(1),gsize(2)) * rad) &
                       & + SPV * cos(glon(1:gsize(1),gsize(2)) * rad)
  end if
end subroutine TL_SPInterpolation2D_Wind
!************************************************************************
subroutine AD_SPInterpolation2D_Wind(DSValueU, DSValueV, DgpvU, DgpvV)
!  === get interpolated vaule (Wind) at the South Pole ===
!      for GAUSSIAN GRID only
!      original Qc/Src/Cda2/poleges.f
  real(R_KIND), intent(in)    :: DSValueU(:), DSValueV(:)
  real(R_KIND), intent(inout) :: DgpvU(:), DgpvV(:)
!
  real(R_KIND) :: SPU, SPV
!DgpvU DgpvV
  if (proj == "GS  ") then
    SPU = sum( DSValueU(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) &
           & + DSValueV(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
    SPV = sum(-DSValueU(1:gsize(1)) * sin(glon(1:gsize(1),gsize(2)) * rad) &
           & + DSValueV(1:gsize(1)) * cos(glon(1:gsize(1),gsize(2)) * rad) ) &
        & / gsize(1)
    DgpvU(1:gsize(1)) =  SPU * cos(glon(1:gsize(1),gsize(2)) * rad) &
                    & - SPV * sin(glon(1:gsize(1),gsize(2)) * rad) &
                    & + DgpvU(1:gsize(1))
    DgpvV(1:gsize(1)) =  SPU * sin(glon(1:gsize(1),gsize(2)) * rad) &
                    & + SPV * cos(glon(1:gsize(1),gsize(2)) * rad) &
                    & + DgpvV(1:gsize(1))
  end if
end subroutine AD_SPInterpolation2D_Wind
!************************************************************************
subroutine TL_NPInterpolation3D(DNValue, Dgpv)
  real(R_KIND), intent(out)   :: DNValue(:,:)
  real(R_KIND), intent(in)    :: Dgpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call TL_PolarInterpolation2D(DNValue(1:gsize(1),k), &
                                 & Dgpv(1:gsize(1),1,k))
    end do
  end if
end subroutine TL_NPInterpolation3D
!************************************************************************
subroutine AD_NPInterpolation3D(DNValue, Dgpv)
  real(R_KIND), intent(in)    :: DNValue(:,:)
  real(R_KIND), intent(inout) :: Dgpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call AD_PolarInterpolation2D(DNValue(1:gsize(1),k), &
                                 & Dgpv(1:gsize(1),1,k))
    end do
  end if
end subroutine AD_NPInterpolation3D
!************************************************************************
subroutine TL_SPInterpolation3D(DSValue, Dgpv)
  real(R_KIND), intent(out)   :: DSValue(:,:)
  real(R_KIND), intent(in)    :: Dgpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call TL_PolarInterpolation2D(DSValue(1:gsize(1),k), &
                                 & Dgpv(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine TL_SPInterpolation3D
!************************************************************************
subroutine AD_SPInterpolation3D(DSValue, Dgpv)
  real(R_KIND), intent(in)    :: DSValue(:,:)
  real(R_KIND), intent(inout) :: Dgpv(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call AD_PolarInterpolation2D(DSValue(1:gsize(1),k), &
                                 & Dgpv(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine AD_SPInterpolation3D
!************************************************************************
subroutine TL_NPInterpolation3D_Wind(DNValueU, DNValueV, DgpvU, DgpvV)
  real(R_KIND), intent(out)   :: DNValueU(:,:), DNValueV(:,:)
  real(R_KIND), intent(in)    :: DgpvU(:,:,:), DgpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call TL_NPInterpolation2D_Wind(DNValueU(1:gsize(1),k), &
                                  & DNValueV(1:gsize(1),k), &
                                  & DgpvU(1:gsize(1),1,k), &
                                  & DgpvV(1:gsize(1),1,k))
    end do
  end if
end subroutine TL_NPInterpolation3D_Wind
!************************************************************************
subroutine AD_NPInterpolation3D_Wind(DNValueU, DNValueV, DgpvU, DgpvV)
  real(R_KIND), intent(in)    :: DNValueU(:,:), DNValueV(:,:)
  real(R_KIND), intent(inout) :: DgpvU(:,:,:), DgpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call AD_NPInterpolation2D_Wind(DNValueU(1:gsize(1),k), &
                                  & DNValueV(1:gsize(1),k), &
                                  & DgpvU(1:gsize(1),1,k), &
                                  & DgpvV(1:gsize(1),1,k))
    end do
  end if
end subroutine AD_NPInterpolation3D_Wind
!************************************************************************
subroutine TL_SPInterpolation3D_Wind(DSValueU, DSValueV, DgpvU, DgpvV)
  real(R_KIND), intent(out)   :: DSValueU(:,:), DSValueV(:,:)
  real(R_KIND), intent(in)    :: DgpvU(:,:,:), DgpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call TL_SPInterpolation2D_Wind(DSValueU(1:gsize(1),k), &
                                  & DSValueV(1:gsize(1),k), &
                                  & DgpvU(1:gsize(1),gsize(2),k), &
                                  & DgpvV(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine TL_SPInterpolation3D_Wind
!************************************************************************
subroutine AD_SPInterpolation3D_Wind(DSValueU, DSValueV, DgpvU, DgpvV)
  real(R_KIND), intent(in)    :: DSValueU(:,:), DSValueV(:,:)
  real(R_KIND), intent(inout) :: DgpvU(:,:,:), DgpvV(:,:,:)
!
  integer(4) :: k
!
  if (proj == "GS  ") then
    do k = 1, gsize(3)
      call AD_SPInterpolation2D_Wind(DSValueU(1:gsize(1),k), &
                                  & DSValueV(1:gsize(1),k), &
                                  & DgpvU(1:gsize(1),gsize(2),k), &
                                  & DgpvV(1:gsize(1),gsize(2),k))
    end do
  end if
end subroutine AD_SPInterpolation3D_Wind
!************************************************************************


end module interpolate_TLAD
