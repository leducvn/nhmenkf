module NodeHField_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl
   use interpolate, only : PosGrid, convert_LatLon_to_GridPos, convert_GridPos_offset
   use NodeInfo_class
   use NodeObsSpaceField_class
   use NodeObsValidField_class
   use NodeObsField_class
   use NodeObsControl_class
   use NodeHFieldCNV_class
   use NodeHFieldTC_class
   use NodeHFieldGNSS_class
   use NodeHFieldRAD_class
   use NodeControl_class
   use NodeProfileControl_class
   use NodeMPI
   implicit none
   !
   type NodeHField
      character(len=10) :: obstype, name
      integer :: nobs
      type(NodeHFieldCNV) :: cnv
      type(NodeHFieldTC) :: tc
      type(NodeHFieldGNSS) :: gnss
      type(NodeHFieldRAD) :: rad
   end type NodeHField
   !
   interface new
      module procedure new_NodeHField
   end interface
   interface destroy
      module procedure destroy_NodeHField
   end interface
   interface display
      module procedure display_NodeHField
   end interface
   interface get_name
      module procedure get_name_NodeHField
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHField
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHField1
      module procedure get_xyloc_NodeHField2
   end interface
   interface apply_Hlogp
      module procedure apply_Hlogp_NodeHField
   end interface
   interface apply_H
      module procedure apply_H_NodeHField1
      module procedure apply_H_NodeHField2
      module procedure apply_H_NodeHField3
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHField
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHField
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHField
   end interface
   !
contains
   !
   subroutine new_NodeHField(self, info, obsspace)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace
      !
      self%obstype = obsspace%obstype
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      if (trim(self%obstype) == 'CNV') then
         call new(self%cnv, info, obsspace%cnv)
      else if (trim(self%obstype) == 'TC') then
         call new(self%tc, info, obsspace%tc)
      else if (trim(self%obstype) == 'GNSS') then
         call new(self%gnss, info, obsspace%gnss)
      else if (trim(self%obstype) == 'RAD') then
         call new(self%rad, info, obsspace%rad)
      end if
      !
      return
   end subroutine new_NodeHField
   !
   !
   !
   subroutine destroy_NodeHField(self)
      implicit none
      type(NodeHField), intent(inout) :: self
      !
      if (trim(self%obstype) == 'CNV') then
         call destroy(self%cnv)
      else if (trim(self%obstype) == 'TC') then
         call destroy(self%tc)
      else if (trim(self%obstype) == 'GNSS') then
         call destroy(self%gnss)
      else if (trim(self%obstype) == 'RAD') then
         call destroy(self%rad)
      end if
      !
      return
   end subroutine destroy_NodeHField
   !
   !
   !
   subroutine display_NodeHField(self)
      implicit none
      type(NodeHField), intent(in) :: self
      !
      if (trim(self%obstype) == 'CNV') then
         call display(self%cnv)
      else if (trim(self%obstype) == 'TC') then
         call display(self%tc)
      else if (trim(self%obstype) == 'GNSS') then
         call display(self%gnss)
      else if (trim(self%obstype) == 'RAD') then
         call display(self%rad)
      end if
      !
      return
   end subroutine display_NodeHField
   !
   !
   !
   subroutine get_name_NodeHField(self, name)
      implicit none
      type(NodeHField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      if (trim(self%obstype) == 'CNV') then
         call get_name(self%cnv, name)
      else if (trim(self%obstype) == 'TC') then
         call get_name(self%tc, name)
      else if (trim(self%obstype) == 'GNSS') then
         call get_name(self%gnss, name)
      else if (trim(self%obstype) == 'RAD') then
         call get_name(self%rad, name)
      end if
      !
      return
   end subroutine get_name_NodeHField
   !
   !
   !
   subroutine get_nobs_NodeHField(self, nobs)
      implicit none
      type(NodeHField), intent(in) :: self
      integer, intent(out) :: nobs
      !
      if (trim(self%obstype) == 'CNV') then
         call get_nobs(self%cnv, nobs)
      else if (trim(self%obstype) == 'TC') then
         call get_nobs(self%tc, nobs)
      else if (trim(self%obstype) == 'GNSS') then
         call get_nobs(self%gnss, nobs)
      else if (trim(self%obstype) == 'RAD') then
         call get_nobs(self%rad, nobs)
      end if
      !
      return
   end subroutine get_nobs_NodeHField
   !
   !
   !
   subroutine get_xyloc_NodeHField1(self, xyloc)
      implicit none
      type(NodeHField), intent(in) :: self
      type(NodeObsField), intent(inout) :: xyloc
      !
      if (trim(self%obstype) == 'CNV') then
         call get_xyloc(self%cnv, xyloc)
      else if (trim(self%obstype) == 'TC') then
         call get_xyloc(self%tc, xyloc)
      else if (trim(self%obstype) == 'GNSS') then
         call get_xyloc(self%gnss, xyloc)
      else if (trim(self%obstype) == 'RAD') then
         call get_xyloc(self%rad, xyloc)
      end if
      !
      return
   end subroutine get_xyloc_NodeHField1
   !
   !
   !
   subroutine get_xyloc_NodeHField2(self, info, xyloc)
      implicit none
      type(NodeHField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsField), intent(inout) :: xyloc
      !
      if (trim(self%obstype) == 'CNV') then
         call get_xyloc(self%cnv, info, xyloc)
      else if (trim(self%obstype) == 'TC') then
         call get_xyloc(self%tc, info, xyloc)
      else if (trim(self%obstype) == 'GNSS') then
         call get_xyloc(self%gnss, info, xyloc)
      else if (trim(self%obstype) == 'RAD') then
         call get_xyloc(self%rad, info, xyloc)
      end if
      !
      return
   end subroutine get_xyloc_NodeHField2
   !
   !
   !
   subroutine apply_Hlogp_NodeHField(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_Hlogp(self%cnv, info, x, obsspace%cnv, valid, y)
      else if (trim(self%obstype) == 'TC') then
         !call apply_Hlogp(self%tc, info, x, obsspace%tc, valid, y)
      else if (trim(self%obstype) == 'GNSS') then
         !call apply_Hlogp(self%gnss, info, x, obsspace%gnss, valid, y)
      else if (trim(self%obstype) == 'RAD') then
         !call apply_Hlogp(self%rad, info, x, obsspace%rad, valid, y)
      end if
      !
      return
   end subroutine apply_Hlogp_NodeHField
   !
   !
   !
   subroutine apply_H_NodeHField1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_H(self%cnv, info, x, obsspace%cnv, valid, y)
      else if (trim(self%obstype) == 'TC') then
         call apply_H(self%tc, info, x, obsspace%tc, valid, y)
      else if (trim(self%obstype) == 'GNSS') then
         call apply_H(self%gnss, info, x, obsspace%gnss, valid, y)
      else if (trim(self%obstype) == 'RAD') then
         call apply_H(self%rad, info, x, obsspace%rad, valid, y)
      end if
      !
      return
   end subroutine apply_H_NodeHField1
   !
   !
   !
   subroutine apply_H_NodeHField2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_H(self%cnv, info, processed, k2ijt, x, obsspace%cnv, valid, y, nx, ny, nxyt)
      else if (trim(self%obstype) == 'TC') then
         call apply_H(self%tc, info, processed, k2ijt, x, obsspace%tc, valid, y, nx, ny, nxyt)
      else if (trim(self%obstype) == 'GNSS') then
         call apply_H(self%gnss, info, processed, k2ijt, x, obsspace%gnss, valid, y, nx, ny, nxyt)
      else if (trim(self%obstype) == 'RAD') then
         call apply_H(self%rad, info, processed, k2ijt, x, obsspace%rad, valid, y, nx, ny, nxyt)
      end if
      !
      return
   end subroutine apply_H_NodeHField2
   !
   !
   !
   subroutine apply_H_NodeHField3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_H(self%cnv, info, ip, jp, ijt2k, x, obsspace%cnv, valid, y, nt)
      else if (trim(self%obstype) == 'TC') then
         call apply_H(self%tc, info, ip, jp, ijt2k, x, obsspace%tc, valid, y, nt)
      else if (trim(self%obstype) == 'GNSS') then
         call apply_H(self%gnss, info, ip, jp, ijt2k, x, obsspace%gnss, valid, y, nt)
      else if (trim(self%obstype) == 'RAD') then
         call apply_H(self%rad, info, ip, jp, ijt2k, x, obsspace%rad, valid, y, nt)
      end if
      !
      return
   end subroutine apply_H_NodeHField3
   
   !
   !
   !
   subroutine initialize_DH_NodeHField(self, info, x, obsspace, valid, y, ne)
      implicit none
      type(NodeHField), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         !call initialize_DH(self%cnv, info, x, obsspace%cnv, valid, y, ne)
      else if (trim(self%obstype) == 'TC') then
         call initialize_DH(self%tc, info, x, obsspace%tc, valid, y, ne)
      else if (trim(self%obstype) == 'GNSS') then
         !call initialize_DH(self%gnss, info, x, obsspace%gnss, valid, y, ne)
      else if (trim(self%obstype) == 'RAD') then
         !call initialize_DH(self%rad, info, x, obsspace%rad, valid, y, ne)
      end if
      !
      return
   end subroutine initialize_DH_NodeHField
   !
   !
   !
   subroutine apply_DH_NodeHField(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_DH(self%cnv, info, xbck, x, obsspace%cnv, valid, y)
      else if (trim(self%obstype) == 'TC') then
         call apply_DH(self%tc, info, xbck, x, obsspace%tc, valid, y)
      else if (trim(self%obstype) == 'GNSS') then
         call apply_DH(self%gnss, info, xbck, x, obsspace%gnss, valid, y)
      else if (trim(self%obstype) == 'RAD') then
         call apply_DH(self%rad, info, xbck, x, obsspace%rad, valid, y)
      end if
      !
      return
   end subroutine apply_DH_NodeHField
   !
   !
   !
   subroutine apply_DHT_NodeHField(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      if (trim(self%obstype) == 'CNV') then
         call apply_DHT(self%cnv, info, xbck, obsspace%cnv, valid, y, x)
      else if (trim(self%obstype) == 'TC') then
         call apply_DHT(self%tc, info, xbck, obsspace%tc, valid, y, x)
      else if (trim(self%obstype) == 'GNSS') then
         call apply_DHT(self%gnss, info, xbck, obsspace%gnss, valid, y, x)
      else if (trim(self%obstype) == 'RAD') then
         call apply_DHT(self%rad, info, xbck, obsspace%rad, valid, y, x)
      end if
      !
      return
   end subroutine apply_DHT_NodeHField
   !
   !
   !
end module NodeHField_class
