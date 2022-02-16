module NodeObsSpaceField_class
! Author: Le Duc
! Created date: 13 Mar 2016
! Dec 10 2020: modify as a prototype for different kinds of observations
!              Now, polymorphism is applied for this class
   use variable, only : r_size, r_sngl
   use NodeInfo_class
   use NodeObsSpaceFieldCNV_class
   use NodeObsSpaceFieldTC_class
   use NodeObsSpaceFieldGNSS_class
   use NodeObsSpaceFieldRAD_class
   implicit none
   !
   type NodeObsSpaceField
      character(len=10) :: obstype, name
      integer :: nobs, nsubobs, nx, ny, nt
      integer, dimension(:,:,:), allocatable :: iobs, mobs
      type(NodeObsSpaceFieldCNV) :: cnv
      type(NodeObsSpaceFieldTC) :: tc
      type(NodeObsSpaceFieldGNSS) :: gnss
      type(NodeObsSpaceFieldRAD) :: rad
   end type NodeObsSpaceField
   !
   interface new
      module procedure new_NodeObsSpaceField
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceField
   end interface
   interface display
      module procedure display_NodeObsSpaceField
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceField
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceField
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceField
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceField1
      module procedure get_mobs_NodeObsSpaceField2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceField
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceField
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceField
   end interface
   !
contains
   !
   subroutine new_NodeObsSpaceField(self, name, filename, nx, ny, nt)
      implicit none
      type(NodeObsSpaceField), intent(inout) :: self
      character(len=*), intent(in) :: name, filename
      integer, intent(in) :: nx, ny, nt
      !
      self%name = name
      self%nx = nx
      self%ny = ny
      self%nt = nt
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	&     trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv' .or. &
	&     trim(self%name) == 'rvl') then
         self%obstype = 'CNV'
         call new(self%cnv, name, filename, nx, ny, nt)
      else if (trim(self%name) == 'tc') then
         self%obstype = 'TC'
         call new(self%tc, name, filename, nx, ny, nt)
      else if (trim(self%name) == 're' .or. trim(self%name) == 'ba') then
         self%obstype = 'GNSS'
         call new(self%gnss, name, filename, nx, ny, nt)
      else if (trim(self%name) == 'rad') then
         self%obstype = 'RAD'
         call new(self%rad, name, filename, nx, ny, nt)
      end if
      !
      return
   end subroutine new_NodeObsSpaceField
   !
   !
   !
   subroutine destroy_NodeObsSpaceField(self)
      implicit none
      type(NodeObsSpaceField), intent(inout) :: self
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
      if (allocated(self%iobs)) deallocate(self%iobs)
      if (allocated(self%mobs)) deallocate(self%mobs)
      !
      return
   end subroutine destroy_NodeObsSpaceField
   !
   !
   !
   subroutine display_NodeObsSpaceField(self)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
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
   end subroutine display_NodeObsSpaceField
   !
   !
   !
   subroutine set_metadata(self)
      implicit none
      type(NodeObsSpaceField), intent(inout) :: self
      !
      allocate(self%iobs(self%nx,self%ny,self%nt))
      allocate(self%mobs(self%nx,self%ny,self%nt))
      if (trim(self%obstype) == 'CNV') then
         self%nobs = self%cnv%nobs
	 self%nsubobs = self%cnv%nsubobs
         self%iobs(:,:,:) = self%cnv%iobs(:,:,:)
	 self%mobs(:,:,:) = self%cnv%mobs(:,:,:)
      else if (trim(self%obstype) == 'TC') then
         self%nobs = self%tc%nobs
	 self%nsubobs = self%tc%nsubobs
         self%iobs(:,:,:) = self%tc%iobs(:,:,:)
	 self%mobs(:,:,:) = self%tc%mobs(:,:,:)
      else if (trim(self%obstype) == 'GNSS') then
         self%nobs = self%gnss%nobs
	 self%nsubobs = self%gnss%nsubobs
         self%iobs(:,:,:) = self%gnss%iobs(:,:,:)
	 self%mobs(:,:,:) = self%gnss%mobs(:,:,:)
      else if (trim(self%obstype) == 'RAD') then
         self%nobs = self%rad%nobs
	 self%nsubobs = self%rad%nsubobs
         self%iobs(:,:,:) = self%rad%iobs(:,:,:)
	 self%mobs(:,:,:) = self%rad%mobs(:,:,:)
      end if
      !
      return
   end subroutine set_metadata
   !
   !
   !
   subroutine get_name_NodeObsSpaceField(self, name)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
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
   end subroutine get_name_NodeObsSpaceField
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceField(self, nobs)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
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
   end subroutine get_nobs_NodeObsSpaceField
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceField(self, nsubobs)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
      integer, intent(out) :: nsubobs
      !
      if (trim(self%obstype) == 'CNV') then
	 call get_nsubobs(self%cnv, nsubobs)
      else if (trim(self%obstype) == 'TC') then
	 call get_nsubobs(self%tc, nsubobs)
      else if (trim(self%obstype) == 'GNSS') then
	 call get_nsubobs(self%gnss, nsubobs)
      else if (trim(self%obstype) == 'RAD') then
	 call get_nsubobs(self%rad, nsubobs)
      end if
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceField
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceField1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      !
      if (trim(self%obstype) == 'CNV') then
	 call get_mobs(self%cnv, i, j, mobs)
      else if (trim(self%obstype) == 'TC') then
	 call get_mobs(self%tc, i, j, mobs)
      else if (trim(self%obstype) == 'GNSS') then
	 call get_mobs(self%gnss, i, j, mobs)
      else if (trim(self%obstype) == 'RAD') then
	 call get_mobs(self%rad, i, j, mobs)
      end if
      !
      return
   end subroutine get_mobs_NodeObsSpaceField1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceField2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      !
      if (trim(self%obstype) == 'CNV') then
	 call get_mobs(self%cnv, i, j, it, mobs)
      else if (trim(self%obstype) == 'TC') then
	 call get_mobs(self%tc, i, j, it, mobs)
      else if (trim(self%obstype) == 'GNSS') then
	 call get_mobs(self%gnss, i, j, it, mobs)
      else if (trim(self%obstype) == 'RAD') then
	 call get_mobs(self%rad, i, j, it, mobs)
      end if
      !
      return
   end subroutine get_mobs_NodeObsSpaceField2
   !
   !
   !
   subroutine read_obs_NodeObsSpaceField(self, myid)
      implicit none
      type(NodeObsSpaceField), intent(inout) :: self
      integer, intent(in) :: myid
      !
      if (trim(self%obstype) == 'CNV') then
         call read_obs(self%cnv, myid)
      else if (trim(self%obstype) == 'TC') then
         call read_obs(self%tc, myid)
      else if (trim(self%obstype) == 'GNSS') then
         call read_obs(self%gnss, myid)
      else if (trim(self%obstype) == 'RAD') then
         call read_obs(self%rad, myid)
      end if
      !
      return
   end subroutine read_obs_NodeObsSpaceField
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceField(self, info, local_object)
      ! scatter a big obsspace from Node 0 to all nodes.
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(inout) :: local_object
      !
      if (trim(self%obstype) == 'CNV') then
         call scatter_obs(self%cnv, info, local_object%cnv)
      else if (trim(self%obstype) == 'TC') then
         call scatter_obs(self%tc, info, local_object%tc)
      else if (trim(self%obstype) == 'GNSS') then
	 call scatter_obs(self%gnss, info, local_object%gnss)
      else if (trim(self%obstype) == 'RAD') then
         call scatter_obs(self%rad, info, local_object%rad)
      end if
      call set_metadata(local_object)
      !
      return
   end subroutine scatter_obs_NodeObsSpaceField
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceField(self, info, global_object)
      ! this subroutine is only for completion and we realy do not use this
      implicit none
      type(NodeObsSpaceField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(inout) :: global_object
      !
      if (trim(self%obstype) == 'CNV') then
         call gather_obs(self%cnv, info, global_object%cnv)
      else if (trim(self%obstype) == 'TC') then
         call gather_obs(self%tc, info, global_object%tc)
      else if (trim(self%obstype) == 'GNSS') then
         call gather_obs(self%gnss, info, global_object%gnss)
      else if (trim(self%obstype) == 'RAD') then
         call gather_obs(self%rad, info, global_object%rad)
      end if
      !
      return
   end subroutine gather_obs_NodeObsSpaceField
   !
   !
   !
end module NodeObsSpaceField_class
