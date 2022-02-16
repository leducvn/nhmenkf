module NodeObsField_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsSpaceField_class
   use NodeObsValidField_class
   use NodeMPI
   implicit none
   !
   type NodeObsField
      character(len=10) :: name
      integer :: nobs, nsubobs, nx, ny
      real(r_size), dimension(:,:), allocatable :: field
   end type NodeObsField
   !
   interface new
      module procedure new_NodeObsField1
      module procedure new_NodeObsField2
   end interface
   interface destroy
      module procedure destroy_NodeObsField
   end interface
   interface display
      module procedure display_NodeObsField
   end interface
   interface assignment(=)
      module procedure copy_NodeObsField
   end interface
   interface get_name
      module procedure get_name_NodeObsField
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsField
   end interface
   interface extract_y
      module procedure extract_y_NodeObsField0
      module procedure extract_y_NodeObsField1
      module procedure extract_y_NodeObsField2
      module procedure extract_y_NodeObsField3
      module procedure extract_y_NodeObsField4
   end interface
   interface set_obs
      module procedure set_obs_NodeObsField
   end interface
   interface set_error
      module procedure set_error_NodeObsField
   end interface
   interface set_field
      module procedure set_field_NodeObsField1
      module procedure set_field_NodeObsField2
      module procedure set_field_NodeObsField3
      module procedure set_field_NodeObsField4
      module procedure set_field_NodeObsField5
      module procedure set_field_NodeObsField6
   end interface
   interface add_obs
      module procedure add_obs_NodeObsField1
      module procedure add_obs_NodeObsField2
      module procedure add_obs_NodeObsField3
   end interface
   interface subtract_obs
      module procedure subtract_obs_NodeObsField1
      module procedure subtract_obs_NodeObsField2
      module procedure subtract_obs_NodeObsField3
   end interface
   interface power_obs
      module procedure power_obs_NodeObsField
   end interface
   interface loga10_obs
      module procedure loga10_obs_NodeObsField
   end interface
   interface multiply_obs
      module procedure multiply_obs_NodeObsField1
      module procedure multiply_obs_NodeObsField2
   end interface
   interface divide_obs
      module procedure divide_obs_NodeObsField1
      module procedure divide_obs_NodeObsField2
      module procedure divide_obs_NodeObsField3
      module procedure divide_obs_NodeObsField4
      module procedure divide_obs_NodeObsField5
      module procedure divide_obs_NodeObsField6
   end interface
   interface compute_normsquare
      module procedure compute_normsquare_NodeObsField
   end interface
   interface innerproduct
      module procedure innerproduct_NodeObsField1
      module procedure innerproduct_NodeObsField2
      module procedure innerproduct_NodeObsField3
   end interface
   interface qccheck
      module procedure qccheck_NodeObsField1
      module procedure qccheck_NodeObsField2
   end interface
   interface stdcheck
      module procedure stdcheck_NodeObsField
   end interface
   interface random_obs
      module procedure random_obs_NodeObsField1
      module procedure random_obs_NodeObsField2
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsField
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsField
   end interface
   interface broadcast_ensobs
      module procedure broadcast_ensobs_NodeObsField1
      module procedure broadcast_ensobs_NodeObsField2
   end interface
   interface allreduce_ensobs
      module procedure allreduce_ensobs_NodeObsField1
      module procedure allreduce_ensobs_NodeObsField2
      module procedure allreduce_ensobs_NodeObsField3
   end interface
   !
contains
   !
   !
   !
   subroutine new_NodeObsField1(self, obsspace)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: obsspace
      !
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      self%nsubobs = obsspace%nsubobs
      self%nx = obsspace%nx
      self%ny = obsspace%ny
      if (self%nobs > 0) then
	 allocate(self%field(self%nobs,self%nsubobs))
	 self%field = 0.d0
      end if
      !
      return
   end subroutine new_NodeObsField1
   !
   !
   !
   subroutine new_NodeObsField2(self, obsspace, nsubobs)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: nsubobs
      !
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      self%nsubobs = nsubobs
      self%nx = obsspace%nx
      self%ny = obsspace%ny
      if (self%nobs > 0) then
         allocate(self%field(self%nobs,self%nsubobs))
         self%field = 0.d0
      end if
      !
      return
   end subroutine new_NodeObsField2
   !
   !
   !
   subroutine destroy_NodeObsField(self)
      implicit none
      type(NodeObsField), intent(inout) :: self
      !
      if (allocated(self%field)) deallocate(self%field)
      !
      return
   end subroutine destroy_NodeObsField
   !
   !
   !
   subroutine display_NodeObsField(self)
      implicit none
      type(NodeObsField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsField
   !
   !
   !
   subroutine copy_NodeObsField(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = object%field(:,:)
      !
      return
   end subroutine copy_NodeObsField
   !
   !
   !
   subroutine get_name_NodeObsField(self, name)
      implicit none
      type(NodeObsField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsField
   !
   !
   !
   subroutine get_nobs_NodeObsField(self, nobs)
      implicit none
      type(NodeObsField), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsField
   !
   !
   !
   subroutine extract_y_NodeObsField0(self, obsspace, valid, mobs, y, mmax)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: mmax
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      integer, intent(inout) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: iobs, jobs, nsubobs
      !
      do iobs = 1, self%nobs
         if (trim(obsspace%obstype) == 'RAD') then
            nsubobs = obsspace%rad%nchannel(iobs)
         else
            nsubobs = self%nsubobs
         end if
         do jobs = 1, nsubobs
            if (valid%field(iobs,jobs) == 0) cycle
            mobs = mobs + 1
            y(mobs) = self%field(iobs,jobs)
         end do
      end do
      !
      return
   end subroutine extract_y_NodeObsField0
   !
   !
   !
   subroutine extract_y_NodeObsField1(self, i, j, obsspace, valid, mobs, y, mmax)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j, mmax
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      integer, intent(inout) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: it, iobs1, iobs2, iobs, jobs, nsubobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    do jobs = 1, nsubobs
	       if (valid%field(iobs,jobs) == 0) cycle
	       mobs = mobs + 1
	       y(mobs) = self%field(iobs,jobs)
	    end do
	 end do
      end do
      !
      return
   end subroutine extract_y_NodeObsField1
   !
   !
   !
   subroutine extract_y_NodeObsField2(self, i, j, ip, jp, obsspace, valid, mobs, y, np, mmax)
      use variable, only : hscale
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      integer, intent(inout) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: k, ix, jy, it, iobs1, iobs2, iobs, jobs, nsubobs
      real(r_size) :: dh, hweight
      !
      do k = 1, np
         ix = ip(k); jy = jp(k)
         dh = sqrt(1.d0*(ix-i)**2+1.d0*(jy-j)**2)
         !if (dh > hscale) cycle
	 call GaspariCohn(dh, hscale, hweight)
	 do it = 1, obsspace%nt
	    if (obsspace%mobs(ix,jy,it) == 0) cycle
	    iobs1 = obsspace%iobs(ix,jy,it)
	    iobs2 = obsspace%iobs(ix,jy,it) + obsspace%mobs(ix,jy,it) - 1
	    do iobs = iobs1, iobs2
	       if (trim(obsspace%obstype) == 'RAD') then
	          nsubobs = obsspace%rad%nchannel(iobs)
	       else
	          nsubobs = self%nsubobs
	       end if
	       do jobs = 1, nsubobs
		  if (valid%field(iobs,jobs) == 0) cycle
		  mobs = mobs + 1
		  y(mobs) = hweight*self%field(iobs,jobs)
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine extract_y_NodeObsField2
   !
   !
   !
   subroutine extract_y_NodeObsField3(self, i, j, ip, jp, obsspace, valid, xyloc, mobs, y, np, mmax)
      use variable, only : hscale
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: xyloc
      integer, intent(inout) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: k, ix, jy, it, iobs1, iobs2, iobs, jobs, nsubobs
      real(r_size) :: dh, hweight
      !
      do k = 1, np
         ix = ip(k); jy = jp(k)
	 do it = 1, obsspace%nt
	    if (obsspace%mobs(ix,jy,it) == 0) cycle
	    iobs1 = obsspace%iobs(ix,jy,it)
	    iobs2 = obsspace%iobs(ix,jy,it) + obsspace%mobs(ix,jy,it) - 1
	    do iobs = iobs1, iobs2
	       dh = sqrt(1.d0*(xyloc%field(iobs,1)-i)**2+1.d0*(xyloc%field(iobs,2)-j)**2)
	       call GaspariCohn(dh, hscale, hweight)
	       if (trim(obsspace%obstype) == 'RAD') then
	          nsubobs = obsspace%rad%nchannel(iobs)
	       else
	          nsubobs = self%nsubobs
	       end if
	       do jobs = 1, nsubobs
		  if (valid%field(iobs,jobs) == 0) cycle
		  mobs = mobs + 1
		  y(mobs) = hweight*self%field(iobs,jobs)
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine extract_y_NodeObsField3
   !
   !
   !
   subroutine extract_y_NodeObsField4(self, i, j, ip, jp, obsspace, valid, xyloc, vloc, mobs, y, np, mmax)
      use variable, only : hscale
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: xyloc
      type(NodeObsField), intent(in) :: vloc
      integer, intent(inout) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: k, ix, jy, it, iobs1, iobs2, iobs, jobs, nsubobs
      real(r_size) :: dh, hweight
      !
      do k = 1, np
         ix = ip(k); jy = jp(k)
	 do it = 1, obsspace%nt
	    if (obsspace%mobs(ix,jy,it) == 0) cycle
	    iobs1 = obsspace%iobs(ix,jy,it)
	    iobs2 = obsspace%iobs(ix,jy,it) + obsspace%mobs(ix,jy,it) - 1
	    do iobs = iobs1, iobs2
	       dh = sqrt(1.d0*(xyloc%field(iobs,1)-i)**2+1.d0*(xyloc%field(iobs,2)-j)**2)
	       call GaspariCohn(dh, hscale, hweight)
	       if (trim(obsspace%obstype) == 'RAD') then
	          nsubobs = obsspace%rad%nchannel(iobs)
	       else
	          nsubobs = self%nsubobs
	       end if
	       do jobs = 1, nsubobs
		  if (valid%field(iobs,jobs) == 0) cycle
		  mobs = mobs + 1
		  y(mobs) = hweight*vloc%field(iobs,jobs)*self%field(iobs,jobs)
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine extract_y_NodeObsField4
   !
   !
   !
   subroutine set_obs_NodeObsField(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: object
      !
      if (self%nobs == 0) return
      if (trim(object%obstype) == 'CNV') then
         self%field(:,:) = object%cnv%obs(:,:)
      else if (trim(object%obstype) == 'TC') then
         self%field(:,:) = object%tc%obs(:,:)
      else if (trim(object%obstype) == 'GNSS') then
         self%field(:,:) = object%gnss%obs(:,:)
      else if (trim(object%obstype) == 'RAD') then
         self%field(:,:) = object%rad%obs(:,:)
      end if
      !
      return
   end subroutine set_obs_NodeObsField
   !
   !
   !
   subroutine set_error_NodeObsField(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: object
      !
      if (self%nobs == 0) return
      if (trim(object%obstype) == 'CNV') then
         self%field(:,:) = object%cnv%error(:,:)
      else if (trim(object%obstype) == 'TC') then
         self%field(:,:) = object%tc%error(:,:)
      else if (trim(object%obstype) == 'GNSS') then
         self%field(:,:) = object%gnss%error(:,:)
      else if (trim(object%obstype) == 'RAD') then
         self%field(:,:) = object%rad%error(:,:)
      end if
      !
      return
   end subroutine set_error_NodeObsField
   !
   !
   !
   subroutine set_field_NodeObsField1(self, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:) = const
      !
      return
   end subroutine set_field_NodeObsField1
   !
   !
   !
   subroutine set_field_NodeObsField2(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = object%field(:,:)
      !
      return
   end subroutine set_field_NodeObsField2
   !
   !
   !
   subroutine set_field_NodeObsField3(self, i, j, obsspace, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: it, iobs1, iobs2, iobs, nsubobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = const
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsField3
   !
   !
   !
   subroutine set_field_NodeObsField4(self, i, j, obsspace, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: it, iobs1, iobs2, iobs, nsubobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsField4
   !
   !
   !
   subroutine set_field_NodeObsField5(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = const
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsField5
   !
   !
   !
   subroutine set_field_NodeObsField6(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsField6
   !
   !
   !
   subroutine add_obs_NodeObsField1(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = self%field(:,:) + object%field(:,:)
      !
      return
   end subroutine add_obs_NodeObsField1
   !
   !
   !
   subroutine add_obs_NodeObsField2(self, i, j, obsspace, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: it, iobs1, iobs2, iobs, nsubobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs) + object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine add_obs_NodeObsField2
   !
   !
   !
   subroutine add_obs_NodeObsField3(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs) + object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine add_obs_NodeObsField3
   !
   !
   !
   subroutine subtract_obs_NodeObsField1(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = self%field(:,:) - object%field(:,:)
      !
      return
   end subroutine subtract_obs_NodeObsField1
   !
   !
   !
   subroutine subtract_obs_NodeObsField2(self, i, j, obsspace, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: nt, it, iobs1, iobs2, iobs, nsubobs
      !
      nt = obsspace%nt
      do it = 1, nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs) - object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine subtract_obs_NodeObsField2
   !
   !
   !
   subroutine subtract_obs_NodeObsField3(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs) - object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine subtract_obs_NodeObsField3
   !
   !
   !
   subroutine power_obs_NodeObsField(self, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:) = self%field(:,:)**const
      !
      return
   end subroutine power_obs_NodeObsField
   !
   !
   !
   subroutine loga10_obs_NodeObsField(self)
      implicit none
      type(NodeObsField), intent(inout) :: self
      !
      if (self%nobs == 0) return
      self%field(:,:) = log10(self%field(:,:))
      !
      return
   end subroutine loga10_obs_NodeObsField
   !
   !
   !
   subroutine multiply_obs_NodeObsField1(self, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:) = const*self%field(:,:)
      !
      return
   end subroutine multiply_obs_NodeObsField1
   !
   !
   !
   subroutine multiply_obs_NodeObsField2(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = object%field(:,:)*self%field(:,:)
      !
      return
   end subroutine multiply_obs_NodeObsField2
   !
   !
   !
   subroutine divide_obs_NodeObsField1(self, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:) = self%field(:,:)/const
      !
      return
   end subroutine divide_obs_NodeObsField1
   !
   !
   !
   subroutine divide_obs_NodeObsField2(self, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      !
      if (self%nobs == 0) return
      where (object%field < 1.e-12)
         self%field = 0.d0
      else where
         self%field = self%field/object%field
      end where
      !
      return
   end subroutine divide_obs_NodeObsField2
   !
   !
   !
   subroutine divide_obs_NodeObsField3(self, i, j, obsspace, const)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: nt, it, iobs1, iobs2, iobs, nsubobs
      !
      nt = obsspace%nt
      do it = 1, nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs)/const
	 end do
      end do
      !
      return
   end subroutine divide_obs_NodeObsField3
   !
   !
   !
   subroutine divide_obs_NodeObsField4(self, i, j, obsspace, object)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: nt, it, iobs1, iobs2, iobs, nsubobs
      !
      nt = obsspace%nt
      do it = 1, nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    where (object%field(iobs,1:nsubobs) < 1.e-12)
	       self%field(iobs,1:nsubobs) = 0.d0
	    else where
	       self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs)/object%field(iobs,1:nsubobs)
	    end where
	 end do	 
      end do
      !
      return
   end subroutine divide_obs_NodeObsField4
   !
   !
   !
   subroutine divide_obs_NodeObsField5(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs)/const
	 end do
      end do
      !
      return
   end subroutine divide_obs_NodeObsField5
   !
   !
   !
   subroutine divide_obs_NodeObsField6(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    where (object%field(iobs,1:nsubobs) < 1.e-12)
	       self%field(iobs,1:nsubobs) = 0.d0
	    else where
	       self%field(iobs,1:nsubobs) = self%field(iobs,1:nsubobs)/object%field(iobs,1:nsubobs)
	    end where
	 end do	 
      end do
      !
      return
   end subroutine divide_obs_NodeObsField6
   !
   !
   !
   subroutine compute_normsquare_NodeObsField(self, normsquare)
      implicit none
      type(NodeObsField), intent(in) :: self
      !type(NodeObsValidField), intent(in) :: valid
      real(r_size), intent(out) :: normsquare
      !integer :: iobs, j
      !
      normsquare = 0.d0
      if (self%nobs == 0) return
      !normsquare = 0.d0
      !do iobs = 1, self%nobs
         !do j = 1, self%nsubobs
            !if (valid%field(iobs,j) == 1) normsquare = normsquare + self%field(iobs,j)**2
	 !end do
      !end do
      ! all invalid entries are zero
      normsquare = sum(self%field(:,:)**2)
      !
      return
   end subroutine compute_normsquare_NodeObsField
   !
   !
   !
   subroutine innerproduct_NodeObsField1(self, object, product)
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeObsField), intent(in) :: object
      real(r_size), intent(out) :: product
      !
      product = 0.d0
      if (self%nobs == 0) return
      product = sum(self%field(:,:)*object%field(:,:))
      !
      return
   end subroutine innerproduct_NodeObsField1
   !
   !
   !
   subroutine innerproduct_NodeObsField2(self, object, valid, product)
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeObsField), intent(in) :: object
      type(NodeObsValidField), intent(in) :: valid
      real(r_size), intent(out) :: product
      integer :: iobs, j
      !
      product = 0.d0
      if (self%nobs == 0) return
      do iobs = 1, self%nobs
         do j = 1, self%nsubobs
            if (valid%field(iobs,j) == 1) product = product + self%field(iobs,j)*object%field(iobs,j)
         end do
      end do
      !
      return
   end subroutine innerproduct_NodeObsField2
   !
   !
   !
   subroutine innerproduct_NodeObsField3(self, object, obsspace, valid, product)
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeObsField), intent(in) :: object
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      real(r_size), intent(out) :: product
      integer :: iobs, jobs, nsubobs
      !
      product = 0.d0
      if (self%nobs == 0) return
      do iobs = 1, self%nobs
         if (trim(obsspace%obstype) == 'RAD') then
	    nsubobs = obsspace%rad%nchannel(iobs)
	 else
	    nsubobs = self%nsubobs
	 end if
         do jobs = 1, nsubobs
            if (valid%field(iobs,jobs) == 1) product = product + self%field(iobs,jobs)*object%field(iobs,jobs)
         end do
      end do
      !
      return
   end subroutine innerproduct_NodeObsField3
   !
   !
   !
   subroutine qccheck_cnv1(self, threshold, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidField), intent(inout) :: valid
      integer :: iobs
      !
      do iobs = 1, self%nobs
	 if (abs(self%field(iobs,1)) > threshold) then
	    valid%field(iobs,1) = 0
	 else
	    valid%field(iobs,1) = 1
	 end if
      end do
      !
      return
   end subroutine qccheck_cnv1
   !
   !
   !
   subroutine qccheck_cnv2(self, threshold, i, j, obsspace, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      integer :: it, iobs1, iobs2, iobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (abs(self%field(iobs,1)) > threshold) then
	       valid%field(iobs,1) = 0
	    else
	       valid%field(iobs,1) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine qccheck_cnv2
   !
   !
   !
   subroutine qccheck_tc1(self, threshold, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidField), intent(inout) :: valid
      integer :: iobs
      !
      do iobs = 1, self%nobs
         !if (sqrt(self%field(iobs,1)**2+self%field(iobs,2)**2) > threshold) then
         if (sqrt(self%field(iobs,1)**2+self%field(iobs,2)**2) > 20.) then
	    valid%field(iobs,:) = 0
	 else
	    valid%field(iobs,1:2) = 1
	    if (abs(self%field(iobs,3)) > threshold) then
	       valid%field(iobs,3) = 1
	    else
	       valid%field(iobs,3) = 1
	    end if
	 end if
      end do
      !
      return
   end subroutine qccheck_tc1
   !
   !
   !
   subroutine qccheck_tc2(self, threshold, i, j, obsspace, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      integer :: it, iobs1, iobs2, iobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (sqrt(self%field(iobs,1)**2+self%field(iobs,2)**2) > 20.) then
	       valid%field(iobs,:) = 0
	    else
	       valid%field(iobs,1:2) = 1
	       if (abs(self%field(iobs,3)) > threshold) then
	          valid%field(iobs,3) = 1
	       else
	          valid%field(iobs,3) = 1
	       end if
	    end if
	 end do
      end do
      !
      return
   end subroutine qccheck_tc2
   !
   !
   !
   subroutine qccheck_gnss1(self, threshold, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidField), intent(inout) :: valid
      integer :: iobs
      !
      do iobs = 1, self%nobs
	 if (abs(self%field(iobs,1)) > threshold) then
	    valid%field(iobs,1) = 0
	 else
	    valid%field(iobs,1) = 1
	 end if
      end do
      !
      return
   end subroutine qccheck_gnss1
   !
   !
   !
   subroutine qccheck_gnss2(self, threshold, i, j, obsspace, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      integer :: it, iobs1, iobs2, iobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (abs(self%field(iobs,1)) > threshold) then
	       valid%field(iobs,1) = 0
	    else
	       valid%field(iobs,1) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine qccheck_gnss2
   !
   !
   !
   subroutine qccheck_rad1(self, threshold, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidField), intent(inout) :: valid
      integer :: iobs, jobs
      !
      do iobs = 1, self%nobs
         do jobs = 1, self%nsubobs
	    if (abs(self%field(iobs,jobs)) > threshold) then
	       valid%field(iobs,jobs) = 0
	    else
	       valid%field(iobs,jobs) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine qccheck_rad1
   !
   !
   !
   subroutine qccheck_rad2(self, threshold, i, j, obsspace, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      integer :: it, iobs1, iobs2, iobs, jobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    do jobs = 1, obsspace%rad%nchannel(iobs)
	       if (abs(self%field(iobs,jobs)) > threshold) then
	          valid%field(iobs,jobs) = 0
	       else
	          valid%field(iobs,jobs) = 1
	       end if
	    end do
	 end do
      end do
      !
      return
   end subroutine qccheck_rad2
   !
   !
   !
   subroutine qccheck_NodeObsField1(self, threshold, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidField), intent(inout) :: valid
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv' .or. &
	& trim(self%name) == 'rvl') then
	 call qccheck_cnv1(self, threshold, valid)
      else if (trim(self%name) == 'tc') then
	 call qccheck_tc1(self, threshold, valid)
      else if (trim(self%name) == 're' .or. trim(self%name) == 'ba') then
	 call qccheck_gnss1(self, threshold, valid)
      else if (trim(self%name) == 'rad') then
	 call qccheck_rad1(self, threshold, valid)
      end if
      !
      return
   end subroutine qccheck_NodeObsField1
   !
   !
   !
   subroutine qccheck_NodeObsField2(self, threshold, i, j, obsspace, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv' .or. &
	& trim(self%name) == 'rvl') then
	 call qccheck_cnv2(self, threshold, i, j, obsspace, valid)
      else if (trim(self%name) == 'tc') then
	 call qccheck_tc2(self, threshold, i, j, obsspace, valid)
      else if (trim(self%name) == 're' .or. trim(self%name) == 'ba') then
	 call qccheck_gnss2(self, threshold, i, j, obsspace, valid)
      else if (trim(self%name) == 'rad') then
	 call qccheck_rad2(self, threshold, i, j, obsspace, valid)
      end if
      !
      return
   end subroutine qccheck_NodeObsField2
   !
   !
   !
   subroutine stdcheck_NodeObsField(self, valid)
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeObsValidField), intent(inout) :: valid
      integer :: iobs, jobs
      !
      do iobs = 1, self%nobs
         do jobs = 1, self%nsubobs
	    if (self%field(iobs,jobs) < 1.e-12) then
	       valid%field(iobs,jobs) = 0
	    else
	       valid%field(iobs,jobs) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine stdcheck_NodeObsField
   !
   !
   !
   subroutine random_obs_NodeObsField1(self)
      use random, only : random_normal
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer :: iobs, jobs
      !
      do iobs = 1, self%nobs
         do jobs = 1, self%nsubobs
            self%field(iobs,jobs) = 1.d0*random_normal()
         end do
      end do
      !
      return
   end subroutine random_obs_NodeObsField1
   !
   !
   !
   subroutine random_obs_NodeObsField2(self, valid)
      use random, only : random_normal
      implicit none
      type(NodeObsField), intent(inout) :: self
      type(NodeObsValidField), intent(in) :: valid
      integer :: iobs, jobs
      !
      do iobs = 1, self%nobs
         do jobs = 1, self%nsubobs
            if (valid%field(iobs,jobs) == 1) then
	       self%field(iobs,jobs) = 1.d0*random_normal()
	    else
	       self%field(iobs,jobs) = 0.d0
	    end if
         end do
      end do
      !
      return
   end subroutine random_obs_NodeObsField2
   !
   !
   !
   subroutine scatter_obs_NodeObsField(self, info, obsspace, local_obsspace, local_object)
      ! scatter a big obs from Node 0 to all nodes.
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace, local_obsspace
      type(NodeObsField), intent(inout) :: local_object
      !
      integer :: myid, myidx, myidy, nproc, nxpe, nype, nx0, ny0, nx, ny, nt
      integer :: iobs1, iobs2, iobs3, iobs4, i, j, it
      integer :: is, ie, js, je, dis, die, djs, dje
      integer, dimension(:,:,:), allocatable :: iobs
      real(r_size), dimension(:,:), allocatable :: field
      !
      call get_nxpe(info, nxpe); call get_nype(info, nype)
      call get_myidx(info, myidx); call get_myidy(info, myidy)
      myid = myidy*nxpe+myidx
      nx0 = obsspace%nx
      ny0 = obsspace%ny
      nx = local_obsspace%nx
      ny = local_obsspace%ny
      nt = local_obsspace%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      allocate(iobs(nx0,ny0,nt))
      allocate(field(obsspace%nobs,obsspace%nsubobs))
      !
      if (myid == 0) then
	 iobs(:,:,:) = obsspace%iobs
         field(:,:) = self%field(:,:)
      end if
      call int_scatter(info, 'xy', iobs, 0, nx0, ny0, nt)
      call broadcast2D('xy', field, 0, obsspace%nobs, obsspace%nsubobs)
      if (local_obsspace%nobs > 0) then
	 do it = 1, nt
	    do j = js, je
	       do i = is, ie
	          if (local_obsspace%mobs(i,j,it) > 0) then
		     iobs1 = local_obsspace%iobs(i,j,it)
		     iobs2 = iobs1 + local_obsspace%mobs(i,j,it) - 1
		     iobs3 = iobs(i,j,it)
		     iobs4 = iobs3 + local_obsspace%mobs(i,j,it) - 1
		     local_object%field(iobs1:iobs2,:) = field(iobs3:iobs4,:)
		  end if
	       end do
	    end do
	 end do
      end if
      deallocate(iobs, field)
      !
      return
   end subroutine scatter_obs_NodeObsField
   !
   !
   !
   subroutine gather_obs_NodeObsField(self, info, obsspace, global_obsspace, global_object)
      ! gather all obs from all nodes into one big obs. This big obs exists in all nodes.
      implicit none
      type(NodeObsField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace, global_obsspace
      type(NodeObsField), intent(inout) :: global_object
      !
      integer :: istat(MPI_STATUS_SIZE)
      integer :: nproc, myid, nx0, ny0, nx, ny, nt, nobs, nmax
      integer :: id, ierror, iobs1, iobs2, iobs3, iobs4, i, j, it
      integer :: idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myidx, myidy, myide, nxpe, nype, nepe, ide, source, destination
      integer, dimension(:,:,:), allocatable :: iobs
      real(kind=r_size), dimension(:,:), allocatable :: field
      !
      call get_myid(info, myid, nproc)
      call get_myidx(info, myidx); call get_nxpe(info, nxpe)
      call get_myidy(info, myidy); call get_nype(info, nype)
      call get_myide(info, myide); call get_nepe(info, nepe)
      nx0 = global_obsspace%nx
      ny0 = global_obsspace%ny
      nx = obsspace%nx
      ny = obsspace%ny
      nt = obsspace%nt
      !
      allocate(iobs(nx0,ny0,nt))
      iobs(1:nx,1:ny,:) = obsspace%iobs
      call int_gather(info, 'xy', iobs, 0, nx0, ny0, nt)
      if (myidx == 0 .and. myidy == 0) then
         if (myide == 0) then
            do ide = 1, nepe-1
               destination = ide*nxpe*nype
               call MPI_SEND(iobs, nx0*ny0*nt, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
            end do
         else
            source = 0
            call MPI_RECV(iobs, nx0*ny0*nt, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
         end if
      end if
      !
      if (myidx == 0 .and. myidy == 0) then
	 nmax = 0
	 do id = 0, nxpe*nype-1
	    call get_info(info, id, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nmax = max(nmax, sum(global_obsspace%mobs(is+dis:ie-die,js+djs:je-dje,:)))
	 end do
	 allocate(field(nmax,self%nsubobs))
	 !
	 do id = 0, nxpe*nype-1
	    source = myide*nxpe*nype + id
	    if (id == 0) then
	       nobs = self%nobs
	    else
	       call MPI_RECV(nobs, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
	    end if
	    if (nobs == 0) cycle
	    if (id == 0) then
	       field(1:nobs,:) = self%field
	    else
	       call MPI_RECV(field(1:nobs,:), nobs*self%nsubobs, r_type, source, 0, MPI_COMM_WORLD, istat, ierror)
	    end if
	    !
	    call get_info(info, id, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    do it = 1, nt
	       do j = js+djs, je-dje
		  do i = is+dis, ie-die
		     if (global_obsspace%mobs(i,j,it) > 0) then
			iobs1 = global_obsspace%iobs(i,j,it)
		        iobs2 = iobs1 + global_obsspace%mobs(i,j,it) - 1
			iobs3 = iobs(i,j,it)
			iobs4 = iobs3 + global_obsspace%mobs(i,j,it) - 1
			global_object%field(iobs1:iobs2,:) = field(iobs3:iobs4,:)
		     end if
		  end do
	       end do
	    end do
         end do
	 deallocate(field)
      else
         destination = myide*nxpe*nype
	 call MPI_SEND(self%nobs, 1, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
	 if (self%nobs > 0) call MPI_SEND(self%field, self%nobs*self%nsubobs, r_type, destination, 0, MPI_COMM_WORLD, ierror)
      end if
      !
      if (myidx == 0 .and. myidy == 0) then
         do id = 1, nxpe*nype-1
            destination = myide*nxpe*nype + id
            call MPI_SEND(global_object%field, global_object%nobs*global_object%nsubobs, r_type, destination, 0, MPI_COMM_WORLD, ierror)
         end do
      else
         source = myide*nxpe*nype
         call MPI_RECV(global_object%field, global_object%nobs*global_object%nsubobs, r_type, source, 0, MPI_COMM_WORLD, istat, ierror)
      end if
      deallocate(iobs)
      !
      return
   end subroutine gather_obs_NodeObsField
   !
   !
   !
   subroutine broadcast_ensobs_NodeObsField1(self, source)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: source
      !
      if (self%nobs == 0) return
      call broadcast2D('e', self%field, source, self%nobs, self%nsubobs)
      !
      return
   end subroutine broadcast_ensobs_NodeObsField1
   !
   !
   !
   subroutine broadcast_ensobs_NodeObsField2(self, i, j, obsspace, source)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: source
      integer :: it, iobs1, iobs2
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call broadcast2D('e', self%field(iobs1:iobs2,:), source, obsspace%mobs(i,j,it), self%nsubobs)
      end do
      !
      return
   end subroutine broadcast_ensobs_NodeObsField2
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsField1(self)
      implicit none
      type(NodeObsField), intent(inout) :: self
      !
      if (self%nobs == 0) return
      call allreduce2D('e', self%field, self%nobs, self%nsubobs)
      !
      return
   end subroutine allreduce_ensobs_NodeObsField1
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsField2(self, i, j, obsspace)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer :: it, iobs1, iobs2
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call allreduce2D('e', self%field(iobs1:iobs2,:), obsspace%mobs(i,j,it), self%nsubobs)
      end do
      !
      return
   end subroutine allreduce_ensobs_NodeObsField2
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsField3(self, processed, k2ijt, obsspace, nxyt)
      implicit none
      type(NodeObsField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer :: ixyt, i, j, it, iobs1, iobs2
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call allreduce2D('e', self%field(iobs1:iobs2,:), obsspace%mobs(i,j,it), self%nsubobs)
      end do
      !
      return
   end subroutine allreduce_ensobs_NodeObsField3
   !
   !
   !
end module NodeObsField_class
