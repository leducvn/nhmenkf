module NodeObsVLocField_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsSpaceField_class
   use NodeObsField_class
   use NodeProfileControl_class
   use NodeControl_class
   use NodeMPI
   implicit none
   !
   type NodeObsVLocField
      character(len=10) :: name, varname(2)
      integer :: nobs, nsubobs, nz, nvar
      real(r_size), dimension(:,:,:), allocatable :: field
   end type NodeObsVLocField
   !
   interface new
      module procedure new_NodeObsVLocField1
      module procedure new_NodeObsVLocField2
   end interface
   interface destroy
      module procedure destroy_NodeObsVLocField
   end interface
   interface display
      module procedure display_NodeObsVLocField
   end interface
   interface assignment(=)
      module procedure copy_NodeObsVLocField
   end interface
   interface get_name
      module procedure get_name_NodeObsVLocField
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsVLocField
   end interface
   interface get_field
      module procedure get_field_NodeObsVLocField
   end interface
   interface set_field
      module procedure set_field_NodeObsVLocField1
      module procedure set_field_NodeObsVLocField2
   end interface
   interface add_vloc
      module procedure add_vloc_NodeObsVLocField
   end interface
   interface subtract_vloc
      module procedure subtract_vloc_NodeObsVLocField
   end interface
   interface power_vloc
      module procedure power_vloc_NodeObsVLocField
   end interface
   interface multiply_vloc
      module procedure multiply_vloc_NodeObsVLocField1
      module procedure multiply_vloc_NodeObsVLocField2
      module procedure multiply_vloc_NodeObsVLocField3
   end interface
   interface divide_vloc
      module procedure divide_vloc_NodeObsVLocField1
      module procedure divide_vloc_NodeObsVLocField2
   end interface
   interface set_logp
      module procedure set_logp_NodeObsVLocField1
      module procedure set_logp_NodeObsVLocField2
   end interface
   interface set_profile
      module procedure set_profile_NodeObsVLocField1
      module procedure set_profile_NodeObsVLocField2
   end interface
   interface ecorap
      module procedure ecorap_NodeObsVLocField
   end interface
   interface allreduce_ensvloc
      module procedure allreduce_ensvloc_NodeObsVLocField
   end interface
   !
contains
   !
   !
   !
   subroutine new_NodeObsVLocField1(self, obsspace, nz)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: nz
      !
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      self%nsubobs = obsspace%nsubobs
      self%nz = nz
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. &
	& trim(self%name) == 't' .or. trim(self%name) == 'rh') then
	 self%nvar = 1
	 self%varname(1) = self%name
      else if (trim(self%name) == 'p') then
	 self%nvar = 1
	 self%varname(1) = 'p'
      else if (trim(self%name) == 'pwv') then
	 self%nvar = 1
	 self%varname(1) = 'qv'
      else if (trim(self%name) == 'rvl') then
	 self%nvar = 2
	 self%varname(1) = 'u'
	 self%varname(2) = 'v'
      else if (trim(self%name) == 'tc') then
	 self%nvar = 2
	 self%varname(1) = 't'
	 self%varname(2) = 'qv'
      else if (trim(self%name) == 're' .or. trim(self%name) == 'ba') then
	 self%nvar = 2
	 self%varname(1) = 't'
	 self%varname(2) = 'qv'
      else if (trim(self%name) == 'rad') then
	 self%nvar = 2
	 self%varname(1) = 't'
	 self%varname(2) = 'qv'
      end if
      !
      if (self%nobs > 0) then
	 allocate(self%field(self%nobs,self%nsubobs,self%nz*self%nvar))
         self%field = 0.d0
      end if
      !
      return
   end subroutine new_NodeObsVLocField1
   !
   !
   !
   subroutine new_NodeObsVLocField2(self, obsspace, nz, nvar)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: nz, nvar
      !
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      self%nsubobs = obsspace%nsubobs
      self%nz = nz
      self%nvar = nvar
      if (self%nobs > 0) then
	 allocate(self%field(self%nobs,self%nsubobs,self%nz*self%nvar))
         self%field = 0.d0
      end if
      !
      return
   end subroutine new_NodeObsVLocField2
   !
   !
   !
   subroutine destroy_NodeObsVLocField(self)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      !
      if (allocated(self%field)) deallocate(self%field)
      !
      return
   end subroutine destroy_NodeObsVLocField
   !
   !
   !
   subroutine display_NodeObsVLocField(self)
      implicit none
      type(NodeObsVLocField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs, self%nsubobs, self%nz, self%nvar
      !
      return
   end subroutine display_NodeObsVLocField
   !
   !
   !
   subroutine copy_NodeObsVLocField(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = object%field(:,:,:)
      !
      return
   end subroutine copy_NodeObsVLocField
   !
   !
   !
   subroutine get_name_NodeObsVLocField(self, name)
      implicit none
      type(NodeObsVLocField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsVLocField
   !
   !
   !
   subroutine get_nobs_NodeObsVLocField(self, nobs)
      implicit none
      type(NodeObsVLocField), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsVLocField
   !
   !
   !
   subroutine get_field_NodeObsVLocField(self, k, object)
      implicit none
      type(NodeObsVLocField), intent(in) :: self
      integer, intent(in) :: k
      type(NodeObsField), intent(inout) :: object
      !
      if (self%nobs == 0) return
      object%field(:,:) = self%field(:,:,k)
      !
      return
   end subroutine get_field_NodeObsVLocField
   !
   !
   !
   subroutine set_field_NodeObsVLocField1(self, const)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = const
      !
      return
   end subroutine set_field_NodeObsVLocField1
   !
   !
   !
   subroutine set_field_NodeObsVLocField2(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = object%field(:,:,:)
      !
      return
   end subroutine set_field_NodeObsVLocField2
   
   !
   !
   !
   subroutine add_vloc_NodeObsVLocField(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = self%field(:,:,:) + object%field(:,:,:)
      !
      return
   end subroutine add_vloc_NodeObsVLocField
   !
   !
   !
   subroutine subtract_vloc_NodeObsVLocField(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = self%field(:,:,:) - object%field(:,:,:)
      !
      return
   end subroutine subtract_vloc_NodeObsVLocField
   !
   !
   !
   subroutine power_vloc_NodeObsVLocField(self, const)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = self%field(:,:,:)**const
      !
      return
   end subroutine power_vloc_NodeObsVLocField
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocField1(self, const)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = const*self%field(:,:,:)
      !
      return
   end subroutine multiply_vloc_NodeObsVLocField1
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocField2(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = object%field(:,:,:)*self%field(:,:,:)
      !
      return
   end subroutine multiply_vloc_NodeObsVLocField2
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocField3(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsField), intent(in) :: object
      integer :: iobs, j
      !
      do iobs = 1, self%nobs
         do j = 1, self%nsubobs
            self%field(iobs,j,:) = object%field(iobs,j)*self%field(iobs,j,:)
         end do
      end do
      !
      return
   end subroutine multiply_vloc_NodeObsVLocField3
   !
   !
   !
   subroutine divide_vloc_NodeObsVLocField1(self, const)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:,:) = self%field(:,:,:)/const
      !
      return
   end subroutine divide_vloc_NodeObsVLocField1
   !
   !
   !
   subroutine divide_vloc_NodeObsVLocField2(self, object)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeObsVLocField), intent(in) :: object
      !
      if (self%nobs == 0) return
      where (object%field < 1.e-12)
         self%field = 0.d0
      else where
         self%field = self%field/object%field
      end where
      !
      return
   end subroutine divide_vloc_NodeObsVLocField2
   !
   !
   !
   subroutine set_logp_NodeObsVLocField1(self, info, obsspace, x)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeControl), intent(in) :: x
      !
      integer :: nx, ny, nz, nt, halo
      integer :: i, j, k, it, iobs, iobs1, iobs2
      integer :: dis, die, djs, dje, is, ie, js, je, istart, iend, jstart, jend
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: field
      !
      nx = x%nx
      ny = x%ny
      nz = x%nlev
      nt = x%nt
      call get_halo(info, halo)
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'logp', it, field, nz)
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
		  istart = i - halo + 1
		  iend = i + halo
		  jstart = j - halo + 1
		  jend = j + halo
		  do k = 1, nz
		     self%field(iobs,:,k) = sum(field(istart:iend,jstart:jend,k))/(4*halo**2)
		  end do
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine set_logp_NodeObsVLocField1
   !
   !
   !
   subroutine set_logp_NodeObsVLocField2(self, info, processed, k2ijt, obsspace, x, nx, ny, nxyt)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeProfileControl), intent(in) :: x
      !
      logical :: update
      integer :: nz, nt
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: field
      !
      nz = x%nlev; nt = obsspace%nt
      if (.not. allocated(field)) allocate(field(nx,ny,nz))
      do it = 1, nt
	 do ixyt = 1, nxyt
	    if (it /= k2ijt(ixyt,3)) cycle
	    i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	    update = .False.
	    imin = max(1,i-1); jmin = max(1,j-1)
	    do ii = imin, i
	    do jj = jmin, j
	       if (.not. processed(ii,jj)) update = .True.
	    end do
	    end do
	    if (.not. update) cycle
	    call get_field(x, 'logp', ixyt, field(i,j,1:nz), nz)
	 end do
	 !
	 do ixyt = 1, nxyt
	    if (it /= k2ijt(ixyt,3)) cycle
	    i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	    if (processed(i,j)) cycle
	    if (obsspace%mobs(i,j,it) == 0) cycle
	    !
	    iobs1 = obsspace%iobs(i,j,it)
	    iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	    do iobs = iobs1, iobs2
	       do k = 1, nz
		  self%field(iobs,:,k) = sum(field(i:i+1,j:j+1,k))/4
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine set_logp_NodeObsVLocField2
   !
   !
   !
   subroutine set_profile_NodeObsVLocField1(self, info, obsspace, x)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeControl), intent(in) :: x
      !
      integer :: nx, ny, nz, nt, nvar, halo
      integer :: i, j, k, it, ivar, iobs, iobs1, iobs2
      integer :: dis, die, djs, dje, is, ie, js, je, istart, iend, jstart, jend
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: field
      !
      nx = x%nx
      ny = x%ny
      nz = x%nlev
      nt = x%nt
      nvar = self%nvar
      call get_halo(info, halo)
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 do ivar = 1, nvar
	    call get_field(x, self%varname(ivar), it, field, nz)
	    !call update_halo(info, field, nx, ny, nz)
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     istart = i - halo + 1
		     iend = i + halo
		     jstart = j - halo + 1
		     jend = j + halo
		     do k = 1, nz
		        self%field(iobs,:,(ivar-1)*nz+k) = sum(field(istart:iend,jstart:jend,k))/(4*halo**2)
		     end do
                  end do
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine set_profile_NodeObsVLocField1
   !
   !
   !
   subroutine set_profile_NodeObsVLocField2(self, info, processed, k2ijt, obsspace, x, nx, ny, nxyt)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeProfileControl), intent(in) :: x
      !
      logical :: update
      integer :: nz, nt, nvar
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it, ivar, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: field
      !
      nz = x%nlev; nt = obsspace%nt; nvar = self%nvar
      if (.not. allocated(field)) allocate(field(nx,ny,nz))
      do it = 1, nt
         do ivar = 1, nvar
  	    do ixyt = 1, nxyt
	       if (it /= k2ijt(ixyt,3)) cycle
	       i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	       update = .False.
	       imin = max(1,i-1); jmin = max(1,j-1)
	       do ii = imin, i
	       do jj = jmin, j
	          if (.not. processed(ii,jj)) update = .True.
	       end do
	       end do
	       if (.not. update) cycle
	       call get_field(x, self%varname(ivar), ixyt, field(i,j,1:nz), nz)
	    end do
	    !
	    do ixyt = 1, nxyt
	       if (it /= k2ijt(ixyt,3)) cycle
	       i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	       if (processed(i,j)) cycle
	       if (obsspace%mobs(i,j,it) == 0) cycle
	       !
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          do k = 1, nz
		     self%field(iobs,:,(ivar-1)*nz+k) = sum(field(i:i+1,j:j+1,k))/4
	          end do
	       end do
	    end do
         end do
      end do
      !
      return
   end subroutine set_profile_NodeObsVLocField2
   !
   !
   !
   subroutine ecorap_NodeObsVLocField(self, scale, logpfield, vlocfield)
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeObsVLocField), intent(in) :: self
      real(r_size), intent(in) :: scale
      type(NodeObsVLocField), intent(in) :: logpfield
      type(NodeObsVLocField), intent(inout) :: vlocfield
      !
      integer :: nz, nvar
      integer :: iobs, j, ivar, k, k1, k2, kk
      real(r_size) :: vlocmax, d, vweight
      real(r_size), dimension(self%nz) :: corr, logp, vloc
      !
      nz = self%nz
      nvar = self%nvar
      do iobs = 1, self%nobs
         do j = 1, self%nsubobs
	    logp(:) = logpfield%field(iobs,j,:)
	    corr(:) = 0.d0
	    vloc(:) = 0.d0
	    do ivar = 1, nvar
	       k1 = (ivar-1)*nz + 1
	       k2 = k1 + nz - 1
               corr(:) = corr(:) + abs(self%field(iobs,j,k1:k2))
	    end do
	    corr(:) = corr(:)/nvar
	    corr(:) = corr(:)**6
	    !
	    do k = 1, nz
	       do kk = 1, nz
	          !vloc(k) = vloc(k) + exp(-0.5*(logp(kk)-logp(k))/scale**2)*corr(kk)
		  d = logp(kk)-logp(k)
		  call GaspariCohn(d, scale, vweight)
		  vloc(k) = vloc(k) + vweight*corr(kk)
	       end do
	    end do
	    vlocmax = maxval(vloc(:))
	    if (vlocmax < 1.e-12) then
	       vloc(:) = 0.d0
	    else
	       vloc(:) = vloc(:)/vlocmax
	    end if
	    vlocfield%field(iobs,j,:) = vloc(:)
         end do
      end do
      !
      return
   end subroutine ecorap_NodeObsVLocField
   !
   !
   !
   subroutine allreduce_ensvloc_NodeObsVLocField(self)
      implicit none
      type(NodeObsVLocField), intent(inout) :: self
      !
      if (self%nobs == 0) return
      call allreduce3D('e', self%field, self%nobs, self%nsubobs, self%nz*self%nvar)
      !
      return
   end subroutine allreduce_ensvloc_NodeObsVLocField
   !
   !
   !
end module NodeObsVLocField_class
