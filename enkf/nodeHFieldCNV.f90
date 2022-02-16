module NodeHFieldCNV_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl
   use interpolate, only : PosGrid, convert_LatLon_to_GridPos, convert_GridPos_offset, &
                         & FlagOverTop, interpolation2D, interpolation3D, interpolation3D_T
   use interpolate_TLAD, only : TL_interpolation2D, TL_interpolation3D, TL_interpolation3D_T, &
			      & AD_interpolation2D, AD_interpolation3D, AD_interpolation3D_T
   use NodeInfo_class
   use NodeObsSpaceFieldCNV_class
   use NodeObsValidField_class
   use NodeObsField_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeProfileControl_class
   use NodeMPI
   implicit none
   !
   type NodeHFieldCNV
      character(len=10) :: name
      integer :: nobs
      type(PosGrid), dimension(:), allocatable :: pg
      ! For ensemble linear tangent H
      integer :: ne
      integer, dimension(:), allocatable :: nrank
      real(r_size), dimension(:,:), allocatable :: eigval
      real(r_size), dimension(:,:,:), allocatable :: Y
      real(r_size), dimension(:,:,:,:,:), allocatable :: X, eigvec
   end type NodeHFieldCNV
   !
   interface new
      module procedure new_NodeHFieldCNV
   end interface
   interface destroy
      module procedure destroy_NodeHFieldCNV
   end interface
   interface display
      module procedure display_NodeHFieldCNV
   end interface
   interface get_name
      module procedure get_name_NodeHFieldCNV
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHFieldCNV
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHFieldCNV1
      module procedure get_xyloc_NodeHFieldCNV2
   end interface
   interface apply_Hlogp
      module procedure apply_Hlogp_NodeHFieldCNV
   end interface
   interface apply_H
      module procedure apply_H_NodeHFieldCNV1
      module procedure apply_H_NodeHFieldCNV2
      module procedure apply_H_NodeHFieldCNV3
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHFieldCNV
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHFieldCNV
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHFieldCNV
   end interface
   !
contains
   !
   subroutine new_NodeHFieldCNV(self, info, obsspace)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      if (self%nobs > 0) then
         allocate(self%pg(self%nobs))
         do iobs = 1, self%nobs
            call convert_LatLon_to_GridPos(self%pg(iobs), obsspace%lat(iobs), obsspace%lon(iobs), 1)
            call convert_GridPos_offset(self%pg(iobs), is-1, js-1)
         end do
      end if
      !
      return
   end subroutine new_NodeHFieldCNV
   !
   !
   !
   subroutine destroy_NodeHFieldCNV(self)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      !
      if (allocated(self%pg)) deallocate(self%pg)
      if (allocated(self%nrank)) then
	 deallocate(self%Y, self%X, self%nrank, self%eigval, self%eigvec)
      end if
      !
      return
   end subroutine destroy_NodeHFieldCNV
   !
   !
   !
   subroutine display_NodeHFieldCNV(self)
      implicit none
      type(NodeHFieldCNV), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs
      !
      return
   end subroutine display_NodeHFieldCNV
   !
   !
   !
   subroutine get_name_NodeHFieldCNV(self, name)
      implicit none
      type(NodeHFieldCNV), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeHFieldCNV
   !
   !
   !
   subroutine get_nobs_NodeHFieldCNV(self, nobs)
      implicit none
      type(NodeHFieldCNV), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeHFieldCNV
   !
   !
   !
   subroutine get_xyloc_NodeHFieldCNV1(self, xyloc)
      implicit none
      type(NodeHFieldCNV), intent(in) :: self
      type(NodeObsField), intent(inout) :: xyloc
      integer :: iobs
      !
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px
         xyloc%field(iobs,2) = self%pg(iobs)%py
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldCNV1
   !
   !
   !
   subroutine get_xyloc_NodeHFieldCNV2(self, info, xyloc)
      implicit none
      type(NodeHFieldCNV), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsField), intent(inout) :: xyloc
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px - is + 1
         xyloc%field(iobs,2) = self%pg(iobs)%py - js + 1
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldCNV2
   !
   !
   !
   subroutine interpolate_logp_conv(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: logp, z, g, field
      !
      nx = x%nx
      ny = x%ny
      nz = x%nlev
      nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(x, 'z', 1, z, nz)
      call get_field(x, 'g', 1, g, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'logp', it, logp, nz)
	 call get_field(x, 'logps', it, logp(:,:,1:1), 1)
	 !call update_halo(info, logp, nx, ny, nz)
         field(:,:,1:nz) = logp(:,:,1:nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          select case(obsspace%ztype(iobs))
		  case(IVITP_HEIGHT) ! vertical interpolation with height
		     if (trim(self%name) == 't') then
			call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
		     else
			call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z)
		     end if
		  case(IVITP_LOGP)   ! vertical interpolation with log_p
		     if (trim(self%name) == 't') then
			call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
		     else
			call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp)
		     end if
		  case(IVITP_GRID)   ! vertical interpolation with grid coordinate
		     call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, g)
		  case(IVITP_SURF)   ! surface data: 2D interpolation
		     call interpolation2D(y%field(iobs,1), self%pg(iobs), field(:,:,1))
		  end select
		  if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                    & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
                     valid%field(iobs,1) = 0
                     y%field(iobs,1) = 0.d0
                  else
                     valid%field(iobs,1) = 1
                  end if
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_logp_conv
   !
   !
   !
   subroutine interpolate_logp_rvl(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: z, logp
      !
      nx = x%nx
      ny = x%ny
      nz = x%nlev
      nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(x, 'z', 1, z, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'logp', it, logp, nz)
	 call get_field(x, 'logps', it, logp(:,:,1:1), 1)
	 !call update_halo(info, logp, nx, ny, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), logp, z)
	          if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                    & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  valid%field(iobs,1) = 1
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_logp_rvl
   !
   !
   !
   subroutine interpolate_conv1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: logp, z, g, field
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(x, 'z', 1, z, nz)
      call get_field(x, 'g', 1, g, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'logp', it, logp, nz)
	 call get_field(x, 'logps', it, logp(:,:,1:1), 1)
	 call get_nlev(x, trim(self%name), nlev)
         call get_field(x, trim(self%name), it, field(:,:,1:nlev), nlev)
	 !call update_halo(info, logp, nx, ny, nz)
	 !call update_halo(info, field(:,:,1:nlev), nx, ny, nlev)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          select case(obsspace%ztype(iobs))
		  case(IVITP_HEIGHT) ! vertical interpolation with height
		     if (trim(self%name) == 't') then
			call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
		     else
			call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z)
		     end if
		  case(IVITP_LOGP)   ! vertical interpolation with log_p
		     if (trim(self%name) == 't') then
			call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
		     else
			call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp)
		     end if
		  case(IVITP_GRID)   ! vertical interpolation with grid coordinate
		     call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, g)
		  case(IVITP_SURF)   ! surface data: 2D interpolation
		     call interpolation2D(y%field(iobs,1), self%pg(iobs), field(:,:,1))
		  end select
		  if (obsspace%ztype(iobs) == IVITP_SURF) then
	             valid%field(iobs,1) = 1
	          else if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                         & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
                     valid%field(iobs,1) = 0
                     y%field(iobs,1) = 0.d0
                  else
                     valid%field(iobs,1) = 1
                  end if
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_conv1
   !
   !
   !
   subroutine interpolate_conv2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      logical :: update
      integer :: nz, nt, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, it, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: z, g, logp, field
      !
      nz = x%nlev; nt = obsspace%nt
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(g)) allocate(g(nx,ny,nz))
      if (.not. allocated(logp)) allocate(logp(nx,ny,nz))
      if (.not. allocated(field)) allocate(field(nx,ny,nz))
      !
      call get_nlev(x, trim(self%name), nlev)
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
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'g', ixyt, g(i,j,1:nz), nz)
	    call get_field(x, 'logp', ixyt, logp(i,j,1:nz), nz)
	    call get_field(x, 'logps', ixyt, logp(i,j,1:1), 1)
	    call get_field(x, trim(self%name), ixyt, field(i,j,1:nlev), nlev)
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
	       select case(obsspace%ztype(iobs))
	       case(IVITP_HEIGHT) ! vertical interpolation with height
		  if (trim(self%name) == 't') then
		     call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
		  else
		     call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z)
		  end if
	       case(IVITP_LOGP)   ! vertical interpolation with log_p
		  if (trim(self%name) == 't') then
		     call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
		  else
		     call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp)
		  end if
	       case(IVITP_GRID)   ! vertical interpolation with grid coordinate
		  call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, g)
	       case(IVITP_SURF)   ! surface data: 2D interpolation
		  call interpolation2D(y%field(iobs,1), self%pg(iobs), field(:,:,1))
	       end select
	       if (obsspace%ztype(iobs) == IVITP_SURF) then
	          valid%field(iobs,1) = 1
	       else if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		      & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
		  valid%field(iobs,1) = 0
		  y%field(iobs,1) = 0.d0
	       else
		  valid%field(iobs,1) = 1
	       end if
	    end do
         end do
      end do
      !
      return
   end subroutine interpolate_conv2
   !
   !
   !
   subroutine interpolate_conv3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, it
      real(r_size), dimension(:,:,:), allocatable, save :: z, g, logp, field
      !
      nx = obsspace%nx; ny = obsspace%ny; nz = x%nlev
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(g)) allocate(g(nx,ny,nz))
      if (.not. allocated(logp)) allocate(logp(nx,ny,nz))
      if (.not. allocated(field)) allocate(field(nx,ny,nz))
      !
      call get_nlev(x, trim(self%name), nlev)
      do it = 1, nt
         do i = ip, ip+1
	 do j = jp, jp+1
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'g', ixyt, g(i,j,1:nz), nz)
	    call get_field(x, 'logp', ixyt, logp(i,j,1:nz), nz)
	    call get_field(x, 'logps', ixyt, logp(i,j,1:1), 1)
	    call get_field(x, trim(self%name), ixyt, field(i,j,1:nlev), nlev)
	 end do
	 end do
	 !
	 i = ip; j = jp
	 ixyt = ijt2k(1,1,it)
	 if (ixyt == 0) cycle
	 if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    select case(obsspace%ztype(iobs))
	    case(IVITP_HEIGHT) ! vertical interpolation with height
	       if (trim(self%name) == 't') then
		  call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
	       else
		  call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, z)
	       end if
	    case(IVITP_LOGP)   ! vertical interpolation with log_p
	       if (trim(self%name) == 't') then
		  call interpolation3D_T(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
	       else
		  call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, logp)
	       end if
	    case(IVITP_GRID)   ! vertical interpolation with grid coordinate
	       call interpolation3D(y%field(iobs,1), self%pg(iobs), obsspace%height(iobs), field, g)
	    case(IVITP_SURF)   ! surface data: 2D interpolation
	       call interpolation2D(y%field(iobs,1), self%pg(iobs), field(:,:,1))
	    end select
	    if (obsspace%ztype(iobs) == IVITP_SURF) then
	       valid%field(iobs,1) = 1
	    else if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
	           & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	    else
	       valid%field(iobs,1) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine interpolate_conv3
   !
   !
   !
   subroutine interpolate_Dconv(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, g, logp, field, dz, dg, dlogp, dfield
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(xbck, 'z', 1, z, nz)
      call get_field(xbck, 'g', 1, g, nz)
      call get_field(x, 'z', 1, dz, nz)
      call get_field(x, 'g', 1, dg, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(xbck, 'logp', it, logp, nz)
	 call get_field(x, 'logp', it, dlogp, nz)
	 call get_field(xbck, 'logps', it, logp(:,:,1:1), 1)
	 call get_field(x, 'logps', it, dlogp(:,:,1:1), 1)
	 call get_nlev(xbck, trim(self%name), nlev)
         call get_field(xbck, trim(self%name), it, field(:,:,1:nlev), nlev)
	 call get_field(x, trim(self%name), it, dfield(:,:,1:nlev), nlev)
	 !call update_halo(info, logp, nx, ny, nz)
	 !call update_halo(info, field(:,:,1:nlev), nx, ny, nlev)
	 !call update_halo(info, dlogp, nx, ny, nz)
	 !call update_halo(info, dfield(:,:,1:nlev), nx, ny, nlev)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          select case(obsspace%ztype(iobs))
		  case(IVITP_HEIGHT) ! vertical interpolation with height
		     if (trim(self%name) == 't') then
			call TL_interpolation3D_T(y%field(iobs,1), dfield, dz, self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
		     else
			call TL_interpolation3D(y%field(iobs,1), dfield, dz, self%pg(iobs), obsspace%height(iobs), field, z)
		     end if
		  case(IVITP_LOGP)   ! vertical interpolation with log_p
		     if (trim(self%name) == 't') then
			call TL_interpolation3D_T(y%field(iobs,1), dfield, dlogp, self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
		     else
			call TL_interpolation3D(y%field(iobs,1), dfield, dlogp, self%pg(iobs), obsspace%height(iobs), field, logp)
		     end if
		  case(IVITP_GRID)   ! vertical interpolation with grid coordinate
		     call TL_interpolation3D(y%field(iobs,1), dfield, dg, self%pg(iobs), obsspace%height(iobs), field, g)
		  case(IVITP_SURF)   ! surface data: 2D interpolation
		     call TL_interpolation2D(y%field(iobs,1), dfield(:,:,1), self%pg(iobs))
		  end select
		  if (obsspace%ztype(iobs) == IVITP_SURF) then
	             valid%field(iobs,1) = 1
	          else if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                         & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
                     valid%field(iobs,1) = 0
                     y%field(iobs,1) = 0.d0
                  else
                     valid%field(iobs,1) = 1
                  end if
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_Dconv
   !
   !
   !
   subroutine interpolate_DconvT(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, g, logp, field, dz, dg, dlogp, dfield
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(xbck, 'z', 1, z, nz)
      call get_field(xbck, 'g', 1, g, nz)
      !
      do it = 1, nt
         dz = 0.0d0
	 dg = 0.0d0
         dlogp = 0.0d0
	 dfield = 0.d0
	 call get_nlev(xbck, trim(self%name), nlev)
	 !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
	    call get_field(xbck, 'logp', it, logp, nz)
	    call get_field(xbck, 'logps', it, logp(:,:,1:1), 1)
	    call get_field(xbck, trim(self%name), it, field(:,:,1:nlev), nlev)
	    !call update_halo(info, logp, nx, ny, nz)
	    !call update_halo(info, field(:,:,1:nlev), nx, ny, nlev)
	    !
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     if (obsspace%ztype(iobs) /= IVITP_SURF .and . &
		      & (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		       & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop)) cycle
		     if (valid%field(iobs,1) == 0) cycle
		     select case(obsspace%ztype(iobs))
		     case(IVITP_HEIGHT) ! vertical interpolation with height
			if (trim(self%name) == 't') then
			   call AD_interpolation3D_T(y%field(iobs,1), dfield, dz, self%pg(iobs), obsspace%height(iobs), field, z, 'Z     ')
			else
			   call AD_interpolation3D(y%field(iobs,1), dfield, dz, self%pg(iobs), obsspace%height(iobs), field, z)
			end if
		     case(IVITP_LOGP)   ! vertical interpolation with log_p
			if (trim(self%name) == 't') then
			   call AD_interpolation3D_T(y%field(iobs,1), dfield, dlogp, self%pg(iobs), obsspace%height(iobs), field, logp, 'logP  ')
			else
			   call AD_interpolation3D(y%field(iobs,1), dfield, dlogp, self%pg(iobs), obsspace%height(iobs), field, logp)
			end if
		     case(IVITP_GRID)   ! vertical interpolation with grid coordinate
			call AD_interpolation3D(y%field(iobs,1), dfield, dg, self%pg(iobs), obsspace%height(iobs), field, g)
		     case(IVITP_SURF)   ! surface data: 2D interpolation
			call AD_interpolation2D(y%field(iobs,1), dfield(:,:,1), self%pg(iobs))
		     end select
		  end do
	       end do
	    end do
	 end if
	 !
	 call add_halo(info, dlogp, nx, ny, nz)
	 call add_halo(info, dfield, nx, ny, nlev)
	 call add_field(x, 'logps', it, dlogp(:,:,1:1), 1)
	 dlogp(:,:,1) = 0.d0
	 call add_field(x, 'logp', it, dlogp, nz)
	 call add_field(x, trim(self%name), it, dfield(:,:,1:nlev), nlev)
      end do
      !
      return
   end subroutine interpolate_DconvT
   !
   !
   !
   subroutine smooth_beam_range(pg, height, range, lgpvz, dat, smth, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      type(PosGrid), intent(in) :: pg  ! information of position in the grid
      real(kind = r_sngl), intent(in) ::  height
      real(kind = r_size), intent(in) :: range
      real(kind = r_size), dimension(nx,ny,nz), intent(in) :: lgpvz, dat
      real(kind = r_size), dimension(nx,ny), intent(out) :: smth
      integer :: n, ix, jy, kz, l, lb, le
      real(kind = r_size) :: wgt(nz), wgt_sum, dist
      !
      do n = 1, 4
	 select case(n)
	 case(1)
	    ix = pg%i1; jy = pg%j1; kz = pg%k11
	 case(2)
	    ix = pg%i2; jy = pg%j1; kz = pg%k21
	 case(3)
	    ix = pg%i1; jy = pg%j2; kz = pg%k12
	 case(4)
	    ix = pg%i2; jy = pg%j2; kz = pg%k22
	 end select
	 !
	 lb = max(kz-2, 1)
	 le = min(kz+3, nz-1)
	 wgt_sum = 0.d0
	 wgt(1:nz) = 0.d0
	 do l = lb, le
	    dist = dabs(lgpvz(ix,jy,l)-height)/range
	    wgt(l) = max(exp(-dist*dist),1.0D-6)
	    wgt_sum = wgt_sum + wgt(l)
	 end do
         !
	 if (lgpvz(ix,jy,lb) < height .and. height < lgpvz(ix,jy,le)) then
	    smth(ix,jy) = 0.0d0
	    do l = lb, le
	       smth(ix,jy) = smth(ix,jy) + dat(ix,jy,l)*wgt(l)
	    end do
	    smth(ix,jy) = smth(ix,jy)/wgt_sum
	 end if
      end do
      !
      return
   end subroutine smooth_beam_range
   !
   !
   !
   subroutine AD_smooth_beam_range(pg, height, range, lgpvz, smth, dat, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      type(PosGrid), intent(in) :: pg  ! information of position in the grid
      real(kind = r_sngl), intent(in) :: height
      real(kind = r_size), intent(in) :: range
      real(kind = r_size), dimension(nx,ny,nz), intent(in) :: lgpvz
      real(kind = r_size), dimension(nx,ny), intent(inout) :: smth
      real(kind = r_size), dimension(nx,ny,nz), intent(inout) :: dat
      integer :: n, ix, jy, kz, l, lb, le
      real(kind = r_size) :: wgt(nz), wgt_sum, dist
      !
      do n = 1, 4
	 select case(n)
	 case(1)
	    ix = pg%i1; jy = pg%j1; kz = pg%k11
	 case(2)
	    ix = pg%i2; jy = pg%j1; kz = pg%k21
	 case(3)
	    ix = pg%i1; jy = pg%j2; kz = pg%k12
	 case(4)
	    ix = pg%i2; jy = pg%j2; kz = pg%k22
	 end select
         !
	 lb = max(kz-2, 1)
	 le = min(kz+3, nz-1)
	 wgt_sum = 0.d0
	 wgt(1:nz) = 0.d0
	 do l = lb, le
	    dist = dabs(lgpvz(ix,jy,l)-height)/range
	    wgt(l) = max(exp(-dist*dist),1.0D-6)
	    wgt_sum = wgt_sum + wgt(l)
	 end do
         !
	 if (lgpvz(ix,jy,lb) < height .and. height < lgpvz(ix,jy,le)) then
	    smth(ix,jy) = smth(ix,jy)/wgt_sum
	    do l = lb, le
	       dat(ix,jy,l) = dat(ix,jy,l) + wgt(l)*smth(ix,jy)
	    end do
	    smth(ix,jy) = 0.0d0
	 end if
      end do
      !
      return
   end subroutine AD_smooth_beam_range
   !
   !
   !
   subroutine interpolate_rvl1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size) :: ploc, utmp, vtmp
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: z, logp, u, v
      real(r_size), dimension(x%nx,x%ny) :: urvl, vrvl
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(x, 'z', 1, z, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'logp', it, logp, nz)
	 call get_field(x, 'logps', it, logp(:,:,1:1), 1)
	 call get_field(x, 'u', it, u, nz)
	 call get_field(x, 'v', it, v, nz)
	 !call update_halo(info, logp, nx, ny, nz)
	 !call update_halo(info, u, nx, ny, nz)
	 !call update_halo(info, v, nx, ny, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          call interpolation3D(ploc, self%pg(iobs), obsspace%height(iobs), logp, z)
	          if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                    & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, u, urvl, nx, ny, nz)
		  call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, v, vrvl, nx, ny, nz)
		  call interpolation2D(utmp, self%pg(iobs), urvl)
		  call interpolation2D(vtmp, self%pg(iobs), vrvl)
		  valid%field(iobs,1) = 1
		  y%field(iobs,1) = utmp*sin(obsspace%azimuth(iobs)) + vtmp*cos(obsspace%azimuth(iobs))
		end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_rvl1
   !
   !
   !
   subroutine interpolate_rvl2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      logical :: update
      integer :: nz, nt
      integer :: ixyt, iobs, iobs1, iobs2, i, j, it, imin, jmin, ii, jj
      real(r_size) :: ploc, utmp, vtmp
      real(r_size), dimension(:,:), allocatable, save :: urvl, vrvl
      real(r_size), dimension(:,:,:), allocatable, save :: z, logp, u, v
      !
      nz = x%nlev; nt = obsspace%nt
      if (.not. allocated(urvl)) allocate(urvl(nx,ny))
      if (.not. allocated(vrvl)) allocate(vrvl(nx,ny))
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(logp)) allocate(logp(nx,ny,nz))
      if (.not. allocated(u)) allocate(u(nx,ny,nz))
      if (.not. allocated(v)) allocate(v(nx,ny,nz))
      !
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
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'logp', ixyt, logp(i,j,1:nz), nz)
	    call get_field(x, 'logps', ixyt, logp(i,j,1:1), 1)
	    call get_field(x, 'u', ixyt, u(i,j,1:nz), nz)
	    call get_field(x, 'v', ixyt, v(i,j,1:nz), nz)
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
	       call interpolation3D(ploc, self%pg(iobs), obsspace%height(iobs), logp, z)
	       if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		 & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
		  valid%field(iobs,1) = 0
		  y%field(iobs,1) = 0.d0
		  cycle
	       end if
	       if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
		 & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
		 & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
		 & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) then
		  valid%field(iobs,1) = 0
		  y%field(iobs,1) = 0.d0
		  cycle
	       end if
	       call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, u, urvl, nx, ny, nz)
	       call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, v, vrvl, nx, ny, nz)
	       call interpolation2D(utmp, self%pg(iobs), urvl)
	       call interpolation2D(vtmp, self%pg(iobs), vrvl)
	       valid%field(iobs,1) = 1
	       y%field(iobs,1) = utmp*sin(obsspace%azimuth(iobs)) + vtmp*cos(obsspace%azimuth(iobs))
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_rvl2
   !
   !
   !
   subroutine interpolate_rvl3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz
      integer :: ixyt, iobs, iobs1, iobs2, i, j, it
      real(r_size) :: ploc, utmp, vtmp
      real(r_size), dimension(:,:), allocatable, save :: urvl, vrvl
      real(r_size), dimension(:,:,:), allocatable, save :: z, logp, u, v
      !
      nx = obsspace%nx; ny = obsspace%ny; nz = x%nlev
      if (.not. allocated(urvl)) allocate(urvl(nx,ny))
      if (.not. allocated(vrvl)) allocate(vrvl(nx,ny))
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(logp)) allocate(logp(nx,ny,nz))
      if (.not. allocated(u)) allocate(u(nx,ny,nz))
      if (.not. allocated(v)) allocate(v(nx,ny,nz))
      !
      do it = 1, nt
         do i = ip, ip+1
	 do j = jp, jp+1
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'logp', ixyt, logp(i,j,1:nz), nz)
	    call get_field(x, 'logps', ixyt, logp(i,j,1:1), 1)
	    call get_field(x, 'u', ixyt, u(i,j,1:nz), nz)
	    call get_field(x, 'v', ixyt, v(i,j,1:nz), nz)
	 end do
	 end do
         !
	 i = ip; j = jp
	 ixyt = ijt2k(1,1,it)
	 if (ixyt == 0) cycle
	 if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    call interpolation3D(ploc, self%pg(iobs), obsspace%height(iobs), logp, z)
	    if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
	      & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	       cycle
	    end if
	    if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
	      & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
	      & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
	      & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	       cycle
	    end if
	    call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, u, urvl, nx, ny, nz)
	    call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, v, vrvl, nx, ny, nz)
	    call interpolation2D(utmp, self%pg(iobs), urvl)
	    call interpolation2D(vtmp, self%pg(iobs), vrvl)
	    valid%field(iobs,1) = 1
	    y%field(iobs,1) = utmp*sin(obsspace%azimuth(iobs)) + vtmp*cos(obsspace%azimuth(iobs))
	 end do
      end do
      !
      return
   end subroutine interpolate_rvl3
   !
   !
   !
   subroutine interpolate_Drvl(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size) :: ploc, dutmp, dvtmp
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, du, dv
      real(r_size), dimension(xbck%nx,xbck%ny) :: durvl, dvrvl
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(xbck, 'z', 1, z, nz)
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'u', it, du, nz)
	 call get_field(x, 'v', it, dv, nz)
	 !call update_halo(info, du, nx, ny, nz)
	 !call update_halo(info, dv, nx, ny, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
                    & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
		    & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, du, durvl, nx, ny, nz)
		  call smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, dv, dvrvl, nx, ny, nz)
		  call TL_interpolation2D(dutmp, durvl, self%pg(iobs))
		  call TL_interpolation2D(dvtmp, dvrvl, self%pg(iobs))
		  valid%field(iobs,1) = 1
		  y%field(iobs,1) = dutmp*sin(obsspace%azimuth(iobs)) + dvtmp*cos(obsspace%azimuth(iobs))
		end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_Drvl
   !
   !
   !
   subroutine interpolate_DrvlT(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: nx, ny, nz, nt
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size) :: ploc, dutmp, dvtmp
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, du, dv
      real(r_size), dimension(xbck%nx,xbck%ny) :: durvl, dvrvl
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      call get_field(xbck, 'z', 1, z, nz)
      durvl = 0.d0
      dvrvl = 0.d0
      !
      do it = 1, nt
         du = 0.d0
         dv = 0.d0
         !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		       & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) cycle
		     if (valid%field(iobs,1) == 0) cycle
		     if (obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j1,4) .or. &
		       & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j1,4) .or. &
		       & obsspace%height(iobs) < z(self%pg(iobs)%i1,self%pg(iobs)%j2,4) .or. &
		       & obsspace%height(iobs) < z(self%pg(iobs)%i2,self%pg(iobs)%j2,4)) cycle
		     dutmp = y%field(iobs,1)*sin(obsspace%azimuth(iobs))
		     dvtmp = y%field(iobs,1)*cos(obsspace%azimuth(iobs))
		     call AD_interpolation2D(dutmp, durvl, self%pg(iobs))
		     call AD_interpolation2D(dvtmp, dvrvl, self%pg(iobs))
		     call AD_smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, durvl, du, nx, ny, nz)
		     call AD_smooth_beam_range(self%pg(iobs), obsspace%height(iobs), obsspace%range(iobs), z, dvrvl, dv, nx, ny, nz)
		   end do
	       end do
	    end do
	 end if
	 !
	 call add_halo(info, du, nx, ny, nz)
	 call add_halo(info, dv, nx, ny, nz)
	 call add_field(x, 'u', it, du, nz)
	 call add_field(x, 'v', it, dv, nz)
      end do
      !
      return
   end subroutine interpolate_DrvlT
   !
   !
   !
   subroutine apply_Hlogp_NodeHFieldCNV(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_logp_conv(self, info, x, obsspace, valid, y)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_logp_rvl(self, info, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine apply_Hlogp_NodeHFieldCNV
   !
   !
   !
   subroutine apply_H_NodeHFieldCNV1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_conv1(self, info, x, obsspace, valid, y)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_rvl1(self, info, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldCNV1
   !
   !
   !
   subroutine apply_H_NodeHFieldCNV2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_conv2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_rvl2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldCNV2
   !
   !
   !
   subroutine apply_H_NodeHFieldCNV3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_conv3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_rvl3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldCNV3
   
   !
   !
   !
   subroutine initialize_DH_NodeHFieldCNV(self, info, x, obsspace, valid, y, ne)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 !call initialize_Dconv(self, info, x, obsspace, valid, y)
      else if (trim(self%name) == 'rvl') then
	 !call initialize_Drvl(self, info, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine initialize_DH_NodeHFieldCNV
   !
   !
   !
   subroutine apply_DH_NodeHFieldCNV(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_Dconv(self, info, xbck, x, obsspace, valid, y)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_Drvl(self, info, xbck, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine apply_DH_NodeHFieldCNV
   !
   !
   !
   subroutine apply_DHT_NodeHFieldCNV(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldCNV), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldCNV), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      if (trim(self%name) == 'u' .or. trim(self%name) == 'v' .or. trim(self%name) == 't' .or. &
	& trim(self%name) == 'p' .or. trim(self%name) == 'rh' .or. trim(self%name) == 'pwv') then
	 call interpolate_DconvT(self, info, xbck, obsspace, valid, y, x)
      else if (trim(self%name) == 'rvl') then
	 call interpolate_DrvlT(self, info, xbck, obsspace, valid, y, x)
      end if
      !
      return
   end subroutine apply_DHT_NodeHFieldCNV
   !
   !
   !
end module NodeHFieldCNV_class
