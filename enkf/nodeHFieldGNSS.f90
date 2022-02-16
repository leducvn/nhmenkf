module NodeHFieldGNSS_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl
   use interpolate, only : PosGrid, convert_LatLon_to_GridPos, convert_GridPos_offset, &
                         & FlagOverTop, interpolation2D, interpolation3D
   use interpolate_TLAD, only : TL_interpolation2D, AD_interpolation2D
   ! --- ROPP modules
   use ropp_fm, only : ropp_fm_refrac_1d, ropp_fm_refrac_1d_tl, ropp_fm_refrac_1d_ad, &
                     & ropp_fm_bangle_1d, ropp_fm_bangle_1d_tl, ropp_fm_bangle_1d_ad
   use ropp_fm_types, only : State1dFM, Obs1drefrac, Obs1dbangle
   use geodesy, only : geometric2geopotential, gravity, r_eff
   use NodeInfo_class
   use NodeObsSpaceFieldGNSS_class
   use NodeObsValidField_class
   use NodeObsField_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeProfileControl_class
   use NodeMPI
   implicit none
   !
   type NodeHFieldGNSS
      character(len=10) :: name
      integer :: nobs
      type(PosGrid), dimension(:), allocatable :: pg
      ! For ensemble linear tangent H
      integer :: ne
      integer, dimension(:), allocatable :: nrank
      real(r_size), dimension(:,:), allocatable :: eigval
      real(r_size), dimension(:,:,:), allocatable :: Y
      real(r_size), dimension(:,:,:,:,:), allocatable :: X, eigvec
   end type NodeHFieldGNSS
   !
   interface new
      module procedure new_NodeHFieldGNSS
   end interface
   interface destroy
      module procedure destroy_NodeHFieldGNSS
   end interface
   interface display
      module procedure display_NodeHFieldGNSS
   end interface
   interface get_name
      module procedure get_name_NodeHFieldGNSS
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHFieldGNSS
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHFieldGNSS1
      module procedure get_xyloc_NodeHFieldGNSS2
   end interface
   interface apply_H
      module procedure apply_H_NodeHFieldGNSS1
      module procedure apply_H_NodeHFieldGNSS2
      module procedure apply_H_NodeHFieldGNSS3
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHFieldGNSS
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHFieldGNSS
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHFieldGNSS
   end interface
   !
contains
   !
   subroutine new_NodeHFieldGNSS(self, info, obsspace)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      if (self%nobs > 0) then
         allocate(self%pg(self%nobs))
         do iobs = 1, self%nobs
            call convert_LatLon_to_GridPos(self%pg(iobs), real(obsspace%lat(iobs),r_sngl), real(obsspace%lon(iobs),r_sngl), 1)
	    call convert_GridPos_offset(self%pg(iobs), is-1, js-1)
         end do
      end if
      !
      return
   end subroutine new_NodeHFieldGNSS
   !
   !
   !
   subroutine destroy_NodeHFieldGNSS(self)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      !
      if (allocated(self%pg)) deallocate(self%pg)
      if (allocated(self%nrank)) then
	 deallocate(self%Y, self%X, self%nrank, self%eigval, self%eigvec)
      end if
      !
      return
   end subroutine destroy_NodeHFieldGNSS
   !
   !
   !
   subroutine display_NodeHFieldGNSS(self)
      implicit none
      type(NodeHFieldGNSS), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs
      !
      return
   end subroutine display_NodeHFieldGNSS
   !
   !
   !
   subroutine get_name_NodeHFieldGNSS(self, name)
      implicit none
      type(NodeHFieldGNSS), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeHFieldGNSS
   !
   !
   !
   subroutine get_nobs_NodeHFieldGNSS(self, nobs)
      implicit none
      type(NodeHFieldGNSS), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeHFieldGNSS
   !
   !
   !
   subroutine get_xyloc_NodeHFieldGNSS1(self, xyloc)
      implicit none
      type(NodeHFieldGNSS), intent(in) :: self
      type(NodeObsField), intent(inout) :: xyloc
      integer :: iobs
      !
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px
	 xyloc%field(iobs,2) = self%pg(iobs)%py
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldGNSS1
   !
   !
   !
   subroutine get_xyloc_NodeHFieldGNSS2(self, info, xyloc)
      implicit none
      type(NodeHFieldGNSS), intent(in) :: self
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
   end subroutine get_xyloc_NodeHFieldGNSS2
   !
   !
   !
   subroutine interpolate_re1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1drefrac) :: y1D
      ! --- ROPP variables
      real(r_size) :: ptmp
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%refrac(1), y1D%geop(1), y1D%weights(1))
      !
      call get_field(x, 'z', 1, z, nz)
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'p', it, p, nz)
	 call get_field(x, 't', it, t, nz)
         call get_field(x, 'qv', it, qv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          y1D%lon = obsspace%lon(iobs)
                  y1D%lat = obsspace%lat(iobs)
                  y1D%refrac(1) = obsspace%obs(iobs,1)
                  y1D%geop(1) = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
                  y1D%weights(1) = 1.0d0
                  do k = 2, nz - 1
                     call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                     call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                     call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                     call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                     if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
                  end do
                  x1D%non_ideal = .false.
                  x1D%lon = y1D%lon
                  x1D%lat = y1D%lat
                  call ropp_fm_refrac_1d(x1D, y1D)
                  y%field(iobs,1) = y1D%refrac(1)
		  !
		  if (y1D%refrac(1) < 0.0d0 ) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  call interpolation3D(ptmp, self%pg(iobs), real(y1D%geop(1),r_sngl), p, z)
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%refrac, y1D%geop, y1D%weights)
      !
      return
   end subroutine interpolate_re1
   !
   !
   !
   subroutine interpolate_re2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      logical :: update
      integer :: nz, nt, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1drefrac) :: y1D
      ! --- ROPP variables
      real(r_size) :: ptmp
      !
      nz = x%nlev; nt = obsspace%nt
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      x1D%n_lev = nz-2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%refrac(1), y1D%geop(1), y1D%weights(1))
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
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
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
	       y1D%lon = obsspace%lon(iobs)
               y1D%lat = obsspace%lat(iobs)
               y1D%refrac(1) = obsspace%obs(iobs,1)
               y1D%geop(1) = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
               y1D%weights(1) = 1.0d0
               do k = 2, nz-1
                  call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                  call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                  call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                  call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                  if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
               end do
               x1D%non_ideal = .false.
               x1D%lon = y1D%lon
               x1D%lat = y1D%lat
               call ropp_fm_refrac_1d(x1D, y1D)
               y%field(iobs,1) = y1D%refrac(1)
	       !
	       if (y1D%refrac(1) < 0.0d0 ) then
		  valid%field(iobs,1) = 0
		  y%field(iobs,1) = 0.d0
		  cycle
	       end if
	       call interpolation3D(ptmp, self%pg(iobs), real(y1D%geop(1),r_sngl), p, z)
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%refrac, y1D%geop, y1D%weights)
      !
      return
   end subroutine interpolate_re2
   !
   !
   !
   subroutine interpolate_re3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it
      real(r_size), dimension(:,:,:), allocatable, save :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1drefrac) :: y1D
      ! --- ROPP variables
      real(r_size) :: ptmp
      !
      nx = obsspace%nx; ny = obsspace%ny; nz = x%nlev
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      x1D%n_lev = nz-2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%refrac(1), y1D%geop(1), y1D%weights(1))
      !
      do it = 1, nt
         do i = ip, ip+1
	 do j = jp, jp+1
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
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
	    y1D%lon = obsspace%lon(iobs)
            y1D%lat = obsspace%lat(iobs)
            y1D%refrac(1) = obsspace%obs(iobs,1)
            y1D%geop(1) = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
            y1D%weights(1) = 1.0d0
            do k = 2, nz-1
               call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
               call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
               call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
               call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
               if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
            end do
            x1D%non_ideal = .false.
            x1D%lon = y1D%lon
            x1D%lat = y1D%lat
            call ropp_fm_refrac_1d(x1D, y1D)
            y%field(iobs,1) = y1D%refrac(1)
	    !
	    if (y1D%refrac(1) < 0.0d0 ) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	       cycle
	    end if
	    call interpolation3D(ptmp, self%pg(iobs), real(y1D%geop(1),r_sngl), p, z)
	    if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
	      & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	    else
	       valid%field(iobs,1) = 1
	    end if
	 end do
      end do
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%refrac, y1D%geop, y1D%weights)
      !
      return
   end subroutine interpolate_re3
   !
   !
   !
   subroutine interpolate_Dre(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, p, t, qv, dz, dp, dt, dqv
      ! --- ROPP types
      type(State1dFM) :: x1D, dx1D
      type(Obs1drefrac) :: y1D
      ! --- ROPP variables
      real(r_size) :: dy1D(1)
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      dx1D%n_lev = nz - 2
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(dx1D%temp(dx1D%n_lev), dx1D%shum(dx1D%n_lev), dx1D%pres(dx1D%n_lev), dx1D%geop(dx1D%n_lev))
      allocate(y1D%refrac(1), y1D%geop(1), y1D%weights(1))
      !
      call get_field(xbck, 'z', 1, z, nz)
      !call get_field(x, 'z', 1, dz, nz)
      dz(:,:,:) = 0.d0
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(xbck, 'p', it, p, nz)
	 call get_field(xbck, 't', it, t, nz)
         call get_field(xbck, 'qv', it, qv, nz)
	 call get_field(x, 'p', it, dp, nz)
	 call get_field(x, 't', it, dt, nz)
         call get_field(x, 'qv', it, dqv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          y1D%lon = obsspace%lon(iobs)
                  y1D%lat = obsspace%lat(iobs)
                  y1D%refrac(1) = obsspace%obs(iobs,1)
                  y1D%geop(1) = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
                  y1D%weights(1) = 1.0d0
                  do k = 2, nz - 1
                     call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                     call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                     call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                     call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                     if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
		     call TL_interpolation2D(dx1D%temp(k-1), dt(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%geop(k-1), dz(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%pres(k-1), dp(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%shum(k-1), dqv(:,:,k), self%pg(iobs))
                  end do
                  x1D%non_ideal = .false.
                  x1D%lon = y1D%lon
                  x1D%lat = y1D%lat
		  call ropp_fm_refrac_1d_tl(x1D, dx1D, y1D, dy1D)
                  y%field(iobs,1) = dy1D(1)
	          !
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(dx1D%temp, dx1D%shum, dx1D%pres, dx1D%geop)
      deallocate(y1D%refrac, y1D%geop, y1D%weights)
      !
      return
   end subroutine interpolate_Dre
   !
   !
   !
   subroutine interpolate_DreT(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, p, t, qv, dz, dp, dt, dqv
      ! --- ROPP types
      type(State1dFM) :: x1D, dx1D
      type(Obs1drefrac) :: y1D
      ! --- ROPP variables
      real(r_size) :: dy1D(1)
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      dx1D%n_lev = nz - 2
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(dx1D%temp(dx1D%n_lev), dx1D%shum(dx1D%n_lev), dx1D%pres(dx1D%n_lev), dx1D%geop(dx1D%n_lev))
      allocate(y1D%refrac(1), y1D%geop(1), y1D%weights(1))
      !
      call get_field(xbck, 'z', 1, z, nz)
      do it = 1, nt
         dz = 0.0d0
	 dp = 0.0d0
         dt = 0.0d0
	 dqv = 0.d0
	 !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
	    call get_field(xbck, 'p', it, p, nz)
	    call get_field(xbck, 't', it, t, nz)
            call get_field(xbck, 'qv', it, qv, nz)
	    !
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		       & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) cycle
		     if (valid%field(iobs,1) == 0) cycle
		     !
		     y1D%lon = obsspace%lon(iobs)
                     y1D%lat = obsspace%lat(iobs)
                     y1D%refrac(1) = obsspace%obs(iobs,1)
                     y1D%geop(1) = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
                     y1D%weights(1) = 1.0d0
                     do k = 2, nz - 1
                        call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                        call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                        call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                        call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                        if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
                     end do
                     x1D%non_ideal = .false.
                     x1D%lon = y1D%lon
                     x1D%lat = y1D%lat
                     dy1D(1) = y%field(iobs,1)
		     dx1D%temp(:) = 0.0d0
                     dx1D%shum(:) = 0.0d0
                     dx1D%pres(:) = 0.0d0
                     dx1D%geop(:) = 0.0d0 
		     call ropp_fm_refrac_1d_ad(x1D, dx1D, y1D, dy1D)
                     do k = nz-1, 2, -1
                        call AD_interpolation2D(dx1D%geop(k-1), dz(:,:,k), self%pg(iobs))
                        call AD_interpolation2D(dx1D%pres(k-1), dp(:,:,k), self%pg(iobs))
                        call AD_interpolation2D(dx1D%temp(k-1), dt(:,:,k), self%pg(iobs))
			call AD_interpolation2D(dx1D%shum(k-1), dqv(:,:,k), self%pg(iobs))
                     end do
		  end do
	       end do
	    end do
	 end if
	 !
	 call add_halo(info, dp, nx, ny, nz)
	 call add_halo(info, dt, nx, ny, nz)
	 call add_halo(info, dqv, nx, ny, nz)
	 call add_field(x, 'p', it, dp, nz)
	 call add_field(x, 't', it, dt, nz)
	 call add_field(x, 'qv', it, dqv, nz)
      end do
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(dx1D%temp, dx1D%shum, dx1D%pres, dx1D%geop)
      deallocate(y1D%refrac, y1D%geop, y1D%weights)
      !
      return
   end subroutine interpolate_DreT
   !
   !
   !
   subroutine interpolate_ba1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1dbangle) :: y1D
      ! --- ROPP variables
      real(r_size) :: ztmp, ptmp
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%bangle(1), y1D%impact(1), y1D%weights(1))
      !
      call get_field(x, 'z', 1, z, nz)
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'p', it, p, nz)
	 call get_field(x, 't', it, t, nz)
         call get_field(x, 'qv', it, qv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          y1D%nobs = 1
	          y1D%lon = obsspace%lon(iobs)
                  y1D%lat = obsspace%lat(iobs)
		  y1D%r_curve = obsspace%rcurve(iobs)
                  y1D%undulation = obsspace%undulation(iobs)
                  y1D%azimuth = obsspace%azimuth(iobs)
                  y1D%bangle(1) = obsspace%obs(iobs,1)
		  y1D%impact(1) = obsspace%impact(iobs)
                  ztmp = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
                  y1D%weights(1) = 1.0d0
                  do k = 2, nz - 1
                     call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                     call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                     call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                     call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                     if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
                  end do
                  x1D%non_ideal = .false.
                  x1D%lon = y1D%lon
                  x1D%lat = y1D%lat
		  x1D%new_bangle_op = .true.
		  call ropp_fm_bangle_1d(x1D, y1D)
                  y%field(iobs,1) = y1D%bangle(1)
		  !
		  if (y1D%bangle(1) < 0.0d0 ) then
		     valid%field(iobs,1) = 0
		     y%field(iobs,1) = 0.d0
		     cycle
		  end if
		  call interpolation3D(ptmp, self%pg(iobs), real(ztmp,r_sngl), p, z)
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%bangle, y1D%impact, y1D%weights)
      !
      return
   end subroutine interpolate_ba1
   !
   !
   !
   subroutine interpolate_ba2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      logical :: update
      integer :: nz, nt, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1dbangle) :: y1D
      ! --- ROPP variables
      real(r_size) :: ztmp, ptmp
      !
      nz = x%nlev; nt = obsspace%nt
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%bangle(1), y1D%impact(1), y1D%weights(1))
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
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
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
	       y1D%nobs = 1
	       y1D%lon = obsspace%lon(iobs)
               y1D%lat = obsspace%lat(iobs)
	       y1D%r_curve = obsspace%rcurve(iobs)
               y1D%undulation = obsspace%undulation(iobs)
               y1D%azimuth = obsspace%azimuth(iobs)
               y1D%bangle(1) = obsspace%obs(iobs,1)
	       y1D%impact(1) = obsspace%impact(iobs)
               ztmp = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
               y1D%weights(1) = 1.0d0
               do k = 2, nz - 1
                  call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                  call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                  call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                  call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                  if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
               end do
               x1D%non_ideal = .false.
               x1D%lon = y1D%lon
               x1D%lat = y1D%lat
	       x1D%new_bangle_op = .true.
	       call ropp_fm_bangle_1d(x1D, y1D)
               y%field(iobs,1) = y1D%bangle(1)
	       !
	       if (y1D%bangle(1) < 0.0d0 ) then
		  valid%field(iobs,1) = 0
		  y%field(iobs,1) = 0.d0
		  cycle
	       end if
	       call interpolation3D(ptmp, self%pg(iobs), real(ztmp,r_sngl), p, z)
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%bangle, y1D%impact, y1D%weights)
      !
      return
   end subroutine interpolate_ba2
   !
   !
   !
   subroutine interpolate_ba3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it
      real(r_size), dimension(:,:,:), allocatable, save :: z, p, t, qv
      ! --- ROPP types
      type(State1dFM) :: x1D
      type(Obs1dbangle) :: y1D
      ! --- ROPP variables
      real(r_size) :: ztmp, ptmp
      !
      nx = obsspace%nx; ny = obsspace%ny; nz = x%nlev
      if (.not. allocated(z)) allocate(z(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(y1D%bangle(1), y1D%impact(1), y1D%weights(1))
      !
      do it = 1, nt
         do i = ip, ip+1
	 do j = jp, jp+1
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    call get_field(x, 'z', ixyt, z(i,j,1:nz), nz)
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
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
	    y1D%nobs = 1
	    y1D%lon = obsspace%lon(iobs)
            y1D%lat = obsspace%lat(iobs)
	    y1D%r_curve = obsspace%rcurve(iobs)
            y1D%undulation = obsspace%undulation(iobs)
            y1D%azimuth = obsspace%azimuth(iobs)
            y1D%bangle(1) = obsspace%obs(iobs,1)
	    y1D%impact(1) = obsspace%impact(iobs)
            ztmp = geometric2geopotential(obsspace%lat(iobs), obsspace%height(iobs))
            y1D%weights(1) = 1.0d0
            do k = 2, nz - 1
               call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
               call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
               call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
               call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
               if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
            end do
            x1D%non_ideal = .false.
            x1D%lon = y1D%lon
            x1D%lat = y1D%lat
	    x1D%new_bangle_op = .true.
	    call ropp_fm_bangle_1d(x1D, y1D)
            y%field(iobs,1) = y1D%bangle(1)
	    !
	    if (y1D%bangle(1) < 0.0d0 ) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	       cycle
	    end if
	    call interpolation3D(ptmp, self%pg(iobs), real(ztmp,r_sngl), p, z)
	    if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
	      & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) then
	       valid%field(iobs,1) = 0
	       y%field(iobs,1) = 0.d0
	    else
	       valid%field(iobs,1) = 1
	    end if
	 end do
      end do
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(y1D%bangle, y1D%impact, y1D%weights)
      !
      return
   end subroutine interpolate_ba3
   !
   !
   !
   subroutine interpolate_Dba(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, p, t, qv, dz, dp, dt, dqv
      ! --- ROPP types
      type(State1dFM) :: x1D, dx1D
      type(Obs1dbangle) :: y1D
      ! --- ROPP variables
      real(r_size) :: dy1D(1)
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      dx1D%n_lev = nz - 2
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(dx1D%temp(dx1D%n_lev), dx1D%shum(dx1D%n_lev), dx1D%pres(dx1D%n_lev), dx1D%geop(dx1D%n_lev))
      allocate(y1D%bangle(1), y1D%impact(1), y1D%weights(1))
      !
      call get_field(xbck, 'z', 1, z, nz)
      !call get_field(x, 'z', 1, dz, nz)
      dz(:,:,:) = 0.0d0
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(xbck, 'p', it, p, nz)
	 call get_field(xbck, 't', it, t, nz)
         call get_field(xbck, 'qv', it, qv, nz)
	 call get_field(x, 'p', it, dp, nz)
	 call get_field(x, 't', it, dt, nz)
         call get_field(x, 'qv', it, dqv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
		  y1D%nobs = 1
	          y1D%lon = obsspace%lon(iobs)
                  y1D%lat = obsspace%lat(iobs)
                  y1D%r_curve = obsspace%rcurve(iobs)
                  y1D%undulation = obsspace%undulation(iobs)
                  y1D%azimuth = obsspace%azimuth(iobs)
		  y1D%g_sfc = gravity(obsspace%lat(iobs))
                  y1D%r_earth = r_eff(obsspace%lat(iobs))
                  y1D%bangle(1) = obsspace%obs(iobs,1)
		  y1D%impact(1) = obsspace%impact(iobs)
                  y1D%weights(1) = 1.0d0
                  do k = 2, nz - 1
                     call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                     call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                     call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                     call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                     if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
		     call TL_interpolation2D(dx1D%temp(k-1), dt(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%geop(k-1), dz(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%pres(k-1), dp(:,:,k), self%pg(iobs))
                     call TL_interpolation2D(dx1D%shum(k-1), dqv(:,:,k), self%pg(iobs))
                  end do
                  x1D%non_ideal = .false.
                  x1D%lon = y1D%lon
                  x1D%lat = y1D%lat
		  x1D%new_bangle_op = .true.
		  call ropp_fm_bangle_1d_tl(x1D, dx1D, y1D, dy1D)
                  y%field(iobs,1) = dy1D(1)
	          !
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
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(dx1D%temp, dx1D%shum, dx1D%pres, dx1D%geop)
      deallocate(y1D%bangle, y1D%impact, y1D%weights)
      !
      return
   end subroutine interpolate_Dba
   !
   !
   !
   subroutine interpolate_DbaT(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: z, p, t, qv, dz, dp, dt, dqv
      ! --- ROPP types
      type(State1dFM) :: x1D, dx1D
      type(Obs1dbangle) :: y1D
      ! --- ROPP variables
      real(r_size) :: dy1D(1)
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      x1D%n_lev = nz - 2 ! Except layers at kz=1 and kz=nz
      dx1D%n_lev = nz - 2
      allocate(x1D%temp(x1D%n_lev), x1D%shum(x1D%n_lev), x1D%pres(x1D%n_lev), x1D%geop(x1D%n_lev))
      allocate(dx1D%temp(dx1D%n_lev), dx1D%shum(dx1D%n_lev), dx1D%pres(dx1D%n_lev), dx1D%geop(dx1D%n_lev))
      allocate(y1D%bangle(1), y1D%impact(1), y1D%weights(1))
      !
      call get_field(xbck, 'z', 1, z, nz)
      do it = 1, nt
         dz = 0.0d0
	 dp = 0.0d0
         dt = 0.0d0
	 dqv = 0.d0
	 !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
	    call get_field(xbck, 'p', it, p, nz)
	    call get_field(xbck, 't', it, t, nz)
            call get_field(xbck, 'qv', it, qv, nz)
	    !
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     if (self%pg(iobs)%range11 == FlagOverTop .or. self%pg(iobs)%range21 == FlagOverTop .or. &
		       & self%pg(iobs)%range12 == FlagOverTop .or. self%pg(iobs)%range12 == FlagOverTop) cycle
		     if (valid%field(iobs,1) == 0) cycle
		     !
		     y1D%nobs = 1
		     y1D%lon = obsspace%lon(iobs)
                     y1D%lat = obsspace%lat(iobs)
                     y1D%r_curve = obsspace%rcurve(iobs)
                     y1D%undulation = obsspace%undulation(iobs)
                     y1D%azimuth = obsspace%azimuth(iobs)
		     y1D%g_sfc = gravity(obsspace%lat(iobs))
                     y1D%r_earth = r_eff(obsspace%lat(iobs))
                     y1D%bangle(1) = obsspace%obs(iobs,1)
		     y1D%impact(1) = obsspace%impact(iobs)
                     y1D%weights(1) = 1.0d0
                     do k = 2, nz - 1
                        call interpolation2D(x1D%temp(k-1), self%pg(iobs), t(:,:,k))
                        call interpolation2D(x1D%geop(k-1), self%pg(iobs), z(:,:,k))
                        call interpolation2D(x1D%pres(k-1), self%pg(iobs), p(:,:,k))
                        call interpolation2D(x1D%shum(k-1), self%pg(iobs), qv(:,:,k))
                        if(x1D%shum(k-1) <= 0.0d0) x1D%shum(k-1) = 0.000001d0
                     end do
                     x1D%non_ideal = .false.
                     x1D%lon = y1D%lon
                     x1D%lat = y1D%lat
                     x1D%new_bangle_op = .true.
		     dy1D(1) = y%field(iobs,1)
		     dx1D%temp(:) = 0.0d0
                     dx1D%shum(:) = 0.0d0
                     dx1D%pres(:) = 0.0d0
                     dx1D%geop(:) = 0.0d0 
		     call ropp_fm_bangle_1d_ad(x1D, dx1D, y1D, dy1D)
                     do k = nz-1, 2, -1
                        call AD_interpolation2D(dx1D%geop(k-1), dz(:,:,k), self%pg(iobs))
                        call AD_interpolation2D(dx1D%pres(k-1), dp(:,:,k), self%pg(iobs))
                        call AD_interpolation2D(dx1D%temp(k-1), dt(:,:,k), self%pg(iobs))
			call AD_interpolation2D(dx1D%shum(k-1), dqv(:,:,k), self%pg(iobs))
                     end do
		  end do
	       end do
	    end do
	 end if
	 !
	 call add_halo(info, dp, nx, ny, nz)
	 call add_halo(info, dt, nx, ny, nz)
	 call add_halo(info, dqv, nx, ny, nz)
	 call add_field(x, 'p', it, dp, nz)
	 call add_field(x, 't', it, dt, nz)
	 call add_field(x, 'qv', it, dqv, nz)
      end do
      deallocate(x1D%temp, x1D%shum, x1D%pres, x1D%geop)
      deallocate(dx1D%temp, dx1D%shum, dx1D%pres, dx1D%geop)
      deallocate(y1D%bangle, y1D%impact, y1D%weights)
      !
      return
   end subroutine interpolate_DbaT
   !
   !
   !
   subroutine apply_H_NodeHFieldGNSS1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 're') then
	 call interpolate_re1(self, info, x, obsspace, valid, y)
      else if (trim(self%name) == 'ba') then
	 call interpolate_ba1(self, info, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldGNSS1
   !
   !
   !
   subroutine apply_H_NodeHFieldGNSS2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 're') then
	 call interpolate_re2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      else if (trim(self%name) == 'ba') then
	 call interpolate_ba2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldGNSS2
   !
   !
   !
   subroutine apply_H_NodeHFieldGNSS3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 're') then
	 call interpolate_re3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      else if (trim(self%name) == 'ba') then
	 call interpolate_ba3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      end if
      !
      return
   end subroutine apply_H_NodeHFieldGNSS3
   !
   !
   !
   subroutine initialize_DH_NodeHFieldGNSS(self, info, x, obsspace, valid, y, ne)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      !
      if (trim(self%name) == 're') then
	 !call initialize_Dre(self, info, x, obsspace, valid, y)
      else if (trim(self%name) == 'ba') then
	 !call initialize_Dba(self, info, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine initialize_DH_NodeHFieldGNSS
   !
   !
   !
   subroutine apply_DH_NodeHFieldGNSS(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      if (trim(self%name) == 're') then
	 call interpolate_Dre(self, info, xbck, x, obsspace, valid, y)
      else if (trim(self%name) == 'ba') then
	 call interpolate_Dba(self, info, xbck, x, obsspace, valid, y)
      end if
      !
      return
   end subroutine apply_DH_NodeHFieldGNSS
   !
   !
   !
   subroutine apply_DHT_NodeHFieldGNSS(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldGNSS), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldGNSS), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      if (trim(self%name) == 're') then
	 call interpolate_DreT(self, info, xbck, obsspace, valid, y, x)
      else if (trim(self%name) == 'ba') then
	 call interpolate_DbaT(self, info, xbck, obsspace, valid, y, x)
      end if
      !
      return
   end subroutine apply_DHT_NodeHFieldGNSS
   !
   !
   !
end module NodeHFieldGNSS_class
