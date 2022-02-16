module NodeField_class
! Author: Le Duc
! Created date: 13 Jul 2014
   use variable, only : r_size
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   type NodeField
      character(len=10) :: name
      integer :: nx, ny, nlev, nt
      real(r_size), dimension(:,:,:,:), allocatable :: field
   end type NodeField
   !
   interface new
      module procedure new_NodeField
   end interface
   interface destroy
      module procedure destroy_NodeField
   end interface
   interface display
      module procedure display_NodeField
   end interface
   interface assignment(=)
      module procedure copy_NodeField
   end interface
   interface get_name
      module procedure get_name_NodeField
   end interface
   interface get_nx
      module procedure get_nx_NodeField
   end interface
   interface get_ny
      module procedure get_ny_NodeField
   end interface
   interface get_nlev
      module procedure get_nlev_NodeField
   end interface
   interface get_nt
      module procedure get_nt_NodeField
   end interface
   interface get_ndim
      module procedure get_ndim_NodeField
   end interface
   interface get_field
      module procedure get_field_NodeField1
      module procedure get_field_NodeField2
      module procedure get_field_NodeField3
      module procedure get_field_NodeField4
   end interface
   interface set_field
      module procedure set_field_NodeField1
      module procedure set_field_NodeField2
      module procedure set_field_NodeField3
      module procedure set_field_NodeField4
      module procedure set_field_NodeField5
   end interface
   interface add_field
      module procedure add_field_NodeField1
      module procedure add_field_NodeField2
   end interface
   interface add
      module procedure add_NodeField1
      module procedure add_NodeField2
   end interface
   interface subtract
      module procedure subtract_NodeField
   end interface
   interface power
      module procedure power_NodeField
   end interface
   interface ratransform
      module procedure ratransform_NodeField
   end interface
   interface multiply
      module procedure multiply_NodeField1
      module procedure multiply_NodeField2
   end interface
   interface divide
      module procedure divide_NodeField1
      module procedure divide_NodeField2
   end interface
   interface compute_normsquare
      module procedure compute_normsquare_NodeField1
      module procedure compute_normsquare_NodeField2
   end interface
   interface innerproduct
      module procedure innerproduct_NodeField1
      module procedure innerproduct_NodeField2
   end interface
   interface localize
      module procedure localize_NodeField1
      module procedure localize_NodeField2
   end interface
   interface broadcast_bck
      module procedure broadcast_bck_NodeField
   end interface
   interface allreduce_ens
      module procedure allreduce_ens_NodeField
   end interface
   interface randomize
      module procedure randomize_NodeField1
      module procedure randomize_NodeField2
   end interface
   !
contains
   !
   subroutine new_NodeField(self, name, nx, ny, nlev, nt)
      implicit none
      type(NodeField), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: nx, ny, nlev, nt
      !
      self%name = name
      self%nx = nx
      self%ny = ny
      self%nlev = nlev
      self%nt = nt
      allocate(self%field(nx,ny,nlev,nt))
      self%field(:,:,:,:) = 0.d0
      !
      return
   end subroutine new_NodeField
   !
   !
   !
   subroutine destroy_NodeField(self)
      implicit none
      type(NodeField), intent(inout) :: self
      !
      if (allocated(self%field)) deallocate(self%field)
      !
      return
   end subroutine destroy_NodeField
   !
   !
   !
   subroutine display_NodeField(self)
      implicit none
      type(NodeField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nx, self%ny, self%nlev, self%nt
      !
      return
   end subroutine display_NodeField
   !
   !
   !
   subroutine copy_NodeField(self, f)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeField), intent(in) :: f
      !
      self%field(:,:,:,:) = f%field(:,:,:,:)
      !
      return
   end subroutine copy_NodeField
   !
   !
   !
   subroutine get_name_NodeField(self, name)
      implicit none
      type(NodeField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeField
   !
   !
   !
   subroutine get_nx_NodeField(self, nx)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(out) :: nx
      !
      nx = self%nx
      !
      return
   end subroutine get_nx_NodeField
   !
   !
   !
   subroutine get_ny_NodeField(self, ny)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(out) :: ny
      !
      ny = self%ny
      !
      return
   end subroutine get_ny_NodeField
   !
   !
   !
   subroutine get_nlev_NodeField(self, nlev)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(out) :: nlev
      !
      nlev = self%nlev
      !
      return
   end subroutine get_nlev_NodeField
   !
   !
   !
   subroutine get_nt_NodeField(self, nt)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(out) :: nt
      !
      nt = self%nt
      !
      return
   end subroutine get_nt_NodeField
   !
   !
   !
   subroutine get_ndim_NodeField(self, ndim)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(out) :: ndim
      !
      ndim = self%nx*self%ny*self%nlev*self%nt
      !
      return
   end subroutine get_ndim_NodeField
   !
   !
   !
   subroutine get_field_NodeField1(self, it, field)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(inout) :: field
      !
      field(:,:,:) = self%field(:,:,:,it)
      !
      return
   end subroutine get_field_NodeField1
   !
   !
   !
   subroutine get_field_NodeField2(self, field)
      implicit none
      type(NodeField), intent(in) :: self
      real(r_size), dimension(self%nx,self%ny,self%nlev,self%nt), intent(out) :: field
      !
      field(:,:,:,:) = self%field(:,:,:,:)
      !
      return
   end subroutine get_field_NodeField2
   !
   !
   !
   subroutine get_field_NodeField3(self, itout, object)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(in) :: itout
      type(NodeField), intent(inout) :: object
      !
      object%field(:,:,:,1) = self%field(:,:,:,itout)
      !
      return
   end subroutine get_field_NodeField3
   !
   !
   !
   subroutine get_field_NodeField4(self, is, ie, js, je, it, field)
      implicit none
      type(NodeField), intent(in) :: self
      integer, intent(in) :: is, ie, js, je, it
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(inout) :: field
      !
      field(is:ie,js:je,:) = self%field(is:ie,js:je,:,it)
      !
      return
   end subroutine get_field_NodeField4
   !
   !
   !
   subroutine set_field_NodeField1(self, it, field)
      implicit none
      type(NodeField), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), intent(in) :: field
      !
      self%field(:,:,:,it) = field
      !
      return
   end subroutine set_field_NodeField1
   !
   !
   !
   subroutine set_field_NodeField2(self, it, field)
      implicit none
      type(NodeField), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(in) :: field
      integer :: i, j, k
      !
      !do k = 1, self%nlev
         !do j = 1, self%ny
	    !do i = 1, self%nx
	       !if (field(i,j,k) /= field(i,j,k)) then ! NaN
	          !self%field(i,j,k,it) = 0.d0
	       !else
		  !self%field(i,j,k,it) = field(i,j,k)
	       !end if
	    !end do
	 !end do
      !end do
      self%field(:,:,:,it) = field(:,:,:)
      !
      return
   end subroutine set_field_NodeField2
   !
   !
   !
   subroutine set_field_NodeField3(self, field)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: field
      !
      self%field(:,:,:,:) = field
      !
      return
   end subroutine set_field_NodeField3
   !
   !
   !
   subroutine set_field_NodeField4(self, field)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), dimension(self%nx,self%ny,self%nlev,self%nt), intent(in) :: field
      !
      self%field(:,:,:,:) = field
      !
      return
   end subroutine set_field_NodeField4
   !
   !
   !
   subroutine set_field_NodeField5(self, is, ie, js, je, it, field)
      implicit none
      type(NodeField), intent(inout) :: self
      integer, intent(in) :: is, ie, js, je, it
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(in) :: field
      !
      self%field(is:ie,js:je,:,it) = field(is:ie,js:je,:)
      !
      return
   end subroutine set_field_NodeField5
   !
   !
   !
   subroutine add_field_NodeField1(self, it, field)
      implicit none
      type(NodeField), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), intent(in) :: field
      !
      self%field(:,:,:,it) = self%field(:,:,:,it) + field
      !
      return
   end subroutine add_field_NodeField1
   !
   !
   !
   subroutine add_field_NodeField2(self, it, field)
      implicit none
      type(NodeField), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(in) :: field
      integer :: i, j, k
      !
      !do k = 1, self%nlev
         !do j = 1, self%ny
	    !do i = 1, self%nx
	       !if (field(i,j,k) /= field(i,j,k)) then ! NaN
	          !self%field(i,j,k,it) = self%field(i,j,k,it) + 0.d0
	       !else
		  !self%field(i,j,k,it) = self%field(i,j,k,it) + field(i,j,k)
	       !end if
	    !end do
	 !end do
      !end do
      self%field(:,:,:,it) = self%field(:,:,:,it) + field(:,:,:)
      !
      return
   end subroutine add_field_NodeField2
   !
   !
   !
   subroutine add_NodeField1(self, const)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:) + const
      self%field(:,:,:,:) = self%field(:,:,:,:) + const
      !
      return
   end subroutine add_NodeField1
   !
   !
   !
   subroutine add_NodeField2(self, object)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeField), intent(in) :: object
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:) + object%field(:,:,1:nlev,:)
      self%field(:,:,:,:) = self%field(:,:,:,:) + object%field(:,:,:,:)
      !
      return
   end subroutine add_NodeField2
   !
   !
   !
   subroutine subtract_NodeField(self, object)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeField), intent(in) :: object
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:) - object%field(:,:,1:nlev,:)
      self%field(:,:,:,:) = self%field(:,:,:,:) - object%field(:,:,:,:)
      !
      return
   end subroutine subtract_NodeField
   !
   !
   !
   subroutine power_NodeField(self, const)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:)**const
      self%field(:,:,:,:) = self%field(:,:,:,:)**const
      !
      return
   end subroutine power_NodeField
   !
   !
   !
   subroutine ratransform_NodeField(self, threshold)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: threshold
      real(r_size) :: a, b
      !
      a = 3.d0*threshold**(2.d0/3.d0)
      b = -2.d0*threshold
      where (self%field >= threshold)
         self%field = a*self%field**(1.d0/3.d0) + b
      end where
      !
      return
   end subroutine ratransform_NodeField
   !
   !
   !
   subroutine multiply_NodeField1(self, const)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = const*self%field(:,:,1:nlev,:)
      self%field(:,:,:,:) = const*self%field(:,:,:,:)
      !
      return
   end subroutine multiply_NodeField1
   !
   !
   !
   subroutine multiply_NodeField2(self, object)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeField), intent(in) :: object
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = object%field(:,:,1:nlev,:)*self%field(:,:,1:nlev,:)
      self%field(:,:,:,:) = object%field(:,:,:,:)*self%field(:,:,:,:)
      !
      return
   end subroutine multiply_NodeField2
   !
   !
   !
   subroutine divide_NodeField1(self, const)
      implicit none
      type(NodeField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:)/const
      self%field(:,:,:,:) = self%field(:,:,:,:)/const
      !
      return
   end subroutine divide_NodeField1
   !
   !
   !
   subroutine divide_NodeField2(self, object)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeField), intent(in) :: object
      !integer :: nlev
      !
      !nlev = self%nlev
      !if (nlev > 1) nlev = nlev - 1
      !where (object%field(:,:,1:nlev,:) < 1.e-12)
         !self%field(:,:,1:nlev,:) = 0.d0
      !else where
         !self%field(:,:,1:nlev,:) = self%field(:,:,1:nlev,:)/object%field(:,:,1:nlev,:)
      !end where
      where (object%field(:,:,:,:) < 1.e-12)
         self%field(:,:,:,:) = 0.d0
      else where
         self%field(:,:,:,:) = self%field(:,:,:,:)/object%field(:,:,:,:)
      end where
      !
      return
   end subroutine divide_NodeField2
   !
   !
   !
   subroutine compute_normsquare_NodeField1(self, info, normsquare)
   ! excluding top boundaries
      implicit none
      type(NodeField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      real(r_size), intent(out) :: normsquare
      integer :: nx, ny, nlev, dis, die, djs, dje, is, ie, js, je
      !
      nx = self%nx; ny = self%ny; nlev = self%nlev
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !if (nlev > 1) nlev = nlev - 1
      normsquare = sum(self%field(is:ie,js:je,1:nlev,:)**2)
      !
      return
   end subroutine compute_normsquare_NodeField1
   !
   !
   !
   subroutine compute_normsquare_NodeField2(self, info, klev, normsquare)
      implicit none
      type(NodeField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: klev
      real(r_size), intent(out) :: normsquare
      integer :: nx, ny, dis, die, djs, dje, is, ie, js, je
      !
      nx = self%nx; ny = self%ny
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      normsquare = sum(self%field(is:ie,js:je,1:klev,:)**2)
      !
      return
   end subroutine compute_normsquare_NodeField2
   !
   !
   !
   subroutine innerproduct_NodeField1(self, object, info, product)
   ! excluding bottom and top boundaries
      implicit none
      type(NodeField), intent(in) :: self
      type(NodeField), intent(in) :: object
      type(NodeInfo), intent(in) :: info
      real(r_size), intent(out) :: product
      integer :: nx, ny, nlev, dis, die, djs, dje, is, ie, js, je
      !
      nx = self%nx; ny = self%ny; nlev = self%nlev
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !if (nlev > 1) nlev = nlev - 1
      product = sum(self%field(is:ie,js:je,1:nlev,:)*object%field(is:ie,js:je,1:nlev,:))
      !
      return
   end subroutine innerproduct_NodeField1
   !
   !
   !
   subroutine innerproduct_NodeField2(self, object, info, nlev, product)
   ! excluding bottom and top boundaries
      implicit none
      type(NodeField), intent(in) :: self
      type(NodeField), intent(in) :: object
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev
      real(r_size), intent(out) :: product
      integer :: nx, ny, dis, die, djs, dje, is, ie, js, je
      !
      nx = self%nx; ny = self%ny
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      product = sum(self%field(is:ie,js:je,1:nlev,:)*object%field(is:ie,js:je,1:nlev,:))
      !
      return
   end subroutine innerproduct_NodeField2
   !
   !
   !
   subroutine localize_NodeField1(self, info, i0, j0)
      use variable, only : hscale
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: i0, j0
      integer :: nx, ny, i, j, ig, jg
      real(r_size) :: dh, hweight
      !
      nx = self%nx; ny = self%ny
      do i = 1, nx
         call get_xindex(info, i, ig)
         do j = 1, ny
	    call get_yindex(info, j, jg)
	    dh = sqrt(1.d0*(ig-i0)**2+1.d0*(jg-j0)**2)
	    if (dh > 0.d0) then
	       call GaspariCohn(dh, hscale, hweight)
	       self%field(i,j,:,:) = hweight*self%field(i,j,:,:)
	    end if
	 end do
      end do
      !
      return
   end subroutine localize_NodeField1
   !
   !
   !
   subroutine localize_NodeField2(self, info, i0, j0, R0)
      use variable, only : hscale
      use enkflib, only : GaspariCohn
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: i0, j0
      real(r_size), intent(in) :: R0
      integer :: nx, ny, i, j, ig, jg
      real(r_size) :: dh, hweight
      !
      nx = self%nx; ny = self%ny
      do i = 1, nx
         call get_xindex(info, i, ig)
         do j = 1, ny
	    call get_yindex(info, j, jg)
	    dh = sqrt(1.d0*(ig-i0)**2+1.d0*(jg-j0)**2) - R0
	    if (dh > 0.d0) then
	       call GaspariCohn(dh, hscale, hweight)
	       self%field(i,j,:,:) = hweight*self%field(i,j,:,:)
	    end if
	 end do
      end do
      !
      return
   end subroutine localize_NodeField2
   !
   !
   !
   subroutine broadcast_bck_NodeField(self, info, source_ide)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: source_ide
      integer :: istat(MPI_STATUS_SIZE)
      integer :: myidx, myidy, myide, nxpe, nype, nepe
      integer :: ierror, ide, destination, source
      !
      call get_myidx(info, myidx); call get_nxpe(info, nxpe)
      call get_myidy(info, myidy); call get_nype(info, nype)
      call get_myide(info, myide); call get_nepe(info, nepe)
      if (myide == source_ide) then
	 do ide = 0, nepe-1
            if (ide == myide) cycle
            destination = ide*nxpe*nype + myidy*nxpe + myidx
	    call MPI_SEND(self%field, self%nx*self%ny*self%nlev*self%nt, r_type, destination, 0, MPI_COMM_WORLD, ierror)
         end do
      else
	 source = source_ide*nxpe*nype + myidy*nxpe + myidx
         call MPI_RECV(self%field, self%nx*self%ny*self%nlev*self%nt, r_type, source, 0, MPI_COMM_WORLD, istat, ierror)
      end if
      !
      return
   end subroutine broadcast_bck_NodeField
   !
   !
   !
   subroutine allreduce_ens_NodeField(self)
      implicit none
      type(NodeField), intent(inout) :: self
      integer :: it
      !
      do it = 1, self%nt
         call allreduce3D('e', self%field(:,:,:,it), self%nx, self%ny, self%nlev)
      end do
      !
      return
   end subroutine allreduce_ens_NodeField
   !
   !
   !
   subroutine dirac(self, info, i0, j0, k0, it0)
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: i0, j0, k0, it0
      integer :: nproc, myid, myide, nx0, ny0, nz0, nx, ny, i, j, k
      real(r_size), dimension(:,:,:), allocatable :: datum
      !
      call get_myid(info, myid, nproc)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nz0 = self%nlev
      call get_myide(info, myide)
      call get_nx(info, nx)
      call get_ny(info, ny)
      allocate(datum(nx0,ny0,nz0))
      !
      if (myid == 0) then
	 datum(:,:,:) = 0.d0
         datum(i0,j0,k0) = 1.d0
      end if
      if (myide == 0) call scatter(info, 'xy', datum, 0, nx0, ny0, nz0)
      call broadcast3D('e', datum(1:nx,1:ny,:), 0, nx, ny, nz0)
      call set_field(self, it0, datum(1:nx,1:ny,:))
      !do i = 1, nx
         !do j = 1, ny
	    !do k = 1, nz0
	       !if (datum(i,j,k) > 0.5) print*, 'SCATTER:', myid, i, j, k
	    !end do
	 !end do
      !end do
      deallocate(datum)
      !
      return
   end subroutine dirac
   !
   !
   !
   subroutine randomize_NodeField1(self, info)
      use random, only : random_normal
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer :: nproc, myid, myide, nx0, ny0, nz0, nt0, nx, ny, i, j, k, it
      real(r_size), dimension(:,:,:), allocatable :: datum
      !
      call get_myid(info, myid, nproc)
      call get_myide(info, myide)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nz0 = self%nlev
      nt0 = self%nt
      call get_nx(info, nx)
      call get_ny(info, ny)
      allocate(datum(nx0,ny0,nz0))
      !
      do it = 1, nt0
         if (myid == 0) then
	    do k = 1, nz0
	       do j = 1, ny0
	          do i = 1, nx0
		     datum(i,j,k) = 1.d0*random_normal()
	          end do
	       end do
	    end do 
         end if
         if (myide == 0) call scatter(info, 'xy', datum, 0, nx0, ny0, nz0)
	 call broadcast3D('e', datum(1:nx,1:ny,:), 0, nx, ny, nz0)
         call set_field(self, it, datum(1:nx,1:ny,:))
      end do
      deallocate(datum)
      !
      return
   end subroutine randomize_NodeField1
   !
   !
   !
   subroutine randomize_NodeField2(self, info, it0)
      use random, only : random_normal
      implicit none
      type(NodeField), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: it0
      integer :: nproc, myid, myide, nx0, ny0, nz0, nx, ny, i, j, k
      real(r_size), dimension(:,:,:), allocatable :: datum
      !
      call get_myid(info, myid, nproc)
      call get_myide(info, myide)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nz0 = self%nlev
      call get_nx(info, nx)
      call get_ny(info, ny)
      allocate(datum(nx0,ny0,nz0))
      !
      call set_field(self, 0.d0)
      if (myid == 0) then
	 do k = 1, nz0
	    do j = 1, ny0
	       do i = 1, nx0
		  datum(i,j,k) = 1.d0*random_normal()
	       end do
	    end do
	 end do 
      end if
      if (myide == 0) call scatter(info, 'xy', datum, 0, nx0, ny0, nz0)
      call broadcast3D('e', datum(1:nx,1:ny,:), 0, nx, ny, nz0)
      call set_field(self, it0, datum(1:nx,1:ny,:))
      deallocate(datum)
      !
      return
   end subroutine randomize_NodeField2
   !
   !
   !
end module NodeField_class
   
