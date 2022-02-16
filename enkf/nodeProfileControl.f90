module NodeProfileControl_class
! Author: Le Duc
! Created date: 27 Jul 2014
   use variable, only : nsoil, r_size, r_sngl
   use NodeInfo_class
   use NodeObsSpaceControl_class
   use NodeProfileField_class
   use NodeControl_class
   use NodeMPI
   implicit none
   !
   type NodeProfileFieldPointer
      type(NodeProfileField), pointer :: p
   end type NodeProfileFieldPointer
   !
   type NodeProfileControl
      integer :: ncontrol
      integer :: nxyt, nlev, nx, ny
      type(NodeProfileField), pointer :: u, v, w, t, p, qv, qc, qi, qr, qs, qg, tsoil, &
                                       & logp, rh, pwv, z, g, ps, logps, pmsl, rain, tgrd
      type(NodeProfileFieldPointer), dimension(:), allocatable :: control
   end type NodeProfileControl
   !
   interface new
      module procedure new_NodeProfileControl
   end interface
   interface destroy
      module procedure destroy_NodeProfileControl
   end interface
   interface display
      module procedure display_NodeProfileControl
   end interface
   interface assignment(=)
      module procedure copy_NodeProfileControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeProfileControl
   end interface
   interface get_nxyt
      module procedure get_nxyt_NodeProfileControl
   end interface
   interface get_name
      module procedure get_name_NodeProfileControl
   end interface
   interface get_nlev
      module procedure get_nlev_NodeProfileControl1
      module procedure get_nlev_NodeProfileControl2
      module procedure get_nlev_NodeProfileControl3
   end interface
   interface get_field
      module procedure get_field_NodeProfileControl1
      module procedure get_field_NodeProfileControl2
      module procedure get_field_NodeProfileControl3
      module procedure get_field_NodeProfileControl4
   end interface
   interface get_control
      module procedure get_control_NodeProfileControl
   end interface
   interface set_field
      module procedure set_field_NodeProfileControl1
      module procedure set_field_NodeProfileControl2
      module procedure set_field_NodeProfileControl3
      module procedure set_field_NodeProfileControl4
      module procedure set_field_NodeProfileControl5
   end interface
   interface set_control
      module procedure set_control_NodeProfileControl
   end interface
   interface add
      module procedure add_NodeProfileControl
   end interface
   interface subtract
      module procedure subtract_NodeProfileControl
   end interface
   interface multiply
      module procedure multiply_NodeProfileControl1
      module procedure multiply_NodeProfileControl2
   end interface
   interface divide
      module procedure divide_NodeProfileControl1
      module procedure divide_NodeProfileControl2
   end interface
   interface allreduce_ens
      module procedure allreduce_ens_NodeProfileControl
   end interface
   interface read_ensemble
      module procedure read_ensemble_NodeProfileControl1
      module procedure read_ensemble_NodeProfileControl2
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_nhm_NodeProfileControl(self, nxyt, nlev, nx, ny)
   ! the surface fields are imposed at the bottom level.
   ! the p field is non-hydrostatic pressure
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: nxyt, nlev, nx, ny
      integer :: ivar
      !
      self%ncontrol = 12
      self%nxyt = nxyt
      self%nlev = nlev
      self%nx = nx
      self%ny = ny
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! w
      allocate(self%w)
      call new(self%w, 'w', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%w
      ! t
      allocate(self%t)
      call new(self%t, 't', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      ! qc
      allocate(self%qc)
      call new(self%qc, 'qc', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qc
      ! qi
      allocate(self%qi)
      call new(self%qi, 'qi', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qi
      ! qr
      allocate(self%qr)
      call new(self%qr, 'qr', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qr
      ! qs
      allocate(self%qs)
      call new(self%qs, 'qs', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qs
      ! qg
      allocate(self%qg)
      call new(self%qg, 'qg', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qg
      ! tsoil
      allocate(self%tsoil)
      call new(self%tsoil, 'tsoil', nxyt, nsoil, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%tsoil
      !
      return
   end subroutine new_nhm_NodeProfileControl
   !
   !
   !
   subroutine new_obs_NodeProfileControl(self, nxyt, nlev, nx, ny)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: nxyt, nlev, nx, ny
      integer :: ivar
      !
      self%ncontrol = 14
      self%nxyt = nxyt
      self%nlev = nlev
      self%nx = nx
      self%ny = ny
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! logp
      allocate(self%logp)
      call new(self%logp, 'logp', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%logp
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      ! rh
      allocate(self%rh)
      call new(self%rh, 'rh', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%rh
      ! pwv
      allocate(self%pwv)
      call new(self%pwv, 'pwv', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%pwv
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      ! logps
      allocate(self%logps)
      call new(self%logps, 'logps', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%logps
      ! pmsl
      allocate(self%pmsl)
      call new(self%pmsl, 'pmsl', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%pmsl
      ! rain
      allocate(self%rain)
      call new(self%rain, 'rain', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%rain
      ! z
      allocate(self%z)
      call new(self%z, 'z', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%z
      ! g
      allocate(self%g)
      call new(self%g, 'g', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%g
      !
      return
   end subroutine new_obs_NodeProfileControl
   !
   !
   !
   subroutine new_enkf5_NodeProfileControl(self, nxyt, nlev, nx, ny)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: nxyt, nlev, nx, ny
      integer :: ivar
      !
      self%ncontrol = 6
      self%nxyt = nxyt
      self%nlev = nlev
      self%nx = nx
      self%ny = ny
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nxyt, nlev, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nxyt, 1, nx, ny)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      !
      return
   end subroutine new_enkf5_NodeProfileControl
   !
   !
   !
   subroutine new_NodeProfileControl(self, control_mode, nxyt, nlev, nx, ny)
      use variable, only : enkf0, enkf1, enkf2, enkf3, enkf4, enkf5
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: control_mode, nxyt, nlev, nx, ny
      !
      if (control_mode == enkf0) then
	 !call new_enkf0_NodeProfileControl(self, nxyt, nlev, nx, ny)
      else if (control_mode == enkf1) then
	 !call new_enkf1_NodeProfileControl(self, nxyt, nlev, nx, ny)
      else if (control_mode == enkf2) then
	 !call new_enkf2_NodeProfileControl(self, nxyt, nlev, nx, ny)
      else if (control_mode == enkf3) then
	 !call new_enkf3_NodeProfileControl(self, nxyt, nlev, nx, ny)
      else if (control_mode == enkf4) then
	 !call new_enkf4_NodeProfileControl(self, nxyt, nlev, nx, ny)
      else if (control_mode == enkf5) then
	 call new_enkf5_NodeProfileControl(self, nxyt, nlev, nx, ny)
      end if
      !
      return
   end subroutine new_NodeProfileControl
   !
   !
   !
   subroutine destroy_NodeProfileControl(self)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeProfileControl
   !
   !
   !
   subroutine display_NodeProfileControl(self)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      !
      print*, 'Control: ', self%ncontrol
      print*, 'Dimension: ', self%nxyt, self%nlev
      !
      return
   end subroutine display_NodeProfileControl
   !
   !
   !
   subroutine get_ncontrol_NodeProfileControl(self, ncontrol)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeProfileControl
   !
   !
   !
   subroutine get_nxyt_NodeProfileControl(self, nxyt)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(out) :: nxyt
      !
      nxyt = self%nxyt
      !
      return
   end subroutine get_nxyt_NodeProfileControl
   !
   !
   !
   subroutine get_name_NodeProfileControl(self, ivar, name)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeProfileControl
   !
   !
   !
   subroutine get_nlev_NodeProfileControl1(self, ivar, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nlev
      !
      call get_nlev(self%control(ivar)%p, nlev)
      !
      return
   end subroutine get_nlev_NodeProfileControl1
   !
   !
   !
   subroutine get_nlev_NodeProfileControl2(self, name, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: nlev
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_nlev(self%control(ivar)%p, nlev)
      !
      return
   end subroutine get_nlev_NodeProfileControl2
   !
   !
   !
   subroutine get_nlev_NodeProfileControl3(self, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(out) :: nlev
      integer :: ivar, nz
      !
      nlev = 0
      do ivar = 1, self%ncontrol
	 call get_nlev(self%control(ivar)%p, nz)
	 nlev = nlev + nz
      end do
      !
      return
   end subroutine get_nlev_NodeProfileControl3
   !
   !
   !
   subroutine copy_NodeProfileControl(self, f)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeProfileControl), intent(in) :: f
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p = f%control(ivar)%p
      end do
      !
      return
   end subroutine copy_NodeProfileControl
   !
   !
   !
   subroutine get_field_NodeProfileControl1(self, name, it, k2ijt, field, nx, ny, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it, nx, ny, nlev
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      real(r_size), dimension(nx,ny,nlev), intent(inout) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_field(self%control(ivar)%p, it, k2ijt, field, nx, ny)
      !
      return
   end subroutine get_field_NodeProfileControl1
   !
   !
   !
   subroutine get_field_NodeProfileControl2(self, name, ixyt, field, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: ixyt, nlev
      real(r_size), dimension(nlev), intent(inout) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_field(self%control(ivar)%p, ixyt, field)
      !
      return
   end subroutine get_field_NodeProfileControl2
   !
   !
   !
   subroutine get_field_NodeProfileControl3(self, name, ixyt, i, j, field, nx, ny, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: ixyt, i, j, nx, ny, nlev
      real(r_size), dimension(nx,ny,nlev), intent(inout) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_field(self%control(ivar)%p, ixyt, i, j, field, nx, ny)
      !
      return
   end subroutine get_field_NodeProfileControl3
   !
   !
   !
   subroutine get_field_NodeProfileControl4(self, field, nlev)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      integer, intent(in) :: nlev
      real(r_size), dimension(self%nxyt,nlev), intent(inout) :: field
      integer :: ivar, nz, k
      !
      if (self%nxyt == 0) return
      k = 1
      do ivar = 1, self%ncontrol
	 call get_nlev(self%control(ivar)%p, nz)
	 call get_field(self%control(ivar)%p, field(:,k:k+nz-1))
	 k = k + nz
      end do
      !
      return
   end subroutine get_field_NodeProfileControl4
   !
   !
   !
   subroutine get_control_NodeProfileControl(self, w, zw)
      implicit none
      type(NodeProfileControl), intent(in) :: self
      type(NodeControl), intent(in) :: w
      type(NodeProfileControl), intent(inout) :: zw
      integer :: nlev, ivar, iweight, ivarnhm
      !
      if (self%nxyt == 0) return
      do iweight = 1, w%ncontrol
         call get_nlev(w, iweight, nlev)
         !if (nlev > 1) nlev = nlev - 1
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            zw%control(iweight)%p%field(:,1:nlev) = self%control(ivarnhm)%p%field(:,1:nlev)
         end do
      end do
      !
      return
   end subroutine get_control_NodeProfileControl
   !
   !
   !
   subroutine set_field_NodeProfileControl1(self, const)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine set_field_NodeProfileControl1
   !
   !
   !
   subroutine set_field_NodeProfileControl2(self, name, ixyt, field, nlev)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: ixyt, nlev
      real(r_size), dimension(nlev), intent(in) :: field
      character(len=10) :: varname
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call set_field(self%control(ivar)%p, ixyt, field)
      !
      return
   end subroutine set_field_NodeProfileControl2
   !
   !
   !
   subroutine set_field_NodeProfileControl3(self, ivar, it, k2ijt, field, nx, ny, nlev)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: ivar, it, nx, ny, nlev
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      real(r_size), dimension(nx,ny,nlev), intent(in) :: field
      !
      call set_field(self%control(ivar)%p, it, k2ijt, field, nx, ny)
      !
      return
   end subroutine set_field_NodeProfileControl3
   !
   !
   !
   subroutine set_field_NodeProfileControl4(self, processed, k2ijt, obsspace, object)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeProfileControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_field_NodeProfileControl4
   !
   !
   !
   subroutine set_field_NodeProfileControl5(self, field, nlev)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer, intent(in) :: nlev
      real(r_size), dimension(self%nxyt,nlev), intent(in) :: field
      integer :: ivar, nz, k
      !
      if (self%nxyt == 0) return
      k = 1
      do ivar = 1, self%ncontrol
	 call get_nlev(self%control(ivar)%p, nz)
	 call set_field(self%control(ivar)%p, field(:,k:k+nz-1))
	 k = k + nz
      end do
      !
      return
   end subroutine set_field_NodeProfileControl5
   !
   !
   !
   subroutine set_control_NodeProfileControl(self, w, zw)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeControl), intent(in) :: w
      type(NodeProfileControl), intent(in) :: zw
      integer :: nlev, ivar, iweight, ivarnhm
      !
      if (self%nxyt == 0) return
      do iweight = 1, w%ncontrol
         call get_nlev(w, iweight, nlev)
         !if (nlev > 1) nlev = nlev - 1
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            self%control(ivarnhm)%p%field(:,1:nlev) = zw%control(iweight)%p%field(:,1:nlev)
         end do
      end do
      !
      return
   end subroutine set_control_NodeProfileControl
   !
   !
   !
   subroutine add_NodeProfileControl(self, object)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeProfileControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine add_NodeProfileControl
   !
   !
   !
   subroutine subtract_NodeProfileControl(self, object)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeProfileControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine subtract_NodeProfileControl
   !
   !
   !
   subroutine multiply_NodeProfileControl1(self, const)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine multiply_NodeProfileControl1
   !
   !
   !
   subroutine multiply_NodeProfileControl2(self, object)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeProfileControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine multiply_NodeProfileControl2
   !
   !
   !
   subroutine divide_NodeProfileControl1(self, const)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine divide_NodeProfileControl1
   !
   !
   !
   subroutine divide_NodeProfileControl2(self, object)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      type(NodeProfileControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine divide_NodeProfileControl2
   !
   !
   !
   subroutine allreduce_ens_NodeProfileControl(self)
      implicit none
      type(NodeProfileControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call allreduce_ens(self%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ens_NodeProfileControl
   !
   !
   !
   subroutine read_ensemble_NodeProfileControl1(self, info, ie, k2ijt, xpert, header, nt, ne, nxyt)
      use variable, only : itout
      implicit none
      integer, intent(in) :: ie, nt, ne, nxyt
      type(NodeProfileControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeControl), intent(inout) :: xpert
      character(*), intent(in) :: header
      !
      character(6) :: infile = '000000'
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, nx0, ny0, nx, ny, nz, nvar, ngroup
      integer :: igroup, it, imember, ivar, id, ierror, dis, die, djs, dje
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_r
      !
      call get_nx(info, nx); call get_ny(info, ny)
      call get_di(info, dis, die); call get_dj(info, djs, dje)
      call get_nx0(info, nx0); call get_ny0(info, ny0)
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype); call get_nepe(info, nepe)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      nvar = self%ncontrol
      ngroup = nt/nproc + 1
      if (mod(nt,nproc) == 0) ngroup = ngroup - 1
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self, ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum(nx0,ny0,nz), datum_r(nx0,ny0,nz,nvar))
      !
      imember = myide*ne + ie
      do igroup = 1, ngroup
	 it = (igroup-1)*nproc + myid+1
	 if (it <= nt) then
	    write(infile(1:2),'(I2.2)') it
	    write(infile(3:6),'(I4.4)') imember
	    !print*, 'Read ', trim(header)//infile
	    open(90, file=trim(header)//infile, form='unformatted')
	    do ivar = 1, nvar
	       read(90) datum_r(:,:,1:nlev(ivar),ivar)
	    end do
	    close(90)
	 end if
	 !
	 do id = 1, nproc
	    it = (igroup-1)*nproc + id
	    if (it > nt) cycle
	    do ivar = 1, nvar
	       if (myid == id-1) datum(:,:,1:nlev(ivar)) = datum_r(:,:,1:nlev(ivar),ivar)
	       call scatter(info, 'xy', datum(:,:,1:nlev(ivar)), id-1, nx0, ny0, nlev(ivar))
	       call set_field(self, ivar, it, k2ijt, datum(1:nx,1:ny,1:nlev(ivar)), nx, ny, nlev(ivar))
	       if (it == nt) then
		  call set_field(xpert, ivar, 2, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)), nlev(ivar))
	       else
		  if (it == itout) call set_field(xpert, ivar, 1, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)), nlev(ivar))
	       end if
	    end do
	 end do
      end do
      deallocate(nlev, datum_r, datum)
      !
      return
   end subroutine read_ensemble_NodeProfileControl1
   !
   !
   !
   subroutine read_ensemble_NodeProfileControl2(self, info, k2ijt, xpert, header, nt, ne, nxyt)
      use variable, only : itout
      implicit none
      integer, intent(in) :: nt, ne, nxyt
      type(NodeProfileControl), dimension(ne), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeControl), dimension(ne), intent(inout) :: xpert
      character(*), intent(in) :: header
      !
      character(6) :: infile = '000000'
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, nx0, ny0, ne0, nx, ny, nz, nvar, ngroup
      integer :: igroup, itm, it, imember, ivar, id, ide, ierror, ks, ke, dis, die, djs, dje
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_r
      !
      call get_nx(info, nx); call get_ny(info, ny)
      call get_di(info, dis, die); call get_dj(info, djs, dje)
      call get_nx0(info, nx0); call get_ny0(info, ny0); call get_ne0(info, ne0)
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype); call get_nepe(info, nepe)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      nvar = self(1)%ncontrol
      ngroup = nt*ne/nproc + 1
      if (mod(nt*ne,nproc) == 0) ngroup = ngroup - 1
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self(1), ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum_r(nx0,ny0,nz,nvar), datum(nx0,ny0,nz))
      !
      call set_index_nohalo(ne0, nepe, myide, ks, ke)
      do igroup = 1, ngroup
	 itm = (igroup-1)*nproc + myid
	 if (itm < nt*ne) then
	    it = itm/ne + 1
	    imember = mod(itm,ne) + ks
	    write(infile(1:2),'(I2.2)') it
	    write(infile(3:6),'(I4.4)') imember
	    !print*, 'Read ', trim(header)//infile
	    open(90, file=trim(header)//infile, form='unformatted')
	    do ivar = 1, nvar
	       read(90) datum_r(:,:,1:nlev(ivar),ivar)
	    end do
	    close(90)
	 end if
	 !
	 do id = 1, nproc
	    itm = (igroup-1)*nproc + id-1
	    if (itm >= nt*ne) cycle
	    it = itm/ne + 1
	    imember = mod(itm,ne) + 1
	    do ivar = 1, nvar
	       if (myid == id-1) datum(:,:,1:nlev(ivar)) = datum_r(:,:,1:nlev(ivar),ivar)
	       call scatter(info, 'xy', datum(:,:,1:nlev(ivar)), id-1, nx0, ny0, nlev(ivar))
	       call set_field(self(imember), ivar, it, k2ijt, datum(1:nx,1:ny,1:nlev(ivar)), nx, ny, nlev(ivar))
	       if (it == nt) then
		  call set_field(xpert(imember), ivar, 2, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)), nlev(ivar))
	       else
		  if (it == itout) call set_field(xpert(imember), ivar, 1, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)), nlev(ivar))
	       end if
	    end do
	 end do
      end do
      deallocate(nlev, datum_r, datum)
      !
      return
   end subroutine read_ensemble_NodeProfileControl2
   !
   !
   !
end module NodeProfileControl_class
   
