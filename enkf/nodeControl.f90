module NodeControl_class
! Author: Le Duc
! Created date: 27 Jul 2014
   use variable, only : nsoil, enkf0, enkf1, enkf2, enkf3, enkf4, enkf5, &
                      & r_size, r_sngl, weight_clim, weight_ens
   use NodeInfo_class
   use NodeField_class
   use NodeMPI
   implicit none
   !
   type NodeFieldPointer
      type(NodeField), pointer :: p
   end type NodeFieldPointer
   !
   type NodeControl
      logical :: nohalo
      integer :: ncontrol
      integer :: nx, ny, nlev, nt
      type(NodeField), pointer :: u, v, w, t, p, qv, qc, qi, qr, qs, qg, tsoil, &
                                & logp, rh, pwv, z, g, uvtqp, uvtq, uvt, uv, ps, logps, pmsl, tgrd, rain
      type(NodeFieldPointer), dimension(:), allocatable :: control
      integer, dimension(:), allocatable :: nvarnhm
      integer, dimension(:,:), allocatable :: ivarnhm, nlevnhm
   end type NodeControl
   !
   integer :: iunhm, ivnhm, iwnhm, itnhm, ipnhm, iqvnhm
   integer :: iqcnhm, iqinhm, iqrnhm, iqsnhm, iqgnhm, itsoilnhm
   !
   interface new
      module procedure new_NodeControl
   end interface
   interface destroy
      module procedure destroy_NodeControl
   end interface
   interface display
      module procedure display_NodeControl
   end interface
   interface assignment(=)
      module procedure copy_NodeControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeControl
   end interface
   interface get_nx
      module procedure get_nx_NodeControl
   end interface
   interface get_ny
      module procedure get_ny_NodeControl
   end interface
   interface get_nt
      module procedure get_nt_NodeControl1
      module procedure get_nt_NodeControl2
   end interface
   interface get_ndim
      module procedure get_ndim_NodeControl
   end interface
   interface get_name
      module procedure get_name_NodeControl
   end interface
   interface get_nlev
      module procedure get_nlev_NodeControl1
      module procedure get_nlev_NodeControl2
      module procedure get_nlev_NodeControl3
   end interface
   interface set_const
      module procedure set_const_NodeControl
   end interface
   interface get_field
      module procedure get_field_NodeControl1
      module procedure get_field_NodeControl2
      module procedure get_field_NodeControl3
      module procedure get_field_NodeControl4
      module procedure get_field_NodeControl5
      module procedure get_field_NodeControl6
      module procedure get_field_NodeControl7
   end interface
   interface get_control
      module procedure get_control_NodeControl
   end interface
   interface set_field
      module procedure set_field_NodeControl1
      module procedure set_field_NodeControl2
      module procedure set_field_NodeControl3
      module procedure set_field_NodeControl4
      module procedure set_field_NodeControl5
      module procedure set_field_NodeControl6
      module procedure set_field_NodeControl7
      module procedure set_field_NodeControl8
   end interface
   interface set_control
      module procedure set_control_NodeControl
   end interface
   interface add_field
      module procedure add_field_NodeControl1
      module procedure add_field_NodeControl2
      module procedure add_field_NodeControl3
      module procedure add_field_NodeControl4
   end interface
   interface add
      module procedure add_NodeControl1
      module procedure add_NodeControl2
      module procedure add_NodeControl3
   end interface
   interface subtract
      module procedure subtract_NodeControl1
      module procedure subtract_NodeControl2
   end interface
   interface power
      module procedure power_NodeControl
   end interface
   interface ratransform
      module procedure ratransform_NodeControl
   end interface
   interface multiply
      module procedure multiply_NodeControl1
      module procedure multiply_NodeControl2
      module procedure multiply_NodeControl3
      module procedure multiply_NodeControl4
   end interface
   interface divide
      module procedure divide_NodeControl1
      module procedure divide_NodeControl2
      module procedure divide_NodeControl3
      module procedure divide_NodeControl4
   end interface
   interface compute_normsquare
      module procedure compute_normsquare_NodeControl
   end interface
   interface compute_moistnorm
      module procedure compute_moistnorm_NodeControl
   end interface
   interface innerproduct
      module procedure innerproduct_NodeControl1
      module procedure innerproduct_NodeControl2
   end interface
   interface localize
      module procedure localize_NodeControl1
      module procedure localize_NodeControl2
   end interface
   interface randomize
      module procedure randomize_NodeControl1
      module procedure randomize_NodeControl2
   end interface
   interface broadcast_bck
      module procedure broadcast_bck_NodeControl
   end interface
   interface allreduce_ens
      module procedure allreduce_ens_NodeControl
   end interface
   interface read_control
      module procedure read_control_NodeControl
   end interface
   interface read_ensemble
      module procedure read_ensemble_NodeControl
   end interface
   interface write_control
      module procedure write_control_NodeControl1
      module procedure write_control_NodeControl2
   end interface
   interface write_ensemble
      module procedure write_ensemble_NodeControl
   end interface
   interface update
      module procedure update_NodeControl
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_nhm_NodeControl(self, info, nohalo, nlev, nt)
   ! the surface fields are imposed at the bottom level.
   ! the p field is non-hydrostatic pressure
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 12
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      iunhm = ivar
      ! v
      allocate(self%v)
      call new(self%v, 'v', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ivnhm = ivar
      ! w
      allocate(self%w)
      call new(self%w, 'w', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%w
      iwnhm = ivar
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      itnhm = ivar
      ! p
      allocate(self%p)
      call new(self%p, 'p', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ipnhm = ivar
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      iqvnhm = ivar
      ! qc
      allocate(self%qc)
      call new(self%qc, 'qc', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qc
      iqcnhm = ivar
      ! qi
      allocate(self%qi)
      call new(self%qi, 'qi', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qi
      iqinhm = ivar
      ! qr
      allocate(self%qr)
      call new(self%qr, 'qr', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qr
      iqrnhm = ivar
      ! qs
      allocate(self%qs)
      call new(self%qs, 'qs', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qs
      iqsnhm = ivar
      ! qg
      allocate(self%qg)
      call new(self%qg, 'qg', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qg
      iqgnhm = ivar
      ! tsoil
      allocate(self%tsoil)
      call new(self%tsoil, 'tsoil', nx, ny, nsoil, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tsoil
      itsoilnhm = ivar
      !
      return
   end subroutine new_nhm_NodeControl
   !
   !
   !
   subroutine new_obs_NodeControl(self, info, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      integer :: nx, ny, ivar
      !
      self%ncontrol = 14
      self%nohalo = .False.
      call get_nx(info, nx); call get_ny(info, ny)
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! logp
      allocate(self%logp)
      call new(self%logp, 'logp', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%logp
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      ! rh
      allocate(self%rh)
      call new(self%rh, 'rh', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rh
      ! pwv
      allocate(self%pwv)
      call new(self%pwv, 'pwv', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pwv
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      ! logps
      allocate(self%logps)
      call new(self%logps, 'logps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%logps
      ! pmsl
      allocate(self%pmsl)
      call new(self%pmsl, 'pmsl', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pmsl
      ! rain
      allocate(self%rain)
      call new(self%rain, 'rain', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rain
      ! z
      allocate(self%z)
      call new(self%z, 'z', nx, ny, nlev, 1)
      ivar = ivar + 1
      self%control(ivar)%p => self%z
      ! g
      allocate(self%g)
      call new(self%g, 'g', nx, ny, nlev, 1)
      ivar = ivar + 1
      self%control(ivar)%p => self%g
      !
      return
   end subroutine new_obs_NodeControl
   !
   !
   !
   subroutine new_sens0_NodeControl(self, info, nlev, nt, nhm)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      type(NodeControl), intent(in) :: nhm
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = .True.
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      self%ps%field(:,:,1,:) = nhm%p%field(:,:,1,:)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      !
      return
   end subroutine new_sens0_NodeControl
   !
   !
   !
   subroutine new_sens1_NodeControl(self, info, nlev, nt, nhm)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      type(NodeControl), intent(in) :: nhm
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 5
      self%nohalo = .True.
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nx, ny, nlev, nt)
      self%u = nhm%u
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', nx, ny, nlev, nt)
      self%v = nhm%v
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      self%t = nhm%t
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      self%qv = nhm%qv
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      self%ps%field(:,:,1,:) = nhm%p%field(:,:,1,:)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      !
      return
   end subroutine new_sens1_NodeControl
   !
   !
   !
   subroutine new_logp_NodeControl(self, info, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      integer :: nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = .False.
      call get_nx(info, nx); call get_ny(info, ny)
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! logp
      allocate(self%logp)
      call new(self%logp, 'logp', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%logp
      !
      return
   end subroutine new_logp_NodeControl
   !
   !
   !
   subroutine new_tqv_NodeControl(self, info, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      integer :: nx, ny, ivar
      !
      self%ncontrol = 2
      self%nohalo = .False.
      call get_nx(info, nx); call get_ny(info, ny)
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      !
      return
   end subroutine new_tqv_NodeControl
   !
   !
   !
   subroutine new_rain_NodeControl(self, info, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      integer :: nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = .False.
      call get_nx(info, nx); call get_ny(info, ny)
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! rain
      allocate(self%rain)
      call new(self%rain, 'rain', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rain
      !
      return
   end subroutine new_rain_NodeControl
   !
   !
   !
   subroutine new_pmsl_NodeControl(self, info, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nlev, nt
      integer :: nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = .False.
      call get_nx(info, nx); call get_ny(info, ny)
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! pmsl
      allocate(self%pmsl)
      call new(self%pmsl, 'pmsl', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pmsl
      !
      return
   end subroutine new_pmsl_NodeControl
   !
   !
   !
   subroutine new_enkf0_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = 1; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uvtq
      allocate(self%uvtqp)
      call new(self%uvtqp, 'uvtqp', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uvtqp
      self%nvarnhm(ivar) = 6
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      self%ivarnhm(ivar,3) = itnhm;  self%nlevnhm(ivar,3) = nlevnhm
      self%ivarnhm(ivar,4) = iqvnhm; self%nlevnhm(ivar,4) = nlevnhm
      self%ivarnhm(ivar,5) = ipnhm; self%nlevnhm(ivar,5) = 1
      self%ivarnhm(ivar,6) = itsoilnhm; self%nlevnhm(ivar,6) = 1
      !
      return
   end subroutine new_enkf0_NodeControl
   !
   !
   !
   subroutine new_sens0old_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 1
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = 1; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uvtq
      allocate(self%uvtqp)
      call new(self%uvtqp, 'uvtqp', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uvtqp
      self%nvarnhm(ivar) = 12
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      self%ivarnhm(ivar,3) = iwnhm;  self%nlevnhm(ivar,3) = nlevnhm
      self%ivarnhm(ivar,4) = itnhm;  self%nlevnhm(ivar,4) = nlevnhm
      self%ivarnhm(ivar,5) = ipnhm; self%nlevnhm(ivar,5) = nlevnhm
      self%ivarnhm(ivar,6) = iqvnhm; self%nlevnhm(ivar,6) = nlevnhm
      self%ivarnhm(ivar,7) = iqcnhm; self%nlevnhm(ivar,7) = nlevnhm
      self%ivarnhm(ivar,8) = iqinhm; self%nlevnhm(ivar,8) = nlevnhm
      self%ivarnhm(ivar,9) = iqrnhm; self%nlevnhm(ivar,9) = nlevnhm
      self%ivarnhm(ivar,10) = iqsnhm; self%nlevnhm(ivar,10) = nlevnhm
      self%ivarnhm(ivar,11) = iqgnhm; self%nlevnhm(ivar,11) = nlevnhm
      self%ivarnhm(ivar,12) = itsoilnhm; self%nlevnhm(ivar,12) = nsoil
      !
      return
   end subroutine new_sens0old_NodeControl
   !
   !
   !
   subroutine new_enkf1_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 2
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uvtq
      allocate(self%uvtqp)
      call new(self%uvtqp, 'uvtqp', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uvtqp
      self%nvarnhm(ivar) = 5
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      self%ivarnhm(ivar,3) = itnhm;  self%nlevnhm(ivar,3) = nlevnhm
      self%ivarnhm(ivar,4) = iqvnhm; self%nlevnhm(ivar,4) = nlevnhm
      self%ivarnhm(ivar,5) = ipnhm;  self%nlevnhm(ivar,5) = 1
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itsoilnhm; self%nlevnhm(ivar,1) = 1
      !
      return
   end subroutine new_enkf1_NodeControl
   !
   !
   !
   subroutine new_enkf2_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 3
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uvtq
      allocate(self%uvtq)
      call new(self%uvtq, 'uvtq', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uvtq
      self%nvarnhm(ivar) = 4
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      self%ivarnhm(ivar,3) = itnhm;  self%nlevnhm(ivar,3) = nlevnhm
      self%ivarnhm(ivar,4) = iqvnhm; self%nlevnhm(ivar,4) = nlevnhm
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = ipnhm; self%nlevnhm(ivar,1) = 1
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itsoilnhm; self%nlevnhm(ivar,1) = 1
      !
      return
   end subroutine new_enkf2_NodeControl
   !
   !
   !
   subroutine new_enkf3_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 4
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uvt
      allocate(self%uvt)
      call new(self%uvt, 'uvt', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uvt
      self%nvarnhm(ivar) = 3
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      self%ivarnhm(ivar,3) = itnhm;  self%nlevnhm(ivar,3) = nlevnhm
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = iqvnhm; self%nlevnhm(ivar,1) = nlevnhm
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = ipnhm; self%nlevnhm(ivar,1) = 1
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itsoilnhm; self%nlevnhm(ivar,1) = 1
      !
      return
   end subroutine new_enkf3_NodeControl
   !
   !
   !
   subroutine new_enkf4_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 5
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! uv
      allocate(self%uv)
      call new(self%uv, 'uv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%uv
      self%nvarnhm(ivar) = 2
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      self%ivarnhm(ivar,2) = ivnhm;  self%nlevnhm(ivar,2) = nlevnhm
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itnhm;  self%nlevnhm(ivar,1) = nlevnhm
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = iqvnhm; self%nlevnhm(ivar,1) = nlevnhm
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = ipnhm; self%nlevnhm(ivar,1) = 1
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itsoilnhm; self%nlevnhm(ivar,1) = 1
      !
      return
   end subroutine new_enkf4_NodeControl
   !
   !
   !
   subroutine new_enkf5_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: nlevnhm, nlev, nt
      integer :: dis, die, djs, dje, nx, ny, ivar
      !
      self%ncontrol = 6
      self%nohalo = nohalo
      call get_nx(info, nx); call get_ny(info, ny)
      if (self%nohalo) then
	 call get_di(info, dis, die); call get_dj(info, djs, dje)
	 nx = nx - dis - die; ny = ny - djs - dje
      end if
      self%nx = nx; self%ny = ny
      self%nlev = nlev; self%nt = nt
      allocate(self%control(self%ncontrol))
      allocate(self%nvarnhm(self%ncontrol))
      allocate(self%ivarnhm(self%ncontrol,20))
      allocate(self%nlevnhm(self%ncontrol,20))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = iunhm;  self%nlevnhm(ivar,1) = nlevnhm
      ! v
      allocate(self%v)
      call new(self%v, 'v', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = ivnhm;  self%nlevnhm(ivar,1) = nlevnhm
      ! t
      allocate(self%t)
      call new(self%t, 't', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itnhm;  self%nlevnhm(ivar,1) = nlevnhm
      ! qv
      allocate(self%qv)
      call new(self%qv, 'qv', nx, ny, nlev, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%qv
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = iqvnhm; self%nlevnhm(ivar,1) = nlevnhm
      ! ps
      allocate(self%ps)
      call new(self%ps, 'ps', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%ps
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = ipnhm; self%nlevnhm(ivar,1) = 1
      ! tgrd
      allocate(self%tgrd)
      call new(self%tgrd, 'tgrd', nx, ny, 1, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%tgrd
      self%nvarnhm(ivar) = 1
      self%ivarnhm(ivar,1) = itsoilnhm; self%nlevnhm(ivar,1) = 1
      !
      return
   end subroutine new_enkf5_NodeControl
   !
   !
   !
   subroutine new_NodeControl(self, control_mode, info, nohalo, nlevnhm, nlev, nt)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      logical, intent(in) :: nohalo
      integer, intent(in) :: control_mode, nlevnhm, nlev, nt
      !
      if (control_mode == enkf0) then
	 call new_enkf0_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      else if (control_mode == enkf1) then
	 call new_enkf1_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      else if (control_mode == enkf2) then
	 call new_enkf2_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      else if (control_mode == enkf3) then
	 call new_enkf3_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      else if (control_mode == enkf4) then
	 call new_enkf4_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      else if (control_mode == enkf5) then
	 call new_enkf5_NodeControl(self, info, nohalo, nlevnhm, nlev, nt)
      end if
      !
      return
   end subroutine new_NodeControl
   !
   !
   !
   subroutine destroy_NodeControl(self)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      if (allocated(self%nvarnhm)) deallocate(self%nvarnhm)
      if (allocated(self%ivarnhm)) deallocate(self%ivarnhm)
      !
      return
   end subroutine destroy_NodeControl
   !
   !
   !
   subroutine display_NodeControl(self)
      implicit none
      type(NodeControl), intent(in) :: self
      !
      print*, 'Control: ', self%ncontrol
      print*, 'Dimension: ', self%nx, self%ny, self%nlev, self%nt
      !
      return
   end subroutine display_NodeControl
   !
   !
   !
   subroutine get_ncontrol_NodeControl(self, ncontrol)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeControl
   !
   !
   !
   subroutine get_nx_NodeControl(self, nx)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(out) :: nx
      !
      nx = self%nx
      !
      return
   end subroutine get_nx_NodeControl
   !
   !
   !
   subroutine get_ny_NodeControl(self, ny)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(out) :: ny
      !
      ny = self%ny
      !
      return
   end subroutine get_ny_NodeControl
   !
   !
   !
   subroutine get_nt_NodeControl1(self, ivar, nt)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nt
      !
      call get_nt(self%control(ivar)%p, nt)
      !
      return
   end subroutine get_nt_NodeControl1
   !
   !
   !
   subroutine get_nt_NodeControl2(self, nt)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(out) :: nt
      !
      nt = self%nt
      !
      return
   end subroutine get_nt_NodeControl2
   !
   !
   !
   subroutine get_ndim_NodeControl(self, ndim)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(out) :: ndim
      integer :: ivar, n
      !
      ndim = 0
      do ivar = 1, self%ncontrol
	 call get_ndim(self%control(ivar)%p, n)
	 ndim = ndim + n
      end do
      !
      return
   end subroutine get_ndim_NodeControl
   !
   !
   !
   subroutine get_name_NodeControl(self, ivar, name)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeControl
   !
   !
   !
   subroutine get_nlev_NodeControl1(self, ivar, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nlev
      !
      call get_nlev(self%control(ivar)%p, nlev)
      !
      return
   end subroutine get_nlev_NodeControl1
   !
   !
   !
   subroutine get_nlev_NodeControl2(self, name, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
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
   end subroutine get_nlev_NodeControl2
   !
   !
   !
   subroutine get_nlev_NodeControl3(self, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
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
   end subroutine get_nlev_NodeControl3
   !
   !
   !
   subroutine copy_NodeControl(self, f)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: f
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p = f%control(ivar)%p
      end do
      !
      return
   end subroutine copy_NodeControl
   !
   !
   !
   subroutine set_const_NodeControl(self, const)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar, nt, it
      !
      do ivar = 1, self%ncontrol
	 call get_nt(self, ivar, nt)
         do it = 1, nt
	    call set_field(self%control(ivar)%p, it, const)
	 end do
      end do
      !
      return
   end subroutine set_const_NodeControl
   !
   !
   !
   subroutine get_field_NodeControl1(self, ivar, it, field, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar, it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(inout) :: field
      !
      call get_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine get_field_NodeControl1
   !
   !
   !
   subroutine get_field_NodeControl2(self, name, it, field, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(inout) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine get_field_NodeControl2
   !
   !
   !
   subroutine get_field_NodeControl3(self, ivar, field, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar, nlev
      real(r_size), dimension(self%nx,self%ny,nlev,self%nt), intent(out) :: field
      !
      call get_field(self%control(ivar)%p, field)
      !
      return
   end subroutine get_field_NodeControl3
   !
   !
   !
   subroutine get_field_NodeControl4(self, ivar, field)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: ivar
      type(NodeField), intent(out) :: field
      !
      field = self%control(ivar)%p
      !
      return
   end subroutine get_field_NodeControl4
   !
   !
   !
   subroutine get_field_NodeControl5(self, itout, object)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: itout
      type(NodeControl), intent(inout) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_field(self%control(ivar)%p, itout, object%control(ivar)%p)
      end do
      !
      return
   end subroutine get_field_NodeControl5
   !
   !
   !
   subroutine get_field_NodeControl6(self, field, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: nlev
      real(r_size), dimension(self%nx,self%ny,nlev,self%nt), intent(inout) :: field
      integer :: ivar, nz, k
      !
      k = 1
      do ivar = 1, self%ncontrol
	 call get_nlev(self%control(ivar)%p, nz)
	 call get_field(self%control(ivar)%p, field(:,:,k:k+nz-1,:))
	 k = k + nz
      end do
      !
      return
   end subroutine get_field_NodeControl6
   !
   !
   !
   subroutine get_field_NodeControl7(self, name, is, ie, js, je, it, field, nlev)
      implicit none
      type(NodeControl), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: is, ie, js, je, it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(inout) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call get_field(self%control(ivar)%p, is, ie, js, je, it, field)
      !
      return
   end subroutine get_field_NodeControl7
   !
   !
   !
   subroutine get_control_NodeControl(self, it, w)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: it
      type(NodeControl), intent(inout) :: w
      integer :: nlev, ivar, iweight, ivarnhm
      !
      do iweight = 1, w%ncontrol
         call get_nlev(w, iweight, nlev)
         !if (nlev > 1) nlev = nlev - 1
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            w%control(iweight)%p%field(:,:,1:nlev,1) = self%control(ivarnhm)%p%field(:,:,1:nlev,it)
         end do
      end do
      !
      return
   end subroutine get_control_NodeControl
   !
   !
   !
   subroutine set_field_NodeControl1(self, ivar, it, field)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar, it
      real(r_size), intent(in) :: field
      !
      call set_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine set_field_NodeControl1
   !
   !
   !
   subroutine set_field_NodeControl2(self, ivar, it, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar, it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(in) :: field
      !
      call set_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine set_field_NodeControl2
   !
   !
   !
   subroutine set_field_NodeControl3(self, name, it, field)
      implicit none
      type(NodeControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it
      real(r_size), intent(in) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call set_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine set_field_NodeControl3
   !
   !
   !
   subroutine set_field_NodeControl4(self, name, it, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(in) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call set_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine set_field_NodeControl4
   !
   !
   !
   subroutine set_field_NodeControl5(self, ivar, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar, nlev
      real(r_size), dimension(self%nx,self%ny,nlev,self%nt), intent(in) :: field
      !
      call set_field(self%control(ivar)%p, field)
      !
      return
   end subroutine set_field_NodeControl5
   !
   !
   !
   subroutine set_field_NodeControl6(self, ivar, field)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar
      type(NodeField), intent(in) :: field
      !
      self%control(ivar)%p = field
      !
      return
   end subroutine set_field_NodeControl6
   !
   !
   !
   subroutine set_field_NodeControl7(self, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: nlev
      real(r_size), dimension(self%nx,self%ny,nlev,self%nt), intent(in) :: field
      integer :: ivar, nz, k
      !
      k = 1
      do ivar = 1, self%ncontrol
	 call get_nlev(self%control(ivar)%p, nz)
	 call set_field(self%control(ivar)%p, field(:,:,k:k+nz-1,:))
	 k = k + nz
      end do
      !
      return
   end subroutine set_field_NodeControl7
   !
   !
   !
   subroutine set_field_NodeControl8(self, name, is, ie, js, je, it, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: is, ie, js, je, it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(in) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call set_field(self%control(ivar)%p, is, ie, js, je, it, field)
      !
      return
   end subroutine set_field_NodeControl8
   !
   !
   !
   subroutine set_control_NodeControl(self, it, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: it
      type(NodeControl), intent(in) :: w
      integer :: nlev, ivar, iweight, ivarnhm
      !
      do iweight = 1, w%ncontrol
         call get_nlev(w, iweight, nlev)
         !if (nlev > 1) nlev = nlev - 1
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            self%control(ivarnhm)%p%field(:,:,1:nlev,it) = w%control(iweight)%p%field(:,:,1:nlev,1)
         end do
      end do
      !
      return
   end subroutine set_control_NodeControl
   !
   !
   !
   subroutine add_field_NodeControl1(self, ivar, it, field)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar, it
      real(r_size), intent(in) :: field
      !
      call add_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine add_field_NodeControl1
   !
   !
   !
   subroutine add_field_NodeControl2(self, ivar, it, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ivar, it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(in) :: field
      !
      call add_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine add_field_NodeControl2
   !
   !
   !
   subroutine add_field_NodeControl3(self, name, it, field)
      implicit none
      type(NodeControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it
      real(r_size), intent(in) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call add_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine add_field_NodeControl3
   !
   !
   !
   subroutine add_field_NodeControl4(self, name, it, field, nlev)
      implicit none
      type(NodeControl), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: it, nlev
      real(r_size), dimension(self%nx,self%ny,nlev), intent(in) :: field
      character(len=10) :: varname
      integer :: nvar, ivar
      !
      nvar = self%ncontrol
      do ivar = 1, nvar
         call get_name(self%control(ivar)%p, varname)
	 if (trim(varname) == trim(name)) exit
      end do
      call add_field(self%control(ivar)%p, it, field)
      !
      return
   end subroutine add_field_NodeControl4
   !
   !
   !
   subroutine add_NodeControl1(self, const)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine add_NodeControl1
   !
   !
   !
   subroutine add_NodeControl2(self, object)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine add_NodeControl2
   !
   !
   !
   subroutine add_NodeControl3(self, object, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call add(self%control(ivarnhm)%p, object%control(ivarnhm)%p)
	 end do
      end do
      !
      return
   end subroutine add_NodeControl3
   !
   !
   !
   subroutine subtract_NodeControl1(self, object)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine subtract_NodeControl1
   !
   !
   !
   subroutine subtract_NodeControl2(self, object, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call subtract(self%control(ivarnhm)%p, object%control(ivarnhm)%p)
	 end do
      end do
      !
      return
   end subroutine subtract_NodeControl2
   !
   !
   !
   subroutine power_NodeControl(self, const)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call power(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine power_NodeControl
   !
   !
   !
   subroutine ratransform_NodeControl(self, threshold)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: threshold
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call ratransform(self%control(ivar)%p, threshold)
      end do
      !
      return
   end subroutine ratransform_NodeControl
   !
   !
   !
   subroutine multiply_NodeControl1(self, const)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine multiply_NodeControl1
   !
   !
   !
   subroutine multiply_NodeControl2(self, object)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine multiply_NodeControl2
   !
   !
   !
   subroutine multiply_NodeControl3(self, const, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call multiply(self%control(ivarnhm)%p, const)
	 end do
      end do
      !
      return
   end subroutine multiply_NodeControl3
   !
   !
   !
   subroutine multiply_NodeControl4(self, object, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call multiply(self%control(ivarnhm)%p, object%control(ivarnhm)%p)
	 end do
      end do
      !
      return
   end subroutine multiply_NodeControl4
   !
   !
   !
   subroutine divide_NodeControl1(self, const)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine divide_NodeControl1
   !
   !
   !
   subroutine divide_NodeControl2(self, object)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine divide_NodeControl2
   !
   !
   !
   subroutine divide_NodeControl3(self, const, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call divide(self%control(ivarnhm)%p, const)
	 end do
      end do
      !
      return
   end subroutine divide_NodeControl3
   !
   !
   !
   subroutine divide_NodeControl4(self, object, w)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: object
      type(NodeControl), intent(in) :: w
      integer :: iweight, ivar, ivarnhm
      !
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            call divide(self%control(ivarnhm)%p, object%control(ivar)%p)
	 end do
      end do
      !
      return
   end subroutine divide_NodeControl4
   !
   !
   !
   subroutine compute_normsquare_NodeControl(self, info, normsquare)
      implicit none
      type(NodeControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      real(r_size), intent(out) :: normsquare
      integer :: ivar
      real(r_size) :: tmp
      !
      normsquare = 0.d0
      do ivar = 1, self%ncontrol
	 call compute_normsquare(self%control(ivar)%p, info, tmp)
	 normsquare = normsquare + tmp
      end do
      !
      return
   end subroutine compute_normsquare_NodeControl
   !
   !
   !
   subroutine compute_moistnorm_NodeControl(self, info, klev, tref, psref, wt, wqv, normsquare)
      use variable, only : rd, cp, lqv
      implicit none
      type(NodeControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: klev
      real(r_size), intent(in) :: tref, psref, wt, wqv
      real(r_size), intent(out) :: normsquare
      integer :: ivar
      real(r_size) :: tmp
      !
      normsquare = 0.d0
      do ivar = 1, self%ncontrol
         if (ivar == iunhm .or. ivar == ivnhm .or. ivar == itnhm .or. ivar == ipnhm .or. ivar == iqvnhm) then
	    if (ivar == ipnhm) then
	       call compute_normsquare(self%control(ivar)%p, info, 1, tmp)
	    else
	       call compute_normsquare(self%control(ivar)%p, info, klev, tmp)
	    end if
	    if (ivar == iunhm .or. ivar == ivnhm) then
	       normsquare = normsquare + tmp
	    else if (ivar == itnhm) then
	       normsquare = normsquare + wt*cp/tref*tmp
	    else if (ivar == ipnhm) then
	       normsquare = normsquare + rd*tref/psref**2*tmp
	    else if (ivar == iqvnhm) then
	       normsquare = normsquare + wqv*lqv**2/(cp*tref)*tmp
	    end if
	 end if
      end do
      !
      return
   end subroutine compute_moistnorm_NodeControl
   !
   !
   !
   subroutine innerproduct_NodeControl1(self, object, info, product)
      implicit none
      type(NodeControl), intent(in) :: self
      type(NodeControl), intent(in) :: object
      type(NodeInfo), intent(in) :: info
      real(r_size), intent(out) :: product
      integer :: ivar
      real(r_size) :: tmp
      !
      product = 0.d0
      do ivar = 1, self%ncontrol
	 call innerproduct(self%control(ivar)%p, object%control(ivar)%p, info, tmp)
	 product = product + tmp
      end do
      !
      return
   end subroutine innerproduct_NodeControl1
   !
   !
   !
   subroutine innerproduct_NodeControl2(self, object, info, w, product)
      implicit none
      type(NodeControl), intent(in) :: self
      type(NodeControl), intent(in) :: object
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: w
      real(r_size), intent(out) :: product
      integer :: nlevnhm, iweight, ivar, ivarnhm
      real(r_size) :: tmp
      !
      product = 0.d0
      do iweight = 1, w%ncontrol
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
	    nlevnhm = w%nlevnhm(iweight,ivar)
            call innerproduct(self%control(ivarnhm)%p, object%control(ivarnhm)%p, info, nlevnhm, tmp)
	    product = product + tmp
	 end do
      end do
      !
      return
   end subroutine innerproduct_NodeControl2
   !
   !
   !
   subroutine localize_NodeControl1(self, info, i0, j0)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: i0, j0
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call localize(self%control(ivar)%p, info, i0, j0)
      end do
      !
      return
   end subroutine localize_NodeControl1
   !
   !
   !
   subroutine localize_NodeControl2(self, info, i0, j0, R0)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: i0, j0
      real(r_size), intent(in) :: R0
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call localize(self%control(ivar)%p, info, i0, j0, R0)
      end do
      !
      return
   end subroutine localize_NodeControl2
   !
   !
   !
   subroutine randomize_NodeControl1(self, info)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call randomize(self%control(ivar)%p, info)
      end do
      !
      return
   end subroutine randomize_NodeControl1
   !
   !
   !
   subroutine randomize_NodeControl2(self, info, it0)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: it0
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call randomize(self%control(ivar)%p, info, it0)
      end do
      !
      return
   end subroutine randomize_NodeControl2
   !
   !
   !
   subroutine broadcast_bck_NodeControl(self, info, source_ide)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: source_ide
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call broadcast_bck(self%control(ivar)%p, info, source_ide)
      end do
      !
      return
   end subroutine broadcast_bck_NodeControl
   !
   !
   !
   subroutine allreduce_ens_NodeControl(self)
      implicit none
      type(NodeControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call allreduce_ens(self%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ens_NodeControl
   !
   !
   !
   subroutine read_control_NodeControl(self, info, header)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: header
      !
      character(2) :: infile = '00'
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, nx0, ny0, nx, ny, nz, nt, nvar, ngroup
      integer :: igroup, it, ivar, id, ierror
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_r
      !
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype); call get_nepe(info, nepe)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      call get_nx0(info, nx0); call get_ny0(info, ny0)
      nx = self%nx; ny = self%ny
      nt = self%nt
      nvar = self%ncontrol
      ngroup = nt/nproc + 1
      if (mod(nt,nproc) == 0) ngroup = ngroup - 1
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self, ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum(nx0,ny0,nz))
      if (myide == 0) allocate(datum_r(nx0,ny0,nz,nvar))
      !
      do igroup = 1, ngroup
         if (myide == 0) then
	    it = (igroup-1)*nproc + myid+1
	    if (it <= nt) then
	       write(infile(1:2),'(I2.2)') it
	       !print*, 'Read ', trim(header)//infile
	       open(90, file=trim(header)//infile, form='unformatted')
	       do ivar = 1, nvar
	          read(90) datum_r(:,:,1:nlev(ivar),ivar)
	       end do
	       close(90)
	    end if
	 end if
	 !
	 do id = 1, nproc
	    it = (igroup-1)*nproc + id
	    if (it > nt) cycle
	    do ivar = 1, nvar
	       if (myide == 0) then
	          if (myid == id-1) datum(:,:,1:nlev(ivar)) = datum_r(:,:,1:nlev(ivar),ivar)
	          call scatter(info, 'xy', datum(:,:,1:nlev(ivar)), id-1, nx0, ny0, nlev(ivar))
	       end if
	       call broadcast3D('e', datum(1:nx,1:ny,1:nlev(ivar)), 0, nx, ny, nlev(ivar))
	       call set_field(self, ivar, it, datum(1:nx,1:ny,1:nlev(ivar)), nlev(ivar))
	    end do
	 end do
      end do
      deallocate(nlev, datum)
      if (myide == 0) deallocate(datum_r)
      !
      return
   end subroutine read_control_NodeControl
   !
   !
   !
   subroutine read_ensemble_NodeControl(self, info, header, ne)
      implicit none
      integer, intent(in) :: ne
      type(NodeControl), dimension(ne), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: header
      !
      character(6) :: infile = '000000'
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, nx0, ny0, ne0, nx, ny, nz, nt, nvar, ngroup
      integer :: igroup, itm, it, imember, ivar, id, ide, ierror, ks, ke
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_r
      !
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype); call get_nepe(info, nepe)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      call get_ne0(info, ne0)
      nx = self(1)%nx
      ny = self(1)%ny
      nt = self(1)%nt
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
	       call set_field(self(imember), ivar, it, datum(1:nx,1:ny,1:nlev(ivar)), nlev(ivar))
	    end do
	 end do
      end do
      deallocate(nlev, datum_r, datum)
      !
      return
   end subroutine read_ensemble_NodeControl
   !
   !
   !
   subroutine write_control_NodeControl1(self, it, info, header)
      implicit none
      type(NodeControl), intent(in) :: self
      integer, intent(in) :: it
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: header
      !
      character(2) :: outfile
      integer :: myid, myidx, myidy, myide, nxpe, nype, nproc, nx0, ny0, nx, ny, nz, nvar
      integer :: id, ivar, dis, die, djs, dje
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_w
      !
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      if (myide > 0) return
      call get_nxpe(info, nxpe); call get_nype(info, nype)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      call get_nx0(info, nx0); call get_ny0(info, ny0)
      call get_nx(info, nx); call get_ny(info, ny)
      nvar = self%ncontrol
      if (self%nohalo) then ! EnSRF
         call get_di(info, dis, die)
         call get_dj(info, djs, dje)
      end if
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self, ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum(nx0,ny0,nz))
      if (myid == 0) allocate(datum_w(nx0,ny0,nz,nvar))
      !
      do ivar = 1, nvar
         if (self%nohalo) then ! EnSRF
	    call get_field(self%control(ivar)%p, it, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)))
	 else
	    call get_field(self%control(ivar)%p, it, datum(1:nx,1:ny,1:nlev(ivar)))
         end if
         call gather(info, 'xy', datum, 0, nx0, ny0, nlev(ivar))
         if (myid == 0) datum_w(:,:,1:nlev(ivar),ivar) = datum(:,:,1:nlev(ivar))
      end do
      !
      if (myid == 0) then
	 write(outfile(1:2),'(I2.2)') it
	 !print*, 'Write ', trim(header)//outfile
	 open(90, file=trim(header)//outfile, form='unformatted')
	 do ivar = 1, nvar
	    write(90) datum_w(1:nx0,1:ny0,1:nlev(ivar),ivar)
	 end do
	 close(90)
	 deallocate(datum_w)
      end if
      deallocate(datum)
      !
      return
   end subroutine write_control_NodeControl1
   !
   !
   !
   subroutine write_control_NodeControl2(self, info, header)
      implicit none
      type(NodeControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: header
      !
      character(2) :: outfile
      integer :: myid, myidx, myidy, myide, nxpe, nype, nproc, nx0, ny0, nx, ny, nz, nt, nvar, ngroup
      integer :: igroup, it, id, ivar, dis, die, djs, dje
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_w
      !
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      if (myide > 0) return
      call get_nxpe(info, nxpe); call get_nype(info, nype)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      call get_nx0(info, nx0); call get_ny0(info, ny0)
      call get_nx(info, nx); call get_ny(info, ny)
      nt = self%nt
      nvar = self%ncontrol
      ngroup = nt/nproc + 1
      if (mod(nt,nproc) == 0) ngroup = ngroup - 1
      if (self%nohalo) then ! EnSRF
         call get_di(info, dis, die)
         call get_dj(info, djs, dje)
      end if
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self, ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum(nx0,ny0,nz))
      if (myid < nt) allocate(datum_w(nx0,ny0,nz,nvar))
      !
      do igroup = 1, ngroup
         do id = 1, nproc
	    it = (igroup-1)*nproc + id
	    if (it > nt) cycle
	    do ivar = 1, nvar
	       if (self%nohalo) then ! EnSRF
	          call get_field(self%control(ivar)%p, it, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)))
	       else
	          call get_field(self%control(ivar)%p, it, datum(1:nx,1:ny,1:nlev(ivar)))
               end if
               call gather(info, 'xy', datum, id-1, nx0, ny0, nlev(ivar))
	       if (myid == id-1) datum_w(:,:,1:nlev(ivar),ivar) = datum(:,:,1:nlev(ivar))
	    end do
	 end do
	 !
	 it = (igroup-1)*nproc + myid+1
	 if (it <= nt) then
	    write(outfile(1:2),'(I2.2)') it
	    !print*, 'Write ', trim(header)//outfile
	    open(90, file=trim(header)//outfile, form='unformatted')
	    do ivar = 1, nvar
	       write(90) datum_w(1:nx0,1:ny0,1:nlev(ivar),ivar)
	    end do
	    close(90)
	 end if
      end do
      deallocate(datum, datum_w)
      !
      return
   end subroutine write_control_NodeControl2
   !
   !
   !
   subroutine write_ensemble_NodeControl(self, it, imember1, imember2, info, header, ne)
      implicit none
      type(NodeControl), dimension(ne), intent(inout) :: self
      integer, intent(in) :: it, imember1, imember2, ne
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: header
      !
      character(6) :: outfile
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, nx0, ny0, ne0, nx, ny, nz, nt, nvar, ngroup
      integer :: igroup, imember, id, ide, ivar, ks, ke, dis, die, djs, dje
      integer, dimension(:), allocatable :: nlev
      real(r_size), dimension(:,:,:), allocatable :: datum
      real(r_sngl), dimension(:,:,:,:), allocatable :: datum_w
      !
      call get_myidx(info, myidx); call get_myidy(info, myidy); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype); call get_nepe(info, nepe)
      myid = myidy*nxpe+myidx
      nproc = nxpe*nype
      call get_nx0(info, nx0); call get_ny0(info, ny0)
      call get_nx(info, nx); call get_ny(info, ny)
      call get_ne0(info, ne0)
      if (self(1)%nohalo) then ! EnSRF
         call get_di(info, dis, die)
         call get_dj(info, djs, dje)
      end if
      nt = self(1)%nt
      nvar = self(1)%ncontrol
      ngroup = ne/nproc + 1
      if (mod(ne,nproc) == 0) ngroup = ngroup - 1
      !
      allocate(nlev(nvar))
      do ivar = 1, nvar
         call get_nlev(self(1), ivar, nlev(ivar))
      end do
      nz = maxval(nlev)
      allocate(datum(nx0,ny0,nz), datum_w(nx0,ny0,nz,nvar))
      !
      call set_index_nohalo(ne0, nepe, myide, ks, ke)
      do igroup = 1, ngroup
         do id = 1, nproc
	    imember = (igroup-1)*nproc + id
	    if (imember > ne .or. imember+ks-1 < imember1 .or. imember+ks-1 > imember2) cycle
	    do ivar = 1, nvar
	       if (self(1)%nohalo) then ! EnSRF
	          call get_field(self(imember)%control(ivar)%p, it, datum(1+dis:nx-die,1+djs:ny-dje,1:nlev(ivar)))
	       else
	          call get_field(self(imember)%control(ivar)%p, it, datum(1:nx,1:ny,1:nlev(ivar)))
               end if
               call gather(info, 'xy', datum, id-1, nx0, ny0, nlev(ivar))
	       if (myid == id-1) datum_w(:,:,1:nlev(ivar),ivar) = datum(:,:,1:nlev(ivar))
	    end do
	 end do
	 !
	 imember = (igroup-1)*nproc + myid + 1
	 if (imember <= ne .and. imember1 <= imember+ks-1 .and. imember+ks-1 <= imember2) then
	    write(outfile(1:2),'(I2.2)') it
	    write(outfile(3:6),'(I4.4)') imember+ks-1
	    !print*, 'Write ', trim(header)//outfile
	    open(90, file=trim(header)//outfile, form='unformatted')
	    do ivar = 1, nvar
	       write(90) datum_w(1:nx0,1:ny0,1:nlev(ivar),ivar)
	    end do
	    close(90)
	 end if
      end do
      deallocate(datum, datum_w)
      !
      return
   end subroutine write_ensemble_NodeControl
   !
   !
   !
   subroutine apply_USIens(self, ne, w, xpert)
   ! no temporal localization
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeControl), dimension(1:ne), intent(in) :: w
      type(NodeControl), dimension(ne), intent(in) :: xpert
      integer :: nlev, nvarnhm, nlevnhm, ivar, iweight, ivarnhm, k, it, ie
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p%field(:,:,:,:) = 0.d0
      end do
      do iweight = 1, w(1)%ncontrol
         call get_nlev(w(1), iweight, nlev)
	 nvarnhm = w(1)%nvarnhm(iweight)
         do ivar = 1, nvarnhm
            ivarnhm = w(1)%ivarnhm(iweight,ivar)
	    nlevnhm = w(1)%nlevnhm(iweight,ivar)
	    if (nvarnhm < 12 .and. nlevnhm > 1) nlevnhm = nlevnhm - 1 ! nvarnhm = 12 is for sensitivity test
	    if (nlev < nlevnhm) then ! no vertical localization
	       do it = 1, self%nt
		  do k = 1, nlevnhm
		     do ie = 1, ne
			self%control(ivarnhm)%p%field(:,:,k,it) = self%control(ivarnhm)%p%field(:,:,k,it) + &
								& w(ie)%control(iweight)%p%field(:,:,1,1)*xpert(ie)%control(ivarnhm)%p%field(:,:,k,it)
		     end do
		     self%control(ivarnhm)%p%field(:,:,k,it) = weight_ens*self%control(ivarnhm)%p%field(:,:,k,it)
		  end do
	       end do
	    else
	       do it = 1, self%nt
		  do ie = 1, ne
		     self%control(ivarnhm)%p%field(:,:,1:nlevnhm,it) = self%control(ivarnhm)%p%field(:,:,1:nlevnhm,it) + &
								     & w(ie)%control(iweight)%p%field(:,:,1:nlevnhm,1)*xpert(ie)%control(ivarnhm)%p%field(:,:,1:nlevnhm,it)
		  end do
		  self%control(ivarnhm)%p%field(:,:,1:nlevnhm,it) = weight_ens*self%control(ivarnhm)%p%field(:,:,1:nlevnhm,it)
	       end do
	    end if
         end do
      end do
      !
      return
   end subroutine apply_USIens
   !
   !
   !
   subroutine apply_USIclim(self, sigma, w)
   ! no temporal localization
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: sigma
      type(NodeControl), intent(in) :: w
      integer :: nlev, ivar, iweight, ivarnhm, it
      !
      do iweight = 1, w%ncontrol
         call get_nlev(w, iweight, nlev)
         if (nlev > 1) nlev = nlev - 1
         do ivar = 1, w%nvarnhm(iweight)
            ivarnhm = w%ivarnhm(iweight,ivar)
            do it = 1, self%nt
               self%control(ivarnhm)%p%field(:,:,1:nlev,it) = self%control(ivarnhm)%p%field(:,:,1:nlev,it) + weight_clim* &
                                                            & w%control(iweight)%p%field(:,:,1:nlev,1)*sigma%control(iweight)%p%field(:,:,1:nlev,1)
            end do
         end do
      end do
      !
      return
   end subroutine apply_USIclim
   !
   !
   !
   subroutine apply_USITens(self, ne, x, xpert)
   ! no temporal localization
      implicit none
      type(NodeControl), dimension(ne), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeControl), intent(in) :: x
      type(NodeControl), dimension(ne), intent(in) :: xpert
      integer :: nlev, nlevnhm, ivar, iweight, ivarnhm, k, it, ie
      !
      do iweight = 1, self(1)%ncontrol
         call get_nlev(self(1), iweight, nlev)
         do ie = 1, ne
            self(ie)%control(iweight)%p%field(:,:,1:nlev,1) = 0.d0
            do ivar = 1, self(1)%nvarnhm(iweight)
               ivarnhm = self(1)%ivarnhm(iweight,ivar)
	       nlevnhm = self(1)%nlevnhm(iweight,ivar)
	       if (nlev < nlevnhm) then ! no vertical localization
		  do it = 1, x%nt
		     do k = 1, nlevnhm
			self(ie)%control(iweight)%p%field(:,:,1,1) = self(ie)%control(iweight)%p%field(:,:,1,1) + &
								   & x%control(ivarnhm)%p%field(:,:,k,it)*xpert(ie)%control(ivarnhm)%p%field(:,:,k,it)
		     end do 
		  end do
	       else
		  do it = 1, x%nt
		     self(ie)%control(iweight)%p%field(:,:,1:nlev,1) = self(ie)%control(iweight)%p%field(:,:,1:nlev,1) + &
								     & x%control(ivarnhm)%p%field(:,:,1:nlev,it)*xpert(ie)%control(ivarnhm)%p%field(:,:,1:nlev,it) 
		  end do
	       end if
            end do
            self(ie)%control(iweight)%p%field(:,:,1:nlev,1) = weight_ens*self(ie)%control(iweight)%p%field(:,:,1:nlev,1)
         end do
      end do
      !
      return
   end subroutine apply_USITens
   !
   !
   !
   subroutine apply_USITclim(self, sigma, x)
   ! no temporal localization
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeControl), intent(in) :: sigma, x
      integer :: nlev, ivar, iweight, ivarnhm, it
      !
      do iweight = 1, self%ncontrol
         call get_nlev(self, iweight, nlev)
         self%control(iweight)%p%field(:,:,1:nlev,1) = 0.d0
         do ivar = 1, self%nvarnhm(iweight)
            ivarnhm = self%ivarnhm(iweight,ivar)
            do it = 1, x%nt
               self%control(iweight)%p%field(:,:,1:nlev,1) = self%control(iweight)%p%field(:,:,1:nlev,1) + &
                                                           & x%control(ivarnhm)%p%field(:,:,1:nlev,it)*sigma%control(iweight)%p%field(:,:,1:nlev,1)
            end do
         end do
         self%control(iweight)%p%field(:,:,1:nlev,1) = weight_clim*self%control(iweight)%p%field(:,:,1:nlev,1)
      end do
      !
      return
   end subroutine apply_USITclim
   !
   !
   !
   subroutine apply_U1000(self, ne, w, xpert)
   ! similar to apply_USIens but for sensitivity test
      implicit none
      type(NodeControl), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeControl), dimension(1:ne), intent(in) :: w
      type(NodeControl), dimension(ne), intent(in) :: xpert
      integer :: nlev, ivar, k, it, ie
      !
      do ivar = 1, self%ncontrol
	 call get_nlev(self, ivar, nlev)
	 do it = 1, self%nt
	    do k = 1, nlev
	       self%control(ivar)%p%field(:,:,k,it) = 0.d0
	       do ie = 1, ne
		  self%control(ivar)%p%field(:,:,k,it) = self%control(ivar)%p%field(:,:,k,it) + &
						       & w(ie)%control(1)%p%field(:,:,1,1)*xpert(ie)%control(ivar)%p%field(:,:,k,it)
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine apply_U1000
   !
   !
   !
   subroutine apply_UT1000(self, ne, x, xpert)
   ! similar to apply_USITens but for sensitivity test
      implicit none
      type(NodeControl), dimension(ne), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeControl), intent(in) :: x
      type(NodeControl), dimension(ne), intent(in) :: xpert
      integer :: nlev, ivar, k, it, ie
      !
      do ie = 1, ne
         self(ie)%control(1)%p%field(:,:,1,1) = 0.d0
         do ivar = 1, x%ncontrol
	    call get_nlev(x, ivar, nlev)
	    do it = 1, x%nt
	       do k = 1, nlev
	          self(ie)%control(1)%p%field(:,:,1,1) = self(ie)%control(1)%p%field(:,:,1,1) + &
						       & x%control(ivar)%p%field(:,:,k,it)*xpert(ie)%control(ivar)%p%field(:,:,k,it)
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine apply_UT1000
   !
   !
   !
   subroutine update_NodeControl(self, info)
      implicit none
      type(NodeControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer :: nx, ny, nlev, nt, ivar, it
      !
      nx = self%nx
      ny = self%ny
      do ivar = 1, self%ncontrol
         call get_nlev(self, ivar, nlev)
	 call get_nt(self, ivar, nt)
         do it = 1, nt
            call update_halo(info, self%control(ivar)%p%field(:,:,:,it), nx, ny, nlev)
         end do
      end do
      !
      return
   end subroutine update_NodeControl
   !
   !
   !
   subroutine apply_ensemble_isometry(self, ne, e, v)
      use enkflib, only : apply_isometry
      implicit none
      type(NodeControl), dimension(ne), intent(inout) :: self
      integer, intent(in) :: ne
      real(r_size), dimension(ne), intent(in) :: e
      real(r_size), dimension(ne,ne), intent(in) :: v
      integer :: nlev, ivar, it, k, j, i, ie
      real(r_size), dimension(ne) :: x
      !
      do ivar = 1, self(1)%ncontrol
         call get_nlev(self(1), ivar, nlev)
         if (nlev > 1) nlev = nlev - 1
         do it = 1, self(1)%nt
            do k = 1, nlev
               do j = 1, self(1)%ny
                  do i = 1, self(1)%nx
                     do ie = 1, ne
                        x(ie) = self(ie)%control(ivar)%p%field(i,j,k,it)
                     end do
                     call apply_isometry(ne, e, v, x)
                     do ie = 1, ne
                        self(ie)%control(ivar)%p%field(i,j,k,it) = x(ie)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !
      return
   end subroutine apply_ensemble_isometry
   !
   !
   !
end module NodeControl_class
   
