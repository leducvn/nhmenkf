module Localization_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsControl_class
   use NodeControl_class
   implicit none
   !
   type Localization
   ! hscale is the parameter 2c in Gaspari and Cohn function: distance where influence is 0
      integer :: nvar, nobsvar
      character(len=10), dimension(:), allocatable :: varname, obsname
      logical, dimension(:,:), allocatable :: correlated
      real(r_size), dimension(:,:), allocatable :: hscale
   end type Localization
   !
   interface new
      module procedure new_Localization
   end interface
   interface destroy
      module procedure destroy_Localization
   end interface
   interface display
      module procedure display_Localization
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_Localization(self, hscale, hscaleqv, x, y)
      implicit none
      type(Localization), intent(inout) :: self
      type(NodeControl), intent(in) :: x
      type(NodeObsControl), intent(in) :: y
      character(len=10) :: varname
      integer :: ivar, nvar, nobsvar
      real(r_size) :: hscale, hscaleqv
      !
      call get_ncontrol(x, nvar)
      call get_ncontrol(y, nobsvar)
      self%nvar = nvar
      self%nobsvar = nobsvar
      allocate(self%varname(nvar), self%obsname(nobsvar))
      do ivar = 1, nvar
         call get_name(x, ivar, varname)
         self%varname(ivar) = varname
      end do
      do ivar = 1, nobsvar
         call get_name(y, ivar, varname)
         self%obsname(ivar) = varname
      end do
      !
      allocate(self%correlated(nvar,nobsvar))
      allocate(self%hscale(nvar,nobsvar))
      self%correlated(:,:) = .True.
      self%hscale(:,:) = hscale
      !
      call set_correlated1(self, 1, 'tgrd', .False.)
      call set_correlated2(self, 'tgrd', 't', .True.)
      call set_correlated2(self, 'tgrd', 'rh', .True.)
      !
      call set_correlated1(self, 1, 'ps', .False.)
      call set_correlated2(self, 'ps', 'u', .True.)
      call set_correlated2(self, 'ps', 'v', .True.)
      call set_correlated2(self, 'ps', 't', .True.)
      call set_correlated2(self, 'ps', 'p', .True.)
      call set_correlated2(self, 'ps', 'tc', .True.)
      call set_correlated2(self, 'ps', 'rvl', .True.)
      !
      call set_hscale1(self, 1, 'qv', hscaleqv)
      call set_hscale1(self, 2, 'pwv', hscaleqv)
      call set_hscale1(self, 2, 'rh', hscaleqv)
      !
      call set_correlated2(self, 'qv', 'p', .False.)
      !call set_hscale1(self, 2, 'p', 0.7*hscale)
      !call set_hscale2(self, 'ps', 'p', 1.3*hscale)
      !
      return
   end subroutine new_Localization
   !
   !
   !
   subroutine destroy_Localization(self)
      implicit none
      type(Localization), intent(inout) :: self
      !
      if (allocated(self%varname)) deallocate(self%varname)
      if (allocated(self%obsname)) deallocate(self%obsname)
      if (allocated(self%correlated)) deallocate(self%correlated)
      if (allocated(self%hscale)) deallocate(self%hscale)
      !
      return
   end subroutine destroy_Localization
   !
   !
   !
   subroutine display_Localization(self)
      implicit none
      type(Localization), intent(in) :: self
      !
      print*, 'Localization: ', self%nvar, self%nobsvar
      !
      return
   end subroutine display_Localization
   !
   !
   !
   subroutine get_correlated(self, ivar, iobsvar, correlated)
      implicit none
      type(Localization), intent(in) :: self
      integer, intent(in) :: ivar, iobsvar
      logical, intent(out) :: correlated
      !
      correlated = self%correlated(ivar,iobsvar)
      !
      return
   end subroutine get_correlated
   !
   !
   !
   subroutine set_correlated1(self, varmode, varname, correlated)
      implicit none
      type(Localization), intent(inout) :: self
      integer, intent(in) :: varmode
      character(len=*), intent(in) :: varname
      logical, intent(in) :: correlated
      logical :: found
      integer :: ivar
      !
      found = .False.
      if (varmode == 1) then
         do ivar = 1, self%nvar
	    if (trim(varname) == trim(self%varname(ivar))) then
	       found = .True.
	       exit
	    end if
         end do
      else if (varmode == 2) then
         do ivar = 1, self%nobsvar
	    if (trim(varname) == trim(self%obsname(ivar))) then
	       found = .True.
	       exit
	    end if
         end do
      else
         return
      end if
      if (.not. found) return
      if (varmode == 1) then
         self%correlated(ivar,:) = correlated
      else
         self%correlated(:,ivar) = correlated
      end if
      !
      return
   end subroutine set_correlated1
   !
   !
   !
   subroutine set_correlated2(self, varname, obsname, correlated)
      implicit none
      type(Localization), intent(inout) :: self
      character(len=*), intent(in) :: varname, obsname
      logical, intent(in) :: correlated
      logical :: found
      integer :: ivar, iobs
      !
      found = .False.
      do ivar = 1, self%nvar
	 if (trim(varname) == trim(self%varname(ivar))) then
	    found = .True.
	    exit
	 end if
      end do
      if (.not. found) return
      found = .False.
      do iobs = 1, self%nobsvar
	 if (trim(obsname) == trim(self%obsname(iobs))) then
	    found = .True.
	    exit
	 end if
      end do
      if (.not. found) return
      self%correlated(ivar,iobs) = correlated
      !
      return
   end subroutine set_correlated2
   !
   !
   !
   subroutine get_hscale(self, ivar, iobsvar, hscale)
      implicit none
      type(Localization), intent(in) :: self
      integer, intent(in) :: ivar, iobsvar
      real(r_size), intent(out) :: hscale
      !
      hscale = self%hscale(ivar,iobsvar)
      !
      return
   end subroutine get_hscale
   !
   !
   !
   subroutine set_hscale1(self, varmode, varname, hscale)
      implicit none
      type(Localization), intent(inout) :: self
      integer, intent(in) :: varmode
      character(len=*), intent(in) :: varname
      real(r_size), intent(in) :: hscale
      logical :: found
      integer :: ivar
      !
      found = .False.
      if (varmode == 1) then
         do ivar = 1, self%nvar
	    if (trim(varname) == trim(self%varname(ivar))) then
	       found = .True.
	       exit
	    end if
         end do
      else if (varmode == 2) then
         do ivar = 1, self%nobsvar
	    if (trim(varname) == trim(self%obsname(ivar))) then
	       found = .True.
	       exit
	    end if
         end do
      else
         return
      end if
      if (.not. found) return
      if (varmode == 1) then
         self%hscale(ivar,:) = hscale
      else
         self%hscale(:,ivar) = hscale
      end if
      !
      return
   end subroutine set_hscale1
   !
   !
   !
   subroutine set_hscale2(self, varname, obsname, hscale)
      implicit none
      type(Localization), intent(inout) :: self
      character(len=*), intent(in) :: varname, obsname
      real(r_size), intent(in) :: hscale
      logical :: found
      integer :: ivar, iobs
      !
      found = .False.
      do ivar = 1, self%nvar
	 if (trim(varname) == trim(self%varname(ivar))) then
	    found = .True.
	    exit
	 end if
      end do
      if (.not. found) return
      found = .False.
      do iobs = 1, self%nobsvar
	 if (trim(obsname) == trim(self%obsname(iobs))) then
	    found = .True.
	    exit
	 end if
      end do
      if (.not. found) return
      self%hscale(ivar,iobs) = hscale
      !
      return
   end subroutine set_hscale2
   !
   !
   !
end module Localization_class