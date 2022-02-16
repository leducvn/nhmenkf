module lbfgsb_driver
   use variable, only : r_size, r_sngl, r_dble
   implicit none
   !
   integer(4), parameter :: m = 10, iprint = -1
   real(kind=r_dble), parameter :: factr = 1.0d+7, pgtol = 1.0d-5, eps = 1.0d-10
   character(len=60), private :: task, csave
   logical, private :: lsave(4)
   integer, private :: isave(44)
   real(kind=r_dble) :: dsave(29)
   integer(4), private, allocatable :: nbd(:)
   integer(4), private, allocatable :: iworka(:) ! iwa in lbfgsb
   real(kind=r_dble), private, allocatable :: lbd(:), ubd(:) ! l, u in lbfgsb
   real(kind=r_dble), private, allocatable :: worka(:) ! wa in lbfgsb
   real(kind=r_dble), private, allocatable :: fwork, xwork(:), gwork(:)
   !
   !
   !
contains
   !
   !
   !
   subroutine Initialize_Minimizer(n)
      implicit none
      integer(4), intent(in) :: n
      !
      allocate(nbd(n), iworka(3*n))
      allocate(lbd(n), ubd(n))
      allocate(worka(2*m*n+5*n+11*m*m+8*m))
      allocate(xwork(n), gwork(n))
      nbd(:) = 0
      lbd(:) = 0.
      ubd(:) = 0.
      xwork(:) = 0.
      task = 'START'
      call setulb(n, m, xwork, lbd, ubd, nbd, fwork, gwork, factr, pgtol, &
                & worka, iworka, task, iprint, csave, lsave, isave, dsave)
      !
      return
   end subroutine Initialize_Minimizer
   !
   !
   !
   subroutine Minimize(n, x, f, g, iflag)
   !      iflag = 1      : continue to minimize
   !      iflag =   0    : succeed to find the minimum point
   !            <   0    : abnormal termination by various reason.
      implicit none
      integer(4), intent(in) :: n
      integer(4), intent(out) :: iflag
      real(kind=r_size), intent(inout) :: x(n)
      real(kind=r_size), intent(inout) :: f
      real(kind=r_size), intent(inout) :: g(n)
      integer :: i
      !
      if(r_sngl == r_size) then
         fwork = real(f, r_dble)
         do i = 1, n
            gwork(i) = real(g(i), r_dble)
         end do
         call setulb(n, m, xwork, lbd, ubd, nbd, fwork, gwork, factr, pgtol, &
                   & worka, iworka, task, iprint, csave, lsave, isave, dsave)
	 do while (task(1:5) == 'NEW_X')
	    if (dsave(13) <= eps*(1.0d0+abs(f))) then
	       task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
	       exit
	    end if
	    call setulb(n, m, xwork, lbd, ubd, nbd, fwork, gwork, factr, pgtol, &
                      & worka, iworka, task, iprint, csave, lsave, isave, dsave)
	 end do
         do i = 1, n
            x(i) = real(xwork(i), r_sngl)
         end do
      else
         call setulb(n, m, xwork, lbd, ubd, nbd, f, g, factr, pgtol, &
                   & worka, iworka, task, iprint, csave, lsave, isave, dsave)
         do while (task(1:5) == 'NEW_X')
	    if (dsave(13) <= eps*(1.0d0+abs(f))) then
	       task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
	       exit
	    end if
	    call setulb(n, m, xwork, lbd, ubd, nbd, f, g, factr, pgtol, &
                      & worka, iworka, task, iprint, csave, lsave, isave, dsave)
	 end do
	 x(:) = xwork(:)
      end if
      !
      if (task(1:2) == 'FG') then
	 iflag = 1
      else if (task(1:4) == 'STOP') then
	 iflag = 0
      else if (task(1:4) == 'CONV') then
	 iflag = 0
      else if (task(1:4) == 'ABNO') then
	 iflag = -1
      else if (task(1:5) == 'ERROR') then
	 iflag = -2
      else
	 iflag = -99
      end if
      !
      return
   end subroutine Minimize
   !
   !
   !
   subroutine Terminate_Minimizer
      implicit none
      !
      deallocate(nbd, iworka)
      deallocate(lbd, ubd)
      deallocate(worka)
      deallocate(xwork, gwork)
      !
      return
   end subroutine Terminate_Minimizer
   !
   !
   !
end module lbfgsb_driver
