module enkflib
! Author: Le Duc
! Created date: 10 Apr 2016
   use variable, only : r_size, r_dble
contains
   !
   !
   !
   subroutine GaspariCohn(z, c2, weight)
   ! c2 is 2c parameter
      implicit none
      real(r_size), intent(in) :: z, c2
      real(r_size), intent(out) :: weight
      real(r_size) :: x
      !
      x = 2.*abs(z)/c2
      if (x <= 1.) then
	 weight = -0.25*x**5  + 0.5*x**4 + 0.625*x**3 - 5./3.*x**2 + 1.
      else if (x <= 2.) then
	 weight = 1./12.*x**5 - 0.5*x**4 + 0.625*x**3 + 5./3.*x**2 - 5.*x + 4. - 2./3.*x**(-1)
      else
	 weight = 0.
      end if
      !
      return
   end subroutine GaspariCohn
   !
   !
   !
   !=======================================================================
   !  Eigenvalue decomposition using subroutine rs
   !    INPUT
   !      INTEGER :: imode           : mode switch (0: only eiven values)
   !      INTEGER :: n               : dimension of matrix
   !      REAL(r_size) :: a(n,n)     : input matrix
   !    OUTPUT
   !      REAL(r_size) :: eival(n)   : eiven values in decending order
   !                                   i.e. eival(1) is the largest
   !      REAL(r_size) :: eivec(n,n) : eiven vectors
   !      INTEGER :: nrank_eff       : number of positive eivenvalues
   !=======================================================================
   SUBROUTINE mtx_eigen(imode,n,a,eival,eivec,nrank_eff)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: imode ! 0: calculate only eigen values
      INTEGER,INTENT(IN) :: n
      REAL(r_size),INTENT(IN) :: a(1:n,1:n)
      REAL(r_size),INTENT(OUT) :: eival(1:n)
      REAL(r_size),INTENT(OUT) :: eivec(1:n,1:n)
      INTEGER,INTENT(OUT) :: nrank_eff
      !
      REAL(r_dble) :: a8(n,n)
      REAL(r_dble) :: eival8(n)
      REAL(r_dble) :: eivec8(n,n)
      REAL(r_dble) :: wrk1(n)
      REAL(r_dble) :: wrk2(n)
      INTEGER :: ierr,i,j
      !
      a8 = a
      eivec8 = 0.0d0
      CALL rs(n,n,a8,eival8,imode,eivec8,wrk1,wrk2,ierr)
      IF (ierr /= 0) THEN
         PRINT *,'!!! ERROR (mtx_eigen): rs error code is ',ierr
         STOP 2
      END IF
      !
      nrank_eff = n
      IF( eival8(n) > 0 ) THEN
         DO i=1,n-1
	    IF( eival8(i) < ABS(eival8(n))*SQRT(EPSILON(eival8(n))) ) THEN
	       nrank_eff = nrank_eff - 1
	       eival8(i) = 0.0d0
	    END IF
         END DO
      ELSE
         PRINT *,'!!! ERROR (mtx_eigen): All Eigenvalues are below 0'
         STOP 2
      END IF
      !
      DO i=1,n
         eival(i) = eival8(n+1-i)
         eivec(:,i) = eivec8(:,n+1-i)
      END DO
      !
      RETURN
   END SUBROUTINE mtx_eigen
   !
   !
   !
   subroutine random_ensemble_isometry(n, e, v)
      use variable, only : r_size, pi
      implicit none
      integer, intent(in) :: n
      real(kind=r_size), dimension(n), intent(out) :: e
      real(kind=r_size), dimension(n,n), intent(out) :: v
      integer :: i, j
      real(kind=r_size) :: norm, product, phi
      !
      ! First vector: unit vector
      call random_number(phi)
      if (phi < 0.5) then
	 e(1) = -1.
      else
	 e(1) = 1.
      end if
      v(:,1) = 1.d0
      !
      do j = 2, n
         norm = 0.d0
	 do while (norm < 1.e-12)
	    call random_number(v(:,j))
	    norm = sqrt(sum(v(:,j)**2))
         end do
      end do
      !
      ! Modified Gram-Schmidt
      do j = 1, n
         norm = sqrt(sum(v(:,j)**2))
         do while (norm < 1.e-12)
	    call random_number(v(:,j))
	    do i = 1, j-1
	       product = sum(v(:,j)*v(:,i))
	       v(:,j) = v(:,j) - product*v(:,i)
	    end do
	    norm = sqrt(sum(v(:,j)**2))
	 end do
         v(:,j) = v(:,j)/norm
         do i = j+1, n
	    product = sum(v(:,j)*v(:,i))
	    v(:,i) = v(:,i) - product*v(:,j)
	 end do
      end do
      !
      ! Random rotation: e contains cos and sin at e(j) and e(j+1)
      do j = 2, n-1, 2
         call random_number(phi)
	 phi = 2.d0*pi*phi - pi
	 e(j) = cos(phi)
	 e(j+1) = sin(phi)
      end do
      ! Reflexion for the last item
      if (mod(n,2) == 0) e(n) = -1
      !
      return
   end subroutine random_ensemble_isometry
   !
   !
   !
   subroutine apply_isometry(n, e, v, x)
      use variable, only : r_size
      implicit none
      integer, intent(in) :: n
      real(kind=r_size), dimension(n), intent(in) :: e
      real(kind=r_size), dimension(n,n), intent(in) :: v
      real(kind=r_size), dimension(n), intent(inout) :: x
      integer :: i
      real(kind=r_size), dimension(n) :: y
      !
      ! Coordinates in V
      do i = 1, n
	 y(i) = sum(x(:)*v(:,i))
      end do
      x(:) = y(:)
      !
      ! Modify coordinates by E
      y(1) = x(1)*e(1)
      do i = 2, n-1, 2
         y(i) = e(i)*x(i) - e(i+1)*x(i+1)
	 y(i+1) = e(i+1)*x(i) + e(i)*x(i+1)
      end do
      ! Reflexion for the last item
      if (mod(n,2) == 0) y(n) = x(n)*e(n)
      !
      ! Add V
      x(:) = 0.d0
      do i = 1, n
         x(:) = x(:) + y(i)*v(:,i)
      end do
      !
      return
   end subroutine apply_isometry
   !
   !
   !
   subroutine test_isometry(n, e, v)
      use variable, only : r_size
      implicit none
      integer, intent(in) :: n
      real(kind=r_size), dimension(n), intent(in) :: e
      real(kind=r_size), dimension(n,n), intent(in) :: v
      integer :: i, j
      real(kind=r_size) :: product
      real(kind=r_size), dimension(n,n) :: x
      !
      x(:,:) = 0.d0
      do i = 1, n
	 x(i,i) = 1.d0
      end do
      do i = 1, n
         call apply_isometry(n, e, v, x(:,i))
      end do
      !
      !do i = 1, n
         !do j = i, n
            !product = sum(v(:,j)*v(:,i))
            !print*, i, j, product
         !end do
      !end do
      do i = 1, n
         do j = i, n
            product = sum(x(:,j)*x(:,i))
            print*, i, j, product
         end do
      end do
      !
      return
   end subroutine test_isometry
   !
   !
   !
end module enkflib
