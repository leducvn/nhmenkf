! $Id: matrix_svd.f90 3551 2013-02-25 09:51:28Z idculv $

!****f* Matrix/matrix_svd *
!
! NAME
!    matrix_svd - Compute singular value decomposition (SVD) of a matrix.
!
! SYNOPSIS
!    call matrix_svd(A, W, U, V)
! 
! DESCRIPTION
!    This subroutine calculates the SVD of a matrix.
!               A = U diag(W) V^T
!    where U and V are matrices containing eigenvectors of A and
!    the diagonal elements of W give the singular values.
!
! INPUTS
!    A   Matrix to be used for calculating the SVD
!
! OUTPUT
!    W   Diagonal elements of singular values of A
!    U   Left-hand singular vector of A
!    V   Right-hand singular vector of A. (Note routine returns V and NOT its 
!        transpose V^T).
!
! NOTES
!   The routine is based on the Numerical Recipes (Press et al. 1992) algorithm
!   for computing the SVD. Only full matrices are supported. Matrices in packed
!   form should be converted to full type before calling matrix_svd.
!
! EXAMPLE
!
! SEE ALSO
!      matrix_sqrt
!
! REFERENCES
!   W.H. Press, S.A. Teukolsjy, W.T. Vetterling and B.P. Flannery,
!   Numerical Recipes in C - The Art of Scientific Computing.
!   2nd Ed., Cambridge University Press, 1992.
!
!    G.H. Golub, C.H. van Loan
!    Matrix computations. 
!    3rd Ed., The Johns Hopkins University Press, 1996.
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

subroutine matrix_svd(A, W, U, V)

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  use typesizes, only: wp => EightByteReal

  implicit none
  
  real(wp), dimension(:,:)  ,               intent(in)  ::  A
  real(wp), dimension(size(A,1)),           intent(out) ::  W
  real(wp), dimension(size(A,1),size(A,2)), intent(out) ::  U
  real(wp), dimension(size(A,2),size(A,2)), intent(out) ::  V

  integer :: n, m
  
  integer  :: flag, i, its, j, jj, k, l, nm
  real(wp) :: anorm, c, f, g, h, s, scale, x, y, z
  real(wp), dimension(size(A,2)) :: rv1
  
  m = size(A,1)
  n = size(A,2)
  
!-------------------------------------------------------------------------------
! 2. Householder reduction to bidiagonal form
!-------------------------------------------------------------------------------

  rv1(:) = 0.0_wp
  g = 0.0_wp
  scale = 0.0_wp
  anorm = 0.0_wp
  
  U(:,:) = A(:,:)

  do i=1,n
     l = i+1
     rv1(i) = scale*g
     g = 0.0_wp
     s = 0.0_wp
     scale = 0.0_wp

     if(i <= m)then
        
        do k=i,m
           scale = scale + abs(U(k,i))
        enddo

        if(scale /= 0.0_wp)then
           do k=i,m
              U(k,i) = U(k,i)/scale
              s = s + U(k,i)*U(k,i)
           enddo

           f = U(i,i)
           g = -SIGN(sqrt(s),f)           
           h = f*g - s
           U(i,i) = f - g

           do j=l,n
              s = 0.0_wp
              do k=i,m
                 s = s + U(k,i)*U(k,j)
              enddo
              f = s/h
              do k=i,m
                 U(k,j) = U(k,j) + f*U(k,i)
              enddo
           enddo

           do k=i,m
              U(k,i) = scale*U(k,i)
           enddo

        endif
     endif

     w(i) = scale*g
     
     g = 0.0_wp
     s = 0.0_wp
     scale = 0.0_wp

     if(i <= m .and. i /= n)then
        
        do k=l,n
           scale = scale + abs(U(i,k))
        enddo

        if(scale /= 0.0_wp)then
           
           do k=l,n
              U(i,k) = U(i,k)/scale
              s = s + U(i,k)*U(i,k)
           enddo
           
           f = U(i,l)           
           g = -SIGN(sqrt(s),f)          
           h = f*g - s          
           U(i,l) = f - g
           
           do k=l,n
              rv1(k) = U(i,k)/h
           enddo
           
           do j=l,m
              s = 0.0_wp
              do k=l,n
                 s = s + U(j,k)*U(i,k)
              enddo
              
              do k=l,n
                 U(j,k) = U(j,k) + s*rv1(k)
              enddo
           enddo
           
           do k=l,n
              U(i,k) = scale*U(i,k)
           enddo
        endif
     endif

     anorm = MAX(anorm, (abs(w(i))+abs(rv1(i))))

  enddo

!-------------------------------------------------------------------------------
! 3. Accumulation of right-hand transformations
!-------------------------------------------------------------------------------

  do i=n,1,-1
     
     if(i < n)then

        if(g /= 0.0_wp)then
           
           ! Double division to avoid possible underflow
           do j=l,n
              V(j,i) = (U(i,j)/U(i,l))/g
           enddo

           do j=l,n
              s = 0.0_wp
              do k=l,n
                 s = s + U(i,k)*V(k,j)
              enddo
              do k=l,n
                 V(k,j) = V(k,j) + s*V(k,i)
              enddo
           enddo
        endif
        
        do j=l,n
           V(i,j) = 0.0_wp
           V(j,i) = 0.0_wp
        enddo

     endif

     V(i,i) = 1.0_wp
     g = rv1(i)
     l = i
  enddo

!-------------------------------------------------------------------------------
! 4. Accumulation of left-hand transformations
!-------------------------------------------------------------------------------

  do i=MIN(m,n),1,-1
     l = i+1
     g = W(i)

     do j=l,n
        U(i,j) = 0.0_wp
     enddo
     
     if(g /= 0.0_wp)then
        g = 1.0_wp/g
        
        do j=l,n
           s = 0.0_wp
           do k=l,m
              s = s + U(k,i)*U(k,j)
           enddo
           f = (s/U(i,i))*g
           do k=i,m
              U(k,j) = U(k,j) + f*U(k,i)
           enddo
        enddo
        
        do j=i,m
           U(j,i) = U(j,i)*g
        enddo
        
     else
        do j=i,m
           U(j,i) = 0.0_wp
        enddo
     endif
     
     U(i,i) = U(i,i) + 1.0_wp
     
  enddo

!-------------------------------------------------------------------------------
! 5. Diagonalisation of the bidiagonal form
!-------------------------------------------------------------------------------

! 5.1 Loop over singular values, and over allowed iterations.
  do k=n,1,-1

     do its=1,30

        flag = 1

        ! Test for splitting (note rv1(1) always zero)
        do l=k,1,-1
           nm = l-1
           if(abs(rv1(l)) + anorm == anorm)then
              flag = 0
              exit
           endif
           
           if(abs(W(nm)) + anorm == anorm)then
              exit
           endif
        enddo

! 5.2 Cancellation of rv1(l) if l>1
       if(flag > 0)then
          
           c = 0.0_wp
           s = 1.0_wp

           do i=l,k
              f = s*rv1(i)
              rv1(i) = c*rv1(i)

              if(abs(f)+anorm == anorm)then
                 exit
              endif

              g = W(i) 
              h = PYTHAG(f, g)
              W(i) = h
              h = 1.0_wp/h
              c = g*h
              s = -f*h

              do j=1,m
                 y = U(j,nm)
                 z = U(j,i)
                 U(j,nm) = y*c + z*s
                 U(j,i) = z*c - y*s
              enddo

           enddo
        endif

! 5.3 Check for convergence

        z = W(k)
        if(l == k)then
           
           if(z < 0.0_wp)then  !! Ensure positive singular values
              W(k) = -z
              do j=1,n
                 V(j,k) = -V(j,k)
              enddo
           endif

           exit
        endif

        if(its == 30)then
           print*, 'No convergence in 30 iterations'
        endif

! 5.4 Shift from bottom 2x2 minor
        x = W(l)
        nm = k-1

        y = W(nm)
        g = rv1(nm)
        h = rv1(k)

        f = ((y-z)*(x+z)+(g-h)*(g+h))/(2.0_wp*h*y)
        g = PYTHAG(f, 1.0_wp)
        f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x

! 5.5 Next QR transformation
        c = 1.0_wp
        s = 1.0_wp

        do j=l,nm
           i = j+1

           g = rv1(i)
           y = W(i)
           h = s*g
           g = c*g

           z = PYTHAG(f,h)
           rv1(j) = z
         
           c = f/z
           s = h/z
           f = x*c + g*s
           g = g*c - x*s
           h = y*s
           y = y*c

           do jj=1,n
              x = V(jj,j)
              z = V(jj,i)
              V(jj,j) = x*c + z*s
              V(jj,i) = z*c - x*s
           enddo

           z = PYTHAG(f,h)
           
           W(j) = z     !! Rotation can be arbitary if z=0
           if(z > 0)then
              z = 1.0_wp/z
              c = f*z
              s = h*z
           endif

           f = c*g + s*y
           x = c*y - s*g

           do jj=1,m
              y = U(jj,j)
              z = U(jj,i)
              U(jj,j) = y*c + z*s
              U(jj,i) = z*c - y*s
           enddo

        enddo

        rv1(l) = 0.0_wp
        rv1(k) = f
        W(k) = x

     enddo

  enddo

contains

!-------------------------------------------------------------------------------
! 6. Compute (a2+b2)^0.5 (without underflow or overflow)
!-------------------------------------------------------------------------------
  function pythag(a_pyth, b_pyth) result (out_pyth)

    implicit none

    real(wp) :: a_pyth
    real(wp) :: b_pyth
    real(wp) :: out_pyth
    real(wp) :: absa, absb

    absa = abs(a_pyth)
    absb = abs(b_pyth)

    if(absa > absb)then
       out_pyth = absa*sqrt(1.0_wp + (absb/absa)**2)
    else if(absb > 0.0_wp)then
       out_pyth = absb*sqrt(1.0_wp + (absa/absb)**2)
    else
       out_pyth = 0.0_wp
    endif
    
  end function pythag

      
end subroutine matrix_svd
