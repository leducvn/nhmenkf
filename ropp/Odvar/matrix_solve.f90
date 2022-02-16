! $Id: matrix_solve.f90 4010 2014-01-10 11:07:40Z idculv $

!****f* Matrix/matrix_solve *
!
! NAME
!    matrix_solve - Solve a linear matrix equation
!
! SYNOPSIS
!    use matrix
!      ...
!    x = matrix_solve(A, b)
! 
! DESCRIPTION
!    This function calls matrix_solve to solve a linear equation of the form
!
!       A x = b
!
!    for arbitrary matrices A and arbitrary right hand side vectors b. 
!    It is assumed that the matrix A is positive 
!    definite, and the solution is obtained using a Cholesky decomposition. 
!
! INPUTS
!    A     Matrix defining the linear equation system.
!    b     Right hand side of the linear equation system.
!
! OUTPUT
!    x     Solution of the linear equation system.
!
! NOTES
!    The matrix A may be in either full or packed form. If the matrix is not 
!    positive definite, the attempt to solve will fail with an error message.
!    These routines are currently available in double precision only.
!
! EXAMPLE
!    The solution of a linear equation with right hand side vector b is
!
!       x = solve(A, b)
!
! SEE ALSO
!    matrix_types
!
! REFERENCES
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

!-------------------------------------------------------------------------------
! 1. Full matrices, double precision (1D RHS)
!-------------------------------------------------------------------------------
function matrix_solve_gen_1d(A, b) result(x)

! 1.1 Declarations

  use typesizes, only: wp => EightByteReal
  use messages
  use matrix_types
  use ropp_utils, only: ropp_MDFV

  implicit none

  real(wp), dimension(:,:), intent(inout)  :: A
  real(wp), dimension(:),   intent(in)  :: b
  real(wp), dimension(size(b))          :: x

  real(wp), dimension(size(b), size(b)) :: G
  real(wp), dimension(size(b))          :: temp
  real(wp), dimension(size(b))          :: y
  integer                               :: i, j, n
  character(len = 256)                  :: routine

! 1.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_solve')
  
! 1.3 Check that matrix is a square matrix

  n = size(A,1)

  if(n /= size(A,2))then
     call message(msg_error, 'LHS matrix to be solved for ' // trim(routine) // ' is not square - aboriting calculation')
     A(:,:) = ropp_MDFV
     x(:) = ropp_MDFV
     RETURN
  endif

! 1.4 Check that matrix and RHS vector have consistent dimensions

  if(n /= size(b))then
     call message(msg_error, &
          'RHS matrix has different dimensions from LHS matrix to be solved for ' &
          // trim(routine) // ' - aboriting calculation')
     A(:,:) = ropp_MDFV
     x(:) = ropp_MDFV
     RETURN
  endif

! 1.5 Determine Cholesky matrix G

   G(:,:) = 0.0_wp 
  
   do i = 1, n
      
      temp(i:n) = A(i:n, i)      ! store lower triangle matrix
      
      if(i /= 1)then         
         do j = 1, i-1
            temp(i:n) = temp(i:n) - G(i,j)*G(i:n, j)
         enddo
      endif
      
!  1.6 Check that matrix is positive definite
      
      if(temp(i) <= TINY(0.0_wp)*100.0_wp)then
         call message(msg_error, 'Matrix to be inverted in '// trim(routine) // ' is not positive definite - aborting')
         A(:,:) = ropp_MDFV
         x(:) = ropp_MDFV
         return
      endif
      
!  1.7 Compute elements of Cholesky matrix G

      G(i:n, i) = temp(i:n) / sqrt(temp(i))
      
   enddo

!  1.8 Solve G.y = b for y by forward substitution
   
   y = b
   
   y(1) = y(1) / G(1,1)
   do j = 2, n 
      y(j) = (y(j) - dot_product( G(j,1:j-1),y(1:j-1))) / G(j,j)
   enddo
   
!  1.9 Solve G^T.x = y for x by backward substitution

   x = y
   x(n) = x(n) / G(n,n)
   do j = n-1, 1, -1
      x(j) = (x(j) - dot_product(G(j+1:n,j),x(j+1:n))) / G(j,j)
   enddo
   
! 1.10 Clean up
   
   call message_set_routine(routine)
   
end function matrix_solve_gen_1d
 

!-------------------------------------------------------------------------------
! 2. Full matrices, double precision (2D RHS)
!-------------------------------------------------------------------------------
function matrix_solve_gen_2d(A, b) result(x)

! 2.1 Declarations

  use typesizes, only: wp => EightByteReal
  use messages
  use matrix_types
  use ropp_utils, only: ropp_MDFV

  implicit none

  real(wp), dimension(:,:), intent(inout)  :: A
  real(wp), dimension(:,:),   intent(in)   :: b
  real(wp), dimension(size(b,1),size(b,2)) :: x

  real(wp), dimension(size(b,1),size(b,1)) :: G
  real(wp), dimension(size(b,1))           :: temp
  real(wp), dimension(size(b,1))           :: y
  integer                                  :: i, j, k, n
  character(len = 256)                     :: routine
  
! 2.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_solve')
  
! 2.3 Check that matrix is a square matrix

  n = size(A,1)

  if(n /= size(A,2))then
     call message(msg_error, 'LHS matrix to be solved for ' // trim(routine) // ' is not square - aboriting calculation')
     A(:,:) = ropp_MDFV
     x(:,:) = ropp_MDFV
     RETURN
  endif

! 2.4 Check that matrix and RHS vector have consistent dimensions

  if(n /= size(b,1))then
     call message(msg_error, &
          'RHS matrix has different dimensions from LHS matrix to be solved for ' &
          // trim(routine) // ' - aboriting calculation')
     A(:,:) = ropp_MDFV
     x(:,:) = ropp_MDFV
     RETURN
  endif

! 2.5 Determine Cholesky matrix G

   G(:,:) = 0.0_wp 
  
   do i = 1, n
      
      temp(i:n) = A(i:n, i)      ! store lower triangle matrix
      
      if(i /= 1)then         
         do j = 1, i-1
            temp(i:n) = temp(i:n) - G(i,j)*G(i:n, j)
         enddo
      endif
      
!  2.6 Check that matrix is positive definite
      
      if(temp(i) <= TINY(0.)*100.0_wp)then
         call message(msg_error, 'Matrix to be inverted in '// trim(routine) // ' is not positive definite - aborting')
         A(:,:) = ropp_MDFV
         x(:,:) = ropp_MDFV
         return
      endif
       
!  2.7 Compute elements of Cholesky matrix G

      G(i:n, i) = temp(i:n) / sqrt(temp(i))
      
   enddo
   
!  2.8 Solve G.y = b for y by forward substitution

   do k = 1, size(b,2)
      
      y = b(k,:)
      
      y(1) = y(1) / G(1,1)
      do j = 2, n 
         y(j) = (y(j) - dot_product( G(j,1:j-1),y(1:j-1))) / G(j,j)
      enddo 
      
!  2.9 Solve G^T.x = y for x by backward substitution

      x(:,k) = y
      x(n,k) = x(n,k) / G(n,n)
      do j = n-1, 1, -1
         x(j,k) = (x(j,k) - dot_product(G(j+1:n,j),x(j+1:n,k))) / G(j,j)
      enddo
      
   enddo

! 2.10 Clean up
   
   call message_set_routine(routine)
   
end function matrix_solve_gen_2d

!-------------------------------------------------------------------------------
! 3. Positive definite packed matrix, double precision (1D RHS)
!-------------------------------------------------------------------------------

function matrix_solve_packed_1d(A, b) result(x)

! 3.1 Declarations

  use typesizes, only: wp => EightByteReal
  use ropp_utils, only: ropp_MDFV
  use messages
  use matrix_types
! use matrix, not_this => matrix_solve_packed_1d 

  implicit none

  type(matrix_pp),        intent(inout) :: A
  real(wp), dimension(:), intent(in)    :: b
  real(wp), dimension(size(b))          :: x

  real(wp), dimension(:,:), pointer     :: A_full => null() 
  real(wp), dimension(size(b), size(b)) :: G
  real(wp), dimension(size(b))          :: temp
  real(wp), dimension(size(b))          :: y
  integer                               :: i, j, n
  character(len = 256)                  :: routine
  
! 3.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_solve')

! 3.3 Check for consistent matrix and RHS vector dimensions

  n = (int(sqrt(8.0_wp * size(A % d)) + 1.0_wp) - 1.0_wp) / 2

  if(n /= size(b))then
     call message(msg_error, &
          'RHS matrix has different dimensions from LHS matrix to be solved for ' &
          // trim(routine) // ' - aboriting calculation')
     A%d(:) = ropp_MDFV
     x(:) = ropp_MDFV
     RETURN
  endif

! 3.4 Check if initialisation is needed

  if (.not. associated(A % f)) then
     A % fact_chol = .false.
  endif

! 3.5 Initialisation - compute Cholesky matrix (if not available from previous call to Cholesky)

  if(.not. A%fact_chol) then
     
     if (associated(A%e)) deallocate(A%e)
     if (associated(A%f)) deallocate(A%f)
     if (associated(A%s)) deallocate(A%s)
     
     allocate(A%f(size(A%d)))
     A%f(:) = 0.0_wp

     ! 3.5.1 Compute elements of full matrix
     
     call matrix_pp2full_alloc(A, A_full)
     
     ! 3.5.2 Check that matrix is a square matrix
     
     if(n /= size(A_full,2))then
        call message(msg_error, 'LHS matrix to be solved for ' // trim(routine) // ' is not square - aboriting calculation')
        A%d(:) = ropp_MDFV
        x(:) = ropp_MDFV
        RETURN
     endif
               
     ! 3.5.3 Determine Cholesky matrix G
     
     G = 0.0_wp 
     
     do i = 1, n
        
        temp(i:n) = A_full(i:n, i)      ! store lower triangle matrix
        
        if(i /= 1)then         
           do j = 1, i-1
              temp(i:n) = temp(i:n) - G(i,j)*G(i:n, j)
           enddo
        endif
        
        !  3.5.4 Check that matrix is positive definite
        
        if(temp(i) <= TINY(0.)*100.0_wp)then
           call message(msg_error, 'Matrix to be inverted in '// trim(routine) // ' is not positive definite - aborting')
           A%d(:) = ropp_MDFV
           x(:) = ropp_MDFV
           RETURN
        endif
        
        !  3.5.5 Compute elements of Cholesky matrix G

        G(i:n, i) = temp(i:n) / sqrt(temp(i))
        
     enddo
     
     ! 3.5.6 Store Cholesky matrix as packed type 
     call matrix_full2pp(G, A%f, 'L')
     
     ! 3.5.7 Set initialised flag in current matrix structure
     
     A % fact_chol = .true.
     deallocate(A_full)
     
  endif

! 3.7 Solve G.y = b for y by forward substitution

  y = b
  y(1) = y(1) / A%f(1)
  do j = 2, n 

      do i=1,j-1
        temp(i) = A%f(j+(2*n-i)*(i-1)/2)
     enddo
     
     y(j) = (y(j) - dot_product( temp(1:j-1),y(1:j-1))) / A%f(j+(2*n-j)*(j-1)/2)
  enddo
  
! 3.8 Solve G^T.x = y for x by backward substitution

  x = y
  x(n) = x(n) / A%f(n+n*(n-1)/2)
  do j = n-1, 1, -1

     do i=j+1,n
        temp(i) = A%f(i+(2*n-j)*(j-1)/2)
     enddo

     x(j) = (x(j) - dot_product(temp(j+1:n),x(j+1:n))) / A%f(j+(2*n-j)*(j-1)/2)
  enddo
  
! 3.9 Clean up

  call message_set_routine(routine)

end function matrix_solve_packed_1d

!-------------------------------------------------------------------------------
! 4. Positive definite packed matrix, double precision (2D RHS)
!-------------------------------------------------------------------------------

function matrix_solve_packed_2d(A, b) result(x)

! 4.1 Declarations

  use typesizes, only: wp => EightByteReal
  use ropp_utils, only: ropp_MDFV
  use messages
  use matrix_types
! use matrix, not_this => matrix_solve_packed_2d 

  implicit none

  type(matrix_pp),        intent(inout)     :: A
  real(wp), dimension(:,:), intent(in)      :: b
  real(wp), dimension(size(b,1),size(b,2))  :: x

  real(wp), dimension(:,:), pointer         :: A_full => null() 
  real(wp), dimension(size(b,1), size(b,1)) :: G
  real(wp), dimension(size(b,1))            :: temp
  real(wp), dimension(size(b,1))            :: y
  integer                                   :: i, j, k, n
  character(len = 256)                      :: routine

! 4.2 Error messages

  call message_get_routine(routine)
  call message_set_routine('matrix_solve')

! 4.3 Check for consistent matrix and RHS vector dimensions

  n = (int(sqrt(8.0_wp * size(A % d)) + 1.0_wp) - 1.0_wp) / 2

  if(n /= size(b,1))then
     call message(msg_error, &
          'RHS matrix has different dimensions from LHS matrix to be solved for ' &
          // trim(routine) // ' - aborting calculation')
     A%d(:) = ropp_MDFV
     x(:,:) = ropp_MDFV
     RETURN
  endif

! 4.4 Check if initialisation is needed

  if (.not. associated(A % f)) then
     A % fact_chol = .false.
  endif

! 4.5 Initialisation - compute Cholesky matrix (if not available from previous call to Cholesky)

  if(.not. A%fact_chol) then
     
     if (associated(A%e)) deallocate(A%e)
     if (associated(A%f)) deallocate(A%f)
     if (associated(A%s)) deallocate(A%s)
     
     allocate(A%f(size(A%d)))
     A%f(:) = 0.0_wp

     ! 4.5.1 Compute elements of full matrix
     
     call matrix_pp2full_alloc(A, A_full)
     
     ! 4.5.2 Check that matrix is a square matrix
     
     if(n /= size(A_full,2))then
        call message(msg_error, 'LHS matrix to be solved for ' // trim(routine) // ' is not square - aboriting calculation')
        A%d(:) = ropp_MDFV
        x(:,:) = ropp_MDFV
        RETURN
     endif
               
     ! 4.5.3 Determine Cholesky matrix G
     
     G = 0.0_wp 
     
     do i = 1, n
        
        temp(i:n) = A_full(i:n, i)      ! store lower triangle matrix
        
        if(i /= 1)then         
           do j = 1, i-1
              temp(i:n) = temp(i:n) - G(i,j)*G(i:n, j)
           enddo
        endif
        
        !  4.5.4 Check that matrix is positive definite
        
        if(temp(i) <= TINY(0.)*100.0_wp)then
           call message(msg_error, 'Matrix to be inverted in '// trim(routine) // ' is not positive definite - aborting')
           A%d(:) = ropp_MDFV
           x(:,:) = ropp_MDFV
           RETURN
        endif
        
        !  4.5.5 Compute elements of Cholesky matrix G

        G(i:n, i) = temp(i:n) / sqrt(temp(i))
        
     enddo
     
     ! 4.5.6 Store Cholesky matrix as packed type 
     call matrix_full2pp(G, A%f, 'L')
     
     ! 4.5.7 Set initialised flag in current matrix structure
     
     A % fact_chol = .true.
     deallocate(A_full)
     
  endif

! 4.7 Solve G.y = b for y by forward substitution

  do k = 1, size(b,2)

     y = b(k,:)
     y(1) = y(1) / A%f(1)
     do j = 2, n 

        do i=1,j-1
           temp(i) = A%f(j+(2*n-i)*(i-1)/2)
        enddo
     
        y(j) = (y(j) - dot_product( temp(1:j-1),y(1:j-1))) / A%f(j+(2*n-j)*(j-1)/2)
     enddo
     
! 4.8 Solve G^T.x = y for x by backward substitution

     x(:,k) = y
     x(n,k) = x(n,k) / A%f(n+n*(n-1)/2)
     do j = n-1, 1, -1
        
        do i=j+1,n
           temp(i) = A%f(i+(2*n-j)*(j-1)/2)
        enddo
        
        x(j,k) = (x(j,k) - dot_product(temp(j+1:n),x(j+1:n,k))) / A%f(j+(2*n-j)*(j-1)/2)
     enddo
     
  enddo

! 4.9 Clean up

  call message_set_routine(routine)

end function matrix_solve_packed_2d

