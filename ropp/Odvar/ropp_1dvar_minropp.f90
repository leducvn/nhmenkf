! $Id: ropp_1dvar_minropp.f90 3551 2013-02-25 09:51:28Z idculv $

!****p* 1DVar/ropp_1dvar_minropp *
!
! NAME
!    ropp_1dvar_minropp - Minimisation routine using Quasi-Newton method
!
! SYNOPSIS
!    ropp_1dvar_minropp(state, J_grad, J_dir, dJ, gconv, niter, indic, 
!                       miter, maxstore)
! 
! DESCRIPTION
!
! INPUTS
!    state      - control vector (1st guess)
!    J_grad     - penalty function gradient
!    J_dir      - quasi-Newton method direction
!    dJ         - expected decrease of cost function
!    gconv      - decrease in norm of penalty function for convergence
!    niter      - iteration counter
!    indic      - communication flag
!    miter      - maximum number of iterations allowed
!    maxstore   - size of storage available
!
! OUTPUT
!    state      - control vector solution
!    J_grad     - penalty function gradient
!    J_dir      - quasi-Newton method direction
!    niter      - iteration counter
!    indic      - communication flag
!
! NOTES
!       Based on limited-memory quasi-Newton method with diagonal scaling.
!       Further details on algorithm methods in  
!       Gilbert, J. C., Lemarechal, C. (1989) Some numerical experiments with 
!          variable storage quasi-Newton algorithms. Mathematical Programming, 
!          45, 407-435  
!      Nocedal, J. (1980) Updating quasi-Newton matrices with limited
!          storage. Mathematics of Computation, 35, 773-782
!      Based on original code written by Don Allan.
!
!      Cost function and gradient must be computed external to this routine.
!      Communication with external routines achieved using indic flag
!            indic = 0   : convergence achieved, no further minimisation
!            indic = 1   : first call in loop
!            indic = 4   : calling program has to compute new values for J and 
!                          gradient of J based for updated control vector
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

SUBROUTINE ropp_1dvar_minropp(x, g, p, dJ, gconv, niter, indic, miter, maxstore)
 
!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------
   
  USE typesizes, ONLY: wp => EightByteReal
  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:), INTENT(inout) :: x         ! Control vector 
                                                     ! IN=guess OUT=solution
  REAL(wp), DIMENSION(:), INTENT(inout) :: g         ! Penalty func gradient
  REAL(wp), DIMENSION(:), INTENT(inout) :: p         ! Quasi-Newton direction 
  REAL(wp),               INTENT(in)    :: dJ        ! Expected decrease of cost
  REAL(wp),               INTENT(in)    :: gconv     ! Decrease in |g| for conv
  INTEGER,                INTENT(inout) :: niter     ! Iteration counter
  INTEGER,                INTENT(inout) :: indic     ! Communication flag
  INTEGER,                INTENT(in)    :: miter     ! Maximum no. of iterations
  INTEGER,                INTENT(in)    :: maxstore  ! Size of storage available

  REAL(wp), DIMENSION(:),   ALLOCATABLE, SAVE :: xold  ! Control vectors
  REAL(wp), DIMENSION(:),   ALLOCATABLE, SAVE :: gold  ! Gradient of cost
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: D     ! Diagonal scaling vector
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: xnew  ! Updated control
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: store ! QN s, y  vectors
  REAL(wp)                              :: gnbar ! Norm of the gradient  
  INTEGER                               :: nuke  ! Update to be (over)written
  INTEGER                               :: icg   ! Nocedals method: =0 intially,
                                                 !     =1 when overwriting store
  INTEGER                               :: con   ! =0 initially, 
                                                 ! =1 when convergence achieved
  INTEGER                               :: restart ! =0 at restart, =1 otherwise
  INTEGER                               :: n     ! Size of control vector

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  n = SIZE(x)
  con = 0
  restart = 0  
  IF(indic == 1)THEN      ! First call of routine in minimisation loop
     niter = 1
     
     IF (ALLOCATED(xold)) DEALLOCATE(xold)
     ALLOCATE(xold(n))
     xold(:) = 0.0_wp
     
     IF (ALLOCATED(gold)) DEALLOCATE(gold)
     ALLOCATE(gold(n))
     gold(:) = 0.0_wp

  ENDIF

  ALLOCATE(D(n))
  ALLOCATE(xnew(n))
  ALLOCATE(store(n,maxstore*2))

!-------------------------------------------------------------------------------
! 3. Initial convergence checks
!-------------------------------------------------------------------------------

  ! Convergence achieved in external routine?
  IF(indic==0)THEN
     DEALLOCATE(xold)
     DEALLOCATE(gold)
     DEALLOCATE(D)
     DEALLOCATE(xnew)
     DEALLOCATE(store)
     RETURN
  ENDIF
  
  ! Gradient convergence check
  gnbar = dsqrt(DOT_PRODUCT(g,g))     ! Find norm(g)
  IF(gnbar < gconv)THEN
     ! Exit if convergence achieved
     con=1
     indic = 0
     DEALLOCATE(xold)
     DEALLOCATE(gold)
     DEALLOCATE(D)
     DEALLOCATE(xnew)
     DEALLOCATE(store)
     RETURN
  ENDIF

!-------------------------------------------------------------------------------
! 4. Main loop
!-------------------------------------------------------------------------------

  DO WHILE (con == 0 .AND. niter <= miter)     

     DO

        IF (niter >= miter) THEN
          indic = 0
          RETURN
        ENDIF
        
! 4.1 Clear store on restart and set initial search direction
! -----------------------------------------------------------
        
        IF (restart == 0) THEN
           nuke = 1
           store = 0.0_wp
           icg = 0
           
           ! Find initial value for direction vector using Fletcher's scaling
           IF(indic == 1)THEN
              gnbar = dsqrt(DOT_PRODUCT(g,g))
              p = -g * (2.0_wp * dJ/gnbar**2)   
           ENDIF
        ENDIF

! 4.2 Compute new value for state vector
! -------------------------------------- 

        CALL ropp_1dvar_minropp_linesearch(p, x, xnew, indic)

! 4.3 Update x,g vectors and store (s and y)
! ------------------------------------------

        IF(indic /= 1)THEN

           niter = niter + 1
           
           xold = x
           x = xnew     
           gold = g

           ! Cleanup if convergence achieved
           IF(indic == 0)THEN
              DEALLOCATE(xold)
              DEALLOCATE(gold)
           ENDIF

           DEALLOCATE(D)
           DEALLOCATE(xnew)
           DEALLOCATE(store)
           
           ! Exit routine if new x value found ok
           RETURN
           
        ENDIF
        
        store(:,nuke) = x - xold
        store(:,nuke+maxstore) = g - gold 

! 4.4 Generate new value of descent direction p using Nocedal's Method 
! --------------------------------------------------------------------        

        ! Initialise diagonal scaling vector
        IF (restart == 0) THEN    
           D = 1.0_wp     
        ENDIF
        
        ! Find diagonal scaling vector D from s and y 
        CALL ropp_1dvar_minropp_dscale(store(:,nuke), store(:,nuke+maxstore), &
                                       D, n) 
        
        ! Find search direction p(n) used in quasi-Newton method

        CALL ropp_1dvar_minropp_hmult(g,store,D,nuke,icg,p,n,maxstore)
        p = - p

        ! If max storage reached write over oldest s and y vectors
        IF (nuke == maxstore) THEN
           nuke = 1
           icg = 1
        ELSE
           nuke = nuke + 1
           restart = 1
        ENDIF
                        
     END DO
     
  END DO

  DEALLOCATE(D)
  DEALLOCATE(xnew)
  DEALLOCATE(store)

CONTAINS

!-------------------------------------------------------------------------------
! 5. Find search direction p(n) used in quasi-Newton method
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_minropp_hmult(g,store,D,nuke,icg,p,n,maxstore)

! 5.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  IMPLICIT NONE

  INTEGER,                  INTENT(in)  :: n         ! Dimension
  INTEGER,                  INTENT(in)  :: maxstore  ! Size of storage available
  REAL(wp), DIMENSION(:),   INTENT(in)  :: g     ! Penalty function gradient
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: store ! Quasi-Newton s, y vectors
  REAL(wp), DIMENSION(:),   INTENT(in)  :: D     ! Diagonal scaling vector
  INTEGER,                  INTENT(in)  :: nuke  ! No. of update to (over)write
  INTEGER,                  INTENT(in)  :: icg   ! =1 when overwriting begins
  REAL(wp), DIMENSION(:),   INTENT(out) :: p     ! QN direction vector
  
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: alpha
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: q
  REAL(wp)                              :: rho
  REAL(wp)                              :: beta
  INTEGER                               :: i,k       ! Loop counters

! 5.2 Initialise variables
! ------------------------

  ALLOCATE(q(n))
  ALLOCATE(alpha(maxstore))

  alpha = 0.0_wp
  q = g
  i = nuke

! 5.3 Compute first matrix product in BFGS method 
! -----------------------------------------------
  
  DO k = 1,maxstore

     ! Find rho, alpha and q
     rho = 1.0_wp / DOT_PRODUCT(store(:,i+maxstore),store(:,i))
     alpha(i) = rho * DOT_PRODUCT(store(:,i),q)
     q = q - alpha(i) * store(:,i+maxstore)
     
     i = i - 1
     
     ! If store has started overwriting reset counter to maxstore   
     ! otherwise exit loop     
     IF (i == 0 .AND. icg == 1) THEN
        i = maxstore
     ENDIF
     
     IF (i == 0) EXIT
     
  END DO
  
! 5.4 Add diagonal scaling
! ------------------------

  DO i=1,n
     p(i) = D(i)*q(i)
  END DO

! 5.5 Check if store has started overwriting
! ------------------------------------------

  ! If store has started overwriting set counter to start from i=nuke+1
  ! otherwise start from i=1  
  IF (icg == 0 .OR. nuke == maxstore) THEN
     i = 1
  ELSE
     i = nuke + 1
  ENDIF
  
! 5.6 Compute second matrix product in BFGS method
! ------------------------------------------------
  
  DO k = 1,maxstore

     rho = 1.0_wp / DOT_PRODUCT(store(:,i+maxstore),store(:,i))
     beta = rho * DOT_PRODUCT(store(:,i+maxstore),p)
     p = p + (store(:,i) * (alpha(i) - beta))
     
     ! if store has not started overwriting exit loop. Otherwise reset 
     ! loop when i=maxstore 
     
     IF (i == nuke) EXIT   
     
     IF (i == maxstore) THEN
        i = 1
     ELSE
        i = i + 1
     ENDIF
     
  END DO

  DEALLOCATE(q)
  DEALLOCATE(alpha)
  
END SUBROUTINE ropp_1dvar_minropp_hmult

!-------------------------------------------------------------------------------
! 6. Find diagonal scaling vector used in quasi-Newton method
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_minropp_dscale(s,y,D,n)
  
! 6.1 Declarations
! ----------------
 
  USE typesizes, ONLY: wp => EightByteReal
  IMPLICIT NONE

  INTEGER,                INTENT(in)  :: n     ! Dimension
  REAL(wp), DIMENSION(:), INTENT(in)  :: s     ! Quasi-Newton method s vector
  REAL(wp), DIMENSION(:), INTENT(in)  :: y     ! Quasi-Newton method y vector
  REAL(wp), DIMENSION(:), INTENT(out) :: D     ! Diagonal scaling vector

  REAL(wp), DIMENSION(:), ALLOCATABLE :: Dinvs
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Dy
  REAL(wp)                            :: elk
  REAL(wp)                            :: enu
  REAL(wp)                            :: edd
  INTEGER                             :: i 

  ALLOCATE(Dinvs(n))
  ALLOCATE(Dy(n))

! 6.2 Find D^-1*s
! ---------------

  DO i=1,n  
     Dinvs(i)=s(i)/D(i)
  END DO
  
! 6.3 Find <D^-1*s,s>
! -------------------

  elk = DOT_PRODUCT(Dinvs,s)

! 6.4 Find <y,s>
! --------------

  enu = DOT_PRODUCT(y,s)

! 6.5 Find D*y
! ------------

  DO i=1,n
     Dy(i) = D(i)*y(i) 
  END DO

! 6.6 Find <D*y,y>
! ---------------- 

  edd = DOT_PRODUCT(Dy,y)

! 6.7 Compute D
! -------------

  DO i=1,n
     
     D(i) = 1.0_wp / &
          ( edd/(enu*D(i)) + y(i)**2/enu - (edd*s(i)**2)/(enu*elk*D(i)**2) )
     
  END DO

  DEALLOCATE(Dinvs)
  DEALLOCATE(Dy)
  
END SUBROUTINE ropp_1dvar_minropp_dscale

!-------------------------------------------------------------------------------
! 7. Compute new value for control vector x
!-------------------------------------------------------------------------------
!
!      Communication with external routines achieved using indic flag
!            indic = 0   : convergence achieved, no further minimisation
!            indic = 1   : first call in loop
!            indic = 4   : calling program has to compute new values for J and 
!                          gradient of J based for updated control vector

SUBROUTINE ropp_1dvar_minropp_linesearch(p,x,xnew,indic)

! 7.1 Declarations
! ----------------
  
  USE typesizes, ONLY: wp => EightByteReal
  IMPLICIT NONE
  
  REAL(wp), DIMENSION(:), INTENT(in)    :: p      ! Search direction vector
  REAL(wp), DIMENSION(:), INTENT(in)    :: x      ! Initial control vector
  REAL(wp), DIMENSION(:), INTENT(out)   :: xnew   ! New control vector 
  INTEGER,                INTENT(inout) :: indic  ! Communication flag

  REAL(wp), PARAMETER :: t = 1.0_wp  ! step length

! 7.2 If convergence already achieved, nothing to be done
! -------------------------------------------------------

  IF (indic == 0) THEN     

     RETURN 

! 7.3 Compute new value of control vector and set flag to find new J, J_grad
! --------------------------------------------------------------------------
  ELSE IF (indic == 1) THEN     
     xnew = x + t*p

    WHERE (ABS(t*p) > ABS(x/2.0))
      xnew = x
    ENDWHERE

     indic = 4            
     RETURN
 
! 7.4 Set flag to estimate new search direction vector p
! ------------------------------------------------------
  ELSE

     indic=1
     
  ENDIF
    
END SUBROUTINE ropp_1dvar_minropp_linesearch
 
END SUBROUTINE ropp_1dvar_minropp
