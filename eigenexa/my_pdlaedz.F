!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: my_pdlaedz_mod
!
! Copyright(C) 2012-2018 RIKEN.
! Copyright(C) 2011-2012 Toshiyuki Imamura
!                        Graduate School of Informatics and Engineering,
!                        The University of Electro-Communications.
! Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
! 
! Redistribution  and  use  in  source and binary forms, with or without
! modification,  are  permitted  provided  that the following conditions
! are met:
! 
! * Redistributions  of  source  code  must  retain  the above copyright
!   notice,  this  list  of  conditions  and  the  following  disclaimer.
! * Redistributions  in  binary  form must reproduce the above copyright
!   notice,  this list of conditions and the following disclaimer in the
!   documentation  and/or other materials provided with the distribution.
! 
! THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
! LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE MY_PDLAEDZ_MOD

      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MY_PDLAEDZ(ISTAGE, HBW, N, N1, ID, Q, IQ, JQ, LDQ, 
     &     DESCQ, Z, WORK, V )
!     
!     This routine is a modified routine of PDLAEDZ
!     which is contained in the ScaLAPACK
!     Modification is done by T.Imamura
!     
!     .. Scalar Arguments ..
      INTEGER            ISTAGE
      INTEGER            ID, IQ, JQ, LDQ, N, N1, HBW
!     ..
!     .. Array Arguments ..
      INTEGER            DESCQ( * )
      DOUBLE PRECISION   Q( LDQ, * ), WORK( * ), Z( * ),
     &     V( * )
!     ..
!     
!     Purpose
!     =======
!     
!     MY_PDLAEDZ Form the z-vector which consists of the last row of Q_1
!     and the first row of Q_2.
!     =====================================================================
!     
!     .. Parameters ..
!     
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     &     MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
!     
      INTEGER            COL, I, IBUF, ICTXT, IIQ, IIZ1, IIZ2, IQCOL,
     &     IQROW, IZ, IZ1, IZ1COL, IZ1ROW, IZ2, IZ2COL,
     &     IZ2ROW, J, JJQ, JJZ1, JJZ2, MYCOL, MYROW, N2,
     &     NB, NBLOC, NPCOL, NPROW, NQ1, NQ2, ZSIZ, 
     &     K
      
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     &     DGERV2D, DGESD2D, INFOG2L, DAXPY
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
!     ..
      integer :: ERR_MPI
!     .. Executable Statements ..
!     

!     This is just to keep ftnchek and toolpack/1 happy
      

      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     &     RSRC_.LT.0 )RETURN
!     
      ICTXT = DESCQ( CTXT_ )
      NB = DESCQ( NB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( ID, ID, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     &     IQROW, IQCOL )
      N2 = N - N1

!     
!$OMP DO
      DO J = 1, N
         WORK(J) = ZERO
      END DO
!$OMP END DO
!     
      IF (ISTAGE == 1) THEN
!     
         DO K = 1, HBW
            CALL INFOG2L( IQ-1+(ID+N1-1-HBW+K), JQ-1+ID, DESCQ, NPROW, 
     &           NPCOL, MYROW, MYCOL, IIZ1, JJZ1, IZ1ROW, IZ1COL )
            NQ1 = NUMROC( N1, NB, MYCOL, IZ1COL, NPCOL )
            IF( ( MYROW.EQ.IZ1ROW ) .AND. ( NQ1.NE.0 ) ) THEN
               I = MOD( MYCOL + NPCOL - MOD( IZ1COL, NPCOL ), NPCOL )
               NBLOC = ( NQ1-1 ) / NB + 1
!$OMP DO
               DO J = 1, NBLOC
                  IZ1 = (J-1)*NB
                  IZ = I*NB + 1 + (J-1)*NB*NPCOL
                  ZSIZ = MIN( NB, NQ1-IZ1 )
                  IF ( ZSIZ > 0 ) THEN
                     IF ( K .EQ. 1 ) THEN
!     CALL DCOPY( ZSIZ, Q(IIZ1,JJZ1+IZ1), LDQ,
!     &                       WORK(IZ), 1)
!     CALL DSCAL( ZSIZ, V(K), WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       V(K) * Q(IIZ1,JJZ1+IZ1:JJZ1+IZ1+ZSIZ-1)
                     ELSE
!     CALL DAXPY( ZSIZ, V(K),
!     &                       Q(IIZ1,JJZ1+IZ1), LDQ, WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       WORK(IZ:IZ+ZSIZ-1) +
     &                       V(K) * Q(IIZ1,JJZ1+IZ1:JJZ1+IZ1+ZSIZ-1)
                     END IF
                  END IF
               END DO
!$OMP END DO
            END IF
         END DO
!     
         DO K = 1, HBW
            CALL INFOG2L( IQ-1+(ID+N1)+K-1, JQ-1+(ID+N1), DESCQ, 
     &           NPROW, NPCOL, MYROW, MYCOL, IIZ2, JJZ2, IZ2ROW, IZ2COL)
            NQ2 = NUMROC( N2, NB, MYCOL, IZ2COL, NPCOL )
            IF( ( MYROW.EQ.IZ2ROW ) .AND. ( NQ2.NE.0 ) ) THEN
               I = MOD( MYCOL + NPCOL - MOD( IZ2COL, NPCOL ), NPCOL )
               NBLOC = ( NQ2-1 ) / NB + 1
!$OMP DO
               DO J = 1, NBLOC
                  IZ2 = (J-1)*NB
                  IZ = NB*I + N1 + 1 + (J-1)*NB*NPCOL
                  ZSIZ = MIN( NB, NQ2-IZ2 )
                  IF ( ZSIZ > 0 ) THEN
                     IF ( K .EQ. 1 ) THEN
!     CALL DCOPY( ZSIZ, Q(IIZ2,JJZ2+IZ2), LDQ,
!     &                       WORK(IZ), 1)
!     CALL DSCAL( ZSIZ, V(K+HBW), WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       V(K+HBW) * Q(IIZ2,JJZ2+IZ2:JJZ2+IZ2+ZSIZ-1)
                     ELSE
!     CALL DAXPY( ZSIZ, V(K+HBW),
!     &                       Q(IIZ2,JJZ2+IZ2), LDQ, WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       WORK(IZ:IZ+ZSIZ-1) +
     &                       V(K+HBW) * Q(IIZ2,JJZ2+IZ2:JJZ2+IZ2+ZSIZ-1)
                     END IF
                  END IF
               END DO
!$OMP END DO
            END IF
         END DO
!     
      ELSE
!     
         DO K = 1, HBW+HBW
            CALL INFOG2L(IQ-1+(ID+N1-1-HBW+K), JQ-1+ID, DESCQ, NPROW,
     &           NPCOL, MYROW, MYCOL, IIZ2, JJZ2, IZ2ROW, IZ2COL )
            NQ1 = NUMROC( N, NB, MYCOL, IZ2COL, NPCOL )
            IF( ( MYROW.EQ.IZ2ROW ) .AND. ( NQ1.NE.0 ) ) THEN
               I = MOD( MYCOL + NPCOL - MOD( IZ2COL, NPCOL ), NPCOL )
               NBLOC = ( NQ1-1 ) / NB + 1
!$OMP DO
               DO J = 1, NBLOC
                  IZ1 = (J-1)*NB
                  IZ = I*NB + 1 + (J-1)*NB*NPCOL
                  ZSIZ = MIN( NB, NQ1-IZ1 )
                  IF ( ZSIZ > 0 ) THEN
                     IF ( K .EQ. 1 ) THEN
!     CALL DCOPY( ZSIZ, Q(IIZ2,JJZ2+IZ1), LDQ,
!     &                       WORK(IZ), 1)
!     CALL DSCAL( ZSIZ, V(K), WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       V(K) * Q(IIZ2,JJZ2+IZ1:JJZ2+IZ1+ZSIZ-1)
                     ELSE
!     CALL DAXPY( ZSIZ, V(K),
!     &                       Q(IIZ2,JJZ2+IZ1), LDQ, WORK(IZ), 1)
                        WORK(IZ:IZ+ZSIZ-1) =
     &                       WORK(IZ:IZ+ZSIZ-1) +
     &                       V(K) * Q(IIZ2,JJZ2+IZ1:JJZ2+IZ1+ZSIZ-1)
                     ENDIF
                  END IF
               END DO
!$OMP END DO
            END IF
         END DO
!     
      END IF
!     
!$OMP BARRIER
!     
!     HERE, WE NEED FLOATING-SUMUP TO RETRIVE THE VECTOR WORK ONTO Z
!$OMP MASTER
      call MPI_ALLREDUCE( WORK( 1 ), Z( 1 ),
     &     N, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, TRD_COMM_WORLD, ERR_MPI )
!$OMP END MASTER
!     
!$OMP BARRIER
!     
      RETURN
!     

!     End of MY_PDLAEDZ
!     
      END SUBROUTINE MY_PDLAEDZ

      END MODULE MY_PDLAEDZ_MOD
