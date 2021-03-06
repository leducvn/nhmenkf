!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: my_pdsxedc_mod
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

!     This subroutine  is modified from scalapack's PDSTEDC to calculate eigenvalues
!     and eigenvectors of a pentadiagonal matrix in parallel using divide and conquer
!     
!     This routine is a modified routine of PSTEDC
!     which is contained in the ScaLAPACK
!     Modification is done by T.Imamura and P.H.Pham
!     
      MODULE MY_PDSxEDC_MOD

      USE MY_PDLASRT_MOD
      USE MY_PDLAED0_MOD
      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MY_PDSxEDC( COMPZ, HBW, N, D, E, LDE, Q, IQ, JQ,  
     &     DESCQ, WORK, LWORK,IWORK, LIWORK, INFO )
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, IQ, JQ, LIWORK, LWORK, N, HBW, LDE
!     ..
!     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
      DOUBLE PRECISION   E( LDE, HBW )
!     ..
!     
!     =====================================================================
!     Purpose
!     =======
!     MY_PDSxEDC computes all eigenvalues and eigenvectors of a
!     symmetric pentadiagonal matrix in parallel, using the divide and
!     conquer algorithm.
!     
!     This code makes very mild assumptions about floating point
!     arithmetic. It will work on machines with a guard digit in
!     add/subtract, or on those binary machines without guard digits
!     which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!     It could conceivably fail on hexadecimal or decimal machines
!     without guard digits, but we know of none.  See DLAED3 for details.
!     
!     Arguments
!     =========
!     
!     COMPZ   (input) CHARACTER*1
!     = 'N':  Compute eigenvalues only.    (NOT IMPLEMENTED YET)
!     = 'I':  Compute eigenvectors of tridiagonal matrix also.
!     = 'V':  Compute eigenvectors of original dense symmetric
!     matrix also.  On entry, Z contains the orthogonal
!     matrix used to reduce the original matrix to
!     pentadiagonal form.            (NOT IMPLEMENTED YET)
!     
!     N       (global input) INTEGER
!     The order of the tridiagonal matrix T.  N >= 0.
!     
!     D       (global input/output) DOUBLE PRECISION array, dimension (N)
!     On entry, the diagonal elements of the pentadiagonal matrix.
!     On exit, if INFO = 0, the eigenvalues in descending order.
!     
!     E1      (global input/output) DOUBLE PRECISION array, dimension (N-1)
!     On entry, the subdiagonal elements of the pentadiagonal matrix.
!     On exit, E has been destroyed.
!     
!     
!     E2      (global input/output) DOUBLE PRECISION array, dimension (N-2)
!     On entry, the second subdiagonal elements of the pentadiagonal 
!     matrix. On exit, E has been destroyed.
!     

!     Q       (local output) DOUBLE PRECISION array,
!     local dimension ( LLD_Q, LOCc(JQ+N-1))
!     Q  contains the orthonormal eigenvectors of the symmetric
!     tridiagonal matrix.
!     On output, Q is distributed across the P processes in block
!     cyclic format.
!     
!     IQ      (global input) INTEGER
!     Q's global row index, which points to the beginning of the
!     submatrix which is to be operated on.
!     
!     JQ      (global input) INTEGER
!     Q's global column index, which points to the beginning of
!     the submatrix which is to be operated on.
!     
!     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
!     The array descriptor for the distributed matrix Z.
!     
!     
!     WORK    (local workspace/output) DOUBLE PRECISION array,
!     dimension (LWORK)
!     On output, WORK(1) returns the workspace needed.
!     
!     LWORK   (local input/output) INTEGER,
!     the dimension of the array WORK.
!     LWORK = 6*N + 2*NP*NQ
!     NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
!     NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
!     
!     If LWORK = -1, the LWORK is global input and a workspace
!     query is assumed; the routine only calculates the minimum
!     size for the WORK array.  The required workspace is returned
!     as the first element of WORK and no error message is issued
!     by PXERBLA.
!     
!     IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
!     On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
!     
!     LIWORK  (input) INTEGER
!     The dimension of the array IWORK.
!     LIWORK = 2 + 7*N + 8*NPCOL
!     
!     INFO    (global output) INTEGER
!     = 0:  successful exit
!     < 0:  If the i-th argument is an array and the j-entry had
!     an illegal value, then INFO = -(i*100+j), if the i-th
!     argument is a scalar and had an illegal value, then
!     INFO = -i.
!     > 0:  The algorithm failed to compute the INFO/(N+1) th
!     eigenvalue while working on the submatrix lying in
!     global rows and columns mod(INFO,N+1).
!     
!     =====================================================================
!     
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     &     MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     DOUBLE PRECISION   ZERO, ONE
!     PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICOFFQ, IIQ, IQCOL, IQROW, IROFFQ, JJQ,
     &     LDQ, LIWMIN, LWMIN, MYCOL, MYROW, NB, NP,
     &     NPCOL, NPROW, NQ, I, J
      DOUBLE PRECISION   ORGNRM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   DLANST
      EXTERNAL           INDXG2P, LSAME, NUMROC, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DLASCL, DSTEDC,
     &     INFOG2L, PXERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD, MAX
!     ..
!     .. Executable Statements ..
!     
!     This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     &     RSRC_.LT.0 )RETURN
!     
!     Test the input parameters.
!     
!     PRINT *, COMPZ, HBW, N, D(1:6) 
!     PRINT *, E(1:6,1), E(1:6,2), IQ, JQ
!     PRINT *, DESCQ(1:9), LWORK, LIWORK

      
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      LDQ = DESCQ( LLD_ )
      NB = DESCQ( NB_ )
!     PRINT *, "DESCQX==", DESCQ(1:9)     
!     PRINT *, "HERE", HBW, N, D(1), E(1,1)
      NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
      NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
!     PRINT *, "HERE", N, NB, MYROW, MYCOL
      INFO = 0
      
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IQ, JQ, DESCQ, 8, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            IROFFQ = MOD( IQ-1, DESCQ( MB_ ) )
            ICOFFQ = MOD( JQ-1, DESCQ( NB_ ) )
            IQROW = INDXG2P( IQ, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            IQCOL = INDXG2P( JQ, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = 6*N + 2*NP*NQ
            LIWMIN = 2 + 7*N + 8*NPCOL
!     
            WORK( 1 ) = DBLE( LWMIN )
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
            IF( .NOT.LSAME( COMPZ, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( IROFFQ.NE.ICOFFQ .OR. ICOFFQ.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -12
            END IF
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCQ( CTXT_ ), 'MY_PDSxEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF


!     
!     Quick return
!     
      IF( N.EQ.0 )
     &     GO TO 10
      CALL INFOG2L( IQ, JQ, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     &     IQROW, IQCOL )
      IF( N.EQ.1 ) THEN
         IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL )
     &        Q( 1 ) = ONE
         GO TO 10
      END IF
!     
!     If N is smaller than the minimum divide size NB, then
!     solve the problem with the serial divide and conquer
!     code locally.
!     Since currently dspedc is not implemented, below conditional
!     cases is not enable. Everything will be solved with pdspedc
!     
!     IF( N.LE.NB ) THEN
!     IF( ( MYROW.EQ.IQROW ) .AND. ( MYCOL.EQ.IQCOL ) ) THEN
!     IPQ = IIQ + ( JJQ-1 )*LDQ
!     CALL DSPEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK,
!     $                   IWORK, LIWORK, INFO )
!     IF( INFO.NE.0 ) THEN
!     INFO = ( N+1 ) + N
!     GO TO 10
!     END IF
!     END IF
!     GO TO 10
!     END IF
!     
!     If P=NPROW*NPCOL=1, solve the problem with DSTEDC.
!     
!     IF( NPCOL*NPROW.EQ.1 ) THEN
!     IPQ = IIQ + ( JJQ-1 )*LDQ
!     CALL DSTEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK, IWORK,
!     $                LIWORK, INFO )
!     GO TO 10
!     END IF
!     
!     Scale matrix to allowable range, if necessary.
!     

!     ORGNRM = DLANST( 'M', N, D, E )
      ORGNRM = ABS ( D(N) )
      DO I = 1, N - 1
         ORGNRM = MAX( ORGNRM, ABS( D(I) ) )
      END DO

      DO J = 1, HBW
         DO  I = 1, N - HBW
            ORGNRM = MAX( ORGNRM, ABS( E(I,J) ) )
         END DO
      END DO

      IF( ORGNRM.NE.ZERO ) THEN
         CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
         if ( N-1 >= 1 ) then
            DO 50 I = 1, HBW
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N-I, 1, E(1,I), 
     &              N-I, INFO )
               
 50         CONTINUE
         end if
      END IF
!     
      CALL MY_PDLAED0( N, HBW, D, E, LDE, Q, IQ, JQ, DESCQ, WORK, IWORK,
     &     NP, NQ, INFO )
!     
!     Sort eigenvalues and corresponding eigenvectors
!     
      CALL MY_PDLASRT( 'I', N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, IWORK,
     &     LIWORK, INFO )
!     
!     Scale back.
!     
      IF( ORGNRM.NE.ZERO )
     &     CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
!     
 10   CONTINUE
!     

      IF( LWORK.GT.0 )
     &     WORK( 1 ) = DBLE( LWMIN )
      IF( LIWORK.GT.0 )
     &     IWORK( 1 ) = LIWMIN
      RETURN

!     
!     End of MY_PDSxEDC
!     
      END SUBROUTINE MY_PDSxEDC

      END MODULE MY_PDSXEDC_MOD
