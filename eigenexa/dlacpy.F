!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: dlacpy_mod
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

      MODULE DLACPY_MOD

      IMPLICIT NONE
      PUBLIC

      CONTAINS

!     
!     ** OpenMP thread parallelized version **
!     
!     Definition:
!     ===========
!     
!     SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!     
!     .. Scalar Arguments ..
!     CHARACTER          UPLO
!     INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
!     DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!     
!     
!     > \par Purpose:
!     =============
!     >
!     > \verbatim
!     >
!     > DLACPY copies all or part of a two-dimensional matrix A to another
!     > matrix B.
!     > \endverbatim
!     
!     Arguments:
!     ==========
!     
!     > \param[in] UPLO
!     > \verbatim
!     >          UPLO is CHARACTER*1
!     >          Specifies the part of the matrix A to be copied to B.
!     >          = 'U':      Upper triangular part
!     >          = 'L':      Lower triangular part
!     >          Otherwise:  All of the matrix A
!     > \endverbatim
!     >
!     > \param[in] M
!     > \verbatim
!     >          M is INTEGER
!     >          The number of rows of the matrix A.  M >= 0.
!     > \endverbatim
!     >
!     > \param[in] N
!     > \verbatim
!     >          N is INTEGER
!     >          The number of columns of the matrix A.  N >= 0.
!     > \endverbatim
!     >
!     > \param[in] A
!     > \verbatim
!     >          A is DOUBLE PRECISION array, dimension (LDA,N)
!     >          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!     >          or trapezoid is accessed; if UPLO = 'L', only the lower
!     >          triangle or trapezoid is accessed.
!     > \endverbatim
!     >
!     > \param[in] LDA
!     > \verbatim
!     >          LDA is INTEGER
!     >          The leading dimension of the array A.  LDA >= max(1,M).
!     > \endverbatim
!     >
!     > \param[out] B
!     > \verbatim
!     >          B is DOUBLE PRECISION array, dimension (LDB,N)
!     >          On exit, B = A in the locations specified by UPLO.
!     > \endverbatim
!     >
!     > \param[in] LDB
!     > \verbatim
!     >          LDB is INTEGER
!     >          The leading dimension of the array B.  LDB >= max(1,M).
!     > \endverbatim
!     
!     =====================================================================
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!     
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!     
!     =====================================================================
!     
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!     
      IF( LSAME( UPLO, 'U' ) ) THEN
!$OMP PARALLEL DO
!$OMP+         PRIVATE(I,J)
!$OMP+         SCHEDULE(STATIC,1)
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
 10         CONTINUE
 20      CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!$OMP PARALLEL DO
!$OMP+         PRIVATE(I,J)
!$OMP+         SCHEDULE(STATIC,1)
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
 30         CONTINUE
 40      CONTINUE
      ELSE
!$OMP PARALLEL DO
!$OMP+         PRIVATE(I,J)
!$OMP+         SCHEDULE(STATIC,1)
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
 50         CONTINUE
 60      CONTINUE
      END IF
      RETURN
!     
!     End of DLACPY
!     
      END SUBROUTINE DLACPY

      END MODULE DLACPY_MOD
