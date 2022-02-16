!     Compute the number of streams,the area coverage of each stream and
!     the cloud distribution in each stream.
SUBROUTINE rttov_cldstr_tl( &
            & profiles,    &
            & profiles_tl, &
            & ircld,       &
            & ircld_tl)
!     Description:
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
!     Method:
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           30/7/2004    Marco Matricardi. ECMWF.
!     2           24/03/2006   Marco Matricardi. ECMWF.
!                              A bug has been fixed tha caused
!                              the computation of a wrong number of
!                              streams in some situation.Also, a new feature
!                              has been introduced that allow to reduce the
!                              number of streams by considering only those
!                              streams whose weight is larger than
!                              cldstr_threshold. cldstr_threshold is pecified
!                              in the module rttov_const. By setting
!                              cldstr_threshold to a negative number,
!                              all the streams will be processed.
!                              This feature is used with caution.
!                              In fact, the sum of the weights of all streams
!                              (including the clear one) must be equal to 1.
!                              If some streams are not considered,
!                              this means that the weight of the clear stream
!                              has to be adjusted. As a consequence, if a too
!                              large value is used for cldstr_threshold,
!                              this can have serious implications for the
!                              accuracy of the results. In conclusion, only
!                              very small values should be used for
!                              cldstr_threshold just to remove the streams
!                              with a very small weight.
!     3           15/07/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     4           02/12/2009   All variables stored as layer arrays (Marco Matricardi)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY : profile_Type, ircld_type
!INTF_OFF
  USE rttov_const, ONLY : ncldtyp
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!       Arguments with intent in:
  TYPE(profile_type), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type), INTENT(IN)    :: profiles_tl(size(profiles))
!       Arguments with intent inout:
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld_tl
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: I, J, ISTR, IJSTR, ILAY
  INTEGER(KIND=jpim) :: nprofiles                    ! Number of profiles
!       Local tangent linear arrays:
  REAL   (KIND=jprb) :: NTOT  (profiles(1)%nlevels)
!       Local direct arrays arrays:
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!      ircld_tl(:)%XSTRCLR=0._jprb
!---------Compute number of streams and cloud distribution in each stream-------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(profiles)
  DO J = 1, NPROFILES
    ircld_tl%CLDCFR(:, j)  = 0._jprb
    ircld_tl%XSTRMIN(:, j) = 0._jprb
    ircld_tl%XSTRMAX(:, j) = 0._jprb
    ntot = 0._jprb
    DO ilay = 1, profiles(1)%nlayers
      DO i = 1, ncldtyp
        IF (profiles(j)%cfrac(I, ilay) /= 0._jprb) THEN
          ircld_tl%cldcfr(ilay, J) = profiles_tl(j)%cfrac(i, ilay)
        ENDIF
      ENDDO
    ENDDO
    ircld%ICOUNT(J) = 0_jpim
!---------Compute the cumulative cloud coverage usin the maximum-random---------
!         overlap assumption
    NTOT(1)         =  - ircld_tl%CLDCFR(1, J)
    DO ILAY = 2, profiles(1)%nlayers
      IF ((ircld%CLDCFR(ILAY - 1, J) > ircld%CLDCFR(ILAY, J))) THEN
        ircld_tl%MAXCOV(ILAY, J) =  - ircld_tl%CLDCFR(ILAY - 1, J)
      ELSE IF ((ircld%CLDCFR(ILAY - 1, J) < ircld%CLDCFR(ILAY, J))) THEN
        ircld_tl%MAXCOV(ILAY, J) =  - ircld_tl%CLDCFR(ILAY, J)
      ELSE IF ((ircld%CLDCFR(ILAY - 1, J) == ircld%CLDCFR(ILAY, J))) THEN
        ircld_tl%MAXCOV(ILAY, J) =  - ircld_tl%CLDCFR(ILAY, J)
      ENDIF
      IF (ircld%NTOTREF(ILAY - 1, J) == 0._jprb .AND. (1. - ircld%CLDCFR(ILAY - 1, J)) == 0._jprb) THEN
        NTOT(ILAY) = ircld_tl%MAXCOV(ILAY, J)
      ELSE
        NTOT(ILAY) = NTOT(ILAY - 1) * ircld%MAXCOV(ILAY, J) / (1._jprb - ircld%CLDCFR(ILAY - 1, J)) +      &
          & ircld_tl%CLDCFR(ILAY - 1, J) * ircld%NTOTREF(ILAY - 1, J) * ircld%MAXCOV(ILAY, J) /            &
          & (1._jprb - ircld%CLDCFR(ILAY - 1, J)) ** 2 +                                                   &
          & ircld_tl%MAXCOV(ILAY, J) * ircld%NTOTREF(ILAY - 1, J) / (1._jprb - ircld%CLDCFR(ILAY - 1, J))
      ENDIF
    ENDDO
    DO ILAY = 1, profiles(1)%nlayers
      IF (ircld%CLDCFR(ILAY, J) /= 0._jprb) THEN
        NTOT(ILAY)   =  - NTOT(ILAY)
      ELSE
        NTOT(ILAY)   = 0._jprb
      ENDIF
    ENDDO
!---------Determine the limits of each stream----------------------------------
    ircld%ICOUNT(J) = 1_jpim
    DO ILAY = 1, profiles(1)%nlayers
      ircld_tl%XSTRMAX(ILAY, J) = NTOT(ILAY)
      ircld_tl%XSTRMIN(ILAY, J) = NTOT(ILAY) - ircld_tl%CLDCFR(ILAY, J)
      IF (ircld%XSTRMINREF(ILAY, J) < 0._JPRB) ircld_tl%XSTRMIN(ILAY, J) = 0._JPRB
      IF (ircld%XSTRMAX(ILAY, J) /= 0._JPRB) THEN
        ircld_tl%XSTR(ircld%ICOUNT(J), J) = ircld_tl%XSTRMIN(ILAY, J)
        ircld%ICOUNT(J)                   = ircld%ICOUNT(J) + 1
        ircld_tl%XSTR(ircld%ICOUNT(J), J) = ircld_tl%XSTRMAX(ILAY, J)
        ircld%ICOUNT(J)                   = ircld%ICOUNT(J) + 1
      ENDIF
    ENDDO
!---------Re-arrange the limits of each stream in ascending order---------------
    LOOP1 : DO ISTR = 2, ircld%ICOUNT(J) - 1
      ircld%A(ISTR, J)    = ircld%XSTRREF(ISTR, J)
      ircld_tl%A(ISTR, J) = ircld_tl%XSTR(ISTR, J)
      DO I = ISTR - 1, 1,  - 1
        IF (ircld%XSTRREF(I, J) <= ircld%A(ISTR, J)) THEN
          ircld%XSTRREF(I + 1, J) = ircld%A(ISTR, J)
          ircld_tl%XSTR(I + 1, J) = ircld_tl%A(ISTR, J)
          CYCLE LOOP1
        ELSE
          ircld%XSTRREF(I + 1, J) = ircld%XSTRREF(I, J)
          ircld_tl%XSTR(I + 1, J) = ircld_tl%XSTR(I, J)
        ENDIF
      ENDDO
      I = 0
      ircld%XSTRREF(I + 1, J) = ircld%A(ISTR, J)
      ircld_tl%XSTR(I + 1, J) = ircld_tl%A(ISTR, J)
    ENDDO LOOP1
    OUTER : DO I = 1, ircld%ILOOP(j)
      INNER : DO ISTR = 1, ircld%ILOOPIN(I, j)
        IF (ircld%XSTRREF2(I, ISTR, j) == ircld%XSTRREF2(I, ISTR + 1, j)) THEN
          IF (ircld%XSTRREF2(I, ISTR, j) /= 1._jprb) THEN
            ircld%ICOUNT1REF(I, j) = ircld%ICOUNT1REF(I, j) - 1
            DO IJSTR = ISTR, ircld%ICOUNT1REF(I, j)
              ircld_tl%XSTR(IJSTR, j) = ircld_tl%XSTR(IJSTR + 1, j)
            ENDDO
          ENDIF
        ENDIF
      ENDDO INNER
    ENDDO OUTER
!---------Comppute the weight of the clear stream------------------------------
    ircld_tl%XSTRCLR(J) =  - ircld_tl%XSTR(ircld%NSTREAMREF(J) + 1, J) + ircld_tl%XSTR(1, J)
!---------Consider only the streams whose weight is greater than cldstr_threshold-------
    IF (ircld%NSTREAMREF(J) /= 0_jpim) THEN
      IF (ircld%ICOUNSTR(J) /= 0_jpim) THEN
        DO ISTR = 1, ircld%ICOUNSTR(J) + 1
          IF (ISTR == 1) THEN
            ircld_tl%XSTR(ISTR, J) = ircld_tl%XSTR(ircld%INDEXSTR(ISTR, J), J)
          ELSE
            ircld_tl%XSTR(ISTR, J) = ircld_tl%XSTR(ISTR - 1, J) +      &
              & (ircld_tl%XSTR(ircld%INDEXSTR(ISTR - 1, J) + 1, J) - ircld_tl%XSTR(ircld%INDEXSTR(ISTR - 1, J), J))
          ENDIF
        ENDDO
        ircld_tl%XSTRCLR(J) =  - ircld_tl%XSTR(ircld%ICOUNSTR(J) + 1, J) + ircld_tl%XSTR(1, J)
      ELSE
        ircld_tl%XSTRCLR(J) = 0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_cldstr_tl
