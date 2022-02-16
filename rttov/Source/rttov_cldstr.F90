!     Compute the number of streams,the area coverage of each stream and
!     the cloud distribution in each stream.
SUBROUTINE rttov_cldstr( &
            & profiles,         &
            & cldstr_threshold, &
            & ircld,            &
            & nstreams)
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
!     3           27/02/2009   Profile levels to include ToA. Cloud inputs on
!                              layers but stored as a level array in profiles.
!                              Distinguish between layer arrays and level
!                              arrays for index labels and looping (P. Rayer)
!     4           02/12/2009   All variables stored as layer arrays (Marco Matricardi)
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY : profile_Type, ircld_type
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : ncldtyp
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!       Arguments with intent out:
  INTEGER(KIND=jpim), INTENT(OUT)   :: nstreams
!       Arguments with intent in:
  TYPE(PROFILE_TYPE), INTENT(IN)    :: PROFILES(:)
  REAL(KIND=jprb),    INTENT(IN)    :: cldstr_threshold
!       Arguments with intent inout:
  TYPE(IRCLD_TYPE  ), INTENT(INOUT) :: IRCLD
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: I, J, ISTR, IJSTR, ILAY
  REAL   (KIND=jprb) :: delta_cfrac
!       Local arrays:
  INTEGER(KIND=jpim) :: nprofiles                        ! Number of profiles
  REAL   (KIND=jprb) :: NTOT(profiles(1)%nlevels)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!      ircld(:)%NSTREAM=0_jpim
!      ircld(:)%XSTRCLR=1._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(profiles)
  nstreams  = 0_jpim
  delta_cfrac = 10.0_jprb*epsilon(1.0_jprb)
  DO J = 1, NPROFILES
!---------Compute number of streams and cloud distribution in each stream-------
    ircld%iflag(:, j)      = 1_jpim
    ircld%ICLDARR(:, :, j) = 0_jpim
    ircld%CLDCFR(:, j)     = 0._jprb
    ircld%CLDTYP(:, :, j)  = 0._jprb
    ircld%XSTR(:, J)       = 0._jprb
    ircld%XSTRMIN(:, j)    = 0._jprb
    ircld%XSTRMAX(:, j)    = 0._jprb
    ircld%FLAG(:, j)       = .FALSE.
    NTOT = 0._jprb
    DO ILAY = 1, profiles(1)%nlayers
      DO I = 1, ncldtyp
        IF (profiles(j)%cfrac(I, Ilay) /= 0.) THEN
          ircld%CLDCFR(ILAY, J) = profiles(j)%cfrac(I, ILay)
        ENDIF
        IF (profiles(j)%cloud(I, Ilay) /= 0.) THEN
          ircld%CLDTYP(I, ILAY, J) = I
        ENDIF
      ENDDO
    ENDDO
    
    ! Check for overcast layers and identical cfrac values on consecutive layers
    ! These need adjusting to ensure TL/AD/K models are correct
    DO ILAY = 2, profiles(1)%nlayers
      IF (ircld%CLDCFR(ILAY, J) > 0.) THEN
        
        ! Check for overcast layers
        IF (ircld%CLDCFR(ILAY, J) == 1.0_jprb) THEN
          ircld%CLDCFR(ILAY, J) = ircld%CLDCFR(ILAY, J) - ilay*delta_cfrac
        ENDIF
        
        ! Check for identical adjacent cfrac (note that this won't always work if cldstr_threshold is +ve
        IF (ircld%CLDCFR(ILAY, J) == ircld%CLDCFR(ILAY-1, J)) THEN
          ircld%CLDCFR(ILAY, J) = ircld%CLDCFR(ILAY, J) - SIGN(delta_cfrac, ircld%CLDCFR(ILAY, J)-0.5_jprb)
        ENDIF
        
      ENDIF        
    ENDDO
    
    ircld%ICOUNT(J)     = 0_jpim
    ircld%ICOUNSTR(J)   = 0_jpim
!---------Compute the cumulative cloud coverage usin the maximum-random---------
!         overlap assumption
    NTOT(1)             = 1. - ircld%CLDCFR(1, J)
    ircld%NTOTREF(1, J) = NTOT(1)
    DO ILAY = 2, profiles(1)%nlayers
      ircld%MAXCOV(ILAY, J) = (1. - MAX(ircld%CLDCFR(ILAY - 1, J), ircld%CLDCFR(ILAY, J)))
      IF (NTOT(ILAY - 1) == 0 .AND. (1. - ircld%CLDCFR(ILAY - 1, J)) == 0._JPRB) THEN
        NTOT(ILAY)             = ircld%MAXCOV(ILAY, J)
        ircld%NTOTREF(ILAY, J) = NTOT(ILAY)
      ELSE
        NTOT(ILAY)             = NTOT(ILAY - 1) * ircld%MAXCOV(ILAY, J) / (1. - ircld%CLDCFR(ILAY - 1, J))
        ircld%NTOTREF(ILAY, J) = NTOT(ILAY)
      ENDIF
    ENDDO
    DO ILAY = 1, profiles(1)%nlayers
      IF (ircld%CLDCFR(ILAY, J) /= 0._JPRB) THEN
        NTOT(ILAY) = 1._JPRB - NTOT(ILAY)
      ELSE
        NTOT(ILAY) = 0._JPRB
      ENDIF
    ENDDO
!---------Determine the limits of each stream----------------------------------
    ircld%ICOUNT(J) = 1_jpim
    DO ILAY = 1, profiles(1)%nlayers
      ircld%XSTRMAX(ILAY, J)    = NTOT(ILAY)
      ircld%XSTRMIN(ILAY, J)    = NTOT(ILAY) - ircld%CLDCFR(ILAY, J)
      ircld%XSTRMINREF(ILAY, J) = ircld%XSTRMIN(ILAY, J)
      IF (ircld%XSTRMIN(ILAY, J) < 0._JPRB) ircld%XSTRMIN(ILAY, J) = 0._JPRB
      IF (ircld%XSTRMAX(ILAY, J) /= 0._JPRB) THEN
        ircld%XSTR(ircld%ICOUNT(J), J) = ircld%XSTRMIN(ILAY, J)
        ircld%ICOUNT(J)                = ircld%ICOUNT(J) + 1
        ircld%XSTR(ircld%ICOUNT(J), J) = ircld%XSTRMAX(ILAY, J)
        ircld%ICOUNT(J)                = ircld%ICOUNT(J) + 1
      ENDIF
    ENDDO
    ircld%XSTRREF(:, J) = ircld%XSTR(:, J)
!---------Re-arrange the limits of each stream in ascending order---------------
    LOOP1 : DO ISTR = 2, ircld%ICOUNT(J) - 1
      ircld%A(ISTR, J) = ircld%XSTR(ISTR, J)
      DO I = ISTR - 1, 1,  - 1
        IF (ircld%XSTR(I, J) <= ircld%A(ISTR, J)) THEN
          ircld%XSTR(I + 1, J) = ircld%A(ISTR, J)
          ircld%IFLAG(ISTR, J) = I
          ircld%FLAG(ISTR, J)  = .TRUE.
          CYCLE LOOP1
        ELSE
          ircld%XSTR(I + 1, J) = ircld%XSTR(I, J)
        ENDIF
      ENDDO
      I = 0
      ircld%XSTR(I + 1, J) = ircld%A(ISTR, J)
    ENDDO LOOP1
    ircld%ICOUNT1(J)    = ircld%ICOUNT(J) - 1
    ircld%ILOOP(J)      = 0_jpim
    ircld%ILOOPIN(:, J) = 0_jpim
    OUTER : DO
      ircld%ILOOP(J)                       = ircld%ILOOP(J) + 1
      ircld%ICOUNT1REF(ircld%ILOOP(J), J)  = ircld%ICOUNT1(J)
      ircld%XSTRREF2(ircld%ILOOP(J), :, J) = ircld%XSTR(:, J)
      INNER : DO ISTR = 1, ircld%ICOUNT1(J) - 1
        ircld%ILOOPIN(ircld%ILOOP(J), J) = ircld%ILOOPIN(ircld%ILOOP(J), J) + 1
        IF (ircld%XSTR(ISTR, J) == ircld%XSTR(ISTR + 1, J)) THEN
          IF (ircld%XSTR(ISTR, J) /= 1._jprb) THEN
            ircld%ICOUNT1(J) = ircld%ICOUNT1(J) - 1
            DO IJSTR = ISTR, ircld%ICOUNT1(J)
              ircld%XSTR(IJSTR, J) = ircld%XSTR(IJSTR + 1, J)
            ENDDO
            CYCLE OUTER
          ELSE
            EXIT OUTER
          ENDIF
        ENDIF
      ENDDO INNER
      EXIT OUTER
    ENDDO OUTER
    ircld%NSTREAM(J) = ircld%ICOUNT1(J) - 1
!---------Comppute the weight of the clear stream------------------------------
    DO ISTR = 1, ircld%NSTREAM(J)
      DO ILAY = 1, profiles(1)%nlayers
        IF (ircld%XSTRMIN(ILAY, J) <= ircld%XSTR(ISTR, J) .AND. ircld%XSTRMAX(ILAY, J) >= ircld%XSTR(ISTR + 1, J)) THEN
          ircld%ICLDARR(ISTR, ILAY, J) = 1_jpim
        ENDIF
      ENDDO
    ENDDO
    IF (ircld%NSTREAM(J) ==  - 1_jpim) THEN
      ircld%NSTREAM(J) = 0_jpim
    ENDIF
    ircld%XSTRCLR(J)    = 1._jprb - (ircld%XSTR(ircld%NSTREAM(J) + 1, J) - ircld%XSTR(1, J))
    ircld%NSTREAMREF(J) = ircld%NSTREAM(J)
!---------Consider only the streams whose weight is greater than cldstr_threshold-------
    IF (ircld%NSTREAM(J) /= 0_jpim) THEN
      DO ISTR = 1, ircld%NSTREAM(J)
        IF ((ircld%XSTR(ISTR + 1, J) - ircld%XSTR(ISTR, J)) >= cldstr_threshold) THEN
          ircld%ICOUNSTR(J)                    = ircld%ICOUNSTR(J) + 1
          ircld%INDEXSTR(ircld%ICOUNSTR(J), J) = ISTR
        ENDIF
      ENDDO
      IF (ircld%ICOUNSTR(J) /= 0_jpim) THEN
        DO ISTR = 1, ircld%ICOUNSTR(J)
          DO ILAY = 1, profiles(1)%nlayers
            ircld%ICLDARR(ISTR, ILAY, J) = ircld%ICLDARR(ircld%INDEXSTR(ISTR, J), ILAY, J)
          ENDDO
        ENDDO
        DO ISTR = 1, ircld%ICOUNSTR(J) + 1
          IF (ISTR == 1_jpim) THEN
            ircld%XSTR(ISTR, J) = ircld%XSTR(ircld%INDEXSTR(ISTR, J), J)
          ELSE
            ircld%XSTR(ISTR, J) = ircld%XSTR(ISTR - 1, J) +      &
              & (ircld%XSTR(ircld%INDEXSTR(ISTR - 1, J) + 1, J) - ircld%XSTR(ircld%INDEXSTR(ISTR - 1, J), J))
          ENDIF
        ENDDO
        ircld%XSTRCLR(J) = 1._jprb - (ircld%XSTR(ircld%ICOUNSTR(J) + 1, J) - ircld%XSTR(1, J))
        ircld%NSTREAM(J) = ircld%ICOUNSTR(J)
      ELSE IF (ircld%ICOUNSTR(J) == 0_jpim) THEN
        ircld%XSTRCLR(J) = 1._jprb
        ircld%NSTREAM(J) = 0_jpim
      ENDIF
    ELSE
      ircld%NSTREAM(J) = 0_jpim
    ENDIF
    nstreams = max(nstreams, ircld%NSTREAM(J))
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLDSTR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_cldstr
