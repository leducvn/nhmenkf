#ifndef _THROW_H
#define _THROW_H

#define TRY INTEGER::L_; CHARACTER(LEN=256)::M_,F_,N_;N_="";M_=N_;F_=__FILE__;ERR=0

#define THROW_L(T,LAB_) IF(T)THEN;L_=__LINE__;GOTO LAB_;ENDIF
#define THROWM_L(T,M,LAB_) IF(T)THEN;L_=__LINE__;M_=M;GOTO LAB_;ENDIF

#define THROW(T) IF(T)THEN;L_=__LINE__;GOTO 999;ENDIF
#define THROWM(T,M) IF(T)THEN;L_=__LINE__;M_=M;GOTO 999;ENDIF

#define CATCH_L(LAB_) LAB_ WRITE(F_,"(A,':',I4.4)")TRIM(F_),L_;CALL RTTOV_ERRORREPORT(ef_,M_,F_)

#define CATCH RETURN;999 WRITE(F_,"(A,':',I4.4)")TRIM(F_),L_;CALL RTTOV_ERRORREPORT(ef_,M_,F_)
#define PCATCH STOP;999 WRITE(F_,"(A,':',I4.4)")TRIM(F_),L_;CALL RTTOV_ERRORREPORT(ef_,M_,F_);CALL RTTOV_EXIT(1_JPIM)

#define INFO(M) WRITE(N_,"(A,':',I4.4)")TRIM(F_),__LINE__;CALL RTTOV_ERRORREPORT(ei_,M,N_) 
#define WARN(M) WRITE(N_,"(A,':',I4.4)")TRIM(F_),__LINE__;CALL RTTOV_ERRORREPORT(ew_,M,N_)

      USE RTTOV_CONST, ONLY : ERRORSTATUS_FATAL
      USE RTTOV_CONST, ONLY : ERRORSTATUS_INFO 
      USE RTTOV_CONST, ONLY : ERRORSTATUS_WARNING
      USE RTTOV_CONST, ONLY : ERRORSTATUS_SUCCESS

      USE RTTOV_CONST, ONLY : ef_ => ERRORSTATUS_FATAL
      USE RTTOV_CONST, ONLY : ei_ => ERRORSTATUS_INFO 
      USE RTTOV_CONST, ONLY : ew_ => ERRORSTATUS_WARNING
      USE RTTOV_CONST, ONLY : es_ => ERRORSTATUS_SUCCESS

#endif