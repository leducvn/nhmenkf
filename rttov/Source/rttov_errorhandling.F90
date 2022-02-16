!
Subroutine rttov_errorhandling (err_unit, verbosity_level)
  ! Description:
  ! Handling of error messages.
  ! Error messages will be sent on the optional unit number errunit.
  !      Default is the value defined in the module for constants.
  !
  ! The levels of verbosity are 
  !  0 = no error messages output 
  !  1 = FATAL errors only printed. these are errors which 
  !      mean that profile should be aborted (e.g. unphysical 
  !      profile input) 
  !  2 = WARNING errors only printed. Errors which can allow 
  !      the computation to continue but the results may be 
  !      suspect (e.g. profile outside basis profile limits) 
  !  3 = INFORMATION messages which inform the user about 
  !      the computation
  !
  ! Copyright:
  !
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  !
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version    Date       Comment
  !  1.0    10/03/2003   Original code (P Brunel)
  !  1.1    07/12/2007   add argument for control of checkinput warnings P.B.
  !
  ! Code Description:
  !   FORTRAN 90, following AAPP standards
  !
  ! Declarations
  !
  ! Global variables:
  ! Modules used:
  !

  Use parkind1, Only : jpim

!INTF_OFF
Use rttov_const, Only : & 
    errorstatus_info,& 
    errorstatus_success,& 
    errorstatus_fatal,& 
    errorstatus_warning,& 
    default_err_unit 
Use rttov_global, Only : & 
    verbose_message,& 
    err_init,& 
    error_unit
!INTF_ON

  Implicit None
  !
  ! Subroutine arguments
  !   Scalar arguments with intent(in):
  Integer(Kind=jpim), Intent (in) :: err_unit       ! Logical error unit
  Integer(Kind=jpim), Intent (in) :: verbosity_level

!INTF_END

  !- End of header --------------------------------------------------------

  ! According to the user defined verbosity level 
  ! defines the verbose flag for each error level 
  ! Default is to output all messages 
  Select Case ( verbosity_level ) 
  Case ( 0_jpim ) 
     verbose_message( errorstatus_success ) = .false. 
     verbose_message( errorstatus_warning ) = .false. 
     verbose_message( errorstatus_fatal   ) = .false. 
     verbose_message( errorstatus_info    ) = .false. 
  Case ( 1_jpim ) 
     verbose_message( errorstatus_success ) = .false. 
     verbose_message( errorstatus_warning ) = .false. 
     verbose_message( errorstatus_fatal   ) = .true. 
     verbose_message( errorstatus_info    ) = .false. 
  Case ( 2_jpim ) 
     verbose_message( errorstatus_success ) = .false. 
     verbose_message( errorstatus_warning ) = .true. 
     verbose_message( errorstatus_fatal   ) = .true. 
     verbose_message( errorstatus_info    ) = .false. 
  Case ( 3_jpim ) 
     verbose_message( errorstatus_success ) = .true. 
     verbose_message( errorstatus_warning ) = .true. 
     verbose_message( errorstatus_fatal   ) = .true. 
     verbose_message( errorstatus_info    ) = .true. 
  Case Default 
     verbose_message( errorstatus_success ) = .true. 
     verbose_message( errorstatus_warning ) = .true. 
     verbose_message( errorstatus_fatal   ) = .true. 
     verbose_message( errorstatus_info    ) = .true. 
  End Select

  ! Definition of the error message logical unit
  ! default is taken from the module for constants
  If( err_unit >= 0 ) then
     error_unit = err_unit
  Else
     error_unit = default_err_unit
  EndIf

  ! Setup initialisation flag
  ! This flag is tested by the rttov_errorreport subroutine
  err_init = .true.

End Subroutine rttov_errorhandling
