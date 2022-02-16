function rttov_opts_eq (opts1, opts2)
! Description:
!   Check equality of options structures
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  Use parkind1, Only : jplm
  Use rttov_types, Only : rttov_options
  Implicit None

  Logical(Kind=jplm) :: rttov_opts_eq
  
  Type(rttov_options), Intent(in) :: opts1
  Type(rttov_options), Intent(in) :: opts2

!INTF_END

! Only options meaningful for calculations
  rttov_opts_eq = &
   ( opts1%addinterp        .eqv. opts2%addinterp        ) .and. &
   ( opts1%addpc            .eqv. opts2%addpc            ) .and. &
   ( opts1%addradrec        .eqv. opts2%addradrec        ) .and. &
   ( opts1%addsolar         .eqv. opts2%addsolar         ) .and. &
   ( opts1%addaerosl        .eqv. opts2%addaerosl        ) .and. &
   ( opts1%addclouds        .eqv. opts2%addclouds        ) .and. &
   ( opts1%switchrad        .eqv. opts2%switchrad        ) .and. &
   ( opts1%spacetop         .eqv. opts2%spacetop         ) .and. &
   ( opts1%lgradp           .eqv. opts2%lgradp           ) .and. &
   ( opts1%use_q2m          .eqv. opts2%use_q2m          ) .and. &
   ( opts1%apply_reg_limits .eqv. opts2%apply_reg_limits ) .and. &
   ( opts1%ozone_data       .eqv. opts2%ozone_data       ) .and. &
   ( opts1%co2_data         .eqv. opts2%co2_data         ) .and. &
   ( opts1%n2o_data         .eqv. opts2%n2o_data         ) .and. &
   ( opts1%co_data          .eqv. opts2%co_data          ) .and. &
   ( opts1%ch4_data         .eqv. opts2%ch4_data         ) .and. &
   ( opts1%clw_data         .eqv. opts2%clw_data         ) .and. &
   ( opts1%addrefrac        .eqv. opts2%addrefrac        ) .and. &
   ( opts1%do_checkinput    .eqv. opts2%do_checkinput    )


End function
