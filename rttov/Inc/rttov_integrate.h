Interface
subroutine rttov_integrate(addcosmic, opts, maxnstreams , chanprof, &! in
                           emissivity, reflectivity, fresnrefl, sunglint, &! in
                           sun, transmission_aux, transmission_scatt_ir_stream, &!in
                           profiles, aux_prof, coef, raytracing, ircld, &! in
                           rad, auxrad, auxrad_stream)      ! inout
  use parkind1, Only : jpim, jprb, jplm
  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, transmission_type_aux, &
                          transmission_scatt_ir_type, sunglint_type, radiance_type, ircld_type, raytracing_type, radiance_aux
  Implicit None
  logical(jplm),   intent(in)                   :: addcosmic   ! switch for adding cosmic background
  type(rttov_options),  intent(in)              :: opts        ! options structure
  integer(jpim),   intent(in)                   :: maxnstreams                !
  type(rttov_chanprof), intent(in)              :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)              :: profiles(:) ! Profiles
  real(jprb),      intent(in)                   :: emissivity(size(chanprof))   ! surface emissivity
  real(jprb),      intent(in)                   :: reflectivity(size(chanprof)) ! surface reflectivity
  real(jprb),      intent(in)                   :: fresnrefl(size(chanprof))
  logical(jplm),   intent(in)                   :: sun(size(chanprof))
  type(ircld_type),     intent(in)              :: ircld
  type(raytracing_type),intent(in)              :: raytracing
  type(sunglint_type),  intent(in)              :: sunglint
  type(transmission_type_aux), intent(inout)    :: transmission_aux            ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir_stream
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux),   intent(inout)           :: auxrad_stream
  type(radiance_type),  intent(inout)           :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance_aux),   intent(inout)           :: auxrad ! auxillary radiances
End Subroutine
End Interface
