Interface
Subroutine rttov_integrate_tl( &
   addcosmic, opts, maxnstreams, chanprof,      &! in
   emissivity, emissivity_tl,    &! in
   reflectivity, reflectivity_tl,  &! in
   fresnrefl, fresnrefl_tl,     &! in
   sunglint, sunglint_tl,      &! in
   sun,              &! in
   transmission_aux, transmission_aux_tl,  &! in
   transmission_scatt_ir_stream, transmission_scatt_ir_stream_tl,&
   profiles, profiles_tl,      &! in
   aux_prof, aux_prof_tl,      &! in
   coef,     &! in
   raytracing, raytracing_tl,    &! in
   ircld, ircld_tl,         &! in
   rad , &! in
   auxrad ,          &! in
   auxrad_stream, auxrad_stream_tl, &! in
   rad_tl           ) ! inout
  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, transmission_type_aux, & 
       transmission_scatt_ir_type, sunglint_type, radiance_type, ircld_type, raytracing_type, radiance_aux
  use parkind1, Only : jpim, jprb, jplm
  use rttov_const, Only : sensor_id_po
  use yomhook, Only : LHOOK, DR_HOOK
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
  type(transmission_type_aux), intent(in)       :: transmission_aux            ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir_stream
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux),   intent(in)              :: auxrad_stream
  type(radiance_type),  intent(in)              :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance_aux),   intent(in)              :: auxrad ! auxillary radiances
  Real(jprb),                Intent(in)         :: emissivity_tl(size(chanprof))
  Real(jprb),                Intent(in)         :: reflectivity_tl(size(chanprof))
  Real(jprb),                Intent(in)         :: fresnrefl_tl(size(chanprof))
  Type(profile_Type),  Intent(in)               :: profiles_tl(size(profiles))
  Type(ircld_type)                 ,intent(in)  :: ircld_tl
  Type(raytracing_type),intent(in)              :: raytracing_tl
  Type(sunglint_Type), Intent(in)               :: sunglint_tl
  Type(transmission_Type_aux), Intent(in)       :: transmission_aux_tl
  type(transmission_scatt_ir_type) ,intent(in)  :: transmission_scatt_ir_stream_tl
  Type(profile_aux) ,  Intent(in)               :: aux_prof_tl
  Type(radiance_Type), Intent(inout)            :: rad_tl ! in because of mem allocation
  Type(radiance_aux),  Intent(inout)            :: auxrad_stream_tl
End Subroutine
End Interface
