Interface
Subroutine rttov_integrate_k( &
   addcosmic, opts, maxnstreams, chanprof,                       &
   emissivity,                   emissivity_k,                  &
   reflectivity,                 reflectivity_k,                &
   fresnrefl,                    fresnrefl_k,                   &
   sunglint,                     sunglint_k,                    &
   sun,                                                          &
   transmission_aux,             transmission_aux_k,            &
   transmission_scatt_ir_stream, transmission_scatt_ir_stream_k,&
   profiles,                     profiles_k,                    &
   aux_prof,                     aux_prof_k,                    &
   coef,                                                         &
   raytracing,                   raytracing_k,                  &
   ircld,                        ircld_k,                       &
   rad,                                                          &
   auxrad,                                                       &
                                 auxrad_stream,                  &
                                 auxrad_stream_k,               &
                                 rad_k)
Use rttov_types, Only : &
   rttov_chanprof, rttov_coef, profile_Type, profile_aux, transmission_type_aux, transmission_scatt_ir_type, sunglint_type, &
   radiance_Type, rttov_options, ircld_type, raytracing_type, radiance_aux
Use parkind1, Only : jpim, jprb, jplm
Implicit None
Logical(jplm), Intent(in)                       :: addcosmic
Type(rttov_options), Intent(in)                 :: opts
Integer(jpim), Intent(in)                       :: maxnstreams
Type(rttov_chanprof), Intent(in)                :: chanprof(:)
Type(profile_Type), Intent(in)                  :: profiles(:)
Real(jprb), Intent(in)                          :: emissivity(size(chanprof))
Real(jprb), Intent(in)                          :: reflectivity(size(chanprof))
Real(jprb), Intent(in)                          :: fresnrefl(size(chanprof))
Logical(jplm), Intent(in)                       :: sun(size(chanprof))
Type(rttov_coef), Intent(in)                    :: coef
Type(profile_aux) , Intent(in)                  :: aux_prof
Type(transmission_Type_aux), Intent(in)         :: transmission_aux
type(transmission_scatt_ir_type), intent(in)    :: transmission_scatt_ir_stream
Type(ircld_type), intent(in)                    :: ircld
Type(raytracing_type), intent(in)               :: raytracing
Type(sunglint_type), Intent(in)                 :: sunglint
Type(radiance_Type), Intent(in)                 :: rad
Type(radiance_aux), Intent(in)                  :: auxrad
Type(radiance_aux), Intent(in)                  :: auxrad_stream
Type(sunglint_type), Intent(inout)              :: sunglint_k
Type(radiance_aux), Intent(inout)               :: auxrad_stream_k
Real(jprb), Intent(inout)                       :: emissivity_k(size(chanprof))
Real(jprb), Intent(inout)                       :: reflectivity_k(size(chanprof))
Real(jprb), Intent(inout)                       :: fresnrefl_k(size(chanprof))
Type(profile_Type), Intent(inout)               :: profiles_k(size(chanprof))
Type(profile_aux), Intent(inout)                :: aux_prof_k
Type(transmission_Type_aux), Intent(inout)      :: transmission_aux_k
type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_stream_k
Type(ircld_type), intent(inout)                 :: ircld_k
Type(raytracing_type), intent(inout)            :: raytracing_k
Type(radiance_Type), Intent(inout)              :: rad_k
End Subroutine
End Interface
