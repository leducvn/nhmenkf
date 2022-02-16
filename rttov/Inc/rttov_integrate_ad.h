Interface
Subroutine rttov_integrate_ad( &
   addcosmic, opts, maxnstreams, chanprof,                       &
   emissivity,                   emissivity_ad,                  &
   reflectivity,                 reflectivity_ad,                &
   fresnrefl,                    fresnrefl_ad,                   &
   sunglint,                     sunglint_ad,                    &
   sun,                                                          &
   transmission_aux,             transmission_aux_ad,            &
   transmission_scatt_ir_stream, transmission_scatt_ir_stream_ad,&
   profiles,                     profiles_ad,                    &
   aux_prof,                     aux_prof_ad,                    &
   coef,                                                         &
   raytracing,                   raytracing_ad,                  &
   ircld,                        ircld_ad,                       &
   rad,                                                          &
   auxrad,                                                       &
                                 auxrad_stream,                  &
                                 auxrad_stream_ad,               &
                                 rad_ad)
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
Type(sunglint_type), Intent(inout)              :: sunglint_ad
Type(radiance_aux), Intent(inout)               :: auxrad_stream_ad
Real(jprb), Intent(inout)                       :: emissivity_ad(size(chanprof))
Real(jprb), Intent(inout)                       :: reflectivity_ad(size(chanprof))
Real(jprb), Intent(inout)                       :: fresnrefl_ad(size(chanprof))
Type(profile_Type), Intent(inout)               :: profiles_ad(size(profiles))
Type(profile_aux), Intent(inout)                :: aux_prof_ad
Type(transmission_Type_aux), Intent(inout)      :: transmission_aux_ad
type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_stream_ad
Type(ircld_type), intent(inout)                 :: ircld_ad
Type(raytracing_type), intent(inout)            :: raytracing_ad
Type(radiance_Type), Intent(inout)              :: rad_ad
End Subroutine
End Interface
