Interface
function rttov_opts_eq (opts1, opts2)
  Use parkind1, Only : jplm
  Use rttov_types, Only : rttov_options
  Implicit None
  Logical(Kind=jplm) :: rttov_opts_eq
  Type(rttov_options), Intent(in) :: opts1
  Type(rttov_options), Intent(in) :: opts2
End Function
End Interface
