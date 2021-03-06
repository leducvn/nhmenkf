Interface
Subroutine rttov_INTEX(    &
       & klevi ,  &! in
       & klevf ,  &! in
       & presi ,  &! in
       & presf ,  &! in
       & veci  ,  &! in
       & vecf  )   ! out
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer(Kind=jpim), Intent(in) :: klevi      ! number of levels of the initial grid
  Integer(Kind=jpim), Intent(in) :: klevf      ! number of levels of the final grid
  Real(Kind=jprb), Intent(in), Dimension(klevi)  :: presi ! initial grid
  Real(Kind=jprb), Intent(in), Dimension(klevf)  :: presf ! final grid
  Real(Kind=jprb), Intent(in), Dimension(klevi)  :: veci  ! initial vec array
  Real(Kind=jprb), Intent(out), Dimension(klevf) :: vecf  ! final vec array
End Subroutine
End Interface
