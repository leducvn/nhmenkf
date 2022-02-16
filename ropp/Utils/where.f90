! $Id: where.f90 1882 2008-10-27 15:45:52Z frhl $

!****f* Arrays/where *
!
! NAME
!    where - Return the subscripts for which a condition is fullfilled.
!
! SYNOPSIS
!    idx => where(mask [, n])
! 
! DESCRIPTION
!    This function returns the subscripts of the true elements of mask.
!
! INPUTS
!    logical      :: mask(:)   The array to be scanned.
!
! OUTPUT
!    int, pointer :: idx(:)    An array of indices where mask == .true.
!
! REFERENCES
!    This is a Fortran 90 implementation of IDL's where() function.
!
! AUTHOR
!    C. Marquardt, West Hill, UK     <christian@marquardt.fsnet.co.uk>
!
!**** 

function where(mask, n) result (indices)

  logical, dimension(:), intent(in)       :: mask
  integer, dimension(:), pointer          :: indices

  integer,               optional         :: n
  integer                                 :: nc, i, j

  nc = count(mask)
  if (present(n)) then
     n = nc
  endif

  nullify(indices)
  allocate(indices(nc))

  j = 1
  do i = 1, size(mask)
     if (mask(i)) then
        indices(j) = i
        j = j + 1
     end if
  end do

end function where
