!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: CSTAB_mod
!
! Purpose
! =======
!
! CSTAB_mod mnages the modules for Cache stabilization routines.
!
!
! Copyright(C) 2012-2018 RIKEN.
! Copyright(C) 2011-2012 Toshiyuki Imamura
!                        Graduate School of Informatics and Engineering,
!                        The University of Electro-Communications.
! Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
! 
! Redistribution  and  use  in  source and binary forms, with or without
! modification,  are  permitted  provided  that the following conditions
! are met:
! 
! * Redistributions  of  source  code  must  retain  the above copyright
!   notice,  this  list  of  conditions  and  the  following  disclaimer.
! * Redistributions  in  binary  form must reproduce the above copyright
!   notice,  this list of conditions and the following disclaimer in the
!   documentation  and/or other materials provided with the distribution.
! 
! THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
! LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module CSTAB_mod

      implicit none
      private


      public :: CSTAB_get_optdim
      public :: CSTAB_adjust_base
      public :: CSTAB_adjust_page
      public :: CSTAB_round_offset

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_get_optdim
!
! Purpose
! =======
!
! Return the recommendational dimension of the array to minimize the
! risk of Cache thrashing
!
!
! Arguments
! =========
!
! n_min   (input) integer
!         minimum dimension of the array
!
! n_unroll (input) integer
!          possible unroll factor
!
! delta_L1 (input) integer
!          window size of L1
!
! delta_L2 (input) integer
!          window size of L2
!
! n_opt   (output) integer
!         recommendational dimension of the array
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_get_optdim(n_min, n_unroll,
     &     delta_L1, delta_L2, n_opt)

      integer, intent(in)    :: n_min
      integer, intent(in)    :: n_unroll
      integer, intent(in)    :: delta_L1
      integer, intent(in)    :: delta_L2
      integer, intent(out)   :: n_opt

      include 'CSTAB.h'

      real(8), parameter     :: ONE = 1d0

      integer                :: n_delta
      integer                :: n_opt2
      integer                :: n_delta2
      integer                :: i,j,k


      n_opt=n_min

      do

        n_opt   = (n_opt-1)/L1_WINDOW+1
        n_opt   = (n_opt/2)*2+1
        n_opt   = n_opt*L1_WINDOW

        n_delta = 0

        do i=1,int((n_unroll*1.2d0-ONE)/L1_WAY+1)

          k=mod(i*n_opt+L1_LSIZE/2,L1_LSIZE)-L1_LSIZE/2
          if (abs(k)<=delta_L1/2) then
            n_delta=(delta_L1/2-k-1)/i+1
            goto 10000
          end if

        end do

        do i=1,int((n_unroll*1.2d0-ONE)/L2_WAY+1)

          k=mod(i*n_opt+L2_LSIZE/2,L2_LSIZE)-L2_LSIZE/2
          if (abs(k)<=delta_L2/2) then
            n_delta=(delta_L2/2-k-1)/i+1
            goto 10000
          end if

        end do

10000   continue

        if (n_delta==0) exit
        n_opt = n_opt + n_delta

      end do

      return

      end subroutine CSTAB_get_optdim

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_adjust_base
!
! Purpose
! =======
!
! Adjust the cache line position
!
!
! Arguments
! =========
!
! a       (input/output) real array, dimension(*)
!         top position of the array a
!
! b       (input/output) real array, dimension(*)
!         top position of the array b
!
! offset  (input/output) integer
!         offset to adjust arrays a and b
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_adjust_base(a,b,offset)

      real(8)                :: a(*)
      real(8)                :: b(*)
      integer                :: offset

      include 'CSTAB.h'


      call get_delta(a(1),b(1),offset)
      offset=(offset/8)

      if (offset>0) then
        offset=mod(L2_LSIZE-mod(+offset,L2_LSIZE),L2_LSIZE)
      else
        offset=mod(L2_LSIZE+mod(-offset,L2_LSIZE),L2_LSIZE)
      end if

      return

      end subroutine CSTAB_adjust_base

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_adjust_page
!
! Purpose
! =======
!
! Adjust the page position
!
!
! Arguments
! =========
!
! a       (input/output) real array, dimension(*)
!         top position of the array a
!
! b       (input/output) real array, dimension(*)
!         top position of the array b
!
! offset  (input/output) integer
!         offset to adjust arrays a and b
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_adjust_page(a,b,offset)

      real(8)                :: a(*)
      real(8)                :: b(*)
      integer                :: offset

      include 'CSTAB.h'


      call get_delta(a(1),b(1),offset)
      offset=(offset/8)
      if (offset>0) then
        offset = mod(PAGE_LSIZE-mod(+offset,PAGE_LSIZE),PAGE_LSIZE)
      else
        offset = mod(PAGE_LSIZE+mod(-offset,PAGE_LSIZE),PAGE_LSIZE)
      end if

      return

      end subroutine CSTAB_adjust_page

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_round_offset
!
! Purpose
! =======
!
! Round up the offset by pagesize
!
!
! Arguments
! =========
!
! offset  (input/output) integer
!         offset
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_round_offset(offset)

      integer, intent(inout) :: offset

      include 'CSTAB.h'


      offset = mod(offset,L2_LSIZE)

      return

      end subroutine CSTAB_round_offset

      end module CSTAB_mod
