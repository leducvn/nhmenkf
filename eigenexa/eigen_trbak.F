!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trbak_mod
!
! Purpose
! =======
!
! eigen_trbak_mod manages the modules of the back-transformation.
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

      module eigen_trbak_mod

!$    use omp_lib

      implicit none
      public

      real(8) :: trbk_time_bcast
      real(8) :: trbk_time_reduc
      real(8) :: trbk_time_fr
      real(8) :: trbk_time_trbk1
      real(8) :: trbk_time_trbk1_
      real(8) :: trbk_time_trbk1x
      real(8) :: trbk_time_trbk1x_
      real(8) :: trbk_time_trbk2

      real(8) :: trbk_time_reduc_overhead_x

!$    integer(kind=omp_lock_kind) :: trbk_lock
!$    integer(8)                  :: trbk_mask(1:2)

      integer :: do_overlap_bcast_level
      integer :: trbk_time_counter
      integer :: trbk_time_interval
      integer :: trbk_time_next
      integer :: trbk_switched

      real(8), pointer :: trbk_buf(:)

      end module eigen_trbak_mod
