!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_house_mod
!
! Purpose
! =======
!
! eigen_house_mod manages the modules for householder transformation.
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

      module eigen_house_mod

      implicit none

      real(8), allocatable   ::  u0_z(:), v0_z(:)
      real(8), allocatable   ::  u1_z(:), v1_z(:)
      real(8), allocatable   ::  tsave(:, :)
      integer                ::  offset1, offset2
      integer                ::  offset3, offset4

      integer, parameter     ::  LOOP_INFO_X_START = 1
      integer, parameter     ::  LOOP_INFO_X_END   = 2
      integer, parameter     ::  LOOP_INFO_X_OINDX = 3
      integer, parameter     ::  LOOP_INFO_X_ROOT  = 4
      integer, parameter     ::  LOOP_INFO_X_POS   = 5

      integer, parameter     ::  LOOP_INFO_Y_START = 6
      integer, parameter     ::  LOOP_INFO_Y_END   = 7
      integer, parameter     ::  LOOP_INFO_Y_OINDX = 8
      integer, parameter     ::  LOOP_INFO_Y_ROOT  = 9
      integer, parameter     ::  LOOP_INFO_Y_POS   = 10

      integer, parameter     ::  LOOP_INFO_V_WIDTH = 11
      integer, parameter     ::  LOOP_INFO_V_CURNT = 12
      integer, parameter     ::  LOOP_INFO_I_BASE  = 13

      integer, parameter     ::  LOOP_INFO_L_RANK  = 15
      integer, parameter     ::  LOOP_INFO_L_SIZE  = 16

      end module eigen_house_mod
