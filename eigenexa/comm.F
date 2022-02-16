!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: comm_mod
!
! Purpose
! =======
!
! comm_mod manages the modules for MPI-communication wrapper routines
! and controls internal algorithms.
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

      module comm_mod

      use eigen_devel_mod
      use mpi

      implicit none
      private


      public  :: pack_dbl
      public  :: pack1_dbl
      public  :: unpack_dbl
      public  :: unpack1_dbl
      public  :: send_dbl
      public  :: send_dblt
      public  :: isend_dbl
      public  :: isend_dblt
      public  :: recv_dbl
      public  :: irecv_dbl
      public  :: irecv_dblt
      public  :: wait_dbl
      public  :: waitall_dbl
      public  :: barrier
      public  :: bcast_dbl
      public  :: bcastw_dbl
      public  :: reduce_dbl
      public  :: allgather_dbl
      public  :: datacast_dbl
      public  :: datacast_dbl2
      public  :: datacast_dblx
      public  :: print_bcast_algorithm
      public  :: print_reduce_algorithm
      public  :: print_gather_algorithm

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: pack_dbl
!
! Purpose
! =======
!
! Packs the buffer
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         original buffer
!
! n       (input) integer
!         number of elements in original buffer
!
! buffer  (output) real(8) array, dimension(ptr+n)
!         packed buffer
!
! ptr     (input/output) integer
!         pointer to the buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

!DIR$ ATTRIBUTES FORCEINLINE :: pack_dbl
      subroutine  pack_dbl(buf, n, buffer, ptr)
      
      real(8), intent(in)    :: buf(*)
      integer, intent(in)    :: n
      real(8), intent(out)   :: buffer(*)
      integer, intent(inout) :: ptr


      buffer(ptr+1:ptr+n) = buf(1:n)
      ptr = ptr + n

      return

      end subroutine pack_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: pack1_dbl
!
! Purpose
! =======
!
! Packs the buffer
!
!
! Arguments
! =========
!
! buf     (input/output) real(8)
!         original buffer
!
! buffer  (output) real(8) array, dimension(ptr+1)
!         packed buffer
!
! ptr     (input/output) integer
!         pointer to the buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

!DIR$ ATTRIBUTES FORCEINLINE :: pack1_dbl
      subroutine  pack1_dbl(buf, buffer, ptr)
      
      real(8), intent(in)    :: buf
      real(8), intent(out)   :: buffer(*)
      integer, intent(inout) :: ptr


      buffer(ptr+1:ptr+1) = buf
      ptr = ptr + 1

      return

      end subroutine pack1_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: unpack_dbl
!
! Purpose
! =======
!
! Un-Packs the buffer
!
!
! Arguments
! =========
!
! buf     (output) real(8) array, dimension(n)
!         original buffer
!
! n       (input) integer
!         number of elements in original buffer
!
! buffer  (input) real(8) array, dimension(ptr+n)
!         packed buffer
!
! ptr     (input/output) integer
!         pointer to the buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

!DIR$ ATTRIBUTES FORCEINLINE :: unpack_dbl
      subroutine unpack_dbl(buf, n, buffer, ptr)
      
      real(8), intent(out)   :: buf(*)
      integer, intent(in)    :: n
      real(8), intent(in)    :: buffer(*)
      integer, intent(inout) :: ptr


      buf(1:n) = buffer(ptr+1:ptr+n)
      ptr = ptr + n

      return

      end subroutine unpack_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: unpack1_dbl
!
! Purpose
! =======
!
! Un-Packs the buffer
!
!
! Arguments
! =========
!
! buf     (output) real(8)
!         original buffer
!
! buffer  (input) real(8) array, dimension(ptr+n)
!         packed buffer
!
! ptr     (input/output) integer
!         pointer to the buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

!DIR$ ATTRIBUTES FORCEINLINE :: unpack1_dbl
      subroutine unpack1_dbl(buf, buffer, ptr)
      
      real(8), intent(out)   :: buf
      real(8), intent(in)    :: buffer(*)
      integer, intent(inout) :: ptr


      buf = buffer(ptr+1)
      ptr = ptr + 1

      return

      end subroutine unpack1_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: send_dbl
!
! Purpose
! =======
!
! Performs a blocking send
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send buffer
!
! n       (input) integer
!         number of elements in send buffer
!
! idest   (input) integer
!         rank of destination
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine send_dbl(buf, n, idest, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: idest
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Send(buf, n, MPI_DOUBLE_PRECISION,
     &     idest-1, 1, icom, ierr)

      return

      end subroutine send_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: send_dblt
!
! Purpose
! =======
!
! Performs a blocking send with message tag
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send buffer
!
! n       (input) integer
!         number of elements in send buffer
!
! idest   (input) integer
!         rank of destination
!
! itag    (input) integer
!         message tag
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine send_dblt(buf, n, idest, itag, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: idest
      integer, intent(in)    :: itag
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Send(buf, n, MPI_DOUBLE_PRECISION,
     &     idest-1, itag, icom, ierr)

      return

      end subroutine send_dblt

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: isend_dbl
!
! Purpose
! =======
!
! Begins a nonblocking send
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send buffer
!
! n       (input) integer
!         number of elements in send buffer
!
! idest   (input) integer
!         rank of destination
!
! ireq    (input/output) integer
!         communication request
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine isend_dbl(buf, n, idest, ireq, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: idest
      integer, intent(inout) :: ireq
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Isend(buf, n, MPI_DOUBLE_PRECISION,
     &     idest-1, 1, icom, ireq, ierr)

      return

      end subroutine isend_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: isend_dblt
!
! Purpose
! =======
!
! Begins a nonblocking send with message tag
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send buffer
!
! n       (input) integer
!         number of elements in send buffer
!
! idest   (input) integer
!         rank of destination
!
! itag    (input) integer
!         message tag
!
! ireq    (input/output) integer
!         communication request
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine isend_dblt(buf, n, idest, itag, ireq, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: idest
      integer, intent(in)    :: itag
      integer, intent(inout) :: ireq
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Isend(buf, n, MPI_DOUBLE_PRECISION,
     &     idest-1, itag, icom, ireq, ierr)

      return

      end subroutine isend_dblt

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: recv_dbl
!
! Purpose
! =======
!
! Blocking receive for a message
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         receive buffer
!
! n       (input) integer
!         maximum number of elements in receive buffer
!
! isrc    (input) integer
!         rank of source
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine recv_dbl(buf, n, isrc, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: isrc
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION,
     &     isrc-1, 1, icom, MPI_STATUS_IGNORE, ierr)

      return

      end subroutine recv_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: recv_dblt
!
! Purpose
! =======
!
! Blocking receive for a message with message tag
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         receive buffer
!
! n       (input) integer
!         maximum number of elements in receive buffer
!
! isrc    (input) integer
!         rank of source
!
! itag    (input) integer
!         message tag
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine recv_dblt(buf, n, isrc, itag, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: isrc
      integer, intent(in)    :: itag
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION,
     &     isrc-1, itag, icom, MPI_STATUS_IGNORE, ierr)

      return

      end subroutine recv_dblt

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: irecv_dbl
!
! Purpose
! =======
!
! Begins a nonblocking receive
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         receive buffer
!
! n       (input) integer
!         maximum number of elements in receive buffer
!
! isrc    (input) integer
!         rank of source
!
! ireq    (input/output) integer
!         communication request
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine irecv_dbl(buf, n, isrc, ireq, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: isrc
      integer, intent(inout) :: ireq
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Irecv(buf, n, MPI_DOUBLE_PRECISION,
     &     isrc-1, 1, icom, ireq, ierr)

      return

      end subroutine irecv_dbl
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: irecv_dblt
!
! Purpose
! =======
!
! Begins a nonblocking receive with message tag
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         receive buffer
!
! n       (input) integer
!         maximum number of elements in receive buffer
!
! isrc    (input) integer
!         rank of source
!
! itag    (input) integer
!         message tag
!
! ireq    (input/output) integer
!         communication request
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine irecv_dblt(buf, n, isrc, itag, ireq, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: isrc
      integer, intent(in)    :: itag
      integer, intent(inout) :: ireq
      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Irecv(buf, n, MPI_DOUBLE_PRECISION,
     &     isrc-1, itag, icom, ireq, ierr)

      return

      end subroutine irecv_dblt

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: wait_dbl
!
! Purpose
! =======
!
! Waits for an MPI request to complete 
!
!
! Arguments
! =========
!
! ireq    (input/output) integer
!         request
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine wait_dbl(ireq)

      integer, intent(inout) :: ireq

      integer                :: ierr


      call MPI_Wait(ireq, MPI_STATUS_IGNORE, ierr)

      return

      end subroutine wait_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: waitall_dbl
!
! Purpose
! =======
!
! Waits for all given MPI Requests to complete
!
!
! Arguments
! =========
!
! n       (input) integer
!         list length
!
! ireq    (input/output) integer array, dimension(n)
!         requests
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine waitall_dbl(n, ireq)

      integer, intent(in)    :: n
      integer, intent(inout) :: ireq(n)

      integer                :: ierr


      call MPI_Waitall(n, ireq, MPI_STATUSES_IGNORE, ierr)

      return

      end subroutine waitall_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: barrier
!
! Purpose
! =======
!
! Sync all processes
!
!
! Arguments
! =========
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine barrier(icom)

      integer, intent(in)    :: icom

      integer                :: ierr


      call MPI_Barrier(icom, ierr)

      return

      end subroutine barrier

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: bcast_dbl
!
! Purpose
! =======
!
! Broadcast data
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send/receive buffer
!
! n       (input) integer
!         number of elements in buffer
!
! iroot   (input) integer
!         rank of broadcast root
!  
! col_id  (input) integer
!         color identity
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine bcast_dbl(buf, n, iroot, col_id, icom)

      real(8), intent(inout) :: buf(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: iroot
      integer, intent(in)    :: col_id
      integer, intent(in)    :: icom

      integer, parameter     :: BCAST_TAG = 100000

      integer                :: ierr, i, j, k, L, Glog, Ghalf, tag

#if TIMER_PRINT
      real(8)                :: timer
#endif

      integer                :: my_rank, world_size, local_rank
      integer                :: rank_from, rank_to
      integer, pointer       :: ireq(:)


#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif


#define BCAST_MPI                       (-1)
#define BCAST_MPI_SEGMENTED             (1)
#define BCAST_SEQUENTIAL                (2)
#define BCAST_BINOMIAL                  (3)
#define BCAST_TRINOMIAL                 (4)
#define BCAST_RELAY_ONEDIR              (5)
#define BCAST_RELAY_TWODIRS             (6)
#define BCAST_UNKNOWN                   (0)

#define BCAST_SEGMENT_SIZE              (4096)

#if defined(__FUJITSU)
!     #define BCAST_ALGORITHM         BCAST_MPI
!     #define BCAST_ALGORITHM         BCAST_BINOMIAL
#define BCAST_ALGORITHM         BCAST_TRINOMIAL
!     #define BCAST_ALGORITHM         BCAST_RELAY_ONEDIR
!     #define BCAST_ALGORITHM         BCAST_RELAY_TWODIRS
#else
#define BCAST_ALGORITHM         BCAST_MPI
#endif

#if BCAST_ALGORITHM == BCAST_UNKNOWN || !defined(BCAST_ALGORITHM)
      error
#endif


      if (n < 1) return

      call MPI_Comm_size(icom, world_size, ierr)
      if (world_size == 1) then
        return
      end if

#if defined(__FUJITSU)
!     On Fujitsu K and FX10 computers, special algorithm may be
!     applied  when communicator is COMM_WORLD.
!     In case n > 512, MPI performs best, thus, is forcibly selected.

      if (icom == TRD_COMM_WORLD 
     &     .or. n > max(512, world_size)) then
!        do i=1,n,BCAST_SEGMENT_SIZE
!          j = min(n-i+1, BCAST_SEGMENT_SIZE)
!          call MPI_Bcast(buf(i), j, MPI_DOUBLE_PRECISION,
!     &         iroot-1, icom, ierr)
!        end do
        call MPI_Bcast(buf, n, MPI_DOUBLE_PRECISION,
     &       iroot-1, icom, ierr)
        goto 99999
      end if
#endif

      call MPI_Comm_rank(icom, my_rank, ierr)
      my_rank = my_rank + 1

#if BCAST_ALGORITHM == BCAST_MPI
      call MPI_Bcast(buf, n, MPI_DOUBLE_PRECISION,
     &     iroot-1, icom, ierr)
#endif

#if BCAST_ALGORITHM == BCAST_MPI_SEGMENTED

!     On Fujitsu K and FX10 computers, the following message dividing
!     on broadcasting reduces the unexpecetd communication overhead
!     when the message size exeeds 4096, somehow.

      do i=1,n,BCAST_SEGMENT_SIZE
        j = min(n-i+1, BCAST_SEGMENT_SIZE)
        call MPI_Bcast(buf(i), j, MPI_DOUBLE_PRECISION,
     &       iroot-1, icom, ierr)
      end do
#endif

#if BCAST_ALGORITHM == BCAST_SEQUENTIAL
!     
!     Sequential algorithm
!     
      if (iroot == my_rank) then
        allocate(ireq(1:world_size))
        j = 0
        do i=1, world_size
          if (i /= iroot) then
            j = j + 1
            call isend_dblt(buf, n, i, i, ireq(j), icom)
          end if
        end do
        if (j > 0) then
          call waitall_dbl(j, ireq)
        end if
        deallocate(ireq)
      else
        call recv_dblt(buf, n, iroot, my_rank, icom)
      end if
#endif

#if BCAST_ALGORITHM == BCAST_BINOMIAL
!     
!     Binomial algorithm
!     
      local_rank = mod(world_size+my_rank-iroot, world_size) ! relative 
                                !rank
      i=1; Glog=0
      do
        if (i >= world_size) exit
        Glog = Glog + 1
        i = i * 2
      end do
      allocate(ireq(1:Glog))
      i=1; j=0
      do
        if (i >= world_size) exit
        if (local_rank < i) then
          if (local_rank+i < world_size) then
            k = mod(my_rank+i-1, world_size)+1
            j = j + 1
            tag = BCAST_TAG + k
            call isend_dblt(buf, n, k, tag, ireq(j), icom)
          end if
        else if (local_rank < 2*i) then
          k = mod(world_size+my_rank-i-1, world_size)+1
          tag = BCAST_TAG + my_rank
          call recv_dblt(buf, n, k, tag, icom)
        end if
        i = i * 2
      end do
      if (j > 0) then
        call waitall_dbl(j, ireq)
      end if
      deallocate(ireq)
#endif

#if BCAST_ALGORITHM == BCAST_TRINOMIAL
!     
!     Binomial algorithm
!     
      local_rank = mod(world_size+my_rank-iroot, world_size) ! relativ
                                !e rank
      i=1; Glog=0
      do
        if (i >= world_size) exit
        Glog = Glog + 2
        i = i * 3
      end do
      allocate(ireq(1:Glog))
      i=1; j=0
      do
        if (i >= world_size) exit
        if (local_rank < i) then
          if (local_rank+1*i < world_size) then
            k = mod(my_rank+1*i-1, world_size)+1
            j = j + 1
            tag = BCAST_TAG + k
            call isend_dblt(buf, n, k, tag, ireq(j), icom)
          end if
          if (local_rank+2*i < world_size) then
            k = mod(my_rank+2*i-1, world_size)+1
            j = j + 1
            tag = BCAST_TAG + k
            call isend_dblt(buf, n, k, tag, ireq(j), icom)
          end if
        else if (local_rank < 2*i) then
          k = mod(world_size+my_rank-1*i-1, world_size)+1
          tag = BCAST_TAG + my_rank
          call recv_dblt(buf, n, k, tag, icom)
        else if (local_rank < 3*i) then
          k = mod(world_size+my_rank-2*i-1, world_size)+1
          tag = BCAST_TAG + my_rank
          call recv_dblt(buf, n, k, tag, icom)
        end if
        i = i * 3
      end do
      if (j > 0) then
        call waitall_dbl(j, ireq)
      end if
      deallocate(ireq)
#endif

#if BCAST_ALGORITHM == BCAST_RELAY_ONEDIR
!     
!     Relay algorithm
!     
      tag = BCAST_TAG
      rank_to   = mod((my_rank-1)+1,            world_size) + 1
      rank_from = mod((my_rank-1)-1+world_size, world_size) + 1

      if (my_rank /= iroot) then
        call recv_dblt(buf, n, rank_from, tag, icom)
      end if
      if (rank_to /= iroot) then
        call send_dblt(buf, n, rank_to,   tag, icom)
      end if
#endif

#if BCAST_ALGORITHM == BCAST_RELAY_TWODIRS
!     
      tag = BCAST_TAG
      Ghalf = world_size/2
      local_rank = mod(world_size+my_rank-iroot, world_size) ! relative
                                ! rank
      if (local_rank >= Ghalf) then
        local_rank = local_rank - world_size
      end if

      if (my_rank == iroot) then
        allocate(ireq(1:2))
        i   = mod((my_rank-1)+1,            world_size)+1
        j   = mod((my_rank-1)-1+world_size, world_size)+1
        if (i /= j) then
          call isend_dblt(buf, n, i,   tag, ireq(1), icom)
          call isend_dblt(buf, n, j,   tag, ireq(2), icom)
          call waitall_dbl(2, ireq)
        else
          call isend_dblt(buf, n, i,   tag, ireq(1), icom)
          call wait_dbl(ireq(1))
        end if
        deallocate(ireq)
      else
        if (local_rank > 0) then
          rank_to   = mod((my_rank-1)+1,            world_size)+1
          rank_from = mod((my_rank-1)-1+world_size, world_size)+1
        else
          rank_to   = mod((my_rank-1)-1+world_size, world_size)+1
          rank_from = mod((my_rank-1)+1,            world_size)+1
        end if

        call recv_dblt(buf, n, rank_from, tag, icom)

        k = mod(world_size+rank_to-iroot, world_size) ! relative rank
        if (k >= Ghalf) k = k - world_size
        if (k * local_rank > 0) then
          call send_dblt(buf, n, rank_to,   tag, icom)
        end if
      end if
#endif


99999 continue

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      time_bcast = time_bcast + timer
      if (col_id > 0) then
        if (col_id >= 1 .and. col_id <= items_bcast) then
          time_bcast_(col_id)     = time_bcast_(col_id) + timer
          counter_bcast_(col_id)  = counter_bcast_(col_id) + n
          messages_bcast_(col_id) = messages_bcast_(col_id) + 1
        end if
      else
        if (-col_id >= 1 .and.
     &      -col_id <= items_bcast+items_redist) then
          time_bcast_(-col_id)     = time_bcast_(-col_id) + timer
          counter_bcast_(-col_id)  = counter_bcast_(-col_id) + n
          messages_bcast_(-col_id) = messages_bcast_(-col_id) + 1
        end if
      end if
#endif

      return

      end subroutine bcast_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: bcastw_dbl
!
! Purpose
! =======
!
! Do the sequence of broadcast ops.
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send/receive buffer
!
! n       (input) integer
!         number of elements in buffer
!
! iroot   (input) integer
!         rank of broadcast root
!  
! lda     (input) integer
!         leading dimension of the buffer
!
! lpx     (input) integer
!         the number of vectors in the buffer
!
! buffer  (output) real(8) array, dimension(n,*)
!         working buffer
!
! col_id  (input) integer
!         color identity
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine bcastw_dbl(buf, n, iroot, lda, lpx, buffer,
     &     col_id, icom)

      real(8), intent(inout) :: buf(1:lda,1:lpx)
      integer, intent(in)    :: n
      integer, intent(in)    :: iroot
      integer, intent(in)    :: lda
      integer, intent(in)    :: lpx
      real(8), intent(out)   :: buffer(1:n, *)
      integer, intent(in)    :: col_id
      integer, intent(in)    :: icom

      integer, parameter     :: BCASTW_TAG = 200000

      integer                :: ierr, i, j, k, tag
      integer                :: my_rank, world_size
      integer                :: local_rank, new_type
      integer, pointer       :: ireq(:)
      integer, pointer       :: rcounts(:)
      integer, pointer       :: displ(:)

#if TIMER_PRINT
      real(8)                :: timer
#endif


      if (lpx == 1) then
        call bcast_dbl(buf, n, iroot, col_id, icom)
        return
      end if

      call MPI_Comm_size(icom, world_size, ierr)
      call MPI_Comm_rank(icom, my_rank, ierr)
      my_rank = my_rank + 1

      if (lpx <= 4) then
        do i=0,lpx-1
          k = mod((iroot)-1+(i), world_size)+1
          call bcast_dbl(buf(1,i+1), n, k, col_id, icom)
        end do
        return
      end if


#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif


      allocate(rcounts(world_size))
      allocate(displ(world_size))

      do i=1,world_size
        k = mod(i-iroot+world_size, world_size)+1
        if (k <= lpx) then
          rcounts(i) = n
          displ(i)   = (k-1)*lda
        else
          rcounts(i) = 0
          displ(i)   = 0
        end if
      end do 

      k = mod(my_rank-iroot+world_size, world_size)+1
      if (k <= lpx) then
        j = n
      else
        j = 0
        k = 1
      end if

      buffer(1:n, 1) = buf(1:n, k)

      call MPI_Allgatherv(buffer, j, MPI_DOUBLE_PRECISION,
     &                    buf, rcounts, displ,
     &                    MPI_DOUBLE_PRECISION,
     &                    icom, ierr)

      deallocate(rcounts)
      deallocate(displ)

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      time_bcast = time_bcast + timer
      if (col_id >= 1 .and. col_id <= items_bcast) then
        time_bcast_(col_id)     = time_bcast_(col_id) + timer
        counter_bcast_(col_id)  = counter_bcast_(col_id) + n*lpx
        messages_bcast_(col_id) = messages_bcast_(col_id) + 1
      end if
#endif

      return

      end subroutine bcastw_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: reduce_dbl
!
! Purpose
! =======
!
! Sum up array data from all processes
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send/receive buffer
!
! wrk     (input/output) real(8) array, dimension(n)
!         working buffer
!
! n       (input) integer
!         number of elements in buffer
!  
! col_id  (input) integer
!         color identity
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine reduce_dbl(buf, wrk, n, col_id, icom)

      real(8), intent(inout) :: buf(1:n)
      real(8), intent(inout) :: wrk(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: col_id
      integer, intent(in)    :: icom

      integer                :: ierr, i, j

#if TIMER_PRINT
      real(8)                :: timer
#endif


#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif

      call MPI_Allreduce(buf, wrk, n, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, icom, ierr)
      buf(1:n) = wrk(1:n)

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      time_reduce = time_reduce + timer
      if (col_id >= 1 .and. col_id <= items_reduce) then
        time_reduce_(col_id)     = time_reduce_(col_id) + timer
        counter_reduce_(col_id)  = counter_reduce_(col_id) + n
        messages_reduce_(col_id) = messages_reduce_(col_id) + 1
      end if
#endif

      return

      end subroutine reduce_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: allgather_dbl
!
! Purpose
! =======
!
! Gather the data onto all processes
!
!
! Arguments
! =========
!
! buf     (input/output) real(8) array, dimension(n)
!         send/receive buffer
!
! wrk     (input/output) real(8) array, dimension(n)
!         working buffer
!
! n       (input) integer
!         number of elements in buffer
!  
! col_id  (input) integer
!         color identity
!
! icom    (input) integer
!         communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine allgather_dbl(buf, wrk, n, col_id, icom)

      real(8), intent(inout) :: buf(1:n)
      real(8), intent(inout) :: wrk(1:n)
      integer, intent(in)    :: n
      integer, intent(in)    :: col_id
      integer, intent(in)    :: icom

      integer                :: ierr, i, j
#if __IBM_REGISTER_VARS
      integer                :: my_rank, world_size
#endif

#if TIMER_PRINT
      real(8)                :: timer
#endif


#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif

#if __IBM_REGISTER_VARS && 0
      call MPI_Comm_size(icom, world_size, ierr)
      call MPI_Comm_rank(icom, my_rank, ierr)
      if (my_rank == 0) then
        wrk(1:n) = buf(1:n)
        do i=1,world_size-1
          call MPI_Recv(wrk(1+i*n), n, MPI_DOUBLE_PRECISION,
     &         i, 1, icom, MPI_STATUS_IGNORE, ierr)
        end do
      else
        call MPI_Send(buf(1), n, MPI_DOUBLE_PRECISION,
     &       0, 1, icom, ierr)
      end if
      call MPI_Bcast(wrk, world_size*n, MPI_DOUBLE_PRECISION,
     &     0, icom, ierr)
#else
      call MPI_Allgather(buf, n, MPI_DOUBLE_PRECISION,
     &     wrk, n, MPI_DOUBLE_PRECISION,
     &     icom, ierr)
#endif

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      time_gather = time_gather + timer
      if (col_id > 0) then
        if (col_id >= 1 .and. col_id <= items_gather) then
          time_gather_(col_id)     = time_gather_(col_id) + timer
          counter_gather_(col_id)  = counter_gather_(col_id) + n
          messages_gather_(col_id) = messages_gather_(col_id) + 1
        end if
      else
        if (-col_id >= 1 .and.
     &      -col_id <= items_gather+items_redist) then
          time_gather_(-col_id)     = time_gather_(-col_id) + timer
          counter_gather_(-col_id)  = counter_gather_(-col_id) + n
          messages_gather_(-col_id) = messages_gather_(-col_id) + 1
        end if
      end if
#endif

      return

      end subroutine allgather_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: datacast_dbl
!
! Purpose
! =======
!
! Redistribute a row-vector onto a column-vector
!
!
! Arguments
! =========
!
! u_y     (input/output) real(8) array, dimension(n)
!         recv buffer
!
! u_x     (input) real(8) array, dimension(n)
!         send buffer
!
! u_t     (input/output) real(8) array, dimension(n)
!         working buffer
!
! u_s     (input/output) real(8) array, dimension(n)
!         working buffer
!
! n       (input) integer
!         length of send/recv vector
!  
! col_id  (input) integer
!         color identity
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine datacast_dbl(u_y, u_x, u_t, u_s, n, col_id)

      real(8), intent(inout) :: u_y(1)
      real(8), intent(in)    :: u_x(1)
      real(8), intent(inout) :: u_t(1)
      real(8), intent(inout) :: u_s(1)
      integer, intent(in)    :: n
      integer, intent(in)    :: col_id

      integer                :: n_x, n_y
      integer                :: nx, ny, nz, ic, i, j, k
      integer, pointer       :: req(:)
      integer                :: reqr(2), reqs(2), x_snod, y_snod
      integer                :: his_rank, her_rank, his_local

#if TIMER_PRINT
      real(8)                :: timer
      real(8)                :: timer_excl1
      real(8)                :: timer_excl2
#endif


#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif

      if (x_nnod == 1) then
        if (y_nnod == 1) then
          u_y(1:n) = u_x(1:n)
          return
        end if
      end if

      n_x = (n-1)/x_nnod+1
      n_y = (n-1)/y_nnod+1

      if (x_nnod == 1) then
!OCL UNROLL(8)
!OCL UXSIMD
        do i=1,n_y
          j = y_inod+y_nnod*(i-1)
          u_y(i) = u_x(j)
        end do                ! i
        return
      end if

      if (x_nnod == y_nnod) then
        if (x_inod == y_inod) then
          u_y(1:n_y) = u_x(1:n_y)
        end if
        k = col_id + items_bcast
        call bcast_dbl(u_y(1), n_y, y_inod, -k, x_COMM_WORLD)
        return
      end if

      if (p0_(x_inod) > 0) then

        x_snod = x_nnod/n_common
        y_snod = y_nnod/n_common

        nx = (n_x-1)/y_snod+1
        ny = (n_x-p0_(x_inod))/y_snod+1
        nz = (n_y-q0_(x_inod))/x_snod+1

!OCL UNROLL(8)
!OCL UXSIMD
        do i=1,ny
          j = p0_(x_inod)+y_snod*(i-1)
          u_t(i) = u_x(j)
        end do
        if ( nx > ny ) then
           u_t(ny+1:nx) = ZERO
        end if

#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime()
#endif
        k = col_id + items_gather
        call allgather_dbl(u_t(1), u_s(1), nx, -k, w_COMM_WORLD)
#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime() - timer_excl2
#endif

        do ic=0,x_snod-1
          his_rank  = mod(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
          if ( n_x >= p0_(his_rank) .AND.
     &         n_y >= q0_(his_rank) ) then
            his_local = (his_rank-1)/n_common+1
            ny = (n_x-p0_(his_rank))/y_snod+1
            nz = (n_y-q0_(his_rank))/x_snod+1
!OCL UNROLL(8)
!OCL UXSIMD
            do i=1,MIN(ny,nz)
              k = q0_(his_rank)+x_snod*(i-1)
              j = (his_local-1)*nx+i
              u_y(k) = u_s(j)
            end do                ! i
          end if
        end do                  ! ic

        if ( y_inod > MOD(n-1,y_nnod)+1 ) then
          u_y(n_y) = ZERO
        end if

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        if ( z_nnod > 1 ) then
          k = col_id + items_bcast
          call bcast_dbl(u_y(1), n_y, z_inod, -k, z_COMM_WORLD)
        end if
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1
#endif

      else

        i = mod(y_inod-1,n_common)
        j = mod(x_inod-1,n_common)
        ic = mod(j-i+n_common*x_nnod*y_nnod,n_common)
        his_rank  = mod(x_inod-1 +x_nnod*y_nnod -ic,x_nnod)+1
        his_local = mod(his_rank-1, n_common)+1

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        k = col_id + items_bcast
        call bcast_dbl(u_y(1), n_y, his_local, -k, z_COMM_WORLD)
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1

        timer_excl2 = ZERO
#endif

      end if

 9999 continue

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      timer = timer - timer_excl1 - timer_excl2
      time_redist = time_redist + timer
      if (col_id >= 1 .and. col_id <= items_redist) then
        time_redist_(col_id)     = time_redist_(col_id) + timer
        counter_redist_(col_id)  = counter_redist_(col_id) + n
        messages_redist_(col_id) = messages_redist_(col_id) + 1
      end if
#endif

      return

      end subroutine datacast_dbl

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: datacast_dbl2
!
! Purpose
! =======
!
! Redistribute a couple of row-vectors onto column-vector
!
!
! Arguments
! =========
!
! ur_y     (input/output) real(8) array, dimension(n)
!          recv bufffer 1
!
! ui_y     (input/output) real(8) array, dimension(n)
!          recev buffer 2
!
! ur_x     (input) real(8) array, dimension(n)
!          send buffer 1
!
! ui_x     (input) real(8) array, dimension(n)
!          send buffer 2
!
! u_t      (input/output) real(8) array, dimension(2*n)
!          working buffer
!
! u_s      (input/output) real(8) array, dimension(2*n)
!          working buffer
!
! n        (input) integer
!          length of send/recv vector
!  
! col_id   (input) integer
!          color identity
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine datacast_dbl2(ur_y,ui_y,ur_x,ui_x,u_t,u_s, n, col_id)

      real(8), intent(inout) :: ur_y(1:n)
      real(8), intent(inout) :: ui_y(1:n)
      real(8), intent(in)    :: ur_x(1:n)
      real(8), intent(in)    :: ui_x(1:n)
      real(8), intent(inout) :: u_t(1:2*n)
      real(8), intent(inout) :: u_s(1:2*n)
      integer, intent(in)    :: n
      integer, intent(in)    :: col_id

      integer                :: n_x, n_y
      integer                :: nx, ny, nz, ic, i, j, k
      integer, pointer       :: req(:)
      integer                :: reqr(2), reqs(2), x_snod, y_snod
      integer                :: his_rank, her_rank, his_local

#if TIMER_PRINT
      real(8)                :: timer
      real(8)                :: timer_excl1
      real(8)                :: timer_excl2
#endif

#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif

      if (x_nnod == 1) then
        if (y_nnod == 1) then
          ur_y(1:n)=ur_x(1:n)
          ui_y(1:n)=ui_x(1:n)
          return
        end if
      end if

      n_x=(n-1)/x_nnod+1
      n_y=(n-1)/y_nnod+1

      if (x_nnod == 1) then
        do i=1,n_y
          j=y_inod+y_nnod*(i-1)
          ur_y(i)=ur_x(j)
          ui_y(i)=ui_x(j)
        end do                ! i
        return
      end if

      if ( x_nnod == y_nnod ) then
        if ( x_inod == y_inod ) then
          u_t(    1:    n_y) = ur_x(1:n_y)
          u_t(n_y+1:n_y+n_y) = ui_x(1:n_y)
        end if
        k = col_id + items_bcast
        call bcast_dbl(u_t(1), 2*n_y, y_inod, -k, x_COMM_WORLD)
        ur_y(1:n_y) = u_t(    1:    n_y)
        ui_y(1:n_y) = u_t(n_y+1:n_y+n_y)
        return
      end if

      if (p0_(x_inod) > 0) then

        x_snod = x_nnod/n_common
        y_snod = y_nnod/n_common

        nx = (n_x-1)/y_snod+1
        ny = (n_x-p0_(x_inod))/y_snod+1
        nz = (n_y-q0_(x_inod))/x_snod+1

        do i=1,ny
          j = p0_(x_inod)+y_snod*(i-1)
          u_t(   i) = ur_x(j)
          u_t(nx+i) = ui_x(j)
        end do
        if ( nx > ny ) then
          u_t(ny+1   :nx   ) = ZERO
          u_t(nx+ny+1:nx+nx) = ZERO
        end if

#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime()
#endif
        k = col_id + items_gather
        call allgather_dbl(u_t(1), u_s(1), 2*nx,-k, w_COMM_WORLD)
#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime() - timer_excl2
#endif

        do ic=0,x_snod-1
          his_rank  = mod(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
          if ( n_x >= p0_(his_rank) .AND.
     &         n_y >= q0_(his_rank) ) then
            his_local = (his_rank-1)/n_common+1
            ny = (n_x-p0_(his_rank))/y_snod+1
            nz = (n_y-q0_(his_rank))/x_snod+1
            do i=1,MIN(ny,nz)
              k = q0_(his_rank)+x_snod*(i-1)
              j = 2*(his_local-1)*nx+i
              ur_y(k) = u_s(   j)
              ui_y(k) = u_s(nx+j)
            end do                ! i
          end if
        end do                  ! ic

        if ( y_inod > MOD(n-1,y_nnod)+1 ) then
          ur_y(n_y) = ZERO
          ui_y(n_y) = ZERO
        end if

        u_s(    1:    n_y) = ur_y(1:n_y)
        u_s(n_y+1:n_y+n_y) = ui_y(1:n_y)

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        if ( z_nnod > 1 ) then
          k = col_id + items_bcast
          call bcast_dbl(u_s(1), 2*n_y, z_inod, -k, z_COMM_WORLD)
        end if
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1
#endif

      else

        i = mod(y_inod-1,n_common)
        j = mod(x_inod-1,n_common)
        ic = mod(j-i+n_common*x_nnod*y_nnod,n_common)
        his_rank  = mod(x_inod-1 +x_nnod*y_nnod -ic,x_nnod)+1
        his_local = mod(his_rank-1, n_common)+1

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        k = col_id + items_bcast
        call bcast_dbl(u_s(1), 2*n_y, his_local, -k, z_COMM_WORLD)
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1

        timer_excl2 = ZERO
#endif

        ur_y(1:n_y) = u_s(    1:    n_y)
        ui_y(1:n_y) = u_s(n_y+1:n_y+n_y)

      end if

 9999 continue

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      timer = timer - timer_excl1 - timer_excl2
      time_redist = time_redist + timer
      if (col_id >= 1 .and. col_id <= items_redist) then
        time_redist_(col_id)     = time_redist_(col_id) + timer
        counter_redist_(col_id)  = counter_redist_(col_id) + 2*n
        messages_redist_(col_id) = messages_redist_(col_id) + 1
      end if
#endif

      return

      end subroutine datacast_dbl2

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: datacast_dblx
!
! Purpose
! =======
!
! Redistribute multiple of row-vectors onto column-vectors
!
!
! Arguments
! =========
!
! nk       (input) integer
!          The number of buffers
!
! u_y      (input/output) real(8) array, dimension(ldv,*)
!          recv bufffer
!
! u_x      (input/output) real(8) array, dimension(ldv,*)
!          send bufffer
!
! ldv      (input) integer
!          leading dimension
!  
! u_t      (input/output) real(8) array, dimension(*)
!          working bufffer
!
! u_s      (input/output) real(8) array, dimension(*)
!          working bufffer
!
! n        (input) integer
!          length of send/recv vector
!  
! col_id   (input) integer
!          color identity
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine datacast_dblx(nk, u_y, u_x, ldv, u_t, u_s, n, col_id)

      integer, intent(in)    :: nk
      real(8), intent(inout) :: u_y(1:ldv, *)
      real(8), intent(inout) :: u_x(1:ldv, *)
      integer, intent(in)    :: ldv
      real(8), intent(inout) :: u_t(*)
      real(8), intent(inout) :: u_s(*)
      integer, intent(in)    :: n
      integer, intent(in)    :: col_id

      integer                :: n_x, n_y
      integer                :: nx, ny, nz, ic, i, j, k, kk
      integer, pointer :: req(:)
      integer :: reqr(2), reqs(2), x_snod, y_snod
      integer :: his_rank, her_rank, his_local

#if TIMER_PRINT
      real(8) :: timer
      real(8) :: timer_excl1
      real(8) :: timer_excl2
#endif

#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif

      if ( x_nnod == 1 ) then
        if ( y_nnod == 1 ) then
          do kk=1,nk
            u_y(1:n, kk)=u_x(1:n, kk)
          end do ! kk
          return
        end if
      end if

      n_x = (n-1)/x_nnod+1
      n_y = (n-1)/y_nnod+1

      if ( x_nnod == 1 ) then
        do kk=1,nk
          do i=1,n_y
            j=y_inod+y_nnod*(i-1)
            u_y(i, kk)=u_x(j, kk)
          end do              ! i
        end do              ! kk
        return
      end if

      if ( x_nnod == y_nnod ) then
        if ( x_inod == y_inod ) then
          do kk=1,nk
            j=(kk-1)*n_y
            u_t(1+j:n_y+j) = u_x(1:n_y, kk)
          end do
        end if
        k = col_id + items_bcast
        call bcast_dbl(u_t(1), n_y*nk, y_inod, -k, x_COMM_WORLD)
        do kk=1,nk
          j=(kk-1)*n_y
          u_y(1:n_y, kk) = u_t(1+j:n_y+j)
        end do ! kk
        return
      end if

      if (p0_(x_inod) > 0) then

        x_snod = x_nnod/n_common
        y_snod = y_nnod/n_common

        nx = (n_x-1)/y_snod+1
        ny = (n_x-p0_(x_inod))/y_snod+1
        nz = (n_y-q0_(x_inod))/x_snod+1

        do kk=1,nk
        do i=1,ny
          j = p0_(x_inod)+y_snod*(i-1)
          u_t(i+((kk-1)*nx)) = u_x(j, kk)
        end do
        if ( nx > ny ) then
          j = (kk-1)*nx
          u_t(ny+1+j:nx+j) = ZERO
        end if
        end do

#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime()
#endif
        k = col_id + items_gather
        call allgather_dbl(u_t(1), u_s(1), nk*nx, -k, w_COMM_WORLD)
#if TIMER_PRINT
        timer_excl2 = eigen_get_wtime() - timer_excl2
#endif

        do ic=0,x_snod-1
          his_rank  = MOD(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
          if ( n_x >= p0_(his_rank) .AND.
     &         n_y >= q0_(his_rank) ) then
            his_local = (his_rank-1)/n_common+1
            ny = (n_x-p0_(his_rank))/y_snod+1
            nz = (n_y-q0_(his_rank))/x_snod+1
            do kk=1,nk
            do i=1,MIN(ny,nz)
              k = q0_(his_rank)+x_snod*(i-1)
              j = (his_local-1)*(nk*nx)+i   +(kk-1)*nx
              u_y(k, kk) = u_s(j)
            end do              ! i
            end do              ! kk
          end if
        end do                 ! ic

        if ( y_inod > MOD(n-1,y_nnod)+1 ) then
          do kk=1,nk
            u_y(n_y,kk) = ZERO
          enddo
        end if

        if ( z_nnod > 1 ) then
          do kk=1,nk
            j = (kk-1)*n_y
            u_s(1+j:n_y+j) = u_y(1:n_y, kk)
          end do
        end if

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        if ( z_nnod > 1 ) then
          k = col_id + items_bcast
          call bcast_dbl(u_s(1), nk*n_y, z_inod, -k, z_COMM_WORLD)
        end if
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1
#endif

      else

        i = mod(y_inod-1,n_common)
        j = mod(x_inod-1,n_common)
        ic = mod(j-i+n_common*x_nnod*y_nnod,n_common)
        his_rank  = mod(x_inod-1 +x_nnod*y_nnod -ic,x_nnod)+1
        his_local = mod(his_rank-1, n_common)+1

#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime()
#endif
        k = col_id + items_bcast
        call bcast_dbl(u_s(1), nk*n_y, his_local, -k, z_COMM_WORLD)
#if TIMER_PRINT
        timer_excl1 = eigen_get_wtime() - timer_excl1

        timer_excl2 = ZERO
#endif
        do kk=1,nk
          j = (kk-1)*n_y
          u_y(1:n_y, kk) = u_s(1+j:n_y+j)
        end do

      end if

 9999 continue

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      timer = timer - timer_excl1 - timer_excl2
      time_redist = time_redist + timer
      if (col_id >= 1 .and. col_id <= items_redist) then
        time_redist_(col_id)     = time_redist_(col_id) + timer
        counter_redist_(col_id)  = counter_redist_(col_id) + 2*n
        messages_redist_(col_id) = messages_redist_(col_id) + 1
      end if
#endif

      return

      end subroutine datacast_dblx

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: print_bcast_algorithm
!
! Purpose
! =======
!
! < purpose of this subroutine ... >
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine print_bcast_algorithm()

      if (BCAST_ALGORITHM == BCAST_MPI) then
        print*, "Bcast algorithm is MPI"
      end if
      if (BCAST_ALGORITHM == BCAST_MPI_SEGMENTED) then
        print*, "Bcast algorithm is MPI_SEGMENTED"
      end if
      if (BCAST_ALGORITHM == BCAST_SEQUENTIAL) then
        print*, "Bcast algorithm is SEQUENTIAL"
      end if
      if (BCAST_ALGORITHM == BCAST_BINOMIAL) then
        print*, "Bcast algorithm is BINOMIAL"
      end if
      if (BCAST_ALGORITHM == BCAST_TRINOMIAL) then
        print*, "Bcast algorithm is TRINOMIAL"
      end if
      if (BCAST_ALGORITHM == BCAST_RELAY_ONEDIR) then
        print*, "Bcast algorithm is RELAY_ONE_DIR"
      end if
      if (BCAST_ALGORITHM == BCAST_RELAY_TWODIRS) then
        print*, "Bcast algorithm is RELAY_TWO_DIRS"
      end if

      end subroutine print_bcast_algorithm

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: print_reduce_algorithm
!
! Purpose
! =======
!
! < purpose of this subroutine ... >
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine print_reduce_algorithm()

      print*, "Reduce algorithm is MPI"

      end subroutine print_reduce_algorithm

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: print_gather_algorithm
!
! Purpose
! =======
!
! < purpose of this subroutine ... >
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine print_gather_algorithm()

      print*, "Gather algorithm is MPI"

      end subroutine print_gather_algorithm

      end module comm_mod
