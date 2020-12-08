!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!! @brief unit test for mpp_alltoall
!! @author MiKyung Lee
!! @description
!! unit test for mpp_alltoall, mpp_alltoallv, mpp_alltoallw
program test_mpp_alltoall

  use platform_mod
  use mpp_mod, only : mpp_init, mpp_init_test_requests_allocated, mpp_error, FATAL
  use mpp_mod, only : mpp_pe, mpp_npes, mpp_alltoall
  use mpp_mod, only : mpp_type_create, mpp_type, mpp_byte

  implicit none

  integer :: npes, ierr

    !> initialize MPI
    call mpp_init( test_level=mpp_init_test_requests_allocated )

    !> get total number of pe's
    npes = mpp_npes()

    !> test mpp_alltoall
    call test_mpp_alltoall_real4(npes)
    call test_mpp_alltoall_real8(npes)
    call test_mpp_alltoall_int4(npes)
    call test_mpp_alltoall_int8(npes)
    !> test mpp_alltoallv
    call test_mpp_alltoallv_real4(npes)
    call test_mpp_alltoallv_real8(npes)
    call test_mpp_alltoallv_int4(npes)
    call test_mpp_alltoallv_int8(npes)
    !> test mpp_alltoallw
    call test_mpp_alltoallw_real4(npes)
    call test_mpp_alltoallw_real8(npes)
    call test_mpp_alltoallw_int4(npes)
    call test_mpp_alltoallw_int8(npes)

    call MPI_FINALIZE(ierr)


  contains

  !>
  !> test mpp_alltoall for REAL4
  !>

  subroutine test_mpp_alltoall_real4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: srcount = 2
    real(r4_kind), parameter :: zero = 0., one = 1.

    integer :: pe, i, ii, j, N, scount, rcount
    real(r4_kind), allocatable :: sbuf(:), rbuf(:)

    !> test sending/receiving 2 elements so that, for example, for npes=4,
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoall--> [0,1,10,11,20,21,30,31]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoall--> [2,3,12,13,22,23,32,33]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoall--> [4,5,14,15,24,25,34,35]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoall--> [6,7,16,17,26,27,36,37]

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving scount/rcount array elements
    scount = srcount  ;  rcount = srcount

    !> allocate sbuf, rbuf
    N = scount*npes - 1
    allocate( sbuf(0:N), rbuf(0:N) )

    !> initialize receiving array
    rbuf = -one

    !> intialize sending array
    do i=0, N
       sbuf(i) = real( 100*pe+i, kind=r4_kind )
    end do

    !> call mpp_alltoall to send/receive one element
    call mpp_alltoall( sbuf, scount, rbuf, rcount )

    !> check
    ii = 0
    do i=0, (npes-1)
       do j=0, scount-1
          if( rbuf(ii) .ne. real( 100*i+scount*pe+j, kind=r4_kind ) ) then
             write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, ii, real(100*i+scount*pe+j), rbuf(ii)
             call mpp_error(FATAL, 'test_mpp_alltoall_real4 failed')
          end if
          ii = ii + 1
       end do
    end do

    deallocate( sbuf, rbuf )


  end subroutine test_mpp_alltoall_real4

  !>
  !> test mpp_alltoall for REAL8
  !>

  subroutine test_mpp_alltoall_real8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: srcount = 2
    real(r8_kind), parameter :: zero = 0., one = 1.

    integer :: pe, i, ii, j, N, scount, rcount
    real(r8_kind), allocatable :: sbuf(:), rbuf(:)

    !> test sending/receiving 2 elements so that, for example, for npes=4,
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoall--> [0,1,10,11,20,21,30,31]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoall--> [2,3,12,13,22,23,32,33]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoall--> [4,5,14,15,24,25,34,35]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoall--> [6,7,16,17,26,27,36,37]

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving scount/rcount array elements
    scount = srcount  ;  rcount = srcount

    !> allocate sbuf, rbuf
    N = scount*npes - 1
    allocate( sbuf(0:N), rbuf(0:N) )

    !> initialize receiving array
    rbuf = -one

    !> intialize sending array
    do i=0, N
       sbuf(i) = real( 100*pe+i, kind=r8_kind )
    end do

    !> call mpp_alltoall to send/receive one element
    call mpp_alltoall( sbuf, scount, rbuf, rcount )

    !> check
    ii = 0
    do i=0, (npes-1)
       do j=0, scount-1
          if( rbuf(ii) .ne. real( 100*i+scount*pe+j, kind=r8_kind ) ) then
             write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, ii, real(100*i+scount*pe+j), rbuf(ii)
             call mpp_error(FATAL, 'test_mpp_alltoall_real8 failed')
          end if
          ii = ii + 1
       end do
    end do

    deallocate( sbuf, rbuf )


  end subroutine test_mpp_alltoall_real8

  !>
  !> test mpp_alltoall for INT4
  !>

  subroutine test_mpp_alltoall_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: srcount = 2
    integer(i4_kind), parameter :: zero = 0, one = 1

    integer :: pe, i, ii, j, N, scount, rcount
    integer(i4_kind), allocatable :: sbuf(:), rbuf(:)

    !> test sending/receiving 2 elements so that, for example, for npes=4,
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoall--> [0,1,10,11,20,21,30,31]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoall--> [2,3,12,13,22,23,32,33]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoall--> [4,5,14,15,24,25,34,35]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoall--> [6,7,16,17,26,27,36,37]

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving scount/rcount array elements
    scount = srcount  ;  rcount = srcount

    !> allocate sbuf, rbuf
    N = scount*npes - 1
    allocate( sbuf(0:N), rbuf(0:N) )

    !> initialize receiving array
    rbuf = -one

    !> intialize sending array
    do i=0, N
       sbuf(i) = int( 100*pe+i, kind=i4_kind )
    end do

    !> call mpp_alltoall to send/receive one element
    call mpp_alltoall( sbuf, scount, rbuf, rcount )

    !> check
    ii = 0
    do i=0, (npes-1)
       do j=0, scount-1
          if( rbuf(ii) .ne. int( 100*i+scount*pe+j, kind=i4_kind ) ) then
             write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, ii, int(100*i+scount*pe+j), rbuf(ii)
             call mpp_error(FATAL, 'test_mpp_alltoall_int4 failed')
          end if
          ii = ii + 1
       end do
    end do

    deallocate( sbuf, rbuf )

  end subroutine test_mpp_alltoall_int4

  !>
  !> test mpp_alltoall for INT8
  !>

  subroutine test_mpp_alltoall_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: srcount = 2
    integer(i8_kind), parameter :: zero = 0, one = 1

    integer :: pe, i, ii, j, N, scount, rcount
    integer(i8_kind), allocatable :: sbuf(:), rbuf(:)

    !> test sending/receiving 2 elements so that, for example, for npes=4,
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoall--> [0,1,10,11,20,21,30,31]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoall--> [2,3,12,13,22,23,32,33]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoall--> [4,5,14,15,24,25,34,35]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoall--> [6,7,16,17,26,27,36,37]

    !> get pe
    pe = mpp_pe()

    !> test sending/receiving scount/rcount array elements
    scount = srcount  ;  rcount = srcount

    !> allocate sbuf, rbuf
    N = scount*npes - 1
    allocate( sbuf(0:N), rbuf(0:N) )

    !> initialize receiving array
    rbuf = -one

    !> intialize sending array
    do i=0, N
       sbuf(i) = int( 100*pe+i, kind=i8_kind )
    end do

    !> call mpp_alltoall to send/receive one element
    call mpp_alltoall( sbuf, scount, rbuf, rcount )

    !> check
    ii = 0
    do i=0, (npes-1)
       do j=0, scount-1
          if( rbuf(ii) .ne. int( 100*i+scount*pe+j, kind=i8_kind ) ) then
             write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, ii, int(100*i+scount*pe+j), rbuf(ii)
             call mpp_error(FATAL, 'test_mpp_alltoall_int8 failed')
          end if
          ii = ii + 1
       end do
    end do

    deallocate( sbuf, rbuf )

  end subroutine test_mpp_alltoall_int8

  !>
  !> test mpp_alltoallv for REAL4
  !>

  subroutine test_mpp_alltoallv_real4(npes)

    implicit none

    integer, intent(in) :: npes

    real(r4_kind) :: zero = 0., one = 1.

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(r4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N  = npes - 1

    !> allocate arrays required for mpp_alltoallv
    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    !> send/receive 1 element
    ssize = 1  ;  rsize = 1

    !> contiguous send/receive
    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    !> assign sbuf/rbuf to send/receive
    do i=0, N
       sbuf(i) = real( 100*pe+i, kind=r4_kind )
    end do
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(100*i+pe, kind=r4_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, i, real(100*i+pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real4 fail in sending/receiving 1 element' )
       end if
    end do

    deallocate( sbuf, rbuf )

    !>send every other element and receive
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    !> send/receive one element
    ssize = 1  ;  rsize = 1

    !> send/receive for every even indices
    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    !> assign sbuf to send
    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i,   kind=r4_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=r4_kind )
    end do

    !> initialize rbuf to receive
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=r4_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, i, real(10*i+2*pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real4 fail sending/receiving every other element' )
       end if
       if ( rbuf(2*i+1).ne.-one ) then
          write(*,'("PE #",i3,"element",i4,"Expected",f6.0,"but received",f6.0)') pe, i,-one, rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real4 fail sending/receiving every other element' )
       end if
    end do

    deallocate( sbuf, rbuf)

  end subroutine test_mpp_alltoallv_real4

  !>
  !> test mpp_alltoallv for REAL8
  !>

  subroutine test_mpp_alltoallv_real8(npes)

    implicit none

    integer, intent(in) :: npes

    real(r8_kind) :: zero = 0., one = 1.

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    real(r8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N  = npes - 1

    !> allocate arrays required for mpp_alltoallv
    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    !> send/receive 1 element
    ssize = 1  ;  rsize = 1

    !> contiguous send/receive
    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    !> assign sbuf/rbuf to send/receive
    do i=0, N
       sbuf(i) = real( 100*pe+i, kind=r8_kind )
    end do
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(i).ne.real(100*i+pe, kind=r8_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, i, real(100*i+pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real8 fail in sending/receiving 1 element' )
       end if
    end do

    deallocate( sbuf, rbuf )

    !>send every other element and receive
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    !> send/receive one element
    ssize = 1  ;  rsize = 1

    !> send/receive for every even indices
    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    !> assign sbuf to send
    do i=0, N
       sbuf(2*i)   = real( 10*pe+2*i,   kind=r8_kind )
       sbuf(2*i+1) = real( 10*pe+2*i+1, kind=r8_kind )
    end do

    !> initialize rbuf to receive
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.real(10*i+2*pe, kind=r8_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, i, real(10*i+2*pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real8 fail sending/receiving every other element' )
       end if
       if ( rbuf(2*i+1).ne.-one ) then
          write(*,'("PE # ",i3," element",i4," Expected",f6.0," but received",f6.0)') pe, i,-one, rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_real8 fail sending/receiving every other element' )
       end if
    end do

    deallocate( sbuf, rbuf)

  end subroutine test_mpp_alltoallv_real8

  !>
  !> test mpp_alltoallv for INT4
  !>

  subroutine test_mpp_alltoallv_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i4_kind) :: zero = 0, one = 1

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    integer(i4_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N  = npes - 1

    !> allocate arrays required for mpp_alltoallv
    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    !> send/receive 1 element
    ssize = 1  ;  rsize = 1

    !> contiguous send/receive
    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    !> assign sbuf/rbuf to send/receive
    do i=0, N
       sbuf(i) = int( 100*pe+i, kind=i4_kind )
    end do
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(i).ne.int(100*i+pe, kind=i4_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i, int(100*i+pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int4 fail in sending/receiving 1 element' )
       end if
    end do

    deallocate( sbuf, rbuf )

    !>send every other element and receive
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    !> send/receive one element
    ssize = 1  ;  rsize = 1

    !> send/receive for every even indices
    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    !> assign sbuf to send
    do i=0, N
       sbuf(2*i)   = int( 10*pe+2*i,   kind=i4_kind )
       sbuf(2*i+1) = int( 10*pe+2*i+1, kind=i4_kind )
    end do

    !> initialize rbuf to receive
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.int(10*i+2*pe, kind=i4_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i, int(10*i+2*pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int4 fail sending/receiving every other element' )
       end if
       if ( rbuf(2*i+1).ne.-one ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i,-one, rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int4 fail sending/receiving every other element' )
       end if
    end do

    deallocate( sbuf, rbuf)

  end subroutine test_mpp_alltoallv_int4

  !>
  !> test mpp_alltoallv for INT8
  !>

  subroutine test_mpp_alltoallv_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer(i8_kind) :: zero = 0, one = 1

    integer :: pe, ierr, i, ii, N
    integer, allocatable :: ssize(:), rsize(:), sdispl(:), rdispl(:)
    integer(i8_kind), allocatable :: sbuf(:), rbuf(:)

    !> get pe
    pe = mpp_pe()
    N  = npes - 1

    !> allocate arrays required for mpp_alltoallv
    allocate( sbuf(0:N),   rbuf(0:N) )
    allocate( ssize(0:N),  rsize(0:N) )
    allocate( sdispl(0:N), rdispl(0:N) )

    !>send one, receive one
    !! process0: [ 0, 1, 2, 3]  --alltoallv--> [0,10,20,30]
    !! process1: [10,11,12,13]  --alltoallv--> [1,11,21,31]
    !! process2: [20,21,22,23]  --alltoallv--> [2,12,22,32]
    !! process3: [30,31,32,33]  --alltoallv--> [3,13,23,33]

    !> send/receive 1 element
    ssize = 1  ;  rsize = 1

    !> contiguous send/receive
    do i=0, N
       sdispl(i) = i  ;  rdispl(i) = i
    end do

    !> assign sbuf/rbuf to send/receive
    do i=0, N
       sbuf(i) = int( 100*pe+i, kind=i8_kind )
    end do
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(i).ne.int(100*i+pe, kind=i8_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i, int(100*i+pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int8 fail in sending/receiving 1 element' )
       end if
    end do

    deallocate( sbuf, rbuf )

    !>send every other element and receive
    !! process0: [ 0, 1, 2, 3, 4, 5, 6, 7]  --alltoallv--> [0,-1,10,-1,20,-1,30,-1]
    !! process1: [10,11,12,13,14,15,16,17]  --alltoallv--> [2,-1,12,-1,22,-1,32,-1]
    !! process2: [20,21,22,23,24,25,26,27]  --alltoallv--> [4,-1,14,-1,24,-1,34,-1]
    !! process3: [30,31,32,33,34,35,36,37]  --alltoallv--> [6,-1,16,-1,26,-1,36,-1]

    allocate( sbuf(0:(2*npes-1)), rbuf(0:(2*npes-1)) )

    !> send/receive one element
    ssize = 1  ;  rsize = 1

    !> send/receive for every even indices
    do i=0, N
       sdispl(i) = 2*i  ;  rdispl(i) = 2*i
    end do

    !> assign sbuf to send
    do i=0, N
       sbuf(2*i)   = int( 10*pe+2*i,   kind=i8_kind )
       sbuf(2*i+1) = int( 10*pe+2*i+1, kind=i8_kind )
    end do

    !> initialize rbuf to receive
    rbuf = -one

    !> call mpp_alltoallv
    call mpp_alltoall(sbuf, ssize, sdispl, rbuf, rsize, rdispl)

    !> check
    do i=0, N
       if ( rbuf(2*i).ne.int(10*i+2*pe, kind=i8_kind) ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i, int(10*i+2*pe), rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int8 fail sending/receiving every other element' )
       end if
       if ( rbuf(2*i+1).ne.-one ) then
          write(*,'("PE # ",i3," element",i4," Expected",i6," but received",i6)') pe, i,-one, rbuf(i)
          call mpp_error( FATAL, 'test_mpp_alltoallv_int8 fail sending/receiving every other element' )
       end if
    end do

    deallocate( sbuf, rbuf)

  end subroutine test_mpp_alltoallv_int8

  !>
  !> test mpp_alltoallw REAL4
  !>

  subroutine test_mpp_alltoallw_real4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n = 9
    integer, parameter :: byte4 = 4
    real(r4_kind), parameter :: zero = 0. , one = 1.

    integer :: pe, i, j, jj, jjj, k, kk
    real(r4_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    real(r4_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    real(r4_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)

    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = real( pe*1000 + i*100 + j*10 + k, kind=r4_kind )
          end do
       end do
    end do
    rbuf = - one

    !> each PE send/receive every other column of length subsize_k=2, and starting from sbuf(start_k=2,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, j=1, kth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, j=1 : [ -1, -1,  112,  113 ]  ; i=1, j=2 : [-1,-1,-1,-1]
    !! i=1, j=3 : [ -1, -1, 1112, 1113 ]  ; i=1, j=4 : [-1,-1,-1,-1]
    !! i=1, j=5 : [ -1, -1, 2112, 2113 ]  ; i=1, j=6 : [-1,-1,-1,-1]
    !! i=1, j=7 : [ -1, -1, 3112, 3114 ]  ; i=1, j=8 : [-1,-1,-1,-1]

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    !> send every other kth column (hence the 2)
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    !> receive in every other kth column (hence the 2)
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=r4_kind )
             if ( i.gt.subsize_i )                           answer=-one !! rbuf = -1 for i>1
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one !! rbuf(k',2pe,1)=rbuf(k',j>2npes,1)=-1
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one !! rbuf(k<start_k,2pe+1,1)=rbuf(k>start_k+subsize_k,2pe+1,1)=-1
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real4 with columns' )
             end if
          end do
       end do
    end do

    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, k=1, jth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, k=1 : [ -1, -1,  112,  113 ]
    !! i=1, k=2 : [ -1, -1, 1112, 1113 ]
    !! i=1, k=3 : [ -1, -1, 2112, 2113 ]
    !! i=1, k=4 : [ -1, -1, 3112, 3114 ]

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=real( (k-1)*1000 + i*100 + j*10 + pe+1, kind=r4_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real4 with rows' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r4_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real4 with rank 2 subarrays' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r4_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real4 with rank 3 subarrays' )
             end if
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_real4

  !>
  !> test mpp_alltoallw REAL8
  !>

  subroutine test_mpp_alltoallw_real8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n = 9
    integer, parameter :: byte8 = 8
    real(r8_kind), parameter :: zero = 0. , one = 1.

    integer :: pe, i, j, jj, jjj, k, kk
    real(r8_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    real(r8_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    real(r8_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)

    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = real( pe*1000 + i*100 + j*10 + k, kind=r8_kind )
          end do
       end do
    end do
    rbuf = - one

    !> each PE send/receive every other column of length subsize_k=2, and starting from sbuf(start_k=2,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, j=1, kth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, j=1 : [ -1, -1,  112,  113 ]  ; i=1, j=2 : [-1,-1,-1,-1]
    !! i=1, j=3 : [ -1, -1, 1112, 1113 ]  ; i=1, j=4 : [-1,-1,-1,-1]
    !! i=1, j=5 : [ -1, -1, 2112, 2113 ]  ; i=1, j=6 : [-1,-1,-1,-1]
    !! i=1, j=7 : [ -1, -1, 3112, 3114 ]  ; i=1, j=8 : [-1,-1,-1,-1]

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    !> send every other kth column (hence the 2)
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    !> receive in every other kth column (hence the 2)
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=real( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=r8_kind )
             if ( i.gt.subsize_i )                           answer=-one !! rbuf = -1 for i>1
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one !! rbuf(k',2pe,1)=rbuf(k',j>2npes,1)=-1
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one !! rbuf(k<start_k,2pe+1,1)=rbuf(k>start_k+subsize_k,2pe+1,1)=-1
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real8 with columns' )
             end if
          end do
       end do
    end do

    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, k=1, jth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, k=1 : [ -1, -1,  112,  113 ]
    !! i=1, k=2 : [ -1, -1, 1112, 1113 ]
    !! i=1, k=3 : [ -1, -1, 2112, 2113 ]
    !! i=1, k=4 : [ -1, -1, 3112, 3114 ]

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=real( (k-1)*1000 + i*100 + j*10 + pe+1, kind=r8_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real8 with rows' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r8_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real8 with rank 2 subarrays' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=real( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=r8_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",f6.0," but received ",f6.0)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_real8 with rank 3 subarrays' )
             end if
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_real8

  !>
  !> test mpp_alltoallw INT4
  !>

  subroutine test_mpp_alltoallw_int4(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n = 9
    integer, parameter :: byte4 = 4
    integer(i4_kind), parameter :: zero = 0 , one = 1

    integer :: pe, i, j, jj, jjj, k, kk
    integer(i4_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    integer(i4_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    integer(i4_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)

    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = int( pe*1000 + i*100 + j*10 + k, kind=i4_kind )
          end do
       end do
    end do
    rbuf = - one

    !> each PE send/receive every other column of length subsize_k=2, and starting from sbuf(start_k=2,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, j=1, kth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, j=1 : [ -1, -1,  112,  113 ]  ; i=1, j=2 : [-1,-1,-1,-1]
    !! i=1, j=3 : [ -1, -1, 1112, 1113 ]  ; i=1, j=4 : [-1,-1,-1,-1]
    !! i=1, j=5 : [ -1, -1, 2112, 2113 ]  ; i=1, j=6 : [-1,-1,-1,-1]
    !! i=1, j=7 : [ -1, -1, 3112, 3114 ]  ; i=1, j=8 : [-1,-1,-1,-1]

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    !> send every other kth column (hence the 2)
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    !> receive in every other kth column (hence the 2)
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=int( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=i4_kind )
             if ( i.gt.subsize_i )                           answer=-one !! rbuf = -1 for i>1
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one !! rbuf(k',2pe,1)=rbuf(k',j>2npes,1)=-1
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one !! rbuf(k<start_k,2pe+1,1)=rbuf(k>start_k+subsize_k,2pe+1,1)=-1
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE #",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int4 with columns' )
             end if
          end do
       end do
    end do

    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, k=1, jth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, k=1 : [ -1, -1,  112,  113 ]
    !! i=1, k=2 : [ -1, -1, 1112, 1113 ]
    !! i=1, k=3 : [ -1, -1, 2112, 2113 ]
    !! i=1, k=4 : [ -1, -1, 3112, 3114 ]

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=int( (k-1)*1000 + i*100 + j*10 + pe+1, kind=i4_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE #",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int4 with rows' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i4_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int4 with rank 2 subarrays' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte4
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte4
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i4_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int4 with rank 3 subarrays' )
             end if
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_int4

  !>
  !> test mpp_alltoallw INT8
  !>

  subroutine test_mpp_alltoallw_int8(npes)

    implicit none

    integer, intent(in) :: npes

    integer, parameter :: n = 9
    integer, parameter :: byte8 = 8
    integer(i8_kind), parameter :: zero = 0 , one = 1

    integer :: pe, i, j, jj, jjj, k, kk
    integer(i8_kind) :: answer

    integer :: array_of_subsizes(3), array_of_starts(3)
    integer :: subsize_i, subsize_j, subsize_k
    integer :: start_i, start_j, start_k
    integer :: ssize(0:npes-1), rsize(0:npes-1), sdispl(0:npes-1), rdispl(0:npes-1)
    integer(i8_kind), target :: sbuf(n,n,n), rbuf(n,n,n)

    integer(i8_kind), dimension(:), pointer :: psbuf, prbuf
    type(mpp_type) :: stype(0:npes-1), rtype(0:npes-1)

    !> get pe
    pe = mpp_pe()

    !> assign sbuf and rbuf
    do i=1, n
       do j=1, n
          do k=1, n
             sbuf(k,j,i) = int( pe*1000 + i*100 + j*10 + k, kind=i8_kind )
          end do
       end do
    end do
    rbuf = - one

    !> each PE send/receive every other column of length subsize_k=2, and starting from sbuf(start_k=2,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, j=1, kth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, j=1 : [ -1, -1,  112,  113 ]  ; i=1, j=2 : [-1,-1,-1,-1]
    !! i=1, j=3 : [ -1, -1, 1112, 1113 ]  ; i=1, j=4 : [-1,-1,-1,-1]
    !! i=1, j=5 : [ -1, -1, 2112, 2113 ]  ; i=1, j=6 : [-1,-1,-1,-1]
    !! i=1, j=7 : [ -1, -1, 3112, 3114 ]  ; i=1, j=8 : [-1,-1,-1,-1]

    !> subarray dimensions
    subsize_k = 5  ;  subsize_j = 1  ;  subsize_i = 1
    start_k   = 3  ;  start_j   = 0  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    !> send every other kth column (hence the 2)
    do i=0, npes-1
       sdispl(i) = 2 * i * n * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    !> receive in every other kth column (hence the 2)
    do i=0, npes-1
       rdispl(i) = 2 * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    do i=1, n
       do j=1, n
          jj = int( (j-1)/2 )
          do k=1, n
             answer=int( jj*1000 + i*100 + (2*pe+1)*10 + k, kind=i8_kind )
             if ( i.gt.subsize_i )                           answer=-one !! rbuf = -1 for i>1
             if ( mod(j,2).eq.0 .or. j.gt.2*npes )           answer=-one !! rbuf(k',2pe,1)=rbuf(k',j>2npes,1)=-1
             if ( k.le.start_k .or. k.gt.subsize_k+start_k ) answer=-one !! rbuf(k<start_k,2pe+1,1)=rbuf(k>start_k+subsize_k,2pe+1,1)=-1
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE #",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int8 with columns' )
             end if
          end do
       end do
    end do

    !> test each PE sending a row of length subsize_i, and starting from sbuf(start_k,start_j,start_i)
    !> EXAMPLE
    !! sbuf for i=1, k=1, jth array
    !! process0: [ 110, 111, 112, 113]
    !! process1: [1110,1111,1112,1113]
    !! process2: [2110,2111,2112,2113]
    !! process3: [3110,3111,3112,3113]
    !! rbuf for pe 0 :
    !! i=1, k=1 : [ -1, -1,  112,  113 ]
    !! i=1, k=2 : [ -1, -1, 1112, 1113 ]
    !! i=1, k=3 : [ -1, -1, 2112, 2113 ]
    !! i=1, k=4 : [ -1, -1, 3112, 3114 ]

    rbuf = - one

    !> subarray dimensions
    subsize_k = 1  ;  subsize_j = 5  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 2  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) =  i * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k, start_j, start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          do k=1, n
             answer=int( (k-1)*1000 + i*100 + j*10 + pe+1, kind=i8_kind )
             if ( i .gt. subsize_i )                         answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j+start_j ) answer=-one
             if ( k .gt. npes )                              answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE #",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int8 with rows' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 2
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 1
    start_k   = 0  ;  start_j   = 1  ;  start_i   = 0

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i8_kind )
             if ( i .gt. subsize_i )                                answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )   answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)         answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int8 with rank 2 subarrays' )
             end if
          end do
       end do
    end do

    !>
    !> send and receive subarray of rank 3
    !>

    rbuf = -one

    !> subarray dimensions
    subsize_k = 2  ;  subsize_j = 2  ;  subsize_i = 2
    start_k   = 1  ;  start_j   = 1  ;  start_i   = 1

    !> send one group to each PE
    ssize = 1
    do i=0, npes-1
       sdispl(i) = i * subsize_k * byte8
    end do

    !> receive one group from each PE
    rsize = 1
    do i=0, npes-1
       rdispl(i) = subsize_j * i * n * byte8
    end do

    !> subarrays (portion of data) in sbuf/rbuf to send/receive
    array_of_subsizes=(/subsize_k, subsize_j, subsize_i/)
    array_of_starts=(/start_k,start_j,start_i/)

    !> initialize mpp_type datatype
    stype(:) = mpp_byte  ;  rtype(:) = mpp_byte

    !> create mpp_type datatype
    do i=0, npes-1
       call mpp_type_create( sbuf, array_of_subsizes, array_of_starts, stype(i) )
       call mpp_type_create( rbuf, array_of_subsizes, array_of_starts, rtype(i) )
    end do

    !> mpp_alltoallW
    psbuf(1:size(sbuf)) => sbuf  ;  prbuf(1:size(rbuf)) => rbuf
    call mpp_alltoall( psbuf, ssize, sdispl, stype, prbuf, rsize, rdispl, stype )

    !> check
    do i=1, n
       do j=1, n
          jj  = int( (j-1-start_j)/subsize_j )
          jjj = mod( (j-1-start_j), subsize_j ) + 1 + start_j
          do k=1, n
             answer=int( jj*1000 + i*100 + jjj*10 + subsize_k*pe+k, kind=i8_kind )
             if ( i.le.start_i .or. i.gt.subsize_i+start_i )       answer=-one
             if ( j.le.start_j .or. j.gt.subsize_j*npes+start_j )  answer=-one
             if ( k.le.start_k .or. k.gt.subsize_k+start_k)        answer=-one
             if ( rbuf(k,j,i) .ne. answer ) then
                write(*,'("PE # ",i3," element",i2,",",i2,",",i2, " Expected",i6," but received ",i6)') pe, k,j,i, answer, rbuf(k,j,i)
                call mpp_error( FATAL, 'error in test_mpp_alltoallw_int8 with rank 3 subarrays' )
             end if
          end do
       end do
    end do

  end subroutine test_mpp_alltoallw_int8

end program test_mpp_alltoall
