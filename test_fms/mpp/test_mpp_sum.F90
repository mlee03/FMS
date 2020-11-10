
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
!! @author Lauren Chilutti
!! @brief Test program for the mpp_sum interface.
!! @email gfdl.climate.model.info@noaa.gov
!! @description This test program is for testing the mpp_sum interface.

program test_mpp_sum

  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod, only : mpp_sync
  use mpp_mod, only : mpp_set_stack_size, mpp_init_test_requests_allocated
  use mpp_mod, only : mpp_sum
  use mpp_mod, only : mpp_error, FATAL
  use platform_mod

  implicit none

  integer                                       :: n, ierr
  real(kind=r4_kind), allocatable, dimension(:) :: a4, b4, c4
  real(kind=r8_kind), allocatable, dimension(:) :: a8, b8, c8
  integer                                       :: i, pe, npes, root
  integer, allocatable, dimension(:)            :: pelist

  call mpp_init(mpp_init_test_requests_allocated)
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()

  if( pe.EQ.root ) print *, '------------------> Calling test_mpp_sum <------------------'
    call test_mpp_sum_scalar(pe,npes,root)
    call test_mpp_sum_2D(pe,npes,root)
    call test_mpp_sum_3D(pe,npes,root)
    call test_mpp_sum_4D(pe,npes,root)
    call test_mpp_sum_5D(pe,npes,root)
    call test_mpp_sum_pelist(pe,npes,root)
    call test_mpp_sum_large(pe,npes,root)
  if( pe.EQ.root ) print *, '------------------> Finished test_mpp_sum <------------------'

  call MPI_FINALIZE(ierr)

contains

  subroutine test_mpp_sum_scalar(pe,npes,root)
    integer, intent(in) :: npes, pe, root

    call test_mpp_sum_scalar_r4(pe,npes,root)
    call test_mpp_sum_scalar_r8(pe,npes,root)
    call test_mpp_sum_scalar_i4(pe,npes,root)
    call test_mpp_sum_scalar_i8(pe,npes,root)

  end subroutine test_mpp_sum_scalar

  !> Test the functionality of mpp_transmit for an r4_scalar.
  subroutine test_mpp_sum_scalar_r4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=1
    real(kind=r4_kind) :: a4, c4
    real(kind=r4_kind), dimension(npes) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    call mpp_sum(a4)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (a4 .ne. c4) call mpp_error(FATAL, "Scalar_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_r4

  !> Test the functionality of mpp_transmit for an r8_scalar.
  subroutine test_mpp_sum_scalar_r8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=1
    real(kind=r8_kind) :: a8, c8
    real(kind=r8_kind), dimension(npes) :: b8

    a8 = real(pe+1, kind=r8_kind)
    call mpp_sync()
    call mpp_sum(a8)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (a8 .ne. c8) call mpp_error(FATAL, "Scalar_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_r8

  !> Test the functionality of mpp_transmit for an i4_scalar.
  subroutine test_mpp_sum_scalar_i4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=1
    integer(kind=i4_kind) :: a4, c4
    integer(kind=i4_kind), dimension(npes) :: b4

    a4 = int(pe+1, kind=i4_kind)
    call mpp_sync()
    call mpp_sum(a4)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (a4 .ne. c4) call mpp_error(FATAL, "Scalar_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_i4

  !> Test the functionality of mpp_transmit for an i8_scalar.
  subroutine test_mpp_sum_scalar_i8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=1
    integer(kind=i8_kind) :: a8, c8
    integer(kind=i8_kind), dimension(npes) :: b8

    a8 = int(pe+1, kind=i8_kind)
    call mpp_sync()
    call mpp_sum(a8)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (a8 .ne. c8) call mpp_error(FATAL, "Scalar_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_i8

  subroutine test_mpp_sum_2D(pe,npes,root)
    integer, intent(in) :: npes, pe, root

    call test_mpp_sum_2D_r4(pe,npes,root)
    call test_mpp_sum_2D_r8(pe,npes,root)
    call test_mpp_sum_2D_i4(pe,npes,root)
    call test_mpp_sum_2D_i8(pe,npes,root)

  end subroutine test_mpp_sum_2D

  !> Test the functionality of mpp_transmit for an r4_2D.
  subroutine test_mpp_sum_2D_r4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=4
    real(kind=r4_kind), dimension(2,2) :: a4, c4
    real(kind=r4_kind), dimension(npes) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "2D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_r4

  !> Test the functionality of mpp_transmit for an r8_2D.
  subroutine test_mpp_sum_2D_r8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=4
    real(kind=r8_kind), dimension(2,2) :: a8, c8
    real(kind=r8_kind), dimension(npes) :: b8

    a8 = real(pe+1, kind=r8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "2D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_r8

  !> Test the functionality of mpp_transmit for an i4_2D.
  subroutine test_mpp_sum_2D_i4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=4
    integer(kind=i4_kind), dimension(2,2) :: a4, c4
    integer(kind=i4_kind), dimension(npes) :: b4

    a4 = int(pe+1, kind=i4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "2D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_i4

  !> Test the functionality of mpp_transmit for an i8_2D.
  subroutine test_mpp_sum_2D_i8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=4
    integer(kind=i8_kind), dimension(2,2) :: a8, c8
    integer(kind=i8_kind), dimension(npes) :: b8

    a8 = int(pe+1, kind=i8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "2D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_i8

  subroutine test_mpp_sum_3D(pe,npes,root)
    integer, intent(in) :: npes, pe, root

    call test_mpp_sum_3D_r4(pe,npes,root)
    call test_mpp_sum_3D_r8(pe,npes,root)
    call test_mpp_sum_3D_i4(pe,npes,root)
    call test_mpp_sum_3D_i8(pe,npes,root)

  end subroutine test_mpp_sum_3D

  !> Test the functionality of mpp_transmit for an r4_3D.
  subroutine test_mpp_sum_3D_r4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=8
    real(kind=r4_kind), dimension(2,2,2) :: a4, c4
    real(kind=r4_kind), dimension(npes) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "3D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_r4

  !> Test the functionality of mpp_transmit for an r8_3D.
  subroutine test_mpp_sum_3D_r8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=8
    real(kind=r8_kind), dimension(2,2,2) :: a8, c8
    real(kind=r8_kind), dimension(npes) :: b8

    a8 = real(pe+1, kind=r8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "3D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_r8

  !> Test the functionality of mpp_transmit for an i4_3D.
  subroutine test_mpp_sum_3D_i4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=8
    integer(kind=i4_kind), dimension(2,2,2) :: a4, c4
    integer(kind=i4_kind), dimension(npes) :: b4

    a4 = int(pe+1, kind=i4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "3D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_i4

  !> Test the functionality of mpp_transmit for an i8_3D.
  subroutine test_mpp_sum_3D_i8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=8
    integer(kind=i8_kind), dimension(2,2,2) :: a8, c8
    integer(kind=i8_kind), dimension(npes) :: b8

    a8 = int(pe+1, kind=i8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "3D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_i8

  subroutine test_mpp_sum_4D(pe,npes,root)
    integer, intent(in) :: npes, pe, root

    call test_mpp_sum_4D_r4(pe,npes,root)
    call test_mpp_sum_4D_r8(pe,npes,root)
    call test_mpp_sum_4D_i4(pe,npes,root)
    call test_mpp_sum_4D_i8(pe,npes,root)

  end subroutine test_mpp_sum_4D

  !> Test the functionality of mpp_transmit for an r4_4D.
  subroutine test_mpp_sum_4D_r4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=16
    real(kind=r4_kind), dimension(2,2,2,2) :: a4, c4
    real(kind=r4_kind), dimension(npes) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "4D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_r4

  !> Test the functionality of mpp_transmit for an r8_4D.
  subroutine test_mpp_sum_4D_r8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=16
    real(kind=r8_kind), dimension(2,2,2,2) :: a8, c8
    real(kind=r8_kind), dimension(npes) :: b8

    a8 = real(pe+1, kind=r8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "4D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_r8

  !> Test the functionality of mpp_transmit for an i4_4D.
  subroutine test_mpp_sum_4D_i4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=16
    integer(kind=i4_kind), dimension(2,2,2,2) :: a4, c4
    integer(kind=i4_kind), dimension(npes) :: b4

    a4 = int(pe+1, kind=i4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "4D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_i4

  !> Test the functionality of mpp_transmit for an i8_4D.
  subroutine test_mpp_sum_4D_i8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=16
    integer(kind=i8_kind), dimension(2,2,2,2) :: a8, c8
    integer(kind=i8_kind), dimension(npes) :: b8

    a8 = int(pe+1, kind=i8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "4D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_i8

  subroutine test_mpp_sum_5D(pe,npes,root)
    integer, intent(in) :: npes, pe, root

    call test_mpp_sum_5D_r4(pe,npes,root)
    call test_mpp_sum_5D_r8(pe,npes,root)
    call test_mpp_sum_5D_i4(pe,npes,root)
    call test_mpp_sum_5D_i8(pe,npes,root)

  end subroutine test_mpp_sum_5D

  !> Test the functionality of mpp_transmit for an r4_5D.
  subroutine test_mpp_sum_5D_r4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=32
    real(kind=r4_kind), dimension(2,2,2,2,2) :: a4, c4
    real(kind=r4_kind), dimension(npes) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "5D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_r4

  !> Test the functionality of mpp_transmit for an r8_5D.
  subroutine test_mpp_sum_5D_r8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=32
    real(kind=r8_kind), dimension(2,2,2,2,2) :: a8, c8
    real(kind=r8_kind), dimension(npes) :: b8

    a8 = real(pe+1, kind=r8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "5D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_r8

  !> Test the functionality of mpp_transmit for an i4_5D.
  subroutine test_mpp_sum_5D_i4(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=32
    integer(kind=i4_kind), dimension(2,2,2,2,2) :: a4, c4
    integer(kind=i4_kind), dimension(npes) :: b4

    a4 = int(pe+1, kind=i4_kind)
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes, 1)/)
      c4 = sum(b4)
      if (all(a4 .ne. c4)) call mpp_error(FATAL, "5D_r4: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_i4

  !> Test the functionality of mpp_transmit for an i8_5D.
  subroutine test_mpp_sum_5D_i8(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=32
    integer(kind=i8_kind), dimension(2,2,2,2,2) :: a8, c8
    integer(kind=i8_kind), dimension(npes) :: b8

    a8 = int(pe+1, kind=i8_kind)
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .EQ. root) then
      b8 = (/(i, i=1,npes, 1)/)
      c8 = sum(b8)
      if (all(a8 .ne. c8)) call mpp_error(FATAL, "5D_r8: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_i8

  !> Test mpp_sum for a select list of pes with r4
  subroutine test_mpp_sum_pelist(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=1
    integer, dimension(0:npes-2) :: pelist
    real(kind=r4_kind), dimension(1) :: a4, c4
    real(kind=r4_kind), dimension(0:npes-2) :: b4

    a4 = real(pe+1, kind=r4_kind)
    call mpp_sync()
    if (pe .LE. npes-2) then
      pelist = (/(i,i=0,npes-2)/)
      call mpp_sum(a4(1:n),n, pelist=pelist)
    endif
    if (pe .EQ. root) then
      b4 = (/(i, i=1,npes-1, 1)/)
      c4 = sum(b4)
      if (a4(1) .ne. c4(1)) call mpp_error(FATAL, "Real4 with pelist: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_pelist

  !> Test mpp_sum for large length variable with r4
  subroutine test_mpp_sum_large(pe,npes,root)
    integer, intent(in) :: npes, pe, root
    integer, parameter :: n=3000
    real(kind=r4_kind), dimension(n) :: a4
    ! Verify mpp_sum works for large array dim with random values
    call random_number(a4)
    call mpp_sync()
    call mpp_sum(a4(1:n),n)

    ! Verify mpp_sum works for large array dim with random values and mpp_sum length less than array dim
    call random_number(a4)
    call mpp_sync()
    call mpp_sum(a4(1:n),n/2)

  end subroutine test_mpp_sum_large

end program test_mpp_sum
